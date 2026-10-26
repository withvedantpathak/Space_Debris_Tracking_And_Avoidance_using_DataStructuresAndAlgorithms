#include <math.h>
#include <stdlib.h>
#include "screen.h"
#include "ds/kdtree.h"
#include "ds/grid.h"

#define PARKED_X 1.0e8     /* prop_state parks unusable objects at x >= 1e9 */

typedef struct {
    const Prop *props; const Vec3 *cur, *nxt; const double *vm; const int *orig;
    const ScreenParams *sp; double t0, t1; int i; int k; int linear;   /* 0 bisection, 1 linear interpolation, 2 hybrid */
    DArray *out; ScreenStats *st; int err;
} Ctx;

/* Relative speed below which the hybrid refinement falls back to bisection (same value as the 'crossing' list of `analyze`). */
#define HYBRID_VREL_KMS 1.0

/* The per-pair test and refinement, identical to conj_bruteforce. */
static void test_pair(Ctx *c, int i, int j) {
    const double h = c->sp->step_s;
    const double thr = c->sp->infl_km ? c->sp->threshold_km + c->sp->infl_km[i] + c->sp->infl_km[j] : c->sp->threshold_km;
    c->st->pair_checks++;
    double dx = c->cur[i].x - c->cur[j].x, dy = c->cur[i].y - c->cur[j].y, dz = c->cur[i].z - c->cur[j].z;
    double L = c->vm[i] + c->vm[j], reach = thr + L * h;
    double d0sq = dx * dx + dy * dy + dz * dz;
    if (c->k == 0 && d0sq < thr * thr) c->st->close_at_start++;
    if (d0sq > reach * reach) return;
    double ex = c->nxt[i].x - c->nxt[j].x, ey = c->nxt[i].y - c->nxt[j].y, ez = c->nxt[i].z - c->nxt[j].z;
    double d0 = sqrt(d0sq), d1 = sqrt(ex * ex + ey * ey + ez * ez);
    if (d0 + d1 - L * h > 2.0 * thr) return;
    c->st->candidates++;
    double tm, dm, vr;
    c->st->refinements++;
    int found;
    if (c->linear) {          /* straight-line relative motion between the two step positions */
        int hybrid_slow = 0;
        double rx = dx, ry = dy, rz = dz, qx = ex - dx, qy = ey - dy, qz = ez - dz;
        double qq = qx * qx + qy * qy + qz * qz, sstar = qq > 0.0 ? -(rx * qx + ry * qy + rz * qz) / qq : -1.0;
        found = sstar > 0.0 && sstar <= 1.0;
        if (found) { double mx = rx + sstar * qx, my = ry + sstar * qy, mz = rz + sstar * qz; dm = sqrt(mx * mx + my * my + mz * mz); tm = c->t0 + sstar * h; vr = sqrt(qq) / h; }
        /* hybrid: the chord of a slow (co-moving) pair has flat minima at many steps (extra events, L33), so slow pairs
           (relative speed < HYBRID_VREL_KMS over the step) are refined by bisection instead */
        hybrid_slow = c->linear == 2 && sqrt(qq) / h < HYBRID_VREL_KMS;
        if (hybrid_slow) found = conj_refine(&c->props[i], &c->props[j], c->t0, c->t1, &tm, &dm, &vr);
    } else found = conj_refine(&c->props[i], &c->props[j], c->t0, c->t1, &tm, &dm, &vr);
    if (found && dm < thr && tm <= c->sp->window_s) {
        Conjunction cj = {i, j, tm, dm, vr};
        if (da_push(c->out, &cj) != 0) c->err = -3;
        c->st->events++;
    }
}

static void visit(int idx, double d2, void *vctx) {      /* idx = index into the compacted array */
    (void)d2;
    Ctx *c = vctx;
    int j = c->orig[idx];
    if (j > c->i) test_pair(c, c->i, j);
}

typedef struct { double lo, hi; int idx; } AltItem;
static int cmp_alt(const void *a, const void *b) { double x = ((const AltItem *)a)->lo, y = ((const AltItem *)b)->lo; return (x > y) - (x < y); }

int conj_fast(const Prop *props, int n, const ScreenParams *sp, FastScreen *fs, DArray *out, ScreenStats *st) {
    ScreenStats local = {0, 0, 0, 0, 0, 0};
    if (conj_check_params(sp, n) != 0) return -4;
    if (fs->cell_scale != 0.0 && !(fs->cell_scale >= 0.05 && fs->cell_scale <= 100.0)) return -4;   /* 0 = default; experiments use 0.25..4 */
    const double h = sp->step_s, thr = sp->threshold_km;
    const long K = (long)ceil(sp->window_s / h);
    size_t N = (size_t)(n > 0 ? n : 1);
    if (fs->method != SCREEN_KDTREE && fs->method != SCREEN_GRID && fs->method != SCREEN_ALT) return -4;
    if (fs->refine_linear < 0 || fs->refine_linear > 2) return -4;
    if (fs->method == SCREEN_ALT && (!fs->rmin_km || !fs->rmax_km)) return -4;
    Vec3 *cur = malloc(sizeof(Vec3) * N), *nxt = malloc(sizeof(Vec3) * N), *pts = malloc(sizeof(Vec3) * N);
    double *vm = malloc(sizeof(double) * N);
    int *orig = malloc(sizeof(int) * N);
    AltItem *alt = NULL;
    int rc = 0, m_alt = 0;
    if (!cur || !nxt || !pts || !vm || !orig) { rc = -3; goto done; }
    double vmax_all = 0.0;
    for (int i = 0; i < n; i++) { vm[i] = conj_vmax(&props[i]); if (prop_ok(&props[i]) && vm[i] > vmax_all) vmax_all = vm[i]; }
    fs->alt_violation_km = 0.0; fs->index_pairs = 0;
    double infl_max = 0.0;
    if (sp->infl_km) for (int i = 0; i < n; i++) if (sp->infl_km[i] > infl_max) infl_max = sp->infl_km[i];
    if (fs->method == SCREEN_ALT) {
        alt = malloc(sizeof(AltItem) * N);
        if (!alt) { rc = -3; goto done; }
        int m = 0;
        for (int i = 0; i < n; i++) if (prop_ok(&props[i])) { alt[m].lo = fs->rmin_km[i] - fs->alt_pad_km; alt[m].hi = fs->rmax_km[i] + fs->alt_pad_km; alt[m].idx = i; m++; }
        qsort(alt, (size_t)m, sizeof(AltItem), cmp_alt);
        m_alt = m;
    }
    prop_positions(props, n, 0.0, cur);
    Ctx c = {props, cur, nxt, vm, orig, sp, 0, 0, 0, 0, fs->refine_linear, out, &local, 0};
    size_t first = out->n;
    for (long k = 0; k < K; k++) {
        c.t0 = (double)k * h; c.t1 = (double)(k + 1) * h; c.k = (int)k; c.cur = cur; c.nxt = nxt;
        prop_positions(props, n, c.t1, nxt);
        local.steps++;
        if (fs->method == SCREEN_ALT) {
            for (int a = 0; a < m_alt && !c.err; a++) {
                int i = alt[a].idx;
                double r = sqrt(cur[i].x * cur[i].x + cur[i].y * cur[i].y + cur[i].z * cur[i].z);
                if (r < 1.0e7) {                                   /* ignore parked */
                    double v1 = fs->rmin_km[i] - r, v2 = r - fs->rmax_km[i];
                    if (v1 > fs->alt_violation_km) fs->alt_violation_km = v1;
                    if (v2 > fs->alt_violation_km) fs->alt_violation_km = v2;
                }
                for (int b = a + 1; b < m_alt && alt[b].lo <= alt[a].hi + thr + 2.0 * infl_max; b++) {
                    int j = alt[b].idx;
                    int lo = i < j ? i : j, hi = i < j ? j : i;
                    if (cur[lo].x >= PARKED_X || cur[hi].x >= PARKED_X) continue;
                    test_pair(&c, lo, hi);
                }
            }
        } else {
            int m = 0;
            for (int i = 0; i < n; i++) if (prop_ok(&props[i]) && cur[i].x < PARKED_X) { pts[m] = cur[i]; orig[m] = i; m++; }
            double Rmax = thr + 2.0 * infl_max + 2.0 * vmax_all * h, cs = fs->cell_scale > 0.0 ? fs->cell_scale : 1.0;
            KDTree kd; Grid gr; int built;
            if (fs->method == SCREEN_KDTREE) built = kd_build(&kd, pts, m); else built = grid_build(&gr, pts, m, fmax(Rmax * cs, 1.0));   /* floor: Rmax can be 0 (empty catalog, zero threshold); exactness does not depend on the cell size */
            if (built != 0) { rc = built == -1 ? -3 : -4; break; }
            for (int a = 0; a < m && !c.err; a++) {
                int i = orig[a]; c.i = i;
                double R = thr + (sp->infl_km ? sp->infl_km[i] : 0.0) + infl_max + (vm[i] + vmax_all) * h;
                if (fs->method == SCREEN_KDTREE) kd_radius(&kd, pts[a], R, visit, &c); else grid_radius(&gr, pts[a], R, visit, &c);
            }
            if (fs->method == SCREEN_KDTREE) kd_free(&kd); else grid_free(&gr);
        }
        if (c.err) { rc = c.err; break; }
        Vec3 *tmp = cur; cur = nxt; nxt = tmp;
    }
    if (rc == 0 && out->n > first) conj_sort_range(out, first);
    fs->index_pairs = local.pair_checks;
    if (st) *st = local;
done:
    free(cur); free(nxt); free(pts); free(vm); free(orig); free(alt);
    return rc;
}
