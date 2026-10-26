#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "conjunction.h"

double conj_vmax(const Prop *p) {
    return 1.01 * sqrt(EARTH_MU_KM3S2 / p->a) * sqrt((1.0 + p->e) / (1.0 - p->e));
}

static double gfun(const Prop *pa, const Prop *pb, double t, double *d, double *vr) {
    Vec3 ra, va, rb, vb;
    prop_state(pa, t, &ra, &va);
    prop_state(pb, t, &rb, &vb);
    double dx = ra.x - rb.x, dy = ra.y - rb.y, dz = ra.z - rb.z;
    double ux = va.x - vb.x, uy = va.y - vb.y, uz = va.z - vb.z;
    if (d) *d = sqrt(dx * dx + dy * dy + dz * dz);
    if (vr) *vr = sqrt(ux * ux + uy * uy + uz * uz);
    return dx * ux + dy * uy + dz * uz;                        /* = 0.5 d(d^2)/dt */
}

int conj_refine(const Prop *pa, const Prop *pb, double t0, double t1, double *tmin, double *dmin, double *vrel) {
    double g0 = gfun(pa, pb, t0, NULL, NULL), g1 = gfun(pa, pb, t1, NULL, NULL);
    if (!(g0 < 0.0 && g1 >= 0.0)) return 0;
    double lo = t0, hi = t1;
    for (int i = 0; i < 60; i++) {                              /* bisection on the sign of g */
        double mid = 0.5 * (lo + hi);
        if (gfun(pa, pb, mid, NULL, NULL) < 0.0) lo = mid; else hi = mid;
    }
    *tmin = 0.5 * (lo + hi);
    gfun(pa, pb, *tmin, dmin, vrel);
    return 1;
}

static int cmp_conj(const void *x, const void *y) {
    const Conjunction *p = x, *q = y;
    if (p->a != q->a) return p->a < q->a ? -1 : 1;
    if (p->b != q->b) return p->b < q->b ? -1 : 1;
    return (p->tca > q->tca) - (p->tca < q->tca);
}
void conj_sort(DArray *ev) { if (ev->n > 1) qsort(ev->data, ev->n, sizeof(Conjunction), cmp_conj); }

void conj_sort_range(DArray *ev, size_t first) {
    if (ev->n > first + 1) qsort((Conjunction *)ev->data + first, ev->n - first, sizeof(Conjunction), cmp_conj);
}

int conj_check_params(const ScreenParams *sp, int n) {
    if (!sp || n < 0) return -4;
    if (!(sp->step_s > 0.0) || !isfinite(sp->step_s)) return -4;                    /* also rejects NaN */
    if (!(sp->window_s >= 0.0) || !isfinite(sp->window_s)) return -4;
    if (!(sp->threshold_km >= 0.0) || !isfinite(sp->threshold_km)) return -4;
    if (!(sp->window_s / sp->step_s <= 2.0e9)) return -4;                           /* step count must fit a long on every platform */
    if (sp->infl_km)
        for (int i = 0; i < n; i++)
            if (!(sp->infl_km[i] >= 0.0) || !isfinite(sp->infl_km[i])) return -4;
    return 0;
}

int conj_bruteforce(const Prop *props, int n, const ScreenParams *sp, DArray *out, ScreenStats *st) {
    ScreenStats local = {0, 0, 0, 0, 0, 0};
    if (conj_check_params(sp, n) != 0) return -4;
    double h = sp->step_s, thr = sp->threshold_km;
    long K = (long)ceil(sp->window_s / h);
    Vec3 *cur = malloc(sizeof(Vec3) * (size_t)(n > 0 ? n : 1)), *nxt = malloc(sizeof(Vec3) * (size_t)(n > 0 ? n : 1));
    double *vm = malloc(sizeof(double) * (size_t)(n > 0 ? n : 1));
    if (!cur || !nxt || !vm) { free(cur); free(nxt); free(vm); return -3; }
    for (int i = 0; i < n; i++) vm[i] = conj_vmax(&props[i]);
    prop_positions(props, n, 0.0, cur);
    size_t first = out->n;
    for (long k = 0; k < K; k++) {
        double t0 = (double)k * h, t1 = (double)(k + 1) * h;
        prop_positions(props, n, t1, nxt);
        local.steps++;
        for (int i = 0; i < n; i++) for (int j = i + 1; j < n; j++) {
            local.pair_checks++;
            double dx = cur[i].x - cur[j].x, dy = cur[i].y - cur[j].y, dz = cur[i].z - cur[j].z;
            double tij = sp->infl_km ? thr + sp->infl_km[i] + sp->infl_km[j] : thr;
            double L = vm[i] + vm[j], reach = tij + L * h;      /* d1 >= d0 - L h, so d0 > thr + L h => skip */
            double d0sq = dx * dx + dy * dy + dz * dz;
            if (k == 0 && d0sq < tij * tij) local.close_at_start++;
            if (d0sq > reach * reach) continue;
            double ex = nxt[i].x - nxt[j].x, ey = nxt[i].y - nxt[j].y, ez = nxt[i].z - nxt[j].z;
            double d0 = sqrt(d0sq), d1 = sqrt(ex * ex + ey * ey + ez * ez);
            if (d0 + d1 - L * h > 2.0 * tij) continue;         /* proven: min distance in interval > thr */
            local.candidates++;
            double tm, dm, vr;
            local.refinements++;
            if (conj_refine(&props[i], &props[j], t0, t1, &tm, &dm, &vr) && dm < tij && tm <= sp->window_s) {
                Conjunction c = {i, j, tm, dm, vr};
                if (da_push(out, &c) != 0) { free(cur); free(nxt); free(vm); return -3; }
                local.events++;
            }
        }
        Vec3 *tmp = cur; cur = nxt; nxt = tmp;
    }
    free(cur); free(nxt); free(vm);
    if (out->n > first) {                                       /* sort only the newly added block */
        qsort((Conjunction *)out->data + first, out->n - first, sizeof(Conjunction), cmp_conj);
    }
    if (st) *st = local;
    return 0;
}

void conj_compare(const DArray *truth, const DArray *test, double tol, CompareResult *res) {
    res->truth = (long)truth->n; res->test = (long)test->n; res->matched = 0;
    const Conjunction *T = truth->data, *E = test->data;
    unsigned char *used = calloc(test->n ? test->n : 1, 1);
    size_t j0 = 0;
    for (size_t i = 0; i < truth->n; i++) {
        while (j0 < test->n && (E[j0].a < T[i].a || (E[j0].a == T[i].a && E[j0].b < T[i].b))) j0++;
        for (size_t j = j0; j < test->n && E[j].a == T[i].a && E[j].b == T[i].b; j++) {
            if (!used[j] && fabs(E[j].tca - T[i].tca) <= tol) { used[j] = 1; res->matched++; break; }
        }
    }
    res->missed = res->truth - res->matched;
    res->extra = res->test - res->matched;
    free(used);
}

int conj_best_per_pair(const DArray *in, DArray *out) {
    const Conjunction *c = in->data;
    for (size_t i = 0; i < in->n; ) {
        size_t best = i, j = i;
        while (j < in->n && c[j].a == c[i].a && c[j].b == c[i].b) { if (c[j].miss_km < c[best].miss_km) best = j; j++; }
        if (da_push(out, &c[best]) != 0) return -3;
        i = j;
    }
    return 0;
}

int conj_write_csv(const DArray *ev, const Object *objs, const char *path) {
    FILE *f = fopen(path, "w");
    if (!f) return -1;
    fprintf(f, "norad_id_a,norad_id_b,tca_s,miss_km,rel_speed_kms\n");
    const Conjunction *c = ev->data;
    for (size_t i = 0; i < ev->n; i++)
        fprintf(f, "%d,%d,%.6f,%.9f,%.6f\n", objs[c[i].a].norad_id, objs[c[i].b].norad_id, c[i].tca, c[i].miss_km, c[i].rel_speed_kms);
    return fclose(f) == 0 ? 0 : -1;
}
