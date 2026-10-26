#include <stdlib.h>
#include "kdtree.h"

static double coord(const Vec3 *p, int ax) { return ax == 0 ? p->x : ax == 1 ? p->y : p->z; }
static double dist2(Vec3 a, Vec3 b) {
    double dx = a.x - b.x, dy = a.y - b.y, dz = a.z - b.z;
    return dx * dx + dy * dy + dz * dz;
}

/* Hoare quickselect: afterwards idx[k] holds the k-th smallest along ax, smaller-or-equal
 * entries are left of k and larger-or-equal right of k. lo,hi inclusive. Expected O(n). */
static void select_nth(const Vec3 *p, int *idx, int lo, int hi, int k, int ax) {
    while (lo < hi) {
        double a = coord(&p[idx[lo]], ax), b = coord(&p[idx[lo + (hi - lo) / 2]], ax), c = coord(&p[idx[hi]], ax);
        double pv = a < b ? (b < c ? b : (a < c ? c : a)) : (a < c ? a : (b < c ? c : b)); /* median of 3 */
        int i = lo, j = hi;
        while (i <= j) {
            while (coord(&p[idx[i]], ax) < pv) i++;
            while (coord(&p[idx[j]], ax) > pv) j--;
            if (i <= j) { int t = idx[i]; idx[i] = idx[j]; idx[j] = t; i++; j--; }
        }
        if (k <= j) hi = j;
        else if (k >= i) lo = i;
        else return;
    }
}

static void build_rec(KDTree *t, int lo, int hi) {   /* range [lo,hi) */
    if (hi - lo <= 1) { if (hi - lo == 1) t->axis[lo] = 0; return; }
    Vec3 mn = t->pts[t->idx[lo]], mx = mn;
    for (int i = lo + 1; i < hi; i++) {
        Vec3 p = t->pts[t->idx[i]];
        if (p.x < mn.x) mn.x = p.x;
        if (p.x > mx.x) mx.x = p.x;
        if (p.y < mn.y) mn.y = p.y;
        if (p.y > mx.y) mx.y = p.y;
        if (p.z < mn.z) mn.z = p.z;
        if (p.z > mx.z) mx.z = p.z;
    }
    double ex = mx.x - mn.x, ey = mx.y - mn.y, ez = mx.z - mn.z;
    int ax = (ex >= ey && ex >= ez) ? 0 : (ey >= ez ? 1 : 2);
    int mid = lo + (hi - lo) / 2;
    select_nth(t->pts, t->idx, lo, hi - 1, mid, ax);
    t->axis[mid] = (unsigned char)ax;
    build_rec(t, lo, mid);
    build_rec(t, mid + 1, hi);
}

int kd_build(KDTree *t, const Vec3 *pts, int n) {
    t->pts = pts; t->n = n; t->idx = NULL; t->axis = NULL;
    if (n <= 0) { t->n = 0; return 0; }
    t->idx = malloc(sizeof(int) * (size_t)n);
    t->axis = malloc((size_t)n);
    if (!t->idx || !t->axis) { kd_free(t); return -1; }
    for (int i = 0; i < n; i++) t->idx[i] = i;
    build_rec(t, 0, n);
    return 0;
}

static void rad_rec(const KDTree *t, int lo, int hi, Vec3 q, double r, double r2, SpatialVisit fn, void *ctx) {
    if (hi <= lo) return;
    int mid = lo + (hi - lo) / 2;
    int pi = t->idx[mid];
    Vec3 p = t->pts[pi];
    double d2 = dist2(p, q);
    if (d2 <= r2) fn(pi, d2, ctx);
    int ax = t->axis[mid];
    double qc = coord(&q, ax), pc = coord(&p, ax);
    if (qc - r <= pc) rad_rec(t, lo, mid, q, r, r2, fn, ctx);       /* left holds coords <= pc */
    if (qc + r >= pc) rad_rec(t, mid + 1, hi, q, r, r2, fn, ctx);   /* right holds coords >= pc */
}

void kd_radius(const KDTree *t, Vec3 q, double r, SpatialVisit fn, void *ctx) {
    if (t->n > 0 && r >= 0) rad_rec(t, 0, t->n, q, r, r * r, fn, ctx);
}

static void nn_rec(const KDTree *t, int lo, int hi, Vec3 q, int exclude, int *best, double *bd2) {
    if (hi <= lo) return;
    int mid = lo + (hi - lo) / 2;
    int pi = t->idx[mid];
    Vec3 p = t->pts[pi];
    if (pi != exclude) {
        double d2 = dist2(p, q);
        if (*best < 0 || d2 < *bd2) { *best = pi; *bd2 = d2; }
    }
    int ax = t->axis[mid];
    double diff = coord(&q, ax) - coord(&p, ax);
    if (diff <= 0) {
        nn_rec(t, lo, mid, q, exclude, best, bd2);
        if (*best < 0 || diff * diff <= *bd2) nn_rec(t, mid + 1, hi, q, exclude, best, bd2);
    } else {
        nn_rec(t, mid + 1, hi, q, exclude, best, bd2);
        if (*best < 0 || diff * diff <= *bd2) nn_rec(t, lo, mid, q, exclude, best, bd2);
    }
}

int kd_nearest(const KDTree *t, Vec3 q, int exclude, double *d2out) {
    int best = -1; double bd2 = 0;
    nn_rec(t, 0, t->n, q, exclude, &best, &bd2);
    if (best >= 0 && d2out) *d2out = bd2;
    return best;
}

size_t kd_memory(const KDTree *t) { return (size_t)t->n * (sizeof(int) + 1); }
void kd_free(KDTree *t) { free(t->idx); free(t->axis); t->idx = NULL; t->axis = NULL; t->n = 0; }
