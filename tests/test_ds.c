/* Unit tests for src/ds: every structure is compared with a naive reference on seeded
 * random data. Plain C, no framework. Exit code 0 only if every check passes. */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "ds/darray.h"
#include "ds/heap.h"
#include "ds/queue.h"
#include "ds/hash.h"
#include "ds/avl.h"
#include "ds/kdtree.h"
#include "ds/grid.h"
#include "ds/dsu.h"

static long checks = 0, fails = 0;
#define CHECK(c) do { checks++; if (!(c)) { fails++; printf("FAIL %s:%d: %s\n", __FILE__, __LINE__, #c); } } while (0)

static unsigned long long rs = 88172645463325252ULL;           /* fixed seed, xorshift64 */
static unsigned long long rnd(void) { rs ^= rs << 13; rs ^= rs >> 7; rs ^= rs << 17; return rs; }
static int rint_(int n) { return (int)(rnd() % (unsigned long long)n); }
static double runif(double lo, double hi) { return lo + (hi - lo) * ((double)(rnd() >> 11) / 9007199254740992.0); }

/* ---------- dynamic array ---------- */
static void test_darray(void) {
    DArray a; da_init(&a, sizeof(int));
    for (int i = 0; i < 10000; i++) CHECK(da_push(&a, &i) == 0);
    CHECK(a.n == 10000);
    int ok = 1;
    for (int i = 0; i < 10000; i++) if (*(int *)da_at(&a, (size_t)i) != i) ok = 0;
    CHECK(ok);
    da_clear(&a); CHECK(a.n == 0);
    da_free(&a);
}

/* ---------- heap ---------- */
static int cmp_d(const void *a, const void *b) { double x = *(const double *)a, y = *(const double *)b; return (x > y) - (x < y); }
static void test_heap(void) {
    Heap h; CHECK(heap_init(&h, 2) == 0);
    HeapItem it;
    CHECK(heap_pop(&h, &it) == 0 && heap_peek(&h, &it) == 0);
    enum { N = 5000 };
    double *ref = malloc(sizeof(double) * N);
    for (int i = 0; i < N; i++) { ref[i] = (double)rint_(300) / 7.0 - 20.0; heap_push(&h, ref[i], i); }
    qsort(ref, N, sizeof(double), cmp_d);
    int ok = 1;
    for (int i = 0; i < N; i++) { if (!heap_pop(&h, &it) || it.key != ref[i]) ok = 0; }
    CHECK(ok); CHECK(h.n == 0);
    /* interleaved push/pop vs naive minimum scan */
    double naive[600]; int m = 0; ok = 1;
    for (int step = 0; step < 4000; step++) {
        if (m == 0 || (m < 600 && rint_(100) < 55)) {
            double k = runif(-1000, 1000); naive[m++] = k; heap_push(&h, k, step);
        } else {
            int bi = 0; for (int i = 1; i < m; i++) if (naive[i] < naive[bi]) bi = i;
            heap_peek(&h, &it); double pk = it.key;
            heap_pop(&h, &it);
            if (it.key != naive[bi] || pk != it.key) ok = 0;
            naive[bi] = naive[--m];
        }
    }
    CHECK(ok);
    free(ref); heap_free(&h);
}

/* ---------- queue ---------- */
static void test_queue(void) {
    Queue q; queue_init(&q);
    int ref[2000], head = 0, tail = 0, ok = 1, v;
    CHECK(queue_empty(&q) && queue_pop(&q, &v) == 0);
    for (int step = 0; step < 20000; step++) {
        if (tail - head < 1500 && rint_(100) < 52) { ref[tail++ % 2000] = step; queue_push(&q, step); }
        else if (tail > head) { if (!queue_pop(&q, &v) || v != ref[head++ % 2000]) ok = 0; }
        else if (queue_pop(&q, &v)) ok = 0;
        if (q.n != (size_t)(tail - head)) ok = 0;
    }
    CHECK(ok);
    queue_free(&q); CHECK(queue_empty(&q));
}

/* ---------- hash ---------- */
static void test_hash(void) {
    enum { D = 3000 };
    Hash h; CHECK(hash_init(&h, 1) == 0);           /* tiny start forces several resizes */
    static int64_t keys[D], val[D]; static char present[D];
    for (int i = 0; i < D; i++) keys[i] = ((int64_t)rint_(1000000) - 500000) * 4294967311LL + i; /* wide, negative, unique */
    int ok = 1; size_t cnt = 0; int64_t got;
    for (int step = 0; step < 60000; step++) {
        int i = rint_(D), op = rint_(10);
        if (op < 5) { val[i] = (int64_t)rnd(); if (hash_put(&h, keys[i], val[i]) != 0) ok = 0; if (!present[i]) { present[i] = 1; cnt++; } }
        else if (op < 8) { if (hash_remove(&h, keys[i]) != present[i]) ok = 0; if (present[i]) { present[i] = 0; cnt--; } }
        else { int f = hash_get(&h, keys[i], &got); if (f != present[i] || (f && got != val[i])) ok = 0; }
        if (hash_count(&h) != cnt) ok = 0;
    }
    CHECK(ok);
    for (int i = 0; i < D; i++) { int f = hash_get(&h, keys[i], &got); if (f != present[i] || (f && got != val[i])) ok = 0; }
    CHECK(ok); CHECK(h.nb >= 8 && h.n <= h.nb);
    hash_free(&h);
}

/* ---------- AVL ---------- */
static void collect(int key, void *val, void *ctx) { (void)val; int *p = ctx; p[++p[0]] = key; }
static void test_avl(void) {
    enum { K = 1000 };
    Avl t; avl_init(&t);
    static char present[K]; static int stored[K]; int ok = 1; void *v;
    for (int step = 0; step < 20000; step++) {
        int k = rint_(K), op = rint_(10);
        if (op < 5) { stored[k] = step; if (avl_insert(&t, k, &stored[k]) != 0) ok = 0; present[k] = 1; }
        else if (op < 8) { if (avl_remove(&t, k) != present[k]) ok = 0; present[k] = 0; }
        else { int f = avl_find(&t, k, &v); if (f != present[k] || (f && v != &stored[k])) ok = 0; }
        if (step % 50 == 0 && !avl_check(&t)) ok = 0;
    }
    CHECK(ok); CHECK(avl_check(&t));
    static int out[K + 1]; out[0] = 0; avl_foreach(&t, collect, out);
    int idx = 0, sorted = 1, prev = -1;
    for (int k = 0; k < K; k++) if (present[k]) { idx++; if (out[idx] != k || k <= prev) sorted = 0; prev = k; }
    CHECK(sorted && out[0] == idx && t.n == (size_t)idx);
    avl_free(&t);
    /* worst case for plain BSTs: sorted insert; height must stay near log2(n) */
    avl_init(&t);
    for (int i = 1; i <= 100000; i++) avl_insert(&t, i, NULL);
    CHECK(avl_check(&t)); CHECK(avl_height(&t) <= 25);
    for (int i = 1; i <= 100000; i += 2) avl_remove(&t, i);
    CHECK(avl_check(&t)); CHECK(t.n == 50000);
    avl_free(&t);
}

/* ---------- spatial structures ---------- */
typedef struct { unsigned char *seen; int dup; int count; } Mark;
static void mark_visit(int idx, double d2, void *ctx) { (void)d2; Mark *m = ctx; if (m->seen[idx]) m->dup = 1; m->seen[idx] = 1; m->count++; }

static double bd2(Vec3 a, Vec3 b) { double dx = a.x - b.x, dy = a.y - b.y, dz = a.z - b.z; return dx * dx + dy * dy + dz * dz; }

/* mode: 0 spread cube, 1 many exact duplicates, 2 thin shell (orbit-like) */
static void gen_points(Vec3 *p, int n, int mode) {
    for (int i = 0; i < n; i++) {
        if (mode == 1) { p[i].x = (double)rint_(5) * 10; p[i].y = (double)rint_(5) * 10; p[i].z = (double)rint_(5) * 10; }
        else if (mode == 2) {                       /* points scaled onto a thin shell, radius 6700-7300 */
            double r = runif(6700, 7300);
            double dx = runif(-1, 1), dy = runif(-1, 1), dz = runif(-1, 1);
            double n2 = dx * dx + dy * dy + dz * dz + 1e-12, g = 1;
            for (int k = 0; k < 40; k++) g = 0.5 * (g + n2 / g);   /* Newton sqrt, avoids libm */
            p[i].x = r * dx / g; p[i].y = r * dy / g; p[i].z = r * dz / g;
        } else { p[i].x = runif(-7000, 7000); p[i].y = runif(-7000, 7000); p[i].z = runif(-7000, 7000); }
    }
}

static void test_spatial(void) {
    int sizes[] = {0, 1, 2, 10, 500, 3000};
    for (int mode = 0; mode < 3; mode++) for (int si = 0; si < 6; si++) {
        int n = sizes[si];
        Vec3 *p = malloc(sizeof(Vec3) * (size_t)(n ? n : 1));
        gen_points(p, n, mode);
        KDTree t; CHECK(kd_build(&t, p, n) == 0);
        double cells[] = {5.0, 40.0, 400.0};
        Grid g[3]; for (int c = 0; c < 3; c++) CHECK(grid_build(&g[c], p, n, cells[c]) == 0);
        unsigned char *seen = malloc((size_t)(n ? n : 1));
        int ok = 1;
        for (int qn = 0; qn < 60; qn++) {
            Vec3 q; if (n > 0 && rint_(2)) { q = p[rint_(n)]; q.x += runif(-3, 3); } else gen_points(&q, 1, mode);
            double rs_[] = {0.0, 1.0, 5.0, 40.0, 120.0, 900.0};
            double r = rs_[rint_(6)];
            int expect = 0; for (int i = 0; i < n; i++) if (bd2(p[i], q) <= r * r) expect++;
            Mark m; m.seen = seen; m.dup = 0; m.count = 0; memset(seen, 0, (size_t)(n ? n : 1));
            kd_radius(&t, q, r, mark_visit, &m);
            if (m.count != expect || m.dup) ok = 0;
            for (int i = 0; i < n; i++) if ((bd2(p[i], q) <= r * r) != (seen[i] != 0)) ok = 0;
            for (int c = 0; c < 3; c++) {
                if (r / cells[c] > 30) continue;      /* absurd cell/radius ratio: millions of empty cells */
                memset(seen, 0, (size_t)(n ? n : 1)); m.dup = 0; m.count = 0;
                grid_radius(&g[c], q, r, mark_visit, &m);
                if (m.count != expect || m.dup) ok = 0;
                for (int i = 0; i < n; i++) if ((bd2(p[i], q) <= r * r) != (seen[i] != 0)) ok = 0;
            }
            /* nearest neighbour, with and without exclusion */
            int ex = (n > 0 && rint_(2)) ? rint_(n) : -1;
            double best = -1; int bi = -1;
            for (int i = 0; i < n; i++) if (i != ex) { double d = bd2(p[i], q); if (bi < 0 || d < best) { best = d; bi = i; } }
            double kd2 = 0; int ki = kd_nearest(&t, q, ex, &kd2);
            if ((ki < 0) != (bi < 0) || (bi >= 0 && (kd2 != best || ki == ex))) ok = 0;
        }
        CHECK(ok);
        CHECK(kd_memory(&t) == (size_t)n * (sizeof(int) + 1));
        for (int c = 0; c < 3; c++) { CHECK(grid_memory(&g[c]) > 0 || n == 0); grid_free(&g[c]); }
        kd_free(&t); free(seen); free(p);
    }
    /* grid rejects bad cell sizes and out-of-range coordinates */
    Vec3 far = {1e12, 0, 0}; Grid g;
    CHECK(grid_build(&g, &far, 1, 1.0) == -2);
    CHECK(grid_build(&g, &far, 1, 0.0) == -2);
}

/* ---------- DSU ---------- */
static void test_dsu(void) {
    enum { N = 400 };
    Dsu d; CHECK(dsu_init(&d, N) == 0);
    int label[N]; for (int i = 0; i < N; i++) label[i] = i;
    int ok = 1, comps = N;
    for (int step = 0; step < 600; step++) {
        int a = rint_(N), b = rint_(N);
        int merged = dsu_union(&d, a, b);
        if (label[a] != label[b]) { int from = label[b], to = label[a]; for (int i = 0; i < N; i++) if (label[i] == from) label[i] = to; comps--; if (!merged) ok = 0; }
        else if (merged) ok = 0;
        if (d.comps != comps) ok = 0;
    }
    for (int i = 0; i < N; i++) {
        int sz = 0; for (int j = 0; j < N; j++) if (label[j] == label[i]) sz++;
        if (dsu_size(&d, i) != sz) ok = 0;
        int j = rint_(N); if (dsu_same(&d, i, j) != (label[i] == label[j])) ok = 0;
    }
    CHECK(ok);
    dsu_free(&d);
}

int main(void) {
    test_darray(); test_heap(); test_queue(); test_hash(); test_avl(); test_spatial(); test_dsu();
    printf("%ld checks, %ld failed\n", checks, fails);
    return fails ? 1 : 0;
}
