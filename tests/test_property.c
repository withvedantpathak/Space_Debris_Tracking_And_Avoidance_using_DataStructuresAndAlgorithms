/* C8 property tests. Complements the example-based unit tests with
 *   A. structural invariants of the hand-written data structures (checked directly on the internal layout),
 *   B. CSV-parser fuzzing with an accounting identity, scope-filter and synthetic-generator properties,
 *   C. randomised differential tests (k-d tree / grid / altitude filter versus the brute-force baseline) over many
 *      configurations, plus metamorphic relations (permutation, threshold monotonicity, window prefix, determinism),
 *   D. mitigation and graph properties (exact cover versus subset enumeration, with and without reductions).
 * All randomness is seeded; every failure prints the trial number so it can be replayed. Plain C99, no frameworks. */
#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "catalog.h"
#include "synth.h"
#include "propagate.h"
#include "conjunction.h"
#include "screen.h"
#include "graph.h"
#include "mitigate.h"
#include "rank.h"
#include "analyze.h"
#include "ds/darray.h"
#include "ds/heap.h"
#include "ds/hash.h"
#include "ds/avl.h"
#include "ds/kdtree.h"
#include "ds/grid.h"
#include "ds/dsu.h"

static long checks = 0, fails = 0;
#define CHECK(c) do { checks++; if (!(c)) { fails++; printf("FAIL %s:%d: %s\n", __FILE__, __LINE__, #c); } } while (0)
#define CHECKF(c, ...) do { checks++; if (!(c)) { fails++; printf("FAIL %s:%d: %s | ", __FILE__, __LINE__, #c); printf(__VA_ARGS__); printf("\n"); } } while (0)

static unsigned long long rs = 0x2545F4914F6CDD1DULL;                /* xorshift64*, reseeded per trial */
static void seed_rng(unsigned long long s) { rs = s * 0x9E3779B97F4A7C15ULL + 1ULL; }
static unsigned long long rnd(void) { rs ^= rs >> 12; rs ^= rs << 25; rs ^= rs >> 27; return rs * 2685821657736338717ULL; }
static int rint_(int n) { return (int)(rnd() % (unsigned long long)n); }
static double runif(double lo, double hi) { return lo + (hi - lo) * ((double)(rnd() >> 11) / 9007199254740992.0); }

/* ======================================================================================================== A */
static int heap_ok(const Heap *h) {
    for (int i = 1; i < h->n; i++) if (h->a[(i - 1) / 2].key > h->a[i].key) return 0;
    return 1;
}
static void test_heap_invariant(void) {
    for (int trial = 0; trial < 60; trial++) {
        seed_rng(1000 + (unsigned)trial);
        Heap h; CHECK(heap_init(&h, 1 + rint_(8)) == 0);
        int ok = 1, cnt = 0; long long sum = 0;
        for (int step = 0; step < 400; step++) {
            if (h.n == 0 || rint_(100) < 60) { int k = rint_(50) - 25; heap_push(&h, (double)k, step); sum += k; cnt++; }
            else {
                HeapItem pk, it; heap_peek(&h, &pk); heap_pop(&h, &it);
                if (pk.key != it.key) ok = 0;
                for (int i = 0; i < h.n; i++) if (h.a[i].key < it.key) ok = 0;      /* popped key is a minimum */
                sum -= (long long)it.key; cnt--;
            }
            if (!heap_ok(&h) || h.n != cnt) ok = 0;
        }
        double prev = -1e300; long long s2 = 0; HeapItem it;                       /* drain: ascending, same multiset (sum, count) */
        while (heap_pop(&h, &it)) { if (it.key < prev) ok = 0; prev = it.key; s2 += (long long)it.key; cnt--; }
        CHECKF(ok && s2 == sum && cnt == 0, "heap trial %d", trial);
        heap_free(&h);
    }
}

static double cv(Vec3 p, int ax) { return ax == 0 ? p.x : ax == 1 ? p.y : p.z; }
static int kd_rec_ok(const KDTree *t, int lo, int hi) {                             /* node of [lo,hi) sits at mid */
    if (hi - lo <= 1) return 1;
    int mid = lo + (hi - lo) / 2, ax = t->axis[mid];
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
    int widest = (ex >= ey && ex >= ez) ? 0 : (ey >= ez ? 1 : 2);
    if (ax != widest) return 0;                                                     /* documented split rule */
    double pv = cv(t->pts[t->idx[mid]], ax);
    for (int i = lo; i < mid; i++) if (cv(t->pts[t->idx[i]], ax) > pv) return 0;    /* left <= pivot */
    for (int i = mid + 1; i < hi; i++) if (cv(t->pts[t->idx[i]], ax) < pv) return 0; /* right >= pivot */
    return kd_rec_ok(t, lo, mid) && kd_rec_ok(t, mid + 1, hi);
}
/* point sets: 0 cube, 1 lattice with many duplicates, 2 thin shell, 3 collinear, 4 all identical, 5 two near-equal clusters */
static void gen_pts(Vec3 *p, int n, int mode) {
    for (int i = 0; i < n; i++) {
        switch (mode) {
        case 1: p[i].x = (double)rint_(4) * 25; p[i].y = (double)rint_(4) * 25; p[i].z = (double)rint_(4) * 25; break;
        case 2: { double th = runif(0, 6.283185307179586), z = runif(-1, 1), r = runif(6700, 7300), s = sqrt(1 - z * z);
                  p[i].x = r * s * cos(th); p[i].y = r * s * sin(th); p[i].z = r * z; break; }
        case 3: p[i].x = 5; p[i].y = -3; p[i].z = runif(-100, 100); break;
        case 4: p[i].x = 1; p[i].y = 2; p[i].z = 3; break;
        case 5: p[i].x = 7000 + (rint_(2) ? 1e-9 * (double)rint_(5) : 0); p[i].y = runif(-1, 1); p[i].z = runif(-1, 1); break;
        default: p[i].x = runif(-7000, 7000); p[i].y = runif(-7000, 7000); p[i].z = runif(-7000, 7000);
        }
    }
}
typedef struct { int target, hit; } SelfCtx;
static void self_visit(int idx, double d2, void *ctx) { (void)d2; SelfCtx *c = ctx; if (idx == c->target) c->hit = 1; }
static void test_kd_invariant(void) {
    int sizes[] = {0, 1, 2, 3, 7, 64, 513}, trials = 0;
    for (int mode = 0; mode < 6; mode++) for (int si = 0; si < 7; si++) for (int rep = 0; rep < 3; rep++) {
        seed_rng(2000 + (unsigned)(mode * 100 + si * 10 + rep));
        int n = sizes[si];
        Vec3 *p = malloc(sizeof(Vec3) * (size_t)(n ? n : 1)); gen_pts(p, n, mode);
        KDTree t; CHECK(kd_build(&t, p, n) == 0);
        int ok = 1; unsigned char *seen = calloc((size_t)(n ? n : 1), 1);
        for (int i = 0; i < n; i++) { if (t.idx[i] < 0 || t.idx[i] >= n || seen[t.idx[i]]) ok = 0; else seen[t.idx[i]] = 1; }   /* permutation */
        for (int i = 0; i < n; i++) if (!seen[i]) ok = 0;
        if (!kd_rec_ok(&t, 0, n)) ok = 0;
        for (int i = 0; i < n && ok; i++) { SelfCtx cx = {i, 0}; kd_radius(&t, p[i], 0.0, self_visit, &cx); if (!cx.hit) ok = 0; }   /* a zero-radius query at a stored point finds it */
        CHECKF(ok, "kd mode %d n %d rep %d", mode, n, rep);
        trials++;
        kd_free(&t); free(seen); free(p);
    }
    CHECK(trials == 6 * 7 * 3);
}

static int grid_ok(const Grid *g, const Vec3 *p, int n, double cell) {
    int ok = 1; int *cnt = calloc((size_t)(n > 0 ? n : 1), sizeof(int)); size_t nodes = 0;
    for (size_t b = 0; b < g->cells.nb; b++) for (HNode *h = g->cells.b[b]; h; h = h->next) {
        nodes++;
        int ix = (int)((h->key >> 42) & 0x1FFFFF) - (1 << 20), iy = (int)((h->key >> 21) & 0x1FFFFF) - (1 << 20), iz = (int)(h->key & 0x1FFFFF) - (1 << 20);
        for (int i = (int)h->val; i >= 0; i = g->next[i]) {
            if (i >= n) { ok = 0; break; }
            if (++cnt[i] > 1) { ok = 0; break; }                                    /* a point on two chains or a cycle */
            if ((int)floor(p[i].x / cell) != ix || (int)floor(p[i].y / cell) != iy || (int)floor(p[i].z / cell) != iz) ok = 0;
        }
    }
    for (int i = 0; i < n; i++) if (cnt[i] != 1) ok = 0;                            /* every point stored exactly once, in its own cell */
    if (nodes != hash_count(&g->cells)) ok = 0;
    free(cnt);
    return ok;
}
static void test_grid_invariant(void) {
    int sizes[] = {0, 1, 5, 200, 1500}; double cells[] = {0.5, 7.0, 90.0, 1000.0};
    for (int mode = 0; mode < 6; mode++) for (int si = 0; si < 5; si++) for (int ci = 0; ci < 4; ci++) {
        seed_rng(3000 + (unsigned)(mode * 100 + si * 10 + ci));
        int n = sizes[si]; double cell = cells[ci];
        Vec3 *p = malloc(sizeof(Vec3) * (size_t)(n ? n : 1)); gen_pts(p, n, mode);
        if (mode == 1) for (int i = 0; i < n; i++) { p[i].x = cell * (double)(rint_(9) - 4); p[i].y = cell * (double)(rint_(9) - 4); p[i].z = cell * (double)(rint_(9) - 4); }   /* exactly on cell boundaries, both signs */
        Grid g; int rc = grid_build(&g, p, n, cell);
        int range_ok = 1;                                                           /* the shell at cell 0.5 needs more than 21 bits: must be refused, not corrupted */
        for (int i = 0; i < n; i++) if (fabs(p[i].x / cell) >= 1048576.0 || fabs(p[i].y / cell) >= 1048576.0 || fabs(p[i].z / cell) >= 1048576.0) range_ok = 0;
        if (!range_ok) { CHECKF(rc == -2, "grid mode %d cell %g should refuse out-of-range", mode, cell); free(p); continue; }
        CHECKF(rc == 0, "grid_build rc %d (mode %d n %d cell %g)", rc, mode, n, cell);
        if (rc == 0) {
            CHECKF(grid_ok(&g, p, n, cell), "grid mode %d n %d cell %g", mode, n, cell);
            int found_self = 1;                                                      /* a zero-radius query at each stored point finds that point */
            for (int i = 0; i < n && found_self; i++) { SelfCtx cx = {i, 0}; grid_radius(&g, p[i], 0.0, self_visit, &cx); if (!cx.hit) found_self = 0; }
            CHECKF(found_self, "grid self-query mode %d n %d cell %g", mode, n, cell);
            grid_free(&g);
        }
        free(p);
    }
}

static int dsu_ok(Dsu *d) {
    int ok = 1, roots = 0; int *cnt = calloc((size_t)(d->n ? d->n : 1), sizeof(int)); long total = 0;
    for (int x = 0; x < d->n; x++) {
        int y = x, steps = 0;
        while (d->parent[y] != y && steps <= d->n) { y = d->parent[y]; steps++; }
        if (steps > d->n) { ok = 0; break; }                                         /* cycle */
        cnt[y]++;
    }
    for (int x = 0; x < d->n && ok; x++) if (d->parent[x] == x) { roots++; if (d->size[x] != cnt[x]) ok = 0; total += d->size[x]; }
    if (ok && (roots != d->comps || total != d->n)) ok = 0;
    free(cnt);
    return ok;
}
static void test_dsu_invariant(void) {
    for (int trial = 0; trial < 80; trial++) {
        seed_rng(4000 + (unsigned)trial);
        int n = 1 + rint_(60); Dsu d; CHECK(dsu_init(&d, n) == 0);
        int ok = dsu_ok(&d);
        for (int step = 0; step < 4 * n; step++) {
            int a = rint_(n), b = rint_(n), before = d.comps, same = dsu_same(&d, a, b);
            int ra = dsu_find(&d, a), rb = dsu_find(&d, b), sa = d.size[ra], sb = d.size[rb];
            int merged = dsu_union(&d, a, b);
            if (merged != !same || d.comps != before - merged || !dsu_same(&d, a, b)) ok = 0;
            if (merged && sa != sb && dsu_find(&d, a) != (sa > sb ? ra : rb)) ok = 0;      /* union by size: the larger set's root survives */
            if (merged && d.size[dsu_find(&d, a)] != sa + sb) ok = 0;
            if (!dsu_ok(&d)) ok = 0;
            int labels_before[64], labels_after[64];                                 /* find() with path halving must not change the partition */
            for (int x = 0; x < n; x++) labels_before[x] = dsu_find(&d, x);
            for (int x = n - 1; x >= 0; x--) labels_after[x] = dsu_find(&d, x);
            for (int x = 0; x < n; x++) if (labels_before[x] != labels_after[x]) ok = 0;
        }
        CHECKF(ok, "dsu trial %d", trial);
        dsu_free(&d);
    }
}

static void test_avl_invariant(void) {
    int keys[64];
    for (int i = 0; i < 64; i++) keys[i] = (i - 32) * 1000;                          /* distinct multiples of 1000 ... */
    keys[0] = -2147483647 - 1; keys[1] = 2147483647; keys[2] = -2147483647; keys[3] = 2147483646;   /* ... and the int extremes */
    keys[4] = 1; keys[5] = -1; keys[6] = 7;                                            /* all distinct from each other */
    for (int trial = 0; trial < 120; trial++) {
        seed_rng(5000 + (unsigned)trial);
        Avl t; avl_init(&t); int ok = 1; char present[64]; memset(present, 0, sizeof present);
        for (int step = 0; step < 300; step++) {
            int i = rint_(64), op = rint_(10);
            if (op < 6) { if (avl_insert(&t, keys[i], &keys[i]) != 0) ok = 0; present[i] = 1; }
            else {
                if (avl_remove(&t, keys[i]) != present[i]) ok = 0;
                present[i] = 0;
            }
            if (!avl_check(&t)) ok = 0;
            if (t.n > 0 && avl_height(&t) > 9) ok = 0;                              /* 1.44 log2(65) < 9 */
            void *v; if (avl_find(&t, keys[i], &v) != present[i]) ok = 0;
        }
        for (int i = 0; i < 64; i++) { avl_remove(&t, keys[i]); if (!avl_check(&t)) ok = 0; }
        for (int i = 0; i < 64; i++) if (avl_remove(&t, keys[i]) != 0) ok = 0;       /* removing absent keys changes nothing */
        CHECKF(ok && t.n == 0 && t.root == NULL && avl_check(&t), "avl trial %d", trial);
        avl_free(&t);
    }
}

/* ======================================================================================================== B */
typedef struct { char *p; size_t n, cap; } Str;
static void str_add(Str *s, const char *b, size_t len) {
    if (!s->p || s->n + len + 1 > s->cap) { s->cap = (s->n + len + 1) * 2 + 16; s->p = realloc(s->p, s->cap); }
    if (len) memcpy(s->p + s->n, b, len);
    s->n += len; s->p[s->n] = 0;
}
static void str_set(Str *s, const char *b, size_t len) { s->n = 0; str_add(s, b, len); }
static void str_free(Str *s) { free(s->p); s->p = NULL; s->n = s->cap = 0; }
static int write_file(const char *path, const char *buf, size_t len) {
    FILE *f = fopen(path, "wb"); if (!f) return -1;
    size_t w = len ? fwrite(buf, 1, len, f) : 0; fclose(f); return w == len ? 0 : -1;
}
static char *read_file(const char *path, size_t *len) {
    FILE *f = fopen(path, "rb"); if (!f) return NULL;
    fseek(f, 0, SEEK_END); long L = ftell(f); fseek(f, 0, SEEK_SET);
    char *b = malloc((size_t)L + 1); size_t r = L > 0 ? fread(b, 1, (size_t)L, f) : 0; fclose(f);
    b[r] = 0; *len = r; return b;
}

#define BASEFILE "build/_prop_base.csv"
#define FUZZFILE "build/_prop_fuzz.csv"
#define MAXLINES 160

/* independent expectation: non-blank lines after the first non-blank line (the header) */
static long count_rows(const char *b, size_t len) {
    size_t i = 0; long rows = 0; int header = 0;
    while (i < len) {
        size_t j = i; while (j < len && b[j] != '\n') j++;
        size_t e = j; if (e > i && b[e - 1] == '\r') e--;
        if (e > i) { if (!header) header = 1; else rows++; }
        i = j + 1;
    }
    return rows;
}
static int obj_valid(const Object *o) {
    return o->norad_id >= 1 && o->norad_id <= 999999999 && isfinite(o->epoch_jd) && o->ecc >= 0 && o->ecc < 1 && o->incl >= 0 && o->incl <= 180
        && o->mean_motion > 0 && o->mean_motion <= 20 && o->raan >= 0 && o->raan < 360 && o->argp >= 0 && o->argp < 360 && o->mean_anom >= 0 && o->mean_anom < 360
        && isfinite(o->bstar) && isfinite(o->mean_motion_dot) && o->perigee_r_km >= EARTH_RADIUS_KM && o->apogee_r_km >= o->perigee_r_km
        && memchr(o->name, 0, sizeof o->name) != NULL && memchr(o->intl_id, 0, sizeof o->intl_id) != NULL;
}
static int catalog_consistent(const Catalog *c) {
    const Object *o = cat_objects(c);
    for (int i = 0; i < cat_count(c); i++) if (!obj_valid(&o[i]) || cat_find(c, o[i].norad_id) != i) return 0;
    return 1;
}

static const char *TOK[] = {"nan", "inf", "-inf", "1e999", "-1", "0", "", "0x1p3", "  ", "1,2", "\"q\"", "99999999999999999999", "1e-999", "2026-13-45T25:61:61", "+", "-", ".", "1e308", "-0"};
#define NTOK ((int)(sizeof TOK / sizeof TOK[0]))

/* apply one random mutation to the line array L[0..*nl) (line 0 is the header) */
static void mutate_once(Str *L, int *nl) {
    int idx = rint_(20) == 0 ? 0 : 1 + rint_(*nl - 1), op = rint_(11);
    Str *s = &L[idx];
    size_t pos = s->n ? (size_t)rint_((int)s->n) : 0;
    switch (op) {
    case 0: { size_t del = (size_t)rint_(12); if (pos + del > s->n) del = s->n - pos; memmove(s->p + pos, s->p + pos + del, s->n - pos - del + 1); s->n -= del; break; }
    case 1: { char ch = (char)(32 + rint_(95)); Str t = {0}; str_add(&t, s->p, pos); str_add(&t, &ch, 1); str_add(&t, s->p + pos, s->n - pos); str_free(s); *s = t; break; }
    case 2: { static const char *sp[] = {",", "\"", "\r", " ", "\t", "\\", ";", "\xff", "\xef\xbb\xbf", "\"\"", ",,,"}; const char *x = sp[rint_(11)];
              Str t = {0}; str_add(&t, s->p, pos); str_add(&t, x, strlen(x)); str_add(&t, s->p + pos, s->n - pos); str_free(s); *s = t; break; }
    case 3: s->n = pos; s->p[pos] = 0; break;
    case 4: { char big[5200]; memset(big, rint_(2) ? 'A' : '7', sizeof big); Str t = {0}; str_add(&t, s->p, pos); str_add(&t, big, sizeof big - 1 - (size_t)rint_(3)); str_add(&t, s->p + pos, s->n - pos); str_free(s); *s = t; break; }
    case 5: if (*nl < MAXLINES) { memmove(&L[idx + 1], &L[idx], sizeof(Str) * (size_t)(*nl - idx)); L[idx + 1].p = NULL; L[idx + 1].n = L[idx + 1].cap = 0; str_add(&L[idx + 1], L[idx].p, L[idx].n); (*nl)++; } break;
    case 6: {                                                                          /* replace field number `field` by a hostile token */
        int field = rint_(17), cur = 0, found = 0; size_t a = 0, e = s->n;
        for (size_t i = 0; i <= s->n; i++)
            if (i == s->n || s->p[i] == ',') { if (cur == field) { e = i; found = 1; break; } cur++; a = i + 1; }
        if (found) { const char *x = TOK[rint_(NTOK)]; Str t = {0}; str_add(&t, s->p, a); str_add(&t, x, strlen(x)); str_add(&t, s->p + e, s->n - e); str_free(s); *s = t; }
        break; }
    case 7: { Str t = {0}; str_add(&t, s->p, s->n); for (int i = 0; i < 70; i++) str_add(&t, ",", 1); str_free(s); *s = t; break; }
    case 8: if (*nl < MAXLINES) { memmove(&L[idx + 1], &L[idx], sizeof(Str) * (size_t)(*nl - idx)); L[idx].p = NULL; L[idx].n = L[idx].cap = 0; str_add(&L[idx], "", 0); (*nl)++; } break;
    case 9: { int j = 1 + rint_(*nl - 1); Str t = L[idx]; L[idx] = L[j]; L[j] = t; break; }
    default: if (*nl > 2) { str_free(&L[idx]); memmove(&L[idx], &L[idx + 1], sizeof(Str) * (size_t)(*nl - idx - 1)); (*nl)--; } break;
    }
}


static int split_lines(const char *b, size_t len, Str *L, int maxl) {
    int nl = 0; size_t i = 0;
    while (i < len && nl < maxl) {
        size_t j = i; while (j < len && b[j] != '\n') j++;
        L[nl].p = NULL; L[nl].n = L[nl].cap = 0; str_add(&L[nl], b + i, j - i); nl++;
        i = j + 1;
    }
    return nl;
}
static void join_lines(const Str *L, int nl, const char *sep, int final_newline, Str *out) {
    str_set(out, "", 0);
    for (int i = 0; i < nl; i++) { str_add(out, L[i].p, L[i].n); if (i + 1 < nl || final_newline) str_add(out, sep, strlen(sep)); }
}

/* Random mutations of a valid file: the loader must never crash, must return 0 or -2, must keep its invariants, and must
 * account for every non-blank data line exactly once (added, replaced, duplicate, suborbital or bad). */
static void test_parser_fuzz(void) {
    Catalog base; cat_init(&base);
    CHECK(synth_generate(&base, NULL, 30, 4242, 3) == 0);
    CHECK(cat_write_csv(&base, BASEFILE) == 0);
    cat_free(&base);
    size_t blen = 0; char *bb = read_file(BASEFILE, &blen); CHECK(bb != NULL && blen > 100);
    if (!bb) return;
    Str L[MAXLINES + 8], M[MAXLINES + 8]; int nl0 = split_lines(bb, blen, L, MAXLINES);
    CHECK(nl0 == 31);
    long loaded = 0, header_err = 0, bad_rows = 0, suborb = 0, dups = 0;
    for (int trial = 0; trial < 600; trial++) {
        seed_rng(6000 + (unsigned)trial);
        for (int i = 0; i < nl0; i++) { M[i].p = NULL; M[i].n = M[i].cap = 0; str_add(&M[i], L[i].p, L[i].n); }
        int nl = nl0, nm = 1 + rint_(3);
        for (int k = 0; k < nm; k++) mutate_once(M, &nl);
        Str out = {0}; join_lines(M, nl, rint_(4) == 0 ? "\r\n" : "\n", rint_(5) != 0, &out);
        CHECK(write_file(FUZZFILE, out.p, out.n) == 0);
        long expect = count_rows(out.p, out.n);
        Catalog c; cat_init(&c); LoadReport rep; report_init(&rep);
        int rc = cat_load_csv(&c, FUZZFILE, &rep);
        CHECKF(rc == 0 || rc == -2, "fuzz trial %d rc %d", trial, rc);
        if (rc == 0) {
            long acc = rep.rows_added + rep.rows_replaced + rep.rows_dup_ignored + rep.rows_suborbital + rep.rows_bad;
            CHECKF(acc == expect, "fuzz trial %d: %ld rows accounted, %ld data lines in file", trial, acc, expect);
            CHECKF(rep.errors_total == rep.rows_bad, "fuzz trial %d errors_total %ld rows_bad %ld", trial, rep.errors_total, rep.rows_bad);
            CHECKF((long)rep.errors.n == (rep.errors_total < 100 ? rep.errors_total : 100), "fuzz trial %d error list length", trial);
            CHECKF(cat_count(&c) == rep.rows_added, "fuzz trial %d count %d added %ld", trial, cat_count(&c), rep.rows_added);
            CHECKF(catalog_consistent(&c), "fuzz trial %d object invariants", trial);
            loaded++; bad_rows += rep.rows_bad; suborb += rep.rows_suborbital; dups += rep.rows_dup_ignored + rep.rows_replaced;
            for (size_t k = 0; k < rep.errors.n; k++) { const ParseError *e = da_at(&rep.errors, k); CHECKF(e->line >= 1 && e->msg[0], "fuzz trial %d error entry %zu", trial, k); }
        } else header_err++;
        report_free(&rep); cat_free(&c); str_free(&out);
        for (int i = 0; i < nl; i++) str_free(&M[i]);
    }
    CHECKF(loaded > 300 && header_err > 0 && bad_rows > 100 && dups > 0, "fuzz mix: loaded %ld header errors %ld bad rows %ld dup/replace %ld suborbital %ld", loaded, header_err, bad_rows, dups, suborb);
    for (int i = 0; i < nl0; i++) str_free(&L[i]);
    free(bb);
}

/* Columns are found by header name, so any column order must give the same catalogue. */
static void test_parser_column_permutation(void) {
    Catalog ref; cat_init(&ref); LoadReport r0; report_init(&r0);
    CHECK(cat_load_csv(&ref, BASEFILE, &r0) == 0 && cat_count(&ref) == 30 && r0.rows_bad == 0);
    size_t blen; char *bb = read_file(BASEFILE, &blen); if (!bb) return;
    Str L[MAXLINES]; int nl = split_lines(bb, blen, L, MAXLINES);
    for (int trial = 0; trial < 30; trial++) {
        seed_rng(8000 + (unsigned)trial);
        int perm[17]; for (int i = 0; i < 17; i++) perm[i] = i;
        for (int i = 16; i > 0; i--) { int j = rint_(i + 1), t = perm[i]; perm[i] = perm[j]; perm[j] = t; }
        Str out = {0}; str_add(&out, "", 0);
        for (int li = 0; li < nl; li++) {
            const char *f[17]; size_t fl[17]; int nf = 0; size_t a = 0;
            for (size_t i = 0; i <= L[li].n && nf < 17; i++) if (i == L[li].n || L[li].p[i] == ',') { f[nf] = L[li].p + a; fl[nf] = i - a; nf++; a = i + 1; }
            if (nf != 17) { CHECKF(0, "line %d has %d fields", li, nf); break; }
            for (int k = 0; k < 17; k++) { str_add(&out, f[perm[k]], fl[perm[k]]); str_add(&out, k < 16 ? "," : "\n", 1); }
        }
        CHECK(write_file(FUZZFILE, out.p, out.n) == 0);
        Catalog c; cat_init(&c); LoadReport rep; report_init(&rep);
        int rc = cat_load_csv(&c, FUZZFILE, &rep);
        int same = rc == 0 && cat_count(&c) == cat_count(&ref) && rep.rows_bad == 0;
        for (int i = 0; same && i < cat_count(&ref); i++) {
            const Object *a1 = &cat_objects(&ref)[i]; int j = cat_find(&c, a1->norad_id);
            if (j < 0) { same = 0; break; }
            const Object *b1 = &cat_objects(&c)[j];
            same = a1->epoch_jd == b1->epoch_jd && a1->mean_motion == b1->mean_motion && a1->ecc == b1->ecc && a1->incl == b1->incl && a1->raan == b1->raan
                && a1->argp == b1->argp && a1->mean_anom == b1->mean_anom && a1->bstar == b1->bstar && a1->mean_motion_dot == b1->mean_motion_dot
                && strcmp(a1->name, b1->name) == 0 && strcmp(a1->intl_id, b1->intl_id) == 0;
        }
        CHECKF(same, "column permutation trial %d", trial);
        report_free(&rep); cat_free(&c); str_free(&out);
    }
    for (int i = 0; i < nl; i++) str_free(&L[i]);
    free(bb); cat_free(&ref); report_free(&r0);
}

/* Binary garbage: a NUL byte inside one row must reject that row only (the next row must not be swallowed);
 * a UTF-16 file (what a Windows shell redirect can produce) must give a header error, not a crash. */
static void test_parser_binary_input(void) {
    Catalog ref; cat_init(&ref); LoadReport r0; report_init(&r0);
    CHECK(cat_load_csv(&ref, BASEFILE, &r0) == 0 && cat_count(&ref) == 30);
    size_t blen; char *bb = read_file(BASEFILE, &blen); if (!bb) return;
    Str L[MAXLINES]; int nl = split_lines(bb, blen, L, MAXLINES);
    for (int victim = 3; victim <= 27; victim += 12) {                      /* data line `victim` is object victim-1 */
        Str out = {0}; str_add(&out, "", 0);
        for (int li = 0; li < nl; li++) {
            if (li == victim) { str_add(&out, L[li].p, 3); str_add(&out, "\0", 1); str_add(&out, L[li].p + 3, L[li].n - 3); }
            else str_add(&out, L[li].p, L[li].n);
            str_add(&out, "\n", 1);
        }
        CHECK(write_file(FUZZFILE, out.p, out.n) == 0);
        Catalog c; cat_init(&c); LoadReport rep; report_init(&rep);
        int rc = cat_load_csv(&c, FUZZFILE, &rep);
        long acc = rep.rows_added + rep.rows_replaced + rep.rows_dup_ignored + rep.rows_suborbital + rep.rows_bad;
        CHECKF(rc == 0 && acc == 30, "NUL in line %d: rc %d, %ld of 30 rows accounted", victim, rc, acc);
        int next_id = cat_objects(&ref)[victim].norad_id;                          /* the row AFTER the damaged one */
        CHECKF(cat_find(&c, next_id) >= 0, "NUL in line %d swallowed the following row (id %d)", victim, next_id);
        CHECKF(cat_count(&c) >= 29, "NUL in line %d: only %d objects loaded", victim, cat_count(&c));
        report_free(&rep); cat_free(&c); str_free(&out);
    }
    { /* UTF-16LE with BOM */
        Str out = {0}; str_add(&out, "\xff\xfe", 2);
        for (size_t i = 0; i < blen; i++) { str_add(&out, bb + i, 1); str_add(&out, "\0", 1); }
        CHECK(write_file(FUZZFILE, out.p, out.n) == 0);
        Catalog c; cat_init(&c); LoadReport rep; report_init(&rep);
        CHECK(cat_load_csv(&c, FUZZFILE, &rep) == -2 && cat_count(&c) == 0);
        report_free(&rep); cat_free(&c); str_free(&out);
    }
    { /* degenerate files: empty, only NULs, header only (with and without newline), 200 kB of random bytes, one 100 kB line */
        const char *hdr = L[0].p; size_t hl = L[0].n;
        struct { const char *name; Str s; int want_rc; } cases[6]; memset(cases, 0, sizeof cases);
        cases[0].name = "empty"; str_add(&cases[0].s, "", 0); cases[0].want_rc = -2;
        cases[1].name = "nuls"; { char z[100]; memset(z, 0, sizeof z); str_add(&cases[1].s, z, sizeof z); } cases[1].want_rc = -2;
        cases[2].name = "header+nl"; str_add(&cases[2].s, hdr, hl); str_add(&cases[2].s, "\n", 1); cases[2].want_rc = 0;
        cases[3].name = "header"; str_add(&cases[3].s, hdr, hl); cases[3].want_rc = 0;
        cases[4].name = "random"; seed_rng(9); for (int i = 0; i < 200000; i++) { char ch = (char)rnd(); str_add(&cases[4].s, &ch, 1); } cases[4].want_rc = 99;     /* 0 or -2 */
        cases[5].name = "longline"; str_add(&cases[5].s, hdr, hl); str_add(&cases[5].s, "\n", 1); { char big[1000]; memset(big, 'x', sizeof big); for (int i = 0; i < 100; i++) str_add(&cases[5].s, big, sizeof big); } cases[5].want_rc = 0;
        for (int k = 0; k < 6; k++) {
            CHECK(write_file(FUZZFILE, cases[k].s.p, cases[k].s.n) == 0);
            Catalog c; cat_init(&c); LoadReport rep; report_init(&rep);
            int rc = cat_load_csv(&c, FUZZFILE, &rep);
            CHECKF(cases[k].want_rc == 99 ? (rc == 0 || rc == -2) : rc == cases[k].want_rc, "degenerate file '%s' rc %d", cases[k].name, rc);
            if (k == 5) CHECK(rep.rows_bad == 1 && cat_count(&c) == 0);
            if (k == 2 || k == 3) CHECK(rep.rows_bad == 0 && cat_count(&c) == 0);
            report_free(&rep); cat_free(&c); str_free(&cases[k].s);
        }
    }
    CHECK(cat_load_csv(&ref, "build/definitely_not_here.csv", &r0) == -1);
    for (int i = 0; i < nl; i++) str_free(&L[i]);
    free(bb); cat_free(&ref); report_free(&r0);
}

/* ---- scope filters (the LEO scope of docs/scope.md rests on these) ---- */
typedef struct { double a, b; } Par;
typedef int (*Pred)(const Object *, const Par *);
static int p_epoch(const Object *o, const Par *p) { return o->epoch_jd >= p->a && o->epoch_jd <= p->b; }
static int p_mm(const Object *o, const Par *p) { return o->mean_motion >= p->a; }
static int p_alt(const Object *o, const Par *p) { return (o->perigee_r_km - EARTH_RADIUS_KM) >= p->a && (o->apogee_r_km - EARTH_RADIUS_KM) <= p->b; }
static int filter_matches(const Catalog *src, const Catalog *dst, int ret, Pred pred, Par par) {
    const Object *s = cat_objects(src), *d = cat_objects(dst); int k = 0;
    for (int i = 0; i < cat_count(src); i++) if (pred(&s[i], &par)) { if (k >= cat_count(dst) || d[k].norad_id != s[i].norad_id) return 0; k++; }
    return k == cat_count(dst) && ret == k;
}
static void check_filters_on(const Catalog *src, const char *label) {
    int n = cat_count(src); const Object *o = cat_objects(src);
    double emin = 1e30, emax = 0; for (int i = 0; i < n; i++) { if (o[i].epoch_jd < emin) emin = o[i].epoch_jd; if (o[i].epoch_jd > emax) emax = o[i].epoch_jd; }
    for (int t = 0; t < 40; t++) {
        seed_rng(10000 + (unsigned)t);
        Par pe = {runif(emin - 1, emax + 1), 0}; pe.b = pe.a + runif(0, emax - emin + 2);
        Catalog d; cat_init(&d); int r = cat_filter_epoch(src, &d, pe.a, pe.b);
        CHECKF(filter_matches(src, &d, r, p_epoch, pe), "%s epoch filter %d", label, t); cat_free(&d);
        static const double mms[] = {0, 1, 6.4, 11.25, 14.9, 15.6, 20, 25};
        Par pm = {mms[rint_(8)], 0};
        cat_init(&d); r = cat_filter_mean_motion(src, &d, pm.a);
        CHECKF(filter_matches(src, &d, r, p_mm, pm), "%s mean-motion filter %d (min %g)", label, t, pm.a); cat_free(&d);
        static const double lo[] = {-1e9, 0, 175, 300, 550, 2000}, hi[] = {100, 800, 2000, 20000, 1e9};
        Par pa = {lo[rint_(6)], hi[rint_(5)]};
        cat_init(&d); r = cat_filter_altitude(src, &d, pa.a, pa.b);
        CHECKF(filter_matches(src, &d, r, p_alt, pa), "%s altitude filter %d (%g..%g)", label, t, pa.a, pa.b);
        Catalog d2; cat_init(&d2); int r2 = cat_filter_altitude(&d, &d2, pa.a, pa.b);            /* idempotent */
        CHECKF(r2 == r && cat_count(&d2) == cat_count(&d), "%s altitude filter idempotent %d", label, t);
        cat_free(&d2); cat_free(&d);
        /* mean-motion then altitude equals altitude then mean-motion (as sets and in order) */
        Catalog x1, x2, y1, y2; cat_init(&x1); cat_init(&x2); cat_init(&y1); cat_init(&y2);
        cat_filter_mean_motion(src, &x1, pm.a); cat_filter_altitude(&x1, &x2, pa.a, pa.b);
        cat_filter_altitude(src, &y1, pa.a, pa.b); cat_filter_mean_motion(&y1, &y2, pm.a);
        int same = cat_count(&x2) == cat_count(&y2);
        for (int i = 0; same && i < cat_count(&x2); i++) same = cat_objects(&x2)[i].norad_id == cat_objects(&y2)[i].norad_id;
        CHECKF(same, "%s filters commute %d", label, t);
        cat_free(&x1); cat_free(&x2); cat_free(&y1); cat_free(&y2);
    }
    /* exact boundaries are inclusive, one ulp beyond excludes */
    for (int t = 0; t < 12 && n > 0; t++) {
        const Object *b = &o[(t * 977) % n]; Catalog d;
        double palt = b->perigee_r_km - EARTH_RADIUS_KM, aalt = b->apogee_r_km - EARTH_RADIUS_KM;
        cat_init(&d); cat_filter_altitude(src, &d, palt, 1e9); CHECKF(cat_find(&d, b->norad_id) >= 0, "%s perigee boundary kept (%d)", label, t); cat_free(&d);
        cat_init(&d); cat_filter_altitude(src, &d, nextafter(palt, 1e18), 1e9); CHECKF(cat_find(&d, b->norad_id) < 0, "%s perigee ulp excluded (%d)", label, t); cat_free(&d);
        cat_init(&d); cat_filter_altitude(src, &d, -1e9, aalt); CHECKF(cat_find(&d, b->norad_id) >= 0, "%s apogee boundary kept (%d)", label, t); cat_free(&d);
        cat_init(&d); cat_filter_altitude(src, &d, -1e9, nextafter(aalt, -1e18)); CHECKF(cat_find(&d, b->norad_id) < 0, "%s apogee ulp excluded (%d)", label, t); cat_free(&d);
        cat_init(&d); cat_filter_mean_motion(src, &d, b->mean_motion); CHECKF(cat_find(&d, b->norad_id) >= 0, "%s mean motion boundary kept (%d)", label, t); cat_free(&d);
        cat_init(&d); cat_filter_mean_motion(src, &d, nextafter(b->mean_motion, 1e18)); CHECKF(cat_find(&d, b->norad_id) < 0, "%s mean motion ulp excluded (%d)", label, t); cat_free(&d);
        cat_init(&d); cat_filter_epoch(src, &d, b->epoch_jd, b->epoch_jd); CHECKF(cat_find(&d, b->norad_id) >= 0, "%s epoch boundary kept (%d)", label, t); cat_free(&d);
    }
}
static void test_filters(void) {
    Catalog all; cat_init(&all);
    CHECK(synth_generate(&all, NULL, 700, 77, 0) == 0);
    Catalog sp; cat_init(&sp); const Object *o = cat_objects(&all);
    for (int i = 0; i < cat_count(&all); i++) { Object c = o[i]; c.epoch_jd += 0.25 * (double)(i % 17); CHECK(cat_add(&sp, &c) == 1); }
    check_filters_on(&sp, "synthetic");
    { /* the project scope itself: nothing outside 175-2000 km survives, nothing inside is lost */
      Catalog leo; cat_init(&leo); int k = cat_filter_altitude(&sp, &leo, SCOPE_LEO_MIN_ALT_KM, SCOPE_LEO_MAX_ALT_KM);
      int inside = 0; for (int i = 0; i < cat_count(&sp); i++) { const Object *x = &cat_objects(&sp)[i]; if (x->perigee_r_km - EARTH_RADIUS_KM >= 175.0 && x->apogee_r_km - EARTH_RADIUS_KM <= 2000.0) inside++; }
      CHECK(k == inside && k > 100 && k < cat_count(&sp));                                  /* synthetic mix contains MEO/GEO objects, which must go */
      for (int i = 0; i < cat_count(&leo); i++) { const Object *x = &cat_objects(&leo)[i]; CHECK(x->perigee_r_km >= EARTH_RADIUS_KM + 175.0 && x->apogee_r_km <= EARTH_RADIUS_KM + 2000.0); }
      cat_free(&leo); }
    cat_free(&sp); cat_free(&all);
    { /* real catalogue, when the data files are present (they are git-ignored) */
      static const char *files[] = {"data/raw/active.csv", "data/raw/fengyun-1c-debris.csv"};
      Catalog real; cat_init(&real); int loaded = 0;
      for (int i = 0; i < 2; i++) { LoadReport rep; report_init(&rep); if (cat_load_csv(&real, files[i], &rep) == 0) loaded++; report_free(&rep); }
      if (loaded > 0 && cat_count(&real) > 50) check_filters_on(&real, "real");
      else printf("  (real data not present: skipping real-catalogue filter checks)\n");
      cat_free(&real); }
}

/* ---- synthetic generator ---- */
static void test_synth_properties(void) {
    static const double TARGETS[9] = {0.5, 1, 2, 4, 6, 8, 12, 15, 25};
    Catalog z; cat_init(&z); CHECK(synth_generate(&z, NULL, 0, 1, 0) == -1 && synth_generate(&z, NULL, 5, 1, -1) == -1 && synth_generate(&z, NULL, 5, 1, 3) == -1 && cat_count(&z) == 0); cat_free(&z);
    for (int trial = 0; trial < 40; trial++) {
        seed_rng(11000 + (unsigned)trial);
        int n = 1 + rint_(300), pairs = rint_(n / 2 + 1); if (pairs > 20) pairs = 20;
        unsigned seed = (unsigned)rnd();
        Catalog a, b; DArray ta, tb; cat_init(&a); cat_init(&b); da_init(&ta, sizeof(PlantedPair)); da_init(&tb, sizeof(PlantedPair));
        int ra = synth_generate(&a, &ta, n, seed, pairs), rb = synth_generate(&b, &tb, n, seed, pairs);
        CHECKF(ra == 0 && rb == 0 && cat_count(&a) == n && cat_count(&b) == n && (int)ta.n == pairs && ta.n == tb.n, "synth trial %d n %d pairs %d rc %d %d", trial, n, pairs, ra, rb);
        if (ra || rb) { cat_free(&a); cat_free(&b); da_free(&ta); da_free(&tb); continue; }
        int ok = catalog_consistent(&a), same = 1;
        for (int i = 0; i < n; i++) {
            const Object *x = &cat_objects(&a)[i], *y = &cat_objects(&b)[i];
            if (x->norad_id < SYNTH_FIRST_ID || x->norad_id >= SYNTH_FIRST_ID + n) ok = 0;
            if (x->perigee_r_km - EARTH_RADIUS_KM < 200.0 - 1e-6) ok = 0;                          /* generator promises perigee >= 200 km */
            if (x->norad_id != y->norad_id || x->mean_motion != y->mean_motion || x->ecc != y->ecc || x->incl != y->incl || x->raan != y->raan || x->argp != y->argp || x->mean_anom != y->mean_anom || x->epoch_jd != y->epoch_jd) same = 0;
        }
        for (size_t k = 0; k < ta.n; k++) {
            const PlantedPair *p = da_at(&ta, k), *q = da_at(&tb, k);
            if (cat_find(&a, p->id_a) < 0 || cat_find(&a, p->id_b) < 0 || p->id_a == p->id_b) ok = 0;
            int match = 0; for (int m = 0; m < 9; m++) if (fabs(p->miss_km - TARGETS[m]) < 1e-5) match = 1;
            if (!match || !(p->tca_s > 0 && p->tca_s < 72100) || !(p->rel_speed_kms > 0)) ok = 0;
            { const Object *oa = &cat_objects(&a)[cat_find(&a, p->id_a)]; if (fabs(p->period_s - 86400.0 / oa->mean_motion) > 1e-3) ok = 0; }
            if (p->tca_s != q->tca_s || p->miss_km != q->miss_km || p->id_a != q->id_a) same = 0;
        }
        CHECKF(ok, "synth trial %d object/truth invariants", trial);
        CHECKF(same, "synth trial %d not deterministic", trial);
        if (n >= 30) {                                                                          /* another seed gives another catalogue */
            Catalog c; cat_init(&c); synth_generate(&c, NULL, n, seed + 1u, pairs);
            int diff = 0; for (int i = 0; i < n && !diff; i++) { int j = cat_find(&c, cat_objects(&a)[i].norad_id); if (j < 0 || cat_objects(&c)[j].mean_anom != cat_objects(&a)[i].mean_anom) diff = 1; }
            CHECKF(diff, "synth trial %d seed+1 gives the same catalogue", trial); cat_free(&c);
        }
        cat_free(&a); cat_free(&b); da_free(&ta); da_free(&tb);
    }
}

/* ---- Kepler solver and two-body conservation laws ---- */
static void test_kepler_and_two_body(void) {
    static const double es[] = {0, 1e-9, 0.001, 0.1, 0.5, 0.79, 0.8, 0.9, 0.97, 0.99};
    for (int ei = 0; ei < 10; ei++) {
        double e = es[ei], worst = 0, prevE = -1; int mono = 1;
        seed_rng(12000 + (unsigned)ei);
        for (int k = 0; k < 3000; k++) {
            double M = k < 6 ? (double)k * 3.141592653589793 : runif(-20, 20);
            double E = kepler_E(M, e), Mw = fmod(M, 6.283185307179586); if (Mw < 0) Mw += 6.283185307179586;
            double res = fabs(E - e * sin(E) - Mw); if (res > worst) worst = res;
        }
        for (int k = 0; k <= 1000; k++) { double M = 6.283185307179586 * (double)k / 1001.0, E = kepler_E(M, e); if (E < prevE - 1e-12) mono = 0; prevE = E; }
        CHECKF(worst < 1e-11, "Kepler residual %.3g for e = %g", worst, e);
        CHECKF(mono, "E(M) not monotone for e = %g", e);
    }
    /* two-body (flags 0): vis-viva, constant angular momentum and the speed bound used by the screening proofs */
    Catalog c; cat_init(&c); CHECK(synth_generate(&c, NULL, 200, 31, 5) == 0);
    double ref; iso_to_jd(SYNTH_EPOCH, &ref);
    int ok = 1;
    for (int i = 0; i < cat_count(&c); i++) {
        Prop p; prop_init(&p, &cat_objects(&c)[i], 0, ref);
        double vmax = conj_vmax(&p), h0 = -1;
        for (int k = 0; k < 40; k++) {
            Vec3 r, v; prop_state(&p, (double)k * 733.0, &r, &v);
            double rr = sqrt(r.x * r.x + r.y * r.y + r.z * r.z), vv2 = v.x * v.x + v.y * v.y + v.z * v.z;
            double vis = EARTH_MU_KM3S2 * (2.0 / rr - 1.0 / p.a);
            if (fabs(vv2 - vis) > 1e-9 * vis) ok = 0;
            double hx = r.y * v.z - r.z * v.y, hy = r.z * v.x - r.x * v.z, hz = r.x * v.y - r.y * v.x, h = sqrt(hx * hx + hy * hy + hz * hz);
            if (h0 < 0) h0 = h; else if (fabs(h - h0) > 1e-9 * h0) ok = 0;
            if (sqrt(vv2) > vmax) ok = 0;
            if (rr < p.a * (1 - p.e) * (1 - 1e-9) || rr > p.a * (1 + p.e) * (1 + 1e-9)) ok = 0;
        }
        Vec3 r0, r1; prop_state(&p, 0.0, &r0, NULL); prop_state(&p, 6.283185307179586 / p.n, &r1, NULL);       /* periodic */
        if (sqrt((r0.x - r1.x) * (r0.x - r1.x) + (r0.y - r1.y) * (r0.y - r1.y) + (r0.z - r1.z) * (r0.z - r1.z)) > 1e-5) ok = 0;
    }
    CHECK(ok);
    cat_free(&c);
}

/* ======================================================================================================== C */
static Prop *mk_props(const Catalog *c, int flags, double ref) {
    int n = cat_count(c);
    Prop *p = malloc(sizeof(Prop) * (size_t)(n > 0 ? n : 1));
    for (int i = 0; i < n; i++) prop_init(&p[i], &cat_objects(c)[i], flags, ref);
    return p;
}
static int events_identical(const DArray *x, const DArray *y) {               /* bit-for-bit: same pair, same tca, same miss, same speed */
    if (x->n != y->n) return 0;
    const Conjunction *a = x->data, *b = y->data;
    for (size_t i = 0; i < x->n; i++)
        if (a[i].a != b[i].a || a[i].b != b[i].b || a[i].tca != b[i].tca || a[i].miss_km != b[i].miss_km || a[i].rel_speed_kms != b[i].rel_speed_kms) return 0;
    return 1;
}
static void ev_init(DArray *d) { da_init(d, sizeof(Conjunction)); }

/* Random configurations: the fast methods must return EXACTLY the events of the brute-force baseline. */
static void test_screening_sweep(void) {
    static const int flagsets[4] = {0, PROP_J2, PROP_J2 | PROP_DRAG, PROP_SGP4};
    static const int ns[] = {0, 1, 2, 3, 5, 17, 60, 120, 200}; static const int nsw[] = {2, 3, 5, 17, 60, 120, 120, 200, 200};
    static const double wins[] = {0, 45, 300, 1200, 3600}, steps[] = {1, 7.5, 10, 30, 60, 90, 1000}, thrs[] = {0, 0.5, 5, 25, 100, 400}, scales[] = {0.3, 1, 3};
    long total_events = 0, trials_with_events = 0, trials_flags[4] = {0, 0, 0, 0}, alt_exact = 0;
    for (int trial = 0; trial < 120; trial++) {
        seed_rng(13000 + (unsigned)trial);
        int n = (trial % 3 == 0 ? ns : nsw)[rint_(9)], pairs = n >= 4 ? rint_(3) : 0, fi = rint_(4), flags = flagsets[fi];
        unsigned seed = (unsigned)rnd();
        double window = wins[rint_(5)], step = steps[rint_(7)], thr = thrs[rint_(6)], infl_amp = rint_(3) == 0 ? 6.0 : 0.0;
        double work = 0.5 * (double)n * (double)n * window / step;
        while (work > 2.0e7 && window > 45) { window /= 2; work /= 2; }          /* keep each trial around a second */
        Catalog c; cat_init(&c);
        if (n > 0) CHECK(synth_generate(&c, NULL, n, seed, pairs) == 0);
        double ref; iso_to_jd(SYNTH_EPOCH, &ref);
        int nn = cat_count(&c); trials_flags[fi]++;
        Prop *p = mk_props(&c, flags, ref);
        double *infl = NULL;
        if (infl_amp > 0 && nn > 0) { infl = malloc(sizeof(double) * (size_t)nn); for (int i = 0; i < nn; i++) infl[i] = infl_amp * (double)((i * 7919) % 101) / 100.0; }
        double *rmin = malloc(sizeof(double) * (size_t)(nn ? nn : 1)), *rmax = malloc(sizeof(double) * (size_t)(nn ? nn : 1));
        for (int i = 0; i < nn; i++) { rmin[i] = cat_objects(&c)[i].perigee_r_km; rmax[i] = cat_objects(&c)[i].apogee_r_km; }
        ScreenParams sp = {window, step, thr, infl};
        DArray base; ev_init(&base); ScreenStats sb;
        int rcb = conj_bruteforce(p, nn, &sp, &base, &sb);
        CHECKF(rcb == 0, "sweep %d brute rc %d", trial, rcb);
        CHECKF(sb.steps == (long)ceil(window / step) && sb.pair_checks == sb.steps * (long)nn * (long)(nn - 1) / 2, "sweep %d baseline step/pair count", trial);
        total_events += (long)base.n; if (base.n) trials_with_events++;
        for (int m = 0; m < 4; m++) {
            FastScreen fs; memset(&fs, 0, sizeof fs); fs.rmin_km = rmin; fs.rmax_km = rmax;
            const char *name = m == 0 ? "kd" : m == 1 ? "grid" : m == 2 ? "grid(default cell)" : "alt";
            fs.method = m == 0 ? SCREEN_KDTREE : m == 3 ? SCREEN_ALT : SCREEN_GRID;
            fs.cell_scale = m == 1 ? scales[rint_(3)] : 0.0;
            fs.alt_pad_km = flags == 0 ? (rint_(2) ? 0.0 : 5.0) : 250.0;
            DArray ev; ev_init(&ev); ScreenStats sf;
            int rc = conj_fast(p, nn, &sp, &fs, &ev, &sf);
            CHECKF(rc == 0, "sweep %d %s rc %d (n %d flags %d window %g step %g thr %g infl %g scale %g)", trial, name, rc, nn, flags, window, step, thr, infl_amp, fs.cell_scale);
            CHECKF(sf.steps == sb.steps && sf.close_at_start == sb.close_at_start, "sweep %d %s steps/close_at_start", trial, name);
            CHECKF(fs.index_pairs == sf.pair_checks && sf.pair_checks <= sb.pair_checks, "sweep %d %s index pairs %ld vs brute %ld", trial, name, sf.pair_checks, sb.pair_checks);
            if (m != 3) {                                                        /* proven-bound methods: identical output, identical candidates */
                CHECKF(events_identical(&base, &ev), "sweep %d %s: %zu events vs baseline %zu (n %d flags %d window %g step %g thr %g infl %g)", trial, name, ev.n, base.n, nn, flags, window, step, thr, infl_amp);
                CHECKF(sf.events == sb.events && sf.candidates == sb.candidates && sf.refinements == sb.refinements, "sweep %d %s counters", trial, name);
            } else {                                                             /* altitude filter: never invents events; exact when its assumption holds */
                CompareResult r; conj_compare(&base, &ev, 1e-9, &r);
                CHECKF(r.extra == 0, "sweep %d alt invented %ld events", trial, r.extra);
                if (flags == 0 || fs.alt_violation_km <= fs.alt_pad_km) { CHECKF(events_identical(&base, &ev), "sweep %d alt differs (%zu vs %zu, pad %g, violation %g)", trial, ev.n, base.n, fs.alt_pad_km, fs.alt_violation_km); alt_exact++; }
                CHECK(sf.candidates <= sb.candidates);
            }
            da_free(&ev);
        }
        da_free(&base); free(rmin); free(rmax); free(infl); free(p); cat_free(&c);
    }
    CHECKF(total_events > 150 && trials_with_events >= 15, "sweep is vacuous: %ld events in %ld trials", total_events, trials_with_events);
    CHECKF(trials_flags[0] > 5 && trials_flags[1] > 5 && trials_flags[2] > 5 && trials_flags[3] > 5 && alt_exact > 40, "sweep did not cover every model: %ld %ld %ld %ld alt-exact %ld", trials_flags[0], trials_flags[1], trials_flags[2], trials_flags[3], alt_exact);
}

typedef struct { int lo, hi; double tca, miss, rel; } IdEv;
static int cmp_idev(const void *x, const void *y) {
    const IdEv *p = x, *q = y;
    if (p->lo != q->lo) return p->lo < q->lo ? -1 : 1;
    if (p->hi != q->hi) return p->hi < q->hi ? -1 : 1;
    return (p->tca > q->tca) - (p->tca < q->tca);
}
static IdEv *to_idev(const DArray *ev, const Catalog *c) {
    IdEv *r = malloc(sizeof(IdEv) * (ev->n ? ev->n : 1)); const Conjunction *e = ev->data;
    for (size_t i = 0; i < ev->n; i++) {
        int ia = cat_objects(c)[e[i].a].norad_id, ib = cat_objects(c)[e[i].b].norad_id;
        r[i].lo = ia < ib ? ia : ib; r[i].hi = ia < ib ? ib : ia; r[i].tca = e[i].tca; r[i].miss = e[i].miss_km; r[i].rel = e[i].rel_speed_kms;
    }
    qsort(r, ev->n, sizeof(IdEv), cmp_idev);
    return r;
}

static void run_kd(const Prop *p, int n, const ScreenParams *sp, DArray *out) {
    FastScreen fs; memset(&fs, 0, sizeof fs); fs.method = SCREEN_KDTREE;
    CHECK(conj_fast(p, n, sp, &fs, out, NULL) == 0);
}

static void test_screening_metamorphic(void) {
    double ref; iso_to_jd(SYNTH_EPOCH, &ref);
    Catalog A; cat_init(&A); CHECK(synth_generate(&A, NULL, 150, 2024, 3) == 0);
    int n = cat_count(&A); Prop *pa = mk_props(&A, PROP_J2 | PROP_DRAG, ref);
    ScreenParams sp = {3600.0, 10.0, 60.0, NULL};
    DArray ea, ea2; ev_init(&ea); ev_init(&ea2);
    run_kd(pa, n, &sp, &ea); run_kd(pa, n, &sp, &ea2);
    CHECKF(ea.n > 10, "metamorphic base run has only %zu events", ea.n);
    CHECK(events_identical(&ea, &ea2));                                              /* determinism */

    /* 1. permutation invariance: inserting the objects in another order changes indices, not events */
    for (int trial = 0; trial < 3; trial++) {
        seed_rng(14000 + (unsigned)trial);
        int *perm = malloc(sizeof(int) * (size_t)n); for (int i = 0; i < n; i++) perm[i] = i;
        for (int i = n - 1; i > 0; i--) { int j = rint_(i + 1), t = perm[i]; perm[i] = perm[j]; perm[j] = t; }
        Catalog B; cat_init(&B); for (int i = 0; i < n; i++) cat_add(&B, &cat_objects(&A)[perm[i]]);
        Prop *pb = mk_props(&B, PROP_J2 | PROP_DRAG, ref);
        DArray eb, ebb; ev_init(&eb); ev_init(&ebb);
        run_kd(pb, n, &sp, &eb); CHECK(conj_bruteforce(pb, n, &sp, &ebb, NULL) == 0);
        IdEv *x = to_idev(&ea, &A), *y = to_idev(&eb, &B), *z = to_idev(&ebb, &B);
        int same = ea.n == eb.n && ea.n == ebb.n;
        for (size_t i = 0; same && i < ea.n; i++)
            same = x[i].lo == y[i].lo && x[i].hi == y[i].hi && fabs(x[i].tca - y[i].tca) < 1e-9 && fabs(x[i].miss - y[i].miss) < 1e-9
                && y[i].lo == z[i].lo && y[i].hi == z[i].hi && fabs(y[i].tca - z[i].tca) < 1e-9 && fabs(y[i].miss - z[i].miss) < 1e-9 && fabs(y[i].rel - z[i].rel) < 1e-9;
        CHECKF(same, "permutation trial %d: %zu vs %zu events", trial, ea.n, eb.n);
        free(x); free(y); free(z); free(perm); free(pb); da_free(&eb); da_free(&ebb); cat_free(&B);
    }

    /* 2. threshold monotonicity: raising the threshold only adds events */
    { double ths[4] = {5, 15, 40, 90}; DArray prev; ev_init(&prev); size_t last = 0;
      for (int k = 0; k < 4; k++) {
          ScreenParams s2 = {3600.0, 10.0, ths[k], NULL}; DArray cur; ev_init(&cur); run_kd(pa, n, &s2, &cur);
          if (k > 0) { CompareResult r; conj_compare(&cur, &prev, 1e-6, &r); CHECKF(r.extra == 0 && cur.n >= last, "threshold %g: %ld of the lower-threshold events missing", ths[k], r.extra); }
          last = cur.n; da_free(&prev); prev = cur;
      }
      CHECK(last > 20); da_free(&prev); }

    /* 3. window prefix: events of a shorter window are the early events of the longer one (also when the window is not a multiple of the step) */
    { double wins[3] = {1800.0, 1805.0, 3000.0}; DArray full; ev_init(&full); run_kd(pa, n, &sp, &full);
      for (int k = 0; k < 3; k++) {
          ScreenParams s2 = {wins[k], 10.0, 60.0, NULL}; DArray part, want; ev_init(&part); ev_init(&want); run_kd(pa, n, &s2, &part);
          for (size_t i = 0; i < full.n; i++) { const Conjunction *e = da_at(&full, i); if (e->tca <= wins[k]) da_push(&want, e); }
          CHECKF(events_identical(&part, &want), "window %g: %zu events, prefix of the full run has %zu", wins[k], part.n, want.n);
          da_free(&part); da_free(&want);
      }
      da_free(&full); }

    /* 4. step independence of the fast pipeline (bisection result does not depend on the screening step) */
    { double st[3] = {5.0, 10.0, 30.0}; DArray e[3];
      for (int k = 0; k < 3; k++) { ScreenParams s2 = {3600.0, st[k], 60.0, NULL}; ev_init(&e[k]); run_kd(pa, n, &s2, &e[k]); }
      for (int k = 1; k < 3; k++) { CompareResult r; conj_compare(&e[0], &e[k], 0.01, &r); CHECKF(r.missed == 0 && r.extra == 0, "step %g vs 5: missed %ld extra %ld", st[k], r.missed, r.extra); }
      for (int k = 0; k < 3; k++) da_free(&e[k]); }

    /* 5. conj_best_per_pair and conj_compare */
    { DArray best; ev_init(&best); CHECK(conj_best_per_pair(&ea, &best) == 0);
      const Conjunction *E = ea.data, *B = best.data; int ok = best.n > 0 && best.n <= ea.n; size_t groups = 0;
      for (size_t i = 0; i < ea.n; ) { size_t j = i; double mn = E[i].miss_km; while (j < ea.n && E[j].a == E[i].a && E[j].b == E[i].b) { if (E[j].miss_km < mn) mn = E[j].miss_km; j++; }
          if (groups >= best.n || B[groups].a != E[i].a || B[groups].b != E[i].b || B[groups].miss_km != mn) ok = 0;
          groups++; i = j; }
      CHECK(ok && groups == best.n);
      DArray half; ev_init(&half); for (size_t i = 0; i < ea.n; i += 2) da_push(&half, da_at(&ea, i));
      CompareResult r; conj_compare(&ea, &ea, 0.0, &r); CHECK(r.matched == (long)ea.n && r.missed == 0 && r.extra == 0);
      conj_compare(&ea, &half, 0.0, &r); CHECK(r.extra == 0 && r.matched == (long)half.n && r.missed == (long)(ea.n - half.n));
      DArray shifted; ev_init(&shifted); for (size_t i = 0; i < ea.n; i++) { Conjunction e = E[i]; e.tca += 1.0e6; da_push(&shifted, &e); }
      conj_compare(&ea, &shifted, 1.0, &r); CHECK(r.matched == 0 && r.missed == (long)ea.n && r.extra == (long)ea.n);
      da_free(&best); da_free(&half); da_free(&shifted); }

    da_free(&ea); da_free(&ea2); free(pa); cat_free(&A);
}

/* Screening parameters that make no sense must be refused with -4, not looped on, wrapped or cast out of range. */
static void test_screen_parameter_validation(void) {
    Catalog c; cat_init(&c); CHECK(synth_generate(&c, NULL, 20, 5, 1) == 0);
    double ref; iso_to_jd(SYNTH_EPOCH, &ref); int n = cat_count(&c);
    Prop *p = mk_props(&c, 0, ref);
    double nanv = nan(""), infv = HUGE_VAL;
    struct { double w, h, t; const char *what; } bad[] = {
        {600, 0, 10, "step 0"}, {600, -10, 10, "negative step"}, {600, nanv, 10, "NaN step"}, {600, infv, 10, "infinite step"},
        {-1, 10, 10, "negative window"}, {nanv, 10, 10, "NaN window"}, {infv, 10, 10, "infinite window"},
        {600, 10, -1, "negative threshold"}, {600, 10, nanv, "NaN threshold"}, {600, 10, infv, "infinite threshold"},
        {3600, 1e-12, 10, "more than 2e9 steps"}};
    for (unsigned i = 0; i < sizeof bad / sizeof bad[0]; i++) {
        ScreenParams sp = {bad[i].w, bad[i].h, bad[i].t, NULL}; DArray ev; ev_init(&ev);
        CHECKF(conj_bruteforce(p, n, &sp, &ev, NULL) == -4 && ev.n == 0, "brute force accepted %s", bad[i].what);
        FastScreen fs; memset(&fs, 0, sizeof fs); fs.method = SCREEN_KDTREE;
        CHECKF(conj_fast(p, n, &sp, &fs, &ev, NULL) == -4 && ev.n == 0, "fast screen accepted %s", bad[i].what);
        da_free(&ev);
    }
    { double *infl = malloc(sizeof(double) * (size_t)n); for (int i = 0; i < n; i++) infl[i] = 1.0; infl[3] = -2.0;
      ScreenParams sp = {600, 10, 10, infl}; DArray ev; ev_init(&ev); CHECK(conj_bruteforce(p, n, &sp, &ev, NULL) == -4);
      infl[3] = nan(""); CHECK(conj_bruteforce(p, n, &sp, &ev, NULL) == -4);
      FastScreen fs; memset(&fs, 0, sizeof fs); fs.method = SCREEN_GRID; CHECK(conj_fast(p, n, &sp, &fs, &ev, NULL) == -4);
      da_free(&ev); free(infl); }
    { ScreenParams sp = {600, 10, 10, NULL}; DArray ev; ev_init(&ev); FastScreen fs; memset(&fs, 0, sizeof fs); fs.method = SCREEN_GRID;
      double badscale[] = {1e-9, 0.01, 1e4, nanv, -1.0};
      for (int i = 0; i < 5; i++) { fs.cell_scale = badscale[i]; CHECKF(conj_fast(p, n, &sp, &fs, &ev, NULL) == -4, "grid cell scale %g accepted", badscale[i]); }
      fs.cell_scale = 0.25; CHECK(conj_fast(p, n, &sp, &fs, &ev, NULL) == 0); fs.cell_scale = 4.0; CHECK(conj_fast(p, n, &sp, &fs, &ev, NULL) == 0);   /* the experiment range stays valid */
      da_free(&ev); }
    { ScreenParams sp = {0, 10, 10, NULL}; DArray ev; ev_init(&ev); ScreenStats st;            /* zero-length window is valid: no steps, no events */
      CHECK(conj_bruteforce(p, n, &sp, &ev, &st) == 0 && st.steps == 0 && ev.n == 0);
      FastScreen fs; memset(&fs, 0, sizeof fs); fs.method = SCREEN_KDTREE; CHECK(conj_fast(p, n, &sp, &fs, &ev, &st) == 0 && st.steps == 0 && ev.n == 0);
      ScreenParams z = {600, 10, 0, NULL}; CHECK(conj_bruteforce(p, n, &z, &ev, &st) == 0 && ev.n == 0);     /* zero threshold: nothing is strictly closer than 0 */
      da_free(&ev); }
    free(p); cat_free(&c);
}

/* ======================================================================================================== D */
static int rand_graph(Graph *g, int n, int pct, int **eu_out, int **ev_out) {
    int cap = n * (n - 1) / 2 + 1, m = 0, *eu = malloc(sizeof(int) * (size_t)cap), *ev = malloc(sizeof(int) * (size_t)cap); double *w = malloc(sizeof(double) * (size_t)cap);
    for (int i = 0; i < n; i++) for (int j = i + 1; j < n; j++) if (rint_(100) < pct) { eu[m] = i; ev[m] = j; w[m] = (double)m; m++; }
    int rc = graph_build(g, n, m, eu, ev, w);
    free(w); *eu_out = eu; *ev_out = ev;
    return rc == 0 ? m : -1;
}
static int optimum_by_enumeration(int n, int m, const int *eu, const int *ev, const unsigned char *mov, int *unres) {
    int u = 0, best = 1 << 30;
    for (int e = 0; e < m; e++) if (mov && !mov[eu[e]] && !mov[ev[e]]) u++;
    for (unsigned mask = 0; mask < (1u << n); mask++) {
        int ok = 1, sz = 0;
        for (int v = 0; v < n; v++) if (mask >> v & 1) { if (mov && !mov[v]) { ok = 0; break; } sz++; }
        if (!ok || sz >= best) continue;
        for (int e = 0; e < m && ok; e++) { if (mov && !mov[eu[e]] && !mov[ev[e]]) continue; if (!(mask >> eu[e] & 1) && !(mask >> ev[e] & 1)) ok = 0; }
        if (ok) best = sz;
    }
    *unres = u; return best;
}
static int cover_is_valid(int m, const int *eu, const int *ev, const unsigned char *mov, const unsigned char *cov) {
    for (int e = 0; e < m; e++) {
        int clearable = !mov || mov[eu[e]] || mov[ev[e]];
        if (clearable && !cov[eu[e]] && !cov[ev[e]]) return 0;
    }
    return 1;
}
/* exact cover = subset-enumeration optimum, with the reductions on AND off; greedy valid and never below the optimum; lower bound never above it */
static void test_cover_reductions(void) {
    long compared = 0;
    for (int trial = 0; trial < 300; trial++) {
        seed_rng(15000 + (unsigned)trial);
        int n = 2 + rint_(11), pct = 10 + rint_(70), *eu, *ev; Graph g;
        int m = rand_graph(&g, n, pct, &eu, &ev); CHECK(m >= 0); if (m < 0) { free(eu); free(ev); continue; }
        unsigned char mov[16], cov[16], cg[16]; int use_mov = rint_(2);
        for (int v = 0; v < n; v++) mov[v] = (unsigned char)(rint_(100) < 70);
        const unsigned char *mp = use_mov ? mov : NULL;
        int unres_ref, opt = optimum_by_enumeration(n, m, eu, ev, mp, &unres_ref);
        for (int red = 1; red >= 0; red--) {
            vc_use_reductions = red; long nodes = 0; int complete = 0, unres = -1;
            memset(cov, 0, sizeof cov);
            int sz = vc_exact(&g, mp, 1000000L, cov, &nodes, &complete, &unres);
            int count = 0; for (int v = 0; v < n; v++) count += cov[v];
            CHECKF(complete && sz == opt && count == sz && unres == unres_ref && cover_is_valid(m, eu, ev, mp, cov),
                   "cover trial %d reductions %d: size %d (optimum %d) complete %d unresolved %d/%d", trial, red, sz, opt, complete, unres, unres_ref);
            compared++;
        }
        vc_use_reductions = 1;
        int ug = -1, gs = vc_greedy(&g, mp, cg, &ug), cnt = 0; for (int v = 0; v < n; v++) cnt += cg[v];
        CHECKF(gs == cnt && gs >= opt && ug == unres_ref && cover_is_valid(m, eu, ev, mp, cg), "greedy trial %d: %d vs optimum %d", trial, gs, opt);
        CHECKF(vc_lower_bound(&g, mp) <= opt, "lower bound above optimum, trial %d", trial);
        free(eu); free(ev); graph_free(&g);
    }
    vc_use_reductions = 1;
    CHECK(compared == 600);
}

static void test_topk_extremes(void) {
    double key[6] = {5, 1, 1, 9, 0.5, 1}; int out[8];
    CHECK(topk_smallest(key, 6, 2147483647, out) == 6 && out[0] == 4 && out[1] == 1 && out[2] == 2 && out[3] == 5 && out[4] == 0 && out[5] == 3);   /* K far above n: all, ascending */
    CHECK(topk_smallest(key, 6, 0, out) == 0 && topk_smallest(key, 6, -3, out) == 0 && topk_smallest(key, 0, 4, out) == 0);
    CHECK(topk_smallest(key, 6, 6, out) == 6 && topk_smallest(key, 6, 1, out) == 1 && out[0] == 4);
    double ties[5] = {2, 2, 2, 2, 2};
    CHECK(topk_smallest(ties, 5, 3, out) == 3 && out[0] == 0 && out[1] == 1 && out[2] == 2);                                              /* ties: lowest indices */
}

/* ---- analyze_events: end-to-end driver on hand-made, broken and extreme event files ---- */
static FILE *wf(const char *path) { FILE *f = fopen(path, "w"); CHECK(f != NULL); return f; }
static int slurp_ids(const char *path, int *ids, int cap, int col_cover) {          /* ids listed in a *_cover.csv with the "in_best_cover" flag set */
    FILE *f = fopen(path, "r"); if (!f) return -1;
    char line[256]; int n = 0; if (!fgets(line, sizeof line, f)) { fclose(f); return -1; }
    while (fgets(line, sizeof line, f)) { int id, g, b; if (sscanf(line, "%d,%d,%d", &id, &g, &b) == 3 && (col_cover ? b : g) && n < cap) ids[n++] = id; }
    fclose(f); return n;
}
static int has_id(const int *ids, int n, int id) { for (int i = 0; i < n; i++) if (ids[i] == id) return 1; return 0; }

static void test_analyze_robustness(void) {
    const char *EV = "build/_p_events.csv", *MOV = "build/_p_movable.csv", *PRE = "build/_p_out";
    FILE *f;
    /* 1. header only, and a file with garbage lines only: no objects, no crash, success */
    f = wf(EV); if (f) { fprintf(f, "norad_id_a,norad_id_b,tca_s,miss_km,rel_speed_kms\n"); fclose(f); }
    CHECK(analyze_events(EV, NULL, 5, PRE, 1000) == 0);
    f = wf(EV); if (f) { fprintf(f, "header\nnot,an,event\n1,2,3\nabc\n\n,,,,\n1,2,3.0,4.0\n99999999999999999999,2,3,4,5\n"); fclose(f); }
    CHECK(analyze_events(EV, NULL, 5, PRE, 1000) == 0);
    /* 2. missing files and bad arguments: non-zero, nothing leaked (checked by LeakSanitizer in the sanitizer build) */
    CHECK(analyze_events("build/_p_no_such_events.csv", NULL, 5, PRE, 10) != 0);
    f = wf(EV); if (f) { fprintf(f, "h\n10,11,100,3,7\n11,12,200,2,7\n"); fclose(f); }
    CHECK(analyze_events(EV, "build/_p_no_such_movable.csv", 5, PRE, 10) != 0);
    /* 3. self pairs (an object with itself) are ignored; extreme K values do not overflow */
    f = wf(EV); if (f) { fprintf(f, "h\n10,10,100,0,0\n10,11,100,3,7\n11,12,200,2,7\n"); fclose(f); }
    CHECK(analyze_events(EV, NULL, 5, PRE, 1000) == 0);
    CHECK(analyze_events(EV, NULL, 0, PRE, 1000) == 0);
    CHECK(analyze_events(EV, NULL, -4, PRE, 1000) == 0);
    CHECK(analyze_events(EV, NULL, 2147483647, PRE, 1000) == 0);
    /* 4. times a 32-bit millisecond key cannot hold: negative, 40 days, 1e15, NaN, infinity */
    f = wf(EV); if (f) { fprintf(f, "h\n1,2,-50,3,7\n1,3,3456000,2,7\n1,4,1e15,1,7\n1,5,nan,4,7\n1,6,inf,5,7\n1,7,1e300,6,7\n1,8,600,7,7\n"); fclose(f); }
    CHECK(analyze_events(EV, NULL, 3, PRE, 1000) == 0);
    /* 5. many events at one instant (the AVL key bump) */
    f = wf(EV); if (f) { fprintf(f, "h\n"); for (int i = 0; i < 400; i++) fprintf(f, "%d,%d,1234.5,%d,7\n", 1000, 2000 + i, 1 + i % 9); fclose(f); }
    CHECK(analyze_events(EV, NULL, 3, PRE, 1000) == 0);
    { FILE *t = fopen("build/_p_out_timeline.csv", "r"); int total = 0, bins = 0; char line[128];
      CHECK(t != NULL); if (t) { if (!fgets(line, sizeof line, t)) line[0] = 0; while (fgets(line, sizeof line, t)) { int mn, c; if (sscanf(line, "%d,%d", &mn, &c) == 2) { total += c; bins++; } } fclose(t); }
      CHECKF(total == 400 && bins >= 3, "timeline lost conjunctions: total %d bins %d", total, bins); }
    /* 6. movable subset: the chosen cover only contains movable objects and clears every clearable conjunction */
    {
        Catalog all; cat_init(&all); CHECK(synth_generate(&all, NULL, 40, 77, 0) == 0);
        Catalog mv; cat_init(&mv);
        for (int i = 0; i < 40; i++) if (i % 3 != 0) cat_add(&mv, &cat_objects(&all)[i]);
        CHECK(cat_write_csv(&mv, MOV) == 0);
        int eu[200], ev[200], m = 0; seed_rng(16000);
        f = wf(EV); if (f) { fprintf(f, "h\n");
            for (int i = 0; i < 40; i++) for (int j = i + 1; j < 40; j++) if (rint_(100) < 9) { eu[m] = cat_objects(&all)[i].norad_id; ev[m] = cat_objects(&all)[j].norad_id; fprintf(f, "%d,%d,%d.0,%d.5,7.0\n", eu[m], ev[m], 100 + m, 1 + m % 7); m++; if (m == 200) break; }
            fclose(f); }
        CHECK(m > 20);
        CHECK(analyze_events(EV, MOV, 5, PRE, 200000) == 0);
        int cover[64], nc = slurp_ids("build/_p_out_cover.csv", cover, 64, 1), ok = nc >= 0;
        for (int i = 0; i < nc; i++) if (cat_find(&mv, cover[i]) < 0) ok = 0;                               /* only movable objects */
        for (int e = 0; e < m; e++) {
            int mu = cat_find(&mv, eu[e]) >= 0, mvv = cat_find(&mv, ev[e]) >= 0;
            if ((mu || mvv) && !has_id(cover, nc, eu[e]) && !has_id(cover, nc, ev[e])) ok = 0;               /* clearable conjunction left uncovered */
        }
        CHECKF(ok, "movable-subset cover invalid (%d chosen, %d events)", nc, m);
        int gcover[64], ng = slurp_ids("build/_p_out_cover.csv", gcover, 64, 0); CHECK(ng >= nc);            /* greedy is never smaller than the best found */
        CHECK(analyze_events(EV, NULL, 5, PRE, 200000) == 0);                                                /* all movable */
        nc = slurp_ids("build/_p_out_cover.csv", cover, 64, 1); ok = nc > 0;
        for (int e = 0; e < m; e++) if (!has_id(cover, nc, eu[e]) && !has_id(cover, nc, ev[e])) ok = 0;
        CHECK(ok);
        cat_free(&mv); cat_free(&all);
    }
}

/* grid_build must leave a Grid that grid_free can release after every kind of failure (bad cell size, cell count out of range) */
static void test_grid_failure_safe(void) {
    Vec3 pts[3] = {{1, 2, 3}, {4, 5, 6}, {1e12, 0, 0}};
    double cells[] = {0.0, -1.0, nan(""), 1e-9};
    for (unsigned i = 0; i < sizeof cells / sizeof cells[0]; i++) {
        Grid g; memset(&g, 0xA5, sizeof g);                                      /* garbage in */
        int rc = grid_build(&g, pts, 3, cells[i]);
        CHECKF(rc == -2, "grid_build with cell %g returned %d", cells[i], rc);
        grid_free(&g); grid_free(&g);                                           /* safe, and safe twice */
    }
    Grid g; memset(&g, 0xA5, sizeof g);
    CHECK(grid_build(&g, pts, 2, 10.0) == 0); grid_free(&g); grid_free(&g);
    CHECK(grid_build(&g, pts, 3, 10.0) == -2); grid_free(&g);                    /* the third point is 1e11 cells away */
}

/* Range checks of the CSV loader at the exact boundaries: one value just inside and one just outside each limit. */
static void test_parser_boundaries(void) {
    struct { const char *ecc, *inc; int bad; } rows[] = {
        {"0", "0", 0}, {"0.05", "180", 0}, {"0.00001", "90", 0}, {"0.1", "179.9999", 0}, {"0.9999999", "10", 0},   /* the last one is in range but suborbital: not a bad row */
        {"1", "10", 1}, {"1.0000001", "10", 1}, {"-0.0000001", "10", 1}, {"2", "10", 1},
        {"0.1", "180.0001", 1}, {"0.1", "180.5", 1}, {"0.1", "181", 1}, {"0.1", "-0.0001", 1}, {"0.1", "360", 1}};
    const unsigned nrows = sizeof rows / sizeof rows[0];
    FILE *f = fopen("build/_p_bound.csv", "w"); CHECK(f != NULL); if (!f) return;
    fprintf(f, "OBJECT_NAME,OBJECT_ID,EPOCH,MEAN_MOTION,ECCENTRICITY,INCLINATION,RA_OF_ASC_NODE,ARG_OF_PERICENTER,MEAN_ANOMALY,NORAD_CAT_ID,BSTAR,MEAN_MOTION_DOT\n");
    int want_bad = 0;
    for (unsigned i = 0; i < nrows; i++) { fprintf(f, "B%u,2026-001A,2026-10-20T00:00:00,12.0,%s,%s,10,20,30,%u,0.0001,0.0\n", i, rows[i].ecc, rows[i].inc, 5000 + i); want_bad += rows[i].bad; }
    fclose(f);
    Catalog c; LoadReport rep; cat_init(&c); report_init(&rep);
    CHECK(cat_load_csv(&c, "build/_p_bound.csv", &rep) == 0);
    CHECKF(rep.rows_bad == want_bad && rep.rows_added + rep.rows_suborbital == (long)nrows - want_bad, "boundary rows: %ld bad (want %d), %ld added + %ld suborbital (want %ld)", rep.rows_bad, want_bad, rep.rows_added, rep.rows_suborbital, (long)nrows - want_bad);
    for (unsigned i = 0; i < nrows; i++) {
        int present = cat_find(&c, (int)(5000 + i)) >= 0;
        if (rows[i].bad) CHECKF(!present, "row %u (ecc %s, inc %s) was accepted", i, rows[i].ecc, rows[i].inc);
        else if (i < 4) CHECKF(present, "row %u (ecc %s, inc %s) was refused", i, rows[i].ecc, rows[i].inc);
    }
    cat_free(&c); report_free(&rep);
}

/* Hybrid refinement (refine_linear = 2): linear interpolation for crossing pairs, bisection for slow pairs. Against the bisection baseline
   it must invent nothing, lose at most 1 % of the events, and reproduce slow events bit for bit; invalid modes are refused. */
static void test_hybrid_refinement(void) {
    double ref; iso_to_jd(SYNTH_EPOCH, &ref);
    long slow_checked = 0, base_total = 0, lin_extra = 0;
    for (int trial = 0; trial < 6; trial++) {
        Catalog c; cat_init(&c); CHECK(synth_generate(&c, NULL, 160, 31000 + (unsigned)trial, 4) == 0);
        int n = cat_count(&c); Prop *p = mk_props(&c, trial % 2 ? (PROP_J2 | PROP_DRAG) : PROP_J2, ref);
        double step = trial < 3 ? 10.0 : 30.0; ScreenParams sp = {3600.0, step, 80.0, NULL};
        DArray base, hyb, lin; ev_init(&base); ev_init(&hyb); ev_init(&lin);
        CHECK(conj_bruteforce(p, n, &sp, &base, NULL) == 0);
        FastScreen fs; memset(&fs, 0, sizeof fs); fs.method = SCREEN_KDTREE; fs.refine_linear = 2;
        CHECK(conj_fast(p, n, &sp, &fs, &hyb, NULL) == 0);
        fs.refine_linear = 1; CHECK(conj_fast(p, n, &sp, &fs, &lin, NULL) == 0);
        CompareResult r, rl; conj_compare(&base, &hyb, 2.0, &r); conj_compare(&base, &lin, 2.0, &rl);
        CHECKF(r.extra == 0 && r.missed * 100 <= (long)base.n, "hybrid trial %d: %ld missed, %ld extra of %zu baseline events", trial, r.missed, r.extra, base.n);
        lin_extra += rl.extra; base_total += (long)base.n;
        const Conjunction *H = hyb.data, *B = base.data;                           /* slow hybrid events are bisection results, bit for bit */
        for (size_t i = 0; i < hyb.n; i++) if (H[i].rel_speed_kms < 1.0) {
            int found = 0; for (size_t j = 0; j < base.n && !found; j++) found = B[j].a == H[i].a && B[j].b == H[i].b && B[j].tca == H[i].tca && B[j].miss_km == H[i].miss_km;
            CHECKF(found, "hybrid slow event (%d,%d) has no bit-identical bisection twin", H[i].a, H[i].b); slow_checked++;
        }
        da_free(&base); da_free(&hyb); da_free(&lin); free(p); cat_free(&c);
    }
    CHECKF(base_total > 100, "hybrid test is vacuous: %ld baseline events", base_total);
    /* co-moving pairs (same orbit, mean anomaly shifted by 0.004-0.03 degrees = 0.5-3.7 km apart): the case that gives pure linear
       interpolation its extra events (L33). Hybrid must report exactly what bisection reports (here: none). */
    {
        Catalog g; cat_init(&g); CHECK(synth_generate(&g, NULL, 12, 77, 0) == 0);
        Catalog w; cat_init(&w); double shift[4] = {0.004, 0.01, 0.02, 0.03};
        for (int i = 0; i < 4; i++) { Object o = cat_objects(&g)[i]; cat_add(&w, &o); o.norad_id += 500000; o.mean_anom += shift[i]; cat_add(&w, &o); }
        for (int i = 4; i < 12; i++) cat_add(&w, &cat_objects(&g)[i]);
        int nw = cat_count(&w); Prop *pw = mk_props(&w, PROP_J2, ref); ScreenParams sw = {1800.0, 10.0, 5.0, NULL};
        DArray b, h, l; ev_init(&b); ev_init(&h); ev_init(&l);
        CHECK(conj_bruteforce(pw, nw, &sw, &b, NULL) == 0);
        FastScreen f2; memset(&f2, 0, sizeof f2); f2.method = SCREEN_KDTREE; f2.refine_linear = 2; CHECK(conj_fast(pw, nw, &sw, &f2, &h, NULL) == 0);
        f2.refine_linear = 1; CHECK(conj_fast(pw, nw, &sw, &f2, &l, NULL) == 0);
        CompareResult rh; conj_compare(&b, &h, 2.0, &rh);
        CHECKF(rh.missed == 0 && rh.extra == 0 && h.n == b.n, "co-moving: bisection %zu events, hybrid %zu (missed %ld, extra %ld)", b.n, h.n, rh.missed, rh.extra);
        CHECKF(l.n > b.n + 100, "co-moving: pure linear gave only %zu events (bisection %zu): the test no longer shows the L33 effect", l.n, b.n);   /* bisection has no interior minimum for a co-moving pair; linear invents hundreds */
        printf("co-moving pairs: bisection %zu events, hybrid %zu, pure linear %zu\n", b.n, h.n, l.n);
        da_free(&b); da_free(&h); da_free(&l); free(pw); cat_free(&w); cat_free(&g);
    }
    Catalog c; cat_init(&c); CHECK(synth_generate(&c, NULL, 20, 3, 1) == 0); Prop *p = mk_props(&c, 0, ref);
    ScreenParams sp = {600, 10, 10, NULL}; DArray ev; ev_init(&ev); FastScreen fs; memset(&fs, 0, sizeof fs); fs.method = SCREEN_KDTREE;
    for (int bad = -1; bad <= 3; bad += 4) { fs.refine_linear = bad; CHECKF(conj_fast(p, cat_count(&c), &sp, &fs, &ev, NULL) == -4, "refine mode %d accepted", bad); }
    da_free(&ev); free(p); cat_free(&c);
    printf("hybrid: %ld baseline events, %ld slow hybrid events checked bit for bit; pure linear had %ld extra events in the same runs\n", base_total, slow_checked, lin_extra);
}

/* ==== MAIN ==== */
/* Optional argument selects one group (so each can be timed alone): basic | screen | meta | valid | graph | analyze | hybrid */
int main(int argc, char **argv) {
    const char *only = argc > 1 ? argv[1] : "all";
    int all = !strcmp(only, "all");
    if (all || !strcmp(only, "basic")) {
        test_heap_invariant(); test_kd_invariant(); test_grid_invariant(); test_dsu_invariant(); test_avl_invariant();
        test_parser_fuzz(); test_parser_column_permutation(); test_parser_binary_input(); test_filters(); test_synth_properties(); test_kepler_and_two_body();
        test_grid_failure_safe(); test_parser_boundaries();
    }
    if (all || !strcmp(only, "screen")) test_screening_sweep();
    if (all || !strcmp(only, "meta")) test_screening_metamorphic();
    if (all || !strcmp(only, "valid")) test_screen_parameter_validation();
    if (all || !strcmp(only, "hybrid")) test_hybrid_refinement();
    if (all || !strcmp(only, "graph")) { test_cover_reductions(); test_topk_extremes(); }
    if (all || !strcmp(only, "analyze")) test_analyze_robustness();
    printf("%ld checks, %ld failed\n", checks, fails);
    return fails != 0;
}
