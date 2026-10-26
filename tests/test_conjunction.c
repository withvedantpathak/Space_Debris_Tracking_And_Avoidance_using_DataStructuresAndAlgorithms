/* Tests for conjunction.c: refinement, the rigorous brute-force baseline, planted-pair recovery,
 * independence from the step size, and agreement with an independent dense-sampling search. */
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "conjunction.h"
#include "synth.h"

static long checks = 0, fails = 0;
#define CHECK(c) do { checks++; if (!(c)) { fails++; printf("FAIL %s:%d: %s\n", __FILE__, __LINE__, #c); } } while (0)
#define NEAR(a, b, tol) do { checks++; if (!(fabs((a) - (b)) <= (tol))) { fails++; printf("FAIL %s:%d: %s=%.12g vs %.12g (tol %g)\n", __FILE__, __LINE__, #a, (double)(a), (double)(b), (double)(tol)); } } while (0)

static Prop *make_props(const Catalog *c, int flags, double ref) {
    Prop *p = malloc(sizeof(Prop) * (size_t)cat_count(c));
    for (int i = 0; i < cat_count(c); i++) prop_init(&p[i], &cat_objects(c)[i], flags, ref);
    return p;
}
static double dist_at(const Prop *a, const Prop *b, double t) {
    Vec3 p, q; prop_state(a, t, &p, NULL); prop_state(b, t, &q, NULL);
    return sqrt((p.x - q.x) * (p.x - q.x) + (p.y - q.y) * (p.y - q.y) + (p.z - q.z) * (p.z - q.z));
}

/* d(t) must be Lipschitz with constant vmax_a + vmax_b for every model variant: this is what makes the pruning rigorous */
static void test_lipschitz(void) {
    Catalog c; cat_init(&c); LoadReport r; report_init(&r);
    int rc = cat_load_csv(&c, "data/raw/fengyun-1c-debris.csv", &r); report_free(&r);
    if (rc != 0) { printf("  (skipping Lipschitz check on real data)\n"); cat_free(&c); return; }
    double ref = cat_median_epoch(&c);
    int ok = 1, n = cat_count(&c);
    for (int fl = 0; fl <= 3; fl++) {
        Prop *p = make_props(&c, fl, ref);
        for (int k = 0; k < 4000; k++) {
            int i = (k * 7919) % n, j = (k * 104729 + 13) % n;
            if (i == j) continue;
            double t = (k % 97) * 977.0, h = 1.0 + (k % 7) * 20.0;
            double L = conj_vmax(&p[i]) + conj_vmax(&p[j]);
            if (fabs(dist_at(&p[i], &p[j], t + h) - dist_at(&p[i], &p[j], t)) > L * h) ok = 0;
        }
        free(p);
    }
    CHECK(ok);
    cat_free(&c);
}

static void test_refine(void) {
    Catalog c; DArray t; cat_init(&c); da_init(&t, sizeof(PlantedPair));
    CHECK(synth_generate(&c, &t, 10, 11, 2) == 0);
    double ref; iso_to_jd(SYNTH_EPOCH, &ref);
    Prop *p = make_props(&c, 0, ref);
    const PlantedPair *pp = da_at(&t, 0);
    int ia = cat_find(&c, pp->id_a), ib = cat_find(&c, pp->id_b);
    double tm, dm, vr;
    CHECK(conj_refine(&p[ia], &p[ib], pp->tca_s - 5, pp->tca_s + 5, &tm, &dm, &vr) == 1);
    NEAR(tm, pp->tca_s, 1e-4); NEAR(dm, pp->miss_km, 1e-5); NEAR(vr, pp->rel_speed_kms, 1e-3);
    CHECK(conj_refine(&p[ia], &p[ib], pp->tca_s + 1, pp->tca_s + 11, &tm, &dm, &vr) == 0);     /* already separating */
    CHECK(conj_refine(&p[ia], &p[ib], pp->tca_s - 40, pp->tca_s - 30, &tm, &dm, &vr) == 0);    /* still approaching: no minimum inside */
    free(p); cat_free(&c); da_free(&t);
}

static int find_event(const DArray *ev, int a, int b, double t, double tol, double *miss) {
    const Conjunction *e = ev->data;
    for (size_t i = 0; i < ev->n; i++) if (e[i].a == a && e[i].b == b && fabs(e[i].tca - t) <= tol) { *miss = e[i].miss_km; return 1; }
    return 0;
}

static void test_planted(void) {
    Catalog c; DArray t, ev; cat_init(&c); da_init(&t, sizeof(PlantedPair)); da_init(&ev, sizeof(Conjunction));
    CHECK(synth_generate(&c, &t, 200, 3, 18) == 0);
    double ref; iso_to_jd(SYNTH_EPOCH, &ref);
    int n = cat_count(&c);
    Prop *p = make_props(&c, 0, ref);
    ScreenParams sp = {86400.0, 10.0, 10.0, NULL}; ScreenStats st;
    CHECK(conj_bruteforce(p, n, &sp, &ev, &st) == 0);
    CHECK(st.steps == 8640 && st.pair_checks == 8640L * n * (n - 1) / 2 && (double)st.candidates < 0.01 * (double)st.pair_checks);
    int ok = 1, hits = 0, decoys = 0, nev = 0;
    for (size_t k = 0; k < t.n; k++) {
        const PlantedPair *pp = da_at(&t, k);
        int ia = cat_find(&c, pp->id_a), ib = cat_find(&c, pp->id_b);
        int a = ia < ib ? ia : ib, b = ia < ib ? ib : ia;
        double half = 0.5 * pp->period_s;
        /* every recurrence tca + k*period/2 inside (0, window] must be an event (hits) or absent (decoys at 10 km) */
        for (long m = -200; m <= 200; m++) {
            double tc = pp->tca_s + (double)m * half, miss;
            if (tc <= 1.0 || tc > sp.window_s - 1.0) continue;
            int found = find_event(&ev, a, b, tc, 0.05, &miss);
            if (pp->miss_km < sp.threshold_km) { if (!found || fabs(miss - pp->miss_km) > 1e-3) { ok = 0; printf("  pair %zu m=%ld t=%.2f found=%d\n", k, m, tc, found); } nev++; }
            else if (found) ok = 0;
        }
        if (pp->miss_km < sp.threshold_km) hits++; else decoys++;
    }
    CHECK(ok); CHECK(hits == 12 && decoys == 6 && nev > 100);
    /* a larger threshold recovers the decoys too */
    DArray ev30; da_init(&ev30, sizeof(Conjunction));
    ScreenParams sp30 = {86400.0, 10.0, 30.0, NULL};
    CHECK(conj_bruteforce(p, n, &sp30, &ev30, NULL) == 0);
    int all = 1;
    for (size_t k = 0; k < t.n; k++) {
        const PlantedPair *pp = da_at(&t, k);
        int ia = cat_find(&c, pp->id_a), ib = cat_find(&c, pp->id_b); double miss;
        if (!find_event(&ev30, ia < ib ? ia : ib, ia < ib ? ib : ia, pp->tca_s, 0.05, &miss) || fabs(miss - pp->miss_km) > 1e-3) all = 0;
    }
    CHECK(all);
    /* CSV and per-pair reduction */
    DArray best; da_init(&best, sizeof(Conjunction));
    CHECK(conj_best_per_pair(&ev30, &best) == 0 && best.n >= t.n);
    for (size_t i = 1; i < best.n; i++) { const Conjunction *e = best.data; if (e[i].a == e[i - 1].a && e[i].b == e[i - 1].b) ok = 0; }
    CHECK(ok);
    CHECK(conj_write_csv(&ev30, cat_objects(&c), "build/t_events.csv") == 0);
    FILE *f = fopen("build/t_events.csv", "r"); char line[200]; int lines = 0; while (fgets(line, sizeof line, f)) lines++; fclose(f);
    CHECK((size_t)lines == ev30.n + 1);
    da_free(&best); da_free(&ev30); free(p); cat_free(&c); da_free(&t); da_free(&ev);
}

/* A catalogue dense enough to have many events at a 150 km threshold */
static void test_step_independence_and_dense_check(void) {
    Catalog c; cat_init(&c);
    CHECK(synth_generate(&c, NULL, 80, 21, 0) == 0);
    double ref; iso_to_jd(SYNTH_EPOCH, &ref);
    int n = cat_count(&c);
    int flags[2] = {0, 3};
    for (int fi = 0; fi < 2; fi++) {
        Prop *p = make_props(&c, flags[fi], ref);
        double steps[3] = {5.0, 10.0, 30.0};
        DArray ev[3];
        for (int s = 0; s < 3; s++) { da_init(&ev[s], sizeof(Conjunction)); ScreenParams sp = {7200.0, steps[s], 150.0, NULL}; CHECK(conj_bruteforce(p, n, &sp, &ev[s], NULL) == 0); }
        CHECK(ev[0].n > 40);
        CompareResult r;
        conj_compare(&ev[0], &ev[1], 0.01, &r); CHECK(r.missed == 0 && r.extra == 0);
        conj_compare(&ev[0], &ev[2], 0.01, &r); CHECK(r.missed == 0 && r.extra == 0);
        conj_compare(&ev[2], &ev[0], 0.01, &r); CHECK(r.missed == 0 && r.extra == 0);
        /* independent dense sampling at 0.5 s: discrete local minima below 150 km */
        long K = 14400;
        double *dprev = malloc(sizeof(double) * (size_t)n * (size_t)n), *dcur = malloc(sizeof(double) * (size_t)n * (size_t)n), *dnext = malloc(sizeof(double) * (size_t)n * (size_t)n);
        Vec3 *pos = malloc(sizeof(Vec3) * (size_t)n);
        DArray dense; da_init(&dense, sizeof(Conjunction));
        for (long k = 0; k <= K; k++) {
            prop_positions(p, n, 0.5 * (double)k, pos);
            for (int i = 0; i < n; i++) for (int j = i + 1; j < n; j++) {
                double dx = pos[i].x - pos[j].x, dy = pos[i].y - pos[j].y, dz = pos[i].z - pos[j].z;
                dnext[i * n + j] = sqrt(dx * dx + dy * dy + dz * dz);
            }
            if (k >= 2) for (int i = 0; i < n; i++) for (int j = i + 1; j < n; j++) {
                double a = dprev[i * n + j], b = dcur[i * n + j], cc = dnext[i * n + j];
                if (b < 150.0 && b < a && b <= cc) { Conjunction e = {i, j, 0.5 * (double)(k - 1), b, 0}; da_push(&dense, &e); }
            }
            double *tmp = dprev; dprev = dcur; dcur = dnext; dnext = tmp;
        }
        conj_sort(&dense);
        /* ignore events within 2 s of the window edges: the dense search cannot see their outer neighbour */
        DArray d2, b2; da_init(&d2, sizeof(Conjunction)); da_init(&b2, sizeof(Conjunction));
        for (size_t i = 0; i < dense.n; i++) { const Conjunction *e = da_at(&dense, i); if (e->tca > 2 && e->tca < 7198) da_push(&d2, e); }
        for (size_t i = 0; i < ev[0].n; i++) { const Conjunction *e = da_at(&ev[0], i); if (e->tca > 2 && e->tca < 7198) da_push(&b2, e); }
        conj_compare(&b2, &d2, 0.6, &r);
        CHECK(r.missed == 0 && r.extra == 0);
        if (r.missed || r.extra) printf("  flags %d: baseline %ld dense %ld matched %ld\n", flags[fi], r.truth, r.test, r.matched);
        da_free(&d2); da_free(&b2); da_free(&dense); free(dprev); free(dcur); free(dnext); free(pos);
        for (int s = 0; s < 3; s++) da_free(&ev[s]);
        free(p);
    }
    cat_free(&c);
}

int main(void) {
    test_lipschitz(); test_refine(); test_planted(); test_step_independence_and_dense_check();
    printf("%ld checks, %ld failed\n", checks, fails);
    return fails ? 1 : 0;
}
