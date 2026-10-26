/* Tests for screen.c: every fast method must reproduce the brute-force baseline exactly. */
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "screen.h"
#include "synth.h"

static long checks = 0, fails = 0;
#define CHECK(c) do { checks++; if (!(c)) { fails++; printf("FAIL %s:%d: %s\n", __FILE__, __LINE__, #c); } } while (0)

static Prop *make_props(const Catalog *c, int flags, double ref) {
    Prop *p = malloc(sizeof(Prop) * (size_t)cat_count(c));
    for (int i = 0; i < cat_count(c); i++) prop_init(&p[i], &cat_objects(c)[i], flags, ref);
    return p;
}

/* Compare one method against the baseline on a catalogue. Returns the baseline event count. */
static long run_case(int n_obj, unsigned seed, int pairs, int flags, double window, double step, double thr, double ev_pad_check, double infl_amp) {
    Catalog c; cat_init(&c);
    CHECK(synth_generate(&c, NULL, n_obj, seed, pairs) == 0);
    double ref; iso_to_jd(SYNTH_EPOCH, &ref);
    int n = cat_count(&c);
    Prop *p = make_props(&c, flags, ref);
    double *infl = NULL;
    if (infl_amp > 0.0) { infl = malloc(sizeof(double) * (size_t)n); for (int i = 0; i < n; i++) infl[i] = infl_amp * (double)((i * 7919) % 101) / 100.0; }
    ScreenParams sp = {window, step, thr, infl};
    DArray base; da_init(&base, sizeof(Conjunction)); ScreenStats sb;
    CHECK(conj_bruteforce(p, n, &sp, &base, &sb) == 0);
    double *rmin = malloc(sizeof(double) * (size_t)n), *rmax = malloc(sizeof(double) * (size_t)n);
    for (int i = 0; i < n; i++) { rmin[i] = cat_objects(&c)[i].perigee_r_km; rmax[i] = cat_objects(&c)[i].apogee_r_km; }
    ScreenMethod ms[5] = {SCREEN_KDTREE, SCREEN_GRID, SCREEN_ALT, SCREEN_GRID, SCREEN_GRID};
    double scale[5] = {1, 1, 1, 0.4, 2.5};
    for (int m = 0; m < 5; m++) {
        DArray ev; da_init(&ev, sizeof(Conjunction)); ScreenStats sf;
        FastScreen fs; memset(&fs, 0, sizeof fs); fs.method = ms[m]; fs.rmin_km = rmin; fs.rmax_km = rmax; fs.alt_pad_km = 30.0; fs.cell_scale = scale[m];
        CHECK(conj_fast(p, n, &sp, &fs, &ev, &sf) == 0);
        CompareResult r; conj_compare(&base, &ev, 1e-9, &r);
        CHECK(r.missed == 0 && r.extra == 0 && (long)ev.n == (long)base.n);
        CHECK(sf.events == sb.events);
        if (ms[m] == SCREEN_ALT) CHECK(sf.candidates <= sb.candidates);        /* altitude filter may drop candidates that never become events */
        else CHECK(sf.candidates == sb.candidates);                            /* proven-bound methods: exactly the baseline's candidates */
        CHECK(sf.steps == sb.steps);
        CHECK(sf.close_at_start == sb.close_at_start);
        if (ms[m] != SCREEN_ALT) CHECK(sf.pair_checks < sb.pair_checks);
        if (ms[m] == SCREEN_ALT) CHECK(fs.alt_violation_km <= ev_pad_check);
        da_free(&ev);
    }
    long nev = (long)base.n;
    da_free(&base); free(infl); free(rmin); free(rmax); free(p); cat_free(&c);
    return nev;
}

int main(void) {
    long e1 = run_case(300, 5, 20, 0, 21600.0, 10.0, 10.0, 1e-6, 0.0);       /* two-body, planted pairs */
    CHECK(e1 > 20);
    long e2 = run_case(250, 9, 10, PROP_J2 | PROP_DRAG, 7200.0, 10.0, 50.0, 30.0, 0.0);   /* J2+drag, larger threshold: many events */
    CHECK(e2 > 50);
    long e3 = run_case(200, 13, 0, PROP_SGP4, 7200.0, 10.0, 25.0, 30.0, 0.0);   /* SGP4 mode: GEO fillers are parked and must be ignored */
    CHECK(e3 >= 0);
    long e4 = run_case(120, 17, 0, 0, 3600.0, 30.0, 400.0, 1e-6, 0.0);           /* big threshold and coarse step: lots of candidates */
    CHECK(e4 > 100);
    long e5 = run_case(220, 23, 12, PROP_J2 | PROP_DRAG, 7200.0, 10.0, 40.0, 30.0, 8.0);   /* per-object distance inflation: fast == brute */
    CHECK(e5 > 50);
    { /* inflation only adds events: every event of the plain run is present in the inflated run */
      Catalog c; cat_init(&c); synth_generate(&c, NULL, 150, 31, 0); double ref; iso_to_jd(SYNTH_EPOCH, &ref); int n = cat_count(&c);
      Prop *p = make_props(&c, 0, ref); double *infl = malloc(sizeof(double) * (size_t)n); for (int i = 0; i < n; i++) infl[i] = 5.0 * (double)(i % 3);
      ScreenParams s0 = {3600.0, 10.0, 15.0, NULL}, s1 = {3600.0, 10.0, 15.0, infl}; DArray a, b; da_init(&a, sizeof(Conjunction)); da_init(&b, sizeof(Conjunction));
      CHECK(conj_bruteforce(p, n, &s0, &a, NULL) == 0 && conj_bruteforce(p, n, &s1, &b, NULL) == 0);
      CompareResult r; conj_compare(&b, &a, 1e-9, &r); CHECK(r.extra == 0 && b.n >= a.n);            /* truth = inflated run, test = plain run */
      da_free(&a); da_free(&b); free(infl); free(p); cat_free(&c); }
    { /* linear interpolation (Gomez et al.) vs bisection: with a short step the two agree closely */
      Catalog c; cat_init(&c); synth_generate(&c, NULL, 200, 5, 20); double ref; iso_to_jd(SYNTH_EPOCH, &ref); int n = cat_count(&c);
      Prop *p = make_props(&c, 0, ref); ScreenParams sp = {7200.0, 1.0, 10.0, NULL}; DArray a, b; da_init(&a, sizeof(Conjunction)); da_init(&b, sizeof(Conjunction));
      FastScreen f1, f2; memset(&f1, 0, sizeof f1); memset(&f2, 0, sizeof f2); f1.method = f2.method = SCREEN_KDTREE; f2.refine_linear = 1;
      CHECK(conj_fast(p, n, &sp, &f1, &a, NULL) == 0 && conj_fast(p, n, &sp, &f2, &b, NULL) == 0);
      CompareResult r; conj_compare(&a, &b, 1.0, &r); CHECK(a.n > 5 && r.matched >= (long)a.n - 1);
      const Conjunction *A = a.data, *B = b.data; double worst = 0; long pairs = 0;
      for (size_t i = 0; i < a.n; i++) for (size_t j = 0; j < b.n; j++) if (A[i].a == B[j].a && A[i].b == B[j].b && fabs(A[i].tca - B[j].tca) < 1.0) { double d = fabs(A[i].miss_km - B[j].miss_km); if (d > worst) worst = d; pairs++; break; }
      CHECK(pairs > 5 && worst < 0.01);
      da_free(&a); da_free(&b); free(p); cat_free(&c); }
    /* bad option */
    { Catalog c; cat_init(&c); synth_generate(&c, NULL, 10, 1, 0); double ref; iso_to_jd(SYNTH_EPOCH, &ref);
      Prop *p = make_props(&c, 0, ref); ScreenParams sp = {100, 10, 10, NULL}; DArray ev; da_init(&ev, sizeof(Conjunction)); FastScreen fs; memset(&fs, 0, sizeof fs);
      fs.method = SCREEN_ALT; CHECK(conj_fast(p, cat_count(&c), &sp, &fs, &ev, NULL) == -4);       /* no radius arrays */
      fs.method = (ScreenMethod)9; CHECK(conj_fast(p, cat_count(&c), &sp, &fs, &ev, NULL) == -4);
      da_free(&ev); free(p); cat_free(&c); }
    printf("%ld checks, %ld failed (events %ld %ld %ld %ld %ld)\n", checks, fails, e1, e2, e3, e4, e5);
    return fails != 0;
}
