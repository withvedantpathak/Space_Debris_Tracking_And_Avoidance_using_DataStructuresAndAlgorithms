/* Tests for catalog.c (time helpers, OMM CSV parser, merge) and synth.c (generator, planted pairs).
 * Run from the repository root (uses build/ for temp files). */
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "catalog.h"
#include "synth.h"

static long checks = 0, fails = 0;
#define CHECK(c) do { checks++; if (!(c)) { fails++; printf("FAIL %s:%d: %s\n", __FILE__, __LINE__, #c); } } while (0)
#define NEAR(a, b, tol) CHECK(fabs((a) - (b)) <= (tol))

static const char *TMP = "build/test_catalog_tmp.csv";
static void write_file(const char *path, const char *text) { FILE *f = fopen(path, "wb"); fputs(text, f); fclose(f); }
static char *slurp(const char *path) {
    FILE *f = fopen(path, "rb"); fseek(f, 0, SEEK_END); long n = ftell(f); fseek(f, 0, SEEK_SET);
    char *b = malloc((size_t)n + 1); size_t got = fread(b, 1, (size_t)n, f); b[got] = 0; fclose(f); return b;
}

static void test_time(void) {
    double jd;
    CHECK(iso_to_jd("2000-01-01T12:00:00", &jd) == 0); NEAR(jd, 2451545.0, 1e-9);
    CHECK(iso_to_jd("2000-01-01", &jd) == 0);          NEAR(jd, 2451544.5, 1e-9);
    CHECK(iso_to_jd("2024-02-29T00:00:00.5", &jd) == 0);
    CHECK(iso_to_jd("2023-02-29T00:00:00", &jd) == -1);          /* not a leap year */
    CHECK(iso_to_jd("2026-13-01T00:00:00", &jd) == -1);
    CHECK(iso_to_jd("2026-10-20T24:00:00", &jd) == -1);
    CHECK(iso_to_jd("hello", &jd) == -1);
    CHECK(iso_to_jd("", &jd) == -1);
    double a, b;
    iso_to_jd("2026-10-20T00:00:00", &a); iso_to_jd("2026-10-21T00:00:00", &b);
    NEAR(b - a, 1.0, 1e-9);
    char buf[40];
    jd_to_iso(a, buf, sizeof buf); CHECK(strcmp(buf, "2026-10-20T00:00:00.000000") == 0);
    iso_to_jd("2024-02-29T23:59:59.999999", &a); jd_to_iso(a, buf, sizeof buf); CHECK(strncmp(buf, "2024-02-29T23:59:5", 18) == 0 || strncmp(buf, "2024-03-01T00:00:0", 18) == 0);
    /* round trip over many dates, including century and leap boundaries */
    int ok = 1; unsigned x = 12345;
    for (int i = 0; i < 20000; i++) {
        x = x * 1664525u + 1013904223u;
        int Y = 1960 + (int)((x >> 8) % 240), M = 1 + (int)((x >> 4) % 12), D = 1 + (int)((x >> 12) % 28);
        int h = (int)((x >> 3) % 24), mi = (int)((x >> 9) % 60), s = (int)((x >> 15) % 60);
        char in[40], out[40];
        snprintf(in, sizeof in, "%04d-%02d-%02dT%02d:%02d:%02d.000000", Y, M, D, h, mi, s);
        if (iso_to_jd(in, &a)) { ok = 0; break; }
        jd_to_iso(a, out, sizeof out);
        /* JD in a double has ~20 us resolution: compare to the second */
        if (strncmp(in, out, 19) != 0) {
            double a2; iso_to_jd(out, &a2);
            if (fabs(a2 - a) > 1e-9) { ok = 0; printf("  roundtrip %s -> %s\n", in, out); break; }
        }
    }
    CHECK(ok);
}

static void test_parser(void) {
    /* physical line numbers: 1 header ... */
    const char *csv =
        "\xEF\xBB\xBFNORAD_CAT_ID,OBJECT_NAME,EPOCH,MEAN_MOTION,ECCENTRICITY,INCLINATION,RA_OF_ASC_NODE,ARG_OF_PERICENTER,MEAN_ANOMALY,OBJECT_ID,BSTAR\r\n" /* 1 */
        "25544,ISS (ZARYA),2026-10-20T10:00:00.000000,15.5,.0006,51.64,100.0,50.0,10.0,1998-067A,.0001\r\n"           /* 2 ok */
        "100,\"FOO, BAR\",2026-10-20T10:00:00,14.0,0.01,98.0,370,-20,30,2000-001A,0\r\n"                                /* 3 ok, quoted comma, angles wrap */
        "\r\n"                                                                                                           /* 4 blank */
        "101,BADECC,2026-10-20T10:00:00,14.0,1.2,98.0,10,20,30,2000-002A,0\r\n"                                         /* 5 bad */
        "102,SHORT,2026-10-20T10:00:00,14.0\r\n"                                                                        /* 6 bad: too few */
        "103,NOTNUM,2026-10-20T10:00:00,abc,0.01,98.0,10,20,30,2000-003A,0\r\n"                                         /* 7 bad */
        "25544,ISS (ZARYA),2026-10-21T10:00:00.000000,15.6,.0006,51.64,100.0,50.0,10.0,1998-067A,.0001\r\n"           /* 8 newer duplicate -> replace */
        "100,FOO OLD,2026-10-19T10:00:00,14.1,0.01,98.0,10,20,30,2000-001A,0\r\n"                                       /* 9 older duplicate -> ignored */
        "104,UNDERGROUND,2026-10-20T10:00:00,17.5,0,50,10,20,30,2000-004A,0\r\n"                                        /* 10 suborbital */
        "105,BADDATE,2026-13-45T00:00:00,14.0,0.01,98.0,10,20,30,2000-005A,0\r\n"                                       /* 11 bad */
        "abc,BADID,2026-10-20T10:00:00,14.0,0.01,98.0,10,20,30,2000-006A,0\r\n"                                         /* 12 bad */
        "106,LASTONE,2026-10-20T10:00:00,2.0,0.7,63.4,10,270,30,2000-007A,0";                                           /* 13 ok, no newline at EOF */
    write_file(TMP, csv);
    Catalog c; CHECK(cat_init(&c) == 0);
    LoadReport r; report_init(&r);
    CHECK(cat_load_csv(&c, TMP, &r) == 0);
    CHECK(r.rows_added == 3 && r.rows_replaced == 1 && r.rows_dup_ignored == 1 && r.rows_suborbital == 1 && r.rows_bad == 5);
    CHECK(r.lines == 13 && r.errors_total == 5 && r.errors.n == 5);
    int expect_lines[5] = {5, 6, 7, 11, 12};
    for (int i = 0; i < 5 && i < (int)r.errors.n; i++) CHECK(((ParseError *)da_at(&r.errors, (size_t)i))->line == expect_lines[i]);
    CHECK(cat_count(&c) == 3);
    int i25 = cat_find(&c, 25544), i100 = cat_find(&c, 100), i106 = cat_find(&c, 106);
    CHECK(i25 >= 0 && i100 >= 0 && i106 >= 0 && cat_find(&c, 104) < 0 && cat_find(&c, 101) < 0 && cat_find(&c, 7) < 0);
    const Object *o = cat_objects(&c);
    NEAR(o[i25].mean_motion, 15.6, 1e-12);                       /* newer element set kept */
    CHECK(strcmp(o[i100].name, "FOO, BAR") == 0);                /* older duplicate ignored; quoted comma kept */
    NEAR(o[i100].raan, 10.0, 1e-9); NEAR(o[i100].argp, 340.0, 1e-9);   /* 370 -> 10, -20 -> 340 */
    NEAR(o[i25].sma_km, 6378.137 + 0, 1e9);                      /* sanity only */
    /* derived values for ISS-like orbit: n=15.6 rev/day -> a ~ 6.7e3 km */
    CHECK(o[i25].sma_km > 6700 && o[i25].sma_km < 6800 && o[i25].perigee_r_km < o[i25].sma_km && o[i25].apogee_r_km > o[i25].sma_km);
    CHECK(o[i106].apogee_r_km > 40000);                          /* Molniya-like */
    report_free(&r);

    /* header problems and missing file */
    write_file(TMP, "OBJECT_NAME,EPOCH,MEAN_MOTION\n");
    report_init(&r); CHECK(cat_load_csv(&c, TMP, &r) == -2 && r.errors_total >= 1); report_free(&r);
    write_file(TMP, "");
    report_init(&r); CHECK(cat_load_csv(&c, TMP, &r) == -2); report_free(&r);
    report_init(&r); CHECK(cat_load_csv(&c, "build/does_not_exist.csv", &r) == -1); report_free(&r);
    cat_free(&c);
}

static void test_mmdot_and_filter(void) {
    write_file(TMP,
        "OBJECT_NAME,OBJECT_ID,EPOCH,MEAN_MOTION,ECCENTRICITY,INCLINATION,RA_OF_ASC_NODE,ARG_OF_PERICENTER,MEAN_ANOMALY,NORAD_CAT_ID,BSTAR,MEAN_MOTION_DOT\n"
        "A,2020-001A,2026-10-18T00:00:00,15.0,0.001,50,1,2,3,201,.0001,.00012345\n"
        "B,2020-001B,2026-10-20T00:00:00,15.1,0.001,50,1,2,3,202,.0001,-.00000050\n"
        "C,2020-001C,2026-10-22T00:00:00,15.2,0.001,50,1,2,3,203,.0001,0\n");
    Catalog c; Catalog f; LoadReport r; cat_init(&c); cat_init(&f); report_init(&r);
    CHECK(cat_load_csv(&c, TMP, &r) == 0 && r.rows_bad == 0 && cat_count(&c) == 3);
    NEAR(cat_objects(&c)[cat_find(&c, 201)].mean_motion_dot, 0.00012345, 1e-12);
    NEAR(cat_objects(&c)[cat_find(&c, 202)].mean_motion_dot, -5.0e-7, 1e-15);
    double lo, hi; iso_to_jd("2026-10-19T00:00:00", &lo); iso_to_jd("2026-10-21T00:00:00", &hi);
    CHECK(cat_filter_epoch(&c, &f, lo, hi) == 1 && cat_find(&f, 202) == 0 && cat_find(&f, 201) < 0);
    double med; iso_to_jd("2026-10-20T00:00:00", &med);
    NEAR(cat_median_epoch(&c), med, 1e-9);
    Catalog empty; cat_init(&empty); CHECK(cat_median_epoch(&empty) == 0.0); cat_free(&empty);
    /* write -> read keeps the field */
    CHECK(cat_write_csv(&c, "build/t_mm.csv") == 0);
    Catalog c2; LoadReport r2; cat_init(&c2); report_init(&r2);
    CHECK(cat_load_csv(&c2, "build/t_mm.csv", &r2) == 0 && cat_count(&c2) == 3);
    NEAR(cat_objects(&c2)[cat_find(&c2, 201)].mean_motion_dot, 0.00012345, 1e-9);
    report_free(&r); report_free(&r2); cat_free(&c); cat_free(&f); cat_free(&c2);
}

static void test_real_data(void) {
    FILE *f = fopen("data/raw/active.csv", "r");
    if (!f) { printf("  (skipping real-data check: data/raw/active.csv absent)\n"); return; }
    fclose(f);
    Catalog c; cat_init(&c); LoadReport r; report_init(&r);
    CHECK(cat_load_csv(&c, "data/raw/active.csv", &r) == 0);
    CHECK(r.rows_bad == 0 && r.rows_added > 10000 && (long)cat_count(&c) == r.rows_added);
    const Object *o = cat_objects(&c); int ok = 1;
    for (int i = 0; i < cat_count(&c); i++) {
        if (!(o[i].sma_km > 6400 && o[i].perigee_r_km >= EARTH_RADIUS_KM) || cat_find(&c, o[i].norad_id) != i) ok = 0;
    }
    CHECK(ok);
    report_free(&r); cat_free(&c);
}

/* ---- independent circular-orbit propagator for verifying planted pairs ---- */
static void pos(const Object *o, double t, double p[3]) {
    double n = o->mean_motion * TWO_PI_ / 86400.0, d2r = TWO_PI_ / 360.0;
    double u = o->mean_anom * d2r + n * t, O = o->raan * d2r, i = o->incl * d2r, a = o->sma_km;
    p[0] = a * (cos(O) * cos(u) - sin(O) * cos(i) * sin(u));
    p[1] = a * (sin(O) * cos(u) + cos(O) * cos(i) * sin(u));
    p[2] = a * sin(i) * sin(u);
}
static double dist(const Object *a, const Object *b, double t) {
    double p[3], q[3]; pos(a, t, p); pos(b, t, q);
    return sqrt((p[0] - q[0]) * (p[0] - q[0]) + (p[1] - q[1]) * (p[1] - q[1]) + (p[2] - q[2]) * (p[2] - q[2]));
}

static void test_synth(void) {
    Catalog c1, c2, c3; DArray t1, t2;
    da_init(&t1, sizeof(PlantedPair)); da_init(&t2, sizeof(PlantedPair));
    cat_init(&c1); cat_init(&c2); cat_init(&c3);
    CHECK(synth_generate(&c1, &t1, 3000, 42, 30) == 0);
    CHECK(synth_generate(&c2, &t2, 3000, 42, 30) == 0);
    CHECK(synth_generate(&c3, NULL, 3000, 43, 30) == 0);
    CHECK(cat_count(&c1) == 3000 && t1.n == 30);
    CHECK(cat_write_csv(&c1, "build/t_a.csv") == 0 && cat_write_csv(&c2, "build/t_b.csv") == 0 && cat_write_csv(&c3, "build/t_c.csv") == 0);
    char *a = slurp("build/t_a.csv"), *b = slurp("build/t_b.csv"), *d = slurp("build/t_c.csv");
    CHECK(strcmp(a, b) == 0);                                        /* same seed -> identical file */
    CHECK(strcmp(a, d) != 0);                                        /* different seed -> different file */
    free(a); free(b); free(d);
    CHECK(memcmp(t1.data, t2.data, t1.n * sizeof(PlantedPair)) == 0);

    /* defensive fallbacks: all-zero weights fall through to the last index; a single attempt per pair often gives up */
    { const double zero[3] = {0, 0, 0}, one[3] = {0, 1, 0};
      CHECK(synth_test_pick(zero, 3, 1) == 2);
      CHECK(synth_test_pick(one, 3, 1) == 1);
      int gave_up = 0, placed = 0;
      for (unsigned sd = 1; sd <= 200; sd++) { if (synth_test_plant(sd, (int)(sd % 9), 1) == 0) placed++; else gave_up++; }
      CHECK(gave_up > 0 && placed > 0);
      CHECK(synth_test_plant(7, 3, 1000) == 0); }

    /* argument checks */
    Catalog bad; cat_init(&bad);
    CHECK(synth_generate(&bad, NULL, 0, 1, 0) == -1 && synth_generate(&bad, NULL, 10, 1, 6) == -1 && synth_generate(&bad, NULL, 10, 1, -1) == -1);
    cat_free(&bad);

    /* reload through the real parser: same objects, valid ranges, ids unique */
    Catalog rl; cat_init(&rl); LoadReport r; report_init(&r);
    CHECK(cat_load_csv(&rl, "build/t_a.csv", &r) == 0 && r.rows_bad == 0 && r.rows_added == 3000 && r.rows_dup_ignored == 0);
    CHECK(cat_count(&rl) == 3000);
    const Object *o = cat_objects(&rl); int ok = 1, nonleo = 0;
    for (int i = 0; i < 3000; i++) {
        if (o[i].ecc < 0 || o[i].ecc >= 1 || o[i].incl < 0 || o[i].incl > 180 || o[i].perigee_r_km < EARTH_RADIUS_KM + 190) ok = 0;
        if (o[i].norad_id < SYNTH_FIRST_ID || o[i].norad_id >= SYNTH_FIRST_ID + 3000) ok = 0;
        if (o[i].sma_km - EARTH_RADIUS_KM >= 2000) nonleo++;
        const Object *m = &cat_objects(&c1)[cat_find(&c1, o[i].norad_id)];
        if (fabs(m->mean_motion - o[i].mean_motion) > 1e-11 || fabs(m->incl - o[i].incl) > 1e-9) ok = 0;
    }
    CHECK(ok);
    CHECK(nonleo > 40 && nonleo < 250);                              /* expected ~4.2% of 3000 = 126 */
    report_free(&r);

    /* planted pairs: verify with an independent propagator built from the RELOADED elements */
    int pok = 1, nhit = 0, ndecoy = 0;
    for (size_t k = 0; k < t1.n; k++) {
        const PlantedPair *p = (const PlantedPair *)da_at(&t1, k);
        int ia = cat_find(&rl, p->id_a), ib = cat_find(&rl, p->id_b);
        if (ia < 0 || ib < 0) { pok = 0; continue; }
        const Object *A = &o[ia], *B = &o[ib];
        if (A->ecc != 0 || B->ecc != 0) pok = 0;
        if (fabs(dist(A, B, p->tca_s) - p->miss_km) > 1e-4) { pok = 0; printf("  pair %zu: dist at tca %.6f vs %.6f\n", k, dist(A, B, p->tca_s), p->miss_km); }
        double best = 1e30;
        for (double t = p->tca_s - 30; t <= p->tca_s + 30; t += 0.02) { double dd = dist(A, B, t); if (dd < best) best = dd; }
        if (best < p->miss_km - 1e-4) { pok = 0; printf("  pair %zu: scan found %.6f < truth %.6f\n", k, best, p->miss_km); }
        if (p->miss_km < 0.4999 || p->miss_km > 25.0001 || p->tca_s < 1000 || p->tca_s > 73000) pok = 0;
        if (p->miss_km < 10) nhit++; else ndecoy++;
        if (p->rel_speed_kms < 1 || p->rel_speed_kms > 15.5) pok = 0;
    }
    CHECK(pok);
    CHECK(nhit == 21 && ndecoy == 9);                                /* targets cycle through 9 values, 6 of them below 10 km */
    cat_free(&rl); cat_free(&c1); cat_free(&c2); cat_free(&c3); da_free(&t1); da_free(&t2);

    /* truth file */
    DArray t; da_init(&t, sizeof(PlantedPair)); Catalog c; cat_init(&c);
    CHECK(synth_generate(&c, &t, 100, 7, 3) == 0 && synth_write_truth(&t, "build/t_truth.csv") == 0);
    char *tt = slurp("build/t_truth.csv");
    CHECK(strncmp(tt, "id_a,id_b,tca_s_after_epoch,miss_km,rel_speed_kms,period_s\n", 58) == 0);
    int lines = 0; for (char *s = tt; *s; s++) if (*s == '\n') lines++;
    CHECK(lines == 4);
    free(tt); da_free(&t); cat_free(&c);
}

int main(void) {
    test_time(); test_parser(); test_mmdot_and_filter(); test_real_data(); test_synth();
    printf("%ld checks, %ld failed\n", checks, fails);
    return fails ? 1 : 0;
}
