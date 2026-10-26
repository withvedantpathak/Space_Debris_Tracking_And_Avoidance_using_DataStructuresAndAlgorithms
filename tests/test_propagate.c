/* Tests for propagate.c against analytic cases, internal consistency and the C2 planted pairs. */
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "propagate.h"
#include "synth.h"

static long checks = 0, fails = 0;
#define CHECK(c) do { checks++; if (!(c)) { fails++; printf("FAIL %s:%d: %s\n", __FILE__, __LINE__, #c); } } while (0)
#define NEAR(a, b, tol) do { checks++; if (!(fabs((a) - (b)) <= (tol))) { fails++; printf("FAIL %s:%d: %s=%.12g vs %.12g (tol %g)\n", __FILE__, __LINE__, #a, (double)(a), (double)(b), (double)(tol)); } } while (0)

static double nrm(Vec3 a) { return sqrt(a.x * a.x + a.y * a.y + a.z * a.z); }
static Vec3 cross(Vec3 a, Vec3 b) { Vec3 c = {a.y * b.z - a.z * b.y, a.z * b.x - a.x * b.z, a.x * b.y - a.y * b.x}; return c; }
static double ddist(Vec3 a, Vec3 b) { Vec3 d = {a.x - b.x, a.y - b.y, a.z - b.z}; return nrm(d); }

static Object make(double a_km, double e, double inc, double raan, double argp, double ma) {
    Object o; memset(&o, 0, sizeof o);
    o.norad_id = 1; iso_to_jd("2026-10-20T00:00:00", &o.epoch_jd);
    o.mean_motion = sqrt(EARTH_MU_KM3S2 / (a_km * a_km * a_km)) * 86400.0 / TWO_PI_;
    o.ecc = e; o.incl = inc; o.raan = raan; o.argp = argp; o.mean_anom = ma;
    obj_derive(&o);
    return o;
}

static void test_kepler(void) {
    double worst = 0;
    for (int ie = 0; ie <= 99; ie++) {
        double e = ie * 0.0099 + (ie == 99 ? 0.0 : 0.0);
        for (int im = -50; im <= 800; im++) {
            double M = im * 0.01;
            double E = kepler_E(M, e), res = E - e * sin(E) - fmod(M, TWO_PI_);
            double Mm = fmod(M, TWO_PI_); if (Mm < 0) Mm += TWO_PI_;
            res = E - e * sin(E) - Mm;
            if (fabs(res) > worst) worst = fabs(res);
        }
    }
    CHECK(worst < 1e-12);
    NEAR(kepler_E(0.0, 0.5), 0.0, 1e-14);
    NEAR(kepler_E(3.141592653589793, 0.9), 3.141592653589793, 1e-12);
    NEAR(kepler_E(1.234, 0.0), 1.234, 1e-14);                     /* circular: E = M */
}

static void test_circular(void) {
    Object o = make(7000.0, 0.0, 51.6, 40.0, 0.0, 10.0);
    Prop p; prop_init(&p, &o, 0, o.epoch_jd);
    double period = 86400.0 / o.mean_motion;
    NEAR(period, TWO_PI_ * sqrt(7000.0 * 7000.0 * 7000.0 / EARTH_MU_KM3S2), 1e-6);   /* Kepler's third law */
    Vec3 r0, v0, r1, v1, r2;
    prop_state(&p, 0, &r0, &v0);
    int ok = 1; double zmax = 0;
    for (int i = 0; i <= 200; i++) {
        double t = period * i / 200.0;
        prop_state(&p, t, &r1, &v1);
        if (fabs(nrm(r1) - 7000.0) > 1e-8) ok = 0;
        if (fabs(nrm(v1) - sqrt(EARTH_MU_KM3S2 / 7000.0)) > 1e-10) ok = 0;
        if (fabs(r1.x * v1.x + r1.y * v1.y + r1.z * v1.z) > 1e-6) ok = 0;   /* r perpendicular to v */
        if (fabs(r1.z) > zmax) zmax = fabs(r1.z);
    }
    CHECK(ok);
    NEAR(zmax, 7000.0 * sin(51.6 * TWO_PI_ / 360.0), 1.0);            /* max height = a sin i (sampled) */
    prop_state(&p, period, &r2, NULL);
    NEAR(ddist(r0, r2), 0.0, 1e-6);                                    /* back to the start after one period */
    prop_state(&p, 5 * period, &r2, NULL);
    NEAR(ddist(r0, r2), 0.0, 1e-5);
    /* t=0 position matches the textbook formula for u = M */
    double u = 10.0 * TWO_PI_ / 360.0, O = 40.0 * TWO_PI_ / 360.0, i = 51.6 * TWO_PI_ / 360.0;
    NEAR(r0.x, 7000.0 * (cos(O) * cos(u) - sin(O) * cos(i) * sin(u)), 1e-9);
    NEAR(r0.z, 7000.0 * sin(i) * sin(u), 1e-9);
}

static void test_elliptical(void) {
    double a = 26000.0, e = 0.72;
    Object o = make(a, e, 63.4, 120.0, 270.0, 0.0);
    Prop p; prop_init(&p, &o, 0, o.epoch_jd);
    Vec3 r, v;
    prop_state(&p, 0, &r, &v);
    NEAR(nrm(r), a * (1 - e), 1e-8);                                   /* M=0 is perigee */
    double h0 = nrm(cross(r, v)), E0 = 0.5 * nrm(v) * nrm(v) - EARTH_MU_KM3S2 / nrm(r);
    NEAR(E0, -EARTH_MU_KM3S2 / (2 * a), 1e-9);                         /* specific energy = -mu/2a */
    NEAR(h0, sqrt(EARTH_MU_KM3S2 * a * (1 - e * e)), 1e-6);
    double period = 86400.0 / o.mean_motion;
    int ok = 1;
    for (int k = 1; k <= 300; k++) {
        double t = period * k / 97.0;
        prop_state(&p, t, &r, &v);
        double rr = nrm(r);
        if (rr < a * (1 - e) - 1e-6 || rr > a * (1 + e) + 1e-6) ok = 0;
        if (fabs(0.5 * nrm(v) * nrm(v) - EARTH_MU_KM3S2 / rr - E0) > 1e-9) ok = 0;   /* vis-viva */
        if (fabs(nrm(cross(r, v)) - h0) > 1e-5) ok = 0;                               /* angular momentum */
        /* finite-difference velocity */
        Vec3 rp, rm, vd; double h = 0.01;
        prop_state(&p, t + h, &rp, NULL); prop_state(&p, t - h, &rm, NULL);
        vd.x = (rp.x - rm.x) / (2 * h); vd.y = (rp.y - rm.y) / (2 * h); vd.z = (rp.z - rm.z) / (2 * h);
        if (ddist(vd, v) > 1e-7) ok = 0;
    }
    CHECK(ok);
    Vec3 ra; prop_state(&p, period / 2, &ra, NULL);
    NEAR(nrm(ra), a * (1 + e), 1e-6);                                  /* half a period later: apogee */
}

static void test_j2_and_epoch(void) {
    /* sun-synchronous: ~98.2 deg at 700 km altitude drifts RAAN by +0.9856 deg/day (360 deg per tropical year) */
    Object ss = make(EARTH_RADIUS_KM + 700.0, 0.0, 98.19, 0, 0, 0);
    double rd, ad, md;
    prop_j2_rates(&ss, &rd, &ad, &md);
    NEAR(rd * 86400.0 * 360.0 / TWO_PI_, 0.9856, 0.01);
    /* prograde 51.6 deg orbit regresses westward */
    Object iss = make(EARTH_RADIUS_KM + 420.0, 0.0005, 51.6, 0, 0, 0);
    prop_j2_rates(&iss, &rd, &ad, &md);
    CHECK(rd < 0);
    NEAR(rd * 86400.0 * 360.0 / TWO_PI_, -5.0, 0.4);
    /* critical inclination 63.435: apsidal drift vanishes */
    Object cr = make(26000.0, 0.7, 63.4349, 0, 0, 0);
    prop_j2_rates(&cr, &rd, &ad, &md);
    NEAR(ad, 0.0, 1e-12);
    /* J2 off: no secular rates */
    Prop p; prop_init(&p, &iss, 0, iss.epoch_jd);
    CHECK(p.raan_rate == 0 && p.argp_rate == 0 && p.m_rate == p.n);
    /* J2 on moves the plane: after one day the orbit normal has turned by about the RAAN rate */
    Prop q; prop_init(&q, &iss, 1, iss.epoch_jd);
    Vec3 r1, v1, r2, v2;
    prop_state(&q, 0, &r1, &v1); prop_state(&q, 86400.0, &r2, &v2);
    Vec3 h1 = cross(r1, v1), h2 = cross(r2, v2);
    double ang = acos((h1.x * h2.x + h1.y * h2.y + h1.z * h2.z) / (nrm(h1) * nrm(h2))) * 360.0 / TWO_PI_;
    NEAR(ang, 5.0 * sin(51.6 * TWO_PI_ / 360.0), 0.4);   /* node turns ~5 deg; the orbit normal turns by that times sin(i) */
    /* changing the reference epoch gives the same trajectory (with and without J2) */
    for (int j2 = 0; j2 <= 1; j2++) {
        Object o = make(7200.0, 0.01, 70.0, 33.0, 77.0, 120.0);
        Prop a, b; Vec3 ra, rb, va, vb;
        prop_init(&a, &o, j2, o.epoch_jd);
        prop_init(&b, &o, j2, o.epoch_jd + 3600.0 / 86400.0);
        prop_state(&a, 3600.0 + 1234.5, &ra, &va); prop_state(&b, 1234.5, &rb, &vb);
        /* limit: a Julian date in a double resolves ~40 us, i.e. ~0.3 m of along-track position */
        NEAR(ddist(ra, rb), 0.0, 5e-4); NEAR(ddist(va, vb), 0.0, 1e-6);
    }
}

static void test_drag(void) {
    Object o = make(6900.0, 0.001, 53.0, 10.0, 20.0, 30.0);
    o.mean_motion_dot = 2.0e-4;                                   /* rev/day^2 (TLE ndot/2 field) */
    Prop nd, dr, dj; prop_init(&nd, &o, 0, o.epoch_jd); prop_init(&dr, &o, PROP_DRAG, o.epoch_jd); prop_init(&dj, &o, PROP_J2 | PROP_DRAG, o.epoch_jd);
    CHECK(nd.m_quad == 0 && dr.m_quad > 0);
    NEAR(dr.m_quad, 2.0e-4 * TWO_PI_ / (86400.0 * 86400.0), 1e-20);
    /* mean anomaly gains ndot/2 * T^2: after one day, 2e-4 rev = 0.072 deg along the orbit (~0.14 km at r=6900 km * 1.26e-3 rad) */
    Vec3 a, b; prop_state(&nd, 86400.0, &a, NULL); prop_state(&dr, 86400.0, &b, NULL);
    NEAR(ddist(a, b), 6900.0 * 2.0e-4 * TWO_PI_, 0.01);
    prop_state(&nd, 0.0, &a, NULL); prop_state(&dr, 0.0, &b, NULL);
    NEAR(ddist(a, b), 0.0, 1e-9);                                  /* no effect at epoch */
    /* quadratic growth: 3 days gives 9x the one-day shift */
    Vec3 c, d; prop_state(&nd, 3 * 86400.0, &c, NULL); prop_state(&dr, 3 * 86400.0, &d, NULL);
    NEAR(ddist(c, d), 9.0 * 6900.0 * 2.0e-4 * TWO_PI_, 0.1);
    /* shifting the reference epoch must not change the trajectory, with drag */
    Prop e1, e2; Vec3 r1, r2, v1, v2;
    prop_init(&e1, &o, PROP_J2 | PROP_DRAG, o.epoch_jd); prop_init(&e2, &o, PROP_J2 | PROP_DRAG, o.epoch_jd + 2.0);
    prop_state(&e1, 2 * 86400.0 + 777.0, &r1, &v1); prop_state(&e2, 777.0, &r2, &v2);
    NEAR(ddist(r1, r2), 0.0, 5e-4); NEAR(ddist(v1, v2), 0.0, 1e-6);
    (void)dj;
}

/* v must be the exact derivative of r for every model variant (J2 drift, drag), on real orbits */
static void test_velocity_consistency(void) {
    Object objs[4] = {make(6900.0, 0.001, 53.0, 10.0, 20.0, 30.0), make(7100.0, 0.02, 98.2, 200.0, 100.0, 300.0),
                      make(26000.0, 0.72, 63.4, 120.0, 270.0, 10.0), make(42164.0, 0.0003, 0.1, 80.0, 10.0, 200.0)};
    objs[0].mean_motion_dot = 3e-4; objs[1].mean_motion_dot = 1e-3;
    int ok = 1;
    for (int k = 0; k < 4; k++) for (int fl = 0; fl <= 3; fl++) for (int j = 0; j < 25; j++) {
        Prop p; Vec3 r, v, rp, rm;
        prop_init(&p, &objs[k], fl, objs[k].epoch_jd + 0.37);
        double t = j * 3777.0 + 12.3, h = 0.05;
        prop_state(&p, t, &r, &v); prop_state(&p, t + h, &rp, NULL); prop_state(&p, t - h, &rm, NULL);
        Vec3 fd = {(rp.x - rm.x) / (2 * h), (rp.y - rm.y) / (2 * h), (rp.z - rm.z) / (2 * h)};
        if (ddist(fd, v) > 1e-7) { ok = 0; printf("  obj %d flags %d t=%g: |fd-v|=%g\n", k, fl, t, ddist(fd, v)); }
    }
    CHECK(ok);
}

static void test_planted(void) {
    Catalog c; DArray t; cat_init(&c); da_init(&t, sizeof(PlantedPair));
    CHECK(synth_generate(&c, &t, 2000, 5, 27) == 0);
    const Object *o = cat_objects(&c); int n = cat_count(&c);
    Prop *props = malloc(sizeof(Prop) * (size_t)n);
    double ref; iso_to_jd(SYNTH_EPOCH, &ref);
    for (int i = 0; i < n; i++) prop_init(&props[i], &o[i], 0, ref);
    int ok = 1;
    for (size_t k = 0; k < t.n; k++) {
        const PlantedPair *pp = (const PlantedPair *)da_at(&t, k);
        const Prop *A = &props[cat_find(&c, pp->id_a)], *B = &props[cat_find(&c, pp->id_b)];
        Vec3 ra, rb;
        prop_state(A, pp->tca_s, &ra, NULL); prop_state(B, pp->tca_s, &rb, NULL);
        if (fabs(ddist(ra, rb) - pp->miss_km) > 1e-4) { ok = 0; printf("  pair %zu: %.6f vs %.6f\n", k, ddist(ra, rb), pp->miss_km); }
        double best = 1e30;
        for (double s = pp->tca_s - 30; s <= pp->tca_s + 30; s += 0.02) {
            prop_state(A, s, &ra, NULL); prop_state(B, s, &rb, NULL);
            double d = ddist(ra, rb); if (d < best) best = d;
        }
        if (best < pp->miss_km - 1e-4) ok = 0;
    }
    CHECK(ok);
    free(props); cat_free(&c); da_free(&t);
}

static void test_real(void) {
    Catalog c; cat_init(&c); LoadReport r; report_init(&r);
    int rc = cat_load_csv(&c, "data/raw/active.csv", &r);
    report_free(&r);
    if (rc != 0) { printf("  (skipping real-data propagation check)\n"); cat_free(&c); return; }
    int n = cat_count(&c), ok = 1;
    const Object *o = cat_objects(&c);
    for (int j2 = 0; j2 <= 3; j2 += (j2 == 1 ? 2 : 1)) for (int i = 0; i < n; i++) {
        Prop p; Vec3 rr, vv; prop_init(&p, &o[i], j2, o[i].epoch_jd);
        double times[3] = {0.0, 86400.0, 7 * 86400.0};
        for (int k = 0; k < 3; k++) {
            prop_state(&p, times[k], &rr, &vv);
            double rad = nrm(rr);
            if (!isfinite(rad) || rad < p.a * (1 - p.e) - 1e-6 || rad > p.a * (1 + p.e) + 1e-6) ok = 0;   /* bounds use the propagator's own a */
            if (fabs(p.a - o[i].sma_km) > 8.0) ok = 0;       /* J2 'a' differs from the no-J2 value by up to ~6 km (equatorial orbits) */
        }
    }
    CHECK(ok);
    cat_free(&c);
}

int main(void) {
    test_kepler(); test_circular(); test_elliptical(); test_j2_and_epoch(); test_drag(); test_velocity_consistency(); test_planted(); test_real();
    printf("%ld checks, %ld failed\n", checks, fails);
    return fails ? 1 : 0;
}
