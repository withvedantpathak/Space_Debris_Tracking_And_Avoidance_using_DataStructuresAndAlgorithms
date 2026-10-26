#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "synth.h"
#include "ds/vec3.h"

/* ---------------- RNG: xorshift64*, seeded through splitmix64 ---------------- */
typedef struct { unsigned long long s; } Rng;
static void rng_seed(Rng *r, unsigned seed) {
    unsigned long long z = (unsigned long long)seed + 0x9e3779b97f4a7c15ULL;
    z = (z ^ (z >> 30)) * 0xbf58476d1ce4e5b9ULL;
    z = (z ^ (z >> 27)) * 0x94d049bb133111ebULL;
    r->s = (z ^ (z >> 31)) | 1ULL;
}
static unsigned long long rng_next(Rng *r) {
    r->s ^= r->s >> 12; r->s ^= r->s << 25; r->s ^= r->s >> 27;
    return r->s * 2685821657736338717ULL;
}
static double rng_u(Rng *r) { return (double)(rng_next(r) >> 11) / 9007199254740992.0; }
static double rng_range(Rng *r, double lo, double hi) { return lo + (hi - lo) * rng_u(r); }
static int rng_int(Rng *r, int n) { return (int)(rng_next(r) % (unsigned long long)n); }
static double rng_norm(Rng *r) {                       /* Box-Muller */
    double u1 = rng_u(r), u2 = rng_u(r);
    if (u1 < 1e-300) u1 = 1e-300;
    return sqrt(-2.0 * log(u1)) * cos(TWO_PI_ * u2);
}
static int pick(Rng *r, const double *w, int n) {
    double tot = 0, x;
    for (int i = 0; i < n; i++) tot += w[i];
    x = rng_u(r) * tot;
    for (int i = 0; i < n; i++) { if (x < w[i]) return i; x -= w[i]; }
    return n - 1;
}

/* Counts per bin from the real snapshot (19,440 objects, LEO = a-R < 2000 km). */
static const double LEO_ALT_W[20] = {0, 3, 79, 1066, 10618, 1776, 1149, 965, 1067, 471, 377, 595, 346, 8, 91, 6, 2, 5, 1, 1}; /* 100 km bins from 0 km */
static const double INC_W[18] = {10, 6, 20, 158, 3924, 5757, 580, 1010, 1417, 5654, 85, 0, 2, 1, 2, 0, 0, 0};              /* 10 deg bins from 0 deg */
#define P_NON_LEO (814.0 / 19440.0)        /* share of objects with a-R >= 2000 km */
#define P_GEO_OF_NON_LEO 0.7               /* coarse split of the non-LEO share */

/* ---------------- vector helpers ---------------- */
static Vec3 v3(double x, double y, double z) { Vec3 v = {x, y, z}; return v; }
static double dot(Vec3 a, Vec3 b) { return a.x * b.x + a.y * b.y + a.z * b.z; }
static Vec3 cross(Vec3 a, Vec3 b) { return v3(a.y * b.z - a.z * b.y, a.z * b.x - a.x * b.z, a.x * b.y - a.y * b.x); }
static double vnorm(Vec3 a) { return sqrt(dot(a, a)); }
static Vec3 vscale(Vec3 a, double s) { return v3(a.x * s, a.y * s, a.z * s); }

/* ---------------- exact circular-orbit motion (e = 0, argp = 0) ---------------- */
typedef struct { double a, inc, raan, m0; } Circ;         /* radians; m0 = argument of latitude at epoch */
static Vec3 circ_pos(const Circ *c, double n, double t) {
    double u = c->m0 + n * t;
    double cO = cos(c->raan), sO = sin(c->raan), ci = cos(c->inc), si = sin(c->inc), cu = cos(u), su = sin(u);
    return v3(c->a * (cO * cu - sO * ci * su), c->a * (sO * cu + cO * ci * su), c->a * si * su);
}
static double circ_dist(const Circ *A, const Circ *B, double n, double t) {
    Vec3 p = circ_pos(A, n, t), q = circ_pos(B, n, t);
    return vnorm(v3(p.x - q.x, p.y - q.y, p.z - q.z));
}
/* Golden-section minimum of the pair distance on [lo,hi] (unimodal on these short windows). */
static double circ_min(const Circ *A, const Circ *B, double n, double lo, double hi, double *tmin) {
    const double g = 0.6180339887498949;
    double x1 = hi - g * (hi - lo), x2 = lo + g * (hi - lo);
    double f1 = circ_dist(A, B, n, x1), f2 = circ_dist(A, B, n, x2);
    for (int i = 0; i < 80; i++) {
        if (f1 < f2) { hi = x2; x2 = x1; f2 = f1; x1 = hi - g * (hi - lo); f1 = circ_dist(A, B, n, x1); }
        else         { lo = x1; x1 = x2; f1 = f2; x2 = lo + g * (hi - lo); f2 = circ_dist(A, B, n, x2); }
    }
    *tmin = 0.5 * (lo + hi);
    return circ_dist(A, B, n, *tmin);
}
static void plane_axes(double raan, double inc, Vec3 *ahat, Vec3 *bhat) {
    *ahat = v3(cos(raan), sin(raan), 0);
    *bhat = v3(-sin(raan) * cos(inc), cos(raan) * cos(inc), sin(inc));
}
static double wrap_deg(double rad) {
    double d = fmod(rad * 180.0 / (TWO_PI_ / 2), 360.0);
    return d < 0 ? d + 360.0 : d;
}

static void fill_object(Object *o, int id, const char *name, double epoch, double a, double ecc,
                        double inc_deg, double raan_deg, double argp_deg, double ma_deg) {
    memset(o, 0, sizeof *o);
    o->norad_id = id;
    snprintf(o->name, sizeof o->name, "%s", name);
    snprintf(o->intl_id, sizeof o->intl_id, "SYN-%d", id - SYNTH_FIRST_ID);
    o->epoch_jd = epoch;
    o->mean_motion = sqrt(EARTH_MU_KM3S2 / (a * a * a)) * 86400.0 / TWO_PI_;
    o->ecc = ecc; o->incl = inc_deg; o->raan = raan_deg; o->argp = argp_deg; o->mean_anom = ma_deg;
    obj_derive(o);
}

static void make_filler(Rng *r, Object *o, int id, double epoch) {
    double alt;
    if (rng_u(r) < P_NON_LEO) {
        alt = rng_u(r) < P_GEO_OF_NON_LEO ? 35786.0 + rng_range(r, -150, 150) : rng_range(r, 2000, 25000);
    } else {
        alt = (pick(r, LEO_ALT_W, 20) + rng_u(r)) * 100.0;
        if (alt < 200) alt = 200;
    }
    double a = EARTH_RADIUS_KM + alt;
    double inc = (pick(r, INC_W, 18) + rng_u(r)) * 10.0;
    double ecc = 2e-4 * exp(2.56 * rng_norm(r));               /* lognormal: median 2e-4, p90 ~ 5e-3 */
    double emax = 1.0 - (EARTH_RADIUS_KM + 200.0) / a;         /* keep perigee >= 200 km altitude */
    if (ecc > 0.05) ecc = 0.05;
    if (ecc > emax) ecc = emax;
    char nm[32];
    snprintf(nm, sizeof nm, "SYN-%06d", id - SYNTH_FIRST_ID);
    fill_object(o, id, nm, epoch, a, ecc, inc, rng_range(r, 0, 360), rng_range(r, 0, 360), rng_range(r, 0, 360));
}

static const double MISS_TARGETS[9] = {0.5, 1, 2, 4, 6, 8, 12, 15, 25};   /* km: first six below 10 km, last three above */

#define SYNTH_MAX_TRIES 1000   /* attempts to place one planted pair before giving up */

static int make_pair(Rng *r, Object *oa, Object *ob, int ida, int idb, int k, double epoch, PlantedPair *t, int max_tries) {
    for (int tries = 0; tries < max_tries; tries++) {
        double alt = rng_range(r, 450, 1000), a = EARTH_RADIUS_KM + alt;
        double ia = rng_range(r, 40, 140) * (TWO_PI_ / 360.0);
        double di = rng_range(r, 15, 60) * (TWO_PI_ / 360.0) * (rng_u(r) < 0.5 ? -1 : 1);
        double ib = ia + di;
        if (ib < 0.17 || ib > 2.97) continue;
        Circ A = {a, ia, rng_range(r, 0, TWO_PI_), 0}, B = {a, ib, rng_range(r, 0, TWO_PI_), 0};
        Vec3 aA, bA, aB, bB;
        plane_axes(A.raan, A.inc, &aA, &bA); plane_axes(B.raan, B.inc, &aB, &bB);
        Vec3 L = cross(cross(aA, bA), cross(aB, bB));
        if (vnorm(L) < 0.3) continue;
        L = vscale(L, 1.0 / vnorm(L));                          /* direction of the planes' intersection line */
        double uA = atan2(dot(L, bA), dot(L, aA)), uB = atan2(dot(L, bB), dot(L, aB));
        Vec3 vA = v3(-sin(uA) * aA.x + cos(uA) * bA.x, -sin(uA) * aA.y + cos(uA) * bA.y, -sin(uA) * aA.z + cos(uA) * bA.z);
        Vec3 vB = v3(-sin(uB) * aB.x + cos(uB) * bB.x, -sin(uB) * aB.y + cos(uB) * bB.y, -sin(uB) * aB.z + cos(uB) * bB.z);
        double cth = dot(vA, vB), theta = acos(cth);
        if (theta < 20.0 * TWO_PI_ / 360.0 || theta > 150.0 * TWO_PI_ / 360.0) continue;

        double n = sqrt(EARTH_MU_KM3S2 / (a * a * a));
        double tstar = rng_range(r, 1800, 72000);               /* crossing time: 0.5 h .. 20 h after epoch */
        double target = MISS_TARGETS[k % 9];
        A.m0 = uA - n * tstar;
        /* Offset B along its track; start from the straight-line estimate, then correct by iteration. */
        double delta = target / (a * cos(0.5 * theta)), m, tc;
        B.m0 = uB + delta - n * tstar;
        m = circ_min(&A, &B, n, tstar - 60, tstar + 60, &tc);
        for (int it = 0; it < 12 && fabs(m - target) > 1e-7; it++) {
            delta *= target / m;
            B.m0 = uB + delta - n * tstar;
            m = circ_min(&A, &B, n, tstar - 60, tstar + 60, &tc);
        }
        if (fabs(m - target) > 1e-6) continue;

        char nm[32];
        snprintf(nm, sizeof nm, "SYN-PLANT-%03d-A", k);
        fill_object(oa, ida, nm, epoch, a, 0.0, wrap_deg(A.inc), wrap_deg(A.raan), 0.0, wrap_deg(A.m0));
        snprintf(nm, sizeof nm, "SYN-PLANT-%03d-B", k);
        fill_object(ob, idb, nm, epoch, a, 0.0, wrap_deg(B.inc), wrap_deg(B.raan), 0.0, wrap_deg(B.m0));
        t->id_a = ida; t->id_b = idb; t->tca_s = tc; t->miss_km = m;
        t->rel_speed_kms = sqrt(EARTH_MU_KM3S2 / a) * 2.0 * sin(0.5 * theta);
        t->period_s = TWO_PI_ / n;
        return 0;
    }
    return -1;
}

int synth_generate(Catalog *cat, DArray *truth, int n, unsigned seed, int pairs) {
    if (n < 1 || pairs < 0 || 2 * pairs > n) return -1;
    double epoch;
    iso_to_jd(SYNTH_EPOCH, &epoch);
    Rng r;
    rng_seed(&r, seed);
    Object *all = malloc(sizeof(Object) * (size_t)n);
    if (!all) return -3;
    int rc = 0;
    for (int k = 0; k < pairs && !rc; k++) {
        PlantedPair t;
        if (make_pair(&r, &all[2 * k], &all[2 * k + 1], SYNTH_FIRST_ID + 2 * k, SYNTH_FIRST_ID + 2 * k + 1, k, epoch, &t, SYNTH_MAX_TRIES)) rc = -1;
        else if (truth && da_push(truth, &t) != 0) rc = -3;
    }
    for (int i = 2 * pairs; i < n && !rc; i++) make_filler(&r, &all[i], SYNTH_FIRST_ID + i, epoch);
    if (!rc) {
        for (int i = n - 1; i > 0; i--) { int j = rng_int(&r, i + 1); Object t = all[i]; all[i] = all[j]; all[j] = t; }   /* Fisher-Yates */
        for (int i = 0; i < n && !rc; i++) if (cat_add(cat, &all[i]) < 0) rc = -3;
    }
    free(all);
    return rc;
}

int synth_write_truth(const DArray *truth, const char *path) {
    FILE *f = fopen(path, "w");
    if (!f) return -1;
    fprintf(f, "id_a,id_b,tca_s_after_epoch,miss_km,rel_speed_kms,period_s\n");
    for (size_t i = 0; i < truth->n; i++) {
        const PlantedPair *p = (const PlantedPair *)da_at(truth, i);
        fprintf(f, "%d,%d,%.6f,%.9f,%.6f,%.6f\n", p->id_a, p->id_b, p->tca_s, p->miss_km, p->rel_speed_kms, p->period_s);
    }
    return fclose(f) == 0 ? 0 : -1;
}

/* Test hooks: reach the defensive fallbacks (pick() rounding fallback, give-up in make_pair) that normal seeds never hit. */
int synth_test_pick(const double *w, int n, unsigned seed) {
    Rng r;
    rng_seed(&r, seed);
    return pick(&r, w, n);
}

int synth_test_plant(unsigned seed, int k, int max_tries) {
    Rng r; Object a, b; PlantedPair t; double epoch;
    rng_seed(&r, seed);
    iso_to_jd(SYNTH_EPOCH, &epoch);
    return make_pair(&r, &a, &b, SYNTH_FIRST_ID, SYNTH_FIRST_ID + 1, k, epoch, &t, max_tries);
}
