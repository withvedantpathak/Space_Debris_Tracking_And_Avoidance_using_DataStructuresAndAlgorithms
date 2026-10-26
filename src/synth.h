#ifndef SYNTH_H
#define SYNTH_H
#include "catalog.h"

/* Seeded synthetic catalogue with planted close pairs of known closest approach.
 *
 * Filler objects: altitude, inclination and orbit-class mix are sampled from histograms of
 * the CelesTrak snapshot in data/raw (see docs/synthetic_model.md); all other elements are
 * simple parametric draws. All objects share one epoch (SYNTH_EPOCH). NORAD ids start at 900000.
 *
 * Planted pairs: two circular orbits (e=0) of equal radius whose planes intersect; B is
 * offset along its track until the TWO-BODY closest approach equals the designed miss
 * distance. The truth (tca_s, miss_km) is obtained numerically from the exact circular-orbit
 * motion, so a correct propagator + refinement must reproduce it. BOTH orbits repeat every period_s and
 * their planes cross at two antipodal points, so the same encounter (same miss distance, by symmetry)
 * recurs at tca_s + k * period_s / 2 for every integer k (also negative k) that falls inside the window.
 * Valid for two-body propagation (flags = 0); with J2 the planes drift and the truth no longer applies. Other (random) objects may
 * of course also have conjunctions; the truth lists only the planted ones. */
#define SYNTH_EPOCH "2026-10-20T00:00:00"
#define SYNTH_FIRST_ID 900000

typedef struct {
    int    id_a, id_b;
    double tca_s;           /* time of closest approach, seconds after epoch */
    double miss_km;         /* minimum distance between A and B (two-body, exact to ~1e-7 km) */
    double rel_speed_kms;   /* relative speed at closest approach */
    double period_s;        /* orbital period of both objects */
} PlantedPair;

/* n objects in total (2*pairs of them planted). truth receives PlantedPair records
 * (may be NULL). Returns 0, -1 bad arguments, -3 out of memory. Deterministic per seed. */
int synth_generate(Catalog *cat, DArray *truth, int n, unsigned seed, int pairs);
int synth_write_truth(const DArray *truth, const char *path);

/* Test hooks (used only by tests/test_catalog.c). */
int synth_test_pick(const double *w, int n, unsigned seed);             /* weighted index pick; n-1 if the weights are exhausted */
int synth_test_plant(unsigned seed, int k, int max_tries);              /* place one pair; 0 ok, -1 gave up after max_tries */
#endif
