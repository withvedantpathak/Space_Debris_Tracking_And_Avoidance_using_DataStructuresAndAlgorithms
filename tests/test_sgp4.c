/* Tests for sgp4.c: reference vectors from the Python sgp4 package (tests/sgp4_vectors.h), error paths, and
 * consistency between velocity and the numerical derivative of position. */
#include <math.h>
#include <stdio.h>
#include <string.h>
#include "sgp4.h"
#include "propagate.h"
#include "sgp4_vectors.h"

static long checks = 0, fails = 0;
#define CHECK(c) do { checks++; if (!(c)) { fails++; printf("FAIL %s:%d: %s\n", __FILE__, __LINE__, #c); } } while (0)
#define NEAR(a, b, tol) do { checks++; if (!(fabs((a) - (b)) <= (tol))) { fails++; printf("FAIL %s:%d: %s=%.12g vs %.12g (tol %g)\n", __FILE__, __LINE__, #a, (double)(a), (double)(b), (double)(tol)); } } while (0)

static Object make(const Sgp4Case *c) {
    Object o; memset(&o, 0, sizeof o);
    o.norad_id = 1; o.mean_motion = c->mm; o.ecc = c->ecc; o.incl = c->inc; o.raan = c->raan; o.argp = c->argp; o.mean_anom = c->ma; o.bstar = c->bstar;
    return o;
}

int main(void) {
    int ncase = (int)(sizeof SGP4_CASES / sizeof SGP4_CASES[0]);
    for (int k = 0; k < ncase; k++) {
        const Sgp4Case *c = &SGP4_CASES[k];
        Object o = make(c); Sgp4 s;
        CHECK(sgp4_init(&s, &o) == SGP4_OK);
        for (int j = 0; j < c->nt; j++) {
            Vec3 r, v;
            int err = sgp4_state(&s, c->t[j], &r, &v);
            CHECK(err == SGP4_OK);
            NEAR(r.x, c->r[j][0], 1e-7); NEAR(r.y, c->r[j][1], 1e-7); NEAR(r.z, c->r[j][2], 1e-7);
            NEAR(v.x, c->v[j][0], 1e-10); NEAR(v.y, c->v[j][1], 1e-10); NEAR(v.z, c->v[j][2], 1e-10);
            /* velocity vs central difference of position (SGP4's velocity is analytic, not exactly the derivative of r,
               so the bound is loose: SGP4 omits derivatives of the short-period terms; measured gap up to ~1.5 m/s) -- catches sign/axis mistakes, not rounding */
            if (k != 2 || j < 2) {
                Vec3 ra, rb; double h = 0.01;   /* minutes */
                if (c->t[j] > h && sgp4_state(&s, c->t[j] - h, &ra, NULL) == SGP4_OK && sgp4_state(&s, c->t[j] + h, &rb, NULL) == SGP4_OK) {
                    NEAR(v.x, (rb.x - ra.x) / (2 * h * 60.0), 3e-3); NEAR(v.y, (rb.y - ra.y) / (2 * h * 60.0), 3e-3); NEAR(v.z, (rb.z - ra.z) / (2 * h * 60.0), 3e-3);
                }
            }
        }
    }
    /* v may be NULL */
    { Object o = make(&SGP4_CASES[0]); Sgp4 s; Vec3 r; CHECK(sgp4_init(&s, &o) == SGP4_OK); CHECK(sgp4_state(&s, 10.0, &r, NULL) == SGP4_OK); }
    /* deep space refused: GEO (1.0027 rev/day) and Molniya (2 rev/day) */
    { Object o = make(&SGP4_CASES[0]); Sgp4 s; o.mean_motion = 1.0027; CHECK(sgp4_init(&s, &o) == SGP4_ERR_DEEP_SPACE); o.mean_motion = 2.0; o.ecc = 0.7; CHECK(sgp4_init(&s, &o) == SGP4_ERR_DEEP_SPACE); }
    /* boundary: period just under / over 225 min is 6.4 rev/day */
    { Object o = make(&SGP4_CASES[0]); Sgp4 s; o.mean_motion = 6.41; CHECK(sgp4_init(&s, &o) == SGP4_OK); o.mean_motion = 6.39; CHECK(sgp4_init(&s, &o) == SGP4_ERR_DEEP_SPACE); }
    /* invalid input */
    { Object o = make(&SGP4_CASES[0]); Sgp4 s; o.ecc = 1.0; CHECK(sgp4_init(&s, &o) == SGP4_ERR_INPUT); o.ecc = -0.1; CHECK(sgp4_init(&s, &o) == SGP4_ERR_INPUT);
      o = make(&SGP4_CASES[0]); o.mean_motion = 0.0; CHECK(sgp4_init(&s, &o) == SGP4_ERR_INPUT); }
    /* decayed: strong drag, long time -> must return an error code, never garbage */
    { Object o = make(&SGP4_CASES[2]); Sgp4 s; Vec3 r; o.bstar = 0.05; CHECK(sgp4_init(&s, &o) == SGP4_OK);
      int bad = 0; for (double t = 0; t < 20000; t += 50) if (sgp4_state(&s, t, &r, NULL) != SGP4_OK) bad = 1;
      CHECK(bad == 1); }
    /* integration with Prop: PROP_SGP4 matches sgp4_state; failed objects are parked and distinct */
    { Object o = make(&SGP4_CASES[0]); o.epoch_jd = 2461000.5; Prop p; Sgp4 s; Vec3 r1, r2;
      prop_init(&p, &o, PROP_SGP4, o.epoch_jd + 0.25);
      CHECK(prop_ok(&p)); sgp4_init(&s, &o); prop_state(&p, 600.0, &r1, NULL); sgp4_state(&s, 0.25 * 1440.0 + 10.0, &r2, NULL);
      NEAR(r1.x, r2.x, 1e-6); NEAR(r1.y, r2.y, 1e-6); NEAR(r1.z, r2.z, 1e-6);
      Object g1 = o, g2 = o; g1.mean_motion = g2.mean_motion = 1.0; g1.norad_id = 5; g2.norad_id = 6; Prop a, b; Vec3 ra, rb;
      prop_init(&a, &g1, PROP_SGP4, o.epoch_jd); prop_init(&b, &g2, PROP_SGP4, o.epoch_jd);
      CHECK(!prop_ok(&a)); prop_state(&a, 0, &ra, NULL); prop_state(&b, 0, &rb, NULL); CHECK(fabs(ra.x - rb.x) >= 999.0 && ra.x > 1e8); }
    printf("%ld checks, %ld failed\n", checks, fails);
    return fails != 0;
}
