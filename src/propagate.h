#ifndef PROPAGATE_H
#define PROPAGATE_H
#include "catalog.h"
#include "ds/vec3.h"
#include "sgp4.h"

/* Orbit propagation from mean elements: Kepler two-body motion, optionally with J2 secular
 * drift of RAAN, argument of perigee and mean anomaly. Output is an inertial Earth-centred
 * frame (the TLE/OMM frame; Earth rotation is not applied, screening needs relative geometry
 * only). Units: km, km/s, seconds.
 *
 * With j2=1 the OMM mean motion is treated as the Kozai mean motion (as SGP4 does): it is the
 * mean-anomaly rate, and the semi-major axis comes from the un-Kozai'd Keplerian mean motion.
 * Object.sma_km (from mean motion, no J2) can therefore differ from the propagated a by up to ~6 km (largest for equatorial orbits; measured on the real catalogue).
 *
 * Not modelled: drag, solar pressure, third bodies, short-period J2 terms, SGP4's Brouwer/Kozai
 * mean-element conversion. Errors against SGP4 are measured in results/propagation_error.csv. */
#define EARTH_J2 1.08262668e-3
#define PROP_J2   1       /* secular J2 drift of RAAN, argp, mean anomaly */
#define PROP_DRAG 2       /* quadratic mean-anomaly term from the element set's MEAN_MOTION_DOT */
#define PROP_SGP4 8       /* use src/sgp4.c (near-earth SGP4, TEME) for position and velocity; the Kepler fields stay
                             filled (J2+drag) because the screening speed bound uses a and e. Objects SGP4 cannot
                             propagate (deep space, decayed) are parked far away and never pair up: see prop_ok(). */

typedef struct {
    double a, e, sqrt1me2, cosi, sini;
    double n;                         /* mean motion (rad/s), two-body */
    double M0, raan0, argp0;          /* radians at the reference time */
    double m_rate, raan_rate, argp_rate;   /* rad/s; m_rate = n when J2 is off */
    int sgp4_on, sgp4_ok;                  /* PROP_SGP4 requested / init succeeded */
    Sgp4 sg; double sg_t0_min;             /* SGP4 state and minutes from element epoch to ref_jd */
    double park;                           /* x coordinate (km) used when SGP4 fails */
    double m_quad;                         /* rad/s^2: drag term, M += m_quad * T^2 (T = seconds since element epoch) */
} Prop;

/* Solve Kepler's equation M = E - e sin E for E (radians). Newton iteration, quadratic
 * convergence, <= 50 iterations, tolerance 1e-14. Works for 0 <= e < 1. */
double kepler_E(double M, double e);

/* Secular J2 rates in rad/s for the object (exposed for tests). m_dot is the Kozai mean motion. */
void prop_j2_rates(const Object *o, double *raan_dot, double *argp_dot, double *m_dot);

/* Prepare an object for propagation. Time t=0 is the Julian date ref_jd (may differ from the
 * object's own epoch: the element set is first advanced to ref_jd). flags = 0 (two-body),
 * PROP_J2, or PROP_J2 | PROP_DRAG (PROP_DRAG alone is allowed). */
void prop_init(Prop *p, const Object *o, int flags, double ref_jd);

/* 1 if the object can be propagated by the requested model (always 1 unless PROP_SGP4 init failed). */
int prop_ok(const Prop *p);

/* Position (km) and velocity (km/s) t seconds after ref_jd. O(1): one Kepler solve. v may be NULL. */
void prop_state(const Prop *p, double t, Vec3 *r, Vec3 *v);

/* Positions of all objects at time t (seconds after ref_jd); out[i] for object i. O(n). */
void prop_positions(const Prop *props, int n, double t, Vec3 *out);
#endif
