#ifndef SGP4_H
#define SGP4_H
#include "catalog.h"
#include "ds/vec3.h"

/* Near-earth SGP4, hand-written from Spacetrack Report #3 as revised in Vallado, Crawford, Hujsak,
 * Kelso, "Revisiting Spacetrack Report #3" (AIAA 2006-6753). WGS-72 constants, as used for TLE/OMM
 * element sets. Output: TEME position (km) and velocity (km/s).
 *
 * Deep-space objects (period >= 225 min, SDP4) are NOT implemented: sgp4_init() returns
 * SGP4_ERR_DEEP_SPACE for them. Validated against the Python sgp4 package: tests/test_sgp4.c
 * (fixed cases) and scripts/validate_sgp4.py (whole real catalogue: positions agree to < 1 mm).
 *
 * Note: like the reference code, the velocity is SGP4's own analytic output, which omits the time derivative of the
 * short-period terms; it differs from d(position)/dt by up to ~1.5 m/s. Harmless for km-scale screening (speed bounds carry a 1 %
 * margin, ~70 m/s), but it is not an exact derivative as the Kepler-model velocity is. */
#define SGP4_OK              0
#define SGP4_ERR_MEAN_ECC    1   /* mean eccentricity out of range during propagation */
#define SGP4_ERR_MEAN_MOTION 2   /* mean motion <= 0 */
#define SGP4_ERR_PERT_ECC    3   /* perturbed eccentricity out of range */
#define SGP4_ERR_SEMILATUS   4   /* semi-latus rectum < 0 */
#define SGP4_ERR_DECAYED     6   /* satellite has decayed (r < 1 Earth radius) */
#define SGP4_ERR_DEEP_SPACE  7   /* period >= 225 min: not supported here */
#define SGP4_ERR_INPUT       8   /* eccentricity or mean motion invalid at init */

typedef struct {
    int isimp;
    double ecco, inclo, nodeo, argpo, mo, no, bstar;   /* radians, rad/min; no = un-Kozai'd mean motion */
    double mdot, argpdot, nodedot, nodecf, omgcof, xmcof, xlcof, aycof;
    double cc1, cc4, cc5, d2, d3, d4, t2cof, t3cof, t4cof, t5cof;
    double eta, delmo, sinmao, x1mth2, x7thm1, con41;
} Sgp4;

/* Initialise from an object's mean elements (OMM: rev/day, degrees). Returns SGP4_OK or an error. */
int sgp4_init(Sgp4 *s, const Object *o);

/* State tsince_min minutes after the element epoch. v may be NULL. Returns SGP4_OK or an error;
 * on error r and v are left untouched. */
int sgp4_state(const Sgp4 *s, double tsince_min, Vec3 *r, Vec3 *v);
#endif
