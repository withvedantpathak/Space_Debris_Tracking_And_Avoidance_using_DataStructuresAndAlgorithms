#ifndef SCREEN_H
#define SCREEN_H
#include "conjunction.h"

/* C5: fast screening. Every method uses the same per-pair test and refinement as the brute-force baseline
 * (conj_bruteforce), so for the proven-bound methods the result is IDENTICAL to the baseline; they only differ in
 * how candidate pairs are found at each time step instead of testing all n(n-1)/2 pairs.
 *
 *  SCREEN_KDTREE  rebuild a k-d tree each step; for each object a radius query with reach thr + (vmax_i + vmax_max) h.
 *  SCREEN_GRID    same with the uniform hashed grid (cell = largest reach).
 *  SCREEN_ALT     altitude filter: objects sorted by perigee radius; a sweep skips pairs whose radial ranges
 *                 [rmin - pad, rmax + pad] are further apart than the threshold. NOT a proven bound under SGP4
 *                 (the radius leaves the mean-element range by a few km); pad is empirical and the largest
 *                 violation seen is reported in ScreenStats.alt_violation_km.
 *
 * Proof for KDTREE/GRID: a pair with an event in [t0,t1] has d0 <= thr + (vmax_i + vmax_j) h (Lipschitz bound, as in the
 * baseline), and vmax_j <= vmax_max, so it is returned by the query of the lower-numbered object.
 * Objects that cannot be propagated (prop_ok == 0) or are parked at a step are left out of the index. */
typedef enum { SCREEN_KDTREE = 1, SCREEN_GRID = 2, SCREEN_ALT = 3 } ScreenMethod;

typedef struct {
    ScreenMethod method;
    const double *rmin_km, *rmax_km;   /* SCREEN_ALT: per-object perigee/apogee radius (km from Earth centre) */
    double alt_pad_km;                 /* SCREEN_ALT: padding on each side */
    double alt_violation_km;           /* out: largest distance by which an object left [rmin, rmax] (unpadded) over the window */
    int refine_linear;                 /* 1 = closest approach by linear interpolation of the positions between steps (Gomez et al.); 2 = hybrid: linear, but bisection for pairs with relative speed < 1 km/s; 0 = bisection (default) */
    double cell_scale;                 /* SCREEN_GRID: cell size = cell_scale * largest reach (0 means 1; otherwise must be in [0.05, 100], else conj_fast returns -4) */
    long index_pairs;                  /* out: candidate pairs returned by the index (the pair_checks of this method) */
} FastScreen;

/* Same contract as conj_bruteforce: events appended to out (sorted by a, b, tca); stats->pair_checks counts the
 * pairs the index handed over. Returns 0, -3 on out of memory, -4 bad option. */
int conj_fast(const Prop *props, int n, const ScreenParams *sp, FastScreen *fs, DArray *out, ScreenStats *st);
#endif
