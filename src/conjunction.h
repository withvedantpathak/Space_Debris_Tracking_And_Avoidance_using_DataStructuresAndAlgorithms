#ifndef CONJUNCTION_H
#define CONJUNCTION_H
#include "propagate.h"
#include "ds/darray.h"

/* Conjunction events and the exact brute-force baseline (ground truth for C5).
 *
 * Definition. For a pair of objects let d(t) be their distance under the propagation model.
 * A CONJUNCTION is a local minimum of d(t) strictly inside the screening window (0, T] with
 * d_min < threshold. Minima at the window start (pairs already separating) and pairs that never
 * approach (no minimum) are not events. A pair can have several events (one per close pass).
 * Pairs that stay within the threshold for the whole window (co-orbiting, e.g. docked or
 * launch-train objects) have no interior minimum, so they do not appear; see docs/known_limitations.md.
 *
 * Time is in seconds after the Prop reference epoch (prop_init ref_jd). */
typedef struct { int a, b; double tca, miss_km, rel_speed_kms; } Conjunction;   /* a < b: indices into the object array */

/* infl_km (optional, may be NULL): per-object extra distance, so the threshold of the pair (i, j) is threshold_km + infl_km[i] + infl_km[j].
 * Used for the epoch-age-aware screening distance (uncertainty grows with the age of the element set). */
typedef struct { double window_s, step_s, threshold_km; const double *infl_km; } ScreenParams;
typedef struct { long steps, pair_checks, candidates, refinements, events, close_at_start; } ScreenStats;   /* close_at_start: pairs already closer than the threshold at t = 0 (co-orbiting or about to pass) */

/* Conservative upper bound on the speed of the model (km/s): perigee speed with 1% margin. */
double conj_vmax(const Prop *p);

/* Closest approach of one pair inside [t0,t1]. Finds the root of g(t) = r_rel . v_rel by bisection
 * (60 iterations), valid when g changes sign from negative to non-negative on the interval, i.e. the
 * distance has a minimum there. Returns 1 and fills *tmin, *dmin (km), *vrel (km/s) if so, else 0.
 * Assumes one minimum per interval, which holds for steps well below the relative orbital period. */
int conj_refine(const Prop *pa, const Prop *pb, double t0, double t1, double *tmin, double *dmin, double *vrel);

/* Exact all-pairs screening, O(n^2 * T/step) pair checks. A pair/interval is discarded only when a
 * proven bound shows d(t) cannot drop below the threshold: d is Lipschitz with constant L = vmax_a +
 * vmax_b, so min over an interval of width h with endpoint distances d0,d1 is >= (d0 + d1 - L h)/2.
 * Surviving intervals are refined with conj_refine. Result: Conjunction records appended to out
 * (sorted by a, b, tca). Returns 0, or -3 on out of memory. */
/* 0 if the parameters are usable, -4 otherwise (non-positive or non-finite step, negative or non-finite window/threshold/inflation, more than 2e9 steps). */
int conj_check_params(const ScreenParams *sp, int n);
int conj_bruteforce(const Prop *props, int n, const ScreenParams *sp, DArray *out, ScreenStats *st);

void conj_sort(DArray *ev);
void conj_sort_range(DArray *ev, size_t first);                 /* sort only ev[first..) */                                    /* by (a, b, tca) */

typedef struct { long truth, test, matched, missed, extra; } CompareResult;
/* Match events of `test` against `truth`: same (a,b) and |tca difference| <= tca_tol_s.
 * missed = truth events not found (recall loss), extra = test events with no truth partner.
 * Both arrays must be sorted with conj_sort. */
void conj_compare(const DArray *truth, const DArray *test, double tca_tol_s, CompareResult *res);

/* Keep only the closest event of each pair. Input sorted by conj_sort; output sorted too. */
int conj_best_per_pair(const DArray *in, DArray *out);

/* CSV: norad_id_a,norad_id_b,tca_s,miss_km,rel_speed_kms. objs maps indices to NORAD ids. */
int conj_write_csv(const DArray *ev, const Object *objs, const char *path);
#endif
