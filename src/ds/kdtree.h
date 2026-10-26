#ifndef KDTREE_H
#define KDTREE_H
#include <stddef.h>
#include "vec3.h"

/* Static 3-d tree over an external point array (not copied; must outlive the tree).
 * Implicit layout: the node of a range [lo,hi) sits at mid=(lo+hi)/2 of a permutation
 * array, split on the axis with the widest extent.
 *   build:   O(n log n) expected (quickselect median at each level)
 *   radius:  O(n^(2/3) + k) worst case, ~O(log n + k) on spread-out data
 *   nearest: O(log n) expected
 * Rebuild per time step; no insert/delete. */
typedef struct { const Vec3 *pts; int n; int *idx; unsigned char *axis; } KDTree;

int    kd_build(KDTree *t, const Vec3 *pts, int n);   /* 0 ok, -1 oom */
void   kd_radius(const KDTree *t, Vec3 q, double r, SpatialVisit fn, void *ctx); /* dist <= r */
int    kd_nearest(const KDTree *t, Vec3 q, int exclude, double *dist2); /* index or -1; exclude = point index to skip, or -1 */
size_t kd_memory(const KDTree *t);
void   kd_free(KDTree *t);
#endif
