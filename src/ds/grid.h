#ifndef GRID_H
#define GRID_H
#include <stddef.h>
#include "vec3.h"
#include "hash.h"

/* Uniform grid ("linked cells") over an external point array (not copied).
 * Cell (ix,iy,iz) is packed into one int64 key and looked up in a hash table; points in
 * a cell are chained through next[]. Choose cell size close to the query radius.
 *   build:  O(n) expected
 *   radius: O(27 + k_cells) when r <= cell; (2*ceil(r/cell)+1)^3 cells otherwise
 * Cell indices must fit in 21 bits (|coord/cell| < 2^20). */
typedef struct { const Vec3 *pts; int n; double cell; Hash cells; int *next; } Grid;

int    grid_build(Grid *g, const Vec3 *pts, int n, double cell); /* 0 ok, -1 oom, -2 bad cell/range */
void   grid_radius(const Grid *g, Vec3 q, double r, SpatialVisit fn, void *ctx); /* dist <= r */
size_t grid_memory(const Grid *g);
void   grid_free(Grid *g);
#endif
