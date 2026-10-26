#ifndef VEC3_H
#define VEC3_H
typedef struct { double x, y, z; } Vec3;

/* Called once per point found by a radius query: point index, squared distance. */
typedef void (*SpatialVisit)(int idx, double dist2, void *ctx);
#endif
