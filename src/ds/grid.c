#include <stdlib.h>
#include "grid.h"

#define OFF   (1 << 20)
#define LIMIT (1 << 21)

static int cell_of(double v, double cell) {          /* floor(v/cell) without libm */
    double s = v / cell;
    if (s > 1e9 || s < -1e9) return -OFF - 1;         /* flagged as out of range */
    int i = (int)s;
    return (double)i > s ? i - 1 : i;
}
static int in_range(int c) { return c + OFF >= 0 && c + OFF < LIMIT; }
static int64_t key_of(int ix, int iy, int iz) {
    return ((int64_t)(ix + OFF) << 42) | ((int64_t)(iy + OFF) << 21) | (int64_t)(iz + OFF);
}

int grid_build(Grid *g, const Vec3 *pts, int n, double cell) {
    g->pts = pts; g->n = n; g->cell = cell; g->next = NULL;
    g->cells.b = NULL; g->cells.nb = g->cells.n = 0;      /* so grid_free is safe after ANY failed build */
    if (!(cell > 0)) return -2;
    if (hash_init(&g->cells, (size_t)(n > 8 ? n : 8)) != 0) return -1;
    g->next = malloc(sizeof(int) * (size_t)(n > 0 ? n : 1));
    if (!g->next) { hash_free(&g->cells); return -1; }
    for (int i = 0; i < n; i++) {
        int ix = cell_of(pts[i].x, cell), iy = cell_of(pts[i].y, cell), iz = cell_of(pts[i].z, cell);
        if (!in_range(ix) || !in_range(iy) || !in_range(iz)) { grid_free(g); return -2; }
        int64_t k = key_of(ix, iy, iz), head;
        g->next[i] = hash_get(&g->cells, k, &head) ? (int)head : -1;
        if (hash_put(&g->cells, k, i) != 0) { grid_free(g); return -1; }
    }
    return 0;
}

void grid_radius(const Grid *g, Vec3 q, double r, SpatialVisit fn, void *ctx) {
    if (g->n <= 0 || r < 0) return;
    double r2 = r * r;
    int span = (int)(r / g->cell) + 1;                /* safe superset of ceil(r/cell) */
    int cx = cell_of(q.x, g->cell), cy = cell_of(q.y, g->cell), cz = cell_of(q.z, g->cell);
    for (int ix = cx - span; ix <= cx + span; ix++) {
        if (!in_range(ix)) continue;
        for (int iy = cy - span; iy <= cy + span; iy++) {
            if (!in_range(iy)) continue;
            for (int iz = cz - span; iz <= cz + span; iz++) {
                if (!in_range(iz)) continue;
                int64_t head;
                if (!hash_get(&g->cells, key_of(ix, iy, iz), &head)) continue;
                for (int i = (int)head; i >= 0; i = g->next[i]) {
                    double dx = g->pts[i].x - q.x, dy = g->pts[i].y - q.y, dz = g->pts[i].z - q.z;
                    double d2 = dx * dx + dy * dy + dz * dz;
                    if (d2 <= r2) fn(i, d2, ctx);
                }
            }
        }
    }
}

size_t grid_memory(const Grid *g) { return hash_memory(&g->cells) + (size_t)g->n * sizeof(int); }
void grid_free(Grid *g) {
    if (g->cells.b) hash_free(&g->cells);
    free(g->next); g->next = NULL; g->n = 0;
}
