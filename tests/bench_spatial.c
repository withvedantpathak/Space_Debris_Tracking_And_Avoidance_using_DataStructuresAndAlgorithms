/* C1 micro-benchmark: fixed-radius self-join (find all pairs closer than R) on synthetic
 * points in a thin shell (radius 6600-8000, orbit-like), using brute force, uniform grid
 * and k-d tree. Output: CSV on stdout. NOT the conjunction pipeline: no propagation yet.
 * Usage: bench_spatial [radius_km]   (default 50) */
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include "ds/kdtree.h"
#include "ds/grid.h"

static unsigned long long rs;
static double runif(double lo, double hi) { rs ^= rs << 13; rs ^= rs >> 7; rs ^= rs << 17; return lo + (hi - lo) * ((double)(rs >> 11) / 9007199254740992.0); }
static double ms(clock_t a, clock_t b) { return 1000.0 * (double)(b - a) / CLOCKS_PER_SEC; }

typedef struct { long pairs; int self; } Cnt;
static void count_visit(int idx, double d2, void *ctx) { (void)d2; Cnt *c = ctx; if (idx > c->self) c->pairs++; }

static void gen(Vec3 *p, int n) {
    for (int i = 0; i < n; i++) {
        double dx, dy, dz, n2;
        do { dx = runif(-1, 1); dy = runif(-1, 1); dz = runif(-1, 1); n2 = dx * dx + dy * dy + dz * dz; } while (n2 > 1 || n2 < 1e-6);
        double g = 1; for (int k = 0; k < 40; k++) g = 0.5 * (g + n2 / g);
        double r = runif(6600, 8000);
        p[i].x = r * dx / g; p[i].y = r * dy / g; p[i].z = r * dz / g;
    }
}

int main(int argc, char **argv) {
    double R = argc > 1 ? atof(argv[1]) : 50.0;
    int sizes[] = {1000, 2000, 5000, 10000, 20000, 50000};
    printf("radius_km,n,seed,method,build_ms,query_ms,pairs,recall,mem_bytes\n");
    for (int si = 0; si < 6; si++) for (int seed = 1; seed <= 3; seed++) {
        int n = sizes[si];
        rs = 88172645463325252ULL + (unsigned long long)seed * 7919ULL;
        Vec3 *p = malloc(sizeof(Vec3) * (size_t)n);
        gen(p, n);
        clock_t t0 = clock();
        long bf = 0; double r2 = R * R;
        for (int i = 0; i < n; i++) for (int j = i + 1; j < n; j++) {
            double dx = p[i].x - p[j].x, dy = p[i].y - p[j].y, dz = p[i].z - p[j].z;
            if (dx * dx + dy * dy + dz * dz <= r2) bf++;
        }
        clock_t t1 = clock();
        printf("%g,%d,%d,brute,0,%.3f,%ld,1.0,0\n", R, n, seed, ms(t0, t1), bf);

        Grid g; t0 = clock(); if (grid_build(&g, p, n, R) != 0) return 1; t1 = clock();
        Cnt c; c.pairs = 0;
        for (int i = 0; i < n; i++) { c.self = i; grid_radius(&g, p[i], R, count_visit, &c); }
        clock_t t2 = clock();
        printf("%g,%d,%d,grid,%.3f,%.3f,%ld,%.6f,%zu\n", R, n, seed, ms(t0, t1), ms(t1, t2), c.pairs, bf ? (double)c.pairs / (double)bf : 1.0, grid_memory(&g));
        grid_free(&g);

        KDTree t; t0 = clock(); if (kd_build(&t, p, n) != 0) return 1; t1 = clock();
        c.pairs = 0;
        for (int i = 0; i < n; i++) { c.self = i; kd_radius(&t, p[i], R, count_visit, &c); }
        t2 = clock();
        printf("%g,%d,%d,kdtree,%.3f,%.3f,%ld,%.6f,%zu\n", R, n, seed, ms(t0, t1), ms(t1, t2), c.pairs, bf ? (double)c.pairs / (double)bf : 1.0, kd_memory(&t));
        kd_free(&t); free(p);
    }
    return 0;
}
