#include <stdlib.h>
#include "dsu.h"

int dsu_init(Dsu *d, int n) {
    d->n = n; d->comps = n;
    d->parent = malloc(sizeof(int) * (size_t)(n > 0 ? n : 1));
    d->size = malloc(sizeof(int) * (size_t)(n > 0 ? n : 1));
    if (!d->parent || !d->size) { dsu_free(d); return -1; }
    for (int i = 0; i < n; i++) { d->parent[i] = i; d->size[i] = 1; }
    return 0;
}

int dsu_find(Dsu *d, int x) {
    while (d->parent[x] != x) { d->parent[x] = d->parent[d->parent[x]]; x = d->parent[x]; }
    return x;
}

int dsu_union(Dsu *d, int a, int b) {
    a = dsu_find(d, a); b = dsu_find(d, b);
    if (a == b) return 0;
    if (d->size[a] < d->size[b]) { int t = a; a = b; b = t; }
    d->parent[b] = a; d->size[a] += d->size[b];
    d->comps--;
    return 1;
}

int dsu_same(Dsu *d, int a, int b) { return dsu_find(d, a) == dsu_find(d, b); }
int dsu_size(Dsu *d, int x) { return d->size[dsu_find(d, x)]; }
void dsu_free(Dsu *d) { free(d->parent); free(d->size); d->parent = d->size = NULL; d->n = d->comps = 0; }
