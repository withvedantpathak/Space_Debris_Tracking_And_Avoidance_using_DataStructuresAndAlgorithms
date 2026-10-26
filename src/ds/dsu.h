#ifndef DSU_H
#define DSU_H

/* Disjoint-set union: union by size + path halving, near-constant amortised time
 * (inverse Ackermann). Used to cluster objects linked by conjunctions. */
typedef struct { int *parent, *size; int n, comps; } Dsu;

int  dsu_init(Dsu *d, int n);            /* 0 ok, -1 oom */
int  dsu_find(Dsu *d, int x);
int  dsu_union(Dsu *d, int a, int b);    /* 1 if two sets were merged, 0 if already joined */
int  dsu_same(Dsu *d, int a, int b);
int  dsu_size(Dsu *d, int x);            /* size of x's set */
void dsu_free(Dsu *d);
#endif
