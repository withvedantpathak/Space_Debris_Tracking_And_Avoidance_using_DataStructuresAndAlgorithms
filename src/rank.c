#include <stdlib.h>
#include "rank.h"
#include "ds/heap.h"

static const double *g_key;
static int cmp_idx(const void *a, const void *b) {
    int i = *(const int *)a, j = *(const int *)b;
    if (g_key[i] != g_key[j]) return g_key[i] < g_key[j] ? -1 : 1;
    return (i > j) - (i < j);
}

int topk_smallest(const double *key, int n, int K, int *out) {
    if (K <= 0 || n <= 0) return 0;
    if (K > n) K = n;                                  /* never keep more than there are (also avoids K + 1 overflowing for K = INT_MAX) */
    Heap h; if (heap_init(&h, K + 1)) return -1;      /* min-heap on -key: the root is the LARGEST key kept */
    for (int i = 0; i < n; i++) {
        if (h.n == K) {                                    /* full: a new key must beat the largest kept one */
            HeapItem top; heap_peek(&h, &top);
            if (key[i] >= -top.key) continue;              /* equal keys: the lower index (already kept) wins */
            heap_pop(&h, &top);
        }
        if (heap_push(&h, -key[i], i)) { heap_free(&h); return -1; }
    }
    /* the heap holds a set of K smallest keys but may have kept the wrong ones among keys equal to the K-th smallest:
       fix that with one pass (all keys below the threshold, then the lowest indices equal to it) */
    HeapItem top; heap_peek(&h, &top);
    double thr = -top.key; int below = 0, c = 0;
    for (int i = 0; i < n; i++) if (key[i] < thr) below++;
    int ties = (h.n < K ? h.n : K) - below;
    for (int i = 0; i < n; i++) {
        if (key[i] < thr) out[c++] = i;
        else if (key[i] == thr && ties > 0) { out[c++] = i; ties--; }
    }
    heap_free(&h);
    g_key = key;
    qsort(out, (size_t)c, sizeof(int), cmp_idx);       /* deterministic order, ties by index */
    return c;
}
