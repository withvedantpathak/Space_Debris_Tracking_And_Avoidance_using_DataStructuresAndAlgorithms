#include <stdlib.h>
#include "heap.h"

int heap_init(Heap *h, int cap) {
    if (cap < 4) cap = 4;
    h->a = malloc(sizeof(HeapItem) * (size_t)cap);
    h->n = 0;
    h->cap = h->a ? cap : 0;
    return h->a ? 0 : -1;
}

int heap_push(Heap *h, double key, int val) {
    if (h->n == h->cap) {
        int nc = h->cap ? h->cap * 2 : 4;
        HeapItem *p = realloc(h->a, sizeof(HeapItem) * (size_t)nc);
        if (!p) return -1;
        h->a = p; h->cap = nc;
    }
    int i = h->n++;                      /* sift up */
    while (i > 0) {
        int p = (i - 1) / 2;
        if (h->a[p].key <= key) break;
        h->a[i] = h->a[p];
        i = p;
    }
    h->a[i].key = key; h->a[i].val = val;
    return 0;
}

int heap_pop(Heap *h, HeapItem *out) {
    if (h->n == 0) return 0;
    HeapItem top = h->a[0], last = h->a[--h->n];
    int i = 0;                           /* sift down */
    for (;;) {
        int c = 2 * i + 1;
        if (c >= h->n) break;
        if (c + 1 < h->n && h->a[c + 1].key < h->a[c].key) c++;
        if (h->a[c].key >= last.key) break;
        h->a[i] = h->a[c];
        i = c;
    }
    if (h->n > 0) h->a[i] = last;
    if (out) *out = top;
    return 1;
}

int heap_peek(const Heap *h, HeapItem *out) {
    if (h->n == 0) return 0;
    *out = h->a[0];
    return 1;
}

void heap_free(Heap *h) { free(h->a); h->a = NULL; h->n = h->cap = 0; }
