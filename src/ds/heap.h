#ifndef HEAP_H
#define HEAP_H

/* Binary min-heap of (double key, int value) pairs, array based.
 * Decrease-key is done by lazy deletion: push the better entry and let the
 * caller discard stale entries when they are popped. For max-heap behaviour
 * push the negated key. */
typedef struct { double key; int val; } HeapItem;
typedef struct { HeapItem *a; int n, cap; } Heap;

int  heap_init(Heap *h, int cap);                 /* 0 ok, -1 out of memory */
int  heap_push(Heap *h, double key, int val);     /* O(log n); 0 ok, -1 oom */
int  heap_pop(Heap *h, HeapItem *out);            /* O(log n); 1 popped, 0 empty */
int  heap_peek(const Heap *h, HeapItem *out);     /* O(1);     1 ok, 0 empty */
void heap_free(Heap *h);
#endif
