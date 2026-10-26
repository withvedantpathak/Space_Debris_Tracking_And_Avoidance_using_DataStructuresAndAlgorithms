#ifndef DARRAY_H
#define DARRAY_H
#include <stddef.h>

/* Growable array of fixed-size elements (amortised O(1) push). */
typedef struct { void *data; size_t n, cap, esz; } DArray;

void   da_init(DArray *a, size_t elem_size);
int    da_push(DArray *a, const void *elem);   /* 0 ok, -1 out of memory */
void  *da_at(const DArray *a, size_t i);       /* no bounds check; i < n */
void   da_clear(DArray *a);                    /* keeps capacity */
void   da_free(DArray *a);
#endif
