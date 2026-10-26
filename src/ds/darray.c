#include <stdlib.h>
#include <string.h>
#include "darray.h"

void da_init(DArray *a, size_t elem_size) {
    a->data = NULL; a->n = 0; a->cap = 0; a->esz = elem_size;
}

/* Amortised O(1): capacity doubles when full. */
int da_push(DArray *a, const void *elem) {
    if (a->n == a->cap) {
        size_t nc = a->cap ? a->cap * 2 : 8;
        void *p = realloc(a->data, nc * a->esz);
        if (!p) return -1;
        a->data = p; a->cap = nc;
    }
    memcpy((char *)a->data + a->n * a->esz, elem, a->esz);
    a->n++;
    return 0;
}

void *da_at(const DArray *a, size_t i) { return (char *)a->data + i * a->esz; }
void da_clear(DArray *a) { a->n = 0; }
void da_free(DArray *a) { free(a->data); a->data = NULL; a->n = a->cap = 0; }
