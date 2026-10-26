#include <stdlib.h>
#include "stack.h"

int stack_init(Stack *s, int cap) {
    if (cap < 4) cap = 4;
    s->a = malloc(sizeof(int) * (size_t)cap);
    s->n = 0; s->cap = s->a ? cap : 0;
    return s->a ? 0 : -1;
}
int stack_push(Stack *s, int v) {
    if (s->n == s->cap) {
        int nc = s->cap ? s->cap * 2 : 8;
        int *p = realloc(s->a, sizeof(int) * (size_t)nc);
        if (!p) return -1;
        s->a = p; s->cap = nc;
    }
    s->a[s->n++] = v;
    return 0;
}
int stack_pop(Stack *s, int *out) { if (!s->n) return 0; *out = s->a[--s->n]; return 1; }
int stack_top(const Stack *s, int *out) { if (!s->n) return 0; *out = s->a[s->n - 1]; return 1; }
int stack_empty(const Stack *s) { return s->n == 0; }
void stack_free(Stack *s) { free(s->a); s->a = NULL; s->n = s->cap = 0; }
