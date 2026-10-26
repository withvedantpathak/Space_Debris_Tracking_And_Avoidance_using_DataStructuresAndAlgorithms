#ifndef STACK_H
#define STACK_H
/* LIFO stack of ints, array based with doubling. Used for iterative DFS and for the undo stack of the
 * exact mitigation search. Zero-initialise or call stack_init. */
typedef struct { int *a; int n, cap; } Stack;
int  stack_init(Stack *s, int cap);     /* 0 ok, -1 oom */
int  stack_push(Stack *s, int v);       /* O(1) amortised; 0 ok, -1 oom */
int  stack_pop(Stack *s, int *out);     /* 1 popped, 0 empty */
int  stack_top(const Stack *s, int *out);   /* 1 ok, 0 empty */
int  stack_empty(const Stack *s);
void stack_free(Stack *s);
#endif
