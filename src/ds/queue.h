#ifndef QUEUE_H
#define QUEUE_H
#include <stddef.h>

/* FIFO queue of ints, singly linked list. Zero-initialise or call queue_init. */
typedef struct QNode { int v; struct QNode *next; } QNode;
typedef struct { QNode *head, *tail; size_t n; } Queue;

void queue_init(Queue *q);
int  queue_push(Queue *q, int v);      /* O(1); 0 ok, -1 oom */
int  queue_pop(Queue *q, int *out);    /* O(1); 1 popped, 0 empty */
int  queue_empty(const Queue *q);
void queue_free(Queue *q);
#endif
