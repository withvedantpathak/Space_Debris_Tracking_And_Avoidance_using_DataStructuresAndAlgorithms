#include <stdlib.h>
#include "queue.h"

void queue_init(Queue *q) { q->head = q->tail = NULL; q->n = 0; }

int queue_push(Queue *q, int v) {
    QNode *n = malloc(sizeof(QNode));
    if (!n) return -1;
    n->v = v; n->next = NULL;
    if (q->tail) q->tail->next = n; else q->head = n;
    q->tail = n;
    q->n++;
    return 0;
}

int queue_pop(Queue *q, int *out) {
    QNode *n = q->head;
    if (!n) return 0;
    if (out) *out = n->v;
    q->head = n->next;
    if (!q->head) q->tail = NULL;
    free(n);
    q->n--;
    return 1;
}

int queue_empty(const Queue *q) { return q->head == NULL; }

void queue_free(Queue *q) { while (queue_pop(q, NULL)) {} }
