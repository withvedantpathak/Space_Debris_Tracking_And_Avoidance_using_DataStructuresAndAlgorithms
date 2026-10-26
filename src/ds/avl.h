#ifndef AVL_H
#define AVL_H
#include <stddef.h>

/* Self-balancing AVL tree, int key -> void* value (the tree does not own values).
 * Height <= 1.44 log2(n); insert/find/remove are O(log n). */
typedef struct AvlNode { int key, h; void *val; struct AvlNode *l, *r; } AvlNode;
typedef struct { AvlNode *root; size_t n; int oom, removed; } Avl;

void   avl_init(Avl *t);
int    avl_insert(Avl *t, int key, void *val);   /* overwrites existing key; 0 ok, -1 oom */
int    avl_find(const Avl *t, int key, void **val); /* 1 found */
int    avl_remove(Avl *t, int key);              /* 1 removed, 0 missing */
int    avl_height(const Avl *t);
int    avl_check(const Avl *t);                  /* 1 if BST order, heights, balance, count all valid */
void   avl_foreach(const Avl *t, void (*fn)(int key, void *val, void *ctx), void *ctx); /* in order */
void   avl_free(Avl *t);
#endif
