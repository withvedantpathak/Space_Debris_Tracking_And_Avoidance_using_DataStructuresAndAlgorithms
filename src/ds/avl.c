#include <stdlib.h>
#include "avl.h"

static int H(const AvlNode *n) { return n ? n->h : 0; }
static void upd(AvlNode *n) { int a = H(n->l), b = H(n->r); n->h = 1 + (a > b ? a : b); }
static int bal(const AvlNode *n) { return H(n->l) - H(n->r); }

static AvlNode *rot_right(AvlNode *y) { AvlNode *x = y->l; y->l = x->r; x->r = y; upd(y); upd(x); return x; }
static AvlNode *rot_left(AvlNode *x)  { AvlNode *y = x->r; x->r = y->l; y->l = x; upd(x); upd(y); return y; }

/* Restore the AVL property at n after an insert or delete below it. O(1). */
static AvlNode *rebalance(AvlNode *n) {
    upd(n);
    int b = bal(n);
    if (b > 1)  { if (bal(n->l) < 0) n->l = rot_left(n->l); return rot_right(n); }
    if (b < -1) { if (bal(n->r) > 0) n->r = rot_right(n->r); return rot_left(n); }
    return n;
}

void avl_init(Avl *t) { t->root = NULL; t->n = 0; t->oom = 0; t->removed = 0; }

static AvlNode *ins(Avl *t, AvlNode *n, int key, void *val) {
    if (!n) {
        AvlNode *x = calloc(1, sizeof(AvlNode));
        if (!x) { t->oom = 1; return NULL; }
        x->key = key; x->val = val; x->h = 1;
        t->n++;
        return x;
    }
    if (key < n->key)      n->l = ins(t, n->l, key, val);
    else if (key > n->key) n->r = ins(t, n->r, key, val);
    else { n->val = val; return n; }
    return rebalance(n);
}

int avl_insert(Avl *t, int key, void *val) {
    t->oom = 0;
    t->root = ins(t, t->root, key, val);
    return t->oom ? -1 : 0;
}

int avl_find(const Avl *t, int key, void **val) {
    const AvlNode *n = t->root;
    while (n) {
        if (key == n->key) { if (val) *val = n->val; return 1; }
        n = key < n->key ? n->l : n->r;
    }
    return 0;
}

static AvlNode *del(Avl *t, AvlNode *n, int key) {
    if (!n) return NULL;
    if (key < n->key)      n->l = del(t, n->l, key);
    else if (key > n->key) n->r = del(t, n->r, key);
    else {
        t->removed = 1;
        if (!n->l || !n->r) {
            AvlNode *c = n->l ? n->l : n->r;
            free(n);
            return c;
        }
        AvlNode *m = n->r;                 /* in-order successor */
        while (m->l) m = m->l;
        n->key = m->key; n->val = m->val;
        n->r = del(t, n->r, m->key);
    }
    return rebalance(n);
}

int avl_remove(Avl *t, int key) {
    t->removed = 0;
    t->root = del(t, t->root, key);
    if (t->removed) t->n--;
    return t->removed;
}

int avl_height(const Avl *t) { return H(t->root); }

/* Returns node count, or -1 if any invariant fails. long long, not long: on Windows (LLP64) long is 32 bits and the
 * exclusive bounds INT_MIN-1 / INT_MAX+1 below would overflow and make every non-empty tree fail the check. */
static long long check(const AvlNode *n, long long lo, long long hi) {
    if (!n) return 0;
    if (n->key <= lo || n->key >= hi) return -1;
    long long a = check(n->l, lo, n->key), b = check(n->r, n->key, hi);
    if (a < 0 || b < 0) return -1;
    int hl = H(n->l), hr = H(n->r);
    if (n->h != 1 + (hl > hr ? hl : hr)) return -1;
    if (hl - hr > 1 || hr - hl > 1) return -1;
    return a + b + 1;
}

int avl_check(const Avl *t) {
    long long c = check(t->root, -2147483647LL - 1 - 1, 2147483647LL + 1);
    return c >= 0 && (size_t)c == t->n;
}

static void walk(const AvlNode *n, void (*fn)(int, void *, void *), void *ctx) {
    if (!n) return;
    walk(n->l, fn, ctx); fn(n->key, n->val, ctx); walk(n->r, fn, ctx);
}
void avl_foreach(const Avl *t, void (*fn)(int, void *, void *), void *ctx) { walk(t->root, fn, ctx); }

static void destroy(AvlNode *n) { if (!n) return; destroy(n->l); destroy(n->r); free(n); }
void avl_free(Avl *t) { destroy(t->root); t->root = NULL; t->n = 0; }
