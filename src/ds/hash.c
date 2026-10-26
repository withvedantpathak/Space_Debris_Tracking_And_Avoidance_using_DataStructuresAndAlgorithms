#include <stdlib.h>
#include "hash.h"

/* splitmix64 finaliser: good avalanche, so power-of-two masking is safe. */
static uint64_t mix(uint64_t x) {
    x += 0x9e3779b97f4a7c15ULL;
    x = (x ^ (x >> 30)) * 0xbf58476d1ce4e5b9ULL;
    x = (x ^ (x >> 27)) * 0x94d049bb133111ebULL;
    return x ^ (x >> 31);
}

int hash_init(Hash *h, size_t buckets) {
    size_t nb = 8;
    while (nb < buckets) nb <<= 1;
    h->b = calloc(nb, sizeof(HNode *));
    h->nb = h->b ? nb : 0;
    h->n = 0;
    return h->b ? 0 : -1;
}

/* O(n) rehash, called rarely (doubling) so puts stay amortised O(1). */
static void grow(Hash *h) {
    size_t nb = h->nb * 2;
    HNode **nbk = calloc(nb, sizeof(HNode *));
    if (!nbk) return;                     /* keep working with longer chains */
    for (size_t i = 0; i < h->nb; i++) {
        HNode *n = h->b[i];
        while (n) {
            HNode *nx = n->next;
            size_t j = mix((uint64_t)n->key) & (nb - 1);
            n->next = nbk[j]; nbk[j] = n;
            n = nx;
        }
    }
    free(h->b);
    h->b = nbk; h->nb = nb;
}

int hash_put(Hash *h, int64_t key, int64_t val) {
    if (!h->b) return -1;                 /* hash_init failed (or was never called): report out of memory instead of crashing */
    size_t i = mix((uint64_t)key) & (h->nb - 1);
    for (HNode *n = h->b[i]; n; n = n->next)
        if (n->key == key) { n->val = val; return 0; }
    HNode *n = malloc(sizeof(HNode));
    if (!n) return -1;
    n->key = key; n->val = val; n->next = h->b[i]; h->b[i] = n;
    h->n++;
    if (h->n > h->nb) grow(h);
    return 0;
}

int hash_get(const Hash *h, int64_t key, int64_t *val) {
    if (!h->b) return 0;
    for (HNode *n = h->b[mix((uint64_t)key) & (h->nb - 1)]; n; n = n->next)
        if (n->key == key) { if (val) *val = n->val; return 1; }
    return 0;
}

int hash_remove(Hash *h, int64_t key) {
    if (!h->b) return 0;
    HNode **pp = &h->b[mix((uint64_t)key) & (h->nb - 1)];
    while (*pp) {
        if ((*pp)->key == key) {
            HNode *dead = *pp;
            *pp = dead->next;
            free(dead);
            h->n--;
            return 1;
        }
        pp = &(*pp)->next;
    }
    return 0;
}

size_t hash_count(const Hash *h) { return h->n; }
size_t hash_memory(const Hash *h) { return h->nb * sizeof(HNode *) + h->n * sizeof(HNode); }

void hash_free(Hash *h) {
    for (size_t i = 0; i < h->nb; i++) {
        HNode *n = h->b[i];
        while (n) { HNode *nx = n->next; free(n); n = nx; }
    }
    free(h->b);
    h->b = NULL; h->nb = h->n = 0;
}
