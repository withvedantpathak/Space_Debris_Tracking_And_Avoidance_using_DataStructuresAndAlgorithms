#ifndef HASH_H
#define HASH_H
#include <stddef.h>
#include <stdint.h>

/* Hash table, separate chaining, int64 key -> int64 value.
 * Average O(1) put/get/remove; table doubles when count exceeds bucket count,
 * keeping the load factor <= 1. Use for NORAD-id lookup and grid cells. */
typedef struct HNode { int64_t key, val; struct HNode *next; } HNode;
typedef struct { HNode **b; size_t nb, n; } Hash;

int    hash_init(Hash *h, size_t buckets);          /* rounded up to a power of 2 */
int    hash_put(Hash *h, int64_t key, int64_t val); /* insert or overwrite; 0 ok, -1 oom */
int    hash_get(const Hash *h, int64_t key, int64_t *val); /* 1 found, 0 missing */
int    hash_remove(Hash *h, int64_t key);           /* 1 removed, 0 missing */
size_t hash_count(const Hash *h);
size_t hash_memory(const Hash *h);                  /* approximate bytes in use */
void   hash_free(Hash *h);
#endif
