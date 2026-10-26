#ifndef RANK_H
#define RANK_H
/* Top-K by smallest key with a bounded binary heap: O(n log K) time, O(K) space (a full sort is O(n log n)).
 * out[0..count) = indices of the K smallest keys, ascending by (key, index). Returns count (<= K), or -1 on oom. */
int topk_smallest(const double *key, int n, int K, int *out);
#endif
