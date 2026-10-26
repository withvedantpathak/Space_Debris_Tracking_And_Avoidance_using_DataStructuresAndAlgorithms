# Conjunction graph, clusters, Top-K, mitigation (C6)

Code: `src/graph.[ch]`, `src/mitigate.[ch]`, `src/rank.[ch]`, `src/analyze.[ch]`, `src/ds/stack.[ch]`. CLI: `build/dlscms analyze <events.csv> <movable.csv|-> <K> <out_prefix> <exact_node_limit> [noreduce]`.
Tests: `tests/test_graph.c` (22,452 checks). Evidence: `scripts/run_c6.sh` -> `results/c6_summary.csv`, `c6_<w>_{components,topk,topk_crossing,cover,timeline,verify}.csv`; `scripts/verify_cover.py` (independent MILP); `scripts/plot_c6.py` -> `results/c6_summary.png`.

## Data and definitions
Input: events found by `fastscreen kd` on the real catalogue (SGP4, LEO scope, 10 km, windows of 1, 3 and 6 h from the median epoch). Graph: nodes = objects with at least one conjunction, edge = a pair with a conjunction (the closest event of the pair), weight = miss distance. NORAD id -> node via the hash table.
Movable objects = those in `data/raw/active.csv` (active satellites). Debris cannot manoeuvre. A pair with no movable object cannot be cleared; it is reported as unresolvable and excluded.

## What each structure does
| Task | Structure / algorithm |
|---|---|
| NORAD id -> node | hash table |
| Adjacency | compressed arrays (CSR) |
| Components | BFS (queue), iterative DFS (explicit stack), union-find: all three must agree (tests and every run) |
| Cascade reach | BFS hop distance from the busiest object |
| Top-K closest | bounded binary heap, O(n log K); equals a full sort (tests); second list for crossing encounters only (relative speed > 1 km/s) |
| Conjunction timeline | AVL tree keyed by TCA (ms), in-order walk -> events per 10 minutes |
| Mitigation | minimum vertex cover: greedy (max-degree via heap, lazy deletion, plus a redundancy pass) vs exact (reductions, independent groups, branch and bound with matching lower bound, undo stack) |

## Results (`results/c6_summary.csv`)
| Window | pairs | objects | clusters | largest cluster | greedy | exact (ours) | MILP optimum | lower bound |
|---|---|---|---|---|---|---|---|---|
| 1 h | 11,914 | 8,923 | 717 | 7,197 objects | 3,737 | **3,687** | 3,687 | 3,310 |
| 3 h | 29,229 | 12,005 | 950 | 9,416 | 4,842 | 4,839 (cut off) | 4,832 | 4,342 |
| 6 h | 52,870 | 14,306 | 895 | 11,045 | 5,614 | **5,595** | 5,595 | 4,921 |

- The conjunction graph is mostly one giant cluster (mean degree 2.7 at 1 h and 7.4 at 6 h, both above the threshold of 1 at which a giant component appears in random graphs) plus hundreds of isolated pairs. Cascade reach from the busiest object: the whole giant cluster within 15-18 hops.
- The exact search matches the independent integer-programming optimum on 1 h and 6 h. On 3 h the giant cluster exceeds the node limit (300,000 nodes): the cover is valid and within 7 objects (0.14%) of the optimum, the optimum is not proven by our search.
- Greedy is 0.2-1.4% above the optimum: it is a good heuristic here, but strictly worse than optimal on 1 / 2 / 8 clusters.
- Reductions matter: pure branch and bound (`noreduce`) finishes only clusters up to ~13 objects at 1 h and takes 2 s for what the reduced search does in 0.6 ms (`variant` rows in the summary).
- Top-K by miss distance is dominated by objects drifting together at ~0 km/s (launch batches); `*_topk_crossing.csv` lists real crossing encounters.

## Limits (see `docs/known_limitations.md` L27-L31)
- Vertex cover assumes a manoeuvre removes all conjunctions of that object; real planning needs fuel, timing, new conjunctions created by the move and operator coordination.
- The giant cluster is the hard case for exact search. Viable fix: a stronger lower bound (LP / Nemhauser-Trotter kernel) and per-node reductions; the MILP check shows an exact solution exists and is cheap for a solver, so it is used here as an independent verifier, not as part of the core.
- `exact_us` timings are single runs (`clock()`).
