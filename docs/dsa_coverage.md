# DSA concept coverage (where each concept is used and tested)

| Concept | Where it is used | Evidence |
|---|---|---|
| Arrays / dynamic arrays | `src/ds/darray`, event lists, CSR graph | `tests/test_ds.c` |
| Linked lists | queue nodes; hash-table chaining; grid cell chains | test_ds |
| Stack | iterative DFS components; undo stack of the exact vertex-cover search | `tests/test_graph.c` |
| Queue | BFS components, hop distances | test_graph |
| Binary heap | Top-K closest (bounded), greedy max-degree with lazy deletion | test_ds, test_graph |
| Hash table | NORAD id lookup, grid cells | test_ds, test_catalog |
| BST / AVL | time-ordered conjunction index (`analyze`); k-d tree is a spatial BST | test_ds (AVL), test_graph (analyze) |
| Sorting | quickselect (k-d tree median), qsort of events | test_ds, test_conjunction |
| Searching | binary-style k-d tree and grid radius queries | test_ds |
| Graphs | conjunction graph (CSR), subgraphs | test_graph |
| BFS / DFS | components, hop distance | test_graph |
| Union-find | components cross-check | test_ds, test_graph |
| Greedy | mitigation heuristic | test_graph, `results/c6_summary.csv` |
| Backtracking / branch and bound | exact vertex cover | test_graph (vs brute force), MILP check |
| Complexity analysis | C1 benchmark, C5 benchmark (O(n^2) vs ~O(n log n)), C6 exact search cost | `results/bench_c1.csv`, `bench_c5.csv`, `c6_summary.csv` |
| Dynamic programming / shortest path | not used (hop distances use BFS, unweighted) | - |
