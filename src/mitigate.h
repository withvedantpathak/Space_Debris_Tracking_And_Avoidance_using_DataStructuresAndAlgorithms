#ifndef MITIGATE_H
#define MITIGATE_H
#include "graph.h"

/* Mitigation planning = minimum vertex cover on the conjunction graph: choose the fewest objects to manoeuvre so that
 * every conjunction involves at least one manoeuvred object (moving an object is assumed to remove all its conjunctions).
 * movable[v] != 0 marks objects that can manoeuvre (active satellites); NULL = all movable. An edge with no movable
 * endpoint cannot be cleared; it is skipped and counted in *unresolved. Minimum vertex cover is NP-hard.
 *
 * in_cover[v] (caller-allocated, g->n bytes) is filled with 1 for chosen objects. All return the cover size, or -1 on oom. */

/* Greedy: forced picks first (edges with exactly one movable endpoint), then repeatedly the movable object that clears
 * the most remaining conjunctions (max-degree, binary heap with lazy deletion), then a redundancy pass that drops
 * any chosen object whose conjunctions are all cleared by the other object too. Not optimal: worst case ratio ~ ln(n). */
int vc_greedy(const Graph *g, const unsigned char *movable, unsigned char *in_cover, int *unresolved);

/* Size of a maximal matching of the clearable edges: a LOWER BOUND on the optimum (disjoint edges need distinct objects). */
int vc_lower_bound(const Graph *g, const unsigned char *movable);

/* Exact: backtracking branch and bound. Branch on the first uncovered edge (try each movable endpoint), undo via an
 * explicit stack, prune with the matching lower bound, start from the greedy solution as upper bound. Stops after
 * node_limit search nodes: *complete = 0 then means in_cover holds the best found so far (still a valid cover), not
 * a proven optimum. *nodes = search nodes visited. */
int vc_exact(const Graph *g, const unsigned char *movable, long node_limit, unsigned char *in_cover, long *nodes, int *complete, int *unresolved);
/* Ablation switch for experiments: 0 disables the leaf/forced reductions in vc_exact (pure branch and bound). Default 1. */
extern int vc_use_reductions;
#endif
