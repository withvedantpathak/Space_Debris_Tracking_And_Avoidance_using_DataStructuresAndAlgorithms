#ifndef GRAPH_H
#define GRAPH_H
/* Undirected conjunction graph: nodes are objects, an edge joins two objects that have a conjunction.
 * Compressed adjacency (CSR): the edges incident to node v are adj_edge[adj_start[v] .. adj_start[v+1]).
 * Edge e joins eu[e] and ev[e] and carries a weight w[e] (the miss distance, km). Edges must join two
 * different nodes; the caller removes duplicate pairs (use conj_best_per_pair). */
typedef struct { int n, m; int *eu, *ev; double *w; int *adj_start, *adj_edge; } Graph;

int  graph_build(Graph *g, int n, int m, const int *eu, const int *ev, const double *w);  /* copies; 0 ok, -1 oom, -2 bad edge */
int  graph_degree(const Graph *g, int v);
/* Connected components, three independent ways. comp[v] = label 0..k-1 numbered by the lowest node of each
 * component; return k, or -1 on out of memory. They must agree (tests). */
int  graph_components_bfs(const Graph *g, int *comp);   /* queue */
int  graph_components_dfs(const Graph *g, int *comp);   /* iterative, explicit stack */
int  graph_components_dsu(const Graph *g, int *comp);   /* disjoint-set union */
/* Hop distance from src to every node by BFS (-1 = unreachable); returns the number of nodes reached (incl. src), -1 on oom.
 * Models how far a cascade of conjunctions can spread. */
int  graph_hops(const Graph *g, int src, int *dist);
/* Subgraph of component c: nodes relabelled 0..k-1 (orig_node[i] = node of g), edges kept with weights. 0 ok, -1 oom. */
int  graph_subgraph(const Graph *g, const int *comp, int c, Graph *out, int *orig_node);
void graph_free(Graph *g);
#endif
