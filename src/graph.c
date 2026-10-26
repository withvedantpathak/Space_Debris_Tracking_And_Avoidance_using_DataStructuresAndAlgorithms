#include <stdlib.h>
#include <string.h>
#include "graph.h"
#include "ds/queue.h"
#include "ds/stack.h"
#include "ds/dsu.h"

int graph_build(Graph *g, int n, int m, const int *eu, const int *ev, const double *w) {
    size_t N = (size_t)(n > 0 ? n : 1), M = (size_t)(m > 0 ? m : 1);
    g->n = n; g->m = m;
    g->eu = malloc(sizeof(int) * M); g->ev = malloc(sizeof(int) * M); g->w = malloc(sizeof(double) * M);
    g->adj_start = calloc(N + 1, sizeof(int)); g->adj_edge = malloc(sizeof(int) * 2 * M);
    if (!g->eu || !g->ev || !g->w || !g->adj_start || !g->adj_edge) { graph_free(g); return -1; }
    for (int e = 0; e < m; e++) {
        if (eu[e] < 0 || eu[e] >= n || ev[e] < 0 || ev[e] >= n || eu[e] == ev[e]) { graph_free(g); return -2; }
        g->eu[e] = eu[e]; g->ev[e] = ev[e]; g->w[e] = w ? w[e] : 0.0;
        g->adj_start[eu[e] + 1]++; g->adj_start[ev[e] + 1]++;
    }
    for (int v = 0; v < n; v++) g->adj_start[v + 1] += g->adj_start[v];
    int *pos = malloc(sizeof(int) * N);
    if (!pos) { graph_free(g); return -1; }
    memcpy(pos, g->adj_start, sizeof(int) * (size_t)n);
    for (int e = 0; e < m; e++) { g->adj_edge[pos[eu[e]]++] = e; g->adj_edge[pos[ev[e]]++] = e; }
    free(pos);
    return 0;
}
int graph_degree(const Graph *g, int v) { return g->adj_start[v + 1] - g->adj_start[v]; }
static int other(const Graph *g, int e, int v) { return g->eu[e] == v ? g->ev[e] : g->eu[e]; }

int graph_components_bfs(const Graph *g, int *comp) {
    for (int v = 0; v < g->n; v++) comp[v] = -1;
    Queue q; queue_init(&q);
    int k = 0;
    for (int s = 0; s < g->n; s++) {
        if (comp[s] != -1) continue;
        comp[s] = k;
        if (queue_push(&q, s)) { queue_free(&q); return -1; }
        int v;
        while (queue_pop(&q, &v)) {
            for (int i = g->adj_start[v]; i < g->adj_start[v + 1]; i++) {
                int w = other(g, g->adj_edge[i], v);
                if (comp[w] == -1) { comp[w] = k; if (queue_push(&q, w)) { queue_free(&q); return -1; } }
            }
        }
        k++;
    }
    queue_free(&q);
    return k;
}

int graph_components_dfs(const Graph *g, int *comp) {
    for (int v = 0; v < g->n; v++) comp[v] = -1;
    Stack st; if (stack_init(&st, 64)) return -1;
    int k = 0;
    for (int s = 0; s < g->n; s++) {
        if (comp[s] != -1) continue;
        comp[s] = k;
        if (stack_push(&st, s)) { stack_free(&st); return -1; }
        int v;
        while (stack_pop(&st, &v)) {
            for (int i = g->adj_start[v]; i < g->adj_start[v + 1]; i++) {
                int w = other(g, g->adj_edge[i], v);
                if (comp[w] == -1) { comp[w] = k; if (stack_push(&st, w)) { stack_free(&st); return -1; } }
            }
        }
        k++;
    }
    stack_free(&st);
    return k;
}

int graph_components_dsu(const Graph *g, int *comp) {
    Dsu d; if (dsu_init(&d, g->n)) return -1;
    for (int e = 0; e < g->m; e++) dsu_union(&d, g->eu[e], g->ev[e]);
    int *label = malloc(sizeof(int) * (size_t)(g->n > 0 ? g->n : 1));
    if (!label) { dsu_free(&d); return -1; }
    for (int v = 0; v < g->n; v++) label[v] = -1;
    int k = 0;
    for (int v = 0; v < g->n; v++) { int r = dsu_find(&d, v); if (label[r] == -1) label[r] = k++; comp[v] = label[r]; }
    free(label); dsu_free(&d);
    return k;
}

int graph_hops(const Graph *g, int src, int *dist) {
    for (int v = 0; v < g->n; v++) dist[v] = -1;
    Queue q; queue_init(&q);
    dist[src] = 0; int reached = 1;
    if (queue_push(&q, src)) return -1;
    int v;
    while (queue_pop(&q, &v)) {
        for (int i = g->adj_start[v]; i < g->adj_start[v + 1]; i++) {
            int w = other(g, g->adj_edge[i], v);
            if (dist[w] == -1) { dist[w] = dist[v] + 1; reached++; if (queue_push(&q, w)) { queue_free(&q); return -1; } }
        }
    }
    queue_free(&q);
    return reached;
}

int graph_subgraph(const Graph *g, const int *comp, int c, Graph *out, int *orig_node) {
    int *local = malloc(sizeof(int) * (size_t)(g->n > 0 ? g->n : 1));
    int *eu = malloc(sizeof(int) * (size_t)(g->m > 0 ? g->m : 1)), *ev = malloc(sizeof(int) * (size_t)(g->m > 0 ? g->m : 1));
    double *w = malloc(sizeof(double) * (size_t)(g->m > 0 ? g->m : 1));
    if (!local || !eu || !ev || !w) { free(local); free(eu); free(ev); free(w); return -1; }
    int k = 0, m = 0;
    for (int v = 0; v < g->n; v++) if (comp[v] == c) { local[v] = k; if (orig_node) orig_node[k] = v; k++; }
    for (int e = 0; e < g->m; e++) if (comp[g->eu[e]] == c) { eu[m] = local[g->eu[e]]; ev[m] = local[g->ev[e]]; w[m] = g->w[e]; m++; }
    int rc = graph_build(out, k, m, eu, ev, w);
    free(local); free(eu); free(ev); free(w);
    return rc;
}

void graph_free(Graph *g) {
    free(g->eu); free(g->ev); free(g->w); free(g->adj_start); free(g->adj_edge);
    g->eu = g->ev = g->adj_start = g->adj_edge = NULL; g->w = NULL; g->n = g->m = 0;
}
