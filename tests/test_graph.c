/* Tests for C6: stack, conjunction graph, components, hop distances, Top-K and vertex-cover mitigation. */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "graph.h"
#include "mitigate.h"
#include "rank.h"
#include "ds/stack.h"
#include "analyze.h"

static long checks = 0, fails = 0;
#define CHECK(c) do { checks++; if (!(c)) { fails++; printf("FAIL %s:%d: %s\n", __FILE__, __LINE__, #c); } } while (0)

static unsigned long long rs = 88172645463325252ULL;
static unsigned rnd(void) { rs ^= rs << 13; rs ^= rs >> 7; rs ^= rs << 17; return (unsigned)(rs >> 11); }

static void test_stack(void) {
    Stack s; CHECK(stack_init(&s, 2) == 0);
    int v;
    CHECK(stack_empty(&s) && !stack_pop(&s, &v) && !stack_top(&s, &v));
    for (int i = 0; i < 1000; i++) CHECK(stack_push(&s, i * 3) == 0);
    CHECK(s.n == 1000 && stack_top(&s, &v) && v == 2997);
    for (int i = 999; i >= 0; i--) { CHECK(stack_pop(&s, &v) && v == i * 3); }
    CHECK(stack_empty(&s));
    stack_free(&s);
    Stack z; memset(&z, 0, sizeof z); CHECK(stack_push(&z, 7) == 0 && stack_pop(&z, &v) && v == 7); stack_free(&z);   /* zero-initialised works */
}

/* random simple graph with n nodes, edge probability p/100 */
static int rand_graph(Graph *g, int n, int p, int **eu_out, int **ev_out) {
    int *eu = malloc(sizeof(int) * (size_t)(n * n + 1)), *ev = malloc(sizeof(int) * (size_t)(n * n + 1)); double *w = malloc(sizeof(double) * (size_t)(n * n + 1)); int m = 0;
    for (int a = 0; a < n; a++) for (int b = a + 1; b < n; b++) if ((int)(rnd() % 100) < p) { eu[m] = a; ev[m] = b; w[m] = (double)(rnd() % 1000) / 10.0; m++; }
    CHECK(graph_build(g, n, m, eu, ev, w) == 0);
    if (eu_out) *eu_out = eu; else free(eu);
    if (ev_out) *ev_out = ev; else free(ev);
    free(w);
    return m;
}

static void test_graph_basic(void) {
    Graph g; int eu[3] = {0, 1, 3}, ev[3] = {1, 2, 4};
    CHECK(graph_build(&g, 6, 3, eu, ev, NULL) == 0);
    CHECK(graph_degree(&g, 1) == 2 && graph_degree(&g, 5) == 0 && graph_degree(&g, 0) == 1);
    int comp[6]; CHECK(graph_components_bfs(&g, comp) == 3 && comp[0] == 0 && comp[1] == 0 && comp[2] == 0 && comp[3] == 1 && comp[4] == 1 && comp[5] == 2);
    graph_free(&g);
    int bad_u[1] = {0}, bad_v[1] = {0}; CHECK(graph_build(&g, 3, 1, bad_u, bad_v, NULL) == -2);       /* self loop */
    bad_v[0] = 5; CHECK(graph_build(&g, 3, 1, bad_u, bad_v, NULL) == -2);                              /* node out of range */
    CHECK(graph_build(&g, 0, 0, NULL, NULL, NULL) == 0); CHECK(graph_components_bfs(&g, comp) == 0); graph_free(&g);
}

static void test_components_and_hops(void) {
    for (int t = 0; t < 300; t++) {
        int n = 1 + (int)(rnd() % 14), p = (int)(rnd() % 40);
        Graph g; int *eu, *ev; int m = rand_graph(&g, n, p, &eu, &ev);
        int c1[16], c2[16], c3[16];
        int k1 = graph_components_bfs(&g, c1), k2 = graph_components_dfs(&g, c2), k3 = graph_components_dsu(&g, c3);
        CHECK(k1 == k2 && k2 == k3 && memcmp(c1, c2, sizeof(int) * (size_t)n) == 0 && memcmp(c1, c3, sizeof(int) * (size_t)n) == 0);
        /* Floyd-Warshall reference */
        int D[16][16];
        for (int a = 0; a < n; a++) for (int b = 0; b < n; b++) D[a][b] = a == b ? 0 : 1000;
        for (int e = 0; e < m; e++) D[eu[e]][ev[e]] = D[ev[e]][eu[e]] = 1;
        for (int x = 0; x < n; x++) for (int a = 0; a < n; a++) for (int b = 0; b < n; b++) if (D[a][x] + D[x][b] < D[a][b]) D[a][b] = D[a][x] + D[x][b];
        for (int a = 0; a < n; a++) { CHECK((D[a][0] < 1000) == (c1[a] == c1[0])); }
        int src = (int)(rnd() % (unsigned)n), dist[16], reached = graph_hops(&g, src, dist), cnt = 0;
        for (int b = 0; b < n; b++) { int ref = D[src][b] >= 1000 ? -1 : D[src][b]; CHECK(dist[b] == ref); if (ref >= 0) cnt++; }
        CHECK(reached == cnt);
        /* subgraph keeps exactly the edges of the component */
        int comp_edges = 0; for (int e = 0; e < m; e++) if (c1[eu[e]] == c1[src]) comp_edges++;
        Graph sg; int orig[16]; CHECK(graph_subgraph(&g, c1, c1[src], &sg, orig) == 0); CHECK(sg.m == comp_edges && sg.n == cnt);
        for (int e = 0; e < sg.m; e++) CHECK(c1[orig[sg.eu[e]]] == c1[src]);
        graph_free(&sg); graph_free(&g); free(eu); free(ev);
    }
}

static void test_topk(void) {
    for (int t = 0; t < 200; t++) {
        int n = (int)(rnd() % 60), K = (int)(rnd() % 20);
        double *key = malloc(sizeof(double) * (size_t)(n + 1)); int *out = malloc(sizeof(int) * (size_t)(K + 1));
        for (int i = 0; i < n; i++) key[i] = (double)(rnd() % 1000) / 8.0;                       /* many ties on purpose */
        int c = topk_smallest(key, n, K, out);
        int want = K < n ? K : n; CHECK(c == want);
        /* reference: selection by repeated minimum, ties by index */
        unsigned char *used = calloc((size_t)(n + 1), 1);
        for (int r = 0; r < want; r++) { int b = -1; for (int i = 0; i < n; i++) if (!used[i] && (b < 0 || key[i] < key[b])) b = i; used[b] = 1; CHECK(out[r] == b); }
        free(used); free(key); free(out);
    }
    int o[1]; CHECK(topk_smallest(NULL, 0, 5, o) == 0);
}

static int valid_cover(const Graph *g, const unsigned char *mov, const unsigned char *cov) {
    for (int e = 0; e < g->m; e++) {
        int cl = (!mov || mov[g->eu[e]] || mov[g->ev[e]]);
        if (cl && !cov[g->eu[e]] && !cov[g->ev[e]]) return 0;
    }
    if (mov) for (int v = 0; v < g->n; v++) if (cov[v] && !mov[v]) return 0;       /* never picks an immovable object */
    return 1;
}
static int csize(const Graph *g, const unsigned char *c) { int s = 0; for (int v = 0; v < g->n; v++) s += c[v] != 0; return s; }
/* brute-force optimum by subset enumeration over movable vertices */
static int brute_vc(const Graph *g, const unsigned char *mov) {
    int best = 1000;
    for (unsigned mask = 0; mask < (1u << g->n); mask++) {
        int ok = 1, sz = 0;
        for (int v = 0; v < g->n; v++) if (mask >> v & 1) { if (mov && !mov[v]) { ok = 0; break; } sz++; }
        if (!ok || sz >= best) continue;
        for (int e = 0; e < g->m && ok; e++) {
            int cl = (!mov || mov[g->eu[e]] || mov[g->ev[e]]);
            if (cl && !(mask >> g->eu[e] & 1) && !(mask >> g->ev[e] & 1)) ok = 0;
        }
        if (ok) best = sz;
    }
    return best;
}

static void test_cover_known(void) {
    struct { int n, m; int eu[12], ev[12]; int opt; } K[] = {
        {3, 3, {0, 1, 0}, {1, 2, 2}, 2},                                      /* triangle */
        {5, 4, {0, 0, 0, 0}, {1, 2, 3, 4}, 1},                                /* star */
        {5, 4, {0, 1, 2, 3}, {1, 2, 3, 4}, 2},                                /* path of 5 */
        {4, 6, {0, 0, 0, 1, 1, 2}, {1, 2, 3, 2, 3, 3}, 3},                    /* K4 */
        {6, 9, {0, 0, 0, 1, 1, 1, 2, 2, 2}, {3, 4, 5, 3, 4, 5, 3, 4, 5}, 3},  /* K3,3 */
        {5, 5, {0, 1, 2, 3, 4}, {1, 2, 3, 4, 0}, 3},                          /* 5-cycle */
    };
    for (unsigned i = 0; i < sizeof K / sizeof K[0]; i++) {
        Graph g; CHECK(graph_build(&g, K[i].n, K[i].m, K[i].eu, K[i].ev, NULL) == 0);
        unsigned char cg[16], cx[16]; long nodes; int complete, unres;
        int gs = vc_greedy(&g, NULL, cg, &unres), xs = vc_exact(&g, NULL, 100000, cx, &nodes, &complete, &unres);
        CHECK(xs == K[i].opt && complete && valid_cover(&g, NULL, cx) && csize(&g, cx) == xs);
        CHECK(valid_cover(&g, NULL, cg) && csize(&g, cg) == gs && gs >= xs);
        CHECK(vc_lower_bound(&g, NULL) <= xs);
        graph_free(&g);
    }
    /* empty graph */
    Graph g; graph_build(&g, 4, 0, NULL, NULL, NULL); unsigned char c[4]; long nd; int comp, un;
    CHECK(vc_greedy(&g, NULL, c, &un) == 0 && vc_exact(&g, NULL, 10, c, &nd, &comp, &un) == 0 && vc_lower_bound(&g, NULL) == 0); graph_free(&g);
}

static void test_cover_random(void) {
    int greedy_worse = 0, incomplete_seen = 0;
    for (int t = 0; t < 1500; t++) {
        int n = 2 + (int)(rnd() % 11), p = 10 + (int)(rnd() % 50);
        Graph g; rand_graph(&g, n, p, NULL, NULL);
        unsigned char mov[16], *movp = NULL;
        if (rnd() % 2) { for (int v = 0; v < n; v++) mov[v] = rnd() % 4 != 0; movp = mov; }
        int unres_ref = 0; for (int e = 0; e < g.m; e++) if (movp && !movp[g.eu[e]] && !movp[g.ev[e]]) unres_ref++;
        int opt = brute_vc(&g, movp);
        unsigned char cg[16], cx[16]; long nodes; int complete, u1 = -1, u2 = -1;
        int gs = vc_greedy(&g, movp, cg, &u1), xs = vc_exact(&g, movp, 1000000, cx, &nodes, &complete, &u2);
        CHECK(u1 == unres_ref && u2 == unres_ref);
        CHECK(valid_cover(&g, movp, cg) && csize(&g, cg) == gs);
        CHECK(valid_cover(&g, movp, cx) && csize(&g, cx) == xs);
        CHECK(complete && xs == opt);
        CHECK(gs >= opt && vc_lower_bound(&g, movp) <= opt);
        if (gs > opt) greedy_worse++;
        /* with a tiny budget the answer must still be a valid cover no smaller than the optimum */
        int xs2 = vc_exact(&g, movp, 2, cx, &nodes, &complete, &u2);
        CHECK(valid_cover(&g, movp, cx) && csize(&g, cx) == xs2 && xs2 >= opt);
        if (!complete) incomplete_seen++;
        graph_free(&g);
    }
    CHECK(greedy_worse > 0);        /* the greedy heuristic is demonstrably suboptimal on some graphs */
    CHECK(incomplete_seen > 0);     /* and the node limit really cuts searches short */
    /* larger sparse graph (reductions + groups): exact must still be valid and never worse than greedy */
    Graph g; rand_graph(&g, 400, 1, NULL, NULL);
    unsigned char *cg = malloc(400), *cx = malloc(400); long nodes; int complete, u;
    int gs = vc_greedy(&g, NULL, cg, &u), xs = vc_exact(&g, NULL, 200000, cx, &nodes, &complete, &u);
    CHECK(valid_cover(&g, NULL, cg) && valid_cover(&g, NULL, cx) && xs <= gs && xs >= vc_lower_bound(&g, NULL) && csize(&g, cx) == xs);
    graph_free(&g); free(cg); free(cx);
}

/* end-to-end driver on a hand-made event file: triangle {10,11,12}, star {20 centre; 21,22,23}, single pair {30,31}, duplicate events of a pair */
static void rd(char *l, size_t n, FILE *f) { if (!fgets(l, (int)n, f)) l[0] = 0; }
static void test_analyze(void) {
    FILE *f = fopen("build/_t_events.csv", "w"); CHECK(f != NULL); if (!f) return;
    fprintf(f, "norad_id_a,norad_id_b,tca_s,miss_km,rel_speed_kms\n");
    fprintf(f, "10,11,100.0,3.0,7.0\n10,12,200.0,2.0,7.0\n11,12,300.0,1.0,7.0\n20,21,400.0,4.0,0.5\n20,22,500.0,5.0,7.0\n20,23,600.0,6.0,7.0\n30,31,700.0,0.5,7.0\n30,31,900.0,0.7,7.0\n");
    fclose(f);
    CHECK(analyze_events("build/_t_events.csv", NULL, 3, "build/_t_c6", 100000) == 0);
    FILE *g = fopen("build/_t_c6_cover.csv", "r"); CHECK(g != NULL);
    int total = 0, in_best = 0, in_greedy = 0; char line[256];
    if (g) { rd(line, sizeof line, g); while (fgets(line, sizeof line, g)) { int id, a, b; if (sscanf(line, "%d,%d,%d", &id, &a, &b) == 3) { total++; in_best += b; in_greedy += a; } } fclose(g); }
    CHECK(in_best == 2 + 1 + 1 && in_greedy == 4);                       /* triangle 2, star 1 (the centre), pair 1 */
    g = fopen("build/_t_c6_topk.csv", "r"); CHECK(g != NULL);
    if (g) { rd(line, sizeof line, g); rd(line, sizeof line, g); CHECK(strstr(line, "1,30,31,700.000,0.5000") != NULL); fclose(g); }   /* best event of the repeated pair wins rank 1 */
    g = fopen("build/_t_c6_topk_crossing.csv", "r"); CHECK(g != NULL);
    if (g) { int rows = 0; rd(line, sizeof line, g); while (fgets(line, sizeof line, g)) { rows++; CHECK(strstr(line, ",0.500,") == NULL || strstr(line, ",7.000,") != NULL); } CHECK(rows == 3); fclose(g); }
    g = fopen("build/_t_c6_timeline.csv", "r"); CHECK(g != NULL);
    if (g) { int bins = 0, sum = 0; rd(line, sizeof line, g); while (fgets(line, sizeof line, g)) { int mn, c; if (sscanf(line, "%d,%d", &mn, &c) == 2) { bins++; sum += c; } } CHECK(bins >= 1 && sum == 7); fclose(g); }   /* 7 pairs */
    CHECK(analyze_events("build/does_not_exist.csv", NULL, 3, "build/_t_c6", 10) != 0);
}

int main(void) {
    test_stack(); test_graph_basic(); test_components_and_hops(); test_topk(); test_cover_known(); test_cover_random(); test_analyze();
    printf("%ld checks, %ld failed\n", checks, fails);
    return fails != 0;
}
