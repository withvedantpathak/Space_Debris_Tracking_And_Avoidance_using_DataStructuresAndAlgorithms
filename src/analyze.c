#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include "analyze.h"
#include "catalog.h"
#include "graph.h"
#include "mitigate.h"
#include "rank.h"
#include "ds/hash.h"
#include "ds/darray.h"
#include "ds/avl.h"

typedef struct { int id_a, id_b; double tca, miss, vrel; } Ev;

#define TIMELINE_MAX_S 2.0e6
typedef struct { int *bins; int nb; } Timeline;
static void timeline_visit(int key, void *val, void *ctx) {     /* in-order walk of the AVL tree: keys arrive sorted by TCA */
    (void)val;
    Timeline *t = ctx; int b = key / 600000;                     /* key = TCA in ms; 10-minute bins */
    if (b >= 0 && b < t->nb) t->bins[b]++;
}
static double now_s(void) { return (double)clock() / CLOCKS_PER_SEC; }

static int read_events(const char *path, DArray *ev) {
    FILE *f = fopen(path, "r"); if (!f) return -1;
    char line[512]; int first = 1;
    while (fgets(line, sizeof line, f)) {
        if (first) { first = 0; continue; }
        Ev e;
        if (sscanf(line, "%d,%d,%lf,%lf,%lf", &e.id_a, &e.id_b, &e.tca, &e.miss, &e.vrel) != 5 || e.id_a == e.id_b) continue;   /* malformed line or an object paired with itself: ignored */
        if (da_push(ev, &e)) { fclose(f); return -3; }
    }
    fclose(f);
    return 0;
}

int analyze_events(const char *events_csv, const char *movable_csv, int K, const char *prefix, long exact_node_limit) {
    DArray raw; da_init(&raw, sizeof(Ev));
    if (read_events(events_csv, &raw)) { fprintf(stderr, "cannot read %s\n", events_csv); da_free(&raw); return 1; }
    /* best (smallest miss) event of each pair; the file is sorted by (a, b, tca) */
    const Ev *R = raw.data; DArray best; da_init(&best, sizeof(Ev));
    for (size_t i = 0; i < raw.n; ) {
        size_t b = i, j = i;
        while (j < raw.n && R[j].id_a == R[i].id_a && R[j].id_b == R[i].id_b) { if (R[j].miss < R[b].miss) b = j; j++; }
        if (da_push(&best, &R[b])) { fprintf(stderr, "out of memory\n"); da_free(&raw); da_free(&best); return 1; }
        i = j;
    }
    const Ev *B = best.data; int m = (int)best.n;
    if (K < 0) K = 0;
    if (K > m) K = m;                                   /* Top-K lists cannot be longer than the number of pairs */
    /* NORAD id -> node index */
    Hash idx; if (hash_init(&idx, 1024) != 0) { da_free(&raw); da_free(&best); return 1; }
    int *ids = malloc(sizeof(int) * (size_t)(2 * m + 1)); int n = 0;
    int *eu = malloc(sizeof(int) * (size_t)(m + 1)), *ev = malloc(sizeof(int) * (size_t)(m + 1)); double *w = malloc(sizeof(double) * (size_t)(m + 1));
    if (!ids || !eu || !ev || !w) { fprintf(stderr, "out of memory\n"); free(ids); free(eu); free(ev); free(w); hash_free(&idx); da_free(&raw); da_free(&best); return 1; }
    int idx_fail = 0;
    for (int e = 0; e < m && !idx_fail; e++) {
        int64_t v;
        if (!hash_get(&idx, B[e].id_a, &v)) { v = n; ids[n++] = B[e].id_a; if (hash_put(&idx, B[e].id_a, v)) { idx_fail = 1; break; } }
        eu[e] = (int)v;
        if (!hash_get(&idx, B[e].id_b, &v)) { v = n; ids[n++] = B[e].id_b; if (hash_put(&idx, B[e].id_b, v)) { idx_fail = 1; break; } }
        ev[e] = (int)v; w[e] = B[e].miss;
    }
    if (idx_fail) { fprintf(stderr, "out of memory\n"); free(ids); free(eu); free(ev); free(w); hash_free(&idx); da_free(&raw); da_free(&best); return 1; }
    Graph g; if (graph_build(&g, n, m, eu, ev, w)) { fprintf(stderr, "graph build failed\n"); free(ids); free(eu); free(ev); free(w); hash_free(&idx); da_free(&raw); da_free(&best); return 1; }
    unsigned char *mov = NULL; int nmov = n;
    if (movable_csv) {
        Catalog mc; LoadReport rep; cat_init(&mc); report_init(&rep);
        if (cat_load_csv(&mc, movable_csv, &rep) != 0) {
            fprintf(stderr, "cannot load %s\n", movable_csv);
            report_free(&rep); cat_free(&mc); graph_free(&g); free(ids); free(eu); free(ev); free(w); hash_free(&idx); da_free(&raw); da_free(&best); return 1;
        }
        mov = calloc((size_t)(n > 0 ? n : 1), 1); nmov = 0;
        if (!mov) { fprintf(stderr, "out of memory\n"); report_free(&rep); cat_free(&mc); graph_free(&g); free(ids); free(eu); free(ev); free(w); hash_free(&idx); da_free(&raw); da_free(&best); return 1; }
        for (int v = 0; v < n; v++) if (cat_find(&mc, ids[v]) >= 0) { mov[v] = 1; nmov++; }
        report_free(&rep); cat_free(&mc);
    }
    /* From here on every allocation is checked; BAIL releases everything allocated so far and reports failure (tests/test_oom.c). */
    int *c1 = malloc(sizeof(int) * (size_t)(n + 1)), *c2 = malloc(sizeof(int) * (size_t)(n + 1)), *c3 = malloc(sizeof(int) * (size_t)(n + 1));
    int *csz = NULL, *cm = NULL, *top = NULL, *dist = NULL; unsigned char *gcov = NULL, *xcov = NULL;
#define BAIL() do { fprintf(stderr, "out of memory\n"); free(c1); free(c2); free(c3); free(csz); free(cm); free(top); free(dist); free(gcov); free(xcov); free(mov); \
                    graph_free(&g); hash_free(&idx); free(ids); free(eu); free(ev); free(w); da_free(&raw); da_free(&best); return 1; } while (0)
    if (!c1 || !c2 || !c3) BAIL();
    double t0 = now_s(); int k1 = graph_components_bfs(&g, c1); double tb = now_s() - t0;
    t0 = now_s(); int k2 = graph_components_dfs(&g, c2); double td = now_s() - t0;
    t0 = now_s(); int k3 = graph_components_dsu(&g, c3); double tu = now_s() - t0;
    int agree = (k1 == k2 && k2 == k3 && memcmp(c1, c2, sizeof(int) * (size_t)n) == 0 && memcmp(c1, c3, sizeof(int) * (size_t)n) == 0);
    if (k1 < 0 || k2 < 0 || k3 < 0) BAIL();
    int k = k1;
    csz = calloc((size_t)(k + 1), sizeof(int)); cm = calloc((size_t)(k + 1), sizeof(int));
    if (!csz || !cm) BAIL();
    for (int v = 0; v < n; v++) csz[c1[v]]++;
    for (int e = 0; e < m; e++) cm[c1[eu[e]]]++;
    int largest = 0, pairs_only = 0; for (int c = 0; c < k; c++) { if (csz[c] > csz[largest]) largest = c; if (csz[c] == 2) pairs_only++; }
    /* Top-K closest approaches among pairs */
    top = malloc(sizeof(int) * (size_t)(K > 0 ? K : 1));
    if (!top) BAIL();
    int nt = topk_smallest(w, m, K, top);
    if (nt < 0) BAIL();
    char path[1024];
    snprintf(path, sizeof path, "%s_topk.csv", prefix);
    FILE *f = fopen(path, "w");
    if (f) {
        fprintf(f, "rank,norad_id_a,norad_id_b,tca_s,miss_km,rel_speed_kms,degree_a,degree_b,cluster_size\n");
        for (int r = 0; r < nt; r++) { int e = top[r]; fprintf(f, "%d,%d,%d,%.3f,%.4f,%.3f,%d,%d,%d\n", r + 1, B[e].id_a, B[e].id_b, B[e].tca, B[e].miss, B[e].vrel, graph_degree(&g, eu[e]), graph_degree(&g, ev[e]), csz[c1[eu[e]]]); }
        fclose(f);
    }
    /* Top-K again, crossing encounters only (relative speed > 1 km/s): objects drifting together at ~0 km/s dominate the plain list */
    {
        double *kc = malloc(sizeof(double) * (size_t)(m + 1)); int ncr = 0;
        int *tc = malloc(sizeof(int) * (size_t)(K > 0 ? K : 1));
        if (!kc || !tc) { free(kc); free(tc); BAIL(); }
        for (int e = 0; e < m; e++) { kc[e] = B[e].vrel > 1.0 ? B[e].miss : 1e30; if (B[e].vrel > 1.0) ncr++; }
        int ntc = topk_smallest(kc, m, K, tc);
        if (ntc < 0) { free(kc); free(tc); BAIL(); }
        snprintf(path, sizeof path, "%s_topk_crossing.csv", prefix);
        FILE *fx = fopen(path, "w");
        if (fx) { fprintf(fx, "rank,norad_id_a,norad_id_b,tca_s,miss_km,rel_speed_kms,degree_a,degree_b,cluster_size\n");
            for (int r = 0; r < ntc; r++) { int e = tc[r]; if (kc[e] >= 1e29) break; fprintf(fx, "%d,%d,%d,%.3f,%.4f,%.3f,%d,%d,%d\n", r + 1, B[e].id_a, B[e].id_b, B[e].tca, B[e].miss, B[e].vrel, graph_degree(&g, eu[e]), graph_degree(&g, ev[e]), csz[c1[eu[e]]]); }
            fclose(fx); }
        printf("pairs with relative speed > 1 km/s (crossing encounters): %d of %d\n", ncr, m);
        free(kc); free(tc);
    }
    /* time-ordered index of the conjunctions: AVL tree keyed by TCA in milliseconds (a key already taken is bumped by 1 ms) */
    {
        Avl tree; avl_init(&tree); Timeline tl; tl.nb = 1;
        /* only TCAs in [0, TIMELINE_MAX_S] are indexed: the key is a 32-bit millisecond count (2.0e6 s = 23 days fits; the
           casts below would otherwise be undefined for negative, NaN, infinite or huge times in a damaged file) */
        long skipped = 0;
        for (int e = 0; e < m; e++) { if (!(B[e].tca >= 0.0 && B[e].tca <= TIMELINE_MAX_S)) continue; int b = (int)(B[e].tca / 600.0) + 1; if (b > tl.nb) tl.nb = b; }
        tl.bins = calloc((size_t)tl.nb, sizeof(int));
        if (!tl.bins) { avl_free(&tree); BAIL(); }
        for (int e = 0; e < m; e++) {
            if (!(B[e].tca >= 0.0 && B[e].tca <= TIMELINE_MAX_S)) { skipped++; continue; }
            int key = (int)(B[e].tca * 1000.0); void *dummy;
            while (avl_find(&tree, key, &dummy)) key++;
            if (avl_insert(&tree, key, (void *)(size_t)(e + 1))) { free(tl.bins); avl_free(&tree); BAIL(); }
        }
        if (skipped) printf("WARNING: %ld conjunction(s) with a TCA outside [0, %.0f] s are not in the time index\n", skipped, (double)TIMELINE_MAX_S);
        avl_foreach(&tree, timeline_visit, &tl);
        snprintf(path, sizeof path, "%s_timeline.csv", prefix);
        FILE *ft = fopen(path, "w");
        if (ft) { fprintf(ft, "minutes_from_start,conjunction_pairs_with_tca_in_bin\n"); for (int b = 0; b < tl.nb; b++) fprintf(ft, "%d,%d\n", b * 10, tl.bins[b]); fclose(ft); }
        int mx = 0; for (int b = 0; b < tl.nb; b++) if (tl.bins[b] > mx) mx = tl.bins[b];
        printf("time index (AVL, height %d for %d conjunctions): busiest 10-minute bin has %d conjunction pairs\n", avl_height(&tree), m, mx);
        free(tl.bins); avl_free(&tree);
    }
    /* mitigation per component */
    snprintf(path, sizeof path, "%s_components.csv", prefix);
    FILE *fc = fopen(path, "w");
    if (fc) fprintf(fc, "component,nodes,edges,movable_nodes,lower_bound,greedy,exact,exact_complete,exact_nodes,greedy_us,exact_us,unresolved_edges\n");
    gcov = calloc((size_t)(n + 1), 1); xcov = calloc((size_t)(n + 1), 1);
    if (!gcov || !xcov) { if (fc) fclose(fc); BAIL(); }
    long tot_greedy = 0, tot_exact = 0, tot_lb = 0, tot_unres = 0, n_exact_done = 0, n_nontrivial = 0, n_gap = 0; double t_greedy = 0, t_exact = 0;
    int max_nontrivial_done = 0;
    for (int c = 0; c < k; c++) {
        Graph sg; int *orig = malloc(sizeof(int) * (size_t)csz[c]);
        if (!orig) { if (fc) fclose(fc); BAIL(); }
        if (graph_subgraph(&g, c1, c, &sg, orig)) { free(orig); if (fc) fclose(fc); BAIL(); }
        unsigned char *smov = NULL, *sc1 = calloc((size_t)sg.n + 1, 1), *sc2 = calloc((size_t)sg.n + 1, 1);
        if (mov) smov = malloc((size_t)sg.n + 1);
        if (!sc1 || !sc2 || (mov && !smov)) { graph_free(&sg); free(orig); free(smov); free(sc1); free(sc2); if (fc) fclose(fc); BAIL(); }
        if (mov) for (int i = 0; i < sg.n; i++) smov[i] = mov[orig[i]];
        int unres = 0, nm = 0; for (int i = 0; i < sg.n; i++) nm += smov ? smov[i] : 1;
        double a = now_s(); int gs = vc_greedy(&sg, smov, sc1, &unres); double tg = now_s() - a;
        int lb = vc_lower_bound(&sg, smov);
        if (gs < 0 || lb < 0) { graph_free(&sg); free(orig); free(smov); free(sc1); free(sc2); if (fc) fclose(fc); BAIL(); }
        long nodes = 0; int complete = 1, xs; double tx;
        if (sg.n <= 2) { xs = gs; tx = 0; memcpy(sc2, sc1, (size_t)sg.n); }          /* a single pair: trivially optimal */
        else { a = now_s(); xs = vc_exact(&sg, smov, exact_node_limit, sc2, &nodes, &complete, &unres); tx = now_s() - a; if (xs < 0) { graph_free(&sg); free(orig); free(smov); free(sc1); free(sc2); if (fc) fclose(fc); BAIL(); } n_nontrivial++; if (complete) { n_exact_done++; if (sg.n > max_nontrivial_done) max_nontrivial_done = sg.n; } if (complete && xs < gs) n_gap++; }
        tot_greedy += gs; tot_exact += xs; tot_lb += lb; tot_unres += unres; t_greedy += tg; t_exact += tx;
        for (int i = 0; i < sg.n; i++) { if (sc1[i]) gcov[orig[i]] = 1; if (sc2[i]) xcov[orig[i]] = 1; }
        if (fc && sg.n > 2) fprintf(fc, "%d,%d,%d,%d,%d,%d,%d,%d,%ld,%.1f,%.1f,%d\n", c, sg.n, sg.m, nm, lb, gs, xs, complete, nodes, tg * 1e6, tx * 1e6, unres);
        graph_free(&sg); free(orig); free(smov); free(sc1); free(sc2);
    }
    if (fc) fclose(fc);
    snprintf(path, sizeof path, "%s_cover.csv", prefix);
    f = fopen(path, "w");
    if (f) { fprintf(f, "norad_id,in_greedy_cover,in_best_cover\n"); for (int v = 0; v < n; v++) if (gcov[v] || xcov[v]) fprintf(f, "%d,%d,%d\n", ids[v], gcov[v], xcov[v]); fclose(f); }
    /* cascade reach: hop distances from the highest-degree node of the largest component */
    int hub = -1; for (int v = 0; v < n; v++) if (c1[v] == largest && (hub < 0 || graph_degree(&g, v) > graph_degree(&g, hub))) hub = v;
    dist = malloc(sizeof(int) * (size_t)(n + 1));
    if (!dist) BAIL();
    int maxhop = 0, reached = 0;
    if (hub >= 0) { reached = graph_hops(&g, hub, dist); if (reached < 0) BAIL(); for (int v = 0; v < n; v++) if (dist[v] > maxhop) maxhop = dist[v]; }
    int maxdeg = 0; for (int v = 0; v < n; v++) if (graph_degree(&g, v) > maxdeg) maxdeg = graph_degree(&g, v);

    printf("events %zu, conjunction pairs (edges) %d, objects involved (nodes) %d, movable among them %d\n", raw.n, m, n, nmov);
    printf("components %d (%d are single pairs), largest %d objects / %d pairs, max degree %d; BFS/DFS/DSU agree: %s (%.1f / %.1f / %.1f ms)\n", k, pairs_only, csz[largest], cm[largest], maxdeg, agree ? "yes" : "NO", tb * 1e3, td * 1e3, tu * 1e3);
    printf("cascade reach from the busiest object of the largest cluster: %d objects within %d hops\n", reached, maxhop);
    printf("mitigation (min vertex cover): greedy %ld objects (%.1f ms), exact/best %ld (%.1f ms), lower bound %ld, unresolvable pairs (no movable object) %ld\n", tot_greedy, t_greedy * 1e3, tot_exact, t_exact * 1e3, tot_lb, tot_unres);
    printf("exact search finished on %ld of %ld clusters with 3+ objects (largest finished: %d objects); greedy strictly worse than optimal on %ld\n", n_exact_done, n_nontrivial, max_nontrivial_done, n_gap);
    printf("C6SUMMARY,%s,%zu,%d,%d,%d,%d,%d,%d,%d,%d,%d,%ld,%ld,%ld,%ld,%ld,%ld,%ld,%.6f,%.6f\n", events_csv, raw.n, m, n, nmov, k, pairs_only, csz[largest], cm[largest], agree, nt,
           tot_greedy, tot_exact, tot_lb, tot_unres, n_nontrivial, n_exact_done, n_gap, t_greedy, t_exact);
    free(dist); free(gcov); free(xcov); free(top); free(csz); free(cm); free(c1); free(c2); free(c3); free(mov);
    graph_free(&g); hash_free(&idx); free(ids); free(eu); free(ev); free(w); da_free(&raw); da_free(&best);
    return 0;
}
