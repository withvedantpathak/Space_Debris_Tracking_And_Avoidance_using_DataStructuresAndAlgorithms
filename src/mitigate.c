#include <stdlib.h>
#include <string.h>
#include "mitigate.h"
#include "ds/heap.h"
#include "ds/stack.h"
#include "ds/dsu.h"

int vc_use_reductions = 1;

static int is_mov(const unsigned char *mov, int v) { return mov ? mov[v] != 0 : 1; }
static int clearable(const Graph *g, const unsigned char *mov, int e) { return is_mov(mov, g->eu[e]) || is_mov(mov, g->ev[e]); }
static int oth(const Graph *g, int e, int v) { return g->eu[e] == v ? g->ev[e] : g->eu[e]; }

/* ---- greedy over an edge list -------------------------------------------------------------------------------- */
/* Adds picks to in_cover. el[0..ne) = edge ids still to be covered (all clearable). scratch: covered[m] zeroed by us, deg[n]. */
static int greedy_core(const Graph *g, const unsigned char *mov, const int *el, int ne, unsigned char *in_cover, unsigned char *covered, int *deg) {
    int size = 0;
    Heap h; if (heap_init(&h, 64)) return -1;
    for (int i = 0; i < ne; i++) { int e = el[i]; covered[e] = 0; deg[g->eu[e]] = 0; deg[g->ev[e]] = 0; }
    for (int i = 0; i < ne; i++) { int e = el[i]; deg[g->eu[e]]++; deg[g->ev[e]]++; }
#define GPICK(v) do { in_cover[v] = 1; size++; \
        for (int i_ = g->adj_start[v]; i_ < g->adj_start[(v) + 1]; i_++) { int e_ = g->adj_edge[i_]; if (covered[e_] != 0 || !in_list[e_]) continue; \
            covered[e_] = 1; int w_ = oth(g, e_, (v)); deg[w_]--; deg[v]--; \
            if (is_mov(mov, w_) && !in_cover[w_] && deg[w_] > 0 && heap_push(&h, -(double)deg[w_], w_)) { heap_free(&h); free(in_list); return -1; } } } while (0)
    unsigned char *in_list = calloc((size_t)(g->m > 0 ? g->m : 1), 1);          /* membership of el, so picks do not touch other edges */
    if (!in_list) { heap_free(&h); return -1; }
    for (int i = 0; i < ne; i++) { in_list[el[i]] = 1; covered[el[i]] = 0; }
    for (int i = 0; i < ne; i++) {                                              /* forced picks */
        int e = el[i]; if (covered[e]) continue;
        int u = g->eu[e], v = g->ev[e];
        if (is_mov(mov, u) && is_mov(mov, v)) continue;
        int p = is_mov(mov, u) ? u : v;
        if (!in_cover[p]) GPICK(p);
    }
    for (int i = 0; i < ne; i++) { int e = el[i]; int ends[2] = {g->eu[e], g->ev[e]};
        for (int s = 0; s < 2; s++) { int v = ends[s]; if (is_mov(mov, v) && !in_cover[v] && deg[v] > 0 && heap_push(&h, -(double)deg[v], v)) { heap_free(&h); free(in_list); return -1; } } }
    HeapItem it;
    while (heap_pop(&h, &it)) {
        int v = it.val;
        if (in_cover[v] || deg[v] <= 0 || (int)(-it.key) != deg[v]) continue;    /* stale entry (lazy decrease-key) */
        GPICK(v);
    }
#undef GPICK
    /* redundancy pass: drop a chosen object when every conjunction it was clearing is also cleared by the other object */
    for (int i = 0; i < ne; i++) {
        int ends[2] = {g->eu[el[i]], g->ev[el[i]]};
        for (int s2 = 0; s2 < 2; s2++) {
            int v = ends[s2];
            if (!in_cover[v]) continue;
            int needed = 0;
            for (int j = g->adj_start[v]; j < g->adj_start[v + 1] && !needed; j++) {
                int e = g->adj_edge[j]; if (!in_list[e]) continue;
                if (!in_cover[oth(g, e, v)]) needed = 1;
            }
            if (!needed) { in_cover[v] = 0; size--; }
        }
    }
    heap_free(&h); free(in_list);
    return size;
}

int vc_greedy(const Graph *g, const unsigned char *mov, unsigned char *in_cover, int *unresolved) {
    int m = g->m, ne = 0, unres = 0;
    memset(in_cover, 0, (size_t)(g->n > 0 ? g->n : 1));
    int *el = malloc(sizeof(int) * (size_t)(m > 0 ? m : 1)), *deg = calloc((size_t)(g->n > 0 ? g->n : 1), sizeof(int));
    unsigned char *cov = calloc((size_t)(m > 0 ? m : 1), 1);
    if (!el || !deg || !cov) { free(el); free(deg); free(cov); return -1; }
    for (int e = 0; e < m; e++) { if (clearable(g, mov, e)) el[ne++] = e; else unres++; }
    int s = greedy_core(g, mov, el, ne, in_cover, cov, deg);
    free(el); free(deg); free(cov);
    if (unresolved) *unresolved = unres;
    return s;
}

/* ---- lower bound --------------------------------------------------------------------------------------------- */
static int matching_bound(const Graph *g, const int *el, int ne, const int *cnt, unsigned char *mark) {
    for (int i = 0; i < ne; i++) { mark[g->eu[el[i]]] = 0; mark[g->ev[el[i]]] = 0; }
    int lb = 0;
    for (int i = 0; i < ne; i++) {
        int e = el[i];
        if (cnt[e] > 0 || mark[g->eu[e]] || mark[g->ev[e]]) continue;
        mark[g->eu[e]] = mark[g->ev[e]] = 1; lb++;
    }
    return lb;
}
int vc_lower_bound(const Graph *g, const unsigned char *mov) {
    int *el = malloc(sizeof(int) * (size_t)(g->m > 0 ? g->m : 1)), *cnt = calloc((size_t)(g->m > 0 ? g->m : 1), sizeof(int));
    unsigned char *mark = calloc((size_t)(g->n > 0 ? g->n : 1), 1);
    if (!el || !cnt || !mark) { free(el); free(cnt); free(mark); return -1; }
    int ne = 0; for (int e = 0; e < g->m; e++) if (clearable(g, mov, e)) el[ne++] = e;
    int lb = matching_bound(g, el, ne, cnt, mark);
    free(el); free(cnt); free(mark);
    return lb;
}

/* ---- exact --------------------------------------------------------------------------------------------------- */
typedef struct { const Graph *g; const unsigned char *mov; int *cnt; unsigned char *chosen, *mark; int *cand, *ncand, *opt; } Bb;

/* Branch and bound on the edges el[0..ne) (one connected group of the reduced problem). Returns the optimum size for the
 * group and fills best[] (n bytes, vertices of the group only), or the best found if the budget ran out (*complete = 0). */
static int bb_group(Bb *b, const int *el, int ne, unsigned char *best_mark, int ub, long *budget, int *complete, long *nodes) {
    const Graph *g = b->g; const unsigned char *mov = b->mov; int *cnt = b->cnt; unsigned char *chosen = b->chosen;
    int best = ub, cur = 0, d = 0, enter = 1, ok = 1;
    Stack undo; if (stack_init(&undo, 64)) return -1;
#define CHOOSE(v) do { chosen[v] = 1; cur++; for (int i_ = g->adj_start[v]; i_ < g->adj_start[(v) + 1]; i_++) cnt[g->adj_edge[i_]]++; } while (0)
#define UNCHOOSE(v) do { chosen[v] = 0; cur--; for (int i_ = g->adj_start[v]; i_ < g->adj_start[(v) + 1]; i_++) cnt[g->adj_edge[i_]]--; } while (0)
    for (;;) {
        if (enter) {
            if (*budget <= 0) { ok = 0; break; }
            (*budget)--; (*nodes)++;
            int e = -1;
            for (int i = 0; i < ne; i++) if (cnt[el[i]] == 0) { e = el[i]; break; }
            if (e < 0) {
                if (cur < best) { best = cur; for (int i = 0; i < ne; i++) { best_mark[g->eu[el[i]]] = chosen[g->eu[el[i]]]; best_mark[g->ev[el[i]]] = chosen[g->ev[el[i]]]; } }
                enter = 0;
            } else if (cur + matching_bound(g, el, ne, cnt, b->mark) >= best) {
                enter = 0;
            } else {
                int k = 0;
                if (is_mov(mov, g->eu[e])) b->cand[2 * d + k++] = g->eu[e];
                if (is_mov(mov, g->ev[e])) b->cand[2 * d + k++] = g->ev[e];
                b->ncand[d] = k; b->opt[d] = 0;
                int v = b->cand[2 * d + b->opt[d]++];
                CHOOSE(v);
                if (stack_push(&undo, v)) { stack_free(&undo); return -1; }
                d++;
                continue;
            }
        }
        int again = 0, v;
        while (d > 0) {                                          /* backtrack: pop the undo stack */
            d--; stack_pop(&undo, &v); UNCHOOSE(v);
            if (b->opt[d] < b->ncand[d]) {
                v = b->cand[2 * d + b->opt[d]++];
                CHOOSE(v);
                if (stack_push(&undo, v)) { stack_free(&undo); return -1; }
                d++; again = 1; break;
            }
        }
        if (!again) break;
        enter = 1;
    }
    /* leave the shared state clean even after a cut-off */
    { int v; while (stack_pop(&undo, &v)) UNCHOOSE(v); }
#undef CHOOSE
#undef UNCHOOSE
    stack_free(&undo);
    if (!ok) *complete = 0;
    return best;
}

int vc_exact(const Graph *g, const unsigned char *mov, long node_limit, unsigned char *in_cover, long *nodes_out, int *complete, int *unresolved) {
    int n = g->n, m = g->m, unres = 0;
    size_t N = (size_t)(n > 0 ? n : 1), M = (size_t)(m > 0 ? m : 1);
    memset(in_cover, 0, N);
    int *cnt = calloc(M, sizeof(int)), *deg = calloc(N, sizeof(int)), *el = malloc(sizeof(int) * M), *grp = malloc(sizeof(int) * M), *gcode = malloc(sizeof(int) * M);
    unsigned char *mark = calloc(N, 1), *chosen = calloc(N, 1), *bestm = calloc(N, 1), *scov = calloc(M, 1), *gcov = calloc(N, 1);
    Bb b = {g, mov, cnt, chosen, mark, malloc(sizeof(int) * 2 * (N + M + 2)), malloc(sizeof(int) * (N + M + 2)), malloc(sizeof(int) * (N + M + 2))};
    Dsu ds; int dsu_ok = dsu_init(&ds, n);
    Stack work; int work_ok = stack_init(&work, 64);
    if (!cnt || !deg || !el || !grp || !gcode || !mark || !chosen || !bestm || !scov || !gcov || !b.cand || !b.ncand || !b.opt || dsu_ok || work_ok) {
        free(cnt); free(deg); free(el); free(grp); free(gcode); free(mark); free(chosen); free(bestm); free(scov); free(gcov); free(b.cand); free(b.ncand); free(b.opt);
        if (!dsu_ok) dsu_free(&ds);
        if (!work_ok) stack_free(&work);
        return -1;
    }
    int ne = 0, cur = 0;
    for (int e = 0; e < m; e++) { if (clearable(g, mov, e)) { el[ne++] = e; deg[g->eu[e]]++; deg[g->ev[e]]++; } else unres++; }
    /* 1. reductions that never lose optimality: (a) an edge with one movable endpoint forces that endpoint; (b) a leaf
          (one uncovered edge) is dominated by its neighbour, so pick the neighbour. Applied until nothing changes. */
#define XPICK(x) do { chosen[x] = 1; cur++; for (int i_ = g->adj_start[x]; i_ < g->adj_start[(x) + 1]; i_++) { int e_ = g->adj_edge[i_]; \
            if (cnt[e_]++ == 0 && clearable(g, mov, e_)) { int w_ = oth(g, e_, (x)); deg[w_]--; deg[x]--; if (stack_push(&work, w_)) goto fail; } } } while (0)
    if (vc_use_reductions) {
    for (int i = 0; i < ne; i++) { int e = el[i], u = g->eu[e], v = g->ev[e];
        if (is_mov(mov, u) && is_mov(mov, v)) continue;
        int p = is_mov(mov, u) ? u : v; if (!chosen[p]) XPICK(p); }
    for (int i = 0; i < ne; i++) { if (stack_push(&work, g->eu[el[i]]) || stack_push(&work, g->ev[el[i]])) goto fail; }
    { int v;
      while (stack_pop(&work, &v)) {
        if (chosen[v] || deg[v] != 1) continue;
        int e = -1; for (int i = g->adj_start[v]; i < g->adj_start[v + 1]; i++) { int x = g->adj_edge[i]; if (cnt[x] == 0 && clearable(g, mov, x)) { e = x; break; } }
        if (e < 0) continue;
        int u = oth(g, e, v), p = is_mov(mov, u) ? u : v;      /* neighbour dominates the leaf (it covers e and maybe more) */
        if (!chosen[p]) XPICK(p);
      } }
    }
#undef XPICK
    memcpy(in_cover, chosen, N);
    /* 2. what is left (edges still uncovered) falls apart into independent groups: solve each one separately, smallest first */
    int nrem = 0;
    for (int i = 0; i < ne; i++) if (cnt[el[i]] == 0) { el[nrem++] = el[i]; dsu_union(&ds, g->eu[el[i]], g->ev[el[i]]); }
    int ng = 0;
    for (int i = 0; i < nrem; i++) gcode[i] = dsu_find(&ds, g->eu[el[i]]);
    /* group edges by DSU root: counting sort via root index */
    int *rootcount = calloc(N, sizeof(int)), *order = malloc(sizeof(int) * (size_t)(nrem + 1)), *roots = malloc(sizeof(int) * (size_t)(nrem + 1));
    if (!rootcount || !order || !roots) { free(rootcount); free(order); free(roots); goto fail; }
    for (int i = 0; i < nrem; i++) { if (rootcount[gcode[i]]++ == 0) roots[ng++] = gcode[i]; }
    int *start = malloc(sizeof(int) * (size_t)(ng + 2)), *gsize = malloc(sizeof(int) * (size_t)(ng + 1)), *gperm = malloc(sizeof(int) * (size_t)(ng + 1));
    if (!start || !gsize || !gperm) { free(rootcount); free(order); free(roots); free(start); free(gsize); free(gperm); goto fail; }
    for (int gi = 0; gi < ng; gi++) { gsize[gi] = rootcount[roots[gi]]; gperm[gi] = gi; }
    for (int i = 1; i < ng; i++) { int x = gperm[i], j = i - 1; while (j >= 0 && gsize[gperm[j]] > gsize[x]) { gperm[j + 1] = gperm[j]; j--; } gperm[j + 1] = x; }   /* insertion sort by size */
    { int off = 0; for (int gi = 0; gi < ng; gi++) { start[gi] = off; off += gsize[gi]; rootcount[roots[gi]] = 0; } start[ng] = off;
      int *fill = calloc((size_t)(ng + 1), sizeof(int)); int *pos_of = malloc(sizeof(int) * N);
      if (!fill || !pos_of) { free(fill); free(pos_of); free(rootcount); free(order); free(roots); free(start); free(gsize); free(gperm); goto fail; }
      for (int gi = 0; gi < ng; gi++) pos_of[roots[gi]] = gi;
      for (int i = 0; i < nrem; i++) { int gi = pos_of[gcode[i]]; order[start[gi] + fill[gi]++] = el[i]; }
      free(fill); free(pos_of); }
    long budget = node_limit, nodes = 0; int all_complete = 1;
    for (int q = 0; q < ng; q++) {
        int gi = gperm[q]; const int *ge = order + start[gi]; int gn = gsize[gi];
        int gc = 1;
        for (int i = 0; i < gn; i++) { gcov[g->eu[ge[i]]] = 0; gcov[g->ev[ge[i]]] = 0; }
        int ub = greedy_core(g, mov, ge, gn, gcov, scov, deg);        /* upper bound, and the answer if the budget is gone */
        if (ub < 0) { free(rootcount); free(order); free(roots); free(start); free(gsize); free(gperm); goto fail; }
        for (int i = 0; i < gn; i++) { bestm[g->eu[ge[i]]] = gcov[g->eu[ge[i]]]; bestm[g->ev[ge[i]]] = gcov[g->ev[ge[i]]]; }
        int bsz = ub;
        if (ub > 1) { bsz = bb_group(&b, ge, gn, bestm, ub, &budget, &gc, &nodes); if (bsz < 0) { free(rootcount); free(order); free(roots); free(start); free(gsize); free(gperm); goto fail; } }
        if (!gc) all_complete = 0;
        for (int i = 0; i < gn; i++) { if (bestm[g->eu[ge[i]]]) in_cover[g->eu[ge[i]]] = 1; if (bestm[g->ev[ge[i]]]) in_cover[g->ev[ge[i]]] = 1; }
        cur += bsz;
    }
    if (nodes_out) *nodes_out = nodes;
    if (complete) *complete = all_complete;
    if (unresolved) *unresolved = unres;
    free(rootcount); free(order); free(roots); free(start); free(gsize); free(gperm);
    free(cnt); free(deg); free(el); free(grp); free(gcode); free(mark); free(chosen); free(bestm); free(scov); free(gcov); free(b.cand); free(b.ncand); free(b.opt);
    dsu_free(&ds); stack_free(&work);
    return cur;
fail:
    free(cnt); free(deg); free(el); free(grp); free(gcode); free(mark); free(chosen); free(bestm); free(scov); free(gcov); free(b.cand); free(b.ncand); free(b.opt);
    dsu_free(&ds); stack_free(&work);
    return -1;
}
