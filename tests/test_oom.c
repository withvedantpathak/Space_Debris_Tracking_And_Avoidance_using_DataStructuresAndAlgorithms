/* C8: out-of-memory robustness. malloc/calloc/realloc are wrapped (link with -Wl,--wrap=malloc,--wrap=calloc,--wrap=realloc; GNU ld,
 * so Linux/MinGW only). For every scenario the harness first counts the allocations, then re-runs it once per allocation with that
 * allocation (and, in a second sweep, EVERY allocation from there on) failing. Required of the code under test:
 *   - never crashes, never reads or writes freed/uninitialised memory (ASan/UBSan build), never leaks (LeakSanitizer at exit);
 *   - either reports an error (non-zero) or returns EXACTLY the result of the unfailed run (no silently wrong answer).
 * Library calls made by libc itself (fopen, printf buffers) are not affected: only the project's own calls go through the wrappers. */
#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include "catalog.h"
#include "synth.h"
#include "propagate.h"
#include "conjunction.h"
#include "screen.h"
#include "graph.h"
#include "mitigate.h"
#include "rank.h"
#include "analyze.h"
#include "ds/darray.h"
#include "ds/heap.h"
#include "ds/hash.h"
#include "ds/queue.h"
#include "ds/stack.h"
#include "ds/avl.h"
#include "ds/kdtree.h"
#include "ds/grid.h"
#include "ds/dsu.h"

void *__real_malloc(size_t n);
void *__real_calloc(size_t a, size_t b);
void *__real_realloc(void *p, size_t n);
static long alloc_count = 0, fail_at = 0;      /* fail_at: 1-based index of the allocation to fail, 0 = never */
static int sticky = 0, armed = 0;
static int should_fail(void) {
    if (!armed) return 0;
    alloc_count++;
    return fail_at > 0 && (alloc_count == fail_at || (sticky && alloc_count > fail_at));
}
void *__wrap_malloc(size_t n) { return should_fail() ? NULL : __real_malloc(n); }
void *__wrap_calloc(size_t a, size_t b) { return should_fail() ? NULL : __real_calloc(a, b); }
void *__wrap_realloc(void *p, size_t n) { return should_fail() ? NULL : __real_realloc(p, n); }

static long checks = 0, fails = 0;
#define CHECKF(c, ...) do { checks++; if (!(c)) { fails++; printf("FAIL %s:%d: %s | ", __FILE__, __LINE__, #c); printf(__VA_ARGS__); printf("\n"); } } while (0)

/* ---- fixtures built once, outside the injection ---- */
static Catalog gcat; static Prop *gprops; static double *grmin, *grmax; static int gn;
static const char *CSV = "build/_oom_cat.csv", *EVF = "build/_oom_events.csv", *MOVF = "build/_oom_mov.csv";

typedef int (*Scenario)(long *sig);          /* returns 0 on success (sig = result) or non-zero when it reports an error */

/* ---- scenarios ---- */
static int sc_structures(long *sig) {
    long s = 0; int rc = -1;
    DArray da; da_init(&da, sizeof(int)); Hash h; int hh = 0; Heap hp; int hpi = 0; Queue q; queue_init(&q); Stack st; int sti = 0;
    Dsu d; int di = 0; Avl t; avl_init(&t);
    for (int i = 0; i < 300; i++) if (da_push(&da, &i)) goto out;
    if (hash_init(&h, 4)) goto out;
    hh = 1;
    for (int i = 0; i < 300; i++) if (hash_put(&h, (int64_t)i * 7919, i)) goto out;
    if (heap_init(&hp, 2)) goto out;
    hpi = 1;
    for (int i = 0; i < 300; i++) if (heap_push(&hp, (double)((i * 37) % 101), i)) goto out;
    for (int i = 0; i < 300; i++) if (queue_push(&q, i)) goto out;
    if (stack_init(&st, 1)) goto out;
    sti = 1;
    for (int i = 0; i < 300; i++) if (stack_push(&st, i)) goto out;
    if (dsu_init(&d, 300)) goto out;
    di = 1;
    for (int i = 1; i < 300; i += 2) dsu_union(&d, i - 1, i);
    for (int i = 0; i < 200; i++) if (avl_insert(&t, (i * 7919) % 1009, NULL)) goto out;
    { HeapItem it; int v; while (heap_pop(&hp, &it)) s += it.val; while (queue_pop(&q, &v)) s += v; while (stack_pop(&st, &v)) s += v; }
    s += (long)da.n + (long)hash_count(&h) + dsu_size(&d, 0) + avl_height(&t);
    rc = 0;
out:
    da_free(&da); if (hh) hash_free(&h); if (hpi) heap_free(&hp); queue_free(&q); if (sti) stack_free(&st); if (di) dsu_free(&d); avl_free(&t);
    *sig = s; return rc;
}
static int sc_spatial(long *sig) {
    Vec3 pts[200]; for (int i = 0; i < 200; i++) { pts[i].x = (i * 37) % 101; pts[i].y = (i * 53) % 97; pts[i].z = (i * 11) % 89; }
    KDTree kd; long s = 0; int rc = kd_build(&kd, pts, 200);
    if (rc) return -1;
    s += kd_nearest(&kd, pts[3], 3, &(double){0}); kd_free(&kd);
    Grid g; rc = grid_build(&g, pts, 200, 12.0);
    if (rc) return -1;
    grid_free(&g); *sig = s; return 0;
}
static int sc_catalog(long *sig) {
    Catalog c; LoadReport rep; cat_init(&c); report_init(&rep);
    int rc = cat_load_csv(&c, CSV, &rep), r2 = 0; long s = 0;   /* the filters return the number kept, or < 0 on error */
    if (rc == 0) {
        Catalog f; cat_init(&f); r2 = cat_filter_mean_motion(&c, &f, 11.0);
        if (r2 >= 0) { Catalog f2; cat_init(&f2); r2 = cat_filter_altitude(&f, &f2, 100.0, 2000.0); if (r2 >= 0) s += cat_count(&f2); cat_free(&f2); }
        cat_free(&f);
        s += cat_count(&c) * 1000;
        for (int i = 0; i < cat_count(&c); i++) s += cat_objects(&c)[i].norad_id % 1000;
    }
    cat_free(&c); report_free(&rep);
    *sig = s; return rc != 0 || r2 < 0;
}
static int sc_synth(long *sig) {
    Catalog c; cat_init(&c); DArray truth; da_init(&truth, sizeof(PlantedPair));
    int rc = synth_generate(&c, &truth, 40, 5, 2); long s = 0;
    if (rc == 0) { s = cat_count(&c) * 100000L + (long)truth.n; for (int i = 0; i < cat_count(&c); i++) s += cat_objects(&c)[i].norad_id % 997; }
    cat_free(&c); da_free(&truth); *sig = s; return rc;
}
static int sc_screen(long *sig) {
    ScreenParams sp = {1200.0, 20.0, 40.0, NULL}; long s = 0; int bad = 0;
    for (int m = 0; m < 4; m++) {
        DArray ev; da_init(&ev, sizeof(Conjunction)); int rc;
        if (m == 0) rc = conj_bruteforce(gprops, gn, &sp, &ev, NULL);
        else { FastScreen fs; memset(&fs, 0, sizeof fs); fs.method = m == 1 ? SCREEN_KDTREE : m == 2 ? SCREEN_GRID : SCREEN_ALT; fs.rmin_km = grmin; fs.rmax_km = grmax; fs.alt_pad_km = 50.0; rc = conj_fast(gprops, gn, &sp, &fs, &ev, NULL); }
        if (rc) bad = 1; else s = s * 31 + (long)ev.n;
        da_free(&ev);
    }
    *sig = s; return bad;
}
static int sc_graph(long *sig) {
    enum { N = 24 };
    int eu[200], ev[200], m = 0; double w[200]; unsigned r = 12345;
    for (int i = 0; i < N; i++) for (int j = i + 1; j < N; j++) { r = r * 1103515245u + 12345u; if ((r >> 16) % 100 < 14) { eu[m] = i; ev[m] = j; w[m] = (double)((r >> 8) % 50); m++; } }
    Graph g; if (graph_build(&g, N, m, eu, ev, w)) return -1;
    int c1[N], c2[N], c3[N], dist[N], bad = 0; long s = 0;
    int k1 = graph_components_bfs(&g, c1), k2 = graph_components_dfs(&g, c2), k3 = graph_components_dsu(&g, c3);
    if (k1 < 0 || k2 < 0 || k3 < 0) bad = 1; else s += k1 + k2 + k3;
    int reach = graph_hops(&g, 0, dist); if (reach < 0) bad = 1; else s += reach;
    if (!bad) { Graph sg; int orig[N]; if (graph_subgraph(&g, c1, c1[0], &sg, orig)) bad = 1; else { s += sg.m; graph_free(&sg); } }
    unsigned char cov[N], cov2[N]; int un = 0; long nodes = 0; int complete = 0;
    int gs = vc_greedy(&g, NULL, cov, &un); if (gs < 0) bad = 1; else s += gs * 3;
    int xs = vc_exact(&g, NULL, 200000, cov2, &nodes, &complete, &un); if (xs < 0) bad = 1; else s += xs * 7 + complete;
    int lb = vc_lower_bound(&g, NULL); if (lb < 0) bad = 1; else s += lb;
    int out[5], nt = topk_smallest(w, m, 5, out); if (nt < 0) bad = 1; else s += nt;
    graph_free(&g); *sig = s; return bad;
}
static long file_sig(const char *path) {                       /* position-weighted byte sum of a file: any change in content changes it */
    FILE *f = fopen(path, "rb"); if (!f) return -1;
    long s = 0, i = 1; int ch; while ((ch = fgetc(f)) != EOF) { s += ch * (i % 251 + 1); i++; }
    fclose(f); return s;
}
static int sc_analyze(long *sig) {
    int rc = analyze_events(EVF, MOVF, 3, "build/_oom_out", 100000);
    *sig = rc == 0 ? file_sig("build/_oom_out_topk.csv") * 7 + file_sig("build/_oom_out_cover.csv") * 13 + file_sig("build/_oom_out_components.csv") * 17 + file_sig("build/_oom_out_timeline.csv") : 0;
    return rc;
}
/* a catalogue must stay consistent (array and index agree) after a failed cat_add, and keep working afterwards */
static int sc_catalog_add(long *sig) {
    Catalog c; if (cat_init(&c) != 0) { cat_free(&c); return -1; }
    int fails_seen = 0; long s = 0;
    for (int i = 0; i < gn; i++) if (cat_add(&c, &cat_objects(&gcat)[i]) < 0) fails_seen = 1;
    for (int i = 0; i < gn; i++) cat_add(&c, &cat_objects(&gcat)[i]);                 /* second pass: duplicates, or the objects that failed */
    int consistent = 1;
    for (int i = 0; i < cat_count(&c); i++) { int id = cat_objects(&c)[i].norad_id; if (cat_find(&c, id) != i) consistent = 0; }
    for (int i = 0; i < cat_count(&c); i++) for (int j = i + 1; j < cat_count(&c); j++) if (cat_objects(&c)[i].norad_id == cat_objects(&c)[j].norad_id) consistent = 0;   /* no duplicate ids */
    s = consistent ? (long)cat_count(&c) : -999;
    cat_free(&c);
    *sig = s;
    return fails_seen && s != -999 && s != (long)gn ? 1 : 0;   /* an error was reported and the catalogue is merely incomplete: fine */
}

/* ---- harness ---- */
static void sweep(const char *name, Scenario fn, long max_points) {
    long s0 = 0; int quiet = !strcmp(name, "analyze");
    if (quiet) { fflush(stdout); if (!freopen("/dev/null", "w", stdout)) return; }
    alloc_count = 0; fail_at = 0; sticky = 0; armed = 1;
    int rc0 = fn(&s0); armed = 0; long N = alloc_count;
    long stride = N > max_points ? (N + max_points - 1) / max_points : 1, runs = 0, errs = 0, same = 0, wrong = 0;
    for (int mode = 0; mode < 2; mode++)
        for (long k = 1; k <= N; k += stride) {
            long s = 0; alloc_count = 0; fail_at = k; sticky = mode; armed = 1;
            int rc = fn(&s); armed = 0; runs++;
            if (rc != 0) errs++;
            else if (s == s0) same++;
            else { wrong++; if (wrong <= 3) { if (quiet) { fflush(stdout); } fprintf(stderr, "  WRONG RESULT %s: failing allocation %ld (%s) gave %ld, expected %ld\n", name, k, mode ? "sticky" : "single", s, s0); } }
        }
    fail_at = 0; sticky = 0;
    if (quiet) { fflush(stdout); if (!freopen("/dev/tty", "w", stdout)) { stdout = fdopen(dup(2), "w"); } }
    printf("%-12s %6ld allocations, %6ld failure runs: %6ld reported an error, %5ld unaffected, %ld silently wrong\n", name, N, runs, errs, same, wrong);
    CHECKF(rc0 == 0, "%s: unfailed run returned %d", name, rc0);
    CHECKF(wrong == 0, "%s: %ld runs gave a wrong answer without reporting an error", name, wrong);
    CHECKF(N > 0, "%s: no allocations were seen: the wrappers are not linked", name);
}

int main(int argc, char **argv) {
    const char *only = argc > 1 ? argv[1] : "all";
    /* fixtures */
    cat_init(&gcat); if (synth_generate(&gcat, NULL, 40, 11, 2) != 0) { printf("fixture failed\n"); return 2; }
    gn = cat_count(&gcat); double ref; iso_to_jd(SYNTH_EPOCH, &ref);
    gprops = malloc(sizeof(Prop) * (size_t)gn); grmin = malloc(sizeof(double) * (size_t)gn); grmax = malloc(sizeof(double) * (size_t)gn);
    for (int i = 0; i < gn; i++) { prop_init(&gprops[i], &cat_objects(&gcat)[i], PROP_J2, ref); grmin[i] = cat_objects(&gcat)[i].perigee_r_km; grmax[i] = cat_objects(&gcat)[i].apogee_r_km; }
    if (cat_write_csv(&gcat, CSV) != 0) { printf("cannot write fixture\n"); return 2; }
    { FILE *f = fopen(EVF, "w"); if (!f) return 2; fprintf(f, "h\n");
      for (int i = 0; i < 30; i++) fprintf(f, "%d,%d,%d.0,%d.5,7.0\n", 900000 + i, 900000 + (i * 7 + 1) % 40 + (i == (i * 7 + 1) % 40), 100 + 37 * i, 1 + i % 9);
      fclose(f); }
    { Catalog mv; cat_init(&mv); for (int i = 0; i < gn; i += 2) cat_add(&mv, &cat_objects(&gcat)[i]); cat_write_csv(&mv, MOVF); cat_free(&mv); }
    struct { const char *name; Scenario fn; long pts; } all[] = {
        {"structures", sc_structures, 4000}, {"spatial", sc_spatial, 4000}, {"catalog", sc_catalog, 1500}, {"synth", sc_synth, 1500},
        {"screen", sc_screen, 1500}, {"graph", sc_graph, 4000}, {"analyze", sc_analyze, 1500}, {"catalog_add", sc_catalog_add, 1500}};
    for (unsigned i = 0; i < sizeof all / sizeof all[0]; i++) if (!strcmp(only, "all") || !strcmp(only, all[i].name)) sweep(all[i].name, all[i].fn, all[i].pts);
    cat_free(&gcat); free(gprops); free(grmin); free(grmax);
    printf("%ld checks, %ld failed\n", checks, fails);
    return fails != 0;
}
