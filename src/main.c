/*
 * Space-Debris Conjunction Screening System (working title) - DSA project.
 *
 *   dlscms info <catalog.csv> [more.csv ...]     load CelesTrak OMM CSV file(s), merge by NORAD id, summarise
 *   dlscms gen <n> <seed> <out.csv> [pairs]      write a synthetic catalogue + <out>_truth.csv (planted pairs)
 *   dlscms prop <catalog.csv> <dt_s> [flags] [common]  state of every object dt_s seconds after its own epoch (CSV); flags 0 two-body, 1 +J2, 3 +J2+drag
 *
 *   dlscms screen <out.csv> <window_s> <step_s> <threshold_km> <flags (1=J2, 2=drag, 4=LEO scope, 8=SGP4)> <fresh_days> <min_mean_motion> <catalog.csv> [more.csv ...]
 *                                                exact brute-force conjunction screening (baseline); reference epoch = median epoch;
 *                                                fresh_days > 0 keeps only objects whose epoch is within that many days of it;
 *                                                min_mean_motion 6.4 keeps near-earth objects only (period < 225 min), 0 keeps all
 *
 * Fast screening commands arrive in C5. Tests: scripts/test.sh
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "catalog.h"
#include "synth.h"
#include "propagate.h"
#include "conjunction.h"
#include "screen.h"
#include "analyze.h"
#include "mitigate.h"
#include <time.h>

static const char *orbit_class(const Object *o) {
    double alt_p = o->perigee_r_km - EARTH_RADIUS_KM, alt_a = o->apogee_r_km - EARTH_RADIUS_KM;
    if (o->ecc >= 0.25 && alt_a >= 2000) return "HEO";
    if (alt_a < 2000) return "LEO";
    if (o->sma_km > 41664 && o->sma_km < 42664 && o->ecc < 0.1) return "GEO";
    (void)alt_p;
    return "MEO/other";
}

static int cmd_info(int argc, char **argv) {
    Catalog cat;
    if (cat_init(&cat)) return 1;
    int status = 0;                                  /* 1 if any file could not be opened or has an unusable header */
    for (int i = 0; i < argc; i++) {
        LoadReport rep;
        report_init(&rep);
        int rc = cat_load_csv(&cat, argv[i], &rep);
        if (rc != 0) status = 1;
        if (rc == -1) { printf("%s: cannot open file\n", argv[i]); report_free(&rep); continue; }
        printf("%s: %ld lines | added %ld, replaced %ld, duplicate ignored %ld, suborbital skipped %ld, bad rows %ld%s\n",
               argv[i], rep.lines, rep.rows_added, rep.rows_replaced, rep.rows_dup_ignored, rep.rows_suborbital,
               rep.rows_bad, rc == -2 ? " [HEADER ERROR]" : "");
        for (size_t k = 0; k < rep.errors.n && k < 10; k++) {
            const ParseError *e = (const ParseError *)da_at(&rep.errors, k);
            printf("    line %d: %s\n", e->line, e->msg);
        }
        if (rep.errors_total > 10) printf("    ... %ld errors in total\n", rep.errors_total);
        report_free(&rep);
    }
    int n = cat_count(&cat);
    const Object *o = cat_objects(&cat);
    int cls[4] = {0, 0, 0, 0};
    const char *names[4] = {"LEO", "MEO/other", "GEO", "HEO"};
    double pmin = 1e30, pmax = 0, ejd_min = 1e30, ejd_max = 0;
    for (int i = 0; i < n; i++) {
        const char *c = orbit_class(&o[i]);
        for (int k = 0; k < 4; k++) if (strcmp(c, names[k]) == 0) cls[k]++;
        double p = o[i].perigee_r_km - EARTH_RADIUS_KM;
        if (p < pmin) pmin = p;
        if (o[i].apogee_r_km - EARTH_RADIUS_KM > pmax) pmax = o[i].apogee_r_km - EARTH_RADIUS_KM;
        if (o[i].epoch_jd < ejd_min) ejd_min = o[i].epoch_jd;
        if (o[i].epoch_jd > ejd_max) ejd_max = o[i].epoch_jd;
    }
    printf("merged catalogue: %d objects", n);
    for (int k = 0; k < 4; k++) printf(" | %s %d", names[k], cls[k]);
    printf("\n");
    if (n) {
        char a[40], b[40];
        jd_to_iso(ejd_min, a, sizeof a); jd_to_iso(ejd_max, b, sizeof b);
        printf("lowest perigee altitude %.0f km, highest apogee altitude %.0f km, epochs %s .. %s\n", pmin, pmax, a, b);
    }
    cat_free(&cat);
    return status;
}

static int cmd_gen(int argc, char **argv) {
    if (argc < 3) { fprintf(stderr, "usage: gen <n> <seed> <out.csv> [pairs]\n"); return 2; }
    int n = atoi(argv[0]);
    unsigned seed = (unsigned)strtoul(argv[1], NULL, 10);
    int pairs = argc > 3 ? atoi(argv[3]) : (n / 100 > 100 ? 100 : n / 100);
    if (pairs < 1 && n >= 2) pairs = 1;
    Catalog cat; DArray truth;
    da_init(&truth, sizeof(PlantedPair));
    if (cat_init(&cat)) return 1;
    int rc = synth_generate(&cat, &truth, n, seed, pairs);
    if (rc) { fprintf(stderr, "generation failed (%d)\n", rc); return 1; }
    char tpath[1024];
    size_t L = strlen(argv[2]);
    if (L > 4 && strcmp(argv[2] + L - 4, ".csv") == 0) snprintf(tpath, sizeof tpath, "%.*s_truth.csv", (int)(L - 4), argv[2]);
    else snprintf(tpath, sizeof tpath, "%s_truth.csv", argv[2]);
    if (cat_write_csv(&cat, argv[2]) || synth_write_truth(&truth, tpath)) { fprintf(stderr, "cannot write output\n"); return 1; }
    printf("wrote %d objects (%d planted pairs) to %s and %s\n", cat_count(&cat), pairs, argv[2], tpath);
    cat_free(&cat); da_free(&truth);
    return 0;
}

#define ALT_PAD_KM 50.0   /* empirical; 26.7 km excursion measured under SGP4 on the real catalogue (10 min), see docs/known_limitations.md L23 */
static int cmd_screen_m(int argc, char **argv, int method, int linear, double cell_scale) {
    if (argc < 8) { fprintf(stderr, "usage: screen <out.csv> <window_s> <step_s> <threshold_km> <flags (1=J2, 2=drag, 4=LEO scope, 8=SGP4)> <fresh_days> <min_mean_motion> <catalog.csv> [more.csv ...]\n"); return 2; }
    ScreenParams sp = {atof(argv[1]), atof(argv[2]), atof(argv[3]), NULL};
    int flags = atoi(argv[4]);
    int scope_leo = (flags & 4) != 0;   /* flag 4 = restrict to the project scope (LEO 175-2000 km) */
    flags &= ~4;
    double fresh = atof(argv[5]), min_mm = atof(argv[6]);
    Catalog all, cat, tmp;
    if (cat_init(&all) || cat_init(&cat) || cat_init(&tmp)) return 1;
    for (int i = 7; i < argc; i++) {
        LoadReport rep; report_init(&rep);
        if (cat_load_csv(&all, argv[i], &rep) != 0) { fprintf(stderr, "cannot load %s\n", argv[i]); return 1; }
        report_free(&rep);
    }
    double ref = cat_median_epoch(&all);
    if (fresh > 0) cat_filter_epoch(&all, &tmp, ref - fresh, ref + fresh);
    else cat_filter_epoch(&all, &tmp, 0.0, 1e12);
    if (scope_leo) {
        Catalog t2; if (cat_init(&t2)) return 1;
        cat_filter_mean_motion(&tmp, &t2, min_mm);
        cat_filter_altitude(&t2, &cat, SCOPE_LEO_MIN_ALT_KM, SCOPE_LEO_MAX_ALT_KM);
        cat_free(&t2);
    } else cat_filter_mean_motion(&tmp, &cat, min_mm);
    int n = cat_count(&cat);
    char iso[40]; jd_to_iso(ref, iso, sizeof iso);
    printf("reference epoch %s, %d of %d objects used (fresh_days=%g, min mean motion %g), flags=%d%s, window %.0f s, step %.1f s, threshold %.2f km\n",
           iso, n, cat_count(&all), fresh, min_mm, flags, scope_leo ? "+LEO scope 175-2000 km" : "", sp.window_s, sp.step_s, sp.threshold_km);
    Prop *props = malloc(sizeof(Prop) * (size_t)(n > 0 ? n : 1));
    if (!props) { fprintf(stderr, "out of memory\n"); return 1; }
    for (int i = 0; i < n; i++) prop_init(&props[i], &cat_objects(&cat)[i], flags, ref);
    double *infl = NULL; const char *ag = getenv("DLSCMS_AGE_GROWTH");     /* km per day of element age, added per object to the screening distance */
    if (ag && atof(ag) > 0.0) {
        infl = malloc(sizeof(double) * (size_t)(n > 0 ? n : 1));
        if (!infl) { fprintf(stderr, "out of memory\n"); return 1; }
        for (int i = 0; i < n; i++) { double age = ref - cat_objects(&cat)[i].epoch_jd; infl[i] = atof(ag) * (age < 0 ? -age : age); }
        sp.infl_km = infl;
        printf("age-aware screening distance: +%.3g km per day of element age per object (DLSCMS_AGE_GROWTH)\n", atof(ag));
    }
    int parked = 0; for (int i = 0; i < n; i++) if (!prop_ok(&props[i])) parked++;
    if (parked) printf("WARNING: %d objects cannot be propagated by SGP4 (deep space or invalid) and are excluded from pairing\n", parked);
    DArray ev; da_init(&ev, sizeof(Conjunction));
    ScreenStats st;
    clock_t t0 = clock();
    FastScreen fs; memset(&fs, 0, sizeof fs);
    double *rmin = NULL, *rmax = NULL;
    int rc;
    if (method == 0) rc = conj_bruteforce(props, n, &sp, &ev, &st);
    else {
        fs.method = (ScreenMethod)method; fs.refine_linear = linear; fs.cell_scale = cell_scale;
        if (method == SCREEN_ALT) {
            rmin = malloc(sizeof(double) * (size_t)(n > 0 ? n : 1)); rmax = malloc(sizeof(double) * (size_t)(n > 0 ? n : 1));
            if (!rmin || !rmax) { fprintf(stderr, "out of memory\n"); return 1; }
            for (int i = 0; i < n; i++) { rmin[i] = cat_objects(&cat)[i].perigee_r_km; rmax[i] = cat_objects(&cat)[i].apogee_r_km; }
            fs.rmin_km = rmin; fs.rmax_km = rmax; fs.alt_pad_km = ALT_PAD_KM;
        }
        rc = conj_fast(props, n, &sp, &fs, &ev, &st);
    }
    if (rc == -4) { fprintf(stderr, "screening: invalid parameters (window >= 0, step > 0, threshold >= 0 and at most 2e9 steps; grid cell scale 0.05..100)\n"); return 1; }
    if (rc) { fprintf(stderr, "screening failed (%d)\n", rc); return 1; }
    if (method == SCREEN_ALT) printf("altitude filter: pad %.0f km, largest excursion outside [perigee, apogee] radius seen at step starts: %.2f km\n", fs.alt_pad_km, fs.alt_violation_km);
    free(rmin); free(rmax); free(infl);
    double secs = (double)(clock() - t0) / CLOCKS_PER_SEC;
    DArray best; da_init(&best, sizeof(Conjunction));
    conj_best_per_pair(&ev, &best);
    if (conj_write_csv(&ev, cat_objects(&cat), argv[0])) { fprintf(stderr, "cannot write %s\n", argv[0]); return 1; }
    printf("pair checks %ld, candidates after proven bound %ld, events %ld on %zu pairs, %ld pairs already within threshold at t=0, %.2f s CPU\n", st.pair_checks, st.candidates, st.events, best.n, st.close_at_start, secs);
    {   /* peak resident memory of this process (Linux only; 0 elsewhere) */
        long kb = 0; FILE *sf = fopen("/proc/self/status", "r");
        if (sf) { char ln[256]; while (fgets(ln, sizeof ln, sf)) if (sscanf(ln, "VmHWM: %ld kB", &kb) == 1) break; fclose(sf); }
        printf("PEAKRSS_KB,%ld\n", kb);
    }
    printf("SUMMARY,%s,%d,%.0f,%.1f,%.2f,%d,%ld,%ld,%ld,%ld,%zu,%ld,%.3f\n", argv[7], n, sp.window_s, sp.step_s, sp.threshold_km, flags,
           st.pair_checks, st.candidates, st.refinements, st.events, best.n, st.close_at_start, secs);
    da_free(&best); da_free(&ev); free(props); cat_free(&cat); cat_free(&tmp); cat_free(&all);
    return 0;
}

static int cmd_screen(int argc, char **argv) { return cmd_screen_m(argc, argv, 0, 0, 1.0); }

static int cmd_prop(int argc, char **argv) {
    if (argc < 2) { fprintf(stderr, "usage: prop <catalog.csv> <dt_s> [flags: 0|1|2|3] [common]\n"); return 2; }
    Catalog cat; LoadReport rep;
    if (cat_init(&cat)) return 1;
    report_init(&rep);
    if (cat_load_csv(&cat, argv[0], &rep) != 0) { fprintf(stderr, "cannot load %s\n", argv[0]); return 1; }
    double dt = atof(argv[1]);
    int flags = argc > 2 ? atoi(argv[2]) : 0;
    int common = argc > 3 && strcmp(argv[3], "common") == 0;        /* time measured from the median epoch instead of each object's own */
    double ref = common ? cat_median_epoch(&cat) : 0.0;
    const Object *o = cat_objects(&cat);
    printf("norad_id,x_km,y_km,z_km,vx_kms,vy_kms,vz_kms\n");
    for (int i = 0; i < cat_count(&cat); i++) {
        Prop p; Vec3 r, v;
        prop_init(&p, &o[i], flags, common ? ref : o[i].epoch_jd);
        prop_state(&p, dt, &r, &v);
        printf("%d,%.6f,%.6f,%.6f,%.9f,%.9f,%.9f\n", o[i].norad_id, r.x, r.y, r.z, v.x, v.y, v.z);
    }
    report_free(&rep); cat_free(&cat);
    return 0;
}

int main(int argc, char **argv) {
    if (argc >= 3 && strcmp(argv[1], "info") == 0) return cmd_info(argc - 2, argv + 2);
    if (argc >= 2 && strcmp(argv[1], "gen") == 0) return cmd_gen(argc - 2, argv + 2);
    if (argc >= 3 && strcmp(argv[1], "prop") == 0) return cmd_prop(argc - 2, argv + 2);
    if (argc >= 10 && strcmp(argv[1], "screen") == 0) return cmd_screen(argc - 2, argv + 2);
    if (argc >= 7 && strcmp(argv[1], "analyze") == 0)                 /* analyze <events.csv> <movable.csv|-> <K> <out_prefix> <exact_node_limit> */
        { if (argc >= 8 && strcmp(argv[7], "noreduce") == 0) vc_use_reductions = 0;
          return analyze_events(argv[2], strcmp(argv[3], "-") == 0 ? NULL : argv[3], atoi(argv[4]), argv[5], atol(argv[6])); }
    if (argc >= 11 && strcmp(argv[1], "fastscreen") == 0) {       /* fastscreen <kd|grid|alt> <same arguments as screen> */
        /* method string: kd | grid | alt, optionally followed by ":linear" or ":hybrid" and/or ":<cell scale>" (e.g. grid:0.5, kd:linear) */
        char spec[64]; strncpy(spec, argv[2], sizeof spec - 1); spec[sizeof spec - 1] = 0;
        char *tok = strtok(spec, ":"); int m = tok && strcmp(tok, "kd") == 0 ? SCREEN_KDTREE : tok && strcmp(tok, "grid") == 0 ? SCREEN_GRID : tok && strcmp(tok, "alt") == 0 ? SCREEN_ALT : 0;
        int linear = 0; double cs = 1.0;
        while (m && (tok = strtok(NULL, ":")) != NULL) { if (strcmp(tok, "linear") == 0) linear = 1; else if (strcmp(tok, "hybrid") == 0) linear = 2; else if (atof(tok) > 0.0) cs = atof(tok); else m = 0; }
        if (m) return cmd_screen_m(argc - 3, argv + 3, m, linear, cs);
    }
    /* No arguments: print "build OK" + usage and succeed (scripts/build.sh run uses this as a smoke check). Anything else that
       reached this point is an unknown or incomplete command: usage on stderr, exit status 2. */
    FILE *uo = argc == 1 ? stdout : stderr;
    fprintf(uo, "Space-debris conjunction screening: build OK\n"
           "usage:\n  %s info <catalog.csv> [more.csv ...]\n  %s gen <n> <seed> <out.csv> [pairs]\n  %s prop <catalog.csv> <dt_s> [flags]\n  %s analyze <events.csv> <movable.csv|-> <K> <out_prefix> <exact_node_limit> [noreduce]\n  %s fastscreen <kd|grid|alt> <same arguments as screen>\n  %s screen <out.csv> <window_s> <step_s> <threshold_km> <flags (1=J2, 2=drag, 4=LEO scope, 8=SGP4)> <fresh_days> <min_mean_motion> <catalog.csv> [more.csv ...]\n", argv[0], argv[0], argv[0], argv[0], argv[0], argv[0]);
    return argc == 1 ? 0 : 2;
}
