#include <math.h>
#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "catalog.h"

/* ---------------- time ---------------- */
static int days_in_month(int y, int m) {
    static const int d[] = {31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31};
    if (m == 2 && ((y % 4 == 0 && y % 100 != 0) || y % 400 == 0)) return 29;
    return d[m - 1];
}

int iso_to_jd(const char *s, double *jd) {
    int Y, M, D, h = 0, mi = 0;
    double sec = 0;
    int n = sscanf(s, "%d-%d-%d%*1[T ]%d:%d:%lf", &Y, &M, &D, &h, &mi, &sec);
    if (n != 6 && n != 3) return -1;
    if (Y < 1957 || Y > 2200 || M < 1 || M > 12 || D < 1 || D > days_in_month(Y, M)) return -1;
    if (h < 0 || h > 23 || mi < 0 || mi > 59 || !(sec >= 0 && sec < 60)) return -1;
    long long a = (14 - M) / 12, y2 = Y + 4800 - a, m2 = M + 12 * a - 3;
    long long jdn = D + (153 * m2 + 2) / 5 + 365 * y2 + y2 / 4 - y2 / 100 + y2 / 400 - 32045;
    *jd = (double)jdn - 0.5 + (h * 3600.0 + mi * 60.0 + sec) / 86400.0;
    return 0;
}

void jd_to_iso(double jd, char *buf, size_t size) {
    double x = jd + 0.5;
    long long J = (long long)floor(x);
    long long us = llround((x - (double)J) * 86400e6);
    if (us >= 86400000000LL) { us -= 86400000000LL; J++; }
    long long f = J + 1401 + (((4 * J + 274277) / 146097) * 3) / 4 - 38;
    long long e = 4 * f + 3, g = (e % 1461) / 4, h = 5 * g + 2;
    int D = (int)((h % 153) / 5 + 1), M = (int)((h / 153 + 2) % 12 + 1);
    int Y = (int)(e / 1461 - 4716 + (12 + 2 - M) / 12);
    long long secs = us / 1000000;
    snprintf(buf, size, "%04d-%02d-%02dT%02d:%02d:%02d.%06d", Y, M, D,
             (int)(secs / 3600), (int)(secs / 60 % 60), (int)(secs % 60), (int)(us % 1000000));
}

void obj_derive(Object *o) {
    double n = o->mean_motion * TWO_PI_ / 86400.0;                 /* rad/s */
    o->sma_km = pow(EARTH_MU_KM3S2 / (n * n), 1.0 / 3.0);          /* Kepler's third law */
    o->perigee_r_km = o->sma_km * (1.0 - o->ecc);
    o->apogee_r_km = o->sma_km * (1.0 + o->ecc);
}

/* ---------------- catalogue ---------------- */
int cat_init(Catalog *c) {
    da_init(&c->objs, sizeof(Object));
    return hash_init(&c->by_id, 1024);
}
int cat_count(const Catalog *c) { return (int)c->objs.n; }
const Object *cat_objects(const Catalog *c) { return (const Object *)c->objs.data; }
int cat_find(const Catalog *c, int id) {
    int64_t v;
    return hash_get(&c->by_id, id, &v) ? (int)v : -1;
}
int cat_add(Catalog *c, const Object *o) {
    int64_t v;
    if (hash_get(&c->by_id, o->norad_id, &v)) {
        Object *old = (Object *)da_at(&c->objs, (size_t)v);
        if (o->epoch_jd > old->epoch_jd) { *old = *o; return 2; }   /* keep the newer element set */
        return 0;
    }
    if (da_push(&c->objs, o) != 0) return -1;
    if (hash_put(&c->by_id, o->norad_id, (int64_t)c->objs.n - 1) != 0) { c->objs.n--; return -1; }   /* keep array and index consistent */
    return 1;
}
void cat_free(Catalog *c) { da_free(&c->objs); hash_free(&c->by_id); }

int cat_filter_epoch(const Catalog *src, Catalog *dst, double jd_lo, double jd_hi) {
    const Object *o = cat_objects(src);
    int k = 0;
    for (int i = 0; i < cat_count(src); i++) {
        if (o[i].epoch_jd < jd_lo || o[i].epoch_jd > jd_hi) continue;
        if (cat_add(dst, &o[i]) < 0) return -1;
        k++;
    }
    return k;
}

int cat_filter_mean_motion(const Catalog *src, Catalog *dst, double mm_min) {
    const Object *o = cat_objects(src);
    int k = 0;
    for (int i = 0; i < cat_count(src); i++) {
        if (o[i].mean_motion < mm_min) continue;
        if (cat_add(dst, &o[i]) < 0) return -1;
        k++;
    }
    return k;
}

int cat_filter_altitude(const Catalog *src, Catalog *dst, double min_alt_km, double max_alt_km) {
    const Object *o = cat_objects(src);
    int k = 0;
    for (int i = 0; i < cat_count(src); i++) {
        if (o[i].perigee_r_km - EARTH_RADIUS_KM < min_alt_km) continue;
        if (o[i].apogee_r_km - EARTH_RADIUS_KM > max_alt_km) continue;
        if (cat_add(dst, &o[i]) < 0) return -1;
        k++;
    }
    return k;
}

static int cmp_dbl(const void *a, const void *b) { double x = *(const double *)a, y = *(const double *)b; return (x > y) - (x < y); }
double cat_median_epoch(const Catalog *c) {
    int n = cat_count(c);
    if (n == 0) return 0.0;
    double *e = malloc(sizeof(double) * (size_t)n);
    if (!e) return 0.0;
    const Object *o = cat_objects(c);
    for (int i = 0; i < n; i++) e[i] = o[i].epoch_jd;
    qsort(e, (size_t)n, sizeof(double), cmp_dbl);
    double m = e[n / 2];
    free(e);
    return m;
}

void report_init(LoadReport *r) { memset(r, 0, sizeof *r); da_init(&r->errors, sizeof(ParseError)); }
void report_free(LoadReport *r) { da_free(&r->errors); }

static void rep_err(LoadReport *r, int line, const char *fmt, ...) {
    r->errors_total++;
    r->rows_bad++;
    if (r->errors.n >= 100) return;
    ParseError e;
    e.line = line;
    va_list ap; va_start(ap, fmt);
    vsnprintf(e.msg, sizeof e.msg, fmt, ap);
    va_end(ap);
    da_push(&r->errors, &e);
}

/* ---------------- CSV parsing ---------------- */
/* Split a line in place at commas; supports "quoted, fields" and "" escapes. */
static int split_csv(char *line, char **f, int max) {
    int n = 0;
    char *p = line;
    for (;;) {
        if (n >= max) return n;
        if (*p == '"') {
            char *w = ++p;
            f[n++] = w;
            for (;;) {
                if (*p == '"' && p[1] == '"') { *w++ = '"'; p += 2; }
                else if (*p == '"') { p++; break; }
                else if (*p == '\0') break;
                else *w++ = *p++;
            }
            *w = '\0';
            while (*p && *p != ',') p++;
        } else {
            f[n++] = p;
            while (*p && *p != ',') p++;
        }
        if (*p == ',') { *p++ = '\0'; continue; }
        return n;
    }
}

static int parse_double(const char *s, double *out) {
    char *end;
    if (!*s) return -1;
    double v = strtod(s, &end);
    while (*end == ' ' || *end == '\t') end++;
    if (*end || !isfinite(v)) return -1;
    *out = v;
    return 0;
}
static double wrap360(double a) { a = fmod(a, 360.0); return a < 0 ? a + 360.0 : a; }

enum { C_NAME, C_OID, C_EPOCH, C_MM, C_ECC, C_INC, C_RAAN, C_ARGP, C_MA, C_ID, C_BSTAR, C_MMDOT, NCOLS };
static const char *COLNAMES[NCOLS] = {"OBJECT_NAME", "OBJECT_ID", "EPOCH", "MEAN_MOTION", "ECCENTRICITY",
    "INCLINATION", "RA_OF_ASC_NODE", "ARG_OF_PERICENTER", "MEAN_ANOMALY", "NORAD_CAT_ID", "BSTAR", "MEAN_MOTION_DOT"};
static const int OPTIONAL_COL[NCOLS] = {1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1};

/* Reads one line into buf (NUL-terminated, at most size-1 characters kept; '\n' and a trailing '\r' are removed).
 * Returns 1 ok, 0 at end of file. *too_long is set if the line did not fit (the rest of it is discarded).
 * *has_nul is set if the line contained a NUL byte; each NUL is replaced by a space so the whole line is rejected as
 * one bad row. (An earlier fgets/strlen version mistook such a line for an over-long one and swallowed the NEXT line.) */
static int read_line(FILE *f, char *buf, size_t size, int *too_long, int *has_nul) {
    size_t n = 0;
    int ch, got = 0;
    *too_long = 0; *has_nul = 0;
    while ((ch = fgetc(f)) != EOF) {
        got = 1;
        if (ch == '\n') break;
        if (ch == '\0') { *has_nul = 1; ch = ' '; }
        if (n + 1 < size) buf[n++] = (char)ch; else *too_long = 1;
    }
    if (!got) return 0;
    buf[n] = '\0';
    if (n && buf[n - 1] == '\r') buf[--n] = '\0';
    return 1;
}

int cat_load_csv(Catalog *c, const char *path, LoadReport *rep) {
    FILE *f = fopen(path, "r");
    if (!f) return -1;
    char line[4096], *fld[64];
    int too_long, has_nul, lineno = 0, col[NCOLS], ncols_needed = 0;
    for (int i = 0; i < NCOLS; i++) col[i] = -1;

    /* header: first non-blank line */
    int have_header = 0;
    while (!have_header && read_line(f, line, sizeof line, &too_long, &has_nul)) {
        lineno++;
        char *h = line;
        if ((unsigned char)h[0] == 0xEF && (unsigned char)h[1] == 0xBB && (unsigned char)h[2] == 0xBF) h += 3;
        if (!*h) continue;
        int nf = split_csv(h, fld, 64);
        for (int k = 0; k < nf; k++) for (int i = 0; i < NCOLS; i++) if (strcmp(fld[k], COLNAMES[i]) == 0) col[i] = k;
        have_header = 1;
    }
    rep->lines = lineno;
    int missing = 0;
    for (int i = 0; i < NCOLS; i++) {
        if (col[i] < 0 && !OPTIONAL_COL[i]) { rep_err(rep, lineno ? lineno : 1, "header lacks column %s", COLNAMES[i]); rep->rows_bad--; missing = 1; }
        if (col[i] + 1 > ncols_needed) ncols_needed = col[i] + 1;
    }
    if (!have_header || missing) { fclose(f); return -2; }

    while (read_line(f, line, sizeof line, &too_long, &has_nul)) {
        lineno++;
        rep->lines = lineno;
        if (!line[0]) continue;
        if (too_long) { rep_err(rep, lineno, "line longer than %d characters", (int)sizeof line - 1); continue; }
        if (has_nul) { rep_err(rep, lineno, "line contains a NUL byte"); continue; }
        int nf = split_csv(line, fld, 64);
        if (nf < ncols_needed) { rep_err(rep, lineno, "too few fields (%d, need %d)", nf, ncols_needed); continue; }

        Object o;
        memset(&o, 0, sizeof o);
        double id, mm, ecc, inc, raan, argp, ma, bstar = 0, mmdot = 0;
        int bad = 0;
        if (parse_double(fld[col[C_ID]], &id) || id < 1 || id > 999999999 || id != floor(id)) { rep_err(rep, lineno, "NORAD_CAT_ID invalid (\"%.20s\")", fld[col[C_ID]]); continue; }
        if (iso_to_jd(fld[col[C_EPOCH]], &o.epoch_jd)) { rep_err(rep, lineno, "EPOCH invalid (\"%.30s\")", fld[col[C_EPOCH]]); continue; }
        const char *nm[8] = {"MEAN_MOTION", "ECCENTRICITY", "INCLINATION", "RA_OF_ASC_NODE", "ARG_OF_PERICENTER", "MEAN_ANOMALY", "BSTAR", "MEAN_MOTION_DOT"};
        double *dst[8] = {&mm, &ecc, &inc, &raan, &argp, &ma, &bstar, &mmdot};
        int src[8] = {C_MM, C_ECC, C_INC, C_RAAN, C_ARGP, C_MA, C_BSTAR, C_MMDOT};
        for (int k = 0; k < 8 && !bad; k++) {
            if (col[src[k]] < 0) continue;                         /* optional column absent */
            if (parse_double(fld[col[src[k]]], dst[k])) { rep_err(rep, lineno, "%s not a number (\"%.20s\")", nm[k], fld[col[src[k]]]); bad = 1; }
        }
        if (bad) continue;
        if (ecc < 0 || ecc >= 1) { rep_err(rep, lineno, "ECCENTRICITY out of range [0,1): %g", ecc); continue; }
        if (inc < 0 || inc > 180) { rep_err(rep, lineno, "INCLINATION out of range [0,180]: %g", inc); continue; }
        if (mm <= 0 || mm > 20) { rep_err(rep, lineno, "MEAN_MOTION out of range (0,20]: %g", mm); continue; }

        o.norad_id = (int)id;
        o.mean_motion = mm; o.ecc = ecc; o.incl = inc; o.bstar = bstar; o.mean_motion_dot = mmdot;
        o.raan = wrap360(raan); o.argp = wrap360(argp); o.mean_anom = wrap360(ma);
        if (col[C_NAME] >= 0) snprintf(o.name, sizeof o.name, "%s", fld[col[C_NAME]]);
        if (col[C_OID] >= 0) snprintf(o.intl_id, sizeof o.intl_id, "%s", fld[col[C_OID]]);
        obj_derive(&o);
        if (o.perigee_r_km < EARTH_RADIUS_KM) { rep->rows_suborbital++; continue; }
        int r = cat_add(c, &o);
        if (r < 0) { fclose(f); return -3; }
        if (r == 1) rep->rows_added++; else if (r == 2) rep->rows_replaced++; else rep->rows_dup_ignored++;
    }
    fclose(f);
    return 0;
}

int cat_write_csv(const Catalog *c, const char *path) {
    FILE *f = fopen(path, "w");
    if (!f) return -1;
    fputs("OBJECT_NAME,OBJECT_ID,EPOCH,MEAN_MOTION,ECCENTRICITY,INCLINATION,RA_OF_ASC_NODE,ARG_OF_PERICENTER,"
          "MEAN_ANOMALY,EPHEMERIS_TYPE,CLASSIFICATION_TYPE,NORAD_CAT_ID,ELEMENT_SET_NO,REV_AT_EPOCH,BSTAR,"
          "MEAN_MOTION_DOT,MEAN_MOTION_DDOT\n", f);
    const Object *o = cat_objects(c);
    for (int i = 0; i < cat_count(c); i++) {
        char ep[40];
        jd_to_iso(o[i].epoch_jd, ep, sizeof ep);
        const char *nmq = strpbrk(o[i].name, ",\"") ? "\"" : "";   /* quote names containing commas */
        fprintf(f, "%s%s%s,%s,%s,%.12f,%.10f,%.10f,%.10f,%.10f,%.10f,0,U,%d,999,0,%.8g,%.10g,0\n",
                nmq, o[i].name, nmq, o[i].intl_id, ep, o[i].mean_motion, o[i].ecc, o[i].incl, o[i].raan, o[i].argp,
                o[i].mean_anom, o[i].norad_id, o[i].bstar, o[i].mean_motion_dot);
    }
    return fclose(f) == 0 ? 0 : -1;
}
