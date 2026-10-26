#ifndef CATALOG_H
#define CATALOG_H
#include <stddef.h>
#include "ds/darray.h"
#include "ds/hash.h"

/* Orbital-element catalogue: object record, CelesTrak OMM-CSV parser with per-line error
 * reporting, merge by NORAD id. Elements follow the OMM/TLE convention (mean motion in
 * rev/day, angles in degrees). Derived radii are two-body values from the mean motion. */
#define EARTH_RADIUS_KM 6378.137
#define EARTH_MU_KM3S2  398600.4418
#define TWO_PI_ 6.283185307179586

typedef struct {
    int    norad_id;
    char   name[32];
    char   intl_id[16];
    double epoch_jd;                 /* Julian date (UTC) of the element set */
    double mean_motion;              /* rev/day */
    double ecc, incl, raan, argp, mean_anom;   /* ecc dimensionless; angles in degrees, [0,360) */
    double bstar;
    double mean_motion_dot;          /* rev/day^2, the TLE/OMM "ndot/2" field (drag proxy) */
    double sma_km, perigee_r_km, apogee_r_km;  /* derived; radii from Earth's centre */
} Object;

typedef struct { int line; char msg[96]; } ParseError;
typedef struct {
    long lines, rows_added, rows_replaced, rows_dup_ignored, rows_suborbital, rows_bad;
    long errors_total;               /* all errors seen */
    DArray errors;                   /* first 100 as ParseError */
} LoadReport;

typedef struct { DArray objs; Hash by_id; } Catalog;

/* time helpers */
int  iso_to_jd(const char *s, double *jd);                 /* "YYYY-MM-DDTHH:MM:SS[.f]" -> JD; 0 ok, -1 bad */
void jd_to_iso(double jd, char *buf, size_t size);         /* microsecond resolution */

void obj_derive(Object *o);                                /* fills sma_km, perigee_r_km, apogee_r_km */

int  cat_init(Catalog *c);
int  cat_add(Catalog *c, const Object *o);                /* 1 added, 2 replaced older duplicate, 0 duplicate ignored, -1 oom */
int  cat_count(const Catalog *c);
const Object *cat_objects(const Catalog *c);               /* contiguous array of cat_count() objects */
int  cat_find(const Catalog *c, int norad_id);             /* index or -1; O(1) average */
/* Copy objects whose epoch lies in [jd_lo, jd_hi] into dst (dst must be initialised). Used to
 * keep only fresh element sets for a screening window. Returns number copied, -1 on oom. */
int  cat_filter_epoch(const Catalog *src, Catalog *dst, double jd_lo, double jd_hi);
/* Copy objects with mean motion >= mm_min rev/day (mm_min = 6.4 keeps near-earth objects, orbital period < 225 min,
 * which is the SGP4 near-earth / deep-space boundary). Returns number copied, -1 on oom. */
int  cat_filter_mean_motion(const Catalog *src, Catalog *dst, double mm_min);
/* Project scope (follows Gomez et al. 2022, sec 3.1): keep objects whose whole orbit lies between
   min_alt_km and max_alt_km above the spherical Earth (perigee >= min, apogee <= max). Returns count kept, -1 on oom. */
#define SCOPE_LEO_MIN_ALT_KM 175.0
#define SCOPE_LEO_MAX_ALT_KM 2000.0
int  cat_filter_altitude(const Catalog *src, Catalog *dst, double min_alt_km, double max_alt_km);
double cat_median_epoch(const Catalog *c);                 /* 0 if empty */
void cat_free(Catalog *c);

void report_init(LoadReport *r);
void report_free(LoadReport *r);

/* Load a CelesTrak OMM CSV. Columns are found by header name, any order. Bad rows are
 * skipped and recorded with their line number. Objects whose perigee is below the Earth's
 * surface are skipped (counted as suborbital). Returns 0, -1 file not found, -2 header lacks
 * a required column, -3 out of memory. O(rows). */
int  cat_load_csv(Catalog *c, const char *path, LoadReport *rep);

/* Write the catalogue as OMM CSV (same 17-column header as CelesTrak). */
int  cat_write_csv(const Catalog *c, const char *path);
#endif
