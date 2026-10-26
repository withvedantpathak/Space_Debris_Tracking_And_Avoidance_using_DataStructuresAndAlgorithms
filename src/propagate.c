#include <math.h>
#include "propagate.h"

double kepler_E(double M, double e) {
    M = fmod(M, TWO_PI_);
    if (M < 0) M += TWO_PI_;
    double E = e < 0.8 ? M + e * sin(M) : 3.141592653589793;     /* standard starting guesses */
    for (int i = 0; i < 50; i++) {
        double f = E - e * sin(E) - M, d = 1.0 - e * cos(E);
        double dE = f / d;
        E -= dE;
        if (fabs(dE) < 1e-14) break;
    }
    return E;
}

/* Undo the Kozai mean-motion convention (Spacetrack Report #3 / SGP4 initialisation).
 * TLE/OMM mean motion n'' already contains the secular J2 effect on the mean anomaly; the
 * Keplerian motion underneath is n0 = n''/(1+d0) with semi-major axis a0 = (mu/n0^2)^(1/3).
 * Returns a0 (km) and n0 (rad/s). */
static void unkozai(const Object *o, double *a0, double *n0) {
    double d2r = TWO_PI_ / 360.0, nk = o->mean_motion * TWO_PI_ / 86400.0;
    double ci = cos(o->incl * d2r), c3 = 3.0 * ci * ci - 1.0, e2 = o->ecc * o->ecc;
    double q = 0.75 * EARTH_J2 * c3 / pow(1.0 - e2, 1.5);
    double a1 = pow(EARTH_MU_KM3S2 / (nk * nk), 1.0 / 3.0);
    double d1 = q * (EARTH_RADIUS_KM / a1) * (EARTH_RADIUS_KM / a1);
    double a = a1 * (1.0 - d1 / 3.0 - d1 * d1 - 134.0 / 81.0 * d1 * d1 * d1);
    double d0 = q * (EARTH_RADIUS_KM / a) * (EARTH_RADIUS_KM / a);
    *n0 = nk / (1.0 + d0);
    *a0 = pow(EARTH_MU_KM3S2 / (*n0 * *n0), 1.0 / 3.0);
}

void prop_j2_rates(const Object *o, double *raan_dot, double *argp_dot, double *m_dot) {
    double d2r = TWO_PI_ / 360.0, a0, n0;
    unkozai(o, &a0, &n0);
    double e2 = o->ecc * o->ecc, p = a0 * (1.0 - e2), ci = cos(o->incl * d2r);
    double k = EARTH_J2 * (EARTH_RADIUS_KM / p) * (EARTH_RADIUS_KM / p);
    *raan_dot = -1.5 * n0 * k * ci;
    *argp_dot = 0.75 * n0 * k * (5.0 * ci * ci - 1.0);
    *m_dot = o->mean_motion * TWO_PI_ / 86400.0;      /* Kozai mean motion IS the mean-anomaly rate */
}

void prop_init(Prop *p, const Object *o, int flags, double ref_jd) {
    p->sgp4_on = (flags & PROP_SGP4) != 0;
    if (p->sgp4_on) flags |= PROP_J2 | PROP_DRAG;      /* keep the Kepler fields (a, e) meaningful for the speed bound */
    p->sgp4_ok = 0; p->sg_t0_min = 0; p->park = 1.0e9;
    double d2r = TWO_PI_ / 360.0, dt0 = (ref_jd - o->epoch_jd) * 86400.0, a0, n0;
    p->a = o->sma_km; p->e = o->ecc; p->sqrt1me2 = sqrt(1.0 - o->ecc * o->ecc);
    p->cosi = cos(o->incl * d2r); p->sini = sin(o->incl * d2r);
    p->n = o->mean_motion * TWO_PI_ / 86400.0;
    p->m_rate = p->n; p->raan_rate = 0; p->argp_rate = 0;
    if (flags & PROP_J2) {
        unkozai(o, &a0, &n0);
        p->a = a0;
        prop_j2_rates(o, &p->raan_rate, &p->argp_rate, &p->m_rate);
    }
    /* Drag: TLE convention M(T) = M0 + n T + (ndot/2) T^2 with MEAN_MOTION_DOT = the (ndot/2) field.
     * Checked against SGP4 on the real catalogue: reading the field as ndot instead (factor 2) is much worse. */
    p->m_quad = (flags & PROP_DRAG) ? o->mean_motion_dot * TWO_PI_ / (86400.0 * 86400.0) : 0.0;
    p->M0 = o->mean_anom * d2r + p->m_rate * dt0 + p->m_quad * dt0 * dt0;   /* advance element set to ref_jd */
    p->m_rate += 2.0 * p->m_quad * dt0;                         /* d/dt of the quadratic at ref time */
    p->raan0 = o->raan * d2r + p->raan_rate * dt0;
    p->argp0 = o->argp * d2r + p->argp_rate * dt0;
    if (p->sgp4_on) {
        p->sgp4_ok = (sgp4_init(&p->sg, o) == SGP4_OK);
        p->sg_t0_min = (ref_jd - o->epoch_jd) * 1440.0;
        p->park = 1.0e9 + 1000.0 * (double)o->norad_id;
    }
}

int prop_ok(const Prop *p) { return !p->sgp4_on || p->sgp4_ok; }

void prop_state(const Prop *p, double t, Vec3 *r, Vec3 *v) {
    if (p->sgp4_on) {
        if (p->sgp4_ok && sgp4_state(&p->sg, p->sg_t0_min + t / 60.0, r, v) == SGP4_OK) return;
        r->x = p->park; r->y = 0; r->z = 0;            /* unusable object: far away, at rest, distinct per object */
        if (v) { v->x = 0; v->y = 0; v->z = 0; }
        return;
    }
    double E = kepler_E(p->M0 + p->m_rate * t + p->m_quad * t * t, p->e);
    double cE = cos(E), sE = sin(E);
    double xp = p->a * (cE - p->e), yp = p->a * p->sqrt1me2 * sE;       /* perifocal position */
    double O = p->raan0 + p->raan_rate * t, w = p->argp0 + p->argp_rate * t;
    double cO = cos(O), sO = sin(O), cw = cos(w), sw = sin(w), ci = p->cosi, si = p->sini;
    /* rotation perifocal -> inertial: Rz(O) Rx(i) Rz(w) */
    double r11 = cO * cw - sO * sw * ci, r12 = -cO * sw - sO * cw * ci;
    double r21 = sO * cw + cO * sw * ci, r22 = -sO * sw + cO * cw * ci;
    double r31 = sw * si,                r32 = cw * si;
    r->x = r11 * xp + r12 * yp; r->y = r21 * xp + r22 * yp; r->z = r31 * xp + r32 * yp;
    if (v) {
        /* Exact time derivative of the position model above: eccentric-anomaly rate from the
         * (possibly drag-accelerated) mean-anomaly rate, plus the rotation of the orbit by the
         * secular RAAN and argument-of-perigee drift. Keeps v consistent with r (needed for
         * rigorous closest-approach bracketing). */
        double Edot = (p->m_rate + 2.0 * p->m_quad * t) / (1.0 - p->e * cE);
        double vx = -p->a * sE * Edot, vy = p->a * p->sqrt1me2 * cE * Edot;
        v->x = r11 * vx + r12 * vy; v->y = r21 * vx + r22 * vy; v->z = r31 * vx + r32 * vy;
        if (p->raan_rate != 0.0 || p->argp_rate != 0.0) {
            double hx = sO * si, hy = -cO * si, hz = ci;                  /* orbit normal */
            v->x += p->raan_rate * (-r->y) + p->argp_rate * (hy * r->z - hz * r->y);
            v->y += p->raan_rate * (r->x)  + p->argp_rate * (hz * r->x - hx * r->z);
            v->z += p->argp_rate * (hx * r->y - hy * r->x);
        }
    }
}

void prop_positions(const Prop *props, int n, double t, Vec3 *out) {
    for (int i = 0; i < n; i++) prop_state(&props[i], t, &out[i], NULL);
}
