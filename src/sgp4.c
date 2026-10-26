#include <math.h>
#include <stddef.h>
#include "sgp4.h"

/* WGS-72 */
#define RE   6378.135
#define MUE  398600.8
#define J2   0.001082616
#define J3   (-0.00000253881)
#define J4   (-0.00000165597)
#define J3OJ2 (J3 / J2)
#define TWOPI (2.0 * 3.14159265358979323846)
#define X2O3 (2.0 / 3.0)

static double xke_(void) { return 60.0 / sqrt(RE * RE * RE / MUE); }

int sgp4_init(Sgp4 *s, const Object *o) {
    const double xke = xke_();
    const double deg = 3.14159265358979323846 / 180.0;
    double no_kozai = o->mean_motion * TWOPI / 1440.0;
    double ecco = o->ecc, inclo = o->incl * deg, nodeo = o->raan * deg, argpo = o->argp * deg, mo = o->mean_anom * deg;
    if (!(ecco >= 0.0 && ecco < 1.0) || !(no_kozai > 0.0)) return SGP4_ERR_INPUT;

    /* initl: un-Kozai the mean motion */
    double eccsq = ecco * ecco, omeosq = 1.0 - eccsq, rteosq = sqrt(omeosq);
    double cosio = cos(inclo), cosio2 = cosio * cosio, sinio = sin(inclo);
    double ak = pow(xke / no_kozai, X2O3);
    double d1 = 0.75 * J2 * (3.0 * cosio2 - 1.0) / (rteosq * omeosq);
    double del = d1 / (ak * ak);
    double adel = ak * (1.0 - del * del - del * (1.0 / 3.0 + 134.0 * del * del / 81.0));
    del = d1 / (adel * adel);
    double no = no_kozai / (1.0 + del);
    double ao = pow(xke / no, X2O3);
    double po = ao * omeosq;
    double con42 = 1.0 - 5.0 * cosio2;
    double con41 = -con42 - cosio2 - cosio2;
    double posq = po * po;
    double rp = ao * (1.0 - ecco);
    if (TWOPI / no >= 225.0) return SGP4_ERR_DEEP_SPACE;   /* period in minutes (un-Kozai'd, as SGP4 decides) */

    s->isimp = (rp < 220.0 / RE + 1.0) ? 1 : 0;
    double ss = 78.0 / RE + 1.0, qzms2t = pow((120.0 - 78.0) / RE, 4.0);
    double sfour = ss, qzms24 = qzms2t;
    double perige = (rp - 1.0) * RE;
    if (perige < 156.0) {
        sfour = perige - 78.0;
        if (perige < 98.0) sfour = 20.0;
        qzms24 = pow((120.0 - sfour) / RE, 4.0);
        sfour = sfour / RE + 1.0;
    }
    double pinvsq = 1.0 / posq;
    double tsi = 1.0 / (ao - sfour);
    s->eta = ao * ecco * tsi;
    double etasq = s->eta * s->eta, eeta = ecco * s->eta, psisq = fabs(1.0 - etasq);
    double coef = qzms24 * pow(tsi, 4.0);
    double coef1 = coef / pow(psisq, 3.5);
    double cc2 = coef1 * no * (ao * (1.0 + 1.5 * etasq + eeta * (4.0 + etasq)) +
                 0.375 * J2 * tsi / psisq * con41 * (8.0 + 3.0 * etasq * (8.0 + etasq)));
    s->cc1 = o->bstar * cc2;
    double cc3 = 0.0;
    if (ecco > 1.0e-4) cc3 = -2.0 * coef * tsi * J3OJ2 * no * sinio / ecco;
    s->x1mth2 = 1.0 - cosio2;
    s->cc4 = 2.0 * no * coef1 * ao * omeosq *
             (s->eta * (2.0 + 0.5 * etasq) + ecco * (0.5 + 2.0 * etasq) -
              J2 * tsi / (ao * psisq) * (-3.0 * con41 * (1.0 - 2.0 * eeta + etasq * (1.5 - 0.5 * eeta)) +
              0.75 * s->x1mth2 * (2.0 * etasq - eeta * (1.0 + etasq)) * cos(2.0 * argpo)));
    s->cc5 = 2.0 * coef1 * ao * omeosq * (1.0 + 2.75 * (etasq + eeta) + eeta * etasq);
    double cosio4 = cosio2 * cosio2;
    double temp1 = 1.5 * J2 * pinvsq * no;
    double temp2 = 0.5 * temp1 * J2 * pinvsq;
    double temp3 = -0.46875 * J4 * pinvsq * pinvsq * no;
    s->mdot = no + 0.5 * temp1 * rteosq * con41 + 0.0625 * temp2 * rteosq * (13.0 - 78.0 * cosio2 + 137.0 * cosio4);
    s->argpdot = -0.5 * temp1 * con42 + 0.0625 * temp2 * (7.0 - 114.0 * cosio2 + 395.0 * cosio4) +
                 temp3 * (3.0 - 36.0 * cosio2 + 49.0 * cosio4);
    double xhdot1 = -temp1 * cosio;
    s->nodedot = xhdot1 + (0.5 * temp2 * (4.0 - 19.0 * cosio2) + 2.0 * temp3 * (3.0 - 7.0 * cosio2)) * cosio;
    s->omgcof = o->bstar * cc3 * cos(argpo);
    s->xmcof = (ecco > 1.0e-4) ? -X2O3 * coef * o->bstar / eeta : 0.0;
    s->nodecf = 3.5 * omeosq * xhdot1 * s->cc1;
    s->t2cof = 1.5 * s->cc1;
    if (fabs(cosio + 1.0) > 1.5e-12) s->xlcof = -0.25 * J3OJ2 * sinio * (3.0 + 5.0 * cosio) / (1.0 + cosio);
    else s->xlcof = -0.25 * J3OJ2 * sinio * (3.0 + 5.0 * cosio) / 1.5e-12;
    s->aycof = -0.5 * J3OJ2 * sinio;
    double delmotemp = 1.0 + s->eta * cos(mo);
    s->delmo = delmotemp * delmotemp * delmotemp;
    s->sinmao = sin(mo);
    s->x7thm1 = 7.0 * cosio2 - 1.0;
    s->d2 = s->d3 = s->d4 = s->t3cof = s->t4cof = s->t5cof = 0.0;
    if (s->isimp != 1) {
        double cc1sq = s->cc1 * s->cc1;
        s->d2 = 4.0 * ao * tsi * cc1sq;
        double temp = s->d2 * tsi * s->cc1 / 3.0;
        s->d3 = (17.0 * ao + sfour) * temp;
        s->d4 = 0.5 * temp * ao * tsi * (221.0 * ao + 31.0 * sfour) * s->cc1;
        s->t3cof = s->d2 + 2.0 * cc1sq;
        s->t4cof = 0.25 * (3.0 * s->d3 + s->cc1 * (12.0 * s->d2 + 10.0 * cc1sq));
        s->t5cof = 0.2 * (3.0 * s->d4 + 12.0 * s->cc1 * s->d3 + 6.0 * s->d2 * s->d2 +
                   15.0 * cc1sq * (2.0 * s->d2 + cc1sq));
    }
    s->ecco = ecco; s->inclo = inclo; s->nodeo = nodeo; s->argpo = argpo; s->mo = mo; s->no = no;
    s->bstar = o->bstar; s->con41 = con41;
    return SGP4_OK;
}

int sgp4_state(const Sgp4 *s, double t, Vec3 *r, Vec3 *v) {
    const double xke = xke_(), vkmpersec = RE * xke / 60.0;
    double xmdf = s->mo + s->mdot * t, argpdf = s->argpo + s->argpdot * t, nodedf = s->nodeo + s->nodedot * t;
    double argpm = argpdf, mm = xmdf, t2 = t * t;
    double nodem = nodedf + s->nodecf * t2;
    double tempa = 1.0 - s->cc1 * t, tempe = s->bstar * s->cc4 * t, templ = s->t2cof * t2;
    if (s->isimp != 1) {
        double delomg = s->omgcof * t;
        double delmtemp = 1.0 + s->eta * cos(xmdf);
        double delm = s->xmcof * (delmtemp * delmtemp * delmtemp - s->delmo);
        double temp = delomg + delm;
        mm = xmdf + temp;
        argpm = argpdf - temp;
        double t3 = t2 * t, t4 = t3 * t;
        tempa = tempa - s->d2 * t2 - s->d3 * t3 - s->d4 * t4;
        tempe = tempe + s->bstar * s->cc5 * (sin(mm) - s->sinmao);
        templ = templ + s->t3cof * t3 + t4 * (s->t4cof + t * s->t5cof);
    }
    double nm = s->no, em = s->ecco, inclm = s->inclo;
    if (nm <= 0.0) return SGP4_ERR_MEAN_MOTION;
    double am = pow(xke / nm, X2O3) * tempa * tempa;
    nm = xke / pow(am, 1.5);
    em = em - tempe;
    if (em >= 1.0 || em < -0.001) return SGP4_ERR_MEAN_ECC;
    if (em < 1.0e-6) em = 1.0e-6;
    mm = mm + s->no * templ;
    double xlm = mm + argpm + nodem;
    nodem = fmod(nodem, TWOPI); argpm = fmod(argpm, TWOPI); xlm = fmod(xlm, TWOPI);
    mm = fmod(xlm - argpm - nodem, TWOPI);
    double sinim = sin(inclm), cosim = cos(inclm);
    double ep = em, xincp = inclm, argpp = argpm, nodep = nodem, mp = mm, sinip = sinim, cosip = cosim;
    if (ep >= 1.0 || ep < 0.0) return SGP4_ERR_PERT_ECC;
    double axnl = ep * cos(argpp);
    double temp = 1.0 / (am * (1.0 - ep * ep));
    double aynl = ep * sin(argpp) + temp * s->aycof;
    double xl = mp + argpp + nodep + temp * s->xlcof * axnl;
    double u = fmod(xl - nodep, TWOPI);
    double eo1 = u, tem5 = 9999.9, sineo1 = 0, coseo1 = 0;
    for (int ktr = 1; fabs(tem5) >= 1.0e-12 && ktr <= 10; ktr++) {
        sineo1 = sin(eo1); coseo1 = cos(eo1);
        tem5 = 1.0 - coseo1 * axnl - sineo1 * aynl;
        tem5 = (u - aynl * coseo1 + axnl * sineo1 - eo1) / tem5;
        if (fabs(tem5) >= 0.95) tem5 = tem5 > 0.0 ? 0.95 : -0.95;
        eo1 += tem5;
    }
    double ecose = axnl * coseo1 + aynl * sineo1, esine = axnl * sineo1 - aynl * coseo1;
    double el2 = axnl * axnl + aynl * aynl, pl = am * (1.0 - el2);
    if (pl < 0.0) return SGP4_ERR_SEMILATUS;
    double rl = am * (1.0 - ecose), rdotl = sqrt(am) * esine / rl, rvdotl = sqrt(pl) / rl;
    double betal = sqrt(1.0 - el2);
    temp = esine / (1.0 + betal);
    double sinu = am / rl * (sineo1 - aynl - axnl * temp), cosu = am / rl * (coseo1 - axnl + aynl * temp);
    double su = atan2(sinu, cosu), sin2u = (cosu + cosu) * sinu, cos2u = 1.0 - 2.0 * sinu * sinu;
    temp = 1.0 / pl;
    double temp1 = 0.5 * J2 * temp, temp2 = temp1 * temp;
    double mrt = rl * (1.0 - 1.5 * temp2 * betal * s->con41) + 0.5 * temp1 * s->x1mth2 * cos2u;
    su = su - 0.25 * temp2 * s->x7thm1 * sin2u;
    double xnode = nodep + 1.5 * temp2 * cosip * sin2u;
    double xinc = xincp + 1.5 * temp2 * cosip * sinip * cos2u;
    double mvt = rdotl - nm * temp1 * s->x1mth2 * sin2u / xke;
    double rvdot = rvdotl + nm * temp1 * (s->x1mth2 * cos2u + 1.5 * s->con41) / xke;
    double sinsu = sin(su), cossu = cos(su), snod = sin(xnode), cnod = cos(xnode), sini = sin(xinc), cosi = cos(xinc);
    double xmx = -snod * cosi, xmy = cnod * cosi;
    double ux = xmx * sinsu + cnod * cossu, uy = xmy * sinsu + snod * cossu, uz = sini * sinsu;
    double vx = xmx * cossu - cnod * sinsu, vy = xmy * cossu - snod * sinsu, vz = sini * cossu;
    if (mrt < 1.0) return SGP4_ERR_DECAYED;
    double mr = mrt * RE;
    r->x = mr * ux; r->y = mr * uy; r->z = mr * uz;
    if (v) {
        v->x = (mvt * ux + rvdot * vx) * vkmpersec;
        v->y = (mvt * uy + rvdot * vy) * vkmpersec;
        v->z = (mvt * uz + rvdot * vz) * vkmpersec;
    }
    return SGP4_OK;
}
