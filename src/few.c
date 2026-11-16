/* ==============================================================
 * HENDOS - HErmite N-body Dynamics Of Systems - few.c
 * FEW-BODY TOOLS: a small integrator for an isolated group of a few
 * stars, and the stability test for a hierarchical triple. These know
 * nothing about the rest of the code. Internal units (G = 1).
 * (Archived Fortran equivalent: archive/fortran/src/nbfew.f90.)
 * ============================================================== */
#include "hendos.h"

/* accelerations a and jerks jk of the nm stars of a group from each
 * other, the square of the shortest pairwise free-fall or flyby time,
 * and the closest pair ic, jc with its squared separation */
static void fewacc(int nm, const double *mm, vec3 *x, vec3 *vv, double e2, vec3 *a, vec3 *jk,
                   double *tmin2, double *r2min, int *ic, int *jc)
{
    memset(a, 0, sizeof(vec3) * (size_t)nm);
    memset(jk, 0, sizeof(vec3) * (size_t)nm);
    *tmin2 = HUGE_VAL;
    *r2min = HUGE_VAL;
    *ic = *jc = -1;
    for (int i = 0; i < nm - 1; i++)
        for (int j = i + 1; j < nm; j++) {
            double d[3], dv[3];
            for (int k = 0; k < 3; k++) {
                d[k] = x[j][k] - x[i][k];
                dv[k] = vv[j][k] - vv[i][k];
            }
            double d2 = dot3(d, d);
            if (d2 < *r2min) {
                *r2min = d2;
                *ic = i;
                *jc = j;
            }
            double rad2 = d2 + e2, rinv = 1.0 / sqrt(rad2), rinv3 = rinv / rad2;
            double rv = 3.0 * dot3(d, dv) / rad2;
            for (int k = 0; k < 3; k++) {
                double w = (dv[k] - rv * d[k]) * rinv3;
                a[i][k] += mm[j] * d[k] * rinv3;
                a[j][k] -= mm[i] * d[k] * rinv3;
                jk[i][k] += mm[j] * w;
                jk[j][k] -= mm[i] * w;
            }
            double tff = rad2 * rad2 * rinv / (mm[i] + mm[j]);
            double tfl = rad2 / (dot3(dv, dv) + 1e-300);
            if (tff < *tmin2) *tmin2 = tff;
            if (tfl < *tmin2) *tmin2 = tfl;
        }
}

/* integrate the isolated group of nm stars (masses mm, positions x and
 * velocities vv relative to its centre of mass) for a time tint.
 * 4th-order Hermite (predict, evaluate, correct, evaluate) with one
 * timestep shared by all the stars, dt = etaf * (shortest pairwise
 * free-fall or flyby time). It stops early if two stars come within
 * sqrt(rc2) of each other (rc2 > 0; ic, jc are then the pair and tdone
 * the time reached), or after NSTMAX steps (ic = -2). Otherwise
 * ic = -1 and tdone = tint. nst returns the number of steps. */
void fewint(int nm, const double *mm, vec3 *x, vec3 *vv, double tint, double etaf, double e2,
            double rc2, double *tdone, int *ic, int *jc, int64_t *nst)
{
    const int64_t NSTMAX = 20000000;
    vec3 a0[KFMAX], j0[KFMAX], a1[KFMAX], j1[KFMAX], xp[KFMAX], vq[KFMAX];
    double t = 0.0, tmin2, r2min;
    int k1, k2;
    *nst = 0;
    *tdone = 0.0;
    fewacc(nm, mm, x, vv, e2, a0, j0, &tmin2, &r2min, &k1, &k2);
    for (;;) {
        if (rc2 > 0.0 && r2min < rc2) {
            *ic = k1;
            *jc = k2;
            *tdone = t;
            return;
        }
        if (*nst >= NSTMAX) {
            *ic = -2;
            *jc = -1;
            *tdone = t;
            return;
        }
        if (t >= tint) break;
        double dt = etaf * sqrt(tmin2);
        int last = 0;
        if (t + dt >= tint) {
            dt = tint - t;
            last = 1;
        }
        for (int i = 0; i < nm; i++)
            for (int k = 0; k < 3; k++) {
                xp[i][k] = x[i][k] + vv[i][k] * dt + a0[i][k] * dt * dt / 2.0
                           + j0[i][k] * dt * dt * dt / 6.0;
                vq[i][k] = vv[i][k] + a0[i][k] * dt + j0[i][k] * dt * dt / 2.0;
            }
        fewacc(nm, mm, xp, vq, e2, a1, j1, &tmin2, &r2min, &k1, &k2);
        for (int i = 0; i < nm; i++)
            for (int k = 0; k < 3; k++) {
                double v1 = vv[i][k] + (a0[i][k] + a1[i][k]) * dt / 2.0
                            + (j0[i][k] - j1[i][k]) * dt * dt / 12.0;
                x[i][k] += (vv[i][k] + v1) * dt / 2.0 + (a0[i][k] - a1[i][k]) * dt * dt / 12.0;
                vv[i][k] = v1;
            }
        /* evaluate at the corrected positions: the start of the next step */
        fewacc(nm, mm, x, vv, e2, a0, j0, &tmin2, &r2min, &k1, &k2);
        (*nst)++;
        t = last ? tint : t + dt;
    }
    *ic = *jc = -1;
    *tdone = tint;
    /* keep the centre of mass exactly at the origin */
    double mt = 0.0, cx[3] = {0, 0, 0}, cvv[3] = {0, 0, 0};
    for (int i = 0; i < nm; i++) {
        mt += mm[i];
        for (int k = 0; k < 3; k++) {
            cx[k] += mm[i] * x[i][k];
            cvv[k] += mm[i] * vv[i][k];
        }
    }
    for (int i = 0; i < nm; i++)
        for (int k = 0; k < 3; k++) {
            x[i][k] -= cx[k] / mt;
            vv[i][k] -= cvv[k] / mt;
        }
}

/* semi-major axis a (negative if unbound), eccentricity and specific
 * angular momentum h of the relative orbit rr, vr with G*(mass) = mu */
void elements(double mu, const double rr[3], const double vr[3], double *a, double *ecc,
              double h[3])
{
    double e = 0.5 * dot3(vr, vr) - mu / norm3(rr);
    *a = -0.5 * mu / e;
    h[0] = rr[1] * vr[2] - rr[2] * vr[1];
    h[1] = rr[2] * vr[0] - rr[0] * vr[2];
    h[2] = rr[0] * vr[1] - rr[1] * vr[0];
    *ecc = sqrt(fmax(0.0, 1.0 + 2.0 * e * dot3(h, h) / (mu * mu)));
}

/* is the triple a stable hierarchy? The inner binary is the bound pair
 * with the smallest semi-major axis; the third star must be bound to
 * it and the outer pericentre must satisfy Mardling & Aarseth (2001):
 *   q_out/a_in > 2.8 [(1 + m3/(m1+m2)) (1 + e_out)/sqrt(1 - e_out)]^(2/5)
 *                    (1 - 0.3 i/pi)
 * Returns the inner pair i1, i2 and outer star i3, the inner orbit rel
 * (star i1 relative to i2), the outer orbit rel2 (star i3 relative to
 * the inner c.o.m.) and the two pericentres. */
void hierarchy(const double mm[3], vec3 x[3], vec3 vv[3], int *stable, int *i1, int *i2,
               int *i3, double rel[6], double rel2[6], double *qin, double *qout)
{
    static const int pa[3] = {0, 0, 1}, pb[3] = {1, 2, 2}, pc[3] = {2, 1, 0};
    double a, ecc, h[3], abest = HUGE_VAL, d[3], dv[3];
    *stable = 0;
    *i1 = *i2 = *i3 = -1;
    *qin = *qout = 0.0;
    memset(rel, 0, 6 * sizeof(double));
    memset(rel2, 0, 6 * sizeof(double));
    for (int p = 0; p < 3; p++) {
        for (int k = 0; k < 3; k++) {
            d[k] = x[pa[p]][k] - x[pb[p]][k];
            dv[k] = vv[pa[p]][k] - vv[pb[p]][k];
        }
        elements(mm[pa[p]] + mm[pb[p]], d, dv, &a, &ecc, h);
        if (a > 0.0 && a < abest) {
            abest = a;
            *i1 = pa[p];
            *i2 = pb[p];
            *i3 = pc[p];
        }
    }
    if (*i1 < 0) return;
    double mb = mm[*i1] + mm[*i2], mt = mb + mm[*i3];
    double ain, ein, hin[3], aout, eout, hout[3];
    for (int k = 0; k < 3; k++) {
        rel[k] = x[*i1][k] - x[*i2][k];
        rel[k + 3] = vv[*i1][k] - vv[*i2][k];
    }
    elements(mb, rel, rel + 3, &ain, &ein, hin);
    for (int k = 0; k < 3; k++) {
        rel2[k] = x[*i3][k] - (mm[*i1] * x[*i1][k] + mm[*i2] * x[*i2][k]) / mb;
        rel2[k + 3] = vv[*i3][k] - (mm[*i1] * vv[*i1][k] + mm[*i2] * vv[*i2][k]) / mb;
    }
    elements(mt, rel2, rel2 + 3, &aout, &eout, hout);
    if (aout <= 0.0 || eout >= 1.0) return;
    *qin = ain * (1.0 - ein);
    *qout = aout * (1.0 - eout);
    double cosi = dot3(hin, hout) / sqrt(dot3(hin, hin) * dot3(hout, hout) + 1e-300);
    double inc = acos(fmax(-1.0, fmin(1.0, cosi)));
    double crit = 2.8 * pow((1.0 + mm[*i3] / mb) * (1.0 + eout) / sqrt(1.0 - eout), 0.4)
                  * (1.0 - 0.3 * inc / 3.141592653589793);
    *stable = (*qout > crit * ain);
}
