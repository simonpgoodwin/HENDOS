/* ==============================================================
 * HENDOS - HErmite N-body Dynamics Of Systems - few.c
 * FEW-BODY TOOLS: orbital elements and the stability test for a
 * hierarchical triple. These know nothing about the rest of the code.
 * Internal units (G = 1). (The small Hermite integrator that used to be
 * here for few-body groups has been replaced by the regularised
 * integrator in ar.c; it is still in the archived Fortran,
 * archive/fortran/src/nbfew.f90.)
 * ============================================================== */
#include "hendos.h"

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
