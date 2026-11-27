/* ==============================================================
 * HENDOS - HErmite N-body Dynamics Of Systems - groups.c
 * BINARIES, FEW-BODY GROUPS, COLLISIONS AND MERGERS
 *
 * GROUPS. A compact set of stars is replaced in the main integration by
 * a single centre-of-mass (c.o.m.) body, and its internal motion is
 * followed separately. The stars then no longer force tiny timesteps on
 * the main integration (where every step costs a prediction of all N
 * bodies). There are three kinds of group (bmode):
 *   1  unperturbed binary   an exact Kepler orbit (costs nothing)
 *   2  regularised group    2 to KFMAX stars - a perturbed binary, or a
 *                           compact multiple in a close encounter -
 *                           integrated by algorithmic regularisation
 *                           with the tidal pull of its perturbers, and
 *                           the perturbers feel its stars (ar.c)
 *   3  frozen triple        a hierarchical triple that passes the
 *                           Mardling & Aarseth (2001) stability test and
 *                           is effectively unperturbed, as two nested
 *                           Kepler orbits (no secular evolution)
 * When a group becomes too perturbed it is 'opened' back into its stars,
 * which are then integrated directly like any other stars.
 *
 * relative perturbation of a group of mass mb and size s by all other
 * bodies k at distance d_k from its c.o.m.:
 *    gamma = 2 s^3 / mb * sum_k m_k / d_k^3
 * (tidal acceleration across the group / internal acceleration). Two
 * sizes are used: g_orb, with s the apocentre of a bound pair (is it
 * unperturbed ALL round its orbit? - for Kepler binaries), and g_cur,
 * with s its current separation, or the diameter of a larger group (is
 * it isolated NOW? - for regularised groups, as in KS regularisation).
 *
 * Star i, after its step, is a candidate if its step < dtclose and its
 * nearest neighbour j is on the same block. Then:
 * CLOSE   i, j single, bound, g_orb < gclose              -> 1 Kepler
 *         i, j single, bound or plunging (unbound and coming together,
 *         pericentre < separation / 1000), g_cur < gfew   -> 2 (a pair)
 * FEW     otherwise the pair and the next nearest body are tried as a
 *         group; so is a star and a nearby group. Made if its own
 *         g_cur < gfew and it pays (see tryfew)            -> 2
 * Groups, after each step of their c.o.m.:
 *   1 Kepler:  g_orb > gopen      -> a regularised pair (2)
 *   2 pair:    g_orb < gclose     -> Kepler (1), if bound
 *   2 triple:  g_cur < gclose and a stable hierarchy -> frozen (3)
 *   3 frozen:  g_cur > gopen      -> regularised (2)
 *   any:       g_cur > 10 gfew    -> takes its nearest body in, if the
 *                                    result is a group, or is opened
 *   2 pair:    unbound, moving apart, and twice as far apart as when it
 *              was made, or g_cur > gfew -> opened
 *   2 any:     perturber list full and its weakest perturber still strong
 *                                 -> opened (too big to regularise)
 *   2 any:     its nearest body is a group perturbing it by more than
 *              gfew on its own   -> takes it in, if the result is a group,
 *                                   or else is opened
 * (gopen > gclose and the factor 10 give hysteresis, so groups do not
 * flicker.) Binaries in the initial conditions are closed at the start.
 * A star that leaves a group makes the group large, so its gamma rises
 * and it is opened; what is left is closed again if it is a binary.
 * With fewbody off there are no regularised groups: perturbed binaries
 * are integrated directly, and only Kepler binaries are made.
 *
 * COLLISIONS AND MERGERS. Two stars that pass within rcoll of each
 * other have collided. With mergers on they are merged into one star
 * (mass and momentum conserved) and the run carries on; a merger is
 * announced on the screen and written to <runid>_mergers.dat. The
 * orbital energy of the pair is lost from the stellar dynamics: it is
 * added up in emerged so that the energy check still works. With
 * mergers off the run stops with an alert, as before. Collisions inside
 * regularised groups are treated in the same way (ar.c).
 *
 * The bodies are kept contiguous in slots 0..n-1 of the particle arrays
 * (n = number of bodies, ntot = number of stars), so the force loops
 * are unchanged: closing puts the c.o.m. in the lowest slot and moves
 * the last bodies into the others; opening uses the c.o.m. slot and
 * slots n onwards. Only one open and one close are done per block.
 * bcm[k] is the slot of group k's c.o.m., and the perturber lists of the
 * regularised groups are kept pointing at the right slots (premap).
 *
 * All of this is done serially between blocks, outside any OpenMP
 * parallel region; the O(N) loops inside are parallelised. With the
 * neighbour scheme on (ac.c) the neighbour lists are kept correct
 * through every change.
 *
 * Groups are invisible in the output: snapshots contain all ntot
 * stars. The types of the stars are kept in btype while they are in a
 * group (the c.o.m. body has itype = 0) and restored when it is opened.
 *
 * Here a group index is 0..nbin-1, and ib[i] = -1 means body i is a
 * single star.
 * ============================================================== */
#include <omp.h>
#include "hendos.h"

#define NFEWFAC 2    /* a multiple of nm stars needs NFEWFAC*nm*(nm-1) bodies to pay */
#define PLUNGE 1.0e-3  /* >>> CHANGE HERE: an unbound pair is regularised only if its
                          pericentre < PLUNGE * its separation (see pairstate) */

static void opengroup(int i, int64_t t8, double tnow, int *nm, int *list);
static int bmode_last_open = 0;


/* ---- debugging: total energy at integer time t8 from the predicted
 * positions (valid between the steps of a block, when every body is
 * at t8). Used only when the environment variable HENDOS_DEBUG is set,
 * to check that no group event changes the energy. */
static int dbgon = -1;
double dbgenergy(double tnow)
{
    double ke = 0.0, pe = 0.0;
#pragma omp parallel for reduction(+:ke,pe) schedule(dynamic, 16)
    for (int i = 0; i < n; i++) {
        ke += 0.5 * m[i] * dot3(vp[i], vp[i]);
        for (int j = i + 1; j < n; j++) {
            double d[3] = {rp[i][0] - rp[j][0], rp[i][1] - rp[j][1], rp[i][2] - rp[j][2]};
            pe -= m[i] * m[j] / sqrt(dot3(d, d) + eps2);
        }
    }
    for (int i = 0; i < n; i++) {
        int k = ib[i], nm;
        vec3 xm[KFMAX], um[KFMAX];
        if (k < 0) continue;
        grpcomp(k, rp[i], vp[i], tnow, &nm, xm, um);
        for (int l = 0; l < nm; l++) {
            double u[3] = {um[l][0] - vp[i][0], um[l][1] - vp[i][1], um[l][2] - vp[i][2]};
            ke += 0.5 * bmass[k][l] * dot3(u, u);
            for (int l2 = l + 1; l2 < nm; l2++) {
                double d[3] = {xm[l][0] - xm[l2][0], xm[l][1] - xm[l2][1], xm[l][2] - xm[l2][2]};
                pe -= bmass[k][l] * bmass[k][l2] / norm3(d);
            }
        }
    }
    return ke + pe + emerged + arpotcorr(tnow, 1);
}
static double dbge1;
static void dbgbefore(double tnow)
{
    if (dbgon < 0) dbgon = getenv("HENDOS_DEBUG") != NULL;
    if (dbgon) dbge1 = dbgenergy(tnow);
}
static void dbgafter(const char *what, double tnow)
{
    if (!dbgon) return;
    double e2 = dbgenergy(tnow);
    printf("  DEBUG t = %.6f Myr %-10s dE/|E| = %10.3E  (groups %d)\n", tnow * tunit, what,
           (e2 - dbge1) / fabs(dbge1), nbin);
}

/* internal energy (kinetic about the c.o.m. + mutual potential) of nm
 * stars of masses mm at positions xm, velocities um */
double grpeint(int nm, const double *mm, vec3 xm[KFMAX], vec3 um[KFMAX])
{
    double mt = 0.0, vcm[3] = {0, 0, 0}, e = 0.0;
    for (int l = 0; l < nm; l++) {
        mt += mm[l];
        for (int c = 0; c < 3; c++) vcm[c] += mm[l] * um[l][c];
    }
    for (int c = 0; c < 3; c++) vcm[c] /= mt;
    for (int l = 0; l < nm; l++) {
        double u[3] = {um[l][0] - vcm[0], um[l][1] - vcm[1], um[l][2] - vcm[2]};
        e += 0.5 * mm[l] * dot3(u, u);
        for (int l2 = l + 1; l2 < nm; l2++) {
            double d[3] = {xm[l][0] - xm[l2][0], xm[l][1] - xm[l2][1], xm[l][2] - xm[l2][2]};
            e -= mm[l] * mm[l2] / norm3(d);
        }
    }
    return e;
}

/* advance a bound two-body relative orbit (G(m1+m2) = mu) by dt with
 * the f and g functions and the eccentric anomaly difference (Danby
 * 1988). Exact; dt is reduced modulo the period. */
void kepler(double mu, const double r0[3], const double v0[3], double dt, double r1[3],
            double v1[3])
{
    const double twopi = 6.283185307179586;
    double r0n = norm3(r0);
    double alpha = 2.0 / r0n - dot3(v0, v0) / mu;          /* 1/a */
    double a = 1.0 / alpha, nmo = sqrt(mu * alpha * alpha * alpha);
    double sig0 = dot3(r0, v0) / sqrt(mu);
    double c1 = 1.0 - r0n / a, c2 = sig0 / sqrt(a);
    double dm = fmod(nmo * dt, twopi);
    if (dm < 0.0) dm += twopi;
    double dtr = dm / nmo;
    /* Kepler's equation in difference form, dm = x - c1 sin x + c2 (1 - cos x).
     * x - dm = e (sin(x + E0) - sin E0), so the root is within 2e < 2 of dm;
     * Newton's method kept inside the bracket */
    double lo = dm - 2.0, hi = dm + 2.0, x = dm;
    for (int it = 0; it < 100; it++) {
        double ff = x - c1 * sin(x) + c2 * (1.0 - cos(x)) - dm;
        double fp = 1.0 - c1 * cos(x) + c2 * sin(x);
        if (ff > 0.0) hi = x; else lo = x;
        double xn = x - ff / fp;
        if (xn <= lo || xn >= hi) xn = 0.5 * (lo + hi);
        if (fabs(xn - x) < 1.e-15) { x = xn; break; }
        x = xn;
    }
    double ce = cos(x), se = sin(x);
    double f = 1.0 - a / r0n * (1.0 - ce), g = dtr + (se - x) / nmo;
    for (int k = 0; k < 3; k++) r1[k] = f * r0[k] + g * v0[k];
    double r1n = norm3(r1);
    double fd = -sqrt(mu * a) * se / (r1n * r0n), gd = 1.0 - a / r1n * (1.0 - ce);
    for (int k = 0; k < 3; k++) v1[k] = fd * r0[k] + gd * v0[k];
}

/* specific orbital energy and apocentre of a relative orbit (the
 * apocentre is huge if unbound) */
static void pairorbit(double mb, const double rr[3], const double vr[3], double *ebind,
                      double *rapo)
{
    double rn = norm3(rr), h[3];
    *ebind = 0.5 * dot3(vr, vr) - mb / rn;
    if (*ebind >= 0.0) {
        *rapo = HUGE_VAL;
        return;
    }
    double a = -0.5 * mb / *ebind;
    h[0] = rr[1] * vr[2] - rr[2] * vr[1];
    h[1] = rr[2] * vr[0] - rr[0] * vr[2];
    h[2] = rr[0] * vr[1] - rr[1] * vr[0];
    double ecc2 = fmax(0.0, 1.0 - dot3(h, h) / (mb * a));
    *rapo = a * (1.0 + sqrt(ecc2));
}

/* two-body pericentre distance of a relative orbit, q = h^2/(mb (1+e)) */
static double pericentre(double mb, const double rr[3], const double vr[3])
{
    double h[3];
    h[0] = rr[1] * vr[2] - rr[2] * vr[1];
    h[1] = rr[2] * vr[0] - rr[0] * vr[2];
    h[2] = rr[0] * vr[1] - rr[1] * vr[0];
    double h2 = dot3(h, h), ebind = 0.5 * dot3(vr, vr) - mb / norm3(rr);
    double ecc = sqrt(fmax(0.0, 1.0 + 2.0 * ebind * h2 / (mb * mb)));
    return h2 / (mb * (1.0 + ecc));
}

/* type of the star made by merging types t1 and t2: the more compact
 * wins. >>> CHANGE HERE for different rules */
int groupmergetype(int t1, int t2)
{
    if (t1 == TBHOLE || t2 == TBHOLE) return TBHOLE;
    if (t1 == TNSTAR || t2 == TNSTAR) return TNSTAR;
    if (t1 == TWDWARF || t2 == TWDWARF) return TWDWARF;
    if (t1 == TSTAR || t2 == TSTAR) return TSTAR;
    return t1 > t2 ? t1 : t2;
}

/* announce a merger on the screen and add a line to <runid>_mergers.dat */
void groupmergelog(double tnow, int id1, double m1, int t1, int id2, double m2, int t2,
                   int idn, double mn, int tn, const double xn[3], const double vn[3],
                   double rsep)
{
    char fname[2 * NSTR];
    double rsun = 6.957e8 / 3.0856775814913673e16;
    printf("  ** MERGER t = %.6f Myr: stars %d and %d (%.4f and %.4f Msun) are now star %d, "
           "%.4f Msun, type %s\n", tnow * tunit, id1, id2, m1, m2, idn, mn, typename_[tn]);
    snprintf(fname, sizeof fname, "%s%s_mergers.dat", datadir, runid);
    FILE *f = fopen(fname, "r");
    int exists = (f != NULL);
    if (f) fclose(f);
    f = fopen(fname, "a");
    if (!f) return;
    if (!exists) {
        fprintf(f, "# format                        = NBODY-MERGERS 1\n");
        fprintf(f, "# runid                         = %s\n", runid);
        fprintf(f, "# units                         = mass Msun, length pc, velocity km/s, "
                   "time Myr\n");
        fprintf(f, "# columns                       = time id1 m1 type1 id2 m2 type2 id m "
                   "type x y z vx vy vz separation\n");
        fprintf(f, "# column_units                  = Myr - Msun - - Msun - - Msun - pc pc pc "
                   "km/s km/s km/s Rsun\n");
        fprintf(f, "# column_notes                  = stars id1 and id2 merged to make star "
                   "id (the id of the more massive); x..vz are the new star\n");
        fprintf(f, "# end_header\n");
    }
    fprintf(f, "%24.16E%10d%24.16E %-12s%10d%24.16E %-12s%10d%24.16E %-12s", tnow * tunit, id1,
            m1, typename_[t1], id2, m2, typename_[t2], idn, mn, typename_[tn]);
    for (int k = 0; k < 3; k++) fprintf(f, "%24.16E", xn[k]);
    for (int k = 0; k < 3; k++) fprintf(f, "%24.16E", vn[k] * vunit);
    fprintf(f, "%16.8E\n", rsep / rsun);
    fclose(f);
}

/* positions xm and velocities um of the nm stars of group k at internal
 * time tnow, given its c.o.m. position xc and velocity vcm. A
 * regularised group whose c.o.m. has just reached the end of its
 * integration (tnow = bend) is first brought up to date (arshift: the
 * ghosts of merged stars are removed). */
void grpcomp(int k, const double xc[3], const double vcm[3], double tnow, int *nm,
             vec3 xm[KFMAX], vec3 um[KFMAX])
{
    double rr[3], vr[3], ro[3], vo[3];
    memset(xm, 0, sizeof(vec3) * KFMAX);
    memset(um, 0, sizeof(vec3) * KFMAX);
    if (bmode[k] == 1) {
        double mb = bmass[k][0] + bmass[k][1];
        kepler(mb, brel[k], brel[k] + 3, tnow - bepoch[k], rr, vr);
        for (int c = 0; c < 3; c++) {
            xm[0][c] = xc[c] + bmass[k][1] / mb * rr[c];
            um[0][c] = vcm[c] + bmass[k][1] / mb * vr[c];
            xm[1][c] = xc[c] - bmass[k][0] / mb * rr[c];
            um[1][c] = vcm[c] - bmass[k][0] / mb * vr[c];
        }
    } else if (bmode[k] == 3) {
        /* inner orbit (stars 0, 1) and outer orbit (star 2 about their c.o.m.) */
        double mb = bmass[k][0] + bmass[k][1], mt = mb + bmass[k][2];
        kepler(mb, brel[k], brel[k] + 3, tnow - bepoch[k], rr, vr);
        kepler(mt, brel2[k], brel2[k] + 3, tnow - bepoch[k], ro, vo);
        for (int c = 0; c < 3; c++) {
            double xi = xc[c] - bmass[k][2] / mt * ro[c], vi = vcm[c] - bmass[k][2] / mt * vo[c];
            xm[0][c] = xi + bmass[k][1] / mb * rr[c];
            um[0][c] = vi + bmass[k][1] / mb * vr[c];
            xm[1][c] = xi - bmass[k][0] / mb * rr[c];
            um[1][c] = vi - bmass[k][0] / mb * vr[c];
            xm[2][c] = xc[c] + mb / mt * ro[c];
            um[2][c] = vcm[c] + mb / mt * vo[c];
        }
    } else {
        double mm[KFMAX];
        vec3 xr[KFMAX], ur[KFMAX];
        if (fabs(tnow - bend[k]) <= 0.5 * tick) arshift(k);
        int nr = arstate(k, tnow, mm, xr, ur);
        for (int l = 0; l < nr; l++)
            for (int c = 0; c < 3; c++) {
                xm[l][c] = xc[c] + xr[l][c];
                um[l][c] = vcm[c] + ur[l][c];
            }
        *nm = nr;
        return;
    }
    *nm = bnum[k];
}

/* two stars have come within rcoll and mergers are off: stop the run */
static void collision(int i, int j, double rsep, double q, double tnow)
{
    double rsun = 6.957e8 / 3.0856775814913673e16;
    int two[2] = {i, j};
    printf("\n******************************************************************************\n");
    printf("  STELLAR COLLISION at t = %.6f Myr\n", tnow * tunit);
    printf("******************************************************************************\n");
    for (int c = 0; c < 2; c++) {
        int l = two[c];
        if (ib[l] < 0)
            printf("  star  id %12d   mass %12.5f Msun   type %s\n", id[l], m[l],
                   typename_[itype[l]]);
        else
            printf("  group ids %12d%12d ...   mass %12.5f Msun\n", bident[ib[l]][0],
                   bident[ib[l]][1], m[l]);
    }
    printf("  separation now          %12.4E Rsun\n", rsep / rsun);
    printf("  two-body pericentre     %12.4E Rsun\n", q / rsun);
    printf("  collision radius rcoll  %12.4E Rsun\n", rcoll);
    printf("  position                %14.6f%14.6f%14.6f pc\n", rp[i][0], rp[i][1], rp[i][2]);
    printf("  The run has been stopped: mergers are switched off.\n");
    printf("******************************************************************************\n");
    exit(1);
}

/* perturbation gamma on a (possible) group of mass mb and size gsize
 * at xc from all bodies except the nex in ex, and the squared shortest
 * pairwise free-fall/flyby time of its c.o.m. Uses rp, vp, and the
 * tree when it is on and up to date. */
static void gperturb(const double xc[3], const double vcm[3], double mb, double gsize, int nex,
                     const int *ex, double *gam, double *tp2)
{
    if (usetree && treevalid) {
        double sum3t;
        treeperturb(xc, vcm, mb, nex, ex, &sum3t, tp2);
        *gam = 2.0 * gsize * gsize * gsize / mb * sum3t;
        return;
    }
    double sum3 = 0.0, t2 = HUGE_VAL;
#pragma omp parallel for reduction(+:sum3) reduction(min:t2) schedule(static) if (n >= nompj)
    for (int j = 0; j < n; j++) {
        int skip = 0;
        for (int e = 0; e < nex; e++)
            if (ex[e] == j) skip = 1;
        if (skip) continue;
        double dx = rp[j][0] - xc[0], dy = rp[j][1] - xc[1], dz = rp[j][2] - xc[2];
        double d2 = dx * dx + dy * dy + dz * dz;
        double dv2 = (vp[j][0] - vcm[0]) * (vp[j][0] - vcm[0])
                     + (vp[j][1] - vcm[1]) * (vp[j][1] - vcm[1])
                     + (vp[j][2] - vcm[2]) * (vp[j][2] - vcm[2]);
        double d3 = d2 * sqrt(d2);
        sum3 += m[j] / d3;
        double tff = d3 / (mb + m[j]), tfl = d2 / (dv2 + 1e-300);
        if (tff < t2) t2 = tff;
        if (tfl < t2) t2 = tfl;
    }
    *gam = 2.0 * gsize * gsize * gsize / mb * sum3;
    *tp2 = t2;
}

/* the body nearest to body i, not counting body jex (-1: any body);
 * -1 if there is none */
static int nearestx(int i, int jex)
{
    int jn = -1;
    double gbest = HUGE_VAL;
#pragma omp parallel if (n >= nompj)
    {
        int jbest = -1;
        double dbest = HUGE_VAL;
#pragma omp for schedule(static) nowait
        for (int j = 0; j < n; j++) {
            if (j == i || j == jex) continue;
            double dx = rp[j][0] - rp[i][0], dy = rp[j][1] - rp[i][1], dz = rp[j][2] - rp[i][2];
            double d2 = dx * dx + dy * dy + dz * dz;
            if (d2 < dbest) {
                dbest = d2;
                jbest = j;
            }
        }
#pragma omp critical(nearestcrit)
        if (jbest >= 0 && (dbest < gbest || (dbest == gbest && jbest < jn))) {
            gbest = dbest;
            jn = jbest;
        }
    }
    return jn;
}

/* copy body src into slot dst */
static void moverow(int src, int dst)
{
#define CPV(a) memcpy(a[dst], a[src], sizeof(vec3))
    CPV(r); CPV(v); CPV(rp); CPV(vp); CPV(f0); CPV(f1); CPV(f2); CPV(f3); CPV(f2e); CPV(f3e);
#undef CPV
    m[dst] = m[src];
    id[dst] = id[src];
    itype[dst] = itype[src];
    ib[dst] = ib[src];
    lev[dst] = lev[src];
    tk[dst] = tk[src];
    dtk[dst] = dtk[src];
    tpair2[dst] = tpair2[src];
    rmin2[dst] = rmin2[src];
    phi[dst] = phi[src];
    if (ib[dst] >= 0) bcm[ib[dst]] = dst;
    premap(src, dst);
    if (useac) acmoverow(src, dst);
}

/* start off the nl new bodies in list at integer time t8: forces from
 * all other bodies (whose rp, vp are at t8), then the first level of
 * each, which is at least levmin and divides t8 */
static void startbody(int nl, const int *list, int64_t t8, int levmin)
{
    tcur8 = t8;
    aclast = 1;
    nact = nl;
    for (int ia = 0; ia < nl; ia++) ilist[ia] = list[ia];
    if (usetree) treebuild();
    forces();
    for (int ia = 0; ia < nl; ia++) {
        int l = list[ia];
        for (int k = 0; k < 3; k++) {
            f0[l][k] = f0p[l][k];
            f1[l][k] = f1p[l][k];
            f2[l][k] = f3[l][k] = f2e[l][k] = f3e[l][k] = 0.0;
        }
        tk[l] = t8;
        double amod = norm3(f0[l]), jmod = norm3(f1[l]);
        double dtw = etap * sqrt(tpair2[l]);
        if (amod > 0.0 && jmod > 0.0) dtw = fmin(dtw, etas * amod / jmod);
        int knew = getlevel(dtw);
        lev[l] = commlev(knew > levmin ? knew : levmin, t8);
        dtk[l] = (int64_t)1 << (kmax - lev[l]);
        if (useac) acstart(l);
    }
}

/* merge single stars i and j into one star at their centre of mass
 * (their predicted positions rp, vp are used, which for the active
 * star are its new position). Mass and momentum are conserved; the
 * star keeps the id of the more massive. The orbital energy of the
 * pair is added to emerged. */
static void mergestars(int i, int j, int64_t t8, double tnow)
{
    double mi = m[i], mj = m[j], mn = mi + mj, rr[3], vr[3], xc[3], vcm[3];
    for (int k = 0; k < 3; k++) {
        rr[k] = rp[i][k] - rp[j][k];
        vr[k] = vp[i][k] - vp[j][k];
        xc[k] = (mi * rp[i][k] + mj * rp[j][k]) / mn;
        vcm[k] = (mi * vp[i][k] + mj * vp[j][k]) / mn;
    }
    double rsep = norm3(rr);
    emerged += 0.5 * mi * mj / mn * dot3(vr, vr) - mi * mj / sqrt(rsep * rsep + eps2);
    int idn = mj > mi ? id[j] : id[i];
    int tn = groupmergetype(itype[i], itype[j]);
    groupmergelog(tnow, id[i], mi, itype[i], id[j], mj, itype[j], idn, mn, tn, xc, vcm, rsep);

    int lo = i < j ? i : j, hi = i < j ? j : i;
    for (int k = 0; k < 3; k++) {
        r[lo][k] = rp[lo][k] = xc[k];
        v[lo][k] = vp[lo][k] = vcm[k];
    }
    m[lo] = mn;
    id[lo] = idn;
    itype[lo] = tn;
    ib[lo] = -1;
    tk[lo] = t8;
    /* remove slot hi by moving the last body into it */
    premap(hi, lo);
    if (useac) {
        acremoved(hi, lo);
        acok[lo] = 0;
    }
    if (hi != n - 1) {
        moverow(n - 1, hi);
        if (useac) acmoved(n - 1, hi);
    }
    n--;
    ntot--;
    nmerged++;
    treevalid = 0;
    startbody(1, &lo, t8, 0);
}

/* after each block: for each active body that came within 2 rcoll of
 * another, find that neighbour and test the pair: a collision if they
 * are within rcoll, or approaching with a two-body pericentre inside
 * rcoll. The stars are then merged, or the run stopped if mergers are
 * off. If one of the two is a group it is opened instead. At most one
 * merger per block. */
void collcheck(int64_t t8, double tnow)
{
    int nact0 = nact, nm, list[KFMAX];
    int *alist = xalloc(sizeof(int) * (size_t)(nact0 > 0 ? nact0 : 1));
    memcpy(alist, ilist, sizeof(int) * (size_t)nact0);
    for (int ia = 0; ia < nact0; ia++) {
        int i = alist[ia];
        if (i >= n) continue;
        if (rmin2[i] >= 4.0 * rcolli * rcolli) continue;
        int j = nearestx(i, -1);
        if (j < 0) continue;
        double rr[3], vr[3];
        for (int k = 0; k < 3; k++) {
            rr[k] = rp[i][k] - rp[j][k];
            vr[k] = vp[i][k] - vp[j][k];
        }
        double rsep = norm3(rr), q = pericentre(m[i] + m[j], rr, vr);
        if (rsep < rcolli || (q < rcolli && dot3(rr, vr) < 0.0)) {
            if (!mergers) collision(i, j, rsep, q, tnow);
            if (ib[i] >= 0) {
                if (tk[i] == t8) opengroup(i, t8, tnow, &nm, list);
            } else if (ib[j] >= 0) {
                if (tk[j] == t8) opengroup(j, t8, tnow, &nm, list);
            } else {
                mergestars(i, j, t8, tnow);
            }
            break;
        }
    }
    free(alist);
}

/* if group k (few-body, three stars, up to date) is a stable hierarchy
 * whose inner binary does not touch, freeze it as two nested Kepler
 * orbits (mode 3), with its stars in the order inner, inner, outer */
static void tryfreeze(int k)
{
    int stable, i1, i2, i3;
    double rel[6], rel2[6], qin, qout;
    hierarchy(bmass[k], bx[k], bv[k], &stable, &i1, &i2, &i3, rel, rel2, &qin, &qout);
    if (!stable) return;
    if (rcoll > 0.0 && qin < rcolli) return;
    /* its true internal energy now, at which it is counted while frozen
     * (see diagnose in io.c) - before the stars are put in a new order */
    bein[k] = grpeint(3, bmass[k], bx[k], bv[k]);
    int idx[3] = {i1, i2, i3}, ti[3], tt[3];
    double tm[3];
    for (int c = 0; c < 3; c++) {
        tm[c] = bmass[k][idx[c]];
        ti[c] = bident[k][idx[c]];
        tt[c] = btype[k][idx[c]];
    }
    for (int c = 0; c < 3; c++) {
        bmass[k][c] = tm[c];
        bident[k][c] = ti[c];
        btype[k][c] = tt[c];
    }
    memcpy(brel[k], rel, sizeof rel);
    memcpy(brel2[k], rel2, sizeof rel2);
    bmode[k] = 3;
    arreset(k);
    nfrozen++;
}

/* replace the nm single stars in list (all at integer time t8 =
 * internal time tnow) by one c.o.m. body, as a group of kind mode
 * (1 Kepler binary, 2 regularised) */
static void closegroup(int nm, const int *list, int mode, int64_t t8, double tnow)
{
    int k = nbin++, lo = list[0];
    double mb = 0.0, xc[3] = {0, 0, 0}, vcm[3] = {0, 0, 0}, fc[4][3], gsize = 0.0, ebind;
    double gam, tp2;
    memset(fc, 0, sizeof fc);
    bnum[k] = nm;
    bmode[k] = mode;
    bepoch[k] = tnow;
    for (int l = 0; l < nm; l++) {
        int i = list[l];
        if (i < lo) lo = i;
        bmass[k][l] = m[i];
        bident[k][l] = id[i];
        btype[k][l] = itype[i];
        mb += m[i];
        for (int c = 0; c < 3; c++) {
            xc[c] += m[i] * r[i][c];
            vcm[c] += m[i] * v[i][c];
            /* c.o.m. derivatives: the mass-weighted mean, in which the
             * mutual forces of the stars cancel exactly */
            fc[0][c] += m[i] * f0[i][c];
            fc[1][c] += m[i] * f1[i][c];
            fc[2][c] += m[i] * f2e[i][c];
            fc[3][c] += m[i] * f3e[i][c];
        }
    }
    for (int c = 0; c < 3; c++) {
        xc[c] /= mb;
        vcm[c] /= mb;
        for (int q = 0; q < 4; q++) fc[q][c] /= mb;
    }
    if (mode == 1) {
        for (int c = 0; c < 3; c++) {
            brel[k][c] = r[list[0]][c] - r[list[1]][c];
            brel[k][c + 3] = v[list[0]][c] - v[list[1]][c];
        }
        pairorbit(mb, brel[k], brel[k] + 3, &ebind, &gsize);
        nclosed++;
    } else {
        for (int l = 0; l < nm; l++) {
            int i = list[l];
            double d[3];
            for (int c = 0; c < 3; c++) {
                bx[k][l][c] = d[c] = r[i][c] - xc[c];
                bv[k][l][c] = v[i][c] - vcm[c];
            }
            gsize = fmax(gsize, 2.0 * norm3(d));
        }
        nfewformed++;
    }
    /* timestep of the c.o.m. (before the stars are removed) */
    gperturb(xc, vcm, mb, gsize, nm, list, &gam, &tp2);
    for (int c = 0; c < 3; c++) {
        f0[lo][c] = fc[0][c];
        f1[lo][c] = fc[1][c];
        f2e[lo][c] = fc[2][c];
        f3e[lo][c] = fc[3][c];
        r[lo][c] = rp[lo][c] = xc[c];
        v[lo][c] = vp[lo][c] = vcm[c];
    }
    m[lo] = mb;
    id[lo] = 0;
    itype[lo] = 0;
    ib[lo] = k;
    bcm[k] = lo;
    npl[k] = 0;
    tk[lo] = t8;
    tpair2[lo] = tp2;
    double dtw = etap * sqrt(tp2), amod = norm3(f0[lo]), jmod = norm3(f1[lo]);
    if (amod > 0.0 && jmod > 0.0) dtw = fmin(dtw, etas * amod / jmod);
    lev[lo] = commlev(getlevel(dtw), t8);
    dtk[lo] = (int64_t)1 << (kmax - lev[lo]);
    /* remove the other slots, highest first, moving the last body into each */
    int rest[KFMAX], nrest = 0;
    for (int l = 0; l < nm; l++)
        if (list[l] != lo) rest[nrest++] = list[l];
    while (nrest > 0) {
        int lmax = 0;
        for (int l = 1; l < nrest; l++)
            if (rest[l] > rest[lmax]) lmax = l;
        int h = rest[lmax];
        rest[lmax] = rest[--nrest];
        premap(h, lo);
        if (useac) acremoved(h, lo);
        if (h != n - 1) {
            moverow(n - 1, h);
            if (useac) acmoved(n - 1, h);
        }
        n--;
    }
    treevalid = 0;
    /* with the neighbour scheme the c.o.m. body needs its own list and
     * force split, so it is started with a full force sum */
    if (useac) {
        acok[lo] = 0;
        startbody(1, &lo, t8, 0);
    }
    if (mode == 2) {
        bein[k] = 0.0;                 /* a pair: its separation when made (see bincheck) */
        if (nm == 2) {
            double d[3] = {bx[k][0][0] - bx[k][1][0], bx[k][0][1] - bx[k][1][1],
                           bx[k][0][2] - bx[k][1][2]};
            bein[k] = norm3(d);
        }
        arstart(k, t8);
    }
    ardirty();
}

/* split the c.o.m. body i (at integer time t8 = internal time tnow)
 * back into its nm stars, in slot i and slots n onwards, and start
 * them off. list returns the slots of the stars. */
static void opengroup(int i, int64_t t8, double tnow, int *nm, int *list)
{
    int k = ib[i], levc = lev[i], mode = bmode[k];
    bmode_last_open = mode;
    vec3 xm[KFMAX], um[KFMAX];
    double xc[3], vcm[3];
    memcpy(xc, r[i], sizeof xc);
    memcpy(vcm, v[i], sizeof vcm);
    grpcomp(k, xc, vcm, tnow, nm, xm, um);
    /* a frozen triple: the change of its true energy while it was frozen */
    if (mode == 3) efrozen += bein[k] - grpeint(*nm, bmass[k], xm, um);
    for (int l = 0; l < *nm; l++) {
        int j = i;
        if (l > 0) {
            j = n++;
            if (useac) {
                nnb[j] = 0;
                acok[j] = 0;
                acadded(j, i);
            }
        }
        for (int c = 0; c < 3; c++) {
            r[j][c] = rp[j][c] = xm[l][c];
            v[j][c] = vp[j][c] = um[l][c];
        }
        m[j] = bmass[k][l];
        id[j] = bident[k][l];
        itype[j] = btype[k][l];
        ib[j] = -1;
        phi[j] = phi[i];
        rmin2[j] = HUGE_VAL;
        list[l] = j;
    }
    if (useac) acok[i] = 0;
    arreset(k);
    /* remove group k from the table, moving the last group into its place */
    int last = nbin - 1;
    if (k != last) {
        bnum[k] = bnum[last];
        bmode[k] = bmode[last];
        memcpy(bmass[k], bmass[last], sizeof bmass[k]);
        memcpy(bident[k], bident[last], sizeof bident[k]);
        memcpy(btype[k], btype[last], sizeof btype[k]);
        memcpy(brel[k], brel[last], sizeof brel[k]);
        memcpy(brel2[k], brel2[last], sizeof brel2[k]);
        bepoch[k] = bepoch[last];
        bein[k] = bein[last];
        memcpy(bx[k], bx[last], sizeof bx[k]);
        memcpy(bv[k], bv[last], sizeof bv[k]);
        bcm[k] = bcm[last];
        ib[bcm[k]] = k;
        armovegroup(last, k);
    }
    nbin--;
    if (mode == 1) nopened++; else nfewopened++;
    treevalid = 0;
    /* forces on the stars and their first timesteps (at least the c.o.m.
     * level, so the step divides t8) */
    startbody(*nm, list, t8, levc);
}

/* what single stars i and j are to each other:
 *   0  nothing to do
 *   1  a bound pair that is unperturbed round its whole orbit: gamma with
 *      its apocentre as the size < gclose -> a Kepler binary
 *   2  a pair that is isolated now: gamma with its CURRENT SEPARATION as
 *      the size < gfew, and either bound or coming together -> a
 *      regularised pair (only with fewbody on)
 * Judging a regularised pair by its current separation is what KS
 * regularisation does (Aarseth 2003): a very eccentric binary whose
 * apocentre reaches out among other stars is still regularised while it
 * plunges through pericentre, and let go when it swings out again; an
 * unbound pair is regularised for a close hyperbolic encounter (only a
 * deep plunge, pericentre < PLUNGE = 1/1000 of its separation). Either,
 * integrated directly, can need steps shorter than the shortest level
 * (kmax) and ruin the energy; regularised, any pericentre is easy. With
 * the perturbers treated properly (ar.c) it does not matter that the
 * pair is perturbed at its apocentre.
 * A bound pair whose pericentre is inside rcoll is never closed as a
 * Kepler binary, which would miss the collision: it is made a
 * regularised pair, which watches for it, or with fewbody off it is left
 * in the main integration, and collcheck merges the stars (or stops the
 * run) when they actually meet - not before. */
/* would a group of mass mb, size sz at xc have more than NPMAX bodies
 * with gamma > 100 gpert? Then its list fills and it is let go at once
 * (the rule in bincheck), so it is not worth making a group at all */
static int overfull(const double xc[3], double mb, double sz, int i, int j)
{
    int cnt = 0;
    double f = 2.0 * sz * sz * sz / (mb * 100.0 * gpert);
    for (int l = 0; l < n; l++) {
        if (l == i || l == j || m[l] <= 0.0) continue;
        double d[3] = {r[l][0] - xc[0], r[l][1] - xc[1], r[l][2] - xc[2]};
        double d2 = dot3(d, d);
        if (f * m[l] > d2 * sqrt(d2) && ++cnt >= NPMAX) return 1;
    }
    return 0;
}

static int pairstate(int i, int j)
{
    if (ib[i] >= 0 || ib[j] >= 0) return 0;
    double mb = m[i] + m[j], rr[3], vr[3], xc[3], vcm[3], ebind, rapo, gam, tp2;
    for (int k = 0; k < 3; k++) {
        rr[k] = r[i][k] - r[j][k];
        vr[k] = v[i][k] - v[j][k];
        xc[k] = (m[i] * r[i][k] + m[j] * r[j][k]) / mb;
        vcm[k] = (m[i] * v[i][k] + m[j] * v[j][k]) / mb;
    }
    pairorbit(mb, rr, vr, &ebind, &rapo);
    double rn = norm3(rr);
    int two[2] = {i, j}, bound = ebind < 0.0;
    if (!bound && (!fewbody || dot3(rr, vr) >= 0.0)) return 0;
    /* an unbound pair only for a deep plunge: its pericentre far inside its
     * separation now (PLUNGE). An ordinary close passage the main
     * integrator follows easily (to 1e-9 in a Plummer sphere), and
     * regularised it only leaves a little error behind each time it is
     * made a group and let go (found: with PLUNGE = 0.1, a 1000-star
     * Plummer sphere made 57 such groups in 10 Myr and the energy error
     * went from 1e-8 to 2e-7; the plunges that wreck the energy, in cold
     * fractals, have pericentres a millionth of the separation) */
    if (!bound && pericentre(mb, rr, vr) > PLUNGE * rn) return 0;
    /* gamma scales as size^3, so one sum gives it for both sizes */
    gperturb(xc, vcm, mb, rn, 2, two, &gam, &tp2);
    if (bound) {
        double gapo = gam * pow(rapo / rn, 3.0);
        int collide = rcoll > 0.0 && pericentre(mb, rr, vr) < rcolli;
        if (gapo < gclose && !collide) return 1;
    }
    if (!fewbody || gam >= gfew) return 0;
    /* a wide bound pair whose list would fill at apocentre: it is let go
     * as soon as it swings out, and made a group again at the next
     * pericentre (found: a 0.03 pc binary in a 1000-star Plummer sphere,
     * made a group and let go every orbit, 5400 times in 10 Myr, put the
     * energy error up from 1e-8 to 2e-7). The main integrator follows it
     * well, so leave it there - unless it plunges deeply */
    if (bound && pericentre(mb, rr, vr) > PLUNGE * rn && overfull(xc, mb, rapo, i, j)) return 0;
    return 2;
}

/* perturbation of the group whose c.o.m. body is i (a regularised group
 * is brought up to date first): gcur with its current size (separation
 * of a pair, diameter of a larger group), gorb with the size of its
 * orbit (apocentre of a bound pair; otherwise the same as gcur) */
static void grpgamma(int i, double tnow, double *gcur, double *gorb, double *gs)
{
    int k = ib[i], nm;
    vec3 xm[KFMAX], um[KFMAX];
    double gsize = 0.0, ebind, rapo = 0.0, tp2, xc[3], vcm[3];
    memcpy(xc, r[i], sizeof xc);
    memcpy(vcm, v[i], sizeof vcm);
    grpcomp(k, xc, vcm, tnow, &nm, xm, um);
    if (nm == 2) {
        double rr[3], vr[3];
        for (int c = 0; c < 3; c++) {
            rr[c] = xm[0][c] - xm[1][c];
            vr[c] = um[0][c] - um[1][c];
        }
        pairorbit(m[i], rr, vr, &ebind, &rapo);
        gsize = norm3(rr);
    } else {
        for (int l = 0; l < nm; l++) {
            double d[3];
            for (int c = 0; c < 3; c++) d[c] = xm[l][c] - xc[c];
            gsize = fmax(gsize, 2.0 * norm3(d));
        }
    }
    gperturb(xc, vcm, m[i], gsize, 1, &i, gcur, &tp2);
    *gs = gsize;
    *gorb = *gcur;
    if (nm == 2) *gorb = rapo == HUGE_VAL ? HUGE_VAL : *gcur * pow(rapo / gsize, 3.0);
}

/* group k, c.o.m. in slot i, at internal time tnow: a Kepler binary or a
 * frozen triple becomes a regularised group (no slots change) */
static void toar(int i, int64_t t8, double tnow)
{
    int k = ib[i], nm;
    vec3 xm[KFMAX], um[KFMAX];
    grpcomp(k, r[i], v[i], tnow, &nm, xm, um);
    if (bmode[k] == 3) efrozen += bein[k] - grpeint(nm, bmass[k], xm, um);
    for (int l = 0; l < nm; l++)
        for (int c = 0; c < 3; c++) {
            bx[k][l][c] = xm[l][c] - r[i][c];
            bv[k][l][c] = um[l][c] - v[i][c];
        }
    bmode[k] = 2;
    bein[k] = 0.0;
    if (nm == 2) {
        double d[3] = {bx[k][0][0] - bx[k][1][0], bx[k][0][1] - bx[k][1][1],
                       bx[k][0][2] - bx[k][1][2]};
        bein[k] = norm3(d);
    }
    nfewformed++;
    arstart(k, t8);
}

/* a regularised pair (c.o.m. in slot i, up to date at tnow) has become
 * unperturbed: it becomes a Kepler binary */
static void tokepler(int i, double tnow)
{
    int k = ib[i];
    for (int c = 0; c < 3; c++) {
        brel[k][c] = bx[k][0][c] - bx[k][1][c];
        brel[k][c + 3] = bv[k][0][c] - bv[k][1][c];
    }
    bepoch[k] = tnow;
    bmode[k] = 1;
    arreset(k);
    nclosed++;
}

/* try to make a regularised group out of the nb bodies in bod (single
 * stars or groups, all at integer time t8). It is made if it has 3 to
 * KFMAX stars, its perturbation by the rest of the cluster is below
 * gfew, and it would pay: the main integration must have at least
 * NFEWFAC*nm*(nm-1) bodies (see below for why), unless it is an
 * unperturbed stable triple, which costs nothing once frozen. Returns 1
 * if the group was made. (Pairs are made in bincheck and binstart.)
 *
 * WHEN A GROUP PAYS. Taking nm stars out saves work only if the main
 * integration is big enough. There, each step of a close pair costs
 * about 2.5 n pair interactions (two force sums over the n bodies and
 * the prediction of all of them). The rule n >= NFEWFAC nm (nm-1),
 * NFEWFAC = 2,
 *      stars in the group        3     4     5     6
 *      bodies needed            12    24    40    60
 * was measured for the old Hermite few-body integrator, which took 2 to 4
 * times as many steps per orbit as the main integrator. The regularised
 * integrator takes 10 to 20 steps per orbit whatever the eccentricity,
 * each of about 40 force evaluations of nm(nm-1)/2 pairs plus nm times
 * the number of perturbers, so the rule is still about right (if
 * anything a little cautious). Accuracy does not depend on n at all: it
 * is set by gfew and artol. An unperturbed stable triple is always worth
 * it, whatever n is, because frozen as Kepler orbits it costs nothing.
 *
 * >>> CHANGE HERE (NFEWFAC, at the top of this file) once the break-even has
 *     been measured with the compiled code: benchmark/bench_kira.py
 *     --skip-kira --modes nofew direct
 */
static int tryfew(int nb, const int *bod, int64_t t8, double tnow)
{
    int nm = 0;
    vec3 xm[KFMAX], um[KFMAX], xg[KFMAX], ug[KFMAX];
    double mm[KFMAX], mt = 0.0, xc[3] = {0, 0, 0}, vcm[3] = {0, 0, 0}, gsize = 0.0, gam, tp2;
    for (int b = 0; b < nb; b++) {
        int i = bod[b], k = ib[i], nk = (k >= 0) ? bnum[k] : 1;
        if (nm + nk > KFMAX) return 0;
        if (k < 0) {
            memcpy(xm[nm], r[i], sizeof(vec3));
            memcpy(um[nm], v[i], sizeof(vec3));
            mm[nm] = m[i];
        } else {
            grpcomp(k, r[i], v[i], tnow, &nk, xg, ug);
            if (nm + nk > KFMAX) return 0;
            for (int l = 0; l < nk; l++) {
                memcpy(xm[nm + l], xg[l], sizeof(vec3));
                memcpy(um[nm + l], ug[l], sizeof(vec3));
                mm[nm + l] = bmass[k][l];
            }
        }
        nm += nk;
        mt += m[i];
        for (int c = 0; c < 3; c++) {
            xc[c] += m[i] * r[i][c];
            vcm[c] += m[i] * v[i][c];
        }
    }
    if (nm < 3) return 0;
    for (int c = 0; c < 3; c++) {
        xc[c] /= mt;
        vcm[c] /= mt;
    }
    for (int l = 0; l < nm; l++) {
        double d[3];
        for (int c = 0; c < 3; c++) d[c] = xm[l][c] - xc[c];
        gsize = fmax(gsize, 2.0 * norm3(d));
    }
    gperturb(xc, vcm, mt, gsize, nb, bod, &gam, &tp2);
    if (gam >= gfew) return 0;
    /* would it pay? */
    if (n < NFEWFAC * nm * (nm - 1)) {
        int stable, i1, i2, i3;
        double rel[6], rel2[6], qin, qout;
        if (nm != 3 || gam >= gclose) return 0;
        hierarchy(mm, xm, um, &stable, &i1, &i2, &i3, rel, rel2, &qin, &qout);
        if (!stable) return 0;
        if (rcoll > 0.0 && qin < rcolli) return 0;
    }
    /* make it: open the groups, then close all the stars together */
    int nl = 0, list[KFMAX], lk[KFMAX], nk;
    for (int b = 0; b < nb; b++) {
        int i = bod[b];
        if (ib[i] >= 0) {
            opengroup(i, t8, tnow, &nk, lk);
            for (int l = 0; l < nk; l++) list[nl++] = lk[l];
        } else {
            list[nl++] = i;
        }
    }
    if (nl < 3) return 0;
    closegroup(nl, list, 2, t8, tnow);
    /* an unperturbed stable triple is frozen at once (the new group is the
     * last one) */
    if (nl == 3 && gam < gclose) tryfreeze(nbin - 1);
    return 1;
}

/* after each block (all active bodies at integer time t8 = internal
 * time tnow), among the bodies that have just stepped:
 *  - groups: change the kind of one that has become more or less
 *    perturbed (Kepler <-> regularised, frozen <-> regularised); open one
 *    that has become too perturbed (or take its perturber into a
 *    regularised group, if that would be isolated)
 *  - single stars with short steps: close one bound pair (Kepler or
 *    regularised), or make a regularised group of a perturbed pair and
 *    the next nearest body
 * (see the top of this file) */
void bincheck(int64_t t8, double tnow)
{
    int nact0 = nact, nm, list[KFMAX];
    int *alist = xalloc(sizeof(int) * (size_t)(nact0 > 0 ? nact0 : 1));
    memcpy(alist, ilist, sizeof(int) * (size_t)nact0);

    for (int ia = 0; ia < nact0; ia++) {
        int i = alist[ia];
        if (i >= n || ib[i] < 0 || tk[i] != t8) continue;
        double gam, gorb, gs;
        grpgamma(i, tnow, &gam, &gorb, &gs);
        int k = ib[i], mode = bmode[k];
        if (mode == 2 && bnum[k] < 2) {          /* mergers have left a single star */
            dbgbefore(tnow);
            opengroup(i, t8, tnow, &nm, list);
            dbgafter("open<2", tnow);
            break;
        }
        /* a regularised group that should simply be let go (no stars move
         * slots, so it counts as the one open of this block):
         *  - an unbound pair moving apart: the encounter is over, and it is
         *    only growing - opened once it is twice as far apart as when it
         *    was made a group (kept in bein, which only frozen triples
         *    otherwise use), or at gamma > gfew; a receding unbound pair is
         *    never made a group, so it cannot flicker
         *  - a group whose perturber list is full while the weakest perturber
         *    kept is still strong: it has grown too big for its perturbers to
         *    be counted, and the ones left off are missing from the energy
         *    (found: unbound pairs that had drifted to 0.05 pc in a 200-star
         *    cluster each put 2e-7 into the energy when finally opened) */
        if (fewbody && mode == 2) {
            int letgo = 0;
            if (bnum[k] == 2 && gorb == HUGE_VAL) {
                double rr[3], vr[3];
                for (int c = 0; c < 3; c++) {
                    rr[c] = bx[k][0][c] - bx[k][1][c];
                    vr[c] = bv[k][0][c] - bv[k][1][c];
                }
                if (dot3(rr, vr) > 0.0 && (gam > gfew || norm3(rr) > 2.0 * bein[k])) letgo = 1;
            }
            if (npl[k] == NPMAX) {
                int j = plist[k][NPMAX - 1];
                double d[3] = {r[j][0] - r[i][0], r[j][1] - r[i][1], r[j][2] - r[i][2]};
                double dn = norm3(d), gj = 2.0 * gs * gs * gs * m[j] / (m[i] * dn * dn * dn);
                if (gj > 100.0 * gpert) letgo = 1;
            }
            if (letgo) {
                dbgbefore(tnow);
                opengroup(i, t8, tnow, &nm, list);
                dbgafter("let go", tnow);
                break;
            }
        }
        if (gam > 10.0 * gfew || (!fewbody && gorb > gopen)) {
            int done = 0;
            if (fewbody) {
                int j = nearestx(i, -1);
                if (j >= 0 && tk[j] == t8) {
                    int two[2] = {i, j};
                    dbgbefore(tnow);
                    done = tryfew(2, two, t8, tnow);
                    if (done) dbgafter("few(g+1)", tnow);
                }
            }
            if (!done) {
                dbgbefore(tnow);
                opengroup(i, t8, tnow, &nm, list);
                dbgafter(bmode_last_open == 1 ? "open-bin" : "open-grp", tnow);
            }
            break;
        }
        /* two groups meeting: each sees the other as a point mass, so the
         * interaction of their internal motions (quadrupole with quadrupole)
         * is left out. So a group whose nearest body is another group takes
         * it in as soon as that group alone perturbs it by more than gfew,
         * if the two together are isolated - not only when its whole
         * gamma passes 10 gfew - and if they are not, it is opened. (Found
         * in cold fractals: two regularised binaries that met as point
         * masses until gamma = 0.1 put errors of 2e-4 to 5e-4 into the
         * energy.) Its stars cannot be made a group again while the other
         * group is that close, because its gamma alone is then above gfew. */
        if (fewbody && mode == 2) {
            int j = nearestx(i, -1);
            if (j >= 0 && ib[j] >= 0 && tk[j] == t8) {
                double d[3] = {r[j][0] - r[i][0], r[j][1] - r[i][1], r[j][2] - r[i][2]};
                double dn = norm3(d), gj = 2.0 * gs * gs * gs * m[j] / (m[i] * dn * dn * dn);
                if (gj > gfew) {
                    int two[2] = {i, j};
                    dbgbefore(tnow);
                    if (tryfew(2, two, t8, tnow)) {
                        dbgafter("few(g+g)", tnow);
                        break;
                    }
                    /* they cannot be one group: this one is opened, so that
                     * the other sees its stars, exactly, as single stars */
                    dbgbefore(tnow);
                    opengroup(i, t8, tnow, &nm, list);
                    dbgafter("open-gg", tnow);
                    break;
                }
            }
        }
        /* changes of kind: no slots change, so several may be done */
        if ((mode == 1 && gorb > gopen) || (mode == 3 && gam > gopen)) {
            dbgbefore(tnow);
            toar(i, t8, tnow);
            dbgafter(mode == 1 ? "kep->AR" : "unfreeze", tnow);
        } else if (mode == 2 && gorb < gclose) {
            double vr[3] = {bv[k][0][0] - bv[k][1][0], bv[k][0][1] - bv[k][1][1],
                            bv[k][0][2] - bv[k][1][2]};
            double rr[3] = {bx[k][0][0] - bx[k][1][0], bx[k][0][1] - bx[k][1][1],
                            bx[k][0][2] - bx[k][1][2]};
            if (bnum[k] == 2 && 0.5 * dot3(vr, vr) < m[i] / norm3(rr)) {
                dbgbefore(tnow);
                tokepler(i, tnow);
                dbgafter("AR->kep", tnow);
            } else if (bnum[k] == 3) {
                dbgbefore(tnow);
                int64_t nf0 = nfrozen;
                tryfreeze(k);
                if (nfrozen > nf0) dbgafter("freeze", tnow);
            }
        }
    }

    /* single stars - only stars with short steps are candidates */
    for (int ia = 0; ia < nact0; ia++) {
        int i = alist[ia], done = 0;
        if (i >= n || ib[i] >= 0 || tk[i] != t8) continue;
        if ((double)dtk[i] * tick >= dtclosei) continue;
        int j = nearestx(i, -1);
        if (j < 0 || tk[j] != t8) continue;
        int two[2] = {i, j};
        if (ib[j] >= 0) {
            if (fewbody) {                               /* the nearest body is a group */
                dbgbefore(tnow);
                done = tryfew(2, two, t8, tnow);
                if (done) dbgafter("few(1+g)", tnow);
            }
        } else {
            int st = pairstate(i, j);
            if (st > 0) {
                dbgbefore(tnow);
                closegroup(2, two, st, t8, tnow);
                dbgafter(st == 1 ? "close" : "close-AR", tnow);
                done = 1;
            } else if (fewbody) {
                int l = nearestx(i, j);
                if (l >= 0 && tk[l] == t8) {
                    int three[3] = {i, j, l};
                    dbgbefore(tnow);
                    done = tryfew(3, three, t8, tnow);
                    if (done) dbgafter("few(3)", tnow);
                }
            }
        }
        if (done) break;
    }
    free(alist);
}

/* at the start (everyone at integer time 0 = internal time tnow):
 * close all bound pairs in the initial conditions (Kepler or
 * regularised), and make regularised groups of perturbed bound pairs
 * with their nearest body where those are isolated enough. After a
 * change slot i may hold a different body, so i only advances when
 * nothing was done. */
void binstart(double tnow)
{
    int i = 0;
    while (i < n) {
        int done = 0;
        if (ib[i] < 0) {
            int j = nearestx(i, -1);
            if (j >= 0 && ib[j] < 0) {
                int two[2] = {i, j}, st = pairstate(i, j);
                if (st > 0) {
                    closegroup(2, two, st, 0, tnow);
                    done = 1;
                } else if (fewbody) {
                    /* only bound pairs: a perturbed binary and its perturber */
                    double d[3], dv[3];
                    for (int k = 0; k < 3; k++) {
                        d[k] = r[i][k] - r[j][k];
                        dv[k] = v[i][k] - v[j][k];
                    }
                    if (0.5 * dot3(dv, dv) < (m[i] + m[j]) / norm3(d)) {
                        int l = nearestx(i, j);
                        if (l >= 0) {
                            int three[3] = {i, j, l};
                            done = tryfew(3, three, 0, tnow);
                        }
                    }
                }
            }
        }
        if (!done) i++;
    }
}
