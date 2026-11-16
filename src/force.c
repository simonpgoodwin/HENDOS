/* ==============================================================
 * HENDOS - HErmite N-body Dynamics Of Systems - force.c
 * GRAVITY AND TIMESTEPS: accelerations and jerks on the active bodies
 * by direct summation (with the Galactic tide), energies, and the
 * choice of each body's block timestep level.
 * (Archived Fortran equivalents: nbforce.f90 and nbstep.f90.)
 * ============================================================== */
#include <omp.h>
#include "hendos.h"

/* add the Galactic tide to the total acceleration a and jerk jk of
 * body i: a = T.x + 2 omega (vy,-vx,0); jerk = T.v + 2 omega (ay,-ax,0)
 * with the TOTAL acceleration, so a must be complete first */
void tideadd(int i, double a[3], double jk[3])
{
    if (!ltide) return;
    for (int k = 0; k < 3; k++) a[k] += dot3(ttid[k], rp[i]);
    a[0] += 2.0 * omega * vp[i][1];
    a[1] -= 2.0 * omega * vp[i][0];
    for (int k = 0; k < 3; k++) jk[k] += dot3(ttid[k], vp[i]);
    jk[0] += 2.0 * omega * a[1];
    jk[1] -= 2.0 * omega * a[0];
}

/* acceleration and jerk on body i from bodies j1..j2-1, with the
 * shortest pairwise free-fall / flyby time squared and the nearest
 * distance squared. The pair interaction of the whole code. */
static inline void sumrange(int i, int j1, int j2, double a[3], double jk[3], double *tmin2,
                            double *rm2)
{
    double xi = rp[i][0], yi = rp[i][1], zi = rp[i][2];
    double vxi = vp[i][0], vyi = vp[i][1], vzi = vp[i][2], mi = m[i];
    double ax = 0, ay = 0, az = 0, jx = 0, jy = 0, jz = 0, t2 = HUGE_VAL, r2 = HUGE_VAL;
    for (int j = j1; j < j2; j++) {
        if (j == i) continue;
        /* separations are j - i, so the force on i is +m(j)*dr/r^3 */
        double dx = rp[j][0] - xi, dy = rp[j][1] - yi, dz = rp[j][2] - zi;
        double dvx = vp[j][0] - vxi, dvy = vp[j][1] - vyi, dvz = vp[j][2] - vzi;
        double d2u = dx * dx + dy * dy + dz * dz;
        if (d2u < r2) r2 = d2u;
        double rad2 = d2u + eps2, rinv = 1.0 / sqrt(rad2), rinv2 = rinv * rinv;
        double mor = m[j] * rinv * rinv2;
        double drdv = 3.0 * (dx * dvx + dy * dvy + dz * dvz) * rinv2;
        ax += mor * dx;
        ay += mor * dy;
        az += mor * dz;
        jx += mor * (dvx - drdv * dx);
        jy += mor * (dvy - drdv * dy);
        jz += mor * (dvz - drdv * dz);
        /* pairwise timescales (squared): free-fall r^3/M and flyby r^2/v^2 */
        double dv2 = dvx * dvx + dvy * dvy + dvz * dvz;
        double tff = rad2 * rad2 * rinv / (mi + m[j]), tfl = rad2 / (dv2 + 1e-300);
        if (tff < t2) t2 = tff;
        if (tfl < t2) t2 = tfl;
    }
    a[0] = ax; a[1] = ay; a[2] = az;
    jk[0] = jx; jk[1] = jy; jk[2] = jz;
    *tmin2 = t2;
    *rm2 = r2;
}

/* force on body i from all bodies; the sum is shared over the threads
 * if par is set (used for small blocks and by the neighbour scheme) */
void fullsum(int i, int par, double a[3], double jk[3], double *tmin2, double *rm2)
{
    if (!par) {
        sumrange(i, 0, n, a, jk, tmin2, rm2);
        return;
    }
    double ax = 0, ay = 0, az = 0, jx = 0, jy = 0, jz = 0, t2 = HUGE_VAL, r2 = HUGE_VAL;
#pragma omp parallel reduction(+:ax,ay,az,jx,jy,jz) reduction(min:t2,r2)
    {
        int nt = omp_get_num_threads(), it = omp_get_thread_num();
        int j1 = (int)((int64_t)n * it / nt), j2 = (int)((int64_t)n * (it + 1) / nt);
        double pa[3], pj[3], pt, pr;
        sumrange(i, j1, j2, pa, pj, &pt, &pr);
        ax += pa[0]; ay += pa[1]; az += pa[2];
        jx += pj[0]; jy += pj[1]; jz += pj[2];
        if (pt < t2) t2 = pt;
        if (pr < r2) r2 = pr;
    }
    a[0] = ax; a[1] = ay; a[2] = az;
    jk[0] = jx; jk[1] = jy; jk[2] = jz;
    *tmin2 = t2;
    *rm2 = r2;
}

static void direct(int i, int par)
{
    double a[3], jk[3], t2, r2;
    fullsum(i, par, a, jk, &t2, &r2);
    tideadd(i, a, jk);
    for (int k = 0; k < 3; k++) {
        f0p[i][k] = a[k];
        f1p[i][k] = jk[k];
    }
    tpair2[i] = t2;
    rmin2[i] = r2;
}

/* accelerations f0p and jerks f1p on the active bodies ilist[0..nact-1]
 * from all bodies at the predicted positions rp, vp; also tpair2 and
 * rmin2. Big blocks (nact >= nompmin): the threads share out the
 * active bodies. Small blocks: one body at a time, the threads share
 * the sum over j (if n >= nompj). Newton's 3rd law is not used. */
void forces(void)
{
    if (useac) {
        if (nact >= nompmin) {
#pragma omp parallel for schedule(dynamic, 1)
            for (int ia = 0; ia < nact; ia++) acforce(ilist[ia], 1);
        } else {
            for (int ia = 0; ia < nact; ia++) acforce(ilist[ia], 0);
        }
    } else if (nact >= nompmin) {
#pragma omp parallel for schedule(static)
        for (int ia = 0; ia < nact; ia++) direct(ilist[ia], 0);
    } else {
        for (int ia = 0; ia < nact; ia++) direct(ilist[ia], n >= nompj);
    }
}

/* kinetic, gravitational potential (negative) and tidal potential
 * energy in internal units; also sets phi[i]. Only meaningful when all
 * bodies are synchronised. With a tide ke + pot + ptid is the Jacobi
 * energy. */
void energy(double *ke, double *pot, double *ptid)
{
    double k = 0.0, pt = 0.0, p = 0.0;
#pragma omp parallel for reduction(+:k,pt) schedule(static) if (n >= nompj)
    for (int i = 0; i < n; i++) {
        k += 0.5 * m[i] * dot3(v[i], v[i]);
        if (ltide)
            for (int c = 0; c < 3; c++) pt -= 0.5 * m[i] * r[i][c] * dot3(ttid[c], r[i]);
    }
#pragma omp parallel for schedule(static) reduction(+:p)
    for (int i = 0; i < n; i++) {
        double poti = 0.0;
        for (int j = 0; j < n; j++) {
            if (j == i) continue;
            double dx = r[j][0] - r[i][0], dy = r[j][1] - r[i][1], dz = r[j][2] - r[i][2];
            poti -= m[j] / sqrt(dx * dx + dy * dy + dz * dz + eps2);
        }
        phi[i] = poti;
        p += 0.5 * m[i] * poti;
    }
    *ke = k;
    *pot = p;
    *ptid = pt;
}

/* ---- timestep levels ------------------------------------------- */

/* smallest level k with dtmax/2^k <= dtw, limited to 0 <= k <= kmax */
int getlevel(double dtw)
{
    int k;
    if (dtw >= dtmax) k = 0;
    else {
        double x = ceil(log(dtmax / dtw) / log(2.0));
        k = (x > 1000.0) ? 1000 : (int)x;
    }
    if (k > kmax) {
        /* steps that want to be shorter than the shortest level are held
         * there and counted */
        k = kmax;
#pragma omp atomic
        nclamp++;
    }
    return k;
}

/* smallest level >= k whose step divides the integer time t8 */
int commlev(int k, int64_t t8)
{
    if (k < 0) k = 0;
    if (k > kmax) k = kmax;
    while (k < kmax && t8 % ((int64_t)1 << (kmax - k)) != 0) k++;
    return k;
}

/* forces on everyone at the start and each body's first level:
 * dt = min(etas*|a|/|j|, etap*tpair) */
void firststep(void)
{
    nact = n;
    for (int i = 0; i < n; i++) {
        ilist[i] = i;
        for (int k = 0; k < 3; k++) {
            rp[i][k] = r[i][k];
            vp[i][k] = v[i][k];
        }
    }
    tcur8 = 0;
    aclast = 1;
    if (useac) acsetrs0();
    forces();
#pragma omp parallel for schedule(static) if (n >= nompj)
    for (int i = 0; i < n; i++) {
        for (int k = 0; k < 3; k++) {
            f0[i][k] = f0p[i][k];
            f1[i][k] = f1p[i][k];
            f2[i][k] = f3[i][k] = f2e[i][k] = f3e[i][k] = 0.0;
        }
        tk[i] = 0;
        double amod = norm3(f0[i]), jmod = norm3(f1[i]);
        double dtw = etap * sqrt(tpair2[i]);
        if (amod > 0.0 && jmod > 0.0) dtw = fmin(dtw, etas * amod / jmod);
        lev[i] = getlevel(dtw);
        dtk[i] = (int64_t)1 << (kmax - lev[i]);
        if (useac) acstart(i);
    }
}

/* new level of body i after its step has been corrected:
 *   dt_A = sqrt(eta (|a||a2| + |a1|^2)/(|a1||a3| + |a2|^2))  (Aarseth),
 *          with a2, a3 at the END of the step
 *   dt_P = etap min_j(sqrt(r^3/(mi+mj)), r/v)                 (pairwise)
 * levels can go down (smaller step) by any amount, up by one level and
 * only when the current time is commensurate with 2 dt */
void newlevel(int i)
{
    double dt = (double)dtk[i] * tick, s[3];
    for (int k = 0; k < 3; k++) s[k] = f2[i][k] + f3[i][k] * dt;
    double amod = norm3(f0[i]), jmod = norm3(f1[i]), smod = norm3(s), cmod = norm3(f3[i]);
    double dtw = etap * sqrt(tpair2[i]);
    double denom = jmod * cmod + smod * smod;
    if (denom > 0.0) dtw = fmin(dtw, sqrt(eta * (amod * smod + jmod * jmod) / denom));
    int knew = getlevel(dtw);
    if (knew > lev[i]) lev[i] = knew;
    else if (knew < lev[i] && tk[i] % (2 * dtk[i]) == 0) lev[i]--;
    dtk[i] = (int64_t)1 << (kmax - lev[i]);
}
