/* ==============================================================
 * HENDOS - HErmite N-body Dynamics Of Systems - ac.c
 * AHMAD-COHEN NEIGHBOUR SCHEME (switch useac) for direct summation
 * (Ahmad & Cohen 1973, J. Comp. Phys., 12, 389; Hermite form after
 * Makino & Aarseth 1992, PASJ, 44, 141)
 *
 * The force on a star is split in two:
 *   irregular force   from its neighbours (the stars inside its
 *                     neighbour sphere, radius rs), which changes fast.
 *                     The Galactic tide is included in it.
 *   regular force     from all the other stars, which changes slowly
 * Each star has two timesteps. On most steps (irregular steps) only the
 * short sum over the neighbour list is done and the regular force is
 * taken from its last value and derivative (linear extrapolation). On
 * regular steps, which are 2, 4, 8... irregular steps apart, the full
 * sum over all stars is done, a new neighbour list is made and the
 * regular force is corrected.
 *
 * How it fits the Hermite scheme. The arrays f0, f1 are still the total
 * force and its derivative, so the predictor and the corrector in the
 * main program are unchanged:
 *   irregular step   f0p = (neighbour sum + tide) + fr + frd*(t - treg)
 *                    f1p = (its derivative) + frd
 *                    the regular part is linear in time, so it drops
 *                    out of the corrector's 2nd and 3rd derivatives,
 *                    which are those of the irregular force alone
 *   regular step     the same, with the OLD neighbour list, then the
 *                    Hermite corrector is applied to the regular force
 *                    over the whole regular step: its 2nd and 3rd
 *                    derivatives r2, r3 give the extra corrections
 *                      dreg  = r2 dtr^4/24 + r3 dtr^5/120   to the position
 *                      dvreg = r2 dtr^3/6  + r3 dtr^4/24    to the velocity
 *                    which the main program adds. Then the new list
 *                    replaces the old one and the force is split again.
 * The regular step is chosen from the regular force and its derivatives
 * with the Aarseth formula and the parameter etar. It is a block step
 * like the irregular one, never shorter than it, and regular times are
 * always whole multiples of the regular step, so every regular time is
 * also one of the star's irregular times.
 *
 * Neighbour list: stars inside rs, plus stars inside 3 rs that are
 * approaching and will be inside rs within one regular step. rs is
 * changed a little at every regular step to keep the number of
 * neighbours near nnbopt. kfmax places are kept free in every list for
 * stars that appear when a group is opened (see below).
 *
 * When bodies are created, removed or moved by the binary and merger
 * routines the lists are kept correct by acmoved, acremoved and
 * acadded; a new body has acok = 0 and gets a full force sum and
 * a new list the first time forces is called for it.
 * (Tested with a Python model of this scheme before it was written:
 * a 150-star Plummer sphere over 4 N-body times gave |dE/E| = 3e-8 with
 * etar = 0.005 and a quarter to a half of the full force sums.)
 *
 * OpenMP: for big blocks the active stars are shared over the threads;
 * for small blocks the full sums and list searches are shared instead.
 *
 * Bodies are 0..n-1 here, and the list of body i is
 * nbl[i*nnbmax .. i*nnbmax + nnb[i] - 1].
 * ============================================================== */
#include <omp.h>
#include "hendos.h"

#define NACMIN 200   /* below this many stars the scheme is switched off */

/* the neighbour scheme is switched off, and the user told, if there are
 * fewer than NACMIN stars: a neighbour list is then a large part of the
 * whole system, so the scheme costs as much as plain direct summation.
 * Called once the number of stars is known, before acalloc.
 *
 * Where NACMIN comes from. A Python model of this code was run on
 * equal-mass Plummer spheres of 30 to 800 stars with and without the
 * scheme, with etar = 0.005 so that the energy errors were the same,
 * and the pair interactions counted. Cost of the scheme / cost of
 * direct summation, if a pair in a neighbour sum costs w times a pair
 * in a full sum (neighbour sums jump about in memory and do not
 * vectorise, so w is more than 1; 2 to 3 is a guess):
 *        N      30    60   100   150   200   300   400   800
 *      w = 1   1.36  0.98  0.87  0.80  0.72  0.65  0.60  0.49
 *      w = 2   1.95  1.26  1.04  0.91  0.81  0.73  0.67  0.55
 *      w = 3   2.54  1.55  1.21  1.01  0.89  0.80  0.74  0.60
 * so the scheme breaks even at about 60 stars if w = 1, 110 if w = 2
 * and 150 if w = 3. w has not been measured (the code had not been
 * compiled when this was written).
 *
 * NACMIN is 200, a margin above the largest of those break-evens, at the
 * user's request (27 Nov 2026), to be safe.
 *
 * >>> CHANGE HERE (NACMIN) once the break-even has been measured with
 *     the compiled code: benchmark/bench_kira.py --skip-kira --modes
 *     direct neighbours at several N, compared at equal energy error
 *
 */
void accheck(int nstars)
{
    if (!useac || nstars >= NACMIN) return;
    useac = 0;
    printf("\n NOTE: the neighbour scheme (useac) has been SWITCHED OFF: there are only "
           "%d stars (fewer than %d).\n", nstars, NACMIN);
    printf("       For so few stars it is no faster than plain direct summation and is less\n"
           "       accurate. The run carries on with plain direct summation.\n\n");
}

void acalloc(void)
{
    if (nnbopt == 0) nnbopt = (int)fmax(10.0, round(pow(nalloc / 8.0, 0.75)));
    if (nnbopt > nalloc - 1) nnbopt = nalloc - 1;
    if (nnbmax == 0) nnbmax = 3 * nnbopt + KFMAX + 10;
    if (nnbmax < nnbopt + KFMAX + 2) die("nnbmax must be at least nnbopt + 8");
    size_t nv = sizeof(vec3) * (size_t)nalloc;
    fi = xalloc(nv); fid = xalloc(nv); fr = xalloc(nv); frd = xalloc(nv);
    fin = xalloc(nv); fidn = xalloc(nv); frn = xalloc(nv); frdn = xalloc(nv);
    dreg = xalloc(nv); dvreg = xalloc(nv);
    rs = xalloc(sizeof(double) * (size_t)nalloc);
    dtwr = xalloc(sizeof(double) * (size_t)nalloc);
    nnb = xalloc(sizeof(int) * (size_t)nalloc);
    levr = xalloc(sizeof(int) * (size_t)nalloc);
    nbl = xalloc(sizeof(int) * (size_t)nalloc * (size_t)nnbmax);
    treg = xalloc(sizeof(int64_t) * (size_t)nalloc);
    dtreg = xalloc(sizeof(int64_t) * (size_t)nalloc);
    lreg = xalloc((size_t)nalloc);
    acok = xalloc((size_t)nalloc);
}

/* first guess of the neighbour radius: a sphere holding nnbopt bodies
 * at the mean density inside the rms radius */
void acsetrs0(void)
{
    double s = 0.0, mt = 0.0;
    for (int i = 0; i < n; i++) {
        s += m[i] * dot3(r[i], r[i]);
        mt += m[i];
    }
    rs0 = sqrt(s / mt) * cbrt((double)nnbopt / (double)n);
}

/* force on body i from the nl bodies in list (same pair interaction
 * as sumrange in force.c - keep the two identical) */
static void nbsum(int i, int nl, const int *list, double a[3], double jk[3], double *tmin2,
                  double *rm2)
{
    double xi = rp[i][0], yi = rp[i][1], zi = rp[i][2];
    double vxi = vp[i][0], vyi = vp[i][1], vzi = vp[i][2], mi = m[i];
    double ax = 0, ay = 0, az = 0, jx = 0, jy = 0, jz = 0, t2 = HUGE_VAL, r2 = HUGE_VAL;
    for (int k = 0; k < nl; k++) {
        int j = list[k];
        if (j == i) continue;
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

/* new neighbour list of body i: the bodies inside rs[i], and those
 * inside 3 rs that are approaching and will be inside rs within one
 * regular step. If there are too many rs is reduced and the search
 * done again; for a new body (init) rs is also increased if there are
 * too few. Afterwards rs is changed by at most 25% towards the radius
 * that would hold nnbopt neighbours. The search is shared over the
 * threads if par is set. Returns the number of neighbours. */
static int aclist(int i, int par, int init, int *newl)
{
    int nlim = nnbmax - KFMAX, nfound = 0, nt = par ? nthreads : 1;
    int *scr = xalloc(sizeof(int) * (size_t)n), *cnt = xalloc(sizeof(int) * (size_t)nt);
    double dtr = init ? 0.0 : (double)dtreg[i] * tick;
    if (rs[i] <= 0.0) rs[i] = rs0;
    for (int it = 1; it <= 30; it++) {
        double rs2 = rs[i] * rs[i];
#pragma omp parallel for schedule(static, 1) if (par) num_threads(nt)
        for (int ic = 0; ic < nt; ic++) {
            int j1 = (int)((int64_t)n * ic / nt), j2 = (int)((int64_t)n * (ic + 1) / nt);
            int k = j1;
            for (int j = j1; j < j2; j++) {
                if (j == i) continue;
                double dx = rp[j][0] - rp[i][0], dy = rp[j][1] - rp[i][1],
                       dz = rp[j][2] - rp[i][2];
                double d2 = dx * dx + dy * dy + dz * dz;
                if (d2 >= rs2) {
                    if (init || d2 >= 9.0 * rs2) continue;
                    double dvx = vp[j][0] - vp[i][0], dvy = vp[j][1] - vp[i][1],
                           dvz = vp[j][2] - vp[i][2];
                    if (dx * dvx + dy * dvy + dz * dvz >= 0.0) continue;
                    double fx = dx + dvx * dtr, fy = dy + dvy * dtr, fz = dz + dvz * dtr;
                    if (fx * fx + fy * fy + fz * fz >= rs2) continue;
                }
                scr[k++] = j;
            }
            cnt[ic] = k - j1;
        }
        nfound = 0;
        for (int ic = 0; ic < nt; ic++) nfound += cnt[ic];
        if (it == 30) break;
        if (nfound > nlim) {
            rs[i] *= fmax(0.5, cbrt((double)nnbopt / nfound));
            continue;
        }
        if (init && nfound < nnbopt / 2 && nfound < n - 1) {
            rs[i] *= fmin(2.0, cbrt((double)nnbopt / (nfound > 1 ? nfound : 1)));
            continue;
        }
        break;
    }
    /* join the pieces up (never more than nlim) */
    int nl = 0;
    for (int ic = 0; ic < nt; ic++) {
        int j1 = (int)((int64_t)n * ic / nt);
        for (int k = 0; k < cnt[ic] && nl < nlim; k++) newl[nl++] = scr[j1 + k];
    }
    free(scr);
    free(cnt);
    rs[i] *= fmin(1.25, fmax(0.8, cbrt((double)nnbopt / (nl > 1 ? nl : 1))));
    return nl;
}

/* force and its derivative on active body i at the block time tcur8,
 * into f0p, f1p, by an irregular step, a regular step or (for a body
 * with acok unset) a first full sum. inpar is set when this is called
 * from inside a parallel loop over the active bodies. Thread safe:
 * only body i's entries are written. */
void acforce(int i, int inpar)
{
    int par = !inpar && n >= nompj, init = !acok[i];
    double an[3], jn[3], at[3], jt[3], ao[3], jo[3], a[3], j[3], tn2, rn2, ta2, ra2, dtr;
    int *list = nbl + (size_t)i * (size_t)nnbmax;

    if (!init && tcur8 != treg[i] + dtreg[i]) {
        /* ---- irregular step: neighbours only */
        nbsum(i, nnb[i], list, an, jn, &tn2, &rn2);
        dtr = (double)(tcur8 - treg[i]) * tick;
        for (int k = 0; k < 3; k++) {
            a[k] = an[k] + fr[i][k] + frd[i][k] * dtr;
            j[k] = jn[k] + frd[i][k];
        }
        tideadd(i, a, j);
        for (int k = 0; k < 3; k++) {
            fin[i][k] = a[k] - fr[i][k] - frd[i][k] * dtr;
            fidn[i][k] = j[k] - frd[i][k];
            f0p[i][k] = a[k];
            f1p[i][k] = j[k];
        }
        tpair2[i] = tn2;
        rmin2[i] = rn2;
        lreg[i] = 0;
#pragma omp atomic
        nirregular++;
        return;
    }

    /* ---- regular step or first sum: all the bodies */
    fullsum(i, par, at, jt, &ta2, &ra2);
    tpair2[i] = ta2;
    rmin2[i] = ra2;
#pragma omp atomic
    nregular++;

    if (!init) {
        /* the end of the irregular step, with the old list */
        double fro[3], frdo[3], r2[3], r3[3], r2e[3];
        dtr = (double)dtreg[i] * tick;
        nbsum(i, nnb[i], list, ao, jo, &tn2, &rn2);
        for (int k = 0; k < 3; k++) {
            a[k] = ao[k] + fr[i][k] + frd[i][k] * dtr;
            j[k] = jo[k] + frd[i][k];
        }
        tideadd(i, a, j);
        for (int k = 0; k < 3; k++) {
            f0p[i][k] = a[k];
            f1p[i][k] = j[k];
            /* Hermite corrector of the regular force (old list) over the regular step */
            fro[k] = at[k] - ao[k];
            frdo[k] = jt[k] - jo[k];
            double a1 = fr[i][k] - fro[k], a2 = frd[i][k] + frdo[k];
            r2[k] = 2.0 * (-3.0 * a1 - (frd[i][k] + a2) * dtr) / (dtr * dtr);
            r3[k] = 6.0 * (2.0 * a1 + a2 * dtr) / (dtr * dtr * dtr);
            dreg[i][k] = r2[k] * pow(dtr, 4) / 24.0 + r3[k] * pow(dtr, 5) / 120.0;
            dvreg[i][k] = r2[k] * pow(dtr, 3) / 6.0 + r3[k] * pow(dtr, 4) / 24.0;
            r2e[k] = r2[k] + r3[k] * dtr;
        }
        /* regular timestep wanted (Aarseth formula on the regular force) */
        double den = norm3(frdo) * norm3(r3) + dot3(r2e, r2e);
        dtwr[i] = den > 0.0
                      ? sqrt(etar * (norm3(fro) * norm3(r2e) + dot3(frdo, frdo)) / den)
                      : dtmax;
        lreg[i] = 1;
        /* with corrector iterations the new list is only made on the last one */
        if (!aclast) return;
    }

    /* new neighbour list, and the force split with it */
    int *newl = xalloc(sizeof(int) * (size_t)nnbmax);
    int nl = aclist(i, par, init, newl);
    nbsum(i, nl, newl, an, jn, &tn2, &rn2);
    for (int k = 0; k < 3; k++) {
        a[k] = at[k];
        j[k] = jt[k];
    }
    tideadd(i, a, j);
    for (int k = 0; k < 3; k++) {
        frn[i][k] = at[k] - an[k];
        frdn[i][k] = jt[k] - jn[k];
        fin[i][k] = a[k] - frn[i][k];
        fidn[i][k] = j[k] - frdn[i][k];
    }
    nnb[i] = nl;
    memcpy(list, newl, sizeof(int) * (size_t)nl);
    free(newl);

    if (init) {
        /* a new body: everything is set now (levr, dtreg are set by acstart) */
        for (int k = 0; k < 3; k++) {
            fi[i][k] = fin[i][k];
            fid[i][k] = fidn[i][k];
            fr[i][k] = frn[i][k];
            frd[i][k] = frdn[i][k];
            f0p[i][k] = a[k];
            f1p[i][k] = j[k];
        }
        treg[i] = tcur8;
        lreg[i] = 0;
        acok[i] = 1;
    }
}

/* start the regular step of a new body i, whose irregular level and
 * time have just been set: the first regular step is the irregular step */
void acstart(int i)
{
    levr[i] = lev[i];
    dtreg[i] = dtk[i];
    treg[i] = tk[i];
}

/* the step of body i has been accepted (f0, f1 set from f0p, f1p and
 * tk[i] to the new time): keep the new irregular force, and after a
 * regular step the new regular force too; f0, f1 are then the exact
 * totals from the full sum. Called before newlevel. */
void acaccept(int i)
{
    for (int k = 0; k < 3; k++) {
        fi[i][k] = fin[i][k];
        fid[i][k] = fidn[i][k];
    }
    if (lreg[i]) {
        for (int k = 0; k < 3; k++) {
            fr[i][k] = frn[i][k];
            frd[i][k] = frdn[i][k];
            f0[i][k] = fin[i][k] + frn[i][k];
            f1[i][k] = fidn[i][k] + frdn[i][k];
        }
        treg[i] = tk[i];
    }
}

/* after newlevel: choose the regular level if this was a regular step
 * (down by any amount, up one level at a time, and its step must
 * divide the current time), and keep the irregular step no longer
 * than the regular step */
void aclevel(int i)
{
    if (lreg[i]) {
        int kr = getlevel(dtwr[i]);
        if (kr < levr[i]) kr = levr[i] - 1;
        levr[i] = commlev(kr, tk[i]);
        dtreg[i] = (int64_t)1 << (kmax - levr[i]);
        lreg[i] = 0;
    }
    if (lev[i] < levr[i]) {
        lev[i] = levr[i];
        dtk[i] = (int64_t)1 << (kmax - lev[i]);
    }
}

/* ---- keeping the lists right when bodies are moved, removed, made -- */

void acmoverow(int src, int dst)
{
    memcpy(fi[dst], fi[src], sizeof(vec3));
    memcpy(fid[dst], fid[src], sizeof(vec3));
    memcpy(fr[dst], fr[src], sizeof(vec3));
    memcpy(frd[dst], frd[src], sizeof(vec3));
    rs[dst] = rs[src];
    nnb[dst] = nnb[src];
    memcpy(nbl + (size_t)dst * (size_t)nnbmax, nbl + (size_t)src * (size_t)nnbmax,
           sizeof(int) * (size_t)nnb[src]);
    levr[dst] = levr[src];
    treg[dst] = treg[src];
    dtreg[dst] = dtreg[src];
    acok[dst] = acok[src];
    lreg[dst] = 0;
}

/* body src has been moved to slot dst: change its index in every list */
void acmoved(int src, int dst)
{
#pragma omp parallel for schedule(static) if (n >= nompj)
    for (int i = 0; i < n; i++) {
        int *list = nbl + (size_t)i * (size_t)nnbmax;
        for (int k = 0; k < nnb[i]; k++)
            if (list[k] == src) list[k] = dst;
    }
}

/* body j is being removed (it has joined a group or merged, and body
 * irep now stands for it; irep = -1 if nothing does): take j out of
 * every list. In a list that held j but not irep, irep takes its place. */
void acremoved(int j, int irep)
{
#pragma omp parallel for schedule(static) if (n >= nompj)
    for (int i = 0; i < n; i++) {
        int *list = nbl + (size_t)i * (size_t)nnbmax, kj = -1, has = 0;
        for (int k = 0; k < nnb[i]; k++) {
            if (list[k] == j) kj = k;
            if (list[k] == irep) has = 1;
        }
        if (kj < 0) continue;
        if (irep >= 0 && irep != i && !has) list[kj] = irep;
        else list[kj] = list[--nnb[i]];
    }
}

/* body jn has appeared beside body i0 (a group was opened): add jn to
 * every list that holds i0. Lists keep KFMAX places free for this. */
void acadded(int jn, int i0)
{
#pragma omp parallel for schedule(static) if (n >= nompj)
    for (int i = 0; i < n; i++) {
        if (i == jn || i == i0) continue;
        int *list = nbl + (size_t)i * (size_t)nnbmax, has = 0;
        for (int k = 0; k < nnb[i]; k++)
            if (list[k] == i0) has = 1;
        if (has && nnb[i] < nnbmax) list[nnb[i]++] = jn;
    }
}

/* the nm bodies in list (all at integer time t8, positions rp, vp) are
 * about to become one body: a group's c.o.m., or a merged star, in the
 * slot of one of them. Every list that holds any of them will then hold
 * that one body, which carries the mass of all of them (see acremoved).
 * So in a list that held only some of them, the others move from the
 * regular force, where they still are, into the irregular force; their
 * pull is taken out of the regular force here, and put into the
 * irregular force kept for the corrector, so that it is counted once.
 * (Found in cold fractals: stars in sub-clumps whose lists held one star
 * of a pair being made a group felt the other star twice until their
 * next regular step - 0.1% of regular steps found the regular force
 * wrong by more than 10% of the total force.) */
void acjoining(int nm, const int *list, int64_t t8)
{
#pragma omp parallel for schedule(dynamic, 16) if (n >= nompj)
    for (int i = 0; i < n; i++) {
        int *li = nbl + (size_t)i * (size_t)nnbmax, inl[KFMAX], nin = 0, skip = 0;
        if (!acok[i]) continue;
        for (int l = 0; l < nm; l++) {
            if (list[l] == i) skip = 1;
            inl[l] = 0;
            for (int k = 0; k < nnb[i]; k++)
                if (li[k] == list[l]) inl[l] = 1;
            nin += inl[l];
        }
        if (skip || nin == 0 || nin == nm) continue;
        double tr = (double)(t8 - treg[i]) * tick, ti = (double)(t8 - tk[i]) * tick;
        for (int l = 0; l < nm; l++) {
            if (inl[l]) continue;
            int j = list[l];
            double d[3], dv[3];
            for (int q = 0; q < 3; q++) {
                d[q] = rp[j][q] - rp[i][q];
                dv[q] = vp[j][q] - vp[i][q];
            }
            double r2 = dot3(d, d) + eps2, rinv = 1.0 / sqrt(r2), mr3 = m[j] * rinv * rinv * rinv;
            double rv = 3.0 * dot3(d, dv) / r2;
            for (int q = 0; q < 3; q++) {
                double a = mr3 * d[q], jk = mr3 * (dv[q] - rv * d[q]);
                /* the regular force is fr + frd (t - treg), linear */
                fr[i][q] -= a - jk * tr;
                frd[i][q] -= jk;
                /* the irregular force at the star's last step, tk */
                fi[i][q] += a - jk * ti;
                fid[i][q] += jk;
            }
        }
    }
}
