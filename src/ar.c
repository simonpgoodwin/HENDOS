/* ==============================================================
 * HENDOS - HErmite N-body Dynamics Of Systems - ar.c
 * REGULARISED FEW-BODY GROUPS (bmode 2): perturbed binaries and
 * compact multiples of 2 to KFMAX stars
 *
 * WHAT IS IN HERE. The internal motion of a group of bmode 2 (see
 * groups.c for when groups are made and broken up) is integrated here
 * by ALGORITHMIC REGULARISATION (AR): the leapfrog of Mikkola &
 * Tanikawa (1999) and Preto & Tremaine (1999) in its logarithmic-
 * Hamiltonian form, made accurate by Gragg-Bulirsch-Stoer (GBS)
 * extrapolation (Mikkola & Aarseth 2002; it is the integrator of AR-CHAIN
 * and, in a slightly different form, of TSUNAMI). With a time
 * transformation dt = ds / U the leapfrog integrates two-body collisions
 * and close approaches exactly in the limit of small eccentricity
 * changes, and it is regular: a pericentre of 1e-6 of the semi-major axis
 * costs no more steps than a circular orbit. So unlike the Hermite
 * integrator it does not need tiny steps at pericentre, and unlike KS it
 * works for any number of stars with no coordinate changes.
 *
 * The leapfrog, with a step h in the new independent variable s:
 *    drift D(h):  dt = h / (T + B);  t += dt;  x += dt v
 *    kick  K(h):  dt = h / U;        v += dt (a_int + a_ext)
 *                                    B -= dt sum m v.a_ext  (v the mean of old and new)
 * T kinetic, U = sum m_i m_j / r_ij the (positive) potential, and B = U - T
 * at the start, i.e. minus the internal energy, carried along with the
 * work done by the outside forces. One step of s is
 *    D(h/2) K(h) D(h) K(h) ... D(h) K(h) D(h/2)    (nsub kicks, h = H/nsub)
 * and GBS does it with nsub = 2, 4, 6, ..., 16 and extrapolates to
 * h -> 0 (the leapfrog error is even in h). A step is accepted when the
 * last two extrapolations agree to artol (relative). The extrapolation
 * can reach 16th order, so the steps are long: about 10 to 20 per orbit.
 *
 * THE OUTSIDE WORLD. The main integration sees a group as one body at its
 * centre of mass (c.o.m.). Here the stars feel, besides each other, the
 * tidal pull of a list of PERTURBERS: every body whose own perturbation
 *    gamma_j = 2 s^3 m_j / (M d_j^3)          (s group size, M its mass,
 *                                              d_j distance from the c.o.m.)
 * is above gpert (standard 1e-7; at most NPMAX bodies, the strongest). The
 * perturbers and the c.o.m. are predicted at any time from their own
 * Taylor series (the same 5th-order prediction as the main integrator),
 * and each star feels a_ext(x_i) minus the mass-weighted mean of a_ext
 * over the group (the mean moves the c.o.m., which the main integration
 * does). The list is made by a search of all the bodies each time the
 * c.o.m. starts a new step; at the start of the other stretches of
 * integration (see WHEN) only the bodies on it are looked at again.
 * A STRONG perturber (gamma_j > GHORFAC gpert, 1e-5 as standard) is
 * followed exactly, star by star. The weak ones only through their summed
 * TIDAL TENSOR T (a = T.x for a star at x from the c.o.m.), worked out at
 * the start and end of the stretch and interpolated linearly in time: in a
 * clumpy cluster a group can have dozens of perturbers, and summing them
 * exactly was 85% of the time. What the tensor leaves out is of order
 * gamma_j (s/d_j), below 1e-5 * 0.05; making the tensor take perturbers up
 * to gamma_j = 1e-4 or 1e-3 was tried, and the energy error of a cold
 * fractal grew in proportion. The perturbers feel the stars exactly
 * whichever they are (see THE OTHER WAY ROUND).
 * A Taylor series is only good over its own step: the series of a star in
 * a tight orbit, with steps of hours, is nonsense (parsecs out) after a
 * few hundred years. So, as in the main integrator, no prediction is used
 * past the end of the body's step: a group is never integrated past the
 * end of the step of its c.o.m. or of a strong perturber, and a weak
 * perturber whose step ends first is held where it is at the start, for
 * the whole stretch (its true motion, of at most the size of its own orbit
 * or its drift, hardly matters because it is weak - and its force must
 * stay smooth, or the GBS extrapolation fails).
 *
 * THE OTHER WAY ROUND. The perturbers must feel the group's stars, not
 * just a point mass at its c.o.m., and the c.o.m. must feel the mean
 * pull on its stars, not the pull at its c.o.m. - or the energy and
 * momentum exchanged with the perturbers are wrong (by about (s/d)^2 of
 * the force: a few per cent for a strong perturber). So arcorrect, called
 * at the end of every force calculation, adds to the force and jerk on
 *   - each active perturber j of group k:
 *         sum_i m_i g(x_i - x_j) - M g(X - x_j)
 *   - each active c.o.m. X of a group:
 *         (1/M) sum_i m_i sum_j m_j g(x_j - x_i) - sum_j m_j g(x_j - X)
 * (g(d) = d/|d|^3 and its time derivative for the jerk). Then the
 * Hamiltonian of the whole system is exactly
 *     H = (point masses) + sum_groups (internal) + sum_groups sum_j [U_resolved - U_point]
 * and diagnose adds the last sum (arpotcorr) so that the energy check is
 * honest.
 *
 * WHEN. At the start of each block, a group that has not been integrated
 * as far as the block time is integrated AHEAD (aradvance), from where it
 * got to, to its HORIZON: the end of the current step of its c.o.m. or of
 * one of its strong perturbers, whichever comes first (all known at that
 * moment, in the integer block clock). The states at the ends of the AR
 * steps (positions, velocities, accelerations) are kept,
 * and the stars at any time between are found by quintic Hermite
 * interpolation (arstate) - that is what arcorrect uses for the
 * perturbers' forces. The steps are kept below ARCAP times the shortest
 * pair dynamical time so that the interpolation is good to about 1e-6
 * (it only enters the corrections above, which are themselves small).
 *
 * LANDING. The end of each stretch must be hit exactly in time, but the
 * steps are in s, not t. The last step is aimed at the end (dt/ds is
 * about 1/U) and then corrected by secant iteration in s until it ends
 * within LANDTOL (1e-6) of the length of that last step from the end; the
 * rest is covered by a second-order Taylor step in time, whose error (of
 * order LANDTOL^3) is far below anything else. Iterating the secant to
 * 1e-14 instead took six extrapolated steps per landing, and with short
 * stretches (strong perturbers with short steps) that was 84% of all the
 * work.
 *
 * TIMESTEPS. The corrections vary on the group's internal timescale. Once
 * a body has taken a step with them, the Aarseth criterion sees how fast
 * its force changes and sets its steps. But the step a body is on when a
 * group is made was chosen without them, and can be far too long: found
 * with an e = 0.9999 binary made a group at the start, sitting at
 * apocentre and about to plunge, with a star falling on to it from six
 * binary sizes away - the star's first step was half the binary's orbit,
 * and the energy error 2e-6 whatever eta was. So when a group is made
 * (arcapnew), the bodies that feel it and are at the start of a step, if
 * their correction is more than CAPFRAC of their force, have that step
 * cut to below
 *    etap sqrt(s^3/M)      the group's own dynamical time (s its diameter),
 *                          just as etap limits the step against any close
 *                          body - this is what a binary at apocentre needs,
 *                          which no derivative shows yet
 *    sqrt(eta) |da|/|dj|   (da, dj the correction and its derivative)
 * which brought the error to 2e-10. Applying the same cap at every step
 * as well was tried: it made a cold fractal 15 times slower for no gain.
 * What IS needed at every step is the pairwise limit (etap times the
 * free-fall or flyby time) measured against the group's STARS: the main
 * integration's own pair sum only sees the c.o.m., so a star flying
 * through or close past a group took steps set by its distance to the
 * c.o.m. and was integrated far too coarsely (found in a cold fractal, an
 * energy jump of 1.7e-4 in one passage). arcorrect therefore lowers
 * tpair2 (and rmin2, for the collision check) of every active body from
 * the distances it has just worked out for the corrections.
 *
 * COLLISIONS inside a group: the smallest separation of each pair is
 * watched during the kicks; if it falls below rcoll the two stars are
 * merged at the end of that AR step (mass and momentum kept, orbital
 * energy to emerged, a line in <runid>_mergers.dat). The merged star
 * keeps the lower slot; the other becomes a massless ghost until the
 * group is next brought up to date (arshift), when it is removed.
 *
 * Every group here has its own work and its own memory, so the groups
 * are integrated in parallel (aradvance); mergers are done one at a time.
 * ============================================================== */
#include <omp.h>
#include "hendos.h"

#define NSEQ 8           /* GBS: leapfrogs with 2, 4, ..., 2*NSEQ kicks */
#define ARCAP 0.5        /* >>> CHANGE HERE: longest AR step, in units of the shortest
                            pair dynamical time sqrt(r^3/(m1+m2)) (about 12 steps per orbit
                            at 0.5); larger is faster, but the interpolation for the
                            perturbers' forces gets worse as ARCAP^6 */
#define GHORFAC 100.0    /* >>> CHANGE HERE: a perturber is 'strong' if its gamma_j >
                            GHORFAC gpert; the group is never integrated past the end of a
                            strong perturber's step (see THE OUTSIDE WORLD above) */
#define LANDTOL 1.0e-6   /* a landing stops within this fraction of the last stretch
                            of the end (see LANDING above) */
#define CAPFRAC 1.0e-4   /* >>> CHANGE HERE: a correction larger than this fraction of
                            a body's force limits its step when a group is made (see
                            TIMESTEPS above) */
#define NYMAX (2 + 6 * KFMAX)

/* the dense output of a group: node q at time t[q] has, for each of its
 * nm stars, m, x[3], v[3], a[3] (relative to the c.o.m.) */
typedef struct {
    int nn, cap, nm;
    size_t caps;                  /* doubles allocated in s */
    double *t, *s;
} dense_t;
static dense_t *dn = NULL;
static int ndn = 0;

/* one group, set up for integration */
typedef struct {
    int nm, np;
    double m[KFMAX], mt;
    double pm[NPMAX], pc[NPMAX][6][3], pt0[NPMAX];    /* strong perturbers: mass, Taylor series */
    int nweak;                                          /* weak perturbers: */
    double tt[2][3][3], tw0, tw1;                       /* their tidal tensor at tw0, tw1 */
    double cc[6][3], ct0;                               /* the c.o.m.: Taylor series */
    double t0;                                          /* internal time at s = 0 */
} arsys;

/* perturber index: which groups each body perturbs (for arcorrect) */
static int *phead = NULL, *pk = NULL, *pnext = NULL, *pmarked = NULL;
/* when each group's perturber list was last made by a search of all the
 * bodies: valid, and the c.o.m.'s step start then */
static int *plok = NULL;
static int64_t *plcm8 = NULL;
static double *plbase = NULL;
static int npool = 0, nmarked = 0, pdirty = 1, phalloc = 0;

static double tolt(void) { return 0.5 * tick; }
static void arcapnew(int k, int64_t t8);

void ardirty(void) { pdirty = 1; }

void aralloc(int ng)
{
    dn = xalloc(sizeof(dense_t) * (size_t)ng);
    plok = xalloc(sizeof(int) * (size_t)ng);
    plcm8 = xalloc(sizeof(int64_t) * (size_t)ng);
    plbase = xalloc(sizeof(double) * (size_t)ng);
    for (int k = 0; k < ng; k++) {
        dn[k].nn = dn[k].cap = dn[k].nm = 0;
        dn[k].caps = 0;
        plok[k] = 0;
        dn[k].t = dn[k].s = NULL;
        npl[k] = 0;
    }
    ndn = ng;
}

/* forget the dense output and perturbers of group k */
void arreset(int k)
{
    dn[k].nn = 0;
    npl[k] = 0;
    plok[k] = 0;
    pdirty = 1;
}

/* group 'from' is being moved to index 'to' of the group table */
void armovegroup(int from, int to)
{
    dense_t tmp = dn[to];
    dn[to] = dn[from];
    dn[from] = tmp;               /* (swapped, so no memory is lost) */
    dn[from].nn = 0;
    npl[to] = npl[from];
    memcpy(plist[to], plist[from], sizeof(int) * (size_t)npl[from]);
    bend[to] = bend[from];
    bend8[to] = bend8[from];
    plok[to] = plok[from];
    plcm8[to] = plcm8[from];
    plbase[to] = plbase[from];
    plok[from] = 0;
    npl[from] = 0;
    pdirty = 1;
}

/* body 'from' is now in slot 'to' (to = -1: it has gone): change the
 * perturber lists to match */
void premap(int from, int to)
{
    for (int k = 0; k < nbin; k++) {
        if (bmode[k] != 2) continue;
        int q = 0;
        for (int e = 0; e < npl[k]; e++) {
            int j = plist[k][e];
            if (j == from) j = to;
            if (j < 0 || j == bcm[k]) continue;
            int dup = 0;
            for (int e2 = 0; e2 < q; e2++)
                if (plist[k][e2] == j) dup = 1;
            if (!dup) plist[k][q++] = j;
        }
        npl[k] = q;
    }
    pdirty = 1;
}

/* ---- Taylor series of body j about its last step ---------------- */
static void tset(int j, double c[6][3], double *t0)
{
    for (int q = 0; q < 3; q++) {
        c[0][q] = r[j][q];
        c[1][q] = v[j][q];
        c[2][q] = f0[j][q] / 2.0;
        c[3][q] = f1[j][q] / 6.0;
        c[4][q] = f2e[j][q] / 24.0;
        c[5][q] = f3e[j][q] / 120.0;
    }
    *t0 = tbase + (double)tk[j] * tick;
}
static inline void tpos(const double c[6][3], double dt, double x[3])
{
    for (int q = 0; q < 3; q++)
        x[q] = c[0][q] + dt * (c[1][q] + dt * (c[2][q] + dt * (c[3][q] + dt * (c[4][q]
               + dt * c[5][q]))));
}

/* ---- forces ------------------------------------------------------ */
/* mutual accelerations a of the stars at x; returns U = sum m m / r, the
 * smallest separation and its pair */
static double accint(const arsys *s, const double *x, double *a, double *rmin, int *ic, int *jc)
{
    int nm = s->nm;
    double u = 0.0, rm = HUGE_VAL;
    memset(a, 0, sizeof(double) * 3 * (size_t)nm);
    for (int l = 0; l < nm; l++) {
        if (s->m[l] == 0.0) continue;
        for (int l2 = l + 1; l2 < nm; l2++) {
            if (s->m[l2] == 0.0) continue;
            double d[3] = {x[3 * l2] - x[3 * l], x[3 * l2 + 1] - x[3 * l + 1],
                           x[3 * l2 + 2] - x[3 * l + 2]};
            double r2 = dot3(d, d), rr = sqrt(r2), ir3 = 1.0 / (r2 * rr);
            for (int q = 0; q < 3; q++) {
                a[3 * l + q] += s->m[l2] * d[q] * ir3;
                a[3 * l2 + q] -= s->m[l] * d[q] * ir3;
            }
            u += s->m[l] * s->m[l2] / rr;
            if (rr < rm) {
                rm = rr;
                if (ic) { *ic = l; *jc = l2; }
            }
        }
    }
    if (rmin) *rmin = rm;
    return u;
}

/* tidal accelerations of the perturbers on the stars at x, at time
 * s->t0 + tau, less their mass-weighted mean */
static void accext(const arsys *s, double tau, const double *x, double *a)
{
    int nm = s->nm;
    memset(a, 0, sizeof(double) * 3 * (size_t)nm);
    if (s->np == 0 && s->nweak == 0) return;
    double t = s->t0 + tau, xc[3], abar[3] = {0, 0, 0};
    if (s->nweak > 0) {
        double w = s->tw1 > s->tw0 ? (t - s->tw0) / (s->tw1 - s->tw0) : 0.0, tt[3][3];
        for (int q = 0; q < 3; q++)
            for (int q2 = 0; q2 < 3; q2++)
                tt[q][q2] = s->tt[0][q][q2] + w * (s->tt[1][q][q2] - s->tt[0][q][q2]);
        for (int l = 0; l < nm; l++)
            for (int q = 0; q < 3; q++) a[3 * l + q] = dot3(tt[q], x + 3 * l);
    }
    tpos(s->cc, t - s->ct0, xc);
    for (int j = 0; j < s->np; j++) {
        double pj[3];
        tpos(s->pc[j], t - s->pt0[j], pj);
        for (int q = 0; q < 3; q++) pj[q] -= xc[q];
        for (int l = 0; l < nm; l++) {
            double d[3] = {pj[0] - x[3 * l], pj[1] - x[3 * l + 1], pj[2] - x[3 * l + 2]};
            double r2 = dot3(d, d), f = s->pm[j] / (r2 * sqrt(r2));
            for (int q = 0; q < 3; q++) a[3 * l + q] += f * d[q];
        }
    }
    for (int l = 0; l < nm; l++)
        for (int q = 0; q < 3; q++) abar[q] += s->m[l] * a[3 * l + q];
    for (int l = 0; l < nm; l++)
        for (int q = 0; q < 3; q++) a[3 * l + q] -= abar[q] / s->mt;
}

static double kinetic(const arsys *s, const double *vv)
{
    double t = 0.0;
    for (int l = 0; l < s->nm; l++)
        t += 0.5 * s->m[l] * (vv[3 * l] * vv[3 * l] + vv[3 * l + 1] * vv[3 * l + 1]
                              + vv[3 * l + 2] * vv[3 * l + 2]);
    return t;
}

/* ---- the leapfrog and GBS ------------------------------------------
 * state y: y[0] time since s->t0, y[1] B, y[2..] x (3 nm), then v (3 nm) */
static void drift(const arsys *s, double *y, double h)
{
    int nm = s->nm;
    double *x = y + 2, *vv = y + 2 + 3 * nm;
    double dt = h / (kinetic(s, vv) + y[1]);
    y[0] += dt;
    for (int q = 0; q < 3 * nm; q++) x[q] += dt * vv[q];
}
static void kick(const arsys *s, double *y, double h, double *rmin, int *ic, int *jc)
{
    int nm = s->nm;
    double *x = y + 2, *vv = y + 2 + 3 * nm, ai[3 * KFMAX], ae[3 * KFMAX], rm;
    int i1 = 0, j1 = 0;
    double u = accint(s, x, ai, &rm, &i1, &j1);
    accext(s, y[0], x, ae);
    if (rm < *rmin) {
        *rmin = rm;
        *ic = i1;
        *jc = j1;
    }
    double dt = h / u, db = 0.0;
    for (int l = 0; l < nm; l++)
        for (int q = 0; q < 3; q++) {
            int c = 3 * l + q;
            double vn = vv[c] + dt * (ai[c] + ae[c]);
            db += s->m[l] * 0.5 * (vv[c] + vn) * ae[c];
            vv[c] = vn;
        }
    y[1] -= dt * db;
}
static void leapfrog(const arsys *s, const double *y0, double hs, int nsub, double *y,
                     double *rmin, int *ic, int *jc)
{
    int ny = 2 + 6 * s->nm;
    double h = hs / nsub;
    memcpy(y, y0, sizeof(double) * (size_t)ny);
    drift(s, y, 0.5 * h);
    for (int q = 0; q < nsub; q++) {
        kick(s, y, h, rmin, ic, jc);
        drift(s, y, q == nsub - 1 ? 0.5 * h : h);
    }
}

/* one GBS step of length hs in s from y0: returns the column at which it
 * converged (>= 2), or -1; yout the result */
static int gbs(const arsys *s, const double *y0, double hs, double *yout, double *rmin,
               int *ic, int *jc)
{
    static const int nseq[NSEQ] = {2, 4, 6, 8, 10, 12, 14, 16};
    int nm = s->nm, ny = 2 + 6 * nm;
    double tab[NSEQ][NSEQ][NYMAX], rs = 0.0, vs = 0.0, ai[3 * KFMAX];
    for (int l = 0; l < nm; l++) {
        if (s->m[l] == 0.0) continue;
        rs = fmax(rs, norm3(y0 + 2 + 3 * l));
        vs = fmax(vs, norm3(y0 + 2 + 3 * nm + 3 * l));
    }
    double u0 = accint(s, y0 + 2, ai, NULL, NULL, NULL);
    for (int kk = 0; kk < NSEQ; kk++) {
        double rm = HUGE_VAL;
        int i1 = 0, j1 = 0;
        leapfrog(s, y0, hs, nseq[kk], tab[kk][0], &rm, &i1, &j1);
        for (int j = 1; j <= kk; j++) {
            double rat = (double)nseq[kk] / nseq[kk - j], den = rat * rat - 1.0;
            for (int c = 0; c < ny; c++)
                tab[kk][j][c] = tab[kk][j - 1][c] + (tab[kk][j - 1][c] - tab[kk - 1][j - 1][c]) / den;
        }
        if (kk < 2) continue;
        const double *a = tab[kk][kk], *b = tab[kk][kk - 1];
        double err = fabs(a[0] - b[0]) / fmax(fabs(a[0]), 1e-300);
        err = fmax(err, fabs(a[1] - b[1]) / u0);
        for (int l = 0; l < nm; l++) {
            if (s->m[l] == 0.0) continue;
            for (int q = 0; q < 3; q++) {
                err = fmax(err, fabs(a[2 + 3 * l + q] - b[2 + 3 * l + q]) / rs);
                err = fmax(err, fabs(a[2 + 3 * nm + 3 * l + q] - b[2 + 3 * nm + 3 * l + q]) / vs);
            }
        }
        if (err < artol) {
            memcpy(yout, a, sizeof(double) * (size_t)ny);
            *rmin = rm;
            *ic = i1;
            *jc = j1;
            return kk;
        }
    }
    return -1;
}

/* advance the state y by a tiny time dt (positive or negative) by a
 * second-order Taylor step, carrying B along: the last bit of a landing */
static void landfix(const arsys *s, double *y, double dt)
{
    int nm = s->nm;
    double *x = y + 2, *vv = y + 2 + 3 * nm, ai[3 * KFMAX], ae[3 * KFMAX], db = 0.0;
    accint(s, x, ai, NULL, NULL, NULL);
    accext(s, y[0], x, ae);
    for (int c = 0; c < 3 * nm; c++) {
        double a = ai[c] + ae[c], vn = vv[c] + a * dt;
        x[c] += vv[c] * dt + 0.5 * a * dt * dt;
        db += s->m[c / 3] * 0.5 * (vv[c] + vn) * ae[c];
        vv[c] = vn;
    }
    y[1] -= dt * db;
    y[0] += dt;
}

/* the longest step in s that keeps the interpolation good (see ARCAP) */
static double hcap(const arsys *s, const double *y)
{
    int nm = s->nm;
    const double *x = y + 2;
    double u = 0.0, tmin = HUGE_VAL;
    for (int l = 0; l < nm; l++) {
        if (s->m[l] == 0.0) continue;
        for (int l2 = l + 1; l2 < nm; l2++) {
            if (s->m[l2] == 0.0) continue;
            double d[3] = {x[3 * l2] - x[3 * l], x[3 * l2 + 1] - x[3 * l + 1],
                           x[3 * l2 + 2] - x[3 * l + 2]};
            double rr = norm3(d);
            u += s->m[l] * s->m[l2] / rr;
            tmin = fmin(tmin, sqrt(rr * rr * rr / (s->m[l] + s->m[l2])));
        }
    }
    return ARCAP * u * tmin;
}

/* ---- dense output -------------------------------------------------- */
static void addnode(int k, const arsys *s, const double *y)
{
    dense_t *d = &dn[k];
    int nm = s->nm;
    if (d->nn == 0) d->nm = nm;
    if (d->nn >= d->cap) {
        d->cap = d->cap ? 2 * d->cap : 64;
        d->t = realloc(d->t, sizeof(double) * (size_t)d->cap);
        if (!d->t) die("out of memory (AR dense output)");
    }
    size_t need = (size_t)d->cap * 10 * (size_t)nm;
    if (need > d->caps) {
        d->s = realloc(d->s, sizeof(double) * need);
        if (!d->s) die("out of memory (AR dense output)");
        d->caps = need;
    }
    double ai[3 * KFMAX], ae[3 * KFMAX];
    accint(s, y + 2, ai, NULL, NULL, NULL);
    accext(s, y[0], y + 2, ae);
    double *p = d->s + (size_t)d->nn * 10 * (size_t)nm;
    for (int l = 0; l < nm; l++) {
        p[10 * l] = s->m[l];
        for (int q = 0; q < 3; q++) {
            p[10 * l + 1 + q] = y[2 + 3 * l + q];
            p[10 * l + 4 + q] = y[2 + 3 * nm + 3 * l + q];
            p[10 * l + 7 + q] = ai[3 * l + q] + ae[3 * l + q];
        }
    }
    d->t[d->nn++] = s->t0 + y[0];
}

/* masses mm, positions x and velocities u RELATIVE TO THE C.O.M. of the
 * stars of group k (bmode 2) at internal time t; returns their number.
 * From the start state if there is no dense output, otherwise
 * interpolated. Changes nothing. Massless ghosts of merged stars may be
 * among them (mass 0) until the next arshift. */
int arstate(int k, double t, double *mm, vec3 *x, vec3 *u)
{
    dense_t *d = &dn[k];
    if (d->nn < 2 || fabs(t - bepoch[k]) <= tolt()) {
        for (int l = 0; l < bnum[k]; l++) {
            mm[l] = bmass[k][l];
            memcpy(x[l], bx[k][l], sizeof(vec3));
            memcpy(u[l], bv[k][l], sizeof(vec3));
        }
        return bnum[k];
    }
    int nm = d->nm, lo = 0, hi = d->nn - 1;
    if (t <= d->t[0]) hi = 1;
    else if (t >= d->t[hi]) lo = hi - 1;
    else
        while (hi - lo > 1) {
            int mid = (lo + hi) / 2;
            if (d->t[mid] <= t) lo = mid; else hi = mid;
        }
    double dt = d->t[hi] - d->t[lo], tau = (t - d->t[lo]) / dt;
    tau = fmin(1.0, fmax(0.0, tau));
    double t2 = tau * tau, t3 = t2 * tau, t4 = t3 * tau, t5 = t4 * tau;
    double h0 = 1 - 10 * t3 + 15 * t4 - 6 * t5, h1 = tau - 6 * t3 + 8 * t4 - 3 * t5;
    double h2 = 0.5 * (t2 - 3 * t3 + 3 * t4 - t5), h3 = 10 * t3 - 15 * t4 + 6 * t5;
    double h4 = -4 * t3 + 7 * t4 - 3 * t5, h5 = 0.5 * (t3 - 2 * t4 + t5);
    double g0 = -30 * t2 + 60 * t3 - 30 * t4, g1 = 1 - 18 * t2 + 32 * t3 - 15 * t4;
    double g2 = 0.5 * (2 * tau - 9 * t2 + 12 * t3 - 5 * t4), g3 = -g0;
    double g4 = -12 * t2 + 28 * t3 - 15 * t4, g5 = 0.5 * (3 * t2 - 8 * t3 + 5 * t4);
    const double *pa = d->s + (size_t)lo * 10 * (size_t)nm, *pb = pa + 10 * (size_t)nm;
    for (int l = 0; l < nm; l++) {
        const double *a = pa + 10 * l, *b = pb + 10 * l;
        mm[l] = a[0];
        for (int q = 0; q < 3; q++) {
            x[l][q] = h0 * a[1 + q] + dt * h1 * a[4 + q] + dt * dt * h2 * a[7 + q]
                      + h3 * b[1 + q] + dt * h4 * b[4 + q] + dt * dt * h5 * b[7 + q];
            u[l][q] = (g0 * a[1 + q] + g3 * b[1 + q]) / dt + g1 * a[4 + q] + g4 * b[4 + q]
                      + dt * (g2 * a[7 + q] + g5 * b[7 + q]);
        }
    }
    return nm;
}

/* group k has been integrated to bend and its c.o.m. has now reached
 * bend: make the end state the start state (bx, bv at bepoch = bend),
 * removing the ghosts of merged stars, and drop the dense output */
void arshift(int k)
{
    dense_t *d = &dn[k];
    if (d->nn < 2) return;
    int nm = d->nm, q = 0;
    const double *p = d->s + (size_t)(d->nn - 1) * 10 * (size_t)nm;
    int id2[KFMAX], ty2[KFMAX];
    memcpy(id2, bident[k], sizeof id2);
    memcpy(ty2, btype[k], sizeof ty2);
    for (int l = 0; l < nm; l++) {
        if (p[10 * l] == 0.0) continue;
        bmass[k][q] = p[10 * l];
        bident[k][q] = id2[l];
        btype[k][q] = ty2[l];
        for (int c = 0; c < 3; c++) {
            bx[k][q][c] = p[10 * l + 1 + c];
            bv[k][q][c] = p[10 * l + 4 + c];
        }
        q++;
    }
    bnum[k] = q;
    bepoch[k] = bend[k];
    d->nn = 0;
}

/* ---- the perturbers --------------------------------------------- */
/* size of the group for gamma: apocentre of a bound pair, otherwise the
 * diameter */
static double arsize(const arsys *s, const double *y)
{
    int nm = s->nm, live[KFMAX], nl = 0;
    for (int l = 0; l < nm; l++)
        if (s->m[l] > 0.0) live[nl++] = l;
    if (nl == 2) {
        int a = live[0], b = live[1];
        double rr[3], vr[3], mb = s->m[a] + s->m[b];
        for (int q = 0; q < 3; q++) {
            rr[q] = y[2 + 3 * a + q] - y[2 + 3 * b + q];
            vr[q] = y[2 + 3 * nm + 3 * a + q] - y[2 + 3 * nm + 3 * b + q];
        }
        double rn = norm3(rr), eb = 0.5 * dot3(vr, vr) - mb / rn;
        if (eb >= 0.0) return rn;
        double h[3] = {rr[1] * vr[2] - rr[2] * vr[1], rr[2] * vr[0] - rr[0] * vr[2],
                       rr[0] * vr[1] - rr[1] * vr[0]};
        double sa = -0.5 * mb / eb, e2 = fmax(0.0, 1.0 - dot3(h, h) / (mb * sa));
        return sa * (1.0 + sqrt(e2));
    }
    double g = 0.0;
    for (int l = 0; l < nm; l++)
        if (s->m[l] > 0.0) g = fmax(g, 2.0 * norm3(y + 2 + 3 * l));
    return g;
}

/* make the perturber list of group k (c.o.m. in slot c) at internal time
 * t0, and set up s (s->nm, masses and the c.o.m. already set). Returns
 * the horizon: the end of the step of the c.o.m. or of a strong
 * perturber, whichever is first (integer time) */
static int64_t mkplist(int k, arsys *s, const double *y, double t0)
{
    int c = bcm[k], np = 0, nst = 0;
    int64_t h8 = tk[c] + dtk[c];
    double gs = arsize(s, y), xc[3], gam[NPMAX];
    int idx[NPMAX];
    tset(c, s->cc, &s->ct0);
    tpos(s->cc, t0 - s->ct0, xc);
    double fac = 2.0 * gs * gs * gs / s->mt;
    /* search all the bodies only when the c.o.m. has started a new step;
     * otherwise just look again at the ones on the list */
    int full = !(plok[k] && plcm8[k] == tk[c] && plbase[k] == tbase);
    int ncand = full ? n : npl[k];
    for (int q = 0; q < ncand; q++) {
        int j = full ? q : plist[k][q];
        if (j == c || j < 0 || j >= n) continue;
        double dt = t0 - (tbase + (double)tk[j] * tick), d2 = dt * dt / 2.0, d3 = d2 * dt / 3.0,
               d4 = d3 * dt / 4.0, d5 = d4 * dt / 5.0, dd[3];
        for (int q2 = 0; q2 < 3; q2++)
            dd[q2] = r[j][q2] + v[j][q2] * dt + f0[j][q2] * d2 + f1[j][q2] * d3
                     + f2e[j][q2] * d4 + f3e[j][q2] * d5 - xc[q2];
        double r2 = dot3(dd, dd), gj = fac * m[j] / (r2 * sqrt(r2));
        if (gj < gpert) continue;
        if (np == NPMAX && gj <= gam[np - 1]) continue;
        int p = np < NPMAX ? np++ : NPMAX - 1;      /* insertion, strongest first */
        while (p > 0 && gam[p - 1] < gj) {
            gam[p] = gam[p - 1];
            idx[p] = idx[p - 1];
            p--;
        }
        gam[p] = gj;
        idx[p] = j;
    }
    if (full) {
        plok[k] = 1;
        plcm8[k] = tk[c];
        plbase[k] = tbase;
    }
    /* the strong perturbers (first in the list) are followed exactly, and
     * set the horizon */
    for (int e = 0; e < np; e++) {
        int j = idx[e];
        plist[k][e] = j;
        if (gam[e] <= GHORFAC * gpert) continue;
        s->pm[nst] = m[j];
        tset(j, s->pc[nst], &s->pt0[nst]);
        nst++;
        if (tk[j] + dtk[j] < h8) h8 = tk[j] + dtk[j];
    }
    s->np = nst;
    npl[k] = np;
    /* the weak ones only through their summed tidal tensor, at the start
     * and end of the stretch (a weak perturber whose step ends first is
     * held still) */
    double t1 = tbase + (double)h8 * tick, xc1[3];
    tpos(s->cc, t1 - s->ct0, xc1);
    memset(s->tt, 0, sizeof s->tt);
    s->tw0 = t0;
    s->tw1 = t1;
    s->nweak = np - nst;
    for (int e = nst; e < np; e++) {
        int j = idx[e];
        double cj[6][3], tj0, p0[3], p1[3];
        tset(j, cj, &tj0);
        tpos(cj, t0 - tj0, p0);
        if (tk[j] + dtk[j] >= h8) tpos(cj, t1 - tj0, p1);
        else memcpy(p1, p0, sizeof p0);
        for (int w = 0; w < 2; w++) {
            double d[3];
            for (int q = 0; q < 3; q++) d[q] = (w ? p1[q] - xc1[q] : p0[q] - xc[q]);
            double d2 = dot3(d, d), id5 = m[j] / (d2 * d2 * sqrt(d2));
            for (int q = 0; q < 3; q++)
                for (int q2 = 0; q2 < 3; q2++)
                    s->tt[w][q][q2] += id5 * (3.0 * d[q] * d[q2] - (q == q2 ? d2 : 0.0));
        }
    }
    return h8;
}

/* set up s and y for group k from its start state at bepoch; returns the
 * horizon (see mkplist) */
static int64_t arsetup(int k, arsys *s, double *y)
{
    int nm = bnum[k];
    s->nm = nm;
    s->mt = 0.0;
    s->t0 = bepoch[k];
    y[0] = 0.0;
    for (int l = 0; l < nm; l++) {
        s->m[l] = bmass[k][l];
        s->mt += s->m[l];
        for (int q = 0; q < 3; q++) {
            y[2 + 3 * l + q] = bx[k][l][q];
            y[2 + 3 * nm + 3 * l + q] = bv[k][l][q];
        }
    }
    double ai[3 * KFMAX];
    y[1] = accint(s, y + 2, ai, NULL, NULL, NULL) - kinetic(s, y + 2 + 3 * nm);
    return mkplist(k, s, y, s->t0);
}

/* a group of bmode 2 has been made (or another kind turned into one) at
 * integer time t8, its c.o.m. in slot bcm[k] and its start state in
 * bx, bv: no dense output yet, and its first perturber list */
void arstart(int k, int64_t t8)
{
    arsys s;
    double y[NYMAX];
    bend8[k] = t8;
    bepoch[k] = bend[k] = tbase + (double)t8 * tick;
    dn[k].nn = 0;
    plok[k] = 0;
    arsetup(k, &s, y);
    pdirty = 1;
    arcapnew(k, t8);
}

/* stars ic and jc of group k have come within rcoll at the state y */
static void armerge(int k, arsys *s, double *y, int ic, int jc)
{
    int nm = s->nm;
    double m1 = s->m[ic], m2 = s->m[jc], mn = m1 + m2, rr[3], vr[3], xn[3], vn[3], xc[3],
           xa[3], va[3], vcm[3];
    double *x = y + 2, *vv = y + 2 + 3 * nm;
    for (int q = 0; q < 3; q++) {
        rr[q] = x[3 * ic + q] - x[3 * jc + q];
        vr[q] = vv[3 * ic + q] - vv[3 * jc + q];
        xn[q] = (m1 * x[3 * ic + q] + m2 * x[3 * jc + q]) / mn;
        vn[q] = (m1 * vv[3 * ic + q] + m2 * vv[3 * jc + q]) / mn;
    }
    double t = s->t0 + y[0], rsep = norm3(rr), dtc = t - s->ct0;
    tpos(s->cc, dtc, xc);
    for (int q = 0; q < 3; q++)       /* c.o.m. velocity from the Taylor series */
        vcm[q] = s->cc[1][q] + dtc * (2 * s->cc[2][q] + dtc * (3 * s->cc[3][q]
                 + dtc * (4 * s->cc[4][q] + dtc * 5 * s->cc[5][q])));
    for (int q = 0; q < 3; q++) {
        xa[q] = xc[q] + xn[q];
        va[q] = vcm[q] + vn[q];
    }
#pragma omp critical(armergecrit)
    {
        if (!mergers) {
            printf("\n******************************************************************************\n");
            printf("  STELLAR COLLISION at t = %.6f Myr (inside a few-body group)\n", t * tunit);
            printf("  star ids %12d%12d   masses %12.5f%12.5f Msun\n", bident[k][ic],
                   bident[k][jc], m1, m2);
            printf("  The run has been stopped: mergers are switched off.\n");
            printf("******************************************************************************\n");
            exit(1);
        }
        emerged += 0.5 * m1 * m2 / mn * dot3(vr, vr) - m1 * m2 / rsep;
        int idn = m2 > m1 ? bident[k][jc] : bident[k][ic];
        int tn = groupmergetype(btype[k][ic], btype[k][jc]);
        groupmergelog(t, bident[k][ic], m1, btype[k][ic], bident[k][jc], m2, btype[k][jc], idn,
                      mn, tn, xa, va, rsep);
        bident[k][ic] = idn;
        btype[k][ic] = tn;
        bmass[k][ic] = mn;
        bmass[k][jc] = 0.0;
        ntot--;
        nmerged++;
    }
    s->m[ic] = mn;
    s->m[jc] = 0.0;
    for (int q = 0; q < 3; q++) {
        x[3 * ic + q] = x[3 * jc + q] = xn[q];
        vv[3 * ic + q] = vv[3 * jc + q] = vn[q];
    }
    double ai[3 * KFMAX];
    y[1] = accint(s, x, ai, NULL, NULL, NULL) - kinetic(s, vv);
}

/* integrate group k from bepoch to its horizon (see mkplist), keeping
 * the dense output */
static void arintegrate(int k)
{
    static int warned = 0;
    arsys s;
    double y[NYMAX], yn[NYMAX], yb[NYMAX];
    int64_t h8 = arsetup(k, &s, y);
    double tgt = tbase + (double)h8 * tick;
    bend8[k] = h8;
    dn[k].nn = 0;
    addnode(k, &s, y);
    double dtot = tgt - s.t0, hs = hcap(&s, y), rm;
    int64_t nst = 0;
    int nfail = 0, ic = 0, jc = 0;
    if (dtot <= 0.0) {
        bend[k] = tgt;
        return;
    }
    for (;;) {
        double ai[3 * KFMAX], u0 = accint(&s, y + 2, ai, NULL, NULL, NULL), rem = dtot - y[0];
        hs = fmin(hs, hcap(&s, y));
        /* the last step is aimed at the end (dt/ds is about 1/U) */
        int aim = hs >= rem * u0;
        double hstep = aim ? rem * u0 : hs;
        int kc = gbs(&s, y, hstep, yn, &rm, &ic, &jc);
        nst++;
        if (kc < 0) {
            hs = 0.3 * hstep;
            if (++nfail > 100) die("the regularised integration of a group failed (step too small)");
            continue;
        }
        nfail = 0;
        int last = 0;
        if (aim || yn[0] >= dtot) {
            /* land on the end: secant iteration on t(s) until within LANDTOL
             * of it, then the rest (a tiny fraction of the step) by a Taylor
             * step in time (see LANDING above) */
            double h0 = 0.0, t0 = y[0], h1 = hstep, t1 = yn[0];
            int ok = 0;
            memcpy(yb, yn, sizeof yb);
            for (int it = 0; it < 30; it++) {
                if (fabs(t1 - dtot) <= LANDTOL * rem) {
                    ok = 1;
                    break;
                }
                double hm = h1 + (dtot - t1) * (h1 - h0) / (t1 - t0);
                if (!(hm > 0.0) || hm > 2.0 * h1) hm = (t1 > dtot ? 0.5 : 1.5) * h1;
                if (gbs(&s, y, hm, yn, &rm, &ic, &jc) < 0) break;
                nst++;
                h0 = h1;
                t0 = t1;
                h1 = hm;
                t1 = yn[0];
                memcpy(yb, yn, sizeof yb);
            }
            if (ok) {
                memcpy(yn, yb, sizeof yb);
                landfix(&s, yn, dtot - yn[0]);
                yn[0] = dtot;
                last = 1;
            } else if (t1 < dtot) {
                memcpy(yn, yb, sizeof yb);      /* not there yet: an ordinary step */
            } else {
                hs = 0.5 * hstep;               /* try again from where we are */
                continue;
            }
        }
        memcpy(y, yn, sizeof y);
        if (rcoll > 0.0 && rm < rcolli && s.m[ic] > 0.0 && s.m[jc] > 0.0) {
            armerge(k, &s, y, ic, jc);
            int live = 0;
            for (int l = 0; l < s.nm; l++)
                if (s.m[l] > 0.0) live++;
            if (live < 2 && !last) {
                /* one star left: it just moves with the c.o.m. (and the
                 * group is opened at the end of the c.o.m. step) */
                addnode(k, &s, y);
                double dt = dtot - y[0];
                for (int q = 0; q < 3 * s.nm; q++) y[2 + q] += dt * y[2 + 3 * s.nm + q];
                y[0] = dtot;
                last = 1;
            }
        }
        addnode(k, &s, y);
        if (last) break;
        if (kc <= 3) hs *= 1.6;
        else if (kc >= 6) hs *= 0.7;
        if (nst > 20000000 && !warned) {
            printf("  ** NOTE: a regularised group needs more than 2e7 steps in one go\n");
            warned = 1;
        }
    }
    bend[k] = tgt;
#pragma omp atomic
    narsteps += nst;
}

/* at the start of the block at integer time tb8: every group of bmode 2
 * that has not been integrated as far as tb8 is integrated ahead, from
 * where it got to, to its new horizon (which is at least tb8, because
 * every step ends at tb8 or later) */
void aradvance(int64_t tb8)
{
    int nk = 0;
    int *todo = xalloc(sizeof(int) * (size_t)(nbin + 1));
    for (int k = 0; k < nbin; k++) {
        if (bmode[k] != 2 || bend8[k] >= tb8) continue;
        arshift(k);
        bepoch[k] = tbase + (double)bend8[k] * tick;
        todo[nk++] = k;
    }
    if (nk > 0) {
#pragma omp parallel for schedule(dynamic, 1) if (nk > 1)
        for (int q = 0; q < nk; q++) arintegrate(todo[q]);
        pdirty = 1;
    }
    free(todo);
}

/* everyone is synchronised at internal time tsync (the end of an output
 * interval): make the end states the start states, at exactly tsync */
void arsync(double tsync)
{
    for (int k = 0; k < nbin; k++) {
        if (bmode[k] != 2) continue;
        arshift(k);
        bepoch[k] = bend[k] = tsync;
        bend8[k] = 0;
    }
}

/* ---- the corrections to the main integration ---------------------- */
static void pindex(void)
{
    int need = 0;
    if (phalloc < nalloc) {
        free(phead);
        free(pmarked);
        phead = xalloc(sizeof(int) * (size_t)nalloc);
        pmarked = xalloc(sizeof(int) * (size_t)nalloc);
        for (int j = 0; j < nalloc; j++) phead[j] = -1;
        phalloc = nalloc;
        nmarked = 0;
    }
    for (int q = 0; q < nmarked; q++) phead[pmarked[q]] = -1;
    nmarked = 0;
    for (int k = 0; k < nbin; k++)
        if (bmode[k] == 2) need += npl[k];
    if (need > npool) {
        npool = need + 256;
        free(pk);
        free(pnext);
        pk = xalloc(sizeof(int) * (size_t)npool);
        pnext = xalloc(sizeof(int) * (size_t)npool);
    }
    int e = 0;
    for (int k = 0; k < nbin; k++) {
        if (bmode[k] != 2) continue;
        for (int q = 0; q < npl[k]; q++) {
            int j = plist[k][q];
            if (j < 0 || j >= n) continue;
            if (phead[j] < 0) pmarked[nmarked++] = j;
            pk[e] = k;
            pnext[e] = phead[j];
            phead[j] = e++;
        }
    }
    pdirty = 0;
}

/* g(d) = m d/|d|^3 and its time derivative, added to a and jk */
static inline void gpair(double mj, const double d[3], const double dv[3], double sgn,
                         double a[3], double jk[3])
{
    double r2 = dot3(d, d), ir2 = 1.0 / r2, mor = sgn * mj * ir2 * sqrt(ir2);
    double rv = 3.0 * dot3(d, dv) * ir2;
    for (int q = 0; q < 3; q++) {
        a[q] += mor * d[q];
        jk[q] += mor * (dv[q] - rv * d[q]);
    }
}

/* internal dynamical time sqrt(s^3/M) of a group of nm stars, s its
 * diameter */
static double tdyngrp(int nm, const double *mm, vec3 *x)
{
    double s = 0.0, mt = 0.0;
    for (int l = 0; l < nm; l++) {
        if (mm[l] == 0.0) continue;
        s = fmax(s, 2.0 * norm3(x[l]));
        mt += mm[l];
    }
    return sqrt(s * s * s / mt);
}

/* pairwise free-fall and flyby times squared of a star of mass ms at
 * separation d, relative velocity dv, from a body of mass mi; and its
 * distance squared */
static inline void pairtime(double mi, double ms, const double d[3], const double dv[3],
                            double *tp2, double *r2m)
{
    double r2 = dot3(d, d), v2 = dot3(dv, dv);
    double tff = r2 * sqrt(r2) / (mi + ms), tfl = r2 / (v2 + 1e-300);
    if (tff < *tp2) *tp2 = tff;
    if (tfl < *tp2) *tp2 = tfl;
    if (r2 < *r2m) *r2m = r2;
}

/* the correction (see the top of this file) to the force a and jerk jk
 * of body i at internal time tb, from the predicted rp, vp: as the c.o.m.
 * of a regularised group, and as a perturber of others; also the shortest
 * dynamical time of those groups (tdg), and the shortest pairwise time
 * squared (tp2) and distance squared (r2m) to their stars. Returns 0 if
 * there is none. The perturber index must be up to date. */
static int corrfor(int i, double tb, double a[3], double jk[3], double *tdg, double *tp2,
                   double *r2m)
{
    int k = ib[i], nm, has = 0;
    double mm[KFMAX];
    vec3 x[KFMAX], u[KFMAX];
    a[0] = a[1] = a[2] = jk[0] = jk[1] = jk[2] = 0.0;
    *tdg = *tp2 = *r2m = HUGE_VAL;
    if (k >= 0 && bmode[k] == 2 && npl[k] > 0) {
        nm = arstate(k, tb, mm, x, u);
        *tdg = fmin(*tdg, tdyngrp(nm, mm, x));
        for (int e = 0; e < npl[k]; e++) {
            int j = plist[k][e];
            double d[3], dv[3];
            for (int l = 0; l < nm; l++) {
                if (mm[l] == 0.0) continue;
                for (int q = 0; q < 3; q++) {
                    d[q] = rp[j][q] - rp[i][q] - x[l][q];
                    dv[q] = vp[j][q] - vp[i][q] - u[l][q];
                }
                gpair(m[j] * mm[l] / m[i], d, dv, 1.0, a, jk);
                pairtime(mm[l], m[j], d, dv, tp2, r2m);
            }
            for (int q = 0; q < 3; q++) {
                d[q] = rp[j][q] - rp[i][q];
                dv[q] = vp[j][q] - vp[i][q];
            }
            gpair(m[j], d, dv, -1.0, a, jk);
        }
        has = 1;
    }
    if (i < phalloc)
        for (int e = phead[i]; e >= 0; e = pnext[e]) {
            int kg = pk[e], c = bcm[kg];
            double d[3], dv[3];
            nm = arstate(kg, tb, mm, x, u);
            *tdg = fmin(*tdg, tdyngrp(nm, mm, x));
            for (int l = 0; l < nm; l++) {
                if (mm[l] == 0.0) continue;
                for (int q = 0; q < 3; q++) {
                    d[q] = rp[c][q] + x[l][q] - rp[i][q];
                    dv[q] = vp[c][q] + u[l][q] - vp[i][q];
                }
                gpair(mm[l], d, dv, 1.0, a, jk);
                pairtime(m[i], mm[l], d, dv, tp2, r2m);
            }
            for (int q = 0; q < 3; q++) {
                d[q] = rp[c][q] - rp[i][q];
                dv[q] = vp[c][q] - vp[i][q];
            }
            gpair(m[c], d, dv, -1.0, a, jk);
            has = 1;
        }
    return has;
}

/* the longest step a correction a, jk allows a body whose force is f;
 * tdg the shortest dynamical time of the groups it comes from (see
 * TIMESTEPS above) */
static double capfor(const double a[3], const double jk[3], const double f[3], double tdg)
{
    double da = norm3(a), dj = norm3(jk);
    if (da <= CAPFRAC * norm3(f)) return HUGE_VAL;
    double cap = etap * tdg;
    if (dj > 0.0) cap = fmin(cap, sqrt(eta) * da / dj);
    return cap;
}

/* the corrections to the force and jerk of the active bodies, at the
 * block time, from the predicted rp, vp */
void arcorrect(void)
{
    int any = 0;
    for (int k = 0; k < nbin && !any; k++)
        if (bmode[k] == 2) any = 1;
    if (!any) return;
    if (pdirty) pindex();
    double tb = tbase + (double)tcur8 * tick;
#pragma omp parallel for schedule(dynamic, 8) if (nact >= nompmin)
    for (int ia = 0; ia < nact; ia++) {
        int i = ilist[ia];
        double a[3], jk[3], tdg, tp2, r2m;
        if (!corrfor(i, tb, a, jk, &tdg, &tp2, &r2m)) continue;
        /* the pairwise step limit (etap) and closest approach measured
         * against the group's stars, not just its c.o.m. (see TIMESTEPS) */
        if (tp2 < tpair2[i]) tpair2[i] = tp2;
        if (r2m < rmin2[i]) rmin2[i] = r2m;
        for (int q = 0; q < 3; q++) {
            f0p[i][q] += a[q];
            f1p[i][q] += jk[q];
            if (useac) {
                fin[i][q] += a[q];
                fidn[i][q] += jk[q];
            }
        }
    }
}

/* group k has just been made at integer time t8: the bodies that feel it
 * and have just started a step (its c.o.m. and its perturbers that are
 * at t8) took that step without knowing of it, so it is shortened now if
 * the cap asks for it (a shorter step from t8 always fits the block
 * clock). Without this, the first step of a star near a binary made at
 * the start can be half the binary's orbit (see TIMESTEPS above). */
static void arcapnew(int k, int64_t t8)
{
    double tb = tbase + (double)t8 * tick, a[3], jk[3], tdg, tp2, r2m;
    int c = bcm[k];
    if (pdirty) pindex();
    for (int e = -1; e < npl[k]; e++) {
        int j = e < 0 ? c : plist[k][e];
        if (tk[j] != t8 || !corrfor(j, tb, a, jk, &tdg, &tp2, &r2m)) continue;
        double cap = fmin(capfor(a, jk, f0[j], tdg), etap * sqrt(tp2));
        if ((double)dtk[j] * tick <= cap) continue;
        int lv = commlev(getlevel(cap), t8);
        if (lv > lev[j]) {
            lev[j] = lv;
            dtk[j] = (int64_t)1 << (kmax - lv);
            if (useac) aclevel(j);
        }
    }
}

/* the energy of the resolved interaction of each regularised group with
 * its perturbers less that of the point mass the main integration uses
 * (see the top of this file), at internal time t, from r (pred = 0, all
 * bodies synchronised) or rp (pred = 1) */
double arpotcorr(double t, int pred)
{
    double e = 0.0, mm[KFMAX];
    vec3 x[KFMAX], u[KFMAX];
    vec3 *pos = pred ? rp : r;
    for (int k = 0; k < nbin; k++) {
        if (bmode[k] != 2 || npl[k] == 0) continue;
        int c = bcm[k], nm = arstate(k, t, mm, x, u);
        for (int e2 = 0; e2 < npl[k]; e2++) {
            int j = plist[k][e2];
            if (j < 0 || j >= n) continue;
            double d[3];
            for (int q = 0; q < 3; q++) d[q] = pos[c][q] - pos[j][q];
            e += m[c] * m[j] / norm3(d);
            for (int l = 0; l < nm; l++) {
                if (mm[l] == 0.0) continue;
                for (int q = 0; q < 3; q++) d[q] = pos[c][q] + x[l][q] - pos[j][q];
                e -= mm[l] * m[j] / norm3(d);
            }
        }
    }
    return e;
}
