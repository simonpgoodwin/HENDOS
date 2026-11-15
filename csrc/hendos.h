/* ==============================================================
 * HENDOS - HErmite N-body Dynamics Of Systems
 * THE C VERSION of the main simulation code (for comparison with the
 * Fortran version in ../src). It reads the same parameter file and
 * initial conditions and writes the same output files, so makeics, the
 * analysis tools, the regression tests and the benchmark all work
 * with either. See CSRC_README.md.
 *
 * hendos.h - everything shared between the C files: the parameters,
 * the particle arrays, the groups, the neighbour scheme and the
 * routines. (The Fortran equivalent is the module nbvars.)
 *
 * Differences from the Fortran to keep in mind when reading both:
 *   - arrays start at 0, so body i is 0..n-1, a group index is
 *     0..nbin-1 and 'no group' / 'no body' is -1 (0 in the Fortran)
 *   - a vector is double[3]; r[i][0] is x of body i (r(1,i) in Fortran)
 * Internal units as in the Fortran: G = 1, Msun, pc.
 * ============================================================== */
#ifndef HENDOS_H
#define HENDOS_H
#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#ifdef HENDOS_MAIN
#define EXT
#else
#define EXT extern
#endif

#define KFMAX 6            /* most stars in a group */
#define NTYPES 5           /* particle types */
#define NSTR 512
typedef double vec3[3];

/* ---- particle types: code 1..NTYPES, 0 = the c.o.m. body of a group */
enum { TSTAR = 1, TWDWARF, TNSTAR, TBHOLE, TPLANET };
extern const char *typename_[NTYPES + 1];

/* ---- parameters (see run/PARAMS_README.md) ---------------------- */
EXT char runid[64], icfile[NSTR], datadir[NSTR], screen[32], tidal[32];
EXT char pfilename[NSTR];
EXT double tend, tout, eps, eps2, tsnap, dewarn, trestart, tana;
EXT int lagrangian, lverbose, restart_;
EXT double oorta, oortb, rho0, tij[3][3], omegat;
EXT double ttid[3][3], omega;          /* tide in internal units */
EXT int ltide;
EXT double eta, etap, etas;
EXT int niter, kmax;
EXT int closebinaries, fewbody;
EXT double dtclose, dtclosei, gclose, gopen, gfew, etafew;
EXT double rcoll, rcolli;
EXT int mergers;
EXT int nensemble, nparallel, ensfirst, ensmember;
EXT char ensbase[64];
EXT int useac, nnbopt, nnbmax;
EXT double etar;
EXT int usetree, nleaf, ntreebuild;    /* read, but the tree is not in the C version */
EXT double theta;
EXT int nompmin, nompj, nthreads;

/* ---- units and the block-timestep clock ------------------------- */
EXT double vunit, tunit;               /* km/s and Myr per internal unit */
EXT double dtmax, tick;                /* output interval, one tick (internal) */
EXT int64_t tend8;                     /* ticks in one output interval */

/* ---- the bodies ------------------------------------------------- */
EXT int n, ntot, nalloc, nbin, nact;   /* bodies, stars, array size, groups, active */
EXT double time0, comshift[3], covshift[3];
EXT char icformat[64];
EXT vec3 *r, *v, *rp, *vp, *rc, *vc, *f0, *f1, *f0p, *f1p, *f2, *f3, *f2e, *f3e;
EXT double *m, *tpair2, *rmin2, *phi;
EXT int *id, *itype, *ib, *lev, *ilist;
EXT int64_t *tk, *dtk;

/* ---- groups: stars treated as one c.o.m. body (see groups.c) ---- */
EXT int *bnum, *bmode;                 /* stars in the group; 1 binary, 2 few-body, 3 triple */
EXT double (*bmass)[KFMAX], (*brel)[6], (*brel2)[6], *bepoch;
EXT double *bein;                     /* frozen triple: its internal energy when frozen */
EXT int (*bident)[KFMAX], (*btype)[KFMAX];
EXT vec3 (*bx)[KFMAX], (*bv)[KFMAX];

/* ---- Ahmad-Cohen neighbour scheme (see ac.c) --------------------- */
EXT vec3 *fi, *fid, *fr, *frd, *fin, *fidn, *frn, *frdn, *dreg, *dvreg;
EXT double *rs, *dtwr, rs0;
EXT int *nnb, *levr, *nbl;             /* list of body i is nbl[i*nnbmax + k] */
EXT int64_t *treg, *dtreg, tcur8;
EXT char *lreg, *acok;
EXT int aclast;

/* ---- counters ---------------------------------------------------- */
EXT int64_t nclamp, nclosed, nopened, nfewformed, nfewopened, nfrozen;
EXT int64_t nfewsteps, nmerged, nregular, nirregular;
EXT double emerged;                   /* orbital energy removed by mergers */
EXT double efrozen;                   /* energy change of the triples unfrozen so far */

/* ---- diagnostics (set by diagnose) ------------------------------- */
EXT double ekin, epot, etid, etot, etot0, eprev, demax;
EXT double cm[3], cv[3], sig3d, rlag[3];
EXT int64_t nmergedsnap;
EXT int ioutstart;

/* ---- routines ---------------------------------------------------- */
/* params.c */
void setunits(void);
void readparams(const char *pfile);
void ensname(const char *name, int k, char *out);
/* io.c */
void checkdatadir(void);
void allocarrays(void);
void readic(void);
void diagnose(double time);
void lagradii(int nf, const double *frac, double *rl);
void startinfo(int nout);
void statusout(int iout, int nout, double time, int64_t nblock, int64_t nforce);
void endinfo(double time, int64_t nforce);
void snapname(double time, char *fname);
void snapshot(double time, int64_t nforce);
void hdrcommon(FILE *f);
void ana_init(void);
void ana_do(double time);
void ana_end(void);
void checkrunid(int nsnap, int nana);
double walltime(void);
void startclock(void);
void restartcounters(int64_t nblock, int64_t nforce);
/* force.c */
void forces(void);
void tideadd(int i, double a[3], double jk[3]);
void fullsum(int i, int par, double a[3], double jk[3], double *tmin2, double *rm2);
void energy(double *ke, double *pot, double *ptid);
void firststep(void);
void newlevel(int i);
int getlevel(double dtw);
int commlev(int k, int64_t t8);
/* ac.c */
void accheck(int nstars);
void acalloc(void);
void acsetrs0(void);
void acforce(int i, int inpar);
void acstart(int i);
void acaccept(int i);
void aclevel(int i);
void acmoverow(int src, int dst);
void acmoved(int src, int dst);
void acremoved(int j, int irep);
void acadded(int jn, int i0);
/* few.c */
void fewint(int nm, const double *mm, vec3 *x, vec3 *vv, double tint, double etaf,
            double e2, double rc2, double *tdone, int *ic, int *jc, int64_t *nst);
void elements(double mu, const double rr[3], const double vr[3], double *a, double *ecc,
              double h[3]);
void hierarchy(const double mm[3], vec3 x[3], vec3 vv[3], int *stable, int *i1, int *i2,
               int *i3, double rel[6], double rel2[6], double *qin, double *qout);
/* groups.c */
void kepler(double mu, const double r0[3], const double v0[3], double dt, double r1[3],
            double v1[3]);
void grpcomp(int k, const double xc[3], const double vc_[3], double tnow, int *nm,
             vec3 xm[KFMAX], vec3 um[KFMAX]);
double grpeint(int nm, const double *mm, vec3 xm[KFMAX], vec3 um[KFMAX]);
void collcheck(int64_t t8, double tnow);
void bincheck(int64_t t8, double tnow);
void binstart(double tnow);
/* restart.c */
void writerestart(double time, int iout, int64_t nblock, int64_t nforce);
void readrestart(double *time, int *iout, int64_t *nblock, int64_t *nforce);

/* small helpers */
static inline double dot3(const double a[3], const double b[3])
{
    return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
}
static inline double norm3(const double a[3]) { return sqrt(dot3(a, a)); }
void die(const char *msg);
void *xalloc(size_t nbytes);
#endif
