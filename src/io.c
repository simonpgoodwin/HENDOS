/* ==============================================================
 * HENDOS - HErmite N-body Dynamics Of Systems - io.c
 * INPUT AND OUTPUT: the data directory, the initial conditions, the
 * screen, diagnostics, snapshots and the in-run analysis file.
 * The formats are described in README.md. They are the same as those
 * of the archived Fortran version (nbsetup.f90, nbout.f90, nbana.f90),
 * so every tool reads files from either; the creator line in each
 * header says which code wrote the file.
 * ============================================================== */
#include <omp.h>
#include <sys/stat.h>
#include <time.h>
#include "hendos.h"

static double clock0, wlast;
static int64_t nblock0, nforce0, nclamp0, nclosed0, nopened0, nmerged0, nfewformed0,
    nfewopened0, nregular0, nirregular0;
static FILE *flagr = NULL;
static double efrozennow = 0.0;   /* efrozen + what unfreezing every triple now would add */
static const char *DLINE =
    " ==============================================================================";
static const char *SLINE =
    " ------------------------------------------------------------------------------";

double walltime(void) { return omp_get_wtime() - clock0; }
void startclock(void)
{
    clock0 = omp_get_wtime();
    wlast = 0.0;
}

/* after a restart: the counters the output intervals are measured from */
void restartcounters(int64_t nblock, int64_t nforce)
{
    nblock0 = nblock; nforce0 = nforce; nclamp0 = nclamp; nclosed0 = nclosed;
    nopened0 = nopened; nmerged0 = nmerged; nfewformed0 = nfewformed;
    nfewopened0 = nfewopened; nregular0 = nregular; nirregular0 = nirregular;
}

static const char *hms(double sec, char *buf)
{
    long s = (long)(sec + 0.5);
    sprintf(buf, "%6ld:%02ld:%02ld", s / 3600, (s / 60) % 60, s % 60);
    return buf;
}

/* make the data directory if it does not exist, and stop now if files
 * cannot be written there */
void checkdatadir(void)
{
    char cmd[2 * NSTR], test[2 * NSTR];
    if (!datadir[0]) return;
    snprintf(cmd, sizeof cmd, "mkdir -p %s", datadir);
    if (system(cmd) != 0) { /* the write test below reports the problem */ }
    snprintf(test, sizeof test, "%s.hendos_write_test", datadir);
    FILE *f = fopen(test, "w");
    if (!f) {
        printf(" ERROR: cannot write to the data directory %s\n", datadir);
        printf("        Check datadir in the &output section of the parameter file.\n");
        exit(1);
    }
    fclose(f);
    remove(test);
}

void allocarrays(void)
{
    size_t nn = (size_t)n, nv = sizeof(vec3) * nn, ng = nn / 2 + 1;
    r = xalloc(nv); v = xalloc(nv); rp = xalloc(nv); vp = xalloc(nv);
    rc = xalloc(nv); vc = xalloc(nv);
    f0 = xalloc(nv); f1 = xalloc(nv); f0p = xalloc(nv); f1p = xalloc(nv);
    f2 = xalloc(nv); f3 = xalloc(nv); f2e = xalloc(nv); f3e = xalloc(nv);
    m = xalloc(sizeof(double) * nn); tpair2 = xalloc(sizeof(double) * nn);
    rmin2 = xalloc(sizeof(double) * nn); phi = xalloc(sizeof(double) * nn);
    id = xalloc(sizeof(int) * nn); itype = xalloc(sizeof(int) * nn);
    ib = xalloc(sizeof(int) * nn); lev = xalloc(sizeof(int) * nn);
    ilist = xalloc(sizeof(int) * nn);
    tk = xalloc(sizeof(int64_t) * nn); dtk = xalloc(sizeof(int64_t) * nn);
    bnum = xalloc(sizeof(int) * ng); bmode = xalloc(sizeof(int) * ng);
    bmass = xalloc(sizeof(*bmass) * ng); bident = xalloc(sizeof(*bident) * ng);
    btype = xalloc(sizeof(*btype) * ng); brel = xalloc(sizeof(*brel) * ng);
    brel2 = xalloc(sizeof(*brel2) * ng); bepoch = xalloc(sizeof(double) * ng);
    bein = xalloc(sizeof(double) * ng);
    bx = xalloc(sizeof(*bx) * ng); bv = xalloc(sizeof(*bv) * ng);
    nalloc = n;
    for (int i = 0; i < n; i++) {
        ib[i] = -1;
        rmin2[i] = HUGE_VAL;
    }
}

static int typecode(const char *name)
{
    for (int k = 1; k <= NTYPES; k++)
        if (!strcasecmp(name, typename_[k])) return k;
    return 0;
}

/* read the initial conditions file icfile. Two formats:
 * NBODY-SNAPSHOT (written by makeics and as snapshots): header lines
 *   '# key = value' ending with '# end_header', then one line per star
 *     version 2:  id m x y z vx vy vz type [anything else is ignored]
 *     version 1:  id m x y z vx vy vz      [anything else is ignored]
 *   the start time is 'time_Myr' in the header
 * legacy: lines  m x y z vx vy vz, optional # lines at the top
 * Units Msun, pc, km/s. Velocities are converted to internal units and
 * everything is moved to the centre-of-mass frame. */
void readic(void)
{
    FILE *f = fopen(icfile, "r");
    char line[2048], key[128], tname[64];
    int snapfmt = 0, ifmt = 1, nhdr = -1, nhead = 0;
    if (!f) {
        printf(" ERROR: cannot open initial conditions file %s\n", icfile);
        exit(1);
    }
    time0 = 0.0;
    long datastart = 0;
    while (fgets(line, sizeof line, f)) {
        char *p = line;
        while (*p == ' ' || *p == '\t') p++;
        if (*p != '#') break;
        nhead++;
        datastart = ftell(f);
        char *eq = strchr(p, '=');
        if (!eq) continue;
        int len = (int)(eq - (p + 1));
        snprintf(key, sizeof key, "%.*s", len, p + 1);
        char *k = key;
        while (*k == ' ') k++;
        for (char *e = k + strlen(k) - 1; e >= k && (*e == ' ' || *e == '\t'); e--) *e = 0;
        if (!strcmp(k, "format")) {
            char *s = strstr(eq + 1, "NBODY-SNAPSHOT");
            if (s) {
                snapfmt = 1;
                if (sscanf(s + 14, "%d", &ifmt) != 1) ifmt = 1;
            }
        } else if (!strcmp(k, "time_Myr")) {
            sscanf(eq + 1, "%lf", &time0);
        } else if (!strcmp(k, "n")) {
            sscanf(eq + 1, "%d", &nhdr);
        }
    }
    if (snapfmt) snprintf(icformat, sizeof icformat, "NBODY-SNAPSHOT %d", ifmt);
    else {
        strcpy(icformat, "legacy (m x y z vx vy vz)");
        time0 = 0.0;
    }
    /* count the stars: lines that are not blank */
    fseek(f, datastart, SEEK_SET);
    n = 0;
    while (fgets(line, sizeof line, f))
        if (strspn(line, " \t\r\n") != strlen(line)) n++;
    if (n < 2) {
        printf(" ERROR: fewer than 2 stars read from %s\n", icfile);
        exit(1);
    }
    if (nhdr >= 0 && nhdr != n)
        printf(" WARNING: header of initial conditions says n = %d but %d stars were read\n",
               nhdr, n);
    allocarrays();
    fseek(f, datastart, SEEK_SET);
    for (int i = 0; i < n;) {
        if (!fgets(line, sizeof line, f)) break;
        if (strspn(line, " \t\r\n") == strlen(line)) continue;
        int got;
        if (snapfmt && ifmt >= 2) {
            got = sscanf(line, "%d %lf %lf %lf %lf %lf %lf %lf %63s", &id[i], &m[i], &r[i][0],
                         &r[i][1], &r[i][2], &v[i][0], &v[i][1], &v[i][2], tname);
            if (got != 9) {
                printf(" ERROR: cannot read particle %d of %s - each line needs  "
                       "id m x y z vx vy vz type\n", i + 1, icfile);
                exit(1);
            }
            itype[i] = typecode(tname);
            if (itype[i] == 0) {
                printf(" ERROR: particle %d of %s has unknown type '%s'\n", i + 1, icfile, tname);
                printf("        known types:");
                for (int k = 1; k <= NTYPES; k++) printf(" %s", typename_[k]);
                printf("\n");
                exit(1);
            }
        } else if (snapfmt) {
            got = sscanf(line, "%d %lf %lf %lf %lf %lf %lf %lf", &id[i], &m[i], &r[i][0],
                         &r[i][1], &r[i][2], &v[i][0], &v[i][1], &v[i][2]);
            if (got != 8) die("cannot read a particle of the initial conditions");
            itype[i] = TSTAR;
        } else {
            got = sscanf(line, "%lf %lf %lf %lf %lf %lf %lf", &m[i], &r[i][0], &r[i][1],
                         &r[i][2], &v[i][0], &v[i][1], &v[i][2]);
            if (got != 7) die("cannot read a particle of the initial conditions");
            id[i] = i + 1;
            itype[i] = TSTAR;
        }
        i++;
    }
    fclose(f);
    /* km/s to internal units, then to the centre-of-mass frame */
    double com[3] = {0, 0, 0}, cov[3] = {0, 0, 0}, mtot = 0.0;
    for (int i = 0; i < n; i++) {
        mtot += m[i];
        for (int k = 0; k < 3; k++) {
            v[i][k] /= vunit;
            com[k] += m[i] * r[i][k];
            cov[k] += m[i] * v[i][k];
        }
    }
    for (int k = 0; k < 3; k++) {
        com[k] /= mtot;
        cov[k] /= mtot;
        comshift[k] = com[k];
        covshift[k] = cov[k] * vunit;
    }
    for (int i = 0; i < n; i++)
        for (int k = 0; k < 3; k++) {
            r[i][k] -= com[k];
            v[i][k] -= cov[k];
        }
    ntot = n;
    nbin = 0;
}

/* index that sorts a[0..nn-1] into increasing order */
static const double *sortkey;
static int cmpidx(const void *a, const void *b)
{
    double x = sortkey[*(const int *)a], y = sortkey[*(const int *)b];
    return (x > y) - (x < y);
}
static void indexsort(int nn, const double *a, int *idx)
{
    for (int i = 0; i < nn; i++) idx[i] = i;
    sortkey = a;
    qsort(idx, (size_t)nn, sizeof(int), cmpidx);
}

/* Lagrangian radii (pc) about the centre of mass cm for the mass
 * fractions frac[0..nf-1], which must be in increasing order */
void lagradii(int nf, const double *frac, double *rl)
{
    double *r2 = xalloc(sizeof(double) * (size_t)n), mt = 0.0, mc = 0.0;
    int *idx = xalloc(sizeof(int) * (size_t)n), k = 0;
    for (int i = 0; i < n; i++) {
        double d[3] = {r[i][0] - cm[0], r[i][1] - cm[1], r[i][2] - cm[2]};
        r2[i] = dot3(d, d);
        mt += m[i];
    }
    indexsort(n, r2, idx);
    for (int q = 0; q < nf; q++) rl[q] = sqrt(r2[idx[n - 1]]);
    for (int i = 0; i < n && k < nf; i++) {
        mc += m[idx[i]];
        while (k < nf && mc >= frac[k] * mt) rl[k++] = sqrt(r2[idx[i]]);
    }
    free(r2);
    free(idx);
}

/* energies (and phi), centre of mass, velocity dispersion and
 * 10/50/90% Lagrangian radii. The internal energy of the groups is
 * added at the current phase of their motion, and emerged (the orbital
 * energy removed by mergers) so that the total is still conserved.
 * FROZEN TRIPLES: two nested Kepler orbits keep the energies of the
 * inner and outer orbits fixed, but leave out the (small) interaction
 * of the outer star with the two inner stars separately, so the true
 * energy of a frozen triple wanders a little with the phases of its
 * orbits - and in a cluster rich in hard binaries that is a large
 * part of the total. So the total counts each frozen triple at its
 * energy when it was frozen (bein), and when a triple is unfrozen the
 * difference from its true energy then is added to efrozen. The total
 * is then conserved to the accuracy of the main integration, and
 * efrozen (header key energy_change_frozen_triples) shows how much the
 * freezing has changed the energy. */
void diagnose(double time)
{
    static const double frac[3] = {0.1, 0.5, 0.9};
    double ecorr = 0.0;
    energy(&ekin, &epot, &etid);
    for (int i = 0; i < n; i++) {
        int k = ib[i], nm;
        vec3 xm[KFMAX], um[KFMAX];
        if (k < 0) continue;
        grpcomp(k, r[i], v[i], time / tunit, &nm, xm, um);
        if (bmode[k] == 3) ecorr += bein[k] - grpeint(nm, bmass[k], xm, um);
        for (int l = 0; l < nm; l++) {
            double u[3] = {um[l][0] - v[i][0], um[l][1] - v[i][1], um[l][2] - v[i][2]};
            ekin += 0.5 * bmass[k][l] * dot3(u, u);
            for (int l2 = l + 1; l2 < nm; l2++) {
                double d[3] = {xm[l][0] - xm[l2][0], xm[l][1] - xm[l2][1], xm[l][2] - xm[l2][2]};
                epot -= bmass[k][l] * bmass[k][l2] / norm3(d);
            }
        }
    }
    efrozennow = efrozen + ecorr;
    etot = ekin + epot + etid + emerged + efrozennow;
    double mt = 0.0;
    for (int k = 0; k < 3; k++) cm[k] = cv[k] = 0.0;
    for (int i = 0; i < n; i++) {
        mt += m[i];
        for (int k = 0; k < 3; k++) {
            cm[k] += m[i] * r[i][k];
            cv[k] += m[i] * v[i][k];
        }
    }
    for (int k = 0; k < 3; k++) {
        cm[k] /= mt;
        cv[k] /= mt;
    }
    sig3d = 0.0;
    for (int i = 0; i < n; i++) {
        double u[3] = {v[i][0] - cv[0], v[i][1] - cv[1], v[i][2] - cv[2]};
        sig3d += m[i] * dot3(u, u);
    }
    sig3d = sqrt(sig3d / mt);
    lagradii(3, frac, rlag);
}

/* text giving the number of each particle type, e.g. 'star 998 black_hole 2' */
static void typecounts(char *out)
{
    int cnt[NTYPES + 1] = {0};
    for (int i = 0; i < n; i++) {
        if (ib[i] < 0) cnt[itype[i]]++;
        else
            for (int l = 0; l < bnum[ib[i]]; l++) cnt[btype[ib[i]][l]]++;
    }
    out[0] = 0;
    for (int k = 1; k <= NTYPES; k++)
        if (cnt[k] > 0) sprintf(out + strlen(out), "%s%s %d", out[0] ? " " : "", typename_[k], cnt[k]);
}

/* table of the occupied timestep levels (verbose screen output) */
static void writelevels(void)
{
    int hist[64] = {0};
    for (int i = 0; i < n; i++) hist[lev[i]]++;
    printf("   timestep levels       level      step (yr)       stars\n");
    for (int k = 0; k <= kmax; k++)
        if (hist[k] > 0) printf("%22s%6d%15.4E%12d\n", "", k, tout * 1.e6 / ldexp(1.0, k), hist[k]);
}

/* screen output at the start of the run; sets the energy the error is
 * measured from */
void startinfo(int nout)
{
    char tc[256];
    diagnose(time0);
    etot0 = etot;
    eprev = etot;
    demax = 0.0;
    startclock();
    restartcounters(0, 0);
    typecounts(tc);
    printf("%s\n", DLINE);
    printf("  HENDOS - HErmite N-body Dynamics Of Systems     run: %s\n", runid);
    printf("%s\n", DLINE);
    printf("  initial conditions     %s  (%s)\n", icfile, icformat);
    printf("  particles              %14d\n", ntot);
    if (ensmember > 0)
        printf("  ensemble member        %14d of %d of ensemble %s\n", ensmember, nensemble, ensbase);
    printf("  OpenMP threads         %14d\n", nthreads);
    printf("  particle types         %s\n", tc);
    printf("  start, end time        %14.4f%14.4f Myr\n", time0, tend);
    printf("  output every           %14.4f Myr  (%d outputs)\n", tout, nout);
    if (usetree) printf("  forces                 %14.2f  (tree opening angle; energies direct)\n", theta);
    else if (useac) printf("  forces                 %14d  (direct, neighbour scheme: neighbours wanted)\n", nnbopt);
    else printf("  forces                        direct\n");
    printf("  tide                   %s\n", tidal);
    if (closebinaries)
        printf("  groups closed at start %14d  (binaries %lld, few-body %lld)\n", nbin,
               (long long)nclosed, (long long)nfewformed);
    printf("  total energy           %14.6E Msun km2/s2\n", etot * vunit * vunit);
    printf("  virial ratio           %14.4f\n", -ekin / epot);
    printf("  half-mass radius       %14.4f pc\n", rlag[1]);
    printf("  velocity dispersion    %14.4f km/s\n", sig3d * vunit);
    if (lverbose) {
        double mt = 0.0, mmin = HUGE_VAL, mmax = 0.0, v2 = vunit * vunit;
        for (int i = 0; i < n; i++) {
            mt += m[i];
            int k = ib[i];
            for (int l = 0; l < (k < 0 ? 1 : bnum[k]); l++) {
                double ml = k < 0 ? m[i] : bmass[k][l];
                mmin = fmin(mmin, ml);
                mmax = fmax(mmax, ml);
            }
        }
        printf("%s\n", SLINE);
        printf("  PARAMETERS from %s\n", pfilename);
        printf("    eta %g   etap %g   etas %g   niter %d   kmax %d (shortest step %.3E yr)\n",
               eta, etap, etas, niter, kmax, tout * 1.e6 / ldexp(1.0, kmax));
        printf("    eps %g pc   nompmin %d   nompj %d   dewarn %g\n", eps, nompmin, nompj, dewarn);
        printf("    closebinaries %d   dtclose %g yr   gclose %g   gopen %g\n", closebinaries,
               dtclose, gclose, gopen);
        printf("    fewbody %d   gfew %g   etafew %g\n", fewbody, gfew, etafew);
        printf("    rcoll %g Rsun   mergers %d\n", rcoll, mergers);
        if (usetree) printf("    usetree 1   theta %g   nleaf %d   ntreebuild %d   etatree %g\n",
                            theta, nleaf, ntreebuild, etatree);
        if (useac) printf("    useac 1   nnbopt %d   nnbmax %d   etar %g\n", nnbopt, nnbmax, etar);
        else printf("    useac 0\n");
        printf("    trestart %g Myr (0 = no restart files)\n", trestart);
        if (ltide) {
            double fconv = 1.0 / (vunit * 1.e3);
            printf("    Omega %.4f km/s/kpc   T diagonal %.2f %.2f %.2f (km/s/kpc)^2\n",
                   omega / fconv, ttid[0][0] / (fconv * fconv), ttid[1][1] / (fconv * fconv),
                   ttid[2][2] / (fconv * fconv));
        }
        printf("%s\n", SLINE);
        printf("  INITIAL CONDITIONS\n");
        printf("    c.o.m. shift            %12.4E%12.4E%12.4E pc\n", comshift[0], comshift[1],
               comshift[2]);
        printf("    c.o.m. velocity shift   %12.4E%12.4E%12.4E km/s\n", covshift[0], covshift[1],
               covshift[2]);
        printf("    mass min, mean, max     %12.4E%12.4E%12.4E Msun\n", mmin, mt / ntot, mmax);
        printf("    total mass              %12.4f Msun\n", mt);
        printf("    Lagrangian 10,50,90%%    %12.4f%12.4f%12.4f pc\n", rlag[0], rlag[1], rlag[2]);
        printf("    kinetic, potential      %12.4E%12.4E Msun km2/s2\n", ekin * v2, epot * v2);
        if (ltide) printf("    tidal energy            %12.4E Msun km2/s2\n", etid * v2);
        /* crossing time 2 r_h / sigma_3D, Spitzer half-mass relaxation time */
        printf("    crossing time 2r_h/sig  %12.4f Myr\n", 2.0 * rlag[1] / sig3d * tunit);
        if (ntot > 10)
            printf("    half-mass relaxation    %12.4f Myr\n",
                   0.138 * sqrt((double)ntot) * pow(rlag[1], 1.5) / sqrt(mt / ntot)
                       / log(0.11 * ntot) * tunit);
        writelevels();
    }
    printf("%s\n", DLINE);
    fflush(stdout);
}

/* screen output at each output time: one line (standard) or a block
 * (verbose); warnings in both */
void statusout(int iout, int nout, double time, int64_t nblock, int64_t nforce)
{
    char b1[32], b2[32], b3[32];
    diagnose(time);
    double de = (etot - etot0) / fabs(etot0), dei = (etot - eprev) / fabs(etot0);
    demax = fmax(demax, fabs(de));
    double w = walltime(), wi = w - wlast;
    int done = iout - ioutstart;
    double eta_s = w * (double)(nout - iout) / (double)(done > 1 ? done : 1);
    int64_t nbi = nblock - nblock0, nfi = nforce - nforce0;
    if (!lverbose) {
        if ((iout - 1) % 50 == 0)
            printf("   time(Myr)     done       dE/E       Q  r_h (pc)   wall time   time left\n");
        printf("%12.4f%8.1f%%%11.2E%8.3f%10.4f%12s%12s\n", time, 100.0 * iout / nout, de,
               -ekin / epot, rlag[1], hms(w, b1), hms(eta_s, b2));
    } else {
        printf("\n ---- t = %.4f Myr    (output %d of %d) ----------------------------------------\n",
               time, iout, nout);
        printf("   energy (Msun km2/s2)     kinetic %12.4E   potential %12.4E\n",
               ekin * vunit * vunit, epot * vunit * vunit);
        if (ltide)
            printf("                            tidal   %12.4E   total     %12.4E\n",
                   etid * vunit * vunit, etot * vunit * vunit);
        else printf("                            total   %12.4E\n", etot * vunit * vunit);
        printf("   energy error   since start%11.3E   this interval%11.3E   max%10.2E\n", de, dei,
               demax);
        printf("   virial ratio            %10.4f      sigma_3D%10.4f km/s\n", -ekin / epot,
               sig3d * vunit);
        printf("   Lagrangian radii (pc)   10/50/90%%  %10.4f%10.4f%10.4f\n", rlag[0], rlag[1],
               rlag[2]);
        printf("   centre of mass          %11.3E%11.3E%11.3E pc\n", cm[0], cm[1], cm[2]);
        printf("   c.o.m. velocity         %11.3E%11.3E%11.3E km/s\n", cv[0] * vunit,
               cv[1] * vunit, cv[2] * vunit);
        printf("   this interval   blocks  %12lld   star steps%14lld\n", (long long)nbi,
               (long long)nfi);
        if (nbi > 0) printf("                   mean stars per block%12.1f\n", (double)nfi / nbi);
        writelevels();
        printf("   steps held at kmax this interval %lld\n", (long long)(nclamp - nclamp0));
        if (closebinaries)
            printf("   groups as one body      %d   (binaries closed %lld, opened %lld this interval)\n",
                   nbin, (long long)(nclosed - nclosed0), (long long)(nopened - nopened0));
        if (fewbody)
            printf("   few-body groups         formed %lld, opened %lld this interval;  triples "
                   "frozen so far %lld;  steps %lld\n", (long long)(nfewformed - nfewformed0),
                   (long long)(nfewopened - nfewopened0), (long long)nfrozen, (long long)nfewsteps);
        if (usetree)
            printf("   tree nodes              %12d   builds so far %lld\n", nnode,
                   (long long)ntreebuilt);
        if (useac)
            printf("   neighbour scheme        %lld full sums, %lld neighbour sums this interval\n",
                   (long long)(nregular - nregular0), (long long)(nirregular - nirregular0));
        if (rcoll > 0.0)
            printf("   mergers                 %lld this interval, %lld so far;  stars left %d\n",
                   (long long)(nmerged - nmerged0), (long long)nmerged, ntot);
        printf("   wall time   interval%12s   total%12s   left%12s\n", hms(wi, b1), hms(w, b2),
               hms(eta_s, b3));
        if (wi > 0.0) printf("   speed      %11.3E star steps per second\n", (double)nfi / wi);
    }
    if (fabs(de) > dewarn)
        printf("  ** WARNING t = %.4f Myr: |dE/E| = %9.2E is above dewarn =%9.2E\n", time, fabs(de),
               dewarn);
    if (nclamp > nclamp0)
        printf("  ** NOTE t = %.4f Myr: %lld steps were held at the shortest timestep (kmax)\n",
               time, (long long)(nclamp - nclamp0));
    if (nmerged > nmerged0 && !lverbose)
        printf("  ** NOTE t = %.4f Myr: %lld merger(s) since the last output (%lld so far) - see "
               "%s_mergers.dat\n", time, (long long)(nmerged - nmerged0), (long long)nmerged, runid);
    eprev = etot;
    wlast = w;
    restartcounters(nblock, nforce);
    fflush(stdout);
}

void endinfo(double time, int64_t nforce)
{
    char b[32];
    printf("%s\n", DLINE);
    printf("  run %s finished at t = %.4f Myr\n", runid, time);
    printf("  star steps             %14lld\n", (long long)nforce);
    printf("  largest energy error   %14.3E\n", demax);
    printf("  wall time              %14s\n", hms(walltime(), b));
    if (closebinaries)
        printf("  binaries closed %lld, opened %lld, groups now %d\n", (long long)nclosed,
               (long long)nopened, nbin);
    if (fewbody)
        printf("  few-body groups formed %lld, opened %lld, triples frozen %lld\n",
               (long long)nfewformed, (long long)nfewopened, (long long)nfrozen);
    if (nmerged > 0) printf("  stellar mergers        %14lld\n", (long long)nmerged);
    printf("%s\n", DLINE);
}

/* ---- header lines:  # key (30 wide) = value -------------------- */
static void hdrs(FILE *f, const char *key, const char *val) { fprintf(f, "# %-30s= %s\n", key, val); }
static void hdrd(FILE *f, const char *key, double val) { fprintf(f, "# %-30s= %24.16E\n", key, val); }
static void hdri(FILE *f, const char *key, long long val) { fprintf(f, "# %-30s= %lld\n", key, val); }
static void hdrv(FILE *f, const char *key, const double *val, int nv)
{
    fprintf(f, "# %-30s=", key);
    for (int k = 0; k < nv; k++) fprintf(f, "%s%24.16E", k ? "" : " ", val[k]);
    fprintf(f, "\n");
}

/* the header lines common to every output file */
void hdrcommon(FILE *f)
{
    char date[64];
    time_t now = time(NULL);
    double mt = 0.0;
    strftime(date, sizeof date, "%Y-%m-%d %H:%M:%S", localtime(&now));
    for (int i = 0; i < n; i++) mt += m[i];
    hdrs(f, "creator", "HENDOS (C, src/main.c)");
    hdrs(f, "date", date);
    hdrs(f, "runid", runid);
    if (ensmember > 0) {
        hdrs(f, "ensemble", ensbase);
        hdri(f, "ensemble_member", ensmember);
        hdri(f, "ensemble_size", nensemble);
    }
    hdrs(f, "parameter_file", pfilename);
    hdrs(f, "initial_conditions", icfile);
    hdri(f, "n", ntot);
    hdrd(f, "total_mass_Msun", mt);
    hdrs(f, "units", "mass Msun, length pc, velocity km/s, time Myr, energy Msun (km/s)^2");
    hdrs(f, "frame", ltide ? "rotating with circular Galactic orbit: x away from Galactic "
                             "centre, y along rotation"
                           : "centre of mass of the initial conditions, non-rotating");
}

/* the header lines that depend on which physics is switched on: each
 * block is written only when its physics is on, so readers must treat
 * every key here as optional. >>> NEW PHYSICS ADDS ITS BLOCK HERE */
static void hdrphysics(FILE *f)
{
    char phys[256] = "gravity";
    double fconv = 1.0 / (vunit * 1.e3);
    if (ltide) strcat(phys, " tide");
    if (closebinaries) strcat(phys, " binaries");
    if (fewbody) strcat(phys, " few_body");
    if (rcoll > 0.0) strcat(phys, mergers ? " mergers" : " collision_stop");
    if (usetree) strcat(phys, " tree");
    if (useac) strcat(phys, " neighbour_scheme");
    hdrs(f, "physics", phys);
    hdrs(f, "tidal", tidal);
    if (ltide) {
        double tv[9];
        hdrd(f, "tidal_omega_km_s_kpc", omega / fconv);
        for (int k = 0; k < 9; k++) tv[k] = ttid[k / 3][k % 3] / (fconv * fconv);
        hdrv(f, "tidal_tensor_km2_s2_kpc2", tv, 9);
        if (!strcmp(tidal, "solar")) {
            hdrd(f, "oort_A_km_s_kpc", oorta);
            hdrd(f, "oort_B_km_s_kpc", oortb);
            hdrd(f, "rho0_Msun_pc3", rho0);
        }
    }
    if (closebinaries) {
        hdri(f, "groups_as_one_body", nbin);
        hdri(f, "binaries_closed", nclosed);
        hdri(f, "binaries_opened", nopened);
    }
    if (fewbody) {
        hdri(f, "few_body_groups_formed", nfewformed);
        hdri(f, "few_body_groups_opened", nfewopened);
        hdri(f, "triples_frozen", nfrozen);
        hdrd(f, "energy_change_frozen_triples", efrozennow * vunit * vunit);
    }
    if (rcoll > 0.0) {
        hdrd(f, "collision_radius_Rsun", rcoll);
        if (mergers) {
            hdri(f, "mergers", nmerged);
            hdri(f, "mergers_since_last_snapshot", nmerged - nmergedsnap);
            hdrd(f, "energy_removed_by_mergers", emerged * vunit * vunit);
        }
    }
    if (useac) {
        hdri(f, "neighbours_wanted", nnbopt);
        hdrd(f, "etar", etar);
    }
    if (usetree) {
        hdrd(f, "tree_opening_angle", theta);
        hdri(f, "tree_nleaf", nleaf);
        hdri(f, "tree_rebuild_interval", ntreebuild);
        hdrd(f, "etatree", etatree);
    }
}

/* snapshot file name  <datadir><runid>_snap_<time>Myr.dat */
void snapname(double time, char *fname)
{
    long long it = llround(time * 1.e4);
    sprintf(fname, "%s%s_snap_%06lld.%04lldMyr.dat", datadir, runid, it / 10000, it % 10000);
}

/* write a snapshot (NBODY-SNAPSHOT 2): header, then one line per star
 *   id m x y z vx vy vz type phi dt
 * all ntot stars, groups split into their stars, in id order. Must be
 * called just after diagnose (statusout or startinfo). */
void snapshot(double time, int64_t nforce)
{
    char fname[2 * NSTR], tc[256];
    snapname(time, fname);
    FILE *f = fopen(fname, "w");
    if (!f) {
        printf(" ERROR: cannot write %s\n", fname);
        exit(1);
    }
    double v2 = vunit * vunit, cvk[3] = {cv[0] * vunit, cv[1] * vunit, cv[2] * vunit};
    typecounts(tc);
    hdrs(f, "format", "NBODY-SNAPSHOT 2");
    hdrcommon(f);
    hdrs(f, "particle_types", tc);
    hdrd(f, "time_Myr", time);
    hdri(f, "bodies", n);
    hdrphysics(f);
    hdrd(f, "softening_pc", eps);
    hdrd(f, "eta", eta);
    hdrd(f, "etap", etap);
    hdri(f, "niter", niter);
    hdri(f, "kmax", kmax);
    hdrd(f, "dtmax_yr", tout * 1.e6);
    hdrd(f, "energy_kinetic", ekin * v2);
    hdrd(f, "energy_potential", epot * v2);
    hdrd(f, "energy_tidal", etid * v2);
    hdrd(f, "energy_total", etot * v2);
    hdrd(f, "energy_initial", etot0 * v2);
    hdrd(f, "energy_error", (etot - etot0) / fabs(etot0));
    hdrd(f, "energy_error_max", demax);
    hdrd(f, "virial_ratio", -ekin / epot);
    hdrv(f, "lagrangian_radii_10_50_90_pc", rlag, 3);
    hdrv(f, "centre_of_mass_pc", cm, 3);
    hdrv(f, "centre_of_mass_velocity_km_s", cvk, 3);
    hdrd(f, "wall_time_s", walltime());
    hdri(f, "star_steps", nforce);
    hdrs(f, "columns", "id m x y z vx vy vz type phi dt");
    hdrs(f, "column_units", "- Msun pc pc pc km/s km/s km/s - (km/s)^2 yr");
    hdrs(f, "column_notes", "type = particle type; phi = gravitational + tidal potential per "
                            "unit mass; dt = current timestep");
    fprintf(f, "# end_header\n");

    int ns = 0, *sid = xalloc(sizeof(int) * (size_t)ntot), *styp = xalloc(sizeof(int) * (size_t)ntot);
    int *sord = xalloc(sizeof(int) * (size_t)ntot);
    double (*sdat)[9] = xalloc(sizeof(*sdat) * (size_t)ntot);
    double *skey = xalloc(sizeof(double) * (size_t)ntot);
    for (int i = 0; i < n; i++) {
        double phit = phi[i], dtyr = (double)dtk[i] * tick * tunit * 1.e6;
        if (ltide)
            for (int c = 0; c < 3; c++) phit -= 0.5 * r[i][c] * dot3(ttid[c], r[i]);
        int k = ib[i];
        if (k < 0) {
            double row[9] = {m[i], r[i][0], r[i][1], r[i][2], v[i][0] * vunit, v[i][1] * vunit,
                             v[i][2] * vunit, phit * v2, dtyr};
            sid[ns] = id[i];
            styp[ns] = itype[i];
            memcpy(sdat[ns++], row, sizeof row);
        } else {
            int nm;
            vec3 xm[KFMAX], um[KFMAX];
            grpcomp(k, r[i], v[i], time / tunit, &nm, xm, um);
            for (int l = 0; l < nm; l++) {
                double phim = phit;
                for (int l2 = 0; l2 < nm; l2++) {
                    if (l2 == l) continue;
                    double d[3] = {xm[l][0] - xm[l2][0], xm[l][1] - xm[l2][1], xm[l][2] - xm[l2][2]};
                    phim -= bmass[k][l2] / norm3(d);
                }
                double row[9] = {bmass[k][l], xm[l][0], xm[l][1], xm[l][2], um[l][0] * vunit,
                                 um[l][1] * vunit, um[l][2] * vunit, phim * v2, dtyr};
                sid[ns] = bident[k][l];
                styp[ns] = btype[k][l];
                memcpy(sdat[ns++], row, sizeof row);
            }
        }
    }
    for (int i = 0; i < ns; i++) skey[i] = (double)sid[i];
    indexsort(ns, skey, sord);
    for (int q = 0; q < ns; q++) {
        int i = sord[q];
        fprintf(f, "%10d", sid[i]);
        for (int c = 0; c < 7; c++) fprintf(f, "%24.16E", sdat[i][c]);
        fprintf(f, " %-12s%16.8E%16.8E\n", typename_[styp[i]], sdat[i][7], sdat[i][8]);
    }
    fclose(f);
    free(sid); free(styp); free(sord); free(sdat); free(skey);
    nmergedsnap = nmerged;
    if (lverbose) printf("   snapshot written: %s\n", fname);
}

/* ---- IN-RUN ANALYSIS ----------------------------------------------
 * Each analysis has an on/off switch in the &analysis section of the
 * parameter file and writes its own file  <runid>_ana_<name>.dat (in
 * the data directory) with a metadata header (same style as the
 * snapshots) followed by one line per analysis time. Analyses are done
 * every tana Myr, a multiple of tout, so all stars are synchronised.
 * On a restart existing files are appended to rather than overwritten.
 * ana_do is called just after statusout, so the diagnostics (cm, rlag,
 * sig3d) are up to date and can be used without recomputing them.
 * Energies (direct-N, expensive) are not computed by analyses.
 *
 * >>> TO ADD A NEW ANALYSIS <name>:
 *   1. add an int switch <name> to hendos.h, a line for it in the
 *      params[] table in params.c with its default in readparams, and
 *      a line in run/params.nml
 *   2. add a FILE * for it next to flagr at the top of this file
 *   3. open it and write its header in ana_init (copy the lagrangian
 *      block)
 *   4. write its line in ana_do
 *   5. close it in ana_end
 *   6. add its file to the check in checkrunid (bottom of this file)
 * ------------------------------------------------------------------ */
static void ananame(const char *name, char *fname)
{
    sprintf(fname, "%s%s_ana_%s.dat", datadir, runid, name);
}

void ana_init(void)
{
    char fname[2 * NSTR];
    struct stat st;
    if (!lagrangian) return;
    ananame("lagrangian", fname);
    if ((time0 > 0.0 || restart_) && stat(fname, &st) == 0) {
        flagr = fopen(fname, "a");            /* carry on with the existing file */
    } else {
        flagr = fopen(fname, "w");
        if (!flagr) die("cannot write the analysis file");
        hdrs(flagr, "format", "NBODY-ANALYSIS 1");
        hdrs(flagr, "analysis", "lagrangian");
        hdrcommon(flagr);
        hdrd(flagr, "analysis_interval_Myr", tana);
        hdrs(flagr, "columns", "time r1% r2% r5% r10% r20% r30% r50% r70% r90% sigma3d");
        hdrs(flagr, "column_units", "Myr pc pc pc pc pc pc pc pc pc km/s");
        fprintf(flagr, "# end_header\n");
    }
    if (lverbose) printf("   analysis file: %s\n", fname);
}

void ana_do(double time)
{
    static const double fl[9] = {0.01, 0.02, 0.05, 0.1, 0.2, 0.3, 0.5, 0.7, 0.9};
    double rl[9];
    if (!lagrangian || !flagr) return;
    lagradii(9, fl, rl);
    fprintf(flagr, "%16.6f", time);
    for (int k = 0; k < 9; k++) fprintf(flagr, "%15.6E", rl[k]);
    fprintf(flagr, "%15.6E\n", sig3d * vunit);
    fflush(flagr);
}

void ana_end(void)
{
    if (flagr) fclose(flagr);
}

/* for a new run stop if output files with this runid already exist, so
 * an earlier run is never overwritten by accident */
void checkrunid(int nsnap, int nana)
{
    char fname[2 * NSTR];
    struct stat st;
    int clash = 0;
    if (nsnap > 0) {
        snapname(0.0, fname);
        if (stat(fname, &st) == 0 && strcmp(fname, icfile)) clash = 1;
    }
    if (nana > 0 && lagrangian) {
        ananame("lagrangian", fname);
        if (stat(fname, &st) == 0) clash = 1;
    }
    sprintf(fname, "%s%s_mergers.dat", datadir, runid);
    if (stat(fname, &st) == 0) clash = 1;
    sprintf(fname, "%s%s_restart.bin", datadir, runid);
    if (stat(fname, &st) == 0) clash = 1;
    if (clash) {
        printf(" ERROR: output files for runid '%s' already exist here.\n", runid);
        printf("        Choose a new runid in the parameter file, or move or delete the old files.\n");
        exit(1);
    }
}
