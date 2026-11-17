/* ==============================================================
 * HENDOS - HErmite N-body Dynamics Of Systems
 * The main program - main.c
 *
 * 4th-order Hermite (Makino & Aarseth
 * 1992) with block timesteps (Makino 1991), timestep = min(Aarseth,
 * pairwise free-fall/flyby), higher-order predictor, P(EC)^n, Galactic
 * tide, Ahmad-Cohen neighbour scheme, binaries / few-body groups /
 * frozen triples as one body, mergers, restart files and ensembles.
 * Forces by direct summation, the neighbour scheme (ac.c) or the
 * Barnes-Hut tree (tree.c).
 *
 * usage:   ./hermite params.nml  [ member ]
 *
 * This file is the main program only: set-up, the ensemble launcher,
 * then the loop over output intervals and blocks (predict, forces,
 * correct, new timesteps). Block times are 8-byte integers in units of
 * tick = dtmax/2^kmax, so they are exact; everyone is synchronised
 * every dtmax (= tout), when the integer clock is reset to zero.
 * ============================================================== */
#define HENDOS_MAIN
#include <omp.h>
#include <sys/stat.h>
#include "hendos.h"
double dbgenergy(double tnow);

static char ensicfile[NSTR];

/* this run is member k of an ensemble: its runid and its initial
 * conditions file get the member number, so every file it writes has
 * it too */
static void setmember(int k)
{
    char tmp[NSTR];
    if (k < 1 || k > nensemble) {
        printf(" ERROR: ensemble member %d asked for, but nensemble in the parameter file "
               "is %d\n", k, nensemble);
        exit(1);
    }
    ensmember = k;
    strcpy(ensbase, runid);
    strcpy(ensicfile, icfile);
    ensname(ensbase, k, tmp);
    snprintf(runid, sizeof runid, "%s", tmp);
    ensname(ensicfile, k, tmp);
    snprintf(icfile, sizeof icfile, "%s", tmp);
}

/* run every member of the ensemble (nensemble > 1 in the parameter file
 * and no member number on the command line).
 *
 * Each member is a separate run of this same program, started here as
 *     hermite <parameter file> <member number>
 * so the members cannot interfere with each other and nothing in the
 * code has to be reset between them. Member k reads the initial
 * conditions ensname(icfile,k) (as written by makeics with the same
 * nensemble) and has the runid ensname(runid,k), so its files are
 *     <runid>_e00k_snap_..., <runid>_e00k_ana_..., <runid>_e00k_mergers.dat
 * and its screen output goes to <datadir>/<runid>_e00k.log.
 * nparallel members are run at the same time; each uses OMP_NUM_THREADS
 * threads, so nparallel x OMP_NUM_THREADS should not be more than the
 * number of cores. For small N, many members with one thread each is
 * by far the best use of a machine.
 * A member can also be run on its own in the same way by hand (or as
 * one job of a job array on a cluster):   ./hermite params.nml 7
 *
 * ADDING MEMBERS LATER. An ensemble can be made bigger at any time:
 * make the extra initial conditions (makeics with the larger nensemble
 * and ensfirst = the first new member), raise nensemble here and run
 * again. Members that have been run already - those whose log file
 * exists - are left alone and only the new ones are run, so nothing
 * else has to be changed. (ensfirst in &run does the same by hand:
 * only members ensfirst to nensemble are looked at.) The members of an
 * ensemble need not all have been run with the same nensemble; the
 * analysis tools take the ensemble to be all the members they find.
 * To run a member again, delete its files first.
 */
static void runensemble(const char *exe)
{
    char name[NSTR], path[3 * NSTR], ic1[NSTR], ic2[NSTR];
    int *todo = xalloc(sizeof(int) * (size_t)nensemble), ntodo = 0, ndone = 0, nmiss = 0,
        nfail = 0;
    struct stat st;
    const char *redir = restart_ ? ">>" : ">";
    ensname(icfile, 1, ic1);
    ensname(icfile, nensemble, ic2);
    printf(" ==============================================================================\n");
    printf("  HENDOS - HErmite N-body Dynamics Of Systems     ensemble: %s\n", runid);
    printf(" ==============================================================================\n");
    printf("  members                %14d\n", nensemble);
    if (ensfirst > 1) printf("  starting at member     %14d\n", ensfirst);
    printf("  run at the same time   %14d\n", nparallel);
    printf("  threads for each       %14d\n", nthreads);
    printf("  initial conditions     %s to %s\n", ic1, ic2);
    printf("  screen output of each  %s<runid>_e<member>.log\n", datadir);
    if (restart_) printf("  every member is restarted from its restart file\n");
    for (int k = ensfirst; k <= nensemble; k++) {
        ensname(runid, k, name);
        snprintf(path, sizeof path, "%s%s.log", datadir, name);
        if (!restart_ && stat(path, &st) == 0) ndone++;
        else todo[ntodo++] = k;
    }
    if (ndone > 0) {
        printf("  already run            %14d  (their log files exist: left alone)\n", ndone);
        printf("  to run now             %14d\n", ntodo);
    }
    if (ntodo == 0) {
        printf("  Every member has been run already. To add members raise nensemble; to run one\n");
        printf("  again delete its files (including its .log file) first.\n");
        return;
    }
    if (!restart_)
        for (int q = 0; q < ntodo; q++) {
            ensname(icfile, todo[q], name);
            if (stat(name, &st) != 0 && ++nmiss <= 5) printf(" ERROR: cannot find %s\n", name);
        }
    if (nmiss > 0) {
        printf(" ERROR: %d of the initial conditions files are missing. Make them with makeics "
               "and the same nensemble\n", nmiss);
        printf("        (and ensfirst = the first new member, if members are being added).\n");
        exit(1);
    }
    printf(" ------------------------------------------------------------------------------\n");
    fflush(stdout);
    char *cmd = xalloc((size_t)nparallel * 4 * NSTR + 16);
    for (int k1 = 0; k1 < ntodo; k1 += nparallel) {
        int k2 = k1 + nparallel < ntodo ? k1 + nparallel : ntodo;
        cmd[0] = 0;
        for (int q = k1; q < k2; q++) {
            ensname(runid, todo[q], name);
            sprintf(cmd + strlen(cmd), "( %s %s %d %s %s%s.log 2>&1 ; echo $? > %s%s.status ) & ",
                    exe, pfilename, todo[q], redir, datadir, name, datadir, name);
        }
        strcat(cmd, "wait");
        double w0 = omp_get_wtime();
        if (system(cmd) != 0) { /* the status files say what happened */ }
        for (int q = k1; q < k2; q++) {
            int ist = -1;
            ensname(runid, todo[q], name);
            snprintf(path, sizeof path, "%s%s.status", datadir, name);
            FILE *f = fopen(path, "r");
            if (f) {
                if (fscanf(f, "%d", &ist) != 1) ist = -1;
                fclose(f);
                remove(path);
            }
            if (ist == 0) printf("  member %6d   finished   %s\n", todo[q], name);
            else {
                nfail++;
                printf("  member %6d   STOPPED EARLY - see %s%s.log\n", todo[q], datadir, name);
            }
        }
        printf("    (members %d to %d took %.1f s)\n", todo[k1], todo[k2 - 1], omp_get_wtime() - w0);
        fflush(stdout);
    }
    printf(" ==============================================================================\n");
    printf("  ensemble finished: %d of %d members ran to the end\n", ntodo - nfail, ntodo);
    if (ndone > 0) printf("  (%d members had been run before)\n", ndone);
    if (nfail > 0) printf("  ** %d stopped early: look at their log files\n", nfail);
    printf(" ==============================================================================\n");
    if (nfail > 0) exit(1);
}

int main(int argc, char **argv)
{
    int nsince = 0, iout0 = 0, nout, nsnap = 0, nana = 0, nrest = 0, member = 0, iout;
    int64_t nblock = 0, nforce = 0, tblock;
    double time = 0.0, tbase;
    char fname[2 * NSTR];
    double tracet1 = 0.0, tracet2 = 0.0, tracee = 0.0;
    if (getenv("HENDOS_TRACE")) sscanf(getenv("HENDOS_TRACE"), "%lf %lf", &tracet1, &tracet2);

    snprintf(pfilename, sizeof pfilename, "%s", argc >= 2 ? argv[1] : "params.nml");
    setunits();
    readparams(pfilename);
    checkdatadir();

    /* ensembles: with a member number this run is that member; without
     * one, start a run for every member in turn and stop */
    if (argc >= 3) {
        char *end;
        member = (int)strtol(argv[2], &end, 10);
        if (*end || member < 1) {
            printf(" ERROR: the ensemble member number is not a whole number > 0: %s\n", argv[2]);
            return 1;
        }
    }
    if (member > 0) setmember(member);
    else if (nensemble > 1) {
        runensemble(argv[0]);
        return 0;
    }
    eps2 = eps * eps;

    /* initial conditions, or the complete state of an earlier run */
    if (restart_) {
        readrestart(&time, &iout0, &nblock, &nforce);
    } else {
        readic();
        accheck(ntot);          /* (switched off if there are too few stars) */
        if (useac) acalloc();
    }
    if (usetree) treealloc();

    nout = (int)lround((tend - time0) / tout);
    if (nout < 1) {
        printf(" ERROR: tend must be at least one tout after the start time of the initial "
               "conditions (%.4f Myr)\n", time0);
        return 1;
    }
    if (restart_ && iout0 >= nout) {
        printf(" ERROR: the restart file is already at %.4f Myr - increase tend to carry the "
               "run on\n", time);
        return 1;
    }
    if (tsnap > 0.0) nsnap = (int)fmax(1.0, round(tsnap / tout));
    if (tana > 0.0) nana = (int)fmax(1.0, round(tana / tout));
    if (trestart > 0.0) nrest = (int)fmax(1.0, round(trestart / tout));

    /* a new run must not overwrite the files of an earlier run */
    if (time0 == 0.0 && !restart_) checkrunid(nsnap, nana);

    /* block timestep set-up */
    dtmax = tout / tunit;
    tick = dtmax / ldexp(1.0, kmax);
    tend8 = (int64_t)1 << kmax;

    if (restart_) {
        if (nana > 0) ana_init();
    } else {
        nclamp = 0;
        ioutstart = 0;
        firststep();
        if (rcoll > 0.0) collcheck(0, time0 / tunit);
        if (closebinaries) binstart(time0 / tunit);
        startinfo(nout);
        time = time0;
        if (nsnap > 0) {
            snapname(time, fname);       /* (not when starting from that same snapshot) */
            if (strcmp(fname, icfile)) snapshot(time, nforce);
        }
        if (nana > 0) {
            ana_init();
            if (time0 == 0.0) ana_do(time);
        }
    }

    /* ================================================================
     * main loop - over output intervals, then over blocks */
    for (iout = iout0 + 1; iout <= nout; iout++) {
        tbase = (time0 + (double)(iout - 1) * tout) / tunit;
        /* the tree is rebuilt at the first block of every output interval,
         * so that a run restarted at an output time is exactly the same as
         * an unbroken one */
        nsince = ntreebuild;
        for (;;) {
            /* next block time and the active particles (in order of index) */
            tblock = INT64_MAX;
#pragma omp parallel for reduction(min:tblock) schedule(static) if (n >= nompj)
            for (int i = 0; i < n; i++)
                if (tk[i] + dtk[i] < tblock) tblock = tk[i] + dtk[i];
            nact = 0;
            for (int i = 0; i < n; i++)
                if (tk[i] + dtk[i] == tblock) ilist[nact++] = i;

            /* predictor - ALL particles to tblock for the force calculation,
             * with f2e, f3e (2nd and 3rd derivatives from each star's last
             * step): essential with the tide, whose Coriolis force depends on
             * velocity */
#pragma omp parallel for schedule(static) if (n >= nompj)
            for (int i = 0; i < n; i++) {
                double dtp = (double)(tblock - tk[i]) * tick;
                double d2 = dtp * dtp / 2.0, d3 = d2 * dtp / 3.0, d4 = d3 * dtp / 4.0,
                       d5 = d4 * dtp / 5.0;
                for (int k = 0; k < 3; k++) {
                    rp[i][k] = r[i][k] + v[i][k] * dtp + f0[i][k] * d2 + f1[i][k] * d3
                               + f2e[i][k] * d4 + f3e[i][k] * d5;
                    vp[i][k] = v[i][k] + f0[i][k] * dtp + f1[i][k] * d2 + f2e[i][k] * d3
                               + f3e[i][k] * d4;
                }
            }

            /* tree: rebuild it every ntreebuild block steps (or after a group
             * was closed or opened, or stars merged); on the steps in between
             * only restock it - recompute the moments from the new predicted
             * positions */
            if (usetree) {
                if (!treevalid || nsince >= ntreebuild) {
                    treebuild();
                    nsince = 0;
                } else {
                    treemoments();
                }
                nsince++;
            }

            /* base for the corrector of the active stars: the 3rd-order
             * prediction from F0 and F1 only */
#pragma omp parallel for schedule(static) if (nact >= nompmin)
            for (int ia = 0; ia < nact; ia++) {
                int i = ilist[ia];
                double dt = (double)dtk[i] * tick;
                for (int k = 0; k < 3; k++) {
                    rc[i][k] = r[i][k] + v[i][k] * dt + f0[i][k] * dt * dt / 2.0
                               + f1[i][k] * dt * dt * dt / 6.0;
                    vc[i][k] = v[i][k] + f0[i][k] * dt + f1[i][k] * dt * dt / 2.0;
                }
            }

            /* evaluate and correct the active particles (niter times) */
            tcur8 = tblock;
            for (int it = 1; it <= niter; it++) {
                aclast = (it == niter);
                forces();
#pragma omp parallel for schedule(static) if (nact >= nompmin)
                for (int ia = 0; ia < nact; ia++) {
                    int i = ilist[ia];
                    double dt = (double)dtk[i] * tick, dt2 = dt * dt, dt3 = dt2 * dt,
                           dt4 = dt3 * dt, dt5 = dt4 * dt;
                    for (int k = 0; k < 3; k++) {
                        double a1 = f0[i][k] - f0p[i][k], a2 = f1[i][k] + f1p[i][k];
                        f2[i][k] = 2.0 * (-3.0 * a1 - (f1[i][k] + a2) * dt) / dt2;
                        f3[i][k] = 6.0 * (2.0 * a1 + a2 * dt) / dt3;
                        rp[i][k] = rc[i][k] + f2[i][k] * dt4 / 24.0 + f3[i][k] * dt5 / 120.0;
                        vp[i][k] = vc[i][k] + f2[i][k] * dt3 / 6.0 + f3[i][k] * dt4 / 24.0;
                        /* neighbour scheme: on a regular step, the correction from
                         * the regular force over the whole regular step as well */
                        if (useac && lreg[i]) {
                            rp[i][k] += dreg[i][k];
                            vp[i][k] += dvreg[i][k];
                        }
                    }
                }
            }

            /* accept the step and get new levels */
#pragma omp parallel for schedule(static) if (nact >= nompmin)
            for (int ia = 0; ia < nact; ia++) {
                int i = ilist[ia];
                double dt = (double)dtk[i] * tick;
                for (int k = 0; k < 3; k++) {
                    r[i][k] = rp[i][k];
                    v[i][k] = vp[i][k];
                    f0[i][k] = f0p[i][k];
                    f1[i][k] = f1p[i][k];
                    f2e[i][k] = f2[i][k] + f3[i][k] * dt;
                    f3e[i][k] = f3[i][k];
                }
                tk[i] = tblock;
                if (useac) acaccept(i);
                newlevel(i);
                if (useac) aclevel(i);
            }
            nblock++;
            nforce += (int64_t)nact * niter;

            /* stellar collisions, then open, close and freeze groups */
            if (rcoll > 0.0) collcheck(tblock, tbase + (double)tblock * tick);
            if (closebinaries) bincheck(tblock, tbase + (double)tblock * tick);

            /* debugging: HENDOS_TRACE="t1 t2" (Myr) prints every block between
             * t1 and t2 in which the total energy changed by more than 1e-7 */
            if (tracet1 < tracet2) {
                double tn = (tbase + (double)tblock * tick) * tunit;
                if (tn >= tracet1 && tn <= tracet2) {
                    double e = dbgenergy(tbase + (double)tblock * tick);
                    if (tracee != 0.0 && fabs(e - tracee) > 1e-7 * fabs(tracee)) {
                        printf("  TRACE t = %.8f Myr  dE/|E| = %10.3E  active %d:", tn,
                               (e - tracee) / fabs(tracee), nact);
                        for (int q = 0; q < nact && q < 6; q++)
                            printf(" id %d (group %d, level %d)", id[ilist[q]], ib[ilist[q]],
                                   lev[ilist[q]]);
                        printf("\n");
                    }
                    tracee = e;
                }
            }
            if (tblock == tend8) break;          /* everyone is synchronised here */
        }

        /* reset the integer clock (regular times too) */
        for (int i = 0; i < n; i++) {
            tk[i] = 0;
            if (useac) treg[i] -= tend8;
        }
        time = time0 + (double)iout * tout;

        /* output: screen first - it also updates the diagnostics */
        statusout(iout, nout, time, nblock, nforce);
        if (nsnap > 0 && iout % nsnap == 0) snapshot(time, nforce);
        if (nana > 0 && iout % nana == 0) ana_do(time);
        if (nrest > 0 && (iout % nrest == 0 || iout == nout)) writerestart(time, iout, nblock, nforce);
    }
    if (nana > 0) ana_end();
    endinfo(time, nforce);
    return 0;
}
