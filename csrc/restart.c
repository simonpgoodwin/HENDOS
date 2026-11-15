/* ==============================================================
 * HENDOS (C version) - restart.c
 * RESTART FILES (parameters trestart and restart) - for long runs.
 * The counterpart of nbrestart.f90: the complete state of the run at
 * an output time, in binary at full precision, so that a run carried
 * on from it does what the original would have done.
 *
 * NOTE: the file has the same name as the Fortran one
 * (<runid>_restart.bin) but it is NOT the same file: a run must be
 * restarted by the version of the code that wrote the file. The C file
 * starts with the words HENDOS-C-RESTART, and the C code refuses
 * anything else.
 *
 * FILE FORMAT - made to be added to, as in the Fortran. After the
 * 16 opening characters the file is a list of sections, each
 *     name (16 characters), number of bytes (8-byte integer), the bytes
 * A reader that does not know a section skips its bytes, so a section
 * for new physics can be added without breaking anything.
 * >>> TO ADD A SECTION: write it in writerestart between secbegin and
 *     secend, and add a branch for its name in readrestart
 * ============================================================== */
#include "hendos.h"

#define RSTVERSION 1
static const char MAGIC[16] = "HENDOS-C-RESTART";

static void restartname(char *fname) { sprintf(fname, "%s%s_restart.bin", datadir, runid); }

static long secpos;
static void secbegin(FILE *f, const char *name)
{
    char tag[16] = {0};
    int64_t nb = 0;
    strncpy(tag, name, 15);
    fwrite(tag, 1, 16, f);
    secpos = ftell(f);
    fwrite(&nb, sizeof nb, 1, f);
}
static void secend(FILE *f)
{
    long here = ftell(f);
    int64_t nb = here - secpos - (long)sizeof(int64_t);
    fseek(f, secpos, SEEK_SET);
    fwrite(&nb, sizeof nb, 1, f);
    fseek(f, here, SEEK_SET);
}
#define WR(f, a, cnt) fwrite((a), sizeof(*(a)), (size_t)(cnt), (f))
#define WR1(f, x) fwrite(&(x), sizeof(x), 1, (f))
static FILE *rf;
static void rd(void *p, size_t size, size_t cnt)
{
    if (fread(p, size, cnt, rf) != cnt) die("the restart file is incomplete");
}
#define RD(a, cnt) rd((a), sizeof(*(a)), (size_t)(cnt))
#define RD1(x) rd(&(x), sizeof(x), 1)

/* write the restart file. Must be called at an output time, just after
 * statusout (everyone is synchronised and the integer clock is zero).
 * Written under a temporary name first, so a crash while writing never
 * destroys the previous file. */
void writerestart(double time, int iout, int64_t nblock, int64_t nforce)
{
    char fname[2 * NSTR], tmp[2 * NSTR + 8];
    int iver = RSTVERSION, kf = KFMAX, nnbm = useac ? nnbmax : 0;
    restartname(fname);
    sprintf(tmp, "%s.tmp", fname);
    FILE *f = fopen(tmp, "wb");
    if (!f) {
        printf(" ERROR: cannot write the restart file %s\n", tmp);
        exit(1);
    }
    fwrite(MAGIC, 1, 16, f);

    secbegin(f, "HEADER");
    WR1(f, iver); WR1(f, nalloc); WR1(f, n); WR1(f, ntot); WR1(f, nbin); WR1(f, kf); WR1(f, nnbm);
    secend(f);

    secbegin(f, "RUN");
    WR1(f, time); WR1(f, time0); WR1(f, iout); WR1(f, nblock); WR1(f, nforce);
    WR(f, comshift, 3); WR(f, covshift, 3); fwrite(icformat, 1, sizeof icformat, f);
    secend(f);

    secbegin(f, "COUNTERS");
    WR1(f, nclamp); WR1(f, nclosed); WR1(f, nopened); WR1(f, nfewformed); WR1(f, nfewopened);
    WR1(f, nfrozen); WR1(f, nfewsteps); WR1(f, nmerged); WR1(f, emerged); WR1(f, nregular);
    WR1(f, nirregular);
    secend(f);

    secbegin(f, "ENERGY");
    WR1(f, etot0); WR1(f, eprev); WR1(f, demax); WR1(f, nmergedsnap);
    secend(f);

    secbegin(f, "PARTICLES");
    WR(f, r, n); WR(f, v, n); WR(f, f0, n); WR(f, f1, n); WR(f, f2, n); WR(f, f3, n);
    WR(f, f2e, n); WR(f, f3e, n); WR(f, m, n); WR(f, tpair2, n); WR(f, rmin2, n); WR(f, phi, n);
    WR(f, id, n); WR(f, itype, n); WR(f, ib, n); WR(f, lev, n); WR(f, dtk, n);
    secend(f);

    if (nbin > 0) {
        secbegin(f, "GROUPS");
        WR(f, bnum, nbin); WR(f, bmode, nbin); WR(f, bmass, nbin); WR(f, bident, nbin);
        WR(f, btype, nbin); WR(f, brel, nbin); WR(f, brel2, nbin); WR(f, bepoch, nbin);
        WR(f, bx, nbin); WR(f, bv, nbin);
        secend(f);
    }

    if (useac) {
        secbegin(f, "NEIGHBOURS");
        WR1(f, rs0); WR1(f, nnbopt);
        WR(f, fi, n); WR(f, fid, n); WR(f, fr, n); WR(f, frd, n); WR(f, rs, n); WR(f, nnb, n);
        WR(f, nbl, (size_t)n * (size_t)nnbmax); WR(f, levr, n); WR(f, dtreg, n); WR(f, acok, n);
        secend(f);
    }

    if (nbin > 0 || efrozen != 0.0) {
        secbegin(f, "FROZEN");
        WR1(f, efrozen); WR(f, bein, nbin);
        secend(f);
    }

    /* >>> NEW SECTIONS GO HERE */

    secbegin(f, "END");
    secend(f);
    fclose(f);
    if (rename(tmp, fname) != 0) printf(" WARNING: could not rename %s\n", tmp);
    if (lverbose) printf("   restart file written: %s\n", fname);
}

/* read the restart file, allocate the arrays and set everything up to
 * carry on from the time it was written */
void readrestart(double *time, int *iout, int64_t *nblock, int64_t *nforce)
{
    char fname[2 * NSTR], tag[17] = {0}, magic[16];
    int gothead = 0, gotpart = 0, gotac = 0, iver, nn, nt, nb, kf, nnbm = 0;
    int64_t nbytes;
    restartname(fname);
    rf = fopen(fname, "rb");
    if (!rf) {
        printf(" ERROR: restart = .true. but there is no restart file %s\n", fname);
        exit(1);
    }
    if (fread(magic, 1, 16, rf) != 16 || memcmp(magic, MAGIC, 16)) {
        printf(" ERROR: %s was not written by the C version of HENDOS.\n", fname);
        printf("        A run must be restarted by the version (Fortran or C) that wrote "
               "the file.\n");
        exit(1);
    }
    for (;;) {
        rd(tag, 1, 16);
        RD1(nbytes);
        if (!strcmp(tag, "END")) break;
        if (strcmp(tag, "HEADER") && !gothead) die("the restart file does not start with its header");
        if (!strcmp(tag, "HEADER")) {
            RD1(iver); RD1(nalloc); RD1(nn); RD1(nt); RD1(nb); RD1(kf); RD1(nnbm);
            if (iver > RSTVERSION || kf != KFMAX)
                die("the restart file was written by a different version of HENDOS");
            /* (the neighbour scheme is switched off if the run began with
             * too few stars, as it was when the file was written) */
            accheck(nalloc);
            if ((nnbm > 0) != (useac != 0))
                die("useac must be the same as when the restart file was written");
            /* the arrays are allocated for the number of stars the run began with */
            n = nalloc;
            allocarrays();
            if (useac) {
                nnbmax = nnbm;
                if (nnbopt < 1) nnbopt = 1;
                acalloc();
            }
            n = nn;
            ntot = nt;
            nbin = nb;
            gothead = 1;
        } else if (!strcmp(tag, "RUN")) {
            RD1(*time); RD1(time0); RD1(*iout); RD1(*nblock); RD1(*nforce);
            RD(comshift, 3); RD(covshift, 3); rd(icformat, 1, sizeof icformat);
        } else if (!strcmp(tag, "COUNTERS")) {
            RD1(nclamp); RD1(nclosed); RD1(nopened); RD1(nfewformed); RD1(nfewopened);
            RD1(nfrozen); RD1(nfewsteps); RD1(nmerged); RD1(emerged); RD1(nregular);
            RD1(nirregular);
        } else if (!strcmp(tag, "ENERGY")) {
            RD1(etot0); RD1(eprev); RD1(demax); RD1(nmergedsnap);
        } else if (!strcmp(tag, "PARTICLES")) {
            RD(r, n); RD(v, n); RD(f0, n); RD(f1, n); RD(f2, n); RD(f3, n);
            RD(f2e, n); RD(f3e, n); RD(m, n); RD(tpair2, n); RD(rmin2, n); RD(phi, n);
            RD(id, n); RD(itype, n); RD(ib, n); RD(lev, n); RD(dtk, n);
            gotpart = 1;
        } else if (!strcmp(tag, "GROUPS")) {
            RD(bnum, nbin); RD(bmode, nbin); RD(bmass, nbin); RD(bident, nbin);
            RD(btype, nbin); RD(brel, nbin); RD(brel2, nbin); RD(bepoch, nbin);
            RD(bx, nbin); RD(bv, nbin);
        } else if (!strcmp(tag, "FROZEN")) {
            RD1(efrozen); RD(bein, nbin);
        } else if (!strcmp(tag, "NEIGHBOURS")) {
            if (!useac) die("useac must be the same as when the restart file was written");
            RD1(rs0); RD1(nnbopt);
            RD(fi, n); RD(fid, n); RD(fr, n); RD(frd, n); RD(rs, n); RD(nnb, n);
            RD(nbl, (size_t)n * (size_t)nnbmax); RD(levr, n); RD(dtreg, n); RD(acok, n);
            gotac = 1;
        } else {
            printf(" NOTE: section %s of the restart file is not known to this version of "
                   "HENDOS and was skipped\n", tag);
            if (fseek(rf, (long)nbytes, SEEK_CUR)) die("the restart file is incomplete");
        }
    }
    fclose(rf);
    if (!gotpart || (useac && !gotac)) die("the restart file is incomplete");

    /* everyone is synchronised at the time of the file, integer clock zero */
    for (int i = 0; i < n; i++) {
        memcpy(rp[i], r[i], sizeof(vec3));
        memcpy(vp[i], v[i], sizeof(vec3));
        tk[i] = 0;
        if (useac) {
            treg[i] = 0;
            lreg[i] = 0;
        }
    }
    printf(" ==============================================================================\n");
    printf("  HENDOS - HErmite N-body Dynamics Of Systems (C version)     run: %s\n", runid);
    printf(" ==============================================================================\n");
    printf("  restarted from         %s\n", fname);
    printf("  time                   %14.4f Myr\n", *time);
    printf("  end time               %14.4f Myr\n", tend);
    printf("  particles              %14d\n", ntot);
    printf("  OpenMP threads         %14d\n", nthreads);
    printf("  energy error so far    %14.3E\n", demax);
    printf(" ==============================================================================\n");
    startclock();
    ioutstart = *iout;
    restartcounters(*nblock, *nforce);
}
