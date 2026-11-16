/* ==============================================================
 * HENDOS - HErmite N-body Dynamics Of Systems - params.c
 * SET-UP: units and the parameter file.
 *
 * The parameter file is written as a Fortran namelist (run/params.nml),
 * the format of the original Fortran version (now in archive/fortran). C has no namelists, so there is a
 * small reader here that understands what those files contain:
 *     &section
 *       name = value      ! comment
 *     /
 * with numbers (1.0, 1.0e-4, 1.0d-4), switches (.true. / .false.),
 * text in quotes, and  tij(1,1) = ...  for the tidal tensor. Sections
 * can be in any order; a parameter left out keeps its standard value.
 * A name that is not a parameter is an error, as it is in Fortran.
 * (Archived Fortran equivalent: readparams in nbsetup.f90.)
 * ============================================================== */
#include <ctype.h>
#include <omp.h>
#include "hendos.h"

const char *typename_[NTYPES + 1] = {"group", "star", "white_dwarf", "neutron_star",
                                     "black_hole", "planet"};

void die(const char *msg)
{
    printf(" ERROR: %s\n", msg);
    exit(1);
}

void *xalloc(size_t nbytes)
{
    void *p = calloc(1, nbytes ? nbytes : 1);
    if (!p) die("out of memory");
    return p;
}

void setunits(void)
{
    /* internal units G=1, M=1 Msun, R=1 pc: velocity unit 0.06558 km/s,
     * time unit 14.91 Myr */
    double gmsun = 1.32712440018e20;    /* G*Msun in m^3 s^-2 (IAU) */
    double pc = 3.0856775814913673e16;  /* m */
    double myr = 3.15576e13;            /* s (Julian) */
    vunit = sqrt(gmsun / pc) / 1.e3;
    tunit = pc / (vunit * 1.e3) / myr;
}

/* the name of member k of an ensemble: _e and the member number (3
 * digits, more if needed) before the file extension if there is one,
 * otherwise at the end. The same rule as ensname in the Fortran. */
void ensname(const char *name, int k, char *out)
{
    char num[32];
    const char *dot = strrchr(name, '.'), *slash = strrchr(name, '/');
    snprintf(num, sizeof num, k < 1000 ? "_e%03d" : "_e%d", k);
    if (dot && dot > name && (!slash || dot > slash + 1)) {
        int len = (int)(dot - name);
        sprintf(out, "%.*s%s%s", len, name, num, dot);
    } else {
        sprintf(out, "%s%s", name, num);
    }
}

/* ---- the table of parameters: section, name, kind, where it lives -- */
typedef struct {
    const char *sec, *name;
    char kind;          /* d real, i whole number, l switch, s text */
    void *ptr;
    int len;            /* room for text */
} Param;

static Param table[] = {
    {"run", "runid", 's', runid, sizeof runid},
    {"run", "icfile", 's', icfile, sizeof icfile},
    {"run", "tend", 'd', &tend, 0},
    {"run", "tout", 'd', &tout, 0},
    {"run", "eps", 'd', &eps, 0},
    {"run", "restart", 'l', &restart_, 0},
    {"run", "nensemble", 'i', &nensemble, 0},
    {"run", "nparallel", 'i', &nparallel, 0},
    {"run", "ensfirst", 'i', &ensfirst, 0},
    {"output", "datadir", 's', datadir, sizeof datadir},
    {"output", "screen", 's', screen, sizeof screen},
    {"output", "tsnap", 'd', &tsnap, 0},
    {"output", "dewarn", 'd', &dewarn, 0},
    {"output", "trestart", 'd', &trestart, 0},
    {"analysis", "tana", 'd', &tana, 0},
    {"analysis", "lagrangian", 'l', &lagrangian, 0},
    {"tides", "tidal", 's', tidal, sizeof tidal},
    {"tides", "oorta", 'd', &oorta, 0},
    {"tides", "oortb", 'd', &oortb, 0},
    {"tides", "rho0", 'd', &rho0, 0},
    {"tides", "omegat", 'd', &omegat, 0},
    {"tides", "tij", 't', tij, 0},
    {"integrator", "eta", 'd', &eta, 0},
    {"integrator", "etap", 'd', &etap, 0},
    {"integrator", "etas", 'd', &etas, 0},
    {"integrator", "niter", 'i', &niter, 0},
    {"integrator", "kmax", 'i', &kmax, 0},
    {"binaries", "closebinaries", 'l', &closebinaries, 0},
    {"binaries", "dtclose", 'd', &dtclose, 0},
    {"binaries", "gclose", 'd', &gclose, 0},
    {"binaries", "gopen", 'd', &gopen, 0},
    {"binaries", "fewbody", 'l', &fewbody, 0},
    {"binaries", "gfew", 'd', &gfew, 0},
    {"binaries", "etafew", 'd', &etafew, 0},
    {"collisions", "rcoll", 'd', &rcoll, 0},
    {"collisions", "mergers", 'l', &mergers, 0},
    {"neighbours", "useac", 'l', &useac, 0},
    {"neighbours", "nnbopt", 'i', &nnbopt, 0},
    {"neighbours", "nnbmax", 'i', &nnbmax, 0},
    {"neighbours", "etar", 'd', &etar, 0},
    {"tree", "usetree", 'l', &usetree, 0},
    {"tree", "theta", 'd', &theta, 0},
    {"tree", "nleaf", 'i', &nleaf, 0},
    {"tree", "ntreebuild", 'i', &ntreebuild, 0},
    {"computing", "nompmin", 'i', &nompmin, 0},
    {"computing", "nompj", 'i', &nompj, 0},
    {NULL, NULL, 0, NULL, 0}};

static void lower(char *s)
{
    for (; *s; s++) *s = (char)tolower((unsigned char)*s);
}

static double getreal(const char *tok, const char *sec, const char *key)
{
    char buf[128], *end;
    snprintf(buf, sizeof buf, "%s", tok);
    for (char *p = buf; *p; p++)
        if (*p == 'd' || *p == 'D') *p = 'e';     /* Fortran 1.0d-4 */
    double x = strtod(buf, &end);
    if (end == buf || *end) {
        printf(" ERROR in &%s section of %s:\n   %s = %s is not a number\n", sec,
               pfilename, key, tok);
        exit(1);
    }
    return x;
}

static int getlogical(const char *tok, const char *sec, const char *key)
{
    char buf[32];
    snprintf(buf, sizeof buf, "%s", tok);
    lower(buf);
    if (!strcmp(buf, ".true.") || !strcmp(buf, "t") || !strcmp(buf, "true")) return 1;
    if (!strcmp(buf, ".false.") || !strcmp(buf, "f") || !strcmp(buf, "false")) return 0;
    printf(" ERROR in &%s section of %s:\n   %s = %s is not .true. or .false.\n", sec,
           pfilename, key, tok);
    exit(1);
}

/* set one parameter from its values (nv of them) */
static void setparam(const char *sec, char *key, char **val, int nv)
{
    int i1 = 0, i2 = 0;
    char name[64];
    snprintf(name, sizeof name, "%s", key);
    char *paren = strchr(name, '(');
    if (paren) {                                   /* tij(1,1) */
        sscanf(paren, "(%d,%d)", &i1, &i2);
        *paren = 0;
    }
    lower(name);
    for (Param *p = table; p->sec; p++) {
        if (strcmp(p->sec, sec) || strcmp(p->name, name)) continue;
        if (nv < 1) {
            printf(" ERROR in &%s section of %s:\n   %s has no value\n", sec, pfilename, key);
            exit(1);
        }
        switch (p->kind) {
        case 'd': *(double *)p->ptr = getreal(val[0], sec, key); break;
        case 'i': *(int *)p->ptr = (int)lround(getreal(val[0], sec, key)); break;
        case 'l': *(int *)p->ptr = getlogical(val[0], sec, key); break;
        case 's': snprintf((char *)p->ptr, (size_t)p->len, "%s", val[0]); break;
        case 't':
            if (paren) {
                if (i1 < 1 || i1 > 3 || i2 < 1 || i2 > 3) die("tij index must be 1 to 3");
                tij[i1 - 1][i2 - 1] = getreal(val[0], sec, key);
            } else {
                /* tij = nine numbers, in Fortran order (first index fastest) */
                for (int k = 0; k < nv && k < 9; k++)
                    tij[k % 3][k / 3] = getreal(val[k], sec, key);
            }
            break;
        }
        return;
    }
    printf(" ERROR in &%s section of %s:\n   %s is not a parameter of this section\n", sec,
           pfilename, key);
    exit(1);
}

/* read the whole file, find each &section ... / and set its parameters.
 * Returns 1 if the section 'run' was there. */
static int parsefile(const char *pfile)
{
    FILE *f = fopen(pfile, "r");
    if (!f) {
        printf(" ERROR: cannot open parameter file %s\n", pfile);
        exit(1);
    }
    fseek(f, 0, SEEK_END);
    long size = ftell(f);
    rewind(f);
    char *text = xalloc((size_t)size + 2);
    if (fread(text, 1, (size_t)size, f) != (size_t)size) die("cannot read the parameter file");
    fclose(f);

    /* tokens: text in quotes (kept with its opening quote as a marker),
     * '=', '/', or anything up to a space, comma, '=', '/' or '!'.
     * '!' starts a comment that runs to the end of the line. */
    int maxtok = (int)size + 8, ntok = 0, gotrun = 0;
    char **tok = xalloc(sizeof(char *) * (size_t)maxtok);
    char *p = text;
    while (*p) {
        char *start;
        size_t len;
        if (isspace((unsigned char)*p) || *p == ',') { p++; continue; }
        if (*p == '!') { while (*p && *p != '\n') p++; continue; }
        if (*p == '\'' || *p == '"') {
            char q = *p;
            start = p++;
            while (*p && *p != q) p++;
            len = (size_t)(p - start);
            if (*p) p++;
        } else if (*p == '=' || *p == '/') {
            start = p++;
            len = 1;
        } else {
            start = p;
            while (*p && !isspace((unsigned char)*p) && *p != ',' && *p != '=' && *p != '!'
                   && *p != '/') {
                if (*p == '(')                     /* tij(1,1) is one token */
                    while (*p && *p != ')') p++;
                if (*p) p++;
            }
            len = (size_t)(p - start);
        }
        tok[ntok] = xalloc(len + 1);
        memcpy(tok[ntok], start, len);
        ntok++;
    }
    char sec[64] = "";
    for (int t = 0; t < ntok;) {
        if (tok[t][0] == '&') {
            snprintf(sec, sizeof sec, "%s", tok[t] + 1);
            lower(sec);
            if (!strcmp(sec, "run")) gotrun = 1;
            t++;
        } else if (!strcmp(tok[t], "/")) {
            sec[0] = 0;
            t++;
        } else if (!sec[0]) {
            t++;                                   /* text outside any section */
        } else {
            char *key = tok[t];
            if (t + 1 >= ntok || strcmp(tok[t + 1], "=")) {
                printf(" ERROR in &%s section of %s:\n   expected  name = value  at '%s'\n",
                       sec, pfile, key);
                exit(1);
            }
            int t2 = t + 2, nv = 0;
            char *vals[16];
            while (t2 < ntok && strcmp(tok[t2], "/") && tok[t2][0] != '&'
                   && !(t2 + 1 < ntok && !strcmp(tok[t2 + 1], "="))) {
                if (nv < 16)
                    vals[nv++] = (tok[t2][0] == '\'' || tok[t2][0] == '"') ? tok[t2] + 1
                                                                           : tok[t2];
                t2++;
            }
            setparam(sec, key, vals, nv);
            t = t2;
        }
    }
    for (int t = 0; t < ntok; t++) free(tok[t]);
    free(tok);
    free(text);
    return gotrun;
}

static int safename(const char *s)
{
    for (; *s; s++)
        if (!isalnum((unsigned char)*s) && *s != '-' && *s != '_' && *s != '.') return 0;
    return 1;
}

void readparams(const char *pfile)
{
    int bad = 0;
    /* standard values -------------------------------------------- */
    runid[0] = icfile[0] = 0;
    tend = 1.0; tout = 0.1; eps = 0.0;
    strcpy(datadir, "data"); strcpy(screen, "standard");
    tsnap = 0.0; dewarn = 1.e-4; tana = 0.0; lagrangian = 1;
    strcpy(tidal, "none");
    oorta = 15.3; oortb = -11.9; rho0 = 0.097; omegat = 0.0;   /* Bovy 2017; McKee+ 2015 */
    memset(tij, 0, sizeof tij);
    eta = 0.01; etap = 0.2; etas = 0.01; niter = 1; kmax = 40;
    closebinaries = 1; fewbody = 1; gfew = 1.e-5; etafew = 0.02;
    rcoll = 10.0; mergers = 1; restart_ = 0; trestart = 0.0;
    nensemble = 1; nparallel = 1; ensfirst = 1;
    useac = 0; nnbopt = 0; nnbmax = 0; etar = 0.01;
    usetree = 0; theta = 0.5; nleaf = 8; ntreebuild = 10;
    dtclose = 100.0; gclose = 1.e-6; gopen = 1.e-5;
    nompmin = 32; nompj = 1024;

    if (!parsefile(pfile)) {
        printf(" ERROR: no &run section found in %s\n", pfile);
        exit(1);
    }
    nthreads = omp_get_max_threads();
    if (nthreads < 1) nthreads = 1;

    /* data directory: held with one trailing /, blank = the run directory */
    if (datadir[0] == '~' && datadir[1] == '/') {
        char tmp[NSTR];
        snprintf(tmp, sizeof tmp, "%s%s", getenv("HOME") ? getenv("HOME") : "", datadir + 1);
        strcpy(datadir, tmp);
    }
    size_t l = strlen(datadir);
    while (l > 0 && datadir[l - 1] == '/') datadir[--l] = 0;
    if (!strcmp(datadir, ".")) datadir[0] = 0;
    if (datadir[0]) strcat(datadir, "/");

    /* checks ------------------------------------------------------ */
#define BAD(msg) do { printf(" ERROR: %s\n", msg); bad = 1; } while (0)
    if (strchr(datadir, ' ')) BAD("datadir name must not contain spaces");
    if (!runid[0]) BAD("runid must be set in &run - it is used in all output file names");
    else if (!safename(runid)) BAD("runid may only contain letters, digits, - _ and . (no spaces)");
    if (!icfile[0]) BAD("icfile must be set in &run");
    else if (strchr(icfile, ' ')) BAD("icfile name must not contain spaces");
    if (tout <= 0.0) BAD("tout must be > 0");
    if (eps < 0.0) BAD("eps must be >= 0");
    if (!strcmp(screen, "verbose")) lverbose = 1;
    else if (!strcmp(screen, "standard")) lverbose = 0;
    else BAD("screen must be 'standard' or 'verbose'");
    if (tsnap < 0.0 || tana < 0.0) BAD("tsnap and tana must be >= 0 (0 = off)");
    if (eta <= 0.0 || etap <= 0.0 || etas <= 0.0) BAD("eta, etap and etas must be > 0");
    if (niter < 1) BAD("niter must be >= 1");
    if (kmax < 1 || kmax > 60) BAD("kmax must be between 1 and 60");

    /* tides - to internal units, 1 km/s/kpc = fconv internal.
     * In a frame rotating with the cluster's circular orbit the extra
     * acceleration is a = T.x - 2 Omega x v. 'solar' builds T from the
     * Oort constants and local density:
     *   Omega = A - B, T = diag(4A(A-B), 0, -(4 pi G rho0 + 2(A^2-B^2))) */
    double fconv = 1.0 / (vunit * 1.e3), pi = 3.141592653589793;
    memset(ttid, 0, sizeof ttid);
    omega = 0.0;
    ltide = 0;
    if (!strcmp(tidal, "none")) {
        ltide = 0;
    } else if (!strcmp(tidal, "solar")) {
        ltide = 1;
        if (rho0 < 0.0) BAD("rho0 must be >= 0");
        omega = (oorta - oortb) * fconv;
        ttid[0][0] = 4.0 * oorta * (oorta - oortb) * fconv * fconv;
        ttid[2][2] = -(4.0 * pi * rho0 + 2.0 * (oorta * oorta - oortb * oortb) * fconv * fconv);
    } else if (!strcmp(tidal, "tensor")) {
        ltide = 1;
        omega = omegat * fconv;
        double asym = 0.0, big = 0.0;
        for (int k = 0; k < 3; k++)
            for (int j = 0; j < 3; j++) {
                ttid[k][j] = 0.5 * (tij[k][j] + tij[j][k]) * fconv * fconv;
                asym = fmax(asym, fabs(tij[k][j] - tij[j][k]));
                big = fmax(big, fabs(tij[k][j]));
            }
        if (asym > 1.e-10 * big)
            printf(" WARNING: tij is not symmetric - using its symmetric part\n");
    } else {
        BAD("tidal must be 'none', 'solar' or 'tensor'");
    }
    if (closebinaries) {
        if (gclose <= 0.0 || gopen <= gclose) BAD("need 0 < gclose < gopen");
        if (dtclose <= 0.0) BAD("dtclose must be > 0");
        if (eps > 0.0) {
            printf(" WARNING: binary treatment switched off because eps > 0 "
                   "(Kepler orbits are unsoftened)\n");
            closebinaries = 0;
        }
    }
    if (!closebinaries) fewbody = 0;
    if (fewbody) {
        if (gfew <= 0.0) BAD("gfew must be > 0");
        if (etafew <= 0.0 || etafew > 0.2) BAD("etafew must be between 0 and 0.2");
    }
    if (usetree) BAD("the Barnes-Hut tree is not in this version of HENDOS - "
                     "set usetree = .false. (only the archived Fortran version, "
                     "archive/fortran, has it)");
    if (useac) {
        if (nnbopt < 0 || nnbmax < 0) BAD("nnbopt and nnbmax must be >= 0 (0 = automatic)");
        if (etar <= 0.0) BAD("etar must be > 0");
    }
    if (nensemble < 1 || nparallel < 1 || nparallel > 64)
        BAD("nensemble must be at least 1 and nparallel between 1 and 64");
    if (ensfirst < 1 || ensfirst > nensemble) BAD("ensfirst must be between 1 and nensemble");
    if (nensemble > 1 && strlen(runid) > 26)
        BAD("runid must be 26 characters or fewer for an ensemble");
    if (trestart < 0.0) BAD("trestart must be >= 0 (0 = no restart files)");
    dtclosei = dtclose * 1.e-6 / tunit;
    if (rcoll < 0.0) BAD("rcoll must be >= 0 (0 = no collision check)");
    rcolli = rcoll * 6.957e8 / 3.0856775814913673e16;      /* Rsun in pc */
    if (bad) exit(1);
#undef BAD

    if (tsnap > 0.0 && (fabs(tsnap / tout - round(tsnap / tout)) > 1.e-6 || tsnap < tout))
        printf(" WARNING: tsnap is not a multiple of tout, snapshots every %.4f Myr\n",
               fmax(1.0, round(tsnap / tout)) * tout);
    if (tana > 0.0 && (fabs(tana / tout - round(tana / tout)) > 1.e-6 || tana < tout))
        printf(" WARNING: tana is not a multiple of tout, analysis every %.4f Myr\n",
               fmax(1.0, round(tana / tout)) * tout);
}
