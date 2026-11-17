/* ==============================================================
 * HENDOS - HErmite N-body Dynamics Of Systems - tree.c
 * BARNES-HUT TREE (switch usetree) for large associations
 *
 * Octree over the bodies' predicted positions rp, with up to nleaf
 * bodies per leaf. Each node has monopole, quadrupole and octupole
 * moments about its centre of mass (raw moments S2_ij = sum m x_i x_j,
 * S3_ijk = sum m x_i x_j x_k, x relative to the node c.o.m.), its c.o.m.
 * velocity, and bmax, the radius about the c.o.m. containing its bodies.
 *
 * Opening criterion: a node is used as a whole if bmax/d < theta, d the
 * distance from the field star to the node c.o.m.; otherwise it is
 * opened. Leaves that are opened are summed directly (exactly as in the
 * direct code, including softening).
 * Acceleration from an accepted node (r from node c.o.m. to the star):
 *   a = -M r/r^3
 *       + (3 S2.r - s r)/r^5 - 5 (3 q - r^2 s) r / (2 r^7)
 *       + (15 C - 6 r (T.r) - 3 r^2 T)/(2 r^7) - 7 (15 c - 9 r^2 T.r) r/(6 r^9)
 *   s = tr S2, q = r.S2.r, C_i = S3_ijk r_j r_k, c = C.r, T_k = S3_iik
 * Jerk from an accepted node: monopole (with the node c.o.m. velocity)
 * plus the quadrupole term from the internal motion of the node's
 * bodies, i.e. the quadrupole formula with S2 replaced by
 * dS2/dt = sum m (x u + u x), u velocities relative to the node c.o.m.
 * (tested: roughly halves the jerk error of the monopole alone).
 *
 * Nodes are numbered in depth-first order from 0 (the root), so the
 * first child of an internal node is node+1, and tnext[node] is the
 * next node once the node's subtree is finished - the tree is walked
 * without a stack. The bodies of each node are a contiguous range of
 * tperm.
 *
 * The tree only needs to be built every so often; it is restocked every
 * step. The topology is rebuilt (treebuild) every ntreebuild block steps
 * (parameter, standard value 10), at the start of every output interval
 * (so that a restarted run is exactly the same as an unbroken one), and
 * after any group is closed or opened or stars merge. On every block
 * step in between only the moments and bmax are recomputed from the new
 * predicted positions (treemoments, O(N)). As bmax is recomputed, the
 * opening criterion stays valid as bodies move: an old tree is never
 * wrong, only slower, because its boxes overlap more and more nodes
 * have to be opened.
 *
 * OpenMP: the build is shared out at the top level. The root's eight
 * octants are independent, so each is given to a thread (2, 4 or 8
 * threads share them evenly; other numbers work, more than 8 are not
 * used): the thread looks through all the particles for those in its
 * octant and builds that octant's whole subtree. Each octant builds into
 * its own region of the node arrays, and the regions are then closed up
 * so that the numbering is exactly that of a serial depth-first build.
 * The moments are computed in parallel over the leaves, then over the
 * same top-level subtrees for the internal nodes. The force walks are
 * shared out over the active stars (in force.c).
 *
 * The Galactic tide is added to every star after the force sum, so it is
 * included whether the tree is used or not. Energies are always direct.
 * The perturbation of groups (gperturb in groups.c) uses the tree when it
 * is up to date; the nearest-neighbour searches do not, because they
 * decide collisions and must not miss a star that has just moved.
 * (Archived Fortran equivalent: nbtree.f90.)
 * ============================================================== */
#include <omp.h>
#include "hendos.h"

#define MAXDEPTH 40
#define NSPARE 125              /* spare nodes per top-level octant */

static int nnodemax;
static int *tperm, *toct, *tfirst, *tcount, *tnext;
static char *tleaf;
static double *tm, (*tcom)[3], (*tcv)[3], *tbmax, (*ts2)[9], (*ts3)[10], *ttr, (*tt3)[3];
static double (*tsd)[9], *tsdtr;

/* node arrays: the root plus 3 nodes per body and NSPARE for each of the
 * 8 top-level octants is generous; overflow stops with a message */
void treealloc(void)
{
    size_t nm;
    nnodemax = 3 * nalloc + 8 * NSPARE + 1;
    nm = (size_t)nnodemax;
    tperm = xalloc(sizeof(int) * (size_t)nalloc);
    toct = xalloc(sizeof(int) * (size_t)nalloc);
    tfirst = xalloc(sizeof(int) * nm);
    tcount = xalloc(sizeof(int) * nm);
    tnext = xalloc(sizeof(int) * nm);
    tleaf = xalloc(nm);
    tm = xalloc(sizeof(double) * nm);
    tcom = xalloc(sizeof(*tcom) * nm);
    tcv = xalloc(sizeof(*tcv) * nm);
    tbmax = xalloc(sizeof(double) * nm);
    ts2 = xalloc(sizeof(*ts2) * nm);
    ts3 = xalloc(sizeof(*ts3) * nm);
    ttr = xalloc(sizeof(double) * nm);
    tt3 = xalloc(sizeof(*tt3) * nm);
    tsd = xalloc(sizeof(*tsd) * nm);
    tsdtr = xalloc(sizeof(double) * nm);
    nnode = 0;
    treevalid = 0;
}

static inline int octant(const double x[3], const double c[3])
{
    return (x[0] >= c[0]) + 2 * (x[1] >= c[1]) + 4 * (x[2] >= c[2]);
}

static void childcentre(const double c[3], double half, int oct, double cc[3])
{
    cc[0] = c[0] + ((oct & 1) ? 0.5 * half : -0.5 * half);
    cc[1] = c[1] + ((oct & 2) ? 0.5 * half : -0.5 * half);
    cc[2] = c[2] + ((oct & 4) ? 0.5 * half : -0.5 * half);
}

/* node holds bodies tperm[ifirst..ilast] in the cube centre +- half.
 * Split into octants in place; children are created depth first. *nn is
 * the last node number used so far in this subtree's region of the node
 * arrays and nlim the last one it may use. Thread safe: only this
 * subtree's nodes and its own part of tperm are written. */
static void buildnode(int node, int ifirst, int ilast, const double centre[3], double half,
                      int depth, int *nn, int nlim)
{
    int cnt[8] = {0}, start[8], pos[8];
    tfirst[node] = ifirst;
    tcount[node] = ilast - ifirst + 1;
    if (tcount[node] <= nleaf || depth >= MAXDEPTH) {
        tleaf[node] = 1;
        tnext[node] = *nn + 1;
        return;
    }
    tleaf[node] = 0;
    /* sort the node's bodies by octant (counting sort) */
    for (int k = ifirst; k <= ilast; k++) cnt[octant(rp[tperm[k]], centre)]++;
    start[0] = ifirst;
    for (int o = 1; o < 8; o++) start[o] = start[o - 1] + cnt[o - 1];
    memcpy(pos, start, sizeof pos);
    int *tmp = xalloc(sizeof(int) * (size_t)tcount[node]);
    for (int k = ifirst; k <= ilast; k++) {
        int ip = tperm[k];
        tmp[pos[octant(rp[ip], centre)]++ - ifirst] = ip;
    }
    memcpy(tperm + ifirst, tmp, sizeof(int) * (size_t)tcount[node]);
    free(tmp);
    /* children (non-empty octants only), depth first */
    for (int o = 0; o < 8; o++) {
        if (cnt[o] == 0) continue;
        if (++(*nn) > nlim) die("tree node overflow - try a larger nleaf");
        double cc[3];
        childcentre(centre, half, o, cc);
        buildnode(*nn, start[o], start[o] + cnt[o] - 1, cc, 0.5 * half, depth + 1, nn, nlim);
    }
    tnext[node] = *nn + 1;
}

/* add a point mass mp at x (relative to the c.o.m.) to the raw moments.
 * s2 is the full 3x3 matrix (row major), s3 packed:
 * xxx xxy xxz xyy xyz xzz yyy yyz yzz zzz */
static inline void addpoint(double mp, const double x[3], double s2[9], double s3[10])
{
    for (int i = 0; i < 3; i++)
        for (int j = 0; j < 3; j++) s2[3 * i + j] += mp * x[i] * x[j];
    s3[0] += mp * x[0] * x[0] * x[0];
    s3[1] += mp * x[0] * x[0] * x[1];
    s3[2] += mp * x[0] * x[0] * x[2];
    s3[3] += mp * x[0] * x[1] * x[1];
    s3[4] += mp * x[0] * x[1] * x[2];
    s3[5] += mp * x[0] * x[2] * x[2];
    s3[6] += mp * x[1] * x[1] * x[1];
    s3[7] += mp * x[1] * x[1] * x[2];
    s3[8] += mp * x[1] * x[2] * x[2];
    s3[9] += mp * x[2] * x[2] * x[2];
}

/* add a child's moments (about its own c.o.m., which is at d from the
 * parent c.o.m.) to the parent's: the parallel-axis theorem
 *   S2_ij  += S2c_ij + mc d_i d_j
 *   S3_ijk += S3c_ijk + d_i S2c_jk + d_j S2c_ik + d_k S2c_ij + mc d_i d_j d_k */
static void shiftadd(double mc, const double d[3], const double s2c[9], const double s3c[10],
                     double s2[9], double s3[10])
{
    static const int ia[10] = {0, 0, 0, 0, 0, 0, 1, 1, 1, 2};
    static const int ja[10] = {0, 0, 0, 1, 1, 2, 1, 1, 2, 2};
    static const int ka[10] = {0, 1, 2, 1, 2, 2, 1, 2, 2, 2};
    for (int i = 0; i < 3; i++)
        for (int j = 0; j < 3; j++) s2[3 * i + j] += s2c[3 * i + j] + mc * d[i] * d[j];
    for (int l = 0; l < 10; l++) {
        int i = ia[l], j = ja[l], k = ka[l];
        s3[l] += s3c[l] + d[i] * s2c[3 * j + k] + d[j] * s2c[3 * i + k] + d[k] * s2c[3 * i + j]
                 + mc * d[i] * d[j] * d[k];
    }
}

/* moments, c.o.m. velocity and bmax of internal node from its children,
 * which must already be done. Thread safe: only this node is written. */
static void nodemoments(int node)
{
    double mt = 0.0, b = 0.0;
    for (int k = 0; k < 3; k++) tcom[node][k] = tcv[node][k] = 0.0;
    for (int c = node + 1; c < tnext[node]; c = tnext[c]) {
        mt += tm[c];
        for (int k = 0; k < 3; k++) {
            tcom[node][k] += tm[c] * tcom[c][k];
            tcv[node][k] += tm[c] * tcv[c][k];
        }
    }
    tm[node] = mt;
    for (int k = 0; k < 3; k++) {
        tcom[node][k] /= mt;
        tcv[node][k] /= mt;
    }
    memset(ts2[node], 0, sizeof ts2[node]);
    memset(ts3[node], 0, sizeof ts3[node]);
    memset(tsd[node], 0, sizeof tsd[node]);
    for (int c = node + 1; c < tnext[node]; c = tnext[c]) {
        double d[3], w[3];
        for (int k = 0; k < 3; k++) {
            d[k] = tcom[c][k] - tcom[node][k];
            w[k] = tcv[c][k] - tcv[node][k];
        }
        shiftadd(tm[c], d, ts2[c], ts3[c], ts2[node], ts3[node]);
        for (int i = 0; i < 3; i++)
            for (int j = 0; j < 3; j++)
                tsd[node][3 * i + j] += tsd[c][3 * i + j] + tm[c] * (d[i] * w[j] + w[i] * d[j]);
        b = fmax(b, norm3(d) + tbmax[c]);
    }
    tbmax[node] = b;
}

/* moments, c.o.m. velocity and bmax of every node from the current rp,
 * vp: leaves directly from their bodies, internal nodes from their
 * children, deepest first (children always have larger node numbers
 * than their parent). The leaves are shared out over the threads; the
 * internal nodes of each top-level subtree only need that subtree, so
 * the (up to 8) subtrees are shared out too, and the root is done last. */
void treemoments(void)
{
    int nsub = 0, sfirst[8], slast[8];
#pragma omp parallel for schedule(dynamic, 64) if (n >= nompj)
    for (int node = 0; node < nnode; node++) {
        if (!tleaf[node]) continue;
        double mt = 0.0, b = 0.0, com[3] = {0, 0, 0}, cv[3] = {0, 0, 0};
        int k1 = tfirst[node], k2 = tfirst[node] + tcount[node];
        for (int k = k1; k < k2; k++) {
            int ip = tperm[k];
            mt += m[ip];
            for (int c = 0; c < 3; c++) {
                com[c] += m[ip] * rp[ip][c];
                cv[c] += m[ip] * vp[ip][c];
            }
        }
        for (int c = 0; c < 3; c++) {
            com[c] /= mt;
            cv[c] /= mt;
        }
        tm[node] = mt;
        memcpy(tcom[node], com, sizeof com);
        memcpy(tcv[node], cv, sizeof cv);
        memset(ts2[node], 0, sizeof ts2[node]);
        memset(ts3[node], 0, sizeof ts3[node]);
        memset(tsd[node], 0, sizeof tsd[node]);
        for (int k = k1; k < k2; k++) {
            int ip = tperm[k];
            double x[3], u[3];
            for (int c = 0; c < 3; c++) {
                x[c] = rp[ip][c] - com[c];
                u[c] = vp[ip][c] - cv[c];
            }
            addpoint(m[ip], x, ts2[node], ts3[node]);
            for (int i = 0; i < 3; i++)
                for (int j = 0; j < 3; j++) tsd[node][3 * i + j] += m[ip] * (x[i] * u[j] + u[i] * x[j]);
            b = fmax(b, dot3(x, x));
        }
        tbmax[node] = sqrt(b);
    }
    /* internal nodes: the top-level subtrees (the children of the root) */
    if (!tleaf[0])
        for (int c = 1; c < tnext[0]; c = tnext[c]) {
            sfirst[nsub] = c;
            slast[nsub++] = tnext[c] - 1;
        }
#pragma omp parallel for schedule(dynamic, 1) if (n >= nompj)
    for (int is = 0; is < nsub; is++)
        for (int node = slast[is]; node >= sfirst[is]; node--)
            if (!tleaf[node]) nodemoments(node);
    if (!tleaf[0]) nodemoments(0);
    /* derived quantities used in every evaluation: trace of S2 and of
     * dS2/dt, and T_k = S3_iik */
#pragma omp parallel for schedule(static) if (n >= nompj)
    for (int node = 0; node < nnode; node++) {
        ttr[node] = ts2[node][0] + ts2[node][4] + ts2[node][8];
        tsdtr[node] = tsd[node][0] + tsd[node][4] + tsd[node][8];
        tt3[node][0] = ts3[node][0] + ts3[node][3] + ts3[node][5];
        tt3[node][1] = ts3[node][1] + ts3[node][6] + ts3[node][8];
        tt3[node][2] = ts3[node][2] + ts3[node][7] + ts3[node][9];
    }
}

/* build the tree topology for bodies 0..n-1 at positions rp, then the
 * moments. The eight top-level octants are built in parallel (see the
 * top of this file). */
void treebuild(void)
{
    double xmin[3] = {HUGE_VAL, HUGE_VAL, HUGE_VAL}, xmax[3] = {-HUGE_VAL, -HUGE_VAL, -HUGE_VAL};
    double x0 = HUGE_VAL, y0 = HUGE_VAL, z0 = HUGE_VAL, x1 = -HUGE_VAL, y1 = -HUGE_VAL, z1 = -HUGE_VAL;
#pragma omp parallel for schedule(static) if (n >= nompj) reduction(min:x0,y0,z0) reduction(max:x1,y1,z1)
    for (int i = 0; i < n; i++) {
        x0 = fmin(x0, rp[i][0]); y0 = fmin(y0, rp[i][1]); z0 = fmin(z0, rp[i][2]);
        x1 = fmax(x1, rp[i][0]); y1 = fmax(y1, rp[i][1]); z1 = fmax(z1, rp[i][2]);
    }
    xmin[0] = x0; xmin[1] = y0; xmin[2] = z0;
    xmax[0] = x1; xmax[1] = y1; xmax[2] = z1;
    double centre[3], half = 0.0;
    for (int k = 0; k < 3; k++) {
        centre[k] = 0.5 * (xmin[k] + xmax[k]);
        half = fmax(half, xmax[k] - xmin[k]);
    }
    half = 0.5 * half * 1.0001 + 1e-300;

    /* the root (node 0) holds everything */
    tfirst[0] = 0;
    tcount[0] = n;
    if (n <= nleaf) {
        /* so few bodies that the root is a leaf */
        for (int i = 0; i < n; i++) tperm[i] = i;
        tleaf[0] = 1;
        nnode = 1;
    } else {
        int cnt[8] = {0}, start[8], base[8], last[8];
        int c0 = 0, c1 = 0, c2 = 0, c3 = 0, c4 = 0, c5 = 0, c6 = 0, c7 = 0;
        tleaf[0] = 0;
        /* top-level octant of every body, and the number in each octant */
#pragma omp parallel for schedule(static) if (n >= nompj) reduction(+:c0,c1,c2,c3,c4,c5,c6,c7)
        for (int i = 0; i < n; i++) {
            int o = octant(rp[i], centre);
            toct[i] = o;
            switch (o) {
            case 0: c0++; break; case 1: c1++; break; case 2: c2++; break;
            case 3: c3++; break; case 4: c4++; break; case 5: c5++; break;
            case 6: c6++; break; default: c7++; break;
            }
        }
        cnt[0] = c0; cnt[1] = c1; cnt[2] = c2; cnt[3] = c3;
        cnt[4] = c4; cnt[5] = c5; cnt[6] = c6; cnt[7] = c7;
        /* each octant's bodies are tperm[start..start+cnt-1], and its
         * subtree is built in its own region of the node arrays */
        start[0] = 0;
        for (int o = 1; o < 8; o++) start[o] = start[o - 1] + cnt[o - 1];
        for (int o = 0; o < 8; o++) base[o] = 1 + 3 * start[o] + NSPARE * o;
        /* one octant per thread (up to 8 threads): look through all the
         * bodies for those in this octant, then build its subtree.
         * Nothing is shared between octants, so there are no races. */
        int nt8 = nthreads < 8 ? nthreads : 8;
#pragma omp parallel for schedule(dynamic, 1) num_threads(nt8) if (n >= nompj)
        for (int o = 0; o < 8; o++) {
            last[o] = base[o] - 1;
            if (cnt[o] == 0) continue;
            int k = start[o];
            for (int i = 0; i < n; i++)
                if (toct[i] == o) tperm[k++] = i;
            double cc[3];
            childcentre(centre, half, o, cc);
            int nn = base[o];
            buildnode(base[o], start[o], start[o] + cnt[o] - 1, cc, 0.5 * half, 1, &nn,
                      base[o] + 3 * cnt[o] + NSPARE - 1);
            last[o] = nn;
        }
        /* close up the gaps between the regions, in octant order, so the
         * nodes are numbered depth first without gaps. A region only ever
         * moves down, so copying forwards is safe. Node numbers inside a
         * region (tnext) move with it. */
        nnode = 1;
        for (int o = 0; o < 8; o++) {
            if (cnt[o] == 0) continue;
            int nused = last[o] - base[o] + 1, delta = base[o] - nnode;
            if (delta > 0)
                for (int k = nnode; k < nnode + nused; k++) {
                    tfirst[k] = tfirst[k + delta];
                    tcount[k] = tcount[k + delta];
                    tleaf[k] = tleaf[k + delta];
                    tnext[k] = tnext[k + delta] - delta;
                }
            nnode += nused;
        }
    }
    tnext[0] = nnode;
    treevalid = 1;
    ntreebuilt++;
    treemoments();
}

/* C_i = S3_ijk r_j r_k from the packed octupole moments */
static inline void octc(const double s3[10], const double r[3], double cv[3])
{
    double xx = r[0] * r[0], yy = r[1] * r[1], zz = r[2] * r[2];
    double xy = r[0] * r[1], xz = r[0] * r[2], yz = r[1] * r[2];
    cv[0] = s3[0] * xx + 2 * s3[1] * xy + 2 * s3[2] * xz + s3[3] * yy + 2 * s3[4] * yz + s3[5] * zz;
    cv[1] = s3[1] * xx + 2 * s3[3] * xy + 2 * s3[4] * xz + s3[6] * yy + 2 * s3[7] * yz + s3[8] * zz;
    cv[2] = s3[2] * xx + 2 * s3[4] * xy + 2 * s3[5] * xz + s3[7] * yy + 2 * s3[8] * yz + s3[9] * zz;
}

static inline void matvec(const double s[9], const double r[3], double out[3])
{
    for (int i = 0; i < 3; i++) out[i] = s[3 * i] * r[0] + s[3 * i + 1] * r[1] + s[3 * i + 2] * r[2];
}

/* acceleration, jerk, pairwise timescale and nearest distance for body i
 * at rp[i], vp[i] from the tree. Thread safe (results in arguments). */
void treeforce(int i, double a[3], double jk[3], double *tmin2, double *rm2)
{
    const double *xi = rp[i], *vi = vp[i];
    double mi = m[i], theta2 = theta * theta, t2 = HUGE_VAL, r2 = HUGE_VAL;
    double ax = 0, ay = 0, az = 0, jx = 0, jy = 0, jz = 0;
    int node = 0;
    while (node < nnode) {
        double rr[3] = {xi[0] - tcom[node][0], xi[1] - tcom[node][1], xi[2] - tcom[node][2]};
        double d2 = dot3(rr, rr);
        if (tbmax[node] * tbmax[node] < theta2 * d2) {
            /* ---- accepted node: multipoles for a, monopole + internal
             * motion quadrupole for the jerk */
            double d = sqrt(d2), rinv = 1.0 / d, rinv2 = rinv * rinv;
            double r5 = rinv2 * rinv2 * rinv, r7 = r5 * rinv2, r9 = r7 * rinv2;
            double sr[3], cv[3], s = ttr[node], q, c, tr, f;
            matvec(ts2[node], rr, sr);
            q = dot3(rr, sr);
            octc(ts3[node], rr, cv);
            c = dot3(cv, rr);
            tr = dot3(tt3[node], rr);
            double mono = -tm[node] * rinv2 * rinv, quad = -2.5 * (3.0 * q - d2 * s) * r7;
            double oct2 = -7.0 / 6.0 * (15.0 * c - 9.0 * d2 * tr) * r9;
            double acc[3];
            for (int k = 0; k < 3; k++)
                acc[k] = mono * rr[k] + (3.0 * sr[k] - s * rr[k]) * r5 + quad * rr[k]
                         + 0.5 * (15.0 * cv[k] - 6.0 * rr[k] * tr - 3.0 * d2 * tt3[node][k]) * r7
                         + oct2 * rr[k];
            ax += acc[0]; ay += acc[1]; az += acc[2];
            double dv[3] = {vi[0] - tcv[node][0], vi[1] - tcv[node][1], vi[2] - tcv[node][2]};
            double mor = tm[node] * rinv2 * rinv, drdv = 3.0 * dot3(rr, dv) * rinv2;
            jx -= mor * (dv[0] - drdv * rr[0]);
            jy -= mor * (dv[1] - drdv * rr[1]);
            jz -= mor * (dv[2] - drdv * rr[2]);
            /* internal motion of the node: quadrupole formula with dS2/dt */
            s = tsdtr[node];
            matvec(tsd[node], rr, sr);
            q = dot3(rr, sr);
            f = -2.5 * (3.0 * q - d2 * s) * r7;
            jx += (3.0 * sr[0] - s * rr[0]) * r5 + f * rr[0];
            jy += (3.0 * sr[1] - s * rr[1]) * r5 + f * rr[1];
            jz += (3.0 * sr[2] - s * rr[2]) * r5 + f * rr[2];
            double dv2 = dot3(dv, dv), tff = d2 * d / (mi + tm[node]), tfl = d2 / (dv2 + 1e-300);
            if (tff < t2) t2 = tff;
            if (tfl < t2) t2 = tfl;
            double dm = d - tbmax[node];
            if (dm * dm < r2) r2 = dm * dm;
            node = tnext[node];
        } else if (tleaf[node]) {
            /* ---- opened leaf: direct sum, exactly as in the direct code
             * (keep the same as sumrange in force.c) */
            int k2 = tfirst[node] + tcount[node];
            for (int k = tfirst[node]; k < k2; k++) {
                int j = tperm[k];
                if (j == i) continue;
                double dx = rp[j][0] - xi[0], dy = rp[j][1] - xi[1], dz = rp[j][2] - xi[2];
                double dvx = vp[j][0] - vi[0], dvy = vp[j][1] - vi[1], dvz = vp[j][2] - vi[2];
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
            node = tnext[node];
        } else {
            node++;                         /* open the node: its first child */
        }
    }
    a[0] = ax; a[1] = ay; a[2] = az;
    jk[0] = jx; jk[1] = jy; jk[2] = jz;
    *tmin2 = t2;
    *rm2 = r2;
}

/* sum over bodies (except those in ex[0..nex-1]) of m/d^3 at xc, and the
 * squared shortest pairwise free-fall/flyby time, using the monopoles of
 * accepted nodes - for the perturbation gamma of a group. The excluded
 * bodies are next to xc, so any node holding one of them is opened. */
void treeperturb(const double xc[3], const double vc[3], double mb, int nex, const int *ex,
                 double *sum3, double *tp2)
{
    double s = 0.0, t2 = HUGE_VAL, theta2 = theta * theta;
    int node = 0;
    while (node < nnode) {
        double rr[3] = {tcom[node][0] - xc[0], tcom[node][1] - xc[1], tcom[node][2] - xc[2]};
        double d2 = dot3(rr, rr);
        if (tbmax[node] * tbmax[node] < theta2 * d2) {
            double dv[3] = {tcv[node][0] - vc[0], tcv[node][1] - vc[1], tcv[node][2] - vc[2]};
            double d3 = d2 * sqrt(d2);
            s += tm[node] / d3;
            t2 = fmin(t2, fmin(d3 / (mb + tm[node]), d2 / (dot3(dv, dv) + 1e-300)));
            node = tnext[node];
        } else if (tleaf[node]) {
            int k2 = tfirst[node] + tcount[node];
            for (int k = tfirst[node]; k < k2; k++) {
                int j = tperm[k], skip = 0;
                for (int e = 0; e < nex; e++)
                    if (ex[e] == j) skip = 1;
                if (skip) continue;
                double r2v[3] = {rp[j][0] - xc[0], rp[j][1] - xc[1], rp[j][2] - xc[2]};
                double dv[3] = {vp[j][0] - vc[0], vp[j][1] - vc[1], vp[j][2] - vc[2]};
                double dd2 = dot3(r2v, r2v), d3 = dd2 * sqrt(dd2);
                s += m[j] / d3;
                t2 = fmin(t2, fmin(d3 / (mb + m[j]), dd2 / (dot3(dv, dv) + 1e-300)));
            }
            node = tnext[node];
        } else {
            node++;
        }
    }
    *sum3 = s;
    *tp2 = t2;
}
