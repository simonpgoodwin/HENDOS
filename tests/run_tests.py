#!/usr/bin/env python3
"""
HENDOS - HErmite N-body Dynamics Of Systems
Regression tests.

Runs the hermite executable on a set of small problems with known answers
and checks each result against a tolerance. Run it after every change to
the code, before trusting a production run:

    cd HENDOS/tests
    python3 run_tests.py                 all the tests
    python3 run_tests.py --list          what the tests are
    python3 run_tests.py --only merger   tests whose name contains 'merger'
    python3 run_tests.py --keep          keep the working files (in ./work)
    python3 run_tests.py --dry           write the input files, run nothing

It needs numpy and the executable HENDOS/run/hermite (make it first; another
one can be given with --hermite). Each test makes its own initial conditions
here, so makeics is not needed. The exit status is 0 if every test passed.

TO ADD A TEST: write a function test_<name>(t) below. It gets a Tester, makes
initial conditions with t.write_ics, runs the code with t.run (giving any
parameters that differ from DEFAULTS) and reports numbers with t.check.
It is found and run automatically. See tests/README.md.
"""
import argparse
import glob
import os
import re
import shutil
import subprocess
import sys
import time

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, '..', 'analysis'))
from hendosio import read_snapshot, read_header  # noqa: E402

# the units of the code: G = 1 with Msun and pc
GMSUN = 1.32712440018e20        # m^3 s^-2
PC = 3.0856775814913673e16      # m
MYR = 3.15576e13                # s
VUNIT = np.sqrt(GMSUN / PC) / 1e3          # km/s
TUNIT = PC / (VUNIT * 1e3) / MYR           # Myr
G = VUNIT ** 2                             # pc (km/s)^2 / Msun

# parameters every test starts from; a test changes what it needs to.
# Each key is 'section.name'.
DEFAULTS = {
    'run.tend': 1.0, 'run.tout': 0.25, 'run.eps': 0.0, 'run.restart': False,
    'output.datadir': '.', 'output.screen': 'verbose', 'output.tsnap': 0.25,
    'output.trestart': 0.0,
    'analysis.tana': 0.0,
    'tides.tidal': 'none',
    'tree.usetree': False,
    'neighbours.useac': False,
    'binaries.closebinaries': True, 'binaries.fewbody': True,
    'collisions.rcoll': 0.0, 'collisions.mergers': True,
    'integrator.eta': 0.02, 'integrator.niter': 1,
}


def fortran(value):
    """A Python value as it is written in a namelist."""
    if isinstance(value, bool):
        return '.true.' if value else '.false.'
    if isinstance(value, str):
        return "'%s'" % value
    if isinstance(value, float):
        return repr(value).replace('e', 'd') if 'e' in repr(value) else repr(value)
    return str(value)


class Failed(Exception):
    pass


class Tester:
    """What a test function uses: makes files, runs the code, checks numbers."""

    def __init__(self, name, hermite, workdir, dry):
        self.name = name
        self.hermite = hermite
        self.dir = os.path.join(workdir, name)
        self.dry = dry
        self.lines = []          # one line per check
        self.ok = True
        if os.path.isdir(self.dir):
            shutil.rmtree(self.dir)
        os.makedirs(self.dir)

    # -- input files ---------------------------------------------------
    def write_ics(self, fname, m, pos, vel, types=None):
        """Initial conditions in snapshot format (Msun, pc, km/s)."""
        m = np.asarray(m, float)
        n = len(m)
        with open(os.path.join(self.dir, fname), 'w') as f:
            f.write('# format                        = NBODY-SNAPSHOT 2\n')
            f.write('# creator                       = HENDOS tests/run_tests.py\n')
            f.write('# n                             = %d\n' % n)
            f.write('# time_Myr                      = 0.0\n')
            f.write('# columns                       = id m x y z vx vy vz type\n')
            f.write('# end_header\n')
            for i in range(n):
                typ = 'star' if types is None else types[i]
                f.write('%10d %24.16e %24.16e %24.16e %24.16e %24.16e %24.16e %24.16e %s\n'
                        % ((i + 1, m[i]) + tuple(pos[i]) + tuple(vel[i]) + (typ,)))

    def write_params(self, runid, icfile, **changes):
        params = dict(DEFAULTS)
        params['run.runid'] = runid
        params['run.icfile'] = icfile
        for key, value in changes.items():
            params[key.replace('__', '.')] = value
        sections = {}
        for key, value in params.items():
            sec, name = key.split('.')
            sections.setdefault(sec, []).append((name, value))
        fname = os.path.join(self.dir, runid + '.nml')
        with open(fname, 'w') as f:
            for sec, items in sections.items():
                f.write('&%s\n' % sec)
                for name, value in items:
                    f.write('  %s = %s\n' % (name, fortran(value)))
                f.write('/\n\n')
        return fname

    # -- running -------------------------------------------------------
    def run(self, runid, icfile, need_snaps=True, **changes):
        """Run hermite. Parameters are given as section__name=value.

        Returns the Run result: .log (screen output) and .snaps (the
        snapshot file names in time order); .last is the last snapshot.
        """
        pfile = self.write_params(runid, icfile, **changes)
        res = RunResult(self.dir, runid)
        if self.dry:
            raise Failed('dry run: input files written, nothing run')
        with open(os.path.join(self.dir, runid + '.log'), 'a') as log:
            proc = subprocess.run([self.hermite, os.path.basename(pfile)], cwd=self.dir,
                                  stdout=log, stderr=subprocess.STDOUT)
        res.read()
        if proc.returncode != 0:
            tail = '\n'.join(res.log.strip().split('\n')[-6:])
            raise Failed('hermite stopped with status %d:\n%s' % (proc.returncode, tail))
        if need_snaps and not res.snaps:
            raise Failed('the run wrote no snapshots')
        return res

    # -- checking ------------------------------------------------------
    def check(self, what, value, limit, kind='below'):
        """Record one check. kind: 'below' (|value| < limit), 'equal'."""
        if kind == 'below':
            good = abs(value) < limit
            text = '%-44s %11.3e   (limit %.1e)' % (what, value, limit)
        else:
            good = value == limit
            text = '%-44s %11s   (should be %s)' % (what, value, limit)
        self.lines.append(('ok  ' if good else 'FAIL') + '  ' + text)
        self.ok = self.ok and good
        return good


class RunResult:
    def __init__(self, directory, runid):
        self.dir = directory
        self.runid = runid
        self.log = ''
        self.snaps = []

    def read(self):
        logfile = os.path.join(self.dir, self.runid + '.log')
        if os.path.exists(logfile):
            with open(logfile) as f:
                self.log = f.read()
        self.snaps = sorted(glob.glob(os.path.join(self.dir, self.runid + '_snap_*Myr.dat')))

    @property
    def last(self):
        return read_snapshot(self.snaps[-1])

    def file(self, suffix):
        return os.path.join(self.dir, self.runid + suffix)


# ----------------------------------------------------------------------
# initial conditions used by several tests
def binary(m1, m2, a, ecc):
    """Two stars at apocentre in the x-y plane, centre of mass at rest.
    Returns m, pos, vel and the period in Myr."""
    mt = m1 + m2
    rr = np.array([a * (1 + ecc), 0.0, 0.0])
    vr = np.array([0.0, np.sqrt(G * mt / a * (1 - ecc) / (1 + ecc)), 0.0])
    pos = np.array([m2 / mt * rr, -m1 / mt * rr])
    vel = np.array([m2 / mt * vr, -m1 / mt * vr])
    period = 2 * np.pi * np.sqrt(a ** 3 / mt) * TUNIT
    return np.array([m1, m2]), pos, vel, period


def plummer(n, mstar=1.0, rhalf=1.0, seed=1):
    """An equal-mass Plummer sphere in virial equilibrium (Aarseth, Henon &
    Wielen 1974), half-mass radius rhalf pc, cut at 10 Plummer radii."""
    rng = np.random.default_rng(seed)
    a = rhalf / 1.288
    mtot = n * mstar
    pos = np.zeros((n, 3))
    vel = np.zeros((n, 3))
    i = 0
    while i < n:
        rad = 1.0 / np.sqrt(rng.random() ** (-2.0 / 3.0) - 1.0)
        if rad > 10.0:
            continue
        pos[i] = rad * sphere_point(rng)
        while True:
            q = rng.random()
            if 0.1 * rng.random() < q * q * (1 - q * q) ** 3.5:
                break
        vel[i] = q * np.sqrt(2.0) * (1 + rad * rad) ** -0.25 * sphere_point(rng)
        i += 1
    pos *= a
    vel *= np.sqrt(G * mtot / a)
    m = np.full(n, mstar)
    pos -= pos.mean(axis=0)
    vel -= vel.mean(axis=0)
    return m, pos, vel


def sphere_point(rng):
    z = 1 - 2 * rng.random()
    phi = 2 * np.pi * rng.random()
    s = np.sqrt(1 - z * z)
    return np.array([s * np.cos(phi), s * np.sin(phi), z])


def energy_error(snap):
    return snap.value('energy_error')


# ----------------------------------------------------------------------
# the tests
def test_kepler_circular(t):
    """Two stars on a circular orbit, integrated directly for 10 periods."""
    m, pos, vel, period = binary(1.0, 1.0, 1.0e-3, 0.0)
    t.write_ics('ics.dat', m, pos, vel)
    r = t.run('kepc', 'ics.dat', run__tend=10 * period, run__tout=period,
              output__tsnap=period, binaries__closebinaries=False)
    s = r.last
    t.check('energy error', energy_error(s), 1e-7)
    t.check('position after 10 periods / a', np.abs(s.pos - pos).max() / 1e-3, 1e-3)


def test_kepler_eccentric(t):
    """Two stars, eccentricity 0.9, integrated directly for 10 periods."""
    m, pos, vel, period = binary(1.0, 0.5, 1.0e-3, 0.9)
    t.write_ics('ics.dat', m, pos, vel)
    r = t.run('kepe', 'ics.dat', run__tend=10 * period, run__tout=period,
              output__tsnap=period, binaries__closebinaries=False)
    t.check('energy error', energy_error(r.last), 1e-5)


def test_binary_closed(t):
    """The same eccentric binary treated as one body: an exact Kepler orbit."""
    m, pos, vel, period = binary(1.0, 0.5, 1.0e-3, 0.9)
    t.write_ics('ics.dat', m, pos, vel)
    r = t.run('bin', 'ics.dat', run__tend=10 * period, run__tout=period,
              output__tsnap=period)
    s = r.last
    t.check('energy error', energy_error(s), 1e-10)
    t.check('position after 10 periods / a', np.abs(s.pos - pos).max() / 1e-3, 1e-6)
    t.check('stars in the snapshot', s.n, 2, 'equal')
    t.check('groups as one body', int(s.value('groups_as_one_body')), 1, 'equal')


def pythagorean():
    m = np.array([3.0, 4.0, 5.0])
    pos = np.array([[1.0, 3.0, 0.0], [-2.0, -1.0, 0.0], [1.0, -1.0, 0.0]])
    return m, pos, np.zeros((3, 3))


def test_pythagorean_direct(t):
    """Burrau's three-body problem to t = 10, integrated directly."""
    m, pos, vel = pythagorean()
    t.write_ics('ics.dat', m, pos, vel)
    r = t.run('pyth', 'ics.dat', run__tend=10 * TUNIT, run__tout=TUNIT,
              output__tsnap=TUNIT, binaries__closebinaries=False)
    t.check('energy error', energy_error(r.last), 1e-4)


def test_pythagorean_fewbody(t):
    """The same as a few-body group, with nine light, distant onlookers
    (a group of three is only taken out if there are 12 bodies or more)."""
    m, pos, vel = pythagorean()
    rng = np.random.default_rng(7)
    far = np.array([200.0 * sphere_point(rng) for _ in range(9)])
    m = np.concatenate([m, np.full(9, 1.0e-3)])
    pos = np.vstack([pos, far])
    vel = np.zeros((12, 3))
    t.write_ics('ics.dat', m, pos, vel)
    r = t.run('pythf', 'ics.dat', run__tend=10 * TUNIT, run__tout=TUNIT,
              output__tsnap=TUNIT)
    s = r.last
    t.check('energy error', energy_error(s), 1e-5)
    t.check('a few-body group was formed', int(s.value('few_body_groups_formed')) >= 1,
            True, 'equal')
    t.check('stars in the snapshot', s.n, 12, 'equal')


def test_fewbody_small_n(t):
    """With only three bodies an unstable triple is NOT taken out as a
    group (it would not pay), but it is still integrated correctly."""
    m, pos, vel = pythagorean()
    t.write_ics('ics.dat', m, pos, vel)
    r = t.run('pyth3', 'ics.dat', run__tend=10 * TUNIT, run__tout=TUNIT,
              output__tsnap=TUNIT)
    s = r.last
    t.check('few-body groups formed', int(s.value('few_body_groups_formed')), 0, 'equal')
    t.check('energy error', energy_error(s), 1e-4)


def test_triple_frozen(t):
    """A stable hierarchical triple: frozen as two Kepler orbits (even
    with only three bodies, because a frozen triple always pays)."""
    m1, p1, v1, pin = binary(1.0, 1.0, 1.0e-4, 0.0)
    # outer star on a circular orbit of radius 5e-3 pc about the binary
    mt = 3.0
    rout = 5.0e-3
    vout = np.sqrt(G * mt / rout)
    pos = np.vstack([p1 - np.array([rout / 3, 0, 0]), [[2 * rout / 3, 0, 0]]])
    vel = np.vstack([v1 - np.array([0, vout / 3, 0]), [[0, 2 * vout / 3, 0]]])
    pout = 2 * np.pi * np.sqrt(rout ** 3 / mt) * TUNIT
    t.write_ics('ics.dat', [1.0, 1.0, 1.0], pos, vel)
    r = t.run('trip', 'ics.dat', run__tend=4 * pout, run__tout=pout / 4,
              output__tsnap=pout)
    s = r.last
    t.check('triples frozen', int(s.value('triples_frozen')), 1, 'equal')
    t.check('energy error', energy_error(s), 1e-3)
    sep = np.linalg.norm(s.pos[0] - s.pos[1])
    t.check('inner separation / a - 1', sep / 1.0e-4 - 1, 1e-3)


def test_plummer_direct(t):
    """256-star Plummer sphere for 2 Myr by direct summation."""
    m, pos, vel = plummer(256)
    t.write_ics('ics.dat', m, pos, vel)
    r = t.run('plum', 'ics.dat', run__tend=2.0)
    s = r.last
    t.check('energy error', energy_error(s), 1e-5)
    t.check('virial ratio - 0.5', s.value('virial_ratio') - 0.5, 0.1)
    t.check('stars in the snapshot', s.n, 256, 'equal')


def test_plummer_neighbours(t):
    """The same with the Ahmad-Cohen neighbour scheme."""
    m, pos, vel = plummer(256)
    t.write_ics('ics.dat', m, pos, vel)
    r = t.run('plac', 'ics.dat', run__tend=2.0, neighbours__useac=True)
    s = r.last
    t.check('energy error', energy_error(s), 1e-4)
    t.check('neighbour_scheme in physics', 'neighbour_scheme' in s.header['physics'],
            True, 'equal')


def test_plummer_tree(t):
    """The same with the Barnes-Hut tree."""
    m, pos, vel = plummer(256)
    t.write_ics('ics.dat', m, pos, vel)
    r = t.run('pltr', 'ics.dat', run__tend=2.0, tree__usetree=True)
    t.check('energy error', energy_error(r.last), 5e-3)


def test_plummer_tide(t):
    """The same in the solar-neighbourhood tidal field: Jacobi energy."""
    m, pos, vel = plummer(256)
    t.write_ics('ics.dat', m, pos, vel)
    r = t.run('pltd', 'ics.dat', run__tend=2.0, tides__tidal='solar')
    s = r.last
    t.check('Jacobi energy error', energy_error(s), 1e-5)
    t.check('tidal keys in the header', 'tidal_omega_km_s_kpc' in s.header, True, 'equal')


def test_header_flexible(t):
    """Header keys come and go with the physics; unknown keys and extra
    columns in an initial conditions file are ignored."""
    m, pos, vel = plummer(64)
    t.write_ics('ics.dat', m, pos, vel)
    r = t.run('hdr', 'ics.dat', run__tend=0.25)
    s = r.last
    t.check('no tidal keys without a tide', 'tidal_omega_km_s_kpc' in s.header,
            False, 'equal')
    t.check("'physics' key present", 'physics' in s.header, True, 'equal')
    # a snapshot (11 columns, many keys) with a made-up key and column added
    lines = open(r.snaps[-1]).read().split('\n')
    out = []
    for line in lines:
        if line.startswith('# columns'):
            line += ' luminosity'
        elif line.startswith('# end_header'):
            out.append('# stellar_evolution             = made_up_for_the_test')
        elif line and not line.startswith('#'):
            line += '  1.0'
        out.append(line)
    with open(os.path.join(t.dir, 'ics2.dat'), 'w') as f:
        f.write('\n'.join(out))
    r2 = t.run('hdr2', 'ics2.dat', run__tend=0.5)
    s2 = r2.last
    t.check('run from the changed snapshot: stars', s2.n, 64, 'equal')
    t.check('run from the changed snapshot: time', s2.time - 0.5, 1e-9)


def test_merger(t):
    """Two stars collide head on and merge; a third star looks on."""
    m = np.array([1.0, 2.0, 0.5])
    pos = np.array([[2.0e-4, 0, 0], [-1.0e-4, 0, 0], [0, 0.05, 0]])
    vel = np.array([[-2.0, 0, 0], [1.0, 0, 0], [0, 0, 0.0]])
    t.write_ics('ics.dat', m, pos, vel)
    r = t.run('mrg', 'ics.dat', run__tend=0.01, run__tout=0.0025,
              output__tsnap=0.0025, collisions__rcoll=10.0)
    s = r.last
    t.check('stars in the snapshot', s.n, 2, 'equal')
    t.check('total mass - 3.5', s.m.sum() - 3.5, 1e-12)
    t.check('mergers in the header', int(s.value('mergers')), 1, 'equal')
    t.check('energy error (with merger energy)', energy_error(s), 1e-5)
    mom = (s.m[:, None] * s.vel).sum(axis=0)
    t.check('total momentum (Msun km/s)', np.abs(mom).max(), 1e-9)
    mfile = r.file('_mergers.dat')
    nlines = 0
    if os.path.exists(mfile):
        nlines = sum(1 for line in open(mfile) if line.strip() and not line.startswith('#'))
    t.check('lines in the merger file', nlines, 1, 'equal')
    t.check("'** MERGER' on the screen", '** MERGER' in r.log, True, 'equal')
    t.check('merged star keeps id 2', 2 in list(s.id), True, 'equal')


def test_merger_off(t):
    """With mergers off a collision stops the run with an alert."""
    m = np.array([1.0, 2.0, 0.5])
    pos = np.array([[2.0e-4, 0, 0], [-1.0e-4, 0, 0], [0, 0.05, 0]])
    vel = np.array([[-2.0, 0, 0], [1.0, 0, 0], [0, 0, 0.0]])
    t.write_ics('ics.dat', m, pos, vel)
    try:
        t.run('mrgoff', 'ics.dat', run__tend=0.01, run__tout=0.0025,
              output__tsnap=0.0025, collisions__rcoll=10.0, collisions__mergers=False)
        stopped = False
    except Failed as err:
        stopped = 'STELLAR COLLISION' in str(err) or 'stopped' in str(err)
    t.check('run stopped at the collision', stopped, True, 'equal')


def test_restart(t):
    """A run stopped half way and restarted ends where the unbroken run does."""
    m, pos, vel = plummer(128, seed=5)
    t.write_ics('ics.dat', m, pos, vel)
    ref = t.run('ref', 'ics.dat', run__tend=1.0, output__trestart=0.5).last
    t.run('rst', 'ics.dat', run__tend=0.5, output__trestart=0.5)
    r = t.run('rst', 'ics.dat', run__tend=1.0, output__trestart=0.5, run__restart=True)
    s = r.last
    t.check('time reached - 1 Myr', s.time - 1.0, 1e-9)
    t.check('largest position difference (pc)', np.abs(s.pos - ref.pos).max(), 1e-9)
    t.check('energy error difference', energy_error(s) - energy_error(ref), 1e-9)


def test_restart_neighbours(t):
    """The same with the neighbour scheme (lists are in the restart file)."""
    m, pos, vel = plummer(192, seed=5)
    t.write_ics('ics.dat', m, pos, vel)
    ac = dict(neighbours__useac=True, output__trestart=0.5)
    ref = t.run('ref', 'ics.dat', run__tend=1.0, **ac).last
    t.run('rst', 'ics.dat', run__tend=0.5, **ac)
    s = t.run('rst', 'ics.dat', run__tend=1.0, run__restart=True, **ac).last
    t.check('largest position difference (pc)', np.abs(s.pos - ref.pos).max(), 1e-9)


def test_ensemble(t):
    """An ensemble of three runs, each with its own output files."""
    for k in (1, 2, 3):
        m, pos, vel = plummer(32, seed=10 + k)
        t.write_ics('ics_e%03d.dat' % k, m, pos, vel)
    r = t.run('ens', 'ics.dat', need_snaps=False, run__tend=0.5, run__nensemble=3,
              run__nparallel=3)
    t.check("'3 of 3 members ran to the end'", '3 of 3 members ran to the end' in r.log,
            True, 'equal')
    ends = []
    for k in (1, 2, 3):
        member = RunResult(t.dir, 'ens_e%03d' % k)
        member.read()
        t.check('member %d: snapshots' % k, len(member.snaps), 3, 'equal')
        if member.snaps:
            s = member.last
            t.check('member %d: ensemble_member in header' % k,
                    int(s.value('ensemble_member')), k, 'equal')
            t.check('member %d: energy error' % k, energy_error(s), 1e-4)
            ends.append(s.pos[0, 0])
        t.check('member %d: its own log file' % k,
                os.path.exists(member.file('.log')), True, 'equal')
    t.check('members differ from each other', len(set(ends)), 3, 'equal')


# ----------------------------------------------------------------------
def all_tests():
    tests = [(name[5:], fn) for name, fn in globals().items()
             if name.startswith('test_') and callable(fn)]
    return tests


def main():
    ap = argparse.ArgumentParser(description='HENDOS regression tests')
    ap.add_argument('--hermite', default=os.path.join(HERE, '..', 'run', 'hermite'),
                    help='the executable to test (default ../run/hermite)')
    ap.add_argument('--only', default='', help='only tests whose name contains this')
    ap.add_argument('--list', action='store_true', help='list the tests and stop')
    ap.add_argument('--keep', action='store_true', help='keep the working files')
    ap.add_argument('--dry', action='store_true',
                    help='write the input files of each test but run nothing')
    args = ap.parse_args()

    tests = [(n, f) for n, f in all_tests() if args.only in n]
    if args.list:
        for name, fn in tests:
            print('%-22s %s' % (name, ' '.join(fn.__doc__.split())))
        return 0
    hermite = os.path.abspath(args.hermite)
    if not args.dry and not os.path.exists(hermite):
        print('cannot find the executable %s - run make in HENDOS/run first' % hermite)
        return 2
    workdir = os.path.join(HERE, 'work')
    os.makedirs(workdir, exist_ok=True)

    npass = 0
    failed = []
    for name, fn in tests:
        t = Tester(name, hermite, workdir, args.dry)
        t0 = time.time()
        try:
            fn(t)
        except Failed as err:
            t.ok = False
            t.lines.append('FAIL  ' + str(err))
        except Exception as err:                       # a broken test or output file
            t.ok = False
            t.lines.append('FAIL  %s: %s' % (type(err).__name__, err))
        print('%-22s %s   %.1f s' % (name, 'PASS' if t.ok else 'FAIL', time.time() - t0))
        for line in t.lines:
            print('    ' + line)
        if t.ok:
            npass += 1
        else:
            failed.append(name)
    print('-' * 70)
    print('%d of %d tests passed' % (npass, len(tests)))
    if failed:
        print('failed: ' + ', '.join(failed))
        print('working files are in %s' % workdir)
    elif not args.keep and not args.dry:
        shutil.rmtree(workdir)
    return 0 if not failed else 1


if __name__ == '__main__':
    sys.exit(main())
