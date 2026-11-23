"""
HENDOS - HErmite N-body Dynamics Of Systems
Benchmark HENDOS against Starlab's kira on the same initial conditions.

    python bench_kira.py <HENDOS ic file> --tend 5 [options]

HENDOS can do the same problem in several ways, so it is run in one or more
MODES (--modes), each of which is compared with kira separately:
    direct       direct summation, binaries and few-body groups on [standard]
    neighbours   the same with the Ahmad-Cohen neighbour scheme
    tree         Barnes-Hut tree (approximate forces)
    nofew        direct, binaries on but few-body groups off
    nobin        direct, every star integrated directly (no binaries/groups)
The modes are the table HENDOS_MODES below - add your own there.

For kira and for each HENDOS mode it runs the same stars to the same time
at several settings of the accuracy parameter, and for every run records
    wall-clock time
    energy error  |E_end - E_start| / |E_start|
Both are worked out here, in the same way for both codes: the time by
timing the run, the energy by direct summation over the stars of the first
and last snapshots. Neither code's own report is used for the comparison.

The fair comparison is at EQUAL ACCURACY, not equal parameters (eta means
different things in the two codes), so the result is a curve for kira and
for each HENDOS mode: time taken against energy error reached. The lower
curve is faster at that accuracy.

It also plots the half-mass radius against time for the most accurate run
of each, as a check that they are doing the same physics, and for every
HENDOS run it records what the binary and few-body machinery did (groups
closed, few-body groups formed, triples frozen), HENDOS's own energy error
and whether the number of stars changed.

Mergers are switched off in HENDOS for the comparison (rcoll = 0): kira
without stellar radii does not merge stars, and a merger removes energy,
which would spoil the energy comparison. Restart files are off too.

Everything goes in a new directory (--out, default bench_<ic file name>):
    ic.dyn, ic.dyn.units.json     the stars in Starlab format and the units
    kira_eta<..>.dyn / .log       kira output and log for each run
    hendos_<mode>_eta<..>/        HENDOS parameter file, log and snapshots
    results.txt, results.json     the table
    bench_time_error.png          wall-clock time against energy error
    bench_halfmass.png            half-mass radius against time

BEFORE THE FIRST RUN check the two command lines just below against your
installation (kira --help). They are the only places that know how the
codes are run.

Needs numpy and matplotlib, and the HENDOS analysis directory next to this
one.
"""
import argparse
import json
import os
import shlex
import subprocess
import sys
import time

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, '..', 'analysis'))
sys.path.insert(0, HERE)
import hendosanalysis as ha                                   # noqa: E402
from hendosio import Run                                      # noqa: E402
from starlabio import (G_PHYS, Units, hendos_to_dyn, kinetic_energy,   # noqa: E402
                       potential_energy, read_dyn)

# ----------------------------------------------------------------------
# HOW THE CODES ARE RUN - check these against your installation
# ----------------------------------------------------------------------
# kira reads a dyn file on standard input, writes snapshots to standard
# output and its log to standard error. All in N-body units. The options,
# as kira_init.C of Starlab 4.4 reads them:
#   -t  length of the run (delta_t)       -d  interval between log outputs
#   -D  interval between snapshots (a number; 'all', 'x<n>' etc. mean
#       other things)                     -a  accuracy parameter eta
#   -e  softening (0, as HENDOS)          -u  switch unperturbed binary
#                                             motion on (HENDOS does the same)
# Run  python bench_kira.py <ic file> --check  first: it runs kira for one
# output interval on the stars and checks that what comes back makes sense.
KIRA_CMD = '{kira} -t {tend} -d {dt} -D {dt} -a {eta} -e 0 -u'

# HENDOS is run in its own directory with the parameter file written below.
HERMITE_CMD = '{hermite} params.nml'

HENDOS_PARAMS = """&run
  runid  = 'bench'
  icfile = '{icfile}'
  tend   = {tend}
  tout   = {dt}
  eps    = 0.0
  restart = .false.
/
&output
  datadir = 'data'
  screen  = 'verbose'
  tsnap   = {dt}
  trestart = 0.0
/
&analysis
  tana = 0.0
/
&integrator
  eta  = {eta}
  etap = {etap}
/
&binaries
  closebinaries = {closebinaries}
  fewbody = {fewbody}
/
&neighbours
  useac = {useac}
  etar  = {etar}
/
&tree
  usetree = {usetree}
/
&collisions
  rcoll = 0.0
/
"""
# rcoll = 0: kira (without stellar radii) does not merge or stop for
# collisions, so HENDOS must not either - and a merger removes energy.

# THE WAYS HENDOS IS RUN. Each mode is one curve on the plots. To add one,
# add a line here (and a parameter to HENDOS_PARAMS above if it needs one).
#   label           name on the plots and in the table
#   colour, marker  how it is drawn (colours checked for colour-blind readers;
#                   the markers and the labels on the curves back them up)
#   the rest        values for the {..} places in HENDOS_PARAMS
# HENDOS has two timestep parameters as well as eta: etap, the pairwise
# free-fall/flyby limit, and (neighbour scheme) etar. Both are tied to eta in
# every run so that all three tighten together - with etap held fixed the
# close encounters set a floor to the accuracy and smaller eta buys nothing
# (found in the first benchmarks: with etap = 0.2 cutting eta from 0.04 to
# 0.01 only reduced the error 3 times; with etap = 5 eta it is 100 times).
HENDOS_MODES = {
    'direct':     dict(label='HENDOS direct', colour='#2a78d6', marker='o',
                       closebinaries=True, fewbody=True, useac=False, usetree=False),
    'neighbours': dict(label='HENDOS neighbours', colour='#1baf7a', marker='^',
                       closebinaries=True, fewbody=True, useac=True, usetree=False),
    'tree':       dict(label='HENDOS tree', colour='#eda100', marker='D',
                       closebinaries=True, fewbody=True, useac=False, usetree=True),
    'nofew':      dict(label='HENDOS no few-body', colour='#e87ba4', marker='v',
                       closebinaries=True, fewbody=False, useac=False, usetree=False),
    'nobin':      dict(label='HENDOS no binaries', colour='#6250d6', marker='P',
                       closebinaries=False, fewbody=False, useac=False, usetree=False),
}
ETAR_OVER_ETA = 0.5
ETAP_OVER_ETA = 5.0
KIRA_STYLE = dict(label='kira', colour='#eb6834', marker='s')

# header keys of the last HENDOS snapshot kept in the results (if present)
HENDOS_KEYS = ('groups_as_one_body', 'binaries_closed', 'binaries_opened',
               'few_body_groups_formed', 'few_body_groups_opened', 'triples_frozen',
               'mergers')

TEXT, TEXT2, GRID, SURFACE = '#0b0b0b', '#52514e', '#dddcd6', '#fcfcfb'


def fbool(value):
    return '.true.' if value else '.false.'


def run_timed(command, cwd=None, stdin=None, stdout=None, stderr=None, env=None):
    """Run a command and return (wall-clock seconds, return code)."""
    t0 = time.perf_counter()
    result = subprocess.run(shlex.split(command), cwd=cwd, stdin=stdin, stdout=stdout,
                            stderr=stderr, env=env)
    return time.perf_counter() - t0, result.returncode


def total_energy(m, pos, vel, g):
    return kinetic_energy(m, vel) + potential_energy(m, pos, g)


def half_mass_curve(snaps):
    return (np.array([s.time for s in snaps]),
            np.array([ha.half_mass_radius(s) for s in snaps]))


# ----------------------------------------------------------------------
def run_kira(out, units, tend, dt, eta, kira='kira'):
    """One kira run. Returns a dictionary of results, and the snapshots in physical units."""
    tag = 'kira_eta%g' % eta
    command = KIRA_CMD.format(kira=kira, tend=repr(tend / units.time), dt=repr(dt / units.time),
                              eta=eta)
    dyn_out = os.path.join(out, tag + '.dyn')
    with open(os.path.join(out, 'ic.dyn')) as fin, open(dyn_out, 'w') as fout, \
            open(os.path.join(out, tag + '.log'), 'w') as flog:
        wall, code = run_timed(command, stdin=fin, stdout=fout, stderr=flog)
    result = dict(code='kira', series='kira', eta=eta, wall_s=wall, ok=(code == 0),
                  command=command)
    snaps = [s.to_physical(units) for s in read_dyn(dyn_out)] if code == 0 else []
    return result, snaps


def run_hendos(out, icfile, hermite, tend, dt, eta, threads, mode):
    """One HENDOS run in the given mode. Returns results and the snapshots."""
    opts = HENDOS_MODES[mode]
    tag = 'hendos_%s_eta%g' % (mode, eta)
    rundir = os.path.join(out, tag)
    datadir = os.path.join(rundir, 'data')
    os.makedirs(datadir, exist_ok=True)
    for name in os.listdir(datadir):
        os.remove(os.path.join(datadir, name))            # HENDOS will not overwrite a run
    with open(os.path.join(rundir, 'params.nml'), 'w') as f:
        f.write(HENDOS_PARAMS.format(
            icfile=os.path.abspath(icfile), tend=tend, dt=dt, eta=eta,
            etar=ETAR_OVER_ETA * eta, etap=ETAP_OVER_ETA * eta, closebinaries=fbool(opts['closebinaries']),
            fewbody=fbool(opts['fewbody']), useac=fbool(opts['useac']),
            usetree=fbool(opts['usetree'])))
    env = dict(os.environ, OMP_NUM_THREADS=str(threads))
    command = HERMITE_CMD.format(hermite=os.path.abspath(hermite))
    with open(os.path.join(rundir, 'hermite.log'), 'w') as flog:
        wall, code = run_timed(command, cwd=rundir, stdout=flog, stderr=subprocess.STDOUT,
                               env=env)
    result = dict(code='HENDOS', series=mode, eta=eta, etap=ETAP_OVER_ETA * eta, wall_s=wall,
                  ok=(code == 0), command=command, threads=threads)
    snaps = []
    if code == 0:
        run = Run(datadir, 'bench')
        snaps = [run.snapshot(i) for i in range(run.nsnap)]
        if snaps:
            last = snaps[-1]
            result['own_energy_error'] = last.value('energy_error')
            result['physics'] = last.header.get('physics', '')
            for key in HENDOS_KEYS:
                if key in last.header:
                    result[key] = int(float(last.header[key]))
    return result, snaps


def finish(result, snaps, e0, tend, n0):
    """Add the energy errors and end time to a result, from the snapshots.

    energy_error_end  |E - E0| / |E0| at the last snapshot
    energy_error      the largest |E - E0| / |E0| at any snapshot - the one
                      used in the plots and table. The error at the end alone
                      is noisy: it can change sign during a run, so a single
                      run can by chance end with a tiny error (both codes did
                      in the first benchmarks).
    """
    if not snaps:
        result['ok'] = False
        return
    errs = [abs(total_energy(s.m, s.pos, s.vel, G_PHYS) - e0) / abs(e0) for s in snaps]
    last = snaps[-1]
    result['energy_error_end'] = errs[-1]
    result['energy_error'] = max(errs)
    result['t_end_Myr'] = float(last.time)
    result['n_end'] = int(last.n)
    result['n_changed'] = bool(last.n != n0)
    result['reached_end'] = bool(abs(last.time - tend) < 1e-3 * tend)


def series_style(series):
    return KIRA_STYLE if series == 'kira' else HENDOS_MODES[series]


# ----------------------------------------------------------------------
def make_plots(out, results, curves, title):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt

    def style(ax):
        ax.set_facecolor(SURFACE)
        ax.grid(True, color=GRID, lw=0.6, which='both')
        ax.set_axisbelow(True)
        for side in ('top', 'right'):
            ax.spines[side].set_visible(False)
        for side in ('left', 'bottom'):
            ax.spines[side].set_color(GRID)
        ax.tick_params(colors=TEXT2, which='both')

    order = []
    for r in results:
        if r['series'] not in order:
            order.append(r['series'])

    # wall-clock time against energy error
    fig, ax = plt.subplots(figsize=(8.0, 5.8), facecolor=SURFACE)
    for series in order:
        st = series_style(series)
        rows = sorted((r for r in results if r['series'] == series and r['ok']
                       and r.get('energy_error', 0) > 0), key=lambda r: r['energy_error'])
        if not rows:
            continue
        x = [r['energy_error'] for r in rows]
        y = [r['wall_s'] for r in rows]
        ax.plot(x, y, color=st['colour'], lw=2, marker=st['marker'], ms=8, mec=SURFACE,
                mew=1.5, label=st['label'])
        ax.annotate('  ' + st['label'], (x[-1], y[-1]), color=TEXT2, va='center',
                    fontsize=9)
        for r in rows:
            ax.annotate('%g' % r['eta'], (r['energy_error'], r['wall_s']),
                        textcoords='offset points', xytext=(0, 9), ha='center',
                        fontsize=8, color=TEXT2)
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.margins(x=0.22)
    ax.set_xlabel('largest relative energy error during the run   (numbers on the points: eta)',
                  color=TEXT2)
    ax.set_ylabel('wall-clock time (s)', color=TEXT2)
    ax.set_title(title, color=TEXT, loc='left', fontsize=11)
    ax.legend(frameon=False, labelcolor=TEXT)
    style(ax)
    fig.tight_layout()
    fig.savefig(os.path.join(out, 'bench_time_error.png'), dpi=150, facecolor=SURFACE)

    # half-mass radius against time, every run: the most accurate run of each
    # code or mode thick, the others thin. After a cold collapse, or once
    # binaries form, runs diverge chaotically - so compare the spread of the
    # curves of one code with that of the other, not one curve with one curve.
    fig, ax = plt.subplots(figsize=(8.0, 5.0), facecolor=SURFACE)
    for series in order:
        runs = sorted((key for key in curves if key[0] == series), key=lambda key: key[1])
        st = series_style(series)
        for i, key in enumerate(runs):
            t, rh = curves[key]
            best = (i == 0)
            ax.plot(t, rh, color=st['colour'], lw=2 if best else 0.9,
                    alpha=1.0 if best else 0.55,
                    label=('%s (eta=%s)' % (st['label'], ', '.join('%g' % k[1] for k in runs))
                           if best else None))
            if best:
                ax.annotate('  ' + st['label'], (t[-1], rh[-1]), color=TEXT2, va='center',
                            fontsize=9)
    ax.margins(x=0.18)
    ax.set_xlabel('time (Myr)', color=TEXT2)
    ax.set_ylabel('half-mass radius (pc)', color=TEXT2)
    ax.set_title(title, color=TEXT, loc='left', fontsize=11)
    ax.legend(frameon=False, labelcolor=TEXT)
    style(ax)
    fig.tight_layout()
    fig.savefig(os.path.join(out, 'bench_halfmass.png'), dpi=150, facecolor=SURFACE)


def check(a, out, units, ic, e0, dt):
    """Run kira and HENDOS (direct) for one output interval on the stars and
    check that both read them, reach the end, keep every star and conserve
    energy, and that the two agree on where the stars are. Returns 0 if all
    is well."""
    import shutil
    good = True

    def report(ok, text):
        print('   %s  %s' % ('ok  ' if ok else 'FAIL', text))
        return ok

    print('\nCHECK: one output interval (%.4g Myr) with each code' % dt)
    if shutil.which(a.kira) is None and not os.path.exists(a.kira):
        report(False, 'cannot find kira (%s) - put it on your path or give --kira' % a.kira)
        good = False
        ksnaps = []
    else:
        kres, ksnaps = run_kira(out, units, dt, dt / 2, a.eta_kira[-1], a.kira)
        print('   kira command:  %s' % kres['command'])
        good &= report(kres['ok'], 'kira ran (log: %s)' % os.path.join(out, 'kira_eta%g.log'
                                                                     % a.eta_kira[-1]))
        good &= report(len(ksnaps) > 0, 'kira wrote snapshots that could be read (%d)'
                       % len(ksnaps))
        if ksnaps:
            last = ksnaps[-1]
            good &= report(last.n == ic.n, 'kira kept every star (%d of %d)' % (last.n, ic.n))
            good &= report(abs(last.time - dt) < 1e-3 * dt,
                           'kira reached the end (%.6g of %.6g Myr) - so its time units '
                           'are being converted correctly' % (last.time, dt))
            e1 = total_energy(last.m, last.pos, last.vel, G_PHYS)
            good &= report(abs(e1 - e0) < 1e-3 * abs(e0),
                           'kira energy error %.2e (should be small)' % (abs(e1 - e0) / abs(e0)))
    hres, hsnaps = run_hendos(out, a.icfile, a.hermite, dt, dt / 2, a.eta_hendos[-1],
                              a.threads, 'direct')
    good &= report(hres['ok'] and len(hsnaps) > 0, 'HENDOS ran (log: %s)'
                   % os.path.join(out, 'hendos_direct_eta%g' % a.eta_hendos[-1], 'hermite.log'))
    if hsnaps:
        last = hsnaps[-1]
        e1 = total_energy(last.m, last.pos, last.vel, G_PHYS)
        good &= report(abs(e1 - e0) < 1e-3 * abs(e0),
                       'HENDOS energy error %.2e (should be small)' % (abs(e1 - e0) / abs(e0)))
        if ksnaps:
            # the same stars in the same order: compare positions, relative
            # to the size of the system
            k, h = ksnaps[-1], last
            if k.n == h.n:
                size = np.sqrt((ic.pos ** 2).sum(1).mean())
                ik = np.argsort(k.id) if hasattr(k, 'id') else np.arange(k.n)
                ih = np.argsort(h.id)
                kp = k.pos[ik] - (k.m[:, None] * k.pos).sum(0) / k.m.sum()
                hp = h.pos[ih] - (h.m[:, None] * h.pos).sum(0) / h.m.sum()
                diff = np.median(np.sqrt(((kp - hp) ** 2).sum(1))) / size
                good &= report(diff < 1e-2, 'kira and HENDOS agree on the positions after '
                               'one interval (median difference %.1e of the rms radius)' % diff)
    print('\nAll is well: run the benchmark.' if good else
          '\nSomething is wrong: look at the logs above before running the benchmark.')
    return 0 if good else 1


def main():
    p = argparse.ArgumentParser(description='Benchmark HENDOS against Starlab kira.')
    p.add_argument('icfile', help='HENDOS initial conditions file (from makeics)')
    p.add_argument('--tend', type=float, required=True, help='end time (Myr)')
    p.add_argument('--nout', type=int, default=20,
                   help='number of outputs and snapshots over the run [20]')
    p.add_argument('--modes', nargs='+', default=['direct'], choices=sorted(HENDOS_MODES),
                   help='the ways HENDOS is run, one curve each [direct]')
    p.add_argument('--eta-hendos', type=float, nargs='+', default=[0.04, 0.02, 0.01],
                   help='HENDOS eta values to run [0.04 0.02 0.01]')
    p.add_argument('--eta-kira', type=float, nargs='+', default=[0.2, 0.1, 0.05],
                   help='kira -a values to run [0.2 0.1 0.05]')
    p.add_argument('--hermite', default=os.path.join(HERE, '..', 'run', 'hermite'),
                   help='the HENDOS executable [../run/hermite]')
    p.add_argument('--threads', type=int, default=1,
                   help='OpenMP threads for HENDOS [1, as kira is serial]')
    p.add_argument('--out', default=None, help='directory for everything [bench_<ic name>]')
    p.add_argument('--skip-kira', action='store_true', help='do not run kira')
    p.add_argument('--kira', default='kira', help='the kira executable [kira, on your path]')
    p.add_argument('--check', action='store_true',
                   help='only check that kira and HENDOS run on these stars and agree '
                        'for one output interval, then stop')
    p.add_argument('--skip-hendos', action='store_true', help='do not run HENDOS')
    a = p.parse_args()

    out = a.out or 'bench_' + os.path.splitext(os.path.basename(a.icfile))[0]
    os.makedirs(out, exist_ok=True)
    dt = a.tend / a.nout

    units, ic = hendos_to_dyn(a.icfile, os.path.join(out, 'ic.dyn'))
    units.save(os.path.join(out, 'ic.dyn.units.json'))
    com_v = (ic.m[:, None] * ic.vel).sum(0) / ic.m.sum()
    e0 = total_energy(ic.m, ic.pos, ic.vel - com_v, G_PHYS)
    print('%d stars, %.1f Msun, virial radius %.4f pc' % (ic.n, units.mass, units.length))
    print('one N-body time unit = %.5f Myr, so tend = %.3f N-body units'
          % (units.time, a.tend / units.time))

    if a.check:
        return check(a, out, units, ic, e0, dt)

    results, curves = [], {}
    jobs = []
    if not a.skip_hendos:
        jobs += [(mode, eta) for mode in a.modes for eta in a.eta_hendos]
    if not a.skip_kira:
        jobs += [('kira', eta) for eta in a.eta_kira]
    for series, eta in jobs:
        label = series_style(series)['label']
        print('running %s with eta = %g ...' % (label, eta), flush=True)
        if series == 'kira':
            result, snaps = run_kira(out, units, a.tend, dt, eta, a.kira)
        else:
            result, snaps = run_hendos(out, a.icfile, a.hermite, a.tend, dt, eta, a.threads,
                                       series)
        finish(result, snaps, e0, a.tend, int(ic.n))
        results.append(result)
        if result['ok']:
            print('   %.2f s, largest energy error %.2e' % (result['wall_s'],
                                                          result['energy_error']))
            curves[(series, eta)] = half_mass_curve(snaps)
        else:
            print('   FAILED - see the log in %s' % out)

    title = '%s   N = %d   to %g Myr' % (os.path.basename(a.icfile), ic.n, a.tend)
    with open(os.path.join(out, 'results.json'), 'w') as f:
        json.dump(dict(icfile=a.icfile, n=int(ic.n), tend_Myr=a.tend, results=results), f, indent=1)
    lines = [title, '',
             '%-19s %7s %10s %13s %13s %10s  %s' % ('run', 'eta', 'wall (s)', 'max |dE/E|',
                                                    'end |dE/E|', 'end (Myr)', 'notes')]
    for r in results:
        label = series_style(r['series'])['label']
        if r['ok']:
            notes = []
            if not r.get('reached_end'):
                notes.append('DID NOT REACH tend')
            if r.get('n_changed'):
                notes.append('NUMBER OF STARS CHANGED to %d' % r['n_end'])
            if r['code'] == 'HENDOS':
                notes.append('%d thread(s)' % r['threads'])
                if r.get('etap'):
                    notes.append('etap %g' % r['etap'])
                if 'own_energy_error' in r:
                    notes.append('own dE/E %.1e' % r['own_energy_error'])
                if 'binaries_closed' in r:
                    notes.append('binaries closed %d' % r['binaries_closed'])
                if 'few_body_groups_formed' in r:
                    notes.append('few-body %d, frozen %d' % (r['few_body_groups_formed'],
                                                            r.get('triples_frozen', 0)))
            lines.append('%-19s %7g %10.2f %13.3e %13.3e %10.4f  %s'
                         % (label, r['eta'], r['wall_s'], r['energy_error'],
                            r['energy_error_end'], r['t_end_Myr'], ', '.join(notes)))
        else:
            lines.append('%-19s %7g %10s %13s %13s %10s  FAILED'
                         % (label, r['eta'], '-', '-', '-', '-'))
    text = '\n'.join(lines)
    with open(os.path.join(out, 'results.txt'), 'w') as f:
        f.write(text + '\n')
    print('\n' + text)
    if any(r['ok'] for r in results):
        make_plots(out, results, curves, title)
        print('\nplots and tables are in %s' % out)


if __name__ == '__main__':
    sys.exit(main())
