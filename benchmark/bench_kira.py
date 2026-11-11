"""
HENDOS - HErmite N-body Dynamics Of Systems
Benchmark HENDOS against Starlab's kira on the same initial conditions.

    python bench_kira.py <HENDOS ic file> --tend 5 [options]

For each code it runs the same stars to the same time at several settings
of the code's accuracy parameter, and for every run records
    wall-clock time
    energy error  |E_end - E_start| / |E_start|
Both are worked out here, in the same way for both codes: the time by
timing the run, the energy by direct summation over the stars of the first
and last snapshots. Neither code's own report is used for the comparison.

The fair comparison is at EQUAL ACCURACY, not equal parameters (eta means
different things in the two codes), so the result is a curve for each
code: time taken against energy error reached. The code whose curve is
lower is faster at that accuracy.

It also plots the half-mass radius against time for the most accurate run
of each code, as a check that they are doing the same physics.

Everything goes in a new directory (--out, default bench_<ic file name>):
    ic.dyn, ic.dyn.units.json     the stars in Starlab format and the units
    kira_eta<..>.dyn / .log       kira output and log for each run
    hendos_eta<..>/               HENDOS parameter file, log and snapshots
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
# output and its log to standard error. All in N-body units.
#   -t  end time            -d  interval between log outputs
#   -D  interval between snapshots        -a  accuracy parameter
# Add anything else your kira needs here (for example the option that
# switches on unperturbed binary motion, or a softening).
KIRA_CMD = 'kira -t {tend} -d {dt} -D {dt} -a {eta}'

# HENDOS is run in its own directory with the parameter file written below.
HERMITE_CMD = '{hermite} params.nml'

HENDOS_PARAMS = """&run
  runid  = 'bench'
  icfile = '{icfile}'
  tend   = {tend}
  tout   = {dt}
  eps    = 0.0
/
&output
  datadir = 'data'
  screen  = 'standard'
  tsnap   = {dt}
/
&analysis
  tana = 0.0
/
&integrator
  eta = {eta}
/
&collisions
  rcoll = 0.0
/
"""
# rcoll = 0: kira (without stellar evolution) does not stop for collisions,
# so HENDOS must not either.

# colours: two codes, two hues that can be told apart by colour-blind readers
BLUE, ORANGE = '#2a78d6', '#eb6834'
TEXT, TEXT2, GRID, SURFACE = '#0b0b0b', '#52514e', '#dddcd6', '#fcfcfb'


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
def run_kira(out, units, tend, dt, eta):
    """One kira run. Returns a dictionary of results, and the snapshots in physical units."""
    tag = 'kira_eta%g' % eta
    command = KIRA_CMD.format(tend=repr(tend / units.time), dt=repr(dt / units.time), eta=eta)
    dyn_out = os.path.join(out, tag + '.dyn')
    with open(os.path.join(out, 'ic.dyn')) as fin, open(dyn_out, 'w') as fout, \
            open(os.path.join(out, tag + '.log'), 'w') as flog:
        wall, code = run_timed(command, stdin=fin, stdout=fout, stderr=flog)
    result = dict(code='kira', eta=eta, wall_s=wall, ok=(code == 0), command=command)
    snaps = [s.to_physical(units) for s in read_dyn(dyn_out)] if code == 0 else []
    return result, snaps


def run_hendos(out, icfile, hermite, tend, dt, eta, threads):
    """One HENDOS run. Returns a dictionary of results, and the snapshots."""
    tag = 'hendos_eta%g' % eta
    rundir = os.path.join(out, tag)
    os.makedirs(rundir, exist_ok=True)
    for name in os.listdir(os.path.join(rundir, 'data')) if os.path.isdir(
            os.path.join(rundir, 'data')) else []:
        os.remove(os.path.join(rundir, 'data', name))     # HENDOS will not overwrite a run
    with open(os.path.join(rundir, 'params.nml'), 'w') as f:
        f.write(HENDOS_PARAMS.format(icfile=os.path.abspath(icfile), tend=tend, dt=dt, eta=eta))
    env = dict(os.environ, OMP_NUM_THREADS=str(threads))
    command = HERMITE_CMD.format(hermite=os.path.abspath(hermite))
    with open(os.path.join(rundir, 'hermite.log'), 'w') as flog:
        wall, code = run_timed(command, cwd=rundir, stdout=flog, stderr=subprocess.STDOUT,
                               env=env)
    result = dict(code='HENDOS', eta=eta, wall_s=wall, ok=(code == 0), command=command,
                  threads=threads)
    snaps = []
    if code == 0:
        run = Run(os.path.join(rundir, 'data'), 'bench')
        snaps = [run.snapshot(i) for i in range(run.nsnap)]
        if snaps:
            result['own_energy_error'] = snaps[-1].value('energy_error')
    return result, snaps


def finish(result, snaps, e0, tend):
    """Add the energy error and end time to a result, from the last snapshot."""
    if not snaps:
        result['ok'] = False
        return
    last = snaps[-1]
    e1 = total_energy(last.m, last.pos, last.vel, G_PHYS)
    result['energy_error'] = abs(e1 - e0) / abs(e0)
    result['t_end_Myr'] = float(last.time)
    result['n_end'] = int(last.n)
    result['reached_end'] = bool(abs(last.time - tend) < 1e-3 * tend)


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

    # wall-clock time against energy error
    fig, ax = plt.subplots(figsize=(7.5, 5.5), facecolor=SURFACE)
    for code, colour, marker in (('HENDOS', BLUE, 'o'), ('kira', ORANGE, 's')):
        rows = sorted((r for r in results if r['code'] == code and r['ok']
                       and r.get('energy_error', 0) > 0), key=lambda r: r['energy_error'])
        if not rows:
            continue
        x = [r['energy_error'] for r in rows]
        y = [r['wall_s'] for r in rows]
        ax.plot(x, y, color=colour, lw=2, marker=marker, ms=8, mec=SURFACE, mew=1.5,
                label=code)
        ax.annotate('  ' + code, (x[-1], y[-1]), color=TEXT2, va='center', fontsize=10)
        for r in rows:
            ax.annotate('eta=%g' % r['eta'], (r['energy_error'], r['wall_s']),
                        textcoords='offset points', xytext=(0, 9), ha='center',
                        fontsize=8, color=TEXT2)
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_xlabel('relative energy error at the end', color=TEXT2)
    ax.set_ylabel('wall-clock time (s)', color=TEXT2)
    ax.set_title(title, color=TEXT, loc='left', fontsize=11)
    ax.legend(frameon=False, labelcolor=TEXT)
    style(ax)
    fig.tight_layout()
    fig.savefig(os.path.join(out, 'bench_time_error.png'), dpi=150, facecolor=SURFACE)

    # half-mass radius against time, most accurate run of each code
    fig, ax = plt.subplots(figsize=(7.5, 5.0), facecolor=SURFACE)
    for code, colour in (('HENDOS', BLUE), ('kira', ORANGE)):
        if code in curves:
            t, rh, eta = curves[code]
            ax.plot(t, rh, color=colour, lw=2, label='%s (eta=%g)' % (code, eta))
            ax.annotate('  ' + code, (t[-1], rh[-1]), color=TEXT2, va='center', fontsize=10)
    ax.set_xlabel('time (Myr)', color=TEXT2)
    ax.set_ylabel('half-mass radius (pc)', color=TEXT2)
    ax.set_title(title, color=TEXT, loc='left', fontsize=11)
    ax.legend(frameon=False, labelcolor=TEXT)
    style(ax)
    fig.tight_layout()
    fig.savefig(os.path.join(out, 'bench_halfmass.png'), dpi=150, facecolor=SURFACE)


def main():
    p = argparse.ArgumentParser(description='Benchmark HENDOS against Starlab kira.')
    p.add_argument('icfile', help='HENDOS initial conditions file (from makeics)')
    p.add_argument('--tend', type=float, required=True, help='end time (Myr)')
    p.add_argument('--nout', type=int, default=20,
                   help='number of outputs and snapshots over the run [20]')
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

    results, curves = [], {}
    jobs = []
    if not a.skip_hendos:
        jobs += [('HENDOS', eta) for eta in a.eta_hendos]
    if not a.skip_kira:
        jobs += [('kira', eta) for eta in a.eta_kira]
    for code, eta in jobs:
        print('running %s with eta = %g ...' % (code, eta), flush=True)
        if code == 'kira':
            result, snaps = run_kira(out, units, a.tend, dt, eta)
        else:
            result, snaps = run_hendos(out, a.icfile, a.hermite, a.tend, dt, eta, a.threads)
        finish(result, snaps, e0, a.tend)
        results.append(result)
        if result['ok']:
            print('   %.2f s, energy error %.2e' % (result['wall_s'], result['energy_error']))
            best = curves.get(code)
            if best is None or eta < best[2]:
                t, rh = half_mass_curve(snaps)
                curves[code] = (t, rh, eta)
        else:
            print('   FAILED - see the log in %s' % out)

    title = '%s   N = %d   to %g Myr' % (os.path.basename(a.icfile), ic.n, a.tend)
    with open(os.path.join(out, 'results.json'), 'w') as f:
        json.dump(dict(icfile=a.icfile, n=int(ic.n), tend_Myr=a.tend, results=results), f, indent=1)
    lines = [title, '',
             '%-8s %8s %12s %14s %12s  %s' % ('code', 'eta', 'wall (s)', 'energy error',
                                             'end (Myr)', 'notes')]
    for r in results:
        if r['ok']:
            notes = '' if r.get('reached_end') else 'DID NOT REACH tend'
            if r['code'] == 'HENDOS':
                notes += ' %d thread(s)' % r['threads']
            lines.append('%-8s %8g %12.2f %14.3e %12.4f  %s'
                         % (r['code'], r['eta'], r['wall_s'], r['energy_error'],
                            r['t_end_Myr'], notes.strip()))
        else:
            lines.append('%-8s %8g %12s %14s %12s  FAILED' % (r['code'], r['eta'], '-', '-', '-'))
    text = '\n'.join(lines)
    with open(os.path.join(out, 'results.txt'), 'w') as f:
        f.write(text + '\n')
    print('\n' + text)
    if any(r['ok'] for r in results):
        make_plots(out, results, curves, title)
        print('\nplots and tables are in %s' % out)


if __name__ == '__main__':
    main()
