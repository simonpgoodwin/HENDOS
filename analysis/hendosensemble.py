#!/usr/bin/env python3
"""
HENDOS - HErmite N-body Dynamics Of Systems
Collate an ensemble: put together the results of runs that differ only in
the random number seed of their initial conditions.

    python hendosensemble.py <data directory>              list the ensembles there
    python hendosensemble.py <data directory> <ensemble>   collate one
    python hendosensemble.py <data directory> --all        collate every one

An ensemble is found from the names of its files: hermite, run with
nensemble > 1, gives the members the runids <ensemble>_e001, <ensemble>_e002,
... so no list of members is needed, and members that are still running or
missing are simply left out (and reported).

For the ensemble it writes, in the data directory (or --out):

  <ensemble>_ens_summary.dat       one line per member: end time, number of
                                   stars, energy error, mergers, binaries and
                                   groups, and every analysis at the end
  <ensemble>_ens_<analysis>.dat    for every analysis in hendosanalysis.py,
                                   against time: number of members, mean,
                                   standard deviation, median, 16th and 84th
                                   percentiles, minimum and maximum
  <ensemble>_ens_<analysis>.png    the same as a plot: every member faintly,
                                   the median, and the band holding the
                                   middle 68% of members
  <ensemble>_ens_mergers.dat       the merger files of all members in one,
                                   with the member number as the first column
  <ensemble>_ens_lagrangian.dat    the same statistics for the run-time
                                   Lagrangian radii files, if there are any

Any analysis you add to hendosanalysis.py with @quantity is collated too,
with no change here. In a script:

    from hendosio import Ensemble
    import hendosanalysis as ha
    ens = Ensemble('data', 'test01')
    times, members, stats = ha.ensemble_timeseries(ens, 'half-mass radius')
    stats['']['median']          # median half-mass radius against time

Needs numpy, and matplotlib for the plots (--noplots to do without).
"""
import argparse
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import hendosanalysis as ha                                   # noqa: E402
from hendosio import Ensemble, find_ensembles, read_header   # noqa: E402

# header keys of a member's last snapshot that go in the summary (if present)
SUMMARY_KEYS = ('energy_error', 'mergers', 'groups_as_one_body', 'binaries_closed',
                'few_body_groups_formed', 'triples_frozen', 'wall_time_s')

# one hue for the ensemble: the median, a band, and the members in grey
BLUE, BAND, MEMBER = '#2a78d6', '#cde2fb', '#b0afa8'
TEXT, TEXT2, GRID, SURFACE = '#0b0b0b', '#52514e', '#dddcd6', '#fcfcfb'


def slug(name):
    """An analysis name as part of a file name (no spaces)."""
    return ''.join(c if c.isalnum() else '_' for c in name).strip('_').lower()


def write_table(path, keys, columns, units, rows):
    """A text table with a '# key = value' header, like every HENDOS file."""
    with open(path, 'w') as f:
        for key, value in keys:
            f.write('# %-29s = %s\n' % (key, value))
        f.write('# %-29s = %s\n' % ('columns', ' '.join(columns)))
        f.write('# %-29s = %s\n' % ('column_units', ' '.join(units)))
        f.write('# end_header\n')
        for row in rows:
            f.write(' '.join(row) + '\n')


def fmt(x):
    return '%16.8e' % x


def common_keys(ens):
    return [('format', 'NBODY-ENSEMBLE 1'), ('creator', 'HENDOS (hendosensemble.py)'),
            ('ensemble', ens.name), ('members_found', len(ens)),
            ('ensemble_size', ens.size if ens.size else 'unknown'),
            ('members_missing', ' '.join(str(k) for k in ens.missing) or 'none'),
            ('units', 'mass Msun, length pc, velocity km/s, time Myr')]


def collate_summary(ens, out):
    """One line per member, from its last snapshot."""
    names = list(ha.QUANTITIES)
    columns = ['member', 'snapshots', 't_end', 'n'] + list(SUMMARY_KEYS)
    units = ['-', '-', 'Myr', '-'] + ['-'] * len(SUMMARY_KEYS)
    for name in names:
        func, label, unit = ha.QUANTITIES[name]
        result = func(ens.runs[0].snapshot(ens.runs[0].nsnap - 1))
        labels = list(result) if isinstance(result, dict) else ['']
        for lab in labels:
            columns.append(slug(name) + ('_' + slug(lab) if lab else ''))
            units.append(unit or '-')
    rows = []
    for member, run in zip(ens.members, ens.runs):
        last = run.snapshot(run.nsnap - 1)
        row = ['%6d' % member, '%6d' % run.nsnap, fmt(last.time), '%8d' % last.n]
        for key in SUMMARY_KEYS:
            try:
                row.append(fmt(float(last.header[key])))
            except (KeyError, ValueError):
                row.append(fmt(np.nan))
        for name in names:
            result = ha.QUANTITIES[name][0](last)
            values = result.values() if isinstance(result, dict) else [result]
            row += [fmt(np.nan if v is None else v) for v in values]
        rows.append(row)
    path = os.path.join(out, ens.name + '_ens_summary.dat')
    write_table(path, common_keys(ens) + [('contents', 'each member at its last snapshot')],
                columns, units, rows)
    return path


def collate_quantity(ens, name, out, plots):
    """Statistics against time of one analysis; returns the files written."""
    func, label, unit = ha.QUANTITIES[name]
    times, members, stats = ha.ensemble_timeseries(ens, name)
    columns, units, data = ['time'], ['Myr'], [times]
    for lab in stats:
        pre = slug(lab) + '_' if lab else ''
        for stat in ha.ENSEMBLE_STATS:
            columns.append(pre + stat)
            units.append('-' if stat == 'n' else (unit or '-'))
            data.append(stats[lab][stat])
    rows = [[fmt(col[i]) for col in data] for i in range(len(times))]
    base = os.path.join(out, '%s_ens_%s' % (ens.name, slug(name)))
    write_table(base + '.dat', common_keys(ens) + [
        ('analysis', name),
        ('statistics', 'over the members at each time: n = members with a value, std = '
                       'standard deviation between members, p16 and p84 = 16th and 84th '
                       'percentiles')], columns, units, rows)
    written = [base + '.dat']
    if plots:
        plot_quantity(ens, name, label, unit, times, members, stats, base + '.png')
        written.append(base + '.png')
    return written


def plot_quantity(ens, name, label, unit, times, members, stats, path):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt

    labels = list(stats)
    fig, axes = plt.subplots(1, len(labels), figsize=(6.4 * len(labels), 4.8),
                             facecolor=SURFACE, squeeze=False, sharey=True)
    for ax, lab in zip(axes[0], labels):
        for row in members[lab]:
            ax.plot(times, row, color=MEMBER, lw=0.7, alpha=0.7, zorder=1)
        st = stats[lab]
        ax.fill_between(times, st['p16'], st['p84'], color=BAND, zorder=2, lw=0)
        ax.plot(times, st['median'], color=BLUE, lw=2.2, zorder=3)
        ax.set_facecolor(SURFACE)
        ax.grid(True, color=GRID, lw=0.6)
        ax.set_axisbelow(True)
        for side in ('top', 'right'):
            ax.spines[side].set_visible(False)
        for side in ('left', 'bottom'):
            ax.spines[side].set_color(GRID)
        ax.tick_params(colors=TEXT2)
        ax.set_xlabel('time (Myr)', color=TEXT2)
        if lab:
            ax.set_title(lab, color=TEXT2, fontsize=10, loc='left')
    axes[0][0].set_ylabel(label + (' (%s)' % unit if unit else ''), color=TEXT2)
    fig.suptitle('%s   ensemble of %d   blue line: median, band: middle 68%% of members, '
                 'grey: each member' % (ens.name, len(ens)), color=TEXT, fontsize=10,
                 x=0.01, ha='left')
    fig.tight_layout(rect=(0, 0, 1, 0.95))
    fig.savefig(path, dpi=140, facecolor=SURFACE)
    plt.close(fig)


def collate_mergers(ens, out):
    """All the members' merger files in one, member number first."""
    lines, header = [], None
    for member, run in zip(ens.members, ens.runs):
        path = os.path.join(ens.datadir, run.runid + '_mergers.dat')
        if not os.path.exists(path):
            continue
        if header is None:
            header = read_header(path)[0]
        for line in open(path):
            if line.strip() and not line.startswith('#'):
                lines.append('%6d %s' % (member, line.rstrip('\n')))
    if header is None:
        return None
    path = os.path.join(out, ens.name + '_ens_mergers.dat')
    with open(path, 'w') as f:
        for key, value in common_keys(ens) + [('mergers', len(lines))]:
            f.write('# %-29s = %s\n' % (key, value))
        f.write('# %-29s = member %s\n' % ('columns', header.get('columns', '')))
        f.write('# %-29s = - %s\n' % ('column_units', header.get('column_units', '')))
        f.write('# end_header\n')
        f.write('\n'.join(lines) + ('\n' if lines else ''))
    return path


def collate_runtime(ens, out):
    """Statistics of the run-time analysis files (e.g. Lagrangian radii)."""
    written = []
    names = sorted(set(n for run in ens.runs for n in run.analysis_files))
    for name in names:
        files = [run.analysis(name) for run in ens.runs if name in run.analysis_files]
        cols = [c for c in files[0].columns if c != 'time']
        times = max((f.column('time') for f in files), key=len)
        columns, units, data = ['time'], ['Myr'], [times]
        for col in cols:
            values = np.full((len(files), len(times)), np.nan)
            for k, f in enumerate(files):
                t = f.column('time')
                idx = np.searchsorted(times, t - 1e-9)
                ok = idx < len(times)
                ok[ok] &= np.abs(times[idx[ok]] - t[ok]) < 1e-6
                values[k, idx[ok]] = f.column(col)[ok]
            st = ha.ensemble_stats(values)
            for stat in ha.ENSEMBLE_STATS:
                columns.append('%s_%s' % (col, stat))
                units.append('-')
                data.append(st[stat])
        rows = [[fmt(c[i]) for c in data] for i in range(len(times))]
        path = os.path.join(out, '%s_ens_%s.dat' % (ens.name, slug(name)))
        write_table(path, common_keys(ens) + [('run_time_analysis', name)], columns, units, rows)
        written.append(path)
    return written


def collate(datadir, name, out, plots=True, say=print):
    ens = Ensemble(datadir, name)
    if not len(ens):
        say('no ensemble called %s in %s' % (name, datadir))
        return []
    ens_runs = [(k, r) for k, r in zip(ens.members, ens.runs) if r.nsnap]
    say('ensemble %s: %d members found%s' % (
        name, len(ens), ' of %d' % ens.size if ens.size else ''))
    if ens.missing:
        say('   missing members: ' + ' '.join(str(k) for k in ens.missing))
    ends = np.array([r.times[-1] for _, r in ens_runs])
    if ends.max() - ends.min() > 1e-6:
        short = [k for (k, r), e in zip(ens_runs, ends) if e < ends.max() - 1e-6]
        say('   members that have not reached %.4f Myr: %s'
            % (ends.max(), ' '.join(str(k) for k in short)))
    os.makedirs(out, exist_ok=True)
    written = [collate_summary(ens, out)]
    for qname in ha.QUANTITIES:
        written += collate_quantity(ens, qname, out, plots)
    mergers = collate_mergers(ens, out)
    if mergers:
        written.append(mergers)
    written += collate_runtime(ens, out)
    for path in written:
        say('   wrote ' + path)
    return written


def main():
    p = argparse.ArgumentParser(description='Collate a HENDOS ensemble.')
    p.add_argument('datadir', help='the data directory of the runs')
    p.add_argument('ensemble', nargs='?', help='name of the ensemble (the runid without _e001)')
    p.add_argument('--all', action='store_true', help='collate every ensemble found')
    p.add_argument('--out', default=None, help='directory for the files [the data directory]')
    p.add_argument('--noplots', action='store_true', help='tables only')
    a = p.parse_args()

    found = find_ensembles(a.datadir)
    if not found:
        print('no ensembles found in %s (members have runids ending _e001, _e002, ...)'
              % a.datadir)
        return 1
    if not a.ensemble and not a.all:
        print('ensembles in %s:' % a.datadir)
        for name, members in found.items():
            print('   %-28s %d members (%d to %d)' % (name, len(members), members[0],
                                                     members[-1]))
        print('collate one with:  python hendosensemble.py %s <name>' % a.datadir)
        return 0
    names = list(found) if a.all else [a.ensemble]
    for name in names:
        collate(a.datadir, name, a.out or a.datadir, plots=not a.noplots)
    return 0


if __name__ == '__main__':
    sys.exit(main())
