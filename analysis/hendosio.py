"""
HENDOS - HErmite N-body Dynamics Of Systems
Reading HENDOS output files.

    read_snapshot(path)    one snapshot (or initial conditions) file
    read_analysis(path)    one run-time analysis file  <runid>_ana_<name>.dat
    find_runs(datadir)     the runs that have snapshots in a directory
    Run(datadir, runid)    all the files of one run, snapshots read when asked for

Needs only numpy. Use it on its own in your own scripts:

    from hendosio import Run
    run = Run('data', 'test01')
    for i in range(run.nsnap):
        s = run.snapshot(i)
        print(s.time, s.m.sum(), s.pos.shape)

Everything is in the units of the files: Msun, pc, km/s, Myr
(timesteps in yr, potentials and energies in (km/s)^2).
"""
import glob
import os
import re

import numpy as np

SNAP_RE = re.compile(r'^(?P<runid>.+)_snap_(?P<time>\d+\.\d+)Myr\.dat$')
ANA_RE = re.compile(r'^(?P<runid>.+)_ana_(?P<name>.+)\.dat$')


# ----------------------------------------------------------------------
def read_header(path):
    """The '# key = value' header of a HENDOS file.

    Returns (header, nlines): a dictionary of text values, and the number
    of header lines (so the data start on line nlines).
    """
    header = {}
    nlines = 0
    with open(path) as f:
        for line in f:
            if not line.lstrip().startswith('#'):
                break
            nlines += 1
            body = line.lstrip()[1:]
            if '=' in body:
                key, value = body.split('=', 1)
                header[key.strip()] = value.strip()
    return header, nlines


def _floats(header, key):
    """The numbers in header[key], or None if it is not there."""
    if key not in header:
        return None
    try:
        return np.array([float(x) for x in header[key].split()])
    except ValueError:
        return None


# ----------------------------------------------------------------------
class Snapshot:
    """One snapshot. Arrays have one entry per star, in id order.

    header   dictionary of the header lines (text)
    time     time (Myr)
    n        number of stars
    id       star identifiers
    m        masses (Msun)
    pos      positions, shape (n, 3) (pc)
    vel      velocities, shape (n, 3) (km/s)
    type     particle types (text: 'star', 'black_hole', ...)
    phi      potential per unit mass ((km/s)^2), or None if not in the file
    dt       current timestep (yr), or None if not in the file
    """

    def __init__(self, path):
        self.path = path
        self.header, nhead = read_header(path)
        columns = self.header.get('columns', '').split()
        with open(path) as f:
            for _ in range(nhead):
                f.readline()
            tokens = f.read().split()
        if not columns:
            # legacy file: m x y z vx vy vz
            columns = ['m', 'x', 'y', 'z', 'vx', 'vy', 'vz']
        ncol = len(columns)
        if len(tokens) % ncol != 0:
            raise ValueError('%s: the data do not fit the %d columns in the header'
                             % (path, ncol))
        table = np.array(tokens, dtype=object).reshape(-1, ncol)
        col = {name: table[:, k] for k, name in enumerate(columns)}

        def num(name):
            return col[name].astype(float) if name in col else None

        self.n = table.shape[0]
        self.m = num('m')
        self.pos = np.column_stack([num('x'), num('y'), num('z')])
        self.vel = np.column_stack([num('vx'), num('vy'), num('vz')])
        self.id = col['id'].astype(int) if 'id' in col else np.arange(1, self.n + 1)
        self.type = (col['type'].astype(str) if 'type' in col
                     else np.full(self.n, 'star'))
        self.phi = num('phi')
        self.dt = num('dt')
        t = _floats(self.header, 'time_Myr')
        self.time = float(t[0]) if t is not None else 0.0

    def value(self, key):
        """A number (or array of numbers) from the header, or None."""
        x = _floats(self.header, key)
        if x is None:
            return None
        return float(x[0]) if len(x) == 1 else x


def read_snapshot(path):
    """Read one snapshot or initial conditions file."""
    return Snapshot(path)


# ----------------------------------------------------------------------
class AnalysisFile:
    """One run-time analysis file.

    header    dictionary of the header lines (text)
    name      the analysis name (e.g. 'lagrangian')
    columns   column names
    units     column units (same length as columns, '' if not given)
    data      the numbers, shape (nrows, ncolumns)
    """

    def __init__(self, path):
        self.path = path
        self.header, nhead = read_header(path)
        self.name = self.header.get('analysis', os.path.basename(path))
        self.data = np.atleast_2d(np.loadtxt(path, comments='#'))
        ncol = self.data.shape[1]
        self.columns = self.header.get('columns', '').split()
        if len(self.columns) != ncol:
            self.columns = ['column%d' % (k + 1) for k in range(ncol)]
        self.units = self.header.get('column_units', '').split()
        if len(self.units) != ncol:
            self.units = [''] * ncol

    def column(self, name):
        """One column by name."""
        return self.data[:, self.columns.index(name)]


def read_analysis(path):
    """Read one run-time analysis file."""
    return AnalysisFile(path)


# ----------------------------------------------------------------------
def find_runs(datadir):
    """The runs with snapshots in datadir: a dictionary runid -> number of snapshots."""
    runs = {}
    for path in glob.glob(os.path.join(datadir, '*_snap_*Myr.dat')):
        match = SNAP_RE.match(os.path.basename(path))
        if match:
            runs[match.group('runid')] = runs.get(match.group('runid'), 0) + 1
    return dict(sorted(runs.items()))


class Run:
    """All the output files of one run.

    datadir, runid
    nsnap           number of snapshots
    times           their times (Myr), from the file names, in order
    paths           their file names, in the same order
    snapshot(i)     snapshot i (read the first time it is asked for, then kept)
    analysis_files  dictionary  analysis name -> file name
    analysis(name)  one run-time analysis file
    """

    def __init__(self, datadir, runid):
        self.datadir = datadir
        self.runid = runid
        found = []
        for path in glob.glob(os.path.join(datadir, runid + '_snap_*Myr.dat')):
            match = SNAP_RE.match(os.path.basename(path))
            if match and match.group('runid') == runid:
                found.append((float(match.group('time')), path))
        found.sort()
        self.times = np.array([t for t, _ in found])
        self.paths = [p for _, p in found]
        self.nsnap = len(found)
        self._cache = {}
        self.analysis_files = {}
        for path in sorted(glob.glob(os.path.join(datadir, runid + '_ana_*.dat'))):
            match = ANA_RE.match(os.path.basename(path))
            if match and match.group('runid') == runid:
                self.analysis_files[match.group('name')] = path
        self._analysis = {}

    def snapshot(self, i):
        if i not in self._cache:
            self._cache[i] = Snapshot(self.paths[i])
        return self._cache[i]

    def analysis(self, name):
        if name not in self._analysis:
            self._analysis[name] = AnalysisFile(self.analysis_files[name])
        return self._analysis[name]
