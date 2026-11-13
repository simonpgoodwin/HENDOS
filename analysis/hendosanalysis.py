"""
HENDOS - HErmite N-body Dynamics Of Systems
Analysis of snapshots.

Each analysis is a function that takes one Snapshot (see hendosio.py) and
returns one number, or a few numbers with labels. timeseries() runs one
over every snapshot of a run. The viewer (hendosview.py) lists everything
registered here, so a new analysis appears in it with no other change.

TO ADD YOUR OWN ANALYSIS: write a function below (or in your own file that
imports this one) and put @quantity(...) above it:

    @quantity('mean mass', 'Mean stellar mass', 'Msun')
    def mean_mass(snap):
        return snap.m.mean()

For several curves on one plot return a dictionary, {label: number}:

    @quantity('mass in halves', 'Mass inside and outside 1 pc', 'Msun')
    def halves(snap):
        r = radius(snap)
        return {'inside': snap.m[r < 1].sum(), 'outside': snap.m[r >= 1].sum()}

Return None (or nan) if the snapshot cannot give the number.

Units are those of the files: Msun, pc, km/s, Myr.
"""
import numpy as np

QUANTITIES = {}      # name -> (function, axis label, unit)


def quantity(name, label, unit=''):
    """Register an analysis so that timeseries() and the viewer can use it."""
    def register(func):
        QUANTITIES[name] = (func, label, unit)
        return func
    return register


# ----------------------------------------------------------------------
# tools
# ----------------------------------------------------------------------
def centre_of_mass(snap):
    """Centre of mass position and velocity."""
    mtot = snap.m.sum()
    return (snap.m[:, None] * snap.pos).sum(0) / mtot, \
           (snap.m[:, None] * snap.vel).sum(0) / mtot


def radius(snap):
    """Distance of each star from the centre of mass (pc)."""
    com, _ = centre_of_mass(snap)
    return np.sqrt(((snap.pos - com) ** 2).sum(1))


def speed(snap):
    """Speed of each star relative to the centre of mass (km/s)."""
    _, cov = centre_of_mass(snap)
    return np.sqrt(((snap.vel - cov) ** 2).sum(1))


def lagrangian_radii(snap, fractions):
    """Radii about the centre of mass holding the given fractions of the mass."""
    r = radius(snap)
    order = np.argsort(r)
    cum = np.cumsum(snap.m[order]) / snap.m.sum()
    return np.array([r[order][min(np.searchsorted(cum, f), snap.n - 1)]
                     for f in fractions])


def specific_energy(snap):
    """Energy per unit mass of each star, 0.5 v^2 + phi ((km/s)^2).

    v is relative to the centre of mass. None if the file has no phi
    column (initial conditions files do not).
    """
    if snap.phi is None:
        return None
    return 0.5 * speed(snap) ** 2 + snap.phi


# ----------------------------------------------------------------------
# the analyses
# ----------------------------------------------------------------------
@quantity('half-mass radius', 'Half-mass radius', 'pc')
def half_mass_radius(snap):
    return lagrangian_radii(snap, [0.5])[0]


@quantity('Lagrangian radii', 'Lagrangian radius', 'pc')
def lagrangian_10_50_90(snap):
    r10, r50, r90 = lagrangian_radii(snap, [0.1, 0.5, 0.9])
    return {'10%': r10, '50%': r50, '90%': r90}


@quantity('velocity dispersion', '3D velocity dispersion', 'km/s')
def velocity_dispersion(snap):
    return np.sqrt((snap.m * speed(snap) ** 2).sum() / snap.m.sum())


@quantity('virial ratio', 'Virial ratio', '')
def virial_ratio(snap):
    return snap.value('virial_ratio')


@quantity('energy error', 'Relative energy error', '')
def energy_error(snap):
    return snap.value('energy_error')


@quantity('bound fraction', 'Fraction of stars that are bound', '')
def bound_fraction(snap):
    e = specific_energy(snap)
    if e is None:
        return None
    return np.mean(e < 0)


@quantity('bound mass', 'Mass in bound stars', 'Msun')
def bound_mass(snap):
    e = specific_energy(snap)
    if e is None:
        return None
    return snap.m[e < 0].sum()


@quantity('furthest star', 'Distance of the furthest star', 'pc')
def furthest_star(snap):
    return radius(snap).max()


@quantity('Q parameter', 'Q parameter (x-y projection)', '')
def q_parameter(snap):
    """Q of Cartwright & Whitworth (2004, MNRAS, 348, 589) in the x-y plane.

    Q = mbar / sbar:  mbar is the mean edge length of the minimum spanning
    tree normalised by sqrt(N A)/(N-1), sbar the mean separation of the
    stars divided by the radius R, where R is the distance of the furthest
    star from the mean position and A = pi R^2.
    Needs scipy; returns None without it.
    """
    try:
        from scipy.sparse.csgraph import minimum_spanning_tree
        from scipy.spatial.distance import pdist, squareform
    except ImportError:
        return None
    xy = snap.pos[:, :2]
    n = len(xy)
    if n < 3:
        return None
    d = pdist(xy)
    rmax = np.sqrt(((xy - xy.mean(0)) ** 2).sum(1)).max()
    area = np.pi * rmax ** 2
    mst_length = minimum_spanning_tree(squareform(d)).sum()
    mbar = mst_length / np.sqrt(n * area)
    sbar = d.mean() / rmax
    return mbar / sbar


# ----------------------------------------------------------------------
def timeseries(run, name, progress=None):
    """Run the analysis called name over every snapshot of a run.

    Returns (times, curves): times in Myr, and curves a dictionary
    {label: array of values} with one entry per curve (the label is ''
    for an analysis that returns a single number).
    progress, if given, is called as progress(i, nsnap) after each snapshot.
    """
    func = QUANTITIES[name][0]
    curves = {}
    for i in range(run.nsnap):
        result = func(run.snapshot(i))
        if not isinstance(result, dict):
            result = {'': result}
        for label, value in result.items():
            if label not in curves:
                curves[label] = np.full(run.nsnap, np.nan)
            if value is not None:
                curves[label][i] = value
        if progress is not None:
            progress(i + 1, run.nsnap)
    return run.times.copy(), curves


# ----------------------------------------------------------------------
# Ensembles (runs that differ only in the seed of their initial conditions)
ENSEMBLE_STATS = ('n', 'mean', 'std', 'median', 'p16', 'p84', 'min', 'max')


def ensemble_timeseries(ensemble, name, progress=None):
    """Run the analysis called name over every snapshot of every member of
    an ensemble (see hendosio.Ensemble) and put the members together.

    Returns (times, members, stats):
      times     the output times (Myr) - those of the member with the most
                snapshots; a member that stopped early, or is missing a
                snapshot, has nan at the times it does not have
      members   {label: array [number of members, number of times]} with
                one row per member, in the order of ensemble.members
      stats     {label: {statistic: array over time}} for the statistics in
                ENSEMBLE_STATS: n (members with a value at that time), mean,
                std (standard deviation between members), median, p16 and
                p84 (16th and 84th percentiles: the middle 68%), min, max
    progress, if given, is called as progress(i, number of members).
    """
    series = []
    for i, run in enumerate(ensemble.runs):
        series.append(timeseries(run, name))
        if progress is not None:
            progress(i + 1, len(ensemble.runs))
    if not series:
        return np.array([]), {}, {}
    times = max((t for t, _ in series), key=len)
    labels = []
    for _, curves in series:
        for label in curves:
            if label not in labels:
                labels.append(label)
    members = {label: np.full((len(series), len(times)), np.nan) for label in labels}
    for k, (t, curves) in enumerate(series):
        # match each of this member's times to the common list of times
        idx = np.searchsorted(times, t - 1e-9)
        ok = (idx < len(times))
        ok[ok] &= np.abs(times[idx[ok]] - t[ok]) < 1e-6
        for label, values in curves.items():
            members[label][k, idx[ok]] = values[ok]
    stats = {label: ensemble_stats(values) for label, values in members.items()}
    return times, members, stats


def ensemble_stats(values):
    """Statistics over the members (rows) of an array [members, times],
    ignoring nan. Returns a dictionary of arrays over time."""
    import warnings
    values = np.asarray(values, float)
    count = np.sum(np.isfinite(values), axis=0)
    with warnings.catch_warnings():
        warnings.simplefilter('ignore', RuntimeWarning)     # times with no values
        out = {
            'n': count.astype(float),
            'mean': np.nanmean(values, axis=0),
            'std': np.nanstd(values, axis=0, ddof=0),
            'median': np.nanmedian(values, axis=0),
            'p16': np.nanpercentile(values, 16, axis=0),
            'p84': np.nanpercentile(values, 84, axis=0),
            'min': np.nanmin(values, axis=0),
            'max': np.nanmax(values, axis=0),
        }
    # the scatter between members is the sample standard deviation
    with np.errstate(invalid='ignore', divide='ignore'):
        out['std'] = np.where(count > 1, out['std'] * np.sqrt(count / np.maximum(count - 1, 1)),
                              np.nan)
    return out
