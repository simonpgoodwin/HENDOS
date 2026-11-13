# HENDOS analysis and viewer

Python tools for looking at a run: a viewer with buttons, and the modules it is
built from, which can be used on their own in scripts.

**Keep this file up to date: if the viewer's controls, the analyses or the
functions listed here change, change this in the same edit.**

| File | What it is |
|---|---|
| `hendosview.py` | the viewer: 3D view, plots against time, images and movies |
| `hendosanalysis.py` | the analyses done on snapshots - add your own here |
| `hendosio.py` | reads snapshots and run-time analysis files; finds runs and ensembles |
| `hendosensemble.py` | puts the results of an ensemble of runs together |
| `maketestdata.py` | makes a few test files to try the viewer on (not a HENDOS run); with a number, a test ensemble too |

Needs Python 3 with `numpy` and `matplotlib`. `scipy` is needed only for the Q
parameter, and `ffmpeg` only for mp4 movies (without it an animated GIF is
written). Nothing else to install: the window uses matplotlib's own buttons.

## Starting the viewer
```
cd HENDOS/run
python ../analysis/hendosview.py data            (first run found in data/)
python ../analysis/hendosview.py data test01     (the run test01)
```
The first argument is the data directory (`datadir` in `params.nml`), the
second the `runid`. To try it before you have a run:
```
python ../analysis/maketestdata.py testdata
python ../analysis/hendosview.py testdata
```
The test files are a crude toy integration, marked as test data in their
headers. They are only there to have something to look at.

## The window
Controls are down the left; the picture is on the right; the slider along the
bottom moves through the snapshots.

### VIEW
- **3D view** - the stars of one snapshot. Drag with the mouse to rotate.
- **Snapshot analysis** - a quantity worked out from every snapshot, plotted
  against time. The first time a quantity is chosen every snapshot is read,
  which takes a moment for a big run; after that it is kept.
- **Run-time files** - the files the code wrote while it ran
  (`<runid>_ana_<name>.dat`). Choose the file on the left and tick the
  columns on the right. Columns with different units are drawn in separate
  panels, one above the other, so each has its own axis.

In both plot views the dashed line marks the time of the current snapshot.

### COLOUR BY (3D view)
| Choice | Stars are coloured by | Key |
|---|---|---|
| mass | mass, on a logarithmic scale | colour bar |
| mass bins | mass in four bins: below 0.5, 0.5-2, 2-8, above 8 Msun | legend with numbers |
| type | particle type; each type also has its own marker shape | legend with numbers |
| speed | speed relative to the centre of mass | colour bar |
| bound / unbound | sign of kinetic + potential energy per unit mass | legend with numbers |
| timestep | current timestep, on a logarithmic scale | colour bar |

Mass and type can be shown together: colour by one and use the options for
the other (colour by mass with **shape by type**, or colour by type with
**size by mass**).

### OPTIONS
- **size by mass** - bigger points for more massive stars.
- **shape by type** - a different marker for each particle type (always on
  when colouring by type).
- **ring 10 most massive** - draws a ring round the ten most massive stars.
- **fixed scale** - the box stays the same size for the whole run (four times
  the initial half-mass radius; stars outside it are not drawn). Off: the box
  is refitted to each snapshot to hold 98 per cent of the stars.
- **centre on c.o.m.** - keeps the centre of mass in the middle.
- **rotate when playing** - turns the view slowly during Play and in movies.
- **dark background**.
- **log y axis (plots)** - logarithmic y axis in the two plot views.

### Buttons
- **Play / Pause** - steps through the snapshots.
- **Next run** - the next run in the same data directory.
- **Zoom in / Zoom out**.
- **Save image** - what is showing, as a PNG in the data directory, without
  the buttons: `<runid>_view_<time>Myr.png` or `<runid>_plot_<name>.png`.
- **Save movie** - the 3D view through every snapshot with the present colours,
  options and viewing angle: `<runid>_movie.mp4` (or `.gif`). The window does
  not respond while it is being made; progress is shown at the bottom left.

## The analyses
In `hendosanalysis.py`. Each takes one snapshot and returns a number (or a few
labelled numbers), and all of them appear in the viewer's QUANTITY list.

| Name | What |
|---|---|
| half-mass radius | radius about the centre of mass holding half the mass |
| Lagrangian radii | the 10, 50 and 90 per cent radii |
| velocity dispersion | mass-weighted 3D velocity dispersion |
| virial ratio | from the snapshot header |
| energy error | relative energy error, from the snapshot header |
| bound fraction | fraction of stars with negative energy (uses the `phi` column) |
| bound mass | mass in those stars |
| furthest star | distance of the furthest star from the centre of mass |
| Q parameter | Q of Cartwright & Whitworth (2004) in the x-y projection (needs scipy) |

### Adding your own
Write a function in `hendosanalysis.py` with `@quantity(...)` above it. It
then appears in the viewer with nothing else to change.
```python
@quantity('mean mass', 'Mean stellar mass', 'Msun')
def mean_mass(snap):
    return snap.m.mean()
```
`snap` has `time`, `n`, `id`, `m`, `pos` (n by 3), `vel` (n by 3), `type`,
`phi`, `dt` and `header`. Return a dictionary `{label: number}` to get several
curves on one plot. Useful tools already there: `centre_of_mass`, `radius`,
`speed`, `lagrangian_radii`, `specific_energy`.

## Ensembles
An ensemble is a set of runs that differ only in the random number seed of
their initial conditions (`nensemble` in `ics.nml` and `params.nml`). hermite
gives the members the runids `<name>_e001`, `<name>_e002`, ... and the tools
here recognise an ensemble from those names: nothing has to be listed, and
members that are missing or still running are left out and reported.

```
python hendosensemble.py data              lists the ensembles in data/
python hendosensemble.py data test01       collates the ensemble test01
python hendosensemble.py data --all        collates every ensemble
```
Options: `--out <directory>` for the files (the data directory otherwise) and
`--noplots` for tables only. For the ensemble `test01` it writes:

| File | What |
|---|---|
| `test01_ens_summary.dat` | one line per member: end time, number of stars, energy error, mergers, binaries and groups, and the value of every analysis at the end |
| `test01_ens_<analysis>.dat` | for each analysis below, against time: the number of members, their mean, the standard deviation between them, median, 16th and 84th percentiles, minimum and maximum |
| `test01_ens_<analysis>.png` | the same as a plot: each member in grey, the median, and the band that holds the middle 68% of the members |
| `test01_ens_mergers.dat` | the merger files of all the members in one, with the member number as the first column |
| `test01_ens_lagrangian.dat` | the same statistics for the run-time Lagrangian radii files |

The tables have the usual `# key = value` header, so `hendosio.read_analysis`
reads them. Every analysis registered with `@quantity` is collated, your own
included, with no change to `hendosensemble.py`. Statistics at each time use
the members that have a snapshot at that time, and the column `n` says how
many that was, so an ensemble can be collated while it is still running.

In the viewer the members appear in the list of runs as separate runs
(`test01_e001`, ...) and can be looked at one at a time; the viewer does not
yet show ensemble statistics.

In a script:
```python
from hendosio import Ensemble, find_ensembles
import hendosanalysis as ha

print(find_ensembles('data'))                 # {'test01': [1, 2, 3, ...]}
ens = Ensemble('data', 'test01')
print(len(ens), ens.size, ens.missing)        # found, asked for, missing
times, members, stats = ha.ensemble_timeseries(ens, 'half-mass radius')
stats['']['median'], stats['']['p16'], stats['']['p84']    # against time
members['']                                   # array [member, time]
run7 = ens.run(7)                             # one member, an ordinary Run
```
For an analysis with several curves (`'Lagrangian radii'`) the dictionaries
have one entry per curve (`stats['50%']['mean']`).

## Using the pieces in a script
```python
from hendosio import Run
import hendosanalysis as ha

run = Run('data', 'test01')
for i in range(run.nsnap):
    snap = run.snapshot(i)
    print(snap.time, ha.half_mass_radius(snap))

t, curves = ha.timeseries(run, 'Lagrangian radii')     # every snapshot
lag = run.analysis('lagrangian')                       # a run-time file
print(lag.columns, lag.data.shape)
```
`hendosview.py` also has two drawing functions that need no window:
`draw_particles(ax, snap, colour=..., ...)` draws one snapshot on a matplotlib
3D axis, and `plot_curves(fig, rect, groups, ...)` draws curves against time.

## Changing how things look
At the top of `hendosview.py`:
- `TYPE_STYLE` - the colour and marker of each particle type (add a line for a
  new type);
- `MASS_BINS` - the edges of the mass bins;
- `LIGHT` and `DARK` - the colours for the two backgrounds.

Only three hues are used for separate things, with neutral grey for ordinary
stars, because more than that cannot be told apart reliably in a scatter of
points by colour-blind readers; that is why types have marker shapes as well.

## What has and has not been tested
The reader, the analyses, every view, the options, saving images and saving a
movie were run on the test files without a screen. `hendosensemble.py` was run
on a test ensemble of five made by `maketestdata.py` (one member cut short,
one with a merger file), not on a real ensemble. The window has not been
used interactively (dragging to rotate, the Play timer, clicking), because
there was no display to try it on.
