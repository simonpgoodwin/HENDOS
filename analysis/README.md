# HENDOS analysis and viewer

Python tools for looking at a run: a viewer with buttons, and the modules it is
built from, which can be used on their own in scripts.

**Keep this file up to date: if the viewer's controls, the analyses or the
functions listed here change, change this in the same edit.**

| File | What it is |
|---|---|
| `hendosview.py` | the viewer: 3D view, plots against time, images and movies |
| `hendosanalysis.py` | the analyses done on snapshots - add your own here |
| `hendosio.py` | reads snapshots and run-time analysis files |
| `maketestdata.py` | makes a few test files to try the viewer on (not a HENDOS run) |

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
movie were run on the test files without a screen. The window has not been
used interactively (dragging to rotate, the Play timer, clicking), because
there was no display to try it on.
