# HENDOS parameters

A full description of every parameter in the two parameter files in this
directory:

- `params.nml`, read by the N-body code `hermite` (first part of this file);
- `ics.nml`, read by the initial conditions generator `makeics` (second part,
  from "Initial conditions" onwards).

**Keep this file up to date: when a parameter is added, removed, renamed or its
standard value or meaning changes, change it here (and in `params.nml`) in the
same edit.**

## How the parameter file works
- Run with `./hermite params.nml` from `HENDOS/run`. With no file name the
  code looks for `params.nml`. Keep one file per run, e.g. `params_run01.nml`.
- The file is a set of sections. Each starts with `&name` and ends with `/`.
- Inside a section each line is `name = value`. Anything after `!` is a comment.
- Text values need quotes (`'solar'`); switches are `.true.` or `.false.`.
- The order of sections, and of parameters inside a section, does not matter.
- Any parameter left out takes its **standard value** (given below). Any whole
  section left out takes all its standard values, except `&run`, which must be
  there. In verbose mode the code lists the sections it did not find.
- A misspelt parameter name stops the run with a message naming the section.
- Units everywhere: masses in Msun, lengths in pc, velocities in km/s, times
  in Myr, unless a parameter says otherwise.
- The values in the `params.nml` supplied are an example run (a 1000-star
  Plummer sphere for 10 Myr); they are not all the standard values.

---

## `&run` - run set-up (must be present; change for every run)

### `runid` (text, no standard value - must be set)
The name of the run. Every output file name starts with it, e.g.
`test01_snap_000001.0000Myr.dat`. Only letters, digits, `-`, `_` and `.` are
allowed (no spaces), and at most 32 characters. A new run (one starting at time
0) will not start if files with this `runid` already exist in the data
directory, so an earlier run cannot be overwritten by accident: choose a new
`runid`, or move or delete the old files.

### `icfile` (text, no standard value - must be set)
The initial conditions file, with its path, e.g. `'data/plum1k.dat'`. No
spaces. Three kinds of file are understood:
- an NBODY-SNAPSHOT 2 file, as written by `makeics` and by `hermite`
  itself: a header, then lines `id m x y z vx vy vz type`;
- an older NBODY-SNAPSHOT 1 file (the same without the type): every particle
  is taken to be a star;
- a legacy file with lines `m x y z vx vy vz` (e.g. from bigfractals): ids are
  set to 1..N, every particle is a star, and the run starts at time 0.

If the file is a snapshot of an earlier run, the run **restarts** from the time
of that snapshot (the `time_Myr` in its header), and analysis files are added
to rather than replaced. The initial conditions are always moved to their
centre-of-mass frame.

### `tend` (Myr, standard value 1.0)
The time at which the run ends. It is a time, not a duration: a run restarted
from a snapshot at 5 Myr with `tend = 10` runs for another 5 Myr. It must be at
least one `tout` after the start time. The run length is rounded to a whole
number of `tout`.

### `tout` (Myr, standard value 0.1)
The interval between lines of screen output. It also sets two other things:
- it is the **longest timestep** any particle can have, and all particles are
  brought to the same time every `tout`;
- the energy is calculated every `tout`, which costs about as much as one
  force calculation on every particle, so a very small `tout` slows the run.

A sensible value is a fraction of a crossing time.

### `eps` (pc, standard value 0.0)
The softening length: the force between two particles is calculated as if
their separation were sqrt(r^2 + eps^2). 0 means no softening, which is what a
collisional calculation wants. If `eps` is greater than 0 the binary treatment
(`&binaries`) is switched off, because the Kepler orbits it uses are unsoftened.

---

### `nensemble` (whole number, standard value 1)
The number of runs to do. With more than 1 the code runs an **ensemble**: a
set of simulations that are the same in every way except the random number
seed of their initial conditions. Make the initial conditions first with
`makeics` and the same `nensemble` in `ics.nml`; that writes one file per
member, with `_e001`, `_e002`, ... put before the extension of `outfile`.
Then `./hermite params.nml`, with `icfile` set to the name *without* a member
number, runs them all:

| | set in the files | member 1 | member 2 |
|---|---|---|---|
| initial conditions | `data/plum1k.dat` | `data/plum1k_e001.dat` | `data/plum1k_e002.dat` |
| runid | `test01` | `test01_e001` | `test01_e002` |
| snapshots | | `test01_e001_snap_...Myr.dat` | `test01_e002_snap_...Myr.dat` |
| screen output | | `data/test01_e001.log` | `data/test01_e002.log` |

Every output file of a member carries its member number, because the number
is part of its runid: snapshots, analysis files, the merger file and the
restart file. The headers of its files also say which ensemble it belongs to
(`ensemble`, `ensemble_member`, `ensemble_size`). While the ensemble runs the
screen shows one line as each member finishes; what each member would have
written to the screen is in its `.log` file. `runid` must be 26 characters or
fewer.

Each member is a separate run of `hermite`, so one member stopping early (a
collision with `mergers` off, say) does not stop the others; the summary at
the end says how many finished. A single member can be run, or run again, on
its own by giving its number after the parameter file, which is also the way
to run an ensemble as a job array on a cluster:
```
./hermite params.nml 7
```
With `restart = .true.` every member is carried on from its own restart file.
To put the results of an ensemble together afterwards use
`analysis/hendosensemble.py`.

### `nparallel` (whole number, standard value 1)
How many members of an ensemble are run at the same time. Each of them uses
`OMP_NUM_THREADS` cores, so `nparallel` times `OMP_NUM_THREADS` should not be
more than the number of cores of the machine. For ensembles of small systems
the fastest way by far is one thread each and as many members at once as
there are cores:
```
export OMP_NUM_THREADS=1        (and nparallel = 16 on a 16-core machine)
```
because the members do not have to talk to each other at all, whereas the
threads of a single small run spend most of their time waiting for each other.

### `restart` (switch, standard value `.false.`)
`.true.` carries the run on from its restart file, `<runid>_restart.bin` in
the data directory, instead of starting again from `icfile`. Use the same
parameter file as the original run; `tend` may be increased to extend a run
that has finished. The run then continues exactly as it would have done:
the restart file holds everything, including timesteps, binaries and the
energy the run started with, so the energy error keeps counting from the
original start. The switches `usetree`, `useac` and `fewbody` must be the same
as when the file was written. See `trestart` for writing the file.

## `&output` - where output goes and how much

### `datadir` (text, standard value `'data'`)
The directory the snapshots and analysis files are written to. The standard
value `'data'` is `HENDOS/run/data`. Any path can be given instead:
- relative to the run directory, e.g. `'data/run01'`;
- in full, e.g. `'/scratch/me/run01'`;
- starting from your home directory, e.g. `'~/nbody/run01'`;
- `'.'` for the run directory itself.

No spaces, and at most about 250 characters. The directory is made if it does
not exist, and the run stops straight away if it cannot be written to. Note
that `datadir` is only for output: `icfile` is given with its own path, and
the screen output goes wherever you send it (e.g. `| tee data/test01.log`).

### `screen` (text, standard value `'standard'`)
How much is written to the screen.
- `'standard'`: a short header, then one line every `tout` (time, percentage
  done, energy error, virial ratio, half-mass radius, wall-clock time used and
  left), and a summary at the end. For production runs.
- `'verbose'`: also every parameter, a description of the initial conditions,
  and at every `tout` a block with the energies, energy errors, Lagrangian
  radii, centre-of-mass drift, numbers of steps, the table of timestep levels,
  binaries, tree and speed. For testing.

Warnings are written in both.

### `tsnap` (Myr, standard value 0.0)
The interval between snapshots. 0 means no snapshots. It should be a whole
multiple of `tout`; if it is not it is rounded to one, with a warning. A
snapshot is also written at the start of a new run. Each snapshot holds every
particle (`id m x y z vx vy vz type phi dt`) with a full header, and can be
used as the `icfile` of a later run.

### `dewarn` (standard value 1.0e-4)
A warning is written to the screen whenever the size of the relative energy
error since the start of the run, |dE/E|, is larger than this. It does not stop
or change the run.

---

### `trestart` (Myr, standard value 0.0)
How often a restart file is written; 0 means never. It should be a whole
number of `tout`. The file is also written when the run ends, so the run can
be extended. Each file replaces the one before, and it is written under a
temporary name first, so a crash while it is being written never loses the
previous one. Restart files are only worth having for long runs: they are in
binary, can only be read by HENDOS, and are about 250 bytes per star (more
with the neighbour scheme). To continue a run from one, set `restart`.

## `&analysis` - analysis done while the code runs

### `tana` (Myr, standard value 0.0)
The interval between analyses. 0 means none. It should be a whole multiple of
`tout` (rounded with a warning if not). Each analysis that is switched on
writes one line to its own file `<runid>_ana_<name>.dat` every `tana`.

### `lagrangian` (switch, standard value `.true.`)
Writes the radii containing 1, 2, 5, 10, 20, 30, 50, 70 and 90 per cent of the
mass, measured from the centre of mass, and the 3D velocity dispersion, to
`<runid>_ana_lagrangian.dat`. Only done if `tana` is greater than 0.

---

## `&tides` - the Galactic tidal field
The cluster is taken to be on a circular orbit in the Galaxy and the run is
done in a frame rotating with that orbit: x points away from the Galactic
centre, y in the direction of Galactic rotation, z out of the disc. Positions
and velocities, in and out, are in this frame. With a tide the energy the code
reports is the Jacobi energy, which is the conserved quantity.

### `tidal` (text, standard value `'none'`)
- `'none'`: an isolated cluster. The other parameters in this section are ignored.
- `'solar'`: the tidal field of the solar neighbourhood, made from `oorta`,
  `oortb` and `rho0`.
- `'tensor'`: your own tidal tensor `tij` and rotation rate `omegat`.

### `oorta` (km/s/kpc, standard value 15.3)
Oort's constant A (Bovy 2017). Used only with `tidal = 'solar'`.

### `oortb` (km/s/kpc, standard value -11.9)
Oort's constant B (Bovy 2017). Used only with `tidal = 'solar'`. The frame
rotates at A - B = 27.2 km/s/kpc with the standard values.

### `rho0` (Msun/pc^3, standard value 0.097)
The total local mass density of the Galactic disc - stars, gas and dark matter
(McKee, Parravano & Hollenbach 2015). It sets the strength of the tide
perpendicular to the disc. Used only with `tidal = 'solar'`.

### `tij(k,l)` ((km/s/kpc)^2, standard value 0.0 for all nine)
The components of your own tidal tensor in the rotating frame, used only with
`tidal = 'tensor'`. The tidal acceleration on a particle at position x is
T.x, and T must include the centrifugal term. Components not given are zero.
The tensor should be symmetric; if it is not, its symmetric part is used and a
warning is written. For reference the solar values are `tij(1,1) = 1664.6`,
`tij(2,2) = 0`, `tij(3,3) = -5427.5`.

### `omegat` (km/s/kpc, standard value 0.0)
The rotation rate of the frame about the z axis, which sets the Coriolis
force. Used only with `tidal = 'tensor'`. The solar value is 27.2.

---

## `&neighbours` - the Ahmad-Cohen neighbour scheme
A faster way of doing direct summation. The force on a star is split into the
part from its nearest neighbours, which changes quickly, and the part from all
the other stars, which changes slowly. The neighbour part is summed every
step; the full sum over all stars is done only every few steps, and in between
the slow part is extrapolated. It is for runs without the tree and is worth
using above a thousand stars or so. The saving grows with N.

### `useac` (switch, standard value `.false.`)
`.true.` uses the neighbour scheme. It is switched off, with a warning, if
`usetree` is on.
With fewer than 150 stars the code switches the neighbour scheme off, says so
on the screen at the start, and carries on with plain direct summation: for so
few stars a neighbour list is a large part of the whole system and the scheme
saves nothing. The value 150 is an estimate of the break-even from a model of
the code, compared at equal energy error: counting pair interactions the
scheme costs 1.36 times direct summation at 30 stars, 0.98 at 60, 0.87 at 100,
0.72 at 200, 0.60 at 400 and 0.49 at 800, and a neighbour sum costs more per
pair than a full sum, which moves the break-even from 60 up to somewhere
between 110 and 150. It should be measured with the compiled code and changed
if need be: it is `nacmin` in `accheck` in `src/nbac.f90`.

### `nnbopt` (whole number, standard value 0)
The number of neighbours wanted for each star. 0 sets it automatically to
(N/8)^(3/4): 37 for 1000 stars, 210 for 10,000. The size of each star's
neighbour sphere is adjusted as the run goes to keep about this many.

### `nnbmax` (whole number, standard value 0)
The longest neighbour list a star can have. 0 sets it automatically to three
times `nnbopt` plus 16. It sets the memory used (4 bytes times `nnbmax` per
star) and must be at least `nnbopt` + 8.

### `etar` (standard value 0.01)
The accuracy parameter of the timestep on which the slow (regular) force is
recomputed, used in the same formula as `eta`. Smaller is more accurate and
slower. If the energy error with the neighbour scheme is worse than without
it, reduce `etar` first: in tests halving it twice (0.02 to 0.005) reduced the
energy error by about 200.

## `&tree` - Barnes-Hut tree forces
For large numbers of particles. Instead of summing the force from every other
particle, distant groups of particles are replaced by a few terms describing
the group. Energies are still calculated exactly.

### `usetree` (switch, standard value `.false.`)
`.true.` uses the tree for the forces; `.false.` uses direct summation, which
is exact. The tree is worth using above a few thousand particles.

### `theta` (standard value 0.5)
The opening angle, between 0 and 1. A group of particles is treated as a whole
if its size divided by its distance is less than `theta`. Smaller is more
accurate and slower: force errors are about 0.02 per cent at 0.3, 0.2 per cent
at 0.5 and 1 per cent at 0.7.

### `nleaf` (standard value 8)
The largest number of particles in the smallest boxes of the tree. Particles
in the same box as each other are always summed directly. Rarely needs changing.

### `ntreebuild` (whole number, standard value 10)
How often the tree is rebuilt, in block steps (a block step is one pass of the
integrator, in which every star whose timestep ends at that moment is moved).
The tree is rebuilt from scratch every `ntreebuild` block steps. On the steps
in between it is only restocked: the stars stay in the boxes they were put in,
and the mass, centre and size of every box are recomputed from where the stars
are now. Restocking is cheaper than rebuilding.

It does not change the accuracy, because the size of each box is recomputed
and the opening angle is applied to that. It changes only the speed: a tree
that has not been rebuilt for a long time has boxes that have spread out and
overlap, so more of them have to be opened. Too small a value wastes time
rebuilding; too large a value slows the tree walks. It must be at least 1
(rebuild every step). The tree is also rebuilt whenever a binary is closed or
opened.

---

## Standard settings - these should rarely need changing

## `&binaries` - unperturbed binaries and few-body encounters
There is no regularisation. Instead a bound pair that nothing else is
disturbing is replaced by a single body at its centre of mass, and its orbit is
followed exactly as a Kepler orbit. This stops a tight binary forcing tiny
timesteps. As soon as the pair is disturbed it is put back as two stars and
integrated normally. How disturbed a binary is is measured by gamma: the tidal
force of all other particles across the binary's orbit, divided by the
binary's own internal force. Binaries never appear as single bodies in the
output.

### `closebinaries` (switch, standard value `.true.`)
Switches this treatment on. It is switched off automatically, with a warning,
if `eps` is greater than 0.

### `dtclose` (yr, standard value 100.0)
Only stars whose timestep is shorter than this are tested to see if they are
in an unperturbed binary. It simply saves testing stars that cannot be in a
tight binary. Binaries already in the initial conditions are all tested at the
start whatever their timesteps.

### `gclose` (standard value 1.0e-6)
A bound pair is treated as one body when its gamma is below this.

### `gopen` (standard value 1.0e-5)
A binary being treated as one body is split back into two stars when its gamma
rises above this. It must be larger than `gclose`; the gap between them stops
a binary being closed and opened over and over.

### `fewbody` (switch, standard value `.true.`)
`.true.` treats close encounters of a few stars separately. When a binary
meets another star or binary, or three stars come together, and the group of
them (3 to 6 stars) is isolated from the rest of the cluster, the group is
taken out of the main calculation and replaced by a single body at its centre
of mass. Its internal motion is integrated on its own, at the same time as the
cluster, until the encounter is over: when a star is thrown out the group
grows, stops being isolated, and its stars are put back. This is much faster
than leaving them in, because each of the very short steps of a close
encounter otherwise costs a prediction of every star in the cluster. A group
of three that settles into a stable hierarchical triple (an inner binary and
a distant third star, by the Mardling & Aarseth 2001 criterion) is then
treated as two exact Kepler orbits. It needs `closebinaries` to be on.

### `gfew` (standard value 1.0e-5)
A few-body group is taken out if its perturbation gamma by the rest of the
cluster (see `gclose`) is less than `gfew`, and put back if it rises above
10 times `gfew`. Larger values take groups out sooner and keep them longer,
which is faster but ignores more of the cluster's pull on them.

### `etafew` (standard value 0.02)
The timestep accuracy parameter of the few-body integrator: its step is
`etafew` times the shortest free-fall or flyby time of any pair in the group
(about 500 steps per orbit of an eccentric binary at 0.02, with an energy
error of a few parts in 10^9 per orbit). Smaller is more accurate and slower.

## `&collisions` - stellar collisions

### `rcoll` (Rsun, standard value 10.0)
Two stars that pass within this distance of each other have collided. What
happens then is set by `mergers`. 0 switches the check off. The same distance
is used for every type of particle at present. A bound pair whose orbit would
bring it within `rcoll` counts as well, as do stars inside a few-body group.

### `mergers` (switch, standard value `.true.`)
`.true.`: the two stars are merged into one and the run carries on. The new
star has the total mass, sits at the centre of mass with the centre-of-mass
velocity (so mass and momentum are conserved), keeps the id of the more
massive of the two, and takes the type of the more compact (a star and a
black hole make a black hole). No mass is lost. Each merger is

- announced on the screen with a line starting `** MERGER`,
- written as one line to `<runid>_mergers.dat` in the data directory (time,
  ids, masses and types of the two stars and of the new one, its position and
  velocity),
- counted in the header of every snapshot (`mergers` and
  `mergers_since_last_snapshot`).

The number of stars in the snapshots goes down by one for each merger. The
orbital energy of the pair is lost from the stellar dynamics; it is kept track
of, so the energy error on the screen is still a true check of the
integration.

`.false.`: the run stops with an alert giving the two stars, as it did before
mergers were added.

## `&integrator` - the Hermite integrator and its timesteps

### `eta` (standard value 0.02)
The accuracy parameter of the Aarseth timestep criterion, which sets each
particle's timestep from its acceleration and its derivatives. Smaller means
shorter steps: more accurate and slower. 0.01 to 0.02 is usual.

### `etap` (standard value 0.2)
A safety net. No particle's timestep may be longer than `etap` times the
shortest time it would take to fall on to, or fly past, any other particle. It
catches fast close encounters that the Aarseth criterion sees too late.

### `etas` (standard value 0.01)
Sets the very first timestep of each particle (as `etas` times acceleration
divided by its rate of change), before there is enough information for the
Aarseth criterion. Also used when a binary is closed or opened.

### `niter` (standard value 1)
The number of times the corrector is applied each step. 1 is the standard
Hermite scheme. 2 gives roughly ten times better energy conservation for twice
the cost.

### `kmax` (standard value 40, between 1 and 60)
Timesteps are `tout` divided by a power of two; `kmax` is the largest power
allowed, so the shortest possible timestep is `tout`/2^`kmax`. A particle that
wants a shorter step is held at this one, and the code reports how many times
that happened. With the standard value and `tout` = 0.1 Myr the shortest step
is about 9e-8 yr, so the limit should never be reached in practice.

## `&computing` - performance only (these do not change the results)
Set the number of cores before running with `export OMP_NUM_THREADS=8`.
The number in use is printed at the start of the run. All the expensive parts
of the code use every core, except building the tree, which is shared out
between at most 8 cores (one for each of the eight boxes the whole system is
first cut into) - 2, 4 or 8 cores share that work evenly.

### `nompmin` (standard value 32)
When at least this many particles are being moved in one step, the cores share
out the particles between them.

### `nompj` (standard value 1024)
When fewer than `nompmin` particles are being moved, the cores instead share
out the sum over all the other particles - but only if the total number of
particles is at least `nompj`, as for small N sharing costs more than it saves.
Also used for the other loops over all particles.

---

# Initial conditions: `ics.nml`

Run with `./makeics ics.nml` from `HENDOS/run`. With no file name the program
looks for `ics.nml`. The file works in the same way as `params.nml`: sections
that start with `&name` and end with `/`, `name = value` lines, comments after
`!`, and standard values for anything left out. `&ics` must be there; `&fractal`
and `&masses` can be left out.

What the program does with them, in order:
1. makes the structure (Plummer sphere or fractal) with every mass set to 1;
2. gives each star its mass;
3. scales the positions to `radius`;
4. moves to the centre of mass;
5. scales the velocities to the virial ratio `qvir` - or, for a Larson
   fractal only, sets them to the size-linewidth relation;
6. makes any changes of your own (the routine `userfinal`);
7. prints a summary and writes the file.

All stars are single and have type `star`.

These parameters cover the standard choices. For anything else the code is
laid out to be changed: see `src/ics/ICS_README.md`.

---

## `&ics` - what to make (must be present)

### `outfile` (text, no standard value - must be set)
The file to write, with its path, e.g. `'data/plum1k.dat'`. No spaces. The
directory must already exist. This is the name to give as `icfile` in
`params.nml`. An existing file of the same name is replaced.

### `n` (whole number, standard value 1000)
The number of stars. The total mass is not an input: it is whatever `n` stars
drawn from the mass function add up to, and is printed at the end.

### `nensemble` (whole number, standard value 1)
The number of sets of initial conditions to make. With more than 1 `makeics`
makes an **ensemble**: sets that are the same in every way except the random
number seed, which is `seed` for the first, `seed` + 1 for the second, and so
on. Each is written to its own file, with the member number put before the
extension of `outfile`: `outfile = 'data/plum1k.dat'` with `nensemble = 3`
writes `data/plum1k_e001.dat`, `data/plum1k_e002.dat` and
`data/plum1k_e003.dat`. Member k is exactly what a single run of `makeics`
with seed `seed` + k - 1 would make, so one member can be remade on its own.
The header of each file records the ensemble, the member number and the seed
used. To run them all, give `hermite` the same `nensemble` (in `&run` of
`params.nml`) and the `outfile` name as its `icfile`.

Remember that members differ in everything the random numbers decide: the
positions and velocities, and also the masses, so the total mass differs a
little from member to member.

### `seed` (whole number, standard value 12345)
The random number seed. The same seed, with the same other parameters, gives
exactly the same stars again, on any machine. Change it to get a different
realisation of the same system.

### `distribution` (text, standard value `'plummer'`)
- `'plummer'`: a Plummer sphere (Aarseth, Henon & Wielen 1974), cut off at 10
  Plummer radii, with velocities drawn from its distribution function.
- `'fractal'`: a box fractal (Goodwin & Whitworth 2004), set up by the
  `&fractal` section.
- `'user'`: your own structure, from the routine `userstructure` in
  `src/ics/usersubs.f90`. As supplied that makes a uniform sphere.

### `radius` (pc, standard value 1.0)
The size of the system. What it means depends on `distribution`:
- `plummer`: the **half-mass radius**. (The Plummer radius is `radius`/1.288
  for the sphere cut at 10 Plummer radii.) This is the half-mass radius of the
  model; the value for the stars actually made is printed, and differs a
  little because the masses are random.
- `fractal`: the **total radius** - the radius of the sphere the fractal is
  cut from. After moving to the centre of mass a few stars can be slightly
  further out than this.
- `user`: whatever length is 1 in the positions `userstructure` returns (the
  total radius, for the uniform sphere supplied).

### `qvir` (standard value 0.5)
The virial ratio, kinetic energy divided by the size of the potential energy.
0.5 is virial equilibrium, less is cold (the system will collapse), more is hot
(it will expand), 0 is at rest. Once the stars have their masses and the
system has its size, all the velocities are multiplied by one factor to give
exactly this value. **Not used when `velocities = 'larson'`** (see below), or
when `setvirial` is off.

### `setvirial` (switch, standard value `.true.`)
`.true.` scales the velocities to `qvir` as just described. `.false.` leaves
the velocities exactly as the structure routine made them, taking them to be
in km/s already. This is for your own structures (`distribution = 'user'`)
that set their own velocities: for a Plummer sphere or a standard fractal the
velocities made are in arbitrary units and must be scaled. The virial ratio
that results is printed.

---

## `&fractal` - the fractal (only used when `distribution = 'fractal'`)
How the fractal is built is explained in `src/ics/FRACTAL_README.md`.

### `fdim` (standard value 1.6)
The fractal dimension, between 0 and 3. 3.0 is a uniform sphere, 1.6 is very
clumpy. Each box keeps the nearest whole number to `ndiv`^`fdim` of its
children, so only certain dimensions are made exactly. With `ndiv` = 2 these
are 1.58, 2.0, 2.32, 2.58, 2.81 and 3.0; the program prints the one you get.

### `velocities` (text, standard value `'standard'`)
How the velocities of nearby stars are related.
- `'standard'`: the velocities of the original `makefractal`. Children inherit
  their parent's velocity plus a random part that is the same size at every
  level, so nearby stars have similar velocities but the velocity difference
  grows only slowly with separation (power about 0.1).
- `'inherit'`: as `standard`, but the random part is proportional to the size
  of the box, so velocity differences grow with separation (power about 0.8).
- `'random'`: an independent Gaussian velocity for every star, so no velocity
  structure at all.
- `'larson'`: correlated velocities that follow Larson's size-linewidth
  relation: the velocity difference between two stars grows as the square
  root of their separation (see `larson.tex`). This option is different from
  the others in two ways:
  - **the velocities are set by `vlarson`, in km/s, and `qvir` is not used.**
    The virial ratio is whatever results, and is printed. It rises as
    `radius` squared over the total mass: for 1000 stars from the standard
    mass function (about 370 Msun) and `fdim` = 1.6 it is about 0.1 for
    `radius` = 0.5 pc, 0.4 for 1 pc and 1.5 for 2 pc, and it varies by a third
    or more from one seed to the next.
  - it needs an `n` by `n` matrix, so memory is 8 `n`^2 bytes (0.8 GB for
    10,000 stars) and the time grows as `n`^3: fine for a few thousand stars.

### `vlarson` (km/s, standard value 1.0)
Only used with `velocities = 'larson'`. It sets the size of the velocities:
two stars `d` pc apart differ, in each velocity component, by `vlarson` x
sqrt(`d`) km/s (rms). A region filled evenly with stars to a radius R pc then
has a one-dimensional velocity dispersion of 0.72 x `vlarson` x sqrt(R) km/s,
so the standard value 1.0 is the observed relation of Solomon et al. (1987),
0.72 (R/pc)^0.5 km/s. Larson's (1981) own relation, 1.1 (L/pc)^0.38 km/s, has
a different slope and cannot be matched by changing `vlarson`.

### `ndiv` (whole number, standard value 2)
How many pieces each box is cut into along each side at every level, so a box
has `ndiv`^3 children. The standard value 2 gives 8 children. It must be
between 2 and 10. A larger value lets other fractal dimensions be made
exactly but gives fewer, coarser levels.

### `coresize` (pc, standard value 0.0)
The smallest length scale in the fractal (it is not a cluster core radius).
- 0: boxes are divided until there are more than 3.5 times `n` of them, as in
  the original `makefractal`. The smallest scale then depends on `n`. This is
  the usual choice.
- greater than 0: boxes are divided until the children of the next level
  would be closer together than `coresize`, and the positions are blurred by
  `coresize`. It must be less than `radius`
  and small enough to leave at least `n` stars inside the sphere; the number
  made rises quickly as it shrinks. As a guide, for 1000 stars with `fdim` =
  1.6 and `radius` = 1 pc, 0.003 pc works.

If the program stops with `nsmax<nstar`, too few stars were left inside the
sphere: change `seed` (with `coresize` = 0 this happens for roughly one seed
in ten at `fdim` = 1.6) or use a smaller `coresize`.

---

## `&masses` - stellar masses (standard settings)

### `imf` (text, standard value `'maschberger'`)
- `'maschberger'`: each star's mass is drawn from the initial mass function
  of Maschberger (2013, MNRAS, 429, 1725): a power law of slope `alpha` at
  high masses that turns over smoothly around `mu`, between `mlo` and `mup`.
- `'equal'`: every star has the mass `mequal`.
- `'user'`: your own masses, from the routine `usermasses` in
  `src/ics/usersubs.f90`. As supplied that makes a Salpeter mass function
  from 0.5 to 50 Msun.

### `mlo` and `mup` (Msun, standard values 0.1 and 2.0)
The lowest and highest mass a star can have. The standard values are those of
the Sheffield sampling code and give a mean mass of 0.38 Msun and a median of
0.24 Msun. (Maschberger's own limits, 0.01 and 150 Msun, give a mean of 0.36
and a median of 0.10 Msun.)

### `alpha` (standard value 2.3)
The slope of the mass function at high masses (Salpeter is 2.35). Must be
greater than 1.

### `beta` (standard value 1.4)
Sets the shape of the mass function at low masses. Must be greater than 1.

### `mu` (Msun, standard value 0.2)
The mass scale at which the mass function turns over.

### `mequal` (Msun, standard value 1.0)
The mass of every star when `imf = 'equal'`. Not used otherwise.
