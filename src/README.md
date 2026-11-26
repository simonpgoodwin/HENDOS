# HENDOS source code - what is where

The N-body code `hermite` is written in C, in the files below. The initial
conditions generator `makeics` is written in Fortran, in `ics/`. The
original Fortran N-body code is kept in `../archive/fortran` (see
`../archive/fortran/README.md`).

**Keep this file up to date: when a function is added, moved, renamed or changes
what it does, change its entry here in the same edit.**

## Files

| File | Job |
|---|---|
| `hendos.h` | everything shared: constants, parameters, particle arrays, groups, neighbour-scheme arrays, prototypes |
| `main.c` | the main program: set-up, the ensemble launcher, the integration loop |
| `params.c` | units, the parameter-file reader, parameter checks, small helpers |
| `io.c` | input and output: initial conditions, diagnostics, screen, snapshots, headers, in-run analysis |
| `force.c` | gravity: forces, jerks, tide, energies; timestep levels |
| `ac.c` | the Ahmad-Cohen neighbour scheme |
| `tree.c` | the Barnes-Hut tree |
| `few.c` | orbital elements and the triple stability test (stand-alone) |
| `groups.c` | binaries, regularised groups, frozen triples, collisions and mergers: when groups are made, changed and broken up |
| `ar.c` | the regularised integration of groups (algorithmic regularisation), their perturbers and their corrections to the main integration |
| `restart.c` | binary restart files |
| `ics/makeics.f90` and the five files with it | the initial conditions generator (Fortran) |

They are built by `run/Makefile`. A new C file has to be added to `CFILES`
there.

## Conventions

- **Arrays start at 0.** Body `i` runs from 0 to `n`−1, and a group index from 0
  to `nbin`−1. "No group" or "no body" is −1.
- **Vectors are `double[3]`** (`typedef double vec3[3]`), so `r[i][0]` is x of
  body `i`. The neighbour list of body `i` is `nbl[i*nnbmax + k]`,
  k = 0 to `nnb[i]`−1.
- **Every shared variable is declared once, in `hendos.h`.** The `EXT` macro
  makes each declaration a definition in `main.c` (which has
  `#define HENDOS_MAIN`) and `extern` everywhere else.
- **Units.** Internally G = 1, mass in Msun and length in pc, so the velocity
  unit is 0.06558 km/s and the time unit 14.91 Myr. All input and output is in
  Msun, pc, km/s and Myr. The conversions (`vunit`, `tunit`) are only used
  when reading and writing.
- **Block times** are 8-byte integers (`int64_t`) in units of
  `tick = dtmax/2^kmax`, so they are exact. The integer clock goes back to zero
  at every output time, when everyone is synchronised.
- Places you are likely to want to change are marked `>>> CHANGE HERE`.

## hendos.h - shared data
No calculations, only declarations of everything shared between the files:
- the particle type table (`NTYPES`, `typename_`, codes `TSTAR`, `TWDWARF`,
  `TNSTAR`, `TBHOLE`, `TPLANET`). To add a type, increase `NTYPES` and add its
  name to `typename_` in `params.c`. Nothing else changes.
- the parameters read from the parameter file (all described in
  `run/PARAMS_README.md`), derived values, counters and the block clock
- the particle arrays: positions, velocities, force derivatives, masses, ids,
  types (`itype`), group index (`ib`), timestep levels
- the table of groups (Kepler binaries, regularised groups, frozen triples):
  `bnum`, `bmode`, `bmass`, `bident`, `btype`, `brel`, `brel2`, `bx`, `bv`,
  `bepoch`, `bein`, `bcm` (the slot of its c.o.m. body), and for regularised
  groups `bend`, `bend8` (how far it has been integrated) and `npl`, `plist`
  (its perturbers), with up to `KFMAX` = 6 stars and `NPMAX` = 256 perturbers
  each
- the neighbour-scheme arrays (split forces, neighbour lists, regular steps)
- the prototype of every function used outside its own file, and the
  helpers `dot3`, `norm3`

## main.c - the main program
- `setmember(k)` - makes this run member `k` of an ensemble. Its runid and its
  initial conditions file get the member number.
- `runensemble(exe)` - runs every member of an ensemble, `nparallel` at a
  time. Each is a separate run of `hermite`, with its screen output in its
  own log file. Members that have been run before are left alone, so an
  ensemble grows when `nensemble` is raised.
- `main` - reads the parameter file (`setunits`, `readparams`,
  `checkdatadir`). For an ensemble it either calls `runensemble` and stops,
  or, given a member number on the command line, calls `setmember`. Then it
  calls `readic` (or `readrestart`), `accheck`, `firststep`, `collcheck`,
  `binstart`, `startinfo` and writes the first snapshot and analysis. Then
  comes the loop over output intervals, and inside it the loop over blocks:
  1. find the next block time and the active bodies, and integrate ahead the
     regularised groups that need it (`aradvance`)
  2. predict every body to that time (with each body's 2nd and 3rd force
     derivatives from its last step)
  3. `forces`, then correct the active bodies (`niter` times). A body on a
     regular step of the neighbour scheme also gets its regular correction.
  4. accept the step (`acaccept`) and choose new levels (`newlevel`,
     `aclevel`)
  (with the tree, between 2 and 3: rebuild it every `ntreebuild` blocks, at
  the first block of every output interval and after any group change, and
  otherwise restock it with `treemoments`)
  5. `collcheck` and `bincheck`

  At the end of each output interval it calls `arsync`, `statusout`,
  `snapshot`, `ana_do` and `writerestart`, and at the end `ana_end` and `endinfo`.
  With the environment variable `HENDOS_TRACE="t1 t2"` it prints every
  block between t1 and t2 Myr in which the total energy changed by more than
  1e-7 (debugging only, because it costs an O(N²) sum per block).

## params.c - set-up and the parameter file
- `die(msg)`, `xalloc(n)` - stop with a message; `malloc` that stops if it
  fails.
- `setunits` - the conversions between internal units and km/s, Myr.
- `ensname(name,k,out)` - the name of member `k` of an ensemble: `_e007` goes
  before the extension, or at the end if there is none. `makeics` uses the
  same rule, so keep the two the same.
- the parameter table `params[]` - one line per parameter: section, name,
  type, where it is stored. **>>> A new parameter is one line here, plus its
  standard value and checks in `readparams`.**
- `parsefile(pfile)` - reads the parameter file. It understands what a
  Fortran namelist file holds: `&section`, `name = value` separated by commas
  or new lines, `!` comments, `/` at the end of a section, `.true./.false.`
  (or `T/F`), `1.0d-3`, text in ' or ", and `tij(i,j)`. A name that is not a
  parameter of its section stops the run with an error.
  - `lower`, `getreal`, `getlogical`, `setparam` (internal) - its pieces.
- `safename(s)` - true if a name has only letters, digits, - _ and .
- `readparams(pfile)` - sets the standard value of every parameter, reads the
  file, checks the values and converts the tide and a few others to internal
  units. With the tree on, the neighbour scheme is switched off.

## io.c - input and output
Set-up
- `checkdatadir` - makes the data directory if needed, and stops if files
  cannot be written there.
- `allocarrays` - allocates the particle and group arrays.
- `readic` - reads the initial conditions: NBODY-SNAPSHOT 2 with types,
  version 1 without, or the legacy `m x y z vx vy vz`. It sets ids and types,
  converts velocities to internal units and moves everything to the
  centre-of-mass frame.
- `typecode(name)` - the type code of a type name (any case), 0 if unknown.

Diagnostics
- `diagnose(time)` - energies (including the internal energy of every group,
  `emerged` for mergers, and the bookkeeping for frozen triples), centre of
  mass, velocity dispersion and the 10/50/90% Lagrangian radii.
- `lagradii(nf,frac,rl)` - Lagrangian radii for any list of mass fractions.
- `indexsort`, `cmpidx` - the index that sorts an array.
- `typecounts(out)` - text giving the number of each particle type.

Screen
- `startclock`, `walltime`, `restartcounters`, `hms` - clocks and the counters
  each output interval is measured from.
- `startinfo(nout)` - the header at the start of a run.
- `statusout(...)` - one line (standard) or a block (verbose) at every output
  time, plus the warnings.
- `endinfo(time,nforce)` - the summary at the end.

Files
- `hdrs`, `hdrd`, `hdri`, `hdrv` - write one `# key = value` header line.
- `hdrcommon(f)` - the header lines common to every output file.
- `hdrphysics(f)` - the header lines that depend on which physics is on,
  each block written only when its physics is on. **>>> New physics adds
  its word to `physics` and its block here.**
- `snapname(time,fname)` - the snapshot file name, with the data directory.
- `snapshot(time,nforce)` - writes a snapshot: all stars, groups split into
  their stars, in id order, with types.
- `ananame`, `ana_init`, `ana_do`, `ana_end` - the in-run analysis files.
  **>>> How to add an analysis is described above `ananame`.**
- `checkrunid(nsnap,nana)` - stops a new run if output files with its runid
  already exist.

## force.c - gravity and timesteps
- `tideadd(i,a,jk)` - adds the Galactic tide and the Coriolis force to an
  acceleration and jerk.
- `sumrange(...)` (internal, inline) - the pair kernel: acceleration and jerk
  on one body from a range of others, with the shortest pairwise
  free-fall/flyby time and the nearest distance. `nbsum` in `ac.c` is a copy
  of it for neighbour lists, so keep the two identical.
- `fullsum(i,par,...)` - force on `i` from all bodies, shared over the threads
  if `par`.
- `direct(i,par)` (internal) - force on one active body, with the tide, into
  `f0p`, `f1p`.
- `forces` - forces on all the active bodies: the neighbour scheme if
  `useac`, the tree if `usetree` (walks shared over the active bodies),
  otherwise parallel over the active bodies (big blocks) or over the j-sum
  (small blocks).
- `energy(ke,pot,ptid)` - kinetic, potential and tidal energy (direct sum),
  and the potential `phi` of every body.
- `getlevel(dtw)` - the timestep level of a wanted step, limited by `kmax`.
- `commlev(k,t8)` - the nearest level that fits the block structure at the
  current time.
- `firststep` - forces on everyone at the start (building the tree first if
  it is on), and the first levels.
- `newlevel(i)` - the new level of body `i` after its step: the smaller of
  the Aarseth and pairwise criteria (with the tree the Aarseth step is never
  shorter than `etatree` |a|/|j|; why is explained above the function), down by any amount, up one level only
  when the time allows.

## ac.c - Ahmad-Cohen neighbour scheme (only when `useac` is on)
The scheme is explained at the top of the file.
- `accheck(nstars)` - switches the scheme off, and says so, below 150 stars.
  Where 150 comes from is explained above it.
- `acalloc`, `acsetrs0` - sizes of the lists and allocation; first guess of
  the neighbour-sphere radius.
- `nbsum(...)` (internal) - force from the stars in a list.
- `aclist(i,par,init,newl)` (internal) - a new neighbour list, with the
  sphere radius moved towards `nnbopt` neighbours.
- `acforce(i,inpar)` - the force on an active body: an irregular step, a
  regular step (full sum, regular corrector, new list) or the first full sum
  of a new body.
- `acstart(i)`, `acaccept(i)`, `aclevel(i)` - start the regular step of a new
  body; keep the new split forces when a step is accepted; choose the
  regular level after `newlevel`.
- `acmoverow(src,dst)`, `acmoved(src,dst)`, `acremoved(j,irep)`,
  `acadded(jn,i0)` - keep every list right when a body is moved, removed or
  created.

## tree.c - Barnes-Hut tree (only when `usetree` is on)
How the tree works, the force formulae and how it is shared over the threads
are explained at the top of the file.
- `treealloc` - allocates the node arrays.
- `treebuild` - builds the octree for the current predicted positions, then
  calls `treemoments`. The root's eight octants are built in parallel, one per
  thread (2, 4 or 8 threads share them evenly; no more than 8 are used): each
  thread looks through all the bodies for those in its octant and builds that
  subtree in its own region of the node arrays; the regions are then closed
  up, so the tree is identical to a serial build.
- `buildnode(...)` (internal) - recursive: sorts a node's bodies into octants
  and makes its children.
- `octant`, `childcentre` (internal) - which octant a position is in; the
  centre of a child box.
- `treemoments` - mass, centre of mass, quadrupole and octupole moments, the
  time derivative of the quadrupole, velocity and size of every node. Called
  every block between rebuilds (restocking).
- `nodemoments(node)`, `addpoint`, `shiftadd` (internal) - an internal node's
  moments from its children; adding a body; the parallel-axis shift.
- `treeforce(i,a,jk,tmin2,rm2)` - acceleration, jerk, shortest pairwise time
  and nearest distance for body `i` by walking the tree. Opened leaves are
  summed exactly as in `sumrange` in `force.c` (keep the two the same).
- `octc`, `matvec` (internal) - the octupole contraction; a 3x3 product.
- `treeperturb(...)` - the perturbation sum on a group, used by `gperturb`.

## few.c - few-body tools (uses nothing else)
- `elements(mu,rr,vr,a,ecc,h)` - semi-major axis, eccentricity and angular
  momentum of a two-body orbit.
- `hierarchy(...)` - is a triple a stable hierarchy (Mardling & Aarseth
  2001)? It returns the inner pair, the outer star and the two orbits.
(The old Hermite few-body integrator `fewint` has gone; it is in the archived
Fortran.)

## ar.c - regularised groups
How the integration works, how the perturbers are found, predicted and
corrected for, when a group is integrated and how its steps are limited is
explained at the top of the file.

The integrator (internal, one group at a time, thread safe)
- `accint`, `accext` - the stars' accelerations from each other (and U), and
  the tidal accelerations of the perturbers.
- `drift`, `kick`, `leapfrog` - the logarithmic-Hamiltonian leapfrog.
- `gbs` - one Gragg-Bulirsch-Stoer step (2, 4, ..., 16 kicks, extrapolated).
- `hcap` - the longest step that keeps the interpolation good (`ARCAP`).
- `mkplist` - the perturber list, the Taylor series of the perturbers and of
  the c.o.m., and the horizon (end of the c.o.m.'s or a strong perturber's
  step). Weak perturbers whose steps end sooner are held still.
- `arsetup`, `arintegrate` - set a group up and integrate it to its horizon,
  landing on it exactly (secant iteration in s), merging stars that collide
  (`armerge`) and keeping the dense output (`addnode`).

Used by the rest of the code
- `aralloc`, `arreset`, `armovegroup`, `premap`, `ardirty` - bookkeeping:
  memory, a group leaving the table or moving in it, a body changing slot.
- `arstart(k,t8)` - a group has just become regularised.
- `aradvance(tb8)` - at the start of each block, integrates the groups that
  are behind (in parallel).
- `arstate(k,t,...)` - a group's stars at any time (interpolated).
- `arshift(k)` - makes the end state the start state, removing merged ghosts.
- `arsync(t)` - at the end of an output interval.
- `arcorrect()` - at the end of every force calculation, the corrections to
  the forces and jerks of the active bodies.
- `arcapnew(k,t8)` (internal) - shortens the steps of bodies near a group
  that has just been made.
- `arpotcorr(t,pred)` - the energy of the corrections, for the energy check.

## groups.c - binaries, regularised groups, collisions and mergers
How groups are closed, opened and frozen is explained at the top of the file.

Kepler orbits and groups
- `grpeint(nm,mm,xm,um)` - internal energy of a few stars.
- `kepler(mu,r0,v0,dt,r1,v1)` - advances a bound two-body orbit by `dt`
  exactly. Kepler's equation is solved by safeguarded Newton iteration in a
  bracket of ±2 about the mean-anomaly change, which always holds the root.
- `grpcomp(k,xc,vcm,tnow,nm,xm,um)` - positions and velocities now of the
  stars of any group (binary, regularised group or frozen triple).
- `toar(i,t8,tnow)`, `tokepler(i,tnow)` - a Kepler binary or frozen triple
  becomes a regularised group; a regularised pair becomes a Kepler binary.
- `pairorbit`, `pericentre` - binding energy and apocentre of a pair; its
  two-body pericentre.

Collisions and mergers
- `collision(i,j,...)` - the collision alert, and stop (when `mergers` is off).
- `groupmergetype(t1,t2)` - the type of a merger product. **>>> CHANGE HERE**
  for other rules.
- `groupmergelog(...)` - the `** MERGER` alert and a line in
  `<runid>_mergers.dat` (both are used by `ar.c` too).
- `mergestars(i,j,t8,tnow)` - merges two single stars.
- `collcheck(t8,tnow)` - after each block, tests bodies that came close to
  another and merges them (or stops).

Neighbours and perturbations
- `gperturb(...)` - the perturbation gamma on a (possible) group (from the
  tree when it is on and up to date).
- `grpgamma(i,tnow,gcur,gorb)` - gamma of the group whose c.o.m. body is `i`,
  with its current size and with its orbit's size (apocentre of a pair).
- `nearestx(i,jex)` - the nearest body to `i`, not counting `jex`.

Closing and opening
- `moverow(src,dst)` - copies all the data of one body to another slot.
- `startbody(nl,list,t8,levmin)` - forces and first steps of new bodies
  (rebuilding the tree first if it is on).
- `closegroup(nm,list,mode,t8,tnow)` - replaces stars by one c.o.m. body.
- `opengroup(i,t8,tnow,nm,list)` - splits a c.o.m. body back into its stars,
  booking the energy change of a frozen triple in `efrozen`.
- `pairstate(i,j)` - 1 if two single stars can be a Kepler binary, 2 if
  they can be a regularised pair (bound, or unbound and coming together),
  else 0. A pair whose pericentre is inside `rcoll` is never a Kepler binary.
- `tryfreeze(k)` - freezes a regularised group of three if it is a stable
  triple, and records its energy then (`bein`).
- `tryfew(nb,bod,t8,tnow)` - makes a regularised group of 3 or more stars if
  it would be isolated and would pay. When it pays is explained above it;
  `NFEWFAC` is at the top of the file.
- `bincheck(t8,tnow)` - after each block: changes, freezes, opens and closes
  groups among the bodies that have just stepped.
- `binstart(tnow)` - closes binaries and makes regularised groups in the
  initial conditions.

Debugging (only with the environment variable `HENDOS_DEBUG` set)
- `dbgenergy`, `dbgbefore`, `dbgafter` - print the change in total energy
  across every close, open, freeze and change of kind of a group. Each change should be
  about 1e-7 or smaller.

## restart.c - restart files
The file format (a list of named sections that a reader can skip) is
explained at the top of the file.
- `restartname`, `secbegin`, `secend`, `rd` (internal) - file name, sections,
  reading.
- `writerestart(time,iout,nblock,nforce)` - writes the complete state. It
  writes to a temporary name first, so a crash never destroys the previous
  file.
- `readrestart(time,iout,nblock,nforce)` - reads it and sets the run up to
  carry on. Sections it does not know are skipped.

## ics/ - the initial conditions generator (Fortran)
One stand-alone program, `makeics`, separate from the N-body code. How to
change it is in `ics/ICS_README.md`, and the fractals are described in
`ics/FRACTAL_README.md`. Its routines:

### ics/makeics.f90 - the program
Written to be changed: the places to change are marked `! >>> CHANGE HERE`
and `! >>> NEW PARAMETER`. The main program has seven steps, one call each:
`makestructure`, `makemasses`, scale the positions to pc, `tocom`,
`virialise` (or `larsonise` for a Larson fractal, or neither if `setvirial`
is off), `userfinal`, then `report` and `writeics`. The seven steps are
inside a loop over the members of an ensemble. Its internal routines:
- `makestructure` - positions and velocities in model units
  (`makeplummer`, `makefractal`, `larsonfractal` or `userstructure`), with
  every mass set to unity, and `rscale`.
- `makemasses` - masses (`maschberger`, all equal, or `usermasses`) and types.
- `ensname(name,k)` - the file name of member `k` of an ensemble, by the same
  rule as `ensname` in `params.c`.
- `tocom` - to the centre of mass and centre of velocity.
- `virialise` - scales the velocities to the virial ratio asked for.
- `larsonise` - Larson fractals: scales the velocities to the size-linewidth
  relation instead.
- `energies`, `report`, `writeics` - energies, the summary on the screen,
  the NBODY-SNAPSHOT 2 file.
- `readinput`, `nmlerror`, `bad` - reads and checks `&ics`, `&masses` and
  `&fractal`.

### ics/usersubs.f90 - the user's own code
- `userstructure` - your own positions and velocities (`distribution =
  'user'`). As supplied it makes a uniform sphere.
- `usermasses` - your own masses (`imf = 'user'`). As supplied it gives a
  Salpeter mass function.
- `userfinal` - any last changes to anything; always called. As supplied it
  does nothing, and has three examples commented out.

### ics/plummersubs.f90, masssubs.f90, icrandom.f90, fractalsubs.f90
- `makeplummer` - a Plummer sphere cut at 10 Plummer radii, in model units.
- `maschberger` - masses from the L3 IMF of Maschberger (2013).
- `ran3`, `gasdev` - uniform (Knuth's subtractive method) and Gaussian random
  numbers.
- `makefractal`, `larsonfractal` - the Sheffield box fractal (Goodwin &
  Whitworth 2004) and its velocity fields, including Larson's law (see
  `ics/FRACTAL_README.md`).
