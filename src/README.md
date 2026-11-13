# HENDOS source code - what is where

Each file is one Fortran module that does one job. The main program
(`hermite.f90`) only contains the integration loop; everything else is in the
modules. `(internal)` marks a routine contained inside the one above it.

**Keep this file up to date: when a routine is added, moved, renamed or changes
what it does, change its entry here in the same edit.**

## Files, in the order they are compiled

| File | Module | Job | Uses |
|---|---|---|---|
| `nbvars.f90`  | `nbvars`  | shared data: parameters, particle arrays, particle types | - |
| `nbsetup.f90` | `nbsetup` | set-up and input: units, parameter file, data directory, initial conditions | nbvars |
| `nbtree.f90`  | `nbtree`  | Barnes-Hut tree | nbvars |
| `nbac.f90`    | `nbac`    | Ahmad-Cohen neighbour scheme | nbvars |
| `nbforce.f90` | `nbforce` | gravity: forces, jerks, tide, energies | nbvars, nbtree, nbac |
| `nbstep.f90`  | `nbstep`  | timesteps: starting steps and block timestep levels | nbvars, nbtree, nbac, nbforce |
| `nbfew.f90`   | `nbfew`   | few-body integrator and triple stability test (stand-alone) | - |
| `nbbin.f90`   | `nbbin`   | binaries, few-body groups, collisions and mergers | nbvars, nbtree, nbac, nbforce, nbstep, nbfew |
| `nbout.f90`   | `nbout`   | output: screen, diagnostics, snapshots | nbvars, nbforce, nbbin |
| `nbrestart.f90` | `nbrestart` | binary restart files | nbvars, nbsetup, nbac, nbout |
| `nbana.f90`  | `nbana`  | in-run analysis files | nbvars, nbout |
| `hermite.f90` | program `hermite` | main program: the integration loop | all of the above |
| `ics/makeics.f90` | program `makeics` | the initial conditions generator | the five files below |
| `ics/plummersubs.f90` | (no module) | Plummer sphere: `makeplummer` | icrandom |
| `ics/fractalsubs.f90` | (no module) | box fractals: `makefractal`, `larsonfractal` | icrandom |
| `ics/masssubs.f90` | (no module) | stellar masses: `maschberger` | icrandom |
| `ics/usersubs.f90` | (no module) | the user's own code: `userstructure`, `usermasses`, `userfinal` | icrandom |
| `ics/icrandom.f90` | (no module) | random numbers: `ran3`, `gasdev` | - |

The compile order and dependencies are in `run/Makefile`; a new module must be
added to `MODS` there (in the right place) with a line saying what it uses.

## nbvars.f90 - shared data
No calculations: declarations of everything shared between modules.
- **particle type table**: `ntypes`, `typename` (the names used in files) and
  the named codes `tstar`, `twdwarf`, `tnstar`, `tbhole`, `tplanet`. Adding a
  type means changing only this table.
- **parameters** read from the parameter file (all described in
  `run/PARAMS_README.md`).
- **derived and internal values**: unit conversions, tide in internal units,
  counters, block-timestep clock.
- **particle arrays**: positions, velocities, force derivatives, masses, ids,
  types (`itype`), timestep levels; and the table of closed binaries
  (`bmass`, `bident`, `btype`, `brel`, `bepoch`), now a table of *groups*:
  binaries, few-body groups and frozen triples (`bnum`, `bmode`, `bx`, `bv`,
  `brel2`), with up to `kfmax` = 6 stars each; and the arrays of the neighbour
  scheme (split forces, neighbour lists, regular steps).
- `typecode(name)` - function: the integer code of the particle type called
  `name` (any case), or 0 if there is no such type.

## nbsetup.f90 - set-up and input (everything done once at the start)
- `setunits` - sets the conversions between internal units (G = 1, Msun, pc)
  and km/s and Myr.
- `readparams(pfile)` - sets the standard value of every parameter, reads the
  namelist sections of the parameter file, checks the values, and converts the
  tide and a few others to internal units.
  - `nmlerror(name)` (internal) - stops with the compiler's message if a
    section is wrong; notes sections that are missing.
  - `safename(str)` (internal) - true if a name has only letters, digits, - _ .
- `checkdatadir` - makes the data directory if needed and stops if files
  cannot be written there.
- `allocarrays` - allocates all the particle and binary arrays.
- `readic` - reads the initial conditions (NBODY-SNAPSHOT 2 with types,
  version 1 without, or legacy `m x y z vx vy vz`), sets ids and types,
  converts velocities to internal units and moves to the centre-of-mass frame.

## nbtree.f90 - Barnes-Hut tree (only used when `usetree` is on)
- `treealloc` - allocates the node arrays.
- `treebuild` - builds the octree for the current predicted positions, then
  calls `treemoments`. The root's eight octants are built in parallel, one
  per thread (2, 4 or 8 threads share them evenly; no more than 8 are used):
  each thread looks through all the bodies for those in its octant and builds
  that subtree in its own region of the node arrays; the regions are then
  closed up, so the tree is identical to a serial build.
- `buildnode(...)` - recursive: sorts a node's bodies into octants and makes
  its children. Thread safe (its node counter is an argument).
- `octant(x,centre)` - function: which of the 8 octants a position is in.
- `treemoments` - mass, centre of mass, quadrupole and octupole moments,
  velocity and size of every node (leaves from bodies, then parents from
  children). Called every block; cheaper than a rebuild. Parallel over the
  leaves, then over the top-level subtrees for the internal nodes.
- `nodemoments(node)` - the moments of one internal node from its children.
- `addpoint(...)` - adds one body to a leaf's moments.
- `shiftadd(...)` - adds a child's moments to its parent's (parallel-axis shift).
- `treeforce(i,...)` - acceleration, jerk, shortest pairwise timescale and
  nearest-neighbour distance for body `i` by walking the tree.
- `octc(...)` - the octupole contraction used by `treeforce`.
- `treeperturb(...)` - tidal perturbation sum on a binary from the tree.
- `treenearest(i,jn)` - nearest body to body `i` from the tree.

## nbac.f90 - Ahmad-Cohen neighbour scheme (only used when `useac` is on)
The scheme is explained at the top of the file.
- `accheck(nstars)` - switches the scheme off, and says so, if there are
  fewer than 150 stars.
- `acalloc` - chooses the list sizes if they were left automatic and
  allocates the arrays.
- `acsetrs0` - first guess of the neighbour-sphere radius.
- `acforce(i,inpar)` - the force on active star `i`: an irregular step
  (neighbour sum plus the extrapolated regular force), a regular step (full
  sum, regular-force corrector, new neighbour list) or the first full sum of
  a new body.
- `pairsum(...)` - adds the force of one star on another (the same pair
  interaction as in `forces`).
- `nbsum(i,nl,list,...)` - force on `i` from the stars in a list.
- `fullsum(i,par,...)` - force on `i` from all stars, shared over the threads
  if `par`.
- `actide(i,a,jk)` - adds the Galactic tide to a total force and jerk.
- `aclist(i,par,init,nl,newl)` - makes a new neighbour list and adjusts the
  sphere radius towards `nnbopt` neighbours.
- `acmoverow(src,dst)` - copies a body's neighbour-scheme data to another slot.
- `acmoved(src,dst)`, `acremoved(j,irep)`, `acadded(jn,i0)` - keep every
  list correct when a body is moved, removed or created.

## nbforce.f90 - gravity
- `forces` - acceleration and jerk on every active particle from all the
  others at their predicted positions: the tree if `usetree`, otherwise direct
  summation (two OpenMP versions, for big and small blocks). Also returns the
  shortest pairwise free-fall/flyby time and nearest-neighbour distance.
  - `finish(...)` (internal) - adds the Galactic tide and Coriolis terms and
    stores the results for one particle.
- `energy(ke,pot,ptid)` - kinetic, potential and tidal energy of the whole
  system (direct sum), and the potential `phi` of every particle.

## nbstep.f90 - timesteps
- `firststep` - forces on everyone at the start and each particle's first
  timestep level.
- `newlevel(i)` - new timestep level of particle `i` after its step: the
  smaller of the Aarseth criterion and the pairwise criterion, with the block
  rules (down by any amount, up one level only when the time allows).
- `commlev(k,t8)` - function: nearest timestep level that fits the block
  structure at the current time.
- `acstart(i)`, `acaccept(i)`, `aclevel(i)` - neighbour scheme: start the
  regular step of a new body; keep the new split forces when a step is
  accepted; choose the regular timestep level after `newlevel`.
- `getlevel(dtw,k)` - turns a wanted timestep into a level (power of two of
  `tout`), limited by `kmax`; counts steps held at the limit.

## nbfew.f90 - few-body tools (stand-alone: no other module is used)
- `fewacc(...)` - accelerations and jerks of the stars of a group from each
  other, the shortest pairwise timescale and the closest pair.
- `fewint(nm,mm,x,v,tint,...)` - integrates an isolated group for a given
  time: 4th-order Hermite, one shared timestep. Stops early if two stars
  collide.
- `elements(mu,rr,vr,a,ecc,h)` - semi-major axis, eccentricity and angular
  momentum of a two-body orbit.
- `hierarchy(...)` - is a triple a stable hierarchy (Mardling & Aarseth
  2001)? Returns the inner pair, the outer star and the two orbits.

## nbbin.f90 - binaries, few-body groups, collisions and mergers
How groups are closed, opened and frozen is explained at the top of the file.

Kepler orbits and groups
- `kepler(mu,r0,v0,dt,r1,v1)` - advances a bound two-body orbit by `dt` exactly.
- `bincomp(k,...)` - positions and velocities of the two stars of binary `k`.
- `grpcomp(k,xc,vc,tnow,nm,xm,um)` - positions and velocities of the stars
  of any group `k` now (binary, few-body group or frozen triple).
- `fewadvance(k,xc,vc,tnow)` - brings the integration of few-body group `k`
  up to date, merging stars inside it that collide.
- `pairorbit(...)` - binding energy and apocentre of a pair.
- `pericentre(...)` - function: two-body pericentre distance of a pair.

Collisions and mergers
- `collision(i,j,...)` - prints the collision alert and stops the run (used
  when `mergers` is off).
  - `collstar(l)` (internal) - prints the id, mass and type of one body.
- `mergetype(t1,t2)` - function: the type of the star made by a merger.
- `mergelog(...)` - announces a merger on the screen and adds a line to
  `<runid>_mergers.dat`.
- `mergestars(i,j,t8,tnow)` - merges two single stars into one.
- `collcheck(t8,tnow)` - after each block, tests active bodies that came close
  to another for a collision (within `rcoll`) and merges them or stops.
- `paircollide(i,j,t8,tnow)` - a bound pair whose orbit touches: merge or stop.

Finding neighbours and perturbations
- `gperturb(...)` - the perturbation gamma on a group from all other bodies.
- `grpgamma(i,tnow,gam)` - gamma of the existing group whose c.o.m. is `i`.
- `nearest(i,jn)` - nearest body to body `i`.
- `nearestx(i,jex,jn)` - the same, not counting body `jex`.

Closing and opening
- `moverow(src,dst)` - copies all the data of one body to another slot.
- `startbody(nl,list,t8,levmin)` - forces and first timesteps of new bodies.
- `closegroup(nm,list,mode,...)` - replaces `nm` stars by one c.o.m. body.
- `opengroup(i,t8,tnow,nm,list)` - splits a c.o.m. body back into its stars.
- `pairstate(i,j)` - function: 0 nothing, 1 the pair can be closed as a
  binary, 2 the pair will collide.
- `tryfew(nb,bod,t8,tnow,formed)` - makes a few-body group of some bodies if
  it would be isolated.
- `tryfreeze(k)` - freezes a few-body group of three if it is a stable triple.
- `bincheck(t8,tnow)` - after each block: freezes, opens and closes groups
  among the bodies that have just stepped.
- `binstart(tnow)` - closes binaries, merges touching pairs and makes
  few-body groups in the initial conditions.

## nbout.f90 - output
Small tools
- `hms(sec)` - function: seconds as h:mm:ss.
- `walltime()` - function: wall-clock seconds since the integration started.
- `typecounts()` - function: text giving the number of each particle type.
- `indexsort(nn,a,idx)` - index that sorts an array (heapsort).
- `datestamp()` - function: date and time as text.
- `physlist()` - function: the physics switched on, for the header key
  `physics`. **New physics adds its word here.**

Diagnostics
- `diagnose(time)` - energies (including closed binaries), centre of mass,
  velocity dispersion and 10/50/90% Lagrangian radii of the current state.
- `lagradii(nf,frac,rl)` - Lagrangian radii for any list of mass fractions.
- `levelhist(hist)` - number of bodies on each timestep level.

Screen
- `writelevels` - prints the table of timestep levels (verbose).
- `startinfo(nout)` - the header at the start of the run.
- `statusout(...)` - the line (standard) or block (verbose) at every output
  time, and the warnings.
- `endinfo(time,nforce)` - the summary at the end.

Files
- `hdrs`, `hdrd`, `hdrv`, `hdri` - write one `# key = value` header line (text,
  real, several reals, integer).
- `hdrcommon(u)` - the header lines common to every output file.
- `hdrphysics(u)` - the header lines that depend on which physics is on (tide,
  binaries, few-body, mergers, neighbour scheme), each block written only when
  its physics is on. **New physics adds its header block here.**
- `snapname(time)` - function: the snapshot file name, including the data
  directory.
- `snapshot(time,nforce)` - writes a snapshot (all stars, closed binaries
  split into their two stars, in id order, with types).

## nbrestart.f90 - restart files
The file format (a list of named sections) is explained at the top of the file.
- `restartname()` - function: the restart file name.
- `secthead(u,name,nrec)` - starts a section of the file.
- `writerestart(time,iout,nblock,nforce)` - writes the complete state.
- `readrestart(time,iout,nblock,nforce)` - reads it, allocates the arrays and
  sets up the run to carry on. Sections it does not know are skipped.

## nbana.f90 - in-run analysis
How to add a new analysis is described at the top of the file.
- `ananame(name)` - function: the analysis file name, including the data
  directory.
- `ana_open(u,name,columns,units)` - opens an analysis file (appends on a
  restart) and writes its header.
- `ana_init` - opens the files of the analyses that are switched on.
- `ana_do(time)` - does all the analyses that are switched on.
- `ana_end` - closes the files.
- `ana_lagrangian(time)` - the example analysis: 1-90% Lagrangian radii and
  velocity dispersion.
- `checkrunid(nsnap,nana)` - stops a new run if output files with its `runid`
  already exist.

## hermite.f90 - main program
No subroutines. In order: read the parameter file name, `setunits`,
`readparams`, `checkdatadir`; then either `readic`, set up the block clock,
`firststep`, `collcheck`, `binstart`, `startinfo`, first snapshot and analysis,
or (with `restart`) `readrestart`. Then the loop over output intervals; inside
it the loop over blocks:
1. find the next block time and the active particles (the list is made in
   parallel: one piece of the particles per thread, joined up in order)
2. predict all particles to that time
3. rebuild the tree or refresh its moments (if `usetree`)
4. `forces`, then correct the active particles (`niter` times); with the
   neighbour scheme a star on a regular step also gets its regular correction
5. accept the step (`acaccept`) and choose new timestep levels (`newlevel`,
   `aclevel`)
6. `collcheck` and `bincheck`

At the end of each output interval: `statusout`, `snapshot`, `ana_do`,
`writerestart`. Finally `ana_end` and `endinfo`.

## ics/ - the initial conditions generator
One stand-alone program, `makeics`, separate from the N-body code (it uses none
of its modules). The program does the steps common to every set of initial
conditions; the routines it calls each make one thing and know nothing about
the others. All random numbers come from `ran3`.

### ics/makeics.f90 - the program
Written to be changed: the places to change are marked `! >>> CHANGE HERE`
and `! >>> NEW PARAMETER`, and `ics/ICS_README.md` gives recipes. The main
program is seven steps, one call each: `makestructure`, `makemasses`, scale
the positions to pc, `tocom`, `virialise` (or `larsonise` for a Larson
fractal, or neither if `setvirial` is off), `userfinal`, then `report` and
`writeics`. Its internal routines:
- `makestructure` - step 1. Calls `makeplummer`, `makefractal`,
  `larsonfractal` or `userstructure` for the positions and velocities in
  model units, sets every mass to unity, and sets `rscale`, the factor that
  takes the positions to pc (from the half-mass radius for a Plummer sphere,
  the total radius for a fractal). A new structure is a new `CASE` here.
- `makemasses` - step 2. Gives the stars their masses (`maschberger`, all
  equal, or `usermasses`) and their types. A new mass function is a new
  `CASE` here.
- `tocom` - step 4. Moves to the centre of mass and centre of velocity.
- `virialise` - step 5. Scales all velocities by one factor to give exactly
  the virial ratio asked for.
- `larsonise` - step 5, Larson fractals only. Scales the velocities to the
  observed size-linewidth relation in km/s instead.
- `energies` - kinetic and potential energy of the stars as they are.
- `report` - step 7. Summary on the screen, measured from the stars as
  written (masses, half-mass radius, virial ratio, velocity dispersion).
- `writeics` - step 7. Writes the NBODY-SNAPSHOT 2 file.
- `readinput` - sets the standard values, reads the namelist sections `&ics`,
  `&masses` and `&fractal` of the parameter file and checks them.
- `nmlerror(name,ios,msg)` - stops with the compiler's message if a section
  is wrong.
- `bad(text)` - stops with an error message.

### ics/usersubs.f90 - the user's own code
The file to edit to make initial conditions your own way. Each routine is
called by `makeics` at one point.
- `userstructure(nstar,rstar,vstar,idum)` - your own positions and
  velocities; used when `distribution = 'user'`. As supplied: a uniform
  sphere with Gaussian velocities.
- `usermasses(nstar,m,rstar,idum)` - your own masses; used when `imf =
  'user'`. As supplied: a Salpeter mass function from 0.5 to 50 Msun.
- `userfinal(nstar,m,rstar,vstar,ptype,idum)` - any last changes to anything;
  always called, just before the output. As supplied: does nothing, with
  three examples commented out (black holes above a mass, most massive star
  to the centre, rotation).

### ics/plummersubs.f90 - Plummer sphere
- `makeplummer(nstar,rstar,vstar,idum)` - positions and velocities of a
  Plummer sphere cut at 10 Plummer radii, in model units (G = 1, total mass
  1, Plummer radius 1), for equal-mass stars. Double precision.

### ics/masssubs.f90 - stellar masses
- `maschberger(nstar,m,idum,alpha,beta,mu,mlo,mup)` - masses of single stars
  from the L3 IMF of Maschberger (2013), sampled exactly with one random
  number per star (the Sheffield sampling code).

### ics/icrandom.f90 - random numbers
- `ran3(idum)` - function: uniform random number between 0 and 1, by Knuth's
  subtractive method (the method of `ran3` in Numerical Recipes, same calling
  convention, written for HENDOS from the published algorithm). Started by
  calling it with a negative `idum`.
- `gasdev(idum)` - function: Gaussian random number (zero mean, unit
  variance) made from `ran3`.

### ics/fractalsubs.f90 - box-fractal routine
How the fractal is built, the core size and the velocity options are explained
in `ics/FRACTAL_README.md`.
The Sheffield box fractal (Goodwin & Whitworth 2004), from the routines
`makefractal` and `makefractal2`; lines that differ from the routines as
supplied are marked `! HENDOS:`.
- `makefractal(nstar,rstar,vstar,fdim,ndiv,idum,radius,vscale,coresize,vtype)` -
  makes the positions, and every velocity field except Larson's.
  - positions: a box is divided into `ndiv`^3 (8 as standard), the nearest
    whole number to `ndiv`^`fdim` of them are kept and divided again, and a
    sphere of radius `radius` is cut out. With `coresize` = 0 dividing stops
    when there are more than 3.5 times `nstar` boxes (the original
    `makefractal`); with `coresize` > 0 it stops when the boxes are smaller
    than `coresize` (`makefractal2`).
  - velocities, `vtype`: `'f'` inherited from the parent plus a random part
    proportional to box size; `'s'` inherited plus a random part the same
    size at every level (the original `makefractal`); `'u'` uncorrelated;
    `'n'` none (zero).
- `larsonfractal(nstar,rstar,vstar,fdim,ndiv,idum,radius,vscale,coresize)` -
  a fractal whose velocities follow Larson's law (velocity difference
  proportional to the square root of separation). It calls `makefractal` for
  the positions (with `vtype` `'n'`), then draws all the velocities together
  as correlated Gaussian numbers. Needs an `nstar` by `nstar` matrix.
