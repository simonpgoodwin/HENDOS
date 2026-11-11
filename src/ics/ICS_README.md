# Changing how the initial conditions are made

The generator `makeics` is laid out so that you can go in and change it. This
file says where everything is and gives recipes for the usual changes. The
parameters themselves are described in `run/PARAMS_README.md`, and the
fractals in `FRACTAL_README.md`.

**Keep this file up to date: if the steps, the hooks or the markers in the
code change, change this in the same edit.**

## Three levels of control
1. **The parameter file `run/ics.nml`** - no programming. Structure, size,
   virial ratio, fractal dimension and velocities, mass function and its
   limits.
2. **`src/ics/usersubs.f90`** - the file meant for your own code. Three
   routines that `makeics` calls at fixed points. You edit these and nothing
   else.
3. **`src/ics/makeics.f90`** - the program itself. Every place you are likely
   to want to change is marked with `! >>> CHANGE HERE` and says what to do.

After any change to the code, type `make` in `HENDOS/run`.

## What the program does
The main program is seven steps, one `CALL` each, so the order is plain to
see and easy to alter.

| Step | What | Routine | Where |
|---|---|---|---|
| 1 | structure: positions and velocities, every mass 1 | `makestructure` | `makeics.f90` |
| 2 | masses and particle types | `makemasses` | `makeics.f90` |
| 3 | size: positions to pc | (one line) | main program |
| 4 | frame: centre of mass to the origin, at rest | `tocom` | `makeics.f90` |
| 5 | velocities to km/s: virial ratio, or Larson's relation | `virialise`, `larsonise` | `makeics.f90` |
| 6 | your own last changes | `userfinal` | `usersubs.f90` |
| 7 | summary and file | `report`, `writeics` | `makeics.f90` |

The comment at the top of each routine says what it **needs** and what it
**changes**, which is what you need to know to move it or leave it out.

## The stars
Everything is in four arrays, for stars `i` = 1 to `n`:

| Array | Holds |
|---|---|
| `m(i)` | mass |
| `r(1:3,i)` | position x, y, z |
| `v(1:3,i)` | velocity vx, vy, vz |
| `ptype(i)` | type: `star`, `white_dwarf`, `neutron_star`, `black_hole`, `planet` |

Their units change as the steps go by:

| After step | `m` | `r` | `v` |
|---|---|---|---|
| 1 | 1 | model units | arbitrary |
| 2 | Msun | model units | arbitrary |
| 3 | Msun | pc | arbitrary |
| 5 | Msun | pc | km/s |

Model units are whatever the structure routine works in (Plummer radius = 1
for a Plummer sphere, total radius = 1 for a fractal). Step 1 also sets
`rscale`, the number that turns them into pc in step 3.

## The routines in `usersubs.f90`

| Routine | Called | Switched on by | As supplied |
|---|---|---|---|
| `userstructure(nstar,rstar,vstar,idum)` | step 1 | `distribution = 'user'` | uniform sphere, Gaussian velocities |
| `usermasses(nstar,m,rstar,idum)` | step 2 | `imf = 'user'` | Salpeter mass function, 0.5 to 50 Msun |
| `userfinal(nstar,m,rstar,vstar,ptype,idum)` | step 6 | always called | does nothing (three examples commented out) |

- **`userstructure`** returns a position and velocity for every star. Positions
  are in units in which the length `radius` refers to is 1 (you decide what
  that length is: total radius, half-mass radius, scale length); `makeics`
  multiplies by `radius`. Velocities can be in any units, because step 5
  scales them - unless you set `setvirial = .false.`, when they must be km/s
  and are left alone.
- **`usermasses`** returns a mass in Msun for every star. It is given the
  positions in case the mass should depend on where the star is.
- **`userfinal`** is given everything in final units (Msun, pc, km/s) and may
  change anything. Nothing is rescaled after it, so what it leaves is what is
  written; the summary on the screen is measured afterwards and shows what
  you ended up with.

## Recipes

### Use different mass limits, or a different slope
Parameter file only: `mlo`, `mup`, `alpha`, `beta`, `mu` in `&masses`.

### Make all stars the same mass
Parameter file only: `imf = 'equal'` and `mequal`.

### Make my own distribution of stars
Set `distribution = 'user'` and write it in `userstructure`. The uniform
sphere there shows the pattern: a loop over stars setting `rstar(1:3,i)` and
`vstar(1:3,i)`.

### Give my structure its own velocities, in km/s, and leave them alone
Write them in `userstructure` and set `setvirial = .false.` in `&ics`.

### Use my own mass function
Set `imf = 'user'` and write it in `usermasses`. The Salpeter example there
shows how to turn one random number into a mass.

### Mass segregation
Either make the masses depend on position in `usermasses` (the positions are
passed in), or swap masses between stars in `userfinal`: example 2 there puts
the most massive star at the centre.

### Change particle types
For a rule that depends on mass, the marked place at the end of `makemasses`
in `makeics.f90`, or example 1 in `userfinal`. The names must be ones the
N-body code knows (`typename` in `src/nbvars.f90`).

### Add rotation, a bulk motion, or anything else to the finished system
`userfinal`. Example 3 there adds solid-body rotation.

### Add a new named structure (so it can be chosen in `ics.nml`)
This is the way to add something permanent, such as a King model.
1. Write a routine that returns positions and velocities, in its own file
   (copy `plummersubs.f90`).
2. In `makestructure` in `makeics.f90`, add a `CASE ('name')` that calls it
   and sets `rscale` (how to get from its units to pc, given `radius`).
3. Add the name to the error message in `CASE DEFAULT` just below.
4. Add the file to `ICSRC` in `run/Makefile`.

Nothing else needs to know about it: steps 2 to 7 work for any structure.

### Add a new named mass function
The same four steps, in `makemasses` (copy `masssubs.f90`).

### Add a new input parameter
Three places in `makeics.f90`, each marked `! >>> NEW PARAMETER`: declare it,
add it to the namelist of its section, and give it a standard value. Then
put it in `run/ics.nml` and describe it in `run/PARAMS_README.md`. To pass it
to a routine in another file, add it to that routine's arguments.

### Change what `radius` means
`rscale` in `makestructure`. For example, for `radius` to be the Plummer
radius rather than the half-mass radius, set `rscale=radius` in the Plummer
case.

### Leave out a step, or change the order
Edit the main program in `makeics.f90`. To keep the system out of its
centre-of-mass frame, put a `!` in front of `CALL tocom`. To do something
between two steps, add your own `CALL` there. Check the "needs" and "changes"
of the routines you move: the masses must be set before `tocom` and
`virialise`, for example, because both use them.

### Set the velocities in another way
Step 5 in the main program chooses between `larsonise` and `virialise`. Two
common changes:
- *A Larson fractal with a virial ratio of my choosing:* call `virialise`
  after `larsonise`. The velocities keep the Larson slope and are scaled to
  `qvir`, so the coefficient is no longer `vlarson`.
- *Scale to something other than the virial ratio:* in `virialise`, `fac` is
  the one number every velocity is multiplied by; work it out your own way.

### Record something extra in the output file
Add a header line in `writeics`, like the ones there. The N-body code ignores
header lines it does not know.

## Rules to keep
- **Do not change the columns of the output file** (`id m x y z vx vy vz
  type`): the N-body code reads them.
- **Take every random number from `ran3(idum)` or `gasdev(idum)`**, pass `idum`
  on unchanged, and never set it. That is what makes one seed give the same
  stars every time.
- **Check the summary on the screen.** It is measured from the stars that are
  written - total mass, half-mass radius, virial ratio, velocity dispersion -
  so it shows what your changes actually did.
- All stars are single. Binaries and multiple systems are not made yet.
