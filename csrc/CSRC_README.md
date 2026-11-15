# HENDOS in C (`csrc/`)

This directory holds a C version of the main simulation code (`hermite`), so
the Fortran can be compared against it. The initial conditions (`makeics`)
and the analysis tools are not duplicated: they stay in Fortran and Python
and work with both versions.

The C version reads **the same parameter file** (`run/params.nml`, written as a
Fortran namelist) and **the same initial conditions**. It writes **the same
output files**, with the same names, header keys and columns: snapshots, the
in-run analysis file, the mergers file and ensemble logs. That means the
viewer, `hendosanalysis.py`, `hendosensemble.py`, the regression tests and
the kira benchmark all work on it unchanged. The `creator` header line tells
you which code wrote a file (`HENDOS C version (csrc)`).

## Building and running

```
cd HENDOS/run
make hermite_c                     # needs only gcc (or any C11 compiler with OpenMP)
./hermite_c params.nml             # exactly like ./hermite
./hermite_c params.nml 7           # ensemble member 7 on its own
OMP_NUM_THREADS=4 ./hermite_c params.nml
```

`make` on its own still builds the Fortran `hermite` and `makeics`.
`make hermite_c` builds the C version and does not need a Fortran compiler.
The compiler and flags are `CC` and `CFLAGS` near the top of `run/Makefile`.

To run the regression tests on the C version:

```
cd HENDOS/tests
python3 run_tests.py --hermite ../run/hermite_c
```

## What is in it

The C version uses the same algorithms as the Fortran, routine for routine:

- 4th-order Hermite with block timesteps (integer clock), the Aarseth and
  pairwise criteria, a higher-order predictor, and P(EC)^n
- direct summation with OpenMP
- the Galactic tide (`solar` or `tensor`, rotating frame)
- the Ahmad-Cohen neighbour scheme, including the automatic switch-off below
  150 stars
- binaries as exact Kepler orbits, few-body groups (with the rule that a group
  of nm stars needs at least 2·nm(nm−1) bodies in the main integration), and
  frozen stable triples
- collisions and mergers, with the screen alert and `<runid>_mergers.dat`
- restart files
- ensembles: runs every member, or one member if you give its number, can add
  members later, and keeps the `.log`/`.status` files
- the in-run Lagrangian-radius analysis
- the same screen output, standard or verbose

**Not in the C version: the Barnes-Hut tree.** If `usetree = .true.` it stops
and says so. That means the `plummer_tree` regression test fails on the C
version by design.

**Restart files are not shared between the two versions.** Both are called
`<runid>_restart.bin` and are built the same way, as named sections that a
reader can skip, so new physics can add a section. But the C file is a plain
binary stream that starts with `HENDOS-C-RESTART`, while the Fortran file is
Fortran unformatted. A run has to be restarted by the version that wrote the
file. The C version checks for this and refuses a file it did not write.

## Files

| File | What is in it | Fortran counterpart |
|---|---|---|
| `hendos.h` | all shared variables, constants and prototypes | `nbvars.f90` |
| `main.c` | set-up, the ensemble launcher, the main loop | `hermite.f90`, `runensemble`/`setmember` in `nbsetup.f90` |
| `params.c` | units, the namelist reader, parameter checks, `ensname` | `nbsetup.f90` |
| `io.c` | data directory, initial conditions, diagnostics, screen, snapshots, headers, in-run analysis | `nbsetup.f90`, `nbout.f90`, `nbana.f90` |
| `force.c` | forces, tide, energy, timestep levels | `nbforce.f90`, `nbstep.f90` |
| `ac.c` | the Ahmad-Cohen neighbour scheme | `nbac.f90` |
| `few.c` | the few-body integrator, orbital elements, the triple stability test | `nbfew.f90` |
| `groups.c` | binaries, few-body groups, frozen triples, collisions, mergers | `nbbin.f90` |
| `restart.c` | restart files | `nbrestart.f90` |

The function names match the Fortran routine names, so `src/README.md` works
as a guide to both. The comments point out where the two versions differ.
The differences to keep in mind:

- **Arrays start at 0.** Body `i` runs from 0 to n−1, and a group index from 0
  to nbin−1. "No group" or "no body" is −1 in C and 0 in Fortran.
- **Vectors are `double[3]`,** so `r[i][0]` is x of body `i` (`r(1,i)` in
  Fortran). The list of neighbours of body `i` is `nbl[i*nnbmax + k]`.
- **The namelist reader** (`parsefile` in `params.c`) accepts what the
  Fortran reads from these files: `name = value` separated by commas or new
  lines, `!` comments, `.true./.false.` (also `T/F`), `1.0d-3`, strings in ' or
  ", and `tij(i,j)`. It is stricter than Fortran in one way: it stops with an
  error naming any parameter it does not know.
  **>>> CHANGE HERE: a new parameter is added in two places: one line in the
  `params[]` table in `params.c`, and its default and checks in
  `readparams`.**

## Debugging aids (C only)

- `HENDOS_DEBUG=1 ./hermite_c params.nml` prints the change in total energy
  across every binary or few-body event: close, open, freeze, a few-body group
  formed. Each line should be tiny (about 1e-7 or less). A big one points
  straight at the event that is wrong.
- `HENDOS_TRACE="t1 t2" ./hermite_c params.nml` (times in Myr) prints every
  block between t1 and t2 in which the total energy changed by more than 1e-7,
  with the stars that were active. This is how the Kepler-solver bug below was
  found. It costs an O(N²) energy sum per block, so keep the window short.

## What comparing the two found

The C version has been compiled and run here (gcc 13, 2 cores). The Fortran
still has not, because there is no Fortran compiler in this environment.
Running the C version found these problems in the shared algorithms, and
each was fixed in **both** versions:

1. **The Kepler solver could return the wrong orbit.** It looked for the
   eccentric anomaly within ±1 of the mean-anomaly change. The root can be up
   to 2e away (x − ΔM = e[sin(x + E0) − sin E0]), so for an eccentric binary
   at some phases the solver converged to the wrong place. A closed binary
   could then come out of `grpcomp` with the wrong position and velocity,
   which wrecked the energy when it was opened. This is the worst of the four.
   In 200,000 random orbits (e < 0.999) the old solver's worst energy error
   was 6e2, and the new one's is 1e-11. In a 420-star cluster with 120
   primordial binaries the run's energy error went from 1.2e-2 to 5e-5.
2. **Stars on a colliding orbit were merged too early.** A bound pair whose
   pericentre was inside `rcoll` was merged as soon as it was found, however
   far apart the stars were at the time. Now such a pair is simply never
   closed as a binary, and `collcheck` merges the stars when they actually
   meet.
3. **Frozen-triple energy bookkeeping.** Two nested Kepler orbits leave out
   part of the interaction between the outer star and the inner pair, so the
   true energy of a frozen triple wanders. Each frozen triple is now counted
   at its energy when it was frozen. When it is unfrozen, the difference is
   added to `efrozen`. The new header key `energy_change_frozen_triples`
   reports this, and restart files carry it in a `FROZEN` section.
4. **The test harness** wrote numpy numbers as `np.float64(...)` in the
   parameter files. It now writes plain numbers.

## Speed

The two versions have not been timed against each other, because the Fortran
cannot be built here. With the same compiler family (`gfortran`/`gcc -O3
-march=native`) they should be within a few tens of per cent of each other.
The inner loops are the same and both vectorise.
