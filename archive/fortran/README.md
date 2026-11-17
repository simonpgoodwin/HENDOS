# HENDOS - the archived Fortran version of the N-body code

This is the original Fortran version of the HENDOS N-body code. **It is no
longer the version in use.** That is the C code in `HENDOS/src`, built as
`run/hermite`. The Fortran is kept as the reference that the C code was
translated from, routine by routine: `src/README.md` here describes every
Fortran routine, and the C functions have the same names.

The initial conditions generator `makeics` is still Fortran, but it is not
archived. It is in `HENDOS/src/ics` and is built by `HENDOS/run/Makefile`.

## Status: never compiled

None of this Fortran has been compiled or run. It was written without a
Fortran compiler. The C translation of it has been compiled, run and tested,
and every bug that testing found was fixed here too (see below). So the
algorithms are tested, but this Fortran text is not. The first compilation
will probably need small fixes.

## Building and running

```
cd HENDOS/archive/fortran/run
make                               # builds hermite_fortran here (gfortran + OpenMP)
cd ../../../run
../archive/fortran/run/hermite_fortran params.nml
```

It reads the same parameter files and initial conditions as the C code, and
it writes the same snapshots, headers, analysis, merger and ensemble files.
The analysis tools therefore work with either version. The `creator` header
line says which code wrote a file: `HENDOS (hermite.f90)` for this one,
`HENDOS (C, src/main.c)` for the C code.

**Restart files are not shared.** Both versions call the file
`<runid>_restart.bin`, but this one writes Fortran unformatted records, while
the C file is a plain stream that starts with `HENDOS-C-RESTART`. A run must
be restarted by the code that wrote the file.

To run the regression tests on it:
```
cd HENDOS/tests
python3 run_tests.py --hermite ../archive/fortran/run/hermite_fortran
```

## What the C translation found (fixed in both versions)

1. **Kepler solver.** It looked for the eccentric anomaly within ±1 of the
   mean-anomaly change, but the root can be up to 2e away
   (x − ΔM = e[sin(x + E0) − sin E0]). So for an eccentric binary at some
   phases it converged to the wrong place. A closed binary could then come
   back from `grpcomp` with the wrong position and velocity, which wrecked the
   energy when the binary was opened. Over 200,000 random orbits (e < 0.999),
   the worst energy error was 6e2 with the old solver and 1e-11 with the new
   one. In a 420-star cluster with 120 primordial binaries, the run's energy
   error went from 1.2e-2 to 5e-5.
2. **Stars on a colliding orbit were merged too early.** A bound pair whose
   pericentre was inside `rcoll` was merged as soon as it was found, however
   far apart the stars were at the time. Now such a pair is never closed as a
   binary, and `collcheck` merges the stars when they actually meet.
3. **Frozen-triple energy bookkeeping.** Two nested Kepler orbits leave out
   part of the interaction between the outer star and the inner pair, so the
   true energy of a frozen triple wanders. Each frozen triple is now counted
   at its energy when it was frozen. When it is unfrozen, the difference is
   added to `efrozen`. The header key `energy_change_frozen_triples` reports
   it, and restart files carry it in a `FROZEN` section.
4. **The standard `eta` is now 0.01** (it was 0.02). There is no
   regularisation, so perturbed hard binaries are integrated directly. In
   binary-rich clusters `eta` = 0.02 gave energy errors of 1e-2 or worse,
   while 0.01 gave about 1e-4.
5. **The tree and the timesteps** (found when the tree was put into the C
   code). The tree force on a star jumps a little whenever a box switches
   between being used whole and being opened. The Hermite corrector reads
   the jumps as large higher derivatives, and the Aarseth criterion then makes
   the steps far too short: a 4,000-star tree run took 9 times as many steps
   as direct summation. The C code fixes this with `etatree` (no step shorter
   than `etatree` |a|/|j|); **this Fortran does not have that fix**. The C
   code also rebuilds the tree at every output time, so restarts are exact,
   and does not use the tree for nearest-neighbour searches.
