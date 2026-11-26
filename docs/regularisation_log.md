# HENDOS - regularised binaries and multiples: where things stand

Log written 24 November 2026, when the work was paused. Everything below is
committed and pushed (GitHub `simonpgoodwin/HENDOS`, commit `7474e9a`).

## What was asked
1. Fix the kira benchmark (tighten `etap` with `eta`, use the largest energy
   error over a run, plot every run's r_h, add neighbours to the fractal test,
   several seeds).
2. Make `etap` = 0.1 the standard value.
3. Treat perturbed binaries properly (KS or kira-style), and find accurate and
   fast ways of dealing with higher-order multiples.

Items 1 and 2 are done (commit `c450089`). Item 3 is done and working, with
some loose ends, which are listed at the end.

## 1 and 2: the benchmark and `etap`
- `bench_kira.py` sets `etap` = 5 `eta` and `etar` = 0.5 `eta` in every HENDOS
  run (`ETAP_OVER_ETA`, `ETAR_OVER_ETA`).
- The energy error plotted is now the largest |dE/E| over all the snapshots;
  the error at the end is also given in the table.
- The r_h plot shows every run: the most accurate of each code thick, the
  rest thin.
- `run_benchmarks.sh` has neighbours in the fractal test and uses three seeds
  (`SEEDS="1 2 3"`) by default.
- The standard `etap` is 0.1 (it was 0.2). PARAMS_README has the table that
  shows why: with `etap` held at 0.2 the close encounters set a floor to the
  accuracy, whatever `eta` is.

## 3: what was built
**A new file, `src/ar.c`**: algorithmic regularisation for groups of 2 to 6
stars (`bmode` 2). It uses the logarithmic-Hamiltonian leapfrog (Mikkola &
Tanikawa 1999) with Gragg-Bulirsch-Stoer extrapolation (as in AR-CHAIN). It
is regular at any pericentre, takes about 10-20 steps per orbit, and each step
is accurate to `artol` = 1e-10.

How a group is coupled to the cluster:
- **Perturbers.** These are bodies whose own gamma on the group is above
  `gpert` (1e-7); at most 64 were kept (256 since 26 November, see the addendum).
  - *Strong* perturbers (gamma above 1e-5) are followed exactly.
  - *Weak* ones act through a summed tidal tensor, interpolated in time.
  - Perturbers are predicted from their Taylor series, never past the end of
    their own step.
- **Back-reaction.** The perturbers, and the group's centre of mass, feel the
  group's *stars*, not a point mass (`arcorrect`, called at the end of every
  force calculation). So there is one Hamiltonian, and the energy check adds
  the correction energy (`arpotcorr`).
- **When a group is integrated.** At the start of a block, a group is
  integrated ahead to its *horizon*: the end of the step of its centre of mass
  or of a strong perturber. In between, its stars are found by quintic
  interpolation. It lands exactly on the horizon by a loose secant iteration
  followed by a tiny Taylor step.
- **Timesteps.**
  - The pairwise step limit (`etap`) and the closest approach are measured
    against the group's *stars*.
  - When a group is made, the steps of the bodies that feel it are cut if
    necessary.
- **Collisions** inside groups are merged as before.
- **Restarts** are exact, and the restart file format is unchanged.

**The group logic in `groups.c`** (summarised at the top of the file):
- A pair becomes a **Kepler binary** if it is unperturbed all round its orbit
  (gamma at apocentre < `gclose`).
- It becomes a **regularised pair** if it is isolated *now* (gamma with its
  current separation < `gfew`). That includes unbound pairs coming together
  (close hyperbolic encounters), and very eccentric binaries plunging through
  pericentre. Judging by the current separation is how KS regularisation does
  it.
- Changes of kind happen in place:
  - Kepler binary -> regularised when gamma > `gopen`, and back when it is
    below `gclose`.
  - Frozen triple <-> regularised, on the same thresholds.
  - Triples are frozen only when they are unperturbed.
- When gamma > 10 `gfew`, the group takes its nearest body in, or else is
  opened.
- **Groups meeting groups.** A group whose nearest body is another group
  perturbing it by more than `gfew` on its own takes it in, or, if the
  combined group would not be isolated, is opened. Two groups see each other
  only as point masses, and this was the largest remaining source of error.

**Parameters.**
- `gfew` is now 1e-2 (it was 1e-5).
- New: `artol` (1e-10) and `gpert` (1e-7).
- `etafew` is no longer used; it is accepted, with a note.
- The old Hermite few-body integrator (`fewint`) has been removed from the C
  code. It is still in the archived Fortran.

**Tests.** There are 25 regression tests (5 new); all pass, also under the
address and undefined-behaviour sanitisers. The new ones:
- **A flyby** past a binary, against a very accurate direct run (with no
  binary treatment): the binary's semi-major axis agrees to 2e-10.
- **A regularised triple** over 20 outer orbits: energy error 2e-10.
- **An e = 0.9999 binary** plunging past a third star: 2e-10.
- **A hyperbolic near miss** at 2 Rsun.
- **A restart** with groups present: identical to the last digit.

**Docs updated:** README, `src/README`, PARAMS_README, `params.nml`, the tests
and benchmark READMEs, and the archive README.

## Measurements
These were done here, on 1000 stars with one thread and `rcoll` = 0, against
the previous commit. The test systems were made in Python, because there is
no Fortran compiler in this workspace (generator in `/tmp/cmp/mkics.py`, not
kept).

| System | Old code | New code |
|---|---|---|
| Plummer, 500 binaries (a = 10-1000 AU), 1 Myr | 385 s, max dE/E 2e-6 | 13 s, 3e-7 |
| Cold fractal (D = 1.6, Q = 0.3), 2 Myr, seed 1 | 470 s, dE/E 2.7e2 (meaningless) | 106 s, 4e-5 |
| seed 2 | 2.0e1 by 0.7 Myr (stopped) | 211 s, 1.1e-4 |
| seed 3 | - | 89 s, 1.4e-4 |

On the cold fractals the old code fails outright. Close, very eccentric or
unbound encounters between massive stars hit the shortest timestep level and
wreck the energy.

How the speed was found (from profiling):
- The weak perturbers' tidal tensor made it 3.5 times faster.
- The cheap landing made it a further 7 times faster.
- A step cap applied at every step was tried and dropped: it was 15 times
  slower for no gain.

## Loose ends, in order of importance
1. **Fractal energy errors of about 1e-4.** Every one of these looked at so
   far was a real, systematic hole (unbound plunges, eccentric plunges,
   stars passing through groups, groups meeting groups), and each has been
   fixed. What is left may be more of the same, or ordinary truncation
   error. The test that tells them apart was about to run when work
   paused: fractal 3 with `eta` = 0.005 and `etap` = 0.05. If the error
   drops, it is truncation.
2. **Groups as perturbers of other groups are point masses.** The
   merge-or-open rule bounds the damage. The full fix is to resolve a group
   perturber into its stars inside another group's integration. That needs
   the two groups' integrations to be coordinated, because they run side by
   side.
3. **Kepler binaries ignore perturbations below `gopen`** (1e-5). In a flyby
   that gave errors of about 1e-5 in the binary's elements, against 1e-10
   when the binary is regularised throughout. Lowering `gclose` and `gopen`
   to 1e-7 and 1e-6 (kira's level) would cost some speed; this has not been
   measured.
4. **Slow-down** (Mikkola & Aarseth 1996) for long-lived, weakly perturbed
   hierarchies. At present these are integrated orbit by orbit, which is
   accurate but can be slow. Frozen triples have no secular evolution.
5. **The perturber search is O(N)** each time a group's centre of mass
   starts a new step; with the tree at large N it should use the tree.
6. **The kira benchmark has not been rerun** with the new code. Next is
   `run_benchmarks.sh` with the user's real `makeics` initial conditions and
   kira.
7. **Tuning, unmeasured:** `NFEWFAC` (when a multiple pays), and the
   constants `ARCAP`, `GHORFAC`, `CAPFRAC` and `LANDTOL` in `ar.c`.
8. **Not yet copied to the claude.ai project:** only `src/ar.c`,
   `src/hendos.h`, `src/few.c` and `src/groups.c` have been updated there.
   The other changed files are up to date on GitHub and in the tarball.

## Debugging aids that found the problems above
- `HENDOS_DEBUG=1` prints the energy change at every group event; the value
  is field 9 of each `DEBUG` line.
- `HENDOS_TRACE="t1 t2"` gives the energy block by block. Use it from a
  restart file just before the problem, because a long window is very slow.

## Addendum, 26 November 2026: before the kira benchmark
A dry run of the benchmark on a smooth 200-star Plummer sphere showed
regularisation making HENDOS *less* accurate: largest dE/E 7e-7, against
2e-9 with `fewbody` off. Each group cost about 1e-7 when it was opened. The
cause was the perturber list overflowing. Mild unbound flybys at about
2e-3 pc had been regularised, had drifted apart within one centre-of-mass
step, and by then had more perturbers than the 64 the list held. Fixes:
- **Unbound pairs are regularised only for a plunge**: two-body pericentre
  less than a tenth of the current separation.
- **The centre of mass of an unbound pair steps at least every `etap` times
  the pair's flyby time r/v.**
- **A receding unbound pair is let go** once it is twice as far apart as when
  it was made a group (its separation then is kept in `bein`), or at
  gamma > `gfew`.
- **A group whose perturber list is full** while the weakest perturber kept
  is still strong is let go.
- **The perturber list holds 256** (`NPMAX`); each list is allocated only
  when a group is regularised.

Results (largest dE/E), before -> after:

| Run | Before | After |
|---|---|---|
| 200-star Plummer, no binaries | 7e-7 | 3.5e-8 |
| fractal 3 | 1.4e-4 | 8e-6 |
| Plummer with 500 binaries | 3e-7 | 4e-8 |

The 25 tests pass, also under the sanitisers.
