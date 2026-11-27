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
  less than a tenth of the current separation (now a thousandth - see the
  27 Nov addendum).
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


## Addendum, 27 Nov: the Plummer floor

The user's kira benchmark showed HENDOS direct stuck at an energy error of
1e-7 to 1e-6 on smooth Plummer spheres whatever eta was, where kira got to
1e-8. On a 1000-star equal-mass Plummer sphere (10 Myr, eta 0.01) the
error was 2.1e-7 with `fewbody` on and 1.4e-8 with it off, and lowering
`gpert` made no difference.

The cause was flicker. With HENDOS_DEBUG, one wide eccentric binary
(apocentre 0.03 pc, pericentre 0.004 pc, period about 0.14 Myr) was being
made a regularised group at every pericentre, where its gamma dropped
below `gfew`, and let go on the way out because its perturber list was
full (at 0.03 pc, every star in the cluster is a perturber above 100
`gpert`). That happened 5400 times in 10 Myr, and each close and let-go
left 1e-9 to 5e-8 behind.

Fixes (groups.c):
- **`overfull()`**: a bound pair is not regularised if more than `NPMAX`
  bodies would have gamma > 100 `gpert` at its apocentre size (so the
  list-full rule would let it go at once), unless it plunges deeply. The
  Hermite integrator follows such pairs perfectly well.
- **`PLUNGE` = 1e-3** (a `>>> CHANGE HERE` define): an unbound pair (or an
  overfull bound one) is regularised only if its pericentre is under a
  thousandth of its separation. On its own this made no difference to the
  Plummer sphere, but the plunges that matter in cold fractals have
  pericentres of 1e-6 of the separation, so nothing is lost.

Results:

| run | before | after |
|---|---|---|
| 1000-star Plummer (fewbody off: 1.4e-8) | 2.0e-7 | 5.3e-8 (2.9e-9 in a second build) |
| Plummer with 500 binaries | 4.2e-8 | 4.2e-8 |
| fractals 1-3 | 1e-5 to 1e-4 | 1.7e-5 to 3.2e-5 |

Only about ten group events are left in the Plummer run. The largest
cost is opening a group (up to 1.6e-8), which is the next thing to look at
if the floor matters.

The 25 tests pass.


## Addendum, 27 Nov (later): the neighbour scheme in sub-clumps

Asked: does the Ahmad-Cohen scheme go wrong in the sub-clumps (10 to 100
stars) of a fractal? In the benchmark the neighbour runs were 2 to 10
times worse than direct on all nine fractal runs.

Test: fractal seed 1 (1000 stars), with the regular-force prediction
error (how far the regular force had moved from its linear prediction by
each regular step, as a fraction of the total force) written out, and
each star matched to its clump by friends-of-friends (linking length
0.1 pc) on the snapshots.

- **The clumps themselves are fine.** Stars in clumps of 10-30 and 31-100
  have median errors of 1e-6 to 3e-6 (99% below 6e-4), no worse than
  stars in big clumps. Stars between the clumps have the largest relative
  errors (median 3e-4), but their forces are weak and smooth.
- **There was a bug in the list bookkeeping.** 0.1% of regular steps found
  the regular force wrong by more than 10% of the total force. When stars
  are made a group (or merged), the body that takes their places in the
  neighbour lists carries the mass of all of them. A star whose list held
  only some of them then felt the others twice - once in the irregular
  force (through the c.o.m.), once in the regular force - until its next
  regular step. In sub-clumps, where groups are made all the time among
  stars with about 37 neighbours, that happened constantly. Fixed by
  `acjoining` (ac.c), called from `closegroup` and `mergestars`: the
  pull of the stars that were not in the list is moved from the regular to
  the irregular force. Regular steps wrong by more than 10%: 0.1% -> none;
  by more than 1%: 0.28% -> 0.001%.

Results (fractal, 1 Myr, `eta` 0.01, `etap` 0.05, `etar` 0.005):

| seed | direct | neighbours, fixed |
|---|---|---|
| 1 | 7.0e-5, 181 s | 1.3e-5, 165 s (2e-4 before the fix) |
| 2 | 1.8e-4, 348 s | 1.7e-5, 303 s |
| 3 | 6.8e-6, 177 s | 1.0e-4, 184 s |

The occasional 1e-4 jumps, which turn up in either mode, are not the
neighbour scheme: they are group-group bookkeeping (next point).

Also: `NACMIN` raised from 150 to 200 (user's request); the neighbour
restart test moved from 192 to 256 stars so that it still tests the
scheme; the benchmark's `etar` is now `eta`/4 (Plummer 1k: 2e-6 at
`eta`/2, 6e-8 at `eta`/4, direct 5e-8; 24 s against 37 s for direct).

**Lead, not fixed: groups that grow to 1 pc.** The large jumps at group
openings in fractals come from regularised groups that have swung out to
about 1 pc (wide pairs can be bound in a cold fractal, and gamma stays
small in its empty parts). Such a group sees the other groups (dozens, a
few pc away) only as point masses, and the missing terms made the
energy wrong by up to 6e-4 until it was opened. A wider group-group test
(any group, with both sizes, 2(s_i+s_j)^3 m_j/(m_i d^3) > gfew, used for
opening and making groups alike) was tried: fractal seed 1 went from
7e-5 to 2e-5 but seed 2 from 2e-5 to 2e-4, with hundreds of openings
with jumps of up to 5.6e-4. It was taken out; the attempt is kept in
`/tmp/groups_ggwide.c` for this session only. A cleaner fix is probably
a size limit on regularised groups, or resolving group perturbers into
their stars.

The 25 tests pass.


## Addendum, 27 Nov (evening): multiples and sub-groups in clumpy clusters

Asked: make HENDOS deal properly with multiples and sub-groups before the
next benchmark, since the use cases have many multiple-star interactions
in sub-clumps.

### How it was tested
Single fractal runs are too chaotic to compare one version with another:
any change (even a diagnostic that rebuilds a perturber list) gives a
different realisation, and one bad event in one realisation decides the
energy error. So each version was run on **nine** cold fractals (1000
stars, D = 1.6, Q = 0.3, seeds 1-9; generator as before) for 0.2 Myr, the
violent early phase, then on seeds 1-3 for 1 Myr, plus the two Plummer
checks. Bad events were found by rerunning the identical run with a
read-only energy check (`grpcomp` must not `arshift` inside the check, or
the trajectory changes), localising the loss in time, and auditing what was
mispriced just before it.

### What was wrong, and the fixes (all in groups.c and ar.c)
1. **The two-groups test was made only if both centres of mass were on the
   same block** - seldom, with block steps. Two hard binaries could then
   come within 2.5 sizes of each other as point masses (up to 9e-4 of the
   energy each time). Now made on any block; so are the tests with a
   group's nearest single star, and the pair test. A body between steps is
   taken into a group at its predicted state (`syncbody`).
2. **Except a regularised group between steps**: it is known only through
   its interpolated dense output, which for a hard binary can be wrong in
   energy by 1e-4 of its binding energy - a 38 Msun binary bound more
   tightly than the whole cluster cost 2e-4. Such a group is taken in only
   at the end of its own step.
3. **The two-groups test looked only at the nearest body, and only at this
   group's size.** Now every group is looked at, with both sizes:
   g = 2 (s_i + s_j)^3 m_j / (m_i d^3) > `gfew` (`grpworst`). The same test
   is used when opening a group and when making a pair or a multiple, so a
   group it opens is not made again at once.
4. **When two groups cannot be merged, the wider one is opened**, never a
   hard binary because a wide pair came by. (Opening whichever group
   noticed first put hard binaries back into the Hermite integrator: 0.2.)
5. **Merging needed the merged group to be isolated to `gfew`.** In dense
   knots two binaries meeting never are. A group meeting another body is now
   accepted below `GMULT` `gfew` = 0.05 (its strong perturbers are followed
   exactly, so this costs time, not accuracy). Three single stars still need
   `gfew` (soft triples cost 2e-7 in a Plummer sphere otherwise).
6. **Perturber lists grow as needed** (any number of weak perturbers, at
   most `NPMAX` = 256 followed exactly). With a fixed list of 256, groups in
   dense regions were let go when it filled, and `overfull` (this morning)
   refused tight pairs in dense knots altogether (one at gamma 4e-6, 18 AU
   apart). Both rules are gone. The 1k Plummer sphere improved too.
7. A loop made and opened the same triple 22,000 times in 0.2 Myr; fixed by
   (3).

Tried and parked (in `git stash`): growing groups up to 16 stars until
isolated. In these fractals no set of up to 16 stars round a knot is
isolated (clumps within clumps), so it found nothing and cost 2-3 times as
much.

### Results (largest |dE/E|)
Nine fractals, 0.2 Myr (`eta` 0.01, `etap` 0.05):

| seed | 1 | 2 | 3 | 4 | 5 | 6 | 7 | 8 | 9 |
|---|---|---|---|---|---|---|---|---|---|
| before | 2.4e-5 | **1.6e-4** | 4.5e-7 | 3.6e-6 | 1.7e-6 | 4.0e-6 | 1.9e-5 | 2.2e-6 | 9.8e-6 |
| after | 2.6e-5 | 2.6e-5 | 6.5e-7 | 6.2e-6 | 9.2e-6 | 2.0e-6 | 6.0e-6 | 1.4e-6 | 4.8e-6 |

Fractals 1-3, 1 Myr: before 7.0e-5, 1.8e-4, 6.8e-6 (181, 348, 177 s);
after 2.6e-5, 4.0e-5, 1.5e-5 (211, 362, 165 s). Seed 1 with neighbours
(`etar` 0.0025): 2.8e-5 in 167 s.

Plummer 1k (10 Myr): 5.3e-8 -> 6.3e-9. Plummer with 500 binaries: 4e-8 ->
1.3e-8.

The worst case is now a few 1e-5 instead of 2e-4. What is left, as far as
is known, is the point-mass treatment of groups for each other below the
`gfew` threshold, and sub-clumps where every group is perturbed. The full
fix for the first is still to resolve group perturbers into their stars
inside each other's integration.

The 25 tests pass.
