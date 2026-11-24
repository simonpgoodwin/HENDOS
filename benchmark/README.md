# Benchmarking HENDOS against Starlab's kira

Runs the same stars through both codes and compares how long each takes to
reach a given accuracy. HENDOS can do the same problem in several ways
(direct sums, the neighbour scheme, the tree, with or without regularised
groups), so it is run in one or more *modes*, and each mode is compared with
kira as a curve of its own.

**Keep this file up to date: if the scripts or what they measure change,
change this in the same edit.**

| File | What it is |
|---|---|
| `run_benchmarks.sh` | **start here**: builds everything, makes the stars, checks kira, runs the whole set of tests |
| `ics/*.nml` | the initial conditions of the tests, for `makeics` |
| `bench_kira.py` | runs both codes on one set of stars, times them, measures the energy errors, makes the table and plots |
| `starlabio.py` | converts a HENDOS file to Starlab's `dyn` format, and reads `dyn` files (kira output) |

Needs Python 3 with `numpy` and `matplotlib`, the `analysis` directory next to
this one, a C compiler with OpenMP and gfortran (to build `hermite` and
`makeics`), and Starlab's `kira`.

## The quick way
```
cd HENDOS/benchmark
./run_benchmarks.sh quick        one small test (a few minutes), to see that all is well
./run_benchmarks.sh              the standard set (an hour or so)
./run_benchmarks.sh big          the standard set and 16,000 stars (many hours)
KIRA=/path/to/kira SEEDS="1 2 3" ./run_benchmarks.sh
```
The script builds `hermite` and `makeics`, makes the initial conditions of
every test, runs the check described below, and stops there if the check
fails. Then it runs `bench_kira.py` on every test. Each test gets its own
directory, `results/<test>_s<seed>/`, and all the tables are collected in
`results/summary.txt`. Settings (kira's path, the seeds, threads) are given in
front of the command and are listed at the top of the script, as are the
tests:

| Test | Stars | To (Myr) | HENDOS modes | What it tests |
|---|---|---|---|---|
| `plummer_eq_1k` | 1000 equal-mass, Plummer, r_h = 1 pc | 10 | direct, neighbours | the basic integrators |
| `plummer_imf_1k` | 1000 from the Maschberger IMF, Plummer | 10 | direct, neighbours | the same with a mass spectrum |
| `fractal_cool_1k` | 1000, box fractal D = 1.6, Q = 0.3 | 5 | direct, neighbours, nofew, nobin | binaries and multiples forming, and what the binary and few-body treatments are worth |
| `plummer_eq_4k` | 4000 equal-mass, Plummer | 4 | direct, neighbours | where the neighbour scheme starts to pay |
| `plummer_eq_16k` (big only) | 16,000 equal-mass, Plummer | 1 | direct, neighbours, tree | large N, and the tree |

Each runs for a few crossing times. Single runs vary a lot (the first
benchmarks showed a factor of ten between runs that differed only in `eta`),
so the script runs **three seeds by default** (`SEEDS="1 2 3"`); give more for
anything you quote, or one (`SEEDS=12345`) for a quick look. The fractal fails to build for
roughly one seed in ten (see `src/ics/FRACTAL_README.md`); the script reports
that and skips the test for that seed.

## The kira command line and the check
The kira command is one line near the top of `bench_kira.py`:
```
KIRA_CMD = '{kira} -t {tend} -d {dt} -D {dt} -a {eta} -e 0 -u'
```
The options are as `kira_init.C` of Starlab 4.4 reads them: `-t` the length of
the run, `-d` the interval between log outputs, `-D` (with a number) the
interval between snapshots, `-a` the accuracy parameter eta, `-e` the
softening (0, as in HENDOS) and `-u` to switch on kira's unperturbed binary
motion (HENDOS treats unperturbed binaries specially too, so this is the like
for like comparison). Times are in N-body units, converted by the script.

The options were read from the Starlab source, but kira itself has not been
run with them here (there was no kira where this was written). So **run the
check first** - the script does it for you, or by hand:
```
cd HENDOS/run
python ../benchmark/bench_kira.py data/plum1k.dat --tend 5 --check
```
It runs kira and HENDOS for one output interval on the stars and checks that
kira ran, wrote snapshots that can be read, kept every star, reached the end
time (which shows the units are converted correctly) and conserved energy;
that HENDOS did the same; and that the two codes agree on where the stars are
after the interval. It ends with "All is well" or says what failed.

The scripts, the check and the driver were tested end to end here with a
stand-in for kira (a small leapfrog integrator that reads and writes dyn files
and takes kira's options), including all five HENDOS modes, two seeds, the
tables, the plots and the summary.

## Running bench_kira.py by hand
For your own initial conditions:
```
cd HENDOS/run
./makeics ics.nml                                   (writes data/plum1k.dat)
python ../benchmark/bench_kira.py data/plum1k.dat --tend 5 --check
python ../benchmark/bench_kira.py data/plum1k.dat --tend 5
python ../benchmark/bench_kira.py data/plum1k.dat --tend 5 --modes direct neighbours
```
Options:

| Option | Meaning | Standard |
|---|---|---|
| `--tend` | end time (Myr) | must be given |
| `--nout` | number of outputs and snapshots over the run | 20 |
| `--modes` | the ways HENDOS is run, one curve each (see below) | direct |
| `--eta-hendos` | HENDOS `eta` values to run in every mode | 0.04 0.02 0.01 |
| `--eta-kira` | kira `-a` values to run | 0.2 0.1 0.05 |
| `--hermite` | the HENDOS executable | `../run/hermite` |
| `--threads` | OpenMP threads for HENDOS | 1 |
| `--out` | directory for everything | `bench_<ic file name>` |
| `--skip-kira`, `--skip-hendos` | run only one of the codes | |
| `--kira` | the kira executable | `kira` (on your path) |
| `--check` | only run the check described above, then stop | |

## The HENDOS modes
| Mode | Forces | Kepler binaries | Regularised groups | Compare with kira? |
|---|---|---|---|---|
| `direct` | direct sums | on | on | yes - the like-for-like comparison |
| `neighbours` | direct sums with the Ahmad-Cohen neighbour scheme | on | on | yes |
| `tree` | Barnes-Hut tree | on | on | not like for like: the forces are approximate |
| `nofew` | direct sums | on | off | shows what regularisation gains (perturbed binaries and close encounters integrated directly) |
| `nobin` | direct sums | off | off | shows what the whole binary treatment gains |

kira treats unperturbed binaries and multiples specially too, so `direct` and
`neighbours` are the fair comparisons with it. `nofew` and `nobin` are there
to compare HENDOS with itself. In `neighbours` mode the regular-force
parameter `etar` is set to half of `eta` in every run, so the two are
tightened together (`ETAR_OVER_ETA` in the script).

The modes are the table `HENDOS_MODES` near the top of `bench_kira.py`; to add
one (a different `theta` for the tree, say), add a line there and, if it needs
a parameter that is not yet in it, a line to `HENDOS_PARAMS` just above.

## What it does
1. Converts the HENDOS initial conditions to a `dyn` file in standard N-body
   units (G = 1, total mass 1, virial radius 1), which is what kira expects,
   and saves the units so that kira's output can be put back into Msun, pc,
   km/s and Myr.
2. Runs kira, and HENDOS in each mode asked for, once for every accuracy value
   given, timing each run.
3. Works out the energy error of every run **in the same way for both
   codes**: by direct summation over the stars of the initial conditions and
   of every snapshot. Neither code's own energy report is used. The figure
   plotted is the **largest** |dE/E| at any snapshot; the error at the end is
   given too, but on its own it is a poor measure, because the error can
   change sign during a run and a run can end with a tiny error by luck
   (both codes did so in the first benchmarks).
4. For every HENDOS run, takes from the header of its last snapshot what the
   binary and few-body machinery did (binaries closed, few-body groups formed,
   triples frozen) and HENDOS's own energy error, and checks that the number
   of stars has not changed.
5. Writes a table and two plots.

HENDOS has two timestep parameters besides `eta`: `etap` (the pairwise
free-fall and flyby limit) and, with neighbours, `etar`. The script ties both
to `eta` (`etap` = 5 `eta`, `etar` = 0.5 `eta`; `ETAP_OVER_ETA` and
`ETAR_OVER_ETA` in `bench_kira.py`) so that all three tighten together. In the
first benchmarks `etap` was left at 0.2, and the close encounters then set a
floor to the accuracy: cutting `eta` from 0.04 to 0.01 reduced the error only
three times, where it should have been a hundred.

HENDOS is run with the collision check off (`rcoll = 0`), so there are **no
mergers**: kira without stellar radii does not merge stars, and a merger
removes orbital energy from the stars, which would spoil an energy comparison
made from the snapshots. Restart files are off. It uses one thread unless told
otherwise, because kira is a serial code.

## What you get
In the output directory:

| File | What |
|---|---|
| `results.txt`, `results.json` | wall-clock time, largest and final energy error of every run |
| `bench_time_error.png` | wall-clock time against largest energy error, one curve for kira and one for each HENDOS mode |
| `bench_halfmass.png` | half-mass radius against time for every run: the most accurate of each code or mode thick, the others thin |
| `ic.dyn`, `ic.dyn.units.json` | the stars in Starlab format, and the units |
| `kira_eta<..>.dyn`, `.log` | kira's snapshots and log for each run |
| `hendos_<mode>_eta<..>/` | HENDOS parameter file, log (verbose) and snapshots for each run |

## Reading the result
- **Compare at equal accuracy, not equal parameters.** `eta` does not mean the
  same thing in the two codes. On `bench_time_error.png` pick an energy error
  and read off the time for each code: the lower curve is faster at that
  accuracy. That is why several accuracy values are run.
- **Check the physics agrees.** The half-mass radius curves should lie on
  top of each other at first and then drift apart, as any two N-body runs of
  the same stars do. After a cold collapse, or once hard binaries form, the
  drift is chaotic and fast, so compare the *spread* of one code's curves (its
  runs at different `eta`) with the spread of the other's, and the averages
  over seeds - not one curve with one curve. A difference larger than both
  spreads, in every seed, would be worth looking into; a large difference
  from the very start means the two codes were not given the same problem.
- **Check every run reached the end.** The table says so if one did not.
- **Compare the two energy errors of each HENDOS run.** The table gives the
  error measured here from the snapshots and, in the notes, HENDOS's own
  (`own dE/E`). They should be close: HENDOS's own figure includes the
  internal energy of every group and its interaction with its perturbers,
  which the snapshots (all the stars) include too. A Kepler binary ignores
  perturbations below `gopen`, which shows up here as a small difference.
- **Look at what the binary machinery did.** The notes give the number of
  Kepler binaries closed, regularised groups formed and triples frozen. A run
  that is slow with none of these has not found its binaries. A cool fractal
  makes thousands of regularised groups in a few Myr (each close encounter
  is one, briefly), which is as it should be.

## A fair set of tests
- An equal-mass Plummer sphere (`imf = 'equal'`) at N = 1000, 4000 and 16000
  with `--modes direct neighbours`: no binaries to speak of, so this compares
  the basic integrators, and shows where the neighbour scheme starts to pay.
- The same with the standard mass function.
- A cool fractal (`distribution = 'fractal'`, `qvir = 0.3`), which forms
  binaries and multiples quickly, with `--modes direct nofew nobin`: this
  tests how each code copes with them, and shows what the binary treatment
  and regularisation in HENDOS are worth.
- When the initial conditions generator can make primordial binaries and
  triples, repeat the Plummer test with them: that is the case the
  regularised groups and frozen triples were written for.
- For HENDOS also try `--threads` equal to your number of cores, and
  `--modes tree` at large N, noting that neither is like for like with serial
  kira.

Run each test a few times with different seeds: single runs vary.

## The converter on its own
```
python starlabio.py <HENDOS file> <dyn file>
```
writes the dyn file and `<dyn file>.units.json`. In a script:
```python
from starlabio import hendos_to_dyn, read_dyn, Units
units, snap = hendos_to_dyn('data/plum1k.dat', 'plum1k.dyn')
out = [s.to_physical(units) for s in read_dyn('kira_output.dyn')]
```
`read_dyn` returns single stars in the frame of the whole system even where
kira has written binaries and multiples as nested groups.
