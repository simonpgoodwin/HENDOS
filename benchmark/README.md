# Benchmarking HENDOS against Starlab's kira

Runs the same stars through both codes and compares how long each takes to
reach a given accuracy.

**Keep this file up to date: if the scripts or what they measure change,
change this in the same edit.**

| File | What it is |
|---|---|
| `bench_kira.py` | runs both codes, times them, measures the energy errors, makes the table and plots |
| `starlabio.py` | converts a HENDOS file to Starlab's `dyn` format, and reads `dyn` files (kira output) |

Needs Python 3 with `numpy` and `matplotlib`, the `analysis` directory next to
this one, HENDOS built (`make` in `run/`), and `kira` on your path.

## Before the first run
Two things could not be checked where this was written, because neither kira
nor a compiled HENDOS was available there:

1. **The kira command line.** It is one line near the top of `bench_kira.py`:
   ```
   KIRA_CMD = 'kira -t {tend} -d {dt} -D {dt} -a {eta}'
   ```
   Check the options against `kira --help` for your version, and add anything
   your runs need (for example the option for unperturbed binary motion).
2. **The dyn files.** Try the converter on its own and let kira read the
   result for a short run:
   ```
   python starlabio.py ../run/data/plum1k.dat test.dyn
   kira -t 1 < test.dyn > /dev/null
   ```

The scripts themselves were tested end to end with stand-in programs in place
of the two codes, including kira-style output with a binary written as a
nested pair.

## Running it
```
cd HENDOS/run
./makeics ics.nml                                   (writes data/plum1k.dat)
python ../benchmark/bench_kira.py data/plum1k.dat --tend 5
```
Options:

| Option | Meaning | Standard |
|---|---|---|
| `--tend` | end time (Myr) | must be given |
| `--nout` | number of outputs and snapshots over the run | 20 |
| `--eta-hendos` | HENDOS `eta` values to run | 0.04 0.02 0.01 |
| `--eta-kira` | kira `-a` values to run | 0.2 0.1 0.05 |
| `--hermite` | the HENDOS executable | `../run/hermite` |
| `--threads` | OpenMP threads for HENDOS | 1 |
| `--out` | directory for everything | `bench_<ic file name>` |
| `--skip-kira`, `--skip-hendos` | run only one of the codes | |

## What it does
1. Converts the HENDOS initial conditions to a `dyn` file in standard N-body
   units (G = 1, total mass 1, virial radius 1), which is what kira expects,
   and saves the units so that kira's output can be put back into Msun, pc,
   km/s and Myr.
2. Runs each code once for every accuracy value given, timing each run.
3. Works out the energy error of every run **in the same way for both
   codes**: by direct summation over the stars of the initial conditions and
   of the last snapshot. Neither code's own energy report is used.
4. Writes a table and two plots.

HENDOS is run with the collision check off (`rcoll = 0`), because kira does
not stop for collisions, and with one thread unless told otherwise, because
kira is a serial code.

## What you get
In the output directory:

| File | What |
|---|---|
| `results.txt`, `results.json` | wall-clock time and energy error of every run |
| `bench_time_error.png` | wall-clock time against energy error, one curve for each code |
| `bench_halfmass.png` | half-mass radius against time for the most accurate run of each |
| `ic.dyn`, `ic.dyn.units.json` | the stars in Starlab format, and the units |
| `kira_eta<..>.dyn`, `.log` | kira's snapshots and log for each run |
| `hendos_eta<..>/` | HENDOS parameter file, log and snapshots for each run |

## Reading the result
- **Compare at equal accuracy, not equal parameters.** `eta` does not mean the
  same thing in the two codes. On `bench_time_error.png` pick an energy error
  and read off the time for each code: the lower curve is faster at that
  accuracy. That is why several accuracy values are run.
- **Check the physics agrees.** The two half-mass radius curves should lie on
  top of each other at first and then drift apart slowly, as any two N-body
  runs of the same stars do. A large difference from the start means the two
  codes were not given the same problem.
- **Check every run reached the end.** The table says so if one did not.

## A fair set of tests
- An equal-mass Plummer sphere (`imf = 'equal'`) at N = 1000, 4000 and 16000:
  no binaries, so this compares the basic integrators.
- The same with the standard mass function.
- A cool fractal (`distribution = 'fractal'`, `qvir = 0.3`), which forms
  binaries quickly: this tests how each code copes with them.
- For HENDOS also try `--threads` equal to your number of cores, and the tree
  (edit `HENDOS_PARAMS` in `bench_kira.py` to add a `&tree` section), noting
  that neither is like for like with serial kira.

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
