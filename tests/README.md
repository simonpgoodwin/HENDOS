# HENDOS regression tests

`run_tests.py` runs the `hermite` executable on small problems whose answers
are known and checks each result against a tolerance. Run it after every
change to the code and before trusting a production run.

```
cd HENDOS/run && make          # the tests use run/hermite
cd ../tests
python3 run_tests.py                 # all the tests
python3 run_tests.py --list          # what the tests are
python3 run_tests.py --only merger   # tests whose name contains 'merger'
python3 run_tests.py --keep          # keep the working files (in tests/work)
python3 run_tests.py --dry           # write the input files, run nothing
python3 run_tests.py --hermite /path/to/another/hermite
```

It needs Python 3 with numpy, and uses `analysis/hendosio.py` to read the
snapshots. Each test makes its own initial conditions, so `makeics` is not
needed. Set `OMP_NUM_THREADS` as for a normal run. The whole set should take
a few minutes. The exit status is 0 if everything passed, so it can be used in
a script.

Each test prints PASS or FAIL and, underneath, every number it checked with
its limit:

```
kepler_eccentric       PASS   0.4 s
    ok    energy error                                   -3.1e-07   (limit 1.0e-05)
```

If a test fails its working files (parameter file, initial conditions, screen
log, snapshots) are left in `tests/work/<test name>/` to look at.

## The tests

| Test | What it runs | What it checks |
|---|---|---|
| `kepler_circular` | two stars, circular orbit, 10 periods, integrated directly | energy; stars back where they started |
| `kepler_eccentric` | two stars, eccentricity 0.9, 10 periods, directly | energy |
| `binary_closed` | the same binary treated as one body | exact Kepler orbit: energy to 1e-10, position, one group |
| `pythagorean_direct` | Burrau's three-body problem to t = 10, directly | energy |
| `pythagorean_fewbody` | the same as a few-body group | energy; one group formed |
| `triple_frozen` | a stable hierarchical triple | frozen as two Kepler orbits; energy; inner orbit unchanged |
| `plummer_direct` | 256-star Plummer sphere, 2 Myr, direct | energy; virial ratio; number of stars |
| `plummer_neighbours` | the same with the neighbour scheme | energy |
| `plummer_tree` | the same with the tree | energy (looser: the tree is approximate) |
| `plummer_tide` | the same in the solar tidal field | Jacobi energy; tidal keys in the header |
| `header_flexible` | a run started from a snapshot with an unknown header key and an extra column added | keys follow the physics; the changed file is read correctly |
| `merger` | two stars collide head on beside a third | one star fewer; mass and momentum conserved; merger file, screen alert and header count; id of the more massive kept |
| `merger_off` | the same with `mergers = .false.` | the run stops with the collision alert |
| `restart` | a run stopped half way and restarted | ends where the unbroken run ends (positions to 1e-9 pc) |
| `restart_neighbours` | the same with the neighbour scheme | the same |

The limits are deliberately loose enough not to fail for harmless reasons
(a different compiler or number of threads) and tight enough to catch a real
mistake. **They were set from tests of the algorithms in Python, before the
Fortran had been compiled, so the first real run may show that some need
adjusting** - tighten any that pass with a lot to spare, and note the values
your build gives so that later changes can be compared with them.

## Adding a test

Write a function `test_<name>(t)` in `run_tests.py`; it is found and run
automatically. `t` is a `Tester`:

- `t.write_ics(fname, m, pos, vel, types=None)` writes initial conditions
  (Msun, pc, km/s).
- `t.run(runid, icfile, section__name=value, ...)` runs the code. Parameters
  not given take the values in `DEFAULTS` at the top of the file. It returns
  the run: `.last` is the last snapshot (a `hendosio.Snapshot`), `.snaps` the
  snapshot files, `.log` the screen output, `.file('_mergers.dat')` any other
  file of the run.
- `t.check(what, value, limit)` passes if |value| < limit;
  `t.check(what, value, expected, 'equal')` passes if they are equal.

Helpers for initial conditions: `binary(m1, m2, a, ecc)`, `plummer(n)`,
`pythagorean()`. The constants `G`, `TUNIT` and `VUNIT` are the code's units.

Add a test whenever something new goes into the code (stellar evolution,
primordial multiples), and whenever a bug is fixed: a test that fails with the
bug and passes without it stops the bug coming back.
