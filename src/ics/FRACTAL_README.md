# The HENDOS box-fractal generator

How `makefractal` builds a fractal, what the core size is, and what the
velocity options do. There are two routines, both in `fractalsubs.f90`:
`makefractal` builds the fractal and makes the inherited and random
velocities, and `larsonfractal` makes a fractal with velocities that follow
Larson's law (it calls `makefractal` for the positions, so the two always
build the fractal in the same way). The program that calls them is `makeics`
(`makeics.f90`). The inputs are in the `&fractal` section of `run/ics.nml`
and are listed in `run/PARAMS_README.md`, and the Larson velocities are derived in `larson.tex`.

**Keep this file up to date: if either routine changes, change this in the same
edit.**

## How the fractal is built
The method is that of Goodwin & Whitworth (2004, A&A, 413, 929).

1. Start with one cube of half-side `radius`, with a single "parent" at its
   centre.
2. Divide each cube into `ndiv` along each side, giving `ndiv`^3 smaller
   cubes with a "child" at the centre of each. `ndiv` is 2 as standard, which
   gives 8 children.
3. Keep some of the children and delete the rest. The number kept is the
   nearest whole number to `ndiv`^`fdim`, chosen at random among the children.
4. Move each surviving child a little away from the centre of its cube (the
   "noise"), so the stars do not lie on a regular grid.
5. The children become the parents of the next level. Go back to step 2 until
   the stopping rule is met.
6. Cut out the sphere of radius `radius` and throw away what is outside.
7. Remove stars at random until exactly `n` are left.

Each pass through steps 2 to 5 is one **level**. The cubes shrink by a factor
`ndiv` at each level, so after L levels they are `ndiv`^L times smaller than
the first cube.

### Fractal dimension
Because a whole number of children must be kept, only some fractal dimensions
are made exactly: those equal to log(k)/log(`ndiv`) for a whole number k. The
program prints the one you will get. With `ndiv` = 2:

| `fdim` asked for | children kept (of 8) | dimension made |
|---|---|---|
| about 1.6 | 3 | 1.58 |
| 2.0 | 4 | 2.00 |
| about 2.3 | 5 | 2.32 |
| about 2.6 | 6 | 2.58 |
| about 2.8 | 7 | 2.81 |
| 3.0 | 8 | 3.00 (uniform) |

### The number of boxes per side, `ndiv`
`ndiv` is the second number on the fractal dimension line of the input file.
2 is the standard value and the one the group's results use. A larger value
changes how the fractal looks without changing its dimension:
- more dimensions can be made exactly. With `ndiv` = 3 a box has 27 children
  and the dimension is log(k)/log(3), so 1.63 (6 kept), 2.00 (9), 2.10 (10),
  2.52 (16) and so on;
- each level is a bigger step in scale (a factor 3 instead of 2), so there
  are fewer levels and the structure is coarser;
- the number of stars made jumps by more at each level, so more are thrown
  away at the end.

The tables below are for `ndiv` = 2.

## The core size
`coresize` is the **smallest length scale in the fractal**. It is not the core
radius of a cluster. It controls two things: when the dividing stops (step 5)
and how big the noise is (step 4). There are two modes.

### `coresize` = 0 - the number of stars sets the smallest scale
This is the original `makefractal`, and the standard setting.

- **Stopping:** dividing stops as soon as there are more than 3.5 times `n`
  children. The factor 3.5 leaves enough to survive cutting out the sphere.
- **Noise:** each child is moved by a random amount, up to a quarter of the
  width of its cube, in each coordinate. The noise shrinks with the cubes.
- **Smallest scale:** whatever the cubes have shrunk to when dividing stops,
  so it depends on `n`. For `n` = 1000, `fdim` = 1.6 and `radius` = 1 pc
  there are 8 levels and the smallest cubes are about 0.008 pc across.

### `coresize` > 0 - you set the smallest scale
This is `makefractal2`.

- **Stopping:** dividing stops when the children of the next level would be
  closer together than `coresize`. The fractal has structure from `radius`
  down to about `coresize` and none below it.
- **Noise:** each child is moved by a Gaussian random amount of dispersion
  `coresize` in each coordinate, the same at every level. Because it is added
  at every level, the total blurring after L levels is about `coresize` times
  the square root of L.
- **Number of stars made:** set by `coresize`, not by `n`. Each level
  multiplies the number by the number of children kept, and the extra stars
  are thrown away at the end. So `coresize` must be small enough to leave at
  least `n` stars in the sphere, but making it much smaller than needed wastes
  memory and time.

Example for `fdim` = 1.6 (3 children kept) and `radius` = 1 pc:

| `coresize` (pc) | levels | stars made before the sphere is cut |
|---|---|---|
| 0.01  | 7 | 2,187 |
| 0.005 | 8 | 6,561 |
| 0.003 | 9 | 19,683 |
| 0.001 | 10 | 59,049 |

Cutting the sphere from a very clumpy fractal can remove most of the stars, so
aim for several times `n`: for `n` = 1000 use 0.003 pc, not 0.005 pc.

### If the program stops with `nsmax<nstar`
Too few stars were left inside the sphere. With `coresize` = 0 change the
random number seed (this happens for roughly one seed in ten at `fdim` = 1.6,
and less often for larger `fdim`). With `coresize` > 0 use a smaller
`coresize`.

## The velocity options
Chosen by `velocities` in the input file, and independent of the core size.

For `standard`, `inherit` and `random` the program scales all the velocities
by one factor afterwards to give the virial ratio asked for, so the choice
sets the *pattern* of the velocities, not their size.

`larson` is the exception. Its velocities are set to the observed
size-linewidth relation in km/s: two stars `d` pc apart differ by
`vlarson` x sqrt(`d`) km/s in each velocity component (rms), and the standard
`vlarson` = 1 is the relation of Solomon et al. (1987), 0.72 (R/pc)^0.5 km/s
for the 1D velocity dispersion of a region of radius R. The virial ratio is
then not something you choose: the program reports what it comes out as. It
rises as radius squared over mass. For 1000 stars from the standard IMF (about
370 Msun) with `fdim` = 1.6 it is about 0.1 for a radius of 0.5 pc, 0.4 for
1 pc and 1.5 for 2 pc, varying by a third or more from one fractal to the next.

| `velocities` | routine | how the velocities are made | velocity difference against separation |
|---|---|---|---|
| `larson` | `larsonfractal` | drawn all together as correlated Gaussian numbers once the positions are fixed | power 0.5 (Larson's relation) |
| `inherit` | `makefractal`, `vtype` `f` | each child takes its parent's velocity plus a random part proportional to the size of its cube | power about 0.8 |
| `standard` | `makefractal`, `vtype` `s` | each child takes its parent's velocity plus a random part of the same size at every level (the original `makefractal`) | power about 0.1 |
| `random` | `makefractal`, `vtype` `u` | an independent Gaussian velocity for every star | none |

The powers are from tests of a Python copy of the routine (`n` = 1000, `fdim`
= 1.6, ten fractals each); a single fractal can differ by about 0.1.

`larson` needs an `n` by `n` matrix: 8 `n`^2 bytes of memory (0.8 GB for
10,000 stars) and a time that grows as `n`^3. The others are fast for any `n`.

## What the program adds
The routines return positions in a sphere of radius 1 and velocities in
arbitrary units. `makeics` then
- sets every mass to unity, then gives each star its mass from the IMF,
- scales the positions so the sphere has the total radius asked for,
- moves the stars to their centre-of-mass frame,
- scales the velocities to the virial ratio (or, for `larson`, to the
  size-linewidth relation),
- writes the file in the HENDOS snapshot format.

Random numbers are from `ran3` (`icrandom.f90`), so one seed gives the same
fractal on any machine.

## Differences from the routines as supplied
The two routines supplied, `makefractal` and `makefractal2`, are rearranged
into `makefractal` (positions, and the inherited and random velocities of
both) and `larsonfractal` (the correlated velocities of `makefractal2`).
Every line that differs is marked `! HENDOS:` in `fractalsubs.f90`.

- The sphere cut keeps positions and velocities together and no longer loses
  the last star. (The original `makefractal` moved positions without their
  velocities.)
- With a core size, the position noise no longer shifts every star by
  -`coresize`/2 at each level.
- `coresize` = 0 no longer divides by zero, and gives the original
  `makefractal` positions.
- The Larson arrays are allocated rather than automatic, so large `n` does not
  overflow the stack.
- `REAL(start)**fdim` (an undeclared variable) is `REAL(ndiv)**fdim`.
- The unused starting level (`lstart`, `nstart`) is removed.
- Random numbers are from `ran3` instead of `ran2`.
- `ndiv` is an argument and an input (it was set to 2 inside the routines),
  and the arrays are sized for any `ndiv`, with a check that they are not
  overrun.
