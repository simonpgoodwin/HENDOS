"""
HENDOS - HErmite N-body Dynamics Of Systems
Make a small set of TEST files for trying out the viewer.

    python maketestdata.py [directory]        (default: testdata)
    python maketestdata.py [directory] 5      a test ENSEMBLE of 5 runs as well

THIS IS NOT A HENDOS RUN. It is a crude softened leapfrog integration of a
few hundred stars, written only so that there is something in the right
file format to look at before a real run exists. The files are marked as
test data in their headers. Do not use them for science.

It writes  demo_snap_<time>Myr.dat  and  demo_ana_lagrangian.dat  in the
HENDOS formats (NBODY-SNAPSHOT 2 and NBODY-ANALYSIS 1). With a number after
the directory it also writes that many smaller runs, demoens_e001, _e002, ...
which differ only in their random number seed, named and labelled as hermite
names and labels the members of an ensemble, for trying hendosensemble.py.
"""
import os
import sys

import numpy as np

G = 4.300917270619619e-3          # pc (km/s)^2 / Msun
KMS_TO_PC_MYR = 1.0227            # 1 km/s in pc/Myr


def header_line(key, value):
    return '# ' + key.ljust(30) + '= ' + value + '\n'


def accel_phi(m, r, eps):
    d = r[None] - r[:, None]
    d2 = (d ** 2).sum(2) + eps ** 2
    np.fill_diagonal(d2, np.inf)
    inv = d2 ** -0.5
    a = G * (m[None, :, None] * d * (inv ** 3)[:, :, None]).sum(1)
    phi = -G * (m[None] * inv).sum(1)
    return a, phi


def main(outdir, runid='demo', seed=42, n=300, nsnap=41, ensemble=None):
    """One test run. ensemble = (name, member, size) marks it as a member."""
    rng = np.random.default_rng(seed)
    tend, eps = 10.0, 0.02
    os.makedirs(outdir, exist_ok=True)

    # Plummer-like positions, IMF-like masses, a cold start so it evolves
    x = rng.random(n)
    rad = 0.6 / np.sqrt(x ** (-2 / 3) - 1)
    rad = np.minimum(rad, 5.0)
    ct = 2 * rng.random(n) - 1
    ph = 2 * np.pi * rng.random(n)
    st = np.sqrt(1 - ct ** 2)
    r = rad[:, None] * np.column_stack([st * np.cos(ph), st * np.sin(ph), ct])
    m = 0.1 * (1 - rng.random(n)) ** -0.75
    m = np.minimum(m, 20.0)
    v = rng.standard_normal((n, 3))
    types = np.full(n, 'star', dtype=object)
    heavy = np.argsort(m)[::-1]
    types[heavy[:2]] = 'black_hole'
    types[heavy[2:5]] = 'neutron_star'
    types[heavy[5:11]] = 'white_dwarf'
    types[np.argsort(m)[:8]] = 'planet'
    mtot = m.sum()
    r -= (m[:, None] * r).sum(0) / mtot
    v -= (m[:, None] * v).sum(0) / mtot
    a, phi = accel_phi(m, r, eps)
    ekin = 0.5 * (m * (v ** 2).sum(1)).sum()
    epot = 0.5 * (m * phi).sum()
    v *= np.sqrt(0.3 * abs(epot) / ekin)
    e0 = 0.5 * (m * (v ** 2).sum(1)).sum() + epot

    dt = 2.0e-3                      # Myr
    nsub = int(round(tend / (nsnap - 1) / dt))
    lag_rows = []
    for isnap in range(nsnap):
        time = isnap * nsub * dt
        a, phi = accel_phi(m, r, eps)
        ekin = 0.5 * (m * (v ** 2).sum(1)).sum()
        epot = 0.5 * (m * phi).sum()
        rr = np.sqrt((r ** 2).sum(1))
        order = np.argsort(rr)
        cum = np.cumsum(m[order]) / mtot
        fr = [0.01, 0.02, 0.05, 0.1, 0.2, 0.3, 0.5, 0.7, 0.9]
        lag = [rr[order][min(np.searchsorted(cum, f), n - 1)] for f in fr]
        sig = np.sqrt((m * (v ** 2).sum(1)).sum() / mtot)
        lag_rows.append([time] + lag + [sig])

        it = int(round(time * 1e4))
        name = os.path.join(outdir, '%s_snap_%06d.%04dMyr.dat' % (runid, it // 10000,
                                                                 it % 10000))
        with open(name, 'w') as f:
            f.write(header_line('format', 'NBODY-SNAPSHOT 2'))
            f.write(header_line('creator', 'maketestdata.py - TEST DATA, NOT A HENDOS RUN'))
            f.write(header_line('runid', runid))
            if ensemble:
                f.write(header_line('ensemble', ensemble[0]))
                f.write(header_line('ensemble_member', str(ensemble[1])))
                f.write(header_line('ensemble_size', str(ensemble[2])))
            f.write(header_line('n', str(n)))
            f.write(header_line('total_mass_Msun', '%.16E' % mtot))
            f.write(header_line('units', 'mass Msun, length pc, velocity km/s, '
                                         'time Myr, energy Msun (km/s)^2'))
            f.write(header_line('time_Myr', '%.16E' % time))
            f.write(header_line('energy_kinetic', '%.16E' % ekin))
            f.write(header_line('energy_potential', '%.16E' % epot))
            f.write(header_line('energy_total', '%.16E' % (ekin + epot)))
            f.write(header_line('energy_error', '%.16E' % ((ekin + epot - e0) / abs(e0))))
            f.write(header_line('virial_ratio', '%.16E' % (-ekin / epot)))
            f.write(header_line('columns', 'id m x y z vx vy vz type phi dt'))
            f.write(header_line('column_units',
                                '- Msun pc pc pc km/s km/s km/s - (km/s)^2 yr'))
            f.write('# end_header\n')
            step = 1e6 * 0.02 * np.sqrt((eps ** 2 + rr ** 2) ** 1.5 / (G * mtot))
            for i in range(n):
                f.write('%10d' % (i + 1)
                        + ''.join('%24.16E' % q for q in (m[i], *r[i], *v[i]))
                        + ' ' + types[i].ljust(12)
                        + '%16.8E%16.8E\n' % (phi[i], step[i]))
        # leapfrog to the next snapshot
        if isnap < nsnap - 1:
            for _ in range(nsub):
                v += 0.5 * dt * a * KMS_TO_PC_MYR
                r += dt * v * KMS_TO_PC_MYR
                a, _ = accel_phi(m, r, eps)
                v += 0.5 * dt * a * KMS_TO_PC_MYR

    with open(os.path.join(outdir, runid + '_ana_lagrangian.dat'), 'w') as f:
        f.write(header_line('format', 'NBODY-ANALYSIS 1'))
        f.write(header_line('analysis', 'lagrangian'))
        f.write(header_line('creator', 'maketestdata.py - TEST DATA, NOT A HENDOS RUN'))
        f.write(header_line('runid', runid))
        f.write(header_line('columns',
                            'time r1% r2% r5% r10% r20% r30% r50% r70% r90% sigma3d'))
        f.write(header_line('column_units', 'Myr pc pc pc pc pc pc pc pc pc km/s'))
        f.write('# end_header\n')
        for row in lag_rows:
            f.write('%16.6f' % row[0] + ''.join('%15.6E' % q for q in row[1:]) + '\n')
    print('wrote %d test snapshots and one analysis file of %s to %s'
          % (nsnap, runid, outdir))


if __name__ == '__main__':
    directory = sys.argv[1] if len(sys.argv) > 1 else 'testdata'
    main(directory)
    if len(sys.argv) > 2:
        size = int(sys.argv[2])
        for member in range(1, size + 1):
            main(directory, runid='demoens_e%03d' % member, seed=100 + member, n=120,
                 nsnap=21, ensemble=('demoens', member, size))
