"""
HENDOS - HErmite N-body Dynamics Of Systems
Moving stars between HENDOS files and Starlab's 'dyn' format.

    hendos_to_dyn(icfile, dynfile)   write a HENDOS initial conditions file
                                     as a Starlab dyn file, in N-body units
    read_dyn(path)                   read every snapshot in a dyn file
                                     (e.g. the output of kira)
    Units                            the conversion between the two

Starlab works in N-body units: G = 1, and kira also assumes the total mass
is 1 and the virial radius is 1 when it sets some of its internal scales.
hendos_to_dyn therefore uses the standard units
    mass unit    = total mass of the system
    length unit  = virial radius,  G M^2 / (2 |W|),  W = potential energy
    G            = 1
which fixes the velocity and time units. They are returned (and can be
saved) so that kira's output can be put back into Msun, pc, km/s and Myr.

The dyn format is nested text:
    (Particle                  the root: the whole system
      N = 3
    (Log ... )Log
    (Dynamics
      system_time  =  0
      m  =  1
      r  =  0 0 0
      v  =  0 0 0
    )Dynamics
    (Hydro )Hydro  (Star )Star
    (Particle                  one of these for each star
      i = 1
      N = 1
      ...
    )Particle
    )Particle
kira writes binaries as a Particle holding two Particles, with the
positions and velocities of the two measured from their parent. read_dyn
follows that, so it always returns single stars in the frame of the root.

Needs only numpy.
"""
import json
import os
import sys

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'analysis'))
from hendosio import read_snapshot              # noqa: E402

G_PHYS = 4.300917270619619e-3      # G in pc (km/s)^2 / Msun
PC_PER_KMS_IN_MYR = 0.9777922217   # 1 pc / (1 km/s) in Myr


def potential_energy(m, pos, g=1.0, chunk=512):
    """Potential energy of point masses (negative), by direct summation."""
    n = len(m)
    w = 0.0
    for a in range(0, n, chunk):
        d = pos[a:a + chunk, None, :] - pos[None, :, :]
        r = np.sqrt((d ** 2).sum(2))
        rows = np.arange(a, min(a + chunk, n))
        r[rows - a, rows] = np.inf
        w -= 0.5 * g * (m[a:a + chunk, None] * m[None, :] / r).sum()
    return w


def kinetic_energy(m, vel):
    return 0.5 * (m * (vel ** 2).sum(1)).sum()


class Units:
    """N-body units (G = 1) in physical units.

    mass (Msun), length (pc), velocity (km/s), time (Myr): the size of one
    N-body unit of each.
    """

    def __init__(self, mass, length):
        self.mass = float(mass)
        self.length = float(length)
        self.velocity = float(np.sqrt(G_PHYS * mass / length))
        self.time = float(length / self.velocity * PC_PER_KMS_IN_MYR)

    def save(self, path):
        with open(path, 'w') as f:
            json.dump(dict(mass_Msun=self.mass, length_pc=self.length,
                           velocity_km_s=self.velocity, time_Myr=self.time), f, indent=1)

    @classmethod
    def load(cls, path):
        with open(path) as f:
            d = json.load(f)
        return cls(d['mass_Msun'], d['length_pc'])


def hendos_to_dyn(icfile, dynfile):
    """Write a HENDOS initial conditions (or snapshot) file as a Starlab dyn file.

    The stars keep their order: star i in the dyn file is line i of the
    HENDOS file. Returns (units, snapshot): the Units used and the HENDOS
    snapshot that was read.
    """
    snap = read_snapshot(icfile)
    mtot = snap.m.sum()
    com = (snap.m[:, None] * snap.pos).sum(0) / mtot
    cov = (snap.m[:, None] * snap.vel).sum(0) / mtot
    pos = snap.pos - com
    vel = snap.vel - cov
    w = potential_energy(snap.m, pos, G_PHYS)
    units = Units(mtot, G_PHYS * mtot ** 2 / (2.0 * abs(w)))
    m = snap.m / units.mass
    pos = pos / units.length
    vel = vel / units.velocity
    with open(dynfile, 'w') as f:
        f.write('(Particle\n  N = %d\n(Log\n  made from the HENDOS file %s\n)Log\n'
                % (snap.n, os.path.basename(icfile)))
        f.write('(Dynamics\n  system_time  =  0\n  m  =  1\n'
                '  r  =  0 0 0\n  v  =  0 0 0\n)Dynamics\n(Hydro\n)Hydro\n(Star\n)Star\n')
        for i in range(snap.n):
            f.write('(Particle\n  i = %d\n  N = 1\n(Log\n)Log\n(Dynamics\n' % (i + 1))
            f.write('  m  =  %.17g\n' % m[i])
            f.write('  r  =  %.17g %.17g %.17g\n' % tuple(pos[i]))
            f.write('  v  =  %.17g %.17g %.17g\n' % tuple(vel[i]))
            f.write(')Dynamics\n(Hydro\n)Hydro\n(Star\n)Star\n)Particle\n')
        f.write(')Particle\n')
    return units, snap


class DynSnapshot:
    """One snapshot from a dyn file: single stars in the frame of the root.

    time, n, m, pos (n by 3), vel (n by 3), id (the 'i' label, 0 if none),
    in the units of the file (N-body units for kira output).
    to_physical(units) gives a copy in Msun, pc, km/s, Myr.
    """

    def __init__(self, time, ident, m, pos, vel):
        self.time = time
        self.id = np.array(ident)
        self.m = np.array(m)
        self.pos = np.array(pos).reshape(-1, 3)
        self.vel = np.array(vel).reshape(-1, 3)
        self.n = len(self.m)
        self.phi = None
        self.dt = None
        self.type = np.full(self.n, 'star')
        self.header = {}

    def value(self, key):
        return None

    def to_physical(self, units):
        return DynSnapshot(self.time * units.time, self.id, self.m * units.mass,
                           self.pos * units.length, self.vel * units.velocity)


def read_dyn(path):
    """Read every snapshot in a Starlab dyn file. Returns a list of DynSnapshot.

    Stars come out in the order of their 'i' labels where they all have
    one, so they line up with the HENDOS file they were made from.
    """
    snapshots = []
    stack = []            # the Particles we are inside: dictionaries
    story = None          # which story we are in, inside a Particle
    with open(path) as f:
        for raw in f:
            line = raw.strip()
            if not line:
                continue
            if line.startswith('(Particle'):
                stack.append(dict(m=0.0, r=np.zeros(3), v=np.zeros(3), i=0, time=None,
                                  children=[]))
                story = None
            elif line.startswith(')Particle'):
                node = stack.pop()
                if stack:
                    stack[-1]['children'].append(node)
                else:
                    snapshots.append(_flatten(node))
            elif line.startswith('(') and stack:
                story = line[1:]
            elif line.startswith(')') and stack:
                story = None
            elif stack and '=' in line and story in (None, 'Dynamics'):
                key, value = [x.strip() for x in line.split('=', 1)]
                node = stack[-1]
                try:
                    if story is None and key == 'i':
                        node['i'] = int(value.split()[0])
                    elif story == 'Dynamics' and key == 'm':
                        node['m'] = float(value.split()[0])
                    elif story == 'Dynamics' and key == 'r':
                        node['r'] = np.array([float(x) for x in value.split()[:3]])
                    elif story == 'Dynamics' and key == 'v':
                        node['v'] = np.array([float(x) for x in value.split()[:3]])
                    elif story == 'Dynamics' and key == 'real_system_time':
                        node['time'] = float(value.split()[0])
                    elif story == 'Dynamics' and key in ('system_time', 't') \
                            and node['time'] is None:
                        node['time'] = float(value.split()[0])
                except ValueError:
                    pass          # e.g. a time written in Starlab's extended format
    return snapshots


def _flatten(root):
    """The single stars under a root Particle, with absolute positions and velocities."""
    ident, m, pos, vel = [], [], [], []

    def walk(node, r0, v0):
        r = r0 + node['r']
        v = v0 + node['v']
        if node['children']:
            for child in node['children']:
                walk(child, r, v)
        else:
            ident.append(node['i'])
            m.append(node['m'])
            pos.append(r)
            vel.append(v)

    walk(root, np.zeros(3), np.zeros(3))
    snap = DynSnapshot(root['time'] if root['time'] is not None else 0.0, ident, m, pos, vel)
    if len(set(ident)) == len(ident) and min(ident) > 0:
        order = np.argsort(snap.id)
        snap = DynSnapshot(snap.time, snap.id[order], snap.m[order], snap.pos[order],
                           snap.vel[order])
    return snap


if __name__ == '__main__':
    if len(sys.argv) != 3:
        raise SystemExit('usage: python starlabio.py <HENDOS ic file> <dyn file to write>')
    u, s = hendos_to_dyn(sys.argv[1], sys.argv[2])
    u.save(sys.argv[2] + '.units.json')
    print('%d stars written to %s' % (s.n, sys.argv[2]))
    print('one N-body unit is  %.6g Msun   %.6g pc   %.6g km/s   %.6g Myr'
          % (u.mass, u.length, u.velocity, u.time))
