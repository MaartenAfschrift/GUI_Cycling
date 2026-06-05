from __future__ import annotations

import numpy as np


def trapezoidal_integrator(x, x1, xd, xd1, dt):
    return (x1 - x) - (0.5 * dt * (xd + xd1))


def get_joint_positions(seg_lengths, fi):
    seg_lengths = np.asarray(seg_lengths, dtype=float).reshape(-1)
    fi = np.asarray(fi, dtype=float)

    nseg, nframes = fi.shape
    if nseg != 4 and nframes == 4:
        fi = fi.T
        nseg, nframes = fi.shape

    r_joint = np.zeros((2, nseg + 1, nframes), dtype=float)
    for frame in range(nframes):
        for k in range(nseg):
            r_prox = r_joint[:, k, frame]
            r_joint[:, k + 1, frame] = r_prox + seg_lengths[k] * np.array(
                [np.cos(fi[k, frame]), np.sin(fi[k, frame])],
                dtype=float,
            )
    return r_joint


