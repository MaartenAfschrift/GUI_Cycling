from __future__ import annotations

from dataclasses import dataclass
from typing import Any

import casadi as ca
import numpy as np

from .parameters import SegmentParameters


def _is_casadi_value(value: Any) -> bool:
    return isinstance(value, (ca.SX, ca.MX, ca.DM))


def _xp(*values: Any):
    return ca if any(_is_casadi_value(v) for v in values) else np


@dataclass(slots=True)
class RigidBodyModel:
    segparms: SegmentParameters

    def skeldyn(self, fi, fid, fidd, Fr5x, Fr5y, M, Tmus):
        xp = _xp(fi, fid, fidd, Fr5x, Fr5y, M, Tmus)
        fi1, fi2, fi3, fi4 = fi[0], fi[1], fi[2], fi[3]
        fid1, fid2, fid3, fid4 = fid[0], fid[1], fid[2], fid[3]
        fidd1, fidd2, fidd3, fidd4 = fidd[0], fidd[1], fidd[2], fidd[3]
        tor1, tor2, tor3 = Tmus[0], Tmus[1], Tmus[2]

        l1, l2, l3, l4 = self.segparms.L
        d1, d2, d3, d4 = self.segparms.d
        m1, m2, m3, m4 = self.segparms.m
        j1, j2, j3, j4 = self.segparms.j
        g = self.segparms.g

        c1, c2, c3, c4 = xp.cos(fi1), xp.cos(fi2), xp.cos(fi3), xp.cos(fi4)
        s1, s2, s3, s4 = xp.sin(fi1), xp.sin(fi2), xp.sin(fi3), xp.sin(fi4)

        A11 = -j1 - d1**2 * m1 - l1**2 * (m2 + m3 + m4)
        A12 = -l1 * (d2 * m2 + l2 * (m3 + m4)) * xp.cos(fi1 - fi2)
        A13 = -l1 * (d3 * m3 + l3 * m4) * xp.cos(fi1 - fi3)
        A14 = -l1 * d4 * m4 * xp.cos(fi1 - fi4)

        A21 = A12
        A22 = -j2 - d2**2 * m2 - l2**2 * (m3 + m4)
        A23 = -l2 * (d3 * m3 + l3 * m4) * xp.cos(fi2 - fi3)
        A24 = -l2 * d4 * m4 * xp.cos(fi2 - fi4)

        A31 = A13
        A32 = A23
        A33 = -j3 - d3**2 * m3 - l3**2 * m4
        A34 = -l3 * d4 * m4 * xp.cos(fi3 - fi4)

        A41 = A14
        A42 = A24
        A43 = A34
        A44 = -j4 - d4**2 * m4

        As = xp.vertcat(
            xp.horzcat(A11, A12, A13, A14),
            xp.horzcat(A21, A22, A23, A24),
            xp.horzcat(A31, A32, A33, A34),
            xp.horzcat(A41, A42, A43, A44),
        )

        bs1 = (l1 * s1) * (l1 * m2 * c1 * fid1**2 + d2 * m2 * c2 * fid2**2) - M
        bs1 -= (l1 * c1) * (l1 * m4 * s1 * fid1**2 + l2 * m4 * s2 * fid2**2 + l3 * m4 * s3 * fid3**2 + d4 * m4 * s4 * fid4**2 - Fr5y - g * m4)
        bs1 += (l1 * s1) * (l1 * m4 * c1 * fid1**2 + l2 * m4 * c2 * fid2**2 + l3 * m4 * c3 * fid3**2 + d4 * m4 * c4 * fid4**2 - Fr5x)
        bs1 -= (l1 * c1) * (l1 * m3 * s1 * fid1**2 + l2 * m3 * s2 * fid2**2 + d3 * m3 * s3 * fid3**2 - g * m3)
        bs1 += (l1 * s1) * (l1 * m3 * c1 * fid1**2 + l2 * m3 * c2 * fid2**2 + d3 * m3 * c3 * fid3**2)
        bs1 -= (l1 * c1) * (l1 * m2 * s1 * fid1**2 + d2 * m2 * s2 * fid2**2 - g * m2)
        bs1 += d1 * c1 * (-d1 * m1 * s1 * fid1**2 + g * m1) + d1**2 * fid1**2 * m1 * c1 * s1

        bs2 = tor1
        bs2 -= (l2 * c2) * (l1 * m4 * s1 * fid1**2 + l2 * m4 * s2 * fid2**2 + l3 * m4 * s3 * fid3**2 + d4 * m4 * s4 * fid4**2 - Fr5y - g * m4)
        bs2 += (l2 * s2) * (l1 * m4 * c1 * fid1**2 + l2 * m4 * c2 * fid2**2 + l3 * m4 * c3 * fid3**2 + d4 * m4 * c4 * fid4**2 - Fr5x)
        bs2 -= (l2 * c2) * (l1 * m3 * s1 * fid1**2 + l2 * m3 * s2 * fid2**2 + d3 * m3 * s3 * fid3**2 - g * m3)
        bs2 += (l2 * s2) * (l1 * m3 * c1 * fid1**2 + l2 * m3 * c2 * fid2**2 + d3 * m3 * c3 * fid3**2)
        bs2 -= d2 * c2 * (l1 * m2 * s1 * fid1**2 + d2 * m2 * s2 * fid2**2 - g * m2)
        bs2 += d2 * s2 * (l1 * m2 * c1 * fid1**2 + d2 * m2 * c2 * fid2**2)

        bs3 = tor2 - tor1
        bs3 -= (l3 * c3) * (l1 * m4 * s1 * fid1**2 + l2 * m4 * s2 * fid2**2 + l3 * m4 * s3 * fid3**2 + d4 * m4 * s4 * fid4**2 - Fr5y - g * m4)
        bs3 += (l3 * s3) * (l1 * m4 * c1 * fid1**2 + l2 * m4 * c2 * fid2**2 + l3 * m4 * c3 * fid3**2 + d4 * m4 * c4 * fid4**2 - Fr5x)
        bs3 -= d3 * c3 * (l1 * m3 * s1 * fid1**2 + l2 * m3 * s2 * fid2**2 + d3 * m3 * s3 * fid3**2 - g * m3)
        bs3 += d3 * s3 * (l1 * m3 * c1 * fid1**2 + l2 * m3 * c2 * fid2**2 + d3 * m3 * c3 * fid3**2)

        e4 = d4 - l4
        bs4 = tor3 - tor2
        bs4 += Fr5x * s4 * e4
        bs4 -= d4 * c4 * (l1 * m4 * s1 * fid1**2 + l2 * m4 * s2 * fid2**2 + l3 * m4 * s3 * fid3**2 + d4 * m4 * s4 * fid4**2 - Fr5y - g * m4)
        bs4 += d4 * s4 * (l1 * m4 * c1 * fid1**2 + l2 * m4 * c2 * fid2**2 + l3 * m4 * c3 * fid3**2 + d4 * m4 * c4 * fid4**2 - Fr5x)
        bs4 -= Fr5y * c4 * e4

        bs = xp.vertcat(bs1, bs2, bs3, bs4)
        return xp.mtimes(As, fidd) - bs

    def hip_kin(self, fi, fid):
        xp = _xp(fi, fid)
        l1, l2, l3, l4 = self.segparms.L
        if xp is ca:
            rh = ca.vertcat(0, 0)
            vh = ca.vertcat(0, 0)
        else:
            rh = np.zeros(2, dtype=float)
            vh = np.zeros(2, dtype=float)

        for k, Lk in enumerate((l1, l2, l3, l4)):
            rh = rh + Lk * xp.vertcat(xp.cos(fi[k]), xp.sin(fi[k]))
            vh = vh + Lk * fid[k] * xp.vertcat(-xp.sin(fi[k]), xp.cos(fi[k]))
        return rh, vh


def build_rigid_body_model(crank_length: float) -> RigidBodyModel:
    from .parameters import build_segment_parameters

    return RigidBodyModel(segparms=build_segment_parameters(crank_length))

