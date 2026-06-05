from __future__ import annotations

from typing import Any

import casadi as ca
import numpy as np


def _is_casadi_value(value: Any) -> bool:
    return isinstance(value, (ca.SX, ca.MX, ca.DM))


def _xp(*values: Any):
    return ca if any(_is_casadi_value(v) for v in values) else np


def muscle_force_length_vu(lcerel):
    xp = _xp(lcerel)
    fmin = 1e-2
    c1_gaus = 0.407742573856005
    c2_gaus = 0.447996880165858
    return ((1 - fmin) / 2.0) * (
        xp.exp(-((lcerel - 1.0) / c1_gaus) ** 2)
        + xp.exp(-((lcerel - 1.0) / c2_gaus) ** 4)
    ) + fmin


def muscle_force_velocity_vu(
    fcerel,
    q,
    lcerel,
    parms,
    bool_activation_dependent: bool = True,
    bool_length_dependent: bool = True,
):
    xp = _xp(fcerel, q, lcerel)

    def sigma(x, x0, w):
        return 1.0 / (1.0 + xp.exp(-(x - x0) / w))

    if bool_length_dependent:
        fisomrel = muscle_force_length_vu(lcerel)
    else:
        fisomrel = 1.0 + 0.0 * fcerel
        lcerel = 1.0 + 0.0 * fcerel

    q0_b = (xp.log(1.0 / parms.vfactmin - 1.0) + parms.q0 * 22.0) / 22.0
    if bool_activation_dependent:
        brel = parms.brel_c / (1.0 + xp.exp(-22.0 * (q - q0_b)))
    else:
        brel = parms.brel_c

    arel = parms.arel_c * fisomrel

    dvdf_isom_con = brel / (q * (fisomrel + arel))
    dvdf_isom_ecc = dvdf_isom_con / parms.slopfac
    dFdvcon0 = 1.0 / dvdf_isom_con
    s_as = 1.0 / parms.sloplin
    p1 = -(fisomrel * q * (parms.fasymp - 1.0)) / (s_as - dFdvcon0 * parms.slopfac)
    p3 = -parms.fasymp * fisomrel * q
    p2 = (fisomrel**2 * q**2 * (parms.fasymp - 1.0) ** 2) / (s_as - dFdvcon0 * parms.slopfac)
    p4 = -s_as

    sig_c1 = sigma(fcerel, fisomrel * q, 0.01)
    lcereld_c = (1.0 - sig_c1) * brel * (fcerel - q * fisomrel) / (fcerel + q * arel)
    sig_c2 = sigma(dvdf_isom_con, parms.sloplin, 0.01)
    lcereld_c = lcereld_c + sig_c2 * parms.sloplin * (fcerel - q * fisomrel)

    sig_e1 = sigma(fcerel, fisomrel * q, 0.01)
    sqrt_term = xp.sqrt(
        fcerel**2
        - 2.0 * fcerel * p1 * p4
        + 2.0 * fcerel * p3
        + p1**2 * p4**2
        - 2.0 * p1 * p3 * p4
        + p3**2
        + 4.0 * p2 * p4
    )
    lcereld_e = sig_e1 * (-(fcerel + p3 + p1 * p4 + sqrt_term) / (2.0 * p4))
    sig_e2 = sigma(dvdf_isom_ecc, parms.sloplin / parms.slopfac, 0.01)
    lcereld_e = lcereld_e + sig_e2 * (parms.sloplin / parms.slopfac) * (fcerel - q * fisomrel)

    return lcereld_c + lcereld_e


def get_force_length_velocity_properties(lMtilde, vMtilde, vMtildemax):
    xp = _xp(lMtilde, vMtilde)

    Fvparam = np.array([-0.318323436899128, -8.14915604347525, -0.37412150864786, 0.885644059915004], dtype=float)
    Faparam = np.array(
        [
            0.814483478343008,
            1.05503342897057,
            0.162384573599574,
            0.0633034484654646,
            0.433004984392647,
            0.71677541339776,
            -0.0299471169706956,
            0.200356847296188,
        ],
        dtype=float,
    )
    Fpparam = np.array([-0.995172050006169, 53.598150033144236], dtype=float)

    e0 = 0.6
    kpe = 4.0
    t5 = xp.exp(kpe * (lMtilde - 1.0) / e0)
    Fpe = ((t5 - 1.0) - Fpparam[0]) / Fpparam[1]

    b11, b21, b31, b41, b12, b22, b32, b42 = Faparam
    b13 = 0.1
    b23 = 1.0
    b33 = 0.5 * np.sqrt(0.5)
    b43 = 0.0

    num3 = lMtilde - b23
    den3 = b33 + b43 * lMtilde
    FMtilde3 = b13 * xp.exp(-0.5 * num3**2 / den3**2)

    num1 = lMtilde - b21
    den1 = b31 + b41 * lMtilde
    FMtilde1 = b11 * xp.exp(-0.5 * num1**2 / den1**2)

    num2 = lMtilde - b22
    den2 = b32 + b42 * lMtilde
    FMtilde2 = b12 * xp.exp(-0.5 * num2**2 / den2**2)

    FMltilde = FMtilde1 + FMtilde2 + FMtilde3

    e1, e2, e3, e4 = Fvparam
    FMvtilde = e1 * xp.log((e2 * vMtilde / vMtildemax + e3) + xp.sqrt((e2 * vMtilde / vMtildemax + e3) ** 2 + 1.0)) + e4
    return Fpe, FMltilde, FMvtilde


