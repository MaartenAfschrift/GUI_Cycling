from __future__ import annotations

from typing import Any

import casadi as ca
import numpy as np


def _is_casadi_value(value: Any) -> bool:
    return isinstance(value, (ca.SX, ca.MX, ca.DM))


def _xp(*values: Any):
    return ca if any(_is_casadi_value(v) for v in values) else np


def bhargava2004_energy(
    exc,
    act,
    lMtilde,
    vM,
    Fce,
    Fpass,
    musclemass,
    pctst,
    Fiso,
    Fmax,
    modelmass,
    b,
    scaleRate=0.0,
):
    xp = _xp(exc, act, lMtilde, vM, Fce, Fpass)

    pctft = 1.0 - pctst

    st_e = pctst * xp.sin(np.pi / 2.0 * exc)
    ft_e = pctft * (1.0 - xp.cos(np.pi / 2.0 * exc))

    decay_function_value = 1.0
    activation_constant_st = 40.0
    activation_constant_ft = 133.0
    Adot = musclemass * decay_function_value * (
        (activation_constant_st * st_e) + (activation_constant_ft * ft_e)
    )

    fiber_length_dep = lMtilde
    maintenance_constant_st = 74.0
    maintenance_constant_ft = 111.0
    Mdot = musclemass * fiber_length_dep * (
        (maintenance_constant_st * st_e) + (maintenance_constant_ft * ft_e)
    )

    F_iso = act * Fiso * Fmax
    fiber_force_total = Fce + Fpass
    alpha = (0.16 * F_iso) + (0.18 * fiber_force_total)
    vM_pos = 0.5 + 0.5 * xp.tanh(b * vM)
    vM_neg = 1.0 - vM_pos
    alpha = alpha + (-alpha + 0.157 * fiber_force_total) * vM_pos
    Sdot = -scaleRate * alpha * vM

    Wdot = -Fce * vM * vM_neg

    Edot_W_beforeClamp = Adot + Mdot + Sdot + Wdot
    Edot_Wkg_beforeClamp_neg = 0.5 + (0.5 * xp.tanh(b * (-Edot_W_beforeClamp)))
    Sdot = Sdot - Edot_W_beforeClamp * Edot_Wkg_beforeClamp_neg

    totalHeatRate = Adot + Mdot + Sdot
    totalHeatRate = totalHeatRate / musclemass
    totalHeatRate = totalHeatRate + (-totalHeatRate + 1.0) * (0.5 + 0.5 * xp.tanh(b * (1.0 - totalHeatRate)))
    totalHeatRate = totalHeatRate * musclemass

    energy_total = totalHeatRate + Wdot

    basal_coef = 1.2
    basal_exp = 1.0
    if _is_casadi_value(energy_total):
        energy_model = basal_coef * (modelmass**basal_exp) + ca.sum1(energy_total)
    else:
        energy_model = basal_coef * (modelmass**basal_exp) + np.sum(energy_total)

    return energy_total, Adot, Mdot, Sdot, Wdot, energy_model

