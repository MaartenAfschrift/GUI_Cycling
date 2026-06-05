from __future__ import annotations

from dataclasses import asdict, dataclass
from typing import Any

import casadi as ca
import numpy as np

from .config import CyclingSettings
from .dynamics import build_rigid_body_model
from .io import export_result_json
from .math_utils import get_joint_positions, trapezoidal_integrator
from .metabolic import bhargava2004_energy
from .muscle import get_force_length_velocity_properties
from .parameters import MuscleParameters, build_muscle_parameters


@dataclass(slots=True)
class SimulationResult:
    fi: np.ndarray
    fid: np.ndarray
    fidd: np.ndarray
    stim: np.ndarray
    Fr5x: np.ndarray
    Fr5y: np.ndarray
    M1: np.ndarray
    t: np.ndarray
    r_joint: np.ndarray
    crank_power: np.ndarray
    average_power: float
    objective_value: float
    stats: dict[str, Any]
    settings: dict[str, Any]
    muscle_params: dict[str, Any] | None = None
    a: np.ndarray | None = None
    lMtilde: np.ndarray | None = None
    vMtilde: np.ndarray | None = None
    Fce: np.ndarray | None = None
    FM: np.ndarray | None = None
    TMus: np.ndarray | None = None
    fse: np.ndarray | None = None
    energy_total: np.ndarray | None = None
    solver_info: dict[str, Any] | None = None

    def summary(self) -> dict[str, Any]:
        summary = {
            "name": self.settings.get("name", "simulation"),
            "average_power": self.average_power,
            "objective_value": self.objective_value,
            "crank_length": self.settings.get("crank_length"),
            "saddle_rx": self.settings.get("saddle_rx"),
            "saddle_ry": self.settings.get("saddle_ry"),
            "freq_rpm": self.settings.get("crank_frequency_rpm"),
        }
        if self.solver_info:
            summary["solver_info"] = self.solver_info
        return summary

    def to_dict(self) -> dict[str, Any]:
        payload = asdict(self)
        payload["summary"] = self.summary()
        return payload


class CyclingSimulator:
    def __init__(self, settings: CyclingSettings | None = None, *, muscle_source: str | None = None):
        self.settings = settings or CyclingSettings()
        self.muscle_params: MuscleParameters = build_muscle_parameters(muscle_source)
        self.rbd = build_rigid_body_model(self.settings.crank_length)

    def _solver_options(self) -> dict[str, Any]:
        return self.settings.solver.as_casadi_options()

    def _cycle_setup(self) -> tuple[int, float, float, float, np.ndarray]:
        s = self.settings
        N = int(s.collocation_points)
        cf = float(s.crank_frequency_hz)
        v_crank = -2.0 * np.pi * cf
        time = 1.0 / cf
        h = time / (N - 1)
        t_vect = np.linspace(0.0, time, N)
        return N, v_crank, time, h, t_vect

    def _muscle_geometry(self, fi_col):
        m = self.muscle_params
        fijo0 = fi_col[2] - fi_col[1]
        fijo1 = fi_col[3] - fi_col[2]
        fijo2 = 1.1 - fi_col[3]

        lmtc = (
            m.A0[:, 0] + m.A1[:, 0] * fijo0 + m.A2[:, 0] * fijo0**2
            + m.A0[:, 1] + m.A1[:, 1] * fijo1 + m.A2[:, 1] * fijo1**2
            + m.A0[:, 2] + m.A1[:, 2] * fijo2 + m.A2[:, 2] * fijo2**2
        )
        lmtc = ca.reshape(lmtc, m.nmus, 1)

        momarm = -ca.hcat(
            [
                m.A1[:, 0] + 2.0 * m.A2[:, 0] * fijo0,
                m.A1[:, 1] + 2.0 * m.A2[:, 1] * fijo1,
                m.A1[:, 2] + 2.0 * m.A2[:, 2] * fijo2,
            ]
        )
        return lmtc, momarm

    def _tendon_force_leuven(self, lcerel_col, fi_col):
        m = self.muscle_params
        lmtc, momarm = self._muscle_geometry(fi_col)
        lse = lmtc - (lcerel_col * ca.DM(m.lce_opt).reshape((m.nmus, 1)))
        lTtilde = lse / ca.DM(m.lse0).reshape((m.nmus, 1))
        fse = ca.exp(35.0 * (lTtilde - 0.995)) / 5.0 - 0.25
        fse_N = fse * ca.DM(m.fmax).reshape((m.nmus, 1))
        return fse_N, momarm

    def _solve_torque_stage(self) -> dict[str, Any]:
        s = self.settings
        model = self.rbd
        N, v_crank, _time, h, _ = self._cycle_setup()
        nseg = 4
        Topt = 300.0

        opti = ca.Opti()
        Sfi = opti.variable(nseg, N)
        Sfid = opti.variable(nseg, N)
        Sfidd = opti.variable(nseg, N)
        stim = opti.variable(3, N)
        SFr5x = opti.variable(1, N)
        SFr5y = opti.variable(1, N)
        SM1 = opti.variable(1, N)

        fi = Sfi * s.scaling.fi
        fid = Sfid * s.scaling.fid
        fidd = Sfidd * s.scaling.fidd
        Fr5x = SFr5x * s.scaling.Fr5x
        Fr5y = SFr5y * s.scaling.Fr5y
        M1 = SM1 * s.scaling.M1

        opti.subject_to(opti.bounded(np.pi / 2.0 + 0.3, fi[1, :], np.pi + 0.3))
        opti.subject_to(opti.bounded(-np.pi, fi[2, :] - fi[1, :], 0.0))
        opti.subject_to(opti.bounded(0.0, fi[3, :] - fi[2, :], np.pi))
        opti.subject_to(opti.bounded(-np.pi, 1.1 - fi[3, :], 0.0))
        opti.subject_to(opti.bounded(-1.0, stim, 1.0))

        fi_guess = np.array([1.5708, 0.9117, 0.9117, 2.6898], dtype=float)
        fid_guess = np.array([v_crank, 1.2534, 3.4792, 2.5262], dtype=float)
        opti.set_initial(Sfi, np.tile(fi_guess.reshape(-1, 1), (1, N)) / s.scaling.fi)
        opti.set_initial(Sfid, np.tile(fid_guess.reshape(-1, 1), (1, N)) / s.scaling.fid)
        opti.set_initial(Sfidd, 0.0)
        opti.set_initial(SFr5x, -200.0 / s.scaling.Fr5x)
        opti.set_initial(SFr5y, -200.0 / s.scaling.Fr5y)
        opti.set_initial(SM1, -50.0 / s.scaling.M1)

        opti.subject_to(fi[0, 0] == 0.5 * np.pi)
        opti.subject_to(fid[0, :] == v_crank)

        x = ca.vertcat(fi[:, : N - 1], fid[:, : N - 1])
        xt1 = ca.vertcat(fi[:, 1:N], fid[:, 1:N])
        xd = ca.vertcat(fid[:, : N - 1], fidd[:, : N - 1])
        xdt1 = ca.vertcat(fid[:, 1:N], fidd[:, 1:N])
        opti.subject_to(trapezoidal_integrator(x, xt1, xd, xdt1, h) == 0)

        Terr = ca.MX.zeros(nseg, N)
        rhc = ca.MX.zeros(2, N)
        for k in range(N):
            Tm = stim[:, k] * Topt
            Terr[:, k] = model.skeldyn(
                fi[:, k],
                fid[:, k],
                fidd[:, k],
                Fr5x[0, k],
                Fr5y[0, k],
                M1[0, k],
                Tm,
            )
            rhc[:, k], _ = model.hip_kin(fi[:, k], fid[:, k])

        opti.subject_to(Terr == 0)
        opti.subject_to(rhc == np.tile(np.array(s.saddle, dtype=float).reshape(2, 1), (1, N)))
        opti.subject_to(fi[1, 0] == fi[1, N - 1])
        opti.subject_to(fid[1, 0] == fid[1, N - 1])

        average_power = ca.sum1(ca.vec(-fid[0, :] * M1)) / N
        if s.fixed_power is not None:
            opti.subject_to(average_power == float(s.fixed_power))

        Jstim = s.objective.w_stim * ca.sumsqr(stim) / N / 3.0
        Jfidd = s.objective.w_qdd * ca.sumsqr(fidd) / N / nseg
        J = Jstim + Jfidd
        opti.minimize(J)

        opti.solver(s.solver.solver, self._solver_options())
        sol = opti.solve()
        return {
            "fi": np.asarray(sol.value(fi), dtype=float),
            "fid": np.asarray(sol.value(fid), dtype=float),
            "fidd": np.asarray(sol.value(fidd), dtype=float),
            "Fr5x": np.asarray(sol.value(Fr5x), dtype=float),
            "Fr5y": np.asarray(sol.value(Fr5y), dtype=float),
            "M1": np.asarray(sol.value(M1), dtype=float),
            "stats": dict(sol.stats()),
        }

    def _solve_muscle_stage(self, torque_guess: dict[str, Any]) -> dict[str, Any]:
        s = self.settings
        m = self.muscle_params
        model = self.rbd
        if s.muscle_model != "Leuven":
            raise ValueError("Only Leuven muscle model is implemented in the Python port.")

        N, v_crank, _time, h, _ = self._cycle_setup()
        nseg = 4
        nmus = m.nmus

        opti = ca.Opti()
        Sfi = opti.variable(nseg, N)
        Sfid = opti.variable(nseg, N)
        Sfidd = opti.variable(nseg, N)
        a = opti.variable(nmus, N)
        lcerel = opti.variable(nmus, N)
        stim = opti.variable(nmus, N)
        SFr5x = opti.variable(1, N)
        SFr5y = opti.variable(1, N)
        SM1 = opti.variable(1, N)
        vcerel_helper = opti.variable(nmus, N)

        fi = Sfi * s.scaling.fi
        fid = Sfid * s.scaling.fid
        fidd = Sfidd * s.scaling.fidd
        Fr5x = SFr5x * s.scaling.Fr5x
        Fr5y = SFr5y * s.scaling.Fr5y
        M1 = SM1 * s.scaling.M1

        opti.subject_to(opti.bounded(0.0, stim, 1.0))
        opti.subject_to(opti.bounded(0.0, a, 1.0))
        opti.subject_to(opti.bounded(0.1, lcerel, 1.7))
        opti.subject_to(opti.bounded(-10.0, vcerel_helper, 10.0))

        opti.set_initial(Sfi, torque_guess["fi"] / s.scaling.fi)
        opti.set_initial(Sfid, torque_guess["fid"] / s.scaling.fid)
        opti.set_initial(Sfidd, torque_guess["fidd"] / s.scaling.fidd)
        opti.set_initial(SFr5x, torque_guess["Fr5x"] / s.scaling.Fr5x)
        opti.set_initial(SFr5y, torque_guess["Fr5y"] / s.scaling.Fr5y)
        opti.set_initial(SM1, torque_guess["M1"] / s.scaling.M1)
        opti.set_initial(a, 1.0)
        opti.set_initial(lcerel, 1.0)
        opti.set_initial(stim, 1.0)
        opti.set_initial(vcerel_helper, 0.0)

        opti.subject_to(fi[0, 0] == 0.5 * np.pi)
        opti.subject_to(fid[0, :] == v_crank)

        adot = ca.MX.zeros(nmus, N)
        T_id_err = ca.MX.zeros(nseg, N)
        rhc = ca.MX.zeros(2, N)
        errFEq = ca.MX.zeros(nmus, N)
        fseM = ca.MX.zeros(nmus, N)
        TMusM = ca.MX.zeros(3, N)
        FceV = ca.MX.zeros(nmus, N)
        FpasV = ca.MX.zeros(nmus, N)
        FMltildeV = ca.MX.zeros(nmus, N)
        FMvtildeV = ca.MX.zeros(nmus, N)
        FMV = ca.MX.zeros(nmus, N)
        energy_total = ca.MX.zeros(nmus, N)

        specific_tension = s.specific_tension
        musclemass = (m.fmax * m.lce_opt) * 1059.7 / (specific_tension * 1e4)
        musclemass_dm = ca.DM(musclemass).reshape((nmus, 1))
        fmax_dm = ca.DM(m.fmax).reshape((nmus, 1))
        lce_opt_dm = ca.DM(m.lce_opt).reshape((nmus, 1))
        pctst = 0.5

        for k in range(N):
            fse_N, momarm = self._tendon_force_leuven(lcerel[:, k], fi[:, k])
            TMus = ca.transpose(ca.mtimes(ca.transpose(fse_N), momarm))
            fseM[:, k] = fse_N
            TMusM[:, k] = TMus

            T_id_err[:, k] = model.skeldyn(
                fi[:, k],
                fid[:, k],
                fidd[:, k],
                Fr5x[0, k],
                Fr5y[0, k],
                M1[0, k],
                TMus,
            )

            rhc[:, k], _ = model.hip_kin(fi[:, k], fid[:, k])
            adot[:, k] = (stim[:, k] - a[:, k]) / 0.03

            Fpe, FMltilde, FMvtilde = get_force_length_velocity_properties(lcerel[:, k], vcerel_helper[:, k], 10.0)
            Fce = a[:, k] * (FMltilde * FMvtilde)
            FM = Fce + Fpe
            errFEq[:, k] = FM - (fse_N / fmax_dm)

            vM = lce_opt_dm * vcerel_helper[:, k]
            Fm = Fce * fmax_dm
            Etot, *_ = bhargava2004_energy(
                stim[:, k],
                a[:, k],
                lcerel[:, k],
                vM,
                Fm,
                Fpe,
                musclemass_dm,
                pctst,
                FMltilde,
                fmax_dm,
                75.0,
                100.0,
                s.metabolism_scale_rate,
            )

            FceV[:, k] = Fce
            FpasV[:, k] = Fpe
            FMltildeV[:, k] = FMltilde
            FMvtildeV[:, k] = FMvtilde
            FMV[:, k] = FM
            energy_total[:, k] = Etot

        opti.subject_to(errFEq == 0)
        opti.subject_to(rhc == np.tile(np.array(s.saddle, dtype=float).reshape(2, 1), (1, N)))
        opti.subject_to(T_id_err == 0)

        x = ca.vertcat(fi[:, : N - 1], fid[:, : N - 1], a[:, : N - 1], lcerel[:, : N - 1])
        xt1 = ca.vertcat(fi[:, 1:N], fid[:, 1:N], a[:, 1:N], lcerel[:, 1:N])
        xd = ca.vertcat(fid[:, : N - 1], fidd[:, : N - 1], adot[:, : N - 1], vcerel_helper[:, : N - 1])
        xdt1 = ca.vertcat(fid[:, 1:N], fidd[:, 1:N], adot[:, 1:N], vcerel_helper[:, 1:N])
        opti.subject_to(trapezoidal_integrator(x, xt1, xd, xdt1, h) == 0)

        crank_power = -fid[0, :] * M1
        average_power = ca.sum1(ca.vec(crank_power)) / N
        if (s.fixed_power is not None) and (s.objective.type != "maxpower"):
            opti.subject_to(average_power == float(s.fixed_power))

        opti.subject_to(fi[1, 0] == fi[1, N - 1])
        opti.subject_to(fid[1, 0] == fid[1, N - 1])
        opti.subject_to(lcerel[:, 0] == lcerel[:, N - 1])
        opti.subject_to(vcerel_helper[:, 0] == vcerel_helper[:, N - 1])

        Jmetab = s.objective.w_metab * ca.sum1(ca.vec(energy_total)) / N / nmus
        Jstim = s.objective.w_stim * ca.sumsqr(stim) / N / nmus
        Jfidd = s.objective.w_qdd * ca.sumsqr(fidd) / N / nseg
        Jvce = s.objective.w_vMtilde * ca.sumsqr(vcerel_helper) / N / nmus
        JM1 = s.objective.w_M1 * ca.sumsqr(M1) / N / nmus
        J_ConF = s.objective.C_forces * (ca.sumsqr(Fr5x) + ca.sumsqr(Fr5y)) / N

        if s.objective.type == "Multi_a_E":
            J = s.objective.scale * (Jmetab + Jstim + Jfidd + Jvce + JM1 + J_ConF)
        elif s.objective.type == "stim":
            J = s.objective.scale * (Jstim + Jvce + JM1 + J_ConF)
        elif s.objective.type == "maxpower":
            JCranckPower = s.objective.w_cranckP * (-average_power + 3000.0)
            J = s.objective.scale * (Jstim + Jvce + JM1 + J_ConF + JCranckPower)
        else:
            J = s.objective.scale * (Jstim + Jvce + JM1 + J_ConF)

        opti.minimize(J)
        opti.solver(s.solver.solver, self._solver_options())
        sol = opti.solve()

        return {
            "fi": np.asarray(sol.value(fi), dtype=float),
            "fid": np.asarray(sol.value(fid), dtype=float),
            "fidd": np.asarray(sol.value(fidd), dtype=float),
            "stim": np.asarray(sol.value(stim), dtype=float),
            "Fr5x": np.asarray(sol.value(Fr5x), dtype=float).reshape(-1),
            "Fr5y": np.asarray(sol.value(Fr5y), dtype=float).reshape(-1),
            "M1": np.asarray(sol.value(M1), dtype=float).reshape(-1),
            "a": np.asarray(sol.value(a), dtype=float),
            "lMtilde": np.asarray(sol.value(lcerel), dtype=float),
            "vMtilde": np.asarray(sol.value(vcerel_helper), dtype=float),
            "Fce": np.asarray(sol.value(FceV), dtype=float),
            "FM": np.asarray(sol.value(FMV), dtype=float),
            "TMus": np.asarray(sol.value(TMusM), dtype=float),
            "fse": np.asarray(sol.value(fseM), dtype=float),
            "energy_total": np.asarray(sol.value(energy_total), dtype=float),
            "objective": float(sol.value(J)),
            "average_power": float(sol.value(average_power)),
            "crank_power": np.asarray(sol.value(crank_power), dtype=float).reshape(-1),
            "stats": dict(sol.stats()),
        }

    def simulate(self) -> SimulationResult:
        s = self.settings
        m = self.muscle_params
        model = self.rbd
        N, _v_crank, _time, _h, t_vect = self._cycle_setup()
        torque_stage = self._solve_torque_stage()
        muscle_stage = self._solve_muscle_stage(torque_stage)

        fi_v = muscle_stage["fi"]
        fid_v = muscle_stage["fid"]
        fidd_v = muscle_stage["fidd"]
        stim_v = muscle_stage["stim"]
        Fr5x_v = muscle_stage["Fr5x"]
        Fr5y_v = muscle_stage["Fr5y"]
        M1_v = muscle_stage["M1"]
        crank_power = muscle_stage["crank_power"]
        average_power_value = muscle_stage["average_power"]
        objective_value = muscle_stage["objective"]
        r_joint = get_joint_positions(model.segparms.L, fi_v)

        stats = muscle_stage["stats"]
        solver_info = {
            "torque_stage": {
                "success": bool(torque_stage["stats"].get("success", False)),
                "return_status": torque_stage["stats"].get("return_status"),
                "iter_count": torque_stage["stats"].get("iter_count"),
            },
            "muscle_stage": {
                "success": bool(muscle_stage["stats"].get("success", False)),
                "return_status": muscle_stage["stats"].get("return_status"),
                "iter_count": muscle_stage["stats"].get("iter_count"),
            },
        }
        return SimulationResult(
            fi=fi_v,
            fid=fid_v,
            fidd=fidd_v,
            stim=stim_v,
            Fr5x=Fr5x_v,
            Fr5y=Fr5y_v,
            M1=M1_v,
            t=t_vect,
            r_joint=r_joint,
            crank_power=crank_power,
            average_power=average_power_value * 2,
            objective_value=objective_value,
            stats=stats,
            settings=s.to_dict(),
            muscle_params={
                "nmus": m.nmus,
                "fmax": m.fmax.tolist(),
                "lce_opt": m.lce_opt.tolist(),
                "muscle_names": list(m.muscle_names),
            },
            a=muscle_stage["a"],
            lMtilde=muscle_stage["lMtilde"],
            vMtilde=muscle_stage["vMtilde"],
            Fce=muscle_stage["Fce"],
            FM=muscle_stage["FM"],
            TMus=muscle_stage["TMus"],
            fse=muscle_stage["fse"],
            energy_total=muscle_stage["energy_total"],
            solver_info=solver_info,
        )

    def simulate_and_export(self, *, name: str | None = None, results_dir: str | None = None) -> tuple[SimulationResult, str]:
        result = self.simulate()
        export_name = name or self.settings.name
        path = export_result_json(result, name=export_name, results_dir=results_dir)
        return result, str(path)


