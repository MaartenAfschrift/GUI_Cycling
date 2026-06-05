from __future__ import annotations

from dataclasses import asdict, dataclass, field
from typing import Any


@dataclass(slots=True)
class Scaling:
    fi: float = 1.0
    fid: float = 10.0
    fidd: float = 1000.0
    Fr5x: float = 100.0
    Fr5y: float = 400.0
    M1: float = 100.0


@dataclass(slots=True)
class ObjectiveSettings:
    type: str = "maxpower"
    w_metab: float = 0.1
    w_stim: float = 100.0
    w_qdd: float = 0.001
    w_vMtilde: float = 1.0e-4
    w_M1: float = 1.0e-4
    scale: float = 1.0e-2
    C_forces: float = 1.0e-5
    w_cranckP: float = 0.1
    w_minNegWork: float = 10.0


@dataclass(slots=True)
class SolverSettings:
    solver: str = "ipopt"
    derivativelevel: str = "second"
    tol: float = 1.0e-4
    max_iter: int = 1000
    linear_solver: str = "mumps"
    nlp_scaling_method: str = "none"
    expand: bool = True

    def as_casadi_options(self) -> dict[str, Any]:
        return {
            "ipopt": {
                "tol": self.tol,
                "max_iter": self.max_iter,
                "linear_solver": self.linear_solver,
                "nlp_scaling_method": self.nlp_scaling_method,
            },
            "expand": self.expand,
            "print_time": False,
        }


@dataclass(slots=True)
class CyclingSettings:
    name: str = "MaartenAfschrift"
    muscle_model: str = "Leuven"
    msk_model: str = "Kistemaker"
    specific_tension: float = 25.0
    crank_frequency_rpm: float = 90.0
    is_isokinetic: bool = True
    saddle_rx: float = -0.05
    saddle_ry: float = 0.87
    crank_length: float = 0.175
    fixed_power: float | None = None
    opt_frequency: bool = False
    collocation_points: int = 50
    metabolism_model: str = "Bhargava2004"
    metabolism_scale_rate: float = 1.0
    bool_plot: bool = False
    scaling: Scaling = field(default_factory=Scaling)
    objective: ObjectiveSettings = field(default_factory=ObjectiveSettings)
    solver: SolverSettings = field(default_factory=SolverSettings)

    @property
    def crank_frequency_hz(self) -> float:
        return self.crank_frequency_rpm / 60.0

    @property
    def saddle(self) -> tuple[float, float]:
        return (self.saddle_rx, self.saddle_ry)

    def to_dict(self) -> dict[str, Any]:
        return asdict(self)

