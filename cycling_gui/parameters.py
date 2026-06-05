from __future__ import annotations

from dataclasses import dataclass
import json
from pathlib import Path
from typing import Any

import numpy as np

PACKAGE_DIR = Path(__file__).resolve().parent
DEFAULT_KVS_JSON = PACKAGE_DIR / "data" / "kvs_parms.json"


@dataclass(slots=True)
class MuscleParameters:
    A0: np.ndarray
    A1: np.ndarray
    A2: np.ndarray
    fmax: np.ndarray
    lce_opt: np.ndarray
    lse0: np.ndarray
    kse: np.ndarray
    q0: float
    rm: float
    gamma_0: float
    kCa: float
    a_act: float
    a1_act: float
    b_act: np.ndarray
    arel_c: float = 0.41
    brel_c: float = 5.2
    fasymp: float = 1.5
    slopfac: float = 2.0
    vfactmin: float = 0.1
    width: float = 0.56
    C_fl: float = -1.0 / (0.56**2)
    b_Fse: float = 1.0
    sloplin: float = 0.0
    nmus: int = 9
    muscle_names: tuple[str, ...] = ()


@dataclass(slots=True)
class SegmentParameters:
    L: np.ndarray
    d: np.ndarray
    m: np.ndarray
    j: np.ndarray
    p: np.ndarray
    g: float = 9.81


@dataclass(slots=True)
class RawKvsData:
    data: dict[str, Any]


def _as_float_array(value: Any) -> np.ndarray:
    return np.asarray(value, dtype=float)


def load_raw_kvs_data(source: str | Path | None = None) -> RawKvsData:
    path = Path(source) if source is not None else DEFAULT_KVS_JSON
    with path.open("r", encoding="utf-8") as fh:
        data = json.load(fh)
    return RawKvsData(data=data)


def build_muscle_parameters(source: str | Path | None = None) -> MuscleParameters:
    raw = load_raw_kvs_data(source).data

    A0 = _as_float_array(raw["A0"])
    A1 = _as_float_array(raw["A1"])
    A2 = _as_float_array(raw["A2"])
    ce_F_max = _as_float_array(raw["ce_F_max"])
    ce_len_opt = _as_float_array(raw["ce_len_opt"])
    se_len_slack = _as_float_array(raw["se_len_slack"])
    hafo_q0 = _as_float_array(raw["hafo_q0"])
    hafo_m = _as_float_array(raw["hafo_m"])

    # MATLAB port of LoadParameters_hamsplit_DAK.m
    A0 = np.vstack([A0[:, 2:5], np.array([0.0, 0.3388, 0.0])])
    A1 = np.vstack([A1[:, 2:5], np.array([0.0, -0.0260, 0.0])])
    A2 = np.vstack([A2[:, 2:5], np.array([0.0, 0.0, 0.0])])

    fmax = np.concatenate([ce_F_max, np.array([ce_F_max[6] * (5.0 / 40.0)])])
    fmax[6] *= 35.0 / 40.0

    lce_opt = np.concatenate([ce_len_opt, np.array([0.11])])
    lse0 = np.concatenate([se_len_slack, np.array([0.20])])
    kse = fmax / ((lse0 * 0.04) ** 2)

    q0 = float(hafo_q0[0])
    rm = float(hafo_m[0])
    gamma_0 = 1.0e-5
    kCa = 0.8e-5
    a_act = -4.587
    a1_act = float(np.log10(np.exp(a_act)))
    b_act = np.array([5.168, 1.081, -0.1909], dtype=float)

    muscle_names = tuple(str(x) for x in raw.get("muscle_names", [])) + ("9 biceps femoris 2",)

    return MuscleParameters(
        A0=A0,
        A1=A1,
        A2=A2,
        fmax=fmax,
        lce_opt=lce_opt,
        lse0=lse0,
        kse=kse,
        q0=q0,
        rm=rm,
        gamma_0=gamma_0,
        kCa=kCa,
        a_act=a_act,
        a1_act=a1_act,
        b_act=b_act,
        sloplin=float((0.1 * 5.2) / (2.0 * 0.005 * 0.0975 * (1.0 + 0.41))),
        muscle_names=muscle_names,
    )


def build_segment_parameters(crank_length: float) -> SegmentParameters:
    L = np.array([crank_length, 0.165, 0.458, 0.4851], dtype=float)
    d = np.array([crank_length / 2.0, 0.120, 0.260, 0.275], dtype=float)
    m = np.array([0.200, 1.234, 3.540, 8.470], dtype=float)
    j = np.array([0.001, 0.010, 0.068, 0.209], dtype=float)
    p = L - d
    return SegmentParameters(L=L, d=d, m=m, j=j, p=p)

