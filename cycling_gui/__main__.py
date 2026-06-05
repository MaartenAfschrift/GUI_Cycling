from __future__ import annotations

import argparse
import sys

from PyQt5.QtWidgets import QApplication

from .config import CyclingSettings, ObjectiveSettings, Scaling, SolverSettings
from .gui import CyclingMainWindow
from .simulation import CyclingSimulator


def _build_settings_from_args(args: argparse.Namespace) -> CyclingSettings:
    return CyclingSettings(
        name=args.name,
        muscle_model="Leuven",
        msk_model="Kistemaker",
        specific_tension=25.0,
        crank_frequency_rpm=args.freq,
        is_isokinetic=True,
        saddle_rx=args.saddle_rx,
        saddle_ry=args.saddle_ry,
        crank_length=args.crank_length,
        fixed_power=None,
        opt_frequency=False,
        collocation_points=args.mesh,
        metabolism_model="Bhargava2004",
        metabolism_scale_rate=1.0,
        bool_plot=False,
        scaling=Scaling(),
        objective=ObjectiveSettings(type="maxpower"),
        solver=SolverSettings(max_iter=args.max_iter, tol=args.tol),
    )


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description="Cycling GUI / simulation")
    parser.add_argument("--headless", action="store_true", help="Run a simulation and exit")
    parser.add_argument("--name", default="MaartenAfschrift", help="Result name")
    parser.add_argument("--freq", type=float, default=90.0, help="Pedalling frequency [rot/min]")
    parser.add_argument("--crank-length", type=float, default=0.175, help="Crank length [m]")
    parser.add_argument("--saddle-rx", type=float, default=-0.05, help="Saddle x position [m]")
    parser.add_argument("--saddle-ry", type=float, default=0.87, help="Saddle y position [m]")
    parser.add_argument("--mesh", type=int, default=50, help="Number of collocation points")
    parser.add_argument("--max-iter", type=int, default=1000, help="IPOPT maximum iterations")
    parser.add_argument("--tol", type=float, default=1.0e-4, help="IPOPT tolerance")

    args = parser.parse_args(argv)

    if args.headless:
        settings = _build_settings_from_args(args)
        simulator = CyclingSimulator(settings)
        result, json_path = simulator.simulate_and_export(name=settings.name)
        print(f"Average power: {result.average_power:.2f} W")
        print(f"Objective: {result.objective_value:.6f}")
        print(f"JSON saved: {json_path}")
        return 0

    app = QApplication.instance()
    if app is None:
        app = QApplication(sys.argv)
    window = CyclingMainWindow()
    window.show()
    return app.exec_()


if __name__ == "__main__":
    raise SystemExit(main())


