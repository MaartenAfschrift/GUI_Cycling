from __future__ import annotations

from pathlib import Path

import numpy as np
from PyQt5.QtCore import QTimer, Qt
from PyQt5.QtWidgets import (
    QApplication,
    QDoubleSpinBox,
    QFormLayout,
    QGridLayout,
    QGroupBox,
    QHBoxLayout,
    QLabel,
    QLineEdit,
    QMainWindow,
    QMessageBox,
    QPushButton,
    QSpinBox,
    QTextEdit,
    QVBoxLayout,
    QWidget,
)
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.figure import Figure

from .config import CyclingSettings, ObjectiveSettings, Scaling, SolverSettings
from .io import ensure_results_dir
from .simulation import CyclingSimulator, SimulationResult


class MplCanvas(FigureCanvas):
    def __init__(self, width=5.0, height=4.0, dpi=100):
        self.figure = Figure(figsize=(width, height), dpi=dpi)
        self.axes = self.figure.add_subplot(111)
        super().__init__(self.figure)


class CyclingMainWindow(QMainWindow):
    def __init__(self):
        super().__init__()
        self.setWindowTitle("Cycling GUI")
        self.resize(1200, 800)

        self.results_dir = ensure_results_dir(Path.cwd())
        self.simulator: CyclingSimulator | None = None
        self.result: SimulationResult | None = None
        self.history: list[float] = []
        self.history_names: list[str] = []
        self._animation_timer = QTimer(self)
        self._animation_timer.timeout.connect(self._advance_animation)
        self._animation_frame = 0
        self._animation_lines = []
        self._animation_labels = []
        self._animation_circle = None
        self._animation_saddle = None

        self._build_ui()
        self._apply_default_state()

    def _build_ui(self):
        central = QWidget(self)
        self.setCentralWidget(central)
        root = QHBoxLayout(central)

        self.visual_canvas = MplCanvas(width=6.5, height=6.0, dpi=100)
        self.visual_ax = self.visual_canvas.axes
        self.visual_ax.set_title("Cycling Model")
        self.visual_ax.set_xlabel("X [m]")
        self.visual_ax.set_ylabel("Y [m]")
        self.visual_ax.set_aspect("equal", adjustable="box")
        self.visual_ax.grid(True, alpha=0.3)
        self.visual_ax.set_xlim(-0.5, 0.5)
        self.visual_ax.set_ylim(-0.3, 1.1)

        left_box = QVBoxLayout()
        left_box.addWidget(self.visual_canvas)

        self.score_canvas = MplCanvas(width=5.0, height=2.8, dpi=100)
        self.score_ax = self.score_canvas.axes
        self.score_ax.set_title("Average mechanical power")
        self.score_ax.set_xlabel("attempts")
        self.score_ax.set_ylabel("power [W]")
        self.score_ax.grid(True, alpha=0.3)

        self.name_edit = QLineEdit()
        self.crank_length_spin = QDoubleSpinBox()
        self.crank_length_spin.setDecimals(3)
        self.crank_length_spin.setRange(0.01, 0.25)
        self.crank_length_spin.setSingleStep(0.005)

        self.saddle_rx_spin = QDoubleSpinBox()
        self.saddle_rx_spin.setDecimals(3)
        self.saddle_rx_spin.setRange(-0.2, 0.2)
        self.saddle_rx_spin.setSingleStep(0.005)

        self.saddle_ry_spin = QDoubleSpinBox()
        self.saddle_ry_spin.setDecimals(3)
        self.saddle_ry_spin.setRange(0.6, 1.1)
        self.saddle_ry_spin.setSingleStep(0.005)

        self.freq_spin = QDoubleSpinBox()
        self.freq_spin.setDecimals(1)
        self.freq_spin.setRange(50.0, 170.0)
        self.freq_spin.setSingleStep(1.0)

        self.objective_edit = QLineEdit()
        self.objective_edit.setReadOnly(True)

        self.ready_lamp = QLabel(" ")
        self.ready_lamp.setFixedSize(18, 18)
        self.ready_lamp.setStyleSheet("background-color: #999999; border-radius: 9px;")

        self.simulate_button = QPushButton("Simulate")
        self.simulate_button.clicked.connect(self.simulate)
        self.visualize_button = QPushButton("Visualise")
        self.visualize_button.clicked.connect(self.visualise)
        self.clear_button = QPushButton("Clear log")
        self.clear_button.clicked.connect(self.clear_log)

        controls_group = QGroupBox("Task settings")
        controls_layout = QGridLayout(controls_group)
        controls_layout.addWidget(QLabel("Name"), 0, 0)
        controls_layout.addWidget(self.name_edit, 0, 1)
        controls_layout.addWidget(QLabel("Crank length [m]"), 1, 0)
        controls_layout.addWidget(self.crank_length_spin, 1, 1)
        controls_layout.addWidget(QLabel("Saddle rx [m]"), 2, 0)
        controls_layout.addWidget(self.saddle_rx_spin, 2, 1)
        controls_layout.addWidget(QLabel("Saddle ry [m]"), 3, 0)
        controls_layout.addWidget(self.saddle_ry_spin, 3, 1)
        controls_layout.addWidget(QLabel("Freq [rot/min]"), 4, 0)
        controls_layout.addWidget(self.freq_spin, 4, 1)
        controls_layout.addWidget(QLabel("Objective"), 5, 0)
        controls_layout.addWidget(self.objective_edit, 5, 1)
        controls_layout.addWidget(self.simulate_button, 6, 0)
        controls_layout.addWidget(self.visualize_button, 6, 1)
        controls_layout.addWidget(self.clear_button, 7, 0)
        controls_layout.addWidget(QLabel("Ready?"), 7, 1)
        controls_layout.addWidget(self.ready_lamp, 7, 1, alignment=Qt.AlignRight)

        self.log_text = QTextEdit()
        self.log_text.setReadOnly(True)
        self.log_text.setMinimumHeight(170)
        self.log_text.setPlaceholderText("Simulation log...")

        right_box = QVBoxLayout()
        right_box.addWidget(self.score_canvas)
        right_box.addWidget(controls_group)
        right_box.addWidget(QLabel("Log"))
        right_box.addWidget(self.log_text, 1)

        root.addLayout(left_box, 2)
        root.addLayout(right_box, 1)

    def _apply_default_state(self):
        self.name_edit.setText("MaartenAfschrift")
        self.crank_length_spin.setValue(0.175)
        self.saddle_rx_spin.setValue(-0.05)
        self.saddle_ry_spin.setValue(0.87)
        self.freq_spin.setValue(90.0)
        self.objective_edit.setText("maxpower")
        self._append_log("Ready.")
        self._set_ready(False)

    def _set_ready(self, ready: bool):
        color = "#39a845" if ready else "#c0392b"
        self.ready_lamp.setStyleSheet(f"background-color: {color}; border-radius: 9px;")

    def _append_log(self, message: str):
        self.log_text.append(message)

    def clear_log(self):
        self.log_text.clear()
        self._append_log("Log cleared.")

    def _build_settings(self) -> CyclingSettings:
        return CyclingSettings(
            name=self.name_edit.text().strip() or "simulation",
            muscle_model="Leuven",
            msk_model="Kistemaker",
            specific_tension=25.0,
            crank_frequency_rpm=float(self.freq_spin.value()),
            is_isokinetic=True,
            saddle_rx=float(self.saddle_rx_spin.value()),
            saddle_ry=float(self.saddle_ry_spin.value()),
            crank_length=float(self.crank_length_spin.value()),
            fixed_power=None,
            opt_frequency=False,
            collocation_points=50,
            metabolism_model="Bhargava2004",
            metabolism_scale_rate=1.0,
            bool_plot=False,
            scaling=Scaling(),
            objective=ObjectiveSettings(type="maxpower"),
            solver=SolverSettings(),
        )

    def simulate(self):
        self._set_ready(False)
        self._append_log("")
        self._append_log("Start simulation with settings:")
        self._append_log(f"  name: {self.name_edit.text().strip() or 'simulation'}")
        self._append_log(f"  crank length: {self.crank_length_spin.value():.3f} m")
        self._append_log(f"  saddle: ({self.saddle_rx_spin.value():.3f}, {self.saddle_ry_spin.value():.3f}) m")
        self._append_log(f"  frequency: {self.freq_spin.value():.1f} rot/min")
        self._append_log("  running...")
        app = QApplication.instance()
        if app is not None:
            app.processEvents()

        try:
            settings = self._build_settings()
            self.simulator = CyclingSimulator(settings)
            self.result, json_path = self.simulator.simulate_and_export(name=settings.name, results_dir=self.results_dir)
            self.history.append(self.result.average_power)
            self.history_names.append(settings.name)
            self._update_score_plot()
            self._draw_initial_visualisation()
            self._append_log("Simulation finished.")
            self._append_log(f"Average mechanical power: {self.result.average_power:.2f} W")
            self._append_log(f"Objective value: {self.result.objective_value:.4f}")
            self._append_log(f"Result saved to: {json_path}")
            self._set_ready(True)
        except Exception as exc:  # pragma: no cover - GUI error handling
            self._set_ready(False)
            self._append_log("Optimizer failed.")
            self._append_log(str(exc))
            QMessageBox.critical(self, "Simulation failed", str(exc))

    def _update_score_plot(self):
        self.score_ax.clear()
        self.score_ax.set_title("Average mechanical power")
        self.score_ax.set_xlabel("attempts")
        self.score_ax.set_ylabel("power [W]")
        self.score_ax.grid(True, alpha=0.3)
        if self.history:
            x = np.arange(1, len(self.history) + 1)
            self.score_ax.plot(x, self.history, "o", color="tab:blue")
            self.score_ax.plot(x, self.history, color="tab:blue", alpha=0.4)
            max_power = max(self.history)
            self.score_ax.axhline(max_power, color="k", linestyle="--", linewidth=1)
            self.score_ax.text(0.5, max_power, f"max {max_power:.1f} W", fontsize=9)
        self.score_canvas.draw_idle()

    def _draw_initial_visualisation(self):
        if self.result is None:
            return
        self.visual_ax.clear()
        self.visual_ax.set_title("Cycling Model")
        self.visual_ax.set_xlabel("X [m]")
        self.visual_ax.set_ylabel("Y [m]")
        self.visual_ax.set_aspect("equal", adjustable="box")
        self.visual_ax.grid(True, alpha=0.3)
        self.visual_ax.set_xlim(-0.5, 0.5)
        self.visual_ax.set_ylim(-0.3, 1.1)

        rj = self.result.r_joint
        segment_colors = ["k", (0, 0.5, 0), "b", "r"]
        segment_names = ["Crank", "Foot", "Lower Leg", "Upper Leg"]
        self._animation_lines = []
        self._animation_labels = []

        for i in range(4):
            line, = self.visual_ax.plot(
                [rj[0, i, 0], rj[0, i + 1, 0]],
                [rj[1, i, 0], rj[1, i + 1, 0]],
                "o-",
                color=segment_colors[i],
                linewidth=3,
                markersize=5,
                markerfacecolor=segment_colors[i],
            )
            self._animation_lines.append(line)
            txt = self.visual_ax.text(
                float(np.mean([rj[0, i, 0], rj[0, i + 1, 0]])),
                float(np.mean([rj[1, i, 0], rj[1, i + 1, 0]])),
                f"  {segment_names[i]}",
                color=segment_colors[i],
                fontsize=9,
                weight="bold",
            )
            self._animation_labels.append(txt)

        crank_radius = float(np.linalg.norm(rj[:, 0, 0] - rj[:, 1, 0]))
        theta = np.linspace(0.0, 2.0 * np.pi, 200)
        self._animation_circle, = self.visual_ax.plot(
            rj[0, 0, 0] + crank_radius * np.cos(theta),
            rj[1, 0, 0] + crank_radius * np.sin(theta),
            "--",
            color=(0.6, 0.6, 0.6),
        )
        self._animation_saddle, = self.visual_ax.plot(
            rj[0, 4, 0],
            rj[1, 4, 0],
            "s",
            markersize=9,
            markeredgecolor="k",
            markerfacecolor="y",
        )
        self.visual_canvas.draw_idle()

    def visualise(self):
        if self.result is None:
            self._append_log("No simulation result available yet.")
            return
        if self._animation_timer.isActive():
            self._animation_timer.stop()
        self._animation_frame = 0
        dt = float(np.mean(np.diff(self.result.t))) if len(self.result.t) > 1 else 0.02
        self._animation_interval_ms = max(1, int(1000.0 * dt * 0.92))
        self._append_log("Starting visualisation...")
        self._animation_timer.start(self._animation_interval_ms)

    def _advance_animation(self):
        if self.result is None:
            self._animation_timer.stop()
            return

        rj = self.result.r_joint
        n_frames = rj.shape[2]
        if self._animation_frame >= n_frames:
            self._animation_timer.stop()
            return

        frame = self._animation_frame
        for i, line in enumerate(self._animation_lines):
            line.set_data([rj[0, i, frame], rj[0, i + 1, frame]], [rj[1, i, frame], rj[1, i + 1, frame]])
            self._animation_labels[i].set_position(
                (
                    float(np.mean([rj[0, i, frame], rj[0, i + 1, frame]])),
                    float(np.mean([rj[1, i, frame], rj[1, i + 1, frame]])),
                )
            )
        self._animation_saddle.set_data([rj[0, 4, frame]], [rj[1, 4, frame]])
        self.visual_canvas.draw_idle()
        self._animation_frame += 1


