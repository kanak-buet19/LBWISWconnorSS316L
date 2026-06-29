"""
Stationary laser absorptivity comparison: experiment vs simulation.

Run from repo root:
    python3 scripts/plot_stationary_abs_exp_sim.py
"""

from pathlib import Path

import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import pandas as pd


ROOT = Path(__file__).parent.parent

EXP_PATH = ROOT / "stationary_laser" / "absorption_data.csv"
SIM_PATH = Path(
    "/home/kanak/coding_ubuntu/LPBF_solvers/LaserbeamFoam_AMR/"
    "tutorials/laserbeamFoam/connor316AMR/validation/cases/"
    "tao_Ti64_stationary_156W_r70um_rho9p2em07/"
    "absorptivity_vs_time/absorptivity_vs_time.csv"
)


plt.rcParams.update({
    "font.family": "DejaVu Sans",
    "font.size": 11,
    "axes.titlesize": 12,
    "axes.labelsize": 11,
    "xtick.labelsize": 9,
    "ytick.labelsize": 9,
    "legend.fontsize": 9,
    "axes.spines.top": False,
    "axes.spines.right": False,
    "axes.grid": True,
    "grid.alpha": 0.35,
    "grid.linewidth": 0.6,
    "lines.linewidth": 1.2,
    "figure.dpi": 150,
})


def load_experiment(path):
    df = pd.read_csv(path, index_col=0, low_memory=False)
    df["Time_ms"] = df["Time"] * 1000.0

    laser_on = df["InputLaser"] > 0
    return df.loc[laser_on, ["Time_ms", "RelativeAbsorption"]].rename(
        columns={"RelativeAbsorption": "Absorptivity_pct"}
    )


def load_simulation(path):
    df = pd.read_csv(path, low_memory=False)
    df["Time_ms"] = df["time_s"] * 1000.0
    df["Absorptivity_pct"] = df["absorptivity"] * 100.0
    return df[["Time_ms", "Absorptivity_pct"]]


exp = load_experiment(EXP_PATH)
sim = load_simulation(SIM_PATH)

fig, ax = plt.subplots(figsize=(8.2, 4.8))

ax.plot(
    exp["Time_ms"],
    exp["Absorptivity_pct"],
    color="#2166ac",
    label="Experiment: NIST stationary laser",
    linewidth=1.0,
)
ax.plot(
    sim["Time_ms"],
    sim["Absorptivity_pct"],
    color="#b2182b",
    label="Simulation: 156 W, r = 70 um",
    linewidth=1.4,
)

ax.set_title("Stationary Ti-6Al-4V Laser Absorptivity")
ax.set_xlabel("Time (ms)")
ax.set_ylabel("Absorptivity (%)")
ax.set_ylim(bottom=0)
ax.yaxis.set_major_locator(ticker.MultipleLocator(10))
ax.yaxis.set_minor_locator(ticker.AutoMinorLocator())
ax.xaxis.set_minor_locator(ticker.AutoMinorLocator())
ax.legend(loc="best", framealpha=0.9)

fig.tight_layout()

for ext in ("png", "pdf"):
    out = ROOT / "scripts" / f"stationary_abs_exp_sim.{ext}"
    fig.savefig(out, dpi=150, bbox_inches="tight")
    print(f"Saved -> {out}")

plt.show()
