"""
Absorption timeseries — both laser experiments, 3-panel figure.
Run from repo root:  python scripts/plot_absorption.py
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
from pathlib import Path

ROOT = Path(__file__).parent.parent

EXPERIMENTS = {
    "Stationary laser (spot, 2 ms pulse)":   ROOT / "stationary_laser" / "absorption_data.csv",
    "Moving laser (scan, 700 mm/s)":         ROOT / "moving_laser"     / "absorption_data.csv",
}
COLORS = {
    "Stationary laser (spot, 2 ms pulse)":  "#2166ac",
    "Moving laser (scan, 700 mm/s)":        "#d6604d",
}

plt.rcParams.update({
    "font.family":      "DejaVu Sans",
    "font.size":        11,
    "axes.titlesize":   12,
    "axes.labelsize":   11,
    "xtick.labelsize":  9,
    "ytick.labelsize":  9,
    "legend.fontsize":  9,
    "axes.spines.top":  False,
    "axes.spines.right":False,
    "axes.grid":        True,
    "grid.alpha":       0.35,
    "grid.linewidth":   0.6,
    "lines.linewidth":  1.0,
    "figure.dpi":       150,
})


def load(path):
    df = pd.read_csv(path, index_col=0, low_memory=False)
    df["AbsAbsorptionUncertainty"] = pd.to_numeric(
        df["AbsAbsorptionUncertainty"], errors="coerce"
    )
    df["Time_ms"] = df["Time"] * 1000
    df["Time_us"] = df["Time"] * 1e6
    return df


def shade_frame_captures(ax, df, color, alpha=0.07):
    """Shade time regions where FrameTrigger == 1."""
    in_frame, t0 = False, None
    for _, row in df.iterrows():
        if row["FrameTrigger"] == 1 and not in_frame:
            t0, in_frame = row["Time_ms"], True
        elif row["FrameTrigger"] == 0 and in_frame:
            ax.axvspan(t0, row["Time_ms"], color=color, alpha=alpha, linewidth=0)
            in_frame = False


datasets = {name: load(p) for name, p in EXPERIMENTS.items()}

fig, axes = plt.subplots(3, 1, figsize=(11, 12))
fig.subplots_adjust(top=0.91, hspace=0.38)
fig.suptitle(
    "NIST AMB2022 — Ti-6Al-4V Laser Absorptance\n"
    "Simultaneous high-speed X-ray & integrating sphere radiometry, APS 32-ID-B",
    fontsize=13, fontweight="bold", y=0.97,
)

# ── Panel 1: Input laser power ────────────────────────────────────────────────
ax = axes[0]
for name, df in datasets.items():
    mask = df["InputLaser"] > 0
    ax.plot(df.loc[mask, "Time_us"], df.loc[mask, "InputLaser"],
            color=COLORS[name], label=name, linewidth=0.9)

ax.set_xlabel("Time (μs)")
ax.set_ylabel("Input Laser Power (W)")
ax.set_title("Input Laser Temporal Profile")
ax.legend(loc="upper right", framealpha=0.9)
ax.set_xlim(left=-100)
ax.yaxis.set_minor_locator(ticker.AutoMinorLocator())

# ── Panel 2: Relative absorption ──────────────────────────────────────────────
ax = axes[1]
for name, df in datasets.items():
    c = COLORS[name]
    laser_on = df["InputLaser"] > 0
    sub = df[laser_on].reset_index(drop=True)
    ax.plot(sub["Time_ms"], sub["RelativeAbsorption"],
            color=c, label=name, linewidth=0.9)
    shade_frame_captures(ax, sub, color=c)

ax.set_xlabel("Time (ms)")
ax.set_ylabel("Relative Absorptance (%)")
ax.set_title("Relative Absorptance  [shaded regions = X-ray frames captured]")
ax.set_ylim(0, 110)
ax.yaxis.set_major_locator(ticker.MultipleLocator(20))
ax.legend(loc="upper right", framealpha=0.9)

# ── Panel 3: Absolute absorbed power + uncertainty ────────────────────────────
ax = axes[2]
for name, df in datasets.items():
    c = COLORS[name]
    sub = df[df["InputLaser"] > 0].copy()
    t = sub["Time_ms"]
    a = sub["AbsoluteAbsorption"]
    u = sub["AbsAbsorptionUncertainty"]
    ax.plot(t, a, color=c, label=name, linewidth=0.9)
    if u.notna().any():
        ax.fill_between(t, a - u, a + u, color=c, alpha=0.22,
                        label=f"{name.split('(')[0].strip()} ± expanded uncertainty")

ax.set_xlabel("Time (ms)")
ax.set_ylabel("Absorbed Power (W)")
ax.set_title("Absolute Absorbed Power with Expanded Measurement Uncertainty")
ax.legend(loc="upper right", framealpha=0.9)
ax.yaxis.set_minor_locator(ticker.AutoMinorLocator())

for ext in ("png", "pdf"):
    out = ROOT / "scripts" / f"absorption_timeseries.{ext}"
    plt.savefig(out, dpi=150, bbox_inches="tight")
    print(f"Saved → {out}")
plt.show()
