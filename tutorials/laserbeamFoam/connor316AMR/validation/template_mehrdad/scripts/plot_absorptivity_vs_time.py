#!/usr/bin/env python3
"""Plot absorptivity history written by laserHeatSource.

Default input:
  absorptivity_vs_time/absorptivity_vs_time.csv

Default outputs:
  post-processing-data/absorptivity_vs_time/absorptivity_vs_time.pdf
  post-processing-data/absorptivity_vs_time/absorptivity_vs_time.png
"""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd


CASE_DIR = Path(__file__).resolve().parents[1]
DEFAULT_CSV = CASE_DIR / "absorptivity_vs_time" / "absorptivity_vs_time.csv"
DEFAULT_OUT_DIR = CASE_DIR / "post-processing-data" / "absorptivity_vs_time"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Plot laser absorptivity versus simulation time."
    )
    parser.add_argument(
        "--csv",
        type=Path,
        default=DEFAULT_CSV,
        help=f"Input CSV file. Default: {DEFAULT_CSV}",
    )
    parser.add_argument(
        "--out-dir",
        type=Path,
        default=DEFAULT_OUT_DIR,
        help=f"Output directory. Default: {DEFAULT_OUT_DIR}",
    )
    parser.add_argument(
        "--name",
        default="absorptivity_vs_time",
        help="Base filename for saved plot.",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    csv_path = args.csv.expanduser().resolve()

    if not csv_path.exists():
        raise FileNotFoundError(f"Cannot find absorptivity CSV: {csv_path}")

    data = pd.read_csv(csv_path)
    required = {"time_s", "laser", "absorptivity"}
    missing = required.difference(data.columns)
    if missing:
        raise ValueError(
            f"CSV missing required columns: {', '.join(sorted(missing))}"
        )

    data = data.sort_values("time_s")
    data["time_us"] = data["time_s"] * 1.0e6

    fig, ax = plt.subplots(figsize=(7.0, 4.2), constrained_layout=True)

    for laser_name, laser_data in data.groupby("laser", sort=False):
        ax.plot(
            laser_data["time_us"],
            laser_data["absorptivity"],
            marker="o",
            markersize=2.5,
            linewidth=1.3,
            label=str(laser_name),
        )

    ax.set_xlabel("Time (microsec)")
    ax.set_ylabel("Absorptivity")
    ax.set_ylim(bottom=0)
    ax.grid(True, alpha=0.3)

    if data["laser"].nunique() > 1:
        ax.legend(frameon=True)

    args.out_dir.mkdir(parents=True, exist_ok=True)
    fig.savefig(args.out_dir / f"{args.name}.pdf")
    fig.savefig(args.out_dir / f"{args.name}.png", dpi=200)
    plt.close(fig)


if __name__ == "__main__":
    main()
