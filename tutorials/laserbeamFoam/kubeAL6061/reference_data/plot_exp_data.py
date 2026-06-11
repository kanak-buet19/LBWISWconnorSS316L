#!/usr/bin/env python3
"""Plot Kube no-sonication experimental digitized data in subplots."""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd


ROOT = Path(__file__).resolve().parent
CSV_DIR = ROOT / "csv"
OUT_DIR = ROOT / "images"
OUT_FILE = OUT_DIR / "kube_exp_data_subplots.png"
EXP_TIME_SHIFT_MS = 2.661


PLOTS = [
    (
        "keyhole_depth_exp.csv",
        "depth_um",
        "Keyhole depth",
        "Depth (um)",
        "#D62728",
    ),
    (
        "keyhole_ar.csv",
        "aspect_ratio_keyholeDepthDividedbyLaserbeamdia80um",
        "Keyhole aspect ratio",
        "Depth / 80 um",
        "#D62728",
    ),
    (
        "meltpool_depth_exp.csv",
        "depth_um",
        "Melt-pool depth",
        "Depth (um)",
        "#F2A900",
    ),
    (
        "meltpool_ar_exp.csv",
        "aspect_ratio_meltpoolWidthbyDepth",
        "Melt-pool aspect ratio",
        "Width / depth",
        "#F2A900",
    ),
]


def load_csv(filename: str, value_column: str) -> pd.DataFrame:
    path = CSV_DIR / filename
    if not path.exists():
        raise FileNotFoundError(f"Missing CSV: {path}")

    data = pd.read_csv(path)
    required = {"time_ms", value_column}
    missing = required.difference(data.columns)
    if missing:
        raise ValueError(f"{path} missing columns: {', '.join(sorted(missing))}")

    return data.sort_values("time_ms")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Plot Kube no-sonication experimental data in shifted time."
    )
    parser.add_argument(
        "--time-shift-ms",
        type=float,
        default=EXP_TIME_SHIFT_MS,
        help="Subtract this value from experimental time_ms.",
    )
    parser.add_argument(
        "--no-shift",
        action="store_true",
        help="Plot original experimental time without shifting.",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    time_shift_ms = 0.0 if args.no_shift else args.time_shift_ms

    fig, axes = plt.subplots(2, 2, figsize=(10.5, 7.0), constrained_layout=True)

    for ax, (filename, value_column, title, ylabel, color) in zip(axes.flat, PLOTS):
        data = load_csv(filename, value_column)
        ax.plot(
            data["time_ms"],
            data[value_column],
            color=color,
            linewidth=1.8,
        )
        ax.set_title(title)
        ax.set_xlabel("Time (ms)")
        ax.set_ylabel(ylabel)
        ax.grid(True, alpha=0.3)

    OUT_DIR.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUT_FILE, dpi=220)
    plt.close(fig)
    print(f"Wrote {OUT_FILE}")


if __name__ == "__main__":
    main()
