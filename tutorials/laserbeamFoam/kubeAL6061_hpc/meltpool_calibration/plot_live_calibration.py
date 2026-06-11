#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
import math
import re
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


def parse_scalar(path: Path, name: str) -> float:
    text = path.read_text()
    match = re.search(rf"\b{name}\s+([^;]+);", text)
    if not match:
        return math.nan
    return float(match.group(1))


def mesh_size_um(case_dir: Path) -> float:
    text = (case_dir / "system" / "blockMeshDict").read_text()
    vertices = [
        tuple(float(match.group(i)) for i in range(1, 4))
        for match in re.finditer(r"\(\s*([0-9.eE+-]+)\s+([0-9.eE+-]+)\s+([0-9.eE+-]+)\s*\)", text)
    ][:8]
    cell_match = re.search(
        r"hex\s*\(0\s+1\s+2\s+3\s+4\s+5\s+6\s+7\)\s*\(\s*(\d+)\s+(\d+)\s+(\d+)\s*\)",
        text,
    )
    if len(vertices) != 8 or not cell_match:
        return math.nan
    xs = [p[0] for p in vertices]
    nx = int(cell_match.group(1))
    return (max(xs) - min(xs)) / nx * 1e6


def exp_onset_ms(exp_dir: Path) -> float:
    data = pd.read_csv(exp_dir / "keyhole_depth_exp.csv")
    depth = pd.to_numeric(data["depth_um"], errors="coerce")
    active = data[depth > 0.0]
    if active.empty:
        return float(data["time_ms"].min())
    return float(active["time_ms"].iloc[0])


def load_exp(exp_dir: Path, filename: str, column: str, onset_ms: float, shift_ms: float) -> tuple[np.ndarray, np.ndarray]:
    data = pd.read_csv(exp_dir / filename)
    data = data[pd.to_numeric(data["time_ms"], errors="coerce") >= onset_ms].copy()
    time_ms = pd.to_numeric(data["time_ms"], errors="coerce").to_numpy(dtype=float) - shift_ms
    values = pd.to_numeric(data[column], errors="coerce").to_numpy(dtype=float)
    finite = np.isfinite(time_ms) & np.isfinite(values)
    return time_ms[finite], values[finite]


def sim_series(case_dir: Path, metric: str, config: dict) -> tuple[np.ndarray, np.ndarray]:
    data = pd.read_csv(case_dir / "post-processing-data" / "vtu_meltpool_geometry.csv")
    time_ms = pd.to_numeric(data["time"], errors="coerce").to_numpy(dtype=float) * 1e3
    time_ms = time_ms + float(config.get("sim_time_shift_ms", 0.0))

    radius_um = parse_scalar(case_dir / "constant" / "LaserProperties", "laserRadius") * 1e6
    beam_diameter_um = 2.0 * radius_um if bool(config.get("keyhole_ar_uses_case_beam_diameter", True)) else 80.0

    if metric == "keyhole_depth":
        values = pd.to_numeric(data["keyholeDepth_um"], errors="coerce").to_numpy(dtype=float)
    elif metric == "keyhole_ar":
        values = pd.to_numeric(data["keyholeDepth_um"], errors="coerce").to_numpy(dtype=float) / beam_diameter_um
    elif metric == "meltpool_depth":
        values = pd.to_numeric(data["meltPoolDepth_um"], errors="coerce").to_numpy(dtype=float)
    elif metric == "meltpool_ar":
        depth = pd.to_numeric(data["meltPoolDepth_um"], errors="coerce").to_numpy(dtype=float)
        width = pd.to_numeric(data["meltPoolWidth_um"], errors="coerce").to_numpy(dtype=float)
        values = np.divide(width, depth, out=np.full_like(width, np.nan), where=depth > 0)
    else:
        raise ValueError(metric)

    finite = np.isfinite(time_ms) & np.isfinite(values)
    return time_ms[finite], values[finite]


def case_label(case_dir: Path, config: dict) -> str:
    radius_um = parse_scalar(case_dir / "constant" / "LaserProperties", "laserRadius") * 1e6
    er = parse_scalar(case_dir / "constant" / "transportProperties", "elec_resistivity")
    mesh = mesh_size_um(case_dir)
    sim_shift = float(config.get("sim_time_shift_ms", 0.0))
    return f"m={mesh:g} um, r={radius_um:g} um, er={er:.2e}, sim +{sim_shift:.3f} ms"


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--base-case", type=Path, required=True)
    parser.add_argument("--runs-dir", type=Path, required=True)
    parser.add_argument("--config", type=Path, required=True)
    parser.add_argument("--output", type=Path, default=None)
    args = parser.parse_args()

    base_case = args.base_case.resolve()
    runs_dir = args.runs_dir.resolve()
    config = json.loads(args.config.read_text())
    output = args.output or base_case / "meltpool_calibration" / "plots" / "live_sim_vs_exp_transient.png"
    output.parent.mkdir(parents=True, exist_ok=True)

    run_cases = sorted(
        case_dir for case_dir in runs_dir.iterdir()
        if (case_dir / "post-processing-data" / "vtu_meltpool_geometry.csv").exists()
    ) if runs_dir.exists() else []

    exp_dir = base_case / "kube_exp_data" / "csv"
    exp_shift = float(config.get("exp_time_shift_ms", 0.0))
    onset_ms = exp_onset_ms(exp_dir) if bool(config.get("plot_from_keyhole_onset", True)) else -math.inf
    x_left = onset_ms - exp_shift if math.isfinite(onset_ms) else math.inf

    metrics = [
        ("keyhole_depth", "keyhole_depth_exp.csv", "depth_um", "Keyhole depth", "Depth (um)"),
        ("keyhole_ar", "keyhole_ar.csv", "aspect_ratio_keyholeDepthDividedbyLaserbeamdia80um", "Keyhole aspect ratio", "Depth / beam diameter"),
        ("meltpool_depth", "meltpool_depth_exp.csv", "depth_um", "Melt-pool depth", "Depth (um)"),
        ("meltpool_ar", "meltpool_ar_exp.csv", "aspect_ratio_meltpoolWidthbyDepth", "Melt-pool aspect ratio", "Width / depth"),
    ]

    fig, axes = plt.subplots(2, 2, figsize=(14, 9), constrained_layout=True)
    axes = axes.ravel()

    for ax, (metric, exp_file, exp_col, title, ylabel) in zip(axes, metrics):
        exp_t, exp_y = load_exp(exp_dir, exp_file, exp_col, onset_ms, exp_shift)
        ax.plot(exp_t, exp_y, color="black", linewidth=2.3, label="Experiment")

        for case_dir in run_cases:
            sim_t, sim_y = sim_series(case_dir, metric, config)
            if sim_t.size:
                x_left = min(x_left, float(np.nanmin(sim_t)))
                ax.plot(sim_t, sim_y, marker="o", markersize=2.2, linewidth=1.1, label=case_label(case_dir, config))

        if math.isfinite(x_left):
            ax.set_xlim(left=x_left)
        ax.set_title(title)
        ax.set_xlabel("Time (ms)")
        ax.set_ylabel(ylabel)
        ax.grid(True, alpha=0.3)

    handles, labels = axes[0].get_legend_handles_labels()
    fig.legend(handles, labels, loc="lower center", bbox_to_anchor=(0.5, -0.08), ncol=3, fontsize=8)
    fig.savefig(output, dpi=160, bbox_inches="tight")
    plt.close(fig)
    print(f"Wrote {output}")


if __name__ == "__main__":
    main()
