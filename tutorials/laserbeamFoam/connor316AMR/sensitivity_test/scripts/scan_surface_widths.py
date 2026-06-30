#!/usr/bin/env python3
"""Compute comparable melt-pool surface widths from final VTU files.

The template analyzer chooses an XY slice from keyhole/interface logic. For
shallow or weak-interface cases that can fall back to an arbitrary z plane.
This script scans a small set of fixed z planes and reports the maximum
surface liquidus width, which is more stable for sensitivity ranking.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
from pathlib import Path
from typing import Any

import numpy as np
import pyvista as pv
from scipy.interpolate import griddata


ROOT = Path(__file__).resolve().parents[1]
CONFIG_PATH = ROOT / "sensitivity_config.json"
Z_PLANES_UM = (515.0, 600.0, 640.0, 679.3, 720.0, 765.0, 840.0)
T_LIQUIDUS_K = 1723.0
NX = 240
NY = 180


def load_config() -> dict[str, Any]:
    return json.loads(CONFIG_PATH.read_text())


def latest_vtu(case_dir: Path) -> Path | None:
    files = list((case_dir / "VTK").glob("*/internal.vtu"))
    if not files:
        return None
    return max(files, key=lambda path: path.stat().st_mtime)


def surface_width_at_z(mesh: pv.DataSet, z_um: float, surface_y_um: float) -> tuple[float, float]:
    slc = mesh.slice(normal=(0.0, 0.0, 1.0), origin=(0.0, 0.0, z_um * 1e-6))
    if slc.n_points < 3 or "T" not in slc.point_data or "alpha.metal" not in slc.point_data:
        return math.nan, math.nan

    points = slc.points
    x_um = points[:, 0] * 1e6
    y_um = points[:, 1] * 1e6
    values = np.column_stack([x_um, y_um])
    x_grid = np.linspace(float(x_um.min()), float(x_um.max()), NX)
    y_grid = np.linspace(float(y_um.min()), float(y_um.max()), NY)
    X, Y = np.meshgrid(x_grid, y_grid)

    temp = np.asarray(slc.point_data["T"])
    alpha = np.asarray(slc.point_data["alpha.metal"])
    temp_grid = griddata(values, temp, (X, Y), method="linear")
    alpha_grid = griddata(values, alpha, (X, Y), method="linear")
    temp_near = griddata(values, temp, (X, Y), method="nearest")
    alpha_near = griddata(values, alpha, (X, Y), method="nearest")
    temp_grid[np.isnan(temp_grid)] = temp_near[np.isnan(temp_grid)]
    alpha_grid[np.isnan(alpha_grid)] = alpha_near[np.isnan(alpha_grid)]

    temp_profile = np.array([np.interp(surface_y_um, y_grid, temp_grid[:, i]) for i in range(NX)])
    alpha_profile = np.array([np.interp(surface_y_um, y_grid, alpha_grid[:, i]) for i in range(NX)])
    melted = (temp_profile >= T_LIQUIDUS_K) & (alpha_profile > 0.1)
    if not np.any(melted):
        return math.nan, float(np.nanmax(temp_profile))

    indices = np.flatnonzero(melted)
    intervals: list[tuple[int, int]] = []
    start = int(indices[0])
    previous = int(indices[0])
    for index in indices[1:]:
        index = int(index)
        if index == previous + 1:
            previous = index
        else:
            intervals.append((start, previous))
            start = previous = index
    intervals.append((start, previous))

    width = max(float(x_grid[right] - x_grid[left]) for left, right in intervals)
    return width, float(np.nanmax(temp_profile))


def scan_case(row: dict[str, Any], surface_y_um: float) -> dict[str, Any]:
    case_dir = ROOT / row["case_dir"]
    out = {
        **row,
        "correctedSurfaceWidth_um": math.nan,
        "correctedWidthZ_um": math.nan,
        "maxSurfaceTemperature_K": math.nan,
        "vtu_file": "",
    }
    vtu = latest_vtu(case_dir)
    if vtu is None:
        return out

    mesh = pv.read(vtu).cell_data_to_point_data(pass_cell_data=True)
    best_width = math.nan
    best_z = math.nan
    best_temp = math.nan
    for z_um in Z_PLANES_UM:
        width, max_temp = surface_width_at_z(mesh, z_um, surface_y_um)
        if math.isfinite(width) and (not math.isfinite(best_width) or width > best_width):
            best_width = width
            best_z = z_um
            best_temp = max_temp

    out["correctedSurfaceWidth_um"] = best_width
    out["correctedWidthZ_um"] = best_z
    out["maxSurfaceTemperature_K"] = best_temp
    out["vtu_file"] = str(vtu.relative_to(ROOT))
    return out


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--quiet", action="store_true")
    args = parser.parse_args()

    config = load_config()
    results_path = ROOT / config["results_root"] / "sensitivity_results.csv"
    if not results_path.exists():
        raise SystemExit("No sensitivity_results.csv found")

    rows = list(csv.DictReader(results_path.open()))
    surface_y_um = float(config["run_settings"]["surface_y_m"]) * 1e6
    scanned = [scan_case(row, surface_y_um) for row in rows]

    out_path = ROOT / config["results_root"] / "corrected_surface_widths.csv"
    fields = [
        "case_id",
        "parameter",
        "label",
        "direction",
        "multiplier",
        "baseline_value",
        "varied_value",
        "unit",
        "status",
        "meltPoolWidth_um",
        "correctedSurfaceWidth_um",
        "correctedWidthZ_um",
        "maxSurfaceTemperature_K",
        "case_dir",
        "vtu_file",
    ]
    with out_path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(scanned)

    if not args.quiet:
        print(f"Wrote {out_path}")


if __name__ == "__main__":
    main()
