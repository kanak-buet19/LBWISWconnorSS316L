#!/usr/bin/env python3
"""Measure melt-pool/keyhole metrics from reconstructed AMR VTU files.

Outputs:
  post-processing-data/vtu_meltpool_geometry.csv
  post-processing-data/vtu_sections/*.png
"""

from __future__ import annotations

import argparse
import csv
import json
import re
from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import numpy as np
import pandas as pd
import pyvista as pv
from scipy.interpolate import griddata


NX, NZ, NY = 180, 180, 180
INTERP_METHOD = "linear"
SURFACE_Y_UM = None
PDF_DPI = 150
PNG_DPI = 180
EXP_TARGET_WIDTH_UM = 436.0
EXP_TARGET_DEPTH_UM = 354.0
EXP_TARGET_DEPTH_WIDTH = 0.8
EXP_TIME_SHIFT_MS = 0.0
CASE_DIR = Path(__file__).resolve().parents[1]
TAOSUN_EXP_CSV_DIR = CASE_DIR / "taosun_exp_data" / "csv"
ABSORPTIVITY_CSV = CASE_DIR / "absorptivity_vs_time" / "absorptivity_vs_time.csv"
LASER_DIAMETER_UM = 120.0


CSV_FIELDS = [
    "time",
    "z_slice_um",
    "surface_y_um",
    "T_liq_K",
    "keyholeDepth_um",
    "keyholeDepthX_um",
    "keyholeDepthY_um",
    "meltPoolDepth_um",
    "meltPoolDepthX_um",
    "meltPoolDepthY_um",
    "keyholeWidth_um",
    "keyholeWidthXLeft_um",
    "keyholeWidthXRight_um",
    "meltPoolWidth_um",
    "meltPoolWidthXLeft_um",
    "meltPoolWidthXRight_um",
]


def parse_scalar(path: Path, name: str) -> float:
    text = path.read_text()
    match = re.search(rf"\b{name}\s+([^;]+);", text)
    if not match:
        raise ValueError(f"Could not find {name} in {path}")
    return float(match.group(1))


def infer_surface_y_um(case: Path) -> float:
    """Infer the initial gas/metal interface from setFieldsDict.

    This case initializes metal with a boxToCell region. Since y is the
    laser-entry/depth direction and the gas sits above the metal, the lower
    y-coordinate of the alpha.metal = 1 box is the initial free surface.
    """
    set_fields = case / "system" / "setFieldsDict"
    text = set_fields.read_text()

    box_match = re.search(
        r"box\s*"
        r"\(\s*([^\s()]+)\s+([^\s()]+)\s+([^\s()]+)\s*\)\s*"
        r"\(\s*([^\s()]+)\s+([^\s()]+)\s+([^\s()]+)\s*\)\s*;",
        text,
    )
    if not box_match:
        raise ValueError(f"Could not infer surface y from {set_fields}")

    y_min = float(box_match.group(2))
    return y_min * 1e6


def vtk_time(path: Path) -> float:
    # Match format: case_time_serial.vtk
    match = re.search(r"_([0-9.eE+-]+)_[0-9]+\.(?:vtk|vtu)$", path.name)
    if match:
        return float(match.group(1))

    # Match format: case_time.vtk or case_time.vtu
    match = re.search(r"_([0-9.eE+-]+)\.(?:vtk|vtu)$", path.name)
    if match:
        return float(match.group(1))

    # Match format: directory name ending with _time
    match = re.search(r"_([0-9.eE+-]+)$", path.parent.name)
    if match:
        return float(match.group(1))

    with path.open("rb") as handle:
        header = handle.read(512).decode("utf-8", errors="ignore")

    match = re.search(r"time='([0-9.eE+-]+)'", header)
    if match:
        return float(match.group(1))

    series_time = vtk_series_time(path)
    if series_time is not None:
        return series_time

    raise ValueError(f"Cannot parse time from {path}")


def vtk_series_time(path: Path) -> float | None:
    for series_path in path.parent.glob("*.vtk.series"):
        try:
            series = json.loads(series_path.read_text())
        except (OSError, json.JSONDecodeError):
            continue

        for item in series.get("files", []):
            if Path(item.get("name", "")).name != path.name:
                continue
            try:
                return float(item["time"])
            except (KeyError, TypeError, ValueError):
                return None

    return None


def converted_vtk_files(vtk_dir: Path) -> list[Path]:
    candidates = [
        *vtk_dir.glob("*/internal.vtu"),
        *vtk_dir.glob("*/internal.vtk"),
        *vtk_dir.glob("*.vtk"),
    ]
    return sorted(candidates, key=vtk_time)


def field(mesh: pv.DataSet, name: str) -> np.ndarray:
    if name in mesh.cell_data:
        return np.asarray(mesh.cell_data[name])
    if name in mesh.point_data:
        return np.asarray(mesh.point_data[name])
    raise KeyError(f"{name} not found in {mesh}")


def slice_to_dataframe(mesh: pv.DataSet, normal: tuple[float, float, float], origin: tuple[float, float, float]) -> pd.DataFrame:
    source = mesh.cell_data_to_point_data(pass_cell_data=True)
    slc = source.slice(normal=normal, origin=origin)
    if slc.n_points < 3:
        return pd.DataFrame()

    if "alpha.metal" not in slc.point_data or "T" not in slc.point_data:
        available = sorted(set(slc.point_data.keys()) | set(slc.cell_data.keys()))
        raise KeyError(f"Slice missing alpha.metal or T. Available fields: {available}")

    xyz = slc.points
    df = pd.DataFrame(
        {
            "x": xyz[:, 0],
            "y": xyz[:, 1],
            "z": xyz[:, 2],
            "alpha.metal": np.asarray(slc.point_data["alpha.metal"]),
            "T": np.asarray(slc.point_data["T"]),
        }
    )
    df["x_um"] = df["x"] * 1e6
    df["y_um"] = df["y"] * 1e6
    df["z_um"] = df["z"] * 1e6
    return df


def interpolate_slice(
    slice_df: pd.DataFrame,
    x_col: str,
    y_col: str,
    nx: int,
    ny: int,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ma.MaskedArray]:
    points = slice_df[[x_col, y_col]].to_numpy()
    xg = np.linspace(slice_df[x_col].min(), slice_df[x_col].max(), nx)
    yg = np.linspace(slice_df[y_col].min(), slice_df[y_col].max(), ny)
    Xg, Yg = np.meshgrid(xg, yg)

    alpha_vals = slice_df["alpha.metal"].to_numpy()
    temp_vals = slice_df["T"].to_numpy()

    Ag = griddata(points, alpha_vals, (Xg, Yg), method=INTERP_METHOD)
    Tg = griddata(points, temp_vals, (Xg, Yg), method=INTERP_METHOD)

    Ag_near = griddata(points, alpha_vals, (Xg, Yg), method="nearest")
    Tg_near = griddata(points, temp_vals, (Xg, Yg), method="nearest")

    Ag[np.isnan(Ag)] = Ag_near[np.isnan(Ag)]
    Tg[np.isnan(Tg)] = Tg_near[np.isnan(Tg)]

    Tg_liq_region = np.ma.masked_where(Ag <= 0.1, Tg)
    return xg, yg, Xg, Yg, Ag, Tg_liq_region


def find_max_depth(contour_set, surface_y_um: float, contour_name: str):
    segments = [seg for seg in contour_set.allsegs[0] if len(seg) > 0]
    if len(segments) == 0:
        print(f"No {contour_name} contour found.")
        return None

    pts = np.vstack(segments)
    x_contour = pts[:, 0]
    y_contour = pts[:, 1]
    mask = y_contour >= surface_y_um

    if not np.any(mask):
        print(f"No {contour_name} contour found below y = {surface_y_um:.1f} um")
        return None

    x_below = x_contour[mask]
    y_below = y_contour[mask]
    i = np.argmax(y_below)
    x_at_max_depth = x_below[i]
    y_at_max_depth = y_below[i]
    max_depth_um = y_at_max_depth - surface_y_um
    return x_at_max_depth, y_at_max_depth, max_depth_um


def find_surface_width(
    contour_set,
    surface_y_um: float,
    contour_name: str,
    reference_x_um: float | None = None,
    prefer_outermost: bool = True,
):
    intersections = []

    for segment in contour_set.allsegs[0]:
        if len(segment) < 2:
            continue

        x_vals = segment[:, 0]
        y_vals = segment[:, 1]
        for i in range(len(segment) - 1):
            x1, y1 = x_vals[i], y_vals[i]
            x2, y2 = x_vals[i + 1], y_vals[i + 1]

            if y1 == y2:
                continue

            if (y1 - surface_y_um) * (y2 - surface_y_um) > 0:
                continue

            t = (surface_y_um - y1) / (y2 - y1)
            if 0.0 <= t <= 1.0:
                intersections.append(x1 + t * (x2 - x1))

    if len(intersections) < 2:
        print(f"Less than 2 {contour_name} intersections with y = {surface_y_um:.1f} um")
        return None

    intersections = np.array(sorted(intersections))
    unique_intersections = []
    for x_int in intersections:
        if len(unique_intersections) == 0 or abs(x_int - unique_intersections[-1]) > 1e-6:
            unique_intersections.append(float(x_int))

    if len(unique_intersections) < 2:
        print(f"Less than 2 unique {contour_name} intersections with y = {surface_y_um:.1f} um")
        return None

    intervals = [
        (unique_intersections[i], unique_intersections[i + 1])
        for i in range(len(unique_intersections) - 1)
    ]

    if reference_x_um is not None and np.isfinite(reference_x_um):
        containing = [
            interval for interval in intervals
            if interval[0] - 1e-6 <= reference_x_um <= interval[1] + 1e-6
        ]
        if containing:
            x_left, x_right = min(containing, key=lambda interval: interval[1] - interval[0])
            return x_left, x_right, x_right - x_left

        x_left, x_right = min(
            intervals,
            key=lambda interval: min(abs(reference_x_um - interval[0]), abs(reference_x_um - interval[1])),
        )
        return x_left, x_right, x_right - x_left

    if prefer_outermost:
        x_left = min(unique_intersections)
        x_right = max(unique_intersections)
    else:
        x_left, x_right = min(intervals, key=lambda interval: interval[1] - interval[0])

    return x_left, x_right, x_right - x_left


def surface_profile(
    Xg: np.ndarray,
    Yg: np.ndarray,
    values: np.ndarray,
    surface_y_um: float,
) -> tuple[np.ndarray, np.ndarray]:
    x_axis = Xg[0, :]
    y_axis = Yg[:, 0]
    profile = np.array(
        [np.interp(surface_y_um, y_axis, values[:, col]) for col in range(values.shape[1])]
    )
    return x_axis, profile


def find_keyhole_surface_width(
    contour_set,
    surface_y_um: float,
    Xg: np.ndarray,
    Yg: np.ndarray,
    alpha_grid: np.ndarray,
    meltpool_width,
    reference_x_um: float | None = None,
    threshold: float = 0.5,
):
    intersections = []

    for segment in contour_set.allsegs[0]:
        if len(segment) < 2:
            continue

        x_vals = segment[:, 0]
        y_vals = segment[:, 1]
        for i in range(len(segment) - 1):
            x1, y1 = x_vals[i], y_vals[i]
            x2, y2 = x_vals[i + 1], y_vals[i + 1]

            if y1 == y2:
                continue

            if (y1 - surface_y_um) * (y2 - surface_y_um) > 0:
                continue

            t = (surface_y_um - y1) / (y2 - y1)
            if 0.0 <= t <= 1.0:
                intersections.append(x1 + t * (x2 - x1))

    if len(intersections) < 2:
        print(f"Less than 2 alpha.metal = {threshold:g} intersections with y = {surface_y_um:.1f} um")
        return None

    intersections = np.array(sorted(intersections))
    unique_intersections = []
    for x_int in intersections:
        if len(unique_intersections) == 0 or abs(x_int - unique_intersections[-1]) > 1e-6:
            unique_intersections.append(float(x_int))

    if len(unique_intersections) < 2:
        print(f"Less than 2 unique alpha.metal = {threshold:g} intersections with y = {surface_y_um:.1f} um")
        return None

    x_surface, alpha_surface = surface_profile(Xg, Yg, alpha_grid, surface_y_um)
    candidates = []

    for x_left, x_right in zip(unique_intersections[:-1], unique_intersections[1:]):
        width = x_right - x_left
        if width <= 0.0:
            continue

        if meltpool_width is not None:
            mp_left, mp_right, mp_width = meltpool_width
            if width >= mp_width:
                continue
            if x_left < mp_left - 1e-6 or x_right > mp_right + 1e-6:
                continue

        sample_x = np.linspace(x_left, x_right, 65)[1:-1]
        alpha_sample = np.interp(sample_x, x_surface, alpha_surface)
        if len(alpha_sample) == 0 or not np.all(np.isfinite(alpha_sample)):
            continue

        if np.max(alpha_sample) <= threshold:
            candidates.append((x_left, x_right, width))

    if len(candidates) == 0:
        print(
            f"No gas-only keyhole interval found at y = {surface_y_um:.1f} um "
            f"inside the melt-pool width"
        )
        return None

    if reference_x_um is not None and np.isfinite(reference_x_um):
        containing = [
            candidate for candidate in candidates
            if candidate[0] - 1e-6 <= reference_x_um <= candidate[1] + 1e-6
        ]
        if containing:
            return min(containing, key=lambda candidate: candidate[2])

        return min(
            candidates,
            key=lambda candidate: min(
                abs(reference_x_um - candidate[0]),
                abs(reference_x_um - candidate[1]),
            ),
        )

    return min(candidates, key=lambda candidate: candidate[2])


def result_value(result, index: int) -> float:
    if result is None:
        return np.nan
    return result[index]


def plot_slice(ax, Xg, Yg, Ag, Tg_liq_region, title: str, xlabel: str, t_liquidus: float, surface_y_um: float):
    pcm = ax.pcolormesh(
        Xg,
        Yg,
        Ag,
        shading="auto",
        cmap="coolwarm",
        vmin=0,
        vmax=1,
        rasterized=True,
    )
    cs_alpha = ax.contour(
        Xg,
        Yg,
        Ag,
        levels=[0.5],
        colors="black",
        linewidths=2,
    )
    if hasattr(cs_alpha, "collections"):
        for collection in cs_alpha.collections:
            collection.set_rasterized(True)

    cs_tliq = ax.contour(
        Xg,
        Yg,
        Tg_liq_region,
        levels=[t_liquidus],
        colors="yellow",
        linewidths=2,
    )
    if hasattr(cs_tliq, "collections"):
        for collection in cs_tliq.collections:
            collection.set_rasterized(True)

    ax.axhline(
        y=surface_y_um,
        linestyle="--",
        linewidth=1.5,
        color="white",
        zorder=10,
    )
    ax.set_title(title)
    ax.set_xlabel(xlabel)
    ax.set_ylabel("Y / depth (um)")
    ax.axis("equal")
    ax.invert_yaxis()
    return pcm, cs_alpha, cs_tliq


def plot_geometry_history(
    ax,
    rows_by_time: dict[float, dict[str, object]],
    current_row: dict[str, object] | None = None,
) -> None:
    rows = {}
    for time_value, row in rows_by_time.items():
        rows[float(time_value)] = row
    if current_row is not None and np.isfinite(float(current_row["time"])):
        rows[float(current_row["time"])] = current_row

    if len(rows) == 0:
        ax.text(
            0.5,
            0.5,
            "No geometry history yet",
            transform=ax.transAxes,
            ha="center",
            va="center",
            color="black",
        )
        ax.set_axis_off()
        return

    times = np.array(sorted(rows), dtype=float)
    times_us = times * 1e6
    series = [
        ("keyholeDepth_um", "Keyhole depth", "black", "-"),
        ("keyholeWidth_um", "Keyhole width", "0.45", "--"),
        ("meltPoolDepth_um", "Melt-pool depth", "#E3A500", "-"),
        ("meltPoolWidth_um", "Melt-pool width", "#4C78A8", "-"),
    ]

    for field_name, label, color, linestyle in series:
        values = np.array([float(rows[time].get(field_name, np.nan)) for time in times])
        finite = np.isfinite(values)
        if np.any(finite):
            ax.plot(
                times_us[finite],
                values[finite],
                color=color,
                linestyle=linestyle,
                linewidth=1.8,
                marker="o",
                markersize=3.5,
                label=label,
            )

    ax.set_title("Geometry history")
    ax.set_xlabel("Time (us)")
    ax.set_ylabel("Length (um)")
    handles, labels = ax.get_legend_handles_labels()
    if handles:
        ax.legend(handles, labels, frameon=True, facecolor="white", edgecolor="black", fontsize=8)


def combined_rows(rows_by_time: dict[float, dict[str, object]], current_row: dict[str, object]) -> dict[float, dict[str, object]]:
    rows = {float(time_value): row for time_value, row in rows_by_time.items()}
    rows[float(current_row["time"])] = current_row
    return rows


def row_float(row: dict[str, object], field_name: str) -> float:
    try:
        return float(row.get(field_name, np.nan))
    except (TypeError, ValueError):
        return np.nan


def sim_series(
    rows_by_time: dict[float, dict[str, object]],
    current_row: dict[str, object],
    metric: str,
) -> tuple[np.ndarray, np.ndarray]:
    rows = combined_rows(rows_by_time, current_row)
    times = np.array(sorted(rows), dtype=float)

    values = []
    for time_value in times:
        row = rows[time_value]
        if metric == "keyhole_depth":
            values.append(row_float(row, "keyholeDepth_um"))
        elif metric == "keyhole_ar":
            values.append(row_float(row, "keyholeDepth_um") / LASER_DIAMETER_UM)
        elif metric == "meltpool_depth":
            values.append(row_float(row, "meltPoolDepth_um"))
        elif metric == "meltpool_ar":
            depth = row_float(row, "meltPoolDepth_um")
            width = row_float(row, "meltPoolWidth_um")
            values.append(width / depth if depth > 0 else np.nan)
        else:
            raise ValueError(f"Unknown metric: {metric}")

    return times * 1.0e3, np.array(values, dtype=float)


def plot_exp_comparison(
    ax,
    rows_by_time: dict[float, dict[str, object]],
    current_row: dict[str, object],
    exp_csv: Path,
    exp_value_col: str,
    metric: str,
    title: str,
    ylabel: str,
    exp_time_shift_ms: float,
) -> None:
    current_time_ms = float(current_row["time"]) * 1.0e3
    ax.set_title(f"{title} at {current_time_ms:.3f} ms")
    ax.set_xlabel("Time (ms)")
    ax.set_ylabel(ylabel)
    ax.grid(True, alpha=0.3)
    ax.axvline(current_time_ms, color="black", linestyle=":", linewidth=1.2, alpha=0.8)

    sim_time_ms, sim_values = sim_series(rows_by_time, current_row, metric)
    finite = np.isfinite(sim_values)
    if np.any(finite):
        ax.plot(
            sim_time_ms[finite],
            sim_values[finite],
            color="#4C78A8",
            marker="o",
            markersize=3.0,
            linewidth=1.5,
            label="Simulation",
        )

    if exp_csv.exists():
        exp_data = pd.read_csv(exp_csv)
        if "time_ms" in exp_data.columns and exp_value_col in exp_data.columns:
            ax.plot(
                exp_data["time_ms"] - exp_time_shift_ms,
                exp_data[exp_value_col],
                color="#E45756",
                linewidth=1.5,
                label=f"Experiment shifted -{exp_time_shift_ms:.3f} ms",
            )
        else:
            ax.text(
                0.5,
                0.5,
                f"Missing columns in {exp_csv.name}",
                transform=ax.transAxes,
                ha="center",
                va="center",
            )
    else:
        ax.text(
            0.5,
            0.5,
            f"Missing {exp_csv.name}",
            transform=ax.transAxes,
            ha="center",
            va="center",
        )

    handles, labels = ax.get_legend_handles_labels()
    if handles:
        ax.legend(handles, labels, frameon=True, facecolor="white", edgecolor="black", fontsize=8)


def plot_absorptivity_history(ax, absorptivity_csv: Path) -> None:
    ax.set_title("Absorptivity")
    ax.set_xlabel("Time (ms)")
    ax.set_ylabel("Absorptivity")
    ax.set_ylim(0.0, 1.0)
    ax.grid(True, alpha=0.3)

    sim_plotted = False
    if absorptivity_csv.exists():
        try:
            data = pd.read_csv(absorptivity_csv)
            required = {"time_s", "laser", "absorptivity"}
            missing = required.difference(data.columns)
            if not missing:
                for laser_name, laser_data in data.sort_values("time_s").groupby("laser", sort=False):
                    ax.plot(
                        laser_data["time_s"] * 1.0e3,
                        laser_data["absorptivity"],
                        marker="o",
                        markersize=2.5,
                        linewidth=1.3,
                        label=f"Simulation ({laser_name})",
                    )
                sim_plotted = True
        except Exception:
            pass

    exp_csv = None
    curr = Path(__file__).resolve()
    for _ in range(5):
        candidate1 = curr / "exp_abs" / "exp_abs.csv"
        candidate2 = curr / "case_original_thermal_prop" / "exp_abs" / "exp_abs.csv"
        if candidate1.exists():
            exp_csv = candidate1
            break
        elif candidate2.exists():
            exp_csv = candidate2
            break
        curr = curr.parent

    exp_plotted = False
    if exp_csv and exp_csv.exists():
        try:
            exp_data = pd.read_csv(exp_csv)
            if "time_ms" in exp_data.columns and "abs" in exp_data.columns:
                ax.plot(
                    exp_data["time_ms"],
                    exp_data["abs"],
                    color="red",
                    linestyle="--",
                    linewidth=1.5,
                    label="Experiment (Tao Sun)",
                    alpha=0.8,
                )
                exp_plotted = True
        except Exception:
            pass

    if not sim_plotted and not exp_plotted:
        ax.text(
            0.5,
            0.5,
            "Absorptivity data not found",
            transform=ax.transAxes,
            ha="center",
            va="center",
        )
        return

    handles, labels = ax.get_legend_handles_labels()
    if handles:
        ax.legend(handles, labels, frameon=True, facecolor="white", edgecolor="black", fontsize=8)


def plot_target_bar(
    ax,
    title: str,
    sim_value: float,
    exp_value: float,
    ylabel: str,
    time_value_s: float,
) -> None:
    time_ms = time_value_s * 1.0e3
    ax.set_title(f"{title} at {time_ms:.3f} ms")
    ax.set_ylabel(ylabel)
    ax.grid(True, axis="y", alpha=0.3)

    values = [sim_value, exp_value]
    colors = ["#4C78A8", "#E45756"]
    bars = ax.bar(["Simulation", "Experiment"], values, color=colors, width=0.58)

    for bar in bars:
        height = bar.get_height()
        if np.isfinite(height):
            ax.annotate(
                f"{height:.3g}",
                xy=(bar.get_x() + bar.get_width() / 2, height),
                xytext=(0, 3),
                textcoords="offset points",
                ha="center",
                va="bottom",
                fontsize=8,
            )


def plot_target_bars_together(
    ax,
    row: dict[str, object],
    time_value_s: float,
) -> None:
    time_ms = time_value_s * 1.0e3
    sim_width = row_float(row, "meltPoolWidth_um")
    sim_depth = row_float(row, "meltPoolDepth_um")
    sim_ratio = sim_depth / sim_width if sim_width > 0 else np.nan

    sim_values = np.array([sim_width, sim_depth, sim_ratio], dtype=float)
    exp_values = np.array(
        [EXP_TARGET_WIDTH_UM, EXP_TARGET_DEPTH_UM, EXP_TARGET_DEPTH_WIDTH],
        dtype=float,
    )
    normalized_sim = sim_values / exp_values
    normalized_exp = np.ones_like(exp_values)

    labels = ["Width", "Depth", "Depth/width"]
    x = np.arange(len(labels))
    bar_width = 0.36

    ax.set_title(f"Target bars at {time_ms:.3f} ms")
    ax.set_ylabel("Value / target")
    ax.grid(True, axis="y", alpha=0.3)
    ax.axhline(1.0, color="black", linestyle=":", linewidth=1.2)

    sim_bars = ax.bar(
        x - bar_width / 2,
        normalized_sim,
        bar_width,
        color="#4C78A8",
        label="Simulation",
    )
    exp_bars = ax.bar(
        x + bar_width / 2,
        normalized_exp,
        bar_width,
        color="#E45756",
        label="Experiment target",
    )

    ax.set_xticks(x)
    ax.set_xticklabels(labels)

    annotations = [
        (sim_bars[0], sim_width, "um"),
        (sim_bars[1], sim_depth, "um"),
        (sim_bars[2], sim_ratio, ""),
        (exp_bars[0], EXP_TARGET_WIDTH_UM, "um"),
        (exp_bars[1], EXP_TARGET_DEPTH_UM, "um"),
        (exp_bars[2], EXP_TARGET_DEPTH_WIDTH, ""),
    ]
    for bar, value, unit in annotations:
        height = bar.get_height()
        if np.isfinite(height) and np.isfinite(value):
            label = f"{value:.3g}" if unit == "" else f"{value:.3g} {unit}"
            ax.annotate(
                label,
                xy=(bar.get_x() + bar.get_width() / 2, height),
                xytext=(0, 3),
                textcoords="offset points",
                ha="center",
                va="bottom",
                fontsize=8,
            )

    ax.legend(frameon=True, facecolor="white", edgecolor="black", fontsize=8)


def mark_surface_width(ax, width_result, surface_y_um: float, color: str) -> None:
    if width_result is None:
        return

    x_left, x_right, _ = width_result
    ax.plot(
        [x_left, x_right],
        [surface_y_um, surface_y_um],
        color=color,
        linewidth=3.0,
        solid_capstyle="round",
        zorder=23,
    )
    ax.scatter(
        [x_left, x_right],
        [surface_y_um, surface_y_um],
        color=color,
        edgecolor="black",
        s=45,
        zorder=24,
    )


def mark_combined_depths(
    ax,
    keyhole_depth,
    meltpool_depth,
    surface_y_um: float,
    x_grid: np.ndarray,
    label_x_name: str,
    extra_text_lines: list[str] | None = None,
    show_contour_positions: bool = True,
) -> None:
    available_depths = [depth for depth in [keyhole_depth, meltpool_depth] if depth is not None]
    if len(available_depths) == 0:
        return

    x_max = x_grid.max()
    x_span = x_max - x_grid.min()
    y_deepest = max(depth[1] for depth in available_depths)
    x_scale = x_max - 0.08 * x_span
    tick_width = 0.025 * x_span

    ax.annotate(
        "",
        xy=(x_scale, y_deepest),
        xytext=(x_scale, surface_y_um),
        arrowprops=dict(
            arrowstyle="<->",
            color="black",
            linewidth=2.5,
            shrinkA=0,
            shrinkB=0,
        ),
        zorder=20,
    )

    text_lines = []
    if keyhole_depth is not None:
        x_keyhole, y_keyhole, keyhole_depth_um = keyhole_depth
        ax.hlines(
            y_keyhole,
            min(x_keyhole, x_scale),
            max(x_keyhole, x_scale),
            colors="black",
            linestyles=":",
            linewidth=1.5,
            zorder=19,
        )
        ax.hlines(
            y_keyhole,
            x_scale - tick_width,
            x_scale + tick_width,
            colors="black",
            linewidth=2.0,
            zorder=21,
        )
        text_lines.append(f"keyhole_depth = {keyhole_depth_um:.1f} um")
        if show_contour_positions:
            text_lines.append(f"keyhole {label_x_name} = {x_keyhole:.1f} um")

    if meltpool_depth is not None:
        x_meltpool, y_meltpool, meltpool_depth_um = meltpool_depth
        ax.hlines(
            y_meltpool,
            min(x_meltpool, x_scale),
            max(x_meltpool, x_scale),
            colors="gold",
            linestyles=":",
            linewidth=1.5,
            zorder=19,
        )
        ax.hlines(
            y_meltpool,
            x_scale - tick_width,
            x_scale + tick_width,
            colors="gold",
            linewidth=2.0,
            zorder=21,
        )
        text_lines.append(f"meltpool_depth = {meltpool_depth_um:.1f} um")
        if show_contour_positions:
            text_lines.append(f"meltpool {label_x_name} = {x_meltpool:.1f} um")

    if extra_text_lines:
        text_lines.extend(extra_text_lines)

    ax.text(
        0.02,
        0.98,
        "\n".join(text_lines),
        transform=ax.transAxes,
        color="black",
        va="top",
        ha="left",
        bbox=dict(facecolor="white", edgecolor="black", alpha=0.85),
        zorder=22,
    )


def measure(
    vtk_file: Path,
    t_meltpool: float,
    surface_y_um: float,
    output_png: Path,
    rows_by_time: dict[float, dict[str, object]],
    exp_csv_dir: Path,
    absorptivity_csv: Path,
    exp_time_shift_ms: float,
) -> dict[str, float]:
    mesh = pv.read(vtk_file)
    x_min, x_max, y_min, y_max, z_min, z_max = mesh.bounds
    time_value = vtk_time(vtk_file)

    fig = plt.figure(figsize=(18, 14), constrained_layout=True)
    gs = fig.add_gridspec(3, 6, height_ratios=[3.0, 1.65, 1.65])
    slice_axes = [fig.add_subplot(gs[0, 0:3]), fig.add_subplot(gs[0, 3:6])]
    comparison_axes = [
        fig.add_subplot(gs[1, 0:2]),
        fig.add_subplot(gs[1, 2:4]),
        fig.add_subplot(gs[1, 4:6]),
        fig.add_subplot(gs[2, 0:2]),
        fig.add_subplot(gs[2, 2:4]),
        fig.add_subplot(gs[2, 4:6]),
    ]

    # TaoSun scans along x and z is the transverse width direction. Use the
    # center-width XY plane only to choose the streamwise x-location where the
    # transverse YZ measurement slice should be taken.
    z_at_max_depth = 0.5 * (z_min + z_max) * 1e6

    z_plane = np.clip(z_at_max_depth * 1e-6, z_min, z_max)
    xy = slice_to_dataframe(mesh, normal=(0.0, 0.0, 1.0), origin=(0.0, 0.0, z_plane))

    if len(xy) < 3:
        plt.close(fig)
        print(f"Skipping {vtk_file}: XY slice has fewer than 3 cells.")
        return empty_row(vtk_file, t_meltpool, surface_y_um, z_plane * 1e6)

    xg, _, Xg, Yg_xy, Ag_xy, Tg_xy_liq_region = interpolate_slice(xy, "x_um", "y_um", NX, NY)

    _, cs_xy, cs_tliq_xy = plot_slice(
        slice_axes[0],
        Xg,
        Yg_xy,
        Ag_xy,
        Tg_xy_liq_region,
        title=f"XY slice at z = {z_plane * 1e6:.1f} um",
        xlabel="X / scan direction (um)",
        t_liquidus=t_meltpool,
        surface_y_um=surface_y_um,
    )

    xy_depth = find_max_depth(cs_xy, surface_y_um, "alpha.metal = 0.5")
    meltpool_xy_depth = find_max_depth(cs_tliq_xy, surface_y_um, f"T = {t_meltpool:.0f} K")
    if xy_depth is not None:
        x_plane = np.clip(xy_depth[0] * 1e-6, x_min, x_max)
    elif meltpool_xy_depth is not None:
        x_plane = np.clip(meltpool_xy_depth[0] * 1e-6, x_min, x_max)
        print(f"{vtk_file}: using melt-pool deepest x for YZ slice = {x_plane * 1e6:.3f} um")
    else:
        x_plane = 0.5 * (x_min + x_max)
        print(f"{vtk_file}: using fallback x for YZ slice = {x_plane * 1e6:.3f} um")

    mark_combined_depths(
        slice_axes[0],
        xy_depth,
        meltpool_xy_depth,
        surface_y_um,
        xg,
        "x",
        extra_text_lines=[f"YZ slice x = {x_plane * 1e6:.1f} um"],
        show_contour_positions=True,
    )

    yz = slice_to_dataframe(mesh, normal=(1.0, 0.0, 0.0), origin=(x_plane, 0.0, 0.0))

    if len(yz) < 3:
        plt.close(fig)
        print(f"Skipping {vtk_file}: YZ slice has fewer than 3 cells.")
        return empty_row(vtk_file, t_meltpool, surface_y_um, z_plane * 1e6)

    zg, _, Zg, Yg, Ag_yz, Tg_yz_liq_region = interpolate_slice(yz, "z_um", "y_um", NZ, NY)

    pcm, cs_yz, cs_tliq_yz = plot_slice(
        slice_axes[1],
        Zg,
        Yg,
        Ag_yz,
        Tg_yz_liq_region,
        title=f"YZ slice at x = {x_plane * 1e6:.1f} um, t = {time_value:.6g} s",
        xlabel="Z / transverse width (um)",
        t_liquidus=t_meltpool,
        surface_y_um=surface_y_um,
    )

    yz_depth = find_max_depth(cs_yz, surface_y_um, "alpha.metal = 0.5")
    meltpool_yz_depth = find_max_depth(cs_tliq_yz, surface_y_um, f"T = {t_meltpool:.0f} K")
    meltpool_width = find_surface_width(cs_tliq_yz, surface_y_um, f"T = {t_meltpool:.0f} K")
    keyhole_width = find_keyhole_surface_width(
        cs_yz,
        surface_y_um,
        Zg,
        Yg,
        Ag_yz,
        meltpool_width,
        reference_x_um=result_value(yz_depth, 0),
    )
    if (
        keyhole_width is not None
        and meltpool_width is not None
        and keyhole_width[2] > meltpool_width[2]
    ):
        print(
            f"{vtk_file}: rejecting keyhole width {keyhole_width[2]:.1f} um "
            f"because it exceeds melt-pool width {meltpool_width[2]:.1f} um"
        )
        keyhole_width = None

    mark_surface_width(slice_axes[1], meltpool_width, surface_y_um, "gold")
    mark_surface_width(slice_axes[1], keyhole_width, surface_y_um, "black")

    yz_summary_lines = []
    if keyhole_width is not None:
        yz_summary_lines.append(f"keyhole_width = {keyhole_width[2]:.1f} um")
    if meltpool_width is not None:
        yz_summary_lines.append(f"meltpool_width = {meltpool_width[2]:.1f} um")

    mark_combined_depths(
        slice_axes[1],
        yz_depth,
        meltpool_yz_depth,
        surface_y_um,
        zg,
        "z",
        extra_text_lines=yz_summary_lines,
        show_contour_positions=False,
    )

    row = {
        "time": time_value,
        "z_slice_um": z_plane * 1e6,
        "surface_y_um": surface_y_um,
        "T_liq_K": t_meltpool,
        "keyholeDepth_um": result_value(yz_depth, 2),
        "keyholeDepthX_um": result_value(yz_depth, 0),
        "keyholeDepthY_um": result_value(yz_depth, 1),
        "meltPoolDepth_um": result_value(meltpool_yz_depth, 2),
        "meltPoolDepthX_um": result_value(meltpool_yz_depth, 0),
        "meltPoolDepthY_um": result_value(meltpool_yz_depth, 1),
        "keyholeWidth_um": result_value(keyhole_width, 2),
        "keyholeWidthXLeft_um": result_value(keyhole_width, 0),
        "keyholeWidthXRight_um": result_value(keyhole_width, 1),
        "meltPoolWidth_um": result_value(meltpool_width, 2),
        "meltPoolWidthXLeft_um": result_value(meltpool_width, 0),
        "meltPoolWidthXRight_um": result_value(meltpool_width, 1),
    }

    plot_exp_comparison(
        comparison_axes[0],
        rows_by_time,
        row,
        exp_csv_dir / "keyhole_depth_exp.csv",
        "depth_um",
        "keyhole_depth",
        "Keyhole depth",
        "Depth (um)",
        exp_time_shift_ms,
    )
    plot_exp_comparison(
        comparison_axes[1],
        rows_by_time,
        row,
        exp_csv_dir / "keyhole_ar.csv",
        "aspect_ratio_keyholeDepthDividedbyLaserbeamdia80um",
        "keyhole_ar",
        "Keyhole aspect ratio",
        f"Depth / {LASER_DIAMETER_UM:g} um",
        exp_time_shift_ms,
    )
    plot_exp_comparison(
        comparison_axes[2],
        rows_by_time,
        row,
        exp_csv_dir / "meltpool_depth_exp.csv",
        "depth_um",
        "meltpool_depth",
        "Melt-pool depth",
        "Depth (um)",
        exp_time_shift_ms,
    )
    plot_exp_comparison(
        comparison_axes[3],
        rows_by_time,
        row,
        exp_csv_dir / "meltpool_ar_exp.csv",
        "aspect_ratio_meltpoolWidthbyDepth",
        "meltpool_ar",
        "Melt-pool aspect ratio",
        "Width / depth",
        exp_time_shift_ms,
    )
    plot_absorptivity_history(comparison_axes[4], absorptivity_csv)
    plot_target_bars_together(comparison_axes[5], row, time_value)

    fig.colorbar(pcm, ax=slice_axes, label="alpha.metal")
    slice_axes[0].legend(
        handles=legend_handles(t_meltpool),
        loc="upper right",
        frameon=True,
        facecolor="white",
        edgecolor="black",
        framealpha=0.9,
    )

    output_png.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_png, format="png", bbox_inches="tight", dpi=PNG_DPI)
    plt.close(fig)

    return row


def empty_row(
    vtk_file: Path,
    t_meltpool: float,
    surface_y_um: float,
    z_slice_um: float = np.nan,
) -> dict[str, float]:
    row = {field_name: np.nan for field_name in CSV_FIELDS}
    row["time"] = vtk_time(vtk_file)
    row["z_slice_um"] = z_slice_um
    row["surface_y_um"] = surface_y_um
    row["T_liq_K"] = t_meltpool
    return row


def legend_handles(t_liquidus: float) -> list[Line2D]:
    return [
        Line2D([0], [0], color="black", linewidth=2, label="alpha.metal = 0.5"),
        Line2D([0], [0], color="yellow", linewidth=2, label=f"T = {t_liquidus:.0f} K"),
        Line2D([0], [0], color="white", linestyle="--", linewidth=1.5, label="substrate surface"),
        Line2D(
            [0],
            [0],
            color="black",
            linewidth=2.5,
            marker=r"$\leftrightarrow$",
            markersize=14,
            label="depth scale",
        ),
        Line2D([0], [0], color="gold", linewidth=2.0, label="melt pool depth tick"),
        Line2D([0], [0], color="black", linewidth=3.0, marker="o", label="keyhole width"),
        Line2D([0], [0], color="gold", linewidth=3.0, marker="o", label="melt pool width"),
    ]


def read_existing_rows(csv_path: Path) -> dict[float, dict[str, str]]:
    if not csv_path.exists():
        return {}

    rows = {}
    with csv_path.open(newline="") as handle:
        reader = csv.DictReader(handle)
        for row in reader:
            try:
                rows[float(row["time"])] = row
            except (KeyError, TypeError, ValueError):
                continue
    return rows


def write_rows(csv_path: Path, rows_by_time: dict[float, dict[str, object]]) -> None:
    with csv_path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=CSV_FIELDS)
        writer.writeheader()
        for time_value in sorted(rows_by_time):
            writer.writerow(rows_by_time[time_value])


def section_index(vtk_file: Path, vtk_dir: Path) -> int:
    vtk_files = converted_vtk_files(vtk_dir)
    vtk_file = vtk_file.resolve()
    for index, candidate in enumerate(vtk_files):
        if candidate.resolve() == vtk_file:
            return index
    return sum(vtk_time(candidate) < vtk_time(vtk_file) for candidate in vtk_files)


def numbered_section_png(image_dir: Path, index: int, vtk_file: Path) -> Path:
    name = vtk_file.parent.name if vtk_file.name in {"internal.vtu", "internal.vtk"} else vtk_file.stem
    return image_dir / f"{index:04d}_{name}.png"


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--case", type=Path, default=Path(__file__).resolve().parents[1])
    parser.add_argument("--vtk-dir", type=Path, default=None)
    parser.add_argument("--vtk-file", type=Path, default=None)
    parser.add_argument("--surface-y-um", type=float, default=SURFACE_Y_UM)
    parser.add_argument("--t-liquidus", type=float, default=None)
    parser.add_argument("--t-meltpool", type=float, default=None)
    parser.add_argument("--exp-csv-dir", type=Path, default=TAOSUN_EXP_CSV_DIR)
    parser.add_argument("--absorptivity-csv", type=Path, default=ABSORPTIVITY_CSV)
    parser.add_argument("--exp-time-shift-ms", type=float, default=EXP_TIME_SHIFT_MS)
    args = parser.parse_args()

    case = args.case.resolve()
    surface_y_um = args.surface_y_um
    if surface_y_um is None:
        surface_y_um = infer_surface_y_um(case)
        print(f"Inferred surface y = {surface_y_um:.6g} um from system/setFieldsDict")

    vtk_dir = args.vtk_dir or case / "VTK"
    out_dir = case / "post-processing-data"
    image_dir = out_dir / "vtu_sections"
    out_dir.mkdir(exist_ok=True)
    image_dir.mkdir(exist_ok=True)

    t_meltpool = args.t_meltpool
    if t_meltpool is None:
        if args.t_liquidus is not None:
            t_meltpool = args.t_liquidus
        else:
            t_meltpool = parse_scalar(case / "constant" / "transportProperties", "Tsolidus")

    if args.vtk_file:
        vtk_files = [args.vtk_file.resolve()]
    else:
        vtk_files = converted_vtk_files(vtk_dir)

    csv_path = out_dir / "vtu_meltpool_geometry.csv"
    rows_by_time = read_existing_rows(csv_path)

    for vtk_file in vtk_files:
        print(f"Analyzing {vtk_file}")
        index = section_index(vtk_file, vtk_dir)
        row = measure(
            vtk_file,
            t_meltpool,
            surface_y_um,
            numbered_section_png(image_dir, index, vtk_file),
            rows_by_time,
            args.exp_csv_dir.resolve(),
            args.absorptivity_csv.resolve(),
            args.exp_time_shift_ms,
        )
        rows_by_time[float(row["time"])] = row

    write_rows(csv_path, rows_by_time)

    print(f"Wrote {csv_path}")
    print(f"Wrote PNG sections to {image_dir}")


if __name__ == "__main__":
    main()
