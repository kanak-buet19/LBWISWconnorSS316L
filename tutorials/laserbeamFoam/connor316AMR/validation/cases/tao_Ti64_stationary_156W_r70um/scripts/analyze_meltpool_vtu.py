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
SURFACE_Y_UM = 128
PNG_DPI = 150
EXP_SUMMARY_CSV = Path(__file__).resolve().parents[1] / "exp_tao_Ti64_stationary_156W_r70um_summary.csv"
EXP_MASK_IMAGE = Path(__file__).resolve().parents[1] / "exp_tao_Ti64_stationary_156W_r70um_mask.png"
COMPARE_DEPTH_FIELD = "keyholeDepth_um"


CSV_FIELDS = [
    "time",
    "z_slice_um",
    "surface_y_um",
    "T_threshold_K",
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


def vtk_time(path: Path) -> float:
    with path.open("rb") as handle:
        header = handle.read(512).decode("utf-8", errors="ignore")

    match = re.search(r"time='([0-9.eE+-]+)'", header)
    if match:
        return float(match.group(1))

    match = re.search(r"_([0-9.eE+-]+)\.(?:vtk|vtu)$", path.name)
    if not match:
        match = re.search(r"_([0-9.eE+-]+)$", path.parent.name)
    if not match:
        raise ValueError(f"Cannot parse time from {path}")
    return float(match.group(1))


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

    Tg_threshold_region = np.ma.masked_where(Ag <= 0.1, Tg)
    return xg, yg, Xg, Yg, Ag, Tg_threshold_region


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


def read_experimental_meltpool(exp_csv: Path) -> dict[str, float] | None:
    if not exp_csv.exists():
        print(f"Experimental melt-pool summary not found: {exp_csv}")
        return None

    df = pd.read_csv(exp_csv)
    if df.empty:
        print(f"Experimental melt-pool summary is empty: {exp_csv}")
        return None

    row = df.iloc[0]
    result = {"depth_um": float(row["depth_um"])}
    if "width_um" in df.columns:
        result["width_um"] = float(row["width_um"])
    for field_name in ("um_per_px", "left_x_px", "right_x_px", "surface_y_px", "bottom_y_px"):
        if field_name in row:
            result[field_name] = float(row[field_name])
    return result


def plot_slice(ax, Xg, Yg, Ag, Tg_threshold_region, title: str, xlabel: str, t_threshold: float, surface_y_um: float):
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

    cs_t = ax.contour(
        Xg,
        Yg,
        Tg_threshold_region,
        levels=[t_threshold],
        colors="yellow",
        linewidths=2,
    )
    if hasattr(cs_t, "collections"):
        for collection in cs_t.collections:
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
    return pcm, cs_alpha, cs_t


def plot_meltpool_comparison(
    ax,
    sim_width_um: float,
    sim_depth_um: float,
    exp_metrics: dict[str, float] | None,
) -> None:
    title = "Keyhole geometry comparison" if "keyhole" in COMPARE_DEPTH_FIELD.lower() else "Melt-pool geometry comparison"
    ax.set_title(title)

    if exp_metrics is None:
        ax.text(
            0.5,
            0.5,
            "Experimental summary CSV not found",
            transform=ax.transAxes,
            ha="center",
            va="center",
            color="black",
        )
        ax.set_axis_off()
        return

    if "width_um" in exp_metrics:
        categories = ["Width", "Depth"]
        sim_values = [sim_width_um, sim_depth_um]
        exp_values = [exp_metrics["width_um"], exp_metrics["depth_um"]]
    else:
        categories = ["Depth"]
        sim_values = [sim_depth_um]
        exp_values = [exp_metrics["depth_um"]]

    x = np.arange(len(categories))
    bar_width = 0.34
    sim_bars = ax.bar(x - bar_width / 2, sim_values, bar_width, label="Simulation", color="#4C78A8")
    exp_bars = ax.bar(x + bar_width / 2, exp_values, bar_width, label="Experiment", color="#F58518")

    ax.set_xticks(x)
    ax.set_xticklabels(categories)
    ax.set_ylabel("Length (um)")
    ax.legend(frameon=True, facecolor="white", edgecolor="black")

    for bars in (sim_bars, exp_bars):
        for bar in bars:
            height = bar.get_height()
            if np.isfinite(height):
                ax.annotate(
                    f"{height:.1f}",
                    xy=(bar.get_x() + bar.get_width() / 2, height),
                    xytext=(0, 3),
                    textcoords="offset points",
                    ha="center",
                    va="bottom",
                    fontsize=8,
                )

def plot_experimental_mask(
    ax,
    exp_metrics: dict[str, float] | None,
    exp_mask_image: Path | None,
    surface_y_um: float,
    x_center_um: float,
    x_limits: tuple[float, float],
    y_limits: tuple[float, float],
) -> None:
    ax.set_title("Experiment mask")
    ax.set_xlabel("X / width (um)")
    ax.set_ylabel("Y / depth (um)")

    if exp_metrics is None or exp_mask_image is None or not exp_mask_image.exists():
        ax.text(
            0.5,
            0.5,
            "Experimental mask not found",
            transform=ax.transAxes,
            ha="center",
            va="center",
            color="black",
        )
        ax.set_axis_off()
        return

    image = plt.imread(exp_mask_image)
    height_px, width_px = image.shape[:2]
    um_per_px = float(exp_metrics.get("um_per_px", 1.0))
    x_mid_px = 0.5 * (
        float(exp_metrics.get("left_x_px", 0.0))
        + float(exp_metrics.get("right_x_px", width_px))
    )
    surface_y_px = float(exp_metrics.get("surface_y_px", 0.0))

    x_min = x_center_um + (0.0 - x_mid_px) * um_per_px
    x_max = x_center_um + (width_px - x_mid_px) * um_per_px
    y_top = surface_y_um + (0.0 - surface_y_px) * um_per_px
    y_bottom = surface_y_um + (height_px - surface_y_px) * um_per_px

    ax.imshow(image, extent=[x_min, x_max, y_bottom, y_top], origin="upper")
    ax.axhline(surface_y_um, linestyle="--", linewidth=1.5, color="white", zorder=10)
    ax.set_xlim(*x_limits)
    ax.set_ylim(*y_limits)
    ax.set_aspect("equal", adjustable="box")


def _read_absorptivity(case_dir: Path) -> pd.DataFrame | None:
    csv = case_dir / "absorptivity_vs_time" / "absorptivity_vs_time.csv"
    if not csv.exists():
        return None
    df = pd.read_csv(csv)
    if not {"time_s", "absorptivity"}.issubset(df.columns):
        return None
    if "laser" not in df.columns:
        df["laser"] = "laser0"
    return df.sort_values("time_s")


def _plot_absorptivity(ax, absorptivity_df: pd.DataFrame | None) -> None:
    ax.set_title("Absorptivity vs time")
    ax.set_xlabel("Time (us)")
    ax.set_ylabel("Absorptivity")
    ax.set_ylim(0, 1)

    if absorptivity_df is None or absorptivity_df.empty:
        ax.text(0.5, 0.5, "No absorptivity data", transform=ax.transAxes, ha="center", va="center")
        ax.set_axis_off()
        return

    single_laser = absorptivity_df["laser"].nunique() == 1
    for laser_name, grp in absorptivity_df.groupby("laser", sort=False):
        lbl = "absorptivity" if single_laser else f"absorptivity ({laser_name})"
        ax.plot(grp["time_s"] * 1e6, grp["absorptivity"], color="#2CA02C", linewidth=1.5, label=lbl)

    ax.legend(frameon=True, facecolor="white", edgecolor="black", fontsize=8)
    ax.grid(True, alpha=0.3)


def plot_keyhole_timeseries(
    ax,
    rows_by_time: dict[float, dict[str, object]],
    current_row: dict[str, object] | None,
    exp_timeseries_df,
    laser_off_s: float | None = None,
) -> None:
    rows = dict(rows_by_time)
    if current_row is not None and np.isfinite(float(current_row["time"])):
        rows[float(current_row["time"])] = current_row

    if exp_timeseries_df is not None:
        ax.plot(
            exp_timeseries_df["t_ms"],
            exp_timeseries_df["keyhole_depth_um"],
            color="#E3A500",
            linewidth=1.5,
            label="Experiment",
            zorder=2,
        )

    if rows:
        times = np.array(sorted(rows), dtype=float)
        values = np.array([float(rows[t].get("keyholeDepth_um", np.nan)) for t in times])
        finite = np.isfinite(values)
        if np.any(finite):
            ax.plot(
                times[finite] * 1e3,
                values[finite],
                color="black",
                linewidth=1.8,
                marker="o",
                markersize=3.5,
                label="Simulation",
                zorder=3,
            )

    if laser_off_s is not None:
        ax.axvline(
            laser_off_s * 1e3,
            color="red",
            linestyle="--",
            linewidth=1.2,
            label=f"Laser off ({laser_off_s * 1e3:.0f} ms)",
        )

    ax.set_xlabel("Time (ms)")
    ax.set_ylabel("Keyhole depth (µm)")
    ax.set_title("Keyhole depth vs time — exp vs sim")
    ax.legend(frameon=True, fontsize=8)
    ax.grid(True, linewidth=0.5, alpha=0.4)


def plot_geometry_history(
    ax,
    rows_by_time: dict[float, dict[str, object]],
    current_row: dict[str, object] | None = None,
    exp_metrics: dict[str, float] | None = None,
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

    if exp_metrics is not None and "depth_um" in exp_metrics:
        ax.axhline(
            exp_metrics["depth_um"],
            color="red",
            linestyle=":",
            linewidth=1.5,
            label=f"Exp. mean depth = {exp_metrics['depth_um']:.0f} µm",
        )
    depth_values = np.array([float(rows[t].get(COMPARE_DEPTH_FIELD, np.nan)) for t in times], dtype=float)
    mean_depth = float(np.nanmean(depth_values)) if np.any(np.isfinite(depth_values)) else np.nan
    if np.isfinite(mean_depth):
        ax.axhline(
            mean_depth,
            color="navy",
            linestyle="--",
            linewidth=1.2,
            label=f"Sim. mean depth = {mean_depth:.0f} µm",
        )
    ax.set_title("Geometry history")
    ax.set_xlabel("Time (us)")
    ax.set_ylabel("Length (um)")

    handles, labels = ax.get_legend_handles_labels()
    if handles:
        ax.legend(handles, labels, frameon=True, facecolor="white", edgecolor="black", fontsize=8)


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
    t_threshold: float,
    surface_y_um: float,
    output_png: Path,
    exp_metrics: dict[str, float] | None,
    exp_mask_image: Path | None,
    rows_by_time: dict[float, dict[str, object]],
    case_dir: Path | None = None,
    exp_timeseries_df=None,
    laser_off_s: float | None = None,
) -> dict[str, float]:
    mesh = pv.read(vtk_file)
    x_min, x_max, y_min, y_max, z_min, z_max = mesh.bounds

    x_plane = 0.5 * (x_min + x_max)
    yz = slice_to_dataframe(mesh, normal=(1.0, 0.0, 0.0), origin=(x_plane, 0.0, 0.0))

    if len(yz) < 3:
        print(f"Skipping {vtk_file}: YZ slice has fewer than 3 cells.")
        return empty_row(vtk_file, t_threshold, surface_y_um)

    zg, _, Zg, Yg, Ag_yz, Tg_yz_liq_region = interpolate_slice(yz, "z_um", "y_um", NZ, NY)

    has_timeseries = exp_timeseries_df is not None
    nrows = 3 if has_timeseries else 2
    height_ratios = [3.0, 1.9, 1.9] if has_timeseries else [3.0, 1.9]
    fig = plt.figure(figsize=(20, 14 if has_timeseries else 11), constrained_layout=True)
    gs = fig.add_gridspec(nrows, 6, height_ratios=height_ratios)
    slice_axes = [fig.add_subplot(gs[0, 0:2]), fig.add_subplot(gs[0, 2:4])]
    absorptivity_ax = fig.add_subplot(gs[0, 4:6])
    comparison_ax = fig.add_subplot(gs[1, 0:2])
    exp_mask_ax = fig.add_subplot(gs[1, 2:4])
    history_ax = fig.add_subplot(gs[1, 4:6])
    timeseries_ax = fig.add_subplot(gs[2, 0:6]) if has_timeseries else None
    time_value = vtk_time(vtk_file)

    pcm, cs_yz, cs_t_yz = plot_slice(
        slice_axes[0],
        Zg,
        Yg,
        Ag_yz,
        Tg_yz_liq_region,
        title=f"YZ slice at x = {x_plane * 1e6:.1f} um, t = {time_value:.6g} s",
        xlabel="Z / scan track (um)",
        t_threshold=t_threshold,
        surface_y_um=surface_y_um,
    )

    yz_depth = find_max_depth(cs_yz, surface_y_um, "alpha.metal = 0.5")
    meltpool_yz_depth = find_max_depth(cs_t_yz, surface_y_um, f"T = {t_threshold:.0f} K")
    mark_combined_depths(slice_axes[0], yz_depth, meltpool_yz_depth, surface_y_um, zg, "z")

    if yz_depth is None:
        z_at_max_depth = yz["z_um"].iloc[len(yz) // 2]
        print(f"{vtk_file}: using fallback z for XY slice = {z_at_max_depth:.3f} um")
    else:
        z_at_max_depth = yz_depth[0]

    z_plane = np.clip(z_at_max_depth * 1e-6, z_min, z_max)
    xy = slice_to_dataframe(mesh, normal=(0.0, 0.0, 1.0), origin=(0.0, 0.0, z_plane))

    if len(xy) < 3:
        plt.close(fig)
        print(f"Skipping {vtk_file}: XY slice has fewer than 3 cells.")
        return empty_row(vtk_file, t_threshold, surface_y_um, z_plane * 1e6)

    xg, _, Xg, Yg_xy, Ag_xy, Tg_xy_liq_region = interpolate_slice(xy, "x_um", "y_um", NX, NY)

    _, cs_xy, cs_t_xy = plot_slice(
        slice_axes[1],
        Xg,
        Yg_xy,
        Ag_xy,
        Tg_xy_liq_region,
        title=f"XY slice at z = {z_plane * 1e6:.1f} um",
        xlabel="X / width (um)",
        t_threshold=t_threshold,
        surface_y_um=surface_y_um,
    )

    xy_depth = find_max_depth(cs_xy, surface_y_um, "alpha.metal = 0.5")
    meltpool_xy_depth = find_max_depth(cs_t_xy, surface_y_um, f"T = {t_threshold:.0f} K")
    meltpool_width = find_surface_width(cs_t_xy, surface_y_um, f"T = {t_threshold:.0f} K")
    keyhole_width = find_keyhole_surface_width(
        cs_xy,
        surface_y_um,
        Xg,
        Yg_xy,
        Ag_xy,
        meltpool_width,
        reference_x_um=result_value(xy_depth, 0),
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

    xy_summary_lines = []
    if keyhole_width is not None:
        xy_summary_lines.append(f"keyhole_width = {keyhole_width[2]:.1f} um")
    if meltpool_width is not None:
        xy_summary_lines.append(f"meltpool_width = {meltpool_width[2]:.1f} um")

    mark_combined_depths(
        slice_axes[1],
        xy_depth,
        meltpool_xy_depth,
        surface_y_um,
        xg,
        "x",
        extra_text_lines=xy_summary_lines,
        show_contour_positions=False,
    )

    row = {
        "time": time_value,
        "z_slice_um": z_plane * 1e6,
        "surface_y_um": surface_y_um,
        "T_threshold_K": t_threshold,
        "keyholeDepth_um": result_value(xy_depth, 2),
        "keyholeDepthX_um": result_value(xy_depth, 0),
        "keyholeDepthY_um": result_value(xy_depth, 1),
        "meltPoolDepth_um": result_value(meltpool_xy_depth, 2),
        "meltPoolDepthX_um": result_value(meltpool_xy_depth, 0),
        "meltPoolDepthY_um": result_value(meltpool_xy_depth, 1),
        "keyholeWidth_um": result_value(keyhole_width, 2),
        "keyholeWidthXLeft_um": result_value(keyhole_width, 0),
        "keyholeWidthXRight_um": result_value(keyhole_width, 1),
        "meltPoolWidth_um": result_value(meltpool_width, 2),
        "meltPoolWidthXLeft_um": result_value(meltpool_width, 0),
        "meltPoolWidthXRight_um": result_value(meltpool_width, 1),
    }

    plot_meltpool_comparison(
        comparison_ax,
        row["meltPoolWidth_um"],
        row[COMPARE_DEPTH_FIELD],
        exp_metrics,
    )
    plot_experimental_mask(
        exp_mask_ax,
        exp_metrics,
        exp_mask_image,
        surface_y_um,
        0.5 * (float(Xg.min()) + float(Xg.max())),
        (float(Xg.min()), float(Xg.max())),
        (float(Yg_xy.max()), float(Yg_xy.min())),
    )
    absorptivity_df = _read_absorptivity(case_dir) if case_dir is not None else None
    _plot_absorptivity(absorptivity_ax, absorptivity_df)
    plot_geometry_history(history_ax, rows_by_time, row, exp_metrics)
    if timeseries_ax is not None:
        plot_keyhole_timeseries(timeseries_ax, rows_by_time, row, exp_timeseries_df, laser_off_s)

    fig.colorbar(pcm, ax=slice_axes, label="alpha.metal")
    slice_axes[1].legend(
        handles=legend_handles(t_threshold),
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
    t_threshold: float,
    surface_y_um: float,
    z_slice_um: float = np.nan,
) -> dict[str, float]:
    row = {field_name: np.nan for field_name in CSV_FIELDS}
    row["time"] = vtk_time(vtk_file)
    row["z_slice_um"] = z_slice_um
    row["surface_y_um"] = surface_y_um
    row["T_threshold_K"] = t_threshold
    return row


def legend_handles(t_threshold: float) -> list[Line2D]:
    return [
        Line2D([0], [0], color="black", linewidth=2, label="alpha.metal = 0.5"),
        Line2D([0], [0], color="yellow", linewidth=2, label=f"T = {t_threshold:.0f} K"),
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
        writer = csv.DictWriter(handle, fieldnames=CSV_FIELDS, extrasaction="ignore")
        writer.writeheader()
        for time_value in sorted(rows_by_time):
            writer.writerow(rows_by_time[time_value])


def section_index(vtk_file: Path, vtk_dir: Path) -> int:
    vtk_files = sorted(vtk_dir.glob("*/internal.vtu"), key=vtk_time)
    vtk_file = vtk_file.resolve()
    for index, candidate in enumerate(vtk_files):
        if candidate.resolve() == vtk_file:
            return index
    return sum(vtk_time(candidate) < vtk_time(vtk_file) for candidate in vtk_files)


def numbered_section_png(png_dir: Path, index: int, vtk_file: Path) -> Path:
    name = vtk_file.parent.name if vtk_file.name == "internal.vtu" else vtk_file.stem
    return png_dir / f"{index:04d}_{name}.png"


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--case", type=Path, default=Path(__file__).resolve().parents[1])
    parser.add_argument("--vtk-dir", type=Path, default=None)
    parser.add_argument("--vtk-file", type=Path, default=None)
    parser.add_argument("--surface-y-um", type=float, default=SURFACE_Y_UM)
    parser.add_argument("--t-threshold", type=float, default=None,
                        help="Temperature threshold for melt-pool contour; defaults to Tsolidus.")
    parser.add_argument("--t-liquidus", type=float, default=None,
                        help="Backward-compatible alias for --t-threshold.")
    parser.add_argument("--exp-summary-csv", type=Path, default=EXP_SUMMARY_CSV)
    parser.add_argument("--exp-mask-image", type=Path, default=EXP_MASK_IMAGE)
    parser.add_argument("--exp-timeseries-csv", type=Path, default=None,
                        help="Time-series exp CSV (t_ms, keyhole_depth_um). "
                             "Falls back to exp_timeseries_csv in case_info.json.")
    args = parser.parse_args()

    case = args.case.resolve()
    vtk_dir = args.vtk_dir or case / "VTK"
    out_dir = case / "post-processing-data"
    png_dir = out_dir / "vtu_sections"
    out_dir.mkdir(exist_ok=True)
    png_dir.mkdir(exist_ok=True)

    t_threshold = args.t_threshold if args.t_threshold is not None else args.t_liquidus
    if t_threshold is None:
        t_threshold = parse_scalar(case / "constant" / "transportProperties", "Tsolidus")

    if args.vtk_file:
        vtk_files = [args.vtk_file.resolve()]
    else:
        vtk_files = sorted(vtk_dir.glob("*/internal.vtu"), key=vtk_time)

    exp_metrics = read_experimental_meltpool(args.exp_summary_csv.resolve())
    exp_mask_image = args.exp_mask_image.resolve()

    case_info = {}
    case_info_path = case / "case_info.json"
    if case_info_path.exists():
        case_info = json.loads(case_info_path.read_text(encoding="utf-8"))

    laser_off_s = case_info.get("laser_off_time_s")

    exp_ts_path = args.exp_timeseries_csv
    if exp_ts_path is None and "exp_timeseries_csv" in case_info:
        exp_ts_path = case.parents[1] / case_info["exp_timeseries_csv"]

    exp_timeseries_df = None
    if exp_ts_path is not None:
        exp_ts_path = Path(exp_ts_path).resolve()
        if exp_ts_path.exists():
            _df = pd.read_csv(exp_ts_path)
            if {"t_ms", "keyhole_depth_um"}.issubset(_df.columns):
                exp_timeseries_df = _df[["t_ms", "keyhole_depth_um"]].dropna().sort_values("t_ms")
            else:
                print(f"exp-timeseries-csv missing required columns (t_ms, keyhole_depth_um): {exp_ts_path}")
        else:
            print(f"exp-timeseries-csv not found: {exp_ts_path}")

    csv_path = out_dir / "vtu_meltpool_geometry.csv"
    rows_by_time = read_existing_rows(csv_path)

    for vtk_file in vtk_files:
        print(f"Analyzing {vtk_file}")
        index = section_index(vtk_file, vtk_dir)
        row = measure(
            vtk_file,
            t_threshold,
            args.surface_y_um,
            numbered_section_png(png_dir, index, vtk_file),
            exp_metrics,
            exp_mask_image,
            rows_by_time,
            case_dir=case,
            exp_timeseries_df=exp_timeseries_df,
            laser_off_s=laser_off_s,
        )
        rows_by_time[float(row["time"])] = row

    write_rows(csv_path, rows_by_time)

    print(f"Wrote {csv_path}")
    print(f"Wrote PNG sections to {png_dir}")


if __name__ == "__main__":
    main()
