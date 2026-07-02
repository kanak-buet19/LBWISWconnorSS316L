#!/usr/bin/env python3
"""Measure melt-pool/keyhole metrics from reconstructed legacy VTK files.

Outputs:
  post-processing-data/vtk_meltpool_geometry.csv
  post-processing-data/vtk_sections/*.png
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
SURFACE_Y_UM = 120
PNG_DPI = 150
EXP_SUMMARY_CSV = Path(__file__).resolve().parents[1] / "exp_hofmann_200W_1200_summary.csv"
EXP_MASK_IMAGE = Path(__file__).resolve().parents[1] / "exp_hofmann_200W_1200_mask.png"
COMPARE_DEPTH_FIELD = "meltPoolDepth_um"
KEYHOLE_DEPRESSION_TOLERANCE_UM = 10.0
SECTION_CANDIDATE_COUNT = 41


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

    match = re.search(r"_([0-9.eE+-]+)\.vtk$", path.name)
    if not match:
        match = re.search(r"_([0-9.eE+-]+)$", path.parent.name)
    if not match:
        raise ValueError(f"Cannot parse time from {path}")
    return float(match.group(1))


def internal_mesh_vtk_files(vtk_dir: Path) -> list[Path]:
    return sorted(vtk_dir.glob("*.vtk"), key=vtk_time)


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


def _plot_absorptivity(
    ax,
    absorptivity_df: pd.DataFrame | None,
    exp_timeseries_df: pd.DataFrame | None = None,
) -> None:
    ax.set_title("Absorptivity vs time")
    ax.set_xlabel("Time (ms)")
    ax.set_ylabel("Absorptivity")
    ax.set_ylim(0, 1)

    has_sim = absorptivity_df is not None and not absorptivity_df.empty
    has_exp = (
        exp_timeseries_df is not None
        and "absorptance" in exp_timeseries_df.columns
        and exp_timeseries_df["absorptance"].notna().any()
    )

    if not has_sim and not has_exp:
        ax.text(0.5, 0.5, "No absorptivity data", transform=ax.transAxes, ha="center", va="center")
        ax.set_axis_off()
        return

    if has_exp:
        exp_abs = exp_timeseries_df.dropna(subset=["t_ms", "absorptance"])
        ax.plot(
            exp_abs["t_ms"],
            exp_abs["absorptance"],
            color="#D62728",
            linewidth=1.5,
            label="Exp. absorptance",
        )

    if has_sim:
        single_laser = absorptivity_df["laser"].nunique() == 1
        for laser_name, grp in absorptivity_df.groupby("laser", sort=False):
            lbl = "Sim. absorptivity" if single_laser else f"Sim. absorptivity ({laser_name})"
            ax.plot(grp["time_s"] * 1e3, grp["absorptivity"], color="#2CA02C", linewidth=1.5, label=lbl)

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
        exp_depth = exp_timeseries_df.dropna(subset=["t_ms", "keyhole_depth_um"])
        if "keyhole_depth_std_um" in exp_depth.columns:
            exp_std = exp_depth["keyhole_depth_std_um"].fillna(0.0)
            ax.fill_between(
                exp_depth["t_ms"],
                exp_depth["keyhole_depth_um"] - exp_std,
                exp_depth["keyhole_depth_um"] + exp_std,
                color="#E3A500",
                alpha=0.2,
                linewidth=0,
                label="Experiment ±1σ",
                zorder=1,
            )
        ax.plot(
            exp_depth["t_ms"],
            exp_depth["keyhole_depth_um"],
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
    mark_combined_depths(slice_axes[0], yz_depth, None, surface_y_um, zg, "z")

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

    mark_surface_width(slice_axes[1], keyhole_width, surface_y_um, "black")

    xy_summary_lines = []
    if keyhole_width is not None:
        xy_summary_lines.append(f"keyhole_width = {keyhole_width[2]:.1f} um")
    mark_combined_depths(
        slice_axes[1],
        xy_depth,
        None,
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
        row["keyholeWidth_um"],
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
    _plot_absorptivity(absorptivity_ax, absorptivity_df, exp_timeseries_df)
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
        Line2D([0], [0], color="black", linewidth=3.0, marker="o", label="keyhole width"),
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
    vtk_files = internal_mesh_vtk_files(vtk_dir)
    vtk_file = vtk_file.resolve()
    for index, candidate in enumerate(vtk_files):
        if candidate.resolve() == vtk_file:
            return index
    return sum(vtk_time(candidate) < vtk_time(vtk_file) for candidate in vtk_files)


def numbered_section_png(png_dir: Path, index: int, vtk_file: Path) -> Path:
    return png_dir / f"{index:04d}_{vtk_file.stem}.png"


def read_laser_z_table(case: Path) -> tuple[np.ndarray, np.ndarray] | None:
    table_path = case / "constant" / "timeVsLaserPosition"
    if not table_path.exists():
        return None

    text = table_path.read_text(encoding="utf-8")
    rows = []
    for match in re.finditer(
        r"\(\s*([0-9.eE+-]+)\s+\(\s*([0-9.eE+-]+)\s+([0-9.eE+-]+)\s+([0-9.eE+-]+)\s*\)\s*\)",
        text,
    ):
        rows.append((float(match.group(1)), float(match.group(4)) * 1e6))

    if len(rows) < 2:
        return None

    rows.sort()
    times = np.array([row[0] for row in rows], dtype=float)
    z_um = np.array([row[1] for row in rows], dtype=float)
    return times, z_um


def scan_z_from_time(case: Path, times: np.ndarray) -> np.ndarray:
    table = read_laser_z_table(case)
    if table is not None:
        table_times, table_z_um = table
        return np.interp(times, table_times, table_z_um)

    try:
        v_scan = parse_scalar(case / "constant" / "transportProperties", "V_scan")
    except (FileNotFoundError, ValueError):
        v_scan = 0.0
    return times * v_scan * 1e6


def scan_time_from_z(case: Path, z_um: np.ndarray) -> np.ndarray:
    table = read_laser_z_table(case)
    if table is not None:
        table_times, table_z_um = table
        order = np.argsort(table_z_um)
        return np.interp(z_um, table_z_um[order], table_times[order])

    try:
        v_scan = parse_scalar(case / "constant" / "transportProperties", "V_scan")
    except (FileNotFoundError, ValueError):
        v_scan = 0.0
    if v_scan <= 0.0:
        return np.zeros_like(z_um, dtype=float)
    return z_um * 1e-6 / v_scan


def stable_keyhole_sections(
    rows_by_time: dict[float, dict[str, object]],
    case: Path,
    tolerance: float = 0.10,
) -> dict[str, object] | None:
    rows = []
    for time_value, row in rows_by_time.items():
        try:
            depth = float(row.get("keyholeDepth_um", np.nan))
        except (TypeError, ValueError):
            depth = np.nan
        if np.isfinite(depth) and depth > 0.0:
            rows.append((float(time_value), depth))

    if len(rows) < 3:
        print("Not enough finite keyhole-depth rows to locate stable melt-pool section zone.")
        return None

    rows.sort()
    times = np.array([row[0] for row in rows], dtype=float)
    depths = np.array([row[1] for row in rows], dtype=float)
    window = min(9, len(depths))
    if window % 2 == 0:
        window -= 1
    window = max(3, window)
    local_mean = (
        pd.Series(depths)
        .rolling(window=window, center=True, min_periods=3)
        .mean()
        .bfill()
        .ffill()
        .to_numpy()
    )
    stable = np.isfinite(local_mean) & (np.abs(depths - local_mean) <= tolerance * local_mean)

    segments = []
    start = None
    for index, is_stable in enumerate(stable):
        if is_stable and start is None:
            start = index
        elif not is_stable and start is not None:
            if index - start >= 3:
                segments.append((start, index - 1))
            start = None
    if start is not None and len(stable) - start >= 3:
        segments.append((start, len(stable) - 1))

    z_um = scan_z_from_time(case, times)
    if segments:
        i0, i1 = max(segments, key=lambda pair: abs(z_um[pair[1]] - z_um[pair[0]]))
    else:
        print("No contiguous stable keyhole segment found; using full finite keyhole-depth history.")
        i0, i1 = 0, len(times) - 1
        segments = [(i0, i1)]

    z0, z1 = sorted((float(z_um[i0]), float(z_um[i1])))
    if not np.isfinite(z0) or not np.isfinite(z1) or abs(z1 - z0) < 1e-9:
        print("Stable keyhole zone has zero scan length; cannot place melt-pool sections.")
        return None

    stable_segments = []
    for start_index, end_index in segments:
        segment_z0, segment_z1 = sorted((float(z_um[start_index]), float(z_um[end_index])))
        if np.isfinite(segment_z0) and np.isfinite(segment_z1) and abs(segment_z1 - segment_z0) > 1e-9:
            stable_segments.append(
                {
                    "time_start_s": float(times[start_index]),
                    "time_end_s": float(times[end_index]),
                    "z_start_um": segment_z0,
                    "z_end_um": segment_z1,
                    "mean_keyhole_depth_um": float(np.mean(depths[start_index : end_index + 1])),
                }
            )

    return {
        "time_start_s": float(times[i0]),
        "time_end_s": float(times[i1]),
        "z_start_um": z0,
        "z_end_um": z1,
        "z_sections_um": np.array([z0 + 0.25 * (z1 - z0), 0.5 * (z0 + z1), z0 + 0.75 * (z1 - z0)]),
        "mean_keyhole_depth_um": float(np.mean(depths[i0 : i1 + 1])),
        "stable_segments": stable_segments,
    }


def alpha_depression_depth_um(
    mesh: pv.DataSet,
    z_um: float,
    surface_y_um: float,
) -> float:
    source = mesh.cell_data_to_point_data(pass_cell_data=True)
    slc = source.slice(normal=(0.0, 0.0, 1.0), origin=(0.0, 0.0, z_um * 1e-6))
    if slc.n_points < 3 or "alpha.metal" not in slc.point_data:
        return 0.0

    xyz = slc.points
    x_um = xyz[:, 0] * 1e6
    y_um = xyz[:, 1] * 1e6
    points = np.column_stack([x_um, y_um])
    xg = np.linspace(float(x_um.min()), float(x_um.max()), NX)
    yg = np.linspace(float(y_um.min()), float(y_um.max()), NY)
    Xg, Yg = np.meshgrid(xg, yg)

    alpha = np.asarray(slc.point_data["alpha.metal"])
    Ag = griddata(points, alpha, (Xg, Yg), method="linear")
    Ag_near = griddata(points, alpha, (Xg, Yg), method="nearest")
    Ag[np.isnan(Ag)] = Ag_near[np.isnan(Ag)]

    fig_tmp, ax_tmp = plt.subplots()
    contour = ax_tmp.contour(Xg, Yg, Ag, levels=[0.5])
    segments = [seg for seg in contour.allsegs[0] if len(seg) > 0]
    plt.close(fig_tmp)
    if not segments:
        return 0.0

    pts = np.vstack(segments)
    below_surface = pts[:, 1] >= surface_y_um
    if not np.any(below_surface):
        return 0.0
    return float(np.max(pts[below_surface, 1]) - surface_y_um)


def contour_surface_width(
    segments: list[np.ndarray],
    surface_y_um: float,
) -> tuple[float, float, float] | None:
    intersections = []
    for segment in segments:
        if len(segment) < 2:
            continue
        for index in range(len(segment) - 1):
            x1, y1 = segment[index]
            x2, y2 = segment[index + 1]
            if y1 == y2:
                continue
            if (y1 - surface_y_um) * (y2 - surface_y_um) > 0:
                continue
            t = (surface_y_um - y1) / (y2 - y1)
            if 0.0 <= t <= 1.0:
                intersections.append(float(x1 + t * (x2 - x1)))

    if len(intersections) < 2:
        return None

    intersections = sorted(intersections)
    unique = []
    for value in intersections:
        if len(unique) == 0 or abs(value - unique[-1]) > 1e-6:
            unique.append(value)

    if len(unique) < 2:
        return None

    return min(unique), max(unique), max(unique) - min(unique)


def select_non_keyhole_section_z(
    stable_zone: dict[str, object],
    case: Path,
    final_mesh: pv.DataSet,
    surface_y_um: float,
    tolerance_um: float = KEYHOLE_DEPRESSION_TOLERANCE_UM,
    candidate_count: int = SECTION_CANDIDATE_COUNT,
) -> tuple[np.ndarray, dict[str, object]]:
    stable_segments = stable_zone.get("stable_segments") or [
        {
            "z_start_um": float(stable_zone["z_start_um"]),
            "z_end_um": float(stable_zone["z_end_um"]),
            "time_start_s": float(stable_zone["time_start_s"]),
            "time_end_s": float(stable_zone["time_end_s"]),
            "mean_keyhole_depth_um": float(stable_zone["mean_keyhole_depth_um"]),
        }
    ]
    source = final_mesh.cell_data_to_point_data(pass_cell_data=True)
    if "alpha.metal" not in source.point_data:
        z0 = float(stable_zone["z_start_um"])
        z1 = float(stable_zone["z_end_um"])
        return np.array([z0 + 0.25 * (z1 - z0), 0.5 * (z0 + z1), z0 + 0.75 * (z1 - z0)]), {
            "accepted_count": 0,
            "rejected_count": 0,
            "tolerance_um": tolerance_um,
            "reason": "final VTU missing alpha.metal for keyhole-footprint filter",
        }

    accepted = []
    rejected = []
    accepted_runs = []

    for segment_index, segment in enumerate(stable_segments):
        z0 = float(segment["z_start_um"])
        z1 = float(segment["z_end_um"])
        candidates = np.linspace(z0, z1, max(candidate_count, 3))
        candidate_times = scan_time_from_z(case, candidates)
        current_run = []

        for z_candidate, time_candidate in zip(candidates, candidate_times):
            depression_um = alpha_depression_depth_um(source, float(z_candidate), surface_y_um)
            item = {
                "z_um": float(z_candidate),
                "time_s": float(time_candidate),
                "vtk_time_s": np.nan,
                "depression_um": depression_um,
                "stable_segment": segment_index + 1,
            }
            if depression_um <= tolerance_um:
                accepted.append(item)
                current_run.append(item)
            else:
                rejected.append(item)
                if len(current_run) >= 3:
                    accepted_runs.append(current_run)
                current_run = []

        if len(current_run) >= 3:
            accepted_runs.append(current_run)

    if accepted_runs:
        selected_run = max(
            accepted_runs,
            key=lambda run: abs(float(run[-1]["z_um"]) - float(run[0]["z_um"])),
        )
        accepted_z = np.array([item["z_um"] for item in selected_run], dtype=float)
        targets = np.quantile(accepted_z, np.array([0.25, 0.5, 0.75]))
        selected = []
        for target in targets:
            ordered = np.argsort(np.abs(accepted_z - target))
            for index in ordered:
                z_value = float(accepted_z[index])
                if not any(np.isclose(z_value, used) for used in selected):
                    selected.append(z_value)
                    break
        selected_z = np.array(sorted(selected), dtype=float)
        selected_z0 = float(min(accepted_z))
        selected_z1 = float(max(accepted_z))
        stable_zone["z_start_um"] = selected_z0
        stable_zone["z_end_um"] = selected_z1
        stable_zone["time_start_s"] = float(scan_time_from_z(case, np.array([selected_z0]))[0])
        stable_zone["time_end_s"] = float(scan_time_from_z(case, np.array([selected_z1]))[0])
    else:
        reason = (
            "keyhole-footprint filter left fewer than 3 accepted candidate sections; "
            "not selecting fallback sections inside the rejected/depressed final alpha region"
        )
        print(reason)
        selected_z = np.array([], dtype=float)

    result = {
        "accepted_count": len(accepted),
        "rejected_count": len(rejected),
        "tolerance_um": tolerance_um,
        "accepted": accepted,
        "rejected": rejected,
        "accepted_run_count": len(accepted_runs),
    }
    if not accepted_runs:
        result["reason"] = reason

    return selected_z, result


def invalid_final_section(reason: str) -> dict[str, object]:
    return {
        "z_um": np.nan,
        "valid": False,
        "reason": reason,
    }


def final_meltpool_slice(
    mesh: pv.DataSet,
    z_um: float,
    surface_y_um: float,
    t_threshold: float,
) -> dict[str, object]:
    source = mesh.cell_data_to_point_data(pass_cell_data=True)
    slc = source.slice(normal=(0.0, 0.0, 1.0), origin=(0.0, 0.0, z_um * 1e-6))
    if slc.n_points < 3:
        return {"z_um": z_um, "valid": False, "reason": "slice has fewer than 3 points"}

    required = ("alpha.metal", "TmaxHistory")
    missing = [name for name in required if name not in slc.point_data]
    if missing:
        available = sorted(set(slc.point_data.keys()) | set(slc.cell_data.keys()))
        return {
            "z_um": z_um,
            "valid": False,
            "reason": f"missing fields {missing}; available fields: {available}",
        }

    xyz = slc.points
    x_um = xyz[:, 0] * 1e6
    y_um = xyz[:, 1] * 1e6
    points = np.column_stack([x_um, y_um])
    xg = np.linspace(float(x_um.min()), float(x_um.max()), NX)
    yg = np.linspace(float(y_um.min()), float(y_um.max()), NY)
    Xg, Yg = np.meshgrid(xg, yg)

    alpha = np.asarray(slc.point_data["alpha.metal"])
    Ag = griddata(points, alpha, (Xg, Yg), method="linear")
    Ag_near = griddata(points, alpha, (Xg, Yg), method="nearest")
    Ag[np.isnan(Ag)] = Ag_near[np.isnan(Ag)]

    history_values = np.asarray(slc.point_data["TmaxHistory"])
    Hg = griddata(points, history_values, (Xg, Yg), method="linear")
    Hg_near = griddata(points, history_values, (Xg, Yg), method="nearest")
    Hg[np.isnan(Hg)] = Hg_near[np.isnan(Hg)]
    contour_values = np.ma.masked_where(Ag < 0.5, Hg)
    contour_level = t_threshold
    field_label = f"TmaxHistory = {t_threshold:.0f} K"
    plot_values = Hg

    fig_tmp, ax_tmp = plt.subplots()
    contour = ax_tmp.contour(Xg, Yg, contour_values, levels=[contour_level])
    segments = [seg for seg in contour.allsegs[0] if len(seg) > 0]
    plt.close(fig_tmp)

    if not segments:
        return {
            "z_um": z_um,
            "valid": False,
            "reason": f"no melted-metal {field_label} contour",
            "Xg": Xg,
            "Yg": Yg,
            "plot_values": plot_values,
            "contour_values": contour_values,
            "contour_level": contour_level,
            "field_label": field_label,
        }

    pts = np.vstack(segments)
    below_surface = pts[:, 1] >= surface_y_um
    if not np.any(below_surface):
        return {
            "z_um": z_um,
            "valid": False,
            "reason": "melted-metal contour is not below the substrate surface",
            "Xg": Xg,
            "Yg": Yg,
            "plot_values": plot_values,
            "contour_values": contour_values,
            "contour_level": contour_level,
            "field_label": field_label,
        }

    pts = pts[below_surface]
    surface_width = contour_surface_width(segments, surface_y_um)
    if surface_width is None:
        return {
            "z_um": z_um,
            "valid": False,
            "reason": "melt boundary has fewer than two intersections with substrate surface",
            "Xg": Xg,
            "Yg": Yg,
            "plot_values": plot_values,
            "contour_values": contour_values,
            "contour_level": contour_level,
            "field_label": field_label,
        }

    bottom_i = int(np.argmax(pts[:, 1]))
    x_left, x_right, width_um = surface_width
    y_bottom = float(pts[bottom_i, 1])

    return {
        "z_um": z_um,
        "valid": True,
        "width_um": width_um,
        "depth_um": y_bottom - surface_y_um,
        "x_left_um": x_left,
        "x_right_um": x_right,
        "x_bottom_um": float(pts[bottom_i, 0]),
        "y_bottom_um": y_bottom,
        "Xg": Xg,
        "Yg": Yg,
        "plot_values": plot_values,
        "contour_values": contour_values,
        "contour_level": contour_level,
        "field_label": field_label,
        "segments": segments,
    }


def write_final_meltpool_summary(csv_path: Path, section_results: list[dict[str, object]]) -> None:
    fields = [
        "section",
        "z_um",
        "depth_um",
        "width_um",
        "x_left_um",
        "x_right_um",
        "x_bottom_um",
        "y_bottom_um",
    ]
    valid = [result for result in section_results if result.get("valid")]
    with csv_path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, extrasaction="ignore")
        writer.writeheader()
        for index, result in enumerate(section_results, start=1):
            row = {field_name: np.nan for field_name in fields}
            row["section"] = f"section_{index}"
            for field_name in fields[1:]:
                row[field_name] = result.get(field_name, np.nan)
            writer.writerow(row)
        if valid:
            writer.writerow(
                {
                    "section": "average",
                    "z_um": np.nan,
                    "depth_um": float(np.mean([float(result["depth_um"]) for result in valid])),
                    "width_um": float(np.mean([float(result["width_um"]) for result in valid])),
                    "x_left_um": np.nan,
                    "x_right_um": np.nan,
                    "x_bottom_um": np.nan,
                    "y_bottom_um": np.nan,
                }
            )


def write_keyhole_filter_summary(csv_path: Path, keyhole_filter: dict[str, object]) -> None:
    fields = ["status", "z_um", "time_s", "vtk_time_s", "depression_um", "tolerance_um"]
    rows = []
    for status in ("accepted", "rejected"):
        for item in keyhole_filter.get(status, []):
            row = {field_name: item.get(field_name, np.nan) for field_name in fields}
            row["status"] = status
            row["tolerance_um"] = keyhole_filter.get("tolerance_um", np.nan)
            rows.append(row)

    with csv_path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def plot_final_meltpool_sections(
    png_path: Path,
    section_results: list[dict[str, object]],
    stable_zone: dict[str, object],
    surface_y_um: float,
    case_info: dict[str, object] | None = None,
) -> None:
    case_info = case_info or {}
    valid = [result for result in section_results if result.get("valid")]
    avg_depth = float(np.mean([float(result["depth_um"]) for result in valid])) if valid else np.nan
    avg_width = float(np.mean([float(result["width_um"]) for result in valid])) if valid else np.nan
    target_depth = float(case_info["target_depth_um"]) if case_info.get("target_depth_um") is not None else np.nan
    target_width = float(case_info["target_width_um"]) if case_info.get("target_width_um") is not None else np.nan

    fig = plt.figure(figsize=(18, 12), constrained_layout=True)
    grid = fig.add_gridspec(3, 6, height_ratios=[1.05, 1.0, 0.85])
    section_axes = [fig.add_subplot(grid[0, 2 * col : 2 * col + 2]) for col in range(3)]
    mask_axes = [fig.add_subplot(grid[1, 2 * col : 2 * col + 2]) for col in range(3)]
    bar_ax = fig.add_subplot(grid[2, :3])
    compare_ax = fig.add_subplot(grid[2, 3:])

    for index, (ax, result) in enumerate(zip(section_axes, section_results), start=1):
        ax.set_title(f"Section {index}: z = {float(result['z_um']):.1f} um")
        ax.set_xlabel("X / width (um)")
        ax.set_ylabel("Y / depth (um)")
        if "Xg" not in result:
            ax.text(0.5, 0.5, result.get("reason", "invalid section"), transform=ax.transAxes, ha="center", va="center")
            ax.set_axis_off()
            continue

        Xg = result["Xg"]
        Yg = result["Yg"]
        plot_values = result["plot_values"]
        contour_values = result["contour_values"]
        contour_level = float(result["contour_level"])
        ax.pcolormesh(Xg, Yg, plot_values, shading="auto", cmap="inferno", rasterized=True)
        ax.contour(Xg, Yg, contour_values, levels=[contour_level], colors="black", linewidths=2)
        ax.axhline(surface_y_um, linestyle="--", linewidth=1.5, color="#D62728", zorder=10)
        ax.axis("equal")
        ax.invert_yaxis()

        if result.get("valid"):
            x_left = float(result["x_left_um"])
            x_right = float(result["x_right_um"])
            y_bottom = float(result["y_bottom_um"])
            x_bottom = float(result["x_bottom_um"])
            ax.plot([x_left, x_right], [surface_y_um, surface_y_um], color="#1F77B4", linewidth=2.5)
            ax.scatter([x_left, x_right, x_bottom], [surface_y_um, surface_y_um, y_bottom], color="#1F77B4", edgecolor="black", zorder=20)
            ax.annotate(
                "",
                xy=(x_bottom, y_bottom),
                xytext=(x_bottom, surface_y_um),
                arrowprops=dict(arrowstyle="<->", color="#D62728", linewidth=2.0),
            )
            text = (
                f"depth = {float(result['depth_um']):.1f} um\n"
                f"width = {float(result['width_um']):.1f} um\n"
                f"{result.get('field_label', '')}"
            )
        else:
            text = str(result.get("reason", "invalid section"))
        ax.text(
            0.02,
            0.98,
            text,
            transform=ax.transAxes,
            va="top",
            ha="left",
            bbox=dict(facecolor="white", edgecolor="black", alpha=0.85),
        )

    for index, (ax, result) in enumerate(zip(mask_axes, section_results), start=1):
        ax.set_title(f"Binary mask {index}")
        ax.set_xlabel("X / width (um)")
        ax.set_ylabel("Y / depth (um)")
        if "Xg" not in result:
            ax.text(0.5, 0.5, result.get("reason", "invalid section"), transform=ax.transAxes, ha="center", va="center")
            ax.set_axis_off()
            continue

        Xg = result["Xg"]
        Yg = result["Yg"]
        contour_values = np.ma.asarray(result["contour_values"])
        contour_level = float(result["contour_level"])
        binary_mask = np.ma.filled(contour_values >= contour_level, False).astype(float)
        ax.pcolormesh(Xg, Yg, binary_mask, shading="auto", cmap="Greys", vmin=0, vmax=1)
        ax.contour(Xg, Yg, contour_values, levels=[contour_level], colors="#0072B2", linewidths=1.8)
        ax.axhline(surface_y_um, linestyle="--", linewidth=1.5, color="#D62728", zorder=10)
        ax.axis("equal")
        ax.invert_yaxis()

        if result.get("valid"):
            x_left = float(result["x_left_um"])
            x_right = float(result["x_right_um"])
            y_bottom = float(result["y_bottom_um"])
            x_bottom = float(result["x_bottom_um"])
            ax.plot([x_left, x_right], [surface_y_um, surface_y_um], color="#E69F00", linewidth=2.5)
            ax.scatter([x_left, x_right, x_bottom], [surface_y_um, surface_y_um, y_bottom], color="#E69F00", edgecolor="black", zorder=20)

    section_labels = [f"S{index}" for index in range(1, len(section_results) + 1)]
    depth_values = [float(result["depth_um"]) if result.get("valid") else np.nan for result in section_results]
    width_values = [float(result["width_um"]) if result.get("valid") else np.nan for result in section_results]
    metric_centers = np.array([0.0, 1.0], dtype=float)
    offsets = np.linspace(-0.22, 0.22, len(section_results))
    bar_width = 0.14
    section_colors = ["#D55E00", "#0072B2", "#009E73"]
    for index, (label, offset, color) in enumerate(zip(section_labels, offsets, section_colors), start=0):
        values = [depth_values[index], width_values[index]]
        bars = bar_ax.bar(metric_centers + offset, values, bar_width, label=label, color=color)
        for bar, value in zip(bars, values):
            if np.isfinite(value):
                bar_ax.text(
                    bar.get_x() + bar.get_width() / 2,
                    value,
                    f"{value:.1f}",
                    ha="center",
                    va="bottom",
                    fontsize=12,
                )
    if np.isfinite(avg_depth):
        bar_ax.hlines(avg_depth, metric_centers[0] - 0.36, metric_centers[0] + 0.36, colors="#D55E00", linestyles="--", linewidth=2.0)
        bar_ax.text(metric_centers[0] + 0.39, avg_depth, f"avg {avg_depth:.1f}", va="center", ha="left", color="#D55E00", fontsize=12)
    if np.isfinite(avg_width):
        bar_ax.hlines(avg_width, metric_centers[1] - 0.36, metric_centers[1] + 0.36, colors="#0072B2", linestyles="--", linewidth=2.0)
        bar_ax.text(metric_centers[1] + 0.39, avg_width, f"avg {avg_width:.1f}", va="center", ha="left", color="#0072B2", fontsize=12)
    bar_ax.set_xticks(metric_centers)
    bar_ax.set_xticklabels(["Depth", "Width"])
    bar_ax.set_ylabel("Dimension (um)")
    bar_ax.set_title("Final melt-pool dimensions")
    bar_ax.grid(axis="y", alpha=0.3)
    bar_ax.legend(title="Section")
    finite_values = [value for value in depth_values + width_values + [avg_depth, avg_width] if np.isfinite(value)]
    if finite_values:
        bar_ax.set_ylim(0.0, max(finite_values) * 1.18)

    compare_metrics = ["Depth", "Width"]
    sim_values = [avg_depth, avg_width]
    exp_values = [target_depth, target_width]
    compare_centers = np.arange(len(compare_metrics), dtype=float)
    compare_width = 0.32
    sim_bars = compare_ax.bar(
        compare_centers - compare_width / 2,
        sim_values,
        compare_width,
        label="Sim avg",
        color="#0072B2",
    )
    exp_bars = compare_ax.bar(
        compare_centers + compare_width / 2,
        exp_values,
        compare_width,
        label="Experiment",
        color="#CC79A7",
    )
    compare_ax.set_xticks(compare_centers)
    compare_ax.set_xticklabels(compare_metrics)
    compare_ax.set_ylabel("Dimension (um)")
    compare_ax.set_title("Average simulation vs experiment")
    compare_ax.grid(axis="y", alpha=0.3)
    compare_ax.legend()
    compare_finite = [value for value in sim_values + exp_values if np.isfinite(value)]
    if compare_finite:
        compare_ax.set_ylim(0.0, max(compare_finite) * 1.18)
    else:
        compare_ax.text(
            0.5,
            0.5,
            "No target_width_um / target_depth_um in case_info.json",
            transform=compare_ax.transAxes,
            ha="center",
            va="center",
        )
    for bars in (sim_bars, exp_bars):
        for bar in bars:
            value = float(bar.get_height())
            if np.isfinite(value):
                compare_ax.text(
                    bar.get_x() + bar.get_width() / 2,
                    value,
                    f"{value:.1f}",
                    ha="center",
                    va="bottom",
                    fontsize=12,
                )

    fig.suptitle(
        "Final melt-pool sections from peak-temperature history when available\n"
        f"stable z = {float(stable_zone['z_start_um']):.1f}-{float(stable_zone['z_end_um']):.1f} um; "
        f"average depth = {avg_depth:.1f} um, average width = {avg_width:.1f} um"
    )
    png_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(png_path, format="png", bbox_inches="tight", dpi=PNG_DPI)
    plt.close(fig)


def plot_final_melttrack_yz(
    png_path: Path,
    mesh: pv.DataSet,
    stable_zone: dict[str, object],
    section_z_um: np.ndarray,
    surface_y_um: float,
    t_threshold: float,
) -> None:
    source = mesh.cell_data_to_point_data(pass_cell_data=True)
    x_mid = 0.5 * (mesh.bounds[0] + mesh.bounds[1])
    slc = source.slice(normal=(1.0, 0.0, 0.0), origin=(x_mid, 0.0, 0.0))

    fig, ax = plt.subplots(figsize=(12, 5.5), constrained_layout=True)
    ax.set_title(f"Final melted-metal track, YZ slice at x = {x_mid * 1e6:.1f} um")
    ax.set_xlabel("Z / scan track (um)")
    ax.set_ylabel("Y / depth (um)")

    if slc.n_points < 3:
        ax.text(0.5, 0.5, "YZ slice has fewer than 3 points", transform=ax.transAxes, ha="center", va="center")
        ax.set_axis_off()
    else:
        required = ("alpha.metal", "TmaxHistory")
        missing = [name for name in required if name not in slc.point_data]
        if missing:
            ax.text(
                0.5,
                0.5,
                f"Missing fields: {missing}",
                transform=ax.transAxes,
                ha="center",
                va="center",
            )
            ax.set_axis_off()
        else:
            xyz = slc.points
            z_um = xyz[:, 2] * 1e6
            y_um = xyz[:, 1] * 1e6
            points = np.column_stack([z_um, y_um])
            zg = np.linspace(float(z_um.min()), float(z_um.max()), NZ)
            yg = np.linspace(float(y_um.min()), float(y_um.max()), NY)
            Zg, Yg = np.meshgrid(zg, yg)

            alpha = np.asarray(slc.point_data["alpha.metal"])
            Ag = griddata(points, alpha, (Zg, Yg), method="linear")
            Ag_near = griddata(points, alpha, (Zg, Yg), method="nearest")
            Ag[np.isnan(Ag)] = Ag_near[np.isnan(Ag)]

            history_values = np.asarray(slc.point_data["TmaxHistory"])
            Hg = griddata(points, history_values, (Zg, Yg), method="linear")
            Hg_near = griddata(points, history_values, (Zg, Yg), method="nearest")
            Hg[np.isnan(Hg)] = Hg_near[np.isnan(Hg)]
            contour_values = np.ma.masked_where(Ag < 0.5, Hg)
            contour_level = t_threshold
            ax.pcolormesh(Zg, Yg, Hg, shading="auto", cmap="inferno", rasterized=True)

            ax.contour(Zg, Yg, contour_values, levels=[contour_level], colors="black", linewidths=2)
            ax.contour(
                Zg,
                Yg,
                Ag,
                levels=[0.5],
                colors="#00BFC4",
                linewidths=1.8,
                linestyles="--",
            )
            ax.axhline(surface_y_um, linestyle="--", linewidth=1.5, color="#D62728", zorder=10)
            ax.axvspan(
                float(stable_zone["z_start_um"]),
                float(stable_zone["z_end_um"]),
                color="#1F77B4",
                alpha=0.12,
                label="stable keyhole zone",
            )
            for index, z_section in enumerate(section_z_um, start=1):
                ax.axvline(float(z_section), color="#1F77B4", linestyle=":", linewidth=1.6)
                ax.text(
                    float(z_section),
                    0.03,
                    f"S{index}",
                    transform=ax.get_xaxis_transform(),
                    ha="center",
                    va="bottom",
                    color="#1F77B4",
                    bbox=dict(facecolor="white", edgecolor="#1F77B4", alpha=0.85),
                )
            ax.axis("equal")
            ax.invert_yaxis()
            handles, labels = ax.get_legend_handles_labels()
            handles.extend(
                [
                    Line2D([0], [0], color="black", linewidth=2, label="melt boundary"),
                    Line2D([0], [0], color="#00BFC4", linestyle="--", linewidth=1.8, label="final alpha.metal = 0.5"),
                ]
            )
            labels.extend(["melt boundary", "final alpha.metal = 0.5"])
            if handles:
                ax.legend(handles, labels, frameon=True, facecolor="white", edgecolor="black")

    png_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(png_path, format="png", bbox_inches="tight", dpi=PNG_DPI)
    plt.close(fig)


def analyze_final_meltpool_sections(
    final_vtk: Path,
    rows_by_time: dict[float, dict[str, object]],
    case: Path,
    surface_y_um: float,
    t_threshold: float,
    out_dir: Path,
    case_info: dict[str, object] | None = None,
) -> None:
    stable_zone = stable_keyhole_sections(rows_by_time, case)
    if stable_zone is None:
        return

    mesh = pv.read(final_vtk)
    z_min_um = mesh.bounds[4] * 1e6
    z_max_um = mesh.bounds[5] * 1e6
    section_z_um, keyhole_filter = select_non_keyhole_section_z(
        stable_zone,
        case,
        mesh,
        surface_y_um,
    )
    stable_zone["z_sections_um"] = section_z_um
    stable_zone["keyhole_filter"] = keyhole_filter
    print(
        "Keyhole-footprint filter accepted "
        f"{keyhole_filter['accepted_count']} candidates and rejected "
        f"{keyhole_filter['rejected_count']} candidates "
        f"(depression tolerance = {keyhole_filter['tolerance_um']:.1f} um)."
    )
    section_z_um = np.clip(np.asarray(section_z_um, dtype=float), z_min_um, z_max_um)
    if len(section_z_um) < 3:
        reason = keyhole_filter.get("reason", "fewer than 3 valid final melt-pool sections")
        section_results = [invalid_final_section(str(reason)) for _ in range(3)]
    else:
        section_results = [
            final_meltpool_slice(mesh, float(z_um), surface_y_um, t_threshold)
            for z_um in section_z_um
        ]

    csv_path = out_dir / "final_meltpool_dimensions.csv"
    filter_csv_path = out_dir / "final_meltpool_keyhole_filter.csv"
    png_path = out_dir / "final_meltpool_sections.png"
    yz_png_path = out_dir / "final_melttrack_yz.png"
    write_final_meltpool_summary(csv_path, section_results)
    write_keyhole_filter_summary(filter_csv_path, keyhole_filter)
    plot_final_meltpool_sections(png_path, section_results, stable_zone, surface_y_um, case_info)
    plot_final_melttrack_yz(yz_png_path, mesh, stable_zone, section_z_um, surface_y_um, t_threshold)
    print(f"Wrote {csv_path}")
    print(f"Wrote {filter_csv_path}")
    print(f"Wrote {png_path}")
    print(f"Wrote {yz_png_path}")


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
    png_dir = out_dir / "vtk_sections"
    out_dir.mkdir(exist_ok=True)
    png_dir.mkdir(exist_ok=True)

    t_threshold = args.t_threshold if args.t_threshold is not None else args.t_liquidus
    if t_threshold is None:
        t_threshold = parse_scalar(case / "constant" / "transportProperties", "Tsolidus")

    if args.vtk_file:
        vtk_files = [args.vtk_file.resolve()]
    else:
        vtk_files = internal_mesh_vtk_files(vtk_dir)

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
                keep_cols = ["t_ms", "keyhole_depth_um"]
                if "keyhole_depth_std_um" in _df.columns:
                    keep_cols.append("keyhole_depth_std_um")
                if "absorptance" in _df.columns:
                    keep_cols.append("absorptance")
                exp_timeseries_df = _df[keep_cols].sort_values("t_ms")
            elif {"time_ms", "depth_mean"}.issubset(_df.columns):
                exp_timeseries_df = pd.DataFrame(
                    {
                        "t_ms": _df["time_ms"],
                        "keyhole_depth_um": _df["depth_mean"],
                    }
                )
                if "depth_std" in _df.columns:
                    exp_timeseries_df["keyhole_depth_std_um"] = _df["depth_std"]
                if "absorptance" in _df.columns:
                    exp_timeseries_df["absorptance"] = _df["absorptance"]
                exp_timeseries_df = exp_timeseries_df.sort_values("t_ms")
            elif {"time_ms", "depth_um"}.issubset(_df.columns):
                exp_timeseries_df = pd.DataFrame(
                    {
                        "t_ms": _df["time_ms"],
                        "keyhole_depth_um": _df["depth_um"],
                    }
                )
                if "absorptance" in _df.columns:
                    exp_timeseries_df["absorptance"] = _df["absorptance"]
                exp_timeseries_df = exp_timeseries_df.sort_values("t_ms")
            else:
                print(
                    "exp-timeseries-csv missing required columns "
                    f"(t_ms, keyhole_depth_um) or (time_ms, depth_mean): {exp_ts_path}"
                )
        else:
            print(f"exp-timeseries-csv not found: {exp_ts_path}")

    csv_path = out_dir / "vtk_meltpool_geometry.csv"
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

    if vtk_files:
        final_vtk = max(vtk_files, key=vtk_time)
        analyze_final_meltpool_sections(
            final_vtk,
            rows_by_time,
            case,
            args.surface_y_um,
            t_threshold,
            out_dir,
            case_info,
        )

    print(f"Wrote {csv_path}")
    print(f"Wrote PNG sections to {png_dir}")


if __name__ == "__main__":
    main()
