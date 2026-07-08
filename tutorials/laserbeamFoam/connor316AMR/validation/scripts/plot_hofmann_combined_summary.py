#!/usr/bin/env python3
"""Plot combined Hofmann final melt-pool section summary."""

from __future__ import annotations

import csv
import json
import re
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pyvista as pv
from matplotlib.lines import Line2D
from scipy.interpolate import griddata


ROOT = Path(__file__).resolve().parents[1]
CASES_DIR = ROOT / "runs"
OUT_PNG = ROOT / "post-processing-data" / "hofmann_combined_summary.png"
OUT_PDF = ROOT / "post-processing-data" / "hofmann_combined_summary.pdf"
HOFMANN_CASES = [
    "hofmann_scantrack_200W_900mms_r25um",
    "hofmann_scantrack_200W_1200mms_r25um",
    "hofmann_scantrack_200W_1500mms_r25um",
    "hofmann_scantrack_250W_600mms_r25um",
]
NX = 180
NY = 180
SURFACE_Y_UM = 128.0
PNG_DPI = 150
TITLE_FONTSIZE = 14
LABEL_FONTSIZE = 12
TICK_FONTSIZE = 10
ANNOTATION_FONTSIZE = 10


def condition_label(case_name: str) -> str:
    match = re.search(r"_(\d+)W_(\d+)mms_", case_name)
    if not match:
        return case_name
    return f"{match.group(1)} W, {match.group(2)} mm/s"


def parse_scalar(path: Path, name: str) -> float:
    pattern = re.compile(rf"^\s*{re.escape(name)}\s+([-+0-9.eE]+)\s*;")
    for line in path.read_text(encoding="utf-8").splitlines():
        match = pattern.match(line)
        if match:
            return float(match.group(1))
    raise KeyError(f"{name} not found in {path}")


def vtk_time(path: Path) -> float:
    match = re.search(r"_([-+0-9.eE]+)\.vtk$", path.name)
    return float(match.group(1)) if match else -1.0


def latest_vtk(case_dir: Path) -> Path:
    files = sorted((case_dir / "VTK").glob("*.vtk"), key=vtk_time)
    if not files:
        raise FileNotFoundError(f"No legacy internal-mesh .vtk found in {case_dir / 'VTK'}")
    return files[-1]


def read_dimensions(case_dir: Path) -> dict[str, dict[str, float | str]]:
    csv_path = case_dir / "post-processing-data" / "final_meltpool_dimensions.csv"
    if not csv_path.exists():
        raise FileNotFoundError(csv_path)
    rows = {}
    with csv_path.open(newline="") as handle:
        for row in csv.DictReader(handle):
            rows[str(row["section"])] = {
                key: (float(value) if key != "section" and value not in ("", "nan") else value)
                for key, value in row.items()
            }
    return rows


def load_case_summary(case_name: str) -> dict[str, object]:
    case_dir = CASES_DIR / case_name
    info = json.loads((case_dir / "case_info.json").read_text(encoding="utf-8"))
    dims = read_dimensions(case_dir)
    section_1 = dims["section_1"]
    average = dims["average"]
    t_solidus = parse_scalar(case_dir / "constant" / "transportProperties", "Tsolidus")
    return {
        "name": case_name,
        "label": condition_label(case_name),
        "dir": case_dir,
        "info": info,
        "vtk": latest_vtk(case_dir),
        "section_1": section_1,
        "average": average,
        "t_solidus": t_solidus,
        "target_depth_um": float(info["target_depth_um"]),
        "target_width_um": float(info["target_width_um"]),
    }


def slice_section(case: dict[str, object]) -> dict[str, object]:
    mesh = pv.read(case["vtk"])
    source = mesh.cell_data_to_point_data(pass_cell_data=True)
    z_um = float(case["section_1"]["z_um"])
    slc = source.slice(normal=(0.0, 0.0, 1.0), origin=(0.0, 0.0, z_um * 1e-6))
    if slc.n_points < 3:
        raise RuntimeError(f"{case['name']}: section 1 slice has fewer than 3 points")

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

    if "TmaxHistory" not in slc.point_data:
        raise RuntimeError(f"{case['name']}: missing TmaxHistory")
    field = np.asarray(slc.point_data["TmaxHistory"])
    field_label = f"TmaxHistory >= {float(case['t_solidus']):.0f} K"
    level = float(case["t_solidus"])

    Fg = griddata(points, field, (Xg, Yg), method="linear")
    Fg_near = griddata(points, field, (Xg, Yg), method="nearest")
    Fg[np.isnan(Fg)] = Fg_near[np.isnan(Fg)]

    contour_values = np.ma.masked_where(Ag < 0.5, Fg)
    plot_values = Fg

    return {
        "Xg": Xg,
        "Yg": Yg,
        "plot_values": plot_values,
        "contour_values": contour_values,
        "level": level,
        "field_label": field_label,
    }


def add_dimension_markers(ax, section: dict[str, float | str], color: str) -> None:
    try:
        x_left = float(section["x_left_um"])
        x_right = float(section["x_right_um"])
        x_bottom = float(section["x_bottom_um"])
        y_bottom = float(section["y_bottom_um"])
    except (KeyError, TypeError, ValueError):
        return
    ax.plot([x_left, x_right], [SURFACE_Y_UM, SURFACE_Y_UM], color=color, linewidth=2.2)
    ax.scatter(
        [x_left, x_right, x_bottom],
        [SURFACE_Y_UM, SURFACE_Y_UM, y_bottom],
        color=color,
        edgecolor="black",
        s=22,
        zorder=20,
    )


def pct_error(sim: float, exp: float) -> float:
    return 100.0 * (sim - exp) / exp if np.isfinite(sim) and np.isfinite(exp) and exp != 0 else np.nan


def style_section_axis(ax, show_ylabel: bool) -> None:
    ax.axhline(SURFACE_Y_UM, linestyle="--", linewidth=1.2, color="#D62728")
    ax.set_xlabel("Width, x (um)", fontsize=LABEL_FONTSIZE)
    if show_ylabel:
        ax.set_ylabel("Depth, y (um)", fontsize=LABEL_FONTSIZE)
    else:
        ax.set_ylabel("")
        ax.set_yticklabels([])
    ax.tick_params(labelsize=TICK_FONTSIZE)
    ax.axis("equal")
    ax.invert_yaxis()


def add_row_label(fig, axes: list[plt.Axes], label: str) -> None:
    box = axes[0].get_position()
    fig.text(
        0.025,
        0.5 * (box.y0 + box.y1),
        label,
        rotation=90,
        va="center",
        ha="center",
        fontsize=TITLE_FONTSIZE,
        fontweight="bold",
    )


def plot_summary(cases: list[dict[str, object]], output_png: Path, output_pdf: Path) -> None:
    fig = plt.figure(figsize=(22, 13.0))
    grid = fig.add_gridspec(
        3,
        4,
        height_ratios=[1.0, 1.0, 0.95],
        left=0.06,
        right=0.99,
        bottom=0.06,
        top=0.84,
        hspace=0.26,
        wspace=0.06,
    )
    temp_axes = [fig.add_subplot(grid[0, i]) for i in range(4)]
    mask_axes = [fig.add_subplot(grid[1, i]) for i in range(4)]
    depth_ax = fig.add_subplot(grid[2, :2])
    width_ax = fig.add_subplot(grid[2, 2:])

    for col, (ax, case) in enumerate(zip(temp_axes, cases)):
        section = case["section_1"]
        data = case["section_data"]
        Xg = data["Xg"]
        Yg = data["Yg"]
        ax.set_title(f"{case['label']}\nSection 1 at z = {float(section['z_um']):.0f} um", fontsize=TITLE_FONTSIZE)
        ax.pcolormesh(Xg, Yg, data["plot_values"], shading="auto", cmap="inferno", rasterized=True)
        ax.contour(Xg, Yg, data["contour_values"], levels=[float(data["level"])], colors="black", linewidths=1.8)
        add_dimension_markers(ax, section, "#0072B2")
        style_section_axis(ax, col == 0)
        ax.text(
            0.02,
            0.98,
            f"D = {float(section['depth_um']):.1f} um\nW = {float(section['width_um']):.1f} um",
            transform=ax.transAxes,
            va="top",
            ha="left",
            fontsize=ANNOTATION_FONTSIZE,
            bbox=dict(facecolor="white", edgecolor="0.3", alpha=0.9),
        )

    for col, (ax, case) in enumerate(zip(mask_axes, cases)):
        section = case["section_1"]
        data = case["section_data"]
        Xg = data["Xg"]
        Yg = data["Yg"]
        contour_values = np.ma.asarray(data["contour_values"])
        level = float(data["level"])
        binary_mask = np.ma.filled(contour_values >= level, False).astype(float)
        ax.set_title(str(case["label"]), fontsize=TITLE_FONTSIZE)
        ax.pcolormesh(Xg, Yg, binary_mask, shading="auto", cmap="Greys", vmin=0, vmax=1)
        ax.contour(Xg, Yg, contour_values, levels=[level], colors="#0072B2", linewidths=1.8)
        add_dimension_markers(ax, section, "#E69F00")
        style_section_axis(ax, col == 0)

    labels = [str(case["label"]) for case in cases]
    x = np.arange(len(labels), dtype=float)
    depth_sim = np.array([float(case["average"]["depth_um"]) for case in cases])
    depth_exp = np.array([float(case["target_depth_um"]) for case in cases])
    width_sim = np.array([float(case["average"]["width_um"]) for case in cases])
    width_exp = np.array([float(case["target_width_um"]) for case in cases])

    def plot_metric(ax, sim_values, exp_values, title, color, legend_loc="upper left"):
        bar_width = 0.34
        sim_bars = ax.bar(x - bar_width / 2, sim_values, bar_width, label="Simulation", color=color)
        exp_bars = ax.bar(x + bar_width / 2, exp_values, bar_width, label="Experiment", color="#8A8A8A")
        ax.set_title(title, fontsize=TITLE_FONTSIZE, fontweight="bold")
        ax.set_xticks(x)
        ax.set_xticklabels(labels, fontsize=TICK_FONTSIZE)
        ax.set_ylabel("Dimension (um)", fontsize=LABEL_FONTSIZE)
        ax.tick_params(axis="y", labelsize=TICK_FONTSIZE)
        ax.grid(axis="y", alpha=0.25)
        finite = np.concatenate([sim_values, exp_values])
        finite = finite[np.isfinite(finite)]
        if finite.size:
            ax.set_ylim(0.0, float(finite.max()) * 1.28)
        for bars in (sim_bars, exp_bars):
            for bar in bars:
                value = float(bar.get_height())
                if np.isfinite(value):
                    ax.text(
                        bar.get_x() + bar.get_width() / 2,
                        value,
                        f"{value:.0f}",
                        ha="center",
                        va="bottom",
                        fontsize=ANNOTATION_FONTSIZE,
                    )
        for index in range(len(labels)):
            err = pct_error(sim_values[index], exp_values[index])
            y_top = max(sim_values[index], exp_values[index])
            ax.text(
                x[index],
                y_top * 1.12,
                f"{err:+.1f}%",
                ha="center",
                va="bottom",
                fontsize=ANNOTATION_FONTSIZE,
                fontweight="bold",
                bbox=dict(facecolor="white", edgecolor="0.75", alpha=0.95),
            )
        ax.legend(loc=legend_loc, frameon=True, facecolor="white", edgecolor="0.8", fontsize=TICK_FONTSIZE)

    plot_metric(depth_ax, depth_sim, depth_exp, "Average melt-pool depth", "#D55E00")
    plot_metric(width_ax, width_sim, width_exp, "Average melt-pool width", "#0072B2", legend_loc="lower left")

    legend_handles = [
        Line2D([0], [0], color="black", linewidth=1.8, label="Melt boundary"),
        Line2D([0], [0], color="#D62728", linestyle="--", linewidth=1.2, label="Initial substrate height"),
        Line2D([0], [0], color="#0072B2", linewidth=2.2, marker="o", label="Measured dimension"),
    ]
    fig.legend(
        handles=legend_handles,
        loc="upper center",
        bbox_to_anchor=(0.5, 0.91),
        ncol=3,
        frameon=True,
        facecolor="white",
        edgecolor="0.8",
        fontsize=LABEL_FONTSIZE,
    )
    fig.suptitle(
        "Hofmann 316L Single-Track Melt-Pool Validation",
        fontsize=20,
        fontweight="bold",
        y=0.975,
    )
    fig.text(
        0.5,
        0.935,
        "Top: peak-temperature melt boundary at section 1. Middle: binary melted-region mask. Bottom: simulation average vs experimental target; labels show percent error.",
        ha="center",
        va="center",
        fontsize=LABEL_FONTSIZE,
    )
    add_row_label(fig, temp_axes, "Temperature field")
    add_row_label(fig, mask_axes, "Binary mask")
    output_png.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_png, format="png", bbox_inches="tight", dpi=PNG_DPI)
    fig.savefig(output_pdf, format="pdf", bbox_inches="tight")
    plt.close(fig)


def main() -> None:
    cases = [load_case_summary(name) for name in HOFMANN_CASES]
    for case in cases:
        case["section_data"] = slice_section(case)
    plot_summary(cases, OUT_PNG, OUT_PDF)
    print(f"Wrote {OUT_PNG}")
    print(f"Wrote {OUT_PDF}")


if __name__ == "__main__":
    main()
