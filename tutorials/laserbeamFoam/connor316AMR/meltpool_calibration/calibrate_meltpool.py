#!/usr/bin/env python3
from __future__ import annotations

import argparse
import concurrent.futures
import csv
import json
import math
import re
import shutil
import subprocess
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pandas as pd


WORK_DIR = Path(__file__).resolve().parent
BASE_CASE = WORK_DIR.parent / "template"
RUNS_DIR = WORK_DIR / "runs"
SUMMARY_CSV = WORK_DIR / "summary.csv"
BEST_TXT = WORK_DIR / "best_case.txt"
PLOTS_DIR = WORK_DIR / "plots"
DEFAULT_CONFIG = WORK_DIR / "calibration_config.json"

DEFAULT_RUN_COMMAND = "source /usr/lib/openfoam/openfoam2506/etc/bashrc; ./Allrun"


@dataclass(frozen=True)
class CaseSpec:
    radius_m: float
    resistivity_ohm_m: float

    @property
    def case_id(self) -> str:
        radius_um = self.radius_m * 1e6
        resistivity_micro = self.resistivity_ohm_m * 1e6
        return f"r{radius_um:06.1f}um_er{resistivity_micro:05.3f}e-6"


def replace_scalar(text: str, name: str, value: float) -> str:
    pattern = re.compile(rf"(^\s*{re.escape(name)}\s+)([^;]+)(;.*$)", re.MULTILINE)
    replacement = rf"\g<1>{value:.12g}\g<3>"
    text, count = pattern.subn(replacement, text, count=1)
    if count != 1:
        raise ValueError(f"Could not patch {name}")
    return text


def effective_config(config: dict[str, object], force_coarse: bool = False) -> dict[str, object]:
    config = json.loads(json.dumps(config))
    coarse = dict(config.get("coarse_test_mode", {}))
    if force_coarse:
        coarse["enabled"] = True
        config["coarse_test_mode"] = coarse
    if not bool(coarse.get("enabled", False)):
        return config

    if "write_interval_s" in coarse:
        config["write_interval_s"] = coarse["write_interval_s"]
    if "simulation_end_s" in coarse:
        config["simulation_end_s"] = coarse["simulation_end_s"]
    if "first_pulse_end_s" in coarse:
        config["first_pulse_end_s"] = coarse["first_pulse_end_s"]
    if "initial_laser_radius_um" in coarse:
        config["laser_radius_um"] = coarse["initial_laser_radius_um"]
    if "initial_elec_resistivity_ohm_m" in coarse:
        config["elec_resistivity_ohm_m"] = coarse["initial_elec_resistivity_ohm_m"]
    if "max_parallel_cases" in coarse:
        config["max_parallel_cases"] = coarse["max_parallel_cases"]
    if "z_base_cell_size_um" in coarse:
        shrink = dict(config.get("shrink_z_geometry", {}))
        shrink["base_cell_size_um"] = coarse["z_base_cell_size_um"]
        config["shrink_z_geometry"] = shrink

    adaptive = dict(config.get("adaptive", {}))
    if "adaptive_max_new_cases" in coarse:
        adaptive["max_new_cases"] = coarse["adaptive_max_new_cases"]
    if "adaptive_laser_radius_search_step_um" in coarse:
        adaptive["laser_radius_search_step_um"] = coarse["adaptive_laser_radius_search_step_um"]
    if "adaptive_elec_resistivity_search_step_ohm_m" in coarse:
        adaptive["elec_resistivity_search_step_ohm_m"] = coarse["adaptive_elec_resistivity_search_step_ohm_m"]
    if "adaptive_candidates_per_iteration" in coarse:
        adaptive["candidates_per_iteration"] = coarse["adaptive_candidates_per_iteration"]
    config["adaptive"] = adaptive
    return config


def vector_literal(values: tuple[float, float, float]) -> str:
    return f"({values[0]:.12g} {values[1]:.12g} {values[2]:.12g})"


def parse_time_vs_laser_position(path: Path) -> list[tuple[float, tuple[float, float, float]]]:
    text = path.read_text()
    rows = []
    pattern = re.compile(
        r"\(\s*([0-9.eE+-]+)\s+\(\s*([0-9.eE+-]+)\s+([0-9.eE+-]+)\s+([0-9.eE+-]+)\s*\)\s*\)"
    )
    for match in pattern.finditer(text):
        rows.append(
            (
                float(match.group(1)),
                (float(match.group(2)), float(match.group(3)), float(match.group(4))),
            )
        )
    if len(rows) < 2:
        raise ValueError(f"Could not parse laser path from {path}")
    return sorted(rows, key=lambda item: item[0])


def interpolate_laser_position(
    rows: list[tuple[float, tuple[float, float, float]]],
    time_s: float,
) -> tuple[float, float, float]:
    if time_s <= rows[0][0]:
        return rows[0][1]
    if time_s >= rows[-1][0]:
        return rows[-1][1]

    for (t0, p0), (t1, p1) in zip(rows[:-1], rows[1:]):
        if t0 <= time_s <= t1:
            f = (time_s - t0) / (t1 - t0)
            return tuple(p0[i] + f * (p1[i] - p0[i]) for i in range(3))
    return rows[-1][1]


def patch_block_mesh_dimensions(case_dir: Path, nx: int | None = None, ny: int | None = None, nz: int | None = None) -> None:
    block_mesh = case_dir / "system" / "blockMeshDict"
    text = block_mesh.read_text()
    pattern = re.compile(
        r"(hex\s*\(0\s+1\s+2\s+3\s+4\s+5\s+6\s+7\)\s*\(\s*)(\d+)(\s+)(\d+)(\s+)(\d+)(\s*\))"
    )
    match = pattern.search(text)
    if not match:
        raise ValueError(f"Could not find block cell counts in {block_mesh}")

    old_nx, old_ny, old_nz = int(match.group(2)), int(match.group(4)), int(match.group(6))
    new_nx = old_nx if nx is None else nx
    new_ny = old_ny if ny is None else ny
    new_nz = old_nz if nz is None else nz
    replacement = rf"\g<1>{new_nx}\g<3>{new_ny}\g<5>{new_nz}\g<7>"
    text = pattern.sub(replacement, text, count=1)
    block_mesh.write_text(text)


def patch_block_mesh_z(case_dir: Path, z_min: float, z_max: float, nz: int) -> None:
    block_mesh = case_dir / "system" / "blockMeshDict"
    text = block_mesh.read_text()

    vertices_match = re.search(r"(vertices\s*\(\s*)(.*?)(\s*\);)", text, flags=re.DOTALL)
    if not vertices_match:
        raise ValueError(f"Could not find vertices block in {block_mesh}")

    tuple_pattern = re.compile(r"\(\s*([0-9.eE+-]+)\s+([0-9.eE+-]+)\s+([0-9.eE+-]+)\s*\)")
    vertices = []
    for match in tuple_pattern.finditer(vertices_match.group(2)):
        vertices.append((float(match.group(1)), float(match.group(2)), float(match.group(3))))
    if len(vertices) != 8:
        raise ValueError(f"Expected 8 block vertices in {block_mesh}, found {len(vertices)}")

    new_vertices = []
    for index, (x, y, _z) in enumerate(vertices):
        new_z = z_min if index < 4 else z_max
        new_vertices.append(f"    {vector_literal((x, y, new_z))}")
    new_vertex_block = vertices_match.group(1) + "\n" + "\n".join(new_vertices) + vertices_match.group(3)
    text = text[: vertices_match.start()] + new_vertex_block + text[vertices_match.end() :]

    block_mesh.write_text(text)
    patch_block_mesh_dimensions(case_dir, nz=nz)


def patch_set_fields_z(case_dir: Path, z_min: float, z_max: float) -> None:
    set_fields = case_dir / "system" / "setFieldsDict"
    text = set_fields.read_text()
    text, count = re.subn(
        r"(box\s*\(\s*0\s+0\.3e-3\s+)([0-9.eE+-]+)(\s*\)\s*\(\s*1\.5e-3\s+1\.0e-3\s+)([0-9.eE+-]+)(\s*\)\s*;)",
        rf"\g<1>{z_min:.12g}\g<3>{z_max:.12g}\g<5>",
        text,
        count=1,
    )
    if count != 1:
        raise ValueError(f"Could not patch setFields box in {set_fields}")
    set_fields.write_text(text)


def patch_shrunk_z_geometry(case_dir: Path, config: dict[str, object]) -> None:
    shrink = dict(config.get("shrink_z_geometry", {}))
    if not bool(shrink.get("enabled", False)):
        return

    laser_path = parse_time_vs_laser_position(case_dir / "constant" / "timeVsLaserPosition")
    start_pos = interpolate_laser_position(laser_path, 0.0)
    end_pos = interpolate_laser_position(laser_path, float(config["simulation_end_s"]))
    z_path_min = min(start_pos[2], end_pos[2])
    z_path_max = max(start_pos[2], end_pos[2])

    offset_before = float(shrink.get("offset_before_um", 300.0)) * 1e-6
    offset_after = float(shrink.get("offset_after_um", 300.0)) * 1e-6
    min_length = float(shrink.get("minimum_length_um", 500.0)) * 1e-6
    cell_size = float(shrink.get("base_cell_size_um", 50.0)) * 1e-6

    z_min = max(0.0, z_path_min - offset_before)
    z_max = z_path_max + offset_after
    if z_max - z_min < min_length:
        pad = 0.5 * (min_length - (z_max - z_min))
        z_min = max(0.0, z_min - pad)
        z_max = z_max + pad

    nz = max(1, int(math.ceil((z_max - z_min) / cell_size)))
    patch_block_mesh_z(case_dir, z_min, z_max, nz)
    patch_set_fields_z(case_dir, z_min, z_max)


def patch_coarse_test_mesh(case_dir: Path, config: dict[str, object]) -> None:
    coarse = dict(config.get("coarse_test_mode", {}))
    if not bool(coarse.get("enabled", False)):
        return

    mesh_cells = coarse.get("mesh_cells", [None, None, None])
    if not isinstance(mesh_cells, list) or len(mesh_cells) != 3:
        raise ValueError("coarse_test_mode.mesh_cells must be [nx, ny, nz_or_null]")

    nx = None if mesh_cells[0] is None else int(mesh_cells[0])
    ny = None if mesh_cells[1] is None else int(mesh_cells[1])
    nz = None if mesh_cells[2] is None else int(mesh_cells[2])
    patch_block_mesh_dimensions(case_dir, nx=nx, ny=ny, nz=nz)


def patch_coarse_test_solver_cost(case_dir: Path, config: dict[str, object]) -> None:
    coarse = dict(config.get("coarse_test_mode", {}))
    if not bool(coarse.get("enabled", False)):
        return

    laser_props = case_dir / "constant" / "LaserProperties"
    text = laser_props.read_text()
    if "nRadial" in coarse:
        text = replace_scalar(text, "nRadial", int(coarse["nRadial"]))
    if "nAngular" in coarse:
        text = replace_scalar(text, "nAngular", int(coarse["nAngular"]))
    laser_props.write_text(text)

    if "numberOfSubdomains" in coarse:
        decompose = case_dir / "system" / "decomposeParDict"
        text = decompose.read_text()
        text = replace_scalar(text, "numberOfSubdomains", int(coarse["numberOfSubdomains"]))
        decompose.write_text(text)


def load_config(path: Path) -> dict[str, object]:
    with path.open() as handle:
        config = json.load(handle)

    required = [
        "target_width_um",
        "target_depth_um",
        "first_pulse_end_s",
        "simulation_end_s",
        "write_interval_s",
        "laser_radius_um",
        "elec_resistivity_ohm_m",
    ]
    missing = [key for key in required if key not in config]
    if missing:
        raise KeyError(f"Missing config keys: {', '.join(missing)}")
    return config


def patch_case(case_dir: Path, spec: CaseSpec, config: dict[str, object]) -> None:
    laser_props = case_dir / "constant" / "LaserProperties"
    transport_props = case_dir / "constant" / "transportProperties"
    control_dict = case_dir / "system" / "controlDict"

    text = laser_props.read_text()
    text = replace_scalar(text, "laserRadius", spec.radius_m)
    laser_props.write_text(text)

    text = transport_props.read_text()
    text = replace_scalar(text, "elec_resistivity", spec.resistivity_ohm_m)
    transport_props.write_text(text)

    text = control_dict.read_text()
    text = replace_scalar(text, "endTime", float(config["simulation_end_s"]))
    text = replace_scalar(text, "writeInterval", float(config["write_interval_s"]))
    text = replace_scalar(text, "purgeWrite", 0)
    control_dict.write_text(text)
    patch_shrunk_z_geometry(case_dir, config)
    patch_coarse_test_mesh(case_dir, config)
    patch_coarse_test_solver_cost(case_dir, config)


def copy_case(spec: CaseSpec, force: bool, config: dict[str, object]) -> Path:
    case_dir = RUNS_DIR / spec.case_id
    if case_dir.exists():
        if not force:
            return case_dir
        shutil.rmtree(case_dir)

    ignore = shutil.ignore_patterns(
        "0",
        "VTK",
        "VTKs",
        "processor*",
        "post-processing-data",
        "absorptivity_vs_time",
        "log.*",
        "*.foam",
        "__pycache__",
        "meltpool_calibration",
    )
    shutil.copytree(BASE_CASE, case_dir, ignore=ignore)
    patch_case(case_dir, spec, config)
    return case_dir


def run_case(case_dir: Path, run: bool, run_command: str, stream_output: bool = True) -> None:
    if not run:
        return
    if stream_output:
        subprocess.run(["bash", "-lc", run_command], cwd=case_dir, check=True, stdin=subprocess.DEVNULL)
        return

    log_path = case_dir / "log.calibrationRunner"
    with log_path.open("w") as log_file:
        subprocess.run(
            ["bash", "-lc", run_command],
            cwd=case_dir,
            check=True,
            stdin=subprocess.DEVNULL,
            stdout=log_file,
            stderr=subprocess.STDOUT,
        )


def peak_metrics(case_dir: Path, config: dict[str, object]) -> dict[str, float | str]:
    csv_path = case_dir / "post-processing-data" / "vtu_meltpool_geometry.csv"
    if not csv_path.exists():
        raise FileNotFoundError(f"Missing analyzer CSV: {csv_path}")

    df = pd.read_csv(csv_path)
    if df.empty:
        raise ValueError(f"Analyzer CSV is empty: {csv_path}")

    first_pulse_end_s = float(config["first_pulse_end_s"])
    first_pulse = df[(df["time"] >= 0.0) & (df["time"] <= first_pulse_end_s + 1e-12)].copy()
    if first_pulse.empty:
        first_pulse = df.copy()

    first_pulse["meltPoolWidth_um"] = pd.to_numeric(first_pulse["meltPoolWidth_um"], errors="coerce")
    first_pulse["meltPoolDepth_um"] = pd.to_numeric(first_pulse["meltPoolDepth_um"], errors="coerce")

    width_peak = float(first_pulse["meltPoolWidth_um"].max(skipna=True))
    depth_peak = float(first_pulse["meltPoolDepth_um"].max(skipna=True))

    width_time = float(first_pulse.loc[first_pulse["meltPoolWidth_um"].idxmax(), "time"]) if math.isfinite(width_peak) else math.nan
    depth_time = float(first_pulse.loc[first_pulse["meltPoolDepth_um"].idxmax(), "time"]) if math.isfinite(depth_peak) else math.nan

    target_width_um = float(config["target_width_um"])
    target_depth_um = float(config["target_depth_um"])
    within_percent = float(config.get("within_percent", 5.0))

    width_error_pct = 100.0 * (width_peak - target_width_um) / target_width_um
    depth_error_pct = 100.0 * (depth_peak - target_depth_um) / target_depth_um
    score = math.hypot(width_error_pct, depth_error_pct)
    if not math.isfinite(score):
        score = math.inf
    within_tolerance = (
        math.isfinite(width_error_pct)
        and math.isfinite(depth_error_pct)
        and abs(width_error_pct) <= within_percent
        and abs(depth_error_pct) <= within_percent
    )

    return {
        "peak_width_um": width_peak,
        "peak_depth_um": depth_peak,
        "peak_width_time_s": width_time,
        "peak_depth_time_s": depth_time,
        "width_error_pct": width_error_pct,
        "depth_error_pct": depth_error_pct,
        "score": score,
        "within_tolerance": str(within_tolerance),
    }


def candidate_specs(config: dict[str, object]) -> list[CaseSpec]:
    radii_um = [float(value) for value in config["laser_radius_um"]]
    resistivities = [float(value) for value in config["elec_resistivity_ohm_m"]]
    return [CaseSpec(radius * 1e-6, er) for radius in radii_um for er in resistivities]


def adaptive_config(config: dict[str, object]) -> dict[str, object]:
    return dict(config.get("adaptive", {}))


def load_existing_summary() -> dict[str, dict[str, str]]:
    if not SUMMARY_CSV.exists():
        return {}
    with SUMMARY_CSV.open(newline="") as handle:
        rows = {}
        for row in csv.DictReader(handle):
            case_id = row.get("case_id", "").strip()
            if not case_id:
                continue
            rows[case_id] = row
        return rows


def numeric_row(row: dict[str, object], key: str, default: float = math.nan) -> float:
    try:
        return float(row[key])
    except (KeyError, TypeError, ValueError):
        return default


def plot_scatter(ax, radius: np.ndarray, resistivity_micro: np.ndarray, values: np.ndarray, title: str, label: str):
    sc = ax.scatter(radius, resistivity_micro, c=values, s=95, cmap="viridis", edgecolor="black", linewidth=0.6)
    ax.set_title(title)
    ax.set_xlabel("Laser radius (um)")
    ax.set_ylabel("Electrical resistivity (x1e-6 ohm m)")
    ax.grid(True, alpha=0.25)
    return sc


def write_plots(rows: list[dict[str, object]], config: dict[str, object]) -> None:
    if not bool(config.get("plot_results", True)):
        return

    finite_rows = []
    for row in rows:
        values = [
            numeric_row(row, "laserRadius_um"),
            numeric_row(row, "elec_resistivity_ohm_m"),
            numeric_row(row, "peak_width_um"),
            numeric_row(row, "peak_depth_um"),
            numeric_row(row, "width_error_pct"),
            numeric_row(row, "depth_error_pct"),
            numeric_row(row, "score"),
        ]
        if all(math.isfinite(value) for value in values):
            finite_rows.append(row)

    if not finite_rows:
        return

    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    PLOTS_DIR.mkdir(exist_ok=True)
    finite_rows = sorted(finite_rows, key=lambda row: numeric_row(row, "score"))
    best = finite_rows[0]

    radius = np.array([numeric_row(row, "laserRadius_um") for row in finite_rows])
    resistivity_micro = np.array([numeric_row(row, "elec_resistivity_ohm_m") * 1e6 for row in finite_rows])
    score = np.array([numeric_row(row, "score") for row in finite_rows])
    width_error = np.array([numeric_row(row, "width_error_pct") for row in finite_rows])
    depth_error = np.array([numeric_row(row, "depth_error_pct") for row in finite_rows])

    fig, ax = plt.subplots(figsize=(8, 6), constrained_layout=True)
    sc = plot_scatter(ax, radius, resistivity_micro, score, "Calibration Score Map", "score")
    ax.scatter(
        [numeric_row(best, "laserRadius_um")],
        [numeric_row(best, "elec_resistivity_ohm_m") * 1e6],
        marker="*",
        s=260,
        color="red",
        edgecolor="black",
        label="Best",
        zorder=5,
    )
    fig.colorbar(sc, ax=ax, label="Score: lower is better")
    ax.legend()
    fig.savefig(PLOTS_DIR / "score_map.png", dpi=180)
    plt.close(fig)

    fig, axes = plt.subplots(1, 2, figsize=(13, 5.5), constrained_layout=True)
    sc0 = plot_scatter(axes[0], radius, resistivity_micro, width_error, "Width Error", "width error")
    sc1 = plot_scatter(axes[1], radius, resistivity_micro, depth_error, "Depth Error", "depth error")
    fig.colorbar(sc0, ax=axes[0], label="Width error (%)")
    fig.colorbar(sc1, ax=axes[1], label="Depth error (%)")
    fig.savefig(PLOTS_DIR / "error_maps.png", dpi=180)
    plt.close(fig)

    target_width = float(config["target_width_um"])
    target_depth = float(config["target_depth_um"])
    best_width = numeric_row(best, "peak_width_um")
    best_depth = numeric_row(best, "peak_depth_um")
    fig, ax = plt.subplots(figsize=(7, 5), constrained_layout=True)
    x = np.arange(2)
    width = 0.36
    ax.bar(x - width / 2, [target_width, target_depth], width, label="Experiment", color="#F58518")
    ax.bar(x + width / 2, [best_width, best_depth], width, label="Best simulation", color="#4C78A8")
    ax.set_xticks(x)
    ax.set_xticklabels(["Width", "Depth"])
    ax.set_ylabel("Length (um)")
    ax.set_title("Best Case vs Experiment")
    ax.legend()
    ax.grid(True, axis="y", alpha=0.25)
    for xpos, value in zip([x[0] - width / 2, x[1] - width / 2, x[0] + width / 2, x[1] + width / 2], [target_width, target_depth, best_width, best_depth]):
        ax.text(xpos, value, f"{value:.1f}", ha="center", va="bottom", fontsize=8)
    fig.savefig(PLOTS_DIR / "best_vs_experiment.png", dpi=180)
    plt.close(fig)

    top = finite_rows[: min(15, len(finite_rows))]
    labels = [str(row["case_id"]) for row in top]
    scores = [numeric_row(row, "score") for row in top]
    fig, ax = plt.subplots(figsize=(10, max(4.5, 0.35 * len(top))), constrained_layout=True)
    y = np.arange(len(top))
    ax.barh(y, scores, color="#4C78A8")
    ax.set_yticks(y)
    ax.set_yticklabels(labels, fontsize=8)
    ax.invert_yaxis()
    ax.set_xlabel("Score: lower is better")
    ax.set_title("Top Calibration Cases")
    ax.grid(True, axis="x", alpha=0.25)
    fig.savefig(PLOTS_DIR / "top_cases.png", dpi=180)
    plt.close(fig)


def write_summary(rows: list[dict[str, object]], config: dict[str, object] | None = None) -> None:
    cleaned_rows = []
    for row in rows:
        if not str(row.get("case_id", "")).strip():
            continue
        score = numeric_row(row, "score", math.inf)
        row["score"] = score if math.isfinite(score) else math.inf
        cleaned_rows.append(row)
    rows = sorted(cleaned_rows, key=lambda row: numeric_row(row, "score", math.inf))
    fieldnames = [
        "case_id",
        "laserRadius_m",
        "laserRadius_um",
        "elec_resistivity_ohm_m",
        "target_width_um",
        "target_depth_um",
        "peak_width_um",
        "peak_depth_um",
        "peak_width_time_s",
        "peak_depth_time_s",
        "width_error_pct",
        "depth_error_pct",
        "score",
        "within_tolerance",
        "case_dir",
    ]
    with SUMMARY_CSV.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)

    if rows:
        best = rows[0]
        BEST_TXT.write_text(
            "\n".join(
                [
                    f"case_id: {best['case_id']}",
                    f"laserRadius_m: {best['laserRadius_m']}",
                    f"laserRadius_um: {best['laserRadius_um']}",
                    f"elec_resistivity_ohm_m: {best['elec_resistivity_ohm_m']}",
                    f"peak_width_um: {best['peak_width_um']}",
                    f"peak_depth_um: {best['peak_depth_um']}",
                    f"width_error_pct: {best['width_error_pct']}",
                    f"depth_error_pct: {best['depth_error_pct']}",
                    f"score: {best['score']}",
                    f"within_tolerance: {best['within_tolerance']}",
                    f"case_dir: {best['case_dir']}",
                ]
            )
            + "\n"
        )
    if config is not None:
        write_plots(rows, config)


def best_completed_row(rows: list[dict[str, object]]) -> dict[str, object] | None:
    finite_rows = []
    for row in rows:
        try:
            score = float(row["score"])
        except (KeyError, TypeError, ValueError):
            continue
        if math.isfinite(score):
            finite_rows.append(row)
    if not finite_rows:
        return None
    return min(finite_rows, key=lambda row: float(row["score"]))


def spec_key(spec: CaseSpec) -> tuple[float, float]:
    return (round(spec.radius_m * 1e6, 6), round(spec.resistivity_ohm_m, 15))


def completed_or_planned_keys(rows: list[dict[str, object]], specs: list[CaseSpec]) -> set[tuple[float, float]]:
    keys = {spec_key(spec) for spec in specs}
    for row in rows:
        try:
            keys.add((round(float(row["laserRadius_um"]), 6), round(float(row["elec_resistivity_ohm_m"]), 15)))
        except (KeyError, TypeError, ValueError):
            continue
    return keys


def feature_matrix(radius_um: np.ndarray, resistivity_ohm_m: np.ndarray, quadratic: bool) -> np.ndarray:
    r = (radius_um - np.mean(radius_um)) / max(np.std(radius_um), 1.0)
    e_micro = resistivity_ohm_m * 1e6
    e = (e_micro - np.mean(e_micro)) / max(np.std(e_micro), 1e-9)
    columns = [np.ones_like(r), r, e]
    if quadratic:
        columns.extend([r * r, e * e, r * e])
    return np.column_stack(columns)


def fit_predictor(rows: list[dict[str, object]]):
    data = []
    for row in rows:
        try:
            item = (
                float(row["laserRadius_um"]),
                float(row["elec_resistivity_ohm_m"]),
                float(row["peak_width_um"]),
                float(row["peak_depth_um"]),
            )
        except (KeyError, TypeError, ValueError):
            continue
        if all(math.isfinite(value) for value in item):
            data.append(item)

    if len(data) < 3:
        return None

    arr = np.asarray(data, dtype=float)
    radius_um = arr[:, 0]
    resistivity = arr[:, 1]
    width = arr[:, 2]
    depth = arr[:, 3]
    quadratic = len(data) >= 8

    x = feature_matrix(radius_um, resistivity, quadratic)
    width_coef, *_ = np.linalg.lstsq(x, width, rcond=None)
    depth_coef, *_ = np.linalg.lstsq(x, depth, rcond=None)
    radius_mean = float(np.mean(radius_um))
    radius_std = float(max(np.std(radius_um), 1.0))
    er_micro = resistivity * 1e6
    er_mean = float(np.mean(er_micro))
    er_std = float(max(np.std(er_micro), 1e-9))

    def predict(candidate_radius_um: np.ndarray, candidate_er: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
        r = (candidate_radius_um - radius_mean) / radius_std
        e = (candidate_er * 1e6 - er_mean) / er_std
        cols = [np.ones_like(r), r, e]
        if quadratic:
            cols.extend([r * r, e * e, r * e])
        x_pred = np.column_stack(cols)
        return x_pred @ width_coef, x_pred @ depth_coef

    return predict


def local_fallback_specs(
    best: dict[str, object],
    rows: list[dict[str, object]],
    planned: list[CaseSpec],
    config: dict[str, object],
) -> list[CaseSpec]:
    adaptive = adaptive_config(config)
    radius_step = float(adaptive.get("laser_radius_search_step_um", 10.0))
    er_step = float(adaptive.get("elec_resistivity_search_step_ohm_m", 1e-7))
    radius_min, radius_max = [float(v) for v in adaptive.get("laser_radius_bounds_um", [200, 460])]
    er_min, er_max = [float(v) for v in adaptive.get("elec_resistivity_bounds_ohm_m", [5e-7, 2.5e-6])]

    center_radius = float(best["laserRadius_um"])
    center_er = float(best["elec_resistivity_ohm_m"])
    used = completed_or_planned_keys(rows, planned)
    candidates = []
    for dr in (-radius_step, 0.0, radius_step):
        for de in (-er_step, 0.0, er_step):
            if dr == 0.0 and de == 0.0:
                continue
            radius = min(max(center_radius + dr, radius_min), radius_max)
            er = min(max(center_er + de, er_min), er_max)
            spec = CaseSpec(radius * 1e-6, er)
            if spec_key(spec) not in used:
                candidates.append(spec)
    return candidates


def propose_adaptive_specs(
    rows: list[dict[str, object]],
    planned: list[CaseSpec],
    config: dict[str, object],
) -> list[CaseSpec]:
    best = best_completed_row(rows)
    if best is None:
        return []

    adaptive = adaptive_config(config)
    count = int(adaptive.get("candidates_per_iteration", 1))
    radius_min, radius_max = [float(v) for v in adaptive.get("laser_radius_bounds_um", [200, 460])]
    er_min, er_max = [float(v) for v in adaptive.get("elec_resistivity_bounds_ohm_m", [5e-7, 2.5e-6])]
    radius_step = float(adaptive.get("laser_radius_search_step_um", 5.0))
    er_step = float(adaptive.get("elec_resistivity_search_step_ohm_m", 5e-8))

    radius_values = np.arange(radius_min, radius_max + 0.5 * radius_step, radius_step)
    er_values = np.arange(er_min, er_max + 0.5 * er_step, er_step)
    radius_grid, er_grid = np.meshgrid(radius_values, er_values, indexing="ij")
    radius_flat = radius_grid.ravel()
    er_flat = er_grid.ravel()

    used = completed_or_planned_keys(rows, planned)
    unused_mask = np.array(
        [spec_key(CaseSpec(r * 1e-6, e)) not in used for r, e in zip(radius_flat, er_flat)],
        dtype=bool,
    )
    if not np.any(unused_mask):
        return []

    predictor = fit_predictor(rows)
    if predictor is None:
        fallback = local_fallback_specs(best, rows, planned, config)
        return fallback[:count]

    pred_width, pred_depth = predictor(radius_flat, er_flat)
    target_width = float(config["target_width_um"])
    target_depth = float(config["target_depth_um"])
    width_error = 100.0 * (pred_width - target_width) / target_width
    depth_error = 100.0 * (pred_depth - target_depth) / target_depth
    predicted_score = np.hypot(width_error, depth_error)
    predicted_score[~unused_mask] = np.inf

    order = np.argsort(predicted_score)
    specs = []
    for index in order:
        if not math.isfinite(float(predicted_score[index])):
            continue
        spec = CaseSpec(float(radius_flat[index]) * 1e-6, float(er_flat[index]))
        if spec_key(spec) in used:
            continue
        specs.append(spec)
        used.add(spec_key(spec))
        if len(specs) >= count:
            break

    if not specs:
        specs = local_fallback_specs(best, rows, planned, config)[:count]
    return specs


def execute_case(
    spec: CaseSpec,
    config: dict[str, object],
    run: bool,
    run_command: str,
    force: bool,
    stream_output: bool = True,
) -> dict[str, object] | None:
    print(f"Preparing {spec.case_id}")
    case_dir = copy_case(spec, force=force, config=config)
    print(f"Running {case_dir}" if run else f"Prepared {case_dir}")
    if run and not stream_output:
        print(f"  log: {case_dir / 'log.calibrationRunner'}")
    run_case(case_dir, run=run, run_command=run_command, stream_output=stream_output)

    if not run:
        return None

    metrics = peak_metrics(case_dir, config)
    return {
        "case_id": spec.case_id,
        "laserRadius_m": spec.radius_m,
        "laserRadius_um": spec.radius_m * 1e6,
        "elec_resistivity_ohm_m": spec.resistivity_ohm_m,
        "target_width_um": float(config["target_width_um"]),
        "target_depth_um": float(config["target_depth_um"]),
        **metrics,
        "case_dir": case_dir,
    }


def add_result_row(
    rows: list[dict[str, object]],
    row: dict[str, object],
    config: dict[str, object],
) -> list[dict[str, object]]:
    rows = [old for old in rows if old["case_id"] != row["case_id"]]
    rows.append(row)
    write_summary(rows, config)
    print(
        "Result "
        f"{row['case_id']} "
        f"width={float(row['peak_width_um']):.3f} um "
        f"depth={float(row['peak_depth_um']):.3f} um "
        f"score={float(row['score']):.3f}"
    )
    return rows


def record_case(
    spec: CaseSpec,
    rows: list[dict[str, object]],
    config: dict[str, object],
    run: bool,
    run_command: str,
    force: bool,
) -> list[dict[str, object]]:
    row = execute_case(spec, config, run, run_command, force, stream_output=True)
    if row is None:
        return rows
    return add_result_row(rows, row, config)


def run_specs(
    specs: list[CaseSpec],
    rows: list[dict[str, object]],
    config: dict[str, object],
    run: bool,
    run_command: str,
    force: bool,
    jobs: int,
) -> list[dict[str, object]]:
    if jobs <= 1 or not run:
        for spec in specs:
            rows = record_case(spec, rows, config, run, run_command, force)
        return rows

    print(f"Running {len(specs)} cases with {jobs} parallel jobs")
    with concurrent.futures.ThreadPoolExecutor(max_workers=jobs) as executor:
        futures = {
            executor.submit(execute_case, spec, config, run, run_command, force, False): spec
            for spec in specs
        }
        for future in concurrent.futures.as_completed(futures):
            spec = futures[future]
            try:
                row = future.result()
            except Exception as exc:
                print(f"ERROR: {spec.case_id} failed. See its run directory logs.")
                raise exc
            if row is not None:
                rows = add_result_row(rows, row, config)
    return rows


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", type=Path, default=DEFAULT_CONFIG)
    parser.add_argument("--run", action="store_true", help="Run OpenFOAM cases after creating them.")
    parser.add_argument(
        "--run-command",
        default=DEFAULT_RUN_COMMAND,
        help="Shell command used inside each case when --run is set.",
    )
    parser.add_argument("--force", action="store_true", help="Delete and recreate case directories.")
    parser.add_argument("--max-cases", type=int, default=None)
    parser.add_argument("--only-missing", action="store_true", help="Skip cases already present in summary.csv.")
    parser.add_argument("--no-adaptive", action="store_true", help="Run only initial config grid.")
    parser.add_argument("--plot-only", action="store_true", help="Regenerate plots from summary.csv and exit.")
    parser.add_argument("--coarse-test", action="store_true", help="Use coarse test settings from config.")
    parser.add_argument("--jobs", type=int, default=None, help="Number of case simulations to run in parallel.")
    args = parser.parse_args()
    config = effective_config(load_config(args.config), force_coarse=args.coarse_test)
    if args.run_command == DEFAULT_RUN_COMMAND:
        args.run_command = str(config.get("run_command", DEFAULT_RUN_COMMAND))
    jobs = args.jobs if args.jobs is not None else int(config.get("max_parallel_cases", 1))
    jobs = max(1, jobs)

    RUNS_DIR.mkdir(parents=True, exist_ok=True)
    existing = load_existing_summary()
    rows: list[dict[str, object]] = list(existing.values())
    if args.plot_only:
        write_plots(rows, config)
        print(f"Wrote plots to {PLOTS_DIR}")
        return

    specs = candidate_specs(config)
    if args.max_cases is not None:
        specs = specs[: args.max_cases]

    specs_to_run = [
        spec for spec in specs
        if not (args.only_missing and spec.case_id in existing)
    ]
    rows = run_specs(specs_to_run, rows, config, args.run, args.run_command, args.force, jobs)

    adaptive = adaptive_config(config)
    if args.run and not args.no_adaptive and bool(adaptive.get("enabled", False)):
        max_new_cases = int(adaptive.get("max_new_cases", 0))
        planned: list[CaseSpec] = []
        adaptive_cases_started = 0
        while adaptive_cases_started < max_new_cases:
            best = best_completed_row(rows)
            if best is not None and str(best.get("within_tolerance", "False")) == "True":
                print("Stopping: best case is within tolerance.")
                break

            proposals = propose_adaptive_specs(rows, planned, config)
            if not proposals:
                print("Stopping: no new adaptive proposal available.")
                break

            remaining = max_new_cases - adaptive_cases_started
            proposals = proposals[:remaining]
            planned.extend(proposals)
            rows = run_specs(proposals, rows, config, True, args.run_command, args.force, jobs)
            adaptive_cases_started += len(proposals)

    if not args.run:
        print("Cases prepared. Re-run with --run to execute and score.")


if __name__ == "__main__":
    main()
