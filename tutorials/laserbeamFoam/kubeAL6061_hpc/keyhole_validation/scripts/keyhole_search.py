#!/usr/bin/env python3
from __future__ import annotations

import argparse
import concurrent.futures
import csv
import json
import math
import random
import re
import shutil
import subprocess
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable

import numpy as np
import pandas as pd


WORK_DIR = Path(__file__).resolve().parents[1]
DEFAULT_CONFIG = WORK_DIR / "config.json"
PLANNED_CASES_JSON = WORK_DIR / "planned_cases.json"
_PYPLOT = None
_PYPLOT_ERROR: Exception | None = None


def pyplot():
    global _PYPLOT, _PYPLOT_ERROR
    if _PYPLOT is not None:
        return _PYPLOT
    if _PYPLOT_ERROR is not None:
        return None
    try:
        import matplotlib

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except Exception as exc:
        _PYPLOT_ERROR = exc
        print(f"WARNING: matplotlib unavailable; plot generation skipped: {exc}")
        return None
    _PYPLOT = plt
    return _PYPLOT


@dataclass(frozen=True)
class CaseSpec:
    name: str
    laser_power_W: float
    laser_radius_um: float
    e_num_density: float
    elec_resistivity_ohm_m: float
    stage: str
    iteration: int = 0

    @property
    def case_id(self) -> str:
        er_micro = self.elec_resistivity_ohm_m * 1e6
        e_scaled = self.e_num_density / 1e29
        return (
            f"{self.stage}_i{self.iteration:02d}_{self.name}"
            f"_p{self.laser_power_W:06.2f}"
            f"_r{self.laser_radius_um:05.2f}"
            f"_ne{e_scaled:05.2f}e29"
            f"_er{er_micro:06.3f}e-6"
        )


def replace_scalar(path: Path, name: str, value: float | int) -> None:
    text = path.read_text()
    pattern = re.compile(rf"(^\s*{re.escape(name)}\s+)([^;]+)(;.*$)", re.MULTILINE)
    text, count = pattern.subn(rf"\g<1>{value:.12g}\g<3>", text, count=1)
    if count != 1:
        raise ValueError(f"Could not patch {name} in {path}")
    path.write_text(text)


def patch_power_history(case_dir: Path, power_W: float, laser_on_s: float, end_s: float) -> None:
    path = case_dir / "constant" / "timeVsLaserPower"
    if laser_on_s >= end_s:
        rows = [
            (0.0, power_W),
            (end_s, power_W),
        ]
    else:
        rows = [
            (0.0, power_W),
            (laser_on_s, power_W),
            (min(laser_on_s + 1.0e-7, end_s), 0.0),
            (end_s, 0.0),
        ]

    path.write_text("(\n" + "".join(f"    ({time:.12g}       {power:.12g})\n" for time, power in rows) + ")\n")


def patch_position_history(case_dir: Path, end_s: float) -> None:
    path = case_dir / "constant" / "timeVsLaserPosition"
    path.write_text(
        "(\n"
        "    (0        (0.4e-3 0.2e-3 0.4e-3))\n"
        f"    ({end_s:.12g}   (0.4e-3 0.2e-3 0.4e-3))\n"
        ")\n"
    )


def patch_block_mesh_size(case_dir: Path, mesh_size_um: float) -> tuple[int, int, int]:
    block_mesh = case_dir / "system" / "blockMeshDict"
    text = block_mesh.read_text()
    vertex_pattern = re.compile(r"\(\s*([0-9.eE+-]+)\s+([0-9.eE+-]+)\s+([0-9.eE+-]+)\s*\)")
    vertices = [
        tuple(float(match.group(i)) for i in range(1, 4))
        for match in vertex_pattern.finditer(text)
    ][:8]
    if len(vertices) != 8:
        raise ValueError(f"Could not parse 8 vertices from {block_mesh}")

    xs = [point[0] for point in vertices]
    ys = [point[1] for point in vertices]
    zs = [point[2] for point in vertices]
    mesh_size_m = mesh_size_um * 1e-6
    nx = max(1, int(round((max(xs) - min(xs)) / mesh_size_m)))
    ny = max(1, int(round((max(ys) - min(ys)) / mesh_size_m)))
    nz = max(1, int(round((max(zs) - min(zs)) / mesh_size_m)))

    block_pattern = re.compile(
        r"(hex\s*\(0\s+1\s+2\s+3\s+4\s+5\s+6\s+7\)\s*\(\s*)(\d+)(\s+)(\d+)(\s+)(\d+)(\s*\))"
    )
    text, count = block_pattern.subn(rf"\g<1>{nx}\g<3>{ny}\g<5>{nz}\g<7>", text, count=1)
    if count != 1:
        raise ValueError(f"Could not patch block cell counts in {block_mesh}")
    text = re.sub(
        r"Base mesh spacing is uniform:.*",
        f"Base mesh spacing is uniform: {mesh_size_um:g} um in x, y, and z.",
        text,
    )
    block_mesh.write_text(text)
    return nx, ny, nz


def patch_dynamic_mesh(case_dir: Path, config: dict) -> None:
    dynamic_mesh = case_dir / "constant" / "dynamicMeshDict"
    text = dynamic_mesh.read_text()
    mesh_cfg = config["mesh"]
    mesh_type = "dynamicRefineFvMesh" if bool(mesh_cfg.get("dynamic_mesh", True)) else "staticFvMesh"
    text = re.sub(
        r"(^\s*dynamicFvMesh\s+)([^;]+)(;.*$)",
        rf"\g<1>{mesh_type}\g<3>",
        text,
        count=1,
        flags=re.MULTILINE,
    )
    text = re.sub(
        r"(^\s*maxRefinement\s+)([^;]+)(;.*$)",
        rf"\g<1>{int(mesh_cfg.get('max_refinement', 1))}\g<3>",
        text,
        count=1,
        flags=re.MULTILINE,
    )
    text = re.sub(
        r"(^\s*maxCells\s+)([^;]+)(;.*$)",
        rf"\g<1>{int(mesh_cfg.get('max_cells', 2000000))}\g<3>",
        text,
        count=1,
        flags=re.MULTILINE,
    )
    dynamic_mesh.write_text(text)


def copy_case(spec: CaseSpec, config: dict, force: bool) -> Path:
    runs_dir = WORK_DIR / str(config.get("runs_dir", "runs"))
    template = WORK_DIR / str(config["template_case"])
    case_dir = runs_dir / spec.case_id

    if case_dir.exists():
        if not force:
            return case_dir
        shutil.rmtree(case_dir)

    ignore = shutil.ignore_patterns(
        "0",
        "0.*",
        "processor*",
        "VTK",
        "VTKs",
        "postProcessing",
        "post-processing-data",
        "absorptivity_vs_time",
        "log.*",
        "*.foam",
        "__pycache__",
    )
    shutil.copytree(template, case_dir, ignore=ignore)
    shutil.copy2(WORK_DIR / "run_single_case.sh", case_dir / "run_single_case.sh")

    paper = config["paper_simulation_setup"]
    mesh = config["mesh"]
    parallel = config["parallel"]

    replace_scalar(case_dir / "constant" / "LaserProperties", "laserRadius", spec.laser_radius_um * 1e-6)
    replace_scalar(case_dir / "constant" / "LaserProperties", "e_num_density", spec.e_num_density)
    replace_scalar(case_dir / "constant" / "transportProperties", "elec_resistivity", spec.elec_resistivity_ohm_m)
    replace_scalar(case_dir / "system" / "decomposeParDict", "numberOfSubdomains", int(parallel["ranks_per_case"]))
    replace_scalar(case_dir / "system" / "controlDict", "endTime", float(paper["simulation_end_s"]))
    replace_scalar(case_dir / "system" / "controlDict", "writeInterval", float(config["write_interval_s"]))
    patch_power_history(case_dir, spec.laser_power_W, float(paper["laser_on_time_s"]), float(paper["simulation_end_s"]))
    patch_position_history(case_dir, float(paper["simulation_end_s"]))
    patch_block_mesh_size(case_dir, float(mesh["base_mesh_um"]))
    patch_dynamic_mesh(case_dir, config)

    metadata = {
        "case_id": spec.case_id,
        "stage": spec.stage,
        "iteration": spec.iteration,
        "laser_power_W": spec.laser_power_W,
        "laser_radius_um": spec.laser_radius_um,
        "e_num_density": spec.e_num_density,
        "elec_resistivity_ohm_m": spec.elec_resistivity_ohm_m,
    }
    (case_dir / "case_parameters.json").write_text(json.dumps(metadata, indent=2) + "\n")
    if "_config_path" in config:
        (case_dir / ".keyhole_config").write_text(str(config["_config_path"]) + "\n")
    return case_dir


def exploration_specs(config: dict) -> list[CaseSpec]:
    specs = []
    for item in config["exploration_cases"]:
        specs.append(
            CaseSpec(
                name=str(item["name"]),
                laser_power_W=float(item["laser_power_W"]),
                laser_radius_um=float(item["laser_radius_um"]),
                e_num_density=float(item["e_num_density"]),
                elec_resistivity_ohm_m=float(item["elec_resistivity_ohm_m"]),
                stage="sensitivity",
                iteration=0,
            )
        )
    return specs


def run_case(case_dir: Path) -> tuple[Path, str]:
    subprocess.run(["bash", "run_single_case.sh"], cwd=case_dir, check=True, stdin=subprocess.DEVNULL)
    return case_dir, "ok"


def read_case_params(case_dir: Path) -> dict:
    path = case_dir / "case_parameters.json"
    if not path.exists():
        return {}
    return json.loads(path.read_text())


def longest_stable_window(time_ms: np.ndarray, stable: np.ndarray, max_gap_ms: float) -> tuple[float, float, float]:
    best_start = math.nan
    best_end = math.nan
    best_duration = 0.0
    start_index = None

    for index, is_stable in enumerate(stable):
        if is_stable and start_index is None:
            start_index = index

        gap_break = (
            start_index is not None
            and index > start_index
            and time_ms[index] - time_ms[index - 1] > max_gap_ms
        )
        end_break = (not is_stable) or gap_break

        if start_index is not None and end_break:
            end_index = index - 1
            duration = float(time_ms[end_index] - time_ms[start_index])
            if duration > best_duration:
                best_duration = duration
                best_start = float(time_ms[start_index])
                best_end = float(time_ms[end_index])
            start_index = index if is_stable and gap_break else None

    if start_index is not None:
        duration = float(time_ms[-1] - time_ms[start_index])
        if duration > best_duration:
            best_duration = duration
            best_start = float(time_ms[start_index])
            best_end = float(time_ms[-1])

    return best_start, best_end, best_duration


def interval_error(value: float, low: float, high: float) -> float:
    if not np.isfinite(value):
        return 1.0e6
    if low <= value <= high:
        return 0.0
    return min(abs(value - low), abs(value - high))


def dimension_gate_config(config: dict) -> dict[str, object]:
    target = config["target_stability"]
    return dict(target.get("dimension_gate", {}))


def dimension_gate_penalty(
    avg_width: float,
    avg_depth: float,
    peak_width: float,
    peak_depth: float,
    peak_keyhole_depth: float,
    config: dict,
) -> dict[str, object]:
    gate = dimension_gate_config(config)
    if not bool(gate.get("enabled", False)):
        return {
            "dimension_gate_pass": True,
            "dimension_gate_penalty": 0.0,
            "hard_runaway": False,
            "hard_runaway_penalty": 0.0,
            "keyhole_depth_clipped": False,
        }

    w_low, w_high = [float(value) for value in gate.get("reasonable_width_um", [0.0, math.inf])]
    d_low, d_high = [float(value) for value in gate.get("reasonable_depth_um", [0.0, math.inf])]
    w_center = max(0.5 * (w_low + w_high), 1.0)
    d_center = max(0.5 * (d_low + d_high), 1.0)
    width_penalty = interval_error(avg_width, w_low, w_high) / w_center
    depth_penalty = interval_error(avg_depth, d_low, d_high) / d_center
    hard_width = float(gate.get("hard_width_um", math.inf))
    hard_depth = float(gate.get("hard_depth_um", math.inf))
    hard_keyhole = float(gate.get("hard_keyhole_depth_um", math.inf))
    hard_runaway = (
        (np.isfinite(peak_width) and peak_width >= hard_width)
        or (np.isfinite(peak_depth) and peak_depth >= hard_depth)
        or (np.isfinite(peak_keyhole_depth) and peak_keyhole_depth >= hard_keyhole)
    )
    clipped = np.isfinite(peak_keyhole_depth) and peak_keyhole_depth >= hard_keyhole
    penalty = float(width_penalty + depth_penalty)
    return {
        "dimension_gate_pass": penalty <= 1.0e-12 and not hard_runaway,
        "dimension_gate_penalty": penalty,
        "hard_runaway": bool(hard_runaway),
        "hard_runaway_penalty": 1.0 if hard_runaway else 0.0,
        "keyhole_depth_clipped": bool(clipped),
    }


def read_experimental_ar(config: dict) -> tuple[np.ndarray, np.ndarray]:
    exp_path = (WORK_DIR / str(config["exp_aspect_ratio_csv"])).resolve()
    if not exp_path.exists():
        return np.array([], dtype=float), np.array([], dtype=float)
    exp = pd.read_csv(exp_path)
    exp.columns = [column.strip() for column in exp.columns]
    exp_time = pd.to_numeric(exp["time_ms"], errors="coerce").to_numpy(dtype=float)
    exp_ar = pd.to_numeric(exp["meltpool_ar_depthByWidth"], errors="coerce").to_numpy(dtype=float)
    finite = np.isfinite(exp_time) & np.isfinite(exp_ar)
    if not np.any(finite):
        return np.array([], dtype=float), np.array([], dtype=float)
    exp_time = exp_time[finite]
    exp_ar = exp_ar[finite]
    return exp_time - float(np.nanmin(exp_time)), exp_ar


def usable_ar_mask(
    time_ms: np.ndarray,
    width: np.ndarray,
    depth: np.ndarray,
    keyhole_depth: np.ndarray,
    config: dict,
) -> np.ndarray:
    ar = np.divide(depth, width, out=np.full_like(depth, np.nan), where=width > 0)
    finite = np.isfinite(time_ms) & np.isfinite(ar)
    gate = dimension_gate_config(config)
    if bool(gate.get("ignore_ar_after_keyhole_clip", False)) and keyhole_depth.size == time_ms.size:
        hard_keyhole = float(gate.get("hard_keyhole_depth_um", math.inf))
        clipped = np.isfinite(keyhole_depth) & (keyhole_depth >= hard_keyhole)
        if np.any(clipped):
            finite &= time_ms < float(np.nanmin(time_ms[clipped]))
    return finite


def candidate_comparison_windows(
    time_ms: np.ndarray,
    width: np.ndarray,
    depth: np.ndarray,
    keyhole_depth: np.ndarray,
    config: dict,
) -> dict[str, object]:
    target = config["target_stability"]
    window_ms = float(target.get("comparison_window_ms", 2.0))
    finite = usable_ar_mask(time_ms, width, depth, keyhole_depth, config)
    if not np.any(finite):
        return []
    finite_time = time_ms[finite]
    min_time = float(np.nanmin(finite_time))
    max_time = float(np.nanmax(finite_time))
    min_coverage = float(target.get("min_comparison_window_coverage", 0.95))
    if max_time - min_time < min_coverage * window_ms:
        start_ms = min_time
        end_ms = max_time
        return [
            {
                "mask": finite & (time_ms >= start_ms) & (time_ms <= end_ms),
                "start_ms": start_ms,
                "end_ms": end_ms,
                "duration_ms": end_ms - start_ms,
            }
        ]

    valid_starts = finite_time[finite_time <= max_time - min_coverage * window_ms]
    starts = np.unique(np.round(valid_starts, decimals=9))
    windows = []
    for start_ms in starts:
        end_ms = min(float(start_ms) + window_ms, max_time)
        mask = finite & (time_ms >= float(start_ms)) & (time_ms <= end_ms)
        if np.count_nonzero(mask) < 2:
            continue
        duration_ms = float(time_ms[mask][-1] - time_ms[mask][0])
        if duration_ms < min_coverage * window_ms:
            continue
        windows.append(
            {
                "mask": mask,
                "start_ms": float(time_ms[mask][0]),
                "end_ms": float(time_ms[mask][-1]),
                "duration_ms": duration_ms,
            }
        )
    return windows


def aspect_ratio_trace_metrics(
    time_ms: np.ndarray,
    width: np.ndarray,
    depth: np.ndarray,
    keyhole_depth: np.ndarray,
    config: dict,
) -> dict[str, float]:
    windows = candidate_comparison_windows(time_ms, width, depth, keyhole_depth, config)
    if not windows:
        return {
            "ar_rmse": math.nan,
            "ar_mae": math.nan,
            "ar_bias": math.nan,
            "ar_slope_per_ms": math.nan,
            "exp_ar_slope_per_ms": math.nan,
            "ar_slope_error_per_ms": math.nan,
            "comparison_start_ms": math.nan,
            "comparison_end_ms": math.nan,
            "comparison_duration_ms": 0.0,
            "comparison_coverage": 0.0,
            "comparison_window_score": math.nan,
            "comparison_window_count": 0,
        }

    target = config["target_stability"]
    window_ms = float(target.get("comparison_window_ms", 2.0))
    exp_time_rel, exp_ar = read_experimental_ar(config)
    if exp_time_rel.size < 2:
        window = windows[-1]
        return {
            "ar_rmse": math.nan,
            "ar_mae": math.nan,
            "ar_bias": math.nan,
            "ar_slope_per_ms": math.nan,
            "exp_ar_slope_per_ms": math.nan,
            "ar_slope_error_per_ms": math.nan,
            "comparison_start_ms": float(window["start_ms"]),
            "comparison_end_ms": float(window["end_ms"]),
            "comparison_duration_ms": float(window["duration_ms"]),
            "comparison_coverage": float(window["duration_ms"]) / max(window_ms, 1.0e-12),
            "comparison_window_score": math.nan,
            "comparison_window_count": len(windows),
        }

    best: dict[str, float] | None = None
    for window in windows:
        mask = np.asarray(window["mask"], dtype=bool)
        sim_time_rel = time_ms[mask] - float(window["start_ms"])
        sim_ar = np.divide(depth[mask], width[mask], out=np.full(np.count_nonzero(mask), np.nan), where=width[mask] > 0)
        finite = np.isfinite(sim_time_rel) & np.isfinite(sim_ar)
        sim_time_rel = sim_time_rel[finite]
        sim_ar = sim_ar[finite]
        if sim_time_rel.size < 2:
            continue
        order = np.argsort(sim_time_rel)
        sim_time_rel = sim_time_rel[order]
        sim_ar = sim_ar[order]
        compare_time = exp_time_rel[(exp_time_rel >= sim_time_rel[0]) & (exp_time_rel <= sim_time_rel[-1])]
        if compare_time.size < 2:
            compare_time = np.linspace(0.0, min(window_ms, sim_time_rel[-1]), num=min(32, sim_time_rel.size))
        if compare_time.size < 2:
            continue

        exp_interp = np.interp(compare_time, exp_time_rel, exp_ar)
        sim_interp = np.interp(compare_time, sim_time_rel, sim_ar)
        diff = sim_interp - exp_interp
        sim_slope = float(np.polyfit(compare_time, sim_interp, 1)[0])
        exp_slope = float(np.polyfit(compare_time, exp_interp, 1)[0])
        rmse = float(np.sqrt(np.nanmean(diff * diff)))
        bias = float(np.nanmean(diff))
        slope_error = float(sim_slope - exp_slope)
        coverage = min(1.0, float(window["duration_ms"]) / max(window_ms, 1.0e-12))
        gate = dimension_gate_penalty(
            float(np.nanmean(width[mask])),
            float(np.nanmean(depth[mask])),
            float(np.nanmax(width[mask])),
            float(np.nanmax(depth[mask])),
            float(np.nanmax(keyhole_depth[mask])) if keyhole_depth.size == time_ms.size and np.any(np.isfinite(keyhole_depth[mask])) else math.nan,
            config,
        )
        score = (
            rmse / max(float(target.get("aspect_ratio_rmse_target", 0.15)), 1.0e-12)
            + 0.35 * abs(bias) / max(float(target.get("aspect_ratio_bias_target", 0.10)), 1.0e-12)
            + 0.35 * abs(slope_error) / max(float(target.get("aspect_ratio_slope_target_per_ms", 0.12)), 1.0e-12)
            + max(0.0, 1.0 - coverage)
            + 2.0 * float(gate["dimension_gate_penalty"])
            + 4.0 * float(gate["hard_runaway_penalty"])
        )
        result = {
            "ar_rmse": rmse,
            "ar_mae": float(np.nanmean(np.abs(diff))),
            "ar_bias": bias,
            "ar_slope_per_ms": sim_slope,
            "exp_ar_slope_per_ms": exp_slope,
            "ar_slope_error_per_ms": slope_error,
            "comparison_start_ms": float(window["start_ms"]),
            "comparison_end_ms": float(window["end_ms"]),
            "comparison_duration_ms": float(window["duration_ms"]),
            "comparison_coverage": coverage,
            "comparison_window_score": float(score),
            "comparison_window_count": len(windows),
        }
        if best is None or score < float(best["comparison_window_score"]):
            best = result

    if best is None:
        return {
            "ar_rmse": math.nan,
            "ar_mae": math.nan,
            "ar_bias": math.nan,
            "ar_slope_per_ms": math.nan,
            "exp_ar_slope_per_ms": math.nan,
            "ar_slope_error_per_ms": math.nan,
            "comparison_start_ms": math.nan,
            "comparison_end_ms": math.nan,
            "comparison_duration_ms": 0.0,
            "comparison_coverage": 0.0,
            "comparison_window_score": math.nan,
            "comparison_window_count": len(windows),
        }
    return best


def write_target_snapshot_artifacts(
    case_dir: Path,
    data: pd.DataFrame,
    config: dict,
) -> dict[str, object]:
    target = config["target_stability"]
    enabled = bool(target.get("save_width_depth_target_snapshots", True))
    result: dict[str, object] = {
        "target_box_hit_count": 0,
        "first_target_box_time_ms": math.nan,
        "last_target_box_time_ms": math.nan,
        "target_box_snapshot_dir": "",
    }
    if not enabled or data.empty:
        return result

    width = pd.to_numeric(data.get("meltPoolWidth_um"), errors="coerce")
    depth = pd.to_numeric(data.get("meltPoolDepth_um"), errors="coerce")
    time_s = pd.to_numeric(data.get("time"), errors="coerce")
    w_low, w_high = [float(value) for value in target["meltpool_width_um"]]
    d_low, d_high = [float(value) for value in target["meltpool_depth_um"]]
    hit_mask = (
        np.isfinite(time_s)
        & np.isfinite(width)
        & np.isfinite(depth)
        & (width >= w_low)
        & (width <= w_high)
        & (depth >= d_low)
        & (depth <= d_high)
    )
    hit_indices = list(np.flatnonzero(hit_mask.to_numpy(dtype=bool)))
    if not hit_indices:
        return result

    out_dir = case_dir / "post-processing-data" / "target_meltpool_snapshots"
    out_dir.mkdir(parents=True, exist_ok=True)
    section_dir = case_dir / "post-processing-data" / "vtu_sections"
    section_pngs = sorted(section_dir.glob("*.png"))
    max_pngs = int(target.get("max_target_snapshot_pngs_per_case", 60))
    selected_indices = set(hit_indices[:max_pngs])
    manifest_rows = []

    for index in hit_indices:
        source_png = section_pngs[index] if index < len(section_pngs) else None
        copied_png = ""
        if source_png is not None and index in selected_indices:
            copied = out_dir / f"target_{len(manifest_rows):04d}_t{float(time_s.iloc[index]) * 1.0e3:.6f}ms_{source_png.name}"
            if not copied.exists():
                shutil.copy2(source_png, copied)
            copied_png = str(copied.relative_to(case_dir))
        manifest_rows.append(
            {
                "case_id": case_dir.name,
                "row_index": index,
                "time_s": float(time_s.iloc[index]),
                "time_ms": float(time_s.iloc[index]) * 1.0e3,
                "meltpool_width_um": float(width.iloc[index]),
                "meltpool_depth_um": float(depth.iloc[index]),
                "depth_over_width": float(depth.iloc[index] / width.iloc[index]) if float(width.iloc[index]) > 0 else math.nan,
                "source_section_png": str(source_png.relative_to(case_dir)) if source_png is not None else "",
                "copied_section_png": copied_png,
            }
        )

    manifest_path = out_dir / "target_meltpool_snapshots.csv"
    pd.DataFrame(manifest_rows).to_csv(manifest_path, index=False)
    times_ms = [row["time_ms"] for row in manifest_rows]
    result.update(
        {
            "target_box_hit_count": len(hit_indices),
            "first_target_box_time_ms": float(min(times_ms)),
            "last_target_box_time_ms": float(max(times_ms)),
            "target_box_snapshot_dir": str(out_dir.relative_to(case_dir)),
        }
    )
    return result


def case_metrics(case_dir: Path, config: dict) -> dict[str, object]:
    params = read_case_params(case_dir)
    early_status_path = case_dir / "early_stop_status.json"
    early_status = json.loads(early_status_path.read_text()) if early_status_path.exists() else {}
    row: dict[str, object] = {
        "case_id": case_dir.name,
        "case_dir": str(case_dir),
        "status": "missing_csv",
        "early_stop_decision": early_status.get("decision", ""),
        "early_stop_reason": early_status.get("reason", ""),
        **params,
    }

    csv_path = case_dir / "post-processing-data" / "vtu_meltpool_geometry.csv"
    if not csv_path.exists():
        row["score"] = float(config["scoring_weights"]["missing_csv_penalty"])
        return row

    data = pd.read_csv(csv_path)
    if data.empty:
        row["status"] = "empty_csv"
        row["score"] = float(config["scoring_weights"]["missing_csv_penalty"])
        return row
    target_snapshot_artifacts = write_target_snapshot_artifacts(case_dir, data, config)

    time_ms = pd.to_numeric(data["time"], errors="coerce").to_numpy(dtype=float) * 1.0e3
    width = pd.to_numeric(data["meltPoolWidth_um"], errors="coerce").to_numpy(dtype=float)
    depth = pd.to_numeric(data["meltPoolDepth_um"], errors="coerce").to_numpy(dtype=float)
    keyhole_depth = pd.to_numeric(data["keyholeDepth_um"], errors="coerce").to_numpy(dtype=float)
    finite = np.isfinite(time_ms) & np.isfinite(width) & np.isfinite(depth)
    time_ms = time_ms[finite]
    width = width[finite]
    depth = depth[finite]
    keyhole_depth = keyhole_depth[finite]

    if time_ms.size == 0:
        keyhole_values = pd.to_numeric(data.get("keyholeDepth_um"), errors="coerce")
        row.update(
            {
                "status": "no_meltpool_contour",
                "final_time_ms": float(pd.to_numeric(data.get("time"), errors="coerce").max() * 1.0e3) if "time" in data else math.nan,
                "stable": False,
                "stable_start_ms": math.nan,
                "stable_end_ms": math.nan,
                "stable_duration_ms": 0.0,
                "avg_stable_width_um": math.nan,
                "avg_stable_depth_um": math.nan,
                "avg_stable_depth_over_width": math.nan,
                "avg_stable_keyhole_depth_um": float(keyhole_values.max()) if not keyhole_values.dropna().empty else math.nan,
                "peak_meltpool_width_um": math.nan,
                "peak_meltpool_depth_um": math.nan,
                "peak_keyhole_depth_um": float(keyhole_values.max()) if not keyhole_values.dropna().empty else math.nan,
                **target_snapshot_artifacts,
                "score": float(config["scoring_weights"]["missing_csv_penalty"]),
            }
        )
        return row

    target = config["target_stability"]
    w_low, w_high = [float(value) for value in target["meltpool_width_um"]]
    d_low, d_high = [float(value) for value in target["meltpool_depth_um"]]
    ar_metrics = aspect_ratio_trace_metrics(time_ms, width, depth, keyhole_depth, config)
    window = {
        "start_ms": ar_metrics["comparison_start_ms"],
        "end_ms": ar_metrics["comparison_end_ms"],
        "duration_ms": ar_metrics["comparison_duration_ms"],
    }
    if np.isfinite(window["start_ms"]) and np.isfinite(window["end_ms"]):
        finite = usable_ar_mask(time_ms, width, depth, keyhole_depth, config)
        window["mask"] = finite & (time_ms >= float(window["start_ms"])) & (time_ms <= float(window["end_ms"]))
    else:
        window["mask"] = np.zeros_like(time_ms, dtype=bool)
    stable_mask = np.asarray(window["mask"], dtype=bool)
    stable_start = float(window["start_ms"])
    stable_end = float(window["end_ms"])
    stable_duration = float(window["duration_ms"])

    if np.any(stable_mask):
        avg_width = float(np.nanmean(width[stable_mask]))
        avg_depth = float(np.nanmean(depth[stable_mask]))
        avg_keyhole_depth = float(np.nanmean(keyhole_depth[stable_mask]))
        avg_ar_depth_width = float(avg_depth / avg_width) if avg_width > 0 else math.nan
    else:
        avg_width = float(np.nanmax(width)) if width.size else math.nan
        avg_depth = float(np.nanmax(depth)) if depth.size else math.nan
        avg_keyhole_depth = float(np.nanmax(keyhole_depth)) if keyhole_depth.size else math.nan
        avg_ar_depth_width = float(avg_depth / avg_width) if np.isfinite(avg_width) and avg_width > 0 else math.nan

    center_width = 0.5 * (w_low + w_high)
    center_depth = 0.5 * (d_low + d_high)
    weights = config["scoring_weights"]
    width_score = interval_error(avg_width, w_low, w_high) / max(1.0, center_width)
    depth_score = interval_error(avg_depth, d_low, d_high) / max(1.0, center_depth)
    finite_keyhole = keyhole_depth[np.isfinite(keyhole_depth)]
    peak_width = float(np.nanmax(width)) if width.size else math.nan
    peak_depth = float(np.nanmax(depth)) if depth.size else math.nan
    peak_keyhole = float(np.nanmax(finite_keyhole)) if finite_keyhole.size else math.nan
    gate_metrics = dimension_gate_penalty(avg_width, avg_depth, peak_width, peak_depth, peak_keyhole, config)
    ar_rmse = ar_metrics["ar_rmse"]
    ar_bias = ar_metrics["ar_bias"]
    ar_slope_error = ar_metrics["ar_slope_error_per_ms"]
    comparison_coverage = ar_metrics["comparison_coverage"]
    stable = (
        comparison_coverage >= 0.95
        and np.isfinite(ar_rmse)
        and ar_rmse <= float(target.get("aspect_ratio_rmse_target", 0.15))
        and abs(ar_bias) <= float(target.get("aspect_ratio_bias_target", 0.10))
        and abs(ar_slope_error) <= float(target.get("aspect_ratio_slope_target_per_ms", 0.12))
        and bool(gate_metrics["dimension_gate_pass"])
    )
    stable_penalty = 0.0 if stable else 1.0
    keyhole_penalty = 0.0 if np.isfinite(peak_keyhole) and peak_keyhole > 40.0 else 1.0
    ar_rmse_score = ar_rmse if np.isfinite(ar_rmse) else 10.0
    ar_bias_score = abs(ar_bias) if np.isfinite(ar_bias) else 10.0
    ar_slope_score = abs(ar_slope_error) if np.isfinite(ar_slope_error) else 10.0
    window_score = max(0.0, 1.0 - comparison_coverage)
    score = (
        float(weights["stable_window"]) * stable_penalty
        + float(weights["aspect_ratio_error"]) * ar_rmse_score
        + float(weights.get("aspect_ratio_bias", 0.0)) * ar_bias_score
        + float(weights.get("aspect_ratio_slope", 0.0)) * ar_slope_score
        + float(weights.get("comparison_window", 0.0)) * window_score
        + float(weights.get("dimension_gate", 0.0)) * float(gate_metrics["dimension_gate_penalty"])
        + float(weights.get("hard_runaway", 0.0)) * float(gate_metrics["hard_runaway_penalty"])
        + float(weights.get("keyhole_missing", 0.0)) * keyhole_penalty
        + float(weights["width_error"]) * width_score
        + float(weights["depth_error"]) * depth_score
    )

    row.update(
        {
            "status": "ok",
            "final_time_ms": float(np.nanmax(time_ms)) if time_ms.size else math.nan,
            "stable": bool(stable),
            "stable_start_ms": stable_start,
            "stable_end_ms": stable_end,
            "stable_duration_ms": stable_duration,
            "comparison_start_ms": ar_metrics["comparison_start_ms"],
            "comparison_end_ms": ar_metrics["comparison_end_ms"],
            "comparison_duration_ms": ar_metrics["comparison_duration_ms"],
            "comparison_coverage": ar_metrics["comparison_coverage"],
            "comparison_window_score": ar_metrics["comparison_window_score"],
            "comparison_window_count": ar_metrics["comparison_window_count"],
            "ar_rmse": ar_metrics["ar_rmse"],
            "ar_mae": ar_metrics["ar_mae"],
            "ar_bias": ar_metrics["ar_bias"],
            "ar_slope_per_ms": ar_metrics["ar_slope_per_ms"],
            "exp_ar_slope_per_ms": ar_metrics["exp_ar_slope_per_ms"],
            "ar_slope_error_per_ms": ar_metrics["ar_slope_error_per_ms"],
            "avg_stable_width_um": avg_width,
            "avg_stable_depth_um": avg_depth,
            "avg_stable_depth_over_width": avg_ar_depth_width,
            "avg_stable_keyhole_depth_um": avg_keyhole_depth,
            "dimension_gate_pass": bool(gate_metrics["dimension_gate_pass"]),
            "dimension_gate_penalty": float(gate_metrics["dimension_gate_penalty"]),
            "hard_runaway": bool(gate_metrics["hard_runaway"]),
            "hard_runaway_penalty": float(gate_metrics["hard_runaway_penalty"]),
            "keyhole_depth_clipped": bool(gate_metrics["keyhole_depth_clipped"]),
            "peak_meltpool_width_um": peak_width,
            "peak_meltpool_depth_um": peak_depth,
            "peak_keyhole_depth_um": peak_keyhole,
            **target_snapshot_artifacts,
            "score": float(score),
        }
    )
    return row


def write_summary(rows: list[dict[str, object]], path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    keys = []
    for row in rows:
        for key in row:
            if key not in keys:
                keys.append(key)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=keys)
        writer.writeheader()
        writer.writerows(rows)


def read_summary(path: Path) -> list[dict[str, object]]:
    if not path.exists():
        return []
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle))


def case_spec_record(spec: CaseSpec) -> dict[str, object]:
    return {
        "case_id": spec.case_id,
        "name": spec.name,
        "stage": spec.stage,
        "iteration": spec.iteration,
        "laser_power_W": spec.laser_power_W,
        "laser_radius_um": spec.laser_radius_um,
        "e_num_density": spec.e_num_density,
        "elec_resistivity_ohm_m": spec.elec_resistivity_ohm_m,
    }


def write_planned_cases(specs: list[CaseSpec], label: str) -> None:
    payload = {
        "label": label,
        "cases": [case_spec_record(spec) for spec in specs],
    }
    PLANNED_CASES_JSON.write_text(json.dumps(payload, indent=2) + "\n")


def read_planned_cases() -> dict[str, object]:
    if not PLANNED_CASES_JSON.exists():
        return {"label": "none", "cases": []}
    return json.loads(PLANNED_CASES_JSON.read_text())


def numeric(row: dict[str, object], key: str, default: float = math.nan) -> float:
    try:
        return float(row.get(key, default))
    except (TypeError, ValueError):
        return default


def best_rows(rows: list[dict[str, object]], limit: int = 5) -> list[dict[str, object]]:
    candidates = [row for row in rows if str(row.get("status")) == "ok"]
    return sorted(candidates, key=lambda row: numeric(row, "score", 1.0e9))[:limit]


def clip(value: float, bounds: list[float]) -> float:
    return max(float(bounds[0]), min(float(bounds[1]), value))


PARAMETER_KEYS = ["laser_power_W", "laser_radius_um", "e_num_density", "elec_resistivity_ohm_m"]


def parameter_matrix(rows: list[dict[str, object]]) -> np.ndarray:
    return np.array([[numeric(row, key) for key in PARAMETER_KEYS] for row in rows], dtype=float)


def row_targets(rows: list[dict[str, object]]) -> dict[str, np.ndarray]:
    stable = np.array([str(row.get("stable")).lower() == "true" for row in rows], dtype=float)
    return {
        "score": np.array([numeric(row, "score", 1.0e6) for row in rows], dtype=float),
        "width": np.array([numeric(row, "avg_stable_width_um") for row in rows], dtype=float),
        "depth": np.array([numeric(row, "avg_stable_depth_um") for row in rows], dtype=float),
        "duration": np.array([numeric(row, "stable_duration_ms", 0.0) for row in rows], dtype=float),
        "ar_rmse": np.array([numeric(row, "ar_rmse", 10.0) for row in rows], dtype=float),
        "ar_bias": np.array([abs(numeric(row, "ar_bias", 10.0)) for row in rows], dtype=float),
        "ar_slope_error": np.array([abs(numeric(row, "ar_slope_error_per_ms", 10.0)) for row in rows], dtype=float),
        "stable": stable,
    }


def random_candidate_pool(config: dict, rng: random.Random, count: int) -> np.ndarray:
    bounds = config["parameter_bounds"]
    pool = np.zeros((count, len(PARAMETER_KEYS)), dtype=float)
    for col, key in enumerate(PARAMETER_KEYS):
        low, high = [float(value) for value in bounds[key]]
        pool[:, col] = [rng.uniform(low, high) for _ in range(count)]
    return pool


def target_error_from_arrays(width: np.ndarray, depth: np.ndarray, duration: np.ndarray, config: dict) -> np.ndarray:
    target = config["target_stability"]
    w_low, w_high = [float(value) for value in target["meltpool_width_um"]]
    d_low, d_high = [float(value) for value in target["meltpool_depth_um"]]
    min_duration = float(target["min_stable_duration_ms"])
    width_center = 0.5 * (w_low + w_high)
    depth_center = 0.5 * (d_low + d_high)

    width_error = np.where(width < w_low, w_low - width, np.where(width > w_high, width - w_high, 0.0))
    depth_error = np.where(depth < d_low, d_low - depth, np.where(depth > d_high, depth - d_high, 0.0))
    duration_error = np.maximum(0.0, min_duration - duration) / max(min_duration, 1.0e-12)
    return width_error / width_center + depth_error / depth_center + duration_error


def dimension_gate_error_from_arrays(width: np.ndarray, depth: np.ndarray, config: dict) -> np.ndarray:
    gate = dimension_gate_config(config)
    if not bool(gate.get("enabled", False)):
        return np.zeros_like(width, dtype=float)
    w_low, w_high = [float(value) for value in gate.get("reasonable_width_um", [0.0, math.inf])]
    d_low, d_high = [float(value) for value in gate.get("reasonable_depth_um", [0.0, math.inf])]
    w_center = max(0.5 * (w_low + w_high), 1.0)
    d_center = max(0.5 * (d_low + d_high), 1.0)
    width_error = np.where(width < w_low, w_low - width, np.where(width > w_high, width - w_high, 0.0))
    depth_error = np.where(depth < d_low, d_low - depth, np.where(depth > d_high, depth - d_high, 0.0))
    return width_error / w_center + depth_error / d_center


def target_error_from_predictions(predictions: dict[str, np.ndarray], config: dict) -> np.ndarray:
    target = config["target_stability"]
    if str(target.get("objective", "")) == "aspect_ratio_trace" and "ar_rmse" in predictions:
        min_duration = float(target["min_stable_duration_ms"])
        duration_error = np.maximum(0.0, min_duration - predictions["duration"]) / max(min_duration, 1.0e-12)
        rmse_error = predictions["ar_rmse"] / max(float(target.get("aspect_ratio_rmse_target", 0.15)), 1.0e-12)
        bias_error = predictions["ar_bias"] / max(float(target.get("aspect_ratio_bias_target", 0.10)), 1.0e-12)
        slope_error = predictions["ar_slope_error"] / max(float(target.get("aspect_ratio_slope_target_per_ms", 0.12)), 1.0e-12)
        dimension_error = dimension_gate_error_from_arrays(predictions["width"], predictions["depth"], config)
        return rmse_error + 0.35 * bias_error + 0.35 * slope_error + duration_error + 2.0 * dimension_error
    return target_error_from_arrays(predictions["width"], predictions["depth"], predictions["duration"], config)


def completed_training_rows(rows: list[dict[str, object]]) -> list[dict[str, object]]:
    usable = []
    for row in rows:
        if str(row.get("status")) != "ok":
            continue
        values = [numeric(row, key) for key in PARAMETER_KEYS]
        values.extend(
            [
                numeric(row, "score"),
                numeric(row, "avg_stable_width_um"),
                numeric(row, "avg_stable_depth_um"),
                numeric(row, "stable_duration_ms"),
                numeric(row, "ar_rmse"),
                numeric(row, "ar_bias"),
                numeric(row, "ar_slope_error_per_ms"),
            ]
        )
        if all(np.isfinite(value) for value in values):
            usable.append(row)
    return usable


def feature_scale(config: dict) -> np.ndarray:
    bounds = config["parameter_bounds"]
    spans = []
    for key in PARAMETER_KEYS:
        low, high = [float(value) for value in bounds[key]]
        spans.append(max(high - low, 1.0e-30))
    return np.array(spans, dtype=float)


def normalized_distance(point: np.ndarray, selected: list[np.ndarray], scale: np.ndarray) -> float:
    if not selected:
        return math.inf
    distances = [float(np.linalg.norm((point - other) / scale)) for other in selected]
    return min(distances)


def select_diverse_batch(
    pool: np.ndarray,
    rank_value: np.ndarray,
    count: int,
    config: dict,
) -> np.ndarray:
    scale = feature_scale(config)
    order = np.argsort(rank_value)
    shortlist_size = max(count * 20, int(len(pool) * float(config["optimization"].get("diversity_fraction", 0.08))))
    shortlist = order[:min(len(order), shortlist_size)]
    selected_indices: list[int] = []
    selected_points: list[np.ndarray] = []

    if len(shortlist) == 0:
        return pool[:0]

    selected_indices.append(int(shortlist[0]))
    selected_points.append(pool[int(shortlist[0])])

    while len(selected_indices) < count and len(selected_indices) < len(shortlist):
        best_index = None
        best_value = -math.inf
        for candidate_index in shortlist:
            candidate_index = int(candidate_index)
            if candidate_index in selected_indices:
                continue
            diversity = normalized_distance(pool[candidate_index], selected_points, scale)
            score = -rank_value[candidate_index] + 0.15 * diversity
            if score > best_value:
                best_value = score
                best_index = candidate_index
        if best_index is None:
            break
        selected_indices.append(best_index)
        selected_points.append(pool[best_index])

    return pool[selected_indices]


def specs_from_pool(pool: np.ndarray, iteration: int) -> list[CaseSpec]:
    specs = []
    for index, values in enumerate(pool):
        specs.append(
            CaseSpec(
                name=f"opt{index:02d}",
                laser_power_W=float(values[0]),
                laser_radius_um=float(values[1]),
                e_num_density=float(values[2]),
                elec_resistivity_ohm_m=float(values[3]),
                stage="opt",
                iteration=iteration,
            )
        )
    return specs


def propose_surrogate_specs(config: dict, rows: list[dict[str, object]], iteration: int) -> list[CaseSpec] | None:
    training_rows = completed_training_rows(rows)
    opt = config["optimization"]
    min_cases = int(opt.get("surrogate_min_cases", 8))
    if len(training_rows) < min_cases:
        return None

    try:
        from sklearn.ensemble import ExtraTreesRegressor
    except Exception as exc:
        print(f"WARNING: sklearn unavailable; falling back to local/random optimizer: {exc}")
        return None

    rng = random.Random(int(opt["random_seed"]) + 1000 * iteration)
    x_train = parameter_matrix(training_rows)
    targets = row_targets(training_rows)
    pool_size = int(opt.get("candidate_pool_size", 20000))
    pool = random_candidate_pool(config, rng, pool_size)

    models = {}
    predictions = {}
    uncertainty = {}
    for name, values in targets.items():
        if name == "stable":
            continue
        model = ExtraTreesRegressor(
            n_estimators=int(opt.get("surrogate_n_estimators", 256)),
            random_state=int(opt["random_seed"]) + iteration,
            min_samples_leaf=1,
            bootstrap=False,
        )
        model.fit(x_train, values)
        tree_predictions = np.vstack([tree.predict(pool) for tree in model.estimators_])
        predictions[name] = np.mean(tree_predictions, axis=0)
        uncertainty[name] = np.std(tree_predictions, axis=0)
        models[name] = model

    predicted_error = target_error_from_predictions(predictions, config)
    target = config["target_stability"]
    w_low, w_high = [float(value) for value in target["meltpool_width_um"]]
    d_low, d_high = [float(value) for value in target["meltpool_depth_um"]]
    min_duration = float(target["min_stable_duration_ms"])
    if str(target.get("objective", "")) == "aspect_ratio_trace" and "ar_rmse" in predictions:
        gate = dimension_gate_config(config)
        if bool(gate.get("enabled", False)):
            gate_w_low, gate_w_high = [float(value) for value in gate.get("reasonable_width_um", [0.0, math.inf])]
            gate_d_low, gate_d_high = [float(value) for value in gate.get("reasonable_depth_um", [0.0, math.inf])]
        else:
            gate_w_low, gate_w_high = -math.inf, math.inf
            gate_d_low, gate_d_high = -math.inf, math.inf
        feasible = (
            (predictions["ar_rmse"] <= float(target.get("aspect_ratio_rmse_target", 0.15)))
            & (predictions["ar_bias"] <= float(target.get("aspect_ratio_bias_target", 0.10)))
            & (predictions["ar_slope_error"] <= float(target.get("aspect_ratio_slope_target_per_ms", 0.12)))
            & (predictions["duration"] >= min_duration)
            & (predictions["width"] >= gate_w_low)
            & (predictions["width"] <= gate_w_high)
            & (predictions["depth"] >= gate_d_low)
            & (predictions["depth"] <= gate_d_high)
        ).astype(float)
    else:
        feasible = (
            (predictions["width"] >= w_low)
            & (predictions["width"] <= w_high)
            & (predictions["depth"] >= d_low)
            & (predictions["depth"] <= d_high)
            & (predictions["duration"] >= min_duration)
        ).astype(float)

    score_scale = max(float(np.nanstd(targets["score"])), 1.0)
    if str(target.get("objective", "")) == "aspect_ratio_trace" and "ar_rmse" in uncertainty:
        uncertainty_bonus = (
            uncertainty["ar_rmse"] / max(float(target.get("aspect_ratio_rmse_target", 0.15)), 1.0e-12)
            + uncertainty["ar_bias"] / max(float(target.get("aspect_ratio_bias_target", 0.10)), 1.0e-12)
            + uncertainty["ar_slope_error"] / max(float(target.get("aspect_ratio_slope_target_per_ms", 0.12)), 1.0e-12)
            + uncertainty["duration"] / max(min_duration, 1.0e-12)
        )
    else:
        uncertainty_bonus = (
            uncertainty["width"] / max(w_high - w_low, 1.0)
            + uncertainty["depth"] / max(d_high - d_low, 1.0)
            + uncertainty["duration"] / max(min_duration, 1.0e-12)
        )
    rank_value = (
        predictions["score"] / score_scale
        + predicted_error
        - float(opt.get("exploration_weight", 0.2)) * uncertainty_bonus
        - float(opt.get("feasibility_weight", 2.0)) * feasible
    )

    selected = select_diverse_batch(pool, rank_value, int(opt["cases_per_iteration"]), config)
    print(
        "Surrogate proposal: "
        f"trained on {len(training_rows)} cases, "
        f"pool={pool_size}, "
        f"best predicted target error={float(np.nanmin(predicted_error)):.4g}"
    )
    return specs_from_pool(selected, iteration)


def propose_local_random_specs(config: dict, rows: list[dict[str, object]], iteration: int) -> list[CaseSpec]:
    opt = config["optimization"]
    bounds = config["parameter_bounds"]
    rng = random.Random(int(opt["random_seed"]) + iteration)
    best = best_rows(rows, limit=3)
    count = int(opt["cases_per_iteration"])
    local_fraction = float(opt["local_search_fraction"])

    if best:
        center = best[0]
    else:
        center = {
            "laser_power_W": np.mean(bounds["laser_power_W"]),
            "laser_radius_um": np.mean(bounds["laser_radius_um"]),
            "e_num_density": np.mean(bounds["e_num_density"]),
            "elec_resistivity_ohm_m": np.mean(bounds["elec_resistivity_ohm_m"]),
        }

    specs: list[CaseSpec] = []
    for index in range(count):
        use_local = index < max(1, int(round(count * 0.75)))
        values = {}
        for key, key_bounds in bounds.items():
            low, high = [float(value) for value in key_bounds]
            span = high - low
            if use_local:
                c = numeric(center, key, 0.5 * (low + high))
                value = rng.uniform(c - local_fraction * span, c + local_fraction * span)
                value = clip(value, key_bounds)
            else:
                value = rng.uniform(low, high)
            values[key] = value

        specs.append(
            CaseSpec(
                name=f"opt{index:02d}",
                laser_power_W=values["laser_power_W"],
                laser_radius_um=values["laser_radius_um"],
                e_num_density=values["e_num_density"],
                elec_resistivity_ohm_m=values["elec_resistivity_ohm_m"],
                stage="opt",
                iteration=iteration,
            )
        )
    return specs


def propose_optimization_specs(config: dict, rows: list[dict[str, object]], iteration: int) -> list[CaseSpec]:
    if str(config["optimization"].get("strategy", "ensemble_surrogate")) == "ensemble_surrogate":
        specs = propose_surrogate_specs(config, rows, iteration)
        if specs is not None:
            return specs
    return propose_local_random_specs(config, rows, iteration)


def run_specs(specs: list[CaseSpec], config: dict, force: bool, run: bool, jobs: int) -> list[Path]:
    runs_dir = WORK_DIR / str(config.get("runs_dir", "runs"))
    runs_dir.mkdir(exist_ok=True)
    write_planned_cases(specs, specs[0].stage if specs else "none")
    case_dirs = [copy_case(spec, config, force=force) for spec in specs]

    if not run:
        update_artifacts(config)
        return case_dirs

    with concurrent.futures.ThreadPoolExecutor(max_workers=jobs) as executor:
        futures = {executor.submit(run_case, case_dir): case_dir for case_dir in case_dirs}
        for future in concurrent.futures.as_completed(futures):
            case_dir = futures[future]
            try:
                future.result()
                print(f"Completed {case_dir.name}")
            except Exception as exc:
                print(f"FAILED {case_dir.name}: {exc}")
    return case_dirs


def collect_and_write_summary(config: dict) -> list[dict[str, object]]:
    runs_dir = WORK_DIR / str(config.get("runs_dir", "runs"))
    rows = [case_metrics(case_dir, config) for case_dir in sorted(runs_dir.iterdir()) if case_dir.is_dir()]
    write_summary(rows, WORK_DIR / str(config["summary_csv"]))
    return rows


def plot_sensitivity(rows: list[dict[str, object]], config: dict) -> None:
    plt = pyplot()
    if plt is None:
        return
    plots_dir = WORK_DIR / str(config.get("plots_dir", "plots"))
    plots_dir.mkdir(exist_ok=True)
    ok_rows = [row for row in rows if str(row.get("status")) == "ok"]
    if not ok_rows:
        return

    params = ["laser_power_W", "laser_radius_um", "e_num_density", "elec_resistivity_ohm_m"]
    fig, axes = plt.subplots(2, 2, figsize=(12, 8), constrained_layout=True)
    axes = axes.ravel()
    score = np.array([numeric(row, "score") for row in ok_rows])

    for ax, param in zip(axes, params):
        x = np.array([numeric(row, param) for row in ok_rows])
        if param == "e_num_density":
            x = x / 1e29
            xlabel = "e_num_density / 1e29"
        elif param == "elec_resistivity_ohm_m":
            x = x * 1e7
            xlabel = "elec_resistivity / 1e-7"
        else:
            xlabel = param
        ax.scatter(x, score, c=score, cmap="viridis_r", edgecolor="black")
        ax.set_xlabel(xlabel)
        ax.set_ylabel("score")
        ax.grid(True, alpha=0.3)

    fig.savefig(plots_dir / "sensitivity_parameter_scores.png", dpi=180)
    plt.close(fig)


def plot_best_aspect_ratio(rows: list[dict[str, object]], config: dict) -> None:
    plt = pyplot()
    if plt is None:
        return
    plots_dir = WORK_DIR / str(config.get("plots_dir", "plots"))
    plots_dir.mkdir(exist_ok=True)
    best = best_rows(rows, limit=1)
    if not best:
        return

    case_dir = Path(str(best[0]["case_dir"]))
    csv_path = case_dir / "post-processing-data" / "vtu_meltpool_geometry.csv"
    if not csv_path.exists():
        return

    data = pd.read_csv(csv_path)
    time_ms = pd.to_numeric(data["time"], errors="coerce").to_numpy(dtype=float) * 1.0e3
    width = pd.to_numeric(data["meltPoolWidth_um"], errors="coerce").to_numpy(dtype=float)
    depth = pd.to_numeric(data["meltPoolDepth_um"], errors="coerce").to_numpy(dtype=float)
    sim_ar = np.divide(depth, width, out=np.full_like(depth, np.nan), where=width > 0)
    stable_start_ms = numeric(best[0], "stable_start_ms")
    comparison_window_ms = float(config.get("target_stability", {}).get("comparison_window_ms", 2.0))
    if not np.isfinite(stable_start_ms):
        stable_start_ms = float(np.nanmin(time_ms[np.isfinite(time_ms)])) if np.any(np.isfinite(time_ms)) else 0.0
    sim_time_rel = time_ms - stable_start_ms
    window_mask = (sim_time_rel >= 0.0) & (sim_time_rel <= comparison_window_ms)

    exp_path = (WORK_DIR / str(config["exp_aspect_ratio_csv"])).resolve()
    fig, ax = plt.subplots(figsize=(8, 5), constrained_layout=True)
    finite = np.isfinite(sim_time_rel) & np.isfinite(sim_ar) & window_mask
    if np.any(finite):
        ax.plot(
            sim_time_rel[finite],
            sim_ar[finite],
            marker="o",
            linewidth=1.5,
            label=f"Simulation stable window: {case_dir.name}",
        )
    else:
        fallback = np.isfinite(time_ms) & np.isfinite(sim_ar)
        ax.plot(
            time_ms[fallback],
            sim_ar[fallback],
            marker="o",
            linewidth=1.5,
            label=f"Simulation full trace: {case_dir.name}",
        )

    if exp_path.exists():
        exp = pd.read_csv(exp_path)
        exp.columns = [column.strip() for column in exp.columns]
        exp_time = pd.to_numeric(exp["time_ms"], errors="coerce").to_numpy(dtype=float)
        exp_ar = pd.to_numeric(exp["meltpool_ar_depthByWidth"], errors="coerce").to_numpy(dtype=float)
        exp_finite = np.isfinite(exp_time) & np.isfinite(exp_ar)
        exp_time_rel = exp_time[exp_finite] - np.nanmin(exp_time[exp_finite])
        ax.plot(exp_time_rel, exp_ar[exp_finite], color="black", linewidth=2.0, label="Experiment, relative time")

    ax.set_xlim(0.0, comparison_window_ms)
    ax.set_xlabel("Time from stable comparison window start (ms)")
    ax.set_ylabel("Melt-pool aspect ratio, depth / width")
    ax.grid(True, alpha=0.3)
    ax.legend()
    fig.savefig(plots_dir / "best_aspect_ratio_vs_exp.png", dpi=180)
    plt.close(fig)


def status_counts(rows: list[dict[str, object]]) -> dict[str, int]:
    counts: dict[str, int] = {}
    for row in rows:
        status = str(row.get("status", "unknown"))
        stable = str(row.get("stable", "")).lower() == "true"
        key = "stable" if stable else status
        counts[key] = counts.get(key, 0) + 1
    return counts


def running_like_cases(rows: list[dict[str, object]]) -> list[dict[str, object]]:
    pending = []
    for row in rows:
        status = str(row.get("status", ""))
        final_time = numeric(row, "final_time_ms")
        if status in {"missing_csv", "empty_csv"} or (np.isfinite(final_time) and final_time < 1.25):
            pending.append(row)
    return pending


def plot_live_dashboard(rows: list[dict[str, object]], config: dict) -> None:
    plt = pyplot()
    if plt is None:
        return
    plots_dir = WORK_DIR / str(config.get("plots_dir", "plots"))
    plots_dir.mkdir(exist_ok=True)
    planned = read_planned_cases()
    planned_cases = list(planned.get("cases", []))
    best = best_rows(rows, limit=1)
    target = config["target_stability"]
    w_low, w_high = [float(value) for value in target["meltpool_width_um"]]
    d_low, d_high = [float(value) for value in target["meltpool_depth_um"]]

    fig = plt.figure(figsize=(18, 12), constrained_layout=True)
    gs = fig.add_gridspec(3, 4)
    ax_status = fig.add_subplot(gs[0, 0])
    ax_score = fig.add_subplot(gs[0, 1])
    ax_targets = fig.add_subplot(gs[0, 2])
    ax_best = fig.add_subplot(gs[0, 3])
    ax_progress = fig.add_subplot(gs[1, 0:2])
    ax_params = fig.add_subplot(gs[1, 2:4])
    ax_table = fig.add_subplot(gs[2, 0:2])
    ax_next = fig.add_subplot(gs[2, 2:4])

    fig.suptitle("Kube Keyhole Validation Progress", fontsize=18, fontweight="bold")

    counts = status_counts(rows)
    if counts:
        labels = list(counts)
        values = [counts[label] for label in labels]
        colors = ["#4C78A8" if label != "stable" else "#54A24B" for label in labels]
        ax_status.bar(labels, values, color=colors)
        ax_status.set_ylabel("cases")
        ax_status.set_title("Case Status")
        ax_status.tick_params(axis="x", rotation=30)
    else:
        ax_status.text(0.5, 0.5, "No cases yet", ha="center", va="center", transform=ax_status.transAxes)
        ax_status.set_axis_off()

    ok_rows = [row for row in rows if str(row.get("status")) == "ok"]
    if ok_rows:
        ordered = sorted(ok_rows, key=lambda row: str(row.get("case_id")))
        scores = np.array([numeric(row, "score") for row in ordered])
        best_so_far = np.minimum.accumulate(scores)
        ax_score.plot(np.arange(len(scores)), scores, "o-", label="case score", alpha=0.55)
        ax_score.plot(np.arange(len(scores)), best_so_far, "k-", linewidth=2.0, label="best so far")
        ax_score.set_title("Score Progress")
        ax_score.set_xlabel("completed case index")
        ax_score.set_ylabel("score")
        ax_score.grid(True, alpha=0.3)
        ax_score.legend(fontsize=8)

        width = np.array([numeric(row, "avg_stable_width_um") for row in ok_rows])
        depth = np.array([numeric(row, "avg_stable_depth_um") for row in ok_rows])
        score = np.array([numeric(row, "score") for row in ok_rows])
        scatter = ax_targets.scatter(width, depth, c=score, cmap="viridis_r", edgecolor="black")
        ax_targets.axvspan(w_low, w_high, color="#54A24B", alpha=0.12)
        ax_targets.axhspan(d_low, d_high, color="#54A24B", alpha=0.12)
        ax_targets.axvline(w_low, color="0.35", linestyle=":")
        ax_targets.axvline(w_high, color="0.35", linestyle=":")
        ax_targets.axhline(d_low, color="0.35", linestyle=":")
        ax_targets.axhline(d_high, color="0.35", linestyle=":")
        ax_targets.set_title("Stable-Window Dimensions")
        ax_targets.set_xlabel("width (um)")
        ax_targets.set_ylabel("depth (um)")
        ax_targets.grid(True, alpha=0.3)
        fig.colorbar(scatter, ax=ax_targets, label="score")
    else:
        ax_score.text(0.5, 0.5, "No completed cases", ha="center", va="center", transform=ax_score.transAxes)
        ax_targets.text(0.5, 0.5, "No dimension data", ha="center", va="center", transform=ax_targets.transAxes)
        ax_score.set_axis_off()
        ax_targets.set_axis_off()

    ax_best.set_title("Best Case")
    ax_best.set_axis_off()
    if best:
        row = best[0]
        lines = [
            f"case: {str(row.get('case_id'))[:42]}",
            f"score: {numeric(row, 'score'):.4g}",
            f"stable: {row.get('stable')}",
            f"stable duration: {numeric(row, 'stable_duration_ms'):.3g} ms",
            f"dimension gate: {row.get('dimension_gate_pass')}",
            f"avg width: {numeric(row, 'avg_stable_width_um'):.3g} um",
            f"avg depth: {numeric(row, 'avg_stable_depth_um'):.3g} um",
            f"depth/width: {numeric(row, 'avg_stable_depth_over_width'):.3g}",
            f"power: {numeric(row, 'laser_power_W'):.3g} W",
            f"radius: {numeric(row, 'laser_radius_um'):.3g} um",
            f"ne: {numeric(row, 'e_num_density') / 1e29:.3g}e29",
            f"er: {numeric(row, 'elec_resistivity_ohm_m'):.3g}",
        ]
        ax_best.text(0.0, 1.0, "\n".join(lines), va="top", ha="left", family="monospace")
    else:
        ax_best.text(0.5, 0.5, "No best case yet", ha="center", va="center", transform=ax_best.transAxes)

    pending = running_like_cases(rows)
    ax_progress.set_title("Current / Incomplete Cases")
    ax_progress.set_axis_off()
    if pending:
        lines = []
        for row in pending[:12]:
            lines.append(
                f"{str(row.get('case_id'))[:52]:52s}  "
                f"{str(row.get('status')):11s}  "
                f"t={numeric(row, 'final_time_ms'):.3g} ms"
            )
        ax_progress.text(0.0, 1.0, "\n".join(lines), va="top", ha="left", family="monospace", fontsize=8)
    else:
        ax_progress.text(0.5, 0.5, "No incomplete cases detected", ha="center", va="center", transform=ax_progress.transAxes)

    ax_params.set_title("Parameter Sensitivity Snapshot")
    if ok_rows:
        params = [
            ("laser_power_W", "Power W", 1.0),
            ("laser_radius_um", "Radius um", 1.0),
            ("e_num_density", "ne / 1e29", 1e-29),
            ("elec_resistivity_ohm_m", "er / 1e-7", 1e7),
        ]
        positions = np.arange(len(params))
        for row in ok_rows:
            values = [numeric(row, key) * scale for key, _, scale in params]
            score = numeric(row, "score")
            alpha = 0.25 if score > np.nanmedian([numeric(r, "score") for r in ok_rows]) else 0.65
            ax_params.plot(positions, values, color="#4C78A8", alpha=alpha, linewidth=1.0)
        if best:
            values = [numeric(best[0], key) * scale for key, _, scale in params]
            ax_params.plot(positions, values, color="#E45756", linewidth=2.5, marker="o", label="best")
            ax_params.legend(fontsize=8)
        ax_params.set_xticks(positions)
        ax_params.set_xticklabels([label for _, label, _ in params])
        ax_params.grid(True, alpha=0.3)
    else:
        ax_params.text(0.5, 0.5, "No sensitivity data yet", ha="center", va="center", transform=ax_params.transAxes)
        ax_params.set_axis_off()

    ax_table.set_title("Top Cases")
    ax_table.set_axis_off()
    top = best_rows(rows, limit=8)
    if top:
        lines = ["score    stable  dur(ms)  width  depth  case"]
        for row in top:
            lines.append(
                f"{numeric(row, 'score'):7.3g}  "
                f"{str(row.get('stable'))[:5]:5s}  "
                f"{numeric(row, 'stable_duration_ms'):7.3g}  "
                f"{numeric(row, 'avg_stable_width_um'):5.0f}  "
                f"{numeric(row, 'avg_stable_depth_um'):5.0f}  "
                f"{str(row.get('case_id'))[:44]}"
            )
        ax_table.text(0.0, 1.0, "\n".join(lines), va="top", ha="left", family="monospace", fontsize=8)
    else:
        ax_table.text(0.5, 0.5, "No ranked cases yet", ha="center", va="center", transform=ax_table.transAxes)

    ax_next.set_title(f"Next / Planned Cases ({planned.get('label', 'none')})")
    ax_next.set_axis_off()
    if planned_cases:
        lines = ["power  radius  ne(e29)  er(e-7)  case"]
        completed_ids = {str(row.get("case_id")) for row in rows}
        for item in planned_cases[:12]:
            marker = "done" if str(item.get("case_id")) in completed_ids else "next"
            lines.append(
                f"{float(item.get('laser_power_W', math.nan)):5.1f}  "
                f"{float(item.get('laser_radius_um', math.nan)):6.2f}  "
                f"{float(item.get('e_num_density', math.nan)) / 1e29:7.2f}  "
                f"{float(item.get('elec_resistivity_ohm_m', math.nan)) * 1e7:7.2f}  "
                f"{marker:4s}  {str(item.get('case_id'))[:32]}"
            )
        ax_next.text(0.0, 1.0, "\n".join(lines), va="top", ha="left", family="monospace", fontsize=8)
    else:
        ax_next.text(0.5, 0.5, "No planned cases recorded", ha="center", va="center", transform=ax_next.transAxes)

    fig.savefig(plots_dir / "live_progress_dashboard.png", dpi=170)
    plt.close(fig)


def save_best_case(rows: list[dict[str, object]], config: dict) -> None:
    best = best_rows(rows, limit=1)
    if not best:
        return
    path = WORK_DIR / "best_case_summary.json"
    path.write_text(json.dumps(best[0], indent=2) + "\n")


def any_stable(rows: list[dict[str, object]]) -> bool:
    return any(str(row.get("stable")).lower() == "true" for row in rows)


def update_artifacts(config: dict) -> list[dict[str, object]]:
    rows = collect_and_write_summary(config)
    plot_sensitivity(rows, config)
    plot_best_aspect_ratio(rows, config)
    plot_live_dashboard(rows, config)
    save_best_case(rows, config)
    return rows


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", type=Path, default=DEFAULT_CONFIG)
    parser.add_argument("--mode", choices=["sensitivity", "optimize", "all", "score"], default="sensitivity")
    parser.add_argument("--run", action="store_true")
    parser.add_argument("--prepare-only", action="store_true")
    parser.add_argument("--force", action="store_true")
    parser.add_argument("--jobs", type=int, default=None)
    parser.add_argument("--ranks-per-case", type=int, default=None)
    args = parser.parse_args()

    config = json.loads(args.config.read_text())
    config["_config_path"] = str(args.config.resolve())
    if args.ranks_per_case is not None:
        config["parallel"]["ranks_per_case"] = int(args.ranks_per_case)
    jobs = int(args.jobs or config["parallel"]["max_parallel_cases"])
    run = bool(args.run and not args.prepare_only)

    if args.mode in ["sensitivity", "all"]:
        specs = exploration_specs(config)
        run_specs(specs, config, force=args.force, run=run, jobs=jobs)
        rows = update_artifacts(config)
        if args.mode == "sensitivity":
            print(f"Wrote {WORK_DIR / str(config['summary_csv'])}")
            return
    else:
        rows = update_artifacts(config)

    if args.mode in ["optimize", "all"]:
        max_iterations = int(config["optimization"]["max_iterations"])
        for iteration in range(1, max_iterations + 1):
            rows = update_artifacts(config)
            if bool(config["optimization"]["stop_when_stable_case_found"]) and any_stable(rows):
                print("Stable target case found. Stopping optimization.")
                break
            if rows and not completed_training_rows(rows):
                print("No successful CFD cases available for optimization. Stopping.")
                break

            specs = propose_optimization_specs(config, rows, iteration)
            run_specs(specs, config, force=args.force, run=run, jobs=jobs)
            rows = update_artifacts(config)
            print(f"Finished optimization iteration {iteration}")

            if not run:
                break

        print(f"Wrote {WORK_DIR / str(config['summary_csv'])}")


if __name__ == "__main__":
    main()
