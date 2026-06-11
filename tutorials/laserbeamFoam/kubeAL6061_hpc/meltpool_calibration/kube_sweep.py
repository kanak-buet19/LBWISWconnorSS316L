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

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


WORK_DIR = Path(__file__).resolve().parent
BASE_CASE = WORK_DIR.parent
RUNS_DIR = WORK_DIR / "runs"
PLOTS_DIR = WORK_DIR / "plots"
SUMMARY_CSV = WORK_DIR / "summary.csv"
DEFAULT_CONFIG = WORK_DIR / "calibration_config.json"
NEXT_PROPOSAL_JSON = WORK_DIR / "next_proposal.json"
ADAPTIVE_STATE_JSON = WORK_DIR / "adaptive_state.json"


@dataclass(frozen=True)
class CaseSpec:
    radius_um: float
    resistivity_ohm_m: float
    mesh_size_um: float

    @property
    def radius_m(self) -> float:
        return self.radius_um * 1e-6

    @property
    def case_id(self) -> str:
        er_micro = self.resistivity_ohm_m * 1e6
        return f"m{self.mesh_size_um:05.1f}um_r{self.radius_um:06.2f}um_er{er_micro:06.3f}e-6"


def value_list(value) -> list[float]:
    if isinstance(value, list):
        return [float(item) for item in value]
    if isinstance(value, dict):
        start = float(value["start"])
        stop = float(value["stop"])
        step = float(value["step"])
        count = int(round((stop - start) / step)) + 1
        return [start + index * step for index in range(count)]
    return [float(value)]


def replace_scalar(path: Path, name: str, value: float | int) -> None:
    text = path.read_text()
    pattern = re.compile(rf"(^\s*{re.escape(name)}\s+)([^;]+)(;.*$)", re.MULTILINE)
    text, count = pattern.subn(rf"\g<1>{value:.12g}\g<3>", text, count=1)
    if count != 1:
        raise ValueError(f"Could not patch {name} in {path}")
    path.write_text(text)


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
    enabled = bool(config.get("dynamic_mesh", True))
    mesh_type = "dynamicRefineFvMesh" if enabled else "staticFvMesh"
    text, count = re.subn(
        r"(^\s*dynamicFvMesh\s+)([^;]+)(;.*$)",
        rf"\g<1>{mesh_type}\g<3>",
        text,
        count=1,
        flags=re.MULTILINE,
    )
    if count != 1:
        raise ValueError(f"Could not patch dynamicFvMesh in {dynamic_mesh}")

    text = re.sub(
        r"(^\s*maxRefinement\s+)([^;]+)(;.*$)",
        rf"\g<1>{int(config.get('max_refinement', 1))}\g<3>",
        text,
        count=1,
        flags=re.MULTILINE,
    )
    text = re.sub(
        r"(^\s*maxCells\s+)([^;]+)(;.*$)",
        rf"\g<1>{int(config.get('max_cells', 1000000))}\g<3>",
        text,
        count=1,
        flags=re.MULTILINE,
    )
    dynamic_mesh.write_text(text)


def copy_case(spec: CaseSpec, config: dict, force: bool) -> Path:
    case_dir = RUNS_DIR / spec.case_id
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
        "output_*.out",
        "output_*.err",
        "*.foam",
        "__pycache__",
        "meltpool_calibration",
        "connor316AMR",
        "kubeAL6061",
    )
    shutil.copytree(BASE_CASE, case_dir, ignore=ignore)

    replace_scalar(case_dir / "constant" / "LaserProperties", "laserRadius", spec.radius_m)
    replace_scalar(case_dir / "constant" / "transportProperties", "elec_resistivity", spec.resistivity_ohm_m)
    replace_scalar(case_dir / "system" / "decomposeParDict", "numberOfSubdomains", int(config.get("ranks_per_case", 8)))
    patch_block_mesh_size(case_dir, spec.mesh_size_um)
    patch_dynamic_mesh(case_dir, config)

    if "simulation_end_s" in config:
        replace_scalar(case_dir / "system" / "controlDict", "endTime", float(config["simulation_end_s"]))
    if "write_interval_s" in config:
        replace_scalar(case_dir / "system" / "controlDict", "writeInterval", float(config["write_interval_s"]))

    return case_dir


def specs_from_config(config: dict) -> list[CaseSpec]:
    if "exploration_cases" in config:
        return [
            CaseSpec(
                float(case["laser_radius_um"]),
                float(case["elec_resistivity_ohm_m"]),
                float(case.get("mesh_size_um", config.get("mesh_size_um", 20.0))),
            )
            for case in config["exploration_cases"]
        ]

    radii = value_list(config["laser_radius_um"])
    resistivities = value_list(config["elec_resistivity_ohm_m"])
    mesh_sizes = value_list(config.get("mesh_size_um", 20.0))
    return [CaseSpec(radius, er, mesh_size) for mesh_size in mesh_sizes for radius in radii for er in resistivities]


def run_case(case_dir: Path, run_command: str) -> tuple[Path, str]:
    subprocess.run(["bash", "-lc", run_command], cwd=case_dir, check=True, stdin=subprocess.DEVNULL)
    return case_dir, "ok"


def last_finite(series: pd.Series) -> float:
    values = pd.to_numeric(series, errors="coerce").dropna()
    if values.empty:
        return math.nan
    return float(values.iloc[-1])


def peak_finite(series: pd.Series) -> float:
    values = pd.to_numeric(series, errors="coerce").dropna()
    if values.empty:
        return math.nan
    return float(values.max())


def keyhole_start_ms(data: pd.DataFrame, threshold_um: float) -> float:
    if "keyholeDepth_um" not in data:
        return math.nan

    depth = pd.to_numeric(data["keyholeDepth_um"], errors="coerce")
    time_s = pd.to_numeric(data["time"], errors="coerce")
    active = data[(depth > threshold_um) & time_s.notna()]
    if active.empty:
        return math.nan

    return float(pd.to_numeric(active["time"], errors="coerce").iloc[0]) * 1e3


def case_metrics(case_dir: Path, spec: CaseSpec, config: dict) -> dict[str, object]:
    csv_path = case_dir / "post-processing-data" / "vtu_meltpool_geometry.csv"
    stop_reason_path = case_dir / "STOP_REASON.json"
    stop_reason = json.loads(stop_reason_path.read_text()) if stop_reason_path.exists() else {}
    if not csv_path.exists():
        return {
            "case_id": spec.case_id,
            "mesh_size_um": spec.mesh_size_um,
            "laser_radius_um": spec.radius_um,
            "elec_resistivity_ohm_m": spec.resistivity_ohm_m,
            "status": "missing_csv",
            "stop_decision": stop_reason.get("decision", ""),
            "stop_reason": stop_reason.get("reason", ""),
            "case_dir": str(case_dir),
        }

    data = pd.read_csv(csv_path)
    if data.empty:
        status = "empty_csv"
    else:
        status = "ok"

    beam_diameter_um = 2.0 * spec.radius_um if bool(config.get("keyhole_ar_uses_case_beam_diameter", True)) else 80.0
    keyhole_depth = pd.to_numeric(data.get("keyholeDepth_um"), errors="coerce")
    melt_depth = pd.to_numeric(data.get("meltPoolDepth_um"), errors="coerce")
    melt_width = pd.to_numeric(data.get("meltPoolWidth_um"), errors="coerce")
    onset_ms = keyhole_start_ms(data, float(config.get("keyhole_onset_threshold_um", 1.0)))
    target_onset_ms = float(config.get("target_keyhole_start_ms", 4.2709))
    sim_shift_ms = float(config.get("sim_time_shift_ms", 0.0))

    return {
        "case_id": spec.case_id,
        "mesh_size_um": spec.mesh_size_um,
        "laser_radius_um": spec.radius_um,
        "elec_resistivity_ohm_m": spec.resistivity_ohm_m,
        "beam_diameter_um": beam_diameter_um,
        "status": status,
        "stop_decision": stop_reason.get("decision", ""),
        "stop_reason": stop_reason.get("reason", ""),
        "final_time_ms": last_finite(data["time"]) * 1e3 if "time" in data else math.nan,
        "keyhole_start_ms": onset_ms,
        "keyhole_start_exp_time_ms": onset_ms + sim_shift_ms if math.isfinite(onset_ms) else math.nan,
        "keyhole_start_error_ms": onset_ms - target_onset_ms if math.isfinite(onset_ms) else math.nan,
        "keyhole_start_abs_error_ms": abs(onset_ms - target_onset_ms) if math.isfinite(onset_ms) else math.nan,
        "peak_keyhole_depth_um": peak_finite(keyhole_depth),
        "peak_keyhole_ar": peak_finite(keyhole_depth / beam_diameter_um),
        "peak_meltpool_depth_um": peak_finite(melt_depth),
        "peak_meltpool_ar": peak_finite(melt_width / melt_depth.replace(0.0, np.nan)),
        "case_dir": str(case_dir),
    }


def write_summary(rows: list[dict[str, object]]) -> None:
    fields = [
        "case_id",
        "mesh_size_um",
        "laser_radius_um",
        "elec_resistivity_ohm_m",
        "beam_diameter_um",
        "status",
        "stop_decision",
        "stop_reason",
        "final_time_ms",
        "keyhole_start_ms",
        "keyhole_start_exp_time_ms",
        "keyhole_start_error_ms",
        "keyhole_start_abs_error_ms",
        "peak_keyhole_depth_um",
        "peak_keyhole_ar",
        "peak_meltpool_depth_um",
        "peak_meltpool_ar",
        "case_dir",
    ]
    with SUMMARY_CSV.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def surrogate_terms(radius_um: np.ndarray, resistivity_ohm_m: np.ndarray) -> np.ndarray:
    radius_center = 41.0
    radius_span = 3.0
    er_center = 3.5e-7
    er_span = 1.0e-7
    r = (radius_um - radius_center) / radius_span
    er = (resistivity_ohm_m - er_center) / er_span
    return np.column_stack([np.ones_like(r), r, er, r * r, r * er, er * er])


def write_next_proposal(rows: list[dict[str, object]], config: dict) -> None:
    target_ms = float(config.get("target_keyhole_start_ms", 4.2709))
    valid_rows = [
        row for row in rows
        if row.get("status") == "ok" and math.isfinite(float(row.get("keyhole_start_ms", math.nan)))
    ]

    payload: dict[str, object] = {
        "target_keyhole_start_ms": target_ms,
        "status": "insufficient_completed_cases",
        "message": "Run the exploration cases first; at least 6 completed onset measurements are needed for the quadratic surrogate.",
    }

    if len(valid_rows) >= 6:
        radius = np.array([float(row["laser_radius_um"]) for row in valid_rows])
        er = np.array([float(row["elec_resistivity_ohm_m"]) for row in valid_rows])
        onset = np.array([float(row["keyhole_start_ms"]) for row in valid_rows])

        coeffs, *_ = np.linalg.lstsq(surrogate_terms(radius, er), onset, rcond=None)

        proposal_cfg = config.get("surrogate_proposal", {})
        radius_min = float(proposal_cfg.get("laser_radius_min_um", 39.0))
        radius_max = float(proposal_cfg.get("laser_radius_max_um", 42.0))
        er_min = float(proposal_cfg.get("elec_resistivity_min_ohm_m", 3.0e-7))
        er_max = float(proposal_cfg.get("elec_resistivity_max_ohm_m", 4.0e-7))
        radius_count = int(proposal_cfg.get("radius_grid_count", 61))
        er_count = int(proposal_cfg.get("elec_resistivity_grid_count", 61))

        radius_grid, er_grid = np.meshgrid(
            np.linspace(radius_min, radius_max, radius_count),
            np.linspace(er_min, er_max, er_count),
            indexing="ij",
        )
        predicted = surrogate_terms(radius_grid.ravel(), er_grid.ravel()).dot(coeffs)
        best_index = int(np.nanargmin(np.abs(predicted - target_ms)))

        payload = {
            "target_keyhole_start_ms": target_ms,
            "status": "ok",
            "model": "quadratic least-squares surrogate for keyhole_start_ms",
            "best_completed_case": min(
                valid_rows,
                key=lambda row: float(row.get("keyhole_start_abs_error_ms", math.inf)),
            ),
            "proposed_next_case": {
                "laser_radius_um": float(radius_grid.ravel()[best_index]),
                "elec_resistivity_ohm_m": float(er_grid.ravel()[best_index]),
                "mesh_size_um": float(config.get("mesh_size_um", 20.0)),
                "predicted_keyhole_start_ms": float(predicted[best_index]),
                "predicted_abs_error_ms": float(abs(predicted[best_index] - target_ms)),
            },
            "surrogate_coefficients": [float(value) for value in coeffs],
        }

    NEXT_PROPOSAL_JSON.write_text(json.dumps(payload, indent=2) + "\n")


def case_from_payload(payload: dict[str, object]) -> CaseSpec:
    return CaseSpec(
        float(payload["laser_radius_um"]),
        float(payload["elec_resistivity_ohm_m"]),
        float(payload.get("mesh_size_um", 20.0)),
    )


def spec_key(spec: CaseSpec) -> tuple[int, int, int]:
    return (
        int(round(spec.radius_um * 1000.0)),
        int(round(spec.resistivity_ohm_m * 1.0e12)),
        int(round(spec.mesh_size_um * 1000.0)),
    )


def load_state(config: dict) -> dict[str, object]:
    if ADAPTIVE_STATE_JSON.exists():
        state = json.loads(ADAPTIVE_STATE_JSON.read_text())
        if state.get("pending_cases"):
            return state

        if not SUMMARY_CSV.exists() and not NEXT_PROPOSAL_JSON.exists():
            state["pending_cases"] = [
                {
                    "laser_radius_um": spec.radius_um,
                    "elec_resistivity_ohm_m": spec.resistivity_ohm_m,
                    "mesh_size_um": spec.mesh_size_um,
                }
                for spec in specs_from_config(config)
            ]
            return state

        if NEXT_PROPOSAL_JSON.exists():
            proposal = json.loads(NEXT_PROPOSAL_JSON.read_text())
            if proposal.get("status") != "ok":
                state["pending_cases"] = [
                    {
                        "laser_radius_um": spec.radius_um,
                        "elec_resistivity_ohm_m": spec.resistivity_ohm_m,
                        "mesh_size_um": spec.mesh_size_um,
                    }
                    for spec in specs_from_config(config)
                ]
                return state

        return state

    return {
        "iteration": 0,
        "pending_cases": [
            {
                "laser_radius_um": spec.radius_um,
                "elec_resistivity_ohm_m": spec.resistivity_ohm_m,
                "mesh_size_um": spec.mesh_size_um,
            }
            for spec in specs_from_config(config)
        ],
        "completed_case_ids": [],
        "history": [],
    }


def save_state(state: dict[str, object]) -> None:
    ADAPTIVE_STATE_JSON.write_text(json.dumps(state, indent=2) + "\n")


def completed_specs_from_summary() -> set[tuple[int, int, int]]:
    if not SUMMARY_CSV.exists():
        return set()

    data = pd.read_csv(SUMMARY_CSV)
    completed = set()
    for _, row in data.iterrows():
        if row.get("status") != "ok":
            continue
        completed.add(
            spec_key(
                CaseSpec(
                    float(row["laser_radius_um"]),
                    float(row["elec_resistivity_ohm_m"]),
                    float(row["mesh_size_um"]),
                )
            )
        )
    return completed


def append_proposal_cases(state: dict[str, object], config: dict, count: int) -> None:
    if not NEXT_PROPOSAL_JSON.exists():
        return

    proposal = json.loads(NEXT_PROPOSAL_JSON.read_text())
    if proposal.get("status") != "ok":
        return

    proposed = proposal.get("proposed_next_case", {})
    if not proposed:
        return

    proposal_cfg = config.get("surrogate_proposal", {})
    radius_step = float(proposal_cfg.get("adaptive_radius_step_um", 0.25))
    er_step = float(proposal_cfg.get("adaptive_elec_resistivity_step_ohm_m", 0.1e-7))
    radius_min = float(proposal_cfg.get("laser_radius_min_um", 39.0))
    radius_max = float(proposal_cfg.get("laser_radius_max_um", 42.0))
    er_min = float(proposal_cfg.get("elec_resistivity_min_ohm_m", 3.0e-7))
    er_max = float(proposal_cfg.get("elec_resistivity_max_ohm_m", 4.0e-7))
    mesh_size = float(proposed.get("mesh_size_um", config.get("mesh_size_um", 20.0)))
    center_radius = float(proposed["laser_radius_um"])
    center_er = float(proposed["elec_resistivity_ohm_m"])

    completed = completed_specs_from_summary()
    pending_payloads = list(state.get("pending_cases", []))
    pending = {spec_key(case_from_payload(item)) for item in pending_payloads}
    candidates: list[CaseSpec] = []
    offsets = [
        (0.0, 0.0),
        (-radius_step, 0.0),
        (radius_step, 0.0),
        (0.0, -er_step),
        (0.0, er_step),
        (-radius_step, -er_step),
        (-radius_step, er_step),
        (radius_step, -er_step),
        (radius_step, er_step),
    ]

    for dr, der in offsets:
        spec = CaseSpec(
            min(max(center_radius + dr, radius_min), radius_max),
            min(max(center_er + der, er_min), er_max),
            mesh_size,
        )
        key = spec_key(spec)
        if key in completed or key in pending:
            continue
        candidates.append(spec)
        pending.add(key)
        if len(candidates) >= count:
            break

    pending_payloads.extend(
        {
            "laser_radius_um": spec.radius_um,
            "elec_resistivity_ohm_m": spec.resistivity_ohm_m,
            "mesh_size_um": spec.mesh_size_um,
        }
        for spec in candidates
    )
    state["pending_cases"] = pending_payloads
    state.setdefault("history", []).append(
        {
            "iteration": state.get("iteration", 0),
            "proposal_source": proposed,
            "added_cases": [
                {
                    "laser_radius_um": spec.radius_um,
                    "elec_resistivity_ohm_m": spec.resistivity_ohm_m,
                    "mesh_size_um": spec.mesh_size_um,
                }
                for spec in candidates
            ],
        }
    )


def has_matched_case(rows: list[dict[str, object]]) -> bool:
    return any(row.get("stop_decision") == "stop_matched" for row in rows)


def run_specs(specs: list[CaseSpec], config: dict, run_command: str, jobs: int, force: bool) -> list[dict[str, object]]:
    RUNS_DIR.mkdir(exist_ok=True)
    case_dirs = [copy_case(spec, config, force) for spec in specs]

    with concurrent.futures.ThreadPoolExecutor(max_workers=jobs) as executor:
        futures = {
            executor.submit(run_case, case_dir, run_command): case_dir
            for case_dir in case_dirs
        }
        for future in concurrent.futures.as_completed(futures):
            case_dir = futures[future]
            try:
                future.result()
                print(f"Completed {case_dir.name}")
            except Exception as exc:
                print(f"FAILED {case_dir.name}: {exc}")

    return [case_metrics(case_dir, spec, config) for case_dir, spec in zip(case_dirs, specs)]


def merge_summary_rows(new_rows: list[dict[str, object]]) -> list[dict[str, object]]:
    if SUMMARY_CSV.exists():
        existing = pd.read_csv(SUMMARY_CSV).to_dict(orient="records")
    else:
        existing = []

    by_case_id = {str(row["case_id"]): row for row in existing}
    for row in new_rows:
        by_case_id[str(row["case_id"])] = row
    return list(by_case_id.values())


def run_adaptive(config: dict, run_command: str, jobs: int, force: bool) -> None:
    adaptive = config.get("adaptive", {})
    batch_size = int(adaptive.get("batch_size", jobs))
    max_iterations = int(adaptive.get("max_iterations", 20))
    stop_on_match = bool(adaptive.get("stop_on_match", True))
    state = load_state(config)

    for _ in range(max_iterations):
        pending_payloads = list(state.get("pending_cases", []))
        if not pending_payloads:
            append_proposal_cases(state, config, batch_size)
            pending_payloads = list(state.get("pending_cases", []))
            if not pending_payloads:
                print("No pending cases and no valid proposal available.")
                break

        batch_payloads = pending_payloads[:batch_size]
        state["pending_cases"] = pending_payloads[batch_size:]
        state["iteration"] = int(state.get("iteration", 0)) + 1
        save_state(state)

        specs = [case_from_payload(item) for item in batch_payloads]
        print(f"=== Adaptive iteration {state['iteration']}: running {len(specs)} cases ===")
        rows = run_specs(specs, config, run_command, jobs, force)
        all_rows = merge_summary_rows(rows)
        write_summary(all_rows)
        write_next_proposal(all_rows, config)
        plot_comparison(all_rows, [CaseSpec(float(row["laser_radius_um"]), float(row["elec_resistivity_ohm_m"]), float(row["mesh_size_um"])) for row in all_rows], config)

        state.setdefault("completed_case_ids", []).extend(str(row["case_id"]) for row in rows)
        save_state(state)

        if stop_on_match and has_matched_case(all_rows):
            print("Matched keyhole onset found; stopping adaptive calibration.")
            break

        append_proposal_cases(state, config, batch_size)
        save_state(state)

    print(f"Wrote {SUMMARY_CSV}")
    print(f"Wrote {NEXT_PROPOSAL_JSON}")
    print(f"Wrote {ADAPTIVE_STATE_JSON}")


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


def sim_series(case_dir: Path, spec: CaseSpec, metric: str, config: dict) -> tuple[np.ndarray, np.ndarray]:
    data = pd.read_csv(case_dir / "post-processing-data" / "vtu_meltpool_geometry.csv")
    time_ms = pd.to_numeric(data["time"], errors="coerce").to_numpy(dtype=float) * 1e3
    time_ms = time_ms + float(config.get("sim_time_shift_ms", 0.0))
    beam_diameter_um = 2.0 * spec.radius_um if bool(config.get("keyhole_ar_uses_case_beam_diameter", True)) else 80.0

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


def plot_comparison(rows: list[dict[str, object]], specs: list[CaseSpec], config: dict) -> None:
    PLOTS_DIR.mkdir(exist_ok=True)
    exp_dir = BASE_CASE / "kube_exp_data" / "csv"
    exp_shift_ms = float(config.get("exp_time_shift_ms", 0.0))
    sim_shift_ms = float(config.get("sim_time_shift_ms", 0.0))
    onset_ms = exp_onset_ms(exp_dir) if bool(config.get("plot_from_keyhole_onset", True)) else -math.inf
    x_min = onset_ms - exp_shift_ms if math.isfinite(onset_ms) else None

    spec_by_id = {spec.case_id: spec for spec in specs}
    metrics = [
        ("keyhole_depth", "keyhole_depth_exp.csv", "depth_um", "Keyhole depth", "Depth (um)"),
        ("keyhole_ar", "keyhole_ar.csv", "aspect_ratio_keyholeDepthDividedbyLaserbeamdia80um", "Keyhole aspect ratio", "Depth / beam diameter"),
        ("meltpool_depth", "meltpool_depth_exp.csv", "depth_um", "Melt-pool depth", "Depth (um)"),
        ("meltpool_ar", "meltpool_ar_exp.csv", "aspect_ratio_meltpoolWidthbyDepth", "Melt-pool aspect ratio", "Width / depth"),
    ]

    fig, axes = plt.subplots(2, 2, figsize=(14, 9), constrained_layout=True)
    axes = axes.ravel()

    for ax, (metric, exp_file, exp_col, title, ylabel) in zip(axes, metrics):
        exp_t, exp_y = load_exp(exp_dir, exp_file, exp_col, onset_ms, exp_shift_ms)
        exp_label = "Experiment" if exp_shift_ms == 0.0 else f"Experiment shifted -{exp_shift_ms:.3f} ms"
        ax.plot(exp_t, exp_y, color="black", linewidth=2.3, label=exp_label)

        for row in rows:
            if row.get("status") != "ok":
                continue
            spec = spec_by_id[str(row["case_id"])]
            case_dir = Path(str(row["case_dir"]))
            if not (case_dir / "post-processing-data" / "vtu_meltpool_geometry.csv").exists():
                continue
            sim_t, sim_y = sim_series(case_dir, spec, metric, config)
            if sim_t.size:
                x_min = float(np.nanmin(sim_t)) if x_min is None else min(x_min, float(np.nanmin(sim_t)))
            label = f"m={spec.mesh_size_um:g} um, r={spec.radius_um:g} um, er={spec.resistivity_ohm_m:.2e}, sim +{sim_shift_ms:.3f} ms"
            ax.plot(sim_t, sim_y, marker="o", markersize=2.5, linewidth=1.2, label=label)

        if x_min is not None:
            ax.set_xlim(left=x_min)
        ax.set_title(title)
        ax.set_xlabel("Time (ms)")
        ax.set_ylabel(ylabel)
        ax.grid(True, alpha=0.3)

    handles, labels = axes[0].get_legend_handles_labels()
    fig.legend(handles, labels, loc="lower center", bbox_to_anchor=(0.5, -0.08), ncol=3, fontsize=8)
    fig.savefig(PLOTS_DIR / "sim_vs_exp_transient.png", dpi=180, bbox_inches="tight")
    plt.close(fig)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", type=Path, default=DEFAULT_CONFIG)
    parser.add_argument("--run-command", required=True)
    parser.add_argument("--jobs", type=int, default=None)
    parser.add_argument("--ranks-per-case", type=int, default=None)
    parser.add_argument("--force", action="store_true")
    parser.add_argument("--adaptive", action="store_true")
    args = parser.parse_args()

    config = json.loads(args.config.read_text())
    if args.ranks_per_case is not None:
        config["ranks_per_case"] = args.ranks_per_case
    jobs = int(args.jobs or config.get("max_parallel_cases", 8))

    if args.adaptive or bool(config.get("adaptive", {}).get("enabled", False)):
        run_adaptive(config, args.run_command, jobs, args.force)
        return

    specs = specs_from_config(config)
    rows = run_specs(specs, config, args.run_command, jobs, args.force)
    write_summary(rows)
    write_next_proposal(rows, config)
    plot_comparison(rows, specs, config)
    print(f"Wrote {SUMMARY_CSV}")
    print(f"Wrote {NEXT_PROPOSAL_JSON}")
    print(f"Wrote {PLOTS_DIR / 'sim_vs_exp_transient.png'}")


if __name__ == "__main__":
    main()
