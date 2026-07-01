#!/usr/bin/env python3
"""Build and run a small CW parametric sweep from the calibration template."""

from __future__ import annotations

import argparse
import csv
import json
import os
import shutil
import subprocess
import sys
import time
from pathlib import Path
from typing import Any


ROOT = Path(__file__).resolve().parents[1]
REPO = ROOT.parent
CONFIG_PATH = ROOT / "parametric_config.json"

CALIB_SCRIPTS = REPO / "calibrate_CW_case" / "scripts"
sys.path.insert(0, str(CALIB_SCRIPTS))

import caselib  # noqa: E402
import evaluate as ev  # noqa: E402


def load_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text())


def baseline_params(calib_config: dict[str, Any]) -> dict[str, float]:
    params: dict[str, float] = {}
    for name, spec in calib_config["parameters"].items():
        if name.startswith("_"):
            continue
        params[name] = float(spec["baseline"])
    return params


def command_exists(command: str) -> bool:
    return shutil.which(command) is not None


def case_surface_y(case_dir: Path) -> float | None:
    try:
        return float(json.loads((case_dir / "case_build.json").read_text())["surface_y_um"])
    except (OSError, KeyError, TypeError, ValueError):
        return None


def sim_progress(case_dir: Path) -> tuple[float, float | None] | None:
    log_file = case_dir / "log.laserbeamFoam"
    if not log_file.exists():
        return None
    try:
        tail = log_file.read_bytes()[-40000:].decode("utf-8", "ignore")
    except OSError:
        return None

    sim_time = None
    delta_t = None
    for line in tail.splitlines():
        if line.startswith("Time = "):
            try:
                sim_time = float(line[7:].strip())
            except ValueError:
                pass
        elif line.startswith("deltaT = "):
            try:
                delta_t = float(line[9:].strip())
            except ValueError:
                pass
    return (sim_time, delta_t) if sim_time is not None else None


def progress_line(case_index: int, total: int, case_id: str, case_dir: Path, end_time: float, status: str) -> str:
    progress = sim_progress(case_dir)
    width = 28
    if progress is None or end_time <= 0:
        bar = "." * ((int(time.monotonic()) % width) + 1)
        bar = bar.ljust(width, "-")
        return f"\r[{bar}] {case_index}/{total} {case_id} {status}"

    sim_time, delta_t = progress
    frac = min(max(sim_time / end_time, 0.0), 1.0)
    filled = int(width * frac)
    bar = "#" * filled + "-" * (width - filled)
    dt_text = "n/a" if delta_t is None else f"{delta_t:.2e}"
    return (
        f"\r[{bar}] {case_index}/{total} {case_id} "
        f"t={sim_time:.3e}/{end_time:.3e}s ({100.0 * frac:5.1f}%) dt={dt_text}"
    )


def read_manifest(results_root: Path) -> list[dict[str, Any]]:
    manifest = results_root / "case_manifest.json"
    if not manifest.exists():
        raise SystemExit("case_manifest.json not found. Run with --generate-only first or run ./Allrun.")
    return json.loads(manifest.read_text())


def generate_cases(config: dict[str, Any], calib_config: dict[str, Any]) -> list[dict[str, Any]]:
    template = (ROOT / config["template_case"]).resolve()
    case_root = ROOT / config["case_root"]
    results_root = ROOT / config["results_root"]
    case_root.mkdir(parents=True, exist_ok=True)
    results_root.mkdir(parents=True, exist_ok=True)

    case_cfg = calib_config["cases"][int(config.get("case_index", 0))]
    geom = dict(calib_config["geometry"])
    geom.update(config.get("geometry_overrides", {}))
    control = calib_config["control"]
    cores = int(os.environ.get("PARAM_CORES_PER_SIM", "4"))
    max_refinement = int(config.get("max_refinement", 2))
    base_params = baseline_params(calib_config)

    active_ids = {item["case_id"] for item in config["sweep_cases"]}
    for old_case in case_root.iterdir():
        if old_case.is_dir() and old_case.name not in active_ids:
            shutil.rmtree(old_case)

    rows: list[dict[str, Any]] = []
    for item in config["sweep_cases"]:
        case_id = item["case_id"]
        params = dict(base_params)
        params[item["parameter"]] = float(item["value"])
        dest = case_root / case_id
        build = caselib.build_case(
            template,
            dest,
            params,
            case_cfg,
            geom,
            control,
            cores,
            max_refinement=max_refinement,
        )
        row = {
            "case_id": case_id,
            "label": item["label"],
            "parameter": item["parameter"],
            "baseline_value": base_params[item["parameter"]],
            "value": float(item["value"]),
            "case_dir": str(dest.relative_to(ROOT)),
            "target_width_um": case_cfg["exp_width_um"],
            "target_depth_um": case_cfg["exp_depth_um"],
            "power_W": case_cfg["P_laser_W"],
            "scan_speed_mm_s": case_cfg["v_scan_mm_s"],
            "laser_radius_um": case_cfg["laserRadius_m"] * 1e6,
            "end_time_s": build["geometry"]["end_time"],
            "max_refinement": max_refinement,
            "cores": cores,
        }
        rows.append(row)

    (results_root / "case_manifest.json").write_text(json.dumps(rows, indent=2) + "\n")
    return rows


def evaluate_case(row: dict[str, Any], calib_config: dict[str, Any], status: str, return_code: int | None) -> dict[str, Any]:
    case_dir = ROOT / row["case_dir"]
    objective = calib_config.get("objective", {})
    result = dict(row)
    result.update({"status": status, "return_code": return_code})
    metrics = ev.evaluate(
        case_dir,
        float(row["target_depth_um"]),
        float(row["target_width_um"]),
        depth_weight=float(objective.get("depth_weight", 1.0)),
        width_weight=float(objective.get("width_weight", 1.0)),
        aspect_ratio_weight=float(objective.get("aspect_ratio_weight", 1.0)),
    )
    for key in (
        "sim_depth_um",
        "sim_width_um",
        "depth_err",
        "width_err",
        "ar_err",
        "case_error",
        "metric_source",
        "n_points",
        "converged",
    ):
        result[key] = metrics.get(key)
    return result


def write_summary(results_root: Path, rows: list[dict[str, Any]]) -> None:
    results_root.mkdir(parents=True, exist_ok=True)
    (results_root / "summary.json").write_text(json.dumps(rows, indent=2) + "\n")
    fields = [
        "case_id",
        "label",
        "parameter",
        "baseline_value",
        "value",
        "status",
        "return_code",
        "sim_width_um",
        "sim_depth_um",
        "target_width_um",
        "target_depth_um",
        "width_err",
        "depth_err",
        "ar_err",
        "case_error",
        "metric_source",
        "n_points",
        "converged",
        "case_dir",
    ]
    with (results_root / "summary.csv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def run_case(
    case_dir: Path,
    log_path: Path,
    env: dict[str, str],
    verbose: bool,
    case_index: int,
    total_cases: int,
    case_id: str,
    end_time: float,
) -> int:
    run_env = dict(env)
    surface_y = case_surface_y(case_dir)
    if surface_y is not None:
        run_env["SURFACE_Y_UM"] = f"{surface_y:.8g}"
    run_env.setdefault("FOAM_SIGFPE", "0")
    run_env.setdefault("MPLCONFIGDIR", str(ROOT / ".mplconfig"))
    run_env.setdefault("DELETE_ANALYZED_VTK", "true")
    log_path.parent.mkdir(parents=True, exist_ok=True)
    if verbose:
        return subprocess.run(["./Allrun_long"], cwd=case_dir, env=run_env, check=False).returncode

    with log_path.open("w") as handle:
        process = subprocess.Popen(
            ["./Allrun_long"],
            cwd=case_dir,
            env=run_env,
            stdout=handle,
            stderr=subprocess.STDOUT,
        )

        while process.poll() is None:
            sys.stdout.write(progress_line(case_index, total_cases, case_id, case_dir, end_time, "starting"))
            sys.stdout.flush()
            time.sleep(5)

    sys.stdout.write(progress_line(case_index, total_cases, case_id, case_dir, end_time, "done"))
    sys.stdout.write("\n")
    sys.stdout.flush()
    return int(process.returncode or 0)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--generate-only", action="store_true", help="Build cases and manifest, then stop.")
    parser.add_argument("--skip-existing", action="store_true", help="Skip cases with final result CSV output.")
    parser.add_argument("--only", nargs="*", help="Optional case_id list to run.")
    parser.add_argument("--verbose", action="store_true", help="Stream case output to terminal.")
    args = parser.parse_args()

    config = load_json(CONFIG_PATH)
    calib_config = load_json((ROOT / config["baseline_config"]).resolve())
    results_root = ROOT / config["results_root"]

    manifest_path = results_root / "case_manifest.json"
    if args.skip_existing and manifest_path.exists():
        manifest = read_manifest(results_root)
    else:
        manifest = generate_cases(config, calib_config)
    if args.generate_only:
        write_summary(results_root, [dict(row, status="generated", return_code=None) for row in manifest])
        print(f"Generated {len(manifest)} cases under {ROOT / config['case_root']}")
        return

    if not command_exists("laserbeamFoam") or not os.environ.get("WM_PROJECT_DIR"):
        raise SystemExit("OpenFOAM/laserbeamFoam is not available. Source OpenFOAM first, or use --generate-only.")

    only = set(args.only or [])
    env = os.environ.copy()
    env.setdefault("PYTHON", os.environ.get("PYTHON", sys.executable))
    rows: list[dict[str, Any]] = []
    run_logs = results_root / "run_logs"
    selected_rows = [row for row in read_manifest(results_root) if not only or row["case_id"] in only]
    total_selected = len(selected_rows)
    selected_index = 0

    for row in read_manifest(results_root):
        case_id = row["case_id"]
        case_dir = ROOT / row["case_dir"]
        final_csv = case_dir / ev.FINAL_CSV_REL
        series_csv = case_dir / ev.CSV_REL

        if only and case_id not in only:
            rows.append(dict(row, status="not_selected", return_code=None))
            continue
        selected_index += 1
        if args.skip_existing and (final_csv.exists() or series_csv.exists()):
            rows.append(evaluate_case(row, calib_config, "skipped_existing", 0))
            write_summary(results_root, rows)
            print(f"[skip] {selected_index}/{total_selected} {case_id}")
            continue

        print(f"[run] {selected_index}/{total_selected} {case_id}: {row['label']}", flush=True)
        return_code = run_case(
            case_dir,
            run_logs / f"{case_id}.log",
            env,
            args.verbose,
            selected_index,
            total_selected,
            case_id,
            float(row["end_time_s"]),
        )
        status = "complete" if return_code == 0 and (final_csv.exists() or series_csv.exists()) else "failed"
        rows.append(evaluate_case(row, calib_config, status, return_code))
        write_summary(results_root, rows)
        print(f"[{status}] {case_id}", flush=True)

    write_summary(results_root, rows)
    print(f"Wrote {results_root / 'summary.csv'}")


if __name__ == "__main__":
    main()
