#!/usr/bin/env python3
"""Run generated sensitivity cases and collect melt-pool width metrics."""

from __future__ import annotations

import argparse
import csv
import json
import os
import shutil
import sys
import subprocess
import time
from pathlib import Path
from typing import Any


ROOT = Path(__file__).resolve().parents[1]
CONFIG_PATH = ROOT / "sensitivity_config.json"


def load_config() -> dict[str, Any]:
    return json.loads(CONFIG_PATH.read_text())


def command_exists(command: str) -> bool:
    return shutil.which(command) is not None


def latest_geometry_row(case_dir: Path) -> dict[str, str] | None:
    csv_path = case_dir / "post-processing-data" / "vtu_meltpool_geometry.csv"
    if not csv_path.exists():
        return None
    with csv_path.open(newline="") as handle:
        rows = list(csv.DictReader(handle))
    if not rows:
        return None
    return max(rows, key=lambda row: float(row.get("time", "0") or 0.0))


def read_manifest(config: dict[str, Any]) -> list[dict[str, Any]]:
    manifest_path = ROOT / config["results_root"] / "case_manifest.json"
    if not manifest_path.exists():
        raise SystemExit("case_manifest.json not found. Run scripts/generate_cases.py first.")
    return json.loads(manifest_path.read_text())


def case_summary(row: dict[str, Any], case_dir: Path, status: str, return_code: int | None) -> dict[str, Any]:
    geometry = latest_geometry_row(case_dir)
    summary = dict(row)
    summary.update(
        {
            "status": status,
            "return_code": return_code,
            "time_s": None,
            "meltPoolWidth_um": None,
            "meltPoolDepth_um": None,
            "keyholeWidth_um": None,
            "keyholeDepth_um": None,
        }
    )
    if geometry is not None:
        for key in ("time", "meltPoolWidth_um", "meltPoolDepth_um", "keyholeWidth_um", "keyholeDepth_um"):
            out_key = "time_s" if key == "time" else key
            value = geometry.get(key, "")
            try:
                summary[out_key] = float(value)
            except (TypeError, ValueError):
                summary[out_key] = None
    return summary


def write_summary(config: dict[str, Any], rows: list[dict[str, Any]]) -> None:
    results_root = ROOT / config["results_root"]
    results_root.mkdir(parents=True, exist_ok=True)
    json_path = results_root / "sensitivity_results.json"
    json_path.write_text(json.dumps(rows, indent=2))

    csv_path = results_root / "sensitivity_results.csv"
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
        "return_code",
        "time_s",
        "meltPoolWidth_um",
        "meltPoolDepth_um",
        "keyholeWidth_um",
        "keyholeDepth_um",
        "case_dir",
    ]
    with csv_path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def run_case(case_dir: Path, env: dict[str, str], log_path: Path, verbose: bool) -> int:
    if verbose:
        process = subprocess.run(["./Allrun"], cwd=case_dir, env=env, check=False)
        return process.returncode

    log_path.parent.mkdir(parents=True, exist_ok=True)
    with log_path.open("w") as log_handle:
        process = subprocess.run(
            ["./Allrun"],
            cwd=case_dir,
            env=env,
            stdout=log_handle,
            stderr=subprocess.STDOUT,
            check=False,
        )
    return process.returncode


def progress_bar(done: int, total: int, label: str, status: str = "") -> None:
    width = 28
    filled = 0 if total == 0 else int(width * done / total)
    bar = "#" * filled + "-" * (width - filled)
    suffix = f" {status}" if status else ""
    line = f"\r[{bar}] {done:>2}/{total:<2} {label}{suffix}"
    sys.stdout.write(line[:120].ljust(120))
    sys.stdout.flush()


def finish_progress() -> None:
    sys.stdout.write("\n")
    sys.stdout.flush()


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--skip-existing", action="store_true", help="Do not rerun cases with geometry CSV output.")
    parser.add_argument("--dry-run", action="store_true", help="Only write planned-case summary, do not call OpenFOAM.")
    parser.add_argument("--only", nargs="*", help="Optional case_id list to run.")
    parser.add_argument("--quiet", action="store_true", help="Keep terminal output minimal.")
    parser.add_argument("--verbose", action="store_true", help="Stream full case output to terminal.")
    args = parser.parse_args()

    config = load_config()
    manifest = read_manifest(config)
    only = set(args.only or [])

    if not args.dry_run and (not command_exists("laserbeamFoam") or not os.environ.get("WM_PROJECT_DIR")):
        raise SystemExit("OpenFOAM/laserbeamFoam is not available. Source OpenFOAM first, or use --dry-run.")

    rows: list[dict[str, Any]] = []
    env = os.environ.copy()
    env.setdefault("PYTHON", sys.executable or config.get("python", "python3"))
    env.setdefault("MPLCONFIGDIR", "/tmp")
    run_logs = ROOT / config["results_root"] / "run_logs"
    selected_total = sum(1 for row in manifest if not only or row["case_id"] in only)
    selected_index = 0
    started_at = time.monotonic()

    for row in manifest:
        case_id = row["case_id"]
        if only and case_id not in only:
            rows.append(case_summary(row, ROOT / row["case_dir"], "not_selected", None))
            continue

        selected_index += 1
        case_dir = ROOT / row["case_dir"]
        if args.dry_run:
            if args.quiet:
                progress_bar(selected_index, selected_total, case_id, "planned")
            rows.append(case_summary(row, case_dir, "planned", None))
            continue

        if args.skip_existing and latest_geometry_row(case_dir) is not None:
            if args.quiet:
                progress_bar(selected_index, selected_total, case_id, "skipped")
            else:
                print(f"[skip] {case_id}: existing geometry CSV")
            rows.append(case_summary(row, case_dir, "skipped_existing", 0))
            continue

        if args.quiet:
            progress_bar(selected_index - 1, selected_total, case_id, "running")
        else:
            print(f"[run] {case_id}")

        log_path = run_logs / f"{case_id}.log"
        return_code = run_case(case_dir, env, log_path, args.verbose)
        status = "complete" if return_code == 0 and latest_geometry_row(case_dir) is not None else "failed"
        if args.quiet:
            progress_bar(selected_index, selected_total, case_id, "ok" if status == "complete" else "failed")
        if status == "failed":
            if args.quiet:
                finish_progress()
            print(f"FAILED: {case_id} (see {log_path})")
        rows.append(case_summary(row, case_dir, status, return_code))
        write_summary(config, rows)

    if args.quiet:
        elapsed = time.monotonic() - started_at
        progress_bar(selected_total, selected_total, "done", f"{elapsed:.0f}s")
        finish_progress()
    write_summary(config, rows)
    if not args.quiet:
        print(f"Wrote {ROOT / config['results_root'] / 'sensitivity_results.csv'}")


if __name__ == "__main__":
    main()
