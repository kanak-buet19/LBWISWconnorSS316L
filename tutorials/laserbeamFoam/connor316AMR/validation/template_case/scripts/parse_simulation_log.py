#!/usr/bin/env python3
import sys
import re
import argparse
import csv
import json
import math
from bisect import bisect_left
from pathlib import Path
from tqdm import tqdm

def parse_control_dict(case_dir: Path) -> float:
    end_time = 1.08e-3
    dict_file = case_dir / "system" / "controlDict"
    if dict_file.exists():
        with open(dict_file, "r") as f:
            content = f.read()
        match = re.search(r"endTime\s+([\d\.e\-+]+);", content)
        if match:
            end_time = float(match.group(1))
    return end_time

def load_case_info(case_dir: Path) -> dict:
    case_info = case_dir / "case_info.json"
    if not case_info.exists():
        return {}
    try:
        with open(case_info, "r") as f:
            return json.load(f)
    except (OSError, json.JSONDecodeError):
        return {}


def resolve_validation_path(case_dir: Path, configured_path: str) -> Path:
    path = Path(configured_path)
    if path.is_absolute():
        return path
    for parent in (case_dir, *case_dir.parents):
        if (parent / "cases.json").is_file():
            return parent / path
    return case_dir / path


def optional_float(row: dict, key: str) -> float | None:
    try:
        value = float(row[key])
        return value if math.isfinite(value) and value > 0.0 else None
    except (KeyError, TypeError, ValueError):
        return None


def get_experimental_metrics(case_dir: Path, case_info: dict) -> tuple[float | None, float | None]:
    summary_files = sorted(case_dir.glob("exp_*_summary.csv"))
    for exp_csv in summary_files:
        try:
            with open(exp_csv, mode="r", newline="") as f:
                row = next(csv.DictReader(f), None)
            if row:
                width = optional_float(row, "width_um")
                depth = optional_float(row, "depth_um")
                if width is not None or depth is not None:
                    return width, depth
        except OSError:
            pass

    return (
        optional_float(case_info, "target_width_um"),
        optional_float(case_info, "target_depth_um"),
    )


def get_experimental_depth_timeseries(case_dir: Path, case_info: dict) -> list[tuple[float, float]]:
    configured_path = case_info.get("exp_timeseries_csv")
    if not configured_path:
        return []
    exp_csv = resolve_validation_path(case_dir, configured_path)
    try:
        with open(exp_csv, mode="r", newline="") as f:
            rows = list(csv.DictReader(f))
    except OSError:
        return []

    if rows and {"t_ms", "keyhole_depth_um"}.issubset(rows[0]):
        time_key, depth_key = "t_ms", "keyhole_depth_um"
    elif rows and {"time_ms", "depth_mean"}.issubset(rows[0]):
        time_key, depth_key = "time_ms", "depth_mean"
    elif rows and {"time_ms", "depth_um"}.issubset(rows[0]):
        time_key, depth_key = "time_ms", "depth_um"
    else:
        return []

    samples = []
    for row in rows:
        try:
            time_ms = float(row[time_key])
            depth_um = float(row[depth_key])
            if math.isfinite(time_ms) and time_ms >= 0.0 and depth_um > 0.0:
                samples.append((time_ms, depth_um))
        except (KeyError, TypeError, ValueError):
            pass
    return sorted(samples)


def interpolated_depth(samples: list[tuple[float, float]], time_ms: float) -> float | None:
    if not samples or time_ms < samples[0][0] or time_ms > samples[-1][0]:
        return None
    index = bisect_left(samples, (time_ms, -math.inf))
    if index == 0:
        return samples[0][1]
    if index == len(samples):
        return samples[-1][1]
    t0, depth0 = samples[index - 1]
    t1, depth1 = samples[index]
    if t1 == t0:
        return depth1
    fraction = (time_ms - t0) / (t1 - t0)
    return depth0 + fraction * (depth1 - depth0)

def get_latest_absorptivity(case_dir: Path) -> float:
    csv_file = case_dir / "absorptivity_vs_time" / "absorptivity_vs_time.csv"
    if not csv_file.exists():
        return 0.0
    try:
        with open(csv_file, "r") as f:
            lines = f.readlines()
            for last_line in reversed(lines):
                last_line = last_line.strip()
                if not last_line:
                    continue
                parts = last_line.split(",")
                if len(parts) >= 5:
                    try:
                        return float(parts[4])
                    except ValueError:
                        pass
    except Exception:
        pass
    return 0.0

def get_latest_meltpool_geometry(case_dir: Path, depth_metric: str) -> tuple[float, float, float]:
    csv_file = case_dir / "post-processing-data" / "vtk_meltpool_geometry.csv"
    if not csv_file.exists():
        return 0.0, 0.0, 0.0
    try:
        with open(csv_file, "r") as f:
            rows = list(csv.DictReader(f))
            depth_field = "keyholeDepth_um" if depth_metric == "keyhole" else "meltPoolDepth_um"
            for row in reversed(rows):
                try:
                    t_val = float(row["time"])
                    depth = float(row[depth_field])
                    width = float(row["meltPoolWidth_um"])
                    if not math.isnan(t_val) and (not math.isnan(depth) or not math.isnan(width)):
                        w_val = 0.0 if math.isnan(width) else width
                        d_val = 0.0 if math.isnan(depth) else depth
                        return t_val, w_val, d_val
                except (KeyError, TypeError, ValueError):
                    pass
    except Exception:
        pass
    return 0.0, 0.0, 0.0

def main():
    parser = argparse.ArgumentParser(description="Clean simulation log parser with tqdm progress bar")
    parser.add_argument("--append", action="store_true", help="Append to log file instead of overwriting")
    args = parser.parse_args()

    case_dir = Path(__file__).resolve().parents[1]
    end_time = parse_control_dict(case_dir)
    
    case_info = load_case_info(case_dir)
    exp_width, exp_depth = get_experimental_metrics(case_dir, case_info)
    exp_depth_timeseries = get_experimental_depth_timeseries(case_dir, case_info)
    depth_metric = case_info.get("depth_metric", "meltpool")
    
    log_file_path = case_dir / "log.laserbeamFoam"
    mode = "a" if args.append else "w"
    
    # Setup tqdm progress bar
    pbar = tqdm(
        total=end_time, 
        unit='s',
        bar_format='{desc} {percentage:3.0f}%|{bar}| {n:.3e}/{total:.3e} s [{elapsed}<{remaining}] {postfix}',
        dynamic_ncols=True
    )
    pbar.set_description("Simulation Progress")
    
    last_time = 0.0
    current_time = 0.0
    current_dt = 0.0
    max_co = 0.0
    power = 0.0
    laser_z = 0.0
    exec_time = 0.0
    abs_val = 0.0
    t_max = 0.0
    p_vap_max_kpa = 0.0
    max_u_metal = 0.0
    max_u_gas = 0.0
    last_printed_mp_time = -1.0
    first_time_seen = True
    dump_errors = False
    
    with open(log_file_path, mode, buffering=1) as log_file:
        for line in sys.stdin:
            # Save the raw output in the background log
            log_file.write(line)
            
            # If we are in error dump mode, write everything to stderr
            if dump_errors:
                sys.stderr.write(line)
                sys.stderr.flush()
                continue
            
            # Check for crash signatures to enter error dump mode
            lower_line = line.lower()
            if any(k in lower_line for k in ["fatal", "segmentation fault", "aborted", "sigsegv", "foam exiting"]):
                dump_errors = True
                pbar.close()
                sys.stderr.write(f"\n[ERROR DETECTED] Simulation crashed! Dumping solver output:\n")
                sys.stderr.write(line)
                sys.stderr.flush()
                continue
            
            # Parse metrics
            if line.startswith("Time = "):
                try:
                    current_time = float(line.split("Time = ")[1].strip())
                    if first_time_seen:
                        pbar.update(current_time)
                        last_time = current_time
                        first_time_seen = False
                    elif current_time > last_time:
                        pbar.update(current_time - last_time)
                        last_time = current_time
                        
                    # Fetch latest absorptivity at every solver step
                    abs_val = get_latest_absorptivity(case_dir)

                    # Check for new melt pool geometry predictions
                    mp_time, mp_width, mp_depth = get_latest_meltpool_geometry(case_dir, depth_metric)
                    if mp_time > 0.0 and mp_time > last_printed_mp_time:
                        comparisons = []
                        if exp_width is not None and mp_width > 0.0:
                            width_err = ((mp_width - exp_width) / exp_width) * 100.0
                            comparisons.append(f"Width Err: {width_err:+.1f}%")
                        depth_target = interpolated_depth(exp_depth_timeseries, mp_time * 1e3)
                        if depth_target is None:
                            depth_target = exp_depth
                        if depth_target is not None and mp_depth > 0.0:
                            depth_err = ((mp_depth - depth_target) / depth_target) * 100.0
                            comparisons.append(f"Depth Err: {depth_err:+.1f}%")
                        if comparisons:
                            pbar.write(
                                f"[ANALYSIS] Time: {mp_time*1e3:.3f} ms | "
                                + " | ".join(comparisons)
                            )
                        last_printed_mp_time = mp_time
                except (ValueError, IndexError):
                    pass
            elif "Courant Number mean:" in line:
                match = re.search(r"max:\s+([\d\.e\-+]+)", line)
                if match:
                    max_co = float(match.group(1))
            elif "deltaT = " in line:
                try:
                    current_dt = float(line.split("deltaT = ")[1].strip())
                except (ValueError, IndexError):
                    pass
            elif "power = " in line:
                try:
                    power = float(line.split("power = ")[1].strip())
                except (ValueError, IndexError):
                    pass
            elif "mean position = " in line:
                match = re.search(r"mean position = \(([\d\.e\-+\s]+)\)", line)
                if match:
                    parts = match.group(1).split()
                    if len(parts) >= 3:
                        laser_z = float(parts[2]) * 1000.0 # convert to mm
            elif "ExecutionTime = " in line:
                match = re.search(r"ExecutionTime = ([\d\.e\-+]+)", line)
                if match:
                    exec_time = float(match.group(1))
            elif "TMax = " in line:
                match = re.search(r"TMax = ([\d\.e\-+]+).*pVapMax = ([\d\.e\-+]+)", line)
                if match:
                    t_max = float(match.group(1))
                    p_vap_max_kpa = float(match.group(2)) / 1000.0
            elif "maxU_metal = " in line:
                match = re.search(r"maxU_metal = ([\d\.e\-+]+).*maxU_gas = ([\d\.e\-+]+)", line)
                if match:
                    max_u_metal = float(match.group(1))
                    max_u_gas = float(match.group(2))

            # Build and update progress bar postfix (keep it minimal and clean)
            postfix = {
                'dt': f"{current_dt:.2e}",
                'maxCo': f"{max_co:.2f}",
            }
            if abs_val > 0.0:
                postfix['Abs'] = f"{abs_val*100.0:.1f}%"
            if t_max > 0.0:
                postfix['TMax'] = f"{t_max:.0f}K"
            if p_vap_max_kpa > 0.0:
                postfix['pVap'] = f"{p_vap_max_kpa:.1f}kPa"
            if max_u_metal > 0.0 or max_u_gas > 0.0:
                postfix['Um'] = f"{max_u_metal:.3g}"
                postfix['Ug'] = f"{max_u_gas:.3g}"
            pbar.set_postfix(postfix)
            
    pbar.close()
    
    if dump_errors:
        print("\n[DEBUG] Simulation failed. Logs saved to log.laserbeamFoam.")
        sys.exit(1)
    else:
        print("\n[DEBUG] Simulation completed successfully. Logs saved to log.laserbeamFoam.")

if __name__ == "__main__":
    main()
