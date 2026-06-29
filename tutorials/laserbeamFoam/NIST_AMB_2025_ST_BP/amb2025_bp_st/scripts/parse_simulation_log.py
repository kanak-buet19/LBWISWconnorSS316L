#!/usr/bin/env python3
import sys
import re
import argparse
import csv
import json
import math
import os
import time
from pathlib import Path
from tqdm import tqdm

EXP_WIDTH_DEFAULT = 114.602
EXP_DEPTH_DEFAULT = 93.953
USE_COLOR = sys.stderr.isatty() and not bool(os.environ.get("NO_COLOR"))

ANSI = {
    "reset": "\033[0m",
    "bold": "\033[1m",
    "dim": "\033[2m",
    "red": "\033[91m",
    "green": "\033[92m",
    "yellow": "\033[93m",
    "blue": "\033[94m",
    "magenta": "\033[95m",
    "cyan": "\033[96m",
    "white": "\033[97m",
}

def paint(text: str, *styles: str) -> str:
    if not USE_COLOR:
        return text
    return "".join(ANSI[item] for item in styles) + text + ANSI["reset"]

def metric(label: str, value: str, color: str = "white") -> str:
    return f"{paint(label, 'bold', color)} {paint(value, color)}"

def courant_color(value: float) -> str:
    if value >= 0.9:
        return "red"
    if value >= 0.5:
        return "yellow"
    return "green"

def temp_color(value: float) -> str:
    if value >= 3000.0:
        return "red"
    if value >= 2000.0:
        return "yellow"
    return "green"

def format_postfix(
    current_dt: float,
    max_co: float,
    power: float,
    laser_z: float,
    abs_val: float,
    t_max: float,
    p_vap_max_kpa: float,
    max_u_metal: float,
    max_u_gas: float,
) -> str:
    fields = [
        metric("dt", f"{current_dt:.2e}", "blue"),
        metric("Co", f"{max_co:.2f}", courant_color(max_co)),
    ]
    if abs_val > 0.0:
        fields.append(metric("Abs", f"{abs_val*100.0:.1f}%", "magenta"))
    if t_max > 0.0:
        fields.append(metric("Tmax", f"{t_max:.0f}K", temp_color(t_max)))
    if p_vap_max_kpa > 0.0:
        fields.append(metric("pVap", f"{p_vap_max_kpa:.1f}kPa", "yellow"))
    if max_u_metal > 0.0 or max_u_gas > 0.0:
        fields.append(metric("Um", f"{max_u_metal:.3g}", "cyan"))
        fields.append(metric("Ug", f"{max_u_gas:.3g}", "blue"))
    return "  ".join(fields)

def format_eta(seconds: float | None) -> str:
    if seconds is None or math.isinf(seconds) or seconds < 0:
        return "--"

    seconds_i = int(round(seconds))
    days, rem = divmod(seconds_i, 86400)
    hours, rem = divmod(rem, 3600)
    minutes, secs = divmod(rem, 60)

    if days:
        return f"{days}D {hours:02d}H {minutes:02d}M"
    if hours:
        return f"{hours}H {minutes:02d}M {secs:02d}S"
    if minutes:
        return f"{minutes}M {secs:02d}S"
    return f"{secs}S"

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

def get_sim_parameters(case_dir: Path) -> tuple[float, float, float, float]:
    p_laser = 200.0
    v_scan = 900.0
    d_laser = 0.05
    t_powder = 0.0

    # 1. Parse Laser Power from constant/timeVsLaserPower
    power_file = case_dir / "constant" / "timeVsLaserPower"
    if power_file.exists():
        with open(power_file, "r") as f:
            content = f.read()
        match = re.search(r"\(0\s+(\d+(?:\.\d+)?)\)", content)
        if match:
            p_laser = float(match.group(1))

    # 2. Parse Laser Radius from constant/LaserProperties
    laser_properties_file = case_dir / "constant" / "LaserProperties"
    if laser_properties_file.exists():
        with open(laser_properties_file, "r") as f:
            content = f.read()
        match_radius = re.search(r"laserRadius\s+(\d+(?:\.\d+)?(?:e-?\d+)?);", content)
        if match_radius:
            radius = float(match_radius.group(1))
            d_laser = round(2.0 * radius * 1000.0, 3) # in mm
        match_powder = re.search(r"PowderSim\s+(true|false);", content)
        if match_powder:
            if match_powder.group(1) == "false":
                t_powder = 0.0

    # 3. Parse Scan Speed from constant/timeVsLaserPosition
    pos_file = case_dir / "constant" / "timeVsLaserPosition"
    if pos_file.exists():
        with open(pos_file, "r") as f:
            content = f.read()
        matches = re.findall(r"\(([\d\.e\-]+)\s+\(([\d\.e\-]+)\s+([\d\.e\-]+)\s+([\d\.e\-]+)\)\)", content)
        if len(matches) >= 2:
            t0, x0, y0, z0 = map(float, matches[0])
            t1, x1, y1, z1 = map(float, matches[-1])
            dt = t1 - t0
            dz = abs(z1 - z0)
            if dt > 0:
                v_scan = round((dz * 1000.0) / dt, 1) # in mm/s

    return p_laser, v_scan, d_laser, t_powder

def get_experimental_metrics(case_dir: Path) -> tuple[float, float]:
    summary_files = sorted(case_dir.glob("exp_hofmann_*_summary.csv"))
    for exp_csv in summary_files:
        try:
            with open(exp_csv, mode="r", newline="") as f:
                row = next(csv.DictReader(f), None)
            if row:
                return float(row["width_um"]), float(row["depth_um"])
        except (KeyError, TypeError, ValueError, OSError):
            pass

    case_info = case_dir / "case_info.json"
    if case_info.exists():
        try:
            with open(case_info, "r") as f:
                info = json.load(f)
            return float(info["target_width_um"]), float(info["target_depth_um"])
        except (KeyError, TypeError, ValueError, OSError, json.JSONDecodeError):
            pass

    exp_csv = Path("/home/kanak/drives/d-drive/work/research/connor_project/papers/2026_Hofmann_Meltpool_data_316L/MeltpoolGeometryData.csv")
    if not exp_csv.exists():
        return EXP_WIDTH_DEFAULT, EXP_DEPTH_DEFAULT

    try:
        p_laser, v_scan, d_laser, t_powder = get_sim_parameters(case_dir)
        with open(exp_csv, mode='r') as f:
            reader = csv.reader(f)
            header = next(reader)
            rows = list(reader)

        matching_widths = []
        matching_depths = []
        tol = 1e-4
        for r in rows:
            if len(r) < 7:
                continue
            try:
                p_val = float(r[1])
                v_val = float(r[2])
                d_val = float(r[3])
                t_val = float(r[4])
                
                if (abs(p_val - p_laser) < tol and
                    abs(v_val - v_scan) < tol and
                    abs(d_val - d_laser) < tol and
                    abs(t_val - t_powder) < tol):
                    matching_widths.append(float(r[5]))
                    matching_depths.append(float(r[6]))
            except ValueError:
                continue

        if matching_widths:
            avg_width = sum(matching_widths) / len(matching_widths)
            avg_depth = sum(matching_depths) / len(matching_depths)
            return avg_width, avg_depth
    except Exception:
        pass
    return EXP_WIDTH_DEFAULT, EXP_DEPTH_DEFAULT

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

def get_latest_meltpool_geometry(case_dir: Path) -> tuple[float, float, float]:
    csv_file = case_dir / "post-processing-data" / "vtu_meltpool_geometry.csv"
    if not csv_file.exists():
        return 0.0, 0.0, 0.0
    try:
        with open(csv_file, "r") as f:
            lines = f.readlines()
            if len(lines) > 1:
                # Find the last line that has valid numbers
                for last_line in reversed(lines):
                    last_line = last_line.strip()
                    if not last_line:
                        continue
                    parts = last_line.split(",")
                    if len(parts) >= 14:
                        try:
                            t_val = float(parts[0])
                            depth = float(parts[7])
                            width = float(parts[13])
                            
                            if not math.isnan(t_val) and (not math.isnan(depth) or not math.isnan(width)):
                                w_val = 0.0 if math.isnan(width) else width
                                d_val = 0.0 if math.isnan(depth) else depth
                                return t_val, w_val, d_val
                        except ValueError:
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
    
    # Load experimental metrics
    exp_width, exp_depth = get_experimental_metrics(case_dir)
    
    log_file_path = case_dir / "log.laserbeamFoam"
    mode = "a" if args.append else "w"
    
    bar_format = (
        paint("{desc}", "bold", "cyan")
        + " {percentage:3.0f}%"
        + " {postfix}"
    )
    # Setup tqdm progress bar
    pbar = tqdm(
        total=end_time, 
        unit='s',
        bar_format=bar_format,
        dynamic_ncols=True,
    )
    pbar.set_description("AMB meltpool")
    
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
    start_wall_time = time.monotonic()
    
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
                    mp_time, mp_width, mp_depth = get_latest_meltpool_geometry(case_dir)
                    if mp_time > 0.0 and mp_time > last_printed_mp_time:
                        # Calculate relative errors compared to experimental values
                        width_err = ((mp_width - exp_width) / exp_width) * 100.0 if mp_width > 0.0 else float('nan')
                        depth_err = ((mp_depth - exp_depth) / exp_depth) * 100.0 if mp_depth > 0.0 else float('nan')
                        
                        # Format errors beautifully
                        width_err_str = f"{width_err:+.1f}%" if not math.isnan(width_err) else "N/A"
                        depth_err_str = f"{depth_err:+.1f}%" if not math.isnan(depth_err) else "N/A"
                        
                        pbar.write(
                            paint("[MELTPOOL]", "bold", "magenta")
                            + f" t={mp_time*1e3:.3f} ms  "
                            + metric("width_err", width_err_str, "yellow")
                            + "  "
                            + metric("depth_err", depth_err_str, "yellow")
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

            elapsed_wall_time = time.monotonic() - start_wall_time
            if pbar.n > 0.0 and pbar.total:
                remaining_wall_time = elapsed_wall_time * max(pbar.total - pbar.n, 0.0) / pbar.n
            else:
                remaining_wall_time = None

            postfix = (
                paint("ETA", "bold", "yellow")
                + " "
                + paint(format_eta(remaining_wall_time), "yellow")
            )
            stats = format_postfix(
                current_dt,
                max_co,
                power,
                laser_z,
                abs_val,
                t_max,
                p_vap_max_kpa,
                max_u_metal,
                max_u_gas,
            )
            if stats:
                postfix += "  " + stats
            pbar.set_postfix_str(postfix)
            
    pbar.close()
    
    if dump_errors:
        print("\n[DEBUG] Simulation failed. Logs saved to log.laserbeamFoam.")
        sys.exit(1)
    else:
        print("\n[DEBUG] Simulation completed successfully. Logs saved to log.laserbeamFoam.")

if __name__ == "__main__":
    main()
