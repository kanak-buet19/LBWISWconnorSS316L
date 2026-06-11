#!/usr/bin/env python3
import os
import json
import re
import sys
from pathlib import Path

# Set non-interactive backend for matplotlib to avoid X11 errors
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

SMOOTHING_WINDOW = 101

def parse_unified_dat(dat_path):
    """Parses unified OpenFOAM fieldMinMax.dat file.
    """
    times = []
    max_Ts = []
    globalMaxPVaps = []
    
    # We will keep track of T and recoilPressure per time step to align them
    data_by_time = {}
    
    with open(dat_path, 'r', errors='ignore') as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith('#'):
                continue
            parts = line.split()
            if len(parts) >= 8:
                try:
                    t_val = float(parts[0])
                    field_name = parts[1]
                    max_val = float(parts[7])
                    
                    if t_val not in data_by_time:
                        data_by_time[t_val] = {'T': None, 'recoilPressure': None}
                    
                    if field_name == 'T':
                        data_by_time[t_val]['T'] = max_val
                    elif field_name == 'recoilPressure':
                        data_by_time[t_val]['recoilPressure'] = max_val
                except ValueError:
                    continue
                    
    # Sort by time and extract lists
    sorted_times = sorted(data_by_time.keys())
    for t in sorted_times:
        t_us = t * 1e6
        t_max = data_by_time[t]['T']
        rp_max = data_by_time[t]['recoilPressure']
        
        # We only append if both are present for that timestep
        if t_max is not None and rp_max is not None:
            times.append(t_us)
            max_Ts.append(t_max)
            globalMaxPVaps.append(rp_max / 1e3)  # Convert Pa to kPa
            
    return times, max_Ts, globalMaxPVaps

def centered_moving_average(values, window_size=SMOOTHING_WINDOW):
    """Return a centered moving average without shortening the series."""
    if not values:
        return []

    window_size = int(window_size)
    if window_size < 2 or len(values) < 3:
        return list(values)

    if window_size % 2 == 0:
        window_size += 1

    half_window = window_size // 2
    smoothed = []

    for i in range(len(values)):
        start = max(0, i - half_window)
        end = min(len(values), i + half_window + 1)
        window = values[start:end]
        smoothed.append(sum(window) / len(window))

    return smoothed

def centered_moving_std(values, window_size=SMOOTHING_WINDOW):
    """Return a centered moving standard deviation without shortening the series."""
    if not values:
        return []

    window_size = int(window_size)
    if window_size < 2 or len(values) < 3:
        return [0.0 for _ in values]

    if window_size % 2 == 0:
        window_size += 1

    half_window = window_size // 2
    moving_std = []

    for i in range(len(values)):
        start = max(0, i - half_window)
        end = min(len(values), i + half_window + 1)
        window = values[start:end]
        mean = sum(window) / len(window)
        variance = sum((value - mean)**2 for value in window) / len(window)
        moving_std.append(variance**0.5)

    return moving_std

def plus_minus(values, spread, multiplier=1.0):
    lower = [value - multiplier*delta for value, delta in zip(values, spread)]
    upper = [value + multiplier*delta for value, delta in zip(values, spread)]
    return lower, upper

def mean_value(values):
    if not values:
        return 0.0

    return sum(values) / len(values)

def min_value(values):
    if not values:
        return 0.0

    return min(values)

def max_value(values):
    if not values:
        return 0.0

    return max(values)

def main():
    case_dir = Path(__file__).resolve().parents[1]
    output_dir = case_dir / "post-processing"
    
    # Try finding unified fieldMinMax dat file first
    dat_files = sorted(case_dir.glob("postProcessing/fieldMinMax/*/fieldMinMax.dat"))
    
    parsed_via_dat = False
    times = []
    max_Ts = []
    globalMaxPVaps = []
    
    if dat_files:
        print(f"Found unified OpenFOAM fieldMinMax output file: {dat_files[-1]}")
        try:
            times, max_Ts, globalMaxPVaps = parse_unified_dat(dat_files[-1])
            if len(times) > 0:
                parsed_via_dat = True
                print(f"Extracted {len(times)} data points from native OpenFOAM function objects.")
        except Exception as e:
            print(f"Warning: Failed to parse .dat file ({e}), falling back to log parsing.")

    # Fallback to log file parsing
    if not parsed_via_dat:
        log_path = case_dir / "log.laserbeamFoam"
        if not log_path.exists():
            print(f"Error: Neither native data files nor log file found at {log_path}", file=sys.stderr)
            print("Please run the simulation first.", file=sys.stderr)
            sys.exit(1)
            
        print(f"Parsing solver log file: {log_path}")
        pattern = re.compile(
            r"\[DBG Recoil\]\s+"
            r"time=(?P<time>[0-9.eE+-]+)\s+"
            r"proc=\d+\s+"
            r"maxPVapCell=\d+\s+"
            r"C=\([^)]+\)\s+"
            r"pVap=(?P<pVap>[0-9.eE+-]+)\s+"
            r"globalMaxPVap=(?P<globalMaxPVap>[0-9.eE+-]+)\s+"
            r"T@pVap=(?P<Tvap>[0-9.eE+-]+)\s+"
            r"alpha@pVap=(?P<alpha>[0-9.eE+-]+)\s+"
            r"magGradAlpha@pVap=[0-9.eE+-]+\s+"
            r"maxT=(?P<maxT>[0-9.eE+-]+)"
        )
        
        times = []
        globalMaxPVaps = []
        max_Ts = []
        
        with open(log_path, 'r', errors='ignore') as f:
            for line in f:
                match = pattern.search(line)
                if match:
                    times.append(float(match.group('time')) * 1e6)  # Convert to microseconds
                    globalMaxPVaps.append(float(match.group('globalMaxPVap')) / 1e3)  # Convert to kPa
                    max_Ts.append(float(match.group('maxT')))
                    
        if not times:
            print("No data points found yet in the log file or postProcessing directories.")
            sys.exit(0)
            
        print(f"Successfully extracted {len(times)} data points from solver log.")

    # Ensure output directory exists
    output_dir.mkdir(parents=True, exist_ok=True)

    smooth_recoil = centered_moving_average(globalMaxPVaps)
    smooth_max_Ts = centered_moving_average(max_Ts)
    recoil_std = centered_moving_std(globalMaxPVaps)
    temp_std = centered_moving_std(max_Ts)
    recoil_lower, recoil_upper = plus_minus(smooth_recoil, recoil_std, multiplier=2.0)
    temp_lower, temp_upper = plus_minus(smooth_max_Ts, temp_std, multiplier=2.0)
    mean_recoil = mean_value(globalMaxPVaps)
    mean_max_T = mean_value(max_Ts)
    summary = {
        "smoothing_window_points": SMOOTHING_WINDOW,
        "std_dev_band_multiplier": 2.0,
        "data_points": len(times),
        "time_us": {
            "start": min_value(times),
            "end": max_value(times),
        },
        "recoil_pressure_kPa": {
            "min": min_value(globalMaxPVaps),
            "max": max_value(globalMaxPVaps),
            "mean": mean_recoil,
            "smoothed_min": min_value(smooth_recoil),
            "smoothed_max": max_value(smooth_recoil),
        },
        "peak_temperature_K": {
            "min": min_value(max_Ts),
            "max": max_value(max_Ts),
            "mean": mean_max_T,
            "smoothed_min": min_value(smooth_max_Ts),
            "smoothed_max": max_value(smooth_max_Ts),
        },
    }
    
    # ----------------------------------------------------
    # Combined plot: Recoil Pressure and Peak Temperature
    # ----------------------------------------------------
    fig, axes = plt.subplots(2, 1, figsize=(10, 8), sharex=True)
    fig.suptitle('Peak Recoil Pressure and Temperature History', fontsize=14, fontweight='bold')

    ax_press = axes[0]
    ax_press.plot(times, globalMaxPVaps, color='#E45756', linewidth=0.7, alpha=0.28, label='Raw Peak Recoil Pressure', zorder=1)
    ax_press.fill_between(times, recoil_lower, recoil_upper, color='#E45756', alpha=0.18, linewidth=0, label='Moving +/-2 std dev', zorder=2)
    ax_press.plot(times, smooth_recoil, color='black', linewidth=3.0, label=f'Smoothed Peak Recoil ({SMOOTHING_WINDOW}-point MA)', zorder=3)
    ax_press.axhline(y=mean_recoil, color='#E45756', linestyle='-.', linewidth=2.0, label=f'Mean ({mean_recoil:.1f} kPa)', zorder=2)
    ax_press.axhline(y=101.325, color='gray', linestyle='--', label='1 atm (101.3 kPa)')
    ax_press.set_title('Vaporization Recoil Pressure', fontsize=12, fontweight='bold')
    ax_press.set_ylabel('Pressure (kPa)', fontsize=10)
    ax_press.grid(True, linestyle=':', alpha=0.6)
    ax_press.legend(frameon=True, facecolor='white', edgecolor='black')

    ax_temp = axes[1]
    ax_temp.plot(times, max_Ts, color='#4C78A8', linewidth=0.7, alpha=0.28, label='Raw Max Domain Temp', zorder=1)
    ax_temp.fill_between(times, temp_lower, temp_upper, color='#4C78A8', alpha=0.18, linewidth=0, label='Moving +/-2 std dev', zorder=2)
    ax_temp.plot(times, smooth_max_Ts, color='black', linewidth=3.0, label=f'Smoothed Max Temp ({SMOOTHING_WINDOW}-point MA)', zorder=3)
    ax_temp.axhline(y=mean_max_T, color='#4C78A8', linestyle='-.', linewidth=2.0, label=f'Mean ({mean_max_T:.1f} K)', zorder=2)
    ax_temp.axhline(y=2792, color='red', linestyle=':', linewidth=2.0, label='Vaporization Temp (2792 K)', zorder=2)
    ax_temp.axhline(y=915, color='green', linestyle=':', linewidth=2.0, label='Liquidus (915 K)', zorder=2)
    ax_temp.set_title('Peak Temperature', fontsize=12, fontweight='bold')
    ax_temp.set_xlabel('Time (us)', fontsize=10)
    ax_temp.set_ylabel('Temperature (K)', fontsize=10)
    ax_temp.grid(True, linestyle=':', alpha=0.6)
    ax_temp.legend(frameon=True, facecolor='white', edgecolor='black')

    plt.tight_layout()

    combined_plot_path = output_dir / "recoil_temperature_history.png"
    plt.savefig(combined_plot_path, dpi=150)
    plt.close()
    print(f"Saved combined history plot to: {combined_plot_path}")

    json_path = output_dir / "recoil_temperature_history_summary.json"
    with open(json_path, "w") as f:
        json.dump(summary, f, indent=2)
    print(f"Saved combined history summary JSON to: {json_path}")
    
    print("\nAll plots updated inside the post-processing directory successfully!")

if __name__ == "__main__":
    main()
