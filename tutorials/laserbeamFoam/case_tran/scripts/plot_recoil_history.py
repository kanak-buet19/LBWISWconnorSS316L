#!/usr/bin/env python3
import os
import re
import sys
from pathlib import Path

# Set non-interactive backend for matplotlib to avoid X11 errors
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

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
    
    # ----------------------------------------------------
    # Plot 1: Recoil Pressure vs Time
    # ----------------------------------------------------
    plt.figure(figsize=(8, 5))
    plt.plot(times, globalMaxPVaps, color='#E45756', linewidth=2, marker='o', markersize=3, label='Peak Recoil Pressure')
    plt.axhline(y=101.325, color='gray', linestyle='--', label='1 atm (101.3 kPa)')
    plt.title('Vaporization Recoil Pressure History', fontsize=12, fontweight='bold')
    plt.xlabel('Time (us)', fontsize=10)
    plt.ylabel('Pressure (kPa)', fontsize=10)
    plt.grid(True, linestyle=':', alpha=0.6)
    plt.legend(frameon=True, facecolor='white', edgecolor='black')
    plt.tight_layout()
    
    press_plot_path = output_dir / "recoil_pressure_history.png"
    plt.savefig(press_plot_path, dpi=150)
    plt.close()
    print(f"Saved pressure history plot to: {press_plot_path}")
    
    # ----------------------------------------------------
    # Plot 2: Peak Temperatures vs Time
    # ----------------------------------------------------
    plt.figure(figsize=(8, 5))
    plt.plot(times, max_Ts, color='#4C78A8', linewidth=2, marker='s', markersize=3, label='Max Domain Temp (maxT)')
    plt.axhline(y=2500, color='red', linestyle=':', label='Vaporization Temp (2500 K)')
    plt.axhline(y=915, color='green', linestyle=':', label='Liquidus (915 K)')
    plt.title('Peak Temperature History', fontsize=12, fontweight='bold')
    plt.xlabel('Time (us)', fontsize=10)
    plt.ylabel('Temperature (K)', fontsize=10)
    plt.grid(True, linestyle=':', alpha=0.6)
    plt.legend(frameon=True, facecolor='white', edgecolor='black')
    plt.tight_layout()
    
    temp_plot_path = output_dir / "temperature_history.png"
    plt.savefig(temp_plot_path, dpi=150)
    plt.close()
    print(f"Saved temperature history plot to: {temp_plot_path}")
    
    print("\nAll plots updated inside the post-processing directory successfully!")

if __name__ == "__main__":
    main()
