#!/usr/bin/env python3
import argparse
import os
import sys
from pathlib import Path

# Set non-interactive backend for matplotlib to avoid X11 errors
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

def parse_unified_dat(dat_path):
    times = []
    max_Ts = []
    globalMaxPVaps = []
    
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
                    
    sorted_times = sorted(data_by_time.keys())
    for t in sorted_times:
        t_us = t * 1e6
        t_max = data_by_time[t]['T']
        rp_max = data_by_time[t]['recoilPressure']
        
        if t_max is not None and rp_max is not None:
            times.append(t_us)
            max_Ts.append(t_max)
            globalMaxPVaps.append(rp_max / 1e3)  # Convert Pa to kPa
            
    return times, max_Ts, globalMaxPVaps

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--case-dir", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--case-name", type=str, required=True)
    args = parser.parse_args()
    
    case_dir = args.case_dir.resolve()
    output_dir = args.output_dir.resolve()
    case_name = args.case_name
    
    dat_files = sorted(case_dir.glob("postProcessing/fieldMinMax/*/fieldMinMax.dat"))
    if not dat_files:
        print(f"Error: No fieldMinMax.dat found under {case_dir}", file=sys.stderr)
        sys.exit(1)
        
    times, max_Ts, globalMaxPVaps = parse_unified_dat(dat_files[-1])
    
    if not times:
        print(f"Error: No synchronized data points extracted from {dat_files[-1]}", file=sys.stderr)
        sys.exit(1)
        
    output_dir.mkdir(parents=True, exist_ok=True)
    
    # Create beautiful, stacked 2-subplot figure
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(9, 7), sharex=True)
    
    # Subplot 1: Recoil Pressure vs Time
    ax1.plot(times, globalMaxPVaps, color='#E45756', linewidth=2, marker='o', markersize=3.5, label='Peak Recoil Pressure')
    ax1.axhline(y=101.325, color='gray', linestyle='--', alpha=0.8, label='1 atm (101.3 kPa)')
    ax1.set_ylabel('Recoil Pressure (kPa)', fontsize=10, fontweight='bold')
    ax1.set_title(f'Keyhole Solver Histories: {case_name}', fontsize=12, fontweight='bold', pad=12)
    ax1.grid(True, linestyle=':', alpha=0.6)
    ax1.legend(frameon=True, facecolor='white', edgecolor='black', loc='upper left')
    
    # Subplot 2: Max Temperature vs Time
    ax2.plot(times, max_Ts, color='#4C78A8', linewidth=2, marker='s', markersize=3.5, label='Max Domain Temp')
    ax2.axhline(y=2792, color='red', linestyle=':', linewidth=1.5, label='Vaporization Temp (2792 K)')
    ax2.axhline(y=915, color='green', linestyle=':', linewidth=1.5, label='Liquidus (915 K)')
    ax2.set_xlabel('Time (us)', fontsize=10, fontweight='bold')
    ax2.set_ylabel('Temperature (K)', fontsize=10, fontweight='bold')
    ax2.grid(True, linestyle=':', alpha=0.6)
    ax2.legend(frameon=True, facecolor='white', edgecolor='black', loc='upper left')
    
    plt.tight_layout()
    
    # Save to central plots/recoil_pressure_plots/
    central_plot_path = output_dir / f"{case_name}.png"
    plt.savefig(central_plot_path, dpi=150)
    plt.close()
    
    # Save copy locally inside case post-processing folder
    local_post_dir = case_dir / "post-processing"
    local_post_dir.mkdir(parents=True, exist_ok=True)
    shutil_copy = True
    try:
        import shutil
        shutil.copy2(central_plot_path, local_post_dir / "recoil_temperature_subplot.png")
    except Exception:
        pass
        
    print(f"Successfully generated recoil pressure and temperature subplot: {central_plot_path}")

if __name__ == "__main__":
    main()
