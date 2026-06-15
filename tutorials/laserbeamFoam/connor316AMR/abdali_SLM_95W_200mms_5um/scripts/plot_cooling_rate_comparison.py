#!/usr/bin/env python3
import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path

# Paths
CASE_DIR = Path(__file__).resolve().parents[1]
CSV_PATH = CASE_DIR / "post-processing-data" / "cooling_rate_comparison.csv"
OUT_DIR = CASE_DIR / "post-processing-data"
OUT_PLOT_PNG = OUT_DIR / "cooling_rate_comparison_plot.png"
OUT_PLOT_PDF = OUT_DIR / "cooling_rate_comparison_plot.pdf"

# Paper data from Abdali et al. (2024)
paper_distances = [20, 80, 140, 200]  # um
paper_cooling_rates = [972400, 384300, 319500, 286400]  # K/s

def main():
    if not CSV_PATH.exists():
        print(f"Error: Simulation CSV file not found at: {CSV_PATH}")
        print("Please run the simulation first to generate the cooling rate data.")
        return

    # Load simulation CSV
    try:
        df = pd.read_csv(CSV_PATH)
    except Exception as e:
        print(f"Error reading CSV: {e}")
        return

    if df.empty:
        print("Error: The CSV file is empty.")
        return

    # Filter for steady state (t >= 2.5 ms)
    t_start_steady = 2.5e-3
    steady_df = df[df["Time_s"] >= t_start_steady]
    
    if steady_df.empty:
        print(f"Warning: No data found for steady-state (t >= {t_start_steady*1000} ms). Using all available data instead.")
        steady_df = df

    # Get average values over the steady-state window for robustness
    avg_row = steady_df.mean()

    # Extract simulation cooling rates (steady-state averages)
    sim_crs_avg = [
        avg_row["CRS_20um"],
        avg_row["CRS_80um"],
        avg_row["CRS_140um"],
        avg_row["CRS_200um"]
    ]
    sim_crt_avg = [
        avg_row["CRT_20um"],
        avg_row["CRT_80um"],
        avg_row["CRT_140um"],
        avg_row["CRT_200um"]
    ]

    # Create a premium 3-panel figure
    fig, axes = plt.subplots(1, 3, figsize=(18, 5.5), constrained_layout=True)
    
    # Define a clean, high-contrast color palette
    colors_dist = {
        20: '#3f51b5',   # Indigo
        80: '#00bcd4',   # Cyan
        140: '#4caf50',  # Green
        200: '#e91e63'   # Pink/Magenta
    }
    
    # ------------------ PANEL 1: COMPARISON VS DISTANCE ------------------
    ax1 = axes[0]
    # Plot paper values
    ax1.plot(paper_distances, paper_cooling_rates, 'ro--', label='Abdali et al. (2024) Paper', linewidth=2.5, markersize=8)
    
    # Plot simulated averages
    ax1.plot(paper_distances, sim_crs_avg, 'bs-', label='Simulated Spatial (G * R_solid) Avg', linewidth=2, markersize=8)
    ax1.plot(paper_distances, sim_crt_avg, 'g^-', label='Simulated Temporal (|dT/dt|) Avg', linewidth=2, markersize=8)
    
    # Plot all transient data points to show distribution/spread
    first_scatter = True
    for dist in paper_distances:
        lbl_crs = 'Simulated Spatial (Transient)' if first_scatter else None
        lbl_crt = 'Simulated Temporal (Transient)' if first_scatter else None
        
        # Spatial transient points
        ax1.scatter([dist] * len(df), df[f"CRS_{dist}um"], color='#82b1ff', alpha=0.25, s=15, edgecolors='none', label=lbl_crs)
        # Temporal transient points
        ax1.scatter([dist] * len(df), df[f"CRT_{dist}um"], color='#b9f6ca', alpha=0.25, s=15, edgecolors='none', label=lbl_crt)
        first_scatter = False

    ax1.set_xlabel('Distance from Melt Pool Bottom (um)', fontsize=11, fontweight='bold')
    ax1.set_ylabel('Solidification Cooling Rate (K/s)', fontsize=11, fontweight='bold')
    ax1.set_title('Cooling Rate vs Distance from Bottom', fontsize=12, fontweight='bold')
    ax1.grid(True, linestyle='--', alpha=0.5)
    ax1.legend(fontsize=9, loc='upper right')
    ax1.ticklabel_format(style='sci', axis='y', scilimits=(0,0))

    # ------------------ PANEL 2: TRANSIENT SPATIAL HISTORY ------------------
    ax2 = axes[1]
    time_ms = df["Time_s"] * 1000.0
    for dist in paper_distances:
        ax2.plot(time_ms, df[f"CRS_{dist}um"], color=colors_dist[dist], linewidth=2, label=f'{dist} um')
        
    ax2.set_xlabel('Time (ms)', fontsize=11, fontweight='bold')
    ax2.set_ylabel('Spatial Cooling Rate (K/s)', fontsize=11, fontweight='bold')
    ax2.set_title('Spatial Cooling Rate (G * R_solid) vs Time', fontsize=12, fontweight='bold')
    ax2.grid(True, linestyle='--', alpha=0.5)
    ax2.legend(fontsize=9, title="Distance")
    ax2.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
    # Add vertical line for start of steady state
    ax2.axvline(x=t_start_steady * 1000.0, color='gray', linestyle=':', linewidth=1.5, label='Steady-state start')

    # ------------------ PANEL 3: TRANSIENT TEMPORAL HISTORY ------------------
    ax3 = axes[2]
    for dist in paper_distances:
        ax3.plot(time_ms, df[f"CRT_{dist}um"], color=colors_dist[dist], linestyle='--', linewidth=2, label=f'{dist} um')
        
    ax3.set_xlabel('Time (ms)', fontsize=11, fontweight='bold')
    ax3.set_ylabel('Temporal Cooling Rate (K/s)', fontsize=11, fontweight='bold')
    ax3.set_title('Temporal Cooling Rate (|dT/dt|) vs Time', fontsize=12, fontweight='bold')
    ax3.grid(True, linestyle='--', alpha=0.5)
    ax3.legend(fontsize=9, title="Distance")
    ax3.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
    # Add vertical line for start of steady state
    ax3.axvline(x=t_start_steady * 1000.0, color='gray', linestyle=':', linewidth=1.5, label='Steady-state start')

    # Save plots
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUT_PLOT_PNG, dpi=200, bbox_inches='tight')
    fig.savefig(OUT_PLOT_PDF, bbox_inches='tight')
    plt.close(fig)

    print(f"Comparison plots saved successfully to:")
    print(f"  - PNG: {OUT_PLOT_PNG}")
    print(f"  - PDF: {OUT_PLOT_PDF}")

if __name__ == "__main__":
    main()
