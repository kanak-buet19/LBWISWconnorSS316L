#!/usr/bin/env python3
import sys
from pathlib import Path
import numpy as np

# Set non-interactive backend for matplotlib
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

def parse_summary_csv(csv_path):
    import csv
    data = []
    with open(csv_path, 'r') as f:
        reader = csv.DictReader(f)
        for row in reader:
            data.append({
                'Resistivity': float(row['Resistivity']),
                'Beta_r': float(row['Beta_r']),
                'e_num_density': float(row['e_num_density']),
                'LatentHeatVap': float(row['LatentHeatVap']),
                'Depth': float(row['Keyhole_Depth_um']),
                'Temp': float(row['Max_Temp_K']),
                'Pressure': float(row['Max_Pressure_kPa'])
            })
    return data

def main():
    study_dir = Path(__file__).resolve().parents[1]
    csv_path = study_dir / "results" / "summary.csv"
    output_dir = study_dir / "plots"
    
    if not csv_path.exists():
        print(f"Error: Summary CSV not found at {csv_path}", file=sys.stderr)
        sys.exit(1)
        
    data = parse_summary_csv(csv_path)
    if not data:
        print("Error: No data rows found in summary CSV", file=sys.stderr)
        sys.exit(1)
        
    output_dir.mkdir(parents=True, exist_ok=True)
    
    # Convert data to numpy arrays for calculation
    depths = np.array([row['Depth'] for row in data])
    res = np.array([row['Resistivity'] for row in data])
    beta = np.array([row['Beta_r'] for row in data])
    enum = np.array([row['e_num_density'] for row in data])
    lh = np.array([row['LatentHeatVap'] for row in data])
    
    # Calculate Pearson Correlation Coefficients
    corrs = {
        'e_num_density': np.corrcoef(enum, depths)[0, 1] if len(np.unique(enum)) > 1 else 0.0,
        'elec_resistivity': np.corrcoef(res, depths)[0, 1] if len(np.unique(res)) > 1 else 0.0,
        'beta_r': np.corrcoef(beta, depths)[0, 1] if len(np.unique(beta)) > 1 else 0.0,
        'LatentHeatVap': np.corrcoef(lh, depths)[0, 1] if len(np.unique(lh)) > 1 else 0.0
    }
    
    # Create the beautiful 4-panel dashboard
    fig, axs = plt.subplots(2, 2, figsize=(14, 11))
    
    # Color palette
    colors = ['#4C78A8', '#F58518', '#E45756', '#72B7B2']
    
    # 1. Top-Left: e_num_density vs Depth
    axs[0, 0].scatter(enum, depths, color=colors[0], s=60, edgecolor='black', zorder=3)
    if len(np.unique(enum)) > 1:
        m, b = np.polyfit(enum, depths, 1)
        axs[0, 0].plot(enum, m*enum + b, color=colors[0], linestyle='--', alpha=0.8)
    axs[0, 0].set_title('Laser Electron Density (e_num_density) vs Depth', fontsize=11, fontweight='bold')
    axs[0, 0].set_xlabel('e_num_density (m^-3)', fontsize=9)
    axs[0, 0].set_ylabel('Keyhole Depth (um)', fontsize=9)
    axs[0, 0].grid(True, linestyle=':', alpha=0.5)
    
    # 2. Top-Right: elec_resistivity vs Depth
    axs[0, 1].scatter(res * 1e7, depths, color=colors[1], s=60, edgecolor='black', zorder=3)
    if len(np.unique(res)) > 1:
        m, b = np.polyfit(res * 1e7, depths, 1)
        axs[0, 1].plot(res * 1e7, m*(res*1e7) + b, color=colors[1], linestyle='--', alpha=0.8)
    axs[0, 1].set_title('Electrical Resistivity (elec_resistivity) vs Depth', fontsize=11, fontweight='bold')
    axs[0, 1].set_xlabel('elec_resistivity (x10^-7 Ohm-m)', fontsize=9)
    axs[0, 1].set_ylabel('Keyhole Depth (um)', fontsize=9)
    axs[0, 1].grid(True, linestyle=':', alpha=0.5)
    
    # 3. Bottom-Left: beta_r vs Depth
    axs[1, 0].scatter(beta, depths, color=colors[2], s=60, edgecolor='black', zorder=3)
    if len(np.unique(beta)) > 1:
        m, b = np.polyfit(beta, depths, 1)
        axs[1, 0].plot(beta, m*beta + b, color=colors[2], linestyle='--', alpha=0.8)
    axs[1, 0].set_title('Recoil Pressure Coefficient (beta_r) vs Depth', fontsize=11, fontweight='bold')
    axs[1, 0].set_xlabel('beta_r', fontsize=9)
    axs[1, 0].set_ylabel('Keyhole Depth (um)', fontsize=9)
    axs[1, 0].grid(True, linestyle=':', alpha=0.5)
    
    # 4. Bottom-Right: LatentHeatVap vs Depth
    axs[1, 1].scatter(lh / 1e6, depths, color=colors[3], s=60, edgecolor='black', zorder=3)
    if len(np.unique(lh)) > 1:
        m, b = np.polyfit(lh / 1e6, depths, 1)
        axs[1, 1].plot(lh / 1e6, m*(lh/1e6) + b, color=colors[3], linestyle='--', alpha=0.8)
    axs[1, 1].set_title('Latent Heat of Vaporization (LatentHeatVap) vs Depth', fontsize=11, fontweight='bold')
    axs[1, 1].set_xlabel('LatentHeatVap (MJ/kg)', fontsize=9)
    axs[1, 1].set_ylabel('Keyhole Depth (um)', fontsize=9)
    axs[1, 1].grid(True, linestyle=':', alpha=0.5)
    
    plt.suptitle('LPBF Keyhole Sensitivity Analysis Dashboard', fontsize=14, fontweight='bold', y=0.98)
    plt.tight_layout(rect=[0, 0.03, 1, 0.95])
    
    dashboard_plot_path = output_dir / "sensitivity_dashboard.png"
    plt.savefig(dashboard_plot_path, dpi=180)
    plt.close()
    
    # 5. Create a standalone relative impact bar chart showing Pearson Correlation
    plt.figure(figsize=(7.5, 5.5))
    features = list(corrs.keys())
    values = [corrs[f] for f in features]
    
    bar_colors = [colors[i] for i in range(len(features))]
    bars = plt.barh(features, values, color=bar_colors, edgecolor='black', height=0.55)
    plt.axvline(x=0.0, color='black', linestyle='-', linewidth=1.1)
    
    # Annotate bar values
    for bar in bars:
        width = bar.get_width()
        x_pos = width + 0.02 if width >= 0 else width - 0.08
        plt.annotate(
            f"{width:+.2f}",
            xy=(x_pos, bar.get_y() + bar.get_height() / 2),
            va='center',
            fontweight='bold',
            fontsize=9
        )
        
    plt.title('Relative Sensitivity of Keyhole Depth (Pearson Correlation)', fontsize=12, fontweight='bold', pad=15)
    plt.xlabel('Correlation Coefficient (Negative = Shallower, Positive = Deeper)', fontsize=10)
    plt.xlim(-1.15, 1.15)
    plt.grid(True, axis='x', linestyle=':', alpha=0.5)
    plt.tight_layout()
    
    corr_plot_path = output_dir / "parameter_importance.png"
    plt.savefig(corr_plot_path, dpi=150)
    plt.close()
    
    print(f"Successfully generated main sensitivity dashboard: {dashboard_plot_path}")
    print(f"Successfully generated parameter correlation chart: {corr_plot_path}")

if __name__ == "__main__":
    main()
