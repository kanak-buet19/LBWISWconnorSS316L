#!/usr/bin/env python3
import os
import re
import shutil
import subprocess
import sys
from pathlib import Path

# targeted test matrix: [elec_resistivity, beta_r, sigma, name]
TEST_CASES = [
    (3.5e-7, 0.2,  0.91, "Baseline_LowAbs"),
    (6.0e-7, 0.2,  0.91, "HighAbs_Runaway"),
    (4.8e-7, 0.05, 0.91, "Goldilocks_LowRecoil"),
    (4.8e-7, 0.2,  1.2,  "Goldilocks_HighSigma"),
    (4.8e-7, 0.05, 1.2,  "Combined_Stabilized")
]

def backup_file(file_path):
    backup_path = file_path.with_suffix(file_path.suffix + ".bak")
    if file_path.exists() and not backup_path.exists():
        shutil.copy2(file_path, backup_path)
    return backup_path

def restore_file(file_path):
    backup_path = file_path.with_suffix(file_path.suffix + ".bak")
    if backup_path.exists():
        shutil.copy2(backup_path, file_path)
        os.remove(backup_path)

def modify_transport_properties(file_path, elec_res, beta_r, sigma):
    content = file_path.read_text()
    
    # 1. Update elec_resistivity
    content = re.sub(
        r"elec_resistivity\s+[0-9.eE+-]+;",
        f"elec_resistivity\t{elec_res:g};",
        content
    )
    
    # 2. Update beta_r
    content = re.sub(
        r"beta_r\s+[0-9.eE+-]+;",
        f"beta_r  {beta_r:g};",
        content
    )
    
    # 3. Update sigma
    content = re.sub(
        r"sigma\s+[0-9.eE+-]+;",
        f"sigma            {sigma:g};",
        content
    )
    
    file_path.write_text(content)

def modify_control_dict(file_path, end_time):
    content = file_path.read_text()
    
    # Update endTime
    content = re.sub(
        r"endTime\s+[0-9.eE+-]+;",
        f"endTime         {end_time:g};",
        content
    )
    
    # Ensure debugNumerics is false to avoid deadlock
    content = re.sub(
        r"debugNumerics\s+\w+;",
        "debugNumerics              false;",
        content
    )
    
    file_path.write_text(content)

def run_command(cmd, cwd):
    print(f"Running command: {' '.join(cmd)}")
    result = subprocess.run(cmd, cwd=cwd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    if result.returncode != 0:
        print(f"Command failed with exit code {result.returncode}")
        print(f"Stderr:\n{result.stderr}")
    return result.returncode == 0

def get_steady_state_results(case_dir):
    # Read keyhole depth from post-processing-data/vtu_meltpool_geometry.csv
    csv_path = case_dir / "post-processing-data" / "vtu_meltpool_geometry.csv"
    depth = 0.0
    if csv_path.exists():
        try:
            with open(csv_path, 'r') as f:
                lines = f.readlines()
                if len(lines) > 1:
                    last_line = lines[-1].strip().split(',')
                    # Index 4 is keyholeDepth_um
                    depth = float(last_line[4])
        except Exception as e:
            print(f"Warning: Failed to parse geometry CSV ({e})")

    # Read peak temp and pressure from postProcessing/fieldMinMax/0/fieldMinMax.dat
    dat_path = sorted(case_dir.glob("postProcessing/fieldMinMax/*/fieldMinMax.dat"))
    max_T = 300.0
    max_P = 101.3
    
    if dat_files := dat_path:
        try:
            with open(dat_files[-1], 'r') as f:
                for line in reversed(f.readlines()):
                    line = line.strip()
                    if not line or line.startswith('#'):
                        continue
                    parts = line.split()
                    if len(parts) >= 8:
                        t_val = float(parts[0])
                        field_name = parts[1]
                        val = float(parts[7])
                        if field_name == 'T' and max_T == 300.0:
                            max_T = val
                        elif field_name == 'recoilPressure' and max_P == 101.3:
                            max_P = val / 1e3  # Convert to kPa
                    if max_T != 300.0 and max_P != 101.3:
                        break
        except Exception as e:
            print(f"Warning: Failed to parse fieldMinMax dat file ({e})")
            
    return depth, max_T, max_P

def main():
    case_dir = Path(__file__).resolve().parents[1]
    
    tp_file = case_dir / "constant" / "transportProperties"
    cd_file = case_dir / "system" / "controlDict"
    
    print("==========================================================")
    print("Starting Automated Keyhole Sensitivity Study")
    print("==========================================================")
    
    # Back up critical dictionaries
    backup_file(tp_file)
    backup_file(cd_file)
    
    # Set fast 0.5 ms run time
    modify_control_dict(cd_file, 0.0005)
    
    results = []
    
    try:
        for elec_res, beta_r, sigma, name in TEST_CASES:
            print("\n----------------------------------------------------------")
            print(f"Configuring Test Case: {name}")
            print(f"Resistivity: {elec_res:g} | Beta: {beta_r:g} | Sigma: {sigma:g}")
            print("----------------------------------------------------------")
            
            modify_transport_properties(tp_file, elec_res, beta_r, sigma)
            
            # Clean previous run
            run_command(["./Allclean"], case_dir)
            
            # Run simulation
            # We temporarily bypass convertOnWrite queue in Allrun for speed
            # by executing OpenFOAM parallel simulation directly
            nCores_proc = subprocess.run(
                ["foamDictionary", "-entry", "numberOfSubdomains", "-value", "system/decomposeParDict"],
                cwd=case_dir, stdout=subprocess.PIPE, text=True
            )
            nCores = nCores_proc.stdout.strip()
            if not nCores:
                nCores = "8"
                
            print("Running preprocessing applications...")
            run_command(["blockMesh"], case_dir)
            run_command(["setFields"], case_dir)
            run_command(["decomposePar"], case_dir)
            
            print(f"Running parallel solver laserbeamFoam on {nCores} cores...")
            solver_cmd = ["mpirun", "--oversubscribe", "-np", nCores, "laserbeamFoam", "-parallel"]
            run_command(solver_cmd, case_dir)
            
            # Reconstruct and post-process
            print("Reconstructing and converting timesteps to legacy VTK...")
            run_command(["./foamVTK.sh", "all"], case_dir)
            
            print("Running meltpool geometry analysis...")
            python_exe = sys.executable
            run_command([python_exe, "scripts/analyze_meltpool_vtu.py"], case_dir)
            
            # Collect results
            depth, max_T, max_P = get_steady_state_results(case_dir)
            print(f"\nResults for {name}:")
            print(f"  Steady-State Keyhole Depth: {depth:.1f} um")
            print(f"  Max Domain Temperature: {max_T:.1f} K")
            print(f"  Max Recoil Pressure: {max_P:.1f} kPa")
            
            results.append((name, elec_res, beta_r, sigma, depth, max_T, max_P))
            
    finally:
        # Restore files to original
        restore_file(tp_file)
        restore_file(cd_file)
        print("\nRestored transportProperties and controlDict back to baseline.")
        
    # Write summary table
    summary_path = case_dir / "post-processing" / "sensitivity_study_results.csv"
    with open(summary_path, 'w') as f:
        f.write("Case_Name,Resistivity,Beta_r,Sigma_N/m,Keyhole_Depth_um,Max_Temp_K,Max_Pressure_kPa\n")
        for r in results:
            f.write(f"{r[0]},{r[1]:g},{r[2]:g},{r[3]:g},{r[4]:.2f},{r[5]:.2f},{r[6]:.2f}\n")
            
    print("\n==========================================================")
    print("Sensitivity Study Complete!")
    print(f"Summary table saved to: {summary_path}")
    print("==========================================================")
    
    # Print the table to terminal
    print(f"{'Case Name':<22} | {'Resist.':<8} | {'Beta':<5} | {'Sigma':<5} | {'Depth (um)':<10} | {'Temp (K)':<8} | {'Pres (kPa)':<10}")
    print("-" * 88)
    for r in results:
        print(f"{r[0]:<22} | {r[1]:<8.1e} | {r[2]:<5.2f} | {r[3]:<5.2f} | {r[4]:<10.1f} | {r[5]:<8.1f} | {r[6]:<10.1f}")

if __name__ == "__main__":
    main()
