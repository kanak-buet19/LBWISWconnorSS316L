#!/usr/bin/env python3
"""
LPBF Solver Sensitivity Study Orchestrator
Uses OpenFOAM-v2412 + laserbeamFoam solver to sweep physical parameters.
Runs 4 concurrent cases with 4 cores each, displaying multi-level tqdm progress bars.
"""

import os
import re
import sys
import json
import queue
import shutil
import itertools
import threading
import csv
import subprocess
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor
from tqdm import tqdm

# OpenFOAM Sourcing Command
OF_SOURCE_CMD = "source /usr/lib/openfoam/openfoam2412/etc/bashrc"

def get_field_maxima(case_dir: Path):
    """
    Parses unified fieldMinMax.dat file to find absolute maximum T and recoilPressure.
    Returns: (max_T, max_P_kPa)
    """
    dat_files = sorted(case_dir.glob("postProcessing/fieldMinMax/*/fieldMinMax.dat"))
    if not dat_files:
        return 300.0, 101.325 # Default ambient values
    
    max_T = 300.0
    max_P = 101325.0 # Pa
    
    with open(dat_files[-1], 'r', errors='ignore') as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith('#'):
                continue
            parts = line.split()
            if len(parts) >= 8:
                try:
                    field_name = parts[1]
                    max_val = float(parts[7])
                    if field_name == 'T':
                        if max_val > max_T:
                            max_T = max_val
                    elif field_name == 'recoilPressure':
                        if max_val > max_P:
                            max_P = max_val
                except ValueError:
                    continue
    return max_T, max_P / 1e3 # Convert Pa to kPa

def setup_case_dir(template_dir: Path, case_dir: Path, params: dict):
    """
    Copies system, constant, initial, scripts, and foamVTK.sh to the sandbox.
    Recreates 0/ folder from initial/, copies .orig fields, and edits parameters.
    """
    if case_dir.exists():
        shutil.rmtree(case_dir)
    case_dir.mkdir(parents=True, exist_ok=True)
    
    # Copy template folders
    for folder in ["system", "constant", "initial", "scripts"]:
        shutil.copytree(template_dir / folder, case_dir / folder, dirs_exist_ok=True)
    shutil.copy2(template_dir / "foamVTK.sh", case_dir / "foamVTK.sh")
    
    # Recreate 0/ folder from initial/
    zeros_dir = case_dir / "0"
    shutil.copytree(case_dir / "initial", zeros_dir, dirs_exist_ok=True)
    
    # Resolve .orig files in 0/
    for orig_file in zeros_dir.glob("**/*.orig"):
        dest_file = orig_file.with_suffix("")
        shutil.copy2(orig_file, dest_file)

    # 1. Edit LaserProperties
    lp_path = case_dir / "constant" / "LaserProperties"
    lp_content = lp_path.read_text()
    lp_content = re.sub(
        r'(e_num_density\s+)[0-9.eE+-]+;',
        rf'\g<1>{params["e_num_density"]:g};',
        lp_content
    )
    lp_path.write_text(lp_content)

    # 2. Edit transportProperties
    tp_path = case_dir / "constant" / "transportProperties"
    tp_content = tp_path.read_text()
    tp_content = re.sub(
        r'(elec_resistivity\s+)[0-9.eE+-]+;',
        rf'\g<1>{params["elec_resistivity"]:g};',
        tp_content
    )
    tp_content = re.sub(
        r'(beta_r\s+)[0-9.eE+-]+;',
        rf'\g<1>{params["beta_r"]:g};',
        tp_content
    )
    tp_content = re.sub(
        r'(LatentHeatVap\s+)[0-9.eE+-]+;',
        rf'\g<1>{params["LatentHeatVap"]:g};',
        tp_content
    )
    tp_path.write_text(tp_content)

    # 3. Edit controlDict (endTime and debugNumerics=false)
    cd_path = case_dir / "system" / "controlDict"
    cd_content = cd_path.read_text()
    cd_content = re.sub(
        r'(endTime\s+)[0-9.eE+-]+;',
        rf'\g<1>{params["endTime"]:g};',
        cd_content
    )
    cd_content = re.sub(
        r'(debugNumerics\s+)\w+;',
        r'\g<1>false;',
        cd_content
    )
    cd_path.write_text(cd_content)

    # 4. Edit decomposeParDict
    dp_path = case_dir / "system" / "decomposeParDict"
    dp_content = dp_path.read_text()
    dp_content = re.sub(
        r'(numberOfSubdomains\s+)\d+;',
        rf'\g<1>{params["nCoresPerRun"]};',
        dp_content
    )
    dp_path.write_text(dp_content)

def run_cmd_failsafe(cmd: str, cwd: Path, log_files=None, desc=""):
    """
    Runs a shell command synchronously. Redirects output to log files.
    If the exit code is non-zero, immediately halts the entire orchestrator script.
    """
    try:
        p = subprocess.Popen(
            cmd,
            shell=True,
            executable="/bin/bash",
            cwd=cwd,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            text=True
        )
        
        stdout_buf = []
        for line in p.stdout:
            stdout_buf.append(line)
            if log_files:
                for f in log_files:
                    f.write(line)
                    f.flush()
                    
        p.wait()
        if p.returncode != 0:
            err_msg = "".join(stdout_buf)
            print(f"\n\n[FAILSAFE] CRITICAL ERROR: Step '{desc}' failed in {cwd} with exit code {p.returncode}!", file=sys.stderr)
            print(f"--- COMMAND ---\n{cmd}\n--- STDOUT/STDERR ---\n{err_msg}", file=sys.stderr)
            # Instantly terminate the orchestrator process to prevent further case runs
            os._exit(1)
    except Exception as e:
        print(f"\n\n[FAILSAFE] CRITICAL EXCEPTION running step '{desc}': {e}", file=sys.stderr)
        os._exit(1)

def run_solver_with_tqdm(case_name: str, cwd: Path, end_time: float, worker_id: int, log_files):
    """
    Launches mpirun solver and parses stdout line-by-line in real-time
    to drive an isolated tqdm progress bar. Halts the orchestrator on crash.
    """
    cmd = f"{OF_SOURCE_CMD} && mpirun --oversubscribe -np 4 laserbeamFoam -parallel"
    
    # Initialize isolated stacked progress bar
    pbar = tqdm(
        total=end_time,
        desc=f"Worker {worker_id} | {case_name[:25]}",
        position=worker_id + 1,
        leave=True,
        unit="s",
        bar_format="{desc}: {percentage:3.0f}%|{bar}| {n:.2e}/{total:.2e} [{elapsed}<{remaining}]"
    )
    
    try:
        p = subprocess.Popen(
            cmd,
            shell=True,
            executable="/bin/bash",
            cwd=cwd,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            text=True
        )
        
        time_pattern = re.compile(r'^Time\s*=\s*([0-9.eE+-]+)')
        stdout_buf = []
        for line in p.stdout:
            stdout_buf.append(line)
            # Save to log files
            for f in log_files:
                f.write(line)
                f.flush()
            
            # Update progress bar
            m = time_pattern.match(line.strip())
            if m:
                t_val = float(m.group(1))
                pbar.n = min(t_val, end_time)
                pbar.refresh()
                
        p.wait()
        pbar.n = end_time
        pbar.refresh()
        pbar.close()
        
        if p.returncode != 0:
            err_msg = "".join(stdout_buf)
            print(f"\n\n[FAILSAFE] CRITICAL ERROR: Solver crashed for {case_name} with exit status {p.returncode}!", file=sys.stderr)
            print(f"--- STDOUT/STDERR ---\n{err_msg}", file=sys.stderr)
            os._exit(1)
            
    except Exception as e:
        pbar.close()
        print(f"\n\n[FAILSAFE] CRITICAL EXCEPTION running solver for {case_name}: {e}", file=sys.stderr)
        os._exit(1)

def execute_case(idx: int, combo: tuple, keys: list, study_dir: Path, template_dir: Path, end_time: float, slot_queue: queue.Queue, summary_csv_path: Path, csv_lock: threading.Lock):
    """
    Task worker executing a single case sweep in full isolation.
    """
    # Get available screen slot and concurrent sandbox index
    worker_id = slot_queue.get()
    
    params = dict(zip(keys, combo))
    params["endTime"] = end_time
    params["nCoresPerRun"] = 4 # Decomposed to exactly 4 cores
    
    # Beautiful unique directory name
    case_name = f"case_dens{params['e_num_density']:.1e}_res{params['elec_resistivity']:.1e}_beta{params['beta_r']:.2f}_lh{params['LatentHeatVap']:.2e}".replace("+", "")
    case_dir = study_dir / "runs" / case_name
    central_log_path = study_dir / "logs" / f"{case_name}.log"
    
    # Setup directories
    case_dir.mkdir(parents=True, exist_ok=True)
    setup_case_dir(template_dir, case_dir, params)
    
    # Open local and central log files
    local_log_file = open(case_dir / "log.laserbeamFoam", "w")
    central_log_file = open(central_log_path, "w")
    log_files = [local_log_file, central_log_file]
    
    try:
        # Preprocessing steps
        run_cmd_failsafe(f"{OF_SOURCE_CMD} && blockMesh", cwd=case_dir, log_files=log_files, desc="blockMesh")
        run_cmd_failsafe(f"{OF_SOURCE_CMD} && setFields", cwd=case_dir, log_files=log_files, desc="setFields")
        run_cmd_failsafe(f"{OF_SOURCE_CMD} && decomposePar", cwd=case_dir, log_files=log_files, desc="decomposePar")
        
        # Run standard parallel simulation with tqdm progress bar
        run_solver_with_tqdm(case_name, case_dir, end_time, worker_id, log_files)
        
        # Post-processing: foamVTK and meltpool analysis
        run_cmd_failsafe(f"{OF_SOURCE_CMD} && ./foamVTK.sh all", cwd=case_dir, log_files=log_files, desc="foamVTK")
        run_cmd_failsafe(f"{OF_SOURCE_CMD} && {sys.executable} scripts/analyze_meltpool_vtu.py --case .", cwd=case_dir, log_files=log_files, desc="analyze_meltpool_vtu")
        
        # Organise meltpool analysis plots
        vtu_sections_dir = case_dir / "post-processing-data" / "vtu_sections"
        png_files = sorted(vtu_sections_dir.glob("*.png"))
        if png_files:
            shutil.copy2(png_files[-1], study_dir / "plots" / "meltpool_analysis" / f"{case_name}.png")
            
        # Run recoil history plots
        recoil_cmd = f"{OF_SOURCE_CMD} && {sys.executable} scripts/plot_recoil_history.py --case-dir . --output-dir ../../plots/recoil_pressure_plots --case-name {case_name}"
        run_cmd_failsafe(recoil_cmd, cwd=case_dir, log_files=log_files, desc="plot_recoil_history")
        
        # Extract steady state results
        depth = 0.0
        geom_csv = case_dir / "post-processing-data" / "vtu_meltpool_geometry.csv"
        if geom_csv.exists():
            try:
                with open(geom_csv, 'r') as f:
                    reader = csv.DictReader(f)
                    rows = list(reader)
                    if rows:
                        depth = float(rows[-1]['keyholeDepth_um'])
            except Exception:
                pass
                
        max_T, max_P_kPa = get_field_maxima(case_dir)
        
        # Thread-safe append to summary CSV
        with csv_lock:
            with open(summary_csv_path, 'a') as f:
                writer = csv.writer(f)
                writer.writerow([
                    params['elec_resistivity'],
                    params['beta_r'],
                    params['e_num_density'],
                    params['LatentHeatVap'],
                    depth,
                    max_T,
                    max_P_kPa
                ])
                
    finally:
        # Close open file handlers
        local_log_file.close()
        central_log_file.close()
        # Free screen worker slot
        slot_queue.put(worker_id)

def main():
    script_dir = Path(__file__).resolve().parent
    study_dir = script_dir.parent
    base_dir = study_dir.parent
    template_dir = base_dir / "case_original_thermal_prop"
    
    config_path = study_dir / "config.json"
    if not config_path.exists():
        print(f"Error: config.json not found at {config_path}", file=sys.stderr)
        sys.exit(1)
        
    with open(config_path, 'r') as f:
        config = json.load(f)
        
    sim_cfg = config["simulation"]
    end_time = sim_cfg["endTime"]
    max_concurrent = sim_cfg["maxConcurrentRuns"]
    
    params = config["parameters"]
    keys = list(params.keys())
    combinations = list(itertools.product(*(params[k] for k in keys)))
    
    print("======================================================================")
    print("Starting Concurrent LPBF Solver Sensitivity Study")
    print(f"Total Sweep Matrix Cases: {len(combinations)}")
    print(f"Concurrency Allocation  : {max_concurrent} concurrent runs @ 4 cores each")
    print("======================================================================")
    
    # Initialize workspace directories
    (study_dir / "runs").mkdir(parents=True, exist_ok=True)
    (study_dir / "logs").mkdir(parents=True, exist_ok=True)
    (study_dir / "plots" / "meltpool_analysis").mkdir(parents=True, exist_ok=True)
    (study_dir / "plots" / "recoil_pressure_plots").mkdir(parents=True, exist_ok=True)
    (study_dir / "results").mkdir(parents=True, exist_ok=True)
    
    # Initialise summary csv
    summary_csv = study_dir / "results" / "summary.csv"
    with open(summary_csv, 'w') as f:
        writer = csv.writer(f)
        writer.writerow([
            "Resistivity",
            "Beta_r",
            "e_num_density",
            "LatentHeatVap",
            "Keyhole_Depth_um",
            "Max_Temp_K",
            "Max_Pressure_kPa"
        ])
        
    # Queue for worker screen positions (0 to max_concurrent-1)
    slot_queue = queue.Queue()
    for i in range(max_concurrent):
        slot_queue.put(i)
        
    csv_lock = threading.Lock()
    
    # Run the concurrent thread pool executor
    with ThreadPoolExecutor(max_workers=max_concurrent) as executor:
        futures = [
            executor.submit(
                execute_case,
                idx,
                combo,
                keys,
                study_dir,
                template_dir,
                end_time,
                slot_queue,
                summary_csv,
                csv_lock
            )
            for idx, combo in enumerate(combinations)
        ]
        
        # Wait for all thread workers to complete (or crash/exit failsafe)
        for fut in futures:
            fut.result()
            
    print("\n======================================================================")
    print("Simulations complete. Running central dashboard script...")
    print("======================================================================")
    
    # Run plot_sensitivity_dashboard.py to generate final publication-grade dashboard
    dashboard_cmd = f"{sys.executable} scripts/plot_sensitivity_dashboard.py"
    subprocess.run(dashboard_cmd, shell=True, cwd=study_dir, check=True)
    
    print("\n======================================================================")
    print("LPBF Parameter Sensitivity Study Complete!")
    print(f"Organised Case Files Saved in  : {study_dir}/runs/")
    print(f"Simulation Central Logs in      : {study_dir}/logs/")
    print(f"Meltpool Geometry Plots in     : {study_dir}/plots/meltpool_analysis/")
    print(f"Recoil Pressure Stacked Plots in: {study_dir}/plots/recoil_pressure_plots/")
    print(f"Summary Comparison Table in     : {summary_csv}")
    print(f"Final Study Dashboard Saved as  : {study_dir}/plots/sensitivity_dashboard.png")
    print("======================================================================")

if __name__ == "__main__":
    main()
