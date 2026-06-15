#!/usr/bin/env python3
import os
import re
import sys
import json
import time
import shutil
import signal
import subprocess
from pathlib import Path
import matplotlib.pyplot as plt

# Paths
BENCHMARK_DIR = Path(__file__).resolve().parent
CASE_DIR = BENCHMARK_DIR.parent
RESULTS_JSON = BENCHMARK_DIR / "benchmark_results.json"

# Experimental parameters
CORE_CONFIGS = [4, 8, 16, 24, 32]
END_TIME = 5.0e-3  # 5.0 ms (proper full end time)

# Global results list to allow access from signal handlers
benchmark_results = []

def get_log_projection(log_path: Path, target_end_time: float) -> tuple:
    """
    Parses OpenFOAM solver log to get:
    - Measured/projected wall-clock solver time (s)
    - Simulated time progress (%)
    - Whether the run completed successfully (bool)
    """
    if not log_path.exists():
        return None, 0.0, False
        
    latest_time = 0.0
    latest_exec_time = 0.0
    is_finished = False
    
    time_pat = re.compile(r"^Time\s*=\s*([\d\.e\-+]+)")
    exec_pat = re.compile(r"ExecutionTime\s*=\s*([\d\.e\-+]+)")
    
    with open(log_path, "r", errors="ignore") as f:
        for line in f:
            if "Finalising parallel run" in line or "End\n" in line or line.strip() == "End":
                is_finished = True
            m_time = time_pat.match(line)
            if m_time:
                latest_time = float(m_time.group(1))
            m_exec = exec_pat.search(line)
            if m_exec:
                latest_exec_time = float(m_exec.group(1))
                
    if is_finished and latest_exec_time > 0.0:
        return latest_exec_time, 100.0, True
        
    if latest_time > 0.0 and latest_exec_time > 0.0:
        progress = (latest_time / target_end_time) * 100.0
        # Prevent division by tiny progress skewing projection initially
        if progress < 0.1:
            return None, progress, False
        projected_time = (target_end_time / latest_time) * latest_exec_time
        return projected_time, progress, False
        
    return None, 0.0, False

def get_updated_results() -> list:
    """
    Combines completed in-memory results with parsed log file projections
    for any currently running or incomplete configurations.
    """
    updated = []
    completed_cores = {r["cores"]: r for r in benchmark_results if r["status"] == "Completed"}
    
    for cores in CORE_CONFIGS:
        if cores in completed_cores:
            updated.append(completed_cores[cores])
            continue
            
        # Parse log file for incomplete/active cases
        log_path = BENCHMARK_DIR / f"run_{cores}W" / "log.laserbeamFoam"
        proj_time, progress, is_finished = get_log_projection(log_path, END_TIME)
        
        if proj_time is not None:
            status = "Completed" if is_finished else "Running (Projected)"
            updated.append({
                "cores": cores,
                "solver_time_s": proj_time,
                "progress_pct": progress,
                "status": status
            })
        else:
            updated.append({
                "cores": cores,
                "solver_time_s": None,
                "progress_pct": progress,
                "status": "Running (Warmup)" if progress > 0.0 else "Not Started"
            })
    return updated

def generate_plots(results: list):
    valid_results = [r for r in results if r["solver_time_s"] is not None]
    if not valid_results:
        return
        
    cores_list = [r["cores"] for r in valid_results]
    times_list = [r["solver_time_s"] for r in valid_results]
    
    base_time = next((r["solver_time_s"] for r in valid_results if r["cores"] == CORE_CONFIGS[0]), valid_results[0]["solver_time_s"])
    
    fig, ax1 = plt.subplots(figsize=(9, 5.5))
    
    # Plot Execution Time (Bars)
    color = '#1f77b4'
    ax1.set_xlabel('Number of Cores', fontweight='bold', fontsize=11)
    ax1.set_ylabel('Solver Execution Time (s)', color=color, fontweight='bold', fontsize=11)
    
    labels = []
    for r in valid_results:
        lbl = str(r["cores"])
        if "projected" in r["status"].lower():
            lbl += "\n(Projected)"
        labels.append(lbl)
        
    bars = ax1.bar(labels, times_list, color=color, alpha=0.7, width=0.4, label='Exec Time')
    ax1.tick_params(axis='y', labelcolor=color)
    ax1.grid(True, linestyle='--', alpha=0.3)
    
    # Add values on top of bars
    for bar in bars:
        height = bar.get_height()
        ax1.annotate(f'{height:.1f}s',
                    xy=(bar.get_x() + bar.get_width() / 2, height),
                    xytext=(0, 3),
                    textcoords="offset points",
                    ha='center', va='bottom', fontweight='bold')
                    
    # Plot Speedup (Line)
    ax2 = ax1.twinx()  
    color = '#d62728'
    ax2.set_ylabel('Speedup Factor', color=color, fontweight='bold', fontsize=11)
    speedups = [base_time / t for t in times_list]
    ax2.plot(labels, speedups, color=color, marker='o', linewidth=2, markersize=8, label='Speedup')
    ax2.tick_params(axis='y', labelcolor=color)
    
    # Theoretical ideal speedup line
    ideal_cores = [c / CORE_CONFIGS[0] for c in cores_list]
    ax2.plot(labels, ideal_cores, color='gray', linestyle=':', label='Ideal Speedup')
    
    plt.title(f'OpenFOAM Parallel Scaling (Full {END_TIME*1e3:.1f} ms SLM track)', fontweight='bold', fontsize=12)
    fig.tight_layout()
    
    plot_path = BENCHMARK_DIR / "benchmark_scaling.png"
    plt.savefig(plot_path, dpi=200)
    plt.close()

def generate_reports(results: list, final=False):
    # Sort by core count
    results = sorted(results, key=lambda x: x["cores"])
    
    # Find base time (usually 4 cores)
    base_time = next((r["solver_time_s"] for r in results if r["cores"] == CORE_CONFIGS[0] and r["solver_time_s"] is not None), None)
    if base_time is None:
        # Fallback to the first available core count
        valid_times = [r["solver_time_s"] for r in results if r["solver_time_s"] is not None]
        if valid_times:
            base_time = valid_times[0]
            
    summary_path = BENCHMARK_DIR / "benchmark_summary.md"
    with open(summary_path, "w") as sf:
        sf.write("# Parallel Scaling Benchmark Summary\n\n")
        status_str = "Completed" if final else "Running (Incremental / Projected)"
        sf.write(f"Status: **{status_str}** | Last Updated: {time.strftime('%Y-%m-%d %H:%M:%S')}\n\n")
        sf.write(f"This scaling test evaluates the performance of the `laserbeamFoam` solver on the autogenous SLM 316L single-track case (Abdali parameters) using core configurations: **{', '.join(map(str, CORE_CONFIGS))} cores** for a full simulation duration of **{END_TIME*1e3:.1f} ms**.\n\n")
        
        sf.write("## Performance Results Table\n\n")
        sf.write("| Cores | Status | Sim Progress | Solver Time (s) | Speedup vs Base | Parallel Efficiency |\n")
        sf.write("| :---: | :---: | :---: | :---: | :---: | :---: |\n")
        
        for r in results:
            cores = r["cores"]
            status = r["status"]
            progress = r["progress_pct"]
            time_val = r["solver_time_s"]
            
            if time_val is not None:
                time_str = f"{time_val:.1f}s"
                if "projected" in status.lower():
                    time_str += " (projected)"
                
                if base_time:
                    speedup = base_time / time_val
                    efficiency = (speedup / (cores / CORE_CONFIGS[0])) * 100.0
                    speedup_str = f"{speedup:.2f}x"
                    eff_str = f"{efficiency:.1f}%"
                else:
                    speedup_str = "1.00x"
                    eff_str = "100.0%"
            else:
                time_str = "N/A"
                speedup_str = "N/A"
                eff_str = "N/A"
                
            sf.write(f"| {cores} | {status} | {progress:.1f}% | {time_str} | {speedup_str} | {eff_str} |\n")
            
        sf.write("\n## Key Findings\n")
        valid_results = [r for r in results if r["solver_time_s"] is not None]
        if valid_results:
            fastest_config = min(valid_results, key=lambda x: x["solver_time_s"])
            sf.write(f"- **Fastest Configuration**: Running on **{fastest_config['cores']} cores** is the fastest, with a measured or projected time of **{fastest_config['solver_time_s']:.2f} seconds**.\n")
        sf.write("\n*Note: Projected times are calculated linearly based on solver progress to account for jobs that are cut off by the cluster time limit before finishing all configurations.*\n")

    # Save incremental JSON
    with open(RESULTS_JSON, "w") as jf:
        json.dump(results, jf, indent=4)

def sigterm_handler(signum, frame):
    print("\n[BENCHMARK] Received termination signal (SIGTERM/SIGINT).")
    print("[BENCHMARK] Writing final scaling reports based on current simulation progress...")
    
    updated_results = get_updated_results()
    generate_reports(updated_results, final=False)
    generate_plots(updated_results)
    
    print("[BENCHMARK] Scaling reports generated successfully. Exiting.")
    sys.exit(0)

# Register signals
signal.signal(signal.SIGTERM, sigterm_handler)
signal.signal(signal.SIGINT, sigterm_handler)

def setup_case(cores: int) -> Path:
    run_dir = BENCHMARK_DIR / f"run_{cores}W"
    print(f"\n[BENCHMARK] Setting up case for {cores} cores in {run_dir.name}...")
    
    if run_dir.exists():
        shutil.rmtree(run_dir)
    run_dir.mkdir(parents=True, exist_ok=True)
    
    shutil.copytree(CASE_DIR / "initial", run_dir / "initial")
    shutil.copytree(CASE_DIR / "constant", run_dir / "constant")
    shutil.copytree(CASE_DIR / "system", run_dir / "system")
    shutil.copytree(run_dir / "initial", run_dir / "0")
    
    # Modify controlDict: set endTime=5.0e-3, writeInterval=1.0 (reduce I/O overhead)
    control_dict_path = run_dir / "system" / "controlDict"
    with open(control_dict_path, "r") as f:
        content = f.read()
    content = re.sub(r"endTime\s+[\d\.e\-+]+;", f"endTime {END_TIME};", content)
    content = re.sub(r"writeInterval\s+[\d\.e\-+]+;", "writeInterval 1.0;", content)
    with open(control_dict_path, "w") as f:
        f.write(content)
        
    # Modify decomposeParDict: set numberOfSubdomains
    decomp_dict_path = run_dir / "system" / "decomposeParDict"
    with open(decomp_dict_path, "r") as f:
        content = f.read()
    content = re.sub(r"numberOfSubdomains\s+\d+;", f"numberOfSubdomains {cores};", content)
    with open(decomp_dict_path, "w") as f:
        f.write(content)
        
    return run_dir

def run_cmd(cmd: list, cwd: Path, log_file=None):
    if log_file:
        with open(log_file, "w") as lf:
            subprocess.run(cmd, cwd=str(cwd), stdout=lf, stderr=subprocess.STDOUT, check=True)
    else:
        subprocess.run(cmd, cwd=str(cwd), stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, check=True)

def main():
    # Load previous partial results if they exist
    global benchmark_results
    if RESULTS_JSON.exists():
        try:
            with open(RESULTS_JSON, "r") as jf:
                benchmark_results = json.load(jf)
            print("[BENCHMARK] Loaded existing results from benchmark_results.json.")
        except Exception:
            pass
            
    completed_cores = {r["cores"] for r in benchmark_results if r["status"] == "Completed"}
    
    for cores in CORE_CONFIGS:
        if cores in completed_cores:
            print(f"[BENCHMARK] Cores {cores} already completed. Skipping.")
            continue
            
        try:
            run_dir = setup_case(cores)
            
            # Setup log path
            log_path = run_dir / "log.laserbeamFoam"
            
            print(f"[BENCHMARK] [{cores} Cores] Running blockMesh...")
            run_cmd(["blockMesh"], run_dir, run_dir / "log.blockMesh")
            
            print(f"[BENCHMARK] [{cores} Cores] Running setFields...")
            run_cmd(["setFields"], run_dir, run_dir / "log.setFields")
            
            print(f"[BENCHMARK] [{cores} Cores] Decomposing domain...")
            run_cmd(["decomposePar"], run_dir, run_dir / "log.decomposePar")
            
            print(f"[BENCHMARK] [{cores} Cores] Running laserbeamFoam solver (Full {END_TIME*1e3:.1f}ms)...")
            solver_cmd = ["mpirun", "--oversubscribe", "-np", str(cores), "laserbeamFoam", "-parallel"]
            
            start_time = time.time()
            run_cmd(solver_cmd, run_dir, log_path)
            solver_time = time.time() - start_time
            
            print(f"[BENCHMARK] [{cores} Cores] Solver execution completed in {solver_time:.2f} seconds.")
            
            # Update results
            # Remove any previous entries for this core configuration
            benchmark_results = [r for r in benchmark_results if r["cores"] != cores]
            benchmark_results.append({
                "cores": cores,
                "solver_time_s": solver_time,
                "progress_pct": 100.0,
                "status": "Completed"
            })
            
            # Regenerate reports and plots incrementally
            updated_results = get_updated_results()
            generate_reports(updated_results, final=(len(completed_cores) + 1 == len(CORE_CONFIGS)))
            generate_plots(updated_results)
            
            # Cleanup mesh/processors to save disk space
            for p in run_dir.iterdir():
                if p.is_dir() and p.name not in ["initial"]:
                    shutil.rmtree(p)
                elif p.is_file() and not p.name.startswith("log."):
                    p.unlink()
                    
        except Exception as e:
            print(f"[ERROR] Benchmark failed for {cores} cores: {e}")
            
    # Final generation
    updated_results = get_updated_results()
    generate_reports(updated_results, final=True)
    generate_plots(updated_results)
    print("\n[BENCHMARK] Scaling benchmark completed successfully.")

if __name__ == "__main__":
    main()
