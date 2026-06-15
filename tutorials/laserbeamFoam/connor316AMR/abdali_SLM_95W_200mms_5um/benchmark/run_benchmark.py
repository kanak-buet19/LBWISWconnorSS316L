#!/usr/bin/env python3
import os
import re
import shutil
import time
import subprocess
from pathlib import Path
import matplotlib.pyplot as plt

# Paths
BENCHMARK_DIR = Path(__file__).resolve().parent
CASE_DIR = BENCHMARK_DIR.parent
PYTHON_EXE = "/home/kanak/.venv/venv314/bin/python"

# Experimental parameters
CORE_CONFIGS = [4, 8, 16, 24, 32]
END_TIME = 9.0e-5  # 0.09 ms

def setup_case(cores: int) -> Path:
    run_dir = BENCHMARK_DIR / f"run_{cores}W"
    print(f"\n[BENCHMARK] Setting up case for {cores} cores in {run_dir.name}...")
    
    # Clean up existing run folder
    if run_dir.exists():
        shutil.rmtree(run_dir)
    run_dir.mkdir(parents=True, exist_ok=True)
    
    # Copy initial, constant, system
    shutil.copytree(CASE_DIR / "initial", run_dir / "initial")
    shutil.copytree(CASE_DIR / "constant", run_dir / "constant")
    shutil.copytree(CASE_DIR / "system", run_dir / "system")
    
    # Prepare "0" directory
    shutil.copytree(run_dir / "initial", run_dir / "0")
    
    # 1. Modify controlDict: set endTime=9.0e-5, writeInterval=1.0 (to avoid I/O overhead skewing benchmark)
    control_dict_path = run_dir / "system" / "controlDict"
    with open(control_dict_path, "r") as f:
        content = f.read()
    content = re.sub(r"endTime\s+[\d\.e\-+]+;", f"endTime {END_TIME};", content)
    content = re.sub(r"writeInterval\s+[\d\.e\-+]+;", "writeInterval 1.0;", content)
    with open(control_dict_path, "w") as f:
        f.write(content)
        
    # 2. Modify decomposeParDict: set numberOfSubdomains
    decomp_dict_path = run_dir / "system" / "decomposeParDict"
    with open(decomp_dict_path, "r") as f:
        content = f.read()
    content = re.sub(r"numberOfSubdomains\s+\d+;", f"numberOfSubdomains {cores};", content)
    with open(decomp_dict_path, "w") as f:
        f.write(content)
        
    return run_dir

def run_cmd(cmd: list, cwd: Path, log_file=None) -> float:
    start_time = time.time()
    if log_file:
        with open(log_file, "w") as lf:
            subprocess.run(cmd, cwd=str(cwd), stdout=lf, stderr=subprocess.STDOUT, check=True)
    else:
        subprocess.run(cmd, cwd=str(cwd), stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, check=True)
    return time.time() - start_time

def benchmark_cores(cores: int) -> dict:
    run_dir = setup_case(cores)
    
    print(f"[BENCHMARK] [{cores} Cores] Running blockMesh...")
    run_cmd(["blockMesh"], run_dir, run_dir / "log.blockMesh")
    
    print(f"[BENCHMARK] [{cores} Cores] Running setFields...")
    run_cmd(["setFields"], run_dir, run_dir / "log.setFields")
    
    print(f"[BENCHMARK] [{cores} Cores] Decomposing domain...")
    run_cmd(["decomposePar"], run_dir, run_dir / "log.decomposePar")
    
    print(f"[BENCHMARK] [{cores} Cores] Running laserbeamFoam solver...")
    solver_cmd = ["mpirun", "--oversubscribe", "-np", str(cores), "laserbeamFoam", "-parallel"]
    
    # Measure execution time
    solver_time = run_cmd(solver_cmd, run_dir, run_dir / "log.laserbeamFoam")
    print(f"[BENCHMARK] [{cores} Cores] Solver execution completed in {solver_time:.2f} seconds.")
    
    return {
        "cores": cores,
        "solver_time_s": solver_time,
        "run_dir": run_dir
    }

def main():
    results = []
    
    # Run benchmark for each core configuration
    for cores in CORE_CONFIGS:
        try:
            res = benchmark_cores(cores)
            results.append(res)
        except Exception as e:
            print(f"[ERROR] Benchmark failed for {cores} cores: {e}")
            
    if not results:
        print("[ERROR] All benchmarks failed.")
        return
        
    print("\n" + "="*50)
    print("                BENCHMARK RESULTS")
    print("="*50)
    print(f"{'Cores':<10}{'Solver Exec Time (s)':<25}{'Speedup vs 4-core':<20}")
    print("-"*50)
    
    base_time = results[0]["solver_time_s"]
    for r in results:
        speedup = base_time / r["solver_time_s"]
        print(f"{r['cores']:<10}{r['solver_time_s']:<25.2f}{speedup:<20.2f}")
    print("="*50)
    
    # Plot results
    cores_list = [r["cores"] for r in results]
    times_list = [r["solver_time_s"] for r in results]
    
    fig, ax1 = plt.subplots(figsize=(8, 5))
    
    # Plot Execution Time
    color = '#1f77b4'
    ax1.set_xlabel('Number of Cores', fontweight='bold', fontsize=11)
    ax1.set_ylabel('Solver Execution Time (s)', color=color, fontweight='bold', fontsize=11)
    bars = ax1.bar([str(c) for c in cores_list], times_list, color=color, alpha=0.7, width=0.4, label='Exec Time')
    ax1.tick_params(axis='y', labelcolor=color)
    ax1.grid(True, linestyle='--', alpha=0.3)
    
    # Add values on top of bars
    for bar in bars:
        height = bar.get_height()
        ax1.annotate(f'{height:.1f}s',
                    xy=(bar.get_x() + bar.get_width() / 2, height),
                    xytext=(0, 3),  # 3 points vertical offset
                    textcoords="offset points",
                    ha='center', va='bottom', fontweight='bold')
                    
    # Plot Speedup/Efficiency on secondary Y axis
    ax2 = ax1.twinx()  
    color = '#d62728'
    ax2.set_ylabel('Speedup Factor', color=color, fontweight='bold', fontsize=11)
    speedups = [base_time / t for t in times_list]
    ax2.plot([str(c) for c in cores_list], speedups, color=color, marker='o', linewidth=2, markersize=8, label='Speedup')
    ax2.tick_params(axis='y', labelcolor=color)
    
    # Theoretical linear speedup line
    ideal_speedup = [c / cores_list[0] for c in cores_list]
    ax2.plot([str(c) for c in cores_list], ideal_speedup, color='gray', linestyle=':', label='Ideal Speedup')
    
    plt.title('OpenFOAM Parallel Scaling Benchmark (0.09 ms SLM Simulation)', fontweight='bold', fontsize=12)
    fig.tight_layout()
    
    # Save the plot
    plot_path = BENCHMARK_DIR / "benchmark_scaling.png"
    plt.savefig(plot_path, dpi=200)
    plt.close()
    
    # Write summary report
    summary_path = BENCHMARK_DIR / "benchmark_summary.md"
    with open(summary_path, "w") as sf:
        sf.write("# Parallel Scaling Benchmark Summary\n\n")
        sf.write(f"This scaling test evaluates the performance of `laserbeamFoam` solver on the autogenous SLM 316L single-track case (Abdali parameters) using **4**, **8**, and **16 cores** for a simulation duration of **{END_TIME*1e3:.2f} ms**.\n\n")
        
        sf.write("## Benchmark Configuration\n")
        sf.write("- **Domain size**: $0.3\\text{ mm} \\times 0.3\\text{ mm} \\times 1.2\\text{ mm}$\n")
        sf.write("- **Mesh grid**: $5\\ \\mathrm{\\mu m}$ fine mesh (with 2-level AMR)\n")
        sf.write(f"- **Simulated duration**: {END_TIME*1e3:.2f} ms\n")
        sf.write("- **Write interval**: Disabled during benchmark to isolate CPU compute performance from I/O bottlenecking.\n\n")
        
        sf.write("## Performance Results Table\n\n")
        sf.write("| Cores | Solver Execution Time (s) | Speedup Factor | Parallel Efficiency |\n")
        sf.write("| :---: | :---: | :---: | :---: |\n")
        for r in results:
            speedup = base_time / r["solver_time_s"]
            efficiency = (speedup / (r["cores"] / cores_list[0])) * 100.0
            sf.write(f"| {r['cores']} | {r['solver_time_s']:.2f}s | {speedup:.2f}x | {efficiency:.1f}% |\n")
            
        sf.write("\n## Key Findings\n")
        fastest_config = min(results, key=lambda x: x["solver_time_s"])
        sf.write(f"- **Optimal Configuration**: Running on **{fastest_config['cores']} cores** reached 0.09 ms the fastest, taking **{fastest_config['solver_time_s']:.2f} seconds**.\n")
        
        # Analyze scaling behavior
        if len(results) > 1:
            ratio_16_8 = results[1]["solver_time_s"] / results[2]["solver_time_s"] if len(results) > 2 else 0
            if ratio_16_8 < 1.1 and ratio_16_8 > 0:
                sf.write("- **Diminishing Returns**: The scaling performance shows significant diminishing returns between 8 and 16 cores. This is typical for small grid domains where inter-processor MPI communication latency outweighs the computational benefits of additional cores.\n")
            else:
                sf.write("- **Scaling Trend**: Increasing core count demonstrates acceleration. However, communication overhead on smaller mesh sizes starts to limit efficiency as core counts increase.\n")
                
        sf.write(f"\nThe scaling visualization has been plotted and saved as `benchmark_scaling.png` inside the `benchmark` directory.\n")
        
    print(f"\n[BENCHMARK] Scaling plot saved to: {plot_path}")
    print(f"[BENCHMARK] Summary report written to: {summary_path}")
    
    # Cleanup big pre-processing and mesh files in subfolders to save space, leaving only log files
    for r in results:
        for p in r["run_dir"].iterdir():
            if p.is_dir() and p.name not in ["initial"]:
                shutil.rmtree(p)
            elif p.is_file() and not p.name.startswith("log."):
                p.unlink()

if __name__ == "__main__":
    main()
