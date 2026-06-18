# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## What this repo does

Bayesian-optimization calibration of 316L SS thermophysical parameters for LaserbeamFoam CFD simulations. Targets are CW single-track melt-pool depth/width from Hofmann (2026) experiments. Calibrated params feed subsequent pulsed-laser gravity-variation sims (Riffel et al. 2026).

## Commands

```bash
# Local run (OpenFOAM must be sourced in shell)
./Allrun

# HPC (Slurm + Apptainer)
sbatch job_calibration.sh

# Override parallelism (local or HPC)
CALIB_TOTAL_CORES=16 CALIB_CORES_PER_SIM=4 ./Allrun

# Wipe all generated data (keeps template_case/, scripts/, config)
./Allclean

# Dry-run sanity check on HPC (skips swarm launch)
DRY_RUN=1 sbatch job_calibration.sh
```

Python deps required: `matplotlib numpy optuna pandas pyvista scipy tqdm vtk`

## Architecture

```
calibration_config.json          ← single source of truth: all knobs
scripts/
  run_calibration.py             ← main: Calibrator orchestrator + Job model
  optimizer.py                   ← BOOptimizer wrapping Optuna TPE + LHS seed
  caselib.py                     ← builds OF cases from template via regex patching
  evaluate.py                    ← melt-pool stability scorer (reads per-sim CSV)
  plots.py                       ← post-run visualization
template_case/                   ← OpenFOAM case template (blockMesh, controlDict, etc.)
  constant/transportProperties   ← patched per candidate (all 12 calibration params)
  system/blockMeshDict           ← patched per case (domain size from geometry config)
  Allrun_long                    ← per-sim runner: blockMesh → setFields → decomposePar → laserbeamFoam
  scripts/analyze_meltpool_vtu.py  ← per-timestep VTU analysis → CSV
  scripts/parse_simulation_log.py  ← tees laserbeamFoam stdout
runs/cand_XX/case_name/          ← generated: one dir per candidate×case
results/                         ← generated: summary.csv, summary.json, best_params.json, plots/
```

## Calibration loop flow

1. `Calibrator.__init__` → loads config, replays prior results into surrogate (`_load_prior`)
2. `BOOptimizer.seed_lhs` → enqueues 16 LHS points as initial design
3. Main loop: `ask()` → `caselib.build_case()` → `Job.start()` (spawns `Allrun_long` subprocess)
4. `Allrun_long` runs OpenFOAM, calls `analyze_meltpool_vtu.py` after each written timestep → appends to `post-processing-data/vtu_meltpool_geometry.csv`
5. `evaluate.py` reads that CSV, finds the longest stable tail (variation < `stabilityTol`), computes `case_error = mean(|depth_err|, |width_err|, |ar_err|)` where `ar_err` = aspect-ratio error (depth/width vs exp) — prevents compensating errors (e.g. depth+10%/width−10%) from scoring well
6. Early abort: if pool has stabilized AND error > 50%, sim is killed immediately
7. `BOOptimizer.tell()` → feeds result back; TPE proposes smarter next candidates
8. On completion: `best_params.json` written, plots generated

## Resume behaviour

Candidates persist via `runs/cand_XX/result.json`. On resubmit, completed candidates are replayed into the surrogate (no re-run); in-flight candidates are re-queued. Safe to resubmit a Slurm job after timeout.

## Key config fields

| Field | Where | Purpose |
|---|---|---|
| `optimizer.nSamples` | config | Total BO budget (default 48) |
| `optimizer.initSamples` | config | LHS exploration phase size (default 16) |
| `execution.coresPerSim` / `totalCores` | config | Parallelism; override with env vars |
| `earlyAbort.errorThreshold` | config | Kill stabilized-but-bad sims at this error fraction |
| `earlyAbort.stabilityTol` | config | Pool is "stable" when (max-min)/mean < this over tail |
| `geometry.*` | config | Domain sizing; drives blockMeshDict patching |
| `control.*` | config | endTime, deltaT, maxCo, writeInterval |

## The 12 calibration parameters

All live in `constant/transportProperties` and are patched by `caselib.patch_transportProperties`:

- `elec_resistivity` — controls laser absorption (ITO model)
- `cp_scale`, `kappa_scale` — multiplicative on solid cp/kappa tables
- `cp_liquid_slope`, `kappa_liquid_slope` — additive dX/dT in liquid phase (T ≥ Tliquidus=1723K, capped at Tvap=3122K)
- `rho`, `nu`, `LatentHeat`, `LatentHeatVap`, `sigma`, `Marangoni_Constant`, `beta_r`

`_sub_entry` patches scalar fields by regex; `_scale_table` / `_apply_liquid_slope` patch table entries. If you add a new parameter, you must add a patcher call in `patch_transportProperties`.

## Environment variables

| Var | Default | Purpose |
|---|---|---|
| `CALIB_TOTAL_CORES` | config value (64) | Total cores available to swarm |
| `CALIB_CORES_PER_SIM` | config value (8) | Cores per laserbeamFoam instance |
| `CALIB_USE_APPTAINER` | `0` | Set to `1` on HPC to run sims inside Apptainer image |
| `OF2506_IMAGE` | `~/openfoam-dev_2506.sif` | Apptainer SIF path |
| `OF2506_USER` | `x-rkanak1` | Username injected into container |
| `PYTHON` | `/home/kanak/.venv/venv314/bin/python` | Python used for orchestrator + analysis |
| `FOAM_SIGFPE` | unset | Set to `0` to suppress FP exceptions in OpenFOAM |

## Watchdog logic

Progress-based (not wall-clock). A sim is killed only if `log.laserbeamFoam` stops advancing:
- No solver output at all for `perSimStartupSec` (3600s) → killed as "no solver output"
- Solver Time frozen for `perSimStallSec` (1800s) → killed as "solver Time stalled"
- Slow-but-progressing sims are never killed.

## Disk cleanup

`cleanup.enabled=true` keeps only the best-so-far candidate's full data (VTK + processor*/). Non-best candidates are stripped after finishing. `strip_heavy()` deletes VTK/, processor*/, and reconstructed time dirs; keeps CSV, logs, `case_build.json`, `result.json`.

## Output files (results/)

| File | Content |
|---|---|
| `best_params.json` | Winning parameter set + per-case errors |
| `summary.json` | All candidates with full series data |
| `summary.csv` | Flat table suitable for quick inspection |
| `status.json` | Live progress (updated every poll interval) |
| `plots/ranked_objective.png` | Candidates ranked by objective |
| `plots/params_vs_objective.png` | Each param vs objective scatter |
| `plots/best_sim_vs_exp.png` | Best candidate sim vs exp bars |
| `plots/best_convergence.png` | Best candidate depth/width stabilization |
