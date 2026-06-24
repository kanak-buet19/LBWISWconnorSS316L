# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## What this workspace is

Tutorial and calibration environment for 316L stainless steel single-track melt-pool studies using `laserbeamFoam`. Houses the master case template, experimental validation runs (Hofmann 2026 CW, Wang Ti-6Al-4V), and a Bayesian parameter calibration framework. All cases follow OpenFOAM directory structure (`initial/`, `constant/`, `system/`, `scripts/`).

For build, solver architecture, and library details, see `../../../CLAUDE.md`.

## Quick commands

```bash
# Run a single case (from case directory)
./Allrun                          # blockMesh → setFields → decomposePar → laserbeamFoam
./Allclean_long                   # remove generated outputs (time dirs, VTK/, post-processing)

# Reconstruct and convert to VTK
./foamVTK.sh all                  # reconstruct all written times, convert to VTK

# Calibration (316L CW melt-pool Bayesian optimization)
cd calibrate_CW_case
./Allrun                          # local run (sequential candidates)
sbatch job_calibration.sh         # HPC (Slurm + Apptainer)
CALIB_TOTAL_CORES=16 CALIB_CORES_PER_SIM=4 ./Allrun   # override parallelism

# Post-calibration
python scripts/plots.py           # generate ranked_objective, params_vs_objective, etc.
watch -n 10 cat results/status.txt    # live progress (HPC runs)

# Analyze one VTU file directly
python template/scripts/analyze_meltpool_vtu.py --vtk-file VTK/.../internal.vtu
```

## Directory structure

| Path | Purpose |
|------|---------|
| `template/` | Master OpenFOAM case; source for copying to new studies |
| `calibrate_CW_case/` | Bayesian optimization workflow (see its CLAUDE.md for full detail) |
| `validation/` | Standalone validation runs; `cases/` holds per-experiment sub-dirs, `cases.json` controls which run |
| `validation/template_case` | Base SS316L template for CW Hofmann validation |
| `validation/template_Ti64` | Base Ti-6Al-4V template for Wang validation |
| `backups/` | Ad-hoc snapshots; not tracked by calibration |
| `reference_paper/`, `docs/`, `exp_figure/` | Literature, notes, experimental images |

## Case structure (template and all derived cases)

```
case_name/
  Allrun, Allclean_long, foamVTK.sh  ← run scripts
  initial/               ← initial fields (T, U, p_rgh, alpha.metal, Laser_boundary)
  constant/
    transportProperties  ← fluid props (rho, nu, cp, kappa, Marangoni, LeeCoeff, etc.)
    LaserProperties      ← laser position/power tables, Fresnel absorption model
    fvSolution, fvSchemes
    dynamicMeshDict      ← AMR config (criteria, dRefine, dUnrefine)
  system/
    blockMeshDict        ← mesh generation (domain size, resolution, boundaries)
    controlDict          ← time stepping, output, max Courant
    decomposeParDict     ← MPI domain decomposition
  scripts/
    analyze_meltpool_vtu.py    ← per-timestep VTU → melt-pool geometry CSV
    parse_simulation_log.py    ← tees laserbeamFoam stdout
  post-processing-data/  ← generated: per-timestep vtu_meltpool_geometry.csv
```

## Calibration targets (Hofmann 2026 CW, SS316L)

Three experimental cases in `calibrate_CW_case/exp_cases.csv`:

| # | Regime | Power (W) | Speed (mm/s) | Beam r (µm) | Exp width (µm) | Exp depth (µm) |
|---|--------|-----------|--------------|-------------|----------------|----------------|
| 1 | conduction | 200 | 1200 | 26 | 94.7 | 73.8 |
| 2 | transition | 200 | 900 | 26 | 114.6 | 94.0 |
| 3 | keyhole | 250 | 600 | 26 | 142.0 | 204.6 |

Objective = mean over cases of `mean(|depth_err|, |width_err|, |AR_err|)`. The AR term prevents compensating depth/width errors.

## Calibration two-phase mesh strategy

Phase 1 (coarse): 40 µm base + maxRefinement=1 → 20 µm fine. Runs up to `nSamples=64`. If best objective < `transitionObjective=0.20`, top-K candidates seed Phase 2.

Phase 2 (fine): 40 µm base + maxRefinement=2 → 10 µm fine. Runs up to `nSamples=160`. Surrogate is warm-started from Phase 1 + LHS seeds.

## Active optimizer knobs (3 tuned params)

| Param | Range | Mechanism |
|---|---|---|
| `elec_resistivity` | 8e-7–1.5e-6 Ω·m | Controls laser absorption via ITO model |
| `LeeCoeff` | 3e5–1e6 1/s | Volumetric evaporation cooling strength |
| `dSigmadT_norm` | maps to -7e-4 to -4.9e-4 N/m/K | Marangoni lever |

All other SS316L thermophysical values are fixed baselines (see `calibrate_CW_case/CLAUDE.md` for full list).

## Naming conventions

Case folders: material/source + power + speed + date or reference, e.g. `cw_scantrack_200W_900mms_2026_06_14`.

Python scripts: use `pathlib`, small helper functions, explicit config keys. Test on existing sample output before launching full calibration.

## Testing & development

1. Run narrowest case first; check `log.blockMesh`, `log.setFields`, `log.laserbeamFoam`.
2. For script edits: test on existing sample output before full calibration run.
3. For case edits: inspect `transportProperties` and `dynamicMeshDict` first.
4. Minimum regression: mesh generation, solver start, VTK conversion, melt-pool CSV produced.

## Commit guidelines

Conventional Commits style, e.g.:
- `feat: add Lee volumetric evaporative cooling model`
- `fix(calibrate): add aspect-ratio penalty to objective function`

Scope commits to one case, script, or model change. Include plots or CSV summaries when results change. Exclude `processor*/`, `VTK/`, time directories, and logs.

See `AGENTS.md` for extended style guidelines. See `calibrate_CW_case/CLAUDE.md` for full calibration architecture, resume logic, watchdog config, and output file descriptions.
