# Repository Guidelines

## Project Structure & Module Organization

This directory is a `laserbeamFoam` tutorial and calibration workspace for 316L single-track melt-pool studies. Reusable case inputs live in `template/`, with OpenFOAM dictionaries under `system/`, material and laser setup under `constant/`, starting fields under `initial/`, and post-processing helpers under `scripts/`. Named study cases such as `cw_scantrack_200W_900mms_2026_06_14/`, `zhang_N01_260W_520mms_BP/`, and `mesh_convergence/mesh_8um/` follow the same layout. Calibration orchestration is in `calibrate_CW_case/`, especially `calibration_config.json` and `scripts/`. Reference data and notes are kept in `reference_paper/`, `docs/`, and `exp_figure/`.

## Build, Test, and Development Commands

Source OpenFOAM before running cases, for example `source /openfoam/bash.rc`.

- `cd template && ./Allrun`: prepare `0/`, run `blockMesh`, `setFields`, then execute `laserbeamFoam` via `Allrun_long`.
- `cd template && ./Allclean_long`: remove generated case outputs such as `0/`, `VTKs/`, and post-processing data.
- `cd template && ./foamVTK.sh all`: reconstruct and convert all written times to VTK.
- `cd calibrate_CW_case && CALIB_TOTAL_CORES=16 CALIB_CORES_PER_SIM=4 ./Allrun`: launch the local Bayesian calibration swarm.
- `python scripts/analyze_meltpool_vtu.py --vtk-file VTK/.../internal.vtu`: analyze one converted VTU file.

## Coding Style & Naming Conventions

Keep shell scripts POSIX/Bash-oriented, with existing `Allrun`, `Allclean`, and `foamVTK.sh` patterns. Use four-space indentation for Python and prefer `pathlib`, small helper functions, and explicit configuration keys. Preserve OpenFOAM dictionary formatting and names such as `LaserProperties`, `transportProperties`, `fvSchemes`, and `fvSolution`. Name new case folders with material or source plus operating point, for example `zhang_N02_260W_870mms_BP`.

## Testing Guidelines

For case edits, run the narrowest relevant case first and inspect `log.blockMesh`, `log.setFields`, and `log.laserbeamFoam`. For script edits, run the specific Python entry point on an existing sample output before launching a full calibration. Treat successful mesh generation, solver start, VTK conversion, and melt-pool CSV production as the minimum regression checks.

## Commit & Pull Request Guidelines

Recent history uses concise imperative subjects, often Conventional Commit style, such as `feat: add Lee volumetric evaporative cooling model` or `fix(calibrate): add aspect-ratio penalty to objective function`. Keep commits scoped to one case, script, or model change. Pull requests should describe changed parameters, affected cases, commands run, and key outputs; include plots or CSV summaries when results change. Do not include bulky generated directories (`processor*`, `VTK/`, `VTKs/`, time directories, or logs) unless explicitly needed for review.
