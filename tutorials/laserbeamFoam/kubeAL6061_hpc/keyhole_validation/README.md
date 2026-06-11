# Keyhole Validation

This folder is a standalone workflow for reproducing the Kube Al6061
paper-style simulation setup, then exploring and optimizing the process
parameters that control stable keyhole/melt-pool dimensions.

The active template is:

```text
reference_template/
```

The runner creates cases under:

```text
runs/
```

## Paper-Style Setup

The reference template starts from the paper simulation setup, but the active
calibration run is longer because the available aspect-ratio experiment file is
a stable-window segment, not the initial heating transient.

- laser power: `350 W`
- beam diameter: `80 um`
- beam radius: `40 um`
- stationary beam
- full calibration laser-on time: `3.0 ms`
- full calibration simulation end time: `3.0 ms`
- smoke-test simulation end time: `1.3 ms`
- base mesh: `20 um`
- dynamic mesh: `maxRefinement 1`
- effective refined cell size: about `10 um`

## Parameters For Exploration

The first sensitivity/optimization variables are:

- laser power: `340-355 W`
- laser radius: `40.5-43 um`
- electron number density: `1e29-4e29`
- electrical resistivity: `3e-7-6.5e-7 ohm m`

## Target

The calibration target is the stable-window melt-pool aspect-ratio trace:

- melt-pool aspect ratio: `depth / width`
- comparison window: best gated `2.0 ms` simulation segment found by sliding
  across the available time history
- experiment: `../kube_exp_data/csv/1ms_meltpool_ar_exp.csv`
- score terms: AR RMSE, AR bias, AR slope mismatch, window coverage,
  dimension-gate penalty, and hard-runaway penalty

Absolute width/depth are not forced to exactly match every instant, but they
gate the aspect-ratio target. A case cannot be marked successful just because
the aspect ratio is stable if the melt-pool width/depth are physically
unreasonable. Keyhole-depth clipping is treated as physical overshoot/censored
data, so AR data after clipping is not trusted for success.

Whenever a case has any written time inside the width/depth diagnostic box, the
workflow saves a manifest and section images for later postprocessing:

```text
runs/<case_id>/post-processing-data/target_meltpool_snapshots/
```

The manifest records every hit time. PNG copying is capped by
`max_target_snapshot_pngs_per_case` in the config to avoid runaway storage.

For the experimental aspect-ratio comparison, the workflow scans possible
simulation windows, uses the best gated window start as `t = 0`, then plots up
to `2.0 ms` of simulation against
`../kube_exp_data/csv/1ms_meltpool_ar_exp.csv`.

## Run On HPC

Submit:

```bash
sbatch keyhole_calibrate_job.sh
```

The job requests 64 MPI tasks for 48 hours, runs 8 cases at a time, and uses 8
MPI ranks per case.

## Adaptive Optimizer

After the initial sensitivity cases, the optimizer uses an ensemble surrogate
model when `scikit-learn` is available:

1. Train `ExtraTreesRegressor` models from completed CFD cases.
2. Predict score, aspect-ratio error metrics, and comparison-window duration for many
   virtual candidates inside the parameter bounds.
3. Rank candidates by predicted AR trace error, dimension-gate feasibility,
   predicted score, uncertainty, and diversity.
4. Select the next batch of 8 diverse candidates.
5. Repeat until a stable target case is found or the iteration limit is reached.

If `scikit-learn` is unavailable, the script falls back to a local/random
trust-region search around the best completed case.

## Early Stop

Each case checks the partial `post-processing-data/vtu_meltpool_geometry.csv`
after every written-time analysis.

The case can stop early as:

- `stop_success`: AR trace matched while width/depth stayed physically reasonable
- `stop_fail`: case is clearly too hot/runaway or too cold
- `continue`: case still might recover

The decision is saved inside each run directory:

```text
early_stop_status.json
```

These decisions are also copied into `summary.csv` and shown in the live
dashboard. Early-stop thresholds are controlled by the `early_stop` block in
`config.json`.

## Run Locally

From this folder, after sourcing OpenFOAM:

```bash
export PYTHON=python3
export MAX_PARALLEL_CASES=1
export RANKS_PER_CASE=8
python3 scripts/keyhole_search.py --mode sensitivity --run
```

For dry generation only:

```bash
python3 scripts/keyhole_search.py --mode sensitivity --prepare-only --force
```

## Smoke Adaptive Pipeline

To verify the adaptive workflow with a coarse static mesh:

```bash
./Allrun
```

This uses:

- `config_smoke.json`
- `40 um` base mesh
- no adaptive mesh refinement
- 6 sensitivity exploration cases
- 2 proposed optimization cases per iteration
- 2 optimization iterations by default
- early-stop checks using the smoke config
- 2 MPI ranks per case
- 2 cases at a time

Outputs:

```text
smoke_summary.csv
smoke_plots/live_progress_dashboard.png
smoke_plots/best_aspect_ratio_vs_exp.png
```

Clean the smoke test with:

```bash
./Allclean
```

## Outputs

- `runs/<case_id>/post-processing-data/vtu_meltpool_geometry.csv`
- `runs/<case_id>/post-processing-data/target_meltpool_snapshots/`
- `summary.csv`
- `plots/live_progress_dashboard.png`
- `plots/sensitivity_parallel_coordinates.png`
- `plots/best_aspect_ratio_vs_exp.png`
- `best_case_summary.json`

## Live Progress Dashboard

The workflow continuously refreshes one combined progress image:

```text
plots/live_progress_dashboard.png
```

It summarizes:

- status counts for all cases
- completed case scores and best-so-far score
- stable-window width/depth target map
- best case parameters and metrics
- incomplete or currently-running cases
- top-ranked cases
- next/planned cases

To manually refresh the dashboard while cases are running:

```bash
python3 scripts/keyhole_search.py --mode score
```

This only scans existing run folders, updates `summary.csv`, and rewrites the
plots. It does not start new simulations.
