# Connor AMR Melt-Pool Calibration

Goal: calibrate `laserRadius` and metal `elec_resistivity` against the Connor/Riffel experimental melt-pool size.

Main settings live in `calibration_config.json`:

- experimental target width/depth
- first-pulse time window
- simulation end time
- write interval
- `laserRadius` candidates
- `elec_resistivity` candidates
- optional z-geometry shrink settings
- optional coarse test mode
- OpenFOAM run command

Method:

1. Copy the baseline `connor316AMR` case into `runs/<case_id>`.
2. Patch `constant/LaserProperties`, `constant/transportProperties`, and `system/controlDict`.
3. If `shrink_z_geometry.enabled` is true, shrink only the copied case's z-domain around the first-pulse laser path plus offsets.
4. Run only the first pulse window.
5. Use the existing `scripts/analyze_meltpool_vtu.py`.
6. Score each simulation using the peak melt-pool dimensions during the first pulse:

   `score = sqrt(width_error_pct^2 + depth_error_pct^2)`
7. If adaptive mode is enabled, fit a surrogate model from completed cases, predict the next best untested radius/resistivity pair, run it, and repeat until either:

   - width and depth are both within `within_percent`
   - `adaptive.max_new_cases` is reached
   - no new candidate is available

Current z shrink:

- first-pulse simulation end: `0.008 s`
- write interval: `0.00025 s`, giving about 33 output samples over `0.008 s`
- laser z travel over that time: about `32 um`
- default z offsets: `700 um` before and after
- resulting z length: about `1432 um`
- resulting base z cells: about `29` instead of `70`

Run:

```bash
./Allrun
```

`Allrun` uses `/home/x-rkanak1/.conda/envs/isw_env/bin/python3` by default on Anvil and forwards extra arguments.

Useful fast test:

```bash
./Allrun --max-cases 1 --no-adaptive
```

Clean generated calibration files:

```bash
./Allclean
```

Test entire pipeline cheaply with coarse mesh/settings:

```bash
./Allclean
./Allrun --coarse-test
```

Coarse test mode currently uses:

- same time window as full calibration
- same initial parameter grid as full calibration
- same adaptive optimizer settings as full calibration
- block mesh x/y cells: `12 x 8`
- coarser z base cell size: `150 um`
- write interval: `0.001 s`
- laser rays: `10 x 48 = 480`
- MPI subdomains: `4`
- parallel cases: `4`
- adaptive proposals per iteration: `4`

Override parallel case count:

```bash
./Allrun --coarse-test --jobs 2
./Allrun --coarse-test --jobs 4
```

With `--jobs 4` and `numberOfSubdomains: 4`, total MPI ranks can reach `16`.
Each parallel case writes its top-level run output to:

```text
runs/<case_id>/log.calibrationRunner
```

The driver sources OpenFOAM v2506 non-interactively by default:

```bash
apptainer exec --cleanenv --env USER=${OF2506_USER:-x-rkanak1} ${OF2506_IMAGE:-$HOME/openfoam-dev_2506.sif} bash -lc 'source /openfoam/bash.rc && export FOAM_SIGFPE=0 && ./Allrun'
```

Override if needed:

```bash
/home/x-rkanak1/.conda/envs/isw_env/bin/python3 calibrate_meltpool.py --run --run-command "apptainer exec --cleanenv --env USER=${OF2506_USER:-x-rkanak1} ${OF2506_IMAGE:-$HOME/openfoam-dev_2506.sif} bash -lc 'source /openfoam/bash.rc && export FOAM_SIGFPE=0 && ./Allrun'"
```

Prepare cases without running:

```bash
/home/x-rkanak1/.conda/envs/isw_env/bin/python3 calibrate_meltpool.py
```

Run only the initial sensitivity grid, with no adaptive proposals:

```bash
/home/x-rkanak1/.conda/envs/isw_env/bin/python3 calibrate_meltpool.py --run --no-adaptive
```

Use another config:

```bash
/home/x-rkanak1/.conda/envs/isw_env/bin/python3 calibrate_meltpool.py --config my_config.json --run
```

Outputs:

- `runs/<case_id>/`: full copied OpenFOAM case and generated data
- `summary.csv`: ranked calibration results
- `best_case.txt`: best parameter set found
- `plots/score_map.png`: parameter map colored by total score
- `plots/error_maps.png`: width/depth error maps
- `plots/best_vs_experiment.png`: best simulation vs experiment bar plot
- `plots/top_cases.png`: ranked top-case plot

Regenerate plots from an existing `summary.csv`:

```bash
./PlotResults
```

Run fine calibration around the coarse-best region:

```bash
./Allclean
./Allrun_fine_zone
```

`fine_zone_config.json` uses full mesh/rays/time settings with:

- radius grid: `335, 345, 355, 365, 375 um`
- resistivity grid: `2.1e-6, 2.25e-6, 2.4e-6, 2.55e-6`
- adaptive bounds: `320-390 um`, `1.9e-6-2.7e-6`
- parallel cases: `2`
