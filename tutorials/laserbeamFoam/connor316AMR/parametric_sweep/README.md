# CW Parametric Sweep

Four-case sweep using `../calibrate_CW_case/template_case` and the baseline
values from `../calibrate_CW_case/calibration_config.json`.

The sweep keeps the Hofmann 200 W, 900 mm/s, 25 um radius case fixed and only
changes one solid-phase table scale at a time:

- `table_cp_solid_scale`: `1.10`, `1.20`
- `table_kappa_solid_scale`: `1.10`, `1.20`

Run after sourcing OpenFOAM:

```bash
source /openfoam/bash.rc
./Allrun
```

Useful commands:

```bash
./Allrun --generate-only
./Allrun --skip-existing
./Allrun --only solid_cp_p10 solid_kappa_p20
./Allclean
```

Outputs are written to `results/summary.csv`, `results/summary.json`, and
per-case logs under `results/run_logs/`.
