# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## What this is

Standalone validation suite for `laserbeamFoam`. Cases are generated from `cases.json` by copying `template_case/` and patching laser/mesh parameters. Results are compared to published melt-pool depth and width measurements.

For solver architecture, build, and library details see `../../../CLAUDE.md`.

## Running cases

```bash
# Interactive menu (select one case or all)
./Allrun

# Non-interactively
./Allrun --case cw_scantrack_200W_900mms_2026_06_14
./Allrun --all                           # include disabled cases too

# Setup only (no run)
./Allrun --setup-only
./Allrun --setup-only --case wang_Ti6Al4V_382W_475mms_r50um

# Rebuild case from template (overwrites case dir, keep post-processing-data)
python3 scripts/setup_cases.py --rebuild --case <name>

# Post-process a specific VTU manually
python3 template_case/scripts/analyze_meltpool_vtu.py \
    --case runs/<name> \
    --vtk-file runs/<name>/VTK/<time>/internal.vtu
```

`Allrun` calls `setup_cases.py` to prepare/sync the case directory, then runs `Allclean_long` + `Allrun` inside each case dir.

## Cases.json — adding or editing cases

`cases.json` is the single source of truth. `setup_cases.py` reads it and templates every OpenFOAM file. Edit `cases.json`; never hand-edit generated case files for things covered by the schema.

### Top-level keys

| Key | Effect |
|---|---|
| `template` | Default template dir (relative to `validation/`) |
| `write_interval_s` | Default write interval, case can override |
| `base_mesh_size_um` | Default base cell size; case can override |
| `amr_levels` | Default AMR refinement levels; case can override |
| `run_enabled_only` | If true, `./Allrun` skips cases with `"enabled": false` |

### Per-case keys

| Key | Required | Notes |
|---|---|---|
| `name` | yes | Directory name under `runs/` |
| `material` | yes | Must match a dir under `materials/` (`SS316L` or `Ti64`) |
| `enabled` | no | Defaults to `true` |
| `power_W` | yes | Laser power |
| `scan_speed_mm_s` | yes | 0 = stationary spot |
| `laser_radius_m` | yes | Beam radius (1/e²) |
| `laser_start_z_m` | yes | Z position at t=0 |
| `laser_end_z_m` | — | Computes `end_time_s` from speed; mutually exclusive with `end_time_s` |
| `end_time_s` | — | Explicit simulation end time |
| `laser_off_time_s` | no | When laser turns off; defaults to `end_time_s` |
| `target_width_um` | no | Expected melt-pool width for comparison plots |
| `target_depth_um` | no | Expected depth; also drives domain Y-sizing |
| `depth_metric` | no | `meltpool` (default) or `keyhole` — selects contour for comparison |
| `base_mesh_size_um` | no | Overrides suite-level default |
| `amr_levels` | no | Overrides suite-level default |
| `domain_x_m` | no | Override X extent |
| `domain_z_m` | no | Override Z extent |
| `domain_depth_m` | no | Override Y (substrate) extent, else auto from `target_depth_um × 1.5` |
| `electric_resistivity_values` | no | Single value or list → generates sweep sub-cases named `<name>_rho<value>` |
| `exp_tag` | no | Prefix for `exp_<tag>_<power>W_<speed>_summary.csv` lookup |
| `exp_timeseries_csv` | no | Path (relative to `validation/`) to `(t_ms, keyhole_depth_um)` CSV for the plot and live depth-error comparison |
| `write_interval_s` | no | Per-case override |
| `template` | no | Per-case template dir override |

### Domain sizing

- Y (depth): auto = gas layer (0.128 mm) + max(1.5 × `target_depth_um`, 224 µm). Override with `domain_depth_m`.
- X, Z: template defaults unless `domain_x_m` / `domain_z_m` set.
- All extents rounded to nearest `base_mesh_size_um` cell.

## Material library

`materials/<MATERIAL>/transportProperties` is the canonical property file. `setup_cases.py` copies it into each case's `constant/` at setup time, overwriting the template's copy.

- `SS316L` — 316L stainless steel (Hofmann 2026 targets)
- `Ti64` — Ti-6Al-4V (Wang / Simonds / Tao targets)

To change a material property for one case only: add `electric_resistivity_values` in `cases.json` (only `elec_resistivity` is patchable this way). For other props, edit the material file and re-run `--setup-only --rebuild`.

## Postprocessing pipeline

`Allrun_long` (called by each case's `Allrun`) runs the solver and pipes stdout through `scripts/parse_simulation_log.py`. After each write it calls `foamVTK.sh <time>` then immediately runs `scripts/analyze_meltpool_vtu.py --vtk-file <latest internal.vtu>`.

`analyze_meltpool_vtu.py` per-timestep outputs:
- `post-processing-data/vtu_meltpool_geometry.csv` — time-series of depth/width metrics
- `post-processing-data/vtu_sections/<N>_<time>.png` — 6-panel diagnostic plot

Metrics measured from `alpha.metal = 0.5` contour (keyhole boundary) and `T = Tsolidus` contour (melt-pool boundary) on the XY and YZ slices at the laser axis.

`COMPARE_DEPTH_FIELD` in the generated `analyze_meltpool_vtu.py` is patched to `keyholeDepth_um` when `depth_metric = "keyhole"`, else `meltPoolDepth_um`.

## Key invariants

- `setup_cases.py` uses regex-replace on known pattern anchors in OpenFOAM files. If you rename a comment header or parameter line in the template, the corresponding `replace_regex` call in `configure_case` will raise `RuntimeError: pattern not found`.
- `timeVsLaserPosition` and `timeVsLaserPower` are always fully regenerated (not patched) from `cases.json` values.
- Cases under `runs/` are generated artifacts. Only commit them if results need to be preserved; otherwise they can be regenerated with `--setup-only`.
