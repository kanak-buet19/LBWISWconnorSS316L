# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## What This Is

NIST AMB 2022 experimental dataset: simultaneous laser absorptance + high-speed X-ray imaging during laser melting of Ti-6Al-4V. DOI: https://doi.org/10.18434/mds2-2525

Two experimental configurations:
- **Spot (stationary laser):** 2.0 ms pulse, laser held fixed
- **Scan (moving laser):** 700 mm/s scan speed, "sky-write" strategy

## Directory Layout

```
beam_profile/     beam profile CSV + visualization notebook
stationary_laser/ spot experiment (2 ms pulse): absorption CSV, X-ray zips, AVI
moving_laser/     scan experiment (700 mm/s): absorption CSV, X-ray zips, AVI
notebooks/        data synchronization notebook
docs/             README, sample diagram PDF
```

## Visualization Scripts

All scripts run from repo root. Require: `numpy`, `pandas`, `matplotlib`, `Pillow`.

```bash
# Absorption timeseries (both experiments, 3-panel figure)
python scripts/plot_absorption.py

# Beam profile heatmap + cross-sections + 1/e² circle
python scripts/plot_beam_profile.py --power 201.3

# Single X-ray frame synced to absorption timeseries
python scripts/plot_xray_sync.py --exp stationary --frame 90
python scripts/plot_xray_sync.py --exp moving --frame 90 --save   # save PNG

# Interactive frame browser (arrow keys or slider)
python scripts/plot_xray_sync.py --exp stationary --browse
```

`--zip` flag selects image set: `raw` (default), `captioned`, or `absorption`.

## Running the Notebooks

```bash
jupyter notebook
```

Both notebooks have hardcoded NIST server paths — update before running.

- [beam_profile/beam_profile_imaging.ipynb](beam_profile/beam_profile_imaging.ipynb) — visualize beam irradiance. Set `LaserPower` (W) and path to `beam_profile/beam_profile_5p5um_normalized.csv`.
- [notebooks/data_sync_example.ipynb](notebooks/data_sync_example.ipynb) — sync X-ray frames to absorption timeseries. Set paths to `{experiment}/absorption_data.csv` and unzipped image directories.

## Data Schema

### Absorption CSVs (`*_Calibrated Absorption Data.csv`)

9 columns (index starts at 1):

| Col | Name | Description |
|-----|------|-------------|
| 1 | index | Row number |
| 2 | `Time` | Seconds |
| 3 | `InputLaser` | Laser pulse profile (W). Ignore sharp artifact in first ~hundreds of ns. |
| 4 | `AbsoluteAbsorption` | Absorbed power (W) |
| 5 | `AbsAbsorptionUncertainty` | Expanded uncertainty in absorbed power (W) |
| 6 | `RelativeAbsorption` | Percent absorption (%) |
| 7 | `CameraTrigger` | Binary: 1 = camera triggered |
| 8 | `FrameTrigger` | Binary: 1 = frame exposure in progress (leading edge = start) |
| 9 | `FrameNumber` | Frame number; matches last 3 digits of image filenames |

A frame is captured only when both `CameraTrigger` AND `FrameTrigger` are 1.

### Image Files

Raw X-ray images are `.tif` files. Last 3 digits of filename = `FrameNumber` (zero-padded). Syncing frames to absorption: find the row index where `FrameNumber` first equals the target frame number — that gives the timestamp.

### Beam Profile CSV

60×60 matrix, pixel spacing = 5.5 μm. Values are normalized so sum×pixel_area = 1. Multiply by laser power (W) to get irradiance (W).

## Key Experimental Parameters

- **Laser wavelength:** 1070 nm (Yb-doped fiber)
- **Beam spot at sample surface:** 122.5 ± 3.0 μm (1/e²) — sample was 2.8 mm below beam waist
- **Beam waist (focus):** 49.5 ± 5 μm — beam profile CSV is measured here, not at sample
- **Angle of incidence:** 7° from surface normal
- **X-ray frame rate:** 50,000 fps, 2.5 μs exposure
- **Sample:** Ti-6Al-4V (NIST SRM 654b), ~300 μm thick, polished to specular finish
- **Atmosphere:** Argon (vacuum pump + backfill)
- **Facility:** APS 32-ID-B beamline, Argonne National Laboratory

## Frame Synchronization Pattern

```python
import pandas as pd
data = pd.read_csv('*_Calibrated Absorption Data.csv', index_col=0, float_precision="high")
frame_indices = {}
for n in range(1, int(data.FrameNumber.max()) + 1):
    row = data.loc[data['FrameNumber'] == n]
    if not row.empty:
        frame_indices[n] = row.index[0]
# frame_indices[N] gives the CSV row index at start of frame N
```

## Reference Papers

Core physics and methodology:
1. Simonds et al. (2021) *Applied Materials Today* 23, 101049 — causal relationship between melt pool geometry and absorption
2. Khairallah, Sun, Simonds (2021) *Additive Manufacturing Letters* 1, 100002 — periodic oscillations as precursor to pore-generating turbulence
3. Simonds et al. (2020) *Procedia CIRP* 94, 775–779 — simultaneous X-ray and absorption measurement method
