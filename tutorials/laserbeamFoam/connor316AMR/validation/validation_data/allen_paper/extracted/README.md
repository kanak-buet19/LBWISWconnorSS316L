# Allen et al. — Digitized Supplementary Data

**Source paper:**
> Troy R. Allen et al., "Energy-coupling mechanisms revealed through simultaneous keyhole depth and absorptance measurements during laser-metal processing," *Phys. Rev. Applied* (2020).

Data digitized from the supplementary PDF (Figs. S1–S13) using `batch_extract_pdf.py`.

---

## Cases

| File prefix       | Figure | Irradiance (MW/cm²) | Power (approx.) | Regime         |
|-------------------|--------|----------------------|-----------------|----------------|
| S1_0p23MWcm2      | S1     | 0.23                 | ~103 W          | Conduction     |
| S2_0p35MWcm2      | S2     | 0.35                 | ~156 W          | Conduction     |
| S3_0p40MWcm2      | S3     | 0.40                 | ~178 W          | Conduction     |
| S4_0p44MWcm2      | S4     | 0.44                 | ~196 W          | Transition     |
| S5_0p45MWcm2      | S5     | 0.45                 | ~200 W          | Transition     |
| S6_0p46MWcm2      | S6     | 0.46                 | ~205 W          | Transition     |
| S7_0p47MWcm2      | S7     | 0.47                 | ~209 W          | Transition     |
| S8_0p49MWcm2      | S8     | 0.49                 | ~218 W          | Keyhole onset  |
| S9_0p52MWcm2      | S9     | 0.52                 | ~231 W          | Keyhole        |
| S10_0p58MWcm2     | S10    | 0.58                 | ~258 W          | Keyhole        |
| S11_0p70MWcm2     | S11    | 0.70                 | ~311 W          | Keyhole        |
| S12_0p81MWcm2     | S12    | 0.81                 | ~360 W          | Keyhole        |
| S13_0p92MWcm2     | S13    | 0.92                 | ~410 W          | Strong keyhole |



---

Beam diameter: 238 µm (1/e², Gaussian)
Pulse duration: 10 ms (stationary)

## Axes (identical across all figures — stated explicitly in the paper)

| Quantity     | Range       | Notes                         |
|--------------|-------------|-------------------------------|
| Time         | 0 – 10 ms   | Uniform 500-bin grid (0.02 ms resolution) |
| Depth        | 0 – 800 µm  | Left axis, black curve        |
| Absorptance  | 0.3 – 0.9   | Right axis, red curve         |

---

## Folder Structure

```
extracted/
├── README.md
├── csv/
│   ├── S{N}_*_scatter.csv   ← all individual depth scatter points
│   └── S{N}_*_stats.csv     ← per-bin depth statistics + absorptance
└── plots/
    └── S{N}_*_verify.png    ← 2×3 verification figure
```

---

## Extraction Methodology

### Step 1 — Plot ROI from PDF vector graphics

The script uses PyMuPDF (`fitz`) to read the PDF as vector data rather than a rasterized image. MATLAB exports figures with the plot interior as a **white-filled rectangle** in the PDF vector stream. This rectangle is extracted directly from the drawing objects (`page.get_drawings()`), giving pixel-exact plot boundaries with no manual cropping.

- Filter: `fill == (1.0, 1.0, 1.0)`, width > 20% page width, height > 10% page height
- The plot interior rectangle was identical across all 13 figures (315.5 × 180.1 pts)
- Pages are rendered at 200 DPI → pixel coordinates computed as `point × (200/72)`

A border margin of **8 px** is trimmed from each side of the cropped plot interior before any pixel analysis, to prevent the frame lines from bleeding into the data region.

### Step 2 — Irradiance label extraction

The irradiance value (e.g. `0.46 MW/cm²`) for each plot is extracted directly from the PDF text layer using regex (`\d+\.\d+\s*MW`), so filenames are labelled automatically.

### Step 3 — Absorptance curve (red)

The red absorptance curve is extracted by HSV colour thresholding:

- Hue: 0 ± 15 (handles MATLAB red, which wraps near hue=0/180)
- Saturation ≥ 150, Value ≥ 80
- Morphological opening (3×3 ellipse kernel) removes isolated noise pixels
- Detected pixels are converted to data coordinates using the fixed axis calibration
- **One median value per x-bin** (500 bins over 0–10 ms) → single time series

The absorptance is a continuous line in the original, so one-median-per-bin is appropriate and accurate.

### Step 4 — Depth curve (black, scatter-aware)

The depth data is extracted using a dark-pixel threshold (`all(R, G, B) < 80`). Unlike the absorptance, the depth is plotted as **scatter dots** in the higher-irradiance cases (S5–S13), where multiple dots appear at different depths at the same time. Standard median binning would discard all but the middle dot.

Instead, within each x-bin the detected dark pixels are clustered vertically:
- Pixels within each x-bin are sorted by y-coordinate
- A **gap > 10 px** in y separates two distinct dots
- Each cluster must contain **≥ 8 pixels** to be counted (filters out tick-mark artifacts, which are 1–3 pixels)
- The median y of each cluster → one depth data point per dot

This recovers the full vertical scatter distribution at each time step.

### Step 5 — Artifact filtering

Two classes of artifacts were identified and removed:

| Artifact | Cause | Fix |
|----------|-------|-----|
| Deep points at t < 0.09 ms | Left frame border bleeding into first x-bins | 8 px border trim |
| Tick-mark bleed (isolated pixels) | Axis tick lines at frame edges | Min cluster size ≥ 8 px |

Points at t = 0.33–0.49 ms with depth 200–300 µm in S11–S13 were inspected and **retained** — they form monotonically ascending sequences consistent with rapid keyhole onset at high irradiance.

---

## Output File Formats

### `*_scatter.csv`
One row per detected depth point. Multiple rows can exist at the same `time_ms` when scatter dots stack at different depths.

| Column    | Unit | Description                    |
|-----------|------|--------------------------------|
| time_ms   | ms   | Time                           |
| depth_um  | µm   | Keyhole depth (individual dot) |

### `*_stats.csv`
One row per time bin (500 rows, uniform 0.02 ms grid). Depth statistics computed across all scatter points in each bin. Absorptance is the per-bin median.

| Column        | Unit | Description                          |
|---------------|------|--------------------------------------|
| time_ms       | ms   | Bin centre time                      |
| depth_mean    | µm   | Mean depth across scatter dots       |
| depth_std     | µm   | Std dev of scatter (spatial spread)  |
| depth_min     | µm   | Shallowest dot in bin                |
| depth_max     | µm   | Deepest dot in bin                   |
| depth_median  | µm   | Median depth                         |
| absorptance   | —    | Median absorptance (0.3–0.9)         |

NaN = no data detected in that time bin.

---

## Uncertainty Interpretation

**Depth std (`depth_std`)** represents the **physical scatter** of the X-ray depth measurement — at any given time, the keyhole depth is not a single value but a distribution. This is the correct uncertainty to use for CFD validation (compare simulated depth against `depth_mean ± depth_std`).

**Absorptance rolling std** (shown in verify plots, not saved to CSV) is **temporal noise** over a 1 ms window — useful for qualitative assessment only, since the raw absorptance already contains the full time-resolved signal.

---

## Verification Plots (`plots/*.png`)

Each figure has a 2-row × 3-column verification layout:

```
Row 0 (Depth):
  [overlay on image] | [all scatter points] | [mean ± 1σ band]

Row 1 (Absorptance):
  [overlay on image] | [raw extracted]      | [rolling mean ± 1σ]
```

Check these plots to confirm extraction quality before using the CSV data.
