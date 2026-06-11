# Pulsed Laser Welding of AISI 316L SS Under Space Conditions
**Riffel et al. (2026) — *Materials & Design*, DOI: 10.1016/j.matdes.2026.116100**

---

## 1. Material

| Property | Value |
|---|---|
| Alloy | AISI 316L Stainless Steel |
| Plate thickness | 3.0 mm (0.118 in) |
| Plate dimensions | 101.6 mm × 101.6 mm (4 in × 4 in) |
| Surface prep | Sanded with 120-µm grit (to increase roughness and laser absorptivity) |

### Chemical Composition (wt%)
| C | Co | Cr | Cu | Mn | Mo | N | Nb | Ni | P | S | Si | Ti | Fe |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| 0.01 | 0.37 | 16.82 | 0.42 | 1.26 | 2.05 | 0.05 | 0.07 | 10.06 | 0.03 | 0.002 | 0.34 | 0.07 | Bal |

### Key Thermophysical Constants Used in Modeling
- Density (ρ): **7966 kg/m³** (constant)
- Creq/Nieq ratio (WRC-1992): **18.9 / 11.8 = 1.59** → falls in austenitic massive transformation region
- Kurz-Fisher constants for AISI 316L: **A = 80, n = 0.33**
- Marangoni number (estimated): **~2.6 × 10⁴** (surface tension dominated flow)
- Rayleigh number: **Ra = 53.9** (Earth), **Ra = 1.6** (low gravity)

---

## 2. Laser System & Process Parameters

| Parameter | Value |
|---|---|
| Laser type | IPG Photonics YLR-150/1500-QCW-MM-AC (pulsed fiber laser) |
| Wavelength | 1070 nm |
| Max peak power | 1.5 kW |
| **Pulse power used** | **1050 W** |
| **Pulse duration** | **7.1 ms** |
| **Pulse frequency** | **10 Hz** |
| **Travel speed** | **4 mm/s** |
| Spot size | 191 µm (measured via PRIMES FocusMonitor FM+) |
| Weld type | Bead-on-plate |
| Mode | Keyhole (finger-like penetration + crown-shaped surface) |

> **Note for simulation:** The highest pulse-power setting was selected as it produces the largest molten zones, maximizing detectable gravity-dependent differences.

---

## 3. Experimental Setup & Environments

### Vacuum Chamber
- Volume: ~0.6 m³
- Vacuum level: **1.33 × 10⁻³ Pa (1 × 10⁻⁶ Torr)**
- Pumps: 2× Varian TV301 turbomolecular + 2× Edwards XDS10 dry scroll roughing pumps
- Motion: axial system moves the part; laser is stationary
- Ambient temperature: ~22 °C (ground lab); ~15 °C (aircraft cabin — may affect penetration)

### Gravity Conditions (Parabolic Flight)
| Condition | Gravity Level | Avg. Acceleration |
|---|---|---|
| On-ground (baseline) | 1 g | ~9.81 m/s² |
| Lunar gravity | ~0.15 g | 0.150 ± 0.003 g (1.47 m/s²) |
| Low gravity | ~0.03 g | 0.033 ± 0.002 g (0.29 m/s²) |

- Gravity recorded via **Summit Instruments 35203A accelerometer** at 25 Hz, mounted on chamber exterior
- No convective heat transfer applied in model (full vacuum — radiation only at boundaries)

### Weld Sequence & Lens Buildup Issue
- Total campaign: 69 welds across 2 days (316L SS, Al 2219-T87, Ti6Al4V)
- Welds analyzed here: **Weld 12** (Earth vs. Lunar) and **Weld 30** (Earth vs. Low-G)
- By Weld 30, metallic vapor had progressively coated the cover glass → reduced effective power → narrower beads
- Ground replicas were made at OSU under identical vacuum/parameters to match optical buildup state

---

## 4. Characterization Methods

| Method | Equipment | Purpose |
|---|---|---|
| Computed Tomography (CT) | NSI microfocus CT, X-RayWorX [P19-702] source, Perkin Elmer [XRD1620/1621] detector | 3D porosity characterization |
| CT software | VGStudio (acquisition), 3D Slicer (analysis) | Pore count, volume, aspect ratio |
| CT voxel size | **3.77 × 10⁻⁵ mm³** | Note: sub-micron pores may be missed |
| Segmentation | Otsu's + Huang's thresholding | Validated against cross-sectional porosity area |
| Optical Microscopy | Leco Olympus GX50, Keyence VHX-S750E | Weld macro/microstructure |
| SEM | Thermo Scientific Quattro | Microstructure, EDS |
| EBSD | EDAX Velocity camera (~4,500 pts/s) | Grain size, aspect ratio, texture |
| EBSD settings | 70° tilt, 20 mm WD, 20 kV, 5 µm spot | — |
| EBSD post-processing | MTEX 5.2 | IPF maps, grain size, aspect ratio per ASTM E2627 |
| Chemical analysis | OES (Optical Emission Spectroscopy) | Composition verification |

---

## 5. Thermal / FEM Model

### Formulation
- Software: **Simufact Welding 2023.b**
- 3D transient heat conduction (Fourier's law + 1st law of thermodynamics)
- **No fluid dynamics / buoyancy modeled** (conduction-only)
- Gravity does not appear explicitly in governing equations

### Governing Equation

$$\frac{\partial}{\partial x}\!\left[k(T)\frac{\partial T}{\partial x}\right] + \frac{\partial}{\partial y}\!\left[k(T)\frac{\partial T}{\partial y}\right] + \frac{\partial}{\partial z}\!\left[k(T)\frac{\partial T}{\partial z}\right] + \dot{Q}(T,t) = \rho(T)\,C_p(T)\,\frac{\partial T}{\partial t}$$

### Heat Source Model
Hybrid: **Conical volumetric** (finger penetration) + **Gaussian surface** (crown profile)

Volumetric heat source:

$$Q(x,y,z) = Q_0 \exp\!\left[-\frac{x^2 + y^2}{r_0(z)^2}\right]$$

$$r_0(z) = r_u + \frac{r_i - r_u}{z_i - z_u}\,(z - z_u)$$

Radiation boundary (Stefan-Boltzmann):

$$\frac{\dot{Q}_r}{A} = -\varepsilon\,\sigma\,(T^4 - T_\infty^4)$$

### Boundary Conditions
- **No convective coefficient** (vacuum conditions)
- Radiation at all external surfaces
- Reference temperature T∞ = **22 °C** (ground and flight)
- Material properties k(T) and Cp(T): temperature-dependent (see paper Supplementary Fig. S2)
- Density: constant at 7966 kg/m³

### Model Calibration
- Heat source parameters tuned by matching simulated vs. experimental weld cross-sections (molten zone shape/size)
- Validated in Supplementary Figs. S4 and S5

---

## 6. Key Quantitative Results

### 6.1 Weld Bead Geometry

| Weld | Gravity | Avg. Width |
|---|---|---|
| Weld 12 | Ground | 0.91 ± 0.05 mm |
| Weld 12 | Lunar (0.15g) | 0.95 ± 0.06 mm |
| Weld 30 | Ground | 0.73 ± 0.08 mm |
| Weld 30 | Low-G (0.03g) | 0.72 ± 0.07 mm |

- All pairs **statistically equivalent** (t-test, 95% CI)
- Fusion zone area: **0.21 mm²** (Weld 12 both conditions), **0.18 mm²** (Weld 30 both conditions)

### 6.2 Porosity (CT Data)

| Weld | Gravity | # Pores | Total Vol (mm³) | Avg. Vol (mm³) | X/Z ar | Y/Z ar |
|---|---|---|---|---|---|---|
| Weld 12 | Ground | 32 | 0.16 | 5.1×10⁻³ ± 0.002 | 0.86 ± 0.09 | 0.84 ± 0.10 |
| Weld 12 | Lunar | 32 | 0.10 | 3.0×10⁻³ ± 0.002 | 0.94 ± 0.06 | 0.85 ± 0.07 |
| Weld 30 | Ground | 34 | 0.11 | 3.2×10⁻³ ± 0.001 | 0.86 ± 0.14 | 0.95 ± 0.13 |
| Weld 30 | Low-G | 13 | 0.03 | 2.6×10⁻³ ± 0.001 | 0.96 ± 0.10 | 0.76 ± 0.09 |

- Porosity reduction: **38.3%** (Earth → Lunar), **69.1%** (Earth → Low-G)
- Reduced gravity → fewer, smaller, more **spherical** pores (X/Z ar closer to 1.0)
- Ground welds → more **elongated** pores in Z-direction (gravity vector direction)

### 6.3 Solidification & Cooling Rates (FEM + Kurz-Fisher)

| Location | FEM dT/dt | Kurz-Fisher dT/dt | FEM R | KF R |
|---|---|---|---|---|
| Weld center (max) | 2.2 × 10⁵ K/s | 5.6 × 10⁵ K/s | 127.5 mm/s | 199.1 mm/s |
| Near fusion boundary (min) | 3.6 × 10⁵ K/s | 1.4 × 10⁵ K/s | 77.3 mm/s | 50.6 mm/s |

- Thermal gradient (G): **2807.3 K/mm** (FEM, within molten zone)
- Solidification complete **before** next pulse (10 Hz → 100 ms between pulses; pulse = 7.1 ms)

### 6.4 Grain Size & Aspect Ratio (EBSD)

| Weld | Gravity | Grain Size (µm) | Aspect Ratio |
|---|---|---|---|
| Weld 12 | Ground | 12.3 ± 0.8 | 4.0 ± 0.3 |
| Weld 12 | Lunar | 12.0 ± 0.7 | 3.8 ± 0.2 |
| Weld 30 | Ground | 10.9 ± 0.8 | 2.7 ± 0.1 |
| Weld 30 | Low-G | 9.8 ± 0.5 | 2.7 ± 0.1 |

- All conditions: **statistically equivalent** (overlapping confidence intervals)
- Columnar grain growth throughout fusion zone in all conditions

### 6.5 Microstructure / Solidification Modes
- **Periphery of fusion zone:** Planar → Cellular dendritic → Fully austenitic (A)
- **Weld center:** Ferrite–Massive Austenite (F/MA) transformation
- Same solidification sequence observed across all gravity levels
- Liquated grain boundaries present in all conditions (not gravity-specific)

---

## 7. Dimensionless Numbers

| Number | Formula | Earth Value | Low-G Value | Significance |
|---|---|---|---|---|
| Marangoni (Ma) | $\displaystyle\frac{\|d\gamma/dT\|\,\Delta T\,w}{\mu\,\alpha}$ | ~2.6 × 10⁴ | ~2.6 × 10⁴ | Surface tension dominates over buoyancy |
| Rayleigh (Ra) | $\displaystyle\frac{g\,\beta\,\Delta T\,L^3}{\nu\,\alpha}$ | 53.9 | 1.6 | Buoyancy convection much weaker in low-G |

> The Marangoni effect is gravity-independent; buoyancy (Ra) scales directly with g. This explains why microstructure is similar across gravities while porosity morphology differs slightly.

### Marangoni Number — Full Calculation (from Supplementary Material)

**Sulfur content:** S = 0.002 wt% = 20 ppm → drives **outward** Marangoni flow (hot center → cold periphery)

| Property | Symbol | Value | Source |
|---|---|---|---|
| Surface tension gradient | dγ/dT | −2.4 × 10⁻⁴ N/(m·K) | Literature (316L in vacuum at this S content) |
| Molten pool ΔT | ΔT | ~1000 K (center ~2400 K, periphery ~1400 K) | FEM simulation |
| Molten pool width (proxy) | w | 9.5 × 10⁻⁴ m (0.95 mm, max bead width — Weld 12 lunar) | Experiment |
| Dynamic viscosity at 2400 °C | μ | 1.98 × 10⁻³ Pa·s (1.98 cP) | Literature |
| Thermal diffusivity at 2400 °C | α | 4.4 × 10⁻⁶ m²/s | Literature |

$$\text{Ma} = \frac{|d\gamma/dT|\,\Delta T\,w}{\mu\,\alpha} = \frac{(2.4\times10^{-4})\times 1000\times(9.5\times10^{-4})}{(1.98\times10^{-3})\times(4.4\times10^{-6})} = 2.6\times10^{4}$$

---

## 8. Hypothesized Mechanisms

### Why pores are more spherical in low gravity
- On Earth: pores formed at keyhole collapse can rise with solidification front due to buoyancy → elongated shape in Z
- In low-G: buoyancy absent → pores stay in place → spherical morphology
- Alternative: gravity affects **keyhole collapse dynamics** (inward collapse + Marangoni) → more vapor entrapment on Earth

### Why microstructure is unaffected by gravity
- Pulsed LBW solidification is complete between pulses (~92.9 ms cooling window at 10 Hz)
- No sustained molten pool → insufficient time for buoyancy-driven convection to alter grain growth
- Rapid solidification (up to 199.1 mm/s) dominates over gravity-driven effects
- Contrasts with CW/steady-state processes where buoyancy has time to act

---

## 9. Unexplained Observations / Open Questions

| Observation | Status |
|---|---|
| Weld 30 has smaller fusion zone area (0.18 vs. 0.21 mm²) and deeper penetration than Weld 12 | Cause not clearly understood; possibly lens focal shift from vapor deposition |
| Slight texture difference (ground: max 2.6, low-G: max 2.1 in ⟨111⟩) | Attributed to small residual buoyancy in ground welds; needs CFD confirmation |
| Low-G Weld 30 has fewer pores (13 vs. 34) but also lower heat input | Cannot fully decouple gravity effect from lens degradation effect |
| Transverse (Y/Z) pore aspect ratio shows no clear trend with gravity | Mechanism unclear; overlapping confidence intervals |

---

## 10. Experimental Confounders (Important for Simulation Validation)

- **Lens buildup:** metallic vapor progressively coated cover glass → reduced effective power over weld sequence
- **Aircraft vibrations:** near parabola peak may influence molten pool flow and pore dynamics
- **Chamber relaxation:** aluminum vacuum chamber may deflect slightly under reduced gravity → focal shift
- **Cabin temperature:** ~15 °C during flight vs. 22 °C in lab → may affect penetration depth
- **Rolling direction:** ground plates perpendicular, flight plates parallel to rolling direction (EBSD confirmed negligible effect)
- **Statistical limitation:** only 1 weld per gravity level per condition analyzed; high experimental cost limits repetitions

---

## 11. Future Work

- [ ] CFD modeling of molten pool under different gravity levels (fluid dynamics + buoyancy + Marangoni)
- [ ] Additional parabolic flight experiments with **CW laser** for comparison with pulsed
- [ ] High-speed imaging to capture transient keyhole and molten pool behavior (for model validation)
- [ ] Higher-resolution sub-micron CT imaging (current voxel size misses sub-micron pores)
- [ ] Studies at **cryogenic temperatures**
- [ ] Experiments with **aluminum alloys** and **titanium alloys** (gravity effects may differ due to lower viscosity)
- [ ] More repetitions to improve statistical confidence

---

## 12. Supplementary Material Details

### FEM Model — Geometry & Thermal Properties (Figs. S2 & S3)
- Temperature-dependent **k(T)** and **Cp(T)** for 316L used as inputs (tabulated in Fig. S2 of paper)
- Heat source geometry parameters (conical radii rᵢ, rᵤ, and conic depth zᵢ−zᵤ) defined in Fig. S3
- Hybrid conical-Gaussian heat source calibrated to reproduce the crown-surface + finger-penetration morphology

### FEM Validation Method (Figs. S4 & S5)
- **Cross-section (Fig. S4):** Simulated molten zone boundary compared to experimental transverse macrograph; heat source tuned until shapes match
- **Longitudinal (Fig. S5):** Isometric view of full pulsed bead including check spot weld prior to bead start; transverse cross-section identifies center and fusion boundary locations; longitudinal cross-section used to extract G, R, and dT/dt along bead

### Base Material Rolling Direction Check (Fig. S6)
- Ground sample: rolled **perpendicular** to weld direction; Flight sample: rolled **parallel**
- EBSD confirmed both base materials have similar equiaxed grain sizes and distribution patterns
- Texture: ground max intensity **2.2**, flight max intensity **2.0** — both weak, similar cluster positions
- Ground EBSD data rotated **90° around TD** to align reference frames for comparison
- **Conclusion:** different rolling directions have negligible effect on weld microstructure

### EDS Maps of F/MA Region (Fig. S7)
- White dots in SEM (weld center) are **enriched in C and Ni**, **depleted in Cr and Fe**
- Identified as **primary dendrite arm spaces** (interdendritic regions), not second-phase particles

### Solidification Mode Diagram (Fig. S8)
- Creq/Nieq = **1.59** with the solidification rate places this work in the **austenitic massive transformation** region
- Confirms the F/MA solidification mode observed experimentally

### Spatter Trajectory Validation (Fig. S1 & Video S1)
- Camera outside vacuum chamber tracks spatter as qualitative gravity confirmation
- **Earth (1g):** spatter follows **parabolic arcs** — gravity clearly evident
- **Low-G (0.03g):** spatter travels in **straight lines** — no gravitational deflection
- Only particles visible ≥3 frames were tracked; single-camera setup limits quantitative analysis
- Future: add perpendicular camera for 3D velocity/position measurement

### Solidification Completeness Between Pulses (Video S2)
- Thermal simulation confirms solidification is **complete before the next pulse fires**
- 10 Hz → 100 ms inter-pulse interval; 7.1 ms pulse → **~92.9 ms available for full solidification**
- This is the fundamental reason gravity cannot influence grain structure in pulsed LBW

---

## 13. Quick Reference for Simulation Setup

```
Material:         AISI 316L SS
                  ρ = 7966 kg/m³ (constant)
                  k(T), Cp(T) = temperature-dependent (see Fig. S2)
                  S content = 0.002 wt% → outward Marangoni flow

Laser:            1070 nm, 1050 W peak, 7.1 ms pulse, 10 Hz, 4 mm/s travel, 191 µm spot
                  Inter-pulse interval = 100 ms → solidification complete between pulses

Heat source:      Conical-Gaussian hybrid (see Fig. S3 for geometric parameters)
                  Target: crown-shaped surface + finger-like penetration
                  Calibrate to: 0.21 mm² fusion zone (Weld 12) or 0.18 mm² (Weld 30)

Boundary:         Radiation only (Stefan-Boltzmann), T∞ = 22°C
                  No convective coefficient — full vacuum
                  ε (emissivity) = from 316L literature

Vacuum:           1.33×10⁻³ Pa (1×10⁻⁶ Torr)

Gravity cases:    1g = 9.81 m/s²  |  0.15g = 1.47 m/s²  |  0.033g = 0.29 m/s²

Fluid props       μ = 1.98×10⁻³ Pa·s at 2400°C
(for CFD):        α = 4.4×10⁻⁶ m²/s at 2400°C
                  dγ/dT = −2.4×10⁻⁴ N/(m·K)
                  Ma ≈ 2.6×10⁴  |  Ra(1g) = 53.9  |  Ra(0.03g) = 1.6

Target outputs:   Molten zone shape, G, R, dT/dt, pore morphology
Validation:       G = 2807.3 K/mm; R = 77.3–127.5 mm/s (FEM); R = 50.6–199.1 mm/s (KF)
                  Spatter: parabolic at 1g, straight lines at 0.03g
```