# AISI 316L Additive Manufacturing — Simulation & Reproduction Reference

Extracted from: Abdali et al., *"Predictive tools for the cooling rate-dependent microstructure evolution of AISI 316L stainless steel in additive manufacturing,"* **Journal of Materials Research and Technology** 29 (2024) 5530–5538. DOI: 10.1016/j.jmrt.2024.03.008

This sheet collects the material data, laser/process parameters, modeling setup, and experimental results needed to reproduce the thermal–fluid simulations and validate against the paper's findings.

---

## 1. Material

| Property | Value |
|---|---|
| Alloy | AISI 316L stainless steel (per ASTM A276) |
| Form | Hot-rolled plate |
| Plate thickness | 12 mm |
| Solidus temperature (used for fusion boundary tracing) | 1650 K |
| Cooling-rate measurement window | between 1700 K and 1650 K (during solidification) |
| Thermophysical properties source | Wei et al., *Prog Mater Sci* 116 (2021) 100703 — Ref. [43] |

### Chemical composition (wt. %)

| C | Si | Mn | Cr | Ni | Mo | Cu | Ti | V | P | S | Fe |
|---|---|---|---|---|---|---|---|---|---|---|---|
| 0.019 | 0.514 | 1.32 | 16.51 | 10.46 | 2.14 | 0.401 | 0.019 | 0.089 | 0.025 | 0.008 | 68.50 |

### Derived equivalents (for constitution-diagram placement)

| Quantity | Value |
|---|---|
| Cr_eq (Hammar-Svensson) | 20.3 |
| Ni_eq (Hammar-Svensson) | 11.7 |
| Cr_eq / Ni_eq | 1.73 |

> Note: the paper plots the steel on multiple diagrams (WRC-1992, Schaeffler, Hammar-Svensson), each using its own equivalent formulas. The Cr_eq/Ni_eq = 1.73 value is the Hammar-Svensson basis used in Fig. 7.

---

## 2. Laser & Process Parameters

### DED (powderless directed energy deposition) — IRB 4600 machine

| Parameter | Value |
|---|---|
| Laser type | Fiber laser |
| Spot size (diameter) | 500 μm |
| Laser power | 1000 W |
| Scan speed | 11 mm/s |
| Melting regime | Conduction mode |

### SLM / LPBF — Mlab cusing R system

| Parameter | Value |
|---|---|
| Laser type | Fiber laser |
| Spot size (diameter) | 40 μm |
| Laser power | 95 W (constant) |
| Scan speeds (all runs) | 200, 300, 500, 600, 700 mm/s |
| Melting regime | Keyhole at low speeds → shallow conduction at highest speed |

**Fig. 4 simulated SLM condition specifically:** scan speed **200 mm/s**, power 95 W, 40 μm spot, keyhole mode. (This is the slowest run; the optical micrographs in Fig. 3 cover the 300–700 mm/s runs.)

---

## 3. Mechanistic Model Setup

General: 3D temperature field + fluid flow; Newtonian, incompressible fluid; laminar flow; solved by Finite Difference Method (FDM). Fusion boundary traced by the 1650 K solidus isotherm. Free surface captured with the Volume of Fluid (VOF) method. Initial temperature = ambient in both cases.

| Setting | DED model | SLM model |
|---|---|---|
| Domain size | 14 × 3 × 2 mm | 1.2 × 0.3 × 0.3 mm |
| Mesh cell size (rectangular cuboid) | 50 μm | 5 μm |
| Beam turn-off time | 0.76 ms | 3.0 ms |
| Candidate cross-section time | 0.47 ms | 2.6 ms |
| Recoil pressure | Neglected (T_max < evaporation T) | Included (keyhole) |

> Discrepancy flag: the Fig. 2a / Fig. 4a image labels read "0.76 s / 0.47 s" and "3.0 ms / 2.6 ms" respectively, while the figure captions state milliseconds in both cases. The caption values (ms) are used above; worth verifying against source data if timing matters to your run.

### Governing equations

| Eq. | Description | Form |
|---|---|---|
| (1) | Mass conservation | ∇·**v** = 0 |
| (2) | Momentum | ∂**v**/∂t + (**v**·∇)**v** = −(1/ρ)∇p + (μ/ρ)∇²**v** − K**v** + f |
| (3) | Energy | ∂h/∂t + **v**·∇h = (1/ρ)∇·(κ∇T) |
| (4) | VOF (free surface) | ∂F/∂t + (**v**·∇)F = 0 |
| (5) | Gaussian heat source | q = (2ηP / πr²)·exp(−2(x²+y²)/r²) |
| (6) | Buoyancy force | F_b = −ρgβ(T − T_l) |
| (7) | Surface tension (Marangoni) | γ = γ₀ + (dγ/dT)(T − T₀) |
| (8) | Recoil pressure | P_r = A·B₀·T_s^(1/2)·exp(−U/T_s) |
| (9) | Surface heat balance (BC) | κ ∂T/∂n = Q − Q_rad − Q_conv − Q_evap |
| (10) | Radiation loss | Q_rad = σε(T⁴ − T₀⁴) |
| (11) | Convection loss | Q_conv = h_c(T − T₀) |
| (12) | Evaporation loss | Q_evap = ρ·V_evap·T |

### Recoil-pressure constants (Eq. 8)

| Constant | Value |
|---|---|
| A (atmospheric pressure) | 0.55 |
| B₀ | 1.78 × 10¹⁰ |

### Symbol reference

P = laser power · η = energy absorption coefficient · r = laser beam radius · **v** = velocity vector · ρ = density · p = hydrodynamic pressure · μ = kinematic viscosity · K = drag coefficient · f = body-force acceleration · h = enthalpy · κ = thermal conductivity · F = volume fraction of fluid · T = temperature · T_s, T_l = solidus/liquidus temps · T₀ = ambient temp · β = thermal expansion coefficient · γ₀ = surface tension of pure metal · U = internal energy · **n** = free-surface normal · σ = Stefan-Boltzmann coefficient · ε = radiation emissivity · h_c = heat transfer coefficient · V_evap = evaporation recession speed

> Note: η (absorption) and the assumption that keyhole multiple reflections are "compensated by increasing nominal power" are stated qualitatively — exact numeric values for η, β, γ₀, dγ/dT, K, h_c, ε, etc. are not tabulated in the paper and must be taken from the 316L property set in Ref. [43].

---

## 4. Experimental & Computed Results (validation targets)

### 4.1 Cooling rates vs. distance from melt pool bottom (Fig. 5)

**DED weld:**

| Distance from melt pool bottom (μm) | Cooling rate (K/s) |
|---|---|
| ~50 | 2650 |
| ~150 | 2250 |
| ~250 | 2210 |
| ~350 | 1360 |

- Range: **1360–2650 K/s** · Average: **2118 K/s**

**SLM weld:**

| Distance from melt pool bottom (μm) | Cooling rate (K/s) |
|---|---|
| ~20 | 972,400 |
| ~80 | 384,300 |
| ~140 | 319,500 |
| ~200 | 286,400 |

- Average: **490,650 K/s**

Key relationship: cooling rate decreases with distance from the melt pool bottom in both cases; SLM is at least **two orders of magnitude** higher than DED at a given distance.

### 4.2 DED residual ferrite (Fig. 1, transverse section)

| Region (toward center) | Residual ferrite δ (vol %) |
|---|---|
| Near melt pool bottom | 17 |
| Middle | 13 |
| Center | 12 |

- Solidification mode: **FA** (primary ferrite) throughout; vermicular/lacy ferrite near fusion boundary, shifting to spherical toward center.

### 4.3 SLM A-mode area fraction (Fig. 3, ∑A = dark/A-mode area %)

| Scan speed V (mm/s) | ∑A (%) |
|---|---|
| 300 | 51 |
| 500 | 52 |
| 600 | 65 |
| 700 | 54 |

- Total residual ferrite: **< 0.5 vol %** across all SLM conditions.
- Dual solidification mode: **A** (austenite) at fusion boundary → **F** (ferrite) in interior.
- ∑A trends up with scan speed, except at the highest speed where the mode shifts to very shallow conduction.

### 4.4 Predicted ferrite (constitution-tool targets)

| Tool | Prediction for this steel |
|---|---|
| WRC-1992 | FA mode, ~6–8 vol % residual ferrite (underestimates DED measurement) |
| Hammar-Svensson (Fig. 7) | FA for DED cooling rate; modification signaled at SLM rate |
| Schaeffler @ 1.9×10³ K/s (Fig. 8b) | Austenite + <10 vol % ferrite (DED-comparable) |
| Schaeffler @ 2×10⁵ K/s (Fig. 9) | Type A′ → fully austenitic (SLM-comparable) |

---

## 5. ORFN Neural-Network Model (ferrite number prediction)

| Aspect | Detail |
|---|---|
| Architecture | 14 input nodes (13 elements + log of cooling rate), 6 hidden nodes, 1 output node (ferrite number) |
| Parameters source | Vitek et al., *Weld J* 82 (2003) — Refs. [32, 33]; normalization, weights, de-normalization tabulated there |
| Implementation | Excel spreadsheet, per formulation in Vitek et al. Ref. [20] |
| Key prediction | Ferrite content rises with cooling rate up to ~30,000 K/s, then declines beyond it |
| DED window result | Over 1360–2650 K/s, ferrite predicted to change by only a few percent (matches observed DED trend qualitatively) |

---

## 6. Reproduction Checklist / Gaps

What the paper gives you directly: full composition, both machines' power/speed/spot, domain & mesh sizes, beam timing, governing equations, recoil constants, and validation data (cooling rates, ferrite fractions, ∑A).

What you must source elsewhere before running:
- Numeric 316L thermophysical properties (κ, ρ, c_p, μ, β, γ₀, dγ/dT, latent heat, liquidus T) → take from Ref. [43] (Wei et al., 2021).
- Absorption coefficient η and the "nominal power increase" factor used to emulate keyhole multiple reflections (qualitative only in paper).
- Boundary-condition constants: h_c, ε, V_evap model.
- ORFN weights/normalization → Refs. [32, 33].

Recommended validation order: (1) match fusion boundary to the 1650 K solidus contour (Figs. 2c, 4b); (2) match the cooling-rate-vs-distance curves in Fig. 5; (3) compare predicted solidification mode/ferrite against Sections 4.2–4.3.