# Zhang et al. (2024) — 316L Bare-Plate Melt Pool: Full Reproduction Reference

**Source:** Zilong Zhang, Tianyu Zhang, Can Sun, Sivaji Karna, Lang Yuan,
*"Understanding Melt Pool Behavior of 316L Stainless Steel in Laser Powder Bed Fusion Additive Manufacturing,"*
**Micromachines** 15(2):170, 2024.
DOI: 10.3390/mi15020170 | Open Access (CC BY) | PMID: 38398900

---

## 1. Experimental Setup

| Item | Value |
|---|---|
| Machine | Aconity MIDI (Germany) |
| Laser type | Single-mode fiber laser |
| Max laser power | 1000 W |
| Laser spot radius (r_laser) | 50 μm (diameter = 100 μm, constant throughout) |
| Atmosphere | Argon |
| O₂ level | < 100 ppm |
| Substrate (bare plate) | 316L stainless steel |
| Powder (powder plate cases) | 316L, gas-atomized, 15–45 μm particle size, mean 30 μm (Carpenter Technology Corp., USA) |
| Powder layer thickness | 50 μm |
| Cross-section etchant | 75 vol% HCl + 25 vol% HNO₃ |
| Cross-section method | Wire-EDM, perpendicular to scan direction, at mid-track |
| Surface topography tool | Keyence VHX-5000 digital optical microscope |

---

## 2. Process Parameters — Single-Track Experiments

All 9 cases run on **both bare plate (BP) and powder plate (PP)**. Each line printed twice, length 20 mm.

| Case | Laser Power (W) | Laser Speed (m/s) | Linear Energy Density (J/m) | Series |
|---|---|---|---|---|
| N01 | 260 | 0.52 | 500 | Constant power |
| N02 | 260 | 0.87 | 300 | Constant power + constant LED |
| N03 | 260 | 1.30 | 200 | Constant power |
| N04 | 260 | 1.47 | 177 | Constant power + constant speed |
| N05 | 260 | 2.20 | 118 | Constant power |
| N06 | 440 | 1.47 | 300 | Constant speed + constant LED |
| N07 | 620 | 2.07 | 300 | Constant LED |
| N08 | 800 | 2.67 | 300 | Constant LED |
| N09 | 620 | 1.47 | 423 | Constant speed |

> Linear energy density = P_laser / V_laser

**Three experimental series:**
- Constant power (260 W): N01–N05 → varies speed and LED
- Constant speed (1.47 m/s): N04, N06, N09 → varies power and LED
- Constant LED (300 J/m): N02, N06, N07, N08 → varies power and speed simultaneously

---

## 3. Cube Sample Parameters (for surface roughness validation)

| Condition | Power (W) | Speed (m/s) |
|---|---|---|
| Low power/speed | 260 | 0.52 |
| High power/speed | 440 | 1.47 |

| Cube parameter | Value |
|---|---|
| Cube size | 10 × 10 × 10 mm |
| Hatch spacing | 100 μm |
| Layer thickness | 30 μm |
| Scan strategy | Back-and-forth with 90° rotation between consecutive layers |

---

## 4. Key Experimental Observations (Validation Targets)

### 4.1 Melting Mode vs. Process Conditions

| Condition | Mode | Criterion |
|---|---|---|
| W/D ratio < 1 | Keyhole | Depth > half-width |
| W/D ratio > 1 | Conduction | Depth < half-width |

**Constant power (260 W), bare plate — mode transitions:**

| Case | Speed (m/s) | LED (J/m) | Mode (BP) |
|---|---|---|---|
| N01 | 0.52 | 500 | Keyhole (W/D = 0.633) |
| N02 | 0.87 | 300 | Conduction |
| N03 | 1.30 | 200 | Conduction (swell-undercut starts) |
| N04 | 1.47 | 177 | Conduction + swell-undercut |
| N05 | 2.20 | 118 | Conduction + balling |

**Constant speed (1.47 m/s), bare plate — mode transitions:**

| Case | Power (W) | Mode (BP) |
|---|---|---|
| N04 | 260 | Conduction |
| N06 | 440 | Transition to keyhole (onset) |
| N09 | 620 | Keyhole |

**Constant LED (300 J/m), bare plate:**

| Case | Power (W) | Speed (m/s) | Mode (BP) |
|---|---|---|---|
| N02 | 260 | 0.87 | Conduction |
| N06 | 440 | 1.47 | Keyhole |
| N07 | 620 | 2.07 | Keyhole |
| N08 | 800 | 2.67 | Keyhole collapses (extreme recoil) |

### 4.2 Surface Morphology Regimes (Bare Plate)

| Morphology | Condition |
|---|---|
| Flat, continuous | 260 W, ≤ 0.87 m/s |
| Swell-undercut | 260 W, ≥ 1.30 m/s; all 1.47 m/s cases |
| Balling (discontinuous) | 260 W, ≥ 1.47 m/s |
| Keyhole void trail | 260 W, 0.52 m/s; high power cases |

> **Swell-undercut:** higher height at melt pool center, lower height at boundaries.
> Begins when recoil pressure expels liquid outward but pool refills only partially.

### 4.3 Bare Plate vs. Powder Plate Differences

| Aspect | Bare Plate | Powder Plate |
|---|---|---|
| Absorptivity | Lower (smooth surface) | Higher (multi-reflection in powder bed) |
| Melt pool height | Lower | Higher (more mass melted) |
| Instability | Lower | Higher (uneven mass addition) |
| Spatter | Less | More (especially high power/speed) |
| Balling height | Lower | Higher under same conditions |

---

## 5. Thermophysical Properties of 316L (Table 2 in paper)

These are the values used to validate the CFD model against bare-plate and powder-plate experiments.

### 5.1 Density

| Phase | Symbol | Value | Unit |
|---|---|---|---|
| Solid | ρ_s | 7950 | kg/m³ |
| Liquid | ρ_l | 8200 − 0.77·T | kg/m³ (T-dependent) |
| Gas | ρ_g | 1.22 | kg/m³ |

### 5.2 Specific Heat

| Phase | Symbol | Value | Unit |
|---|---|---|---|
| Solid | Cp,m | 415 + 0.1838·T | J/kg/K (T-dependent) |
| Liquid | Cp,l | 830 | J/kg/K (constant) |
| Gas | Cp,g | 1006.43 | J/kg/K (constant) |

### 5.3 Thermal Conductivity

| Phase | Symbol | Value | Unit |
|---|---|---|---|
| Solid | k_s | 9.23 + 0.0139·T | W/m/K (T-dependent) |
| Liquid | k_l | 5.5 + 0.0133·T | W/m/K (T-dependent) |
| Gas | k_g | 0.02 | W/m/K (constant) |

**Evaluated at key temperatures (for tabulation in solvers like OpenFOAM):**

| T (K) | k_s (W/m/K) | Cp,s (J/kg/K) | ρ_l (kg/m³) |
|---|---|---|---|
| 300 | 13.40 | 470.1 | 7969 |
| 500 | 16.18 | 506.9 | 7815 |
| 800 | 20.35 | 562.0 | 7584 |
| 1200 | 25.91 | 635.6 | 7276 |
| 1658 | 32.26 | 719.7 | 6924 |
| 1723 (liq) | k_l = 28.42 | 830 (liq) | 6876 |
| 2000 (liq) | k_l = 32.10 | 830 (liq) | 6660 |
| 2500 (liq) | k_l = 38.75 | 830 (liq) | 6275 |
| 3090 (vap) | k_l = 46.60 | 830 (liq) | 5822 |

### 5.4 Phase Change Temperatures

| Property | Symbol | Value | Unit |
|---|---|---|---|
| Solidus temperature | T_s | 1658 | K |
| Liquidus temperature | T_l | 1723 | K |
| Evaporation temperature | T_v | 3090 | K |

### 5.5 Latent Heats

| Property | Symbol | Value | Unit |
|---|---|---|---|
| Latent heat of melting | L_m | 2.6 × 10⁵ | J/kg |
| Latent heat of vaporization | L_v | 7.45 × 10⁶ | J/kg |

### 5.6 Transport & Other Properties

| Property | Symbol | Value | Unit |
|---|---|---|---|
| Viscosity, metallic phase | μ_m | 0.006 | N·s/m² (Pa·s) |
| Viscosity, gas phase | μ_g | 1.85 × 10⁻⁵ | N·s/m² (Pa·s) |
| Molar mass | M | 0.05593 | kg/mol |
| Universal gas constant | R | 8.314 | J/mol/K |
| Stefan–Boltzmann constant | σ_s | 5.67 × 10⁻⁸ | W/m²/K⁴ |
| Emissivity | ε | 0.5 | — |
| Ambient temperature | T₀ | 300 | K |
| Ambient pressure | P₀ | 1.013 × 10⁵ | Pa |

### 5.7 Model Coefficients

| Parameter | Symbol | Value | Unit | Notes |
|---|---|---|---|---|
| Laser absorptivity | α_laser | 0.5 | — | Bare plate value |
| Recoil pressure coefficient | A_recoil | 0.6 | — | |
| Evaporation energy-loss coeff. | A_evap | 0.5 | — | |
| Mushy zone constant | A_mushy | 1 × 10⁹ | kg/m³/s | Selected for small mushy zone; literature range 10³–10¹⁵ |
| Laser spot radius | r_laser | 50 | μm | |

---

## 6. Computational Model Setup

### 6.1 Solver

| Item | Detail |
|---|---|
| Software | ANSYS Fluent 2022 R1 |
| Algorithm | SIMPLE (pressure-based segregated) |
| Flow assumption | Newtonian, incompressible |
| Turbulence | Standard k-ε model |
| Free surface | Volume of Fluid (VOF) |

### 6.2 Simulation Domains

| Domain | Size (μm³) | Mesh size | Purpose |
|---|---|---|---|
| Single-track (core) | 1500 × 90 × 270 | 3.5 μm | High-res melt pool zone |
| Single-track (overall) | 1500 × 250 × 637.5 | 3.5–10 μm | Full domain with graded mesh |
| Multi-track | 2900 × 800 × 530 | 8 μm | Surface roughness prediction |

### 6.3 Governing Equations Summary

| Equation | Conservation law |
|---|---|
| ∂ρ/∂t + ∇·(ρu) = 0 | Mass |
| ∂(ρu)/∂t + ∇·(ρuu) = −∇p + ∇·(μ_eff(∇u+∇uᵀ)) + ρg + F_source | Momentum |
| ∂(ρH)/∂t + ∇·(ρuH) = ∇·(k_eff∇T) + q_source | Energy |
| ∂α_m/∂t + ∇·(α_m·u) = 0 | VOF (volume fraction) |

### 6.4 Momentum Source Terms

| Term | Physical origin |
|---|---|
| F_surface tension = σκn̂\|∇α_m\|(2ρ/(ρ_m+ρ_g)) | Surface tension (normal) |
| F_marangoni = [∇σ − (n̂·∇σ)n̂]\|∇α_m\|(2ρ/(ρ_m+ρ_g)) | Marangoni (tangential) |
| F_recoil = A_recoil·P₀·exp[L_v·M(T−T_v)/(RT·T_v)]·n̂\|∇α_m\|(2ρ/(ρ_m+ρ_g)) | Recoil pressure |
| F_damping = [(1−γ_l)²/(γ_l³+0.001)]·A_mushy·(u−V_laser) | Mushy zone damping (Carman-Kozeny) |

### 6.5 Energy Source Terms

| Term | Physical origin |
|---|---|
| q_laser = α_laser·(0.5·P/πr²)·exp[−(r/2r_laser)²]·\|∇α_m\|·(2ρ/(ρ_m+ρ_g)) | Gaussian laser input (surface deposition) |
| q_evap = −A_evap·L_v·√(M/2πRT)·P₀·exp[L_v·M(T−T_v)/(RT·T_v)]·\|∇α_m\|·(2ρ/(ρ_m+ρ_g)) | Evaporation energy loss |
| q_radiation = −σ_s·ε·(T⁴−T₀⁴)·\|∇α_m\|·(2ρ/(ρ_m+ρ_g)) | Radiation loss |

### 6.6 Liquid Fraction (Mushy Zone)

| Condition | γ_l |
|---|---|
| T ≤ T_s (1658 K) | 0 |
| T_s < T < T_l | (T − T_s)/(T_l − T_s) |
| T ≥ T_l (1723 K) | 1 |

---

## 7. Surface Tension Model

Surface tension σ is **nonlinear with temperature** (Fig. 2 in paper — not a simple linear fit). Key values approximately:

| T (K) | σ (N/m) approx. |
|---|---|
| 1723 (liquidus) | ~1.75 |
| 2000 | ~1.60 |
| 2500 | ~1.35 |
| 3090 (vap) | ~1.00 |

> The Marangoni term dσ/dT is negative (surface tension decreases with temperature), driving flow from hot center to cool periphery. The exact σ(T) curve should be read from Fig. 2 of the paper; a simple linear approximation will introduce error particularly at high temperatures relevant to keyhole mode.

---

## 8. Comparison Notes vs. Other Sources

| Property | Zhang et al. (2024) | King et al. (2014) | Abdali et al. (2024) |
|---|---|---|---|
| Configuration | Bare plate + powder, LPBF | Powder bed only, LPBF | Bare plate, DED + SLM |
| Absorptivity | 0.5 (bare), higher (powder) | 0.4 (powder) | Not tabulated (qualitative) |
| L_v (J/kg) | 7.45 × 10⁶ | Not tabulated | Not tabulated |
| T_vap (K) | 3090 | Not tabulated | Not tabulated |
| ρ_s (kg/m³) | 7950 | 7980 | Not tabulated |
| D (m²/s) | 5.38 × 10⁻⁶ (from King) | 5.38 × 10⁻⁶ | Not tabulated |
| h_s at melting (J/kg) | 2.6 × 10⁵ (L_m) | 1.2 × 10⁶ | Not tabulated |

> Note: King's h_s = 1.2 × 10⁶ J/kg is the enthalpy at melting (sensible + latent heat from 300 K to T_liquidus), not just the latent heat. Zhang's L_m = 2.6 × 10⁵ J/kg is the latent heat of fusion only. These are not directly comparable.

---

## 9. Reproduction Checklist

### What this paper gives you directly
- Full process parameter matrix (9 cases × bare + powder)
- Complete temperature-dependent thermophysical properties (ρ, k, c_p for solid, liquid, gas)
- All model coefficients (α_laser, A_recoil, A_evap, A_mushy, ε)
- Phase-change temperatures (T_s, T_l, T_v) and latent heats (L_m, L_v)
- Simulation domain sizes and mesh resolution
- Governing equations and all source term formulations
- Melt pool width and depth measurements vs. process parameters (Figs. 5–7)
- Mode maps (conduction vs. keyhole) for bare and powder plate
- Surface morphology observations (flat, swell-undercut, balling, humping)

### What you must handle separately
- Exact σ(T) curve: read from Fig. 2 of the paper (nonlinear; cannot be reduced to a single dσ/dT)
- Powder absorption model: they reference Gusarov et al. (2009) for powder-layer laser heating — needed only if reproducing powder-plate cases
- Quantitative melt pool dimension tables: the exact measured W and D values per case are presented in Figs. 5–7 as plots, not as tables — values must be digitized from the figures
- Turbulence model constants: standard k-ε defaults assumed; not stated explicitly

### Key flags when comparing to your OpenFOAM file
| Property | Zhang et al. | Your OpenFOAM file | Flag |
|---|---|---|---|
| k_l (liquid) | 5.5 + 0.0133T (rises with T) | Flat 26.90 W/m/K above liquidus | Diverges at high T — check keyhole depth |
| L_v (J/kg) | 7.45 × 10⁶ | 6.33 × 10⁶ | ~15% lower in your file — affects recoil |
| T_vap (K) | 3090 | 3122 | Minor (~1%) |
| ε (emissivity) | 0.5 | 0.3 (emS/emL) | 40% lower in your file — less radiation loss |
| α_laser (bare) | 0.5 | Not shown in file shown | Confirm in your laser BC |
| A_mushy | 1 × 10⁹ | Not shown in file shown | Confirm in your solver settings |