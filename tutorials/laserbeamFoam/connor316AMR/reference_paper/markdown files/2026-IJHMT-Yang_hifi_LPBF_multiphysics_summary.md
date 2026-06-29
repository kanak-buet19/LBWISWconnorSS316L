# Paper Summary

**Title:** High-fidelity multiphysics modeling of laser powder bed fusion: Coupling laser absorption, vaporization, and powder heat transfer
**Authors:** Jian Yang, Akash Aggarwal, Christian Leinenbach
**Journal:** International Journal of Heat and Mass Transfer 265 (2026) 128801
**DOI:** [https://doi.org/10.1016/j.ijheatmasstransfer.2026.128801](https://doi.org/10.1016/j.ijheatmasstransfer.2026.128801)

---

## Process Parameters — Validation Cases

### Ti-6Al-4V bare plate (Cunningham et al. 2019)

| Case | Power (W) | Speed (mm/s) | Spot w(1/e²) (µm) | Mode |
|------|-----------|--------------|-------------------|------|
| Stationary | 156 | 0 | 70 | Deep unstable keyhole |
| Fast scan | 416 | 1500 | 47.5 | Stable elongated depression |
| Scan A | 230 | 400 | 47.5 | Unstable keyhole |
| Scan B | 208 | 600 | 47.5 | Stable keyhole |

- Wavelength: 1070 nm
- Defocus: adjusted per case to match experimental spot size

### Ti-6Al-4V LPBF with powder bed (Simonds et al. 2021)

| Power (W) | Speed (mm/s) | Spot w(1/e²) (µm) | Powder layer |
|-----------|--------------|-------------------|--------------|
| 201 | 700 | 61.25 | 50 µm (default, not specified in exp) |

- Defocusing: −22 mm below focal plane → spot diameter (1/e²) = 122.5 µm at surface

### SS316L bare plate (Yang et al. 2024, TOMCAT PSI synchrotron)

| Case | Power (W) | Speed (mm/s) | Spot w(1/e²) (µm) |
|------|-----------|--------------|-------------------|
| Sample e | 85 | 50 | 71 |
| Remelting | 70 | 50 | 30 |

- Remelting case: initial substrate temperature = 800 K (10 s dwell between scans)

---

## Geometry

| Parameter | Value |
|-----------|-------|
| Standard CFD domain | 1200 × 300 × 600 µm (scan × width × depth) |
| Substrate height | 400 µm |
| Powder bed thickness | 50 µm |
| Argon gas layer | 200 µm (above powder surface) |
| DEM domain (powder generation) | 1200 × 300 × 400 µm |
| Bare plate Ti64 validation | substrate = 430 µm; no powder |
| Boundary — top/bottom (y) | adiabatic walls |
| Boundary — left/right (z) | adiabatic walls |
| Boundary — front/back (x) | natural convection (swak4Foam), h calibrated to 300 K ambient |
| Argon gas velocity | 0.5 m/s along −z (scan direction) |
| Argon pressure | 1 atm |

---

## Mesh

| Parameter | Value |
|-----------|-------|
| Type | Structured, locally refined |
| Minimum cell size (near melt pool / powder) | 5 µm |
| Mesh strategy | Adaptive local refinement around powder bed; coarser elsewhere |
| Powder bed initialization | Continuous α field from DEM particle geometry (not binary 0/1) |

---

## Timestep

Not explicitly stated. Minimum cell size = 5 µm with Courant-limited explicit scheme implied; sub-µs steps expected near keyhole.

---

## Powder Bed Properties (DEM — LIGGGHTS)

| Parameter | Value |
|-----------|-------|
| Particle size distribution | Log-normal; mean diameter 35 µm, std 15 µm |
| Particle size range | 7.5 – 57.5 µm (5 µm bins) |
| Material (DEM) | SS316L surrogate: density 8000 kg/m³, Young's modulus 196 GPa, Poisson ratio 0.3 |
| Friction coefficient | 0.2 |
| Restitution coefficient | 0.3 |
| Contact model | Johnson–Kendall–Roberts (SJKR) cohesive |
| k_p (unmelted powder conductivity) | Calibrated to match homogenized packed-powder model (Wei et al. 2018) |
| Effective kappa mixing rule | k = α⁴ kₚ + (1 − α⁴) k_gas (non-linear; suppresses over-conduction) |
| Absorptivity enhancement factor | 1.2 – 1.4× above smooth-surface Fresnel (calibrated at P=1 W) |

---

## Key Model Parameters

### Vaporization (multi-element, Hertz–Langmuir)

| Parameter | Value |
|-----------|-------|
| Retro-diffusion coefficient β_R | 0.18 |
| Elements considered | Ti, Al, V (Ti6Al4V); Fe, Cr, Ni (SS316L) |
| Recoil pressure smoothing radius | r_T = 800 K rolling-circle |
| C_scale (recoil smoothing) | 100 kg m⁻¹ s⁻² K⁻¹ |

### Vapor plume attenuation & scattering

| Parameter | Value |
|-----------|-------|
| Max attenuation coefficient µ_alpha^max | 310 m⁻¹ |
| Max plume height H_plume | 1 mm |
| Max scanning speed for plume model u_max | 1.5 m/s |
| Max blurring constant C_max | 0.25 |
| Attenuation law | Beer–Lambert |

---

## Thermophysical Properties

Temperature-dependent (piecewise functions) for all phases. Shown in Appendix Fig. A6; not tabulated in main text. Key features:

**Ti-6Al-4V:**
- Density, cp, kappa: piecewise T-dependent (solid/mushy/liquid)
- Surface tension: linearly decreasing with T
- Absorptivity: modified Fresnel, slightly increasing with T until boiling point

**SS316L:**
- Density, cp, kappa: piecewise T-dependent
- Surface tension: non-monotone (increases then decreases) due to oxygen/sulfur surface-active elements → temperature-dependent Marangoni (CST) is critical
- Absorptivity: temperature-dependent, from Ebrahimi et al. (2022)

---

## Validation Cases

| # | Material | Configuration | Exp. reference | Metrics compared |
|---|----------|---------------|----------------|-----------------|
| 1 | Ti-6Al-4V | Bare plate stationary, 156 W, 0 mm/s | Cunningham et al. (2019) | Melt pool shape, keyhole depth vs time, absorptivity |
| 2 | Ti-6Al-4V | Bare plate scan, 416 W, 1500 mm/s | Cunningham et al. (2019) | Melt pool shape, absorptivity |
| 3 | Ti-6Al-4V | Bare plate scan, 230 W/400 mm/s and 208 W/600 mm/s | Cunningham et al. (2019) | Keyhole morphology (stable/unstable) |
| 4 | Ti-6Al-4V | Single-track LPBF with powder, 201 W, 700 mm/s | Simonds et al. (2021) | Keyhole morphology + laser absorption vs time |
| 5 | SS316L | Bare plate scan (sample e), 85 W, 50 mm/s | Yang et al. (2024) TOMCAT PSI | Melt pool morphology at 2/3/4 ms |
| 6 | SS316L | Bare plate remelting, 70 W, 50 mm/s | Yang et al. (2024) TOMCAT PSI | Keyhole formation and development |

**Sensitivity study (Table 1):** ρ, Cₚ, k each +10% → ~5–10% change in melt pool size. Neglecting T-dependent absorptivity → −22% keyhole depth. Neglecting vapor effects → +62% keyhole depth. Neglecting T-dependent CST (SS316L) → +70% keyhole depth.
