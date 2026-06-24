# Gan et al. 2025 — AMB2022 Ti-6Al-4V Notes

**DOI**: 10.1007/s40964-024-00637-6  
**Journal**: Progress in Additive Manufacturing (2025) 10:491–515  
**Title**: Benchmark study of melt pool and keyhole dynamics, laser absorptance, and porosity in additive manufacturing of Ti-6Al-4V

## Experimental setup

- Laser: 1070 nm Yb-doped fiber, Gaussian, **7° angle of incidence** from normal
- Spot size: **122.5 µm diameter** (1/e²); beam waist 49.5 µm ± 5 µm
- Substrate: Ti-6Al-4V NIST SRM 654b, ~300 µm thick, polished
- Environment: vacuum + argon backfill
- Measurement: APS 32-ID-B synchrotron X-ray (in situ), integrating sphere absorptance
- Laser pulse duration: ~2 ms

## Cases validated in this paper

| Config | P (W) | SS (mm/s) | Substrate |
|--------|-------|-----------|-----------|
| Stationary | 172 | 0 | bare |
| Stationary | 201 | 0 | bare |
| Stationary | 201 | 0 | powder (100 µm) |
| Scanning | 261 | 700 | bare |
| Scanning | 201 | 700 | powder |
| Scanning | 261 | 700 | powder |

> These are a **subset** of the full AMB2022 dataset (Simonds et al. [34]). Original dataset covers wider power range.

## Material properties (Table 1)

| Property | Value |
|----------|-------|
| ρ_s (kg/m³) | 4575 − 0.168T |
| ρ_L (kg/m³) | 4616.5 − 0.254T |
| μ (Pa·s) | exp(−1.6 + 5346/T) × 10⁻³ |
| C_p,s (J/kg·K) | 483.04 − 0.22T (T < 1268 K); 412.70 − 0.18T (1268–1923 K) |
| C_p,L (J/kg·K) | 831 |
| κ_s (W/m·K) | −0.32 + 1.46×10⁻²T |
| κ_L (W/m·K) | −0.66 + 1.83×10⁻²T |
| T_solidus (K) | 1877 |
| T_liquidus (K) | 1923 |
| T_boil (K) | 3533 |
| ΔH_fusion (kJ/kg) | 286 |
| ΔH_evap (kJ/kg) | 9830 |
| h_conv (W/m²·K) | 25 |
| Fresnel ε | 0.15 |

## Laser / absorption model

- Gaussian beam: Q = (2P / π r²_b) exp(−2r²/r²_b)
- Fresnel absorptance with ε = 0.15 (electrical conductance constant)
- Recoil pressure: P_recoil = 0.54 P_a exp(ΔH_lv/R_v T_boil · (1 − T_boil/T))
- Clausius–Clapeyron evaporative cooling

## Simulation domains (FLOW-3D, 7 µm cubical cells)

| Config | x (µm) | y (µm) | z (µm) | Cells |
|--------|---------|---------|---------|-------|
| Stationary | 600 | 400 | 1560 | ~1.09 M |
| Scanning | 1700 | 400 | 740 | ~1.47 M |

Scanning: laser moves along +x from (50, 200) µm to (1650, 200) µm.

## Key results

- Absorptivity: rises rapidly < 0.1 ms → steady ~0.6–0.85 (keyhole mode)
- Stationary: sim overestimates depth, width more accurate
- Powder layer reduces melt pool depth + width vs bare (absorbs/disperses energy before reaching solid)
- Keyhole aspect ratio scales linearly with NEI (R² ≈ 0.95–0.99)
- Keyhole angle → ~90° as NEI increases (stationary); for scanning, angle increases with keyhole number Ke
- Pore size decreases with larger spot size at fixed power (more distributed energy)

## Dimensionless numbers used

NEI = ηP / (κ r (T_1 − T_0))  — normalized energy input (stationary)

Ke = ηP / ((T_1 − T_0) π ρ C_p √(α V_s) r²_0)  — keyhole number (scanning)
