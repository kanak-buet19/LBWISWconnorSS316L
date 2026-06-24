# Keyhole Pore Formation Simulation Parameters
**Material:** Ti-6Al-4V | **Condition:** Atmospheric Pressure (1 atm)
**Reference:** Wang et al., npj Computational Materials (2022) 8:22

---

## 1. Laser Process Parameters

| Parameter | Case 2 | Case 3 | Case 4 |
|---|---|---|---|
| Laser Power (W) | 382 | 382 | 382 |
| Scanning Speed (mm/s) | 525 | 500 | 475 |
| Laser Spot Diameter (μm) | 100 | 100 | 100 |
| Ambient Pressure (atm) | 1.0 | 1.0 | 1.0 |
| Substrate Type | Bare plate | Bare plate | Bare plate |

> **Note:** Case 1 (205 W, 500 mm/s, 1 atm) was used only for instant bubble formation study and is excluded here as it represents a different experimental condition.

---

## 2. Simulation Domain & Numerical Settings

| Parameter | Value |
|---|---|
| Mesh Size | 4 μm |
| Physical Simulation Time | 2000 μs |
| Laser Beam Profile | Parallel (constant diameter) |
| Free Surface Method | Volume of Fluid (VoF) |
| Flow Assumption | Incompressible Newtonian, Laminar |
| Drag Force Model | Blake-Kozeny (Darcy) |
| Primary Dendrite Arm Spacing (λ₁) | 5 μm |
| Darcy Coefficient Constant (180μ/ρλ₁²) | 5.57 × 10⁶ |

---

## 3. Melt Pool & Keyhole Geometry (Simulation Results)

| Parameter | Case 2 (525 mm/s) | Case 3 (500 mm/s) | Case 4 (475 mm/s) |
|---|---|---|---|
| Mean Keyhole Depth — Sim (μm) | 345 | 372 | 388 |
| Mean Keyhole Depth — Exp (μm) | 355 | 388 | 456 |
| Std. Dev. Keyhole Depth — Sim (μm) | 30 | 29 | 35 |
| Std. Dev. Keyhole Depth — Exp (μm) | 10 | 14 | 17 |
| Max Keyhole Depth Fluctuation (μm) | 152 | 170 | 156 |
| Mean Keyhole Pore Size (μm) | 37 | 49 | 82 |
| Max Keyhole Pore Size — Sim (μm) | 76 | ~50 | ~120 |

### Energy Absorptivity on Keyhole Surface

| Parameter | Case 2 | Case 3 | Case 4 |
|---|---|---|---|
| Front Wall — Mean (%) | 25.54 | 26.07 | 27.07 |
| Front Wall — Std. Dev. (%) | 2.29 | 2.51 | 2.63 |
| Rear Wall — Mean (%) | 50.04 | 49.88 | 49.04 |
| Rear Wall — Std. Dev. (%) | 3.02 | 3.26 | 4.21 |

---

## 4. Thermophysical Properties of Ti-6Al-4V

| Property | Symbol | Value | Unit |
|---|---|---|---|
| Solidus Temperature | Tₛ | 1878 | K |
| Liquidus Temperature | T_l | 1928 | K |
| Boiling Temperature | T_b | 3315 | K |
| Density (solidus) | ρ | 4400 | kg/m³ |
| Latent Heat of Melting | L_m | 2.86 × 10⁵ | J/kg |
| Latent Heat of Evaporation | L_v | 9.70 × 10⁶ | J/kg |
| Saturated Vapour Pressure | P_e | 1.013 × 10⁵ | Pa (at T_b = 3315 K) |
| Specific Heat (solidus) | c_s | 570 | J/K/kg |
| Specific Heat (liquidus) | c_l | 831 | J/K/kg |
| Thermal Conductivity (solidus) | k_s | 16 | W/m/K |
| Thermal Conductivity (liquidus) | k_l | 32 | W/m/K |
| Surface Radiation Coefficient | ε | 0.4 | — |
| Surface Tension Coefficient | σ₀ | 1.68 | N/m |
| Temperature Sensitivity of σ | σᵀₛ | 2.6 × 10⁻⁴ | N/m/K |
| Dynamic Viscosity | μ | 0.005 | Pa·s |

---

## 5. Physics & Boundary Conditions

### Governing Physics Included
- Heat transfer (conduction, convection, radiation, evaporation loss)
- Molten pool fluid flow (incompressible Newtonian, laminar)
- Marangoni effect (surface tension gradient-driven flow)
- Recoil pressure from metal evaporation
- Darcy drag force in mushy zone (Blake-Kozeny model)
- Laser ray-tracing with Fresnel equation (multi-reflection)
- Volume of Fluid (VoF) free surface tracking
- Boussinesq approximation for buoyancy

### Thermal Boundary Condition
$$-k\nabla T \cdot n = h(T - T_{env}) + \varepsilon\sigma_s(T^4 - T^4_{env}) + \dot{q}_{evp}$$

| Constant | Value |
|---|---|
| Stefan-Boltzmann constant (σₛ) | 5.6704 × 10⁻⁸ W/m²/K⁴ |
| Reference temperature (T_ref) | T_l = 1928 K |

### Bubble Pressure (Adiabatic after formation)
$$p = p_0\left(\frac{V_0}{V}\right)^\gamma$$

---

## 6. Key Porosity Outcomes Summary

| Case | Speed (mm/s) | Mean Pore Size (μm) | Pore Shape | Distribution |
|---|---|---|---|---|
| Case 2 | 525 | 37 | Spherical | Horizontal |
| Case 3 | 500 | 49 | Slightly irregular | Horizontal |
| Case 4 | 475 | **82** | Irregular, flat/sharp bottom | Non-horizontal |

> **Worst case:** Case 4 (475 mm/s) — lowest speed produces largest, most irregular pores with highest keyhole instability (rear wall absorptivity std. dev. = 4.21%).