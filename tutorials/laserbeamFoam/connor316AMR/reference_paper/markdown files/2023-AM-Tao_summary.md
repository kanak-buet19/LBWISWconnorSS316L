# Paper Summary

**Title:** Quantifying the mechanisms of keyhole pore evolutions and the role of metal-vapor condensation in laser powder bed fusion
**Authors:** Tao Yu, Jidong Zhao
**Journal:** Additive Manufacturing 72 (2023) 103642
**DOI:** [https://doi.org/10.1016/j.addma.2023.103642](https://doi.org/10.1016/j.addma.2023.103642)

---

## Process Parameters

| Parameter | Value |
|-----------|-------|
| Material | Ti-6Al-4V bare plate |
| Laser power | 205 W |
| Scanning speed | 500 mm/s |
| Laser spot radius | 100 µm |
| Laser model | VOF-compatible ray-tracing with Fresnel absorption |
| Ambient gas | Argon (treated as single ambient gas phase) |

---

## Geometry

| Parameter | Value |
|-----------|-------|
| Domain | 720 × 400 × 320 µm (x × y × z) |
| Scanning direction | x (longest axis, 720 µm) |
| Width direction | y (400 µm) |
| Depth direction | z (320 µm) |
| Note | Fig. 1A shows domain cut at y = 200 µm for visualization only; full domain simulated |
| Boundary — bottom | No-slip wall; zero gradient for pressure and temperature |
| Boundary — all others | Fixed value for pressure; zero gradient for velocity and temperature |

---

## Mesh

| Parameter | Value |
|-----------|-------|
| Coarse base mesh size | 20 µm |
| Fine (dense) mesh size | 5 µm |
| Refinement strategy | Dynamic meshing; regions with \|∇α₁\| > 0, T > 1000 K, and Sᵢ > 0 are refined each timestep |
| α₁ | Volume fraction of metal phase |
| Sᵢ | Absorbed laser energy source term |

---

## Timestep

| Parameter | Value |
|-----------|-------|
| Time step | 1 × 10⁻⁷ s (100 ns) |

---

## Thermophysical Properties — Ti-6Al-4V

### Constant parameters

| Parameter | Value |
|-----------|-------|
| Room temperature T₀ | 300 K |
| Solidus temperature Tₛ | 1878 K |
| Liquidus temperature Tₗ | 1923 K |
| Boiling temperature T_LV | 3133 K |
| Molar mass M | 446.07 g/mol |
| Viscosity — liquid µₗ | 0.05 Pa·s |
| Viscosity — solid µₛ | 1.13 Pa·s |
| Latent heat of fusion Lf | 2.88 × 10⁵ m²/s² |
| Latent heat of evaporation Lv | 10⁷ m²/s² |
| Permeability coefficient Kc | 5.56 × 10⁶ kg/(m³·s) |
| Constant to avoid div-by-zero Cₖ | 10⁻³ |
| Surface tension at liquidus σᵢ | 1.68 N/m (kg/s²) |
| dσ/dT (Marangoni coefficient) | −2.6 × 10⁻⁴ N/(m·K) |
| Refractive index e | 3.47 |
| Electrical conductance coefficient ε | 0.2 |
| Reflectivity at room temperature f_R0 | 0.95 |
| Reflectivity at liquidus f_Rl | 0.63 |
| Absorption coefficient at room temperature γ₀ | 47.1 µm⁻¹ |
| Absorption coefficient at liquidus γ₀ | 0.192 µm⁻¹ |
| Convective heat transfer coefficient h | 19 kg·s⁻³·K |

### Ambient gas properties

| Parameter | Value |
|-----------|-------|
| Gas density ρ₂ | 1.87 kg/m³ |
| Gas viscosity µ₂ | 2.5 × 10⁻⁵ Pa·s |
| Gas heat capacity C₂ | 500 J/(kg·K) |
| Gas thermal conductivity k₂ | 0.021 W/(m·K) |

### Temperature-dependent properties

**Density ρ₁ (kg/m³):**
```
4420                          T < 1268 K
4420 − 0.154(T − 298)         1268 K ≤ T < 1923 K
3920 − 0.680(T − 1923)        T ≥ 1923 K
```

**Heat capacity C₁ (J/(kg·K)):**
```
411.5                              T < 1268 K
411.5 + 0.2T + 5×10⁻⁷ T²         1268 K ≤ T < 1923 K
830                                T ≥ 1923 K
```

**Thermal conductivity k₁ (W/(m·K)):**
```
~19.0                                         T < 1268 K
−0.80 + 0.0018T − 2×10⁻⁸ T²                 1268 K ≤ T < 1923 K
33.4                                           1923 K ≤ T < 1973 K
34.6                                           T ≥ 1973 K
```

---

## Simulation Case — Stationary Laser (Own Study)

### Process Parameters

| Parameter | Value |
|-----------|-------|
| Material | Ti-6Al-4V |
| Laser power | 156 W |
| Laser type | Stationary (no scanning) |
| Laser spot diameter | 140 µm (radius 70 µm) |
| Laser on | 0 → 2 ms |
| Cooling period | 2 ms → 2.5 ms (laser off) |
| Total simulation time | 2.5 ms |

### Geometry

| Parameter | Value |
|-----------|-------|
| Domain | 400 × 400 × 350 µm (x × y × z) |
| x | 400 µm |
| y (width) | 400 µm |
| z (depth) | 350 µm |

### Mesh

| Parameter | Value |
|-----------|-------|
| Base mesh size | 20 µm |
| Fine mesh size | 5 µm |
| Refinement strategy | Dynamic (same AMR criteria as Tao et al.) |

---

## Validation Cases

| Case | Description | Reference |
|------|-------------|-----------|
| Benchmark I | Single-phase Marangoni-driven cavity flow; isotherms and Nusselt number vs. analytical/numerical solutions | Bergman (1988), Salid (2012), Sen & Davis (1982) |
| Benchmark II | Two-phase Marangoni flow with free surface; dimensionless surface heights at left/right walls | Sasmal & Hochstein (1994), Francois et al. (2006), Saldi (2012) |
| Benchmark III | Free-surface Marangoni flow with phase change; liquid-solid interface shape | Salid (2012), Tan et al. (2006) |
| Benchmark IV | Stationary laser keyhole depth evolution; 140 µm spot, 364 W power | Cunningham et al. (2019) [Ref. 18] |
| Primary validation | Moving keyhole dynamics (J-like keyhole formation, pore formation, collapse, splitting); qualitative and quantitative shape/speed comparison | Zhao et al. (2020) megahertz X-ray imaging [Ref. 13] |
