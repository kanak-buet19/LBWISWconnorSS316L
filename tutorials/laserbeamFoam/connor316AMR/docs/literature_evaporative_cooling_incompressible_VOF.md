# Evaporative Cooling in Incompressible VOF Solvers for LPBF/Laser Welding

## Problem Statement

In incompressible VOF laser melting CFD (laserbeamFoam), temperature overshoots past the boiling
point (T_vap ≈ 3122 K for 316L → T_max ~5000 K). Root cause:

- Laser flux at interface: ~2.4×10¹⁰ W/m² (200 W, 25 µm radius)
- Evaporative cooling Qv at T_vap: ~3.1×10⁸ W/m²
- Ratio: ~78× gap — T must rise to ~4600–5000 K for the exponential in Qv to balance
- Qv is interface-restricted via `gradAlphaDamper` (CSF approach) — no cooling in bulk liquid
- Result: CFL crash (velocity spike → CoNum spike → Δt → 1e-8), unphysical melt pool

## Approaches in Literature

### 1. Hertz-Knudsen / Anisimov Evaporative Cooling (Current laserbeamFoam Approach)

**Used by**: Flint et al. (2023), Khairallah et al. (2016), Klassen et al. (2014)

**Formulation**:
```
Qv = 0.82 * Lv * M * P0 * exp(Lv*M*(T - Tvap) / (R*T*Tvap)) / sqrt(2π*M*R*T)
```
Applied at interface only via `|∇α|` (CSF method).

**Limitation**: Exponential saturates too slowly for high-power-density LPBF. Accommodation
coefficient 0.82 is for low-speed evaporation — at keyhole-relevant fluxes, the Knudsen layer
departs from this equilibrium assumption.

**References**:
- Flint TF, Robson JD, Parivendhan G, Cardiff P. *laserbeamFoam: Laser ray-tracing and thermally
  induced state transition simulation toolkit.* SoftwareX, 2023.
  DOI: [10.1016/j.softx.2022.101299](https://doi.org/10.1016/j.softx.2022.101299)
- Klassen A, Scharowsky T, Körner C. *Evaporation model for beam based additive manufacturing
  using free surface lattice Boltzmann methods.* J Phys D: Appl Phys, 2014.
  DOI: [10.1088/0022-3727/47/27/275303](https://doi.org/10.1088/0022-3727/47/27/275303)
- Anisimov SI, Khokhlov VA. *Instabilities in Laser-Matter Interaction.* CRC Press, 1995.
  (Original recoil pressure derivation)

---

### 2. Lee Model — Volumetric Phase Change (Most Common Alternative)

**Used by**: OpenFOAM `interCondensatingEvaporatingFoam`, `icoReactingMultiPhaseInterFoam`,
  Bunaziv et al. (2025), Tao & Zhao (2022, 2023)

**Formulation** (evaporation, T > T_sat):
```
m_dot_evap = r * ρ_l * α_l * (T - T_sat) / T_sat
```
Energy sink: `S_E = -m_dot_evap * Lv`

**Key difference from current approach**:
- Acts **volumetrically** wherever α_l > 0 and T > T_sat — not just at interface
- Relaxation form: drives T back toward T_sat at rate proportional to r
- Coefficient `r` can be physically derived from Hertz-Knudsen, not just tuned

**Physical derivation of r** (from kinetic theory):
```
r = (6 * β / d) * sqrt(M / (2π*R*T_sat)) * Lv * (ρ_v * ρ_l / (ρ_l - ρ_v))
```
where β = accommodation coefficient, d = characteristic interface thickness.

**Key insight**: Larger `r` → tighter T clamping near T_sat. Values 10³–10⁶ common for
laser welding after accounting for 1/T_sat scaling in OpenFOAM's implementation.

**Why relevant**: The Lee model is the **simplest volumetric generalization** of the
interface-only Hertz-Knudsen approach. It converts "interface flux" → "volumetric
relaxation" without requiring a separate vapor phase.

**References**:
- Chen G, Nie T, Yan X. *An explicit expression of the empirical factor in a widely used
  phase change model.* Int J Heat Mass Transfer, 2020.
  DOI: [10.1016/j.ijheatmasstransfer.2019.119279](https://doi.org/10.1016/j.ijheatmasstransfer.2019.119279)
- Lee WH. *A pressure iteration scheme for two-phase flow modeling.* In: Veziroglu TN (ed)
  Multiphase Transport Fundamentals, Reactor Safety, Applications. Hemisphere, 1980.
  (Original Lee model)
- Bunaziv I et al. *Numerical modelling of phase changes in laser materials processing using
  OpenFOAM.* IOP Conf Ser: Mater Sci Eng, 2025.
  DOI: [10.1088/1757-899X/1332/1/012006](https://doi.org/10.1088/1757-899X/1332/1/012006)

---

### 3. Extended Enthalpy-Porosity with Vaporization Latent Heat

**Used by**: Verhaeghe et al. (2009), Courtois et al. (2014)

**Concept**: Extend the enthalpy-porosity framework (currently handles melting: T_solidus →
T_liquidus) to also absorb vaporization latent heat above T_boil.

**Implementation**:
```
h(T) = ∫ cp dT + f_l(T) * L_fusion + f_v(T) * L_vaporization
```
where f_v(T) transitions from 0 → 1 over [T_boil - ΔT, T_boil + ΔT].

**Energy equation**: `∂(ρh)/∂t + ∇·(ρuh) = ∇·(k∇T) + S_laser`

**Advantage**: No explicit mass transfer needed. T cannot exceed T_boil + ΔT because all
excess energy goes into latent heat of vaporization. Physics-based: energy conservation is
exact. Smooth formulation with existing enthalpy-porosity iteration.

**Limitation**: No mass loss accounting (vapor mass removed, not tracked). Recoil pressure
still needs separate model. The smearing width ΔT is a numerical parameter — too small = stiff,
too large = unphysical.

**References**:
- Verhaeghe F, Craeghs T, Heulens J, Pandelaers L. *A pragmatic model for selective laser
  melting with evaporation.* Acta Materialia, 2009.
  DOI: [10.1016/j.actamat.2009.08.027](https://doi.org/10.1016/j.actamat.2009.08.027)
  (Key finding: evaporation consumes ~4× more enthalpy than heating from room temp to fully
  molten for Ti6Al4V — cannot be neglected)
- Courtois M, Carin M, Le Masson P, Gaied S, Balabane M. *A complete model of keyhole and
  melt pool dynamics to analyze instabilities and collapse during laser welding.*
  J Laser Applications, 2014.
  DOI: [10.2351/1.4886835](https://doi.org/10.2351/1.4886835)

---

### 4. Tao & Zhao High-Fidelity Framework (2022) — Gold Standard for OpenFOAM VOF

**Reference paper**:
- Yu T, Zhao J. *Quantitative simulation of selective laser melting of metals enabled by new
  high-fidelity multiphase, multiphysics computational tool.*
  Comput Methods Appl Mech Eng, 2022.
  DOI: [10.1016/j.cma.2022.115422](https://doi.org/10.1016/j.cma.2022.115422)
  (Open Access PDF: [jzhao.people.ust.hk](http://jzhao.people.ust.hk/home/PDFs/2022-CMAME-Tao.pdf))

**Key features**:
- Dual interface scheme: isoAdvector (melt pool, CFL<1) + MULES (vapor jet, CFL>1)
- Knudsen layer evaporation model in VOF
- Ray tracing with Fresnel absorption
- Full energy equation: laser + conduction + convection + radiation + fusion + vaporization
- Vaporization sink term: `0.82 * ρ₀*Lv*M / sqrt(2πMRT) * exp(Lv*M*(T-T_LV)/(R*T*T_LV)) * |∇α| * 2Cρ/(C₁ρ₁+C₂ρ₂)`
- CFD-DEM coupling for powder dynamics
- Validated against synchrotron X-ray experiments

**Takeaway**: Even the gold-standard framework uses the same interface-restricted Qv
formulation. The difference is: (a) dual VOF schemes handle vapor dynamics better,
(b) better mesh resolution at interface, (c) CFD-DEM captures powder effects that
modify energy coupling.

**Open-source solver** (simplified, no DEM):
- [`thermocapillaryInterFoam`](https://github.com/pzimbrod/thermocapillaryInterFoam) —
  VOF + laser + Marangoni + melting/evaporation + AMR
- [`lpbfFoam`](https://github.com/hdubbs/lpbfFoam) — fork for PBF

---

### 5. Parameter-Scaled CSF — Fixing Interface Resolution Error

**Reference**:
- Schreter-Fleischhacker M et al. *A consistent diffuse-interface finite element approach to
  rapid melt–vapor dynamics in metal additive manufacturing.*
  Comput Methods Appl Mech Eng, 2025.
  DOI: [10.1016/j.cma.2025.117700](https://doi.org/10.1016/j.cma.2025.117700)
  (arXiv: [2501.18781](https://arxiv.org/abs/2501.18781))

**Key finding**: Errors in recoil pressure are ~10× larger than errors in interface temperature
due to exponential amplification. Classical CSF (Brackbill) underpredicts evaporative cooling
at under-resolved interfaces. A parameter-scaled CSF reduces temperature error by one order of
magnitude.

**Implication**: Part of laserbeamFoam's T overshoot may be numerical — the `gradAlphaDamper`
smearing across 2-3 cells (~12-18 µm) under-resolves the interface for accurate Qv computation.
Sharper interface (isoAdvector vs MULES) or narrower smearing would help, but not solve the
fundamental energy balance gap.

---

### 6. Compressible Mass-of-Fluid (MoF) — Implicit Recoil from Phase Change

**Reference**:
- Zenz C, Buttazzoni M, Florian T, Crespo Armijos KE, Gómez Vázquez R, Liedl G, Otto A.
  *A compressible multiphase Mass-of-Fluid model for the simulation of laser-based
  manufacturing processes.* Computers & Fluids, 2024.
  DOI: [10.1016/j.compfluid.2023.106109](https://doi.org/10.1016/j.compfluid.2023.106109)

**Key innovation**: Evaporation recoil pressure is NOT explicitly modeled — it emerges from
compressible phase change + Tait EOS. Phase mass tracked directly (not volume fraction),
ensuring mass conservation. No empirical recoil coefficient (0.54) needed.

**Relevance**: User explicitly ruled out compressible solver. But the conceptual insight
matters: the `mass_dot` volumetric source term in the energy equation IS the correct
physics-based sink — it's the same concept as the Lee model but with full vapor dynamics.
The question is whether we can use the energy sink without the compressible vapor phase.

---

### 7. Khairallah et al. (2016) — ALE3D Benchmark

**Reference**:
- Khairallah SA, Anderson AT, Rubenchik A, King WE. *Laser powder-bed fusion additive
  manufacturing: Physics of complex melt flow and formation mechanisms of pores, spatter,
  and denudation zones.* Acta Materialia, 2016.
  DOI: [10.1016/j.actamat.2015.12.004](https://doi.org/10.1016/j.actamat.2015.12.004)

**Approach**: ALE3D (LLNL): ray-tracing + recoil pressure + Marangoni + evaporation cooling.
Uses recoil pressure from Clausius-Clapeyron with surface-temperature-dependent Psat.
Showed that recoil pressure dominates melt flow in keyhole regime. Evaporative cooling is
modeled as surface heat flux (not volumetric).

**Limitation noted**: Gas/vapor phase not explicitly modeled — later addressed by Tao & Zhao
and the compressible MoF group.

---

## Summary: Viable Physics-Based Approaches for laserbeamFoam

| Approach | Physics Basis | Implementation Complexity | Requires Vapor Phase? | Solves T Overshoot? |
|---|---|---|---|---|
| **1. Lee volumetric relaxation** | Kinetic theory → linearized HK | Low (add source term to TEqn) | No | Yes (drives T→T_sat) |
| **2. Extended enthalpy-porosity** | Latent heat absorption above T_boil | Low (modify cp/T table) | No | Yes (caps T via enthalpy) |
| **3. Enhanced Qv coefficient** | Accommodation coeff tuning | Minimal | No | Partial (only at interface) |
| **4. Sharper interface (isoAdvector)** | Better-resolved interface → accurate Qv | Medium | No | Partial (numerical fix) |
| **5. Full multiphase (mass_dot)** | Compressible vapor phase | High | Yes | Yes (complete physics) |

### Recommendation

**Approach 1 (Lee model)** is the best fit for the user's constraint ("port mass_dot concept
without vapor model"). It:
- Is physically derived from the same Hertz-Knudsen kinetics as the current Qv
- Acts volumetrically (solves the "only-at-interface" limitation)
- Does NOT require a separate vapor phase
- Has a tunable coefficient with physical bounds (can be derived from kinetic theory)
- Is already implemented in OpenFOAM's phase change framework as reference
- Naturally relaxes T → T_sat with rate proportional to superheating

**Approach 2 (extended enthalpy-porosity)** is the simplest fallback if Approach 1 causes
convergence issues. It reuses the existing epsilon1 iteration machinery.

**Approach 3 (enhanced coefficient)** is worth doing regardless — the 0.82 accommodation
coefficient assumes low-speed evaporation. At LPBF-relevant temperatures, the Knudsen layer
is in a transition regime where effective accommodation may differ.
