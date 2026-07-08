# testrun19 / testrun21 — case summary

Paired comparison study: same laser/geometry/mesh setup, only surface-tension
(Marangoni) sign/magnitude differs, to test effect of sulfur (PLC) content on
melt-pool flow direction.

| | testrun19 | testrun21 |
|---|---|---|
| `sigma` (N/m) | 1.8648 | 1.421 |
| `Marangoni_Constant` (N/m/K) | -2.28e-4 | +1.244e-3 |
| Reference condition | pure 316L, η=0 | 316L + PLC, η=0.1 |
| `numberOfSubdomains` | 8 (scotch) | 32 (hierarchical 2×2×8) |

Negative Marangoni_Constant (testrun19) → surface tension falls with T →
outward (center-to-edge) surface flow, shallow/wide pool.
Positive (testrun21) → surface tension rises with T → inward flow,
deep/narrow pool (typical surfactant/sulfur reversal effect).

## Solver / physics

- Solver: `laserbeamFoam` (incompressible VOF + AMR), run via case `Allrun`
  (not the `application` field in controlDict, which is stale:
  `Flint_multiphaseEulerFoamD` — leftover, unused).
- Phases: `metal` (316L SS), `gas` (ambient).
- Interface tracking: isoAdvector.

### Material properties (`constant/transportProperties`)

**metal (316L SS)**
- ν = 5e-7 m²/s, ρ = 8000 kg/m³
- Tsolidus = 1658 K, Tliquidus = 1723 K, LatentHeat = 2.7e5 J/kg
- β (thermal exp.) = 5e-6 1/K
- κ(T) = 11.5 + 0.00307·T W/m/K (poly_kappa)
- cp = 500 J/kg/K (const)
- elec_resistivity = 1.6384e-7 Ω·m → ~10% absorptivity at normal incidence
  (Drude–Fresnel model)

**gas**
- ν = 1.48e-5 m²/s, ρ = 1 kg/m³, κ = 0.04 W/m/K, cp = 520 J/kg/K
- Tsolidus/Tliquidus = 1 K / 10 K (dummy — keeps gas always "liquid" so
  solidification physics never triggers there)

**Shared (both cases)**
- p0 = 1e5 Pa, Tvap = 3068 K, Mm = 0.05 kg/mol, LatentHeatVap = 7.45e6 J/kg
- Emissivity: solid 0.4, liquid 0.1
- TRef = 300 K

### Laser (`constant/LaserProperties`, `timeVsLaserPosition/Power`)

- Wavelength = 1.064 µm (Yb fiber laser)
- Beam radius = 3.5e-5 m (35 µm), Radius_Flavour = 1.336 (M² = 1.497)
- Power = const 147.3684 W over full run
- Scan: x fixed at 0, y fixed at 10 µm (in the gas headspace above the free
  surface), z moves 100 µm → 400 µm over 1.5 ms
  → scan speed = 300 µm / 1.5 ms = **0.2 m/s (200 mm/s)**
- Volumetric energy density: VED = P / (v · π r²)
  = 147.3684 / (200 · π · 0.035²) = **191.5 J/mm³** (same for both cases —
  P, v, r identical)

### Geometry / mesh (`system/blockMeshDict`)

- Domain: 0.32 mm (x) × 0.5 mm (y) × 0.5 mm (z), `convertToMeters 0.001`
- Base mesh: 32×50×50 cells → 10 µm cubic cells
- y=0 face = `atmosphere` (free surface, laser incidence side);
  y=0.5 mm face = `lowerWall` (substrate base)
- Initial `alpha.metal`: metal fills y ∈ [100, 500] µm; y ∈ [0,100] µm is gas
  headspace above the free surface (room for keyhole/vapor depression)

### AMR (`constant/dynamicMeshDict`)

- `dynamicRefineFvMesh`, refine on `alpha.metal` (generic interface AMR, not
  one of the custom `interfaceGradT`/`interfaceT100`/`interfaceT300` criteria)
- Refine band: alpha ∈ [0.001, 0.999], unrefine below 0.001
- maxRefinement = 1 → 10 µm → 5 µm at the interface
- refineInterval = 5, nBufferLayers = 1, maxCells = 200000

### Time control (`system/controlDict`)

- endTime = 1.5 ms, initial deltaT = 1e-8 s, adjustTimeStep on
- maxCo = 0.1, maxAlphaCo = 0.1, maxDeltaT = 2e-6 s
- writeInterval = 1e-5 s (adjustableRunTime)

### Boundary conditions

- T: `frontAndBack` fixedGradient (uniform -50 K/m, ambient heat loss);
  `lowerWall`/`atmosphere` zeroGradient
- Laser_boundary: `frontAndBack` fixedValue -1 (marks domain edges out of
  laser search region)

### Known gaps fixed this session

- `constant/trackProperties` was missing from both cases (solver requires
  it) — added, `trackDuration` set to 0.0015 s to match `endTime` and the
  `timeVsLaserPosition` table end time.
- `foamVTK.sh` (reconstruct + `foamToVTK`) was missing — added standard
  version from `template/`.
