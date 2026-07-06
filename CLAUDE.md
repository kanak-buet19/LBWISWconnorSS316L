# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Environment setup

```bash
of2412          # source OpenFOAM v2412 (also v2506, v2512 supported)
venv314         # activate Python venv (needed for analysis/calibration scripts)
```

## Build

```bash
./Allwmake -j               # build all libs then all solvers
./Allwmake -j -s            # with reduced output
./Allwclean                 # clean everything

# Build only one layer
(cd src && ./Allwmake -j)          # libraries only
(cd applications && ./Allwmake -j) # solvers only
```

## Run and test

```bash
# Run a tutorial case
cd tutorials/laserbeamFoam/connor316AMR/some_case
./Allrun

# Run all tutorials as smoke test (single-timestep)
cd tutorials && ./Alltest

# Calibration framework (316L SS Bayesian optimization)
cd tutorials/laserbeamFoam/connor316AMR/calibrate_CW_case
./Allrun                          # local run
sbatch job_calibration.sh         # HPC (Slurm + Apptainer)
CALIB_TOTAL_CORES=16 CALIB_CORES_PER_SIM=4 ./Allrun   # override parallelism
```

## Architecture

Three solver applications + five custom libraries + Python calibration framework.

### Solvers (`applications/solvers/`)

| Solver | Type | Key feature |
|---|---|---|
| `laserbeamFoam` | Incompressible, two-phase VOF + AMR | Primary solver: ray-tracing heat source, melting/solidification (enthalpy-porosity), Marangoni, recoil pressure, AMR via `dynamicRefineFvMesh` |
| `compressibleLaserbeamFoam` | Compressible, multi-component | Vaporisation with mass transfer, multi-phase mixture thermo |
| `laserMeltFoam` | Incompressible, three-phase | Legacy solver; isoAdvector only, solid/liquid/gas phases |

### Libraries (`src/`, built in order)

1. **`geometricVoF`** → `libgeometricVoFLMFOAM.so` — isoAdvector + PLIC interface reconstruction
2. **`laserHeatSource`** → `liblaserHeatSource.so` — Lagrangian ray-tracing with Fresnel absorption; core of the laser model
3. **`MTHD`** → `libMTHD.so` — magneto-thermo-hydro-dynamic model (Lorentz force, Joule heating)
4. **`transportModels`** → `libincompressibleTransportModelsLMFOAM.so` — custom viscosity + phase mixture models
5. **`turbulenceModel`** → `libincompressibleTurbulenceModelLMFOAM.so` — runtime-selectable turbulence (laminar only; hook for RAS/LES)

### Solver code structure (all solvers follow this split-header pattern)

```
laserbeamFoam.C       ← main loop (PIMPLE + sub-cycles)
createFields.H        ← field registration and initialisation
UEqn.H                ← momentum equation (Darcy, Marangoni, recoil, Lorentz)
pEqn.H                ← pressure equation (surface tension, recoil, buoyancy)
TEqn.H                ← energy equation with latent heat iteration (epsilon1)
updateProps.H         ← material properties + AMR indicator computation
updateKappaCp.H       ← temperature-dependent cp/kappa via tables/polynomials
```

### AMR implementation

Indicator computed in `updateProps.H`, driven by `dynamicRefineFvMesh`. Four criteria selectable via `constant/dynamicMeshDict`:
- `interfaceGradT` — refine at interface where |grad T| > threshold
- `interfaceOnly` — refine at interface + 1 dilation layer
- `interfaceT100` / `interfaceT300` — refine at interface where T > (Tsolidus - deltaT)

### Calibration framework (`connor316AMR/calibrate_CW_case/`)

Bayesian optimization (Optuna TPE) calibrates 12 thermophysical parameters against Hofmann (2026) CW melt-pool measurements. See `calibrate_CW_case/CLAUDE.md` for full architecture. Key files:
- `calibration_config.json` — all parameter ranges, optimizer settings, geometry
- `scripts/run_calibration.py` — main orchestrator (`Calibrator` + `Job` model)
- `scripts/caselib.py` — patches `constant/transportProperties` per candidate via regex
- `scripts/evaluate.py` — melt-pool stability scoring from per-timestep CSV

## Key conventions

- **Runtime selection tables**: `TypeName("...")` in header, `defineTypeNameAndDebug(...)`, `addToRunTimeSelectionTable(...)` in source
- **Dictionary-driven config**: all physics parameters read from `constant/` dictionaries, never hardcoded
- **Laser config**: `constant/LaserProperties` with `timeVsLaserPosition` and `timeVsLaserPower` tables
- **Ray-tracing**: keep in `src/laserHeatSource/`, not in solver loops
- **Make/files**: every `.C` file must be listed; update when adding/removing source files
- See `.agents/skills/laserbeamfoam/SKILL.md` for detailed C++ style and OpenFOAM convention rules

## Adding a runtime-selectable model

1. Add `TypeName("...")` in header
2. Register with `addToRunTimeSelectionTable(...)` in source
3. Add `.C` to the module's `Make/files`
4. Ensure dictionary `type` string matches the registered class name
