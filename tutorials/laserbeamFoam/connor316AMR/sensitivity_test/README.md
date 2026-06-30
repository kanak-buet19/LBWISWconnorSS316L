# 316L Transport Sensitivity Test

This folder runs a short-track sensitivity study for melt-pool width using the
local `template/` case structure with
`validation/materials/SS316L/transportProperties` overlaid as the material
baseline. Beam diameter, laser power, and scan speed are held fixed; only
transport/material parameters are perturbed.

Study setup:

- base mesh: 40 um, `32 x 32 x 32` cells over a `1.28 mm` cube
- AMR: `maxRefinement 1`
- scan track: `250 um` centered in the domain
- scan speed: `0.9 m/s`
- run time: `0.0002777778 s`
- cooling: convection and radiation disabled for the quick forced-heating comparison
- outputs: final-time VTU analysis, CSV summary, and `results/dashboard.html`

Run from this directory after sourcing OpenFOAM:

```bash
source /openfoam/bash.rc
./Allrun
```

`Allrun` uses `${PYTHON}` when set, otherwise `python3` from the active shell.

Useful partial commands:

```bash
./scripts/generate_cases.py
./scripts/run_sensitivity.py --skip-existing
./scripts/build_dashboard.py
./Allclean
```

The dashboard ranks parameters by the largest absolute melt-pool width shift
from baseline across the low/high perturbations.
