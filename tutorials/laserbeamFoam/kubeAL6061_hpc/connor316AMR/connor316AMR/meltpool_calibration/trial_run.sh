#!/bin/bash

set -Eeo pipefail

fail()
{
    echo "FAIL: $*" >&2
    exit 1
}

check()
{
    echo
    echo "=== $* ==="
}

scriptDir="$(cd "$(dirname "$0")" && pwd)"
caseDir="$(cd "$scriptDir/.." && pwd)"
cd "$scriptDir"

OF2506_IMAGE="${OF2506_IMAGE:-$HOME/openfoam-dev_2506.sif}"
OF2506_USER="${OF2506_USER:-x-rkanak1}"
PYTHON="${PYTHON:-/home/x-rkanak1/.conda/envs/isw_env/bin/python3}"

of2506()
{
    apptainer exec --cleanenv --env USER="$OF2506_USER" "$OF2506_IMAGE" bash -lc "source /openfoam/bash.rc && $*"
}

check "Paths"
echo "trial_run.sh: $0"
echo "meltpool_calibration: $scriptDir"
echo "base case: $caseDir"
echo "OF2506_IMAGE: $OF2506_IMAGE"
echo "PYTHON: $PYTHON"

[ -d "$scriptDir" ] || fail "missing meltpool_calibration directory"
[ -d "$caseDir/system" ] || fail "missing base case system directory"
[ -d "$caseDir/constant" ] || fail "missing base case constant directory"
[ -d "$caseDir/initial" ] || fail "missing base case initial directory"
[ -f "$scriptDir/calibrate_meltpool.py" ] || fail "missing calibrate_meltpool.py"
[ -f "$scriptDir/fine_zone_config.json" ] || fail "missing fine_zone_config.json"
[ -x "$scriptDir/calibrate_fine_job.sh" ] || fail "calibrate_fine_job.sh not executable"
[ -x "$caseDir/Allrun" ] || fail "base case Allrun not executable"
[ -x "$caseDir/Allrun_long" ] || fail "base case Allrun_long not executable"

check "Scheduler tools"
if command -v sbatch >/dev/null 2>&1; then
    echo "sbatch: $(command -v sbatch)"
else
    echo "WARN: sbatch not found in this shell"
fi
if command -v squeue >/dev/null 2>&1; then
    echo "squeue: $(command -v squeue)"
else
    echo "WARN: squeue not found in this shell"
fi

check "Apptainer"
command -v apptainer || fail "apptainer not found"
[ -f "$OF2506_IMAGE" ] || fail "OpenFOAM 2506 image not found: $OF2506_IMAGE"
apptainer --version

check "OpenFOAM 2506 container"
of2506 '
set -e
echo "WM_PROJECT_DIR=$WM_PROJECT_DIR"
echo "WM_PROJECT_VERSION=${WM_PROJECT_VERSION:-unset}"
echo "FOAM_APPBIN=${FOAM_APPBIN:-unset}"
echo "FOAM_USER_APPBIN=${FOAM_USER_APPBIN:-unset}"
command -v mpirun
mpirun --version | head -3
command -v blockMesh
command -v setFields
command -v decomposePar
command -v foamDictionary
command -v reconstructPar
command -v foamToVTK
command -v laserbeamFoam
'

check "Python environment"
[ -x "$PYTHON" ] || fail "PYTHON not executable: $PYTHON"
"$PYTHON" - <<'PY'
import importlib
import sys

required = ["matplotlib", "numpy", "pandas", "pyvista", "scipy", "vtk"]
missing = []

for module in required:
    try:
        imported = importlib.import_module(module)
        version = getattr(imported, "__version__", "unknown")
        print(f"{module}: OK ({version})")
    except Exception as exc:
        missing.append(f"{module}: {exc}")

if missing:
    print("Missing modules:")
    for item in missing:
        print(f"  - {item}")
    sys.exit(1)
PY

check "Calibration config"
"$PYTHON" - <<'PY'
import json
from pathlib import Path

config = json.loads(Path("fine_zone_config.json").read_text())
required = [
    "target_width_um",
    "target_depth_um",
    "first_pulse_end_s",
    "simulation_end_s",
    "write_interval_s",
    "laser_radius_um",
    "elec_resistivity_ohm_m",
]
missing = [key for key in required if key not in config]
if missing:
    raise SystemExit(f"missing fine_zone_config keys: {missing}")

radii = config["laser_radius_um"]
resistivities = config["elec_resistivity_ohm_m"]
grid_cases = len(radii) * len(resistivities)
adaptive = config.get("adaptive", {})
print(f"radius candidates: {radii}")
print(f"resistivity candidates: {resistivities}")
print(f"initial grid cases: {grid_cases}")
print(f"config max_parallel_cases: {config.get('max_parallel_cases')}")
print(f"adaptive enabled: {adaptive.get('enabled', False)}")
print(f"adaptive max_new_cases: {adaptive.get('max_new_cases', 0)}")
PY

check "MPI rank plan"
subdomains="$(of2506 "foamDictionary -entry numberOfSubdomains -value '$caseDir/system/decomposeParDict'")"
echo "numberOfSubdomains: $subdomains"
[ "$subdomains" = "8" ] || fail "expected numberOfSubdomains 8, got $subdomains"
echo "calibrate_fine_job.sh requests: 64 Slurm tasks"
echo "calibrate_fine_job.sh runs: 8 concurrent cases"
echo "planned rank use: 8 cases * $subdomains ranks = 64 ranks"

check "Case preparation dry run"
"$PYTHON" calibrate_meltpool.py \
    --config fine_zone_config.json \
    --max-cases 1 \
    --no-adaptive
firstCase="$(find "$scriptDir/runs" -maxdepth 1 -mindepth 1 -type d ! -name '_trial_preflight_case' | sort | head -1)"
[ -n "$firstCase" ] || fail "calibrate_meltpool.py did not prepare any run case"
echo "prepared case: $firstCase"
[ -f "$firstCase/system/controlDict" ] || fail "prepared case missing controlDict"
[ -f "$firstCase/constant/LaserProperties" ] || fail "prepared case missing LaserProperties"
[ -f "$firstCase/constant/transportProperties" ] || fail "prepared case missing transportProperties"

check "Container can read prepared case"
of2506 "cd '$firstCase' && foamDictionary -entry endTime -value system/controlDict && foamDictionary -entry writeInterval -value system/controlDict && foamDictionary -entry numberOfSubdomains -value system/decomposeParDict"

check "Result"
echo "PASS: environment looks ready for:"
echo "  sbatch calibrate_fine_job.sh"
echo
echo "Note: this preflight prepares one calibration case but does not run solver."
