#!/bin/bash
#SBATCH --job-name=paramSweep
#SBATCH --time=48:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=32
#SBATCH --account=mch250110
#SBATCH --output=param_sweep_%j.out
#SBATCH --error=param_sweep_%j.err

# Slurm driver for the CW parametric sweep.
#
# Submit:
#   sbatch job.sh
#
# Override sizing:
#   sbatch --ntasks-per-node=40 --export=ALL,PARAM_CORES_PER_SIM=8 job.sh
#
# Dry run:
#   sbatch --export=ALL,DRY_RUN=1 job.sh

set -Eeo pipefail

on_error()
{
    local status="$?"
    echo
    echo "ERROR: parametric_sweep/job.sh failed (exit $status) at line ${BASH_LINENO[0]:-?}"
    echo "Command: ${BASH_COMMAND:-?}"
    echo "Working directory: $PWD"
    echo "Time: $(date)"
    exit "$status"
}
trap on_error ERR

if [ -n "${SLURM_SUBMIT_DIR:-}" ]; then
    cd "$SLURM_SUBMIT_DIR"
else
    cd "$(dirname "$0")"
fi
caseDir="$(pwd)"

export OF2506_IMAGE="${OF2506_IMAGE:-$HOME/openfoam-dev_2506.sif}"
export OF2506_USER="${OF2506_USER:-x-rkanak1}"
export PYTHON="${PYTHON:-/home/x-rkanak1/.conda/envs/isw_env/bin/python3}"
export PARAM_TOTAL_CORES="${PARAM_TOTAL_CORES:-${SLURM_NTASKS:-32}}"
export PARAM_CORES_PER_SIM="${PARAM_CORES_PER_SIM:-8}"
export FOAM_SIGFPE=0
export MPLCONFIGDIR="$caseDir/.mplconfig"
mkdir -p "$MPLCONFIGDIR"

echo "=== Slurm CW parametric sweep bootstrap ==="
echo "Sweep dir:           $caseDir"
echo "Job ID:              ${SLURM_JOB_ID:-unknown}"
echo "Node list:           ${SLURM_NODELIST:-unknown}"
echo "Slurm tasks:         ${SLURM_NTASKS:-unset}"
echo "PARAM_TOTAL_CORES:   $PARAM_TOTAL_CORES"
echo "PARAM_CORES_PER_SIM: $PARAM_CORES_PER_SIM"
echo "OF2506_IMAGE:        $OF2506_IMAGE"
echo "PYTHON:              $PYTHON"
echo "==========================================="

if [ "$PARAM_TOTAL_CORES" -lt "$PARAM_CORES_PER_SIM" ]; then
    echo "ERROR: PARAM_TOTAL_CORES must be >= PARAM_CORES_PER_SIM" >&2
    exit 2
fi

echo "=== Environment Check ==="
command -v apptainer
test -f "$OF2506_IMAGE"
"$PYTHON" - <<'PY'
import importlib
import sys

required = ["matplotlib", "numpy", "pandas", "pyvista", "scipy", "vtk"]
missing = []
for module in required:
    try:
        importlib.import_module(module)
    except Exception as exc:
        missing.append(f"{module}: {exc}")

if missing:
    print("ERROR: missing Python modules needed by post-processing:")
    for item in missing:
        print(f"  - {item}")
    sys.exit(1)

print("Host Python environment OK.")
PY
apptainer exec --cleanenv \
    --env USER="$OF2506_USER" \
    --env PYTHON="$PYTHON" \
    --env PARAM_TOTAL_CORES="$PARAM_TOTAL_CORES" \
    --env PARAM_CORES_PER_SIM="$PARAM_CORES_PER_SIM" \
    --env FOAM_SIGFPE=0 \
    --env MPLCONFIGDIR=/tmp \
    "$OF2506_IMAGE" \
    bash -lc 'source /openfoam/bash.rc &&
command -v mpirun &&
command -v blockMesh &&
command -v setFields &&
command -v decomposePar &&
command -v foamDictionary &&
command -v reconstructPar &&
command -v foamToVTK &&
command -v laserbeamFoam &&
"$PYTHON" - <<'"'"'PY'"'"'
import importlib
import sys

required = ["matplotlib", "numpy", "pandas", "pyvista", "scipy", "vtk"]
missing = []
for module in required:
    try:
        importlib.import_module(module)
    except Exception as exc:
        missing.append(f"{module}: {exc}")

if missing:
    print("ERROR: host Python is not usable inside Apptainer:")
    for item in missing:
        print(f"  - {item}")
    sys.exit(1)

print("Container OpenFOAM and analysis Python OK.")
PY'
echo "========================="

if [ "${DRY_RUN:-0}" = "1" ]; then
    echo "=== DRY RUN ==="
    echo "Would run parametric sweep from $caseDir"
    exit 0
fi

echo "Parametric sweep started at: $(date)"
apptainer exec --cleanenv \
    --env USER="$OF2506_USER" \
    --env PYTHON="$PYTHON" \
    --env PARAM_TOTAL_CORES="$PARAM_TOTAL_CORES" \
    --env PARAM_CORES_PER_SIM="$PARAM_CORES_PER_SIM" \
    --env FOAM_SIGFPE=0 \
    --env MPLCONFIGDIR=/tmp \
    --env DELETE_ANALYZED_VTK="${DELETE_ANALYZED_VTK:-true}" \
    "$OF2506_IMAGE" \
    bash -lc 'source /openfoam/bash.rc && cd "$1" && shift && ./Allrun "$@"' \
        bash "$caseDir" "$@"
echo "Parametric sweep completed at: $(date)"
