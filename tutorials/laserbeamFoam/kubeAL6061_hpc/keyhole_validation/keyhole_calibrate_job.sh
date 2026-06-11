#!/bin/bash
#SBATCH --job-name=kubeKeyhole
#SBATCH --time=48:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=64
#SBATCH --account=mch250110
#SBATCH --output=keyhole_%j.out
#SBATCH --error=keyhole_%j.err

set -Eeo pipefail

on_error()
{
    local status="$?"
    local line="${BASH_LINENO[0]:-unknown}"
    local command="${BASH_COMMAND:-unknown}"

    echo
    echo "ERROR: keyhole validation job failed"
    echo "Exit status: $status"
    echo "Line: $line"
    echo "Command: $command"
    echo "Working directory: $PWD"
    echo "Time: $(date)"
    echo

    exit "$status"
}

trap on_error ERR

if [ -n "${SLURM_SUBMIT_DIR:-}" ]; then
    cd "$SLURM_SUBMIT_DIR"
fi

caseDir="$(pwd)"

OF2506_IMAGE="${OF2506_IMAGE:-$HOME/openfoam-dev_2506.sif}"
OF2506_USER="${OF2506_USER:-x-rkanak1}"
PYTHON="${PYTHON:-/home/x-rkanak1/.conda/envs/isw_env/bin/python3}"

export FOAM_SIGFPE="${FOAM_SIGFPE:-0}"
export MPLCONFIGDIR="${MPLCONFIGDIR:-/tmp}"
export PYTHON
KEYHOLE_CONFIG="${KEYHOLE_CONFIG:-config.json}"
case "$KEYHOLE_CONFIG" in
    /*) ;;
    *) KEYHOLE_CONFIG="$caseDir/$KEYHOLE_CONFIG" ;;
esac
export KEYHOLE_CONFIG
export RANKS_PER_CASE="${RANKS_PER_CASE:-8}"
export MAX_PARALLEL_CASES="${MAX_PARALLEL_CASES:-8}"

runKeyholeInsideContainer()
{
    source /openfoam/bash.rc

    echo "=== OpenFOAM environment ==="
    echo "WM_PROJECT_VERSION=${WM_PROJECT_VERSION:-unset}"
    echo "FOAM_USER_APPBIN=${FOAM_USER_APPBIN:-unset}"
    command -v laserbeamFoam
    command -v mpirun
    echo "============================"

    "$PYTHON" - <<'PY'
import importlib
import sys

required = ["matplotlib", "numpy", "pandas", "pyvista", "scipy", "sklearn", "vtk"]
missing = []
for module in required:
    try:
        importlib.import_module(module)
    except Exception as exc:
        missing.append(f"{module}: {exc}")

if missing:
    print("Missing Python modules:")
    for item in missing:
        print("  -", item)
    sys.exit(1)

print("Python environment OK")
PY

    "$PYTHON" scripts/keyhole_search.py \
        --config "$KEYHOLE_CONFIG" \
        --mode all \
        --run \
        --force \
        --jobs "$MAX_PARALLEL_CASES" \
        --ranks-per-case "$RANKS_PER_CASE"
}

if [ "${KEYHOLE_CALIBRATE_IN_CONTAINER:-0}" = "1" ]; then
    runKeyholeInsideContainer
    exit 0
fi

echo "=== Kube keyhole validation job ==="
echo "Case directory: $caseDir"
echo "Job ID: ${SLURM_JOB_ID:-unknown}"
echo "Node list: ${SLURM_NODELIST:-unknown}"
echo "Tasks: ${SLURM_NTASKS:-unknown}"
echo "OF2506_IMAGE: $OF2506_IMAGE"
echo "PYTHON: $PYTHON"
echo "Config: $KEYHOLE_CONFIG"
echo "MAX_PARALLEL_CASES: $MAX_PARALLEL_CASES"
echo "RANKS_PER_CASE: $RANKS_PER_CASE"
echo "Started: $(date)"
echo "==================================="

command -v apptainer
test -f "$OF2506_IMAGE"
test -x "$PYTHON"
test -f ./scripts/keyhole_search.py
test -f "$KEYHOLE_CONFIG"

apptainer exec \
    --cleanenv \
    --env USER="$OF2506_USER" \
    --env PYTHON="$PYTHON" \
    --env MPLCONFIGDIR="$MPLCONFIGDIR" \
    --env FOAM_SIGFPE="$FOAM_SIGFPE" \
    --env KEYHOLE_CONFIG="$KEYHOLE_CONFIG" \
    --env RANKS_PER_CASE="$RANKS_PER_CASE" \
    --env MAX_PARALLEL_CASES="$MAX_PARALLEL_CASES" \
    --env KEYHOLE_CALIBRATE_IN_CONTAINER=1 \
    "$OF2506_IMAGE" \
    bash "$caseDir/keyhole_calibrate_job.sh"

echo "Finished: $(date)"
