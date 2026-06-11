#!/bin/bash
#SBATCH --job-name=connorFineCal
#SBATCH --time=12:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=64
#SBATCH --account=mch250110
#SBATCH --output=output_%j.out
#SBATCH --error=output_%j.err

set -Eeo pipefail

on_error()
{
    local status="$?"
    local line="${BASH_LINENO[0]:-unknown}"
    local command="${BASH_COMMAND:-unknown}"

    echo
    echo "ERROR: calibrate_fine_job.sh failed"
    echo "Exit status: $status"
    echo "Line: $line"
    echo "Command: $command"
    echo "Working directory: $PWD"
    echo "Time: $(date)"
    echo

    exit "$status"
}

trap on_error ERR

echo "=== Slurm fine calibration bootstrap ==="
echo "Script: $0"
echo "Submit directory: ${SLURM_SUBMIT_DIR:-unset}"
echo "Launch directory: $PWD"
echo "Job ID: ${SLURM_JOB_ID:-unknown}"
echo "Node list: ${SLURM_NODELIST:-unknown}"
echo "Number of tasks: ${SLURM_NTASKS:-unknown}"
echo "Forwarded calibration args: $*"
echo "==========================="

if [ -n "${SLURM_SUBMIT_DIR:-}" ]; then
    cd "$SLURM_SUBMIT_DIR"
fi

scriptDir="$(pwd)"
cd "$scriptDir"

export MPLCONFIGDIR="$scriptDir/.mplconfig"
mkdir -p "$MPLCONFIGDIR"

OF2506_IMAGE="${OF2506_IMAGE:-$HOME/openfoam-dev_2506.sif}"
OF2506_USER="${OF2506_USER:-x-rkanak1}"

of2506()
{
    apptainer exec --cleanenv --env USER="$OF2506_USER" "$OF2506_IMAGE" bash -lc "source /openfoam/bash.rc && $*"
}

export FOAM_SIGFPE=0
export PYTHON="${PYTHON:-/home/x-rkanak1/.conda/envs/isw_env/bin/python3}"

echo "=== Environment Check ==="
command -v apptainer
test -f "$OF2506_IMAGE"
echo "PYTHON=$PYTHON"
echo "OF2506_IMAGE=$OF2506_IMAGE"
of2506 '
echo "WM_PROJECT_DIR=$WM_PROJECT_DIR"
echo "WM_PROJECT_VERSION=${WM_PROJECT_VERSION:-unset}"
echo "FOAM_USER_APPBIN=${FOAM_USER_APPBIN:-unset}"
echo "FOAM_APPBIN=${FOAM_APPBIN:-unset}"
command -v mpirun
mpirun --version
command -v blockMesh
command -v setFields
command -v decomposePar
command -v foamDictionary
command -v laserbeamFoam
'
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
    print("ERROR: missing Python modules needed by calibration/post-processing:")
    for item in missing:
        print(f"  - {item}")
    sys.exit(1)

print("Python environment OK: all calibration/post-processing modules imported.")
PY
echo "========================="

echo "Fine calibration started at: $(date)"
echo "Running 8 cases at a time with 8 MPI ranks per case, using up to 64 ranks total."

if [ "${DRY_RUN:-0}" = "1" ]; then
    echo "=== DRY RUN ==="
    echo "scriptDir=$scriptDir"
    echo "MPLCONFIGDIR=$MPLCONFIGDIR"
    echo "SLURM_SUBMIT_DIR=${SLURM_SUBMIT_DIR:-unset}"
    echo "SLURM_NTASKS=${SLURM_NTASKS:-unset}"
    echo "Would run: $PYTHON calibrate_meltpool.py --config fine_zone_config.json --run --jobs 8 --run-command <apptainer Allrun>"
    echo "DRY RUN PASS: no calibration cases started."
    exit 0
fi

runCommand="apptainer exec --cleanenv --env USER=\"$OF2506_USER\" \"$OF2506_IMAGE\" bash -lc 'source /openfoam/bash.rc && export FOAM_SIGFPE=0 && export PYTHON=$PYTHON && export MPLCONFIGDIR=/tmp && ./Allrun'"

"$PYTHON" calibrate_meltpool.py \
    --config fine_zone_config.json \
    --run \
    --jobs 8 \
    --run-command "$runCommand" \
    "$@"

echo "Fine calibration completed at: $(date)"
