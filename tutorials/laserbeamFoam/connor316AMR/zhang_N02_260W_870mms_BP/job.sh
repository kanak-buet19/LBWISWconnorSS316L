#!/bin/bash
#SBATCH --job-name=zhN02_260W_870mms
#SBATCH --time=24:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=16
#SBATCH --account=mch250110
#SBATCH --output=output_%j.out
#SBATCH --error=output_%j.err

# Zhang N02: 260 W / 0.87 m/s — conduction mode, bare plate 316L
# Domain: 800x500x1200 um, 20 um base -> 5 um AMR (2 levels)
# Expected runtime: ~1.7 ms sim time

set -Eeo pipefail

on_error()
{
    local status="$?"
    local line="${BASH_LINENO[0]:-unknown}"
    local command="${BASH_COMMAND:-unknown}"

    echo
    echo "ERROR: job.sh failed"
    echo "Exit status: $status"
    echo "Line: $line"
    echo "Command: $command"
    echo "Working directory: $PWD"
    echo "Time: $(date)"
    echo

    exit "$status"
}

trap on_error ERR

echo "=== Slurm zhang_N01_260W_520mms_BP bootstrap ==="
echo "Script: $0"
echo "Submit directory: ${SLURM_SUBMIT_DIR:-unset}"
echo "Launch directory: $PWD"
echo "Job ID: ${SLURM_JOB_ID:-unknown}"
echo "Node list: ${SLURM_NODELIST:-unknown}"
echo "Number of tasks: ${SLURM_NTASKS:-unknown}"
echo "==========================="

if [ -n "${SLURM_SUBMIT_DIR:-}" ]; then
    cd "$SLURM_SUBMIT_DIR"
fi
caseDir="$(pwd)"

export MPLCONFIGDIR="$caseDir/.mplconfig"
mkdir -p "$MPLCONFIGDIR"

OF2506_IMAGE="${OF2506_IMAGE:-$HOME/openfoam-dev_2506.sif}"
OF2506_USER="${OF2506_USER:-x-rkanak1}"
export PYTHON="${PYTHON:-/home/x-rkanak1/.conda/envs/isw_env/bin/python3}"
export FOAM_SIGFPE=0

of2506()
{
    apptainer exec --cleanenv --env USER="$OF2506_USER" "$OF2506_IMAGE" \
        bash -lc "source /openfoam/bash.rc && $*"
}

echo "=== Environment Check ==="
command -v apptainer
test -f "$OF2506_IMAGE"
echo "PYTHON=$PYTHON"
echo "OF2506_IMAGE=$OF2506_IMAGE"
of2506 '
echo "WM_PROJECT_DIR=$WM_PROJECT_DIR"
command -v laserbeamFoam
'

subdomains="$(of2506 "foamDictionary -entry numberOfSubdomains -value '$caseDir/system/decomposeParDict'")"
echo "numberOfSubdomains: $subdomains"

"$PYTHON" - <<'PY'
import importlib, sys
required = ["matplotlib", "numpy", "pandas", "pyvista", "scipy", "vtk"]
missing = []
for module in required:
    try:
        importlib.import_module(module)
    except Exception as exc:
        missing.append(f"{module}: {exc}")
if missing:
    print("ERROR: missing Python modules:")
    for item in missing:
        print(f"  - {item}")
    sys.exit(1)
print("Python environment OK.")
PY
echo "========================="

if [ "${DRY_RUN:-0}" = "1" ]; then
    echo "DRY RUN — solver not started."
    exit 0
fi

echo "Solver started at: $(date)"

apptainer exec --cleanenv --env USER="$OF2506_USER" "$OF2506_IMAGE" bash -lc "
    source /openfoam/bash.rc &&
    export FOAM_SIGFPE=0 &&
    export PYTHON='$PYTHON' &&
    export MPLCONFIGDIR=/tmp &&
    cd '$caseDir' &&
    ./Allrun_long
"

echo "Solver completed at: $(date)"
