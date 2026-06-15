#!/bin/bash
#SBATCH --job-name=abdali_benchmark
#SBATCH --time=24:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=32
#SBATCH --account=mch250110
#SBATCH --output=output_%j.out
#SBATCH --error=output_%j.err
#SBATCH --mail-user=rakibul.buet19@gmail.com
#SBATCH --mail-type=END,FAIL

# Autogenous SLM 316L single-track parallel scaling benchmark (Abdali parameters).
# Sweeps across 4, 8, 16, 24, and 32 cores for a 0.09ms simulation test.
#
# Submit with:   sbatch job.sh

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

echo "=== Slurm connor316AMR benchmark bootstrap ==="
echo "Script: $0"
echo "Submit directory: ${SLURM_SUBMIT_DIR:-unset}"
echo "Launch directory: $PWD"
echo "Job ID: ${SLURM_JOB_ID:-unknown}"
echo "Node list: ${SLURM_NODELIST:-unknown}"
echo "Number of tasks: ${SLURM_NTASKS:-unknown}"
echo "============================================="

# Run from the directory where this script lives.
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
echo "WM_PROJECT_VERSION=${WM_PROJECT_VERSION:-unset}"
command -v mpirun
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
    print("ERROR: missing Python modules needed by post-processing:")
    for item in missing:
        print(f"  - {item}")
    sys.exit(1)

print("Python environment OK: all post-processing modules imported.")
PY
echo "========================="

echo "Benchmark started at: $(date)"

apptainer exec --cleanenv --env USER="$OF2506_USER" "$OF2506_IMAGE" bash -lc "
    source /openfoam/bash.rc &&
    export FOAM_SIGFPE=0 &&
    export PYTHON='$PYTHON' &&
    export MPLCONFIGDIR=/tmp &&
    cd '$caseDir' &&
    ./Allrun
"

echo "Benchmark completed at: $(date)"
