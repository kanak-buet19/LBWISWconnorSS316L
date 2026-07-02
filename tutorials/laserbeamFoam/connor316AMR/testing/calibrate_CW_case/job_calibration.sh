#!/bin/bash
#SBATCH --job-name=cwCalib
#SBATCH --time=72:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=32
#SBATCH --account=mch250110
#SBATCH --output=calib_%j.out
#SBATCH --error=calib_%j.err
#
# HPC driver for the CW melt-pool calibration swarm (Slurm + Apptainer).
#
# Sizing the swarm:
#   --ntasks-per-node  = total cores the orchestrator may use (= CALIB_TOTAL_CORES)
#   CALIB_CORES_PER_SIM = cores per individual laserbeamFoam sim
#   parallel sims      = CALIB_TOTAL_CORES / CALIB_CORES_PER_SIM
#   e.g. 32 cores / 4 per sim = 8 sims at once:
#        sbatch --ntasks-per-node=32 --export=ALL,CALIB_CORES_PER_SIM=4 job_calibration.sh
#
# Submit:  sbatch job_calibration.sh
#
set -Eeo pipefail

on_error()
{
    local status="$?"
    echo
    echo "ERROR: job_calibration.sh failed (exit $status) at line ${BASH_LINENO[0]:-?}"
    echo "Command: ${BASH_COMMAND:-?}"
    echo "Time: $(date)"
    exit "$status"
}
trap on_error ERR

echo "=== Slurm CW calibration bootstrap ==="
echo "Job ID:    ${SLURM_JOB_ID:-unknown}"
echo "Node list: ${SLURM_NODELIST:-unknown}"
echo "NTasks:    ${SLURM_NTASKS:-unknown}"
echo "======================================"

if [ -n "${SLURM_SUBMIT_DIR:-}" ]; then
    cd "$SLURM_SUBMIT_DIR"
fi
caseDir="$(pwd)"

# --- environment (mirror connor316AMR job.sh) ---
export CALIB_USE_APPTAINER=1   # run each sim inside the OpenFOAM image (HPC)
export OF2506_IMAGE="${OF2506_IMAGE:-$HOME/openfoam-dev_2506.sif}"
export OF2506_USER="${OF2506_USER:-x-rkanak1}"
export PYTHON="${PYTHON:-/home/x-rkanak1/.conda/envs/isw_env/bin/python3}"
export FOAM_SIGFPE=0
export MPLCONFIGDIR="$caseDir/.mplconfig"
mkdir -p "$MPLCONFIGDIR"

# --- swarm sizing (dynamic) ---
export CALIB_TOTAL_CORES="${CALIB_TOTAL_CORES:-${SLURM_NTASKS:-32}}"
export CALIB_CORES_PER_SIM="${CALIB_CORES_PER_SIM:-4}"
echo "CALIB_TOTAL_CORES   = $CALIB_TOTAL_CORES"
echo "CALIB_CORES_PER_SIM = $CALIB_CORES_PER_SIM"
echo "parallel sims       = $(( CALIB_TOTAL_CORES / CALIB_CORES_PER_SIM ))"

# --- sanity checks ---
echo "=== Environment Check ==="
command -v apptainer
test -f "$OF2506_IMAGE"
"$PYTHON" - <<'PY'
import importlib, sys
req = ["matplotlib", "numpy", "optuna", "pandas", "pyvista", "scipy", "tqdm", "vtk"]
miss = []
for m in req:
    try:
        importlib.import_module(m)
    except Exception as e:
        miss.append(f"{m}: {e}")
if miss:
    print("ERROR: missing Python modules:")
    [print("  -", x) for x in miss]
    sys.exit(1)
print("Python environment OK.")
PY
apptainer exec --cleanenv --env USER="$OF2506_USER" "$OF2506_IMAGE" \
    bash -lc 'source /openfoam/bash.rc && command -v laserbeamFoam blockMesh foamToVTK >/dev/null && echo "OpenFOAM image OK."'
echo "========================="

if [ "${DRY_RUN:-0}" = "1" ]; then
    echo "=== DRY RUN: not launching swarm ==="
    exit 0
fi

echo "Calibration started at: $(date)"
"$PYTHON" scripts/run_calibration.py "$@"
echo "Calibration completed at: $(date)"
