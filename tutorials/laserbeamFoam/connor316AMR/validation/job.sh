#!/bin/bash
#SBATCH --job-name=validationAMR
#SBATCH --time=24:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=64
#SBATCH --account=mch250110
#SBATCH --output=validation_%j.out
#SBATCH --error=validation_%j.err

# Slurm driver for validation/cases.json.
#
# Submit default Hofmann validation cases:
#   sbatch job.sh
#
# Submit selected cases:
#   sbatch job.sh --case hofmann_scantrack_200W_900mms_r25um --case zhang_N01_260W_520mms_BP
#   sbatch --export=ALL,VALIDATION_CASES="case_a,case_b" job.sh
#
# Sizing:
#   total cores come from SLURM_NTASKS, or VALIDATION_TOTAL_CORES outside Slurm.
#   each case gets at most VALIDATION_MAX_CORES_PER_CASE cores, capped at 8.
#   default is 8 cores/case, so 64 tasks runs 8 cases at once.

set -Eeo pipefail

on_error()
{
    local status="$?"
    echo
    echo "ERROR: validation/job.sh failed (exit $status) at line ${BASH_LINENO[0]:-?}"
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
suiteDir="$(pwd)"

export OF2506_IMAGE="${OF2506_IMAGE:-$HOME/openfoam-dev_2506.sif}"
export OF2506_USER="${OF2506_USER:-x-rkanak1}"
export PYTHON="${PYTHON:-/home/x-rkanak1/.conda/envs/isw_env/bin/python3}"
export FOAM_SIGFPE=0
export MPLCONFIGDIR="$suiteDir/.mplconfig"
mkdir -p "$MPLCONFIGDIR"

totalCores="${VALIDATION_TOTAL_CORES:-${SLURM_NTASKS:-64}}"
maxCoresPerCase="${VALIDATION_MAX_CORES_PER_CASE:-8}"
if [ "$maxCoresPerCase" -gt 8 ]; then
    maxCoresPerCase=8
fi
if [ "$maxCoresPerCase" -lt 1 ]; then
    maxCoresPerCase=1
fi

DEFAULT_HOFMANN_CASES=(
    hofmann_scantrack_200W_900mms_r25um
    hofmann_scantrack_200W_1200mms_r25um
    hofmann_scantrack_200W_1500mms_r25um
    hofmann_scantrack_250W_600mms_r25um
)

setupArgs=("$@")
if [ "${#setupArgs[@]}" -eq 0 ] && [ -n "${VALIDATION_CASES:-}" ]; then
    caseList="${VALIDATION_CASES//,/ }"
    setupArgs=()
    for caseName in $caseList; do
        setupArgs+=(--case "$caseName")
    done
elif [ "${#setupArgs[@]}" -eq 0 ]; then
    setupArgs=()
    for caseName in "${DEFAULT_HOFMANN_CASES[@]}"; do
        setupArgs+=(--case "$caseName")
    done
fi

of2506()
{
    apptainer exec --cleanenv \
        --env USER="$OF2506_USER" \
        --env PYTHON="$PYTHON" \
        --env FOAM_SIGFPE=0 \
        --env MPLCONFIGDIR=/tmp \
        "$OF2506_IMAGE" \
        bash -lc "source /openfoam/bash.rc && $*"
}

echo "=== Slurm validation bootstrap ==="
echo "Suite directory: $suiteDir"
echo "Job ID:          ${SLURM_JOB_ID:-unknown}"
echo "Node list:       ${SLURM_NODELIST:-unknown}"
echo "Slurm tasks:     ${SLURM_NTASKS:-unset}"
echo "Total cores:     $totalCores"
echo "Max cores/case:  $maxCoresPerCase"
echo "OF2506_IMAGE:    $OF2506_IMAGE"
echo "PYTHON:          $PYTHON"
echo "Selection args:  ${setupArgs[*]:-(enabled cases)}"
echo "=================================="

if [ "${DRY_RUN:-0}" = "1" ]; then
    echo "=== DRY RUN: skipping Apptainer/OpenFOAM environment check ==="
else
    echo "=== Environment Check ==="
    command -v apptainer
    test -f "$OF2506_IMAGE"
    "$PYTHON" - <<'PY'
import importlib
import sys

required = ["matplotlib", "numpy", "pandas", "pyvista", "rich", "scipy", "tqdm", "vtk"]
missing = []
for module in required:
    try:
        importlib.import_module(module)
    except Exception as exc:
        missing.append(f"{module}: {exc}")

if missing:
    print("ERROR: missing Python modules needed by validation post-processing:")
    for item in missing:
        print(f"  - {item}")
    sys.exit(1)

print("Host Python environment OK.")
PY
    of2506 '
echo "WM_PROJECT_DIR=$WM_PROJECT_DIR"
echo "WM_PROJECT_VERSION=${WM_PROJECT_VERSION:-unset}"
command -v mpirun
command -v blockMesh
command -v setFields
command -v decomposePar
command -v foamDictionary
command -v reconstructPar
command -v reconstructParMesh
command -v foamToVTK
command -v laserbeamFoam
"$PYTHON" - <<'"'"'PY'"'"'
import importlib
import sys

required = ["matplotlib", "numpy", "pandas", "pyvista", "scipy", "tqdm", "vtk"]
missing = []
for module in required:
    try:
        importlib.import_module(module)
    except Exception as exc:
        missing.append(f"{module}: {exc}")

if missing:
    print("ERROR: host Python is not usable inside the Apptainer run:")
    for item in missing:
        print(f"  - {item}")
    sys.exit(1)

print("Container OpenFOAM and analysis Python OK.")
PY
'
    echo "========================="
fi

echo "Preparing validation cases..."
"$PYTHON" scripts/setup_cases.py "${setupArgs[@]}"
mapfile -t caseNames < <("$PYTHON" scripts/setup_cases.py --list "${setupArgs[@]}")

totalCases="${#caseNames[@]}"
if [ "$totalCases" -eq 0 ]; then
    echo "No validation cases selected."
    exit 0
fi

if [ "$totalCores" -lt 1 ]; then
    echo "ERROR: total core count must be positive, got $totalCores" >&2
    exit 2
fi

coresPerCase="$maxCoresPerCase"
if [ "$totalCores" -lt "$coresPerCase" ]; then
    coresPerCase="$totalCores"
fi
parallelCases=$((totalCores / coresPerCase))
if [ "$parallelCases" -lt 1 ]; then
    parallelCases=1
fi
if [ "$parallelCases" -gt "$totalCases" ]; then
    parallelCases="$totalCases"
fi

echo "=== Run Plan ==="
echo "Selected cases:     $totalCases"
echo "Cores per case:     $coresPerCase"
echo "Parallel case slots: $parallelCases"
printf '  - %s\n' "${caseNames[@]}"
echo "================"

set_case_cores()
{
    local caseName="$1"
    local cores="$2"
    local dict="cases/${caseName}/system/decomposeParDict"

    "$PYTHON" - "$dict" "$cores" <<'PY'
import re
import sys
from pathlib import Path

path = Path(sys.argv[1])
cores = sys.argv[2]
text = path.read_text()
new_text, count = re.subn(
    r"^numberOfSubdomains\s+\d+\s*;",
    f"numberOfSubdomains  {cores};",
    text,
    count=1,
    flags=re.MULTILINE,
)
if count != 1:
    raise SystemExit(f"numberOfSubdomains entry not found in {path}")
path.write_text(new_text)
PY
}

run_case()
{
    local caseName="$1"
    local cores="$2"
    local caseDir="$suiteDir/cases/$caseName"

    set_case_cores "$caseName" "$cores"
    echo "[$(date)] START $caseName (${cores} cores)"
    apptainer exec --cleanenv \
        --env USER="$OF2506_USER" \
        --env PYTHON="$PYTHON" \
        --env FOAM_SIGFPE=0 \
        --env MPLCONFIGDIR=/tmp \
        "$OF2506_IMAGE" \
        bash -lc "source /openfoam/bash.rc && cd '$caseDir' && ./Allclean_long && ./Allrun"
    echo "[$(date)] DONE  $caseName"
}

run_wave()
{
    local startIndex="$1"
    local waveSize="$2"
    local pids=()
    local names=()
    local results=()
    local monitorArgs=()
    local monitorPid=""
    local status=0
    local i=""
    local caseName=""

    echo
    echo "===== Starting wave: ${waveSize} case(s), ${coresPerCase} cores each ====="

    for ((i = 0; i < waveSize; i++)); do
        caseName="${caseNames[$((startIndex + i))]}"
        names+=("$caseName")
        monitorArgs+=(--case "$caseName")
        echo "  -> $caseName (log: cases/$caseName/log.validationJob)"
        run_case "$caseName" "$coresPerCase" > "cases/$caseName/log.validationJob" 2>&1 &
        pids+=("$!")
    done

    if [ "${VALIDATION_MONITOR:-1}" = "1" ]; then
        "$PYTHON" scripts/monitor_case_progress.py --cores "$coresPerCase" --interval 5 "${monitorArgs[@]}" &
        monitorPid="$!"
    fi

    for ((i = 0; i < waveSize; i++)); do
        if wait "${pids[$i]}"; then
            results[$i]=0
        else
            results[$i]=1
            status=1
        fi
    done

    if [ -n "$monitorPid" ]; then
        kill "$monitorPid" 2>/dev/null || true
        wait "$monitorPid" 2>/dev/null || true
    fi

    for ((i = 0; i < waveSize; i++)); do
        if [ "${results[$i]}" -eq 0 ]; then
            echo "  OK: ${names[$i]}"
        else
            echo "  FAILED: ${names[$i]} (see cases/${names[$i]}/log.validationJob)" >&2
        fi
    done

    if [ "$status" -ne 0 ]; then
        exit "$status"
    fi
}

if [ "${DRY_RUN:-0}" = "1" ]; then
    echo "=== DRY RUN: cases prepared, solver not started ==="
    exit 0
fi

echo "Validation started at: $(date)"
index=0
while [ "$index" -lt "$totalCases" ]; do
    remaining=$((totalCases - index))
    waveSize="$parallelCases"
    if [ "$remaining" -lt "$waveSize" ]; then
        waveSize="$remaining"
    fi
    run_wave "$index" "$waveSize"
    index=$((index + waveSize))
done
echo "Validation completed at: $(date)"
