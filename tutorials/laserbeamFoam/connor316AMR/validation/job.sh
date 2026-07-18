#!/bin/bash
#SBATCH --job-name=validationAMR
#SBATCH --time=48:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=36
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
#   every selected case runs concurrently (one wave), each using its own
#   per-case `cores` from cases.json (default 8; see setup_cases.py set_cores).
#   Size your #SBATCH allocation (or VALIDATION_TOTAL_CORES outside Slurm) to
#   match the SUM of the selected cases' cores yourself -- this script no
#   longer overrides per-case core counts or auto-computes parallelism.

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

# Informational only now -- used to warn if it doesn't match the sum of
# per-case cores below, not to size or cap anything.
totalCores="${VALIDATION_TOTAL_CORES:-${SLURM_NTASKS:-64}}"

DEFAULT_HOFMANN_CASES=(
    scantrack_SS316_200W_200mms_r40um
    scantrack_SS316_200W_200mms_r40um_linearSigma
    scantrack_SS316_200W_200mms_r40um_constantSigmaAtTvap
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
command -v laserbeamFoamISW
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

mapfile -t caseNames < <("$PYTHON" scripts/setup_cases.py --list "${setupArgs[@]}")

totalCases="${#caseNames[@]}"
if [ "$totalCases" -eq 0 ]; then
    echo "No validation cases selected."
    exit 0
fi

# One batch folder per job submission: runs/case_NNN_<date>_<time>/, holding
# every selected case as runs/<batchName>/<caseName>/ (built directly via
# --dest-root, no separate staging copy). NNN is a global counter shared with
# validation/Allrun, so local runs and sbatch submissions never collide.
runTimestamp="$(date +%Y%m%d_%H%M%S)"

nextBatchName()
{
    local maxN=0
    local n
    for d in runs/case_[0-9][0-9][0-9]_*; do
        [ -d "$d" ] || continue
        n="$(basename "$d" | sed -n 's/^case_\([0-9]\{3\}\)_.*/\1/p')"
        [ -n "$n" ] || continue
        n=$((10#$n))
        if [ "$n" -gt "$maxN" ]; then
            maxN="$n"
        fi
    done
    printf 'case_%03d_%s' "$((maxN + 1))" "$runTimestamp"
}

batchName="$(nextBatchName)"
echo "Preparing validation cases into runs/${batchName}/ ..."
"$PYTHON" scripts/setup_cases.py --dest-root "runs/${batchName}" "${setupArgs[@]}"

# Each case's numberOfSubdomains was already set by setup_cases.py (from its
# cases.json `cores`, default 8) -- read it back as the source of truth
# instead of re-deriving or overwriting it here.
caseCores=()
sumCores=0
for caseName in "${caseNames[@]}"; do
    dict="runs/${batchName}/${caseName}/system/decomposeParDict"
    cores="$(sed -n 's/^numberOfSubdomains[[:space:]]\+\([0-9]\+\).*/\1/p' "$dict" | head -1)"
    if [ -z "$cores" ]; then
        echo "ERROR: could not read numberOfSubdomains from $dict" >&2
        exit 2
    fi
    caseCores+=("$cores")
    sumCores=$((sumCores + cores))
done

echo "=== Run Plan ==="
echo "Batch folder:       runs/${batchName}"
echo "Selected cases:     $totalCases"
echo "Allocated cores:    $totalCores"
echo "Sum of case cores:  $sumCores"
if [ "$sumCores" -ne "$totalCores" ]; then
    echo "WARNING: allocated cores ($totalCores) != sum of per-case cores ($sumCores)." >&2
    echo "         Size your #SBATCH allocation / VALIDATION_TOTAL_CORES to match the sum." >&2
fi
for i in "${!caseNames[@]}"; do
    echo "  - ${caseNames[$i]} (${caseCores[$i]} cores)"
done
echo "================"

run_case()
{
    local caseName="$1"
    local cores="$2"
    local caseDir="$suiteDir/runs/${batchName}/${caseName}"

    echo "[$(date)] START $caseName (${cores} cores)"
    apptainer exec --cleanenv \
        --env USER="$OF2506_USER" \
        --env PYTHON="$PYTHON" \
        --env FOAM_SIGFPE=0 \
        --env MPLCONFIGDIR=/tmp \
        "$OF2506_IMAGE" \
        bash -lc "source /openfoam/bash.rc && cd '$caseDir' && ./Allrun"
    echo "[$(date)] DONE  $caseName"
}

if [ "${DRY_RUN:-0}" = "1" ]; then
    echo "=== DRY RUN: cases prepared, solver not started ==="
    exit 0
fi

echo "Validation started at: $(date)"
echo
echo "===== Starting all ${totalCases} case(s) concurrently ====="

pids=()
monitorArgs=()
monitorPid=""
status=0

for i in "${!caseNames[@]}"; do
    caseName="${caseNames[$i]}"
    cores="${caseCores[$i]}"
    monitorArgs+=(--case "${batchName}/${caseName}")
    echo "  -> $caseName (log: runs/${batchName}/$caseName/log.validationJob)"
    run_case "$caseName" "$cores" > "runs/${batchName}/$caseName/log.validationJob" 2>&1 &
    pids+=("$!")
done

if [ "${VALIDATION_MONITOR:-1}" = "1" ]; then
    "$PYTHON" scripts/monitor_case_progress.py --cores "$sumCores" --interval 5 "${monitorArgs[@]}" &
    monitorPid="$!"
fi

for i in "${!caseNames[@]}"; do
    if ! wait "${pids[$i]}"; then
        status=1
        echo "  FAILED: ${caseNames[$i]} (see runs/${batchName}/${caseNames[$i]}/log.validationJob)" >&2
    fi
done

if [ -n "$monitorPid" ]; then
    kill "$monitorPid" 2>/dev/null || true
    wait "$monitorPid" 2>/dev/null || true
fi

echo "Validation completed at: $(date)"
if [ "$status" -ne 0 ]; then
    exit "$status"
fi
