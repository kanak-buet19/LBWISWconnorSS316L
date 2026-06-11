#!/bin/bash
#SBATCH --job-name=kubeCal
#SBATCH --time=12:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=72
#SBATCH --account=mch250110
#SBATCH --output=calibrate_%j.out
#SBATCH --error=calibrate_%j.err

set -Eeo pipefail

scriptDir="$(cd "$(dirname "$0")" && pwd)"
calibrationDir="$scriptDir/meltpool_calibration"
reportDir="$calibrationDir/run_reports/${SLURM_JOB_ID:-manual}"
reportCsv="$reportDir/postprocess_events.csv"

initReport()
{
    mkdir -p "$reportDir"
    ( set -o noclobber; printf 'timestamp,case_id,time,event,detail\n' > "$reportCsv" ) 2>/dev/null || true
}

csvEscape()
{
    printf '%s' "$1" | sed 's/"/""/g'
}

recordPostProcessEvent()
{
    local timeName="$1"
    local event="$2"
    local detail="$3"
    local caseId="${CASE_ID:-$(basename "$PWD")}"
    local lockDir="$reportDir/.postprocess_events.lock"

    initReport

    while ! mkdir "$lockDir" 2>/dev/null; do
        sleep 0.1
    done

    printf '"%s","%s","%s","%s","%s"\n' \
        "$(date -Is)" \
        "$(csvEscape "$caseId")" \
        "$(csvEscape "$timeName")" \
        "$(csvEscape "$event")" \
        "$(csvEscape "$detail")" >> "$reportCsv"

    rm -rf "$lockDir"
}

updateLiveCalibrationPlot()
{
    local lockDir="$calibrationDir/.live_plot.lock"
    local plotLog="${CASE_LOG_DIR:-$reportDir}/log.live_plot"

    mkdir -p "$(dirname "$plotLog")"

    if mkdir "$lockDir" 2>/dev/null; then
        trap 'rm -rf "$lockDir"' RETURN
        "$PYTHON" "$calibrationDir/plot_live_calibration.py" \
            --base-case "$scriptDir" \
            --runs-dir "$calibrationDir/runs" \
            --config "$calibrationDir/calibration_config.json" \
            >> "$plotLog" 2>&1 \
            || {
                echo "WARNING: live calibration plot failed" >&2
                recordPostProcessEvent "-" "live_plot_failed" "plot_live_calibration.py failed"
            }
        rm -rf "$lockDir"
        trap - RETURN
    fi
}

on_error()
{
    local status="$?"
    local line="${BASH_LINENO[0]:-unknown}"
    local command="${BASH_COMMAND:-unknown}"

    echo
    echo "ERROR: kube calibration job failed"
    echo "Exit status: $status"
    echo "Line: $line"
    echo "Command: $command"
    echo "Working directory: $PWD"
    echo "Time: $(date)"
    echo

    exit "$status"
}

trap on_error ERR

latestConvertedVTK()
{
    find VTK \( \
            -path "*/internal.vtu" \
            -o -path "*/internal.vtk" \
            -o -path "VTK/*.vtk" \
        \) -type f -printf '%T@ %p\n' 2>/dev/null \
        | sort -g \
        | tail -1 \
        | cut -d' ' -f2-
}

postProcessWrittenTime()
{
    local timeName="$1"

    if [ "${KUBE_LIVE_VTK:-1}" = "1" ]; then
        analyzeWrittenTime "$timeName"
        updateLiveCalibrationPlot
    fi
}

analyzeWrittenTime()
{
    local timeName="$1"
    local vtkFile=""
    local analysisLog="${CASE_LOG_DIR:-.}/log.analyze_meltpool_vtu"
    local postProcessLog="${CASE_LOG_DIR:-.}/log.livePostProcess"

    if ! "$scriptDir/foamVTK.sh" "$timeName" >> "$postProcessLog" 2>&1; then
        echo "WARNING: foamVTK.sh failed for time ${timeName}" >&2
        recordPostProcessEvent "$timeName" "foamVTK_failed" "foamVTK.sh failed"
        return 0
    fi

    vtkFile="$(latestConvertedVTK)"

    if [ -z "$vtkFile" ]; then
        echo "WARNING: no internal VTK file found after converting ${timeName}" >&2
        recordPostProcessEvent "$timeName" "no_vtk" "No converted VTK file found after conversion"
        return 0
    fi

    if ! MPLCONFIGDIR="${MPLCONFIGDIR:-/tmp}" "$PYTHON" \
        scripts/analyze_meltpool_vtu.py \
        --vtk-file "$vtkFile" \
        >> "$analysisLog" 2>&1
    then
        echo "WARNING: meltpool analysis failed for ${vtkFile}" >&2
        recordPostProcessEvent "$timeName" "analysis_failed" "$vtkFile"
    fi
}

postProcessFinishedCase()
{
    local timeName=""

    timeName="$(
        find processor0 -maxdepth 1 -type d -printf '%f\n' 2>/dev/null \
            | awk '$1 ~ /^[0-9]/ {print}' \
            | sort -g \
            | tail -1
    )"

    if [ -z "$timeName" ]; then
        echo "WARNING: no written time found for final post-processing" >&2
        return 0
    fi

    echo "Final VTK conversion for time ${timeName}"
    analyzeWrittenTime "$timeName"
    updateLiveCalibrationPlot
}

postProcessWorkerLoop()
{
    local lockDir=".postProcessLive.lock"
    local latestFile=".postProcessLatest"
    local queuedTime=""

    trap 'rm -rf "$lockDir"' EXIT

    while [ -s "$latestFile" ]; do
        queuedTime="$(tail -1 "$latestFile" 2>/dev/null || true)"
        : > "$latestFile"

        if [ -n "$queuedTime" ]; then
            postProcessWrittenTime "$queuedTime"
        fi
    done

    rm -rf "$lockDir"
    trap - EXIT
}

queuePostProcessTime()
{
    local timeName="$1"
    local lockDir=".postProcessLive.lock"
    local latestFile=".postProcessLatest"

    if [ "${KUBE_LIVE_VTK:-1}" != "1" ]; then
        return 0
    fi

    printf '%s\n' "$timeName" > "$latestFile"

    if mkdir "$lockDir" 2>/dev/null; then
        postProcessWorkerLoop >> "${CASE_LOG_DIR:-.}/log.livePostProcess" 2>&1 < /dev/null &
    fi
}

waitForLivePostProcessing()
{
    local lockDir=".postProcessLive.lock"
    local attempt=0

    while [ -d "$lockDir" ] && [ "$attempt" -lt 600 ]; do
        sleep 1
        attempt=$((attempt + 1))
    done

    if [ -d "$lockDir" ]; then
        echo "WARNING: live post-processing still running; continuing without waiting longer" >&2
    fi
}

convertOnWrite()
{
    local currentTime=""
    local lastProcessedTime=""

    while IFS= read -r line
    do
        case "$line" in
            "Time = "*)
                currentTime="${line#Time = }"
                ;;
            "Writing "*|"Writing data at time "*)
                if [ -n "$currentTime" ] && [ "$currentTime" != "$lastProcessedTime" ]; then
                    echo "[$(date -Is)] $(basename "$PWD"): queued live VTK for time ${currentTime}"
                    queuePostProcessTime "$currentTime"
                    lastProcessedTime="$currentTime"
                fi
                ;;
        esac
    done
}

runSingleCase()
{
    source /openfoam/bash.rc
    export FOAM_TUTORIALS="${FOAM_TUTORIALS:-$PWD}"
    . "$WM_PROJECT_DIR/bin/tools/RunFunctions"

    export FOAM_SIGFPE="${FOAM_SIGFPE:-0}"
    export MPLCONFIGDIR="${MPLCONFIGDIR:-/tmp}"
    export CASE_ID="$(basename "$PWD")"
    export CASE_LOG_DIR="$calibrationDir/logs/$CASE_ID"
    mkdir -p "$CASE_LOG_DIR"
    rm -f "$CASE_LOG_DIR"/log.*
    initReport

    echo "=== Running calibration case: $PWD ==="
    echo "Case logs: $CASE_LOG_DIR"
    command -v blockMesh
    command -v setFields
    command -v decomposePar
    command -v foamDictionary
    command -v laserbeamFoam
    command -v mpirun

    rm -rf 0 processor* VTK VTKs postProcessing log.* STOP_CASE STOP_REASON.json .solverPid .postProcessLive.lock .postProcessLatest
    cp -r initial 0
    find 0 -name "*.orig" -exec bash -c 'cp "$0" "${0%.orig}"' {} \;

    runApplication blockMesh
    runApplication setFields
    runApplication decomposePar
    mv -f log.blockMesh log.setFields log.decomposePar "$CASE_LOG_DIR"/ 2>/dev/null || true

    local nCores=""
    nCores="$(foamDictionary -entry numberOfSubdomains -value system/decomposeParDict)"

    echo "Running laserbeamFoam on ${nCores} MPI ranks"
    set +e
    (
        set -m
        mpirun --oversubscribe -np "$nCores" \
            laserbeamFoam \
            -parallel \
            2>&1 &
        solverPid="$!"
        printf '%s\n' "$solverPid" > .solverPid
        wait "$solverPid"
    ) | tee "$CASE_LOG_DIR/log.laserbeamFoam" | convertOnWrite
    solverStatus="${PIPESTATUS[0]}"
    set -e

    if [ "$solverStatus" -ne 0 ]; then
        exit "$solverStatus"
    fi

    touch kubeAL6061.foam
    waitForLivePostProcessing
    postProcessFinishedCase
}

runCalibrationInsideContainer()
{
    source /openfoam/bash.rc
    export FOAM_TUTORIALS="${FOAM_TUTORIALS:-$PWD}"

    echo "=== OpenFOAM environment ==="
    echo "WM_PROJECT_VERSION=${WM_PROJECT_VERSION:-unset}"
    echo "FOAM_USER_APPBIN=${FOAM_USER_APPBIN:-unset}"
    command -v laserbeamFoam
    command -v mpirun
    echo "============================"

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
    print("Missing Python modules:")
    for item in missing:
        print("  -", item)
    sys.exit(1)

print("Python environment OK")
PY

    local configPath="${CALIBRATION_CONFIG:-meltpool_calibration/calibration_config.json}"
    local jobs="${MAX_PARALLEL_CASES:-8}"
    local runCommand=""
    local forceArgs=()

    case "$configPath" in
        /*) ;;
        *) configPath="$PWD/$configPath" ;;
    esac
    export CALIBRATION_CONFIG="$configPath"

    runCommand="KUBE_CALIBRATE_RUN_CASE=1 CALIBRATION_CONFIG='$configPath' PYTHON='$PYTHON' MPLCONFIGDIR='${MPLCONFIGDIR:-/tmp}' FOAM_SIGFPE='${FOAM_SIGFPE:-0}' bash '$PWD/calibrate_job.sh'"
    if [ "${FORCE_REBUILD_CASES:-0}" = "1" ]; then
        forceArgs=(--force)
    fi

    "$PYTHON" meltpool_calibration/kube_sweep.py \
        --config "$configPath" \
        --jobs "$jobs" \
        "${forceArgs[@]}" \
        --run-command "$runCommand"
}

if [ "${KUBE_CALIBRATE_RUN_CASE:-0}" = "1" ]; then
    runSingleCase
    exit 0
fi

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

if [ "${KUBE_CALIBRATE_IN_CONTAINER:-0}" = "1" ]; then
    runCalibrationInsideContainer
    exit 0
fi

echo "=== Kube Al6061 calibration job ==="
echo "Case directory: $caseDir"
echo "Job ID: ${SLURM_JOB_ID:-unknown}"
echo "Node list: ${SLURM_NODELIST:-unknown}"
echo "Tasks: ${SLURM_NTASKS:-unknown}"
echo "OF2506_IMAGE: $OF2506_IMAGE"
echo "PYTHON: $PYTHON"
echo "Config: ${CALIBRATION_CONFIG:-meltpool_calibration/calibration_config.json}"
echo "Started: $(date)"
echo "==================================="

command -v apptainer
test -f "$OF2506_IMAGE"
test -x "$PYTHON"
test -f ./meltpool_calibration/kube_sweep.py
test -f "${CALIBRATION_CONFIG:-meltpool_calibration/calibration_config.json}"

apptainer exec \
    --cleanenv \
    --env USER="$OF2506_USER" \
    --env PYTHON="$PYTHON" \
    --env MPLCONFIGDIR="$MPLCONFIGDIR" \
    --env FOAM_SIGFPE="$FOAM_SIGFPE" \
    --env CALIBRATION_CONFIG="${CALIBRATION_CONFIG:-meltpool_calibration/calibration_config.json}" \
    --env MAX_PARALLEL_CASES="${MAX_PARALLEL_CASES:-8}" \
    --env FORCE_REBUILD_CASES="${FORCE_REBUILD_CASES:-0}" \
    --env KUBE_LIVE_VTK="${KUBE_LIVE_VTK:-1}" \
    --env KUBE_CALIBRATE_IN_CONTAINER=1 \
    "$OF2506_IMAGE" \
    bash "$caseDir/calibrate_job.sh"

echo "Finished: $(date)"
