#!/bin/bash

set -Eeo pipefail

scriptDir="$(cd "$(dirname "$0")" && pwd)"
baseCase="$(cd "$scriptDir/.." && pwd)"

stopSolver()
{
    local reason="$1"
    local pid=""

    echo "Stopping case: ${reason}"
    touch STOP_CASE
    if [ -f .solverPid ]; then
        pid="$(cat .solverPid)"
    fi
    if [ -n "$pid" ]; then
        kill -TERM "-$pid" 2>/dev/null || kill -TERM "$pid" 2>/dev/null || true
    fi
}

stopSolverIfRequested()
{
    local decision=""

    decision="$(
        "$PYTHON" "$scriptDir/onset_stop_check.py" \
            --case "$PWD" \
            --config "$CALIBRATION_CONFIG"
    )"

    case "$decision" in
        stop_early|stop_late|stop_matched)
            echo "Early-stop decision: ${decision}"
            stopSolver "$decision"
            ;;
    esac
}

updateLiveCalibrationPlot()
{
    local lockDir="$scriptDir/.live_plot.lock"

    if mkdir "$lockDir" 2>/dev/null; then
        trap 'rm -rf "$lockDir"' RETURN
        "$PYTHON" "$scriptDir/plot_live_calibration.py" \
            --base-case "$baseCase" \
            --runs-dir "$scriptDir/runs" \
            --config "$scriptDir/calibration_config.json" \
            || echo "WARNING: live calibration plot failed" >&2
        rm -rf "$lockDir"
        trap - RETURN
    fi
}

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
    local vtkFile=""

    if ! ./foamVTK.sh "$timeName"; then
        echo "WARNING: foamVTK.sh failed for time ${timeName}" >&2
        return 0
    fi

    vtkFile="$(latestConvertedVTK)"

    if [ -z "$vtkFile" ]; then
        echo "WARNING: no converted VTK file found after converting ${timeName}" >&2
        return 0
    fi

    MPLCONFIGDIR="${MPLCONFIGDIR:-/tmp}" "$PYTHON" \
        scripts/analyze_meltpool_vtu.py \
        --vtk-file "$vtkFile" \
        || echo "WARNING: meltpool analysis failed for ${vtkFile}" >&2

    stopSolverIfRequested
    updateLiveCalibrationPlot
}

postProcessQueueWorker()
{
    local queueFile="$1"
    local queuedTime=""

    while IFS= read -r queuedTime
    do
        if [ "$queuedTime" = "__DONE__" ]; then
            break
        fi

        postProcessWrittenTime "$queuedTime"
    done < "$queueFile"
}

convertOnWrite()
{
    local currentTime=""
    local lastProcessedTime=""
    local queueFile=""
    local workerPid=""

    queueFile="$(mktemp -u ".postProcessQueue.XXXXXX")"
    mkfifo "$queueFile"
    postProcessQueueWorker "$queueFile" &
    workerPid="$!"
    exec 3>"$queueFile"
    rm -f "$queueFile"

    while IFS= read -r line
    do
        echo "$line"

        case "$line" in
            "Time = "*)
                currentTime="${line#Time = }"
                ;;
            "Writing "*|"Writing data at time "*)
                if [ -n "$currentTime" ] && [ "$currentTime" != "$lastProcessedTime" ]; then
                    echo "Queueing VTK conversion for time ${currentTime}"
                    printf '%s\n' "$currentTime" >&3
                    lastProcessedTime="$currentTime"
                fi
                ;;
        esac
    done

    printf '%s\n' "__DONE__" >&3
    exec 3>&-
    wait "$workerPid"
}

if [ -z "${WM_PROJECT_DIR:-}" ]; then
    echo "ERROR: OpenFOAM environment is not sourced. Source OpenFOAM first." >&2
    exit 1
fi

. "$WM_PROJECT_DIR/bin/tools/RunFunctions"

export PYTHON="${PYTHON:-python3}"
export MPLCONFIGDIR="${MPLCONFIGDIR:-/tmp}"
export FOAM_SIGFPE="${FOAM_SIGFPE:-0}"
export CALIBRATION_CONFIG="${CALIBRATION_CONFIG:-$scriptDir/calibration_config.json}"

echo "=== Running local calibration case: $PWD ==="
rm -rf 0 processor* VTK VTKs postProcessing log.* STOP_CASE STOP_REASON.json .solverPid
cp -r initial 0
find 0 -name "*.orig" -exec bash -c 'cp "$0" "${0%.orig}"' {} \;

runApplication blockMesh
runApplication setFields
runApplication decomposePar

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
) | convertOnWrite | tee log.laserbeamFoam
solverStatus="${PIPESTATUS[0]}"
set -e

if [ "$solverStatus" -ne 0 ] && [ ! -f STOP_CASE ]; then
    exit "$solverStatus"
fi

touch kubeAL6061.foam
