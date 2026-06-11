#!/bin/bash
#SBATCH --job-name=kubeAL6061
#SBATCH --time=12:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=8
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
    echo "ERROR: kubeAL6061 job failed"
    echo "Exit status: $status"
    echo "Line: $line"
    echo "Command: $command"
    echo "Working directory: $PWD"
    echo "Time: $(date)"
    echo

    exit "$status"
}

trap on_error ERR

latestInternalVTU()
{
    find VTK -path "*/internal.vtu" -type f -printf '%T@ %p\n' 2>/dev/null \
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

    vtkFile="$(latestInternalVTU)"

    if [ -z "$vtkFile" ]; then
        echo "WARNING: no internal.vtu found after converting ${timeName}" >&2
        return 0
    fi

    MPLCONFIGDIR="${MPLCONFIGDIR:-/tmp}" "$PYTHON" \
        scripts/analyze_meltpool_vtu.py \
        --vtk-file "$vtkFile" \
        || echo "WARNING: meltpool analysis failed for ${vtkFile}" >&2
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

runCaseInsideContainer()
{
    source /openfoam/bash.rc

    . "$WM_PROJECT_DIR/bin/tools/RunFunctions"

    export FOAM_SIGFPE="${FOAM_SIGFPE:-0}"
    export MPLCONFIGDIR="${MPLCONFIGDIR:-/tmp}"

    echo "=== OpenFOAM environment ==="
    echo "WM_PROJECT_VERSION=${WM_PROJECT_VERSION:-unset}"
    echo "FOAM_APPBIN=${FOAM_APPBIN:-unset}"
    echo "FOAM_USER_APPBIN=${FOAM_USER_APPBIN:-unset}"
    command -v blockMesh
    command -v setFields
    command -v decomposePar
    command -v foamDictionary
    command -v laserbeamFoam
    command -v mpirun
    echo "============================"

    echo "Copying 'initial' to 0"
    rm -rf 0 processor* VTK VTKs postProcessing log.*
    cp -r initial 0
    find 0 -name "*.orig" -exec bash -c 'cp "$0" "${0%.orig}"' {} \;

    runApplication blockMesh
    runApplication setFields
    runApplication decomposePar

    local nCores=""
    nCores="$(foamDictionary -entry numberOfSubdomains -value system/decomposeParDict)"

    echo "Running laserbeamFoam on ${nCores} MPI ranks"
    mpirun --oversubscribe -np "$nCores" \
        laserbeamFoam \
        -parallel \
        2>&1 | convertOnWrite | tee log.laserbeamFoam

    touch kubeAL6061.foam
}

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

if [ "${KUBE_JOB_IN_CONTAINER:-0}" = "1" ]; then
    runCaseInsideContainer
    exit 0
fi

echo "=== Kube Al6061 HPC job ==="
echo "Case directory: $caseDir"
echo "Job ID: ${SLURM_JOB_ID:-unknown}"
echo "Node list: ${SLURM_NODELIST:-unknown}"
echo "Tasks: ${SLURM_NTASKS:-unknown}"
echo "OF2506_IMAGE: $OF2506_IMAGE"
echo "PYTHON: $PYTHON"
echo "Started: $(date)"
echo "==========================="

command -v apptainer
test -f "$OF2506_IMAGE"
test -x "$PYTHON"
test -f ./foamVTK.sh
test -f ./scripts/analyze_meltpool_vtu.py

apptainer exec \
    --cleanenv \
    --env USER="$OF2506_USER" \
    --env PYTHON="$PYTHON" \
    --env MPLCONFIGDIR="$MPLCONFIGDIR" \
    --env FOAM_SIGFPE="$FOAM_SIGFPE" \
    --env KUBE_JOB_IN_CONTAINER=1 \
    "$OF2506_IMAGE" \
    bash "$caseDir/job.sh"

echo "Finished: $(date)"
