#!/bin/bash
#SBATCH --job-name=taosunMeshConv
#SBATCH --time=4-00:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=88
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
    echo "ERROR: taosun mesh convergence job failed"
    echo "Exit status: $status"
    echo "Line: $line"
    echo "Command: $command"
    echo "Working directory: $PWD"
    echo "Time: $(date)"
    echo

    exit "$status"
}

trap on_error ERR

runPipelineInsideContainer()
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
    "$PYTHON" --version
    echo "============================"

    runMeshConvergence
}

caseNames()
{
    "$PYTHON" - mesh_convergence.json <<'PY'
import json
cfg = json.load(open("mesh_convergence.json"))
for case in cfg["cases"]:
    print(case["name"])
PY
}

caseValue()
{
    local caseName="$1"
    local key="$2"

    "$PYTHON" - mesh_convergence.json "$caseName" "$key" <<'PY'
import json
import sys
cfg = json.load(open(sys.argv[1]))
case = next(item for item in cfg["cases"] if item["name"] == sys.argv[2])
value = case[sys.argv[3]]
if isinstance(value, list):
    print(" ".join(str(v) for v in value))
else:
    print(value)
PY
}

maxRefinement()
{
    local caseName="$1"

    "$PYTHON" - mesh_convergence.json "$caseName" <<'PY'
import json
import sys
cfg = json.load(open(sys.argv[1]))
case = next(item for item in cfg["cases"] if item["name"] == sys.argv[2])
print(case.get("mesh_refinement", cfg.get("max_refinement", 1)))
PY
}

patchCaseMesh()
{
    local runDir="$1"
    local cells="$2"
    local maxRef="$3"
    local cores="$4"

    "$PYTHON" - "$runDir" "$cells" "$maxRef" "$cores" <<'PY'
from pathlib import Path
import re
import sys

run = Path(sys.argv[1])
cells = tuple(sys.argv[2].split())
max_ref = sys.argv[3]
cores = sys.argv[4]

block_mesh = run / "system" / "blockMeshDict"
text = block_mesh.read_text()
text = re.sub(
    r"hex \(0 1 2 3 4 5 6 7\) \([^)]+\)",
    f"hex (0 1 2 3 4 5 6 7) ({cells[0]} {cells[1]} {cells[2]})",
    text,
)
block_mesh.write_text(text)

dynamic_mesh = run / "constant" / "dynamicMeshDict"
text = dynamic_mesh.read_text()
text = re.sub(r"\bmaxRefinement\s+\d+\s*;", f"maxRefinement   {max_ref};", text)
dynamic_mesh.write_text(text)

decompose = run / "system" / "decomposeParDict"
text = decompose.read_text()
text = re.sub(r"\bnumberOfSubdomains\s+\d+\s*;", f"numberOfSubdomains  {cores};", text)
decompose.write_text(text)
PY
}

prepareRun()
{
    local name="$1"
    local runDir="runs/$name"

    mkdir -p runs

    if [ ! -d "$runDir" ]; then
        echo "Creating $runDir"
        mkdir -p "$runDir"
        cp -a template/. "$runDir/"
    else
        echo "Using existing $runDir"
    fi

    patchCaseMesh \
        "$runDir" \
        "$(caseValue "$name" cells)" \
        "$(maxRefinement "$name")" \
        "$(caseValue "$name" cores)"
}

latestConvertedVTK()
{
    find VTK -maxdepth 1 -type f -name "*.vtk" -printf '%T@ %p\n' 2>/dev/null \
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

    if [ -n "${CURRENT_CASE_NAME:-}" ] && [ -n "${CURRENT_BASE_CELL_UM:-}" ]; then
        collectStats "$CURRENT_CASE_NAME" "$CURRENT_BASE_CELL_UM"
    fi
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

liveStatsLoop()
{
    local name="$1"
    local baseCellUm="$2"
    local interval="${STATS_REFRESH_SECONDS:-60}"

    while true
    do
        collectAbsorptivityStats "$name" "$baseCellUm" || true
        sleep "$interval"
    done
}

runOneCase()
{
    local name="$1"
    local baseCellUm="$2"
    local runDir="runs/$name"
    local nCores=""
    local logMode="-a"
    local statsPid=""

    echo "=== Running $name ==="
    CURRENT_CASE_NAME="$name"
    CURRENT_BASE_CELL_UM="$baseCellUm"
    cd "$runDir"

    if [ "${FRESH:-0}" = "1" ] || [ ! -d processor0 ]; then
        echo "Starting fresh $name"
        rm -rf 0 processor* VTK VTKs postProcessing log.*
        cp -r initial 0
        find 0 -name "*.orig" -exec bash -c 'cp "$0" "${0%.orig}"' {} \;

        runApplication blockMesh
        runApplication setFields
        runApplication decomposePar
        logMode=""
    else
        echo "Resuming $name from latest processor time"
        [ -f log.laserbeamFoam ] || logMode=""
    fi

    nCores="$(foamDictionary -entry numberOfSubdomains -value system/decomposeParDict)"
    echo "Running laserbeamFoam on ${nCores} MPI ranks"

    liveStatsLoop "$name" "$baseCellUm" &
    statsPid="$!"

    mpirun --oversubscribe -np "$nCores" \
        laserbeamFoam \
        -parallel \
        2>&1 | convertOnWrite | tee $logMode log.laserbeamFoam

    kill "$statsPid" 2>/dev/null || true
    wait "$statsPid" 2>/dev/null || true
    collectStats "$name" "$baseCellUm"

    touch taosunAL6061.foam
    cd "$caseDir"
}

collectStats()
{
    local name="$1"
    local baseCellUm="$2"
    local label="${baseCellUm}um"
    local runDir="runs/$name"
    local geometryCsv="$runDir/post-processing-data/vtu_meltpool_geometry.csv"

    (
    cd "$caseDir"

    collectAbsorptivityStats "$name" "$baseCellUm"

    mkdir -p stats/keyhole_depth stats/meltpool_depth

    if [ -f "$geometryCsv" ]; then
        "$PYTHON" - "$geometryCsv" "$label" <<'PY'
import csv
import sys
from pathlib import Path

geometry = Path(sys.argv[1])
label = sys.argv[2]
rows = list(csv.DictReader(geometry.open()))

for folder, field in [
    ("keyhole_depth", "keyholeDepth_um"),
    ("meltpool_depth", "meltPoolDepth_um"),
]:
    out = Path("stats") / folder / f"{folder}_{label}.csv"
    out.parent.mkdir(parents=True, exist_ok=True)
    with out.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=["time_s", field])
        writer.writeheader()
        for row in rows:
            writer.writerow({"time_s": row.get("time", ""), field: row.get(field, "")})
PY
    else
        echo "WARNING: missing $geometryCsv" >&2
    fi
    )
}

collectAbsorptivityStats()
{
    local name="$1"
    local baseCellUm="$2"
    local label="${baseCellUm}um"
    local runDir="runs/$name"
    local absorptivityCsv="$runDir/absorptivity_vs_time/absorptivity_vs_time.csv"

    mkdir -p "$caseDir/stats/absorptivity"

    if [ -f "$caseDir/$absorptivityCsv" ]; then
        cp "$caseDir/$absorptivityCsv" "$caseDir/stats/absorptivity/absorptivity_${label}.csv"
    fi
}

runMeshConvergence()
{
    local name=""
    local baseCellUm=""
    local runSim=""
    local pids=()
    local pidNames=()
    local failures=0
    local i=""

    for name in $(caseNames)
    do
        baseCellUm="$(caseValue "$name" base_cell_um)"
        runSim="$(caseValue "$name" run_sim)"

        prepareRun "$name"

        if [ "${FRESH:-0}" = "1" ]; then
            rm -f "runs/$name/.mesh_convergence_complete"
        fi

        if [ "$runSim" != "True" ] && [ "$runSim" != "true" ] && [ "$runSim" != "1" ]; then
            echo "Skipping $name because run_sim=false"
            collectStats "$name" "$baseCellUm"
            continue
        fi

        if [ -f "runs/$name/.mesh_convergence_complete" ]; then
            echo "Skipping $name because completion marker exists"
            collectStats "$name" "$baseCellUm"
            continue
        fi

        (
            runOneCase "$name" "$baseCellUm"
            collectStats "$name" "$baseCellUm"
            touch "runs/$name/.mesh_convergence_complete"
        ) &
        pids+=("$!")
        pidNames+=("$name")
        echo "Launched $name in background with PID ${pids[-1]}"
    done

    for i in "${!pids[@]}"
    do
        if wait "${pids[$i]}"; then
            echo "Completed ${pidNames[$i]}"
        else
            echo "ERROR: ${pidNames[$i]} failed" >&2
            failures=1
        fi
    done

    if [ "$failures" -ne 0 ]; then
        return 1
    fi

    echo "Mesh convergence pipeline complete."
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

if [ "${TAOSUN_MESH_CONV_IN_CONTAINER:-0}" = "1" ]; then
    runPipelineInsideContainer
    exit 0
fi

echo "=== Tao Sun mesh convergence HPC job ==="
echo "Case directory: $caseDir"
echo "Job ID: ${SLURM_JOB_ID:-unknown}"
echo "Node list: ${SLURM_NODELIST:-unknown}"
echo "Tasks: ${SLURM_NTASKS:-unknown}"
echo "OF2506_IMAGE: $OF2506_IMAGE"
echo "PYTHON: $PYTHON"
echo "Started: $(date)"
echo "Resumable: completed cases have runs/<case>/.mesh_convergence_complete"
echo "==========================="

command -v apptainer
test -f "$OF2506_IMAGE"
test -x "$PYTHON"
test -f ./mesh_convergence.json
test -f ./job.sh
test -d ./template

apptainer exec \
    --cleanenv \
    --env USER="$OF2506_USER" \
    --env PYTHON="$PYTHON" \
    --env MPLCONFIGDIR="$MPLCONFIGDIR" \
    --env FOAM_SIGFPE="$FOAM_SIGFPE" \
    --env TAOSUN_MESH_CONV_IN_CONTAINER=1 \
    "$OF2506_IMAGE" \
    bash "$caseDir/job.sh"

echo "Finished: $(date)"
