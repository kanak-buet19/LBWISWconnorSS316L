#!/bin/bash

set -e

. "$WM_PROJECT_DIR/bin/tools/RunFunctions"

foam_to_vtk()
{
    foamToVTK -legacy "$@"
    if find VTK -type f \( -name "*.vtu" -o -name "*.vtp" -o -name "*.vtm" \) | grep -q .; then
        echo "[ERROR] foamToVTK produced XML VTK files; expected legacy .vtk only." >&2
        return 1
    fi
}

is_numeric_time()
{
    [[ "$1" =~ ^([0-9]+|[0-9]*\.[0-9]+)([eE][-+]?[0-9]+)?$ ]]
}

remove_reconstructed_time()
{
    local timeName="$1"
    if [ "$timeName" = "0" ] || ! is_numeric_time "$timeName"; then
        return
    fi
    if [ -d "$timeName" ]; then
        echo "[DEBUG] Removing reconstructed root time directory ${timeName}..."
        rm -rf -- "$timeName"
    fi
}

remove_reconstructed_times()
{
    find . -maxdepth 1 -type d -printf '%f\n' \
        | awk '$1 != "0" && $1 ~ /^([0-9]+|[0-9]*\.[0-9]+)([eE][-+]?[0-9]+)?$/ {print}' \
        | while IFS= read -r timeName; do
            remove_reconstructed_time "$timeName"
        done
}

if [ "$#" -gt 0 ] && [ "$1" = "all" ]; then
    echo "[DEBUG] Reconstructing all time steps and converting to VTK..."
    rm -f log.reconstructParMesh log.reconstructPar log.foamToVTK
    rm -rf VTK

    if [ ! -d constant/polyMesh ]; then
        echo "[DEBUG] Base mesh constant/polyMesh is missing. Re-running blockMesh..."
        blockMesh >> log.blockMesh 2>&1 || { echo "[ERROR] blockMesh failed! Exiting."; exit 1; }
        echo "[DEBUG] blockMesh completed successfully."
    fi

    echo "[DEBUG] Running reconstructParMesh -constant..."
    reconstructParMesh -constant >> log.reconstructParMesh 2>&1 || { echo "[ERROR] reconstructParMesh -constant failed! Exiting."; exit 1; }
    
    echo "[DEBUG] Running reconstructParMesh -time 0:..."
    reconstructParMesh -time "0:" >> log.reconstructParMesh 2>&1 || { echo "[ERROR] reconstructParMesh -time 0: failed! Exiting."; exit 1; }

    echo "[DEBUG] Running reconstructPar -time 0:..."
    reconstructPar -time "0:" >> log.reconstructPar 2>&1 || { echo "[ERROR] reconstructPar failed! Exiting."; exit 1; }

    echo "[DEBUG] Running foamToVTK..."
    foam_to_vtk -useTimeName >> log.foamToVTK 2>&1 || { echo "[ERROR] foamToVTK failed! Exiting."; exit 1; }

    remove_reconstructed_times

    echo "[DEBUG] All reconstructions and VTK conversions completed successfully."
    exit 0
fi

if [ "$#" -gt 0 ]; then
    timeName="$1"
else
    timeName=$(
        find processor0 -maxdepth 1 -type d -printf '%f\n' \
            | awk '$1 ~ /^[0-9]/ {print}' \
            | sort -g \
            | tail -1
    )
fi

if [ -z "$timeName" ]; then
    echo "No processor0 time directory found" >&2
    exit 1
fi

echo "[DEBUG] Converting time step ${timeName} to VTK..."

if [ ! -d constant/polyMesh ]; then
    echo "[DEBUG] Base mesh constant/polyMesh is missing. Re-running blockMesh..."
    blockMesh >> log.blockMesh 2>&1 || { echo "[ERROR] blockMesh failed! Exiting."; exit 1; }
    echo "[DEBUG] blockMesh completed successfully."
fi

echo "[DEBUG] Running reconstructParMesh for time step ${timeName}..."
reconstructParMesh -time "$timeName" >> log.reconstructParMesh 2>&1 || { echo "[ERROR] reconstructParMesh failed for time ${timeName}! Exiting."; exit 1; }

echo "[DEBUG] Running reconstructPar for time step ${timeName}..."
reconstructPar -time "$timeName" >> log.reconstructPar 2>&1 || { echo "[ERROR] reconstructPar failed for time ${timeName}! Exiting."; exit 1; }

echo "[DEBUG] Running foamToVTK for time step ${timeName}..."
foam_to_vtk -time "$timeName" -useTimeName >> log.foamToVTK 2>&1 || { echo "[ERROR] foamToVTK failed for time ${timeName}! Exiting."; exit 1; }

remove_reconstructed_time "$timeName"
