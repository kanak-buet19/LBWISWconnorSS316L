#!/bin/bash

set -e

if [ "$#" -gt 0 ] && [ "$1" = "all" ]; then
    rm -f log.reconstructParMesh log.reconstructPar log.foamToVTK
    rm -rf VTK

    reconstructParMesh -constant >> log.reconstructParMesh 2>&1
    reconstructParMesh -time "0:" >> log.reconstructParMesh 2>&1
    reconstructPar -time "0:" >> log.reconstructPar 2>&1
    foamToVTK -useTimeName >> log.foamToVTK 2>&1

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

echo "Converting time ${timeName} to VTK"

reconstructParMesh -time "$timeName" >> log.reconstructParMesh 2>&1
reconstructPar -time "$timeName" >> log.reconstructPar 2>&1
foamToVTK -time "$timeName" -useTimeName >> log.foamToVTK 2>&1
