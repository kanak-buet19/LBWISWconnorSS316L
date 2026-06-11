#!/bin/bash

set -e

. "$WM_PROJECT_DIR/bin/tools/RunFunctions"

if [ "$#" -gt 0 ] && [ "$1" = "all" ]; then
    rm -f log.reconstructParMesh log.reconstructPar log.foamToVTK
    rm -rf VTK

    runApplication reconstructParMesh -constant
    runApplication reconstructParMesh -time "0:"
    runApplication reconstructPar -time "0:"
    runApplication foamToVTK -legacy

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

reconstructParMesh -time "0:$timeName" >> log.reconstructParMesh 2>&1 || true
reconstructPar -time "$timeName" >> log.reconstructPar 2>&1
foamToVTK -time "$timeName" -legacy >> log.foamToVTK 2>&1

# Rename the newly created index-based .vtk file to contain the actual timename followed by the serial number (Option 1)
newestVtk=$(find VTK -maxdepth 1 -name "*.vtk" ! -name "*$timeName*" -type f -printf '%T@ %p\n' 2>/dev/null | sort -g | tail -1 | cut -d' ' -f2-)
if [ -n "$newestVtk" ]; then
    dirName=$(dirname "$newestVtk")
    baseName=$(basename "$newestVtk")
    caseName=$(echo "$baseName" | sed -E 's/_[0-9]+\.vtk$//')
    if [ "$caseName" = "$baseName" ]; then
        caseName="${baseName%.vtk}"
    fi

    # Calculate 4-digit padded serial number (index starting from 0000)
    index=$(find processor0 -maxdepth 1 -type d -printf '%f\n' | grep -E '^[0-9]' | sort -g | grep -n -F "$timeName" | cut -d':' -f1)
    if [ -n "$index" ]; then
        serialNum=$(printf "%04d" $((index - 1)))
    else
        serialNum="0000"
    fi

    mv "$newestVtk" "${dirName}/${caseName}_${timeName}_${serialNum}.vtk"
    echo "Renamed $newestVtk to ${dirName}/${caseName}_${timeName}_${serialNum}.vtk"
fi
