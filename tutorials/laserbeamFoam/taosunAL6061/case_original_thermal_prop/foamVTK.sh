#!/bin/bash

set -e

logDir="${CASE_LOG_DIR:-.}"
mkdir -p "$logDir"

resolveTimeName()
{
    local requested="$1"
    local resolved=""

    if [ -d "processor0/$requested" ]; then
        printf '%s\n' "$requested"
        return 0
    fi

    resolved="$(
        find processor0 -maxdepth 1 -type d -printf '%f\n' 2>/dev/null \
            | awk -v target="$requested" '
                $1 ~ /^[0-9]/ {
                    value = $1 + 0.0
                    diff = value - target
                    if (diff < 0) diff = -diff
                    tolerance = target
                    if (tolerance < 0) tolerance = -tolerance
                    tolerance = tolerance * 1.0e-4
                    if (tolerance < 1.0e-12) tolerance = 1.0e-12
                    if (diff <= tolerance && (best == "" || diff < bestDiff)) {
                        best = $1
                        bestDiff = diff
                    }
                }
                END {
                    if (best != "") print best
                }
            '
    )"

    if [ -n "$resolved" ]; then
        printf '%s\n' "$resolved"
        return 0
    fi

    printf '%s\n' "$requested"
}

waitForTimeDir()
{
    local requested="$1"
    local resolved=""
    local attempt=0

    while [ "$attempt" -lt 20 ]; do
        resolved="$(resolveTimeName "$requested")"
        if [ -d "processor0/$resolved" ]; then
            printf '%s\n' "$resolved"
            return 0
        fi
        sleep 0.25
        attempt=$((attempt + 1))
    done

    echo "No processor0 time directory found for requested time ${requested}" >&2
    return 1
}

processorCount()
{
    find . -maxdepth 1 -type d -name 'processor[0-9]*' -printf '%f\n' 2>/dev/null | wc -l
}

caseName()
{
    basename "$PWD"
}

timeDirs()
{
    find processor0 -maxdepth 1 -type d -printf '%f\n' 2>/dev/null \
        | awk '$1 ~ /^[0-9]/ {print}' \
        | sort -g
}

renameLegacyVTKToTimeName()
{
    local timeName="$1"
    local newestVtk=""
    local dirName=""
    local baseName=""
    local caseName=""
    local index=""
    local serialNum=""

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
}

timeIsComplete()
{
    local timeName="$1"
    local nProcs="$2"
    local completeProcs=0

    completeProcs="$(
        find processor*/"$timeName"/polyMesh/boundary -type f 2>/dev/null \
            | wc -l
    )"

    [ "$completeProcs" -eq "$nProcs" ]
}

waitForCompleteProcessorTime()
{
    local timeName="$1"
    local nProcs=""
    local attempt=0

    nProcs="$(processorCount)"
    if [ "$nProcs" -le 0 ]; then
        echo "No processor directories found" >&2
        return 1
    fi

    while [ "$attempt" -lt 120 ]; do
        if timeIsComplete "$timeName" "$nProcs"; then
            sleep 2
            if timeIsComplete "$timeName" "$nProcs"; then
                return 0
            fi
        fi
        sleep 1
        attempt=$((attempt + 1))
    done

    echo "Processor time ${timeName} did not become complete across ${nProcs} processors" >&2
    return 1
}

if [ "$#" -gt 0 ] && [ "$1" = "all" ]; then
    rm -f "$logDir/log.reconstructParMesh" "$logDir/log.reconstructPar" "$logDir/log.foamToVTK"
    rm -rf VTK

    reconstructParMesh -constant >> "$logDir/log.reconstructParMesh" 2>&1 || true
    reconstructParMesh -time "0:" >> "$logDir/log.reconstructParMesh" 2>&1 || true
    reconstructPar -time "0:" >> "$logDir/log.reconstructPar" 2>&1

    while IFS= read -r timeName
    do
        echo "Converting reconstructed time ${timeName} to legacy VTK"
        foamToVTK -legacy -time "$timeName" >> "$logDir/log.foamToVTK" 2>&1
        renameLegacyVTKToTimeName "$timeName"
    done < <(timeDirs)

    exit 0
fi

if [ "$#" -gt 0 ]; then
    timeName="$(waitForTimeDir "$1")"
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

waitForCompleteProcessorTime "$timeName"
reconstructParMesh -time "0:$timeName" >> "$logDir/log.reconstructParMesh" 2>&1 || true
reconstructPar -time "$timeName" >> "$logDir/log.reconstructPar" 2>&1
foamToVTK -legacy -time "$timeName" >> "$logDir/log.foamToVTK" 2>&1
renameLegacyVTKToTimeName "$timeName"
