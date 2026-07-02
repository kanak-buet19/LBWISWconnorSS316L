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

normalize_vtk_outputs()
{
    "${PYTHON:-python3}" <<'PY'
from __future__ import annotations

import json
import re
from collections import defaultdict
from pathlib import Path


def time_suffix(value: float | str) -> str:
    if isinstance(value, str):
        return value
    return f"{value:.12g}"


def rename_logged_foam_vtk() -> None:
    log_path = Path("log.foamToVTK")
    if not log_path.exists():
        return

    current_time = None
    time_re = re.compile(r"^Time:\s+([-+0-9.eE]+)")
    path_re = re.compile(r'^\s+(?:Internal|Boundary)\s+:\s+"([^"]+_)(\d+)(\.vtk)"')

    for line in log_path.read_text(errors="ignore").splitlines():
        time_match = time_re.match(line)
        if time_match:
            current_time = time_match.group(1)
            continue

        path_match = path_re.match(line)
        if not path_match or current_time is None:
            continue

        source = Path("".join(path_match.groups()))
        target = Path(f"{path_match.group(1)}{current_time}{path_match.group(3)}")
        if source == target or not source.exists():
            continue
        if target.exists():
            target.unlink()
        source.rename(target)


def parse_time_from_name(path: Path) -> float | None:
    match = re.search(r"_([-+0-9.eE]+)\.vtk$", path.name)
    if not match:
        return None
    try:
        return float(match.group(1))
    except ValueError:
        return None


def write_series(series_path: Path, entries: list[tuple[float, Path]], base_dir: Path) -> None:
    entries = sorted(entries, key=lambda item: item[0])
    data = {
        "file-series-version": "1.0",
        "files": [
            {"name": path.relative_to(base_dir).as_posix(), "time": time}
            for time, path in entries
        ],
    }
    series_path.write_text(json.dumps(data, indent=2) + "\n")


def write_vtk_series() -> None:
    vtk_dir = Path("VTK")
    if not vtk_dir.exists():
        return

    groups: dict[Path, list[tuple[float, Path]]] = defaultdict(list)
    for path in vtk_dir.rglob("*.vtk"):
        time_value = parse_time_from_name(path)
        if time_value is None:
            continue
        prefix = re.sub(r"_[-+0-9.eE]+\.vtk$", "", path.name)
        groups[path.with_name(prefix + ".vtk.series")].append((time_value, path))

    for series_path, entries in groups.items():
        write_series(series_path, entries, series_path.parent)


def normalize_ray_vtks() -> None:
    rays_dir = Path("VTKs")
    if not rays_dir.exists():
        return

    for series_path in rays_dir.glob("*.vtk.series"):
        try:
            data = json.loads(series_path.read_text())
        except (OSError, json.JSONDecodeError):
            continue

        changed = False
        for entry in data.get("files", []):
            try:
                time_value = float(entry["time"])
            except (KeyError, TypeError, ValueError):
                continue

            source = series_path.parent / str(entry.get("name", ""))
            prefix = re.sub(r"_[-+0-9.eE]+\.vtk$", "", source.name)
            target = series_path.parent / f"{prefix}_{time_suffix(time_value)}.vtk"

            if source.exists() and source != target:
                if target.exists():
                    target.unlink()
                source.rename(target)
                changed = True
            if target.exists() and entry.get("name") != target.name:
                entry["name"] = target.name
                changed = True

        if changed:
            series_path.write_text(json.dumps(data, indent=2) + "\n")


rename_logged_foam_vtk()
write_vtk_series()
normalize_ray_vtks()
PY
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
    foam_to_vtk >> log.foamToVTK 2>&1 || { echo "[ERROR] foamToVTK failed! Exiting."; exit 1; }
    normalize_vtk_outputs

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
foam_to_vtk -time "$timeName" >> log.foamToVTK 2>&1 || { echo "[ERROR] foamToVTK failed for time ${timeName}! Exiting."; exit 1; }
normalize_vtk_outputs

remove_reconstructed_time "$timeName"
