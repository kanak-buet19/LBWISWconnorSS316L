#!/bin/bash
# Reconstruct the parallel dynamic-mesh (AMR / topology-changing) case and
# convert every written time to legacy VTK (.vtk ASCII).
#
# The mesh differs per processor and per time (dynamicRefineFvMesh), so the
# mesh must be reconstructed per time with reconstructParMesh *before* the
# fields are mapped with reconstructPar.
cd "${0%/*}" || exit 1
set -o pipefail
. "$WM_PROJECT_DIR/bin/tools/RunFunctions"

rm -f log.reconstructParMesh.* log.reconstructPar log.foamToVTK
rm -rf VTK

# 1. Reconstruct the dynamic mesh (constant + all time directories)
reconstructParMesh -constant   2>&1 | tee log.reconstructParMesh.constant
reconstructParMesh -time '0:'  2>&1 | tee log.reconstructParMesh.times

# 2. Map the fields onto the reconstructed per-time meshes
reconstructPar -time '0:'       2>&1 | tee log.reconstructPar

# 3. Convert to legacy VTK
foamToVTK -legacy -time '0:'    2>&1 | tee log.foamToVTK
