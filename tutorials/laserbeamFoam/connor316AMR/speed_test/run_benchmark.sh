#!/usr/bin/env bash
#
# Strong-scaling benchmark for laserbeamFoam on a single shared-memory node.
#
# Runs the SAME case (identical mesh, identical initial state) on a range of
# MPI rank counts for a fixed wall-clock budget, then compares how far each
# configuration got. Because every run is stopped by the clock rather than by
# a fixed amount of work, the analysis compares the runs at a COMMON TIMESTEP
# INDEX (see scripts/analyze_speedtest.py).
#
# Host topology assumption (verified on EN4226751L, Xeon w5-2445):
#   10 physical cores, 2 hyperthreads each -> 20 logical CPUs, 1 NUMA node.
#   hwloc logical PU ids: even = first thread of each core, odd = its sibling.
# Ranks are therefore pinned to spread over physical cores first and only then
# start doubling up on hyperthreads.
#
# Usage:
#   ./run_benchmark.sh                       # 8 10 12 16 20, 120 s each
#   BENCH_DURATION=60 ./run_benchmark.sh     # shorter runs
#   BENCH_CORES="4 8 10" ./run_benchmark.sh  # custom rank counts
#
set -u

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
CASE_DIR="${SCRIPT_DIR}/case"
RAW_DIR="${SCRIPT_DIR}/results/raw"

DURATION="${BENCH_DURATION:-120}"      # measured seconds per configuration
WARMUP="${BENCH_WARMUP:-25}"           # discarded warm-up run (page cache, turbo)
REPEATS="${BENCH_REPEATS:-1}"          # repeats of the whole sweep, for noise bars
FIRST_REP="${BENCH_FIRST_REP:-1}"      # start numbering here (to add repeats later)
read -r -a CORES <<< "${BENCH_CORES:-8 10 12 16 20}"

N_PHYS=$(lscpu -e=CORE 2>/dev/null | tail -n +2 | sort -u | wc -l)
N_LOG=$(nproc)

mkdir -p "${RAW_DIR}"

# --- environment -------------------------------------------------------------
if [ -z "${WM_PROJECT_DIR:-}" ]; then
    set +u   # the OpenFOAM bashrc dereferences unset variables
    # shellcheck disable=SC1091
    source /usr/lib/openfoam/openfoam2412/etc/bashrc >/dev/null 2>&1
    set -u
fi
command -v laserbeamFoam >/dev/null || { echo "[ERROR] laserbeamFoam not on PATH"; exit 1; }

# --- helpers -----------------------------------------------------------------

# Build an hwloc-logical PU list that fills distinct physical cores before it
# starts using hyperthread siblings.
cpuset_for() {
    local n=$1 list=() i
    for ((i=0; i<N_PHYS && ${#list[@]}<n; i++)); do list+=($((i*2))); done
    for ((i=0; i<N_PHYS && ${#list[@]}<n; i++)); do list+=($((i*2+1))); done
    local IFS=,
    echo "${list[*]}"
}

# Remove every artefact of a previous run but keep constant/polyMesh and 0/,
# so all configurations start from a byte-identical initial state.
reset_case() {
    ( cd "${CASE_DIR}" || exit 1
      rm -rf processor*
      # drop any written time directory other than 0
      for d in [0-9]*; do
          [ -d "$d" ] && [ "$d" != "0" ] && rm -rf "$d"
      done
      rm -rf postProcessing
    )
}

run_one() {
    local n=$1 secs=$2 tag=$3
    local cs; cs=$(cpuset_for "$n")

    reset_case
    ( cd "${CASE_DIR}" && foamDictionary -entry numberOfSubdomains -set "$n" \
        system/decomposeParDict >/dev/null 2>&1 ) || { echo "[ERROR] could not set subdomains"; return 1; }

    # --- decomposition (measured separately; not part of the solve time) ---
    local d0 d1 dt_dec
    d0=$(date +%s.%N)
    ( cd "${CASE_DIR}" && decomposePar -force > "${RAW_DIR}/log.decomposePar.np${n}${tag}" 2>&1 ) \
        || { echo "[ERROR] decomposePar failed for np=$n"; return 1; }
    d1=$(date +%s.%N)
    dt_dec=$(awk -v a="$d0" -v b="$d1" 'BEGIN{printf "%.2f", b-a}')

    # --- solve (fixed wall-clock budget) ---
    local s0 s1 dt_solve
    s0=$(date +%s.%N)
    ( cd "${CASE_DIR}" && timeout -s TERM "${secs}" \
        mpirun --use-hwthread-cpus --cpu-set "${cs}" --bind-to cpu-list:ordered \
               -np "$n" laserbeamFoam -parallel \
        > "${RAW_DIR}/log.solve.np${n}${tag}" 2>&1 )
    s1=$(date +%s.%N)
    dt_solve=$(awk -v a="$s0" -v b="$s1" 'BEGIN{printf "%.2f", b-a}')

    # make sure nothing survived the timeout
    pkill -f "laserbeamFoam -parallel" 2>/dev/null
    sleep 2

    local steps
    steps=$(grep -c "^Time = " "${RAW_DIR}/log.solve.np${n}${tag}" 2>/dev/null || echo 0)

    # record metadata alongside the log
    cat > "${RAW_DIR}/meta.np${n}${tag}.txt" <<EOF
ranks=${n}
cpuset=${cs}
decompose_s=${dt_dec}
solve_wall_s=${dt_solve}
timesteps=${steps}
budget_s=${secs}
EOF
    printf '  np=%-3s decompose=%6ss  solve=%7ss  timesteps=%s\n' \
        "$n" "$dt_dec" "$dt_solve" "$steps"
}

# --- main --------------------------------------------------------------------
echo "=================================================================="
echo " laserbeamFoam strong-scaling benchmark"
echo "=================================================================="
echo " host           : $(hostname)"
echo " CPU            : $(lscpu | sed -n 's/^Model name: *//p')"
echo " physical cores : ${N_PHYS}   logical CPUs: ${N_LOG}"
echo " case           : ${CASE_DIR}"
echo " base mesh      : $(grep -oE 'nCells: *[0-9]+' "${CASE_DIR}/log.blockMesh" | head -1)"
echo " rank counts    : ${CORES[*]}"
echo " budget per run : ${DURATION} s   (warm-up ${WARMUP} s)"
echo " started        : $(date -Is)"
echo "------------------------------------------------------------------"

if command -v uptime >/dev/null; then echo " load before: $(uptime)"; fi

if [ "${WARMUP}" -gt 0 ]; then
    echo "[warm-up] discarded ${WARMUP}s run at np=${CORES[0]} ..."
    run_one "${CORES[0]}" "${WARMUP}" ".warmup" >/dev/null 2>&1
    rm -f "${RAW_DIR}"/*.warmup* 2>/dev/null
    echo "[warm-up] done"
fi

echo "[measured runs]"
for ((rep=FIRST_REP; rep<FIRST_REP+REPEATS; rep++)); do
    echo " -- repeat ${rep} --"
    for n in "${CORES[@]}"; do
        run_one "$n" "${DURATION}" ".rep${rep}"
    done
done

reset_case

echo "------------------------------------------------------------------"
echo " finished: $(date -Is)"
echo " raw logs: ${RAW_DIR}"
echo " next    : python3 scripts/analyze_speedtest.py"
echo "=================================================================="
