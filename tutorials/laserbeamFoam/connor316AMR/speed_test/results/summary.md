# laserbeamFoam — optimum CPU count

**Answer: run this solver on 10 MPI ranks** — one rank per physical core. It had the fastest mean in the sweep (1.09× the 8-rank baseline), and it is the largest rank count at which every rank still gets its own core.

⚠️ **Read the statistics before quoting a number.** Run-to-run noise on this machine is large (±11% at worst). With 4 repeats, the only pairwise difference that clears p < 0.05 is 10 vs 20. The ordering below is a consistent *trend*, not a set of proven differences.

## Setup

- **Case**: `hofmann_scantrack_200W_900mms_r25um` (SS316L, 200 W, 900 mm/s, r = 25 µm)
- **Mesh**: 51 × 36 × 150 = **275,400 cells** base (10 µm), 1 AMR level
- **Host**: Intel Xeon w5-2445 — **10 physical cores / 20 logical CPUs** (2-way hyperthreading), 1 NUMA node, 26.3 MiB L3
- **Method**: every rank count started from a byte-identical initial state and ran for a fixed wall-clock budget of 120 s, repeated **4×**. Field writing was disabled, so the measurement is pure compute + communication with no I/O.
- **Comparison point**: all runs are compared at **timestep 73** — the furthest step every run reached. Simulated time at that step agrees to within 0.001% across all configurations, so each has done identical physics work and the wall time to get there is a like-for-like measurement.
- **Pinning**: ranks were spread across distinct physical cores first (`--cpu-set` in hwloc order); only counts above 10 use hyperthread siblings. This gives each configuration the most favourable placement available.

## Results

Wall time to reach timestep 73, mean of 4 runs ± 95% interval (2×SEM). Lower is better.

| Ranks | Physical cores | Time to step 73 (s) | 95% interval | Speed-up vs 8 | Efficiency | s / timestep | Steps in budget | p vs best |
|------:|---------------:|----------------------------:|:-------------:|-------------:|-----------:|-------------:|----------------:|----------:|
| **8** | 8 | 99.72 ± 2.6 | 95 – 105 | 1.00× | 100% | 1.631 | 85 | 0.091 |
| **10** ⬅ **best** | 10 | 91.47 ± 3.2 | 85 – 98 | 1.09× | 87% | 1.536 | 91 | — |
| **12** *(HT)* | 10 | 96.89 ± 2.1 | 93 – 101 | 1.03× | 69% | 1.609 | 86 | 0.210 |
| **16** *(HT)* | 10 | 105.12 ± 5.8 | 94 – 117 | 0.95× | 47% | 1.678 | 82 | 0.097 |
| **20** *(HT)* | 10 | 107.18 ± 4.1 | 99 – 116 | 0.93× | 37% | 1.714 | 79 | 0.026&nbsp;✅ |

*(HT) = more ranks than physical cores, so at least some cores run two ranks on their two hyperthreads. Efficiency is relative to the 8-rank baseline, which is 100% by construction. `p vs best` is a Welch t-test against 10 ranks; ✅ marks p < 0.05.*

## Verdict

- **Use 10 ranks.** Fastest mean (91.5 s to step 73, 1.09× the 8-rank baseline) and it is exactly one rank per physical core — the natural hardware boundary on this machine.
- **Never go above 10.** This is the firmest result in the sweep: the best hyperthreaded configuration reached only 1.03×, i.e. it never beat plain one-rank-per-core placement, and 16/20 ranks were slower in every repeat while occupying the entire machine.
- **8 and 12 ranks are acceptable fallbacks.** They are within noise of 10; the honest gap is a few percent, not a factor.
- **The knee is at 10 ranks** — the last increment that still bought at least 5% over the previous count.

### How solid is this?

Each configuration was run **4×**, with the sweep order reversed halfway through so that warm-up and turbo-frequency drift over a sweep average out instead of penalising whichever configurations happen to run last. Even so, run-to-run spread reached ±11% of the mean — this is a shared desktop workstation, not a quiet compute node.

**Only the 10-vs-20 comparison(s) reach p < 0.05.** Differences among the leading configurations (8, 10, 12, 16) are *not* statistically significant at n=4: 10 has the best mean and led in most individual repeats, but this data cannot prove it beats its immediate neighbours.

The t-test is not the whole story though: **16, 20 ranks lost to 10 in *every one* of the 4 repeats**, including the reversed-order ones where they ran first and had the thermal advantage. A clean sweep like that is consistent directional evidence even where the t-test is inconclusive, which is why the advice to stay at or below 10 ranks is firm.

That does not weaken the recommendation: 10 ranks has the best point estimate, sits on the hardware boundary, and costs nothing to choose over its neighbours. It does mean you should not quote the 9% margin over 8 ranks as a precise figure. To separate the leaders properly you would need roughly 10 repeats per configuration on an otherwise idle machine (`BENCH_REPEATS=10 BENCH_CORES="8 10 12" ./run_benchmark.sh`).

### Why more ranks stop helping

1. **There are only 10 physical cores.** Rank counts above 10 add no compute hardware — they place two MPI ranks on the two hyperthreads of one core, where the ranks share execution units, L1/L2 cache and the core's share of memory bandwidth. Finite-volume CFD is memory-bandwidth bound, so the sibling thread mostly adds cache pressure and MPI traffic rather than work done.
2. **Subdomains get too small.** At 10 ranks the case is already down to ~27,540 cells per rank. Halo exchange and the pressure solve's global reductions grow as a share of each timestep as subdomains shrink, so the parallel overhead rises even before hyperthreading bites.
3. **Start-up cost grows with rank count** — MPI init and mesh distribution went from 1.2 s at 8 ranks to 2.6 s at 20 ranks, plus ~2.3 s of `decomposePar` in every case.

### Practical note

Leaving 10 logical CPUs free is not waste — 10 ranks is as fast as anything else measured *and* it leaves the machine responsive for the other user's desktop session. If you want to get more total work done, run two studies side by side rather than one wide job: throughput per rank is much better below the physical-core limit than above it.

## Files

- `timings.csv` — aggregated, one row per rank count (all metrics above)
- `per_run.csv` — one row per individual run, including each repeat
- `per_step.csv` — per-timestep elapsed and incremental wall time
- `plots/dashboard.png` — all four charts in one figure
- `plots/scaling_speedup.png`, `plots/efficiency.png`, `plots/time_per_step.png`, `plots/progress_curves.png`
- `raw/` — unmodified solver logs, `decomposePar` logs and per-run metadata
- `../run_benchmark.sh` — the harness; `BENCH_CORES`, `BENCH_DURATION`, `BENCH_REPEATS` are the knobs
