# speed_test — how many MPI ranks should laserbeamFoam use?

Strong-scaling benchmark of `laserbeamFoam` on this workstation, using one real
validation case. **Result: use 10 ranks.** See `results/summary.md` for the full
write-up, numbers and caveats.

## Layout

```
speed_test/
  run_benchmark.sh              driver — runs the sweep, writes results/raw/
  scripts/analyze_speedtest.py  parses the logs, writes CSVs, plots and summary.md
  case/                         the benchmark case (meshed once, reused by every run)
  results/
    summary.md                  ← read this
    timings.csv                 aggregated, one row per rank count
    per_run.csv                 one row per individual run (every repeat)
    per_step.csv                per-timestep wall clock
    plots/dashboard.png         all four charts in one figure
    plots/*.png                 individual charts
    raw/                        unmodified solver + decomposePar logs, run metadata
    benchmark_run*.log          driver console output
```

## Re-running it

```bash
of2412                                   # OpenFOAM env + venv

./run_benchmark.sh                       # default: 8 10 12 16 20, 120 s each, 1 repeat
python scripts/analyze_speedtest.py      # regenerate CSVs, plots and summary.md
```

Knobs (environment variables):

| Variable | Default | Meaning |
|---|---|---|
| `BENCH_CORES` | `8 10 12 16 20` | rank counts to test, in the order they run |
| `BENCH_DURATION` | `120` | measured wall-clock seconds per run |
| `BENCH_WARMUP` | `25` | discarded warm-up run (page cache, turbo) |
| `BENCH_REPEATS` | `1` | repeats of the whole sweep |
| `BENCH_FIRST_REP` | `1` | starting repeat number — use it to *add* repeats later |

Adding two more repeats to an existing dataset, in reversed order so that
warm-up drift does not bias whichever configuration runs last:

```bash
BENCH_WARMUP=25 BENCH_FIRST_REP=5 BENCH_REPEATS=2 \
  BENCH_CORES="20 16 12 10 8" ./run_benchmark.sh
python scripts/analyze_speedtest.py
```

## How the measurement works

Every run starts from a byte-identical initial state and is stopped by the
clock, not by a fixed amount of work. Runs therefore cannot be compared on
"how long did it take" — instead the analysis picks
`K = min(timesteps completed over all runs)` and compares the wall time each
configuration needed to reach **timestep K**. Simulated time at step K agrees
across configurations to within 0.001%, so every run has done the same physics
work and the comparison is like-for-like.

Other things the harness controls for:

- **Writing is disabled** (`writeInterval` set very large), so the numbers are
  pure compute + MPI with no I/O noise.
- **Rank pinning** spreads ranks over distinct physical cores first and only
  then uses hyperthread siblings, via `--cpu-set` in hwloc logical order. This
  gives every configuration the most favourable placement available.
- **`decomposePar` is timed separately** and excluded from the solve time.
- **Repeats and sweep order.** The host has 10 physical cores and 20 logical
  CPUs; run-to-run noise is large because it is a shared desktop machine.
  Repeats 1–2 ran 8→20 and repeats 3–4 ran 20→8, so warm-up and turbo drift
  average out instead of systematically penalising whichever ranks run last.
  This mattered: 16 ranks scored 17% better when it ran first than when it ran
  last.

## Caveat worth knowing

Run-to-run spread reaches ±11%. With 4 repeats only the 10-vs-20 difference
clears p < 0.05 — the leading configurations (8, 10, 12) are not separable by
this data. 10 ranks has the best mean and sits on the hardware boundary
(one rank per physical core), so it is the right default, but do not quote the
margin over 8 ranks as a precise number. To separate the leaders properly:

```bash
BENCH_REPEATS=10 BENCH_CORES="8 10 12" ./run_benchmark.sh
```
