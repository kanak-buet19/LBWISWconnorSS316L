#!/usr/bin/env python3
"""
Live residual watcher for laserbeamFoam.
Parses log every REFRESH_SEC, saves PNG + CSV each iteration.

Usage:
    python3 residuals/plot_residuals.py              # auto-finds log, loops forever
    python3 residuals/plot_residuals.py /path/to/log
    python3 residuals/plot_residuals.py --once       # parse once and exit (post-run)
"""

import re
import sys
import time
import argparse
from pathlib import Path
from collections import defaultdict

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

LOG_NAME    = "log.laserbeamFoam"
CASE_DIR    = Path(__file__).parent.parent   # scripts/../ = case root
OUT_DIR     = CASE_DIR / "residuals"
PNG_PATH    = OUT_DIR / "residuals.png"
CSV_PATH    = OUT_DIR / "residuals.csv"
REFRESH_SEC = 5

RE_ENDTIME  = re.compile(r'^\s*endTime\s+([\d.eE+\-]+)')

# ── patterns (from real log) ─────────────────────────────────────────────────
RE_TIME    = re.compile(r'^Time\s*=\s*([\d.eE+\-]+)')
RE_SOLVER  = re.compile(r'Solving for (\S+),\s*Initial residual\s*=\s*([\d.eE+\-]+)')
RE_EPSILON = re.compile(r'Correcting epsilon1,\s*mean residual\s*=\s*[\d.eE+\-]+,\s*max residual\s*=\s*([\d.eE+\-]+)')
# anchored to avoid matching "Interface Courant Number"
RE_CO      = re.compile(r'^Courant Number mean:\s*([\d.eE+\-]+)\s+max:\s*([\d.eE+\-]+)')
RE_ICO     = re.compile(r'^Interface Courant Number mean:\s*([\d.eE+\-]+)\s+max:\s*([\d.eE+\-]+)')
RE_DT      = re.compile(r'^deltaT\s*=\s*([\d.eE+\-]+)')
RE_TMAX_PVAP = re.compile(r'^TMax\s*=\s*([\d.eE+\-]+),\s*pVapMax\s*=\s*([\d.eE+\-]+)')
# ─────────────────────────────────────────────────────────────────────────────

STYLE = {
    'T':            dict(color='#e63946', lw=1.2, label='T'),
    'alpha.metal':  dict(color='#457b9d', lw=1.2, label='alpha.metal'),
    'p_rgh':        dict(color='#2d6a4f', lw=1.2, label='p_rgh'),
    'epsilon1_max': dict(color='#9b2226', lw=1.4, ls='--', label='epsilon1 max resid'),
    'TMax':         dict(color='#9d4edd', lw=1.2, label='TMax'),
    'pVapMax':      dict(color='#0077b6', lw=1.2, label='pVapMax (kPa)'),
}


def read_end_time(log_path):
    """Read endTime from system/controlDict relative to case root."""
    ctrl = CASE_DIR / "system" / "controlDict"
    if ctrl.exists():
        with open(ctrl, 'r', errors='replace') as f:
            for line in f:
                m = RE_ENDTIME.match(line)
                if m:
                    return float(m.group(1))
    return None


def parse_log(path):
    data    = defaultdict(list)
    courant = []
    icourant = []
    deltat  = []
    current_time = None
    seen = set()
    eps1_this_step = None   # tracks last (converged) epsilon1 max per timestep

    with open(path, 'r', errors='replace') as f:
        for line in f:
            m = RE_TIME.match(line)
            if m:
                # flush last epsilon1 value from previous timestep
                if current_time is not None and eps1_this_step is not None:
                    data['epsilon1_max'].append((current_time, eps1_this_step))
                current_time = float(m.group(1))
                seen = set()
                eps1_this_step = None
                continue

            if current_time is None:
                continue

            m = RE_DT.match(line.lstrip())
            if m:
                deltat.append((current_time, float(m.group(1))))
                continue

            m = RE_TMAX_PVAP.match(line.lstrip())
            if m:
                data['TMax'].append((current_time, float(m.group(1))))
                data['pVapMax'].append((current_time, float(m.group(2))))
                continue

            m = RE_SOLVER.search(line)
            if m:
                field = m.group(1).rstrip('.,')
                if field not in seen:
                    data[field].append((current_time, float(m.group(2))))
                    seen.add(field)
                continue

            m = RE_EPSILON.search(line)
            if m:
                eps1_this_step = float(m.group(1))   # overwrite → keeps last (converged) value
                continue

            m = RE_CO.match(line.lstrip())
            if m:
                courant.append((current_time, float(m.group(2))))
                continue

            m = RE_ICO.match(line.lstrip())
            if m:
                icourant.append((current_time, float(m.group(2))))

    # flush last timestep's epsilon1 (file may end mid-step)
    if current_time is not None and eps1_this_step is not None:
        data['epsilon1_max'].append((current_time, eps1_this_step))

    return data, courant, icourant, deltat


def axis_limits(data, courant, icourant):
    """Derive consistent y-axis limits from data."""
    # y residuals: floor at min seen, ceil at 2
    all_resid = [r for f in ['T', 'alpha.metal', 'p_rgh']
                 for _, r in data[f] if r > 0]
    ymin_r = min(all_resid) * 0.1 if all_resid else 1e-13
    ymin_r = min(ymin_r, 1e-13)
    ylim_resid = (ymin_r, 2.0)

    # y epsilon1: always show full [1e-16, 2] so reference lines visible
    ylim_eps = (1e-16, 2.0)

    ylim_co = (0.0, 2.0)

    return ylim_resid, ylim_eps, ylim_co


def _us(pairs):
    """Convert list of (t_s, val) → (t_µs, val)."""
    return [(t * 1e6, v) for t, v in pairs]


def save_plot(data, courant, icourant, deltat, end_time):
    ylim_resid, ylim_eps, ylim_co = axis_limits(data, courant, icourant)

    all_times = [t for v in data.values() for t, _ in v]
    last_t = max(all_times) if all_times else 0.0
    now_str = time.strftime('%Y-%m-%d %H:%M:%S')

    fig, axes = plt.subplots(2, 3, figsize=(18, 10))
    fig.suptitle(
        f"laserbeamFoam — live residuals & physical extrema    |    "
        f"sim time: {last_t*1e6:.2f} µs    |    plotted: {now_str}",
        fontsize=11
    )

    # [0,0]: solver residuals
    ax = axes[0, 0]
    for field in ['T', 'alpha.metal', 'p_rgh']:
        if data[field]:
            t, r = zip(*_us(data[field]))
            ax.semilogy(t, r, **STYLE[field])
    ax.set_ylim(ylim_resid)
    ax.set_ylabel('Initial residual')
    ax.set_xlabel('Time (µs)')
    ax.set_title('Solver residuals')
    ax.legend(fontsize=8, loc='upper right')
    ax.grid(True, which='both', ls=':', alpha=0.4)

    # [0,1]: epsilon1 max residual
    ax = axes[0, 1]
    if data['epsilon1_max']:
        t, r = zip(*_us(data['epsilon1_max']))
        r = [max(v, 1e-16) for v in r]   # clamp zeros for semilogy
        ax.semilogy(t, r, **STYLE['epsilon1_max'])
    ax.axhline(1.0,  color='red',   ls=':',  lw=1.2, label='ping-pong (=1)')
    ax.axhline(1e-6, color='green', ls=':',  lw=1.2, label='epsilonTol 1e-6')
    ax.set_ylim(ylim_eps)
    ax.set_ylabel('max residual')
    ax.set_xlabel('Time (µs)')
    ax.set_title('epsilon1 (melt-loop) — must stay below 1e-6')
    ax.legend(fontsize=8, loc='upper right')
    ax.grid(True, which='both', ls=':', alpha=0.4)

    # [0,2]: Maximum Temperature
    ax = axes[0, 2]
    if data['TMax']:
        t, tmax = zip(*_us(data['TMax']))
        ax.plot(t, tmax, **STYLE['TMax'])
    ax.set_ylabel('Temperature (K)')
    ax.set_xlabel('Time (µs)')
    ax.set_title('Maximum Temperature')
    ax.legend(fontsize=8, loc='upper right')
    ax.grid(True, ls=':', alpha=0.4)

    # [1,0]: Courant
    ax = axes[1, 0]
    if courant:
        t, c = zip(*_us(courant))
        ax.plot(t, c, color='#6a0572', lw=1.2, label='Co max')
    if icourant:
        t, c = zip(*_us(icourant))
        ax.plot(t, c, color='#b5838d', lw=0.9, ls='--', label='Interface Co max')
    ax.axhline(0.5,  color='orange', ls='--', lw=1, label='maxCo 0.5')
    ax.axhline(0.25, color='green',  ls=':',  lw=1, label='target 0.25')
    ax.set_ylim(ylim_co)
    ax.set_ylabel('Courant number')
    ax.set_xlabel('Time (µs)')
    ax.set_title('Courant number')
    ax.legend(fontsize=8, loc='upper right')
    ax.grid(True, which='both', ls=':', alpha=0.4)

    # [1,1]: deltaT
    ax = axes[1, 1]
    if deltat:
        t, dt = zip(*_us(deltat))
        ax.semilogy(t, dt, color='#f4a261', lw=1.2, label='deltaT')
    ax.set_ylabel('deltaT (s)')
    ax.set_xlabel('Time (µs)')
    ax.set_title('Timestep size')
    ax.legend(fontsize=8, loc='upper right')
    ax.grid(True, which='both', ls=':', alpha=0.4)

    # [1,2]: Maximum Recoil Pressure
    ax = axes[1, 2]
    if data['pVapMax']:
        t, pmax = zip(*_us(data['pVapMax']))
        pmax_kpa = [p / 1000.0 for p in pmax]
        ax.semilogy(t, pmax_kpa, **STYLE['pVapMax'])
    ax.set_ylabel('Recoil Pressure (kPa)')
    ax.set_xlabel('Time (µs)')
    ax.set_title('Maximum Recoil Pressure')
    ax.legend(fontsize=8, loc='upper right')
    ax.grid(True, which='both', ls=':', alpha=0.4)

    plt.tight_layout(rect=[0, 0, 1, 0.96])
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    fig.savefig(PNG_PATH, dpi=150, bbox_inches='tight')
    plt.close(fig)


def save_csv(data, courant, icourant, deltat):
    lookup = defaultdict(dict)
    for field, vals in data.items():
        for t, r in vals:
            lookup[t][field] = r
    for t, c in courant:
        lookup[t]['Co_max'] = c
    for t, c in icourant:
        lookup[t]['ICo_max'] = c
    for t, dt in deltat:
        lookup[t]['deltaT'] = dt

    fields = sorted(set(k for d in lookup.values() for k in d))
    with open(CSV_PATH, 'w') as f:
        f.write('time,' + ','.join(fields) + '\n')
        for t in sorted(lookup):
            row = [str(t)] + [str(lookup[t].get(k, '')) for k in fields]
            f.write(','.join(row) + '\n')


def find_log():
    candidates = [
        CASE_DIR / LOG_NAME,
        Path(LOG_NAME),
    ]
    return next((p for p in candidates if p.exists()), None)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('logfile', nargs='?', default=None)
    parser.add_argument('--once', action='store_true', help='parse once and exit')
    args = parser.parse_args()

    if args.logfile:
        log_path = Path(args.logfile)
    else:
        log_path = find_log()

    # If log doesn't exist yet (watcher started before solver), wait for it
    if log_path is None or not log_path.exists():
        if args.once:
            print(f"Log not found. Pass path as argument.")
            sys.exit(1)
        # loop mode: wait up to 10 min for log to appear
        print(f"Waiting for {LOG_NAME} ...")
        deadline = time.time() + 600
        while time.time() < deadline:
            time.sleep(REFRESH_SEC)
            log_path = find_log()
            if log_path is not None and log_path.exists():
                break
        else:
            print(f"Timed out waiting for log.")
            sys.exit(1)

    end_time = read_end_time(log_path)

    print(f"Watching: {log_path.resolve()}")
    print(f"PNG: {PNG_PATH.resolve()}")
    print(f"endTime: {end_time}")
    if not args.once:
        print(f"Refreshing every {REFRESH_SEC}s. Ctrl-C to stop.")

    def run_once():
        data, courant, icourant, deltat = parse_log(log_path)
        if data or courant:
            save_plot(data, courant, icourant, deltat, end_time)
            save_csv(data, courant, icourant, deltat)
            n_times = len(set(t for v in data.values() for t, _ in v))
            print(f"[{time.strftime('%H:%M:%S')}] {n_times} timesteps → {PNG_PATH.name}")

    if args.once:
        run_once()
        return

    try:
        while True:
            run_once()
            time.sleep(REFRESH_SEC)
    except KeyboardInterrupt:
        print("\nFinal save...")
        run_once()
        print("Done.")


if __name__ == '__main__':
    main()
