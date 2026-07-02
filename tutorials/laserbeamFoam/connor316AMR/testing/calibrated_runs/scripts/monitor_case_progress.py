#!/usr/bin/env python3
"""Live dashboard for parallel validation case progress."""

from __future__ import annotations

import argparse
import re
import signal
import time
from dataclasses import dataclass
from pathlib import Path

from rich.console import Console
from rich.live import Live
from rich.progress import BarColumn, Progress, TextColumn
from rich.table import Table


ROOT = Path(__file__).resolve().parents[1]
RUNNING = True


def stop_dashboard(_signum, _frame) -> None:
    global RUNNING
    RUNNING = False


@dataclass
class CaseState:
    name: str
    end_time_s: float
    time_s: float = 0.0
    dt_s: float = 0.0
    max_co: float = 0.0
    t_max_k: float = 0.0
    p_vap_kpa: float = 0.0
    u_metal: float = 0.0
    u_gas: float = 0.0
    log_mtime: float = 0.0
    status: str = "starting"


def short_label(name: str) -> str:
    match = re.search(r"hofmann_scantrack_(\d+)W_(\d+)mms", name)
    if match:
        return f"{match.group(1)}W {match.group(2)}mm/s"
    return name


def parse_end_time(case_dir: Path) -> float:
    control = case_dir / "system" / "controlDict"
    if not control.exists():
        return 1.0
    match = re.search(r"^\s*endTime\s+([0-9.eE+-]+)\s*;", control.read_text(), re.MULTILINE)
    return float(match.group(1)) if match else 1.0


def read_tail(path: Path, max_bytes: int = 120_000) -> str:
    if not path.exists():
        return ""
    size = path.stat().st_size
    with path.open("rb") as handle:
        if size > max_bytes:
            handle.seek(size - max_bytes)
        return handle.read().decode("utf-8", errors="replace")


def last_float(pattern: str, text: str) -> float | None:
    matches = re.findall(pattern, text, flags=re.MULTILINE)
    if not matches:
        return None
    value = matches[-1]
    if isinstance(value, tuple):
        value = value[-1]
    try:
        return float(value)
    except ValueError:
        return None


def update_state(state: CaseState, case_dir: Path) -> None:
    log_path = case_dir / "log.laserbeamFoam"
    validation_log = case_dir / "log.validationAllrun"
    text = read_tail(log_path)
    if not text:
        text = read_tail(validation_log)

    if log_path.exists():
        state.log_mtime = log_path.stat().st_mtime
    elif validation_log.exists():
        state.log_mtime = validation_log.stat().st_mtime

    if not text:
        return

    value = last_float(r"^Time =\s+([0-9.eE+-]+)", text)
    if value is not None:
        state.time_s = value
        state.status = "running"

    value = last_float(r"deltaT =\s+([0-9.eE+-]+)", text)
    if value is not None:
        state.dt_s = value

    value = last_float(r"Courant Number mean:.*?max:\s+([0-9.eE+-]+)", text)
    if value is not None:
        state.max_co = value

    value = last_float(r"TMax =\s+([0-9.eE+-]+)", text)
    if value is not None:
        state.t_max_k = value

    value = last_float(r"pVapMax =\s+([0-9.eE+-]+)", text)
    if value is not None:
        state.p_vap_kpa = value / 1000.0

    match = re.findall(r"maxU_metal =\s+([0-9.eE+-]+).*?maxU_gas =\s+([0-9.eE+-]+)", text)
    if match:
        try:
            state.u_metal = float(match[-1][0])
            state.u_gas = float(match[-1][1])
        except ValueError:
            pass

    if "Simulation completed successfully" in text:
        state.status = "done"
    if any(token in text.lower() for token in ("fatal", "segmentation fault", "foam exiting", "aborted")):
        state.status = "failed"


def progress_bar(percent: float) -> Progress:
    progress = Progress(
        TextColumn(""),
        BarColumn(bar_width=24),
        TextColumn("{task.percentage:>3.0f}%"),
        expand=False,
    )
    progress.add_task("", total=100.0, completed=max(0.0, min(100.0, percent)))
    return progress


def render(states: list[CaseState], cores: int) -> Table:
    table = Table(
        title=f"Parallel validation batch - {len(states)} case(s), {cores} cores each",
        expand=True,
    )
    table.add_column("Case", no_wrap=True)
    table.add_column("Progress", justify="center")
    table.add_column("Time", justify="right")
    table.add_column("dt", justify="right")
    table.add_column("maxCo", justify="right")
    table.add_column("TMax", justify="right")
    table.add_column("pVap", justify="right")
    table.add_column("U metal/gas", justify="right")
    table.add_column("Status", justify="center")

    now = time.time()
    for state in states:
        percent = 100.0 * state.time_s / state.end_time_s if state.end_time_s > 0 else 0.0
        status = state.status
        if status == "running" and state.log_mtime and now - state.log_mtime > 90.0:
            status = "quiet"
        style = {
            "starting": "yellow",
            "running": "green",
            "quiet": "yellow",
            "done": "cyan",
            "failed": "red",
        }.get(status, "white")
        table.add_row(
            short_label(state.name),
            progress_bar(percent),
            f"{state.time_s * 1e3:.3f}/{state.end_time_s * 1e3:.3f} ms",
            f"{state.dt_s:.2e}" if state.dt_s > 0 else "-",
            f"{state.max_co:.2f}" if state.max_co > 0 else "-",
            f"{state.t_max_k:.0f} K" if state.t_max_k > 0 else "-",
            f"{state.p_vap_kpa:.1f} kPa" if state.p_vap_kpa > 0 else "-",
            f"{state.u_metal:.2g}/{state.u_gas:.2g}" if state.u_metal > 0 or state.u_gas > 0 else "-",
            f"[{style}]{status}[/{style}]",
        )
    return table


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--case", action="append", required=True, help="Case name under cases/")
    parser.add_argument("--cores", type=int, required=True)
    parser.add_argument("--interval", type=float, default=2.0)
    args = parser.parse_args()
    signal.signal(signal.SIGTERM, stop_dashboard)
    signal.signal(signal.SIGINT, stop_dashboard)

    states = [
        CaseState(name=name, end_time_s=parse_end_time(ROOT / "cases" / name))
        for name in args.case
    ]
    console = Console()
    with Live(render(states, args.cores), console=console, refresh_per_second=4, transient=False) as live:
        while RUNNING:
            for state in states:
                update_state(state, ROOT / "cases" / state.name)
            live.update(render(states, args.cores))
            time.sleep(args.interval)


if __name__ == "__main__":
    main()
