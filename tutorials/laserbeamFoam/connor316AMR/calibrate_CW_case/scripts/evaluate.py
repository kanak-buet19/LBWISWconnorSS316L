"""Melt-pool stability + error evaluation for calibration.

Reads `<case>/post-processing-data/vtu_meltpool_geometry.csv` (written live by
analyze_meltpool_vtu.py during the run) and reduces the time series of
meltPoolDepth_um / meltPoolWidth_um to a single 'stable' value, then computes
relative error vs the experimental target.
"""

from __future__ import annotations

import csv
import math
from pathlib import Path

CSV_REL = "post-processing-data/vtu_meltpool_geometry.csv"


def _finite(x: str | float) -> float | None:
    try:
        v = float(x)
    except (TypeError, ValueError):
        return None
    return v if math.isfinite(v) else None


def read_series(case_dir: Path) -> list[dict]:
    """Return [{time, depth, width}, ...] for rows with finite depth AND width."""
    csv_path = Path(case_dir) / CSV_REL
    if not csv_path.exists():
        return []
    rows = []
    with csv_path.open(newline="") as fh:
        for r in csv.DictReader(fh):
            t = _finite(r.get("time"))
            d = _finite(r.get("meltPoolDepth_um"))
            w = _finite(r.get("meltPoolWidth_um"))
            if t is not None and d is not None and w is not None:
                rows.append({"time": t, "depth": d, "width": w})
    rows.sort(key=lambda x: x["time"])
    return rows


def _stable_tail(values: list[float], tol: float,
                 tail_frac: float = 0.25) -> tuple[float, bool, int]:
    """Check if the recent tail of the series is stable (used for early-abort during partial data).

    Uses the last tail_frac fraction of available points (min 2).
    Returns (mean_of_tail, converged, n_tail).
    """
    n = len(values)
    if n == 0:
        return float("nan"), False, 0
    if n == 1:
        return values[0], False, 1
    tail_n = max(2, int(n * tail_frac))
    tail = values[-tail_n:]
    mean_v = sum(tail) / tail_n
    if mean_v == 0:
        return mean_v, False, tail_n
    var = (max(tail) - min(tail)) / abs(mean_v)
    return mean_v, var < tol, tail_n


def _stable(values: list[float], tol: float) -> tuple[float, bool, int]:
    """Find the most stable window in the middle of the series.

    Skips first and last 20% of points (transient + end effects).
    Slides all possible windows of size >= 2 over the middle 60%.
    Picks the window with lowest (max-min)/mean variation.
    Returns (mean_of_best_window, converged, n_stable).
    converged = True when n_stable >= 2 AND variation < tol.
    """
    n = len(values)
    if n == 0:
        return float("nan"), False, 0
    if n == 1:
        return values[0], False, 1

    skip = max(1, n // 5)          # 20% from each end
    i0 = skip
    i1 = n - skip                  # exclusive

    if i1 - i0 < 2:               # middle too short, use all
        i0, i1 = 0, n

    best_var = float("inf")
    best_mean = sum(values[i0:i1]) / (i1 - i0)
    best_n = i1 - i0

    for size in range(2, i1 - i0 + 1):
        for i in range(i0, i1 - size + 1):
            w = values[i : i + size]
            mean_v = sum(w) / size
            if mean_v == 0:
                continue
            var = (max(w) - min(w)) / abs(mean_v)
            if var < best_var:
                best_var = var
                best_mean = mean_v
                best_n = size

    return best_mean, best_n >= 2 and best_var < tol, best_n


def evaluate(case_dir: Path, exp_depth_um: float, exp_width_um: float,
             stability_tol: float = 0.10) -> dict:
    series = read_series(case_dir)
    n = len(series)
    depths = [r["depth"] for r in series]
    widths = [r["width"] for r in series]

    depth, d_conv, d_n = _stable(depths, stability_tol)
    width, w_conv, w_n = _stable(widths, stability_tol)

    def rel(sim, exp):
        return float("nan") if not exp else (sim - exp) / exp

    depth_err = rel(depth, exp_depth_um)
    width_err = rel(width, exp_width_um)

    # Aspect ratio error prevents compensating errors (e.g., depth +10%,
    # width -10%) from scoring well just because abs errors cancel in mean.
    exp_ar = exp_depth_um / exp_width_um if exp_width_um else float("nan")
    sim_ar = depth / width if width else float("nan")
    ar_err = (abs(sim_ar - exp_ar) / exp_ar
              if exp_ar and sim_ar == sim_ar and exp_ar == exp_ar
              else float("nan"))

    errs = [e for e in (abs(depth_err), abs(width_err), ar_err) if e == e]
    case_error = sum(errs) / len(errs) if errs else float("nan")

    return {
        "n_points": n,
        "sim_depth_um": depth, "sim_width_um": width,
        "exp_depth_um": exp_depth_um, "exp_width_um": exp_width_um,
        "depth_err": depth_err, "width_err": width_err,
        "ar_err": ar_err,
        "depth_converged": d_conv, "width_converged": w_conv,
        "depth_n_stable": d_n, "width_n_stable": w_n,
        "converged": d_conv and w_conv,
        "case_error": case_error,
        "series": series,
    }


def should_abort(case_dir: Path, exp_depth_um: float, exp_width_um: float,
                 error_threshold: float, stability_tol: float,
                 min_points: int, end_time: float = 0.0,
                 min_progress: float = 0.70) -> tuple[bool, str]:
    """Abort only when >= min_progress of sim is done, tail has STABILIZED, and error > threshold.

    Uses tail stability (last 25% of current data) so partial-data polling
    doesn't falsely trigger on a tiny mid-60% window. Final scoring uses
    _stable() on the complete dataset after endTime.
    """
    series = read_series(case_dir)
    n = len(series)
    if n < max(min_points, 2):
        return False, ""
    if end_time > 0 and series[-1]["time"] / end_time < min_progress:
        return False, ""
    depths = [r["depth"] for r in series]
    widths = [r["width"] for r in series]
    depth, d_conv, _ = _stable_tail(depths, stability_tol)
    width, w_conv, _ = _stable_tail(widths, stability_tol)
    if not (d_conv and w_conv):
        return False, ""

    def rel(sim, exp):
        return float("nan") if not exp else (sim - exp) / exp

    depth_err = rel(depth, exp_depth_um)
    width_err = rel(width, exp_width_um)
    exp_ar = exp_depth_um / exp_width_um if exp_width_um else float("nan")
    sim_ar = depth / width if width else float("nan")
    ar_err = (abs(sim_ar - exp_ar) / exp_ar
              if exp_ar and sim_ar == sim_ar and exp_ar == exp_ar
              else float("nan"))

    worst = max((abs(e) for e in (depth_err, width_err, ar_err) if e == e), default=0.0)
    if worst > error_threshold:
        return True, (f"tail stabilized but error {worst*100:.0f}% > "
                      f"{error_threshold*100:.0f}% "
                      f"(d={depth:.1f} w={width:.1f} ar={ar_err:.3f})")
    return False, ""
