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


def _stable(values: list[float], tol: float) -> tuple[float, bool, int]:
    """Find longest stable tail window.

    Walk backward from the end; extend while (max-min)/mean < tol.
    Stop when adding an older point breaks that bound — that's the transient edge.
    Returns (mean_of_stable_window, converged, n_stable).
    converged = True only when n_stable >= 2.
    """
    n = len(values)
    if n == 0:
        return float("nan"), False, 0
    if n == 1:
        return values[0], False, 1

    stable_n = 1
    for size in range(2, n + 1):
        w = values[-size:]
        mean_v = sum(w) / len(w)
        if mean_v == 0:
            break
        if (max(w) - min(w)) / abs(mean_v) < tol:
            stable_n = size
        else:
            break

    w = values[-stable_n:]
    return sum(w) / len(w), stable_n >= 2, stable_n


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


def compute_physics_penalty(params: dict, config: dict | None = None) -> float:
    """Wiedemann-Franz correlation penalty.

    For metals, electronic thermal conductivity scales as κ_e ∝ 1/ρ_elec:
    the same electron-scattering mechanisms that increase electrical resistivity
    MUST decrease electronic thermal conductivity.

    This penalty fires when the optimizer picks BOTH high ρ_elec (reducing
    absorption → shallower pool) AND high kappa_liquid (increasing conduction →
    deeper pool) — a physically contradictory compensation that might fit the
    calibration data but will fail under different process conditions.

    Reference point: nominal 316L values (ρ_ref=7e-7 Ω·m, κ_ref=26.9 W/m/K).
    Allowed κ at a given ρ = κ_ref × (ρ_ref / ρ_elec) + κ_lattice_max.
    """
    if config and not config.get("enabled", True):
        return 0.0
    weight = config.get("weight", 0.05) if config else 0.05

    RHO_REF = 7e-7         # Ω·m, nominal 316L at Tmelt
    K_REF = 26.9            # W/m/K, nominal total κ at Tliquidus
    K_LATTICE_MAX = 10.0    # W/m/K, generous upper bound for phonon contribution

    rho_elec = params.get("elec_resistivity", RHO_REF)
    kappa_liq = params.get("kappa_liquid_value", K_REF)

    # Expected electronic κ from WF scaling about the reference point
    k_e_expected = K_REF * (RHO_REF / max(rho_elec, 1e-12))
    # Maximum physically plausible total κ at this ρ_elec
    k_max_allowed = k_e_expected + K_LATTICE_MAX

    if kappa_liq > k_max_allowed:
        return weight * (kappa_liq - k_max_allowed) / k_max_allowed
    return 0.0


def should_abort(case_dir: Path, exp_depth_um: float, exp_width_um: float,
                 error_threshold: float, stability_tol: float,
                 min_points: int) -> tuple[bool, str]:
    """Abort only when the pool has STABILIZED but is still > threshold off."""
    ev = evaluate(case_dir, exp_depth_um, exp_width_um, stability_tol)
    if ev["n_points"] < max(min_points, 2):
        return False, ""
    if not ev["converged"]:
        return False, ""
    worst = max((abs(e) for e in (ev["depth_err"], ev["width_err"], ev.get("ar_err", float("nan")))
                if e == e), default=0.0)
    if worst > error_threshold:
        return True, (f"stabilized but error {worst*100:.0f}% > "
                      f"{error_threshold*100:.0f}% "
                      f"(d={ev['sim_depth_um']:.1f} w={ev['sim_width_um']:.1f} "
                      f"ar={ev.get('ar_err', float('nan')):.3f})")
    return False, ""
