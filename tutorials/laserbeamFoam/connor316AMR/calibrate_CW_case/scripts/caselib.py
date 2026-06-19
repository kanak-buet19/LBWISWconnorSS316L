"""Case construction for CW melt-pool calibration.

Builds a runnable laserbeamFoam case from `template_case/` by applying:
  - candidate (thermophysical) params  -> shared across BOTH cases of a candidate
  - case (process + geometry) params   -> P, v, laserRadius, domain size, scan path

All edits are plain-text regex patches on the OpenFOAM dictionaries so no
OpenFOAM/foamDictionary call is needed at build time.
"""

from __future__ import annotations

import json
import math
import re
import shutil
from pathlib import Path


def g(x: float) -> str:
    """Format a number for an OpenFOAM dictionary."""
    return f"{x:.8g}"


def n_cells(length: float, base: float) -> int:
    return max(1, int(math.ceil(round(length / base, 6))))


def round_up_to_base(length: float, base: float) -> float:
    return n_cells(length, base) * base


# --------------------------------------------------------------------------- #
# Geometry
# --------------------------------------------------------------------------- #
def build_geometry(case: dict, geom: dict, control: dict) -> dict:
    base = geom["base_cell"]
    gas = geom["gas_thickness"]
    x_extent = geom["x_extent"]
    lead_in = geom["lead_in"]
    tail = geom["tail"]
    end_time = control["endTime"]

    v = case["v_scan_mm_s"] / 1000.0          # m/s
    scan = v * end_time                        # distance laser travels in window
    z0 = lead_in
    z1 = lead_in + scan
    z_extent = round_up_to_base(z1 + tail, base)

    substrate = max(geom["substrate_depth_factor"] * case["exp_depth_um"] * 1e-6,
                    geom["substrate_min"])
    substrate = round_up_to_base(substrate, base)
    y_extent = gas + substrate
    x_extent = round_up_to_base(x_extent, base)

    return {
        "x_extent": x_extent, "y_extent": y_extent, "z_extent": z_extent,
        "nx": n_cells(x_extent, base), "ny": n_cells(y_extent, base),
        "nz": n_cells(z_extent, base),
        "gas_thickness": gas, "x_center": x_extent / 2.0,
        "z0": z0, "z1": z1, "scan": scan, "end_time": end_time,
    }


# --------------------------------------------------------------------------- #
# Dictionary patchers
# --------------------------------------------------------------------------- #
def _sub_entry(text: str, key: str, value: str) -> str:
    new, k = re.subn(rf"(\b{key}\s+)[^;]+;", rf"\g<1>{value};", text, count=1)
    if k == 0:
        raise ValueError(f"entry '{key}' not found")
    return new


def patch_blockMeshDict(path: Path, geo: dict) -> None:
    t = path.read_text()
    X, Y, Z = g(geo["x_extent"]), g(geo["y_extent"]), g(geo["z_extent"])
    verts = (
        "vertices\n(\n"
        f"    (0 0 0)\n"
        f"    ({X} 0 0)\n"
        f"    ({X} {Y} 0)\n"
        f"    (0 {Y} 0)\n"
        f"    (0 0 {Z})\n"
        f"    ({X} 0 {Z})\n"
        f"    ({X} {Y} {Z})\n"
        f"    (0 {Y} {Z})\n"
        ");"
    )
    # consume the WHOLE vertices block up to its terminating ');'
    t = re.sub(r"vertices\s*\(.*?\n\);", verts, t, count=1, flags=re.S)
    t = re.sub(r"hex \(0 1 2 3 4 5 6 7\) \([^)]*\)",
               f"hex (0 1 2 3 4 5 6 7) ({geo['nx']} {geo['ny']} {geo['nz']})", t, count=1)
    path.write_text(t)


def patch_setFieldsDict(path: Path, geo: dict) -> None:
    t = path.read_text()
    box = (f"box (0 {g(geo['gas_thickness'])} 0) "
           f"({g(geo['x_extent'])} {g(geo['y_extent'])} {g(geo['z_extent'])});")
    t = re.sub(r"box \([^)]*\) \([^)]*\);", box, t, count=1)
    path.write_text(t)


def patch_controlDict(path: Path, control: dict) -> None:
    t = path.read_text()
    t = _sub_entry(t, "endTime", g(control["endTime"]))
    t = _sub_entry(t, "writeInterval", g(control["writeInterval"]))
    t = _sub_entry(t, "deltaT", g(control["deltaT"]))
    t = _sub_entry(t, "maxDeltaT", g(control["maxDeltaT"]))
    t = _sub_entry(t, "maxCo", g(control["maxCo"]))
    t = _sub_entry(t, "maxAlphaCo", g(control["maxAlphaCo"]))
    path.write_text(t)


def patch_timeVsLaserPosition(path: Path, geo: dict) -> None:
    xc = g(geo["x_center"])
    path.write_text(
        "(\n"
        f"    (0           ({xc}  0  {g(geo['z0'])}))\n"
        f"    ({g(geo['end_time'])}  ({xc}  0  {g(geo['z1'])}))\n"
        ")\n"
    )


def patch_timeVsLaserPower(path: Path, power_w: float, end_time: float) -> None:
    path.write_text(
        "(\n"
        f"    (0              {g(power_w)})\n"
        f"    ({g(end_time)}     {g(power_w)})\n"
        f"    ({g(end_time + 1e-6)}  0)\n"
        f"    ({g(end_time + 1e-5)}  0)\n"
        ")\n"
    )


def patch_trackProperties(path: Path, end_time: float) -> None:
    path.write_text(_sub_entry(path.read_text(), "trackDuration", g(end_time)))


def patch_LaserProperties(path: Path, radius_m: float) -> None:
    path.write_text(_sub_entry(path.read_text(), "laserRadius", g(radius_m)))


def patch_decomposeParDict(path: Path, cores: int) -> None:
    t, n = re.subn(r"(numberOfSubdomains\s+)\d+;", rf"\g<1>{cores};",
                   path.read_text(), count=1)
    if n == 0:
        raise ValueError(f"entry 'numberOfSubdomains' not found in {path}")
    path.write_text(t)


T_SOLIDUS = 1658.0
T_LIQUIDUS = 1723.0


def _scale_solid_table(text: str, name: str, factor: float) -> str:
    """Scale solid-phase table entries (T <= Tsolidus) by factor.

    Liquid entries (T > Tsolidus) are left untouched — they will be set
    directly by _set_liquid_table afterwards.
    """
    m = re.search(rf"({name}[^\n]*\n\s*\()(.*?)(\n\s*\)\s*;)", text, re.S)
    if not m:
        raise ValueError(f"table '{name}' not found")
    n_rows = 0

    def _row(mm: re.Match) -> str:
        nonlocal n_rows
        T = float(mm.group(1))
        val = float(mm.group(2))
        if T <= T_SOLIDUS:
            val *= factor
            n_rows += 1
        return f"({mm.group(1)}    {g(val)})"

    body = re.sub(r"\(\s*([0-9.eE+-]+)\s+([0-9.eE+-]+)\s*\)", _row, m.group(2))
    if n_rows == 0:
        raise ValueError(f"table '{name}': no solid rows found (T <= {T_SOLIDUS})")
    return text[:m.start()] + m.group(1) + body + m.group(3) + text[m.end():]


def _set_liquid_table(text: str, name: str, value: float,
                       slope: float = 0.0, T_cap: float = 3122.0) -> str:
    """Set liquid-phase table entries (T >= Tliquidus) to value + slope*(T-Tliquidus).

    Values above T_cap (Tvap_atm) are pinned to the T_cap delta so the
    far-field extrapolation stays bounded.
    """
    m = re.search(rf"({name}[^\n]*\n\s*\()(.*?)(\n\s*\)\s*;)", text, re.S)
    if not m:
        raise ValueError(f"table '{name}' not found for liquid-value patch")
    n_rows = 0

    def _row(mm: re.Match) -> str:
        nonlocal n_rows
        T = float(mm.group(1))
        if T >= T_LIQUIDUS:
            n_rows += 1
            T_eff = min(T, T_cap)
            val = value + slope * (T_eff - T_LIQUIDUS)
            return f"({mm.group(1)}    {g(val)})"
        return mm.group(0)

    body = re.sub(r"\(\s*([0-9.eE+-]+)\s+([0-9.eE+-]+)\s*\)", _row, m.group(2))
    if n_rows == 0:
        raise ValueError(f"table '{name}': no liquid rows found (T >= {T_LIQUIDUS})")
    return text[:m.start()] + m.group(1) + body + m.group(3) + text[m.end():]


def patch_transportProperties(path: Path, params: dict,
                               v_scan_mm_s: float = 900.0) -> None:
    """Patch candidate thermophysical parameters into transportProperties.

    Marangoni_Constant is derived from sigma × dSigmadT_norm (physically
    coupled — prevents optimizer from picking unrealistic dσ/dT ratios).
    Solid cp/kappa are scaled by narrow solid_scale factors (well-known from
    DSC). Liquid cp/kappa are set directly at Tliquidus plus a narrow slope.
    V_scan is patched per-case.
    """
    t = path.read_text()

    # -- scalars -----------------------------------------------------------
    t = _sub_entry(t, "elec_resistivity", g(params["elec_resistivity"]))
    t = _sub_entry(t, "beta_r", g(params["beta_r"]))
    t = _sub_entry(t, "sigma", g(params["sigma"]))
    marangoni = params["sigma"] * params["dSigmadT_norm"]
    t = _sub_entry(t, "Marangoni_Constant", g(marangoni))
    t = _sub_entry(t, "LatentHeatVap", g(params["LatentHeatVap"]))
    t = _sub_entry(t, "LeeCoeff", g(params["LeeCoeff"]))
    # metal rho only (gas rho = 1, no collision with baseline 7950)
    t = re.sub(r"(\brho\s+)7950(\.\d+)?\b", rf"\g<1>{g(params['rho'])}", t, count=1)
    # metal nu first (gas nu is 1.48e-05, comes after metal block)
    t = re.sub(r"(\bnu\s+)7e-7\b", rf"\g<1>{g(params['nu'])}", t, count=1)
    # LatentHeat (fusion): requires trailing whitespace to avoid matching LatentHeatVap
    t = _sub_entry(t, "LatentHeat", g(params["LatentHeat"]))

    # -- cp / kappa tables --------------------------------------------------
    # 1. Scale SOLID entries only (narrow range — well-known from DSC)
    t = _scale_solid_table(t, "table_cp",    params["cp_solid_scale"])
    t = _scale_solid_table(t, "table_kappa", params["kappa_solid_scale"])
    # 2. Set LIQUID entries directly (value at Tliquidus + narrow slope)
    t = _set_liquid_table(t, "table_kappa", params["kappa_liquid_value"],
                           params["kappa_liquid_slope"])
    t = _set_liquid_table(t, "table_cp",    params["cp_liquid_value"],
                           params["cp_liquid_slope"])

    # -- case-specific scan speed -------------------------------------------
    t = _sub_entry(t, "V_scan", g(v_scan_mm_s / 1000.0))

    path.write_text(t)


# --------------------------------------------------------------------------- #
# Top-level builder
# --------------------------------------------------------------------------- #
def build_case(template_dir: Path, dest: Path, candidate_params: dict,
               case_cfg: dict, geom: dict, control: dict, cores: int) -> dict:
    if dest.exists():
        shutil.rmtree(dest)
    shutil.copytree(template_dir, dest)
    # strip any stray run artifacts copied from template
    for junk in ("0", "VTK", "VTKs", "post-processing-data", "absorptivity_vs_time",
                 ".mplconfig", "dynamicCode"):
        p = dest / junk
        if p.is_dir():
            shutil.rmtree(p)
    for p in dest.glob("processor*"):
        shutil.rmtree(p)
    for p in dest.glob("log.*"):
        p.unlink()

    geo = build_geometry(case_cfg, geom, control)

    patch_blockMeshDict(dest / "system" / "blockMeshDict", geo)
    patch_setFieldsDict(dest / "system" / "setFieldsDict", geo)
    patch_controlDict(dest / "system" / "controlDict", control)
    patch_decomposeParDict(dest / "system" / "decomposeParDict", cores)
    patch_timeVsLaserPosition(dest / "constant" / "timeVsLaserPosition", geo)
    patch_timeVsLaserPower(dest / "constant" / "timeVsLaserPower",
                           case_cfg["P_laser_W"], geo["end_time"])
    patch_trackProperties(dest / "constant" / "trackProperties", geo["end_time"])
    patch_LaserProperties(dest / "constant" / "LaserProperties", case_cfg["laserRadius_m"])
    patch_transportProperties(dest / "constant" / "transportProperties",
                             candidate_params, case_cfg["v_scan_mm_s"])

    record = {"case": case_cfg, "geometry": geo, "params": candidate_params,
              "cores": cores, "surface_y_um": geom["surface_y_um"]}
    (dest / "case_build.json").write_text(json.dumps(record, indent=2))
    (dest / "exp_summary.csv").write_text(
        "width_um,depth_um\n"
        f"{case_cfg['exp_width_um']},{case_cfg['exp_depth_um']}\n"
    )
    return record
