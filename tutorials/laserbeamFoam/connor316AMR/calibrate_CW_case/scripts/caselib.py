"""Case construction for CW melt-pool calibration.

Builds a runnable laserbeamFoam case from `template_case/` by applying:
  - candidate (thermophysical) params  -> shared across all cases of a candidate
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
    v = case["v_scan_mm_s"] / 1000.0          # m/s
    if "laser_start_z_m" in case and "laser_end_z_m" in case:
        z0 = float(case["laser_start_z_m"])
        z1 = float(case["laser_end_z_m"])
        scan = z1 - z0
        laser_end_time = scan / v if v > 0 else float(case.get("end_time_s", control["endTime"]))
        laser_off_time = float(case.get("laser_off_time_s", laser_end_time))
        end_time = float(case.get("end_time_s", control.get("endTime", laser_off_time)))
    elif "z_track_m" in geom:
        scan = geom["z_track_m"]
        end_time = scan / v
        laser_end_time = end_time
        laser_off_time = end_time
        z0 = lead_in
        z1 = lead_in + scan
    else:
        end_time = control["endTime"]
        scan = v * end_time                    # distance laser travels in window
        laser_end_time = end_time
        laser_off_time = end_time
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
        "laser_end_time": laser_end_time, "laser_off_time": laser_off_time,
    }


# --------------------------------------------------------------------------- #
# Dictionary patchers
# --------------------------------------------------------------------------- #
def _sub_entry(text: str, key: str, value: str) -> str:
    new, k = re.subn(rf"(\b{key}\s+)[^;]+;", rf"\g<1>{value};", text, count=1)
    if k == 0:
        raise ValueError(f"entry '{key}' not found")
    return new


def _sub_entry_in_block(text: str, block: str, key: str, value: str) -> str:
    pattern = rf"({block}\s*\{{)(.*?)(\n\}})"
    m = re.search(pattern, text, flags=re.S)
    if not m:
        raise ValueError(f"block '{block}' not found")
    body = _sub_entry(m.group(2), key, value)
    return text[:m.start()] + m.group(1) + body + m.group(3) + text[m.end():]


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
        f"    ({g(geo['laser_end_time'])}  ({xc}  0  {g(geo['z1'])}))\n"
        ")\n"
    )


def patch_timeVsLaserPower(path: Path, power_w: float, laser_off_time: float,
                           end_time: float) -> None:
    path.write_text(
        "(\n"
        f"    (0              {g(power_w)})\n"
        f"    ({g(laser_off_time)}     {g(power_w)})\n"
        f"    ({g(laser_off_time + 1e-6)}  0)\n"
        f"    ({g(max(end_time, laser_off_time + 5e-6))}  0)\n"
        ")\n"
    )


def patch_trackProperties(path: Path, end_time: float) -> None:
    path.write_text(_sub_entry(path.read_text(), "trackDuration", g(end_time)))


def patch_LaserProperties(path: Path, radius_m: float) -> None:
    path.write_text(_sub_entry(path.read_text(), "laserRadius", g(radius_m)))


def patch_dynamicMeshDict(path: Path, max_refinement: int) -> None:
    t = path.read_text()
    t = re.sub(
        r"(maxRefinement\s+)\d+;.*",
        rf"\g<1>{max_refinement};",
        t,
        count=1,
    )
    path.write_text(t)


def patch_decomposeParDict(path: Path, cores: int) -> None:
    t, n = re.subn(r"(numberOfSubdomains\s+)\d+;", rf"\g<1>{cores};",
                   path.read_text(), count=1)
    if n == 0:
        raise ValueError(f"entry 'numberOfSubdomains' not found in {path}")
    path.write_text(t)


def _scale_table(text: str, name: str, factor: float) -> str:
    """Scale every value in an OpenFOAM two-column property table."""
    m = re.search(rf"({name}[^\n]*\n\s*\()(.*?)(\n\s*\)\s*;)", text, re.S)
    if not m:
        raise ValueError(f"table '{name}' not found")
    n_rows = 0

    def _row(mm: re.Match) -> str:
        nonlocal n_rows
        val = float(mm.group(2))
        val *= factor
        n_rows += 1
        return f"({mm.group(1)}    {g(val)})"

    body = re.sub(r"\(\s*([0-9.eE+-]+)\s+([0-9.eE+-]+)\s*\)", _row, m.group(2))
    if n_rows == 0:
        raise ValueError(f"table '{name}': no rows found")
    return text[:m.start()] + m.group(1) + body + m.group(3) + text[m.end():]


def _metal_phase_temperature(text: str, name: str) -> float:
    m = re.search(rf"(metal\s*\{{.*?\b{name}\s+)([0-9.eE+-]+)\s*;", text, re.S)
    if not m:
        raise ValueError(f"metal '{name}' not found")
    return float(m.group(2))


def _scale_table_by_phase(text: str, name: str, solid_factor: float,
                          liquid_factor: float, t_solidus: float,
                          t_liquidus: float) -> str:
    """Scale a table with solid/liquid factors and linear mushy-zone blending."""
    if t_liquidus <= t_solidus:
        raise ValueError("Tliquidus must be greater than Tsolidus")
    m = re.search(rf"({name}[^\n]*\n\s*\()(.*?)(\n\s*\)\s*;)", text, re.S)
    if not m:
        raise ValueError(f"table '{name}' not found")
    n_rows = 0

    def _factor(temp: float) -> float:
        if temp <= t_solidus:
            return solid_factor
        if temp >= t_liquidus:
            return liquid_factor
        w = (temp - t_solidus) / (t_liquidus - t_solidus)
        return solid_factor * (1.0 - w) + liquid_factor * w

    def _row(mm: re.Match) -> str:
        nonlocal n_rows
        temp = float(mm.group(1))
        val = float(mm.group(2)) * _factor(temp)
        n_rows += 1
        return f"({mm.group(1)}    {g(val)})"

    body = re.sub(r"\(\s*([0-9.eE+-]+)\s+([0-9.eE+-]+)\s*\)", _row, m.group(2))
    if n_rows == 0:
        raise ValueError(f"table '{name}': no rows found")
    return text[:m.start()] + m.group(1) + body + m.group(3) + text[m.end():]


def patch_transportProperties(path: Path, params: dict,
                               v_scan_mm_s: float = 900.0) -> None:
    """Patch candidate thermophysical parameters into transportProperties.

    The cp/kappa tables can be scaled either as whole tables using the legacy
    table_*_scale parameters, or split into solid/liquid factors with linear
    blending through the mushy interval. V_scan is patched per case.
    """
    t = path.read_text()

    # -- scalars -----------------------------------------------------------
    t = _sub_entry(t, "elec_resistivity", g(params["elec_resistivity"]))
    t = _sub_entry(t, "sigma", g(params["sigma"]))
    t = _sub_entry(t, "Marangoni_Constant", g(params["Marangoni_Constant"]))
    t = _sub_entry(t, "LatentHeatVap", g(params["LatentHeatVap"]))
    t = _sub_entry_in_block(t, "metal", "rho", g(params["rho"]))
    t = _sub_entry_in_block(t, "metal", "nu", g(params["nu"]))

    # -- cp / kappa tables --------------------------------------------------
    t_solidus = _metal_phase_temperature(t, "Tsolidus")
    t_liquidus = _metal_phase_temperature(t, "Tliquidus")

    if "table_kappa_solid_scale" in params or "table_kappa_liquid_scale" in params:
        k_solid = params.get("table_kappa_solid_scale", params.get("table_kappa_scale", 1.0))
        k_liquid = params.get("table_kappa_liquid_scale", params.get("table_kappa_scale", 1.0))
        t = _scale_table_by_phase(t, "table_kappa", k_solid, k_liquid, t_solidus, t_liquidus)
    else:
        t = _scale_table(t, "table_kappa", params["table_kappa_scale"])

    if "table_cp_solid_scale" in params or "table_cp_liquid_scale" in params:
        cp_solid = params.get("table_cp_solid_scale", params.get("table_cp_scale", 1.0))
        cp_liquid = params.get("table_cp_liquid_scale", params.get("table_cp_scale", 1.0))
        t = _scale_table_by_phase(t, "table_cp", cp_solid, cp_liquid, t_solidus, t_liquidus)
    else:
        t = _scale_table(t, "table_cp", params["table_cp_scale"])

    # -- case-specific scan speed -------------------------------------------
    t = _sub_entry(t, "V_scan", g(v_scan_mm_s / 1000.0))

    path.write_text(t)


def patch_analyzer_defaults(path: Path) -> None:
    """Point copied analyzer at the per-candidate experimental summary."""
    t = path.read_text()
    t = re.sub(
        r'^EXP_SUMMARY_CSV = Path\(__file__\)\.resolve\(\)\.parents\[1\] / ".*"$',
        'EXP_SUMMARY_CSV = Path(__file__).resolve().parents[1] / "exp_summary.csv"',
        t,
        count=1,
        flags=re.M,
    )
    t = re.sub(
        r'^EXP_MASK_IMAGE = Path\(__file__\)\.resolve\(\)\.parents\[1\] / ".*"$',
        'EXP_MASK_IMAGE = Path(__file__).resolve().parents[1] / "exp_mask.png"',
        t,
        count=1,
        flags=re.M,
    )
    path.write_text(t)


# --------------------------------------------------------------------------- #
# Top-level builder
# --------------------------------------------------------------------------- #
def build_case(template_dir: Path, dest: Path, candidate_params: dict,
               case_cfg: dict, geom: dict, control: dict, cores: int,
               max_refinement: int = 2) -> dict:
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
    ctrl = dict(control)
    ctrl["endTime"] = geo["end_time"]

    patch_blockMeshDict(dest / "system" / "blockMeshDict", geo)
    patch_setFieldsDict(dest / "system" / "setFieldsDict", geo)
    patch_controlDict(dest / "system" / "controlDict", ctrl)
    patch_dynamicMeshDict(dest / "constant" / "dynamicMeshDict", max_refinement)
    patch_decomposeParDict(dest / "system" / "decomposeParDict", cores)
    P_sim = case_cfg["P_laser_W"] * candidate_params.get("power_scale", 1.0)
    r_sim = case_cfg["laserRadius_m"] * candidate_params.get("radius_scale", 1.0)

    patch_timeVsLaserPosition(dest / "constant" / "timeVsLaserPosition", geo)
    patch_timeVsLaserPower(
        dest / "constant" / "timeVsLaserPower",
        P_sim,
        geo["laser_off_time"],
        geo["end_time"],
    )
    patch_trackProperties(dest / "constant" / "trackProperties", geo["laser_end_time"])
    patch_LaserProperties(dest / "constant" / "LaserProperties", r_sim)
    patch_transportProperties(dest / "constant" / "transportProperties",
                             candidate_params, case_cfg["v_scan_mm_s"])
    patch_analyzer_defaults(dest / "scripts" / "analyze_meltpool_vtu.py")

    surface_y_um = float(geom.get("surface_y_um", geo["gas_thickness"] * 1e6))
    record = {"case": case_cfg, "geometry": geo, "params": candidate_params,
              "cores": cores, "surface_y_um": surface_y_um,
              "max_refinement": max_refinement}
    (dest / "case_build.json").write_text(json.dumps(record, indent=2))
    (dest / "exp_summary.csv").write_text(
        "width_um,depth_um\n"
        f"{case_cfg['exp_width_um']},{case_cfg['exp_depth_um']}\n"
    )
    case_info = {
        "name": case_cfg["name"],
        "material": "SS316L",
        "power_W": case_cfg["P_laser_W"],
        "scan_speed_mm_s": case_cfg["v_scan_mm_s"],
        "laser_radius_m": case_cfg["laserRadius_m"],
        "laser_start_z_m": geo["z0"],
        "laser_end_z_m": geo["z1"],
        "laser_off_time_s": geo["laser_off_time"],
        "end_time_s": geo["end_time"],
        "target_width_um": case_cfg["exp_width_um"],
        "target_depth_um": case_cfg["exp_depth_um"],
        "base_mesh_size_um": geom["base_cell"] * 1e6,
        "amr_levels": max_refinement,
        "surface_y_um": surface_y_um,
    }
    (dest / "case_info.json").write_text(json.dumps(case_info, indent=2) + "\n")
    return record
