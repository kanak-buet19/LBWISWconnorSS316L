#!/usr/bin/env python3
"""Prepare Hofmann validation cases from cases.json."""

from __future__ import annotations

import argparse
import json
import math
import re
import shutil
from copy import deepcopy
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
CONFIG = ROOT / "cases.json"
MATERIALS_DIR = ROOT / "materials"


def read_text(path: Path) -> str:
    return path.read_text(encoding="utf-8")


def write_text(path: Path, text: str) -> None:
    path.write_text(text, encoding="utf-8")


def replace_regex(path: Path, pattern: str, repl: str) -> None:
    text = read_text(path)
    new_text, count = re.subn(pattern, repl, text, flags=re.MULTILINE)
    if count == 0:
        raise RuntimeError(f"pattern not found in {path}: {pattern}")
    write_text(path, new_text)


def fmt_m(value: float) -> str:
    return f"{value:.6e}"


def fmt_param(value: float) -> str:
    text = f"{value:.6g}"
    return text.replace("+", "")


def resistivity_tag(value: float) -> str:
    return f"rho{fmt_param(value).replace('-', 'm').replace('.', 'p')}"


def parse_resistivity_values(case: dict) -> list[float]:
    keys = (
        "electric_resistivity",
        "electric_resistivity_sweep",
        "electric_resistivity_values",
        "elec_resistivity",
        "elec_resistivity_sweep",
        "elec_resistivity_values",
    )
    present = [key for key in keys if key in case]
    if not present:
        return []
    if len(present) > 1:
        raise RuntimeError(
            f"case {case['name']}: set only one electric resistivity key, found {', '.join(present)}"
        )

    raw = case[present[0]]
    if isinstance(raw, str):
        parts = [part.strip() for part in raw.split(",")]
        values = [float(part) for part in parts if part]
    elif isinstance(raw, (int, float)):
        values = [float(raw)]
    elif isinstance(raw, list):
        values = [float(value) for value in raw]
    else:
        raise RuntimeError(
            f"case {case['name']}: electric resistivity must be a number, list, or comma-separated string"
        )

    if not values:
        raise RuntimeError(f"case {case['name']}: electric resistivity sweep is empty")
    for value in values:
        if value <= 0:
            raise RuntimeError(f"case {case['name']}: electric resistivity must be positive")
    return values


def set_electric_resistivity(case_dir: Path, value: float) -> None:
    replace_regex(
        case_dir / "constant" / "transportProperties",
        r"^(\s*elec_resistivity\s+)[-+0-9.eE]+(\s*;.*)$",
        rf"\g<1>{fmt_param(value)}\2",
    )


def expand_case_sweeps(cases: list[dict]) -> list[dict]:
    expanded = []
    for case in cases:
        values = parse_resistivity_values(case)
        if len(values) <= 1:
            item = deepcopy(case)
            if values:
                item["_electric_resistivity_value"] = values[0]
                item["electric_resistivity_value"] = values[0]
            expanded.append(item)
            continue

        for value in values:
            item = deepcopy(case)
            item["base_name"] = case["name"]
            item["name"] = f"{case['name']}_{resistivity_tag(value)}"
            item["_electric_resistivity_value"] = value
            item["electric_resistivity_value"] = value
            expanded.append(item)
    return expanded


def target_tag(case: dict) -> str:
    power = int(case["power_W"])
    speed = int(case["scan_speed_mm_s"])
    prefix = case.get("exp_tag", "hofmann")
    return f"exp_{prefix}_{power}W_{speed}_summary.csv"


def mesh_bounds(block_mesh: Path) -> tuple[float, float, float]:
    text = read_text(block_mesh)
    match = re.search(r"vertices\s*\((.*?)\);\s*blocks", text, flags=re.DOTALL)
    if not match:
        raise RuntimeError(f"could not parse vertices from {block_mesh}")
    vertices = re.findall(
        r"\(\s*([-+0-9.eE]+)\s+([-+0-9.eE]+)\s+([-+0-9.eE]+)\s*\)",
        match.group(1),
    )
    if not vertices:
        raise RuntimeError(f"no vertices found in {block_mesh}")
    xs, ys, zs = zip(*((float(x), float(y), float(z)) for x, y, z in vertices))
    return max(xs), max(ys), max(zs)


def mesh_cells(block_mesh: Path) -> tuple[int, int, int]:
    text = read_text(block_mesh)
    match = re.search(r"hex\s*\([^)]+\)\s*\(\s*(\d+)\s+(\d+)\s+(\d+)\s*\)", text)
    if not match:
        raise RuntimeError(f"could not parse block cells from {block_mesh}")
    return tuple(int(v) for v in match.groups())


_GAS_M = 0.128e-3
_SUBSTRATE_FACTOR = 1.5
_SUBSTRATE_MIN_M = 224e-6


def resize_domain_depth(case_dir: Path, case: dict) -> None:
    """Set Y-extent: explicit domain_depth_m overrides auto 1.5×target_depth_um."""
    base_m = float(case.get("base_mesh_size_um", 32.0)) * 1e-6
    if "domain_depth_m" in case:
        ny = max(1, round(float(case["domain_depth_m"]) / base_m))
        new_y_max = ny * base_m
    elif case.get("target_depth_um") is not None:
        depth_um = float(case["target_depth_um"])
        n = max(1, math.ceil(round(_SUBSTRATE_FACTOR * depth_um * 1e-6 / base_m, 6)))
        n_min = max(1, math.ceil(round(_SUBSTRATE_MIN_M / base_m, 6)))
        substrate_m = max(n, n_min) * base_m
        new_y_max = _GAS_M + substrate_m
        ny = max(1, round(new_y_max / base_m))
    else:
        return

    block_mesh = case_dir / "system" / "blockMeshDict"
    x_max, old_y_max, z_max = mesh_bounds(block_mesh)
    if abs(old_y_max - new_y_max) < 1e-10:
        return

    X, Y, Z = fmt_m(x_max), fmt_m(new_y_max), fmt_m(z_max)
    new_verts = (
        "vertices\n(\n"
        f"    (0      0    0    )  // 0\n"
        f"    ({X}  0    0    )  // 1\n"
        f"    ({X}  {Y}  0    )  // 2\n"
        f"    (0      {Y}  0    )  // 3\n"
        f"    (0      0    {Z}  )  // 4\n"
        f"    ({X}  0    {Z}  )  // 5\n"
        f"    ({X}  {Y}  {Z}  )  // 6\n"
        f"    (0      {Y}  {Z}  )  // 7\n"
        ");"
    )
    text = re.sub(r"vertices\s*\(.*?\);", new_verts, read_text(block_mesh), count=1, flags=re.DOTALL)
    text = re.sub(
        r"(hex\s*\([^)]+\)\s*\(\s*)(\d+)(\s+)(\d+)(\s+)(\d+)(\s*\))",
        lambda m: f"{m.group(1)}{m.group(2)}{m.group(3)}{ny}{m.group(5)}{m.group(6)}{m.group(7)}",
        text, count=1,
    )
    write_text(block_mesh, text)


def resize_domain_x(case_dir: Path, case: dict) -> None:
    """Set X extent to domain_x_m, rounded to base cell."""
    domain_x_m = case.get("domain_x_m")
    if domain_x_m is None:
        return
    base_m = float(case.get("base_mesh_size_um", 32.0)) * 1e-6
    nx = max(1, round(float(domain_x_m) / base_m))
    new_x_max = nx * base_m

    block_mesh = case_dir / "system" / "blockMeshDict"
    old_x_max, y_max, z_max = mesh_bounds(block_mesh)
    if abs(old_x_max - new_x_max) < 1e-10:
        return

    X, Y, Z = fmt_m(new_x_max), fmt_m(y_max), fmt_m(z_max)
    new_verts = (
        "vertices\n(\n"
        f"    (0      0    0    )  // 0\n"
        f"    ({X}  0    0    )  // 1\n"
        f"    ({X}  {Y}  0    )  // 2\n"
        f"    (0      {Y}  0    )  // 3\n"
        f"    (0      0    {Z}  )  // 4\n"
        f"    ({X}  0    {Z}  )  // 5\n"
        f"    ({X}  {Y}  {Z}  )  // 6\n"
        f"    (0      {Y}  {Z}  )  // 7\n"
        ");"
    )
    text = re.sub(r"vertices\s*\(.*?\);", new_verts, read_text(block_mesh), count=1, flags=re.DOTALL)
    text = re.sub(
        r"(hex\s*\([^)]+\)\s*\(\s*)(\d+)(\s+)(\d+)(\s+)(\d+)(\s*\))",
        lambda m: f"{m.group(1)}{nx}{m.group(3)}{m.group(4)}{m.group(5)}{m.group(6)}{m.group(7)}",
        text, count=1,
    )
    write_text(block_mesh, text)


def resize_domain_z(case_dir: Path, case: dict) -> None:
    """Set Z extent to domain_z_m, rounded to base cell."""
    domain_z_m = case.get("domain_z_m")
    if domain_z_m is None:
        return
    base_m = float(case.get("base_mesh_size_um", 32.0)) * 1e-6
    nz = max(1, round(domain_z_m / base_m))
    new_z_max = nz * base_m

    block_mesh = case_dir / "system" / "blockMeshDict"
    x_max, y_max, old_z_max = mesh_bounds(block_mesh)
    if abs(old_z_max - new_z_max) < 1e-10:
        return

    X, Y, Z = fmt_m(x_max), fmt_m(y_max), fmt_m(new_z_max)
    new_verts = (
        "vertices\n(\n"
        f"    (0      0    0    )  // 0\n"
        f"    ({X}  0    0    )  // 1\n"
        f"    ({X}  {Y}  0    )  // 2\n"
        f"    (0      {Y}  0    )  // 3\n"
        f"    (0      0    {Z}  )  // 4\n"
        f"    ({X}  0    {Z}  )  // 5\n"
        f"    ({X}  {Y}  {Z}  )  // 6\n"
        f"    (0      {Y}  {Z}  )  // 7\n"
        ");"
    )
    text = re.sub(r"vertices\s*\(.*?\);", new_verts, read_text(block_mesh), count=1, flags=re.DOTALL)
    text = re.sub(
        r"(hex\s*\([^)]+\)\s*\(\s*)(\d+)(\s+)(\d+)(\s+)(\d+)(\s*\))",
        lambda m: f"{m.group(1)}{m.group(2)}{m.group(3)}{m.group(4)}{m.group(5)}{nz}{m.group(7)}",
        text, count=1,
    )
    write_text(block_mesh, text)


def configure_mesh(case_dir: Path, case: dict) -> None:
    if "base_mesh_size_um" not in case and "amr_levels" not in case:
        return

    block_mesh = case_dir / "system" / "blockMeshDict"
    dynamic_mesh = case_dir / "constant" / "dynamicMeshDict"

    base_um = float(case.get("base_mesh_size_um", 32.0))
    amr_levels = int(case.get("amr_levels", 2))
    if base_um <= 0:
        raise RuntimeError(f"base_mesh_size_um must be positive for {case['name']}")
    if amr_levels < 0:
        raise RuntimeError(f"amr_levels must be non-negative for {case['name']}")

    finest_um = base_um / (2 ** amr_levels)
    x_max, y_max, z_max = mesh_bounds(block_mesh)
    old_nx, old_ny, old_nz = mesh_cells(block_mesh)

    # Preserve the existing Hofmann template cell counts for its native 32 um base.
    if abs(base_um - 32.0) < 1e-12:
        nx, ny, nz = old_nx, old_ny, old_nz
    else:
        nx = max(1, round(x_max * 1e6 / base_um))
        ny = max(1, round(y_max * 1e6 / base_um))
        nz = max(1, round(z_max * 1e6 / base_um))

    replace_regex(
        block_mesh,
        r"^// Base .* near melt pool\.$",
        f"// Base {base_um:g} um; {amr_levels}-level AMR -> {finest_um:g} um near melt pool.",
    )
    replace_regex(
        block_mesh,
        r"^//   x: .*\bcells\b.*$",
        f"//   x: {x_max * 1e3:.3f} mm  ({nx} cells)",
    )
    replace_regex(
        block_mesh,
        r"^//   y: .*\bcells\b.*$",
        f"//   y: {y_max * 1e3:.3f} mm  ({ny} cells)",
    )
    replace_regex(
        block_mesh,
        r"^//   z: .*\bcells\b.*$",
        f"//   z: {z_max * 1e3:.3f} mm  ({nz} cells)",
    )
    replace_regex(
        block_mesh,
        r"^(\s*hex\s*\([^)]+\)\s*)\(\s*\d+\s+\d+\s+\d+\s*\)(.*)$",
        rf"\1({nx} {ny} {nz})\2",
    )
    replace_regex(
        dynamic_mesh,
        r"^maxRefinement\s+\d+;.*$",
        f"maxRefinement   {amr_levels};  // {base_um:g} um base -> {amr_levels} levels -> {finest_um:g} um",
    )


def configure_case(case_dir: Path, case: dict) -> None:
    power = int(case["power_W"])
    speed = int(case["scan_speed_mm_s"])
    radius = float(case["laser_radius_m"])
    spot = radius * 2e3  # diameter in mm, for comments only
    width_value = case.get("target_width_um")
    depth_value = case.get("target_depth_um")
    has_target = depth_value is not None
    width = float(width_value) if width_value is not None else None
    depth = float(depth_value) if has_target else None
    start_z = float(case.get("laser_start_z_m", 0.100e-3))
    configured_end_z = case.get("laser_end_z_m")
    configured_end_time = case.get("end_time_s")
    configured_laser_off_time = case.get("laser_off_time_s")
    write_interval = float(case.get("write_interval_s", 1e-5))
    if configured_end_z is not None:
        end_z = float(configured_end_z)
        laser_end_time = (end_z - start_z) / (speed / 1000.0)
        end_time = laser_end_time
        case["end_time_s"] = end_time
    elif configured_end_time is not None:
        end_time = float(configured_end_time)
        laser_end_time = end_time
        end_z = start_z + (speed / 1000.0) * end_time
    else:
        raise RuntimeError(f"case {case['name']}: must set either end_time_s or laser_end_z_m")
    laser_off_time = float(configured_laser_off_time) if configured_laser_off_time is not None else laser_end_time
    if laser_off_time > end_time:
        raise RuntimeError(
            f"case {case['name']}: laser_off_time_s ({laser_off_time:g}) cannot exceed end_time_s ({end_time:g})"
        )
    name = case["name"]
    job_name = f"cw{power}W{speed}mms"
    summary_name = target_tag(case) if has_target else f"exp_{name}_summary.csv"
    mask_name = summary_name.replace("_summary.csv", "_mask.png")
    process_comment = f"P={power}W, v={speed}mm/s, d_laser={spot:.2f}mm."
    if has_target:
        if width is not None:
            process_comment += f" Target: depth~{depth:.0f}um, width~{width:.0f}um."
        else:
            process_comment += f" Target: depth~{depth:.0f}um."

    resize_domain_depth(case_dir, case)
    resize_domain_x(case_dir, case)
    resize_domain_z(case_dir, case)
    configure_mesh(case_dir, case)

    job_sh = case_dir / "job.sh"
    if not job_sh.exists():
        template_job = ROOT / case.get("template", "template_case") / "job.sh"
        if template_job.exists():
            shutil.copy2(template_job, job_sh)
    if job_sh.exists():
        replace_regex(job_sh, r"^#SBATCH --job-name=.*$", f"#SBATCH --job-name={job_name}")
        replace_regex(job_sh, r"^# (P=.*|Zhang .*process point:.*)$", f"# {process_comment}")
    replace_regex(case_dir / "system" / "blockMeshDict", r"^// (P=.*|Zhang .*process point:.*)$", f"// {process_comment}")
    replace_regex(
        case_dir / "constant" / "LaserProperties",
        r"^// (P=.*d_laser=.*|N01 bare plate:.*)$",
        f"// P={power}W, v={speed}mm/s, d_laser={spot:.2f}mm (radius={radius * 1e6:.1f}um).",
    )
    replace_regex(
        case_dir / "constant" / "LaserProperties",
        r"^laserRadius\s+[-+0-9.eE]+;.*$",
        f"laserRadius {radius:.6g};  // {spot * 1000.0:.0f} um spot diameter",
    )
    replace_regex(
        case_dir / "constant" / "transportProperties",
        r"^V_scan\s+[-+0-9.eE]+;\s*//.*$",
        f"V_scan          {speed / 1000.0:g};            // Scan speed ({speed} mm/s)",
    )
    replace_regex(
        case_dir / "constant" / "trackProperties",
        r"^// .*validation track.*$|^// 1\.0 mm / .*m/s = .*$",
        f"// {(end_z - start_z) * 1e3:.3f} mm validation track at {speed} mm/s.",
    )
    replace_regex(
        case_dir / "constant" / "trackProperties",
        r"^trackDuration\s+[-+0-9.eE]+;.*$",
        f"trackDuration      {laser_end_time:.4e};",
    )
    replace_regex(
        case_dir / "system" / "controlDict",
        r"^endTime\s+[-+0-9.eE]+;.*$",
        f"endTime         {end_time:.2e};",
    )
    replace_regex(
        case_dir / "system" / "controlDict",
        r"^writeInterval\s+[-+0-9.eE]+;.*$",
        f"writeInterval   {write_interval:.2e};",
    )
    template_dir = ROOT / case.get("template", "template_case")
    shutil.copy2(
        template_dir / "scripts" / "parse_simulation_log.py",
        case_dir / "scripts" / "parse_simulation_log.py",
    )
    replace_regex(
        case_dir / "scripts" / "parse_simulation_log.py",
        r"^    v_scan = [-+0-9.]+$",
        f"    v_scan = {float(speed):.1f}",
    )
    replace_regex(
        case_dir / "scripts" / "analyze_meltpool_vtu.py",
        r'^EXP_SUMMARY_CSV = Path\(__file__\)\.resolve\(\)\.parents\[1\] / ".*"$',
        f'EXP_SUMMARY_CSV = Path(__file__).resolve().parents[1] / "{summary_name}"',
    )
    replace_regex(
        case_dir / "scripts" / "analyze_meltpool_vtu.py",
        r'^EXP_MASK_IMAGE = Path\(__file__\)\.resolve\(\)\.parents\[1\] / ".*"$',
        f'EXP_MASK_IMAGE = Path(__file__).resolve().parents[1] / "{mask_name}"',
    )
    depth_metric = case.get("depth_metric", "meltpool")
    compare_depth_field = "keyholeDepth_um" if depth_metric == "keyhole" else "meltPoolDepth_um"
    replace_regex(
        case_dir / "scripts" / "analyze_meltpool_vtu.py",
        r'^COMPARE_DEPTH_FIELD = ".*"$',
        f'COMPARE_DEPTH_FIELD = "{compare_depth_field}"',
    )

    x_max, y_max, z_max = mesh_bounds(case_dir / "system" / "blockMeshDict")
    replace_regex(
        case_dir / "system" / "setFieldsDict",
        r"^\s*// metal:.*$",
        f"        // metal:     y=0.128mm to {y_max * 1e3:.3f}mm, full z length",
    )
    replace_regex(
        case_dir / "system" / "setFieldsDict",
        r"^\s*box \(0 0\.128e-3 0\) \([^)]+\);$",
        f"        box (0 0.128e-3 0) ({fmt_m(x_max)} {fmt_m(y_max)} {fmt_m(z_max)});",
    )

    laser_x = x_max / 2.0
    write_text(
        case_dir / "constant" / "timeVsLaserPosition",
        "(\n"
        f"    (0           ({fmt_m(laser_x)}  0  {fmt_m(start_z)}))\n"
        f"    ({laser_end_time:.4e}  ({fmt_m(laser_x)}  0  {fmt_m(end_z)}))\n"
        ")\n",
    )
    write_text(
        case_dir / "constant" / "timeVsLaserPower",
        "(\n"
        f"    (0           {power})\n"
        f"    ({laser_off_time:.4e}      {power})\n"
        f"    ({laser_off_time + 1e-6:.4e} 0)\n"
        f"    ({max(end_time, laser_off_time + 5e-6):.4e} 0)\n"
        ")\n",
    )

    for old_summary in case_dir.glob("exp_*_summary.csv"):
        if not has_target or old_summary.name != summary_name:
            old_summary.unlink()
    if has_target:
        if width is not None:
            write_text(case_dir / summary_name, f"width_um,depth_um\n{width},{depth}\n")
        else:
            write_text(case_dir / summary_name, f"depth_um\n{depth}\n")
    write_text(case_dir / "case_info.json", json.dumps(case, indent=2) + "\n")

    print(f"prepared {name}")


def load_config() -> dict:
    return json.loads(CONFIG.read_text(encoding="utf-8"))


def select_cases(cfg: dict, names: list[str], include_all: bool) -> list[dict]:
    use_all = include_all or not cfg.get("run_enabled_only", True)
    cases = cfg["cases"] if use_all else [c for c in cfg["cases"] if c.get("enabled", True)]
    cases = expand_case_sweeps(cases)
    if not names:
        return cases

    by_name = {c["name"]: c for c in cases}
    by_base_name: dict[str, list[dict]] = {}
    for case in cases:
        by_base_name.setdefault(case.get("base_name", case["name"]), []).append(case)
    selected = []
    missing = []
    for name in names:
        case = by_name.get(name)
        if case is not None:
            selected.append(case)
        elif name in by_base_name:
            selected.extend(by_base_name[name])
        else:
            missing.append(name)
    if missing:
        available = ", ".join(sorted(by_name))
        raise SystemExit(f"unknown case(s): {', '.join(missing)}\navailable: {available}")
    return selected


def _sync_template(template: Path, case_dir: Path, skip_patterns: set[str]) -> None:
    """Copy items from template into case_dir only if they don't exist yet."""
    import fnmatch
    for item in template.iterdir():
        if any(fnmatch.fnmatch(item.name, pat) for pat in skip_patterns):
            continue
        dest = case_dir / item.name
        if not dest.exists():
            if item.is_dir():
                shutil.copytree(item, dest)
            else:
                shutil.copy2(item, dest)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--setup-only", action="store_true", help="prepare cases and exit")
    parser.add_argument("--list-enabled", action="store_true", help="print enabled case names")
    parser.add_argument("--list", action="store_true", help="print selected case names")
    parser.add_argument("--case", action="append", default=[], help="case name to prepare/run; repeatable")
    parser.add_argument("--all", action="store_true", help="include disabled cases too")
    parser.add_argument("--rebuild", action="store_true", help="recreate cases from template")
    args = parser.parse_args()

    cfg = load_config()
    template = ROOT / cfg.get("template", "template_case")
    cases_dir = ROOT / cfg.get("cases_dir", "cases")
    cases_dir.mkdir(parents=True, exist_ok=True)

    selected = select_cases(cfg, args.case, args.all)
    if args.list_enabled:
        for case in [c for c in cfg["cases"] if c.get("enabled", True)]:
            print(case["name"])
        return
    if args.list:
        for case in selected:
            print(case["name"])
        return

    _SKIP = {"__pycache__", "0", "VTK", "VTKs", "processor*", "post-processing-data", "absorptivity_vs_time", "log.*"}

    for case in selected:
        case = dict(case)
        if "write_interval_s" not in case and "write_interval_s" in cfg:
            case["write_interval_s"] = cfg["write_interval_s"]
        if "base_mesh_size_um" not in case and "base_mesh_size_um" in cfg:
            case["base_mesh_size_um"] = cfg["base_mesh_size_um"]
        if "amr_levels" not in case and "amr_levels" in cfg:
            case["amr_levels"] = cfg["amr_levels"]
        case_dir = cases_dir / case["name"]
        case_template = ROOT / case.get("template", cfg.get("template", "template_case"))
        if args.rebuild and case_dir.exists():
            shutil.rmtree(case_dir)
        if not case_dir.exists():
            shutil.copytree(case_template, case_dir, ignore=shutil.ignore_patterns(*_SKIP))
        else:
            _sync_template(case_template, case_dir, _SKIP)
        material = case.get("material")
        if material:
            mat_file = MATERIALS_DIR / material / "transportProperties"
            if not mat_file.exists():
                raise RuntimeError(f"material transportProperties not found: {mat_file}")
            shutil.copy2(mat_file, case_dir / "constant" / "transportProperties")
        if "_electric_resistivity_value" in case:
            set_electric_resistivity(case_dir, case["_electric_resistivity_value"])
        configure_case(case_dir, case)

    if args.setup_only:
        return


if __name__ == "__main__":
    main()
