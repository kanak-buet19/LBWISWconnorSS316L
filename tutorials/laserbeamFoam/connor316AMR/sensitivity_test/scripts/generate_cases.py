#!/usr/bin/env python3
"""Generate short-track transport-property sensitivity cases."""

from __future__ import annotations

import json
import argparse
import re
import shutil
from dataclasses import dataclass
from pathlib import Path
from typing import Any


ROOT = Path(__file__).resolve().parents[1]
CONFIG_PATH = ROOT / "sensitivity_config.json"


@dataclass(frozen=True)
class CaseVariant:
    case_id: str
    parameter: str
    label: str
    direction: str
    multiplier: float


def load_config() -> dict[str, Any]:
    return json.loads(CONFIG_PATH.read_text())


def format_number(value: float) -> str:
    return f"{value:.10g}"


def replace_entry(text: str, entry: str, value: str, required: bool = True) -> str:
    pattern = re.compile(rf"(^\s*{re.escape(entry)}\s+)([^;]+)(;)", re.MULTILINE)
    new_text, count = pattern.subn(rf"\g<1>{value}\3", text, count=1)
    if count == 0 and required:
        raise ValueError(f"Entry {entry} not found")
    return new_text


def replace_block_entries(text: str, block: str, replacements: dict[str, str]) -> str:
    start = re.search(rf"(^\s*{re.escape(block)}\s*\{{)", text, re.MULTILINE)
    if not start:
        raise ValueError(f"Block {block} not found")

    brace_index = text.find("{", start.start())
    depth = 0
    end_index = None
    for index in range(brace_index, len(text)):
        if text[index] == "{":
            depth += 1
        elif text[index] == "}":
            depth -= 1
            if depth == 0:
                end_index = index + 1
                break
    if end_index is None:
        raise ValueError(f"Block {block} is not closed")

    block_text = text[start.start():end_index]
    for entry, value in replacements.items():
        block_text = replace_entry(block_text, entry, value)
    return text[:start.start()] + block_text + text[end_index:]


def scalar_value(text: str, path: str, optional: bool = False) -> float | None:
    if "." in path:
        block, entry = path.split(".", 1)
        match = re.search(
            rf"\b{re.escape(block)}\s*\{{(?P<body>.*?)^\s*\}}",
            text,
            flags=re.MULTILINE | re.DOTALL,
        )
        if not match:
            if optional:
                return None
            raise ValueError(f"Block {block} not found")
        search_text = match.group("body")
    else:
        entry = path
        search_text = text

    match = re.search(rf"^\s*{re.escape(entry)}\s+([^;]+);", search_text, flags=re.MULTILINE)
    if not match:
        if optional:
            return None
        raise ValueError(f"Entry {path} not found")
    return float(match.group(1))


def scale_scalar(text: str, path: str, multiplier: float, optional: bool = False) -> tuple[str, float | None, float | None]:
    baseline = scalar_value(text, path, optional=optional)
    if baseline is None:
        return text, None, None
    varied = baseline * multiplier
    value_text = format_number(varied)
    if "." in path:
        block, entry = path.split(".", 1)
        text = replace_block_entries(text, block, {entry: value_text})
    else:
        text = replace_entry(text, path, value_text)
    return text, baseline, varied


def scale_table_by_temperature(
    text: str,
    path: str,
    multiplier: float,
    temperature_min: float | None,
    temperature_max: float | None,
) -> tuple[str, float, float]:
    if "." not in path:
        raise ValueError(f"Table path must include block.entry: {path}")

    block, table_name = path.split(".", 1)
    block_match = re.search(
        rf"\b{re.escape(block)}\s*\{{(?P<body>.*?)^\s*\}}",
        text,
        flags=re.MULTILINE | re.DOTALL,
    )
    if not block_match:
        raise ValueError(f"Block {block} not found")

    body = block_match.group("body")
    table_match = re.search(
        rf"(?P<head>^\s*{re.escape(table_name)}[^\n]*\n\s*\(\n)(?P<rows>.*?)(?P<tail>^\s*\)\s*;)",
        body,
        flags=re.MULTILINE | re.DOTALL,
    )
    if not table_match:
        raise ValueError(f"Table {path} not found")

    row_pattern = re.compile(r"(?P<prefix>\(\s*)(?P<T>[-+0-9.eE]+)(?P<mid>\s+)(?P<V>[-+0-9.eE]+)(?P<suffix>\s*\))")
    changed_values: list[tuple[float, float]] = []

    def replace_row(match: re.Match[str]) -> str:
        temperature = float(match.group("T"))
        value = float(match.group("V"))
        in_lower = temperature_min is None or temperature >= temperature_min
        in_upper = temperature_max is None or temperature <= temperature_max
        if in_lower and in_upper:
            new_value = value * multiplier
            changed_values.append((value, new_value))
            return (
                f"{match.group('prefix')}{match.group('T')}"
                f"{match.group('mid')}{format_number(new_value)}{match.group('suffix')}"
            )
        return match.group(0)

    rows = row_pattern.sub(replace_row, table_match.group("rows"))
    if not changed_values:
        raise ValueError(f"No rows selected for {path}")

    new_table = table_match.group("head") + rows + table_match.group("tail")
    new_body = body[:table_match.start()] + new_table + body[table_match.end():]
    new_text = text[:block_match.start("body")] + new_body + text[block_match.end("body"):]
    baseline_mean = sum(value for value, _ in changed_values) / len(changed_values)
    varied_mean = sum(value for _, value in changed_values) / len(changed_values)
    return new_text, baseline_mean, varied_mean


def patch_block_mesh(case_dir: Path, settings: dict[str, Any]) -> None:
    path = case_dir / "system" / "blockMeshDict"
    text = path.read_text()
    x_len, y_len, z_len = settings["domain_m"]
    nx, ny, nz = settings["cells"]
    text = re.sub(
        r"Base mesh spacing is uniform: .*",
        f"Base mesh spacing is uniform: {format_number(settings['base_cell_um'])} um in x, y, and z.",
        text,
    )
    vertices = f"""vertices
(
    (0       0       0)
    ({format_number(x_len)}  0       0)
    ({format_number(x_len)}  {format_number(y_len)}  0)
    (0       {format_number(y_len)}  0)
    (0       0       {format_number(z_len)})
    ({format_number(x_len)}  0       {format_number(z_len)})
    ({format_number(x_len)}  {format_number(y_len)}  {format_number(z_len)})
    (0       {format_number(y_len)}  {format_number(z_len)})
);"""
    text = re.sub(r"vertices\s*\(.*?\);", vertices, text, flags=re.DOTALL)
    text = re.sub(
        r"hex \(0 1 2 3 4 5 6 7\) \([^)]+\)",
        f"hex (0 1 2 3 4 5 6 7) ({nx} {ny} {nz})",
        text,
    )
    path.write_text(text)


def patch_set_fields(case_dir: Path, settings: dict[str, Any]) -> None:
    path = case_dir / "system" / "setFieldsDict"
    text = path.read_text()
    x_len, y_len, z_len = settings["domain_m"]
    surface_y = settings["surface_y_m"]
    replacement = f"box (0 {format_number(surface_y)} 0) ({format_number(x_len)} {format_number(y_len)} {format_number(z_len)});"
    text = re.sub(r"box\s+\([^)]+\)\s+\([^)]+\);", replacement, text, count=1)
    path.write_text(text)


def patch_dynamic_mesh(case_dir: Path, settings: dict[str, Any]) -> None:
    path = case_dir / "constant" / "dynamicMeshDict"
    text = path.read_text()
    text = replace_entry(text, "maxRefinement", str(settings["max_refinement"]))
    text = replace_entry(text, "refineInterval", str(settings["refine_interval"]))
    path.write_text(text)


def patch_control_dict(case_dir: Path, settings: dict[str, Any]) -> None:
    path = case_dir / "system" / "controlDict"
    text = path.read_text()
    text = replace_entry(text, "startFrom", "startTime")
    text = replace_entry(text, "endTime", format_number(settings["end_time_s"]))
    text = replace_entry(text, "writeInterval", format_number(settings["write_interval_s"]))
    text = replace_entry(text, "purgeWrite", "1")
    path.write_text(text)


def patch_decompose(case_dir: Path, settings: dict[str, Any]) -> None:
    path = case_dir / "system" / "decomposeParDict"
    text = path.read_text()
    text = replace_entry(text, "numberOfSubdomains", str(settings["number_of_subdomains"]))
    path.write_text(text)


def patch_laser_tables(case_dir: Path, settings: dict[str, Any]) -> None:
    center_x, center_y, center_z = settings["scan_center_m"]
    track = settings["scan_track_m"]
    speed = settings["scan_speed_m_per_s"]
    end_time = settings["end_time_s"]
    z0 = center_z - 0.5 * track
    z1 = center_z + 0.5 * track
    (case_dir / "constant" / "timeVsLaserPosition").write_text(
        "(\n"
        f"    (0       ({format_number(center_x)} {format_number(center_y)} {format_number(z0)}))\n"
        f"    ({format_number(end_time)}     ({format_number(center_x)} {format_number(center_y)} {format_number(z1)}))\n"
        ")\n"
    )

    transport = case_dir / "constant" / "transportProperties"
    t_text = transport.read_text()
    t_text = replace_entry(t_text, "V_scan", format_number(speed), required=False)
    transport.write_text(t_text)


def patch_transport_defaults(case_dir: Path, settings: dict[str, Any]) -> None:
    if not settings.get("disable_convective_and_radiative_cooling", False):
        return
    path = case_dir / "constant" / "transportProperties"
    text = path.read_text()
    text = replace_entry(text, "h_convection", "0", required=False)
    text = replace_entry(text, "emS", "0", required=False)
    text = replace_entry(text, "emL", "0", required=False)
    path.write_text(text)


def patch_allrun(case_dir: Path, settings: dict[str, Any]) -> None:
    path = case_dir / "Allrun_long"
    text = path.read_text()
    text = re.sub(r"^convertVTK=.*$", f"convertVTK={str(settings['convert_vtk_during_run']).lower()}", text, flags=re.MULTILINE)
    text = re.sub(r"^analyzeMeltpool=.*$", f"analyzeMeltpool={str(settings['analyze_meltpool_during_run']).lower()}", text, flags=re.MULTILINE)
    path.write_text(text)


def apply_base_patches(case_dir: Path, settings: dict[str, Any]) -> None:
    patch_block_mesh(case_dir, settings)
    patch_set_fields(case_dir, settings)
    patch_dynamic_mesh(case_dir, settings)
    patch_control_dict(case_dir, settings)
    patch_decompose(case_dir, settings)
    patch_laser_tables(case_dir, settings)
    patch_transport_defaults(case_dir, settings)
    patch_allrun(case_dir, settings)


def variants(config: dict[str, Any]) -> list[CaseVariant]:
    delta = float(config["sweep"]["relative_delta"])
    out = [CaseVariant("baseline", "baseline", "Baseline", "baseline", 1.0)]
    for parameter in config["sweep"]["parameters"]:
        for direction, multiplier in (("low", 1.0 - delta), ("high", 1.0 + delta)):
            out.append(
                CaseVariant(
                    f"{parameter['name']}__{direction}",
                    parameter["name"],
                    parameter["label"],
                    direction,
                    multiplier,
                )
            )
    return out


def copy_template(template_dir: Path, case_dir: Path) -> None:
    if case_dir.exists():
        shutil.rmtree(case_dir)
    ignore = shutil.ignore_patterns(
        "0",
        "[0-9]*",
        "processor*",
        "VTK",
        "VTKs",
        "post-processing-data",
        "log.*",
        "*.foam",
        "__pycache__",
    )
    shutil.copytree(template_dir, case_dir, ignore=ignore)


def overlay_transport_properties(config: dict[str, Any], case_dir: Path) -> None:
    override = config.get("transport_properties_override")
    if not override:
        return
    source = (ROOT / override).resolve()
    if not source.exists():
        raise FileNotFoundError(f"transport_properties_override not found: {source}")
    shutil.copy2(source, case_dir / "constant" / "transportProperties")


def apply_variant(config: dict[str, Any], case_dir: Path, variant: CaseVariant) -> dict[str, Any]:
    metadata: dict[str, Any] = {
        "case_id": variant.case_id,
        "parameter": variant.parameter,
        "label": variant.label,
        "direction": variant.direction,
        "multiplier": variant.multiplier,
        "baseline_value": None,
        "varied_value": None,
        "unit": "",
    }
    if variant.parameter == "baseline":
        return metadata

    parameter = next(item for item in config["sweep"]["parameters"] if item["name"] == variant.parameter)
    path = case_dir / parameter["file"]
    text = path.read_text()
    optional = bool(parameter.get("optional", False))

    if parameter["operation"] == "scale_scalar":
        text, baseline, varied = scale_scalar(text, parameter["path"], variant.multiplier, optional=optional)
    elif parameter["operation"] == "scale_table_by_temperature":
        text, baseline, varied = scale_table_by_temperature(
            text,
            parameter["path"],
            variant.multiplier,
            parameter.get("temperature_min_K"),
            parameter.get("temperature_max_K"),
        )
    else:
        raise ValueError(f"Unsupported operation: {parameter['operation']}")

    path.write_text(text)
    metadata.update(
        {
            "baseline_value": baseline,
            "varied_value": varied,
            "unit": parameter.get("unit", ""),
        }
    )
    return metadata


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--quiet", action="store_true")
    args = parser.parse_args()

    config = load_config()
    template_dir = (ROOT / config["template_case"]).resolve()
    case_root = ROOT / config["case_root"]
    results_root = ROOT / config["results_root"]
    case_root.mkdir(parents=True, exist_ok=True)
    results_root.mkdir(parents=True, exist_ok=True)

    active_variants = variants(config)
    active_case_ids = {variant.case_id for variant in active_variants}
    for old_case in case_root.iterdir():
        if old_case.is_dir() and old_case.name not in active_case_ids:
            shutil.rmtree(old_case)

    metadata_rows = []
    for variant in active_variants:
        case_dir = case_root / variant.case_id
        copy_template(template_dir, case_dir)
        overlay_transport_properties(config, case_dir)
        apply_base_patches(case_dir, config["run_settings"])
        metadata = apply_variant(config, case_dir, variant)
        metadata["case_dir"] = str(case_dir.relative_to(ROOT))
        metadata_rows.append(metadata)

    (results_root / "case_manifest.json").write_text(json.dumps(metadata_rows, indent=2))
    if not args.quiet:
        print(f"Generated {len(metadata_rows)} cases in {case_root}")


if __name__ == "__main__":
    main()
