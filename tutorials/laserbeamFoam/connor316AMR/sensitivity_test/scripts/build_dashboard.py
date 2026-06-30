#!/usr/bin/env python3
"""Build an HTML dashboard for the transport sensitivity study."""

from __future__ import annotations

import csv
import argparse
import html
import json
import math
from collections import defaultdict
from pathlib import Path
from typing import Any


ROOT = Path(__file__).resolve().parents[1]
CONFIG_PATH = ROOT / "sensitivity_config.json"


def load_config() -> dict[str, Any]:
    return json.loads(CONFIG_PATH.read_text())


def load_rows(config: dict[str, Any]) -> list[dict[str, Any]]:
    corrected_path = ROOT / config["results_root"] / "corrected_surface_widths.csv"
    if corrected_path.exists():
        with corrected_path.open(newline="") as handle:
            return list(csv.DictReader(handle))

    results_path = ROOT / config["results_root"] / "sensitivity_results.csv"
    if not results_path.exists():
        manifest_path = ROOT / config["results_root"] / "case_manifest.json"
        if manifest_path.exists():
            return json.loads(manifest_path.read_text())
        return []
    with results_path.open(newline="") as handle:
        return list(csv.DictReader(handle))


def number(value: Any) -> float | None:
    try:
        if value in (None, ""):
            return None
        out = float(value)
    except (TypeError, ValueError):
        return None
    if not math.isfinite(out):
        return None
    return out


def fmt(value: Any, digits: int = 3) -> str:
    numeric = number(value)
    if numeric is None:
        return "n/a"
    return f"{numeric:.{digits}g}"


def compute_rankings(rows: list[dict[str, Any]]) -> tuple[dict[str, Any] | None, list[dict[str, Any]]]:
    baseline = next((row for row in rows if row.get("case_id") == "baseline"), None)
    width_field = "correctedSurfaceWidth_um" if baseline and number(baseline.get("correctedSurfaceWidth_um")) is not None else "meltPoolWidth_um"
    baseline_width = number(baseline.get(width_field)) if baseline else None
    grouped: dict[str, list[dict[str, Any]]] = defaultdict(list)
    for row in rows:
        if row.get("parameter") not in (None, "", "baseline"):
            grouped[row["parameter"]].append(row)

    ranked = []
    for parameter, group in grouped.items():
        shifts = []
        for row in group:
            width = number(row.get(width_field))
            shift = None if baseline_width is None or width is None else width - baseline_width
            shifts.append((row, shift))
        finite = [abs(shift) for _, shift in shifts if shift is not None]
        max_abs_shift = max(finite) if finite else None
        label = group[0].get("label", parameter)
        ranked.append(
            {
                "parameter": parameter,
                "label": label,
                "max_abs_shift_um": max_abs_shift,
                "rows": shifts,
                "width_field": width_field,
            }
        )
    ranked.sort(key=lambda item: -1 if item["max_abs_shift_um"] is None else item["max_abs_shift_um"], reverse=True)
    return baseline, ranked


def status_badge(status: str | None) -> str:
    status = status or "not_run"
    classes = {
        "complete": "good",
        "skipped_existing": "good",
        "planned": "warn",
        "failed": "bad",
    }
    return f"<span class=\"badge {classes.get(status, 'warn')}\">{html.escape(status)}</span>"


def bar_rows(ranked: list[dict[str, Any]]) -> str:
    max_shift = max((item["max_abs_shift_um"] or 0.0 for item in ranked), default=0.0)
    out = []
    for index, item in enumerate(ranked, start=1):
        shift = item["max_abs_shift_um"]
        pct = 0.0 if not max_shift or shift is None else 100.0 * shift / max_shift
        out.append(
            "<tr>"
            f"<td>{index}</td>"
            f"<td>{html.escape(item['label'])}</td>"
            f"<td>{fmt(shift, 4)}</td>"
            f"<td><div class=\"bar\"><span style=\"width:{pct:.1f}%\"></span></div></td>"
            "</tr>"
        )
    return "\n".join(out)


def case_table(rows: list[dict[str, Any]], baseline_width: float | None) -> str:
    width_field = "correctedSurfaceWidth_um" if any(number(row.get("correctedSurfaceWidth_um")) is not None for row in rows) else "meltPoolWidth_um"
    out = []
    for row in rows:
        width = number(row.get(width_field))
        shift = None if baseline_width is None or width is None else width - baseline_width
        out.append(
            "<tr>"
            f"<td>{html.escape(row.get('case_id', ''))}</td>"
            f"<td>{html.escape(row.get('label', ''))}</td>"
            f"<td>{html.escape(row.get('direction', ''))}</td>"
            f"<td>{fmt(row.get('multiplier'))}</td>"
            f"<td>{fmt(row.get('varied_value'), 5)} {html.escape(row.get('unit', '') or '')}</td>"
            f"<td>{fmt(width, 4)}</td>"
            f"<td>{fmt(shift, 4)}</td>"
            f"<td>{status_badge(row.get('status'))}</td>"
            "</tr>"
        )
    return "\n".join(out)


def perturbation_cards(ranked: list[dict[str, Any]], baseline_width: float | None) -> str:
    cards = []
    for item in ranked:
        low = next((pair for pair in item["rows"] if pair[0].get("direction") == "low"), None)
        high = next((pair for pair in item["rows"] if pair[0].get("direction") == "high"), None)
        def cell(pair: tuple[dict[str, Any], float | None] | None) -> str:
            if pair is None:
                return "n/a"
            row, shift = pair
            width_field = item.get("width_field", "meltPoolWidth_um")
            return f"{fmt(row.get(width_field), 4)} um <small>({fmt(shift, 4)} um)</small>"
        cards.append(
            "<article>"
            f"<h3>{html.escape(item['label'])}</h3>"
            f"<p class=\"metric\">{fmt(item['max_abs_shift_um'], 4)} um</p>"
            "<dl>"
            f"<dt>-20%</dt><dd>{cell(low)}</dd>"
            f"<dt>+20%</dt><dd>{cell(high)}</dd>"
            "</dl>"
            "</article>"
        )
    return "\n".join(cards)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--quiet", action="store_true")
    args = parser.parse_args()

    config = load_config()
    rows = load_rows(config)
    baseline, ranked = compute_rankings(rows)
    width_field = "correctedSurfaceWidth_um" if baseline and number(baseline.get("correctedSurfaceWidth_um")) is not None else "meltPoolWidth_um"
    baseline_width = number(baseline.get(width_field)) if baseline else None
    completed = sum(1 for row in rows if row.get("status") in ("complete", "skipped_existing"))
    planned = len(rows)
    strongest = ranked[0] if ranked and ranked[0]["max_abs_shift_um"] is not None else None
    settings = config["run_settings"]

    css = """
    :root { color-scheme: light; --ink:#17202a; --muted:#667085; --line:#d0d5dd; --bg:#f7f8fa; --panel:#fff; --blue:#2563eb; --good:#0f766e; --warn:#b45309; --bad:#b42318; }
    * { box-sizing: border-box; }
    body { margin:0; font-family: Inter, ui-sans-serif, system-ui, -apple-system, BlinkMacSystemFont, "Segoe UI", sans-serif; color:var(--ink); background:var(--bg); }
    header { padding:28px 32px 18px; background:#fff; border-bottom:1px solid var(--line); }
    h1 { margin:0 0 8px; font-size:28px; letter-spacing:0; }
    h2 { margin:0 0 14px; font-size:18px; letter-spacing:0; }
    h3 { margin:0 0 10px; font-size:14px; letter-spacing:0; }
    p { margin:0; color:var(--muted); }
    main { padding:24px 32px 36px; display:grid; gap:22px; }
    section { background:var(--panel); border:1px solid var(--line); border-radius:8px; padding:18px; }
    .kpis { display:grid; grid-template-columns:repeat(4, minmax(0, 1fr)); gap:14px; }
    .kpi { background:#fff; border:1px solid var(--line); border-radius:8px; padding:16px; }
    .kpi b { display:block; font-size:24px; margin-top:6px; }
    .grid { display:grid; grid-template-columns:1fr 1.25fr; gap:22px; align-items:start; }
    table { width:100%; border-collapse:collapse; font-size:13px; }
    th, td { padding:9px 10px; border-bottom:1px solid #eaecf0; text-align:left; vertical-align:middle; }
    th { color:#475467; font-weight:650; background:#f9fafb; }
    .bar { width:100%; height:10px; background:#e5e7eb; border-radius:999px; overflow:hidden; }
    .bar span { display:block; height:100%; background:var(--blue); }
    .cards { display:grid; grid-template-columns:repeat(auto-fit, minmax(220px, 1fr)); gap:12px; }
    article { border:1px solid var(--line); border-radius:8px; padding:14px; background:#fff; }
    .metric { color:var(--ink); font-size:22px; font-weight:700; margin-bottom:10px; }
    dl { display:grid; grid-template-columns:64px 1fr; gap:6px 10px; margin:0; font-size:13px; }
    dt { color:var(--muted); }
    dd { margin:0; }
    small { color:var(--muted); }
    .badge { display:inline-block; border-radius:999px; padding:3px 8px; font-size:12px; font-weight:650; }
    .good { color:var(--good); background:#ccfbf1; }
    .warn { color:var(--warn); background:#fef3c7; }
    .bad { color:var(--bad); background:#fee4e2; }
    .note { margin-top:10px; font-size:13px; color:var(--muted); }
    @media (max-width: 900px) { main, header { padding-left:16px; padding-right:16px; } .kpis, .grid { grid-template-columns:1fr; } }
    """

    html_text = f"""<!doctype html>
<html lang="en">
<head>
  <meta charset="utf-8">
  <meta name="viewport" content="width=device-width, initial-scale=1">
  <title>316L Transport Sensitivity Dashboard</title>
  <style>{css}</style>
</head>
<body>
  <header>
    <h1>316L Transport Sensitivity Dashboard</h1>
    <p>Short 250 um scan, 40 um base mesh, AMR level 1. Beam diameter, power, and speed held fixed.</p>
  </header>
  <main>
    <div class="kpis">
      <div class="kpi"><p>Cases complete</p><b>{completed}/{planned}</b></div>
      <div class="kpi"><p>Baseline width</p><b>{fmt(baseline_width, 4)} um</b></div>
      <div class="kpi"><p>Strongest parameter</p><b>{html.escape(strongest['label']) if strongest else 'n/a'}</b></div>
      <div class="kpi"><p>Largest width shift</p><b>{fmt(strongest['max_abs_shift_um'] if strongest else None, 4)} um</b></div>
    </div>

    <section>
      <h2>Study Setup</h2>
      <table>
        <tr><th>Setting</th><th>Value</th></tr>
        <tr><td>Base mesh</td><td>{fmt(settings['base_cell_um'])} um, {settings['cells'][0]} x {settings['cells'][1]} x {settings['cells'][2]}</td></tr>
        <tr><td>Domain</td><td>{', '.join(fmt(v * 1e3, 4) for v in settings['domain_m'])} mm</td></tr>
        <tr><td>Scan track</td><td>{fmt(settings['scan_track_m'] * 1e6, 4)} um</td></tr>
        <tr><td>Scan speed</td><td>{fmt(settings['scan_speed_m_per_s'])} m/s</td></tr>
        <tr><td>End time</td><td>{fmt(settings['end_time_s'] * 1e6, 4)} us</td></tr>
        <tr><td>Cooling</td><td>{'convection and radiation disabled' if settings['disable_convective_and_radiative_cooling'] else 'template settings'}</td></tr>
      </table>
    </section>

    <div class="grid">
      <section>
        <h2>Ranked Width Sensitivity</h2>
        <table>
          <thead><tr><th>#</th><th>Parameter</th><th>Max |shift| (um)</th><th>Relative bar</th></tr></thead>
          <tbody>{bar_rows(ranked)}</tbody>
        </table>
        <p class="note">Ranking uses the larger absolute width shift from the -20% and +20% cases. Width source: {html.escape(width_field)}.</p>
      </section>
      <section>
        <h2>Parameter Cards</h2>
        <div class="cards">{perturbation_cards(ranked, baseline_width)}</div>
      </section>
    </div>

    <section>
      <h2>Case Results</h2>
      <table>
        <thead><tr><th>Case</th><th>Parameter</th><th>Direction</th><th>Multiplier</th><th>Value</th><th>Width (um)</th><th>Shift (um)</th><th>Status</th></tr></thead>
        <tbody>{case_table(rows, baseline_width)}</tbody>
      </table>
    </section>
  </main>
</body>
</html>
"""
    out_path = ROOT / config["dashboard"]
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text(html_text)
    if not args.quiet:
        print(f"Wrote {out_path}")


if __name__ == "__main__":
    main()
