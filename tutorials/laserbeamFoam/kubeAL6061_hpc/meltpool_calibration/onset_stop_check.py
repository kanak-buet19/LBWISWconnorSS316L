#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

import pandas as pd


def write_reason(case_dir: Path, payload: dict[str, object]) -> None:
    (case_dir / "STOP_REASON.json").write_text(json.dumps(payload, indent=2) + "\n")


def keyhole_start_ms(data: pd.DataFrame, threshold_um: float) -> float:
    if "time" not in data or "keyholeDepth_um" not in data:
        return math.nan

    time_s = pd.to_numeric(data["time"], errors="coerce")
    depth_um = pd.to_numeric(data["keyholeDepth_um"], errors="coerce")
    active = data[(depth_um > threshold_um) & time_s.notna()]
    if active.empty:
        return math.nan

    return float(pd.to_numeric(active["time"], errors="coerce").iloc[0]) * 1.0e3


def latest_time_ms(data: pd.DataFrame) -> float:
    if "time" not in data:
        return math.nan
    values = pd.to_numeric(data["time"], errors="coerce").dropna()
    if values.empty:
        return math.nan
    return float(values.max()) * 1.0e3


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--case", type=Path, required=True)
    parser.add_argument("--config", type=Path, required=True)
    args = parser.parse_args()

    case_dir = args.case.resolve()
    config_path = args.config
    if not config_path.exists() and not config_path.is_absolute():
        fallback = Path(__file__).resolve().parent / config_path.name
        if fallback.exists():
            config_path = fallback

    if not config_path.exists():
        print("continue")
        return

    config = json.loads(config_path.read_text())
    csv_path = case_dir / "post-processing-data" / "vtu_meltpool_geometry.csv"
    if not csv_path.exists():
        print("continue")
        return

    data = pd.read_csv(csv_path)
    if data.empty:
        print("continue")
        return

    target_ms = float(config.get("target_keyhole_start_ms", 4.2709))
    tolerance_ms = float(config.get("keyhole_start_tolerance_ms", 0.05))
    threshold_um = float(config.get("keyhole_onset_threshold_um", 1.0))
    onset_ms = keyhole_start_ms(data, threshold_um)
    latest_ms = latest_time_ms(data)

    if math.isfinite(onset_ms):
        if onset_ms < target_ms - tolerance_ms:
            decision = "stop_early"
            reason = "keyhole_started_too_early"
        elif onset_ms > target_ms + tolerance_ms:
            decision = "stop_late"
            reason = "keyhole_started_too_late"
        else:
            decision = "stop_matched"
            reason = "keyhole_start_within_tolerance"

        write_reason(
            case_dir,
            {
                "decision": decision,
                "reason": reason,
                "target_keyhole_start_ms": target_ms,
                "tolerance_ms": tolerance_ms,
                "keyhole_start_ms": onset_ms,
                "keyhole_start_error_ms": onset_ms - target_ms,
                "latest_time_ms": latest_ms,
            },
        )
        print(decision)
        return

    if math.isfinite(latest_ms) and latest_ms > target_ms + tolerance_ms:
        write_reason(
            case_dir,
            {
                "decision": "stop_late",
                "reason": "no_keyhole_by_late_bound",
                "target_keyhole_start_ms": target_ms,
                "tolerance_ms": tolerance_ms,
                "keyhole_start_ms": math.nan,
                "latest_time_ms": latest_ms,
            },
        )
        print("stop_late")
        return

    print("continue")


if __name__ == "__main__":
    main()
