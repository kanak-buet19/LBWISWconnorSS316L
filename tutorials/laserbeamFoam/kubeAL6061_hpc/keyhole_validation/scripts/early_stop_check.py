#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

import numpy as np
import pandas as pd


def numeric_series(data: pd.DataFrame, column: str) -> np.ndarray:
    if column not in data:
        return np.array([], dtype=float)
    return pd.to_numeric(data[column], errors="coerce").to_numpy(dtype=float)


def longest_stable_duration(time_ms: np.ndarray, stable: np.ndarray, max_gap_ms: float) -> tuple[float, float, float]:
    best = (math.nan, math.nan, 0.0)
    start = None
    for index, ok in enumerate(stable):
        if ok and start is None:
            start = index
        gap_break = start is not None and index > start and time_ms[index] - time_ms[index - 1] > max_gap_ms
        if start is not None and ((not ok) or gap_break):
            end = index - 1
            duration = float(time_ms[end] - time_ms[start])
            if duration > best[2]:
                best = (float(time_ms[start]), float(time_ms[end]), duration)
            start = index if ok and gap_break else None
    if start is not None:
        duration = float(time_ms[-1] - time_ms[start])
        if duration > best[2]:
            best = (float(time_ms[start]), float(time_ms[-1]), duration)
    return best


def read_experimental_ar(config: dict) -> tuple[np.ndarray, np.ndarray]:
    config_dir = Path(config.get("_config_dir", "."))
    exp_path = (config_dir / str(config["exp_aspect_ratio_csv"])).resolve()
    if not exp_path.exists():
        return np.array([], dtype=float), np.array([], dtype=float)
    exp = pd.read_csv(exp_path)
    exp.columns = [column.strip() for column in exp.columns]
    exp_time = pd.to_numeric(exp["time_ms"], errors="coerce").to_numpy(dtype=float)
    exp_ar = pd.to_numeric(exp["meltpool_ar_depthByWidth"], errors="coerce").to_numpy(dtype=float)
    finite = np.isfinite(exp_time) & np.isfinite(exp_ar)
    if not np.any(finite):
        return np.array([], dtype=float), np.array([], dtype=float)
    exp_time = exp_time[finite]
    exp_ar = exp_ar[finite]
    return exp_time - float(np.nanmin(exp_time)), exp_ar


def dimension_gate_config(config: dict) -> dict[str, object]:
    return dict(config["target_stability"].get("dimension_gate", {}))


def dimension_gate_pass(width: np.ndarray, depth: np.ndarray, keyhole: np.ndarray, mask: np.ndarray, config: dict) -> dict[str, object]:
    gate = dimension_gate_config(config)
    if not bool(gate.get("enabled", False)):
        return {"passed": True, "penalty": 0.0, "hard_runaway": False, "keyhole_depth_clipped": False}

    w_low, w_high = [float(value) for value in gate.get("reasonable_width_um", [0.0, math.inf])]
    d_low, d_high = [float(value) for value in gate.get("reasonable_depth_um", [0.0, math.inf])]
    avg_width = float(np.nanmean(width[mask])) if np.any(mask) else math.nan
    avg_depth = float(np.nanmean(depth[mask])) if np.any(mask) else math.nan
    peak_width = float(np.nanmax(width)) if width.size else math.nan
    peak_depth = float(np.nanmax(depth)) if depth.size else math.nan
    peak_keyhole = float(np.nanmax(keyhole[np.isfinite(keyhole)])) if keyhole.size and np.any(np.isfinite(keyhole)) else math.nan
    hard_width = float(gate.get("hard_width_um", math.inf))
    hard_depth = float(gate.get("hard_depth_um", math.inf))
    hard_keyhole = float(gate.get("hard_keyhole_depth_um", math.inf))
    width_error = 0.0 if np.isfinite(avg_width) and w_low <= avg_width <= w_high else 1.0
    depth_error = 0.0 if np.isfinite(avg_depth) and d_low <= avg_depth <= d_high else 1.0
    hard_runaway = (
        (np.isfinite(peak_width) and peak_width >= hard_width)
        or (np.isfinite(peak_depth) and peak_depth >= hard_depth)
        or (np.isfinite(peak_keyhole) and peak_keyhole >= hard_keyhole)
    )
    clipped = np.isfinite(peak_keyhole) and peak_keyhole >= hard_keyhole
    return {
        "passed": width_error == 0.0 and depth_error == 0.0 and not hard_runaway,
        "penalty": width_error + depth_error,
        "hard_runaway": bool(hard_runaway),
        "keyhole_depth_clipped": bool(clipped),
        "avg_width_um": avg_width,
        "avg_depth_um": avg_depth,
        "peak_width_um": peak_width,
        "peak_depth_um": peak_depth,
        "peak_keyhole_depth_um": peak_keyhole,
    }


def aspect_ratio_window_status(time_ms: np.ndarray, width: np.ndarray, depth: np.ndarray, keyhole: np.ndarray, config: dict) -> dict[str, object]:
    target = config["target_stability"]
    window_ms = float(target.get("comparison_window_ms", 2.0))
    ar = np.divide(depth, width, out=np.full_like(depth, np.nan), where=width > 0)
    finite = np.isfinite(time_ms) & np.isfinite(ar)
    gate = dimension_gate_config(config)
    if bool(gate.get("ignore_ar_after_keyhole_clip", False)) and keyhole.size == time_ms.size:
        hard_keyhole = float(gate.get("hard_keyhole_depth_um", math.inf))
        clipped = np.isfinite(keyhole) & (keyhole >= hard_keyhole)
        if np.any(clipped):
            finite &= time_ms < float(np.nanmin(time_ms[clipped]))
    if not np.any(finite):
        return {"ready": False, "reason": "no_ar_rows", "duration_ms": 0.0}

    finite_time = time_ms[finite]
    min_coverage = float(target.get("min_comparison_window_coverage", 0.95))
    max_time = float(np.nanmax(finite_time))
    valid_starts = finite_time[finite_time <= max_time - min_coverage * window_ms]
    if valid_starts.size == 0:
        return {
            "ready": False,
            "reason": "waiting_for_ar_comparison_window",
            "start_ms": float(np.nanmin(finite_time)),
            "end_ms": max_time,
            "duration_ms": max_time - float(np.nanmin(finite_time)),
        }

    exp_time, exp_ar = read_experimental_ar(config)
    if exp_time.size < 2:
        return {"ready": False, "reason": "not_enough_experiment_points", "start_ms": float(valid_starts[0]), "end_ms": max_time, "duration_ms": 0.0}

    best = None
    for start in np.unique(np.round(valid_starts, decimals=9)):
        mask = finite & (time_ms >= float(start)) & (time_ms <= float(start) + window_ms)
        if np.count_nonzero(mask) < 2:
            continue
        duration_ms = float(time_ms[mask][-1] - time_ms[mask][0])
        if duration_ms < min_coverage * window_ms:
            continue
        sim_time = time_ms[mask] - float(time_ms[mask][0])
        sim_ar = ar[mask]
        order = np.argsort(sim_time)
        sim_time = sim_time[order]
        sim_ar = sim_ar[order]
        compare_time = exp_time[(exp_time >= sim_time[0]) & (exp_time <= sim_time[-1])]
        if compare_time.size < 2:
            continue
        sim_interp = np.interp(compare_time, sim_time, sim_ar)
        exp_interp = np.interp(compare_time, exp_time, exp_ar)
        diff = sim_interp - exp_interp
        rmse = float(np.sqrt(np.nanmean(diff * diff)))
        bias = float(np.nanmean(diff))
        slope_error = float(np.polyfit(compare_time, sim_interp, 1)[0] - np.polyfit(compare_time, exp_interp, 1)[0])
        gate_status = dimension_gate_pass(width, depth, keyhole, mask, config)
        score = (
            rmse / max(float(target.get("aspect_ratio_rmse_target", 0.15)), 1.0e-12)
            + 0.35 * abs(bias) / max(float(target.get("aspect_ratio_bias_target", 0.10)), 1.0e-12)
            + 0.35 * abs(slope_error) / max(float(target.get("aspect_ratio_slope_target_per_ms", 0.12)), 1.0e-12)
            + 2.0 * float(gate_status["penalty"])
            + 4.0 * float(gate_status["hard_runaway"])
        )
        candidate = {
            "ready": True,
            "matched": bool(
                rmse <= float(target.get("aspect_ratio_rmse_target", 0.15))
                and abs(bias) <= float(target.get("aspect_ratio_bias_target", 0.10))
                and abs(slope_error) <= float(target.get("aspect_ratio_slope_target_per_ms", 0.12))
                and gate_status["passed"]
            ),
            "start_ms": float(time_ms[mask][0]),
            "end_ms": float(time_ms[mask][-1]),
            "duration_ms": duration_ms,
            "ar_rmse": rmse,
            "ar_bias": bias,
            "ar_slope_error_per_ms": slope_error,
            "dimension_gate_pass": bool(gate_status["passed"]),
            "dimension_gate_penalty": gate_status["penalty"],
            "hard_runaway": bool(gate_status["hard_runaway"]),
            "keyhole_depth_clipped": bool(gate_status["keyhole_depth_clipped"]),
            "comparison_window_score": float(score),
        }
        if best is None or score < float(best["comparison_window_score"]):
            best = candidate

    if best is None:
        return {"ready": False, "reason": "no_overlap_with_experiment", "start_ms": float(valid_starts[0]), "end_ms": max_time, "duration_ms": 0.0}
    return best


def decide(case_dir: Path, config: dict) -> dict[str, object]:
    csv_path = case_dir / "post-processing-data" / "vtu_meltpool_geometry.csv"
    if not csv_path.exists():
        return {"decision": "continue", "reason": "no_csv_yet"}

    data = pd.read_csv(csv_path)
    if data.empty:
        return {"decision": "continue", "reason": "empty_csv"}

    time_ms = numeric_series(data, "time") * 1.0e3
    width = numeric_series(data, "meltPoolWidth_um")
    depth = numeric_series(data, "meltPoolDepth_um")
    keyhole = numeric_series(data, "keyholeDepth_um")
    finite = np.isfinite(time_ms) & np.isfinite(width) & np.isfinite(depth)
    time_ms = time_ms[finite]
    width = width[finite]
    depth = depth[finite]
    keyhole = keyhole[finite]

    if time_ms.size == 0:
        return {"decision": "continue", "reason": "no_finite_rows"}

    latest_time = float(np.nanmax(time_ms))
    latest_width = float(width[-1])
    latest_depth = float(depth[-1])
    latest_keyhole = float(keyhole[-1]) if keyhole.size else math.nan
    early = config.get("early_stop", {})
    target = config["target_stability"]

    ar_status = aspect_ratio_window_status(time_ms, width, depth, keyhole, config)
    if ar_status.get("ready") and ar_status.get("matched"):
        return {
            "decision": "stop_success",
            "reason": "aspect_ratio_trace_matched",
            "latest_time_ms": latest_time,
            "stable_start_ms": ar_status.get("start_ms", math.nan),
            "stable_end_ms": ar_status.get("end_ms", math.nan),
            "stable_duration_ms": ar_status.get("duration_ms", 0.0),
            "ar_rmse": ar_status.get("ar_rmse", math.nan),
            "ar_bias": ar_status.get("ar_bias", math.nan),
            "ar_slope_error_per_ms": ar_status.get("ar_slope_error_per_ms", math.nan),
            "latest_width_um": latest_width,
            "latest_depth_um": latest_depth,
            "latest_keyhole_depth_um": latest_keyhole,
        }

    too_hot = early.get("too_hot", {})
    if latest_time >= float(too_hot.get("min_time_ms", 0.0)):
        if latest_depth > float(too_hot.get("meltpool_depth_um", math.inf)):
            return {
                "decision": "stop_fail",
                "reason": "too_hot_meltpool_depth",
                "latest_time_ms": latest_time,
                "latest_width_um": latest_width,
                "latest_depth_um": latest_depth,
                "latest_keyhole_depth_um": latest_keyhole,
            }
        if np.isfinite(latest_keyhole) and latest_keyhole > float(too_hot.get("keyhole_depth_um", math.inf)):
            return {
                "decision": "stop_fail",
                "reason": "too_hot_keyhole_depth",
                "latest_time_ms": latest_time,
                "latest_width_um": latest_width,
                "latest_depth_um": latest_depth,
                "latest_keyhole_depth_um": latest_keyhole,
            }

    too_cold = early.get("too_cold", {})
    if latest_time >= float(too_cold.get("check_time_ms", math.inf)):
        if latest_width < float(too_cold.get("meltpool_width_um", -math.inf)) and latest_depth < float(too_cold.get("meltpool_depth_um", -math.inf)):
            return {
                "decision": "stop_fail",
                "reason": "too_cold_after_laser_window",
                "latest_time_ms": latest_time,
                "latest_width_um": latest_width,
                "latest_depth_um": latest_depth,
                "latest_keyhole_depth_um": latest_keyhole,
            }

    return {
        "decision": "continue",
        "reason": str(ar_status.get("reason", "still_possible")),
        "latest_time_ms": latest_time,
        "stable_duration_ms": ar_status.get("duration_ms", 0.0),
        "ar_rmse": ar_status.get("ar_rmse", math.nan),
        "ar_bias": ar_status.get("ar_bias", math.nan),
        "ar_slope_error_per_ms": ar_status.get("ar_slope_error_per_ms", math.nan),
        "latest_width_um": latest_width,
        "latest_depth_um": latest_depth,
        "latest_keyhole_depth_um": latest_keyhole,
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--case", type=Path, default=Path.cwd())
    parser.add_argument("--config", type=Path, default=None)
    parser.add_argument("--status-file", type=Path, default=None)
    args = parser.parse_args()

    case_dir = args.case.resolve()
    config_path = args.config or case_dir.parents[1] / "config.json"
    config = json.loads(config_path.read_text())
    config["_config_dir"] = str(config_path.resolve().parent)

    if not bool(config.get("early_stop", {}).get("enabled", True)):
        result = {"decision": "continue", "reason": "early_stop_disabled"}
    else:
        result = decide(case_dir, config)

    status_file = args.status_file or case_dir / "early_stop_status.json"
    status_file.write_text(json.dumps(result, indent=2) + "\n")
    print(result["decision"])


if __name__ == "__main__":
    main()
