#!/usr/bin/env python3
"""Compare melt-pool metrics across mesh convergence cases.

Reads post-processing-data/vtu_meltpool_geometry.csv from each case,
extracts quasi-steady values (last 25% of simulation time), and plots
width + depth vs. fine cell size.

Output: results/mesh_convergence_comparison.pdf + results/mesh_convergence_summary.csv
"""

from __future__ import annotations

import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

SCRIPT_DIR = Path(__file__).resolve().parent
MESH_CONV_DIR = SCRIPT_DIR.parent
CONFIG_FILE = MESH_CONV_DIR / "mesh_convergence_config.json"
RESULTS_DIR = MESH_CONV_DIR / "results"

METRICS = [
    ("meltPoolWidth_um", "Melt-pool width (µm)"),
    ("meltPoolDepth_um", "Melt-pool depth (µm)"),
    ("keyholeWidth_um", "Keyhole width (µm)"),
    ("keyholeDepth_um", "Keyhole depth (µm)"),
]

COLORS = {
    "mesh_4um": "#E45756",
    "mesh_6um": "#F58518",
    "mesh_8um": "#4C78A8",
}

MARKERS = {
    "mesh_4um": "^",
    "mesh_6um": "s",
    "mesh_8um": "o",
}


def load_config() -> dict:
    with CONFIG_FILE.open() as fh:
        return json.load(fh)


def load_case_csv(case_dir: Path) -> pd.DataFrame | None:
    csv_path = case_dir / "post-processing-data" / "vtu_meltpool_geometry.csv"
    if not csv_path.exists():
        print(f"  [WARN] CSV not found: {csv_path}")
        return None
    df = pd.read_csv(csv_path)
    if df.empty:
        print(f"  [WARN] Empty CSV: {csv_path}")
        return None
    df = df.sort_values("time").reset_index(drop=True)
    return df


def quasi_steady_mean(df: pd.DataFrame, col: str, fraction: float = 0.25) -> float:
    """Mean over last `fraction` of time range, ignoring NaN."""
    t_min = df["time"].min()
    t_max = df["time"].max()
    t_cutoff = t_max - fraction * (t_max - t_min)
    subset = df.loc[df["time"] >= t_cutoff, col].dropna()
    if subset.empty:
        return np.nan
    return float(subset.mean())


def quasi_steady_std(df: pd.DataFrame, col: str, fraction: float = 0.25) -> float:
    t_min = df["time"].min()
    t_max = df["time"].max()
    t_cutoff = t_max - fraction * (t_max - t_min)
    subset = df.loc[df["time"] >= t_cutoff, col].dropna()
    if subset.empty or len(subset) < 2:
        return 0.0
    return float(subset.std())


def main() -> None:
    config = load_config()
    RESULTS_DIR.mkdir(exist_ok=True)

    cases_cfg = config["cases"]
    summary_rows = []

    # Gather data
    data: dict[str, dict] = {}
    for case_name, cfg in cases_cfg.items():
        case_dir = MESH_CONV_DIR / case_name
        print(f"Loading {case_name} ...")
        df = load_case_csv(case_dir)
        if df is None:
            print(f"  [SKIP] {case_name} has no data yet.")
            continue

        row = {
            "case": case_name,
            "fine_cell_um": cfg["fine_cell_um"],
            "base_cell_um": cfg["base_cell_um"],
            "amr_levels": cfg["amr_levels"],
            "n_timesteps": len(df),
        }
        for metric, _ in METRICS:
            row[f"{metric}_mean"] = quasi_steady_mean(df, metric)
            row[f"{metric}_std"] = quasi_steady_std(df, metric)

        data[case_name] = {"cfg": cfg, "df": df, "row": row}
        summary_rows.append(row)

    if not summary_rows:
        print("No case data found. Run simulations first.")
        return

    summary_df = pd.DataFrame(summary_rows).sort_values("fine_cell_um")
    summary_csv = RESULTS_DIR / "mesh_convergence_summary.csv"
    summary_df.to_csv(summary_csv, index=False)
    print(f"\nSummary written to {summary_csv}")
    print(summary_df[["case", "fine_cell_um"] + [f"{m[0]}_mean" for m in METRICS]].to_string(index=False))

    # ---- Figures ----
    fig, axes = plt.subplots(2, 2, figsize=(12, 9))
    axes = axes.flatten()

    x_all = sorted({cfg["fine_cell_um"] for cfg in cases_cfg.values()})

    for ax_idx, (metric, ylabel) in enumerate(METRICS):
        ax = axes[ax_idx]

        x_vals, y_means, y_stds, labels = [], [], [], []
        for case_name in sorted(data, key=lambda k: data[k]["cfg"]["fine_cell_um"]):
            row = data[case_name]["row"]
            fine_um = row["fine_cell_um"]
            mean_val = row[f"{metric}_mean"]
            std_val = row[f"{metric}_std"]

            if np.isfinite(mean_val):
                x_vals.append(fine_um)
                y_means.append(mean_val)
                y_stds.append(std_val)
                labels.append(case_name)

                ax.errorbar(
                    fine_um,
                    mean_val,
                    yerr=std_val if std_val > 0 else None,
                    fmt=MARKERS.get(case_name, "o"),
                    color=COLORS.get(case_name, "gray"),
                    markersize=9,
                    capsize=5,
                    linewidth=1.5,
                    label=f"{case_name} ({fine_um}µm)",
                    zorder=5,
                )

        if len(x_vals) >= 2:
            sort_idx = np.argsort(x_vals)
            ax.plot(
                np.array(x_vals)[sort_idx],
                np.array(y_means)[sort_idx],
                color="gray",
                linewidth=1.2,
                linestyle="--",
                zorder=3,
            )

        ax.set_xlabel("Fine cell size (µm)")
        ax.set_ylabel(ylabel)
        ax.set_title(ylabel)
        ax.set_xticks(x_all)
        ax.set_xticklabels([f"{x}µm" for x in x_all])
        ax.legend(frameon=True, facecolor="white", edgecolor="black", fontsize=8)
        ax.grid(True, linestyle=":", alpha=0.5)

        if not x_vals:
            ax.text(0.5, 0.5, "No data", transform=ax.transAxes, ha="center", va="center", color="gray")

    # ---- Time-history overlay ----
    fig2, axes2 = plt.subplots(2, 2, figsize=(14, 9))
    axes2 = axes2.flatten()

    for ax_idx, (metric, ylabel) in enumerate(METRICS):
        ax = axes2[ax_idx]
        has_data = False

        for case_name in sorted(data, key=lambda k: data[k]["cfg"]["fine_cell_um"]):
            df = data[case_name]["df"]
            fine_um = data[case_name]["cfg"]["fine_cell_um"]

            col_vals = df[metric].dropna() if metric in df.columns else pd.Series(dtype=float)
            if col_vals.empty:
                continue

            t_us = df.loc[col_vals.index, "time"] * 1e6
            ax.plot(
                t_us,
                col_vals,
                color=COLORS.get(case_name, "gray"),
                label=f"{case_name} ({fine_um}µm)",
                linewidth=1.5,
                marker=MARKERS.get(case_name, "o"),
                markersize=3,
            )
            has_data = True

        ax.set_xlabel("Time (µs)")
        ax.set_ylabel(ylabel)
        ax.set_title(f"{ylabel} — time history")
        if has_data:
            ax.legend(frameon=True, facecolor="white", edgecolor="black", fontsize=8)
        ax.grid(True, linestyle=":", alpha=0.5)

    fig.suptitle("Mesh convergence — quasi-steady melt-pool metrics\n(P=200W, v=900mm/s, 316L SS)", fontsize=12)
    fig.tight_layout()

    fig2.suptitle("Mesh convergence — melt-pool metric time histories\n(P=200W, v=900mm/s, 316L SS)", fontsize=12)
    fig2.tight_layout()

    out_pdf = RESULTS_DIR / "mesh_convergence_comparison.pdf"
    out_pdf2 = RESULTS_DIR / "mesh_convergence_time_history.pdf"
    fig.savefig(out_pdf, format="pdf", bbox_inches="tight", dpi=150)
    fig2.savefig(out_pdf2, format="pdf", bbox_inches="tight", dpi=150)
    plt.close("all")

    print(f"\nPlots written to:")
    print(f"  {out_pdf}")
    print(f"  {out_pdf2}")


if __name__ == "__main__":
    main()
