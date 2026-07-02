"""Plotting for CW calibration results. All functions are best-effort and skip
gracefully when data is missing."""

from __future__ import annotations

from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


def _finite_pairs(xs, ys):
    out = [(x, y) for x, y in zip(xs, ys) if x == x and y == y]
    return ([p[0] for p in out], [p[1] for p in out])


def plot_ranked_objective(summary: list[dict], out: Path) -> None:
    done = [c for c in summary if c.get("objective") == c.get("objective")]
    if not done:
        return
    done = sorted(done, key=lambda c: c["objective"])
    labels = [f"#{c['id']:02d}" for c in done]
    vals = [c["objective"] for c in done]
    colors = ["#2ca02c" if i == 0 else "#4C78A8" for i in range(len(done))]
    fig, ax = plt.subplots(figsize=(max(6, len(done) * 0.5), 4))
    ax.bar(labels, vals, color=colors)
    ax.set_ylabel("objective (mean abs rel error)")
    ax.set_xlabel("candidate (ranked)")
    ax.set_title("Calibration objective by candidate (green = best)")
    ax.axhline(0.10, ls="--", c="grey", lw=1, label="10% line")
    ax.legend()
    plt.xticks(rotation=90)
    fig.tight_layout()
    fig.savefig(out, dpi=130)
    plt.close(fig)


def plot_params_vs_objective(summary: list[dict], param_names: list[str], out: Path) -> None:
    done = [c for c in summary if c.get("objective") == c.get("objective")]
    if not done:
        return
    obj = [c["objective"] for c in done]
    ncol = min(3, len(param_names))
    nrow = int(np.ceil(len(param_names) / ncol))
    fig, axes = plt.subplots(nrow, ncol, figsize=(4 * ncol, 3.2 * nrow), squeeze=False)
    for i, p in enumerate(param_names):
        ax = axes[i // ncol][i % ncol]
        xs = [c["params"].get(p, float("nan")) for c in done]
        x, y = _finite_pairs(xs, obj)
        sc = ax.scatter(x, y, c=y, cmap="viridis_r", s=40, edgecolor="k", linewidth=0.3)
        ibest = int(np.argmin(obj))
        ax.scatter([xs[ibest]], [obj[ibest]], marker="*", s=220, color="red",
                   edgecolor="k", zorder=5, label="best")
        ax.set_xlabel(p)
        ax.set_ylabel("objective")
        ax.legend(fontsize=8)
    for j in range(len(param_names), nrow * ncol):
        axes[j // ncol][j % ncol].axis("off")
    fig.suptitle("Parameter value vs objective")
    fig.colorbar(sc, ax=axes, shrink=0.6, label="objective")
    fig.savefig(out, dpi=130)
    plt.close(fig)


def plot_best_bars(best: dict, out: Path) -> None:
    cases = best.get("cases", {})
    if not cases:
        return
    names = list(cases.keys())
    fig, axes = plt.subplots(1, 2, figsize=(10, 4))
    for ax, metric, title in ((axes[0], "depth", "Melt-pool depth (um)"),
                              (axes[1], "width", "Melt-pool width (um)")):
        sim = [cases[n].get(f"sim_{metric}_um", float("nan")) for n in names]
        exp = [cases[n].get(f"exp_{metric}_um", float("nan")) for n in names]
        x = np.arange(len(names))
        ax.bar(x - 0.2, sim, 0.4, label="sim", color="#4C78A8")
        ax.bar(x + 0.2, exp, 0.4, label="exp", color="#F58518")
        ax.set_xticks(x)
        ax.set_xticklabels(names, rotation=15, fontsize=8)
        ax.set_title(title)
        ax.legend()
    fig.suptitle(f"Best candidate #{best.get('id')} (objective {best.get('objective'):.3f})")
    fig.tight_layout()
    fig.savefig(out, dpi=130)
    plt.close(fig)


def plot_best_convergence(best: dict, out: Path) -> None:
    cases = best.get("cases", {})
    series_any = any(cases[n].get("series") for n in cases)
    if not series_any:
        return
    fig, axes = plt.subplots(1, 2, figsize=(11, 4))
    for n, c in cases.items():
        s = c.get("series", [])
        if not s:
            continue
        t = [r["time"] * 1e3 for r in s]  # ms
        axes[0].plot(t, [r["depth"] for r in s], "-o", label=n)
        axes[1].plot(t, [r["width"] for r in s], "-o", label=n)
        axes[0].axhline(c.get("exp_depth_um", float("nan")), ls="--", lw=1)
        axes[1].axhline(c.get("exp_width_um", float("nan")), ls="--", lw=1)
    axes[0].set_title("Depth vs time (dashed = exp)")
    axes[1].set_title("Width vs time (dashed = exp)")
    for ax in axes:
        ax.set_xlabel("time (ms)")
        ax.set_ylabel("um")
        ax.legend(fontsize=8)
    fig.suptitle(f"Best candidate #{best.get('id')} melt-pool stabilization")
    fig.tight_layout()
    fig.savefig(out, dpi=130)
    plt.close(fig)


def make_all(summary: list[dict], best: dict | None, param_names: list[str],
             out_dir: Path) -> None:
    out_dir.mkdir(parents=True, exist_ok=True)
    try:
        plot_ranked_objective(summary, out_dir / "ranked_objective.png")
        plot_params_vs_objective(summary, param_names, out_dir / "params_vs_objective.png")
        if best:
            plot_best_bars(best, out_dir / "best_sim_vs_exp.png")
            plot_best_convergence(best, out_dir / "best_convergence.png")
    except Exception as exc:  # never let plotting kill the run
        print(f"[WARN] plotting failed: {exc}")
