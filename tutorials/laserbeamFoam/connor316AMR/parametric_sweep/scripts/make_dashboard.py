#!/usr/bin/env python3
"""Create a central PNG dashboard for the parametric sweep results."""

from __future__ import annotations

import argparse
import csv
import json
import math
import re
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.image as mpimg  # noqa: E402
import matplotlib.pyplot as plt  # noqa: E402


ROOT = Path(__file__).resolve().parents[1]


def finite_float(value: object) -> float | None:
    try:
        number = float(value)
    except (TypeError, ValueError):
        return None
    return number if math.isfinite(number) else None


def read_manifest(path: Path) -> list[dict[str, object]]:
    if not path.exists():
        return []
    return json.loads(path.read_text())


def read_metrics(case_dir: Path) -> tuple[float | None, float | None, float | None]:
    csv_path = case_dir / "post-processing-data" / "vtk_meltpool_geometry.csv"
    if not csv_path.exists():
        return None, None, None

    best_time = None
    best_depth = None
    best_width = None
    with csv_path.open(newline="") as handle:
        for row in csv.DictReader(handle):
            depth = finite_float(row.get("keyholeDepth_um"))
            if depth is None:
                continue
            if best_depth is None or depth > best_depth:
                best_depth = depth
                best_width = finite_float(row.get("keyholeWidth_um"))
                best_time = finite_float(row.get("time"))
    return best_depth, best_width, best_time


def png_time(path: Path) -> float | None:
    match = re.search(r"_([0-9.]+e[+-]?[0-9]+)\.png$", path.name)
    return finite_float(match.group(1)) if match else None


def section_png(case_dir: Path, target_time: float | None) -> Path | None:
    pngs = sorted((case_dir / "post-processing-data" / "vtk_sections").glob("*.png"))
    if not pngs:
        return None
    if target_time is not None:
        timed_pngs = [(path, png_time(path)) for path in pngs]
        timed_pngs = [(path, time_value) for path, time_value in timed_pngs if time_value is not None]
        if timed_pngs:
            return min(timed_pngs, key=lambda item: abs(item[1] - target_time))[0]
    return max(pngs, key=lambda path: path.stat().st_mtime)


def xy_crop(image):
    height, width = image.shape[:2]
    return image[: int(height * 0.58), int(width * 0.255) : int(width * 0.63)]


def build_dashboard(root: Path, output: Path) -> None:
    manifest = read_manifest(root / "results" / "case_manifest.json")
    if not manifest:
        manifest = [{"case_id": path.name, "label": path.name, "case_dir": str(path)} for path in sorted((root / "cases").iterdir()) if path.is_dir()]

    rows = []
    for item in manifest:
        case_dir = root / str(item["case_dir"])
        depth, width, time_value = read_metrics(case_dir)
        rows.append(
            {
                "case_id": str(item["case_id"]),
                "label": str(item.get("label", item["case_id"])),
                "case_dir": case_dir,
                "depth": depth,
                "width": width,
                "time": time_value,
                "png": section_png(case_dir, time_value),
            }
        )

    n_cases = max(1, len(rows))
    fig = plt.figure(figsize=(max(14, 3.1 * n_cases), 8), constrained_layout=True)
    grid = fig.add_gridspec(2, n_cases, height_ratios=[1.05, 1.0])

    for index, row in enumerate(rows):
        ax_img = fig.add_subplot(grid[0, index])
        ax_img.set_title(row["label"], fontsize=10)
        ax_img.set_xticks([])
        ax_img.set_yticks([])
        if row["png"] is None:
            ax_img.text(0.5, 0.5, "No section PNG yet", ha="center", va="center", transform=ax_img.transAxes)
        else:
            image = mpimg.imread(row["png"])
            ax_img.imshow(xy_crop(image))

    labels = [row["label"] for row in rows]
    y_pos = list(range(len(rows)))
    depth_values = [row["depth"] if row["depth"] is not None else 0.0 for row in rows]
    width_values = [row["width"] if row["width"] is not None else 0.0 for row in rows]

    split = max(1, n_cases // 2)
    ax_depth = fig.add_subplot(grid[1, :split])
    ax_width = fig.add_subplot(grid[1, split:])

    for ax, values, title, xlabel, color in (
        (ax_depth, depth_values, "Keyhole depth", "Max depth (um)", "#2878B5"),
        (ax_width, width_values, "Keyhole width", "Width at max-depth section (um)", "#D55E00"),
    ):
        ax.barh(y_pos, values, color=color)
        ax.set_title(title, fontsize=12)
        ax.set_yticks(y_pos, labels=labels, fontsize=8)
        ax.invert_yaxis()
        ax.set_xlabel(xlabel)
        ax.grid(axis="x", alpha=0.25)
        xmax = max(values) * 1.15 if max(values, default=0.0) > 0 else 1.0
        ax.set_xlim(0, xmax)
        for y, value in zip(y_pos, values):
            if value > 0:
                ax.text(value, y, f" {value:.1f}", va="center", ha="left", fontsize=8)
            else:
                ax.text(0, y, " pending", va="center", ha="left", fontsize=8, color="0.45")

    fig.suptitle("Parametric sweep dashboard", fontsize=16, fontweight="bold")
    output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output, dpi=180)
    plt.close(fig)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--root", type=Path, default=ROOT)
    parser.add_argument("--output", type=Path, default=ROOT / "results" / "dashboard.png")
    args = parser.parse_args()
    build_dashboard(args.root.resolve(), args.output.resolve())
    print(f"Wrote {args.output}")


if __name__ == "__main__":
    main()
