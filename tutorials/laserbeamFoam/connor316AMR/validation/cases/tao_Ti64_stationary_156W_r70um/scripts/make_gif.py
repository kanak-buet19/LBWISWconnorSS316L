#!/usr/bin/env python3
"""Create animated GIF from vtu_sections PNGs.

Usage:
    python scripts/make_gif.py                  # auto-detect case from script location
    python scripts/make_gif.py --case /path/to/case
    python scripts/make_gif.py --fps 10 --out custom.gif
"""

from __future__ import annotations

import argparse
from pathlib import Path

from PIL import Image


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--case", type=Path, default=Path(__file__).resolve().parents[1])
    parser.add_argument("--fps", type=float, default=5.0)
    parser.add_argument("--out", type=Path, default=None)
    args = parser.parse_args()

    case = args.case.resolve()
    png_dir = case / "post-processing-data" / "vtu_sections"
    out_gif = args.out or case / "post-processing-data" / "meltpool_sections.gif"

    frames = sorted(png_dir.glob("*.png"))
    if not frames:
        raise SystemExit(f"No PNGs found in {png_dir}")

    print(f"Loading {len(frames)} frames from {png_dir}")
    imgs = [Image.open(f).convert("RGBA") for f in frames]

    duration_ms = int(1000 / args.fps)
    imgs[0].save(
        out_gif,
        save_all=True,
        append_images=imgs[1:],
        duration=duration_ms,
        loop=0,
        optimize=False,
    )
    print(f"Wrote {out_gif}  ({out_gif.stat().st_size / 1e6:.1f} MB, {args.fps} fps)")


if __name__ == "__main__":
    main()
