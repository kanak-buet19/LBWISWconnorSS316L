"""
Browse X-ray frames synced to absorption timeseries.

Usage:
  # Single frame
  python scripts/plot_xray_sync.py --exp stationary --frame 90

  # Interactive browser (arrow keys or slider)
  python scripts/plot_xray_sync.py --exp moving --browse

  # Save frame as PNG
  python scripts/plot_xray_sync.py --exp stationary --frame 90 --save

--exp:  stationary | moving
--zip:  raw | captioned | absorption  (default: raw)
"""

import argparse
import zipfile
import io
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import matplotlib.ticker as ticker
from matplotlib.widgets import Slider
from pathlib import Path

ROOT = Path(__file__).parent.parent

EXP_META = {
    "stationary": {"dir": ROOT / "stationary_laser", "label": "Stationary laser (spot, 2 ms pulse)"},
    "moving":     {"dir": ROOT / "moving_laser",     "label": "Moving laser (scan, 700 mm/s)"},
}
ZIP_KEYS = {
    "raw":        "xray_raw.zip",
    "captioned":  "xray_processed_captioned.zip",
    "absorption": "xray_processed_with_absorption.zip",
}

plt.rcParams.update({
    "font.family":      "DejaVu Sans",
    "font.size":        10,
    "axes.titlesize":   11,
    "axes.labelsize":   10,
    "xtick.labelsize":  9,
    "ytick.labelsize":  9,
    "legend.fontsize":  9,
    "axes.spines.top":  False,
    "axes.spines.right":False,
    "axes.grid":        True,
    "grid.alpha":       0.35,
    "grid.linewidth":   0.6,
    "figure.dpi":       150,
})


def load_absorption(exp_dir):
    df = pd.read_csv(exp_dir / "absorption_data.csv", index_col=0, low_memory=False)
    df["AbsAbsorptionUncertainty"] = pd.to_numeric(
        df["AbsAbsorptionUncertainty"], errors="coerce"
    )
    df["Time_ms"] = df["Time"] * 1000
    return df


def build_frame_index(df):
    """Map frame_number → CSV row index at the start of that frame."""
    index = {}
    for n in range(1, int(df["FrameNumber"].max()) + 1):
        rows = df.loc[df["FrameNumber"] == n]
        if not rows.empty:
            index[n] = rows.index[0]
    return index


def load_image_from_zip(zf, frame_num):
    """Find .tif whose filename ends with zero-padded frame_num."""
    from PIL import Image
    target = str(frame_num).zfill(3)
    for name in zf.namelist():
        stem = name.rsplit(".", 1)[0]
        if stem.endswith(target) and name.lower().endswith(".tif"):
            with zf.open(name) as f:
                return np.array(Image.open(io.BytesIO(f.read())))
    return None


def render(df, frame_index, zf, frame_num, ax_img, ax_rel, ax_abs, exp_label, zip_type):
    for ax in (ax_img, ax_rel, ax_abs):
        ax.cla()

    # ── X-ray image ───────────────────────────────────────────────────────────
    img = load_image_from_zip(zf, frame_num)
    if img is not None:
        vmin, vmax = np.percentile(img, [2, 98])
        ax_img.imshow(img, cmap="gray", vmin=vmin, vmax=vmax)
        ax_img.set_title(f"X-ray Frame {frame_num}  [{zip_type}]", fontsize=10)
    else:
        ax_img.text(0.5, 0.5, f"Frame {frame_num} not found in zip",
                    ha="center", va="center", transform=ax_img.transAxes, color="red")
    ax_img.axis("off")

    # Base data (laser on)
    laser_on = df["InputLaser"] > 0
    sub = df[laser_on]
    t = sub["Time_ms"]

    # Marker at current frame
    t_frame = r_frame = abs_frame = None
    if frame_num in frame_index:
        idx      = frame_index[frame_num]
        t_frame  = df.loc[idx, "Time_ms"]
        r_frame  = df.loc[idx, "RelativeAbsorption"]
        abs_frame = df.loc[idx, "AbsoluteAbsorption"]

    # ── Relative absorptance ──────────────────────────────────────────────────
    ax_rel.plot(t, sub["RelativeAbsorption"], color="#2166ac", linewidth=0.9,
                label="Relative absorptance")
    if t_frame is not None:
        ax_rel.axvline(t_frame, color="#d6604d", linewidth=1.4, linestyle="--", alpha=0.9)
        ax_rel.plot(t_frame, r_frame, "o", color="#d6604d", markersize=6,
                    label=f"Frame {frame_num}:  {r_frame:.1f}%  at  {t_frame:.4f} ms")
    ax_rel.set_xlabel("Time (ms)")
    ax_rel.set_ylabel("Relative Absorptance (%)")
    ax_rel.set_title("Relative Absorptance vs. Time")
    ax_rel.set_ylim(0, 110)
    ax_rel.yaxis.set_major_locator(ticker.MultipleLocator(20))
    ax_rel.legend(loc="upper left", framealpha=0.9)

    # ── Absolute absorbed power ───────────────────────────────────────────────
    a = sub["AbsoluteAbsorption"]
    u = sub["AbsAbsorptionUncertainty"]
    ax_abs.plot(t, a, color="#4d9a6b", linewidth=0.9, label="Absorbed power (W)")
    if u.notna().any():
        ax_abs.fill_between(t, a - u, a + u, color="#4d9a6b", alpha=0.2,
                            label="± expanded uncertainty")
    if t_frame is not None:
        ax_abs.axvline(t_frame, color="#d6604d", linewidth=1.4, linestyle="--", alpha=0.9)
        ax_abs.plot(t_frame, abs_frame, "o", color="#d6604d", markersize=6,
                    label=f"Frame {frame_num}:  {abs_frame:.1f} W")
    ax_abs.set_xlabel("Time (ms)")
    ax_abs.set_ylabel("Absorbed Power (W)")
    ax_abs.set_title("Absolute Absorbed Power vs. Time")
    ax_abs.legend(loc="upper left", framealpha=0.9)
    ax_abs.yaxis.set_minor_locator(ticker.AutoMinorLocator())


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--exp",    default="stationary", choices=["stationary", "moving"])
    parser.add_argument("--zip",    default="raw",        choices=["raw", "captioned", "absorption"])
    parser.add_argument("--frame",  type=int,             default=None)
    parser.add_argument("--browse", action="store_true",  help="Interactive frame browser")
    parser.add_argument("--save",   action="store_true",  help="Save PNG instead of showing")
    args = parser.parse_args()

    meta     = EXP_META[args.exp]
    exp_dir  = meta["dir"]
    exp_label = meta["label"]
    df        = load_absorption(exp_dir)
    frame_index = build_frame_index(df)
    frames    = sorted(frame_index.keys())
    zf        = zipfile.ZipFile(exp_dir / ZIP_KEYS[args.zip])

    if args.browse:
        fig = plt.figure(figsize=(15, 7), constrained_layout=True)
        fig.suptitle(
            f"NIST AMB2022 — {exp_label}  |  Arrow keys or slider to navigate",
            fontsize=11, fontweight="bold",
        )
        gs   = gridspec.GridSpec(2, 3, figure=fig, height_ratios=[10, 0.7], hspace=0.5)
        ax_img = fig.add_subplot(gs[0, 0])
        ax_rel = fig.add_subplot(gs[0, 1])
        ax_abs = fig.add_subplot(gs[0, 2])
        ax_sl  = fig.add_subplot(gs[1, :])

        slider = Slider(ax_sl, "Frame #", frames[0], frames[-1],
                        valinit=frames[0], valstep=1, color="#2166ac")

        def update(val):
            fn = int(slider.val)
            if fn not in frame_index:
                return
            render(df, frame_index, zf, fn, ax_img, ax_rel, ax_abs, exp_label, args.zip)
            fig.canvas.draw_idle()

        slider.on_changed(update)

        def on_key(event):
            fn = int(slider.val)
            if event.key == "right" and fn < frames[-1]:
                slider.set_val(fn + 1)
            elif event.key == "left" and fn > frames[0]:
                slider.set_val(fn - 1)

        fig.canvas.mpl_connect("key_press_event", on_key)
        render(df, frame_index, zf, frames[0], ax_img, ax_rel, ax_abs, exp_label, args.zip)
        plt.show()

    else:
        frame_num = args.frame if args.frame is not None else frames[len(frames) // 2]
        fig = plt.figure(figsize=(15, 5), constrained_layout=True)
        fig.suptitle(
            f"NIST AMB2022 — {exp_label}  |  X-ray sync with absorptance",
            fontsize=12, fontweight="bold",
        )
        gs     = gridspec.GridSpec(1, 3, figure=fig, width_ratios=[1, 1.1, 1.1])
        ax_img = fig.add_subplot(gs[0])
        ax_rel = fig.add_subplot(gs[1])
        ax_abs = fig.add_subplot(gs[2])
        render(df, frame_index, zf, frame_num, ax_img, ax_rel, ax_abs, exp_label, args.zip)

        if args.save:
            for ext in ("png", "pdf"):
                out = ROOT / "scripts" / f"xray_sync_{args.exp}_frame{frame_num:03d}.{ext}"
                plt.savefig(out, dpi=150, bbox_inches="tight")
                print(f"Saved → {out}")
        else:
            plt.show()

    zf.close()


if __name__ == "__main__":
    main()
