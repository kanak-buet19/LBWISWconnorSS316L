"""
Beam profile: 2D irradiance heatmap + H/V cross-sections.
Run from repo root:  python scripts/plot_beam_profile.py [--power WATTS]
Default power: 201.3 W (stationary experiment value).
Note: profile measured at beam waist; sample surface was 2.8 mm below → 122.5 μm 1/e² at sample.
"""

import argparse
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
from pathlib import Path

ROOT = Path(__file__).parent.parent
PIXEL_UM = 5.5
N = 60

plt.rcParams.update({
    "font.family":      "DejaVu Sans",
    "font.size":        11,
    "axes.titlesize":   12,
    "axes.labelsize":   11,
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

parser = argparse.ArgumentParser()
parser.add_argument("--power", type=float, default=201.3,
                    help="Laser power in Watts (default: 201.3)")
args = parser.parse_args()

data = np.genfromtxt(
    ROOT / "beam_profile" / "beam_profile_5p5um_normalized.csv", delimiter=","
)
irr = data * args.power   # W/pixel  (integral ≈ laser power)
x   = np.arange(N) * PIXEL_UM   # μm

# Peak location
pr, pc = np.unravel_index(np.argmax(irr), irr.shape)
peak_x, peak_y = pc * PIXEL_UM, pr * PIXEL_UM

# 1/e² radius from horizontal slice through peak
threshold = irr[pr, pc] / np.e**2
h_slice   = irr[pr, :]
crossings = np.where(h_slice < threshold)[0]
crossings = crossings[crossings > pc]   # right side
r_1e2_um  = (crossings[0] - pc) * PIXEL_UM if len(crossings) else N * PIXEL_UM / 2

fig = plt.figure(figsize=(14, 5), constrained_layout=True)
fig.suptitle(
    f"Laser Beam Profile at Waist  (P = {args.power} W, pixel = {PIXEL_UM} μm)\n"
    f"Sample surface 2.8 mm below waist → calculated spot diameter at sample = 122.5 μm (1/e²)",
    fontsize=12, fontweight="bold",
)

axes = fig.subplot_mosaic([["heatmap", "horiz", "vert"]],
                          width_ratios=[1.1, 1, 1])

# ── 2D heatmap ────────────────────────────────────────────────────────────────
ax = axes["heatmap"]
extent = [0, N * PIXEL_UM, 0, N * PIXEL_UM]
im = ax.imshow(irr, origin="lower", extent=extent, cmap="inferno", aspect="equal")
cb = plt.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
cb.set_label("Irradiance (W per pixel area)", fontsize=10)

circle = plt.Circle((peak_x, peak_y), r_1e2_um, color="white",
                     fill=False, linewidth=1.8, linestyle="--",
                     label=f"1/e² radius = {r_1e2_um:.1f} μm")
ax.add_patch(circle)
ax.axhline(peak_y, color="cyan", linewidth=0.9, alpha=0.6, linestyle=":")
ax.axvline(peak_x, color="cyan", linewidth=0.9, alpha=0.6, linestyle=":")
ax.set_xlabel("x (μm)")
ax.set_ylabel("y (μm)")
ax.set_title("2D Irradiance Map")
ax.legend(loc="upper right", fontsize=8, framealpha=0.7)
ax.grid(False)

# ── Horizontal cross-section ──────────────────────────────────────────────────
ax = axes["horiz"]
ax.plot(x, irr[pr, :], color="#2166ac", linewidth=1.5, label="Horizontal slice")
ax.axhline(threshold, color="#555", linestyle="--", linewidth=1.1, label="1/e² threshold")
ax.axvline(peak_x, color="gray", linestyle=":", linewidth=0.9)
ax.axvline(peak_x + r_1e2_um, color="#aaa", linestyle=":", linewidth=0.9)
ax.axvline(peak_x - r_1e2_um, color="#aaa", linestyle=":", linewidth=0.9)
ax.fill_between(x, irr[pr, :], 0, alpha=0.15, color="#2166ac")
ax.set_xlabel("x (μm)")
ax.set_ylabel("Irradiance (W per pixel area)")
ax.set_title("Horizontal Cross-section Through Peak")
ax.legend(framealpha=0.9)
ax.yaxis.set_minor_locator(ticker.AutoMinorLocator())

# ── Vertical cross-section ────────────────────────────────────────────────────
ax = axes["vert"]
ax.plot(x, irr[:, pc], color="#d6604d", linewidth=1.5, label="Vertical slice")
ax.axhline(threshold, color="#555", linestyle="--", linewidth=1.1, label="1/e² threshold")
ax.axvline(peak_y, color="gray", linestyle=":", linewidth=0.9)
ax.axvline(peak_y + r_1e2_um, color="#aaa", linestyle=":", linewidth=0.9)
ax.axvline(peak_y - r_1e2_um, color="#aaa", linestyle=":", linewidth=0.9)
ax.fill_between(x, irr[:, pc], 0, alpha=0.15, color="#d6604d")
ax.set_xlabel("y (μm)")
ax.set_ylabel("Irradiance (W per pixel area)")
ax.set_title("Vertical Cross-section Through Peak")
ax.legend(framealpha=0.9)
ax.yaxis.set_minor_locator(ticker.AutoMinorLocator())

print(f"1/e² radius at waist : {r_1e2_um:.1f} μm  (diameter {2*r_1e2_um:.1f} μm)")
print(f"Integral check       : {irr.sum():.2f} W  (expected ≈ {args.power} W)")
for ext in ("png", "pdf"):
    out = ROOT / "scripts" / f"beam_profile.{ext}"
    plt.savefig(out, dpi=150, bbox_inches="tight")
    print(f"Saved → {out}")
plt.show()
