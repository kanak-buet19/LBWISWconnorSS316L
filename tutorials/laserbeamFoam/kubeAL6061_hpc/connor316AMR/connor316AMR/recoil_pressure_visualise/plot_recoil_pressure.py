#!/usr/bin/env python3
"""
Plot recoil pressure models from ambient-pressure-dependent boiling point.

Models:
    condensation:
        PR = 0.5*(1 + beta)*Pamb*exp(Lv_mol*(T - Tb_amb)/(R*T*Tb_amb))

    modified vacuum:
        PR = condensation + P0 - Pamb

Below Tb_amb, recoil pressure is not plotted/activated.
"""

from pathlib import Path
import math

import matplotlib.pyplot as plt
import numpy as np


# --------------------------------------------------------------------
# Constants (AH36 steel values from paper)
# --------------------------------------------------------------------
P0 = 101325.0      # Pa, reference pressure (1 atm)
LV = 1.053e7       # J/kg, latent heat of vaporization
MM = 0.027         # kg/mol, molar mass
R = 8.314          # J/mol/K, universal gas constant
TB_REF = 2792.0    # K, boiling temperature at P0

LV_MOL = LV*MM     # J/mol


# --------------------------------------------------------------------
# Cases: (ambient pressure [Pa], condensation coefficient beta)
# --------------------------------------------------------------------
CASES = [
    (1e1, 0.8),
    (1e3, 0.4),
    (1e5, 0.18),
]

T_MAX = 3200.0
N_T = 400
STAT_TEMPERATURES = [2400.0, 2600.0, 2800.0, 3000.0, 3200.0]


def tvap(P2, P1=P0, T1=TB_REF, dHvap=LV_MOL):
    """Boiling temperature at pressure P2 from Clausius-Clapeyron."""
    ln_term = math.log(P2/P1)
    inv_T2 = (1.0/T1) - (R/dHvap)*ln_term
    if inv_T2 <= 0:
        raise ValueError("Invalid input: non-physical T2")
    return 1.0/inv_T2


def pr_condensation(T, Pamb, beta, Tb_amb):
    """Condensation-coefficient recoil pressure [Pa]."""
    return 0.5*(1.0 + beta)*Pamb*np.exp(
        LV_MOL*(T - Tb_amb)/(R*T*Tb_amb)
    )


def pr_modified(T, Pamb, beta, Tb_amb):
    """Modified reduced-pressure recoil pressure [Pa]."""
    return pr_condensation(T, Pamb, beta, Tb_amb) + P0 - Pamb


def active_value(model_func, T, Pamb, beta, Tb_amb):
    """Return zero below ambient boiling temperature."""
    if T < Tb_amb:
        return 0.0
    return model_func(np.array([T]), Pamb, beta, Tb_amb)[0]


def main():
    fig, ax = plt.subplots(figsize=(10, 6))
    stats_rows = []

    for Pamb, beta in CASES:
        Tb_amb = tvap(Pamb)
        T_range = np.linspace(Tb_amb, T_MAX, N_T)

        y_cond = pr_condensation(T_range, Pamb, beta, Tb_amb)/1e3
        y_mod = pr_modified(T_range, Pamb, beta, Tb_amb)/1e3

        line = ax.plot(
            T_range,
            y_cond,
            linewidth=2.5,
            label=f"Condensation: Pamb={Pamb:.0e} Pa, beta={beta}, Tb={Tb_amb:.1f} K",
        )
        color = line[0].get_color()
        ax.plot(
            T_range,
            y_mod,
            "--",
            linewidth=2.5,
            color=color,
            label=f"Modified: Pamb={Pamb:.0e} Pa",
        )

        for T_stat in STAT_TEMPERATURES:
            cond = active_value(pr_condensation, T_stat, Pamb, beta, Tb_amb)
            mod = active_value(pr_modified, T_stat, Pamb, beta, Tb_amb)
            stats_rows.append((Pamb, beta, Tb_amb, T_stat, cond, mod))

    ax.set_xlabel("Surface Temperature (K)", fontsize=14)
    ax.set_ylabel("Recoil Pressure PR (kPa)", fontsize=14)
    ax.set_title("Recoil Pressure Models vs Ambient Pressure", fontsize=15)
    ax.grid(True, alpha=0.3)
    ax.legend(fontsize=9)

    plt.tight_layout()

    output_path = Path(__file__).with_name("recoil_pressure_comparison.png")
    fig.savefig(output_path, dpi=300)

    print(f"Saved plot: {output_path}")
    print()
    print(
        "Pamb_Pa     beta   Tb_amb_K   T_K     "
        "PR_cond_kPa   PR_modified_kPa"
    )
    print("-"*72)
    for Pamb, beta, Tb_amb, T_stat, cond, mod in stats_rows:
        print(
            f"{Pamb:9.3e}  "
            f"{beta:5.2f}  "
            f"{Tb_amb:9.2f}  "
            f"{T_stat:6.0f}  "
            f"{cond/1e3:11.4f}  "
            f"{mod/1e3:15.4f}"
        )

    plt.show()


if __name__ == "__main__":
    main()
