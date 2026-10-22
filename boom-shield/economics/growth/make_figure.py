#!/usr/bin/env python3
"""growth.png from the CSVs written by growth_model.py."""
from __future__ import annotations

import csv
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = Path(__file__).resolve().parent
BLUE, ORANGE, AQUA = "#2a78d6", "#eb6834", "#1baf7a"
INK, INK2, GRID, SURF = "#0b0b0b", "#52514e", "#e4e3df", "#fcfcfb"

# Approximate US domestic enplanements, M [K] (BTS search summaries G7/G8 + background knowledge);
# 2024-25 are T-100 domestic segment totals per coverage.md.
HIST = [(1980, 275), (1990, 424), (2000, 640), (2001, 596), (2007, 679), (2008, 652), (2009, 618),
        (2017, 742), (2019, 811), (2020, 340), (2024, 856), (2025, 845)]


def rd(name):
    with open(HERE / name, newline="") as f:
        return list(csv.DictReader(f))


def main():
    plt.rcParams.update({"font.size": 9, "axes.edgecolor": INK2, "axes.labelcolor": INK, "xtick.color": INK2,
                         "ytick.color": INK2, "axes.spines.top": False, "axes.spines.right": False,
                         "figure.facecolor": SURF, "axes.facecolor": SURF, "axes.titlesize": 10,
                         "axes.titleweight": "bold", "axes.titlelocation": "left"})
    fig, axs = plt.subplots(2, 2, figsize=(13, 9.2))
    traj = rd("trajectories.csv")
    yrs = [int(r["year"]) for r in traj]

    # (a) traffic index
    ax = axs[0][0]
    ax.plot([y for y, _ in HIST], [v / 845 for _, v in HIST], color=INK2, lw=1.5, marker="o", ms=3,
            label="History: domestic enplanements (approx.)")
    for case, col in (("low", ORANGE), ("base", BLUE), ("high", AQUA)):
        ax.plot(yrs, [float(r[f"traffic_{case}"]) for r in traj], color=col, lw=2, label=f"> 1,500 km O&D, {case}")
    for m in (1.5, 2, 3, 5, 10):
        ax.axhline(m, color=GRID, lw=0.8, zorder=0)
        ax.text(2094, m * 1.03, f"{m}×", color=INK2, fontsize=8)
    ax.set_yscale("log")
    ax.set_yticks([0.3, 0.5, 1, 2, 3, 5, 10])
    ax.set_yticklabels(["0.3", "0.5", "1", "2", "3", "5", "10"])
    ax.set_xlim(1978, 2100)
    ax.set_ylim(0.28, 12)
    ax.set_ylabel("Traffic index (2025 = 1)")
    ax.set_title("(a) US long-haul traffic: 10× is not in sight this century")
    ax.legend(frameon=False, fontsize=8, loc="upper left")

    grid = rd("roi_grid.csv")

    def series(vs, cost, col, world="base", capture="1.0", design="M2.0", net="spine", energy="central"):
        rows = [g for g in grid if g["design"] == design and g["network"] == net and g["demand"] == "base"
                and g["vs"] == vs and g["world"] == world and g["energy"] == energy and g["cost"] == cost
                and g["capture"] == capture]
        rows.sort(key=lambda g: int(g["open_year"]))
        return [int(g["open_year"]) for g in rows], [float(g[col]) for g in rows]

    for k, (vs, title) in enumerate((("bl", "(b) Spine M2.0 vs Boomless, 100 % capture"),
                                     ("sub", "(c) Spine M2.0 vs subsonic, 100 % capture"))):
        ax = axs[0][1] if k == 0 else axs[1][0]
        for cost, col, lab in (("integrated", BLUE, "Integrated $14/m²"), ("aggressive", ORANGE, "Aggressive $23/m²"),
                               ("base", AQUA, "Base $106/m²")):
            x, g = series(vs, cost, "bc_gross")
            _, n = series(vs, cost, "bc_net")
            ax.plot(x, g, color=col, lw=2, label=f"{lab}, gross")
            ax.plot(x, n, color=col, lw=2, ls="--", label=f"{lab}, net of fuel + CO₂")
        xl, gl = series(vs, "integrated", "bc_gross", world="low")
        xh, gh = series(vs, "integrated", "bc_gross", world="high")
        ax.fill_between(xl, gl, gh, color=BLUE, alpha=0.12, lw=0, label="Integrated gross, low–high world")
        ax.axhline(1, color=INK, lw=1)
        ax.axhline(0, color=INK2, lw=0.6)
        ax.set_xlim(2030, 2100)
        ax.set_xlabel("Opening year (build starts 10 yr earlier)")
        ax.set_ylabel("Life-cycle B/C at 7 %")
        ax.set_title(title)
        if vs == "bl":
            ax.set_ylim(-0.1, 3.0)
            ax.legend(frameon=False, fontsize=7.5, loc="upper left")
        else:
            ax.set_ylim(-6, 8)
            ax.text(2032, 7.2, "Key as in (b); band = integrated gross, low–high world", fontsize=8, color=INK2)
        ax.grid(axis="y", color=GRID, lw=0.6)

    # (d) flights/day
    ax = axs[1][1]
    base_f = 4716
    for capt, col in ((1.0, BLUE), (0.3, ORANGE)):
        for case, ls in (("base", "-"),):
            ax.plot(yrs, [base_f * capt * float(r[f"traffic_{case}"]) for r in traj], color=col, lw=2,
                    label=f"Spine busiest km, {int(capt * 100)} % capture, base traffic")
        ax.fill_between(yrs, [base_f * capt * float(r["traffic_low"]) for r in traj],
                        [base_f * capt * float(r["traffic_high"]) for r in traj], color=col, alpha=0.12, lw=0)
    ax.plot(yrs, [6440 * float(r["traffic_base"]) for r in traj], color=AQUA, lw=2,
            label="6,989 km shoulder, 100 %, base traffic")
    for cap, lab in ((10000, "Lid cap, base 10,000/day"), (1000, "Lid cap, low 1,000/day")):
        ax.axhline(cap, color=INK, lw=1, ls=":")
        ax.text(2031, cap * 1.07, lab, fontsize=8, color=INK)
    ax.set_yscale("log")
    ax.set_xlim(2025, 2100)
    ax.set_ylim(500, 40000)
    ax.set_ylabel("Supersonic flights/day (65 seats)")
    ax.set_title("(d) Fatigue cap binds ~2045–2060 at full capture")
    ax.legend(frameon=False, fontsize=8, loc="upper left")
    ax.grid(axis="y", color=GRID, lw=0.6)

    fig.tight_layout()
    fig.savefig(HERE / "growth.png", dpi=130)
    print("wrote growth.png")


if __name__ == "__main__":
    main()
