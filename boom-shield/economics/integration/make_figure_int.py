#!/usr/bin/env python3
"""Figure for integration.md -> integration.png (called by integration.py)."""
from __future__ import annotations

import math

BLUE, ORANGE, AQUA, YELLOW = "#2a78d6", "#eb6834", "#1baf7a", "#eda100"
INK, INK2, GRID, SURF = "#0b0b0b", "#52514e", "#e4e3df", "#fcfcfb"
M20 = "M2.0 (Concorde-class, psf baseline)"
CATS = {
    "Lid & envelope": ["Lid mass (PE sheet/laminate)", "Lid UV/weathering skin", "Load tapes / scrim (UHMWPE)",
                       "Floor film", "Partitions and bellows", "Soft ties", "Seams and fabrication"],
    "Tethers, anchors, stations": ["Tethers", "Anchors", "Winch/tension stations"],
    "H2, controls, shuttle, install": ["H2 first fill", "Controls, sensors, transponders", "Diurnal gas shuttle",
                                       "Launch / installation"],
}


def make(here, results, bcr_rows, learn_rows, plant_rows):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    plt.rcParams.update({"font.size": 9, "axes.edgecolor": INK2, "axes.labelcolor": INK, "xtick.color": INK2,
                         "ytick.color": INK2, "axes.spines.top": False, "axes.spines.right": False,
                         "figure.facecolor": SURF, "axes.facecolor": SURF})
    fig, axs = plt.subplots(1, 3, figsize=(15.5, 5.4), gridspec_kw={"width_ratios": [1.1, 1.0, 1.1]})

    # (a) $/m2 by category, M2.0, all with the second tether row
    cases = [("cost-down plausible + 2nd tether row", "Cost-down\nplausible", None),
             ("cost-down aggressive + 2nd tether row", "Cost-down\naggressive", None),
             ("integrated, plausible levers, central", "Integrated,\nplausible levers", "plausible"),
             ("integrated, aggressive levers, central", "Integrated,\naggressive levers", "aggressive")]
    cols = [BLUE, ORANGE, AQUA, YELLOW]
    ax = axs[0]
    for i, (key, lab, lev) in enumerate(cases):
        r = results[(M20, key)]
        vals = [sum(r["lines"][k] for k in ks) for ks in CATS.values()] + [r["eng"] + r["cont"] + r["land"]]
        bottom = 0.0
        for j, (v, c) in enumerate(zip(vals, cols)):
            ax.bar(i, v - 0.12, bottom=bottom + 0.06, width=0.62, color=c, linewidth=0,
                   label=(list(CATS) + ["Engineering, contingency, land"])[j] if i == 0 else None)
            bottom += v
        tot = r["usd_m2"]
        if lev:
            lo = results[(M20, f"integrated, {lev} levers, low")]["usd_m2"]
            hi = results[(M20, f"integrated, {lev} levers, high")]["usd_m2"]
            ax.plot([i, i], [lo, hi], color=INK2, lw=1.2)
            ax.plot([i - 0.08, i + 0.08], [lo, lo], color=INK2, lw=1.2)
            ax.plot([i - 0.08, i + 0.08], [hi, hi], color=INK2, lw=1.2)
            ax.text(i + 0.12, hi, f"{hi:.0f}", color=INK2, fontsize=8, va="center")
            ax.text(i + 0.12, lo, f"{lo:.0f}", color=INK2, fontsize=8, va="center")
        if lev:
            ax.text(i - 0.33, tot + 0.5, f"${tot:.1f}", ha="right", va="bottom", color=INK, fontsize=9,
                    fontweight="bold")
        else:
            ax.text(i, tot + 0.6, f"${tot:.1f}", ha="center", va="bottom", color=INK, fontsize=9, fontweight="bold")
    ax.axhline(19.8, color=INK2, lw=1, ls="--", label="Cost-down aggressive as published (1 row): $19.8")
    ax.set_xticks(range(len(cases)))
    ax.set_xticklabels([c[1] for c in cases], fontsize=8)
    ax.set_ylabel("Capital cost, $/m² (M2.0, two tether rows)")
    ax.set_ylim(0, 62)
    ax.yaxis.grid(True, color=GRID, lw=0.8)
    ax.set_axisbelow(True)
    ax.legend(frameon=False, fontsize=8, loc="upper right")
    ax.set_title("(a) Capital cost per m²; whiskers = integrated low–high", fontsize=10, loc="left", color=INK)

    # (b) Wright's law: unit cost relative to the pilot unit (end of a 1 % pilot) vs cumulative units
    ax = axs[1]
    items = [("modules", 0.15, 14456, BLUE), ("tether stations", 0.12, 1446, ORANGE),
             ("control cells", 0.18, 101184, AQUA), ("seam-metres (×10⁻⁵)", 0.10, 1.83e9 / 1e5, YELLOW)]
    for name, lr, N, c in items:
        x0 = N * 0.01
        xs = [x0 * (N / x0) ** (k / 60) for k in range(61)]
        b = -math.log2(1 - lr)
        ys = [(x / x0) ** (-b) for x in xs]
        ax.plot(xs, ys, color=c, lw=2)
        ax.plot([xs[-1]], [ys[-1]], "o", ms=7, color=c, mec=SURF, mew=2)
        ax.text(xs[-1] * 1.25, ys[-1], f"{name}\nLR {lr:.0%}: last unit ×{ys[-1]:.2f}", fontsize=8, va="center",
                color=INK)
    ax.set_xscale("log")
    ax.set_xlim(8, 2e6)
    ax.set_ylim(0, 1.05)
    ax.set_xlabel("Cumulative units built on the spine (log)")
    ax.set_ylabel("Unit labour/overhead cost ÷ cost at end of 1 % pilot")
    ax.yaxis.grid(True, color=GRID, lw=0.8)
    ax.set_axisbelow(True)
    ax.set_title("(b) Wright's law on the rote units (central rates)", fontsize=10, loc="left", color=INK)

    # (c) B/C at 100 % capture, base demand, M2.0: spine and 6,989 km shoulder
    ax = axs[2]
    order = ["cost-down plausible + 2nd row", "cost-down aggressive + 2nd row",
             "integrated, plausible levers (central)", "integrated, aggressive levers (central)",
             "integrated, aggressive levers (low cost)"]
    labels = ["Cost-down plausible", "Cost-down aggressive", "Integrated, plausible levers",
              "Integrated, aggressive levers", "Integrated, aggr., low cost"]
    get = {(r["network"], r["cost_case"]): r for r in bcr_rows
           if r["design"] == M20 and r["demand"] == "base"}
    for i, case in enumerate(order):
        y = len(order) - 1 - i
        for net, mk, fill in (("spine", "o", True), ("shoulder 6,989 km", "o", False)):
            r = get.get((net, case))
            if not r:
                continue
            for key, c, dy in (("bcr_sub_100", BLUE, 0.12), ("bcr_bl_100", ORANGE, -0.12)):
                v = float(r[key])
                ax.plot([v], [y + dy], mk, ms=8, color=c if fill else SURF, mec=c, mew=2)
                if fill:
                    ax.text(v * 1.08, y + dy, f"{v:.2f}", fontsize=8, va="center", color=INK)
    ax.axvline(1.0, color=INK2, lw=1, ls="--")
    ax.set_xscale("log")
    ax.set_xlim(0.05, 5)
    ax.set_ylim(-0.5, len(order) + 0.5)
    ax.set_yticks(range(len(order)))
    ax.set_yticklabels(list(reversed(labels)), fontsize=8)
    ax.set_xlabel("B/C at 100 % capture, base demand (log)\n10 % capture: divide by 10")
    ax.xaxis.grid(True, color=GRID, lw=0.8)
    ax.set_axisbelow(True)
    from matplotlib.lines import Line2D
    hs = [Line2D([], [], marker="o", ls="", color=BLUE, ms=8, label="vs subsonic"),
          Line2D([], [], marker="o", ls="", color=ORANGE, ms=8, label="vs Boomless"),
          Line2D([], [], marker="o", ls="", color=SURF, mec=INK2, mew=2, ms=8, label="6,989 km network (open)")]
    ax.legend(handles=hs, frameon=False, fontsize=8, loc="upper left", bbox_to_anchor=(0.0, 1.0))
    ax.set_title("(c) M2.0 B/C, spine (filled) and shoulder (open)", fontsize=10, loc="left", color=INK)

    fig.tight_layout()
    fig.savefig(here / "integration.png", dpi=130, facecolor=SURF)
    plt.close(fig)
