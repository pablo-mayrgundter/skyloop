#!/usr/bin/env python3
"""Figure for cost-down.md: reads scenarios.csv, waterfall.csv, capex_categories.csv -> cost_down.png."""
from __future__ import annotations

import csv
from pathlib import Path

HERE = Path(__file__).resolve().parent
BLUE, ORANGE, AQUA, YELLOW = "#2a78d6", "#eb6834", "#1baf7a", "#eda100"
INK, INK2, GRID, SURF = "#0b0b0b", "#52514e", "#e4e3df", "#fcfcfb"


def rd(name):
    with open(HERE / name, newline="") as f:
        return list(csv.DictReader(f))


def make(*_):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    plt.rcParams.update({"font.size": 9, "axes.edgecolor": INK2, "axes.labelcolor": INK, "xtick.color": INK2,
                         "ytick.color": INK2, "axes.spines.top": False, "axes.spines.right": False,
                         "figure.facecolor": SURF, "axes.facecolor": SURF})
    fig, axs = plt.subplots(1, 3, figsize=(15, 5.2), gridspec_kw={"width_ratios": [1.0, 1.25, 1.0]})

    # (a) capex by category, stacked, M2.0 and M1.7 x scenario
    cats = ["Lid & envelope", "Tethers, anchors, stations", "H2, controls, shuttle, install",
            "Engineering, contingency, land"]
    cols = [BLUE, ORANGE, AQUA, YELLOW]
    cc = rd("capex_categories.csv")
    groups = [(d, s) for d in ("M2.0", "M1.7") for s in ("base", "plausible", "aggressive")]
    ax = axs[0]
    x = []
    for i, (d, s) in enumerate(groups):
        xi = i + (0.6 if d == "M1.7" else 0)
        x.append(xi)
        bottom = 0
        for c, col in zip(cats, cols):
            v = next(float(r["usd_m2"]) for r in cc if r["design"].startswith(d) and r["scenario"] == s and r["category"] == c)
            ax.bar(xi, v, 0.72, bottom=bottom, color=col, edgecolor=SURF, linewidth=1.5,
                   label=c if i == 0 else None)
            bottom += v
        ax.text(xi, bottom + 1.5, f"${bottom:.0f}", ha="center", va="bottom", color=INK, fontsize=8)
    ax.set_xticks(x, [f"{d}\n{s[:5]}." if s != "base" else f"{d}\nbase" for d, s in groups], fontsize=8)
    ax.set_ylabel("Capital cost, $/m² of planform")
    ax.set_title("(a) Capex by category", loc="left", fontsize=10, color=INK)
    ax.set_ylim(0, 118)
    ax.legend(frameon=False, fontsize=7.5, loc="upper right", ncol=1)
    ax.yaxis.grid(True, color=GRID, linewidth=0.6)
    ax.set_axisbelow(True)

    # (b) cumulative net annual cost by lever step, M2.0
    wf = rd("waterfall.csv")
    ax = axs[1]
    for level, col, off in (("plausible", BLUE, -0.19), ("aggressive", ORANGE, 0.19)):
        rows = [r for r in wf if r["design"].startswith("M2.0") and r["level"] == level]
        ys = list(range(len(rows)))[::-1]
        vals = [float(r["net_musd_km_yr"]) for r in rows]
        ax.barh([y + off for y in ys], vals, 0.36, color=col, edgecolor=SURF, linewidth=1.0, label=level)
        for y, v in zip(ys, vals):
            ax.text(v + 0.3, y + off, f"{v:.1f}", va="center", fontsize=7.5, color=INK2)
        labels = [r["step"].replace("Base (cost_model)", "Base (economics.md)") for r in rows]
    ax.set_yticks(list(range(len(labels)))[::-1], [lab.replace(": ", ":\n", 1) if len(lab) > 34 else lab for lab in labels], fontsize=8)
    ax.set_xlabel("Annualised cost net of dual-use, $M per km per year (M2.0, cumulative)")
    ax.set_title("(b) Cost-down roadmap, levers stacked in order", loc="left", fontsize=10, color=INK)
    ax.legend(frameon=False, fontsize=8, loc="lower right")
    ax.xaxis.grid(True, color=GRID, linewidth=0.6)
    ax.set_axisbelow(True)

    # (c) B/C at 100 % capture, spine, log scale
    sc = rd("scenarios.csv")
    jt = rd("joint_demand.csv") if (HERE / "joint_demand.csv").exists() else []
    ax = axs[2]
    xs = {"base": 0, "plausible": 1, "aggressive": 2}
    for d, mk, dx in (("M2.0", "o", -0.08), ("M1.7", "s", 0.08)):
        for vs, col in (("bl", BLUE), ("sub", ORANGE)):
            pts = [(xs[r["scenario"]] + dx, float(r[f"bcr_{vs}_100"])) for r in sc if r["design"].startswith(d)]
            ax.plot([p[0] for p in pts], [p[1] for p in pts], color=col, linewidth=2 if d == "M2.0" else 1.2,
                    linestyle="-" if d == "M2.0" else "--", marker=mk, markersize=6, markeredgecolor=SURF,
                    label=f"{d} vs {'Boomless' if vs == 'bl' else 'subsonic'}")
            if jt:
                jp = [(xs[r["cost_scenario"]] + dx, float(r[f"bcr_{vs}_100"])) for r in jt
                      if r["design"].startswith(d) and r["demand_scenario"] == "base + life-avg growth x1.4"]
                ax.scatter([p[0] for p in jp], [p[1] for p in jp], marker=mk, s=28, facecolors="none", edgecolors=col,
                           linewidths=1.0, label="+ induced demand & growth (open)" if (d, vs) == ("M2.0", "bl") else None)
    ax.axhline(1.0, color=INK2, linewidth=1, linestyle=":")
    ax.text(0.0, 1.06, "B/C = 1", fontsize=8, color=INK2, ha="center", va="bottom")
    ax.set_yscale("log")
    ax.set_ylim(0.02, 5)
    ax.set_xticks([0, 1, 2], ["base", "plausible", "aggressive"])
    ax.set_ylabel("Spine B/C at 100 % capture (no fatigue cap)")
    ax.set_title("(c) Trunk B/C after cost-down", loc="left", fontsize=10, color=INK)
    ax.legend(frameon=False, fontsize=7, loc="upper left")
    ax.yaxis.grid(True, color=GRID, linewidth=0.6, which="major")
    ax.set_axisbelow(True)

    fig.tight_layout()
    fig.savefig(HERE / "cost_down.png", dpi=150)
    print("wrote cost_down.png")


if __name__ == "__main__":
    make()
