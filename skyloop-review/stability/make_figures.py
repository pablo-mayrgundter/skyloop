"""Figures for stability.md, from the CSVs written by stability.py:
stability_spectra.png (load and motion spectra vs allowable), stability_map.png (Mach x scheme),
stability_bigbore.png (bore trade)."""
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from matplotlib.colors import ListedColormap  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
INK, INK2, GRID, SURF = "#0b0b0b", "#52514e", "#e4e3df", "#fcfcfb"
PASS, MARG, FAIL = "#2a78d6", "#a9c8ef", "#ecebe7"
S1, S2, S3, S4 = "#2a78d6", "#eb6834", "#1baf7a", "#6250d6"
CMAP = ListedColormap([FAIL, MARG, PASS])
BIND = {"ride": "ride", "clearance": "stroke", "lift": "lift", "crit_speed": "v_crit", "power": "power",
        "support": "flight", "safety": "safety", "economics": "B/C"}


def style(ax, title):
    ax.set_facecolor(SURF)
    ax.set_title(title, loc="left", fontsize=10, color=INK, fontweight="bold")
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    for s in ("left", "bottom"):
        ax.spines[s].set_color(INK2)
    ax.tick_params(colors=INK2, labelsize=8)
    ax.xaxis.label.set_color(INK2)
    ax.yaxis.label.set_color(INK2)
    ax.grid(True, color=GRID, lw=0.6, which="major")
    ax.set_axisbelow(True)


def fig_spectra():
    g = pd.read_csv(os.path.join(HERE, "gust_spectrum.csv"))
    g = g[pd.to_numeric(g["lambda_lo_km"], errors="coerce").notna()].copy()
    for c in g.columns:
        g[c] = pd.to_numeric(g[c], errors="coerce")
    lc = np.sqrt(g["lambda_lo_km"] * g["lambda_hi_km"])
    mo = pd.read_csv(os.path.join(HERE, "motion_octaves.csv"))
    fig, axs = plt.subplots(1, 3, figsize=(15, 4.6), facecolor=SURF)
    ax = axs[0]
    style(ax, "(a) Gusts per octave of wavelength, 20 km")
    for lev, col in (("quiet", S3), ("typical", S1), ("active", S2)):
        ax.plot(lc, g[f"u_rms_{lev}_m_s"], color=col, lw=2, marker="o", ms=4, label=f"{lev}")
    ax.set_xscale("log")
    ax.set_xlabel("wavelength along corridor, km")
    ax.set_ylabel("u' rms per octave, m/s")
    ax.text(0.03, 0.55, "lateral load ~180 N/m per 1 m/s\n(187 m sail, U 18 m/s crosswind)", transform=ax.transAxes,
            fontsize=8, color=INK2)
    ax.legend(fontsize=8, frameon=False, loc="upper left", title="Lindborg fit x", title_fontsize=8)
    for i, M in enumerate((2.0, 4.0)):
        ax = axs[1 + i]
        style(ax, f"({'bc'[i]}) Tube motion per octave vs ride allowance, M{M:g}")
        cases = [("bare tethers, T 10 MN", S2), ("keel EI 1e14 (h 150)", S4), ("taut line T 0.3 GN", S3),
                 ("taut line T 1 GN", S1)]
        for name, col in cases:
            d = mo[(mo.case == name) & (mo.Mach == M)]
            ax.plot(d.lambda_km, d.amp_m, color=col, lw=2, marker="o", ms=4, label=name)
        d = mo[(mo.case == cases[0][0]) & (mo.Mach == M)]
        ax.plot(d.lambda_km, d.allow_m, color=INK, lw=2, ls="--", label="allowed (0.05 g, 1 Hz secondary)")
        ax.set_xscale("log")
        ax.set_yscale("log")
        ax.set_ylim(1e-4, 1e3)
        ax.set_xlabel("wavelength along corridor, km")
        ax.set_ylabel("displacement amplitude per octave, m")
        ax.legend(fontsize=7.5, frameon=False, loc="upper left")
    fig.tight_layout()
    fig.savefig(os.path.join(HERE, "stability_spectra.png"), dpi=140, facecolor=SURF)
    plt.close(fig)


def fig_map():
    df = pd.read_csv(os.path.join(HERE, "feasibility_map.csv"))
    schemes = list(dict.fromkeys(df.scheme))
    machs = sorted(df.Mach.unique())
    code = {"fail": 0, "marginal": 1, "pass": 2}
    fig, axs = plt.subplots(1, 2, figsize=(15, 7.0), facecolor=SURF, gridspec_kw=dict(width_ratios=[1.25, 1]))
    for ax, col, ttl in ((axs[0], "tech_status", "(a) Technical status (ride, stroke, lift, v_crit, power, flight, safety)"),
                         (axs[1], "full_status", "(b) Including B/C (base / low-cost proxy)")):
        Z = np.zeros((len(schemes), len(machs)))
        for i, s in enumerate(schemes):
            for j, M in enumerate(machs):
                r = df[(df.scheme == s) & (df.Mach == M)].iloc[0]
                Z[i, j] = code[r[col]]
                if col == "tech_status":
                    txt = BIND.get(r["binding"], r["binding"]) if r[col] != "pass" else ""
                else:
                    txt = f"{r['BC_base']:.2f}\n{r['BC_low']:.2f}"
                ax.text(j, i, txt, ha="center", va="center", fontsize=6.5, color=INK)
        ax.imshow(Z, cmap=CMAP, vmin=0, vmax=2, aspect="auto")
        ax.set_xticks(range(len(machs)))
        ax.set_xticklabels([f"{m:g}" for m in machs], fontsize=8, color=INK2)
        ax.set_yticks(range(len(schemes)))
        ax.set_yticklabels(schemes if ax is axs[0] else [], fontsize=8, color=INK)
        ax.set_xlabel("Mach (V = M x 295 m/s)", color=INK2)
        ax.set_title(ttl, loc="left", fontsize=10, color=INK, fontweight="bold")
        ax.set_xticks(np.arange(-.5, len(machs)), minor=True)
        ax.set_yticks(np.arange(-.5, len(schemes)), minor=True)
        ax.grid(which="minor", color=SURF, lw=2)
        ax.tick_params(which="minor", length=0)
        for s in ax.spines.values():
            s.set_visible(False)
    from matplotlib.patches import Patch
    axs[1].legend(handles=[Patch(color=PASS, label="pass"), Patch(color=MARG, label="marginal (within ~2x)"),
                           Patch(color=FAIL, label="fail")], fontsize=8, frameon=False,
                  loc="upper left", bbox_to_anchor=(0, -0.08), ncol=3)
    axs[0].text(0, -0.1, "cell text = binding constraint where not a pass", transform=axs[0].transAxes,
                fontsize=8, color=INK2)
    fig.tight_layout()
    fig.savefig(os.path.join(HERE, "stability_map.png"), dpi=140, facecolor=SURF)
    plt.close(fig)


def fig_bigbore():
    b = pd.read_csv(os.path.join(HERE, "bigbore.csv"))
    m = b[b["mode"] == "maglev, active guideway"]
    f = b[b["mode"] == "flying pod"]
    s = pd.read_csv(os.path.join(HERE, "sensitivity.csv"))
    s = s[s.study == "min bore"]
    fig, axs = plt.subplots(1, 3, figsize=(15, 4.4), facecolor=SURF)
    ax = axs[0]
    style(ax, "(a) Tube wall vs its own vacuum lift")
    x = np.arange(len(m))
    ax.bar(x - 0.2, m.wall_kg_m, 0.38, color=S2, label="wall (lattice, vacuum-review closure; D 5 sandwich)")
    ax.bar(x + 0.2, m.vac_lift_kg_m, 0.38, color=S1, label="vacuum lift at 5.5 Pa")
    ax.set_xticks(x)
    ax.set_xticklabels([f"D {d:g}" for d in m.D_m])
    ax.set_ylabel("kg per m")
    ax.legend(fontsize=7.5, frameon=False, loc="upper left")
    ax = axs[1]
    style(ax, "(b) Envelope radius and sail with the bore")
    ax.plot(m.D_m, m.env_R_m, color=S1, lw=2, marker="o", label="H2 envelope radius, m")
    ax.plot(m.D_m, m.sail_W_m, color=S2, lw=2, marker="o", label="total sail width (2R + D), m")
    ax.axhline(100, color=INK2, ls="--", lw=1)
    ax.text(m.D_m.min(), 103, "regime limit R 100 m", fontsize=7.5, color=INK2)
    ax.set_xlabel("bore D, m")
    ax.legend(fontsize=7.5, frameon=False, loc="upper left")
    ax = axs[2]
    style(ax, "(c) Smallest bore that passes ride + 4-sigma stroke")
    cols = {"bare tethers": S2, "taut 0.1 GN": S3, "taut 0.3 GN": S1}
    for (base, kind), d in s.groupby(["base", "kind"]):
        ls = "-" if kind == "actuated guideway" else "--"
        ax.plot(d.Mach, d.D_min_m, color=cols[base], lw=2, ls=ls, marker="o", ms=4, label=f"{base}, {kind}")
    ax.set_xlabel("Mach")
    ax.set_ylabel("minimum bore D, m")
    ax.legend(fontsize=7, frameon=False, loc="upper left")
    fig.tight_layout()
    fig.savefig(os.path.join(HERE, "stability_bigbore.png"), dpi=140, facecolor=SURF)
    plt.close(fig)


def main():
    for fn in (fig_spectra, fig_bigbore, fig_map):
        try:
            fn()
            print(f"figure {fn.__name__} ok")
        except FileNotFoundError as e:
            print(f"figure {fn.__name__} skipped: {e}")


if __name__ == "__main__":
    main()
