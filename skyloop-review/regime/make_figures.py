"""Figures for regime.md: reads regime_grid.csv / regime_pods.csv, writes regime_maps.png,
regime_pods.png and regime_glide.png."""
import math
import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from matplotlib.colors import ListedColormap  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.dont_write_bytecode = True
import regime as R  # noqa: E402

INK, INK2, GRID, SURF = "#0b0b0b", "#52514e", "#e4e3df", "#fcfcfb"
PASS, MARG, FAIL = "#2a78d6", "#a9c8ef", "#ecebe7"
ECON = "#eb6834"
CMAP = ListedColormap([FAIL, MARG, PASS])
LAB = {"support": "support\n(lift)", "drag_energy": "drag /\ncrossovers", "heating": "heating",
       "vacuum": "vacuum\n(leaks)", "structure": "structure", "alignment": "alignment", "route": "curve\nradius",
       "safety": "breach\nsafety", "capacity": "capacity", "economics": "B/C"}
SET_T = {"floating": "Floating, 20 km", "elevated": "Elevated on pylons", "buried": "Buried / tunnel"}


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


def edges(v, log=False):
    v = np.asarray(v, float)
    if log:
        lv = np.log10(v)
        e = np.concatenate([[lv[0] - (lv[1] - lv[0]) / 2], (lv[1:] + lv[:-1]) / 2, [lv[-1] + (lv[-1] - lv[-2]) / 2]])
        return 10 ** e
    return np.concatenate([[v[0] - (v[1] - v[0]) / 2], (v[1:] + v[:-1]) / 2, [v[-1] + (v[-1] - v[-2]) / 2]])


def panel(ax, sub, xcol, ycol, scol, bcol, logy, econ_col=None, label_marginal=True):
    xs = np.sort(sub[xcol].unique())
    ys = np.sort(sub[ycol].unique())
    Z = np.zeros((len(ys), len(xs)))
    B = np.empty((len(ys), len(xs)), dtype=object)
    E = np.zeros((len(ys), len(xs)))
    idx = {(round(r[xcol], 6), round(r[ycol], 6)): r for _, r in sub.iterrows()}
    for j, y in enumerate(ys):
        for i, x in enumerate(xs):
            r = idx[(round(x, 6), round(y, 6))]
            Z[j, i] = {"fail": 0, "marginal": 1, "pass": 2}[r[scol]]
            B[j, i] = r[bcol]
            if econ_col:
                E[j, i] = 1.0 if (r[econ_col] >= 1.0 and r[scol] != "fail") else 0.0
    ax.pcolormesh(edges(xs), edges(ys, logy), Z, cmap=CMAP, vmin=-0.5, vmax=2.5, shading="flat")
    # boundaries between binding constraints inside non-pass cells
    codes = {k: n for n, k in enumerate(sorted(set(B.ravel())))}
    C = np.vectorize(lambda b: codes[b])(B).astype(float)
    C[Z == 2] = np.nan
    if len(codes) > 1:
        ax.contour(xs, ys, np.nan_to_num(C, nan=-1), levels=np.arange(-0.5, len(codes)), colors=INK2,
                   linewidths=0.5, alpha=0.6)
    # label each binding region at its centroid (only sizeable regions)
    for b, n in codes.items():
        mask = (C == n) if label_marginal else ((C == n) & (Z == 0))
        if mask.sum() < 0.04 * mask.size:
            continue
        jj, ii = np.nonzero(mask)
        jm, im = int(np.median(jj)), int(np.median(ii))
        if not mask[jm, im]:
            k = np.argmin((jj - jj.mean()) ** 2 + (ii - ii.mean()) ** 2)
            jm, im = jj[k], ii[k]
        ax.text(xs[im], ys[jm], LAB.get(b, b), fontsize=7.2, color=INK, ha="center", va="center",
                bbox=dict(boxstyle="round,pad=0.15", fc=SURF, ec="none", alpha=0.8))
    if econ_col is not None and E.max() > 0:
        ax.contour(xs, ys, E, levels=[0.5], colors=ECON, linewidths=2.0)
    if logy:
        ax.set_yscale("log")


def fig_maps():
    d = pd.read_csv(os.path.join(HERE, "regime_grid.csv"))
    d = d[d.rule == "moving block 10 s"]
    fig, axs = plt.subplots(2, 3, figsize=(13.5, 8.2), sharex=True, sharey=True)
    fig.patch.set_facecolor(SURF)
    for row, sup in enumerate(["maglev", "glide"]):
        for col, sk in enumerate(R.SETTINGS):
            ax = axs[row, col]
            sub = d[(d.setting == sk) & (d.support == sup)]
            panel(ax, sub, "Mach", "p_Pa", "ops_status_bestbeta", "ops_binding_bestbeta", True, econ_col="BC_low")
            style(ax, f"{SET_T[sk]} · {sup}")
            pa = R.SETTINGS[sk]["p_amb"]
            ax.axhline(1e-3 * pa, color=INK, lw=1.0, ls=(0, (4, 3)))
            ax.text(7.95, 1e-3 * pa * 1.12, "10⁻³ of ambient", fontsize=7, color=INK, ha="right", va="bottom")
            sec = ax.secondary_yaxis("right", functions=(lambda p, pa=pa: p / pa, lambda f, pa=pa: f * pa))
            sec.tick_params(colors=INK2, labelsize=7.5)
            if col == 2:
                sec.set_ylabel("fraction of local ambient", color=INK2, fontsize=8)
            if col == 0:
                ax.set_ylabel("tube pressure, Pa", fontsize=9)
            if row == 1:
                ax.set_xlabel("Mach (V = M × 295 m/s)", fontsize=9)
    from matplotlib.patches import Patch
    from matplotlib.lines import Line2D
    h = [Patch(fc=PASS, label="all technical checks + capacity pass"),
         Patch(fc=MARG, label="marginal (worst check within ~2×)"),
         Patch(fc=FAIL, ec=INK2, lw=0.4, label="fails: label = binding constraint"),
         Line2D([], [], color=ECON, lw=2, label="B/C ≥ 1 at low cost (inside line)"),
         Line2D([], [], color=INK, lw=1, ls=(0, (4, 3)), label="10⁻³ of local ambient")]
    fig.legend(handles=h, loc="lower center", ncol=5, frameon=False, fontsize=8.5)
    fig.suptitle("Skyloop operating regime: pressure × Mach, 50-seat 17 t pod, best bore ≤ 8 m, "
                 "10 s moving block", x=0.01, ha="left", fontsize=11.5, color=INK, fontweight="bold")
    fig.tight_layout(rect=(0, 0.04, 1, 0.97))
    fig.savefig(os.path.join(HERE, "regime_maps.png"), dpi=130, facecolor=SURF)
    print("wrote regime_maps.png")


def fig_pods():
    d = pd.read_csv(os.path.join(HERE, "regime_pods.csv"))
    cases = [("elevated", "maglev", "brick-wall"), ("elevated", "maglev", "moving block 10 s"),
             ("buried", "maglev", "moving block 10 s"), ("elevated", "glide", "moving block 10 s")]
    fig, axs = plt.subplots(1, 4, figsize=(14.5, 4.6), sharey=True)
    fig.patch.set_facecolor(SURF)
    for ax, (sk, sup, rule) in zip(axs, cases):
        sub = d[(d.setting == sk) & (d.support == sup) & (d.rule == rule)].copy()
        sub["pod_t"] = sub["pod_t"].round(2)
        panel(ax, sub, "Mach", "pod_t", "ops_status", "ops_binding", False, econ_col="BC_low", label_marginal=False)
        style(ax, f"{SET_T[sk]} · {sup}\n{rule}")
        ax.set_xlabel("Mach", fontsize=9)
    axs[0].set_ylabel("pod mass, t  (seats ≈ (m − 2 t) / 0.3 t)", fontsize=9)
    fig.suptitle("Pod size × Mach (best pressure and bore per point). Orange: B/C ≥ 1 at low cost",
                 x=0.01, ha="left", fontsize=11, color=INK, fontweight="bold")
    fig.tight_layout(rect=(0, 0, 1, 0.93))
    fig.savefig(os.path.join(HERE, "regime_pods.png"), dpi=130, facecolor=SURF)
    print("wrote regime_pods.png")


def fig_glide():
    """Pressure needed to glide (C_L <= 1) vs pressure allowed by heating, per Mach."""
    fig, axs = plt.subplots(1, 2, figsize=(12, 4.4))
    fig.patch.set_facecolor(SURF)
    Ms = np.round(np.arange(1.2, 8.01, 0.2), 2)
    ps = np.geomspace(0.1, 5000, 160)
    for ax, sk in zip(axs, ["floating", "elevated"]):
        for m_pod, lab, col in [(3500, "docs dart 3.5 t, 5 seats", "#1baf7a"), (R.REF_MASS, "50 seats, 17 t", "#2a78d6")]:
            need, heat = [], []
            for M in Ms:
                pn = ph = np.nan
                for p in ps:
                    best = None
                    for b in R.BETAS:
                        if R.pod(m_pod)["d"] / math.sqrt(b) > R.D_MAX:
                            continue
                        r = R.evaluate(sk, "glide", float(M), float(p), m_pod, b, "moving block 10 s")
                        best = r if best is None or r["m_support"] > best["m_support"] else best
                    if np.isnan(pn) and best["m_support"] >= 0:
                        pn = p
                    if best["m_heating"] >= 0:
                        ph = p
                need.append(pn)
                heat.append(ph)
            ax.plot(Ms, need, color=col, lw=2, label=f"{lab}: lift needs ≥")
            ax.plot(Ms, heat, color=col, lw=2, ls=(0, (2, 2)), label=f"{lab}: heating allows ≤")
        pa = R.SETTINGS[sk]["p_amb"]
        ax.axhline(1e-3 * pa, color=INK, lw=1, ls=(0, (4, 3)))
        ax.text(1.3, 1e-3 * pa * 1.15, f"10⁻³ ambient = {1e-3 * pa:.1f} Pa", fontsize=7.5, color=INK)
        ax.set_yscale("log")
        ax.set_ylim(0.3, 5000)
        style(ax, f"{SET_T[sk]}: aerodynamic glide window (C_L ≤ 1, bore ≤ 8 m)")
        ax.set_xlabel("Mach", fontsize=9)
        ax.set_ylabel("tube pressure, Pa", fontsize=9)
        ax.grid(True, color=GRID, lw=0.6)
        ax.legend(fontsize=7.5, frameon=False, loc="lower left")
    fig.tight_layout()
    fig.savefig(os.path.join(HERE, "regime_glide.png"), dpi=130, facecolor=SURF)
    print("wrote regime_glide.png")


if __name__ == "__main__":
    fig_maps()
    fig_pods()
    fig_glide()
