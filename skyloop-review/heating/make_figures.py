"""Figures for heating.md, read from the CSVs written by heating.py:
heating_maps.png (Mach x pressure feasibility, W/S x bore), heating_pod.png (scaling, temperatures, sink, radius),
heating_tube.png (gas dT per pass, max pods/h)."""
import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from matplotlib.colors import ListedColormap  # noqa: E402
from matplotlib.lines import Line2D  # noqa: E402
from matplotlib.patches import Patch  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.dont_write_bytecode = True
import heating as H  # noqa: E402

INK, INK2, GRID, SURF = "#0b0b0b", "#52514e", "#e4e3df", "#fcfcfb"
PASS, MARG, FAIL = "#2a78d6", "#a9c8ef", "#ecebe7"
S1, S2, S3, S4 = "#2a78d6", "#eb6834", "#1baf7a", "#eda100"
CMAP = ListedColormap([FAIL, MARG, PASS])
SHORT = {"lift": "lift", "span": "span > bore", "trip": "trip", "skin_T": "skin T", "le_T": "nose/LE T",
         "heat_sink": "heat sink", "gas_dT_pass": "gas ΔT\nper pass", "throughput": "gas heat\n→ pods/h",
         "film_T": "film T", "guideway_T": "guideway T"}


def style(ax, title):
    ax.set_facecolor(SURF)
    ax.set_title(title, loc="left", fontsize=9.5, color=INK, fontweight="bold")
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    for s in ("left", "bottom"):
        ax.spines[s].set_color(INK2)
    ax.tick_params(colors=INK2, labelsize=8)
    ax.xaxis.label.set_color(INK2)
    ax.yaxis.label.set_color(INK2)
    ax.grid(color=GRID, lw=0.6)
    ax.set_axisbelow(True)


def edges(v, log=False):
    v = np.asarray(v, float)
    if log:
        lv = np.log10(v)
        e = np.concatenate([[lv[0] - (lv[1] - lv[0]) / 2], (lv[1:] + lv[:-1]) / 2, [lv[-1] + (lv[-1] - lv[-2]) / 2]])
        return 10 ** e
    return np.concatenate([[v[0] - (v[1] - v[0]) / 2], (v[1:] + v[:-1]) / 2, [v[-1] + (v[-1] - v[-2]) / 2]])


def panel(ax, sub, scol="status_ns", bcol="binding_ns"):
    xs = np.sort(sub.Mach.unique())
    ys = np.sort(sub.p_Pa.unique())
    piv_s = sub.pivot_table(index="p_Pa", columns="Mach", values=scol, aggfunc="first").loc[ys, xs]
    piv_b = sub.pivot_table(index="p_Pa", columns="Mach", values=bcol, aggfunc="first").loc[ys, xs]
    Z = piv_s.replace({"fail": 0, "marginal": 1, "pass": 2}).to_numpy(float)
    B = piv_b.to_numpy(object)
    ax.pcolormesh(edges(xs), edges(ys, True), Z, cmap=CMAP, vmin=-0.5, vmax=2.5, shading="flat")
    codes = {k: n for n, k in enumerate(sorted(set(B.ravel())))}
    C = np.vectorize(lambda b: codes[b])(B).astype(float)
    C[Z == 2] = np.nan
    if len(codes) > 1:
        ax.contour(xs, ys, np.nan_to_num(C, nan=-1), levels=np.arange(-0.5, len(codes)), colors=INK2,
                   linewidths=0.5, alpha=0.6)
    for b, n in codes.items():
        mask = C == n
        if mask.sum() < 0.03 * mask.size:
            continue
        jj, ii = np.nonzero(mask)
        k = np.argmin((jj - np.median(jj)) ** 2 + (ii - np.median(ii)) ** 2)
        ax.text(xs[ii[k]], ys[jj[k]], SHORT.get(b, b), fontsize=7, color=INK, ha="center", va="center",
                bbox=dict(boxstyle="round,pad=0.15", fc=SURF, ec="none", alpha=0.85))
    ax.set_yscale("log")


def fig_maps():
    d = pd.read_csv(os.path.join(HERE, "heating_grid.csv"))
    bl = pd.read_csv(os.path.join(HERE, "heating_bestld.csv"))
    wss, ds = [50.0, 150.0, 500.0], H.D_LIST
    fig, axs = plt.subplots(3, 3, figsize=(13.5, 11.5), sharex=True, sharey=True)
    fig.patch.set_facecolor(SURF)
    for i, ws in enumerate(wss):
        for j, D in enumerate(ds):
            ax = axs[i, j]
            sub = d[(d.WS_Pa == ws) & (d.D_bore == D)]
            panel(ax, sub)
            b = bl[bl.WS_Pa == ws]
            ax.plot(b.Mach, b.p_bestLD_Pa, color=INK, lw=1.4, ls=(0, (5, 3)))
            ax.axhline(300, color=INK2, lw=0.8, ls=":")
            ax.axhline(5.5, color=S2, lw=1.2)
            sp = sub.span_m.dropna().iloc[0]
            fit = "fits" if sp <= 0.5 * D else ("tight" if sp <= 0.8 * D else "does NOT fit")
            style(ax, f"W/S {ws:g} Pa · bore {D:g} m · span {sp:.0f} m {fit}")
            if sp > 0.5 * D:
                ax.add_patch(plt.Rectangle((0, 0), 1, 1, transform=ax.transAxes, fill=False, hatch="////",
                                           ec=INK2 if sp > 0.8 * D else GRID, lw=0, alpha=0.5))
            if j == 0:
                ax.set_ylabel("tube pressure, Pa", fontsize=9)
            if i == 2:
                ax.set_xlabel("Mach (V = M × 295 m/s)", fontsize=9)
    h = [Patch(fc=PASS, label="all checks pass"), Patch(fc=MARG, label="marginal (worst check within its tolerance)"),
         Patch(fc=FAIL, ec=INK2, lw=0.4, label="fails: label = binding constraint"),
         Line2D([], [], color=INK, lw=1.4, ls=(0, (5, 3)), label="best-L/D pressure"),
         Line2D([], [], color=S2, lw=1.2, label="5.5 Pa (10⁻³ of ambient)"),
         Line2D([], [], color=INK2, lw=0.8, ls=":", label="300 Pa (study range top)"),
         Patch(fc="none", ec=INK2, hatch="////", lw=0, label="hatched: wing span > 0.5 × bore (dark: > 0.8 ×)")]
    fig.legend(handles=h, loc="lower center", ncol=3, frameon=False, fontsize=8.5)
    fig.suptitle("Heating-feasible region, 3.5 t winged pod in a wide tube: Mach × pressure by wing loading and bore\n"
                 "(span-vs-bore check shown by hatching, not colour; 120 pods/h, 4,000 km, Inconel 1,250 K skin, C/SiC 1,700 K nose R 0.1 m / LE R 3 cm)",
                 x=0.01, ha="left", fontsize=11, color=INK, fontweight="bold")
    fig.tight_layout(rect=(0, 0.05, 1, 0.95))
    fig.savefig(os.path.join(HERE, "heating_maps.png"), dpi=120, facecolor=SURF)
    plt.close(fig)
    print("wrote heating_maps.png")


def fig_pod():
    sc = pd.read_csv(os.path.join(HERE, "heating_scaling.csv"))
    sk = pd.read_csv(os.path.join(HERE, "heating_sink.csv"))
    rd = pd.read_csv(os.path.join(HERE, "heating_radius.csv"))
    fig, axs = plt.subplots(2, 2, figsize=(12.5, 9.0))
    fig.patch.set_facecolor(SURF)
    # (a) trip heat
    ax = axs[0, 0]
    cases = [("q 300 Pa (C_L~0.5)", S2, "q 300 Pa"), ("q 1 kPa (C_L~0.15)", S3, "q 1 kPa"),
             ("q 3 kPa (C_L~0.05)", S4, "q 3 kPa"), ("best L/D", S1, "best L/D")]
    for c, col, lab in cases:
        s = sc[sc.case == c]
        ax.plot(s.Mach, s.Q_in_trip_GJ, color=col, lw=2)
        n = s.exp_Q_in_trip_GJ.iloc[0]
        dy = {"q 3 kPa": 1.12, "best L/D": 0.9}.get(lab, 1.0)
        ax.text(s.Mach.iloc[-1] + 0.08, s.Q_in_trip_GJ.iloc[-1] * dy, f"{lab} (∝V^{n:.1f})", fontsize=7.5, color=INK, va="center")
    s = sc[sc.case == "best L/D"]
    ax.plot(s.Mach, s.drag_energy_trip_GJ, color=INK2, lw=1.2, ls="--")
    ax.text(3.05, s.drag_energy_trip_GJ.iloc[0] * 1.08, "drag energy, best L/D", fontsize=7.5, color=INK2)
    ax.set_yscale("log")
    ax.set_xlim(3, 9.6)
    ax.set_xlabel("Mach")
    ax.set_ylabel("heat into pod per 4,000 km trip, GJ")
    style(ax, "(a) Heat per trip at fixed q: grows slowly, ~5–15 % of drag energy")
    # (b) temperatures on best L/D line
    ax = axs[0, 1]
    s = sc[sc.case == "best L/D"]
    for k, col, lab in (("T_nose", S2, "nose, R 0.1 m"), ("T_le", S4, "wing LE, R 3 cm"),
                        ("T_wing_lower", S1, "wing lower surface"), ("T_wing_upper", S3, "wing upper surface")):
        ax.plot(s.Mach, s[k], color=col, lw=2)
        ax.text(8.05, s[k].iloc[-1], lab, fontsize=7.5, color=INK, va="center")
    for T, lab in ((1700, "C/SiC LE 1,700 K"), (1250, "Inconel TPS 1,250 K"), (700, "Ti 700 K"), (420, "Al / CFRP 420 K")):
        ax.axhline(T, color=INK2, lw=0.8, ls=":")
        ax.text(3.02, T + 15, lab, fontsize=7, color=INK2)
    ax.set_xlim(3, 9.4)
    ax.set_xlabel("Mach")
    ax.set_ylabel("radiative-equilibrium temperature, K")
    style(ax, "(b) Temperatures on the best-L/D line, W/S 150 Pa")
    # (c) sink fraction
    ax = axs[1, 0]
    sub = sk[(sk.L_trip_km == 4000) & (sk.p_Pa == 30)]
    for m, col in zip(sorted(sub.materials.unique(), key=lambda x: -sub[sub.materials == x].f_sink_water.max()),
                      (S2, S4, S1, S3)):
        t = sub[sub.materials == m]
        ax.plot(t.Mach, 100 * t.f_sink_water.clip(lower=1e-4), color=col, lw=2)
        if m.startswith("Al"):
            ax.text(8.05, 100 * t.f_sink_water.iloc[-1], "Al / CFRP skin 420 K", fontsize=7.5, color=INK, va="center")
    ax.text(8.05, 0.1, "Ti 700 K, Inconel 1,250 K,\nC/SiC 1,600 K (identical:\ncabin load only)", fontsize=7.5,
            color=INK, va="center")
    ax.axhline(5, color=INK2, lw=0.8, ls=":")
    ax.text(3.02, 5.5, "5 % pass", fontsize=7, color=INK2)
    ax.axhline(15, color=INK2, lw=0.8, ls=":")
    ax.text(3.02, 16.5, "15 % marginal", fontsize=7, color=INK2)
    ax.set_yscale("log")
    ax.set_xlim(3, 9.6)
    ax.set_xlabel("Mach")
    ax.set_ylabel("vented-water heat sink, % of pod mass")
    style(ax, "(c) Heat sink for NYC–LAX, 30 Pa, W/S 150 Pa")
    # (d) radius
    ax = axs[1, 1]
    for M, col in ((4.0, S3), (6.0, S1), (8.0, S2)):
        t = rd[(rd.Mach == M) & (rd.q_Pa == 1000)]
        ax.plot(t.R_m * 100, t.T_nose_eq, color=col, lw=2)
        ax.plot(t.R_m * 100, t.T_le_eq, color=col, lw=1.2, ls="--")
        ax.text(0.5, t.T_nose_eq.iloc[0] + 20, f"M{M:g}", fontsize=8, color=INK)
    ax.axhline(1700, color=INK2, lw=0.8, ls=":")
    ax.text(30, 1715, "C/SiC 1,700 K", fontsize=7, color=INK2)
    ax.set_xscale("log")
    ax.set_xlabel("nose / leading-edge radius, cm")
    ax.set_ylabel("radiative-equilibrium temperature, K")
    ax.legend(handles=[Line2D([], [], color=INK2, lw=2, label="sphere nose"),
                       Line2D([], [], color=INK2, lw=1.2, ls="--", label="LE swept 63°")], fontsize=7.5, frameon=False)
    style(ax, "(d) Nose and LE vs radius, q = 1 kPa")
    fig.tight_layout()
    fig.savefig(os.path.join(HERE, "heating_pod.png"), dpi=120, facecolor=SURF)
    plt.close(fig)
    print("wrote heating_pod.png")


def fig_tube():
    t = pd.read_csv(os.path.join(HERE, "heating_tube.csv"))
    fig, axs = plt.subplots(1, 2, figsize=(12.5, 4.6))
    fig.patch.set_facecolor(SURF)
    sub = t[(t.WS_Pa == 150) & (t.Mach == 6)]
    for D, col in ((20, S2), (40, S4), (80, S1)):
        s = sub[sub.D_bore == D]
        axs[0].plot(s.p_Pa, s.dT_pass, color=col, lw=2, marker="o", ms=4)
        axs[0].text(s.p_Pa.iloc[-1] * 1.1, s.dT_pass.iloc[-1], f"{D} m", fontsize=8, color=INK, va="center")
        axs[1].plot(s.p_Pa, s.N_max_per_h, color=col, lw=2, marker="o", ms=4)
        axs[1].text(s.p_Pa.iloc[-1] * 1.1, s.N_max_per_h.iloc[-1], f"{D} m", fontsize=8, color=INK, va="center")
    for ax in axs:
        ax.set_xscale("log")
        ax.set_yscale("log")
        ax.set_xlabel("tube pressure, Pa")
    axs[0].axhline(22, color=INK2, lw=0.8, ls=":")
    axs[0].text(1.1, 25, "22 K pass (Δρ/ρ ≈ 10 %)", fontsize=7, color=INK2)
    axs[0].axhline(87, color=INK2, lw=0.8, ls=":")
    axs[0].text(1.1, 95, "87 K marginal", fontsize=7, color=INK2)
    axs[0].set_ylabel("bore-mean gas ΔT per pass, K")
    style(axs[0], "(a) Gas heating per pod pass, M6, W/S 150 Pa")
    axs[1].axhline(120, color=INK2, lw=0.8, ls=":")
    axs[1].text(1.1, 130, "120 pods/h assumed need", fontsize=7, color=INK2)
    axs[1].set_ylabel("max pods/h (gas +50 K mean)")
    style(axs[1], "(b) Gas-to-wall heat rejection caps throughput")
    for ax in axs:
        ax.set_xlim(0.9, 2000)
    fig.tight_layout()
    fig.savefig(os.path.join(HERE, "heating_tube.png"), dpi=120, facecolor=SURF)
    plt.close(fig)
    print("wrote heating_tube.png")


if __name__ == "__main__":
    fig_maps()
    fig_pod()
    fig_tube()
