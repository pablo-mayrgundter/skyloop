#!/usr/bin/env python3
"""Figures for highmach.md (reads this folder's CSVs only)."""
from __future__ import annotations

from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

HERE = Path(__file__).resolve().parent
FIG = HERE / "figures"
FIG.mkdir(exist_ok=True)
C1, C2, C3, C4, INK, INK2, GRID = "#2a78d6", "#eb6834", "#1baf7a", "#6250d6", "#0b0b0b", "#52514e", "#e4e3df"
GREY = "#a9a8a3"
plt.rcParams.update({"font.size": 9, "axes.edgecolor": INK2, "axes.labelcolor": INK, "xtick.color": INK2,
                     "ytick.color": INK2, "font.family": "DejaVu Sans"})


def tidy(a):
    a.grid(color=GRID, lw=0.6)
    a.set_axisbelow(True)
    a.spines[["top", "right"]].set_visible(False)


def fig_time():
    T = pd.read_csv(HERE / "door_to_door.csv")
    t = T[(T.pair == "NYC-LAX") & (T.overhead == "base (130 min)")].set_index("mode")
    order = ["subsonic", "Boomless M1.15", "shield M1.7", "shield M2.2", "shield M3.0", "shield M3.5", "shield M4.0",
             "shield M-inf, M4 thrust", "floor: no cruise, no ramps"]
    labels = ["Subsonic M0.85", "Boomless M1.15", "Shield M1.7", "Shield M2.2", "Shield M3.0", "Shield M3.5",
              "Shield M4.0", "Shield, infinite Mach\n(M4-concept thrust)", "Floor: no cruise,\nno ramps"]
    comp = [("Ground: access, processing, egress", ["access", "processing", "egress"], GREY),
            ("Taxi, climb, descent, approach, merge", ["taxi", "climb", "descent", "approach", "merge"], C4),
            ("Acceleration + deceleration ramps", ["accel", "decel"], C2),
            ("Cruise", ["cruise"], C1)]
    fig, (a1, a2) = plt.subplots(1, 2, figsize=(13.5, 4.8), dpi=150, gridspec_kw=dict(width_ratios=[1.35, 1]))
    y = np.arange(len(order))[::-1]
    left = np.zeros(len(order))
    for lab, cols, col in comp:
        v = np.array([t.loc[o, cols].sum() for o in order]) / 60
        a1.barh(y, v, left=left, height=0.62, color=col, edgecolor="white", linewidth=1.5, label=lab)
        left += v
    for yy, o in zip(y, order):
        a1.text(t.loc[o, "door_to_door_h"] + 0.06, yy, f"{t.loc[o, 'door_to_door_h']:.2f} h", va="center",
                fontsize=8, color=INK)
    a1.axvline(4.0, color=INK, lw=0.9, ls=(0, (4, 3)))
    a1.text(4.05, y[-1] - 0.45, "4 h day-return limit (base)", fontsize=7.5, color=INK2, va="bottom")
    a1.set_yticks(y)
    a1.set_yticklabels(labels, fontsize=8)
    a1.set_xlabel("Door-to-door hours, NYC-LAX (130 min ground overhead)")
    a1.set_xlim(0, 8.0)
    a1.legend(fontsize=7.5, frameon=False, loc="lower right")
    a1.set_title("(a) Where the time goes: cruise is already the minority at Mach 2", loc="left", fontsize=9.5)
    tidy(a1)
    S = pd.read_csv(HERE / "d2d_sensitivity.csv")
    S["Mx"] = S.M.replace(np.inf, 5.0)
    for oh, ls in (("base (130 min)", "-"), ("fast (100 min)", (0, (5, 2))), ("walk-on (70 min)", (0, (1.5, 1.5)))):
        for thr, col in (("Overture thrust (0.23 m/s2)", C2), ("concept thrust", C1), ("2x concept thrust", C3)):
            s = S[(S.overhead == oh) & (S.thrust == thr)].sort_values("Mx")
            a2.plot(s.Mx, s.d2d_h, ls=ls, color=col, lw=1.8, marker="o", ms=3,
                    label=f"{thr}" if oh == "base (130 min)" else None)
        s = S[(S.overhead == oh) & (S.thrust == "concept thrust")].sort_values("Mx")
        a2.text(5.08, s.d2d_h.iloc[-1], oh.split(" (")[0] + f"\n{oh.split('(')[1][:-1]}", fontsize=7, color=INK2, va="center")
    a2.axhline(4.0, color=INK, lw=0.9, ls=(0, (4, 3)))
    bl = pd.read_csv(HERE / "door_to_door.csv")
    blv = bl[(bl.pair == "NYC-LAX") & (bl["mode"] == "Boomless M1.15") & (bl.overhead == "base (130 min)")].door_to_door_h.iloc[0]
    a2.axhline(blv, color=GREY, lw=1.2)
    a2.text(1.72, blv + 0.05, f"Boomless, 130 min: {blv:.2f} h", fontsize=7.5, color=INK2)
    a2.set_xticks([1.7, 2.0, 2.5, 3.0, 3.5, 4.0, 5.0])
    a2.set_xticklabels(["1.7", "2", "2.5", "3", "3.5", "4", "inf"])
    a2.set_xlim(1.6, 5.6)
    a2.set_xlabel("Cruise Mach (shield on own great circle)")
    a2.set_ylabel("Door-to-door hours, NYC-LAX")
    a2.legend(fontsize=7.2, frameon=False, loc="upper right", title="Line colour = thrust; style = overhead",
              title_fontsize=7.2)
    a2.set_title("(b) Ground time matters more than Mach above ~2.5", loc="left", fontsize=9.5)
    tidy(a2)
    fig.tight_layout()
    fig.savefig(FIG / "highmach_time.png", facecolor="white")


def fig_shield():
    CR = pd.read_csv(HERE / "ceiling.csv")
    SM = pd.read_csv(HERE / "shield_vs_mach.csv")
    fig, (a1, a2) = plt.subplots(1, 2, figsize=(12.5, 4.4), dpi=150)
    for x, col, lab in ((0.0, C1, "pure H2"), (0.01, C3, "1 % air"), (0.02, C2, "2 % air"), (0.05, C4, "5 % air")):
        c = CR[(CR.H_km == 20) & (CR.air_frac == x)].sort_values("T_gas")
        a1.plot(c.T_gas, c.ceiling_mach, "-o", ms=3, color=col, lw=1.8)
        a1.text(c.T_gas.iloc[-1] + 3, c.ceiling_mach.iloc[-1], lab, fontsize=7.5, color=INK2, va="center")
    for M, lab in ((3.0, "concept B M3.0"), (3.5, "concept C M3.5"), (4.0, "concept D M4.0")):
        a1.axhline(M, color=GREY, lw=0.9, ls=(0, (4, 3)))
        a1.text(212, M + 0.03, lab, fontsize=7, color=INK2)
    a1.axvline(216.65, color=GREY, lw=0.8)
    a1.text(219, 2.82, "ambient 216.65 K", fontsize=7, color=INK2)
    a1.set_xlabel("H2 gas temperature, K (20 km; heat loss ~1.5 W/m2K per K of superheat)")
    a1.set_ylabel("Evanescent ceiling: Mach at which V = c_gas")
    a1.set_xlim(205, 385)
    a1.set_title("(a) The ceiling: purity costs more than warmth gains", loc="left", fontsize=9.5)
    tidy(a1)
    s = SM[(SM.H_mode == "preferred")]
    for gas_, col in (("ambient", C1), ("+40 K", C3), ("+80 K (solar noon)", C2)):
        for x, ls in ((0.0, "-"), (0.02, (0, (4, 2)))):
            c = s[(s.gas == gas_) & (s.air_frac == x)].sort_values("M")
            c = c[np.isfinite(c.d_min_m)]
            a2.plot(c.M, c.d_min_m, ls=ls, color=col, lw=1.8, label=f"{gas_}, {'pure' if x == 0 else '2 % air'}")
    a2.annotate("TBCC weight step\n(63 -> 86 t mid-cruise)", (3.08, 85), xytext=(3.2, 40), fontsize=7, color=INK2,
                arrowprops=dict(arrowstyle="-", color=INK2, lw=0.8))
    a2.axhline(54, color=GREY, lw=0.9)
    a2.text(1.72, 57, "repo baseline 54-55 m (Concorde 147 t, M2)", fontsize=7.2, color=INK2)
    a2.set_xlabel("Cruise Mach (concept weight; shield at preferred altitude, q = 25 kPa)")
    a2.set_ylabel("Min H2 thickness for 0.11 psf, m")
    a2.legend(fontsize=7, frameon=False, loc="upper left")
    a2.set_title("(b) Required layer: weight x Mach, then a cliff at the ceiling", loc="left", fontsize=9.5)
    tidy(a2)
    fig.tight_layout()
    fig.savefig(FIG / "highmach_shield.png", facecolor="white")


def fig_aircraft():
    A = pd.read_csv(HERE / "aircraft_param.csv")
    A = A[A.closes == True]  # noqa: E712
    fig, ax = plt.subplots(1, 4, figsize=(14, 3.5), dpi=150)
    for kT, col, lab in ((1.22, C2, "fixed-cycle thrust (k_T 1.22)"), (1.40, C1, "adaptive, +15 % (1.40)"),
                         (1.60, C3, "+30 % thrust (1.60)")):
        a = A[A.kT == kT].sort_values("M")
        ax[3].plot(a.M, a.a_acc, color=col, lw=1.8, label=lab)
    a = A[A.kT == 1.40].sort_values("M")
    for i, (col_, yl, ttl) in enumerate((("LD", "Cruise L/D", "(a) L/D"),
                                         ("t_per_seat", "MTOW per seat, t", "(b) Weight per seat"),
                                         ("fuel_vs_subsonic_seat", "Fuel per seat-km / subsonic", "(c) Fuel burn"))):
        ax[i].plot(a.M, a[col_], color=C1, lw=1.8)
        ax[i].set_ylabel(yl)
        ax[i].set_title(ttl, loc="left", fontsize=9.5)
    for M in (2.25, 3.05):
        for i in range(4):
            ax[i].axvline(M, color=GREY, lw=0.8, ls=(0, (3, 3)))
    ax[1].text(2.28, ax[1].get_ylim()[1] * 0.97, "titanium", fontsize=7, color=INK2, va="top")
    ax[1].text(3.08, ax[1].get_ylim()[1] * 0.97, "TBCC", fontsize=7, color=INK2, va="top")
    ax[3].set_ylabel("Mean ramp acceleration, m/s2")
    ax[3].set_title("(d) Ramp acceleration", loc="left", fontsize=9.5)
    ax[3].legend(fontsize=7, frameon=False)
    for a_ in ax:
        a_.set_xlabel("Cruise Mach")
        tidy(a_)
    fig.tight_layout()
    fig.savefig(FIG / "highmach_aircraft.png", facecolor="white")


def fig_bcr():
    R = pd.read_csv(HERE / "bcr_by_mach.csv")
    R = R[(R.network == "spine") & (R.demand == "base") & (R.capture == 1.0) & ~R.design.str.startswith("repo")]
    S = pd.read_csv(HERE / "summary_by_mach.csv")
    S = S[~S.design.str.startswith("repo")].sort_values("M")
    fig, (a1, a2) = plt.subplots(1, 2, figsize=(12.5, 4.4), dpi=150)
    for col_, lab, c in (("bl_base_100", "base cost", C4), ("bl_plausible_100", "plausible cost", C1),
                         ("bl_aggressive_100", "aggressive cost", C3), ("bl_integrated_100", "integrated (sibling)", C2)):
        a1.plot(S.M, S[col_], "-o", color=c, lw=1.8, ms=4)
        a1.text(4.05, S[col_].iloc[-1], lab, fontsize=7.5, color=INK2, va="center")
    a1.axhline(1.0, color=INK, lw=0.9)
    a1.text(1.72, 1.02, "break-even", fontsize=7.5, color=INK)
    a1.set_xlim(1.6, 4.6)
    a1.set_ylim(0, 1.15)
    a1.set_xlabel("Cruise Mach (spine, 100 % capture, base demand, two tether rows)")
    a1.set_ylabel("B/C vs Boomless (time + induced; fuel not charged)")
    a1.set_title("(a) Mach roughly doubles B/C; it stops rising after M2.5", loc="left", fontsize=9.5)
    tidy(a1)
    b = R[(R.vs == "bl") & (R.cost == "base")].sort_values("M")
    a2.plot(b.M, b.benefit_busd, "-o", color=C1, lw=1.8, ms=4)
    a2.plot(b.M, b.extra_fuel_busd, "-o", color=C2, lw=1.8, ms=4)
    a2.plot(b.M, b.extra_fuel_busd + b.carbon_busd, ls=(0, (4, 2)), color=C2, lw=1.6)
    a2.plot(b.M, b.cost_busd, "-", color=GREY, lw=1.6)
    for y, lab in ((b.benefit_busd.iloc[-1], "time + induced benefit"), (b.extra_fuel_busd.iloc[-1], "extra fuel vs Overture-class"),
                   (b.extra_fuel_busd.iloc[-1] + b.carbon_busd.iloc[-1], "extra fuel + CO2 ($190/t)"),
                   (b.cost_busd.iloc[-1], "shield cost, base")):
        a2.text(4.05, y, lab, fontsize=7.5, color=INK2, va="center")
    a2.set_xlim(1.6, 4.9)
    a2.set_xlabel("Cruise Mach (spine, 100 % capture vs Boomless)")
    a2.set_ylabel("$ B per year")
    a2.set_title("(b) Above ~M2.3 the extra fuel outweighs the time saved", loc="left", fontsize=9.5)
    tidy(a2)
    fig.tight_layout()
    fig.savefig(FIG / "highmach_bcr.png", facecolor="white")


if __name__ == "__main__":
    fig_time()
    fig_shield()
    fig_aircraft()
    if (HERE / "bcr_by_mach.csv").exists():
        fig_bcr()
    print("figures written")
