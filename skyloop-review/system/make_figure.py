#!/usr/bin/env python3
"""figures/system.png: (a) NYC-LAX door to door, (b) capacity vs Mach, (c) $/km vs comparators, (d) net B/C by opening."""
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

HERE = Path(__file__).resolve().parent
S = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4", "#008300", "#6250d6", "#e34948"]
INK, INK2, GRID, SURF = "#0b0b0b", "#52514e", "#e4e3df", "#fcfcfb"
plt.rcParams.update({"font.size": 9, "axes.edgecolor": INK2, "axes.labelcolor": INK, "xtick.color": INK2,
                     "ytick.color": INK2, "axes.titlesize": 10, "axes.titleweight": "bold", "figure.facecolor": SURF,
                     "axes.facecolor": SURF, "axes.grid": True, "grid.color": GRID, "grid.linewidth": 0.6,
                     "axes.spines.top": False, "axes.spines.right": False})

fig, ax = plt.subplots(2, 2, figsize=(13, 10.5))

# (a) door to door
T = pd.read_csv(HERE / "d2d.csv")
T = T[T.pair == "NYC-LAX"]
pick = [("subsonic", "Subsonic M0.85"), ("Boomless M1.15", "Boomless M1.15"), ("BoomShield M2.2", "BoomShield M2.2"),
        ("Skyloop M2 (ground portal, 100 min)", "Skyloop M2, ground portal"),
        ("Skyloop M4 (ground portal, 100 min)", "Skyloop M4, ground portal"),
        ("Skyloop M6 (ground portal, 100 min)", "Skyloop M6, ground portal"),
        ("Skyloop M4 (ground portal, walk-on 70 min)", "Skyloop M4, walk-on 70 min"),
        ("Skyloop M4 (20 km station by lift (20.5 m/s))", "Skyloop M4, 20 km station by lift"),
        ("Skyloop M4 (20 km station by aircraft docking)", "Skyloop M4, aircraft docking")]
a = ax[0, 0]
parts = [("ground", "ground: access, processing, egress"), ("taxi_climb_desc", "taxi, climb, descent, approach"),
         ("ramps", "accel/brake ramps"), ("cruise", "cruise"), ("vertical", "vertical access to 20 km")]
y = np.arange(len(pick))[::-1]
for k, (m, lab) in enumerate(pick):
    r = T[T["mode"] == m].iloc[0]
    left = 0
    for j, (c, _) in enumerate(parts):
        w = r[c] / 60
        if w > 0:
            a.barh(y[k], w, left=left, color=S[j], height=0.62, edgecolor=SURF, linewidth=1.5)
        left += w
    a.text(left + 0.08, y[k], f"{r.d2d_h:.1f} h", va="center", color=INK, fontsize=8.5)
a.set_yticks(y, [p[1] for p in pick])
a.axvline(4.0, color=INK2, ls="--", lw=1)
a.text(4.03, y[0] + 0.45, "4 h: day-return limit", color=INK2, fontsize=8)
a.set_xlabel("door to door, hours (NYC-LAX)")
a.set_title("(a) Ground time, not Mach, sets the trip")
a.legend(handles=[plt.Rectangle((0, 0), 1, 1, color=S[j]) for j in range(len(parts))], labels=[p[1] for p in parts],
         fontsize=7.5, loc="upper center", bbox_to_anchor=(0.45, -0.12), ncol=3, frameon=False)
a.grid(axis="y", visible=False)

# (b) capacity
C = pd.read_csv(HERE / "capacity.csv")
b = ax[0, 1]
series = [("brick-wall", 50, "brick-wall, 50-seat pods"), ("brick-wall", 150, "brick-wall, 150-seat pods"),
          ("brick-wall", 4, "brick-wall, 4-seat pods (docs)"), ("moving block 10 s", 50, "moving block 10 s, 50 seats")]
for j, (rule, seats, lab) in enumerate(series):
    s = C[(C.rule == rule) & (C.seats == seats)].sort_values("M")
    b.plot(s.M, s.pax_yr_M, color=S[j], lw=2, marker="o", ms=5, label=lab)
b.axhline(81.4, color=INK2, ls="--", lw=1)
b.text(2.05, 95, "spine demand, 81 M pax/yr", color=INK2, fontsize=8)
b.set_yscale("log")
b.set_xlabel("Mach (V = M x 295 m/s)")
b.set_ylabel("pax per year per tube pair (both ways), M")
b.set_title("(b) Safe headway caps a tube far below demand")
b.legend(fontsize=7.5, frameon=False, loc="lower left")
b.set_xticks([2, 3, 4, 6])

# (c) cost per km
K = pd.read_csv(HERE / "comparators.csv")
c = ax[1, 0]
K = K.iloc[::-1].reset_index(drop=True)
for i, r in K.iterrows():
    col = S[0] if "this study" in r.system else (S[3] if "Skyloop doc" in r.system else S[2] if "BoomShield" in r.system
                                                  else "#9a9993")
    c.plot([r.usd_m_km_low, r.usd_m_km_high], [i, i], color=col, lw=5, solid_capstyle="round")
    c.plot(r.usd_m_km_mid, i, "o", color=SURF, ms=4)
c.set_yticks(range(len(K)), K.system, fontsize=7.5)
c.set_xscale("log")
B = pd.read_csv(HERE / "breakeven.csv")
be = B[(B.M == 4) & (B.opening == 2050) & (B.capture == 1.0) & (B.kind == "net") & (B.vs == "sub")
       & B.network.isin(["spine+LAX all-stop", "NYC-LAX direct", "NYC-MIA direct", "spine express"])]
lo, hi = be.breakeven_avg_tube_usd_m_km.min(), be.breakeven_avg_tube_usd_m_km.max()
c.axvspan(lo, hi, color=S[5], alpha=0.12, lw=0)
c.text(np.sqrt(lo * hi), len(K) - 0.3, f"Skyloop break-even\n(net, 100 % capture, 2050):\n${lo:.0f}-{hi:.0f} M/km",
       ha="center", va="top", fontsize=7.5, color=INK)
c.set_xlabel("$M per km (log), low-high")
c.set_title("(c) Cost per km against comparators")
c.grid(axis="y", visible=False)

# (d) net B/C by opening year
L = pd.read_csv(HERE / "lifecycle.csv")
d = ax[1, 1]
lines = [("spine+LAX all-stop", "base", "sub", "spine+LAX, base cost, vs subsonic"),
         ("spine+LAX all-stop", "low", "sub", "spine+LAX, low cost, vs subsonic"),
         ("NYC-LAX direct", "low", "bl", "NYC-LAX, low cost, vs Boomless"),
         ("NYC-LAX direct", "low", "m22", "NYC-LAX, low cost, vs BoomShield M2.2"),
         ("spine+LAX all-stop", "ground vactrain base", "sub", "same network as a ground vactrain ($62 M/km)")]
for j, (n, cost, vs, lab) in enumerate(lines):
    s = L[(L.network == n) & (L.cost == cost) & (L.vs == vs) & (L.M == 4) & (L.capture == 1.0)].sort_values("opening")
    d.plot(s.opening, s.bc_net, color=S[j], lw=2, marker="o", ms=5, label=lab)
    s3 = L[(L.network == n) & (L.cost == cost) & (L.vs == vs) & (L.M == 4) & (L.capture == 0.3)].sort_values("opening")
    d.plot(s3.opening, s3.bc_net, color=S[j], lw=1.2, ls=":", marker=None)
d.axhline(1, color=INK2, lw=1)
d.set_xlabel("opening year")
d.set_ylabel("life-cycle B/C, net of fuel, CO2 and electricity")
d.set_title("(d) Skyloop M4: pays only at low cost and full capture")
d.legend(fontsize=7.5, frameon=False, loc="upper center", bbox_to_anchor=(0.5, -0.12), ncol=2,
         title="solid 100 % capture, dotted 30 %", title_fontsize=7.5)
d.set_xticks([2045, 2050, 2055, 2060])
fig.tight_layout()
(HERE / "figures").mkdir(exist_ok=True)
fig.savefig(HERE / "figures" / "system.png", dpi=140)
print("ok")
