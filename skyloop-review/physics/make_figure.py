"""Figure for physics.md: reads this folder's CSVs, writes physics.png."""
import csv
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
COL = {2.0: "#2a78d6", 4.0: "#eb6834", 6.0: "#1baf7a"}
INK, INK2, GRID, SURF = "#0b0b0b", "#52514e", "#e4e3df", "#fcfcfb"


def rd(name):
    return list(csv.DictReader(open(os.path.join(HERE, name))))


def style(ax, title):
    ax.set_facecolor(SURF)
    ax.set_title(title, loc="left", fontsize=10.5, color=INK, fontweight="bold")
    ax.grid(True, which="major", color=GRID, lw=0.6)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    for s in ("left", "bottom"):
        ax.spines[s].set_color(INK2)
    ax.tick_params(colors=INK2, labelsize=8.5)
    ax.xaxis.label.set_color(INK2)
    ax.yaxis.label.set_color(INK2)


fig, axs = plt.subplots(2, 2, figsize=(12, 8.6), facecolor=SURF)

# (a) per-metre budget, 20 km, D 3 m
mb = [r for r in rd("tube_mass_budget.csv") if r["alt_km"] == "20" and r["D_m"] == "3" and r["case"] == "central"][0]
sens = {r["case"]: r for r in rd("envelope_sensitivity.csv")}
base = sens["baseline: C_D 0.6, V/H 1.5, film 0.10"]
best = sens["C_D 0.3 and V/H 0.75"]
items = [
    ("tube's own vacuum lift", float(mb["tube_own_vacuum_lift_kg_m"]), "#9085e9"),
    ("tube wall (buckling, 5.5 kPa)", float(mb["wall_kg_m"]), "#2a78d6"),
    ("guideway / levitation rail", float(mb["guideway_kg_m"]), "#2a78d6"),
    ("LSM stator, route average (M4, 1,000 km)", float(mb["stator_avg_kg_m"]), "#2a78d6"),
    ("tube + track total", float(mb["payload_kg_m"]), "#184f95"),
    ("H2 envelope film", float(base["film_kg_m"]), "#eb6834"),
    ("tether hold-down", float(base["holddown_kg_m"]), "#eb6834"),
    ("H2 lift needed, best case", float(best["lift_kg_m"]), "#1baf7a"),
    ("H2 lift needed, baseline", float(base["lift_kg_m"]), "#1baf7a"),
]
ax = axs[0, 0]
y = list(range(len(items)))[::-1]
for yi, (lab, v, c) in zip(y, items):
    ax.barh(yi, v, color=c, height=0.62, left=0.1)
    ax.text(v * 1.12, yi, f"{v:,.1f}" if v < 10 else f"{v:,.0f}", va="center", fontsize=8.5, color=INK)
ax.set_yticks(y)
ax.set_yticklabels([i[0] for i in items], fontsize=8.5)
ax.set_xscale("log")
ax.set_xlim(0.2, 8000)
ax.set_xlabel("kg per metre of line (log)")
style(ax, "a. Mass and lift per metre, D 3 m tube at 20 km")

# (b) aero energy vs pressure
ax = axs[0, 1]
dr = rd("drag_vs_pressure.csv")
for M in (2.0, 4.0, 6.0):
    pts = [(float(r["p_Pa"]), float(r["aero_Wh_per_paxkm"])) for r in dr
           if float(r["Mach"]) == M and r["D_tube_m"] == "3"]
    ax.plot(*zip(*pts), color=COL[M], lw=2, marker="o", ms=4)
    ax.text(pts[-1][0] * 1.15, pts[-1][1], f"M{M:g}", color=INK, fontsize=9, va="center")
ax.axhspan(0.5, 2.0, color="#e4e3df", zorder=0)
ax.text(125, 1.0, "docs: 0.5-2 Wh total", fontsize=8.5, color=INK2, va="center")
ax.axvline(101, color=INK2, lw=1, ls="--")
ax.text(110, 0.012, "docs' '<10^-3 atm'", fontsize=8.5, color=INK2, rotation=90, va="bottom")
ax.set_xscale("log")
ax.set_yscale("log")
ax.set_xlim(0.08, 3000)
ax.set_xlabel("tube pressure, Pa")
ax.set_ylabel("aero drag energy, Wh per passenger-km")
style(ax, "b. Drag energy vs tube pressure (D 3 m, unchoked)")

# (c) alignment
ax = axs[1, 0]
al = rd("alignment_allowable.csv")
for M in (2.0, 4.0, 6.0):
    pts = [(float(r["wavelength_m"]) / 1e3, float(r["A_allow_0p05g_m"])) for r in al if float(r["Mach"]) == M]
    ax.plot(*zip(*pts), color=COL[M], lw=2)
    ax.text(pts[-1][0] * 1.08, pts[-1][1], f"M{M:g}", color=INK, fontsize=9, va="center")
mo = rd("floating_tube_motion.csv")
lo = min(float(r["tether_top_motion_m"]) for r in mo if r["u_gust_m_s"] not in ("", "nan"))
hi = max(float(r["total_m"]) for r in mo if r["u_gust_m_s"] not in ("", "nan"))
ax.fill_between([3, 12], lo, hi, color="#eda100", alpha=0.35, lw=0)
ax.text(3.2, hi * 1.25, f"floating tube under 1-5 m/s gusts: {lo:.0f}-{hi:.0f} m", fontsize=8.5, color=INK)
ax.axhline(0.93, color=INK2, lw=1, ls="--")
ax.text(0.06, 1.1, "pod-to-wall clearance 0.93 m", fontsize=8.5, color=INK2)
ax.axhline(0.05, color=INK2, lw=1, ls=":")
ax.text(0.06, 0.058, "maglev suspension stroke ~5 cm", fontsize=8.5, color=INK2)
ax.set_xscale("log")
ax.set_yscale("log")
ax.set_xlim(0.05, 80)
ax.set_ylim(1e-5, 2000)
ax.set_xlabel("wavelength of tube deviation, km")
ax.set_ylabel("allowed amplitude at 0.05 g ride limit, m")
style(ax, "c. Alignment: allowed vs expected tube deviation")

# (d) accel length
ax = axs[1, 1]
pr = rd("propulsion.csv")
for M in (2.0, 4.0, 6.0):
    pts = [(float(r["a_g"]), float(r["L_acc_km"])) for r in pr if float(r["Mach"]) == M]
    ax.plot(*zip(*pts), color=COL[M], lw=2, marker="o", ms=4)
    ax.text(pts[0][0] * 0.93, pts[0][1], f"M{M:g}", color=INK, fontsize=9, va="center", ha="right")
ax.axvspan(0.2, 0.5, color="#e4e3df", zorder=0)
ax.text(0.21, 12, "passenger range 0.2-0.5 g", fontsize=8.5, color=INK2)
ax.set_xscale("log")
ax.set_yscale("log")
ax.set_xlim(0.06, 1.2)
ax.set_xlabel("sustained acceleration, g")
ax.set_ylabel("length to reach speed, km (same again to stop)")
style(ax, "d. Launch length vs acceleration")

fig.tight_layout()
out = os.path.join(HERE, "physics.png")
fig.savefig(out, dpi=130, facecolor=SURF)
print("wrote", out)
