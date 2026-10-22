"""Figure for lid-fatigue.md: load per pass, lid response, allowable traffic by detail, trunk B/C vs cap."""
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
BLUE, ORANGE, AQUA, YELLOW = "#2a78d6", "#eb6834", "#1baf7a", "#eda100"
INK, INK2, GRID, BG = "#0b0b0b", "#52514e", "#e4e3df", "#fcfcfb"
plt.rcParams.update({"font.size": 8.5, "axes.edgecolor": INK2, "axes.labelcolor": INK, "xtick.color": INK2,
                     "ytick.color": INK2, "axes.titlesize": 9.5, "axes.titleweight": "bold"})


def clean(ax):
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    ax.grid(color=GRID, lw=0.6)
    ax.set_axisbelow(True)


wf = pd.read_csv(os.path.join(HERE, "lid_waveforms.csv"))
w = wf[wf.case.str.startswith("Concorde")].copy()
for k in ("w_top", "w_bot", "e", "x_gas"):          # reference displacements to 2 s before the pass
    w[k] = w[k] - w[k].iloc[0]
t0 = w.t_s[w.p_inc.abs() > 1].min()
t = (w.t_s - t0) * 1e3
fig = plt.figure(figsize=(11, 7.6), facecolor=BG)
gs = fig.add_gridspec(2, 2, hspace=0.42, wspace=0.42)

ax = fig.add_subplot(gs[0, 0], facecolor=BG)
ax.plot(t, w.p_inc, color=INK2, lw=1.4, label="incident (free field)")
ax.plot(t, w.p_top, color=BLUE, lw=2, label="on the lid (reflection-doubled)")
ax.plot(t, w.p_gas, color=ORANGE, lw=2, label="gas just under the lid")
ax.set_xlim(-50, 400)
ax.set_xlabel("time from bow shock, ms")
ax.set_ylabel("overpressure, Pa")
ax.set_title("(a) Load per pass: Concorde M2, h = 250 m, H₂ 55 m", loc="left")
ax.legend(frameon=False, fontsize=7.5, loc="upper right")
clean(ax)

ax = fig.add_subplot(gs[0, 1], facecolor=BG)
ax.plot(t, w.w_top * 100, color=BLUE, lw=2, label="lid heave")
ax.plot(t, w.e * 100, color=AQUA, lw=2, label="layer thickness change (tie stroke)")
ax.plot(t, w.x_gas * 100, color=ORANGE, lw=2, label="sideways gas excursion under lid (partition top)")
ax.set_xlim(-300, 1500)
ax.set_xlabel("time from bow shock, ms")
ax.set_ylabel("displacement, cm")
ax.set_title("(b) Lid response: heave ±10 cm; lid strain ~2×10⁻⁵", loc="left")
ax.set_ylim(-15, 100)
ax.legend(frameon=False, fontsize=7.5, loc="upper right")
clean(ax)

# (c) allowable passes/day by detail, Concorde M2 55 m heavy lid, 30 yr, SF 10, PE/PET laminate lid
a = pd.read_csv(os.path.join(HERE, "allowable_traffic.csv"))
a = a[(a.life == "30 yr") & a.load_case.str.startswith("Concorde M2, H2 55 m, lid 4.45")]
dets = [("Joint flap, light (0.5 kg/m2, 1 m)", "PE/PET laminate", "Light joint flap"),
        ("Slack partition hem/crease (50 um)", "LLDPE balloon film", "Slack film partition crease"),
        ("Tie fitting, 50 mm grommet", "PE/PET laminate", "Tie at 50 mm grommet"),
        ("Tie fitting, 50 mm grommet", "Coated fabric (PU/PVF on polyester)", "Tie grommet, coated fabric"),
        ("Heat-sealed lap seam", "PE/PET laminate", "Heat-sealed seam"),
        ("Tie fitting, 0.5 m spread pad", "PE/PET laminate", "Tie on 0.5 m spread pad"),
        ("Joint flap, mass-matched", "PE/PET laminate", "Mass-matched joint flap"),
        ("Rolling-lobe partition / rim (R 0.5 m)", "LLDPE balloon film", "Rolling-lobe partition/rim"),
        ("Lid field (sheet)", "PE/PET laminate", "Lid field (4.8 mm sheet)")]
ax = fig.add_subplot(gs[1, 0], facecolor=BG)
for i, (d, m, lab) in enumerate(dets):
    x = a[(a.detail == d) & (a.material == m)].set_index("case")
    lo, ba, hi = (max(x.passes_per_day_boom_only[c], 0.03) for c in ("low", "base", "high"))
    ax.plot([lo, hi], [i, i], color=BLUE, lw=5, solid_capstyle="round", alpha=0.35)
    ax.plot([ba], [i], "o", color=BLUE, ms=6, mec=BG, mew=1.5)
ax.set_yticks(range(len(dets)))
ax.set_yticklabels([d[2] for d in dets], fontsize=7.5)
ax.set_xscale("log")
ax.set_xlim(0.03, 2e6)
for xv, lab, ha in ((100, "old cap\n100/day ", "right"), (4373, " busiest trunk km,\n full capture", "left")):
    ax.axvline(xv, color=INK2, lw=1, ls="--")
    ax.text(xv, 6.5, lab, fontsize=7, color=INK2, va="center", ha=ha,
            bbox=dict(fc=BG, ec="none", pad=1))
ax.set_xlabel("allowable passes/day (boom only), 30 yr, life factor 10; ≥10⁶ shown at 10⁶")
ax.set_title("(c) Fatigue cap by detail: low–high bar, base dot", loc="left")
clean(ax)

# (d) trunk B/C vs cap, M2.0 spine
b = pd.read_csv(os.path.join(HERE, "trunk_bcr_revised.csv"))
b = b[b.mach == 2.0]
order = ["old 100/day", "revised low 1,000/day", "revised base 10,000/day", "revised high (no fatigue cap)"]
labs = ["old\n100/day", "low\n1,000/day", "base\n10,000/day", "high\nno cap"]
ax = fig.add_subplot(gs[1, 1], facecolor=BG)
xx = np.arange(len(order))
for c, col, mk in ((0.1, AQUA, "o"), (1.0, BLUE, "s")):
    for key, ls, nm in (("bcr_vs_boomless", "-", "vs Boomless"), ("bcr_vs_subsonic", "--", "vs subsonic")):
        y = [b[(b.cap == o) & (b.capture == c)][key].iloc[0] for o in order]
        ax.plot(xx, y, color=col if key == "bcr_vs_boomless" else (ORANGE if c == 1.0 else YELLOW),
                ls=ls, marker=mk, ms=6, lw=2, mec=BG, label=f"{nm}, {int(c * 100)} % capture")
ax.set_xticks(xx)
ax.set_xticklabels(labs)
ax.set_yscale("log")
ax.set_ylabel("trunk B/C (base cost)")
ax.set_title("(d) M2.0 spine B/C vs lid cap (break-even = 1)", loc="left")
ax.legend(frameon=False, fontsize=7.2, loc="upper left")
clean(ax)

fig.savefig(os.path.join(HERE, "lid_fatigue.png"), dpi=150, bbox_inches="tight", facecolor=BG)
print("ok")
