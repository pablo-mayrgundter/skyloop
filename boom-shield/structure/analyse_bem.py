"""Summarise results/raft_bem.csv -> results/raft_bem_summary.csv and figures/raft_vs_slab.png.

Metrics per config (dB re unshielded incident, power averages):
  under  : phi 0-3 deg on the R = 3 km circle
  circle : phi 0-50 deg
  line   : near-field line 100 m below the shield, |y| <= 300 m
bands: low = 2-3 Hz, audible-relevant = 5-30 Hz (3-D study: > 5 Hz sets PLdB), per frequency too.
"""
import csv
import os
from collections import defaultdict

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
RES = os.path.join(HERE, "results")
FIG = os.path.join(HERE, "figures")
os.makedirs(FIG, exist_ok=True)


def load():
    d = defaultdict(lambda: defaultdict(lambda: defaultdict(list)))
    with open(os.path.join(RES, "raft_bem.csv")) as fh:
        for r in csv.DictReader(fh):
            x = float(r["phi_or_y"])
            if r["where"] == "circle":
                keys = ["circle"] + (["under"] if x <= 3 else [])
            else:
                keys = ["line"]
            for k in keys:
                d[r["config"]][float(r["f"])][k].append(float(r["dB_re_incident"]))
    return d


def pmean(v):
    v = np.asarray(v, float)
    return 10 * np.log10(np.mean(10 ** (v / 10)))


def summary(d):
    rows = []
    for cfg, byf in d.items():
        fs = sorted(byf)
        row = dict(config=cfg)
        for k in ("under", "circle", "line"):
            per = {f: pmean(byf[f][k]) for f in fs}
            for f in fs:
                row[f"{k}_{f:g}Hz"] = round(per[f], 1)
            lo = [per[f] for f in fs if f <= 3]
            hi = [per[f] for f in fs if f >= 5]
            row[f"{k}_2-3Hz"] = round(pmean(lo), 1) if lo else None
            row[f"{k}_5-30Hz"] = round(pmean(hi), 1) if hi else None
        rows.append(row)
    keys = sorted({k for r in rows for k in r}, key=lambda k: (k != "config", k))
    with open(os.path.join(RES, "raft_bem_summary.csv"), "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=keys)
        w.writeheader()
        for r in rows:
            w.writerow(r)
    return rows


def plot(d):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    ink, ink2, grid = "#0b0b0b", "#52514e", "#e4e3df"
    s1, s2, s3, s4 = "#2a78d6", "#eb6834", "#1baf7a", "#6250d6"
    panels = [
        ("Raft D 30 m: hoop-stretch stiffness Et (Δp 300 Pa)", [
            ("slab_H2_30", "slab H₂ 30 m", ink2, "-"),
            ("raft_D30_dp300", "raft, Et 1e5 N/m (Et/R² 444)", s1, "-"),
            ("raft_D30_dp300_Et3e5", "raft, Et 3e5 (1,333)", s2, "-"),
            ("raft_D30_dp300_Et1e6", "raft, Et 1e6 (4,444)", s4, "-")]),
        ("Raft geometry: air gaps and diameter", [
            ("slab_H2_30", "slab H₂ 30 m", ink2, "-"),
            ("raft_D30_dp300_gap2", "raft D30, 2 m air gaps", s2, "-"),
            ("raft_D30_dp300_gap5", "raft D30, 5 m air gaps", s4, "-"),
            ("raft_D60_dp300", "raft D60, touching", s1, "-")]),
        ("Floating lid (heavy top, light floor) and ties", [
            ("slab_H2_30", "slab H₂ 30 m, no ties", ink2, "-"),
            ("slab_H2_30_heavytop_K10", "floating lid, K 10 Pa/m", s1, "-"),
            ("slab_H2_30_light_K10", "light lid 0.15, K 10", s3, "-"),
            ("slab_H2_30_K300", "lift-closed, K 300 Pa/m", s4, "-")]),
    ]
    fig, axs = plt.subplots(1, 3, figsize=(15, 5.6), sharey=True)
    for ax, (title, series) in zip(axs, panels):
        for key, lab, col, ls in series:
            if key not in d:
                continue
            fs = sorted(d[key])
            y = [pmean(d[key][f]["circle"]) for f in fs]
            ax.plot(fs, y, ls, color=col, lw=2, marker="o", ms=5, label=lab)
        ax.set_xscale("log")
        ax.set_xticks([2, 3, 5, 10, 20, 30])
        ax.set_xticklabels(["2", "3", "5", "10", "20", "30"])
        ax.set_title(title, fontsize=11, color=ink, loc="left")
        ax.set_xlabel("frequency [Hz]", color=ink2)
        ax.grid(True, color=grid, lw=0.8)
        for sp_ in ("top", "right"):
            ax.spines[sp_].set_visible(False)
        ax.xaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())
        ax.legend(frameon=False, fontsize=9, loc="upper center", bbox_to_anchor=(0.5, -0.16), ncol=2)
    axs[0].set_ylabel("transmitted, dB re unshielded\n(power mean φ 0–50°, R 3 km)", color=ink2)
    axs[0].set_ylim(-45, 0)
    fig.suptitle("Cross-plane BEM, Concorde M2, h 250 m, W 1 km, H₂ at 20 km (lift-closed faces unless noted)",
                 fontsize=11, color=ink, x=0.01, ha="left")
    fig.tight_layout()
    fig.savefig(os.path.join(FIG, "raft_vs_slab.png"), dpi=130)


if __name__ == "__main__":
    d = load()
    rows = summary(d)
    for r in rows:
        print(f"{r['config']:30s} under 2-3 {r.get('under_2-3Hz')}  5-30 {r.get('under_5-30Hz')} | circle 2-3 "
              f"{r.get('circle_2-3Hz')}  5-30 {r.get('circle_5-30Hz')} | line 5-30 {r.get('line_5-30Hz')}")
    plot(d)
