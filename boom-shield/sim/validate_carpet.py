"""Validation of off-track propagation (ground carpet) and of the matching-circle start.

1. Concorde (CONCORDE preset: M 2.0 at 16.8 km, standard atmosphere, no wind),
   unshielded: ground signature, peak overpressure and PLdB (raw and Gate-0 rise-time
   proxy) versus lateral distance, and the lateral cut-off (carpet half-width), compared
   with published carpet extents:
     * NASA Dryden fact sheet FS-2003-11-016 "Sonic Booms": carpet width ~1 mile per
       1000 ft of altitude (half-width 0.5 mi/kft), independent of the aircraft;
     * Maglieri, Acoustics Today 11(2), 2015 (after NASA SP-2014-622, Maglieri et al.,
       "Sonic Boom: Six Decades of Research"): ray diagram for a Concorde-type aircraft
       at 18 km, primary-carpet edge at ~40 km (25 nmi) lateral.
   Lateral cut-offs are also printed for 52 000 ft, 18 km and 60 000 ft.
2. Linear-inside-R_m error of the 3-D pipeline (Gate-0 geometry: Concorde M 2.0 at
   20.25 km, shield top 20 km, r0 = 1.5 L): (a) direct nonlinear propagation from the
   Whitham near field at r0 versus (b) linear uniform-air near field on a matching
   circle R_m = 1, 2, 4 km then ``propagate_from_circle`` (no correction) and (c) the
   same with the back-propagation aging correction (``age_from = r0``), at
   phi = 0, 20, 40, 50 deg.  A "shield stand-in" (x0.25 and 2 ms Gaussian smoothing
   applied at 20 km; direct: inserted in the nonlinear march; circle: applied to the
   linear near field at the shield's slant distance, then post-shield aging from
   there, and additionally with the nonlinearly aged incident ``aged_near_field``)
   checks the corrections on a wave that the shield has changed.

Outputs results/validation_carpet.png and results/validation_carpet.csv.
Usage:  python boom-shield/sim/validate_carpet.py [--procs 4] [--fs 50000]
"""
from __future__ import annotations

import argparse
import csv
import dataclasses
import math
import os
import sys
import time
from multiprocessing import Pool

import numpy as np
from scipy.ndimage import gaussian_filter1d

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)

from boomsim.propagate import (aged_near_field, ground_carpet,  # noqa: E402
                               ground_signature, lateral_cutoff, propagate_from_circle)
from boomsim.proxy import scored_metrics  # noqa: E402
from boomsim.source import CONCORDE, near_field  # noqa: E402

OUT = os.path.join(HERE, "results")
FT = 0.3048
MILE = 1609.344
H_S = 20000.0
AC20 = dataclasses.replace(CONCORDE, mach=2.0, altitude=H_S + 250.0)   # Gate-0 geometry
R0 = 1.5 * CONCORDE.length
CARPET_DEG = [0, 10, 20, 30, 40, 45, 50, 52, 53.5, 54.5, 54.85, 56]
LZ_DEG = [0, 20, 40, 50]
RMS = [1000.0, 2000.0, 4000.0]
SHIELD_GAIN, SHIELD_SIGMA = 0.25, 2e-3


def shield_standin(t, p):
    fs = 1.0 / (t[1] - t[0])
    return t, SHIELD_GAIN * gaussian_filter1d(p, SHIELD_SIGMA * fs, mode="constant")


def lz_job(job):
    """One linear-zone run -> metrics dict (+ decimated signature)."""
    phi = math.radians(job["deg"])
    fs = job["fs"]
    ac = AC20
    if job["kind"] == "direct":
        ins = (H_S, shield_standin) if job["shield"] else None
        t, p = ground_signature(ac, r0=R0, fs=fs, phi=phi, insert=ins)
    else:
        Rm = job["R_m"]
        if job["shield"]:
            r_sh = (ac.altitude - H_S) / math.cos(phi)       # slant distance to shield top
            if job["kind"] == "circle_preaged":                # nonlinear incident at the shield
                tc, pc = aged_near_field(ac, r_sh, fs, phi=phi, r0=R0)
            else:
                tc, pc = near_field(ac, r_sh, fs, phi=phi)
            tc, pc = shield_standin(tc, pc)
            pc = pc * math.sqrt(r_sh / Rm)                   # linear uniform air to R_m
            age = r_sh
        else:
            tc, pc = near_field(ac, Rm, fs, phi=phi)
            age = R0
        info = {}
        t, p, _ = propagate_from_circle(ac, phi, tc, pc, Rm, z_ref=H_S,
                                        age_from=None if job["kind"] == "circle_linear" else age,
                                        info=info)
        job = {**job, "wkb_factor": info["factor"], "z_start": info["z_start"]}
    m = scored_metrics(t, p)
    i0 = int(np.argmax(np.abs(p) > 1e-3 * np.abs(p).max()))
    sl = slice(max(0, i0 - int(0.02 * fs)), i0 + int(0.40 * fs), 5)
    return {**job, "pldb": m["pldb"], "pldb_proxy": m["pldb_proxy"], "peak_pa": m["peak_pa"],
            "rise_ms": m["rise_time_10_90_ms"], "duration_ms": m["duration_ms"],
            "csel": m["csel"], "sig": (t[sl] - t[i0], p[sl])}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--procs", type=int, default=4)
    ap.add_argument("--fs", type=float, default=50e3)
    args = ap.parse_args()
    os.makedirs(OUT, exist_ok=True)
    t_start = time.time()

    # ---------------------------------------------------------------- 1. carpet
    phis = [math.radians(d) for d in CARPET_DEG]
    car = ground_carpet(CONCORDE, phis, fs=args.fs, procs=args.procs)
    cp = car["carpet"]
    print(f"Concorde M{CONCORDE.mach} at {CONCORDE.altitude:.0f} m: phi_c = "
          f"{math.degrees(cp['phi_cutoff']):.2f} deg, half-width {cp['half_width_m'] / 1e3:.1f} km")
    cut_rows = []
    for name, z in (("52 000 ft (FS-016 cond.)", 52000 * FT), ("CONCORDE preset 16.8 km", CONCORDE.altitude),
                    ("18 km (Maglieri 2015 ex.)", 18000.0), ("60 000 ft", 60000 * FT)):
        ac = dataclasses.replace(CONCORDE, altitude=z)
        pc, yc = lateral_cutoff(z, ac.V)
        rule = 0.5 * MILE * z / FT / 1000.0
        cut_rows.append(dict(case=name, altitude_m=z, phi_c_deg=math.degrees(pc), half_width_km=yc / 1e3,
                             rule_of_thumb_km=rule / 1e3))
        print(f"  cut-off {name:28s}: phi_c {math.degrees(pc):5.2f} deg  y_c {yc / 1e3:5.1f} km"
              f"  (1 mi/kft rule: {rule / 1e3:5.1f} km)")
    for i, d in enumerate(CARPET_DEG):
        print(f"  phi {d:5.2f}  y_g {car['y_g'][i] / 1e3:6.2f} km  peak {car['peak_pa'][i]:6.1f} Pa  "
              f"PLdB {car['pldb'][i]:6.1f} / proxy {car['pldb_proxy'][i]:6.1f}  rise {car['rise_ms'][i]:5.2f} ms")

    # ---------------------------------------------------------------- 2. linear zone
    jobs = []
    for shield in (False, True):
        degs = LZ_DEG if not shield else [0, 40]
        for d in degs:
            jobs.append(dict(kind="direct", deg=d, R_m=float("nan"), shield=shield, fs=args.fs))
            for Rm in RMS:
                kinds = ("circle_linear", "circle_aged") + (("circle_preaged",) if shield else ())
                for kind in kinds:
                    jobs.append(dict(kind=kind, deg=d, R_m=Rm, shield=shield, fs=args.fs))
    with Pool(args.procs) as pool:
        lz = pool.map(lz_job, jobs)
    ref = {(r["deg"], r["shield"]): r for r in lz if r["kind"] == "direct"}
    for r in lz:
        u = ref[(r["deg"], r["shield"])]
        r["dpldb"] = r["pldb"] - u["pldb"]
        r["dpldb_proxy"] = r["pldb_proxy"] - u["pldb_proxy"]
        r["peak_ratio"] = r["peak_pa"] / u["peak_pa"]
    print("\nlinear-inside-R_m error (circle - direct), Concorde M2 at 20.25 km:")
    for r in lz:
        if r["kind"] != "direct":
            print(f"  {'shield-standin' if r['shield'] else 'unshielded':14s} phi {r['deg']:3d}  "
                  f"R_m {r['R_m'] / 1e3:3.0f} km  {r['kind']:13s}  dPLdB {r['dpldb']:+6.2f}  "
                  f"proxy {r['dpldb_proxy']:+6.2f}  peak x{r['peak_ratio']:.3f}")
        else:
            print(f"  {'shield-standin' if r['shield'] else 'unshielded':14s} phi {r['deg']:3d}  direct: "
                  f"peak {r['peak_pa']:.1f} Pa  PLdB {r['pldb']:.2f} / {r['pldb_proxy']:.2f}")

    write_csv(car, cut_rows, lz, os.path.join(OUT, "validation_carpet.csv"))
    plot(car, lz, os.path.join(OUT, "validation_carpet.png"))
    print(f"\nwall time {time.time() - t_start:.0f} s")


def write_csv(car, cut_rows, lz, path):
    with open(path, "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["section", "case", "phi_deg", "R_m", "y_g_m", "delay_s", "peak_pa", "pldb", "pldb_proxy",
                    "rise_ms", "duration_ms", "csel", "dpldb", "dpldb_proxy", "peak_ratio", "note"])
        for i, f in enumerate(car["phi"]):
            m = car["metrics"][i]
            w.writerow(["carpet", "Concorde 16.8 km M2", f"{math.degrees(f):.3f}", "", f"{car['y_g'][i]:.1f}",
                        f"{car['delay'][i]:.3f}", f"{car['peak_pa'][i]:.3f}", f"{car['pldb'][i]:.3f}",
                        f"{car['pldb_proxy'][i]:.3f}", f"{car['rise_ms'][i]:.3f}",
                        f"{m['duration_ms']:.2f}" if m else "", f"{car['csel'][i]:.3f}", "", "", "",
                        "cut off" if car["cutoff"][i] else ""])
        cp = car["carpet"]
        w.writerow(["carpet_summary", "Concorde 16.8 km M2", f"{math.degrees(cp['phi_cutoff']):.3f}", "",
                    f"{cp['half_width_m']:.1f}", "", f"{cp['undertrack_peak_pa']:.3f}",
                    f"{cp['undertrack_pldb']:.3f}", f"{cp['undertrack_pldb_proxy']:.3f}", "", "", "", "", "", "",
                    f"max PLdB {cp['max_pldb']:.2f} / proxy {cp['max_pldb_proxy']:.2f}; y_g = half-width"])
        for r in cut_rows:
            w.writerow(["cutoff", r["case"], f"{r['phi_c_deg']:.3f}", "", f"{1e3 * r['half_width_km']:.1f}",
                        "", "", "", "", "", "", "", "", "", "",
                        f"1 mi per 1000 ft rule half-width {r['rule_of_thumb_km']:.1f} km"])
        for r in lz:
            w.writerow(["linear_zone", ("shield-standin " if r["shield"] else "unshielded ") + r["kind"],
                        r["deg"], "" if r["kind"] == "direct" else f"{r['R_m']:.0f}", "", "",
                        f"{r['peak_pa']:.3f}", f"{r['pldb']:.3f}", f"{r['pldb_proxy']:.3f}",
                        f"{r['rise_ms']:.3f}", f"{r['duration_ms']:.2f}", f"{r['csel']:.3f}",
                        f"{r['dpldb']:+.3f}", f"{r['dpldb_proxy']:+.3f}", f"{r['peak_ratio']:.4f}",
                        "Concorde M2 20.25 km" + (f"; start factor {r['wkb_factor']:.4f}" if "wkb_factor" in r else "")])


def plot(car, lz, path):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    ink, ink2, muted, grid, surf = "#0b0b0b", "#52514e", "#8a8984", "#e4e3df", "#fcfcfb"
    cat = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100"]
    plt.rcParams.update({"font.size": 8.5, "axes.edgecolor": muted, "axes.labelcolor": ink,
                         "xtick.color": ink2, "ytick.color": ink2, "axes.titlecolor": ink,
                         "legend.frameon": False})
    fig, ax = plt.subplots(2, 3, figsize=(15, 8.6))
    fig.patch.set_facecolor(surf)
    for a in ax.flat:
        a.set_facecolor(surf)
        a.grid(True, color=grid, lw=0.6)
        for s in ("top", "right"):
            a.spines[s].set_visible(False)
    ok = ~car["cutoff"]
    y = car["y_g"][ok] / 1e3
    yc = car["carpet"]["half_width_m"] / 1e3
    rule = 0.5 * MILE * CONCORDE.altitude / FT / 1000.0 / 1e3

    def edges(a):
        a.axvline(yc, color=ink2, lw=1.0, ls="--")
        a.text(yc, a.get_ylim()[1], f" ray cut-off {yc:.1f} km", color=ink2, va="top", ha="right",
               rotation=90, fontsize=7.5)
        a.axvline(rule, color=muted, lw=1.0, ls=":")
        a.text(rule, a.get_ylim()[0], f" 1 mi/1000 ft rule {rule:.0f} km", color=muted, va="bottom",
               ha="left", rotation=90, fontsize=7.5)

    a = ax[0, 0]
    a.plot(y, car["peak_pa"][ok], color=cat[0], lw=2, marker="o", ms=4)
    a.set_ylim(0, 1.15 * np.nanmax(car["peak_pa"]))
    a.set_xlim(0, 50)
    edges(a)
    a.set(xlabel="ground lateral offset y_g [km]", ylabel="peak overpressure [Pa]",
          title="Concorde M2.0, 16.8 km: peak overpressure across the carpet")

    a = ax[0, 1]
    a.plot(y, car["pldb"][ok], color=cat[0], lw=2, marker="o", ms=4, label="PLdB raw")
    a.plot(y, car["pldb_proxy"][ok], color=cat[1], lw=2, marker="s", ms=4, label="PLdB rise-time proxy (σ = 1 ms)")
    a.set_xlim(0, 50)
    a.set_ylim(80, 115)
    edges(a)
    a.legend(loc="lower left")
    a.set(xlabel="ground lateral offset y_g [km]", ylabel="PLdB", title="Loudness across the carpet")

    a = ax[0, 2]
    show = [0, 30, 45, 54.5]
    for k, d in enumerate(show):
        i = CARPET_DEG.index(d)
        t, p = car["t"][i], car["p"][i]
        i0 = int(np.argmax(np.abs(p) > 1e-3 * np.abs(p).max()))
        sl = slice(max(0, i0 - 1000), i0 + 20000, 5)
        a.plot(1e3 * (t[sl] - t[i0]), p[sl], color=cat[k], lw=1.4,
               label=f"φ = {d:g}°, y_g = {car['y_g'][i] / 1e3:.1f} km")
    a.legend(loc="upper right")
    a.set(xlabel="time from first arrival [ms]", ylabel="overpressure [Pa]",
          title="Ground signatures (×1.9 reflection)")

    for col, shield in ((0, False), (1, True)):
        a = ax[1, col]
        degs = sorted({r["deg"] for r in lz if r["shield"] == shield})
        for k, d in enumerate(degs):
            for kind, ls, mk, lab in (("circle_linear", "-", "o", "linear inside R_m"),
                                      ("circle_aged", "--", "s", "post-shield aging (age_from)"),
                                      ("circle_preaged", ":", "^", "+ aged incident at shield")):
                rr = sorted([r for r in lz if r["shield"] == shield and r["deg"] == d and r["kind"] == kind],
                            key=lambda r: r["R_m"])
                if not rr:
                    continue
                a.plot([r["R_m"] / 1e3 for r in rr], [r["dpldb"] for r in rr], color=cat[k], ls=ls,
                       marker=mk, ms=5, lw=1.6)
        a.axhline(0, color=ink2, lw=0.8)
        a.set_xticks([1, 2, 4])
        from matplotlib.lines import Line2D
        hs = [Line2D([], [], color=cat[k], lw=2, label=f"φ = {d}°") for k, d in enumerate(degs)]
        styles = (("-", "o", "linear inside R_m"), ("--", "s", "post-shield aging (age_from)"))
        if shield:
            styles += ((":", "^", "+ aged incident at shield"),)
        hs += [Line2D([], [], color=ink2, ls=ls, marker=mk, ms=5, label=lab) for ls, mk, lab in styles]
        a.legend(handles=hs, loc="upper left", fontsize=7.5)
        a.set(xlabel="matching-circle radius R_m [km]", ylabel="ΔPLdB raw (circle − direct)",
              title=("Linear near zone error, unshielded (Concorde M2, 20.25 km)" if not shield else
                     "Same with a shield stand-in at 20 km (×0.25, 2 ms smoothing)"))

    a = ax[1, 2]
    d0 = [r for r in lz if not r["shield"] and r["deg"] == 0]
    sel = [("direct", None, ink, "direct nonlinear from r0 = 1.5 L"),
           ("circle_linear", 4000.0, cat[1], "circle R_m = 4 km, linear inside"),
           ("circle_aged", 4000.0, cat[0], "circle R_m = 4 km, aging correction")]
    for kind, Rm, c, lab in sel:
        r = [x for x in d0 if x["kind"] == kind and (Rm is None or x["R_m"] == Rm)][0]
        t, p = r["sig"]
        a.plot(1e3 * t, p, color=c, lw=1.4 if kind != "circle_aged" else 1.0,
               ls="-" if kind != "circle_aged" else "--", label=f"{lab}: {r['pldb']:.1f} PLdB")
    a.legend(loc="upper right", fontsize=7.5)
    a.set(xlabel="time from first arrival [ms]", ylabel="overpressure [Pa]",
          title="Undertrack ground signature, Concorde M2 at 20.25 km")
    fig.tight_layout()
    fig.savefig(path, dpi=130, facecolor=surf)


if __name__ == "__main__":
    main()
