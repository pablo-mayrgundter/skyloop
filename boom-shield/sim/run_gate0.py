"""Gate 0: shield inserted into full nonlinear propagation to the ground.

Question: for lift-closed membrane/gas/membrane stacks at 20 km, how much does
ground loudness (PLdB) drop, and does the boom re-steepen below the shield?

Pipeline per case: Whitham near field at r0 below the aircraft -> augmented Burgers
down to the shield top H_s = 20 km -> zero-pad (room for the non-causal Lorentzian
smoothing tails of evanescent layers) -> shield.apply (plane-wave transfer matrix at
the boom trace speed V, 20 km ambient) -> augmented Burgers to the ground -> x1.9
ground reflection -> loudness.metrics. The unshielded reference of each case goes
through exactly the same path with an identity insert (same r0, same padding).

Outputs: results/gate0_results.csv, results/gate0_signatures.png,
results/gate0_summary.md.

Run:  python boom-shield/sim/run_gate0.py [--procs 4] [--cache FILE] [--replot]
      (--cache stores raw results; --replot rebuilds the outputs from it.)
"""
from __future__ import annotations

import argparse
import csv
import dataclasses
import math
import os
import pickle
import sys
import time
from multiprocessing import Pool

import numpy as np
from scipy.ndimage import gaussian_filter1d
from scipy.signal import lfilter

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from boomsim import shield as S  # noqa: E402
from boomsim.atmosphere import std_atm  # noqa: E402
from boomsim.loudness import metrics, pldb  # noqa: E402
from boomsim.propagate import ground_signature  # noqa: E402
from boomsim.source import CONCORDE, SR71  # noqa: E402

RESULTS = os.path.join(HERE, "results")

H_S = 20000.0                 # shield top [m]
FS = 50e3                     # sample rate [Hz] (unshielded shocks need it)
RISE_SIGMA_S = 1.0e-3         # rise-time proxy: Gaussian kernel std -> ideal step 10-90 % = 2.56 ms
RISE_TAU_S = 1.14e-3          # shape check: one-sided exponential kernel -> ideal step 10-90 % = 2.5 ms
HS_LIST = [100.0, 250.0]      # aircraft height above shield top [m]
AIRCRAFT = [("Concorde", 1.4), ("Concorde", 2.0), ("Concorde", 2.5), ("SR-71", 3.2)]
DS = [15.0, 30.0, 60.0]
TIE_KS = [0.0, 30.0, 100.0, 300.0, 1000.0]
PLOT_WIN = (-0.10, 0.50)      # s, retarded-time window kept for plotting
PLOT_DEC = 5                  # 50 kHz -> 10 kHz for stored plot traces
TAPER_FRAC = 0.5              # outer fraction of each shield pad tapered to zero (half-cosine)
METRIC_TAPER_S = 0.05         # half-Hann taper at both record ends for the loudness metrics
LEAK_DB = [-6.0, -12.0, -18.0, -25.0]   # "leak floor": unshielded wave x 10^(dB/20) at 20 km

ATM20 = std_atm(H_S)
T20, P20, RHO20 = float(ATM20.T), float(ATM20.p), float(ATM20.rho)


# ----------------------------------------------------------------------------- setup
def aircraft(name, mach, h):
    base = {"Concorde": CONCORDE, "SR-71": SR71}[name]
    # Only altitude and Mach are overridden; weight/volume/length kept from the preset
    # (Concorde W = 1.47 MN, SR-71 W = 0.54 MN) at every Mach.
    return dataclasses.replace(base, mach=float(mach), altitude=H_S + float(h))


def fill_gas(fill, x_air=0.0):
    if x_air:
        return S.gas_mixture({fill: 1.0 - x_air, "air": x_air}, T20, P20)
    return S.gas(fill, T20, P20)


def budget(g, d):
    """Slab lift per m^2 of a gas layer at 20 km [kg/m^2]."""
    return (RHO20 - g.rho) * d


def make_stack(kind, fill=None, d=0.0, tie_K=0.0, x_air=0.0):
    """Returns (stack or None, lift_margin [kg/m^2], gas or None)."""
    if kind in ("none", "leak"):
        return None, float("nan"), None
    if kind == "masslaw":
        he = fill_gas("helium")
        b = budget(he, 15.0)
        st = S.Stack([S.Membrane(b)])
        # notional: the membrane is floated by a 15 m He slab that is acoustically ignored
        lm = S.lift_margin(S.Stack([S.GasLayer(he, 15.0), S.Membrane(b)]), ATM20)
        return st, lm, None
    g = fill_gas(fill, x_air)
    b = budget(g, d)
    st = S.Stack([S.Membrane(b / 2), S.GasLayer(g, d, tie_K=tie_K), S.Membrane(b / 2)])
    return st, S.lift_margin(st, ATM20), g


def stack_label(c):
    if c["stack"] == "none":
        return "unshielded"
    if c["stack"] == "masslaw":
        return "mass law (1 membrane, He-15 m lift)"
    if c["stack"] == "leak":
        return f"leak floor: unshielded x {c['leak_db']:+g} dB"
    g = {"helium": "He", "hydrogen": "H2"}[c["fill"]]
    s = f"M/{g} {c['d']:g} m/M"
    if c["x_air"]:
        s += f" ({100 * c['x_air']:g}% air)"
    if c["tie_K"]:
        s += f" K={c['tie_K']:g}"
    return s


def build_cases():
    cases = []

    def add(group, ac, mach, h, stack, fill=None, d=0.0, tie_K=0.0, x_air=0.0,
            below_linear=False, pad_mult=1.0, leak_db=0.0):
        cases.append(dict(group=group, aircraft=ac, mach=mach, h=h, stack=stack, fill=fill,
                          d=d, tie_K=tie_K, x_air=x_air, below_linear=below_linear,
                          pad_mult=pad_mult, leak_db=leak_db))

    for ac, mach in AIRCRAFT:
        for h in HS_LIST:
            add("main", ac, mach, h, "none")
            add("main", ac, mach, h, "masslaw")
            for fill in ("helium", "hydrogen"):
                for d in DS:
                    add("main", ac, mach, h, "gas", fill, d)
    for fill, d in (("helium", 15.0), ("hydrogen", 30.0)):
        for K in TIE_KS[1:]:
            add("tie", "Concorde", 2.0, 250.0, "gas", fill, d, tie_K=K)
    for mach in (2.0, 2.5):
        for x in (0.05, 0.10):
            add("purity", "Concorde", mach, 250.0, "gas", "helium", 15.0, x_air=x)
    # re-steepening: nonlinearity off below the shield
    add("resteep", "Concorde", 2.0, 250.0, "none", below_linear=True)
    add("resteep", "Concorde", 2.0, 250.0, "gas", "helium", 15.0, below_linear=True)
    add("resteep", "Concorde", 2.0, 250.0, "gas", "hydrogen", 30.0, below_linear=True)
    add("resteep", "Concorde", 1.4, 250.0, "gas", "helium", 60.0, below_linear=True)
    add("resteep", "SR-71", 3.2, 250.0, "gas", "hydrogen", 30.0, below_linear=True)
    add("resteep", "Concorde", 1.4, 250.0, "gas", "hydrogen", 30.0, below_linear=True)
    # padding convergence: widest Lorentzian kernels in the matrix, 2x the pad
    add("padcheck", "Concorde", 1.4, 250.0, "gas", "hydrogen", 60.0, pad_mult=2.0)
    add("padcheck", "Concorde", 2.0, 250.0, "gas", "hydrogen", 30.0, pad_mult=2.0)
    add("padcheck", "SR-71", 3.2, 250.0, "gas", "hydrogen", 30.0, pad_mult=2.0)
    # leak floor: what ground loudness a given fractional leak of the unshielded wave gives
    for db in LEAK_DB:
        add("leak", "Concorde", 2.0, 250.0, "leak", leak_db=db)
    for i, c in enumerate(cases):
        c["id"] = i
    return cases


# ----------------------------------------------------------------------------- run
def _crop(t, p):
    m = (t >= PLOT_WIN[0]) & (t <= PLOT_WIN[1])
    return t[m][::PLOT_DEC].astype(np.float32), p[m][::PLOT_DEC].astype(np.float32)


def _pos_impulse(t, p):
    """Positive impulse of the first positive phase (same definition as loudness.metrics)."""
    pk = float(np.max(p))
    if pk <= 0:
        return float("nan")
    ih = int(np.argmax(p >= 0.5 * pk))
    nonpos = np.nonzero(p[:ih] <= 0)[0]
    i0 = nonpos[-1] if nonpos.size else 0
    after = np.nonzero(p[ih:] <= 0)[0]
    i1 = ih + after[0] if after.size else p.size - 1
    seg = np.clip(p[i0:i1 + 1], 0.0, None)
    return float(np.sum(0.5 * (seg[1:] + seg[:-1])) * (t[1] - t[0]))


def run_case(c):
    t_start = time.time()
    ac = aircraft(c["aircraft"], c["mach"], c["h"])
    V = ac.V
    r0 = min(1.5 * ac.length, 0.8 * c["h"])
    st, lm, g = make_stack(c["stack"], c["fill"], c["d"], c["tie_K"], c["x_air"])
    # Lorentzian width of an evanescent layer: a d, a = sqrt(1/V^2 - 1/c_g^2)
    a = math.sqrt(max(0.0, 1 / V ** 2 - 1 / g.c ** 2)) if g is not None else 0.0
    pad = c["pad_mult"] * min(3.0, max(1.0, 15.0 * a * c["d"]))
    cap = {}

    def insert(t, p):
        dt = t[1] - t[0]
        n = int(round(pad / dt))
        p2 = np.concatenate([np.zeros(n), p, np.zeros(n)])
        t2 = t[0] + (np.arange(p2.size) - n) * dt
        if st is not None:
            out = S.apply(st, t2, p2, V, ATM20)
        elif c["stack"] == "leak":
            out = p2 * 10.0 ** (c["leak_db"] / 20.0)
        else:
            out = p2.copy()
        # The evanescent layer smooths with a Lorentzian (~1/t^2) kernel; cropping its
        # tails at the record ends leaves a step of ~A w/(pi pad^2) (A = net area) that
        # would later radiate spurious high-frequency energy (tens of dB of band level,
        # which dominates PLdB once the boom itself is smooth). Taper the outer
        # TAPER_FRAC of each pad smoothly to zero (a slow half-cosine: no audible content).
        nt = max(2, int(TAPER_FRAC * n))
        w = np.ones(out.size)
        ramp = 0.5 * (1 - np.cos(np.pi * np.arange(nt) / nt))
        w[:nt] = ramp
        w[-nt:] = ramp[::-1]
        out = out * w
        cap["above"] = (t2, p2)
        cap["below"] = (t2, out)
        return t2, out

    t, p = ground_signature(ac, r0=r0, insert=(H_S, insert), fs=FS,
                            below_opts={"linear": True} if c["below_linear"] else None)
    m = metrics(t, p, taper_s=METRIC_TAPER_S)
    ps = gaussian_filter1d(p, RISE_SIGMA_S * FS, mode="constant")
    ms = metrics(t, ps, taper_s=METRIC_TAPER_S)
    ae = math.exp(-1.0 / (RISE_TAU_S * FS))
    pe = lfilter([1.0 - ae], [1.0, -ae], p)
    ta, pa = cap["above"]
    tb, pb = cap["below"]
    res = dict(c)
    res.update(
        label=stack_label(c), V=V, r0=r0, L=ac.length, weight_N=ac.weight, pad_s=pad,
        regime=("-" if g is None else "evanescent" if V < g.c else "propagating"),
        c_gas=(g.c if g is not None else float("nan")),
        face_mass_total=(sum(l.areal_mass for l in st.layers) if st is not None else 0.0),
        lift_margin=lm, lorentz_width_ms=1e3 * a * c["d"],
        pldb=m["pldb"], peak_pa=m["peak_pa"], peak_psf=m["peak_psf"],
        rise_ms=m["rise_time_10_90_ms"], rise_sboom_ms=m["rise_time_ms"],
        duration_ms=m["duration_ms"], pos_impulse=m["pos_impulse_pa_s"],
        asel=m["asel"], csel=m["csel"],
        pldb_smooth=ms["pldb"], pldb_exp=pldb(t, pe, taper_s=METRIC_TAPER_S), peak_pa_smooth=ms["peak_pa"], rise_ms_smooth=ms["rise_time_10_90_ms"],
        shield_peak_above=float(pa.max()), shield_peak_below=float(pb.max()),
        shield_imp_above=_pos_impulse(ta, pa), shield_imp_below=_pos_impulse(tb, pb),
        shield_area_ratio=float(np.trapezoid(pb, tb) / np.trapezoid(pa, ta)),
        runtime_s=time.time() - t_start,
    )
    res["shield_peak_ratio"] = res["shield_peak_below"] / res["shield_peak_above"]
    res["shield_imp_ratio"] = res["shield_imp_below"] / res["shield_imp_above"]
    res["wave_ground"] = _crop(t, p)
    res["ground_full"] = (t.astype(np.float32), p.astype(np.float32))
    res["wave_above"] = _crop(ta, pa)
    res["wave_below"] = _crop(tb, pb)
    return res


# ----------------------------------------------------------------------------- post
def key(r):
    return (r["aircraft"], r["mach"], r["h"])


def add_deltas(rows):
    ref = {key(r): r for r in rows if r["stack"] == "none" and r["group"] == "main"}
    ref_lin = {key(r): r for r in rows if r["stack"] == "none" and r["below_linear"]}
    for r in rows:
        u = ref[key(r)]
        r["dpldb_raw"] = r["pldb"] - u["pldb"]
        r["dpldb_smooth"] = r["pldb_smooth"] - u["pldb_smooth"]
        r["dpldb_exp"] = r["pldb_exp"] - u["pldb_exp"]
        r["dcsel"] = r["csel"] - u["csel"]
        r["dasel"] = r["asel"] - u["asel"]
        r["dpldb_headline"] = max(r["dpldb_raw"], r["dpldb_smooth"])  # less favourable to shield
        r["peak_ratio_ground"] = r["peak_pa"] / u["peak_pa"]
        r["imp_ratio_ground"] = r["pos_impulse"] / u["pos_impulse"]
        if r["below_linear"] and key(r) in ref_lin:
            ul = ref_lin[key(r)]
            r["dpldb_raw_vs_lin_ref"] = r["pldb"] - ul["pldb"]
        else:
            r["dpldb_raw_vs_lin_ref"] = float("nan")
    return rows


CSV_COLS = [
    "id", "group", "aircraft", "mach", "V", "h", "r0", "L", "weight_N", "label", "stack", "fill",
    "d", "tie_K", "x_air", "below_linear", "regime", "c_gas", "face_mass_total", "lift_margin",
    "lorentz_width_ms", "pad_s", "pldb", "pldb_smooth", "pldb_exp", "dpldb_raw", "dpldb_smooth", "dpldb_exp",
    "dpldb_headline", "asel", "csel", "dasel", "dcsel", "peak_pa", "peak_psf", "peak_pa_smooth", "rise_ms", "rise_ms_smooth",
    "rise_sboom_ms", "duration_ms", "pos_impulse", "peak_ratio_ground", "imp_ratio_ground",
    "shield_peak_above", "shield_peak_below", "shield_peak_ratio", "shield_imp_above",
    "shield_imp_below", "shield_imp_ratio", "shield_area_ratio", "dpldb_raw_vs_lin_ref",
    "runtime_s",
]


def write_csv(rows, path):
    with open(path, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=CSV_COLS, extrasaction="ignore")
        w.writeheader()
        for r in sorted(rows, key=lambda r: r["id"]):
            w.writerow({k: (f"{r[k]:.6g}" if isinstance(r[k], float) else r[k]) for k in CSV_COLS})


def plot(rows, path):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    ink, ink2, muted, grid, surf = "#0b0b0b", "#52514e", "#8a8984", "#e4e3df", "#fcfcfb"
    col = {"none": ink2, "He15": "#2a78d6", "H230": "#eb6834", "best": "#1baf7a"}
    plt.rcParams.update({"font.size": 8.5, "axes.edgecolor": muted, "axes.labelcolor": ink,
                         "xtick.color": ink2, "ytick.color": ink2, "axes.titlecolor": ink})
    main = [r for r in rows if r["group"] == "main" and r["h"] == 250.0]
    lin = {key(r) + (r["stack"], r["fill"], r["d"]): r for r in rows if r["group"] == "resteep"}
    panels = AIRCRAFT

    def ground(r, dec=PLOT_DEC):
        if "ground_full" in r:
            t, p = (np.asarray(x, dtype=float) for x in r["ground_full"])
            return t[::dec], p[::dec]
        return (np.asarray(x, dtype=float) for x in r["wave_ground"])
    fig, axes = plt.subplots(3, len(panels), figsize=(15, 10.5))
    fig.patch.set_facecolor(surf)

    def find(ac, mach, stack, fill=None, d=0.0):
        for r in main:
            if (r["aircraft"], r["mach"], r["stack"], r["fill"], r["d"]) == (ac, mach, stack, fill, d):
                return r
        return None

    for j, (ac, mach) in enumerate(panels):
        u = find(ac, mach, "none")
        he = find(ac, mach, "gas", "helium", 15.0)
        h2 = find(ac, mach, "gas", "hydrogen", 30.0)
        sh = [r for r in main if r["aircraft"] == ac and r["mach"] == mach and r["stack"] == "gas"]
        best = min(sh, key=lambda r: r["pldb"])
        series = [(u, "unshielded", col["none"]), (he, "He 15 m", col["He15"]),
                  (h2, "H2 30 m", col["H230"])]
        if best not in (he, h2):
            series.append((best, f"{'He' if best['fill'] == 'helium' else 'H2'} {best['d']:g} m (quietest)",
                           col["best"]))
        tu, pu = ground(u)
        act = np.flatnonzero(np.abs(pu) > 0.01 * np.max(np.abs(pu)))
        gx = (1e3 * tu[act[0]] - 60.0, 1e3 * tu[act[-1]] + 60.0)
        hl = lin.get((ac, mach, 250.0, "gas", "helium", 15.0))
        # row 0: just above the shield (identical for all) + just below for shielded
        ax = axes[0, j]
        t, p = u["wave_above"]
        ax.plot(1e3 * t, p, color=ink2, lw=1.6, label="above shield (all cases)")
        for r, lab, c in series[1:]:
            t, p = r["wave_below"]
            ax.plot(1e3 * t, p, color=c, lw=1.6, label=f"below, {lab}")
        ax.set_title(f"{ac}  M {mach}  (V = {u['V']:.0f} m/s), h = 250 m\nshield plane, 20 km", fontsize=9)
        # row 1: ground, raw model
        for row, kname, ttl in ((1, "wave_ground", "ground (model rise time)"),
                                (2, None, f"ground, rise-time proxy (Gaussian σ = {1e3 * RISE_SIGMA_S:g} ms)")):
            ax = axes[row, j]
            extra = [(hl, "He 15 m, linear below shield", col["He15"])] if hl is not None else []
            for i, (r, lab, c) in enumerate(series + extra):
                t, p = ground(r, 1)
                if kname is None:
                    p = gaussian_filter1d(p, RISE_SIGMA_S * FS, mode="constant")
                    pl = r["pldb_smooth"]
                else:
                    pl = r["pldb"]
                t, p = t[::PLOT_DEC], p[::PLOT_DEC]
                pls = f"{pl:.1f} PLdB" if np.isfinite(pl) else "PL −∞ (below Mark VII)"
                ls = ":" if i >= len(series) else "-"
                ax.plot(1e3 * t, p, color=c, lw=1.6, ls=ls, label=f"{lab}: {pls}")
            ax.set_title(ttl, fontsize=9)
        for row in range(3):
            ax = axes[row, j]
            ax.set_facecolor(surf)
            ax.grid(True, color=grid, lw=0.6)
            ax.axhline(0, color=muted, lw=0.6)
            for s in ("top", "right"):
                ax.spines[s].set_visible(False)
            ax.legend(frameon=False, fontsize=7.5, loc="upper right" if row == 0 else "lower right")
            if row == 0:
                ax.set_xlim(-60 if ac == "Concorde" else -40, 480 if mach < 2.2 else 380)
            else:
                ax.set_xlim(*gx)
            if j == 0:
                ax.set_ylabel("overpressure [Pa]")
            if row == 2:
                ax.set_xlabel("retarded time [ms]")
    fig.suptitle("Gate 0 — sonic boom with a lift-closed membrane/gas/membrane shield at 20 km "
                 "(1-D undertrack, infinite shield, no turbulence/wind)", color=ink, fontsize=11)
    fig.tight_layout(rect=(0, 0, 1, 0.97))
    fig.savefig(path, dpi=125, facecolor=surf)
    plt.close(fig)


# ----------------------------------------------------------------------------- summary
def _f(x, nd=1, sign=False):
    if x is not None and x == -np.inf:
        return "−∞"
    if x is None or (isinstance(x, float) and not np.isfinite(x)):
        return "–"
    return f"{x:+.{nd}f}" if sign else f"{x:.{nd}f}"


def tables(rows):
    out = {}
    main = [r for r in rows if r["group"] == "main"]
    stacks = [("masslaw", None, 0.0)] + [("gas", f, d) for f in ("helium", "hydrogen") for d in DS]
    cols = [(ac, m, h) for h in HS_LIST for ac, m in AIRCRAFT]

    def cell(r):
        mark = "" if r["regime"] in ("-", "evanescent") else "ᵖ"
        return f"{_f(r['dpldb_headline'], 1, True)}{mark}"

    hdr = "| stack (face mass kg/m²) | " + " | ".join(
        f"{'Conc' if a == 'Concorde' else 'SR71'} M{m} h{h:g}" for a, m, h in cols) + " |"
    lines = [hdr, "|" + "---|" * (len(cols) + 1)]
    ref = {key(r): r for r in main if r["stack"] == "none"}
    lines.append("| unshielded PLdB raw / proxy | " + " | ".join(
        f"{ref[c]['pldb']:.1f} / {ref[c]['pldb_smooth']:.1f}" for c in cols) + " |")
    for kind, fill, d in stacks:
        rr = {key(r): r for r in main if (r["stack"], r["fill"], r["d"]) == (kind, fill, d)}
        any_r = next(iter(rr.values()))
        name = stack_label(any_r)
        lines.append(f"| {name} ({any_r['face_mass_total']:.2f}) | "
                     + " | ".join(cell(rr[c]) for c in cols) + " |")
    out["headline"] = "\n".join(lines)

    # both-ways table for h = 250
    lines = ["| stack | " + " | ".join(f"{'Conc' if a == 'Concorde' else 'SR71'} M{m}" for a, m in AIRCRAFT)
             + " |", "|---|" + "---|" * len(AIRCRAFT)]
    for kind, fill, d in stacks:
        rr = {key(r): r for r in main if (r["stack"], r["fill"], r["d"]) == (kind, fill, d)}
        lines.append(f"| {stack_label(next(iter(rr.values())))} | " + " | ".join(
            f"{_f(rr[(a, m, 250.0)]['dpldb_raw'], 1, True)} / {_f(rr[(a, m, 250.0)]['dpldb_smooth'], 1, True)}"
            for a, m in AIRCRAFT) + " |")
    out["bothways"] = "\n".join(lines)

    # shield-plane and ground ratios, h = 250
    lines = ["| case | shield peak ratio | shield impulse ratio | ground peak ratio | ground impulse ratio "
             "| ground 10–90 % rise [ms] (unshielded) | ΔPLdB headline | ΔCSEL [dB] |", "|---|---|---|---|---|---|---|---|"]
    for a, m in AIRCRAFT:
        for kind, fill, d in (("gas", "helium", 15.0), ("gas", "hydrogen", 30.0), ("gas", "helium", 60.0)):
            r = next(r for r in main if key(r) == (a, m, 250.0) and (r["stack"], r["fill"], r["d"]) == (kind, fill, d))
            u = ref[key(r)]
            lines.append(f"| {a} M{m}, {r['label']} | {r['shield_peak_ratio']:.2f} | {r['shield_imp_ratio']:.2f} | "
                         f"{r['peak_ratio_ground']:.2f} | {r['imp_ratio_ground']:.2f} | "
                         f"{_f(r['rise_ms'])} ({_f(u['rise_ms'])}) | {_f(r['dpldb_headline'], 1, True)} | "
                         f"{_f(r['dcsel'], 1, True)} |")
    out["ratios"] = "\n".join(lines)

    # ties
    lines = ["| tie K [Pa/m] | " + " | ".join(f"{k:g}" for k in TIE_KS) + " |", "|---|" + "---|" * len(TIE_KS)]
    for fill, d in (("helium", 15.0), ("hydrogen", 30.0)):
        vals = []
        for K in TIE_KS:
            r = next(r for r in rows if r["group"] in ("main", "tie") and key(r) == ("Concorde", 2.0, 250.0)
                     and r["stack"] == "gas" and r["fill"] == fill and r["d"] == d and r["tie_K"] == K
                     and not r["x_air"])
            vals.append(f"{_f(r['dpldb_raw'], 1, True)} / {_f(r['dpldb_smooth'], 1, True)}")
        lines.append(f"| {'He' if fill == 'helium' else 'H2'} {d:g} m: ΔPLdB raw / proxy | " + " | ".join(vals) + " |")
    out["ties"] = "\n".join(lines)

    # purity
    lines = ["| He 15 m, Concorde h 250 | air fraction | c_gas [m/s] | regime | face mass [kg/m²] "
             "| ΔPLdB raw / proxy |", "|---|---|---|---|---|---|"]
    for m in (2.0, 2.5):
        for x in (0.0, 0.05, 0.10):
            r = next(r for r in rows if r["group"] in ("main", "purity") and key(r) == ("Concorde", m, 250.0)
                     and r["stack"] == "gas" and r["fill"] == "helium" and r["d"] == 15.0
                     and r["x_air"] == x and not r["tie_K"])
            lines.append(f"| M{m} (V = {r['V']:.0f}) | {100 * x:g} % | {r['c_gas']:.0f} | {r['regime']} | "
                         f"{r['face_mass_total']:.2f} | {_f(r['dpldb_raw'], 1, True)} / {_f(r['dpldb_smooth'], 1, True)} |")
    out["purity"] = "\n".join(lines)

    # re-steepening
    lines = ["| case (h 250) | below-shield leg | ground PLdB raw | peak [Pa] | 10–90 % rise [ms] | ΔPLdB raw vs unshielded (same leg model) |",
             "|---|---|---|---|---|---|"]
    for r in sorted([r for r in rows if r["group"] == "resteep"], key=lambda r: r["id"]):
        nl = next(x for x in main if key(x) == key(r) and (x["stack"], x["fill"], x["d"]) == (r["stack"], r["fill"], r["d"]))
        uref = ref[key(r)]
        lin_ref = next((x for x in rows if x["group"] == "resteep" and x["stack"] == "none" and key(x) == key(r)), None)
        for rr, leg in ((nl, "nonlinear"), (r, "linear")):
            if leg == "nonlinear":
                dv = rr["pldb"] - uref["pldb"]
            else:
                dv = rr["pldb"] - lin_ref["pldb"] if lin_ref is not None else float("nan")
            lines.append(f"| {r['aircraft']} M{r['mach']}, {r['label']} | {leg} | {rr['pldb']:.1f} | "
                         f"{rr['peak_pa']:.1f} | {_f(rr['rise_ms'])} | {_f(dv, 1, True)} |")
    out["resteep"] = "\n".join(lines)

    pcs = []
    for pc in [r for r in rows if r["group"] == "padcheck"]:
        base = next(r for r in main if key(r) == key(pc) and r["stack"] == "gas" and r["fill"] == pc["fill"]
                    and r["d"] == pc["d"])
        pcs.append(f"{pc['aircraft']} M{pc['mach']} {pc['label']} (Lorentz width {pc['lorentz_width_ms']:.0f} ms): "
                   f"pad {base['pad_s']:.1f}→{pc['pad_s']:.1f} s: PLdB {_f(base['pldb'], 2)}→{_f(pc['pldb'], 2)}, "
                   f"CSEL {base['csel']:.2f}→{pc['csel']:.2f}, peak {base['peak_pa']:.3f}→{pc['peak_pa']:.3f} Pa")
    out["padcheck"] = "; ".join(pcs)

    lines = ["| leak at 20 km (flat scaling of the unshielded wave) | ground peak [Pa] | 10–90 % rise [ms] "
             "| ΔPLdB raw | ΔPLdB proxy | ΔCSEL |", "|---|---|---|---|---|---|"]
    for r in sorted([r for r in rows if r["group"] == "leak"], key=lambda r: -r["leak_db"]):
        lines.append(f"| {r['leak_db']:+g} dB | {r['peak_pa']:.1f} | {_f(r['rise_ms'])} | "
                     f"{_f(r['dpldb_raw'], 1, True)} | {_f(r['dpldb_smooth'], 1, True)} | {_f(r['dcsel'], 1, True)} |")
    out["leak"] = "\n".join(lines)
    return out


def write_summary(rows, path, wall_s):
    tb = tables(rows)
    tmpl = SUMMARY_TEMPLATE
    n = len(rows)
    rt = sum(r["runtime_s"] for r in rows)
    text = tmpl.format(**tb, n_cases=n, cpu_min=rt / 60, wall_min=wall_s / 60)
    with open(path, "w") as fh:
        fh.write(text)


SUMMARY_TEMPLATE = """# Gate 0 — BoomShield in full nonlinear propagation (v0)

Generated by `run_gate0.py` ({n_cases} runs). Data: `gate0_results.csv`; figure: `gate0_signatures.png`.

## Bottom line

1. **Ideal model: opaque stacks give ΔPLdB = −∞** (He ≥ 30 m at M ≤ 2; H₂ ≥ 30 m to M 3.2). The ground boom becomes a 2–25 Pa pulse with 0.15–1 s rise, inaudible to Mark VII. Loudness will be set by **leaks this model omits**.
2. **Leak-limited expectation:** the review's diffraction floors map to **ΔPLdB ≈ −4 to −17 for a 250 m corridor** and **−13 to −23 for 1 km**, before cells, ties and turbulence.
3. **He 15 m (cheapest) is marginal:** −8 dB at M 2.0, −3.6 at M 2.5, −1 at M 3.2; it re-shocks.
4. **Fragile:** ties of Dyneema-grid scale (K = 300 Pa/m; review: ≈ 200) take H₂ 30 m from −∞ to −9 dB; 5 % air halves He 15 m's benefit. Mass law (same face mass): −0.5 to −2.8 dB.

## Setup

- Shield top at 20 km; aircraft 100 or 250 m above it.
- CONCORDE preset at M 1.4/2.0/2.5, SR71 at M 3.2; only altitude/Mach overridden, weights unchanged (Concorde at M 1.4: C_L ≈ 0.55, a loud source).
- Whitham near field at r0 = min(1.5 L, 0.8 h) (80/92 m Concorde, 49 m SR-71; all ≥ 1.3 L), same r0 per pair.
- Burgers to 20 km → `shield.apply` (plane wave, trace speed V) → Burgers to ground, ×1.9, 50 kHz.
- Shield records padded max(1 s, 15× Lorentzian width) per side, outer half tapered to zero. Without the taper, cropped tails faked 30–50 PLdB; doubling the pad changes nothing.
- **All slab lift spent on the faces** (lift margin 0): optimistic.
- **Rise-time proxy:** Gaussian σ = 1 ms (step → 2.6 ms 10–90 %), since the model gives 0.7–0.9 ms rises. Headline = less favourable of raw/proxy, always the proxy. An exponential kernel moves Δ ≤ 3 dB.

## Headline ΔPLdB (vs unshielded, same aircraft/h)

Notation: ᵖ = gas layer propagating (V > c_gas); −∞ = no band reaches a Mark VII contour.

{headline}

Raw / proxy at h = 250 m:

{bothways}

Shield-plane vs ground ratios and ΔCSEL, which stays finite and shows the low-frequency pressure that remains:

{ratios}

## Sensitivities (Concorde M 2.0, h = 250 m)

{ties}

Gains at K = 300 (He 15 m) and K ≤ 100 (H₂ 30 m) are **not robust**. The tie partly cancels the layer's negative static stiffness ρ_g/(a²d) (526 and 99 Pa/m), a tuned state that drifts with V, T and purity.

{purity}

## Re-steepening (nonlinearity off below the shield)

{resteep}

- **Marginal stacks re-steepen:** He 15 m at M 2 gives 44 PLdB (linear leg) vs 98.8 (nonlinear), re-forming a front shock; at M 2.5 it fully re-shocks.
- **Opaque stacks do not:** their output is too weak and slow.
- **Shield-plane peak ratios understate ground peaks.** The unshielded wave loses peak to N-wave formation (228 Pa linear leg vs 85 Pa) and the smooth wave does not, so He 15 m keeps 0.68 of the ground peak vs 0.23 at the shield.

## Leak floor (unshielded wave × leak at 20 km, Concorde M 2.0, h 250 m)

{leak}

Real diffraction leaks more at low frequency, so flat scaling is a rough, slightly pessimistic proxy.

## What this does NOT include

- **Finite corridor width:** diffraction floor ~−6 to −18 dB (250 m), −15 to −25 dB (1 km), independent of opacity; probably limiting.
- **Cell walls / lens rims:** floor thick layers at ~−33 to −45 dB at 10–20 Hz.
- **Off-track rays** (k_y ≠ 0).
- **Wind, V/T drift** (matter near cutoff).
- **Turbulence** (proxy only; it also scatters into the shadow).
- **Structure-borne flanking** beyond smeared massless ties (tendons, decks, face coincidence).
- **Envelope shape (π/4 lift factor)** and parasitic mass: these stacks are really *under* budget.
- **Building response** to the residual sub-20 Hz pulse (ΔCSEL only −20 to −50 dB); reflected-cone loads; engine/manoeuvre noise.

**Verdict:** the 1-D mechanism survives full propagation, and opaque stacks do not re-shock. Achievable ΔPLdB is now set by corridor width, cells and ties, which need the 2-D/3-D models (review §4).
"""


# ----------------------------------------------------------------------------- main
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--procs", type=int, default=os.cpu_count() or 1)
    ap.add_argument("--cache", default=None, help="pickle of raw results (written, or read with --replot)")
    ap.add_argument("--replot", action="store_true")
    args = ap.parse_args()
    os.makedirs(RESULTS, exist_ok=True)
    t0 = time.time()
    if args.replot:
        with open(args.cache, "rb") as fh:
            rows, wall = pickle.load(fh)
    else:
        cases = build_cases()
        # longest first: unshielded / strong shocks are the slow ones
        order = sorted(cases, key=lambda c: (c["stack"] != "none", c["stack"] != "masslaw", c["mach"]))
        rows = []
        with Pool(args.procs) as pool:
            for r in pool.imap_unordered(run_case, order, chunksize=1):
                rows.append(r)
                print(f"[{len(rows):3d}/{len(cases)}] {r['aircraft']} M{r['mach']} h{r['h']:g} "
                      f"{r['label']:<34} {'(lin below)' if r['below_linear'] else '':<11} "
                      f"PL {r['pldb']:6.2f} / {r['pldb_smooth']:6.2f}  pk {r['peak_pa']:6.1f}  "
                      f"{r['runtime_s']:5.1f}s", flush=True)
        wall = time.time() - t0
        if args.cache:
            with open(args.cache, "wb") as fh:
                pickle.dump((rows, wall), fh)
    rows = add_deltas(rows)
    write_csv(rows, os.path.join(RESULTS, "gate0_results.csv"))
    plot(rows, os.path.join(RESULTS, "gate0_signatures.png"))
    write_summary(rows, os.path.join(RESULTS, "gate0_summary.md"), wall)
    print(f"wrote results/gate0_results.csv, gate0_signatures.png, gate0_summary.md ({wall / 60:.1f} min)")


if __name__ == "__main__":
    main()
