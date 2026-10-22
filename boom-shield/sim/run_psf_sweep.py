"""Peak-overpressure (psf) sweep: BoomShield vs the proposed 0.11 psf overland limit.

Question: which floating-lid stacks (gas, thickness d, shield altitude H_s) bring the
*ground peak overpressure* undertrack to <= 0.11 psf (5.27 Pa) for Concorde-class,
SR-71 and Overture-like aircraft flying h = 250 m above the shield; and what limits it.

Pipeline per case (Gate-0 1-D undertrack chain, run_gate0.py conventions):
Whitham near field at r0 = min(1.5 L, 0.8 h) -> augmented Burgers to H_s -> zero pad
(+ half-cosine taper of the outer half of the pad) -> shield.apply (plane wave, trace
speed V, ambient std_atm(H_s)) -> augmented Burgers to the ground -> x1.9 -> metrics
(raw and Gate-0 rise-time proxy).  The above-shield leg (aircraft -> H_s) is shared by
every stack of one aircraft/Mach/H_s and cached.  Below the shield, smooth transmitted
waves are propagated at FS_LO (decimated); a case is re-run at the full 50 kHz whenever
the reduced-rate result could be audible or re-steepen (see run_case).

Stacks (top to bottom), lift = (rho_air - rho_gas) d at H_s:
  float : Membrane(lid = lift - 0.1) / gas d / Membrane(0.1)    lift-closed floating lid
  fixed : Membrane(2.4)              / gas d / Membrane(0.1)    margin = lift - 2.5
  split : Membrane(lift/2) / gas d / Membrane(lift/2)            Gate-0 / 3-D convention
                                                                  (validation only)
Outputs: results/psf_sweep.csv, results/psf_sweep.png (results/psf_sweep.md is written
by hand from the printed tables).  Cache: results/psf_cache/ (resumable).

Run:  python boom-shield/sim/run_psf_sweep.py --phase validate|coarse|refine|extra|3d|report
      [--procs 2]
"""
from __future__ import annotations

import argparse
import csv
import dataclasses
import json
import math
import os
import pickle
import sys
import time
from multiprocessing import get_context

import numpy as np
from scipy.signal import resample_poly

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from boomsim import propagate as PR  # noqa: E402
from boomsim import shield as S  # noqa: E402
from boomsim.atmosphere import std_atm  # noqa: E402
from boomsim.proxy import scored_metrics  # noqa: E402
from boomsim.source import CONCORDE, SR71, Aircraft  # noqa: E402

RESULTS = os.path.join(HERE, "results")
CACHE = os.path.join(RESULTS, "psf_cache")

PSF = 47.880258            # Pa per psf
LIMIT_PSF = 0.11
LIMIT_PA = LIMIT_PSF * PSF  # 5.267 Pa
FS = 50e3
FS_LO = 5e3                # below-shield rate for smooth (opaque) transmissions
H_AC = 250.0               # aircraft height above the shield top [m]
FLOOR_M = 0.1              # floor film [kg/m^2]
FIXED_LID = 2.4            # fixed-lid variant [kg/m^2]
TAPER_FRAC = 0.5

# Boom Overture-like preset, copied from boom-shield/infographic/make_data.py:
#  length 201 ft = 61.3 m; cruise M 1.7; MTOW ~77 t; mid-cruise ~0.81 MTOW -> 62 t ->
#  W = 0.61 MN; volume ~300 m^3; lift on the wing, ~0.40 L ... 0.90 L.
OVERTURE = Aircraft("Overture", length=61.3, mach=1.7, altitude=18300.0,
                    weight=62000.0 * 9.81, volume=300.0,
                    lift_start=0.40, lift_end=0.90)
PRESETS = {"Concorde": CONCORDE, "SR-71": SR71, "Overture": OVERTURE}
# reference wing areas for the C_L flag [m^2]: Concorde 358.25 (type data), SR-71 167
# (1,800 ft^2), Overture ~240 (ASSUMED: no published figure; Concorde-like wing loading
# at its 77 t MTOW).  Operating limits used for flags: C_L > 0.5, q > 35 kPa
# (~ Concorde/SR-71 equivalent-airspeed limits, 450-530 KEAS).
WING_AREA = {"Concorde": 358.25, "SR-71": 167.0, "Overture": 240.0}
CL_MAX, Q_MAX = 0.5, 35e3

AIRCRAFT = [("Concorde", 1.4), ("Overture", 1.7), ("Concorde", 2.0), ("Concorde", 2.5),
            ("SR-71", 3.2)]
H_LIST = [14e3, 16e3, 18e3, 20e3, 22e3, 25e3]
D_LIST = [15.0, 30.0, 45.0, 60.0, 80.0, 100.0, 130.0, 160.0, 200.0]
D_COARSE = [15.0, 30.0, 60.0, 100.0, 160.0, 200.0]
GASES = ["helium", "hydrogen"]
LIDS = ["float", "fixed"]


# ----------------------------------------------------------------------------- setup
def aircraft(name, mach, H_s):
    return dataclasses.replace(PRESETS[name], mach=float(mach), altitude=float(H_s) + H_AC)


def flight_flags(name, mach, H_s):
    ac = aircraft(name, mach, H_s)
    a = ac.atm
    q = 0.5 * 1.4 * float(a.p) * mach ** 2
    CL = ac.weight / (q * WING_AREA[name])
    return dict(q_kpa=q / 1e3, CL=CL, flag=("CL" if CL > CL_MAX else "") + ("q" if q > Q_MAX else ""))


def make_stack(fill, d, lid, H_s, tie_K=0.0):
    atm = std_atm(H_s)
    g = S.gas(fill, float(atm.T), float(atm.p))
    lift = (float(atm.rho) - g.rho) * d
    if lid == "float":
        top, bot = max(lift - FLOOR_M, 1e-3), FLOOR_M
    elif lid == "fixed":
        top, bot = FIXED_LID, FLOOR_M
    elif lid == "split":
        top = bot = lift / 2
    else:
        raise ValueError(lid)
    st = S.Stack([S.Membrane(top), S.GasLayer(g, d, tie_K=tie_K), S.Membrane(bot)])
    return st, g, lift, top, bot, lift - top - bot


def lf_tau(fill, d, lid, mach, H_s):
    """Low-frequency time constant of |T|: |T|^-2 = 1 + (omega tau)^2 + O(omega^4).

    From the small-kappa*d expansion of the stack matrix (limp faces m1, m2, evanescent
    gas, a^2 = 1/V^2 - 1/c_g^2, Z_a = rho_a/sqrt(1/c^2 - 1/V^2)):
      den/2 = 1 + omega^2 a^2 d (d + m/rho_g)/2 + j omega [(m + rho_g d)/Z_a - Z_a a^2 d/rho_g]/2
      tau^2 = a^2 d (d + m/rho_g) + [Z_a a^2 d/rho_g - (m + rho_g d)/Z_a]^2 / 4,  m = m1 + m2.
    """
    atm = std_atm(H_s)
    st, g, lift, top, bot, _ = make_stack(fill, d, lid, H_s)
    V = mach * float(atm.c)
    a2 = 1 / V ** 2 - 1 / g.c ** 2
    Za = float(atm.rho) / math.sqrt(1 / float(atm.c) ** 2 - 1 / V ** 2)
    m = top + bot
    return math.sqrt(a2 * d * (d + m / g.rho) + (Za * a2 * d / g.rho - (m + g.rho * d) / Za) ** 2 / 4)


def case_key(c):
    if c["stack"] == "none":
        return f"{c['ac']}_M{c['mach']:g}_H{c['H_s'] / 1e3:g}_none"
    k = f"{c['ac']}_M{c['mach']:g}_H{c['H_s'] / 1e3:g}_{c['fill'][:2]}{c['d']:g}_{c['lid']}"
    if c.get("tie_K"):
        k += f"_K{c['tie_K']:g}"
    if c.get("variant"):
        k += "_" + c["variant"]
    return k


def mk(ac, mach, H_s, stack="gas", fill=None, d=0.0, lid="float", **kw):
    c = dict(ac=ac, mach=float(mach), H_s=float(H_s), stack=stack, fill=fill, d=float(d), lid=lid,
             tie_K=0.0, variant="", force_hi=False)
    c.update(kw)
    c["key"] = case_key(c)
    return c


def evanescent(fill, mach, H_s):
    atm = std_atm(H_s)
    g = S.gas(fill, float(atm.T), float(atm.p))
    return mach * float(atm.c) < g.c


# ----------------------------------------------------------------------------- cache
def _p(key, what="res"):
    return os.path.join(CACHE, f"{key}.{what}.pkl")


def load(key, what="res"):
    p = _p(key, what)
    if os.path.exists(p):
        with open(p, "rb") as fh:
            return pickle.load(fh)
    return None


def save(key, obj, what="res"):
    os.makedirs(CACHE, exist_ok=True)
    p = _p(key, what)
    with open(p + ".tmp", "wb") as fh:
        pickle.dump(obj, fh)
    os.replace(p + ".tmp", p)


# ----------------------------------------------------------------------------- legs
def above_leg(name, mach, H_s):
    """Near field -> H_s (nonlinear, 50 kHz). Cached per aircraft/Mach/H_s."""
    key = f"above_{name}_M{mach:g}_H{H_s / 1e3:g}"
    got = load(key, "above")
    if got is not None:
        return got
    ac = aircraft(name, mach, H_s)
    r0 = min(1.5 * ac.length, 0.8 * H_AC)
    t, p = PR.near_field(ac, r0, FS)
    ray = PR.RayPath(ac.altitude, ac.V, z_min=0.0)
    z0 = ac.altitude - r0
    t, p = PR.propagate(t, p, ray, z0, H_s)
    out = dict(t0=float(t[0]), dt=float(t[1] - t[0]), p=p.astype(np.float64), r0=r0)
    save(key, out, "above")
    return out


def _pad_taper(t, p, pad, out_fn):
    dt = t[1] - t[0]
    n = int(round(pad / dt))
    p2 = np.concatenate([np.zeros(n), p, np.zeros(n)])
    t2 = t[0] + (np.arange(p2.size) - n) * dt
    out = out_fn(t2, p2)
    nt = max(2, int(TAPER_FRAC * n))
    w = np.ones(out.size)
    ramp = 0.5 * (1 - np.cos(np.pi * np.arange(nt) / nt))
    w[:nt] = ramp
    w[-nt:] = ramp[::-1]
    return t2, out * w


def _pos_impulse(t, p):
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
    """One sweep case -> result dict (metrics + small stored traces)."""
    got = load(c["key"])
    if got is not None:
        return got
    t_start = time.time()
    ac = aircraft(c["ac"], c["mach"], c["H_s"])
    V = ac.V
    atmS = std_atm(c["H_s"])
    ab = above_leg(c["ac"], c["mach"], c["H_s"])
    t = ab["t0"] + np.arange(ab["p"].size) * ab["dt"]
    p = ab["p"].copy()
    res = dict(c)
    res.update(V=V, weight_N=ac.weight, **flight_flags(c["ac"], c["mach"], c["H_s"]))
    if c["stack"] == "none":
        st, pad, tau = None, 1.0, 0.0
        res.update(lift=float("nan"), lid_m=0.0, floor_m=0.0, margin=float("nan"), c_gas=float("nan"),
                   regime="-", tau_lf=0.0)
    else:
        st, g, lift, top, bot, margin = make_stack(c["fill"], c["d"], c["lid"], c["H_s"], c["tie_K"])
        a = math.sqrt(max(0.0, 1 / V ** 2 - 1 / g.c ** 2))
        tau = lf_tau(c["fill"], c["d"], c["lid"], c["mach"], c["H_s"]) if V < g.c else 0.0
        # pad: Gate-0 rule (15 Lorentz widths, >= 1 s) extended to 8 low-frequency time
        # constants (thick layers), capped at 12 s
        pad = min(30.0 if c["variant"] == "pad30" else 12.0, max(1.0, 15.0 * a * c["d"], 8.0 * tau))
        res.update(lift=lift, lid_m=top, floor_m=bot, margin=margin, c_gas=g.c,
                   regime="evanescent" if V < g.c else "propagating", tau_lf=tau)
    fn = (lambda tt, pp: S.apply(st, tt, pp, V, atmS)) if st is not None else (lambda tt, pp: pp.copy())
    tb, pb = _pad_taper(t, p, pad, fn)
    res["shield_peak_above"] = float(p.max())
    res["shield_peak_below"] = float(pb.max())
    ray = PR.RayPath(ac.altitude, V, z_min=0.0)

    def below(fs_dec):
        if fs_dec == FS:
            tt, pp = tb, pb
        else:
            q = int(round(FS / fs_dec))
            pp = resample_poly(pb, 1, q)
            tt = tb[0] + np.arange(pp.size) * (q * (tb[1] - tb[0]))
        tg, pg = PR.propagate(tt, pp, ray, c["H_s"], 0.0)
        return tg, 1.9 * pg

    # reduced rate first for shielded cases; full rate if anything could be audible/steep
    used = FS
    if st is not None and not c["force_hi"]:
        tg, pg = below(FS_LO)
        m = scored_metrics(tg, pg)
        hf = (np.isfinite(m["pldb"]) and m["pldb"] > 20.0) or not (m["rise_time_10_90_ms"] > 20.0)
        if hf:
            tg, pg = below(FS)
            m = scored_metrics(tg, pg)
        else:
            used = FS_LO
    else:
        tg, pg = below(FS)
        m = scored_metrics(tg, pg)
    res.update(fs_below=used, pad_s=pad, pldb=m["pldb"], pldb_proxy=m["pldb_proxy"],
               peak_pa=m["peak_pa"], peak_psf=m["peak_pa"] / PSF, rise_ms=m["rise_time_10_90_ms"],
               duration_ms=m["duration_ms"], pos_impulse=m["pos_impulse_pa_s"], csel=m["csel"],
               asel=m["asel"], neg_peak_pa=float(pg.min()), runtime_s=time.time() - t_start)
    # stored traces (decimated, peak-centred window)
    ipk = int(np.argmax(pg))
    dtg = tg[1] - tg[0]
    w = int(round(3.0 / dtg))
    s0, s1 = max(0, ipk - w), min(pg.size, ipk + w)
    dec = max(1, int(round(1e-3 / dtg)))      # 1 kHz storage
    res["trace"] = (tg[s0:s1:dec].astype(np.float32), pg[s0:s1:dec].astype(np.float32))
    save(c["key"], res)
    return res


def run_all(cases, procs, tag=""):
    todo = [c for c in cases if load(c["key"]) is None]
    print(f"{tag}: {len(cases)} cases, {len(todo)} to run", flush=True)
    # above legs first (serially per aircraft group, in the pool)
    legs = sorted({(c["ac"], c["mach"], c["H_s"]) for c in todo})
    t0 = time.time()
    rows = []
    with get_context("fork").Pool(procs) as pool:
        for _ in pool.starmap(above_leg, legs, chunksize=1):
            pass
        order = sorted(todo, key=lambda c: (c["stack"] != "none", c["d"]))
        for i, r in enumerate(pool.imap_unordered(run_case, order, chunksize=1)):
            print(f"[{i + 1:4d}/{len(todo)}] {r['key']:<38} pk {r['peak_pa']:7.2f} Pa "
                  f"({r['peak_psf']:.3f} psf) PL {r['pldb']:6.1f}/{r['pldb_proxy']:6.1f} "
                  f"rise {r['rise_ms']:7.1f} ms fs {r['fs_below'] / 1e3:g}k {r['runtime_s']:5.1f}s "
                  f"[{time.time() - t0:6.0f}s]", flush=True)
    return [load(c["key"]) for c in cases]


# ----------------------------------------------------------------------------- case lists
def baseline_cases():
    return [mk(a, m, H, "none") for a, m in AIRCRAFT for H in H_LIST]


def grid_cases(ds, lids=LIDS, hs=H_LIST):
    C = []
    for a, m in AIRCRAFT:
        for H in hs:
            for fill in GASES:
                if not evanescent(fill, m, H):
                    continue          # He at M >= ~2.9 propagates: skip (Gate 0: no benefit)
                for d in ds:
                    for lid in lids:
                        C.append(mk(a, m, H, "gas", fill, d, lid))
    return C


def validation_cases():
    """Gate-0 / 3-D comparison (split faces, H_s 20 km) + reduced-rate checks."""
    C = []
    for a, m, fill, d in (("Concorde", 2.0, "hydrogen", 30.0), ("Concorde", 2.0, "helium", 30.0),
                          ("Concorde", 2.0, "helium", 60.0), ("SR-71", 3.2, "hydrogen", 30.0),
                          ("Concorde", 1.4, "helium", 30.0), ("Concorde", 2.0, "hydrogen", 60.0)):
        C.append(mk(a, m, 20e3, "gas", fill, d, "split", force_hi=True, variant="hi"))
        C.append(mk(a, m, 20e3, "gas", fill, d, "split"))
        C.append(mk(a, m, 20e3, "gas", fill, d, "float", force_hi=True, variant="hi"))
    return C


# ----------------------------------------------------------------------------- report
CSV_COLS = ["key", "ac", "mach", "V", "weight_N", "H_s", "q_kpa", "CL", "flag", "stack", "fill", "d",
            "lid", "tie_K", "variant", "regime", "c_gas", "lift", "lid_m", "floor_m", "margin", "tau_lf",
            "pad_s", "fs_below", "peak_pa", "peak_psf", "pass_0p11", "neg_peak_pa", "pldb", "pldb_proxy",
            "dpldb_headline", "csel", "dcsel", "pos_impulse", "rise_ms", "duration_ms",
            "shield_peak_above", "shield_peak_below", "base_peak_pa", "base_pldb", "runtime_s"]


def collect():
    rows = []
    for f in sorted(os.listdir(CACHE)) if os.path.isdir(CACHE) else []:
        if f.endswith(".res.pkl"):
            with open(os.path.join(CACHE, f), "rb") as fh:
                rows.append(pickle.load(fh))
    base = {(r["ac"], r["mach"], r["H_s"]): r for r in rows if r["stack"] == "none"}
    for r in rows:
        b = base.get((r["ac"], r["mach"], r["H_s"]))
        r["pass_0p11"] = bool(r["peak_pa"] <= LIMIT_PA)
        if b is None:
            continue
        r["base_peak_pa"], r["base_pldb"] = b["peak_pa"], b["pldb"]
        r["dpldb_headline"] = max(r["pldb"] - b["pldb"], r["pldb_proxy"] - b["pldb_proxy"])
        r["dcsel"] = r["csel"] - b["csel"]
    return rows


def write_csv(rows, path):
    with open(path, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=CSV_COLS, extrasaction="ignore")
        w.writeheader()
        for r in sorted(rows, key=lambda r: (r["ac"], r["mach"], r["H_s"], r["stack"] != "none",
                                             r.get("fill") or "", r["lid"], r["d"], r["variant"])):
            w.writerow({k: (f"{r.get(k):.6g}" if isinstance(r.get(k), float) else r.get(k, ""))
                        for k in CSV_COLS})


def min_thickness(rows, ac, mach, H, fill, lid):
    """Smallest d meeting the limit (log-interpolated between grid points); None if none."""
    rr = sorted([r for r in rows if (r["ac"], r["mach"], r["H_s"], r.get("fill"), r["lid"]) ==
                 (ac, mach, H, fill, lid) and r["stack"] == "gas" and not r["variant"] and not r["tie_K"]],
                key=lambda r: r["d"])
    if not rr:
        return None, None
    prev = None
    for r in rr:
        if r["peak_pa"] <= LIMIT_PA:
            if prev is None:
                return r["d"], f"≤{r['d']:g}"
            # log-log interpolation of peak(d)
            x0, x1 = math.log(prev["d"]), math.log(r["d"])
            y0, y1 = math.log(prev["peak_pa"]), math.log(r["peak_pa"])
            x = x0 + (math.log(LIMIT_PA) - y0) * (x1 - x0) / (y1 - y0)
            return math.exp(x), f"{math.exp(x):.0f}"
        prev = r
    return None, f">{rr[-1]['d']:g}"


def report(rows):
    print("\n=== baselines (unshielded): peak Pa / psf, PLdB raw / proxy, q, C_L ===")
    for a, m in AIRCRAFT:
        line = []
        for H in H_LIST:
            b = next((r for r in rows if (r["ac"], r["mach"], r["H_s"], r["stack"]) == (a, m, H, "none")), None)
            if b:
                line.append(f"H{H / 1e3:g}: {b['peak_pa']:.1f}/{b['peak_psf']:.2f} {b['pldb']:.1f}/{b['pldb_proxy']:.1f}"
                            f" q{b['q_kpa']:.0f} CL{b['CL']:.2f}{'!' + b['flag'] if b['flag'] else ''}")
        print(f"{a} M{m}: " + " | ".join(line))
    for lid in LIDS:
        print(f"\n=== min thickness [m] for <= 0.11 psf, lid={lid} (rows aircraft x gas, cols H_s) ===")
        print("| aircraft | gas | " + " | ".join(f"{H / 1e3:g} km" for H in H_LIST) + " |")
        for a, m in AIRCRAFT:
            for fill in GASES:
                cells = []
                for H in H_LIST:
                    if not evanescent(fill, m, H):
                        cells.append("n/a (V>c_gas)")
                        continue
                    _, s = min_thickness(rows, a, m, H, fill, lid)
                    cells.append(s or "–")
                print(f"| {a} M{m} | {'He' if fill == 'helium' else 'H₂'} | " + " | ".join(cells) + " |")
    print("\n=== peak [Pa] vs d at H_s 20 km, float lid ===")
    for a, m in AIRCRAFT:
        for fill in GASES:
            rr = sorted([r for r in rows if (r["ac"], r["mach"], r["H_s"], r.get("fill"), r["lid"]) ==
                         (a, m, 20e3, fill, "float") and not r["variant"] and not r["tie_K"]], key=lambda r: r["d"])
            if rr:
                print(f"{a} M{m} {fill[:2]}: " + "  ".join(f"{r['d']:g}:{r['peak_pa']:.2f}/{r['pldb_proxy']:.0f}"
                                                         for r in rr))


# ----------------------------------------------------------------------------- figure
def plot(rows, path):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    ink, ink2, muted, grid, surf = "#0b0b0b", "#52514e", "#8a8984", "#e4e3df", "#fcfcfb"
    pal = {"Concorde M1.4": "#2a78d6", "Overture M1.7": "#1baf7a", "Concorde M2": "#eb6834",
           "Concorde M2.5": "#8f4fd6", "SR-71 M3.2": "#c43c5a"}
    plt.rcParams.update({"font.size": 8.5, "axes.edgecolor": muted, "axes.labelcolor": ink,
                         "xtick.color": ink2, "ytick.color": ink2, "axes.titlecolor": ink})
    fig, axes = plt.subplots(2, 4, figsize=(20, 9.4))
    fig.patch.set_facecolor(surf)

    def series(a, m, H, fill, lid):
        rr = sorted([r for r in rows if (r["ac"], r["mach"], r["H_s"], r.get("fill"), r["lid"]) ==
                     (a, m, H, fill, lid) and r["stack"] == "gas" and not r["variant"] and not r["tie_K"]],
                    key=lambda r: r["d"])
        return rr

    # (0,0) H2 / (0,1) He: peak vs d at 20 km, float lid; fixed lid dashed
    for j, fill in enumerate(("hydrogen", "helium")):
        ax = axes[0, j]
        for a, m in AIRCRAFT:
            lab = f"{a} M{m:g}"
            for lid, ls in (("float", "-"), ("fixed", "--")):
                rr = series(a, m, 20e3, fill, lid)
                if rr:
                    ax.plot([r["d"] for r in rr], [r["peak_pa"] for r in rr], ls, color=pal[lab], lw=1.6,
                            marker="o" if lid == "float" else None, ms=3,
                            label=lab if lid == "float" else None)
        ax.axhline(LIMIT_PA, color=ink, lw=1.0, ls=":")
        ax.text(16, LIMIT_PA * 1.07, "0.11 psf (5.27 Pa)", color=ink, fontsize=8)
        ax.set_xscale("log")
        ax.set_yscale("log")
        ax.set_xlabel("gas thickness d [m]")
        ax.set_ylabel("ground peak overpressure, undertrack [Pa]")
        ax.set_title(f"{'H₂' if fill == 'hydrogen' else 'He'} layer, shield at 20 km "
                     "(solid: floating lid, dashed: fixed 2.4 kg/m² lid)", fontsize=9)
        ax.legend(frameon=False, fontsize=7.5)
    # (0,2) altitude dependence: H2 float, Concorde M2 & Overture, several d
    ax = axes[0, 2]
    for (a, m), mk_ in ((("Concorde", 2.0), "o"), (("Overture", 1.7), "s")):
        for d, col in ((30.0, "#eb6834"), (60.0, "#2a78d6"), (100.0, "#1baf7a"), (200.0, "#8f4fd6")):
            pts = []
            for H in H_LIST:
                rr = [r for r in series(a, m, H, "hydrogen", "float") if r["d"] == d]
                if rr:
                    pts.append((H / 1e3, rr[0]["peak_pa"]))
            if pts:
                ax.plot(*zip(*pts), "-" if a == "Concorde" else "--", marker=mk_, ms=3, color=col, lw=1.4,
                        label=f"{a} M{m:g}, H₂ {d:g} m")
    ax.axhline(LIMIT_PA, color=ink, lw=1.0, ls=":")
    ax.set_yscale("log")
    ax.set_xlabel("shield altitude H_s [km] (aircraft at H_s + 250 m)")
    ax.set_ylabel("ground peak [Pa]")
    ax.set_title("Altitude dependence (H₂, floating lid)", fontsize=9)
    ax.legend(frameon=False, fontsize=7)
    # (1,0) PLdB proxy vs peak: all shielded cases
    ax = axes[1, 0]
    for a, m in AIRCRAFT:
        lab = f"{a} M{m:g}"
        rr = [r for r in rows if (r["ac"], r["mach"]) == (a, m) and r["stack"] == "gas" and not r["variant"]
              and not r["tie_K"]]
        x = np.array([r["peak_pa"] for r in rr])
        y = np.array([r["pldb_proxy"] if np.isfinite(r["pldb_proxy"]) else 5.0 for r in rr])
        ax.scatter(x, y, s=7, color=pal[lab], alpha=0.7, label=lab)
    ax.axvline(LIMIT_PA, color=ink, lw=1.0, ls=":")
    ax.axhline(75, color=muted, lw=0.8, ls="--")
    ax.text(LIMIT_PA * 1.05, 92, "0.11 psf", fontsize=8)
    ax.text(1.1, 77, "~75 PLdB (X-59 target)", fontsize=7.5, color=ink2)
    ax.text(1.1, 8, "PLdB −∞ (plotted at 5)", fontsize=7.5, color=ink2)
    ax.set_xscale("log")
    ax.set_xlabel("ground peak [Pa]")
    ax.set_ylabel("PLdB, rise-time proxy (1-D, no edge leaks)")
    ax.set_title("Loudness vs peak: all shielded 1-D cases", fontsize=9)
    ax.legend(frameon=False, fontsize=7)
    # (1,1) floor model: peak vs analytic lift-swell estimate
    ax = axes[1, 1]
    fl = load("floor_model", "aux")
    for a, m in AIRCRAFT:
        lab = f"{a} M{m:g}"
        rr = [r for r in rows if (r["ac"], r["mach"]) == (a, m) and r["stack"] == "gas" and not r["variant"]
              and not r["tie_K"] and r["regime"] == "evanescent" and r.get("swell_est")]
        ax.scatter([r["swell_est"] for r in rr], [r["peak_pa"] for r in rr], s=7, color=pal[lab], alpha=0.7,
                   label=lab)
    lim = np.array([0.3, 100])
    ax.plot(lim, lim, color=ink, lw=0.8)
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel("lift-swell estimate  C·1.9·R(H)·βW / (2π √(2βh) (Vτ)^{3/2})  [Pa]")
    ax.set_ylabel("simulated ground peak [Pa]")
    ttl = "Floor physics: peak follows the smoothed-lift estimate"
    if fl:
        ttl += f" (C = {fl['C']:.2f})"
    ax.set_title(ttl, fontsize=9)
    ax.legend(frameon=False, fontsize=7)
    # (1,2) signatures: Concorde M2 20 km H2 float
    ax = axes[1, 2]
    b = next((r for r in rows if (r["ac"], r["mach"], r["H_s"], r["stack"]) == ("Concorde", 2.0, 20e3, "none")), None)
    if b:
        t, p = b["trace"]
        i = int(np.argmax(p))
        ax.plot(t - t[i], p, color=ink2, lw=1.2, label=f"unshielded ({b['peak_pa']:.0f} Pa)")
    cols = ["#eb6834", "#2a78d6", "#1baf7a", "#8f4fd6"]
    for d, col in zip((30.0, 60.0, 100.0, 200.0), cols):
        rr = [r for r in series("Concorde", 2.0, 20e3, "hydrogen", "float") if r["d"] == d]
        if rr:
            t, p = rr[0]["trace"]
            i = int(np.argmax(p))
            ax.plot(t - t[i], p, color=col, lw=1.6, label=f"H₂ {d:g} m ({rr[0]['peak_pa']:.1f} Pa)")
    ax.axhline(LIMIT_PA, color=ink, lw=1.0, ls=":")
    ax.set_xlim(-1.5, 2.0)
    ax.set_ylim(-15, 30)
    ax.set_xlabel("time re peak [s]")
    ax.set_ylabel("ground overpressure [Pa]")
    ax.set_title("Concorde M2, shield 20 km: ground swell (unshielded clipped)", fontsize=9)
    ax.legend(frameon=False, fontsize=7.5)
    # (0,3) 3-D carpets: peak vs ground offset (Concorde M2, h 250 m, shield 20 km)
    ax = axes[0, 3]
    pts3 = {}
    with open(os.path.join(RESULTS, "3d_results.csv")) as fh:
        for r in csv.DictReader(fh):
            pts3.setdefault(r["key"], []).append((float(r["y_g_m"]) / 1e3, float(r["peak_pa"])))
    for k, lab, col, ls in (("C2_h250_none", "unshielded (÷4)", ink2, ":"),
                            ("C2_h250_W1000_H230", "W 1 km, H₂ 30 m (lift-closed split)", "#eb6834", "-"),
                            ("C2_h250_W1000_He60", "W 1 km, He 60 m (split)", "#8f4fd6", "-")):
        if k in pts3:
            y, pk = zip(*sorted(pts3[k]))
            pk = np.array(pk) / (4.0 if "none" in k else 1.0)
            ax.plot(y, pk, ls, color=col, lw=1.4, label=lab + " [3d_results]")
    for W, col in ((1000.0, "#2a78d6"), (2000.0, "#1baf7a")):
        o = load(f"3d_Concorde_M2_hy60_float_W{W:g}_R3000", "aux")
        if o:
            rr = [r for r in o["rows"] if "peak_pa" in r]
            ax.plot([r["y_g"] / 1e3 for r in rr], [r["peak_pa"] for r in rr], "-", marker="o", ms=3.5,
                    color=col, lw=1.8, label=f"W {W / 1e3:g} km, H₂ 60 m floating lid [new 3-D]")
    offs = []
    for phi in (0.0, 24.0, 39.0, 54.0):
        o = load(f"off_Concorde_M2_H20_hy60_float_phi{phi:g}", "aux")
        if o:
            offs.append((phi, o["peak_pa"]))
    if offs:
        ys = {0.0: 0.0, 24.0: 9.6, 39.0: 18.3, 54.0: 39.3}
        ax.plot([ys[p] for p, _ in offs], [v for _, v in offs], "x", color="#2a78d6", ms=7, mew=1.6,
                label="H₂ 60 m, INFINITE shield, off-track 1-D")
    ax.axhline(LIMIT_PA, color=ink, lw=1.0, ls=":")
    ax.set_xlabel("ground offset from track y_g [km]")
    ax.set_ylabel("ground peak [Pa]")
    ax.set_title("3-D check: finite width puts the peak in the outer carpet", fontsize=9)
    ax.legend(frameon=False, fontsize=6.8)
    ax.set_ylim(0, 26)
    # (1,3) 1-D vs 3-D: undertrack peak and carpet-maximum peak
    ax = axes[1, 3]
    for x, yu, ym, lab in compare_1d_3d():
        ax.plot([x, x], [yu, ym], color=muted, lw=0.8)
        ax.plot(x, yu, "o", color="#2a78d6", ms=5)
        ax.plot(x, ym, "^", color="#eb6834", ms=6)
        ax.annotate(lab, (x, ym), textcoords="offset points", xytext=(4, 3), fontsize=6.5, color=ink2)
    ax.plot([], [], "o", color="#2a78d6", label="3-D undertrack (R_m 8 km where run)")
    ax.plot([], [], "^", color="#eb6834", label="3-D maximum over the computed carpet")
    lim = np.array([2, 40])
    ax.plot(lim, lim, color=ink, lw=0.8)
    ax.axhline(LIMIT_PA, color=ink, lw=1.0, ls=":")
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel("1-D undertrack ground peak [Pa]")
    ax.set_ylabel("3-D ground peak [Pa]")
    ax.set_title("1-D screen vs 3-D (W = 1 km unless noted)", fontsize=9)
    ax.legend(frameon=False, fontsize=7, loc="upper left")
    for ax in axes.flat:
        if not ax.axison:
            continue
        ax.set_facecolor(surf)
        ax.grid(True, color=grid, lw=0.6, which="both")
        for s in ("top", "right"):
            ax.spines[s].set_visible(False)
    fig.suptitle("BoomShield vs the 0.11 psf overland limit — 1-D undertrack chain, floating-lid stacks, "
                 "aircraft 250 m above the shield", color=ink, fontsize=11)
    fig.tight_layout(rect=(0, 0, 1, 0.97))
    fig.savefig(path, dpi=120, facecolor=surf)
    plt.close(fig)


def compare_1d_3d():
    """[(1-D peak, 3-D undertrack peak, 3-D carpet max, label)] for the 3-D cases that exist."""
    pts3 = {}
    with open(os.path.join(RESULTS, "3d_results.csv")) as fh:
        for r in csv.DictReader(fh):
            pts3.setdefault(r["key"], []).append((float(r["phi_deg"]), float(r["peak_pa"])))
    out = []
    for k3, k1, lab in (("C2_h250_W1000_H230", "Concorde_M2_H20_hy30_split", "C M2 H₂30"),
                        ("C2_h250_W1000_He30", "Concorde_M2_H20_he30_split", "C M2 He30"),
                        ("C2_h250_W1000_He60", "Concorde_M2_H20_he60_split", "C M2 He60"),
                        ("S3.2_h250_W1000_H230", "SR-71_M3.2_H20_hy30_split", "SR-71 H₂30"),
                        ("C1.4_h250_W1000_He30", "Concorde_M1.4_H20_he30_split", "C M1.4 He30")):
        r1 = load(k1)
        if r1 is None or k3 not in pts3:
            continue
        far = dict(pts3.get(k3 + "_R8000", []))
        yu = far.get(0.0, dict(pts3[k3])[0.0])
        out.append((r1["peak_pa"], yu, max(v for _, v in pts3[k3]), lab))
    for W in (1000.0, 2000.0):
        o = load(f"3d_Concorde_M2_hy60_float_W{W:g}_R3000", "aux")
        r1 = load("Concorde_M2_H20_hy60_float")
        if o and r1:
            rr = [r for r in o["rows"] if "peak_pa" in r]
            out.append((r1["peak_pa"], rr[0]["peak_pa"], max(r["peak_pa"] for r in rr),
                        f"C M2 H₂60 float W{W / 1e3:g}km"))
    return out


# ----------------------------------------------------------------------------- floor model
def ray_gain(name, mach, H_s):
    """Linear undertrack amplitude ratio p(ground)/p(H_s) for a ray from the aircraft
    (no absorption), times the 1.9 reflection."""
    ac = aircraft(name, mach, H_s)
    ray = PR.RayPath(ac.altitude, ac.V, z_min=0.0)
    return 1.9 * float(ray.amplitude_ratio(H_s, 0.0))


def swell_estimate(r):
    """Lift-swell estimate (uncalibrated, C = 1): the transmitted k_x -> 0 part is the
    aircraft's lift spread over l = V tau along track; its undertrack wave at the shield
    plane is that of a lifting line of weight W and length l:
        p_s ~ beta W / (2 pi sqrt(2 beta h) l^{3/2})  (Whitham, F ~ A_L / l^{3/2}),
    then x 1.9 R(H) linear to the ground."""
    ac = aircraft(r["ac"], r["mach"], r["H_s"])
    b = ac.beta
    ell = ac.V * r["tau_lf"]
    if ell <= 0:
        return float("nan")
    return ray_gain(r["ac"], r["mach"], r["H_s"]) * b * ac.weight / (2 * math.pi * math.sqrt(2 * b * H_AC) * ell ** 1.5)


def floor_analysis(rows):
    out = {}
    sh = [r for r in rows if r["stack"] == "gas" and not r["variant"] and not r["tie_K"]
          and r["regime"] == "evanescent"]
    for r in sh:
        r["swell_est"] = swell_estimate(r)
    # calibrate C on thick-ish layers (l >= 2 L), log-mean
    fit = [r for r in sh if r["V"] * r["tau_lf"] >= 2 * aircraft(r["ac"], r["mach"], r["H_s"]).length
           and r["fs_below"] == FS_LO]
    lr = np.array([math.log(r["peak_pa"] / r["swell_est"]) for r in fit])
    C = float(np.exp(lr.mean())) if lr.size else float("nan")
    out.update(C=C, n_fit=len(fit), spread_db=float(20 * lr.std() / math.log(10)) if lr.size else float("nan"),
               ratio_minmax=(float(np.exp(lr.min()) / C), float(np.exp(lr.max()) / C)) if lr.size else None)
    for r in sh:
        r["swell_est"] *= C
    # ray gains and quasi-static footprint
    gains = {f"{a}_M{m:g}_H{H / 1e3:g}": ray_gain(a, m, H) for a, m in AIRCRAFT for H in H_LIST}
    out["gains"] = gains
    qs = {}
    for a, m in AIRCRAFT:
        ac = aircraft(a, m, 20e3)
        # quasi-static lift footprint on the ground: integral = W over ~ pi beta H^2
        qs[f"{a}_M{m:g}"] = ac.weight / (math.pi * ac.beta * (20e3 + H_AC) ** 2)
    out["quasi_static_pa"] = qs
    save("floor_model", out, "aux")
    return out


def lift_scaling_cases():
    """Same stack, aircraft weight x0.5 and volume ~0: is the swell linear in W and
    lift-dominated?"""
    C = []
    for d in (30.0, 100.0, 200.0):
        for v in ("halfW", "novol"):
            C.append(mk("Concorde", 2.0, 20e3, "gas", "hydrogen", d, "float", variant=v))
    return C


def run_variant(c):
    """run_case with a modified preset (variant halfW / novol / ymax)."""
    name = c["ac"]
    orig = PRESETS[name]
    if c["variant"] == "halfW":
        PRESETS[name] = dataclasses.replace(orig, weight=0.5 * orig.weight)
    elif c["variant"] == "novol":
        PRESETS[name] = dataclasses.replace(orig, volume=1e-3)
    try:
        # separate above-leg cache for the modified source
        global above_leg
        f0 = above_leg

        def ab(n, m, H):
            key = f"above_{n}_M{m:g}_H{H / 1e3:g}_{c['variant']}"
            got = load(key, "above")
            if got is not None:
                return got
            ac = aircraft(n, m, H)
            r0 = min(1.5 * ac.length, 0.8 * H_AC)
            t, p = PR.near_field(ac, r0, FS)
            ray = PR.RayPath(ac.altitude, ac.V, z_min=0.0)
            t, p = PR.propagate(t, p, ray, ac.altitude - r0, H)
            out = dict(t0=float(t[0]), dt=float(t[1] - t[0]), p=p, r0=r0)
            save(key, out, "above")
            return out
        above_leg = ab
        return run_case(c)
    finally:
        above_leg = f0
        PRESETS[name] = orig


def ymax_case(d):
    """Near-field record length check: F computed to 12 L instead of 3 L."""
    c = mk("Concorde", 2.0, 20e3, "gas", "hydrogen", d, "float", variant="ymax12")
    got = load(c["key"])
    if got is not None:
        return got
    nf0 = PR.near_field
    PR.near_field = lambda ac, r0, fs, **kw: nf0(ac, r0, fs, y_max_factor=12.0, **kw)
    global above_leg
    f0 = above_leg

    def ab(n, m, H):
        key = f"above_{n}_M{m:g}_H{H / 1e3:g}_ymax12"
        got = load(key, "above")
        if got is not None:
            return got
        ac = aircraft(n, m, H)
        r0 = min(1.5 * ac.length, 0.8 * H_AC)
        t, p = PR.near_field(ac, r0, FS)
        ray = PR.RayPath(ac.altitude, ac.V, z_min=0.0)
        t, p = PR.propagate(t, p, ray, ac.altitude - r0, H)
        out = dict(t0=float(t[0]), dt=float(t[1] - t[0]), p=p, r0=r0)
        save(key, out, "above")
        return out
    above_leg = ab
    try:
        return run_case(c)
    finally:
        above_leg = f0
        PR.near_field = nf0


# ----------------------------------------------------------------------------- 3-D check
def run_3d_confirm(fill, d, lid="float", mach=2.0, acn="Concorde", W=1000.0, R_m=3000.0,
                   phis=(0.0, 12.0, 24.0, 33.0, 39.0, 45.0, 51.0, 54.0), workers=2):
    """Finite-width slab through the run_3d.py pipeline (cross-plane BEM -> matching
    circle -> propagate_from_circle), floating-lid face masses, shield top at 20 km.
    Cached in results/psf_cache."""
    key = f"3d_{acn}_M{mach:g}_{fill[:2]}{d:g}_{lid}_W{W:g}_R{R_m:g}"
    got = load(key, "aux")
    if got is not None:
        return got
    import run_3d as R3
    from boomsim import crossplane as cp
    R3.WORKERS = workers
    ac = aircraft(acn, mach, R3.H_S)
    if W is None:
        cfg = cp.CrossPlaneConfig(ac, W=None, h=H_AC, f_max=R3.F_MAX, R_m=R_m, phis=list(phis),
                                  incident="aged", workers=workers, df=R3.DF_MAX)
        top = bot = tau = 0.0
    else:
        st, g, lift, top, bot, _ = make_stack(fill, d, lid, R3.H_S)
        tau = lf_tau(fill, d, lid, mach, R3.H_S)
        cfg = cp.CrossPlaneConfig(ac, W=W, T=d, profile="slab", gas=g, m_top=top, m_bot=bot,
                                  h=H_AC, f_max=R3.F_MAX, R_m=R_m, phis=list(phis),
                                  incident="aged", workers=workers)
        # record >= 6 s and >= 2 s + 15 low-frequency time constants (thick layers)
        cfg.df = min(R3.DF_MAX, cfg.auto_df(), 1.0 / math.ceil(2.0 + 15.0 * tau))
    t0 = time.time()
    f, P = cp.matching_spectra(cfg)
    wall_sp = time.time() - t0
    print(f"  [{key}] spectra {f.size} freqs df={cfg.df:.4f} Hz: {wall_sp:.0f} s", flush=True)
    R3._G.update(cfg=cfg, f=f, P=P)
    rows = []
    import multiprocessing as mp
    with mp.get_context("fork").Pool(workers) as pool:
        for j, row in pool.imap_unordered(R3._ground_job, [(j, "full", R3.CROP_REL) for j in range(len(phis))]):
            rows.append(row)
    rows.sort(key=lambda r: r["phi"])
    out = dict(key=key, fill=fill, d=d, lid=lid, W=W, R_m=R_m, top=top, bot=bot, wall_spec=wall_sp,
               df=cfg.df, tau=tau,
               rows=[dict(phi=r["phi"], y_g=r["y_g"], circle_peak=r["circle_peak"], sig=r.get("sig"),
                          **(r["metrics"] or {})) for r in rows], wall=time.time() - t0)
    save(key, out, "aux")
    for r in out["rows"]:
        print(f"  phi {r['phi']:4.0f}  y_g {r['y_g'] / 1e3:5.1f} km  peak {r.get('peak_pa', float('nan')):6.2f} Pa  "
              f"PLdB {r.get('pldb', float('nan')):6.1f} / {r.get('pldb_proxy', float('nan')):6.1f}", flush=True)
    print(f"  total {out['wall']:.0f} s", flush=True)
    return out


def offtrack_1d(fill, d, lid, phi_deg, mach=2.0, acn="Concorde", H_s=20e3):
    """1-D chain on the off-track ray at azimuth phi under an INFINITE shield.

    With limp (N = 0) faces the layered transfer depends only on the horizontal
    slowness magnitude, so the oblique incidence is shield.apply with the effective
    trace speed V_eff = 1/sqrt(1/V^2 + s_y^2), s_y = cos(mu0) sin(phi)/c_a (Snell
    invariant of the ray; c = c_a through the isothermal layer)."""
    key = f"off_{acn}_M{mach:g}_H{H_s / 1e3:g}_{(fill or 'none')[:2]}{d:g}_{lid}_phi{phi_deg:g}"
    got = load(key, "aux")
    if got is not None:
        return got
    ac = aircraft(acn, mach, H_s)
    atmS = std_atm(H_s)
    phi = math.radians(phi_deg)
    c_a = float(ac.atm.c)
    sy = math.sqrt(1 - (c_a / ac.V) ** 2) * math.sin(phi) / c_a
    V_eff = 1.0 / math.sqrt(1 / ac.V ** 2 + sy ** 2)
    if fill is None:
        fn = None
        pad = 1.0
    else:
        st, g, *_ = make_stack(fill, d, lid, H_s)
        a = math.sqrt(max(0.0, 1 / V_eff ** 2 - 1 / g.c ** 2))
        pad = min(12.0, max(1.0, 15.0 * a * d, 8.0 * lf_tau(fill, d, lid, mach, H_s)))

        def fn(t, p):
            return _pad_taper(t, p, pad, lambda tt, pp: S.apply(st, tt, pp, V_eff, atmS))
    r0 = min(1.5 * ac.length, 0.8 * H_AC)
    if fn is None:
        t, p = PR.ground_signature(ac, r0=r0, fs=FS, phi=phi)
    else:
        t, p = PR.ground_signature(ac, r0=r0, insert=(H_s, fn), fs=FS, phi=phi)
    m = scored_metrics(t, p)
    out = dict(key=key, phi=phi_deg, V_eff=V_eff, peak_pa=m["peak_pa"], pldb=m["pldb"],
               pldb_proxy=m["pldb_proxy"], rise_ms=m["rise_time_10_90_ms"])
    save(key, out, "aux")
    print(f"  {key}: V_eff {V_eff:.0f} m/s peak {m['peak_pa']:.2f} Pa PLdB {m['pldb']:.1f}", flush=True)
    return out


# ----------------------------------------------------------------------------- main
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--phase", required=True,
                    choices=["validate", "baseline", "coarse", "refine", "extra", "3d", "offtrack", "report",
                             "time"])
    ap.add_argument("--W", type=float, default=1000.0)
    ap.add_argument("--R", type=float, default=3000.0)
    ap.add_argument("--procs", type=int, default=2)
    ap.add_argument("--fill", default="hydrogen")
    ap.add_argument("--d", type=float, default=60.0)
    ap.add_argument("--mach", type=float, default=2.0)
    ap.add_argument("--ac", default="Concorde")
    args = ap.parse_args()
    t0 = time.time()
    if args.phase == "time":
        c = mk(args.ac, args.mach, 20e3, "gas", args.fill, args.d, "float", variant="timing")
        for k in (c["key"],):
            p = _p(k)
            if os.path.exists(p):
                os.remove(p)
        r = run_case(c)
        print({k: r[k] for k in ("peak_pa", "pldb", "pldb_proxy", "rise_ms", "fs_below", "runtime_s", "pad_s")})
    elif args.phase == "validate":
        run_all(validation_cases(), args.procs, "validate")
    elif args.phase == "baseline":
        run_all(baseline_cases(), args.procs, "baseline")
    elif args.phase == "coarse":
        run_all(baseline_cases() + grid_cases(D_COARSE), args.procs, "coarse")
    elif args.phase == "refine":
        run_all(grid_cases(D_LIST), args.procs, "refine")
    elif args.phase == "extra":
        for c in lift_scaling_cases():
            r = run_variant(c)
            print(c["key"], f"{r['peak_pa']:.3f}")
        for d in (30.0, 100.0, 200.0):
            r = ymax_case(d)
            print(r["key"], f"{r['peak_pa']:.3f}")
        tie = [mk("Concorde", 2.0, 20e3, "gas", "hydrogen", d, "float", tie_K=10.0) for d in (30.0, 60.0, 100.0, 200.0)]
        tie += [mk(a, m, 20e3, "gas", "hydrogen", 200.0, "float", variant="pad30")
                for a, m in (("Concorde", 1.4), ("Concorde", 2.0))]
        run_all(tie, args.procs, "ties K=10")
    elif args.phase == "3d":
        W = None if args.W <= 0 else args.W
        run_3d_confirm(args.fill, args.d, "float", args.mach, args.ac, W=W, R_m=args.R, workers=args.procs)
    elif args.phase == "offtrack":
        import itertools
        jobs = list(itertools.product([(None, 0.0), ("hydrogen", 30.0), ("hydrogen", 60.0), ("helium", 60.0)],
                                      [0.0, 24.0, 39.0, 54.0]))
        with get_context("fork").Pool(args.procs) as pool:
            pool.starmap(offtrack_1d, [(f, d, "float", ph) for (f, d), ph in jobs])
    elif args.phase == "report":
        rows = collect()
        fl = floor_analysis(rows)
        print(json.dumps({k: v for k, v in fl.items() if k != "gains"}, indent=1))
        write_csv(rows, os.path.join(RESULTS, "psf_sweep.csv"))
        plot(rows, os.path.join(RESULTS, "psf_sweep.png"))
        report(rows)
    print(f"done in {time.time() - t0:.0f} s")


if __name__ == "__main__":
    main()
