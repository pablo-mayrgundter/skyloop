"""3-D finite-width corridor study: ground carpets under a cross-track-finite BoomShield.

Question: with finite corridor width W (edge diffraction), edge treatments, cell walls,
ties and off-track rays, how much does ground PLdB drop across the whole boom carpet?

Pipeline per configuration (contract: DESIGN-3D.md):
  1. cross-plane solver (boomsim.crossplane): total-field spectra P(f, phi) on a matching
     circle of radius R_m around the flight axis (uniform air at the shield altitude,
     aged incident), for downward azimuths phi >= 0 (all configurations are symmetric);
  2. per phi: band-limited signature synthesised at 50 kHz, cropped to where
     |p| > CROP_REL max|p| (+0.1 s margins, 0.05 s half-Hann ends; same rule for every
     run), then propagate.propagate_from_circle(ac, phi, t, p, R_m, z_ref=20 km,
     age_from=aging_radius(cfg, phi)) -> ground signature and lateral offset y_g;
  3. proxy.scored_metrics (PLdB raw and with the Gate-0 rise-time proxy, CSEL, ...).
Unshielded baselines go through exactly the same path (cfg with W=None).

Extras: R_m sensitivity (W = 1 km He 30 m at R_m = 2, 4, 6, 8 km; 12 and 18 km indicative,
the circle then reaches below the 11 km tropopause), near-track re-runs (phi = 0, 4, 12, 24 deg)
of every shielded config at R_m = 8 km (near-track levels are not converged at 3 km),
line-vs-smeared tie check on the circle, band split (< 5 Hz vs
> 5 Hz parts of the circle signature propagated separately) at phi = 0 and at the
loudest azimuth of each shielded config, direct near-field cross-check of the unshielded
baseline (propagate.ground_carpet start="near_field"), crop-threshold check.

Outputs: results/3d_results.csv (per phi), results/3d_carpet_summary.csv,
results/3d_carpets.png, results/3d_signatures.png; tables printed for 3d_summary.md.
Cache (resumable, per config and per phi): results/3d_cache/.

Run:  python boom-shield/sim/run_3d.py [--only KEY ...] [--report] [--procs 4]
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

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from boomsim import crossplane as cp  # noqa: E402
from boomsim import propagate as PR  # noqa: E402
from boomsim import shield as SH  # noqa: E402
from boomsim.atmosphere import std_atm  # noqa: E402
from boomsim.proxy import scored_metrics  # noqa: E402
from boomsim.source import CONCORDE, SR71  # noqa: E402

RESULTS = os.path.join(HERE, "results")
CACHE = os.path.join(RESULTS, "3d_cache")

H_S = 20000.0             # shield top [m]
FS = 50e3                 # ground-propagation sample rate (as Gate 0)
F_MAX = 100.0             # cross-plane band limit [Hz]
DF_MAX = 1.0 / 6.0        # Hz: record >= 6 s (auto_df of the widest shield); baselines use it
R_MAIN = 3000.0           # matching radius [m] (= 3 W at W = 1 km, 12 W at 250 m)
R_SENS = (2000.0, 4000.0, 6000.0, 8000.0)  # R_m sensitivity (2 W, 4 W, 6 W, 8 W at W = 1 km)
CROP_REL = 1e-3           # crop threshold of the circle signature
CROP_MARGIN, CROP_TAPER = 0.1, 0.05
F_BAND = 5.0              # band split [Hz] (raised-cosine crossover 4..6 Hz)
PHI_GRID = [0, 1, 2, 3, 4, 5, 6, 8, 10, 12, 14, 16, 18, 20, 22, 24, 26, 28, 30,
            33, 36, 39, 42, 45, 48, 51, 54, 57]
PHI_SENS = [0, 4, 8, 12, 16, 20, 24, 28, 33, 39, 45, 51, 54]
PHI_SENS_FAR = [0, 4, 12, 24]   # R_m = 6, 8 km (convergence of the near-track values)
R_BEYOND = (12000.0, 18000.0)   # indicative: circle extends below the 11 km tropopause
R_FAR = 8000.0            # far circle for the near-track re-runs (circle start stays above 11 km)
PLOT_PHIS = (0, 10, 20, 30, 39, 45, 51)
SIG_WIN = (-0.2, 1.6)     # s, stored ground-signature window
SIG_DEC = 5               # 50 kHz -> 10 kHz for stored traces

ATM20 = std_atm(H_S)
T20, P20, RHO20 = float(ATM20.T), float(ATM20.p), float(ATM20.rho)


# ----------------------------------------------------------------------------- setup
def aircraft(name, mach, h):
    base = {"Concorde": CONCORDE, "SR-71": SR71}[name]
    return dataclasses.replace(base, mach=float(mach), altitude=H_S + float(h))


def phis_for(ac, grid=PHI_GRID):
    phc = math.degrees(PR.lateral_cutoff(ac.altitude, ac.V)[0])
    return [float(p) for p in grid if p < phc - 0.5]


def gas(fill):
    return SH.gas(fill, T20, P20)


def build_configs():
    """List of dicts. 'base' = key of the unshielded baseline of the same aircraft/h."""
    C = []

    def add(key, group, acn, mach, h, W=None, fill=None, T=0.0, profile="slab", fringe=0.0,
            cells=0.0, tieK=0.0, R=R_MAIN, phis=None, label=None):
        C.append(dict(key=key, group=group, acn=acn, mach=mach, h=h, W=W, fill=fill, T=T,
                      profile=profile, fringe=fringe, cells=cells, tieK=tieK, R=R,
                      phis=phis, label=label or key))

    def base(acn, mach, h):
        k = f"{'C' if acn == 'Concorde' else 'S'}{mach:g}_h{h:g}"
        add(k + "_none", "baseline", acn, mach, h, label=f"{acn} M{mach:g} h{h:g}: unshielded")
        return k

    k = base("Concorde", 2.0, 250)
    for W in (1000.0, 250.0, 500.0):
        for fill, g in (("helium", "He"), ("hydrogen", "H2")):
            add(f"{k}_W{W:g}_{g}30", "width", "Concorde", 2.0, 250, W, fill, 30.0,
                label=f"W {W:g} m, {g} 30 m slab")
    add(f"{k}_W1000_He60", "stack", "Concorde", 2.0, 250, 1000.0, "helium", 60.0,
        label="W 1000 m, He 60 m slab")
    add(f"{k}_W1000_He60lens", "stack", "Concorde", 2.0, 250, 1000.0, "helium", 60.0, "lens",
        label="W 1000 m, He 60 m lens (rim 6 m)")
    for fill, g in (("helium", "He"), ("hydrogen", "H2")):
        add(f"{k}_W1000_{g}30_fr150", "edge", "Concorde", 2.0, 250, 1000.0, fill, 30.0,
            fringe=150.0, label=f"W 1000 m, {g} 30 m, raised-cosine fringe 150 m")
    add(f"{k}_W1000_He30_cells100", "cells", "Concorde", 2.0, 250, 1000.0, "helium", 30.0,
        cells=100.0, label="W 1000 m, He 30 m, partitions every 100 m")
    for K in (30.0, 300.0):
        add(f"{k}_W1000_H230_ties{K:g}", "ties", "Concorde", 2.0, 250, 1000.0, "hydrogen", 30.0,
            tieK=K, label=f"W 1000 m, H2 30 m, line ties 25 m (K = {K:g} Pa/m)")
    for R in R_SENS:
        add(f"{k}_W1000_He30_R{R:g}", "Rsens", "Concorde", 2.0, 250, 1000.0, "helium", 30.0,
            R=R, phis=PHI_SENS if R <= 4000.0 else PHI_SENS_FAR,
            label=f"W 1000 m, He 30 m, R_m = {R:g} m")
    add(f"{k}_none_R8000", "Rsens", "Concorde", 2.0, 250, R=8000.0, phis=[0, 12, 24],
        label="unshielded, R_m = 8000 m (pipeline R_m check)")
    # beyond the isothermal layer (uniform-air circle below 11 km): indicative only
    for R in R_BEYOND:
        add(f"{k}_W1000_He30_R{R:g}", "Rsens", "Concorde", 2.0, 250, 1000.0, "helium", 30.0,
            R=R, phis=PHI_SENS_FAR[:3], label=f"W 1000 m, He 30 m, R_m = {R:g} m (indicative)")
    add(f"{k}_none_R18000", "Rsens", "Concorde", 2.0, 250, R=18000.0, phis=[0, 12],
        label="unshielded, R_m = 18000 m (pipeline R_m check)")
    # other aircraft / h
    k = base("SR-71", 3.2, 250)
    add(f"{k}_W1000_H230", "aircraft", "SR-71", 3.2, 250, 1000.0, "hydrogen", 30.0,
        label="SR-71 M3.2: W 1000 m, H2 30 m")
    k = base("Concorde", 1.4, 250)
    add(f"{k}_W1000_He30", "aircraft", "Concorde", 1.4, 250, 1000.0, "helium", 30.0,
        label="Concorde M1.4: W 1000 m, He 30 m")
    k = base("Concorde", 2.0, 100)
    add(f"{k}_W1000_He30", "aircraft", "Concorde", 2.0, 100, 1000.0, "helium", 30.0,
        label="Concorde M2 h100: W 1000 m, He 30 m")
    # far matching circle (R_m = 8 km, near-track azimuths): the R_m check shows the
    # near-track ground levels are not converged at 3 km (see 3d_summary.md)
    have = {c["key"] for c in C}
    for c in list(C):
        if c["group"] in ("baseline", "Rsens") or c["key"] + "_R8000" in have:
            continue
        C.append(dict(c, key=c["key"] + "_R8000", group="far", R=R_FAR, phis=PHI_SENS_FAR,
                      label=c["label"] + " [R_m 8 km]"))
    for c in C:
        if c["group"] != "baseline":
            c["base"] = "_".join(c["key"].split("_")[:2]) + "_none"
        else:
            c["base"] = c["key"]
    return C


def lens_mean_depth(T, T_rim, q=0.5):
    """Footprint-mean gas depth of the lens profile T_rim + (T - T_rim)(1 - u^2)^q."""
    from scipy.special import beta as B
    return T_rim + (T - T_rim) * 0.5 * B(0.5, q + 1.0)


def make_cfg(c, phis=None):
    ac = aircraft(c["acn"], c["mach"], c["h"])
    phis = c["phis"] if phis is None else phis
    phis = phis_for(ac) if phis is None else [float(p) for p in phis]
    kw = dict(h=float(c["h"]), f_max=F_MAX, R_m=c["R"], phis=phis, incident="aged", workers=4)
    if c["W"] is None:
        cfg = cp.CrossPlaneConfig(ac, W=None, df=DF_MAX, **kw)
        return cfg
    g = gas(c["fill"])
    T = c["T"]
    if c["profile"] == "lens":
        d_lift = lens_mean_depth(T, 0.1 * T)
    else:
        d_lift = T
    face = (RHO20 - g.rho) * d_lift / 2.0       # all slab lift spent on the two faces
    W = c["W"]
    extra = {}
    if c["cells"]:
        extra["partitions"] = list(np.arange(-W / 2 + c["cells"], W / 2 - 1e-6, c["cells"]))
    if c["tieK"]:
        s = 25.0
        extra["line_ties"] = [(float(y), c["tieK"] * s) for y in np.arange(-W / 2 + s / 2, W / 2, s)]
    cfg = cp.CrossPlaneConfig(ac, W=W, T=T, profile=c["profile"], gas=g, m_top=face, m_bot=face,
                              fringe_width=c["fringe"], **extra, **kw)
    cfg.df = min(DF_MAX, cfg.auto_df())
    return cfg


# ----------------------------------------------------------------------------- caching
def _path(key, what):
    return os.path.join(CACHE, f"{key}_{what}.pkl")


def _load(key, what):
    p = _path(key, what)
    if os.path.exists(p):
        with open(p, "rb") as fh:
            return pickle.load(fh)
    return None


def _save(key, what, obj):
    os.makedirs(CACHE, exist_ok=True)
    p = _path(key, what)
    with open(p + ".tmp", "wb") as fh:
        pickle.dump(obj, fh)
    os.replace(p + ".tmp", p)


# ----------------------------------------------------------------------------- circle
def taper_weights(cfg, f):
    w = np.ones(f.size)
    f0 = cfg.taper_from * cfg.f_max
    m = f > f0
    w[m] = 0.5 * (1 + np.cos(np.pi * (f[m] - f0) / (cfg.f_max - f0)))
    return w


def band_weight(f, band):
    lo = np.clip((F_BAND + 1.0 - f) / 2.0, 0.0, 1.0)
    lo = np.sin(0.5 * np.pi * lo) ** 2           # 1 below 4 Hz, 0 above 6 Hz, 0.5 at 5 Hz
    if band == "lo":
        return lo
    if band == "hi":
        return 1.0 - lo
    return np.ones_like(f)


def lf_fraction(cfg, f, P):
    """Fraction of the circle signature energy (0 < f <= f_max, band-limit taper applied)
    below F_BAND, per circle point."""
    E = np.abs(P) ** 2 * taper_weights(cfg, f)[:, None] ** 2
    E[0] = 0.0
    return E[f < F_BAND].sum(0) / E.sum(0)


def circle(c):
    """Matching spectra (cached)."""
    sp = _load(c["key"], "spec")
    if sp is not None:
        return sp
    cfg = make_cfg(c)
    t0 = time.time()
    f, P = cp.matching_spectra(cfg)
    sp = dict(f=f, P=P, phis=np.asarray(cfg.phis), R=cfg.matching_radius(), df=f[1] - f[0],
              wall=time.time() - t0, f_split=cfg.split_frequency(),
              face=(None if not cfg.has_shield else float(cfg.m_top)))
    _save(c["key"], "spec", sp)
    print(f"  [{c['key']}] circle spectra: {f.size} freqs, df={f[1] - f[0]:.3f} Hz, "
          f"R_m={sp['R']:.0f} m, {sp['wall']:.0f} s", flush=True)
    return sp


# ----------------------------------------------------------------------------- ground
_G = {}


def crop(t, p, rel=CROP_REL):
    a = np.abs(p)
    idx = np.flatnonzero(a > rel * a.max())
    t0 = max(t[0], t[idx[0]] - CROP_MARGIN - CROP_TAPER)
    t1 = min(t[-1], t[idx[-1]] + CROP_MARGIN + CROP_TAPER)
    m = (t >= t0) & (t <= t1)
    tt, pp = t[m], p[m].copy()
    n = int(round(CROP_TAPER / (tt[1] - tt[0])))
    w = 0.5 * (1 - np.cos(np.pi * np.arange(n) / n))
    pp[:n] *= w
    pp[-n:] *= w[::-1]
    return tt, pp


def _ground_job(args):
    j, band, rel = args
    cfg, f, P = _G["cfg"], _G["f"], _G["P"]
    t_start = time.time()
    phi_d = float(cfg.phis[j])
    phi = math.radians(phi_d)
    cfg.fs = FS
    Pj = P[:, [j]] * band_weight(f, band)[:, None]
    t, p = cp.spectra_to_signatures(cfg, f, Pj)
    p = p[:, 0]
    circ_peak = float(p.max())
    tc, pc = crop(t, p, rel)
    tg, pg, yg = PR.propagate_from_circle(cfg.aircraft, phi, tc, pc, cfg.matching_radius(),
                                          z_ref=H_S, age_from=float(cp.aging_radius(cfg, phi)))
    row = dict(phi=phi_d, band=band, crop_rel=rel, y_g=float(yg), circle_peak=circ_peak,
               crop=(float(tc[0]), float(tc[-1])), wall=0.0)
    if tg is None:
        row["metrics"] = None
        return j, row
    m = scored_metrics(tg, pg)
    row["metrics"] = {k: float(v) for k, v in m.items()}
    if band == "full" and rel == CROP_REL and int(round(phi_d)) in PLOT_PHIS and phi_d == round(phi_d):
        s = (tg >= SIG_WIN[0]) & (tg <= SIG_WIN[1])
        row["sig"] = (tg[s][::SIG_DEC].astype(np.float32), pg[s][::SIG_DEC].astype(np.float32))
    row["wall"] = time.time() - t_start
    return j, row


def ground(c, procs, band="full", idx=None, rel=CROP_REL, what="ground"):
    """Ground results for azimuth indices idx (default all), cached per phi."""
    sp = circle(c)
    cfg = make_cfg(c)
    have = _load(c["key"], what) or {}
    idx = list(range(len(sp["phis"]))) if idx is None else list(idx)
    bands = (band,) if isinstance(band, str) else tuple(band)
    jobs = [(j, b, rel) for j in sorted(idx, key=lambda j: sp["phis"][j]) for b in bands
            if (float(sp["phis"][j]), b, rel) not in have]
    todo = jobs
    if todo:
        _G.update(cfg=cfg, f=sp["f"], P=sp["P"])
        t0 = time.time()
        import multiprocessing as mp
        with mp.get_context("fork").Pool(min(procs, len(jobs))) as pool:
            for j, row in pool.imap_unordered(_ground_job, jobs):
                have[(row["phi"], row["band"], row["crop_rel"])] = row
                _save(c["key"], what, have)
        print(f"  [{c['key']}] ground {'+'.join(bands)} rel={rel:g}: {len(todo)} jobs, "
              f"{time.time() - t0:.0f} s", flush=True)
    return have


def direct_check(c, procs):
    d = _load(c["key"], "direct")
    if d is not None:
        return d
    ac = aircraft(c["acn"], c["mach"], c["h"])
    phis = phis_for(ac)
    t0 = time.time()
    res = PR.ground_carpet(ac, np.radians(phis), start="near_field", procs=procs, keep=False)
    d = dict(phi=np.array(phis), y_g=res["y_g"], pldb=res["pldb"], pldb_proxy=res["pldb_proxy"],
             peak_pa=res["peak_pa"], csel=res["csel"], wall=time.time() - t0)
    _save(c["key"], "direct", d)
    print(f"  [{c['key']}] direct near-field carpet: {time.time() - t0:.0f} s", flush=True)
    return d


def tie_check():
    """Line ties vs smeared ties (H2 30 m, W 1 km, K = 300 Pa/m): band levels re unshielded
    on the matching circle (f_max 30 Hz, phi = 0 and 22 deg). Cached."""
    d = _load("tiecheck", "circle")
    if d is not None:
        return d
    C = build_configs()
    base = next(c for c in C if c["key"] == "C2_h250_W1000_H230")
    u = make_cfg(next(c for c in C if c["key"] == "C2_h250_none"), phis=[0.0, 22.0])
    u.f_max = 30.0
    f, P0 = cp.matching_spectra(u)
    var = {"no ties": {},
           "line 25 m": dict(line_ties=[(float(y), 300 * 25.0) for y in np.arange(-487.5, 500, 25.0)]),
           "line 12.5 m": dict(line_ties=[(float(y), 300 * 12.5) for y in np.arange(-493.75, 500, 12.5)]),
           "line 5 m": dict(line_ties=[(float(y), 300 * 5.0) for y in np.arange(-497.5, 500, 5.0)]),
           "smeared": dict(tie_K=300.0)}
    bands = [(0.1, 2), (2, 5), (5, 10), (10, 14), (14, 20), (20, 25)]
    d = dict(bands=bands, phis=[0.0, 22.0], rows={})
    for lab, kw in var.items():
        cfg = make_cfg(dict(base, tieK=0.0), phis=[0.0, 22.0])
        cfg.f_max = 30.0
        for k, v in kw.items():
            setattr(cfg, k, v)
        _, P = cp.matching_spectra(cfg)
        d["rows"][lab] = [[10 * np.log10(np.sum(abs(P[(f >= a) & (f < b), j]) ** 2)
                                         / np.sum(abs(P0[(f >= a) & (f < b), j]) ** 2)) for a, b in bands]
                          for j in range(2)]
    _save("tiecheck", "circle", d)
    return d


# ----------------------------------------------------------------------------- analysis
def carpet_rows(c):
    """Per-phi rows (band full, main crop) sorted by phi; None if not computed."""
    sp = _load(c["key"], "spec")
    g = _load(c["key"], "ground")
    if sp is None or g is None:
        return None
    cfg = make_cfg(c)
    lf = lf_fraction(cfg, sp["f"], sp["P"])
    rows = []
    for j, ph in enumerate(sp["phis"]):
        r = g.get((float(ph), "full", CROP_REL))
        if r is None:
            return None
        rows.append(dict(r, lf=float(lf[j])))
    return rows


NAN = float("nan")


def _m(r, k):
    return NAN if r["metrics"] is None else r["metrics"][k]


def deltas(rows, base):
    bmap = {r["phi"]: r for r in base}
    out = []
    for r in rows:
        b = bmap.get(r["phi"])
        if b is None:
            out.append((NAN, NAN, NAN))
            continue
        dr = _m(r, "pldb") - _m(b, "pldb")
        dp = _m(r, "pldb_proxy") - _m(b, "pldb_proxy")
        out.append((dr, dp, max(dr, dp)))
    return out


def protected_strip(y, d, thr=-10.0):
    """Contiguous half-width [m] from the track where headline delta <= thr; also the
    total protected length (one side) if the region is not contiguous."""
    y = np.asarray(y, float)
    d = np.where(np.isneginf(d), -300.0, np.asarray(d, float))
    ok = d <= thr
    if not ok[0]:
        half = 0.0
    else:
        bad = np.flatnonzero(~ok)
        if bad.size == 0:
            half = float(y[-1])
        else:
            i = bad[0]
            half = float(y[i - 1] + (y[i] - y[i - 1]) * (thr - d[i - 1]) / (d[i] - d[i - 1]))
    # total length: piecewise-linear interpolation on a fine grid
    yy = np.linspace(y[0], y[-1], 4001)
    dd = np.interp(yy, y, d)
    total = float(np.sum(dd <= thr) * (yy[1] - yy[0]))
    return half, total


def summarize(c, rows, base):
    ph = np.array([r["phi"] for r in rows])
    yg = np.array([r["y_g"] for r in rows])
    pl = np.array([_m(r, "pldb") for r in rows])
    pp = np.array([_m(r, "pldb_proxy") for r in rows])
    s = dict(key=c["key"], label=c["label"], group=c["group"], aircraft=c["acn"], mach=c["mach"],
             h=c["h"], W=c["W"] if c["W"] else NAN, R_m=c["R"], n_phi=len(rows),
             phi_max_deg=float(ph.max()), y_max_km=float(yg.max()) / 1e3,
             ut_pldb=pl[0], ut_pldb_proxy=pp[0],
             max_pldb=float(np.max(pl)), max_pldb_proxy=float(np.max(pp)),
             phi_at_max=float(ph[int(np.argmax(pp))]), y_at_max_km=float(yg[int(np.argmax(pp))]) / 1e3,
             lf_frac_phi0=rows[0]["lf"], lf_frac_min=float(min(r["lf"] for r in rows)),
             ut_peak=_m(rows[0], "peak_pa"), ut_rise=_m(rows[0], "rise_time_10_90_ms"),
             ut_csel=_m(rows[0], "csel"), ut_impulse=_m(rows[0], "pos_impulse_pa_s"))
    if base is not None and c["group"] != "baseline":
        b = {r["phi"]: r for r in base}
        d = deltas(rows, base)
        s["d_ut_raw"], s["d_ut_proxy"], s["d_ut"] = d[0]
        bpl = max(_m(b[p], "pldb") for p in ph)
        bpp = max(_m(b[p], "pldb_proxy") for p in ph)
        s["d_max_raw"] = s["max_pldb"] - bpl
        s["d_max_proxy"] = s["max_pldb_proxy"] - bpp
        s["d_max"] = max(s["d_max_raw"], s["d_max_proxy"])
        dh = np.array([x[2] for x in d])
        s["d_worst_local"] = float(np.max(dh))
        s["phi_worst_local"] = float(ph[int(np.argmax(dh))])
        s["strip_half_km"], s["strip_total_km"] = [v / 1e3 for v in protected_strip(yg, dh)]
        s["ut_csel_d"] = _m(rows[0], "csel") - _m(b[0.0], "csel")
    far = carpet_rows(dict(c, key=c["key"] + "_R8000", R=R_FAR, phis=PHI_SENS_FAR)) if c["group"] not in ("Rsens", "far", "baseline") else None
    if far and base is not None:
        fd = deltas(far, base)
        s["ut_pldb_R8"], s["ut_pldb_proxy_R8"] = _m(far[0], "pldb"), _m(far[0], "pldb_proxy")
        s["d_ut_R8"] = fd[0][2]
        s["d_near_R8"] = {r["phi"]: x[2] for r, x in zip(far, fd)}
        s["d_near_R3"] = {r["phi"]: x[2] for r, x in zip(rows, d) if r["phi"] in s["d_near_R8"]}
    return s


# ----------------------------------------------------------------------------- outputs
RES_COLS = ["key", "label", "group", "aircraft", "mach", "h_m", "W_m", "R_m", "phi_deg", "y_g_m",
            "pldb", "pldb_proxy", "csel", "peak_pa", "rise_10_90_ms", "pos_impulse_pa_s",
            "d_pldb_raw", "d_pldb_proxy", "d_pldb_headline", "lf_frac_circle_lt5Hz",
            "circle_peak_pa"]
SUM_COLS = ["key", "label", "group", "aircraft", "mach", "h", "W", "R_m", "n_phi", "phi_max_deg",
            "y_max_km", "ut_pldb", "ut_pldb_proxy", "max_pldb", "max_pldb_proxy", "phi_at_max",
            "y_at_max_km", "d_ut_raw", "d_ut_proxy", "d_ut", "d_max_raw", "d_max_proxy", "d_max",
            "d_worst_local", "phi_worst_local", "strip_half_km", "strip_total_km", "ut_peak",
            "ut_rise", "ut_impulse", "ut_csel", "ut_csel_d", "lf_frac_phi0", "lf_frac_min",
            "ut_pldb_R8", "ut_pldb_proxy_R8", "d_ut_R8"]


def _fmt(v):
    if isinstance(v, float):
        if math.isnan(v):
            return ""
        if math.isinf(v):
            return "-inf" if v < 0 else "inf"
        return f"{v:.4g}"
    return v


def write_outputs(C):
    allrows = {c["key"]: carpet_rows(c) for c in C}
    summ = []
    with open(os.path.join(RESULTS, "3d_results.csv"), "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(RES_COLS)
        for c in C:
            rows = allrows[c["key"]]
            if rows is None:
                continue
            base = allrows.get(c["base"])
            if c["group"] == "Rsens" or base is None:
                d = [(NAN, NAN, NAN)] * len(rows) if base is None else deltas(rows, base)
            else:
                d = deltas(rows, base)
            for r, dd in zip(rows, d):
                w.writerow([_fmt(x) for x in (
                    c["key"], c["label"], c["group"], c["acn"], c["mach"], c["h"],
                    c["W"] if c["W"] else NAN, c["R"], r["phi"], r["y_g"], _m(r, "pldb"),
                    _m(r, "pldb_proxy"), _m(r, "csel"), _m(r, "peak_pa"),
                    _m(r, "rise_time_10_90_ms"), _m(r, "pos_impulse_pa_s"), dd[0], dd[1], dd[2],
                    r["lf"], r["circle_peak"])])
            summ.append(summarize(c, rows, base))
    with open(os.path.join(RESULTS, "3d_carpet_summary.csv"), "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(SUM_COLS)
        for s in summ:
            w.writerow([_fmt(s.get(k, NAN)) for k in SUM_COLS])
    return allrows, summ


def _dB(x, sign=True):
    if x is None or (isinstance(x, float) and math.isnan(x)):
        return "–"
    if math.isinf(x):
        return "−∞"
    return (f"{x:+.1f}" if sign else f"{x:.1f}").replace("-", "−")


def print_tables(C, allrows, summ):
    S = {s["key"]: s for s in summ}
    print("\n### carpet summary")
    print("| config | undertrack PLdB raw/proxy | ΔPLdB undertrack (raw/proxy → headline) | max PLdB raw/proxy (φ) "
          "| ΔPLdB max-over-carpet | worst local Δ (φ) | protected strip half-width km (total) | LF frac φ=0 "
          "| ΔPLdB undertrack, R_m 8 km |")
    print("|---|---|---|---|---|---|---|---|---|")
    for c in C:
        s = S.get(c["key"])
        if s is None or c["group"] == "far":
            continue
        if c["group"] == "baseline":
            print(f"| {s['label']} | {_dB(s['ut_pldb'], 0)} / {_dB(s['ut_pldb_proxy'], 0)} | – | "
                  f"{_dB(s['max_pldb'], 0)} / {_dB(s['max_pldb_proxy'], 0)} ({s['phi_at_max']:g}°) | – | – | – | "
                  f"{s['lf_frac_phi0']:.2f} | – |")
            continue
        print(f"| {s['label']} | {_dB(s['ut_pldb'], 0)} / {_dB(s['ut_pldb_proxy'], 0)} | "
              f"{_dB(s.get('d_ut_raw'))} / {_dB(s.get('d_ut_proxy'))} → **{_dB(s.get('d_ut'))}** | "
              f"{_dB(s['max_pldb'], 0)} / {_dB(s['max_pldb_proxy'], 0)} ({s['phi_at_max']:g}°) | "
              f"**{_dB(s.get('d_max'))}** | {_dB(s.get('d_worst_local'))} ({s.get('phi_worst_local', NAN):g}°) | "
              f"{s.get('strip_half_km', NAN):.1f} ({s.get('strip_total_km', NAN):.1f}) | {s['lf_frac_phi0']:.2f} | "
              f"{_dB(s.get('d_ut_R8'))} |")
    print("\n### near-track headline ΔPLdB: R_m 3 km vs 8 km (φ = 0 / 4 / 12 / 24°)")
    print("| config | R_m 3 km | R_m 8 km |")
    print("|---|---|---|")
    for c in C:
        s = S.get(c["key"])
        if s is None or "d_near_R8" not in s:
            continue
        a, b = s["d_near_R3"], s["d_near_R8"]
        print(f"| {s['label']} | " + " / ".join(_dB(a.get(p)) for p in sorted(b)) + " | "
              + " / ".join(_dB(b[p]) for p in sorted(b)) + " |")
    # per-phi compact table for the width study
    print("\n### per-phi headline ΔPLdB (Concorde M2 h250)")
    keys = [c["key"] for c in C if c["key"].startswith("C2_h250") and c["group"] not in ("Rsens", "far")]
    rows0 = allrows.get("C2_h250_none")
    if rows0:
        print("| φ | y_g km | unsh PLdB raw/proxy | " + " | ".join(k.replace("C2_h250_", "") for k in keys[1:]) + " |")
        print("|" + "---|" * (len(keys) + 2))
        for i, r in enumerate(rows0):
            cells = []
            for k in keys[1:]:
                rr = allrows.get(k)
                if rr is None:
                    cells.append("")
                    continue
                cells.append(_dB(deltas([rr[i]], rows0)[0][2]))
            print(f"| {r['phi']:g} | {r['y_g'] / 1e3:.1f} | {_dB(_m(r, 'pldb'), 0)}/{_dB(_m(r, 'pldb_proxy'), 0)} | "
                  + " | ".join(cells) + " |")
    # ground metrics undertrack
    print("\n### undertrack ground signature metrics")
    print("| config | peak Pa | rise 10-90 ms | pos. impulse Pa s | CSEL | ΔCSEL | circle peak Pa | LF frac |")
    print("|---|---|---|---|---|---|---|---|")
    for c in C:
        s = S.get(c["key"])
        rr = allrows.get(c["key"])
        if s is None:
            continue
        print(f"| {s['label']} | {s['ut_peak']:.1f} | {s['ut_rise']:.1f} | {s['ut_impulse']:.2f} | "
              f"{s['ut_csel']:.1f} | {_dB(s.get('ut_csel_d'))} | {rr[0]['circle_peak']:.1f} | {s['lf_frac_phi0']:.2f} |")


def band_report(C, allrows):
    print("\n### band split (circle signature split at 5 Hz, each part propagated alone)")
    print("| config | φ | full PLdB raw/proxy | < 5 Hz part | > 5 Hz part |")
    print("|---|---|---|---|---|")
    out = []
    for c in C:
        b = _load(c["key"], "bands")
        rows = allrows.get(c["key"])
        if b is None or rows is None:
            continue
        for ph in sorted({k[0] for k in b}):
            full = next(r for r in rows if r["phi"] == ph)
            lo, hi = b.get((ph, "lo", CROP_REL)), b.get((ph, "hi", CROP_REL))
            if lo is None or hi is None:
                continue
            print(f"| {c['label']} | {ph:g} | {_dB(_m(full, 'pldb'), 0)}/{_dB(_m(full, 'pldb_proxy'), 0)} | "
                  f"{_dB(_m(lo, 'pldb'), 0)}/{_dB(_m(lo, 'pldb_proxy'), 0)} | "
                  f"{_dB(_m(hi, 'pldb'), 0)}/{_dB(_m(hi, 'pldb_proxy'), 0)} |")
            out.append((c["key"], ph, full, lo, hi))
    return out


def checks_report(C, allrows):
    print("\n### unshielded: circle pipeline vs direct near-field path (Concorde M2 h250)")
    c = next(c for c in C if c["key"] == "C2_h250_none")
    d = _load(c["key"], "direct")
    rows = allrows.get(c["key"])
    if d is not None and rows is not None:
        print("| φ | y_g km | circle PLdB raw/proxy | direct raw/proxy | Δ raw/proxy | peak circle/direct |")
        print("|---|---|---|---|---|---|")
        dr, dp = [], []
        for i, r in enumerate(rows):
            a, b = _m(r, "pldb"), d["pldb"][i]
            ap_, bp = _m(r, "pldb_proxy"), d["pldb_proxy"][i]
            dr.append(a - b)
            dp.append(ap_ - bp)
            if r["phi"] in (0, 2, 6, 10, 20, 30, 39, 45, 51, 54):
                print(f"| {r['phi']:g} | {r['y_g'] / 1e3:.1f} | {a:.1f}/{ap_:.1f} | {b:.1f}/{bp:.1f} | "
                      f"{a - b:+.2f}/{ap_ - bp:+.2f} | {_m(r, 'peak_pa'):.1f}/{d['peak_pa'][i]:.1f} |")
        print(f"max |Δ| raw {np.max(np.abs(dr)):.2f} dB, proxy {np.max(np.abs(dp)):.2f} dB")
    tc = _load("tiecheck", "circle")
    if tc is not None:
        print("\n### line vs smeared ties (H2 30 m, W 1 km, K = 300): circle band level re unshielded, dB")
        print("| ties | φ | " + " | ".join(f"{a:g}-{b:g} Hz" for a, b in tc["bands"]) + " |")
        print("|---|---|" + "---|" * len(tc["bands"]))
        for lab, rr in tc["rows"].items():
            for ph, r in zip(tc["phis"], rr):
                print(f"| {lab} | {ph:g} | " + " | ".join(f"{x:.1f}" for x in r) + " |")
    print("\n### crop-threshold check (1e-4 vs 1e-3)")
    for c in C:
        g = _load(c["key"], "cropcheck")
        rows = allrows.get(c["key"])
        if not g or rows is None:
            continue
        for (ph, band, rel), r in g.items():
            full = next(x for x in rows if x["phi"] == ph)
            print(f"{c['key']} φ={ph:g}: PLdB raw {_m(full, 'pldb'):.2f} -> {_m(r, 'pldb'):.2f}, proxy "
                  f"{_m(full, 'pldb_proxy'):.2f} -> {_m(r, 'pldb_proxy'):.2f}, CSEL {_m(full, 'csel'):.2f} -> "
                  f"{_m(r, 'csel'):.2f}")
    print("\n### R_m sensitivity (W 1000 He 30): ground PLdB raw/proxy")
    Rs = (2000.0, R_MAIN) + tuple(R_SENS[1:]) + tuple(R_BEYOND)
    keys = ["C2_h250_W1000_He30" if R == R_MAIN else f"C2_h250_W1000_He30_R{R:g}" for R in Rs]
    R = {k: {x["phi"]: x for x in (allrows.get(k) or [])} for k in keys}
    if R[keys[0]]:
        print("| φ | y_g km | " + " | ".join(f"R_m {r / 1e3:g} km" for r in Rs) + " |")
        print("|---|---|" + "---|" * len(Rs))
        for ph in sorted(R[keys[0]]):
            cells = []
            for k in keys:
                x = R[k].get(ph)
                cells.append("" if x is None else f"{_dB(_m(x, 'pldb'), 0)}/{_dB(_m(x, 'pldb_proxy'), 0)}")
            print(f"| {ph:g} | {R[keys[1]][ph]['y_g'] / 1e3:.1f} | " + " | ".join(cells) + " |")
    for ku in ("C2_h250_none_R8000", "C2_h250_none_R18000"):
      u8, u3 = allrows.get(ku), allrows.get("C2_h250_none")
      if u8 and u3:
        u3 = {x["phi"]: x for x in u3}
        for x in u8:
            y = u3[x["phi"]]
            print(f"unshielded φ={x['phi']:g}: {ku[13:]} {_m(x, 'pldb'):.2f}/{_m(x, 'pldb_proxy'):.2f} vs "
                  f"3 km {_m(y, 'pldb'):.2f}/{_m(y, 'pldb_proxy'):.2f} (peak {_m(x, 'peak_pa'):.1f} vs "
                  f"{_m(y, 'peak_pa'):.1f} Pa)")


# ----------------------------------------------------------------------------- figures
COL = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4", "#008300", "#6250d6", "#e34948"]
INK, INK2, GRID = "#0b0b0b", "#52514e", "#dddcd8"


def _style(ax):
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    for s in ("left", "bottom"):
        ax.spines[s].set_color(INK2)
    ax.tick_params(colors=INK2, labelsize=8)
    ax.grid(True, color=GRID, lw=0.6)
    ax.set_axisbelow(True)


def plot_carpets(C, allrows, path, floor=-5.0):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    panels = [
        ("(a) Concorde M2, h 250 m: He 30 m slab vs width", "C2_h250_none",
         ["C2_h250_W250_He30", "C2_h250_W500_He30", "C2_h250_W1000_He30"]),
        ("(b) Concorde M2, h 250 m: H$_2$ 30 m slab vs width", "C2_h250_none",
         ["C2_h250_W250_H230", "C2_h250_W500_H230", "C2_h250_W1000_H230"]),
        ("(c) W 1 km variants (Concorde M2, h 250 m)", "C2_h250_none",
         ["C2_h250_W1000_He60", "C2_h250_W1000_He60lens", "C2_h250_W1000_He30_fr150",
          "C2_h250_W1000_H230_fr150", "C2_h250_W1000_He30_cells100", "C2_h250_W1000_H230_ties30",
          "C2_h250_W1000_H230_ties300"]),
        ("(d) other aircraft / height, W 1 km (dotted: own unshielded)", None,
         [("S3.2_h250_none", "S3.2_h250_W1000_H230"), ("C1.4_h250_none", "C1.4_h250_W1000_He30"),
          ("C2_h100_none", "C2_h100_W1000_He30")]),
    ]
    lab = {c["key"]: c["label"] for c in C}
    fig, axs = plt.subplots(2, 2, figsize=(13, 12), sharey=True)
    for ax, (title, bk, keys) in zip(axs.ravel(), panels):
        _style(ax)

        def line(rows, color, label, ls="-", lw=1.8, mk="o"):
            y = np.array([r["y_g"] for r in rows]) / 1e3
            v = np.array([_m(r, "pldb_proxy") for r in rows])
            v = np.where(np.isfinite(v), v, floor)
            ax.plot(y, v, ls=ls, color=color, lw=lw, marker=mk, ms=3.5, label=label)

        if bk is not None and allrows.get(bk):
            line(allrows[bk], INK, "unshielded", lw=2.2, mk=None)
            for i, k in enumerate(keys):
                if allrows.get(k):
                    line(allrows[k], COL[i], lab[k].replace("W ", "W=").split(": ")[-1])
                far = carpet_rows(next(c for c in C if c["key"] == k + "_R8000"))
                if far:
                    y = np.array([r["y_g"] for r in far]) / 1e3
                    v = np.array([_m(r, "pldb_proxy") for r in far])
                    ax.plot(y, np.where(np.isfinite(v), v, floor), ls="none", marker="o", ms=7,
                            mfc="none", mec=COL[i], mew=1.4)
        else:
            for i, (b, k) in enumerate(keys):
                if allrows.get(b):
                    line(allrows[b], COL[i], None, ls=":", lw=1.4, mk=None)
                if allrows.get(k):
                    line(allrows[k], COL[i], lab[k])
        ax.set_title(title, fontsize=10, color=INK, loc="left")
        ax.set_xlabel("lateral ground offset y_g [km]", fontsize=9, color=INK2)
        ax.legend(fontsize=7.5, frameon=False, loc="upper center", bbox_to_anchor=(0.5, -0.12), ncol=2)
        ax.axhline(floor, color=INK2, lw=0.6, ls="--")
        ax.set_ylim(floor - 4, 102)
        ax.set_yticks([floor] + list(range(0, 101, 20)))
        ax.set_yticklabels(["−∞"] + [str(v) for v in range(0, 101, 20)])
    for ax in axs[:, 0]:
        ax.set_ylabel("PLdB with rise-time proxy", fontsize=9, color=INK2)
    fig.suptitle("Ground carpets, finite-width shield at 20 km (half carpet; symmetric). "
                 "Points on the −∞ line: no Mark VII band audible.\n"
                 "Filled: R_m = 3 km (all azimuths); open circles: same config re-run at R_m = 8 km "
                 "(near-track azimuths; see R_m check)", fontsize=10, color=INK)
    fig.tight_layout(rect=(0, 0, 1, 0.96))
    fig.savefig(path, dpi=130)
    plt.close(fig)


def plot_signatures(C, path):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    lab = {c["key"]: c["label"] for c in C}
    sets = [("φ = 0 (undertrack)", 0.0), ("φ = 30°", 30.0), ("φ = 45°", 45.0)]
    keys = ["C2_h250_W250_He30", "C2_h250_W1000_He30", "C2_h250_W1000_H230", "C2_h250_W1000_He60lens"]
    fig, axs = plt.subplots(2, 3, figsize=(14, 7.5))
    for col, (title, ph) in enumerate(sets):
        for row, scale in enumerate(("full", "zoom")):
            ax = axs[row, col]
            _style(ax)
            g0 = _load("C2_h250_none", "ground") or {}
            r0 = g0.get((ph, "full", CROP_REL))
            if row == 0 and r0 is not None and "sig" in r0:
                ax.plot(r0["sig"][0], r0["sig"][1], color=INK, lw=1.0, label="unshielded")
            for i, k in enumerate(keys):
                if row == 1 and "W250" in k:
                    continue
                g = _load(k, "ground") or {}
                r = g.get((ph, "full", CROP_REL))
                if r is not None and "sig" in r:
                    ax.plot(r["sig"][0], r["sig"][1], color=COL[i], lw=1.2, label=lab[k])
            ax.set_xlim(-0.15, 1.4 if row else 0.4)
            if row == 0:
                ax.set_title(f"Concorde M2, h 250 m, {title}", fontsize=9.5, color=INK, loc="left")
            else:
                ax.set_title("W = 1 km only (zoom); late arrivals = edge-diffracted waves", fontsize=9,
                             color=INK2, loc="left")
            ax.set_xlabel("retarded time [s]", fontsize=9, color=INK2)
            ax.set_ylabel("ground pressure [Pa]", fontsize=9, color=INK2)
            if col == 0:
                ax.legend(fontsize=7.5, frameon=False)
    fig.tight_layout()
    fig.savefig(path, dpi=130)
    plt.close(fig)


# ----------------------------------------------------------------------------- main
def run(C, procs, only=None):
    t0 = time.time()
    order = [c for c in C if c["group"] == "baseline"]
    order += [c for c in C if c["group"] != "baseline"]
    for c in order:
        if only and c["key"] not in only:
            continue
        print(f"[{time.time() - t0:6.0f} s] {c['key']}: {c['label']}", flush=True)
        ground(c, procs)
        if c["key"] == "C2_h250_none":
            direct_check(c, procs)
    # band split at phi = 0 and at the loudest azimuth of each shielded main config
    for c in order:
        if c["group"] in ("baseline", "Rsens", "far") or (only and c["key"] not in only):
            continue
        rows = carpet_rows(c)
        if rows is None:
            continue
        pp = [_m(r, "pldb_proxy") for r in rows]
        idx = sorted({0, int(np.argmax(pp))})
        ground(c, procs, band=("lo", "hi"), idx=idx, what="bands")
    if not only:
        tie_check()
    # crop-threshold check
    for k in ("C2_h250_none", "C2_h250_W1000_H230"):
        c = next(c for c in C if c["key"] == k)
        if only and k not in only:
            continue
        ground(c, procs, idx=[0], rel=1e-4, what="cropcheck")
    print(f"total {time.time() - t0:.0f} s", flush=True)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--procs", type=int, default=4)
    ap.add_argument("--only", nargs="*", default=None)
    ap.add_argument("--report", action="store_true", help="only rebuild outputs from the cache")
    ap.add_argument("--list", action="store_true")
    a = ap.parse_args()
    C = build_configs()
    if a.list:
        for c in C:
            print(c["key"], "|", c["label"], "| base", c["base"])
        return
    if not a.report:
        run(C, a.procs, a.only)
    allrows, summ = write_outputs(C)
    print_tables(C, allrows, summ)
    band_report(C, allrows)
    checks_report(C, allrows)
    plot_carpets(C, allrows, os.path.join(RESULTS, "3d_carpets.png"))
    plot_signatures(C, os.path.join(RESULTS, "3d_signatures.png"))


if __name__ == "__main__":
    main()
