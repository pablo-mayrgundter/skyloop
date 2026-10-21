"""Wind and boundary-layer turbulence acting on the 3-D BoomShield ground carpets.

Re-uses the matching-circle spectra cached by run_3d.py (results/3d_cache/*_spec.pkl;
the cross-plane solver is NOT re-run) for Concorde M2, h = 250 m:
  unshielded (C2_h250_none), W = 1 km H2 30 m (C2_h250_W1000_H230), W = 500 m He 30 m
  (C2_h250_W500_He30), all at R_m = 3 km, and propagates them to the ground with
  boomsim.propagate (wind=...) exactly as run_3d does (age_from = aging radius, 50 kHz).

1. Wind: winter jet stream as head-, tail- and crosswind, QBO-type stratospheric easterly as
   headwind (boomsim.propagate profiles).  Carpet geometry (y_g, cut-off on both sides),
   ground PLdB per phi, Delta PLdB shielded - unshielded (headline = max of raw and Gate-0
   rise-time proxy, as in the 3-D study).  Calm values are the 3-D cache (identical
   pipeline; re-checked on one azimuth).
   Shield frame: the station-kept shield sees the ground-frame trace speed V_g = V + u_along
   (aircraft at fixed AIR Mach).  1-D estimate with shield.transfer at V_g, and an upper
   bound on the ground effect: circle spectra multiplied by T(V_g)/T(V) and propagated.
2. Turbulence: boomsim.turbulence phase screen (von Karman, PBL) at d = 500 m above the
   ground, Monte-Carlo, applied (A) to the deterministic ground signature (linear below the
   screen) and (B) to the signature at 500 m followed by the full nonlinear march to the
   ground (re-steepening), giving a bracket.  Same random kernels for shielded and
   unshielded (paired Delta PLdB).

Outputs: results/wind_turbulence.csv, results/wind_turbulence.png; tables printed for
results/wind_turbulence.md.  Per-job cache (resumable): $WT_CACHE or
results/wind_turbulence_cache/.  Serial (one process).

Run:  python boom-shield/sim/run_wind_turbulence.py [--part wind|frame|turb|report|all]
"""
from __future__ import annotations

import argparse
import csv
import math
import os
import pickle
import sys
import time

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import run_3d as R3  # noqa: E402
from boomsim import crossplane as cp  # noqa: E402
from boomsim import propagate as PR  # noqa: E402
from boomsim import shield as SH  # noqa: E402
from boomsim import turbulence as TU  # noqa: E402
from boomsim.loudness import metrics as wave_metrics  # noqa: E402
from boomsim.proxy import scored_metrics  # noqa: E402

RESULTS = os.path.join(HERE, "results")
CACHE = os.environ.get("WT_CACHE", os.path.join(RESULTS, "wind_turbulence_cache"))
KEYS = {"none": "C2_h250_none", "W1000": "C2_h250_W1000_H230", "W500": "C2_h250_W500_He30"}
LABEL = {"none": "unshielded", "W1000": "W 1 km, H2 30 m", "W500": "W 500 m, He 30 m"}
CASES = {
    "jet_head": lambda: PR.jet_stream_winter(180.0),
    "jet_tail": lambda: PR.jet_stream_winter(0.0),
    "jet_cross": lambda: PR.jet_stream_winter(90.0),
    "qbo_head": lambda: PR.qbo_easterly(180.0),
}
CASE_LABEL = {"calm": "calm", "jet_head": "jet, headwind", "jet_tail": "jet, tailwind",
              "jet_cross": "jet, crosswind (+y)", "qbo_head": "QBO easterly, headwind"}
PHI_SH = [0.0, 6.0, 16.0, 24.0, 36.0, 48.0]                 # shielded configs
PHI_W500 = [0.0, 16.0, 24.0, 48.0]          # W 500 m runs are slowest (cross: -24, -48 added)
PHI_UN = [0.0, 3.0, 6.0, 12.0, 16.0, 24.0, 30.0, 36.0, 42.0, 48.0, 54.0]
TURB_PHI = [0.0, 24.0, 48.0]
Z_SCREEN = 500.0
N_TURB = {"moderate": 200, "weak": 100, "strong": 100}
N_TURB_NL = 100
SEED = 20261021


# ------------------------------------------------------------------- cache helpers
def _p(name):
    return os.path.join(CACHE, name + ".pkl")


def load(name):
    if os.path.exists(_p(name)):
        with open(_p(name), "rb") as fh:
            return pickle.load(fh)
    return None


def save(name, obj):
    os.makedirs(CACHE, exist_ok=True)
    with open(_p(name) + ".tmp", "wb") as fh:
        pickle.dump(obj, fh)
    os.replace(_p(name) + ".tmp", _p(name))


_CFG = {}


def config(cfgkey):
    if cfgkey not in _CFG:
        C = {c["key"]: c for c in R3.build_configs()}
        c = C[KEYS[cfgkey]]
        sp = R3._load(c["key"], "spec")
        g = R3._load(c["key"], "ground")
        cfg = R3.make_cfg(c)
        cfg.fs = R3.FS
        _CFG[cfgkey] = (c, sp, g, cfg)
    return _CFG[cfgkey]


def circle_signature(cfgkey, phi_deg, spec_factor=None):
    """Cropped matching-circle signature (as run_3d).  ``spec_factor(f, j)`` -> (factor,
    additive) modifies the circle spectrum column j: P -> factor P + additive."""
    c, sp, g, cfg = config(cfgkey)
    j = int(np.flatnonzero(np.isclose(sp["phis"], abs(phi_deg)))[0])
    P = sp["P"][:, [j]]
    if spec_factor is not None:
        fac, add = spec_factor(sp["f"], j)
        P = P * fac[:, None] + add[:, None]
    t, p = cp.spectra_to_signatures(cfg, sp["f"], P)
    tc, pc = R3.crop(t, p[:, 0])
    phi = math.radians(phi_deg)
    return cfg, tc, pc, float(cp.aging_radius(cfg, abs(phi)))


def hp_stats(t, p, f_c=20.0):
    """Peak of the > f_c Hz part (the audible, edge/shock content) and max slope [Pa/ms]."""
    P = np.fft.rfft(p)
    f = np.fft.rfftfreq(p.size, t[1] - t[0])
    P[f < f_c] = 0.0
    ph = np.fft.irfft(P, p.size)
    dpdt = np.diff(p) / (t[1] - t[0]) / 1e3
    return float(np.max(np.abs(ph))), float(np.max(dpdt))


def score(t, p):
    m = {k: float(v) for k, v in scored_metrics(t, p).items()}
    m["hp_peak"], m["max_slope"] = hp_stats(t, p)
    return m


# ------------------------------------------------------------------- 1. wind
def wind_job(cfgkey, case, phi_deg, spec_factor=None, tag=None):
    name = f"wind_{cfgkey}_{tag or case}_{phi_deg:+g}"
    row = load(name)
    if row is not None:
        return row
    wind = None if case == "calm" else CASES[case]()
    cfg, tc, pc, age = circle_signature(cfgkey, phi_deg, spec_factor)
    t0 = time.time()
    info = {}
    tg, pg, yg = PR.propagate_from_circle(cfg.aircraft, math.radians(phi_deg), tc, pc,
                                          cfg.matching_radius(), z_ref=R3.H_S, age_from=age,
                                          wind=wind, info=info)
    row = dict(cfg=cfgkey, case=tag or case, phi=phi_deg, y_g=float(yg), wall=0.0)
    if tg is None:
        row.update(cutoff=True, metrics=None)
    else:
        row.update(cutoff=False, metrics=score(tg, pg), arrival=info.get("arrival"))
        if phi_deg in (0.0, 24.0):
            s = (tg >= tg[np.argmax(np.abs(pg))] - 0.4) & (tg <= tg[np.argmax(np.abs(pg))] + 1.4)
            row["sig"] = (tg[s][::5].astype(np.float32), pg[s][::5].astype(np.float32))
    row["wall"] = time.time() - t0
    save(name, row)
    print(f"  wind {name}: {row['wall']:.0f} s, y_g = {yg / 1e3:.2f} km, "
          f"PLdB = {row['metrics']['pldb'] if row['metrics'] else float('nan'):.1f}", flush=True)
    return row


def calm_row(cfgkey, phi_deg):
    c, sp, g, cfg = config(cfgkey)
    r = g[(abs(float(phi_deg)), "full", R3.CROP_REL)]
    return dict(cfg=cfgkey, case="calm", phi=phi_deg, y_g=math.copysign(r["y_g"], phi_deg or 1.0),
                cutoff=r["metrics"] is None, metrics=r["metrics"])


def phis_for(cfgkey, case):
    base = {"none": PHI_UN, "W1000": PHI_SH, "W500": PHI_W500}[cfgkey]
    if cfgkey == "none" and case == "qbo_head":
        base = PHI_SH
    if case == "jet_cross":
        neg = [p for p in base if p > 0] if cfgkey != "W500" else [24.0, 48.0]
        return [-p for p in neg[::-1]] + base
    return base


def run_wind():
    # verify the calm pipeline reproduces the 3-D cache on one azimuth
    chk = wind_job("W1000", "calm", 16.0)
    ref = calm_row("W1000", 16.0)["metrics"]
    print(f"calm check W1000 phi 16: PLdB {chk['metrics']['pldb']:.4f} vs cache {ref['pldb']:.4f}")
    order = [("none", c) for c in ("jet_head", "jet_tail", "jet_cross")] + \
            [("W1000", c) for c in ("jet_head", "jet_tail", "jet_cross")] + \
            [("W500", c) for c in ("jet_head", "jet_tail", "jet_cross")] + \
            [(k, "qbo_head") for k in ("none", "W1000", "W500")]
    for cfgkey, case in order:
        for ph in phis_for(cfgkey, case):
            wind_job(cfgkey, case, ph)


def geometry_table(ac):
    """Carpet geometry per wind case (ray tracing only)."""
    rows = []
    for case in ["calm"] + list(CASES):
        w = None if case == "calm" else CASES[case]()
        if w is None:
            pc, yc = PR.lateral_cutoff(ac.altitude, ac.V)
            pp, yp, pn, yn = pc, yc, -pc, -yc
            Vg = ac.V
        else:
            pp, yp = PR.lateral_cutoff(ac.altitude, ac.V, wind=w, side=1)
            pn, yn = PR.lateral_cutoff(ac.altitude, ac.V, wind=w, side=-1)
            Vg = PR.RayPath(ac.altitude, ac.V, wind=w).V_ground
        yg = {}
        for ph in (0.0, 24.0, 48.0, -24.0, -48.0):
            r = PR.RayPath(ac.altitude, ac.V, phi=math.radians(ph), wind=w, allow_cutoff=True)
            a = r.arrival(0.0)
            yg[ph] = (a["y"], a["x"], a["delay"])
        rows.append(dict(case=case, V_g=Vg, phi_cut_pos=math.degrees(pp), y_cut_pos=yp,
                         phi_cut_neg=math.degrees(pn), y_cut_neg=yn, yg=yg))
    return rows


# ------------------------------------------------------------------- shield frame
def shield_stack(cfgkey):
    c, sp, g, cfg = config(cfgkey)
    gas = cfg.gas
    return SH.Stack([SH.Membrane(float(cfg.m_top)), SH.GasLayer(gas, float(cfg.T)),
                     SH.Membrane(float(cfg.m_bot))]), gas


def frame_transfer(cfgkey, dV, mode="transmitted"):
    """Circle-spectrum modification for a ground-frame trace speed V + dV at the shield.

    mode "transmitted" (used): only the part of the field that goes THROUGH the gas
    changes, P -> P + P_inc (T(V+dV) - T(V)), with P_inc the unshielded circle spectrum
    at the same azimuth and T the 1-D stack transmission (insertion, re ambient air); the
    edge-diffracted part is unchanged.  mode "scale" (diagnostic, not physical): the whole
    shielded field multiplied by T(V+dV)/T(V)."""
    stack, gas = shield_stack(cfgkey)
    c, sp, g, cfg = config(cfgkey)
    V = cfg.aircraft.V
    atm = R3.ATM20
    P_inc = config("none")[1]["P"]

    def fac(f, j):
        f = np.asarray(f, float)
        T0, _ = SH.transfer(stack, np.maximum(f, 1e-3), V, atm, reference="air")
        T1, _ = SH.transfer(stack, np.maximum(f, 1e-3), V + dV, atm, reference="air")
        if mode == "scale":
            return T1 / T0, np.zeros(f.size, complex)
        return np.ones(f.size), P_inc[:, j] * (T1 - T0)
    return fac


def frame_table():
    rows = []
    bands = [(0.2, 2.0), (2.0, 5.0), (5.0, 10.0), (10.0, 20.0), (20.0, 50.0)]
    for cfgkey in ("W1000", "W500"):
        stack, gas = shield_stack(cfgkey)
        c, sp, g, cfg = config(cfgkey)
        V = cfg.aircraft.V
        for dV in (-60.0, -25.0, -15.0, 0.0, 15.0, 25.0, 60.0):
            f = np.linspace(0.2, 50.0, 2000)
            T, _ = SH.transfer(stack, f, V + dV, R3.ATM20, reference="air")
            kap = 2 * np.pi * np.sqrt(max(1 / (V + dV) ** 2 - 1 / gas.c ** 2, 0.0))
            rows.append(dict(cfg=cfgkey, gas=gas.kind, c_gas=gas.c, dV=dV, V_g=V + dV,
                             margin=gas.c / (V + dV), kappa_d_per_Hz=kap * float(cfg.T),  # Np/Hz
                             **{f"T_{a:g}_{b:g}": 10 * np.log10(np.mean(np.abs(T[(f >= a) & (f < b)]) ** 2))
                                for a, b in bands}))
    return rows


def run_frame():
    for cfgkey in ("W1000", "W500"):
        for dV, tag in ((-15.0, "frame_m15"), (15.0, "frame_p15"), (-60.0, "frame_m60"), (60.0, "frame_p60")):
            if cfgkey == "W500" and abs(dV) > 20:
                continue
            for ph in (0.0,):
                wind_job(cfgkey, "calm", ph, spec_factor=frame_transfer(cfgkey, dV, "scale"),
                         tag=tag + "scale")
                wind_job(cfgkey, "calm", ph, spec_factor=frame_transfer(cfgkey, dV), tag=tag)


# ------------------------------------------------------------------- 2. turbulence
def turb_prep(cfgkey, phi_deg):
    """Signature at the screen altitude (free field) + deterministic ground signature."""
    name = f"turbprep_{cfgkey}_{phi_deg:g}"
    d = load(name)
    if d is not None:
        return d
    cfg, tc, pc, age = circle_signature(cfgkey, phi_deg)
    phi = math.radians(phi_deg)
    t0 = time.time()
    ts, ps, _ = PR.propagate_from_circle(cfg.aircraft, phi, tc, pc, cfg.matching_radius(),
                                         z_ref=R3.H_S, age_from=age, z_ground=Z_SCREEN,
                                         reflection=1.0)
    ray = PR.RayPath(cfg.aircraft.altitude, cfg.aircraft.V, phi=phi, allow_cutoff=True)
    tg, pg = PR.propagate(ts, ps, ray, Z_SCREEN, 0.0)
    pg = 1.9 * pg
    d = dict(ts=ts, ps=ps, tg=tg, pg=pg, det=score(tg, pg), ray=None, wall=time.time() - t0)
    save(name, d)
    print(f"  turbprep {name}: {d['wall']:.0f} s, det PLdB {d['det']['pldb']:.2f} "
          f"(3-D cache {calm_row(cfgkey, phi_deg)['metrics']['pldb']:.2f})", flush=True)
    return d


def kernels(preset, n):
    name = f"kernels_{preset}_{n}"
    K = load(name)
    if K is not None:
        return K
    p = TU.PRESETS[preset]
    rng = np.random.default_rng(SEED + hash(preset) % 1000)
    K = []
    t0 = time.time()
    for i in range(n):
        k, i0, info = TU.realisation(p, R3.FS, rng)
        K.append((k.astype(np.float64), i0, info))
    save(name, K)
    print(f"  kernels {preset}: {n} in {time.time() - t0:.0f} s", flush=True)
    return K


def run_turb():
    preps = {(k, ph): turb_prep(k, ph) for k in KEYS for ph in TURB_PHI}
    for preset, n in N_TURB.items():
        K = kernels(preset, n)
        name = f"turbA_{preset}"
        if load(name) is None:
            res = {}
            for (k, ph), d in preps.items():
                res[(k, ph)] = [score(d["tg"], TU.apply_kernel(d["pg"], kk, i0)) for kk, i0, _ in K]
            save(name, res)
            print(f"  turb A {preset} done", flush=True)
    # (B) nonlinear below the screen, moderate turbulence
    K = kernels("moderate", N_TURB["moderate"])[:N_TURB_NL]
    for (k, ph), d in preps.items():
        name = f"turbB_moderate_{k}_{ph:g}"
        if load(name) is not None:
            continue
        ray = PR.RayPath(R3.aircraft("Concorde", 2.0, 250).altitude,
                         R3.aircraft("Concorde", 2.0, 250).V, phi=math.radians(ph), allow_cutoff=True)
        t0 = time.time()
        out = []
        for kk, i0, _ in K:
            ps = TU.apply_kernel(d["ps"], kk, i0)
            tg, pg = PR.propagate(d["ts"], ps, ray, Z_SCREEN, 0.0)
            out.append(score(tg, 1.9 * pg))
        save(name, out)
        print(f"  turb B {k} phi {ph:g}: {time.time() - t0:.0f} s", flush=True)


# ------------------------------------------------------------------- report
def _pc(x, q):
    x = np.asarray(x, float)
    return float(np.nanpercentile(x, q))


def headline(m, b):
    return max(m["pldb"] - b["pldb"], m["pldb_proxy"] - b["pldb_proxy"])


def wind_rows():
    out = {}
    for cfgkey in KEYS:
        for case in ["calm"] + list(CASES):
            for ph in (PHI_UN if case == "calm" else phis_for(cfgkey, case)):
                if case == "calm":
                    out[(cfgkey, case, ph)] = calm_row(cfgkey, ph)
                else:
                    r = load(f"wind_{cfgkey}_{case}_{ph:+g}")
                    if r is not None:
                        out[(cfgkey, case, ph)] = r
    return out


def report():
    ac = R3.aircraft("Concorde", 2.0, 250)
    W = wind_rows()
    csv_rows = []
    print("\n### carpet geometry (ray tracing)")
    print("| case | V_g m/s | cut-off +y: phi, y_c km | cut-off -y: phi, y_c km | half-width km (Δ) "
          "| y_g(24°) km | y_g(48°) km | y_g(-48°) km | x_g(0) km |")
    print("|---|---|---|---|---|---|---|---|---|")
    geo = geometry_table(ac)
    hw0 = 0.5 * (geo[0]["y_cut_pos"] - geo[0]["y_cut_neg"])
    for g in geo:
        hw = 0.5 * (g["y_cut_pos"] - g["y_cut_neg"])
        print(f"| {CASE_LABEL[g['case']]} | {g['V_g']:.0f} | {g['phi_cut_pos']:.1f}°, {g['y_cut_pos'] / 1e3:.1f} | "
              f"{g['phi_cut_neg']:.1f}°, {g['y_cut_neg'] / 1e3:.1f} | {hw / 1e3:.1f} ({(hw - hw0) / 1e3:+.1f}) | "
              f"{g['yg'][24.0][0] / 1e3:.1f} | {g['yg'][48.0][0] / 1e3:.1f} | {g['yg'][-48.0][0] / 1e3:.1f} | "
              f"{g['yg'][0.0][1] / 1e3:.1f} |")
        csv_rows.append(dict(kind="geometry", config="", case=g["case"], phi_deg="", y_g_m="",
                             V_g=g["V_g"], y_cut_pos=g["y_cut_pos"], y_cut_neg=g["y_cut_neg"],
                             phi_cut_pos=g["phi_cut_pos"], phi_cut_neg=g["phi_cut_neg"]))

    print("\n### ground PLdB per azimuth (raw/proxy) and headline Δ vs unshielded in the same wind")
    for case in ["calm"] + list(CASES):
        for cfgkey in KEYS:
            phs = phis_for(cfgkey, case if case != "calm" else "jet_head")
            for ph in phs:
                r = W.get((cfgkey, case, ph))
                if r is None:
                    continue
                b = W.get(("none", case, ph))
                m = r["metrics"]
                d = (headline(m, b["metrics"]) if (m and b and b["metrics"] and cfgkey != "none") else float("nan"))
                csv_rows.append(dict(kind="wind", config=cfgkey, case=case, phi_deg=ph, y_g_m=r["y_g"],
                                     pldb=m["pldb"] if m else "", pldb_proxy=m["pldb_proxy"] if m else "",
                                     peak_pa=m["peak_pa"] if m else "", d_pldb_headline=d,
                                     cutoff=r["cutoff"]))

    def summ(cfgkey, case):
        phs = [p for p in phis_for(cfgkey, "jet_cross" if case == "jet_cross" else "jet_head")
               if (cfgkey, case, p) in W]
        rr = [W[(cfgkey, case, p)] for p in phs]
        ok = [r for r in rr if r["metrics"]]
        if not ok:
            return None
        ut = W.get((cfgkey, case, 0.0))
        mx = max(ok, key=lambda r: r["metrics"]["pldb"])
        return dict(ut=ut["metrics"] if ut else None, mx=mx, n=len(ok), phs=phs)

    print("\n| config | case | undertrack PLdB raw/proxy | ΔPLdB undertrack (headline) | max PLdB over computed φ (y km) "
          "| ΔPLdB max (headline) | Δ at φ = ±24°, ±48° |")
    print("|---|---|---|---|---|---|---|")
    for cfgkey in ("W1000", "W500"):
        for case in ["calm"] + list(CASES):
            s, b = summ(cfgkey, case), summ("none", case)
            if s is None or b is None:
                continue
            common = [p for p in s["phs"] if (("none", case, p) in W)]
            sm = max((W[(cfgkey, case, p)] for p in common if W[(cfgkey, case, p)]["metrics"]),
                     key=lambda r: r["metrics"]["pldb_proxy"])
            bm = max((W[("none", case, p)] for p in common if W[("none", case, p)]["metrics"]),
                     key=lambda r: r["metrics"]["pldb_proxy"])
            dmax = max(sm["metrics"]["pldb"] - bm["metrics"]["pldb"],
                       sm["metrics"]["pldb_proxy"] - bm["metrics"]["pldb_proxy"])
            dut = headline(s["ut"], b["ut"])
            side = []
            for p in (24.0, 48.0, -24.0, -48.0):
                r, rb = W.get((cfgkey, case, p)), W.get(("none", case, p))
                if r and rb and r["metrics"] and rb["metrics"]:
                    side.append(f"{p:+g}: {headline(r['metrics'], rb['metrics']):+.1f}")
            print(f"| {LABEL[cfgkey]} | {CASE_LABEL[case]} | {s['ut']['pldb']:.1f}/{s['ut']['pldb_proxy']:.1f} | "
                  f"{dut:+.1f} | {sm['metrics']['pldb']:.1f} ({sm['y_g'] / 1e3:.1f}) | {dmax:+.1f} | {'; '.join(side)} |")
            csv_rows.append(dict(kind="wind_summary", config=cfgkey, case=case, d_ut=dut, d_max=dmax))
    print("\n| unshielded | case | undertrack PLdB raw/proxy | max PLdB raw (y km) |")
    for case in ["calm"] + list(CASES):
        s = summ("none", case)
        if s:
            print(f"| unshielded | {CASE_LABEL[case]} | {s['ut']['pldb']:.1f}/{s['ut']['pldb_proxy']:.1f} | "
                  f"{s['mx']['metrics']['pldb']:.1f} ({s['mx']['y_g'] / 1e3:.1f}) |")

    print("\n### shield frame (1-D stack at V_g = V + dV), band |T|^2 in dB")
    fr = frame_table()
    for r in fr:
        print(f"| {r['cfg']} | {r['gas']} c={r['c_gas']:.0f} | dV {r['dV']:+.0f} | V_g {r['V_g']:.0f} | "
              f"c_gas/V_g {r['margin']:.2f} | κd/f {r['kappa_d_per_Hz']:.3f} | "
              + " | ".join(f"{r[k]:.1f}" for k in r if k.startswith("T_")) + " |")
        csv_rows.append(dict(kind="frame1d", config=r["cfg"], case=f"dV{r['dV']:+g}",
                             **{k: v for k, v in r.items() if k.startswith("T_")}))
    print("\nground effect of the frame change alone (spectra x T(V_g)/T(V), no wind), Δ vs calm:")
    for cfgkey in ("W1000", "W500"):
        for tag in ("frame_m15", "frame_p15", "frame_m60", "frame_p60", "frame_m15scale",
                    "frame_p15scale", "frame_m60scale", "frame_p60scale"):
            for ph in (0.0, 24.0):
                r = load(f"wind_{cfgkey}_{tag}_{ph:+g}")
                if r is None:
                    continue
                c0 = calm_row(cfgkey, ph)["metrics"]
                print(f"  {cfgkey} {tag} phi {ph:g}: dPLdB raw {r['metrics']['pldb'] - c0['pldb']:+.2f}, "
                      f"proxy {r['metrics']['pldb_proxy'] - c0['pldb_proxy']:+.2f}")
                csv_rows.append(dict(kind="frame_ground", config=cfgkey, case=tag, phi_deg=ph,
                                     pldb=r["metrics"]["pldb"], pldb_proxy=r["metrics"]["pldb_proxy"],
                                     d_pldb_raw=r["metrics"]["pldb"] - c0["pldb"],
                                     d_pldb_proxy=r["metrics"]["pldb_proxy"] - c0["pldb_proxy"]))

    print("\n### turbulence (PLdB median [P10, P90]); det raw / proxy; Δ paired vs unshielded")
    TA = {pre: load(f"turbA_{pre}") for pre in N_TURB}
    for pre in N_TURB:
        if TA[pre] is None:
            continue
        print(f"\n#### {pre} (A: linear below screen)")
        print("| config | φ | det PLdB raw / proxy | turb PLdB median [P10, P90] | ΔPLdB det raw / proxy / headline "
              "| ΔPLdB turb median [P10, P90] | peak ratio median [P10,P90] | rise ms det / median turb "
              "| >20 Hz peak ratio median | max slope ratio median |")
        print("|---|---|---|---|---|---|---|---|---|---|")
        for k in KEYS:
            for ph in TURB_PHI:
                d = load(f"turbprep_{k}_{ph:g}")
                b = load(f"turbprep_none_{ph:g}")
                L = TA[pre][(k, ph)]
                Lb = TA[pre][("none", ph)]
                pl = [x["pldb"] for x in L]
                dd = [x["pldb"] - y["pldb"] for x, y in zip(L, Lb)]
                det, bd = d["det"], b["det"]
                pk = [x["peak_pa"] / det["peak_pa"] for x in L]
                hp = [x["hp_peak"] / det["hp_peak"] for x in L]
                sl = [x["max_slope"] / det["max_slope"] for x in L]
                rise = [x["rise_time_10_90_ms"] for x in L]
                dh = headline(det, bd) if k != "none" else float("nan")
                print(f"| {LABEL[k]} | {ph:g} | {det['pldb']:.1f} / {det['pldb_proxy']:.1f} | "
                      f"{np.median(pl):.1f} [{_pc(pl, 10):.1f}, {_pc(pl, 90):.1f}] | "
                      f"{det['pldb'] - bd['pldb']:+.1f} / {det['pldb_proxy'] - bd['pldb_proxy']:+.1f} / {dh:+.1f} | "
                      f"{np.median(dd):+.1f} [{_pc(dd, 10):+.1f}, {_pc(dd, 90):+.1f}] | "
                      f"{np.median(pk):.2f} [{_pc(pk, 10):.2f}, {_pc(pk, 90):.2f}] | "
                      f"{det['rise_time_10_90_ms']:.2f} / {np.median(rise):.2f} | {np.median(hp):.2f} | {np.median(sl):.2f} |")
                csv_rows.append(dict(kind=f"turbA_{pre}", config=k, phi_deg=ph, pldb=det["pldb"],
                                     pldb_proxy=det["pldb_proxy"], turb_median=np.median(pl),
                                     turb_p10=_pc(pl, 10), turb_p90=_pc(pl, 90),
                                     d_turb_median=np.median(dd), d_turb_p10=_pc(dd, 10),
                                     d_turb_p90=_pc(dd, 90), d_det_raw=det["pldb"] - bd["pldb"],
                                     d_det_proxy=det["pldb_proxy"] - bd["pldb_proxy"]))
    print("\n#### moderate, B: nonlinear below the screen")
    for k in KEYS:
        for ph in TURB_PHI:
            L = load(f"turbB_moderate_{k}_{ph:g}")
            Lb = load(f"turbB_moderate_none_{ph:g}")
            if L is None or Lb is None:
                continue
            d = load(f"turbprep_{k}_{ph:g}")["det"]
            pl = [x["pldb"] for x in L]
            dd = [x["pldb"] - y["pldb"] for x, y in zip(L, Lb)]
            pk = [x["peak_pa"] / d["peak_pa"] for x in L]
            rise = [x["rise_time_10_90_ms"] for x in L]
            print(f"| {LABEL[k]} | {ph:g} | {d['pldb']:.1f} | {np.median(pl):.1f} [{_pc(pl, 10):.1f}, {_pc(pl, 90):.1f}] | "
                  f"{np.median(dd):+.1f} [{_pc(dd, 10):+.1f}, {_pc(dd, 90):+.1f}] | "
                  f"peak {np.median(pk):.2f} [{_pc(pk, 10):.2f}, {_pc(pk, 90):.2f}] | rise {np.median(rise):.2f} |")
            csv_rows.append(dict(kind="turbB_moderate", config=k, phi_deg=ph, turb_median=np.median(pl),
                                 turb_p10=_pc(pl, 10), turb_p90=_pc(pl, 90), d_turb_median=np.median(dd),
                                 d_turb_p10=_pc(dd, 10), d_turb_p90=_pc(dd, 90)))
    for pre in N_TURB:
        K = load(f"kernels_{pre}_{N_TURB[pre]}")
        if K:
            print(f"{pre}: sigma_mu2 {TU.PRESETS[pre].sigma_mu2:.2e}, screen rms delay inside r_turb "
                  f"{1e3 * np.median([i['dT_rms_inner'] for _, _, i in K]):.2f} ms (analytic infinite screen "
                  f"{1e3 * TU.delay_sigma(TU.PRESETS[pre]):.2f} ms)")

    first = ["kind", "config", "case", "phi_deg", "y_g_m"]
    cols = first + sorted({k for r in csv_rows for k in r} - set(first))
    with open(os.path.join(RESULTS, "wind_turbulence.csv"), "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols)
        w.writeheader()
        for r in csv_rows:
            w.writerow({k: (f"{v:.4g}" if isinstance(v, float) else v) for k, v in r.items()})
    plot(W, geo)


# ------------------------------------------------------------------- figure
COL = {"calm": "#52514e", "jet_tail": "#2a78d6", "jet_head": "#eb6834", "jet_cross": "#1baf7a",
       "qbo_head": "#6250d6"}


def plot(W, geo):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    ink, ink2, grid = "#0b0b0b", "#52514e", "#e4e3df"
    plt.rcParams.update({"font.size": 8.5, "axes.edgecolor": ink2, "axes.labelcolor": ink,
                         "xtick.color": ink2, "ytick.color": ink2, "axes.spines.top": False,
                         "axes.spines.right": False})
    fig, ax = plt.subplots(2, 3, figsize=(13.5, 8.2))
    # (a) wind profiles
    a = ax[0, 0]
    z = np.linspace(0, 21000, 400)
    a.plot(np.interp(z, np.array(PR.JET_KM) * 1e3, PR.JET_MS), z / 1e3, color=COL["jet_tail"], lw=2,
           label="winter jet (speed)")
    a.plot(-np.interp(z, np.array(PR.QBO_KM) * 1e3, PR.QBO_MS), z / 1e3, color=COL["qbo_head"], lw=2,
           label="QBO easterly (eastward comp.)")
    a.axhline(20.0, color=ink2, lw=0.8, ls=":")
    a.text(-28, 20.3, "shield 20 km", color=ink2, fontsize=7.5)
    a.set_xlabel("wind speed [m/s]")
    a.set_ylabel("altitude [km]")
    a.set_title("(a) wind profiles", loc="left", fontsize=9)
    a.legend(frameon=False, fontsize=7.5, loc="lower right")
    a.grid(color=grid, lw=0.6)
    # (b, c) ground PLdB (proxy) vs y_g: unshielded and W1000
    for axx, cfgkey, ttl in ((ax[0, 1], "none", "(b) unshielded"), (ax[0, 2], "W1000", "(c) W 1 km, H2 30 m")):
        for case in ["calm"] + list(CASES):
            phs = PHI_UN if case == "calm" else phis_for(cfgkey, case)
            pts = [(W[(cfgkey, case, p)]["y_g"], W[(cfgkey, case, p)]["metrics"]["pldb_proxy"])
                   for p in phs if (cfgkey, case, p) in W and W[(cfgkey, case, p)]["metrics"]]
            if case == "calm":
                pts += [(-y, v) for y, v in pts if y > 0]
            if not pts:
                continue
            pts.sort()
            y, v = np.array(pts).T
            axx.plot(y / 1e3, v, "-o", ms=3.5, lw=1.6, color=COL[case], label=CASE_LABEL[case])
        axx.set_xlabel("ground lateral offset y_g [km]")
        axx.set_ylabel("PLdB (rise-time proxy)")
        axx.set_title(ttl + ": ground PLdB vs y_g", loc="left", fontsize=9)
        axx.grid(color=grid, lw=0.6)
    ax[0, 1].legend(frameon=False, fontsize=7.5)
    # (d) headline delta vs y_g, both shields
    a = ax[1, 0]
    for cfgkey, ls in (("W1000", "-"), ("W500", "--")):
        for case in ["calm"] + list(CASES):
            phs = PHI_UN if case == "calm" else phis_for(cfgkey, case)
            pts = []
            for p in phs:
                r, b = W.get((cfgkey, case, p)), W.get(("none", case, p))
                if r and b and r["metrics"] and b["metrics"]:
                    pts.append((r["y_g"], headline(r["metrics"], b["metrics"])))
            if case == "calm":
                pts += [(-y, v) for y, v in pts if y > 0]
            if not pts:
                continue
            pts.sort()
            y, v = np.array(pts).T
            a.plot(y / 1e3, v, ls, marker="o", ms=3, lw=1.4, color=COL[case],
                   label=f"{'1 km' if cfgkey == 'W1000' else '500 m'}, {CASE_LABEL[case]}" if case in ("calm", "jet_cross") else None)
    a.set_xlabel("ground lateral offset y_g [km]")
    a.set_ylabel("ΔPLdB vs unshielded (headline)")
    a.set_title("(d) ΔPLdB, solid W 1 km / dashed W 500 m", loc="left", fontsize=9)
    a.grid(color=grid, lw=0.6)
    a.legend(frameon=False, fontsize=7)
    # (e) turbulence PLdB distributions
    a = ax[1, 1]
    TA = load("turbA_moderate")
    if TA is not None:
        pos, lab = 0, []
        for k in KEYS:
            for ph in TURB_PHI:
                d = load(f"turbprep_{k}_{ph:g}")["det"]
                L = [x["pldb"] - d["pldb"] for x in TA[(k, ph)]]
                B = load(f"turbB_moderate_{k}_{ph:g}")
                a.vlines(pos, _pc(L, 10), _pc(L, 90), color=COL["jet_tail"], lw=5, alpha=0.55)
                a.plot(pos, np.median(L), "o", color=COL["jet_tail"], ms=5)
                if B is not None:
                    LB = [x["pldb"] - d["pldb"] for x in B]
                    a.vlines(pos + 0.25, _pc(LB, 10), _pc(LB, 90), color=COL["jet_head"], lw=5, alpha=0.55)
                    a.plot(pos + 0.25, np.median(LB), "o", color=COL["jet_head"], ms=5)
                a.plot(pos - 0.25, d["pldb_proxy"] - d["pldb"], "x", color=ink, ms=7, mew=1.5)
                lab.append(f"{ {'none': 'un', 'W1000': '1k', 'W500': '500'}[k]} {ph:g}°")
                pos += 1
        a.set_xticks(range(pos))
        a.set_xticklabels(lab, rotation=60, fontsize=7)
        a.axhline(0.0, color=ink2, lw=0.8)
        a.set_ylabel("PLdB − deterministic raw PLdB [dB]")
        a.set_title("(e) moderate PBL turbulence: P10–P90 bar, median dot", loc="left", fontsize=9)
        from matplotlib.lines import Line2D
        a.legend(handles=[Line2D([], [], color=COL["jet_tail"], lw=5, alpha=.55, label="linear below screen"),
                          Line2D([], [], color=COL["jet_head"], lw=5, alpha=.55, label="nonlinear below screen"),
                          Line2D([], [], color=ink, marker="x", ls="", mew=1.5, label="Gate-0 rise-time proxy")],
                 frameon=False, fontsize=7, loc="lower right")
        a.set_ylim(-13.0, 4.5)
        a.grid(color=grid, lw=0.6, axis="y")
    # (f) example signatures with turbulence
    a = ax[1, 2]
    K = load(f"kernels_moderate_{N_TURB['moderate']}")
    if K is not None:
        for k, off, col in (("none", 0.0, COL["calm"]), ("W1000", 0.0, COL["jet_tail"])):
            d = load(f"turbprep_{k}_0")
            tg, pg = d["tg"], d["pg"]
            if k == "none":
                i = int(np.argmax(np.diff(pg)))
                win = (tg > tg[i] - 0.004) & (tg < tg[i] + 0.012)
                scale, lbl = 1.0, "unshielded front shock"
            else:
                P = np.fft.rfft(pg)
                f = np.fft.rfftfreq(pg.size, tg[1] - tg[0])
                ph_ = np.fft.irfft(np.where(f >= 20, P, 0), pg.size)
                i = int(np.argmax(np.diff(ph_)))
                win = (tg > tg[i] - 0.010) & (tg < tg[i] + 0.012)
                scale, lbl = 10.0, "1 km shield: leak pulse ×10"
            t0 = tg[i]
            for j in range(4):
                kk, i0, _ = K[j]
                pt = TU.apply_kernel(pg, kk, i0)
                a.plot(1e3 * (tg[win] - t0), scale * pt[win], color=col, lw=0.8, alpha=0.5)
            a.plot(1e3 * (tg[win] - t0), scale * pg[win], color=ink if k == "none" else col, lw=2.0, label=lbl)
        a.set_xlabel("time re steepest compression [ms]")
        a.set_ylabel("overpressure [Pa]")
        a.set_title("(f) φ = 0: deterministic (thick) and 4 turbulent realisations", loc="left", fontsize=9)
        a.legend(frameon=False, fontsize=7)
        a.grid(color=grid, lw=0.6)
    fig.tight_layout()
    fig.savefig(os.path.join(RESULTS, "wind_turbulence.png"), dpi=150)
    print("figure written")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--part", default="all", choices=["wind", "frame", "turb", "report", "all"])
    a = ap.parse_args()
    if a.part in ("wind", "all"):
        run_wind()
    if a.part in ("turb", "all"):
        run_turb()
    if a.part in ("frame", "all"):
        run_frame()
    if a.part in ("report", "all"):
        report()


if __name__ == "__main__":
    main()
