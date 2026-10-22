"""Focused sonic booms ("superbooms") from accelerating and turning flight, with BoomShield.

Question: where do caustics from unsteady flight (Mach-cutoff transition during
acceleration, turns, corridor entry/exit) reach the ground, how strong is the focused
boom against the proposed 0.11 psf limit, and does the shield (H2 layer at 20 km, aircraft
250 m above) suppress it?

Pipeline per focus point (boomsim.focus):
  ray caustic search on the ground (J = 0 of the (t_e, phi) -> ground map) ->
  local fold parameters (L_c, R, tilt) from the two traced post-focus arrivals ->
  incident = steady level-flight signature on the same ray (Gate-0 chain: Whitham near
  field -> augmented Burgers; with the shield: -> 20 km -> shield.apply at the ray's
  trace speed 1/|s_h| -> ground) ->
  focus = linear Airy fold solution (raw and Gate-0 rise-time proxy) and nonlinear NKZ
  (nonlinear Tricomi) march -> x1.9 ground reflection -> peak psf, PLdB.

Phases:  python run_focus.py --phase validate|scenarios|points|report [--procs 2]
Cache: results/focus_cache/ (one pickle per scenario and per focus point; resumable).
Outputs: results/focus_booms.csv, results/focus_booms.png (focus_booms.md by hand).
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
import traceback
from multiprocessing import get_context

import numpy as np
from scipy.signal import resample_poly

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from boomsim import focus as F  # noqa: E402
from boomsim import propagate as PR  # noqa: E402
from boomsim import shield as S  # noqa: E402
from boomsim.atmosphere import std_atm  # noqa: E402
from boomsim.loudness import metrics  # noqa: E402
from boomsim.proxy import rise_time_proxy  # noqa: E402
from boomsim.source import CONCORDE, Aircraft  # noqa: E402
import run_psf_sweep as PS  # noqa: E402

RESULTS = os.path.join(HERE, "results")
CACHE = os.path.join(RESULTS, "focus_cache")
PSF = F.PSF
LIMIT_PA = 0.11 * PSF
FS = 50e3
FS_FOCUS = 10e3            # rate for the focus computations (rise times >= 2 ms)
FS_LO = 5e3
H_S = 20000.0              # shield top
Z_AC = H_S + PS.H_AC       # 20.25 km
REFL = 1.9

OVERTURE = PS.OVERTURE
# F-16 (validation against Downing et al. 1998): length 15.03 m; test weight ~9.8 t
# (clean, part fuel; assumed); equivalent volume 25 m^3 (S_max ~ 2.8 m^2, assumed);
# lift on the wing ~0.35 L ... 0.85 L.
F16 = Aircraft("F-16", length=15.03, mach=1.2, altitude=3048.0, weight=9800 * 9.81,
               volume=25.0, lift_start=0.35, lift_end=0.85)
PRESETS = {"Overture": OVERTURE, "Concorde": CONCORDE, "F-16": F16}

# shield stacks (floating lid, PS.make_stack conventions)
STACKS = {"none": None, "H2_55": ("hydrogen", 55.0), "H2_30": ("hydrogen", 30.0)}


# ----------------------------------------------------------------------------- cache
def _path(key):
    return os.path.join(CACHE, f"{key}.pkl")


def load(key):
    p = _path(key)
    if os.path.exists(p):
        with open(p, "rb") as fh:
            return pickle.load(fh)
    return None


def save(key, obj):
    os.makedirs(CACHE, exist_ok=True)
    p = _path(key)
    with open(p + ".tmp", "wb") as fh:
        pickle.dump(obj, fh)
    os.replace(p + ".tmp", p)


# ----------------------------------------------------------------------------- scenarios
def scenarios():
    """Scenario definitions: key -> dict(ac, kind, params, phis_deg (points to evaluate),
    shielded (aircraft above the shield?))."""
    S_ = {}
    for a in (0.02, 0.05, 0.10, 0.15):
        S_[f"OV_acc{a:g}"] = dict(ac="Overture", kind="level", z=Z_AC, m0=0.95, m1=1.7, a=a,
                                  pts=[0, 15, 30], shield=True,
                                  label=f"Overture level accel M0.95-1.7, {a:g} g, 20.25 km")
    S_["OV_climb"] = dict(ac="Overture", kind="climb", z0=12000.0, z1=Z_AC, m0=0.95, m1=1.7,
                          a=0.05, pts=[0, 15, 30], shield=False,
                          label="Overture climb 12->20.25 km + accel M0.95-1.7, 0.05 g")
    S_["OV_climb_hi"] = dict(ac="Overture", kind="climb", z0=17000.0, z1=Z_AC, m0=0.95,
                             m1=1.7, a=0.05, pts=[0], shield=False,
                             label="Overture climb 17->20.25 km + accel, 0.05 g")
    for b in (5, 10, 15, 20, 25):
        S_[f"OV_turn{b}"] = dict(ac="Overture", kind="turn", z=Z_AC, mach=1.7, bank=b,
                                 roll_in=10.0, pts=["steady"], shield=True,
                                 label=f"Overture M1.7 level turn, bank {b} deg")
    for a in (-0.05, -0.15):
        S_[f"OV_dec{-a:g}"] = dict(ac="Overture", kind="level", z=Z_AC, m0=1.7, m1=0.95, a=a,
                                   pts=[], shield=True,
                                   label=f"Overture level decel M1.7-0.95, {a:g} g")
    S_["CC_acc0.05"] = dict(ac="Concorde", kind="level", z=Z_AC, m0=0.95, m1=2.0, a=0.05,
                            pts=[0, 15, 30], shield=True,
                            label="Concorde level accel M0.95-2.0, 0.05 g, 20.25 km")
    S_["CC_acc0.15"] = dict(ac="Concorde", kind="level", z=Z_AC, m0=0.95, m1=2.0, a=0.15,
                            pts=[0], shield=True,
                            label="Concorde level accel M0.95-2.0, 0.15 g, 20.25 km")
    S_["CC_climb"] = dict(ac="Concorde", kind="climb", z0=8500.0, z1=16800.0, m0=0.95,
                          m1=2.0, a=0.04, pts=[0, 15], shield=False,
                          label="Concorde-like departure: climb 8.5->16.8 km + accel, 0.04 g")
    for a in (0.1, 0.2, 0.3):
        S_[f"F16_acc{a:g}"] = dict(ac="F-16", kind="level", z=3048.0, m0=0.95, m1=1.25, a=a,
                                   pts=[0], shield=False,
                                   label=f"F-16 level accel at 10 kft, {a:g} g (validation)")
    return S_


def trajectory(sc):
    if sc["kind"] == "level":
        return F.Trajectory.level_accel(sc["z"], sc["m0"], sc["m1"], sc["a"], pad_s=40.0)
    if sc["kind"] == "climb":
        return F.Trajectory.climb_accel(sc["z0"], sc["z1"], sc["m0"], sc["m1"], sc["a"],
                                        pad_s=40.0)
    if sc["kind"] == "turn":
        return F.Trajectory.level_turn(sc["z"], sc["mach"], sc["bank"], 120.0, pad_s=40.0,
                                       roll_in_s=sc["roll_in"])
    raise ValueError(sc["kind"])


def run_scenario(key, sc):
    got = load("scen_" + key)
    if got is not None:
        return got
    t0 = time.time()
    tr = trajectory(sc)
    col = tr.col
    out = dict(key=key, sc=sc)
    if sc["kind"] == "turn":
        phis = np.radians(np.r_[np.arange(-52.9, -40, 0.1), np.arange(-40, 40.1, 2.0),
                                np.arange(40.1, 53.0, 0.1)])
        te = np.arange(-25.0, 110.0, 1.0)
        res = F.scan(tr, te, phis, along_phi=True)
    else:
        phis = np.radians(np.arange(0.0, 72.1, 2.5) if sc["kind"] != "turn" else [])
        te = np.arange(tr.t0 + 1.0, tr.t1 - 1.0, 0.5)
        res = F.scan(tr, te, phis)
    out["roots"] = res["roots"]
    out["roots_phi"] = res.get("roots_phi", [])
    out["edges"] = res["edges"]
    out["J_min"] = float(np.nanmin(res["J"])) if np.isfinite(res["J"]).any() else np.nan
    out["J_max"] = float(np.nanmax(res["J"])) if np.isfinite(res["J"]).any() else np.nan
    out["n_ok"] = int(np.isfinite(res["J"]).sum())
    # focus line: ground points, emission Mach, shield crossing points
    fl = []
    for (t_, ph), kind in [(r, "te") for r in res["roots"]] + \
                          [(r, "phi") for r in out["roots_phi"]]:
        r0 = F.rays(tr, t_, ph, 0.0)
        st = r0["state"]
        rs = F.rays(tr, t_, ph, H_S) if float(st["z"]) > H_S else None
        fl.append(dict(te=t_, phi=ph, axis=kind, X=float(r0["P"][0]), Y=float(r0["P"][1]),
                       M=float(st["M"]), z_a=float(st["z"]), x_a=float(st["x"]),
                       y_a=float(st["y"]), psi=float(st["psi"]),
                       xs=float(rs["P"][0]) if rs is not None else np.nan,
                       ys=float(rs["P"][1]) if rs is not None else np.nan))
    out["focus_line"] = fl
    # cut-off: first/last time any ray reaches the ground (undertrack) and Mach there
    m_cut = col.cutoff_mach(float(tr.state(tr.t0)["z"]))
    out["M_cut_start"] = m_cut
    ed0 = [e for e in res["edges"] if abs(e[1]) < 1e-9]
    out["edge_phi0"] = [(e[0], float(tr.state(e[0])["M"]), float(tr.state(e[0])["x"]))
                        for e in ed0]
    # ground-reaching rays: earliest shield crossing (all phi) for shield-start rules
    if float(tr.state(tr.t0)["z"]) > H_S and sc["kind"] == "level":
        okm = np.isfinite(res["J"])
        if okm.any():
            T, Pp = np.meshgrid(te, phis)
            rs = F.rays(tr, T[okm], Pp[okm], H_S)
            st = rs["state"]
            out["shield_cross"] = dict(xs=rs["P"][0], ys=rs["P"][1], te=T[okm], phi=Pp[okm],
                                       x_a=st["x"], M=st["M"])
    # undertrack caustic altitude vs emission Mach (acceleration: where the caustic lies
    # relative to the shield)
    if sc["kind"] in ("level", "climb") and sc.get("a", 0) > 0:
        zc = []
        t_m1 = tr.time_at_mach(1.0005)
        t_end = res["roots"][0][0] + 30.0 if res["roots"] else tr.t1 - 1
        for t_ in np.unique(np.r_[t_m1 + np.geomspace(0.05, max(1.0, t_end - t_m1), 40)]):
            try:
                z_ = F.caustic_altitude(tr, t_, 0.0, z_lo=-900.0)
            except Exception:
                z_ = np.nan
            zc.append((float(t_), float(tr.state(t_)["M"]), z_))
        out["caustic_z"] = zc
    # deceleration / turn diagnostics: max geometric amplification on the ground grid
    if sc["kind"] in ("level", "turn") and sc.get("a", 1) < 0:
        okm = np.isfinite(res["J"])
        T, Pp = np.meshgrid(te, phis)
        g = F.focus_geometry(tr, T[okm], Pp[okm], 0.0)
        G = np.where(g["ok"], g["G"], np.nan)
        out["G_max"] = float(np.nanmax(G))
        out["G_max_at"] = (float(T[okm][np.nanargmax(G)]), float(Pp[okm][np.nanargmax(G)]))
    out["runtime_s"] = time.time() - t0
    save("scen_" + key, out)
    return out


def pick_points(key, sc, scen):
    """Focus points to evaluate in a scenario: list of (te, phi)."""
    pts = []
    if sc["kind"] == "turn":
        rp = [r for r in scen["roots_phi"] if r[0] > 15.0]
        if rp:
            # steady turn: fold at fixed phi; take the one at mid time
            rp.sort()
            pts.append(rp[len(rp) // 2])
        return pts
    roots = scen["roots"]
    for ph in sc["pts"]:
        cand = [r for r in roots if abs(math.degrees(r[1]) - ph) < 1e-6]
        if cand:
            pts.append(min(cand))          # first fold in time
    return pts


# ----------------------------------------------------------------------------- incident
def incident(acn, M, z_a, phi_ray, phi_lift, stack_key):
    """Free-field (no ground reflection) steady level-flight ground signature on the ray
    (M, z_a, phi_ray), lift seen at azimuth phi_lift, optionally through the shield.
    Returns dict(t, p) at FS_FOCUS (unshielded) or FS_LO (shielded) and diagnostics."""
    key = (f"inc_{acn}_M{M:.6f}_z{z_a:.1f}_p{math.degrees(phi_ray):.4f}"
           f"_l{math.degrees(phi_lift):.4f}_{stack_key}")
    got = load(key)
    if got is not None:
        return got
    ac = dataclasses.replace(PRESETS[acn], mach=float(M), altitude=float(z_a))
    r0 = min(1.5 * ac.length, 0.8 * PS.H_AC)
    t, p = PR.near_field(ac, r0, FS, phi=float(phi_lift))
    ray = PR.RayPath(ac.altitude, ac.V, z_min=0.0, phi=float(phi_ray))
    if ray.cutoff:
        raise ValueError("incident ray cut off")
    z0 = ac.altitude - r0 * math.cos(phi_ray)
    st = STACKS[stack_key]
    out = dict(key=key, acn=acn, M=M, z_a=z_a, phi=phi_ray, phi_lift=phi_lift,
               stack=stack_key)
    if st is None:
        tg, pg = PR.propagate(t, p, ray, z0, 0.0)
        q = int(round(FS / FS_FOCUS))
        pd = resample_poly(pg, 1, q)
        out.update(t=tg[0] + np.arange(pd.size) * q / FS, p=pd, fs=FS_FOCUS)
        m = metrics(tg, REFL * pg)
        out.update(peak_pa=m["peak_pa"] * 1.0, pldb=m["pldb"],
                   pldb_proxy=metrics(tg, rise_time_proxy(tg, REFL * pg))["pldb"],
                   rise_ms=m["rise_time_10_90_ms"], dur_ms=m["duration_ms"])
    else:
        fill, d = st
        t, p = PR.propagate(t, p, ray, z0, H_S)
        stk, g, lift, top, bot, margin = PS.make_stack(fill, d, "float", H_S)
        s_h = math.hypot(1.0 / ac.V, ray.sy)
        V_tr = 1.0 / s_h
        atmS = std_atm(H_S)
        a = math.sqrt(max(0.0, 1 / V_tr ** 2 - 1 / g.c ** 2))
        tau = PS.lf_tau(fill, d, "float", V_tr / float(atmS.c), H_S) if V_tr < g.c else 0.0
        pad = min(12.0, max(1.0, 15.0 * a * d, 8.0 * tau))
        tb, pb = PS._pad_taper(t, p, pad, lambda tt, pp: S.apply(stk, tt, pp, V_tr, atmS))
        q = int(round(FS / FS_LO))
        pp = resample_poly(pb, 1, q)
        tt = tb[0] + np.arange(pp.size) * q / FS
        tg, pg = PR.propagate(tt, pp, ray, H_S, 0.0)
        out.update(t=tg, p=pg, fs=FS_LO, V_tr=V_tr, c_gas=g.c, tau=tau, pad=pad)
        m = metrics(tg, REFL * pg)
        out.update(peak_pa=m["peak_pa"], pldb=m["pldb"], pldb_proxy=m["pldb"],
                   rise_ms=m["rise_time_10_90_ms"], dur_ms=m["duration_ms"])
    save(key, out)
    return out


# ----------------------------------------------------------------------------- focus point
def focus_amplitude(inc, fp, nkz=True):
    """Airy (linear) and NKZ (nonlinear) focus of a free-field incident signature."""
    t, p = inc["t"], inc["p"]
    res = {}
    lin = F.airy_focus(t, p, fp["L_c"], fp["R"], fp["c"])
    res["lin_peak"] = REFL * lin["peak"]
    res["lin_xi"] = lin["xi_peak"]
    m = metrics(lin["t"], REFL * lin["p"])
    res["lin_pldb"] = m["pldb"]
    if inc["stack"] == "none":
        pr = rise_time_proxy(t, p)
        linp = F.airy_focus(t, pr, fp["L_c"], fp["R"], fp["c"])
        res["linproxy_peak"] = REFL * linp["peak"]
        res["linproxy_pldb"] = metrics(linp["t"], REFL * linp["p"])["pldb"]
    # loudness/width of the amplified zone: ground distance (normal to the focus line)
    # over which the linear peak exceeds the limit / 2x the incident
    sa = math.sin(fp["alpha"]) if fp["alpha"] > 0 else np.nan
    xi, pk = lin["xi"], REFL * lin["peaks"]
    over = xi[pk > LIMIT_PA]
    res["zone_over_limit_m"] = float((over.max() - over.min()) / sa) if over.size else 0.0
    if nkz:
        # nonlinear march on the rise-time-proxy incident (turbulence-like 2-5 ms shocks;
        # identical to the raw one for the smooth shielded swell)
        q = rise_time_proxy(t, p) if inc["stack"] == "none" else p
        n = F.nkz_focus(t, q, fp["L_c"], fp["R"], fp["c"], fp["rho"], X_end=6.0,
                        record_from=5.0)
        res.update(nl_peak=REFL * n["peak_nl"], nl_lin_grid=REFL * n["peak_lin"],
                   nl_ratio=n["ratio"], mu=n["mu"], nl_hist=n["hist"], nl_T0=n["T0"],
                   nl_delta=n["delta"], nl_Gamma=n["Gamma"], rise_ms=1e3 * F.shock_rise(t, q))
        # PLdB of the linear (proxy) waveform rescaled to the nonlinear peak (approximate)
        lw = linp if inc["stack"] == "none" else lin
        sc = n["peak_nl"] / lw["peak"]
        res["nl_pldb"] = metrics(lw["t"], REFL * sc * lw["p"])["pldb"]
        res["focus_peak"] = res["nl_peak"]
        res["focus_pldb"] = res["nl_pldb"]
    else:
        res["focus_peak"] = res["lin_peak"]
        res["focus_pldb"] = res["lin_pldb"]
    return res


def run_point(job):
    key, scen_key, te, phi = job
    got = load(key)
    if got is not None:
        return got
    t0 = time.time()
    sc = scenarios()[scen_key]
    tr = trajectory(sc)
    try:
        fp = F.fold_parameters(tr, te, phi, 0.0)
    except Exception as e:
        out = dict(key=key, scen=scen_key, te=te, phi=phi, error=f"fold: {e}")
        save(key, out)
        return out
    out = dict(key=key, scen=scen_key, label=sc["label"], ac=sc["ac"], fold=fp,
               te=te, phi=phi)
    phi_lift = fp["phi_eff"] + fp["bank"]
    variants = ["none"]
    if sc["shield"] and fp["z_a"] > H_S:
        variants += ["H2_55"] + (["H2_30"] if sc["ac"] == "Overture" else [])
    # steady cruise reference on the same azimuth, for context
    for sk in variants:
        try:
            inc = incident(sc["ac"], fp["M_eff"], fp["z_a"], fp["phi_eff"], phi_lift, sk)
        except Exception as e:
            out[sk] = dict(error=f"incident: {e}")
            continue
        r = dict(ss_peak=inc["peak_pa"], ss_pldb=inc["pldb"], ss_rise_ms=inc["rise_ms"],
                 ss_dur_ms=inc["dur_ms"], V_tr=inc.get("V_tr", np.nan))
        try:
            r.update(focus_amplitude(inc, fp, nkz=True))
        except Exception as e:
            r["error"] = f"focus: {e}\n{traceback.format_exc()}"
        out[sk] = r
    # shield crossing of the focus ray
    if fp["z_a"] > H_S:
        rs = F.rays(tr, te, phi, H_S)
        out["xs"], out["ys"] = float(rs["P"][0]), float(rs["P"][1])
    out["runtime_s"] = time.time() - t0
    save(key, out)
    print(f"  {key}: {out['runtime_s']:.0f} s", flush=True)
    return out


def point_jobs():
    jobs = []
    for key, sc in scenarios().items():
        scen = load("scen_" + key)
        if scen is None:
            continue
        for te, ph in pick_points(key, sc, scen):
            jobs.append((f"pt_{key}_phi{math.degrees(ph):.2f}", key, te, ph))
    return jobs


# ----------------------------------------------------------------------------- validation
def validate():
    """Uniform-atmosphere analytic checks + F-16 steady boom (Downing et al. 1998 quote
    5.0 psf predicted for steady M1.2 at 10,000 ft)."""
    out = {}
    from boomsim.propagate import uniform_atmosphere
    col = F.Column(atmosphere=uniform_atmosphere(0.0))
    c = float(col.c[0])
    rows = []
    for a in (0.5, 1.0, 2.0):
        tr = F.Trajectory.level_accel(15000.0, 0.95, 1.7, a / F.G0, column=col)
        sc = F.scan(tr, np.arange(tr.t0 + 1, tr.t1 - 1, 0.5), [0.0])
        te = sc["roots"][0][0]
        M = float(tr.state(te)["M"])
        rows.append(dict(a=a, M_model=M, M_analytic=math.sqrt(1 + (a * 15000 / c ** 2) ** (2 / 3))))
    out["uniform_Mf"] = rows
    print("uniform atmosphere focus Mach (beta^3 = aH/c^2):", rows)
    ac = dataclasses.replace(F16, mach=1.2, altitude=3048.0)
    t, p = PR.ground_signature(ac)
    m = metrics(t, p)
    out["F16_steady_M1.2"] = dict(peak_psf=m["peak_pa"] / PSF, pldb=m["pldb"])
    print("F-16 steady M1.2 10 kft:", out["F16_steady_M1.2"])
    # Overture / Concorde steady cruise references at 20.25 km
    for acn, M in (("Overture", 1.7), ("Concorde", 2.0)):
        ac = dataclasses.replace(PRESETS[acn], mach=M, altitude=Z_AC)
        t, p = PR.ground_signature(ac, r0=min(1.5 * ac.length, 0.8 * PS.H_AC))
        m = metrics(t, p)
        out[f"{acn}_cruise"] = dict(peak_psf=m["peak_pa"] / PSF, pldb=m["pldb"])
        print(acn, "cruise", out[f"{acn}_cruise"])
    # Guiraud / Gill-Seebass scaling check of the NKZ solver: step-like incident,
    # peak ~ amplitude^(4/5) in the strongly nonlinear limit
    tt = np.arange(0, 0.6, 1e-4)
    rise = 2e-3
    step = (0.5 * (1 + np.tanh((tt - 0.1) / (rise / 2.197)))
            * np.exp(-np.clip(tt - 0.1, 0, None) / 0.1))
    gs = []
    for amp in (1.0, 4.0, 16.0, 64.0, 256.0, 1024.0):
        n = F.nkz_focus(tt, amp * step, 1000.0, 50e3, 340.0, 1.2, X_end=4.0, record_from=3.0)
        lin = F.nkz_focus(tt, amp * step, 1000.0, 50e3, 340.0, 1.2, X_end=0.1, linear=True)
        gs.append(dict(amp=amp, peak=n["peak_nl"], lin=lin["peak_lin"], mu=n["mu"],
                       A=n["amp"] * amp, ratio=n["peak_nl"] / lin["peak_lin"]))
        print("  NKZ step", gs[-1], flush=True)
    a_ = np.log([g["amp"] for g in gs])
    p_ = np.log([g["peak"] for g in gs])
    out["guiraud_slopes"] = list(np.diff(p_) / np.diff(a_))
    print("  d ln(peak)/d ln(amplitude):", np.round(out["guiraud_slopes"], 3))
    out["guiraud"] = gs
    save("validate", out)
    return out


def collect():
    rows = []
    for j in point_jobs():
        r = load(j[0])
        if r is None:
            continue
        fp = r.get("fold")
        base = dict(point=j[0], scenario=r["scen"], label=r.get("label", ""),
                    aircraft=r.get("ac", ""), phi_deg=math.degrees(r["phi"]))
        if fp is None:
            rows.append(dict(base, error=r.get("error", "")))
            continue
        base.update(M_emit=fp["M"], M_eff=fp["M_eff"], z_a_m=fp["z_a"], X_km=fp["X"] / 1e3,
                    Y_km=fp["Y"] / 1e3, x_aircraft_km=fp["x_a"] / 1e3,
                    L_c_m=fp["L_c"], R_km=fp["R"] / 1e3, alpha_deg=math.degrees(fp["alpha"]),
                    xs_shield_km=r.get("xs", np.nan) / 1e3, ys_shield_km=r.get("ys", np.nan) / 1e3)
        for sk in ("none", "H2_55", "H2_30"):
            v = r.get(sk)
            if v is None:
                continue
            row = dict(base, shield=sk, error=v.get("error", "").split("\n")[0])
            for k in ("ss_peak", "ss_pldb", "ss_rise_ms", "ss_dur_ms", "V_tr", "lin_peak",
                      "linproxy_peak", "lin_pldb", "nl_peak", "nl_lin_grid", "nl_ratio", "mu",
                      "nl_pldb", "focus_peak", "focus_pldb", "zone_over_limit_m"):
                if k in v:
                    row[k] = v[k]
            if "focus_peak" in v:
                # unshielded: nonlinear NKZ peak; shielded (mu ~ 1e-4, linear): the full-band
                # Airy peak if larger (NKZ is band-limited); the ground maximum in the focus
                # zone is at least the steady value on the same ray
                fpk = v["focus_peak"] if sk == "none" else max(v["focus_peak"], v["lin_peak"])
                row["focus_peak"] = fpk
                row["focus_psf"] = fpk / PSF
                row["ss_psf"] = v["ss_peak"] / PSF
                row["focus_factor"] = fpk / v["ss_peak"]
                row["ground_max_psf"] = max(fpk, v["ss_peak"]) / PSF
                row["pass_0.11psf"] = bool(max(fpk, v["ss_peak"]) <= LIMIT_PA)
            rows.append(row)
    return rows


def write_csv(rows, path):
    keys = []
    for r in rows:
        for k in r:
            if k not in keys:
                keys.append(k)
    with open(path, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=keys)
        w.writeheader()
        for r in rows:
            w.writerow({k: (f"{v:.6g}" if isinstance(v, float) else v) for k, v in r.items()})


def print_tables(rows):
    print("\nFocus points (ground peak incl. x1.9; psf):")
    print(f"{'scenario':16s} {'phi':>5s} {'M':>6s} {'R km':>7s} {'L_c m':>7s} {'shield':>6s} "
          f"{'steady':>7s} {'lin':>7s} {'NL':>7s} {'F':>5s} {'PLdB ss':>7s} {'PLdB f':>7s}")
    for r in rows:
        if "focus_peak" not in r:
            print(r["scenario"], r.get("shield", ""), r.get("error", ""))
            continue
        print(f"{r['scenario']:16s} {r['phi_deg']:5.1f} {r['M_emit']:6.3f} {r['R_km']:7.1f} "
              f"{r['L_c_m']:7.0f} {r['shield']:>6s} {r['ss_psf']:7.3f} "
              f"{r['lin_peak'] / PSF:7.3f} {r['focus_psf']:7.3f} {r['focus_factor']:5.2f} "
              f"{r['ss_pldb']:7.1f} {r['focus_pldb']:7.1f}")


def report():
    rows = collect()
    write_csv(rows, os.path.join(RESULTS, "focus_booms.csv"))
    print_tables(rows)
    scen = {k: load("scen_" + k) for k in scenarios()}
    print("\nScenario geometry:")
    for k, d in scen.items():
        if d is None:
            continue
        fl = [x for x in d["focus_line"] if x["axis"] == "te"]
        fp_ = [x for x in d["focus_line"] if x["axis"] == "phi"]
        msg = f"{k:14s} roots {len(fl):3d}/{len(fp_):3d}"
        if fl:
            u = min(fl, key=lambda x: abs(x["phi"]))
            phmax = max(abs(math.degrees(x["phi"])) for x in fl)
            ymax = max(abs(x["Y"]) for x in fl)
            msg += (f"  undertrack M_f {u['M']:.4f} X {u['X'] / 1e3:.1f} km (aircraft at "
                    f"{u['x_a'] / 1e3:.1f} km), focus line to |phi| {phmax:.1f} deg, "
                    f"|Y| {ymax / 1e3:.1f} km")
            if np.isfinite(u["xs"]):
                msg += f"; shield crossing {(u['xs'] - u['x_a']) / 1e3:+.2f} km ahead"
                ys = max(abs(x["ys"]) for x in fl)
                msg += f", max |y| at shield {ys:.0f} m"
        if fp_:
            msg += (f"  steady-turn fold phi {np.median([math.degrees(x['phi']) for x in fp_]):.2f}"
                    f" deg, ground |Y| {np.median([math.hypot(x['X'] - x['x_a'], x['Y'] - x['y_a']) for x in fp_]) / 1e3:.1f} km from aircraft")
        if "G_max" in d:
            msg += f"  max G {d['G_max']:.3f}"
        if d.get("edge_phi0"):
            msg += "  cut-off (phi 0) at M " + ", ".join(f"{e[1]:.4f}" for e in d["edge_phi0"])
        print(msg)
    plot(rows, scen, os.path.join(RESULTS, "focus_booms.png"))


def plot(rows, scen, path):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(2, 2, figsize=(12, 9))
    C = dict(none="#c0392b", H2_55="#2471a3", H2_30="#17a589")
    # (a) undertrack caustic altitude vs emission Mach
    a = ax[0, 0]
    for k, ls in (("OV_acc0.02", ":"), ("OV_acc0.05", "-"), ("OV_acc0.15", "--")):
        d = scen.get(k)
        if d and d.get("caustic_z"):
            z = np.array([(m, zz) for _, m, zz in d["caustic_z"]])
            a.plot(z[:, 0], z[:, 1] / 1e3, ls, color="k", label=k.replace("OV_acc", "Overture ") + " g")
    a.axhline(20.0, color=C["H2_55"], lw=3, alpha=0.5, label="shield top (20 km)")
    a.axhline(20.25, color="0.5", lw=0.8, label="aircraft (20.25 km)")
    a.axhline(0.0, color="0.3", lw=0.8)
    mc = F.default_column().cutoff_mach(Z_AC)
    a.axvline(mc, color="0.5", ls=":", lw=0.8)
    a.text(mc + 0.003, 10, f"cut-off M {mc:.3f}", rotation=90, fontsize=8, color="0.4")
    a.set_xlim(1.0, 1.3)
    a.set_ylim(-1, 21)
    a.set_xlabel("emission Mach")
    a.set_ylabel("caustic altitude on the undertrack ray [km]")
    a.set_title("(a) Where the acceleration caustic lies")
    a.legend(fontsize=7, loc="lower left")
    # (b) focus lines on the ground
    b = ax[0, 1]
    for k, col in (("OV_acc0.02", "#7d3c98"), ("OV_acc0.05", "#c0392b"), ("OV_acc0.15", "#d68910"),
                   ("CC_acc0.05", "#2e4053")):
        d = scen.get(k)
        if not d:
            continue
        fl = sorted([x for x in d["focus_line"] if x["axis"] == "te"], key=lambda x: x["phi"])
        X = np.array([x["X"] for x in fl]) / 1e3
        Y = np.array([x["Y"] for x in fl]) / 1e3
        xa = np.array([x["x_a"] for x in fl]) / 1e3
        b.plot(np.r_[X[::-1], X], np.r_[-Y[::-1], Y], "-", color=col, label=k)
        b.plot(xa[0], 0, "o", color=col, ms=4)
    b.set_xlabel("along track x [km] (acceleration starts at x = 0)")
    b.set_ylabel("y [km]")
    b.set_title("(b) Ground focus lines (dot: aircraft at emission, phi = 0)")
    b.legend(fontsize=7)
    b.grid(alpha=0.3)
    # (c) focus peak vs acceleration
    c = ax[1, 0]
    acc = [0.02, 0.05, 0.1, 0.15]
    for sk, lab in (("none", "no shield"), ("H2_55", "H2 55 m"), ("H2_30", "H2 30 m")):
        ys, ysteady = [], []
        for a_ in acc:
            rr = [r for r in rows if r["scenario"] == f"OV_acc{a_:g}" and r.get("shield") == sk
                  and abs(r["phi_deg"]) < 1e-6 and "focus_psf" in r]
            ys.append(rr[0]["focus_psf"] if rr else np.nan)
            ysteady.append(rr[0]["ss_psf"] if rr else np.nan)
        c.semilogy(acc, ys, "o-", color=C[sk], label=f"focus, {lab}")
        c.semilogy(acc, ysteady, "s:", color=C[sk], ms=3, alpha=0.7, label=f"steady at M_f, {lab}")
    c.axhline(0.11, color="k", ls="--", lw=1)
    c.text(0.021, 0.12, "0.11 psf", fontsize=8)
    c.set_xlabel("level acceleration [g] (Overture, 20.25 km, undertrack)")
    c.set_ylabel("ground peak [psf]")
    c.set_title("(c) Focus peak vs acceleration")
    c.legend(fontsize=7)
    c.grid(alpha=0.3, which="both")
    # (d) waveforms: unshielded steady vs focus; shielded focus
    d_ = ax[1, 1]
    r0 = load("pt_OV_acc0.05_phi0.00")
    if r0 is not None and "fold" in r0:
        fp = r0["fold"]
        for sk, lab in (("none", "no shield"), ("H2_55", "H2 55 m")):
            try:
                inc = incident("Overture", fp["M_eff"], fp["z_a"], fp["phi_eff"], fp["phi_eff"] + fp["bank"], sk)
            except Exception:
                continue
            lin = F.airy_focus(inc["t"], inc["p"], fp["L_c"], fp["R"], fp["c"])
            sc_ = r0[sk].get("focus_peak", REFL * lin["peak"]) / (REFL * lin["peak"])
            i0 = int(np.argmax(inc["p"]))
            d_.plot(inc["t"] - inc["t"][i0], REFL * inc["p"], color=C[sk], lw=0.8, alpha=0.6,
                    label=f"steady, {lab}")
            j0 = int(np.argmax(lin["p"]))
            d_.plot(lin["t"] - lin["t"][j0], REFL * sc_ * lin["p"], color=C[sk], lw=1.5,
                    label=f"focus, {lab}")
        d_.axhline(LIMIT_PA, color="k", ls="--", lw=1)
        d_.set_xlim(-0.6, 1.2)
        d_.set_xlabel("time [s]")
        d_.set_ylabel("ground overpressure [Pa]")
        d_.set_title("(d) Overture 0.05 g undertrack: steady vs focus")
        d_.legend(fontsize=7)
        d_.grid(alpha=0.3)
    fig.tight_layout()
    fig.savefig(path, dpi=130)
    print("wrote", path)


# ----------------------------------------------------------------------------- extras
def _turn_has_fold(bank, mach=1.7, z=Z_AC, roll_in=10.0, steady_only=True):
    tr = F.Trajectory.level_turn(z, mach, bank, 60.0, pad_s=30.0, roll_in_s=roll_in)
    phis = np.radians(np.r_[np.arange(-53.0, -30.0, 0.05), np.arange(30.0, 53.0, 0.05)])
    te = np.array([20.0, 40.0]) if steady_only else np.arange(-roll_in - 2, 40.0, 0.5)
    sc = F.scan(tr, te, phis, along_phi=True)
    return bool(sc["roots_phi"]) or (not steady_only and bool(sc["roots"]))


def _levelof_has_fold(level_off_s, acn_z0=17000.0):
    tr = F.Trajectory.climb_accel(acn_z0, Z_AC, 0.95, 1.7, 0.05, pad_s=40.0,
                                  level_off_s=level_off_s)
    te = np.arange(tr.t0 + 1.0, tr.t1 - 1.0, 0.5)
    sc = F.scan(tr, te, np.radians([0.0, 20.0, 40.0]))
    return any(float(tr.state(t)["M"]) > 1.6 for t, _ in sc["roots"])


def extras():
    out = load("extras") or {}
    # 1. turn: bank threshold for a ground caustic (steady turn and with a 10 s roll-in)
    if "turn_threshold" not in out:
        res = {}
        for mach in (1.7, 1.4):
            for steady in (True, False):
                lo, hi = 5.0, 40.0
                if not _turn_has_fold(hi, mach, steady_only=steady):
                    res[(mach, steady)] = float("nan")
                    continue
                for _ in range(10):
                    mid = 0.5 * (lo + hi)
                    if _turn_has_fold(mid, mach, steady_only=steady):
                        hi = mid
                    else:
                        lo = mid
                res[(mach, steady)] = hi
                print(f"turn M{mach} {'steady' if steady else 'with 10 s roll-in'}: "
                      f"ground caustic above bank {hi:.1f} deg", flush=True)
        out["turn_threshold"] = res
        save("extras", out)
    # 2. top-of-climb level-off (push-over): shortest level-off without a ground caustic
    if "leveloff" not in out:
        lo, hi = 0.0, 120.0
        if _levelof_has_fold(hi):
            out["leveloff"] = float("nan")
        else:
            for _ in range(8):
                mid = 0.5 * (lo + hi)
                if _levelof_has_fold(mid):
                    lo = mid
                else:
                    hi = mid
            out["leveloff"] = hi
        print(f"level-off at M1.7: ground caustic for level-off shorter than ~{out['leveloff']:.0f} s "
              f"(Gaussian sigma = T/4)", flush=True)
        save("extras", out)
    # 3. distances along track: M = 1, cut-off, focus emission (acceleration / deceleration)
    if "distances" not in out:
        rows = []
        col = F.default_column()
        mcut = col.cutoff_mach(Z_AC)
        for key, sc in scenarios().items():
            if sc["kind"] != "level" or sc["z"] != Z_AC:
                continue
            tr = trajectory(sc)
            d = load("scen_" + key)
            row = dict(key=key, a=sc["a"], M_cut=mcut)
            for lab, M in (("M1", 1.0), ("Mcut", mcut), ("M1.7", 1.7 - 1e-6)):
                if sc["a"] > 0:
                    t_ = tr.time_at_mach(M)
                else:
                    tt = np.linspace(tr.t0, tr.t1, 200001)
                    mm = tr.state(tt)["M"]
                    i = np.flatnonzero(mm <= M)
                    t_ = float(tt[i[0]]) if i.size else np.nan
                row[f"x_{lab}_km"] = float(tr.state(t_)["x"]) / 1e3 if np.isfinite(t_) else np.nan
            if d and d["roots"]:
                u = [r for r in d["roots"] if abs(r[1]) < 1e-9]
                if u:
                    row["x_focus_emit_km"] = float(tr.state(u[0][0])["x"]) / 1e3
                    row["X_focus_ground_km"] = float(F.rays(tr, u[0][0], 0.0, 0.0)["P"][0]) / 1e3
            if d and "shield_cross" in d:
                scx = d["shield_cross"]
                row["dx_cross_max_km"] = float(np.max(scx["xs"] - scx["x_a"])) / 1e3
                row["dx_cross_min_km"] = float(np.min(scx["xs"] - scx["x_a"])) / 1e3
                row["ys_max_m"] = float(np.max(np.abs(scx["ys"])))
                row["xs_first_km"] = float(np.min(scx["xs"])) / 1e3
                row["xs_last_km"] = float(np.max(scx["xs"])) / 1e3
            fl = [x for x in (d or {}).get("focus_line", []) if x["axis"] == "te"]
            if fl:
                row["focus_ys_max_m"] = float(max(abs(x["ys"]) for x in fl))
                row["focus_Y_max_km"] = float(max(abs(x["Y"]) for x in fl)) / 1e3
                row["focus_phi_max_deg"] = float(max(abs(math.degrees(x["phi"])) for x in fl))
            rows.append(row)
            print(row, flush=True)
        out["distances"] = rows
        save("extras", out)
    return out


# ----------------------------------------------------------------------------- on-ramps
TRUNK = os.path.join(HERE, "..", "network-maps", "trunk")
A_TRUNK = 0.23             # m/s^2, trunk acceleration (network-maps/trunk/build_trunk.py)
A_TRUNK_DEC = 0.46
EDGE_MARGIN_M2 = 450.0     # m inside the edge that passed in 3-D (Concorde M2, H2 60 m, 1.5 km)


def cutoff_with_wind():
    """Effective undertrack cut-off Mach at 20.25 km (ground-frame trace speed must
    exceed c + u everywhere below) for the repo's representative wind profiles."""
    z = np.linspace(0.0, Z_AC, 4000)
    c = np.asarray(std_atm(z).c)
    out = {}
    for name, w in (("calm", None), ("winter jet, tailwind", PR.jet_stream_winter(0.0)),
                    ("winter jet, headwind", PR.jet_stream_winter(180.0)),
                    ("QBO easterly, headwind", PR.qbo_easterly(180.0)),
                    ("QBO easterly, tailwind", PR.qbo_easterly(0.0))):
        u = np.zeros_like(z) if w is None else w(z)[0]
        out[name] = float((np.max(c + u) - u[-1]) / c[-1])
    return out


def lateral_reach(M, z_a=Z_AC):
    """Max |y| at the shield top (20 km) of steady level rays that reach the ground."""
    col = F.default_column()
    tr = F.Trajectory.level_accel(z_a, M, M + 1e-9, 1e-6, pad_s=5.0)
    phis = np.radians(np.arange(0.0, 89.0, 0.05))
    r = F.rays(tr, np.zeros_like(phis), phis, 0.0)
    ok = r["ok"]
    if not ok.any():
        return 0.0, 0.0
    ph = phis[ok].max()
    rs = F.rays(tr, 0.0, ph, H_S)
    return float(abs(rs["P"][1])), float(math.degrees(ph))


def accel_km(m0, m1, a=A_TRUNK):
    c = 295.1
    return ((m1 * c) ** 2 - (m0 * c) ** 2) / (2 * a) / 1e3


def onramp():
    out = load("onramp") or {}
    # 1. minimum H2 thickness vs Mach (undertrack, floating lid, 1-D), Overture and Concorde
    d_list = [5.0, 10.0, 20.0, 30.0, 45.0, 60.0, 80.0]
    for d in d_list:
        STACKS[f"H2_{d:g}"] = ("hydrogen", d)
    if "thick" not in out:
        rows = []
        for acn, Ms in (("Overture", (1.16, 1.2, 1.3, 1.5, 1.7, 2.0)),
                        ("Concorde", (1.16, 1.2, 1.3, 1.5, 2.0))):
            for M in Ms:
                pk = {}
                for d in d_list:
                    inc = incident(acn, M, Z_AC, 0.0, 0.0, f"H2_{d:g}")
                    pk[d] = inc["peak_pa"]
                    if inc["peak_pa"] < 0.3 * LIMIT_PA:
                        break
                ds = np.array(sorted(pk))
                ps = np.array([pk[d] for d in ds])
                def dmin(lim):
                    if ps[0] <= lim:
                        return float(ds[0])
                    i = np.flatnonzero(ps <= lim)
                    if i.size == 0:
                        return float("nan")
                    i = int(i[0])
                    x0, x1 = math.log(ps[i - 1]), math.log(ps[i])
                    return float(ds[i - 1] + (ds[i] - ds[i - 1]) * (math.log(lim) - x0) / (x1 - x0))
                inc55 = incident(acn, M, Z_AC, 0.0, 0.0, "H2_60") if 60.0 in pk else None
                row = dict(ac=acn, M=M, d_min=dmin(LIMIT_PA), d_min_focus=dmin(LIMIT_PA / 1.3),
                           peaks={float(k): float(v) for k, v in pk.items()})
                ys, phc = lateral_reach(M)
                row.update(ys_max=ys, phi_c=phc)
                # swell duration through the min-thickness layer (cross-plane wavelength scale)
                dd = min(pk, key=lambda d: abs(d - (row["d_min"] if np.isfinite(row["d_min"]) else 30)))
                incd = incident(acn, M, Z_AC, 0.0, 0.0, f"H2_{dd:g}")
                row["swell_ms"] = float(incd["dur_ms"])
                rows.append(row)
                print(acn, M, {k: (round(v, 1) if isinstance(v, float) else v) for k, v in row.items()
                               if k != "peaks"}, flush=True)
        out["thick"] = rows
        save("onramp", out)
    if "wind" not in out:
        out["wind"] = cutoff_with_wind()
        save("onramp", out)
    print("cut-off Mach with winds:", {k: round(v, 3) for k, v in out["wind"].items()})
    # 2. width scaling: rays cross within |ys|; edge margin scaled with the cross-plane
    #    wavelength of the swell, lambda_perp = V T_swell / beta (anchored at M2.0: 450 m)
    th = [r for r in out["thick"] if r["ac"] == "Overture"]
    def lam(r):
        V = r["M"] * 295.1
        return V * r["swell_ms"] / 1e3 / math.sqrt(r["M"] ** 2 - 1)
    # the 3-D width runs that passed (psf_sweep: 1.5 km, rays >= 450 m inside the edge)
    # were Concorde M2.0 with H2 60 m: anchor the margin there
    ref = [r for r in out["thick"] if r["ac"] == "Concorde" and abs(r["M"] - 2.0) < 1e-9][0]
    for r in out["thick"]:
        r["lambda_perp_m"] = lam(r)
        r["edge_margin_m"] = EDGE_MARGIN_M2 * lam(r) / lam(ref)
        r["W_req_m"] = 2 * (r["ys_max"] + r["edge_margin_m"])
    # 3. ramps of the M2.0 spine and shield km
    import csv as _csv
    with open(os.path.join(TRUNK, "ramps.csv")) as fh:
        ramps = [r for r in _csv.DictReader(fh) if r["tree"] == "greedy_M2.0"]
    mc = F.default_column().cutoff_mach(Z_AC)
    km = dict(
        accel_cut_to_M2=accel_km(mc, 2.0), accel_115_to_M2=accel_km(1.15, 2.0),
        onset_115_130=accel_km(1.15, 1.30), decel_M2_to_cut=accel_km(mc, 2.0, A_TRUNK_DEC),
        accel_cut_to_Mf=accel_km(mc, 1.165), accel_cut_to_M13=accel_km(mc, 1.3),
        wind_cut_range=(min(out["wind"].values()), max(out["wind"].values())),
        accel_low_to_high_cut=accel_km(min(out["wind"].values()), max(out["wind"].values())))
    out["km"] = km
    out["ramps"] = ramps
    print("km:", {k: (round(v, 1) if isinstance(v, float) else v) for k, v in km.items()})
    save("onramp", out)
    write_onramp_csv(out)
    plot_onramp(out)
    return out


def write_onramp_csv(out):
    rows = []
    for r in out["thick"]:
        rows.append(dict(table="thickness_width", aircraft=r["ac"], mach=r["M"],
                         h2_min_m_0p11psf=r["d_min"], h2_min_m_with_focus_x1p3=r["d_min_focus"],
                         rays_max_y_at_shield_m=r["ys_max"], lateral_cutoff_phi_deg=r["phi_c"],
                         swell_duration_ms=r["swell_ms"], lambda_perp_m=r["lambda_perp_m"],
                         edge_margin_m=r["edge_margin_m"], width_required_m=r["W_req_m"]))
    for k, v in out["wind"].items():
        rows.append(dict(table="cutoff_wind", case=k, cutoff_mach=v))
    for k, v in out["km"].items():
        rows.append(dict(table="distances_km", case=k, value=v if not isinstance(v, tuple) else "-".join(f"{x:.3f}" for x in v)))
    for r in out["ramps"]:
        rows.append(dict(table="ramps_M2_spine", case=r["nearest_metro"], km_from_metro=r["km_from_metro"],
                         lat=r["lat"], lon=r["lon"], joining_pax_M=r["joining_pax_M"],
                         share_of_joins=r["share_of_joins"]))
    write_csv(rows, os.path.join(RESULTS, "onramp_shield.csv"))


def plot_onramp(out):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(1, 3, figsize=(15, 4.6))
    a = ax[0]
    for acn, col in (("Overture", "#2471a3"), ("Concorde", "#c0392b")):
        rr = [r for r in out["thick"] if r["ac"] == acn]
        a.plot([r["M"] for r in rr], [r["d_min"] for r in rr], "o-", color=col, label=f"{acn}: 0.11 psf")
        a.plot([r["M"] for r in rr], [r["d_min_focus"] for r in rr], "s--", color=col, ms=3,
               label=f"{acn}: 0.11 psf / 1.3 (focus margin)")
    a.set_xlabel("Mach (level, 20.25 km)")
    a.set_ylabel("minimum H2 layer [m] (1-D undertrack)")
    a.set_title("(a) Layer needed along the ramp")
    a.grid(alpha=0.3)
    a.legend(fontsize=7)
    b = ax[1]
    for acn, col in (("Overture", "#2471a3"), ("Concorde", "#c0392b")):
        rr = [r for r in out["thick"] if r["ac"] == acn]
        b.plot([r["M"] for r in rr], [r["W_req_m"] / 1e3 for r in rr], "o-", color=col, lw=2,
               label=f"{acn}: required width")
        b.plot([r["M"] for r in rr], [2 * r["ys_max"] / 1e3 for r in rr], ":", color=col,
               label=f"{acn}: ray footprint at the shield (2 max|y|)")
    b.axhline(1.5, color="0.4", ls="--", lw=1)
    b.text(1.17, 1.53, "trunk baseline 1.5 km", fontsize=7, color="0.3")
    b.set_xlabel("Mach")
    b.set_ylabel("km")
    b.set_title("(b) Shield width: footprint + edge margin (scaling)")
    b.grid(alpha=0.3)
    b.legend(fontsize=7)
    c = ax[2]
    # merge schematic (not to scale across track): trunk along x, feeders at 0/30/60 deg
    # meeting the trunk axis at x = 0 at M <= day's cut-off, acceleration on the axis
    km = out["km"]
    W = 1.6
    c.fill_between([-30, 160], -W / 2, W / 2, color="#aed6f1", label="trunk shield, ~1.6 km wide")
    c.fill_between([-10, 0], -W / 2, W / 2, color="#1f618d", label="apron, 10 km before accel start")
    c.fill_between([0, km["accel_cut_to_Mf"]], -W / 2, W / 2, color="#d68910",
                   label=f"focus-emitting stretch (M_cut to M_f, {km['accel_cut_to_Mf']:.0f} km)")
    c.fill_between([0, km["onset_115_130"]], W / 2, W / 2 + 0.6, color="#f5b041", alpha=0.6,
                   label=f"M1.15-1.3 ({km['onset_115_130']:.0f} km)")
    for th_ in (30, 60):
        dx = 4.0 / math.tan(math.radians(th_))
        c.plot([-25 - dx, -25], [-4.0, 0.0], color="0.35", lw=1)
    c.plot([-32, 0], [0, 0], color="0.35", lw=1)
    c.text(-22, -2.6, "feeders at M <= day's cut-off:\nno ground boom, turn freely", fontsize=7)
    c.text(80, 2.6, f"M2.0 reached {km['accel_cut_to_M2']:.0f} km downstream\n(0.023 g; off chart)",
           fontsize=7)
    c.set_xlim(-32, 160)
    c.set_ylim(-4.5, 4.0)
    c.set_xlabel("km along trunk (accel start at 0)")
    c.set_ylabel("km across (exaggerated)")
    c.set_title("(c) Merge at cut-off Mach, then accelerate on axis")
    c.legend(fontsize=6, loc="lower right")
    fig.tight_layout()
    path = os.path.join(RESULTS, "onramp_shield.png")
    fig.savefig(path, dpi=130)
    print("wrote", path)


# ----------------------------------------------------------------------------- driver
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--phase", required=True,
                    choices=["validate", "scenarios", "points", "extras", "onramp", "report", "all"])
    ap.add_argument("--procs", type=int, default=2)
    ap.add_argument("--only", default="")
    args = ap.parse_args()
    t0 = time.time()
    if args.phase in ("scenarios", "all"):
        for key, sc in scenarios().items():
            if args.only and args.only not in key:
                continue
            r = run_scenario(key, sc)
            print(f"{key}: {len(r['roots'])} te-roots, {len(r['roots_phi'])} phi-roots, "
                  f"{r['runtime_s']:.0f} s", flush=True)
    if args.phase in ("points", "all"):
        jobs = [j for j in point_jobs() if load(j[0]) is None and args.only in j[0]]
        print(f"{len(jobs)} focus points to run", flush=True)
        if args.procs > 1:
            with get_context("fork").Pool(args.procs, maxtasksperchild=1) as pool:
                for _ in pool.imap_unordered(run_point, jobs):
                    pass
        else:
            for j in jobs:
                run_point(j)
    if args.phase == "validate":
        validate()
    if args.phase in ("extras", "all"):
        extras()
    if args.phase == "onramp":
        onramp()
    if args.phase in ("report", "all"):
        report()
    print(f"done in {time.time() - t0:.0f} s")


if __name__ == "__main__":
    main()
