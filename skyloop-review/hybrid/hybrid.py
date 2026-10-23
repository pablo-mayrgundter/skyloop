#!/usr/bin/env python3
"""Skyloop maglev-plus-wings hybrid: feasible map (screening level).

Pods ride an actively steered maglev guideway inside a 10-80 m low-pressure bore of the floating (20 km)
tube.  The linear motor gives contactless thrust and fine guidance; wings carry a fraction f of the weight
(0 <= f <= 1), unloading the levitation magnets.  f = 1 is flight.py's free-flying pod over a steered strip;
f = 0 is stability.md's maglev pod on an actuated guideway.

Run:   python3 skyloop-review/hybrid/hybrid.py             (~12 min, one process; grid and levers resume)
       python3 skyloop-review/hybrid/hybrid.py figs        (figures only, from the CSVs)
Out:   hybrid_suspension.csv  magnetic gap / suspension vs f: load, EDS gap, frequency, drag, aero damping,
                              and the guideway stroke each suspension needs (stability spectra)
       hybrid_tolerance.csv   allowable guideway error amplitude vs wavelength for each suspension (fixed guideway)
       hybrid_stroke.csv      guideway stroke and magnetic-gap travel needed vs Mach, bore, line tension
       hybrid_grid.csv        best wing loading per (f, pod, bore, Mach, pressure) (> 2 MB, gitignored)
       hybrid_band.csv        highest passing / marginal Mach per f x pod x bore, binding check above it
       hybrid_bore.csv        smallest passing bore vs Mach per f x pod
       hybrid_points.csv      best operating points (fastest, energy-optimal, capacity-optimal) per f
       hybrid_econ.csv        B/C hook (system.md / flight.md conventions) for the best point per f x bore x Mach
       hybrid_levers.csv      what moves the band
       hybrid_propulsion.csv  side check: onboard propulsion for free-flying pods
       hybrid_maps.png, hybrid_summary.png

Reuses read-only: flight/flight.py (stroke conventions, finish/bands logic, economics data, adders),
heating/heating.py (pod_state, tube_side, plate, q_stag, q_net, thresholds), stability/stability.py (line,
spectra, transfer, constants), regime/regime.py (pod, q_hole, beta_kantrowitz), system/sky_model.py
(headway_s, ride_min, run_cost, conventions).  Mach convention: V = M x 295.07 m/s.  (A) = assumption.
"""
from __future__ import annotations

import csv
import math
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
REV = os.path.dirname(HERE)
sys.dont_write_bytecode = True
for _d in ("flight", "heating", "stability"):
    sys.path.insert(0, os.path.join(REV, _d))
import flight as FL  # noqa: E402  (read-only; imports heating, stability, regime, physics, sky_model)

H, S, R, sm = FL.H, FL.S, FL.R, FL.sm
G, R_E, C_REF = FL.G, FL.R_E, FL.C_REF
RHO_AMB, P_AMB, T_AMB = FL.RHO_AMB, FL.P_AMB, FL.T_AMB
LF, PEAK_FACTOR = FL.LF, FL.PEAK_FACTOR
KX = S.KX

# ----------------------------------------------------------------------------- grids
F_LIST = [0.0, 0.25, 0.5, 0.75, 0.9, 1.0]
M_GRID = np.round(np.arange(2.0, 10.001, 0.25), 2)
P_GRID = np.geomspace(5.0, 1000.0, 20)
D_LIST = [10.0, 15.0, 20.0, 30.0, 40.0, 60.0, 80.0]
PODS = [8000.0, 17000.0, 30000.0, 47000.0]           # 20 / 50 / 93 / 150 seats (regime.pod, 300 kg/seat)
WS_LIST = FL.WS_LIST
M_STROKE = np.round(np.arange(2.0, 10.001, 0.5), 2)

# ----------------------------------------------------------------------------- base assumptions (A)
BASE = dict(
    FL.BASE,
    LD_mag=sm.P["LD_mag"],   # EDS magnetic lift/drag (system: 300; Inductrack ~200-300 at 500 km/h, rises with V)
    g0=0.10,                 # m, SC-EDS nominal gap at full magnetic load (SCMaglev ~0.10 m)
    g_min=0.03,              # m, closest approach kept in service
    tau_p=1.0,               # m, magnet pole pitch: lift ~ exp(-2 pi g / tau_p)
    rough=0.01,              # m, short-wave error of the cable-robot guideway (supports ~25 m) the gap absorbs
    guide_w=1.0,             # m, magnet/guideway width beyond the pod body (lateral room at f = 0)
    f_s=S.F_S, z_s=S.Z_S,    # pod suspension on the guideway: 1 Hz passive, zeta 0.3 (stability/flight)
)

CHECKS = ["lift", "stroke", "gap", "trip", "power", "motor", "breach", "blockage", "reflection", "crossover",
          "gas_dT", "throughput", "skin_T", "le_T", "heat_sink", "film_T", "guideway_T", "energy"]
SHORT = dict(FL.SHORT, gap="magnetic gap")
RANK = FL.RANK


def write_csv(name, rows, mode="w"):
    path = os.path.join(HERE, name)
    keys = []
    for r in rows:
        keys += [k for k in r if k not in keys]
    new = mode == "w" or not os.path.exists(path)
    with open(path, mode, newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=keys, extrasaction="ignore")
        if new:
            w.writeheader()
        for r in rows:
            w.writerow({k: (f"{v:.4g}" if isinstance(v, float) else v) for k, v in r.items()})
    print(f"wrote {name} ({len(rows)} rows)", flush=True)


def read_csv(name):
    path = os.path.join(HERE, name)
    if not os.path.exists(path):
        return []
    out = []
    for r in csv.DictReader(open(path)):
        rr = {}
        for k, v in r.items():
            try:
                rr[k] = float(v)
            except (TypeError, ValueError):
                rr[k] = v
        out.append(rr)
    return out


# ============================================================================= 1. gap and suspension
def eds_gap(f, F=BASE):
    """Gap at fixed magnets when wings carry f of the weight (EDS lift ~ exp(-2 pi g / tau_p))."""
    if f >= 1.0:
        return F["g_mag"]                                        # flight.py: thrust gap held by flight control
    return F["g0"] + F["tau_p"] / (2 * math.pi) * math.log(1.0 / (1.0 - f))


def eds_freq(f, F=BASE):
    """Passive EDS heave frequency: k = (2 pi / tau_p) (1 - f) W  ->  f_n = sqrt(2 pi g (1-f) / tau_p) / 2 pi."""
    return math.sqrt(2 * math.pi * G * max(1.0 - f, 0.0) / F["tau_p"]) / (2 * math.pi)


def gap_usable(f, F=BASE):
    return eds_gap(f, F) - F["g_min"]


def aero_zeta(f, M, CL, fn):
    """Heave damping ratio from the wing: c = L CL_alpha / (CL V); zeta = c / (2 m omega)."""
    if f <= 0 or fn <= 0:
        return 0.0
    V = M * C_REF
    cla = 4.0 / math.sqrt(M * M - 1)                             # supersonic flat plate
    return f * G * cla / (CL * V) / (2 * 2 * math.pi * fn)


def ride2(L, M, fm, zm, kind, s_gap, plane="y", cache=None, fixed=False):
    """Two-stage filter: steered guideway path (skyhook low-pass at fc, stability Z_C) then the pod's magnetic
    suspension (fm, zm, passive or skyhook).  Returns the smallest guideway stroke (4 sigma + scallop) whose
    cabin ride is <= 0.05 g and whose gap travel (3 sigma + scallop) fits s_gap.  fixed=True: no steering."""
    V = M * C_REF
    if cache is not None and plane in cache:
        Sy = cache[plane]
    else:
        Sy, _ = S.spectra(L, plane)
        if cache is not None:
            cache[plane] = Sy
    ks = 2 * math.pi / L["spacing"]
    f_var = L["f_var_y"] if plane == "y" else L["f_var_z"]
    y_st = L["eps"] * f_var / S.K_of(L, ks, plane)
    Ta_m, Tr_m = S.transfer(V, KX, fm, zm, kind)
    Ta_m0, Tr_m0 = S.transfer(V, ks, fm, zm, kind)
    fcs = [None] if fixed else list(np.geomspace(2e-5, 3.0, 110))
    cands = []
    for fc in fcs:
        if fc is None:
            Ta_c, Tr_c, Ta_c0, Tr_c0 = np.ones_like(KX), np.zeros_like(KX), 1.0, 0.0
        else:
            Ta_c, Tr_c = S.transfer(V, KX, fc, S.Z_C, "skyhook")
            Ta_c0, Tr_c0 = S.transfer(V, ks, fc, S.Z_C, "skyhook")
        a_pk = S.PK * math.sqrt(float(np.trapezoid(Sy * (V * KX) ** 4 * Ta_c * Ta_m, KX))) \
            + (V * ks) ** 2 * y_st * math.sqrt(Ta_c0 * Ta_m0)
        st = S.N_SIG_CLEAR * math.sqrt(float(np.trapezoid(Sy * Tr_c, KX))) + y_st * math.sqrt(Tr_c0)
        gp = 3 * math.sqrt(float(np.trapezoid(Sy * Ta_c * Tr_m, KX))) + y_st * math.sqrt(Ta_c0 * Tr_m0)
        cands.append(dict(a_pk=a_pk, stroke=st, gap=gp, fc=fc if fc else float("nan")))
    ok = [c for c in cands if c["a_pk"] <= S.A_LIM and c["gap"] <= s_gap]
    if ok:
        return dict(min(ok, key=lambda c: c["stroke"]), ok=1)
    return dict(min(cands, key=lambda c: c["a_pk"]), ok=0)


def make_line(D, T, eps):
    return S.line(D=D, payload_extra=20.0 + S.tension_mass(T), eps=eps, T=T)     # as flight.stroke_table


SUSP = [  # name, f, fm, zeta, kind, usable gap (m)
    ("flight / stability convention (1 Hz passive)", None, 1.0, 0.3, "passive", 1e9),
    ("rigid maglev: EMS 10 mm gap, 8 Hz", 0.0, 8.0, 0.5, "passive", 0.005),
    ("SC-EDS f=0: 0.10 m gap, 1.25 Hz", 0.0, None, 0.3, "passive", None),
    ("hybrid EDS f=0.5", 0.5, None, 0.3, "passive", None),
    ("hybrid EDS f=0.75", 0.75, None, 0.3, "passive", None),
    ("hybrid EDS f=0.9", 0.9, None, 0.3, "passive", None),
    ("active skyhook 1 Hz, 0.3 m gap", None, 1.0, 0.7, "skyhook", 0.30),
    ("active skyhook 0.3 Hz, 0.3 m gap", None, 0.3, 0.7, "skyhook", 0.30),
    ("active skyhook 0.3 Hz, 1 m gap", None, 0.3, 0.7, "skyhook", 1.0),
    ("soft 0.1 Hz, 3 m travel (pod as its own reference)", None, 0.1, 0.7, "skyhook", 3.0),
]


def _susp_params(s, F=BASE):
    name, f, fm, z, kind, sg = s
    if fm is None:
        fm = eds_freq(f, F)
    if sg is None:
        sg = gap_usable(f, F) - F["rough"]
    return name, f, fm, z, kind, sg


def stage_suspension(F=BASE):
    """Item 1: what the gap and the wing do to load, drag, frequency, damping, and to the stroke needed."""
    rows = []
    m, pax = 17000.0, 40.0
    for f in F_LIST:
        W = m * G
        fn = eds_freq(f, F)
        r = dict(kind="gap vs f (17 t, fixed magnets)", f=f, magnet_load_kN=(1 - f) * W / 1e3,
                 eds_gap_m=eds_gap(f, F), gap_usable_m=gap_usable(f, F), eds_freq_Hz=fn,
                 mag_drag_N=(1 - f) * W / F["LD_mag"],
                 mag_Wh_paxkm=(1 - f) * W / F["LD_mag"] / F["eta_d"] / pax / 3.6,
                 # EMS alternative: holding power ~1.7 kW/t at 10 mm (Transrapid, unverified), ~ load x gap^2
                 ems_hold_MW_at_0p1m=1.7e3 * (1 - f) * m / 1e3 * (0.1 / 0.01) ** 2 / 1e6)
        for M in (2.0, 5.0, 8.0):
            st = H.pod_state(M, 200.0, 400.0, dict(H.BASE, m=m))
            if st.get("lift_ok"):
                r[f"wing_Wh_paxkm_M{M:g}_200Pa"] = f * st["drag_N"] / F["eta_d"] / pax / 3.6
                r[f"aero_zeta_M{M:g}"] = aero_zeta(f, M, st["CL"], fn if fn > 0 else 1.0)
        rows.append(r)
    # guideway stroke needed with each suspension, 40 m bore, 0.1 / 0.3 GN (ride 0.05 g, gap travel within usable)
    for T in (0.1e9, 0.3e9):
        L = make_line(40.0, T, F["eps"])
        cy = {}
        for s in SUSP:
            name, f, fm, z, kind, sg = _susp_params(s, F)
            r = dict(kind="guideway stroke by suspension (40 m bore)", T_GN=T / 1e9, suspension=name, f_n_Hz=fm,
                     zeta=z, damping=kind, gap_usable_m=sg)
            for M in (2.0, 4.0, 6.0, 8.0, 10.0):
                rd = ride2(L, M, fm, z, kind, sg, "y", cy)
                r[f"stroke_M{M:g}_m"] = rd["stroke"] if rd["ok"] else float("nan")
                r[f"gap_M{M:g}_m"] = rd["gap"]
                fx = ride2(L, M, fm, z, kind, 1e9, "y", cy, fixed=True)
                r[f"fixed_guideway_a_g_M{M:g}"] = fx["a_pk"] / G
            rows.append(r)
            print(f"  susp {T/1e9:g} GN {name}: " + ", ".join(f"M{M:g} {r[f'stroke_M{M:g}_m']:.1f}" for M in (2, 6, 10)),
                  flush=True)
    write_csv("hybrid_suspension.csv", rows)


TOL_CASES = [  # name, primary (fm, zeta, kind, usable) or ("eds", f, zeta, kind), secondary (fs, zeta, stroke) or None
    ("rigid maglev: EMS 10 mm gap (30 Hz servo) + 1 Hz secondary", (30.0, 0.7, "skyhook", 0.005), (1.0, 0.3, 0.10)),
    ("hybrid f=0: SC-EDS 0.10 m gap as the suspension", ("eds", 0.0, 0.7, "skyhook"), None),
    ("hybrid f=0.75: EDS gap as the suspension", ("eds", 0.75, 0.7, "skyhook"), None),
    ("hybrid f=0.9: EDS gap as the suspension", ("eds", 0.9, 0.7, "skyhook"), None),
    ("active large gap 0.3 Hz, 0.3 m (any f)", (0.3, 0.7, "skyhook", 0.30), None),
    ("free flight f=1 (0.5 Hz, 0.12 m thrust gap)", (0.5, 0.7, "skyhook", 0.12), None),
]


def allow_amp(V, lam, case, F=BASE):
    """Largest sinusoidal guideway error (amplitude, m) at wavelength lam that keeps cabin <= 0.05 g and every
    stage within its travel (fixed guideway)."""
    name, prim, sec = case
    if prim[0] == "eds":
        _, f, z, kind = prim
        fm, sg = eds_freq(f, F), gap_usable(f, F)
    else:
        fm, z, kind, sg = prim
    k = 2 * math.pi / lam
    Ta_p, Tr_p = S.transfer(V, k, fm, z, kind)
    lims = [sg / max(math.sqrt(Tr_p), 1e-12)]
    Ta = Ta_p
    if sec:
        fs, zs, ss = sec
        Ta_s, Tr_s = S.transfer(V, k, fs, zs, "passive")
        lims.append(ss / max(math.sqrt(Ta_p * Tr_s), 1e-12))
        Ta = Ta_p * Ta_s
    lims.append(S.A_LIM / ((V * k) ** 2 * math.sqrt(Ta)))
    return min(lims)


def stage_tolerance(F=BASE):
    rows = []
    lams = np.geomspace(10.0, 1e5, 41)
    for M in (2.0, 4.0, 6.0, 8.0):
        V = M * C_REF
        for lam in lams:
            r = dict(Mach=M, lambda_m=lam, elevated_achievable_m=0.003 + 0.015 * lam / 1e3)
            for c in TOL_CASES:
                r[c[0]] = allow_amp(V, lam, c, F)
            rows.append(r)
    # elevated ground tube (regime: 3 mm + 15 mm/km): highest Mach each suspension tolerates on a fixed guideway
    for c in TOL_CASES:
        Mmax = float("nan")
        for M in np.arange(1.0, 12.01, 0.1):
            V = M * C_REF
            if all(allow_amp(V, lam, c, F) >= 0.003 + 0.015 * lam / 1e3 for lam in lams):
                Mmax = M
        rows.append(dict(kind="elevated ground tube, fixed guideway: highest Mach", case=c[0], Mach=Mmax))
        print(f"  tolerance {c[0]}: elevated M{Mmax:.1f}", flush=True)
    # floating tube motion per octave (stability.octaves), for the figure
    L = make_line(40.0, 0.1e9, F["eps"])
    for o in S.octaves(L, 2.0):
        rows.append(dict(kind="floating tube motion (0.1 GN, 40 m)", lambda_m=o["lambda_km"] * 1e3, amp_m=o["amp_m"]))
    write_csv("hybrid_tolerance.csv", rows)


# ============================================================================= 2. stroke table
_STROKE = {}


def stroke_table(T, eps, F=BASE, d_list=D_LIST):
    key = (round(T / 1e6, 3), round(eps, 4), F["f_s"], F["z_s"], tuple(d_list))
    if key in _STROKE:
        return _STROKE[key]
    tab, rows = {}, []
    for D in d_list:
        L = make_line(D, T, eps)
        cache = {}
        sy, sz, gy, gz = [], [], [], []
        for M in M_STROKE:
            ry = ride2(L, float(M), F["f_s"], F["z_s"], "passive", 1e9, "y", cache)
            rz = ride2(L, float(M), F["f_s"], F["z_s"], "passive", 1e9, "z", cache)
            sy.append(ry["stroke"])
            sz.append(rz["stroke"])
            gy.append(ry["gap"])
            gz.append(rz["gap"])
            rows.append(dict(T_GN=T / 1e9, eps=eps, D_bore=D, Mach=float(M), stroke_y_m=sy[-1], stroke_z_m=sz[-1],
                             gap_y_m=gy[-1], gap_z_m=gz[-1], fc_y_Hz=ry["fc"], env_R_m=L["Re"]))
        tab[D] = dict(y=np.array(sy), z=np.array(sz), gy=np.array(gy), gz=np.array(gz), L=L)
    _STROKE[key] = (tab, rows)
    return _STROKE[key]


def s_req(tab, D, M):
    return tuple(float(np.interp(M, M_STROKE, tab[D][k])) for k in ("y", "z", "gy", "gz"))


def stage_stroke(F=BASE):
    rows = []
    for T in (0.1e9, 0.3e9, 1.0e9):
        rows += stroke_table(T, F["eps"], F)[1]
    write_csv("hybrid_stroke.csv", rows)
    # consistency with flight.py (same line, same 1 Hz suspension)
    fl = {(r["T_GN"], r["D_bore"], r["Mach"]): r["stroke_y_m"] for r in FL.read_csv("flight_stroke.csv")
          if abs(r["eps"] - 0.1) < 1e-6}
    for r in rows:
        k = (r["T_GN"], r["D_bore"], r["Mach"])
        if k in fl and r["Mach"] in (2.0, 6.0, 10.0):
            print(f"    check vs flight: {k}: hybrid {r['stroke_y_m']:.2f} m, flight {fl[k]:.2f} m")


# ============================================================================= 3. one point
def hd(F, m):
    return FL.hd(F, m)


def body_state(M, p, d):
    """Body-only flight state (no wing): pod_state at a vanishing wing loading, body and nose terms only."""
    st = H.pod_state(M, p, 1.0, d)
    V, q = st["V"], st["q_Pa"]
    pd_ = R.pod(d["m"])
    fric = H.plate(p, H.T_GAS, V, pd_["L"], st["T_body"], d["re_tr"])[1] * pd_["S_wet"]
    D_b = fric + q * 0.15 * pd_["A"] + p * pd_["A"]
    rho, T0 = H.freestream(M, p)[1], st["T0"]
    A_nose = 0.5 * 2 * math.pi * d["R_n"] ** 2
    Q_b = (st["q_body_kWm2"] * pd_["S_wet"] + st["q_nose_kWm2"] * A_nose) * 1e3
    ex = 0.0
    if st["T_body"] > d["T_lim_skin"]:
        ex += H.q_net(lambda Tw: H.plate(p, H.T_GAS, V, pd_["L"], Tw, d["re_tr"])[0], d["T_lim_skin"]) * pd_["S_wet"]
    if st["T_nose"] > d["T_lim_le"]:
        ex += H.q_net(lambda Tw: H.q_stag(rho, V, T0, d["R_n"], Tw), d["T_lim_le"]) * A_nose
    return dict(V=V, q=q, D_body=D_b, Q_body=Q_b, Qex_body=ex, T_body=st["T_body"], T_nose=st["T_nose"],
                q_body=st["q_body_kWm2"], T0=T0, relief=st["relief"], Q_cabin_kW=None)


def hybrid_state(M, p, f, WS, b0, stw, d):
    """Combine the body with a wing carrying f of the weight at wing loading WS (per-area state from stw)."""
    pd_ = R.pod(d["m"])
    st = dict(Mach=M, V=b0["V"], p_Pa=p, q_Pa=b0["q"], relief=b0["relief"], T0=b0["T0"], T_body=b0["T_body"],
              T_nose=b0["T_nose"])
    if f <= 0:
        st.update(lift_ok=1, WS_Pa=0.0, S_m2=0.0, span_m=0.0, chord_m=pd_["d"], alpha_deg=0.0, CL=0.0,
                  drag_N=b0["D_body"], Q_in_kW=b0["Q_body"] / 1e3, Q_excess_kW=b0["Qex_body"] / 1e3,
                  T_le=0.0, T_wing_lower=b0["T_body"], q_wing_lower_kWm2=b0["q_body"])
    else:
        if not stw.get("lift_ok"):
            return dict(st, lift_ok=0, alpha_deg=90.0, WS_Pa=WS, S_m2=0.0, span_m=0.0)
        S_f = f * stw["S_m2"]
        b = math.sqrt(d["AR"] * S_f)
        D_wing_full = stw["drag_N"] - b0["D_body"]
        A_le_full = 0.5 * math.pi * d["R_le"] * stw["span_m"] / math.cos(math.radians(stw["sweep_deg"]))
        Q_wing_full = ((stw["q_wing_lower_kWm2"] + stw["q_wing_upper_kWm2"]) * stw["S_m2"]
                       + stw["q_le_kWm2"] * A_le_full) * 1e3
        Qex_wing_full = max(stw["Q_excess_kW"] * 1e3 - b0["Qex_body"], 0.0)
        st.update(lift_ok=1, WS_Pa=WS, S_m2=S_f, span_m=b, chord_m=S_f / b, alpha_deg=stw["alpha_deg"],
                  CL=stw["CL"], drag_N=b0["D_body"] + f * D_wing_full, Q_in_kW=(b0["Q_body"] + f * Q_wing_full) / 1e3,
                  Q_excess_kW=(b0["Qex_body"] + f * Qex_wing_full) / 1e3, T_le=stw["T_le"],
                  T_wing_lower=stw["T_wing_lower"], q_wing_lower_kWm2=stw["q_wing_lower_kWm2"])
    st["Q_cabin_kW"] = (d["k_ins"] / d["t_ins"] * max(0.0, min(st["T_body"], d["T_lim_skin"]) - d["T_cab"])
                        * pd_["S_wet"] + d["w_seat"] * pd_["n"]) / 1e3
    return st


def _st(mg, tol, k, m, t):
    mg[k], tol[k] = m, t


def evaluate(st, D, m, f, F, tab, d):
    pd_ = R.pod(m)
    n, pax = pd_["n"], pd_["n"] * LF
    r = dict(f=f, pod_t=m / 1e3, seats=n, D_bore=D, Mach=st["Mach"], p_Pa=st["p_Pa"], WS_Pa=st.get("WS_Pa", 0.0),
             span_m=st.get("span_m", 0.0))
    mg, tol = {}, {}
    if not st.get("lift_ok"):
        for k in CHECKS:
            _st(mg, tol, k, -9.0 if k == "lift" else 0.0, 0.3)
        r.update(alpha_deg=90.0)
        return finish(r, mg, tol)
    V, M, p = st["V"], st["Mach"], st["p_Pa"]
    Da = st["drag_N"]                                                # aerodynamic drag (heats the tube gas)
    W = m * G * st["relief"]
    Dm = (1 - f) * W / F["LD_mag"]                                   # magnetic drag (heats the guideway)
    Dt = Da + Dm
    A_b = math.pi * D * D / 4
    r.update(alpha_deg=st["alpha_deg"], CL=st["CL"], drag_aero_kN=Da / 1e3, drag_mag_kN=Dm / 1e3,
             LD_total=m * G / Dt, S_m2=st["S_m2"])
    # 1. lift: wing incidence (heating); f = 0 has no wing
    if f > 0:
        _st(mg, tol, "lift", math.log10(d["alpha_le_ok"] / max(st["alpha_deg"], 1e-3)),
            math.log10(d["alpha_le_marg"] / d["alpha_le_ok"]))
    else:
        _st(mg, tol, "lift", 1.0, 0.3)
    # 2. steering stroke: room beside the wing (or the pod + magnets) at 4 sigma, and below the body
    sy, sz, gy, gz = s_req(tab, D, M)
    w_lat = max(st["span_m"], pd_["d"] + F["guide_w"])
    av_y = (D - w_lat) / 2 - F["clear"]
    av_z = (D - (pd_["d"] + 1.5)) / 2 - F["clear"]
    _st(mg, tol, "stroke", min(math.log10(max(av_y, 1e-3) / sy), math.log10(max(av_z, 1e-3) / sz)), math.log10(1.5))
    r.update(stroke_req_y_m=sy, stroke_av_y_m=av_y, stroke_req_z_m=sz, D_req_m=w_lat + 2 * (sy + F["clear"]))
    # 3. magnetic gap travel: suspension relative motion + short-wave guideway error vs usable gap
    g_need = max(gy, gz) + F["rough"]
    g_av = gap_usable(f, F)
    _st(mg, tol, "gap", math.log10(g_av / g_need), math.log10(1.5))
    r.update(gap_need_m=g_need, gap_usable_m=g_av)
    # 4. trip length
    _st(mg, tol, "trip", math.log10(F["a_acc"] * F["L_trip"] / V ** 2), 1e-6)
    # capacity (sets traffic)
    h_rule = sm.headway_s(M, F["rule"])
    ts0 = H.tube_side(st, D, d, N=F["N_ok"])
    N_heat = ts0["N_max_per_h"]
    N_max = min(3600.0 / h_rule, N_heat)
    cap = N_max * n * LF * PEAK_FACTOR
    dem = F["demand"] * F["capture"]
    n_pairs = max(1.0, dem / cap)
    N_dem = dem / (n * LF * PEAK_FACTOR)
    N_ops = min(N_max, N_dem / n_pairs)
    r.update(N_heat_h=N_heat, N_max_h=N_max, cap_pair_Mpax=cap / 1e6, n_pairs=n_pairs, N_ops_h=N_ops)
    # 5. launch power (flight conventions; total drag)
    KE = 0.5 * m * V * V
    P_pod = (Dt + m * F["a_acc"]) * V / F["eta_d"]
    P_stn = KE / F["eta_d"] * N_ops / 3600.0
    z1 = math.log10(F["P_pod_ok"] / P_pod), math.log10(F["P_pod_marg"] / F["P_pod_ok"])
    z2 = math.log10(F["P_stn_ok"] / P_stn), math.log10(F["P_stn_marg"] / F["P_stn_ok"])
    _st(mg, tol, "power", *(z1 if z1[0] / z1[1] < z2[0] / z2[1] else z2))
    r.update(P_pod_MW=P_pod / 1e6, P_station_GW=P_stn / 1e9, KE_tTNT=KE / 4.184e9)
    # 6. motor: LSM frequency and thrust shear across the gap (the hybrid's gap, not flight's 0.15 m)
    tau = min(F["tau_max"], pd_["L"] / 6)
    fq = V / (2 * tau)
    sig_req = (Dt + m * F["a_acc"]) / (F["w_strip"] * 0.8 * pd_["L"])
    sig_av = F["sig0"] * math.exp(-math.pi * eds_gap(f, F) / tau)
    z1 = math.log10(F["f_ok"] / fq), math.log10(F["f_marg"] / F["f_ok"])
    z2 = math.log10(sig_av / sig_req), math.log10(2)
    _st(mg, tol, "motor", *(z1 if z1[0] / z1[1] < z2[0] / z2[1] else z2))
    r.update(lsm_Hz=fq, shear_req_kPa=sig_req / 1e3, shear_av_kPa=sig_av / 1e3)
    # 7. breach: drag deceleration in 20 km air (wing unloaded); wing overload if not, vs magnet hold-down
    q_amb = 0.5 * RHO_AMB * V * V
    a_br = q_amb * (0.5 * pd_["A"] + 0.006 * st["S_m2"]) / (m * G)
    _st(mg, tol, "breach", math.log10(F["g_breach_ok"] / a_br), math.log10(F["g_breach_marg"] / F["g_breach_ok"]))
    r.update(breach_decel_g=a_br, breach_wing_overload_x_weight=f * P_AMB / p)
    # 8. Kantrowitz blockage
    A_fr = pd_["A"] + st["span_m"] * 0.04 * st["chord_m"] + st["S_m2"] * math.sin(math.radians(st["alpha_deg"]))
    beta = A_fr / A_b
    bK = R.beta_kantrowitz(M)
    _st(mg, tol, "blockage", math.log10(0.5 * bK / beta), math.log10(2))
    r.update(beta=beta)
    # 9. bow shock reflected onto the pod
    d_w = max(0.5, (D - w_lat) / 2 - sy / 2)
    x_r = 2 * d_w * math.sqrt(M * M - 1)
    L_ref = max(pd_["L"], 2 * st["chord_m"])
    _st(mg, tol, "reflection", math.log10(x_r / L_ref), math.log10(2))
    # 10. self-pumping crossovers (aero drag only pushes gas)
    L_seg = p * A_b * V * h_rule / Da / 1e3
    _st(mg, tol, "crossover", math.log10(L_seg / 1.0), math.log10(4.0))
    r.update(crossover_km=L_seg)
    # 11. tube gas per pass and gas-to-wall throughput (aero drag only)
    _st(mg, tol, "gas_dT", math.log10(d["dT_pass_ok"] / ts0["dT_pass"]), math.log10(d["dT_pass_marg"] / d["dT_pass_ok"]))
    _st(mg, tol, "throughput", math.log10(N_heat / F["N_ok"]), math.log10(F["N_ok"] / F["N_marg"]))
    r.update(Pi=ts0["Pi"], dT_pass_K=ts0["dT_pass"])
    # 12. pod heating
    Tsk = max(st["T_wing_lower"], st["T_body"])
    _st(mg, tol, "skin_T", max(4 * math.log10(d["T_lim_skin"] / Tsk), -0.99 * math.log10(2)), math.log10(2))
    Tle = max(st["T_nose"], st["T_le"])
    _st(mg, tol, "le_T", max(4 * math.log10(d["T_lim_le"] / Tle), -0.99 * math.log10(2)), math.log10(2))
    fs, _, _ = H.sink_fraction(st, d)
    if math.isnan(fs):
        _st(mg, tol, "heat_sink", 0.0, 0.3)
    else:
        _st(mg, tol, "heat_sink", math.log10(d["f_sink_ok"] / max(fs, 1e-6)), math.log10(d["f_sink_marg"] / d["f_sink_ok"]))
    r.update(T_skin_K=Tsk, T_le_K=Tle, f_sink=fs)
    # 13. film and strip (strip also takes the magnetic drag as eddy heat)
    ts = H.tube_side(st, D, d, N=max(N_ops, 1e-3))
    E_conv = st["q_wing_lower_kWm2"] * 1e3 * 2 * st["chord_m"] / V
    q_strip = N_ops / 3600 * (ts["E_strip_pass_Jm2"] + (F["imp"] - 1) * E_conv + Dm / F["w_strip"])
    T_strip = (H.T_ENV ** 4 + q_strip / (2 * H.EPS * H.SIGMA)) ** 0.25
    Tf = max(ts["T_film_sun"], ts["T_film_tube"])
    _st(mg, tol, "film_T", 4 * math.log10(d["T_film_lim"] / Tf), 4 * math.log10(d["T_film_marg"] / d["T_film_lim"]))
    _st(mg, tol, "guideway_T", 4 * math.log10(d["T_strip_lim"] / T_strip), 4 * math.log10(453.0 / d["T_strip_lim"]))
    r.update(T_strip_K=T_strip)
    # 14. energy per pax-km
    L = F["L_trip"]
    e_ke = KE * (1 / F["eta_d"] - F["eta_d"] * F["eta_r"]) / (pax * L) / 3.6
    e_dr = Da / F["eta_d"] / pax / 3.6
    e_mg = Dm / F["eta_d"] / pax / 3.6
    e_ht = F["hotel_w_seat"] * n / pax / V / 3.6
    leak = F["leak"] * D / 3.0
    P_pump = 10e3 * leak * R.q_hole(P_AMB, T_AMB, 250.0) / p
    paxh = max(N_ops * pax, 1e-6)
    e_pu = P_pump / paxh
    e_ac = F["act_kw_km"] * 1e3 / paxh
    E = e_ke + e_dr + e_mg + e_ht + e_pu + e_ac
    _st(mg, tol, "energy", math.log10(F["E_ok"] / E), math.log10(F["E_marg"] / F["E_ok"]))
    r.update(E_Wh_paxkm=E, E_kinetic=e_ke, E_aero=e_dr, E_mag=e_mg, E_drag=e_dr + e_mg, E_pump=e_pu,
             E_other=e_ht + e_ac, d2d_NYC_LAX_h=FL.d2d_h(M, F["a_acc"]))
    return finish(r, mg, tol)


def finish(r, mg, tol):
    def grade(skip):
        s, w, z, nf = "pass", None, 1e9, 0
        for k in CHECKS:
            if k in skip:
                continue
            g = "pass" if mg[k] >= 0 else ("marginal" if mg[k] >= -tol[k] else "fail")
            if g == "fail":
                s, nf = "fail", nf + 1
            elif g == "marginal" and s == "pass":
                s = "marginal"
            if mg[k] / tol[k] < z:
                z, w = mg[k] / tol[k], k
        return s, w, z, nf
    for k in CHECKS:
        r[f"m_{k}"] = round(mg[k], 4)
    s, w, z, nf = grade(())
    snt, wnt, znt, _ = grade(("trip",))
    r.update(status=s, binding=w, score=z, nfail=nf, status_nt=snt, binding_nt=wnt, score_nt=znt)
    return r


def better(a, b, key="status"):
    return FL.better(a, b, key)


KEEP = ["f", "pod_t", "seats", "D_bore", "Mach", "p_Pa", "WS_Pa", "span_m", "S_m2", "alpha_deg", "CL",
        "drag_aero_kN", "drag_mag_kN", "LD_total", "stroke_req_y_m", "stroke_av_y_m", "stroke_req_z_m", "D_req_m",
        "gap_need_m", "gap_usable_m", "N_heat_h", "N_max_h", "cap_pair_Mpax", "n_pairs", "N_ops_h", "P_pod_MW",
        "P_station_GW", "KE_tTNT", "lsm_Hz", "shear_req_kPa", "shear_av_kPa", "breach_decel_g",
        "breach_wing_overload_x_weight", "beta", "crossover_km", "Pi", "dT_pass_K", "T_skin_K", "T_le_K", "f_sink",
        "T_strip_K", "E_Wh_paxkm", "E_kinetic", "E_aero", "E_mag", "E_drag", "E_pump", "E_other", "d2d_NYC_LAX_h",
        "status", "binding", "score", "nfail", "status_nt", "binding_nt", "score_nt"]


def run_grid(F, f_list=F_LIST, pods=PODS, d_list=D_LIST, m_grid=M_GRID, p_grid=P_GRID, ws_list=WS_LIST,
             key="status", out=None, done=()):
    tab, _ = stroke_table(F["T_line"], F["eps"], F, tuple(d_list))
    rows = []
    for m in pods:
        d = hd(F, m)
        todo = [f for f in f_list if (f, m) not in done]
        if not todo:
            continue
        prow = {f: [] for f in todo}
        for M in m_grid:
            for p in p_grid:
                b0 = body_state(float(M), float(p), d)
                stws = {WS: H.pod_state(float(M), float(p), WS, d) for WS in ws_list} if any(f > 0 for f in todo) else {}
                for f in todo:
                    best = {D: None for D in d_list}
                    cands = [(0.0, None)] if f == 0 else list(stws.items())
                    for WS, stw in cands:
                        st = hybrid_state(float(M), float(p), f, WS, b0, stw, d)
                        for D in d_list:
                            r = evaluate(st, D, m, f, F, tab, d)
                            if better(r, best[D], key):
                                best[D] = r
                    prow[f] += [{k: best[D].get(k) for k in KEEP} for D in d_list]
        for f in todo:
            rows += prow[f]
            if out:
                write_csv(out, prow[f], mode="a")
            n_ok = sum(r["status"] == "pass" for r in prow[f])
            print(f"  pod {m/1e3:g} t, f {f:g}: {n_ok} pass of {len(prow[f])}", flush=True)
    return rows


# ============================================================================= stages
def stage_grid(F=BASE):
    path = os.path.join(HERE, "hybrid_grid.csv")
    done = set()
    if os.path.exists(path) and "--fresh" not in sys.argv:
        old = read_csv("hybrid_grid.csv")
        cnt = {}
        for r in old:
            cnt[(r["f"], r["pod_t"] * 1e3)] = cnt.get((r["f"], r["pod_t"] * 1e3), 0) + 1
        full = len(M_GRID) * len(P_GRID) * len(D_LIST)
        done = {k for k, v in cnt.items() if v >= full}
        if len(done) < len(cnt):                                     # drop a partial block, rewrite the rest
            keep = [r for r in old if (r["f"], r["pod_t"] * 1e3) in done]
            os.remove(path)
            if keep:
                write_csv("hybrid_grid.csv", keep)
    elif os.path.exists(path):
        os.remove(path)
    run_grid(F, out="hybrid_grid.csv", done=done)
    return read_csv("hybrid_grid.csv")


def max_mach(rows, key="status", ok=("pass",)):
    Ms = [r["Mach"] for r in rows if r[key] in ok]
    return max(Ms) if Ms else float("nan")


def top_binding(rows, Mtop, key="status_nt"):
    bk = "binding" if key == "status" else "binding_nt"
    sk = "score" if key == "status" else "score_nt"
    if math.isnan(Mtop):
        cand = rows
    else:
        nxt = [r["Mach"] for r in rows if r["Mach"] > Mtop + 1e-9]
        if not nxt:
            return "grid edge"
        cand = [r for r in rows if r["Mach"] == min(nxt)]
    best = max(cand, key=lambda r: (RANK[r[key]], -r["nfail"], r[sk]))
    return best[bk]


def stage_band(rows):
    out = []
    for f in F_LIST:
        for m in sorted(set(r["pod_t"] for r in rows)):
            for D in D_LIST:
                sub = [r for r in rows if r["f"] == f and r["pod_t"] == m and r["D_bore"] == D]
                if not sub:
                    continue
                Mp = max_mach(sub, "status_nt")
                Mm = max_mach(sub, "status_nt", ("pass", "marginal"))
                row = dict(f=f, pod_t=m, seats=sub[0]["seats"], D_bore=D, Mmax_pass=Mp, Mmax_marg=Mm,
                           bound_above=top_binding(sub, Mp))
                ps = [r for r in sub if r["status_nt"] == "pass"]
                if ps:
                    top = [r for r in ps if r["Mach"] == Mp]
                    row.update(p_lo_at_Mmax=min(r["p_Pa"] for r in top), p_hi_at_Mmax=max(r["p_Pa"] for r in top),
                               E_at_Mmax=min(r["E_Wh_paxkm"] for r in top), E_min=min(r["E_Wh_paxkm"] for r in ps),
                               cap_max_Mpax=max(r["cap_pair_Mpax"] for r in ps))
                out.append(row)
    write_csv("hybrid_band.csv", out)


def stage_bore(rows):
    """Smallest bore (of D_LIST) with a passing (and a marginal) point at each Mach, per f x pod."""
    out = []
    for f in F_LIST:
        for m in sorted(set(r["pod_t"] for r in rows)):
            for M in M_GRID:
                sub = [r for r in rows if r["f"] == f and r["pod_t"] == m and r["Mach"] == M]
                if not sub:
                    continue
                ps = [r for r in sub if r["status_nt"] == "pass"]
                mk = [r for r in sub if r["status_nt"] in ("pass", "marginal")]
                row = dict(f=f, pod_t=m, Mach=float(M), D_min_pass=min((r["D_bore"] for r in ps), default=float("nan")),
                           D_min_marg=min((r["D_bore"] for r in mk), default=float("nan")))
                if ps:
                    Dm = row["D_min_pass"]
                    b = min([r for r in ps if r["D_bore"] == Dm], key=lambda r: r["E_Wh_paxkm"])
                    row.update(p_Pa=b["p_Pa"], E_Wh_paxkm=b["E_Wh_paxkm"], span_m=b["span_m"],
                               stroke_req_y_m=b["stroke_req_y_m"], D_req_m=b["D_req_m"])
                # geometric bore need (stroke + width), smallest over pressure at 80 m (other checks pass)
                g80 = [r for r in sub if r["D_bore"] == 80.0 and r["status_nt"] != "fail"]
                row["D_geom_m"] = min((r["D_req_m"] for r in g80), default=float("nan"))
                out.append(row)
    write_csv("hybrid_bore.csv", out)


PICK = ["f", "pod_t", "seats", "D_bore", "Mach", "p_Pa", "WS_Pa", "span_m", "LD_total", "drag_aero_kN",
        "drag_mag_kN", "stroke_req_y_m", "stroke_av_y_m", "gap_need_m", "gap_usable_m", "T_le_K", "T_skin_K",
        "dT_pass_K", "N_max_h", "cap_pair_Mpax", "n_pairs", "E_Wh_paxkm", "E_aero", "E_mag", "E_kinetic", "E_pump",
        "P_pod_MW", "P_station_GW", "KE_tTNT", "breach_decel_g", "breach_wing_overload_x_weight", "lsm_Hz",
        "d2d_NYC_LAX_h", "status_nt", "binding_nt"]


def pick(r):
    return {k: r.get(k) for k in PICK}


def stage_points(rows):
    out = []
    ok = [r for r in rows if r["status_nt"] == "pass"]
    okm = [r for r in rows if r["status_nt"] in ("pass", "marginal")]
    for f in F_LIST:
        for D in (20.0, 40.0, 80.0, None):
            c = [r for r in ok if r["f"] == f and (D is None or r["D_bore"] == D)]
            cm = [r for r in okm if r["f"] == f and (D is None or r["D_bore"] == D)]
            allr = [r for r in rows if r["f"] == f and (D is None or r["D_bore"] == D)]
            tag = "any bore" if D is None else f"{D:g} m"
            if cm:
                Mt = max(r["Mach"] for r in cm)
                b = min([r for r in cm if r["Mach"] == Mt], key=lambda r: (r["D_bore"], r["E_Wh_paxkm"]))
                out.append(dict(kind=f"fastest marginal, {tag}", bound_above=top_binding(allr, Mt), **pick(b)))
            if not c:
                out.append(dict(kind=f"fastest passing, {tag}", f=f, Mach=float("nan"),
                                bound_above=top_binding(allr, float("nan"))))
                continue
            Mt = max(r["Mach"] for r in c)
            b = min([r for r in c if r["Mach"] == Mt], key=lambda r: (r["D_bore"], r["E_Wh_paxkm"]))
            out.append(dict(kind=f"fastest passing, {tag}", bound_above=top_binding(allr, Mt), **pick(b)))
            e = min(c, key=lambda r: r["E_Wh_paxkm"])
            out.append(dict(kind=f"energy-optimal, {tag}", **pick(e)))
            k = max(c, key=lambda r: (r["cap_pair_Mpax"], -r["E_Wh_paxkm"]))
            out.append(dict(kind=f"capacity-optimal, {tag}", **pick(k)))
    write_csv("hybrid_points.csv", out)


# ----------------------------------------------------------------------------- economics (flight/system conventions)
_ADD = {}


def adders(D, T, eps):
    """$M/km added by bore wall, tension, envelope/tether growth and guideway actuators (flight.adders)."""
    key = (D, round(T / 1e6, 3), round(eps, 4))
    if key not in _ADD:
        L = make_line(D, T, eps)
        add_kg = {"tension": S.tension_mass(T), "wall": max(0.0, S.tube_props(D)[0] - S.tube_props(S.D_BASE)[0])}
        _, _, cb, cl = S.economics(4.0, L, add_kg, S.ACT_USD_M_KM, 0.0, 0.0, D)
        f4 = R.float_cost(4.0)
        _ADD[key] = (cb - f4, cl - f4 * 0.28)
    return _ADD[key]


def tube_cost(net, M, D, F):
    """As flight.tube_cost, but the run keeps its levitation guideway AND gets a full-length LSM."""
    km = FL.NETS[net]
    if net == "NYC-LAX direct":
        ramp = 2 * min(km / 2, max(sm.accel_km(M), sm.P["ramp_min_km"]))
    else:
        sd = FL.sysd()
        fr = np.interp(min(M, 6.0), [2, 3, 4, 6], [sd[(net, "sub", 1.0, x)]["ramp"] / km for x in (2.0, 3.0, 4.0, 6.0)])
        ramp = km * (1.0 if M >= 6 else fr)
    run = km - ramp
    a_b, a_l = adders(D, F["T_line"], F["eps"])
    out = {}
    for lvl, add in (("base", a_b), ("low", a_l)):
        i = sm.LVL[lvl]
        c_run = sm.run_cost(lvl)
        c_ramp = c_run * sm.RAMP_MULT[i] + sm.LSM[i] * (1 + sm.SOFT[i]) * (1 + sm.CONTING[i])
        out[lvl] = (run * (c_run + FL.LSM_EFF[lvl]) + ramp * c_ramp) / km + add
    return out


def bc_hook(r, net, F):
    """flight.bc_hook with the hybrid tube cost (2050 opening; gross/net x 30/100 % x comparator x cost level)."""
    sd = FL.sysd()
    M = r["Mach"]
    Mx = min(M, 6.0)
    tc = tube_cost(net, M, r["D_bore"], F)
    cap_pair = r["cap_pair_Mpax"] * 1e6
    out = dict(net=net, tube_base_usd_m_km=tc["base"], tube_low_usd_m_km=tc["low"])
    k_run = sm.run_cost("low") / sm.run_cost("base")
    for vs in FL.CMP:
        for cap in (0.3, 1.0):
            v = {k: float(np.interp(Mx, [2, 3, 4, 6], [sd[(net, vs, cap, x)][k] for x in (2.0, 3.0, 4.0, 6.0)]))
                 for k in ("ben", "en", "C0", "dC", "C0l", "c_base", "paxkm", "wh_sys")}
            ben = v["ben"]
            e_sys = v["wh_sys"]
            if net == "NYC-LAX direct" and M > 6:
                sv6, svM = FL.CMP[vs] - FL.d2d_h(6.0), FL.CMP[vs] - FL.d2d_h(M)
                ben *= svM / sv6
                e_sys = float(sm.kwh_per_pax(4045, M)[0]) * 1e3 / 4045
            e_mag = sm.P["kg_per_seat"] / LF * G / sm.P["LD_mag"] / sm.P["eta_drive"] / 3.6
            dE = r["E_drag"] + r["E_pump"] + r["E_other"] - e_mag
            if net == "NYC-LAX direct" and M > 6:
                dE += e_sys - v["wh_sys"]
            dPV = dE * v["paxkm"] * cap * R.USD_KWH * R.PV / 1e3
            pairs = max(1.0, F["demand"] * cap / cap_pair)
            for lvl in ("base", "low"):
                C0 = v["C0"] if lvl == "base" else v["C0l"]
                dC = v["dC"] * (1.0 if lvl == "base" else k_run)
                cb = v["c_base"] * (1.0 if lvl == "base" else k_run)
                C = C0 + dC * tc[lvl] * pairs / cb
                tag = f"{vs}_{int(cap * 100)}_{lvl}"
                out[f"bcg_{tag}"] = ben / C
                out[f"bcn_{tag}"] = (ben + v["en"] - dPV) / C
            out[f"pairs_{int(cap * 100)}"] = pairs
    return out


def stage_econ(rows, F=BASE):
    out = []
    cand = [r for r in rows if r["status_nt"] in ("pass", "marginal")]
    for net in FL.NETS:
        for f in F_LIST:
            for D in (20.0, 40.0, 80.0):
                for M in sorted(set(r["Mach"] for r in cand)):
                    if abs(M * 2 - round(M * 2)) > 1e-9:
                        continue
                    c = [r for r in cand if r["f"] == f and r["D_bore"] == D and r["Mach"] == M]
                    if not c:
                        continue
                    best = None
                    for r in c:
                        h = bc_hook(r, net, F)
                        if best is None or h["bcn_sub_100_base"] > best[0]:           # flight.stage_econ convention
                            best = (h["bcn_sub_100_base"], h, r)
                    _, bh, br = best
                    row = dict(f=f, Mach=M, D_bore=D, pod_t=br["pod_t"], p_Pa=br["p_Pa"], status=br["status_nt"],
                               E_Wh_paxkm=br["E_Wh_paxkm"], cap_pair_Mpax=br["cap_pair_Mpax"],
                               d2d_NYC_LAX_h=br["d2d_NYC_LAX_h"])
                    row.update(bh)
                    out.append(row)
    write_csv("hybrid_econ.csv", out)


LEVERS = [
    ("base (0.1 GN line, LD_mag 300, moving block 10 s, 0.3 g)", {}),
    ("line tension 0.3 GN", dict(T_line=0.3e9)),
    ("line tension 1 GN", dict(T_line=1.0e9)),
    ("magnetic L/D 100 (vs 300)", dict(LD_mag=100.0)),
    ("magnetic L/D 1,000", dict(LD_mag=1000.0)),
    ("low-AR wing (AR 1)", dict(AR=1.0)),
    ("launch power cap 100 MW per pod (vs 200)", dict(P_pod_ok=100e6, P_pod_marg=400e6)),
    ("comfort accel 0.15 g", dict(a_acc=0.15 * G)),
    ("brick-wall headway (0.5 g)", dict(rule="brick-wall")),
    ("30 % capture (24 M pax/yr)", dict(capture=0.3)),
    ("energy pass at 100 Wh/pax-km (vs 340)", dict(E_ok=100.0, E_marg=200.0)),
    ("station launch budget 5 GW pass (vs 1)", dict(P_stn_ok=5e9, P_stn_marg=20e9)),
    ("station 5 GW + line tension 0.3 GN", dict(P_stn_ok=5e9, P_stn_marg=20e9, T_line=0.3e9)),
    ("station 5 GW + line tension 1 GN", dict(P_stn_ok=5e9, P_stn_marg=20e9, T_line=1.0e9)),
]


def stage_levers():
    path = os.path.join(HERE, "hybrid_levers.csv")
    have = {r["lever"] for r in read_csv("hybrid_levers.csv")} if os.path.exists(path) and "--fresh" not in sys.argv else set()
    if not have and os.path.exists(path):
        os.remove(path)
    mg = np.round(np.arange(2.0, 10.001, 0.5), 2)
    pg = np.geomspace(5.0, 1000.0, 10)
    dl = [20.0, 40.0, 80.0]
    for name, mod in LEVERS:
        if name in have:
            continue
        F = dict(BASE, **mod)
        g = run_grid(F, f_list=[0.0, 0.5, 1.0], pods=[17000.0, 47000.0], d_list=dl, m_grid=mg, p_grid=pg, key="status_nt")
        r = dict(lever=name)
        for f in (0.0, 0.5, 1.0):
            for D in dl:
                sub = [x for x in g if x["D_bore"] == D and x["f"] == f]
                Mp = max_mach(sub, "status_nt")
                r[f"f{f:g}_D{D:g}_Mmax_pass"] = Mp
                r[f"f{f:g}_D{D:g}_Mmax_marg"] = max_mach(sub, "status_nt", ("pass", "marginal"))
                r[f"f{f:g}_D{D:g}_bound"] = top_binding(sub, Mp)
                ps = [x for x in sub if x["status_nt"] == "pass"]
                r[f"f{f:g}_D{D:g}_minE"] = min((x["E_Wh_paxkm"] for x in ps), default=float("nan"))
        write_csv("hybrid_levers.csv", [r], mode="a")
        print(f"  lever {name}: " + ", ".join(f"f{f:g}/D40 M{r[f'f{f:g}_D40_Mmax_pass']}" for f in (0.0, 0.5, 1.0)),
              flush=True)


# ----------------------------------------------------------------------------- 5. onboard propulsion side check
def stage_propulsion():
    rows = []
    m, M, pax, L = 17000.0, 5.0, 40.0, 4000e3
    V = M * C_REF
    t = L / V
    d = dict(H.BASE, m=m)
    for f, p, WS in ((1.0, 200.0, 400.0), (0.0, 50.0, 400.0)):
        if f > 0:
            Dr = H.pod_state(M, p, WS, d)["drag_N"]
        else:
            Dr = body_state(M, p, d)["D_body"] + m * G / BASE["LD_mag"]
        q = 0.7 * p * M * M
        I = Dr * t
        A_b = math.pi * 40.0 ** 2 / 4
        rho = p / (287.0 * 216.65)
        for name, isp, kind in (("rocket LOX/LH2 (Isp 450 s)", 450.0, "rocket"),
                                ("rocket kerosene/LOX (Isp 330 s)", 330.0, "rocket"),
                                ("scramjet / ramjet (needs q ~25-100 kPa)", 1000.0, "air"),
                                ("battery + fan or propeller", None, "fan"),
                                ("battery + onboard LIM on a passive reaction rail", None, "lim")):
            r = dict(case=name, f=f, Mach=M, p_Pa=p, drag_kN=Dr / 1e3, cruise_s=t, impulse_MNs=I / 1e6, q_tube_Pa=q)
            if kind == "rocket":
                mp_lin = I / (isp * G)
                mp_eq = m * (math.exp(Dr / m * t / (isp * G)) - 1) if f == 0 else mp_lin   # f=1: drag ~ weight
                r.update(propellant_t=mp_lin / 1e3, propellant_rocket_eq_t=mp_eq / 1e3, propellant_over_pod=mp_lin / m,
                         exhaust_kg_per_km_tube=mp_lin / (L / 1e3), tube_gas_kg_per_km_40m=rho * A_b * 1e3)
                # exhaust (H2O / CO2) to pump at 100 pods/h, as vapour at tube pressure and 250 K
                vdot = mp_lin / (L / 1e3) * 100 / 3600 / (p / (461.5 * 250.0))
                r.update(pump_m3s_per_km_100pods_h=vdot, pump_MW_per_km=vdot * 10e3 / 1e6,
                         pump_vs_regime_installed_x=vdot / R.S_LIM)
            elif kind == "air":
                r.update(q_needed_kPa="25-100", q_shortfall_x=50e3 / q, combustor_p_needed_kPa=">=50",
                         combustor_p_available_kPa=p * 30 / 1e3)
            elif kind == "fan":
                mdot = rho * V * R.pod(m)["A"]                       # all air through the pod's frontal area
                ve = V + Dr / mdot
                r.update(captured_kg_s=mdot, exhaust_velocity_ms=ve, prop_eff=2 * V / (V + ve),
                         power_MW=Dr * V / (2 * V / (V + ve)) / 1e6,
                         battery_t_at_0p3kWh_kg=Dr * L / (2 * V / (V + ve)) / 3.6e6 / 0.3 / 1e3,
                         inlet_T0_K=216.65 * (1 + 0.2 * M * M))
            else:
                r.update(battery_t_at_0p3kWh_kg=Dr * L / 0.85 / 3.6e6 / 0.3 / 1e3, energy_MWh=Dr * L / 0.85 / 3.6e9)
            rows.append(r)
    write_csv("hybrid_propulsion.csv", rows)


# ============================================================================= figures
INK, INK2, GRIDC, SURF = FL.INK, FL.INK2, FL.GRIDC, FL.SURF
PASS, MARG, FAIL = FL.PASS, FL.MARG, FL.FAIL
S1, S2, S3, S4 = FL.S1, FL.S2, FL.S3, FL.S4
F_COL = {0.0: "#6da7ec", 0.25: "#3987e5", 0.5: "#256abf", 0.75: "#184f95", 0.9: "#104281", 1.0: "#0d366b"}


def fig_maps(pod=17.0):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.colors import ListedColormap
    from matplotlib.lines import Line2D
    from matplotlib.patches import Patch
    rows = [r for r in read_csv("hybrid_grid.csv") if r["pod_t"] == pod]
    dsh = [20.0, 40.0, 80.0]
    fig, axs = plt.subplots(len(dsh), len(F_LIST), figsize=(2.75 * len(F_LIST), 2.6 * len(dsh) + 0.9),
                            sharex=True, sharey=True)
    fig.patch.set_facecolor(SURF)
    cmap = ListedColormap([FAIL, MARG, PASS])
    for i, D in enumerate(dsh):
        for j, f in enumerate(F_LIST):
            ax = axs[i, j]
            sub = [r for r in rows if r["D_bore"] == D and r["f"] == f]
            FL._style(ax, f"f = {f:g} · bore {D:g} m")
            if not sub:
                continue
            xs = np.array(sorted(set(r["Mach"] for r in sub)))
            ys = np.array(sorted(set(r["p_Pa"] for r in sub)))
            ix = {x: k for k, x in enumerate(xs)}
            iy = {y: k for k, y in enumerate(ys)}
            Z = np.zeros((len(ys), len(xs)))
            B = np.empty((len(ys), len(xs)), object)
            E = np.full((len(ys), len(xs)), np.nan)
            for r in sub:
                a, b_ = iy[r["p_Pa"]], ix[r["Mach"]]
                Z[a, b_] = RANK[r["status_nt"]]
                B[a, b_] = r["binding_nt"]
                if r["status_nt"] != "fail" and isinstance(r.get("E_Wh_paxkm"), float):
                    E[a, b_] = r["E_Wh_paxkm"]
            ax.pcolormesh(FL._edges(xs), FL._edges(ys, True), Z, cmap=cmap, vmin=-0.5, vmax=2.5, shading="flat")
            codes = {k: n for n, k in enumerate(sorted(set(B.ravel()) - {None, ""}))}
            C = np.vectorize(lambda b: codes.get(b, -1))(B).astype(float)
            C[Z == 2] = np.nan
            for bnd, n_ in codes.items():
                mask = C == n_
                if mask.sum() < 0.05 * mask.size:
                    continue
                jj, ii = np.nonzero(mask)
                k = np.argmin((jj - np.median(jj)) ** 2 + (ii - np.median(ii)) ** 2)
                ax.text(xs[ii[k]], ys[jj[k]], SHORT.get(bnd, bnd), fontsize=6, color=INK, ha="center", va="center",
                        bbox=dict(boxstyle="round,pad=0.12", fc=SURF, ec="none", alpha=0.85))
            with np.errstate(invalid="ignore"):
                if np.isfinite(E).any() and np.nanmin(E) < 150:
                    cs = ax.contour(xs, ys, E, levels=[100, 150], colors=S2, linewidths=0.8, linestyles=["-", "--"])
                    ax.clabel(cs, fmt=lambda v: f"{v:.0f} Wh", fontsize=5.5)
            ax.set_yscale("log")
            if i == len(dsh) - 1:
                ax.set_xlabel("Mach", fontsize=8)
            if j == 0:
                ax.set_ylabel("tube pressure, Pa", fontsize=8)
    handles = [Patch(fc=PASS, label="passes every check"), Patch(fc=MARG, label="marginal"),
               Patch(fc=FAIL, label="fails (label = binding check)"),
               Line2D([], [], color=S2, lw=0.8, label="100 / 150 Wh per pax-km (jet fuel ~340)")]
    fig.legend(handles=handles, loc="lower center", ncol=4, fontsize=7.5, frameon=False)
    fig.suptitle(f"Maglev-plus-wings hybrid, {pod:g} t / 50-seat pod: feasible Mach × pressure by wing-lift fraction f "
                 "and bore (0.1 GN line, moving block 10 s, 4,000 km)", fontsize=10, color=INK, x=0.01, ha="left")
    fig.tight_layout(rect=(0, 0.05, 1, 0.97))
    fig.savefig(os.path.join(HERE, "hybrid_maps.png"), dpi=125, facecolor=SURF)
    plt.close(fig)
    print("wrote hybrid_maps.png")


def fig_summary():
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    tol = read_csv("hybrid_tolerance.csv")
    bore = read_csv("hybrid_bore.csv")
    band = read_csv("hybrid_band.csv")
    ec = read_csv("hybrid_econ.csv")
    grid = [r for r in read_csv("hybrid_grid.csv") if r["D_bore"] == 80.0]
    fig, axs = plt.subplots(1, 4, figsize=(17, 4.1))
    fig.patch.set_facecolor(SURF)
    # (a) allowable error amplitude vs wavelength at M6, fixed guideway, + floating tube motion
    ax = axs[0]
    FL._style(ax, "a. Allowed guideway error, fixed guideway, Mach 6")
    t6 = [r for r in tol if r.get("Mach") == 6.0 and r.get("kind", "") in ("", None)]
    lam = [r["lambda_m"] for r in t6]
    cols = [S2, S4, S3, "#008300", S1, INK2]
    for c, col in zip(TOL_CASES, cols):
        ax.plot(lam, [r[c[0]] for r in t6], color=col, lw=1.6, label=c[0].split(":")[0].split("(")[0].strip())
    ax.plot(lam, [r["elevated_achievable_m"] for r in t6], color=INK, lw=1.0, ls=":", label="elevated ground: 3 mm + 15 mm/km")
    fm = [r for r in tol if r.get("kind") == "floating tube motion (0.1 GN, 40 m)"]
    ax.plot([r["lambda_m"] for r in fm], [r["amp_m"] for r in fm], color=INK, lw=1.4, ls="--", label="floating tube motion")
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_ylim(1e-4, 1e3)
    ax.set_xlabel("wavelength, m")
    ax.set_ylabel("amplitude, m")
    ax.legend(fontsize=5.8, frameon=False, loc="upper left")
    # (b) bore needed for stroke + span (or pod + magnets), 17 t, other checks not failing
    ax = axs[1]
    FL._style(ax, "b. Bore needed for stroke + width, 17 t, 0.1 GN")
    for f in F_LIST:
        s = sorted([r for r in bore if r["f"] == f and r["pod_t"] == 17.0 and r["D_geom_m"] == r["D_geom_m"]],
                   key=lambda r: r["Mach"])
        ax.plot([r["Mach"] for r in s], [r["D_geom_m"] for r in s], color=F_COL[f], lw=1.8, label=f"f = {f:g}")
    for y in (20, 40):
        ax.axhline(y, color=INK2, lw=0.8, ls=":")
    ax.set_ylim(0, 80)
    ax.set_xlabel("Mach")
    ax.set_ylabel("bore diameter, m")
    ax.legend(fontsize=7, frameon=False, ncol=2, loc="lower right")
    # (c) lowest energy among non-failing points, 80 m bore, any pod
    ax = axs[2]
    FL._style(ax, "c. Lowest energy, any pod, 80 m bore (not failing)")
    for f in F_LIST:
        xs, ys = [], []
        for M in M_GRID:
            c = [r["E_Wh_paxkm"] for r in grid if r["f"] == f and r["Mach"] == M and r["D_bore"] == 80.0
                 and r["status_nt"] != "fail"]
            if c:
                xs.append(M)
                ys.append(min(c))
        ax.plot(xs, ys, color=F_COL[f], lw=1.6, label=f"f = {f:g}")
    ax.axhline(340, color=INK2, ls="--", lw=0.9)
    ax.text(2.1, 347, "jet fuel ~340", fontsize=6.5, color=INK2, va="bottom")
    ax.set_ylim(0, 480)
    ax.set_xlabel("Mach")
    ax.set_ylabel("Wh per pax-km")
    ax.legend(fontsize=7, frameon=False, ncol=3, loc="upper left")
    # (d) B/C net vs subsonic, 100 %, low cost, NYC-LAX, 40 m
    ax = axs[3]
    FL._style(ax, "d. Best B/C net vs subsonic, 100 %, low cost, 40 m")
    for f in F_LIST:
        s = sorted([r for r in ec if r["f"] == f and r["D_bore"] == 40.0 and r["net"] == "NYC-LAX direct"],
                   key=lambda r: r["Mach"])
        ax.plot([r["Mach"] for r in s], [r["bcn_sub_100_low"] for r in s], color=F_COL[f], lw=1.6, label=f"f = {f:g}")
    ax.axhline(1.0, color=INK, lw=0.8)
    ax.set_xlabel("Mach")
    ax.set_ylabel("B/C")
    ax.legend(fontsize=7, frameon=False, ncol=2)
    fig.tight_layout()
    fig.savefig(os.path.join(HERE, "hybrid_summary.png"), dpi=125, facecolor=SURF)
    plt.close(fig)
    print("wrote hybrid_summary.png")


# ============================================================================= main
def main():
    args = [a for a in sys.argv[1:] if not a.startswith("--")]
    stages = args or ["susp", "tolerance", "stroke", "propulsion", "grid", "band", "bore", "points", "econ", "levers",
                      "figs"]
    rows = None
    if "susp" in stages:
        stage_suspension()
    if "tolerance" in stages:
        stage_tolerance()
    if "stroke" in stages:
        stage_stroke()
    if "propulsion" in stages:
        stage_propulsion()
    if "grid" in stages:
        rows = stage_grid()
    if any(s in stages for s in ("band", "bore", "points", "econ")) and rows is None:
        rows = read_csv("hybrid_grid.csv")
    if "band" in stages:
        stage_band(rows)
    if "bore" in stages:
        stage_bore(rows)
    if "points" in stages:
        stage_points(rows)
    if "econ" in stages:
        stage_econ(rows)
    if "levers" in stages:
        stage_levers()
    if "figs" in stages:
        fig_maps()
        fig_summary()


if __name__ == "__main__":
    main()
