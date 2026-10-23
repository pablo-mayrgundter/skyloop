#!/usr/bin/env python3
"""Floating-tube stability study (Skyloop review, screening level).

Maps the rough feasible regime of a floating, tethered, H2-lifted tube at 20 km against the
ride/alignment limit, across stabilisation schemes (tethers, bracing, keel truss, taut line,
active trim, big-bore decoupling) and Mach 1-4.

Run:   python3 skyloop-review/stability/stability.py        (~1-3 min, one process)
       python3 skyloop-review/stability/stability.py --figs (figures only, from CSVs)
Out:   gust_spectrum.csv      wind/load spectrum per octave of wavelength, steady loads
       tether_field.csv       horizontal/vertical stiffness per km vs layout, fibre mass
       motion_octaves.csv     tube displacement per octave vs allowable, per scheme and Mach
       keel_sizing.csv        keel EI and mass vs depth, span, Mach
       active_trim.csv        actuator options: power, bandwidth, stroke
       bigbore.csv            D = 5/8/20/40/80 m: wall, self-lift, sail, envelope, pod flight, pumping, breach
       feasibility_map.csv    Mach x scheme (x bore): margins, status, binding constraint, B/C
       stability_spectra.png, stability_map.png, stability_bigbore.png

Model (all (A) unless cited): lateral line model per metre of corridor
   K(k) = k_t g(k s) + T k^2 + EI k^4,   H(k, w) = 1 / (K - m w^2 + i c w)
driven by frozen-advected crosswind turbulence (U across the line), force f = rho U C_D W u' chi.
Wind spectrum: Lindborg (1999) fit to Nastrom-Gage, E(k) = a k^-5/3 + b k^-3 per component.
Pods see the frozen spatial shape at speed V; cabin = secondary suspension (f_s 1 Hz, zeta 0.3);
'decoupled' modes add a low-pass reference path (skyhook, zeta 0.7) whose stroke is the bore clearance.

Reuses read-only: physics/skyloop_physics.py (atmosphere, tether_for, wall_sandwich, envelope),
regime/regime.py (breakeven, float_cost, stanton, pod, LSM/cost conventions), system/sky_model.py.
Mach convention (repo): V = M x 295.1 m/s.
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
sys.path.insert(0, os.path.join(REV, "regime"))
import regime as R  # noqa: E402  (read-only; imports physics + sky_model)

sp, sm = R.sp, R.sm
G = 9.81
A20 = sp.atm(20.0)
RHO, P_AMB, T_AMB = float(A20.rho), float(A20.p), float(A20.T)
RHO_H2 = P_AMB * sp.M_H2 / (sp.R_STAR * T_AMB)
C_REF = sp.C_REF

# ---------------------------------------------------------------- wind environment (A unless cited)
U_MEAN, U_OPS, U_DES = 18.0, 30.0, 40.0     # winter mean / operations limit / design (physics, altitude-trade)
CD = 0.6                                   # envelope and tube, crosswind
N_BV = 0.021                               # stratospheric buoyancy frequency, s^-1
LB_A, LB_B = 9.1e-4, 3.0e-10               # Lindborg (1999) fit, per component (memory, unverified)
LEVELS = {"quiet": 0.3, "typical": 1.0, "active": 5.0}   # variance multipliers (summer / mean / jet+mountain waves)
K_MIN = 2 * math.pi / 5.0e6                # largest scale considered, 5,000 km
W_RATIO = 0.05                             # w'/u' for the drag-rotation term (gravity waves are near-horizontal)
DIURNAL = 0.092                            # BoomShield lift swing (zero-pressure); superpressure cells ~0

# ---------------------------------------------------------------- ride and pod
A_LIM = 0.05 * G                           # peak lateral/vertical ride limit (physics, regime)
PK = 2.5                                   # peak / rms for the random part
F_S, Z_S = 1.0, 0.3                        # passive secondary suspension
STROKE_SEC = 0.10                          # m, secondary + magnetic gap budget (3 sigma)
Z_C = 0.7                                  # skyhook / reference-path filter damping (decoupled modes)
N_SIG_CLEAR = 4.0                          # collision margin in sigmas for decoupled stroke

# ---------------------------------------------------------------- structure and materials (A)
E_UH, RHO_UH, SIG_UH = 116e9, 975.0, 1.2e9         # UHMWPE tether (physics, BoomShield)
E_TM, RHO_TM, SIG_TM = 150e9, 1570.0, 1.1e9        # CFRP pultruded rod tension member (creep-safe)
E_UD, RHO_UD = 130e9, 1600.0                       # keel chords, UD CFRP
KEEL_WEB, KEEL_3D = 1.3, 1.5                       # webs/joints; space truss for both planes
SPACING = 3000.0                                   # tether spacing (BoomShield, physics)
EPS_BRIDLE = 0.10                                  # residual non-uniformity of slowly varying loads with bridles
T_BASE = 10e6                                      # physics: 2-10 MN line tension
D_BASE = 5.0                                       # 17 t / 50-seat maglev pod, beta ~0.2
R_ENV_PASS, R_ENV_MARG = 100.0, 150.0              # regime: R_ENV_LIM 100 m, tolerance 1.5x

# ---------------------------------------------------------------- economics (system/regime conventions)
SOFT_BASE, SOFT_LOW = 1.25 * 1.30, 1.15 * 1.15
COST_KG = {  # $/kg installed: (low, base)
    "tension": (15.0, 30.0), "keel": (30.0, 60.0), "fibre": (15.0, 35.0), "wall": (20.0, 40.0),
}
ACT_USD_M_KM = (1.0, 3.0)                          # internal guideway actuators + controls, $M/km
THRUST_USD_KW = (300.0, 800.0)
USD_KWH = R.USD_KWH
PV = R.PV
PAXKM = R.PAXKM_PER_KM


def write_csv(name, rows):
    path = os.path.join(HERE, name)
    keys = []
    for r in rows:
        keys += [k for k in r if k not in keys]
    with open(path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=keys)
        w.writeheader()
        for r in rows:
            w.writerow({k: (f"{v:.4g}" if isinstance(v, float) else v) for k, v in r.items()})
    print(f"wrote {name} ({len(rows)} rows)", flush=True)


# ================================================================== 1. gust environment
def E1(k, level=1.0):
    """One-sided 1-D spectrum of one horizontal component, m^3/s^2 (Lindborg fit)."""
    k = np.asarray(k, float)
    return level * (LB_A * k ** (-5 / 3) + LB_B * k ** (-3.0)) * (k >= K_MIN)


A2 = LB_A / (2 * 2.2405)          # 2-D isotropic form whose x-projection is E1/2 (two-sided)
B2 = LB_B / math.pi


def Phi2(kx, ky, level=1.0):
    k = np.sqrt(kx * kx + ky * ky)
    return level * (A2 * k ** (-8 / 3) + B2 * k ** (-4.0)) * (k >= K_MIN)


def band_var(l_lo, l_hi, level=1.0, n=200):
    k = np.geomspace(2 * math.pi / l_hi, 2 * math.pi / l_lo, n)
    return float(np.trapezoid(E1(k, level), k))


def f_steady(U, W):
    return 0.5 * RHO * U * U * CD * W


def stage_gust(W=175.0, R_e=86.0):
    rows = []
    edges = [100, 200, 400, 800, 1600, 3200, 6400, 12800, 25600, 51200, 102400, 204800, 409600, 819200, 1638400]
    for lo, hi in zip(edges[:-1], edges[1:]):
        r = dict(lambda_lo_km=lo / 1e3, lambda_hi_km=hi / 1e3)
        for lev, x in LEVELS.items():
            r[f"u_rms_{lev}_m_s"] = math.sqrt(band_var(lo, hi, x))
        k = 2 * math.pi / math.sqrt(lo * hi)
        chi = 1 / math.sqrt(1 + (k * W / 2) ** 2)
        r["admittance"] = chi
        r["f_lat_rms_typical_N_m"] = RHO * U_MEAN * CD * W * r["u_rms_typical_m_s"] * chi
        r["f_vert_rms_typical_N_m"] = (RHO * math.pi * R_e ** 2 * N_BV + 0.5 * RHO * U_MEAN * CD * W * W_RATIO) \
            * r["u_rms_typical_m_s"] * chi
        rows.append(r)
    # cumulative and steady loads as extra rows
    for U in (5.0, U_MEAN, U_OPS, U_DES):
        rows.append(dict(lambda_lo_km="steady", lambda_hi_km=f"U={U:g}", f_lat_rms_typical_N_m=f_steady(U, W)))
    lift = sp.h2_lift(20.0) * math.pi * R_e ** 2
    rows.append(dict(lambda_lo_km="diurnal", lambda_hi_km="9.2% lift",
                     f_vert_rms_typical_N_m=DIURNAL * lift * G))
    rows.append(dict(lambda_lo_km="patch", lambda_hi_km="eps 1e-4..1e-3, <1 km",
                     u_rms_typical_m_s=math.sqrt(0.75 * 1.5 * (1e-4) ** (2 / 3) * (1000 / (2 * math.pi)) ** (2 / 3)),
                     u_rms_active_m_s=math.sqrt(0.75 * 1.5 * (1e-3) ** (2 / 3) * (1000 / (2 * math.pi)) ** (2 / 3))))
    write_csv("gust_spectrum.csv", rows)
    return rows


# ================================================================== 2. tether field
def env_size(payload, D_sail=0.0, n_fibre=2.0):
    """H2 envelope radius for payload [kg/m], with an extra tube sail D_sail [m] and V-pair tethers
    (n_fibre legs per station, each sized for the full design side load). Returns dict."""
    lam = sp.h2_lift(20.0)
    Re = 60.0
    for _ in range(300):
        F = 0.5 * RHO * U_DES ** 2 * CD * (2 * Re + D_sail)
        t = sp.tether_for(F, 20.0)
        extra_leg = 0.5 * (n_fibre - 1) * t["m_t"] / SPACING        # half the extra legs' weight at the top
        need = sp.RESERVE * (payload + 2 * math.pi * Re * sp.FILM) + t["holddown"] + extra_leg
        Rn = math.sqrt(need / (lam * math.pi))
        if abs(Rn - Re) < 1e-5:
            break
        Re = 0.5 * Re + 0.5 * Rn
    th = math.atan(1 / sp.V_OVER_H)
    k_ax = sp.E_F * t["A"] / t["L"]
    return dict(R=Re, W=2 * Re + D_sail, tether=t, lift=lam * math.pi * Re ** 2,
                k_h=n_fibre * k_ax * math.sin(th) ** 2 + t["T"] / t["L"] * math.cos(th) ** 2,
                k_v=n_fibre * k_ax * math.cos(th) ** 2, fibre_kg_m=n_fibre * t["m_t"] / SPACING,
                h2_kg_m=math.pi * Re ** 2 * RHO_H2, film_kg_m=2 * math.pi * Re * sp.FILM)


def stage_tethers(payload):
    rows = []
    e = env_size(payload)
    t = e["tether"]
    H = t["H"]                         # design horizontal load per station
    h = 20e3
    best = E_UH / (4 * RHO_UH * h * h)  # N/m per kg of fibre at 45 deg (upper bound, any passive layout)
    lay = [
        ("vertical (pendulum)", 0.0, 1, SIG_UH),
        ("single inclined V/H 1.5", math.degrees(math.atan(1 / 1.5)), 1, SIG_UH),
        ("V-pair V/H 1.5 (baseline)", math.degrees(math.atan(1 / 1.5)), 2, SIG_UH),
        ("V-pair 45 deg", 45.0, 2, SIG_UH),
        ("V-pair 45 deg, stress 0.6 GPa", 45.0, 2, 0.6e9),
        ("V-pair 45 deg, stress 0.3 GPa", 45.0, 2, 0.3e9),
        ("cross-braced truss (pairs + diagonals)", 45.0, 2, SIG_UH),
        ("V-pair 45 deg, M55J carbon", 45.0, 2, SIG_UH),
    ]
    for name, th_deg, n, sig in lay:
        th = math.radians(th_deg)
        L = h / math.cos(th)
        E, rho = (540e9, 1910.0) if "M55J" in name else (E_UH, RHO_UH)
        if "M55J" in name:
            sig = 1.5e9
        if th_deg == 0:
            T_leg = 45e6                       # held by overlift; lateral only through pendulum T/L
            A = T_leg / sig
            k_h = T_leg / L
            k_v = E * A / L
        else:
            T_leg = H / math.sin(th)           # each leg sized for the full design side load
            A = T_leg / sig
            k_ax = E * A / L
            k_h = n * k_ax * math.sin(th) ** 2 + n * T_leg / 2 / L * math.cos(th) ** 2
            k_v = n * k_ax * math.cos(th) ** 2
        m_f = n * rho * A * L
        if "cross-braced" in name:
            m_f *= 1.4                         # diagonals between neighbouring stations: share load, add no ground path
        for s in (3000.0, 1000.0, 300.0):
            rows.append(dict(layout=name, angle_deg=th_deg, spacing_km=s / 1e3, design_H_MN=H / 1e6,
                             fibre_kg_m=m_f / s, k_h_N_m2=k_h / s, k_v_N_m2=k_v / s,
                             k_h_per_fibre=k_h / m_f, bound_45deg=best,
                             strain_at_work_pct=100 * sig / E,
                             gust_motion_m_1ms=RHO * U_MEAN * CD * e["W"] * 1.0 / (k_h / s),
                             note="spacing changes inter-tether bow only; k per km fixed by fibre mass"))
    write_csv("tether_field.csv", rows)
    return e, best


# ================================================================== 3. line response
KX = np.sort(2 * math.pi / np.geomspace(40.0, 3.0e6, 320))          # rad/m, lambda 40 m - 3,000 km
S_GRID = np.concatenate([[0.0], np.geomspace(1e-4, 1e4, 360)])


def line(payload_extra=0.0, D=D_BASE, T=T_BASE, EI=0.0, k_mult=1.0, spacing=SPACING, eps=1.0,
         level=1.0, U=U_MEAN, wall=None, vac_lift=None, tube_EI=None, p_in=5.5, k_servo=1.0, cancel=1.0):
    """Assemble the line: envelope, tethers, mass, damping and loads, for bore D."""
    if wall is None:
        wall, vac_lift, tube_EI = tube_props(D, p_in)
    payload = wall - vac_lift + 60.0 + 150.0 * 0.3 + 15.0 + payload_extra
    payload = max(payload, 20.0)
    e = env_size(payload, D_sail=D, n_fibre=2.0 * k_mult)
    Re, W = e["R"], e["W"]
    k_h, k_v = k_servo * e["k_h"] / SPACING, k_servo * e["k_v"] / SPACING   # per metre; servo winches raise it
    m = payload + e["film_kg_m"] + e["h2_kg_m"] + e["fibre_kg_m"] / 3 + RHO * math.pi * (Re ** 2 + D * D / 4)
    c = RHO * U * CD * W
    cy = cancel * RHO * U * CD * W                                   # cancel < 1: thrusters null the gust force
    cz = cancel * (RHO * math.pi * Re ** 2 * N_BV + 0.5 * RHO * U * CD * W * W_RATIO)
    eps = eps * cancel
    f_var_y = f_steady(U_OPS, W) - f_steady(U_MEAN, W)              # slowly varying, after pre-camber
    f_var_z = DIURNAL * e["lift"] * G
    return dict(payload=payload, Re=Re, W=W, D=D, k_h=k_h, k_v=k_v, m=m, c=c, cy=cy, cz=cz, T=T,
                EI=EI + tube_EI, spacing=spacing, eps=eps, level=level, U=U, f_var_y=f_var_y, f_var_z=f_var_z,
                fibre=e["fibre_kg_m"], env=e)


def tube_props(D, p_in=5.5):
    """Wall mass (kg/m), gross vacuum lift (kg/m) and own EI for bore D."""
    dp = P_AMB - p_in
    lift = (RHO - p_in / (sp.R_AIR * 250.0)) * math.pi * D * D / 4
    m2, tf, _ = sp.wall_sandwich(dp, D / 2)
    sand = (m2 + sp.LINER) * math.pi * D * sp.JOINTS
    ratio = float(np.interp(D, [10, 30, 100], [1.44, 1.27, 1.23]))     # vacuum review closure (lattice)
    latt = ratio * RHO * math.pi * D * D / 4 + sp.LINER * math.pi * D * sp.JOINTS + 5.0
    if D <= D_BASE:
        wall = sand
        EI = sp.E_QI * math.pi * (D / 2) ** 3 * 2 * tf
    else:
        wall = min(sand, latt + 80.0)    # + track-carrying frames, so D 5 sandwich ~ lattice + 80 (A)
        t_eq = 0.5 * wall / (sp.RHO_CF * math.pi * D)                  # half the wall mass works axially
        EI = sp.E_QI * math.pi * (D / 2) ** 3 * t_eq
    return wall, lift, EI


def K_of(L, k, plane="y"):
    kt = L["k_h"] if plane == "y" else L["k_v"]
    g = 1.0 / (1.0 + (k * L["spacing"] / math.pi) ** 4)
    return kt * g + L["T"] * k * k + L["EI"] * k ** 4


def spectra(L, plane="y"):
    """One-sided spatial displacement spectrum S_y(kx) [m^2/(rad/m)] and point-velocity variance."""
    kx = KX[:, None]
    ky = S_GRID[None, :] * kx
    w = L["U"] * ky
    K = K_of(L, KX, plane)[:, None]
    H2 = 1.0 / ((K - L["m"] * w * w) ** 2 + (L["c"] * w) ** 2)
    k = np.sqrt(kx * kx + ky * ky)
    chi2 = 1.0 / (1.0 + (k * L["W"] / 2) ** 2)
    cf = L["cy"] if plane == "y" else L["cz"]
    F = cf * cf * chi2 * Phi2(kx, ky, L["level"])
    Sy = 2 * 2 * np.trapezoid(H2 * F, ky, axis=1)
    Sv = 2 * 2 * np.trapezoid(w * w * H2 * F, ky, axis=1)
    sig_v = math.sqrt(float(np.trapezoid(Sv, KX)))
    return Sy, sig_v


def transfer(V, k, fs, z, kind="passive"):
    r = V * k / (2 * math.pi * fs)
    den = (1 - r * r) ** 2 + (2 * z * r) ** 2
    if kind == "passive":
        Ta = (1 + (2 * z * r) ** 2) / den
        Tr = r ** 4 / den
    else:   # skyhook
        Ta = 1.0 / den
        Tr = (r ** 4 + 4 * z * z * r * r) / den
    return Ta, Tr


def ride(L, M, mode="coupled", s_avail=None, planes=("y", "z"), cache=None, secondary=True):
    """Peak cabin acceleration and stroke use for one scheme at Mach M.
    mode 'coupled': pod on the tube's guideway (1 Hz secondary, 0.1 m stroke).
    mode 'decoupled': a low-pass reference path (actuated guideway or flying pod) with stroke s_avail,
    optimised cut-off; the cabin then also has the 1 Hz secondary (actuated guideway case)."""
    V = M * C_REF
    out = {}
    worst = None
    for pl in planes:
        Sy, sv = cache[pl] if cache and pl in cache else spectra(L, pl)
        ks = 2 * math.pi / L["spacing"]
        f_var = L["f_var_y"] if pl == "y" else L["f_var_z"]
        y_st = L["eps"] * f_var / K_of(L, ks, pl)      # periodic scallop at lambda = spacing
        Ta_s, Tr_s = transfer(V, KX, F_S, Z_S)
        Ta_s0, Tr_s0 = transfer(V, ks, F_S, Z_S)
        if mode != "coupled" and not secondary:          # flying pod: no secondary suspension
            Ta_s, Ta_s0 = np.ones_like(KX), 1.0
        if mode == "coupled":
            sa = math.sqrt(float(np.trapezoid(Sy * (V * KX) ** 4 * Ta_s, KX)))
            sr = math.sqrt(float(np.trapezoid(Sy * Tr_s, KX)))
            a_st = (V * ks) ** 2 * y_st * math.sqrt(Ta_s0)
            r_st = y_st * math.sqrt(Tr_s0)
            a_pk = PK * sa + a_st
            stroke = 3 * sr + r_st
            s_lim = STROKE_SEC
            fc = float("nan")
        else:
            cands = []
            for fc_ in np.geomspace(2e-5, 1.0, 70):
                Ta_c, Tr_c = transfer(V, KX, fc_, Z_C, "skyhook")
                Ta_c0, Tr_c0 = transfer(V, ks, fc_, Z_C, "skyhook")
                sa = math.sqrt(float(np.trapezoid(Sy * (V * KX) ** 4 * Ta_c * Ta_s, KX)))
                sr = math.sqrt(float(np.trapezoid(Sy * Tr_c, KX)))
                a_pk_ = PK * sa + (V * ks) ** 2 * y_st * math.sqrt(Ta_c0 * Ta_s0)
                st_ = N_SIG_CLEAR * sr + y_st * math.sqrt(Tr_c0)
                cands.append((a_pk_, st_, fc_))
            both = [c for c in cands if c[1] <= s_avail and c[0] <= A_LIM]
            ride_ok = [c for c in cands if c[0] <= A_LIM]
            stroke_ok = [c for c in cands if c[1] <= s_avail]
            if both:                       # feasible: smoothest ride within the stroke
                a_pk, stroke, fc = min(both, key=lambda c: c[0])
            elif ride_ok:                  # ride achievable only with more stroke than the bore gives
                a_pk, stroke, fc = min(ride_ok, key=lambda c: c[1])
            elif stroke_ok:
                a_pk, stroke, fc = min(stroke_ok, key=lambda c: c[0])
            else:
                a_pk, stroke, fc = min(cands, key=lambda c: c[0])
            s_lim = s_avail
            a_st = float("nan")
        o = dict(a_pk=a_pk, stroke=stroke, s_lim=s_lim, fc=fc, sig_v=sv, y_st=y_st,
                 sig_y=math.sqrt(float(np.trapezoid(Sy, KX))))
        out[pl] = o
        m_r = math.log10(A_LIM / a_pk)
        m_s = math.log10(s_lim / max(stroke, 1e-9))
        sc = min(m_r, m_s)
        if worst is None or sc < worst[0]:
            worst = (sc, pl, m_r, m_s)
    return out, worst


def octaves(L, M, plane="y"):
    Sy, _ = spectra(L, plane)
    V = M * C_REF
    rows = []
    edges = np.array([100, 200, 400, 800, 1600, 3200, 6400, 12800, 25600, 51200, 102400, 204800, 409600, 819200])
    for lo, hi in zip(edges[:-1], edges[1:]):
        sel = (KX >= 2 * math.pi / hi) & (KX <= 2 * math.pi / lo)
        amp = math.sqrt(2 * float(np.trapezoid(Sy[sel], KX[sel]))) if sel.sum() > 1 else 0.0
        lc = math.sqrt(lo * hi)
        kc = 2 * math.pi / lc
        Ta, _ = transfer(V, kc, F_S, Z_S)
        rows.append(dict(lambda_km=lc / 1e3, amp_m=amp, allow_m=A_LIM / ((V * kc) ** 2 * math.sqrt(Ta))))
    return rows


# ================================================================== schemes
def solve_param(build, M, mode="coupled", s_avail=None, lo=1e5, hi=1e13, n=16):
    """Smallest stiffness parameter x (log bisection) that passes ride + stroke at Mach M."""
    def ok(x):
        L = build(x)
        _, w = ride(L, M, mode, s_avail, secondary=SECONDARY)
        return w[0] >= 0, L
    good, Lh = ok(hi)
    if not good:
        return None, Lh
    a, b = math.log10(lo), math.log10(hi)
    Lb = Lh
    for _ in range(n):
        mid = 0.5 * (a + b)
        g, Lm = ok(10 ** mid)
        if g:
            b, Lb = mid, Lm
        else:
            a = mid
    return 10 ** b, Lb


SECONDARY = True


def tension_mass(T):
    return T / (SIG_TM / RHO_TM)


def keel_mass(EI, h):
    return EI * 4 * RHO_UD / (E_UD * h * h) * KEEL_WEB * KEEL_3D


def build_taut(T, extra=0.0, D=D_BASE, **kw):
    """Taut line: tension member mass feeds back into the envelope (2 iterations on sail)."""
    return line(payload_extra=tension_mass(T) + extra, D=D, T=T, eps=EPS_BRIDLE, **kw)


def build_keel(EI, h, extra=0.0, D=D_BASE, **kw):
    return line(payload_extra=keel_mass(EI, h) + extra, D=D, EI=EI, eps=EPS_BRIDLE, **kw)


def build_stiff(mult, extra=0.0, **kw):
    return line(payload_extra=extra, k_mult=mult, eps=EPS_BRIDLE, **kw)


# ================================================================== 4. active trim
def active_power(L, M=2.0):
    """Power per 1,000 km for actuator options, and their bandwidth."""
    Sy, sv = spectra(L, "y")
    sig_y = math.sqrt(float(np.trapezoid(Sy, KX)))
    # gust load rms seen by the line, scales up to 100 km
    var_u = band_var(L["W"] * 2, 1e5, L["level"])
    f_rms = L["cy"] * math.sqrt(var_u)
    rows = []
    # (i) envelope thrusters (ideal actuator disc, area 2R per m) cancelling f_rms
    A = L["W"]
    P_ideal = f_rms ** 1.5 / math.sqrt(2 * RHO * A)
    rows.append(dict(option="envelope thrusters / rotors", power_MW_per_1000km=P_ideal * 1e6 / 1e6 / 0.6 * 1.25,
                     bandwidth_Hz="0.001-0.05", stroke_m="n/a",
                     note="ideal disc power / 0.6 figure of merit x 1.25 peak-to-rms; cancels gusts at the source"))
    # (ii) tether winches: reel to cancel elastic stretch dT*L/EA; reactive power T0*v, 15 % net loss
    t = L["env"]["tether"]
    stretch_v = sv   # top velocity ~ line velocity (m/s), winch must match it
    n_st = 2 * 1e6 / SPACING
    P_react = t["T"] * stretch_v * n_st
    rows.append(dict(option="tether winches (servo stiffness)", power_MW_per_1000km=0.15 * P_react / 1e6,
                     bandwidth_Hz="0.001-0.05 (wave transit 2.2 s)", stroke_m=f"{3*sig_y:.0f}",
                     reactive_GW=P_react / 1e9,
                     note="holds tether tops only; between-tether band unchanged"))
    # (iii) internal guideway on long-stroke actuators: weight passively balanced, dynamic part + losses
    m_g = 60.0 + 45.0 + 20.0                                    # guideway + stator avg + actuators, kg/m
    a_path = 0.05                                               # m/s^2 rms, guideway path following
    P_dyn = m_g * a_path * sv + 0.2 * m_g * G * sv * 0.1        # inertia + 10 % unbalanced weight at 20 % loss
    rows.append(dict(option="internal guideway on actuators (moving guideway)",
                     power_MW_per_1000km=P_dyn * 1e6 / 1e6, bandwidth_Hz="0.01-0.3 (tube motion), passive >1 Hz",
                     stroke_m=f"bore-limited; tube sigma_y {sig_y:.1f} m", note="needs the stroke; see bigbore"))
    rows.append(dict(option="(no weight balance) same, actuators lift the guideway",
                     power_MW_per_1000km=0.2 * m_g * G * sv * 1e6 / 1e6 * 2,
                     bandwidth_Hz="same", stroke_m="same", note="20 % loss on reactive F*v"))
    for r in rows:
        r.update(Mach=M, sig_v_m_s=sv, sig_y_m=sig_y, f_rms_N_m=f_rms)
    return rows


# ================================================================== 5. big bore
def pod_flight(D, M, m_pod=R.REF_MASS, seats=50):
    """Free-flying winged pod in a bore D: wing span 0.5 D, AR 2.5 (A)."""
    V = M * C_REF
    b = 0.5 * D
    S = b * b / 2.5
    Mg = V / math.sqrt(1.4 * 287.05 * 250.0)
    LD = 0.7 * 4 * (Mg + 3) / Mg                       # 0.7 x Kuchemann (low Re, tube) (A)
    CL = 0.3
    q = m_pod * G / (S * CL)
    p = q / (0.7 * Mg * Mg)
    rho = p / (287.05 * 250.0)
    pd = R.pod(m_pod)
    Swet = pd["S_wet"] + 2.1 * S
    Re = rho * V * pd["L"] / R.mu(250.0)
    T_aw = 250.0 * (1 + 0.85 * 0.2 * Mg * Mg)
    T_sk = 500.0                                       # Ti / hot structure skin (A)
    Q_in = R.stanton(Re) * rho * V * R.CP * max(0.0, T_aw - T_sk) * Swet + 150.0 * seats
    Q_out = 0.85 * R.SIGMA * (T_sk ** 4 - 250.0 ** 4) * Swet
    pax = seats * 0.8
    E = (m_pod * G / LD) / 0.85 / pax / 3.6           # Wh/pax-km, lift-induced drag only
    a_coast = G / LD
    L_boost = 0.1 * V * V / a_coast / 1e3              # km for a 5 % speed sag
    return dict(span_m=b, wing_m2=S, LD=LD, p_req_Pa=p, p_frac_amb=p / P_AMB, heat_ratio=Q_in / Q_out,
                E_Wh_paxkm=E, boost_spacing_km=L_boost, breach_overload_g=P_AMB / p, power_MW=m_pod * G / LD * V / 0.85 / 1e6)


def stage_bigbore():
    rows = []
    for D in (D_BASE, 8.0, 20.0, 40.0, 80.0):
        for mode, p_in in (("maglev, active guideway", 5.5), ("flying pod", 500.0)):
            wall, lift, EI = tube_props(D, p_in)
            m2, _, _ = sp.wall_sandwich(P_AMB - p_in, D / 2)
            sand = (m2 + sp.LINER) * math.pi * D * sp.JOINTS
            L = line(D=D, p_in=p_in, eps=EPS_BRIDLE)
            vol = math.pi * D * D / 4 * 1000
            leak_cm2 = 1.0 * D / 3.0                          # physics 1 cm2/km at D 3, scaled with wall area
            S_req = leak_cm2 * R.q_hole(P_AMB, T_AMB, 250.0) / p_in
            r = dict(D_m=D, mode=mode, p_in_Pa=p_in, wall_kg_m=wall, wall_sandwich_kg_m=sand,
                     vac_lift_kg_m=lift, self_lift_frac=lift / wall, tube_EI_Nm2=EI,
                     payload_kg_m=L["payload"], env_R_m=L["Re"], sail_W_m=L["W"], m_eff_kg_m=L["m"],
                     k_h_N_m2=L["k_h"], gust_f_N_per_ms=L["cy"],
                     stroke_avail_m=(D - 5.0) / 2 - 0.5 if mode.startswith("maglev") else D / 4 - 0.5,
                     volume_m3_km=vol, pumpdown_days_10m3s=vol / 10.0 * math.log(P_AMB / p_in) / 86400,
                     leak_cm2_km=leak_cm2, pump_kW_km=10 * S_req,
                     breach_pV_tTNT_km=(P_AMB - p_in) * vol / 4.184e9)
            if mode == "flying pod":
                for M in (1.5, 2.0, 3.0, 4.0):
                    pf = pod_flight(D, M)
                    r.update({f"M{M:g}_{k}": v for k, v in pf.items()
                              if k in ("p_req_Pa", "heat_ratio", "E_Wh_paxkm", "boost_spacing_km", "breach_overload_g")})
            rows.append(r)
    write_csv("bigbore.csv", rows)
    return rows


# ================================================================== 6. feasibility map
MACHS = [1.0, 1.25, 1.5, 1.75, 2.0, 2.5, 3.0, 3.5, 4.0]


def status(m, tol):
    return "pass" if m >= 0 else ("marginal" if m >= -tol else "fail")


def crit_speed(L):
    kt = L["k_h"]
    return math.sqrt((L["T"] + 2 * math.sqrt(max(kt * L["EI"], 0.0))) / L["m"])


_E0 = {}


def _e0():
    _E0.update(env_size(line()["payload"], D_sail=D_BASE))
    return _E0


def economics(M, L, add_kg, add_usd_m=(0.0, 0.0), power_MW=0.0, extra_wh=0.0, D=D_BASE):
    """B/C proxy as in regime.py: break-even avg tube cost / (floating cost + adders)."""
    BE = R.breakeven(M)
    e0 = _E0 if _E0 else _e0()
    env_scale = L["Re"] / e0["R"] - 1
    teth_scale = L["W"] / e0["W"] - 1
    out = {}
    for i, lvl in enumerate(("low", "base")):
        base = R.float_cost(M) * (0.28 if lvl == "low" else 1.0)
        soft = SOFT_LOW if lvl == "low" else SOFT_BASE
        adders = sum(kg * COST_KG[k][i] for k, kg in add_kg.items()) / 1e3        # $M/km
        adders += (9.8 if lvl == "base" else 2.6) * env_scale + (11.0 if lvl == "base" else 3.0) * teth_scale
        adders += add_usd_m[i]
        c_en = (power_MW / 1e3 * 8760 * USD_KWH * 1e3 * PV / 1e6) + extra_wh / 1e3 * PAXKM * USD_KWH * PV / 1e6
        out[lvl] = base + adders * soft + c_en
    return BE / out["base"], BE / out["low"], out["base"], out["low"]


def evaluate(scheme, M, D=D_BASE):
    """Return a dict of margins/status for one scheme x Mach (x bore)."""
    V = M * C_REF
    mg, tol, info = {}, {}, {}
    add_kg, add_usd, power, extra_wh = {}, (0.0, 0.0), 0.0, 0.0
    mode, s_av = "coupled", None
    if scheme == "bare tethers":
        L = line(eps=1.0)
    elif scheme == "bracing + bridles":
        L = line(eps=EPS_BRIDLE)
        add_kg["fibre"] = 0.4 * L["fibre"]
    elif scheme == "stiff tethers x4":
        L = build_stiff(4.0)
        add_kg["fibre"] = 3 * line()["fibre"]
    elif scheme in ("keel h80", "keel h150"):
        h = 80.0 if scheme == "keel h80" else 150.0
        EI, L = solve_param(lambda x: build_keel(x, h), M, lo=1e11, hi=1e17)
        if EI is None:
            EI = 1e17
        add_kg["keel"] = keel_mass(EI, h)
        info.update(EI=EI, keel_kg_m=keel_mass(EI, h))
    elif scheme == "taut line":
        T, L = solve_param(lambda x: build_taut(x), M, lo=1e7, hi=3e10)
        if T is None:
            T = 3e10
        add_kg["tension"] = tension_mass(T)
        info.update(T_GN=T / 1e9, tension_kg_m=tension_mass(T))
    elif scheme == "active thrusters":
        L0 = line(eps=EPS_BRIDLE)
        pr = active_power(L0, M)
        f_sc = EPS_BRIDLE * L0["f_var_y"] / math.sqrt(2)                 # residual scallop force to null too
        f_g = pr[0]["f_rms_N_m"]
        power = pr[0]["power_MW_per_1000km"] * ((f_g ** 2 + f_sc ** 2) / f_g ** 2) ** 0.75
        L = line(eps=EPS_BRIDLE, cancel=0.1)                             # 90 % cancellation (A)
        add_usd = tuple(power * 1e3 * c / 1e6 / 1e3 for c in THRUST_USD_KW)    # $M per km
    elif scheme == "active winches":
        L = line(eps=EPS_BRIDLE, k_servo=10.0)                           # servo stiffness 10x at the tether tops
        power = active_power(line(eps=EPS_BRIDLE), M)[1]["power_MW_per_1000km"]
        add_usd = (1.0, 3.0)
    elif scheme.startswith("active guideway") or scheme.startswith("taut + active"):
        mode, s_av = "decoupled", (D - 5.0) / 2 - 0.5
        if scheme.startswith("taut"):
            T, L = solve_param(lambda x: build_taut(x, extra=20.0, D=D), M, "decoupled", s_av, lo=1e7, hi=3e10)
            if T is None:
                T = 3e10
            add_kg["tension"] = tension_mass(T)
            info.update(T_GN=T / 1e9, tension_kg_m=tension_mass(T))
        else:
            L = line(D=D, payload_extra=20.0, eps=EPS_BRIDLE)
        power = active_power(L, M)[2]["power_MW_per_1000km"]
        add_usd = ACT_USD_M_KM
    elif scheme.startswith("flying pod"):
        mode, s_av = "decoupled", D / 4 - 0.5
        L = line(D=D, p_in=500.0, payload_extra=-80.0, eps=EPS_BRIDLE)     # no levitation track, boost stations
        pf = pod_flight(D, M)
        info.update(p_req=pf["p_req_Pa"], E_fly=pf["E_Wh_paxkm"], heat=pf["heat_ratio"])
        extra_wh = max(0.0, pf["E_Wh_paxkm"] - 40.0)
    else:
        raise KeyError(scheme)
    if D > D_BASE:
        w0 = tube_props(D_BASE)[0]
        add_kg["wall"] = add_kg.get("wall", 0.0) + max(0.0, tube_props(D)[0] - w0)

    rd, worst = ride(L, M, mode, s_av, secondary=not scheme.startswith("flying pod"))
    mg["ride"] = worst[2]
    tol["ride"] = math.log10(2)
    mg["clearance"] = worst[3]
    tol["clearance"] = math.log10(1.5)
    mg["lift"] = math.log10(R_ENV_PASS / L["Re"])
    tol["lift"] = math.log10(R_ENV_MARG / R_ENV_PASS)
    vc = crit_speed(L)
    ratio = V / vc
    mg["crit_speed"] = math.log10(ratio / 1.4) if ratio >= 1 else math.log10(0.7 / ratio)
    tol["crit_speed"] = math.log10(1.4 / 1.15) if ratio >= 1 else math.log10(0.85 / 0.7)
    # power (active schemes) and support (flying pods)
    if power > 0:
        mg["power"] = math.log10(50.0 / power)
        tol["power"] = math.log10(6.0)
    if scheme.startswith("flying pod"):
        m_E = math.log10(R.E_LIM / info["E_fly"])
        m_h = math.log10(1 / info["heat"])
        m_p = math.log10(1000.0 / info["p_req"]) if info["p_req"] > 1000 else -0.01   # propulsion unproven
        mg["support"] = min(m_E, m_h, m_p)
        tol["support"] = math.log10(R.E_JET / R.E_LIM)
    # safety: floating evacuation is marginal everywhere (regime); tension stored energy; winged pods in a breach
    m_s = -0.05
    tol_s = math.log10(4)
    if "T_GN" in info:
        e_mj = info["T_GN"] * 1e9 * (SIG_TM / E_TM) / 2 / 1e6        # MJ per m stored
        info["strain_MJ_m"] = e_mj
        m_s = min(m_s, math.log10(0.5 / e_mj) if e_mj > 0.5 else -0.05)     # BoomShield tether ~0.15 MJ/m
    if scheme.startswith("flying pod"):
        m_s = min(m_s, math.log10(5.0 / (P_AMB / info["p_req"])))
    mg["safety"], tol["safety"] = m_s, tol_s
    bc_b, bc_l, c_b, c_l = economics(M, L, add_kg, add_usd, power, extra_wh, D)
    mg["economics"] = math.log10(bc_b)
    tol["economics"] = math.log10(bc_l / bc_b) if bc_l >= 1 else -1e-9

    out = dict(scheme=scheme, Mach=M, D_m=D, env_R_m=L["Re"], sail_W_m=L["W"], payload_kg_m=L["payload"],
               a_pk_y=rd["y"]["a_pk"], a_pk_z=rd["z"]["a_pk"], stroke_y=rd["y"]["stroke"], stroke_z=rd["z"]["stroke"],
               stroke_lim=rd["y"]["s_lim"], fc_Hz=rd["y"]["fc"], sig_y_m=rd["y"]["sig_y"], sig_z_m=rd["z"]["sig_y"],
               sig_v_m_s=rd["y"]["sig_v"], scallop_y_m=rd["y"]["y_st"], scallop_z_m=rd["z"]["y_st"],
               V_over_vcrit=ratio, power_MW_1000km=power, BC_base=bc_b, BC_low=bc_l,
               cost_base_usd_m_km=c_b, cost_low_usd_m_km=c_l, worst_plane=worst[1], **info)
    sts = {}
    for d in mg:
        out[f"m_{d}"] = mg[d]
        sts[d] = status(mg[d], tol[d])
        out[f"s_{d}"] = sts[d]
    tech = [d for d in mg if d != "economics"]
    ts = [sts[d] for d in tech]
    out["tech_status"] = "fail" if "fail" in ts else ("marginal" if "marginal" in ts else "pass")
    score = {d: (mg[d] if mg[d] >= 0 else mg[d] - (0 if mg[d] >= -tol[d] else 10)) for d in tech}
    out["binding"] = min(score, key=score.get)
    allst = ts + [sts["economics"]]
    out["full_status"] = "fail" if "fail" in allst else ("marginal" if "marginal" in allst else "pass")
    if out["tech_status"] != "fail" and sts["economics"] == "fail":
        out["binding"] = "economics"
    return out


SCHEMES = [("bare tethers", None), ("bracing + bridles", None), ("stiff tethers x4", None),
           ("keel h80", None), ("keel h150", None), ("taut line", None), ("active thrusters", None), ("active winches", None),
           ("active guideway D8", 8.0), ("active guideway D20", 20.0), ("active guideway D40", 40.0),
           ("active guideway D80", 80.0), ("taut + active D8", 8.0), ("taut + active D20", 20.0),
           ("flying pod D20", 20.0), ("flying pod D40", 40.0), ("flying pod D80", 80.0)]


def stage_map():
    rows = []
    path = os.path.join(HERE, "feasibility_map.csv")
    done = set()
    if os.path.exists(path) and "--fresh" not in sys.argv:          # resume after a restart
        rows = list(csv.DictReader(open(path)))
        done = {(r["scheme"], float(r["Mach"])) for r in rows}
    for name, D in SCHEMES:
        for M in MACHS:
            if (name, M) in done:
                continue
            r = evaluate(name, M, D or D_BASE)
            rows.append(r)
            print(f"  {name:22s} M{M:<5g} {r['tech_status']:8s} bind={r['binding']:11s} "
                  f"a={max(r['a_pk_y'], r['a_pk_z'])/G:.3g} g  R={r['env_R_m']:.0f}  B/C={r['BC_base']:.2f}/{r['BC_low']:.2f}",
                  flush=True)
            write_csv("feasibility_map.csv", rows)
    return rows


# ================================================================== keel table
def stage_keel(payload):
    rows = []
    L0 = line(eps=EPS_BRIDLE)
    for M in (1.5, 2.0, 3.0, 4.0, 6.0):
        V = M * C_REF
        for s in (500.0, 1000.0, 3000.0, 5000.0):
            # band load at lambda = 2 s (octave), and slowly varying residual at lambda = s
            var = band_var(1.41 * s, 2.83 * s)
            f_g = L0["cy"] * math.sqrt(var)
            lam = 2 * s
            EI_g = V * V * PK * f_g * lam ** 2 / (4 * math.pi ** 2 * A_LIM)
            EI_s = V * V * EPS_BRIDLE * L0["f_var_y"] * s * s / (4 * math.pi ** 2 * A_LIM)
            r_s = V / s / F_S
            if r_s > 2:          # scallop above the suspension: limit to 1/3 of the stroke instead
                EI_s = EPS_BRIDLE * L0["f_var_y"] / (STROKE_SEC / 3) / (2 * math.pi / s) ** 4
            EI = max(EI_g, EI_s)
            r = dict(Mach=M, span_km=s / 1e3, f_gust_N_m=f_g, EI_gust=EI_g, EI_steady=EI_s, EI_req=EI)
            for h in (20.0, 50.0, 80.0, 150.0):
                r[f"mass_t_m_h{h:g}"] = keel_mass(EI, h) / 1e3
            r["T_equiv_GN"] = EI * (2 * math.pi / lam) ** 2 / 1e9
            r["tension_alt_t_m"] = tension_mass(r["T_equiv_GN"] * 1e9) / 1e3
            rows.append(r)
    write_csv("keel_sizing.csv", rows)


def stage_motion():
    rows = []
    cases = [("bare tethers, T 10 MN", line(eps=1.0)),
             ("taut line T 0.3 GN", build_taut(3e8)),
             ("taut line T 1 GN", build_taut(1e9)),
             ("keel EI 1e14 (h 150)", build_keel(1e14, 150.0))]
    for name, L in cases:
        for M in (2.0, 4.0):
            for o in octaves(L, M):
                o.update(case=name, Mach=M)
                rows.append(o)
    # spectrum levels for the bare case
    for lev in ("quiet", "active"):
        L = line(eps=1.0, level=LEVELS[lev])
        for o in octaves(L, 2.0):
            o.update(case=f"bare tethers, {lev} spectrum", Mach=2.0)
            rows.append(o)
    write_csv("motion_octaves.csv", rows)


def stage_sensitivity():
    """Taut-line tension needed vs bridle residual, spectrum level, mean wind; and the minimum bore for
    the decoupled (actuated-guideway / flying-pod) schemes per Mach."""
    rows = []
    for M in (1.5, 2.0, 3.0):
        for eps in (0.03, 0.1, 0.3, 1.0):
            for lev in ("quiet", "typical", "active"):
                T, L = solve_param(lambda x: line(payload_extra=tension_mass(x), T=x, eps=eps, level=LEVELS[lev]),
                                   M, lo=1e7, hi=3e10)
                rows.append(dict(study="taut line T_req", Mach=M, eps=eps, level=lev,
                                 T_GN=(T or float("nan")) / 1e9, env_R_m=L["Re"],
                                 tension_kg_m=tension_mass(T) if T else float("nan")))
    for M in (1.5, 2.0, 3.0, 4.0):
        for base, T in (("bare tethers", T_BASE), ("taut 0.1 GN", 1e8), ("taut 0.3 GN", 3e8)):
            for kind in ("actuated guideway", "flying pod"):
                def ok(D):
                    sa = (D - 5.0) / 2 - 0.5 if kind == "actuated guideway" else D / 4 - 0.5
                    pin = 5.5 if kind == "actuated guideway" else 500.0
                    L = line(D=D, T=T, payload_extra=tension_mass(T) if T > T_BASE else 0.0, eps=EPS_BRIDLE, p_in=pin)
                    _, w = ride(L, M, "decoupled", sa, secondary=(kind == "actuated guideway"))
                    return w[0] >= 0, L
                a, b = 7.0, 160.0
                gb, Lb = ok(b)
                if not gb:
                    rows.append(dict(study="min bore", Mach=M, base=base, kind=kind, D_min_m=float("nan")))
                    continue
                for _ in range(14):
                    mid = 0.5 * (a + b)
                    g, Lm = ok(mid)
                    if g:
                        b, Lb = mid, Lm
                    else:
                        a = mid
                rows.append(dict(study="min bore", Mach=M, base=base, kind=kind, D_min_m=b, env_R_m=Lb["Re"]))
        print(f"  sensitivity M{M} done", flush=True)
        write_csv("sensitivity.csv", rows)
    write_csv("sensitivity.csv", rows)


def main():
    L0 = line()
    print(f"baseline: R_env {L0['Re']:.1f} m, W {L0['W']:.0f} m, k_h {L0['k_h']:.1f} N/m2, k_v {L0['k_v']:.1f}, "
          f"m_eff {L0['m']:.0f} kg/m, payload {L0['payload']:.0f}, fibre {L0['fibre']:.0f} kg/m", flush=True)
    stage_gust(W=L0["W"], R_e=L0["Re"])
    stage_tethers(L0["payload"])
    stage_motion()
    stage_keel(L0["payload"])
    act = []
    for M in (2.0, 4.0):
        act += active_power(line(eps=EPS_BRIDLE), M)
    write_csv("active_trim.csv", act)
    stage_bigbore()
    if "--no-sens" not in sys.argv:
        stage_sensitivity()
    stage_map()


if __name__ == "__main__":
    if "--figs" not in sys.argv:
        main()
    if os.path.exists(os.path.join(HERE, "make_figures.py")):
        import importlib.util  # regime/ also has a make_figures.py on sys.path
        spec = importlib.util.spec_from_file_location("stab_figs", os.path.join(HERE, "make_figures.py"))
        mf = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(mf)
        mf.main()
