#!/usr/bin/env python3
"""Skyloop aerothermal heating study: winged pods flying inside a wide (20-80 m) low-pressure tube, Mach 3-8.

Screening level.  Maps pod heating (stagnation, leading edge, acreage), cooling (radiative equilibrium to the
tube wall, heat sinks), and tube-side heating (gas, wall/films, guideway) over Mach x tube pressure, for
several wing loadings and bores, and marks the binding constraint at each point.

Run:   python3 skyloop-review/heating/heating.py           (~1-3 min, one process)
       python3 skyloop-review/heating/make_figures.py      (figures from the CSVs)
Out:   heating_grid.csv      W/S x bore x Mach x pressure: flight state, temperatures, heat loads, tube-side
                             numbers, margins, status, binding constraint
       heating_bestld.csv    best-L/D pressure per W/S x Mach (the efficient flight line)
       heating_scaling.csv   heat per trip vs Mach at fixed q and on the best-L/D line (the q*V*St test)
       heating_radius.csv    nose / leading-edge temperature vs radius (0.5 cm - 1 m), Knudsen number
       heating_sink.csv      heat-sink mass fraction vs Mach, pressure, trip length, material set
       heating_tube.csv      gas dT per pass, relaxation time, max pods/h, pumping-as-cooling, film, guideway
       heating_levers.csv    how design levers move the feasible Mach band
       heating_summary.csv   feasible band per W/S x bore

Reuses read-only: regime/regime.py (pod geometry, Sutherland viscosity, C_REF).
Mach convention (repo): V = M x 295.07 m/s; tube gas = air at 216.65 K (20 km ambient).  (A) = assumption.
"""
from __future__ import annotations

import csv
import math
import os
import sys

import numpy as np
from scipy.optimize import brentq

HERE = os.path.dirname(os.path.abspath(__file__))
REV = os.path.dirname(HERE)
sys.dont_write_bytecode = True
sys.path.insert(0, os.path.join(REV, "regime"))
import regime as R  # noqa: E402  (read-only import; its main() is not run)

# ----------------------------------------------------------------------------- constants
G, R_E = 9.81, 6.371e6
GAM, R_AIR, CP, SIGMA = 1.4, 287.05, 1005.0, 5.670e-8
CV = CP / GAM
PR = 0.71
T_GAS = 216.65                       # tube gas = 20 km ambient
A_GAS = math.sqrt(GAM * R_AIR * T_GAS)
C_REF = R.C_REF                      # 295.07 m/s
T_ENV = 250.0                        # K, tube wall seen by the pod (inside a sunlit H2 envelope) (A)
EPS = 0.85                           # pod surface emissivity (A)
K_SG = 1.7415e-4                     # Sutton-Graves, air, SI (W/m2, kg/m3, m, m/s) [Sutton & Graves 1971]

# ----------------------------------------------------------------------------- pod design (A)
BASE = dict(
    m=3500.0,            # kg, coordinator's pod
    AR=2.0,              # delta-like wing, LE sweep atan(4/AR) = 63 deg
    R_n=0.10,            # m, nose radius
    R_le=0.03,           # m, wing leading-edge radius
    T_lim_skin=1250.0,   # K, Inconel-617 metallic TPS reuse limit (NASA: 1800-1900 F)
    T_lim_le=1700.0,     # K, C/SiC leading edge (ONERA: 1100-1500 C demonstrated; low-pO2 active oxidation caveat)
    h_sink=2.5e6,        # J/kg, water evaporated and vented (A)
    f_sink_ok=0.05,      # heat-sink mass fraction: pass <= 5 %, marginal <= 15 %
    f_sink_marg=0.15,
    L_trip=4000e3,       # m, NYC-LAX
    a_acc=0.3 * G,       # launch / brake acceleration
    k_ins=0.04, t_ins=0.05, T_cab=295.0,   # cabin insulation (fibrous at 1-300 Pa) (A)
    w_seat=150.0,        # W per seat (people + avionics), as regime.py
    re_tr=1.0e6,         # transition Reynolds number (edge, x) (A)
    dT_gas=50.0,         # K, allowed steady bore-mean gas temperature rise (A)
    dT_pass_ok=22.0,     # K, per-pass bore-mean rise: pass (Pi = 0.25) ...
    dT_pass_marg=87.0,   # ... marginal (Pi = 1)
    f_mix=1.0,           # multiplier on gas-to-wall Nusselt number (wake stirring, fans) (A)
    N_req=120.0,         # pods/h per tube needed for 'pass' (30 s headway; 600 seats/h at 5 seats); marginal 30/h (A)
    N_marg=30.0,
    alpha_le_ok=25.0,    # deg, lower-surface incidence: pass ...
    alpha_le_marg=40.0,  # ... marginal (Newtonian C_L peak ~55 deg)
    span_frac=0.5,       # span <= 0.5 D passes, <= 0.8 D marginal
    eta_lsm=0.90,        # motor-strip efficiency (copper + iron) (A)
    w_strip=2.0,         # m, motor-strip width (A)
    T_strip_lim=423.0,   # K, class F/H coil insulation (155-180 C) (A)
    film_ae=0.5,         # alpha_solar / eps_IR of the envelope film (A)
    T_film_lim=353.0,    # K, PE film (80 C); UHMWPE ~120 C as marginal
    T_film_marg=393.0,
)

WS_LIST = [50.0, 100.0, 150.0, 300.0, 500.0]
D_LIST = [20.0, 40.0, 80.0]
M_GRID = np.round(np.arange(3.0, 8.001, 0.1), 2)
P_GRID = np.geomspace(1.0, 1000.0, 46)        # 1-300 Pa asked; extended to 1 kPa to show the heating ceiling


# ----------------------------------------------------------------------------- gas properties
def mu(T):
    return R.mu(T)


def k_gas(T):
    return 0.0241 * (T / 273.15) ** 0.81


def mfp(p, T):
    return mu(T) / p * math.sqrt(math.pi * R_AIR * T / 2.0)


def freestream(M, p):
    V = M * C_REF
    rho = p / (R_AIR * T_GAS)
    h0 = CP * T_GAS + 0.5 * V * V
    return V, rho, 0.5 * rho * V * V, h0 / CP


# ----------------------------------------------------------------------------- shocks
def beta_max(M):
    g = GAM
    s2 = ((g + 1) * M * M - 4 + math.sqrt((g + 1) * ((g + 1) * M ** 4 + 8 * (g - 1) * M * M + 16))) / (4 * g * M * M)
    return math.asin(math.sqrt(min(1.0, s2)))


def theta_of_beta(M, b):
    return math.atan(2 / math.tan(b) * (M * M * math.sin(b) ** 2 - 1) / (M * M * (GAM + math.cos(2 * b)) + 2))


def compression(M, T, pr):
    """Oblique shock giving static-pressure ratio pr.  Returns (T2, V-ratio, deflection deg) or None if detached."""
    Mn1_2 = 1 + (pr - 1) * (GAM + 1) / (2 * GAM)
    if Mn1_2 > M * M:
        return None
    b = math.asin(math.sqrt(Mn1_2) / M)
    if b > beta_max(M):
        return None
    rr = (GAM + 1) * Mn1_2 / ((GAM - 1) * Mn1_2 + 2)       # rho2/rho1
    T2 = T * pr / rr
    Vt, Vn2 = math.cos(b), math.sin(b) / rr
    return T2, math.hypot(Vt, Vn2), math.degrees(theta_of_beta(M, b))


def pitot_ratio(M):
    g = GAM
    return ((g + 1) ** 2 * M * M / (4 * g * M * M - 2 * (g - 1))) ** (g / (g - 1)) * (1 - g + 2 * g * M * M) / (g + 1)


# ----------------------------------------------------------------------------- heating correlations
def plate(p_e, T_e, V_e, L, T_w, re_tr):
    """Mean heat flux and shear over a flat plate of length L (Eckert reference temperature, laminar ->
    turbulent at Re_x = re_tr, weak viscous-interaction correction on the laminar part)."""
    rho_e = p_e / (R_AIR * T_e)
    M_e = V_e / math.sqrt(GAM * R_AIR * T_e)
    Ts = T_e * (1 + 0.032 * M_e ** 2 + 0.58 * (T_w / T_e - 1))
    rho_s, mu_s = p_e / (R_AIR * Ts), mu(Ts)
    T_aw_l = T_e + math.sqrt(PR) * V_e ** 2 / (2 * CP)
    T_aw_t = T_e + PR ** (1 / 3) * V_e ** 2 / (2 * CP)
    a = rho_s * V_e / mu_s                                   # Re*/x
    x_t = min(L, re_tr * mu(T_e) / (rho_e * V_e))
    # laminar: q(x) = 0.332 Pr^-2/3 rho* V cp (Taw-Tw) (a x)^-1/2 ; mean over [0, x_t] -> x 2
    c_l = 0.332 * PR ** (-2 / 3) * rho_s * V_e * CP * (T_aw_l - T_w) / math.sqrt(a)
    t_l = 0.332 * rho_s * V_e ** 2 / math.sqrt(a)            # tau = Cf*/2 rho* V^2
    Q_l = 2 * c_l * math.sqrt(x_t)
    F_l = 2 * t_l * math.sqrt(x_t)
    # weak viscous interaction (insulated-plate form, p/p_e = 1 + 0.31 chi + 0.05 chi^2), heating ~ sqrt(p)
    C = (p_e / (R_AIR * T_w)) * mu(T_w) / (rho_e * mu(T_e))
    chi = M_e ** 3 * math.sqrt(C / max(rho_e * V_e * 0.25 * x_t / mu(T_e), 1.0))
    f_vi = math.sqrt(1 + 0.31 * chi + 0.05 * chi * chi)
    Q_l *= f_vi
    F_l *= f_vi
    Q_t = F_t = 0.0
    if x_t < L:
        c_t = 0.0296 * PR ** (-2 / 3) * rho_s * V_e * CP * (T_aw_t - T_w) / a ** 0.2
        t_t = 0.0296 * rho_s * V_e ** 2 / a ** 0.2
        Q_t = c_t * (L ** 0.8 - x_t ** 0.8) / 0.8
        F_t = t_t * (L ** 0.8 - x_t ** 0.8) / 0.8
    return (Q_l + Q_t) / L, (F_l + F_t) / L, dict(M_e=M_e, chi=chi, turb_frac=1 - x_t / L, Re_L=rho_e * V_e * L / mu(T_e))


def q_stag(rho, V, h0T, R_n, T_w, M=None, p=None):
    """Sphere stagnation point: Sutton-Graves cold wall -> hot wall enthalpy ratio, bridged to the
    free-molecular limit 0.5 rho V^3 (full accommodation)."""
    qc = K_SG * math.sqrt(rho / R_n) * V ** 3 * max(0.0, h0T - T_w) / (h0T - 300.0)
    qfm = 0.5 * rho * V ** 3 * max(0.0, 1 - T_w / h0T)
    return qc / math.sqrt(1 + (qc / qfm) ** 2) if qfm > 0 else 0.0


def q_fay_riddell(M, p, R_n, T_w):
    """Perfect-gas Fay-Riddell stagnation heating (Le = 1, no dissociation), cross-check only."""
    V, rho, q, T0 = freestream(M, p)
    p_s = p * pitot_ratio(M)
    rho_s = p_s / (R_AIR * T0)
    rho_w = p_s / (R_AIR * T_w)
    due = math.sqrt(2 * (p_s - p) / rho_s) / R_n
    return 0.763 * PR ** -0.6 * (rho_w * mu(T_w)) ** 0.1 * (rho_s * mu(T0)) ** 0.4 * CP * (T0 - T_w) * math.sqrt(due)


def rad_eq(qfun, T_hi):
    """Wall temperature where convective input = grey radiation to the tube wall."""
    f = lambda Tw: qfun(Tw) - EPS * SIGMA * (Tw ** 4 - T_ENV ** 4)
    lo, hi = T_ENV + 1e-3, max(T_hi - 1e-3, T_ENV + 2e-3)
    if f(lo) <= 0:
        return lo
    if f(hi) >= 0:
        return hi
    return brentq(f, lo, hi, xtol=0.05)


def q_net(qfun, Tw):
    """Excess convective flux a surface held at Tw must sink (W/m2), >= 0."""
    return max(0.0, qfun(Tw) - EPS * SIGMA * (Tw ** 4 - T_ENV ** 4))


# ----------------------------------------------------------------------------- the pod
def pod_state(M, p, WS, d=BASE):
    """Flight state and heating for a winged pod at Mach M, tube pressure p, wing loading WS."""
    V, rho, q, T0 = freestream(M, p)
    pd_ = R.pod(d["m"])
    W = d["m"] * G
    relief = 1 - V * V / (G * R_E)
    S = W / WS
    b = math.sqrt(d["AR"] * S)
    cbar = S / b
    lam = math.atan(4.0 / d["AR"])
    WSe = WS * relief
    out = dict(Mach=M, V=V, p_Pa=p, WS_Pa=WS, q_Pa=q, CL=WSe / q, T0=T0, relief=relief, S_m2=S, span_m=b,
               chord_m=cbar, sweep_deg=math.degrees(lam), mfp_mm=mfp(p, T_GAS) * 1e3)
    # surface pressures: upper suction limited to 0.7 p_inf (vacuum limit), rest on the lower surface
    dpu = min(0.5 * WSe, 0.7 * p)
    p_l, p_u = p + (WSe - dpu), p - dpu
    comp = compression(M, T_GAS, p_l / p)
    if comp is None:
        out.update(lift_ok=0, alpha_deg=90.0)
        return out
    T_l, vr, alpha = comp
    V_l = V * vr
    T_u = T_GAS * (p_u / p) ** ((GAM - 1) / GAM)
    V_u = math.sqrt(max(1.0, 2 * CP * (T0 - T_u)))
    out.update(lift_ok=1, alpha_deg=alpha, p_lower_Pa=p_l, T_lower_edge=T_l)
    rt = d["re_tr"]

    def wl(Tw):
        return plate(p_l, T_l, V_l, cbar, Tw, rt)[0]

    def wu(Tw):
        return plate(p_u, T_u, V_u, cbar, Tw, rt)[0]

    def bo(Tw):
        return plate(p, T_GAS, V, pd_["L"], Tw, rt)[0]

    def no(Tw):
        return q_stag(rho, V, T0, d["R_n"], Tw)

    def le(Tw):
        # swept cylinder: sphere of radius R_le / sqrt(2) (velocity gradient), x cos(sweep) (NASA TN, 60 deg delta)
        return q_stag(rho, V, T0, d["R_le"], Tw) / math.sqrt(2) * math.cos(lam)

    L_le = b / math.cos(lam)
    surf = {  # name: (flux fn, area m2, temperature limit)
        "nose": (no, 0.5 * 2 * math.pi * d["R_n"] ** 2, d["T_lim_le"]),
        "le": (le, 0.5 * math.pi * d["R_le"] * L_le, d["T_lim_le"]),
        "wing_lower": (wl, S, d["T_lim_skin"]),
        "wing_upper": (wu, S, d["T_lim_skin"]),
        "body": (bo, pd_["S_wet"], d["T_lim_skin"]),
    }
    Qin = Qcold = Qsink = 0.0
    for k, (fn, A, Tl) in surf.items():
        Teq = rad_eq(fn, T0)
        out[f"T_{k}"] = Teq
        out[f"q_{k}_kWm2"] = fn(Teq) / 1e3
        out[f"q_{k}_cold_kWm2"] = fn(300.0) / 1e3
        Qin += fn(Teq) * A
        Qcold += fn(300.0) * A
        Qsink += q_net(fn, Tl) * A if Teq > Tl else 0.0
    # cabin: insulation leak from the hot body + payload heat; cannot be radiated (skin is hotter than cabin)
    T_body = min(out["T_body"], d["T_lim_skin"])
    leak = d["k_ins"] / d["t_ins"] * max(0.0, T_body - d["T_cab"]) * pd_["S_wet"]
    pay = d["w_seat"] * pd_["n"]
    # drag: flat-plate lift (L tan alpha) + skin friction + wing thickness wave drag + body wave + vacuum base
    F_fric = (plate(p_l, T_l, V_l, cbar, out["T_wing_lower"], rt)[1] + plate(p_u, T_u, V_u, cbar, out["T_wing_upper"], rt)[1]) * S \
        + plate(p, T_GAS, V, pd_["L"], out["T_body"], rt)[1] * pd_["S_wet"]
    D_lift = W * relief * math.tan(math.radians(alpha))
    D_wave = q * (0.002 * S + 0.15 * pd_["A"]) + p * pd_["A"]
    D = D_lift + F_fric + D_wave
    pl = plate(p_l, T_l, V_l, cbar, out["T_wing_lower"], rt)[2]
    out.update(Q_in_kW=Qin / 1e3, Q_in_cold_kW=Qcold / 1e3, Q_excess_kW=Qsink / 1e3, Q_cabin_kW=(leak + pay) / 1e3,
               drag_N=D, drag_fric_N=F_fric, LD=W / D, heat_over_fricpower=Qcold / (F_fric * V),
               Re_chord=pl["Re_L"], chi_lower=pl["chi"], turb_frac_lower=pl["turb_frac"],
               Kn_nose=mfp(p, T_GAS) / d["R_n"], Kn_le=mfp(p, T_GAS) / d["R_le"],
               q_nose_SG_cold_kWm2=K_SG * math.sqrt(rho / d["R_n"]) * V ** 3 / 1e3)
    return out


def trip_time(V, d=BASE):
    """Equivalent full-heat time for a trip: cruise + accel/decel at a_acc (heat ~ V^1.5 while accelerating)."""
    if V * V / d["a_acc"] > d["L_trip"]:
        return float("nan")
    return d["L_trip"] / V - 0.2 * V / d["a_acc"]


def sink_fraction(st, d=BASE):
    t = trip_time(st["V"], d)
    Q = (st["Q_excess_kW"] + st["Q_cabin_kW"]) * 1e3 * t
    return Q / (d["h_sink"] * d["m"]), Q, t


# ----------------------------------------------------------------------------- tube side
def tube_side(st, D_b, d=BASE, N=None):
    p, V = st["p_Pa"], st["V"]
    rho = p / (R_AIR * T_GAS)
    A = math.pi * D_b ** 2 / 4
    drag = st["drag_N"]
    Pi = drag / (p * A)
    dT_pass = drag / (rho * A * CV)
    # gas -> wall: conduction (Nu_D = 8 for uniform heating, mean temperature) or natural convection
    nu, al = mu(T_GAS) / rho, k_gas(T_GAS) / (rho * CP)
    dT = d["dT_gas"]
    Ra = G * dT / T_GAS * D_b ** 3 / (nu * al)
    Nu = max(8.0, 0.1 * Ra ** (1 / 3)) * d["f_mix"]
    h = Nu * k_gas(T_GAS) / D_b
    P_max = math.pi * D_b * h * dT                       # W per m of tube
    E_gas = max(drag - st["Q_in_kW"] * 1e3 / V, 0.3 * drag)  # J/m into the gas (pod-radiated part goes to the wall)
    N_max = 3600 * P_max / E_gas
    tau = rho * A * CV / (math.pi * D_b * h)
    N = d["N_req"] if N is None else N
    P_avg = N / 3600 * drag                               # W/m, everything ends on the wall
    q_wall = P_avg / (math.pi * D_b)
    vdot_km = P_avg * 1e3 / (rho * CP * dT)               # m3/s per km of tube to carry the heat out as gas
    # films: envelope film in sun (alpha/eps), plus traffic flux on the tube film (both faces radiate)
    eps_f = 0.2
    T_film_sun = ((d["film_ae"] * eps_f * 1361 * 1.3 + eps_f * 240) / (2 * eps_f * SIGMA)) ** 0.25
    T_film_tube = (T_ENV ** 4 + q_wall / (2 * eps_f * SIGMA)) ** 0.25
    # guideway / motor strip: LSM loss + pod belly radiation + gap convection, at launch thrust, radiation only
    c_root = 2 * st["chord_m"]
    F_launch = drag + d["m"] * d["a_acc"]
    E_lsm = (1 / d["eta_lsm"] - 1) * F_launch / d["w_strip"]              # J/m2 per pass
    E_rad = 0.5 * EPS * SIGMA * st["T_wing_lower"] ** 4 * c_root / V      # view factor 0.5 (gap ~ m) (A)
    E_conv = st["q_wing_lower_kWm2"] * 1e3 * c_root / V
    q_strip = N / 3600 * (E_lsm + E_rad + E_conv)
    T_strip = (T_ENV ** 4 + q_strip / (2 * EPS * SIGMA)) ** 0.25
    return dict(D_bore=D_b, Pi=Pi, dT_pass=dT_pass, Ra=Ra, Nu=Nu, tau_relax_s=tau, N_max_per_h=N_max,
                q_wall_Wm2=q_wall, vdot_cool_m3s_km=vdot_km, T_film_sun=T_film_sun, T_film_tube=T_film_tube,
                E_strip_pass_Jm2=E_lsm + E_rad + E_conv, T_strip=T_strip)


# ----------------------------------------------------------------------------- checks
CHECKS = ["lift", "span", "trip", "skin_T", "le_T", "heat_sink", "gas_dT_pass", "throughput", "film_T", "guideway_T"]
LABEL = {"lift": "lift (C_L / incidence)", "span": "wing span > bore", "skin_T": "skin temperature",
         "le_T": "nose / LE temperature", "heat_sink": "heat-sink mass", "gas_dT_pass": "gas ΔT per pass",
         "throughput": "gas heat → pods/h", "trip": "trip too short", "film_T": "film temperature", "guideway_T": "guideway temperature",
}


def evaluate(st, D_b, d=BASE):
    r = dict(st)
    mg, tol = {}, {}
    if not st.get("lift_ok"):
        mg = {k: (-9.0 if k == "lift" else 0.0) for k in CHECKS}
        tol = {k: 0.3 for k in CHECKS}
        r.update(D_bore=D_b)
    else:
        mg["lift"] = math.log10(d["alpha_le_ok"] / st["alpha_deg"])
        tol["lift"] = math.log10(d["alpha_le_marg"] / d["alpha_le_ok"])
        mg["trip"] = math.log10(d["L_trip"] * d["a_acc"] / st["V"] ** 2)
        tol["trip"] = 1e-6
        mg["span"] = math.log10(d["span_frac"] * D_b / st["span_m"])
        tol["span"] = math.log10(0.8 / d["span_frac"])
        Tsk = max(st["T_wing_lower"], st["T_body"])
        # above the material limit the surface must be actively cooled: 'marginal', and the heat goes to the sink check
        tol["skin_T"] = math.log10(2)
        mg["skin_T"] = max(4 * math.log10(d["T_lim_skin"] / Tsk), -0.99 * tol["skin_T"])
        Tle = max(st["T_nose"], st["T_le"])
        tol["le_T"] = math.log10(2)
        mg["le_T"] = max(4 * math.log10(d["T_lim_le"] / Tle), -0.99 * tol["le_T"])
        f, Q, t = sink_fraction(st, d)
        r.update(t_trip_s=t, Q_sink_MJ=Q / 1e6, f_sink=f, Q_in_trip_MJ=st["Q_in_kW"] * 1e3 * t / 1e6,
                 Q_in_cold_trip_MJ=st["Q_in_cold_kW"] * 1e3 * t / 1e6)
        if math.isnan(f):                                  # cruise Mach not reachable: flagged by 'trip'
            mg["heat_sink"], tol["heat_sink"] = 0.0, 0.3
        else:
            mg["heat_sink"] = math.log10(d["f_sink_ok"] / max(f, 1e-6))
            tol["heat_sink"] = math.log10(d["f_sink_marg"] / d["f_sink_ok"])
        ts = tube_side(st, D_b, d)
        r.update(ts)
        mg["gas_dT_pass"] = math.log10(d["dT_pass_ok"] / ts["dT_pass"])
        tol["gas_dT_pass"] = math.log10(d["dT_pass_marg"] / d["dT_pass_ok"])
        mg["throughput"] = math.log10(ts["N_max_per_h"] / d["N_req"])
        tol["throughput"] = math.log10(d["N_req"] / d["N_marg"])
        Tf = max(ts["T_film_sun"], ts["T_film_tube"])
        mg["film_T"] = 4 * math.log10(d["T_film_lim"] / Tf)
        tol["film_T"] = 4 * math.log10(d["T_film_marg"] / d["T_film_lim"])
        mg["guideway_T"] = 4 * math.log10(d["T_strip_lim"] / ts["T_strip"])
        tol["guideway_T"] = 4 * math.log10(453.0 / d["T_strip_lim"])
    worst, wm = None, 1e9
    st_all = "pass"
    for k in CHECKS:
        s = "pass" if mg[k] >= 0 else ("marginal" if mg[k] >= -tol[k] else "fail")
        r[f"m_{k}"] = round(mg[k], 4)
        r[f"s_{k}"] = s
        if s == "fail":
            st_all = "fail"
        elif s == "marginal" and st_all == "pass":
            st_all = "marginal"
        # binding = smallest margin in units of its own tolerance
        z = mg[k] / tol[k]
        if z < wm:
            wm, worst = z, k
    r["status"] = st_all
    r["binding"] = worst
    # same, with the span-vs-bore check left out (maps show span separately)
    sn, wn, zn = "pass", None, 1e9
    for k in CHECKS:
        if k == "span":
            continue
        if r[f"s_{k}"] == "fail":
            sn = "fail"
        elif r[f"s_{k}"] == "marginal" and sn == "pass":
            sn = "marginal"
        if mg[k] / tol[k] < zn:
            zn, wn = mg[k] / tol[k], k
    r["status_ns"], r["binding_ns"] = sn, wn
    # heating-only status (flight checks excluded) for the 'heating-feasible band'
    hs = "pass"
    for k in CHECKS[3:]:
        if r[f"s_{k}"] == "fail":
            hs = "fail"
            break
        if r[f"s_{k}"] == "marginal":
            hs = "marginal"
    r["heat_status"] = hs if st.get("lift_ok") else "fail"
    return r


# ----------------------------------------------------------------------------- io
def write_csv(name, rows):
    keys = []
    for r in rows:
        for k in r:
            if k not in keys:
                keys.append(k)
    with open(os.path.join(HERE, name), "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=keys)
        w.writeheader()
        for r in rows:
            w.writerow({k: (f"{v:.5g}" if isinstance(v, float) else v) for k, v in r.items()})
    print(f"wrote {name} ({len(rows)} rows)")


# ----------------------------------------------------------------------------- stages
def stage_grid(d=BASE, ws_list=WS_LIST, d_list=D_LIST, m_grid=M_GRID, p_grid=P_GRID):
    rows = []
    for WS in ws_list:
        for M in m_grid:
            for p in p_grid:
                st = pod_state(float(M), float(p), WS, d)
                for D_b in d_list:
                    rows.append(evaluate(st, D_b, d))
    return rows


def best_ld_p(M, WS, d=BASE):
    """Pressure that maximises L/D at Mach M and wing loading WS (golden search in log p)."""
    from scipy.optimize import minimize_scalar

    def f(lp):
        st = pod_state(M, 10 ** lp, WS, d)
        return -st["LD"] if st.get("lift_ok") else 10.0
    res = minimize_scalar(f, bounds=(-0.5, 4.0), method="bounded", options=dict(xatol=0.01))
    return 10 ** res.x


def stage_bestld(d=BASE):
    rows = []
    for WS in WS_LIST:
        for M in np.arange(3.0, 8.01, 0.5):
            p = best_ld_p(float(M), WS, d)
            st = pod_state(float(M), p, WS, d)
            rows.append(dict(WS_Pa=WS, Mach=M, p_bestLD_Pa=p, q_Pa=st["q_Pa"], CL=st["CL"], alpha_deg=st["alpha_deg"],
                             LD=st["LD"], drag_N=st["drag_N"], T_le=st["T_le"], T_nose=st["T_nose"],
                             T_wing_lower=st["T_wing_lower"], Q_in_kW=st["Q_in_kW"]))
    return rows


def stage_scaling(d=BASE):
    """Heat per trip at fixed q (lift-supported at fixed W/S and C_L) and on the best-L/D line."""
    rows = []
    WS = 150.0
    for label, qfix in (("q 300 Pa (C_L~0.5)", 300.0), ("q 1 kPa (C_L~0.15)", 1000.0), ("q 3 kPa (C_L~0.05)", 3000.0),
                        ("best L/D", None)):
        for M in np.arange(3.0, 8.01, 0.25):
            V = M * C_REF
            p = qfix / (0.7 * M * M) if qfix else best_ld_p(float(M), WS, d)
            st = pod_state(float(M), p, WS, d)
            if not st.get("lift_ok"):
                continue
            t = trip_time(V, d)
            f, Q, _ = sink_fraction(st, d)
            # components of cold-wall input
            comp = {k: st[f"q_{k}_cold_kWm2"] for k in ("nose", "le", "wing_lower", "wing_upper", "body")}
            rows.append(dict(case=label, Mach=M, V=V, p_Pa=p, q_Pa=st["q_Pa"], CL=st["CL"], LD=st["LD"],
                             t_trip_s=t, Q_in_cold_kW=st["Q_in_cold_kW"], Q_in_kW=st["Q_in_kW"],
                             Q_in_cold_trip_GJ=st["Q_in_cold_kW"] * t / 1e6, Q_in_trip_GJ=st["Q_in_kW"] * t / 1e6,
                             Q_sink_trip_MJ=Q / 1e6, f_sink=f, drag_kN=st["drag_N"] / 1e3,
                             drag_energy_trip_GJ=st["drag_N"] * d["L_trip"] / 1e9,
                             heat_over_fricpower=st["heat_over_fricpower"],
                             T_nose=st["T_nose"], T_le=st["T_le"], T_wing_lower=st["T_wing_lower"],
                             T_wing_upper=st["T_wing_upper"], T_body=st["T_body"],
                             **{f"qcold_{k}_kWm2": v for k, v in comp.items()}))
    # power-law exponents of trip heat vs V
    for label in sorted(set(r["case"] for r in rows)):
        sub = [r for r in rows if r["case"] == label and not math.isnan(r["t_trip_s"])]
        lv = np.log([r["V"] for r in sub])
        for key in ("Q_in_cold_kW", "Q_in_cold_trip_GJ", "Q_in_trip_GJ"):
            n = np.polyfit(lv, np.log([r[key] for r in sub]), 1)[0]
            for r in sub:
                r[f"exp_{key}"] = n
    return rows


def stage_radius(d=BASE):
    rows = []
    WS = 150.0
    for M in (3.0, 4.0, 5.0, 6.0, 7.0, 8.0):
        for qfix in (300.0, 1000.0, 3000.0):
            p = qfix / (0.7 * M * M)
            V, rho, q, T0 = freestream(M, p)
            for Rr in (0.005, 0.01, 0.02, 0.03, 0.05, 0.1, 0.2, 0.3, 0.5, 1.0):
                fn = lambda Tw: q_stag(rho, V, T0, Rr, Tw)
                fl = lambda Tw: q_stag(rho, V, T0, Rr, Tw) / math.sqrt(2) * math.cos(math.atan(4 / d["AR"]))
                Tn, Tl = rad_eq(fn, T0), rad_eq(fl, T0)
                qc = K_SG * math.sqrt(rho / Rr) * V ** 3
                rows.append(dict(Mach=M, q_Pa=qfix, p_Pa=p, T0=T0, R_m=Rr, Kn=mfp(p, T_GAS) / Rr,
                                 Re_R=rho * V * Rr / mu(T_GAS), q_SG_cold_kWm2=qc / 1e3,
                                 q_FM_kWm2=0.5 * rho * V ** 3 / 1e3, T_nose_eq=Tn, T_le_eq=Tl,
                                 FR_over_SG=q_fay_riddell(M, p, Rr, 300.0) / qc,
                                 sink_kWm2_at_1700K=q_net(fn, 1700.0) / 1e3))
    return rows


MATERIALS = {
    "Al / CFRP skin (420 K) + C/SiC LE (1,700 K)": dict(T_lim_skin=420.0, T_lim_le=1700.0),
    "Ti hot structure (700 K) + C/SiC LE (1,700 K)": dict(T_lim_skin=700.0, T_lim_le=1700.0),
    "Inconel TPS (1,250 K) + C/SiC LE (1,700 K) [base]": dict(),
    "C/SiC TPS (1,600 K) + UHTC LE (2,200 K)": dict(T_lim_skin=1600.0, T_lim_le=2200.0),
}


def stage_sink():
    rows = []
    WS = 150.0
    for mname, mod in MATERIALS.items():
        for L in (500e3, 1000e3, 4000e3):
            d = dict(BASE, L_trip=L, **mod)
            for M in np.arange(3.0, 8.01, 0.5):
                for p in (3.0, 10.0, 30.0, 100.0, 300.0):
                    st = pod_state(float(M), p, WS, d)
                    if not st.get("lift_ok"):
                        continue
                    f, Q, t = sink_fraction(st, d)
                    rows.append(dict(materials=mname, L_trip_km=L / 1e3, Mach=M, p_Pa=p, CL=st["CL"],
                                     alpha_deg=st["alpha_deg"], t_trip_s=t, Q_excess_kW=st["Q_excess_kW"],
                                     Q_cabin_kW=st["Q_cabin_kW"], Q_sink_MJ=Q / 1e6, f_sink_water=f,
                                     f_sink_PCM=f * d["h_sink"] / 0.25e6, f_sink_LH2sys=f * d["h_sink"] / 2.5e6 * 1.0,
                                     T_skin=max(st["T_wing_lower"], st["T_body"]), T_le=max(st["T_le"], st["T_nose"])))
    return rows


def stage_tube():
    rows = []
    for WS in (50.0, 150.0, 500.0):
        for M in (3.0, 4.0, 5.0, 6.0, 8.0):
            for p in (1.0, 3.0, 5.5, 10.0, 30.0, 100.0, 300.0, 1000.0):
                st = pod_state(M, p, WS)
                if not st.get("lift_ok"):
                    continue
                for D_b in D_LIST:
                    ts = tube_side(st, D_b)
                    P_max = ts["N_max_per_h"] * st["drag_N"] / 3600
                    rows.append(dict(WS_Pa=WS, Mach=M, p_Pa=p, drag_kN=st["drag_N"] / 1e3, LD=st["LD"],
                                     drag_MJ_per_km=st["drag_N"] / 1e3, **ts,
                                     seats_per_h_3p5t=ts["N_max_per_h"] * R.pod(BASE["m"])["n"],
                                     seats_per_h_300kg=3600 * P_max / (300 * G / st["LD"])))
    return rows


LEVERS = [
    ("base: W/S 150 Pa, bore 80 m, R_n 0.1 m, R_LE 3 cm, Inconel 1,250 K / C-SiC 1,700 K", {}, 150.0, 80.0),
    ("sharper nose and LE: R_n 3 cm, R_LE 1 cm", dict(R_n=0.03, R_le=0.01), 150.0, 80.0),
    ("blunter: R_n 0.3 m, R_LE 10 cm", dict(R_n=0.3, R_le=0.10), 150.0, 80.0),
    ("lower wing loading 50 Pa (bigger wings)", {}, 50.0, 80.0),
    ("higher wing loading 500 Pa", {}, 500.0, 80.0),
    ("bore 40 m (W/S 300 to fit span)", {}, 300.0, 40.0),
    ("bore 20 m (W/S 500 to fit span)", {}, 500.0, 20.0),
    ("hotter materials: C/SiC TPS 1,600 K + UHTC LE 2,200 K", dict(T_lim_skin=1600.0, T_lim_le=2200.0), 150.0, 80.0),
    ("Ti hot structure 700 K", dict(T_lim_skin=700.0), 150.0, 80.0),
    ("Al / CFRP skin 420 K (sink takes the rest)", dict(T_lim_skin=420.0), 150.0, 80.0),
    ("Al / CFRP 420 K + 15 % heat sink allowed", dict(T_lim_skin=420.0, f_sink_ok=0.15, f_sink_marg=0.30), 150.0, 80.0),
    ("gas stirring x3 and +150 K allowed gas rise", dict(f_mix=3.0, dT_gas=150.0), 150.0, 80.0),
    ("throughput 30 pods/h instead of 120", dict(N_req=30.0, N_marg=8.0), 150.0, 80.0),
    ("17 t, 50-seat pod (W/S 300 to fit span)", dict(m=17000.0), 300.0, 80.0),
    ("1,000 km trip", dict(L_trip=1000e3), 150.0, 80.0),
]


def window(rows, M, key="status", ok=("pass",)):
    """Pressure window at Mach M where rows[key] is in ok, and the binding check just outside each edge."""
    sub = sorted([r for r in rows if abs(r["Mach"] - M) < 1e-6], key=lambda r: r["p_Pa"])
    good = [k for k, r in enumerate(sub) if r[key] in ok]
    if not good:
        best = max(sub, key=lambda r: min(r[f"m_{c}"] / 1.0 for c in CHECKS))
        return "none", best["binding"], ""
    lo, hi = good[0], good[-1]
    b_lo = sub[lo - 1]["binding"] if lo > 0 else "(grid edge)"
    b_hi = sub[hi + 1]["binding"] if hi < len(sub) - 1 else "(grid edge)"
    return f"{sub[lo]['p_Pa']:.3g}-{sub[hi]['p_Pa']:.3g}", b_lo, b_hi


def band(rows, key="status"):
    out = {}
    for lvl, ok in (("pass", ("pass",)), ("marg", ("pass", "marginal"))):
        Ms = sorted(set(r["Mach"] for r in rows if r[key] in ok))
        out[f"{key}_M_{lvl}"] = f"{Ms[0]:.1f}-{Ms[-1]:.1f}" if Ms else "none"
    return out


def lever_row(g, extra):
    r = dict(extra)
    r.update(band(g))
    r.update(band(g, "heat_status"))
    for M in (3.0, 4.0, 6.0, 8.0):
        w, blo, bhi = window(g, M)
        r[f"M{M:g}_p_pass"], r[f"M{M:g}_floor"], r[f"M{M:g}_ceiling"] = w, blo, bhi
        w, blo, bhi = window(g, M, ok=("pass", "marginal"))
        r[f"M{M:g}_p_marg"] = w
        r[f"M{M:g}_heat_p_pass"] = window(g, M, "heat_status")[0]
    return r


def stage_levers():
    rows = []
    mg = np.array([3.0, 4.0, 6.0, 8.0])
    for name, mod, WS, D_b in LEVERS:
        d = dict(BASE, **mod)
        g = stage_grid(d, [WS], [D_b], mg, P_GRID)
        bl = pod_state(6.0, best_ld_p(6.0, WS, d), WS, d)
        bl8 = pod_state(8.0, best_ld_p(8.0, WS, d), WS, d)
        e = evaluate(bl, D_b, d)
        r = lever_row(g, dict(lever=name, WS_Pa=WS, D_bore=D_b, M6_bestLD_p=bl["p_Pa"],
                                  M8_T_le_bestLD=max(bl8["T_le"], bl8["T_nose"]), M8_T_skin_bestLD=bl8["T_wing_lower"],
                                  M6_f_sink_bestLD=e["f_sink"], M6_Nmax_bestLD=e["N_max_per_h"]))
        rows.append(r)
        print(f"  lever: {name}: M6 pass {r['M6_p_pass']} ({r['M6_floor']} / {r['M6_ceiling']}), "
              f"M8 pass {r['M8_p_pass']} ({r['M8_floor']} / {r['M8_ceiling']})")
    return rows


def stage_summary(grid):
    rows = []
    for WS in WS_LIST:
        for D_b in D_LIST:
            sub = [r for r in grid if r["WS_Pa"] == WS and r["D_bore"] == D_b]
            rows.append(lever_row(sub, dict(WS_Pa=WS, D_bore=D_b)))
    return rows


def stage_check():
    """The coordinator's numbers: 3.5 t, ~100 m2 wing, at 5.5 and 50 Pa."""
    rows = []
    for WS in (343.0, 150.0, 50.0):
        for p in (5.5, 50.0):
            for M in (3.0, 6.0, 8.0):
                st = pod_state(M, p, WS)
                r = dict(WS_Pa=WS, S_m2=BASE["m"] * G / WS, p_Pa=p, Mach=M, q_Pa=st["q_Pa"], CL_req=st["CL"],
                         lift_ok=st.get("lift_ok"), relief=st["relief"])
                if st.get("lift_ok"):
                    r.update(alpha_deg=st["alpha_deg"], LD=st["LD"], drag_kN=st["drag_N"] / 1e3,
                             T_nose=st["T_nose"], T_le=st["T_le"], T_skin=st["T_wing_lower"])
                    for D_b in (20.0, 80.0):
                        ts = tube_side(st, D_b)
                        r[f"dT_pass_D{D_b:g}"] = ts["dT_pass"]
                        r[f"Nmax_D{D_b:g}"] = ts["N_max_per_h"]
                rows.append(r)
    return rows


def main():
    args = sys.argv[1:]
    stages = args or ["grid", "bestld", "scaling", "radius", "sink", "tube", "check", "levers"]
    if "grid" in stages:
        g = stage_grid()
        keep = ["WS_Pa", "D_bore", "Mach", "p_Pa", "q_Pa", "CL", "alpha_deg", "LD", "drag_N", "span_m", "T0",
                "T_nose", "T_le", "T_wing_lower", "T_wing_upper", "T_body", "Q_in_kW", "Q_excess_kW", "Q_cabin_kW",
                "f_sink", "Kn_le", "Re_chord", "Pi", "dT_pass", "tau_relax_s", "N_max_per_h", "T_film_sun",
                "T_film_tube", "T_strip"] + [f"m_{k}" for k in CHECKS] + \
               ["status", "binding", "status_ns", "binding_ns", "heat_status"]
        write_csv("heating_grid.csv", [{k: r[k] for k in keep if k in r} for r in g])
        write_csv("heating_summary.csv", stage_summary(g))
    if "bestld" in stages:
        write_csv("heating_bestld.csv", stage_bestld())
    if "scaling" in stages:
        write_csv("heating_scaling.csv", stage_scaling())
    if "radius" in stages:
        write_csv("heating_radius.csv", stage_radius())
    if "sink" in stages:
        write_csv("heating_sink.csv", stage_sink())
    if "tube" in stages:
        write_csv("heating_tube.csv", stage_tube())
    if "check" in stages:
        write_csv("heating_check.csv", stage_check())
    if "levers" in stages:
        write_csv("heating_levers.csv", stage_levers())


if __name__ == "__main__":
    main()
