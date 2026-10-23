#!/usr/bin/env python3
"""Skyloop operating-regime map (screening level).

Maps pass / marginal / fail per discipline over tube pressure (0.1 Pa - 1 kPa), Mach (1-8),
pod size (2-30 t), blockage ratio beta = (d_pod/D)^2, support (maglev vs aerodynamic glide),
tube setting (floating at 20 km / elevated on pylons / buried) and headway rule.

Run:   python3 skyloop-review/regime/regime.py          (~1-2 min, one process)
Out:   regime_grid.csv     reference pod (50 seats, 15 t), every setting x support x Mach x pressure,
                           best beta per point, margins + status per discipline, both headway rules
       regime_pods.csv     pod mass x Mach (best p and beta), per setting x support x rule
       regime_summary.csv  feasible ranges per setting x support x rule
       glide_check.csv     the coordinator's 5.5 Pa glide check, re-derived
       regime_levers.csv   Mach limits per binding discipline under levers (suspension, survey, tethers, R)
       then make_figures.py -> regime_maps.png, regime_pods.png, regime_glide.png
       (python3 regime.py --levers runs only the lever stage)

Reuses read-only: physics/skyloop_physics.py (atmosphere, wall_sandwich, envelope, kantrowitz_beta,
mfp, C_REF) and system/sky_model.py (cap_pair_yr, headway_s, ride_min, run_cost, P).
Mach convention (repo): V = M x 295.1 m/s whatever the setting.  (A) = assumption.
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
sys.path.insert(0, os.path.join(REV, "physics"))
sys.path.insert(0, os.path.join(REV, "system"))
import skyloop_physics as sp  # noqa: E402  (read-only import; its main() is not run)
import sky_model as sm  # noqa: E402  (read-only import)

G = 9.81
GAM, R_AIR, SIGMA = 1.4, 287.05, 5.670e-8
C_REF = sp.C_REF                        # 295.1 m/s
CP = 1005.0

# ----------------------------------------------------------------------------- grids
P_GRID = np.geomspace(0.1, 1000.0, 41)
M_GRID = np.round(np.arange(1.0, 8.01, 0.2), 2)
BETAS = [0.03, 0.05, 0.08, 0.12, 0.16, 0.22, 0.30]
D_MAX = 8.0                             # m, largest practical bore (Channel Tunnel running bores 7.6 m) (A)
RULES = ["brick-wall", "moving block 10 s"]
SUPPORTS = ["maglev", "glide"]
DEMAND = 81e6                           # spine pax/yr (BoomShield coverage; system capacity table)

# ----------------------------------------------------------------------------- pod family (A)
KG_PER_SEAT = 300.0                     # loaded, incl. 2 t base below; system.md used 500 (conservative)
LF = sm.P["lf"]                         # 0.8


def pod(m_kg):
    """Pressurised pod sized from mass: seats, abreast, diameter, length, areas."""
    n = max(2, int(round((m_kg - 2000.0) / KG_PER_SEAT)) if m_kg > 2600 else 2)
    k = 1 if n <= 6 else 2 if n <= 30 else 3 if n <= 60 else 4
    d = 0.55 * k + 0.6                  # m: seats 0.55 m + aisle/walls (A)
    L = math.ceil(n / k) * 0.85 + 7.0   # m: 0.85 m pitch + nose/tail/lock/equipment (A)
    return dict(m=m_kg, n=n, k=k, d=d, L=L, A=math.pi * d * d / 4, S_wet=0.9 * math.pi * d * L)


REF_MASS = 2000.0 + 50 * KG_PER_SEAT    # 17 t, 50 seats
# ----------------------------------------------------------------------------- settings (A unless noted)
A20 = sp.atm(20.0)
SETTINGS = {
    "floating": dict(label="Floating, 20 km (H2-lifted, tethered)", p_amb=float(A20.p), rho_amb=float(A20.rho),
                     T_amb=float(A20.T), T_gas=250.0,
                     A_float=25.9, lam_float=6000.0,      # physics: 1 m/s gust, 10 MN tension, 6 km wavelength
                     R_avail=1e12, leak_cm2=1.0, outgas=1e-5, evac="marginal"),
    "elevated": dict(label="Elevated on pylons, near ground", p_amb=100.1e3, rho_amb=1.213, T_amb=287.4, T_gas=288.0,
                     A0=3e-3, slope=1.5e-5,               # sun bow between 30 m pylons; 15 mm/km long-wave after realignment
                     R_avail=150e3, leak_cm2=0.1, outgas=1e-7, evac="pass"),
    "buried": dict(label="Buried / deep tunnel", p_amb=101.3e3, rho_amb=1.225, T_amb=288.0, T_gas=288.0,
                   A0=1e-3, slope=5e-6,                   # thermally stable; settlement/seismic only
                   R_avail=1e12, leak_cm2=0.1, outgas=1e-7, evac="pass"),
}

# ride / alignment
A_LIM = 0.05 * G                        # physics.md ride limit
F_S = 1.0                               # Hz, secondary suspension (passive); active 0.5 Hz is a lever
# energy
ETA_D = ETA_R = 0.85                    # system.md
L_TRIP = 1500e3                         # m, representative passenger trip (A)
LD_MAG = 300.0                          # EDS lift/drag (system.md)
HOTEL_W_SEAT = 400.0                    # system.md
E_LIM, E_JET = 100.0, 340.0             # Wh/pax-km: keep a clear advantage / jet-fuel parity
# glide
CL_MAX = 1.0                            # optimistic ram / ground-effect lift coefficient (coordinator's value)
GAP = 0.05                              # m, cushion gap for Knudsen check
# heating
T_SKIN, EPS = 400.0, 0.85
STORE_FRAC, STORE_J_KG = 0.03, 4.0e5    # 3 % of pod mass as ice/PCM heat sink
R_NOSE, T_NOSE_MAX = 0.10, 1900.0
# vacuum
S_LIM = 20.0                            # m3/s per km installed (~200 kW/km, ~5 t of pumps per km)
# structure
R_ENV_LIM = 100.0
WALL_LIM = 6000.0                       # kg/m, steel ring-stiffened tube on pylons
# safety
A_PLUG_OK, A_PLUG_MAX = 5.0, 20.0       # g, meeting ambient air in a full-bore breach
# economics
PAXKM_PER_KM = 40e6                     # spine all-stop: 146 B pax-km / 3,686 km (system networks.csv)
USD_KWH = 0.07 + 0.15 * 190e-3          # electricity + 2045 grid CO2 at $190/t
PV = (1 - 1.07 ** -25) / 0.07           # 11.65
LSM_EFF = {"low": 8.0 * 1.15 * 1.15, "base": 15.0 * 1.25 * 1.30}       # $M/km incl soft + contingency
GUIDE_EFF = {"low": 3.0 * 1.15 * 1.15, "base": 6.0 * 1.25 * 1.30}
STOP_SPACING = 720e3                    # spine all-stop mean leg (3,600 km / 5)


# ----------------------------------------------------------------------------- helpers
def mu(T):
    return 1.458e-6 * T ** 1.5 / (T + 110.4)


def cf(Re):
    Re = max(Re, 50.0)
    return 1.328 / math.sqrt(Re) if Re < 1e6 else 0.074 * Re ** -0.2


def stanton(Re):
    return min(0.5 * cf(Re) * 0.72 ** (-2 / 3), 0.05)


def beta_kantrowitz(Mg):
    if Mg <= 1.0:   # subsonic: choke when annulus reaches sonic (isentropic)
        return max(0.0, 1 - 1 / sp.area_ratio_isentropic(max(Mg, 0.05)))
    return sp.kantrowitz_beta(Mg)[0]


def piston_dp(p, up, c):
    """Pressure rise ahead of a piston moving at up into gas at rest (normal-shock relation)."""
    x = up / c * (GAM + 1) / 2
    Ms = (x + math.sqrt(x * x + 4)) / 2
    return p * 2 * GAM / (GAM + 1) * (Ms * Ms - 1), Ms * c


def q_hole(p0, T0, T_gas):
    """Choked leak through 1 cm2 (C_d 0.6), Pa m3/s at tube temperature."""
    mdot = 0.6e-4 * p0 * math.sqrt(GAM / (R_AIR * T0)) * (2 / (GAM + 1)) ** ((GAM + 1) / (2 * (GAM - 1)))
    return mdot * R_AIR * T_gas


# floating envelope radius: precompute R(payload) once with the physics model (slow iteration)
_PAY = np.array([100, 150, 200, 250, 300, 400, 500, 700, 1000], float)
_RENV = None


def r_env(payload):
    global _RENV
    if _RENV is None:
        _RENV = np.array([sp.envelope(x, 20.0)["R"] for x in _PAY])
    return float(np.interp(payload, _PAY, _RENV))


# breakeven average tube cost, $M/km (system breakeven.csv: spine+LAX all-stop, net vs subsonic, 2050, 100 %)
_BE = {2.0: 48.70, 3.0: 55.14, 4.0: 53.84, 6.0: 48.19}
_FLOAT_COST = {2.0: 98.35, 3.0: 101.9, 4.0: 127.6, 6.0: 157.0}   # same file: base avg tube (run + ramps)


def saved_h(M):
    """NYC-LAX door-to-door saving vs subsonic (7.29 h), 100 min portal (system d2d convention)."""
    return 7.291 - (100 / 60 + float(sm.ride_min(3974.4 * sm.P["path"], M)) / 60)


_k = np.polyfit(list(_BE), [v / saved_h(m) for m, v in _BE.items()], 1)


def breakeven(M):
    if 2.0 <= M <= 6.0:
        return float(np.interp(M, list(_BE), list(_BE.values())))
    return float(np.polyval(_k, M) * saved_h(M))


def float_cost(M):
    return float(np.interp(M, [1.0] + list(_FLOAT_COST), [90.0] + list(_FLOAT_COST.values())))


# ----------------------------------------------------------------------------- core evaluation
DISC = ["support", "drag_energy", "heating", "vacuum", "structure", "alignment", "route", "safety",
        "capacity", "economics"]
TECH = DISC[:8]
OPS = DISC[:9]                          # technical + capacity


def status(m, tol):
    return "pass" if m >= 0 else ("marginal" if m >= -tol else "fail")


def evaluate(set_key, support, M, p, m_pod, beta, rule):
    s = SETTINGS[set_key]
    pd_ = pod(m_pod)
    T = s["T_gas"]
    c = math.sqrt(GAM * R_AIR * T)
    V = M * C_REF
    Mg = V / c
    rho = p / (R_AIR * T)
    q = 0.5 * rho * V * V
    D = pd_["d"] / math.sqrt(beta)
    A_t = math.pi * D * D / 4
    W = m_pod * G
    pax = pd_["n"] * LF
    h = sm.headway_s(M, rule)
    out = dict(setting=set_key, support=support, Mach=M, V=V, p_Pa=p, p_frac_ambient=p / s["p_amb"],
               pod_t=m_pod / 1e3, seats=pd_["n"], d_pod=pd_["d"], beta=beta, D_tube=D, rule=rule, q_Pa=q)
    mg, tol = {}, {}

    # 1. support -------------------------------------------------------------------------------
    S_plan = 0.85 * D * pd_["L"]
    if support == "glide":
        CL_req = W / max(q * S_plan, 1e-9)
        Kn = sp.mfp(p, T) / GAP
        LD = 15.0 if Mg < 1 else 4 * (Mg + 3) / Mg          # subsonic WIG ~15; Kuchemann supersonic
        D_sup = W / LD
        mg["support"] = min(math.log10(CL_MAX / CL_req), math.log10(0.1 / Kn))
        tol["support"] = math.log10(1.5)
        out.update(CL_req=CL_req, LD=LD, S_plan=S_plan, wing_loading=m_pod / S_plan, Kn_gap=Kn)
    else:
        D_sup = W / LD_MAG
        mg["support"], tol["support"] = 1.0, 0.0
    out["LSM_Hz_tau1m"] = V / 2.0

    # 2. drag, choking, self-pumping, energy -------------------------------------------------
    Re = rho * V * pd_["L"] / mu(T)
    CD = 0.15 + cf(Re) * pd_["S_wet"] / pd_["A"]
    bK = beta_kantrowitz(Mg)
    choked = beta > bK
    D_par = q * CD * pd_["A"]
    if choked:
        dp_pist, Vs = piston_dp(p, V, c)
        D_par = max(D_par, 1.8 * q * pd_["A"], dp_pist * A_t * 0.5)   # (A) half the full-bore piston
        out.update(piston_dp_Pa=dp_pist, precursor_shock_m_s=Vs, catch_pod_ahead_s=V * h / max(Vs - V, 1.0))
    D_tot = D_par + D_sup
    KE = 0.5 * m_pod * V * V
    trac = KE * (1 / ETA_D - ETA_D * ETA_R)
    hotel = HOTEL_W_SEAT * pd_["n"] * (L_TRIP / V)
    E = (trac + D_tot * L_TRIP / ETA_D + hotel) / (pax * L_TRIP / 1e3) / 3600      # Wh/pax-km
    L_seg = p * A_t * V * h / D_tot / 1e3                                          # km between crossovers
    m_E = math.log10(E_LIM / E)
    m_seg = math.log10(L_seg / 1.0)       # pass >= 1 km, marginal >= 0.25 km (Channel Tunnel relief ducts every 250 m)
    mg["drag_energy"] = min(m_E, m_seg)        # choking itself is allowed; its cost shows in drag and crossovers
    tol["drag_energy"] = math.log10(4.0) if m_seg < m_E else math.log10(E_JET / E_LIM)
    out.update(Re_pod=Re, CD=CD, beta_kantrowitz=bK, choked=int(choked), drag_par_N=D_par, drag_support_N=D_sup,
               power_MW=D_tot * V / ETA_D / 1e6, E_Wh_paxkm=E, crossover_km=L_seg,
               needs_continuous_LSM=int(support == "glide" or D_tot / m_pod * (L_TRIP / V) > 0.1 * V))

    # 3. heating -------------------------------------------------------------------------------
    T_aw = T * (1 + 0.85 * 0.2 * Mg * Mg)
    Q_body = stanton(Re) * rho * V * CP * max(0.0, T_aw - T_SKIN) * pd_["S_wet"]
    Q_cush = 0.0
    if support == "glide":
        p_c = p + W / S_plan
        rho_c = p_c / (R_AIR * T_aw)
        u_c = 0.5 * V
        Re_c = rho_c * u_c * pd_["L"] / mu(T_aw)
        Q_cush = stanton(Re_c) * rho_c * u_c * CP * max(0.0, T_aw - T_SKIN) * S_plan
    Q_in = Q_body + Q_cush + 150.0 * pd_["n"]
    t_trip = L_TRIP / V + V / (0.3 * G)
    Q_out = EPS * SIGMA * (T_SKIN ** 4 - T ** 4) * pd_["S_wet"] + STORE_FRAC * m_pod * STORE_J_KG / t_trip
    q_n = min(1.7415e-4 * math.sqrt(rho / R_NOSE) * V ** 3, 0.5 * rho * V ** 3)
    T_n = (q_n / (0.85 * SIGMA) + T ** 4) ** 0.25
    mg["heating"] = min(math.log10(Q_out / Q_in), 4 * math.log10(T_NOSE_MAX / T_n))
    tol["heating"] = math.log10(2)
    out.update(T_aw=T_aw, Q_in_kW=Q_in / 1e3, Q_out_kW=Q_out / 1e3, T_nose=T_n)

    # 4. vacuum --------------------------------------------------------------------------------
    Q_gas = s["leak_cm2"] * q_hole(s["p_amb"], s["T_amb"], T) + s["outgas"] * math.pi * D * 1000
    S_req = Q_gas / p
    mg["vacuum"], tol["vacuum"] = math.log10(S_LIM / S_req), math.log10(5)
    S_inst = max(S_req, 10.0)
    out.update(S_req_m3s_km=S_req, pump_kW_km=10 * S_req,
               reevac_h_10km=A_t * 1000 / S_inst * math.log(s["p_amb"] / p) / 3600)

    # 5. structure -----------------------------------------------------------------------------
    dp = s["p_amb"] - p
    if set_key == "floating":
        m2, _, _ = sp.wall_sandwich(dp, D / 2)
        wall = (m2 + sp.LINER) * math.pi * D * sp.JOINTS
        frac = 1.0 if support == "glide" else min(1.0, 2 * (V * V / (0.6 * G)) / STOP_SPACING)
        guide = 10.0 if support == "glide" else 60.0
        payload = wall + guide + 150.0 * frac + 15.0 + m_pod / (V * h)
        R_e = r_env(payload)
        mg["structure"], tol["structure"] = math.log10(R_ENV_LIM / R_e), math.log10(1.5)
        out.update(wall_kg_m=wall, payload_kg_m=payload, R_env_m=R_e)
    elif set_key == "elevated":
        m_mono, t_mono = sp.wall_monocoque(dp, D / 2, E=200e9, rho=7850.0)
        wall = 0.6 * m_mono * math.pi * D            # ring-stiffened steel ~60 % of monocoque (A)
        mg["structure"], tol["structure"] = math.log10(WALL_LIM / wall), math.log10(2)
        out.update(wall_kg_m=wall)
    else:
        mg["structure"], tol["structure"] = 0.5, 0.0
        out.update(wall_kg_m=float("nan"))

    # 6. alignment / ride ----------------------------------------------------------------------
    A_allow_c = A_LIM / (2 * math.pi * F_S) ** 2
    if set_key == "floating":
        lam = s["lam_float"]
        f = V / lam
        A_allow = A_LIM / (2 * math.pi * min(f, F_S)) ** 2
        A_ach = s["A_float"]
    else:
        lam = V / F_S
        A_allow = A_allow_c
        A_ach = s["A0"] + s["slope"] * lam
    mg["alignment"], tol["alignment"] = math.log10(A_allow / A_ach), math.log10(2)
    out.update(align_lambda_m=lam, A_allow_m=A_allow, A_achievable_m=A_ach)

    # 7. route geometry ------------------------------------------------------------------------
    R_min = V * V / (0.5 * G)
    L_acc = V * V / (2 * 0.3 * G)
    mg["route"], tol["route"] = min(math.log10(s["R_avail"] / R_min), 3.0), math.log10(2)
    out.update(R_min_km=R_min / 1e3, L_acc_km=L_acc / 1e3, fits_720km_leg=int(2 * L_acc <= STOP_SPACING))

    # 8. capacity ------------------------------------------------------------------------------
    cap = sm.cap_pair_yr(M, rule, seats=pd_["n"])[0]
    mg["capacity"], tol["capacity"] = math.log10(cap / DEMAND), math.log10(3)
    n_pairs = max(1.0, DEMAND / cap)
    out.update(headway_s=h, cap_pair_Mpax_yr=cap / 1e6, tube_pairs_needed=n_pairs)

    # 9. safety --------------------------------------------------------------------------------
    a_plug = 0.5 * (0.5 * s["rho_amb"] * V * V) * pd_["A"] / m_pod / G
    m_s = math.log10(A_PLUG_OK / a_plug)
    if s["evac"] == "marginal":
        m_s = min(m_s, -0.05)
    mg["safety"], tol["safety"] = m_s, math.log10(A_PLUG_MAX / A_PLUG_OK)
    t_stop = V / (0.5 * G)
    out.update(a_plug_g=a_plug, KE_GJ=KE / 1e9, KE_tTNT=KE / 4.184e9, pods_within_stop=t_stop / h)

    # 10. economics hook -----------------------------------------------------------------------
    BE = breakeven(M)
    lsm_frac = 1.0 if support == "glide" else min(1.0, 2 * L_acc / STOP_SPACING)
    E_ref = (trac + (W / LD_MAG + 0.4 * 3.0 / (R_AIR * T) * V * V * pd_["A"] / 2) * L_TRIP / ETA_D + hotel) \
        / (pax * L_TRIP / 1e3) / 3600
    c_energy = max(0.0, E - E_ref) / 1e3 * PAXKM_PER_KM * USD_KWH * PV / 1e6 \
        + 10 * S_req * 8760 * 0.07 * PV / 1e6
    cost = {}
    for lvl in ("low", "base"):
        if set_key == "floating":
            base = float_cost(M) * (0.28 if lvl == "low" else 1.0) * (0.8 + 0.2 * (D / 4) ** 2)
            extra = (LSM_EFF[lvl] - GUIDE_EFF[lvl]) if support == "glide" else 0.0
        else:
            c0 = {"elevated": (34.0, 62.0), "buried": (100.0, 207.0)}[set_key][0 if lvl == "low" else 1]
            sc = (0.5 + 0.5 * (D / 4) ** 2) if set_key == "elevated" else (0.4 + 0.6 * (D / 4) ** 1.5)
            base = c0 * sc
            extra = lsm_frac * LSM_EFF[lvl] - (GUIDE_EFF[lvl] if support == "glide" else 0.0)
        cost[lvl] = (base + extra) * n_pairs + c_energy
    bc_base, bc_low = BE / cost["base"], BE / cost["low"]
    mg["economics"] = math.log10(bc_base)
    tol["economics"] = math.log10(bc_low / bc_base) if bc_low >= 1 else -1e-9
    out.update(breakeven_usd_m_km=BE, cost_base_usd_m_km=cost["base"], cost_low_usd_m_km=cost["low"],
               BC_base=bc_base, BC_low=bc_low)

    # roll-up ----------------------------------------------------------------------------------
    for d in DISC:
        out[f"m_{d}"] = mg[d]
        out[f"s_{d}"] = status(mg[d], tol[d])
    for lab, ds in (("tech", TECH), ("ops", OPS), ("full", DISC)):
        sts = [out[f"s_{d}"] for d in ds]
        out[f"{lab}_status"] = "fail" if "fail" in sts else ("marginal" if "marginal" in sts else "pass")
        out[f"{lab}_nfail"] = sts.count("fail")
        # binding: worst normalised margin (fails ranked by how far beyond the marginal band)
        score = {d: (mg[d] if mg[d] >= 0 else mg[d] - (0 if mg[d] >= -tol[d] else 10)) for d in ds}
        bind = min(score, key=score.get)
        if lab != "tech" and out["tech_status"] == "fail":
            bind = out["tech_binding"]           # report the technical killer first
        elif lab == "full" and out["ops_status"] == "fail":
            bind = out["ops_binding"]
        out[f"{lab}_binding"] = bind
        out[f"{lab}_score"] = min(score.values())
    return out


def best_over(set_key, support, M, p, m_pod, rule, key="full"):
    """Best blockage ratio for this point (max of the worst margin)."""
    best = None
    for b in BETAS:
        if pod(m_pod)["d"] / math.sqrt(b) > D_MAX:
            continue
        r = evaluate(set_key, support, M, p, m_pod, b, rule)
        rank = ({"pass": 2, "marginal": 1, "fail": 0}[r[f"{key}_status"]], -r[f"{key}_nfail"], r[f"{key}_score"])
        if best is None or rank > best[0]:
            best = (rank, r)
    return best[1]


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


# ----------------------------------------------------------------------------- stages
def stage_glide_check():
    """Coordinator's check: 5.5 Pa, 3.5 t pod, C_L 1."""
    rows = []
    for T, lab in [(216.65, "ambient T 217 K"), (250.0, "tube gas 250 K (physics)")]:
        for p in [5.5, 15, 40]:
            for M in [2, 4, 6]:
                V = M * C_REF
                q = 0.5 * p / (R_AIR * T) * V * V
                S = 3500 * G / q
                Mg = V / math.sqrt(GAM * R_AIR * T)
                LD = 4 * (Mg + 3) / Mg
                for ld in (LD, 2.5):
                    P = 3500 * G * V / ld
                    rows.append(dict(T_gas=lab, p_Pa=p, p_frac_20km=p / float(A20.p), Mach=M, q_Pa=q,
                                     area_m2_CL1_3p5t=S, wing_loading_kg_m2_at_38m2=3500 / 38.25,
                                     q_needed_for_38m2=3500 * G / 38.25, LD=ld, power_MW=P / 1e6,
                                     Wh_paxkm_4pax=P / V / 4 / 3.6, Wh_paxkm_20pax=P / V / 20 / 3.6))
    write_csv("glide_check.csv", rows)


def stage_grid():
    rows = []
    for sk in SETTINGS:
        for sup in SUPPORTS:
            for rule in RULES:
                for M in M_GRID:
                    for p in P_GRID:
                        r = best_over(sk, sup, float(M), float(p), REF_MASS, rule, key="full")
                        rt = best_over(sk, sup, float(M), float(p), REF_MASS, rule, key="tech")
                        ro = best_over(sk, sup, float(M), float(p), REF_MASS, rule, key="ops")
                        r.update(ops_status_bestbeta=ro["ops_status"], ops_binding_bestbeta=ro["ops_binding"])
                        r.update(tech_status_bestbeta=rt["tech_status"], tech_binding_bestbeta=rt["tech_binding"],
                                 tech_score_bestbeta=rt["tech_score"], beta_tech=rt["beta"])
                        rows.append(r)
        print(f"  grid {sk} done", flush=True)
    write_csv("regime_grid.csv", rows)
    return rows


def stage_pods():
    rows = []
    masses = [2000, 3500, 5000, 7000, 10000, 13000, 17000, 20000, 25000, 30000]
    ps = np.geomspace(0.3, 1000, 15)
    for sk in SETTINGS:
        for sup in SUPPORTS:
            for rule in RULES:
                for m in masses:
                    for M in M_GRID[::2]:
                        best = None
                        for p in ps:
                            r = best_over(sk, sup, float(M), float(p), m, rule, key="ops")
                            cat = {"pass": 2, "marginal": 1, "fail": 0}[r["ops_status"]]
                            rank = (cat, min(r["BC_low"], 1.0) if cat else 0.0, -r["ops_nfail"], r["ops_score"])
                            if best is None or rank > best[0]:
                                best = (rank, r)
                        rows.append(best[1])
        print(f"  pods {sk} done", flush=True)
    write_csv("regime_pods.csv", rows)


def stage_summary(rows):
    out = []
    for sk in SETTINGS:
        for sup in SUPPORTS:
            for rule in RULES:
                sub = [r for r in rows if r["setting"] == sk and r["support"] == sup and r["rule"] == rule]
                for lab, key in (("technical", "tech_status_bestbeta"), ("ops", "ops_status_bestbeta"),
                                 ("full", "full_status")):
                    for st in ("pass", "marginal"):
                        ok = [r for r in sub if r[key] == "pass" or (st == "marginal" and r[key] == "marginal")]
                        if ok:
                            Ms = [r["Mach"] for r in ok]
                            ps = [r["p_Pa"] for r in ok]
                            out.append(dict(setting=sk, support=sup, rule=rule, layer=lab, level=st, n_points=len(ok),
                                            M_min=min(Ms), M_max=max(Ms), p_min=min(ps), p_max=max(ps),
                                            BC_base_max=max(r["BC_base"] for r in ok),
                                            BC_low_max=max(r["BC_low"] for r in ok),
                                            E_min=min(r["E_Wh_paxkm"] for r in ok)))
                        else:
                            from collections import Counter
                            bk = {"technical": "tech_binding_bestbeta", "ops": "ops_binding_bestbeta"}.get(lab, "full_binding")
                            bind = Counter(r[bk] for r in sub).most_common(3)
                            out.append(dict(setting=sk, support=sup, rule=rule, layer=lab, level=st, n_points=0,
                                            binding_most_common=";".join(f"{b}:{n}" for b, n in bind)))
    write_csv("regime_summary.csv", out)


def mach_limit(set_key, disc, p=10.0, m_pod=REF_MASS, rule="moving block 10 s", level="pass"):
    """Highest Mach on M_GRID at which discipline `disc` still passes (or is at least marginal)."""
    ok = [M for M in M_GRID if (r := best_over(set_key, "maglev", float(M), p, m_pod, rule, key="tech"))
          and (r[f"s_{disc}"] == "pass" or (level == "marginal" and r[f"s_{disc}"] == "marginal"))]
    return max(ok) if ok else float("nan")


def stage_levers():
    """What opens the regime: Mach limits per binding discipline, base vs lever (maglev, 10 Pa, 50 seats)."""
    rows = []
    levers = [
        ("base", {}, {}),
        ("active suspension 0.5 Hz", {"F_S": 0.5}, {}),
        ("survey/realignment 5 mm per km", {}, {"slope": 5e-6}),
        ("both of the above", {"F_S": 0.5}, {"slope": 5e-6}),
        ("straighter right-of-way R 500 km (elevated)", {}, {"R_avail": 500e3}),
        ("floating: 10x stiffer tethers (2.6 m)", {}, {"A_float": 2.59}),
        ("floating: active hold to 0.3 m", {}, {"A_float": 0.3}),
    ]
    for name, glob, sett in levers:
        saved_g = {k: globals()[k] for k in glob}
        globals().update(glob)
        for sk in SETTINGS:
            keys = {k: v for k, v in sett.items() if k in SETTINGS[sk]}
            if (sett and not keys) or ("(elevated)" in name and sk != "elevated"):
                continue
            saved_s = {k: SETTINGS[sk][k] for k in keys}
            SETTINGS[sk].update(keys)
            row = dict(lever=name, setting=sk)
            for disc in ("alignment", "route", "safety"):
                row[f"Mmax_{disc}_pass"] = mach_limit(sk, disc)
                row[f"Mmax_{disc}_marginal"] = mach_limit(sk, disc, level="marginal")
            rows.append(row)
            SETTINGS[sk].update(saved_s)
        globals().update(saved_g)
    # pod mass for safety: heavier / smaller-frontal pods
    for m in (5000, 10000, 17000, 30000):
        for sk in ("elevated", "buried", "floating"):
            rows.append(dict(lever=f"pod {m/1e3:g} t (safety only)", setting=sk,
                             Mmax_safety_pass=mach_limit(sk, "safety", m_pod=m),
                             Mmax_safety_marginal=mach_limit(sk, "safety", m_pod=m, level="marginal")))
    write_csv("regime_levers.csv", rows)


if __name__ == "__main__":
    if "--levers" in sys.argv:
        stage_levers()
        sys.exit()
    stage_glide_check()
    rows = stage_grid()
    stage_summary(rows)
    if "--no-pods" not in sys.argv:
        stage_pods()
    stage_levers()
