"""Skyloop physics re-review: screening numbers for the buoyant near-vacuum tube.

Order-of-magnitude screening, reusing BoomShield methods where they apply
(US-1976 atmosphere via boomsim, H2 lift density, UHMWPE tether model, assumed
mid-latitude wind profiles, vacuum-structure closure ratios, choked-leak model).

Run:  python3 skyloop-review/physics/skyloop_physics.py
Out:  CSVs in this folder + physics.png (figure made by make_figure.py)

All assumptions are marked (A) in comments; they are engineering estimates, not data.
"""
from __future__ import annotations

import csv
import math
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", ".."))
sys.path.insert(0, os.path.join(ROOT, "boom-shield", "sim"))
from boomsim.atmosphere import std_atm, R_STAR, G0  # noqa: E402  (read-only import)

# ------------------------------------------------------------------ constants
R_AIR, GAM, K_B = 287.05, 1.4, 1.380649e-23
M_AIR, M_H2 = 28.9644e-3, 2.016e-3
C_REF = float(std_atm(20000.0).c)          # 295.1 m/s: "Mach" = V / c at 20 km ambient
MACHS = [2.0, 4.0, 6.0]
T_GAS = 250.0                              # (A) tube gas temperature, K
C_GAS = math.sqrt(GAM * R_AIR * T_GAS)     # 317 m/s

# pod (docs/07: loaded 1500-3500 kg, 2-4 pax)
M_POD = 3000.0       # kg (A, central)
PAX = 3.0            # passengers per pod (A, central)
A_POD = 1.0          # m^2 frontal (A): ~1.13 m diameter needle
L_POD = 15.0         # m (A)
S_POD = 50.0         # m^2 wetted (A)
CD_POD = 0.25        # (A) supersonic slender body, wave + laminar friction, unchoked

# CFRP / materials (A, handbook-class)
E_QI, RHO_CF, NU = 70e9, 1600.0, 0.3       # quasi-isotropic CFRP laminate
RHO_CORE = 50.0                            # honeycomb / foam core
SIG_CF_ALLOW = 300e6                       # face compression allowable incl. SF
SF, KNOCK = 1.5, 0.56                      # safety factor; SP-8007-class knock-down
TF_MIN, H_MIN = 0.4e-3, 0.010              # minimum face / core gauge
LINER = 0.15                               # kg/m^2 metal-foil vacuum liner (A)
JOINTS = 1.20                              # frames, joints, flanges (A)

# guideway, stator, utilities [kg/m] (A): low / central / high
GUIDEWAY = (30.0, 60.0, 120.0)             # EDS conductor plates + brackets + local stiffening
STATOR = (50.0, 150.0, 300.0)              # long-stator LSM in launch/catch zones only
UTIL = 15.0                                # pumps, cables, sensors, comms

# envelopes and tethers (BoomShield methods)
FILM = 0.10           # kg/m^2 envelope laminate incl. seams (A; BoomShield 0.03-0.2)
CD_ENV = 0.6          # crosswind cylinder, supercritical Re (A)
RESERVE = 1.15        # lift reserve (BoomShield structure study)
V_OVER_H = 1.5        # tether top angle (BoomShield)
SIGMA_F, RHO_F, E_F = 1.2e9, 975.0, 116e9   # UHMWPE working stress, density, modulus
TETHER_SPACING = 3000.0                     # m (BoomShield); docs/08 says 30 km

Z_W = np.array([0, 2, 5, 8, 10.5, 12, 14, 16, 18, 20, 22, 25], float)
U_WINTER = np.array([8, 15, 25, 40, 50, 48, 35, 25, 20, 18, 20, 25], float)
U_DESIGN = np.array([20, 30, 50, 75, 90, 85, 65, 50, 42, 40, 42, 50], float)


def atm(z_km):
    return std_atm(np.asarray(z_km, float) * 1e3)


def U(profile, z_km):
    return float(np.interp(z_km, Z_W, profile))


def h2_lift(z_km):
    a = atm(z_km)
    return float(a.rho - a.p * M_H2 / (R_STAR * a.T))


def write_csv(name, rows):
    path = os.path.join(HERE, name)
    with open(path, "w", newline="") as f:
        keys = []
        for r in rows:
            keys += [k for k in r if k not in keys]
        w = csv.DictWriter(f, fieldnames=keys)
        w.writeheader()
        for r in rows:
            w.writerow({k: (f"{v:.4g}" if isinstance(v, float) else v) for k, v in r.items()})
    print(f"wrote {name} ({len(rows)} rows)")


# ================================================================== 1. tube structure
def wall_sandwich(dp, R):
    """Lightest CFRP sandwich wall for a long cylinder under external pressure dp.
    p_cr = 3 D_b / R^3, D_b = E tf h^2 / (2 (1-nu^2)); require p_cr >= dp*SF/KNOCK and
    face hoop compression dp*R/(2 tf) <= allowable. Returns kg/m^2, tf, h."""
    pd = dp * SF / KNOCK
    best = None
    for tf in np.geomspace(TF_MIN, 20e-3, 400):
        if dp * SF * R / (2 * tf) > SIG_CF_ALLOW:
            continue
        h = max(H_MIN, math.sqrt(2 * (1 - NU**2) * pd * R**3 / (3 * E_QI * tf)))
        m = 2 * tf * RHO_CF + h * RHO_CORE
        if best is None or m < best[0]:
            best = (m, tf, h)
    return best


def wall_monocoque(dp, R, E=E_QI, rho=RHO_CF):
    pd = dp * SF / KNOCK
    t = R * (4 * (1 - NU**2) * pd / E) ** (1 / 3)
    return t * rho, t


def tether_for(F_side, z_km, spacing=TETHER_SPACING, voh=None):
    """BoomShield tether model (altitude_trade.tether), side load F_side [N/m]."""
    V_OVER_H = voh or globals()["V_OVER_H"]
    z = np.linspace(0, z_km, 400)
    dz = (z[1] - z[0]) * 1e3
    a = atm(z)
    qz = 0.5 * a.rho * np.interp(z, Z_W, U_DESIGN) ** 2
    L = z_km * 1e3 / math.cos(math.atan(1 / V_OVER_H))
    d = 0.08
    for _ in range(80):
        D = float(np.sum(qz * 1.2 * d) * dz)
        H = F_side * spacing + 0.5 * D
        A = math.hypot(H, V_OVER_H * H) / SIGMA_F
        m_t = RHO_F * A * L
        V = V_OVER_H * H + 0.5 * m_t * G0
        T = math.hypot(H, V)
        A = T / SIGMA_F
        d_new = math.sqrt(4 * A / (math.pi * 0.65))
        if abs(d_new - d) < 1e-7:
            break
        d = d_new
    m_t = RHO_F * A * L
    theta = math.atan(1 / V_OVER_H)                       # from vertical
    k_ax = E_F * A / L
    k_h = k_ax * math.sin(theta) ** 2 + T / L * math.cos(theta) ** 2
    k_v = k_ax * math.cos(theta) ** 2
    holddown = (V + 0.5 * m_t * G0) / G0 / spacing      # kg of lift per m of line
    return dict(T=T, H=H, V=V, A=A, m_t=m_t, L=L, k_h=k_h, k_v=k_v, holddown=holddown, d=d)


def envelope(payload, z_km, film=FILM, cd=CD_ENV, voh=None, tethered=True):
    """H2 cylinder radius R_e lifting payload [kg/m] + film + tether hold-down."""
    lam = h2_lift(z_km)
    Ud = U(U_DESIGN, z_km)
    rho = float(atm(z_km).rho)
    R = 10.0
    for _ in range(200):
        F = 0.5 * rho * Ud**2 * cd * 2 * R
        t = tether_for(F, z_km, voh=voh)
        need = RESERVE * (payload + 2 * math.pi * R * film) + (t["holddown"] if tethered else 0.0)
        R_new = math.sqrt(need / (lam * math.pi))
        if abs(R_new - R) < 1e-4:
            break
        R = 0.5 * R + 0.5 * R_new
    return dict(R=R, F_design=F, lift=lam * math.pi * R**2, film=2 * math.pi * R * film,
                holddown=t["holddown"], tether=t, h2_kg_m=math.pi * R**2 * float(atm(z_km).p) * M_H2 / (R_STAR * float(atm(z_km).T)))


def stator_fraction(M, route_km=1000.0, a_g=0.3):
    L_acc = (M * C_REF) ** 2 / (2 * a_g * G0) / 1e3
    return min(1.0, 2 * L_acc / route_km), L_acc


def section1():
    rows = []
    for z in [0.0, 10.0, 15.0, 20.0, 25.0]:
        dp = float(atm(z).p)                     # near-vacuum inside
        for D in [2.0, 3.0, 4.0, 6.0]:
            R = D / 2
            m2, tf, h = wall_sandwich(dp, R)
            wall = (m2 + LINER) * math.pi * D * JOINTS
            mono, tmono = wall_monocoque(dp, R)
            mono_m = mono * math.pi * D * JOINTS
            fr, _ = stator_fraction(4.0)
            for case, gi in [("low", 0), ("central", 1), ("high", 2)]:
                pay = wall + GUIDEWAY[gi] + fr * STATOR[gi] + UTIL + 3.0   # +3 kg/m pods
                vac_lift = float(atm(z).rho) * math.pi * R**2
                row = dict(alt_km=z, D_m=D, case=case, dp_Pa=dp, face_mm=tf * 1e3, core_mm=h * 1e3,
                           wall_kg_m=wall, monocoque_CFRP_kg_m=mono_m, monocoque_t_mm=tmono * 1e3,
                           guideway_kg_m=GUIDEWAY[gi], stator_avg_kg_m=fr * STATOR[gi],
                           util_pods_kg_m=UTIL + 3.0, payload_kg_m=pay,
                           tube_own_vacuum_lift_kg_m=vac_lift, self_buoyancy_ratio=vac_lift / pay)
                if z >= 10.0:
                    env = envelope(pay, z)
                    row.update(H2_lift_kg_m3=h2_lift(z), env_R_m=env["R"], env_film_kg_m=env["film"],
                               holddown_kg_m=env["holddown"], total_lift_kg_m=env["lift"],
                               side_load_design_kN_m=env["F_design"] / 1e3,
                               tether_T_MN=env["tether"]["T"] / 1e6, tether_t=env["tether"]["m_t"] / 1e3,
                               tether_k_h_N_m=env["tether"]["k_h"], H2_t_per_km=env["h2_kg_m"])
                rows.append(row)
    write_csv("tube_mass_budget.csv", rows)

    # docs/07 invballoon check (vacuum lobes, R_env 20-50 m, 150 kg/m core)
    inv = []
    for Re in [20, 30, 40, 50]:
        gross = float(atm(20).rho) * math.pi * Re**2
        doc = {20: "sinks", 30: "marginal", 40: "+233", 50: "+449"}[Re]
        for ratio, lab in [(0.64, "ideal buckling-free floor"), (1.23, "BoomShield closure, best"),
                           (1.29, "BoomShield closure, D 250")]:
            inv.append(dict(R_env_m=Re, docs_claim_kg_m=doc, vacuum_gross_lift_kg_m=gross,
                            structure_ratio=ratio, basis=lab, structure_kg_m=ratio * gross,
                            net_after_150kg_core=gross * (1 - ratio) - 150.0))
        # the doc's own implied accounting: membrane only, rings 1 t / 50 m
        inv.append(dict(R_env_m=Re, docs_claim_kg_m=doc, vacuum_gross_lift_kg_m=gross,
                        structure_ratio=float("nan"), basis="docs implied: film 0.25 kg/m2 + 1 t ring/50 m",
                        structure_kg_m=2 * math.pi * Re * 0.25 + 20.0,
                        net_after_150kg_core=gross - 2 * math.pi * Re * 0.25 - 20 - 150))
    write_csv("invballoon_check.csv", inv)

    # envelope sensitivity at 20 km, D 3 m central payload
    pay = [r["payload_kg_m"] for r in rows if r["alt_km"] == 20.0 and r["D_m"] == 3.0 and r["case"] == "central"][0]
    sens = []
    for lab, kw in [("baseline: C_D 0.6, V/H 1.5, film 0.10", {}),
                    ("faired/lenticular envelope C_D 0.3", dict(cd=0.3)),
                    ("flatter tethers V/H 0.75", dict(voh=0.75)),
                    ("C_D 0.3 and V/H 0.75", dict(cd=0.3, voh=0.75)),
                    ("light film 0.05 kg/m2", dict(film=0.05)),
                    ("untethered (lift for payload + film only)", dict(tethered=False))]:
        e = envelope(pay, 20.0, **kw)
        sens.append(dict(case=lab, payload_kg_m=pay, env_R_m=e["R"], lift_kg_m=e["lift"], film_kg_m=e["film"],
                         holddown_kg_m=e["holddown"] if kw.get("tethered", True) else 0.0,
                         side_load_design_kN_m=e["F_design"] / 1e3, tether_T_MN=e["tether"]["T"] / 1e6,
                         tether_t=e["tether"]["m_t"] / 1e3, H2_t_per_km=e["h2_kg_m"],
                         lift_over_payload=e["lift"] / pay))
    write_csv("envelope_sensitivity.csv", sens)

    # docs/07 core-tube claim: 150 kg/m, 3 m bore, 10 Hz on 50 m spans
    m, R, L = 150.0, 1.5, 50.0
    t_eq = m / JOINTS / (math.pi * 3.0 * RHO_CF)
    EI = E_QI * math.pi * R**3 * t_eq
    f1 = math.pi / (2 * L**2) * math.sqrt(EI / m)
    EI_need = (10.0 * 2 * L**2 / math.pi) ** 2 * m
    print(f"docs/07 core 150 kg/m: t_eq {t_eq*1e3:.1f} mm, EI {EI:.3g}, f1(50 m, pinned) {f1:.2f} Hz; "
          f"10 Hz needs EI {EI_need:.3g} -> {EI_need/EI:.1f}x")
    return dict(f1=f1, EI=EI, EI_ratio=EI_need / EI)


# ================================================================== 2. aero / vacuum
def area_ratio_isentropic(M):
    g = GAM
    return (1 / M) * ((2 / (g + 1)) * (1 + (g - 1) / 2 * M**2)) ** ((g + 1) / (2 * (g - 1)))   # A/A*


def kantrowitz_beta(M):
    """Max blockage beta = A_pod/A_tube for a self-starting (unchoked) supersonic pod:
    normal shock at the annulus entry, then isentropic to sonic in the gap."""
    g = GAM
    p0r = (((g + 1) * M**2) / ((g - 1) * M**2 + 2)) ** (g / (g - 1)) * \
          ((g + 1) / (2 * g * M**2 - (g - 1))) ** (1 / (g - 1))
    # mass conservation with total-pressure loss: A_t/A_1 = (A*/A)(M1) / p0r
    At_over_A1 = (1 / area_ratio_isentropic(M)) / p0r
    return 1 - min(1.0, At_over_A1), 1 - 1 / area_ratio_isentropic(M)


def mfp(p, T=T_GAS, d=3.7e-10):
    return K_B * T / (math.sqrt(2) * math.pi * d**2 * p)


def section2():
    kr = []
    for M in [1.0, 1.2, 1.5, 2, 2.5, 3, 4, 5, 6, 8]:
        bk, bi = kantrowitz_beta(M)
        kr.append(dict(M_gas=M, V_m_s=M * C_GAS, beta_max_kantrowitz=bk, beta_max_isentropic=bi))
    write_csv("kantrowitz.csv", kr)

    rows = []
    for M in MACHS:
        V = M * C_REF
        for p in [0.1, 0.3, 1, 3, 10, 30, 100, 300, 1000]:
            rho = p / (R_AIR * T_GAS)
            q = 0.5 * rho * V**2
            for D in [2.0, 3.0, 4.0]:
                A_t = math.pi * D**2 / 4
                beta = A_POD / A_t
                bk, _ = kantrowitz_beta(V / C_GAS)
                choked = beta > bk
                D_un = q * A_POD * CD_POD
                D_ch_lo = 1.8 * q * A_POD                     # standing normal shock (C_p0 ~1.8)
                D_ch_hi = 0.5 * (GAM + 1) * rho * V**2 * A_t  # piston on full bore
                drag = D_ch_lo if choked else D_un
                Wh = drag * 1000 / PAX / 3600
                mu = 1.458e-6 * T_GAS**1.5 / (T_GAS + 110.4)
                Re = rho * V * L_POD / mu
                # heating
                q_fm = 0.5 * rho * V**3                       # free-molecular upper bound (alpha=1)
                q_sg = 1.7415e-4 * math.sqrt(rho / 0.05) * V**3   # Sutton-Graves, 5 cm nose
                q_nose = min(q_fm, q_sg) if mfp(p) / 0.05 < 0.1 else min(q_fm, max(q_sg, 0.5 * q_fm))
                T_nose = (q_nose / (0.85 * 5.670e-8) + T_GAS**4) ** 0.25
                rows.append(dict(Mach=M, V_m_s=V, p_Pa=p, D_tube_m=D, beta=beta, beta_kantrowitz=bk,
                                 choked=int(choked), q_Pa=q, mfp_mm=mfp(p) * 1e3, Re_pod=Re,
                                 drag_unchoked_N=D_un, drag_choked_lo_N=D_ch_lo, drag_choked_hi_N=D_ch_hi,
                                 drag_N=drag, aero_Wh_per_paxkm=Wh, drag_power_kW=drag * V / 1e3,
                                 pod_heat_kW_30pct=0.3 * drag * V / 1e3, nose_q_kW_m2=q_nose / 1e3,
                                 nose_T_eq_K=T_nose))
    write_csv("drag_vs_pressure.csv", rows)

    # pumping and leaks per km, D = 3 m tube, 20 km vs ground
    pump = []
    D = 3.0
    area, vol = math.pi * D * 1000, math.pi * D**2 / 4 * 1000
    for z in [0.0, 20.0]:
        a = atm(z)
        p0, T0 = float(a.p), float(a.T)
        mdot_cm2 = 0.6 * 1e-4 * p0 * math.sqrt(GAM / (R_AIR * T0)) * (2 / (GAM + 1)) ** ((GAM + 1) / (2 * (GAM - 1)))
        Q_hole = mdot_cm2 * R_AIR * T_GAS     # Pa m3/s per cm2 at tube temperature
        sources = {
            "outgassing CFRP/epoxy, unlined (1e-5 Pa m3/s/m2)": 1e-5 * area,
            "outgassing metal-lined, unbaked (1e-7)": 1e-7 * area,
            "seal permeation, 20 joints/km x 9.4 m (1e-7 Pa m3/s per m)": 20 * math.pi * D * 1e-7,
            "holes 0.1 cm2/km": 0.1 * Q_hole,
            "holes 1 cm2/km": 1.0 * Q_hole,
            "holes 10 cm2/km": 10.0 * Q_hole,
        }
        for name, Q in sources.items():
            for p in [1.0, 10.0]:
                S = Q / p
                pump.append(dict(alt_km=z, source=name, p_tube_Pa=p, Q_Pa_m3_s_per_km=Q, S_m3_s_per_km=S,
                                 P_pump_kW_per_km=10.0 * S, P_ideal_W_per_km=Q * math.log(p0 / p),
                                 pump_mass_kg_per_km=250 * S))
        # pump-down / re-evacuation of a 10 km section with S = 10 m3/s/km installed
        for S_km in [10.0, 50.0]:
            t = vol / S_km * math.log(p0 / 10.0)
            pump.append(dict(alt_km=z, source=f"re-evacuate to 10 Pa, S={S_km:g} m3/s/km", p_tube_Pa=10.0,
                             Q_Pa_m3_s_per_km=float("nan"), S_m3_s_per_km=S_km,
                             P_pump_kW_per_km=10 * S_km, P_ideal_W_per_km=float("nan"),
                             pump_mass_kg_per_km=250 * S_km, hours=t / 3600,
                             min_energy_kWh_per_km=p0 * vol / 3.6e6))
    write_csv("pumping_per_km.csv", pump)


# ================================================================== 3. propulsion and energy
def section3():
    rows = []
    for M in MACHS:
        V = M * C_REF
        KE = 0.5 * M_POD * V**2
        for a_g in [0.1, 0.2, 0.3, 0.5, 1.0]:
            a = a_g * G0
            L = V**2 / (2 * a)
            t = V / a
            P_peak = M_POD * a * V
            rows.append(dict(Mach=M, V_m_s=V, a_g=a_g, L_acc_km=L / 1e3, t_acc_s=t, P_peak_MW=P_peak / 1e6,
                             KE_pod_GJ=KE / 1e9, KE_pod_kWh=KE / 3.6e6, KE_pax_kWh=KE / PAX / 3.6e6,
                             KE_pod_tTNT=KE / 4.184e9, pods_in_launch_at_4000ph=4000 / 3600 * t,
                             station_avg_GW_4000ph=4000 / 3600 * KE / 1e9,
                             LSM_Hz_tau05=V / (2 * 0.5), LSM_Hz_tau1=V / (2 * 1.0)))
    write_csv("propulsion.csv", rows)

    en = []
    for M in MACHS:
        V = M * C_REF
        KE_pax_Wh = 0.5 * M_POD * V**2 / PAX / 3600
        for trip in [300, 1000, 3000, 9000]:
            for eta in [0.70, 0.80, 0.90]:      # round-trip regen efficiency (A)
                for LD in [300, 1000, 3000]:    # EDS lift/drag (Inductrack ~300 at 500 km/h, rises with v)
                    lev = M_POD * G0 / LD * 1000 / PAX / 3600
                    aero1 = 0.5 * (1.0 / (R_AIR * T_GAS)) * V**2 * A_POD * CD_POD * 1000 / PAX / 3600
                    kin = KE_pax_Wh * (1 - eta) / trip
                    fits = 2 * V**2 / (2 * 0.3 * G0) / 1e3 <= trip
                    en.append(dict(Mach=M, trip_km=trip, regen_eta=eta, EDS_LD=LD,
                                   kinetic_loss_Wh_paxkm=kin, levitation_Wh_paxkm=lev,
                                   aero_at_1Pa_Wh_paxkm=aero1, total_Wh_paxkm=kin + lev + aero1,
                                   no_regen_Wh_paxkm=KE_pax_Wh / trip + lev + aero1,
                                   fits_at_0p3g=int(fits)))
    write_csv("energy_per_paxkm.csv", en)


# ================================================================== 4. ride and alignment
def section4():
    rows = []
    for M in MACHS:
        V = M * C_REF
        r = dict(Mach=M, V_m_s=V)
        for a_g in [0.1, 0.2, 0.5]:
            r[f"R_horiz_km_{a_g}g"] = V**2 / (a_g * G0) / 1e3
        for a_g in [0.05, 0.1]:
            r[f"R_vert_km_{a_g}g"] = V**2 / (a_g * G0) / 1e3
        r["earth_curv_unload_g"] = V**2 / 6.371e6 / G0
        r["aero_lift_wing_m2_at_10Pa_CL0.3"] = M_POD * G0 / (0.5 * 10 / (R_AIR * T_GAS) * V**2 * 0.3)
        # brick-wall headway (safety)
        for b in [0.5, 1.0]:
            s = V**2 / (2 * b * G0)
            r[f"stop_km_{b}g"] = s / 1e3
            r[f"stop_s_{b}g"] = V / (b * G0)
            r[f"brickwall_headway_s_{b}g"] = s / V                      # distance headway = stop distance
            r[f"brickwall_pods_per_h_{b}g"] = 3600 / (s / V)
            r[f"brickwall_pax_per_h_{b}g"] = 3600 / (s / V) * PAX
        rows.append(r)
    write_csv("ride_geometry.csv", rows)

    # allowable tube deviation vs wavelength, and the floating-tube estimate
    al = []
    lams = np.geomspace(50, 50000, 61)
    for M in MACHS:
        V = M * C_REF
        for lam in lams:
            w = 2 * math.pi * V / lam
            al.append(dict(Mach=M, wavelength_m=lam, freq_Hz=V / lam,
                           A_allow_0p05g_m=0.05 * G0 / w**2, A_allow_0p15g_m=0.15 * G0 / w**2))
    write_csv("alignment_allowable.csv", al)

    # floating-tube motion: 20 km, central D=3 m case envelope
    env_rows = []
    pay = None
    for r in csv.DictReader(open(os.path.join(HERE, "tube_mass_budget.csv"))):
        if r["alt_km"] == "20" and r["D_m"] == "3" and r["case"] == "central":
            pay = float(r["payload_kg_m"])
    env = envelope(pay, 20.0)
    Re, t = env["R"], env["tether"]
    rho = float(atm(20).rho)
    Um = U(U_WINTER, 20)
    for up in [1.0, 3.0, 5.0]:
        dF = rho * Um * up * CD_ENV * 2 * Re                    # N/m gust increment
        x_top = dF * TETHER_SPACING / t["k_h"]
        for T_ax in [2e6, 10e6]:
            sag = dF * TETHER_SPACING**2 / (8 * T_ax)
            env_rows.append(dict(u_gust_m_s=up, dF_N_m=dF, tether_k_h_N_m=t["k_h"], tether_top_motion_m=x_top,
                                 axial_tension_MN=T_ax / 1e6, midspan_extra_m=sag, total_m=x_top + sag,
                                 wavelength_m=2 * TETHER_SPACING,
                                 allow_M2_0p05g=0.05 * G0 / (2 * math.pi * 2 * C_REF / (2 * TETHER_SPACING))**2,
                                 allow_M6_0p05g=0.05 * G0 / (2 * math.pi * 6 * C_REF / (2 * TETHER_SPACING))**2))
    # diurnal lift swing (BoomShield 9.2 %/day at 20 km) -> vertical
    dFv = 0.092 * env["lift"] * G0
    env_rows.append(dict(u_gust_m_s=float("nan"), dF_N_m=dFv, tether_k_h_N_m=t["k_v"],
                         tether_top_motion_m=dFv * TETHER_SPACING / t["k_v"], axial_tension_MN=float("nan"),
                         midspan_extra_m=float("nan"), total_m=dFv * TETHER_SPACING / t["k_v"],
                         wavelength_m=float("nan"), allow_M2_0p05g=float("nan"), allow_M6_0p05g=float("nan")))
    write_csv("floating_tube_motion.csv", env_rows)
    # stiffness needed to hold 3.6 cm at M6
    x_allow = 0.05 * G0 / (2 * math.pi * 6 * C_REF / (2 * TETHER_SPACING))**2
    dF3 = rho * Um * 3.0 * CD_ENV * 2 * Re
    k_need = dF3 * TETHER_SPACING / x_allow
    print(f"envelope R {Re:.1f} m, tether k_h {t['k_h']:.3g} N/m, k needed {k_need:.3g} ({k_need/t['k_h']:.0f}x); "
          f"x_allow M6 {x_allow*100:.1f} cm")
    # active station-keeping power to cancel dF with rotors of disc area 2R per metre
    for up in [1.0, 3.0]:
        dF = rho * Um * up * CD_ENV * 2 * Re
        A = 2 * Re          # m2 of disc per m of line (A)
        vi = math.sqrt(dF / (2 * rho * A))
        print(f"active hold u'={up}: {dF:.0f} N/m, ideal rotor power {dF*vi/1e3:.2f} kW/m = {dF*vi:.3g} MW per 1000 km")
    return env


# ================================================================== 5. thermal
def section5():
    rows = []
    for M in MACHS:
        V = M * C_REF
        for dT in [0.1, 1.0, 5.0, 20.0]:
            for alpha, mat in [(2e-6, "CFRP quasi-iso"), (23e-6, "aluminium rail")]:
                kappa = alpha * dT / 3.0
                rows.append(dict(Mach=M, dT_across_D_K=dT, material=mat, kappa_per_m=kappa,
                                 R_bow_km=1 / kappa / 1e3, lat_acc_g_if_free=V**2 * kappa / G0))
    write_csv("thermal_bow.csv", rows)
    for alpha, mat in [(2e-6, "CFRP"), (23e-6, "Al")]:
        print(f"{mat}: 60 K swing -> {alpha*60*50*1e3:.1f} mm per 50 m bay, {alpha*60*1e6:.0f} m per 1000 km")
    # station braking heat if grid lost: energy per pod into eddy rail
    for M in MACHS:
        V = M * C_REF
        KE = 0.5 * M_POD * V**2
        L = V**2 / (2 * 0.5 * G0)
        print(f"M{M:g}: KE {KE/1e9:.2f} GJ, eddy-brake heat {KE/L/1e3:.1f} kJ/m per pod over {L/1e3:.0f} km at 0.5 g")


# ================================================================== 6. safety
def section6():
    a20, a0 = atm(20), atm(0)
    vol_km = math.pi * 9 / 4 * 1000
    print(f"stored pV per km (D 3 m): 20 km {float(a20.p)*vol_km/1e6:.0f} MJ ({float(a20.p)*vol_km/4.184e6:.0f} kg TNT); "
          f"ground {float(a0.p)*vol_km/1e6:.0f} MJ ({float(a0.p)*vol_km/4.184e6:.0f} kg TNT)")
    c0 = float(a20.c)
    cs = c0 * math.sqrt(2 / (GAM + 1))
    print(f"breach inflow: sonic {cs:.0f} m/s at break, rho* {0.634*float(a20.rho):.3f}; "
          f"head of expansion {cs + 2*cs/(GAM-1):.0f} m/s")
    for M in MACHS:
        V = M * C_REF
        rho = float(a20.rho)
        q = 0.5 * rho * V**2
        dec = q * A_POD * CD_POD / M_POD / G0
        dec_ch = 1.8 * q * A_POD / M_POD / G0
        sg = 1.7415e-4 * math.sqrt(rho / 0.05) * V**3
        sg50 = 1.7415e-4 * math.sqrt(rho / 0.5) * V**3
        print(f"M{M:g} into ambient 20 km air: q {q/1e3:.0f} kPa, decel {dec:.1f}-{dec_ch:.1f} g, "
              f"nose heating {sg/1e6:.2f} MW/m2 (5 cm) {sg50/1e6:.2f} (50 cm)")


if __name__ == "__main__":
    s1 = section1()
    section2()
    section3()
    section4()
    section5()
    section6()
