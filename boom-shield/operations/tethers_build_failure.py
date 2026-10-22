"""BoomShield operations screen: faired / buoyant tethers, build-out, H2 make-up,
tether and segment failure, anchor land use.

Screening arithmetic only (numpy + matplotlib + boomsim.atmosphere, and the
assumed wind profiles and tether model of ../altitude/altitude_trade.py, imported
read-only so the numbers stay comparable with altitude-trade.md section 6).

Run:  python3 boom-shield/operations/tethers_build_failure.py
Writes into boom-shield/operations/:
  tether_sizing.csv        bare / ribbon / hair / rigid-faired / jet-only-faired tethers,
                           30 m x 1 km and 54 m x 1.5 km layers, rim C_D 1.0 and 0.3
  buoyant_tether.csv       H2-filled fairing: chord for neutral buoyancy, drag and gust load vs height
  fall_footprint.csv       falling tether / membrane: terminal speed, fall time, drift, impact energy
  h2_makeup.csv            make-up need, PV share, water-capture energy, supply alternatives
  failure_landuse.csv      N-1 redistribution, uplift, compartment loss, anchor and land-use numbers
  tethers_build_failure.png  four-panel figure

All "ASSUMED" constants are engineering guesses, flagged in tethers-build-failure.md.
"""
from __future__ import annotations

import csv
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "altitude"))
sys.path.insert(0, os.path.join(HERE, "..", "sim"))
import altitude_trade as AT  # noqa: E402  (read-only reuse of winds, atmosphere, tether constants)
from boomsim.atmosphere import G0  # noqa: E402

atm, U = AT.atm, AT.U
U_DES, U_WIN = AT.U_DESIGN, AT.U_WINTER
SIG, RHO_F, FILL = AT.SIGMA_F, AT.RHO_F, AT.FILL      # 1.2 GPa, 975 kg/m3, 0.65
VH = AT.V_OVER_H                                       # 1.5
E_F = 116e9            # Pa, SK78-class fibre modulus (~120 N/tex x 0.97; supplier sheets)
H_S = 20.0             # km, shield top
SPACING = 3000.0       # m between tether stations (altitude-trade 6)
M_H2 = AT.M_H2
R_STAR = 8.314462618


def wcsv(name, rows):
    with open(os.path.join(HERE, name), "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        w.writeheader()
        for r in rows:
            w.writerow({k: (f"{v:.4g}" if isinstance(v, float) else v) for k, v in r.items()})


# =============================================================== A. tether sizing
def side_load(T_layer, W, U_s, cd_rim=1.0):
    """Shield side load per m of corridor [N/m], altitude-trade 3 form."""
    a = atm(H_S)
    q = 0.5 * a.rho * U_s ** 2
    cf = AT.skin_cf(a.rho, AT.mu_air(a.T), U_s, L=W)
    return q * (cd_rim * T_layer + 2 * W * cf)


# tether "types": drag coefficient on frontal thickness, thickness / fibre-rope diameter,
# extra linear mass as a function of rope diameter d [m], and the height band faired
def fair_rigid_mass(d):
    t = 1.2 * d                      # fairing thickness encloses the rope (ASSUMED)
    c = t / 0.22                     # t/c 0.22, typical cable-fairing section (ASSUMED)
    return 2.1 * c * 1.0 + 0.5       # 1.0 kg/m2 shell incl. ribs, +0.5 kg/m collars/bearings (ASSUMED)


TYPES = {
    "bare":            dict(cd=1.2,  tf=1.0, extra=lambda d: 0.0, band=None),
    "ribbon":          dict(cd=0.8,  tf=1.0, extra=lambda d: 0.4, band=None),
    "hair":            dict(cd=1.0,  tf=1.0, extra=lambda d: 0.3, band=None),
    "rigid_faired":    dict(cd=0.15, tf=1.2, extra=fair_rigid_mass, band=None),
    "rigid_faired_cd0.3": dict(cd=0.30, tf=1.2, extra=fair_rigid_mass, band=None),
    "rigid_jet_only":  dict(cd=0.15, tf=1.2, extra=fair_rigid_mass, band=(6.0, 16.0)),
}


def tether(F_side, kind="bare", prof=U_DES, spacing=SPACING, W=1000.0):
    """altitude_trade.tether() generalised to fairings (same assumptions:
    drag integrated as for a vertical tether, H_top = F*s + D/2, V/H = 1.5, 34 deg)."""
    ty = TYPES[kind]
    z = np.linspace(0.0, H_S, 400)
    dz = (z[1] - z[0]) * 1e3
    q = 0.5 * atm(z).rho * U(prof, z) ** 2
    infair = np.ones_like(z, bool) if ty["band"] is None else (z >= ty["band"][0]) & (z <= ty["band"][1])
    L = H_S * 1e3 / np.cos(np.arctan(1 / VH))
    frac_f = infair.mean()
    d = 0.15
    for _ in range(200):
        cd = np.where(infair, ty["cd"], TYPES["bare"]["cd"])
        tf = np.where(infair, ty["tf"], 1.0)
        D = np.sum(q * cd * tf * d) * dz
        H = F_side * spacing + 0.5 * D
        A_f = 1.8 * H / SIG
        m_fib = RHO_F * A_f * L
        m_ext = ty["extra"](d) * L * frac_f
        V = VH * H + 0.5 * (m_fib + m_ext) * G0
        T = np.hypot(H, V)
        A_f = T / SIG
        d_new = np.sqrt(4 * A_f / (np.pi * FILL))
        if abs(d_new - d) < 1e-7:
            break
        d = d_new
    m_fib = RHO_F * A_f * L
    m_ext = ty["extra"](d) * L * frac_f
    lift_cost = (V + 0.5 * (m_fib + m_ext) * G0) / G0 / (spacing * W)
    return dict(D_MN=D / 1e6, H_top_MN=H / 1e6, T_top_MN=T / 1e6, d_cm=100 * d,
                fibre_kg_per_m=RHO_F * A_f, fibre_t=m_fib / 1e3, fairing_t=m_ext / 1e3,
                total_t=(m_fib + m_ext) / 1e3, liftcost_kg_m2=lift_cost,
                drag_share_of_H=0.5 * D / H, L_km=L / 1e3)


def section_a():
    rows = []
    layers = [("30m_x_1km", 30.0, 1000.0), ("54m_x_1.5km", 54.0, 1500.0), ("77m_x_2km", 77.0, 2000.0)]
    for lname, T, W in layers:
        for cd_rim in (1.0, 0.3):
            F = side_load(T, W, U(U_DES, H_S), cd_rim)
            for kind in TYPES:
                r = tether(F, kind, W=W)
                rw = tether(side_load(T, W, U(U_WIN, H_S), cd_rim), kind, prof=U_WIN, W=W)
                rows.append(dict(layer=lname, rim_CD=cd_rim, tether=kind,
                                 Fside_design_kN_per_m=F / 1e3, **r,
                                 D_winter_MN=rw["D_MN"], T_top_winter_MN=rw["T_top_MN"]))
    wcsv("tether_sizing.csv", rows)
    return rows


# =============================================================== B. H2-filled buoyant fairing
SKIN_ENV = 0.25        # kg/m2, H2-tight aerostat laminate (ASSUMED, aerostat-class)
CDC = 0.012            # chord-based profile drag, aligned, Re ~1e7 (ASSUMED)


def h2_net(z):
    a = atm(z)
    rho_h2 = a.p * M_H2 / (R_STAR * a.T)
    return a.rho - rho_h2


def section_b(m_line):
    """Chord of a t/c section whose H2 volume floats m_line kg/m plus its own skin."""
    rows = []
    for tc in (0.12, 0.20):
        for z in (0.0, 2.0, 5.0, 8.0, 10.5, 14.0, 18.0):
            dl = h2_net(z)
            k = 0.685 * tc * dl                      # kg lift per m per c^2
            s = 2.05 * SKIN_ENV                      # skin kg/m per m of chord
            c = (s + np.sqrt(s * s + 4 * k * m_line)) / (2 * k)
            a = atm(z)
            q = 0.5 * a.rho * U(U_DES, z) ** 2
            gust = 0.10 if z < 12 else 0.05          # w'/U, ASSUMED (tropospheric / stratospheric)
            rows.append(dict(t_over_c=tc, z_km=z, net_lift_kg_m3=dl, chord_m=c, thickness_m=tc * c,
                             section_m2=0.685 * tc * c * c,
                             drag_aligned_N_per_m=q * CDC * c,
                             drag_bare_N_per_m=q * 1.2 * 0.20,
                             gust_side_N_per_m=q * 2 * np.pi * gust * c))
    # small fairing (t 0.24 m, c 1.1 m) filled with H2: what fraction of the rope weight it floats
    for z in (0.0, 5.0, 10.5, 18.0):
        A = 0.685 * 0.24 * 1.09 - np.pi / 4 * 0.20 ** 2
        rows.append(dict(t_over_c=0.22, z_km=z, net_lift_kg_m3=h2_net(z), chord_m=1.09, thickness_m=0.24,
                         section_m2=A, drag_aligned_N_per_m=np.nan, drag_bare_N_per_m=np.nan,
                         gust_side_N_per_m=np.nan, lift_kg_per_m=A * h2_net(z),
                         share_of_line_weight=A * h2_net(z) / m_line))
    keys = sorted({k for r in rows for k in r}, key=lambda k: list(rows[-1].keys()).index(k)
                  if k in rows[-1] else -1)
    for r in rows:
        for k in rows[-1]:
            r.setdefault(k, np.nan)
    wcsv("buoyant_tether.csv", rows)
    return rows


# =============================================================== C/D. falling objects
def fall(m_lin, width, cd, z0=H_S, prof=U_WIN, n=800):
    """Terminal-velocity descent of a line element (broadside) or membrane.
    Horizontal drift assumes the element moves with the wind (low inertia).
    m_lin [kg/m] and width [m] (or areal mass with width=1)."""
    z = np.linspace(0.0, z0, n)
    rho = atm(z).rho
    vt = np.sqrt(2 * m_lin * G0 / (rho * cd * width))
    t = np.trapezoid(1 / vt, z * 1e3)
    drift = np.trapezoid(U(prof, z) / vt, z * 1e3)
    return vt, t, drift


def section_d_fall(m_fib, d):
    tf = 1.2 * d
    mf = m_fib + fair_rigid_mass(d)
    cases = [
        ("bare rope, broadside", m_fib, d, 1.2),
        ("faired, fairing aligned with fall", mf, tf, 0.15),
        ("faired, flat/tumbling", mf, tf / 0.22, 1.0),
        ("partly buoyant (net 2 kg/m, bluff 1 m)", 2.0, 1.0, 1.0),
        ("deflated lid membrane 3.9 kg/m2", 3.9, 1.0, 1.2),
    ]
    rows = []
    for name, m, w, cd in cases:
        for pname, prof in (("winter_mean", U_WIN), ("design", U_DES)):
            vt, t, dr = fall(m, w, cd, prof=prof)
            rows.append(dict(case=name, wind=pname, m_kg_per_m=m, width_m=w, CD=cd,
                             v_ground_m_s=vt[0], v_10km_m_s=float(np.interp(10, np.linspace(0, H_S, len(vt)), vt)),
                             fall_min=t / 60, drift_km=dr / 1e3,
                             KE_kJ_per_m_ground=0.5 * m * vt[0] ** 2 / 1e3))
    wcsv("fall_footprint.csv", rows)
    return rows


# =============================================================== C. H2 make-up
def section_c():
    rows = []
    a = atm(H_S)
    ppmv = 5.0
    w_mass = ppmv * 1e-6 * 18.015 / 28.964           # kg water / kg air
    rho_w = w_mass * a.rho                           # kg water per m3 air
    air_m3_per_kgH2 = 8.94 / rho_w                    # 8.94 kg water per kg H2
    air_kg_per_kgH2 = air_m3_per_kgH2 * a.rho
    # cold trap: cool air 217 -> 178 K (ice p_sat ~0.2 of ambient pw), COP 0.5, recuperation 0 / 90 %
    cool = 1005 * 40 * air_kg_per_kgH2 / 3.6e6       # kWh_th per kg H2
    fan = air_m3_per_kgH2 * 30 / 0.6 / 3.6e6          # kWh per kg H2 at 30 Pa bed loss
    regen = (8.94 * 2.8e6 + 8.94 / 0.075 * 900 * 150) / 3.6e6
    PV_YIELD = 3000 * 0.12                            # kWh/m2/yr: ~3,000 kWh/m2/yr horizontal at 20 km, 12 % thin film (ASSUMED)
    for inv_Mt, W in ((2.0, 1500.0), (4.0, 2000.0)):
        for leak in (0.05, 0.10, 0.20):
            need_t_km_yr = inv_Mt * 1e9 / 4000 * leak / 1e3
            per_m2 = need_t_km_yr * 1e3 / (1000 * W)                       # kg/m2/yr
            kwh_m2 = per_m2 * 52.5
            flow_air = need_t_km_yr * 1e3 * air_m3_per_kgH2 / 3.156e7      # m3/s per km
            rows.append(dict(inventory_Mt=inv_Mt, width_m=W, leak_pct_yr=100 * leak,
                             makeup_t_per_km_yr=need_t_km_yr, makeup_Mt_yr=inv_Mt * leak,
                             makeup_kg_per_m2_yr=per_m2, electrolysis_kWh_per_m2_yr=kwh_m2,
                             PV_area_share_pct=100 * kwh_m2 / PV_YIELD,
                             air_m3_s_per_km=flow_air, contactor_m2_per_km_at_2ms=flow_air / 2.0,
                             corridor_kg_s=inv_Mt * leak * 1e9 / 3.156e7,
                             per_supply_tether_g_s_every10th=inv_Mt * leak * 1e12 / 3.156e7 / 133,
                             tanker_flights_per_day_6t=inv_Mt * leak * 1e6 / 6.2 / 365))
    consts = dict(water_kg_per_m3_air=rho_w, air_m3_per_kgH2=air_m3_per_kgH2,
                  coldtrap_kWh_e_per_kgH2_no_recup=cool / 0.5, coldtrap_kWh_e_per_kgH2_90pc_recup=0.1 * cool / 0.5,
                  sorbent_fan_kWh_per_kgH2=fan, sorbent_regen_kWh_per_kgH2=regen,
                  electrolysis_kWh_per_kgH2=52.5, water_hose_head_MPa=1000 * G0 * 24e3 / 1e6,
                  water_lift_kWh_per_kgH2=8.94 * G0 * 24e3 / 3.6e6)
    # aircraft exhaust water: Concorde-class ~20 t/h fuel at cruise (ASSUMED, background), EI 1.23
    fuel_per_km = 20000 / (590 * 3.6)                  # kg fuel per km
    for flights in (20, 50, 100):
        h2o = fuel_per_km * 1.23 * flights * 365 / 1e3  # t water per km per yr
        consts[f"exhaust_water_t_km_yr_{flights}fl"] = h2o
        consts[f"exhaust_H2_equiv_t_km_yr_{flights}fl"] = h2o / 8.94
    # H2 feed hose: 24 km, 20 bar inlet, isothermal 250 K, Darcy f 0.015 (ASSUMED)
    for ID, mdot, p0 in ((0.05, 0.024, 20e5), (0.06, 0.048, 20e5), (0.06, 0.095, 40e5), (0.08, 0.19, 40e5)):
        p, Lseg = p0, 24e3 / 400
        for _ in range(400):
            rho = p * M_H2 / (R_STAR * 250)
            v = mdot / (rho * np.pi * ID ** 2 / 4)
            p -= 0.015 * Lseg / ID * 0.5 * rho * v ** 2 + rho * G0 * Lseg * 0.83
            if p <= 0:
                break
        consts[f"hose_ID{int(ID*1000)}mm_{int(mdot*1000)}gs_{int(p0/1e5)}bar_in_p_top_bar"] = max(p, 0) / 1e5
    for k, v in consts.items():
        rows.append(dict(inventory_Mt=k, width_m=v))
    wcsv("h2_makeup.csv", rows)
    return rows, consts


# =============================================================== D/E. failure, land use
def section_e(sz):
    psf = [r for r in sz if r["layer"] == "54m_x_1.5km" and r["rim_CD"] == 1.0 and r["tether"] == "bare"][0]
    T = psf["T_top_MN"] * 1e6
    H = psf["H_top_MN"] * 1e6
    V = VH * H
    out = {}
    # elastic recoil
    eps = SIG / E_F
    out["strain_pct"] = 100 * eps
    out["release_wave_s"] = psf["L_km"] * 1e3 / np.sqrt(E_F / RHO_F)
    out["recoil_v_m_s"] = SIG / np.sqrt(E_F * RHO_F)
    out["stored_energy_GJ"] = 0.5 * SIG * eps * (psf["fibre_t"] * 1e3 / RHO_F) / 1e9
    out["stored_energy_t_TNT"] = out["stored_energy_GJ"] / 4.184
    # N-1: neighbours pick up half each; sudden release DAF 1.3-2.0 on the increment
    for daf in (1.0, 1.5, 2.0):
        out[f"neighbour_load_factor_DAF{daf}"] = 1 + 0.5 * daf
    out["two_adjacent_lost_neighbour_factor"] = 2.0
    # uplift on loss of hold-down: 0.55 kg/m2 net over a 3 km bay
    Wd = 1500.0
    up_pa = 0.55 * G0
    m_add = 4.5 + atm(H_S).rho * np.pi * Wd / 4
    out["uplift_Pa"] = up_pa
    out["bay_uplift_MN"] = up_pa * SPACING * Wd / 1e6
    out["rise_250m_s_inertial"] = np.sqrt(2 * 250 / (up_pa / m_add))
    out["rise_terminal_m_s"] = np.sqrt(2 * up_pa / (atm(H_S).rho * 1.2))
    # along-track rim tendon catenary: span 3 km sag 5 %, then 6 km span same tendon length
    w = side_load(54.0, Wd, U(U_DES, H_S))
    s0 = 0.05 * SPACING
    Lc = 2 * (SPACING + 8 * s0 ** 2 / (3 * SPACING))
    s1 = np.sqrt((Lc - 2 * SPACING) * 3 * 2 * SPACING / 8)
    out["tendon_T_3km_MN"] = w * SPACING ** 2 / (8 * s0) / 1e6
    out["tendon_T_6km_MN"] = w * (2 * SPACING) ** 2 / (8 * s1) / 1e6
    out["bow_6km_m"] = s1
    # compartment: 250 m x 1.5 km x 55 m
    a = atm(H_S)
    rho_h2 = a.p * M_H2 / (R_STAR * a.T)
    vol = 250 * Wd * 55
    out["compartment_H2_t"] = vol * rho_h2 / 1e3
    out["compartment_TNT_kt"] = vol * rho_h2 * 120e6 / 4.184e12
    out["compartment_structure_t"] = 3.96 * 250 * Wd / 1e3
    out["ground_gas_fraction_pct"] = 100 * (a.p / 101325) * (288.15 / a.T)
    out["ground_mattress_m"] = 3.96 / (1.225 - 101325 * M_H2 / (R_STAR * 288.15))
    # flat-mattress ascent: 10 % free lift on 3.96 kg/m2, broadside C_D 1.2
    zz = np.linspace(0, H_S, 400)
    w_up = np.sqrt(2 * 0.1 * 3.96 * G0 / (atm(zz).rho * 1.2))
    out["ascent_w_ground_m_s"] = w_up[0]
    out["ascent_h"] = np.trapezoid(1 / w_up, zz * 1e3) / 3600
    jet = (zz >= 8) & (zz <= 13)
    out["ascent_jet_min"] = np.trapezoid(1 / w_up[jet], zz[jet] * 1e3) / 60
    # side load on a 250 m x 1.5 km module at 10.5 km, mattress thickness 3.96/dlift
    T10 = 3.96 / h2_net(10.5)
    q10 = 0.5 * atm(10.5).rho * np.array([25.0, 35.0, 90.0]) ** 2
    out["module_T_at_10km_m"] = T10
    for qq, uu in zip(q10, (25, 35, 90)):
        out[f"module_force_narrow_end_MN_U{uu}"] = qq * (T10 * 250 + 2 * 0.002 * 250 * Wd) / 1e6
        out[f"module_force_broadside_MN_U{uu}"] = qq * (T10 * Wd + 2 * 0.002 * 250 * Wd) / 1e6
    # tug power to move a 250 m module at 5 m/s relative, 20 km
    q5 = 0.5 * atm(H_S).rho * 25
    out["tug_kW_5ms"] = q5 * (55 * 250 + 2 * 0.002 * 250 * Wd) * 5 / 0.6 / 1e3
    out["compartments_absorbed_per_bay"] = 0.55 / 4.5 * SPACING / 250
    # anchors: gravity block resisting V (SF 1.5) and H_bottom = H + D/2 by friction mu 0.5 (SF 1.5)
    Hb = H + 0.5 * psf["D_MN"] * 1e6
    m_block = max(1.5 * V, 1.5 * Hb / 0.5 + V) / G0
    out["anchor_V_MN"] = V / 1e6
    out["anchor_H_MN"] = Hb / 1e6
    out["gravity_block_t"] = m_block / 1e3
    out["gravity_block_m3"] = m_block / 2400
    for r_ex in (100.0, 150.0, 300.0):
        ha_site = 1.5 + np.pi * r_ex ** 2 / 1e4
        out[f"site_ha_r{int(r_ex)}"] = ha_site
        out[f"ha_per_km_one_row_r{int(r_ex)}"] = ha_site * 1000 / SPACING
        out[f"ha_per_km_two_rows_r{int(r_ex)}"] = 2 * ha_site * 1000 / SPACING
    out["tether_ground_offset_km"] = H_S * np.tan(np.arctan(1 / VH))
    out["build_film_km2"] = 4000 * Wd / 1e3 * 2.3
    out["build_structure_Mt"] = 4000e3 * Wd * 3.96 / 1e9
    for yrs, days in ((10, 100), (10, 200), (20, 100)):
        out[f"bays_per_window_day_{yrs}y_{days}d"] = (4000e3 / SPACING) / (yrs * days)
        out[f"km2_lid_per_window_day_{yrs}y_{days}d"] = 4000 * Wd / 1e3 / (yrs * days)
    wcsv("failure_landuse.csv", [dict(quantity=k, value=v) for k, v in out.items()])
    return out


# =============================================================== figure
C = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100"]
INK, INK2, GRID = "#0b0b0b", "#52514e", "#e4e3df"


def figure(sz, bt, fl, consts):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    plt.rcParams.update({"font.size": 8.5, "axes.edgecolor": INK2, "axes.labelcolor": INK,
                         "xtick.color": INK2, "ytick.color": INK2, "axes.spines.top": False,
                         "axes.spines.right": False})
    fig, ax = plt.subplots(2, 2, figsize=(11, 8.2))

    # (a) top tension, psf layer, by tether type and rim CD
    a = ax[0, 0]
    kinds = ["bare", "ribbon", "hair", "rigid_faired", "rigid_jet_only"]
    labels = ["bare", "ribbon", "hair", "rigid\nfaired", "rigid,\njet only"]
    x = np.arange(len(kinds))
    for j, (rim, col, lab) in enumerate(((1.0, C[0], "bluff rim, C_D 1.0"), (0.3, C[2], "rounded rim, C_D 0.3"))):
        vals = [[r for r in sz if r["layer"] == "54m_x_1.5km" and r["rim_CD"] == rim and r["tether"] == k][0]
                for k in kinds]
        a.bar(x + (j - 0.5) * 0.38, [v["T_top_MN"] for v in vals], 0.36, color=col, label=lab, zorder=3)
        for xi, v in zip(x, vals):
            a.text(xi + (j - 0.5) * 0.38, v["T_top_MN"] + 0.4, f"{v['total_t']:.0f} t", ha="center",
                   fontsize=7, color=INK2)
    a.set_xticks(x, labels)
    a.set_ylabel("tether top tension, design wind [MN]")
    a.set_title("(a) 54 m x 1.5 km layer, 3 km spacing: fairing trims ~10 %,\n     the shield rim matters more (labels: tether mass)",
                loc="left", fontsize=9)
    a.grid(axis="y", color=GRID, zorder=0)
    a.set_ylim(0, 35)
    a.legend(frameon=False, fontsize=7.5, loc="upper center", ncol=2)

    # (b) per-metre lateral load vs height: bare, rigid faired, buoyant fairing aligned and gust
    b = ax[0, 1]
    z = np.linspace(0, 20, 200)
    q = 0.5 * atm(z).rho * U(U_DES, z) ** 2
    psf = [r for r in sz if r["layer"] == "54m_x_1.5km" and r["rim_CD"] == 1.0 and r["tether"] == "bare"][0]
    d = psf["d_cm"] / 100
    b.plot(q * 1.2 * d, z, color=C[0], lw=2, label=f"bare rope, d {100*d:.0f} cm")
    b.plot(q * 0.15 * 1.2 * d, z, color=C[2], lw=2, label="rigid fairing, aligned")
    bb = [r for r in bt if r["t_over_c"] == 0.12]
    zb = np.array([r["z_km"] for r in bb])
    b.plot([r["drag_aligned_N_per_m"] for r in bb], zb, "o-", color=C[1], lw=2, ms=4,
           label="H2-buoyant fairing (t/c 0.12), aligned")
    b.plot([r["gust_side_N_per_m"] for r in bb], zb, "s--", color=C[3], lw=2, ms=4,
           label="same, gust side load (w'/U 5-10 %)")
    b.set_xscale("log")
    b.set_xlabel("lateral load per metre of tether, design wind [N/m]")
    b.set_ylabel("height [km]")
    b.set_title("(b) A fairing big enough to float the rope is a wing:\n     its gust load dwarfs the bare-rope drag",
                loc="left", fontsize=9)
    b.grid(color=GRID, which="both", lw=0.5)
    b.legend(frameon=False, fontsize=7, loc="upper right")

    # (c) fall: drift vs fall time
    c = ax[1, 0]
    names = ["bare rope, broadside", "faired, fairing aligned with fall", "faired, flat/tumbling",
             "partly buoyant (net 2 kg/m, bluff 1 m)", "deflated lid membrane 3.9 kg/m2"]
    short = ["bare rope", "faired, edge-on", "faired, flat", "partly buoyant", "lid membrane"]
    for i, (n, s) in enumerate(zip(names, short)):
        for wind, mk in (("winter_mean", "o"), ("design", "^")):
            r = [r for r in fl if r["case"] == n and r["wind"] == wind][0]
            col = (C + ["#6250d6"])[i if i < 4 else 4]
            c.scatter(r["fall_min"], r["drift_km"], s=46, marker=mk, color=col, edgecolor="white", lw=1.5, zorder=3)
            if wind == "design":
                off = [(-30, 8), (8, -4), (8, -4), (-150, 2), (-140, -3)][i]
                c.annotate(f"{s}: {r['v_ground_m_s']:.0f} m/s, {r['KE_kJ_per_m_ground']:.2g} kJ/m" if i < 4 else
                           f"{s}: {r['v_ground_m_s']:.0f} m/s, {r['KE_kJ_per_m_ground']*1e3:.0f} J/m2",
                           (r["fall_min"], r["drift_km"]), xytext=off, textcoords="offset points",
                           fontsize=7, color=INK2)
    c.set_xscale("log")
    c.set_yscale("log")
    from matplotlib.ticker import FuncFormatter, NullFormatter
    for axis in (c.xaxis, c.yaxis):
        axis.set_major_formatter(FuncFormatter(lambda v, _: f"{v:g}"))
        axis.set_minor_formatter(NullFormatter())
    c.set_xlim(1.2, 80)
    c.set_xticks([2, 5, 10, 20, 50])
    c.set_xlabel("fall time from 20 km [min]")
    c.set_ylabel("downwind drift of the top end [km]")
    c.set_title("(c) What falls, how fast, how far (o winter mean, ^ design wind)", loc="left", fontsize=9)
    c.grid(color=GRID, which="both", lw=0.5)

    # (d) energy per kg H2 by route
    dd = ax[1, 1]
    routes = [("electrolysis only (ground or lid)", 52.5),
              ("+ H2 hose compression, 40 bar", 52.5 + 1.5),
              ("+ water hose lift (24 km head)", 52.5 + consts["water_lift_kWh_per_kgH2"]),
              ("+ lid sorbent, wind-driven", 52.5 + consts["sorbent_regen_kWh_per_kgH2"]),
              ("+ lid sorbent, fan 30 Pa", 52.5 + consts["sorbent_regen_kWh_per_kgH2"] + consts["sorbent_fan_kWh_per_kgH2"]),
              ("+ lid cold trap, 90 % recuperated", 52.5 + consts["coldtrap_kWh_e_per_kgH2_90pc_recup"]),
              ("+ lid cold trap, no recuperation", 52.5 + consts["coldtrap_kWh_e_per_kgH2_no_recup"])]
    y = np.arange(len(routes))[::-1]
    cols = [C[0], C[0], C[3], C[2], C[2], C[1], C[1]]
    dd.barh(y, [r[1] for r in routes], color=cols, height=0.6, zorder=3)
    for yi, r in zip(y, routes):
        dd.text(r[1] * 1.08, yi, f"{r[1]:,.0f}", va="center", fontsize=7.5, color=INK2)
    dd.set_yticks(y, [r[0] for r in routes])
    dd.set_xscale("log")
    dd.set_xlim(30, 3e5)
    dd.set_xlabel("energy per kg of make-up H2 [kWh]")
    dd.set_title("(d) Make-up H2: lid electrolysis is cheap; getting 9 kg of water\n     from 5 ppmv air is not",
                 loc="left", fontsize=9)
    dd.grid(axis="x", color=GRID, which="both", lw=0.5, zorder=0)
    fig.tight_layout()
    fig.savefig(os.path.join(HERE, "tethers_build_failure.png"), dpi=150)


def main():
    sz = section_a()
    psf_bare = [r for r in sz if r["layer"] == "54m_x_1.5km" and r["rim_CD"] == 1.0 and r["tether"] == "bare"][0]
    bt = section_b(psf_bare["fibre_kg_per_m"])
    fl = section_d_fall(psf_bare["fibre_kg_per_m"], psf_bare["d_cm"] / 100)
    mk, consts = section_c()
    ev = section_e(sz)
    figure(sz, bt, fl, consts)
    for r in sz:
        print(f"{r['layer']:12s} rim{r['rim_CD']:.1f} {r['tether']:18s} F {r['Fside_design_kN_per_m']:.2f} kN/m "
              f"D {r['D_MN']:.2f} H {r['H_top_MN']:.2f} T {r['T_top_MN']:.1f} MN d {r['d_cm']:.1f} cm "
              f"fib {r['fibre_kg_per_m']:.1f} kg/m tot {r['total_t']:.0f} t (fair {r['fairing_t']:.0f}) "
              f"lift {r['liftcost_kg_m2']:.3f} dragshare {r['drag_share_of_H']:.2f} Dwin {r['D_winter_MN']:.2f}")
    print()
    for r in bt:
        print({k: (round(v, 3) if isinstance(v, float) else v) for k, v in r.items()})
    print()
    for r in fl:
        print(f"{r['case']:42s} {r['wind']:12s} v0 {r['v_ground_m_s']:.1f} v10 {r['v_10km_m_s']:.1f} "
              f"t {r['fall_min']:.1f} min drift {r['drift_km']:.1f} km KE {r['KE_kJ_per_m_ground']:.3g} kJ/m")
    print()
    for r in mk:
        print({k: (round(v, 4) if isinstance(v, float) else v) for k, v in r.items()})
    print()
    for k, v in ev.items():
        print(f"{k:40s} {v:.4g}")


if __name__ == "__main__":
    main()
