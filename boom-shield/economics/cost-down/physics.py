#!/usr/bin/env python3
"""BoomShield cost-down: physics helpers (lid mass re-optimisation, stratospheric PV, lid thermal).

Run:  python3 physics.py   -> lid_reopt.csv, min_structure.csv, pv_resource.csv, pv_thermal.csv

Read-only use of the trunk: constants follow structure-concept.md / psf_sweep.md / cost_model.py.
"""
from __future__ import annotations

import csv
import math
from pathlib import Path

HERE = Path(__file__).resolve().parent

# 20 km standard atmosphere and H2 at ambient temperature (psf_sweep / structure-concept)
T_AMB = 216.65
P_AMB = 5529.0
RHO_A = 0.0889
C_A = 295.1
RHO_G = P_AMB * 2.016e-3 / (8.314 * T_AMB)          # 0.00619 kg/m3
C_G = math.sqrt(1.405 * 8.314 / 2.016e-3 * T_AMB)    # ~1,120 m/s
LIFT = 0.0827                                         # kg per m3 (structure Table 2)
FLOOR = 0.10                                          # kg/m2


def tau(d, m, V):
    """psf_sweep 'Floor physics': low-pass time constant of a floating-lid layer.
    tau^2 = a^2 d (d + m/rho_g) + 1/4 [Za a^2 d/rho_g - (m + rho_g d)/Za]^2.
    The ground peak goes as W * (V tau)^-3/2, so equal tau = equal psf for a given aircraft."""
    a2 = 1 / V ** 2 - 1 / C_G ** 2
    Za = RHO_A / math.sqrt(1 / C_A ** 2 - 1 / V ** 2)
    t2 = a2 * d * (d + m / RHO_G) + 0.25 * (Za * a2 * d / RHO_G - (m + RHO_G * d) / Za) ** 2
    return math.sqrt(t2)


def d_for_tau(target, V, lid_of_d):
    lo, hi = 1.0, 400.0
    for _ in range(80):
        mid = 0.5 * (lo + hi)
        if tau(mid, lid_of_d(mid), V) >= target:
            hi = mid
        else:
            lo = mid
    return hi


# Reference requirement: psf_sweep minimum thickness at 20 km with the floating lid m = lift - 0.1
AIRCRAFT = {"M1.7 (Overture)": (1.7 * C_A, 20.0), "M2.0 (Concorde)": (2.0 * C_A, 54.0)}


def lid_reopt():
    rows = []
    for name, (V, d_ref) in AIRCRAFT.items():
        t_req = tau(d_ref, LIFT * d_ref - FLOOR, V)
        # peak check helper: peak ratio vs requirement = (tau_req/tau)^1.5
        cases = [
            ("psf_sweep reference (lid = lift - floor, no hold-down)", 0.0, None),
            ("lift-closed, hold-down 0.62 (bare tethers, rim C_D 1)", 0.62, None),
            ("lift-closed, hold-down 0.25 (rounded rim, C_D 0.3)", 0.25, None),
            ("light lid 0.40 kg/m2, excess lift on tethers", 0.62, 0.40),
            ("light lid 1.0 kg/m2, excess lift on tethers", 0.62, 1.0),
        ]
        for label, hd, fixed in cases:
            if fixed is None:
                f = (lambda d, hd=hd: LIFT * d - FLOOR - hd)
            else:
                f = (lambda d, m=fixed: m)
            d = d_for_tau(t_req, V, f)
            m = f(d)
            excess = LIFT * d - FLOOR - m - hd if fixed is not None else 0.0
            rows.append({"aircraft": name, "case": label, "tau_req_s": round(t_req, 3),
                         "H2_thickness_m": round(d, 1), "lid_kg_m2": round(m, 2),
                         "lift_kg_m2": round(LIFT * d, 2), "holddown_kg_m2": round(hd + excess, 2),
                         "excess_lift_to_tethers_kg_m2": round(excess, 2),
                         "H2_kg_m2": round(RHO_G * d, 3)})
        # sensitivity: at the reference thickness, peak rise if the lid is lightened
        for m in (0.4, 1.0, 2.3):
            if m < LIFT * d_ref - FLOOR:
                r = (t_req / tau(d_ref, m, V)) ** 1.5
                rows.append({"aircraft": name, "case": f"same d={d_ref:.0f} m, lid {m} kg/m2: peak x{r:.2f}",
                             "tau_req_s": round(t_req, 3), "H2_thickness_m": d_ref, "lid_kg_m2": m,
                             "lift_kg_m2": round(LIFT * d_ref, 2), "holddown_kg_m2": "", "excess_lift_to_tethers_kg_m2": "",
                             "H2_kg_m2": round(RHO_G * d_ref, 3)})
    return rows


# Minimum structural mass of the lid (per m2), screening.  [K] handbook-class.
MIN_STRUCT = [
    # item, low, base, high kg/m2, basis
    ("Face film (LLDPE, 0.2 kN/m face + boom 0.3 kN/m transient, 5-10 MPa allowable at 217 K)", 0.03, 0.05, 0.09,
     "t = N/sigma = 25-90 um; structure-concept face ~0.2 kN/m"),
    ("Bonded second ply (fatigue: bonded multi-ply, no free flaps)", 0.02, 0.04, 0.06, "lid-fatigue design C"),
    ("UV/ozone top skin (PVF/PVDF or metallised, 12-25 um)", 0.02, 0.035, 0.05, "lid-fatigue S4; economics UV line"),
    ("Load-spreading tapes / scrim (UHMWPE)", 0.08, 0.10, 0.15, "structure-concept section 2; cost_model tape_kg_m2"),
    ("Seams, tie pads, rolling-lobe hinges (+15-25 %)", 0.02, 0.04, 0.08, "structure-concept +15 %; lid-fatigue pads"),
    ("Controls, sensors, transponders", 0.03, 0.05, 0.08, "structure-concept 0.05"),
]


# ---------------------------------------------------------------------------------------------
# Stratospheric PV resource
# ---------------------------------------------------------------------------------------------
S0 = 1361.0


def insolation_20km(lat_deg, tau_am1=0.97, horizontal=True):
    """Annual horizontal insolation (kWh/m2/yr) at 20 km. Beam only: above 20 km lies ~5 % of the
    air mass; Rayleigh + ozone (Chappuis/Hartley) give ~3 % loss at air mass 1 [K]; SoftBank reports
    20 km irradiance within 1-2 % of AM0 at high sun (search summary). Diffuse sky from above ~0."""
    lat = math.radians(lat_deg)
    tot = 0.0
    for day in range(365):
        dec = math.radians(23.44) * math.sin(2 * math.pi * (284 + day + 1) / 365)
        ecc = 1 + 0.033 * math.cos(2 * math.pi * (day + 1) / 365)
        for k in range(288):  # 5-min steps
            h = math.radians((k + 0.5) * 360 / 288 - 180)
            cz = math.sin(lat) * math.sin(dec) + math.cos(lat) * math.cos(dec) * math.cos(h)
            if cz <= 0.01:
                continue
            am = min(1 / cz, 30)
            tot += S0 * ecc * tau_am1 ** am * cz * (24 / 288)
    return tot / 1000.0


def lid_temperature(absorb_solar, eta_pv, emis_top, emis_bot, S=1285.0, ir_up=250.0, ir_down=20.0, h_conv=2.0):
    """Noon radiative-convective balance of a thin lid at 20 km (both faces convect at h to 217 K air).
    S: midsummer noon horizontal irradiance at ~40 N, 20 km. ir_up: earth+troposphere upwelling IR at 20 km [K]."""
    sig = 5.670e-8
    q_in = absorb_solar * S - eta_pv * S + emis_bot * ir_up + emis_top * ir_down
    lo, hi = 150.0, 500.0
    for _ in range(80):
        T = 0.5 * (lo + hi)
        q_out = (emis_top + emis_bot) * sig * T ** 4 + 2 * h_conv * (T - T_AMB)
        if q_out > q_in:
            hi = T
        else:
            lo = T
    return T


def main():
    rows = lid_reopt()
    with open(HERE / "lid_reopt.csv", "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        w.writeheader()
        w.writerows(rows)
    for r in rows:
        print(r)

    with open(HERE / "min_structure.csv", "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["item", "low_kg_m2", "base_kg_m2", "high_kg_m2", "basis"])
        for it in MIN_STRUCT:
            w.writerow(it)
        w.writerow(["TOTAL lid structure", *[round(sum(x[i] for x in MIN_STRUCT), 3) for i in (1, 2, 3)], ""])
    print("min lid structure kg/m2:", [round(sum(x[i] for x in MIN_STRUCT), 3) for i in (1, 2, 3)])

    pv = []
    for lat in (36.0, 40.0, 42.0):
        H = insolation_20km(lat)
        Htoa = insolation_20km(lat, tau_am1=1.0)
        pv.append({"lat_deg": lat, "H_20km_kWh_m2_yr": round(H), "H_TOA_kWh_m2_yr": round(Htoa)})
    with open(HERE / "pv_resource.csv", "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(pv[0].keys()))
        w.writeheader()
        w.writerows(pv)
    print(pv)

    th = []
    for label, a, eta, et, eb in [
        ("white UV skin (alpha 0.15), no PV", 0.15, 0.0, 0.85, 0.85),
        ("clear/metallised lid (alpha 0.25, low-e under)", 0.25, 0.0, 0.85, 0.10),
        ("thin-film PV lid, high-e underside", 0.90, 0.15, 0.85, 0.85),
        ("thin-film PV lid, low-e underside (decoupled from gas)", 0.90, 0.15, 0.85, 0.10),
    ]:
        T = lid_temperature(a, eta, et, eb)
        th.append({"lid": label, "noon_T_K": round(T, 1), "noon_T_C": round(T - 273.15, 1),
                   "excess_over_air_K": round(T - T_AMB, 1)})
    with open(HERE / "pv_thermal.csv", "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(th[0].keys()))
        w.writeheader()
        w.writerows(th)
    for r in th:
        print(r)


if __name__ == "__main__":
    main()
