#!/usr/bin/env python3
"""Higher-Mach study: parametric supersonic airliner across Mach 1.7-4.0 (screening, [U]).

Run:  python3 hm_aircraft.py   -> aircraft_param.csv, concepts.csv

Model (every coefficient is an engineering estimate unless cited in evidence.csv):
  cruise L/D      = k_LD * 4 (M + 3) / M                 Kuchemann barrier (E13), k_LD = 0.85 (modern design)
  overall eff.    eta0(M) tabulated per engine family     TSFC = V / (eta0 * LHV)
  range factor    RF = eta0 * LHV / g * L/D               Breguet
  accel fuel      (V^2/2 + g h) / (eta_acc * LHV) of W    energy to reach cruise speed and altitude
  sizing          W0 = payload / (1 - f_empty(M) - f_fuel)  closes on the design mission
  ramp accel      a = g (k_T - 1) / (L/D) * energy-height factor   k_T = max thrust / cruise drag
Calibration: the M1.7 fixed-cycle point reproduces Overture (65 seats, 4,250 nmi, ~77 t MTOW) and the
repo's 0.23 m/s^2 ramp fit (trunk.md).
"""
from __future__ import annotations

import math
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
G = 9.81
LHV = 43.2e6                 # J/kg, Jet-A / JP-class kerosene
A_SOUND = 295.1              # m/s, isothermal stratosphere (repo convention)
SEATS = 65                   # repo convention (cost_model.SEATS, Overture-class)
M_PAX = 100.0                # kg per seat: passenger + bags
K_LD = 0.85                  # achieved / Kuchemann max L/D (Concorde ~7.3 at M2 = 0.73; HSCT targets ~0.9)
RESERVE_F = 0.06             # reserve + contingency fuel fraction of W0 [U]
FUEL_KG_SEAT_KM_SUB = 0.090 / 3.16 * 0.80   # subsonic narrowbody: 90 g CO2/pax-km (demand.md) at LF 0.8 -> per seat-km
DESIGN_RANGE_KM = 5000.0     # domestic design mission (NYC-LAX flown + feeders + diversion), excl. reserves

# overall propulsive efficiency eta0 = V / (TSFC * LHV), by engine family (installed, cruise) [U]
#   fixed-cycle medium-bypass (Overture/Symphony class, no afterburner): ~0.40 at M1.7 (Concorde Olympus 0.41 at M2)
#   adaptive / variable cycle: third stream CLOSED in supercruise -> low-bypass turbojet-like cycle; credit
#     only ~3 % SFC at supersonic cruise (E5), the 25 % (E1, E2) applies at subsonic/loiter points
#   TBCC above ~M3: turbine to ~M2.5-3 then ramjet (E12); ramjet eta0 ~0.45-0.5 at M3.5-4
ETA0_FIXED = {1.7: 0.400, 2.0: 0.415, 2.2: 0.425, 2.5: 0.435, 3.0: 0.440, 3.5: 0.430, 4.0: 0.420}
ADAPT_SFC_CREDIT = 0.03      # supersonic-cruise SFC credit of an adaptive cycle (E5); subsonic legs 0.15-0.25
ETA0_TBCC = {3.0: 0.440, 3.5: 0.460, 4.0: 0.475}

# empty weight = f_empty(M) * W0 + W_FIXED_SEAT * seats (furnishing, systems, crew: does not shrink with W0).
# Calibrated so the M1.7 fixed-cycle point at Overture's own 7,870 km (4,250 nmi) gives ~77 t MTOW.
# Steps for titanium (> M2.25) and TBCC (> M3.05) structure [U].
W_FIXED_SEAT = 250.0         # kg per seat [U]
F_EMPTY_17 = None            # set by calibrate() below


def f_empty(M, engine="adaptive"):
    f = F_EMPTY_17 + 0.02 * (M - 1.7)
    if M > 2.25:
        f += 0.030        # titanium / high-temperature structure, hot-fuel thermal management (E11)
    if M > 3.05:
        f += 0.050        # dual flowpath TBCC, insulation, ramjet ducts
    if engine == "adaptive":
        f += 0.005        # third-stream ducting, valves
    return f


def interp(tab, M):
    ks = sorted(tab)
    return float(np.interp(M, ks, [tab[k] for k in ks]))


def cruise_alt_km(M, q_target=25e3):
    """Preferred cruise altitude: dynamic pressure ~25 kPa (Concorde ~25 kPa at M2/16.8 km; SR-71 ~21 kPa),
    and not below 20.25 km (shield top 20 km + 250 m), capped at 26 km."""
    from boomsim_atm import pressure_alt
    p = q_target / (0.7 * M ** 2)
    return float(np.clip(pressure_alt(p) / 1e3, 20.25, 26.0))


def ld(M, k=K_LD):
    return k * 4 * (M + 3) / M


def eta0(M, engine):
    if engine == "tbcc":
        return interp(ETA0_TBCC, max(M, 3.0))
    e = interp(ETA0_FIXED, M)
    return e / (1 - ADAPT_SFC_CREDIT) if engine == "adaptive" else e


def tsfc_lb(M, engine):
    V = M * A_SOUND
    return V / (eta0(M, engine) * LHV) * 3600 * G * 1.0    # kg/(N s) -> lb/(lbf h): x 3600 x g


def accel_frac(M, h_km, engine):
    """Fuel fraction to climb and accelerate from rest to cruise (energy / (eta_acc LHV)).
    eta_acc = 0.6 eta0: thrust lapse, transonic drag rise, partial afterburner [U]."""
    V = M * A_SOUND
    E = 0.5 * V ** 2 + G * h_km * 1e3
    return E / (0.6 * eta0(M, engine) * LHV)


def size(M, engine="adaptive", seats=SEATS, R_km=DESIGN_RANGE_KM, k_ld=K_LD, kT=1.40):
    """Close W0 on the design mission. Returns dict or None if it does not close."""
    h = cruise_alt_km(M)
    L_D = ld(M, k_ld)
    e0 = eta0(M, engine)
    RF = e0 * LHV / G * L_D / 1e3                     # km
    fa = accel_frac(M, h, engine)
    fc = (1 - fa) * (1 - math.exp(-R_km / RF))
    ff = fa + fc + RESERVE_F
    fe = f_empty(M, engine) + (0.02 * (kT - 1.40) / 0.2 if kT > 1.40 else 0.0)   # bigger engines for thrust margin
    fp = 1 - fe - ff
    if fp <= 0.10:
        return None
    W0 = seats * (M_PAX + W_FIXED_SEAT) / fp
    W_mid = W0 * (1 - fa - 0.5 * fc)
    return dict(M=M, engine=engine, h_cruise_km=h, LD=L_D, eta0=e0, tsfc_lb_lbf_h=tsfc_lb(M, engine),
                RF_km=RF, f_accel=fa, f_cruise=fc, f_fuel=ff, f_empty=fe, f_payload=fp, W0_t=W0 / 1e3,
                W_mid_t=W_mid / 1e3, t_per_seat=W0 / 1e3 / seats, kT=kT)


def trip_fuel(ac, flown_km):
    """Block fuel (kg) for a trip of flown_km: accel/climb energy + Breguet cruise on the rest."""
    W0 = ac["W0_t"] * 1e3
    # aircraft at trip weight: take-off weight ~ W0 minus unused design fuel
    fa = ac["f_accel"]
    fc = (1 - fa) * (1 - math.exp(-flown_km / ac["RF_km"]))
    return W0 * (fa + fc) * 0.97        # 3 %: TOW below MTOW on a shorter-than-design trip [U]


def ramp_accel(ac, kT=None, h0_km=11.0):
    """Mean along-track acceleration on the ramp (m/s^2): a = g (k_T - 1)/(L/D), reduced by the extra
    climb from ~11 km to cruise altitude done during the ramp (energy height).  kT = 1.22 reproduces the
    repo's 0.23 m/s^2 Overture fit at M1.7 (trunk.md)."""
    kT = ac["kT"] if kT is None else kT
    M = ac["M"]
    V, v0 = M * A_SOUND, 0.95 * A_SOUND
    a_raw = G * (kT - 1) / ac["LD"]
    dKE = 0.5 * (V ** 2 - v0 ** 2)
    dPE = G * max(ac["h_cruise_km"] - 18.0, 0) * 1e3     # the repo's 0.23 fit already includes the climb to ~18 km
    return a_raw * dKE / (dKE + dPE)


def calibrate(W0_ref=77.0e3, R_ref=7870.0):
    """Solve F_EMPTY_17 so the fixed-cycle M1.7 aircraft at Overture's range weighs W0_ref."""
    global F_EMPTY_17
    F_EMPTY_17 = 0.0
    ac = size(1.7, "fixed", R_km=R_ref, kT=1.22)
    F_EMPTY_17 = 1 - ac["f_fuel"] - SEATS * (M_PAX + W_FIXED_SEAT) / W0_ref
    return F_EMPTY_17


calibrate()

CONCEPTS = {
    # name: Mach, engine, thrust ratio kT, structure, note
    "Overture M1.7 (ref)": (1.7, "fixed", 1.22, "aluminium/composite", "repo reference; Symphony-class (E7)"),
    "A: M2.2 Al-limit, adaptive": (2.2, "adaptive", 1.40, "aluminium (130 C skin, E11)", "Concorde-like structure limit"),
    "B: M3.0 Ti, adaptive": (3.0, "adaptive", 1.40, "titanium (~290 C avg, E11)", "XA10x-class turbine to M3, dry/min AB"),
    "C: M3.5 TBCC": (3.5, "tbcc", 1.60, "titanium + insulation", "turbine to ~M2.8, ramjet above (E12)"),
    "D: M4.0 TBCC": (4.0, "tbcc", 1.60, "titanium/Inconel hot structure", "beyond ambient-H2 evanescent ceiling"),
}


def concept_table():
    rows = []
    for name, (M, eng, kT, struct, note) in CONCEPTS.items():
        ac = size(M, eng, kT=kT, R_km=7870.0 if name.startswith("Overture") else DESIGN_RANGE_KM)
        if ac is None:
            rows.append(dict(concept=name, M=M, engine=eng, closes=False))
            continue
        a = ramp_accel(ac)
        if name.startswith("Overture"):
            a = 0.23                                     # repo fit exactly
            ac["W_mid_t"] = 62.0                         # psf_sweep Overture preset (mid-cruise)
        ac.update(concept=name, structure=struct, note=note, closes=True, a_acc=a, a_dec=min(2 * a, 0.9),
                  fuel_kg_seat_km_nyclax=trip_fuel(ac, 4093) / (SEATS * 4093))
        ac["fuel_vs_subsonic_seat"] = ac["fuel_kg_seat_km_nyclax"] / FUEL_KG_SEAT_KM_SUB
        ac["q_kPa_at_20km"] = 0.7 * 5529 * M ** 2 / 1e3
        rows.append(ac)
    return pd.DataFrame(rows)


def sweep():
    rows = []
    for M in np.round(np.arange(1.7, 4.01, 0.1), 2):
        eng = "fixed" if M < 1.75 else ("adaptive" if M <= 3.05 else "tbcc")
        for kT in (1.22, 1.40, 1.60):
            ac = size(M, eng, kT=kT)
            if ac is None:
                rows.append(dict(M=M, engine=eng, kT=kT, closes=False))
                continue
            ac.update(closes=True, a_acc=ramp_accel(ac), fuel_kg_seat_km=trip_fuel(ac, 4093) / (SEATS * 4093))
            ac["fuel_vs_subsonic_seat"] = ac["fuel_kg_seat_km"] / FUEL_KG_SEAT_KM_SUB
            ac["q_kPa_cruise"] = 25.0 if ac["h_cruise_km"] > 20.3 else 0.7 * 5529 * M ** 2 / 1e3
            ac["stag_T_C"] = 216.65 * (1 + 0.2 * 0.89 * M ** 2) - 273.15   # recovery factor 0.89
            rows.append(ac)
    return pd.DataFrame(rows)


if __name__ == "__main__":
    S = sweep()
    S.to_csv(HERE / "aircraft_param.csv", index=False, float_format="%.4g")
    C = concept_table()
    C.to_csv(HERE / "concepts.csv", index=False, float_format="%.4g")
    pd.set_option("display.width", 250)
    print(S[S.kT == 1.40][["M", "engine", "h_cruise_km", "LD", "eta0", "tsfc_lb_lbf_h", "f_fuel", "f_empty",
                          "f_payload", "W0_t", "t_per_seat", "a_acc", "fuel_kg_seat_km", "fuel_vs_subsonic_seat",
                          "stag_T_C"]].round(3).to_string())
    print(C[["concept", "M", "W0_t", "W_mid_t", "LD", "tsfc_lb_lbf_h", "a_acc", "a_dec", "h_cruise_km",
             "fuel_kg_seat_km_nyclax", "fuel_vs_subsonic_seat"]].round(3).to_string())
