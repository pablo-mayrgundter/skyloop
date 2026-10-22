"""Vacuum BoomShield: pumping, leakage, outgassing, puncture and start-up numbers.

Order-of-magnitude screening. Per km of corridor, raft of D = 30 m tubes over W = 250 m
(skin area = pi x planform; volume = pi D/4 x planform). Permeabilities are typical
room-temperature handbook-class values (NOT verified against a primary source here;
see the review); at 216 K they are ~10-100x lower, but sunlit film runs warmer.

Run: python boom-shield/vacuum/vacuum_ops.py   -> prints tables, writes ops_pumping.csv
"""
from __future__ import annotations

import csv
import math
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "sim"))
from boomsim.atmosphere import std_atm  # noqa: E402

A20 = std_atm(20000.0)
P0, RHO0, T0 = float(A20.p), float(A20.rho), float(A20.T)
R_U, R_AIR, GAM = 8.314462618, 287.05, 1.4
BARRER = 3.35e-16          # mol m / (m^2 s Pa)

W, D, LKM = 250.0, 30.0, 1000.0
A_SKIN = math.pi * W * LKM                  # m^2 of skin per km (raft of tubes)
VOL = math.pi * D / 4 * W * LKM             # m^3 per km

# typical O2 / N2 permeabilities at ~23 C [Barrer] (handbook-class, unverified here)
FILMS = {
    "LDPE": (2.9, 1.0),
    "polyimide (Kapton-type)": (0.14, 0.03),
    "PET (biax.)": (0.025, 0.005),
    "PET + metallised/AlOx barrier (~100x)": (2.5e-4, 5e-5),
    "EVOH layer (dry)": (5e-4, 1e-4),
}
T_FILM = 25e-6
# rough-vacuum pump rule of thumb: ~10 kW per m^3/s of speed at 1-50 Pa with a roots
# stack exhausting to 5.5 kPa (incl. backing); ~250 kg per m^3/s of installed speed.
KW_PER_M3S, KG_PER_M3S = 10.0, 250.0


def choked_mdot(area, Cd=0.6):
    return Cd * area * P0 * math.sqrt(GAM / (R_AIR * T0)) * (2 / (GAM + 1)) ** ((GAM + 1) / (2 * (GAM - 1)))


def main():
    rows = []
    print(f"per km: skin {A_SKIN:.3g} m^2, volume {VOL:.3g} m^3, air mass at ambient {RHO0 * VOL / 1e3:.0f} t")
    pO2, pN2 = 0.21 * P0, 0.78 * P0
    print("\nPermeation (25 um film, 23 C values):")
    for name, (PO2, PN2) in FILMS.items():
        n = (PO2 * pO2 + PN2 * pN2) * BARRER / T_FILM * A_SKIN      # mol/s
        Q = n * R_U * T0                                            # Pa m^3/s
        for pr in (1e-3, 1e-2, 1e-1):
            S = Q / (pr * P0)
            rows.append(dict(source=f"permeation {name}", p_ratio=pr, Q_Pa_m3_s=Q, S_m3_s=S,
                             P_kW=S * KW_PER_M3S, pump_kg=S * KG_PER_M3S,
                             P_ideal_W=Q * math.log(1 / pr)))
        print(f"  {name:40s} Q = {Q:9.3g} Pa m^3/s/km -> S = {Q / (1e-3 * P0):8.3g} m^3/s at 1e-3, "
              f"{Q / (1e-1 * P0):8.3g} at 1e-1")
    print("\nHoles (choked inflow), total equivalent hole area per km:")
    for a_cm2 in (0.01, 1.0, 10.0, 100.0):
        md = choked_mdot(a_cm2 * 1e-4)
        Q = md * R_AIR * T0
        for pr in (1e-3, 1e-2, 1e-1):
            S = Q / (pr * P0)
            rows.append(dict(source=f"holes {a_cm2:g} cm2/km", p_ratio=pr, Q_Pa_m3_s=Q, S_m3_s=S,
                             P_kW=S * KW_PER_M3S, pump_kg=S * KG_PER_M3S,
                             P_ideal_W=Q * math.log(1 / pr)))
        print(f"  {a_cm2:7.2f} cm^2: mdot {md * 1e3:8.3g} g/s, S = {Q / (1e-3 * P0):8.3g} m^3/s at 1e-3 "
              f"({Q / (1e-3 * P0) * KW_PER_M3S:8.3g} kW, {Q / (1e-3 * P0) * KG_PER_M3S / 1e3:7.3g} t pumps), "
              f"{Q / (1e-1 * P0):8.3g} m^3/s at 1e-1")
    q_out = 1e-6  # Pa m^3/(s m^2), polymer film after ~100 h at room T (unbaked)
    print(f"\nOutgassing at {q_out:g} Pa m^3/s/m^2: Q = {q_out * A_SKIN:.3g} -> "
          f"S = {q_out * A_SKIN / (1e-3 * P0):.3g} m^3/s at 1e-3")
    lift_t = RHO0 * VOL / 1e3
    print(f"\nLift per km (gross, raft D30, W250): {lift_t:.0f} t; pump mass budget for 1 m^3/s: {KG_PER_M3S} kg")
    print("\nPump-down from ambient at altitude, t = V/S ln(1/p_ratio):")
    for S in (1.0, 10.0, 100.0):
        for pr in (1e-1, 1e-3):
            print(f"  S = {S:5.0f} m^3/s/km, to {pr:g}: {VOL / S * math.log(1 / pr) / 86400:7.1f} days")
    print("\nPuncture: time to lose 10 % of lift for a 30 m x 250 m module (one tube)")
    vmod = math.pi * D ** 2 / 4 * 250
    for a in (1e-4, 1e-2, 1.0, 10.0):
        md = choked_mdot(a)
        print(f"  hole {a:7.4f} m^2: {0.1 * RHO0 * vmod / md:9.3g} s")
    print(f"\nStored work p V of one module: {P0 * vmod / 1e9:.2f} GJ (~{P0 * vmod / 4.184e9:.2f} t TNT eq.)")
    print("\nLift fraction (1 - p_ratio) vs He 0.86 / H2 0.93:")
    for pr in (0.2, 0.14, 0.1, 0.07, 0.05, 1e-2, 1e-3):
        print(f"  p_ratio {pr:6g}: lift fraction {1 - pr:.3f}, p_int {pr * P0:7.1f} Pa")
    print("\nConstant-dp ascent: lift density at dp = 5.5 kPa is dp/(R T):")
    for z in (0, 5000, 10000, 15000, 20000):
        a = std_atm(float(z))
        print(f"  z {z / 1e3:4.0f} km: p_ext {float(a.p):7.0f} Pa, p_int {float(a.p) - P0:7.0f} Pa, "
              f"lift {P0 / (R_AIR * float(a.T)):.4f} kg/m^3 ({P0 / (R_AIR * float(a.T)) / RHO0:.2f} of 20 km)")
    with open(os.path.join(HERE, "ops_pumping.csv"), "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]))
        w.writeheader()
        for r in rows:
            w.writerow({k: (f"{v:.4g}" if isinstance(v, float) else v) for k, v in r.items()})


if __name__ == "__main__":
    main()
