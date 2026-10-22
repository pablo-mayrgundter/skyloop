"""BoomShield structural concept study: lift, mass, superpressure, tie stiffness, compression,
wind, tether and hydrogen-inventory screening numbers (numpy only).

Run:  python3 budget.py   ->  results/budget_*.csv and printed tables.
All values per m^2 of PLANFORM unless stated. 20 km US-1976 (boomsim.atmosphere).
Assumed material values are marked ASSUMED; they are typical handbook-class numbers,
not verified product data.
"""
import csv
import math
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.abspath(os.path.join(HERE, "..", "sim")))
from boomsim.atmosphere import std_atm  # noqa: E402

OUT = os.path.join(HERE, "results")
os.makedirs(OUT, exist_ok=True)

A = std_atm(20000.0)
RHO, P, T, C = float(A.rho), float(A.p), float(A.T), float(A.c)
G = 9.745                    # m/s^2 at 20 km
RU = 8.314462618
M_H2, M_HE = 2.016e-3, 4.0026e-3
GAM_H2 = 1.405


def rho_gas(M, dp=0.0):
    return (P + dp) * M / (RU * T)


RHO_H2 = rho_gas(M_H2)
DRHO = RHO - RHO_H2           # lift density, kg/m^3
V = 2.0 * C                   # Mach 2 trace speed at 20 km

# ---- ASSUMED materials -----------------------------------------------------------------
S_LAM = 5.0e5 / 4.0   # envelope laminate: ~1e5 N/m ultimate at ~0.2 kg/m^2 (HAA class), SF 4 -> N/m per kg/m^2
M_MIN = 0.03          # minimum-gauge gas-barrier film, kg/m^2 (~25-30 um PE/EVOH multilayer)
SEAM = 1.15           # seams, joints, local doublers
S_FIB = 3.0e9 / 970.0 / 3.0   # UHMWPE tendon, SF 3: J/kg = (N/m^2)/(kg/m^3)
E_FIB, SIG_FIB_W = 110e9, 1.0e9   # tendon modulus and working stress (strain ~0.9 %)
RHO_CF, E_CF = 1600.0, 120e9      # CFRP strut (with knock-down)

# ---- operating assumptions --------------------------------------------------------------
SUPERHEAT = {"low (coated)": 10.0, "central": 20.0, "high (clear film, measured He SPB)": 39.0}
DP_NIGHT = 150.0      # Pa, minimum night superpressure (shape + margin)
RESERVE = 0.15        # fraction of gross lift held as reserve (leak, ageing, diurnal, tether pull)
CTRL = 0.05           # controls, valves, sensors, power, kg/m^2
TETHER_SPREAD = 0.10  # load-spreading tendons / patches for tether loads, kg/m^2


def dp_swing(dT):
    """Constant-volume superpressure rise for a gas superheat dT (K)."""
    return P * dT / T


def skin_mass(N):
    """Areal mass of one skin carrying tension N (N/m)."""
    return SEAM * max(N / S_LAM, M_MIN)


def rows_to_csv(name, rows):
    with open(os.path.join(OUT, name), "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
        w.writeheader()
        for r in rows:
            w.writerow(r)


def fmt(rows, keys):
    print(" | ".join(keys))
    for r in rows:
        print(" | ".join(f"{r[k]:.3g}" if isinstance(r[k], float) else str(r[k]) for k in keys))
    print()


def main():
    print(f"20 km: rho_air {RHO:.4f}, rho_H2 {RHO_H2:.5f}, rho_He {rho_gas(M_HE):.5f} kg/m3, "
          f"p {P:.0f} Pa, T {T:.2f} K, lift density H2 {DRHO:.4f} (He {RHO - rho_gas(M_HE):.4f}) kg/m3")
    ag2 = 1 / V ** 2 - GAM_H2 * 0 - 1 / math.sqrt(GAM_H2 * RU * T / M_H2) ** 2
    print(f"evanescent gas must have rho_g <= gamma p / V^2 = {1.4 * P / V**2:.4f} kg/m3 at M2 "
          f"(air {RHO:.4f}): any evanescent filling lifts >= {(RHO - 1.4 * P / V**2) / RHO:.0%} of air")
    for k, dT in SUPERHEAT.items():
        print(f"superheat {k}: {dT:.0f} K -> constant-volume dp swing {dp_swing(dT):.0f} Pa")
    print()

    dT = SUPERHEAT["central"]
    dpmax = DP_NIGHT + dp_swing(dT)
    rows = []
    # ------------------------------------------------------------------ (a) raft of tubes
    for D in (20.0, 30.0, 60.0):
        R = D / 2
        lift = DRHO * math.pi * D / 4          # circle/slab factor pi/4
        for dpm, tag in ((dpmax, "central"), (DP_NIGHT + dp_swing(SUPERHEAT["high (clear film, measured He SPB)"]), "hot")):
            skin = math.pi * skin_mass(dpm * R)      # perimeter pi D over width D
            lacing = 0.01
            mass = skin + lacing + CTRL + TETHER_SPREAD
            K_ov = 3 * dpm / R
            rows.append(dict(arch="(a) raft", size=f"D {D:.0f}", case=tag, lift=lift, skin=skin,
                             ties_or_struts=lacing, other=CTRL + TETHER_SPREAD, mass=mass,
                             reserve=RESERVE * lift, net=lift * (1 - RESERVE) - mass,
                             dp_max=dpm, K_face=K_ov, note="K = 3 dp/R ovalisation (n=2)"))
    # ------------------------------------------------------------------ (b) drop-stitch SP pillow
    for Tt in (30.0, 60.0):
        lift = DRHO * Tt
        Nf = dpmax * Tt / 2
        faces = 2 * skin_mass(Nf)
        ties = 970.0 * Tt * dpmax / SIG_FIB_W * 3.0     # SF 3 on tie area
        mass = faces + ties + CTRL + TETHER_SPREAD
        K = dpmax / ((SIG_FIB_W / E_FIB) * Tt)           # stiff fibre ties at working strain
        rows.append(dict(arch="(b) drop-stitch SP", size=f"T {Tt:.0f}", case="fibre ties", lift=lift, skin=faces,
                         ties_or_struts=ties, other=CTRL + TETHER_SPREAD, mass=mass, reserve=RESERVE * lift,
                         net=lift * (1 - RESERVE) - mass, dp_max=dpmax, K_face=K,
                         note="K = dp/(eps_w T), independent of tie spacing"))
        eps_need = dpmax / (30.0 * Tt)
        rows.append(dict(arch="(b) drop-stitch SP", size=f"T {Tt:.0f}", case="soft ties K=30", lift=lift, skin=faces,
                         ties_or_struts=ties, other=CTRL + TETHER_SPREAD, mass=mass, reserve=RESERVE * lift,
                         net=lift * (1 - RESERVE) - mass, dp_max=dpmax, K_face=30.0,
                         note=f"needs tie working strain {eps_need:.0%}; static stretch {dpmax/30:.0f} m, "
                              f"day-night stroke {dp_swing(dT)/30:.0f} m (> T): not constant-volume"))
    # ------------------------------------------------------------------ (c) lens + cross-width struts
    W = 1000.0
    for Tt in (30.0, 60.0):
        lift = DRHO * math.pi / 4 * Tt                 # elliptic section mean depth (pi/4) T
        Rc = W ** 2 / (4 * Tt)                         # each face sagitta T/2
        for dpm, tag in ((DRHO * G * Tt, "zero-pressure (dp = lift head)"), (dpmax, "superpressure")):
            N = dpm * Rc
            faces = 2 * skin_mass(N)
            Cp = 2 * N                                 # rim compression per m of corridor
            best = None
            for s in (250.0, 500.0, 1000.0, 2000.0):   # strut spacing along corridor
                phi = 200.0                            # r/t
                Pst = Cp * s
                Ast = W * math.sqrt(4 * Pst / (math.pi * E_CF * phi))
                m_st = RHO_CF * Ast * W / (s * W)
                r_st = math.sqrt(Ast * phi / (2 * math.pi))
                if best is None or m_st < best[0]:
                    best = (m_st, s, r_st, Pst / Ast / 1e6)
            mass = faces + best[0] + CTRL + TETHER_SPREAD
            rows.append(dict(arch="(c) lens + struts", size=f"T {Tt:.0f}", case=tag, lift=lift, skin=faces,
                             ties_or_struts=best[0], other=CTRL + TETHER_SPREAD, mass=mass, reserve=RESERVE * lift,
                             net=lift * (1 - RESERVE) - mass, dp_max=dpm, K_face=0.0,
                             note=f"face N {N/1e3:.0f} kN/m; rim compression {Cp/1e3:.0f} kN/m; best strut "
                                  f"spacing {best[1]:.0f} m, tube radius {best[2]:.1f} m, stress {best[3]:.0f} MPa"))
    # ------------------------------------------------------------------ (d) floating-lid ZP pillow
    for Tt in (20.0, 30.0, 60.0):
        lift = DRHO * Tt * 0.99                        # ~1 % rim rounding at W = 1 km
        Nf = DRHO * G * Tt ** 2 / 4                    # rim bulge tension shared by lid and floor
        lid, floor = skin_mass(Nf) + 0.03, skin_mass(Nf)   # lid +0.03 UV/thermal-control layer
        walls = 0.01                                   # side walls, limp along-track partitions, bulkheads
        ties = 0.02                                    # soft tie elements (springs) + lines, K 3-10 Pa/m
        mass = lid + floor + walls + ties + CTRL + TETHER_SPREAD
        rows.append(dict(arch="(d) floating-lid ZP", size=f"T {Tt:.0f}", case="soft ties K 3-10", lift=lift,
                         skin=lid + floor, ties_or_struts=walls + ties, other=CTRL + TETHER_SPREAD, mass=mass,
                         reserve=RESERVE * lift, net=lift * (1 - RESERVE) - mass, dp_max=DRHO * G * Tt,
                         K_face=10.0,
                         note=f"face N {Nf:.0f} N/m; RT needs K > {DRHO*G:.2f} Pa/m; diurnal lift swing "
                              f"+-{dT/T:.0%} = {lift*dT/T:.2f} kg/m2"))
    rows_to_csv("budget_architectures.csv", rows)
    fmt(rows, ["arch", "size", "case", "lift", "skin", "ties_or_struts", "mass", "net", "dp_max", "K_face"])
    for r in rows:
        print(f"  {r['arch']} {r['size']} {r['case']}: {r['note']}")
    print()

    # ------------------------------------------------------------------ wind
    wrows = []
    for U in (10.0, 20.0, 30.0, 40.0, 50.0):
        q = 0.5 * RHO * U ** 2
        for Tt in (30.0, 60.0):
            side = q * (0.8 * Tt + 2 * 0.002 * 1000.0)            # N per m of corridor (crosswind)
            wg = 0.05 * U                                           # gust/wave w' ~ 5 % of U (ASSUMED)
            lift_p = q * 2 * math.pi * wg / U                       # thin-airfoil pressure on 1 km chord
            wrows.append(dict(U=U, T=Tt, q=q, side_N_per_m=side, gust_lift_Pa=lift_p,
                              gust_lift_frac_of_H2_lift=lift_p / (DRHO * G * Tt)))
    rows_to_csv("budget_wind.csv", wrows)
    fmt(wrows, list(wrows[0].keys()))

    # ------------------------------------------------------------------ tethers (station every L km)
    trows = []
    for U, F in ((20.0, None), (30.0, None), (40.0, None)):
        F = 0.5 * RHO * U ** 2 * (0.8 * 30.0 + 4.0)
        for L in (2.0, 5.0, 10.0):
            Ht = F * L * 1e3
            elev = math.radians(45.0)
            Ten = Ht / math.cos(elev)
            A_t = Ten / SIG_FIB_W
            d = math.sqrt(4 * A_t / math.pi)
            m_t = 970.0 * A_t * 20e3 / math.sin(elev)
            down = Ten * math.sin(elev) / G + m_t                       # kg
            jet = 490e3 * d / 0.05                                       # addendum: 350-630 kN for 5 cm
            trows.append(dict(U=U, spacing_km=L, horiz_MN=Ht / 1e6, tether_d_cm=100 * d, tether_t=m_t / 1e3,
                              downforce_kg_per_m2=down / (L * 1e3 * 1000.0), jet_drag_MN=jet / 1e6))
    rows_to_csv("budget_tethers.csv", trows)
    fmt(trows, list(trows[0].keys()))

    # ------------------------------------------------------------------ inventory
    irows = []
    L = 4.0e6
    for W_, Tt, f_sec in ((1000.0, 30.0, 1.0), (1000.0, 30.0, math.pi / 4), (1000.0, 60.0, 1.0), (1000.0, 20.0, 1.0)):
        Vol = L * W_ * Tt * f_sec
        Vstp = Vol * (P / 101325.0) * (273.15 / T)
        m_h2 = Vol * RHO_H2
        irows.append(dict(W=W_, T=Tt, section_factor=f_sec, V_alt_m3=Vol, V_stp_m3=Vstp, H2_Mt=m_h2 / 1e9,
                          H2_frac_world_yr=m_h2 / 97e9, He_years_world=Vstp / 185e6,
                          H2_energy_PJ=m_h2 * 120e6 / 1e15))
    rows_to_csv("budget_inventory.csv", irows)
    fmt(irows, list(irows[0].keys()))
    # permeation (room-temperature LDPE ~ 7 Barrer ASSUMED; cold factor exp(-Ea/R (1/T - 1/296)), Ea 35 kJ/mol ASSUMED)
    barrer = 3.35e-16
    for name, Pb, l in (("LDPE 25 um, 23C value", 7.0, 25e-6), ("LDPE 25 um at 217 K", 7.0 * math.exp(-35e3 / RU * (1 / T - 1 / 296)), 25e-6),
                        ("EVOH-barrier laminate (~0.05 Barrer, ASSUMED)", 0.05, 5e-6)):
        J = Pb * barrer * P / l * 3.156e7 * M_H2      # kg / m^2 skin / yr
        frac = 2 * J / (RHO_H2 * 30.0)                 # two faces over a 30 m slab
        print(f"permeation {name}: {J*1e3:.3g} g/m2/yr per skin -> {frac:.2%} of a 30 m layer per yr")


if __name__ == "__main__":
    main()
