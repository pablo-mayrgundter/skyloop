#!/usr/bin/env python3
"""Higher-Mach study: shield requirements above Mach 2 (screening).

Run:  python3 hm_shield.py   -> ceiling.csv, shield_req.csv, shield_alt_trade.csv, shield_vs_mach.csv

Physics is psf_sweep.md's lift-swell model, re-implemented here (validated against psf_sweep.csv):
  tau^2 = a^2 d (d + m/rho_g) + 1/4 [Za a^2 d/rho_g - (m + rho_g d)/Za]^2,  a^2 = 1/V^2 - 1/c_g^2
  p_ground = C R(H) beta W / (2 pi sqrt(2 beta h) (V tau)^1.5)           (C = 1.51 Concorde M2, 1.40 M2.5,
                                                                           1.87 SR-71 shape; R incl. x1.9)
Evanescent only while V < c_g: the ceiling is V = c_g(T_gas, purity).  Warmer gas raises c_g (sqrt T)
and lowers rho_g; air contamination lowers c_g sharply (review_evanescent.md: 5 % air -> 867 m/s).
Lid is lift-closed: m = lift*d - floor 0.1 - hold-down 0.62 kg/m2 (cost-down / operations), >= 0.32.
"""
from __future__ import annotations

import math
from pathlib import Path

import numpy as np
import pandas as pd

from boomsim_atm import atm, ray_gain

HERE = Path(__file__).resolve().parent
RU = 8.314
M_H2, M_AIR = 2.016e-3, 28.96e-3
G_H2, G_AIR = 1.405, 1.400
LIMIT_PA = 0.11 * 47.88          # 5.27 Pa
H_AC = 250.0
FLOOR, HOLDDOWN, MIN_LID = 0.10, 0.62, 0.32
C_FIT = {"concorde": 1.51, "concorde_M2.5": 1.40, "sr71": 1.87}
A_SOUND = 295.1


def gas(T_g, x_air, p):
    Mm = (1 - x_air) * M_H2 + x_air * M_AIR
    # mixture gamma from cv weighting (both ~diatomic: ~1.40)
    cv = (1 - x_air) * RU / (G_H2 - 1) + x_air * RU / (G_AIR - 1)
    gam = (cv + RU) / cv
    return dict(c=math.sqrt(gam * RU * T_g / Mm), rho=p * Mm / (RU * T_g))


def tau(d, m, V, H, T_g=None, x_air=0.0):
    a_ = atm(H)
    g = gas(T_g or a_["T"], x_air, a_["p"])
    a2 = 1 / V ** 2 - 1 / g["c"] ** 2
    if a2 <= 0:
        return 0.0                                   # propagating: no evanescent low-pass
    Za = a_["rho"] / math.sqrt(1 / a_["c"] ** 2 - 1 / V ** 2)
    t2 = a2 * d * (d + m / g["rho"]) + 0.25 * (Za * a2 * d / g["rho"] - (m + g["rho"] * d) / Za) ** 2
    return math.sqrt(t2)


def lift_per_m(H, T_g=None, x_air=0.0):
    a_ = atm(H)
    return a_["rho"] - gas(T_g or a_["T"], x_air, a_["p"])["rho"]


def lid(d, H, T_g=None, x_air=0.0, holddown=HOLDDOWN):
    return max(lift_per_m(H, T_g, x_air) * d - FLOOR - holddown, MIN_LID)


def peak(d, W_N, M, H, T_g=None, x_air=0.0, C=C_FIT["concorde"], holddown=HOLDDOWN, V=None):
    V = V or M * A_SOUND
    b = math.sqrt(M ** 2 - 1)
    t = tau(d, lid(d, H, T_g, x_air, holddown), V, H, T_g, x_air)
    if t <= 0:
        return float("inf")
    return C * ray_gain(H, V) * b * W_N / (2 * math.pi * math.sqrt(2 * b * H_AC) * (V * t) ** 1.5)


def d_min(W_N, M, H, T_g=None, x_air=0.0, C=C_FIT["concorde"], holddown=HOLDDOWN, dmax=2000.0):
    """Minimum H2 thickness (m) for <= 0.11 psf undertrack, and whether lift (not psf) sets it."""
    if peak(dmax, W_N, M, H, T_g, x_air, C, holddown) > LIMIT_PA:
        return float("inf"), "psf"
    lo, hi = 1.0, dmax
    for _ in range(60):
        mid = 0.5 * (lo + hi)
        if peak(mid, W_N, M, H, T_g, x_air, C, holddown) <= LIMIT_PA:
            hi = mid
        else:
            lo = mid
    d_lift = (MIN_LID + FLOOR + holddown) / lift_per_m(H, T_g, x_air)
    return (hi, "psf") if hi >= d_lift else (d_lift, "lift")


def ceiling_mach(H, T_g, x_air):
    a_ = atm(H)
    return gas(T_g, x_air, a_["p"])["c"] / A_SOUND


def ymax_at_shield(M):
    """Max |y| at the shield of ground-reaching rays (onramp_shield.md Table 1: 320 m at M1.7, 355 m at M2);
    extrapolated +0.25 h per Mach above 2 [U]."""
    return float(np.interp(M, [1.16, 1.2, 1.3, 1.5, 1.7, 2.0], [48, 121, 201, 278, 320, 355])) if M <= 2.0 \
        else 355.0 + 0.25 * H_AC * (M - 2.0)


def width_km(M, d, W_N, H, T_g=None, x_air=0.0):
    """Width scaled as onramp_shield.md: 2 x (max|y| + edge margin), margin 450 m x lambda_perp / ref,
    lambda_perp = V tau / beta; ref = Concorde M2 60 m H2 at 20 km; then x 1.75 km / ref width (cost base)."""
    V = M * A_SOUND
    b = math.sqrt(M ** 2 - 1)
    t = tau(d, lid(d, H, T_g, x_air), V, H, T_g, x_air)
    lam = V * t / b
    Vr = 2.0 * A_SOUND
    lam_ref = Vr * tau(60.0, lid(60.0, 20e3), Vr, 20e3) / math.sqrt(3)
    w = 2 * (ymax_at_shield(M) + 450.0 * lam / lam_ref)
    w_ref = 2 * (355.0 + 450.0)
    return 1.75 * w / w_ref


def validate():
    """Reproduce psf_sweep's minima (lid = lift - floor, no hold-down)."""
    rows = []
    for name, W, M, C, ref in (("Concorde M2.0", 1.47e6, 2.0, 1.51, 54), ("Concorde M2.5", 1.47e6, 2.5, 1.40, 77),
                               ("SR-71 M3.2", 0.54e6, 3.2, 1.87, 71), ("Overture M1.7", 62e3 * 9.81, 1.7, 1.60, 20)):
        d, _ = d_min(W, M, 20e3, C=C, holddown=0.0)
        rows.append(dict(case=name, psf_sweep_m=ref, this_model_m=round(d, 1)))
    return pd.DataFrame(rows)


SUPERHEAT_U = 1.5          # W/m2K heat loss of the warm layer, both faces (convection ~0.3-0.6 + low-e radiation) [U]
SPINE_AREA_M2 = 3614e3 * 1.75e3
US_ELEC_TWH = 4200.0


def heating(dT):
    w = SUPERHEAT_U * dT
    twh = w * SPINE_AREA_M2 * 8760 / 1e12
    return w, twh


if __name__ == "__main__":
    pd.set_option("display.width", 250)
    V = validate()
    V.to_csv(HERE / "shield_validation.csv", index=False)
    print(V.to_string())

    # ---------------- ceiling vs gas temperature and purity
    cr = []
    for H in (20e3, 22e3, 24e3, 25e3):
        a_ = atm(H)
        for T_g in (a_["T"], 240.0, 260.0, 280.0, 300.0, 330.0, 360.0):
            for x in (0.0, 0.01, 0.02, 0.05):
                w, twh = heating(T_g - a_["T"])
                cr.append(dict(H_km=H / 1e3, T_amb=round(a_["T"], 2), T_gas=round(T_g, 1), superheat_K=round(T_g - a_["T"], 1),
                               air_frac=x, c_gas=gas(T_g, x, a_["p"])["c"], ceiling_mach=ceiling_mach(H, T_g, x),
                               lift_kg_m3=lift_per_m(H, T_g, x), heat_W_m2=w, heat_TWh_yr_spine=twh,
                               lift_swing_pct=0.46 * (T_g - a_["T"])))
    CR = pd.DataFrame(cr)
    CR.to_csv(HERE / "ceiling.csv", index=False, float_format="%.4g")
    print(CR[(CR.H_km == 20) | (CR.H_km == 25)].round(3).to_string())

    # ---------------- thickness vs Mach for concept weights (shield at aircraft altitude - 250 m)
    import hm_aircraft as ha
    rows = []
    for M in np.round(np.arange(1.7, 4.01, 0.1), 2):
        eng = "fixed" if M < 1.75 else ("adaptive" if M <= 3.05 else "tbcc")
        ac = ha.size(M, eng, kT=1.40 if eng != "tbcc" else 1.60)
        W = ac["W_mid_t"] * 1e3 * 9.81
        for Hmode in ("20 km", "preferred"):
            H = 20e3 if Hmode == "20 km" else (ac["h_cruise_km"] * 1e3 - H_AC)
            a_ = atm(H)
            for T_lab, T_g in (("ambient", None), ("+40 K", a_["T"] + 40), ("+80 K (solar noon)", a_["T"] + 80)):
                for x in (0.0, 0.02):
                    d, by = d_min(W, M, H, T_g, x)
                    lift = lift_per_m(H, T_g, x)
                    gg = gas(T_g or a_["T"], x, a_["p"])
                    rows.append(dict(M=M, W_mid_t=ac["W_mid_t"], H_mode=Hmode, H_s_km=H / 1e3, gas=T_lab, air_frac=x,
                                     q_kPa=0.7 * atm(H + H_AC)["p"] * M ** 2 / 1e3,
                                     ceiling_mach=ceiling_mach(H, T_g or a_["T"], x), d_min_m=d, set_by=by,
                                     lid_kg_m2=lid(d, H, T_g, x) if np.isfinite(d) else np.nan,
                                     h2_kg_m2=gg["rho"] * d if np.isfinite(d) else np.nan,
                                     width_km=width_km(M, d, W, H, T_g, x) if np.isfinite(d) else np.nan,
                                     T_equiv20_m=d * lift / lift_per_m(20e3) if np.isfinite(d) else np.nan))
    SM = pd.DataFrame(rows)
    SM.to_csv(HERE / "shield_vs_mach.csv", index=False, float_format="%.4g")
    print(SM[(SM.gas == "ambient") & (SM.air_frac == 0)][["M", "W_mid_t", "H_mode", "H_s_km", "q_kPa", "ceiling_mach",
                                                          "d_min_m", "set_by", "lid_kg_m2", "width_km", "T_equiv20_m"]].round(2).to_string())
