"""BoomShield operating-altitude trade (non-acoustic side).

Screening arithmetic only (numpy + the US Std Atmosphere 1976 from boomsim).
Run:  python3 boom-shield/altitude/altitude_trade.py
Writes into boom-shield/altitude/:
  altitude_metrics.csv   one row per shield altitude H_s
  lid_thickness.csv      H2/He thickness for lid masses 1.2-4.8 kg/m^2
  wind_profiles.csv      the assumed climatological wind profiles
  altitude_trade.png     four-panel figure

Every wind number is an ASSUMED representative mid-latitude (30-45 N) profile
built from the summaries cited in altitude-trade.md; it is not an ERA5 product.
"""
from __future__ import annotations

import csv
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "sim"))
from boomsim.atmosphere import std_atm, R_STAR, G0  # noqa: E402  (read-only import)

# ---------------------------------------------------------------- constants
M_AIR, M_H2, M_HE = 28.9644e-3, 2.016e-3, 4.0026e-3
W_SHIELD = 1000.0          # m, shield width
T_ACOUSTIC = 30.0          # m, minimum layer thickness for the acoustic mechanism
M_SYS = 2.46               # kg/m^2, baseline system mass (lift of 30 m H2 at 20 km)
H_GAP = 250.0              # m, aircraft above shield
ALTS = np.array([12, 14, 16, 18, 20, 22, 25], dtype=float)   # km
LID_MASSES = [1.2, 2.4, 3.6, 4.8]                            # kg/m^2

# ---------------------------------------------------------------- winds (ASSUMED)
Z_W = np.array([0, 2, 5, 8, 10.5, 12, 14, 16, 18, 20, 22, 25], float)  # km
U_WINTER = np.array([8, 15, 25, 40, 50, 48, 35, 25, 20, 18, 20, 25], float)
U_SUMMER = np.array([5, 8, 12, 20, 25, 24, 16, 10, 6, 5, 7, 12], float)
U_DESIGN = np.array([20, 30, 50, 75, 90, 85, 65, 50, 42, 40, 42, 50], float)  # ~P99 winter


def U(profile, z_km):
    return np.interp(z_km, Z_W, profile)


# ---------------------------------------------------------------- helpers
def atm(z_km):
    return std_atm(np.asarray(z_km, float) * 1e3)


def lift_density(z_km, M_gas, superheat=0.0):
    """Net lift per m^3 of gas [kg/m^3] at ambient p,T (gas at T+superheat)."""
    a = atm(z_km)
    rho_g = a.p * M_gas / (R_STAR * (a.T + superheat))
    return a.rho - rho_g, rho_g


def skin_cf(rho, mu, Uw, L=W_SHIELD):
    Re = rho * Uw * L / mu
    return 0.455 / np.log10(Re) ** 2.58


def mu_air(T):  # Sutherland
    return 1.458e-6 * T ** 1.5 / (T + 110.4)


def q_flight(z_km, M):
    return 0.5 * 1.4 * atm(z_km).p * M ** 2


# ---------------------------------------------------------------- tethers
SIGMA_F = 1.2e9      # Pa, UHMWPE working stress on fibre area (~3.5 GPa fibre / 3)
RHO_F = 975.0        # kg/m^3
FILL = 0.65          # fibre fill of rope cross-section
CD_TETHER = 1.2
V_OVER_H = 1.5       # vertical/horizontal tension ratio at the top (34 deg from vertical)


def tether(H_s_km, F_side_per_m, spacing=3000.0, prof=U_DESIGN):
    z = np.linspace(0.0, H_s_km, 400)
    dz = (z[1] - z[0]) * 1e3
    a = atm(z)
    qz = 0.5 * a.rho * U(prof, z) ** 2
    d = 0.08
    for _ in range(60):
        D = np.sum(qz * CD_TETHER * d) * dz           # drag, tether taken vertical (conservative)
        H_top = F_side_per_m * spacing + 0.5 * D
        L = H_s_km * 1e3 / np.cos(np.arctan(1 / V_OVER_H))
        A_f = np.sqrt(H_top ** 2 + (V_OVER_H * H_top) ** 2) / SIGMA_F
        m_t = RHO_F * A_f * L
        V_top = V_OVER_H * H_top + m_t * G0 * 0.5
        T_top = np.hypot(H_top, V_top)
        A_f = T_top / SIGMA_F
        d_new = np.sqrt(4 * A_f / (np.pi * FILL))
        if abs(d_new - d) < 1e-6:
            break
        d = d_new
    m_t = RHO_F * A_f * L
    lift_cost = (V_top + 0.5 * m_t * G0) / G0 / (spacing * W_SHIELD)   # kg/m^2 of shield
    # share of drag from the 8-13 km jet band
    jet = (z >= 8) & (z <= 13)
    jet_frac = np.sum(qz[jet]) / np.sum(qz)
    return dict(D=D, H_top=H_top, T_top=T_top, d=d, m_t=m_t, L=L,
                lift_cost=lift_cost, jet_frac=jet_frac)


# ---------------------------------------------------------------- fall time
def fall_time(z_top_km, z_bot_km, m=M_SYS, cd=1.2):
    """Time [min] for a deflated membrane of areal mass m to fall between levels."""
    if z_top_km <= z_bot_km:
        return 0.0
    z = np.linspace(z_bot_km, z_top_km, 400)
    vt = np.sqrt(2 * m * G0 / (atm(z).rho * cd))
    return np.trapezoid(1.0 / vt, z * 1e3) / 60.0


# ---------------------------------------------------------------- superheat (toy)
H_R = 0.30           # W/m^2K, linearised radiative exchange of a low-emissivity film
H_C20 = 0.60         # W/m^2K, film-air convection at 20 km (natural/forced, ~rho^(2/3))
DT20 = 20.0          # K, baseline design superheat at 20 km (structure-concept.md)


def superheat(z_km):
    rho = atm(z_km).rho
    rho20 = atm(20.0).rho
    hc = H_C20 * (rho / rho20) ** (2 / 3)
    q = DT20 * (H_C20 + H_R)
    return q / (hc + H_R)


# ---------------------------------------------------------------- carlson scaling
def boom_ratio(Ha_km, Ha_ref=20.25, p_exp=0.5):
    """Ground-overpressure ratio vs reference aircraft altitude.

    Carlson (NASA TP-1122): dp ~ K_s * sqrt(p_v p_g) * h_e^(-3/4).  Volume-dominated
    (K_s const) -> p_v^0.5; lift-dominated (K_s ~ sqrt(K_L), K_L ~ W/p_v) -> p_v^0.
    h_e approximated by geometric altitude.
    """
    pv, pr = atm(Ha_km).p, atm(Ha_ref).p
    return (pv / pr) ** p_exp * (Ha_km / Ha_ref) ** -0.75


def main():
    rows = []
    for Hs in ALTS:
        a = atm(Hs)
        dl_h2, rho_h2 = lift_density(Hs, M_H2)
        dl_he, _ = lift_density(Hs, M_HE)
        T_lift = M_SYS / dl_h2
        T_use = max(T_ACOUSTIC, T_lift)
        lift_use = dl_h2 * T_use
        surplus = lift_use - M_SYS
        h2_kg_m2 = rho_h2 * T_use
        inv_km = h2_kg_m2 * W_SHIELD * 1e3 / 1e3          # t per km
        inv_corr = inv_km * 4000 / 1e6                    # Mt for 4000 km
        dT = superheat(Hs)
        swing = dT / float(a.T)                           # fractional lift swing
        mu = mu_air(float(a.T))
        out = dict(H_s_km=Hs, T_K=float(a.T), p_Pa=float(a.p), rho=float(a.rho),
                   lift_H2_kg_m3=float(dl_h2), lift_He_kg_m3=float(dl_he),
                   T_for_246_H2_m=float(T_lift), T_used_m=float(T_use),
                   lift_at_T_used=float(lift_use), surplus_lift_kg_m2=float(surplus),
                   H2_t_per_km=float(inv_km), H2_Mt_4000km=float(inv_corr),
                   superheat_K=float(dT), lift_swing_pct=100 * float(swing),
                   lift_swing_kg_m2=float(lift_use * swing),
                   breathing_m=float(T_use * swing))
        for tag, prof in (("mean_winter", U_WINTER), ("mean_summer", U_SUMMER), ("design", U_DESIGN)):
            Uw = U(prof, Hs)
            qd = 0.5 * float(a.rho) * Uw ** 2
            cf = skin_cf(float(a.rho), mu, Uw)
            F = qd * (1.0 * T_use + 2 * W_SHIELD * cf)
            out[f"U_{tag}"] = float(Uw)
            out[f"q_{tag}_Pa"] = qd
            out[f"Fside_{tag}_kN_per_m"] = F / 1e3
            # gust lift (w' = 5 % U, thin-aerofoil 2*pi*alpha) and 1 deg attitude lift
            out[f"gustlift_{tag}_Pa"] = qd * 2 * np.pi * 0.05
        out["lid_weight_Pa"] = M_SYS * G0
        tet = tether(Hs, out["Fside_design_kN_per_m"] * 1e3)
        out.update({f"tether_{k}": float(v) for k, v in tet.items()})
        tet_mean = tether(Hs, out["Fside_mean_winter_kN_per_m"] * 1e3, prof=U_WINTER)
        out["tether_winter_D_kN"] = tet_mean["D"] / 1e3
        out["fall_min_to_FL510"] = fall_time(Hs, 15.54)
        out["fall_min_to_FL410"] = fall_time(Hs, 12.50)
        Ha = Hs + H_GAP / 1e3
        out["H_aircraft_km"] = Ha
        for M in (1.7, 2.0, 3.2):
            out[f"q_M{M}_kPa"] = float(q_flight(Ha, M)) / 1e3
        out["WS_opt_M2_Pa"] = float(q_flight(Ha, 2.0)) * 0.14
        out["near_field_dp_scale"] = float(atm(Ha).p / atm(20.25).p)
        out["boom_vol_dB"] = 20 * np.log10(boom_ratio(Ha, p_exp=0.5))
        out["boom_lift_dB"] = 20 * np.log10(boom_ratio(Ha, p_exp=0.0))
        rows.append(out)

    keys = list(rows[0].keys())
    with open(os.path.join(HERE, "altitude_metrics.csv"), "w", newline="") as f:
        w = csv.DictWriter(f, keys)
        w.writeheader()
        for r in rows:
            w.writerow({k: (f"{v:.4g}" if isinstance(v, float) else v) for k, v in r.items()})

    with open(os.path.join(HERE, "lid_thickness.csv"), "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["H_s_km"] + [f"H2_T_m_lid{m}" for m in LID_MASSES] + [f"He_T_m_lid{m}" for m in LID_MASSES])
        for Hs in ALTS:
            h2 = lift_density(Hs, M_H2)[0]
            he = lift_density(Hs, M_HE)[0]
            # lid + 0.1 floor + 0.28 structure/controls (structure-concept.md table 2)
            w.writerow([Hs] + [f"{(m + 0.38) / h2:.1f}" for m in LID_MASSES]
                       + [f"{(m + 0.38) / he:.1f}" for m in LID_MASSES])

    with open(os.path.join(HERE, "wind_profiles.csv"), "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["z_km", "U_mean_winter", "U_mean_summer", "U_design_P99"])
        for i, z in enumerate(Z_W):
            w.writerow([z, U_WINTER[i], U_SUMMER[i], U_DESIGN[i]])

    make_figure(rows)
    for r in rows:
        print(f"{r['H_s_km']:4.0f} km  T_lift {r['T_for_246_H2_m']:5.1f} m  T_used {r['T_used_m']:4.0f}  "
              f"surplus {r['surplus_lift_kg_m2']:5.2f}  H2 {r['H2_Mt_4000km']:.2f} Mt  dT {r['superheat_K']:4.1f} K "
              f"swing {r['lift_swing_pct']:4.1f}%  Fside_d {r['Fside_design_kN_per_m']:.2f} kN/m  "
              f"gust_d {r['gustlift_design_Pa']:.1f} Pa  tether D {r['tether_D']/1e3:.0f} kN  "
              f"Ttop {r['tether_T_top']/1e6:.2f} MN d {100*r['tether_d']:.1f} cm m {r['tether_m_t']/1e3:.0f} t "
              f"liftcost {r['tether_lift_cost']:.2f}  jetfrac {r['tether_jet_frac']:.2f}  "
              f"fall510 {r['fall_min_to_FL510']:.1f} min  q2 {r['q_M2.0_kPa']:.1f}  WS {r['WS_opt_M2_Pa']:.0f}  "
              f"boom +{r['boom_lift_dB']:.1f}..{r['boom_vol_dB']:.1f} dB  nf x{r['near_field_dp_scale']:.2f}")


def make_figure(rows):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    C = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100"]
    ink, ink2, grid = "#0b0b0b", "#52514e", "#e4e3df"
    plt.rcParams.update({"font.size": 9, "axes.edgecolor": ink2, "axes.labelcolor": ink,
                         "xtick.color": ink2, "ytick.color": ink2, "axes.grid": True,
                         "grid.color": grid, "grid.linewidth": 0.6,
                         "axes.spines.top": False, "axes.spines.right": False,
                         "figure.facecolor": "#fcfcfb", "axes.facecolor": "#fcfcfb"})
    fig, ax = plt.subplots(2, 2, figsize=(10, 7.6))
    H = np.array([r["H_s_km"] for r in rows])
    Hf = np.linspace(12, 25, 131)

    # (a) thickness for lid masses
    a = ax[0, 0]
    for i, m in enumerate([1.2, 2.4, 4.8]):
        h2 = np.array([(m + 0.38) / lift_density(z, M_H2)[0] for z in Hf])
        a.plot(Hf, h2, color=C[i], lw=2, label=f"H$_2$, lid {m} kg/m²")
    he = np.array([(2.4 + 0.38) / lift_density(z, M_HE)[0] for z in Hf])
    a.plot(Hf, he, color=C[1], lw=1.2, ls="--", label="He, lid 2.4 kg/m²")
    a.axhline(30, color=ink2, lw=1, ls=":")
    a.text(12.2, 32, "acoustic minimum ≈ 30 m", color=ink2, fontsize=8)
    a.set_yscale("log")
    a.set_xlabel("shield altitude H$_s$ [km]")
    a.set_ylabel("layer thickness to float the system [m]")
    a.set_title("(a) Thickness needed (lid + 0.38 kg/m² floor/structure)", loc="left", fontsize=9.5)
    a.legend(frameon=False, fontsize=8)

    # (b) wind profiles + rho U^2
    b = ax[0, 1]
    zf = np.linspace(0, 25, 251)
    for i, (lab, prof) in enumerate((("winter mean", U_WINTER), ("summer mean", U_SUMMER),
                                     ("design (~P99)", U_DESIGN))):
        b.plot(U(prof, zf), zf, color=C[i], lw=2, label=lab)
    b.axhspan(18.29, 25, color=C[0], alpha=0.06)
    b.text(70, 21.5, "above FL600\n(upper Class E)", color=ink2, fontsize=8)
    b.axhspan(8.8, 15.5, color=C[3], alpha=0.10)
    b.text(57, 15.7, "airliner/bizjet cruise FL290–FL510", color=ink2, fontsize=8, va="bottom")
    b.set_xlabel("wind speed [m/s] (assumed mid-latitude profiles)")
    b.set_ylabel("altitude [km]")
    b.set_title("(b) Winds the shield and its tethers see", loc="left", fontsize=9.5)
    b.legend(frameon=False, fontsize=8, loc="lower right")

    # (c) tether drag & mass, shield side load (design wind, 3 km spacing)
    c = ax[1, 0]
    D = np.array([r["tether_D"] for r in rows]) / 1e6
    Fs = np.array([r["Fside_design_kN_per_m"] for r in rows]) * 3000 / 1e3
    Tt = np.array([r["tether_T_top"] for r in rows]) / 1e6
    c.plot(H, Tt, "-o", color=C[0], lw=2, ms=5, label="tether top tension")
    c.plot(H, Fs, "-o", color=C[1], lw=2, ms=5, label="shield side load per 3 km")
    c.plot(H, D, "-o", color=C[2], lw=2, ms=5, label="wind drag on the tether itself")
    c.set_yscale("log")
    c.set_xlabel("shield altitude H$_s$ [km]")
    c.set_ylabel("force [MN]")
    c.set_title("(c) Per tether, design wind, 3 km spacing", loc="left", fontsize=9.5)
    c.legend(frameon=False, fontsize=8)

    # (d) aircraft side
    d = ax[1, 1]
    for i, M in enumerate((1.7, 2.0, 3.2)):
        d.plot(Hf, [q_flight(z + 0.25, M) / 1e3 for z in Hf], color=C[i], lw=2, label=f"q at Mach {M}")
    pts = [("Concorde M2.02\n16.8–18 km", 17.2, q_flight(17.2, 2.02) / 1e3),
           ("Overture M1.7\nFL600", 18.29, q_flight(18.29, 1.7) / 1e3),
           ("SR-71 M3.2\n~24 km", 24.0, q_flight(24.0, 3.2) / 1e3)]
    for lab, z, q in pts:
        d.plot([z - 0.25], [q], "o", color=ink, ms=6)
        d.annotate(lab, (z - 0.25, q), xytext=(6, 4), textcoords="offset points", fontsize=7.5, color=ink2)
    d.set_yscale("log")
    d.set_xlabel("shield altitude H$_s$ [km] (aircraft at H$_s$ + 250 m)")
    d.set_ylabel("flight dynamic pressure q [kPa]")
    d.set_title("(d) What H$_s$ asks of the aircraft", loc="left", fontsize=9.5)
    d.legend(frameon=False, fontsize=8, loc="upper right")

    fig.tight_layout()
    fig.savefig(os.path.join(HERE, "altitude_trade.png"), dpi=150)


if __name__ == "__main__":
    main()
