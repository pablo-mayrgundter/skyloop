"""Vacuum ("inverted balloon") BoomShield: structural sizing and mass closure.

Screening model (numpy only). For each architecture we size the members that keep an
evacuated volume open against ambient pressure at 20 km and compare structure mass per
unit *planform* area with the lift per unit planform area.

Architectures
-------------
tube  (i)/(iii)  circular tube of diameter D, axis along-track.  Skin = barrier film
                 spanning b = 0.25 m between longitudinal UHMWPE tendons; tendons sag
                 inward between rings (span s, sag/span 0.1) -> along-track tension
                 N_xx = dp s^2/(8 delta).  Rings (CFRP, lattice or monolithic tube) take
                 the radial line load q = dp s in hoop compression and must not buckle
                 (q_cr = 3 EI/R^3).  Longitudinal CFRP stringers carry the end-cap load
                 dp R/2 per unit circumference plus the tendon pull N_xx.
                 A raft (iii) of touching tubes has exactly the same mass/lift per
                 planform as one tube of the same D; a single tube (i) covering the
                 corridor is the same calculation at D = W.
lens  (ii)       flat evacuated slab of thickness T (lens rims ignored -> optimistic).
                 Faces = film + two-way tendon net spanning a strut grid of pitch a.
                 Vertical CFRP struts (length T, load dp a^2) and horizontal CFRP
                 compression grids in both face planes (resultant N + dp T/2 per unit
                 length in each plane, each direction) close the load path.

Compression members are sized as the lighter of (a) a thin-walled CFRP tube (strength,
Euler with end-fixity 0.8, local shell buckling with the NASA SP-8007 knock-down,
minimum gauge) and (b) a 3-chord lattice column (global Euler with shear knock-down,
L/500 bow imperfection amplified, chords sized as (a) between lacing nodes, +30 % lacing).
Rings use the same two options with q_cr >= SF q / 0.75 and an R/300 out-of-round
imperfection, amplified by 1/(1 - SF q / q_cr).

All loads x SF (default 1.5).  Joints/fittings +15 % on compression structure, +10 % on
tension (terminations, seams).  NOT included: tethers, pumps, power, controls, the
acoustic deck itself, end bulkheads (small when amortised over 250 m modules).

Run:  python boom-shield/vacuum/vacuum_structure.py
Outputs (boom-shield/vacuum/): structure_closure.csv, structure_closure.png
"""
from __future__ import annotations

import csv
import math
import os
import sys
from dataclasses import dataclass

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "sim"))
from boomsim.atmosphere import std_atm  # noqa: E402

ATM20 = std_atm(20000.0)
P20, RHO20, T20 = float(ATM20.p), float(ATM20.rho), float(ATM20.T)
ATM0 = std_atm(0.0)


# ------------------------------------------------------------------ materials
@dataclass(frozen=True)
class Mat:
    rho: float     # kg/m^3
    E: float       # Pa (member axis)
    sig: float     # Pa, usable stress (compression for CFRP, sustained working for fibres)


# CFRP filament-wound / pultruded tube, UD-dominated (IM fibre, ~60 % Vf):
# E ~ 120-140 GPa, UD compressive strength ~1.2-1.5 GPa.  We use 0.75 x 1.2 GPa
# (cold, UV, fatigue, defects) as the usable compressive stress; SF is applied to loads.
CFRP = Mat(rho=1550.0, E=120e9, sig=0.9e9)
E_LOCAL = 70e9           # effective modulus for local shell buckling of a wound tube
# UHMWPE (Dyneema SK75/SK78-class) fibre: ~3.5 GPa break, but creep limits sustained
# load; 0.8 GPa working (~23 % of break) is already aggressive for multi-year loading.
UHMWPE = Mat(rho=975.0, E=110e9, sig=0.8e9)
# Barrier film (PET-class laminate with metallised/AlOx barrier): 200 MPa break, 70 MPa working.
FILM = Mat(rho=1390.0, E=4e9, sig=70e6)
T_FILM_MIN = 25e-6       # m, minimum handleable/pinhole-safe gauge
T_CFRP_MIN = 0.25e-3     # m, minimum wound-tube wall

SF = 1.5
K_EULER = 0.8            # end-fixity / eccentricity knock-down on Euler loads
K_SHEAR = 0.8            # lattice shear-flexibility knock-down
F_JOINT_C = 1.15
F_JOINT_T = 1.10
LACING = 0.30
SAG = 0.10               # sag/span of tendons and film
B_FILM = 0.25            # m, tendon pitch (film span)


def gamma_sp8007(r_over_t):
    """NASA SP-8007 axial-compression knock-down for thin cylinders."""
    phi = np.sqrt(np.maximum(r_over_t, 1.0)) / 16.0
    return 1.0 - 0.901 * (1.0 - np.exp(-phi))


_R_GRID = np.logspace(np.log10(3e-3), np.log10(4.0), 260)


def tube_member(Pd, L, mat=CFRP):
    """Min mass per length of a thin CFRP tube carrying design load Pd over length L.
    Returns (mass_per_length [kg/m], area [m^2], radius [m])."""
    r = _R_GRID
    t_str = Pd / (mat.sig * 2 * np.pi * r)
    t_eul = Pd * L ** 2 / (K_EULER * np.pi ** 3 * mat.E * r ** 3)
    t = np.maximum.reduce([t_str, t_eul, np.full_like(r, T_CFRP_MIN)])
    for _ in range(8):  # local buckling with r/t-dependent knock-down
        g = gamma_sp8007(r / t)
        t_loc = np.sqrt(Pd / (2 * np.pi * 0.605 * g * E_LOCAL))
        t = np.maximum(t, t_loc)
    t = np.where(t > 0.2 * r, np.inf, t)  # not thin-walled any more: reject
    m = mat.rho * 2 * np.pi * r * t
    i = int(np.argmin(m))
    return float(m[i]), float(2 * np.pi * r[i] * t[i]), float(r[i])


_H_GRID = np.logspace(np.log10(0.05), np.log10(20.0), 70)


def lattice_member(Pd, L, mat=CFRP):
    """Min mass per length of a 3-chord lattice column (chord circumradius h)."""
    best = (np.inf, np.inf, np.nan)
    for h in _H_GRID:
        if h > L / 3:
            break
        Ac = Pd / 3 / mat.sig
        ok = False
        for _ in range(30):
            I = 1.5 * Ac * h ** 2
            Pcr = K_SHEAR * K_EULER * np.pi ** 2 * mat.E * I / L ** 2
            if Pd >= 0.95 * Pcr:
                Ac *= 1.5
                continue
            M = Pd * (L / 500.0) / (1 - Pd / Pcr)
            Pc = Pd / 3 + M / (1.5 * h)
            mc, Ac_new, _ = tube_member(Pc, h, mat)
            if abs(Ac_new - Ac) < 1e-3 * Ac:
                Ac = Ac_new
                ok = True
                break
            Ac = Ac_new
        if not ok:
            continue
        m = 3 * mat.rho * Ac * (1 + LACING)
        if m < best[0]:
            best = (m, 3 * Ac, h)
    return best


def column(Pd, L):
    """Lighter of tube and lattice. Returns (mass/len, axial area, kind)."""
    mt, At, _ = tube_member(Pd, L)
    ml, Al, _ = lattice_member(Pd, L)
    return (mt, At, "tube") if mt <= ml else (ml, Al, "lattice")


def ring_member(R, q):
    """Min-mass ring of radius R under radial line load q [N/m] (= dp * spacing).
    Returns (mass per circumferential length, q_cr/q, kind)."""
    best = (np.inf, np.nan, "")
    w0 = R / 300.0
    qd = SF * q
    # (a) monolithic tube ring
    for r in _R_GRID:
        # need EI >= qd R^3 / (3*0.75), plus strength/local with amplified moment
        t = max(T_CFRP_MIN, qd * R ** 3 / (3 * 0.75) / (CFRP.E * np.pi * r ** 3))
        for _ in range(40):
            A = 2 * np.pi * r * t
            I = np.pi * r ** 3 * t
            qcr = 3 * CFRP.E * I / R ** 3
            if qd >= 0.95 * qcr:
                t *= 1.3
                continue
            M = qd * R * w0 / (1 - qd / qcr)
            sig = qd * R / A + M * r / I
            g = gamma_sp8007(r / t)
            sig_loc = 0.605 * g * E_LOCAL * t / r
            if sig <= min(CFRP.sig, sig_loc):
                break
            t *= 1.08
        if t > 0.2 * r:
            continue
        m = CFRP.rho * 2 * np.pi * r * t
        if m < best[0]:
            best = (m, 3 * CFRP.E * np.pi * r ** 3 * t / R ** 3 / q, "tube")
    # (b) 3-chord lattice ring, depth (circumradius) h
    for h in _H_GRID:
        if h > R / 4:
            break
        Ac = qd * R / 3 / CFRP.sig
        ok = False
        for _ in range(40):
            EI = K_SHEAR * CFRP.E * 1.5 * Ac * h ** 2
            qcr = 3 * EI / R ** 3
            if qd >= 0.75 * qcr:
                Ac *= 1.3
                continue
            M = qd * R * w0 / (1 - qd / qcr)
            Pc = qd * R / 3 + M / (1.5 * h)
            _, Ac_new, _ = tube_member(Pc, h)
            if abs(Ac_new - Ac) < 1e-3 * Ac:
                Ac = Ac_new
                ok = True
                break
            Ac = Ac_new
        if not ok:
            continue
        m = 3 * CFRP.rho * Ac * (1 + LACING)
        if m < best[0]:
            EI = K_SHEAR * CFRP.E * 1.5 * Ac * h ** 2
            best = (m, 3 * EI / R ** 3 / q, "lattice")
    return best


def skin(dp, span, b=B_FILM):
    """Film + tendon mass per unit skin area, tendon tension N [N/m] for tendons spanning
    `span` with sag/span SAG, film spanning b."""
    N_f = dp * b / (8 * SAG)
    t_f = max(T_FILM_MIN, SF * N_f / FILM.sig)
    m_film = FILM.rho * t_f * 1.10
    N = dp * span / (8 * SAG)
    m_tend = UHMWPE.rho * SF * N / UHMWPE.sig * (1 + 8 / 3 * SAG ** 2) * F_JOINT_T
    return m_film, m_tend, N


# ------------------------------------------------------------------ architectures
def tube_design(D, dp=P20, rho_air=RHO20, p_ratio=1e-3, s_grid=None):
    R = D / 2
    s_grid = s_grid or [0.75, 1.0, 1.5, 2.0, 3.0, 4.0, 6.0, 8.0, 12.0]
    lift = (1 - p_ratio) * rho_air * np.pi * D / 4
    best = None
    for s in s_grid:
        m_film, m_tend, Nxx = skin(dp, s)
        mr, F_ring, rkind = ring_member(R, dp * s)
        ring_pp = mr * 2 * np.pi * R / s / D * F_JOINT_C
        n_x = dp * R / 2 + Nxx
        st_best = (np.inf, None, None)
        for ns in (3, 4, 6, 8, 12, 24, 48, 96, 192, 384):
            P = SF * n_x * 2 * np.pi * R / ns
            mu, A, kind = column(P, s)
            tot = mu * ns / D * F_JOINT_C
            if tot < st_best[0]:
                st_best = (tot, ns, kind)
        rec = dict(arch="tube", D=D, T_or_D=D, s_or_a=s, lift=lift,
                   film=m_film * np.pi, tendon=m_tend * np.pi, ring=ring_pp,
                   stringer=st_best[0], strut_v=0.0, strut_h=0.0,
                   N_xx=Nxx, F_ring=F_ring, ring_kind=rkind, n_stringers=st_best[1],
                   stringer_kind=st_best[2])
        rec["total"] = rec["film"] + rec["tendon"] + rec["ring"] + rec["stringer"]
        if best is None or rec["total"] < best["total"]:
            best = rec
    best["ratio"] = best["total"] / best["lift"]
    R = D / 2
    # "perfect structure" floor: hoop + axial(+tendon pull) compression at the usable
    # stress with SF and joints, no buckling penalty (ideal hierarchical members)
    s_ = best["s_or_a"]
    comp = CFRP.rho * SF * F_JOINT_C / CFRP.sig * (dp * R * 2 * np.pi * R + (dp * R / 2 + best["N_xx"]) * 2 * np.pi * R) / D
    best["floor"] = best["film"] + best["tendon"] + comp
    best["floor_ratio"] = best["floor"] / best["lift"]
    # ring coupling stiffness per unit planform area, bracket (see review sec. C)
    F = best["F_ring"]
    best["K_ring_lo"] = 1.1 * (F - 1) * dp / R     # diametral squeeze (point loads)
    best["K_ring_hi"] = 3.0 * (F - 1) * dp / R     # n = 2 modal stiffness incl. pressure softening
    # membrane wave speed along-track: tendon tension over smeared skin-area mass
    best["c_m"] = math.sqrt(best["N_xx"] * np.pi / best["total"])
    best["c_m_skin_only"] = math.sqrt(best["N_xx"] / ((best["film"] + best["tendon"]) / np.pi))
    return best


def lens_design(T, dp=P20, rho_air=RHO20, p_ratio=1e-3, a_grid=None):
    a_grid = a_grid or [1.0, 1.5, 2.0, 3.0, 4.0, 6.0, 8.0, 12.0, 16.0, 24.0]
    lift = (1 - p_ratio) * rho_air * T
    best = None
    for a in a_grid:
        m_film, m_tend, N = skin(dp, a)
        faces = 2 * (m_film + m_tend)
        mv, Av, kv = column(SF * dp * a * a, T)
        strut_v = mv * T / a ** 2 * F_JOINT_C
        Ph = SF * (N + dp * T / 2) * a
        mh, Ah, kh = column(Ph, a)
        strut_h = 2 * 2 * mh / a * F_JOINT_C
        rec = dict(arch="lens", D=np.nan, T_or_D=T, s_or_a=a, lift=lift,
                   film=2 * m_film, tendon=2 * m_tend, ring=0.0, stringer=0.0,
                   strut_v=strut_v, strut_h=strut_h, N_xx=N, F_ring=np.nan,
                   ring_kind="", n_stringers=0, stringer_kind=f"{kv}/{kh}",
                   A_strut=Av)
        rec["total"] = faces + strut_v + strut_h
        if best is None or rec["total"] < best["total"]:
            best = rec
    best["ratio"] = best["total"] / best["lift"]
    a = best["s_or_a"]
    comp = CFRP.rho * SF * F_JOINT_C / CFRP.sig * (dp * T + 2 * 2 * (best["N_xx"] + dp * T / 2))
    best["floor"] = best["film"] + best["tendon"] + comp
    best["floor_ratio"] = best["floor"] / best["lift"]
    best["K_strut"] = CFRP.E * best["A_strut"] / (T * a ** 2)
    best["c_m"] = math.sqrt(best["N_xx"] * 2 / best["total"])
    best["c_m_skin_only"] = math.sqrt(best["N_xx"] / ((best["film"] + best["tendon"]) / 2))
    return best


def gas_reference(D_or_T, frac, shape="tube"):
    """He/H2 envelope of the same shape: lift fraction frac, film envelope ~60 g/m^2
    (incl. seams) times skin area per planform; slightly superpressured (no compression)."""
    if shape == "tube":
        lift = frac * RHO20 * np.pi * D_or_T / 4
        env = 0.060 * np.pi
    else:
        lift = frac * RHO20 * D_or_T
        env = 0.060 * 2
    return lift, env


# ------------------------------------------------------------------ main
COLS = ["case", "arch", "T_or_D", "s_or_a", "lift", "film", "tendon", "ring", "stringer",
        "strut_v", "strut_h", "total", "ratio", "floor", "floor_ratio", "net", "N_xx", "F_ring", "ring_kind",
        "n_stringers", "stringer_kind", "K_ring_lo", "K_ring_hi", "K_strut", "c_m",
        "c_m_skin_only", "he_net", "h2_net"]


def run_all():
    rows = []
    for D in (10, 15, 20, 30, 45, 60, 100, 250):
        r = tube_design(D)
        r["case"] = f"tube D={D} m (20 km)"
        rows.append(r)
    for T in (10, 15, 30, 60, 100):
        r = lens_design(T)
        r["case"] = f"lens T={T} m (20 km)"
        rows.append(r)
    # altitude scaling check: same tube at sea level (dp = 101 kPa, rho = 1.225)
    for D in (30, 100):
        r = tube_design(D, dp=float(ATM0.p), rho_air=float(ATM0.rho))
        r["case"] = f"tube D={D} m (sea level)"
        rows.append(r)
    for r in rows:
        r["net"] = r["lift"] - r["total"]
        shape = "tube" if r["arch"] == "tube" else "lens"
        if "sea" in r["case"]:
            r["he_net"] = r["h2_net"] = np.nan
            continue
        for key, frac in (("he_net", 0.86), ("h2_net", 0.93)):
            lift, env = gas_reference(r["T_or_D"], frac, shape)
            r[key] = lift - env
    return rows


def strength_floor(dp=P20, rho_air=RHO20):
    """Pure-strength lower bound: hoop (dp R) + axial (dp R/2) compression at the usable
    CFRP stress, no buckling, no SF: M/M_lift = 3 dp rho_c / (sigma rho_air) for a cylinder;
    sphere bound 1.5 dp rho_c/(sigma rho_air)."""
    k = dp * CFRP.rho / (CFRP.sig * rho_air)
    return dict(cyl=3 * k, cyl_SF=3 * k * SF, sphere=1.5 * k)


def monocoque_sphere(dp=P20, rho_air=RHO20, E=70e9, rho_s=1550.0, knock=0.25, sf=1.5):
    """Homogeneous spherical shell: buckling p_cr = 2E(t/R)^2/sqrt(3(1-nu^2)) * knock >= sf dp.
    Returns shell mass / displaced-air mass."""
    tR = math.sqrt(sf * dp * math.sqrt(3 * (1 - 0.3 ** 2)) / (2 * E * knock))
    return 3 * rho_s * tR / rho_air


def plot(rows, path):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    ink, grid, surf = "#0b0b0b", "#e4e3df", "#fcfcfb"
    comps = [("film", "film", "#9ec5f0"), ("tendon", "UHMWPE tendons", "#2a78d6"),
             ("ring", "CFRP rings", "#eb6834"), ("stringer", "CFRP stringers", "#f2a35e"),
             ("strut_v", "vertical struts", "#7a4fd0"), ("strut_h", "horizontal compression grid", "#b39ae6")]
    sel = [r for r in rows if "sea" not in r["case"]]
    fig, ax = plt.subplots(figsize=(11, 5))
    fig.patch.set_facecolor(surf)
    ax.set_facecolor(surf)
    x = np.arange(len(sel))
    bottom = np.zeros(len(sel))
    for k, lab, c in comps:
        v = np.array([r[k] / r["lift"] for r in sel])
        ax.bar(x, v, bottom=bottom, color=c, label=lab, width=0.7)
        bottom += v
    ax.axhline(1.0, color=ink, lw=1.6, label="vacuum gross lift (1 - p_ratio) = floats")
    ax.scatter(x, [r["floor_ratio"] for r in sel], marker="_", s=700, color="#c0392b", lw=2.2,
               zorder=5, label="ideal floor (no buckling penalty)")
    ax.scatter(x, [r["he_net"] / r["lift"] for r in sel], marker="_", s=700, color="#1baf7a", lw=2.2,
               zorder=5, label="He envelope net lift, same shape (fraction of vacuum gross)")
    ax.set_xticks(x)
    ax.set_xticklabels([r["case"].replace(" (20 km)", "") for r in sel], rotation=30, ha="right")
    ax.set_ylabel("structure mass / vacuum gross lift")
    ax.set_ylim(0, 1.6)
    ax.grid(True, axis="y", color=grid, lw=0.6)
    for s_ in ("top", "right"):
        ax.spines[s_].set_visible(False)
    ax.legend(frameon=False, fontsize=8, ncol=3, loc="upper center", bbox_to_anchor=(0.5, -0.28))
    ax.set_title("Vacuum shield at 20 km: structure mass as a fraction of lift "
                 "(SF 1.5, CFRP compression, UHMWPE tension, film skin)", fontsize=10, color=ink)
    fig.tight_layout()
    fig.savefig(path, dpi=130, facecolor=surf)
    plt.close(fig)


def main():
    rows = run_all()
    out = os.path.join(HERE, "structure_closure.csv")
    with open(out, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=COLS, extrasaction="ignore")
        w.writeheader()
        for r in rows:
            w.writerow({k: (f"{r.get(k):.4g}" if isinstance(r.get(k), float) else r.get(k, ""))
                        for k in COLS})
    plot(rows, os.path.join(HERE, "structure_closure.png"))
    print(f"20 km: p = {P20:.0f} Pa, rho = {RHO20:.4f}, T = {T20:.2f} K")
    sf = strength_floor()
    sl = strength_floor(float(ATM0.p), float(ATM0.rho))
    print(f"strength floor (cylinder, no SF/buckling) M/lift = {sf['cyl']:.2f} at 20 km, "
          f"{sl['cyl']:.2f} at sea level; with SF {sf['cyl_SF']:.2f}; sphere bound {sf['sphere']:.2f}")
    print(f"monocoque CFRP sphere M/lift: 20 km {monocoque_sphere():.1f}, sea level "
          f"{monocoque_sphere(float(ATM0.p), float(ATM0.rho)):.1f}")
    hdr = f"{'case':28s} {'s/a':>5s} {'lift':>6s} {'film':>6s} {'tend':>6s} {'ring':>6s} {'strg':>6s} {'strV':>6s} {'strH':>6s} {'total':>7s} {'ratio':>6s} {'F':>5s} {'Kring':>11s} {'Kstrut':>8s} {'c_m':>5s} {'cm_sk':>6s}"
    print(hdr)
    for r in rows:
        kr = f"{r.get('K_ring_lo', np.nan):.0f}-{r.get('K_ring_hi', np.nan):.0f}" if r["arch"] == "tube" else "-"
        print(f"{r['case']:28s} {r['s_or_a']:5.2f} {r['lift']:6.2f} {r['film']:6.3f} {r['tendon']:6.3f} "
              f"{r['ring']:6.2f} {r['stringer']:6.2f} {r['strut_v']:6.2f} {r['strut_h']:6.2f} "
              f"{r['total']:7.2f} {r['ratio']:6.2f} {r['F_ring']:5.2f} {kr:>11s} "
              f"{r.get('K_strut', np.nan):8.3g} {r['c_m']:5.0f} {r['c_m_skin_only']:6.0f} floor {r['floor_ratio']:.2f} "
              f"[{r['ring_kind']}/{r['stringer_kind']}, n_s={r['n_stringers']}]")


if __name__ == "__main__":
    main()
