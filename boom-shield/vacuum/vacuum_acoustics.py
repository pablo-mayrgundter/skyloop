"""Vacuum ("inverted balloon") BoomShield: acoustics vs the Gate-0 He/H2 stacks.

Parts
-----
1. Gas-state checks for low-pressure air at 20 km: sound speed, impedance, mean free path,
   classical absorption, thermal boundary layer.
2. Shield-plane sweep (plane wave at trace speed V = M c20, boomsim.shield):
   p_ratio in {1e-1, 1e-2, 1e-3} x d in {15, 30, 60} m x M in {1.4, 2.0, 2.5, 3.2},
   lift-closed faces (all slab lift (1-p_ratio) rho d on the two faces, as Gate 0),
   ideal 0.1 s and 0.3 s N-waves -> peak and positive-impulse ratios; He/H2 references.
3. Structural coupling (smeared tie_K) and face tension (Membrane N): sweeps.
4. Ground PLdB (optional, needs boomsim.propagate): same pipeline and rise-time proxy as
   run_gate0.py (logic copied, not imported), h = 250 m, four aircraft cases.

Run:  python boom-shield/vacuum/vacuum_acoustics.py [--ground] [--procs 4]
Outputs (boom-shield/vacuum/): acoustics_shield_plane.csv, acoustics_ties_tension.csv,
acoustics_ground.csv (with --ground), acoustics_transfer.png, acoustics_coincidence.png
"""
from __future__ import annotations

import argparse
import csv
import dataclasses
import math
import os
import sys
import time
from multiprocessing import Pool

import numpy as np
from scipy.ndimage import gaussian_filter1d

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "sim"))
from boomsim import shield as S  # noqa: E402
from boomsim.atmosphere import std_atm  # noqa: E402

H_S = 20000.0
ATM20 = std_atm(H_S)
T20, P20, RHO20, C20 = float(ATM20.T), float(ATM20.p), float(ATM20.rho), float(ATM20.c)
MACHS = [1.4, 2.0, 2.5, 3.2]
P_RATIOS = [1e-1, 1e-2, 1e-3]
DS = [15.0, 30.0, 60.0]
FS = 20e3

# structural coupling estimates from vacuum_structure.py (D = 30 m raft / T = 30 m lens)
K_RING_30 = (575.0, 1569.0)   # Pa/m, ring ovalisation bracket, D = 30 m tube
K_STRUT_30 = 7.66e4           # Pa/m, CFRP struts, T = 30 m lens


# --------------------------------------------------------------------- gas checks
def gas_checks(p_ratio):
    kB = 1.380649e-23
    d_mol = 3.7e-10                       # m, kinetic diameter of air molecules
    p = p_ratio * P20
    mfp = kB * T20 / (math.sqrt(2) * math.pi * d_mol ** 2 * p)
    g = S.gas("air", T20, p)
    mu, kappa, cp, gam = 1.42e-5, 0.0195, 1005.0, 1.4   # air at 217 K
    # classical (viscous + thermal) absorption, Np/m: alpha = w^2/(2 rho c^3) [4/3 mu + (g-1) k/cp]
    def alpha(f):
        w = 2 * math.pi * f
        return w ** 2 / (2 * g.rho * g.c ** 3) * (4 / 3 * mu + (gam - 1) * kappa / cp)
    dth = lambda f: math.sqrt(2 * kappa / (g.rho * cp * 2 * math.pi * f))
    return dict(p_ratio=p_ratio, p_Pa=p, rho=g.rho, c=g.c, Z=g.rho * g.c, mfp_m=mfp,
                alpha_10Hz_dB_per_100m=8.686 * alpha(10) * 100,
                alpha_100Hz_dB_per_100m=8.686 * alpha(100) * 100,
                alpha_1kHz_dB_per_100m=8.686 * alpha(1000) * 100,
                thermal_BL_10Hz_m=dth(10))


# --------------------------------------------------------------------- stacks
def vac_stack(p_ratio, d, tie_K=0.0, N=0.0, eta=0.0, face_mass=None):
    g = S.gas("air", T20, p_ratio * P20)
    b = (RHO20 - g.rho) * d if face_mass is None else 2 * face_mass
    return S.Stack([S.Membrane(b / 2, N=N, eta=eta), S.GasLayer(g, d, tie_K=tie_K),
                    S.Membrane(b / 2, N=N, eta=eta)]), g


def gas_stack(fill, d, tie_K=0.0):
    g = S.gas(fill, T20, P20)
    b = (RHO20 - g.rho) * d
    return S.Stack([S.Membrane(b / 2), S.GasLayer(g, d, tie_K=tie_K), S.Membrane(b / 2)]), g


def f0_mass_air_mass(p_ratio, d, V, m_face):
    """Double-wall resonance with the low-pressure gas spring (propagating layer):
    K_g = rho_g / (d (1/c^2 - 1/V^2)) = gamma p_int / (d cos^2 theta)."""
    g = S.gas("air", T20, p_ratio * P20)
    Kg = g.rho / (d * (1 / g.c ** 2 - 1 / V ** 2))
    return Kg, math.sqrt(Kg * 2 / m_face) / (2 * math.pi)


# --------------------------------------------------------------------- N-waves
def nwave(dur, fs=FS, pad=8.0):
    n = int(round(dur * fs))
    core = np.linspace(1.0, -1.0, n)
    npad = int(round(pad * fs))
    p = np.concatenate([np.zeros(npad), core, np.zeros(npad)])
    t = (np.arange(p.size) - npad) / fs
    return t, p


def pos_impulse(t, p):
    """Positive impulse of the first positive phase (copied from run_gate0._pos_impulse)."""
    pk = float(np.max(p))
    if pk <= 0:
        return float("nan")
    ih = int(np.argmax(p >= 0.5 * pk))
    nonpos = np.nonzero(p[:ih] <= 0)[0]
    i0 = nonpos[-1] if nonpos.size else 0
    after = np.nonzero(p[ih:] <= 0)[0]
    i1 = ih + after[0] if after.size else p.size - 1
    seg = np.clip(p[i0:i1 + 1], 0.0, None)
    return float(np.sum(0.5 * (seg[1:] + seg[:-1])) * (t[1] - t[0]))


_NW = {d: nwave(d) for d in (0.1, 0.3)}
_IMP = {d: pos_impulse(*_NW[d]) for d in _NW}


def shield_ratios(stack, V):
    out = {}
    for dur, (t, p) in _NW.items():
        y = S.apply(stack, t, p, V, ATM20, pad_factor=2)
        out[f"peak_{dur:g}"] = float(np.max(y))
        out[f"imp_{dur:g}"] = pos_impulse(t, y) / _IMP[dur]
        out[f"rms_{dur:g}"] = float(np.sqrt(np.sum(y ** 2) / np.sum(p ** 2)))
    return out


def band_T(stack, V, fs=(2, 5, 10, 20, 50, 100)):
    T, _ = S.transfer(stack, np.array(fs, float), V, ATM20)
    return {f"T_{f}Hz_dB": 20 * math.log10(max(abs(x), 1e-30)) for f, x in zip(fs, T)}


# --------------------------------------------------------------------- part 2/3
def shield_plane_sweep():
    rows = []
    for M in MACHS:
        V = M * C20
        for fill in ("helium", "hydrogen"):
            for d in DS:
                st, g = gas_stack(fill, d)
                r = dict(stack=f"{'He' if fill == 'helium' else 'H2'} {d:g} m", fill=fill,
                         p_ratio=1.0, d=d, M=M, V=V, face_mass=(RHO20 - g.rho) * d / 2,
                         regime="evanescent" if V < g.c else "propagating", f0=float("nan"))
                r.update(shield_ratios(st, V))
                r.update(band_T(st, V))
                rows.append(r)
        for pr in P_RATIOS:
            for d in DS:
                st, g = vac_stack(pr, d)
                mf = (RHO20 - g.rho) * d / 2
                Kg, f0 = f0_mass_air_mass(pr, d, V, mf)
                r = dict(stack=f"vac {pr:g} {d:g} m", fill="air", p_ratio=pr, d=d, M=M, V=V,
                         face_mass=mf, regime="propagating", f0=f0, K_gas=Kg,
                         f_cav1=C20 / (2 * d * math.sqrt(1 - C20 ** 2 / V ** 2)))
                r.update(shield_ratios(st, V))
                r.update(band_T(st, V))
                rows.append(r)
    return rows


def ties_tension_sweep():
    rows = []
    d, pr = 30.0, 1e-3
    for M in MACHS:
        V = M * C20
        for K in (0, 30, 100, 300, 1000, 3000, 1e4, 1e5):
            st, _ = vac_stack(pr, d, tie_K=K)
            r = dict(sweep="tie", stack=f"vac 1e-3 30 m, K={K:g}", M=M, V=V, tie_K=K, N=0.0)
            r.update(shield_ratios(st, V))
            r.update(band_T(st, V))
            rows.append(r)
            st, _ = gas_stack("hydrogen", 30.0, tie_K=K)
            r = dict(sweep="tie", stack=f"H2 30 m, K={K:g}", M=M, V=V, tie_K=K, N=0.0)
            r.update(shield_ratios(st, V))
            r.update(band_T(st, V))
            rows.append(r)
        mf = (RHO20 * (1 - pr)) * d / 2
        for N in (1e4, 3e4, 1e5, 2e5, 3e5, 5e5, 1e6):
            for K in (0.0, K_RING_30[0]):
                st, _ = vac_stack(pr, d, N=N, eta=0.02, tie_K=K)
                r = dict(sweep="tension", stack=f"vac 1e-3 30 m, N={N:g}, K={K:g}", M=M, V=V,
                         tie_K=K, N=N, eta=0.02, c_m=math.sqrt(N / mf))
                r.update(shield_ratios(st, V))
                rows.append(r)
    return rows


def coincidence_map(K=0.0):
    """Peak ratio (0.1 s N-wave) vs Mach and face tension, vacuum 30 m, lift-closed faces,
    smeared tie K [Pa/m] between the faces."""
    pr, d = 1e-3, 30.0
    mf = RHO20 * (1 - pr) * d / 2
    Ms = np.linspace(1.2, 3.6, 21)
    Ns = np.logspace(3, 6.3, 19)
    Z = np.zeros((Ns.size, Ms.size))
    t, p = nwave(0.1, pad=3.0)
    for i, N in enumerate(Ns):
        for j, M in enumerate(Ms):
            st, _ = vac_stack(pr, d, N=N, eta=0.02, tie_K=K)
            y = S.apply(st, t, p, M * C20, ATM20, pad_factor=2)
            Z[i, j] = np.max(y)
    return Ms, Ns, mf, Z


def baselines():
    """Fair baselines for the ring-coupled raft (d_eff = pi/4 * 30 m, lift-closed faces):
    the same faces with ambient air between (an ordinary double wall, no lift) and one
    membrane of the same total mass (mass law)."""
    rows = []
    d = math.pi / 4 * 30
    b = RHO20 * d
    air = S.gas("air", T20, P20)
    for M in MACHS:
        V = M * C20
        cases = [("air-filled double wall (same faces)",
                  S.Stack([S.Membrane(b / 2), S.GasLayer(air, d), S.Membrane(b / 2)])),
                 ("mass law (one membrane, same mass)", S.Stack([S.Membrane(b)])),
                 ("vacuum raft D30, ring K=575", vac_stack(1e-3, d, tie_K=K_RING_30[0])[0]),
                 ("vacuum raft D30, ring K=1569", vac_stack(1e-3, d, tie_K=K_RING_30[1])[0]),
                 ("vacuum raft D30, no structure", vac_stack(1e-3, d)[0]),
                 ("vacuum lens T30, struts K=7.7e4", vac_stack(1e-3, 30.0, tie_K=K_STRUT_30)[0]),
                 ("single tube D100, ring K=166", vac_stack(1e-3, math.pi / 4 * 100, tie_K=166.0)[0]),
                 ("single tube D100, ring K=454", vac_stack(1e-3, math.pi / 4 * 100, tie_K=454.0)[0]),
                 ("single tube D250, ring K=63", vac_stack(1e-3, math.pi / 4 * 250, tie_K=63.0)[0]),
                 ("single tube D250, ring K=173", vac_stack(1e-3, math.pi / 4 * 250, tie_K=173.0)[0])]
        for lab, st in cases:
            r = dict(stack=lab, M=M, V=V)
            r.update(shield_ratios(st, V))
            rows.append(r)
    return rows


# --------------------------------------------------------------------- part 4: ground
H_ABOVE = 250.0
RISE_SIGMA_S = 1.0e-3   # Gate-0 proxy (copied): Gaussian sigma = 1 ms
METRIC_TAPER_S = 0.05
TAPER_FRAC = 0.5
G_FS = 50e3
AIRCRAFT = [("Concorde", 1.4), ("Concorde", 2.0), ("Concorde", 2.5), ("SR-71", 3.2)]


def ground_cases(which="main"):
    if which == "large":
        stacks = [("none", {})]
        for D, Ks in ((100, (166.0, 454.0)), (250, (63.0, 173.0))):
            for K in Ks:
                stacks.append((f"single tube D{D} rings K={K:g}",
                               dict(kind="vac", pr=1e-3, d=math.pi / 4 * D, K=K)))
        cases = []
        for ac, M in AIRCRAFT:
            for lab, spec in stacks:
                cases.append(dict(aircraft=ac, mach=M, label=lab, **spec))
        for i, c in enumerate(cases):
            c["id"] = i
        return cases
    stacks = [
        ("none", {}),
        ("He 15 m", dict(kind="gas", fill="helium", d=15.0)),
        ("H2 30 m", dict(kind="gas", fill="hydrogen", d=30.0)),
        ("vac 1e-3 15 m", dict(kind="vac", pr=1e-3, d=15.0)),
        ("vac 1e-3 30 m", dict(kind="vac", pr=1e-3, d=30.0)),
        ("vac 1e-3 60 m", dict(kind="vac", pr=1e-3, d=60.0)),
        ("vac 1e-2 30 m", dict(kind="vac", pr=1e-2, d=30.0)),
        ("vac 1e-1 30 m", dict(kind="vac", pr=1e-1, d=30.0)),
        ("raft D30 rings K=575", dict(kind="vac", pr=1e-3, d=math.pi / 4 * 30, K=K_RING_30[0])),
        ("raft D30 rings K=1569", dict(kind="vac", pr=1e-3, d=math.pi / 4 * 30, K=K_RING_30[1])),
        ("lens T30 struts K=7.7e4", dict(kind="vac", pr=1e-3, d=30.0, K=K_STRUT_30)),
    ]
    cases = []
    for ac, M in AIRCRAFT:
        for lab, spec in stacks:
            cases.append(dict(aircraft=ac, mach=M, label=lab, **spec))
    for i, c in enumerate(cases):
        c["id"] = i
    return cases


def run_ground_case(c):
    from boomsim.loudness import metrics
    from boomsim.propagate import ground_signature
    from boomsim.source import CONCORDE, SR71
    t0 = time.time()
    base = {"Concorde": CONCORDE, "SR-71": SR71}[c["aircraft"]]
    ac = dataclasses.replace(base, mach=float(c["mach"]), altitude=H_S + H_ABOVE)
    V = ac.V
    r0 = min(1.5 * ac.length, 0.8 * H_ABOVE)
    kind = c.get("kind")
    if kind == "gas":
        st, g = gas_stack(c["fill"], c["d"])
        a = math.sqrt(max(0.0, 1 / V ** 2 - 1 / g.c ** 2))
        pad = min(3.0, max(1.0, 15.0 * a * c["d"]))
    elif kind == "vac":
        st, g = vac_stack(c["pr"], c["d"], tie_K=c.get("K", 0.0))
        pad = 3.0   # propagating layer: long, weak cavity ringing; generous pad
    else:
        st, pad = None, 1.0
    cap = {}

    def insert(t, p):
        dt = t[1] - t[0]
        n = int(round(pad / dt))
        p2 = np.concatenate([np.zeros(n), p, np.zeros(n)])
        t2 = t[0] + (np.arange(p2.size) - n) * dt
        out = S.apply(st, t2, p2, V, ATM20) if st is not None else p2.copy()
        nt = max(2, int(TAPER_FRAC * n))
        w = np.ones(out.size)
        ramp = 0.5 * (1 - np.cos(np.pi * np.arange(nt) / nt))
        w[:nt] = ramp
        w[-nt:] = ramp[::-1]
        out = out * w
        cap["above"], cap["below"] = (t2, p2), (t2, out)
        return t2, out

    t, p = ground_signature(ac, r0=r0, insert=(H_S, insert), fs=G_FS)
    m = metrics(t, p, taper_s=METRIC_TAPER_S)
    ps = gaussian_filter1d(p, RISE_SIGMA_S * G_FS, mode="constant")
    ms = metrics(t, ps, taper_s=METRIC_TAPER_S)
    ta, pa = cap["above"]
    tb, pb = cap["below"]
    res = {k: v for k, v in c.items()}
    res.update(V=V, pldb=m["pldb"], pldb_proxy=ms["pldb"], peak_pa=m["peak_pa"],
               rise_ms=m["rise_time_10_90_ms"], csel=m["csel"], pos_impulse=m["pos_impulse_pa_s"],
               shield_peak_ratio=float(pb.max() / pa.max()),
               shield_imp_ratio=pos_impulse(tb, pb) / pos_impulse(ta, pa),
               runtime_s=time.time() - t0)
    keep = (t > -0.2) & (t < 1.2)
    res["_wave"] = (t[keep][::10].astype(np.float32), p[keep][::10].astype(np.float32))
    return res


def ground_runs(procs, which="main"):
    cases = ground_cases(which)
    with Pool(procs) as pool:
        rows = pool.map(run_ground_case, cases, chunksize=1)
    ref = {(r["aircraft"], r["mach"]): r for r in rows if r["label"] == "none"}
    for r in rows:
        u = ref[(r["aircraft"], r["mach"])]
        r["dpldb_raw"] = r["pldb"] - u["pldb"]
        r["dpldb_proxy"] = r["pldb_proxy"] - u["pldb_proxy"]
        r["dpldb_headline"] = max(r["dpldb_raw"], r["dpldb_proxy"])
        r["dcsel"] = r["csel"] - u["csel"]
        r["ground_peak_ratio"] = r["peak_pa"] / u["peak_pa"]
    return rows


# --------------------------------------------------------------------- output
def write_csv(rows, path, cols=None):
    cols = cols or sorted({k for r in rows for k in r if not k.startswith("_")})
    with open(path, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols, extrasaction="ignore")
        w.writeheader()
        for r in rows:
            w.writerow({k: (f"{r[k]:.5g}" if isinstance(r.get(k), float) else r.get(k, "")) for k in cols})


def plot_transfer(path):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    ink2, muted, grid, surf = "#52514e", "#8a8984", "#e4e3df", "#fcfcfb"
    f = np.logspace(-1, 3, 4000)
    fig, axes = plt.subplots(1, 2, figsize=(12, 4.6))
    fig.patch.set_facecolor(surf)
    V = 2.0 * C20
    ax = axes[0]
    for lab, st, c, ls in [
        ("He 15 m (Gate 0)", gas_stack("helium", 15)[0], "#2a78d6", "-"),
        ("H2 30 m (Gate 0)", gas_stack("hydrogen", 30)[0], "#eb6834", "-"),
        ("vacuum 1e-3, 30 m, no structure", vac_stack(1e-3, 30)[0], "#0b0b0b", "-"),
        ("vacuum 1e-1, 30 m, no structure", vac_stack(1e-1, 30)[0], "#8a8984", "-"),
        ("vacuum raft D30, ring K=575", vac_stack(1e-3, math.pi / 4 * 30, tie_K=K_RING_30[0])[0], "#1baf7a", "--"),
        ("vacuum lens T30, struts K=7.7e4", vac_stack(1e-3, 30, tie_K=K_STRUT_30)[0], "#7a4fd0", "--"),
        ("single membrane, same mass (mass law)", S.Stack([S.Membrane(RHO20 * 30)]), "#c0392b", ":"),
    ]:
        T, _ = S.transfer(st, f, V, ATM20)
        ax.semilogx(f, 20 * np.log10(np.maximum(np.abs(T), 1e-12)), color=c, ls=ls, lw=1.4, label=lab)
    ax.set_ylim(-160, 10)
    ax.set_xlabel("frequency [Hz]")
    ax.set_ylabel("|T| [dB]")
    ax.set_title("Transmission at Mach 2.0 (lift-closed faces)", fontsize=10)
    ax.legend(frameon=False, fontsize=7.5, loc="lower left")
    ax = axes[1]
    for K, c in [(0, "#0b0b0b"), (100, "#2a78d6"), (1000, "#1baf7a"), (1e4, "#eb6834"), (1e5, "#7a4fd0")]:
        T, _ = S.transfer(vac_stack(1e-3, 30, tie_K=K)[0], f, V, ATM20)
        ax.semilogx(f, 20 * np.log10(np.maximum(np.abs(T), 1e-12)), color=c, lw=1.4,
                    label=f"vacuum 1e-3, 30 m, tie K = {K:g} Pa/m")
    T, _ = S.transfer(S.Stack([S.Membrane(RHO20 * 30)]), f, V, ATM20)
    ax.semilogx(f, 20 * np.log10(np.abs(T)), color="#c0392b", ls=":", lw=1.4, label="mass law, same mass")
    ax.set_ylim(-160, 10)
    ax.set_xlabel("frequency [Hz]")
    ax.set_title("Structure-borne coupling (smeared ties), Mach 2.0", fontsize=10)
    ax.legend(frameon=False, fontsize=7.5, loc="lower left")
    for ax in axes:
        ax.set_facecolor(surf)
        ax.grid(True, color=grid, lw=0.6)
        for s in ("top", "right"):
            ax.spines[s].set_visible(False)
    fig.tight_layout()
    fig.savefig(path, dpi=130, facecolor=surf)
    plt.close(fig)


def plot_coincidence(maps, path):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    surf = "#fcfcfb"
    fig, axes = plt.subplots(1, len(maps), figsize=(7.0 * len(maps), 5))
    fig.patch.set_facecolor(surf)
    for ax, (K, (Ms, Ns, mf, Z)) in zip(np.atleast_1d(axes), maps):
        _coinc_panel(fig, ax, K, Ms, Ns, mf, Z)
    fig.tight_layout()
    fig.savefig(path, dpi=130, facecolor=surf)
    plt.close(fig)


def _coinc_panel(fig, ax, K, Ms, Ns, mf, Z):
    im = ax.pcolormesh(Ms, Ns, 20 * np.log10(Z), shading="auto", cmap="viridis", vmin=-60, vmax=0)
    ax.set_yscale("log")
    cm = np.sqrt(Ns / mf)
    ax.plot(cm / C20, Ns, color="white", lw=1.5, ls="--", label="c_m = V (coincidence)")
    for N, lab in [(1.25 * P20 * 3, "tube tendons, ring pitch 3 m"), (1.25 * P20 * 12, "ring pitch 12 m")]:
        ax.axhline(N, color="#eb6834", lw=1.0)
        ax.text(1.22, N * 1.08, lab, color="#eb6834", fontsize=7.5)
    ax.set_xlim(Ms[0], Ms[-1])
    ax.set_xlabel("Mach (trace speed V = M·295 m/s)")
    ax.set_ylabel("along-track face tension N [N/m]")
    ax.set_title(f"Vacuum 1e-3, 30 m, faces m = {mf:.2f} kg/m² each, η = 0.02, tie K = {K:g} Pa/m\n"
                 "peak ratio of 0.1 s N-wave at the shield [dB]", fontsize=9.5)
    fig.colorbar(im, ax=ax, label="20 log10(peak ratio)")
    ax.legend(frameon=False, fontsize=7.5, loc="lower right", labelcolor="white")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--ground", action="store_true")
    ap.add_argument("--procs", type=int, default=4)
    ap.add_argument("--only-baselines", action="store_true")
    ap.add_argument("--ground-large", action="store_true",
                    help="only the large single-tube ground set -> acoustics_ground_large.csv")
    args = ap.parse_args()

    if args.ground_large:
        g = ground_runs(args.procs, "large")
        _write_ground(g, "acoustics_ground_large.csv")
        return
    bl = baselines()
    write_csv(bl, os.path.join(HERE, "acoustics_baselines.csv"))
    print("baselines (0.1 s N-wave at the shield):")
    for r in bl:
        print(f"  {r['stack']:38s} M{r['M']:.1f} peak {r['peak_0.1']:.3f} imp {r['imp_0.1']:.3f}")
    if args.only_baselines:
        return

    print(f"20 km: T {T20:.2f} K, p {P20:.0f} Pa, rho {RHO20:.4f}, c {C20:.1f} m/s")
    for pr in (1e-1, 1e-2, 1e-3, 1e-4):
        g = gas_checks(pr)
        print("  " + ", ".join(f"{k}={v:.3g}" for k, v in g.items()))

    rows = shield_plane_sweep()
    write_csv(rows, os.path.join(HERE, "acoustics_shield_plane.csv"))
    print("\nshield plane (lift-closed faces): peak/impulse ratios, 0.1 s and 0.3 s N-waves")
    print(f"{'stack':18s} {'M':>4s} {'mface':>6s} {'f0':>6s} {'fcav':>6s} {'pk.1':>6s} {'imp.1':>6s} {'pk.3':>6s} {'imp.3':>6s} {'T10':>7s} {'T50':>7s}")
    for r in rows:
        print(f"{r['stack']:18s} {r['M']:4.1f} {r['face_mass']:6.2f} {r.get('f0', np.nan):6.2f} "
              f"{r.get('f_cav1', np.nan):6.1f} {r['peak_0.1']:6.3f} {r['imp_0.1']:6.3f} "
              f"{r['peak_0.3']:6.3f} {r['imp_0.3']:6.3f} {r['T_10Hz_dB']:7.1f} {r['T_50Hz_dB']:7.1f}")

    tt = ties_tension_sweep()
    write_csv(tt, os.path.join(HERE, "acoustics_ties_tension.csv"))
    print("\nties / tension sweep")
    for r in tt:
        print(f"{r['stack']:38s} M{r['M']:.1f} pk.1 {r['peak_0.1']:.3f} imp.1 {r['imp_0.1']:.3f} "
              f"pk.3 {r['peak_0.3']:.3f} T10 {r.get('T_10Hz_dB', np.nan):7.1f} T50 {r.get('T_50Hz_dB', np.nan):7.1f}"
              + (f" c_m {r['c_m']:.0f}" if 'c_m' in r else ""))
    plot_transfer(os.path.join(HERE, "acoustics_transfer.png"))
    maps = [(K, coincidence_map(K)) for K in (0.0, K_RING_30[0])]
    plot_coincidence(maps, os.path.join(HERE, "acoustics_coincidence.png"))

    if args.ground:
        try:
            import boomsim.propagate  # noqa: F401
        except Exception as e:  # noqa: BLE001
            print(f"ground runs skipped: boomsim.propagate import failed ({e!r})")
            return
        t0 = time.time()
        g = ground_runs(args.procs)
        print(f"\nground runs ({time.time() - t0:.0f} s)")
        _write_ground(g, "acoustics_ground.csv")


def _write_ground(g, name):
    cols = ["aircraft", "mach", "V", "label", "pr", "d", "K", "pldb", "pldb_proxy", "dpldb_raw",
            "dpldb_proxy", "dpldb_headline", "dcsel", "peak_pa", "ground_peak_ratio", "rise_ms",
            "shield_peak_ratio", "shield_imp_ratio", "runtime_s"]
    write_csv(g, os.path.join(HERE, name), cols)
    for r in g:
        print(f"{r['aircraft']:8s} M{r['mach']:.1f} {r['label']:26s} PLdB {r['pldb']:6.1f}/{r['pldb_proxy']:6.1f} "
              f"dPL {r['dpldb_headline']:+7.1f} dCSEL {r['dcsel']:+6.1f} shield pk {r['shield_peak_ratio']:.3f} "
              f"imp {r['shield_imp_ratio']:.3f} gnd pk {r['ground_peak_ratio']:.3f} rise {r['rise_ms']:.1f} ms")


if __name__ == "__main__":
    main()
