"""Independent check of discrete line ties vs smeared ties (results/line_ties_check.md).

Semi-analytical Floquet-Bloch model of a laterally infinite membrane / gas / membrane
stack (crossplane conventions: k_x = omega/V, cross-plane incidence k_y = 0, z down
positive velocities, e^{+j omega t}) with line ties periodic in y at pitch P.

Per Floquet order n (k_n = 2 pi n / P, horizontal slowness s_n = sqrt(1/V^2 + k_n^2/w^2)):

    [Zt + Za + Zg11,   -Zg12     ] [v_t]   [2 p_i delta_n0 + f_t,n]
    [-Zg12,       Zb + Za + Zg11 ] [v_b] = [f_b,n                 ]

Z_t,b: face impedance (shield.Membrane(m, N).impedance(w, s_n): isotropic tension, so
the y-tension term N k_n^2 is included), Za = rho w / k_z,air (radiation, both sides),
gas 2-port Zg11 = j rho_g w coth(kappa d)/kappa, Zg12 = j rho_g w /(kappa sinh(kappa d)).
Tie at y = 0 of each cell, line stiffness k_l = K P [N/m per m of tie line], attached
over a strip of width b: force F per unit x-length uniform over the strip (Fourier
weight S_n/P, S_n = sinc(k_n b / 2)), elongation = strip average of w_b - w_t.
Eliminating F (Sherman-Morrison over all orders):

    F = -s_l h_0 / (1 + (s_l / P) sum_n S_n^2 g_n),   s_l = k_l / (j w),
    g_n = u^T A_n^{-1} u,  u = (1, -1),  h_0 = u^T A_0^{-1} (2 p_i, 0).

Only n = 0 in the sum gives exactly the smeared GasLayer(tie_K = K) result.

Run:  python boom-shield/sim/check_lineties.py [floquet|bem|all]
"""
from __future__ import annotations

import csv
import math
import os
import pickle
import sys
import time

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from boomsim import shield as SH  # noqa: E402
from boomsim.atmosphere import std_atm  # noqa: E402

RES = os.path.join(HERE, "results")
CACHE = os.path.join(RES, "3d_cache")
H_S = 20000.0
ATM = std_atm(H_S)
T20, P20, RHO20 = float(ATM.T), float(ATM.p), float(ATM.rho)
C20 = float(ATM.c)


def setup(fill="hydrogen", d=30.0, mach=2.0):
    g = SH.gas(fill, T20, P20)
    m = (RHO20 - g.rho) * d / 2.0          # lift-closed faces (Gate 0 / 3-D study)
    return g, m, mach * C20


# ------------------------------------------------------------------- Floquet model
def _orders(f, V, g, d, m, N, ky, eta=0.0):
    """Per-order matrices for frequency f (scalar) and wavenumbers ky (array).
    Returns (A11, A12, A22, Za) arrays."""
    w = 2 * np.pi * f
    s = np.sqrt(1.0 / V ** 2 + (ky / w) ** 2)
    kza = w * np.sqrt((1.0 / C20 ** 2 - s ** 2) + 0j)
    kza = np.where(kza.imag > 0, -kza, kza)
    Za = RHO20 * w / kza
    a2 = s ** 2 - 1.0 / g.c ** 2
    if np.any(a2 <= 0):
        raise ValueError("gas not evanescent for every order")
    kap = w * np.sqrt(a2)
    x = kap * d
    e = np.exp(-2 * x)
    coth = (1 + e) / (1 - e)
    csch = 2 * np.exp(-x) / (1 - e)
    Zg11 = 1j * g.rho * w * coth / kap
    Zg12 = 1j * g.rho * w * csch / kap
    Zm = SH.Membrane(m, N, eta).impedance(np.full(ky.shape, w), s)
    return Zm + Za + Zg11, -Zg12, Zm + Za + Zg11, Za


def floquet(f, V, g, d, m, K, P, b=None, N=0.0, n_max=None, eta=0.0, smeared=False):
    """Specular T_0 (pressure below / incident, interface reference), total transmitted
    power fraction tau_tot, and the effective tie stiffness ratio K_eff/K (= F_line /
    (K P e_0) with e_0 the n = 0 elongation response). f scalar.
    b: tie attachment width (default P/1000); smeared=True keeps n = 0 only."""
    b = P * 1e-3 if b is None else float(b)
    w = 2 * np.pi * f
    if smeared or K == 0:
        n = np.array([0])
    else:
        if n_max is None:
            n_max = int(min(4e5, max(400, 60 * P / b)))
            if N > 0:   # tension: g_n ~ 1/n^2, sum converges fast once N k^2 >> w^2 m
                n_max = int(min(n_max, max(400, 40 * P * w * math.sqrt(m / N) + 400)))
        n = np.arange(-n_max, n_max + 1)
    ky = 2 * np.pi * n / P
    A11, A12, A22, Za = _orders(f, V, g, d, m, N, ky, eta)
    det = A11 * A22 - A12 ** 2
    # A^{-1} = [[A22, -A12], [-A12, A11]] / det
    gn = (A22 + A11 + 2 * A12) / det                  # u^T A^-1 u, u = (1, -1)
    S = np.sinc(ky * b / (2 * np.pi))                 # sin(k b/2)/(k b/2)
    S = np.ones_like(S) if smeared else S
    i0 = int(np.nonzero(n == 0)[0][0])
    pi2 = 2.0
    v0 = np.array([A22[i0], -A12[i0]]) / det[i0] * pi2   # A0^-1 (2 p_i, 0)
    h0 = v0[0] - v0[1]
    sl = K * P / (1j * w)                              # line stiffness / (j w)
    Ssum = np.sum(S ** 2 * gn)
    Fl = -sl * h0 / (1 + sl / P * Ssum)                # tie force per unit x-length
    # order velocities: v_n = A_n^-1 (2 p_i delta + F S_n/P u)
    fn = Fl * S / P
    vb = (-A12 * fn - A11 * fn) / det                 # bottom row of A^-1 applied to (fn, -fn)
    vb[i0] += v0[1]
    Tn = Za * vb
    prop = np.abs(Za.imag) < 1e-9 * np.abs(Za)
    kz = RHO20 * w / Za
    tau = np.sum(np.abs(Tn[prop]) ** 2 * (kz[prop].real / kz[i0].real))
    Feq = -sl * h0 / (1 + sl / P * gn[i0] * 1.0)       # smeared-equivalent tie force
    return Tn[i0], tau, Fl / Feq


def floquet_band(fs, *a, **k):
    out = np.array([floquet(f, *a, **k) for f in fs])
    return out[:, 0], out[:, 1].real, out[:, 2]


def smeared_transfer(fs, V, g, d, m, K, N=0.0, eta=0.0):
    st = SH.Stack([SH.Membrane(m, N, eta), SH.GasLayer(g, d, tie_K=K), SH.Membrane(m, N, eta)])
    return SH.transfer(st, np.asarray(fs, float), V, ATM)[0]


def band_db(fs, tau, tau_ref, lo=5.0, hi=30.0, weight="nwave"):
    """10 log10 of band power of tau relative to tau_ref, weight 1/f^2 (N-wave spectrum)."""
    m = (fs >= lo) & (fs <= hi)
    wgt = 1.0 / fs[m] ** 2 if weight == "nwave" else np.ones(m.sum())
    return 10 * np.log10(np.sum(wgt * tau[m]) / np.sum(wgt * tau_ref[m]))


# ------------------------------------------------------------------- sweeps
FREQ = np.round(np.arange(1.0, 100.001, 0.1), 3)
PITCHES = (5.0, 12.5, 25.0, 50.0)
KS = (3.0, 30.0, 300.0, 3000.0)
TENSIONS = (0.0, 1e3, 1e4, 1e5)
B_PHYS = 0.5          # nominal tie attachment width [m]
ETA = 0.05            # nominal face loss factor (acts on the tension term only)


def _free(fill, mach, N=0.0, eta=ETA, fs=FREQ):
    g, m, V = setup(fill, 30.0, mach)
    return np.abs(smeared_transfer(fs, V, g, 30.0, m, 0.0, N, eta)) ** 2


def sweep():
    path = os.path.join(CACHE, "lineties_floquet.pkl")
    if os.path.exists(path):
        with open(path, "rb") as fh:
            return pickle.load(fh)
    t0 = time.time()
    out = {}
    for fill in ("helium", "hydrogen"):
        for mach in (2.0, 1.4):
            g, m, V = setup(fill, 30.0, mach)
            for N in TENSIONS:
                out[(fill, mach, N, "free")] = _free(fill, mach, N)
            for K in KS:
                for N in TENSIONS:
                    out[(fill, mach, K, N, "smeared")] = np.abs(smeared_transfer(FREQ, V, g, 30.0, m, K, N, ETA)) ** 2
                    for P in PITCHES:
                        T0, tau, r = floquet_band(FREQ, V, g, 30.0, m, K, P, b=B_PHYS, N=N, eta=ETA)
                        out[(fill, mach, K, N, P)] = (np.abs(T0) ** 2, tau, r)
            print(f"  {fill} M{mach}: {time.time() - t0:.0f} s", flush=True)
    os.makedirs(CACHE, exist_ok=True)
    with open(path, "wb") as fh:
        pickle.dump(out, fh)
    return out


def band_loss(fill, mach, K, P, N, b=B_PHYS, eta=ETA, fs=None):
    """5-30 Hz band loss [dB] re the tie-free stack (same N): 1/f^2-weighted total
    transmitted power. P = 0 -> smeared ties."""
    fs = np.arange(5.0, 30.001, 0.1) if fs is None else fs
    g, m, V = setup(fill, 30.0, mach)
    ref = _free(fill, mach, N, eta, fs)
    if P == 0:
        tau = np.abs(smeared_transfer(fs, V, g, 30.0, m, K, N, eta)) ** 2
    else:
        tau = floquet_band(fs, V, g, 30.0, m, K, P, b=b, N=N, eta=eta)[1]
    return band_db(fs, tau, ref)


def allowable_K(fill, mach, P, N, b=B_PHYS, eta=ETA, loss_db=3.0):
    """Smallest smeared-equivalent K [Pa/m] at which the 5-30 Hz band loss first exceeds
    loss_db (scan upward on a log grid, then bisect); ties at or below it are safe."""
    grid = np.logspace(-0.5, 5.0, 34)
    prev = grid[0]
    for K in grid:
        if band_loss(fill, mach, K, P, N, b, eta) > loss_db:
            lo, hi = math.log10(prev), math.log10(K)
            for _ in range(8):
                mid = 0.5 * (lo + hi)
                if band_loss(fill, mach, 10 ** mid, P, N, b, eta) > loss_db:
                    hi = mid
                else:
                    lo = mid
            return 10 ** lo
        prev = K
    return float("inf")


def rule_table():
    path = os.path.join(CACHE, "lineties_rule.pkl")
    out = {}
    if os.path.exists(path):
        with open(path, "rb") as fh:
            out = pickle.load(fh)
    t0 = time.time()
    for fill in ("hydrogen", "helium"):
        for mach in (2.0, 1.4):
            for N in TENSIONS:
                for P in (0.0,) + PITCHES:
                    key = (fill, mach, N, P)
                    if key in out:
                        continue
                    out[key] = allowable_K(fill, mach, P, N)
                    print(f"  allow {key}: {out[key]:.3g} Pa/m ({time.time() - t0:.0f} s)", flush=True)
                    with open(path + ".tmp", "wb") as fh:
                        pickle.dump(out, fh)
                    os.replace(path + ".tmp", path)
    return out


# ------------------------------------------------------------------- BEM comparison
def bem_T0(f, K=300.0, P=25.0, kind="line", b=None, h=None, N=0.0, eta=0.0, W=1000.0,
           fill="hydrogen", mach=2.0):
    """Cross-plane BEM (boomsim.crossplane) of a W-wide 30 m slab (Concorde geometry,
    h 250 m) under a plane wave k_y = 0 from above. Returns (T_0, b_used, n_panels):
    T_0 = panel-length-weighted mean transmitted pressure on the bottom face over
    |y| <= 250 m (an integer number of periods) / incident pressure at the top face.
    kind: "none" | "smeared" | "line" (line ties at pitch P, k_line = K P, width b).
    N > 0: face tension, x part through shield.Membrane, y part added here as a
    finite-difference N w_yy coupling between neighbouring face panels (w = 0 at the
    face ends)."""
    import dataclasses
    import warnings
    from boomsim import crossplane as cp
    from boomsim.source import CONCORDE
    g, m, V = setup(fill, 30.0, mach)
    ac = dataclasses.replace(CONCORDE, mach=mach, altitude=H_S + 250.0)
    face = SH.Membrane(m, N, eta) if N > 0 else m
    kw = {}
    if kind == "smeared":
        kw = dict(tie_K=K)
    elif kind == "line":
        kw = dict(line_ties=[(float(y), K * P) for y in np.arange(-W / 2 + P / 2, W / 2, P)],
                  line_tie_width=b)
    cfg = cp.CrossPlaneConfig(ac, h=250.0, W=W, T=30.0, gas=g, m_top=face, m_bot=face,
                              f_max=100.0, phis=[0.0], workers=1, **kw)
    cfg.n_dense = 2000
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", RuntimeWarning)
        mesh = cp.build_mesh(cfg, f, h_elem=h)
        if len(mesh["P"]) > cfg.n_dense:
            mesh = cp.build_mesh(cfg, f, h_elem=h, uniform=True)
    Pn, nf = mesh["P"], mesh["nf"]
    yn = np.r_[Pn.a[:nf, 0], Pn.b[nf - 1, 0]]
    if N > 0:   # lateral tension: force N w_yy on each face panel (v along n = down)
        Nc = N * (1 + 1j * eta) / (1j * mesh["w"])
        c = Pn.c[:nf, 0]
        Ls = Pn.L[:nf]
        hp = np.r_[np.diff(c), 0.5 * Ls[-1]]
        hm = np.r_[0.5 * Ls[0], np.diff(c)]
        for off in (0, nf):
            for i in range(nf):
                mesh["ties"].append((off + i, off + i, -Nc * (1 / hp[i] + 1 / hm[i]) / Ls[i]))
                if i + 1 < nf:
                    mesh["ties"].append((off + i, off + i + 1, Nc / hp[i] / Ls[i]))
                if i > 0:
                    mesh["ties"].append((off + i, off + i - 1, Nc / hm[i] / Ls[i]))
    k = mesh["k"]
    inc = lambda Y, Z: (np.exp(1j * k * Z), 0 * Y + 0j, 1j * k * np.exp(1j * k * Z))
    sol = cp.bem_solve(cfg, f, mesh=mesh, inc=inc)
    ib = nf + np.arange(nf)
    sel = np.abs(Pn.c[ib, 0]) <= 250.0
    L = Pn.L[ib][sel]
    T0 = np.sum(L * sol.p_ext[ib][sel]) / np.sum(L) / np.exp(1j * k * (-250.0))
    i = int(np.clip(np.searchsorted(yn, 0.5 * P) - 1, 0, nf - 1))
    b_used = b if b is not None else float(Pn.L[i])
    return T0, b_used, len(Pn)


def bem_compare(freqs=(2, 3, 4, 5, 6, 7, 8, 9, 10, 12, 14, 20), K=300.0, P=25.0):
    """Cached BEM runs: legacy (one panel) at several element sizes, resolved widths,
    and tensioned faces."""
    path = os.path.join(CACHE, "lineties_bem.pkl")
    res = {}
    if os.path.exists(path):
        with open(path, "rb") as fh:
            res = pickle.load(fh)
    runs = [("none", None, None, 0.0), ("smeared", None, None, 0.0)]
    runs += [("line", None, h, 0.0) for h in (None, 6.0, 3.0, 1.5, 0.75)]
    runs += [("line", 1.5, h, 0.0) for h in (None, 3.0)]
    runs += [("line", 1.5, h, 1e4) for h in (3.0, 1.5)]
    runs += [("none", None, 3.0, 1e4), ("smeared", None, 3.0, 1e4)]
    for (kind, b, h, N) in runs:
        for f in freqs:
            key = (kind, b, h, N, float(f), K, P)
            if key in res:
                continue
            if N > 0 and f > 14:
                continue
            t0 = time.time()
            T0, bu, n = bem_T0(f, K, P, kind, b, h, N, ETA if N > 0 else 0.0)
            res[key] = dict(T0=T0, b=bu, n=n)
            print(f"  BEM {key}: {10*np.log10(abs(T0)**2):6.1f} dB b={bu:.2f} {n} panels "
                  f"{time.time()-t0:.1f} s", flush=True)
            with open(path + ".tmp", "wb") as fh:
                pickle.dump(res, fh)
            os.replace(path + ".tmp", path)
    return res


# ------------------------------------------------------------------- report
def report():
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    o = sweep()
    rule = rule_table()
    bem = bem_compare()
    fs = FREQ
    # ---- CSV: band losses and allowable K
    with open(os.path.join(RES, "line_ties_check.csv"), "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["fill", "mach", "N_N_per_m", "K_Pa_per_m", "pitch_m (0 = smeared)",
                    "band_loss_5_30Hz_dB", "allowable_K_Pa_per_m (3 dB)"])
        for fill in ("hydrogen", "helium"):
            for mach in (2.0, 1.4):
                for N in TENSIONS:
                    ref = o[(fill, mach, N, "free")]
                    for K in KS:
                        w.writerow([fill, mach, N, K, 0, round(band_db(fs, o[(fill, mach, K, N, "smeared")], ref), 2),
                                    round(rule[(fill, mach, N, 0.0)], 2)])
                        for P in PITCHES:
                            w.writerow([fill, mach, N, K, P, round(band_db(fs, o[(fill, mach, K, N, P)][1], ref), 2),
                                        round(rule[(fill, mach, N, P)], 2)])
    # ---- figure
    fig, ax = plt.subplots(1, 3, figsize=(16, 4.8))
    db = lambda x: 10 * np.log10(np.abs(x) ** 2)
    fb = sorted({k[4] for k in bem if k[0] == "line" and k[1] is None and k[2] is None})
    g, m, V = setup("hydrogen", 30.0, 2.0)
    ff = np.arange(1.0, 25.01, 0.05)
    ax[0].plot(ff, db(smeared_transfer(ff, V, g, 30.0, m, 0.0)), "k:", label="no ties")
    ax[0].plot(ff, db(smeared_transfer(ff, V, g, 30.0, m, 300.0)), "k-", label="smeared K = 300")
    cols = {None: "C3", 3.0: "C1", 0.75: "C2"}
    for h, c in cols.items():
        pts = [(f, bem[("line", None, h, 0.0, f, 300.0, 25.0)]) for f in fb]
        ax[0].plot([p[0] for p in pts], [db(p[1]["T0"]) for p in pts], "o", color=c, ms=5,
                   label=f"BEM legacy, panel {'34/f' if h is None else h} m")
        ax[0].plot([p[0] for p in pts], [db(floquet(p[0], V, g, 30.0, m, 300.0, 25.0, b=p[1]["b"])[0])
                                          for p in pts], "-", color=c, lw=1)
    pts = [(f, bem[("line", 1.5, None, 0.0, f, 300.0, 25.0)]) for f in fb]
    ax[0].plot([p[0] for p in pts], [db(p[1]["T0"]) for p in pts], "s", mfc="none", color="C0", ms=7,
               label="BEM fixed, b = 1.5 m")
    ax[0].plot(ff, [db(floquet(f, V, g, 30.0, m, 300.0, 25.0, b=1.5)[0]) for f in ff], "-", color="C0", lw=1)
    ax[0].set(xlim=(1, 25), ylim=(-70, 5), xlabel="f [Hz]", ylabel="|T_0|^2 [dB]",
              title="H2 30 m M2, 25 m line ties (K = 300), N = 0\nmarkers BEM, lines Floquet (same b)")
    ax[0].legend(fontsize=7, loc="lower left")
    for N, c in zip(TENSIONS, ("C3", "C1", "C2", "C0")):
        ax[1].plot(fs, 10 * np.log10(o[("hydrogen", 2.0, 300.0, N, 25.0)][1]), color=c, lw=1, label=f"25 m ties, N = {N:g}")
    ax[1].plot(fs, 10 * np.log10(o[("hydrogen", 2.0, 300.0, 1e4, "smeared")]), "k-", label="smeared K = 300")
    ax[1].plot(fs, 10 * np.log10(o[("hydrogen", 2.0, 0.0 if False else 1e4, "free")]), "k:", label="no ties")
    ax[1].set(xlim=(1, 40), ylim=(-100, 5), xlabel="f [Hz]", ylabel="total transmitted power [dB]",
              title="Floquet, b = 0.5 m, eta = 0.05: face tension")
    ax[1].legend(fontsize=7, loc="lower left")
    for N, c in zip(TENSIONS, ("C3", "C1", "C2", "C0")):
        ax[2].semilogy(PITCHES, [rule[("hydrogen", 2.0, N, P)] for P in PITCHES], "o-", color=c, label=f"H2 M2, N = {N:g}")
        ax[2].semilogy(PITCHES, [rule[("helium", 1.4, N, P)] for P in PITCHES], "s--", color=c, mfc="none",
                       label=f"He M1.4, N = {N:g}")
    ax[2].axhline(rule[("hydrogen", 2.0, 0.0, 0.0)], color="k", lw=1, label="smeared, H2 M2")
    ax[2].axhline(rule[("helium", 1.4, 0.0, 0.0)], color="k", ls="--", lw=1, label="smeared, He M1.4")
    ax[2].set(xlabel="tie pitch P [m]", ylabel="allowable K = k_line/P [Pa/m]",
              title="<= 3 dB loss of tie-free 5-30 Hz TL (b = 0.5 m)")
    ax[2].legend(fontsize=6, ncol=2)
    fig.tight_layout()
    fig.savefig(os.path.join(RES, "line_ties_check.png"), dpi=130)


if __name__ == "__main__":
    what = sys.argv[1] if len(sys.argv) > 1 else "all"
    if what in ("floquet", "all"):
        sweep()
    if what in ("rule", "all"):
        rule_table()
    if what in ("bem", "all"):
        bem_compare()
    if what in ("report", "all"):
        report()
