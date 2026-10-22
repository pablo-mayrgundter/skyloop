"""Cross-plane BEM of a raft of parallel superpressure H2 tubes (BoomShield structure study).

Uses boomsim.bem2d / boomsim.crossplane read-only. Same physics and conventions as
crossplane.py (steady flight along x, every field ~ exp(j(wt - k_x x)), k_x = w/V, 2-D
problem in the (y, z) cross-plane per frequency; exterior air Helmholtz, H2 interior
modified Helmholtz because c_H2 > V).

Raft model
----------
N tubes of diameter D (radius R) side by side across the width W, axes along the corridor
(x), tops at depth h below the flight axis. Either

* touching (gap = 0): neighbouring tubes are flattened against each other on a contact
  strip of height b (two skins, no gas): a vertical "chord" panel set between the two gas
  domains (mass 2 m_s, both hoop tensions), or
* gapped (gap > 0): separate circles with an air slit of width ``gap`` between them
  (exterior air connects top and bottom through the slit).

Skin (each tube): areal mass m_s [kg/m^2 of skin], hoop tension N = dp R, along-track
tension N_x = dp R / 2 (membrane impedance j w m - j N_x w / V^2, as shield.Membrane),
plus the curved-membrane geometric stiffness in the cross-plane

    f_out = N (w'' + kappa^2 w) - kappa * (E t / C) * sum_i kappa_i w_i L_i

(w outward normal displacement, '' along the skin, kappa = 1/R on arcs and 0 on contact
chords, C the loop length; the second term is hoop stretching with free tangential
sliding, E t the skin's hoop membrane stiffness). For a closed ring this gives the classical
mode stiffness (n^2 - 1) dp / R (n >= 1; n = 2 ovalisation 3 dp/R) and (E t - N)/R^2 for
breathing (n = 0). It enters bem2d through ``Problem.ties`` (force on panel j along n_j =
sum_i val_ji v_i, membrane law p_A - p_B = Z v - f). Verified against a cylinder-harmonic
series with mode-dependent impedance in ``self_test()``.
"""
from __future__ import annotations

import math
import os
import sys

import numpy as np
from scipy import special as sp

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.abspath(os.path.join(HERE, "..", "sim")))

from boomsim import bem2d  # noqa: E402
from boomsim import crossplane as cp  # noqa: E402
from boomsim import shield as SH  # noqa: E402
from boomsim.atmosphere import std_atm  # noqa: E402
from boomsim.source import CONCORDE  # noqa: E402
import dataclasses  # noqa: E402

Z_SHIELD = 20000.0
ATM = std_atm(Z_SHIELD)
RHO, C0, P0, T0 = float(ATM.rho), float(ATM.c), float(ATM.p), float(ATM.T)
H2 = SH.gas("hydrogen", T0, P0)


def aircraft(mach=2.0, h=250.0):
    return dataclasses.replace(CONCORDE, mach=float(mach), altitude=Z_SHIELD + float(h))


# ----------------------------------------------------------------------------- geometry
def _arc(xc, zc, R, th0, th1, h_e):
    n = max(2, int(math.ceil(R * (th1 - th0) / h_e)))
    th = np.linspace(th0, th1, n + 1)
    pts = np.stack([xc + R * np.cos(th), zc + R * np.sin(th)], 1)
    return pts[:-1], pts[1:]


def build_raft(D, W, h, gap=0.0, b=1.0, dp=300.0, Et=1e5, m_s=None, h_e=2.0, V=590.0,
               gas=H2, rho=RHO, n_tubes=None):
    """Panels + bookkeeping for a raft. Returns dict(P, side_a, side_b, kind, Z(w), loops...).

    m_s None -> lift-closed skin: (rho - rho_g) R / 2 per m^2 of skin (all lift in the skin).
    """
    R = 0.5 * D
    if m_s is None:
        m_s = max((rho - gas.rho) * R / 2.0, 1e-3)
    if gap > 0:
        s = D + gap
        b = 0.0
    else:
        s = 2.0 * math.sqrt(R * R - (0.5 * b) ** 2)
    if n_tubes is None:
        n_tubes = max(1, int(math.floor((W - D) / s + 1e-9)) + 1)
    yc = (np.arange(n_tubes) - 0.5 * (n_tubes - 1)) * s
    zc = -h - R
    thc = math.atan2(0.5 * b, 0.5 * s) if b > 0 else 0.0   # junction half-angle
    A, B, sa, sb, typ = [], [], [], [], []   # typ 0 arc, 1 chord
    loops = [[] for _ in range(n_tubes)]       # (panel index, sigma, kappa)
    idx = 0

    def add(a_, b_, da, db, t):
        nonlocal idx
        A.append(a_); B.append(b_); sa.append(np.full(len(a_), da)); sb.append(np.full(len(a_), db))
        typ.append(np.full(len(a_), t))
        out = np.arange(idx, idx + len(a_))
        idx += len(a_)
        return out

    # chords first (between tube i and i+1), bottom -> top so n = +y (A = left tube, B = right)
    chord = {}
    if b > 0:
        for i in range(n_tubes - 1):
            ym = 0.5 * (yc[i] + yc[i + 1])
            n = max(2, int(math.ceil(b / h_e)))
            zz = np.linspace(zc - 0.5 * b, zc + 0.5 * b, n + 1)
            a_ = np.stack([np.full(n, ym), zz[:-1]], 1)
            b_ = np.stack([np.full(n, ym), zz[1:]], 1)
            chord[i] = add(a_, b_, i + 1, i + 2, 1)
    for i in range(n_tubes):
        hasL, hasR = (i > 0 and b > 0), (i < n_tubes - 1 and b > 0)
        t0 = thc if hasR else 0.0
        t1 = math.pi - thc if hasL else math.pi
        t2 = math.pi + thc if hasL else math.pi
        t3 = 2 * math.pi - thc if hasR else 2 * math.pi
        top = add(*_arc(yc[i], zc, R, t0, t1, h_e), i + 1, 0, 0)
        loops[i] += [(j, 1.0, 1.0 / R) for j in top]
        if hasL:
            loops[i] += [(j, -1.0, 0.0) for j in chord[i - 1][::-1]]   # downward
        bot = add(*_arc(yc[i], zc, R, t2, t3, h_e), i + 1, 0, 0)
        loops[i] += [(j, 1.0, 1.0 / R) for j in bot]
        if hasR:
            loops[i] += [(j, 1.0, 0.0) for j in chord[i]]               # upward
    P = bem2d.Panels(np.concatenate(A), np.concatenate(B))
    typ = np.concatenate(typ)
    nskin = np.where(typ == 1, 2.0, 1.0)
    return dict(P=P, side_a=np.concatenate(sa), side_b=np.concatenate(sb), kind=np.zeros(len(P), int),
                typ=typ, nskin=nskin, loops=loops, ncell=n_tubes, R=R, D=D, m_s=m_s, dp=dp, Et=Et,
                N=dp * R, Nx=dp * R / 2.0, yc=yc, zc=zc, b=b, gap=gap, V=V, gas=gas,
                width=float(yc[-1] - yc[0] + D))


def tension_ties(geo, w):
    """List of (i, j, val) for the hoop geometric stiffness + hoop stretching."""
    P = geo["P"]
    N, Et = geo["N"], geo["Et"]
    acc = {}
    for loop in geo["loops"]:
        idx = np.array([l[0] for l in loop])
        sig = np.array([l[1] for l in loop])
        kap = np.array([l[2] for l in loop])
        L = P.L[idx]
        n = idx.size
        s = np.cumsum(L) - 0.5 * L          # panel-centre arc length
        C = L.sum()
        Kop = np.zeros((n, n))
        for j in range(n):
            jm, jp = (j - 1) % n, (j + 1) % n
            hm = (s[j] - s[jm]) % C
            hp = (s[jp] - s[j]) % C
            Kop[j, jp] += N * 2.0 / (hp * (hp + hm))
            Kop[j, jm] += N * 2.0 / (hm * (hp + hm))
            Kop[j, j] += -N * 2.0 / (hp * hm) + N * kap[j] ** 2
        if Et:
            Kop -= np.outer(kap, kap * L) * (Et / C)
        # f_n,j = sig_j * sum_i Kop_ji sig_i w_n,i,  w = v / (j w)
        M = (sig[:, None] * Kop * sig[None, :]) / (1j * w)
        jj, ii = np.nonzero(Kop)
        for a_, c_ in zip(jj, ii):
            key = (idx[a_], idx[c_])
            acc[key] = acc.get(key, 0.0) + M[a_, c_]
    return [(i, j, v) for (i, j), v in acc.items()]


def raft_problem(geo, f, k, kappa, rho=RHO):
    w = 2 * np.pi * f
    P = geo["P"]
    V = geo["V"]
    Z = geo["nskin"] * (1j * w * geo["m_s"] - 1j * geo["Nx"] * w / V ** 2)
    doms = [bem2d.Domain(bem2d.Helmholtz(k), rho, exterior=True)]
    doms += [bem2d.Domain(bem2d.ModHelmholtz(kappa), geo["gas"].rho) for _ in range(geo["ncell"])]
    ties = tension_ties(geo, w) if (geo["N"] or geo["Et"]) else []
    return bem2d.Problem(w, P, doms, geo["side_a"], geo["side_b"], geo["kind"], Z.astype(complex), ties=ties)


def wavenumbers(f, V, gas=H2):
    w = 2 * np.pi * f
    k = w * math.sqrt(1 / C0 ** 2 - 1 / V ** 2)
    s2 = 1 / V ** 2 - 1 / gas.c ** 2
    assert s2 > 0, "gas must be evanescent (c_g > V)"
    return w, k, w * math.sqrt(s2)


def element_size(f, V, D, gas=H2, npl=10.0, kappa_res=1.0):
    w, k, kap = wavenumbers(f, V, gas)
    return min(2 * np.pi / k / npl, kappa_res / kap, math.pi * D / 48.0)


def incident(ac, f, k):
    Sv, SL = cp.source_strengths(ac, np.array([f]))
    a = np.array([Sv[0], SL[0]])
    return lambda Y, Zc: cp.multipole_field(a, k, Y, Zc, grad=True)


def solve_raft(ac, f, X, D, W, h, **kw):
    V = ac.V
    w, k, kap = wavenumbers(f, V)
    geo = build_raft(D, W, h, h_e=element_size(f, V, D), V=V, **kw)
    prob = raft_problem(geo, f, k, kap)
    inc = incident(ac, f, k)
    P = geo["P"]
    p, py, pz = inc(P.c[:, 0], P.c[:, 1])
    sol = bem2d.solve(prob, p, py * P.n[:, 0] + pz * P.n[:, 1])
    X = np.asarray(X, float)
    p0 = inc(X[:, 0], X[:, 1])[0]
    return sol.field(X, inc_p=p0), p0, geo


def solve_slab(ac, f, X, W, h, T, **kw):
    cfg = cp.CrossPlaneConfig(ac, h=h, W=W, T=T, gas=H2, incident="linear", atm=ATM, **kw)
    sol = cp.bem_solve(cfg, f)
    X = np.asarray(X, float)
    return cp.bem_field(sol, X), sol.inc(X[:, 0], X[:, 1])[0], cfg


# ----------------------------------------------------------------------------- self test
def _series_mode_Z(k, ki, rho_e, rho_i, Zn, a, xs, X, nmax=60):
    """Line source at xs scattered by a K0-gas cylinder with membrane impedance Zn(n)."""
    rs = np.hypot(*xs)
    ths = np.arctan2(xs[1], xs[0])
    r = np.hypot(X[:, 0], X[:, 1])
    th = np.arctan2(X[:, 1], X[:, 0])
    w = None
    p = sp.hankel2(0, k * np.hypot(X[:, 0] - xs[0], X[:, 1] - xs[1]))
    for n in range(-nmax, nmax + 1):
        al = sp.hankel2(n, k * rs)
        J, Jp = sp.jv(n, k * a), sp.jvp(n, k * a)
        H, Hp = sp.hankel2(n, k * a), sp.h2vp(n, k * a)
        Zi, Zip = sp.iv(n, ki * a), sp.ivp(n, ki * a)
        Z, w = Zn(abs(n))
        M = np.array([[k * Hp / rho_e, -ki * Zip / rho_i], [-H + Z / (1j * w * rho_e) * k * Hp, Zi]], complex)
        rhs = np.array([-al * k * Jp / rho_e, al * J - Z / (1j * w * rho_e) * al * k * Jp])
        bb, _ = np.linalg.solve(M, rhs)
        p = p + bb * sp.hankel2(n, k * r) * np.exp(1j * n * (th - ths))
    return p


def self_test(verbose=True):
    """Single tensioned tube (gap model, one tube) vs mode series; returns max rel. error."""
    f, V = 8.0, 590.0
    w, k, kap = wavenumbers(f, V)
    D, dp, Et, m_s = 30.0, 400.0, 3e4, 0.6
    geo = build_raft(D, D, 0.0, gap=1.0, dp=dp, Et=Et, m_s=m_s, h_e=0.4, V=V, n_tubes=1)
    # move centre to the origin for the series
    P0 = geo["P"]
    zc = geo["zc"]
    geo["P"] = bem2d.Panels(P0.a - [0, zc], P0.b - [0, zc])
    prob = raft_problem(geo, f, k, kap)
    xs = np.array([0.0, 120.0])
    P = geo["P"]
    d = P.c - xs
    r = np.hypot(d[:, 0], d[:, 1])
    sol = bem2d.solve(prob, sp.hankel2(0, k * r), -k * sp.hankel2(1, k * r) * (d * P.n).sum(1) / r)
    X = np.array([[0, -60.0], [40, -40], [-50, 10], [0, 80], [100, -20], [0, -400.0]])
    dX = X - xs
    pb = sol.field(X, inc_p=sp.hankel2(0, k * np.hypot(dX[:, 0], dX[:, 1])))
    R, N, Nx = D / 2, dp * D / 2, dp * D / 4

    def Zn(n):
        stiff = N * (n * n - 1) / R ** 2 + (Et / R ** 2 if n == 0 else 0.0)
        return 1j * w * m_s - 1j * Nx * w / V ** 2 + stiff / (1j * w), w

    ps = _series_mode_Z(k, kap, RHO, H2.rho, Zn, R, xs, X)
    # same without tension, to show the tension matters
    def Z0(n):
        return 1j * w * m_s, w
    p0 = _series_mode_Z(k, kap, RHO, H2.rho, Z0, R, xs, X)
    err = np.max(np.abs(pb - ps) / np.abs(ps))
    if verbose:
        print("self-test: max rel err vs series %.2e; tension effect %.2e" %
              (err, np.max(np.abs(ps - p0) / np.abs(ps))))
    return err


if __name__ == "__main__":
    self_test()
