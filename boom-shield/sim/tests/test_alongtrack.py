"""Tests for the along-track-varying (3-D) shield model (boomsim.alongtrack)."""
import math
import os
import sys

import numpy as np
import pytest
from scipy import special as sp

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))

from boomsim import alongtrack as AT  # noqa: E402
from boomsim import crossplane as cp  # noqa: E402
from boomsim import shield as SH  # noqa: E402
from boomsim.atmosphere import std_atm  # noqa: E402
from boomsim.propagate import RayPath, uniform_atmosphere  # noqa: E402
from boomsim.proxy import scored_metrics  # noqa: E402
from boomsim.source import CONCORDE  # noqa: E402


def _cfg(**kw):
    import dataclasses
    ac = dataclasses.replace(CONCORDE, mach=2.0, altitude=20250.0)
    a = std_atm(20000.0)
    g = SH.gas("hydrogen", float(a.T), float(a.p))
    kw.setdefault("W", 1000.0)
    return cp.CrossPlaneConfig(ac, h=250.0, T=30.0, gas=g, incident="linear", f_max=100.0,
                               df=0.5, R_m=3000.0, phis=[0.0], workers=1, **kw)


# ----------------------------------------------------------------------------- masks
@pytest.mark.parametrize("shape", ["sawtooth", "rc_teeth", "random", "linear_sawtooth"])
def test_serration_parseval_and_mean(shape):
    s = AT.Serration(20.0, 150.0, shape, seed=3)
    u = np.linspace(0.01, 0.99, 23)
    n = np.arange(-4000, 4001)
    C = s.coefficients(u, n)
    tau = s.tau(u)
    assert np.allclose(C[n == 0][0].real, tau, atol=1e-12)
    # sum_n |c_n|^2 = (1/L) int chi^2 = tau (tail ~ 1/n: truncation error ~ 1/(pi^2 N))
    assert np.allclose(np.sum(np.abs(C) ** 2, axis=0), tau, atol=2e-4 * s.n_cells)
    # one-sided and symmetric teeth have the same |c_n| (only the phase differs)
    if shape == "sawtooth":
        r = AT.Serration(20.0, 150.0, "rc_teeth").coefficients(u, n[:50])
        assert np.allclose(np.abs(C[:50]), np.abs(r), atol=1e-12)


def test_random_serration_keeps_cell_coverage():
    s = AT.Serration(10.0, 150.0, "random", seed=1)
    u = np.linspace(0, 1, 41)
    a, b = s.intervals(u)
    # covered length in each cell = tau * Lambda (no overlap, no root leak)
    for c in range(s.n_cells):
        lo, hi = c * s.Lam, (c + 1) * s.Lam
        ln = np.sum(np.clip(np.minimum(b, hi) - np.maximum(a, lo), 0, None), axis=1)
        assert np.allclose(ln, s.tau(u) * s.Lam, atol=1e-9)


def test_joint_coefficients():
    j = AT.Joint(250.0, 5.0)
    c = j.coefficients(np.arange(-20000, 20001))
    assert abs(c[20000] - 5.0 / 250.0) < 1e-14
    assert abs(np.sum(np.abs(c) ** 2) - 5.0 / 250.0) < 2e-5


# ----------------------------------------------------------------------------- PO aperture
def test_stack_difference_matches_po_solve():
    cfg = _cfg()
    f = 20.0
    a = cp.incident_coefficients(cfg, [f])[0]
    ap = cp.po_solve(cfg, f, a)
    D = AT.stack_difference(cfg, f, ap["y"], a=a)
    assert np.max(np.abs(D - ap["d"])) < 1e-9 * np.max(np.abs(ap["d"]))
    # an "air" fill with massless faces is transparent: D = 0
    import dataclasses
    c2 = dataclasses.replace(cfg, m_top=1e-6, m_bot=1e-6)
    air = SH.gas("air", float(cfg.atm.T), float(cfg.atm.p))
    D0 = AT.stack_difference(c2, f, ap["y"][::50], a=a, gas=air)
    assert np.max(np.abs(D0)) < 1e-3 * np.max(np.abs(ap["d"]))


# ----------------------------------------------------------------------------- rays
@pytest.mark.parametrize("sx,sy", [(1 / 590.0, 0.0), (1 / 590.0, 0.0012), (0.0, 0.0025),
                                   (1 / 2000.0, -0.002), (0.001, 0.0027)])
def test_ground_rays_match_raypath(sx, sy):
    z_src = 19985.0
    rays = AT.GroundRays(z_src, z_anchor=20250.0)
    dy, G, ok = rays.gain(sx, sy, 50.0)
    assert ok and G > 0
    V = 1.0 / sx if sx > 0 else 1e12
    ca = rays.c_s
    cmu = math.sqrt(1 - (ca / V) ** 2)
    phi = math.asin(sy * ca / cmu)
    rp = RayPath(z_src, V, phi=phi, allow_cutoff=True)
    assert abs(float(dy) - rp.y_ground) < 1e-4 * max(1.0, abs(rp.y_ground)) + 0.1
    # Jacobian dy_g/ds_y = I1 + s_y^2 I3 vs RayPath.dy_dphi
    szg = math.sqrt(1 / rays.c_g ** 2 - sx ** 2 - sy ** 2)
    I13 = np.interp(szg, rays.szg_tab, rays.I1_tab) + sy ** 2 * np.interp(szg, rays.szg_tab, rays.I3s_tab) / szg
    I13_rp = rp.dy_dphi(0.0) * ca / (cmu * math.cos(phi))
    assert abs(I13 / I13_rp - 1) < 1e-4


def test_ground_rays_cutoff():
    rays = AT.GroundRays(19985.0)
    _, G, ok = rays.gain(1 / 330.0, 0.0, 10.0)        # trace speed below c_ground: turns
    assert not ok and G == 0


def test_flux_mapping_uniform_matches_rayleigh():
    """In uniform air the flux mapping equals the far field of the exact 2-D Rayleigh
    integral of the order (incl. the oblique order wavenumber kappa_n)."""
    c = 300.0
    z_src = 3000.0
    rays = AT.GroundRays(z_src, atmosphere=uniform_atmosphere(0.0), dz=10.0, reflection=1.0)
    f = 40.0
    w = 2 * np.pi * f
    rays._abs_cache[float(f)] = np.zeros(rays.szg_tab.size)
    c = rays.c_g
    sx = 0.6 / c
    k0 = w / c
    kap = k0 * math.sqrt(1 - (sx * c) ** 2)
    y = np.arange(-40.0, 40.0, 0.1) + 0.05
    d = np.exp(-(y / 12.0) ** 2) * np.exp(-1j * 0.3 * kap * y)
    ky = np.linspace(-0.95, 0.95, 3001) * kap
    Dt = (np.exp(1j * np.outer(ky, y)) @ d) * 0.1
    dyo, G, ok = rays.gain(sx, ky / w, f)
    yg = np.array([-1500.0, -800.0, -300.0, 200.0, 900.0])
    E_flux = np.interp(yg, dyo[ok], (np.abs(Dt) ** 2 * G)[ok])
    X = np.stack([yg, np.full(yg.size, -z_src)], 1)
    p = AT.rayleigh_2d(d, y, 0.1, 0.0, kap, X)
    assert np.allclose(E_flux, np.abs(p) ** 2, rtol=0.02)


def test_orders_reproduce_3d_rayleigh():
    """Sum of Floquet orders (each a 2-D Rayleigh integral with its own kappa_n) equals
    the brute-force 3-D Rayleigh integral of a periodic aperture (lossy medium for
    convergence) -- checks the sign convention k_x,n = k_x0 + 2 pi n / L."""
    k = 1.0 * (1 - 0.03j)
    kx0 = 0.5
    Lam = 15.0
    s = AT.Serration(Lam, 1.0, "sawtooth")                 # use the interval machinery
    a, b = np.array([[0.0]]), np.array([[6.0]])             # chi = 1 on [0, 6) per period
    gy = lambda yy: np.exp(-(yy / 5.0) ** 2)                # noqa: E731
    P = np.array([3.0, 2.0, -20.0])
    # 3-D: p = (1/2 pi) int d |z|/R (j k + 1/R) e^{-jkR}/R dx dy
    dx = 0.25
    x = np.arange(-300.0, 300.0, dx) + dx / 2
    yv = np.arange(-15.0, 15.0, 0.25) + 0.125
    chi = ((np.mod(x, Lam) >= 0) & (np.mod(x, Lam) < 6.0)).astype(float)
    Xg, Yg = np.meshgrid(x, yv, indexing="ij")
    d = (chi[:, None] * np.exp(-1j * kx0 * Xg)) * gy(Yg)
    Rr = np.sqrt((P[0] - Xg) ** 2 + (P[1] - Yg) ** 2 + P[2] ** 2)
    p3 = np.sum(d * abs(P[2]) / Rr * (1j * k + 1 / Rr) * np.exp(-1j * k * Rr) / Rr) * dx * 0.25 / (2 * np.pi)
    # orders
    n = np.arange(-60, 61)
    cn = AT.interval_coefficients(a, b, Lam, n)[:, 0]
    y2 = np.arange(-15.0, 15.0, 0.05) + 0.025
    p2 = 0.0
    for ni, ci in zip(n, cn):
        kxn = kx0 + 2 * np.pi * ni / Lam
        kap = np.sqrt(k * k - kxn * kxn + 0j)
        if kap.imag > 0:
            kap = -kap
        pn = AT.rayleigh_2d(gy(y2) + 0j, y2, 0.05, 0.0, kap, np.array([[P[1], P[2]]]))[0]
        p2 += ci * np.exp(-1j * kxn * P[0]) * pn
    assert abs(p2 - p3) < 0.02 * abs(p3)


# ----------------------------------------------------------------------------- loudness
def _nwave(fs=20000.0):
    t = np.arange(-0.2, 0.6, 1 / fs)
    p = np.where((t >= 0) & (t < 0.3), 50.0 * (1 - 2 * t / 0.3), 0.0)
    return t, p


def test_combined_pl_equals_pldb_without_scatter():
    t, p = _nwave()
    m = scored_metrics(t, p)
    raw, prx = AT.combined_pl(AT.signature_bands(t, p), np.array([1.0, 2.0]), np.zeros(2))
    assert abs(raw - m["pldb"]) < 1e-9
    assert abs(prx - m["pldb_proxy"]) < 1e-9


def test_band_energies_esd_matches_signal():
    t, p = _nwave()
    eb = AT.band_energies_signal(t, p, taper_s=0.0)
    f = np.concatenate([np.arange(0.25, 20, 0.25), np.arange(20, 100, 0.5), np.arange(100, 400, 1.0)])
    dt = t[1] - t[0]
    P = np.exp(-2j * np.pi * np.outer(f, t)) @ p * dt
    ee = AT.band_energies_esd(f, 2 * np.abs(P) ** 2)
    fc = 10.0 ** (np.arange(42) / 10.0)
    m = (fc >= 3.0) & (fc <= 300.0)
    assert np.allclose(10 * np.log10(ee[m] / eb[m]), 0.0, atol=0.3)


def test_fringe_orders_budget_consistent():
    cfg = _cfg()
    rays = AT.GroundRays(cfg.aircraft.altitude + cfg.z_mid, z_anchor=cfg.aircraft.altitude)
    s = AT.Serration(30.0, 150.0, "sawtooth")
    yg = np.arange(-60000.0, 60001.0, 500.0)
    o = AT.fringe_orders(cfg, s, 40.0, rays, yg)
    assert 0.11 < o["ap_scat"] / o["ap_D"] < 0.14          # ~ int tau(1-tau)/w_f = 1/8
    assert 0 <= o["ap_ground"] <= o["ap_prop"] <= o["ap_scat"] * 1.01
    assert np.all(o["E"] >= 0) and o["E"].max() > 0
    # orders of a short period at low frequency are all evanescent
    s2 = AT.Serration(10.0, 150.0, "sawtooth")
    o2 = AT.fringe_orders(cfg, s2, 15.0, rays, yg)
    assert o2["ap_ground"] == 0.0 and o2["ap_prop"] < 1e-4 * o2["ap_scat"]
