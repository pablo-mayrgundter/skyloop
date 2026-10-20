import os
import sys

import numpy as np
import pytest

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))

from boomsim.atmosphere import saturation_ratio, std_atm  # noqa: E402
from boomsim.propagate import (RayPath, burgers_step, ground_signature,  # noqa: E402
                               iso_absorption, propagate, propagation_exponent,
                               propagation_filter, relaxation_params,
                               uniform_atmosphere)
from boomsim.source import CONCORDE, Aircraft  # noqa: E402

A0 = std_atm(0.0)
C0 = float(A0.c)


def _uniform_ray(altitude=20100.0, mach=2.0):
    return RayPath(altitude, mach * C0, atmosphere=uniform_atmosphere(0.0),
                   humidity=lambda z: np.full_like(np.asarray(z, float), 1.0))


# ------------------------------------------------------------------ absorption
def test_iso_9613_table_values():
    # ISO 9613-2 Table 2 (from ISO 9613-1): 20 C, 70 % RH, 101.325 kPa, dB/km
    T, p = 293.15, 101325.0
    h = 70.0 * saturation_ratio(T) / (p / 101325.0)
    ref = {63: 0.1, 125: 0.3, 250: 1.1, 500: 2.8, 1000: 5.0, 2000: 9.0, 4000: 22.9, 8000: 76.6}
    f = np.array(list(ref))
    a = 8.686 * 1000 * iso_absorption(f * 1.0, T, p, h)
    for fi, ai in zip(f, a):
        assert ai == pytest.approx(ref[fi], abs=0.06, rel=0.02), fi


def test_filter_magnitude_matches_iso():
    f = np.logspace(0, 4.5, 200)
    for z in (0.0, 5000.0, 17000.0):
        a = std_atm(z)
        h = 1.0 if z == 0 else 0.05
        H = propagation_filter(f, a.T, a.p, a.c, h, 1000.0)
        alpha = iso_absorption(f, a.T, a.p, h)
        assert np.allclose(-np.log(np.abs(H)) / 1000.0, alpha, rtol=1e-10)


def test_relaxation_dispersion_is_causal_and_speeds_up_high_f():
    a = std_atm(0.0)
    c, h, ds = float(a.c), 1.0, 2000.0
    r = relaxation_params(a.T, a.p, h)
    # phase speed from the relaxation terms: k(w) = w/c - Im(ln H)/ds
    f = np.logspace(0, 7, 600)
    lnH = propagation_exponent(f, a.T, a.p, c, h, ds, thermoviscous=False)
    assert np.all(lnH.imag >= 0)                   # advance relative to c (rfft convention)
    k = 2 * np.pi * f / c - lnH.imag / ds
    cp = 2 * np.pi * f / k
    assert np.all(np.diff(cp) >= -1e-9)            # monotone increase with frequency
    assert cp[0] == pytest.approx(c, rel=1e-9)
    # high-frequency limit: (c_inf - c)/c = sum eps_v/2 with eps_v = A_v c/pi
    eps = float((r["A_O"] + r["A_N"]) * c / np.pi)
    assert (cp[-1] - c) / c == pytest.approx(eps / 2, rel=2e-3)
    # each process: dispersion step (c_inf - c_0)/c ~ 1e-6..1e-3 (typical air values)
    assert 1e-6 < eps / 2 < 1e-3
    # causality: impulse response has (essentially) nothing before the
    # infinite-frequency arrival tau = -ds*eps/(2c); the zero-phase version does.
    fs, n = 200e3, 2 ** 16
    fr = np.fft.rfftfreq(n, 1 / fs)
    adv = ds * eps / (2 * c)
    for disp, expect_causal in ((True, True), (False, False)):
        Hc = propagation_filter(fr, a.T, a.p, c, h, ds, thermoviscous=False, dispersion=disp)
        g = np.fft.irfft(Hc, n=n)
        tau = np.where(np.arange(n) < n // 2, np.arange(n), np.arange(n) - n) / fs
        e_post = np.sum(g[tau > 1.2 * adv] ** 2)    # dispersive tail (after the front)
        e_pre = np.sum(g[tau < -1.2 * adv] ** 2)    # precursor (must vanish if causal)
        if expect_causal:
            assert e_pre / e_post < 1e-3
        else:
            assert e_pre / e_post > 0.5             # zero-phase: symmetric


# ------------------------------------------------------------------ geometry
def test_raypath_uniform_is_cylindrical():
    ray = _uniform_ray()
    za = ray.altitude
    nz = np.sqrt(1 - 1 / 2.0 ** 2)
    z = np.array([za - 100.0, za - 1000.0, 1000.0])
    assert np.allclose(ray.nz(z), nz)
    assert np.allclose(ray.s(z), (za - z) / nz, rtol=1e-6)
    A = ray.area(z)
    assert np.allclose(A / A[0], (za - z) / (za - z[0]), rtol=1e-6)


def test_raypath_stratified_sanity():
    ray = RayPath(17000.0, 590.0)
    z = np.array([16000.0, 8000.0, 0.0])
    # n_z = sqrt(1 - c^2/V^2): ray flattens where c is larger (lower altitude)
    assert np.all(np.diff(ray.nz(z)) < 0)
    assert np.all(np.diff(ray.area(z)) > 0)
    with pytest.raises(ValueError):
        RayPath(17000.0, 320.0)           # V < c(ground): cut-off


# ------------------------------------------------------------------ propagation
def test_linear_lossless_uniform_spreading():
    ray = _uniform_ray()
    fs = 20e3
    t = np.arange(-0.05, 0.05, 1 / fs)
    p = 100.0 * np.exp(-(t / 2e-3) ** 2)
    za = ray.altitude
    z1, z2 = za - 1000.0, za - 10000.0
    t2, p2 = propagate(t, p, ray, z1, z2, linear=True, lossless=True)
    assert p2.max() / p.max() == pytest.approx(np.sqrt(1000.0 / 10000.0), rel=1e-6)
    # and stratified: matches the Blokhintzev invariant
    ray_s = RayPath(17000.0, 590.0)
    t3, p3 = propagate(t, p, ray_s, 16000.0, 0.0, linear=True, lossless=True)
    assert p3.max() / p.max() == pytest.approx(float(ray_s.amplitude_ratio(16000.0, 0.0)), rel=1e-6)


def test_linear_absorption_of_tone():
    ray = _uniform_ray()
    fs, f0 = 40e3, 2000.0
    t = np.arange(4000) / fs                       # exactly 200 periods
    p = np.sin(2 * np.pi * f0 * t)
    za = ray.altitude
    z1, z2 = za - 5000.0, za - 5000.0 - 2000.0 * ray.nz(za)    # 2 km of ray
    t2, p2 = propagate(t, p, ray, z1, z2, linear=True, auto_extend=False)
    spread = float(ray.amplitude_ratio(z1, z2))
    alpha = iso_absorption(f0, A0.T, A0.p, 1.0)
    amp = np.sqrt(2 * np.mean(p2 ** 2))
    assert amp == pytest.approx(spread * np.exp(-alpha * 2000.0), rel=1e-6)


def test_burgers_step_plane_nwave():
    # plane N-wave, weak-shock theory: P = P0/sqrt(1+X P0/T0), T = T0 sqrt(1+X P0/T0)
    fs, T0, P0, k = 50e3, 0.02, 100.0, 1.0
    dt = 1 / fs
    t = np.arange(-0.15, 0.15, dt)
    p = np.where(np.abs(t) < T0, -P0 * t / T0, 0.0)
    X, Xend = 0.0, 3 * T0 / P0
    while X < Xend - 1e-15:
        ds = min(0.8 * dt / (k * np.abs(p).max()), Xend - X)
        p = burgers_step(p, k, dt, ds)
        X += ds
    fac = np.sqrt(1 + X * P0 / T0)
    assert p.max() == pytest.approx(P0 / fac, rel=2e-3)
    assert p.min() == pytest.approx(-P0 / fac, rel=2e-3)
    Ipos = np.sum(p[p > 0]) * dt
    assert Ipos == pytest.approx(P0 * T0 / 2, rel=2e-3)        # conservation
    assert 2 * Ipos / p.max() == pytest.approx(T0 * fac, rel=3e-3)


def test_nonlinear_lossless_uniform_nwave_cylindrical():
    """N-wave with cylindrical spreading in uniform air vs weak-shock theory;
    asymptotically peak ~ r^-3/4 and duration ~ r^1/4."""
    ray = _uniform_ray()
    za = ray.altitude
    nz = float(ray.nz(za))
    fs = 50e3
    T0, P0 = 0.005, 2000.0
    t = np.arange(-0.02, 0.06, 1 / fs)
    p = np.where(np.abs(t) < T0, -P0 * t / T0, 0.0)
    z1 = za - 100.0
    s1 = float(ray.s(z1))
    k = 1.2 / (float(A0.rho) * C0 ** 3)
    res = {}
    zc, tc, pc = z1, t, p
    for zt in (za - 5000.0, za - 10000.0, za - 20000.0):
        tc, pc = propagate(tc, pc, ray, zc, zt, lossless=True)
        zc = zt
        s = float(ray.s(zt))
        X = 2 * k * np.sqrt(s1) * (np.sqrt(s) - np.sqrt(s1))
        fac = np.sqrt(1 + X * P0 / T0)
        P_th = P0 * np.sqrt(s1 / s) / fac
        T_th = T0 * fac
        dt = tc[1] - tc[0]
        P_num = pc.max()
        T_num = 2 * np.sum(pc[pc > 0]) * dt / P_num
        assert P_num == pytest.approx(P_th, rel=0.02)
        assert T_num == pytest.approx(T_th, rel=0.02)
        res[s] = (P_num, T_num)
    s_a, s_b = sorted(res)[-2:]
    slope_p = np.log(res[s_b][0] / res[s_a][0]) / np.log(s_b / s_a)
    slope_t = np.log(res[s_b][1] / res[s_a][1]) / np.log(s_b / s_a)
    assert -0.82 < slope_p < -0.70
    assert 0.18 < slope_t < 0.30


def test_insert_hook():
    ac = Aircraft("t", length=30.0, mach=2.0, altitude=15000.0, weight=3e5, volume=60.0)
    opts = dict(fs=10e3, linear=True, lossless=True)
    t0, p0 = ground_signature(ac, **opts)
    calls = []

    def ident(t, p):
        calls.append(p.size)
        return p

    t1, p1 = ground_signature(ac, insert=(8000.0, ident), **opts)
    assert len(calls) == 1
    n = min(p0.size, p1.size)
    assert np.allclose(p0[:n], p1[:n], atol=1e-9 * np.abs(p0).max())
    t2, p2 = ground_signature(ac, insert=(8000.0, lambda t, p: (t, 0.5 * p)), **opts)
    assert np.abs(p2).max() == pytest.approx(0.5 * np.abs(p0).max(), rel=1e-9)
    t3, p3 = ground_signature(ac, reflection=1.0, **opts)
    assert np.abs(p3).max() == pytest.approx(np.abs(p0).max() / 1.9, rel=1e-12)


def test_concorde_ground_signature():
    t, p = ground_signature(CONCORDE, fs=25e3)
    pk = p.max()
    assert 80.0 < pk < 110.0                      # ~1.94 psf measured (NASA FS-016)
    dt = t[1] - t[0]
    i_pk = int(np.argmax(p))
    i10 = np.flatnonzero(p[: i_pk + 1] < 0.1 * pk)[-1]
    i90 = np.flatnonzero(p[: i_pk + 1] < 0.9 * pk)[-1]
    rise = (i90 - i10) * dt
    assert 2e-4 < rise < 5e-3
    i_mn = int(np.argmin(p))
    i_end = i_mn + np.flatnonzero(p[i_mn:] > 0.1 * p.min())[0]
    dur = (i_end - i10) * dt
    assert 0.2 < dur < 0.35
    assert -1.1 * pk < p.min() < -0.6 * pk
