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


def test_below_opts_linear_leg():
    """below_opts applies only below the insert: linear below == propagate legs by hand."""
    from boomsim.source import near_field
    ac = CONCORDE
    opts = dict(fs=10e3, r0=90.0)
    z_ins = 10000.0
    t1, p1 = ground_signature(ac, insert=(z_ins, lambda t, p: (t, p)), below_opts={"linear": True},
                              reflection=1.0, **opts)
    t0, p0 = near_field(ac, 90.0, 10e3)
    ray = RayPath(ac.altitude, ac.V)
    t0, p0 = propagate(t0, p0, ray, ac.altitude - 90.0, z_ins)
    t0, p0 = propagate(t0, p0, ray, z_ins, 0.0, linear=True)
    np.testing.assert_allclose(p1, p0, atol=1e-9 * np.max(np.abs(p0)))
    # linear below the insert leaves a much less steep front than the nonlinear run
    t2, p2 = ground_signature(ac, insert=(z_ins, lambda t, p: (t, p)), reflection=1.0, **opts)
    assert np.max(np.diff(p2)) > 1.5 * np.max(np.diff(p1))
    with pytest.raises(ValueError):
        ground_signature(ac, below_opts={"linear": True}, **opts)


# ------------------------------------------------------------------ off-track rays
from boomsim.propagate import (circle_start_factors, ground_carpet,  # noqa: E402
                               lateral_cutoff, propagate_from_circle)
from boomsim.proxy import RISE_SIGMA_S, rise_time_proxy, scored_metrics  # noqa: E402
from boomsim.source import near_field  # noqa: E402


def test_raypath_phi0_is_legacy_formula():
    ray = RayPath(17000.0, 590.0, phi=0.0)
    ray_d = RayPath(17000.0, 590.0)
    for name in ("nz_tab", "s_tab", "lat_tab", "rhoc_tab", "h_tab"):
        assert np.array_equal(getattr(ray, name), getattr(ray_d, name))
    # A = c s_z int_z^{z_a} dz'/s_z, s_z = sqrt(1/c^2 - 1/V^2) (trapezoid)
    from scipy.integrate import cumulative_trapezoid
    z, c = ray.z_tab, ray.c_tab
    sz = np.sqrt(1 / c ** 2 - 1 / 590.0 ** 2)
    I1 = -cumulative_trapezoid(1 / sz[::-1], z[::-1], initial=0.0)[::-1]
    zz = np.linspace(0, 16990, 101)
    assert np.array_equal(ray.area(zz), np.interp(zz, z, c * sz) * np.interp(zz, z, I1))
    # general (off-track) cell rule is continuous as phi -> 0
    ray_e = RayPath(17000.0, 590.0, phi=1e-7)
    assert np.allclose(ray_e.area(zz), ray.area(zz), rtol=1e-6)
    assert np.allclose(ray_e.s(zz), ray.s(zz), rtol=1e-6)
    assert ray.y_ground == 0.0 and not ray.cutoff


@pytest.mark.parametrize("phi", [0.0, 0.3, 0.8, -0.6])
def test_raypath_uniform_offtrack_straight(phi):
    za, V = 20000.0, 2.0 * C0
    ray = RayPath(za, V, atmosphere=uniform_atmosphere(0.0), phi=phi,
                  humidity=lambda z: np.full_like(np.asarray(z, float), 1.0))
    cmu = np.sqrt(1 - (C0 / V) ** 2)
    tmu = (C0 / V) / cmu
    z = np.array([za - 50.0, za - 1000.0, 7000.0, 0.0])
    dz = za - z
    r = dz / np.cos(phi)                           # distance from the flight axis
    assert np.allclose(ray.y(z), dz * np.tan(phi), rtol=1e-9, atol=1e-9)
    assert np.allclose(ray.x(z), r * tmu, rtol=1e-9)
    assert np.allclose(ray.s(z), r / cmu, rtol=1e-9)
    assert np.allclose(ray.travel_time(z), r / cmu / C0, rtol=1e-9)
    assert np.allclose(ray.nz(z), cmu * np.cos(phi), rtol=1e-12)
    assert np.allclose(ray.area(z), C0 * r, rtol=1e-6)          # cylindrical, A -> c r
    assert ray.y_ground == pytest.approx(za * np.tan(phi), rel=1e-9, abs=1e-9)
    arr = ray.arrival()
    assert arr["delay"] == pytest.approx(arr["T"] - arr["x"] / V)
    # straight-line geometry: emission-to-arrival distance = c T
    assert np.hypot(np.hypot(arr["x"], arr["y"]), za) == pytest.approx(C0 * arr["T"], rel=1e-9)
    # linear lossless amplitude ~ r^-1/2 along the off-track ray
    fs = 20e3
    t = np.arange(-0.05, 0.05, 1 / fs)
    p = 100.0 * np.exp(-(t / 2e-3) ** 2)
    _, p2 = propagate(t, p, ray, z[1], z[3], linear=True, lossless=True)
    assert p2.max() / p.max() == pytest.approx(np.sqrt(r[1] / r[3]), rel=1e-6)


@pytest.mark.parametrize("phi", [0.2, 0.6, 0.9])
def test_offtrack_tube_area_matches_neighbouring_rays(phi):
    """Analytic lateral Jacobian = finite difference of traced neighbours (std atm)."""
    za, V = CONCORDE.altitude, CONCORDE.V
    d = 1e-5
    ray = RayPath(za, V, phi=phi)
    rp, rm = RayPath(za, V, phi=phi + d), RayPath(za, V, phi=phi - d)
    z = np.array([za - 500.0, 11000.0, 3000.0, 0.0])
    fd = (rp.y(z) - rm.y(z)) / (2 * d)
    assert np.allclose(ray.dy_dphi(z), fd, rtol=1e-6)
    # tube area = |n_z| dY/dphi (dX/dx0 = 1), normalised by cos(mu0)/c_a
    assert np.allclose(ray.area(z), ray.nz(z) * ray.dy_dphi(z) * ray.c_a / ray.cos_mu0,
                       rtol=1e-12)
    # near the source A -> c_a r
    dz = 20.0
    assert float(ray.area(za - dz)) == pytest.approx(ray.c_a * dz / np.cos(phi), rel=2e-3)


def test_lateral_cutoff_concorde():
    za, V = CONCORDE.altitude, CONCORDE.V
    phi_c, y_c = lateral_cutoff(za, V)
    # analytic: sin(phi_c) = c_a sqrt(1/c_0^2 - 1/V^2)/cos(mu0) (c is largest at the ground)
    c_a, c_0 = float(std_atm(za).c), float(std_atm(0.0).c)
    assert np.sin(phi_c) == pytest.approx(c_a * np.sqrt(1 / c_0 ** 2 - 1 / V ** 2)
                                          / np.sqrt(1 - (c_a / V) ** 2), rel=1e-9)
    assert 35e3 < y_c < 50e3                     # ~1 mile per 1000 ft rule: ~44 km
    below = RayPath(za, V, phi=phi_c - 1e-3)
    above = RayPath(za, V, phi=phi_c + 1e-3)
    assert not below.cutoff and below.y_ground < y_c
    assert above.cutoff and 0.0 < above.z_turn < 2000.0 and np.isnan(above.y_ground)
    with pytest.raises(ValueError):
        RayPath(za, V, phi=phi_c + 1e-3, allow_cutoff=False)
    with pytest.raises(ValueError):
        propagate(np.arange(100) / 1e4, np.zeros(100), above, 5000.0, 0.0)
    with pytest.raises(ValueError):
        ground_signature(CONCORDE, phi=phi_c + 0.01, fs=5e3)
    # uniform air: no cut-off below 90 deg
    assert lateral_cutoff(17000.0, 2 * C0, atmosphere=uniform_atmosphere(0.0))[0] == pytest.approx(np.pi / 2)


def test_offtrack_ground_signature_concorde():
    t0, p0 = ground_signature(CONCORDE, fs=10e3)
    t1, p1 = ground_signature(CONCORDE, fs=10e3, phi=0.0)
    assert np.array_equal(p0, p1)
    t2, p2 = ground_signature(CONCORDE, fs=10e3, phi=np.radians(45.0))
    assert 0.4 * p0.max() < p2.max() < 0.8 * p0.max()


def _uniform_kw():
    return dict(atmosphere=uniform_atmosphere(0.0),
                humidity=lambda z: np.full_like(np.asarray(z, float), 1.0))


def test_circle_start_uniform_air_equals_direct():
    ac = Aircraft("t", length=30.0, mach=2.0, altitude=15000.0, weight=3e5, volume=60.0)
    kw = _uniform_kw()
    fs, r0, Rm = 10e3, 45.0, 1500.0
    for phi in (0.0, 0.5):
        # linear lossless: circle start == direct start (uniform air: factor exactly 1)
        tc, pc = near_field(ac, Rm, fs, phi=phi)
        info = {}
        tg, pg, yg = propagate_from_circle(ac, phi, tc, pc, Rm, info=info, linear=True,
                                           lossless=True, **kw)
        assert info["factor"] == pytest.approx(1.0, abs=1e-6)
        assert yg == pytest.approx(ac.altitude * np.tan(phi), rel=1e-9, abs=1e-9)
        td, pd = ground_signature(ac, r0=r0, fs=fs, phi=phi, linear=True, lossless=True, **kw)
        assert pg.max() == pytest.approx(pd.max(), rel=1e-4)
        # nonlinear: aging correction from r0 reproduces the direct run
        tg, pg, _ = propagate_from_circle(ac, phi, tc, pc, Rm, age_from=r0, **kw)
        td, pd = ground_signature(ac, r0=r0, fs=fs, phi=phi, **kw)
        n = min(pg.size, pd.size)
        assert np.allclose(pg[:n], pd[:n], atol=1e-6 * pd.max())
        # without it the ground boom is too strong (missing nonlinear age inside R_m)
        tg, pg, _ = propagate_from_circle(ac, phi, tc, pc, Rm, **kw)
        assert pg.max() > 1.05 * pd.max()


def test_circle_start_wkb_factor_isothermal():
    ac = Aircraft("t", length=30.0, mach=2.0, altitude=19000.0, weight=3e5, volume=60.0)
    H = 287.0531 * float(std_atm(15000.0).T) / 9.80665
    for phi, Rm in ((0.0, 2000.0), (0.6, 4000.0)):
        f = circle_start_factors(ac, phi, Rm)
        dz = Rm * np.cos(phi)
        assert f["z_start"] == pytest.approx(19000.0 - dz)
        # 11-20 km isothermal: c cancels, rho ~ exp(-z/H); straight rays: geometric = 1
        assert f["impedance"] == pytest.approx(np.exp(dz / (2 * H)), rel=2e-3)
        assert f["geometric"] == pytest.approx(1.0, abs=1e-3)
        f_ref = circle_start_factors(ac, phi, Rm, z_ref=18000.0, anchor="ref")
        assert f_ref["impedance"] == pytest.approx(np.exp((18000.0 - (19000.0 - dz)) / (2 * H)), rel=2e-3)


def test_ground_carpet_structure():
    ac = Aircraft("t", length=30.0, mach=2.0, altitude=15000.0, weight=3e5, volume=60.0)
    phi_c, y_c = lateral_cutoff(ac.altitude, ac.V)
    phis = [0.0, 0.5, phi_c - 0.02, phi_c + 0.05]
    opts = dict(fs=8e3, linear=True, lossless=True)
    res = ground_carpet(ac, phis, procs=2, **opts)
    assert list(res["cutoff"]) == [False, False, False, True]
    assert res["y_g"][0] == 0.0 and np.all(np.diff(res["y_g"][:3]) > 0) and np.isnan(res["y_g"][3])
    assert res["t"][3] is None and res["metrics"][3] is None
    assert np.all(np.diff(res["peak_pa"][:3]) < 0)          # weaker off-track
    c = res["carpet"]
    assert c["half_width_m"] == pytest.approx(y_c) and c["phi_cutoff"] == pytest.approx(phi_c)
    assert c["undertrack_pldb"] == res["pldb"][0] == c["max_pldb"]
    assert res["y_g"][2] < y_c
    # identical to the serial direct call
    t0, p0 = ground_signature(ac, phi=0.5, **opts)
    assert np.array_equal(res["p"][1], p0)
    m = scored_metrics(t0, p0)
    assert res["pldb_proxy"][1] == m["pldb_proxy"]
    # circle start (dict of signatures)
    Rm = 1000.0
    circ = {f: near_field(ac, Rm, 8e3, phi=f) for f in phis[:2]}
    rc = ground_carpet(ac, phis[:2], start="circle", circle=circ, R_m=Rm, procs=1, **opts)
    assert np.allclose(rc["peak_pa"], res["peak_pa"][:2], rtol=0.02)


def test_proxy_matches_gate0_definition():
    from scipy.ndimage import gaussian_filter1d
    fs = 50e3
    t, p = ground_signature(Aircraft("t", length=30.0, mach=2.0, altitude=15000.0,
                                     weight=3e5, volume=60.0), fs=fs, linear=True, lossless=True)
    ps = gaussian_filter1d(p, RISE_SIGMA_S * fs, mode="constant")
    assert np.allclose(rise_time_proxy(t, p), ps, rtol=1e-9, atol=1e-12 * np.abs(p).max())


def test_aged_near_field_reproduces_direct():
    from boomsim.propagate import aged_near_field
    ac = Aircraft("t", length=30.0, mach=2.0, altitude=15000.0, weight=3e5, volume=60.0)
    kw = _uniform_kw()
    fs, r0, Rm = 10e3, 45.0, 1500.0
    for phi in (0.0, 0.5):
        tc, pc = aged_near_field(ac, Rm, fs, phi=phi, r0=r0, **kw)
        tg, pg, _ = propagate_from_circle(ac, phi, tc, pc, Rm, **kw)
        td, pd = ground_signature(ac, r0=r0, fs=fs, phi=phi, **kw)
        n = min(pg.size, pd.size)
        # same physics; only the absorption-chunk boundaries differ (one more stop at R_m)
        assert np.allclose(pg[:n], pd[:n], atol=2e-3 * pd.max())
    t1, p1 = aged_near_field(ac, 30.0, fs, r0=r0)          # inside r0: linear near field
    assert np.array_equal(p1, near_field(ac, 30.0, fs)[1])
