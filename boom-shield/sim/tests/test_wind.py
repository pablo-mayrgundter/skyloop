"""Ray tracing and Burgers propagation in a horizontally stratified moving medium."""
import dataclasses
import os
import sys

import numpy as np
import pytest

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))

from boomsim import propagate as PR  # noqa: E402
from boomsim.source import CONCORDE  # noqa: E402

AC = dataclasses.replace(CONCORDE, mach=2.0, altitude=20250.0)
UATM = PR.uniform_atmosphere(0.0)
C0 = float(UATM(0.0).c)


def test_default_is_windless():
    r = PR.RayPath(AC.altitude, AC.V, phi=0.2, allow_cutoff=True)
    r2 = PR.RayPath(AC.altitude, AC.V, phi=0.2, allow_cutoff=True, wind=None)
    assert r.wind is None and r2.wind is None
    assert np.array_equal(r.area(np.r_[0.0, 5e3]), r2.area(np.r_[0.0, 5e3]))
    assert np.all(r.omega(np.r_[0.0, 1e4]) == 1.0)


@pytest.mark.parametrize("phi", [0.0, 0.3, 0.6, -0.5])
def test_calm_matches_windless(phi):
    r0 = PR.RayPath(AC.altitude, AC.V, phi=phi, allow_cutoff=True)
    r1 = PR.RayPath(AC.altitude, AC.V, phi=phi, allow_cutoff=True, wind=PR.calm())
    z = np.linspace(0.0, 19000.0, 9)
    assert np.allclose(r1.area(z), r0.area(z), rtol=1e-7)
    assert r1.y_ground == pytest.approx(r0.y_ground, abs=1e-6)
    assert r1.travel_time(0.0) == pytest.approx(r0.travel_time(0.0), rel=1e-12)
    assert r1.x(0.0) == pytest.approx(r0.x(0.0), rel=1e-8)
    assert r1.s(0.0) == pytest.approx(r0.s(0.0), rel=1e-8)
    assert r1.dy_dphi(500.0) == pytest.approx(r0.dy_dphi(500.0), rel=1e-7)


@pytest.mark.parametrize("uv", [(30.0, 0.0), (-40.0, 0.0), (0.0, 30.0), (20.0, -25.0)])
@pytest.mark.parametrize("phi", [0.0, 0.4, -0.5])
def test_uniform_wind_is_galilean_shift_geometry(uv, phi):
    """Uniform air + uniform wind: the air-frame ray is the windless one (azimuth phi about
    the air velocity); in the ground frame it is rotated by the crab angle and advected by
    w T; tube area and travel time are unchanged; Omega = c/(c + w.n)."""
    u, v = uv
    w = PR.WindProfile.uniform(u, v)
    V = 2.0 * C0
    rw = PR.RayPath(15000.0, V, atmosphere=UATM, phi=phi, wind=w, allow_cutoff=True)
    rn = PR.RayPath(15000.0, V, atmosphere=UATM, phi=phi, allow_cutoff=True)
    T = rw.travel_time(0.0)
    assert T == pytest.approx(rn.travel_time(0.0), rel=1e-12)
    e = rw.yaw
    assert np.sin(e) == pytest.approx(-v / V, rel=1e-12, abs=1e-15)
    xn, yn = rn.x(0.0), rn.y(0.0)
    assert rw.x(0.0) == pytest.approx(xn * np.cos(e) - yn * np.sin(e) + u * T, abs=1e-6)
    assert rw.y(0.0) == pytest.approx(xn * np.sin(e) + yn * np.cos(e) + v * T, abs=1e-6)
    z = np.linspace(0.0, 14000.0, 8)
    assert np.allclose(rw.area(z), rn.area(z), rtol=1e-10)
    assert rw.s(0.0) == pytest.approx(rn.s(0.0), rel=1e-12)
    n = rw.n_launch
    assert rw.omega(0.0) == pytest.approx(C0 / (C0 + u * n[0] + v * n[1]), rel=1e-12)
    assert rw.sx * rw.V_ground == pytest.approx(1.0, rel=1e-13)
    assert rw.V_ground == pytest.approx(np.sqrt(V ** 2 - v ** 2) + u, rel=1e-13)


@pytest.mark.parametrize("uv", [(40.0, 0.0), (-40.0, 0.0), (0.0, 35.0)])
@pytest.mark.parametrize("mode", ["linear", "nonlinear"])
def test_uniform_wind_is_galilean_shift_signature(uv, mode):
    """Absorbing (and nonlinear) march: the ground signature in uniform wind equals the
    windless one with its time axis compressed by Omega (a ground microphone sees the
    pattern sweep by at the ground speed); amplitude unchanged.  Linear + absorption
    checks the intrinsic-frequency filter (exact); the nonlinear case checks the
    Omega-scaled Burgers step (agreement limited by shock capture: 0.05 % at 200 kHz,
    < 1 % at 50 kHz)."""
    ac = dataclasses.replace(CONCORDE, mach=2.0, altitude=4000.0)
    phi = 0.3
    hum = lambda z: np.full_like(np.asarray(z, float), 1.0)  # noqa: E731
    kw = dict(atmosphere=UATM, humidity=hum, fs=50e3, phi=phi, linear=(mode == "linear"))
    tn, pn = PR.ground_signature(ac, **kw)
    w = PR.WindProfile.uniform(*uv)
    tw, pw = PR.ground_signature(ac, wind=w, **kw)
    om = float(PR.RayPath(ac.altitude, ac.V, atmosphere=UATM, phi=phi, wind=w).omega(0.0))
    assert abs(om - 1.0) > 0.02
    rel = 1e-9 if mode == "linear" else 1e-2
    assert pw.max() == pytest.approx(pn.max(), rel=rel)
    assert pw.min() == pytest.approx(pn.min(), rel=rel)
    span_n = tn[np.argmin(pn)] - tn[np.argmax(pn)]
    span_w = tw[np.argmin(pw)] - tw[np.argmax(pw)]
    assert span_w == pytest.approx(om * span_n, rel=1e-9 if mode == "linear" else 2e-3)
    if mode == "linear":                 # whole waveform, on the scaled time grid
        assert np.allclose(tw, tn * om, rtol=0, atol=1e-12)
        assert np.max(np.abs(pw - pn)) < 1e-9 * pn.max()


@pytest.mark.parametrize("wind", [PR.jet_stream_winter(0.0), PR.jet_stream_winter(180.0),
                                  PR.jet_stream_winter(90.0), PR.qbo_easterly(180.0)],
                         ids=["tail", "head", "cross", "qbo"])
@pytest.mark.parametrize("phi", [0.0, 0.5, -0.7])
def test_wind_jacobian_matches_neighbouring_rays(wind, phi):
    r = PR.RayPath(AC.altitude, AC.V, phi=phi, wind=wind, allow_cutoff=True)
    d = 1e-5
    a = PR.RayPath(AC.altitude, AC.V, phi=phi + d, wind=wind, allow_cutoff=True)
    b = PR.RayPath(AC.altitude, AC.V, phi=phi - d, wind=wind, allow_cutoff=True)
    for z in (0.0, 8000.0):
        fd = (a.y(z) - b.y(z)) / (2 * d)
        assert r.dy_dphi(z) == pytest.approx(fd, rel=1e-6)


def test_blokhintzev_pierce_form():
    """amplitude_ratio equals Pierce's form p^2 |v_ray| A_n/(rho c^2 Omega) = const with the
    ray-normal tube area A_n built independently from neighbouring traced rays."""
    w = PR.jet_stream_winter(60.0)
    phi, d = 0.35, 1e-5
    r = PR.RayPath(AC.altitude, AC.V, phi=phi, wind=w)
    nb = [PR.RayPath(AC.altitude, AC.V, phi=phi + s * d, wind=w) for s in (1, -1)]

    def pierce(z):
        a = r.atm(z)
        u, v = w(z)
        c = float(a.c)
        om = 1.0 - u * r.sx - v * r.sy
        sz = np.sqrt(om ** 2 / c ** 2 - r.sx ** 2 - r.sy ** 2)
        s = np.array([r.sx, r.sy, -sz])
        vr = c ** 2 * s / om + np.array([u, v, 0.0])
        # tube: e1 = launch-time step (dx0 = 1 along x), e2 = d(position)/dphi at fixed z
        e1 = np.array([1.0, 0.0, 0.0])
        e2 = np.array([(nb[0].x(z) - nb[1].x(z)) / (2 * d), (nb[0].y(z) - nb[1].y(z)) / (2 * d), 0.0])
        An = abs(np.dot(np.cross(e1, e2), vr / np.linalg.norm(vr)))
        return np.linalg.norm(vr) * An / (float(a.rho) * c ** 2 * om)

    for z1, z2 in ((19000.0, 11000.0), (11000.0, 0.0), (17000.0, 2000.0)):
        expect = np.sqrt(pierce(z1) / pierce(z2))
        assert float(r.amplitude_ratio(z1, z2)) == pytest.approx(expect, rel=1e-5)


def test_crosswind_mirror_symmetry_and_asymmetric_carpet():
    wp = PR.jet_stream_winter(90.0)
    wm = PR.jet_stream_winter(-90.0)
    for phi in (0.0, 0.4):
        a = PR.RayPath(AC.altitude, AC.V, phi=phi, wind=wp, allow_cutoff=True)
        b = PR.RayPath(AC.altitude, AC.V, phi=-phi, wind=wm, allow_cutoff=True)
        assert a.y_ground == pytest.approx(-b.y_ground, abs=1e-6)
    pp, yp = PR.lateral_cutoff(AC.altitude, AC.V, wind=wp, side=1)
    pn, yn = PR.lateral_cutoff(AC.altitude, AC.V, wind=wp, side=-1)
    assert yp > 0 > yn
    assert yp > -yn                       # carpet pushed downwind
    p0, y0 = PR.lateral_cutoff(AC.altitude, AC.V)
    qp, qy = PR.lateral_cutoff(AC.altitude, AC.V, wind=PR.calm(), side=1)
    assert qp == pytest.approx(p0, abs=1e-6)
    assert qy == pytest.approx(y0, rel=1e-3)


def test_head_and_tail_ground_speed():
    for toward, sgn in ((0.0, 1.0), (180.0, -1.0)):
        r = PR.RayPath(AC.altitude, AC.V, wind=PR.jet_stream_winter(toward))
        ua = float(PR.jet_stream_winter(toward)(AC.altitude)[0])
        assert sgn * ua > 0
        assert r.V_ground == pytest.approx(AC.V + ua, rel=1e-12)
        assert r.arrival(0.0)["delay"] == pytest.approx(r.travel_time(0.0) - r.x(0.0) / r.V_ground)


def test_propagate_wind_kwarg_rebuilds_ray():
    w = PR.WindProfile.uniform(25.0, 10.0)
    ray = PR.RayPath(3000.0, 2 * C0, atmosphere=UATM, phi=0.2)
    rw = ray.with_wind(w)
    t = np.arange(4000) / 50e3
    p = 50.0 * np.exp(-((t - 0.04) / 0.004) ** 2)
    _, a = PR.propagate(t, p, ray, 2900.0, 0.0, wind=w)
    _, b = PR.propagate(t, p, rw, 2900.0, 0.0)
    _, c = PR.propagate(t, p, ray, 2900.0, 0.0)
    assert np.array_equal(a, b)
    assert not np.allclose(a, c)


def test_wind_profiles_shapes():
    jet = PR.jet_stream_winter(0.0)
    u, v = jet(np.array([0.0, 11000.0, 20000.0]))
    assert 50 <= u[1] <= 70 and 10 <= u[2] <= 20 and np.all(v == 0)
    q = PR.qbo_easterly(180.0)
    u, _ = q(np.array([20000.0]))
    assert -30 <= u[0] <= -20
    c = PR.jet_stream_winter(90.0)
    u, v = c(np.array([11000.0]))
    assert abs(u[0]) < 1e-9 and v[0] == pytest.approx(65.0)
