"""Tests for boomsim.focus (unsteady-flight ray caustics and fold-caustic focus amplitude)."""
import math

import numpy as np
import pytest

from boomsim import focus as F
from boomsim.atmosphere import std_atm
from boomsim.propagate import RayPath, uniform_atmosphere


@pytest.fixture(scope="module")
def ucol():
    return F.Column(atmosphere=uniform_atmosphere(0.0))


# ------------------------------------------------------------------ column integrals
@pytest.mark.parametrize("mach", [1.2, 1.7, 2.5])
def test_column_integrals_match_raypath(mach):
    za = 20250.0
    V = mach * float(std_atm(za).c)
    rp = RayPath(za, V, z_min=0.0)
    I = F.default_column().integrals(1.0 / V, 0.0, za)
    assert bool(I["ok"])
    assert float(I["I1"]) == pytest.approx(rp.I1_tab[0], rel=1e-5)
    assert float(I["TT"]) == pytest.approx(rp.tt_tab[0], rel=1e-5)


def test_cutoff_mach_standard_atmosphere():
    col = F.default_column()
    assert col.cutoff_mach(20250.0) == pytest.approx(340.294 / 295.1, rel=2e-3)
    # rays just below cut-off do not reach the ground
    c_a = float(col.c_at(20250.0))
    s = 1.0 / (0.995 * col.cutoff_mach(20250.0) * c_a)
    assert not bool(col.integrals(s, 0.0, 20250.0)["ok"])


# ------------------------------------------------------------------ Jacobian
def test_jacobian_matches_finite_differences():
    tr = F.Trajectory.level_turn(18000.0, 1.6, 20.0, 60.0, pad_s=20.0, roll_in_s=8.0)
    for te, ph in [(-3.0, 0.3), (10.0, -0.5), (25.0, 0.1)]:
        r = F.rays(tr, te, ph, 500.0)
        h_t, h_p = 1e-3, 1e-5
        dPt = (F.rays(tr, te + h_t, ph, 500.0)["P"] - F.rays(tr, te - h_t, ph, 500.0)["P"]) / (2 * h_t)
        dPp = (F.rays(tr, te, ph + h_p, 500.0)["P"] - F.rays(tr, te, ph - h_p, 500.0)["P"]) / (2 * h_p)
        np.testing.assert_allclose(r["dPdt"], dPt, rtol=1e-4, atol=1e-3)
        np.testing.assert_allclose(r["dPdp"], dPp, rtol=1e-4, atol=1e-2)


def test_jacobian_climb_matches_finite_differences():
    tr = F.Trajectory.climb_accel(12000.0, 20000.0, 0.95, 1.7, 0.05, pad_s=10.0)
    te = tr.time_at_mach(1.4)
    r = F.rays(tr, te, 0.2, 0.0)
    h = 1e-3
    dPt = (F.rays(tr, te + h, 0.2, 0.0)["P"] - F.rays(tr, te - h, 0.2, 0.0)["P"]) / (2 * h)
    np.testing.assert_allclose(r["dPdt"], dPt, rtol=1e-4, atol=1e-3)


def test_steady_level_flight_has_unit_gain_and_no_caustic():
    tr = F.Trajectory.level_accel(20250.0, 1.7, 1.7 + 1e-9, 1e-6, pad_s=50.0)
    phis = np.radians([0.0, 20.0, 40.0])
    g = F.focus_geometry(tr, 0.0 * phis - 20.0, phis, 0.0)
    np.testing.assert_allclose(g["G"], 1.0, rtol=1e-6)
    sc = F.scan(tr, np.arange(-40, 0, 2.0), phis)
    assert sc["roots"] == []


# ------------------------------------------------------------------ caustics
@pytest.mark.parametrize("a", [0.5, 2.0])
def test_uniform_atmosphere_focus_mach(ucol, a):
    """Level acceleration in uniform air: undertrack ground focus where
    beta^3 = a H / c^2 (envelope of straight rays)."""
    H = 15000.0
    c = float(ucol.c[0])
    tr = F.Trajectory.level_accel(H, 0.95, 1.7, a / F.G0, column=ucol)
    sc = F.scan(tr, np.arange(tr.t0 + 1, tr.t1 - 1, 0.5), [0.0])
    assert len(sc["roots"]) == 1
    M = float(tr.state(sc["roots"][0][0])["M"])
    assert (M * M - 1) ** 1.5 == pytest.approx(a * H / c ** 2, rel=1e-6)


def test_fold_parameters_match_envelope_curvature(ucol):
    """R from the arrival-time split of the two post-focus rays equals the radius of
    curvature of the analytic ray envelope (straight rays: R = caustic radius)."""
    H, a = 15000.0, 1.0
    tr = F.Trajectory.level_accel(H, 0.95, 1.7, a / F.G0, column=ucol)
    te = F.scan(tr, np.arange(tr.t0 + 1, tr.t1 - 1, 0.5), [0.0])["roots"][0][0]
    fp = F.fold_parameters(tr, te, 0.0)
    ts = np.linspace(te - 5, te + 5, 4001)
    st = tr.state(ts)
    th = np.arccos(1 / st["M"])
    dth = np.gradient(th, ts)
    dxe = np.gradient(st["x"], ts)
    d = np.stack([np.cos(th), -np.sin(th)])
    nd = np.stack([np.cos(th), -np.sin(th)]) * dth
    lam = np.sin(th) * dxe / (nd[0] * d[0] + nd[1] * d[1])
    xc, zc = st["x"] + lam * d[0], H + lam * d[1]
    i = int(np.argmin(np.abs(zc)))
    dx, dz = np.gradient(xc, ts), np.gradient(zc, ts)
    k = (dx * np.gradient(dz, ts) - dz * np.gradient(dx, ts)) / (dx ** 2 + dz ** 2) ** 1.5
    assert fp["X"] == pytest.approx(xc[i], abs=5.0)
    assert fp["R"] == pytest.approx(1 / abs(k[i]), rel=2e-3)
    assert fp["alpha"] == pytest.approx(math.atan2(abs(dz[i]), abs(dx[i])), rel=1e-3)
    # geometric amplitude ~ xi^(-1/4): Lambda constant over the sampled distances
    lams = [b["Lam"] for b in fp["branches"]]
    assert np.ptp(lams) / np.mean(lams) < 0.01


def test_stratified_acceleration_focus_just_above_cutoff():
    tr = F.Trajectory.level_accel(20250.0, 0.95, 1.7, 0.05)
    sc = F.scan(tr, np.arange(tr.t0 + 1, tr.t1 - 1, 0.5), [0.0])
    assert len(sc["roots"]) == 1
    M = float(tr.state(sc["roots"][0][0])["M"])
    m_cut = F.default_column().cutoff_mach(20250.0)
    assert m_cut < M < m_cut + 0.06


def test_gentler_acceleration_focuses_closer_to_cutoff():
    Ms = []
    for a in (0.02, 0.1):
        tr = F.Trajectory.level_accel(20250.0, 0.95, 1.7, a)
        sc = F.scan(tr, np.arange(tr.t0 + 1, tr.t1 - 1, 0.5), [0.0])
        Ms.append(float(tr.state(sc["roots"][0][0])["M"]))
    assert Ms[0] < Ms[1]


def test_deceleration_has_no_caustic():
    tr = F.Trajectory.level_accel(20250.0, 1.7, 0.95, -0.1)
    sc = F.scan(tr, np.arange(tr.t0 + 1, tr.t1 - 1, 0.5), np.radians([0, 20, 40]))
    assert sc["roots"] == []
    assert np.nanmin(sc["J"]) > 0


def test_turn_caustic_inside_only_above_threshold_bank():
    phis = np.radians(np.r_[np.arange(-52.9, -45, 0.2), 0.0, np.arange(45, 53, 0.2)])
    out = {}
    for b in (10.0, 30.0):
        tr = F.Trajectory.level_turn(20250.0, 1.7, b, 60.0, pad_s=30.0, roll_in_s=10.0)
        sc = F.scan(tr, np.arange(10.0, 50.0, 10.0), phis, along_phi=True)
        out[b] = sc["roots_phi"]
    assert out[10.0] == []
    assert out[30.0] and all(ph > 0 for _, ph in out[30.0])     # inside of a left turn


# ------------------------------------------------------------------ Airy focus
def _pulse(fs=4000.0, f0=20.0):
    t = np.arange(-2.0, 2.0, 1.0 / fs)
    p = -np.gradient(np.exp(-(t * f0 * 2) ** 2), t) / (f0 * 4)   # zero-mean pulse
    return t, p / p.max()


def test_airy_far_field_recovers_geometric_incident():
    t, p = _pulse()
    L_c, R, c = 2000.0, 1e5, 340.0
    xi = 6000.0
    res = F.airy_focus(t, p, L_c, R, c, xi=[xi])
    q, tt = res["p"], res["t"]
    # incident arrives (4/3c) sqrt(2/R) xi^1.5 / 2 earlier than the caustic time origin
    dT = (4 / (3 * c)) * math.sqrt(2 / R) * xi ** 1.5
    early = tt < -0.25 * dT
    assert dT > 0.5
    assert q[early].max() == pytest.approx((L_c / xi) ** 0.25, rel=0.02)
    # equal energy in the reflected (Hilbert-shifted) branch
    late = tt > 0.25 * dT
    assert np.sum(q[late] ** 2) == pytest.approx(np.sum(q[early] ** 2), rel=0.02)


def test_airy_caustic_peak_scales_with_frequency_one_sixth():
    """Linear fold focus of a narrow-band pulse: peak/incident ~ (k^2 L_c^3/R)^(1/12)
    x const, so doubling the frequency raises it by 2^(1/6)."""
    L_c, R, c = 2000.0, 1e5, 340.0
    pk = []
    for f0 in (10.0, 20.0):
        t = np.arange(-3, 3, 1 / 8000.0)
        p = np.cos(2 * np.pi * f0 * t) * np.exp(-(t * f0 / 4) ** 2)
        pk.append(F.airy_focus(t, p, L_c, R, c)["peak"])
    assert pk[1] / pk[0] == pytest.approx(2 ** (1 / 6), rel=0.03)


# ------------------------------------------------------------------ NKZ
def _nwave(T=0.1, rise=3e-3, fs=10e3, amp=20.0):
    t = np.arange(-0.3, 0.4, 1 / fs)
    p = np.where((t > 0) & (t < T), amp * (1 - 2 * t / T), 0.0)
    from scipy.ndimage import gaussian_filter1d
    return t, gaussian_filter1d(p, rise / 2.56 * fs)


def test_nkz_linear_is_the_airy_steady_state():
    t, p = _nwave()
    n = F.nkz_focus(t, p, 1500.0, 6e4, 340.0, 1.2, linear=True, X_end=0.5)
    assert n["peak_nl"] == pytest.approx(n["peak_lin"], rel=1e-9)
    a = F.airy_focus(t, p, 1500.0, 6e4, 340.0)
    # band-limited grid vs full-band Airy: within 15 %
    assert n["peak_lin"] == pytest.approx(a["peak"], rel=0.15)


def test_nkz_weak_wave_is_linear_and_strong_wave_saturates():
    t, p = _nwave(amp=1.0)
    weak = F.nkz_focus(t, 0.05 * p, 1500.0, 6e4, 340.0, 1.2, X_end=2.0)
    assert weak["ratio"] == pytest.approx(1.0, abs=0.03)
    strong = F.nkz_focus(t, 200.0 * p, 1500.0, 6e4, 340.0, 1.2, X_end=3.0)
    assert strong["ratio"] < 0.9
