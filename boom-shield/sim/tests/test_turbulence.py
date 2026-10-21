"""Phase-screen / Kirchhoff turbulence kernel (boomsim.turbulence)."""
import os
import sys

import numpy as np
import pytest

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))

from boomsim import turbulence as TU  # noqa: E402

FS = 50e3
SMALL = TU.with_(TU.PRESETS["moderate"], n=512, r_turb=60.0, r_blend=20.0, r_ap_in=90.0,
                 r_ap=120.0)


def _grid(p):
    x = (np.arange(p.n) - p.n // 2) * p.dx
    X, Y = np.meshgrid(x, x, indexing="ij")
    return x, X, Y


def _pulse(sigma=0.5e-3, n=20000):
    t = (np.arange(n) - n // 2) / FS
    return t, np.exp(-0.5 * (t / sigma) ** 2)


def test_no_turbulence_is_identity():
    x, X, _ = _grid(SMALL)
    k, i0 = TU.screen_kernel(np.zeros_like(X), SMALL, FS, x=x)
    assert k[i0] == 1.0
    assert np.all(np.delete(k, i0) == 0.0)


def test_unperturbed_aperture_area_is_exact_ramp():
    """Area of {r^2/(2cd) < t} = 2 pi c d t: the second difference is a unit spike."""
    x, X, Y = _grid(SMALL)
    geo = (X ** 2 + Y ** 2) / (2 * SMALL.c * SMALL.d)
    dt = 1 / FS
    disk = (np.hypot(X, Y) <= 100.0).astype(float)
    A = TU._area_curve(geo, disk, x, dt, -3 * dt, 1700)
    t = -3 * dt + np.arange(1700) * dt
    ramp = 2 * np.pi * SMALL.c * SMALL.d * np.maximum(t, 0)
    # linear interpolation of r^2 on 0.5 m triangles: constant area offset (0.25 m^2)
    assert np.allclose(A[:300], ramp[:300], rtol=0, atol=1e-4 * ramp[299])
    # per-bin areas: grid (piecewise-linear contour) noise of order 0.1 %
    assert np.allclose(np.diff(A[10:300]), np.diff(ramp[10:300]), rtol=2e-3)


@pytest.mark.parametrize("F", [1500.0, 3000.0, -2000.0])
def test_parabolic_lens_focal_gain(F):
    """dT = -r^2/(2cF): effective distance d' = 1/(1/d - 1/F), amplitude gain d'/d."""
    x, X, Y = _grid(SMALL)
    k, i0 = TU.screen_kernel(-(X ** 2 + Y ** 2) / (2 * SMALL.c * F), SMALL, FS, x=x)
    _, p = _pulse()
    out = TU.apply_kernel(p, k, i0)
    assert out.max() == pytest.approx(F / (F - SMALL.d), rel=3e-3)


def test_tilt_is_pure_delay():
    x, X, _ = _grid(SMALL)
    k, i0 = TU.screen_kernel(2e-5 * X, SMALL, FS, x=x)
    t, p = _pulse()
    out = TU.apply_kernel(p, k, i0)
    assert out.max() == pytest.approx(1.0, rel=2e-3)
    assert np.sum(out) == pytest.approx(np.sum(p), rel=1e-6)


def test_screen_variance_and_kernel_dc():
    """Spectral synthesis reproduces the analytic rms when the box >> L0 (here only the
    screen statistics are tested, so a coarse grid is fine); kernels conserve the mean
    (DC gain 1), so low frequencies pass unchanged."""
    p = TU.with_(SMALL, n=1024, dx=4.0, L0=60.0)
    rng = np.random.default_rng(0)
    var = np.mean([np.var(TU.von_karman_screen(p, rng)[1]) for _ in range(8)])
    assert np.sqrt(var) == pytest.approx(TU.delay_sigma(p), rel=0.1)
    k, i0, _ = TU.realisation(SMALL, FS, np.random.default_rng(1))
    assert np.sum(k) == pytest.approx(1.0, abs=1e-6)
    f = np.array([2.0, 5.0])
    j = np.arange(k.size) - i0
    H = np.array([np.sum(k * np.exp(-2j * np.pi * fi * j / FS)) for fi in f])
    assert np.all(np.abs(np.abs(H) - 1) < 0.03)


def test_turbulence_spikes_and_rounds_statistically():
    """A 0.5 ms pulse: some realisations focus (peak > 1), others defocus (< 1); the
    median peak is reduced (scattering spreads energy in delay)."""
    rng = np.random.default_rng(3)
    _, p = _pulse()
    peaks = []
    for _ in range(12):
        k, i0, _ = TU.realisation(SMALL, FS, rng)
        peaks.append(TU.apply_kernel(p, k, i0).max())
    peaks = np.array(peaks)
    assert peaks.min() < 0.9 and np.median(peaks) < 1.0
    assert peaks.max() > peaks.min() * 1.2
