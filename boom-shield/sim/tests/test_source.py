import os
import sys

import numpy as np
import pytest
from scipy.integrate import quad

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))

from boomsim.atmosphere import std_atm  # noqa: E402
from boomsim.source import (CONCORDE, SR71, Aircraft, equivalent_area,  # noqa: E402
                            equivalent_area_slope, f_function, lift_distribution,
                            near_field, volume_area)

AC = Aircraft("test", length=60.0, mach=2.0, altitude=17000.0, weight=1.5e6, volume=450.0)


def test_aircraft_properties():
    assert AC.V == pytest.approx(2.0 * float(std_atm(17000.0).c))
    assert AC.beta == pytest.approx(np.sqrt(3.0))


def test_sears_haack_volume_and_lift_area():
    x = np.linspace(0, AC.length, 200001)
    assert np.trapezoid(volume_area(AC, x), x) == pytest.approx(AC.volume, rel=1e-6)
    assert np.trapezoid(lift_distribution(AC, x), x) == pytest.approx(AC.weight, rel=1e-6)
    a = std_atm(AC.altitude)
    AL_end = AC.beta * AC.weight / (float(a.rho) * AC.V ** 2)
    assert float(equivalent_area(AC, AC.length)) == pytest.approx(AL_end, rel=1e-9)
    assert float(equivalent_area(AC, 2 * AC.length)) == pytest.approx(AL_end, rel=1e-9)
    # analytic slope agrees with numerical derivative of the area
    xm = np.linspace(1.0, AC.length - 1.0, 50)
    num = (equivalent_area(AC, xm + 1e-4) - equivalent_area(AC, xm - 1e-4)) / 2e-4
    assert np.allclose(equivalent_area_slope(AC, xm), num, rtol=1e-5, atol=1e-8)


def _f_reference(ac, y):
    """F via the A' form, integrated by parts with the substitution u = sqrt(y - xi):
    F = (1/2pi) [A'(y)/sqrt(0)... ] -> use F = (1/pi) d/dy int_0^y A'(xi) sqrt(y-xi)^-1 ... ;
    here simply F = (1/2pi) d/dy int_0^y A'(xi)/sqrt(y-xi) * 2 ... computed numerically:
    G(y) = int_0^y A'(xi)/sqrt(y-xi) dxi  =>  F = G'(y)/(2 pi)."""
    def G(yy):
        # substitution xi = yy - u^2 removes the singularity
        val, _ = quad(lambda u: 2.0 * equivalent_area_slope(ac, yy - u * u),
                      0.0, np.sqrt(yy), limit=400, points=None)
        return val
    d = 1e-3
    return (G(y + d) - G(y - d)) / (2 * d) / (2 * np.pi)


@pytest.mark.parametrize("frac", [0.1, 0.3, 0.55, 0.8, 1.4])
def test_f_function_against_quadrature(frac):
    y = frac * AC.length
    assert f_function(AC, y) == pytest.approx(_f_reference(AC, y), rel=2e-3, abs=2e-4)


def test_f_function_nose_limit():
    # Sears-Haack nose: A = 8 S (x/L)^{3/2} => A'' ~ C x^{-1/2}, C = 6 S / L^{3/2},
    # so F(0+) = (1/2pi) C pi = C/2 (finite: the near field starts with a shock).
    ac = Aircraft("sh", length=10.0, mach=2.0, altitude=15000.0, weight=0.0, volume=1.0)
    S = 16 * ac.volume / (3 * np.pi * ac.length)
    C = 6.0 * S / ac.length ** 1.5
    F = f_function(ac, np.array([0.005, 0.02]))
    assert np.allclose(F, C / 2.0, rtol=0.01)


def test_f_function_area_balance_and_tail():
    y = np.linspace(0, 6 * CONCORDE.length, 12001)
    F = f_function(CONCORDE, y)
    # F integrates to ~0 far behind the body (A_e' -> 0): residual ~ A_e(L)/(2 pi sqrt(Y))
    I = np.trapezoid(F, y)
    assert abs(I) < 0.15 * np.max(np.cumsum(F) * (y[1] - y[0]))
    # tail decays
    assert abs(F[-1]) < 0.02 * F.max()


def test_near_field_scaling():
    fs = 20e3
    t1, p1 = near_field(CONCORDE, 60.0, fs)
    t2, p2 = near_field(CONCORDE, 240.0, fs)
    assert t1[1] - t1[0] == pytest.approx(1 / fs)
    assert np.allclose(p2, p1 / 2.0)
    a = CONCORDE.atm
    y = CONCORDE.V * t1
    sel = (y > 0) & (y < 2.9 * CONCORDE.length)       # record holds F up to y = 3 L
    F = f_function(CONCORDE, y[sel])
    pref = float(a.p) * 1.4 * CONCORDE.mach ** 2 * F / np.sqrt(2 * CONCORDE.beta * 60.0)
    assert np.allclose(p1[sel], pref, atol=6e-3 * np.abs(pref).max())
    # zero padding at both ends
    assert p1[0] == 0.0 and p1[-1] == 0.0 and np.abs(p1).max() > 100.0


def test_reference_aircraft_shape_factor():
    # Carlson (NASA TP-1122) shape factor K_S = 2^{1/4} gamma (gamma+1)^{-1/2} sqrt(I)/l^{3/4},
    # I = max_y int_0^y F.  Carlson's chart gives ~0.09 for SST (Concorde-class) configurations.
    for ac, lo, hi in ((CONCORDE, 0.07, 0.11), (SR71, 0.10, 0.19)):
        y = np.linspace(0, 2 * ac.length, 8001)
        F = f_function(ac, y)
        I = np.max(np.cumsum(F) * (y[1] - y[0]))
        ks = 2 ** 0.25 * 1.4 / np.sqrt(2.4) * np.sqrt(I) / ac.length ** 0.75
        assert lo < ks < hi


# ------------------------------------------------------------------ azimuth (off-track)
def test_equivalent_area_azimuth():
    x = np.linspace(-5.0, 1.5 * AC.length, 777)
    from boomsim.source import lift_area
    assert np.array_equal(equivalent_area(AC, x, phi=0.0), equivalent_area(AC, x))
    assert np.array_equal(equivalent_area_slope(AC, x, phi=0.0), equivalent_area_slope(AC, x))
    for phi in (0.3, 0.9, -0.9):
        assert np.allclose(equivalent_area(AC, x, phi),
                           volume_area(AC, x) + np.cos(phi) * lift_area(AC, x), rtol=1e-12)
        xm = np.linspace(1.0, AC.length - 1.0, 40)
        num = (equivalent_area(AC, xm + 1e-4, phi) - equivalent_area(AC, xm - 1e-4, phi)) / 2e-4
        assert np.allclose(equivalent_area_slope(AC, xm, phi), num, rtol=1e-5, atol=1e-8)


def test_f_function_and_near_field_azimuth():
    y = np.linspace(0, 2.5 * AC.length, 501)
    F0 = f_function(AC, y)
    assert np.array_equal(f_function(AC, y, phi=0.0), F0)
    Fv = f_function(AC, y, phi=0.5 * np.pi)            # cos = 6e-17: volume only
    for phi in (0.4, 1.0):
        # F is linear in A_e: F(phi) = F_v + cos(phi) (F(0) - F_v)
        assert np.allclose(f_function(AC, y, phi=phi), Fv + np.cos(phi) * (F0 - Fv),
                           atol=1e-9 * np.abs(F0).max())
    t0, p0 = near_field(CONCORDE, 90.0, 10e3)
    t1, p1 = near_field(CONCORDE, 90.0, 10e3, phi=0.0)
    assert np.array_equal(p0, p1) and np.array_equal(t0, t1)
    # off-track the lift part (which dominates the Concorde F-function) is weaker
    _, p2 = near_field(CONCORDE, 90.0, 10e3, phi=0.8)
    assert p2.max() < p0.max()
