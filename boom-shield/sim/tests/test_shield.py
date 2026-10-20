import os
import sys
from types import SimpleNamespace

import numpy as np
import pytest

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))

from boomsim.shield import (GasLayer, Membrane, Screen, Stack, apply, gas,  # noqa: E402
                            lift_margin, transfer)

# 20 km ambient (US Std 1976): 216.65 K, 5474.9 Pa
T20, P20 = 216.65, 5474.9
_air = gas("air", T20, P20)
ATM = SimpleNamespace(T=T20, p=P20, rho=_air.rho, c=_air.c)
HE = gas("helium", T20, P20)
H2 = gas("hydrogen", T20, P20)


def test_gas_properties():
    assert ATM.rho == pytest.approx(0.0880, rel=2e-3)
    assert ATM.c == pytest.approx(295.07, rel=1e-3)
    assert HE.c == pytest.approx(866, rel=2e-3)
    assert H2.c == pytest.approx(1120, rel=5e-3)
    thin = gas("air", T20, 0.1 * P20)
    assert thin.c == pytest.approx(ATM.c)
    assert thin.rho == pytest.approx(0.1 * ATM.rho)


def test_mass_law_normal_incidence():
    m = 0.5
    f = np.logspace(-1, 4, 200)
    T, R = transfer(Stack([Membrane(m)]), f, np.inf, ATM)
    w = 2 * np.pi * f
    expected = 1.0 / (1.0 + 1j * w * m / (2 * ATM.rho * ATM.c))
    np.testing.assert_allclose(T, expected, rtol=1e-10)
    # very large finite V converges too
    T2, _ = transfer(Stack([Membrane(m)]), f, 1e9, ATM)
    np.testing.assert_allclose(T2, expected, rtol=1e-8)


@pytest.mark.parametrize("mach", [1.2, 2.0, 3.2])
def test_mass_law_oblique(mach):
    m = 1.1
    V = mach * ATM.c
    cos_t = np.sqrt(1 - (ATM.c / V) ** 2)
    f = np.logspace(-1, 3, 100)
    T, _ = transfer(Stack([Membrane(m)]), f, V, ATM)
    w = 2 * np.pi * f
    np.testing.assert_allclose(np.abs(T), 1 / np.abs(1 + 1j * w * m * cos_t / (2 * ATM.rho * ATM.c)),
                               rtol=1e-10)


def test_tensioned_membrane_coincidence():
    # c_m = V -> Z_m = 0 (lossless) -> transparent
    V = 2.0 * ATM.c
    m = 0.2
    T, _ = transfer(Stack([Membrane(m, N=m * V ** 2)]), [1.0, 10.0, 100.0], V, ATM)
    np.testing.assert_allclose(np.abs(T), 1.0, atol=1e-12)


@pytest.mark.parametrize("V", [1.4 * 295.07, 2.0 * 295.07, 3.2 * 295.07, np.inf])
@pytest.mark.parametrize("g", [HE, H2, gas("air", T20, 0.1 * P20)])
def test_energy_conservation(V, g):
    f = np.logspace(-1, 4, 300)
    stack = Stack([Membrane(0.3, N=50.0), GasLayer(g, 30.0), Membrane(0.2),
                   GasLayer(_air, 2.0), Membrane(0.05)])
    T, R = transfer(stack, f, V, ATM)
    assert np.all(np.isfinite(T)) and np.all(np.isfinite(R))
    np.testing.assert_allclose(np.abs(T) ** 2 + np.abs(R) ** 2, 1.0, atol=1e-9)


def test_lossy_elements_dissipate():
    f = np.logspace(0, 3, 50)
    for el in (Screen(200.0), Membrane(0.3, N=500.0, eta=0.1)):
        T, R = transfer(Stack([el, GasLayer(HE, 10.0), Membrane(0.1)]), f, 2 * ATM.c, ATM)
        assert np.all(np.abs(T) ** 2 + np.abs(R) ** 2 < 1.0 + 1e-12)


def test_evanescent_tunnelling_slope():
    V = 2.0 * ATM.c  # 590 m/s < c_He
    assert V < HE.c
    f = 5.0
    w = 2 * np.pi * f
    kappa = w * np.sqrt(1 / V ** 2 - 1 / HE.c ** 2)
    ds = np.linspace(100, 300, 9)  # kappa d from ~3.8 to ~11
    T = np.array([transfer(Stack([Membrane(1e-4), GasLayer(HE, d), Membrane(1e-4)]), f, V, ATM)[0][0]
                  for d in ds])
    slope = np.polyfit(ds, np.log(np.abs(T)), 1)[0]
    assert slope == pytest.approx(-kappa, rel=1e-3)
    # no overflow at absurd kappa d
    T_far, R_far = transfer(Stack([Membrane(1e-4), GasLayer(HE, 300.0), Membrane(1e-4)]),
                            [1e4], V, ATM)
    assert np.isfinite(T_far[0]) and abs(T_far[0]) < 1e-100
    assert abs(R_far[0]) == pytest.approx(1.0, abs=1e-9)


def test_double_wall_resonance():
    m1, m2, d = 0.1, 0.15, 0.1
    rc2 = ATM.rho * ATM.c ** 2
    f0 = np.sqrt(rc2 * (m1 + m2) / (d * m1 * m2)) / (2 * np.pi)
    assert d < 0.1 * ATM.c / f0  # d << lambda
    f = np.linspace(0.5 * f0, 1.5 * f0, 20001)
    T, _ = transfer(Stack([Membrane(m1), GasLayer(_air, d), Membrane(m2)]), f, np.inf, ATM)
    f_peak = f[np.argmax(np.abs(T))]
    assert f_peak == pytest.approx(f0, rel=0.03)
    # well above f0 transmission falls off much faster than single-wall mass law
    w = 2 * np.pi * 4 * f0
    T4, _ = transfer(Stack([Membrane(m1), GasLayer(_air, d), Membrane(m2)]), [4 * f0], np.inf, ATM)
    single = abs(1 / (1 + 1j * w * (m1 + m2) / (2 * ATM.rho * ATM.c)))
    assert abs(T4[0]) < single


def _nwave(fs=20000.0, dur=0.1, pad=0.5):
    t = np.arange(-pad, dur + pad, 1 / fs)
    p = np.where((t >= 0) & (t < dur), 1.0 - 2.0 * t / dur, 0.0)
    return t, p


def test_apply_identity_and_dc():
    t, p = _nwave()
    for stack in (Stack([]), Stack([GasLayer(_air, 25.0)])):
        out = apply(stack, t, p, 2 * ATM.c, ATM)
        assert out.dtype.kind == "f"
        np.testing.assert_allclose(out, p, atol=1e-9)
    # DC gain 1 for membranes and gas layers
    for st in (Stack([Membrane(1.0)]),
               Stack([Membrane(0.5), GasLayer(HE, 15.0), Membrane(0.5)]),
               Stack([Membrane(0.5), GasLayer(H2, 120.0), Membrane(0.5)])):
        T0, R0 = transfer(st, [0.0, 1e-4], 2 * ATM.c, ATM)
        np.testing.assert_allclose(T0[0], 1.0, atol=1e-6)
        np.testing.assert_allclose(T0[1], 1.0, atol=1e-3)
        assert abs(R0[0]) < 1e-6
        # pulse area (DC) preserved. Evanescent layers smooth with a heavy-tailed
        # (Lorentzian ~1/t^2) kernel, so use a long window to keep the tails.
        tt = np.arange(-20.0, 20.0, 2e-3)
        pulse = np.exp(-(tt / 0.05) ** 2)
        out = apply(st, tt, pulse, 2 * ATM.c, ATM)
        assert np.trapezoid(out, tt) == pytest.approx(np.trapezoid(pulse, tt), rel=1e-2)


def test_lift_margin():
    he = gas("helium", T20, P20)
    assert lift_margin(Stack([GasLayer(he, 15.0)]), ATM) == pytest.approx(1.138, abs=0.01)
    st = Stack([Membrane(0.2), GasLayer(he, 15.0), Membrane(0.3), Screen(100.0, m=0.05)])
    assert lift_margin(st, ATM) == pytest.approx((ATM.rho - he.rho) * 15.0 - 0.55)
    assert lift_margin(Stack([GasLayer(_air, 10.0)]), ATM) == pytest.approx(0.0)
