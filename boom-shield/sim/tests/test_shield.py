import os
import sys
from types import SimpleNamespace

import numpy as np
import pytest

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))

from boomsim.shield import (GasLayer, Membrane, Screen, Stack, apply, gas,  # noqa: E402
                            gas_mixture, lift_margin, transfer)

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


# ------------------------------------------------------------------ tie stiffness
def _budget(g, d):
    return (ATM.rho - g.rho) * d


def _tied(g, d, K, m1=None, m2=None):
    b = _budget(g, d)
    m1 = b / 2 if m1 is None else m1
    m2 = b / 2 if m2 is None else m2
    return Stack([Membrane(m1), GasLayer(g, d, tie_K=K), Membrane(m2)])


def test_tie_zero_is_untied():
    f = np.logspace(-1, 4, 300)
    for V in (1.4 * ATM.c, 2.0 * ATM.c, 3.2 * ATM.c, np.inf):
        for g, d in ((HE, 15.0), (H2, 60.0), (HE, 300.0)):
            T0, R0 = transfer(_tied(g, d, 0.0), f, V, ATM)
            b = _budget(g, d)
            T1, R1 = transfer(Stack([Membrane(b / 2), GasLayer(g, d), Membrane(b / 2)]), f, V, ATM)
            np.testing.assert_array_equal(T0, T1)
            np.testing.assert_array_equal(R0, R1)


@pytest.mark.parametrize("K", [1.0, 35.0, 300.0, 1e4, 1e8])
@pytest.mark.parametrize("V", [1.4 * 295.07, 2.0 * 295.07, 3.2 * 295.07, np.inf])
def test_tie_energy_conservation(K, V):
    f = np.logspace(-1, 4, 300)
    stack = Stack([Membrane(0.3, N=50.0), GasLayer(HE, 30.0, tie_K=K), Membrane(0.2),
                   GasLayer(H2, 300.0, tie_K=K / 3), Membrane(0.05)])
    T, R = transfer(stack, f, V, ATM)
    assert np.all(np.isfinite(T)) and np.all(np.isfinite(R))
    np.testing.assert_allclose(np.abs(T) ** 2 + np.abs(R) ** 2, 1.0, atol=1e-9)


@pytest.mark.parametrize("mach", [1.4, 2.0])
@pytest.mark.parametrize("g,d", [(HE, 15.0), (H2, 30.0), (HE, 60.0)])
def test_rigid_tie_mass_law(mach, g, d):
    """K -> infinity: faces move together; at low f the stack is one limp sheet of mass
    m1 + m2 + rho_g d (the gas between rigidly coupled faces is carried along and its
    (negative) compliance is short-circuited)."""
    V = mach * ATM.c
    b = _budget(g, d)
    f = np.array([0.2, 0.5, 1.0, 2.0])
    T, _ = transfer(_tied(g, d, 1e9, 0.3 * b, 0.7 * b), f, V, ATM)
    Tm, _ = transfer(Stack([Membrane(b + g.rho * d)]), f, V, ATM)
    np.testing.assert_allclose(np.abs(T), np.abs(Tm), rtol=2e-2)
    # and it is close to the plain mass law of the faces alone (gas mass is small)
    Tf, _ = transfer(Stack([Membrane(b)]), f, V, ATM)
    assert np.all(np.abs(20 * np.log10(np.abs(T) / np.abs(Tf))) < 1.0)


def test_tie_large_kappa_d_no_overflow():
    V = 2.0 * ATM.c
    for K in (1e-3, 1.0, 1e3):
        T, R = transfer(_tied(HE, 300.0, K), [1e3, 1e4], V, ATM)
        assert np.all(np.isfinite(T)) and np.all(np.isfinite(R))
        np.testing.assert_allclose(np.abs(T) ** 2 + np.abs(R) ** 2, 1.0, atol=1e-9)


# Reference values: independent direct 2x2 face-velocity solve from the adversarial
# review (results/review_evanescent.md section 3.1; scratch misc.py ``tie``), M 2, 10 Hz.
@pytest.mark.parametrize("g,d,K,ref_db", [
    (HE, 15.0, 0.0, -17.8), (HE, 15.0, 100.0, -19.8), (HE, 15.0, 1e3, -10.8),
    (HE, 15.0, 1e4, -3.4), (HE, 60.0, 1e3, -26.3), (HE, 60.0, 1e4, -8.8),
])
def test_tie_matches_review_values(g, d, K, ref_db):
    T, _ = transfer(_tied(g, d, K), [10.0], 2.0 * ATM.c, ATM)
    assert 20 * np.log10(abs(T[0])) == pytest.approx(ref_db, abs=0.06)


@pytest.mark.parametrize("mach,g,d,K6", [
    (2.0, HE, 60.0, 35.0), (2.0, H2, 30.0, 105.0), (1.4, HE, 15.0, 275.0), (2.0, HE, 15.0, 955.0),
])
def test_tie_6db_loss_threshold_matches_review(mach, g, d, K6):
    """Review table 3.1: smallest K at which 10 Hz transmission is 6 dB above K = 0."""
    V = mach * ATM.c
    T0 = 20 * np.log10(abs(transfer(_tied(g, d, 0.0), [10.0], V, ATM)[0][0]))
    Ks = np.logspace(-1, 5, 601)
    TT = np.array([20 * np.log10(abs(transfer(_tied(g, d, K), [10.0], V, ATM)[0][0])) for K in Ks])
    K_found = Ks[np.argmax(TT > T0 + 6)]
    assert K_found == pytest.approx(K6, rel=0.03)


def test_tie_static_cancellation_scale():
    """At K = rho_g/(a^2 d) the tie cancels the evanescent gas layer's negative static
    stiffness: the layer is infinitely compliant at low f and T -> 0 there."""
    V = 2.0 * ATM.c
    d = 30.0
    a2 = 1 / V ** 2 - 1 / HE.c ** 2
    Kc = HE.rho / (a2 * d)
    f = np.array([0.05, 0.1])
    Tc, _ = transfer(_tied(HE, d, Kc), f, V, ATM)
    T0, _ = transfer(_tied(HE, d, 0.0), f, V, ATM)
    assert np.all(np.abs(Tc) < 0.05) and np.all(np.abs(T0) > 0.9)


# ------------------------------------------------------------------- gas mixtures
def test_gas_mixture():
    pure = gas_mixture({"helium": 1.0}, T20, P20)
    assert pure.c == pytest.approx(HE.c, rel=1e-12) and pure.rho == pytest.approx(HE.rho, rel=1e-12)
    airmix = gas_mixture({"air": 1.0}, T20, P20)
    assert airmix.c == pytest.approx(ATM.c, rel=1e-12)
    # review section 3.4: 5 % air drops c_He to ~751 m/s and c_H2 to ~867 m/s
    he5 = gas_mixture({"helium": 0.95, "air": 0.05}, T20, P20)
    h25 = gas_mixture({"hydrogen": 0.95, "air": 0.05}, T20, P20)
    assert he5.c == pytest.approx(751.3, abs=0.5)
    assert h25.c == pytest.approx(867.4, abs=0.5)
    # ideal mixing: M linear in mole fraction, density linear, gamma between components
    he10 = gas_mixture({"helium": 0.9, "air": 0.1}, T20, P20)
    assert he10.M == pytest.approx(0.9 * 4.0026e-3 + 0.1 * 28.964e-3)
    assert he10.rho == pytest.approx(0.9 * HE.rho + 0.1 * ATM.rho, rel=1e-12)
    assert 1.4 < he10.gamma < 5 / 3
    # cv additivity check: monatomic + diatomic 50/50 -> cv = 2 R -> gamma = 1.5
    assert gas_mixture({"helium": 0.5, "air": 0.5}, T20, P20).gamma == pytest.approx(1.5)
    with pytest.raises(ValueError):
        gas_mixture({"helium": 0.9}, T20, P20)
    with pytest.raises(ValueError):
        gas_mixture({"neon": 1.0}, T20, P20)
    # contaminated layer gives less lift
    assert lift_margin(Stack([GasLayer(he5, 15.0)]), ATM) < lift_margin(Stack([GasLayer(HE, 15.0)]), ATM)
