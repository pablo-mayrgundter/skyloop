"""Tests for boomsim.loudness (Stevens Mark VII PLdB and waveform metrics)."""

import os
import sys
import warnings

import numpy as np
import pytest
from scipy.signal import resample_poly

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))

from boomsim.loudness import PSF_TO_PA, band_levels, metrics, pldb  # noqa: E402

DATA = os.path.join(os.path.dirname(__file__), "data")


def nwave(dp_psf, T, tr, fs=50e3, pad=0.5):
    """Symmetric N-wave: peak ``dp_psf``, shock-to-shock duration ``T`` [s],
    linear-ramp shocks of rise time ``tr`` [s] (ramp ends at the nominal shock
    instant), ``pad`` seconds of zeros on each side."""
    dp = dp_psf * PSF_TO_PA
    t = np.arange(0.0, T + 2 * pad, 1.0 / fs)
    x = t - pad
    p = np.where((x >= 0) & (x <= T), dp * (1 - 2 * x / T), 0.0)
    m = (x >= -tr) & (x < 0)
    p[m] = dp * (x[m] + tr) / tr
    m = (x > T) & (x <= T + tr)
    p[m] = -dp * (1 - (x[m] - T) / tr)
    return t, p


# ---------------------------------------------------------------------------
# Reference anchors
# ---------------------------------------------------------------------------
def test_pyldb_reference_signature():
    """PyLdB regression signature (``misc/panair_r1.sig``, a PANAIR/sBOOM-type
    shaped ground signature, 0.38 psf peak) evaluated by PyLdB with
    ``pad_front=6, pad_rear=6, len_window=800`` gives 77.67985 PLdB
    (https://github.com/usuaero/PyLdB, ``verif/test.py``).  PyLdB's authors
    report agreement with NASA's LCASB within 0.2 PLdB (Bolander & Hunsaker,
    AIAA 2019-2091).  The record ends at -0.03 psf, so the same 800-sample
    half-Hann taper is required.  Achieved: 77.677 (diff 0.003 dB)."""
    d = np.load(os.path.join(DATA, "pyldb_panair_r1.npz"))
    t = d["time_ms"] * 1e-3
    p = d["pressure_psf"] * PSF_TO_PA
    dt = (t[-1] - t[0]) / (t.size - 1)
    assert pldb(t, p, taper_s=800 * dt) == pytest.approx(77.67985, abs=0.5)
    assert pldb(t, p, taper_s=800 * dt) == pytest.approx(77.67985, abs=0.05)
    # Without the taper the truncation step is (correctly) heard as a click.
    with pytest.warns(RuntimeWarning):
        assert pldb(t, p) > 79.0


@pytest.mark.parametrize("dp_psf,T,tr,ref", [
    (1.00, 0.3, 1e-3, 100.293),
    (1.00, 0.3, 3e-3, 95.799),
    (1.94, 0.2, 2e-3, 104.246),
    (0.30, 0.1, 5e-3, 81.999),
])
def test_matches_pyldb_on_nwaves(dp_psf, T, tr, ref):
    """Cross-check against PyLdB (usuaero/PyLdB @0e40b88, run with
    ``pad_front=pad_rear=10, len_window=0`` on the identical 50 kHz N-wave).
    Reference values computed offline with that code; achieved |diff| <= 0.012 dB."""
    assert pldb(*nwave(dp_psf, T, tr)) == pytest.approx(ref, abs=0.1)


def test_concorde_class_nwave():
    """NASA quotes Concorde's boom as ~105 PLdB at 1.94 psf
    (https://www.nasa.gov/aeronautics/sonic-boom-heads-for-a-thump/).
    Rise time/duration are not given there; for a 1.94 psf N-wave of 200-300 ms
    and 1-3 ms rise we get 102.4-106.9 PLdB.  Loose check: 105 +- 3 at 2 ms."""
    for T in (0.2, 0.3):
        assert pldb(*nwave(1.94, T, 2e-3)) == pytest.approx(105.0, abs=3.0)


# ---------------------------------------------------------------------------
# Sanity behaviour
# ---------------------------------------------------------------------------
def test_one_psf_nwave_level():
    """~1 psf, 300 ms, 1-3 ms rise lands in the mid/high 90s to ~100 PLdB."""
    for tr in (1e-3, 2e-3, 3e-3):
        assert 93.0 < pldb(*nwave(1.0, 0.3, tr)) < 103.0


def test_overpressure_scaling():
    """At fixed shape, 10x overpressure adds roughly 20 dB (Mark VII
    compresses slightly at high level; band-level growth is exactly 20 dB)."""
    a = pldb(*nwave(0.3, 0.2, 2e-3))
    b = pldb(*nwave(3.0, 0.2, 2e-3))
    assert 17.0 < b - a < 25.0
    _, la = band_levels(*nwave(0.3, 0.2, 2e-3))
    _, lb = band_levels(*nwave(3.0, 0.2, 2e-3))
    np.testing.assert_allclose(lb - la, 20.0, atol=1e-3)


def test_rise_time_dependence():
    """PL decreases monotonically with rise time, a few dB per doubling in the
    ms range (Shepherd & Sullivan 1991: loudness of N-waves controlled by
    overpressure and rise time)."""
    trs = [0.5e-3, 1e-3, 2e-3, 4e-3, 8e-3]
    pl = np.array([pldb(*nwave(1.0, 0.3, tr)) for tr in trs])
    steps = -np.diff(pl)
    assert np.all(steps > 1.0) and np.all(steps < 7.0)


def test_silence_is_minus_inf():
    t = np.arange(0, 1, 1 / 20e3)
    assert pldb(t, np.zeros_like(t)) == -np.inf


# ---------------------------------------------------------------------------
# Invariance
# ---------------------------------------------------------------------------
def test_zero_padding_invariance():
    ref = pldb(*nwave(1.0, 0.3, 2e-3, pad=0.1))
    for pad in (0.5, 2.0):
        assert abs(pldb(*nwave(1.0, 0.3, 2e-3, pad=pad)) - ref) < 0.2
    # asymmetric padding added by hand
    t, p = nwave(1.0, 0.3, 2e-3, pad=0.1)
    p2 = np.concatenate([np.zeros(3 * 50000), p, np.zeros(7919)])
    t2 = np.arange(p2.size) / 50e3
    assert abs(pldb(t2, p2) - ref) < 0.2


def test_sample_rate_invariance():
    ref = pldb(*nwave(1.0, 0.3, 2e-3, fs=100e3))
    for fs in (20e3, 25e3, 50e3, 96e3):
        assert abs(pldb(*nwave(1.0, 0.3, 2e-3, fs=fs)) - ref) < 0.2
    # band-limited resampling of the same record
    t, p = nwave(1.0, 0.3, 2e-3, fs=50e3)
    for up, down in ((2, 1), (3, 5)):
        pr = resample_poly(p, up, down)
        tr = np.arange(pr.size) * (t[1] - t[0]) * down / up
        assert abs(pldb(tr, pr) - pldb(t, p)) < 0.2


# ---------------------------------------------------------------------------
# metrics()
# ---------------------------------------------------------------------------
def test_metrics_nwave():
    T, tr, fs = 0.3, 2e-3, 50e3
    t, p = nwave(1.0, T, tr, fs=fs)
    m = metrics(t, p)
    assert m["peak_pa"] == pytest.approx(PSF_TO_PA, rel=1e-9)
    assert m["peak_psf"] == pytest.approx(1.0, rel=1e-9)
    # linear-ramp shocks (front 0 -> dp, rear -dp -> 0, both in tr):
    # peak / max slope = tr
    assert m["rise_time_ms"] == pytest.approx(1e3 * tr, rel=0.02)
    assert m["rise_time_10_90_ms"] == pytest.approx(0.8e3 * tr, rel=0.03)
    # shock mid-points: front at -tr/2, rear at T + tr/2
    assert m["duration_ms"] == pytest.approx(1e3 * (T + tr), abs=0.05)
    # positive phase: ramp (dp*tr/2) + triangle (dp*T/4)
    assert m["pos_impulse_pa_s"] == pytest.approx(PSF_TO_PA * (tr / 2 + T / 4), rel=1e-3)
    assert m["pldb"] == pytest.approx(pldb(t, p))
    assert 80 < m["csel"] < 130 and m["asel"] < m["csel"]


def test_metrics_front_shock_only_rise_time():
    """With a rear shock no steeper than the front (here: slow rear), the sBOOM
    rise time equals the front-shock ramp time."""
    fs, tr = 50e3, 3e-3
    t = np.arange(0, 1.0, 1 / fs)
    x = t - 0.3
    p = np.interp(x, [-tr, 0, 0.2, 0.2 + 0.05, 0.3], [0, 50, -50, 0, 0], left=0, right=0)
    m = metrics(t, p)
    assert m["rise_time_ms"] == pytest.approx(3.0, rel=0.02)
    assert m["rise_time_10_90_ms"] == pytest.approx(2.4, rel=0.03)


def test_no_warning_for_clean_signal():
    with warnings.catch_warnings():
        warnings.simplefilter("error")
        metrics(*nwave(1.0, 0.3, 2e-3))
