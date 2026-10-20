"""Sonic-boom loudness and waveform metrics.

Stevens Mark VII Perceived Level (PL, in PLdB) of a transient pressure
signature, computed with the procedure NASA uses for sonic booms
(Shepherd & Sullivan 1991), plus simple waveform metrics.

Procedure (``pldb``)
--------------------
1. Energy spectral density (ESD) of the signal, one-sided,
   ``E(f) = 2 |P(f)|^2`` with ``P(f) = dt * rfft(p)``, so that
   ``integral E(f) df = integral p(t)^2 dt``.  The record is zero-padded to at
   least ``_MIN_RECORD_S`` seconds so that even the 1.25 Hz band (0.29 Hz wide)
   spans a few FFT bins; this makes the result independent of the caller's
   padding and sample rate (see tests).  No window is applied by default; a
   record that does not start and end at ambient pressure produces a spurious
   step (PyLdB tapers both ends with a half-Hann window for that reason; the
   same taper is available through ``taper_s``).
2. Band energies ``E_i`` in the 42 one-third-octave bands with nominal centres
   1 Hz ... 12.5 kHz (exact base-10 band edges ``10^((n +- 0.5)/10)``), by
   trapezoidal integration of the ESD with the band edges inserted.
3. Band levels ``L_i = 10 log10(E_i / (0.07 s * p_ref^2)) - 3`` with
   ``p_ref = 20 uPa``: the band energy is treated as if delivered within the
   0.07 s integration time of the ear (Johnson & Robinson 1969; Shepherd &
   Sullivan 1991).  The ``-3 dB`` term (equivalently ``E/(2 * 0.07 * p_ref^2)``)
   follows PyLdB, which was benchmarked against NASA's LCASB loudness code
   (Bolander & Hunsaker 2019); it is kept so that results line up with that
   lineage of codes.
4. Mark VII equal-loudness contours: each band level is transformed into the
   level of the 3150 Hz reference band that is equally loud (``L_eq``), using
   Stevens' piecewise-straight contours in the closed form of Jackson &
   Leventhall (1973), as transcribed in PyLdB.  Below 80 Hz the contours are
   straight lines in (log f, L) converging on 160 dB at 1 Hz, so the 1 Hz band
   only contributes above 160 dB.
5. ``L_eq`` -> perceived magnitude ``S_i`` [sones] from Stevens' Mark VII
   table (Stevens 1972, Table I; 1 sone = 32 dB, doubling every 9 dB at
   moderate levels), linear interpolation.
6. Total ``S_t = S_m + F (sum S_i - S_m)`` with ``S_m = max S_i`` and the Mark
   VII one-third-octave summation factor ``F(S_m)`` (Stevens 1972 table of F
   versus loudest-band magnitude), linear interpolation.
7. ``PL = 32 + 9 log2(S_t)``.

Table provenance
----------------
All three tables below (``_SONES`` vs L_eq = 1..140 dB, ``_SUM_F`` vs
``_SONES[9:104]`` i.e. 0.181..348 sones, and the contour constants in
``_equivalent_level``) were copied digit-for-digit from the open-source PyLdB
(MIT licence), Utah State University AeroLab, Bolander & Hunsaker:
  https://github.com/usuaero/PyLdB  (file ``pyldb/core.py`` at commit 0e40b88,
  module globals ``SONES_L``, ``SUMMATION_FACTORS`` and functions
  ``_equivalent_loudness`` / ``_loud_limits_400``).
PyLdB cites as their origin:
  S. S. Stevens, "Perceived level of noise by Mark VII and decibels (E),"
  J. Acoust. Soc. Am. 51(2B), 575-601 (1972)  [sone & F tables, contours];
  G. M. Jackson & H. G. Leventhall, "Calculation of the perceived level of
  noise (PLdB) using Stevens' method (Mark VII)," Applied Acoustics 6(1),
  23-34 (1973)  [closed-form contours];
  K. P. Shepherd & B. M. Sullivan, "A loudness calculation procedure applied
  to shaped sonic booms," NASA TP-3134 (1991),
  https://ntrs.nasa.gov/citations/19920002547  [application to booms].
The NASA and journal originals could not be retrieved from this environment,
so the tables were not independently re-checked against them; the sone table
was checked to follow 2^((L-32)/9) to table rounding (< 1 %) from 13 dB up.
The whole chain reproduces PyLdB to ~0.01 dB (see tests), and PyLdB's authors
report agreement with NASA LCASB within 0.2 PLdB.
PyLdB's own nominal band edges (rounded ANSI values) were replaced by exact
base-10 edges here (differences < 0.5 % in frequency).

Waveform metrics (``metrics``)
------------------------------
See the ``metrics`` docstring for the definitions of rise time, duration and
positive impulse.
"""

from __future__ import annotations

import warnings

import numpy as np
from scipy.fft import next_fast_len

__all__ = ["pldb", "metrics", "band_levels", "PSF_TO_PA", "P_REF"]

PSF_TO_PA = 47.880259  # 1 lbf/ft^2 in Pa
P_REF = 20e-6  # Pa
T_EAR = 0.07  # s, Mark VII / Johnson-Robinson integration time
_MIN_RECORD_S = 10.0  # s, minimum zero-padded FFT record length (PL converged
                      # to < 0.001 dB already at 4 s for a 300 ms N-wave)

# --- One-third-octave bands: nominal centres 1 Hz .. 12.5 kHz (n = 0..41) ---
_BAND_N = np.arange(42)
BAND_CENTERS = np.array([
    1.0, 1.25, 1.6, 2.0, 2.5, 3.15, 4.0, 5.0, 6.3, 8.0,
    10.0, 12.5, 16.0, 20.0, 25.0, 31.5, 40.0, 50.0, 63.0, 80.0,
    100.0, 125.0, 160.0, 200.0, 250.0, 315.0, 400.0, 500.0, 630.0, 800.0,
    1000.0, 1250.0, 1600.0, 2000.0, 2500.0, 3150.0, 4000.0, 5000.0, 6300.0,
    8000.0, 10000.0, 12500.0])
_BAND_LO = 10.0 ** ((_BAND_N - 0.5) / 10.0)
_BAND_HI = 10.0 ** ((_BAND_N + 0.5) / 10.0)

# --- Stevens Mark VII perceived magnitude [sones] vs L_eq = 1..140 dB -------
# Source: PyLdB (usuaero/PyLdB @0e40b88, SONES_L), after Stevens (1972).
_LEQ_TABLE = np.arange(1.0, 141.0)
_SONES = np.array([
    0.078, 0.087, 0.097, 0.107, 0.118, 0.129, 0.141, 0.153, 0.166, 0.181,
    0.196, 0.212, 0.230, 0.248, 0.269, 0.290, 0.314, 0.339, 0.367, 0.396,
    0.428, 0.463, 0.500, 0.540, 0.583, 0.630, 0.680, 0.735, 0.794, 0.857,
    0.926, 1.000, 1.080, 1.170, 1.260, 1.360, 1.470, 1.590, 1.710, 1.850,
    2.000, 2.160, 2.330, 2.520, 2.720, 2.940, 3.180, 3.430, 3.700, 4.000,
    4.320, 4.670, 5.040, 5.440, 5.880, 6.350, 6.860, 7.410, 8.000, 8.640,
    9.330, 10.10, 10.90, 11.80, 12.70, 13.70, 14.80, 16.00, 17.30, 18.70,
    20.20, 21.80, 23.50, 25.40, 27.40, 29.60, 32.00, 34.60, 37.30, 40.30,
    43.50, 47.00, 50.80, 54.90, 59.30, 64.00, 69.10, 74.70, 80.60, 87.10,
    94.10, 102.0, 110.0, 119.0, 128.0, 138.0, 149.0, 161.0, 174.0, 188.0,
    203.0, 219.0, 237.0, 256.0, 276.0, 299.0, 323.0, 348.0, 376.0, 406.0,
    439.0, 474.0, 512.0, 553.0, 597.0, 645.0, 697.0, 752.0, 813.0, 878.0,
    948.0, 1024., 1106., 1194., 1290., 1393., 1505., 1625., 1756., 1896.,
    2048., 2212., 2389., 2580., 2787., 3010., 3251., 3511., 3792., 4096.])

# --- Mark VII one-third-octave summation factor F vs loudest band [sones] ---
# Abscissa is _SONES[9:104] (0.181 .. 348 sones, i.e. L_eq = 10 .. 104 dB).
# Source: PyLdB (usuaero/PyLdB @0e40b88, SUMMATION_FACTORS), after Stevens (1972).
_SUM_F_SONES = _SONES[9:104]
_SUM_F = np.array([
    0.100, 0.122, 0.140, 0.158, 0.174, 0.187, 0.200, 0.212, 0.222, 0.232,
    0.241, 0.250, 0.259, 0.267, 0.274, 0.281, 0.287, 0.293, 0.298, 0.303,
    0.308, 0.312, 0.316, 0.319, 0.320, 0.322, 0.322, 0.320, 0.319, 0.317,
    0.314, 0.311, 0.308, 0.304, 0.300, 0.296, 0.292, 0.288, 0.284, 0.279,
    0.275, 0.270, 0.266, 0.262, 0.258, 0.253, 0.248, 0.244, 0.240, 0.235,
    0.230, 0.226, 0.222, 0.217, 0.212, 0.208, 0.204, 0.200, 0.197, 0.195,
    0.194, 0.193, 0.192, 0.191, 0.190, 0.190, 0.190, 0.190, 0.190, 0.190,
    0.191, 0.191, 0.192, 0.193, 0.194, 0.195, 0.197, 0.199, 0.201, 0.203,
    0.205, 0.208, 0.210, 0.212, 0.215, 0.217, 0.219, 0.221, 0.223, 0.224,
    0.225, 0.226, 0.227, 0.227, 0.227])
assert _SUM_F.size == _SUM_F_SONES.size == 95


# ---------------------------------------------------------------------------
# Spectral helpers
# ---------------------------------------------------------------------------
def _check_signal(t, p):
    t = np.asarray(t, dtype=float)
    p = np.asarray(p, dtype=float)
    if t.ndim != 1 or p.ndim != 1 or t.size != p.size:
        raise ValueError("t and p must be 1-D arrays of equal length")
    if t.size < 2:
        raise ValueError("need at least two samples")
    dt = (t[-1] - t[0]) / (t.size - 1)
    if not dt > 0:
        raise ValueError("t must be increasing")
    return t, p, dt


def _esd(p, dt, min_record_s=_MIN_RECORD_S):
    """One-sided energy spectral density [Pa^2 s / Hz] on a zero-padded grid."""
    nfft = next_fast_len(max(p.size, int(np.ceil(min_record_s / dt))), real=True)
    P = np.fft.rfft(p, n=nfft) * dt
    esd = 2.0 * np.abs(P) ** 2
    esd[0] *= 0.5
    if nfft % 2 == 0:
        esd[-1] *= 0.5
    f = np.fft.rfftfreq(nfft, d=dt)
    return f, esd


def _band_energies(f, esd, lo, hi):
    """Trapezoidal integral of ``esd`` over each [lo, hi] with edges inserted."""
    df = np.diff(f)
    cum = np.concatenate(([0.0], np.cumsum(0.5 * (esd[1:] + esd[:-1]) * df)))

    def cum_at(fe):
        fe = np.clip(fe, f[0], f[-1])
        i = np.clip(np.searchsorted(f, fe, side="right") - 1, 0, f.size - 2)
        s_e = np.interp(fe, f, esd)
        return cum[i] + 0.5 * (esd[i] + s_e) * (fe - f[i])

    return np.maximum(cum_at(hi) - cum_at(lo), 0.0)


def _taper(p, dt, taper_s):
    """Half-Hann taper over ``taper_s`` seconds at each end (PyLdB ``window``)."""
    n = int(round(taper_s / dt)) if taper_s else 0
    if n <= 0:
        if p.size and max(abs(p[0]), abs(p[-1])) > 1e-3 * np.max(np.abs(p)):
            warnings.warn("signal does not start/end at ambient; the end step "
                          "adds spurious high-frequency energy (use taper_s)",
                          RuntimeWarning, stacklevel=3)
        return p
    if 2 * n > p.size:
        raise ValueError("taper longer than half the record")
    w = np.hanning(2 * n)
    p = p.copy()
    p[:n] *= w[:n]
    p[-n:] *= w[n:]
    return p


def band_levels(t, p, taper_s=0.0):
    """One-third-octave band levels used by Mark VII (step 3 of the module doc).

    Returns ``(fc, L)``: nominal centre frequencies [Hz] (1 Hz ... 12.5 kHz) and
    band levels [dB re 20 uPa, 0.07 s convention].  Bands above fs/2 get -inf.
    """
    t, p, dt = _check_signal(t, p)
    p = _taper(p, dt, taper_s)
    f, esd = _esd(p, dt)
    e = _band_energies(f, esd, _BAND_LO, _BAND_HI)
    e[_BAND_LO >= f[-1]] = 0.0
    with np.errstate(divide="ignore"):
        L = 10.0 * np.log10(e / (T_EAR * P_REF**2)) - 3.0
    return BAND_CENTERS.copy(), L


# ---------------------------------------------------------------------------
# Mark VII
# ---------------------------------------------------------------------------
def _limits_400(fc, lower, upper, L, x):
    """Transform a level at band ``fc`` (100-400 Hz region) to 3150 Hz equivalent.

    Stevens' contours between 1 Hz and 400 Hz are straight in (log f, L):
    below ``lower`` they converge on 115 dB at 1 Hz, above ``upper`` on 160 dB
    at 1 Hz, and in between they are parallel (offset ``x`` dB relative to
    400 Hz).  The final ``-8`` maps 400 Hz to the 3150 Hz reference.
    Constants from PyLdB ``_loud_limits_400`` (after Jackson & Leventhall 1973).
    """
    lg = np.log10(400.0) / np.log10(fc)
    return np.where(
        L <= lower, 115.0 - (115.0 - L) * lg,
        np.where(L <= upper, L - x, 160.0 - (160.0 - L) * lg)) - 8.0


def _equivalent_level(L):
    """Band levels (42 bands, 1 Hz..12.5 kHz) -> equally-loud 3150 Hz levels."""
    L = np.asarray(L, dtype=float)
    Leq = np.full(L.shape, -np.inf)
    i = _BAND_N
    fc = BAND_CENTERS
    fin = np.isfinite(L)

    m = (i >= 40) & fin                     # 10, 12.5 kHz
    Leq[m] = L[m] + 4.0 * (39 - i[m])
    m = (i >= 35) & (i <= 39) & fin         # 3150 .. 8000 Hz: reference region
    Leq[m] = L[m]
    m = (i >= 32) & (i <= 34) & fin         # 1600 .. 2500 Hz
    Leq[m] = L[m] - 2.0 * (35 - i[m])
    m = (i >= 27) & (i <= 31) & fin         # 500 .. 1250 Hz
    Leq[m] = L[m] - 8.0
    m = (i >= 20) & (i <= 26) & fin         # 100 .. 400 Hz
    x = 1.5 * (26 - i[m])
    Leq[m] = _limits_400(fc[m], 76.0 + x, 121.0 + x, L[m], x)
    m = (i >= 1) & (i <= 19) & fin          # 1.25 .. 80 Hz: map to 80 Hz first
    L80 = 160.0 - (160.0 - L[m]) * np.log10(80.0) / np.log10(fc[m])
    Leq[m] = _limits_400(80.0, 86.5, 131.5, L80, 10.5)
    # 1 Hz band (i == 0): contours converge on 160 dB at 1 Hz, so a 1 Hz band
    # below 160 dB is inaudible (L_eq -> -inf) and one above it is infinitely
    # loud; clamp the latter to the top of the table.
    if fin[0] and L[0] > 160.0:
        Leq[0] = _LEQ_TABLE[-1]
    return Leq


def _total_sones(Leq):
    S = np.interp(Leq, _LEQ_TABLE, _SONES, left=0.0, right=_SONES[-1])
    S[~np.isfinite(Leq)] = 0.0
    Sm = S.max()
    F = np.interp(Sm, _SUM_F_SONES, _SUM_F, left=_SUM_F[0], right=_SUM_F[-1])
    return Sm + F * (S.sum() - Sm), S


def pldb(t, p, taper_s=0.0):
    """Stevens Mark VII Perceived Level [PLdB] of a pressure history.

    Parameters
    ----------
    t : array, uniformly sampled time [s].
    p : array, overpressure [Pa].  The record should start and end at ambient
        (p ~ 0); it is zero-padded internally, so any padding the caller adds
        does not matter.  A RuntimeWarning is issued if it does not.
    taper_s : float, optional.  Length [s] of a half-Hann taper applied to each
        end of the record before analysis (0 = none).  PyLdB's default is 800
        samples.

    Returns
    -------
    float, PL in PLdB (``-inf`` for an all-zero signal).
    """
    _, L = band_levels(t, p, taper_s)
    St, _ = _total_sones(_equivalent_level(L))
    if St <= 0:
        return float("-inf")
    return float(32.0 + 9.0 * np.log2(St))


# ---------------------------------------------------------------------------
# Frequency weightings (IEC 61672) and exposure levels
# ---------------------------------------------------------------------------
def _weight_c2(f):
    f2 = np.asarray(f, float) ** 2
    rc = (12194.0**2 * f2) / ((f2 + 20.6**2) * (f2 + 12194.0**2))
    return rc**2


def _weight_a2(f):
    f2 = np.asarray(f, float) ** 2
    ra = (12194.0**2 * f2**2) / ((f2 + 20.6**2)
                                 * np.sqrt((f2 + 107.7**2) * (f2 + 737.9**2))
                                 * (f2 + 12194.0**2))
    return ra**2


_A1K2 = 1.0 / _weight_a2(1000.0)  # normalise to 0 dB at 1 kHz
_C1K2 = 1.0 / _weight_c2(1000.0)


def _sel(f, esd, w2):
    e = np.trapezoid(esd * w2, f) if hasattr(np, "trapezoid") else np.trapz(esd * w2, f)
    return float(10.0 * np.log10(e / (P_REF**2 * 1.0))) if e > 0 else float("-inf")


# ---------------------------------------------------------------------------
# Waveform metrics
# ---------------------------------------------------------------------------
def _crossing(t, p, i, level):
    """Linear-interpolated time where p crosses ``level`` between i and i+1."""
    p0, p1 = p[i], p[i + 1]
    if p1 == p0:
        return t[i]
    return t[i] + (level - p0) / (p1 - p0) * (t[i + 1] - t[i])


def metrics(t, p, taper_s=0.0):
    """Waveform metrics of a boom signature.

    Returns a dict with

    ``peak_pa``, ``peak_psf``
        Maximum overpressure ``max(p)``.
    ``rise_time_ms``
        sBOOM convention: ``max(p) / max(dp/dt)`` with the slope taken by
        first differences over the whole record (for a boom the steepest
        compression is the front or rear shock).
    ``rise_time_10_90_ms``
        Time for the front shock to rise from 10 % to 90 % of ``max(p)``:
        from the first upward crossing of 90 %, back to the last preceding
        sample below 10 % (linear interpolation at both ends).
    ``duration_ms``
        Shock-to-shock duration: from the front-shock mid-point (first upward
        crossing of ``max(p)/2``) to the rear-shock mid-point (last upward
        crossing of ``min(p)/2``).  For an N-wave this is the classical
        duration irrespective of rise time.  NaN if there is no negative phase.
    ``pos_impulse_pa_s``
        Integral of ``p`` over the first positive phase (zero crossing before
        the front shock to the next zero crossing).
    ``pldb``
        Stevens Mark VII PL, see :func:`pldb`.
    ``asel``, ``csel``
        A- and C-weighted sound exposure levels [dB re (20 uPa)^2 s]
        (IEC 61672 weightings applied to the energy spectral density).

    ``taper_s`` is passed to :func:`pldb` (and used for the SELs).
    """
    t, p, dt = _check_signal(t, p)
    out = {}
    pk = float(np.max(p))
    out["peak_pa"] = pk
    out["peak_psf"] = pk / PSF_TO_PA

    slope = np.diff(p) / dt
    smax = float(np.max(slope)) if slope.size else 0.0
    out["rise_time_ms"] = 1e3 * pk / smax if (pk > 0 and smax > 0) else float("nan")

    nan = float("nan")
    rt1090 = dur = imp = nan
    if pk > 0:
        # front shock: first upward crossing of 90 %
        i90 = int(np.argmax(p >= 0.9 * pk))
        below10 = np.nonzero(p[:i90] < 0.1 * pk)[0]
        if i90 > 0 and below10.size:
            j = below10[-1]
            rt1090 = 1e3 * (_crossing(t, p, i90 - 1, 0.9 * pk) - _crossing(t, p, j, 0.1 * pk))
        # duration: front-shock mid-point to rear-shock mid-point
        ih = int(np.argmax(p >= 0.5 * pk))
        t_front = _crossing(t, p, ih - 1, 0.5 * pk) if ih > 0 else t[0]
        pn = float(np.min(p))
        if pn < 0:
            up = np.nonzero((p[:-1] <= 0.5 * pn) & (p[1:] > 0.5 * pn))[0]
            if up.size and t[up[-1]] > t_front:
                dur = 1e3 * (_crossing(t, p, up[-1], 0.5 * pn) - t_front)
        # positive impulse of the first positive phase
        nonpos = np.nonzero(p[:ih] <= 0)[0]
        i0 = nonpos[-1] if nonpos.size else 0
        after = np.nonzero(p[ih:] <= 0)[0]
        i1 = ih + after[0] if after.size else p.size - 1
        seg = np.clip(p[i0:i1 + 1], 0.0, None)
        imp = float(np.sum(0.5 * (seg[1:] + seg[:-1])) * dt)
    out["rise_time_10_90_ms"] = rt1090
    out["duration_ms"] = dur
    out["pos_impulse_pa_s"] = imp
    out["pldb"] = pldb(t, p, taper_s)

    f, esd = _esd(_taper(p, dt, taper_s) if taper_s else p, dt)
    out["asel"] = _sel(f, esd, _weight_a2(f) * _A1K2)
    out["csel"] = _sel(f, esd, _weight_c2(f) * _C1K2)
    return out
