"""Aircraft source model: equivalent area, Whitham F-function, near-field signature.

Public API
----------
Aircraft(name, length, mach, altitude, weight, volume, lift_start=0.4, lift_end=1.0)
equivalent_area(ac, x, phi=0)       A_e(x, phi) = A_v(x) + A_L(x) cos(phi)      [m^2]
equivalent_area_slope(ac, x, phi=0) dA_e/dx (analytic)                          [m]
f_function(ac, y, phi=0)            Whitham F(y) = (1/2pi) int_0^y A_e''(xi)/sqrt(y-xi) dxi  [m^1/2]
near_field(ac, r0, fs, phi=0)       (t, p): dp = p_inf gamma M^2 F(y)/sqrt(2 beta r0), t = y/V
CONCORDE, SR71               reference aircraft used for validation (see docstrings)

Model
-----
* Volume: Sears-Haack body, A_v = (16 Vol/(3 pi L)) [4 xi (1 - xi)]^{3/2}, xi = x/L.
* Lift: A_L(x) = beta/(rho V^2) int_0^x l dx with a sin^2 lift distribution on
  [lift_start L, lift_end L] integrating to the weight (l and l' vanish at both ends).
* Azimuth: phi is the downward azimuth around the flight axis, measured from straight
  down (phi = 0 undertrack, +-90 deg horizontal).  In linearised supersonic theory the
  volume is an axisymmetric line of sources and the lift a line of vertical dipoles,
  whose field carries cos(phi); so the equivalent area seen at azimuth phi is
  A_e(x, phi) = A_v(x) + A_L(x) cos(phi) (e.g. Plotkin 2002, Int. J. Acoust. Vib. 7;
  Carlson NASA TP-1122 lift term).  phi = 0 reproduces the undertrack values bit for
  bit.
* F-function: A_e'' is treated as piecewise constant on a uniform grid of spacing h
  (cell value from the *analytic* slope, (A'(x_{i+1}) - A'(x_i))/h, so the integrable
  x^{-1/2} nose singularity of the Sears-Haack A'' is handled exactly in the mean),
  and the kernel is integrated analytically per cell:
  int_cell dxi/sqrt(y - xi) = 2 [sqrt(y - x_i) - sqrt(y - x_{i+1})].
  On the matching y-grid this is a discrete convolution, evaluated by FFT.
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from scipy.signal import fftconvolve

from .atmosphere import GAMMA, std_atm


@dataclass
class Aircraft:
    name: str
    length: float          # m
    mach: float
    altitude: float        # m
    weight: float          # N (supported by lift)
    volume: float          # m^3 (Sears-Haack distributed)
    lift_start: float = 0.4   # fraction of length
    lift_end: float = 1.0

    @property
    def atm(self):
        return std_atm(self.altitude)

    @property
    def V(self) -> float:
        return float(self.mach * self.atm.c)

    @property
    def beta(self) -> float:
        return float(np.sqrt(self.mach ** 2 - 1.0))


# --- reference aircraft (validation) -------------------------------------------
# Concorde: overall length 61.66 m; cruise M 2.0 at FL520-FL600 (NASA FS-016 quotes
# 1.94 psf at M 2, 52 000 ft).  Mid-cruise mass ~150 t -> W = 1.47 MN.  Volume 450 m^3:
# fuselage ~ 2.9 m diameter x ~50 m effective (~300 m^3) + wing 358 m^2 x ~0.4 m mean
# thickness x ~0.7 fill (~100 m^3) + nacelles (~50 m^3).  Lift on the delta wing,
# root chord ~0.33 L ... ~0.9 L.
CONCORDE = Aircraft("Concorde", length=61.66, mach=2.0, altitude=16800.0,
                    weight=1.47e6, volume=450.0, lift_start=0.35, lift_end=0.90)
# SR-71: length 32.74 m; cruise M 3.0 at ~24.4 km (80 000 ft).  MTOW 78 t, empty
# ~30.6 t: mid-mission mass ~55 t -> W = 0.54 MN.  Volume 150 m^3: fuel ~ 36 t of JP-7
# (~45 m^3) plus fuselage/chines/nacelles (2 nacelles ~ 1.3 m^2 x 12 m; fuselage
# + chines ~ 3 m^2 x 25 m).  Lift on chine + delta wing, ~0.25 L ... 0.95 L.
SR71 = Aircraft("SR-71", length=32.74, mach=3.0, altitude=24400.0,
                weight=0.54e6, volume=150.0, lift_start=0.25, lift_end=0.95)


# --- equivalent area ----------------------------------------------------------
def _lift_scale(ac: Aircraft) -> float:
    a = ac.atm
    return ac.beta / (float(a.rho) * ac.V ** 2)


def lift_distribution(ac: Aircraft, x):
    """l(x) [N/m]: sin^2 bump on [lift_start L, lift_end L] integrating to the weight."""
    x = np.asarray(x, dtype=float)
    a, b = ac.lift_start * ac.length, ac.lift_end * ac.length
    u = (x - a) / (b - a)
    inside = (u > 0) & (u < 1)
    return np.where(inside, 2.0 * ac.weight / (b - a) * np.sin(np.pi * np.clip(u, 0, 1)) ** 2, 0.0)


def _lift_cumulative(ac: Aircraft, x):
    x = np.asarray(x, dtype=float)
    a, b = ac.lift_start * ac.length, ac.lift_end * ac.length
    u = np.clip((x - a) / (b - a), 0.0, 1.0)
    # int_0^u 2 sin^2(pi u') du' = u - sin(2 pi u)/(2 pi)
    return ac.weight * (u - np.sin(2 * np.pi * u) / (2 * np.pi))


def _sears_haack_smax(ac: Aircraft) -> float:
    return 16.0 * ac.volume / (3.0 * np.pi * ac.length)


def volume_area(ac: Aircraft, x):
    xi = np.clip(np.asarray(x, dtype=float) / ac.length, 0.0, 1.0)
    return _sears_haack_smax(ac) * (4.0 * xi * (1.0 - xi)) ** 1.5


def lift_area(ac: Aircraft, x):
    return _lift_scale(ac) * _lift_cumulative(ac, x)


def equivalent_area(ac: Aircraft, x, phi: float = 0.0):
    """A_e(x, phi) = A_v(x) + A_L(x) cos(phi) [m^2]; behind the body the lift part is
    constant (= beta W cos(phi)/(rho V^2))."""
    if phi == 0.0:
        return volume_area(ac, x) + lift_area(ac, x)
    return volume_area(ac, x) + np.cos(phi) * lift_area(ac, x)


def equivalent_area_slope(ac: Aircraft, x, phi: float = 0.0):
    """Analytic dA_e/dx [m] at downward azimuth phi [rad]."""
    x = np.asarray(x, dtype=float)
    xi = x / ac.length
    inside = (xi > 0) & (xi < 1)
    xc = np.clip(xi, 0.0, 1.0)
    dAv = np.where(inside, _sears_haack_smax(ac) * 1.5 * np.sqrt(4 * xc * (1 - xc))
                   * 4.0 * (1 - 2 * xc) / ac.length, 0.0)
    return dAv + (_lift_scale(ac) * float(np.cos(phi))) * lift_distribution(ac, x)


# --- F-function -----------------------------------------------------------------
def _f_on_grid(ac: Aircraft, h: float, n_y: int, oversample: int = 1, phi: float = 0.0):
    """F at y_m = m h, m = 0..n_y-1 (F(0) = 0 by the left-limit convention).

    The cell-constant A'' approximation has an O(m^-3/2) relative error at the m-th
    sample behind the Sears-Haack nose singularity (27 % at m = 1, 0.6 % at m = 8);
    ``oversample`` evaluates on an h/oversample grid and decimates.
    """
    if oversample > 1:
        return _f_on_grid(ac, h / oversample, (n_y - 1) * oversample + 1, phi=phi)[::oversample]
    n_body = int(np.ceil(ac.length / h))
    x_nodes = np.arange(n_body + 1) * h
    dA = equivalent_area_slope(ac, np.minimum(x_nodes, ac.length), phi)
    A2 = np.diff(dA) / h                            # piecewise-constant A'' per cell
    j = np.arange(1, n_y)                           # m - i for cells fully left of y_m
    w = 2.0 * np.sqrt(h) * (np.sqrt(j) - np.sqrt(j - 1))
    conv = fftconvolve(A2, w)[: n_y - 1]            # sum_i A2_i w_{m-i}, m = 1..n_y-1
    F = np.zeros(n_y)
    F[1:] = conv / (2.0 * np.pi)
    return F


def f_function(ac: Aircraft, y, h=None, phi: float = 0.0):
    """Whitham F(y) [m^1/2] at arbitrary y (F = 0 for y <= 0), downward azimuth phi [rad].

    Computed on a uniform grid of spacing ``h`` (default L/20000) and linearly
    interpolated.
    """
    y = np.asarray(y, dtype=float)
    h = ac.length / 20000.0 if h is None else h
    y_max = max(float(np.max(y)), h) if y.size else h
    n_y = int(np.ceil(y_max / h)) + 2
    F = _f_on_grid(ac, h, n_y, oversample=4, phi=phi)
    return np.interp(y, np.arange(n_y) * h, F, left=0.0)


# --- near field ---------------------------------------------------------------------
def near_field(ac: Aircraft, r0: float, fs: float, y_max_factor: float = 3.0,
               pad_before: float | None = None, pad_after: float | None = None,
               phi: float = 0.0):
    """Linear Whitham near-field signature at radial distance r0 from the flight axis
    at downward azimuth ``phi`` [rad] (0 = straight below; uses A_e(x, phi)).

    dp = p_inf gamma M^2 F(y)/sqrt(2 beta r0), t = y/V (retarded time; t = 0 at the
    nose characteristic).  F is computed directly on the sample grid (dy = V/fs) up to
    y = y_max_factor * L (A_e'(L) = 0 so the tail decays ~ y^-3/2) and the record is
    zero padded by ``pad_before``/``pad_after`` seconds (default 2 L/V each).
    """
    a = ac.atm
    V = ac.V
    h = V / fs
    n_y = int(np.ceil(y_max_factor * ac.length / h)) + 1
    F = _f_on_grid(ac, h, n_y, oversample=8, phi=phi)
    p_sig = float(a.p) * GAMMA * ac.mach ** 2 * F / np.sqrt(2.0 * ac.beta * r0)
    tL = ac.length / V
    nb = int(round((2 * tL if pad_before is None else pad_before) * fs))
    na = int(round((2 * tL if pad_after is None else pad_after) * fs))
    p = np.concatenate([np.zeros(nb), p_sig, np.zeros(na)])
    t = (np.arange(p.size) - nb) / fs
    return t, p
