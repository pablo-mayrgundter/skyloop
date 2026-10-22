"""Boundary-layer turbulence acting on a ground boom signature: random phase screen +
Fresnel-Kirchhoff integral, evaluated for all frequencies at once (delay-area method).

Physics
-------
Turbulence in the planetary boundary layer (PBL, lowest ~1 km) perturbs the effective sound
speed by mu = -(dT/(2 T0) + w'/c0) (temperature and along-path, here vertical, velocity
fluctuations; Ostashev & Wilson, Acoustics in Moving Inhomogeneous Media, 2nd ed. 2015).
Following the thin-phase-screen picture used for sonic-boom/N-wave turbulence studies
(Plotkin & George, JFM 54 (1972) 449; Pierce, JASA 44 (1968) 1052 "spikes";
Lipkens & Blackstock, JASA 103 (1998) 148; Locey & Sparrow, Noise Control Eng. J. 55
(2007) 314 "turbulence filters"), the layer is collapsed into one screen at distance ``d``
above the receiver that imposes a travel-time perturbation dT(xi) (xi = 2-D screen
coordinate), the line integral of -mu/c0 across a layer of thickness ``layer``.
The 2-D spectrum of dT for a von Karman index spectrum

    Phi_mu(k) = sigma_mu^2 Gamma(11/6)/(pi^{3/2} Gamma(1/3)) L0^3 / (1 + k^2 L0^2)^{11/6}

is F_T(k_perp) = 2 pi layer Phi_mu(k_perp, 0)/c0^2 (Markov approximation).

Below the screen the wave travels linearly (Fresnel-Kirchhoff, plane incident wave):

    P_out(f) = P_in(f) (j w/(2 pi c d)) iint w(xi) exp(-j w tau(xi)) d^2 xi,
    tau(xi) = dT(xi) + |xi|^2/(2 c d).

The double integral is the Fourier transform of the delay-area density g(t) = dA/dt,
A(t) = area of {tau < t} (weighted), so the impulse response is

    h(t) = A''(t) / (2 pi c d)       (no turbulence: A = 2 pi c d t -> h = delta).

A(t) is computed EXACTLY for tau linearly interpolated on a triangulated screen grid (each
triangle's area-vs-delay is piecewise quadratic), accumulated with kink sums on the output
time grid: cost O(N_triangles + N_t), all frequencies, caustics (spikes) and geometric
spreading of arrivals (rounding) included, diffraction resolved down to the grid scale.
Focusing gives spikes (dA/dt has caustic peaks), defocusing and multipath round the shock.

Aperture.  Turbulence is applied inside radius ``r_turb`` and blended to the screen mean
outside it (cos^2 blend over ``r_blend``); the Kirchhoff aperture is tapered to zero
between ``r_ap_in`` and ``r_ap``.  The same aperture without turbulence gives h0, and the
returned kernel is delta + h - h0: the artificial aperture edge cancels exactly and
dT = 0 returns the identity.  Frequencies whose first Fresnel zone sqrt(lambda d) exceeds
r_turb (below ~10-20 Hz) are left essentially untouched, as physically expected (the
turbulence is weak at those scales anyway).

Limitations: one screen (no amplitude scintillation building up inside the layer),
plane-wave normal incidence (oblique rays: use slant ``d``/``layer``), linear below the
screen (see apply_screen_then_propagate in run_wind_turbulence.py for the nonlinear
bracket), frozen turbulence, ground reflection applied as a factor (both paths see the
same screen).

Public API
----------
TurbulenceParams(sigma_mu2, L0, layer, d, c, dx, n, r_turb, r_blend, r_ap_in, r_ap, name)
PRESETS: dict name -> TurbulenceParams ("weak", "moderate", "strong")
delay_sigma(params) -> analytic rms of the screen travel-time perturbation [s]
von_karman_screen(params, rng) -> (xi, dT) square grid [m], [s]
screen_kernel(dT, params, fs) -> (k, i0): discrete kernel, p_out = conv(p, k) aligned at i0
apply_kernel(p, k, i0) -> p_out (same length as p)
realisation(params, fs, rng) -> (k, i0, info) one random kernel
"""
from __future__ import annotations

from dataclasses import dataclass, replace

import numpy as np
from scipy.signal import fftconvolve
from scipy.special import gamma


@dataclass(frozen=True)
class TurbulenceParams:
    sigma_mu2: float = 1.0e-5      # variance of the effective index fluctuation mu
    L0: float = 150.0              # von Karman length scale (outer scale) [m]
    layer: float = 1000.0          # turbulent path length through the PBL [m]
    d: float = 500.0               # screen -> receiver distance [m]
    c: float = 340.0               # sound speed in the layer [m/s]
    dx: float = 0.5                # screen grid spacing [m]
    n: int = 1024                  # grid points per side (box n*dx; periodic FFT screen)
    r_turb: float = 120.0          # turbulence applied inside this radius ...
    r_blend: float = 40.0          # ... and blended to the screen mean over this width
    r_ap_in: float = 170.0         # Kirchhoff aperture taper start [m]
    r_ap: float = 230.0            # aperture radius [m]
    name: str = "custom"

    @property
    def sigma_mu(self) -> float:
        return float(np.sqrt(self.sigma_mu2))


def _mu2(sigma_T, T0, sigma_w, c0=340.0):
    return (sigma_T / (2.0 * T0)) ** 2 + (sigma_w / c0) ** 2


# Representative PBL states (temperature std sigma_T [K], vertical-velocity std sigma_w
# [m/s]): stable/light-wind, typical daytime (flight-test "moderate"), strongly convective.
PRESETS = {
    "weak": TurbulenceParams(sigma_mu2=_mu2(0.2, 288.0, 0.3), L0=100.0, name="weak"),
    "moderate": TurbulenceParams(sigma_mu2=_mu2(0.5, 288.0, 1.0), L0=150.0, name="moderate"),
    "strong": TurbulenceParams(sigma_mu2=_mu2(1.0, 288.0, 1.7), L0=200.0, name="strong"),
}


def _phi_mu(k, p: TurbulenceParams):
    a = gamma(11.0 / 6.0) / (np.pi ** 1.5 * gamma(1.0 / 3.0))
    return p.sigma_mu2 * a * p.L0 ** 3 / (1.0 + (k * p.L0) ** 2) ** (11.0 / 6.0)


def _F_T(k, p: TurbulenceParams):
    """2-D spectrum of the screen travel-time perturbation [s^2 m^2]."""
    return 2.0 * np.pi * p.layer * _phi_mu(k, p) / p.c ** 2


def delay_sigma(p: TurbulenceParams) -> float:
    """rms travel-time perturbation of an infinite screen: sqrt(int F_T d^2k) [s]
    (closed form: sigma_mu^2 layer L0 * 2 sqrt(pi) Gamma(5/6)/Gamma(1/3) / c^2)."""
    var = p.sigma_mu2 * p.layer * p.L0 * 2.0 * np.sqrt(np.pi) * gamma(5.0 / 6.0) / gamma(1.0 / 3.0)
    return float(np.sqrt(var) / p.c)


def von_karman_screen(p: TurbulenceParams, rng):
    """Random screen dT on an n x n periodic grid (spectral synthesis). Returns (x, dT)
    with x the 1-D grid coordinate centred on 0 [m] and dT [s] (zero mean)."""
    n, dx = int(p.n), float(p.dx)
    kx = 2.0 * np.pi * np.fft.fftfreq(n, dx)
    K = np.hypot(kx[:, None], kx[None, :])
    dk = 2.0 * np.pi / (n * dx)
    amp = np.sqrt(_F_T(K, p)) * dk
    amp[0, 0] = 0.0
    z = rng.standard_normal((n, n)) + 1j * rng.standard_normal((n, n))
    dT = np.real(np.fft.fft2(amp * z))
    x = (np.arange(n) - n // 2) * dx
    return x, np.fft.fftshift(dT)


def _weights(R, p: TurbulenceParams):
    """Kirchhoff aperture weight and turbulence blend at radius R."""
    w_ap = np.where(R <= p.r_ap_in, 1.0,
                    np.where(R >= p.r_ap, 0.0,
                             np.cos(0.5 * np.pi * (R - p.r_ap_in) / (p.r_ap - p.r_ap_in)) ** 2))
    b = np.clip((R - p.r_turb) / p.r_blend, 0.0, 1.0)
    w_tu = np.cos(0.5 * np.pi * b) ** 2
    return w_ap, w_tu


def _area_curve(tau, wgt, x, dt, t0, nt):
    """A(t_j), t_j = t0 + j dt (j = 0..nt-1): weighted area of {tau < t} for tau linear on
    the triangles of the grid (two per cell), exactly (piecewise quadratic per triangle)."""
    dx = float(x[1] - x[0])
    # triangles: (i,j),(i+1,j),(i,j+1) and (i+1,j+1),(i,j+1),(i+1,j)
    a, b_, c_, d_ = tau[:-1, :-1], tau[1:, :-1], tau[:-1, 1:], tau[1:, 1:]
    wc = 0.25 * (wgt[:-1, :-1] + wgt[1:, :-1] + wgt[:-1, 1:] + wgt[1:, 1:])
    keep = wc > 0.0
    T = np.concatenate([np.stack([a[keep], b_[keep], c_[keep]], 1),
                        np.stack([d_[keep], c_[keep], b_[keep]], 1)])
    W = np.concatenate([wc[keep], wc[keep]]) * 0.5 * dx * dx      # triangle area x weight
    T.sort(axis=1)
    t1, t2, t3 = T[:, 0], T[:, 1], T[:, 2]
    eps = 1e-4 * dt                      # vertex values closer than this are separated
    t2 = np.maximum(t2, t1 + eps)
    t3 = np.maximum(t3, t2 + eps)
    s1 = 2.0 * W / ((t2 - t1) * (t3 - t1))
    s2 = 2.0 * W / ((t3 - t1) * (t3 - t2))
    # kinks of A'' (piecewise constant): +s1 at t1, -(s1+s2) at t2, +s2 at t3
    kt = np.concatenate([t1, t2, t3]) - t0
    kc = np.concatenate([s1, -(s1 + s2), s2])
    m = np.floor(kt / dt).astype(np.int64)
    if m.min() < 0 or m.max() >= nt - 1:
        raise ValueError("delay range exceeds the time grid")
    e = (m + 1) * dt - kt                # in (0, dt]
    C0 = np.bincount(m, kc, minlength=nt)
    C1 = np.bincount(m, kc * e, minlength=nt)
    C2 = np.bincount(m, kc * e * e, minlength=nt)
    mm = np.arange(nt, dtype=float)
    P0, P1, P2 = np.cumsum(C0), np.cumsum(C0 * mm), np.cumsum(C0 * mm * mm)
    Q0, Q1 = np.cumsum(C1), np.cumsum(C1 * mm)
    R0 = np.cumsum(C2)
    # A_j = sum_{m <= j-1} [C0 (n dt)^2/2 + C1 n dt + C2/2], n = j-1-m
    A = np.zeros(nt)
    jm = mm[1:] - 1.0
    A[1:] = (0.5 * dt * dt * (jm * jm * P0[:-1] - 2.0 * jm * P1[:-1] + P2[:-1])
             + dt * (jm * Q0[:-1] - Q1[:-1]) + 0.5 * R0[:-1])
    return A


def screen_kernel(dT, p: TurbulenceParams, fs: float, x=None):
    """Discrete turbulence kernel for screen ``dT`` (n x n, centred grid).

    Returns (k, i0): p_out[n] = sum_j k[j] p_in[n - j + i0], i.e. ``apply_kernel``.
    k = delta + h(dT) - h(0) with the same aperture (exact identity for dT = 0).
    """
    n = dT.shape[0]
    if x is None:
        x = (np.arange(n) - n // 2) * p.dx
    X, Y = np.meshgrid(x, x, indexing="ij")
    R = np.hypot(X, Y)
    w_ap, w_tu = _weights(R, p)
    inner = R <= p.r_turb
    mean_in = float(np.mean(dT[inner]))
    dTb = mean_in + (dT - mean_in) * w_tu
    geo = (X * X + Y * Y) / (2.0 * p.c * p.d)
    dt = 1.0 / fs
    tau = dTb - mean_in + geo
    near = R <= p.r_ap + 2.0 * p.dx
    lo = min(float(tau[near].min()), 0.0)
    hi = max(float(tau[near].max()), float(geo[near].max()))
    t0 = (np.floor(lo / dt) - 3) * dt
    nt = int(np.ceil((hi - t0) / dt)) + 6
    A = _area_curve(tau, w_ap, x, dt, t0, nt)
    A0 = _area_curve(geo, w_ap, x, dt, t0, nt)
    norm = 2.0 * np.pi * p.c * p.d * dt
    k = (np.diff(A, 2) - np.diff(A0, 2)) / norm          # k[j] at time t0 + (j+1) dt
    i0 = int(round(-t0 / dt)) - 1                         # index of t = 0
    k[i0] += 1.0
    return k, i0


def apply_kernel(p, k, i0):
    """Convolve signature samples ``p`` with kernel (k, i0); output aligned with input."""
    p = np.asarray(p, dtype=float)
    out = fftconvolve(p, k)
    return out[i0:i0 + p.size]


def realisation(p: TurbulenceParams, fs: float, rng):
    """One random kernel; info has the screen rms inside r_turb and the stationary delay."""
    x, dT = von_karman_screen(p, rng)
    k, i0 = screen_kernel(dT, p, fs, x=x)
    X, Y = np.meshgrid(x, x, indexing="ij")
    inner = np.hypot(X, Y) <= p.r_turb
    return k, i0, dict(dT_rms_inner=float(np.std(dT[inner])), dT_rms=float(np.std(dT)))


def with_(p: TurbulenceParams, **kw) -> TurbulenceParams:
    return replace(p, **kw)
