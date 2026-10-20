"""Undertrack ray geometry and augmented-Burgers propagation of a sonic boom.

Public API
----------
RayPath(altitude, V, atmosphere=std_atm, humidity=None)
    Undertrack ray from an aircraft in steady level flight (no wind).
    .atm(z), .nz(z) (vertical direction cosine), .s(z) (arc length from the aircraft),
    .area(z) (ray-tube area up to a constant), .humidity(z), .amplitude_ratio(z1, z2).
propagate(t, p, ray, z_from, z_to, **opts) -> (t, p)
ground_signature(ac, r0=None, insert=None, fs=50e3, reflection=1.9, **opts) -> (t, p)
iso_absorption(f, T, p, h) -> alpha [Np/m] (ISO 9613-1), relaxation_params(T, p, h)
propagation_filter(f, T, p, c, h, ds, ...) -> complex H(f) for one linear step
propagation_exponent(...) -> ln H(f) (attenuation + dispersive phase)
uniform_atmosphere(z_ref=0.0) -> callable giving a height-independent Atm (tests)

Geometry (no wind, stratified)
------------------------------
Horizontal slowness s_x = 1/V is conserved; s_z(z) = sqrt(1/c^2 - 1/V^2) (requires
V > c everywhere, i.e. no cut-off).  The ray's vertical direction cosine is
n_z = c s_z, so ds = dz/(c s_z).  The ray-tube area (normal to the ray) is the product
of the along-track width (rays launched at successive times: dx * n_z) and the lateral
width from an azimuth increment (s_y = dphi * s_z(z_a) conserved, so
dy = dphi s_z(z_a) int_z^{z_a} dz'/s_z):  A(z) ~ c s_z * int_z^{z_a} dz'/s_z.
In a uniform atmosphere this is A ~ (z_a - z) ~ distance (cylindrical spreading).
Linear amplitude invariant (Blokhintzev, no wind): p sqrt(A/(rho c)) = const.

The near field from :func:`boomsim.source.near_field` is computed at radial distance
r0 from the flight axis; for the undertrack ray that is exactly the depth below the
aircraft, so the march starts at z = z_a - r0 (the near-field theory itself assumes a
uniform atmosphere over r < r0).

Augmented Burgers equation (Cleveland 1995; Rallabhandi 2011, sBOOM)
-------------------------------------------------------------------
In retarded time tau = t - int ds/c0 (c0 = local small-signal/equilibrium sound
speed), per arc length s, operator-split into:
 (a) nonlinearity  dp/ds = beta_n/(rho c^3) p dp/dtau, beta_n = 1.2.  Solved in
     conservative form dp/ds + d/dtau(-k p^2/2) = 0 with a finite-volume scheme:
     MUSCL-Hancock (van Leer limiter, second order, one flux evaluation per step),
     Engquist-Osher flux (entropy-satisfying for this convex flux).  Shocks are
     captured without multivaluedness and obey the equal-area rule; the step is set
     adaptively by the CFL condition k max|p| ds/dt <= cfl (default 0.8).  Only the
     active part of the record (|p| > 1e-10 max|p|) is updated.  (SSP-RK2 variant
     ``burgers_step_rk2`` kept for cross-checks.)
 (b) thermoviscous absorption  H = exp(-a_cl f^2 ds), a_cl = 1.84e-11 (p_r/p)(T/T0)^1/2.
 (c) O2 and N2 relaxation with causal dispersion (ISO 9613-1 frequencies/strengths):
     H_v = exp(-(w/c)(eps_v/2) w tau_v/(1 + j w tau_v) ds), tau_v = 1/(2 pi f_rv),
     eps_v = A_v c/pi, which gives |H_v| = exp(-A_v f^2 f_rv/(f_rv^2 + f^2) ds) and a
     phase speed that rises with frequency (in the rfft e^{-jwt} convention a positive
     phase exp(+jw d) is an advance).
 (d) spreading + stratification: p *= sqrt(A_old/A_new * (rho c)_new/(rho c)_old),
     applied after every nonlinear sub-step.
(b)+(c) are applied together in the frequency domain once per "linear chunk" of arc
length <= ds_lin (default 20 m) using the atmosphere at the chunk mid-point; the
nonlinear sub-steps inside a chunk are as small as the CFL condition requires.

Sampling: default fs = 50 kHz.  Ground rise times are set by N2/O2 relaxation and
are ~1-5 ms, i.e. 50-250 samples; at altitude the physical shock is far thinner than
any affordable grid and the captured shock is ~2-3 samples wide (20-60 us), which
does not matter downstream because the ground rise time is a local balance of
steepening and relaxation that is re-established over the last few km.  The record is
zero padded to an FFT-friendly length and extended automatically if the waveform
approaches either end.
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from scipy import fft as sfft
from scipy.integrate import cumulative_trapezoid

from .atmosphere import P_REF, Atm, molar_humidity, std_atm
from .source import near_field

T0_ISO = 293.15      # ISO 9613-1 reference temperature
BETA_AIR = 1.2       # coefficient of nonlinearity (gamma + 1)/2


# ------------------------------------------------------------------------------
# atmosphere helpers
# ------------------------------------------------------------------------------
def uniform_atmosphere(z_ref: float = 0.0):
    """Return a callable z -> Atm with the std_atm(z_ref) state at every height."""
    a0 = std_atm(z_ref)

    def atm(z):
        z = np.asarray(z, dtype=float)
        one = np.ones_like(z)
        return Atm(z=z, T=a0.T * one, p=a0.p * one, rho=a0.rho * one, c=a0.c * one)

    return atm


# ------------------------------------------------------------------------------
# absorption (ISO 9613-1)
# ------------------------------------------------------------------------------
def relaxation_params(T, p, h):
    """ISO 9613-1 parameters.

    Returns dict with f_rO, f_rN [Hz], A_O, A_N [Np/m] (alpha_v = A_v f^2 f_rv/(f_rv^2+f^2)),
    and a_cl [Np s^2/m] (alpha_cl = a_cl f^2).
    """
    T = np.asarray(T, dtype=float)
    pr = np.asarray(p, dtype=float) / P_REF
    h = np.asarray(h, dtype=float)
    Tr = T / T0_ISO
    f_rO = pr * (24.0 + 4.04e4 * h * (0.02 + h) / (0.391 + h))
    f_rN = pr * Tr ** -0.5 * (9.0 + 280.0 * h * np.exp(-4.170 * (Tr ** (-1.0 / 3.0) - 1.0)))
    A_O = Tr ** -2.5 * 0.01275 * np.exp(-2239.1 / T)
    A_N = Tr ** -2.5 * 0.1068 * np.exp(-3352.0 / T)
    a_cl = 1.84e-11 / pr * Tr ** 0.5
    return dict(f_rO=f_rO, f_rN=f_rN, A_O=A_O, A_N=A_N, a_cl=a_cl)


def iso_absorption(f, T, p, h):
    """ISO 9613-1 pure-tone absorption coefficient alpha(f) [Np/m] (x 8.686 for dB/m)."""
    f = np.asarray(f, dtype=float)
    r = relaxation_params(T, p, h)
    f2 = f * f
    return (r["a_cl"] * f2
            + r["A_O"] * f2 * r["f_rO"] / (r["f_rO"] ** 2 + f2)
            + r["A_N"] * f2 * r["f_rN"] / (r["f_rN"] ** 2 + f2))


def propagation_exponent(f, T, p, c, h, ds, thermoviscous=True, relaxation=True,
                         dispersion=True):
    """ln H(f) for the linear (b)+(c) operators over arc length ds (see
    :func:`propagation_filter`); real part = -alpha ds, imaginary part = phase advance."""
    f = np.asarray(f, dtype=float)
    w = 2.0 * np.pi * f
    r = relaxation_params(T, p, h)
    lnH = np.zeros(f.shape, dtype=complex)
    if thermoviscous:
        lnH -= r["a_cl"] * f * f * ds
    if relaxation:
        for A, fr in ((r["A_O"], r["f_rO"]), (r["A_N"], r["f_rN"])):
            tau = 1.0 / (2.0 * np.pi * fr)
            eps = A * c / np.pi
            wt = w * tau
            term = (w / c) * (eps / 2.0) * wt / (1.0 + 1j * wt) * ds
            lnH -= term if dispersion else term.real
    return lnH


def propagation_filter(f, T, p, c, h, ds, thermoviscous=True, relaxation=True,
                       dispersion=True):
    """Complex transfer function over arc length ds for the linear (b)+(c) operators.

    rfft convention (forward kernel e^{-jwt}); retarded time w.r.t. the low-frequency
    sound speed c.  With dispersion=False the relaxation term is applied as pure
    (zero-phase) attenuation.
    """
    return np.exp(propagation_exponent(f, T, p, c, h, ds, thermoviscous=thermoviscous,
                                       relaxation=relaxation, dispersion=dispersion))


# ------------------------------------------------------------------------------
# ray geometry
# ------------------------------------------------------------------------------
class RayPath:
    """Undertrack ray from an aircraft at ``altitude`` [m] flying level at ``V`` [m/s].

    Quantities are tabulated on a fine z grid (``dz`` m) from ``z_min`` to the
    aircraft altitude and linearly interpolated.  ``atmosphere`` is a callable z -> Atm
    (default std_atm); ``humidity`` a callable z -> h [%] (default molar_humidity
    evaluated with the given atmosphere).
    """

    def __init__(self, altitude, V, atmosphere=None, humidity=None, dz=2.0, z_min=0.0):
        self.altitude = float(altitude)
        self.V = float(V)
        self._atm_fn = std_atm if atmosphere is None else atmosphere
        n = int(np.ceil((self.altitude - z_min) / dz)) + 1
        z = np.linspace(z_min, self.altitude, n)
        a = self._atm_fn(z)
        self.z_tab = z
        self.T_tab = np.asarray(a.T, dtype=float) * np.ones_like(z)
        self.p_tab = np.asarray(a.p, dtype=float) * np.ones_like(z)
        self.rho_tab = np.asarray(a.rho, dtype=float) * np.ones_like(z)
        self.c_tab = np.asarray(a.c, dtype=float) * np.ones_like(z)
        if np.any(self.c_tab >= self.V):
            raise ValueError("ray cut-off: V <= c(z) somewhere on the path")
        if humidity is None:
            self.h_tab = molar_humidity(z, atm=a) * np.ones_like(z)
        else:
            self.h_tab = np.asarray(humidity(z), dtype=float) * np.ones_like(z)
        sz = np.sqrt(1.0 / self.c_tab ** 2 - 1.0 / self.V ** 2)
        self.sz_tab = sz
        self.nz_tab = self.c_tab * sz
        # integrals from z up to the aircraft altitude
        self.s_tab = (cumulative_trapezoid(1.0 / self.nz_tab[::-1], z[::-1], initial=0.0))[::-1] * -1.0
        self.lat_tab = (cumulative_trapezoid(1.0 / sz[::-1], z[::-1], initial=0.0))[::-1] * -1.0
        self.rhoc_tab = self.rho_tab * self.c_tab

    def _i(self, tab, z):
        return np.interp(z, self.z_tab, tab)

    def atm(self, z) -> Atm:
        z = np.asarray(z, dtype=float)
        return Atm(z=z, T=self._i(self.T_tab, z), p=self._i(self.p_tab, z),
                   rho=self._i(self.rho_tab, z), c=self._i(self.c_tab, z))

    def nz(self, z):
        """Vertical direction cosine of the ray, c s_z."""
        return self._i(self.nz_tab, z)

    def s(self, z):
        """Arc length along the ray from the aircraft down to altitude z [m]."""
        return self._i(self.s_tab, z)

    def area(self, z):
        """Ray-tube area up to a constant: c s_z int_z^{z_a} dz'/s_z [m s / m ...]."""
        return self._i(self.nz_tab, z) * self._i(self.lat_tab, z)

    def humidity(self, z):
        return self._i(self.h_tab, z)

    def amplitude_ratio(self, z1, z2):
        """Linear lossless p(z2)/p(z1) = sqrt(A1/A2 * (rho c)_2/(rho c)_1)."""
        return np.sqrt(self.area(z1) / self.area(z2)
                       * self._i(self.rhoc_tab, z2) / self._i(self.rhoc_tab, z1))


# ------------------------------------------------------------------------------
# nonlinear step: conservative finite volume for dp/ds + d/dtau(-k p^2/2) = 0
# ------------------------------------------------------------------------------
def _slopes(u):
    """Cell values with two edge-extended ghost cells each side, and van Leer slopes
    for cells 1..N+2 of that padded array."""
    ue = np.empty(u.size + 4)
    ue[2:-2] = u
    ue[:2] = u[0]
    ue[-2:] = u[-1]
    d = np.diff(ue)                     # N+3
    d1, d2 = d[:-1], d[1:]
    prod = d1 * d2
    s = np.zeros_like(prod)
    np.divide(2.0 * prod, d1 + d2, out=s, where=prod > 0.0)
    return ue[1:-1], s                  # N+2 cells (one ghost each side)


def _eo_flux(uL, uR, k):
    """Engquist-Osher flux for g(u) = -k u^2/2 (entropy-satisfying)."""
    return -0.5 * k * (np.minimum(uL, 0.0) ** 2 + np.maximum(uR, 0.0) ** 2)


def _flux_divergence(u, k):
    """(F_{i+1/2} - F_{i-1/2}) with MUSCL/van Leer reconstruction + EO flux."""
    uc, s = _slopes(u)
    F = _eo_flux(uc[:-1] + 0.5 * s[:-1], uc[1:] - 0.5 * s[1:], k)
    return F[1:] - F[:-1]


def burgers_step_rk2(p, k, dt, ds):
    """SSP-RK2 MUSCL step (TVD for k max|p| ds/dt <= 0.5)."""
    lam = ds / dt
    p1 = p - lam * _flux_divergence(p, k)
    return 0.5 * (p + p1 - lam * _flux_divergence(p1, k))


def burgers_step(p, k, dt, ds):
    """Advance dp/ds = k p dp/dtau by ds with MUSCL-Hancock (2nd order, one flux
    evaluation; stable for CFL k max|p| ds/dt <= 1, default used 0.8)."""
    lam = ds / dt
    uc, s = _slopes(p)
    half = (0.5 * lam * k) * s * uc     # -(lam/2)(g(u+) - g(u-)) for g = -k u^2/2
    um = uc - 0.5 * s + half
    up = uc + 0.5 * s + half
    F = _eo_flux(up[:-1], um[1:], k)
    return p - lam * (F[1:] - F[:-1])


def _active_window(p, rel=1e-10, guard=8):
    a = np.abs(p)
    idx = np.flatnonzero(a > rel * a.max())
    if idx.size == 0:
        return 0, 0
    return max(0, idx[0] - guard), min(p.size, idx[-1] + guard + 1)


# ------------------------------------------------------------------------------
# propagation
# ------------------------------------------------------------------------------
def _extend(t, p, dt, front, back):
    n0 = p.size
    target = sfft.next_fast_len(n0 + front + back, real=True)
    back += target - (n0 + front + back)
    p = np.concatenate([np.zeros(front), p, np.zeros(back)])
    t = t[0] + (np.arange(p.size) - front) * dt
    return t, p


def propagate(t, p, ray: RayPath, z_from, z_to, linear=False, lossless=False,
              thermoviscous=True, relaxation=True, dispersion=True, beta_n=BETA_AIR,
              cfl=0.8, ds_lin=20.0, auto_extend=True, edge_frac=0.03, edge_tol=1e-3,
              stats=None):
    """Augmented-Burgers propagation along ``ray`` from z_from down to z_to (< z_from).

    Options: ``linear`` (skip nonlinearity), ``lossless`` (skip absorption/dispersion),
    ``thermoviscous``/``relaxation``/``dispersion`` (individual linear terms), ``cfl``
    (nonlinear CFL number), ``ds_lin`` (max arc length per linear chunk [m]),
    ``auto_extend`` (zero-pad when the signal reaches the record ends).  If ``stats``
    is a dict it receives step counts.  Returns (t, p) (t may be longer than input).
    """
    if z_to > z_from:
        raise ValueError("propagate marches downward: need z_to <= z_from")
    t = np.asarray(t, dtype=float)
    p = np.array(p, dtype=float)
    dt = float(t[1] - t[0])
    if auto_extend:
        t, p = _extend(t, p, dt, 0, 0)
    n_nl = n_lin = 0
    z = float(z_from)
    eps_z = 1e-9 * max(1.0, abs(z_from))
    while z > z_to + eps_z:
        dz_chunk = min(ds_lin * float(ray.nz(z)), z - z_to)
        z_end = z - dz_chunk
        if z_end - z_to < 1e-6 * ds_lin:
            z_end = z_to
        # --- nonlinear sub-steps + spreading
        zc = z
        while zc > z_end + eps_z:
            dz = zc - z_end
            if not linear:
                pmax = float(np.max(np.abs(p)))
                if pmax > 0.0:
                    a = ray.atm(zc)
                    k = beta_n / (float(a.rho) * float(a.c) ** 3)
                    dz = min(dz, cfl * dt / (k * pmax) * float(ray.nz(zc)))
                    zm = zc - 0.5 * dz
                    am = ray.atm(zm)
                    k = beta_n / (float(am.rho) * float(am.c) ** 3)
                    i0, i1 = _active_window(p)
                    p[i0:i1] = burgers_step(p[i0:i1], k, dt, dz / float(ray.nz(zm)))
                    n_nl += 1
            p *= float(ray.amplitude_ratio(zc, zc - dz))
            zc -= dz
        # --- linear absorption + dispersion over the chunk
        if not lossless and (thermoviscous or relaxation):
            zm = 0.5 * (z + z_end)
            am = ray.atm(zm)
            ds = float(ray.s(z_end) - ray.s(z))
            f = sfft.rfftfreq(p.size, dt)
            H = propagation_filter(f, am.T, am.p, am.c, ray.humidity(zm), ds,
                                   thermoviscous=thermoviscous, relaxation=relaxation,
                                   dispersion=dispersion)
            p = sfft.irfft(sfft.rfft(p) * H, n=p.size)
            n_lin += 1
        z = z_end
        # --- keep the waveform away from the record ends
        if auto_extend:
            m = max(4, int(edge_frac * p.size))
            thr = edge_tol * np.max(np.abs(p))
            front = int(0.25 * p.size) if np.max(np.abs(p[:m])) > thr else 0
            back = int(0.25 * p.size) if np.max(np.abs(p[-m:])) > thr else 0
            if front or back:
                t, p = _extend(t, p, dt, front, back)
    if stats is not None:
        stats.update(nonlinear_steps=n_nl, linear_steps=n_lin, n_samples=p.size)
    return t, p


def ground_signature(ac, r0=None, insert=None, fs=50e3, reflection=1.9, z_ground=0.0,
                     atmosphere=None, humidity=None, below_opts=None, **opts):
    """Near field at r0 below the aircraft -> augmented Burgers to the ground.

    r0 defaults to 1.5 body lengths.  ``insert=(z_ins, fn)``: propagate to z_ins,
    apply ``p = fn(t, p)`` there (fn may also return a new ``(t, p)``), continue to
    the ground.  The ground reflection factor is applied last.  Extra keyword options
    go to :func:`propagate`; ``below_opts`` (dict, needs ``insert``) overrides them
    for the leg below z_ins only, e.g. ``below_opts={"linear": True}`` to switch
    nonlinearity off after the insert (re-steepening diagnostics).
    """
    r0 = 1.5 * ac.length if r0 is None else float(r0)
    t, p = near_field(ac, r0, fs)
    ray = RayPath(ac.altitude, ac.V, atmosphere=atmosphere, humidity=humidity,
                  z_min=min(0.0, z_ground))
    z0 = ac.altitude - r0
    if insert is not None:
        z_ins, fn = insert
        if not (z_ground <= z_ins <= z0):
            raise ValueError(f"insert altitude {z_ins} outside [{z_ground}, {z0}]")
        t, p = propagate(t, p, ray, z0, z_ins, **opts)
        out = fn(t, p)
        if isinstance(out, tuple):
            t, p = out
        else:
            p = np.asarray(out, dtype=float)
        t, p = propagate(t, p, ray, z_ins, z_ground, **{**opts, **(below_opts or {})})
    else:
        if below_opts:
            raise ValueError("below_opts requires insert")
        t, p = propagate(t, p, ray, z0, z_ground, **opts)
    return t, reflection * p
