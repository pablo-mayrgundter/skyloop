"""Ray geometry (undertrack and off-track) and augmented-Burgers propagation of a sonic boom.

Public API
----------
RayPath(altitude, V, atmosphere=std_atm, humidity=None, phi=0.0, allow_cutoff=None)
    Ray from an aircraft in steady level flight (no wind) at downward azimuth phi.
    .atm(z), .nz(z) (vertical direction cosine), .s(z) (arc length from the aircraft),
    .area(z) (ray-tube area up to a constant), .humidity(z), .amplitude_ratio(z1, z2),
    .x(z), .y(z), .travel_time(z), .dy_dphi(z), .arrival(z), .y_ground, .cutoff, .z_turn.
lateral_cutoff(altitude, V, atmosphere=None, z_ground=0.0) -> (phi_c, y_c)
propagate(t, p, ray, z_from, z_to, **opts) -> (t, p)
ground_signature(ac, r0=None, insert=None, fs=50e3, reflection=1.9, phi=0.0, **opts) -> (t, p)
propagate_from_circle(ac, phi, t, p, R_m, z_ref=20000.0, anchor="source", age_from=None,
                      **opts) -> (t_g, p_g, y_g)
    start from a uniform-air matching-circle signature (3-D pipeline), WKB-corrected.
circle_start_factors(ac, phi, R_m, ...) -> dict (the linear start factors alone)
ground_carpet(ac, phis, start="near_field" | "circle", procs=4, **opts) -> dict
    per-phi ground signatures, y_g, loudness (raw + Gate-0 rise-time proxy) and carpet
    metrics (undertrack PLdB, max PLdB, half-width).
iso_absorption(f, T, p, h) -> alpha [Np/m] (ISO 9613-1), relaxation_params(T, p, h)
propagation_filter(f, T, p, c, h, ds, ...) -> complex H(f) for one linear step
propagation_exponent(...) -> ln H(f) (attenuation + dispersive phase)
uniform_atmosphere(z_ref=0.0) -> callable giving a height-independent Atm (tests)

Geometry (no wind, stratified)
------------------------------
Off-track rays: see :class:`RayPath` (horizontal slowness (1/V, s_y) conserved,
s_z = sqrt(1/c^2 - 1/V^2 - s_y^2), lateral cut-off where s_z -> 0 above the ground).
Undertrack (phi = 0, the Gate-0 path, kept bit-identical):
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
from .proxy import METRIC_TAPER_S, RISE_SIGMA_S, scored_metrics
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
    """Ray from an aircraft at ``altitude`` [m] flying level along +x at ``V`` [m/s],
    launched at downward azimuth ``phi`` [rad] (0 = undertrack; |phi| < pi/2).

    Stratified, windless atmosphere.  Launch wave normal (sin mu0, cos mu0 sin phi,
    -cos mu0 cos phi), sin mu0 = c(z_a)/V; conserved horizontal slowness (1/V, s_y),
    s_y = cos mu0 sin phi / c(z_a); vertical slowness s_z = sqrt(1/c^2 - 1/V^2 - s_y^2);
    ray direction = wave normal (no wind), n_z = c s_z.

    Quantities are tabulated on a fine z grid (``dz`` m) from ``z_min`` (or the turning
    altitude, see below) to the aircraft altitude and linearly interpolated.
    ``atmosphere`` is a callable z -> Atm (default std_atm); ``humidity`` a callable
    z -> h [%] (default molar_humidity evaluated with the given atmosphere).

    Lateral cut-off: if s_z^2 <= 0 somewhere above ``z_min`` the ray turns at the
    highest such altitude ``z_turn`` and never reaches the ground; then ``cutoff`` is
    True, the tables stop at z_turn and ``y_ground`` is NaN.  ``allow_cutoff``
    (default: True for phi != 0, False for phi == 0, the legacy behaviour) chooses
    between flagging that and raising ValueError.

    Integrals from z up to z_a: I1 = int dz/s_z, I3 = int dz/s_z^3.  Ray position
    relative to the emission point: x = I1/V, y = s_y I1; travel time int dz/(c^2 s_z).
    Ray-tube area from the neighbouring rays (launch time t0 <-> x0 = V t0, azimuth
    phi): the horizontal footprint has Jacobian |d(X,Y)/d(x0,phi)| = dY/dphi (dX/dx0 = 1,
    dY/dx0 = 0 by steady flight) with dY/dphi = (cos mu0 cos phi/c_a) (I1 + s_y^2 I3);
    the area normal to the ray is that times |n_z|.  ``area`` returns it normalised by
    cos mu0/c_a, i.e. A = c s_z cos(phi) (I1 + s_y^2 I3), which tends to c_a r near the
    source (r = distance from the flight axis) and is c r everywhere in uniform air.
    At phi = 0 this is exactly the undertrack formula A = c s_z I1, and for phi = 0 the
    legacy trapezoidal tables are used unchanged (bit-identical Gate-0 results).  For
    phi != 0 each z-cell is integrated exactly for s_z^2 linear in z
    (int dz/s_z = 2 dz/(r_i + r_{i+1}), int dz/s_z^3 = 2 dz/(r_i r_{i+1} (r_i + r_{i+1})),
    r = s_z), which keeps the integrable 1/sqrt singularity of grazing rays accurate;
    this cell rule is the exact s_y-derivative of the I1 rule, so the analytic Jacobian
    equals the finite difference of neighbouring traced rays (tested).
    """

    def __init__(self, altitude, V, atmosphere=None, humidity=None, dz=2.0, z_min=0.0,
                 phi=0.0, allow_cutoff=None):
        self.altitude = float(altitude)
        self.V = float(V)
        self.phi = float(phi)
        if not abs(self.phi) < 0.5 * np.pi:
            raise ValueError("phi must satisfy |phi| < pi/2 (downward rays)")
        legacy = self.phi == 0.0
        if allow_cutoff is None:
            allow_cutoff = not legacy
        self._atm_fn = std_atm if atmosphere is None else atmosphere
        n = int(np.ceil((self.altitude - z_min) / dz)) + 1
        z = np.linspace(z_min, self.altitude, n)
        a = self._atm_fn(z)
        c = np.asarray(a.c, dtype=float) * np.ones_like(z)
        c_a = float(c[-1])
        if c_a >= self.V:
            raise ValueError("no Mach cone: V <= c at the aircraft altitude")
        self.c_a = c_a
        self.cos_mu0 = float(np.sqrt(1.0 - (c_a / self.V) ** 2))
        self.sy = self.cos_mu0 * float(np.sin(self.phi)) / c_a
        self.cutoff = False
        self.z_turn = float("nan")
        q = 1.0 / c ** 2 - 1.0 / self.V ** 2 - self.sy ** 2
        if np.any(q <= 0.0) if not legacy else np.any(c >= self.V):
            if not allow_cutoff:
                raise ValueError("ray cut-off: V <= c(z) somewhere on the path" if legacy
                                 else f"ray cut-off: turns above z_min (phi={self.phi:g})")
            i = int(np.flatnonzero(q <= 0.0)[-1])
            zt = float(z[i] + (z[i + 1] - z[i]) * (-q[i]) / (q[i + 1] - q[i]))
            self.cutoff = True
            self.z_turn = zt
            z = np.concatenate([[zt], z[i + 1:]])
            a = self._atm_fn(z)
            c = np.asarray(a.c, dtype=float) * np.ones_like(z)
            q = 1.0 / c ** 2 - 1.0 / self.V ** 2 - self.sy ** 2
            q[0] = 0.0
        self.z_tab = z
        self.T_tab = np.asarray(a.T, dtype=float) * np.ones_like(z)
        self.p_tab = np.asarray(a.p, dtype=float) * np.ones_like(z)
        self.rho_tab = np.asarray(a.rho, dtype=float) * np.ones_like(z)
        self.c_tab = c
        if humidity is None:
            self.h_tab = molar_humidity(z, atm=a) * np.ones_like(z)
        else:
            self.h_tab = np.asarray(humidity(z), dtype=float) * np.ones_like(z)
        legacy = legacy and not self.cutoff
        if legacy:
            sz = np.sqrt(1.0 / self.c_tab ** 2 - 1.0 / self.V ** 2)
        else:
            sz = np.sqrt(np.maximum(q, 0.0))
        self.sz_tab = sz
        self.nz_tab = self.c_tab * sz

        def cum(d):                         # integral from z_tab[k] up to z_a
            return np.concatenate([np.cumsum(d[::-1])[::-1], [0.0]])

        dzc = np.diff(z)
        rs = sz[:-1] + sz[1:]
        with np.errstate(divide="ignore", invalid="ignore"):
            dI1 = 2.0 * dzc / rs
            dI3 = 2.0 * dzc / (sz[:-1] * sz[1:] * rs)
        cbar = 0.5 * (c[:-1] + c[1:])
        self.tt_tab = cum(dI1 / cbar ** 2)           # travel time from the aircraft
        if legacy:
            # integrals from z up to the aircraft altitude (legacy trapezoid, unchanged)
            self.s_tab = (cumulative_trapezoid(1.0 / self.nz_tab[::-1], z[::-1], initial=0.0))[::-1] * -1.0
            self.lat_tab = (cumulative_trapezoid(1.0 / sz[::-1], z[::-1], initial=0.0))[::-1] * -1.0
            self.I1_tab = self.lat_tab
            self.I13_tab = self.lat_tab
        else:
            self.s_tab = cum(dI1 / cbar)
            self.I1_tab = cum(dI1)
            self.I13_tab = self.I1_tab + self.sy ** 2 * cum(dI3)
            self.lat_tab = float(np.cos(self.phi)) * self.I13_tab
        self.rhoc_tab = self.rho_tab * self.c_tab

    def _i(self, tab, z):
        return np.interp(z, self.z_tab, tab)

    @property
    def z_bottom(self) -> float:
        """Lowest altitude on the ray table (z_min, or z_turn for a cut-off ray)."""
        return float(self.z_tab[0])

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
        """Ray-tube area up to a constant: c s_z cos(phi) (I1 + s_y^2 I3) (-> c_a r at the
        source; = c s_z int_z^{z_a} dz'/s_z undertrack)."""
        return self._i(self.nz_tab, z) * self._i(self.lat_tab, z)

    def humidity(self, z):
        return self._i(self.h_tab, z)

    def amplitude_ratio(self, z1, z2):
        """Linear lossless p(z2)/p(z1) = sqrt(A1/A2 * (rho c)_2/(rho c)_1)."""
        return np.sqrt(self.area(z1) / self.area(z2)
                       * self._i(self.rhoc_tab, z2) / self._i(self.rhoc_tab, z1))

    # --- ray position / timing -------------------------------------------------
    def x(self, z):
        """Along-track offset of the ray point at altitude z from the emission point [m]."""
        return self._i(self.I1_tab, z) / self.V

    def y(self, z):
        """Lateral offset of the ray point at altitude z [m] (sign of phi)."""
        return self.sy * self._i(self.I1_tab, z)

    def travel_time(self, z):
        """Travel time from emission to altitude z [s]."""
        return self._i(self.tt_tab, z)

    def dy_dphi(self, z):
        """Lateral Jacobian dY/dphi [m/rad] of the ray family at altitude z."""
        return self.cos_mu0 * np.cos(self.phi) / self.c_a * self._i(self.I13_tab, z)

    def arrival(self, z=None):
        """Arrival of the ray at altitude z (default the table bottom z_min).

        Returns dict(x, y, z, T, delay): position relative to the emission point, travel
        time T, and ``delay`` = T - x/V, the time between the aircraft passing over the
        arrival point's x and the boom arriving.  NaN entries for a cut-off ray when z is
        below the turning altitude.
        """
        z = self.z_bottom if z is None else float(z)
        if self.cutoff and z < self.z_turn:
            nan = float("nan")
            return dict(x=nan, y=nan, z=z, T=nan, delay=nan)
        x = float(self.x(z))
        T = float(self.travel_time(z))
        return dict(x=x, y=float(self.y(z)), z=z, T=T, delay=T - x / self.V)

    @property
    def y_ground(self) -> float:
        """Lateral offset y_g where the ray reaches z_min (NaN if cut off)."""
        return float("nan") if self.cutoff else float(self.y(self.z_tab[0]))


def lateral_cutoff(altitude, V, atmosphere=None, z_ground=0.0, dz=2.0):
    """Lateral cut-off of the boom carpet (no wind).

    The ray with azimuth phi reaches the ground iff s_y^2 < min_{z_g<=z<=z_a}(1/c^2 - 1/V^2).
    Returns (phi_c [rad], y_c [m]): limiting azimuth and the ground lateral offset of
    the (grazing) limiting ray, i.e. the carpet half-width.  phi_c = pi/2, y_c = inf if
    no cut-off occurs below 90 deg.
    """
    atm = std_atm if atmosphere is None else atmosphere
    n = int(np.ceil((altitude - z_ground) / dz)) + 1
    z = np.linspace(z_ground, altitude, n)
    c = np.asarray(atm(z).c, dtype=float) * np.ones_like(z)
    c_a = float(c[-1])
    cos_mu0 = np.sqrt(1.0 - (c_a / V) ** 2)
    qmin = float(np.min(1.0 / c ** 2 - 1.0 / V ** 2))
    if qmin <= 0.0:
        return 0.0, 0.0
    sphi = c_a * np.sqrt(qmin) / cos_mu0
    if sphi >= 1.0:
        return 0.5 * np.pi, float("inf")
    phi_c = float(np.arcsin(sphi))
    for fac in (1.0 - 1e-10, 1.0 - 1e-7, 1.0 - 1e-5):
        ray = RayPath(altitude, V, atmosphere=atmosphere, dz=dz, z_min=z_ground,
                      phi=phi_c * fac, allow_cutoff=True)
        if not ray.cutoff:
            return phi_c, ray.y_ground
    return phi_c, float("nan")


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
    if ray.cutoff and z_to < ray.z_turn:
        raise ValueError(f"ray is cut off: turns at z = {ray.z_turn:.0f} m, above z_to = {z_to:g} m")
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
                     atmosphere=None, humidity=None, below_opts=None, phi=0.0, **opts):
    """Near field at r0 from the flight axis -> augmented Burgers to the ground.

    r0 defaults to 1.5 body lengths.  ``phi`` [rad] is the downward launch azimuth
    (0 = undertrack): the near field uses A_e(x, phi) and starts on the off-track ray at
    altitude z_a - r0 cos(phi) (uniform air assumed inside r0, as for the near-field
    theory itself); ValueError if that ray is laterally cut off.  ``insert=(z_ins, fn)``:
    propagate to z_ins, apply ``p = fn(t, p)`` there (fn may also return a new
    ``(t, p)``), continue to the ground.  The ground reflection factor is applied last.
    Extra keyword options go to :func:`propagate`; ``below_opts`` (dict, needs ``insert``)
    overrides them for the leg below z_ins only, e.g. ``below_opts={"linear": True}`` to
    switch nonlinearity off after the insert (re-steepening diagnostics).
    """
    r0 = 1.5 * ac.length if r0 is None else float(r0)
    phi = float(phi)
    t, p = near_field(ac, r0, fs, phi=phi)
    ray = RayPath(ac.altitude, ac.V, atmosphere=atmosphere, humidity=humidity,
                  z_min=min(0.0, z_ground), phi=phi)
    if ray.cutoff and ray.z_turn > z_ground:
        raise ValueError(f"phi = {phi:g} rad is beyond the lateral cut-off (turns at "
                         f"z = {ray.z_turn:.0f} m)")
    z0 = ac.altitude - r0 * float(np.cos(phi))
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


# ------------------------------------------------------------------------------
# start from a matching circle (3-D pipeline)
# ------------------------------------------------------------------------------
def circle_start_factors(ac, phi, R_m, ray=None, z_ref=20000.0, anchor="source",
                         atmosphere=None, humidity=None, z_ground=0.0):
    """Linear factors that map a uniform-air signature on the matching circle onto the
    actual stratified ray (see :func:`propagate_from_circle`).

    Returns dict(ray, z_start, geometric, impedance, factor, z_anchor).
    """
    phi = float(phi)
    if ray is None:
        ray = RayPath(ac.altitude, ac.V, atmosphere=atmosphere, humidity=humidity,
                      z_min=min(0.0, z_ground), phi=phi, allow_cutoff=True)
    z_s = ac.altitude - float(R_m) * float(np.cos(phi))
    if z_s <= ray.z_bottom:
        raise ValueError(f"matching circle start z = {z_s:g} m is below the ray table")
    if anchor == "source":
        z_anchor = ac.altitude
    elif anchor == "ref":
        z_anchor = float(z_ref)
    else:
        z_anchor = float(anchor)
    atm_fn = ray._atm_fn
    a_s, a_0 = atm_fn(np.array([z_s])), atm_fn(np.array([z_anchor]))
    imp = float(np.sqrt((a_s.rho * a_s.c) / (a_0.rho * a_0.c))[0])
    geo = float(np.sqrt(ray.c_a * float(R_m) / float(ray.area(z_s))))
    return dict(ray=ray, z_start=z_s, geometric=geo, impedance=imp, factor=geo * imp,
                z_anchor=z_anchor)


def propagate_from_circle(ac, phi, t, p, R_m, z_ref=20000.0, anchor="source",
                          age_from=None, reflection=1.9, z_ground=0.0, atmosphere=None,
                          humidity=None, info=None, **opts):
    """Propagate a matching-circle signature to the ground along the off-track ray.

    ``(t, p)`` is the signature at radius ``R_m`` [m] from the flight axis at downward
    azimuth ``phi`` [rad], computed (e.g. by the cross-plane solver) linearly in UNIFORM
    air with the ambient state at ``z_ref``.  Returns ``(t_g, p_g, y_g)``: the ground
    signature (reflection factor applied) and its lateral offset; ``(None, None, nan)``
    if the ray is laterally cut off.  ``info`` (dict) receives the start altitude and
    correction factors.  Extra keywords go to :func:`propagate`.

    Start point.  In uniform air the ray at azimuth phi stays in the half-plane of that
    azimuth through the flight axis and is straight, so the circle point is at slant
    depth R_m cos(phi) below the aircraft.  The march starts on the *actual* stratified
    ray at altitude z_s = z_a - R_m cos(phi).  Approximation: the actual ray's curvature
    inside R_m is ignored for the start position (inside 11-20 km c is constant, so the
    actual ray is straight there too; above 20 km dc/dz = +1 K/km gives a negligible
    bend over a few km).

    Amplitude (linear WKB / Blokhintzev invariant p sqrt(A/(rho c)) = const).  The
    uniform-air solution conserves p sqrt(A_u) with straight-ray tubes A_u = c_a r, the
    actual one p sqrt(A/(rho c)); both are anchored where they coincide.  Hence

        p_actual(z_s) = p_circle * sqrt(c_a R_m / A(z_s)) * sqrt((rho c)(z_s) / (rho c)(z_anchor)),

    with A the stratified ray-tube area (normalised so that A -> c_a r at the source).
    (rho c)_ref of the uniform medium cancels: only the anchor matters.  ``anchor``:
    "source" (default; z_anchor = aircraft altitude) is right when the circle field's
    source strength is Whitham's near field with the aircraft's ambient state
    (``source.near_field``, as the cross-plane solver's contract specifies); "ref"
    (z_anchor = z_ref) if the source was normalised with the z_ref state; or an altitude.
    The shield itself is linear, so the factor is the same with or without it.  As
    11-20 km is isothermal, c cancels there and the factor is sqrt(rho(z_s)/rho(z_a)),
    ~ exp((z_a - z_s)/(2 H)), H = R T/g = 6.34 km: +17 % (1.4 dB) for R_m = 2 km.

    Nonlinearity inside R_m.  The circle signature carries no nonlinear distortion
    (or absorption) from the near field out to R_m.  ``age_from=r_f`` [m] applies a cheap
    correction: the corrected circle signature is scaled back linearly (lossless WKB,
    actual atmosphere) to radius r_f on the same ray and the full model is marched from
    there.  For a simple wave this restores exactly the missing nonlinear age
    int k p ds from r_f to R_m (the lossless Burgers solution with spreading depends on
    the path only through that age).  Use r_f = 1.5 L (the near-field start of the
    direct calculation) for the part of the field that never met the shield, or the
    shield's slant distance for the transmitted part; default None = no correction.
    """
    phi = float(phi)
    ray = RayPath(ac.altitude, ac.V, atmosphere=atmosphere, humidity=humidity,
                  z_min=min(0.0, z_ground), phi=phi, allow_cutoff=True)
    if ray.cutoff and ray.z_turn > z_ground:
        if info is not None:
            info.update(cutoff=True, z_turn=ray.z_turn)
        return None, None, float("nan")
    fac = circle_start_factors(ac, phi, R_m, ray=ray, z_ref=z_ref, anchor=anchor)
    z_s = fac["z_start"]
    t = np.asarray(t, dtype=float)
    p = np.asarray(p, dtype=float) * fac["factor"]
    z_from = z_s
    if age_from is not None and float(age_from) < float(R_m):
        z_from = ac.altitude - float(age_from) * float(np.cos(phi))
        p = p / float(ray.amplitude_ratio(z_from, z_s))
    t, p = propagate(t, p, ray, z_from, z_ground, **opts)
    y_g = ray.arrival(z_ground)["y"]
    if info is not None:
        info.update(cutoff=False, z_start=z_s, z_march=z_from, geometric=fac["geometric"],
                    impedance=fac["impedance"], factor=fac["factor"],
                    z_anchor=fac["z_anchor"], y_g=y_g, arrival=ray.arrival(z_ground))
    return t, reflection * p, y_g


# ------------------------------------------------------------------------------
# ground carpet
# ------------------------------------------------------------------------------
def _carpet_job(job):
    phi = job["phi"]
    ac = job["ac"]
    zg = job["z_ground"]
    ray = RayPath(ac.altitude, ac.V, atmosphere=job["atmosphere"], humidity=job["humidity"],
                  z_min=min(0.0, zg), phi=phi, allow_cutoff=True)
    out = dict(phi=phi, cutoff=bool(ray.cutoff and ray.z_turn > zg), z_turn=ray.z_turn)
    arr = ray.arrival(zg)
    out.update(y_g=arr["y"], x_g=arr["x"], delay=arr["delay"])
    if out["cutoff"]:
        out.update(t=None, p=None, metrics=None)
        return out
    if job["start"] == "near_field":
        t, p = ground_signature(ac, r0=job["r0"], fs=job["fs"], reflection=job["reflection"],
                                z_ground=zg, atmosphere=job["atmosphere"],
                                humidity=job["humidity"], phi=phi, **job["opts"])
    else:
        tc, pc = job["circle"]
        t, p, _ = propagate_from_circle(ac, phi, tc, pc, job["R_m"], z_ref=job["z_ref"],
                                        anchor=job["anchor"], age_from=job["age_from"],
                                        reflection=job["reflection"], z_ground=zg,
                                        atmosphere=job["atmosphere"],
                                        humidity=job["humidity"], **job["opts"])
    out["metrics"] = scored_metrics(t, p, taper_s=job["taper_s"], sigma_s=job["sigma_s"])
    if job["keep"]:
        out.update(t=t, p=p)
    else:
        out.update(t=None, p=None)
    return out


def ground_carpet(ac, phis, start="near_field", circle=None, R_m=None, r0=None, fs=50e3,
                  reflection=1.9, z_ground=0.0, atmosphere=None, humidity=None,
                  z_ref=20000.0, anchor="source", age_from=None, procs=4,
                  taper_s=METRIC_TAPER_S, sigma_s=RISE_SIGMA_S, keep=True, **opts):
    """Ground signatures, lateral offsets and loudness across the boom carpet.

    ``phis``: downward launch azimuths [rad] (0 = undertrack; the carpet is symmetric,
    so non-negative values suffice).  ``start``:
      * "near_field": Whitham near field at r0 (default 1.5 L) with A_e(x, phi), direct
        nonlinear propagation (:func:`ground_signature`);
      * "circle": matching-circle signatures from ``circle`` -- a dict {phi: (t, p)} or a
        callable phi -> (t, p) (called in this process) -- at radius ``R_m``, propagated
        with :func:`propagate_from_circle` (``z_ref``, ``anchor``, ``age_from``).
    Runs the azimuths in a multiprocessing pool of ``procs`` workers (serial if 1 or if
    the job cannot be pickled, e.g. a closure atmosphere).  Extra keywords go to
    :func:`propagate`.

    Returns dict with arrays ``phi``, ``y_g`` [m], ``x_g``, ``delay`` [s] (boom arrival
    after the aircraft passes over the point's x), ``cutoff`` (bool), ``pldb``,
    ``pldb_proxy`` (Gate-0 rise-time proxy), ``peak_pa``, ``csel``, ``rise_ms``; lists
    ``metrics`` (dict per phi from :func:`boomsim.proxy.scored_metrics`, None if cut off)
    and ``t``/``p`` (ground signatures, None if cut off or keep=False); and ``carpet``:
    undertrack PLdB raw/proxy (phi = 0 entry, NaN if absent), max PLdB raw/proxy over the
    carpet, ``half_width_m`` (ground offset of the grazing cut-off ray from
    :func:`lateral_cutoff`) and ``phi_cutoff``.
    """
    import pickle
    from multiprocessing import Pool

    phis = [float(f) for f in np.atleast_1d(phis)]
    r0 = 1.5 * ac.length if r0 is None else float(r0)
    if start not in ("near_field", "circle"):
        raise ValueError("start must be 'near_field' or 'circle'")
    if start == "circle" and (circle is None or R_m is None):
        raise ValueError("start='circle' needs circle and R_m")
    jobs = []
    for phi in phis:
        job = dict(phi=phi, ac=ac, start=start, r0=r0, fs=fs, reflection=reflection,
                   z_ground=z_ground, atmosphere=atmosphere, humidity=humidity,
                   R_m=R_m, z_ref=z_ref, anchor=anchor, age_from=age_from,
                   taper_s=taper_s, sigma_s=sigma_s, keep=keep, opts=opts, circle=None)
        if start == "circle":
            sig = circle(phi) if callable(circle) else circle[phi]
            job["circle"] = (np.asarray(sig[0], float), np.asarray(sig[1], float))
        jobs.append(job)
    parallel = procs > 1 and len(jobs) > 1
    if parallel:
        try:
            pickle.dumps(jobs[0])
        except Exception:
            parallel = False
    if parallel:
        with Pool(min(procs, len(jobs))) as pool:
            rows = pool.map(_carpet_job, jobs)
    else:
        rows = [_carpet_job(j) for j in jobs]

    nan = float("nan")

    def col(key):
        return np.array([r["metrics"][key] if r["metrics"] is not None else nan for r in rows])

    res = dict(phi=np.array(phis), y_g=np.array([r["y_g"] for r in rows]),
               x_g=np.array([r["x_g"] for r in rows]),
               delay=np.array([r["delay"] for r in rows]),
               cutoff=np.array([r["cutoff"] for r in rows]),
               metrics=[r["metrics"] for r in rows],
               t=[r["t"] for r in rows], p=[r["p"] for r in rows],
               pldb=col("pldb"), pldb_proxy=col("pldb_proxy"), peak_pa=col("peak_pa"),
               csel=col("csel"), rise_ms=col("rise_time_10_90_ms"))
    phi_c, y_c = lateral_cutoff(ac.altitude, ac.V, atmosphere=atmosphere, z_ground=z_ground)
    i0 = [i for i, f in enumerate(phis) if f == 0.0]
    ok = ~res["cutoff"]
    res["carpet"] = dict(
        undertrack_pldb=res["pldb"][i0[0]] if i0 else nan,
        undertrack_pldb_proxy=res["pldb_proxy"][i0[0]] if i0 else nan,
        undertrack_peak_pa=res["peak_pa"][i0[0]] if i0 else nan,
        max_pldb=float(np.nanmax(res["pldb"])) if ok.any() else nan,
        max_pldb_proxy=float(np.nanmax(res["pldb_proxy"])) if ok.any() else nan,
        half_width_m=y_c, phi_cutoff=phi_c,
        y_max_computed=float(np.nanmax(np.abs(res["y_g"][ok]))) if ok.any() else nan)
    return res
