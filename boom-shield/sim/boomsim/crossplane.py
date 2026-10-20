"""Cross-plane (y, z) solver for a finite-width BoomShield under a supersonic aircraft.

API
---
CrossPlaneConfig(aircraft, h, W | y_left/y_right, T, profile, gas, m_top, m_bot, m_side,
                 m_part, partitions, tie_K, line_ties, T_rim, lens_q, fringe_width,
                 fringe_profile, fringe_level, f_max, df, f_split, R_m, phis, incident, ...)
matching_signatures(cfg) -> {phi_deg: (tau, p)}  total field on the matching circle R_m
                            (retarded time tau = t - beta R_m/V; unshielded if no shield)
matching_spectra(cfg, X=None) -> (f, P[f, point])
field_at(cfg, f, X, method="auto"|"bem"|"po"), field_map(cfg, f, y, z) -> complex p
incident_coefficients(cfg, f), multipole_field(a, k, y, z), aging_radius(cfg, phi)
incident_vs_whitham(ac, r, phi)
Downstream: ``propagate.propagate_from_circle(ac, radians(phi), t, p, cfg.matching_radius(),
age_from=aging_radius(cfg, radians(phi)))`` (the aged incident already contains the
nonlinear distortion up to aging_radius).
Validation: ``python -m boomsim.crossplane_validation`` -> results/crossplane_validation.*

Contract: ``DESIGN-3D.md``. Steady level flight along x at speed V over a shield that is
invariant along x. Every field is ``P(y, z; omega) exp(j(omega t - k_x x))``, ``k_x =
omega/V``, so each frequency is a 2-D problem in the cross-plane (y lateral, z up;
the flight axis is the origin, the shield lies below it).

* air:            lap p + k^2 p = 0,      k = k_x beta = omega sqrt(1/c^2 - 1/V^2)
* gas, c_g > V:   lap p - kappa^2 p = 0,  kappa = omega sqrt(1/V^2 - 1/c_g^2)
* gas, c_g < V:   lap p + k_g^2 p = 0,    k_g = omega sqrt(1/c_g^2 - 1/V^2)

Conventions as ``DESIGN.md``: spectra are continuous Fourier transforms with kernel
e^{-j omega t} (numpy rfft convention), physical e^{+j omega t}; outgoing waves are
H^(2). Signatures are returned in retarded time tau = t - beta R_m / V (the nose
characteristic of the incident wave reaches the matching circle at tau = 0, exactly
like ``source.near_field``).

Incident field (exact linear slender-body solution)
---------------------------------------------------
    p_inc = S_v(omega) H0^(2)(k r) + S_L(omega) cos(phi) H1^(2)(k r)
    S_v = -j (p gamma M^2 / 4V) FT[A_v''](k_x),   S_L = -(p gamma M^2 / 4V) FT[A_L''](k_x)

with FT[g](k) = int g(xi) e^{-j k xi} d xi over the body (xi from the nose), phi the
azimuth from straight down (cos phi = -z/r), A_v and A_L the volume and lift
equivalent areas of :mod:`boomsim.source`. Derivation: the axial Fourier transform of
the supersonic source/doublet line gives 2-D Helmholtz in the cross-plane; matching the
large-k r asymptote H0 ~ sqrt(2/(pi k r)) e^{-j(k r - pi/4)} (H1 ~ j H0) to Whitham,
dp = p gamma M^2 F(y)/sqrt(2 beta r), F = (1/2pi) int A''(xi)/sqrt(y - xi) d xi, whose
transform is FT[A''] sqrt(pi/(j k))/(2 pi), fixes the constants (the lift part is the
z-derivative of the source line, hence cos(phi) H1 with the same far-field level at
phi = 0). At omega -> 0 the lift term tends to the 2-D lift dipole
p = (p gamma M^2/V) A_L(L) cos(phi)/(2 pi beta r) (net impulse W/(2 pi V r)); that
non-asymptotic part is absent from Whitham's formula. f = 0 is set to this limit and the
shield is taken as transparent there (limp faces: T(0) = 1).

Shield geometry (``CrossPlaneConfig``)
--------------------------------------
Body between y_left + w_f/2 and y_right - w_f/2 (w_f = edge-fringe width, 0 by
default), mid-plane at depth h + T/2 below the axis (``h`` = depth of the top face at
its thickest point). ``profile="slab"`` (thickness T) or ``"lens"``:
T(y) = T_rim + (T - T_rim) (1 - u^2)^q, u = (y - y_c)/a, q = ``lens_q`` (0.5 = elliptic),
closed by vertical rim membranes of height T_rim. Top/bottom faces: limp membranes
(``m_top``/``m_bot`` [kg/m^2] or any object with ``impedance(omega, sx)``, e.g.
``shield.Membrane`` with x-tension); sides ``m_side``; internal vertical partitions at
``partitions`` (y list) with ``m_part``. Ties: smeared ``tie_K`` [Pa/m] (float or
callable of y) between top and bottom faces at the same y, and/or ``line_ties`` =
[(y, k_line [Pa = N/m per m of length])] discrete springs.

Edge fringe (apodised / serrated edge), ``fringe_width`` w_f > 0: on each side a
coverage profile tau(y) (fraction of x covered by shield at lateral position y) runs
from 1 (body) to 0 over w_f, centred on the nominal edge (tau = 1/2 at y_left /
y_right, so the covered area equals that of a sharp edge). Physical picture: plan-view
serrations or a random edge profile of period/correlation length Lambda along x.
The model keeps only the x-averaged (zeroth-order) field. Validity:
  (i)   Lambda << sqrt(lambda R), the along-track Fresnel zone at the receiver distance R
        (lambda = c/f), so the receiver averages over several periods;
  (ii)  higher orders k_x +- 2 pi n/Lambda are dropped: they are evanescent (no leak) if
        Lambda < lambda M/(M + 1); otherwise they radiate as separate weak cones whose total
        power is bounded by the coverage variance, ~ tau (1 - tau) |1 - T_stack|^2;
  (iii) the fringe acts as a graded edge only where w_f >~ lambda_perp/2,
        lambda_perp = V/(f beta) (~170 m at 2 Hz, ~34 m at 10 Hz, Mach 2): at the lowest
        frequencies a fringe narrower than that is invisible and the shield behaves like a
        sharp edge at tau = 1/2.
Models:
  * Kirchhoff path: transmitted aperture field tau T_stack p + (1 - tau) p, reflected
    tau R_stack p (exact coverage averaging of the two local responses).
  * BEM path: a zero-thickness sheet (``fringe_level``: body mid-plane, default, or
    top-face level) with impedance Z_eq(y) chosen so that the plane-wave sheet
    transmission 1/(1 + Z_eq k_z/(2 rho omega)) equals 1 - tau(y) (1 - T_stack) at the
    local incidence (k_y = k y/r from the axis). Exact for transmission at that angle
    (BEM and Kirchhoff agree to ~1 dB under a wide fringe); its reflection 1 - T_sheet
    differs in phase from tau R_stack (matters only above the shield). The fringe teeth
    carry the rim thickness. The sheet level changes the deep-shadow benefit by a few dB
    (both variants are reported in the validation).

Lens rim: ``T_rim`` (minimum rim thickness, default 0.1 T) and ``lens_q`` (taper exponent)
shape the rim; a graded rim is itself an apodisation, but thin rims also leak (review
sec. 3.3): both are in the BEM (rim side membranes of height T_rim).

Solvers
-------
BEM (:mod:`boomsim.bem2d`), used for f <= f_split (default: slab -> the whole band up
to f_max, lens -> 25 Hz): exterior Helmholtz with Burton-Miller (alpha = -j/k; no
fictitious resonances), each gas cell its own K0 (or Helmholtz) domain, membranes,
sides and partitions as coupling panels with p_A - p_B = Z v - f_tie, piecewise-constant
collocation, Laplace singular parts integrated analytically (incl. the finite-part
hypersingular term), elements <= lambda_perp/n (n = n_per_lambda = 10 up to 30 Hz,
tapering to n_per_lambda_hf = 7 at >= 60 Hz), <= gas wavelength/n, <= kappa_res/kappa,
<= W/40, <= 1.5 T; outer corner elements graded 1/4, 1/4, 1/2. Small problems use a
dense LU; slabs above ``n_dense`` panels use a uniform mesh whose face blocks are
Toeplitz (FFT matvecs) with GMRES preconditioned by the exact infinite-slab inverse
(Strang block-circulant, 4x4 per k_y mode) plus a near-field LU for sides/partitions.

Physical optics (Kirchhoff) + exact angular spectrum, for f > f_split (and always
available as ``method="po"``): the incident field on the shield mid-plane is written
as its exact plane-wave spectrum (propagating k_y = k sin(theta), evanescent
k_y = +-k cosh(u)). Inside the covered strip each k_y component is multiplied by the
infinite-stack response at horizontal slowness s_h = sqrt(1/V^2 + k_y^2/omega^2) (the
``shield.transfer`` chain matrices, vectorised; local thickness for a lens), referenced
to the mid-plane (an equal air path removed, so the stack is a thin sheet there);
outside the strip the incident field passes. The aperture differences are radiated with
the exact 2-D Rayleigh integral (equivalent to the angular spectrum) downward
(transmitted) and upward (reflected). Approximations: Kirchhoff truncation at the edges
(black-screen edge waves, thickness lumped at the mid-plane), no partitions, discrete
line ties smeared to K = k_line/spacing, lens face curvature ignored. Accuracy: within
~1 dB of the BEM where transmission dominates; in the deep shadow of an opaque stack the
real edge waves depend on the face impedances (rigid faces diffract more than a black
screen, the lift-closed He stack less), so PO is off by 2-10 dB there (validation).

Accuracy and limitations (see results/crossplane_validation.md)
* BEM transmission floor: the extinction residual of the discretisation limits resolvable
  transmitted power to about -52...-56 dB re incident (He 60 m emulation; mesh dependent);
  levels at the circle below
  ~-45 dB re incident are bounds. Deep-shadow levels move by ~1-3 dB with mesh density at
  70 Hz (7-10 elements/lambda); partitions converge as O(h) (~1-2 dB at -40 dB).
* Lens above f_split is PO only (no structured solver for curved faces).
* The matching-circle field is exact; edge waves arriving there are not radial from the
  axis, so ray marching from the circle is approximate unless R_m >> W (or use a larger
  R_m: cost is only more evaluation points).
* The f = 0 bin and the low-frequency tail carry the 2-D static lift dipole (net impulse
  W/(2 pi V r)) that Whitham's asymptote lacks.

Incident field: ``incident="aged"`` (default) adds to the exact linear field a
multipole correction fitted to the nonlinear aging between the aircraft and the shield
(see :func:`incident_coefficients`, :func:`aging_radius`); ``"linear"`` is the exact
slender-body field.
"""
from __future__ import annotations

import math
import os
from dataclasses import dataclass
from typing import Sequence

import numpy as np
from scipy import special as sp

from . import bem2d
from . import shield as SH
from .atmosphere import GAMMA, std_atm
from .source import Aircraft, lift_area, volume_area

# ======================================================================= config


@dataclass
class CrossPlaneConfig:
    aircraft: Aircraft
    h: float = 250.0                  # depth of the shield top (thickest point) below the axis [m]
    W: float | None = None            # symmetric width (sets y_left/right = -+W/2); None = no shield
    T: float = 30.0                   # gas thickness [m] (lens: centre thickness)
    profile: str = "slab"             # "slab" | "lens"
    gas: object = None                # shield.Gas (default: pure helium at ambient T, p)
    m_top: object = None              # kg/m^2 or object with .impedance(omega, sx); None = lift-closed
    m_bot: object = None
    m_side: object = None             # None -> m_bot
    m_part: object = None             # None -> m_bot
    partitions: Sequence = ()         # y positions of internal vertical partitions
    tie_K: object = 0.0               # Pa/m, float or callable(y)
    line_ties: Sequence = ()          # [(y, k_line [Pa])]
    y_left: float | None = None       # asymmetric edges (override W)
    y_right: float | None = None
    T_rim: float | None = None        # lens rim thickness (default 0.1 T)
    lens_q: float = 0.5               # lens taper exponent (0.5 elliptic)
    fringe_width: float = 0.0         # edge fringe width w_f [m]
    fringe_profile: object = "raised_cosine"  # "linear" | "raised_cosine" | callable(u in [0,1]) -> tau
    fringe_level: str = "mid"         # BEM fringe sheet at the body "mid"-plane or "top"-face level
    f_max: float = 100.0              # Hz, band limit (raised-cosine taper from taper_from*f_max)
    df: float | None = None           # Hz; None = from geometry (record covers all edge-wave delays)
    f_split: float | None = None      # Hz, BEM below / Kirchhoff above (None: slab -> f_max, lens -> 25)
    n_dense: int = 700                # slab: dense BEM up to this many panels, structured above
    R_m: float | None = None          # matching-circle radius (default: encloses the shield)
    phis: Sequence | None = None      # degrees; default -70..70 step 10
    atm: object = None                # uniform medium (default std_atm(altitude - h))
    fs: float | None = None           # output sample rate (default 10 f_max)
    t_pre: float = 0.4                # s of record before tau = 0
    taper_from: float = 0.75
    n_per_lambda: float = 10.0
    n_per_lambda_hf: float = 7.0      # used from 60 Hz up (linear taper from 30 Hz)
    kappa_res: float = 1.0            # element <= kappa_res / kappa
    workers: int = 4
    incident: str = "aged"            # "aged" (aged to the shield, see incident_coefficients) | "linear"

    # ---- derived -------------------------------------------------------------
    def __post_init__(self):
        if self.W is not None and self.y_left is None:
            self.y_left, self.y_right = -0.5 * self.W, 0.5 * self.W
        if self.y_left is not None:
            self.W = float(self.y_right - self.y_left)
        if self.atm is None:
            self.atm = std_atm(self.aircraft.altitude - self.h)
        if self.has_shield:
            if self.gas is None:
                self.gas = SH.gas("helium", float(self.atm.T), float(self.atm.p))
            face = max((float(self.atm.rho) - self.gas.rho) * self.T / 2, 1e-3)
            if self.m_top is None:
                self.m_top = face
            if self.m_bot is None:
                self.m_bot = face
            if self.m_side is None:
                self.m_side = self.m_bot
            if self.m_part is None:
                self.m_part = self.m_bot
            if self.T_rim is None:
                self.T_rim = 0.1 * self.T
        if self.phis is None:
            self.phis = np.arange(-70.0, 70.1, 10.0)
        self.phis = np.asarray(self.phis, float)
        if self.fs is None:
            self.fs = 10.0 * self.f_max

    @property
    def has_shield(self) -> bool:
        return self.y_left is not None and self.T > 0 and self.W > 0

    @property
    def V(self) -> float:
        return self.aircraft.V

    @property
    def c(self) -> float:
        return float(self.atm.c)

    @property
    def rho(self) -> float:
        return float(self.atm.rho)

    @property
    def beta(self) -> float:
        """beta in the cross-plane medium: k_perp = (omega/V) beta."""
        return math.sqrt((self.V / self.c) ** 2 - 1.0)

    @property
    def z_mid(self) -> float:
        return -self.h - 0.5 * self.T

    @property
    def body(self):
        """(y0, y1) of the full-coverage body."""
        return self.y_left + 0.5 * self.fringe_width, self.y_right - 0.5 * self.fringe_width

    @property
    def extent(self):
        """(y0, y1) of everything (body + fringes)."""
        return self.y_left - 0.5 * self.fringe_width, self.y_right + 0.5 * self.fringe_width

    def thickness(self, y):
        y = np.asarray(y, float)
        y0, y1 = self.body
        inside = (y >= y0 - 1e-9) & (y <= y1 + 1e-9)
        if self.profile == "slab":
            d = np.full(y.shape, float(self.T))
        elif self.profile == "lens":
            u = np.clip((y - 0.5 * (y0 + y1)) / (0.5 * (y1 - y0)), -1, 1)
            d = self.T_rim + (self.T - self.T_rim) * (1 - u * u) ** self.lens_q
        else:
            raise ValueError("profile must be 'slab' or 'lens'")
        rim = self.T if self.profile == "slab" else self.T_rim
        return np.where(inside, d, rim)   # fringe teeth carry the rim thickness

    def coverage(self, y):
        """tau(y): 1 on the body, fringe profile on each side, 0 outside."""
        y = np.asarray(y, float)
        if not self.has_shield:
            return np.zeros_like(y)
        y0, y1 = self.body
        tau = ((y >= y0) & (y <= y1)).astype(float)
        wf = self.fringe_width
        if wf > 0:
            prof = self.fringe_profile
            if callable(prof):
                g = prof
            elif prof == "linear":
                g = lambda u: 1.0 - u
            elif prof == "raised_cosine":
                g = lambda u: 0.5 * (1.0 + np.cos(np.pi * u))
            else:
                raise ValueError("fringe_profile")
            for ye, sgn in ((y0, -1.0), (y1, 1.0)):
                u = sgn * (y - ye) / wf
                m = (u > 0) & (u <= 1)
                tau[m] = np.clip(g(u[m]), 0.0, 1.0)
        return tau

    def tie_K_at(self, y):
        y = np.asarray(y, float)
        K = self.tie_K(y) if callable(self.tie_K) else np.full(y.shape, float(self.tie_K))
        return np.asarray(K, float) * np.ones_like(y)

    def stack(self, d, K=0.0):
        return SH.Stack([_as_membrane(self.m_top), SH.GasLayer(self.gas, d, tie_K=K),
                         _as_membrane(self.m_bot)])

    def corners(self):
        if not self.has_shield:
            return np.zeros((0, 2))
        e0, e1 = self.extent
        return np.array([[e0, -self.h], [e1, -self.h], [e0, -self.h - self.T], [e1, -self.h - self.T]])

    def matching_radius(self) -> float:
        if self.R_m is not None:
            return float(self.R_m)
        if not self.has_shield:
            return 2.0 * self.h
        return 1.25 * float(np.max(np.hypot(*self.corners().T)))

    def circle_points(self, R=None):
        R = self.matching_radius() if R is None else R
        ph = np.radians(self.phis)
        return np.stack([R * np.sin(ph), -R * np.cos(ph)], axis=1)

    def auto_df(self) -> float:
        L, V = self.aircraft.length, self.V
        R = self.matching_radius()
        dmax = 0.0
        if self.has_shield:
            X = self.circle_points(R)
            for ye in self.extent:
                e = np.array([ye, self.z_mid])
                dmax = max(dmax, float(np.max(np.hypot(*e) + np.hypot(*(X - e).T) - R)))
        T_rec = 1.3 * (self.beta * dmax / V + 3 * L / V) + self.t_pre + 1.0
        T_rec = max(2.0, math.ceil(T_rec * 2) / 2)
        return 1.0 / T_rec

    def split_frequency(self) -> float:
        if self.f_split is not None:
            return float(self.f_split)
        return float(self.f_max) if self.profile == "slab" else 25.0

    def frequencies(self):
        df = self.df if self.df is not None else self.auto_df()
        n = int(math.floor(self.f_max / df + 1e-9))
        return np.arange(n + 1) * df


def _as_membrane(m):
    return m if hasattr(m, "impedance") else SH.Membrane(float(m))


def _impedance(m, omega, sx):
    if hasattr(m, "impedance"):
        return complex(np.asarray(m.impedance(np.array([omega]), sx)).ravel()[0])
    return 1j * omega * float(m)


# ======================================================================= incident field


def _area_second_derivative_cells(ac: Aircraft, n=4000):
    """Cell-constant A_v'' and A_L'' on n cells over [0, L_end] (slopes from central
    differences of the source areas; exact cell means of A'' up to O(delta^2))."""
    xe = ac.length * max(1.0, ac.lift_end)
    h = xe / n
    x = np.arange(n + 1) * h
    dl = 1e-4 * h

    def slope(fun):
        s = (fun(ac, x + dl) - fun(ac, np.maximum(x - dl, 0.0))) / (x + dl - np.maximum(x - dl, 0.0))
        s[0] = (fun(ac, np.array([dl]))[0] - fun(ac, np.array([0.0]))[0]) / dl
        return s

    sv = slope(volume_area)
    sl = slope(lift_area)
    sv[-1] = 0.0
    sl[-1] = 0.0
    return x[:-1], h, np.diff(sv) / h, np.diff(sl) / h, float(lift_area(ac, np.array([xe]))[0])


def source_strengths(ac: Aircraft, f, n_cells=4000):
    """(S_v, S_L) [Pa s] at frequencies f so that p_inc = S_v H0(k r) + S_L cos(phi) H1(k r).
    f = 0 returns 0 (use :func:`incident_field`, which inserts the static lift limit)."""
    f = np.atleast_1d(np.asarray(f, float))
    a = ac.atm
    q = float(a.p) * GAMMA * ac.mach ** 2           # = rho V^2
    V = ac.V
    x0, h, A2v, A2l, _ = _area_second_derivative_cells(ac, n_cells)
    kx = 2 * np.pi * f / V
    out_v = np.zeros(f.size, complex)
    out_l = np.zeros(f.size, complex)
    for i0 in range(0, f.size, 64):
        kk = kx[i0:i0 + 64, None]
        ph = np.exp(-1j * kk * x0[None, :])
        cell = np.where(kk == 0, h, (1 - np.exp(-1j * kk * h)) / np.where(kk == 0, 1, 1j * kk))
        Fv = (ph * A2v[None, :]).sum(1) * cell[:, 0]
        Fl = (ph * A2l[None, :]).sum(1) * cell[:, 0]
        out_v[i0:i0 + 64] = -1j * q / (4 * V) * Fv
        out_l[i0:i0 + 64] = -q / (4 * V) * Fl
    out_v[f == 0] = 0
    out_l[f == 0] = 0
    return out_v, out_l


def static_lift_coefficient(ac: Aircraft, beta_medium=None):
    """f -> 0 limit: p_inc(0) = C cos(phi) / r [Pa s]."""
    a = ac.atm
    q = float(a.p) * GAMMA * ac.mach ** 2
    _, _, _, _, AL = _area_second_derivative_cells(ac, 200)
    b = ac.beta if beta_medium is None else beta_medium
    return q / ac.V * AL / (2 * np.pi * b)


def multipole_field(a, k, Y, Zc, grad=False):
    """p = sum_n a_n cos(n phi) H_n^(2)(k r), phi = atan2(y, -z) (azimuth from straight
    down); with grad=True also returns dp/dy, dp/dz. a = [S_v, S_L] is the linear
    slender-body field."""
    a = np.atleast_1d(np.asarray(a, complex))
    Y = np.asarray(Y, float)
    Zc = np.asarray(Zc, float)
    r = np.hypot(Y, Zc)
    z = k * r
    ph = np.arctan2(Y, -Zc)
    nmax = len(a) - 1
    H = [sp.hankel2(n, z) for n in range(nmax + 2)]
    p = np.zeros(np.broadcast(Y, Zc).shape, complex)
    gy = np.zeros_like(p)
    gz = np.zeros_like(p)
    sph, cph = np.sin(ph), np.cos(ph)
    for n, an in enumerate(a):
        if an == 0:
            continue
        cn = np.cos(n * ph)
        p += an * cn * H[n]
        if grad:
            Hp = -H[1] if n == 0 else 0.5 * (H[n - 1] - H[n + 1])
            dr = an * cn * k * Hp                       # d/dr
            dphi = -an * n * np.sin(n * ph) * H[n] / r  # (1/r) d/dphi
            gy += dr * sph + dphi * cph
            gz += -dr * cph + dphi * sph
    if grad:
        return p, gy, gz
    return p


def incident_field(Sv, SL, k, Y, Zc, grad=False):
    """p (and dp/dy, dp/dz) of S_v H0(k r) + S_L cos(phi) H1(k r)."""
    return multipole_field([Sv, SL], k, Y, Zc, grad)


def _spectral_weight(a, theta):
    """Plane-wave amplitude of sum_n a_n cos(n phi) H_n(k r) in the lower half-space:
    cos(n phi) H_n(k r) = (1/pi) int j^n cos(n theta) e^{-j k (y sin theta + |z| cos theta)} d theta
    (theta = +-(pi/2 + j u) on the evanescent branches, d theta -> j du)."""
    a = np.atleast_1d(np.asarray(a, complex))
    out = np.zeros(np.shape(theta), complex)
    for n, an in enumerate(a):
        if an != 0:
            out += an * (1j ** n) * np.cos(n * theta)
    return out


AGING_PHIS = np.radians(np.arange(0.0, 80.1, 5.0))


def aging_radius(cfg, phi):
    """Radius [m] at which the incident field is matched to the aged near field at
    azimuth phi [rad]: the slant distance to the shield's top plane, h / cos(phi),
    capped at the matching radius. Pass it as ``age_from`` to
    ``propagate.propagate_from_circle`` for consistency."""
    phi = np.asarray(phi, float)
    return np.minimum(cfg.h / np.maximum(np.cos(phi), 1e-6), cfg.matching_radius())


def _signature_spectrum(t, p, f):
    """Continuous FT int p e^{-j 2 pi f t} dt at frequencies f (FFT when f = n df)."""
    t = np.asarray(t, float)
    p = np.asarray(p, float)
    f = np.asarray(f, float)
    dt = t[1] - t[0]
    if f.size > 2 and f[0] == 0 and np.allclose(np.diff(f), f[1]):
        N = int(round(1.0 / (f[1] * dt)))
        if N >= p.size and abs(N * f[1] * dt - 1) < 1e-9:
            X = np.fft.rfft(p, N)[:f.size]
            if X.size == f.size:
                return X * dt * np.exp(-2j * np.pi * f * t[0])
    out = np.empty(len(f), complex)
    for i0 in range(0, len(f), 32):
        out[i0:i0 + 32] = np.exp(-2j * np.pi * np.outer(f[i0:i0 + 32], t)) @ p * dt
    return out


def incident_coefficients(cfg, f):
    """Multipole coefficients A[i, n] of the incident field at frequencies f.

    ``cfg.incident == "linear"``: A = [S_v, S_L] (exact slender-body field).
    ``"aged"``: the linear coefficients plus a least-squares correction dA[i, n],
    n <= min(8, 0.4 k h), chosen so that at radius r_j = aging_radius(phi_j) on
    azimuths phi_j = 0..80 deg (step 5) the field changes by exactly what nonlinear aging
    changes the signature there: FT[propagate.aged_near_field(r_j, phi_j)] -
    FT[source.near_field(r_j, phi_j)] (uniform-air matching convention; converted from
    retarded to absolute phase with exp(-j omega beta r_j / V)). The exact Hankel
    near-zone structure of the linear field is kept; only the aging increment is carried
    by the fitted multipoles (symmetric in phi, as the source is)."""
    f = np.atleast_1d(np.asarray(f, float))
    Sv, SL = source_strengths(cfg.aircraft, f)
    if cfg.incident == "linear":
        return np.stack([Sv, SL], axis=1)
    if cfg.incident != "aged":
        raise ValueError("incident must be 'linear' or 'aged'")
    key = (repr(cfg.aircraft), cfg.h, cfg.matching_radius(), cfg.c, f.tobytes())
    if key in _COEF_CACHE:
        return _COEF_CACHE[key]
    from . import propagate as PR
    from .source import near_field
    ac = cfg.aircraft
    nmax = 8
    A = np.zeros((f.size, nmax + 1), complex)
    A[:, 0], A[:, 1] = Sv, SL
    phis = AGING_PHIS
    rj = aging_radius(cfg, phis)
    fs = max(20000.0, 40 * float(np.max(f)))
    D = np.zeros((phis.size, f.size), complex)
    for j, (ph, r) in enumerate(zip(phis, rj)):
        ta, pa = PR.aged_near_field(ac, float(r), fs, phi=float(ph))
        tl, pl = near_field(ac, float(r), fs, phi=float(ph))
        D[j] = _signature_spectrum(ta, pa, f) - _signature_spectrum(tl, pl, f)
    w = 2 * np.pi * f
    k = w * math.sqrt(1 / cfg.c ** 2 - 1 / cfg.V ** 2)
    for i in range(f.size):
        if f[i] == 0:
            continue
        N = int(np.clip(np.floor(0.4 * k[i] * cfg.h), 1, nmax))
        B = np.stack([np.cos(n * phis) * sp.hankel2(n, k[i] * rj) for n in range(N + 1)], axis=1)
        rhs = D[:, i] * np.exp(-1j * w[i] * cfg.beta * rj / cfg.V)
        sc = np.linalg.norm(B, axis=0)
        x, *_ = np.linalg.lstsq(B / sc, rhs, rcond=1e-8)
        A[i, :N + 1] += x / sc
    _COEF_CACHE[key] = A
    return A


_COEF_CACHE = {}


# ======================================================================= geometry / BEM


def _wavenumbers(cfg: CrossPlaneConfig, f):
    w = 2 * np.pi * f
    k = w * math.sqrt(1 / cfg.c ** 2 - 1 / cfg.V ** 2)
    g = cfg.gas
    if g is None:
        return w, k, ("mod", 0.0)
    s2 = 1 / g.c ** 2 - 1 / cfg.V ** 2
    if s2 > 0:
        return w, k, ("helm", w * math.sqrt(s2))
    return w, k, ("mod", w * math.sqrt(-s2))


def _split_interval(s0, s1, h, grade=(True, True)):
    """Nodes on [s0, s1] with spacing <= h; the end element at a graded end is split
    1/4, 1/4, 1/2 (outer shield corners only: grading at partition junctions slowed
    convergence in tests)."""
    n = max(1, int(math.ceil((s1 - s0) / h - 1e-9)))
    e = np.linspace(s0, s1, n + 1)
    d = e[1] - e[0]
    extra = []
    if grade[0]:
        extra += [s0 + 0.25 * d, s0 + 0.5 * d]
    if grade[1]:
        extra += [s1 - 0.25 * d, s1 - 0.5 * d]
    return np.unique(np.concatenate([e, extra])) if extra else e


def build_mesh(cfg: CrossPlaneConfig, f: float, h_elem: float | None = None, uniform: bool = False):
    """Panels + bookkeeping for the BEM at frequency f.

    uniform=True (slab only): equal face panels (no corner grading), partitions snapped
    to the nearest face node, fringe sheets with the face spacing, so that all face and
    sheet panels form uniform collinear groups (``mesh["groups"]``) for the structured
    FFT solver."""
    w, k, (gkind, gk) = _wavenumbers(cfg, f)
    lam = 2 * np.pi / k
    # elements per exterior wavelength: n_per_lambda up to 30 Hz, tapering linearly to
    # n_per_lambda_hf at 60 Hz and above (deep-shadow levels move by ~1-3 dB, see docs)
    npl = float(np.interp(f, [30.0, 60.0], [cfg.n_per_lambda, min(cfg.n_per_lambda, cfg.n_per_lambda_hf)]))
    h_e = lam / npl
    if gkind == "helm":
        h_e = min(h_e, 2 * np.pi / gk / npl)
    else:
        h_e = min(h_e, cfg.kappa_res / max(gk, 1e-12))
    y0, y1 = cfg.body
    h_e = min(h_e, (y1 - y0) / 40.0, 1.5 * cfg.T)
    if cfg.fringe_width > 0:
        h_e_fr = min(h_e, cfg.fringe_width / 12.0)
    if h_elem is not None:
        h_e = h_elem
    zc = cfg.z_mid
    # ---- face nodes in y, equal arc length on the top curve between breakpoints
    parts = sorted(p for p in cfg.partitions if y0 < p < y1)
    if uniform:
        if cfg.profile != "slab":
            raise ValueError("uniform mesh needs a slab")
        n_u = int(math.ceil((y1 - y0) / h_e))
        h_e = (y1 - y0) / n_u
        parts = sorted(set(float(y0 + round((p - y0) / h_e) * h_e) for p in parts))
        parts = [p for p in parts if y0 < p < y1]
    brk = [y0] + parts + [y1]
    yc, a = 0.5 * (y0 + y1), 0.5 * (y1 - y0)
    th_f = np.linspace(-np.pi / 2, np.pi / 2, 40001)
    yf = yc + a * np.sin(th_f)
    zf = 0.5 * cfg.thickness(yf)
    sf = np.concatenate([[0], np.cumsum(np.hypot(np.diff(yf), np.diff(zf)))])
    ynodes = []
    for b0, b1 in zip(brk[:-1], brk[1:]):
        s0, s1 = np.interp([b0, b1], yf, sf)
        gr = (not uniform and b0 == y0, not uniform and b1 == y1)
        se = _split_interval(s0, s1, h_e * (1 + 1e-9), grade=gr)
        ye = np.interp(se, sf, yf)
        ye[0], ye[-1] = b0, b1
        ynodes.append(ye if not ynodes else ye[1:])
    yn = np.concatenate(ynodes)
    dn = cfg.thickness(yn)
    dn[0] = dn[-1] = cfg.T if cfg.profile == "slab" else cfg.T_rim
    top = np.stack([yn, zc + 0.5 * dn], 1)
    bot = np.stack([yn, zc - 0.5 * dn], 1)
    nf = yn.size - 1
    cell_of_face = np.searchsorted(np.asarray(parts), 0.5 * (yn[:-1] + yn[1:])) + 1  # 1..ncell
    ncell = len(parts) + 1
    A, Bp, sa, sb, kind, Zs, tag = [], [], [], [], [], [], []
    sx = 1.0 / cfg.V
    Zt, Zb = _impedance(cfg.m_top, w, sx), _impedance(cfg.m_bot, w, sx)
    Zsd, Zp = _impedance(cfg.m_side, w, sx), _impedance(cfg.m_part, w, sx)
    # top: left->right, n down, A=ext(0), B=cell
    A.append(top[:-1]); Bp.append(top[1:]); sa.append(np.zeros(nf, int)); sb.append(cell_of_face)
    kind.append(np.zeros(nf, int)); Zs.append(np.full(nf, Zt)); tag.append(np.full(nf, 1))
    # bottom: left->right, n down, A=cell, B=ext
    A.append(bot[:-1]); Bp.append(bot[1:]); sa.append(cell_of_face); sb.append(np.zeros(nf, int))
    kind.append(np.zeros(nf, int)); Zs.append(np.full(nf, Zb)); tag.append(np.full(nf, 2))

    def vertical(yv, zlo, zhi, a_dom, b_dom, Z, tg, even=True):
        n = max(2, int(math.ceil((zhi - zlo) / h_e)))
        if even and n % 2:
            n += 1
        zz = np.linspace(zlo, zhi, n + 1)
        P0 = np.stack([np.full(n, yv), zz[:-1]], 1)
        P1 = np.stack([np.full(n, yv), zz[1:]], 1)
        A.append(P0); Bp.append(P1); sa.append(np.full(n, a_dom)); sb.append(np.full(n, b_dom))
        kind.append(np.zeros(n, int)); Zs.append(np.full(n, Z)); tag.append(np.full(n, tg))

    # sides (bottom -> top => n = +y)
    vertical(y0, bot[0, 1], top[0, 1], 0, 1, Zsd, 3)
    vertical(y1, bot[-1, 1], top[-1, 1], ncell, 0, Zsd, 3)
    for i, yp in enumerate(parts):
        dp = float(cfg.thickness(np.array([yp]))[0])
        vertical(yp, zc - 0.5 * dp, zc + 0.5 * dp, i + 1, i + 2, Zp, 4, even=False)
    # fringe sheets (left->right, n down, exterior both sides) at the body mid-plane
    # (default; matches the coverage-averaged transmission to ~1 dB) or at the top-face
    # level ("top": continues the top face, no exposed corner at the fringe root)
    if cfg.fringe_width > 0:
        wf = cfg.fringe_width
        z_sh = top[0, 1] if cfg.fringe_level == "top" else zc
        for s0, s1 in ((y0 - wf, y0), (y1, y1 + wf)):
            if uniform:
                n = max(4, int(round(wf / h_e)))
                ys = (np.arange(n + 1) * h_e + s0) if s0 >= y1 else (s1 - h_e * np.arange(n, -1, -1))
            else:
                n = max(4, int(math.ceil(wf / h_e_fr)))
                ys = np.linspace(s0, s1, n + 1)
            A.append(np.stack([ys[:-1], np.full(n, z_sh)], 1)); Bp.append(np.stack([ys[1:], np.full(n, z_sh)], 1))
            sa.append(np.zeros(n, int)); sb.append(np.zeros(n, int)); kind.append(np.ones(n, int))
            Zs.append(np.zeros(n, complex)); tag.append(np.full(n, 5))
    sizes = np.cumsum([0] + [len(a_) for a_ in A])
    groups = [np.arange(sizes[0], sizes[1]), np.arange(sizes[1], sizes[2])]
    if cfg.fringe_width > 0:
        groups += [np.arange(sizes[-3], sizes[-2]), np.arange(sizes[-2], sizes[-1])]
    P = bem2d.Panels(np.concatenate(A), np.concatenate(Bp))
    mesh = dict(P=P, side_a=np.concatenate(sa), side_b=np.concatenate(sb), kind=np.concatenate(kind),
                Z=np.concatenate(Zs).astype(complex), tag=np.concatenate(tag), nf=nf, ncell=ncell,
                h_e=h_e, w=w, k=k, gkind=gkind, gk=gk, groups=groups, uniform=uniform)
    # ---- fringe sheet impedances
    if cfg.fringe_width > 0:
        sh = np.nonzero(mesh["kind"] == 1)[0]
        yy = P.c[sh, 0]
        tau = cfg.coverage(yy)
        ky = k * yy / np.hypot(yy, P.c[sh, 1])       # local incidence from the axis
        kz = np.sqrt(k * k - ky * ky)
        rim = cfg.T if cfg.profile == "slab" else cfg.T_rim
        Kt = cfg.tie_K_at(yy)
        Tst = np.array([_stack_TR(cfg, w, np.array([kyi]), rim, Ki)[0][0] for kyi, Ki in zip(ky, Kt)])
        Tt = 1.0 - tau * (1.0 - Tst)
        Tt = np.where(np.abs(Tt) < 1e-12, 1e-12, Tt)
        mesh["Z"][sh] = (2 * cfg.rho * w / kz) * (1.0 / Tt - 1.0)
    # ---- ties
    ties = []
    it, ib = np.arange(nf), nf + np.arange(nf)
    nzt, nzb = P.n[it, 1], P.n[ib, 1]
    dt, db = -nzt, -nzb
    K = cfg.tie_K_at(P.c[it, 0])
    s = 1.0 / (1j * w)
    for i in np.nonzero(K != 0)[0]:
        c_t = K[i] * s * abs(nzt[i]) * dt[i]
        c_b = K[i] * s * abs(nzb[i]) * db[i]
        ties += [(it[i], ib[i], c_t * db[i]), (it[i], it[i], -c_t * dt[i]),
                 (ib[i], ib[i], -c_b * db[i]), (ib[i], it[i], c_b * dt[i])]
    for (yk, kl) in cfg.line_ties:
        i = int(np.clip(np.searchsorted(yn, yk) - 1, 0, nf - 1))
        c_t = kl * s / P.L[it[i]] * dt[i]
        c_b = kl * s / P.L[ib[i]] * db[i]
        ties += [(it[i], ib[i], c_t * db[i]), (it[i], it[i], -c_t * dt[i]),
                 (ib[i], ib[i], -c_b * db[i]), (ib[i], it[i], c_b * dt[i])]
    mesh["ties"] = ties
    return mesh


def bem_solve(cfg: CrossPlaneConfig, f: float, a=None, mesh=None, inc=None, fast_tol=1e-5):
    """BEM solution at one frequency. ``inc(Y, Z) -> (p, py, pz)`` overrides the aircraft
    incident field (used by tests). Slabs with more than ``cfg.n_dense`` panels use the
    structured FFT/GMRES solver on a uniform mesh, everything else the dense solver."""
    if mesh is None:
        mesh = build_mesh(cfg, f)
        if cfg.profile == "slab" and len(mesh["P"]) > cfg.n_dense:
            mesh = build_mesh(cfg, f, uniform=True)
    w, k = mesh["w"], mesh["k"]
    if inc is None:
        if a is None:
            a = incident_coefficients(cfg, [f])[0]
        inc = lambda Y, Z: multipole_field(a, k, Y, Z, grad=True)
    P = mesh["P"]
    gas = cfg.gas
    gk = mesh["gk"]
    gkern = bem2d.Helmholtz(gk) if mesh["gkind"] == "helm" else bem2d.ModHelmholtz(gk)
    doms = [bem2d.Domain(bem2d.Helmholtz(k), cfg.rho, exterior=True)]
    doms += [bem2d.Domain(gkern, gas.rho) for _ in range(mesh["ncell"])]
    prob = bem2d.Problem(w, P, doms, mesh["side_a"], mesh["side_b"], mesh["kind"], mesh["Z"],
                         ties=mesh["ties"])
    p, py, pz = inc(P.c[:, 0], P.c[:, 1])
    dpdn = py * P.n[:, 0] + pz * P.n[:, 1]
    if mesh.get("uniform"):
        sol = bem2d.solve_fast(prob, p, dpdn, mesh["groups"], tol=fast_tol)
    else:
        sol = bem2d.solve(prob, p, dpdn)
    sol.inc = inc
    sol.mesh = mesh
    return sol


def bem_field(sol, X):
    X = np.asarray(X, float).reshape(-1, 2)
    p0 = sol.inc(X[:, 0], X[:, 1])[0]
    return sol.field(X, inc_p=p0)


def bem_grad(sol, X, dirs):
    X = np.asarray(X, float).reshape(-1, 2)
    _, py, pz = sol.inc(X[:, 0], X[:, 1])
    d0 = py * dirs[:, 0] + pz * dirs[:, 1]
    return sol.field(X, inc_p=d0, grad_dirs=dirs)


# ======================================================================= Kirchhoff / PO


def _stack_TR(cfg, omega, ky, d, K=0.0):
    """Infinite-stack (T, R) for horizontal wavenumbers (k_x = omega/V, ky), both
    referenced to the stack mid-plane (equal air path removed: thin-sheet equivalent).
    ky: array (any shape); d: thickness (scalar or array broadcastable to ky)."""
    ky = np.asarray(ky, float)
    d = np.asarray(d, float) * np.ones_like(ky)
    sx = np.sqrt(1.0 / cfg.V ** 2 + (ky / omega) ** 2)
    om = np.full(ky.shape, float(omega))
    st = SH.Stack([_as_membrane(cfg.m_top), SH.GasLayer(cfg.gas, d, tie_K=float(K)),
                   _as_membrane(cfg.m_bot)])
    A, B, C, D, L = SH._stack_matrix(st, om, sx)
    kza = SH._kz(om, cfg.c, sx)
    Za = cfg.rho * om / kza
    den = A + B / Za + Za * C + D
    ph = np.exp(np.minimum((1j * kza * d).real, 700.0) + 1j * (1j * kza * d).imag)  # e^{j kz d}
    T = 2.0 * np.exp(-L) / den * ph
    R = (2.0 * (A + B / Za) / den - 1.0) * ph
    return T, R


def _gl(n, a, b):
    x, w = sp.roots_legendre(int(n))
    return 0.5 * (b - a) * x + 0.5 * (a + b), 0.5 * (b - a) * w


def po_solve(cfg: CrossPlaneConfig, f: float, a=None, theta_factor=1.0, dy_lambda=10.0):
    """Physical-optics aperture fields at frequency f. Returns a dict with the aperture
    grid y, spacing dy, plane z_a, transmitted difference d(y) and reflected r(y)."""
    w, k, _ = _wavenumbers(cfg, f)
    if a is None:
        a = incident_coefficients(cfg, [f])[0]
    za = cfg.z_mid
    e0, e1 = cfg.extent
    n_y = max(16, int(math.ceil((e1 - e0) / (2 * np.pi / k / dy_lambda))))
    dy = (e1 - e0) / n_y
    y = e0 + (np.arange(n_y) + 0.5) * dy
    tau = cfg.coverage(y)
    dthk = cfg.thickness(y)
    K = cfg.tie_K_at(y)
    if len(cfg.line_ties):
        yt = np.sort([t[0] for t in cfg.line_ties])
        sp_ = np.gradient(yt) if yt.size > 1 else np.array([cfg.W])
        for (yk, kl), s in zip(sorted(cfg.line_ties), sp_):
            m = np.abs(y - yk) <= 0.5 * s
            K[m] += kl / s
    rho_max = float(np.max(np.hypot(np.abs(za), np.maximum(abs(e0), abs(e1)))))
    n_th = int(theta_factor * (0.7 * k * rho_max * np.pi / 2 + 60))
    th, wth = _gl(n_th, -np.pi / 2, np.pi / 2)
    ky_p = k * np.sin(th)
    wp = (1 / np.pi) * np.exp(-1j * k * np.cos(th) * abs(za)) * _spectral_weight(a, th) * wth
    u_max = math.asinh(40.0 / (k * abs(za)))
    uu, wu = _gl(48, 0.0, u_max)
    ky_e = k * np.cosh(uu)
    we = (1 / np.pi) * np.exp(-k * np.sinh(uu) * abs(za)) * _spectral_weight(a, np.pi / 2 + 1j * uu) * 1j * wu
    ky_all = np.concatenate([ky_p, ky_e, -ky_e])
    w_all = np.concatenate([wp, we, we])
    E = np.exp(-1j * np.outer(y, ky_all))
    dfield = np.zeros(n_y, complex)
    rfield = np.zeros(n_y, complex)
    groups = {}
    for i, (dd, kk) in enumerate(zip(dthk, K)):
        groups.setdefault((round(float(dd), 6), round(float(kk), 9)), []).append(i)
    if cfg.profile == "slab" and len(groups) <= 4:
        for (dd, kk), idx in groups.items():
            Tm, Rm = _stack_TR(cfg, w, np.abs(ky_all), dd, kk)
            idx = np.asarray(idx)
            dfield[idx] = E[idx] @ ((Tm - 1.0) * w_all)
            rfield[idx] = E[idx] @ (Rm * w_all)
    else:
        for kk in np.unique(np.round(K, 9)):
            idx = np.nonzero(np.round(K, 9) == kk)[0]
            Tm, Rm = _stack_TR(cfg, w, np.abs(ky_all)[None, :], dthk[idx][:, None], kk)
            dfield[idx] = np.sum(E[idx] * (Tm - 1.0) * w_all[None, :], axis=1)
            rfield[idx] = np.sum(E[idx] * Rm * w_all[None, :], axis=1)
    return dict(y=y, dy=dy, za=za, d=tau * dfield, r=tau * rfield, k=k, a=a)


def po_field(ap, X):
    """Total field at points X (M, 2) from the PO aperture fields (Rayleigh integrals)."""
    X = np.asarray(X, float).reshape(-1, 2)
    k, za = ap["k"], ap["za"]
    p = multipole_field(ap["a"], k, X[:, 0], X[:, 1])
    below = X[:, 1] < za
    for sel, src in ((below, ap["d"]), (~below, ap["r"])):
        if not np.any(sel):
            continue
        Xs = X[sel]
        out = np.zeros(Xs.shape[0], complex)
        for i0 in range(0, Xs.shape[0], 512):
            xx = Xs[i0:i0 + 512]
            dyv = xx[:, None, 0] - ap["y"][None, :]
            dz = np.abs(xx[:, None, 1] - za)
            rr = np.hypot(dyv, dz)
            out[i0:i0 + 512] = (-0.5j * k) * np.sum(src[None, :] * sp.hankel2(1, k * rr) * dz / rr, axis=1) * ap["dy"]
        p[sel] += out
    return p


# ======================================================================= drivers


def field_at(cfg: CrossPlaneConfig, f: float, X, method="auto", a=None):
    """Total complex pressure spectrum [Pa s] at points X (M, 2) at frequency f > 0.
    method: "auto" (BEM if f <= f_split), "bem", "po"."""
    X = np.asarray(X, float).reshape(-1, 2)
    w, k, _ = _wavenumbers(cfg, f)
    if a is None:
        a = incident_coefficients(cfg, [f])[0]
    if not cfg.has_shield:
        return multipole_field(a, k, X[:, 0], X[:, 1])
    if method == "auto":
        method = "bem" if f <= cfg.split_frequency() else "po"
    if method == "bem":
        return bem_field(bem_solve(cfg, f, a), X)
    return po_field(po_solve(cfg, f, a), X)


def field_map(cfg: CrossPlaneConfig, f: float, y, z, method="auto"):
    """|p| map helper: complex total field on the grid (y, z) (1-D arrays) at frequency f,
    shape (len(z), len(y)); points inside the shield body are NaN."""
    Yg, Zg = np.meshgrid(np.asarray(y, float), np.asarray(z, float))
    X = np.stack([Yg.ravel(), Zg.ravel()], 1)
    p = field_at(cfg, f, X, method=method)
    if cfg.has_shield:
        y0, y1 = cfg.body
        d = cfg.thickness(X[:, 0])
        inside = (X[:, 0] > y0) & (X[:, 0] < y1) & (np.abs(X[:, 1] - cfg.z_mid) < 0.5 * d)
        p[inside] = np.nan
    return p.reshape(Zg.shape)


_WORK = {}


def _work_one(args):
    i, f, method = args
    cfg, X, a = _WORK["cfg"], _WORK["X"], _WORK["A"][i]
    return i, field_at(cfg, f, X, method=method, a=a)


def matching_spectra(cfg: CrossPlaneConfig, X=None, workers=None, method="auto"):
    """(f, P) with P[n, m] the total-field spectrum [Pa s] at frequency f[n] at point m
    (default: the matching circle points of cfg.phis). f = 0 is the static lift limit."""
    f = cfg.frequencies()
    X = cfg.circle_points() if X is None else np.asarray(X, float).reshape(-1, 2)
    A = incident_coefficients(cfg, f)
    P = np.zeros((f.size, X.shape[0]), complex)
    r = np.hypot(X[:, 0], X[:, 1])
    P[0] = static_lift_coefficient(cfg.aircraft, cfg.beta) * (-X[:, 1] / r) / r
    jobs = [(i, f[i], method) for i in range(1, f.size)]
    # heaviest (high-f BEM) first
    jobs.sort(key=lambda j: -(j[1] if j[1] <= cfg.split_frequency() else 0.05 * j[1]))
    workers = cfg.workers if workers is None else workers
    _WORK.update(cfg=cfg, X=X, A=A)
    if workers > 1 and len(jobs) > 1:
        import multiprocessing as mp
        ctx = mp.get_context("fork")
        with ctx.Pool(workers, initializer=_limit_threads) as pool:
            for i, p in pool.imap_unordered(_work_one, jobs, chunksize=1):
                P[i] = p
    else:
        for job in jobs:
            i, p = _work_one(job)
            P[i] = p
    return f, P


def _limit_threads():
    os.environ["OPENBLAS_NUM_THREADS"] = "1"
    try:  # best effort: OpenBLAS bundled with numpy/scipy
        import ctypes
        import glob
        import numpy as _np
        for lib in glob.glob(os.path.join(os.path.dirname(_np.__file__), "..", "scipy_openblas*", "lib", "*.so*")) + \
                glob.glob(os.path.join(os.path.dirname(_np.__file__), "..", "numpy.libs", "*openblas*.so*")):
            try:
                L = ctypes.CDLL(lib)
            except OSError:
                continue
            for name in ("scipy_openblas_set_num_threads64_", "openblas_set_num_threads64_",
                         "openblas_set_num_threads", "scipy_openblas_set_num_threads"):
                if hasattr(L, name):
                    getattr(L, name)(1)
                    return
    except Exception:
        pass


def spectra_to_signatures(cfg: CrossPlaneConfig, f, P, R=None, retarded=True):
    """Band-limit (raised-cosine taper from taper_from*f_max to f_max) and synthesise
    time signals. Returns t (N,), p (N, M); tau = t - beta R/V if retarded."""
    df = f[1] - f[0]
    T_rec = 1.0 / df
    fs = cfg.fs
    N = int(round(T_rec * fs))
    N += N % 2
    dt = T_rec / N
    nb = N // 2 + 1
    X = np.zeros((nb, P.shape[1]), complex)
    n = min(nb, f.size)
    taper = np.ones(n)
    f0 = cfg.taper_from * cfg.f_max
    m = f[:n] > f0
    taper[m] = 0.5 * (1 + np.cos(np.pi * (f[:n][m] - f0) / (cfg.f_max - f0)))
    shift = cfg.t_pre
    if retarded:
        R = cfg.matching_radius() if R is None else R
        shift -= cfg.beta * R / cfg.V      # multiply by e^{+j w beta R/V} e^{-j w t_pre}
    X[:n] = P[:n] * (taper * np.exp(-2j * np.pi * f[:n] * shift))[:, None] / dt
    p = np.fft.irfft(X, N, axis=0)
    t = np.arange(N) * dt - cfg.t_pre
    return t, p


def matching_signatures(cfg: CrossPlaneConfig, workers=None, return_spectra=False):
    """dict phi_deg -> (t, p): total-field time signatures (retarded time) on the
    matching circle R_m = cfg.matching_radius() for each azimuth in cfg.phis
    (unshielded when cfg has no shield)."""
    f, P = matching_spectra(cfg, workers=workers)
    t, p = spectra_to_signatures(cfg, f, P)
    out = {float(ph): (t, p[:, i]) for i, ph in enumerate(cfg.phis)}
    if return_spectra:
        return out, (f, P)
    return out


def incident_vs_whitham(ac: Aircraft, r: float, phi: float = 0.0, f_max: float = 100.0, df: float = 0.25):
    """(tau, p_exact, p_whitham) at radius r, azimuth phi [rad]: the exact linear Hankel
    incident field (medium = flight altitude, f = 0 static-lift bin dropped) and
    ``source.near_field(ac, r, phi=phi)``, both band-limited identically (raised-cosine
    taper to f_max)."""
    from .source import near_field
    cfg = CrossPlaneConfig(ac, h=r, W=None, f_max=f_max, df=df, R_m=r, phis=[math.degrees(phi)],
                           atm=ac.atm, incident="linear")
    f = cfg.frequencies()
    A = incident_coefficients(cfg, f)
    w = 2 * np.pi * f
    k = w * math.sqrt(1 / cfg.c ** 2 - 1 / cfg.V ** 2)
    X = cfg.circle_points()
    P = np.array([multipole_field(A[i], k[i], X[:, 0], X[:, 1]) if f[i] > 0 else 0 * X[:, 0]
                  for i in range(f.size)])
    t, p = spectra_to_signatures(cfg, f, P)
    tw, pw = near_field(ac, r, 20000.0, y_max_factor=12.0, phi=phi)
    G = _signature_spectrum(tw, pw, f) * np.exp(-1j * w * cfg.beta * r / cfg.V)
    G[0] = 0.0
    _, pW = spectra_to_signatures(cfg, f, G[:, None])
    return t, p[:, 0], pW[:, 0]
