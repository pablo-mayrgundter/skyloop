"""Along-track-varying (fully 3-D) BoomShield: serrated edges and module joints.

The cross-plane model (:mod:`boomsim.crossplane`, contract ``DESIGN-3D.md``) needs a
shield that is invariant along the flight direction x. This module adds shields that
are *periodic* in x (period L_s): plan-view edge serrations (edge position y_e(x), depth
w_f, period Lambda) and transverse module joints (gap of width g every L_m).

API
---
Serration(Lam, w_f=150, shape="sawtooth"|"rc_teeth"|"random"|"linear_sawtooth",
          n_cells=None, seed=0)          plan-view fringe mask; .tau(u), .intervals(u),
                                         .coefficients(u, n) -> c_n(y)
Joint(L_m, g, kind="open"|"membranes")   transverse joint; .coefficients(n)
stack_difference(cfg, f, y, a=None, gas=None) -> D(y)   PO difference field below an
                                         infinite stack (= crossplane.po_solve's dfield)
rayleigh_2d(d, y, dy, za, kappa, X)      2-D Rayleigh integral of an aperture field
GroundRays(z_src, z_anchor, ...)          ray tables from an aperture plane to the ground
    .gain(sx, sy, f) -> (y_offset, G)    |P_ground|^2 = |D~(k_y)|^2 G (flux mapping)
fringe_orders(cfg, serr, f, rays, y_grid, ...) / joint_orders(cfg, joint, f, rays, ...)
    per-frequency x-averaged ground energy spectrum of the non-specular orders and the
    aperture-energy budget (evanescent / propagating / reaching the ground)
band_energies_signal(t, p), band_energies_esd(f, esd), pl_from_bands(e),
combined_pl(bands_spec, f, esd_scat)     Mark VII PL of specular + scattered energy

Method (physical optics + Floquet orders + ray-tube flux mapping)
-----------------------------------------------------------------
Ground-frame fields ~ exp(j(omega t - k_x x)); the incident has k_x0 = omega/V. Kirchhoff
(PO, as crossplane's PO path): just below the shield mid-plane the field is

    p_inc(y) + chi(x, y) D(y),   D(y) = sum_ky (T_stack(k_x0, k_y) - 1) a(k_y) e^{-j k_y y}

with chi the plan-view coverage (1 covered, 0 open) and T_stack the infinite-stack
transmission at the *incident* wavenumbers (local response; ``stack_difference``). For a
joint, (1 - chi_joint) (D_gap - D_stack) is added (D_gap = 0 for an open gap). With chi
periodic in x, chi = sum_n c_n(y) e^{-j 2 pi n x / L_s}, so the aperture field splits into
Floquet orders with k_x,n = omega/V + 2 pi n/L_s, each an x-plane wave and therefore an
exact 2-D problem in its own cross-plane with wavenumber kappa_n = sqrt(k0^2 - k_x,n^2):
propagating iff k_x,n^2 + k_y^2 < k0^2. c_0 = tau(y) (the x-averaged coverage), so the
n = 0 order is exactly crossplane's x-averaged fringe model; sum_{n != 0} |c_n|^2 =
tau (1 - tau) is what that model drops. c_n are exact (sums over covered x-intervals).

Ground mapping. Each order's aperture spectrum D~_n(k_y) = int d_n(y) e^{j k_y y} dy is
a set of plane waves with horizontal slowness (s_x, s_y) = (k_x,n, k_y)/omega, conserved in
the stratified windless atmosphere. A plane-wave band dk_y carries vertical energy flux
|D~|^2 k_z dk_y / (2 pi rho omega) (per unit x, exact angular-spectrum Parseval), conserved
along ray tubes; at the ground it is spread over dy_g = (I1 + s_y^2 I3) dk_y / omega
(y_g - y_c = s_y I1, I1 = int dz/s_z, I3 = int dz/s_z^3 from the aperture down). Hence

    |P_g|^2 = |D~|^2 * omega/(2 pi (I1 + s_y^2 I3)) * (rho_g s_z,src)/(rho_src s_z,g)
              * exp(-2 int alpha ds) * R^2 * rho_src/rho(z_a)

(R = 1.9 ground reflection; the last factor puts the uniform-air aperture field into the
pipeline's "source" anchor convention of ``propagate.circle_start_factors``). In uniform
air this is the stationary-phase far field |D~|^2 k_z^3/(2 pi dz kappa^2) (tested). The
launch point is the aperture centre y_c (fringe centre, or 0 for a joint): ground
positions are uncertain by +- half the aperture width (75 m fringe, 500 m joint). Rays with
s_x^2 + s_y^2 >= 1/c_ground^2 turn before the ground (ducted, or beyond lateral cut-off)
and are counted as "propagating, not reaching the ground". Upward (reflected) orders are
not followed.

Energy at the ground. At a fixed ground point the time signal of order n carries the
factor e^{-j 2 pi n x_g / L_s}; averaged over x_g (one period) the cross terms between
orders (and with the specular n = 0 field) vanish, so the x-averaged energy spectral
density is sum_n |P_g,n|^2. Within a 1/3-octave band the cross terms also average out at
fixed x (different orders arrive from different directions, so their relative phase
varies quickly with frequency); band energies are therefore added incoherently to the
specular (n = 0) ground signature's band energies, and Mark VII PL is computed from the
sum (``combined_pl``; identical to ``loudness.pldb`` when the scattered part is zero).

Limitations: PO (Kirchhoff truncation at tooth edges, local stack response, thickness
lumped at the mid-plane: teeth narrower than the gas-layer lateral decay length ~1/kappa
are less opaque than PO assumes); the scattered field is propagated linearly (see
``resteepening_check``); band limit of the incident spectrum and of the computation.
"""
from __future__ import annotations

import dataclasses
import math
from dataclasses import dataclass

import numpy as np
from scipy import special as sp

from . import crossplane as cp
from . import loudness as LD
from . import shield as SH
from .atmosphere import molar_humidity, std_atm
from .propagate import iso_absorption
from .proxy import METRIC_TAPER_S, RISE_SIGMA_S, rise_time_proxy

REFLECTION = 1.9


# =============================================================================== masks
@dataclass
class Serration:
    """Plan-view serrated edge of depth ``w_f`` and period ``Lam`` along x.

    Local fringe coordinate u in [0, 1]: u = 0 at the root (full coverage, joins the
    body), u = 1 at the tooth tips. All shapes except ``linear_sawtooth`` have the
    raised-cosine x-averaged coverage tau(u) = (1 + cos pi u)/2 of crossplane's
    ``fringe_profile="raised_cosine"`` (so their n = 0 order is the existing x-averaged
    fringe model):

    * ``sawtooth``: one-sided teeth, covered x in [k Lam, k Lam + tau Lam) (one flank
      straight across-track, the other curved so that tau is raised-cosine);
    * ``rc_teeth``: symmetric teeth centred in each period (both flanks curved);
    * ``random``: symmetric teeth, each cyclically shifted within its own period by a
      random fraction theta_k (uniform), over a supercell of ``n_cells`` periods (default
      8; the order spacing is then 2 pi/(n_cells Lam)). The per-period coverage is
      exactly tau, so there is no root leak.
    * ``linear_sawtooth``: straight-flanked one-sided teeth, tau(u) = 1 - u.
    """
    Lam: float
    w_f: float = 150.0
    shape: str = "sawtooth"
    n_cells: int | None = None
    seed: int = 0

    def __post_init__(self):
        if self.shape not in ("sawtooth", "rc_teeth", "random", "linear_sawtooth"):
            raise ValueError(f"unknown serration shape {self.shape!r}")
        if self.n_cells is None:
            self.n_cells = 8 if self.shape == "random" else 1
        rng = np.random.default_rng(self.seed)
        self.theta = rng.random(self.n_cells) if self.shape == "random" else np.zeros(self.n_cells)

    @property
    def L_s(self) -> float:
        """Supercell length (period of the mask) [m]."""
        return self.Lam * self.n_cells

    def tau(self, u):
        u = np.clip(np.asarray(u, float), 0.0, 1.0)
        if self.shape == "linear_sawtooth":
            return 1.0 - u
        return 0.5 * (1.0 + np.cos(np.pi * u))

    def intervals(self, u):
        """Covered x-intervals in [0, L_s): arrays (a, b) of shape (len(u), 2 n_cells)
        (empty pieces have a = b)."""
        u = np.atleast_1d(np.asarray(u, float))
        L = self.Lam
        wdt = self.tau(u)[:, None] * L                       # (Nu, 1)
        k = np.arange(self.n_cells)[None, :] * L             # cell starts (1, N)
        if self.shape in ("sawtooth", "linear_sawtooth"):
            lo, hi = k + 0 * wdt, k + wdt
        elif self.shape == "rc_teeth":
            lo, hi = k + 0.5 * (L - wdt), k + 0.5 * (L + wdt)
        else:
            c = k + self.theta[None, :] * L
            lo, hi = c - 0.5 * wdt, c + 0.5 * wdt
        a1, b1 = np.maximum(lo, k), np.minimum(hi, k + L)
        # wrapped remainder (random shape only)
        a2 = np.where(lo < k, lo + L, np.where(hi > k + L, k, 0.0))
        b2 = np.where(lo < k, k + L, np.where(hi > k + L, hi - L, 0.0))
        a = np.concatenate([a1, a2], axis=1)
        b = np.concatenate([b1, np.maximum(b2, a2)], axis=1)
        return a, b

    def coefficients(self, u, n):
        """c_n(u), shape (len(n), len(u)): chi = sum_n c_n e^{-j 2 pi n x / L_s}."""
        a, b = self.intervals(u)
        return interval_coefficients(a, b, self.L_s, n)


def interval_coefficients(a, b, L_s, n):
    """Exact Fourier coefficients c_n = (1/L_s) int chi e^{+j 2 pi n x/L_s} dx of the
    indicator of the union of intervals [a, b) (arrays (Nu, P)); returns (len(n), Nu)."""
    n = np.atleast_1d(np.asarray(n))
    kap = 2 * np.pi * n / L_s
    out = np.zeros((n.size, a.shape[0]), complex)
    nz = kap != 0
    kz = kap[nz][:, None]
    for p in range(a.shape[1]):
        ap, bp = a[:, p][None, :], b[:, p][None, :]
        out[nz] += (np.exp(1j * kz * bp) - np.exp(1j * kz * ap)) / (1j * kz)
        out[~nz] += (bp - ap)
    return out / L_s


@dataclass
class Joint:
    """Transverse module joint: a strip of width ``g`` [m] every ``L_m`` [m] along x across
    the whole shield width, where the gas layer is interrupted. ``kind``: "open" (nothing)
    or "membranes" (the two face membranes with ambient air between them, same spacing)."""
    L_m: float
    g: float
    kind: str = "open"

    def coefficients(self, n):
        """Fourier coefficients of the gap indicator (gap over [0, g))."""
        a = np.zeros((1, 1))
        b = np.full((1, 1), float(self.g))
        return interval_coefficients(a, b, self.L_m, n)[:, 0]


# =============================================================================== PO aperture
def stack_difference(cfg: cp.CrossPlaneConfig, f: float, y, a=None, gas=None, d=None,
                     theta_factor: float = 1.0):
    """PO difference field D(y) = (field below an infinite stack) - (incident) on the
    stack mid-plane z = cfg.z_mid at the points y [m], frequency f > 0 (the spectral sum
    of ``crossplane.po_solve``: exact plane-wave spectrum of the incident, propagating and
    evanescent branches, times T(k_x0, |k_y|) - 1 of ``crossplane._stack_TR``). ``gas``
    replaces the fill (e.g. ambient air for a "membranes only" joint); ``d`` the
    thickness (default cfg.T)."""
    y = np.atleast_1d(np.asarray(y, float))
    w, k, _ = cp._wavenumbers(cfg, f)
    if a is None:
        a = cp.incident_coefficients(cfg, [f])[0]
    za = cfg.z_mid
    rho_max = float(np.max(np.hypot(abs(za), np.abs(y))))
    n_th = int(theta_factor * (0.7 * k * rho_max * np.pi / 2 + 60))
    th, wth = cp._gl(n_th, -np.pi / 2, np.pi / 2)
    ky_p = k * np.sin(th)
    wp = (1 / np.pi) * np.exp(-1j * k * np.cos(th) * abs(za)) * cp._spectral_weight(a, th) * wth
    u_max = math.asinh(40.0 / (k * abs(za)))
    uu, wu = cp._gl(48, 0.0, u_max)
    ky_e = k * np.cosh(uu)
    we = (1 / np.pi) * np.exp(-k * np.sinh(uu) * abs(za)) * cp._spectral_weight(a, np.pi / 2 + 1j * uu) * 1j * wu
    ky_all = np.concatenate([ky_p, ky_e, -ky_e])
    w_all = np.concatenate([wp, we, we])
    c2 = cfg if gas is None else dataclasses.replace(cfg, gas=gas)
    Tm, _ = cp._stack_TR(c2, w, np.abs(ky_all), cfg.T if d is None else d, 0.0)
    out = np.empty(y.size, complex)
    v = (Tm - 1.0) * w_all
    for i0 in range(0, y.size, 2048):
        out[i0:i0 + 2048] = np.exp(-1j * np.outer(y[i0:i0 + 2048], ky_all)) @ v
    return out


def rayleigh_2d(d, y, dy, za, kappa, X):
    """Field below/above an aperture d(y) (uniform samples y, spacing dy) on the plane z =
    za, at points X (M, 2): (-j kappa/2) sum d H1^(2)(kappa rho) |z - za|/rho dy (the
    radiation of crossplane.po_field without the incident; kappa may be complex with
    Im <= 0)."""
    X = np.asarray(X, float).reshape(-1, 2)
    out = np.zeros(X.shape[0], complex)
    for i0 in range(0, X.shape[0], 256):
        xx = X[i0:i0 + 256]
        dyv = xx[:, None, 0] - y[None, :]
        dz = np.abs(xx[:, None, 1] - za)
        rr = np.hypot(dyv, dz)
        out[i0:i0 + 256] = (-0.5j * kappa) * np.sum(d[None, :] * sp.hankel2(1, kappa * rr) * dz / rr, axis=1) * dy
    return out


# =============================================================================== rays
class GroundRays:
    """Ray tables from an aperture at altitude ``z_src`` down to the ground (windless,
    stratified ``atmosphere``; default std_atm), as functions of the horizontal slowness
    magnitude s_h (the medium is horizontally isotropic, so the horizontal ray offset is
    (s_x, s_y) I1(s_h)). Tabulated on the ground vertical slowness s_zg = sqrt(1/c_g^2 -
    s_h^2) (refined towards grazing); the integrals use RayPath's exact cell rule for s_z^2
    linear in z. Assumes the ground has the largest sound speed below z_src (true for the
    standard atmosphere), so a ray reaches the ground iff s_h < 1/c_g.

    ``gain(sx, sy, f)`` -> (dy, G): lateral ground offset from the launch point and the
    factor G with |P_g|^2 = |D~(k_y)|^2 G (module docstring), including absorption
    (ISO 9613-1 with ``molar_humidity``, as propagate), the ground reflection factor and
    the source-anchor density factor rho(z_src)/rho(z_anchor)."""

    def __init__(self, z_src, z_anchor=None, atmosphere=None, dz=5.0, z_ground=0.0,
                 n_tab=500, reflection=REFLECTION):
        self.z_src = float(z_src)
        atm = std_atm if atmosphere is None else atmosphere
        n = int(math.ceil((self.z_src - z_ground) / dz)) + 1
        z = np.linspace(z_ground, self.z_src, n)
        a = atm(z)
        self.z = z
        self.c = np.asarray(a.c, float) * np.ones_like(z)
        self.rho = np.asarray(a.rho, float) * np.ones_like(z)
        zm = 0.5 * (z[1:] + z[:-1])
        am = atm(zm)
        self.zm = zm
        self.Tm = np.asarray(am.T, float) * np.ones_like(zm)
        self.pm = np.asarray(am.p, float) * np.ones_like(zm)
        self.hm = molar_humidity(zm, atm=am) * np.ones_like(zm)
        self.cm = 0.5 * (self.c[1:] + self.c[:-1])
        if np.max(self.c) > self.c[0] * (1 + 1e-9):
            raise ValueError("GroundRays assumes the ground has the largest sound speed")
        self.c_g, self.c_s = float(self.c[0]), float(self.c[-1])
        self.rho_g, self.rho_s = float(self.rho[0]), float(self.rho[-1])
        za = self.z_src if z_anchor is None else float(z_anchor)
        self.anchor = self.rho_s / float(np.asarray(atm(np.array([za])).rho).ravel()[0])
        self.reflection = float(reflection)
        smax = 1.0 / self.c_g
        e = np.concatenate([np.geomspace(1e-6, 2e-2, n_tab // 3, endpoint=False),
                            np.linspace(2e-2, 1.0, n_tab - n_tab // 3)])
        szg = smax * e                                       # ground vertical slowness
        sh2 = smax ** 2 - szg ** 2
        self.szg_tab = szg
        q = 1.0 / self.c[None, :] ** 2 - sh2[:, None]        # s_z^2 on the z grid
        r = np.sqrt(np.maximum(q, 0.0))
        r[:, 0] = szg
        dzc = np.diff(z)[None, :]
        rs = r[:, :-1] + r[:, 1:]
        dI1 = 2.0 * dzc / rs
        dI3 = 2.0 * dzc / (r[:, :-1] * r[:, 1:] * rs)
        self.I1_tab = dI1.sum(1)
        self.I3s_tab = dI3.sum(1) * szg                      # s_zg * I3 (smooth)
        self.W = dI1 / self.cm[None, :]                      # arc length per cell
        self.szs_tab = np.sqrt(1.0 / self.c_s ** 2 - sh2)
        self._abs_cache = {}

    def _absorb(self, f):
        key = float(f)
        if key not in self._abs_cache:
            al = iso_absorption(f, self.Tm, self.pm, self.hm)
            self._abs_cache[key] = self.W @ al
        return self._abs_cache[key]

    def gain(self, sx, sy, f):
        """(dy, G, ok): arrays broadcast from sx, sy; ok = reaches the ground."""
        sx, sy = np.broadcast_arrays(np.asarray(sx, float), np.asarray(sy, float))
        sh2 = sx ** 2 + sy ** 2
        szg2 = 1.0 / self.c_g ** 2 - sh2
        ok = szg2 > (self.szg_tab[0]) ** 2
        szg = np.sqrt(np.where(ok, szg2, self.szg_tab[0] ** 2))
        I1 = np.interp(szg, self.szg_tab, self.I1_tab)
        I3 = np.interp(szg, self.szg_tab, self.I3s_tab) / szg
        Ab = np.interp(szg, self.szg_tab, self._absorb(f))
        szs = np.sqrt(np.maximum(1.0 / self.c_s ** 2 - sh2, 0.0))
        w = 2 * np.pi * f
        I13 = I1 + sy ** 2 * I3
        G = (w / (2 * np.pi * I13)) * (self.rho_g * szs) / (self.rho_s * szg) \
            * np.exp(-2.0 * Ab) * self.reflection ** 2 * self.anchor
        G = np.where(ok, G, 0.0)
        return sy * I1, G, ok

    def offset(self, sx, sy):
        """Lateral ground offset s_y I1 (NaN if the ray does not reach the ground)."""
        dy, _, ok = self.gain(sx, sy, 1.0)
        return np.where(ok, dy, np.nan)


# =============================================================================== orders
def _deposit(acc, y_grid, y_s, v_s):
    """Linear interpolation of the samples (y_s increasing, values v_s) onto y_grid,
    added to acc (zero outside the sample range)."""
    if y_s.size < 2:
        return
    m = (y_grid >= y_s[0]) & (y_grid <= y_s[-1])
    if np.any(m):
        acc[m] += np.interp(y_grid[m], y_s, v_s)


def _ky_grid(k0, width, oversample=4.0):
    dk = 2 * np.pi / (oversample * width)
    n = int(math.ceil(k0 / dk))
    return np.arange(-n, n + 1) * dk, dk


def fringe_orders(cfg, serr: Serration, f, rays: GroundRays, y_grid, a=None, block=None,
                  ky_oversample=4.0):
    """Non-specular orders of the right-hand serrated fringe at frequency f.

    The fringe is centred on the nominal edge y_r = cfg.y_right (tau = 1/2 there), root
    at y_r - w_f/2. ``block`` = dict(u, y, dy, n, C) (precomputed c_n on a y grid; see
    :func:`fringe_block`), computed here if None. Returns dict with
    ``E``: x-averaged |P_g|^2 [(Pa s)^2] of all orders n != 0 on ``y_grid`` (right fringe
    only; the left fringe is the mirror image, E_left(y) = E(-y)),
    ``E0``: same for the n = 0 order (linear flux mapping; diagnostic only),
    budget (aperture energies int |d|^2 dy [(Pa s)^2 m]): ``ap_spec`` (n = 0, fringe
    only), ``ap_scat`` (all n != 0, exact Parseval: int tau (1 - tau) |D|^2), ``ap_prop``
    (n != 0 propagating at the shield), ``ap_ground`` (n != 0 reaching the ground), and
    ``ap_D`` = int |D|^2 dy over the fringe (reference: a fully covered fringe);
    ``y_mean``, ``sx_mean``, ``sy_mean``: energy-weighted landing offset and slowness of
    the ground-reaching scattered part."""
    w = 2 * np.pi * f
    k0 = w / cfg.c
    kx0 = w / cfg.V
    if block is None:
        block = fringe_block(serr, f, cfg.c, cfg.V)
    u, yl, dy, n, C = block["u"], block["y"], block["dy"], block["n"], block["C"]
    yc = float(cfg.y_right)
    D = stack_difference(cfg, f, yc + yl, a=a)
    tau = serr.tau(u)
    out = dict(f=f, ap_D=float(np.sum(np.abs(D) ** 2) * dy),
               ap_spec=float(np.sum(np.abs(tau * D) ** 2) * dy),
               ap_scat=float(np.sum(tau * (1 - tau) * np.abs(D) ** 2) * dy))
    kxn = kx0 + 2 * np.pi * n / serr.L_s
    sel = (np.abs(kxn) < k0)
    ky, dk = _ky_grid(k0, serr.w_f, ky_oversample)
    E = np.exp(1j * np.outer(yl, ky)) * dy                  # (Ny, Nky)
    M = C[sel] * D[None, :]
    S = M @ E                                                # (Nn_sel, Nky)
    nn = n[sel]
    kx = kxn[sel]
    prop = (kx[:, None] ** 2 + ky[None, :] ** 2) < k0 ** 2
    P2 = np.abs(S) ** 2
    scat = (nn != 0)[:, None]
    out["ap_prop"] = float(np.sum(P2 * (prop & scat)) * dk / (2 * np.pi))
    dyo, G, ok = rays.gain(kx[:, None] / w, ky[None, :] / w, f)
    gnd = ok & prop
    out["ap_ground"] = float(np.sum(P2 * (gnd & scat)) * dk / (2 * np.pi))
    Es = np.zeros(y_grid.size)
    E0 = np.zeros(y_grid.size)
    V = P2 * G
    for i in range(nn.size):
        m = gnd[i]
        if np.count_nonzero(m) < 2:
            continue
        _deposit(E0 if nn[i] == 0 else Es, y_grid, yc + dyo[i, m], V[i, m])
    out["E"] = Es
    out["E0"] = E0
    # energy-weighted mean landing offset and slowness of the ground-reaching scattered part
    wgt = P2 * (gnd & scat) * dk / (2 * np.pi)
    tot = wgt.sum()
    if tot > 0:
        out["y_mean"] = float(np.sum(wgt * (yc + np.where(gnd, dyo, 0.0))) / tot)
        out["sx_mean"] = float(np.sum(wgt * (kx[:, None] / w)) / tot)
        out["sy_mean"] = float(np.sum(wgt * (ky[None, :] / w)) / tot)
    else:
        out["y_mean"] = out["sx_mean"] = out["sy_mean"] = float("nan")
    return out


def fringe_block(serr: Serration, f_hi, c, V, dy=None):
    """c_n(y) on a fringe grid adequate for frequencies up to f_hi: all orders with
    |k_x,n| < k0(f_hi), y spacing resolving the orders' lateral phase rate (2 pi n_eff
    tau' for one-sided teeth, pi n_eff tau' for symmetric/random ones, n_eff = n/n_cells)
    plus k0, with a 1.4x margin (no aliasing into the propagating k_y range)."""
    w = 2 * np.pi * f_hi
    k0 = w / c
    kx0 = w / V
    Ls = serr.L_s
    n_lo = int(math.floor((-k0 - kx0) * Ls / (2 * np.pi)))
    n_hi = int(math.ceil((k0 - kx0) * Ls / (2 * np.pi)))
    n = np.arange(n_lo, n_hi + 1)
    n_eff = max(abs(n_lo), abs(n_hi)) / serr.n_cells
    # lateral phase rate of c_n(y): one-sided teeth 2 pi n tau', symmetric/random pi n tau'
    tau_p = 1.0 / serr.w_f if serr.shape == "linear_sawtooth" else np.pi / (2 * serr.w_f)
    fac = 2.0 if serr.shape in ("sawtooth", "linear_sawtooth") else 1.0
    k_content = k0 + fac * np.pi * n_eff * tau_p
    if dy is None:
        dy = min(np.pi / (1.4 * k_content), serr.w_f / 64)
    ny = int(math.ceil(serr.w_f / dy))
    dy = serr.w_f / ny
    u = (np.arange(ny) + 0.5) / ny
    yl = (u - 0.5) * serr.w_f                                # local y (fringe centre 0)
    C = serr.coefficients(u, n)
    return dict(u=u, y=yl, dy=dy, n=n, C=C, f_hi=f_hi)


def joint_orders(cfg, joints, f, rays: GroundRays, y_grid, a=None, dy_lambda=8.0,
                 ky_oversample=2.0, air=None):
    """Non-specular orders of transverse joints (gap across the whole sharp-edged shield
    [y_left, y_right]) at frequency f. ``joints``: a Joint or a list of Joints of the same
    ``kind`` (they share the aperture difference dD = D_gap - D_stack). Returns a dict (or
    a list of dicts) with ``E`` (x-averaged |P_g|^2 of all n != 0 on y_grid; launch point
    y = 0, symmetric) and the aperture-energy budget as :func:`fringe_orders`
    (``ap_scat`` = (g/L)(1 - g/L) int |dD|^2, ``ap_spec`` = (g/L)^2 int |dD|^2: the n = 0
    coherent leak, which is (g/L) dD, see :func:`joint_aperture`)."""
    single = isinstance(joints, Joint)
    joints = [joints] if single else list(joints)
    kinds = {j.kind for j in joints}
    if len(kinds) != 1:
        raise ValueError("joints must share kind")
    w = 2 * np.pi * f
    k0 = w / cfg.c
    kx0 = w / cfg.V
    y, dy, dD = joint_aperture(cfg, f, kinds.pop(), a=a, dy_lambda=dy_lambda, air=air)
    apD = float(np.sum(np.abs(dD) ** 2) * dy)
    ky, dk = _ky_grid(k0, cfg.y_right - cfg.y_left, ky_oversample)
    S0 = np.empty(ky.size, complex)
    for i0 in range(0, ky.size, 1024):
        S0[i0:i0 + 1024] = np.exp(1j * np.outer(ky[i0:i0 + 1024], y)) @ dD * dy
    P2 = np.abs(S0) ** 2                                     # |S0(k_y)|^2
    res = [None] * len(joints)
    for L in sorted({jt.L_m for jt in joints}):
        n_lo = int(math.floor((-k0 - kx0) * L / (2 * np.pi)))
        n_hi = int(math.ceil((k0 - kx0) * L / (2 * np.pi)))
        n = np.arange(n_lo, n_hi + 1)
        n = n[n != 0]
        kx = kx0 + 2 * np.pi * n / L
        base = np.zeros((n.size, y_grid.size))              # per-order ground |P|^2 (c_n = 1)
        prop_e = np.zeros(n.size)
        gnd_e = np.zeros(n.size)
        for i in range(n.size):
            prop = (kx[i] ** 2 + ky ** 2) < k0 ** 2
            if not np.any(prop):
                continue
            prop_e[i] = np.sum(P2[prop]) * dk / (2 * np.pi)
            dyo, G, ok = rays.gain(kx[i] / w, ky[prop] / w, f)
            if np.count_nonzero(ok) < 2:
                continue
            gnd_e[i] = np.sum(P2[prop][ok]) * dk / (2 * np.pi)
            _deposit(base[i], y_grid, dyo[ok], (P2[prop] * G)[ok])
        for jj, jt in enumerate(joints):
            if jt.L_m != L:
                continue
            q = jt.g / jt.L_m
            c2 = np.abs(jt.coefficients(n)) ** 2
            res[jj] = dict(f=f, ap_D=apD, ap_spec=q ** 2 * apD, ap_scat=q * (1 - q) * apD,
                           ap_prop=float(c2 @ prop_e), ap_ground=float(c2 @ gnd_e), E=c2 @ base)
    return res[0] if single else res


def joint_aperture(cfg, f, kind, a=None, dy_lambda=8.0, air=None):
    """(y, dy, dD): PO aperture difference D_gap - D_stack across the sharp shield at f
    (D_gap = 0 for an "open" joint; ambient-air fill between the face membranes for
    "membranes"). The n = 0 (x-averaged) leak of a joint of open fraction q = g/L_m is
    q dD, coherent with the shield's own field."""
    k0 = 2 * np.pi * f / cfg.c
    y0, y1 = cfg.y_left, cfg.y_right
    lam = 2 * np.pi / k0
    ny = max(64, int(math.ceil((y1 - y0) / (lam / dy_lambda))))
    dy = (y1 - y0) / ny
    y = y0 + (np.arange(ny) + 0.5) * dy
    Ds = stack_difference(cfg, f, y, a=a)
    if kind == "open":
        dD = -Ds
    elif kind == "membranes":
        if air is None:
            air = SH.gas("air", float(cfg.atm.T), float(cfg.atm.p))
        dD = stack_difference(cfg, f, y, a=a, gas=air) - Ds
    else:
        raise ValueError("joint kind must be 'open' or 'membranes'")
    return y, dy, dD


# =============================================================================== loudness
def band_energies_signal(t, p, taper_s=METRIC_TAPER_S):
    """Mark VII 1/3-octave band energies [Pa^2 s] of a signature (loudness.pldb steps 1-2)."""
    t, p, dt = LD._check_signal(t, p)
    p = LD._taper(p, dt, taper_s)
    f, esd = LD._esd(p, dt)
    e = LD._band_energies(f, esd, LD._BAND_LO, LD._BAND_HI)
    e[LD._BAND_LO >= f[-1]] = 0.0
    return e


def band_energies_esd(f, esd):
    """Band energies of a one-sided ESD [Pa^2 s/Hz] sampled at increasing f (any
    spacing; zero outside [f[0], f[-1]])."""
    f = np.asarray(f, float)
    esd = np.asarray(esd, float)
    e = LD._band_energies(f, esd, LD._BAND_LO, LD._BAND_HI)
    e[LD._BAND_HI <= f[0]] = 0.0
    e[LD._BAND_LO >= f[-1]] = 0.0
    return e


def pl_from_bands(e):
    """Mark VII PL [PLdB] from 42 band energies (loudness.pldb steps 3-7)."""
    with np.errstate(divide="ignore"):
        L = 10.0 * np.log10(np.asarray(e, float) / (LD.T_EAR * LD.P_REF ** 2)) - 3.0
    St, _ = LD._total_sones(LD._equivalent_level(L))
    if St <= 0:
        return float("-inf")
    return float(32.0 + 9.0 * np.log2(St))


def signature_bands(t, p, taper_s=METRIC_TAPER_S, sigma_s=RISE_SIGMA_S):
    """(raw, proxy) band energies of a ground signature (proxy = rise-time proxy)."""
    return (band_energies_signal(t, p, taper_s),
            band_energies_signal(t, rise_time_proxy(t, p, sigma_s), taper_s))


def combined_pl(bands_spec, f, E_scat, sigma_s=RISE_SIGMA_S):
    """(PL raw, PL proxy) of specular band energies (raw, proxy) plus the x-averaged
    scattered ground spectrum E_scat(f) = sum_n |P_g,n|^2 [(Pa s)^2] (one-sided ESD
    2 E_scat; the proxy applies the Gaussian smoothing's |H|^2 = exp(-(2 pi f sigma)^2))."""
    esd = 2.0 * np.asarray(E_scat, float)
    eb = band_energies_esd(f, esd)
    ebp = band_energies_esd(f, esd * np.exp(-(2 * np.pi * np.asarray(f) * sigma_s) ** 2))
    return pl_from_bands(bands_spec[0] + eb), pl_from_bands(bands_spec[1] + ebp)


def crossplane_taper(f, f_max=100.0, taper_from=0.75):
    """crossplane's band-limit taper (1 below taper_from*f_max, raised cosine to f_max)."""
    f = np.asarray(f, float)
    f0 = taper_from * f_max
    w = np.where(f > f0, 0.5 * (1 + np.cos(np.pi * np.clip((f - f0) / (f_max - f0), 0, 1))), 1.0)
    return np.where(f >= f_max, 0.0, w)


# =============================================================================== nonlinearity check
def resteepening_check(rays_atm_alt, sx, sy, f_b, p_ground, z_src, z_start, beta=1.2,
                       reflection=REFLECTION):
    """Nonlinear distortion parameter Sigma = int beta omega p(s) /(rho c^3) ds of a tone of
    frequency f_b whose free-field amplitude at the ground is p_ground/reflection, along the
    ray (s_x, s_y) from altitude z_start (below the aperture z_src) to the ground, amplitude
    from the linear WKB invariant (propagate.RayPath with altitude z_src, V = 1/s_x).
    Sigma < ~0.3: linear propagation is adequate (harmonic generation < ~2 %); Sigma >= 1:
    a shock forms. Returns Sigma."""
    from .propagate import RayPath
    V = 1.0 / max(abs(sx), 1e-12)
    c_a = float(std_atm(np.array([z_src])).c[0])
    cmu = math.sqrt(max(1.0 - (c_a / V) ** 2, 1e-12))
    sphi = float(np.clip(sy * c_a / cmu, -0.999999, 0.999999))
    ray = RayPath(z_src, V, phi=math.asin(sphi), allow_cutoff=True)
    z = np.linspace(max(ray.z_bottom, 0.0) + 1.0, z_start, 2000)
    amp = (p_ground / reflection) * ray.amplitude_ratio(z[0], z)
    a = ray.atm(z)
    ds = np.abs(np.gradient(z)) / np.maximum(ray.nz(z), 1e-6)
    return float(np.sum(beta * 2 * np.pi * f_b * amp / (a.rho * a.c ** 3) * ds))
