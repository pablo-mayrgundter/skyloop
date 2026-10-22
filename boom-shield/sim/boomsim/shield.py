"""Layered transfer-matrix model of a BoomShield + lift budget.

Geometry and conventions
------------------------
The shield is a horizontally infinite stack of layers, listed **top to bottom**
(the boom arrives from above). Ambient air (``atm.rho``, ``atm.c``) bounds the
stack on both sides. Linear acoustics, harmonic fields

    p(x, zeta, t) = P exp(j(omega t - k_x x - k_z zeta)),

where ``zeta`` is depth measured *downward* (direction of incidence) and
``k_x = omega / V`` is fixed in every layer (Snell; V = aircraft speed = trace
velocity of the boom on any horizontal surface).

Normal-velocity (positive downward) impedance of a wave travelling in +zeta is
``Z = rho omega / k_z``. The state vector is ``(p, v_zeta)``; every element has a
2x2 matrix ``M`` with ``[p, v]_top = M [p, v]_bottom``:

* Membrane / Screen (series impedance Z_s):  [[1, Z_s], [0, 1]]
* Gas layer of thickness d:  [[cos k_z d, j Z sin k_z d], [j sin(k_z d)/Z, cos k_z d]]
  (optionally with face-to-face tie springs ``tie_K`` in parallel, see
  :func:`_tie_layer_scaled`)

Branch choice. ``k_z = sqrt((omega/c)^2 - k_x^2)`` with Im(k_z) <= 0. For a gas
"subsonic" relative to the trace speed (V < c_gas) this gives ``k_z = -j kappa``,
``kappa = omega sqrt(1/V^2 - 1/c_gas^2)``, so ``exp(-j k_z zeta) = exp(-kappa zeta)``:
the field decays away from the incident side (with our e^{+j omega t} convention
this is the physical, bounded branch). Note the layer matrix depends only on
``k_z^2`` (cos is even; ``Z sin = rho omega sin(k_z d)/k_z`` and
``sin/Z = k_z sin(k_z d)/(rho omega)`` are even), so inside a finite layer the
branch is immaterial; it matters only for the semi-infinite bounding air, where
V > c makes k_z real and positive (outgoing, downward-travelling wave).

Numerical stability. cos/sin of ``k_z d`` grow like exp(|Im k_z d|)/2, which
overflows for kappa d > ~700 (d up to 300 m, f up to 10 kHz gives kappa d ~ 1e4).
Each layer matrix is therefore computed pre-scaled by exp(-|Im k_z d|) and the
scale is accumulated as a log factor ``L``; T carries exp(-L) explicitly while R
is scale-invariant.

Evanescent layers and causality. With V < c_gas and negligible membrane mass,
all impedances are frequency independent, so T(omega) ~ C exp(-|omega| a d),
a = sqrt(1/V^2 - 1/c_gas^2): a Poisson-kernel (Lorentzian, width a d) smoothing
of the signature. It is non-causal in t at a fixed x - this is physical: for a
steady supersonic flight the field is stationary in the aircraft frame, so "time"
here is really the coordinate x - V t, and a subsonic (elliptic) layer smooths
in that coordinate both ways.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from typing import Sequence

import numpy as np

R_UNIVERSAL = 8.314462618  # J/(mol K)

_GASES = {
    # kind: (molar mass [kg/mol], gamma)
    "air": (28.964e-3, 1.4),
    "helium": (4.0026e-3, 5.0 / 3.0),
    "hydrogen": (2.016e-3, 1.405),
}


# --------------------------------------------------------------------------- gas
@dataclass(frozen=True)
class Gas:
    kind: str
    T: float      # K
    p: float      # Pa (absolute pressure inside the layer)
    M: float      # kg/mol
    gamma: float

    @property
    def rho(self) -> float:
        return self.p * self.M / (R_UNIVERSAL * self.T)

    @property
    def c(self) -> float:
        return float(np.sqrt(self.gamma * R_UNIVERSAL * self.T / self.M))


def gas(kind: str, T: float, p: float) -> Gas:
    """Ideal gas ``kind`` ("air", "helium", "hydrogen") at temperature T [K] and
    absolute pressure p [Pa]. For partially evacuated air pass a p below ambient:
    density scales with p, sound speed (depends on T only) is unchanged."""
    kind = kind.lower()
    if kind not in _GASES:
        raise ValueError(f"unknown gas {kind!r}; choose from {sorted(_GASES)}")
    M, gamma = _GASES[kind]
    return Gas(kind, float(T), float(p), M, gamma)


def gas_mixture(fractions: dict, T: float, p: float) -> Gas:
    """Ideal-gas mixture at T [K], absolute pressure p [Pa].

    ``fractions`` maps gas kind (see :func:`gas`) to **mole** (= volume) fraction,
    e.g. ``{"helium": 0.95, "air": 0.05}``; fractions are normalised to sum to 1
    (they must already sum to 1 within 1e-6). Ideal mixing:

        M = sum x_i M_i,  c_v = sum x_i R/(gamma_i - 1),  gamma = (c_v + R)/c_v

    (molar heat capacities add; frozen composition, no relaxation of the mixture).
    """
    if not fractions:
        raise ValueError("empty mixture")
    xs = {k.lower(): float(v) for k, v in fractions.items()}
    tot = sum(xs.values())
    if any(v < 0 for v in xs.values()) or abs(tot - 1.0) > 1e-6:
        raise ValueError(f"mole fractions must be >= 0 and sum to 1 (got {tot})")
    M = cv = 0.0
    for k, x in xs.items():
        if k not in _GASES:
            raise ValueError(f"unknown gas {k!r}; choose from {sorted(_GASES)}")
        Mi, gi = _GASES[k]
        M += x / tot * Mi
        cv += x / tot * R_UNIVERSAL / (gi - 1.0)
    label = "+".join(f"{k}{x / tot:.3g}" for k, x in xs.items() if x > 0)
    return Gas(label, float(T), float(p), M, (cv + R_UNIVERSAL) / cv)


# ---------------------------------------------------------------------- elements
@dataclass
class Membrane:
    """Limp sheet: areal mass m [kg/m^2], in-plane tension N [N/m], loss factor eta.

    m w_tt - N(1+j eta) w_xx = delta_p  =>
    Z_m = j omega m (1 - (1+j eta) N k_x^2 / (m omega^2)) = j omega m (1 - (1+j eta) c_m^2/V^2).
    """
    m: float
    N: float = 0.0
    eta: float = 0.0

    @property
    def areal_mass(self) -> float:
        return self.m

    def impedance(self, omega, sx):
        omega = np.asarray(omega, dtype=float)
        # j omega m - j (1+j eta) N k_x^2 / omega, with k_x = omega sx
        return 1j * omega * self.m - 1j * (1 + 1j * self.eta) * self.N * omega * sx ** 2


@dataclass
class Screen:
    """Flow-resistive sheet: series impedance R [Pa s/m] (real) plus optional mass m."""
    R: float
    m: float = 0.0

    @property
    def areal_mass(self) -> float:
        return self.m

    def impedance(self, omega, sx):
        omega = np.asarray(omega, dtype=float)
        return self.R + 1j * omega * self.m + 0 * omega


@dataclass
class GasLayer:
    """Slab of gas ``gas`` (see :func:`gas`, :func:`gas_mixture`) of thickness d [m].

    ``tie_K`` [Pa/m]: smeared spring stiffness per unit area tying the layer's top
    and bottom boundaries (tendons/drop-threads/webs between the two faces, smeared
    over the footprint: tension per tie / tie spacing^2), acting in parallel with the
    gas. The spring is massless and lossless; only the thickness change w_bot - w_top
    loads it. See :func:`_tie_layer_scaled` for the derivation.
    """
    gas: Gas
    d: float
    tie_K: float = 0.0

    @property
    def areal_mass(self) -> float:  # buoyancy is handled in lift_margin
        return 0.0


@dataclass
class Stack:
    """Layers listed top (incident side) to bottom."""
    layers: Sequence = field(default_factory=list)

    @property
    def thickness(self) -> float:
        return float(sum(l.d for l in self.layers if isinstance(l, GasLayer)))


# ----------------------------------------------------------------------- physics
def _kz(omega, c, sx):
    """Normal wavenumber with Im(k_z) <= 0 (decaying/outgoing in +zeta)."""
    k2 = (omega ** 2) * (1.0 / c ** 2 - sx ** 2) + 0j
    kz = np.sqrt(k2)
    return np.where(kz.imag > 0, -kz, kz)


def _slowness(V) -> float:
    return 0.0 if np.isinf(V) else 1.0 / float(V)


def _layer_matrix_scaled(layer: GasLayer, omega, sx):
    """Gas-layer matrix scaled by exp(-|Im k_z d|). Returns (a, b, c, d, logscale)."""
    g = layer.gas
    kz = _kz(omega, g.c, sx)
    phi = kz * layer.d
    B = np.abs(phi.imag)
    ep = np.exp(1j * phi - B)    # e^{+j phi} e^{-B}, |.| <= 1
    em = np.exp(-1j * phi - B)   # e^{-j phi} e^{-B}, |.| <= 1
    cos_s = 0.5 * (ep + em)
    sin_s = (ep - em) / 2j
    # sin(phi)/kz, robust as kz -> 0 (limit d)
    small = np.abs(phi) < 1e-8
    kz_safe = np.where(small, 1.0, kz)
    sinc_s = np.where(small, layer.d * np.exp(-B), sin_s / kz_safe)
    rw = g.rho * omega
    a = cos_s
    b = 1j * rw * sinc_s                             # j Z sin(phi), Z = rho omega / kz
    rw_safe = np.where(rw == 0, 1.0, rw)
    c = np.where(rw == 0, 0.0, 1j * kz * sin_s / rw_safe)  # j sin(phi) / Z
    d = cos_s
    if layer.tie_K:
        return _tie_layer_scaled(a, b, c, d, B, layer.tie_K, omega)
    return a, b, c, d, B


def _tie_layer_scaled(a, b, c, d, B, K, omega):
    """Add a smeared tie spring K [Pa/m] in parallel with a gas layer.

    Derivation (true, unscaled quantities; state (p, v), v positive downward,
    1 = top face, 2 = bottom face of the layer). Impedance (2-port) form of the gas
    layer, from p1 = a p2 + b v2, v1 = c p2 + d v2 (ad - bc = 1):

        p1 = (a/c) v1 - (1/c) v2,   p2 = (1/c) v1 - (d/c) v2.

    The tie carries tension K e, e = w2 - w1 = (v2 - v1)/(j w) (thickness change).
    It pulls the top face down and the bottom face up by K e, so in the face
    equations it acts exactly like replacing the gas pressures by
    P_i = p_i - K e = p_i + s (v1 - v2), s = K/(j w):

        [P1, P2] = (Z_gas + s [[1, -1], [1, -1]]) [v1, v2]

    (= the textbook (K/jw)[[1,-1],[-1,1]] once v2 is taken positive *into* the layer
    from below, i.e. v2 -> -v2 and P2 -> -P2). Converting back to chain form gives

        A' = (a + s c)/(1 + s c),            C' = c/(1 + s c),
        D' = (d + s c)/(1 + s c),            B' = (b + s(a + d - 2))/(1 + s c),

    (B' simplified with ad - bc = 1; A'D' - B'C' = 1 still, so reciprocity and the
    lossless energy balance survive). We use this direct chain form rather than
    building Z explicitly because Z has 1/c, singular wherever c = 0.

    Scaling. The inputs are scaled: true = e^B * (a_s, b_s, c_s, d_s). With
    u = e^{-B} (may underflow to 0 for huge kappa d):

        M' = N / (u + s c_s),  N = [[a_s + s c_s, b_s + s(a_s + d_s - 2u)],
                                    [c_s,         d_s + s c_s]],

    every entry of N is O(1 + |s|) and the denominator stays finite, so this is
    overflow-safe. Physically: once e^{kappa d} is large any finite tie dominates the
    face-to-face coupling (M' ~ N/(s c_s), no exponential factor), i.e. the ties
    are the transmission path, as expected. We return N/n_max * |den|/den with log
    scale ln(n_max/|den|) so the caller's (scaled entries, log scale) bookkeeping
    is unchanged.

    Note 1 + s c -> 1 - K d (1/V^2 - 1/c_g^2)/rho_g as w -> 0 (evanescent layer:
    the gas has *negative* static stiffness rho_g/(a^2 d) under the travelling
    load); at K = rho_g/(a^2 d) the two cancel, C' blows up and T -> 0 at low f.
    """
    s = K / (1j * omega)
    u = np.exp(-B)
    sc = s * c
    den = u + sc
    den = np.where(den == 0, 1e-300, den)
    n11 = a + sc
    n12 = b + s * (a + d - 2.0 * u)
    n21 = c
    n22 = d + sc
    nmax = np.maximum.reduce([np.abs(n11), np.abs(n12), np.abs(n21), np.abs(n22)])
    nmax = np.where(nmax == 0, 1.0, nmax)
    ph = np.abs(den) / den / nmax
    logscale = np.log(nmax) - np.log(np.abs(den))
    return n11 * ph, n12 * ph, n21 * ph, n22 * ph, logscale


def _stack_matrix(stack: Stack, omega, sx):
    """Scaled total matrix (A, B, C, D) and log scale L (true M = exp(L) * scaled)."""
    one = np.ones_like(omega, dtype=complex)
    A, B, C, D = one.copy(), 0 * one, 0 * one, one.copy()
    L = np.zeros_like(omega, dtype=float)
    for layer in stack.layers:
        if isinstance(layer, GasLayer):
            a, b, c, d, s = _layer_matrix_scaled(layer, omega, sx)
            L = L + s
        elif hasattr(layer, "impedance"):
            a, b, c, d = one, layer.impedance(omega, sx) * one, 0 * one, one
        else:
            raise TypeError(f"unsupported layer {layer!r}")
        A, B, C, D = A * a + B * c, A * b + B * d, C * a + D * c, C * b + D * d
    return A, B, C, D, L


def transfer(stack: Stack, f, V, atm, reference: str = "interface"):
    """Complex pressure transmission T and reflection R at frequencies f [Hz].

    T = transmitted pressure amplitude at the bottom face / incident amplitude at the
    top face (same x). R = reflected / incident at the top face. Bounding media are
    ambient air (atm.rho, atm.c) so power transmission is |T|^2.
    ``reference="air"`` multiplies T by exp(+j Re(k_z,air) D) (D = total gas thickness),
    i.e. the insertion transfer relative to the same thickness of ambient air
    (a pure time shift; an all-ambient-air stack then gives T = 1).
    f = 0 is evaluated as the f -> 0 limit (T -> 1 for any stack without screens).
    """
    f = np.atleast_1d(np.asarray(f, dtype=float))
    f_eval = np.where(f == 0, 1e-9, np.abs(f))
    omega = 2 * np.pi * f_eval
    sx = _slowness(V)
    A, B, C, D, L = _stack_matrix(stack, omega, sx)
    kza = _kz(omega, atm.c, sx)
    Za = atm.rho * omega / kza
    den = A + B / Za + Za * C + D
    T = 2.0 * np.exp(-L) / den
    R = 2.0 * (A + B / Za) / den - 1.0
    if reference == "air":
        T = T * np.exp(1j * kza.real * stack.thickness)
    elif reference != "interface":
        raise ValueError("reference must be 'interface' or 'air'")
    return T, R


def apply(stack: Stack, t, p, V, atm, pad_factor: int = 4, reference: str = "air"):
    """Transmitted waveform (same length/time base as ``p``) in the air below.

    Uses rfft with the record zero-padded to >= pad_factor * len(p) and centred in the
    padded buffer, so both causal tails and the (physical, see module doc) non-causal
    precursor of evanescent layers have room before being cropped back. Inputs should
    already carry some quiet margin on both sides; evanescent layers smooth with a
    heavy-tailed (~1/t^2, Lorentzian) kernel of width ~ d sqrt(1/V^2 - 1/c_gas^2), whose
    far tails are cropped with the record. By default the transfer is
    referenced to an equal thickness of ambient air (pure delay removed).
    """
    t = np.asarray(t, dtype=float)
    p = np.asarray(p, dtype=float)
    n = p.size
    dt = t[1] - t[0]
    nfft = 1 << int(np.ceil(np.log2(max(2, pad_factor * n))))
    off = (nfft - n) // 2
    buf = np.zeros(nfft)
    buf[off:off + n] = p
    P = np.fft.rfft(buf)
    f = np.fft.rfftfreq(nfft, dt)
    T, _ = transfer(stack, f, V, atm, reference=reference)
    out = np.fft.irfft(P * T, nfft)
    return out[off:off + n]


def lift_margin(stack: Stack, atm) -> float:
    """Net static lift per unit area [kg/m^2] of the stack treated as uniform slabs:

        sigma_net = sum_gas (rho_air - rho_gas) d  -  sum_elements areal_mass.

    Positive = floats with that much spare per m^2. This slab model is the
    optimistic upper bound: real envelopes are lens/ellipse shaped, so the mean gas
    depth of a lens of peak depth d is ~ (pi/4) d (elliptic section) or less, while
    the membrane area (hence mass) per m^2 of footprint is >= 1.
    """
    sigma = 0.0
    for layer in stack.layers:
        if isinstance(layer, GasLayer):
            sigma += (atm.rho - layer.gas.rho) * layer.d
        else:
            sigma -= layer.areal_mass
    return float(sigma)
