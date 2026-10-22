"""Secondary sonic booms: upper atmosphere to 200 km, 3-D ray tracing of returning rays,
nonlinear + absorbing propagation along them, ground binning.

Public API
----------
UpperAtmosphere(z_max=200e3, dz=25.0)        tables of T, p, rho, c, M, gamma, x_O, mu, h
    .at(z) -> dict, .T(z), .p(z), .rho(z), .c(z) ... (vectorised, linear interpolation)
upper_atm()                                  cached default instance
USSA76_TABLE                                 reference values above 86 km (validation)
absorption_exponent(f, ua, z, ds)            ln H(f) over air-path ds at altitude z
upper_wind(season, toward_deg)               representative 0-150 km wind profiles
Medium(wind=None, smooth_m=1000.0)           smoothed c(z), wind(z) and their gradients for rays
launch_slowness(V, z_a, phi, medium)         ground-frame slowness of a Mach-cone ray
trace(phis, V, z_a, medium, ...)             vectorised RK4 ray tracer with ground/shield mirrors
fan(ac, medium, phis, shield=None)           landing table of a ray fan (secondary hits classified)
ray_path(ac, phi, medium, shield=None)       recorded central path + ray-tube Q (3 rays)
march_family(ac, phi, medium, family, ...)   near field -> nonlinear march -> descent -> h(t)
bin_ground(...)                              caustic-regularised ground levels vs lateral range

Model summary (see results/secondary_booms.md for the discussion)
-----------------------------------------------------------------
Atmosphere.  0-86 km: boomsim.atmosphere.std_atm (US Standard Atmosphere 1976, USSA76).
Above 86 km USSA76 changes form (NOAA-S/T 76-1562, 1976, Part 1 sec. 1.2-1.3): kinetic
temperature in four segments,
    86-91 km   T = 186.8673 K (isothermal)
    91-110 km  T = Tc + A sqrt(1 - ((z - 91)/a)^2),  Tc = 263.1905 K, A = -76.3232 K,
               a = -19.9429 km (elliptical)
    110-120 km T = 240 + 12 (z - 110)   [K, z in km]
    > 120 km   T = T_inf - (T_inf - 360) exp(-lambda xi),  T_inf = 1000 K,
               lambda = 0.01875 /km, xi = (z - 120)(r0 + 120)/(r0 + z), r0 = 6356.766 km
and diffusive separation of N2, O, O2, Ar, He above ~100 km, which lowers the mean molar
mass M.  Rather than integrating the species diffusion equations we take M(z) from the
USSA76 tables (PCHIP through the tabulated values at 86, 90, 100, 110, 120, 130, 150,
200 km, ``USSA76_TABLE``) and integrate hydrostatic balance dp/dz = -rho g(z) with
rho = p M/(R* T), g = g0 (r0/(r0 + z))^2.  This reproduces the tabulated p and rho to
<~2 % (tests).  Atomic-oxygen mole fraction x_O = (M_air - M)/(M_air - M_O) (N2/O2/O
mixture approximation; 0.37 at 150 km vs ~0.35 in USSA76), gamma from additive molar
heat capacities (c_v = R (2.5 (1 - x_O) + 1.5 x_O)), c = sqrt(gamma R* T/M).
Water vapour: ISO-profile humidity (boomsim.atmosphere.molar_humidity) to 30 km,
then a constant 5 ppmv mixing ratio (h = 5e-4 %), typical of the middle atmosphere.

Absorption.  Classical + rotational (bulk-viscosity) absorption
    alpha_cl = (2 pi^2 f^2 / (gamma p c)) mu [4/3 + mu_B/mu + (gamma - 1)/Pr],
    mu_B = 0.6 mu (1 - x_O) (rotational relaxation of N2/O2; atomic O has none),
    Pr = 0.72, mu = Sutherland's law (1.458e-6 T^1.5/(T + 110.4)),
in the form of Sutherland & Bass, "Atmospheric absorption in the atmosphere up to
160 km", J. Acoust. Soc. Am. 115 (2004) 1012-1032.  At 20 C, 1 atm this equals the ISO
9613-1 classical+rotational coefficient (1.84e-11 f^2/p_r) to 1 %.  It scales as 1/p:
at 100 km (p = 0.032 Pa) alpha(1 Hz) ~ 0.4 dB/km, at 130 km ~ 15 dB/km.  O2/N2 vibrational
relaxation (ISO 9613-1 with causal dispersion, boomsim.propagate.propagation_exponent)
is kept; it is negligible above ~60 km.  Continuum acoustics needs mean free path <<
wavelength: ~3 m at 120 km and ~30 m at 150 km (USSA76) against 300-1000 m for the
0.3-1 Hz waves that survive there; rays turning above ~140 km carry no energy anyway.

Rays.  Horizontally stratified moving medium w = (u(z), v(z)), Hamiltonian ray equations
in ground-frame slowness s (dispersion c|s| + w.s = 1):  s_x, s_y conserved,
    dx/dt = c^2 s_x/Omega + u,  dy/dt = c^2 s_y/Omega + v,  dz/dt = c^2 s_z/Omega,
    ds_z/dt = -(Omega c'/c + u' s_x + v' s_y),  Omega = 1 - w.s_h,
integrated in time (fixed-step RK4) so turning points are regular.  c and w are smoothed
with a 1 km Gaussian for the rays (piecewise-linear USSA76 gradients otherwise make
spurious micro-caustics).  Launch normals lie on the Mach cone about the air velocity
(aircraft crabbing into a cross wind as in boomsim.propagate).  The ground (z = 0) and,
optionally, the shield top are mirrors.

Amplitude.  Blokhintzev wave-action invariant along the tube (Pierce 1989, sec. 8-6):
    p^2 |Q| / (rho c^2 Omega) = const,   Q = (e_x x dr/dphi) . v_ray,
the (signed) flux area of the tube spanned by the steady-flight translation e_x and the
launch-azimuth derivative dr/dphi (computed from two neighbour rays phi +- dphi at equal
travel time; the triple product is unchanged by sliding along the rays, so Q is finite at
turning points).  At the ground Q = (dY/dphi) v_z.  Mirror reflections flip the sign of Q
and are removed; a genuine sign change is a caustic (-pi/2 phase shift = Hilbert
transform, applied to the waveform).

Waveform.  Augmented Burgers along the ray (boomsim.propagate.burgers_step, MUSCL/
Engquist-Osher, CFL-limited; spreading by the invariant; absorption/dispersion in the
frequency domain per <= 200 m chunk of air path, at the intrinsic frequency Omega f).
Delta-p/p_ambient grows by ~1e3 between 20 and 100 km, so the march is strongly
nonlinear on the way up; shocks are captured conservatively (equal-area rule), so the
N-wave decay and lengthening do not depend on resolving the shock thickness, and the
record is decimated as the waveform lengthens (kept at >= ``n_span`` samples per active
span).  ``max_dp_over_p`` is reported (weak-shock theory needs it << 1).
After the apex the march continues nonlinearly to a reference point on the descending
leg (``dz_ref`` below the apex, or earlier if the tube starts to collapse); below it
Delta-p/p falls by orders of magnitude and the rest is linear: the waveform is stored as
h(t) = p sqrt(|Q|/(rho c^2 Omega)) with the descent absorption applied.

Ground levels and caustics.  Ray amplitude at the ground is p = h sqrt(rho c^2 Omega/
|Y_phi v_z|), singular at the fold caustic that bounds every secondary carpet.  Following
energy conservation we bin the landing energy in lateral range Y:
    p_bin^2 ~ (rho c^2 Omega / |v_z|)_g  sum_{rays in bin} E_h(phi) dphi / dY,
which is finite at caustics (dY/dphi cancels).  dY = 3 km by default ~ the Airy/caustic
layer width (lambda^2 L)^(1/3) for 0.3-2 Hz waves and L ~ 100 km (sensitivity 1.5-6 km).
Bins are formed per branch (a contiguous phi-interval with one turning region); the
default peak is the largest single-branch value, the "coherent" bound adds branch peaks.
The ground observer sees the incident wave times 1.9 (same convention as the primary).
"""
from __future__ import annotations

import functools
import math
from dataclasses import dataclass, field

import numpy as np
from scipy import fft as sfft
from scipy.interpolate import PchipInterpolator
from scipy.ndimage import gaussian_filter1d

from .atmosphere import P_REF, molar_humidity, std_atm
from .propagate import WindProfile, burgers_step, propagation_exponent
from .source import near_field

R_STAR = 8.31432           # J/(mol K), USSA76
G0 = 9.80665
R0_EARTH = 6356766.0       # m (USSA76 r0)
M_AIR = 28.9644e-3
M_O = 15.9994e-3
BETA_AIR = 1.2
PR_NUMBER = 0.72
BULK_RATIO = 0.6           # mu_B/mu for N2/O2 (rotational relaxation)
H2O_MIDDLE = 5e-4          # % molar (5 ppmv) above 30 km
REFLECT_UP = 0.9           # ground reflection of the incident wave (1.9 = 1 + 0.9)
GROUND_FACTOR = 1.9
PSF = 47.880258

# USSA76 tabulated values: z [km] -> (T [K], p [Pa], rho [kg/m^3], M [g/mol]).
# NOAA-S/T 76-1562 (1976), Table I/II, geometric altitude.
USSA76_TABLE = {
    86.0: (186.87, 3.7338e-1, 6.958e-6, 28.95),
    90.0: (186.87, 1.8359e-1, 3.416e-6, 28.91),
    100.0: (195.08, 3.2011e-2, 5.604e-7, 28.40),
    110.0: (240.00, 7.1042e-3, 9.708e-8, 27.27),
    120.0: (360.00, 2.5382e-3, 2.222e-8, 26.20),
    130.0: (469.27, 1.2505e-3, 8.152e-9, 25.44),
    150.0: (634.39, 4.5422e-4, 2.076e-9, 24.10),
    200.0: (854.56, 8.4736e-5, 2.541e-10, 21.30),
}


# ==========================================================================================
# atmosphere
# ==========================================================================================
def ussa76_temperature_above86(z):
    """USSA76 kinetic temperature [K] for 86 km <= z [m] (formulas in module docstring)."""
    zk = np.asarray(z, dtype=float) / 1e3
    T = np.full(zk.shape, 186.8673)
    m2 = (zk >= 91.0) & (zk < 110.0)
    T = np.where(m2, 263.1905 - 76.3232 * np.sqrt(np.clip(1.0 - ((zk - 91.0) / -19.9429) ** 2,
                                                          0.0, None)), T)
    m3 = (zk >= 110.0) & (zk < 120.0)
    T = np.where(m3, 240.0 + 12.0 * (zk - 110.0), T)
    m4 = zk >= 120.0
    r0 = R0_EARTH / 1e3
    xi = (zk - 120.0) * (r0 + 120.0) / (r0 + zk)
    T = np.where(m4, 1000.0 - 640.0 * np.exp(-0.01875 * xi), T)
    return T


class UpperAtmosphere:
    """Tabulated atmosphere 0..z_max (USSA76 below 86 km, extended model above)."""

    def __init__(self, z_max=200e3, dz=25.0):
        self.dz = float(dz)
        z = np.arange(0.0, z_max + 0.5 * dz, dz)
        self.z = z
        lo = z <= 86000.0
        a = std_atm(z[lo])
        T = np.empty_like(z)
        p = np.empty_like(z)
        M = np.full_like(z, M_AIR)
        T[lo] = a.T
        p[lo] = a.p
        hi = ~lo
        zt = np.array(sorted(USSA76_TABLE)) * 1e3
        Mt = np.array([USSA76_TABLE[k][3] for k in sorted(USSA76_TABLE)]) * 1e-3
        M[hi] = PchipInterpolator(zt, Mt)(z[hi])
        T[hi] = ussa76_temperature_above86(z[hi])
        # hydrostatic integration from the 86 km value (trapezoid in ln p)
        i86 = int(np.flatnonzero(lo)[-1])
        g = G0 * (R0_EARTH / (R0_EARTH + z)) ** 2
        k = M * g / (R_STAR * T)
        lnp = np.log(p[i86]) - np.concatenate(
            [[0.0], np.cumsum(0.5 * (k[i86:-1] + k[i86 + 1:]) * dz)])
        # tie to the tabulated 86 km pressure (std_atm gives 0.3734 Pa there)
        p[i86:] = np.exp(lnp)
        self.T, self.p, self.M = T, p, M
        self.rho = p * M / (R_STAR * T)
        xO = np.clip((M_AIR - M) / (M_AIR - M_O), 0.0, 1.0)
        self.xO = xO
        cv = 2.5 * (1.0 - xO) + 1.5 * xO
        self.gamma = (cv + 1.0) / cv
        self.c = np.sqrt(self.gamma * R_STAR * T / M)
        self.mu = 1.458e-6 * T ** 1.5 / (T + 110.4)
        h = np.full_like(z, H2O_MIDDLE)
        m30 = z <= 30000.0
        h[m30] = molar_humidity(z[m30])
        # blend smoothly into the constant mixing ratio
        h[m30] = np.maximum(h[m30], 0.0)
        self.h = h

    # fast linear interpolation on the uniform grid
    def _w(self, z):
        x = np.clip(np.asarray(z, dtype=float) / self.dz, 0.0, self.z.size - 1.000001)
        i = x.astype(int)
        return i, x - i

    def interp(self, name, z):
        arr = getattr(self, name)
        i, w = self._w(z)
        return arr[i] * (1.0 - w) + arr[i + 1] * w

    def interp_log(self, name, z):
        arr = np.log(getattr(self, name))
        i, w = self._w(z)
        return np.exp(arr[i] * (1.0 - w) + arr[i + 1] * w)

    def state(self, z):
        """dict of T, p, rho, c, gamma, mu, xO, h at z (p and rho log-interpolated)."""
        return dict(T=self.interp("T", z), p=self.interp_log("p", z),
                    rho=self.interp_log("rho", z), c=self.interp("c", z),
                    gamma=self.interp("gamma", z), mu=self.interp("mu", z),
                    xO=self.interp("xO", z), h=self.interp("h", z))


@functools.lru_cache(maxsize=4)
def upper_atm(z_max=200e3, dz=25.0) -> UpperAtmosphere:
    return UpperAtmosphere(z_max, dz)


def classical_coefficient(st):
    """a such that alpha_cl+rot = a f^2 [Np s^2/m] (Sutherland & Bass 2004 form)."""
    mu = st["mu"]
    g = st["gamma"]
    brack = 4.0 / 3.0 + BULK_RATIO * (1.0 - st["xO"]) + (g - 1.0) / PR_NUMBER
    return 2.0 * np.pi ** 2 * mu * brack / (g * st["p"] * st["c"])


def absorption_exponent(f, st, ds, relaxation=True, dispersion=True):
    """ln H(f) over air path ds [m] at the state ``st`` (dict from UpperAtmosphere.state)."""
    f = np.asarray(f, dtype=float)
    lnH = -classical_coefficient(st) * f * f * ds + 0j
    if relaxation:
        lnH = lnH + propagation_exponent(f, st["T"], st["p"], st["c"], st["h"], ds,
                                         thermoviscous=False, relaxation=True,
                                         dispersion=dispersion)
    return lnH


# ==========================================================================================
# winds
# ==========================================================================================
# Representative mid-latitude (~40 N) zonal-mean zonal wind [m/s], westerly positive,
# vs altitude [km].  Troposphere/lower stratosphere as boomsim.propagate.jet_stream_winter
# (winter) / a weaker summer jet; middle atmosphere from the CIRA-86 climatology (Fleming,
# Chandra, Barnett & Corney, Adv. Space Res. 10(12) (1990) 11-59; Andrews, Holton & Leovy,
# Middle Atmosphere Dynamics, 1987, Fig. 1.3): winter westerly jet ~50-65 m/s near 60-65 km
# at 40 N (stronger, 70-80 m/s, poleward), summer easterly jet ~-50 to -60 m/s near 65-70 km;
# lower thermosphere (90-130 km) dominated by tides of +-20-50 m/s (set to zero here and
# treated as a sensitivity).  Values are rounded from the published cross-sections (not a
# digitisation) -- representative, not a climatology.
UPPER_KM = [0, 1, 3, 6, 9, 10.5, 12, 15, 17, 20, 22, 26, 30, 35, 40, 45, 50, 55, 60, 65,
            70, 75, 80, 85, 90, 95, 100, 150]
WINTER_MS = [3, 10, 18, 35, 55, 65, 65, 35, 22, 15, 15, 15, 18, 25, 33, 42, 52, 60, 64,
             62, 55, 45, 32, 20, 10, 5, 0, 0]
SUMMER_MS = [2, 5, 8, 15, 22, 25, 25, 12, 5, 0, -5, -10, -15, -22, -30, -38, -46, -53,
             -58, -60, -55, -45, -30, -15, -5, 0, 0, 0]


def upper_wind(season="winter", toward_deg=0.0, tide=0.0):
    """WindProfile in track axes for the zonal-mean profile ``season`` ('winter' westerly,
    'summer' easterly) blowing toward ``toward_deg`` relative to the flight direction for
    positive speeds (0 = tailwind, 180 = headwind, 90 = toward +y).  ``tide`` [m/s] adds a
    uniform lower-thermosphere wind between 95 and 140 km (same direction), a crude stand-in
    for the tidal winds that dominate there."""
    sp = np.array(WINTER_MS if season == "winter" else SUMMER_MS, dtype=float)
    zk = np.array(UPPER_KM, dtype=float)
    if tide:
        zk2 = np.arange(0.0, 151.0, 1.0)
        sp2 = np.interp(zk2, zk, sp)
        bump = tide * np.clip(np.minimum((zk2 - 95) / 10, (140 - zk2) / 10), 0, 1)
        zk, sp = zk2, sp2 + bump
    return WindProfile.from_speed(zk * 1e3, sp, toward_deg,
                                  name=f"{season}{'+tide' if tide else ''} -> {toward_deg:g}")


# ==========================================================================================
# ray medium and tracer
# ==========================================================================================
class Medium:
    """Ray medium: smoothed c(z) (from UpperAtmosphere) and wind (u, v)(z) + gradients."""

    def __init__(self, wind=None, smooth_m=1000.0, ua=None):
        self.ua = upper_atm() if ua is None else ua
        z = self.ua.z
        self.dz = self.ua.dz
        sig = smooth_m / self.dz
        self.c = gaussian_filter1d(self.ua.c, sig, mode="nearest")
        self.dc = np.gradient(self.c, self.dz)
        self.wind = wind
        if wind is None:
            self.u = np.zeros_like(z)
            self.v = np.zeros_like(z)
        else:
            u, v = wind(z)
            self.u = gaussian_filter1d(np.asarray(u, float), sig, mode="nearest")
            self.v = gaussian_filter1d(np.asarray(v, float), sig, mode="nearest")
        self.du = np.gradient(self.u, self.dz)
        self.dv = np.gradient(self.v, self.dz)
        self.calm = wind is None or (not np.any(self.u) and not np.any(self.v))

    def fields(self, z):
        x = np.clip(z / self.dz, 0.0, self.c.size - 1.000001)
        i = x.astype(int)
        w = x - i
        w1 = 1.0 - w

        def g(a):
            return a[i] * w1 + a[i + 1] * w

        return g(self.c), g(self.dc), g(self.u), g(self.du), g(self.v), g(self.dv)


def launch_slowness(V, z_a, phi, medium):
    """Ground-frame slowness (s_x, s_y, s_z) of Mach-cone rays at azimuth(s) ``phi`` [rad]
    (0 = straight down, +pi/2 = toward +y horizontal, pi = straight up) for an aircraft at
    air speed V [m/s] holding its ground track along +x (crab into a cross wind)."""
    phi = np.atleast_1d(np.asarray(phi, dtype=float))
    c_a, _, u_a, _, v_a, _ = (float(q[0]) for q in medium.fields(np.array([z_a])))
    # use the unsmoothed sound speed at flight level for the Mach angle (V = M c(z_a))
    sin_mu = min(1.0, c_a / V)
    cos_mu = math.sqrt(1.0 - sin_mu ** 2)
    eps = -math.asin(max(-1.0, min(1.0, v_a / V)))
    ea = np.array([math.cos(eps), math.sin(eps), 0.0])
    ep = np.array([-math.sin(eps), math.cos(eps), 0.0])
    n = (sin_mu * ea[:, None] + cos_mu * (np.sin(phi)[None, :] * ep[:, None]
                                          - np.cos(phi)[None, :] * np.array([0, 0, 1.0])[:, None]))
    den = c_a + u_a * n[0] + v_a * n[1]
    s = n / den
    return s[0], s[1], s[2]


def _rhs(z, sx, sy, sz, med):
    c, dc, u, du, v, dv = med.fields(z)
    om = 1.0 - u * sx - v * sy
    c2 = c * c / om
    return (c2 * sx + u, c2 * sy + v, c2 * sz,
            -(om * dc / c + du * sx + dv * sy), c)


@dataclass
class Shield:
    """Horizontal mirror at z (top face) for rays with |phi| < phi_max (rays launched from
    the aircraft that hit the strip).  mode: 'reflect' (mirror) or 'transmit' (an event is
    recorded but the ray continues)."""
    z: float = 20000.0
    phi_max: float = math.atan2(500.0, 250.0)
    mode: str = "reflect"


def trace(phis, V, z_a, medium, dt=0.25, t_max=2400.0, z_top=160e3, shield=None,
          max_ground=2, record=False, z0=None):
    """Vectorised RK4 tracer.  Returns dict with per-ray hit tables (ground/shield events)
    and, if ``record``, full time histories (t, x, y, z, sz, sigma) [n_steps, n_rays].

    Each ray stops after ``max_ground`` ground hits, above ``z_top`` (status 'escaped') or
    at t_max ('timeout').  Hits: kind (0 ground, 1 shield), t, x, y, sz, vz, apex of the leg
    that ends at the hit (max z since the previous event or launch), sigma.
    """
    phis = np.atleast_1d(np.asarray(phis, dtype=float))
    n = phis.size
    sx, sy, sz = launch_slowness(V, z_a, phis, medium)
    x = np.zeros(n)
    y = np.zeros(n)
    z = np.full(n, float(z_a if z0 is None else z0))
    sig = np.zeros(n)
    alive = np.ones(n, bool)
    status = np.zeros(n, int)          # 0 alive/done-ground, 1 escaped, 2 timeout
    nground = np.zeros(n, int)
    shield_ok = np.zeros(n, bool) if shield is None else (np.abs(phis) < shield.phi_max)
    shield_done = np.zeros(n, bool)
    apex = z.copy()
    H = 2 + max_ground
    hk = np.full((n, H), -1, int)
    ht = np.full((n, H), np.nan)
    hx = np.full((n, H), np.nan)
    hy = np.full((n, H), np.nan)
    hvz = np.full((n, H), np.nan)
    hsz = np.full((n, H), np.nan)
    hap = np.full((n, H), np.nan)
    hsig = np.full((n, H), np.nan)
    hstep = np.full((n, H), -1, int)
    nh = np.zeros(n, int)
    rec = [] if record else None
    t = 0.0
    step = 0
    if record:
        rec.append(np.stack([x, y, z, sz, sig]))

    def add_hit(idx, kind, frac, xs, ys, zs_, szs, sigs, vz):
        for j, i in enumerate(idx):
            k = nh[i]
            if k >= H:
                continue
            hk[i, k] = kind
            ht[i, k] = t + frac[j] * dt
            hx[i, k] = xs[j]
            hy[i, k] = ys[j]
            hsz[i, k] = szs[j]
            hvz[i, k] = vz[j]
            hap[i, k] = apex[i]
            hsig[i, k] = sigs[j]
            hstep[i, k] = step
            nh[i] = k + 1

    while np.any(alive) and t < t_max:
        a = np.flatnonzero(alive)
        X, Y, Z, SZ, SG = x[a], y[a], z[a], sz[a], sig[a]
        SX, SY = sx[a], sy[a]
        k1 = _rhs(Z, SX, SY, SZ, medium)
        k2 = _rhs(Z + 0.5 * dt * k1[2], SX, SY, SZ + 0.5 * dt * k1[3], medium)
        k3 = _rhs(Z + 0.5 * dt * k2[2], SX, SY, SZ + 0.5 * dt * k2[3], medium)
        k4 = _rhs(Z + dt * k3[2], SX, SY, SZ + dt * k3[3], medium)
        inc = [dt / 6.0 * (k1[j] + 2 * k2[j] + 2 * k3[j] + k4[j]) for j in range(5)]
        Xn, Yn, Zn, SZn, SGn = X + inc[0], Y + inc[1], Z + inc[2], SZ + inc[3], SG + inc[4]
        # shield crossing (downward, first time)
        if shield is not None:
            m = shield_ok[a] & ~shield_done[a] & (Z >= shield.z) & (Zn < shield.z)
            if np.any(m):
                idx = a[m]
                fr = (Z[m] - shield.z) / (Z[m] - Zn[m])
                xs = X[m] + fr * (Xn[m] - X[m])
                ys = Y[m] + fr * (Yn[m] - Y[m])
                szs = SZ[m] + fr * (SZn[m] - SZ[m])
                sgs = SG[m] + fr * (SGn[m] - SG[m])
                c, _, u, _, v, _ = medium.fields(np.full(idx.size, shield.z))
                om = 1 - u * sx[idx] - v * sy[idx]
                vz = c * c * szs / om
                add_hit(idx, 1, fr, xs, ys, None, szs, sgs, vz)
                shield_done[idx] = True
                apex[idx] = shield.z
                if shield.mode == "reflect":
                    Zn[m] = 2 * shield.z - Zn[m]
                    SZn[m] = -SZn[m]
        # ground crossing
        m = (Zn < 0.0) & (Z >= 0.0)
        if np.any(m):
            idx = a[m]
            fr = Z[m] / (Z[m] - Zn[m])
            xs = X[m] + fr * (Xn[m] - X[m])
            ys = Y[m] + fr * (Yn[m] - Y[m])
            szs = SZ[m] + fr * (SZn[m] - SZ[m])
            sgs = SG[m] + fr * (SGn[m] - SG[m])
            c, _, u, _, v, _ = medium.fields(np.zeros(idx.size))
            om = 1 - u * sx[idx] - v * sy[idx]
            vz = c * c * szs / om
            add_hit(idx, 0, fr, xs, ys, None, szs, sgs, vz)
            apex[idx] = 0.0
            nground[idx] += 1
            Zn[m] = -Zn[m]
            SZn[m] = -SZn[m]
            stop = nground[idx] >= max_ground
            alive[idx[stop]] = False
        x[a], y[a], z[a], sz[a], sig[a] = Xn, Yn, Zn, SZn, SGn
        apex[a] = np.maximum(apex[a], Zn)
        esc = alive & (z > z_top)
        status[esc] = 1
        alive[esc] = False
        t += dt
        step += 1
        if record:
            rec.append(np.stack([x, y, z, sz, sig]))
    status[alive] = 2
    out = dict(phi=phis, sx=sx, sy=sy, status=status, n_hits=nh, kind=hk, t=ht, x=hx, y=hy,
               vz=hvz, sz=hsz, apex=hap, sigma=hsig, step=hstep, dt=dt,
               shield_mode=None if shield is None else shield.mode)
    if record:
        R = np.array(rec)              # [steps, 5, n]
        out["rec"] = R
        out["t_rec"] = np.arange(R.shape[0]) * dt
    return out


# ==========================================================================================
# fan: landing tables and family classification
# ==========================================================================================
SECONDARY_APEX = 25000.0     # a landing whose last leg peaked above this is "secondary"


def classify_hits(tr, z_a):
    """For each ray, list of (hit index, family) for secondary landings.  family:
    'direct' (no earlier ground/shield event), 'indirect' (after a ground reflection),
    'shield' (after a shield reflection), 'transmit' (passed a transmit-mode shield event;
    also ground-reflected or not)."""
    out = []
    for i in range(tr["phi"].size):
        fams = []
        prior_ground = prior_shield = prior_secondary = False
        for k in range(tr["n_hits"][i]):
            kind = tr["kind"][i, k]
            if kind == 1:
                prior_shield = prior_shield or tr.get("shield_mode") == "reflect"
                continue
            if tr["apex"][i, k] > SECONDARY_APEX:
                if prior_secondary:
                    fam = "multihop"          # second bounce of a secondary: not scored
                elif prior_shield:
                    fam = "shield"
                elif prior_ground:
                    fam = "indirect"
                else:
                    fam = "direct"
                fams.append((k, fam))
                prior_secondary = True
            prior_ground = True
        out.append(fams)
    return out


def fan(ac, medium, phis, shield=None, dt=0.25, **kw):
    """Trace a fan and tabulate secondary landings.  Returns dict with the trace and a
    list of records (phi, family, Y, X, t, vz, apex, hit index)."""
    tr = trace(phis, ac.V, ac.altitude, medium, dt=dt, shield=shield, **kw)
    cls = classify_hits(tr, ac.altitude)
    rows = []
    for i, fams in enumerate(cls):
        for k, fam in fams:
            rows.append(dict(i=i, phi=tr["phi"][i], family=fam, k=k, Y=tr["y"][i, k],
                             X=tr["x"][i, k], t=tr["t"][i, k], vz=tr["vz"][i, k],
                             apex=tr["apex"][i, k], sx=tr["sx"][i], sy=tr["sy"][i]))
    return tr, rows


def lateral_cutoff_phi(ac, medium, side=1, tol=1e-4):
    """Largest downward |phi| whose ray reaches the ground directly (bisection)."""
    lo, hi = 0.0, math.pi / 2
    while hi - lo > tol:
        mid = 0.5 * (lo + hi)
        tr = trace([side * mid], ac.V, ac.altitude, medium, dt=0.5, t_max=600.0,
                   max_ground=1, z_top=ac.altitude + 2000.0)
        ok = tr["n_hits"][0] > 0 and tr["apex"][0, 0] <= ac.altitude + 1.0
        lo, hi = (mid, hi) if ok else (lo, mid)
    return lo


# ==========================================================================================
# recorded path + ray tube
# ==========================================================================================
@dataclass
class Path:
    phi: float
    t: np.ndarray
    x: np.ndarray
    y: np.ndarray
    z: np.ndarray
    sigma: np.ndarray       # air-frame path length (c dt integrated)
    omega: np.ndarray
    Q: np.ndarray           # mirror-corrected signed tube flux area
    Q_ok: np.ndarray
    events: list            # (step index, kind) kind: 'ground', 'shield'
    i_end: int              # index of the family landing step (path ends here)
    landing: dict
    sx: float
    sy: float
    caustics: list = field(default_factory=list)


def _fields_full(medium, z, sx, sy):
    c, dc, u, du, v, dv = medium.fields(z)
    om = 1.0 - u * sx - v * sy
    return c, u, v, om


def ray_path(ac, phi, medium, family, shield=None, dt=0.2, dphi=2e-4, t_max=2400.0,
             z_top=160e3):
    """Record the central ray at ``phi`` up to its ``family`` landing (see classify_hits)
    plus the tube flux area Q from neighbours phi +- dphi.  Returns Path or None if the
    ray has no such landing."""
    phis = np.array([phi - dphi, phi, phi + dphi])
    tr = trace(phis, ac.V, ac.altitude, medium, dt=dt, t_max=t_max, z_top=z_top,
               shield=shield, max_ground=2, record=True)
    if family == "primary":
        hit = [k for k in range(tr["n_hits"][1])
               if tr["kind"][1, k] == 0 and tr["apex"][1, k] <= SECONDARY_APEX][:1]
    else:
        cls = classify_hits(tr, ac.altitude)[1]
        hit = [k for k, fam in cls if fam == family]
    if not hit:
        return None
    k = hit[0]
    i_end = int(tr["step"][1, k])
    R = tr["rec"]
    # central-ray events up to the landing
    events = []
    for j in range(tr["n_hits"][1]):
        if tr["step"][1, j] < i_end:
            events.append((int(tr["step"][1, j]), "ground" if tr["kind"][1, j] == 0 else "shield"))
    sx, sy = float(tr["sx"][1]), float(tr["sy"][1])
    x, y, z, sz, sg = (R[: i_end + 1, j, 1] for j in range(5))
    c, u, v, om = _fields_full(medium, z, sx, sy)
    vx = c * c * sx / om + u
    vy = c * c * sy / om + v
    vz = c * c * sz / om
    rp = (R[: i_end + 1, 0:3, 2] - R[: i_end + 1, 0:3, 0]) / (2 * dphi)
    # Q = (e_x x r_phi) . v = -r_phi_z v_y + r_phi_y v_z
    Q = -rp[:, 2] * vy + rp[:, 1] * vz
    # mirror sign correction (central ray events; reflect-mode shield flips too)
    sign = np.ones_like(Q)
    bad = np.zeros(Q.size, bool)
    all_steps = []
    for r in range(3):
        for j in range(tr["n_hits"][r]):
            all_steps.append(int(tr["step"][r, j]))
    for st, kind in events:
        if kind == "ground" or (shield is not None and shield.mode == "reflect"):
            sign[st:] *= -1.0
    for st in all_steps:
        bad[max(0, st - 3): st + 4] = True
    Qc = Q * sign
    ok = ~bad
    # repair the samples next to events (log-|Q| interpolation, keep sign of neighbours)
    if np.any(bad) and np.sum(ok) > 2:
        idx = np.arange(Q.size)
        mag = np.exp(np.interp(idx, idx[ok], np.log(np.abs(Qc[ok]) + 1e-300)))
        sgn = np.sign(np.interp(idx, idx[ok], Qc[ok]))
        sgn[sgn == 0] = 1.0
        Qc = np.where(bad, sgn * mag, Qc)
    caustics = [int(i) for i in np.flatnonzero(np.diff(np.sign(Qc)) != 0) if ok[i] and ok[i + 1]]
    landing = dict(t=tr["t"][1, k], x=tr["x"][1, k], y=tr["y"][1, k], vz=tr["vz"][1, k],
                   apex=tr["apex"][1, k])
    return Path(phi=phi, t=np.arange(i_end + 1) * dt, x=x, y=y, z=z, sigma=sg, omega=om, Q=Qc,
                Q_ok=ok, events=events, i_end=i_end, landing=landing, sx=sx, sy=sy,
                caustics=caustics)


# ==========================================================================================
# waveform march
# ==========================================================================================
def _extend(t, p, front, back):
    dt = t[1] - t[0]
    n = p.size + front + back
    nn = sfft.next_fast_len(n, real=True)
    back += nn - n
    p = np.concatenate([np.zeros(front), p, np.zeros(back)])
    t = t[0] + (np.arange(p.size) - front) * dt
    return t, p


def _span(p, rel=1e-4):
    a = np.abs(p)
    idx = np.flatnonzero(a > rel * a.max())
    if idx.size == 0:
        return 0, p.size
    return idx[0], idx[-1] + 1


def _manage(t, p, n_span, fs_min, edge=0.05):
    """Keep the record padded (extend if signal near the ends), crop excess zeros and
    decimate by 2 (spectral) while the active span holds > 2 n_span samples."""
    dt = t[1] - t[0]
    i0, i1 = _span(p)
    span = i1 - i0
    m = max(8, int(edge * p.size))
    if i0 < m or i1 > p.size - m:
        pad = max(span // 2, 64)
        t, p = _extend(t, p, pad, pad)
        i0, i1 = _span(p)
        span = i1 - i0
    # decimate
    while span > 2 * n_span and 1.0 / dt > 2 * fs_min:
        n = p.size - (p.size % 2)
        P = sfft.rfft(p[:n])
        half = n // 2
        P2 = P[: half // 2 + 1].copy()
        P2[-1] = P2[-1].real
        p = sfft.irfft(P2, n=half) * 0.5
        t = t[0] + np.arange(half) * 2 * dt
        dt *= 2
        i0, i1 = _span(p)
        span = i1 - i0
    # crop if the record is > 4x the span
    if p.size > 4 * span + 256:
        pad = max(span, 128)
        a, b = max(0, i0 - pad), min(p.size, i1 + pad)
        t, p = t[a:b], p[a:b]
        t, p = _extend(t, p, 0, 0)
    return t, p


def _filter(t, p, H_of_f, pad_frac=1.0):
    """Apply a frequency response (callable f -> H, rfft convention) with zero padding."""
    n = p.size
    pad = int(pad_frac * n)
    nn = sfft.next_fast_len(n + 2 * pad, real=True)
    buf = np.zeros(nn)
    buf[pad:pad + n] = p
    dt = t[1] - t[0]
    f = sfft.rfftfreq(nn, dt)
    out = sfft.irfft(sfft.rfft(buf) * H_of_f(f), n=nn)
    t2 = t[0] + (np.arange(nn) - pad) * dt
    return t2, out


def hilbert_shift(p):
    """-pi/2 phase shift (caustic passage), rfft convention: positive frequencies x (-j)."""
    P = sfft.rfft(p)
    P = P * (-1j)
    P[0] = 0.0
    return sfft.irfft(P, n=p.size)


@dataclass
class MarchResult:
    t: np.ndarray            # stored waveform time axis (ground-frame retarded time)
    h: np.ndarray            # normalised waveform at the ground (see module docstring)
    E_h: float
    max_dp_over_p: float
    z_ref: float
    apex: float
    n_caustic_up: int
    n_caustic_down: int
    p_ref_peak: float
    info: dict


def march(t, p, path: Path, i_start, ops=None, ua=None, linear=False, cfl=0.8,
          ds_lin=200.0, n_span=4096, fs_min=500.0, dz_ref=15000.0, q_floor=0.1,
          nonlinear_to_ground=False, beta=BETA_AIR, stats=None):
    """March (t, p) along ``path`` from sample i_start.  ``ops``: dict step index ->
    callable(t, p) -> (t, p) applied when the march reaches that sample (shield R/T,
    ground reflection).  Returns MarchResult (h normalised at the ground)."""
    ua = upper_atm() if ua is None else ua
    ops = ops or {}
    z = path.z
    st_all = ua.state(z)
    rho, c = st_all["rho"], st_all["c"]
    pa = st_all["p"]
    om = path.omega
    Q = path.Q
    inv = rho * c * c * om / np.maximum(np.abs(Q), 1e-300)     # p^2 ~ inv  (Blokhintzev)
    # aerial caustics (thermospheric returns fold near the apex): |Q| floored at q_floor x
    # its running maximum in the march (caps the focusing gain at q_floor^-1/2)
    i_end = path.i_end
    # apex and reference point on the final leg
    last_event = max([s for s, _ in path.events] + [i_start])
    i_apex = last_event + int(np.argmax(z[last_event:i_end + 1]))
    if nonlinear_to_ground:
        i_ref = i_end
    else:
        i_ref = i_end
        zap = z[i_apex]
        for i in range(i_apex, i_end + 1):
            if z[i] <= zap - dz_ref:
                i_ref = i
                break
    t = np.asarray(t, float).copy()
    p = np.asarray(p, float).copy()
    t, p = _extend(t, p, p.size // 2, p.size // 2)
    qmax = abs(Q[i_start])
    max_dpp = float(np.max(np.abs(p)) / pa[i_start])
    n_caus_up = 0
    acc_ds = 0.0
    acc_lnH_pts = []
    i = i_start
    n_nl = 0
    sQ = np.sign(Q[i_start])
    while i < i_ref:
        j = i + 1
        dsig = float(path.sigma[j] - path.sigma[i])
        # nonlinear sub-steps over this interval (effective length Omega dsigma)
        if not linear and dsig > 0:
            dse = dsig * float(om[i])
            k = beta / (float(rho[i]) * float(c[i]) ** 3)
            dt = t[1] - t[0]
            done = 0.0
            while done < dse:
                pm = float(np.max(np.abs(p)))
                if pm <= 0:
                    break
                ds = min(dse - done, cfl * dt / (k * pm))
                a0, a1 = _span(p, 1e-10)
                a0, a1 = max(0, a0 - 8), min(p.size, a1 + 8)
                p[a0:a1] = burgers_step(p[a0:a1], k, dt, ds)
                done += ds
                n_nl += 1
        # spreading (Blokhintzev), with aerial-caustic floor and Hilbert at sign change
        qmax = max(qmax, abs(Q[j]))
        qi = max(abs(Q[i]), q_floor * qmax)
        qj = max(abs(Q[j]), q_floor * qmax)
        p *= math.sqrt((rho[j] * c[j] ** 2 * om[j] / qj) / (rho[i] * c[i] ** 2 * om[i] / qi))
        if np.sign(Q[j]) != sQ and path.Q_ok[j]:
            p = hilbert_shift(p)
            n_caus_up += 1
            sQ = np.sign(Q[j])
        # absorption accumulated over <= ds_lin of air path
        acc_lnH_pts.append((0.5 * (z[i] + z[j]), dsig, float(0.5 * (om[i] + om[j]))))
        acc_ds += dsig
        if acc_ds >= ds_lin or j == i_ref or j in ops:
            dt = t[1] - t[0]
            f = sfft.rfftfreq(p.size, dt)
            arr = np.array(acc_lnH_pts)
            wsum = max(float(arr[:, 1].sum()), 1e-12)
            zm = float(np.sum(arr[:, 0] * arr[:, 1]) / wsum)
            omm = float(np.sum(arr[:, 2] * arr[:, 1]) / wsum)
            # p (hence alpha) varies exponentially: use the log-mean of alpha over the chunk
            z_lo, z_hi = float(arr[:, 0].min()), float(arr[:, 0].max())
            st = ua.state(np.array(zm))
            lnH = absorption_exponent(f * omm, st, acc_ds)
            if z_hi - z_lo > 50.0:
                pl, ph = ua.interp_log("p", np.array(z_lo)), ua.interp_log("p", np.array(z_hi))
                ratio = float(np.log(pl / ph) / (pl / ph - 1.0) * pl / st["p"]) if pl != ph else 1.0
                # exact mean of 1/p over an exponential segment relative to 1/p(mid)
                lnH = lnH.real * (1.0 / ratio) + 1j * lnH.imag
            p = sfft.irfft(sfft.rfft(p) * np.exp(lnH), n=p.size)
            acc_lnH_pts = []
            acc_ds = 0.0
            t, p = _manage(t, p, n_span, fs_min)
        if j in ops:
            out = ops[j](t, p)
            t, p = out
            t, p = _manage(t, p, n_span, fs_min)
        max_dpp = max(max_dpp, float(np.max(np.abs(p)) / pa[j]))
        i = j
    # linear descent: normalise and filter
    p_ref_peak = float(np.max(np.abs(p)))
    h = p / math.sqrt(inv[i_ref])
    dt = t[1] - t[0]
    t, h = _extend(t, h, h.size // 2, h.size // 2)
    f = sfft.rfftfreq(h.size, dt)
    lnH = np.zeros(f.size, complex)
    # group the descent into ~ds_lin pieces
    seg_z = []
    acc = 0.0
    zs = []
    for ii in range(i_ref, i_end):
        dsig = float(path.sigma[ii + 1] - path.sigma[ii])
        zs.append((0.5 * (z[ii] + z[ii + 1]), dsig, float(om[ii])))
    if zs:
        zz = np.array([a for a, _, _ in zs])
        dd = np.array([b for _, b, _ in zs])
        oo = np.array([o for _, _, o in zs])
        nb = max(1, int(np.ceil(dd.sum() / ds_lin)))
        edges = np.linspace(0, len(zs), nb + 1).astype(int)
        for a, b in zip(edges[:-1], edges[1:]):
            if b <= a:
                continue
            w = dd[a:b]
            zm = float(np.sum(zz[a:b] * w) / max(w.sum(), 1e-12))
            omm = float(np.sum(oo[a:b] * w) / max(w.sum(), 1e-12))
            st = ua.state(np.array(zm))
            lnH += absorption_exponent(f * omm, st, float(w.sum()))
    h = sfft.irfft(sfft.rfft(h) * np.exp(lnH), n=h.size)
    # caustics on the descent (Hilbert per crossing)
    sgn = np.sign(Q[i_ref:i_end + 1])
    okd = path.Q_ok[i_ref:i_end + 1]
    n_down = int(np.sum((np.diff(sgn) != 0) & okd[1:] & okd[:-1]))
    for _ in range(n_down % 4):
        h = hilbert_shift(h)
    t, h = _manage(t, h, n_span, fs_min)
    E = float(np.sum(h * h) * (t[1] - t[0]))
    if stats is not None:
        stats.update(n_nl=n_nl)
    return MarchResult(t=t, h=h, E_h=E, max_dp_over_p=max_dpp, z_ref=float(z[i_ref]),
                       apex=float(z[i_apex]), n_caustic_up=n_caus_up, n_caustic_down=n_down,
                       p_ref_peak=p_ref_peak, info=dict(i_ref=i_ref, i_apex=i_apex,
                                                        n_nl=n_nl, fs_end=1.0 / (t[1] - t[0])))


def primary_signature(ac, phi, medium=None, fs=50e3, r0=None, **kw):
    """Primary boom (incident x 1.9) at the ground from this module's tracer + march
    (validation against boomsim.propagate.ground_signature)."""
    medium = Medium() if medium is None else medium
    path = ray_path(ac, phi, medium, "primary", dt=0.05)
    path, res = march_family(ac, phi, medium, "direct", path=path, fs=fs, r0=r0,
                             nonlinear_to_ground=True, **kw)
    ua = upper_atm()
    st = ua.state(np.array(0.0))
    q = abs(path.Q[path.i_end])
    return res.t, GROUND_FACTOR * res.h * math.sqrt(float(st["rho"]) * float(st["c"]) ** 2
                                                    * float(path.omega[path.i_end]) / q)


def ground_value(path: Path, res: MarchResult, medium, dY_dphi, ua=None):
    """Plain ray-theory ground waveform (incident, no 1.9) for a single ray given
    dY/dphi at landing (diverges at caustics; diagnostic only)."""
    ua = upper_atm() if ua is None else ua
    st = ua.state(np.array(0.0))
    c, u, v, om = _fields_full(medium, np.array([0.0]), path.sx, path.sy)
    Qg = abs(dY_dphi * path.landing["vz"])
    return res.t, res.h * math.sqrt(float(st["rho"]) * float(st["c"]) ** 2 * float(om[0]) / Qg)


# ==========================================================================================
# families: start waveform and operators
# ==========================================================================================
def start_index(path: Path, z_a, r0):
    r = np.hypot(path.y, path.z - z_a)
    i = int(np.argmax(r >= r0))
    return max(i, 1)


def shield_operator(stack, atm20, V_trace, which):
    """Return callable(t, p) applying R (which='R') or T (which='T', reference 'air')."""
    from . import shield as SH

    def H(f):
        T, R = SH.transfer(stack, f, V_trace, atm20, reference="air")
        return R if which == "R" else T

    def op(t, p):
        return _filter(t, p, H)

    return op


def scale_operator(k):
    def op(t, p):
        return t, k * p
    return op


def march_family(ac, phi, medium, family, shield=None, stack=None, fs=20e3, r0=None,
                 path=None, **kw):
    """Near field at azimuth phi -> march along the ``family`` path.  family in
    {'direct', 'indirect', 'shield', 'transmit'}:
      direct   unshielded geometry, first secondary landing without earlier events
      indirect unshielded geometry, ground reflection (x REFLECT_UP) then secondary landing
      shield   reflect-mode shield geometry, R(f) at the shield top
      transmit transmit-mode geometry, T(f) at the shield top (+ ground reflection if any)
    Returns (path, MarchResult) or (None, None) if no such landing."""
    from .atmosphere import std_atm
    r0 = 1.5 * ac.length if r0 is None else r0
    if family in ("direct", "indirect"):
        sh = None
        fam = family
    elif family == "shield":
        sh = Shield(z=shield.z, phi_max=shield.phi_max, mode="reflect")
        fam = "shield"
    elif family == "transmit":
        sh = Shield(z=shield.z, phi_max=shield.phi_max, mode="transmit")
        fam = None
    else:
        raise ValueError(family)
    if path is None:
        if fam is None:
            # transmit: the landing is classified as 'direct' (cut-off rays) or 'indirect'
            path = ray_path(ac, phi, medium, "indirect", shield=sh)
            if path is None:
                path = ray_path(ac, phi, medium, "direct", shield=sh)
        else:
            path = ray_path(ac, phi, medium, fam, shield=sh)
    if path is None:
        return None, None
    i0 = start_index(path, ac.altitude, r0)
    r_start = float(np.hypot(path.y[i0], path.z[i0] - ac.altitude))   # >= r0 (path sample)
    t, p = near_field(ac, r_start, fs, phi=abs(float(phi)))
    if phi < 0:
        pass  # symmetric source (lift dipole depends on cos phi only)
    t = t * float(path.omega[i0])
    ops = {}
    V_tr = 1.0 / math.hypot(path.sx, path.sy)
    atm20 = std_atm(shield.z) if shield is not None else None
    for st, kind in path.events:
        if kind == "ground":
            ops[st + 1] = scale_operator(REFLECT_UP)
        elif kind == "shield":
            ops[st + 1] = shield_operator(stack, atm20, V_tr, "R" if family == "shield" else "T")
    res = march(t, p, path, i0, ops=ops, **kw)
    return path, res


# ==========================================================================================
# ground binning
# ==========================================================================================
def branches(phi, Y, apex, gap=None):
    """Split a fine fan's landings of one family into branches: contiguous runs in phi
    (gap > 1.5 fan steps breaks), same apex class (< or >= 70 km), and monotonic Y(phi)
    (a sign change of dY/dphi = fold caustic starts a new branch)."""
    order = np.argsort(phi)
    phi, Y, apex = phi[order], Y[order], apex[order]
    dphi = np.median(np.diff(phi)) if phi.size > 1 else 1.0
    gap = 1.5 * dphi if gap is None else gap
    lab = np.zeros(phi.size, int)
    b = 0
    for i in range(1, phi.size):
        new = (phi[i] - phi[i - 1] > gap) or ((apex[i] >= 70e3) != (apex[i - 1] >= 70e3))
        if not new and i >= 2:
            d1 = Y[i - 1] - Y[i - 2]
            d2 = Y[i] - Y[i - 1]
            if d1 * d2 < 0 and (phi[i - 1] - phi[i - 2] <= gap):
                new = True
        if new:
            b += 1
        lab[i] = b
    out = np.empty_like(lab)
    out[order] = lab
    return out


def bin_ground(rows, marched, medium, dphi_fine, dY=3000.0, Y_edges=None, ua=None):
    """Caustic-regularised ground levels for one family.

    rows: fine-fan landing records (dicts with phi, Y, vz, apex, sx, sy, t, X).
    marched: dict phi -> MarchResult for marched rays of this family.
    Returns list of per-bin dicts (Y_center, branch, energy-equivalent waveform)."""
    ua = upper_atm() if ua is None else ua
    if not rows or not marched:
        return []
    phi = np.array([r["phi"] for r in rows])
    Y = np.array([r["Y"] for r in rows])
    ap = np.array([r["apex"] for r in rows])
    vz = np.array([r["vz"] for r in rows])
    sx = np.array([r["sx"] for r in rows])
    sy = np.array([r["sy"] for r in rows])
    tt = np.array([r["t"] for r in rows])
    XX = np.array([r["X"] for r in rows])
    st0 = ua.state(np.array(0.0))
    c0, _, u0, _, v0, _ = medium.fields(np.zeros(phi.size))
    om = 1 - u0 * sx - v0 * sy
    G = float(st0["rho"]) * float(st0["c"]) ** 2 * om / np.abs(vz)
    lab = branches(phi, Y, ap)
    mphi = np.array(sorted(marched))
    out = []
    if Y_edges is None:
        lo = math.floor(Y.min() / dY) * dY
        hi = math.ceil(Y.max() / dY) * dY + dY
        Y_edges = np.arange(lo, hi + 0.5 * dY, dY)
    for b in np.unique(lab):
        m = lab == b
        pb = phi[m]
        # marched rays inside this branch's phi-range (+- one fan step)
        lo, hi = pb.min() - 1.01 * dphi_fine, pb.max() + 1.01 * dphi_fine
        inb = mphi[(mphi >= lo) & (mphi <= hi)]
        # also accept marched rays of the same apex class within 2 deg if none inside
        if inb.size == 0:
            near = mphi[np.argsort(np.abs(mphi - pb.mean()))[:1]]
            if near.size and abs(near[0] - pb.mean()) < math.radians(2.0) + 0.5 * (hi - lo):
                inb = near
            else:
                continue
        logE = np.log(np.array([max(marched[q].E_h, 1e-300) for q in inb]))
        # upsample the branch (smooth within a branch) so that consecutive landings are
        # < dY/8 apart: avoids empty/double-counted bins when dY/dphi is large
        o = np.argsort(pb)
        pbs = pb[o]
        Yb, Gb = Y[m][o], G[m][o]
        apb = ap[m][o]
        dlb = (tt[m] - XX[m] * sx[m])[o]
        if pbs.size > 1:
            n_up = int(min(400, max(1, math.ceil(np.max(np.abs(np.diff(Yb))) / (dY / 8.0)))))
        else:
            n_up = 1
        dq = dphi_fine / n_up
        pu = pbs[0] - 0.5 * dphi_fine + (np.arange(pbs.size * n_up) + 0.5) * dq
        Yu = np.interp(pu, pbs, Yb)
        Gu = np.interp(pu, pbs, Gb)
        apu = np.interp(pu, pbs, apb)
        dlu = np.interp(pu, pbs, dlb)
        Ei = np.exp(np.interp(pu, inb, logE))
        w = Gu * Ei * dq                     # energy * width contribution [Pa^2 s m]
        ib = np.digitize(Yu, Y_edges) - 1
        for j in np.unique(ib):
            if j < 0 or j >= Y_edges.size - 1:
                continue
            mm = ib == j
            eps = float(np.sum(w[mm])) / dY                    # Pa^2 s per unit (incident)
            # representative waveform: nearest marched ray
            pc = float(np.mean(pu[mm]))
            q = inb[int(np.argmin(np.abs(inb - pc)))]
            r = marched[q]
            scale = math.sqrt(eps / r.E_h) if r.E_h > 0 else 0.0
            out.append(dict(Y=0.5 * (Y_edges[j] + Y_edges[j + 1]), branch=int(b),
                            phi_deg=math.degrees(pc), energy=eps * GROUND_FACTOR ** 2,
                            t=r.t, p=GROUND_FACTOR * scale * r.h, apex=float(np.mean(apu[mm])),
                            delay_s=float(np.mean(dlu[mm])), n_rays=int(mm.sum()),
                            phi_marched=q))
    return out
