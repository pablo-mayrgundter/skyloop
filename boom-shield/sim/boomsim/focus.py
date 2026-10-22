"""Focused sonic booms ("superbooms") from unsteady flight: ray caustics and focus amplitude.

Public API
----------
Column(atmosphere=None, z_min=-1000, z_max=30000, dz=10)
    Tabulated windless stratified column; ``integrals(s, z_lo, z_hi)`` returns the ray
    integrals I1 = int dz/s_z, I3 = int dz/s_z^3 and the travel time for horizontal
    slowness s (exact per cell for s_z^2 linear in z).
Trajectory(t, mach, z, psi=0.0, bank=0.0, column=None)
    Aircraft path from Mach(t), altitude z(t) and heading psi(t); horizontal position
    integrated.  Constructors: ``level_accel``, ``climb_accel``, ``level_turn``.
rays(traj, te, phi, z_obs)                 ray landing points, arrival times, Jacobian
focus_geometry(traj, te, phi, z_obs)       J and the geometric amplification G (finite
                                           away from caustics) relative to steady level
                                           flight on the same ray
scan(traj, te, phis, z_obs=0)              J on a (te, phi) grid; fold roots per phi
fold_parameters(traj, te, phi, z_obs=0)    local fold-caustic parameters at a ground focus:
                                           Lambda, Theta (from the two post-focus arrivals),
                                           caustic tilt alpha, L_c, R
caustic_altitude(traj, te, phi)            height where the ray touches the caustic
airy_focus(t, p, L_c, R, c)                linear fold-caustic (Airy) focus waveform and peak
nkz_focus(t, p, L_c, R, c, rho)            nonlinear (NKZ / nonlinear Tricomi) focus peak

Geometry
--------
Windless, horizontally stratified atmosphere.  At emission time t_e the aircraft at
r_a(t_e) moves with velocity v (|v| = V = M c(z_a), flight-path angle gamma, heading psi).
Wave normals lie on the Mach cone about v:
    n = (1/M) e_v + sqrt(1 - 1/M^2) (cos(phi) e_1 + sin(phi) e_2),
e_1 = the downward unit vector normal to e_v in the vertical plane, e_2 = e_v x e_1
(phi = 0 straight down; phi > 0 toward +y, the left of a +x heading).  The horizontal
slowness s_h = n_h/c(z_a) is conserved; the ray reaches altitude z after a horizontal
displacement s_h I1(|s_h|; z, z_a) and travel time int dz/(c^2 s_z).

Ray-tube area on a horizontal plane: J = det d(X, Y)/d(t_e, phi), with
    dX/dq = dr_h/dq + (ds_h/dq) I1 + s_h (|s_h| I3 d|s_h|/dq + (dz_a/dq)/s_z(z_a))
(I1 and I3 analytic; ds_h/dt_e by central differences of the smooth cone normal).  The
area normal to the ray is J |n_z|.  A caustic is J = 0; on the ground it is the focus
line.  The geometric amplitude relative to steady level flight on the same ray (same
horizontal slowness; M_eff, phi_eff) is, by the Blokhintzev invariant,
    G^2 = [J_ref(z)/J_ref(z_near)] / [J(z)/J(z_near)]
with both tube areas normalised just below the aircraft (z_near = z_a - 30 m).

Fold caustic amplitude
----------------------
Near a smooth fold the incident geometric field is p_ss (L_c/xi)^(1/4) at distance xi
from the caustic (lit side) and the two branches arrive Delta T = (4/3c) sqrt(2/R) xi^(3/2)
apart, R = relative radius of curvature of rays and caustic.  On the ground (tilt alpha
between caustic surface and ground, ground distance d from the focus line) this reads
G^4 d = Lambda = L_c/sin(alpha), Delta T = Theta d^(3/2); ``fold_parameters`` measures
Lambda and Theta from the two traced post-focus arrivals and alpha from the caustic
surface.  Linear theory (Airy uniform solution, e.g. Pierce 1989 sec. 9-4; Coulouvrat
2000; the time-domain form used by Plotkin's PCBoom) gives in the rfft convention
    P(xi, w) = P_ss(w) L_c^(1/4) 2 sqrt(pi) (2k^2/R)^(1/12) Ai(-(2k^2/R)^(1/3) xi) e^(+j pi/4),
which tends to the incident + (Hilbert-shifted) reflected arrivals far from the caustic.
Nonlinearity (Guiraud 1965; Gill & Seebass 1973; Marchiano, Coulouvrat & Grenon 2003)
enters through the nonlinear Tricomi / NKZ equation, solved by ``nkz_focus`` in the
normalised form P_XT = P_ZZ - Z P_TT + mu (P^2/2)_TT + Gamma P_TTT with
xi = delta Z, delta = (R lambda^2/2)^(1/3), x = ell X, ell = (2 R^2 lambda)^(1/3),
lambda = c T0, mu = beta ell/(rho c^2 lambda) per Pa.
"""
from __future__ import annotations

import math
from dataclasses import dataclass

import numpy as np
from scipy import fft as sfft
from scipy.interpolate import CubicSpline
from scipy.optimize import brentq
from scipy.special import airy

from .atmosphere import std_atm

PSF = 47.880258          # Pa per psf
G0 = 9.80665
BETA_AIR = 1.2


# =====================================================================================
# column integrals
# =====================================================================================
class Column:
    """Windless stratified column tabulated on a uniform z grid (``dz`` m).

    ``integrals(s, z_lo, z_hi)`` (broadcast arrays) returns dict(I1, I3, TT, ok, r_lo,
    r_hi): I1 = int_{z_lo}^{z_hi} dz/s_z, I3 = int dz/s_z^3, TT = int dz/(c^2 s_z) with
    s_z = sqrt(1/c^2 - s^2), each cell integrated exactly for s_z^2 linear in z (keeps the
    integrable singularity of grazing rays accurate; the I3 rule is the exact s-derivative
    of the I1 rule divided by s).  ``ok`` is False where s_z^2 <= 0 somewhere in
    [z_lo, z_hi] (ray turns before z_lo); r_lo, r_hi = s_z at the ends.
    """

    def __init__(self, atmosphere=None, z_min=-1000.0, z_max=30000.0, dz=10.0):
        self.atm_fn = std_atm if atmosphere is None else atmosphere
        n = int(round((z_max - z_min) / dz)) + 1
        self.z = np.linspace(z_min, z_max, n)
        self.dz = float(self.z[1] - self.z[0])
        a = self.atm_fn(self.z)
        self.c = np.asarray(a.c, dtype=float) * np.ones(n)
        self.rho = np.asarray(a.rho, dtype=float) * np.ones(n)
        self.inv_c2 = 1.0 / self.c ** 2

    def c_at(self, z):
        return np.interp(z, self.z, self.c)

    def rho_at(self, z):
        return np.interp(z, self.z, self.rho)

    def cmax_between(self, z_lo, z_hi):
        """max c over [z_lo, z_hi] (broadcast)."""
        z_lo, z_hi = np.broadcast_arrays(np.asarray(z_lo, float), np.asarray(z_hi, float))
        out = np.empty(z_lo.shape)
        # running maxima from the bottom and suffix structure: c is piecewise monotone, so
        # evaluate on the grid inside the interval plus both ends
        for idx in np.ndindex(z_lo.shape):
            a, b = z_lo[idx], z_hi[idx]
            m = (self.z > a) & (self.z < b)
            v = [self.c_at(a), self.c_at(b)]
            if m.any():
                v.append(self.c[m].max())
            out[idx] = max(v)
        return out

    def _q(self, s2, z):
        """s_z^2 at z (linear interpolation of 1/c^2, consistent with the cell rule)."""
        return np.interp(z, self.z, self.inv_c2) - s2

    def integrals(self, s, z_lo, z_hi, chunk=400):
        s, z_lo, z_hi = np.broadcast_arrays(np.asarray(s, float), np.asarray(z_lo, float),
                                            np.asarray(z_hi, float))
        shp = s.shape
        s, z_lo, z_hi = s.ravel(), z_lo.ravel(), z_hi.ravel()
        out = {k: np.full(s.size, np.nan) for k in ("I1", "I3", "TT", "r_lo", "r_hi")}
        ok_all = np.zeros(s.size, bool)
        zg, dz = self.z, self.dz
        for i0 in range(0, s.size, chunk):
            sl = slice(i0, min(s.size, i0 + chunk))
            ss, a, b = s[sl], z_lo[sl], z_hi[sl]
            s2 = ss[:, None] ** 2
            k_lo = np.clip(np.floor((a - zg[0]) / dz).astype(int), 0, zg.size - 2)
            k_hi = np.clip(np.floor((b - zg[0]) / dz).astype(int), 0, zg.size - 2)
            kmin, kmax = int(k_lo.min()), int(k_hi.max()) + 1
            zz = zg[kmin:kmax + 1]
            q = self.inv_c2[kmin:kmax + 1][None, :] - s2          # (n, K)
            # cells strictly inside: from k_lo+1 .. k_hi-1 (node index), partial at ends
            r = np.sqrt(np.maximum(q, 0.0))
            rs = r[:, :-1] + r[:, 1:]
            with np.errstate(divide="ignore", invalid="ignore"):
                c1 = 2.0 * dz / rs
                c3 = 2.0 * dz / (r[:, :-1] * r[:, 1:] * rs)
                cbar2 = (0.5 * (self.c[kmin:kmax] + self.c[kmin + 1:kmax + 1])) ** 2
                ct = c1 / cbar2[None, :]
            # cells touching a turning point are invalid (the ok flag rejects any ray whose
            # interval contains one); zero them so they cannot poison the cumulative sums
            bad = (q[:, :-1] <= 0.0) | (q[:, 1:] <= 0.0)
            c1 = np.where(bad, 0.0, c1)
            c3 = np.where(bad, 0.0, c3)
            ct = np.where(bad, 0.0, ct)
            cum1 = np.concatenate([np.zeros((ss.size, 1)), np.cumsum(c1, axis=1)], axis=1)
            cum3 = np.concatenate([np.zeros((ss.size, 1)), np.cumsum(c3, axis=1)], axis=1)
            cumt = np.concatenate([np.zeros((ss.size, 1)), np.cumsum(ct, axis=1)], axis=1)
            rows = np.arange(ss.size)
            jl, jh = k_lo - kmin, k_hi - kmin
            # partial cell [z_k, z] for z in cell k (q linear in the cell)
            ql = self._q(ss ** 2, a)
            qh = self._q(ss ** 2, b)
            rl, rh = np.sqrt(np.maximum(ql, 0.0)), np.sqrt(np.maximum(qh, 0.0))
            rkl, rkh = r[rows, jl], r[rows, jh]
            cl = 0.5 * (self.c_at(a) + self.c[k_lo])
            ch = 0.5 * (self.c_at(b) + self.c[k_hi])
            with np.errstate(divide="ignore", invalid="ignore"):
                pl1 = 2.0 * (a - zz[jl]) / (rkl + rl)
                ph1 = 2.0 * (b - zz[jh]) / (rkh + rh)
                pl3 = np.where(a > zz[jl], 2.0 * (a - zz[jl]) / (rkl * rl * (rkl + rl)), 0.0)
                ph3 = np.where(b > zz[jh], 2.0 * (b - zz[jh]) / (rkh * rh * (rkh + rh)), 0.0)
            I1 = (cum1[rows, jh] + ph1) - (cum1[rows, jl] + pl1)
            I3 = (cum3[rows, jh] + ph3) - (cum3[rows, jl] + pl3)
            TT = (cumt[rows, jh] + ph1 / ch ** 2) - (cumt[rows, jl] + pl1 / cl ** 2)
            # turning check: q must stay > 0 over [a, b]
            idx = np.arange(zz.size)[None, :]
            inside = (idx > jl[:, None]) & (idx <= jh[:, None])
            qmin = np.where(inside, q, np.inf).min(axis=1)
            ok = (np.minimum(np.minimum(ql, qh), qmin) > 0.0) & (b >= a)
            out["I1"][sl] = np.where(ok, I1, np.nan)
            out["I3"][sl] = np.where(ok, I3, np.nan)
            out["TT"][sl] = np.where(ok, TT, np.nan)
            out["r_lo"][sl] = rl
            out["r_hi"][sl] = rh
            ok_all[sl] = ok
        res = {k: v.reshape(shp) for k, v in out.items()}
        res["ok"] = ok_all.reshape(shp)
        return res

    def cutoff_mach(self, z_a, z_ground=0.0):
        """Undertrack cut-off Mach for level flight: max c over [z_ground, z_a] / c(z_a)."""
        return float(self.cmax_between(z_ground, z_a) / self.c_at(z_a))


_COLUMN = None


def default_column():
    global _COLUMN
    if _COLUMN is None:
        _COLUMN = Column()
    return _COLUMN


# =====================================================================================
# trajectory
# =====================================================================================
class Trajectory:
    """Aircraft path from tabulated Mach(t), altitude z(t) [m], heading psi(t) [rad]
    (and bank(t) [rad], positive = left wing down = left turn, used only to tilt the lift
    in the source).  Velocity |v| = M c(z); vertical speed dz/dt; horizontal speed
    sqrt(V^2 - vz^2) along psi; horizontal position integrated from (x0, y0).
    All quantities are cubic splines of t."""

    def __init__(self, t, mach, z, psi=0.0, bank=0.0, column=None, x0=0.0, y0=0.0, name=""):
        self.col = default_column() if column is None else column
        t = np.asarray(t, float)
        one = np.ones_like(t)
        mach = np.asarray(mach, float) * one
        z = np.asarray(z, float) * one
        psi = np.asarray(psi, float) * one
        bank = np.asarray(bank, float) * one
        self.name = name
        self.t = t
        c = self.col.c_at(z)
        V = mach * c
        zs = CubicSpline(t, z)
        vz = zs(t, 1)
        if np.any(np.abs(vz) >= V):
            raise ValueError("vertical speed exceeds airspeed")
        vh = np.sqrt(V ** 2 - vz ** 2)
        vx, vy = vh * np.cos(psi), vh * np.sin(psi)
        from scipy.integrate import cumulative_trapezoid
        x = x0 + cumulative_trapezoid(vx, t, initial=0.0)
        y = y0 + cumulative_trapezoid(vy, t, initial=0.0)
        self._sp = {k: CubicSpline(t, v) for k, v in
                    dict(x=x, y=y, z=z, mach=mach, psi=psi, bank=bank).items()}
        self.t0, self.t1 = float(t[0]), float(t[-1])

    # --- state ------------------------------------------------------------------
    def state(self, te):
        te = np.asarray(te, float)
        sp = self._sp
        z = sp["z"](te)
        c = self.col.c_at(z)
        M = sp["mach"](te)
        V = M * c
        vz = sp["z"](te, 1)
        psi = sp["psi"](te)
        vh = np.sqrt(np.maximum(V ** 2 - vz ** 2, 0.0))
        return dict(t=te, x=sp["x"](te), y=sp["y"](te), z=z, c=c, M=M, V=V, vz=vz,
                    vx=vh * np.cos(psi), vy=vh * np.sin(psi), psi=psi,
                    gamma=np.arcsin(np.clip(vz / V, -1, 1)), bank=sp["bank"](te))

    def time_at_mach(self, M):
        """First time the Mach number reaches M (NaN if never)."""
        tt = np.linspace(self.t0, self.t1, 20001)
        mm = self._sp["mach"](tt)
        i = np.flatnonzero((mm[:-1] - M) * (mm[1:] - M) <= 0)
        if i.size == 0:
            return float("nan")
        i = int(i[0])
        return float(brentq(lambda u: float(self._sp["mach"](u)) - M, tt[i], tt[i + 1]))

    # --- constructors -----------------------------------------------------------
    @classmethod
    def level_accel(cls, z, m0, m1, accel_g, column=None, pad_s=60.0, dt=0.25, psi=0.0,
                    name=""):
        """Level flight at altitude z: constant along-track acceleration accel_g [g]
        (negative = deceleration) from Mach m0 to m1, then pad_s of steady flight at m1
        (and pad_s steady at m0 before)."""
        col = default_column() if column is None else column
        c = float(col.c_at(z))
        a = accel_g * G0
        T = (m1 - m0) * c / a
        t = np.arange(-pad_s, T + pad_s + dt, dt)
        mach = np.clip(m0 + a * np.clip(t, 0, T) / c, min(m0, m1), max(m0, m1))
        # smooth the corners slightly (C1) so the spline velocity derivative is bounded
        return cls(t, mach, z, psi=psi, column=col, name=name or f"level {accel_g:+g} g")

    @classmethod
    def climb_accel(cls, z0, z1, m0, m1, accel_g, column=None, pad_s=60.0, dt=0.25,
                    level_off_s=60.0, name=""):
        """Climb from z0 to z1 while accelerating at accel_g [g] (airspeed rate) from
        Mach m0 to m1; altitude rises linearly in Mach (a constant-dM/dz schedule).  The
        altitude history is smoothed with a Gaussian of standard deviation level_off_s/4
        so that the pitch-over into the climb and the level-off at the top are gradual
        (an instantaneous level-off is an infinitely fast push-over and would focus)."""
        from scipy.ndimage import gaussian_filter1d
        col = default_column() if column is None else column
        zz = np.linspace(z0, z1, 2001)
        mm = np.linspace(m0, m1, zz.size)
        V = mm * col.c_at(zz)
        tt = (V - V[0]) / (accel_g * G0)
        pad = max(pad_s, level_off_s)
        t = np.arange(-pad, tt[-1] + pad + dt, dt)
        tc = np.clip(t, 0, tt[-1])
        z = np.interp(tc, tt, zz)
        if level_off_s > 0:
            z = gaussian_filter1d(z, level_off_s / 4.0 / dt, mode="nearest")
        mach = np.interp(tc, tt, mm)
        return cls(t, mach, z, column=col, name=name or f"climb {accel_g:g} g")

    @classmethod
    def level_turn(cls, z, mach, bank_deg, duration_s, column=None, pad_s=30.0, dt=0.25,
                   roll_in_s=5.0, name=""):
        """Level coordinated turn at constant Mach: turn rate g tan(bank)/V (left turn,
        heading increasing), entered after pad_s of straight flight with a linear roll-in
        of roll_in_s seconds ending at t = 0."""
        col = default_column() if column is None else column
        V = mach * float(col.c_at(z))
        b = math.radians(bank_deg)
        w = G0 * math.tan(b) / V
        t = np.arange(-pad_s, duration_s + dt, dt)
        # roll in over 5 s (smooth heading rate)
        ramp = (np.clip((t + roll_in_s) / roll_in_s, 0.0, 1.0) if roll_in_s > 0
                else (t >= 0).astype(float))
        rate = w * ramp
        from scipy.integrate import cumulative_trapezoid
        psi = cumulative_trapezoid(rate, t, initial=0.0)
        bank = b * ramp
        return cls(t, mach, z, psi=psi, bank=bank, column=col,
                   name=name or f"turn {bank_deg:g} deg")


# =====================================================================================
# rays and Jacobians
# =====================================================================================
def _frame(gamma, psi):
    cg, sg, cp, sp_ = np.cos(gamma), np.sin(gamma), np.cos(psi), np.sin(psi)
    ev = np.stack([cg * cp, cg * sp_, sg])
    e1 = np.stack([sg * cp, sg * sp_, -cg])
    e2 = np.cross(ev, e1, axis=0)
    return ev, e1, e2


def cone_normal(M, gamma, psi, phi):
    """Mach-cone wave normal (3, ...) and its phi-derivative."""
    ev, e1, e2 = _frame(gamma, psi)
    sb = np.sqrt(np.maximum(1.0 - 1.0 / M ** 2, 0.0))
    n = ev / M + sb * (np.cos(phi) * e1 + np.sin(phi) * e2)
    dn = sb * (-np.sin(phi) * e1 + np.cos(phi) * e2)
    return n, dn


def _slowness(traj, te, phi):
    st = traj.state(te)
    n, dn = cone_normal(st["M"], st["gamma"], st["psi"], phi)
    return st, n[:2] / st["c"], dn[:2] / st["c"], n[2]


def rays(traj: Trajectory, te, phi, z_obs=0.0, h=1e-3):
    """Rays emitted at te [s] with azimuth phi [rad] (broadcast) evaluated at altitude
    z_obs (broadcast).  Returns dict with P (2, ...) landing point, t_arr, J, dPdt, dPdp,
    s_h (2, ...), s, n_z (launch), ok, state."""
    te, phi, z_obs = np.broadcast_arrays(np.asarray(te, float), np.asarray(phi, float),
                                         np.asarray(z_obs, float))
    st, sh, dsh_dp, nz = _slowness(traj, te, phi)
    _, sh_p, _, _ = _slowness(traj, te + h, phi)
    _, sh_m, _, _ = _slowness(traj, te - h, phi)
    dsh_dt = (sh_p - sh_m) / (2 * h)
    s = np.hypot(sh[0], sh[1])
    ok0 = (st["M"] > 1.0) & (nz < 0.0) & (z_obs <= st["z"])
    I = traj.col.integrals(np.where(ok0, s, 0.0), z_obs, np.where(ok0, st["z"], z_obs + 1.0))
    ok = ok0 & I["ok"]
    I1, I3 = I["I1"], I["I3"]
    ds_dt = (sh[0] * dsh_dt[0] + sh[1] * dsh_dt[1]) / s
    ds_dp = (sh[0] * dsh_dp[0] + sh[1] * dsh_dp[1]) / s
    vr = np.stack([st["vx"], st["vy"]])
    P = np.stack([st["x"], st["y"]]) + sh * I1
    dPdt = vr + dsh_dt * I1 + sh * (s * I3 * ds_dt + st["vz"] / I["r_hi"])
    dPdp = dsh_dp * I1 + sh * (s * I3 * ds_dp)
    J = dPdt[0] * dPdp[1] - dPdt[1] * dPdp[0]
    return dict(P=P, t_arr=te + I["TT"], J=J, dPdt=dPdt, dPdp=dPdp, s_h=sh, s=s, I1=I1,
                I3=I3, ok=ok, state=st, nz_launch=nz)


def _reference_jacobian(st, sh, z_obs, col):
    """J of steady LEVEL flight (heading psi, same horizontal slowness s_h, i.e. Mach
    M_eff = 1/(c_a s_par) and azimuth phi_eff) at altitude z_obs."""
    psi = st["psi"]
    e = np.stack([np.cos(psi), np.sin(psi)])
    ep = np.stack([-np.sin(psi), np.cos(psi)])
    c = st["c"]
    s_par = sh[0] * e[0] + sh[1] * e[1]
    s_perp = sh[0] * ep[0] + sh[1] * ep[1]
    M_eff = 1.0 / (c * s_par)
    cmu = np.sqrt(np.maximum(1.0 - 1.0 / M_eff ** 2, 0.0))
    sphi = np.clip(s_perp * c / np.where(cmu > 0, cmu, 1.0), -1.0, 1.0)
    phi_eff = np.arcsin(sphi)
    # level cone derivative wrt phi_eff: d n_h/dphi = cmu * (0, cos phi) in body axes
    # (e_1 = -z has no horizontal part; e_2 = left)
    dsh = (cmu * np.cos(phi_eff) / c) * ep
    s = np.hypot(sh[0], sh[1])
    I = col.integrals(s, z_obs, st["z"])
    ds_dp = (sh[0] * dsh[0] + sh[1] * dsh[1]) / s
    dPdp = dsh * I["I1"] + sh * (s * I["I3"] * ds_dp)
    V = M_eff * c
    dPdt = V * e
    J = dPdt[0] * dPdp[1] - dPdt[1] * dPdp[0]
    return J, M_eff, phi_eff


def focus_geometry(traj: Trajectory, te, phi, z_obs=0.0, z_near=30.0):
    """J at z_obs and the geometric amplification G = p/p_ref (Blokhintzev; p_ref =
    steady level flight on the same ray), plus M_eff, phi_eff of that reference."""
    r = rays(traj, te, phi, z_obs)
    st = r["state"]
    zn = st["z"] - z_near
    rn = rays(traj, te, phi, zn)
    J_ref, M_eff, phi_eff = _reference_jacobian(st, r["s_h"], z_obs, traj.col)
    J_refn, _, _ = _reference_jacobian(st, r["s_h"], zn, traj.col)
    with np.errstate(divide="ignore", invalid="ignore"):
        G2 = (J_ref / J_refn) / (r["J"] / rn["J"])
    r.update(J_near=rn["J"], J_ref=J_ref, J_ref_near=J_refn, G2=G2,
             G=np.sqrt(np.abs(G2)), M_eff=M_eff, phi_eff=phi_eff)
    return r


# =====================================================================================
# caustic search
# =====================================================================================
def _edge(traj, ph, z_obs, t_bad, t_good, tol=1e-7):
    """Bisect the ray-cutoff boundary between t_bad (no ground ray) and t_good."""
    for _ in range(60):
        if abs(t_good - t_bad) < tol:
            break
        m = 0.5 * (t_bad + t_good)
        if bool(rays(traj, m, ph, z_obs)["ok"]):
            t_good = m
        else:
            t_bad = m
    return t_good


def scan(traj: Trajectory, te, phis, z_obs=0.0, refine_edges=True, along_phi=False):
    """J on the grid te x phis at altitude z_obs and the sign-change roots in te for
    each phi.  Near every cut-off boundary (rays start/stop reaching z_obs) the grid is
    refined geometrically, because for gentle accelerations the fold lies within a
    fraction of a second of the cut-off.  With ``along_phi`` sign changes along phi at
    each grid te are also returned (``roots_phi``, folds running along the track as in a
    steady turn).  Returns dict(te, phis, J (n_phi, n_te), ok, roots: list of
    (te*, phi), roots_phi, edges: list of (t_edge, phi))."""
    te = np.asarray(te, float)
    phis = np.asarray(phis, float)
    T, Pp = np.meshgrid(te, phis)
    r = rays(traj, T, Pp, z_obs)
    J = np.where(r["ok"], r["J"], np.nan)
    roots, edges = [], []
    for i, ph in enumerate(phis):
        tt, jj = te, J[i]
        if refine_edges:
            okr = np.isfinite(jj)
            extra = []
            for k in np.flatnonzero(okr[:-1] != okr[1:]):
                tb, tg = (te[k], te[k + 1]) if not okr[k] else (te[k + 1], te[k])
                tc = _edge(traj, ph, z_obs, tb, tg)
                edges.append((tc, float(ph)))
                span = abs(te[1] - te[0]) * 4
                sg = 1.0 if tg > tb else -1.0
                extra.append(tc + sg * np.geomspace(1e-6, span, 40))
            if extra:
                ex = np.concatenate(extra)
                rx = rays(traj, ex, ph, z_obs)
                jx = np.where(rx["ok"], rx["J"], np.nan)
                tt = np.concatenate([te, ex])
                jj = np.concatenate([J[i], jx])
                o = np.argsort(tt)
                tt, jj = tt[o], jj[o]
        good = np.isfinite(jj[:-1]) & np.isfinite(jj[1:])
        idx = np.flatnonzero(good & (np.sign(jj[:-1]) != np.sign(jj[1:])))
        for k in idx:
            f = lambda u: float(rays(traj, u, ph, z_obs)["J"])
            try:
                roots.append((float(brentq(f, tt[k], tt[k + 1], xtol=1e-10)), float(ph)))
            except ValueError:
                pass
    # sign changes along phi at fixed te (folds that run along the track, e.g. turns)
    roots_phi = []
    if along_phi:
        for k, t_ in enumerate(te):
            jj = J[:, k]
            good = np.isfinite(jj[:-1]) & np.isfinite(jj[1:])
            idx = np.flatnonzero(good & (np.sign(jj[:-1]) != np.sign(jj[1:])))
            for i in idx:
                f = lambda u: float(rays(traj, t_, u, z_obs)["J"])
                try:
                    roots_phi.append((float(t_), float(brentq(f, phis[i], phis[i + 1],
                                                              xtol=1e-10))))
                except ValueError:
                    pass
    return dict(te=te, phis=phis, J=J, ok=r["ok"], roots=roots, roots_phi=roots_phi,
                edges=edges, P=r["P"])


def caustic_altitude(traj: Trajectory, te, phi, z_lo=-900.0, z_hi=None, n=400):
    """Altitude z_c where the ray (te, phi) touches the caustic (J(z) = 0), searching
    from just below the aircraft down to z_lo (default -900 m: extends a little below
    the ground so that rays whose caustic is 'virtual' still have a root).  NaN if J
    keeps its sign.  Returns the highest root (first caustic met)."""
    st = traj.state(te)
    zt = float(st["z"]) - 5.0 if z_hi is None else float(z_hi)
    zz = np.linspace(zt, z_lo, n)
    r = rays(traj, te, phi, zz)
    J = np.where(r["ok"], r["J"], np.nan)
    good = np.isfinite(J[:-1]) & np.isfinite(J[1:])
    idx = np.flatnonzero(good & (np.sign(J[:-1]) != np.sign(J[1:])))
    if idx.size == 0:
        return float("nan")
    k = int(idx[0])
    f = lambda z: float(rays(traj, te, phi, z)["J"])
    return float(brentq(f, zz[k + 1], zz[k], xtol=1e-6))


def _newton_land(traj, q0, target, z_obs, tol=1e-7, it=40):
    q = np.array(q0, float)
    for _ in range(it):
        r = rays(traj, q[0], q[1], z_obs)
        F = r["P"] - target
        if not r["ok"]:
            return None
        if np.hypot(*F) < tol:
            return q, r
        A = np.array([[r["dPdt"][0], r["dPdp"][0]], [r["dPdt"][1], r["dPdp"][1]]])
        try:
            dq = np.linalg.solve(A, -F)
        except np.linalg.LinAlgError:
            return None
        # damp large steps
        lim = np.array([2.0, 0.05])
        sc = np.max(np.abs(dq) / lim)
        if sc > 1:
            dq /= sc
        q = q + dq
    r = rays(traj, q[0], q[1], z_obs)
    if r["ok"] and np.hypot(*(r["P"] - target)) < 1e-4:
        return q, r
    return None


def fold_parameters(traj: Trajectory, te, phi, z_obs=0.0, d_list=(5.0, 10.0, 20.0, 40.0)):
    """Local fold-caustic parameters at the focus point (te, phi) on the plane z_obs.

    Returns dict(X, Y, N_g (lit-side unit normal of the focus line), Lambda [m],
    Theta [s/m^1.5], alpha [rad], L_c [m], R [m], c, rho, M, M_eff, phi_eff,
    branches = per-d diagnostics, t_arr, s_h)."""
    col = traj.col
    r0 = rays(traj, te, phi, z_obs)
    A = np.array([[r0["dPdt"][0], r0["dPdp"][0]], [r0["dPdt"][1], r0["dPdp"][1]]],
                 dtype=float)
    U, S, Vt = np.linalg.svd(A)
    e_null = Vt[1]                      # in (te, phi) space
    Ng = U[:, 1]                        # ground normal of the focus line (perp. to range)
    # second derivative along the null direction (finite differences of the map)
    scale_t = 1.0
    hq = 1e-2 * scale_t
    Pp = rays(traj, te + hq * e_null[0], phi + hq * e_null[1], z_obs)["P"]
    Pm = rays(traj, te - hq * e_null[0], phi - hq * e_null[1], z_obs)["P"]
    D2 = (Pp + Pm - 2 * r0["P"]) / hq ** 2
    curv = float(Ng @ D2)
    if curv < 0:
        Ng, curv = -Ng, -curv
    P0 = r0["P"]
    out_b = []
    for d in d_list:
        eps = math.sqrt(2 * d / curv)
        target = P0 + d * Ng
        sols = []
        for sgn in (-1, 1):
            q0 = (te + sgn * eps * e_null[0], phi + sgn * eps * e_null[1])
            res = _newton_land(traj, q0, target, z_obs)
            if res is None:
                break
            sols.append(res)
        if len(sols) < 2 or abs(sols[0][0][0] - sols[1][0][0]) < 1e-9:
            continue
        g1 = focus_geometry(traj, sols[0][0][0], sols[0][0][1], z_obs)
        g2 = focus_geometry(traj, sols[1][0][0], sols[1][0][1], z_obs)
        dT = abs(float(g1["t_arr"] - g2["t_arr"]))
        G1, G2 = float(g1["G"]), float(g2["G"])
        out_b.append(dict(d=d, dT=dT, G1=G1, G2=G2, Lam=d * G1 ** 2 * G2 ** 2,
                          Theta=dT / d ** 1.5))
    if not out_b:
        raise RuntimeError("fold_parameters: could not find the two post-focus arrivals")
    # extrapolate to d -> 0 (linear in sqrt(d)) if more than one point
    ds = np.array([b["d"] for b in out_b])
    Lam = np.array([b["Lam"] for b in out_b])
    The = np.array([b["Theta"] for b in out_b])
    if ds.size >= 2:
        x = np.sqrt(ds)
        Lam0 = float(np.polyval(np.polyfit(x, Lam, 1), 0.0))
        The0 = float(np.polyval(np.polyfit(x, The, 1), 0.0))
        if not (Lam0 > 0 and The0 > 0):
            Lam0, The0 = float(Lam[0]), float(The[0])
    else:
        Lam0, The0 = float(Lam[0]), float(The[0])
    # caustic surface tilt: caustic points of neighbouring rays
    hs = [(2e-2, 0.0), (0.0, 2e-3)]
    pts = {}
    for key, (dt_, dp_) in zip("tp", hs):
        pp = []
        for sg in (-1, 1):
            tt, ph = te + sg * dt_, phi + sg * dp_
            zc = caustic_altitude(traj, tt, ph, z_lo=z_obs - 900.0)
            if not np.isfinite(zc):
                pp = None
                break
            rr = rays(traj, tt, ph, zc)
            pp.append(np.array([rr["P"][0], rr["P"][1], zc]))
        pts[key] = pp
    alpha = float("nan")
    vecs = [pts[k][1] - pts[k][0] for k in "tp" if pts[k] is not None]
    if len(vecs) == 2:
        Nv = np.cross(vecs[0], vecs[1])
        Nv /= np.linalg.norm(Nv)
        alpha = float(math.acos(min(1.0, abs(Nv[2]))))
    elif len(vecs) == 1:
        # 2-D fallback: caustic curve along the t-family, ground normal direction
        v = vecs[0]
        alpha = float(math.atan2(abs(v[2]), abs(float(Ng @ v[:2]))))
    c = float(col.c_at(z_obs))
    rho = float(col.rho_at(z_obs))
    sa = math.sin(alpha)
    L_c = Lam0 * sa
    R = 2.0 * (4.0 * sa ** 1.5 / (3.0 * c * The0)) ** 2
    g0 = focus_geometry(traj, te, phi, z_obs)
    st = r0["state"]
    return dict(te=float(te), phi=float(phi), X=float(P0[0]), Y=float(P0[1]), N_g=Ng,
                Lambda=Lam0, Theta=The0, alpha=alpha, L_c=L_c, R=R, c=c, rho=rho,
                M=float(st["M"]), M_eff=float(g0["M_eff"]), phi_eff=float(g0["phi_eff"]),
                z_a=float(st["z"]), t_arr=float(r0["t_arr"]), s=float(r0["s"]),
                s_h=r0["s_h"].astype(float), branches=out_b,
                x_a=float(st["x"]), y_a=float(st["y"]), bank=float(st["bank"]),
                psi=float(st["psi"]))


# =====================================================================================
# focus amplitude: linear Airy (fold caustic)
# =====================================================================================
def _airy_kernel(f, xi, L_c, R, c):
    """Complex fold-caustic transfer (rfft convention) at distance xi [m] from the
    caustic (xi > 0 lit side), relative to the incident steady waveform."""
    k = 2 * np.pi * np.asarray(f, float) / c
    kk = np.where(k > 0, k, 1.0)
    a = (2 * kk ** 2 / R) ** (1.0 / 3.0)
    ai = airy(-a * xi)[0]
    H = L_c ** 0.25 * 2 * np.sqrt(np.pi) * a ** 0.25 * ai * np.exp(1j * np.pi / 4)
    # f = 0: the caustic does not amplify a static pressure (limit is singular in
    # theory; the waveforms here carry no mean) -> 0
    return np.where(k > 0, H, 0.0)


def airy_focus(t, p, L_c, R, c, xi=None, pad=4, return_waveform=True):
    """Linear focus of the free-field incident waveform (t, p) at a fold caustic.

    Scans distances xi [m] (default: log-spaced over the Airy layers of the waveform's
    band, both sides of the caustic) and returns dict(peak, xi_peak, t, p (waveform at
    the peak), xi, peaks (max p vs xi), neg (min p vs xi))."""
    t = np.asarray(t, float)
    p = np.asarray(p, float)
    dt = t[1] - t[0]
    n = p.size
    nfft = sfft.next_fast_len(pad * n, real=True)
    off = (nfft - n) // 2
    buf = np.zeros(nfft)
    buf[off:off + n] = p
    Pw = sfft.rfft(buf)
    f = sfft.rfftfreq(nfft, dt)
    if xi is None:
        # Airy layer thickness for 1 Hz ... 500 Hz
        lam = c / np.array([500.0, 1.0])
        dl = (R * lam ** 2 / 2) ** (1 / 3)
        xi = np.concatenate([-np.geomspace(dl[0] * 0.05, dl[1] * 1.5, 40)[::-1],
                             [0.0], np.geomspace(dl[0] * 0.05, dl[1] * 8.0, 160)])
    xi = np.asarray(xi, float)
    peaks = np.empty(xi.size)
    negs = np.empty(xi.size)
    best = (-np.inf, None, None)
    for i, x in enumerate(xi):
        q = sfft.irfft(Pw * _airy_kernel(f, x, L_c, R, c), nfft)
        peaks[i] = q.max()
        negs[i] = q.min()
        if peaks[i] > best[0]:
            best = (peaks[i], x, q)
    out = dict(peak=float(best[0]), xi_peak=float(best[1]), xi=xi, peaks=peaks, neg=negs)
    if return_waveform:
        q = best[2]
        tt = t[0] + (np.arange(nfft) - off) * dt
        out.update(t=tt, p=q)
    return out


def incident_at(t, p, L_c, xi):
    """Geometric-acoustics incident amplitude at distance xi (diagnostic)."""
    return (L_c / xi) ** 0.25 * np.asarray(p)


# =====================================================================================
# focus amplitude: nonlinear (NKZ / nonlinear Tricomi), normalised marching solver
# =====================================================================================
def _muscl_burgers(P, k, dT, dX):
    """MUSCL-Hancock + Engquist-Osher step of P_X = k P P_T along the last axis
    (same scheme as propagate.burgers_step, vectorised over rows)."""
    lam = dX / dT
    ue = np.concatenate([P[:, :1], P[:, :1], P, P[:, -1:], P[:, -1:]], axis=1)
    d = np.diff(ue, axis=1)
    d1, d2 = d[:, :-1], d[:, 1:]
    prod = d1 * d2
    s = np.zeros_like(prod)
    np.divide(2.0 * prod, d1 + d2, out=s, where=prod > 0.0)
    uc = ue[:, 1:-1]
    half = (0.5 * lam * k) * s * uc
    um = uc - 0.5 * s + half
    up = uc + 0.5 * s + half
    uL, uR = up[:, :-1], um[:, 1:]
    F = -0.5 * k * (np.minimum(uL, 0.0) ** 2 + np.maximum(uR, 0.0) ** 2)
    return P - lam * (F[:, 1:] - F[:, :-1])


def nkz_focus(t, p, L_c, R, c, rho, beta=BETA_AIR, T0=None, n_rise=4.0, omega_cut=4.0,
              dT=None, dZ=None, Z_top=None, Z_bot=6.0, X_end=4.0, rise_s=None, linear=False,
              cfl=0.4, dX_max=0.02, record_from=None, verbose=False, workers=1):
    """Nonlinear fold-caustic focus by marching the normalised NKZ equation
    P_XT = P_ZZ - Z P_TT + mu (P^2/2)_TT + Gamma P_TTT in X from the linear (Airy)
    steady solution of the incident waveform (t, p) [free field, Pa].

    Units: T0 [s] (default: incident 10-90 % rise time / n_rise), lambda = c T0,
    delta = (R lambda^2/2)^(1/3), ell = (2 R^2 lambda)^(1/3), mu = beta ell/(rho c^2 lambda)
    per Pa.  Gamma keeps the incident shock at its rise time (Taylor balance at the
    incident amplitude at Z = 1).  Returns dict(peak_lin, peak_nl, ratio, Z_peak, xi_peak,
    hist (X, max), t, p (waveform at the peak), mu, Gamma, delta, grid)."""
    t = np.asarray(t, float)
    p = np.asarray(p, float)
    dt = t[1] - t[0]
    pk = p.max()
    if rise_s is None:
        # shock rise from the steepest compression (0.8 peak / max slope = the 10-90 %
        # time of a linear ramp); robust to slow precursors that fool a 10-90 % search
        smax = float(np.max(np.diff(p))) / dt
        rise_s = 0.8 * pk / smax if smax > 0 else 1e-3
    if T0 is None:
        T0 = max(5e-5, rise_s / n_rise)
    if dT is None:
        dT = math.pi / (1.3 * omega_cut)
    lam = c * T0
    delta = (R * lam ** 2 / 2.0) ** (1.0 / 3.0)
    ell = (2.0 * R ** 2 * lam) ** (1.0 / 3.0)
    mu = 0.0 if linear else beta * ell / (rho * c ** 2 * lam)
    amp = (L_c / delta) ** 0.25          # incident at Z = 1 is amp * p
    # resample incident onto the T grid (T in units of T0)
    Tn = (t - t[0]) / T0
    dur = Tn[-1]
    # find signal support
    sup = np.flatnonzero(np.abs(p) > 1e-3 * np.abs(p).max())
    T_a, T_b = Tn[sup[0]], Tn[sup[-1]]
    if Z_top is None:
        # low-frequency Airy layer: Omega_min^(-2/3), Omega_min ~ 2 pi/(signal length)
        Om = 2 * np.pi / max(T_b - T_a, 10.0)
        Z_top = float(np.clip(4.0 * Om ** (-2.0 / 3.0) + 12.0, 20.0, 60.0))
    spread = (2.0 / 3.0) * Z_top ** 1.5
    NT = sfft.next_fast_len(int(np.ceil(((T_b - T_a) + 2 * spread + 60.0) / dT)), real=True)
    Tg = np.arange(NT) * dT
    T_c = 0.5 * (Tg[-1]) - 0.5 * (T_b - T_a)       # place incident centred
    f_in = np.interp(Tg, Tn - T_a + T_c, p, left=0.0, right=0.0) * amp
    Om = 2 * np.pi * sfft.rfftfreq(NT, dT)
    sponge = 6.0
    if dZ is None:
        dZ = 0.8 * math.pi / (omega_cut * math.sqrt(Z_top + sponge))
    Zg = np.arange(-Z_bot - sponge, Z_top + sponge, dZ)
    NZ = sfft.next_fast_len(Zg.size)
    Zg = -Z_bot - sponge + np.arange(NZ) * dZ
    kz = 2 * np.pi * sfft.fftfreq(NZ, dZ)
    # initial linear steady solution: P(Z, Om) = F(Om) 2 sqrt(pi) Om^(1/6) Ai(-Om^(2/3) Z) e^{j pi/4}
    Fw = sfft.rfft(f_in)
    lp = np.clip((Om - 0.6 * omega_cut) / (0.4 * omega_cut), 0.0, 1.0)
    Fw = Fw * np.cos(0.5 * np.pi * lp) ** 2
    Omp = np.where(Om > 0, Om, 1.0)
    ai = airy(-np.outer(Zg, Omp ** (2.0 / 3.0)))[0]
    H = 2 * np.sqrt(np.pi) * Omp ** (1.0 / 6.0) * np.exp(1j * np.pi / 4)
    Pw = ai * (H * Fw)[None, :]
    Pw[:, 0] = 0.0
    # The linear steady field P_ref (X-independent) is the incident supply: only the
    # nonlinear departure D = P - P_ref is marched by the linear operator and damped in
    # the sponges (so the incident keeps arriving from beyond Z_top, and nothing wraps).
    def ramp(x):
        return 0.5 * (1 + np.cos(np.pi * np.clip(x, 0, 1)))
    top_w = ramp((Zg - Z_top) / sponge)
    bot_w = ramp((-Z_bot - Zg) / sponge)
    sig = 3.0 * ((1 - top_w) + (1 - bot_w))
    P_ref = sfft.irfft(Pw, NT, axis=1, workers=workers)
    del Pw, ai
    # absorption: Gamma from Taylor balance at incident half-jump a (at Z = 1)
    a_half = 0.5 * amp * pk
    rise_n = rise_s / T0
    Gamma = 0.0 if linear else mu * a_half * rise_n / 4.39
    # Taylor balance at every height: the incident amplitude falls as Z^(-1/4) away from
    # the caustic, so Gamma(Z) = Gamma min(1, Z^(-1/4)) keeps the incident shock at its
    # rise time all the way in (a uniform Gamma would over-smooth the far incident)
    GZ = Gamma * np.minimum(1.0, np.maximum(Zg, 1.0) ** -0.25)
    P_ref_w = sfft.rfft(P_ref, axis=1, workers=workers) if Gamma > 0 else None

    def lin_factors(dX):
        diff = np.exp(1j * np.outer(kz ** 2, 1.0 / Omp) * dX)
        diff[:, 0] = 0.0
        refr_half = np.exp(-1j * np.outer(Zg, Om) * dX / 2)
        absn = np.exp(-np.outer(GZ, Om ** 2) * dX)
        damp = np.exp(-sig * dX)[:, None]
        dref = (sfft.irfft(P_ref_w * (absn - 1.0), NT, axis=1, workers=workers)
                if Gamma > 0 else 0.0)
        return diff, refr_half, absn, damp, dref
    lin_peak0 = float(P_ref[(Zg >= -Z_bot) & (Zg <= 6.0)].max())
    D = np.zeros_like(P_ref)
    X = 0.0
    hist = []
    cache = {}
    rec_from = X_end - 1.0 if record_from is None else record_from
    best = (-np.inf, None, None)
    zsel = (Zg >= -2.0) & (Zg <= 6.0)
    pmax0 = float(np.abs(P_ref).max())
    while X < X_end - 1e-12:
        dX = dX_max
        if mu > 0:
            pmax = max(pmax0, float(np.abs(P_ref + D).max()))
            lim = cfl * dT / (mu * pmax)
            # step ladder dX_max 2^(-k/2) so the linear factors are reused (cache hits)
            k = max(0, int(math.ceil(2.0 * math.log2(dX_max / lim)))) if lim < dX_max else 0
            dX = dX_max * 2.0 ** (-0.5 * k)
        dX = min(dX, X_end - X)
        key = round(dX, 12)
        if key not in cache:
            if len(cache) > 2:
                cache.clear()
            cache[key] = lin_factors(dX)
        diff, refr_half, absn, damp, dref = cache[key]
        Dw = sfft.rfft(D, axis=1, workers=workers)
        Dw *= refr_half
        Dk = sfft.fft(Dw, axis=0, workers=workers)
        Dk *= diff
        Dw = sfft.ifft(Dk, axis=0, workers=workers)
        Dw *= refr_half
        Dw *= absn
        D = sfft.irfft(Dw, NT, axis=1, workers=workers)
        D += dref
        D *= damp
        if mu > 0:
            D = _muscl_burgers(P_ref + D, mu, dT, dX) - P_ref
        X += dX
        sub = P_ref[zsel] + D[zsel]
        m = float(sub.max())
        hist.append((X, m))
        if X >= rec_from and m > best[0]:
            iz = int(np.argmax(sub.max(axis=1)))
            best = (m, float(Zg[zsel][iz]), sub[iz].copy())
        if verbose and len(hist) % 50 == 0:
            print(f"  X={X:.3f} max={m:.3f} dX={dX:.4f}", flush=True)
    peak = best[0]
    tt = Tg * T0
    return dict(peak_lin=lin_peak0, peak_nl=float(peak), ratio=float(peak / lin_peak0),
                Z_peak=best[1], xi_peak=best[1] * delta if best[1] is not None else np.nan,
                hist=np.array(hist), t=tt, p=best[2], mu=mu, Gamma=Gamma, delta=delta,
                ell=ell, T0=T0, amp=amp,
                grid=dict(NZ=NZ, NT=NT, dZ=dZ, dT=dT, Z_top=Z_top))
