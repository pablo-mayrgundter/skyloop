"""Generic 2-D multi-domain boundary element method (collocation, piecewise constant).

Used by :mod:`boomsim.crossplane` for the low-frequency cross-plane solve, but kept
free of any BoomShield specifics.

Problem class
-------------
Domains: one unbounded exterior (Helmholtz, wavenumber k, outgoing; may contain
the incident-field sources) and any number of bounded interior domains, each either
Helmholtz (wavenumber k_d) or modified Helmholtz (decay constant kappa_d)::

    Helmholtz           (lap + k^2) p = 0,      G = -(j/4) H0^(2)(k r)
    modified Helmholtz  (lap - kappa^2) p = 0,  G = K0(kappa r) / (2 pi)

(time factor e^{+j omega t}; both G satisfy lap G +- ... = -delta and share the
Laplace singularity -ln(r)/2pi). Normal velocity in domain d: v = -(1/(j omega rho_d))
dp/dn.

Boundaries are straight panels with unit normal ``n`` pointing from side A into side B.
Two panel types:

* membrane panel (A, B different domains): normal velocity v continuous,
  ``p_A - p_B = Z v - f_tie(v)`` (Z = series impedance, e.g. j omega m; f_tie = extra
  normal force per unit area that is linear in the panel velocities, e.g. ties).
* sheet panel (exterior air on both sides, zero thickness): v continuous,
  ``p_A - p_B = Z v``.  Represented by its pressure jump mu = Z v (double layer).

Formulation
-----------
Direct boundary integral equations, collocated at panel midpoints:

    1/2 p + K[p] - S[q] = p_inc      (CBIE, every domain; q = dp/d(nu), nu outward)

and for the exterior the Burton-Miller combination CBIE + alpha * HBIE with
``alpha = -j/k`` (``-j/k`` is the e^{+j omega t} image of the classical i/k) which removes
the fictitious interior resonances of the exterior equation:

    1/2 q - A[q] + H[p] = dp_inc/d(nu)  (HBIE)

Sheet panels add -K[mu] to the exterior representation and get the hypersingular
equation ``-j omega rho v = dp_inc/dn + (A[q] - H[p])_body - H[mu]_sheet``.
Unknowns per membrane panel: p_A and v (p_B eliminated with the membrane law); per
sheet panel one scaled velocity.

Integration: every kernel is split into its Laplace part, integrated analytically
over the straight panel for any target point (log, double layer, its gradient and
the finite-part hypersingular integral), plus the bounded/log-singular difference
kernel, integrated with Gauss-Legendre (16 points; on the panel itself with a
quadratic map that removes the log singularity). Far pairs (centre distance > 3 L)
use 2-point and intermediate pairs 6-point Gauss on the full kernel.
"""
from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np
from scipy import special as sp

TWO_PI = 2.0 * np.pi


# ------------------------------------------------------------------- kernels
class Helmholtz:
    """G = -(j/4) H0^(2)(k r): outgoing for e^{+j omega t}."""

    def __init__(self, k: float):
        self.k = float(k)

    def radial(self, r):
        k = self.k
        z = k * r
        H0, H1 = hankel2_01(z)
        G = -0.25j * H0
        G1 = 0.25j * k * H1
        G2 = 0.25j * k * k * (H0 - H1 / z)
        return G, G1, G2


class ModHelmholtz:
    """G = K0(kappa r) / (2 pi)."""

    def __init__(self, kappa: float):
        self.kappa = max(float(kappa), 1e-12)
        self.cutoff = 45.0 / self.kappa      # K0(45) ~ 6e-21

    def radial(self, r):
        a = self.kappa
        z = a * r
        K0, K1 = sp.k0(z), sp.k1(z)
        G = K0 / TWO_PI
        G1 = -a * K1 / TWO_PI
        G2 = a * a * (K0 + K1 / z) / TWO_PI
        return G + 0j, G1 + 0j, G2 + 0j


def hankel2_01(z, z_asym=25.0):
    """H0^(2)(z), H1^(2)(z) for real z > 0: Bessel functions below z_asym, Hankel's
    asymptotic expansion (5 terms, relative error < 3e-9 for z >= 25) above."""
    z = np.asarray(z, float)
    H0 = np.empty(z.shape, complex)
    H1 = np.empty(z.shape, complex)
    lo = z < z_asym
    zl = z[lo]
    H0[lo] = sp.j0(zl) - 1j * sp.y0(zl)
    H1[lo] = sp.j1(zl) - 1j * sp.y1(zl)
    zh = z[~lo]
    if zh.size:
        u = 1.0 / (8.0 * zh)
        ph = np.exp(-1j * (zh - 0.25 * np.pi)) * np.sqrt(2.0 / (np.pi * zh))
        for n, H in ((0, H0), (1, H1)):
            mu = 4.0 * n * n
            # sum_k (-j)^k a_k(n) / z^k with a_k = prod_{m=1..k} (mu - (2m-1)^2) / (k! 8^k)
            s = np.ones_like(zh, dtype=complex)
            term = np.ones_like(zh, dtype=complex)
            for kk in range(1, 6):
                term = term * (mu - (2 * kk - 1) ** 2) * u / kk * (-1j)
                s = s + term
            H[~lo] = ph * s * (1j if n == 1 else 1.0)
    return H0, H1


def _laplace_radial(r):
    return -np.log(r) / TWO_PI, -1.0 / (TWO_PI * r), 1.0 / (TWO_PI * r * r)


# ------------------------------------------------------------------- panels
@dataclass
class Panels:
    a: np.ndarray   # (N, 2) start points
    b: np.ndarray   # (N, 2) end points

    def __post_init__(self):
        self.a = np.asarray(self.a, float).reshape(-1, 2)
        self.b = np.asarray(self.b, float).reshape(-1, 2)
        d = self.b - self.a
        self.L = np.hypot(d[:, 0], d[:, 1])
        if np.any(self.L <= 0):
            raise ValueError("zero-length panel")
        self.t = d / self.L[:, None]
        # normal = tangent rotated by -90 deg: (t_y, -t_x)  (caller orients panels)
        self.n = np.stack([self.t[:, 1], -self.t[:, 0]], axis=1)
        self.c = 0.5 * (self.a + self.b)

    def __len__(self):
        return self.L.size


_GL = {m: np.polynomial.legendre.leggauss(m) for m in (2, 6, 16)}


def _kernel_values(kern, rvec, nx, ny, which, diff=False):
    """Kernel values at separation rvec = x - y (..., 2). nx/ny broadcastable (..., 2)."""
    r = np.hypot(rvec[..., 0], rvec[..., 1])
    G, G1, G2 = kern.radial(r)
    if diff:
        gL, g1L, g2L = _laplace_radial(r)
        G, G1, G2 = G - gL, G1 - g1L, G2 - g2L
    out = {}
    rhx = rhy = None
    if "K" in which or "H" in which:
        rhy = (rvec[..., 0] * ny[..., 0] + rvec[..., 1] * ny[..., 1]) / r
    if "A" in which or "H" in which:
        rhx = (rvec[..., 0] * nx[..., 0] + rvec[..., 1] * nx[..., 1]) / r
    if "S" in which:
        out["S"] = G
    if "K" in which:
        out["K"] = -G1 * rhy
    if "A" in which:
        out["A"] = G1 * rhx
    if "H" in which:
        nn = nx[..., 0] * ny[..., 0] + nx[..., 1] * ny[..., 1]
        out["H"] = -G2 * rhx * rhy - G1 * (nn - rhx * rhy) / r
    return out


def _laplace_panel(x, nx, c, t, n, a, which, self_mask):
    """Analytic integrals of the Laplace kernels over straight panels.

    x, nx: (P, 2) targets/target normals; c, t, n: (P, 2) panel centre/tangent/normal;
    a: (P,) half-lengths; self_mask: (P,) target is the panel's own midpoint
    (principal value / finite part)."""
    d = x - c
    xi = d[:, 0] * t[:, 0] + d[:, 1] * t[:, 1]
    eta = d[:, 0] * n[:, 0] + d[:, 1] * n[:, 1]
    eta = np.where(self_mask, 0.0, eta)
    xi = np.where(self_mask, 0.0, xi)
    r1s = (xi + a) ** 2 + eta ** 2   # to the s = -a end
    r2s = (xi - a) ** 2 + eta ** 2   # to the s = +a end
    out = {}
    theta = np.arctan2(2 * a * eta, xi ** 2 + eta ** 2 - a ** 2)
    theta = np.where(self_mask, 0.0, theta)  # principal value on the panel itself
    if "S" in which:
        def F(u):
            u2 = u * u + eta * eta
            lg = np.where(u2 > 0, u * np.log(np.where(u2 > 0, u2, 1.0)), 0.0)
            return lg - 2 * u
        # 2 eta atan(u/eta) terms combine to eta * theta
        val = F(a - xi) - F(-a - xi) + 2 * eta * theta
        out["S"] = -val / (4 * np.pi)
    if "K" in which:
        out["K"] = theta / TWO_PI
    if "A" in which or "H" in which:
        tx = nx[:, 0] * t[:, 0] + nx[:, 1] * t[:, 1]
        nnx = nx[:, 0] * n[:, 0] + nx[:, 1] * n[:, 1]
    if "A" in which:
        dS_dxi = -np.log(r1s / r2s) / (4 * np.pi)
        dS_deta = -theta / TWO_PI
        out["A"] = tx * dS_dxi + nnx * dS_deta
    if "H" in which:
        dD_dxi = (eta / r1s - eta / r2s) / TWO_PI
        dD_deta = -((a - xi) / r2s + (a + xi) / r1s) / TWO_PI
        out["H"] = tx * dD_dxi + nnx * dD_deta
    return out


def pair_integrals(kern, X, NX, P: Panels, ii, jj, which, selfm=None, near=1.5, mid=3.0):
    """Kernel integrals for the target/panel pairs (X[ii], panel jj): dict name -> (npairs,)."""
    ii = np.asarray(ii, int)
    jj = np.asarray(jj, int)
    if selfm is None:
        selfm = np.zeros(ii.size, bool)
    out = {w: np.empty(ii.size, complex) for w in which}
    if ii.size == 0:
        return out
    x = X[ii]
    nx = NX[ii]
    dist = np.hypot(x[:, 0] - P.c[jj, 0], x[:, 1] - P.c[jj, 1])
    ratio = dist / P.L[jj]
    for m, sel in ((2, ratio >= mid), (6, (ratio < mid) & (ratio >= near))):
        k = np.nonzero(sel)[0]
        if k.size == 0:
            continue
        gg, ww = _GL[m]
        j = jj[k]
        ys = P.c[j][:, None, :] + (0.5 * P.L[j])[:, None, None] * gg[None, :, None] * P.t[j][:, None, :]
        kv = _kernel_values(kern, x[k][:, None, :] - ys, nx[k][:, None, :], P.n[j][:, None, :], which)
        wts = (0.5 * P.L[j])[:, None] * ww[None, :]
        for w in which:
            out[w][k] = np.sum(kv[w] * wts, axis=1)
    k = np.nonzero(ratio < near)[0]
    if k.size:
        j = jj[k]
        sm = selfm[k]
        lap = _laplace_panel(x[k], nx[k], P.c[j], P.t[j], P.n[j], 0.5 * P.L[j], which, sm)
        g16, w16 = _GL[16]
        ys = P.c[j][:, None, :] + (0.5 * P.L[j])[:, None, None] * g16[None, :, None] * P.t[j][:, None, :]
        kv = _kernel_values(kern, x[k][:, None, :] - ys, nx[k][:, None, :], P.n[j][:, None, :], which, diff=True)
        wts = (0.5 * P.L[j])[:, None] * w16[None, :]
        diff = {w: np.sum(kv[w] * wts, axis=1) for w in which}
        if np.any(sm):
            # self panels: s = +-a u^2 on each half (removes the log singularity at s=0)
            si = np.nonzero(sm)[0]
            js = j[si]
            a = 0.5 * P.L[js]
            u = 0.5 * (g16 + 1.0)
            wu = 0.5 * w16
            s = a[:, None] * u[None, :] ** 2
            ds = 2 * a[:, None] * u[None, :] * wu[None, :]
            acc = {w: np.zeros(si.size, complex) for w in which}
            for sign in (1.0, -1.0):
                ys = P.c[js][:, None, :] + sign * s[:, :, None] * P.t[js][:, None, :]
                kvs = _kernel_values(kern, x[k][si][:, None, :] - ys, nx[k][si][:, None, :],
                                     P.n[js][:, None, :], which, diff=True)
                for w in which:
                    acc[w] += np.sum(kvs[w] * ds, axis=1)
            for w in which:
                diff[w][si] = acc[w]
            for w in ("K", "A"):
                if w in which:   # vanish identically on a flat panel
                    diff[w][si] = 0.0
                    lap[w][si] = 0.0
        for w in which:
            out[w][k] = lap[w] + diff[w]
    return out


def panel_integrals(kern, X, NX, P: Panels, which=("S", "K"), self_of=None,
                    near=1.5, mid=3.0, chunk_pairs=400_000):
    """Integrals over every panel of P of the kernels in ``which`` for targets X.

    Returns dict name -> (M, N) complex array with
      S = int G ds, K = int dG/dn_y ds, A = int dG/dn_x ds, H = int d2G/dn_x dn_y ds
    (finite part / principal value when target i is the midpoint of panel
    ``self_of[i]``; ``self_of[i] = -1`` for none).
    """
    X = np.asarray(X, float).reshape(-1, 2)
    M, N = X.shape[0], len(P)
    NX = np.zeros_like(X) if NX is None else np.asarray(NX, float).reshape(-1, 2)
    if self_of is None:
        self_of = -np.ones(M, int)
    self_of = np.asarray(self_of, int)
    out = {w: np.zeros((M, N), complex) for w in which}
    rows = max(1, chunk_pairs // max(1, N))
    cutoff = getattr(kern, "cutoff", np.inf)
    for i0 in range(0, M, rows):
        i1 = min(M, i0 + rows)
        ii, jj = np.meshgrid(np.arange(i0, i1), np.arange(N), indexing="ij")
        ii, jj = ii.ravel(), jj.ravel()
        if np.isfinite(cutoff):   # exponentially small interactions (K0) are dropped
            d = np.hypot(X[ii, 0] - P.c[jj, 0], X[ii, 1] - P.c[jj, 1]) - P.L[jj]
            keep = d < cutoff
            ii, jj = ii[keep], jj[keep]
        vals = pair_integrals(kern, X, NX, P, ii, jj, which, self_of[ii] == jj, near, mid)
        for w in which:
            out[w][ii, jj] = vals[w]
    return out


# ------------------------------------------------------------------- problem
@dataclass
class Domain:
    kernel: object          # Helmholtz | ModHelmholtz
    rho: float
    exterior: bool = False


@dataclass
class Problem:
    """Multi-domain BEM problem.

    panels: all panels. side_a/side_b: domain index per panel (-1 unused for sheets).
    kind: 0 = membrane, 1 = sheet (exterior both sides). Z: (N,) complex impedance.
    ties: list of (i, j, value): extra normal force on panel i = value * v_j
    (per unit area; added as p_A - p_B = Z v - F v).
    """
    omega: float
    panels: Panels
    domains: list
    side_a: np.ndarray
    side_b: np.ndarray
    kind: np.ndarray
    Z: np.ndarray
    ties: list = field(default_factory=list)
    burton_miller: bool = True
    alpha: complex | None = None


@dataclass
class Solution:
    prob: Problem
    p_ext: np.ndarray      # exterior pressure on exterior-facing membrane panels (all N; nan elsewhere)
    q_ext: np.ndarray      # exterior dp/dnu (nu outward from the exterior) on those panels
    s_ext: np.ndarray      # orientation sign of nu relative to panel n
    mu: np.ndarray         # pressure jump on sheet panels (0 elsewhere)
    pA: np.ndarray
    pB: np.ndarray
    v: np.ndarray
    cond: float = np.nan

    def field(self, X, inc_p=None, grad_dirs=None):
        """Exterior pressure at points X (M,2) (+ incident ``inc_p`` (M,) if given).
        With grad_dirs=(M,2) unit vectors, returns the directional derivative instead
        (``inc_p`` then is the incident directional derivative)."""
        prob = self.prob
        ext = [d for d, dom in enumerate(prob.domains) if dom.exterior][0]
        kern = prob.domains[ext].kernel
        P = prob.panels
        body = np.nonzero(np.isfinite(self.p_ext))[0]
        sheet = np.nonzero(prob.kind == 1)[0]
        X = np.asarray(X, float).reshape(-1, 2)
        out = np.zeros(X.shape[0], complex) if inc_p is None else np.asarray(inc_p, complex).copy()
        if grad_dirs is None:
            I = panel_integrals(kern, X, None, P, ("S", "K"))
            out += I["S"][:, body] @ self.q_ext[body] - I["K"][:, body] @ (self.s_ext[body] * self.p_ext[body])
            if sheet.size:
                out -= I["K"][:, sheet] @ self.mu[sheet]
        else:
            I = panel_integrals(kern, X, grad_dirs, P, ("A", "H"))
            out += I["A"][:, body] @ self.q_ext[body] - I["H"][:, body] @ (self.s_ext[body] * self.p_ext[body])
            if sheet.size:
                out -= I["H"][:, sheet] @ self.mu[sheet]
        return out


def solve(prob: Problem, inc_p, inc_dpdn, return_cond=False) -> Solution:
    """Solve. inc_p/inc_dpdn: incident pressure and dp/dn (along panel n) at panel
    midpoints (only used on exterior-facing panels)."""
    P = prob.panels
    N = len(P)
    w = prob.omega
    nd = len(prob.domains)
    ext = [d for d, dom in enumerate(prob.domains) if dom.exterior]
    assert len(ext) == 1
    ext = ext[0]
    mem = np.nonzero(prob.kind == 0)[0]
    sheet = np.nonzero(prob.kind == 1)[0]
    Nm, Ns = mem.size, sheet.size
    # global unknowns: [pA (Nm), wv (Nm), ws (Ns)]
    # membrane:  v = D w (D real diag scaling), pB = pA - (Z - F) v
    # sheet:     v = Z0 ws / (Z + Z0), mu = Z v
    Zm = prob.Z[mem]
    rho_ext = prob.domains[ext].rho
    k_ext = prob.domains[ext].kernel.k
    Z0 = rho_ext * w / k_ext
    Dm = Z0 / (np.abs(Zm) + Z0)
    loc = -np.ones(N, int)
    loc[mem] = np.arange(Nm)
    locs = -np.ones(N, int)
    locs[sheet] = np.arange(Ns)
    # (Z - F) restricted to membranes, as dense Nm x Nm (ties are sparse, few)
    ZF = np.diag(Zm.astype(complex))
    for (i, j, val) in prob.ties:
        ZF[loc[i], loc[j]] -= val
    ZFD = ZF * Dm[None, :]            # pB = pA - ZFD @ w
    Zs = prob.Z[sheet]
    vs_scale = Z0 / (Zs + Z0)         # v = vs_scale * ws
    mu_scale = Zs * vs_scale          # mu = mu_scale * ws

    nunk = 2 * Nm + Ns
    rows = []
    rhs = []
    alpha = prob.alpha if prob.alpha is not None else -1j / k_ext

    for d, dom in enumerate(prob.domains):
        on_a = mem[prob.side_a[mem] == d]
        on_b = mem[prob.side_b[mem] == d]
        idx = np.concatenate([on_a, on_b])
        if idx.size == 0:
            continue
        s = np.concatenate([np.ones(on_a.size), -np.ones(on_b.size)])
        isB = np.concatenate([np.zeros(on_a.size, bool), np.ones(on_b.size, bool)])
        Pd = Panels(P.a[idx], P.b[idx])
        # orientation: use panel normals as given (Panels recomputes the same normals)
        Xc = P.c[idx]
        NXc = P.n[idx] * s[:, None]           # outward normal of domain d at targets
        self_of = np.arange(idx.size)
        which = ("S", "K", "A", "H") if (dom.exterior and prob.burton_miller) else ("S", "K")
        I = panel_integrals(dom.kernel, Xc, NXc, Pd, which, self_of=self_of)
        nloc = idx.size
        # p_d(j) = pA(j) [A side]  or  pA(j) - ZFD w  [B side];  q_d(j) = s_j (-j w rho) v_j
        # CBIE:  1/2 p_d + sum K_nu p_d - sum S q_d = p_inc
        Kmat = I["K"] * s[None, :]            # K wrt nu_y = s_j n_j
        Cp = 0.5 * np.eye(nloc) + Kmat        # coefficient on p_d(j)
        Cq = -I["S"]                          # coefficient on q_d(j)
        if dom.exterior and prob.burton_miller:
            # HBIE (x-normal is nu_x = NXc already): 1/2 q - A q + H_nu p = dpinc/dnu
            Hmat = I["H"] * s[None, :]
            Cp = Cp + alpha * Hmat
            Cq = Cq + alpha * (0.5 * np.eye(nloc) - I["A"])
        R = np.zeros((nloc, nunk), complex)
        li = loc[idx]
        # p_d columns
        R[:, li] += Cp
        if np.any(isB):
            R[:, Nm:2 * Nm] -= Cp[:, isB] @ ZFD[li[isB], :]
        # q_d = s (-j w rho) v, v = D w
        R[:, Nm + li] += Cq * (s * (-1j * w * dom.rho) * Dm[li])[None, :]
        b = np.zeros(nloc, complex)
        if dom.exterior:
            b += inc_p[idx]
            if prob.burton_miller:
                b += alpha * inc_dpdn[idx] * s
            if Ns:
                Psh = Panels(P.a[sheet], P.b[sheet])
                wh = ("K", "H") if prob.burton_miller else ("K",)
                Ish = panel_integrals(dom.kernel, Xc, NXc, Psh, wh)
                Csh = Ish["K"]
                if prob.burton_miller:
                    Csh = Csh + alpha * Ish["H"]
                R[:, 2 * Nm:] += Csh * mu_scale[None, :]
        rows.append(R)
        rhs.append(b)

    if Ns:
        # sheet rows: -j w rho v - sum_body [A q - H_nu p] + sum_sh H mu = dpinc/dn
        dom = prob.domains[ext]
        on_a = mem[prob.side_a[mem] == ext]
        on_b = mem[prob.side_b[mem] == ext]
        idx = np.concatenate([on_a, on_b])
        s = np.concatenate([np.ones(on_a.size), -np.ones(on_b.size)])
        isB = np.concatenate([np.zeros(on_a.size, bool), np.ones(on_b.size, bool)])
        Xs, NXs = P.c[sheet], P.n[sheet]
        Pd = Panels(P.a[idx], P.b[idx])
        I = panel_integrals(dom.kernel, Xs, NXs, Pd, ("A", "H"))
        Psh = Panels(P.a[sheet], P.b[sheet])
        Ish = panel_integrals(dom.kernel, Xs, NXs, Psh, ("H",), self_of=np.arange(Ns))
        R = np.zeros((Ns, nunk), complex)
        li = loc[idx]
        Cp = I["H"] * s[None, :]
        R[:, li] += Cp
        if np.any(isB):
            R[:, Nm:2 * Nm] -= Cp[:, isB] @ ZFD[li[isB], :]
        R[:, Nm + li] += -I["A"] * (s * (-1j * w * dom.rho) * Dm[li])[None, :]
        R[:, 2 * Nm:] += Ish["H"] * mu_scale[None, :]
        R[np.arange(Ns), 2 * Nm + np.arange(Ns)] += -1j * w * dom.rho * vs_scale
        rows.append(R)
        rhs.append(np.asarray(inc_dpdn)[sheet].astype(complex))

    A = np.vstack(rows)
    b = np.concatenate(rhs)
    assert A.shape == (nunk, nunk), A.shape
    x = np.linalg.solve(A, b)
    cond = np.linalg.cond(A) if return_cond else np.nan
    pAm = x[:Nm]
    wv = x[Nm:2 * Nm]
    ws = x[2 * Nm:]
    pA = np.full(N, np.nan + 0j)
    pB = np.full(N, np.nan + 0j)
    v = np.zeros(N, complex)
    mu = np.zeros(N, complex)
    pA[mem] = pAm
    pB[mem] = pAm - ZFD @ wv
    v[mem] = Dm * wv
    v[sheet] = vs_scale * ws
    mu[sheet] = mu_scale * ws
    # exterior traces (nu = outward normal of the exterior domain = s * n)
    p_ext = np.full(N, np.nan + 0j)
    q_ext = np.zeros(N, complex)
    s_ext = np.zeros(N)
    on_a = mem[prob.side_a[mem] == ext]
    on_b = mem[prob.side_b[mem] == ext]
    p_ext[on_a] = pA[on_a]
    p_ext[on_b] = pB[on_b]
    s_ext[on_a] = 1.0
    s_ext[on_b] = -1.0
    q_ext[on_a] = -1j * w * rho_ext * v[on_a]
    q_ext[on_b] = 1j * w * rho_ext * v[on_b]
    return Solution(prob, p_ext, q_ext, s_ext, mu, pA, pB, v, cond)


# ------------------------------------------------------------------- fast (structured) solver
#
# For long flat bodies (slab shields) the face panels form uniform collinear runs
# ("groups": equal length, equal direction, consecutive). Kernel blocks between two such
# runs are Toeplitz, so they are generated from O(n) kernel integrals and applied with
# FFTs. Remaining panels ("irregular": sides, partitions) are handled with dense blocks.
# The system is solved with GMRES, preconditioned by a sparse LU of the near-field part
# (all pairs closer than ``r_near``) of exactly the same operator.

import scipy.sparse as sps
from scipy.fft import fft, ifft, next_fast_len
from scipy.sparse.linalg import LinearOperator, gmres, splu


class _ToepBlock:
    def __init__(self, kern, P, g1, g2, which):
        n1, n2 = len(g1), len(g2)
        h = P.L[g1[0]]
        t = P.t[g1[0]]
        if not (np.allclose(P.L[g2], h, rtol=1e-9) and np.allclose(P.t[g2], t, atol=1e-12)):
            raise ValueError("groups must share spacing and direction")
        m = np.arange(-(n2 - 1), n1)
        X = P.c[g1[0]][None, :] + (m * h)[:, None] * t[None, :]   # target i=m relative to source j=0
        NX = np.repeat(P.n[g1[0]][None, :], m.size, 0)
        selfm = (g1[0] == g2[0]) & (m == 0)
        vals = pair_integrals(kern, X, NX, P, np.arange(m.size), np.full(m.size, g2[0]), which, selfm)
        self.n1, self.n2, self.m = n1, n2, m
        self.vals = vals
        self.nfft = next_fast_len(m.size + n2 - 1)
        self.F = {w: fft(vals[w], self.nfft) for w in which}
        self.g1, self.g2 = np.asarray(g1), np.asarray(g2)
        self.dist = np.hypot(*(X - P.c[g2[0]][None, :]).T)

    def apply(self, w, xf):
        """xf = fft(x[g2], nfft)."""
        return ifft(self.F[w] * xf)[self.n2 - 1:self.n2 - 1 + self.n1]

    def near(self, w, r_near):
        sel = np.nonzero(self.dist < r_near)[0]
        rows, cols, vals = [], [], []
        for k in sel:
            m = self.m[k]
            i = np.arange(max(0, m), min(self.n1, self.n2 + m))
            rows.append(self.g1[i])
            cols.append(self.g2[i - m])
            vals.append(np.full(i.size, self.vals[w][k]))
        if not rows:
            return np.zeros(0, int), np.zeros(0, int), np.zeros(0, complex)
        return np.concatenate(rows), np.concatenate(cols), np.concatenate(vals)


class FastOp:
    """Kernel integrals restricted to one domain: y = I_w (masked) @ x for w in which."""

    def __init__(self, kern, P, groups, in_dom, which, cache, r_near):
        self.N = len(P)
        self.which = which
        self.in_dom = in_dom
        gset = [np.asarray(g) for g in groups]
        ingroup = np.zeros(self.N, bool)
        for g in gset:
            ingroup[g] = True
        self.irr = np.nonzero(in_dom & ~ingroup)[0]
        self.dom = np.nonzero(in_dom)[0]
        self.blocks = []
        for a, g1 in enumerate(gset):
            if not np.any(in_dom[g1]):
                continue
            for b, g2 in enumerate(gset):
                if not np.any(in_dom[g2]):
                    continue
                key = (id(kern), a, b, which)
                if key not in cache:
                    cache[key] = _ToepBlock(kern, P, g1, g2, which)
                blk = cache[key]
                self.blocks.append((blk, in_dom[g1], in_dom[g2]))
        # dense blocks with irregular panels
        grp_dom = np.setdiff1d(self.dom, self.irr)
        self.grp_dom = grp_dom
        if self.irr.size:
            self.B1 = panel_integrals(kern, P.c[self.irr], P.n[self.irr], _sub(P, self.dom), which,
                                      self_of=np.searchsorted(self.dom, self.irr))
            self.B2 = panel_integrals(kern, P.c[grp_dom], P.n[grp_dom], _sub(P, self.irr), which)
        self.P = P
        self.r_near = r_near
        self._mat = {}

    def matvec(self, w, x):
        y = np.zeros(self.N, complex)
        xm = np.where(self.in_dom, x, 0.0)
        if getattr(self, "_fx", None) is not None and self._fx[0] is x:
            fcache = self._fx[1]
        else:
            fcache = {}
            self._fx = (x, fcache)
        for blk, m1, m2 in self.blocks:
            key = (id(blk.g2), blk.nfft)
            if key not in fcache:
                fcache[key] = fft(xm[blk.g2], blk.nfft)
            y[blk.g1] += np.where(m1, blk.apply(w, fcache[key]), 0.0)
        if self.irr.size:
            y[self.irr] += self.B1[w] @ xm[self.dom]
            y[self.grp_dom] += self.B2[w] @ xm[self.irr]
        return y

    def near_matrix(self, w):
        if w in self._mat:
            return self._mat[w]
        R, C, V = [], [], []
        for blk, m1, m2 in self.blocks:
            r, c, v = blk.near(w, self.r_near)
            ok = self.in_dom[r] & self.in_dom[c]
            R.append(r[ok]); C.append(c[ok]); V.append(v[ok])
        P = self.P
        if self.irr.size:
            for rows, cols, B in ((self.irr, self.dom, self.B1[w]), (self.grp_dom, self.irr, self.B2[w])):
                d = np.hypot(P.c[rows][:, None, 0] - P.c[cols][None, :, 0], P.c[rows][:, None, 1] - P.c[cols][None, :, 1])
                ii, jj = np.nonzero(d < self.r_near)
                R.append(rows[ii]); C.append(cols[jj]); V.append(B[ii, jj])
        M = sps.csr_matrix((np.concatenate(V), (np.concatenate(R), np.concatenate(C))), shape=(self.N, self.N))
        self._mat[w] = M
        return M


def _sub(P, idx):
    return Panels(P.a[idx], P.b[idx])


class _Lin:
    """Wraps 'apply kernel w of op' so the same algebra works on vectors and sparse matrices."""

    def __init__(self, op, w):
        self.op, self.w = op, w

    def __matmul__(self, x):
        if sps.issparse(x):
            return self.op.near_matrix(self.w) @ x
        return self.op.matvec(self.w, x)


def solve_fast(prob: Problem, inc_p, inc_dpdn, groups, tol=1e-7, r_near_factor=3.0, maxiter=2000,
               verbose=False, r_near=None, precond="circulant") -> Solution:
    """Structured solve (same equations as :func:`solve`). ``groups``: list of index arrays
    of uniform collinear panel runs (all with the same spacing and direction)."""
    P = prob.panels
    N = len(P)
    w = prob.omega
    ext = [d for d, dom in enumerate(prob.domains) if dom.exterior][0]
    mem = np.nonzero(prob.kind == 0)[0]
    sheet = np.nonzero(prob.kind == 1)[0]
    Nm, Ns = mem.size, sheet.size
    nunk = 2 * Nm + Ns
    rho_ext = prob.domains[ext].rho
    k_ext = prob.domains[ext].kernel.k
    Z0 = rho_ext * w / k_ext
    Zm = np.zeros(N, complex)
    Zm[mem] = prob.Z[mem]
    Dm = np.zeros(N)
    Dm[mem] = Z0 / (np.abs(prob.Z[mem]) + Z0)
    vs = np.zeros(N, complex)
    mus = np.zeros(N, complex)
    vs[sheet] = Z0 / (prob.Z[sheet] + Z0)
    mus[sheet] = prob.Z[sheet] * vs[sheet]
    ZF = sps.diags(Zm).tolil()
    for (i, j, val) in prob.ties:
        ZF[i, j] -= val
    ZF = ZF.tocsr()
    alpha = prob.alpha if prob.alpha is not None else -1j / k_ext
    h = float(np.median(P.L))
    import os as _os
    if r_near is None:
        r_near = r_near_factor * h
    r_near = float(_os.environ.get("BEM_RNEAR", r_near))
    # maps from unknowns to full-N vectors (sparse N x nunk)
    E_pA = sps.csr_matrix((np.ones(Nm), (mem, np.arange(Nm))), shape=(N, nunk))
    E_v = sps.csr_matrix((Dm[mem], (mem, Nm + np.arange(Nm))), shape=(N, nunk))
    if Ns:
        E_v = E_v + sps.csr_matrix((vs[sheet], (sheet, 2 * Nm + np.arange(Ns))), shape=(N, nunk))
        E_mu = sps.csr_matrix((mus[sheet], (sheet, 2 * Nm + np.arange(Ns))), shape=(N, nunk))
    else:
        E_mu = sps.csr_matrix((N, nunk))
    cache = {}
    doms = []
    for d, dom in enumerate(prob.domains):
        onA = np.zeros(N, bool)
        onB = np.zeros(N, bool)
        onA[mem] = prob.side_a[mem] == d
        onB[mem] = prob.side_b[mem] == d
        in_dom = onA | onB
        if dom.exterior:
            in_dom = in_dom.copy()
            in_dom[sheet] = True
        if not np.any(onA | onB):
            continue
        which = ("S", "K", "A", "H") if dom.exterior else ("S", "K")
        op = FastOp(dom.kernel, P, groups, in_dom, which, cache, r_near)
        doms.append((d, dom, onA, onB, op))
    D = sps.diags

    def rows(x_pA, x_v, x_mu):
        """Residual rows given full-N pA, v, mu (vectors or sparse N x nunk)."""
        pB = x_pA - ZF @ x_v
        out = []
        for d, dom, onA, onB, op in doms:
            sd = onA.astype(float) - onB.astype(float)
            pd = D(onA.astype(float)) @ x_pA + D(onB.astype(float)) @ pB
            qd = D(sd * (-1j * w * dom.rho)) @ x_v
            xK = D(sd) @ pd
            if dom.exterior:
                xK = xK + x_mu
            r = 0.5 * pd + _Lin(op, "K") @ xK - _Lin(op, "S") @ qd
            if dom.exterior:
                Hx = _Lin(op, "H") @ xK
                Ax = _Lin(op, "A") @ qd
                r = r + alpha * (D(sd) @ Hx + 0.5 * qd - D(sd) @ Ax)
            out.append(_take(r, np.nonzero(onA | onB)[0]))
            if dom.exterior and Ns:
                rs = Hx - Ax + D(np.where(np.isin(np.arange(N), sheet), -1j * w * dom.rho, 0.0)) @ x_v
                out.append(_take(rs, sheet))
        return out

    def matvec(u):
        u = np.asarray(u).ravel()
        res = rows(E_pA @ u, E_v @ u, E_mu @ u)
        return np.concatenate(res)

    Mnear = sps.vstack(rows(E_pA, E_v, E_mu)).tocsr()
    # row bookkeeping: (panel, domain) of every residual row
    row_panel, row_dom = [], []
    for d, dom, onA, onB, op in doms:
        idx = np.nonzero(onA | onB)[0]
        row_panel.append(idx)
        row_dom.append(np.full(idx.size, d))
        if dom.exterior and Ns:
            row_panel.append(sheet)
            row_dom.append(np.full(Ns, -2))
    row_panel = np.concatenate(row_panel)
    row_dom = np.concatenate(row_dom)
    # rhs
    b = []
    for d, dom, onA, onB, op in doms:
        idx = np.nonzero(onA | onB)[0]
        bb = np.zeros(idx.size, complex)
        if dom.exterior:
            sd = (onA.astype(float) - onB.astype(float))[idx]
            bb = inc_p[idx] + alpha * inc_dpdn[idx] * sd
        b.append(bb)
        if dom.exterior and Ns:
            b.append(np.asarray(inc_dpdn)[sheet].astype(complex))
    b = np.concatenate(b)
    assert Mnear.shape == (nunk, nunk), Mnear.shape
    Aop = LinearOperator((nunk, nunk), matvec=matvec, dtype=complex)
    pre = None
    if precond == "circulant" and len(groups) >= 2 and len(groups[0]) == len(groups[1]):
        pre = _circulant_precond(matvec, Mnear, groups[0], groups[1], prob, ext, mem, loc_of(mem, N),
                                 row_panel, row_dom, nunk, Nm)
    if pre is None:
        lu = splu(Mnear.tocsc())
        pre = lu.solve
    Mop = LinearOperator((nunk, nunk), matvec=pre, dtype=complex)
    it = [0]

    def cb(_):
        it[0] += 1
    x, info = gmres(Aop, b, M=Mop, rtol=tol, restart=300, maxiter=maxiter, callback=cb,
                    callback_type="pr_norm")
    if info != 0:
        raise RuntimeError(f"GMRES did not converge (info={info})")
    if verbose:
        print(f"  GMRES iterations: {it[0]}, N={N}")
    pAm = x[:Nm]
    wv = x[Nm:2 * Nm]
    ws = x[2 * Nm:]
    full_pA = E_pA @ x
    full_v = E_v @ x
    pA = np.full(N, np.nan + 0j)
    pB = np.full(N, np.nan + 0j)
    pA[mem] = pAm
    pB[mem] = (full_pA - ZF @ full_v)[mem]
    v = full_v
    mu = E_mu @ x
    p_ext = np.full(N, np.nan + 0j)
    q_ext = np.zeros(N, complex)
    s_ext = np.zeros(N)
    on_a = mem[prob.side_a[mem] == ext]
    on_b = mem[prob.side_b[mem] == ext]
    p_ext[on_a] = pA[on_a]
    p_ext[on_b] = pB[on_b]
    s_ext[on_a] = 1.0
    s_ext[on_b] = -1.0
    q_ext[on_a] = -1j * w * rho_ext * v[on_a]
    q_ext[on_b] = 1j * w * rho_ext * v[on_b]
    sol = Solution(prob, p_ext, q_ext, s_ext, mu, pA, pB, v)
    sol.iterations = it[0]
    return sol


def _take(r, idx):
    if sps.issparse(r):
        return r.tocsr()[idx]
    return r[idx]


def loc_of(mem, N):
    loc = -np.ones(N, int)
    loc[mem] = np.arange(mem.size)
    return loc


def _circulant_precond(matvec, Mnear, gtop, gbot, prob, ext, mem, loc, row_panel, row_dom, nunk, Nm):
    """Block-circulant (Strang) preconditioner for two aligned face groups (top/bottom,
    same y nodes): per face index i the 4 unknowns (pA, w) x (top, bottom) and the 4 rows
    (exterior, interior) x (top, bottom). The face-to-face operator is block Toeplitz;
    its generating sequences are probed with 4 matvecs at the middle index, wrapped into
    circulants and inverted per Fourier mode (4x4). This is the exact inverse for an
    infinite slab, so it captures the global top/bottom coupling. Remaining unknowns
    (sides, partitions, sheets) use a sparse LU of their near-field diagonal block."""
    n = len(gtop)
    if np.any(prob.kind[gtop] != 0) or np.any(prob.kind[gbot] != 0):
        return None
    # rows: for each group panel, its exterior row and its interior row
    key = {(p, d): r for r, (p, d) in enumerate(zip(row_panel, row_dom))}
    R = np.zeros((n, 4), int)
    U = np.zeros((n, 4), int)
    for gi, g in enumerate((gtop, gbot)):
        for i, p in enumerate(g):
            other = prob.side_b[p] if prob.side_a[p] == ext else prob.side_a[p]
            if (p, ext) not in key or (p, other) not in key:
                return None
            R[i, 2 * gi] = key[(p, ext)]
            R[i, 2 * gi + 1] = key[(p, other)]
            U[i, 2 * gi] = loc[p]
            U[i, 2 * gi + 1] = Nm + loc[p]
    ic = n // 2
    c = np.zeros((n, 4, 4), complex)
    shift = (np.arange(n) - ic) % n
    for u in range(4):
        e = np.zeros(nunk, complex)
        e[U[ic, u]] = 1.0
        y = matvec(e)
        for r in range(4):
            c[shift, r, u] = y[R[:, r]]
    Chat = np.fft.fft(c, axis=0)
    try:
        Cinv = np.linalg.inv(Chat)
    except np.linalg.LinAlgError:
        return None
    grp_rows = R.ravel()
    grp_unk = U.ravel()
    rest_rows = np.setdiff1d(np.arange(nunk), grp_rows)
    rest_unk = np.setdiff1d(np.arange(nunk), grp_unk)
    lu_rest = None
    if rest_rows.size:
        if rest_rows.size != rest_unk.size:
            return None
        sub = Mnear[rest_rows][:, rest_unk].tocsc()
        lu_rest = splu(sub)

    def apply(r):
        r = np.asarray(r).ravel()
        x = np.zeros(nunk, complex)
        rh = np.fft.fft(r[R], axis=0)                     # (n, 4)
        xh = np.einsum("kij,kj->ki", Cinv, rh)
        x[U] = np.fft.ifft(xh, axis=0)
        if lu_rest is not None:
            x[rest_unk] = lu_rest.solve(r[rest_rows])
        return x

    return apply
