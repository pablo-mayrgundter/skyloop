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
        j0, y0, j1, y1 = sp.j0(z), sp.y0(z), sp.j1(z), sp.y1(z)
        H0 = j0 - 1j * y0
        H1 = j1 - 1j * y1
        G = -0.25j * H0
        G1 = 0.25j * k * H1
        G2 = 0.25j * k * k * (H0 - H1 / z)
        return G, G1, G2


class ModHelmholtz:
    """G = K0(kappa r) / (2 pi)."""

    def __init__(self, kappa: float):
        self.kappa = max(float(kappa), 1e-12)

    def radial(self, r):
        a = self.kappa
        z = a * r
        K0, K1 = sp.k0(z), sp.k1(z)
        G = K0 / TWO_PI
        G1 = -a * K1 / TWO_PI
        G2 = a * a * (K0 + K1 / z) / TWO_PI
        return G + 0j, G1 + 0j, G2 + 0j


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


def panel_integrals(kern, X, NX, P: Panels, which=("S", "K"), self_of=None,
                    near=1.5, mid=3.0, chunk_elems=3_000_000):
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
    out = {w: np.empty((M, N), complex) for w in which}
    # --- regular quadrature (2 or 6 points by distance) in row chunks
    rows = max(1, chunk_elems // max(1, N * 6))
    g2, w2 = _GL[2]
    g6, w6 = _GL[6]
    near_pairs = []
    for i0 in range(0, M, rows):
        i1 = min(M, i0 + rows)
        x = X[i0:i1]
        dist = np.hypot(x[:, None, 0] - P.c[None, :, 0], x[:, None, 1] - P.c[None, :, 1])
        ratio = dist / P.L[None, :]
        for gg, ww, sel in ((g2, w2, ratio >= mid), (g6, w6, (ratio < mid) & (ratio >= near))):
            ii, jj = np.nonzero(sel)
            if ii.size == 0:
                continue
            ys = P.c[jj][:, None, :] + (0.5 * P.L[jj])[:, None, None] * gg[None, :, None] * P.t[jj][:, None, :]
            rvec = x[ii][:, None, :] - ys
            kv = _kernel_values(kern, rvec, NX[i0 + ii][:, None, :], P.n[jj][:, None, :], which)
            wts = (0.5 * P.L[jj])[:, None] * ww[None, :]
            for w in which:
                out[w][i0 + ii, jj] = np.sum(kv[w] * wts, axis=1)
        ii, jj = np.nonzero(ratio < near)
        near_pairs.append((i0 + ii, jj))
    ii = np.concatenate([p[0] for p in near_pairs]) if near_pairs else np.zeros(0, int)
    jj = np.concatenate([p[1] for p in near_pairs]) if near_pairs else np.zeros(0, int)
    if ii.size:
        selfm = self_of[ii] == jj
        lap = _laplace_panel(X[ii], NX[ii], P.c[jj], P.t[jj], P.n[jj], 0.5 * P.L[jj], which, selfm)
        # difference kernel, 16-point Gauss (non-self)
        g16, w16 = _GL[16]
        ys = P.c[jj][:, None, :] + (0.5 * P.L[jj])[:, None, None] * g16[None, :, None] * P.t[jj][:, None, :]
        rvec = X[ii][:, None, :] - ys
        kv = _kernel_values(kern, rvec, NX[ii][:, None, :], P.n[jj][:, None, :], which, diff=True)
        wts = (0.5 * P.L[jj])[:, None] * w16[None, :]
        diff = {w: np.sum(kv[w] * wts, axis=1) for w in which}
        if np.any(selfm):
            # self panels: s = +-a u^2 on each half (removes the log singularity at s=0)
            si = np.nonzero(selfm)[0]
            a = 0.5 * P.L[jj[si]]
            u = 0.5 * (g16 + 1.0)
            wu = 0.5 * w16
            s = a[:, None] * u[None, :] ** 2
            ds = 2 * a[:, None] * u[None, :] * wu[None, :]
            for sign in (1.0, -1.0):
                ys = P.c[jj[si]][:, None, :] + sign * s[:, :, None] * P.t[jj[si]][:, None, :]
                rvec = X[ii[si]][:, None, :] - ys
                kvs = _kernel_values(kern, rvec, NX[ii[si]][:, None, :], P.n[jj[si]][:, None, :],
                                     which, diff=True)
                for w in which:
                    if sign > 0:
                        diff[w][si] = 0.0
                    diff[w][si] += np.sum(kvs[w] * ds, axis=1)
            for w in ("K", "A"):
                if w in which:   # vanish identically on a flat panel
                    diff[w][si] = 0.0
                    lap[w][si] = 0.0
        for w in which:
            out[w][ii, jj] = lap[w] + diff[w]
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
