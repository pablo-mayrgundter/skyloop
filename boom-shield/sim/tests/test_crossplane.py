"""Tests for the cross-plane finite-width shield solver (boomsim.crossplane, boomsim.bem2d).

Fast tests run by default; ``@pytest.mark.slow`` ones (a few seconds to ~1 min each)
are skipped unless ``--runslow`` is given or ``RUN_SLOW=1`` is set.
"""
import os
import sys

import numpy as np
import pytest
from scipy import special as sp

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))

from boomsim import bem2d  # noqa: E402
from boomsim import crossplane as cp  # noqa: E402
from boomsim import shield as SH  # noqa: E402
from boomsim.source import CONCORDE, near_field  # noqa: E402

slow = pytest.mark.skipif(not os.environ.get("RUN_SLOW"), reason="slow; set RUN_SLOW=1")


def _lin(**kw):
    kw.setdefault("incident", "linear")
    return cp.CrossPlaneConfig(CONCORDE, **kw)


# ----------------------------------------------------------------------------- incident
def _whitham_vs_exact(r, phi=0.0):
    return cp.incident_vs_whitham(CONCORDE, r, phi=phi)


@pytest.mark.parametrize("r_over_L,tol_peak,tol_max", [(1.5, 0.025, 0.05), (5.0, 0.008, 0.015)])
def test_incident_matches_whitham(r_over_L, tol_peak, tol_max):
    """Exact Hankel field at phi = 0 vs Whitham near field, both band-limited to 100 Hz
    (measured: 1.5 L peak +1.5 %, max error 3.5 %; 5 L +0.4 %, 1.0 %; the difference is
    the O(1/(k r)) near-field correction that Whitham's asymptote omits)."""
    t, p, pw = _whitham_vs_exact(r_over_L * CONCORDE.length)
    err = np.max(np.abs(p - pw)) / np.max(np.abs(pw))
    assert abs(p.max() / pw.max() - 1) < tol_peak
    assert err < tol_max, err


def test_incident_azimuth_lift_cos_phi():
    """Off-track: exact field matches near_field(phi) (A_v + A_L cos phi) at r = 5 L."""
    t, p, pw = _whitham_vs_exact(5 * CONCORDE.length, phi=np.radians(50.0))
    assert np.max(np.abs(p - pw)) / np.max(np.abs(pw)) < 0.06
    # lift part alone scales as cos(phi)
    k = 0.2
    a = [0.0, 1.0 + 0.5j]
    for ph in (0.3, 0.9):
        y, z = 300 * np.sin(ph), -300 * np.cos(ph)
        assert cp.multipole_field(a, k, y, z) == pytest.approx(np.cos(ph) * cp.multipole_field(a, k, 0.0, -300.0))


def test_multipole_plane_wave_spectrum_and_gradient():
    k = 0.3
    y = np.array([37.0, -80.0, 5.0])
    z = np.array([-120.0, -60.0, -200.0])
    th, w = sp.roots_legendre(3000)
    th, w = th * np.pi / 2, w * np.pi / 2
    u, wu = np.polynomial.legendre.leggauss(300)
    u, wu = (u + 1) * 3, wu * 3
    for n in range(4):
        a = np.zeros(n + 1)
        a[n] = 1.0
        ex, gy, gz = cp.multipole_field(a, k, y, z, grad=True)
        num = (1 / np.pi) * (np.exp(-1j * k * (np.outer(y, np.sin(th)) + np.outer(abs(z), np.cos(th))))
                             * cp._spectral_weight(a, th)) @ w
        for sgn in (1, -1):
            num += (1 / np.pi) * (np.exp(-1j * np.outer(y, sgn * k * np.cosh(u)) - np.outer(abs(z), k * np.sinh(u)))
                                  * cp._spectral_weight(a, np.pi / 2 + 1j * u) * 1j) @ wu
        np.testing.assert_allclose(num, ex, rtol=1e-8)
        e = 1e-4
        fy = (cp.multipole_field(a, k, y + e, z) - cp.multipole_field(a, k, y - e, z)) / (2 * e)
        fz = (cp.multipole_field(a, k, y, z + e) - cp.multipole_field(a, k, y, z - e)) / (2 * e)
        np.testing.assert_allclose(gy, fy, rtol=1e-6, atol=1e-9 * np.abs(ex).max())
        np.testing.assert_allclose(gz, fz, rtol=1e-6, atol=1e-9 * np.abs(ex).max())


# ----------------------------------------------------------------------------- BEM core
def _cyl_series(k, ki, modified, rho_e, rho_i, Z, a, xs, X, w, nmax=80):
    rs = np.hypot(*xs)
    ths = np.arctan2(xs[1], xs[0])
    r = np.hypot(X[:, 0], X[:, 1])
    th = np.arctan2(X[:, 1], X[:, 0])
    p = sp.hankel2(0, k * np.hypot(X[:, 0] - xs[0], X[:, 1] - xs[1]))
    for n in range(-nmax, nmax + 1):
        al = sp.hankel2(n, k * rs)
        J, Jp = sp.jv(n, k * a), sp.jvp(n, k * a)
        H, Hp = sp.hankel2(n, k * a), sp.h2vp(n, k * a)
        Zi, Zip = (sp.iv(n, ki * a), sp.ivp(n, ki * a)) if modified else (sp.jv(n, ki * a), sp.jvp(n, ki * a))
        M = np.array([[k * Hp / rho_e, -ki * Zip / rho_i], [-H + Z / (1j * w * rho_e) * k * Hp, Zi]], complex)
        rhs = np.array([-al * k * Jp / rho_e, al * J - Z / (1j * w * rho_e) * al * k * Jp])
        b, _ = np.linalg.solve(M, rhs)
        p = p + b * sp.hankel2(n, k * r) * np.exp(1j * n * (th - ths))
    return p


@pytest.mark.parametrize("ki,modified,rho_i,mass", [(0.15, True, 0.012, 1.0), (0.9, False, 0.012, 0.5),
                                                    (0.3, False, 0.088, 1e9)])
def test_bem_penetrable_cylinder_vs_series(ki, modified, rho_i, mass):
    """Membrane cylinder (K0 gas, Helmholtz gas, rigid) vs exact cylinder-harmonic series."""
    w, rho_e, k, a, N = 2 * np.pi * 10.0, 0.088, 0.3, 30.0, 400
    xs = np.array([0.0, 120.0])
    X = np.array([[0, -60.0], [40, -40], [-50, 10], [0, 80], [100, -20]])
    Z = 1j * w * mass
    th = np.linspace(0, 2 * np.pi, N + 1)
    pts = a * np.stack([np.cos(th), np.sin(th)], 1)
    P = bem2d.Panels(pts[:-1], pts[1:])
    inner = bem2d.ModHelmholtz(ki) if modified else bem2d.Helmholtz(ki)
    prob = bem2d.Problem(w, P, [bem2d.Domain(bem2d.Helmholtz(k), rho_e, True), bem2d.Domain(inner, rho_i)],
                         np.ones(N, int), np.zeros(N, int), np.zeros(N, int), np.full(N, Z))
    d = P.c - xs
    r = np.hypot(d[:, 0], d[:, 1])
    sol = bem2d.solve(prob, sp.hankel2(0, k * r), -k * sp.hankel2(1, k * r) * (d * P.n).sum(1) / r)
    dX = X - xs
    pb = sol.field(X, inc_p=sp.hankel2(0, k * np.hypot(dX[:, 0], dX[:, 1])))
    ps = _cyl_series(k, ki, modified, rho_e, rho_i, Z, a, xs, X, w)
    assert np.max(np.abs(pb - ps) / np.abs(ps)) < 0.012


def test_fast_structured_solver_matches_dense():
    cfg = _lin(h=150.0, W=300.0, T=20.0, partitions=[-50.0, 60.0], tie_K=50.0, fringe_width=40.0)
    f = 8.0
    m = cp.build_mesh(cfg, f, uniform=True)
    md = dict(m)
    md["uniform"] = False
    X = cfg.circle_points()
    p1 = cp.bem_field(cp.bem_solve(cfg, f, mesh=m), X)
    p2 = cp.bem_field(cp.bem_solve(cfg, f, mesh=md), X)
    np.testing.assert_allclose(p1, p2, rtol=2e-3, atol=2e-3 * np.abs(p2).max())


def test_no_fictitious_resonances():
    """Fine frequency scan of a small He slab (evanescent interior: no physical cavity
    modes) through the exterior BIE's fictitious frequencies (interior Dirichlet modes of
    a 60 x 20 m air box: (4,1) 14.2, (5,1) 16.5, (1,2) 17.3 Hz): the response must stay
    smooth (Burton-Miller)."""
    cfg = _lin(h=100.0, W=60.0, T=20.0, m_top=0.3, m_bot=0.3, m_side=0.3)
    X = np.array([[0.0, -200.0], [60.0, -180.0]])
    fs = np.linspace(13.5, 18.0, 19)
    vals = np.array([cp.field_at(cfg, f, X, method="bem") for f in fs])
    inc = np.array([cp.field_at(_lin(h=100.0, W=None), f, X) for f in fs])
    db = 20 * np.log10(np.abs(vals / inc))
    assert np.all(np.isfinite(db))
    assert np.abs(np.diff(db, 2, axis=0)).max() < 1.0


# ----------------------------------------------------------------------------- validation a)
def test_unshielded_matching_signatures_equal_incident_and_whitham():
    cfg = _lin(h=250.0, W=None, f_max=100.0, df=0.25, R_m=5 * CONCORDE.length, phis=[0.0, 30.0],
               atm=CONCORDE.atm, workers=1)
    sig, (f, P) = cp.matching_signatures(cfg, return_spectra=True)
    A = cp.incident_coefficients(cfg, f)
    k = 2 * np.pi * f * np.sqrt(1 / cfg.c ** 2 - 1 / cfg.V ** 2)
    X = cfg.circle_points()
    for i in (5, 50, 300):
        np.testing.assert_allclose(P[i], cp.multipole_field(A[i], k[i], X[:, 0], X[:, 1]), rtol=1e-12)
    t, p = sig[0.0]
    tw, pw = near_field(CONCORDE, cfg.matching_radius(), 20000.0, y_max_factor=12.0)
    pwi = np.interp(t, tw, pw)
    # the exact field carries the static lift footprint (f=0) that Whitham lacks; compare peaks
    assert p.max() == pytest.approx(pw.max(), rel=0.08)   # band-limited (100 Hz) vs. exact shock


# ----------------------------------------------------------------------------- validation c)
def _flux_through_rectangle(cfg, f, y0, y1, z0, z1, n=600):
    sol = cp.bem_solve(cfg, f)
    s = np.linspace(0, 1, n, endpoint=False) + 0.5 / n
    segs = [((y0, z1), (y1, z1), (0, 1)), ((y1, z1), (y1, z0), (1, 0)),
            ((y1, z0), (y0, z0), (0, -1)), ((y0, z0), (y0, z1), (-1, 0))]
    w = 2 * np.pi * f
    tot = inc_top = 0.0
    for (a, b, nrm) in segs:
        a, b = np.array(a, float), np.array(b, float)
        X = a + s[:, None] * (b - a)
        ds = np.linalg.norm(b - a) / n
        dirs = np.repeat(np.array([nrm], float), n, 0)
        p = cp.bem_field(sol, X)
        dp = cp.bem_grad(sol, X, dirs)
        vn = -dp / (1j * w * cfg.rho)
        tot += np.sum(0.5 * np.real(p * np.conj(vn))) * ds
        if nrm == (0, 1):
            pi, gy, gz = sol.inc(X[:, 0], X[:, 1])
            inc_top += np.sum(0.5 * np.real(pi * np.conj(-gz / (1j * w * cfg.rho)))) * ds
    return tot, abs(inc_top)


def test_energy_conservation_lossless_and_lossy():
    """Net outward power of the total field through a closed contour around the shield:
    zero for a lossless stack (incident = reflected + transmitted), negative (absorbed)
    with a flow-resistive face."""
    kw = dict(h=150.0, W=200.0, T=20.0, partitions=[30.0], tie_K=80.0)
    cfg = _lin(**kw)
    tot, inc = _flux_through_rectangle(cfg, 6.0, -130.0, 130.0, -200.0, -120.0)
    assert abs(tot) / inc < 5e-3
    lossy = _lin(m_top=SH.Screen(30.0, 1.0), **kw)
    tot2, inc2 = _flux_through_rectangle(lossy, 6.0, -130.0, 130.0, -200.0, -120.0)
    assert tot2 / inc2 < -0.05


# ----------------------------------------------------------------------------- symmetry / edges
def test_mirror_symmetry_asymmetric_edges_and_fringe():
    a = _lin(h=150.0, y_left=-80.0, y_right=170.0, T=20.0, fringe_width=50.0, phis=[-40.0, 0.0, 25.0])
    b = _lin(h=150.0, y_left=-170.0, y_right=80.0, T=20.0, fringe_width=50.0, phis=[40.0, 0.0, -25.0])
    for f, meth in ((6.0, "bem"), (30.0, "po")):
        pa = cp.field_at(a, f, a.circle_points(), method=meth)
        pb = cp.field_at(b, f, b.circle_points(), method=meth)
        np.testing.assert_allclose(pa, pb, rtol=2e-3, atol=1e-3 * np.abs(pa).max())


def test_fringe_sheet_matches_coverage_model():
    """Under a wide fringe the BEM impedance sheet reproduces the Kirchhoff coverage
    model tau T_stack + (1 - tau) locally (the defining property of the equivalence)."""
    c = _lin(h=250.0, W=600.0, T=30.0, fringe_width=400.0)
    yy = np.array([450.0, 500.0, 550.0])           # tau = 0.69, 0.5, 0.31
    X = np.stack([yy, np.full(3, c.z_mid - 60.0)], 1)
    p0 = cp.field_at(_lin(h=250.0, W=None), 10.0, X)
    db_b = 20 * np.log10(np.abs(cp.field_at(c, 10.0, X, method="bem") / p0))
    db_p = 20 * np.log10(np.abs(cp.field_at(c, 10.0, X, method="po") / p0))
    assert np.max(np.abs(db_b - db_p)) < 1.5, (db_b, db_p)


def test_coverage_profile_area_preserving():
    c = _lin(h=150.0, W=1000.0, fringe_width=200.0)
    y = np.linspace(-800, 800, 160001)
    assert np.trapezoid(c.coverage(y), y) == pytest.approx(1000.0, rel=1e-4)
    assert c.coverage(np.array([-500.0, 500.0])) == pytest.approx([0.5, 0.5])


# ----------------------------------------------------------------------------- validation b)
@slow
def test_wide_shield_approaches_1d_transfer():
    """W >> h: BEM (W = 5 km) and the spectral PO (W = 20 km, exact for an infinite slab
    away from the edges) at the centre vs shield.transfer at the trace speed V."""
    g = SH.gas("helium", 216.65, 5474.9)
    X = np.array([[0.0, -250.0 - 15.0 - 300.0]])
    for f, W, meth in ((2.0, 5000.0, "bem"), (3.0, 5000.0, "bem"), (10.0, 20000.0, "po"), (20.0, 20000.0, "po")):
        cfg = _lin(h=250.0, W=W, T=15.0, gas=g)
        p = cp.field_at(cfg, f, X, method=meth)[0]
        p0 = cp.field_at(_lin(h=250.0, W=None), f, X)[0]
        T1, _ = SH.transfer(cfg.stack(15.0), [f], cfg.V, cfg.atm, reference="air")
        assert 20 * np.log10(abs(p / p0)) == pytest.approx(20 * np.log10(abs(T1[0])), abs=1.0)


# ----------------------------------------------------------------------------- validation d)
def _halfplane_hard(k, r, th, r0, th0, mmax=None):
    """Exact total field of H0^(2)(k|x - x0|) with a rigid half-plane on theta = 0
    (wedge eigenfunction series, nu = m/2)."""
    r = np.asarray(r, float)
    th = np.asarray(th, float)
    rl, rg = np.minimum(r, r0), np.maximum(r, r0)
    mmax = int(2 * (k * max(r0, r.max()) + 40)) if mmax is None else mmax
    out = np.zeros(r.shape, complex)
    for m in range(mmax + 1):
        nu = m / 2.0
        eps = 1.0 if m == 0 else 2.0
        out += eps / 2 * sp.jv(nu, k * rl) * sp.hankel2(nu, k * rg) * np.cos(nu * th) * np.cos(nu * th0)
    return 4j * (-1j * np.pi / (2 * np.pi)) * out


@slow
def test_rigid_strip_edge_vs_halfplane():
    """Thin rigid strip (1 m, 1e6 kg/m^2 faces), line source above-left of its left edge;
    field on a circle around that edge vs the exact rigid half-plane solution."""
    f = 5.0
    cfg = _lin(h=300.0, y_left=0.0, y_right=3000.0, T=1.0, m_top=1e6, m_bot=1e6, m_side=1e6)
    k = 2 * np.pi * f * np.sqrt(1 / cfg.c ** 2 - 1 / cfg.V ** 2)
    E = np.array([0.0, cfg.z_mid])
    xs = E + np.array([-150.0, 200.0])

    def inc(Y, Zc):
        dy, dz = Y - xs[0], Zc - xs[1]
        r = np.hypot(dy, dz)
        H1 = sp.hankel2(1, k * r)
        return sp.hankel2(0, k * r), -k * H1 * dy / r, -k * H1 * dz / r

    sol = cp.bem_solve(cfg, f, inc=inc)
    th = np.radians(np.arange(20.0, 341.0, 20.0))
    R = 120.0
    X = E + R * np.stack([np.cos(th), np.sin(th)], 1)     # theta from +y (along the strip)
    pb = cp.bem_field(sol, X)
    r0 = np.hypot(*(xs - E))
    th0 = np.arctan2(xs[1] - E[1], xs[0] - E[0]) % (2 * np.pi)
    ph = _halfplane_hard(k, np.full(th.size, R), th, r0, th0)
    err = np.abs(np.abs(pb) - np.abs(ph)) / np.abs(ph).max()
    assert err.max() < 0.05, err


# ----------------------------------------------------------------------------- validation e)
@slow
def test_bem_vs_kirchhoff_overlap_transparent_stack():
    """Partially transparent stack (double wall with air): below the shield the field is
    dominated by transmission, where PO + angular spectrum must agree with the BEM."""
    air = SH.gas("air", 216.65, 5474.9)
    cfg = _lin(h=250.0, W=1000.0, T=15.0, gas=air, m_top=0.4, m_bot=0.4, m_side=0.4,
               phis=np.arange(-40.0, 41.0, 10.0))
    X = cfg.circle_points()
    for f in (15.0, 20.0, 25.0):
        pb = cp.field_at(cfg, f, X, method="bem")
        pp = cp.field_at(cfg, f, X, method="po")
        assert np.max(np.abs(20 * np.log10(np.abs(pb / pp)))) < 1.5


@slow
def test_aged_incident_fit():
    from boomsim import propagate as PR
    cfg = cp.CrossPlaneConfig(CONCORDE, h=250.0, W=None, f_max=100.0, df=0.25, incident="aged")
    f = cfg.frequencies()
    A = cp.incident_coefficients(cfg, f)
    w = 2 * np.pi * f
    k = w * np.sqrt(1 / cfg.c ** 2 - 1 / cfg.V ** 2)
    for phd in (0.0, 22.5, 47.5):     # between and on the fit azimuths
        ph = np.radians(phd)
        r = float(cp.aging_radius(cfg, ph))
        X = np.array([[r * np.sin(ph), -r * np.cos(ph)]])
        P = np.array([cp.multipole_field(A[i], k[i], X[:, 0], X[:, 1])[0] if f[i] > 0 else 0 for i in range(f.size)])
        _, p = cp.spectra_to_signatures(cfg, f, P[:, None], R=r)
        ta, pa = PR.aged_near_field(CONCORDE, r, 4000.0, phi=ph)
        G = cp._signature_spectrum(ta, pa, f) * np.exp(-1j * w * cfg.beta * r / cfg.V)
        G[0] = 0
        _, pa = cp.spectra_to_signatures(cfg, f, G[:, None], R=r)
        assert np.max(np.abs(p - pa)) / np.max(np.abs(pa)) < 0.03
