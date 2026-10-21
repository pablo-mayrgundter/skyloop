"""Discrete line ties: Floquet-Bloch reference model (check_lineties.py) vs the smeared
GasLayer tie and vs the cross-plane BEM line ties (results/line_ties_check.md).

Regression for the line-tie fix in boomsim.crossplane: without a physical attachment
width the tie acts on one boundary element, and on the tension-free BEM faces the
result then depends on the element size (lambda_perp/10). With ``line_tie_width`` the
strip is resolved and the BEM matches the Floquet model independent of the mesh.
"""
import os
import sys
import warnings

import numpy as np
import pytest

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))

import check_lineties as CL  # noqa: E402
from boomsim import crossplane as cp  # noqa: E402
from boomsim.source import CONCORDE  # noqa: E402

FS = np.array([1.0, 3.0, 5.0, 8.0, 12.0, 20.0, 50.0])


def _db(x):
    return 20 * np.log10(np.abs(x))


@pytest.mark.parametrize("fill,mach", [("hydrogen", 2.0), ("helium", 1.4)])
@pytest.mark.parametrize("K,N", [(30.0, 0.0), (300.0, 0.0), (3000.0, 1e4)])
def test_floquet_order0_is_smeared_gaslayer(fill, mach, K, N):
    """Keeping only n = 0 reproduces shield.transfer with GasLayer(tie_K = K) exactly."""
    g, m, V = CL.setup(fill, 30.0, mach)
    Ts = CL.smeared_transfer(FS, V, g, 30.0, m, K, N)
    Tf = np.array([CL.floquet(f, V, g, 30.0, m, K, 25.0, N=N, smeared=True)[0] for f in FS])
    assert np.max(np.abs(Tf - Ts) / np.abs(Ts)) < 1e-9


@pytest.mark.parametrize("N", [0.0, 1e4])
def test_floquet_full_coverage_is_smeared(N):
    """Ties attached over the whole pitch (b = P) load only n = 0: smeared result."""
    g, m, V = CL.setup("hydrogen", 30.0, 2.0)
    Ts = CL.smeared_transfer(FS, V, g, 30.0, m, 300.0, N)
    Tf = np.array([CL.floquet(f, V, g, 30.0, m, 300.0, 25.0, b=25.0, N=N)[0] for f in FS])
    assert np.max(np.abs(Tf - Ts) / np.abs(Ts)) < 1e-9


def test_floquet_small_pitch_limit():
    """P -> 0 recovers the smeared tie when the faces carry tension (the n != 0 face
    compliance ~ P^2/N vanishes), but NOT for limp faces with b/P fixed (the tie then
    sees the strip mass m b: correction -2 K (P/b - 1)/(w^2 m), independent of P)."""
    g, m, V = CL.setup("hydrogen", 30.0, 2.0)
    for f in (5.0, 10.0, 20.0):
        Ts = CL.smeared_transfer([f], V, g, 30.0, m, 300.0, 1e4)[0]
        Tf = CL.floquet(f, V, g, 30.0, m, 300.0, 0.25, b=0.025, N=1e4)[0]
        assert abs(_db(Tf) - _db(Ts)) < 0.3
    T1 = CL.floquet(5.0, V, g, 30.0, m, 300.0, 0.5, b=0.25)[0]
    T2 = CL.floquet(5.0, V, g, 30.0, m, 300.0, 0.05, b=0.025)[0]
    Ts = CL.smeared_transfer([5.0], V, g, 30.0, m, 300.0)[0]
    assert abs(_db(T1) - _db(T2)) < 0.3 and abs(_db(T1) - _db(Ts)) > 5.0


def test_limp_line_tie_decouples_as_width_vanishes():
    """N = 0: a line tie of vanishing attachment width leaves the stack tie-free
    (a limp face has no line stiffness to carry the tie force)."""
    g, m, V = CL.setup("hydrogen", 30.0, 2.0)
    for f in (5.0, 10.0):
        T0 = CL.smeared_transfer([f], V, g, 30.0, m, 0.0)[0]
        Tl = CL.floquet(f, V, g, 30.0, m, 300.0, 25.0, b=0.005)[0]
        assert abs(_db(Tl) - _db(T0)) < 1.0


def test_line_tie_legacy_warns():
    cfg = cp.CrossPlaneConfig(CONCORDE, W=200.0, T=30.0, line_ties=[(0.0, 100.0)], phis=[0.0])
    with pytest.warns(RuntimeWarning, match="element size"):
        cp.build_mesh(cfg, 5.0)
    cfg.line_tie_width = 1.0
    with warnings.catch_warnings():
        warnings.simplefilter("error", RuntimeWarning)
        mesh = cp.build_mesh(cfg, 5.0)
    yn = np.r_[mesh["P"].a[:mesh["nf"], 0], mesh["P"].b[mesh["nf"] - 1, 0]]
    assert np.min(np.abs(yn - 0.5)) < 1e-9 and np.min(np.abs(yn + 0.5)) < 1e-9


def test_bem_line_tie_width_resolved_matches_floquet():
    """W = 1 km H2 30 m slab, K = 300 Pa/m equivalent at 25 m pitch, 5 Hz plane wave.
    Legacy (one panel of ~6.8 m): BEM equals Floquet with b = panel size and is nearly
    transparent (the strip resonance). With b = 1.5 m resolved: mesh-independent and
    equal to Floquet(b = 1.5 m)."""
    g, m, V = CL.setup("hydrogen", 30.0, 2.0)
    f = 5.0
    T_leg, b_leg, _ = CL.bem_T0(f, 300.0, 25.0, "line")
    T_fl = CL.floquet(f, V, g, 30.0, m, 300.0, 25.0, b=b_leg)[0]
    assert abs(_db(T_leg) - _db(T_fl)) < 2.0
    T_a, _, _ = CL.bem_T0(f, 300.0, 25.0, "line", b=1.5)
    T_b, _, _ = CL.bem_T0(f, 300.0, 25.0, "line", b=1.5, h=3.0)
    T_fl15 = CL.floquet(f, V, g, 30.0, m, 300.0, 25.0, b=1.5)[0]
    assert abs(_db(T_a) - _db(T_b)) < 0.5
    assert abs(_db(T_a) - _db(T_fl15)) < 1.0
    assert _db(T_a) < _db(T_leg) - 10.0
