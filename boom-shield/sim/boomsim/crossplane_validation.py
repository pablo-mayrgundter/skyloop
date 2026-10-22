"""Validation runs and figure for the cross-plane solver (boomsim.crossplane).

    python -m boomsim.crossplane_validation [--cache FILE] [--only NAME ...] [--workers 4]

Writes results/crossplane_validation.png and results/crossplane_validation.md.
Each block is cached (pickle) so the figure can be rebuilt without recomputing.
"""
from __future__ import annotations

import argparse
import math
import os
import pickle
import time
from types import SimpleNamespace

import numpy as np

from . import crossplane as cp
from . import shield as SH
from .source import CONCORDE, Aircraft

HERE = os.path.dirname(os.path.abspath(__file__))
RESULTS = os.path.join(os.path.dirname(HERE), "results")
T20, P20 = 216.65, 5474.9


def _ratio_db(p, p0):
    return 20 * np.log10(np.abs(p) / np.abs(p0))


# ------------------------------------------------------------------------- blocks
def block_incident():
    """Exact incident vs Whitham at phi = 0 (r = 1.5 L, 5 L) and phi = 50 deg (5 L)."""
    out = {}
    for lab, r, phd in (("1.5L", 1.5 * CONCORDE.length, 0.0), ("5L", 5 * CONCORDE.length, 0.0),
                        ("5L_phi50", 5 * CONCORDE.length, 50.0)):
        t, p, pw = cp.incident_vs_whitham(CONCORDE, r, phi=math.radians(phd))
        out[lab] = dict(t=t, p=p, pw=pw, peak_ratio=p.max() / pw.max(),
                        max_err=np.max(np.abs(p - pw)) / np.max(np.abs(pw)),
                        rms_err=np.sqrt(np.mean((p - pw) ** 2)) / np.sqrt(np.mean(pw ** 2)))
    return out


def _he30(W, **kw):
    kw.setdefault("incident", "linear")
    return cp.CrossPlaneConfig(CONCORDE, h=250.0, W=W, T=30.0, **kw)


def block_maps():
    y = np.linspace(-1000, 1000, 251)
    z = np.linspace(-1000, -40, 121)
    out = dict(y=y, z=z)
    for lab, cfg in (("W250", _he30(250.0)), ("W1000", _he30(1000.0)),
                     ("W1000_asym", _he30(None, y_left=-400.0, y_right=600.0)),
                     ("W1000_fringe", _he30(1000.0, fringe_width=150.0))):
        p = cp.field_map(cfg, 10.0, y, z)
        p0 = cp.field_map(cp.CrossPlaneConfig(CONCORDE, h=250.0, W=None, incident="linear"), 10.0, y, z)
        out[lab] = _ratio_db(p, p0)
        out[lab + "_cfg"] = (cfg.body, cfg.extent, cfg.z_mid, cfg.T)
    return out


def block_overlap():
    """BEM vs PO at circle points, 2-40 Hz: opaque He 30 m and a transparent double wall."""
    f = np.array([2, 3, 5, 7, 10, 14, 17, 20, 25, 30, 40.0])
    phis = np.array([0.0, -30.0, -55.0])
    air = SH.gas("air", T20, P20)
    cases = {"He 30 m": _he30(1000.0, phis=phis),
             "double wall (air 15 m, 0.4 kg/m2 faces)": cp.CrossPlaneConfig(
                 CONCORDE, h=250.0, W=1000.0, T=15.0, gas=air, m_top=0.4, m_bot=0.4, m_side=0.4,
                 phis=phis, incident="linear")}
    out = dict(f=f, phis=phis)
    for lab, cfg in cases.items():
        X = cfg.circle_points(716.0)
        p0 = np.array([cp.field_at(cp.CrossPlaneConfig(CONCORDE, h=250.0, W=None, incident="linear"), fi, X)
                       for fi in f])
        pb = np.array([cp.field_at(cfg, fi, X, method="bem") for fi in f])
        pp = np.array([cp.field_at(cfg, fi, X, method="po") for fi in f])
        out[lab] = dict(bem=_ratio_db(pb, p0), po=_ratio_db(pp, p0))
    return out


def _review_strip(W, h, f, mach):
    from scipy.special import fresnel
    lam = 295.07 / f / np.sqrt(1 - 1 / mach ** 2)
    v = (W / 2) * np.sqrt(2 / (lam * h))
    S, C = fresnel(v)
    return 20 * np.log10(abs(1 - 2 * (C + 1j * S) / (1 + 1j)))


def block_floor():
    """Review sec. 3.2: undertrack two-edge diffraction floor, observer 20 km below."""
    rows = []
    for W, h in ((250.0, 100.0), (250.0, 250.0), (1000.0, 250.0)):
        for f in (2.0, 10.0, 20.0):
            X = np.array([[0.0, -20000.0]])
            p0 = cp.field_at(cp.CrossPlaneConfig(CONCORDE, h=h, W=None, incident="linear"), f, X)[0]
            he = cp.CrossPlaneConfig(CONCORDE, h=h, W=W, T=30.0, incident="linear")
            rig = cp.CrossPlaneConfig(CONCORDE, h=h, W=W, T=2.0, m_top=1e6, m_bot=1e6, m_side=1e6,
                                      incident="linear")
            r = dict(W=W, h=h, f=f,
                     review_M2=_review_strip(W, h, f, 2.0), review_M14=_review_strip(W, h, f, 1.4),
                     po_black=_ratio_db(cp.field_at(rig, f, X, method="po")[0], p0),
                     bem_rigid=_ratio_db(cp.field_at(rig, f, X, method="bem")[0], p0),
                     bem_he30=_ratio_db(cp.field_at(he, f, X, method="bem")[0], p0))
            rows.append(r)
    return rows


def _beam(k, ky0, w0, y0, z0, nk=241):
    """Exact down-going Gaussian beam (plane-wave superposition) with trace wavenumber ky0."""
    dk = 2.0 / w0
    ky = ky0 + np.linspace(-6, 6, nk) * dk
    ky = ky[np.abs(ky) < k]
    kz = np.sqrt(k * k - ky * ky)
    g = np.exp(-((ky - ky0) * w0 / 2) ** 2) * (ky[1] - ky[0])

    def inc(Y, Z):
        Y = np.asarray(Y, float)[..., None]
        Z = np.asarray(Z, float)[..., None]
        e = g * np.exp(-1j * ky * (Y - y0) + 1j * kz * (Z - z0))
        return e.sum(-1), (-1j * ky * e).sum(-1), (1j * kz * e).sum(-1)
    return inc


def _flux_down(cfg, f, sol, inc, z_line, Y, n):
    y = np.linspace(-Y, Y, n)
    X = np.stack([y, np.full(n, z_line)], 1)
    dirs = np.repeat([[0.0, -1.0]], n, 0)
    w = 2 * np.pi * f
    p = cp.bem_field(sol, X)
    dp = cp.bem_grad(sol, X, dirs)
    pt = 0.5 * np.real(p * np.conj(-dp / (1j * w * cfg.rho)))
    pi, gy, gz = inc(X[:, 0], X[:, 1])
    pin = 0.5 * np.real(pi * np.conj(gz / (1j * w * cfg.rho)))
    dy = y[1] - y[0]
    return pt.sum() * dy, pin.sum() * dy


def block_cells():
    """Review sec. 3.3 emulation: V -> infinity makes the cross-plane a plain 2-D acoustic
    plane; a Gaussian beam with trace wavenumber omega/V_true (M 2) on a 3 km shield with
    rigid partitions every L = 250 m reproduces the review's x-periodic cell model (walls
    block lateral sloshing, faces free). Total transmitted power ratio, walls vs none."""
    atm = SimpleNamespace(T=T20, p=P20, rho=SH.gas("air", T20, P20).rho, c=SH.gas("air", T20, P20).c)
    he = SH.gas("helium", T20, P20)
    fake = Aircraft("2-D", length=60.0, mach=1e6, altitude=20250.0, weight=1.0, volume=1.0)
    V_true = 2.0 * atm.c
    rows = []
    for d in (15.0, 60.0):
        for f in (2.0, 10.0, 20.0):
            res = dict(d=d, f=f)
            w = 2 * np.pi * f
            T1, _ = SH.transfer(SH.Stack([SH.Membrane((atm.rho - he.rho) * d / 2), SH.GasLayer(he, d),
                                          SH.Membrane((atm.rho - he.rho) * d / 2)]), [f], V_true, atm)
            res["tmm"] = 20 * np.log10(abs(T1[0]))
            for lab, parts in (("nowall", ()), ("walls", tuple(np.arange(-1250.0, 1251.0, 250.0)))):
                cfg = cp.CrossPlaneConfig(fake, h=250.0, W=3000.0, T=d, gas=he, atm=atm,
                                          partitions=parts, m_part=1e6, incident="linear")
                k = w / atm.c
                inc = _beam(k, w / V_true, 600.0, 125.0, -250.0)
                sol = cp.bem_solve(cfg, f, inc=inc)
                lam = 2 * np.pi / k
                Pt, Pi = _flux_down(cfg, f, sol, inc, -250.0 - d - 30.0, 2600.0, int(5200 / (lam / 8)))
                res[lab] = 10 * np.log10(abs(Pt / Pi))
            rows.append(res)
    return rows


def block_edges():
    """Edge treatments (coordinator request): W = 1 km He 30 m, M 2, h = 250 m.
    (a) sharp symmetric, (b) asymmetric (-400..600 m), (c) raised-cosine fringe 150 m."""
    f = np.array([2, 3, 4, 5, 6, 8, 10, 12, 14, 17, 20, 24, 28, 34, 40.0])
    phis = np.arange(-70.0, 70.1, 10.0)
    cfgs = {"a sharp": _he30(1000.0, phis=phis, R_m=716.0),
            "b asymmetric": _he30(None, y_left=-400.0, y_right=600.0, phis=phis, R_m=716.0),
            "c fringe 150 m": _he30(1000.0, fringe_width=150.0, phis=phis, R_m=716.0),
            "c' fringe 150 m (top-level sheet)": _he30(1000.0, fringe_width=150.0, fringe_level="top",
                                                       phis=phis, R_m=716.0)}
    out = dict(f=f, phis=phis)
    X = cfgs["a sharp"].circle_points()
    p0 = np.array([cp.field_at(cp.CrossPlaneConfig(CONCORDE, h=250.0, W=None, incident="linear"), fi, X)
                   for fi in f])
    out["p0"] = p0
    for lab, cfg in cfgs.items():
        out[lab] = np.array([cp.field_at(cfg, fi, X, method="bem") for fi in f])
        out[lab + " PO"] = np.array([cp.field_at(cfg, fi, X, method="po") for fi in f])
    return out


def block_full(workers=4):
    """Matching signatures (aged incident), W = 1 km He 30 m slab vs unshielded; f_max and
    df convergence on W = 250 m."""
    out = {}
    runs = {
        "W1000": _he30(1000.0, incident="aged", f_max=100.0),
        "none": cp.CrossPlaneConfig(CONCORDE, h=250.0, W=None, incident="aged", f_max=100.0,
                                    R_m=_he30(1000.0).matching_radius()),
        "W250_f100": _he30(250.0, incident="aged", f_max=100.0),
        "W250_f200": _he30(250.0, incident="aged", f_max=200.0),
    }
    df250 = runs["W250_f100"].auto_df()
    runs["W250_f100_halfdf"] = _he30(250.0, incident="aged", f_max=100.0, df=df250 / 2)
    for lab, cfg in runs.items():
        t0 = time.time()
        sig, (f, P) = cp.matching_signatures(cfg, workers=workers, return_spectra=True)
        out[lab] = dict(sig=sig, f=f, P=P, R_m=cfg.matching_radius(), time=time.time() - t0,
                        nf=f.size, df=f[1])
        print(f"  {lab}: {f.size} freqs, {time.time() - t0:.0f} s")
    return out


def block_edges_full(workers=4):
    """Full matching runs (aged incident, 100 Hz) for the edge-treatment variants."""
    out = {}
    R = _he30(1000.0).matching_radius()
    for lab, cfg in (("b asymmetric", _he30(None, y_left=-400.0, y_right=600.0, incident="aged", R_m=R)),
                     ("c fringe 150 m", _he30(1000.0, fringe_width=150.0, incident="aged", R_m=R,
                                              df=_he30(1000.0).auto_df()))):
        t0 = time.time()
        sig, (f, P) = cp.matching_signatures(cfg, workers=workers, return_spectra=True)
        out[lab] = dict(sig=sig, f=f, P=P, time=time.time() - t0, df=f[1], nf=f.size)
        print(f"  {lab}: {f.size} freqs, {time.time() - t0:.0f} s")
    return out


BLOCKS = dict(incident=block_incident, maps=block_maps, overlap=block_overlap, floor=block_floor,
              cells=block_cells, edges=block_edges, full=block_full, edges_full=block_edges_full)


# ------------------------------------------------------------------------- figure / summary
def _band_db(P, P0, f, lo, hi):
    m = (f >= lo) & (f < hi)
    return 10 * np.log10(np.sum(np.abs(P[m]) ** 2, axis=0) / np.sum(np.abs(P0[m]) ** 2, axis=0))


def make_figure(D, path):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig = plt.figure(figsize=(18, 20))
    gs = fig.add_gridspec(5, 3, hspace=0.42, wspace=0.28)
    # 1 incident vs Whitham
    if "incident" in D:
        ax = fig.add_subplot(gs[0, 0])
        for lab, c in (("1.5L", "C0"), ("5L", "C1")):
            d = D["incident"][lab]
            m = (d["t"] > -0.05) & (d["t"] < 0.2)
            ax.plot(d["t"][m] * 1e3, d["pw"][m], c, lw=2, alpha=0.4, label=f"Whitham r={lab}")
            ax.plot(d["t"][m] * 1e3, d["p"][m], c + "--", lw=1,
                    label=f"exact Hankel r={lab} (peak {d['peak_ratio']:.3f}, max err {100 * d['max_err']:.1f}%)")
        ax.set_xlabel("retarded time [ms]"); ax.set_ylabel("p [Pa]")
        ax.set_title("Incident field vs Whitham, phi=0, band-limited 100 Hz"); ax.legend(fontsize=7)
    # 2-3 maps
    if "maps" in D:
        M = D["maps"]
        for j, lab in enumerate(("W250", "W1000")):
            ax = fig.add_subplot(gs[0, 1 + j])
            im = ax.pcolormesh(M["y"], M["z"], np.clip(M[lab], -50, 10), cmap="viridis", shading="auto",
                               vmin=-50, vmax=10)
            ax.set_aspect("equal"); ax.set_title(f"10 Hz |p/p_inc| dB, He 30 m, h 250 m, {lab}")
            ax.set_xlabel("y [m]"); ax.set_ylabel("z [m]")
            fig.colorbar(im, ax=ax, shrink=0.7)
    # 4 overlap
    if "overlap" in D:
        O = D["overlap"]
        for j, lab in enumerate([k for k in O if k not in ("f", "phis")]):
            ax = fig.add_subplot(gs[1, j])
            for i, ph in enumerate(O["phis"]):
                ax.semilogx(O["f"], O[lab]["bem"][:, i], f"C{i}-", label=f"BEM phi={ph:.0f}")
                ax.semilogx(O["f"], O[lab]["po"][:, i], f"C{i}--", label=f"PO  phi={ph:.0f}")
            ax.set_title(f"BEM vs Kirchhoff/PO on R_m=716 m: {lab}", fontsize=9)
            ax.set_xlabel("f [Hz]"); ax.set_ylabel("|p/p_inc| [dB]"); ax.legend(fontsize=7); ax.grid(alpha=0.3)
    # floor table as plot
    if "floor" in D:
        ax = fig.add_subplot(gs[1, 2])
        rows = D["floor"]
        lab = [f"{r['W']:.0f}/{r['h']:.0f}\n{r['f']:.0f}Hz" for r in rows]
        x = np.arange(len(rows))
        for key, mk in (("review_M2", "k^"), ("po_black", "C0o"), ("bem_rigid", "C1s"), ("bem_he30", "C2D")):
            ax.plot(x, [r[key] for r in rows], mk, label=key)
        ax.set_xticks(x); ax.set_xticklabels(lab, fontsize=7)
        ax.set_ylabel("undertrack |p/p_inc| at 20 km [dB]"); ax.set_title("Diffraction floor (review 3.2)")
        ax.legend(fontsize=7); ax.grid(alpha=0.3)
    # full signatures
    if "full" in D:
        F = D["full"]
        ax = fig.add_subplot(gs[2, 0])
        for i, ph in enumerate((0.0, -30.0, -60.0)):
            t, p = F["W1000"]["sig"][ph]
            t0, p0 = F["none"]["sig"][ph]
            m = (t > -0.3) & (t < 2.0)
            m0 = (t0 > -0.3) & (t0 < 2.0)
            ax.plot(t0[m0], p0[m0] / 10, f"C{i}", alpha=0.35, lw=1, label=f"unshielded/10 phi={ph:.0f}")
            ax.plot(t[m], p[m], f"C{i}", lw=1, label=f"W=1 km He 30 phi={ph:.0f}")
        ax.set_xlabel("retarded time [s]"); ax.set_ylabel("p [Pa]"); ax.legend(fontsize=7)
        ax.set_title(f"Matching signatures, R_m={F['W1000']['R_m']:.0f} m (aged incident, 100 Hz)")
        ax = fig.add_subplot(gs[2, 1])
        for lab, st in (("W250_f100", "C0-"), ("W250_f200", "C1--"), ("W250_f100_halfdf", "C2:")):
            t, p = F[lab]["sig"][0.0]
            m = (t > -0.2) & (t < 1.5)
            ax.plot(t[m], p[m], st, lw=1, label=f"{lab} (df={F[lab]['df']:.3f})")
        ax.set_title("Convergence, W=250 m He 30 m, phi=0"); ax.set_xlabel("retarded time [s]")
        ax.legend(fontsize=7)
    # cells
    if "cells" in D:
        ax = fig.add_subplot(gs[2, 2])
        rows = D["cells"]
        x = np.arange(len(rows))
        ax.plot(x, [r["tmm"] for r in rows], "k^", label="1-D TMM (ideal)")
        ax.plot(x, [r["nowall"] for r in rows], "C0o", label="BEM, no walls")
        ax.plot(x, [r["walls"] for r in rows], "C3s", label="BEM, rigid walls L=250 m")
        rev = {(15.0, 10.0): -13, (60.0, 10.0): -33, (60.0, 20.0): -44, (60.0, 2.0): -3}
        ax.plot([i for i, r in enumerate(rows) if (r["d"], r["f"]) in rev],
                [rev[(r["d"], r["f"])] for r in rows if (r["d"], r["f"]) in rev], "kx", ms=10, label="review 3.3")
        ax.set_xticks(x); ax.set_xticklabels([f"He {r['d']:.0f}\n{r['f']:.0f}Hz" for r in rows], fontsize=7)
        ax.set_ylim(-130, 5); ax.set_ylabel("transmitted power [dB]"); ax.legend(fontsize=7); ax.grid(alpha=0.3)
        ax.set_title("Cell-wall floor emulation (review 3.3), M 2")
    # edges
    if "maps" in D:
        M = D["maps"]
        for j, lab in enumerate(("W1000", "W1000_asym", "W1000_fringe")):
            ax = fig.add_subplot(gs[3, j])
            im = ax.pcolormesh(M["y"], M["z"], np.clip(M[lab], -50, 10), cmap="viridis", shading="auto",
                               vmin=-50, vmax=10)
            ax.set_aspect("equal"); ax.set_title(f"edge treatment, 10 Hz |p/p_inc| dB: {lab}")
            fig.colorbar(im, ax=ax, shrink=0.7)
    if "edges" in D:
        E = D["edges"]
        ax = fig.add_subplot(gs[4, 0])
        i0 = list(E["phis"]).index(0.0)
        for j, lab in enumerate(("a sharp", "b asymmetric", "c fringe 150 m")):
            ax.semilogx(E["f"], _ratio_db(E[lab][:, i0], E["p0"][:, i0]), f"C{j}-", label=lab + " (BEM)")
            ax.semilogx(E["f"], _ratio_db(E[lab + " PO"][:, i0], E["p0"][:, i0]), f"C{j}:", lw=1,
                        label=lab + " (PO)")
        ax.set_title("Edge treatments: centre track phi=0"); ax.set_xlabel("f [Hz]")
        ax.set_ylabel("|p/p_inc| [dB]"); ax.legend(fontsize=7); ax.grid(alpha=0.3)
        ax = fig.add_subplot(gs[4, 1])
        under = np.abs(E["phis"]) <= 50.0
        for j, lab in enumerate(("a sharp", "b asymmetric", "c fringe 150 m")):
            r = _ratio_db(E[lab], E["p0"])[:, under]
            ax.semilogx(E["f"], r.max(axis=1), f"C{j}-", label=lab + ": max over |phi|<=50")
        ax.set_xlabel("f [Hz]"); ax.set_ylabel("dB"); ax.legend(fontsize=7); ax.grid(alpha=0.3)
        ax.set_title("Edge treatments: worst azimuth under the shield")
    if "edges_full" in D and "full" in D:
        ax = fig.add_subplot(gs[4, 2])
        runs = [("a sharp", D["full"]["W1000"])] + list(D["edges_full"].items())
        for j, (lab, d) in enumerate(runs):
            t, p = d["sig"][0.0]
            m = (t > -0.3) & (t < 2.5)
            ax.plot(t[m], p[m], f"C{j}", lw=1, label=f"{lab}, phi=0")
        ax.set_xlabel("retarded time [s]"); ax.set_ylabel("p [Pa]"); ax.legend(fontsize=7)
        ax.set_title("Edge treatments: matching signature phi=0 (aged, 100 Hz)")
    fig.savefig(path, dpi=110, bbox_inches="tight")


def summary_md(D):
    L = ["# Cross-plane solver validation (`boomsim/crossplane.py`)", "",
         "Generated by `python -m boomsim.crossplane_validation`; figure `crossplane_validation.png`.", ""]
    if "incident" in D:
        L += ["## a) Incident field vs Whitham (phi = 0, band-limited 100 Hz)", "",
              "| r | peak ratio exact/Whitham | max abs error / peak | rms error |", "|---|---|---|---|"]
        for lab, d in D["incident"].items():
            L.append(f"| {lab} | {d['peak_ratio']:.3f} | {100 * d['max_err']:.1f} % | {100 * d['rms_err']:.1f} % |")
        L.append("")
    if "floor" in D:
        L += ["## f) Diffraction floor (review 3.2), undertrack observer 20 km below", "",
              "| W/h | f | review formula M2 | review M1.4 | PO black strip | BEM rigid strip | BEM He 30 m |",
              "|---|---|---|---|---|---|---|"]
        for r in D["floor"]:
            L.append(f"| {r['W']:.0f}/{r['h']:.0f} | {r['f']:.0f} | {r['review_M2']:.1f} | {r['review_M14']:.1f} | "
                     f"{r['po_black']:.1f} | {r['bem_rigid']:.1f} | {r['bem_he30']:.1f} |")
        L.append("")
    if "cells" in D:
        L += ["## f) Cell-wall floor emulation (review 3.3), He, M 2, L = 250 m", "",
              "| d | f | 1-D TMM | BEM no walls | BEM rigid walls |", "|---|---|---|---|---|"]
        for r in D["cells"]:
            L.append(f"| {r['d']:.0f} | {r['f']:.0f} | {r['tmm']:.1f} | {r['nowall']:.1f} | {r['walls']:.1f} |")
        L.append("")
    if "overlap" in D:
        O = D["overlap"]
        L += ["## e) BEM vs PO at the matching circle (dB, |p/p_inc|)", ""]
        for lab in [k for k in O if k not in ("f", "phis")]:
            d = O[lab]["bem"] - O[lab]["po"]
            L.append(f"- {lab}: BEM - PO, 14-25 Hz, phi = {[float(x) for x in O['phis']]}: "
                     + ", ".join(f"{np.abs(d[(O['f'] >= 14) & (O['f'] <= 25), i]).max():.1f}" for i in range(d.shape[1]))
                     + " dB (max |diff|)")
        L.append("")
    if "edges_full" in D and "full" in D:
        A = D["full"]["W1000"]
        phis = np.array(sorted(A["sig"].keys()))
        under = np.abs(phis) <= 50.0
        i0 = list(phis).index(0.0)
        L += ["## Edge treatments (W = 1 km He 30 m slab, M 2, h = 250 m, aged incident, 100 Hz)", "",
              "Change vs the sharp symmetric edge (a), from the full matching runs: band-mean "
              "|P|^2 on the matching circle (R_m = 716 m) at phi = 0 and for the worst azimuth "
              "under the shield (|phi| <= 50 deg), and signature peaks.", "",
              "| variant | band | centre track phi=0 | worst |phi|<=50 |", "|---|---|---|---|"]

        def band(d, lo, hi):
            m = (d["f"] >= lo) & (d["f"] < hi)
            return 10 * np.log10(np.mean(np.abs(d["P"][m]) ** 2, axis=0))
        for lab, d in D["edges_full"].items():
            for lo, hi in ((2, 5), (5, 10), (10, 20), (20, 50), (50, 100)):
                ba, bb = band(A, lo, hi), band(d, lo, hi)
                L.append(f"| {lab} | {lo}-{hi} Hz | {bb[i0] - ba[i0]:+.1f} dB | "
                         f"{bb[under].max() - ba[under].max():+.1f} dB |")
        L += ["", "| variant | peak phi=0 [Pa] | max peak over |phi|<=50 [Pa] |", "|---|---|---|"]
        for lab, d in [("a sharp", A)] + list(D["edges_full"].items()):
            pk = np.array([d["sig"][p][1].max() for p in phis])
            L.append(f"| {lab} | {pk[i0]:.2f} | {pk[under].max():.2f} |")
        L.append("")
    if "full" in D:
        F = D["full"]
        L += ["## Matching signatures and convergence", ""]
        for lab, d in F.items():
            t, p = d["sig"][0.0]
            L.append(f"- {lab}: {d['nf']} frequencies (df = {d['df']:.3f} Hz), {d['time']:.0f} s wall; "
                     f"phi=0 peak {p.max():.2f} Pa, min {p.min():.2f} Pa")
        for a, b in (("W250_f100", "W250_f200"), ("W250_f100", "W250_f100_halfdf")):
            t1, p1 = F[a]["sig"][0.0]
            t2, p2 = F[b]["sig"][0.0]
            p2i = np.interp(t1, t2, p2)
            m = (t1 > -0.2) & (t1 < 2.0)
            L.append(f"- {a} vs {b}, phi=0: max |diff| / peak = "
                     f"{np.max(np.abs(p1[m] - p2i[m])) / np.max(np.abs(p1[m])):.3f}")
        L.append("")
    return "\n".join(L)


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--cache", default=None)
    ap.add_argument("--only", nargs="*", default=None)
    ap.add_argument("--workers", type=int, default=4)
    ap.add_argument("--plot-only", action="store_true")
    a = ap.parse_args(argv)
    D = {}
    if a.cache and os.path.exists(a.cache):
        D = pickle.load(open(a.cache, "rb"))
    if not a.plot_only:
        for name, fn in BLOCKS.items():
            if a.only is not None and name not in a.only:
                continue
            t0 = time.time()
            print(f"block {name} ...", flush=True)
            D[name] = fn(a.workers) if name in ("full", "edges_full") else fn()
            print(f"block {name}: {time.time() - t0:.0f} s", flush=True)
            if a.cache:
                pickle.dump(D, open(a.cache, "wb"))
    os.makedirs(RESULTS, exist_ok=True)
    make_figure(D, os.path.join(RESULTS, "crossplane_validation.png"))
    open(os.path.join(RESULTS, "crossplane_validation.md"), "w").write(summary_md(D))
    print("written", RESULTS)


if __name__ == "__main__":
    main()
