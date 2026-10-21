"""Along-track-varying (3-D) shield study: serrated edges and module joints.

Questions (see results/alongtrack_summary.md):
  Q1. Plan-view serrated edges (period Lambda, depth w_f = 150 m; sawtooth, symmetric
      raised-cosine teeth, random-phase teeth): how does the true 3-D field compare with
      the x-averaged coverage ("fringe") model, how much energy goes into non-specular
      orders k_x,n = omega/V + 2 pi n/Lambda, and which of them reach the ground?
  Q2. Transverse module joints every L_m = 250 / 500 m, gap g = 1, 5, 20 m, open or
      face-membranes-only: how much leaks and where does it land?

Case: Concorde M2, h = 250 m (shield top 20 km), W = 1 km, H2 30 m slab (the cached
run_3d configuration C2_h250_W1000_H230 and its _fr150 / _R8000 variants).

Pipeline (method: boomsim/alongtrack.py docstring):
  specular / n = 0 part: the cached cross-plane BEM matching-circle spectra of run_3d
      (sharp edge; x-averaged raised-cosine 150 m fringe), plus, for joints, the PO n = 0
      leak (g/L) * Rayleigh[D_gap - D_stack] added coherently on the circle; each
      signature then goes through the unchanged nonlinear pipeline
      (spectra_to_signatures, crop, propagate_from_circle with age_from) -> ground
      signature -> 1/3-octave band energies (raw and rise-time proxy);
  scattered part (n != 0): PO aperture orders -> angular spectrum -> ray-tube flux mapping
      to the ground (linear, with absorption) -> x-averaged ground ESD, added to the band
      energies -> Mark VII PL; headline delta = less favourable of raw / proxy vs the
      unshielded run_3d baseline (as run_3d).

Stages (cached under results/alongtrack_cache/, resumable):
  valid     validations (PO vs crossplane PO/BEM on the circle, flux mapping vs the
            linear circle pipeline, n = 0 of a serration vs the x-averaged fringe)
  scatter   non-specular orders of the serrated fringes (all cases)
  joints    joint orders + n = 0 leak spectra on the matching circles
  specular  nonlinear ground propagation of the n = 0 parts
  report    tables (stdout), results/alongtrack_results.csv, alongtrack_budget.csv,
            figures results/alongtrack_*.png

Run:  python boom-shield/sim/run_alongtrack.py [stage ...] [--procs 2]
Uses at most --procs worker processes (default 2) with single-threaded BLAS.
"""
from __future__ import annotations

import os

os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
os.environ.setdefault("OMP_NUM_THREADS", "1")
os.environ.setdefault("MKL_NUM_THREADS", "1")

import argparse  # noqa: E402
import csv  # noqa: E402
import math  # noqa: E402
import pickle  # noqa: E402
import sys  # noqa: E402
import time  # noqa: E402

import numpy as np  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import run_3d as R3  # noqa: E402  (read-only use: configs, cache, crop, constants)
from boomsim import alongtrack as AT  # noqa: E402
from boomsim import crossplane as cp  # noqa: E402
from boomsim import propagate as PR  # noqa: E402
from boomsim.proxy import scored_metrics  # noqa: E402

RESULTS = os.path.join(HERE, "results")
CACHE = os.path.join(RESULTS, "alongtrack_cache")
BASE = "C2_h250_W1000_H230"
FRINGE = BASE + "_fr150"
NONE = "C2_h250_none"
PROCS = 2

# scattered-field frequency grid (non-uniform; 1/3-octave band integration is trapezoidal)
F_S = np.concatenate([np.arange(0.25, 20.0, 0.25), np.arange(20.0, 100.0, 0.5),
                      np.arange(100.0, 200.0, 1.0), np.arange(200.0, 400.001, 2.0)])
BLOCK_EDGES = [5.0, 10.0, 20.0, 40.0, 70.0, 100.0, 150.0, 200.0, 300.0, 400.0]
Y_GRID = np.arange(-60000.0, 60000.1, 250.0)
CARPET_EDGE = None   # set in main (lateral cut-off half-width of the unshielded carpet)

LAMS = [1.0, 3.0, 10.0, 30.0, 100.0, 300.0]
SHAPES = ["sawtooth", "rc_teeth", "random"]
SCAT_CASES = [(s, L, 150.0, 0) for s in SHAPES for L in LAMS]
SCAT_CASES += [("random", 30.0, 150.0, 1), ("random", 300.0, 150.0, 1),      # 2nd realisation
               ("linear_sawtooth", 10.0, 150.0, 0), ("linear_sawtooth", 100.0, 150.0, 0),
               ("sawtooth", 3.0, 75.0, 0), ("sawtooth", 3.0, 300.0, 0),
               ("sawtooth", 10.0, 75.0, 0), ("sawtooth", 10.0, 300.0, 0)]
SCORED_LAMS = [10.0, 30.0, 100.0, 300.0, 1.0, 3.0]

JOINT_LG = [(250.0, 1.0), (250.0, 5.0), (250.0, 20.0), (500.0, 1.0), (500.0, 5.0), (500.0, 20.0)]
JOINT_KINDS = ["open", "membranes"]
PHI_J = [0.0, 8.0, 16.0, 24.0, 33.0, 42.0, 51.0]
BANDS_REPORT = [(2, 5), (5, 10), (10, 20), (20, 50), (50, 100), (100, 200), (200, 400)]


# ----------------------------------------------------------------------------- helpers
def _p(name):
    return os.path.join(CACHE, name)


def load(name):
    p = _p(name)
    if os.path.exists(p):
        with open(p, "rb") as fh:
            return pickle.load(fh)
    return None


def save(name, obj):
    os.makedirs(CACHE, exist_ok=True)
    p = _p(name)
    with open(p + ".tmp", "wb") as fh:
        pickle.dump(obj, fh)
    os.replace(p + ".tmp", p)


_CONF = None


def confs():
    global _CONF
    if _CONF is None:
        _CONF = {c["key"]: c for c in R3.build_configs()}
    return _CONF


def base_cfg(key=BASE):
    cfg = R3.make_cfg(confs()[key])
    cfg.workers = 1
    return cfg


def rays_for(cfg):
    return AT.GroundRays(cfg.aircraft.altitude + cfg.z_mid, z_anchor=cfg.aircraft.altitude)


def incident_S():
    A = load("incident_FS.pkl")
    if A is None:
        t0 = time.time()
        A = cp.incident_coefficients(base_cfg(), F_S)
        save("incident_FS.pkl", A)
        print(f"  incident coefficients on F_S: {time.time() - t0:.0f} s", flush=True)
    return A


def pool(n):
    import multiprocessing as mp
    return mp.get_context("fork").Pool(n, initializer=cp._limit_threads)


def case_name(shape, Lam, w_f, seed):
    return f"scat_{shape}_L{Lam:g}_w{w_f:g}_s{seed}"


def case_label(shape, Lam, w_f=150.0):
    nm = {"sawtooth": "sawtooth", "rc_teeth": "raised-cosine teeth", "random": "random-phase teeth",
          "linear_sawtooth": "straight sawtooth"}[shape]
    return f"{nm}, Lambda {Lam:g} m, w_f {w_f:g} m"


# ----------------------------------------------------------------------------- scatter
_W = {}
_BLK = {}


F_HF = np.concatenate([F_S, np.arange(404.0, 1000.001, 4.0)])
HF_EDGES = BLOCK_EDGES + [500.0, 600.0, 700.0, 800.0, 900.0, 1000.0]
HF_CASES = [("sawtooth", 1.0, 150.0, 0), ("sawtooth", 2.0, 150.0, 0), ("sawtooth", 3.0, 150.0, 0),
            ("random", 1.0, 150.0, 0), ("random", 3.0, 150.0, 0)]


def _fr_job(i):
    F = _W.get("F", F_S)
    f = float(F[i])
    serr, cfg = _W["serr"], _W["cfg"]
    e = next(x for x in _W.get("edges", BLOCK_EDGES) if x >= f - 1e-9)
    if e not in _BLK:
        _BLK.clear()           # one block at a time (jobs arrive in descending f)
        _BLK[e] = AT.fringe_block(serr, e, cfg.c, cfg.V)
    o = AT.fringe_orders(cfg, serr, f, _W["rays"], Y_GRID, a=_W["A"][i], block=_BLK[e])
    return i, o


def incident_HF():
    A = load("incident_FHF.pkl")
    if A is None:
        A = cp.incident_coefficients(base_cfg(), F_HF)
        save("incident_FHF.pkl", A)
    return A


def scatter_case(shape, Lam, w_f=150.0, seed=0, procs=PROCS, hf=False):
    name = case_name(shape, Lam, w_f, seed) + ("_hf" if hf else "") + ".pkl"
    d = load(name)
    if d is not None:
        return d
    cfg = base_cfg()
    serr = AT.Serration(Lam, w_f, shape, seed=seed)
    F = F_HF if hf else F_S
    _W.update(serr=serr, cfg=cfg, rays=rays_for(cfg), A=incident_HF() if hf else incident_S(),
              F=F, edges=HF_EDGES if hf else BLOCK_EDGES)
    _BLK.clear()
    t0 = time.time()
    res = [None] * F.size
    jobs = list(range(F.size))[::-1]
    with pool(procs) as pl:
        for i, o in pl.imap_unordered(_fr_job, jobs, chunksize=4):
            res[i] = o
    keys = ["ap_D", "ap_spec", "ap_scat", "ap_prop", "ap_ground", "y_mean", "sx_mean", "sy_mean"]
    d = dict(shape=shape, Lam=Lam, w_f=w_f, seed=seed, f=F, y=Y_GRID,
             E=np.array([o["E"] for o in res]), E0=np.array([o["E0"] for o in res]),
             **{k: np.array([o[k] for o in res]) for k in keys}, wall=time.time() - t0)
    save(name, d)
    print(f"  [{name}] {time.time() - t0:.0f} s", flush=True)
    return d


# ----------------------------------------------------------------------------- joints
def _jt_job(i):
    f = float(F_S[i])
    cfg = _W["cfg"]
    out = {}
    for kind in JOINT_KINDS:
        js = [AT.Joint(L, g, kind) for (L, g) in JOINT_LG]
        out[kind] = AT.joint_orders(cfg, js, f, _W["rays"], Y_GRID, a=_W["A"][i])
    return i, out


def joint_scatter(procs=PROCS):
    d = load("joints_scat.pkl")
    if d is not None:
        return d
    cfg = base_cfg()
    _W.update(cfg=cfg, rays=rays_for(cfg), A=incident_S())
    t0 = time.time()
    res = [None] * F_S.size
    with pool(procs) as pl:
        for i, o in pl.imap_unordered(_jt_job, list(range(F_S.size))[::-1], chunksize=2):
            res[i] = o
    d = dict(f=F_S, y=Y_GRID, wall=time.time() - t0)
    keys = ["ap_D", "ap_spec", "ap_scat", "ap_prop", "ap_ground"]
    for kind in JOINT_KINDS:
        for jj, (L, g) in enumerate(JOINT_LG):
            rr = [res[i][kind][jj] for i in range(F_S.size)]
            d[(kind, L, g)] = dict(E=np.array([o["E"] for o in rr]),
                                   **{k: np.array([o[k] for o in rr]) for k in keys})
    save("joints_scat.pkl", d)
    print(f"  joint orders: {time.time() - t0:.0f} s", flush=True)
    return d


def _jl_job(i):
    f = float(_W["f"][i])
    cfg = _W["cfg"]
    k = 2 * np.pi * f * math.sqrt(1 / cfg.c ** 2 - 1 / cfg.V ** 2)
    out = {}
    for kind in JOINT_KINDS:
        y, dy, dD = AT.joint_aperture(cfg, f, kind, a=_W["A"][i], dy_lambda=10.0)
        out[kind] = {R: AT.rayleigh_2d(dD, y, dy, cfg.z_mid, k, X) for R, X in _W["X"].items()}
    return i, out


def joint_leak(procs=PROCS):
    """n = 0 joint leak per unit open fraction, Q = Rayleigh[D_gap - D_stack], on the
    matching circles of the cached run_3d spectra (R_m = 3 km, all phi; 8 km, near track)."""
    d = load("joints_leak.pkl")
    if d is not None:
        return d
    s3, s8 = R3._load(BASE, "spec"), R3._load(BASE + "_R8000", "spec")
    cfg = base_cfg()
    f = s3["f"]
    assert np.allclose(f, s8["f"])
    X = {}
    for s in (s3, s8):
        ph = np.radians(s["phis"])
        X[float(s["R"])] = np.stack([s["R"] * np.sin(ph), -s["R"] * np.cos(ph)], 1)
    A = cp.incident_coefficients(cfg, f)
    _W.update(cfg=cfg, f=f, X=X, A=A)
    t0 = time.time()
    Q = {kind: {R: np.zeros((f.size, X[R].shape[0]), complex) for R in X} for kind in JOINT_KINDS}
    with pool(procs) as pl:
        for i, o in pl.imap_unordered(_jl_job, list(range(1, f.size))[::-1], chunksize=4):
            for kind in JOINT_KINDS:
                for R in X:
                    Q[kind][R][i] = o[kind][R]
    d = dict(f=f, Q=Q, phis={float(s3["R"]): s3["phis"], float(s8["R"]): s8["phis"]},
             wall=time.time() - t0)
    save("joints_leak.pkl", d)
    print(f"  joint n=0 leak on circles: {time.time() - t0:.0f} s", flush=True)
    return d


# ----------------------------------------------------------------------------- specular
_PW = {}


def _prop_job(args):
    tag, phi_d, R, Pcol = args
    t_start = time.time()
    cfg = _PW["cfg"][R]
    f = _PW["f"]
    phi = math.radians(phi_d)
    t, p = cp.spectra_to_signatures(cfg, f, Pcol[:, None])
    tc, pc = R3.crop(t, p[:, 0])
    tg, pg, yg = PR.propagate_from_circle(cfg.aircraft, phi, tc, pc, R, z_ref=R3.H_S,
                                          age_from=float(cp.aging_radius(cfg, phi)))
    row = dict(tag=tag, phi=phi_d, R=R, y_g=float(yg), wall=0.0)
    if tg is None:
        row.update(metrics=None, bands=None)
        return row
    row["metrics"] = {k: float(v) for k, v in scored_metrics(tg, pg).items()}
    row["bands"] = AT.signature_bands(tg, pg)
    if phi_d in (0.0, 24.0, 39.0):
        s = (tg >= -0.2) & (tg <= 1.6)
        row["sig"] = (tg[s][::5].astype(np.float32), pg[s][::5].astype(np.float32))
    row["wall"] = time.time() - t_start
    return row


def specular_jobs():
    """All (tag, phi, R, spectrum) n = 0 propagation jobs."""
    jobs = []
    s3, s8 = R3._load(BASE, "spec"), R3._load(BASE + "_R8000", "spec")
    f3, f8 = R3._load(FRINGE, "spec"), R3._load(FRINGE + "_R8000", "spec")
    for s, tag in ((f3, "fringe"), (f8, "fringe")):
        for j, ph in enumerate(s["phis"]):
            jobs.append((tag, float(ph), float(s["R"]), s["P"][:, j]))
    for s, phs in ((s3, (0.0, 12.0, 24.0, 39.0)), (s8, (0.0,))):           # reproduction check
        for j, ph in enumerate(s["phis"]):
            if float(ph) in phs:
                jobs.append(("sharp", float(ph), float(s["R"]), s["P"][:, j]))
    lk = joint_leak()
    for kind in JOINT_KINDS:
        for (L, g) in JOINT_LG:
            tag = f"joint_{kind}_L{L:g}_g{g:g}"
            for s, phs in ((s3, PHI_J), (s8, (0.0,))):
                Q = lk["Q"][kind][float(s["R"])]
                for j, ph in enumerate(s["phis"]):
                    if float(ph) in phs:
                        jobs.append((tag, float(ph), float(s["R"]), s["P"][:, j] + (g / L) * Q[:, j]))
    return jobs, s3["f"]


def specular(procs=PROCS, only=None):
    have = load("specular.pkl") or {}
    jobs, f = specular_jobs()
    todo = [j for j in jobs if (j[0], j[1], j[2]) not in have and (only is None or j[0].startswith(only))]
    if not todo:
        return have
    _PW.update(f=f, cfg={3000.0: base_cfg(BASE), 8000.0: base_cfg(BASE + "_R8000")})
    for c in _PW["cfg"].values():
        c.fs = R3.FS
    todo.sort(key=lambda j: (-j[1], j[0]))      # large phi (slow) first
    t0 = time.time()
    with pool(procs) as pl:
        for k, row in enumerate(pl.imap_unordered(_prop_job, todo)):
            have[(row["tag"], row["phi"], row["R"])] = row
            if k % 4 == 3 or k == len(todo) - 1:
                save("specular.pkl", have)
            print(f"    {row['tag']} phi={row['phi']:g} R={row['R']:g}: {row['wall']:.0f} s "
                  f"[{k + 1}/{len(todo)}, {time.time() - t0:.0f} s]", flush=True)
    save("specular.pkl", have)
    return have


# ----------------------------------------------------------------------------- validation
def validate(procs=PROCS):
    v = load("valid.pkl")
    if v is not None:
        return v
    v = {}
    cfg = base_cfg()
    # (a) PO difference field = crossplane.po_solve
    f = 30.0
    a = cp.incident_coefficients(cfg, [f])[0]
    ap = cp.po_solve(cfg, f, a)
    D = AT.stack_difference(cfg, f, ap["y"], a=a)
    v["po_match"] = float(np.max(np.abs(D - ap["d"])) / np.max(np.abs(ap["d"])))
    # (b) n = 0 of a true serration = crossplane's x-averaged fringe PO (circle, R 3 km)
    cfr = base_cfg(FRINGE)
    s3 = R3._load(BASE, "spec")
    ph = np.radians(s3["phis"])
    X = np.stack([3000 * np.sin(ph), -3000 * np.cos(ph)], 1)
    rows = []
    for f in (12.0, 30.0, 60.0):
        a = cp.incident_coefficients(cfr, [f])[0]
        k = 2 * np.pi * f * math.sqrt(1 / cfr.c ** 2 - 1 / cfr.V ** 2)
        p_cp = cp.field_at(cfr, f, X, method="po", a=a) - cp.multipole_field(a, k, X[:, 0], X[:, 1])
        # 3-D PO: body (tau = 1) + both fringes with c_0(y) of a sawtooth (exact intervals)
        apf = cp.po_solve(cfr, f, a)
        y = apf["y"]
        Dy = AT.stack_difference(cfr, f, y, a=a)
        s = AT.Serration(30.0, cfr.fringe_width, "sawtooth")
        y0, y1 = cfr.body
        u = np.where(y > y1, (y - y1) / s.w_f, np.where(y < y0, (y0 - y) / s.w_f, -1.0))
        c0 = np.ones_like(y)
        m = u >= 0
        c0[m] = s.coefficients(np.clip(u[m], 0, 1), [0])[0].real
        p_at = AT.rayleigh_2d(c0 * Dy, y, apf["dy"], cfr.z_mid, k, X)
        rows.append((f, float(np.max(np.abs(p_at - p_cp)) / np.max(np.abs(p_cp)))))
    v["n0_vs_fringe_po"] = rows
    # (c) sharp edge and x-averaged fringe: PO vs BEM (cached) on the circle, band levels
    sf3 = R3._load(FRINGE, "spec")
    fc = s3["f"]
    sel = np.nonzero((fc >= 5.0) & (fc <= 100.0))[0][::2]
    Ab = cp.incident_coefficients(cfg, fc[sel])
    po = {"sharp": np.zeros((sel.size, X.shape[0]), complex), "fringe": np.zeros((sel.size, X.shape[0]), complex)}
    inc = np.zeros((sel.size, X.shape[0]), complex)
    for ii, i in enumerate(sel):
        f = float(fc[i])
        k = 2 * np.pi * f * math.sqrt(1 / cfg.c ** 2 - 1 / cfg.V ** 2)
        inc[ii] = cp.multipole_field(Ab[ii], k, X[:, 0], X[:, 1])
        po["sharp"][ii] = cp.field_at(cfg, f, X, method="po", a=Ab[ii])
        po["fringe"][ii] = cp.field_at(cfr, f, X, method="po", a=Ab[ii])
    bem = {"sharp": s3["P"][sel], "fringe": sf3["P"][sel]}
    v["po_bem"] = dict(f=fc[sel], phis=s3["phis"], po=po, bem=bem, inc=inc)
    # (d) flux mapping vs the linear circle pipeline (unshielded incident)
    v["flux_vs_pipeline"] = flux_check()
    save("valid.pkl", v)
    return v


def flux_check():
    cfg = base_cfg()
    cfg0 = base_cfg(NONE)
    cfg0.fs = R3.FS
    sp_ = R3._load(NONE, "spec")
    f, P = sp_["f"], sp_["P"]
    rays = rays_for(cfg)
    fsel = [10.0, 30.0, 60.0]
    idx = [int(round(x / sp_["df"])) for x in fsel]
    A = cp.incident_coefficients(cfg, f[idx])
    rows = []
    for j, phd in enumerate(sp_["phis"]):
        if phd not in (0, 12, 24, 39):
            continue
        t, p = cp.spectra_to_signatures(cfg0, f, P[:, [j]])
        tc, pc = R3.crop(t, p[:, 0])
        tg, pg, yg = PR.propagate_from_circle(cfg0.aircraft, math.radians(phd), tc, pc, sp_["R"],
                                              z_ref=R3.H_S, linear=True)
        dt = tg[1] - tg[0]
        out = []
        for ii, fi in enumerate(fsel):
            Pg = np.sum(pg * np.exp(-2j * np.pi * fi * tg)) * dt
            w = 2 * np.pi * fi
            k = w * math.sqrt(1 / cfg.c ** 2 - 1 / cfg.V ** 2)
            th = np.linspace(-1.4, 1.4, 20001)
            ky, kz = k * np.sin(th), k * np.cos(th)
            Dt2 = 4 * np.abs(cp._spectral_weight(A[ii], th)) ** 2 / kz ** 2
            dy, G, ok = rays.gain(1 / cfg.V, ky / w, fi)
            E = np.interp(yg, (dy + abs(cfg.z_mid) * np.tan(th))[ok], (Dt2 * G)[ok])
            out.append(10 * np.log10(E / abs(Pg) ** 2))
        rows.append((float(phd), float(yg), out))
    return dict(f=fsel, rows=rows)


# ----------------------------------------------------------------------------- analysis
def unshielded_ref():
    g = R3._load(NONE, "ground")
    return {float(k[0]): v["metrics"] for k, v in g.items() if k[1] == "full" and k[2] == R3.CROP_REL}


def smooth_y(E, n=3):
    """Boxcar over n ground-grid points (750 m) along the last axis."""
    if n <= 1:
        return E
    ker = np.ones(n) / n
    return np.apply_along_axis(lambda r: np.convolve(r, ker, mode="same"), -1, E)


def scat_E(d, mirror=True):
    """x-averaged scattered |P_g|^2 (f, y) of both fringes (left = mirror of right)."""
    E = d["E"]
    if mirror:
        E = E + E[:, ::-1]
    return smooth_y(E)


def E_at(E, y):
    """E(f, y) interpolated at ground offset y (|y| symmetric)."""
    yy = abs(y)
    j = np.searchsorted(Y_GRID, yy)
    j = int(np.clip(j, 1, Y_GRID.size - 1))
    w = (yy - Y_GRID[j - 1]) / (Y_GRID[j] - Y_GRID[j - 1])
    return (1 - w) * E[:, j - 1] + w * E[:, j]


BL_VARIANTS = {
    "full": lambda f: np.ones_like(f),                          # scattered up to 400 Hz
    "200": lambda f: (f <= 200.0).astype(float),
    "100": lambda f: AT.crossplane_taper(f),                    # band limit of the circle spectra
}


def score(spec_row, E_f, ref, variant="full"):
    """Combined metrics at one ground point: specular (n = 0) row + scattered E(f)."""
    if spec_row["bands"] is None:
        return None
    Ew = E_f * BL_VARIANTS[variant](F_S)
    raw, prx = AT.combined_pl(spec_row["bands"], F_S, Ew)
    dr, dp = raw - ref["pldb"], prx - ref["pldb_proxy"]
    return dict(pldb=raw, pldb_proxy=prx, d_raw=dr, d_proxy=dp, d=max(dr, dp))


def scat_only_pl(E_f, variant="full"):
    z = [np.zeros(42), np.zeros(42)]
    return AT.combined_pl(z, F_S, E_f * BL_VARIANTS[variant](F_S))


def band_sum(f, x, lo, hi):
    """int x df over [lo, hi) (trapezoid on the non-uniform grid)."""
    m = (f >= lo) & (f < hi)
    if np.count_nonzero(m) < 2:
        return 0.0
    return float(np.trapezoid(x[m], f[m]))


# ----------------------------------------------------------------------------- report
def report():
    os.makedirs(RESULTS, exist_ok=True)
    spec = load("specular.pkl") or {}
    ref = unshielded_ref()
    g_sharp = R3._load(BASE, "ground")
    g_fr = R3._load(FRINGE, "ground")
    g_sharp8 = R3._load(BASE + "_R8000", "ground")
    g_fr8 = R3._load(FRINGE + "_R8000", "ground")
    phc, ycut = PR.lateral_cutoff(base_cfg().aircraft.altitude, base_cfg().V)
    print(f"\nunshielded carpet lateral cut-off: y = {ycut / 1e3:.1f} km (phi {math.degrees(phc):.1f} deg)")

    # ---- reproduction checks
    print("\n## Reproduction of run_3d ground metrics (re-propagated cached spectra)")
    print("| case | R_m | phi | PLdB raw: here / 3d_results | proxy: here / 3d_results |")
    print("|---|---|---|---|---|")
    for (tag, ph, R), row in sorted(spec.items()):
        if tag not in ("sharp", "fringe") or row["metrics"] is None:
            continue
        g = {("sharp", 3000.0): g_sharp, ("sharp", 8000.0): g_sharp8,
             ("fringe", 3000.0): g_fr, ("fringe", 8000.0): g_fr8}[(tag, R)]
        r0 = g.get((ph, "full", R3.CROP_REL))
        if r0 is None or (tag == "fringe" and ph not in (0.0, 24.0, 39.0, 54.0) and R == 3000.0):
            continue
        m, m0 = row["metrics"], r0["metrics"]
        print(f"| {tag} | {R / 1e3:g} km | {ph:g} | {m['pldb']:.2f} / {m0['pldb']:.2f} | "
              f"{m['pldb_proxy']:.2f} / {m0['pldb_proxy']:.2f} |")

    rows_csv = []

    def spec_row(tag, ph, R):
        return spec.get((tag, ph, R))

    # ---- serrations: per-phi combined
    print("\n## Serrated edges: ground Delta PLdB (headline = less favourable of raw/proxy), R_m 3 km")
    fr_phis = sorted(ph for (t, ph, R) in spec if t == "fringe" and R == 3000.0)
    cases = {}
    for (shape, Lam, w_f, seed) in SCAT_CASES:
        d = load(case_name(shape, Lam, w_f, seed) + ".pkl")
        if d is not None:
            cases[(shape, Lam, w_f, seed)] = d
    summ = []

    def add_rows(label, key, tag, E, phis, R, extra):
        out = []
        for ph in phis:
            sr = spec_row(tag, ph, R)
            if sr is None or sr["metrics"] is None:
                continue
            Ef = E_at(E, sr["y_g"]) if E is not None else np.zeros(F_S.size)
            sc = score(sr, Ef, ref[ph])
            sc100 = score(sr, Ef, ref[ph], "100")
            sc200 = score(sr, Ef, ref[ph], "200")
            so = scat_only_pl(Ef)
            row = dict(key=key, label=label, R_m=R, phi_deg=ph, y_g_m=sr["y_g"],
                       pldb_n0=sr["metrics"]["pldb"], pldb_proxy_n0=sr["metrics"]["pldb_proxy"],
                       pldb=sc["pldb"], pldb_proxy=sc["pldb_proxy"], d_raw=sc["d_raw"],
                       d_proxy=sc["d_proxy"], d_headline=sc["d"], d_headline_scat100=sc100["d"],
                       d_headline_scat200=sc200["d"], pldb_scat_only=so[0], pldb_proxy_scat_only=so[1],
                       **extra)
            out.append(row)
        return out

    # references
    for tag, label in (("sharp_cached", "sharp edge (run_3d)"), ("fringe", "x-averaged raised-cosine fringe 150 m")):
        pass
    ref_rows = []
    for (ph, b, cr), r in sorted(g_sharp.items()):
        if b != "full" or cr != R3.CROP_REL or r["metrics"] is None:
            continue
        dr = r["metrics"]["pldb"] - ref[ph]["pldb"]
        dp = r["metrics"]["pldb_proxy"] - ref[ph]["pldb_proxy"]
        ref_rows.append(dict(key="sharp", label="sharp edge (run_3d cache)", R_m=3000.0, phi_deg=ph,
                             y_g_m=r["y_g"], pldb_n0=r["metrics"]["pldb"],
                             pldb_proxy_n0=r["metrics"]["pldb_proxy"], pldb=r["metrics"]["pldb"],
                             pldb_proxy=r["metrics"]["pldb_proxy"], d_raw=dr, d_proxy=dp,
                             d_headline=max(dr, dp), d_headline_scat100=max(dr, dp),
                             d_headline_scat200=max(dr, dp), pldb_scat_only=float("-inf"),
                             pldb_proxy_scat_only=float("-inf"), shape="", Lam_m="", w_f_m="",
                             seed="", joint_kind="", L_m="", g_m=""))
    for (ph, b, cr), r in sorted(g_sharp8.items()):
        if b != "full" or cr != R3.CROP_REL or r["metrics"] is None:
            continue
        dr = r["metrics"]["pldb"] - ref[ph]["pldb"]
        dp = r["metrics"]["pldb_proxy"] - ref[ph]["pldb_proxy"]
        ref_rows.append(dict(ref_rows[0], R_m=8000.0, phi_deg=ph, y_g_m=r["y_g"],
                             pldb_n0=r["metrics"]["pldb"], pldb_proxy_n0=r["metrics"]["pldb_proxy"],
                             pldb=r["metrics"]["pldb"], pldb_proxy=r["metrics"]["pldb_proxy"],
                             d_raw=dr, d_proxy=dp, d_headline=max(dr, dp),
                             d_headline_scat100=max(dr, dp), d_headline_scat200=max(dr, dp)))
    rows_csv += ref_rows
    blank = dict(shape="", Lam_m="", w_f_m="", seed="", joint_kind="", L_m="", g_m="")
    fr_rows = add_rows("x-averaged raised-cosine fringe 150 m (n = 0 only)", "fringe_xavg", "fringe",
                       None, fr_phis, 3000.0, dict(blank, shape="x-averaged", w_f_m=150.0))
    fr_rows += add_rows("x-averaged raised-cosine fringe 150 m (n = 0 only)", "fringe_xavg", "fringe",
                        None, [0.0, 4.0, 12.0, 24.0], 8000.0, dict(blank, shape="x-averaged", w_f_m=150.0))
    rows_csv += fr_rows
    for (shape, Lam, w_f, seed), d in sorted(cases.items()):
        if w_f != 150.0 or shape == "linear_sawtooth":
            continue
        E = scat_E(d)
        key = case_name(shape, Lam, w_f, seed)[5:]
        ex = dict(blank, shape=shape, Lam_m=Lam, w_f_m=w_f, seed=seed)
        rows_csv += add_rows(case_label(shape, Lam, w_f), key, "fringe", E, fr_phis, 3000.0, ex)
        rows_csv += add_rows(case_label(shape, Lam, w_f), key, "fringe", E, [0.0, 4.0, 12.0, 24.0], 8000.0, ex)

    # ---- joints
    jsc = load("joints_scat.pkl")
    if jsc is not None:
        for kind in JOINT_KINDS:
            for (L, g) in JOINT_LG:
                tag = f"joint_{kind}_L{L:g}_g{g:g}"
                E = smooth_y(jsc[(kind, L, g)]["E"])
                ex = dict(blank, joint_kind=kind, L_m=L, g_m=g)
                lab = f"joint {kind}, L_m {L:g} m, g {g:g} m"
                rows_csv += add_rows(lab, tag, tag, E, PHI_J, 3000.0, ex)
                rows_csv += add_rows(lab, tag, tag, E, [0.0], 8000.0, ex)
                rows_csv += add_rows(lab + " (n = 0 only)", tag + "_n0", tag, None, PHI_J, 3000.0, ex)
                rows_csv += add_rows(lab + " (n = 0 only)", tag + "_n0", tag, None, [0.0], 8000.0, ex)

    # ---- CSV
    cols = ["key", "label", "shape", "Lam_m", "w_f_m", "seed", "joint_kind", "L_m", "g_m", "R_m",
            "phi_deg", "y_g_m", "pldb_n0", "pldb_proxy_n0", "pldb", "pldb_proxy", "d_raw", "d_proxy",
            "d_headline", "d_headline_scat100", "d_headline_scat200", "pldb_scat_only",
            "pldb_proxy_scat_only"]
    with open(os.path.join(RESULTS, "alongtrack_results.csv"), "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols)
        w.writeheader()
        for r in rows_csv:
            w.writerow({k: (f"{r[k]:.4g}" if isinstance(r[k], float) else r[k]) for k in cols})

    # ---- summary tables
    by = {}
    for r in rows_csv:
        by.setdefault((r["key"], r["R_m"]), []).append(r)

    def summarize(key):
        r3 = sorted(by.get((key, 3000.0), []), key=lambda r: r["phi_deg"])
        r8 = {r["phi_deg"]: r for r in by.get((key, 8000.0), [])}
        if not r3:
            return None
        ut3 = next((r for r in r3 if r["phi_deg"] == 0.0), None)
        ut8 = r8.get(0.0)
        mx = max(r3, key=lambda r: r["d_headline"])
        prot = protected(r3)
        return dict(label=r3[0]["label"], ut3=ut3["d_headline"] if ut3 else np.nan,
                    ut8=ut8["d_headline"] if ut8 else np.nan, mx=mx["d_headline"], y_mx=mx["y_g_m"],
                    mx100=max(r["d_headline_scat100"] for r in r3),
                    mx200=max(r["d_headline_scat200"] for r in r3), prot=prot,
                    ut8_100=ut8["d_headline_scat100"] if ut8 else np.nan)

    print("\n| case | Delta undertrack R3 / R8 | Delta max over carpet (y) | max, scattered <= 200 / <= 100 Hz | protected (>= 10 dB) |")
    print("|---|---|---|---|---|")
    order = ["sharp", "fringe_xavg"] + [case_name(s, L, 150.0, 0)[5:] for s in ("sawtooth", "random", "rc_teeth")
                                        for L in SCORED_LAMS] + \
            [case_name("random", L, 150.0, 1)[5:] for L in (30.0, 300.0)]
    table = {}
    for key in order:
        s = summarize(key)
        if s is None:
            continue
        table[key] = s
        print(f"| {s['label']} | {s['ut3']:.1f} / {s['ut8']:.1f} | {s['mx']:.1f} ({s['y_mx'] / 1e3:.1f} km) | "
              f"{s['mx200']:.1f} / {s['mx100']:.1f} | {s['prot'] / 1e3:.1f} km |")
    print("\n| joint | Delta undertrack R3 / R8 | Delta max (y) | n = 0 only: undertrack R8 / max | max, scattered <= 100 Hz |")
    print("|---|---|---|---|---|")
    for kind in JOINT_KINDS:
        for (L, g) in JOINT_LG:
            tag = f"joint_{kind}_L{L:g}_g{g:g}"
            s, s0 = summarize(tag), summarize(tag + "_n0")
            if s is None:
                continue
            table[tag] = s
            table[tag + "_n0"] = s0
            print(f"| {kind}, L {L:g}, g {g:g} | {s['ut3']:.1f} / {s['ut8']:.1f} | {s['mx']:.1f} ({s['y_mx'] / 1e3:.1f} km) | "
                  f"{s0['ut8']:.1f} / {s0['mx']:.1f} | {s['mx100']:.1f} |")

    # ---- scattered-only maxima (incl. beyond the carpet edge)
    print("\n## Scattered orders alone (x-averaged), max PL over |y| <= 60 km, and where")
    print("| case | max PLdB raw / proxy (y) | beyond carpet edge |y| > %.1f km: max raw |" % (ycut / 1e3))
    print("|---|---|---|")
    so_tab = {}
    for (shape, Lam, w_f, seed), d in sorted(cases.items()):
        E = scat_E(d)
        ys = np.arange(0.0, 60001.0, 1000.0)
        pl = np.array([scat_only_pl(E_at(E, y)) for y in ys])
        i = int(np.argmax(pl[:, 0]))
        beyond = pl[ys > ycut, 0]
        so_tab[(shape, Lam, w_f, seed)] = (ys, pl)
        print(f"| {case_label(shape, Lam, w_f)} s{seed} | {pl[i, 0]:.1f} / {pl[i, 1]:.1f} ({ys[i] / 1e3:.0f} km) | "
              f"{np.max(beyond):.1f} |")
    if jsc is not None:
        for kind in JOINT_KINDS:
            for (L, g) in JOINT_LG:
                E = smooth_y(jsc[(kind, L, g)]["E"])
                ys = np.arange(0.0, 60001.0, 1000.0)
                pl = np.array([scat_only_pl(E_at(E, y)) for y in ys])
                i = int(np.argmax(pl[:, 0]))
                so_tab[(kind, L, g)] = (ys, pl)
                print(f"| joint {kind} L {L:g} g {g:g} | {pl[i, 0]:.1f} / {pl[i, 1]:.1f} ({ys[i] / 1e3:.0f} km) | "
                      f"{np.max(pl[ys > ycut, 0]):.1f} |")

    # ---- budgets
    budget_rows = budgets(cases, jsc)
    resteep = resteep_table(cases, jsc, spec, ref)
    v = load("valid.pkl")
    if v is not None:
        print_validation(v)
    figures(cases, jsc, rows_csv, so_tab, ycut, v)
    save("report_tables.pkl", dict(table=table, budget=budget_rows, resteep=resteep))
    return table


def protected(r3, thr=-10.0):
    """Contiguous half-width [m] from the track where the headline delta <= thr."""
    yy = [r["y_g_m"] for r in r3]
    dd = [r["d_headline"] for r in r3]
    w = 0.0
    for y, d in zip(yy, dd):
        if d <= thr:
            w = y
        else:
            break
    return w


def budgets(cases, jsc):
    """Aperture-energy budget per band: fractions of the scattered (n != 0) aperture energy
    that are propagating at the shield / reach the ground; ground-reaching scattered energy
    relative to the specular n = 0 fringe aperture energy and to a fully covered fringe's
    |D|^2. Writes results/alongtrack_budget.csv."""
    rows = []
    print("\n## Scattered-order budget (sum over band of aperture energy; fractions of n != 0 total)")
    print("| case | band Hz | scat / n0 (fringe) | propagating / scat | reaching ground / scat | ground / int|D|^2 [dB] |")
    print("|---|---|---|---|---|---|")
    items = [((s, L, w, sd), d) for (s, L, w, sd), d in sorted(cases.items())]
    if jsc is not None:
        for kind in JOINT_KINDS:
            for (L, g) in JOINT_LG:
                items.append(((f"joint_{kind}", L, g, 0), jsc[(kind, L, g)]))
    for key, d in items:
        f = F_S
        for lo, hi in BANDS_REPORT:
            sc = band_sum(f, d["ap_scat"], lo, hi)
            if sc <= 0:
                continue
            spn = band_sum(f, d["ap_spec"], lo, hi)
            pr = band_sum(f, d["ap_prop"], lo, hi)
            gr = band_sum(f, d["ap_ground"], lo, hi)
            dd = band_sum(f, d["ap_D"], lo, hi)
            r = dict(case=f"{key[0]}_L{key[1]:g}_{'w' if not str(key[0]).startswith('joint') else 'g'}{key[2]:g}_s{key[3]}",
                     band_lo=lo, band_hi=hi, scat_over_n0=sc / spn if spn > 0 else np.nan,
                     prop_frac=pr / sc, ground_frac=gr / sc,
                     ground_rel_D_dB=10 * np.log10(gr / dd) if gr > 0 else -np.inf)
            rows.append(r)
            if (key[0] in ("sawtooth", "random") and key[2] == 150.0 and key[3] == 0) or str(key[0]).startswith("joint_open"):
                print(f"| {r['case']} | {lo}-{hi} | {r['scat_over_n0']:.3g} | {r['prop_frac']:.3f} | "
                      f"{r['ground_frac']:.3f} | {r['ground_rel_D_dB']:.1f} |")
    cols = ["case", "band_lo", "band_hi", "scat_over_n0", "prop_frac", "ground_frac", "ground_rel_D_dB"]
    with open(os.path.join(RESULTS, "alongtrack_budget.csv"), "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols)
        w.writeheader()
        for r in rows:
            w.writerow({k: (f"{r[k]:.4g}" if isinstance(r[k], float) else r[k]) for k in cols})
    return rows


def resteep_table(cases, jsc, spec, ref):
    """Re-steepening check at the ground point/band where the scattered field is loudest:
    Sigma = int beta omega p /(rho c^3) ds along the dominant ray, from 1 km below the
    shield to the ground, with p from the WKB invariant and the band energy delivered
    within the 0.07 s ear window (conservative amplitude)."""
    print("\n## Re-steepening check of the scattered field (dominant band at its loudest ground point)")
    print("| case | y | band | p_eq [Pa] (free field) | Sigma (from 1 km below shield) |")
    print("|---|---|---|---|---|")
    out = []
    cfg = base_cfg()
    z_src = cfg.aircraft.altitude + cfg.z_mid
    fc = 10.0 ** (np.arange(42) / 10.0)
    for key, d in sorted(cases.items()):
        if key[2] != 150.0 or key[3] != 0 or key[1] not in (10.0, 30.0, 100.0, 300.0) or key[0] == "linear_sawtooth":
            continue
        E = scat_E(d)
        ys = np.arange(0.0, 60001.0, 1000.0)
        best = None
        for y in ys:
            Ef = E_at(E, y)
            eb = AT.band_energies_esd(F_S, 2 * Ef)
            ib = int(np.argmax(eb * fc ** 1.0))          # loudness-weighted dominant band (approx.)
            if best is None or eb[ib] > best[2]:
                best = (y, ib, eb[ib])
        y, ib, e = best
        m = (F_S >= LD_lo(ib)) & (F_S < LD_hi(ib))
        wgt = d["ap_ground"][m]
        if wgt.sum() <= 0:
            continue
        sx = float(np.nansum(d["sx_mean"][m] * wgt) / wgt.sum())
        sy = float(np.nansum(d["sy_mean"][m] * wgt) / wgt.sum())
        p_eq = math.sqrt(2.0) * math.sqrt(e / 0.07)
        try:
            S = AT.resteepening_check(None, sx, abs(sy), fc[ib], p_eq, z_src, z_src - 1000.0)
        except Exception as exc:   # noqa: BLE001
            S = float("nan")
            print("   resteep failed:", exc)
        out.append(dict(case=key, y=y, band=fc[ib], p_eq=p_eq / AT.REFLECTION, Sigma=S))
        print(f"| {case_label(key[0], key[1], key[2])} | {y / 1e3:.0f} km | {fc[ib]:.0f} Hz | "
              f"{p_eq / AT.REFLECTION:.3f} | {S:.3f} |")
    return out


def LD_lo(i):
    return 10.0 ** ((i - 0.5) / 10.0)


def LD_hi(i):
    return 10.0 ** ((i + 0.5) / 10.0)


def print_validation(v):
    print("\n## Validation")
    print(f"- PO difference field vs crossplane.po_solve: max rel. diff {v['po_match']:.1e}")
    for f, e in v["n0_vs_fringe_po"]:
        print(f"- n = 0 order of a sawtooth (exact c_0) vs crossplane fringe PO, f = {f:g} Hz: max rel. diff {e:.1e}")
    fv = v["flux_vs_pipeline"]
    print(f"- flux mapping vs linear propagate_from_circle (unshielded), dB at f = {fv['f']}:")
    for ph, yg, o in fv["rows"]:
        print(f"    phi {ph:g} (y {yg / 1e3:.1f} km): " + ", ".join(f"{x:+.2f}" for x in o))
    pb = v["po_bem"]
    f = pb["f"]
    print("\n| circle R 3 km | band | " + " | ".join(f"phi {p:g}" for p in pb["phis"] if p in (0, 12, 24, 39, 54)) + " |")
    print("|---|---|" + "---|" * 5)
    for case in ("sharp", "fringe"):
        for lo, hi in ((5, 10), (10, 20), (20, 50), (50, 100)):
            m = (f >= lo) & (f < hi)
            cells = []
            for j, p in enumerate(pb["phis"]):
                if p not in (0, 12, 24, 39, 54):
                    continue
                inc = np.sum(np.abs(pb["inc"][m, j]) ** 2)
                Lb = 10 * np.log10(np.sum(np.abs(pb["bem"][case][m, j]) ** 2) / inc)
                Lp = 10 * np.log10(np.sum(np.abs(pb["po"][case][m, j]) ** 2) / inc)
                cells.append(f"{Lb:.1f} / {Lp:.1f}")
            print(f"| {case} BEM / PO dB re inc | {lo}-{hi} | " + " | ".join(cells) + " |")


# ----------------------------------------------------------------------------- figures
SERIES = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4", "#008300", "#6250d6", "#e34948"]
INK, INK2, GRID = "#0b0b0b", "#52514e", "#e4e3df"


def _style(ax):
    ax.grid(True, color=GRID, lw=0.6)
    for sp_ in ("top", "right"):
        ax.spines[sp_].set_visible(False)
    for sp_ in ("left", "bottom"):
        ax.spines[sp_].set_color(INK2)
    ax.tick_params(colors=INK2, labelsize=8)


def figures(cases, jsc, rows_csv, so_tab, ycut, v):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    plt.rcParams.update({"font.size": 9, "axes.titlesize": 10, "axes.labelcolor": INK,
                         "text.color": INK, "figure.facecolor": "#fcfcfb", "axes.facecolor": "#fcfcfb"})
    by = {}
    for r in rows_csv:
        if r["R_m"] == 3000.0:
            by.setdefault(r["key"], []).append(r)

    def curve(key):
        rr = sorted(by.get(key, []), key=lambda r: r["y_g_m"])
        return np.array([r["y_g_m"] for r in rr]) / 1e3, np.array([r["d_headline"] for r in rr])

    # ---- 1: carpets
    fig, axs = plt.subplots(1, 3, figsize=(13.5, 4.2), sharey=True)
    for ax, shape, title in ((axs[0], "sawtooth", "Sawtooth serrations, w_f 150 m"),
                             (axs[1], "random", "Random-phase serrations, w_f 150 m")):
        for key, lab, col, ls in (("sharp", "sharp edge", INK2, "-"),
                                  ("fringe_xavg", "x-averaged fringe (n = 0 only)", INK, "--")):
            y, d = curve(key)
            ax.plot(y, np.maximum(d, -60), color=col, ls=ls, lw=1.6, label=lab)
        for k, L in enumerate([10.0, 30.0, 100.0, 300.0]):
            y, d = curve(case_name(shape, L, 150.0, 0)[5:])
            if y.size:
                ax.plot(y, np.maximum(d, -60), color=SERIES[k], lw=2, marker="o", ms=3,
                        label=f"Lambda = {L:g} m")
        ax.set_title(title, loc="left")
        ax.set_xlabel("ground lateral offset y_g [km]")
        _style(ax)
    axs[0].set_ylabel("Delta PLdB vs unshielded (headline)")
    axs[0].legend(frameon=False, fontsize=8, loc="lower right")
    ax = axs[2]
    y, d = curve("sharp")
    ax.plot(y, np.maximum(d, -60), color=INK2, lw=1.6, label="no joints (sharp edge)")
    k = 0
    for (L, g) in JOINT_LG:
        y, d = curve(f"joint_open_L{L:g}_g{g:g}")
        if y.size:
            ax.plot(y, d, color=SERIES[k], lw=2, marker="o", ms=3,
                    ls="-" if L == 250.0 else ":", label=f"open joint g {g:g} m / {L:g} m")
        k += 1
    ax.set_title("Open module joints (W 1 km, H2 30 m)", loc="left")
    ax.set_xlabel("ground lateral offset y_g [km]")
    ax.legend(frameon=False, fontsize=7.5, loc="lower right")
    _style(ax)
    axs[0].set_ylim(-62, 2)
    fig.text(0.01, 0.005, "Values below -60 dB (inaudible) are drawn at -60. Concorde M2, h 250 m, "
             "W 1 km, H2 30 m slab, R_m 3 km.", color=INK2, fontsize=7.5)
    fig.tight_layout(rect=(0, 0.03, 1, 1))
    fig.savefig(os.path.join(RESULTS, "alongtrack_carpets.png"), dpi=150)
    plt.close(fig)

    # ---- 2: orders: ground-reaching fraction vs f, and scattered-only PL vs y
    fig, axs = plt.subplots(1, 2, figsize=(11.5, 4.2))
    ax = axs[0]
    for k, L in enumerate(LAMS):
        d = cases.get(("sawtooth", L, 150.0, 0))
        if d is None:
            continue
        frac = np.where(d["ap_scat"] > 0, d["ap_ground"] / d["ap_scat"], 0.0)
        ax.plot(F_S, frac, color=SERIES[k], lw=1.6, label=f"Lambda = {L:g} m")
        V, cg = base_cfg().V, 340.29
        f_th = 1.0 / (L * (1 / V + 1 / cg))
        if f_th < 400:
            ax.axvline(f_th, color=SERIES[k], lw=0.8, ls=":")
    ax.set_xscale("log")
    ax.set_xlim(1, 400)
    ax.set_xlabel("frequency [Hz]")
    ax.set_ylabel("fraction of n != 0 aperture energy reaching the ground")
    ax.set_title("Sawtooth: scattered orders that reach the ground\n(dotted: first order can reach the ground, f = 1/(Lambda(1/V + 1/c_g)))",
                 loc="left", fontsize=9)
    ax.legend(frameon=False, fontsize=8)
    _style(ax)
    ax = axs[1]
    k = 0
    for L in (3.0, 10.0, 30.0, 100.0, 300.0):
        for shape, ls in (("sawtooth", "-"), ("random", "--")):
            t = so_tab.get((shape, L, 150.0, 0))
            if t is None:
                continue
            ys, pl = t
            ax.plot(ys / 1e3, np.maximum(pl[:, 0], 0), color=SERIES[LAMS.index(L)], ls=ls, lw=1.6,
                    label=f"{shape}, Lambda {L:g} m")
    ax.axvline(ycut / 1e3, color=INK2, lw=0.8)
    ax.text(ycut / 1e3 + 0.5, 2, "unshielded\ncarpet edge", color=INK2, fontsize=7.5)
    ax.set_xlabel("ground lateral offset y_g [km]")
    ax.set_ylabel("PLdB of the scattered orders alone (raw)")
    ax.set_title("Where the scattered orders land (x-averaged, both fringes)", loc="left", fontsize=9)
    ax.legend(frameon=False, fontsize=7, ncol=2)
    _style(ax)
    fig.tight_layout()
    fig.savefig(os.path.join(RESULTS, "alongtrack_orders.png"), dpi=150)
    plt.close(fig)

    # ---- 3: joints, open vs membranes
    if jsc is not None:
        fig, axs = plt.subplots(1, 2, figsize=(11.5, 4.0), sharey=True)
        for ax, kind in zip(axs, JOINT_KINDS):
            y, d = curve("sharp")
            ax.plot(y, d, color=INK2, lw=1.6, label="no joints")
            for k, (L, g) in enumerate(JOINT_LG):
                y, d = curve(f"joint_{kind}_L{L:g}_g{g:g}")
                y0, d0 = curve(f"joint_{kind}_L{L:g}_g{g:g}_n0")
                if y.size:
                    ax.plot(y, d, color=SERIES[k], lw=2, ls="-" if L == 250.0 else ":", marker="o", ms=3,
                            label=f"g {g:g} m every {L:g} m")
            ax.set_title(f"Joints: {'open gap' if kind == 'open' else 'face membranes only (air between)'}", loc="left")
            ax.set_xlabel("ground lateral offset y_g [km]")
            _style(ax)
        axs[0].set_ylabel("Delta PLdB vs unshielded (headline)")
        axs[0].legend(frameon=False, fontsize=7.5)
        fig.tight_layout()
        fig.savefig(os.path.join(RESULTS, "alongtrack_joints.png"), dpi=150)
        plt.close(fig)


# ----------------------------------------------------------------------------- main
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("stages", nargs="*", default=["valid", "scatter", "joints", "specular", "report"])
    ap.add_argument("--procs", type=int, default=PROCS)
    ap.add_argument("--only", default=None, help="specular: only tags starting with this")
    a = ap.parse_args()
    procs = min(a.procs, 2)
    t0 = time.time()
    if "valid" in a.stages:
        validate(procs)
        print(f"valid done {time.time() - t0:.0f} s", flush=True)
    if "scatter" in a.stages:
        incident_S()
        for c in sorted(SCAT_CASES, key=lambda c: (c[1] * (8 if c[0] == "random" else 1), c[0])):
            scatter_case(*c, procs=procs)
        print(f"scatter done {time.time() - t0:.0f} s", flush=True)
    if "hf" in a.stages:
        for c in HF_CASES:
            scatter_case(*c, procs=procs, hf=True)
        print(f"hf done {time.time() - t0:.0f} s", flush=True)
    if "joints" in a.stages:
        incident_S()
        joint_scatter(procs)
        joint_leak(procs)
        print(f"joints done {time.time() - t0:.0f} s", flush=True)
    if "specular" in a.stages:
        specular(procs, only=a.only)
        print(f"specular done {time.time() - t0:.0f} s", flush=True)
    if "report" in a.stages:
        report()


if __name__ == "__main__":
    main()
