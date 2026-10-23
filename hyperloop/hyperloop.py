#!/usr/bin/env python3
"""Hyperloop (our version): an elevated, ground-level low-pressure tube with compliant-gap maglev pods.

Screening-level feasibility map across disciplines (as with BoomShield): Mach 1-12 x setting (elevated on
pylons / at-grade / buried) x suspension (rigid EMS, passive SC-EDS gap, active large gap) x long-wave trim.

Run:   python3 -I hyperloop/hyperloop.py            all stages (terrain tiles cached in hyperloop/.cache)
       python3 -I hyperloop/hyperloop.py figs       figures only, from the CSVs
Out:   hl_terrain.csv        sampled ground profile along every spine leg and lateral offset (AWS Terrarium DEM)
       hl_alignment.csv      guideway error spectrum per setting vs allowable amplitude per suspension
       hl_align_limits.csv   highest Mach per setting x suspension x trim (fixed-guideway ride, 0.05 g)
       hl_route.csv          vertical alignment per leg x setting x Mach: tunnel / viaduct / elevated km, cost
       hl_speed.csv          terrain-following alternative: local speed caps and leg times
       hl_grid.csv           Mach x setting x suspension x trim x pressure x blockage: margins per discipline
                             (hl_maps_grid.csv is the > 2 MB variant, gitignored)
       hl_map.csv            best point per Mach x setting x suspension x trim (status, binding)
       hl_cost.csv           bottom-up $M/km per tube pair, low / base / high, per element and setting
       hl_econ.csv           door to door, capacity, B/C (gross / net, 30 / 100 %, 2045-2060, vs sub / Boomless /
                             BoomShield M2.2), break-even $M/km
       hl_levers.csv         what moves the frontier
       hl_maps.png, hl_route.png, hl_econ.png

Reuses read-only (imported, never modified on disk): regime/regime.py (evaluate: drag, Kantrowitz, pod heating,
vacuum, structure, breach, capacity), hybrid/hybrid.py (allow_amp, TOL_CASES, suspension transfer),
heating/heating.py (gas properties for the per-pass tube-gas check), physics/skyloop_physics.py (Kantrowitz,
wall), system/sky_model.py (Net, evaluate, ride_min, cap_pair_yr, lifecycle, costs, comparators),
flight/flight.py (imported through hybrid).  In-process overrides of their parameter dicts are restored after use.
Mach convention (repo): V = M x 295.1 m/s.  (A) = assumption; web items are unverified search summaries.
"""
from __future__ import annotations

import csv
import io
import math
import os
import sys
import urllib.request

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
REV = os.path.join(ROOT, "skyloop-review")
CACHE = os.path.join(HERE, ".cache")
os.makedirs(CACHE, exist_ok=True)
sys.dont_write_bytecode = True
for _d in ("regime", "hybrid", "flight", "heating", "stability", "physics", "system"):
    sys.path.insert(0, os.path.join(REV, _d))

G = 9.81
R_E = 6371.0e3


def write_csv(name, rows):
    path = os.path.join(HERE, name)
    keys = []
    for r in rows:
        keys += [k for k in r if k not in keys]
    with open(path, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=keys)
        w.writeheader()
        for r in rows:
            w.writerow({k: (f"{v:.4g}" if isinstance(v, (float, np.floating)) else v) for k, v in r.items()})
    print(f"wrote {name} ({len(rows)} rows, {os.path.getsize(path) / 1e6:.2f} MB)", flush=True)


def read_csv(name):
    path = os.path.join(HERE, name)
    if not os.path.exists(path):
        return []
    out = []
    for r in csv.DictReader(open(path)):
        rr = {}
        for k, v in r.items():
            try:
                rr[k] = float(v)
            except (TypeError, ValueError):
                rr[k] = v
        out.append(rr)
    return out


# ============================================================================= 0. trunk and terrain
STN = {"LGA": (40.777, -73.874), "CLE": (41.412, -81.850), "ORD": (41.974, -87.907), "OMA": (41.303, -95.894),
       "DEN": (39.856, -104.674), "LAS": (36.084, -115.154), "LAX": (33.942, -118.409)}   # sky_model df (airports)
LEGS = [("LGA", "CLE"), ("CLE", "ORD"), ("ORD", "OMA"), ("OMA", "DEN"), ("DEN", "LAS"), ("LAS", "LAX")]
OFFSETS = [-50e3, -25e3, 0.0, 25e3, 50e3]     # m, lateral mid-leg offset of a circular-arc alignment
DX = 500.0                                    # m, profile sampling
L_URB = 20e3                                  # m, urban approach at each end (system.md convention: tunnel/viaduct)
ZOOM = 9                                      # Terrarium z9: ~235 m pixels at 40 N
TILE_URL = "https://s3.amazonaws.com/elevation-tiles-prod/terrarium/{z}/{x}/{y}.png"   # AWS open data (Mapzen)


def _xyz(lat, lon):
    la, lo = np.radians(lat), np.radians(lon)
    return np.stack([np.cos(la) * np.cos(lo), np.cos(la) * np.sin(lo), np.sin(la)], -1)


def gc_dist(p, q):
    a, b = _xyz(*p), _xyz(*q)
    return R_E * 2 * math.asin(min(1.0, float(np.linalg.norm(a - b)) / 2))


def leg_path(a, b, off):
    """Points every DX along a great circle a->b displaced laterally by off*4f(1-f) (an arc of R = L^2/8|off|)."""
    p, q = STN[a], STN[b]
    A, B = _xyz(*p), _xyz(*q)
    om = math.acos(float(np.clip(A @ B, -1, 1)))
    L = R_E * om
    n = int(L // DX) + 1
    f = np.linspace(0, 1, n)
    P = (np.sin((1 - f) * om)[:, None] * A + np.sin(f * om)[:, None] * B) / math.sin(om)
    nrm = np.cross(A, B)
    nrm /= np.linalg.norm(nrm)                         # unit normal of the great-circle plane (left of travel)
    y = off * 4 * f * (1 - f) / R_E                    # angular displacement
    P = P * np.cos(y)[:, None] + nrm[None, :] * np.sin(y)[:, None]
    lat = np.degrees(np.arcsin(np.clip(P[:, 2], -1, 1)))
    lon = np.degrees(np.arctan2(P[:, 1], P[:, 0]))
    s = np.concatenate([[0], np.cumsum(R_E * np.arccos(np.clip(np.sum(P[1:] * P[:-1], 1), -1, 1)))])
    return lat, lon, s, L


def _tile_xy(lat, lon, z):
    n = 2 ** z
    x = (lon + 180.0) / 360.0 * n
    la = np.radians(lat)
    y = (1.0 - np.log(np.tan(la) + 1 / np.cos(la)) / math.pi) / 2.0 * n
    return x, y


_TILES = {}


def _tile(z, tx, ty):
    key = (z, tx, ty)
    if key in _TILES:
        return _TILES[key]
    fn = os.path.join(CACHE, f"t{z}_{tx}_{ty}.png")
    if not os.path.exists(fn):
        req = urllib.request.Request(TILE_URL.format(z=z, x=tx, y=ty), headers={"User-Agent": "hyperloop-study"})
        for k in range(8):
            try:
                data = urllib.request.urlopen(req, timeout=60).read()
                break
            except Exception as e:                     # proxy resets: retry with back-off
                print(f"    tile {key} retry {k}: {e}", flush=True)
                import time
                time.sleep(2 + 3 * k)
        else:
            raise RuntimeError(f"tile {key} failed")
        with open(fn + ".tmp", "wb") as fh:
            fh.write(data)
        os.replace(fn + ".tmp", fn)
    import matplotlib.image as mpimg
    im = mpimg.imread(fn)                              # float 0-1, RGB(A)
    rgb = np.round(im[..., :3] * 255.0)
    elev = rgb[..., 0] * 256 + rgb[..., 1] + rgb[..., 2] / 256 - 32768.0
    _TILES[key] = elev
    return elev


def elevation(lat, lon, z=ZOOM):
    x, y = _tile_xy(lat, lon, z)
    out = np.empty(len(lat))
    tx, ty = np.floor(x).astype(int), np.floor(y).astype(int)
    for (a, b) in sorted(set(zip(tx, ty))):
        m = (tx == a) & (ty == b)
        e = _tile(z, a, b)
        px = np.clip(((x[m] - a) * 256).astype(int), 0, 255)
        py = np.clip(((y[m] - b) * 256).astype(int), 0, 255)
        out[m] = e[py, px]
    return out


def stage_terrain():
    """Ground profile along each spine leg for five lateral offsets (cached CSV; tiles cached in .cache)."""
    if os.path.exists(os.path.join(HERE, "hl_terrain.csv")):
        return load_terrain()
    rows = []
    for a, b in LEGS:
        for off in OFFSETS:
            lat, lon, s, L = leg_path(a, b, off)
            z = elevation(lat, lon)
            print(f"  terrain {a}-{b} offset {off / 1e3:+.0f} km: {s[-1] / 1e3:.0f} km, "
                  f"z {z.min():.0f}-{z.max():.0f} m, tiles cached {len(_TILES)}", flush=True)
            for i in range(len(s)):
                rows.append(dict(leg=f"{a}-{b}", off_km=off / 1e3, s_km=s[i] / 1e3, lat=lat[i], lon=lon[i], z_m=z[i]))
    write_csv("hl_terrain.csv", rows)
    return load_terrain()


def load_terrain():
    T = {}
    for r in read_csv("hl_terrain.csv"):
        T.setdefault((r["leg"], r["off_km"]), []).append((r["s_km"], r["z_m"], r["lat"], r["lon"]))
    return {k: np.array(v) for k, v in T.items()}




# ----------------------------------------------------------------------------- reused modules (read-only)
import hybrid as HY  # noqa: E402  (imports flight, heating, stability, regime, physics, sky_model)
R, sm, H, S = HY.R, HY.sm, HY.H, HY.S
C_REF = R.C_REF                                   # 295.1 m/s
SETTINGS = ["elevated", "at-grade", "buried"]
M_GRID = np.round(np.arange(1.0, 12.001, 0.25), 2)

# ============================================================================= 1. alignment and ride
# Guideway error as a sinusoid amplitude A(lambda) per source; sources add as RSS (A).  Each band rolls off as
# (lambda/lo)^2 below and (hi/lambda) above.  'trim' = slow active long-wave trim (jacks at each pylon / slab
# adjusters driven by a survey every few weeks): it removes settlement, creep and ground movement, not the daily
# thermal or wind parts.  Seismic is a stop-and-resurvey event, not a ride input (section 6).
SPAN = 30.0                                       # m, pylon spacing (Transrapid 24.8-31 m girders)
ERR = {  # (name, A [m] untrimmed, A trimmed, lam_lo, lam_hi) or ('slope', per m untrimmed, trimmed)
    "elevated": [
        ("survey + as-built", 1.0e-3, 1.0e-3, 5.0, 1e9),
        ("survey drift (long wave)", ("slope", 2e-6, 2e-6)),
        ("thermal camber, 10 K top-bottom, 30 m span", 3.0e-3, 3.0e-3, 25.0, 70.0),
        ("pylon differential settlement (Shanghai re-align at 5 mm)", 5.0e-3, 1.0e-3, 40.0, 300.0),
        ("wind on pylons, 25 m/s operating limit", 2.0e-3, 2.0e-3, 60.0, 2000.0),
        ("creep / long-wave ground movement", ("slope", 15e-6, 3e-6)),
    ],
    "at-grade": [
        ("survey + as-built", 1.0e-3, 1.0e-3, 5.0, 1e9),
        ("survey drift (long wave)", ("slope", 2e-6, 2e-6)),
        ("embankment settlement (slab-track limit 15 mm)", 8.0e-3, 2.0e-3, 20.0, 200.0),
        ("frost heave (Midwest frost depth 1-1.5 m)", 5.0e-3, 2.0e-3, 10.0, 100.0),
        ("long-wave ground movement", ("slope", 20e-6, 5e-6)),
    ],
    "buried": [
        ("survey + as-built", 1.0e-3, 1.0e-3, 5.0, 1e9),
        ("survey drift (long wave)", ("slope", 2e-6, 2e-6)),
        ("lining settlement / convergence", 2.0e-3, 1.0e-3, 20.0, 500.0),
        ("long-wave ground movement", ("slope", 5e-6, 2e-6)),
    ],
}
SUSP = [  # hybrid.TOL_CASES format: name, primary (f_n, zeta, kind, usable gap) or ('eds', f, zeta, kind), secondary
    ("rigid EMS 10 mm + 1 Hz secondary (Transrapid-like)",) + tuple(HY.TOL_CASES[0][1:]),
    ("passive SC-EDS gap 0.10 m (SCMaglev-like)",) + tuple(HY.TOL_CASES[1][1:]),
    ("active gap 0.5 Hz, 0.15 m", (0.5, 0.7, "skyhook", 0.15), None),
    ("active large gap 0.3 Hz, 0.3 m",) + tuple(HY.TOL_CASES[4][1:]),
]
SUSP_SHORT = ["rigid EMS", "passive EDS", "active 0.5 Hz", "active 0.3 Hz"]
LAMS = np.geomspace(10.0, 1e5, 61)


def _band(lam, A, lo, hi):
    if lam < lo:
        return A * (lam / lo) ** 2
    if lam > hi:
        return A * hi / lam
    return A


def err_amp(lam, setting, trim, parts=False):
    out = {}
    for e in ERR[setting]:
        if isinstance(e[1], tuple):
            out[e[0]] = (e[1][2] if trim else e[1][1]) * lam
        else:
            out[e[0]] = _band(lam, e[2] if trim else e[1], e[3], e[4])
    tot = math.sqrt(sum(v * v for v in out.values()))
    return (tot, out) if parts else tot


def align_ratio(M, setting, trim, case):
    """min over wavelength of allowable / achievable amplitude (>= 1 pass, >= 0.5 marginal) and the binding lambda."""
    V = M * C_REF
    best = (1e9, None)
    for lam in LAMS:
        r = HY.allow_amp(V, lam, case) / err_amp(lam, setting, trim)
        if r < best[0]:
            best = (r, lam)
    return best


def stage_alignment():
    rows = []
    for st in SETTINGS:
        for trim in (0, 1):
            for lam in LAMS:
                tot, parts = err_amp(lam, st, trim, parts=True)
                r = dict(kind="error spectrum", setting=st, trim=trim, lambda_m=lam, A_total_m=tot)
                r.update({f"A[{k}]": v for k, v in parts.items()})
                rows.append(r)
    for M in (2.0, 4.0, 6.0, 8.0, 10.0, 12.0):
        for lam in LAMS:
            r = dict(kind="allowable", Mach=M, lambda_m=lam)
            for c, sh in zip(SUSP, SUSP_SHORT):
                r[sh] = HY.allow_amp(M * C_REF, lam, c)
            rows.append(r)
    write_csv("hl_alignment.csv", rows)
    lim = []
    for st in SETTINGS:
        for trim in (0, 1):
            for c, sh in zip(SUSP, SUSP_SHORT):
                mp = mm = 0.0
                lam_b = None
                for M in np.arange(0.2, 20.01, 0.1):
                    rr, lb = align_ratio(M, st, trim, c)
                    if rr >= 1:
                        mp, lam_b = M, lb
                    if rr >= 0.5:
                        mm = M
                lim.append(dict(setting=st, trim=trim, suspension=sh, Mach_pass=round(mp, 1), Mach_marginal=round(mm, 1),
                                binding_lambda_m=lam_b if lam_b else float("nan")))
                print(f"  alignment {st:9s} trim {trim} {sh:14s}: pass M{mp:.1f}, marginal M{mm:.1f}", flush=True)
    write_csv("hl_align_limits.csv", lim)
    return lim


# ============================================================================= 2. route geometry
A_V = 0.05 * G            # vertical-curve acceleration: HSR min vertical radii 10-24 km at 300-350 km/h ~ 0.05 g
A_V_LEVER = 0.10 * G      # seated, belted pods (lever)
A_LAT = {"rail-like (0.1 g felt + 8 deg cant)": 0.24 * G, "banked tube (regime base)": 0.5 * G,
         "rolling pod (0.15 g felt + 25 deg)": 0.62 * G}
A_LAT_BASE = 0.5 * G
GRADE_MAX = 0.05          # (A) LSM-driven pods climb 5 % easily (Transrapid 10 %)
H_OBST, L_OBST = 2.0e3, 50e3   # dodge a 2 km obstacle (town, lake, reserve) every ~50 km of route (A)
BAND = {  # guideway height above ground that the setting builds without a special structure: (h_lo, h_hi), m
    "elevated": (6.0, 25.0), "at-grade": (-2.0, 4.0)}
TUNNEL_BELOW = {"elevated": -10.0, "at-grade": -12.0}   # deeper than this -> bored / mined tunnel
M_ROUTE = [0.3, 0.6, 1.0, 1.5, 2.0, 3.0, 4.0, 5.0, 6.0, 8.0, 10.0, 12.0]   # 0.3 ~ HSR 320 km/h


def r_min_h(M, a=A_LAT_BASE):
    return (M * C_REF) ** 2 / a


def dodge_km(M, a=A_LAT_BASE):
    """Reverse-curve length to side-step an obstacle of H_OBST at R_min (two arcs)."""
    return 2 * math.sqrt(r_min_h(M, a) * H_OBST) / 1e3


def profile_lp(t, kappa, setting):
    """Vertical alignment z over ground t (DX spacing) with |z''| <= kappa and |z'| <= GRADE_MAX.
    elevated / at-grade: min sum 3*below-band + 1*above-band + 5*above 100 m (prefer viaduct to tunnel);
    buried: z <= t - 30 m cover, minimise total depth."""
    from scipy.optimize import linprog
    from scipy.sparse import diags, eye, hstack, vstack, csr_matrix
    n = len(t)
    D2 = diags([np.ones(n - 2), -2 * np.ones(n - 2), np.ones(n - 2)], [0, 1, 2], shape=(n - 2, n))
    D1 = diags([-np.ones(n - 1), np.ones(n - 1)], [0, 1], shape=(n - 1, n))
    c2, c1 = kappa * DX * DX, GRADE_MAX * DX
    if setting == "buried":
        A = vstack([D2, -D2, D1, -D1]).tocsr()
        b = np.concatenate([np.full(n - 2, c2), np.full(n - 2, c2), np.full(n - 1, c1), np.full(n - 1, c1)])
        res = linprog(-np.ones(n), A_ub=A, b_ub=b, bounds=[(None, ti - 30.0) for ti in t], method="highs")
        return res.x if res.status == 0 else None
    h_lo, h_hi = BAND[setting]
    I = eye(n, format="csr")
    Z = csr_matrix((n, n))
    Zs2, Zs1 = csr_matrix((n - 2, n)), csr_matrix((n - 1, n))
    A = vstack([hstack([D2, Zs2, Zs2, Zs2]), hstack([-D2, Zs2, Zs2, Zs2]),
                hstack([D1, Zs1, Zs1, Zs1]), hstack([-D1, Zs1, Zs1, Zs1]),
                hstack([-I, -I, Z, Z]),            # u >= t + h_lo - z
                hstack([I, Z, -I, Z]),             # v >= z - t - h_hi
                hstack([I, Z, Z, -I])]).tocsr()    # w >= z - t - 100
    b = np.concatenate([np.full(n - 2, c2), np.full(n - 2, c2), np.full(n - 1, c1), np.full(n - 1, c1),
                        -(t + h_lo), t + h_hi, t + 100.0])
    c = np.concatenate([np.zeros(n), 3 * np.ones(n), np.ones(n), 5 * np.ones(n)])
    bnds = [(None, None)] * n + [(0, None)] * (3 * n)
    res = linprog(c, A_ub=A, b_ub=b, bounds=bnds, method="highs")
    return res.x[:n] if res.status == 0 else None


def classify(z, t, setting):
    """km by structure class and the height-weighted viaduct factor."""
    h = z - t
    km = DX / 1e3
    out = dict(tunnel_km=0.0, deep_tunnel_km=0.0, cut_km=0.0, std_km=0.0, tall_km=0.0, tall_factor_km=0.0,
               max_h=float(h.max()), max_depth=float(-h.min()))
    if setting == "buried":
        out["tunnel_km"] = len(h) * km
        out["deep_tunnel_km"] = float(np.sum(h < -800)) * km
        return out
    tb = TUNNEL_BELOW[setting]
    h_lo, h_hi = BAND[setting]
    out["tunnel_km"] = float(np.sum(h < tb)) * km
    out["deep_tunnel_km"] = float(np.sum(h < -800)) * km
    out["cut_km"] = float(np.sum((h >= tb) & (h < h_lo - 2))) * km
    out["std_km"] = float(np.sum((h >= h_lo - 2) & (h <= h_hi))) * km
    tall = h > h_hi
    out["tall_km"] = float(np.sum(tall)) * km
    out["tall_factor_km"] = float(np.sum(np.minimum(1 + (h[tall] - h_hi) / 30.0, 8.0))) * km   # +1x per 30 m (A)
    return out


def stage_route(terrain, a_v=A_V, tag="base"):
    """Vertical alignment per leg x lateral offset x setting x Mach; best offset with R_arc >= R_min."""
    fn = f"hl_route_{tag}.csv" if tag != "base" else "hl_route.csv"
    done = {(r["leg"], r["off_km"], r["setting"], r["Mach"]): r for r in read_csv(fn)}
    rows = list(done.values())
    hfn = os.path.join(CACHE, f"route_h_{tag}.npz")
    H0 = dict(np.load(hfn)) if os.path.exists(hfn) else {}
    for (leg, off), arr in terrain.items():
        s, t = arr[:, 0] * 1e3, arr[:, 1]
        L = s[-1]
        keep = (s >= L_URB) & (s <= L - L_URB)
        t_in = np.maximum(t[keep], 0.0)               # below sea level only at the coast (water -> bridge cost)
        R_arc = L * L / (8 * abs(off * 1e3)) if off else 1e12
        for st in SETTINGS:
            for M in M_ROUTE:
                hkey = f"{leg}|{st}|{M:g}"
                if (leg, off, st, M) in done and (off != 0 or st == "buried" or hkey in H0):
                    continue
                V = M * C_REF
                z = profile_lp(t_in, a_v / V ** 2, st)
                if z is None:
                    print(f"  LP failed {leg} {off} {st} M{M}", flush=True)
                    continue
                c = classify(z, t_in, st)
                if off == 0 and st != "buried":
                    H0[hkey] = (z - t_in).astype(np.float32)
                r = dict(leg=leg, off_km=off, setting=st, Mach=M, a_v_g=a_v / G, L_km=L / 1e3,
                         rural_km=float(keep.sum()) * DX / 1e3, R_arc_km=R_arc / 1e3, **c)
                rows = [x for x in rows if (x["leg"], x["off_km"], x["setting"], x["Mach"]) != (leg, off, st, M)]
                rows.append(r)
                done[(leg, off, st, M)] = r
                if st == "elevated" and off == 0:
                    print(f"  route {leg} {st} M{M:g}: tunnel {c['tunnel_km']:.0f} km, tall {c['tall_km']:.0f} km, "
                          f"max depth {c['max_depth']:.0f} m", flush=True)
            write_csv(fn, rows)
            np.savez(hfn, **H0)
    return rows, H0


def route_summary(rows, a_v_tag="base"):
    """Spine totals per setting x Mach for the full-speed strategy: per leg, cheapest offset with R_arc >= R_min.
    Cheapness proxy: tunnel km x 7 + tall-factor km + std km (relative structure cost, base)."""
    out = {}
    by = {}
    for r in rows:
        by.setdefault((r["setting"], r["Mach"], r["leg"]), []).append(r)
    for (st, M, leg), rs in by.items():
        ok = [r for r in rs if r["R_arc_km"] * 1e3 >= r_min_h(M)] or [r for r in rs if r["off_km"] == 0]
        if st == "buried":
            best = min(ok, key=lambda r: r["deep_tunnel_km"] + r["max_depth"] / 100)
        else:
            best = min(ok, key=lambda r: 7 * r["tunnel_km"] + r["tall_factor_km"] + r["std_km"] + 1.5 * r["cut_km"])
        out.setdefault((st, M), {})[leg] = best
    summ = {}
    for (st, M), legs in out.items():
        tot = {k: sum(l[k] for l in legs.values()) for k in ("rural_km", "tunnel_km", "deep_tunnel_km", "cut_km",
                                                               "std_km", "tall_km", "tall_factor_km", "L_km")}
        tot["max_depth"] = max(l["max_depth"] for l in legs.values())
        tot["max_h"] = max(l["max_h"] for l in legs.values())
        tot["tunnel_share"] = tot["tunnel_km"] / tot["rural_km"]
        tot["legs"] = legs
        summ[(st, M)] = tot
    return summ


def route_at(summ, st, M, key):
    Ms = sorted({m for (s, m) in summ if s == st})
    return float(np.interp(M, Ms, [summ[(st, m)][key] for m in Ms]))


# ============================================================================= 3-6. physics grid and map
P_GRID = np.geomspace(1.0, 1000.0, 13)
PODS = [8000.0, 17000.0, 30000.0]                 # 20 / 50 / 93 seats (regime.pod, 300 kg/seat)
REF_POD = 17000.0
RULE = "moving block 10 s"
DEMAND = R.DEMAND                                  # 81 M pax/yr on the spine
PAX_H_DIR = DEMAND / 365 * sm.P["peak_share"] / 2  # 9,250 pax/h per direction in the peak
K_GAS = H.k_gas(288.0)
CV = 1005.0 / 1.4
DT_GAS = 50.0                                      # K, steady gas rise allowed (heating.py)
ETA = sm.P["eta_drive"]
P_STN_OK, P_STN_MARG = 1e9, 5e9                    # W, station launch power per direction (hybrid convention)
L_NONSTOP = 4074e3                                 # m, NYC-LAX along the spine alignment
CHECKS = ["alignment", "route", "drag_energy", "heating", "gas_heat", "vacuum", "structure", "safety", "capacity",
          "power", "trip"]
LABEL = {"alignment": "alignment / ride", "route": "route (terrain, curves)", "drag_energy": "energy / choking",
         "heating": "pod heating", "gas_heat": "tube-gas heat", "vacuum": "vacuum (leaks)", "structure": "structure",
         "safety": "breach deceleration", "capacity": "capacity", "power": "station launch power",
         "trip": "NYC-LAX trip length"}


def st_of(m, tol):
    return "pass" if m >= 0 else ("marginal" if m >= -tol else "fail")


def physics_point(setting, M, p, m_pod, beta, F=None):
    """Margins (log10 limit/demand) and tolerances for every check except alignment and route."""
    F = F or {}
    r = R.evaluate("buried" if setting == "buried" else "elevated", "maglev", M, p, m_pod, beta, F.get("rule", RULE))
    mg, tol = {}, {}
    for k in ("drag_energy", "heating", "vacuum", "structure", "capacity"):
        mg[k], tol[k] = r[f"m_{k}"], (math.log10(2) if k != "capacity" else math.log10(3))
    tol["drag_energy"] = math.log10(R.E_JET / R.E_LIM)
    tol["vacuum"] = math.log10(5)
    if setting == "buried":
        mg["structure"] = 0.5
    # breach: regime's full-bore plug; lever divides the air density met (fast valves, pods braking first)
    a_plug = r["a_plug_g"] / F.get("breach_rho_div", 1.0)
    mg["safety"], tol["safety"] = math.log10(R.A_PLUG_OK / a_plug), math.log10(R.A_PLUG_MAX / R.A_PLUG_OK)
    # tube gas: per-pass rise and gas->wall throughput (heating.tube_side physics at 288 K, Nu_D = 8)
    pd_ = R.pod(m_pod)
    V = M * C_REF
    rho = p / (R.R_AIR * 288.0)
    A_t = math.pi * r["D_tube"] ** 2 / 4
    drag = r["drag_par_N"]
    dT_pass = drag / (rho * A_t * CV)
    pods_h = PAX_H_DIR / (pd_["n"] * sm.P["lf"])
    N_max = 3600 * 8 * math.pi * K_GAS * DT_GAS * F.get("f_mix", 1.0) / drag
    m_pass = math.log10(H.BASE["dT_pass_ok"] / dT_pass)
    m_thr = math.log10(N_max / pods_h)
    mg["gas_heat"] = min(m_pass, m_thr)
    tol["gas_heat"] = math.log10(H.BASE["dT_pass_marg"] / H.BASE["dT_pass_ok"]) if m_pass < m_thr else math.log10(3)
    # station launch power, gross, one direction, at peak demand (pods depart every 3600/pods_h s)
    KE = 0.5 * m_pod * V * V
    P_stn = KE * min(pods_h, 3600 / sm.headway_s(M, F.get("rule", RULE))) / 3600 / ETA
    P_lim = F.get("P_stn", P_STN_OK)
    mg["power"], tol["power"] = math.log10(P_lim / P_stn), math.log10(P_STN_MARG / P_STN_OK)
    # trip: a nonstop NYC-LAX run must fit launch + brake at 0.3 g
    mg["trip"], tol["trip"] = math.log10(L_NONSTOP * sm.P["a_acc"] / V ** 2), 1e-6
    out = dict(setting=setting, Mach=M, p_Pa=p, pod_t=m_pod / 1e3, seats=pd_["n"], beta=beta, D_tube=r["D_tube"],
               choked=r["choked"], beta_kantrowitz=r["beta_kantrowitz"], CD=r["CD"], drag_N=drag,
               E_Wh_paxkm=r["E_Wh_paxkm"], crossover_km=r["crossover_km"], T_nose=r["T_nose"],
               Q_in_kW=r["Q_in_kW"], Q_out_kW=r["Q_out_kW"], S_req_m3s_km=r["S_req_m3s_km"],
               pump_kW_km=r["pump_kW_km"], reevac_h_10km=r["reevac_h_10km"], wall_kg_m=r["wall_kg_m"],
               a_plug_g=a_plug, KE_tTNT=r["KE_tTNT"], pods_within_stop=r["pods_within_stop"],
               cap_pair_Mpax=r["cap_pair_Mpax_yr"], pods_h_peak=pods_h, gas_dT_pass_K=dT_pass, gas_Nmax_h=N_max,
               P_station_GW=P_stn / 1e9, P_station_net_GW=P_stn / 1e9 * (1 - ETA * ETA * sm.P["eta_regen"]),
               L_acc_km=V * V / (2 * sm.P["a_acc"]) / 1e3)
    return out, mg, tol


def route_margin(summ, setting, M, F=None):
    F = F or {}
    a_lat = F.get("a_lat", A_LAT_BASE)
    m_h = math.log10(L_OBST / 1e3 / dodge_km(M, a_lat))
    t_h = math.log10(3)
    if setting == "buried":
        dmax = route_at(summ, "buried", M, "max_depth")
        m_v, t_v = math.log10(1200.0 / dmax), math.log10(2300.0 / 1200.0)      # Gotthard ~2,300 m overburden
        info = dict(max_overburden_m=dmax, tunnel_share=1.0)
    else:
        sh = route_at(summ, setting, M, "tunnel_share")
        m_v, t_v = math.log10(0.15 / max(sh, 1e-3)), math.log10(0.35 / 0.15)
        info = dict(tunnel_share=sh, max_overburden_m=route_at(summ, setting, M, "max_depth"))
    info["dodge_km"] = dodge_km(M, a_lat)
    info["R_min_km"] = r_min_h(M, a_lat) / 1e3
    if m_h / t_h < m_v / t_v:
        return m_h, t_h, info
    return m_v, t_v, info


def combine(mg, tol, checks=CHECKS):
    sts = {k: st_of(mg[k], tol[k]) for k in checks}
    s = "fail" if "fail" in sts.values() else ("marginal" if "marginal" in sts.values() else "pass")
    z = {k: (mg[k] / tol[k] if tol[k] > 0 else (1e3 if mg[k] >= 0 else -1e3)) for k in checks}
    bind = min(z, key=z.get)
    return s, bind, sts, min(z.values())


RANKS = {"pass": 2, "marginal": 1, "fail": 0}


def stage_grid(summ, F=None, write=True, settings=SETTINGS, m_grid=M_GRID, pods=PODS, tag=""):
    """Best (p, beta, pod) per setting x suspension x trim x Mach; full grid written to hl_grid.csv (gitignored)."""
    F = F or {}
    lims = {}
    for st in settings:
        for trim in (0, 1):
            for c, sh in zip(SUSP, SUSP_SHORT):
                lims[(st, trim, sh)] = c
    grid, best = [], {}
    align_cache = {}
    for st in settings:
        for M in m_grid:
            m_r, t_r, rinfo = route_margin(summ, st, M, F)
            for m_pod in pods:
                for b in R.BETAS:
                    if R.pod(m_pod)["d"] / math.sqrt(b) > R.D_MAX:
                        continue
                    for p in P_GRID:
                        out, mg, tol = physics_point(st, M, p, m_pod, b, F)
                        mg["route"], tol["route"] = m_r, t_r
                        out.update(rinfo)
                        for trim in (0, 1):
                            for sh in SUSP_SHORT:
                                key = (st, trim, sh, M)
                                if key not in align_cache:
                                    a, lam = align_ratio(M, st, trim, lims[(st, trim, sh)])
                                    align_cache[key] = (math.log10(a), lam)
                                mg["alignment"], tol["alignment"] = align_cache[key][0], math.log10(2)
                                s, bind, sts, zmin = combine(mg, tol)
                                rank = (RANKS[s], -list(sts.values()).count("fail"), zmin, -out["E_Wh_paxkm"])
                                if key not in best or rank > best[key][0]:
                                    best[key] = (rank, dict(out, trim=trim, suspension=sh, status=s, binding=bind,
                                                            **{f"m_{k}": mg[k] for k in CHECKS},
                                                            **{f"s_{k}": sts[k] for k in CHECKS}))
                                if write and trim == 1 and sh == "active 0.3 Hz" and m_pod == REF_POD:
                                    grid.append(dict(out, status=s, binding=bind, **{f"m_{k}": mg[k] for k in CHECKS}))
        print(f"  grid {tag}{st} done", flush=True)
    rows = [v[1] for v in best.values()]
    if write:
        write_csv("hl_grid.csv", grid)
        write_csv("hl_map.csv", rows)
    return rows


def frontier(rows, key="status", ok=("pass",)):
    """Highest Mach with status in ok, per setting x trim x suspension, and the binding check just above it."""
    out = {}
    by = {}
    for r in rows:
        by.setdefault((r["setting"], r["trim"], r["suspension"]), []).append(r)
    for k, rs in by.items():
        rs = sorted(rs, key=lambda r: r["Mach"])
        good = [r for r in rs if r[key] in ok]
        mx = max((r["Mach"] for r in good), default=float("nan"))
        above = [r for r in rs if r["Mach"] > (mx if mx == mx else 0)]
        out[k] = (mx, above[0]["binding"] if above else "grid edge", good[-1] if good else None)
    return out


# ============================================================================= 7. cost, door to door, B/C
LV = {"low": 0, "base": 1, "high": 2}
COST = {  # $M per km of tube PAIR, installed, before soft costs and contingency (low, base, high) (A)
    "levitation / guidance track, SC-EDS ground coils, both bores": (4.0, 8.0, 16.0),
    "vacuum: pumps, sector gate valves every ~5 km, seals": (1.5, 3.0, 6.0),
    "power, moving-block signalling, comms": (2.0, 4.0, 8.0),
}
STEEL_USD_KG = (2.0, 3.0, 4.5)            # fabricated, coated, erected steel tube (pipeline / bridge steel) (A)
EGRESS = {"elevated": (1.5, 3.0, 6.0), "at-grade": (1.0, 2.0, 4.0), "buried": (0.0, 0.0, 0.0)}   # tunnels: in tunnel $
TRIM = {"elevated": (0.5, 1.5, 3.0), "at-grade": (0.5, 1.0, 2.0), "buried": (0.3, 0.8, 1.5)}
PYLONS = (8.0, 15.0, 30.0)                # 33 pylons/km + piles carrying both tubes (China HSR viaduct ~$9-11 M/km)
AT_GRADE = (6.0, 13.0, 26.0)              # slab / embankment 4-8-16 + road and wildlife crossings 2-5-10
TUNNEL5 = (50.0, 110.0, 250.0)            # twin TBM bores ~5 m ID, rural / mountain, incl. cross-passages
DEEP_X = 1.4                              # overburden > 800 m: heat, squeezing rock, access shafts
LAND = {"elevated": (0.5, 2.0, 6.0), "at-grade": (0.8, 3.0, 9.0), "buried": (0.2, 0.5, 1.0)}
LSM_C = sm.LSM                            # 8 / 15 / 30 on launch, brake and booster sections
BOOST = 0.05                              # share of cruise km with booster LSM (drag + EDS loss ~0.004 g)
POWER_B_GW = (0.3, 0.5, 1.0)              # $B per GW of station launch power above the 1 GW in sm.STATION
SUBURB_X = 2.0                            # express line through an intermediate metro (40 km) vs rural


def soft(i):
    return (1 + sm.SOFT[i]) * (1 + sm.CONTING[i])


def tube_steel(D, setting, i):
    """$M/km for the two steel tubes (regime ring-stiffened wall model at this bore); a tunnel needs only a liner."""
    from skyloop_physics import wall_monocoque
    m_mono, _ = wall_monocoque(100e3, D / 2, E=200e9, rho=7850.0)
    kg_m = 0.6 * m_mono * math.pi * D
    return 2 * kg_m * STEEL_USD_KG[i] / 1e3 * (0.3 if setting == "buried" else 1.0), kg_m


def class_cost(setting, D, lvl, trim=1):
    """$M/km (incl. soft + contingency) for each structure class of this setting."""
    i = LV[lvl]
    common = sum(v[i] for v in COST.values()) + (TRIM[setting][i] if trim else 0.0)
    steel, _ = tube_steel(D, setting, i)
    liner, _ = tube_steel(D, "buried", i)
    tun = TUNNEL5[i] * (D / 5.0) ** 1.5
    base_std = PYLONS[i] if setting == "elevated" else AT_GRADE[i]
    k = soft(i)
    c = dict(
        std=(common + steel + base_std + EGRESS[setting][i] + LAND[setting][i]) * k,
        cut=(common + steel + 1.5 * AT_GRADE[i] + EGRESS[setting][i] + LAND[setting][i]) * k,
        tall_unit=(PYLONS[i]) * k,                      # x (factor - 1) on top of std per km
        tunnel=(common + liner + tun + LAND["buried"][i]) * k,
        lsm=LSM_C[i] * k)
    c["deep_extra"] = (DEX := (DEEP_X - 1) * tun * k)
    return c


def leg_capex(rt, setting, D, lvl, trim=1):
    """$M for the rural part of one leg from its route classification (one tube pair, no LSM)."""
    c = class_cost(setting, D, lvl, trim)
    if setting == "buried":
        return rt["tunnel_km"] * c["tunnel"] + rt["deep_tunnel_km"] * c["deep_extra"]
    return (rt["std_km"] * c["std"] + rt["cut_km"] * c["cut"] + rt["tall_km"] * c["std"]
            + (rt["tall_factor_km"] - rt["tall_km"]) * c["tall_unit"] + rt["tunnel_km"] * c["tunnel"]
            + rt["deep_tunnel_km"] * c["deep_extra"])


def stage_cost():
    rows = []
    for st in SETTINGS:
        for D in (4.0, 5.0, 6.0, 8.0):
            for lvl in ("low", "base", "high"):
                i = LV[lvl]
                steel, kg = tube_steel(D, st, i)
                c = class_cost(st, D, lvl)
                r = dict(setting=st, bore_m=D, level=lvl, tube_steel_t_m=kg / 1e3, tube_or_liner=steel,
                         **{k: v[i] for k, v in COST.items()}, trim=TRIM[st][i], egress=EGRESS[st][i],
                         land=LAND[st][i], civil=(PYLONS[i] if st == "elevated" else AT_GRADE[i] if st == "at-grade"
                                                  else TUNNEL5[i] * (D / 5) ** 1.5),
                         soft_x_contingency=soft(i), total_std_usd_m_km=c["std"] if st != "buried" else c["tunnel"],
                         tunnel_usd_m_km=c["tunnel"], lsm_usd_m_km=c["lsm"])
                rows.append(r)
    write_csv("hl_cost.csv", rows)
    return rows


# ----------------------------------------------------------------------------- speed profiles and leg times
A_ACC = sm.P["a_acc"]
_RIDE0 = sm.ride_min           # original sky_model.ride_min; patched in-process only inside net_eval, then restored


def profile_time(L, vcap, V0=0.0, V1=0.0, a=A_ACC):
    """Minimum time over a leg of length L (m) sampled at DX with speed caps vcap; returns s, km of speed change."""
    n = len(vcap)
    v = vcap.copy()
    v[0], v[-1] = min(v[0], V0), min(v[-1], V1)
    for i in range(1, n):
        v[i] = min(v[i], math.sqrt(v[i - 1] ** 2 + 2 * a * DX))
    for i in range(n - 2, -1, -1):
        v[i] = min(v[i], math.sqrt(v[i + 1] ** 2 + 2 * a * DX))
    vm = np.maximum(0.5 * (v[1:] + v[:-1]), 1e-3)
    t = float(np.sum(DX / vm))
    dv = np.abs(np.diff(v)) > 1e-6
    dke = float(np.sum(np.maximum(v[:-1] ** 2 - v[1:] ** 2, 0)) / 2)     # J/kg braked in route (stations excluded)
    return t, float(dv.sum()) * DX / 1e3, v, dke


def follow_caps(H0, leg, setting, M, n_full):
    """Local Mach cap along the offset-0 alignment: the highest route Mach whose LP solution needs no tunnel there.
    Where every route Mach needs a tunnel, the tunnel is built and the cap is the cruise Mach.  The cap is a
    10 km running minimum so the pod does not hunt.  Urban ends are uncapped (costed as urban approaches)."""
    Ms = [m for m in M_ROUTE if m <= M + 1e-9]
    hs = [H0[f"{leg}|{setting}|{m:g}"] for m in Ms]
    n = len(hs[0])
    cap = np.full(n, M)
    which = np.full(n, len(Ms) - 1)
    any_ok = np.zeros(n, bool)
    for k, (m, h) in enumerate(zip(Ms, hs)):
        ok = h >= TUNNEL_BELOW[setting]
        first = ok & ~any_ok
        cap = np.where(first, m, np.where(ok, m, cap))
        which = np.where(ok, k, which)
        any_ok |= ok
    cap = np.where(any_ok, cap, M)
    w = int(5e3 / DX)
    capm = np.array([cap[max(0, i - w):i + w + 1].min() for i in range(n)])
    h_used = np.array([hs[which[i]][i] for i in range(n)])
    pad = (n_full - n) // 2
    full = np.concatenate([np.full(pad, M), capm, np.full(n_full - n - pad, M)])
    return full, h_used, capm


def classify_h(h, setting):
    """classify() for a stitched height array."""
    return classify(h, np.zeros_like(h), setting)


def leg_plan(summ, H0, terrain, setting, M, strategy):
    """Per leg: path length, ride time (s), LSM km, rural classification for the chosen strategy."""
    Mr = min(M_ROUTE, key=lambda m: abs(m - M)) if M not in M_ROUTE else M
    plans = {}
    V = M * C_REF
    for a, b in LEGS:
        leg = f"{a}-{b}"
        if strategy == "full speed" or setting == "buried":
            # interpolate classification between bracketing route Machs
            lo = max([m for m in M_ROUTE if m <= M] or [M_ROUTE[0]])
            hi = min([m for m in M_ROUTE if m >= M] or [M_ROUTE[-1]])
            rl, rh = summ[(setting, lo)]["legs"][leg], summ[(setting, hi)]["legs"][leg]
            f = 0.0 if hi == lo else (M - lo) / (hi - lo)
            rt = {k: (1 - f) * rl[k] + f * rh[k] for k in ("tunnel_km", "deep_tunnel_km", "cut_km", "std_km", "tall_km",
                                                         "tall_factor_km", "rural_km", "L_km")}
            L = rt["L_km"] * 1e3
            t = float(_RIDE0(L / 1e3, M)) * 60
            lsm_km = min(L, V * V / A_ACC) / 1e3 + BOOST * max(0.0, L - V * V / A_ACC) / 1e3
            vcap = None
        else:
            arr = terrain[(leg, 0.0)]
            L = arr[-1, 0] * 1e3
            n_full = len(arr)
            cap, h_used, _ = follow_caps(H0, leg, setting, Mr, n_full)
            vcap = np.minimum(cap, M) * C_REF
            t, km_dv, v, dke = profile_time(L, vcap)
            rt = classify_h(h_used, setting)
            rt["rural_km"] = len(h_used) * DX / 1e3
            rt["L_km"] = L / 1e3
            lsm_km = km_dv + BOOST * L / 1e3
            rt["dKE_J_kg"] = dke
        plans[leg] = dict(rt=rt, L=L, t=t, lsm_km=lsm_km, vcap=vcap)
    return plans


def express_time(plans, M, strategy):
    if strategy == "full speed":
        return float(_RIDE0(sum(p["L"] for p in plans.values()) / 1e3, M)) * 60
    vc = np.concatenate([p["vcap"] for p in plans.values()])
    t, _, _, _ = profile_time(len(vc) * DX, vc)
    return t




# ----------------------------------------------------------------------------- network evaluation
NETS = {"spine+LAX all-stop": [(sm.SPINE + ["LAX"], "allstop")], "NYC-LAX nonstop": [(["LGA", "LAX"], "express")]}
M_ECON = [1.0, 1.5, 2.0, 3.0, 4.0, 5.0, 6.0, 8.0, 10.0, 12.0]
STRATS = [("elevated", "full speed"), ("elevated", "terrain-following"), ("at-grade", "full speed"),
          ("at-grade", "terrain-following"), ("buried", "full speed")]
OPENINGS = (2045, 2050, 2055, 2060)
CMPS = {"sub": "subsonic", "bl": "Boomless", "m22": "BoomShield M2.2"}
POD_P = dict(seats=50, kg_per_seat=340.0, p_tube=10.0, cd_A=1.05, path=1.0)   # regime 17 t pod at 10 Pa (A)
_INP = None


def inp():
    global _INP
    if _INP is None:
        _INP = sm.load_inputs()
    return _INP


def _gc_km(a, b):
    return sm.gc(inp()["df"], a, b)


def net_eval(name, setting, strategy, M, summ, H0, terrain, D=5.0, F=None):
    """sky_model demand / benefit / energy for this tube network with hyperloop leg times; capex bottom-up."""
    F = F or {}
    plans = leg_plan(summ, H0, terrain, setting, M, strategy)
    table = {}
    for a, b in LEGS:
        table[round(_gc_km(a, b), 3)] = plans[f"{a}-{b}"]["t"] / 60
    table[round(_gc_km("LGA", "LAX"), 3)] = express_time(plans, M, strategy) / 60
    path_eff = sum(p["L"] for p in plans.values()) / 1e3 / sum(_gc_km(a, b) for a, b in LEGS)

    def ride(L_km, MM, a=None):
        key = round(float(L_km), 3)
        if key in table and abs(MM - M) < 1e-9:
            return table[key]
        return float(_RIDE0(float(L_km) * path_eff, MM, a))

    saved = {k: sm.P[k] for k in POD_P}
    sm.P.update(POD_P)
    sm.ride_min = ride
    try:
        net = sm.Net(name, NETS[name], M, inp()["df"])
        evs = {c: sm.evaluate(inp(), net, c) for c in CMPS}
        pods = {lvl: sm.pods_capex(evs["sub"], net, lvl) for lvl in LV}
        Rm = net.R.copy()
        k = net.k
    finally:
        sm.ride_min = _RIDE0
        sm.P.update(saved)
    # capacity: tube pairs on the busiest segment at a 10 s moving block (system convention, 1.3 peaking)
    cap_pair = sm.cap_pair_yr(M, F.get("rule", RULE), seats=POD_P["seats"])[0]
    pairs = int(min(6, max(1, math.ceil(evs["sub"]["peak"] * 1.3 / cap_pair))))
    # capex
    legs_used = LEGS if True else []
    n_st = len(net.st)
    ends = 2 * len(LEGS) if name.endswith("all-stop") else 2
    V = M * C_REF
    pods_h = PAX_H_DIR / (POD_P["seats"] * sm.P["lf"])
    P_dir = 0.5 * REF_POD * V * V * min(pods_h, 3600 / sm.headway_s(M, RULE)) / 3600 / ETA / 1e9
    capex = {}
    for lvl in LV:
        i = LV[lvl]
        c = class_cost(setting, D, lvl)
        rural = sum(leg_capex(plans[f"{a}-{b}"]["rt"], setting, D, lvl) for a, b in LEGS)
        lsm = sum(plans[f"{a}-{b}"]["lsm_km"] for a, b in LEGS) * c["lsm"]
        if name.endswith("all-stop"):
            suburb = 0.0
        else:   # nonstop line passes 5 intermediate metros: 40 km each at the suburban premium
            suburb = 5 * 2 * L_URB / 1e3 * (c["std"] * SUBURB_X if setting != "buried" else c["tunnel"])
            lsm = (V * V / A_ACC / 1e3 + BOOST * sum(p["L"] for p in plans.values()) / 1e3) * c["lsm"] \
                if strategy == "full speed" else lsm
        tube = (rural + lsm + suburb) * pairs / 1e3
        urban = ends * L_URB / 1e3 * sm.URBAN[i] * (1 if pairs == 1 else 0.6 + 0.4 * pairs) / 1e3
        n_launch = 2 * (n_st - 2) + 2 if name.endswith("all-stop") else 2
        stations = n_st * sm.STATION[i] + n_launch * max(0.0, P_dir - 0.5) * POWER_B_GW[i]
        capex[lvl] = dict(tube=tube, urban=urban, stations=stations, pods=pods[lvl],
                          total=tube + urban + stations + pods[lvl], rural=rural * pairs / 1e3)
    rural_km = sum(plans[f"{a}-{b}"]["rt"]["rural_km"] for a, b in LEGS)
    tun_km = sum(plans[f"{a}-{b}"]["rt"]["tunnel_km"] for a, b in LEGS)
    # energy correction: in-route braking on the terrain-following profile
    dke = sum(plans[f"{a}-{b}"]["rt"].get("dKE_J_kg", 0.0) for a, b in LEGS)
    e_scale = 1 + 0.7 * dke / (len(LEGS) * 0.5 * V * V)
    t_nyc_lax = Rm[k["LGA"], k["LAX"]]
    t_den_chi = Rm[k["DEN"], k["ORD"]] if "DEN" in k else float("nan")
    out = dict(network=name, setting=setting, strategy=strategy, Mach=M, bore_m=D, pairs=pairs,
               cap_pair_Mpax=cap_pair / 1e6, route_km=sum(p["L"] for p in plans.values()) / 1e3, rural_km=rural_km,
               tunnel_km=tun_km, tunnel_share=tun_km / rural_km,
               nyc_lax_d2d_h=(100 + t_nyc_lax) / 60, den_chi_d2d_h=(100 + t_den_chi) / 60,
               station_launch_GW_dir=P_dir, e_scale=e_scale)
    for lvl in LV:
        for kk, v in capex[lvl].items():
            out[f"capex_{kk}_{lvl}_busd"] = v
        out[f"avg_tube_usd_m_km_{lvl}"] = capex[lvl]["tube"] * 1e3 / (out["route_km"] * pairs)
    for cmp, ev in evs.items():
        t = ev["tot"]
        out.update({f"{cmp}_pax_M": t["captured_pax_M"] + t["net_induced_M"], f"{cmp}_benefit_busd": t["benefit_busd"],
                    f"{cmp}_fuel_Mt": t["fuel_saved_Mt"], f"{cmp}_kwh_TWh": t["kwh_TWh"] * e_scale,
                    f"{cmp}_wh_paxkm": t["wh_per_pax_km"] * e_scale, f"{cmp}_paxkm_B": t["pax_km_B"],
                    f"{cmp}_peak_M": ev["peak"] / 1e6, f"{cmp}_day_return_M": t["day_return_pax_M"]})
    return out, evs, capex


def bc_rows(o, capex, F=None):
    """Life-cycle B/C (sky_model.lifecycle): gross / net x capture x opening x comparator x cost level, + break-even."""
    F = F or {}
    rows, be = [], []
    cap_pax = o["pairs"] * o["cap_pair_Mpax"] * 1e6
    for cmp in CMPS:
        ann = dict(benefit_busd=o[f"{cmp}_benefit_busd"], fuel_saved_Mt=o[f"{cmp}_fuel_Mt"],
                   kwh_TWh=o[f"{cmp}_kwh_TWh"], pax_M=o[f"{cmp}_pax_M"], pax_km_B=o[f"{cmp}_paxkm_B"],
                   peak=o[f"{cmp}_peak_M"] * 1e6)
        build = 10 if o["route_km"] > 2000 else 7
        for oy in OPENINGS:
            for capture in (0.3, 1.0):
                r = dict(network=o["network"], setting=o["setting"], strategy=o["strategy"], Mach=o["Mach"],
                         vs=cmp, opening=oy, capture=capture)
                for lvl in LV:
                    lc = sm.lifecycle(ann, capex[lvl]["total"], sm.OM_FRAC[LV[lvl]], oy, capture, build=build,
                                      cap_pax=cap_pax)
                    r[f"bcg_{lvl}"], r[f"bcn_{lvl}"] = lc["bc_gross"], lc["bc_net"]
                rows.append(r)
                # break-even average tube cost ($M/km of pair), base level, other items fixed
                fixed = capex["base"]["total"] - capex["base"]["tube"]
                var = capex["base"]["tube"]
                lc1 = sm.lifecycle(ann, fixed + var, sm.OM_FRAC[1], oy, capture, build=build, cap_pax=cap_pax)
                lc0 = sm.lifecycle(ann, fixed, sm.OM_FRAC[1], oy, capture, build=build, cap_pax=cap_pax)
                for kind, num in (("gross", lc1["pv_ben"]), ("net", lc1["pv_ben"] + lc1["pv_energy"])):
                    kk = (num - lc0["pv_cost"]) / max(lc1["pv_cost"] - lc0["pv_cost"], 1e-9)
                    r[f"be_{kind}_usd_m_km"] = kk * o["avg_tube_usd_m_km_base"]
    return rows


def stage_econ(summ, H0, terrain, bore):
    pts, bcs = [], []
    for name in NETS:
        for st, strat in STRATS:
            for M in M_ECON:
                o, evs, cx = net_eval(name, st, strat, M, summ, H0, terrain, D=bore(st, M))
                pts.append(o)
                bcs += bc_rows(o, cx)
            print(f"  econ {name} {st} {strat} done", flush=True)
    write_csv("hl_points.csv", pts)
    write_csv("hl_econ.csv", bcs)
    return pts, bcs


# ============================================================================= levers
def _with(mod, **kw):
    old = {k: getattr(mod, k) for k in kw}
    for k, v in kw.items():
        setattr(mod, k, v)
    return old


def stage_levers(summ, summ_av, H0, terrain, bore):
    """Frontier (elevated, active 0.3 Hz gap + trim, 17 t pods) and spine B/C under single levers and all together."""
    out = []
    feas = [("base", {}, None),
            ("breach: fast valves + early braking (air met at 1/10 density)", {"breach_rho_div": 10.0}, None),
            ("stations 5 GW per direction", {"P_stn": 5e9}, None),
            ("regeneration 0.95 each way", {}, "regen"),
            ("vertical curves at 0.10 g (seated, belted)", {}, "av"),
            ("rolling pod 0.62 g in curves", {"a_lat": A_LAT["rolling pod (0.15 g felt + 25 deg)"]}, None),
            ("rail-like 0.24 g in curves", {"a_lat": A_LAT["rail-like (0.1 g felt + 8 deg cant)"]}, None),
            ("brick-wall headway (no moving-block certification)", {"rule": "brick-wall"}, None),
            ("all positive levers", {"breach_rho_div": 10.0, "P_stn": 5e9,
                                     "a_lat": A_LAT["rolling pod (0.15 g felt + 25 deg)"]}, "regen+av")]
    for name, F, special in feas:
        old = None
        sm_ = summ
        if special and "regen" in special:
            old = _with(R, ETA_D=0.95, ETA_R=0.95)
        if special and "av" in special:
            sm_ = summ_av
        try:
            rows = stage_grid(sm_, F, write=False, settings=["elevated"], pods=PODS, tag="lever ")
        finally:
            if old:
                _with(R, **old)
        fr = frontier([r for r in rows if r["trim"] == 1 and r["suspension"] == "active 0.3 Hz"])
        fm = frontier([r for r in rows if r["trim"] == 1 and r["suspension"] == "active 0.3 Hz"],
                      ok=("pass", "marginal"))
        k = ("elevated", 1, "active 0.3 Hz")
        out.append(dict(kind="feasibility", lever=name, Mach_pass=fr[k][0], binding_above=fr[k][1],
                        Mach_marginal=fm[k][0], binding_above_marginal=fm[k][1]))
        print(f"  lever {name}: pass M{fr[k][0]} ({fr[k][1]}), marginal M{fm[k][0]} ({fm[k][1]})", flush=True)
    # economics: spine+LAX all-stop, elevated full speed, net vs subsonic, 2050, 100 %, at M2 / M4 / M6
    econ = [("base", {}, None), ("tunnels at $50 M/km base (TBM cost-down)", {}, "tun"),
            ("elevated civil + tube at the low-cost level", {}, "low"),
            ("vertical curves at 0.10 g", {}, "av"), ("70 min walk-on stations", {}, "walk"),
            ("brick-wall headway", {"rule": "brick-wall"}, None)]
    global TUNNEL5
    for name, F, special in econ:
        for M in (2.0, 4.0, 6.0):
            t5 = TUNNEL5
            o_sky = sm.P["O_sky"]
            try:
                if special == "tun":
                    TUNNEL5 = (25.0, 50.0 / soft(1) * 1.0 * 1.0 if False else 50.0 / 1.625 * 1.0, 100.0)
                if special == "walk":
                    sm.P["O_sky"] = (20, 30, 20)
                o, evs, cx = net_eval("spine+LAX all-stop", "elevated", "full speed", M,
                                      summ_av if special == "av" else summ, H0, terrain, D=bore("elevated", M), F=F)
                if special == "low":
                    cx = {k: dict(v) for k, v in cx.items()}
                    low_t = cx["low"]["tube"]
                    cx["base"]["total"] += low_t - cx["base"]["tube"]
                    cx["base"]["tube"] = low_t
                b = [r for r in bc_rows(o, cx) if r["vs"] == "sub" and r["opening"] == 2050 and r["capture"] == 1.0][0]
                bl = [r for r in bc_rows(o, cx) if r["vs"] == "bl" and r["opening"] == 2050 and r["capture"] == 1.0][0]
            finally:
                TUNNEL5 = t5
                sm.P["O_sky"] = o_sky
            out.append(dict(kind="B/C spine net vs subsonic 2050 100 %", lever=name, Mach=M, pairs=o["pairs"],
                            avg_tube_usd_m_km=o["avg_tube_usd_m_km_base"], nyc_lax_d2d_h=o["nyc_lax_d2d_h"],
                            bcn_base=b["bcn_base"], bcn_low=b["bcn_low"], bcg_base=b["bcg_base"],
                            bcn_base_vs_boomless=bl["bcn_base"]))
            print(f"  lever {name} M{M}: net B/C base {b['bcn_base']:.2f} low {b['bcn_low']:.2f}", flush=True)
    write_csv("hl_levers.csv", out)
    return out
