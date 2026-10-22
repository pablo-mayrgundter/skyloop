#!/usr/bin/env python3
"""
BoomShield traffic coverage: how much US long-haul O&D does the trunk address?

Stages (each cached in .cache/ so a restarted run resumes):
  1. Demand.  Real O&D from DOT Consumer Airfare Report Table 1 (top-1,000 contiguous-US
     city-pair markets, DB1B 10 % ticket sample scaled to 100 %; 2025 Q1-Q4 average,
     one-way trips/day, both directions combined), mirrored in a public GitHub repo
     (Kaichung862/mist6380_project, Consumer_Airfare.csv; extract saved as
     dot_table1_2024_2025.csv).  Pairs below the Table 1 cutoff (~200-240 pax/day) among
     the same 129 city markets are filled with a log-linear gravity model (city fixed
     effects + distance-band effects) fitted to Table 1 and clipped at the cutoff.
     National scale from BTS T-100 Domestic Segment totals (845.2 M enplanements 2025,
     856.4 M 2024; vaibhavsoni85-spec/US-Domestic-Passenger-Traffic-Analysis-SQL report).
  2. Baseline times.  routes.md / build_trunk.py time model (imported).  Boomless
     baseline = min(Boomless great circle, Boomless feeders + M1.7 over open-water
     tracks >= 40 km offshore: Atlantic, Florida Straits + Gulf, Pacific).
  3. Addressable traffic: pairs whose own great-circle shield (M1.7 / M2.0, same ramps,
     5 min merge) saves >= 10 min vs that baseline.  Upper-bound pax-h saved.
  4. Spine coverage (LAS-DEN-OMA-CHI-CLE-NYC) and the biggest uncaptured flows with a
     reason (ramp overhead / ocean baseline / geometry-detour / partial).
  5. Network expansions: named scenarios + greedy growth from the spine on a
     metro + Steiner-grid candidate graph; coverage curve and B/C vs trunk km with
     ../../economics/cost_model.py (base cost, with and without the lid fatigue cap).
  6. Figures.

Usage: python3 build_coverage.py [--stage N]     (one process; ~20-40 min cold)
"""
from __future__ import annotations

import argparse
import json
import math
import pickle
import sys
from pathlib import Path

import networkx as nx
import numpy as np
import pandas as pd
from scipy.spatial import Delaunay

HERE = Path(__file__).resolve().parent
NM = HERE.parent
CACHE = HERE / ".cache"
FIG = HERE / "figures"
sys.path.insert(0, str(NM / "trunk"))
import build_trunk as bt  # noqa: E402  (also imports ../../economics/cost_model.py as bt.cm)

cm = bt.cm
KM_PER_MI = 1.609344

# ------------------------------------------------------------------------------------
# National totals (verified where marked)
# ------------------------------------------------------------------------------------
T100_ENPL = {2024: 856.42e6, 2025: 845.17e6}   # BTS T-100 Domestic Segment, US carriers, class F
ENPL_PER_TRIP = (1.20, 1.30, 1.40)              # low / central / high; [U] DB1B coupons per directional trip
DOT_TOP1000_SHARE = 0.47                        # DOT report text (2017-2024 editions): top 1,000 = 47 % of all domestic pax

# City market -> airport whose coordinates represent it (main airport of the DOT city market)
CITY_AP = {
    "Albany, NY": "ALB", "Albuquerque, NM": "ABQ", "Allentown/Bethlehem/Easton, PA": "ABE", "Amarillo, TX": "AMA",
    "Appleton, WI": "ATW", "Asheville, NC": "AVL", "Aspen, CO": "ASE", "Atlanta, GA (Metropolitan Area)": "ATL",
    "Atlantic City, NJ": "ACY", "Austin, TX": "AUS", "Bangor, ME": "BGR", "Bend/Redmond, OR": "RDM",
    "Birmingham, AL": "BHM", "Bismarck/Mandan, ND": "BIS", "Boise, ID": "BOI", "Boston, MA (Metropolitan Area)": "BOS",
    "Bozeman, MT": "BZN", "Buffalo, NY": "BUF", "Burlington, VT": "BTV", "Cedar Rapids/Iowa City, IA": "CID",
    "Charleston, SC": "CHS", "Charlotte, NC": "CLT", "Charlottesville, VA": "CHO", "Chicago, IL": "ORD",
    "Cincinnati, OH": "CVG", "Cleveland, OH (Metropolitan Area)": "CLE", "Colorado Springs, CO": "COS",
    "Columbia, SC": "CAE", "Columbus, OH": "CMH", "Dallas/Fort Worth, TX": "DFW", "Dayton, OH": "DAY",
    "Denver, CO": "DEN", "Des Moines, IA": "DSM", "Detroit, MI": "DTW", "Eagle, CO": "EGE", "El Paso, TX": "ELP",
    "Eugene, OR": "EUG", "Everett, WA": "PAE", "Fargo, ND": "FAR", "Fayetteville, AR": "XNA", "Flint, MI": "FNT",
    "Fort Myers, FL": "RSW", "Fort Wayne, IN": "FWA", "Fresno, CA": "FAT", "Grand Rapids, MI": "GRR",
    "Greensboro/High Point, NC": "GSO", "Greenville/Spartanburg, SC": "GSP", "Harlingen/San Benito, TX": "HRL",
    "Harrisburg, PA": "MDT", "Hartford, CT": "BDL", "Houston, TX": "IAH", "Huntsville, AL": "HSV",
    "Indianapolis, IN": "IND", "Jackson, WY": "JAC", "Jackson/Vicksburg, MS": "JAN", "Jacksonville, FL": "JAX",
    "Kalispell, MT": "FCA", "Kansas City, MO": "MCI", "Key West, FL": "EYW", "Knoxville, TN": "TYS",
    "Las Vegas, NV": "LAS", "Little Rock, AR": "LIT", "Los Angeles, CA (Metropolitan Area)": "LAX",
    "Louisville, KY": "SDF", "Lubbock, TX": "LBB", "Madison, WI": "MSN", "Martha's Vineyard, MA": "MVY",
    "Medford, OR": "MFR", "Memphis, TN": "MEM", "Miami, FL (Metropolitan Area)": "MIA", "Midland/Odessa, TX": "MAF",
    "Milwaukee, WI": "MKE", "Minneapolis/St. Paul, MN": "MSP", "Mission/McAllen/Edinburg, TX": "MFE",
    "Myrtle Beach, SC": "MYR", "Nantucket, MA": "ACK", "Nashville, TN": "BNA", "New Haven, CT": "HVN",
    "New Orleans, LA": "MSY", "New York City, NY (Metropolitan Area)": "LGA", "Norfolk, VA (Metropolitan Area)": "ORF",
    "Oklahoma City, OK": "OKC", "Omaha, NE": "OMA", "Orlando, FL": "MCO", "Palm Springs, CA": "PSP",
    "Panama City, FL": "ECP", "Pasco/Kennewick/Richland, WA": "PSC", "Pensacola, FL": "PNS", "Peoria, IL": "PIA",
    "Philadelphia, PA": "PHL", "Phoenix, AZ": "PHX", "Pittsburgh, PA": "PIT", "Portland, ME": "PWM",
    "Portland, OR": "PDX", "Provo, UT": "PVU", "Punta Gorda, FL": "PGD", "Raleigh/Durham, NC": "RDU",
    "Reno, NV": "RNO", "Richmond, VA": "RIC", "Rochester, NY": "ROC", "Sacramento, CA": "SMF",
    "Salt Lake City, UT": "SLC", "San Antonio, TX": "SAT", "San Diego, CA": "SAN",
    "San Francisco, CA (Metropolitan Area)": "SFO", "Sanford, FL": "SFB", "Santa Barbara, CA": "SBA",
    "Santa Rosa, CA": "STS", "Sarasota/Bradenton, FL": "SRQ", "Savannah, GA": "SAV", "Seattle, WA": "SEA",
    "Sioux Falls, SD": "FSD", "Spokane, WA": "GEG", "St. Louis, MO": "STL", "Steamboat Springs, CO": "HDN",
    "Syracuse, NY": "SYR", "Tallahassee, FL": "TLH", "Tampa, FL (Metropolitan Area)": "TPA", "Trenton, NJ": "TTN",
    "Tucson, AZ": "TUS", "Tulsa, OK": "TUL", "Valparaiso, FL": "VPS", "Vero Beach, FL": "VRB",
    "Washington, DC (Metropolitan Area)": "DCA", "West Palm Beach/Palm Beach, FL": "PBI", "Wichita, KS": "ICT",
    "Wilmington, NC": "ILM", "Worcester, MA": "ORH",
}
# trunk.md metro code -> code used here (same coordinates as build_trunk.METROS for these)
TRUNK2HERE = {"NYC": "LGA", "CHI": "ORD", "WAS": "DCA", "SFO": "SFO", "FLL": "MIA", "HOU": "IAH"}

SPINE = [("LGA", "CLE"), ("CLE", "ORD"), ("ORD", "OMA"), ("OMA", "DEN"), ("DEN", "LAS")]
CUT_MIN = 10.0
MACHS = (1.7, 2.0)
CFG = dict(bt.BASE_CFG, R=40)    # same aircraft, best feeder, 5 min merge, cut 10 min; 40 nearest ramp points
R_GREEDY = 10                    # ramp points per metro inside the greedy search (speed); final evals use CFG


def cached(name, fn):
    CACHE.mkdir(exist_ok=True)
    f = CACHE / f"{name}.pkl"
    if f.exists():
        return pickle.loads(f.read_bytes())
    out = fn()
    f.write_bytes(pickle.dumps(out))
    return out


# =====================================================================================
# 1. Demand
# =====================================================================================
def load_table1(year=2025):
    d = pd.read_csv(HERE / "dot_table1_2024_2025.csv")
    d = d[d.Year == year]
    nq = d.quarter.nunique()
    g = d.groupby(["city1", "city2"]).agg(nsmiles=("nsmiles", "first"), pax_day=("passengers", "sum"),
                                         fare=("fare", "mean")).reset_index()
    g["pax_day"] /= nq                  # markets absent in a quarter count as below cutoff there (~0)
    g["pax_yr"] = g.pax_day * 365
    cut_day = float(d.groupby("quarter").passengers.min().mean())
    return g, cut_day, nq


def city_table(g):
    import airportsdata
    ap = airportsdata.load("IATA")
    cities = sorted(set(g.city1) | set(g.city2))
    miss = [c for c in cities if c not in CITY_AP]
    assert not miss, miss
    rows = []
    for c in cities:
        a = ap[CITY_AP[c]]
        rows.append(dict(code=CITY_AP[c], city=c, lat=a["lat"], lon=a["lon"]))
    df = pd.DataFrame(rows)
    # keep trunk.md coordinates for the metros the trunk study used (identical geometry)
    tm = bt.metros_df().set_index("code")
    for tcode, here in list(TRUNK2HERE.items()) + [(c, c) for c in ("CLE", "OMA", "DEN", "LAS", "LAX", "SEA", "ATL",
                                                                     "BOS", "DFW", "PHX", "MSP", "MCO", "TPA")]:
        if tcode in tm.index and here in set(df.code):
            df.loc[df.code == here, ["lat", "lon"]] = tm.loc[tcode, ["lat", "lon"]].values
    return df.reset_index(drop=True)


DBANDS = [0, 500, 1000, 1500, 2000, 2500, 3000, 3500, 5000]


def gravity_fill(df, D, obs, cut_yr):
    """log T_ij = a_i + a_j + f(distance band) fitted on observed Table 1 pairs (OLS in logs).
    Predict unobserved pairs, clip at the Table 1 cutoff.  Returns symmetric matrix of fills."""
    n = len(df)
    bands = np.array([300, 500, 700, 1000, 1300, 1600, 2000, 2500, 3000, 3600, 6000])
    ii, jj, y = obs
    def X_of(i, j):
        d = D[i, j]
        b = np.clip(np.searchsorted(bands, d), 0, len(bands) - 1)
        X = np.zeros((len(i), n + len(bands)))
        X[np.arange(len(i)), i] += 1
        X[np.arange(len(i)), j] += 1
        X[np.arange(len(i)), n + b] = 1
        return X
    X = X_of(ii, jj)
    coef, *_ = np.linalg.lstsq(X, np.log(y), rcond=None)
    res = np.log(y) - X @ coef
    s2 = float(res.var())
    iu, ju = np.triu_indices(n, 1)
    pred = np.exp(X_of(iu, ju) @ coef + s2 / 2)
    F = np.zeros((n, n))
    F[iu, ju] = np.minimum(pred, cut_yr)
    F = F + F.T
    F[ii, jj] = 0
    F[jj, ii] = 0
    return F, float(np.sqrt(s2))


def stage_demand():
    g, cut_day, nq = load_table1(2025)
    df = city_table(g)
    idx = {c: k for k, c in zip(df.code, df.city)}
    pos = {c: i for i, c in enumerate(df.code)}
    D = bt.dist_matrix(df)
    n = len(df)
    T = np.zeros((n, n))
    ii = np.array([pos[idx[c]] for c in g.city1])
    jj = np.array([pos[idx[c]] for c in g.city2])
    T[ii, jj] = g.pax_yr.values
    T[jj, ii] = g.pax_yr.values
    cut_yr = cut_day * 365
    F, sig = gravity_fill(df, D, (ii, jj, g.pax_yr.values), cut_yr)
    print(f"Table 1 2025: {len(g)} markets, {g.pax_yr.sum() / 1e6:.1f} M O&D trips/yr; cutoff {cut_day:.0f}/day; "
          f"fill among {n} cities: {F[np.triu_indices(n, 1)].sum() / 1e6:.1f} M (rms log {sig:.2f})")
    # distance check: DOT nonstop miles vs our great circle
    dd = g.nsmiles.values * KM_PER_MI / D[ii, jj]
    print(f"DOT nsmiles / GC ratio median {np.median(dd):.3f}, 5-95 % {np.percentile(dd, 5):.2f}-{np.percentile(dd, 95):.2f}")
    return dict(df=df, D=D, T_obs=T, T_fill=F, cut_yr=cut_yr, sig=sig, n_markets=len(g))


def band_table(D, T_obs, T_fill, nat_trips):
    iu, ju = np.triu_indices(len(D), 1)
    d = D[iu, ju]
    o, f = T_obs[iu, ju], T_fill[iu, ju]
    rows = []
    loc_tot = (o + f).sum()
    for lo in (0, 1500, 2500, 3500):
        m = d > lo
        rows.append(dict(band=f">{lo} km" if lo else "all", pax_table1_M=o[m].sum() / 1e6,
                         pax_fill_M=f[m].sum() / 1e6, pax_located_M=(o + f)[m].sum() / 1e6,
                         pax_km_located_bn=((o + f) * d)[m].sum() / 1e9,
                         share_of_located_pax=(o + f)[m].sum() / loc_tot,
                         share_of_located_paxkm=((o + f) * d)[m].sum() / ((o + f) * d).sum(),
                         national_scaled_M=(o + f)[m].sum() / loc_tot * nat_trips / 1e6))
    return pd.DataFrame(rows)


def gravity_vs_dot(dem):
    """Compare trunk.md's gravity model (od_matrix.csv) with Table 1 for shared pairs."""
    od = pd.read_csv(NM / "trunk" / "od_matrix.csv", usecols=["a", "b", "distance_km", "pax_yr"])
    df = dem["df"]
    pos = {c: i for i, c in enumerate(df.code)}
    rows = []
    for a, b, d, p in od.itertuples(index=False):
        ca, cb = TRUNK2HERE.get(a, a), TRUNK2HERE.get(b, b)
        if ca in pos and cb in pos:
            i, j = pos[ca], pos[cb]
            rows.append(dict(a=a, b=b, distance_km=d, gravity_pax=p, dot_pax=dem["T_obs"][i, j],
                             dot_or_fill_pax=dem["T_obs"][i, j] + dem["T_fill"][i, j]))
    r = pd.DataFrame(rows)
    r["ratio_gravity_to_dot"] = r.gravity_pax / r.dot_or_fill_pax.clip(lower=1)
    return r


# =====================================================================================
# 2. Baseline (Boomless + open-water M1.7) and dedicated-shield upper bound
# =====================================================================================
# offshore M1.7 tracks (>= ~45 km from the coast; author's waypoints)
OCEAN = {
    "atl_gulf": [("A0", 42.0, -69.3), ("A1", 40.1, -72.0), ("A2", 38.6, -73.9), ("A3", 36.5, -74.8),
                 ("A4", 35.0, -74.9), ("A5", 33.3, -77.3), ("A6", 31.8, -79.6), ("A7", 30.0, -80.7),
                 ("A8", 28.3, -79.9), ("A9", 26.6, -79.5), ("A10", 25.2, -79.6), ("F1", 24.2, -81.0),
                 ("F2", 24.4, -82.8), ("G1", 26.3, -83.4), ("G2", 27.8, -83.6), ("G3", 28.7, -86.5),
                 ("G4", 28.2, -89.6), ("G5", 28.3, -92.5), ("G6", 28.3, -94.6), ("G7", 27.0, -96.7)],
    "pacific": [("P1", 48.0, -125.5), ("P2", 46.0, -124.8), ("P3", 43.5, -125.2), ("P4", 41.0, -125.0),
                ("P5", 38.4, -124.2), ("P6", 37.3, -123.2), ("P7", 35.5, -122.0), ("P8", 34.0, -121.0),
                ("P9", 32.6, -118.8), ("P10", 32.0, -117.8)],
}


def ext_df(df, extra):
    """Cities followed by extra (code, lat, lon) points (zero demand)."""
    e = pd.DataFrame([dict(code=c, city=c, lat=la, lon=lo) for c, la, lo in extra])
    return pd.concat([df, e], ignore_index=True)


def eval_on(trunk, D, mach, cfg, pairs):
    """Trunk routing: returns t_tr, on, off, and geometry (feeder + trunk km)."""
    iu, ju = pairs
    idx, dist = trunk.nearest(cfg["R"])
    a_idx, a_d = idx[iu][:, :, None], dist[iu][:, :, None]
    b_idx, b_d = idx[ju][:, None, :], dist[ju][:, None, :]
    L = trunk.Dt[a_idx, b_idx]
    best = np.full(L.shape, np.inf)
    for f in (("sub",) if cfg["feeder"] == "sub" else ("sub", "bl")):
        t, _ = bt.t_trunk(a_d, L, b_d, mach, f, cfg["merge_min"], cfg["a_acc"], cfg["a_dec"])
        best = np.minimum(best, t)
    flat = best.reshape(len(iu), -1)
    k = np.argmin(flat, 1)
    R = idx.shape[1]
    r = np.arange(len(iu))
    on, off = idx[iu, k // R], idx[ju, k % R]
    return dict(t=flat[r, k], on=on, off=off, d_on=dist[iu, k // R], d_off=dist[ju, k % R],
                L=trunk.Dt[on, off])


def stage_baseline(dem):
    df, D = dem["df"], dem["D"]
    n = len(df)
    iu, ju = np.triu_indices(n, 1)
    d = D[iu, ju]
    t_sub = bt.t_direct(d, "sub")
    t_bl_gc = bt.t_direct(d, "bl")
    t_ocean = np.full(len(iu), np.inf)
    ocean_name = np.array([""] * len(iu), dtype=object)
    for name, pts in OCEAN.items():
        dfx = ext_df(df, pts)
        edges = [(pts[k][0], pts[k + 1][0]) for k in range(len(pts) - 1)]
        tr = bt.Trunk(dfx, edges, step=100.0)
        r = eval_on(tr, None, 1.7, dict(CFG, merge_min=0.0), (iu, ju))
        upd = r["t"] < t_ocean
        t_ocean = np.where(upd, r["t"], t_ocean)
        ocean_name = np.where(upd & np.isfinite(r["t"]), name, ocean_name)
    t_bl = np.minimum(t_bl_gc, t_ocean)
    out = dict(iu=iu, ju=ju, d=d, t_sub=t_sub, t_bl_gc=t_bl_gc, t_ocean=t_ocean, ocean=ocean_name, t_bl=t_bl)
    for M in MACHS:
        ub = np.minimum(*(bt.t_trunk(0.0, d, 0.0, M, f, CFG["merge_min"])[0] for f in ("sub", "bl")))
        out[f"t_ub_M{M}"] = ub
    return out


# =====================================================================================
# 3-5. Network evaluation, economics, greedy growth
# =====================================================================================
BIG = 1e7


def make_trunk(dfx, edges):
    tr = bt.Trunk(dfx, edges)
    comp = np.full(len(tr.ll), -1)
    for c, nodes in enumerate(nx.connected_components(tr.G)):
        comp[list(nodes)] = c
    tr.Dt = np.where(comp[:, None] == comp[None, :], tr.Dt, BIG)
    return tr


def net_eval(dfx, base, T, edges, mach, sel):
    """Evaluate a trunk network for pairs `sel` (indices into base arrays)."""
    iu, ju = base["iu"][sel], base["ju"][sel]
    tr = make_trunk(dfx, edges)
    r = eval_on(tr, None, mach, CFG, (iu, ju))
    t = r["t"]
    pax = T[iu, ju]
    s_bl = base["t_bl"][sel] - t
    s_bl = np.where(s_bl >= CUT_MIN, s_bl, 0.0)
    s_sub = base["t_sub"][sel] - t
    s_sub = np.where(s_sub >= CUT_MIN, s_sub, 0.0)
    return dict(tr=tr, sel=sel, t=t, pax=pax, s_bl=s_bl, s_sub=s_sub, on=r["on"], off=r["off"],
                d_on=r["d_on"], d_off=r["d_off"], L=r["L"], B_bl=float((pax * s_bl).sum() / 60),
                B_sub=float((pax * s_sub).sum() / 60), km=tr.km)


def econ(ne, mach):
    """B/C at 100 % / 10 % capture, base cost; lid-capped (100 flights/day per km) and uncapped."""
    ann = bt.annual_cost_per_km(mach, "base")[0]
    vt = cm.VTTS["base"] * cm.v("availability", "base")
    km = ne["km"]
    ev = dict(pax=ne["pax"], sav=ne["s_bl"], on=ne["on"], off=ne["off"], use=ne["s_bl"] > 0)
    lb = bt.segment_loads(ne["tr"], ev, "sav", "use")
    seg_pax = max([v[0] for v in lb.values()] or [0])
    out = dict(cost_musd_km_yr=ann / 1e6, peak_seg_pax_M=seg_pax / 1e6,
               peak_flights_day=seg_pax / (bt.SEATS * bt.LOAD_FACTOR * 365))
    per_flight = bt.SEATS * bt.LOAD_FACTOR * 365
    for c in (1.0, 0.1):
        ben = ne["B_bl"] * c * vt
        out[f"bcr_bl_cap{int(c * 100)}_nofatigue"] = ben / (ann * km)
        # lid fatigue caps: economics.md 100 flights/day; ../../fatigue revised low 1,000 / base 10,000 per day
        for cap, tag in ((100, "lidcapped"), (1000, "lid1000"), (10000, "lid10000")):
            ben_cap = sum(v[1] * c * min(1.0, cap * per_flight / max(v[0] * c, 1)) for v in lb.values()) * vt
            out[f"bcr_bl_cap{int(c * 100)}_{tag}"] = ben_cap / (ann * km)
        out[f"bcr_sub_cap{int(c * 100)}_nofatigue"] = ne["B_sub"] * c * vt / (ann * km)
    return out, lb


def steiner_grid(dlat=3.0, dlon=3.5):
    from matplotlib.path import Path as MPath
    gj = json.loads((NM / "routes" / ".cache" / "ne_50m_land.geojson").read_text())
    polys = []
    for ft in gj["features"]:
        g = ft["geometry"]
        for poly in (g["coordinates"] if g["type"] == "MultiPolygon" else [g["coordinates"]]):
            ring = np.array(poly[0])
            if ring[:, 0].min() < -60 and ring[:, 0].max() > -130 and ring[:, 1].max() > 25 and len(ring) > 500:
                polys.append(MPath(ring))
    pts = []
    k = 0
    for la in np.arange(30.5, 48.6, dlat):
        for lo in np.arange(-121.0, -69.0, dlon):
            if any(P.contains_point((lo, la)) for P in polys):
                # keep >= ~80 km inland: all 4 neighbours at 0.75 deg also on land
                if all(any(P.contains_point((lo + dx, la + dy)) for P in polys)
                       for dx, dy in ((0.9, 0), (-0.9, 0), (0, 0.75), (0, -0.75))):
                    pts.append((f"S{k:02d}", float(la), float(lo)))
                    k += 1
    return pts


def candidate_graph(dfx, n_city, node_ok, max_edge=950.0, direct_km=1150.0):
    lat0 = math.radians(dfx.lat.mean())
    keep = np.where(node_ok)[0]
    P = np.c_[dfx.lon.values[keep] * math.cos(lat0), dfx.lat.values[keep]]
    Dx = bt.dist_matrix(dfx)
    G = nx.Graph()
    for s in Delaunay(P).simplices:
        for a, b in ((s[0], s[1]), (s[1], s[2]), (s[0], s[2])):
            a, b = int(keep[a]), int(keep[b])
            if Dx[a, b] <= max_edge:
                G.add_edge(a, b, length=float(Dx[a, b]))
    for a in keep:
        for b in keep:
            if a < b and a < n_city and b < n_city and Dx[a, b] <= direct_km:
                G.add_edge(int(a), int(b), length=float(Dx[a, b]))
    return G


def greedy(dfx, base, T, mach, sel, start_edges, cand_nodes, G, max_steps=30, log=print, tag="",
           resume=(), km_max=14000.0, stop_frac=0.10):
    """Grow by the attachment path (tree -> candidate city, shortest in G) with the highest marginal
    pax-h/km.  `resume`: node-sequence strings of steps already taken (from a log)."""
    R_final = CFG["R"]
    CFG["R"] = R_GREEDY
    try:
        return _greedy(dfx, base, T, mach, sel, start_edges, cand_nodes, G, max_steps, log, tag, resume, km_max,
                       stop_frac)
    finally:
        CFG["R"] = R_final


def _greedy(dfx, base, T, mach, sel, start_edges, cand_nodes, G, max_steps, log, tag, resume, km_max, stop_frac):
    code = list(dfx.code)
    edges = list(start_edges)
    ne = net_eval(dfx, base, T, edges, mach, sel)
    rows = [dict(step=0, added="start:" + " ".join(f"{a}-{b}" for a, b in edges), km=ne["km"], B_bl=ne["B_bl"],
                 B_sub=ne["B_sub"], marginal_per_km=ne["B_bl"] / ne["km"], edges=list(edges))]
    for k, path in enumerate(resume, 1):
        nodes = path.split("-")
        ce = list(zip(nodes[:-1], nodes[1:]))
        nn = net_eval(dfx, base, T, edges + ce, mach, sel)
        r = (nn["B_bl"] - ne["B_bl"]) / max(nn["km"] - ne["km"], 1)
        edges, ne = edges + ce, nn
        rows.append(dict(step=k, added=path, km=ne["km"], B_bl=ne["B_bl"], B_sub=ne["B_sub"], marginal_per_km=r,
                         edges=list(edges)))
    if resume and (rows[-1]["marginal_per_km"] < stop_frac * ne["B_bl"] / ne["km"] or ne["km"] > km_max):
        return rows
    for step in range(len(resume) + 1, max_steps + 1):
        tree_nodes = {code.index(x) for e in edges for x in e}
        dist, paths = nx.multi_source_dijkstra(G, tree_nodes, weight="length")
        res = []
        for v in cand_nodes:
            if v in tree_nodes or v not in paths:
                continue
            p = paths[v]
            ce = [(code[a], code[b]) for a, b in zip(p[:-1], p[1:])]
            nn = net_eval(dfx, base, T, edges + ce, mach, sel)
            res.append(((nn["B_bl"] - ne["B_bl"]) / max(nn["km"] - ne["km"], 1), nn, ce))
        if not res:
            break
        res.sort(key=lambda r: -r[0])
        r, nn, ce = res[0]
        if r <= 0:
            break
        edges = edges + ce
        ne = nn
        rows.append(dict(step=step, added="-".join([ce[0][0]] + [b for _, b in ce]), km=ne["km"], B_bl=ne["B_bl"],
                         B_sub=ne["B_sub"], marginal_per_km=r, edges=list(edges)))
        log(f"  {tag} M{mach} step {step}: +{rows[-1]['added']} km={ne['km']:.0f} B={ne['B_bl'] / 1e6:.2f}M "
            f"marg={r:.0f}/km avg={ne['B_bl'] / ne['km']:.0f}/km")
        if r < stop_frac * ne["B_bl"] / ne["km"] or ne["km"] > km_max:
            break
    return rows


SCEN = {
    "spine": SPINE,
    "spine+LAX": SPINE + [("LAS", "LAX")],
    "spine+north SEA-MSP-CHI": SPINE + [("SEA", "MSP"), ("MSP", "ORD")],
    "spine+Florida CHI-ATL-MCO-MIA": SPINE + [("ORD", "ATL"), ("ATL", "MCO"), ("MCO", "MIA")],
    "spine+east coast BOS-NYC-WAS-CLT-ATL-MCO-MIA": SPINE + [("BOS", "LGA"), ("LGA", "DCA"), ("DCA", "CLT"),
                                                            ("CLT", "ATL"), ("ATL", "MCO"), ("MCO", "MIA")],
    "south only LAX-PHX-DFW-ATL-WAS": [("LAX", "PHX"), ("PHX", "DFW"), ("DFW", "ATL"), ("ATL", "DCA")],
    "two-spine (spine+LAX+south LAX-PHX-DFW-ATL-WAS)": SPINE + [("LAS", "LAX"), ("LAX", "PHX"), ("PHX", "DFW"),
                                                                 ("DFW", "ATL"), ("ATL", "DCA")],
    "two-spine + ATL-MCO-MIA": SPINE + [("LAS", "LAX"), ("LAX", "PHX"), ("PHX", "DFW"), ("DFW", "ATL"),
                                        ("ATL", "DCA"), ("ATL", "MCO"), ("MCO", "MIA")],
    "X at Denver (spine+LAX + SEA-SLC-DEN-DFW-IAH)": SPINE + [("LAS", "LAX"), ("SEA", "SLC"), ("SLC", "DEN"),
                                                              ("DEN", "DFW"), ("DFW", "IAH")],
    "X + north + south (spine+LAX+SEA-MSP-CHI+south+ATL-MIA)": SPINE + [
        ("LAS", "LAX"), ("SEA", "MSP"), ("MSP", "ORD"), ("LAX", "PHX"), ("PHX", "DFW"), ("DFW", "ATL"),
        ("ATL", "DCA"), ("ATL", "MCO"), ("MCO", "MIA")],
}


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--stage", type=int, default=99)
    args = ap.parse_args()
    dem = cached("demand", stage_demand)
    df, D = dem["df"], dem["D"]
    T = dem["T_obs"] + dem["T_fill"]
    nat = {k: T100_ENPL[2025] / k for k in ENPL_PER_TRIP}
    bands = band_table(D, dem["T_obs"], dem["T_fill"], nat[1.30])
    print(bands.round(3).to_string())
    gv = gravity_vs_dot(dem)
    gv.round(3).to_csv(HERE / "gravity_vs_dot.csv", index=False)
    lh = gv[gv.distance_km > 2500]
    print(f"gravity/DOT: all pairs sum {gv.gravity_pax.sum() / gv.dot_or_fill_pax.sum():.2f}; "
          f">2500 km {lh.gravity_pax.sum() / lh.dot_or_fill_pax.sum():.2f}; median ratio {gv.ratio_gravity_to_dot.median():.2f}")
    if args.stage <= 1:
        sys.exit()
    base = cached("baseline", lambda: stage_baseline(dem))
    iu, ju = base["iu"], base["ju"]
    code = df.code.values
    pos = {c: i for i, c in enumerate(code)}
    pax = T[iu, ju]
    if args.stage == 2:
        for a, b in [("LGA", "MIA"), ("BOS", "MIA"), ("LAX", "MIA"), ("IAH", "MIA"), ("SEA", "SAN"), ("LGA", "LAX"),
                     ("SEA", "MIA"), ("BOS", "LAX"), ("LGA", "IAH"), ("BOS", "TPA"), ("SFO", "SEA"), ("LAX", "ATL")]:
            i, j = sorted((pos[a], pos[b]))
            k = np.where((iu == i) & (ju == j))[0][0]
            print(a, b, round(base["d"][k]), "sub", round(base["t_sub"][k]), "blgc", round(base["t_bl_gc"][k]),
                  "ocean", round(base["t_ocean"][k], 1), base["ocean"][k], "ub1.7", round(base["t_ub_M1.7"][k], 1),
                  "ub2.0", round(base["t_ub_M2.0"][k], 1), "paxM", round(pax[k] / 1e6, 2))
        sys.exit()

    # ---------------------------------------------------------------- 3. addressable
    n = len(df)
    grid = cached("grid", steiner_grid)
    dfx = ext_df(df, grid)
    d = base["d"]
    t_ub_min = np.minimum(base["t_ub_M1.7"], base["t_ub_M2.0"])
    sel = np.where(base["t_sub"] - t_ub_min >= CUT_MIN)[0]          # superset of anything a shield can help
    add_rows = []
    for M in MACHS:
        t_ub = base[f"t_ub_M{M}"]
        for basis, tb in (("Boomless+ocean", base["t_bl"]), ("Boomless great circle", base["t_bl_gc"]),
                          ("subsonic", base["t_sub"])):
            s = tb - t_ub
            for thr in (10, 20, 30):
                m = s >= thr
                add_rows.append(dict(mach=M, baseline=basis, min_saving=thr, pairs=int(m.sum()),
                                     pax_M=pax[m].sum() / 1e6, pax_km_bn=(pax * d)[m].sum() / 1e9,
                                     pax_h_M=(pax * s)[m].sum() / 60 / 1e6,
                                     mean_saving_min=(pax * s)[m].sum() / max(pax[m].sum(), 1),
                                     share_of_located_pax=pax[m].sum() / pax.sum(),
                                     share_of_pax_over_1500km=pax[m].sum() / pax[d > 1500].sum()))
    add = pd.DataFrame(add_rows)
    add.to_csv(HERE / "addressable_summary.csv", index=False, float_format="%.4g")
    print(add[add.min_saving == 10].round(3).to_string())

    # ---------------------------------------------------------------- 4. spine coverage
    def spine_all():
        return {M: net_eval(dfx, base, T, SPINE, M, sel) for M in MACHS}
    sp = cached("spine_eval", spine_all)
    cov_rows = []
    pair_cols = {}
    for M in MACHS:
        ne = sp[M]
        t_ub = base[f"t_ub_M{M}"]
        ub = base["t_bl"] - t_ub
        A = ub >= CUT_MIN
        ub_s = np.where(A, ub, 0.0)
        s_net = np.zeros(len(iu))
        s_net[sel] = ne["s_bl"]
        s_sub_net = np.zeros(len(iu))
        s_sub_net[sel] = ne["s_sub"]
        ub_sub = base["t_sub"] - t_ub
        ub_sub = np.where(ub_sub >= CUT_MIN, ub_sub, 0)
        C = s_net > 0
        cov_rows.append(dict(
            mach=M, trunk_km=ne["km"], addressable_pax_M=pax[A].sum() / 1e6, captured_pax_M=pax[C].sum() / 1e6,
            share_pax=pax[C & A].sum() / pax[A].sum(),
            addressable_pax_km_bn=(pax * d)[A].sum() / 1e9, share_pax_km=(pax * d)[C & A].sum() / (pax * d)[A].sum(),
            ub_pax_h_M=(pax * ub_s).sum() / 60e6, captured_pax_h_M=(pax * s_net).sum() / 60e6,
            share_pax_h=(pax * s_net).sum() / (pax * ub_s).sum(),
            ub_pax_h_vs_sub_M=(pax * ub_sub).sum() / 60e6, captured_pax_h_vs_sub_M=(pax * s_sub_net).sum() / 60e6,
            share_pax_h_vs_sub=(pax * s_sub_net).sum() / (pax * ub_sub).sum(),
            mean_saving_captured_min=(pax * s_net).sum() / max(pax[C].sum(), 1)))
        tt = np.full(len(iu), np.nan)
        tt[sel] = ne["t"]
        pair_cols[M] = dict(t_ub=t_ub, ub=ub_s, s_net=s_net, t_net=tt)
        if M == 2.0:
            geo = {k: np.full(len(iu), np.nan) for k in ("d_on", "d_off", "L")}
            for k in geo:
                geo[k][sel] = ne[k]
    cov = pd.DataFrame(cov_rows)
    cov.to_csv(HERE / "spine_coverage.csv", index=False, float_format="%.4g")
    print(cov.round(3).T.to_string())

    # per-pair table + reasons (M2.0)
    pc = pair_cols[2.0]
    lost = pax * (pc["ub"] - pc["s_net"]).clip(min=0) / 60
    ub_gc = base["t_bl_gc"] - base["t_ub_M2.0"]
    detour = (geo["d_on"] * bt.PATH + geo["L"] + geo["d_off"] * bt.PATH) / d - 1
    reason = np.full(len(iu), "", dtype=object)
    notadd = pc["ub"] <= 0
    reason[notadd & (ub_gc < CUT_MIN)] = "ramp overhead (too short even on own shield)"
    reason[notadd & (ub_gc >= CUT_MIN)] = "ocean baseline (Boomless + M1.7 offshore already close)"
    miss = (~notadd) & (pc["s_net"] <= 0)
    reason[miss & (detour > 0.15)] = "detour (feeders + trunk > 1.15x great circle)"
    reason[miss & (detour <= 0.15) & (geo["L"] < 1500)] = "geometry: trunk run too short for ramps"
    reason[miss & (detour <= 0.15) & (geo["L"] >= 1500)] = "geometry: off-axis feeders eat the gain"
    part = (~notadd) & (pc["s_net"] > 0)
    reason[part] = "partially captured"
    P = pd.DataFrame(dict(a=code[iu], b=code[ju], distance_km=d.round(0), pax_yr=pax.round(0),
                          table1=dem["T_obs"][iu, ju] > 0,
                          t_sub=base["t_sub"].round(1), t_boomless_gc=base["t_bl_gc"].round(1),
                          t_boomless_ocean=np.where(np.isfinite(base["t_ocean"]), base["t_ocean"], np.nan).round(1),
                          t_baseline=base["t_bl"].round(1),
                          t_ub_M17=base["t_ub_M1.7"].round(1), t_ub_M20=base["t_ub_M2.0"].round(1),
                          t_spine_M17=pair_cols[1.7]["t_net"].round(1), t_spine_M20=pc["t_net"].round(1),
                          save_ub_M20=pc["ub"].round(1), save_spine_M20=pc["s_net"].round(1),
                          ub_pax_h=(pax * pc["ub"] / 60).round(0), spine_pax_h=(pax * pc["s_net"] / 60).round(0),
                          lost_pax_h=lost.round(0), detour_via_spine=detour.round(3),
                          spine_km_used=geo["L"].round(0), reason=reason))
    P = P[(P.ub_pax_h > 0) | (P.spine_pax_h > 0) | ((P.distance_km > 1500) & (P.pax_yr > 2e5))]
    P.sort_values("lost_pax_h", ascending=False).to_csv(HERE / "pairs.csv", index=False)
    top = P.sort_values("lost_pax_h", ascending=False).head(25)
    print(top[["a", "b", "distance_km", "pax_yr", "save_ub_M20", "save_spine_M20", "lost_pax_h", "detour_via_spine",
               "spine_km_used", "reason"]].to_string())
    rs = P.groupby("reason").agg(pax_M=("pax_yr", lambda x: x.sum() / 1e6), lost_pax_h_M=("lost_pax_h", lambda x: x.sum() / 1e6))
    print(rs)
    rs.to_csv(HERE / "uncaptured_by_reason.csv", float_format="%.4g")

    # named flow groups
    GROUPS = {
        "SEA-MIA": (["SEA"], ["MIA"]), "LAX-BOS": (["LAX"], ["BOS"]), "SFO-BOS": (["SFO"], ["BOS"]),
        "LAX-MIA": (["LAX"], ["MIA"]), "SEA-ATL": (["SEA"], ["ATL"]), "LAX-NYC": (["LAX"], ["LGA"]),
        "SFO-NYC": (["SFO"], ["LGA"]),
        "Texas-Northeast": (["DFW", "IAH", "AUS", "SAT"], ["LGA", "BOS", "PHL", "DCA", "BDL", "PVD", "ORH", "HVN"]),
        "Texas-California": (["DFW", "IAH", "AUS", "SAT"], ["LAX", "SFO", "SAN", "SMF"]),
        "Florida-West": (["MIA", "MCO", "TPA", "PBI", "RSW", "JAX", "SFB", "SRQ", "PGD"],
                         ["LAX", "SFO", "SAN", "SEA", "PDX", "SMF", "LAS", "PHX", "DEN", "SLC"]),
        "Pacific NW-East": (["SEA", "PDX", "PAE", "GEG", "BOI"],
                            ["LGA", "BOS", "PHL", "DCA", "ATL", "CLT", "RDU", "MIA", "MCO", "TPA", "BDL", "PIT"]),
        "Southeast-West": (["ATL", "CLT", "BNA", "RDU", "MSY"], ["LAX", "SFO", "SAN", "SEA", "PDX", "LAS", "PHX"]),
        "Northeast-Florida": (["LGA", "BOS", "PHL", "DCA", "BDL"], ["MIA", "MCO", "TPA", "PBI", "RSW", "JAX"]),
        "California-Northeast": (["LAX", "SFO", "SAN", "SMF"], ["LGA", "BOS", "PHL", "DCA", "BDL"]),
        "California-Midwest": (["LAX", "SFO", "SAN", "SMF"], ["ORD", "DTW", "MSP", "CLE", "CMH", "IND", "STL", "MCI"]),
    }
    grows = []
    for gname, (A_, B_) in GROUPS.items():
        m = (np.isin(code[iu], A_) & np.isin(code[ju], B_)) | (np.isin(code[iu], B_) & np.isin(code[ju], A_))
        r = dict(group=gname, pairs=int(m.sum()), pax_M=pax[m].sum() / 1e6,
                 mean_km=(pax * d)[m].sum() / max(pax[m].sum(), 1))
        for M in MACHS:
            pcm = pair_cols[M]
            r[f"ub_pax_h_M{M}_M"] = (pax * pcm["ub"])[m].sum() / 60e6
            r[f"spine_pax_h_M{M}_M"] = (pax * pcm["s_net"])[m].sum() / 60e6
            r[f"share_M{M}"] = r[f"spine_pax_h_M{M}_M"] / max(r[f"ub_pax_h_M{M}_M"], 1e-9)
        sub = P[(P.a.isin(A_) & P.b.isin(B_)) | (P.a.isin(B_) & P.b.isin(A_))]
        r["main_reason"] = (sub.groupby("reason").lost_pax_h.sum().idxmax() if len(sub) and sub.lost_pax_h.sum() > 0 else "")
        grows.append(r)
    gdf = pd.DataFrame(grows)
    gdf.to_csv(HERE / "flow_groups.csv", index=False, float_format="%.4g")
    print(gdf.round(3).to_string())
    if args.stage <= 4:
        sys.exit()

    # ---------------------------------------------------------------- 5. expansions
    def scen_all():
        out = []
        for name, edges in SCEN.items():
            for M in MACHS:
                ne = net_eval(dfx, base, T, edges, M, sel)
                ec, _ = econ(ne, M)
                out.append(dict(scenario=name, mach=M, trunk_km=ne["km"], pax_using_M=ne["pax"][ne["s_bl"] > 0].sum() / 1e6,
                                pax_h_vs_bl_M=ne["B_bl"] / 1e6, pax_h_vs_sub_M=ne["B_sub"] / 1e6,
                                share_of_ub_pax_h=ne["B_bl"] / ((pax * np.where(base["t_bl"] - base[f"t_ub_M{M}"] >= CUT_MIN, base["t_bl"] - base[f"t_ub_M{M}"], 0)).sum() / 60),
                                pax_h_per_km=ne["B_bl"] / ne["km"], **ec))
                print(f"  {name} M{M}: {ne['km']:.0f} km, {ne['B_bl'] / 1e6:.1f} M pax-h, B/C {ec['bcr_bl_cap100_nofatigue']:.4f}")
        return out
    scen = pd.DataFrame(cached("scenarios", scen_all))
    scen.to_csv(HERE / "scenarios.csv", index=False, float_format="%.4g")

    node_mass = T.sum(1)
    ok = np.r_[node_mass >= 1.0e6, np.ones(len(grid), bool)]
    G = candidate_graph(dfx, n, ok)
    cand = [i for i in range(n) if node_mass[i] >= 1.5e6 and i in G]
    print(f"candidate graph: {G.number_of_nodes()} nodes, {G.number_of_edges()} edges; {len(cand)} candidate endpoints")
    import re
    logf = CACHE / "run5_part1.log"
    resume = re.findall(r"spine M2\.0 step \d+: \+(\S+) km=", logf.read_text()) if logf.exists() else []
    growth = cached("greedy_spine_M2.0", lambda: greedy(dfx, base, T, 2.0, sel, SPINE, cand, G, tag="spine",
                                                        resume=resume))

    def scratch():
        code_l = list(dfx.code)
        ubv = base["t_bl"] - base["t_ub_M2.0"]
        w = pax * np.where(ubv >= CUT_MIN, ubv, 0)
        best = None
        for k in np.argsort(-w)[:15]:
            if iu[k] not in G or ju[k] not in G:
                continue
            p = nx.shortest_path(G, int(iu[k]), int(ju[k]), weight="length")
            e = [(code_l[a], code_l[b]) for a, b in zip(p[:-1], p[1:])]
            ne = net_eval(dfx, base, T, e, 2.0, sel)
            if best is None or ne["B_bl"] / 1 > best[0]:
                best = (ne["B_bl"], e)
        return greedy(dfx, base, T, 2.0, sel, best[1], cand, G, tag="scratch", max_steps=14, km_max=11000)
    growth0 = cached("greedy_scratch_M2.0", scratch)

    def curve(rows, label):
        out = []
        for r in rows:
            for M in MACHS:
                ne = net_eval(dfx, base, T, r["edges"], M, sel)
                ec, _ = econ(ne, M)
                out.append(dict(path=label, step=r["step"], added=r["added"], mach=M, trunk_km=ne["km"],
                                pax_using_M=ne["pax"][ne["s_bl"] > 0].sum() / 1e6, pax_h_vs_bl_M=ne["B_bl"] / 1e6,
                                pax_h_vs_sub_M=ne["B_sub"] / 1e6, marginal_pax_h_per_km=r["marginal_per_km"],
                                avg_pax_h_per_km=ne["B_bl"] / ne["km"], **ec,
                                edges=" ".join(f"{a}-{b}" for a, b in r["edges"])))
        return out
    cc = pd.DataFrame(cached("curve_spine", lambda: curve(growth, "greedy from spine"))
                      + cached("curve_scratch", lambda: curve(growth0, "greedy from scratch")))
    for M in MACHS:
        ubp = (pax * np.where(base["t_bl"] - base[f"t_ub_M{M}"] >= CUT_MIN, base["t_bl"] - base[f"t_ub_M{M}"], 0)).sum() / 60e6
        cc.loc[cc.mach == M, "share_of_ub_pax_h"] = cc.loc[cc.mach == M, "pax_h_vs_bl_M"] / ubp
    cc.to_csv(HERE / "coverage_curve.csv", index=False, float_format="%.4g")
    print(cc[cc.mach == 2.0].drop(columns="edges").round(4).to_string())

    # ---------------------------------------------------------------- knee
    def knee(c):
        """Kneedle on cumulative pax-h vs km with (0, 0) prepended; plus B/C-max and 0.5x-average rule."""
        x = np.r_[0, c.trunk_km.values]
        y = np.r_[0, c.pax_h_vs_bl_M.values]
        xn, yn = x / x.max(), y / y.max()
        k = int(np.argmax(yn - xn)) - 1
        half = c.index[(c.marginal_pax_h_per_km < 0.5 * c.avg_pax_h_per_km.shift(1).fillna(np.inf))]
        return max(k, 0), (c.index.get_loc(half[0]) - 1 if len(half) else len(c) - 1)
    ks = cc[(cc.path == "greedy from spine") & (cc.mach == 2.0)].reset_index(drop=True)
    k_knee, k_half = knee(ks)
    # 'shoulder' (map network): last step before the marginal pax-h/km first drops below 0.3x the spine's average
    sp_avg = ks.avg_pax_h_per_km.iloc[0]
    low = ks.index[ks.marginal_pax_h_per_km < 0.3 * sp_avg]
    k_map = int(low[0] - 1) if len(low) else len(ks) - 1
    print(f"knee (kneedle) step {k_knee} at {ks.trunk_km[k_knee]:.0f} km; 0.5x-average rule stops at step {max(k_half, 0)}; "
          f"shoulder/map network step {k_map} at {ks.trunk_km[k_map]:.0f} km")
    json.dump(dict(knee_step=k_knee, knee_km=float(ks.trunk_km[k_knee]), half_rule_step=max(k_half, 0),
                   shoulder_step=k_map, shoulder_km=float(ks.trunk_km[k_map]), shoulder_edges=ks.edges[k_map]),
              open(HERE / "knee.json", "w"), indent=1)
    if args.stage <= 5:
        sys.exit()

    # ---------------------------------------------------------------- 6. figures
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D
    FIG.mkdir(exist_ok=True)
    C1, C2, C3, C4, INK, INK2, GRID = "#2a78d6", "#eb6834", "#1baf7a", "#6250d6", "#0b0b0b", "#52514e", "#e4e3df"
    plt.rcParams.update({"font.size": 9, "axes.edgecolor": INK2, "axes.labelcolor": INK, "xtick.color": INK2,
                         "ytick.color": INK2, "font.family": "DejaVu Sans"})

    # ---- map
    exp_edges = [tuple(e.split("-")) for e in ks.edges[k_map].split()]
    nets = {"spine": SPINE, "exp": [e for e in exp_edges if e not in SPINE and e[::-1] not in SPINE]}
    ne_full = net_eval(dfx, base, T, exp_edges, 2.0, sel)
    _, lb_full = econ(ne_full, 2.0)
    ne_sp = sp[2.0]
    _, lb_sp = econ(ne_sp, 2.0)
    fig, ax = plt.subplots(figsize=(11, 7.0), dpi=150)
    ext = (-2600.0, 2400.0, -1750.0, 1650.0)
    bt._basemap(ax, ext)
    # offshore M1.7 tracks
    for name, pts in OCEAN.items():
        x, y = bt._albers([p[1] for p in pts], [p[2] for p in pts])
        ax.plot(x, y, color="#8fb3de", lw=1.2, ls=(0, (1, 1.6)), zorder=2)
    trf = ne_full["tr"]
    spine_codes = {c for e in SPINE for c in e}
    def seg_is_spine(tr, u, v):
        ei = tr.G[u][v]["edge"]
        a, b = tr.edges[ei]
        return (a, b) in SPINE or (b, a) in SPINE
    WS = 0.12   # line width per M pax/yr
    for (u, v), (p, ph) in lb_full.items():
        x, y = bt._albers(trf.ll[[u, v], 0], trf.ll[[u, v], 1])
        col = C1 if seg_is_spine(trf, u, v) else C2
        ax.plot(x, y, color=col, lw=max(0.8, p / 1e6 * WS), solid_capstyle="round", zorder=5)
    for u, v in trf.G.edges():
        if (min(u, v), max(u, v)) not in lb_full:
            x, y = bt._albers(trf.ll[[u, v], 0], trf.ll[[u, v], 1])
            ax.plot(x, y, color=C2, lw=0.8, zorder=4)
    # top uncaptured flows (by the spine), dashed great-circle arcs
    Pt = pd.read_csv(HERE / "pairs.csv").head(12)
    for r in Pt.itertuples():
        ia, ib = pos[r.a], pos[r.b]
        pa, pb = (df.lat[ia], df.lon[ia]), (df.lat[ib], df.lon[ib])
        arc = np.array([bt.gc_interp(pa, pb, f) for f in np.linspace(0, 1, 40)])
        x, y = bt._albers(arc[:, 0], arc[:, 1])
        ax.plot(x, y, color=C4, lw=0.6 + r.lost_pax_h / 1e6 * 2.2, ls=(0, (4, 2.5)), alpha=0.85, zorder=3)
    lab = set(spine_codes) | {c for e in exp_edges for c in e} | set(Pt.a) | set(Pt.b)
    is_st = lambda c: len(c) == 3 and c[0] == "S" and c[1:].isdigit()
    for c in lab:
        if is_st(c):
            continue
        i = pos[c]
        x, y = bt._albers(df.lat[i], df.lon[i])
        ax.plot(x, y, "o", ms=3.2, color=INK, zorder=6)
        ax.annotate(c, (x, y), xytext=(3, 3), textcoords="offset points", fontsize=7.5, color=INK, zorder=7)
    for c in {c for e in exp_edges for c in e if is_st(c)}:
        i = list(dfx.code).index(c)
        x, y = bt._albers(dfx.lat[i], dfx.lon[i])
        ax.plot(x, y, "s", ms=2.5, color=INK2, zorder=6)
    km_sp = ne_sp["km"]
    km_ex = ne_full["km"]
    leg = [Line2D([], [], color=C1, lw=4, label=f"Current spine, {km_sp:,.0f} km (width = pax/yr through each km)"),
           Line2D([], [], color=C2, lw=3, label=f"Greedy additions to {km_ex:,.0f} km"),
           Line2D([], [], color=C4, lw=1.6, ls=(0, (4, 2.5)), label="Top 12 flows the spine misses (width = pax-h lost)"),
           Line2D([], [], color="#8fb3de", lw=1.4, ls=(0, (1, 1.6)), label="Offshore M1.7 tracks (Boomless baseline)"),
           ]
    for w, lab_ in ((10, "10 M"), (40, "40 M"), (80, "80 M pax/yr")):
        leg.append(Line2D([], [], color=C1, lw=w * WS, label=lab_))
    ax.legend(handles=leg, loc="lower left", fontsize=7.5, frameon=False)
    ax.set_title(f"BoomShield at Mach 2.0: the spine and its greedy expansion (100 % capture, DOT 2025 O&D)",
                 fontsize=10.5, color=INK, loc="left")
    fig.tight_layout()
    fig.savefig(FIG / "coverage_map.png", facecolor="white")
    plt.close(fig)

    # ---- coverage curve
    sc = scen
    fig, (a1, a2) = plt.subplots(1, 2, figsize=(11, 4.4), dpi=150)
    for M, col in ((2.0, C1), (1.7, C2)):
        c = cc[(cc.path == "greedy from spine") & (cc.mach == M)]
        a1.plot(np.r_[0, c.trunk_km] / 1e3, np.r_[0, c.pax_h_vs_bl_M], "-o", color=col, lw=2, ms=4,
                label=f"Greedy from spine, M{M}")
        ubp = add[(add.mach == M) & (add.baseline == "Boomless+ocean") & (add.min_saving == 10)].pax_h_M.iloc[0]
        a1.axhline(ubp, color=col, lw=1, ls=(0, (4, 3)))
        a1.text(0.3, ubp + 2, f"Upper bound M{M}: own-great-circle shields, {ubp:.0f} M pax-h", color=INK2, fontsize=7.5)
        c0 = cc[(cc.path == "greedy from scratch") & (cc.mach == M)]
        if M == 2.0:
            a1.plot(np.r_[0, c0.trunk_km] / 1e3, np.r_[0, c0.pax_h_vs_bl_M], "-", color=C3, lw=1.5,
                    label="Greedy from scratch, M2.0")
        s2 = sc[sc.mach == M]
        a1.plot(s2.trunk_km / 1e3, s2.pax_h_vs_bl_M, "D", ms=5, mfc="white", mec=col, mew=1.4, ls="none")
        a2.plot(c.trunk_km / 1e3, c.bcr_bl_cap100_nofatigue, "-o", color=col, lw=2, ms=4,
                label=f"M{M}, no fatigue cap")
        a2.plot(c.trunk_km / 1e3, c.bcr_bl_cap100_lidcapped, "--", color=col, lw=1.5,
                label=f"M{M}, lid cap 100 flights/day")
        if M == 2.0:
            a2.plot(c.trunk_km / 1e3, c.bcr_bl_cap100_lid1000, ":", color=col, lw=1.8,
                    label="M2.0, lid cap 1,000 flights/day (fatigue/ revised low)")
        a2.plot(s2.trunk_km / 1e3, s2.bcr_bl_cap100_nofatigue, "D", ms=5, mfc="white", mec=col, mew=1.4, ls="none")
    kk = ks.iloc[k_knee]
    a1.annotate(f"knee: {kk.trunk_km:,.0f} km", (kk.trunk_km / 1e3, kk.pax_h_vs_bl_M), xytext=(40, -45),
                textcoords="offset points", fontsize=8, color=INK, arrowprops=dict(arrowstyle="-", color=INK2, lw=0.8))
    a1.plot([], [], "D", ms=5, mfc="white", mec=INK2, ls="none", label="Named scenarios")
    a1.set_xlabel("Total shield length, 1,000 km")
    a1.set_ylabel("Pax-hours saved per year vs Boomless, M (100 % capture)")
    a1.set_xlim(0, None)
    a1.set_ylim(0, None)
    a1.legend(fontsize=7.5, frameon=False, loc="lower right")
    a2.set_yscale("log")
    a2.set_xlabel("Total shield length, 1,000 km")
    a2.set_ylabel("B/C vs Boomless (100 % capture, base cost, log scale)")
    a2.text(0.99, 0.01, "10,000 flights/day (fatigue/ revised base) = no cap at these loads", transform=a2.transAxes,
            ha="right", va="bottom", fontsize=7, color=INK2)
    a2.legend(fontsize=7.5, frameon=False, loc="center left", bbox_to_anchor=(0.0, 0.42))
    for a in (a1, a2):
        a.grid(color=GRID, lw=0.6)
        a.spines[["top", "right"]].set_visible(False)
    a1.set_title("Captured time savings vs network size", fontsize=10, loc="left", color=INK)
    a2.set_title("Benefit/cost vs network size", fontsize=10, loc="left", color=INK)
    fig.tight_layout()
    fig.savefig(FIG / "coverage_curve.png", facecolor="white")
    plt.close(fig)
    print("figures written")
