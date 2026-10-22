#!/usr/bin/env python3
"""
BoomShield shared cross-country trunk: network analysis (screening level).

Thesis tested: a shared trunk carrying ALL long-haul North American O&D flows that
cross the country (feeders flown conventionally to/from on-ramps) raises traffic per
shield-km far above the pairwise route economics in ../routes and ../../economics.

Pipeline
  1. O&D demand matrix among 57 North American city markets: gravity model calibrated
     to the OAG 2025 airport-pair seats in ../routes/routes.csv (BTS DB1B / T-100 and
     DOT Consumer Airfare tables were unreachable from this sandbox).
  2. Door-to-door time per pair: subsonic direct, Boomless (M1.15) direct, and via a
     trunk tree (feeder at M0.85 or M1.15 to an on-ramp, accelerate on the shield to
     Mach M, cruise, decelerate, feeder off). Same-aircraft (merge penalty) or
     hub-transfer (separate flights + connection penalty).
  3. Trunk trees: (a) the user's heuristic tree; (b) a greedy tree per Mach grown on a
     Delaunay graph of the metros, maximising pax-hours saved vs Boomless per km.
  4. Per-segment traffic, pax-h saved per km, capture 10/30/100 %, fare premium.
  5. Economics with ../../economics/cost_model.py (capex/opex per km, per-Mach design).
  6. On-ramp acceleration zones (Mach-onset focus) for the focus-boom sibling study.
  7. Map figure.

Usage: python3 build_trunk.py [--quick]     (~15 min cold, ~6 min with cached trees; one process)
Outputs go next to this file; intermediate results are cached in .cache/ so a
restarted run resumes.
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
BS = HERE.parent.parent
CACHE = HERE / ".cache"
FIG = HERE / "figures"
R_E = 6371.0

# =====================================================================================
# 1. City markets.  Coordinates = principal airport.  mass = domestic/North-American
#    O&D enplanements (originating passengers, excludes connections), M/yr, ~2024.
#    Author estimates from airport totals x O&D share x domestic share; UNVERIFIED.
# =====================================================================================
METROS = [
    # code, name, lat, lon, country, mass (M O&D enpl./yr), airports in market
    ("NYC", "New York", 40.777, -73.874, "US", 42.0, "JFK LGA EWR HPN ISP"),
    ("LAX", "Los Angeles", 33.942, -118.409, "US", 36.0, "LAX SNA BUR ONT LGB"),
    ("SFO", "SF Bay Area", 37.621, -122.379, "US", 25.0, "SFO OAK SJC"),
    ("CHI", "Chicago", 41.974, -87.907, "US", 21.5, "ORD MDW"),
    ("WAS", "Washington/Baltimore", 38.852, -77.037, "US", 27.0, "DCA IAD BWI"),
    ("LAS", "Las Vegas", 36.084, -115.154, "US", 25.0, "LAS"),
    ("FLL", "South Florida", 25.796, -80.287, "US", 25.0, "MIA FLL PBI"),
    ("MCO", "Orlando", 28.431, -81.308, "US", 22.0, "MCO SFB"),
    ("DEN", "Denver", 39.856, -104.674, "US", 21.0, "DEN"),
    ("DFW", "Dallas/Fort Worth", 32.900, -97.040, "US", 19.0, "DFW DAL"),
    ("PHX", "Phoenix", 33.434, -112.012, "US", 17.5, "PHX AZA"),
    ("ATL", "Atlanta", 33.641, -84.428, "US", 16.0, "ATL"),
    ("SEA", "Seattle", 47.450, -122.309, "US", 16.0, "SEA PAE"),
    ("BOS", "Boston", 42.366, -71.010, "US", 15.0, "BOS"),
    ("HOU", "Houston", 29.990, -95.337, "US", 14.5, "IAH HOU"),
    ("SAN", "San Diego", 32.734, -117.190, "US", 12.0, "SAN"),
    ("TPA", "Tampa", 27.976, -82.533, "US", 11.5, "TPA"),
    ("BNA", "Nashville", 36.124, -86.678, "US", 11.0, "BNA"),
    ("AUS", "Austin", 30.197, -97.666, "US", 10.5, "AUS"),
    ("MSP", "Minneapolis", 44.885, -93.222, "US", 9.5, "MSP"),
    ("PHL", "Philadelphia", 39.874, -75.243, "US", 8.5, "PHL"),
    ("DTW", "Detroit", 42.212, -83.353, "US", 8.0, "DTW"),
    ("PDX", "Portland", 45.589, -122.597, "US", 7.8, "PDX"),
    ("CLT", "Charlotte", 35.214, -80.943, "US", 7.0, "CLT"),
    ("RDU", "Raleigh-Durham", 35.878, -78.788, "US", 7.3, "RDU"),
    ("STL", "St. Louis", 38.749, -90.370, "US", 6.8, "STL"),
    ("SMF", "Sacramento", 38.695, -121.591, "US", 6.7, "SMF"),
    ("MSY", "New Orleans", 29.993, -90.258, "US", 6.6, "MSY"),
    ("SLC", "Salt Lake City", 40.790, -111.979, "US", 6.5, "SLC"),
    ("MCI", "Kansas City", 39.298, -94.714, "US", 5.4, "MCI"),
    ("RSW", "Fort Myers", 26.536, -81.755, "US", 5.2, "RSW"),
    ("SAT", "San Antonio", 29.534, -98.470, "US", 5.0, "SAT"),
    ("CLE", "Cleveland", 41.412, -81.850, "US", 4.9, "CLE"),
    ("PIT", "Pittsburgh", 40.492, -80.233, "US", 4.7, "PIT"),
    ("IND", "Indianapolis", 39.717, -86.294, "US", 4.6, "IND"),
    ("CMH", "Columbus", 39.998, -82.892, "US", 4.4, "CMH"),
    ("CVG", "Cincinnati", 39.049, -84.668, "US", 4.0, "CVG"),
    ("JAX", "Jacksonville", 30.494, -81.688, "US", 3.5, "JAX"),
    ("MKE", "Milwaukee", 42.947, -87.897, "US", 3.4, "MKE"),
    ("BDL", "Hartford", 41.939, -72.683, "US", 3.3, "BDL"),
    ("ABQ", "Albuquerque", 35.040, -106.609, "US", 2.6, "ABQ"),
    ("OMA", "Omaha", 41.303, -95.894, "US", 2.6, "OMA"),
    ("BUF", "Buffalo", 42.940, -78.732, "US", 2.5, "BUF"),
    ("BOI", "Boise", 43.564, -116.223, "US", 2.5, "BOI"),
    ("RIC", "Richmond", 37.505, -77.320, "US", 2.4, "RIC"),
    ("MEM", "Memphis", 35.042, -89.977, "US", 2.4, "MEM"),
    ("OKC", "Oklahoma City", 35.393, -97.601, "US", 2.3, "OKC"),
    ("RNO", "Reno", 39.499, -119.768, "US", 2.3, "RNO"),
    ("ELP", "El Paso", 31.807, -106.378, "US", 2.0, "ELP"),
    ("TUS", "Tucson", 32.116, -110.941, "US", 1.9, "TUS"),
    # Canada / Mexico: mass = all North-American O&D (domestic + transborder);
    # cross-border pairs are scaled by BORDER below.
    ("YYZ", "Toronto", 43.678, -79.625, "CA", 13.0, "YYZ YTZ"),
    ("YVR", "Vancouver", 49.197, -123.182, "CA", 8.0, "YVR"),
    ("YUL", "Montreal", 45.470, -73.741, "CA", 6.5, "YUL"),
    ("YYC", "Calgary", 51.131, -114.010, "CA", 5.5, "YYC"),
    ("YOW", "Ottawa", 45.322, -75.667, "CA", 2.4, "YOW"),
    ("MEX", "Mexico City", 19.436, -99.072, "MX", 14.0, "MEX NLU"),
]
BORDER = {frozenset(("US", "CA")): 0.35, frozenset(("US", "MX")): 0.25,
          frozenset(("CA", "MX")): 0.15}            # border/visa/currency deterrence [U]
US_DOMESTIC_ENPL_2024 = 0.87e9   # BTS monthly releases 2024 summed (Jun-Dec quoted, Jan-May est.) [U]
ENPL_PER_OD_TRIP = 1.33          # connections: ~1.33 enplanements per O&D one-way trip [U]
US_OD_TRIPS = US_DOMESTIC_ENPL_2024 / ENPL_PER_OD_TRIP

# Calibration anchors: OAG 2025 US airport-pair seats (routes.csv source OAG25US),
# converted to city-market O&D pax = seats x LF x local share / airport share-of-market.
# local share = fraction of seats carrying local (non-connecting) O&D pax  [U]
# airport share = this airport pair's share of all nonstop seats in the city market [U]
LF_ANCHOR = 0.83
D0_GROUND = 400.0   # km: below ~d0 car/rail take the market (fixed; short anchors cannot fit it)
LOW_DEMAND = 0.75   # 'low' case: scales T so each metro sends ~65 % (not ~85 %) of its O&D to the 56 others
ANCHORS = [
    # metro_a, metro_b, airport pair, OAG seats, local share, airport share
    ("NYC", "LAX", "JFK-LAX", 3431159, 0.70, 0.50),
    ("NYC", "CHI", "LGA-ORD", 3334011, 0.65, 0.50),
    ("LAX", "SFO", "LAX-SFO", 3306223, 0.70, 0.35),
    ("LAS", "LAX", "LAS-LAX", 3283510, 0.85, 0.45),
    ("ATL", "MCO", "ATL-MCO", 3226190, 0.40, 1.00),
    ("DEN", "PHX", "DEN-PHX", 3118880, 0.40, 1.00),
    ("LAX", "CHI", "LAX-ORD", 2800000, 0.55, 0.70),
    ("ATL", "NYC", "ATL-LGA", 2763776, 0.55, 0.55),
    ("ATL", "FLL", "ATL-FLL", 2600000, 0.45, 0.45),
]

# =====================================================================================
# 2. Time model (same constants as ../routes/build_routes.py; Overture-class)
# =====================================================================================
A_MS = 295.1                     # m/s, isothermal stratosphere
A_KMH = A_MS * 3.6
V085, V094, V115 = 0.85 * A_KMH, 0.94 * A_KMH, 1.15 * A_KMH
TAXI = 22.0                      # min, taxi out + in
APPR = (10.0, 60.0)              # min, km terminal approach
PATH = 1.03                      # flown / great-circle on conventional legs
CLIMB = (18.0, 230.0)            # to FL350-390
ACC_BL = (7.0, 130.0)            # M0.9 -> 1.15
DESC_SUB = (18.0, 200.0)
DEC_BL = (4.0, 75.0)             # M1.15 -> subsonic before descent
FL600_EXTRA = (4.0, 60.0)        # extra descent from shield flight level (routes: 22/260 vs 18/200)
A_ACC = 0.23                     # m/s^2 mean along-track accel (fits M0.95->1.7 in 16 min/377 km)
A_DEC = 0.46                     # m/s^2 (fits 1.7->0.95 in 8 min/188 km)
V_SUB_START = 0.95 * A_MS        # accel onto the trunk starts from M0.95 (sub feeder)
V_BL_START = 1.15 * A_MS


def t_direct(d_gc, mode):
    """Block minutes, all-overland direct flight, as routes.md (sub / bl / ov094)."""
    L = np.asarray(d_gc, float) * PATH
    sub = np.where(L >= 490, TAXI + 46 + (L - 490) / V085 * 60, TAXI + 46 * L / 490)
    if mode == "sub":
        return sub
    o94 = np.where(L >= 490, TAXI + 46 + (L - 490) / V094 * 60, sub)
    bl = np.where(L >= 695, TAXI + 57 + (L - 695) / V115 * 60, np.inf)
    return np.minimum(np.minimum(bl, o94), sub)


def t_trunk(d_on, L_tr, d_off, mach, feeder, merge_min=5.0, a_acc=A_ACC, a_dec=A_DEC):
    """Door-to-door block minutes for ONE aircraft: feeder leg (great-circle d_on km) to
    the on-ramp, L_tr km along the trunk (accelerate, cruise at `mach`, decelerate, all
    over the shield), feeder leg d_off km to destination.  Peak Mach falls below
    `mach` if the trunk run is too short for the ramps.  Returns inf where the trunk
    gives no supersonic run.  Arrays broadcast."""
    d_on, L_tr, d_off = (np.asarray(x, float) for x in (d_on, L_tr, d_off))
    if feeder == "sub":
        s1, t_climb, v_f, v0 = CLIMB[1], CLIMB[0], V085, V_SUB_START
        desc_d, desc_t = DESC_SUB[1], DESC_SUB[0]
    else:
        s1, t_climb, v_f, v0 = CLIMB[1] + ACC_BL[1], CLIMB[0] + ACC_BL[0], V115, V_BL_START
        desc_d, desc_t = DESC_SUB[1] + DEC_BL[1], DESC_SUB[0] + DEC_BL[0]
    s_in = d_on * PATH
    s_out = s_in + L_tr
    P = s_out + d_off * PATH
    s_fdesc = P - APPR[1] - desc_d                      # start of final descent
    s0 = np.maximum(s1, s_in)                           # accel can start
    s9 = np.minimum(s_out, s_fdesc)                     # decel must end
    A = (s9 - s0 - FL600_EXTRA[1]) * 1e3                # m available for accel+decel+cruise
    VM = mach * A_MS
    k = 1 / (2 * a_acc) + 1 / (2 * a_dec)
    need = (VM ** 2 - v0 ** 2) * k
    Vp = np.where(A >= need, VM, np.sqrt(np.maximum(v0 ** 2 + np.maximum(A, 0) / k, 0)))
    cruise = np.maximum(A - (Vp ** 2 - v0 ** 2) * k, 0)  # m
    t_sup = ((Vp - v0) / a_acc + (Vp - v0) / a_dec + np.where(Vp > 0, cruise / np.maximum(Vp, 1), 0)) / 60
    feed_km = np.maximum(s0 - s1, 0) + np.maximum(s_fdesc - s9, 0)
    t = (TAXI + t_climb + feed_km / v_f * 60 + t_sup + FL600_EXTRA[0] + desc_t + APPR[0]
         + merge_min)
    ok = (A > 0) & (Vp > v0 * 1.02) & (s_fdesc >= s1)
    return np.where(ok, t, np.inf), np.where(ok, Vp / A_MS, 0.0)


def accel_km(m_from, m_to, a=A_ACC):
    return ((m_to * A_MS) ** 2 - (m_from * A_MS) ** 2) / (2 * a) / 1e3


# =====================================================================================
# Geometry helpers
# =====================================================================================
def xyz(lat, lon):
    la, lo = np.radians(lat), np.radians(lon)
    return np.stack([np.cos(la) * np.cos(lo), np.cos(la) * np.sin(lo), np.sin(la)], -1)


def gc(lat1, lon1, lat2, lon2):
    a, b = xyz(lat1, lon1), xyz(lat2, lon2)
    return R_E * 2 * np.arcsin(np.clip(np.linalg.norm(a - b, axis=-1) / 2, 0, 1))


def gc_interp(p, q, f):
    a, b = xyz(*p), xyz(*q)
    om = math.acos(float(np.clip(a @ b, -1, 1)))
    if om < 1e-9:
        return p
    v = (math.sin((1 - f) * om) * a + math.sin(f * om) * b) / math.sin(om)
    return float(np.degrees(np.arcsin(v[2]))), float(np.degrees(np.arctan2(v[1], v[0])))


# =====================================================================================
# 3. Demand: gravity model  T_ij = k m_i m_j d^-beta (1 - exp(-(d/d0)^2)) b_ij
# =====================================================================================
def metros_df():
    df = pd.DataFrame(METROS, columns=["code", "name", "lat", "lon", "country", "mass_M", "airports"])
    return df


def dist_matrix(df):
    la, lo = df.lat.values, df.lon.values
    return gc(la[:, None], lo[:, None], la[None, :], lo[None, :])


def deterrence(d, beta, d0):
    d = np.maximum(d, 1.0)
    return d ** (-beta) * (1 - np.exp(-(d / d0) ** 2))


def border_matrix(df):
    c = df.country.values
    B = np.ones((len(c), len(c)))
    for i in range(len(c)):
        for j in range(len(c)):
            if c[i] != c[j]:
                B[i, j] = BORDER[frozenset((c[i], c[j]))]
    return B


def calibrate(df, D):
    idx = {c: i for i, c in enumerate(df.code)}
    m = df.mass_M.values * 1e6
    rows = []
    for a, b, ap, seats, loc, share in ANCHORS:
        rows.append(dict(metro_a=a, metro_b=b, airport_pair=ap, oag_seats_2025=seats, local_share=loc,
                         airport_share=share, od_pax_est=seats * LF_ANCHOR * loc / share,
                         i=idx[a], j=idx[b]))
    an = pd.DataFrame(rows)
    best = None
    for beta in np.arange(0.0, 2.01, 0.05):
        for d0 in (D0_GROUND,):
            g = m[an.i] * m[an.j] * deterrence(D[an.i, an.j], beta, d0)
            lk = np.mean(np.log(an.od_pax_est) - np.log(g))
            res = np.log(an.od_pax_est) - np.log(g) - lk
            sse = float((res ** 2).sum())
            if best is None or sse < best[0]:
                best = (sse, beta, d0, math.exp(lk))
    sse, beta, d0, k = best
    g = k * m[an.i] * m[an.j] * deterrence(D[an.i, an.j], beta, d0)
    an["model_pax"] = g.round(0)
    an["ratio_model_to_anchor"] = (g / an.od_pax_est).round(2)
    an["distance_km"] = D[an.i, an.j].round(0)
    return dict(beta=beta, d0=d0, k=k, rms_log=math.sqrt(sse / len(an))), an.drop(columns=["i", "j"])


def demand(df, D, cal, beta=None):
    m = df.mass_M.values * 1e6
    beta = cal["beta"] if beta is None else beta
    k = cal["k"]
    if beta != cal["beta"]:   # re-scale k so the anchor total is preserved
        an = [(df.code.tolist().index(a), df.code.tolist().index(b)) for a, b, *_ in ANCHORS]
        ii, jj = np.array(an).T
        g0 = (m[ii] * m[jj] * deterrence(D[ii, jj], cal["beta"], cal["d0"])).sum() * cal["k"]
        g1 = (m[ii] * m[jj] * deterrence(D[ii, jj], beta, cal["d0"])).sum()
        k = g0 / g1
    T = k * m[:, None] * m[None, :] * deterrence(D, beta, cal["d0"]) * border_matrix(df)
    np.fill_diagonal(T, 0)
    return T   # annual passengers, both directions, per unordered pair (symmetric matrix)


# =====================================================================================
# 4. Trunk trees and pair routing
# =====================================================================================
class Trunk:
    """A tree of metro nodes (straight great-circle edges), subdivided every ~step km
    into ramp points.  Pairs join/leave at any ramp point (same aircraft) or at a
    metro node on the tree (hub transfer)."""

    def __init__(self, df, edges, step=100.0):
        self.df = df
        self.code = list(df.code)
        self.edges = [tuple(e) for e in edges]
        self.G = nx.Graph()
        pts = {}   # (lat, lon) -> id
        self.pt_ll = []
        self.pt_node = []   # metro index or -1

        def pid(i_or_ll, node=-1):
            key = i_or_ll
            if key not in pts:
                pts[key] = len(self.pt_ll)
                self.pt_ll.append(key)
                self.pt_node.append(node)
            return pts[key]
        self.seg = []   # (edge index, p, q, length)
        for ei, (a, b) in enumerate(self.edges):
            ia, ib = self.code.index(a), self.code.index(b)
            pa = (df.lat[ia], df.lon[ia])
            pb = (df.lat[ib], df.lon[ib])
            L = float(gc(*pa, *pb))
            n = max(1, int(round(L / step)))
            prev = pid(("N", ia), ia)
            for s in range(1, n + 1):
                if s == n:
                    cur = pid(("N", ib), ib)
                else:
                    cur = pid(("E", ei, s))
                self.G.add_edge(prev, cur, length=L / n, edge=ei)
                prev = cur
        # coordinates for each point
        ll = []
        for key in self.pt_ll:
            if key[0] == "N":
                ll.append((df.lat[key[1]], df.lon[key[1]]))
            else:
                _, ei, s = key
                a, b = self.edges[ei]
                ia, ib = self.code.index(a), self.code.index(b)
                L = float(gc(df.lat[ia], df.lon[ia], df.lat[ib], df.lon[ib]))
                n = max(1, int(round(L / step)))
                ll.append(gc_interp((df.lat[ia], df.lon[ia]), (df.lat[ib], df.lon[ib]), s / n))
        self.ll = np.array(ll)
        self.pt_node = np.array(self.pt_node)
        self.km = sum(d["length"] for *_, d in self.G.edges(data=True))
        n = len(self.ll)
        self.Dt = np.zeros((n, n))
        for s, lens in nx.all_pairs_dijkstra_path_length(self.G, weight="length"):
            for t, v in lens.items():
                self.Dt[s, t] = v
        # metro -> ramp point distances
        self.Dm = gc(df.lat.values[:, None], df.lon.values[:, None], self.ll[None, :, 0], self.ll[None, :, 1])
        self.node_pts = np.where(self.pt_node >= 0)[0]

    def nearest(self, R):
        R = min(R, len(self.ll))
        idx = np.argsort(self.Dm, axis=1)[:, :R]
        return idx, np.take_along_axis(self.Dm, idx, 1)


def evaluate(trunk, D, T, mach, cfg, pairs=None, detail=False):
    """Route every pair over the trunk at `mach`.  Returns per-pair best trunk time,
    savings vs subsonic and Boomless direct, and (if detail) ramp choices."""
    n = len(D)
    if pairs is None:
        iu, ju = np.triu_indices(n, 1)
    else:
        iu, ju = pairs
    pax = T[iu, ju]
    d = D[iu, ju]
    t_sub = t_direct(d, "sub")
    t_bl = t_direct(d, "bl")
    feeders = ("sub",) if cfg["feeder"] == "sub" else ("sub", "bl")
    if cfg["connection"] == "same":
        idx, dist = trunk.nearest(cfg["R"])
        a_idx, a_d = idx[iu][:, :, None], dist[iu][:, :, None]
        b_idx, b_d = idx[ju][:, None, :], dist[ju][:, None, :]
        L = trunk.Dt[a_idx, b_idx]
        best = np.full(L.shape, np.inf)
        mp = np.zeros(L.shape)
        for f in feeders:
            t, m_ = t_trunk(a_d, L, b_d, mach, f, cfg["merge_min"], cfg["a_acc"], cfg["a_dec"])
            upd = t < best
            best = np.where(upd, t, best)
            mp = np.where(upd, m_, mp)
        flat = best.reshape(len(iu), -1)
        k = np.argmin(flat, 1)
        t_tr = flat[np.arange(len(iu)), k]
        R = idx.shape[1]
        on = idx[iu, k // R]
        off = idx[ju, k % R]
        peak = mp.reshape(len(iu), -1)[np.arange(len(iu)), k]
    else:   # hub transfer: separate feeder flights to/from metro nodes on the trunk
        nodes = trunk.node_pts
        met = trunk.pt_node[nodes]
        X = cfg["transfer_min"]
        Dn = D[:, met]                               # metro -> hub metro
        tf = np.minimum(t_direct(Dn, "sub"), t_direct(Dn, "bl")) if cfg["feeder"] != "sub" else t_direct(Dn, "sub")
        tf = np.where(Dn < 1, 0.0, tf + X)           # no feeder (or transfer) if already at the hub
        Lnn = trunk.Dt[np.ix_(nodes, nodes)]
        ttr = np.full(Lnn.shape, np.inf)
        pk = np.zeros(Lnn.shape)
        for f in feeders:
            t, m_ = t_trunk(0.0, Lnn, 0.0, mach, f, 0.0, cfg["a_acc"], cfg["a_dec"])
            upd = t < ttr
            ttr, pk = np.where(upd, t, ttr), np.where(upd, m_, pk)
        tot = tf[iu][:, :, None] + ttr[None, :, :] + tf[ju][:, None, :]
        flat = tot.reshape(len(iu), -1)
        k = np.argmin(flat, 1)
        t_tr = flat[np.arange(len(iu)), k]
        nn = len(nodes)
        on, off = nodes[k // nn], nodes[k % nn]
        peak = pk[k // nn, k % nn]
    sav_bl = np.where(t_bl - t_tr >= cfg["cut_min"], t_bl - t_tr, 0.0)
    use = sav_bl > 0
    sav_sub = np.where(use, t_sub - t_tr, 0.0)   # trunk users (Boomless world) vs subsonic
    # in a world with no Boomless Cruise the trunk is used whenever it beats subsonic
    sav_sub_world = np.where(t_sub - t_tr >= cfg["cut_min"], t_sub - t_tr, 0.0)
    out = dict(i=iu, j=ju, pax=pax, d=d, t_sub=t_sub, t_bl=t_bl, t_tr=t_tr, use=use,
               sav_bl=sav_bl, sav_sub=sav_sub, sav_sub_world=sav_sub_world,
               on=on, off=off, peak=peak)
    out["B_bl"] = float((pax * sav_bl).sum() / 60)            # pax-h/yr at 100 % capture
    out["B_sub"] = float((pax * sav_sub_world).sum() / 60)
    return out


def segment_loads(trunk, ev, key="sav_bl", use_key="use"):
    """Pax/yr and pax-h/yr (attributed pro rata by trunk km) on every trunk sub-segment."""
    load = {}
    use = ev[use_key]
    path_cache = {}
    for n_, (p, s, a, b) in enumerate(zip(ev["pax"], ev[key], ev["on"], ev["off"])):
        if not use[n_] or a == b:
            continue
        kk = (min(a, b), max(a, b))
        if kk not in path_cache:
            path_cache[kk] = nx.shortest_path(trunk.G, kk[0], kk[1], weight="length")
        path = path_cache[kk]
        Ltot = trunk.Dt[a, b]
        for u, v in zip(path[:-1], path[1:]):
            e = (min(u, v), max(u, v))
            ln = trunk.G[u][v]["length"]
            rec = load.setdefault(e, [0.0, 0.0])
            rec[0] += p
            rec[1] += p * s / 60 * ln / Ltot
    return load


# =====================================================================================
# 5. Tree construction
# =====================================================================================
HEURISTIC = {
    # user's eyeballed shape: West {SFO, LAX, PHX} - LAS - DEN - ... - CHI - East {YYZ, BOS, NYC, PHL, WAS}
    "heuristic_star": [("SFO", "LAS"), ("LAX", "LAS"), ("PHX", "LAS"), ("LAS", "DEN"), ("DEN", "CHI"),
                       ("CHI", "YYZ"), ("CHI", "BOS"), ("CHI", "NYC"), ("CHI", "PHL"), ("CHI", "WAS")],
    # same, but the eastern fan splits at Cleveland instead of five parallel spokes from Chicago
    "heuristic_cle": [("SFO", "LAS"), ("LAX", "LAS"), ("PHX", "LAS"), ("LAS", "DEN"), ("DEN", "CHI"),
                      ("CHI", "CLE"), ("CLE", "YYZ"), ("CLE", "BOS"), ("CLE", "NYC"), ("CLE", "PHL"),
                      ("CLE", "WAS")],
}

BASE_CFG = dict(connection="same", feeder="best", merge_min=5.0, transfer_min=60.0, cut_min=10.0,
                R=10, a_acc=A_ACC, a_dec=A_DEC)


def candidate_graph(df, D, max_km=1600.0, direct_km=1100.0):
    lat0 = np.radians(df.lat.mean())
    P = np.c_[df.lon.values * math.cos(lat0), df.lat.values]
    tri = Delaunay(P)
    G = nx.Graph()
    for s in tri.simplices:
        for a, b in ((s[0], s[1]), (s[1], s[2]), (s[0], s[2])):
            if D[a, b] <= max_km:
                G.add_edge(int(a), int(b), length=float(D[a, b]))
    n = len(D)
    for a in range(n):          # plus direct edges up to direct_km (straighter trunks)
        for b in range(a + 1, n):
            if D[a, b] <= direct_km:
                G.add_edge(a, b, length=float(D[a, b]))
    return G


def tree_from_nodes_edges(df, idx_edges):
    return [(df.code[a], df.code[b]) for a, b in idx_edges]


def score(df, D, T, edges, mach, cfg, pairs):
    tr = Trunk(df, edges)
    ev = evaluate(tr, D, T, mach, cfg, pairs)
    return ev["B_bl"], tr.km


def greedy_tree(df, D, T, mach, cfg, theta=0.5, n_seed=25, log=print):
    n = len(df)
    G = candidate_graph(df, D)
    iu, ju = np.triu_indices(n, 1)
    # cut 1: pairs that gain < cut_min even on a dedicated great-circle shield
    ub, _ = t_trunk(0.0, D[iu, ju], 0.0, mach, "sub", cfg["merge_min"], cfg["a_acc"], cfg["a_dec"])
    ub_bl, _ = t_trunk(0.0, D[iu, ju], 0.0, mach, "bl", cfg["merge_min"], cfg["a_acc"], cfg["a_dec"])
    gain = t_direct(D[iu, ju], "bl") - np.minimum(ub, ub_bl)
    keep = gain >= cfg["cut_min"]
    pairs = (iu[keep], ju[keep])
    info = dict(pairs_total=len(iu), pairs_kept_ub=int(keep.sum()),
                pax_total=float(T[iu, ju].sum()), pax_kept_ub=float(T[iu, ju][keep].sum()))
    # seeds: shortest candidate-graph paths of the biggest (pax x gain) pairs
    w = T[iu, ju] * np.where(keep, gain, 0)
    order = np.argsort(-w)[:n_seed]
    best = None
    for k in order:
        path = nx.shortest_path(G, int(iu[k]), int(ju[k]), weight="length")
        e = list(zip(path[:-1], path[1:]))
        B, km = score(df, D, T, tree_from_nodes_edges(df, e), mach, cfg, pairs)
        if best is None or B > best[0]:
            best = (B, km, e)
    B, km, edges = best
    tree_nodes = set(a for e in edges for a in e)
    log(f"  M{mach}: seed {[df.code[a] for a in tree_nodes]} B={B / 1e6:.2f}M pax-h, {km:.0f} km")
    growth = [dict(step=0, action="seed", added="-".join(df.code[a] for a in _path_order(edges)),
                   km=km, B_bl=B, marginal_per_km=B / km, avg_per_km=B / km)]
    accepted = True
    step = 0
    while True:
        step += 1
        dist, paths = nx.multi_source_dijkstra(G, tree_nodes, weight="length")
        cands = {}
        for v, p in paths.items():
            if v in tree_nodes or len(p) < 2:
                continue
            cands[v] = list(zip(p[:-1], p[1:]))
        if not cands:
            break
        res = []
        for v, ce in cands.items():
            Bn, kmn = score(df, D, T, tree_from_nodes_edges(df, edges + ce), mach, cfg, pairs)
            res.append(((Bn - B) / max(kmn - km, 1), Bn, kmn, ce, v))
        res.sort(key=lambda r: -r[0])
        r, Bn, kmn, ce, v = res[0]
        if r <= 0:
            break
        ok = r >= theta * B / km
        growth.append(dict(step=step, action="add" if (ok and accepted) else "beyond_stop",
                           added="-".join(df.code[a] for a in [ce[0][0]] + [b for _, b in ce]),
                           km=kmn, B_bl=Bn, marginal_per_km=r, avg_per_km=Bn / kmn))
        if not ok and accepted:
            accepted = False
            final = (list(edges), B, km)
        if not accepted and (r < 0.1 * B / km or step > 60):
            break
        edges = edges + ce
        tree_nodes |= {a for e in ce for a in e}
        B, km = Bn, kmn
    if accepted:
        final = (list(edges), B, km)
    edges, B, km = final
    # prune: drop leaf edges whose loss per km < theta x average
    while True:
        deg = {}
        for a, b in edges:
            deg[a] = deg.get(a, 0) + 1
            deg[b] = deg.get(b, 0) + 1
        best_rm = None
        for e in edges:
            if deg[e[0]] == 1 or deg[e[1]] == 1:
                rest = [x for x in edges if x != e]
                if not rest:
                    continue
                Bn, kmn = score(df, D, T, tree_from_nodes_edges(df, rest), mach, cfg, pairs)
                loss = (B - Bn) / max(km - kmn, 1)
                if best_rm is None or loss < best_rm[0]:
                    best_rm = (loss, e, Bn, kmn)
        if best_rm and best_rm[0] < theta * B / km:
            loss, e, Bn, kmn = best_rm
            edges = [x for x in edges if x != e]
            growth.append(dict(step=len(growth), action="prune", added=f"{df.code[e[0]]}-{df.code[e[1]]}",
                               km=kmn, B_bl=Bn, marginal_per_km=loss, avg_per_km=Bn / kmn))
            B, km = Bn, kmn
        else:
            break
    log(f"  M{mach}: final {km:.0f} km, B={B / 1e6:.2f}M pax-h/yr vs Boomless (100 % capture)")
    return tree_from_nodes_edges(df, edges), growth, info


def _path_order(edges):
    G = nx.Graph(edges)
    ends = [n for n in G if G.degree(n) == 1]
    return nx.shortest_path(G, ends[0], ends[1]) if len(ends) == 2 else list(G.nodes)


# =====================================================================================
# 6. Ramps: acceleration and Mach-onset (focus) zones
# =====================================================================================
ONSET = (1.10, 1.30)   # Mach band where accelerating aircraft focus a superboom (focus_booms.md: 1.16-1.28)


def ramp_zones(trunk, ev, mach, min_share=0.0):
    """Trunk sub-segments inside acceleration zones and Mach-onset zones (pax-weighted).
    Conservative: subsonic feeder (accel from M0.95); if the pair joins at its own
    airport the accel starts after the 230 km climb."""
    acc = {}
    ons = {}
    ramps = {}
    on0, on1 = accel_km(0.95, ONSET[0]), accel_km(0.95, ONSET[1])
    Dm = trunk.Dm
    for n_, (p, a, b, i, pk) in enumerate(zip(ev["pax"], ev["on"], ev["off"], ev["i"], ev["peak"])):
        if not ev["use"][n_] or a == b:
            continue
        for (src, dst, org) in ((a, b, ev["i"][n_]), (b, a, ev["j"][n_])):
            # both directions of the O&D flow (half the pax each way)
            off0 = max(0.0, CLIMB[1] - Dm[org, src] * PATH)
            path = nx.shortest_path(trunk.G, src, dst, weight="length")
            s = 0.0
            r = ramps.setdefault(src, [0.0, 0.0])
            r[0] += p / 2
            for u, v in zip(path[:-1], path[1:]):
                ln = trunk.G[u][v]["length"]
                e = (min(u, v), max(u, v))
                lo, hi = s, s + ln
                d_acc_p = accel_km(0.95, pk) if pk > 0 else 0
                if hi > off0 and lo < off0 + d_acc_p:
                    acc[e] = acc.get(e, 0) + p / 2
                if hi > off0 + on0 and lo < off0 + on1:
                    ons[e] = ons.get(e, 0) + p / 2
                s = hi
                if s > off0 + max(d_acc_p, on1):
                    break
    return acc, ons, ramps


# =====================================================================================
# 7. Economics (../../economics/cost_model.py)
# =====================================================================================
sys.path.insert(0, str(BS / "economics"))
import cost_model as cm  # noqa: E402

cm.load_operations()
MACH_DESIGN = {
    # Mach: (H2 thickness m, widths low/base/high km, note)  -- STATUS.md §4/§6
    1.4: (36.0, (1.5, 1.5, 2.0), "lift-set floor (psf alone needs ~28 m, Concorde M1.4)"),
    1.7: (36.0, (1.5, 1.5, 2.0), "economics.md M1.7 design (Overture, lift-set)"),
    2.0: (55.0, (1.5, 1.75, 2.0), "economics.md M2.0 psf baseline"),
    2.5: (77.0, (1.5, 1.75, 2.0), "psf_sweep: 77 m at 20 km (Concorde M2.5)"),
    3.2: (75.0, (1.5, 1.75, 2.0), "psf_sweep: 71 m (SR-71 M3.2; q > 35 kPa at 20 km, out of envelope)"),
}
for M_, (T_, W_, _) in MACH_DESIGN.items():
    cm.DESIGNS[f"trunk M{M_}"] = {"T": T_, "W": W_}


def annual_cost_per_km(mach, case="base"):
    d = f"trunk M{mach}"
    cx = cm.capex(d, case)
    ox = cm.opex(cx, case)
    return cm.annual_cost_km(cx, ox, case), cx["musd_km"], cx["usd_m2"]


SEATS = cm.SEATS
LOAD_FACTOR = cm.LOAD_FACTOR
LID_CAP_FLIGHTS_DAY = cm.FLIGHTS_CAP          # ~1e6 lid cycles / 30 yr ~ 90/day; 100 used
CAP_PAX_YR = LID_CAP_FLIGHTS_DAY * SEATS * LOAD_FACTOR * 365   # through any one km
PAIRWISE_BCR = 1.5e-3       # economics.md: DEN-ORD vs Boomless, M1.7 base, 10 % capture
FOCUS_KM = None             # set from ramp analysis


def _albers(lat, lon, lat1=29.5, lat2=45.5, lat0=37.5, lon0=-96.0):
    p1, p2, p0 = map(math.radians, (lat1, lat2, lat0))
    n = (math.sin(p1) + math.sin(p2)) / 2
    C = math.cos(p1) ** 2 + 2 * n * math.sin(p1)
    r0 = math.sqrt(C - 2 * n * math.sin(p0)) / n
    lat, lon = np.radians(np.asarray(lat, float)), np.radians(np.asarray(lon, float) - lon0)
    r = np.sqrt(C - 2 * n * np.sin(lat)) / n
    th = n * lon
    return R_E * r * np.sin(th), R_E * (r0 - r * np.cos(th))


def _basemap(ax, extent):
    land = BS / "network-maps" / "routes" / ".cache" / "ne_50m_land.geojson"
    layers = [(land, "fill", "#efeee9", None), (CACHE / "ne_50m_lakes.geojson", "fill", "#fcfcfb", None),
              (CACHE / "ne_50m_admin_1_states_provinces_lines.geojson", "line", "#d9d7d0", 0.4),
              (CACHE / "ne_50m_admin_0_boundary_lines_land.geojson", "line", "#b9b6ad", 0.8)]
    for f, kind, col, lw in layers:
        if not f.exists():
            continue
        gj = json.loads(f.read_text())
        for ft in gj["features"]:
            g = ft["geometry"]
            if g is None:
                continue
            if g["type"] in ("Polygon", "MultiPolygon"):
                polys = g["coordinates"] if g["type"] == "MultiPolygon" else [g["coordinates"]]
                for poly in polys:
                    ring = np.array(poly[0])
                    if ring[:, 0].max() < -170 or ring[:, 0].min() > -45 or ring[:, 1].max() < 10 or ring[:, 1].min() > 65:
                        continue
                    ring = ring[(ring[:, 1] > 5)]
                    if len(ring) < 3:
                        continue
                    x, y = _albers(ring[:, 1], np.clip(ring[:, 0], -170, -40))
                    ax.fill(x, y, color=col, lw=0, zorder=0)
            else:
                lines = g["coordinates"] if g["type"] == "MultiLineString" else [g["coordinates"]]
                for ln in lines:
                    a = np.array(ln)
                    if a[:, 0].max() < -170 or a[:, 0].min() > -45 or a[:, 1].min() > 65 or a[:, 1].max() < 10:
                        continue
                    x, y = _albers(a[:, 1], a[:, 0])
                    ax.plot(x, y, color=col, lw=lw, zorder=1)
    ax.set_xlim(extent[0], extent[1])
    ax.set_ylim(extent[2], extent[3])
    ax.set_aspect("equal")
    ax.axis("off")


def _draw_trunk(ax, data, wscale, color, onset_color=None, min_w=0.6):
    G = data["G"]
    ll = data["ll"]
    for (u, v), (p, ph) in data["lb"].items():
        x, y = _albers(ll[[u, v], 0], ll[[u, v], 1])
        ax.plot(x, y, color=color, lw=max(min_w, p / 1e6 * wscale), solid_capstyle="round", zorder=4)
    for u in G:
        for v in G[u]:
            if (min(u, v), max(u, v)) not in data["lb"]:
                x, y = _albers(ll[[u, v], 0], ll[[u, v], 1])
                ax.plot(x, y, color=color, lw=min_w, alpha=0.5, zorder=3)
    if onset_color:
        flow = {e: v[0] for e, v in data["lb"].items()}
        for e, p in data["ons"].items():
            if p >= 0.05 * flow.get(e, 1e18):
                x, y = _albers(ll[list(e), 0], ll[list(e), 1])
                ax.plot(x, y, color=onset_color, lw=2.2, zorder=5, solid_capstyle="butt")


def make_figures(df):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D
    blue, orange, ink, ink2, grid, bg = "#2a78d6", "#eb6834", "#0b0b0b", "#52514e", "#e4e3df", "#fcfcfb"
    plt.rcParams.update({"font.size": 9, "axes.edgecolor": grid, "axes.labelcolor": ink2,
                         "xtick.color": ink2, "ytick.color": ink2})
    summ = pd.read_csv(HERE / "trunk_summary.csv")
    head_M = 2.0
    x0, y0 = _albers(np.array([24, 24, 50, 50]), np.array([-125, -66, -125, -66]))
    ext = (x0.min() - 50, x0.max() + 120, y0.min() - 150, y0.max() + 60)
    mx, my = _albers(df.lat.values, df.lon.values)

    # ---- Figure 1: headline map --------------------------------------------------------
    fig, ax = plt.subplots(figsize=(12, 6.5))
    fig.patch.set_facecolor(bg)
    _basemap(ax, ext)
    data = pickle.loads((CACHE / f"ev_greedy_M{head_M}_M{head_M}.pkl").read_bytes())
    peak = max(p for p, _ in data["lb"].values()) / 1e6
    wscale = 9.0 / peak
    hc = pickle.loads((CACHE / f"ev_heuristic_cle_M{head_M}.pkl").read_bytes())
    for (u, v) in hc["lb"]:
        x, y = _albers(hc["ll"][[u, v], 0], hc["ll"][[u, v], 1])
        ax.plot(x, y, color=ink2, lw=0.9, ls=(0, (3, 2)), zorder=2)
    _draw_trunk(ax, data, wscale, blue, onset_color=orange)
    ax.scatter(mx, my, s=df.mass_M.values * 1.6, color=ink, zorder=6, lw=0.6, edgecolor=bg)
    rp = pd.read_csv(HERE / "ramps.csv")
    rp = rp[(rp.tree == f"greedy_M{head_M}") & (rp.mach == head_M)]
    rx, ry = _albers(rp.lat.values, rp.lon.values)
    ax.scatter(rx, ry, s=rp.joining_pax_M.values * 14, facecolor="none", edgecolor=orange, lw=1.4, zorder=7)
    lab = set(df.code[df.mass_M >= 6.5]) | {"SLC", "OMA", "CLE", "PIT", "IND", "BOI", "YYZ", "YUL", "YVR", "MEX"}
    for c, x, y in zip(df.code, mx, my):
        if c in lab and ext[0] < x < ext[1] and ext[2] < y < ext[3]:
            ax.text(x + 25, y + 25, c, fontsize=7.5, color=ink2, zorder=7)
    hs = [Line2D([], [], color=blue, lw=max(0.6, v * wscale), label=f"{v:g} M pax/yr")
          for v in (20, 40, 80) if v < peak * 1.3]
    hs += [plt.scatter([], [], s=10 * 14, facecolor="none", edgecolor=orange, lw=1.4,
                       label="on-ramp, area = pax/yr joining (10 M shown)")]
    hs += [Line2D([], [], color=orange, lw=2.2, label="Mach-onset (focus) zone, M1.1-1.3"),
           Line2D([], [], color=ink2, lw=0.9, ls=(0, (3, 2)), label="user heuristic tree (Cleveland fan)")]
    ax.legend(handles=hs, loc="lower left", frameon=False, fontsize=8)
    row = summ[(summ.tree == f"greedy_M{head_M}") & (summ.mach == head_M) &
               (summ.variant == "same_aircraft") & (summ.demand == "base")].iloc[0]
    ax.set_title(f"Greedy BoomShield trunk at Mach {head_M:g}: {row.trunk_km:,.0f} km, "
                 f"{row.pax_using_M:.0f} M pax/yr routed (100 % capture), busiest km {row.busiest_seg_pax_M:.1f} M pax/yr\n"
                 "line width = O&D passengers per year through each km; gravity-model demand (unverified)",
                 loc="left", color=ink, fontsize=10)
    fig.tight_layout()
    fig.savefig(FIG / "trunk_map.png", dpi=150, bbox_inches="tight", pad_inches=0.15, facecolor=bg)
    plt.close(fig)

    # ---- Figure 2: trees by Mach (small multiples, same width scale) -------------------
    fig, axs = plt.subplots(2, 3, figsize=(13, 6.4))
    fig.patch.set_facecolor(bg)
    panels = [("greedy_M1.4", 1.4, "greedy spine, Mach 1.4"), ("greedy_M2.0", 2.0, "greedy spine, Mach 2.0"),
              ("greedy_M3.2", 3.2, "greedy, Mach 3.2"), ("greedy_branched_M2.0", 2.0, "greedy branched, Mach 2.0"),
              ("greedy_branched_M3.2", 3.2, "greedy branched, Mach 3.2"), ("heuristic_cle", 2.0, "user heuristic, Mach 2.0")]
    for ax, (name, mm, t) in zip(axs.flat, panels):
        _basemap(ax, ext)
        fn = CACHE / f"ev_{name}_M{mm}.pkl"
        d = pickle.loads(fn.read_bytes())
        _draw_trunk(ax, d, wscale * 0.6, blue, onset_color=orange)
        ax.scatter(mx, my, s=df.mass_M.values * 0.6, color=ink, zorder=6, lw=0)
        r = summ[(summ.tree == name) & (summ.mach == mm) & (summ.variant == "same_aircraft") &
                 (summ.demand == "base")].iloc[0]
        ax.set_title(f"{t}: {r.trunk_km:,.0f} km, peak {r.busiest_seg_pax_M:.0f} M pax/yr,\n"
                     f"B/C vs Boomless {r.bcr_vs_bl_cap100:.3f} (100 % capture)", loc="left", fontsize=8.5, color=ink)
    fig.tight_layout()
    fig.savefig(FIG / "trunk_trees_by_mach.png", dpi=130, bbox_inches="tight", pad_inches=0.15, facecolor=bg)
    plt.close(fig)

    # ---- Figure 3: benefit/cost vs Mach --------------------------------------------------
    g = summ[(summ.variant == "same_aircraft") & (summ.demand == "base") &
             (summ.tree == "greedy_M" + summ.mach.astype(str))].sort_values("mach")
    fig, (a1, a2) = plt.subplots(1, 2, figsize=(13.5, 4.8))
    fig.patch.set_facecolor(bg)
    a1.plot(g.mach, g.densest_km_pax_h_vs_bl / 1e3, color=blue, lw=2, marker="o", ms=5, label="densest km, vs Boomless")
    a1.plot(g.mach, g.pax_h_per_km_vs_bl / 1e3, color=blue, lw=2, ls="--", marker="o", ms=5, label="trunk average, vs Boomless")
    a1.plot(g.mach, g.densest_km_pax_h_vs_sub / 1e3, color=orange, lw=2, marker="s", ms=5, label="densest km, vs subsonic")
    a1.plot(g.mach, g.pax_h_per_km_vs_sub / 1e3, color=orange, lw=2, ls="--", marker="s", ms=5, label="trunk average, vs subsonic")
    a1.set_xlabel("trunk cruise Mach")
    a1.set_ylabel("thousand pax-hours saved per km per year (100 % capture)")
    a1.set_title("Benefit per trunk km", loc="left", color=ink, fontsize=10)
    a1.legend(frameon=False, fontsize=8, loc="upper left")
    a1.grid(color=grid, lw=0.6)
    a1.set_ylim(0, 95)
    for c, ls in zip(CAPTURES, (":", "--", "-")):
        a2.plot(g.mach, g[f"bcr_vs_bl_cap{int(c * 100)}"], color=blue, lw=2, ls=ls, marker="o", ms=4,
                label=f"vs Boomless, {int(c * 100)} % capture")
    a2.plot(g.mach, g["bcr_vs_bl_cap100_lidcapped"], color=ink2, lw=1.6, marker="v", ms=4,
            label=f"vs Boomless, 100 %, lid cap {LID_CAP_FLIGHTS_DAY}/day")
    a2.plot(g.mach, g["bcr_vs_sub_cap100"], color=orange, lw=2, marker="s", ms=4, label="vs subsonic, 100 % capture")
    a2.axhline(PAIRWISE_BCR, color=ink2, lw=0.8, ls=(0, (2, 2)))
    a2.text(1.42, PAIRWISE_BCR * 1.15, "pairwise DEN-ORD (economics.md)", fontsize=7.5, color=ink2)
    a2.axhline(1.0, color=ink, lw=0.8)
    a2.text(1.42, 1.15, "B/C = 1", fontsize=7.5, color=ink)
    a2.set_yscale("log")
    a2.set_ylim(5e-4, 2)
    a2.set_xlabel("trunk cruise Mach")
    a2.set_ylabel("benefit / annualised cost, whole trunk (base cost)")
    a2.set_title("Benefit/cost, greedy trunk per Mach", loc="left", color=ink, fontsize=10)
    a2.grid(color=grid, lw=0.6, which="both")
    a2.legend(frameon=False, fontsize=7.5, loc="upper left", bbox_to_anchor=(1.01, 1.0))
    for a in (a1, a2):
        for sp in ("top", "right"):
            a.spines[sp].set_visible(False)
    fig.tight_layout()
    fig.savefig(FIG / "trunk_bcr_vs_mach.png", dpi=150, bbox_inches="tight", pad_inches=0.15, facecolor=bg)
    plt.close(fig)


# =====================================================================================
# 8. Driver
# =====================================================================================
MACHS = [1.4, 1.7, 2.0, 2.5, 3.2]
THETA = 0.5            # stop growing when marginal pax-h/km < THETA x trunk average
THETA_BRANCHED = 0.25
CAPTURES = [0.10, 0.30, 1.00]
VARIANTS = {
    "same_aircraft": dict(BASE_CFG),
    "same_aircraft_subsonic_feeder": dict(BASE_CFG, feeder="sub"),
    "hub_transfer_60min": dict(BASE_CFG, connection="hub", transfer_min=60.0),
    "hub_transfer_30min": dict(BASE_CFG, connection="hub", transfer_min=30.0),
    "same_aircraft_accel_x2": dict(BASE_CFG, a_acc=2 * A_ACC, a_dec=2 * A_DEC),
}


def edge_name(trunk, e):
    """Metro-to-metro tree edge that sub-segment e belongs to."""
    ei = trunk.G[e[0]][e[1]]["edge"]
    return ei, f"{trunk.edges[ei][0]}-{trunk.edges[ei][1]}"


def summarise(tree_name, trunk, ev, mach, variant, dem, cost):
    use = ev["use"]
    lb = segment_loads(trunk, ev, "sav_bl", "use")
    ev_s = dict(ev, use=ev["sav_sub_world"] > 0)
    ls = segment_loads(trunk, ev_s, "sav_sub_world", "use")
    seg_pax = np.array([v[0] for v in lb.values()]) if lb else np.zeros(1)
    km = trunk.km
    ann = cost["ann_base"]
    vt = cm.VTTS["base"] * cm.v("availability", "base")
    row = dict(tree=tree_name, mach=mach, variant=variant, demand=dem, trunk_km=round(km),
               pairs_using=int(use.sum()), pax_using_M=ev["pax"][use].sum() / 1e6,
               share_US_OD_trips=ev["pax"][use].sum() / US_OD_TRIPS,
               mean_saving_vs_bl_min=float((ev["pax"] * ev["sav_bl"]).sum() / max(ev["pax"][use].sum(), 1)),
               pax_h_vs_bl_M=ev["B_bl"] / 1e6, pax_h_vs_sub_M=ev["B_sub"] / 1e6,
               pax_using_vs_sub_world_M=ev["pax"][ev["sav_sub_world"] > 0].sum() / 1e6,
               busiest_seg_pax_M=seg_pax.max() / 1e6,
               pax_h_per_km_vs_bl=ev["B_bl"] / km, pax_h_per_km_vs_sub=ev["B_sub"] / km,
               densest_km_pax_h_vs_bl=max([v[1] / trunk.G[e[0]][e[1]]["length"] for e, v in lb.items()] or [0]),
               densest_km_pax_h_vs_sub=max([v[1] / trunk.G[e[0]][e[1]]["length"] for e, v in ls.items()] or [0]),
               cost_musd_km_yr_base=ann / 1e6)
    for c in CAPTURES:
        # whole-trunk BCR; capacity-capped version limits captured pax on every km to the lid cap
        ben = sum(v[1] * c for v in lb.values()) * vt
        ben_cap = sum(v[1] * c * min(1.0, CAP_PAX_YR / max(v[0] * c, 1)) for v in lb.values()) * vt
        ben_s = sum(v[1] * c for v in ls.values()) * vt
        row[f"bcr_vs_bl_cap{int(c * 100)}"] = ben / (ann * km)
        row[f"bcr_vs_sub_cap{int(c * 100)}"] = ben_s / (ann * km)
        row[f"bcr_vs_bl_cap{int(c * 100)}_lidcapped"] = ben_cap / (ann * km)
        row[f"bcr_densest_km_vs_bl_cap{int(c * 100)}"] = row["densest_km_pax_h_vs_bl"] * c * vt / ann
        row[f"toll_usd_per_trip_cap{int(c * 100)}"] = ann * km / max(ev["pax"][use].sum() * c, 1)
    for cs in ("low", "high"):
        row[f"bcr_vs_bl_cap100_cost_{cs}"] = sum(v[1] for v in lb.values()) * vt / (cost[f"ann_{cs}"] * km)
    row["flights_day_busiest_cap100"] = seg_pax.max() / (SEATS * LOAD_FACTOR * 365)
    return row, lb, ls


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--quick", action="store_true", help="skip demand sensitivities")
    args = ap.parse_args()
    CACHE.mkdir(exist_ok=True)
    FIG.mkdir(exist_ok=True)
    df = metros_df()
    D = dist_matrix(df)
    cal, an = calibrate(df, D)
    T = demand(df, D, cal)
    print(f"gravity: beta={cal['beta']:.2f} d0={cal['d0']:.0f} km rms(log)={cal['rms_log']:.2f}; "
          f"modelled O&D {T.sum() / 2e6:.0f} M pax/yr")
    an.to_csv(HERE / "calibration_anchors.csv", index=False)
    m = df.mass_M.values * 1e6
    out = df.copy()
    out["modelled_od_to_other_metros_M"] = (T.sum(1) / 2 / 1e6).round(2)
    out["share_of_mass_to_modelled_metros"] = (T.sum(1) / (2 * m)).round(2)
    out.to_csv(HERE / "metros.csv", index=False)
    demands = {"base": T}
    if not args.quick:
        demands["beta0.25"] = demand(df, D, cal, beta=0.25)
        demands["low_x0.75"] = T * LOW_DEMAND

    # ---- trees ----------------------------------------------------------------------
    trees = dict(HEURISTIC)
    growth_all = []
    for M in MACHS:
        f = CACHE / f"greedy_M{M}.pkl"
        if f.exists():
            edges, growth, info = pickle.loads(f.read_bytes())
        else:
            edges, growth, info = greedy_tree(df, D, T, M, BASE_CFG)
            f.write_bytes(pickle.dumps((edges, growth, info)))
        trees[f"greedy_M{M}"] = edges
        for g in growth:
            growth_all.append(dict(mach=M, tree="greedy", **g))
        f = CACHE / f"greedy_branched_M{M}.pkl"     # looser stop rule: more branches
        if f.exists():
            e2, g2, _ = pickle.loads(f.read_bytes())
        else:
            e2, g2, i2 = greedy_tree(df, D, T, M, BASE_CFG, theta=THETA_BRANCHED)
            f.write_bytes(pickle.dumps((e2, g2, i2)))
        trees[f"greedy_branched_M{M}"] = e2
        for g in g2:
            growth_all.append(dict(mach=M, tree="greedy_branched", **g))
        print(f"M{M}: pairs kept after upper-bound cut {info['pairs_kept_ub']}/{info['pairs_total']}, "
              f"{info['pax_kept_ub'] / 1e6:.0f}/{info['pax_total'] / 1e6:.0f} M pax")
    pd.DataFrame(growth_all).round(1).to_csv(HERE / "greedy_growth.csv", index=False)
    pd.DataFrame([dict(tree=k, edges=" ".join(f"{a}-{b}" for a, b in v)) for k, v in trees.items()]) \
        .to_csv(HERE / "trees.csv", index=False)

    # ---- costs ----------------------------------------------------------------------
    costs = {}
    crow = []
    for M in MACHS:
        c = {}
        for case in cm.CASES:
            ann, musd, usd_m2 = annual_cost_per_km(M, case)
            c[f"ann_{case}"] = ann
            crow.append(dict(mach=M, case=case, h2_m=MACH_DESIGN[M][0],
                             width_km=MACH_DESIGN[M][1][cm.CASES.index(case)], usd_m2=round(usd_m2, 1),
                             capex_musd_km=round(musd, 1), annual_musd_km_yr=round(ann / 1e6, 2),
                             note=MACH_DESIGN[M][2]))
        c["ann_base"] = c["ann_base"]
        costs[M] = c
    pd.DataFrame(crow).to_csv(HERE / "cost_per_km_by_mach.csv", index=False)

    # ---- evaluations ----------------------------------------------------------------
    summ, segs, pairs_out, ramp_rows = [], [], [], []
    iu, ju = np.triu_indices(len(df), 1)
    od = pd.DataFrame(dict(a=df.code.values[iu], b=df.code.values[ju], distance_km=D[iu, ju].round(0),
                           pax_yr=T[iu, ju].round(0), t_sub_min=t_direct(D[iu, ju], "sub").round(1),
                           t_boomless_min=t_direct(D[iu, ju], "bl").round(1)))
    for tname, edges in trees.items():
        trunk = Trunk(df, edges)
        for M in MACHS:
            for vname, cfg in VARIANTS.items():
                for dname, Tx in demands.items():
                    if dname != "base" and vname != "same_aircraft":
                        continue
                    ev = evaluate(trunk, D, Tx, M, cfg)
                    row, lb, ls = summarise(tname, trunk, ev, M, vname, dname, costs[M])
                    summ.append(row)
                    if vname != "same_aircraft" or dname != "base":
                        continue
                    # per tree edge: traffic and benefit per km
                    agg = {}
                    for e, (p, ph) in lb.items():
                        ei, nm = edge_name(trunk, e)
                        ln = trunk.G[e[0]][e[1]]["length"]
                        r = agg.setdefault(ei, dict(edge=nm, km=0.0, pax_max=0.0, pax_km=0.0, ph_bl=0.0, ph_sub=0.0))
                        r["pax_max"] = max(r["pax_max"], p)
                        r["pax_km"] += p * ln
                        r["ph_bl"] += ph
                    for e, (p, ph) in ls.items():
                        ei, nm = edge_name(trunk, e)
                        r = agg.setdefault(ei, dict(edge=nm, km=0.0, pax_max=0.0, pax_km=0.0, ph_bl=0.0, ph_sub=0.0))
                        r["ph_sub"] += ph
                    for ei, (a, b) in enumerate(trunk.edges):
                        L = float(gc(*df.loc[df.code == a, ["lat", "lon"]].values[0],
                                     *df.loc[df.code == b, ["lat", "lon"]].values[0]))
                        r = agg.get(ei, dict(edge=f"{a}-{b}", km=0.0, pax_max=0.0, pax_km=0.0, ph_bl=0.0, ph_sub=0.0))
                        vt = cm.VTTS["base"] * cm.v("availability", "base")
                        ann = costs[M]["ann_base"]
                        segs.append(dict(tree=tname, mach=M, edge=r["edge"], km=round(L),
                                         pax_yr_mean_M=r["pax_km"] / L / 1e6, pax_yr_max_M=r["pax_max"] / 1e6,
                                         pax_h_per_km_vs_bl=r["ph_bl"] / L, pax_h_per_km_vs_sub=r["ph_sub"] / L,
                                         **{f"bcr_vs_bl_cap{int(c * 100)}": r["ph_bl"] / L * c * vt / ann for c in CAPTURES},
                                         bcr_vs_sub_cap100=r["ph_sub"] / L * vt / ann,
                                         flights_day_cap100=r["pax_km"] / L / (SEATS * LOAD_FACTOR * 365)))
                    if tname.startswith("greedy") and float(tname.split("M")[1]) == M or tname == "heuristic_cle":
                        od[f"t_trunk_{tname}_M{M}"] = np.where(np.isfinite(ev["t_tr"]), ev["t_tr"].round(1), np.nan)
                        acc, ons, ramps = ramp_zones(trunk, ev, M)
                        flow = {e: v[0] for e, v in lb.items()}
                        acc_km = sum(trunk.G[e[0]][e[1]]["length"] for e, p in acc.items() if p >= 0.05 * flow.get(e, 1e18))
                        ons_km = sum(trunk.G[e[0]][e[1]]["length"] for e, p in ons.items() if p >= 0.05 * flow.get(e, 1e18))
                        ons_any = sum(trunk.G[e[0]][e[1]]["length"] for e in ons)
                        tot_on = sum(v[0] for v in ramps.values())
                        main_r = sorted(ramps.items(), key=lambda kv: -kv[1][0])
                        for pid_, (p, _) in main_r:
                            if p < 0.02 * tot_on:
                                continue
                            la, lo = trunk.ll[pid_]
                            dd = gc(la, lo, df.lat.values, df.lon.values)
                            k = int(np.argmin(dd))
                            ramp_rows.append(dict(tree=tname, mach=M, lat=round(la, 3), lon=round(lo, 3),
                                                  nearest_metro=df.code[k], km_from_metro=round(float(dd[k])),
                                                  joining_pax_M=p / 1e6, share_of_joins=p / tot_on))
                        summ[-1].update(accel_zone_km=acc_km, onset_zone_km=ons_km, onset_zone_km_any=ons_any,
                                        n_ramps_2pct=sum(1 for _, (p, _) in main_r if p >= 0.02 * tot_on))
                        with open(CACHE / f"ev_{tname}_M{M}.pkl", "wb") as fh:
                            pickle.dump(dict(edges=edges, lb=lb, ls=ls, acc=acc, ons=ons,
                                             ll=trunk.ll, G=nx.to_dict_of_dicts(trunk.G)), fh)
        print(f"  evaluated {tname}")
        pd.DataFrame(summ).to_csv(HERE / "trunk_summary.csv", index=False, float_format="%.4g")
        pd.DataFrame(segs).to_csv(HERE / "trunk_segments.csv", index=False, float_format="%.4g")
    pd.DataFrame(ramp_rows).to_csv(HERE / "ramps.csv", index=False, float_format="%.4g")
    od.to_csv(HERE / "od_matrix.csv", index=False)
    print("done; figures next")
    make_figures(df)


if __name__ == "__main__":
    main()
