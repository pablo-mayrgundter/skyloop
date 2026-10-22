#!/usr/bin/env python3
"""
BoomShield route ranking and block-time model for Overture-class aircraft.

Builds on ../america-north (seed airports, top-pairs template, trunk hypothesis).

For each city pair it computes
  * great-circle distance, land fraction and boom-restricted ("overland") fraction
    (Natural Earth 1:50m land, plus a coastal buffer for the boom carpet),
  * the user's ranking metric  overland_km x annual seats,
  * block times for a Mach 0.85 subsonic jet, Overture with Boomless Cruise
    (Mach 1.15 effective over restricted ground, Mach 1.7 over open water) and
    Overture at Mach 1.7 everywhere on a BoomShield corridor,
  * time saved per passenger and passenger-hours saved per year,
  * a coastal-offset check (would a seaward detour at Mach 1.7 beat the shield?),
  * passenger-hours saved per km of the LAX-DEN-ORD-JFK trunk hypothesis.

Usage:  python3 build_routes.py [--capture 0.10] [--load-factor 0.8]
        [--buffer-km 40] [--fill-template]

Needs numpy, pandas, shapely, scipy, matplotlib. The Natural Earth land file is
downloaded once into .cache/ (public domain, naturalearthdata.com).
Light CPU: < 1 min.
"""

from __future__ import annotations

import argparse
import json
import math
import urllib.request
from pathlib import Path

import numpy as np
import pandas as pd
import shapely
from scipy.spatial import cKDTree

HERE = Path(__file__).resolve().parent
NA = HERE.parent / "america-north"
CACHE = HERE / ".cache"
FIG = HERE / "figures"
NE_URL = ("https://raw.githubusercontent.com/nvkelso/natural-earth-vector/"
          "master/geojson/ne_50m_land.geojson")
R_E = 6371.0

# ---------------------------------------------------------------------------
# Airports (IATA -> lat, lon).  North American seeds come from the scaffold.
# ---------------------------------------------------------------------------
EXTRA_AIRPORTS = {
    "LGA": (40.7769, -73.8740), "EWR": (40.6895, -74.1745), "BOS": (42.3656, -71.0096),
    "PHX": (33.4342, -112.0116), "FLL": (26.0742, -80.1506), "HNL": (21.3187, -157.9225),
    "OGG": (20.8986, -156.4305), "CJU": (33.5113, 126.4930), "GMP": (37.5583, 126.7906),
    "CTS": (42.7752, 141.6923), "HND": (35.5494, 139.7798), "FUK": (33.5859, 130.4511),
    "HAN": (21.2212, 105.8072), "SGN": (10.8188, 106.6520), "MEL": (-37.6690, 144.8410),
    "SYD": (-33.9399, 151.1753), "JED": (21.6796, 39.1565), "RUH": (24.9576, 46.6988),
    "OKA": (26.1958, 127.6459), "BOM": (19.0896, 72.8656), "DEL": (28.5562, 77.1000),
    "PEK": (40.0799, 116.6031), "SHA": (31.1979, 121.3363), "CAN": (23.3924, 113.2988),
    "SZX": (22.6393, 113.8107), "HGH": (30.2295, 120.4344), "CTU": (30.5785, 103.9471),
    "HKG": (22.3080, 113.9185), "TPE": (25.0797, 121.2342), "CAI": (30.1219, 31.4056),
    "ICN": (37.4602, 126.4407), "NRT": (35.7720, 140.3929), "KUL": (2.7456, 101.7072),
    "SIN": (1.3644, 103.9915), "KIX": (34.4347, 135.2440), "DXB": (25.2532, 55.3657),
    "BKK": (13.6900, 100.7501), "CGK": (-6.1256, 106.6559), "LHR": (51.4700, -0.4543),
    "BLR": (13.1986, 77.7066), "HYD": (17.2403, 78.4294), "CCU": (22.6547, 88.4467),
}

# ---------------------------------------------------------------------------
# Traffic.  value = annual seats (both directions) unless metric says otherwise.
# Pages could not be fetched from this sandbox (Wikipedia, OAG, BTS, DOT all
# blocked); figures come from search-result summaries of the cited pages and
# are marked accordingly.  'est' rows are order-of-magnitude placeholders.
# ---------------------------------------------------------------------------
SOURCES = {
    "OAG24W": "OAG Busiest Routes 2024 (Dec 2024), world domestic top 10, via "
              "airlineratings.com/articles/oag-reveals-the-busiest-airline-routes-of-2024 "
              "and airwaysmag.com/new-post/busiest-flight-routes-of-2024 (search summary; unverified)",
    "OAG24CN": "OAG Busiest Routes 2024, China list, via airwaysmag.com and "
               "aviationa2z.com 2024-12-18 (search summary; unverified)",
    "OAG24I": "OAG press release 2024 international top 10 "
              "(oag.com/pressroom/hong-kong-taipei-worlds-busiest-international-airline-route-2024; "
              "gulf-times.com) (search summary; unverified)",
    "OAG25US": "OAG Busiest Routes 2025, US domestic top 10, via onemileatatime.com/guides/"
               "busiest-airline-routes and travelweekly.com JFK-LAX 2025 (search summary; unverified)",
    "WIKI_IN": "Wikipedia 'List of busiest passenger flight routes', India table, FY2023-24 "
               "passengers (page blocked here; values from search summary; unverified)",
    "EST": "Author placeholder estimate (daily frequencies x gauge, order of magnitude); "
           "replace with BTS DB1B/T-100 or OAG schedules",
}

TRAFFIC = [
    # world domestic, OAG 2024 seats
    ("CJU", "GMP", 14.2e6, "seats", 2024, "OAG24W", "world #1 domestic"),
    ("CTS", "HND", 11.93e6, "seats", 2024, "OAG24W", ""),
    ("FUK", "HND", 11.34e6, "seats", 2024, "OAG24W", ""),
    ("HAN", "SGN", 10.63e6, "seats", 2024, "OAG24W", ""),
    ("MEL", "SYD", 9.22e6, "seats", 2024, "OAG24W", ""),
    ("JED", "RUH", 8.70e6, "seats", 2024, "OAG24W", "9.8M seats in 2025 (OAG)"),
    ("HND", "OKA", 8.03e6, "seats", 2024, "OAG24W", ""),
    ("BOM", "DEL", 7.96e6, "seats", 2024, "OAG24W", ""),
    ("PEK", "SHA", 7.71e6, "seats", 2024, "OAG24W", "7.5M in 2025"),
    ("CAN", "SHA", 7.01e6, "seats", 2024, "OAG24W", ""),
    ("SHA", "SZX", 6.86e6, "seats", 2024, "OAG24CN", ""),
    ("PEK", "SZX", 5.60e6, "seats", 2024, "OAG24CN", ""),
    ("CAN", "HGH", 4.64e6, "seats", 2024, "OAG24CN", ""),
    ("CTU", "PEK", 4.60e6, "seats", 2024, "OAG24CN", "CTU coords used; much traffic now TFU/PKX"),
    ("CAN", "PEK", 4.27e6, "seats", 2024, "OAG24CN", ""),
    # world international, OAG 2024 seats
    ("HKG", "TPE", 6.8e6, "seats", 2024, "OAG24I", ""),
    ("CAI", "JED", 5.5e6, "seats", 2024, "OAG24I", ""),
    ("ICN", "NRT", 5.4e6, "seats", 2024, "OAG24I", ""),
    ("KUL", "SIN", 5.4e6, "seats", 2024, "OAG24I", ""),
    ("ICN", "KIX", 4.98e6, "seats", 2024, "OAG24I", ""),
    ("DXB", "RUH", 4.3e6, "seats", 2024, "OAG24I", ""),
    ("BKK", "HKG", 4.2e6, "seats", 2024, "OAG24I", ""),
    ("CGK", "SIN", 4.07e6, "seats", 2024, "OAG24I", ""),
    ("BKK", "SIN", 4.03e6, "seats", 2024, "OAG24I", ""),
    ("JFK", "LHR", 4.01e6, "seats", 2024, "OAG24I", "transatlantic: Overture's design case"),
    # India, passengers
    ("BLR", "DEL", 4.78e6, "passengers", 2024, "WIKI_IN", "FY2023-24"),
    ("BLR", "BOM", 4.43e6, "passengers", 2024, "WIKI_IN", "FY2023-24"),
    ("DEL", "HYD", 2.92e6, "passengers", 2024, "WIKI_IN", "FY2023-24"),
    ("CCU", "DEL", 2.82e6, "passengers", 2024, "WIKI_IN", "FY2023-24"),
    # US domestic, OAG 2025 seats
    ("JFK", "LAX", 3431159, "seats", 2025, "OAG25US", "US #1 2025; IATA: 2.2M pax 2024"),
    ("LGA", "ORD", 3334011, "seats", 2025, "OAG25US", "NYC-Chicago trunk"),
    ("LAX", "SFO", 3306223, "seats", 2025, "OAG25US", ""),
    ("LAS", "LAX", 3283510, "seats", 2025, "OAG25US", ""),
    ("ATL", "MCO", 3226190, "seats", 2025, "OAG25US", ""),
    ("DEN", "PHX", 3118880, "seats", 2025, "OAG25US", ""),
    ("HNL", "OGG", 3.0e6, "seats", 2025, "OAG25US", ""),
    ("LAX", "ORD", 2.8e6, "seats", 2025, "OAG25US", ""),
    ("ATL", "LGA", 2763776, "seats", 2025, "OAG25US", ""),
    ("ATL", "FLL", 2.6e6, "seats", 2025, "OAG25US", ""),
    # North American trunks without a published figure found: placeholders
    ("SFO", "JFK", 2.2e6, "seats_est", 2025, "EST", ""),
    ("LAX", "ATL", 2.0e6, "seats_est", 2025, "EST", ""),
    ("SFO", "ORD", 1.9e6, "seats_est", 2025, "EST", ""),
    ("LAX", "DFW", 2.4e6, "seats_est", 2025, "EST", ""),
    ("LAX", "DEN", 2.2e6, "seats_est", 2025, "EST", ""),
    ("DEN", "ORD", 2.1e6, "seats_est", 2025, "EST", ""),
    ("LAX", "SEA", 2.3e6, "seats_est", 2025, "EST", ""),
    ("SEA", "SFO", 1.9e6, "seats_est", 2025, "EST", ""),
    ("LAX", "EWR", 1.5e6, "seats_est", 2025, "EST", ""),
    ("SFO", "EWR", 1.8e6, "seats_est", 2025, "EST", ""),
    ("BOS", "LAX", 1.3e6, "seats_est", 2025, "EST", ""),
    ("BOS", "SFO", 1.4e6, "seats_est", 2025, "EST", ""),
    ("SEA", "JFK", 1.0e6, "seats_est", 2025, "EST", ""),
    ("LAX", "MIA", 1.1e6, "seats_est", 2025, "EST", ""),
    ("DEN", "JFK", 0.9e6, "seats_est", 2025, "EST", ""),
    ("DFW", "JFK", 0.9e6, "seats_est", 2025, "EST", ""),
    ("DFW", "LGA", 1.6e6, "seats_est", 2025, "EST", ""),
    ("IAH", "JFK", 0.4e6, "seats_est", 2025, "EST", "Houston-NYC mostly via EWR"),
    ("IAH", "EWR", 1.5e6, "seats_est", 2025, "EST", ""),
    ("MIA", "JFK", 1.6e6, "seats_est", 2025, "EST", ""),
    ("ORD", "JFK", 0.7e6, "seats_est", 2025, "EST", "NYC-Chicago mostly LGA-ORD"),
    ("DFW", "ORD", 2.0e6, "seats_est", 2025, "EST", ""),
    ("DEN", "LAS", 2.6e6, "seats_est", 2025, "EST", ""),
    ("LAX", "PHX", 1.9e6, "seats_est", 2025, "EST", ""),
    ("DEN", "DFW", 1.6e6, "seats_est", 2025, "EST", ""),
    ("SEA", "ORD", 1.1e6, "seats_est", 2025, "EST", ""),
    ("ATL", "ORD", 1.5e6, "seats_est", 2025, "EST", ""),
]

# ---------------------------------------------------------------------------
# Flight-profile assumptions (screening values; see routes.md).
# Speeds: stratosphere a = 295.1 m/s (ISA, 11-20 km).
# ---------------------------------------------------------------------------
A_STRAT = 295.1 * 3.6  # km/h
V085 = 0.85 * A_STRAT  # 903 km/h
V094 = 0.94 * A_STRAT  # 998
V115 = 1.15 * A_STRAT  # 1222
V17 = 1.70 * A_STRAT   # 1806
TAXI_MIN = 15 + 7      # taxi-out + taxi-in, all types
APPROACH = (10.0, 60.0)  # (min, km) terminal area / approach, all types
PATH_FACTOR = 1.03       # flown / great-circle distance (airways, SIDs/STARs)

# (minutes, km) for ramps; ramps exclude the approach constant
CLIMB_SUB = (18.0, 230.0)    # brake release -> top of climb FL350-390
DESC_SUB = (18.0, 200.0)     # top of descent -> terminal entry
ACC_BL = (7.0, 130.0)        # M0.9 -> M1.15 at cutoff altitude (~FL450-500)
DEC_BL = (4.0, 75.0)
ACC_17 = (16.0, 390.0)       # M0.95 -> M1.7 with climb to ~FL600, no afterburner
DEC_17 = (8.0, 185.0)
DESC_17 = (22.0, 260.0)      # FL600 -> terminal entry
# Mach 1.15 -> 1.7 transition (over water in Boomless mode, or onto the shield)
ACC_BL17 = (ACC_17[0] - ACC_BL[0], ACC_17[1] - ACC_BL[1])     # 9 min, 260 km
DEC_17BL = (DEC_17[0] - DEC_BL[0], DEC_17[1] - DEC_BL[1])     # 4 min, 110 km

PROFILES = {
    # name: (up ramp, down ramp, restricted speed, can use M1.7 over open water)
    "sub":  ((CLIMB_SUB[0], CLIMB_SUB[1]), (DESC_SUB[0], DESC_SUB[1]), V085, False),
    "ov094": ((CLIMB_SUB[0], CLIMB_SUB[1]), (DESC_SUB[0], DESC_SUB[1]), V094, True),
    "bl":   ((CLIMB_SUB[0] + ACC_BL[0], CLIMB_SUB[1] + ACC_BL[1]),
             (DEC_BL[0] + DESC_SUB[0], DEC_BL[1] + DESC_SUB[1]), V115, True),
    "m17":  ((CLIMB_SUB[0] + ACC_17[0], CLIMB_SUB[1] + ACC_17[1]),
             (DEC_17[0] + DESC_17[0], DEC_17[1] + DESC_17[1]), V17, True),
}


# ---------------------------------------------------------------------------
# Geometry helpers
# ---------------------------------------------------------------------------
def to_xyz(lat, lon):
    la, lo = np.radians(lat), np.radians(lon)
    return np.stack([np.cos(la) * np.cos(lo), np.cos(la) * np.sin(lo), np.sin(la)], -1)


def to_latlon(v):
    v = v / np.linalg.norm(v, axis=-1, keepdims=True)
    return np.degrees(np.arcsin(v[..., 2])), np.degrees(np.arctan2(v[..., 1], v[..., 0]))


def gc_km(p, q):
    return R_E * 2 * math.asin(min(1.0, np.linalg.norm(to_xyz(*p) - to_xyz(*q)) / 2))


def gc_points(p, q, step_km=5.0):
    """Great circle p->q sampled at ~step_km (midpoints of equal sub-steps)."""
    a, b = to_xyz(*p), to_xyz(*q)
    om = math.acos(float(np.clip(a @ b, -1, 1)))
    d = om * R_E
    n = max(2, int(math.ceil(d / step_km)))
    t = (np.arange(n) + 0.5) / n
    v = (np.sin((1 - t) * om)[:, None] * a + np.sin(t * om)[:, None] * b) / math.sin(om)
    lat, lon = to_latlon(v)
    return lat, lon, np.full(n, d / n)


def destination(lat, lon, bearing_deg, dist_km):
    la, lo, br, dr = map(math.radians, (lat, lon, bearing_deg, 0.0))
    dr = dist_km / R_E
    la2 = math.asin(math.sin(la) * math.cos(dr) + math.cos(la) * math.sin(dr) * math.cos(br))
    lo2 = lo + math.atan2(math.sin(br) * math.sin(dr) * math.cos(la),
                          math.cos(dr) - math.sin(la) * math.sin(la2))
    return math.degrees(la2), (math.degrees(lo2) + 540) % 360 - 180


def bearing(p, q):
    la1, lo1, la2, lo2 = map(math.radians, (*p, *q))
    y = math.sin(lo2 - lo1) * math.cos(la2)
    x = math.cos(la1) * math.sin(la2) - math.sin(la1) * math.cos(la2) * math.cos(lo2 - lo1)
    return math.degrees(math.atan2(y, x))


def interp_gc(p, q, f):
    a, b = to_xyz(*p), to_xyz(*q)
    om = math.acos(float(np.clip(a @ b, -1, 1)))
    v = (math.sin((1 - f) * om) * a + math.sin(f * om) * b) / math.sin(om)
    la, lo = to_latlon(v)
    return float(la), float(lo)


class LandMask:
    """Natural Earth 1:50m land: point-in-land test plus distance-to-coast."""

    def __init__(self, buffer_km: float):
        CACHE.mkdir(exist_ok=True)
        f = CACHE / "ne_50m_land.geojson"
        if not f.exists():
            print("downloading Natural Earth 50m land ...")
            urllib.request.urlretrieve(NE_URL, f)
        gj = json.loads(f.read_text())
        geoms = [shapely.geometry.shape(ft["geometry"]) for ft in gj["features"]]
        polys = []
        for g in geoms:
            polys.extend(getattr(g, "geoms", [g]))
        self.tree = shapely.STRtree(polys)
        rings = []
        for pg in polys:
            rings.append(pg.exterior)
            rings.extend(pg.interiors)
        coast = shapely.segmentize(shapely.MultiLineString(rings), 0.02)  # ~2 km
        xy = shapely.get_coordinates(coast)
        self.kd = cKDTree(to_xyz(xy[:, 1], xy[:, 0]))
        self.buffer_km = buffer_km

    def classify(self, lat, lon):
        pts = shapely.points(lon, lat)
        idx = self.tree.query(pts, predicate="intersects")[0]
        land = np.zeros(len(lat), bool)
        land[idx] = True
        dist, _ = self.kd.query(to_xyz(lat, lon))
        coast_km = 2 * R_E * np.arcsin(np.clip(dist / 2, 0, 1))
        restricted = land | (coast_km <= self.buffer_km)
        return land, restricted


# ---------------------------------------------------------------------------
# Time model
# ---------------------------------------------------------------------------
def runs(flags, ds):
    """Contiguous runs of flags -> list of (flag, length_km)."""
    out = []
    for f, d in zip(flags, ds):
        if out and out[-1][0] == f:
            out[-1][1] += d
        else:
            out.append([bool(f), float(d)])
    return out


def water_time_h(length_km, from_v):
    """Time over an open-water run entered at speed from_v: go to M1.7 if long enough."""
    if from_v >= V17:
        return length_km / V17
    if from_v == V115:
        acc, dec = ACC_BL17, DEC_17BL
    else:  # from subsonic cruise
        acc, dec = ACC_17, DEC_17
    need = acc[1] + dec[1]
    if length_km < need:
        return length_km / from_v
    return (acc[0] + dec[0]) / 60 + (length_km - need) / V17


def block_time_h(profile, s_km, ds_km, restricted, shield=False):
    """Block time (h) for one profile along a sampled path.

    s_km: cumulative distance at sample centres (already x PATH_FACTOR).
    restricted: True where a boom must not reach the ground (land + buffer).
    shield: if True, restricted ground is covered by BoomShield (only for m17).
    Returns inf if the profile cannot fit its ramps into the route.
    """
    up, dn, v_res, water17 = PROFILES[profile]
    L = s_km[-1] + ds_km[-1] / 2
    ramp_d = up[1] + dn[1] + APPROACH[1]
    ramp_t = (up[0] + dn[0] + APPROACH[0]) / 60
    if L < ramp_d:
        if profile in ("sub", "ov094"):
            # truncated profile: same mean speed as the full ramps
            return TAXI_MIN / 60 + ramp_t * L / ramp_d
        return math.inf
    win = (s_km >= up[1]) & (s_km <= L - dn[1] - APPROACH[1])
    t = TAXI_MIN / 60 + ramp_t
    if profile == "m17" and not shield:
        # above the cutoff Mach (end of Boomless accel .. start of Boomless decel)
        # the boom reaches the ground, so that stretch must be open water
        sup = ((s_km >= PROFILES["bl"][0][1])
               & (s_km <= L - PROFILES["bl"][1][1] - APPROACH[1]))
        if restricted[sup].any():
            return math.inf
    if profile == "m17":
        return t + ds_km[win].sum() / V17
    for res, length in runs(restricted[win], ds_km[win]):
        if res or not water17:
            t += length / v_res
        else:
            t += water_time_h(length, v_res)
    return t


def overture_times(s, ds, restricted):
    t_sub = block_time_h("sub", s, ds, restricted)
    t_094 = block_time_h("ov094", s, ds, restricted)
    t_bl = min(block_time_h("bl", s, ds, restricted), t_094)
    t_m17_noshield = min(block_time_h("m17", s, ds, restricted, shield=False), t_bl)
    t_m17 = min(block_time_h("m17", s, ds, restricted, shield=True), t_bl)
    return t_sub, t_094, t_bl, t_m17, t_m17_noshield


def path_samples(waypoints, step_km=5.0):
    lat, lon, ds = [], [], []
    for p, q in zip(waypoints[:-1], waypoints[1:]):
        a, b, d = gc_points(p, q, step_km)
        lat.append(a), lon.append(b), ds.append(d)
    lat, lon, ds = map(np.concatenate, (lat, lon, ds))
    ds = ds * PATH_FACTOR
    s = np.cumsum(ds) - ds / 2
    return lat, lon, s, ds


def coastal_offset(mask, p, q, t_bl_gc):
    """Best seaward detour (2 offset waypoints) flown boomless/M1.7, no shield."""
    best = None
    brg = bearing(p, q)
    for f in (0.2, 0.3):
        w1, w2 = interp_gc(p, q, f), interp_gc(p, q, 1 - f)
        for side in (90, -90):
            for off in (40, 80, 120, 180, 250, 350):
                a = destination(*w1, brg + side, off)
                b = destination(*w2, brg + side, off)
                lat, lon, s, ds = path_samples([p, a, b, q], 10.0)
                _, res = mask.classify(lat, lon)
                t = min(block_time_h("bl", s, ds, res), block_time_h("ov094", s, ds, res))
                if best is None or t < best[0]:
                    best = (t, s[-1] + ds[-1] / 2, res.mean(), off, f)
    if best is None or best[0] > t_bl_gc - 2 / 60:  # must beat boomless GC by 2 min
        return None
    return best


# ---------------------------------------------------------------------------
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--capture", type=float, default=0.10,
                    help="fraction of the route's passengers flying Overture")
    ap.add_argument("--load-factor", type=float, default=0.80)
    ap.add_argument("--buffer-km", type=float, default=40.0,
                    help="coastal buffer counted as boom-restricted")
    ap.add_argument("--trunk-buffer-km", type=float, default=200.0)
    ap.add_argument("--fill-template", action="store_true")
    args = ap.parse_args()

    seeds = pd.read_csv(NA / "north_america_seed_airports.csv")
    apt = {r.iata: (r.lat_deg, r.lon_deg) for r in seeds.itertuples()}
    apt.update(EXTRA_AIRPORTS)
    mask = LandMask(args.buffer_km)

    # trunk hypothesis -> densified KD-tree for corridor allocation
    trunks = pd.read_csv(NA / "north_america_shared_trunks_hypothesis.csv")
    seg_pts, seg_id = [], []
    for i, r in trunks.iterrows():
        la, lo, _ = gc_points(apt[r.from_iata], apt[r.to_iata], 5.0)
        seg_pts.append(to_xyz(la, lo))
        seg_id += [i] * len(la)
    seg_kd = cKDTree(np.concatenate(seg_pts))
    seg_id = np.array(seg_id)
    seg_alloc = {i: np.zeros(2) for i in trunks.index}  # [vs_sub, vs_boomless] pax-h

    rows = []
    for o, d, val, metric, year, src, note in TRAFFIC:
        p, q = apt[o], apt[d]
        dist = gc_km(p, q)
        lat, lon, s, ds = path_samples([p, q], 5.0)
        land, res = mask.classify(lat, lon)
        w = ds / ds.sum()
        land_frac, res_frac = float((w * land).sum()), float((w * res).sum())
        t_sub, t_094, t_bl, t_m17, t_m17_ns = overture_times(s, ds, res)
        seats = val / args.load_factor if metric == "passengers" else val
        pax = seats * args.load_factor * args.capture
        sav_sub = t_sub - t_m17
        sav_bl = t_bl - t_m17
        # shield length actually flown supersonic over restricted ground
        L = s[-1] + ds[-1] / 2
        up_bl = PROFILES["bl"][0][1]
        dn_bl = PROFILES["bl"][1][1] + APPROACH[1]
        win = (s >= up_bl) & (s <= L - dn_bl) & res
        shield_km = float(ds[win].sum() / PATH_FACTOR) if sav_bl > 1e-6 else 0.0
        notes = [note] if note else []
        if land_frac > 0.3 and res_frac - land_frac > 0.1:
            notes.append("coastal: GC within %d km of shore" % args.buffer_km)
        off = coastal_offset(mask, p, q, t_bl)
        off_txt = ""
        if off:
            t_off, L_off, rf_off, offkm, f = off
            off_txt = (f"seaward detour {offkm} km offset: +{L_off - L:.0f} km, "
                       f"{rf_off:.0%} restricted, block {t_off * 60:.0f} min "
                       f"(shield M1.7 {t_m17 * 60:.0f}, boomless GC {t_bl * 60:.0f})")
            notes.append(off_txt)
        if t_m17_ns < t_bl - 1e-6:
            notes.append("M1.7 possible without shield on GC (open water)")
        if sav_bl <= 1e-6 and L < sum(PROFILES["m17"][0][1:] + PROFILES["m17"][1][1:]) + APPROACH[1]:
            notes.append("too short to reach M1.7")
        # corridor allocation (North American routes only)
        if sav_bl > 1e-6 or sav_sub > 1e-6:
            sw = (s >= PROFILES["m17"][0][1]) & (s <= L - PROFILES["m17"][1][1] - APPROACH[1]) & res
            if sw.any():
                dd, ii = seg_kd.query(to_xyz(lat[sw], lon[sw]))
                near = 2 * R_E * np.arcsin(np.clip(dd / 2, 0, 1)) <= args.trunk_buffer_km
                tot = ds[sw].sum()
                for k in np.unique(seg_id[ii[near]]):
                    frac = ds[sw][near & (seg_id[ii] == k)].sum() / tot
                    seg_alloc[k] += frac * pax * np.array([sav_sub, sav_bl])
        rows.append(dict(
            origin=o, dest=d, distance_km=round(dist), land_fraction=round(land_frac, 3),
            overland_fraction=round(res_frac, 3), overland_km=round(res_frac * dist),
            annual_seats_or_pax=int(val), metric=metric, year=year, source=src,
            seats_equiv=int(seats),
            rank_metric_bn_seat_km=round(res_frac * dist * seats / 1e9, 2),
            time_subsonic_h=round(t_sub, 2), time_overture_M094_h=round(t_094, 2),
            time_boomless_h=round(t_bl, 2), time_M17_h=round(t_m17, 2),
            save_vs_subsonic_min=round(sav_sub * 60, 1),
            save_boomless_vs_subsonic_min=round((t_sub - t_bl) * 60, 1),
            save_vs_boomless_min=round(sav_bl * 60, 1),
            shield_km=round(shield_km),
            capture=args.capture, load_factor=args.load_factor,
            pax_hours_saved_vs_subsonic=round(pax * sav_sub),
            pax_hours_saved_vs_boomless=round(pax * sav_bl),
            notes="; ".join(notes)))

    df = pd.DataFrame(rows).sort_values("rank_metric_bn_seat_km", ascending=False)
    df.insert(0, "rank", range(1, len(df) + 1))
    df.to_csv(HERE / "routes.csv", index=False)

    # corridor table
    seg = trunks[["segment_order", "segment_type", "from_iata", "to_iata"]].copy()
    seg["distance_km"] = [round(gc_km(apt[a], apt[b])) for a, b in zip(trunks.from_iata, trunks.to_iata)]
    seg["pax_h_vs_subsonic"] = [round(seg_alloc[i][0]) for i in trunks.index]
    seg["pax_h_vs_boomless"] = [round(seg_alloc[i][1]) for i in trunks.index]
    seg["pax_h_per_km_vs_sub"] = (seg.pax_h_vs_subsonic / seg.distance_km).round(0)
    seg["pax_h_per_km_vs_bl"] = (seg.pax_h_vs_boomless / seg.distance_km).round(0)
    seg.to_csv(HERE / "corridor_segments.csv", index=False)

    # savings vs distance for an all-overland route (min useful length)
    curve = []
    for D in range(200, 4601, 50):
        dd = np.full(int(D / 5), 5.0 * PATH_FACTOR)
        ss = np.cumsum(dd) - dd / 2
        res = np.ones_like(dd, bool)
        t_sub, t_094, t_bl, t_m17, _ = overture_times(ss, dd, res)
        curve.append(dict(distance_km=D, t_sub=t_sub, t_094=t_094, t_bl=t_bl, t_m17=t_m17))
    cv = pd.DataFrame(curve)
    cv.to_csv(HERE / "savings_vs_distance.csv", index=False, float_format="%.4f")

    for col, lab in (("t_bl", "Boomless vs subsonic"), ("t_m17", "M1.7 vs subsonic")):
        sv = (cv.t_sub - cv[col]) * 60
        for thr in (15, 30, 60):
            hit = cv.distance_km[sv >= thr]
            print(f"{lab}: >= {thr} min from {hit.min() if len(hit) else 'n/a'} km")
    inc = (cv.t_bl - cv.t_m17) * 60
    for thr in (0.5, 10, 20, 30):
        hit = cv.distance_km[inc >= thr]
        print(f"M1.7 vs Boomless: >= {thr} min from {hit.min() if len(hit) else 'n/a'} km")

    if args.fill_template:
        fill_template(df)
    make_figure(df, cv)
    pd.set_option("display.width", 250)
    print(df[["rank", "origin", "dest", "distance_km", "land_fraction", "overland_fraction",
              "seats_equiv", "rank_metric_bn_seat_km", "time_subsonic_h", "time_boomless_h",
              "time_M17_h", "save_vs_subsonic_min", "save_vs_boomless_min", "shield_km",
              "pax_hours_saved_vs_boomless"]].to_string(index=False))
    print(seg.to_string(index=False))


def fill_template(df):
    f = NA / "north_america_top_pairs_template.csv"
    t = pd.read_csv(f)
    key = {frozenset((r.origin, r.dest)): r for r in df.itertuples()}
    for i, r in t.iterrows():
        k = key.get(frozenset((r.origin, r.dest)))
        if k is None:
            continue
        seats = k.seats_equiv
        t.loc[i, "annual_passengers"] = int(round(seats * 0.8, -3))
        t.loc[i, "reference_year"] = k.year
        tag = "OAG seats" if k.source != "EST" else "PLACEHOLDER estimate"
        t.loc[i, "notes"] = (f"{tag} {seats / 1e6:.2f}M x LF 0.8 (source {k.source}, "
                             f"see ../routes/routes.csv); replace with DB1B O&D")
    t["annual_passengers"] = t["annual_passengers"].astype("Int64")
    t.to_csv(f, index=False)


def make_figure(df, cv):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    blue, orange, ink, ink2, grid = "#2a78d6", "#eb6834", "#0b0b0b", "#52514e", "#e4e3df"
    plt.rcParams.update({"font.size": 9, "axes.edgecolor": grid, "axes.labelcolor": ink2,
                         "xtick.color": ink2, "ytick.color": ink2})
    fig, (a1, a2) = plt.subplots(1, 2, figsize=(12, 6.2), gridspec_kw={"width_ratios": [1.25, 1]})
    fig.patch.set_facecolor("#fcfcfb")
    top = df.head(20).iloc[::-1]
    lab = [f"{o}–{d}" + (" *" if s == "EST" else "") for o, d, s in zip(top.origin, top.dest, top.source)]
    y = np.arange(len(top))
    a1.barh(y, top.rank_metric_bn_seat_km, color=blue, height=0.7, edgecolor="#fcfcfb", linewidth=1)
    for yi, v, ph in zip(y, top.rank_metric_bn_seat_km, top.pax_hours_saved_vs_boomless):
        a1.text(v + 0.1, yi, f"{v:.1f}  ({ph / 1e3:.0f}k pax-h)", va="center", color=ink2, fontsize=8)
    a1.set_yticks(y, lab)
    a1.set_xlabel("overland km × annual seats (billion seat-km)")
    a1.set_title("Top 20 by overland seat-km\n(label: pax-h/yr a shield adds vs Boomless, 10% capture)",
                 loc="left", color=ink, fontsize=10)
    a1.grid(axis="x", color=grid, lw=0.6)
    a1.set_axisbelow(True)
    a1.set_xlim(0, top.rank_metric_bn_seat_km.max() * 1.35)
    for sp in ("top", "right"):
        a1.spines[sp].set_visible(False)

    x = cv.distance_km
    s_bl = (cv.t_sub - cv.t_bl) * 60
    s_17 = (cv.t_sub - cv.t_m17) * 60
    a2.plot(x, s_17, color=blue, lw=2, label="Overture M1.7 on shield")
    a2.plot(x, s_bl, color=orange, lw=2, label="Overture Boomless M1.15")
    a2.fill_between(x, s_bl, s_17, color=blue, alpha=0.12, lw=0)
    a2.text(x.iloc[-1], s_17.iloc[-1] + 3, "M1.7 shield", color=ink, ha="right", fontsize=8)
    a2.text(x.iloc[-1], s_bl.iloc[-1] - 9, "Boomless M1.15", color=ink, ha="right", fontsize=8)
    a2.text(3300, (s_bl.iloc[-25] + s_17.iloc[-25]) / 2, "what the\nshield buys", color=ink2, fontsize=8)
    for name, D in (("SFO–LAX", 543), ("LAX–ORD", 2802), ("LAX–JFK", 3974)):
        a2.axvline(D, color=grid, lw=0.8, zorder=0)
        a2.text(D + 30, 2, name, color=ink2, fontsize=7.5, rotation=90, va="bottom")
    a2.set_xlabel("all-overland great-circle distance (km)")
    a2.set_ylabel("block minutes saved vs Mach 0.85 jet")
    a2.set_title("Minutes saved vs route length", loc="left", color=ink, fontsize=10)
    a2.grid(color=grid, lw=0.6)
    a2.legend(frameon=False, loc="upper left")
    for sp in ("top", "right"):
        a2.spines[sp].set_visible(False)
    fig.text(0.01, 0.01, "* placeholder seat estimate. Screening model; see routes.md for assumptions.",
             color=ink2, fontsize=7.5)
    fig.tight_layout(rect=(0, 0.03, 1, 1))
    FIG.mkdir(exist_ok=True)
    fig.savefig(FIG / "routes_ranking.png", dpi=150, facecolor=fig.get_facecolor())


if __name__ == "__main__":
    main()
