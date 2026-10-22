#!/usr/bin/env python3
"""Skyloop system & economics re-review (screening model).

Run:  python3 sky_model.py            (one process, ~1-2 min; writes CSVs here, caches in .cache/)

Stages (each writes its CSV as soon as it is done, so a restart loses little):
  1. access.csv        access options to a 20 km tube: time per end, capacity, comment
  2. d2d.csv           door-to-door decomposition NYC-LAX, DEN-CHI at M2/3/4/6 vs subsonic, Boomless, BoomShield M2.2
  3. capacity.csv      headway (brick-wall vs moving block), pods/h, pax/yr per tube pair, vs spine demand
  4. cost_perkm.csv    bottom-up $/km (stratospheric run, troposphere ramp, urban approach, stations) low/base/high
     comparators.csv   HSR, Hyperloop studies, maglev, LIGO ($/km; evidence.csv)
  5. networks.csv      network evaluation on the coverage/demand O&D set (DOT 2025 + gravity fill), each comparator
     p2p_rank.csv      every point-to-point tube among the top-40 metros, ranked by benefit per $ of cost
  6. lifecycle.csv     life-cycle B/C (growth/ conventions: 7 %, 25 yr, 5-yr ramp-up, traffic + VOT growth,
                       central energy path) for openings 2045-2060, 30 % / 100 % capture, gross / net of energy+CO2
     breakeven.csv     $/km (run + ramps) at which B/C = 1

Reads read-only: boom-shield/network-maps/demand/.cache/inputs.pkl (O&D, block times, fares, time zones),
highmach/.cache/net_evals.pkl (BoomShield M2.2 spine per-pair block times), highmach/hm_aircraft.py (fuel),
highmach/hm_time.py (door-to-door phases), demand/build_demand.py (day-return rule, scenarios),
economics/growth/growth_model.py (traffic, VOT, energy paths).
[U] = engineering assumption.  Web figures: evidence.csv (search summaries, unverified).
"""
from __future__ import annotations

import importlib.util
import math
import pickle
import sys
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
BS = ROOT / "boom-shield"
NM = BS / "network-maps"
CACHE = HERE / ".cache"
CACHE.mkdir(exist_ok=True)
sys.dont_write_bytecode = True
for p in ("trunk", "coverage", "demand", "highmach"):
    sys.path.insert(0, str(NM / p))
import build_trunk as bt  # noqa: E402
import hm_aircraft as ha  # noqa: E402

_spec = importlib.util.spec_from_file_location("growth_model", BS / "economics" / "growth" / "growth_model.py")
gm = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(gm)

G = 9.81
A_MS = bt.A_MS                      # 295.1 m/s: "Mach" = nominal speed label inside the tube (repo convention)
V_TIME = 86.0                       # $/h (economics base)
AVAIL = 0.95
CUT = 10.0                          # min: smallest saving that counts (coverage convention)
MACHS = [2.0, 3.0, 4.0, 6.0]

# ----------------------------------------------------------------------------------------- base assumptions
P = dict(
    a_acc=0.3 * G,          # longitudinal accel/decel, seated & belted (physics study uses 0.3 g) [U]
    a_emerg=0.5 * G,        # emergency braking for brick-wall headway [U]
    t_react=2.0,            # s, detection + actuation [U]
    path=1.02,              # tube length / great circle (a straight tube; bends cost R = v^2/a_lat) [U]
    dwell=4.0,              # min, intermediate stop with doors sealed to a station lock [U]
    transfer_tube=20.0,     # min, change of pod at a shared station [U]
    transfer_air=45.0,      # min, airport feeder -> tube station (co-located) [U]
    drive_kmh=80.0, drive_detour=1.2, drive_max_km=250.0,
    # ground overhead at a ground portal: access / processing / egress (min)
    O_sky=(35, 40, 25),     # base 100 min: aviation-style screening, no gates/taxi; pods board in a lock [U]
    O_air=(40, 55, 35),     # highmach/demand base 130 min
    ramp_min_km=120.0,      # km: shortest climb 0 -> 20 km at a_n 0.1 g while accelerating (h ~ S a_n / 2a_t) [U]
    L_urb=20.0,             # km of urban approach (tunnel / viaduct) per station [U]
    seats=50, lf=0.80,      # pod seats; load factor [U]
    kg_per_seat=500.0,      # loaded pod mass per seat (Hyperloop-class pods 15-20 t / 28 pax) [U]
    eta_drive=0.85,         # grid -> pod (LSM + power electronics) [U]
    eta_regen=0.85,         # pod -> grid on braking (round trip ~0.72) [U]
    p_tube=10.0,            # Pa residual pressure (Hyperloop Alpha 100 Pa; physics study sets the real value) [U]
    cd_A=4.0,               # m^2, C_D x frontal area with blockage [U]
    LD_mag=300.0,           # magnetic lift/drag at speed (EDS) [U]
    hotel_kw_seat=0.4,      # kW per seat, life support + HVAC [U]
    pump_kw_km=20.0,        # kW per km of tube pair, steady leak pumping (outside p only 5.5 kPa) [U]
    peak_share=1 / 12,      # peak hour = 1/12 of the day's flow (8 %) [U]
    elec_usd_kwh=0.07,      # industrial power [U]
    airline_nonfuel=0.06,   # $/pax-km, avoided airline non-fuel cost (CASM ex-fuel ~9-10 c/seat-mi) [K][U]
    sky_ops_pax=15.0,       # $/pax, station staff + security + pod operations [U]
)
FUEL_SUB_PAX_KM = 0.090 / 3.16      # kg fuel per pax-km, subsonic (demand.md 90 g CO2)
CO2_PER_KG = 3.16

# grid carbon intensity, kg CO2/kWh (Ember 2024: 0.384; S14).  central: IRA-era decline [U]
GRID = {"stalled": [(2024, 0.384), (2050, 0.30), (2100, 0.25)],
        "central": [(2024, 0.384), (2035, 0.25), (2045, 0.15), (2060, 0.08), (2100, 0.05)],
        "clean": [(2024, 0.384), (2035, 0.10), (2040, 0.02), (2100, 0.02)]}


def grid_ci(case, y):
    pts = GRID[case]
    for (y0, v0), (y1, v1) in zip(pts, pts[1:]):
        if y0 <= y <= y1:
            return v0 + (v1 - v0) * (y - y0) / (y1 - y0)
    return pts[-1][1]


# ----------------------------------------------------------------------------------------- kinematics
def ride_min(L_km, M, a=None):
    """Station-to-station time on one tube (accelerate at a, cruise at M, brake at a), minutes."""
    a = a or P["a_acc"]
    V = M * A_MS
    L = np.asarray(L_km, dtype=float) * 1e3
    d2 = V * V / a                                   # accel + decel distance
    t = np.where(L >= d2, L / V + V / a, 2 * np.sqrt(np.maximum(L, 0) / a))
    return t / 60.0


def accel_km(M, a=None):
    a = a or P["a_acc"]
    return (M * A_MS) ** 2 / (2 * a) / 1e3


def ramp_km(M, a=None):
    return max(accel_km(M, a), P["ramp_min_km"])


def headway_s(M, rule):
    V = M * A_MS
    if rule == "brick-wall":
        return 1.1 * (V / (2 * P["a_emerg"]) + P["t_react"])
    if rule == "moving block 10 s":
        return 10.0
    if rule == "virtual coupling 3 s":
        return 3.0
    raise KeyError(rule)


def cap_pair_yr(M, rule, seats=None):
    """pax/yr a tube pair (one bore each way) can carry, both directions."""
    seats = seats or P["seats"]
    pods_h = 3600.0 / headway_s(M, rule)
    pax_h_dir = pods_h * seats * P["lf"]
    return pax_h_dir * 2 / P["peak_share"] * 365, pods_h, pax_h_dir


# ----------------------------------------------------------------------------------------- energy per pax
def kwh_per_pax(L_km, M, n_cycles=1):
    """Grid kWh per passenger for a tube trip of L_km with n_cycles accel/brake cycles."""
    V = M * A_MS
    m_pp = P["kg_per_seat"] / P["lf"]
    ke = 0.5 * m_pp * V * V
    trac = n_cycles * ke * (1 / P["eta_drive"] - P["eta_drive"] * P["eta_regen"])
    rho = P["p_tube"] / (287.0 * 250.0)
    drag_pod = 0.5 * rho * V * V * P["cd_A"]
    pax_pod = P["seats"] * P["lf"]
    aero = drag_pod * np.asarray(L_km) * 1e3 / pax_pod / P["eta_drive"]
    mag = m_pp * G / P["LD_mag"] * np.asarray(L_km) * 1e3 / P["eta_drive"]
    hotel = P["hotel_kw_seat"] / P["lf"] * 1e3 * (np.asarray(L_km) * 1e3 / V + 1800) / 1.0
    return (trac + aero + mag + hotel) / 3.6e6, dict(traction=trac / 3.6e6, aero=aero / 3.6e6, magnetic=mag / 3.6e6,
                                                     hotel=hotel / 3.6e6)


# ----------------------------------------------------------------------------------------- inputs
def load_inputs():
    f = CACHE / "inputs_sky.pkl"
    if f.exists():
        return pickle.loads(f.read_bytes())
    inp = pickle.loads((NM / "demand" / ".cache" / "inputs.pkl").read_bytes())
    df = inp["df"]
    n = len(df)
    iu, ju = inp["iu"], inp["ju"]
    b = inp["base"]
    TS = np.zeros((n, n))
    TS[iu, ju] = TS[ju, iu] = b["t_sub"]
    D = inp["D"]
    hm = pickle.loads((NM / "highmach" / ".cache" / "net_evals.pkl").read_bytes())
    e22 = hm["A: M2.2 Al-limit|spine"]
    t22 = b["t_bl"].copy()
    t22[e22["sel"]] = np.minimum(t22[e22["sel"]], e22["t"])
    out = dict(df=df, iu=iu, ju=ju, d=b["d"], t_sub=b["t_sub"], t_bl=b["t_bl"], t_m22=t22, pax=inp["pax"],
               fare=inp["fare"], tz=inp["tz"], TS=TS, D=D)
    f.write_bytes(pickle.dumps(out))
    return out


# fuel per passenger for each comparator aircraft (hm_network.fuel_per_pax conventions, LF 0.75)
_AC_REF = ha.size(1.7, "fixed", kT=1.22)          # Overture-class, flies Boomless at the reference
_AC_22 = ha.size(2.2, "adaptive", kT=1.40)


def fuel_ss(ac, km):
    L = np.asarray(km) * 1.03
    fa = ac["f_accel"]
    f = ac["W0_t"] * 1e3 * (fa + (1 - fa) * (1 - np.exp(-L / ac["RF_km"]))) * 0.97
    return f / (ha.SEATS * 0.75)


def comparator_world(inp, cmp):
    """Door-to-door minutes and fuel kg/pax per pair for comparator 'sub' | 'bl' | 'm22'."""
    O = sum(P["O_air"])
    d = inp["d"]
    f_sub = FUEL_SUB_PAX_KM * d * 1.03
    if cmp == "sub":
        return inp["t_sub"] + O, f_sub
    bl_used = inp["t_bl"] < inp["t_sub"] - 1.0
    f_bl = np.where(bl_used, fuel_ss(_AC_REF, d), f_sub)
    if cmp == "bl":
        return inp["t_bl"] + O, f_bl
    sh = inp["t_m22"] < inp["t_bl"] - CUT
    return np.where(sh, inp["t_m22"], inp["t_bl"]) + O, np.where(sh, fuel_ss(_AC_22, d), f_bl)


# ----------------------------------------------------------------------------------------- networks
def gc(df, a, b):
    r = df.set_index("code")
    return float(bt.gc(r.lat[a], r.lon[a], r.lat[b], r.lon[b]))


class Net:
    """lines: list of (station codes, mode) with mode 'allstop' | 'express'."""

    def __init__(self, name, lines, M, df):
        self.name, self.lines, self.M = name, lines, M
        self.st = sorted({s for ln, _ in lines for s in ln})
        self.k = {s: i for i, s in enumerate(self.st)}
        n = len(self.st)
        self.cum = []
        W = np.full((n, n), np.inf)
        how = {}
        for li, (ln, mode) in enumerate(lines):
            seg = [gc(df, a, b) * P["path"] for a, b in zip(ln, ln[1:])]
            cum = np.concatenate([[0], np.cumsum(seg)])
            self.cum.append((cum, seg))
            for i in range(len(ln)):
                for j in range(len(ln)):
                    if i == j:
                        continue
                    lo, hi = min(i, j), max(i, j)
                    if mode == "allstop":
                        t = float(sum(ride_min(seg[k], M) for k in range(lo, hi))) + P["dwell"] * (hi - lo - 1)
                        cyc = hi - lo
                    else:
                        t = float(ride_min(cum[hi] - cum[lo], M))
                        cyc = 1
                    a, b = self.k[ln[i]], self.k[ln[j]]
                    if t + P["transfer_tube"] < W[a, b]:
                        W[a, b] = t + P["transfer_tube"]
                        how[(a, b)] = (li, lo, hi, cum[hi] - cum[lo], cyc)
        # Floyd-Warshall with next hop (transfer penalty on every boarding; one removed at the end)
        R = W.copy()
        np.fill_diagonal(R, 0)
        nxt = {(a, b): b for (a, b) in how}
        for m in range(n):
            for a in range(n):
                for b in range(n):
                    if R[a, m] + R[m, b] < R[a, b]:
                        R[a, b] = R[a, m] + R[m, b]
                        nxt[(a, b)] = nxt[(a, m)]
        self.R = R - P["transfer_tube"]
        np.fill_diagonal(self.R, np.inf)
        self.how, self.nxt = how, nxt

    def legs(self, a, b):
        out = []
        while a != b:
            c = self.nxt[(a, b)]
            out.append(self.how[(a, c)])
            a = c
        return out

    def km(self):
        return sum(c[-1] for c, _ in self.cum)


def evaluate(inp, net, cmp, cap=1.0, O_sky=None, scen="base"):
    """Per-pair Skyloop response vs comparator.  Returns dict of arrays + totals (2025 demand)."""
    import build_demand as bd
    sc = bd.SCEN[scen]
    df, D, TS = inp["df"], inp["D"], inp["TS"]
    code = df.code.values
    O_sky = O_sky or P["O_sky"]
    O_air = P["O_air"]
    ns, nm = len(net.st), len(df)
    sidx = [int(np.where(code == s)[0][0]) for s in net.st]
    # access A[i, s] and egress E[s, j] (minutes, including terminal processing)
    A = np.full((nm, ns), np.inf)
    E = np.full((ns, nm), np.inf)
    for k, si in enumerate(sidx):
        dk = D[:, si]
        drive = 60 * dk * P["drive_detour"] / P["drive_kmh"]
        air_in = O_air[0] + O_air[1] + TS[:, si] + P["transfer_air"]
        air_out = P["transfer_air"] + TS[si, :] + O_air[2]
        A[:, k] = np.where(dk <= P["drive_max_km"], np.minimum(drive + O_sky[1] + 15, air_in), air_in)
        E[k, :] = np.where(dk <= P["drive_max_km"], np.minimum(drive + 10, air_out), air_out)
        A[si, k] = O_sky[0] + O_sky[1]
        E[k, si] = O_sky[2]
    iu, ju = inp["iu"], inp["ju"]
    best = np.full(len(iu), np.inf)
    bs = np.zeros(len(iu), int)
    bt_ = np.zeros(len(iu), int)
    for a in range(ns):
        for b in range(ns):
            if a == b or not np.isfinite(net.R[a, b]):
                continue
            # both directions: i->j uses (a,b); the reverse uses (b,a) with A/E swapped - symmetric enough
            t1 = A[iu, a] + net.R[a, b] + E[b, ju]
            t2 = A[ju, b] + net.R[b, a] + E[a, iu]
            t = 0.5 * (t1 + t2)
            upd = t < best
            best[upd], bs[upd], bt_[upd] = t[upd], a, b
    T0, fuel0 = comparator_world(inp, cmp)
    s = T0 - best
    use = s >= CUT
    s = np.where(use, s, 0.0)
    Q0 = inp["pax"]
    dz = inp["tz"][ju] - inp["tz"][iu]
    day, eve = bd.feas_gain(T0, np.where(use, best, T0), 0.0, dz, sc["W"], sc["Dmax"])
    day, eve = np.where(use, day, 0), np.where(use, eve, 0)
    credit = sc["d_day"] * day + sc["d_eve"] * eve
    F = inp["fare"]
    GC0 = F + V_TIME * T0 / 60
    GC1 = F + V_TIME * np.where(use, best, T0) / 60 - credit
    ratio = np.where(use, (GC1 / GC0) ** sc["eps"], 1.0)
    gross = np.where(use, cap * Q0 * np.maximum(ratio - 1, 0), 0.0)
    net_new = sc["n"] * gross
    capt = cap * Q0 * use
    time_usd = capt * s / 60 * V_TIME
    cred_usd = capt * credit
    ind_usd = 0.5 * gross * np.maximum(GC0 - GC1, 0)
    # tube km, cycles and loads
    tube_km = np.zeros(len(iu))
    cyc = np.zeros(len(iu))
    seg_load = [np.zeros(len(c[1])) for c in net.cum]
    pairs_used = np.where(use)[0]
    legs_cache = {}
    for p_ in pairs_used:
        key = (bs[p_], bt_[p_])
        if key not in legs_cache:
            legs_cache[key] = net.legs(*key)
        trips = capt[p_] + net_new[p_]
        for (li, lo, hi, km, cy) in legs_cache[key]:
            tube_km[p_] += km
            cyc[p_] += cy
            seg_load[li][lo:hi] += trips
    kwh = np.where(use, kwh_per_pax(tube_km, net.M, cyc)[0], 0.0)
    # feeder flights remain: fuel of the feeder legs (subsonic) for off-line ends
    d_feed = np.zeros(len(iu))
    for p_ in pairs_used:
        i, j = iu[p_], ju[p_]
        si, sj = sidx[bs[p_]], sidx[bt_[p_]]
        for x, y in ((i, si), (j, sj)):
            if x != y and D[x, y] > P["drive_max_km"]:
                d_feed[p_] += D[x, y]
    fuel_saved = capt * (fuel0 - FUEL_SUB_PAX_KM * d_feed * 1.03)          # kg/yr
    fuel_new = net_new * FUEL_SUB_PAX_KM * d_feed * 1.03                    # induced trips' feeder fuel
    kwh_tot = (capt + net_new) * kwh
    peak = max(float(x.max()) if len(x) else 0.0 for x in seg_load)
    return dict(T0=T0, T1=best, s=s, use=use, capt=capt, net_new=net_new, day=day, eve=eve,
                time_usd=time_usd, cred_usd=cred_usd, ind_usd=ind_usd, kwh=kwh_tot, fuel_saved=fuel_saved - fuel_new,
                tube_km=tube_km, seg_load=seg_load, peak=peak,
                tot=dict(captured_pax_M=capt.sum() / 1e6, net_induced_M=net_new.sum() / 1e6,
                         pax_h_M=(capt * s).sum() / 60e6, mean_save_min=(capt * s).sum() / max(capt.sum(), 1),
                         day_return_pax_M=(capt * day).sum() / 1e6,
                         benefit_busd=(time_usd.sum() + cred_usd.sum() + ind_usd.sum()) * AVAIL / 1e9,
                         time_only_busd=time_usd.sum() * AVAIL / 1e9,
                         fuel_saved_Mt=fuel_saved.sum() / 1e9 - fuel_new.sum() / 1e9,
                         kwh_TWh=kwh_tot.sum() / 1e9, pax_km_B=((capt + net_new) * tube_km).sum() / 1e9,
                         wh_per_pax_km=kwh_tot.sum() / max(((capt + net_new) * tube_km).sum(), 1) * 1e3,
                         peak_seg_pax_M=peak / 1e6))


# ----------------------------------------------------------------------------------------- costs
# $M per km of tube PAIR (two bores, one envelope), installed, before soft costs [U]; see system.md section 4
COST = {
    #                                low   base  high      (physics study: wall 35-60 kg/m per 3-4 m bore at 20 km)
    "core tubes (2 bores ~4 m, CFRP sandwich + joints/stiffening, 120-300 kg/m at $30-80/kg)": (3.0, 8.0, 25.0),
    "guideway, levitation track, alignment actuators": (3.0, 6.0, 12.0),
    "vacuum: pumps, sector gate valves, seals and joints": (1.0, 2.0, 5.0),
    "buoyancy: H2 envelope R 43/82/88 m (physics), film 0.32-0.66 M m2/km at $8/15/30 + gas": (2.6, 9.8, 20.5),
    "tethers, anchors, winches (physics: 45 MN every 3 km ~ 70 % of BoomShield t/km)": (3.0, 11.0, 22.0),
    "power and control along the line": (1.0, 2.0, 4.0),
    "safety: sector refuges, rescue vehicles, egress": (1.0, 2.0, 4.0),
}
SOFT = (0.15, 0.25, 0.35)            # engineering / PM
CONTING = (0.15, 0.30, 0.50)         # contingency (first of a kind)
RAMP_MULT = (1.5, 2.0, 3.0)          # troposphere crossing: jet-stream loads ~40x q, weather, icing, airspace [U]
LSM = (8.0, 15.0, 30.0)              # $M/km, linear motor + feeders on accel/brake sections (both bores) [U]
URBAN = (60.0, 150.0, 300.0)         # $M/km urban approach (tunnel / viaduct; Chuo ~$165-260 M/km) [U]
STATION = (1.0, 2.5, 5.0)            # $B per station: terminal, locks, power (~1 GW peak) [U]
SWITCH = (0.2, 0.5, 1.5)             # $B per high-speed switch (express intermediate station) [U]
POD_SEAT = (0.08, 0.15, 0.30)        # $M per pod seat [U]
OM_FRAC = (0.015, 0.025, 0.04)       # O&M per yr as share of capex (rail 1-2 %; BoomShield 4-6 %) [U]
LVL = {"low": 0, "base": 1, "high": 2}
# ground vactrain $M/km per tube pair: Great Lakes Hyperloop 34 / Hyperloop One own 62 / Chuo-like tunnels 207 (S1-S5)
GROUND = {"low": 34.0, "base": 62.0, "high": 207.0}


def run_cost(lvl):
    i = LVL[lvl]
    raw = sum(v[i] for v in COST.values())
    return raw * (1 + SOFT[i]) * (1 + CONTING[i])


def net_capex(net, lvl="base", bores=1, km_scale=1.0):
    """$B capex for a network: stratospheric run, ramps (+LSM), urban approaches, stations, switches, pods."""
    i = LVL[lvl]
    c_run = run_cost(lvl) * km_scale
    c_ramp = (run_cost(lvl) * RAMP_MULT[i] + LSM[i] * (1 + SOFT[i]) * (1 + CONTING[i])) * km_scale
    Lr = ramp_km(net.M)
    run_km = ramp_tot = 0.0
    n_int = 0
    for (ln, mode), (cum, seg) in zip(net.lines, net.cum):
        if mode == "allstop":
            r = sum(2 * min(s / 2, Lr) for s in seg)
        else:
            r = 2 * min(cum[-1] / 2, Lr)
            n_int += len(ln) - 2
        ramp_tot += r
        run_km += cum[-1] - r
    # express: each intermediate station needs 4 single-bore ramp tubes (= 2 pair-equivalents) and 2 switches
    ex_ramp = n_int * 2 * Lr
    tube = (run_km * c_run + (ramp_tot + ex_ramp) * c_ramp) * bores / 1e3
    urb = len(net.st) * P["L_urb"] * URBAN[i] / 1e3
    stn = len(net.st) * STATION[i] + n_int * 2 * SWITCH[i]
    return dict(tube_busd=tube, urban_busd=urb, station_busd=stn, total_ex_pods=tube + urb + stn,
                run_km=run_km, ramp_km=ramp_tot + ex_ramp, c_run=c_run, c_ramp=c_ramp, n_int=n_int)


def pods_capex(ev, net, lvl="base"):
    """Fleet sized from annual pax: peak-hour pods x round-trip cycle."""
    pax_yr = (ev["capt"] + ev["net_new"]).sum()
    if pax_yr <= 0:
        return 0.0
    mean_km = float(((ev["capt"] + ev["net_new"]) * ev["tube_km"]).sum() / pax_yr)
    pax_peak_h = pax_yr / 365 * P["peak_share"]
    pods_h = pax_peak_h / (P["seats"] * P["lf"])
    cycle_h = 2 * (float(ride_min(mean_km, net.M)) / 60 + 0.33)
    fleet = pods_h * cycle_h * 1.15
    return fleet * P["seats"] * POD_SEAT[LVL[lvl]] / 1e3


# ----------------------------------------------------------------------------------------- life cycle
def lifecycle(ann, capex_b, om_frac, open_year, capture=1.0, world="base", energy_case="central",
              grid_case="central", r=0.07, life=25, build=10, cap_pax=np.inf, resource=False):
    """ann: 2025 totals at 100 % capture (benefit_busd, fuel_saved_Mt, kwh_TWh, pax_M, pax_km_B, peak).
    Returns gross and net (of fuel, CO2, electricity) B/C and NPV(2026)."""
    tcase, vcase = gm.WORLDS[world]
    flows = []
    for y in range(open_year - build, open_year):
        flows.append((y, capex_b / build, 0.0, 0.0))
    for k, y in enumerate(range(open_year, open_year + life)):
        ramp = min(1.0, 0.5 + 0.5 * k / 5)
        T = gm.traffic_index(tcase, y)
        vol = capture * ramp * T
        capf = min(1.0, cap_pax / (ann["peak"] * vol)) if ann["peak"] > 0 else 1.0
        vol *= capf
        ben = ann["benefit_busd"] * vol * gm.vot_index(vcase, y)
        pf, scc, co2f, eff = gm.energy(energy_case, y)
        fuel_kg = ann["fuel_saved_Mt"] * 1e9 * vol * eff
        e_fuel = fuel_kg * (pf + CO2_PER_KG * co2f * scc / 1e3) / 1e9
        kwh = ann["kwh_TWh"] * 1e9 * vol
        e_el = kwh * (P["elec_usd_kwh"] + grid_ci(grid_case, y) * scc / 1e3) / 1e9
        extra = 0.0
        if resource:
            extra = (ann["pax_km_B"] * 1e9 * P["airline_nonfuel"] - ann["pax_M"] * 1e6 * P["sky_ops_pax"]) * vol / 1e9
        flows.append((y, capex_b * om_frac, ben, e_fuel - e_el + extra))
    C = sum(f[1] / (1 + r) ** (f[0] - open_year) for f in flows)
    B = sum(f[2] / (1 + r) ** (f[0] - open_year) for f in flows)
    N = sum(f[3] / (1 + r) ** (f[0] - open_year) for f in flows)
    d26 = (1 + r) ** (open_year - 2026)
    return dict(bc_gross=B / C, bc_net=(B + N) / C, npv26_gross=(B - C) / d26, npv26_net=(B + N - C) / d26,
                pv_cost=C, pv_ben=B, pv_energy=N)


# ----------------------------------------------------------------------------------------- stages
def stage_access():
    """Four ways to reach a tube at 20 km; minutes added per trip end beyond a ground-portal station."""
    rows = []
    h = 20e3
    for v, lab in ((10.0, "tether climber 10 m/s"), (20.5, "elevator at record 20.5 m/s (S12)"),
                   (50.0, "hypothetical 50 m/s pressurised lift")):
        ride = h / v / 60
        rows.append(dict(option="lift on a 20 km tower/tether", variant=lab, ride_min=ride,
                         per_end_min=ride + 8.0, both_ends_min=2 * (ride + 8.0),
                         capacity_pax_h=f"{60 / (2 * ride + 6) * 60:.0f} per 60-pax car per shaft",
                         blocker="no structure carries a guided car to 20 km: tallest building 0.83 km, rope "
                                 "lifts ~0.5-1 km (S12), 20 km inflatable tower is a patent (S13); tether sways "
                                 "~10-20 m in gusts and crosses the jet stream; cabin must be fully pressurised "
                                 "above 19 km (Armstrong limit)"))
    rows.append(dict(option="aircraft / VTOL docking at a 20 km station", variant="climb 25 min + dock 15 min",
                     ride_min=40.0, per_end_min=22 + 40 + 15, both_ends_min=2 * (22 + 40 + 15),
                     capacity_pax_h="one aircraft per dock cycle", blocker="keeps the whole airport overhead (130 min) "
                     "and adds taxi, climb to 20 km and an unprecedented passenger transfer at altitude; helicopters "
                     "and eVTOL ceilings are ~6 km, so it needs a high-altitude airliner, i.e. BoomShield's aircraft"))
    for M in MACHS:
        rows.append(dict(option="ground portal + climbing accel ramp", variant=f"M{M:g} at 0.3 g",
                         ride_min=0.0, per_end_min=0.0, both_ends_min=0.0,
                         capacity_pax_h="= line capacity",
                         blocker=f"ramp {ramp_km(M):.0f} km per end through the troposphere (accel {accel_km(M):.0f} km); "
                                 "ramp, LSM and urban approach dominate cost"))
    rows.append(dict(option="ground-level vactrain (no altitude)", variant="tube at grade/tunnel", ride_min=0.0,
                     per_end_min=0.0, both_ends_min=0.0, capacity_pax_h="= line capacity",
                     blocker="right-of-way, curvature (R = v^2/a: 1,400 km at M4, 0.1 g), tunnels; Hyperloop's problems"))
    A = pd.DataFrame(rows)
    A.to_csv(HERE / "access.csv", index=False, float_format="%.1f")
    return A


def stage_d2d(inp):
    import hm_time as hmt
    PAIRS = {"NYC-LAX": ("LGA", "LAX", (40.64, -73.78), (33.94, -118.41)),
             "DEN-CHI": ("DEN", "ORD", (39.86, -104.67), (41.98, -87.90))}
    rows = []
    O = P["O_air"]
    for pair, (a, b, p, q) in PAIRS.items():
        d = hmt.gc_km(p, q)
        for name, ph in (("subsonic", hmt.phases_conv(d, "sub")), ("Boomless M1.15", hmt.phases_conv(d, "bl"))):
            blk = sum(ph[k] for k in hmt.BLOCK)
            rows.append(dict(pair=pair, gc_km=d, mode=name, ground=sum(O), taxi_climb_desc=blk - ph["cruise"],
                             ramps=0.0, cruise=ph["cruise"], vertical=0.0, d2d_h=(blk + sum(O)) / 60))
        acc, dec, h = hmt.accel_for(2.2)
        ph = hmt.phases_shield(d, 2.2, acc, dec, h)
        blk = sum(ph[k] for k in hmt.BLOCK)
        rows.append(dict(pair=pair, gc_km=d, mode="BoomShield M2.2", ground=sum(O),
                         taxi_climb_desc=blk - ph["cruise"] - ph["accel"] - ph["decel"], ramps=ph["accel"] + ph["decel"],
                         cruise=ph["cruise"], vertical=0.0, d2d_h=(blk + sum(O)) / 60))
        L = d * P["path"]
        for M in MACHS:
            t = float(ride_min(L, M))
            rmp = min(2 * (M * A_MS) / P["a_acc"] / 60, t)
            for acc_name, Og, vert in (("ground portal, 100 min", sum(P["O_sky"]), 0.0),
                                       ("ground portal, walk-on 70 min", 70.0, 0.0),
                                       ("ground portal, airport-like 130 min", 130.0, 0.0),
                                       ("20 km station by lift (20.5 m/s)", sum(P["O_sky"]), 2 * (20e3 / 20.5 / 60 + 8)),
                                       ("20 km station by aircraft docking", 130.0, 2 * (22 + 40 + 15))):
                D2 = Og + t + vert
                rows.append(dict(pair=pair, gc_km=d, mode=f"Skyloop M{M:g} ({acc_name})", ground=Og,
                                 taxi_climb_desc=0.0, ramps=rmp, cruise=t - rmp, vertical=vert, d2d_h=D2 / 60,
                                 peak_mach=min(M, math.sqrt(P["a_acc"] * L * 1e3) / A_MS)))
    T = pd.DataFrame(rows)
    import build_demand as bd
    dz = {"NYC-LAX": -3, "DEN-CHI": +1}
    for sc in ("base", "high"):
        s = bd.SCEN[sc]
        T[f"day_return_{sc}"] = [bool(bd.day_window(np.array([r.d2d_h]), dz[r.pair], s["W"], s["Dmax"])[0]
                                      and bd.day_window(np.array([r.d2d_h]), -dz[r.pair], s["W"], s["Dmax"])[0])
                                 for r in T.itertuples()]
    T.to_csv(HERE / "d2d.csv", index=False, float_format="%.3f")
    return T


def stage_capacity():
    rows = []
    for M in MACHS:
        for rule in ("brick-wall", "moving block 10 s", "virtual coupling 3 s"):
            for seats in (4, 28, 50, 150):
                cap, pods_h, pax_h = cap_pair_yr(M, rule, seats)
                rows.append(dict(M=M, rule=rule, seats=seats, headway_s=headway_s(M, rule),
                                 stop_km=(M * A_MS) ** 2 / (2 * P["a_emerg"]) / 1e3, pods_h_dir=pods_h,
                                 pax_h_dir=pax_h, pax_yr_M=cap / 1e6, tube_pairs_for_spine_81M=81.4e6 / cap))
    C = pd.DataFrame(rows)
    C.to_csv(HERE / "capacity.csv", index=False, float_format="%.3g")
    return C


def stage_costs():
    rows = []
    for k, v in COST.items():
        rows.append(dict(item=k, low=v[0], base=v[1], high=v[2]))
    raw = [sum(v[i] for v in COST.values()) for i in range(3)]
    rows.append(dict(item="subtotal, stratospheric run", low=raw[0], base=raw[1], high=raw[2]))
    rows.append(dict(item="x soft costs x contingency = run $M/km", low=run_cost("low"), base=run_cost("base"),
                     high=run_cost("high")))
    for lvl in ("low", "base", "high"):
        pass
    rr = [run_cost(l) * RAMP_MULT[i] + LSM[i] * (1 + SOFT[i]) * (1 + CONTING[i]) for l, i in LVL.items()]
    rows.append(dict(item="troposphere ramp incl. LSM, $M/km", low=rr[0], base=rr[1], high=rr[2]))
    rows.append(dict(item="urban approach, $M/km", low=URBAN[0], base=URBAN[1], high=URBAN[2]))
    rows.append(dict(item="station, $B each", low=STATION[0], base=STATION[1], high=STATION[2]))
    rows.append(dict(item="pod, $M per seat", low=POD_SEAT[0], base=POD_SEAT[1], high=POD_SEAT[2]))
    rows.append(dict(item="O&M, share of capex per yr", low=OM_FRAC[0], base=OM_FRAC[1], high=OM_FRAC[2]))
    pd.DataFrame(rows).to_csv(HERE / "cost_perkm.csv", index=False, float_format="%.3g")
    comp = pd.DataFrame([
        dict(system="Hyperloop Alpha LA-SF (2013 claim)", usd_m_km_low=11, usd_m_km_high=13, src="S2"),
        dict(system="Virgin Hyperloop Missouri (2018, spokeswoman)", usd_m_km_low=19, usd_m_km_high=25, src="S3"),
        dict(system="Great Lakes Hyperloop (TEMS 2019)", usd_m_km_low=31, usd_m_km_high=37, src="S4"),
        dict(system="Hyperloop One Bay Area (own 2016 estimate)", usd_m_km_low=52, usd_m_km_high=75, src="S1"),
        dict(system="HSR China / Europe (World Bank)", usd_m_km_low=17, usd_m_km_high=39, src="S7"),
        dict(system="HSR California Merced-Bakersfield", usd_m_km_low=130, usd_m_km_high=175, src="S7"),
        dict(system="Shanghai Transrapid maglev (2004)", usd_m_km_low=40, usd_m_km_high=45, src="K1"),
        dict(system="Chuo Shinkansen maglev, 86 % tunnel", usd_m_km_low=165, usd_m_km_high=260, src="S5"),
        dict(system="LIGO beam tube, 1.2 m, 1e-9 torr (tube only)", usd_m_km_low=6.3, usd_m_km_high=8, src="S6"),
        dict(system="BoomShield spine integrated / base (capex per km)", usd_m_km_low=25, usd_m_km_high=186,
             src="integration"),
        dict(system="Skyloop doc 03 claim (tube + lift)", usd_m_km_low=5, usd_m_km_high=10, src="docs/03"),
        dict(system="Skyloop this study: stratospheric run", usd_m_km_low=run_cost("low"),
             usd_m_km_high=run_cost("high"), src="cost_perkm.csv"),
        dict(system="Skyloop this study: troposphere ramp", usd_m_km_low=rr[0], usd_m_km_high=rr[2],
             src="cost_perkm.csv"),
    ])
    comp["usd_m_km_mid"] = np.sqrt(comp.usd_m_km_low * comp.usd_m_km_high)
    comp.loc[comp.system.str.contains("this study: strat"), "usd_m_km_mid"] = run_cost("base")
    comp.loc[comp.system.str.contains("this study: tropo"), "usd_m_km_mid"] = rr[1]
    comp.to_csv(HERE / "comparators.csv", index=False, float_format="%.3g")
    return comp


SPINE = ["LGA", "CLE", "ORD", "OMA", "DEN", "LAS"]


def network_defs():
    return {
        "spine all-stop": [(SPINE, "allstop")],
        "spine express": [(SPINE, "express")],
        "spine+LAX all-stop": [(SPINE + ["LAX"], "allstop")],
        "NYC-LAX direct": [(["LGA", "LAX"], "express")],
        "NYC-CHI direct": [(["LGA", "ORD"], "express")],
        "DEN-CHI direct": [(["DEN", "ORD"], "express")],
        "LAX-SFO direct": [(["LAX", "SFO"], "express")],
        "SFO-LAX-LAS all-stop": [(["SFO", "LAX", "LAS"], "allstop")],
        "NYC-SFO direct": [(["LGA", "SFO"], "express")],
        "NYC-MIA direct": [(["LGA", "MIA"], "express")],
        "transcon pair NYC-LAX + NYC-SFO": [(["LGA", "LAX"], "express"), (["LGA", "SFO"], "express")],
        "transcon star NYC-LAX, NYC-SFO, CHI-LAX, NYC-CHI": [(["LGA", "LAX"], "express"), (["LGA", "SFO"], "express"),
                                                             (["ORD", "LAX"], "express"), (["LGA", "ORD"], "express")],
    }


def stage_networks(inp):
    f = CACHE / "net_rows.pkl"
    done = pickle.loads(f.read_bytes()) if f.exists() else {}
    df = inp["df"]
    for name, lines in network_defs().items():
        for M in MACHS:
            net = Net(name, lines, M, df)
            for cmp in ("sub", "bl", "m22"):
                key = (name, M, cmp)
                if key in done:
                    continue
                ev = evaluate(inp, net, cmp)
                cx = net_capex(net)
                pods = pods_capex(ev, net)
                capyr = cap_pair_yr(M, "moving block 10 s")[0]
                bores = int(min(6, max(1, math.ceil(ev["peak"] * 1.3 / capyr))))
                cx_b = net_capex(net, bores=bores)
                done[key] = dict(network=name, M=M, vs=cmp, km=net.km(), stations=len(net.st), **ev["tot"],
                                 bores_movblock=bores, cap_pair_movblock_M=capyr / 1e6,
                                 cap_pair_brickwall_M=cap_pair_yr(M, "brick-wall")[0] / 1e6,
                                 run_km=cx["run_km"], ramp_km=cx["ramp_km"], capex_1pair_busd=cx["total_ex_pods"] + pods,
                                 capex_busd=cx_b["total_ex_pods"] + pods, pods_busd=pods,
                                 tube_busd=cx_b["tube_busd"], urban_busd=cx_b["urban_busd"],
                                 station_busd=cx_b["station_busd"])
                f.write_bytes(pickle.dumps(done))
            print(f"  {name} M{M:g}: vs sub {done[(name, M, 'sub')]['benefit_busd']:.1f} $B/yr, "
                  f"capex {done[(name, M, 'sub')]['capex_busd']:.0f} $B", flush=True)
    N = pd.DataFrame(list(done.values()))
    # single-year annualised B/C (CRF 7 %, 25 yr) for orientation; life-cycle in stage_lifecycle
    crf = 0.07 / (1 - 1.07 ** -25)
    N["annual_cost_busd"] = N.capex_busd * (crf + OM_FRAC[1])
    N["bc_annual_gross"] = N.benefit_busd / N.annual_cost_busd
    N.to_csv(HERE / "networks.csv", index=False, float_format="%.4g")
    return N


def stage_p2p(inp, M=4.0, top=40):
    """Every direct tube among the top metros by O&D volume, standalone, vs subsonic and vs Boomless."""
    f = CACHE / f"p2p_M{M:g}.pkl"
    if f.exists():
        R = pickle.loads(f.read_bytes())
    else:
        R = []
    seen = {(r["a"], r["b"]) for r in R}
    df = inp["df"]
    vol = np.zeros(len(df))
    np.add.at(vol, inp["iu"], inp["pax"])
    np.add.at(vol, inp["ju"], inp["pax"])
    codes = df.code.values[np.argsort(-vol)[:top]]
    for x in range(len(codes)):
        for y in range(x + 1, len(codes)):
            a, b = codes[x], codes[y]
            if (a, b) in seen:
                continue
            d = gc(df, a, b)
            if d < 400:
                continue
            net = Net(f"{a}-{b}", [([a, b], "express")], M, df)
            row = dict(a=a, b=b, km=d * P["path"])
            for cmp in ("sub", "bl"):
                ev = evaluate(inp, net, cmp)
                row[f"ben_{cmp}"] = ev["tot"]["benefit_busd"]
                row[f"pax_{cmp}_M"] = ev["tot"]["captured_pax_M"]
                row[f"fuel_{cmp}_Mt"] = ev["tot"]["fuel_saved_Mt"]
                row[f"peak_{cmp}_M"] = ev["peak"] / 1e6
            cx = net_capex(net)
            row["capex_busd"] = cx["total_ex_pods"]
            R.append(row)
            if len(R) % 50 == 0:
                f.write_bytes(pickle.dumps(R))
    f.write_bytes(pickle.dumps(R))
    T = pd.DataFrame(R)
    crf = 0.07 / (1 - 1.07 ** -25)
    T["bc_sub"] = T.ben_sub / (T.capex_busd * (crf + OM_FRAC[1]))
    T["bc_bl"] = T.ben_bl / (T.capex_busd * (crf + OM_FRAC[1]))
    T = T.sort_values("bc_sub", ascending=False)
    T.to_csv(HERE / f"p2p_rank_M{M:g}.csv", index=False, float_format="%.4g")
    return T


def stage_lifecycle(N):
    rows, be = [], []
    for r in N.itertuples():
        ann = dict(benefit_busd=r.benefit_busd, fuel_saved_Mt=r.fuel_saved_Mt, kwh_TWh=r.kwh_TWh,
                   pax_M=r.captured_pax_M + r.net_induced_M, pax_km_B=r.pax_km_B, peak=r.peak_seg_pax_M * 1e6)
        build = 10 if r.km > 2000 else 7
        cap_pax = r.bores_movblock * r.cap_pair_movblock_M * 1e6
        for oy in (2045, 2050, 2055, 2060):
            for capture in (0.3, 1.0):
                for lvl in ("low", "base", "high"):
                    i = LVL[lvl]
                    # rescale capex by cost level: tube + urban + station scale with their tables; pods scale
                    scale = None
                    capex = capex_level(r, lvl)
                    lc = lifecycle(ann, capex, OM_FRAC[i], oy, capture, build=build, cap_pax=cap_pax)
                    lcr = lifecycle(ann, capex, OM_FRAC[i], oy, capture, build=build, cap_pax=cap_pax, resource=True)
                    rows.append(dict(network=r.network, M=r.M, vs=r.vs, opening=oy, capture=capture, cost=lvl,
                                     capex_busd=capex, bc_gross=lc["bc_gross"], bc_net=lc["bc_net"],
                                     bc_net_plus_airline_cost=lcr["bc_net"], npv26_net_busd=lc["npv26_net"],
                                     pv_energy_busd=lc["pv_energy"], pv_ben_busd=lc["pv_ben"], pv_cost_busd=lc["pv_cost"]))
                for gl, cg in GROUND.items():
                    capex = (r.run_km + r.ramp_km) * cg * r.bores_movblock / 1e3 + r.urban_busd + r.station_busd \
                        + r.pods_busd
                    lc = lifecycle(ann, capex, OM_FRAC[1], oy, capture, build=build, cap_pax=cap_pax)
                    rows.append(dict(network=r.network, M=r.M, vs=r.vs, opening=oy, capture=capture,
                                     cost=f"ground vactrain {gl}", capex_busd=capex, bc_gross=lc["bc_gross"],
                                     bc_net=lc["bc_net"], bc_net_plus_airline_cost=np.nan,
                                     npv26_net_busd=lc["npv26_net"], pv_energy_busd=lc["pv_energy"],
                                     pv_ben_busd=lc["pv_ben"], pv_cost_busd=lc["pv_cost"]))
                # break-even: scale the per-km (run + ramp) cost so that net B/C = 1, base level, stations fixed
                fixed = r.urban_busd + r.station_busd + r.pods_busd
                var = r.tube_busd
                lc1 = lifecycle(ann, fixed + var, OM_FRAC[1], oy, capture, build=build, cap_pax=cap_pax)
                lc0 = lifecycle(ann, fixed, OM_FRAC[1], oy, capture, build=build, cap_pax=cap_pax)
                for kind, numer in (("gross", lc1["pv_ben"]), ("net", lc1["pv_ben"] + lc1["pv_energy"])):
                    # PV cost is linear in capex: C(k) = C0 + k (C1 - C0)
                    k = (numer - lc0["pv_cost"]) / max(lc1["pv_cost"] - lc0["pv_cost"], 1e-9)
                    km_tot = r.run_km + r.ramp_km
                    be.append(dict(network=r.network, M=r.M, vs=r.vs, opening=oy, capture=capture, kind=kind,
                                   scale_vs_base=k, breakeven_run_usd_m_km=k * run_cost("base"),
                                   breakeven_avg_tube_usd_m_km=k * var * 1e3 / max(km_tot * r.bores_movblock, 1),
                                   base_avg_tube_usd_m_km=var * 1e3 / max(km_tot * r.bores_movblock, 1)))
    L = pd.DataFrame(rows)
    L.to_csv(HERE / "lifecycle.csv", index=False, float_format="%.4g")
    B = pd.DataFrame(be)
    B.to_csv(HERE / "breakeven.csv", index=False, float_format="%.4g")
    return L, B


def capex_level(r, lvl):
    i = LVL[lvl]
    k_run = run_cost(lvl) / run_cost("base")
    return (r.tube_busd * k_run + r.urban_busd * URBAN[i] / URBAN[1] + r.station_busd * STATION[i] / STATION[1]
            + r.pods_busd * POD_SEAT[i] / POD_SEAT[1])


if __name__ == "__main__":
    pd.set_option("display.width", 250)
    pd.set_option("display.max_columns", 40)
    inp = load_inputs()
    print(stage_access()[["option", "variant", "both_ends_min"]].to_string())
    T = stage_d2d(inp)
    print(T[["pair", "mode", "ground", "taxi_climb_desc", "ramps", "cruise", "vertical", "d2d_h", "day_return_base"]]
          .round(2).to_string())
    C = stage_capacity()
    print(C[C.seats == 50][["M", "rule", "headway_s", "pods_h_dir", "pax_h_dir", "pax_yr_M",
                             "tube_pairs_for_spine_81M"]].round(2).to_string())
    print(stage_costs().round(1).to_string())
    print("kWh/pax NYC-LAX:", {M: round(float(kwh_per_pax(4050, M)[0]), 1) for M in MACHS},
          kwh_per_pax(4050, 4.0)[1])
    N = stage_networks(inp)
    print(N[["network", "M", "vs", "km", "captured_pax_M", "pax_h_M", "mean_save_min", "day_return_pax_M",
             "benefit_busd", "fuel_saved_Mt", "kwh_TWh", "wh_per_pax_km", "peak_seg_pax_M", "bores_movblock",
             "capex_busd", "bc_annual_gross"]].round(2).to_string())
    L, B = stage_lifecycle(N)
    print("lifecycle written", len(L), len(B))
    if "--p2p" in sys.argv:
        print(stage_p2p(inp).head(25).round(3).to_string())
