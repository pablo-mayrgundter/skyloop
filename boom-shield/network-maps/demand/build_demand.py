#!/usr/bin/env python3
"""
BoomShield induced demand: how much new traffic does a shielded supersonic trunk create,
and what does it do to B/C?

Reads (read-only) ../coverage (cached demand matrix, Boomless/subsonic baselines, network
evaluator), ../trunk/build_trunk.py (time model, segment loads) and
../../economics/cost_model.py (value of time, annual cost per km).  Writes CSVs + figures here.

Model, per O&D pair k with base demand Q0 (DOT 2025 Table 1 + gravity fill):
  generalised cost  GC = fare + V * (block + O) / 60       (V = $86/h, O = access/egress overhead)
  reference world   r = Boomless Cruise (fare F_X) or subsonic (fare F_S, DOT Table 1 fare)
  shield world      captured share c of the pair rides the trunk: block t_net, fare F_X (+ toll)
                    minus threshold credits where a same-day return / after-work arrival
                    becomes feasible (time-zone aware, per direction)
  gross route response   dQ = c Q0 [(GC1/GC0)^eps - 1]   (constant elasticity, route level, >= 0)
  net new trips          n * dQ ; the rest (1-n) is diverted from other city pairs
  benefit                c Q0 V dt (as coverage.md, time only) + c Q0 credits + 1/2 dQ (GC0 - GC1)
                         (rule of half; diversion does not change it to first order)
Scenarios low/base/high vary eps, n, supersonic fare multiple phi, overhead O and the
threshold parameters; 'TTB bound' is a travel-time-budget limit (trips ~ 1/door-to-door time).

Usage: python3 build_demand.py   (one process, ~4 min cold; caches in .cache/)
"""
from __future__ import annotations

import json
import math
import pickle
import sys
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
NM = HERE.parent
CACHE = HERE / ".cache"
FIG = HERE / "figures"
sys.path.insert(0, str(NM / "coverage"))
import build_coverage as bc  # noqa: E402  (imports trunk/build_trunk.py as bc.bt and cost_model as bc.cm)

bt, cm = bc.bt, bc.cm
CACHE.mkdir(exist_ok=True)
FIG.mkdir(exist_ok=True)

V = cm.VTTS["base"]                         # $86/h (economics.md base)
AVAIL = cm.v("availability", "base")        # 0.95
SEAT_YR = bt.SEATS * bt.LOAD_FACTOR * 365   # pax/yr per daily Overture flight
LID_CAP = {"low": 1000, "base": 10000}      # fatigue/lid-fatigue.md revised caps, flights/day

SCEN = {
    #        eps = route-level GC elasticity; n = net-new share; phi = supersonic fare / today's fare
    #        O = door-to-door overhead (min); W = hours at destination for a day return;
    #        d_day / d_eve = $ credit per trip where a day return / after-work arrival becomes feasible
    #        Dmax = longest door-to-door leg (h) people accept for a day return
    "low":  dict(eps=-0.8, n=0.50, phi=5.0, O=150, W=6.0, Dmax=3.5, d_day=0.0, d_eve=0.0),
    "base": dict(eps=-1.3, n=0.60, phi=3.0, O=130, W=5.0, Dmax=4.0, d_day=60.0, d_eve=20.0),
    "high": dict(eps=-2.0, n=0.75, phi=1.5, O=100, W=4.0, Dmax=4.5, d_day=100.0, d_eve=30.0),
}
GROWTH = 1.40      # background traffic over asset life (~2 %/yr, life-average) [U] - separate from induced
CO2_SUB = 0.090    # kg CO2 per pax-km, US domestic narrowbody [K]
CO2_MULT = {"low": 3.0, "base": 5.0, "high": 7.0}   # ICCT 2018 supersonic / subsonic per pax
SCC = 190.0        # $/t CO2 (EPA 2023, 2 %) [K]


def cached(name, fn):
    f = CACHE / f"{name}.pkl"
    if f.exists():
        return pickle.loads(f.read_bytes())
    out = fn()
    f.write_bytes(pickle.dumps(out))
    return out


def cov_cached(name):
    return pickle.loads((bc.CACHE / f"{name}.pkl").read_bytes())


# =====================================================================================
# 1. Inputs: demand, baselines, fares, time zones
# =====================================================================================
def tz_of(lon):
    """Standard-time offset from longitude bands (rough US zones; DST and Arizona ignored)."""
    lon = np.asarray(lon)
    return np.select([lon > -86.5, lon > -102.0, lon > -114.5], [-5, -6, -7], -8).astype(float)


def stage_inputs():
    dem, base, grid = cov_cached("demand"), cov_cached("baseline"), cov_cached("grid")
    df, D = dem["df"], dem["D"]
    T = dem["T_obs"] + dem["T_fill"]
    iu, ju = base["iu"], base["ju"]
    # fares: DOT Table 1 average one-way fare; fill pairs from log-log fit on distance
    g, _, _ = bc.load_table1(2025)
    pos = {c: i for i, c in enumerate(df.code)}
    city2code = dict(zip(df.city, df.code))
    Fm = np.full(D.shape, np.nan)
    for c1, c2, f in zip(g.city1, g.city2, g.fare):
        i, j = pos[city2code[c1]], pos[city2code[c2]]
        Fm[i, j] = Fm[j, i] = f
    fare = Fm[iu, ju]
    obs = np.isfinite(fare)
    b1, b0 = np.polyfit(np.log(base["d"][obs]), np.log(fare[obs]), 1)
    fit = np.exp(b0 + b1 * np.log(base["d"]))
    fare = np.where(obs, fare, fit)
    print(f"fares: {obs.sum()} Table 1 pairs, fit fare = {math.exp(b0):.1f} * km^{b1:.3f}; "
          f"pax-weighted mean {np.average(fare, weights=T[iu, ju]):.0f} $")
    tz = tz_of(df.lon.values)
    return dict(df=df, D=D, T=T, base=base, grid=grid, fare=fare, tz=tz, iu=iu, ju=ju,
                pax=T[iu, ju], fare_fit=(b0, b1))


# =====================================================================================
# 2. Networks (re-evaluated with coverage's evaluator; cached)
# =====================================================================================
def networks():
    cc = pd.read_csv(NM / "coverage" / "coverage_curve.csv")
    cc = cc[(cc.path == "greedy from spine") & (cc.mach == 2.0)].sort_values("step")
    nets = [("spine", 1.7, bc.SPINE), ("spine", 2.0, bc.SPINE)]
    for name, edges in bc.SCEN.items():
        if name != "spine":
            nets.append((name, 2.0, edges))
    for r in cc.itertuples():
        if r.step > 0:
            nets.append((f"greedy step {r.step}", 2.0, [tuple(e.split("-")) for e in r.edges.split()]))
    return nets


def stage_net_evals(inp):
    dfx = bc.ext_df(inp["df"], inp["grid"])
    sel = np.where(inp["base"]["t_sub"] - np.minimum(inp["base"]["t_ub_M1.7"], inp["base"]["t_ub_M2.0"])
                   >= bc.CUT_MIN)[0]
    f = CACHE / "net_evals.pkl"
    out = pickle.loads(f.read_bytes()) if f.exists() else {}
    for name, M, edges in networks():
        key = f"{name}|M{M}"
        if key in out:
            continue
        ne = bc.net_eval(dfx, inp["base"], inp["T"], edges, M, sel)
        out[key] = dict(ne, edges=edges, mach=M, name=name)
        f.write_bytes(pickle.dumps(out))          # save as we go
        print(f"  evaluated {key}: {ne['km']:.0f} km, {ne['B_bl'] / 1e6:.1f} M pax-h")
    return out


# =====================================================================================
# 3. Demand response
# =====================================================================================
def day_window(D_h, dz, W, Dmax):
    """Hours available at destination on a same-day return (leave home 06:00, home by 23:00
    origin local; business window 08:30-18:00 destination local).  dz = dest - origin tz."""
    start = np.maximum(6.0 + D_h + dz, 8.5)
    end = np.minimum(23.0 - D_h + dz, 18.0)
    return ((end - start) >= W) & (D_h <= Dmax)


def evening_ok(D_h, dz):
    """Leave work 17:30 origin local and arrive at the destination door by 00:30 local."""
    return 17.5 + D_h + dz <= 24.5


def feas_gain(T0, T1, O, dz, W, Dmax):
    """Fraction of the pair's directional trips (0, 0.5, 1) newly gaining day-return /
    after-work-arrival feasibility."""
    d0, d1 = (T0 + O) / 60, (T1 + O) / 60
    day = sum((day_window(d1, s * dz, W, Dmax) & ~day_window(d0, s * dz, W, Dmax)).astype(float) for s in (1, -1)) / 2
    eve = sum((evening_ok(d1, s * dz) & ~evening_ok(d0, s * dz)).astype(float) for s in (1, -1)) / 2
    return day, eve


def respond(inp, ne, cmp, sc, cap, growth=1.0, ttb=False, toll=0.0):
    """Per-pair response for network ne vs reference cmp ('bl' or 'sub') at capture cap."""
    base, sel = inp["base"], ne["sel"]
    iu, ju = inp["iu"][sel], inp["ju"][sel]
    Q0 = ne["pax"] * growth
    s = ne["s_bl"] if cmp == "bl" else ne["s_sub"]
    use = s > 0
    T0 = base["t_bl"][sel] if cmp == "bl" else base["t_sub"][sel]
    T1 = np.where(use, T0 - s, T0)
    FS = inp["fare"][sel]
    FX = sc["phi"] * FS
    F0 = FX if cmp == "bl" else FS
    O = sc["O"]
    dz = inp["tz"][ju] - inp["tz"][iu]
    day, eve = feas_gain(T0, T1, O, dz, sc["W"], sc["Dmax"])
    day, eve = np.where(use, day, 0), np.where(use, eve, 0)
    credit = sc["d_day"] * day + sc["d_eve"] * eve
    GC0 = F0 + V * (T0 + O) / 60
    GC1 = FX + toll + V * (T1 + O) / 60 - credit
    if ttb:      # travel-time budget bound: trips ~ 1 / door-to-door time, fares ignored, all net new
        ratio_q = (T0 + O) / (T1 + O)
        dGC = V * (T0 - T1) / 60
    else:
        ratio_q = np.where(use, (GC1 / GC0) ** sc["eps"], 1.0)
        dGC = GC0 - GC1
    gross = np.where(use & (dGC > 0), cap * Q0 * np.maximum(ratio_q - 1, 0), 0.0)
    n = 1.0 if ttb else sc["n"]
    ex_h = cap * Q0 * s / 60                                  # existing captured pax-h (coverage.md metric)
    cred_usd = cap * Q0 * credit
    ind_usd = 0.5 * gross * np.maximum(dGC, 0)                # rule of half
    return dict(sel=sel, iu=iu, ju=ju, Q0=Q0, use=use, s=s, gross=gross, net=n * gross,
                ex_h=ex_h, cred_usd=cred_usd, ind_usd=ind_usd, day=day, eve=eve, GC0=GC0, GC1=GC1,
                T0=T0, T1=T1, cap_pax=cap * Q0 * use)


def summarise(inp, ne, r, cap):
    ann = bt.annual_cost_per_km(ne["mach"], "base")[0]
    cost = ann * ne["km"]
    ex_usd = r["ex_h"].sum() * V
    ben = (ex_usd + r["cred_usd"].sum() + r["ind_usd"].sum()) * AVAIL
    fixed = ex_usd * AVAIL
    return dict(km=ne["km"], captured_pax_M=r["cap_pax"].sum() / 1e6, gross_induced_M=r["gross"].sum() / 1e6,
                net_induced_M=r["net"].sum() / 1e6, induced_pct=r["gross"].sum() / max(r["cap_pax"].sum(), 1) * 100,
                pax_h_fixed_M=r["ex_h"].sum() / 1e6,
                pax_h_equiv_total_M=(ex_usd + r["cred_usd"].sum() + r["ind_usd"].sum()) / V / 1e6,
                credit_share=r["cred_usd"].sum() / max(ben / AVAIL, 1), induced_share=r["ind_usd"].sum() / max(ben / AVAIL, 1),
                bcr_fixed=fixed / cost, bcr_plus_induced=(ex_usd + r["ind_usd"].sum()) * AVAIL / cost, bcr=ben / cost, uplift=ben / max(fixed, 1e-9),
                pax_day_return_M=(r["cap_pax"] * r["day"]).sum() / 1e6,
                pax_evening_M=(r["cap_pax"] * r["eve"]).sum() / 1e6)


# =====================================================================================
# main
# =====================================================================================
if __name__ == "__main__":
    inp = cached("inputs", stage_inputs)
    nets = stage_net_evals(inp)
    code = inp["df"].code.values
    variants = [(k, SCEN[k], {}) for k in SCEN] + [("TTB bound", SCEN["base"], dict(ttb=True)),
                                                    ("base + life-avg growth x1.4", SCEN["base"], dict(growth=GROWTH)),
                                                    ("base, toll $25", SCEN["base"], dict(toll=25.0))]

    # ---------------------------------------------------------------- 3a. scenario table, all networks
    rows = []
    for key, ne in nets.items():
        for vname, sc, kw in variants:
            for cmp in ("bl", "sub"):
                for cap in (0.1, 1.0):
                    r = respond(inp, ne, cmp, sc, cap, **kw)
                    rows.append(dict(network=ne["name"], mach=ne["mach"], scenario=vname, vs=cmp, capture=cap,
                                     **summarise(inp, ne, r, cap)))
    R = pd.DataFrame(rows)
    R.to_csv(HERE / "scenario_results.csv", index=False, float_format="%.5g")
    sp = R[(R.network == "spine") & (R.mach == 2.0)]
    print(sp[["scenario", "vs", "capture", "captured_pax_M", "gross_induced_M", "net_induced_M", "induced_pct",
              "pax_h_fixed_M", "pax_h_equiv_total_M", "bcr_fixed", "bcr", "uplift"]].round(4).to_string())

    # ---------------------------------------------------------------- 3b. per-pair detail, spine M2.0 base
    ne_sp = nets["spine|M2.0"]
    rb = respond(inp, ne_sp, "bl", SCEN["base"], 1.0)
    rs = respond(inp, ne_sp, "sub", SCEN["base"], 1.0)
    rh = respond(inp, ne_sp, "bl", SCEN["high"], 1.0)
    P = pd.DataFrame(dict(a=code[rb["iu"]], b=code[rb["ju"]], km=inp["base"]["d"][rb["sel"]].round(0),
                          pax_yr=rb["Q0"].round(0), fare_sub=inp["fare"][rb["sel"]].round(0),
                          t_sub=inp["base"]["t_sub"][rb["sel"]].round(1), t_bl=rb["T0"].round(1),
                          t_spine=np.where(rb["use"], rb["T1"], np.nan).round(1),
                          save_bl=rb["s"].round(1), save_sub=rs["s"].round(1),
                          gc_bl=rb["GC0"].round(0), gc_shield=rb["GC1"].round(0),
                          induced_bl_base=rb["gross"].round(0), induced_bl_high=rh["gross"].round(0),
                          induced_sub_base=rs["gross"].round(0),
                          day_return_gain=rb["day"], evening_gain=rb["eve"]))
    P = P[rb["use"] | rs["use"]].sort_values("pax_yr", ascending=False)
    P.to_csv(HERE / "pairs_induced.csv", index=False)
    print(P.head(15).to_string())

    # threshold check on named pairs (door-to-door hours, base overhead)
    O = SCEN["base"]["O"]
    chk = []
    for a, b in [("LGA", "LAX"), ("LGA", "SFO"), ("LGA", "LAS"), ("LGA", "DEN"), ("ORD", "LAX"), ("ORD", "LAS"),
                 ("BOS", "LAX"), ("LGA", "PHX"), ("DCA", "LAX"), ("LGA", "SEA")]:
        m = ((P.a == a) & (P.b == b)) | ((P.a == b) & (P.b == a))
        if m.any():
            q = P[m].iloc[0]
            chk.append(dict(pair=f"{a}-{b}", km=q.km, d2d_sub_h=(q.t_sub + O) / 60, d2d_bl_h=(q.t_bl + O) / 60,
                            d2d_spine_h=(q.t_spine + O) / 60, day_return_gain=q.day_return_gain,
                            evening_gain=q.evening_gain))
    CH = pd.DataFrame(chk)
    CH.round(2).to_csv(HERE / "threshold_check.csv", index=False)
    print(CH.round(2).to_string())

    # ---------------------------------------------------------------- 3c. segment flows (spine, shoulder)
    knee = json.load(open(NM / "coverage" / "knee.json"))
    shoulder_key = f"greedy step {knee['shoulder_step']}|M2.0"
    seg_rows = []
    flights = []
    for key in ("spine|M2.0", shoulder_key):
        ne = nets[key]
        tr = ne["tr"]
        for vname in ("fixed", "low", "base", "high", "base + growth"):
            sc = SCEN["base" if vname in ("fixed", "base + growth") else vname]
            r = respond(inp, ne, "bl", sc, 1.0, growth=GROWTH if vname == "base + growth" else 1.0)
            pax = r["cap_pax"] + (0 if vname == "fixed" else r["gross"])
            ev = dict(pax=pax, sav=r["s"], on=ne["on"], off=ne["off"], use=r["use"])
            lb = bt.segment_loads(tr, ev, "sav", "use")
            agg = {}
            for (u, v), (p, ph) in lb.items():
                ei = tr.G[u][v]["edge"]
                a_, b_ = tr.edges[ei]
                rec = agg.setdefault(f"{a_}-{b_}", [0.0, 0.0])
                rec[0] = max(rec[0], p)
            for eg, (pmax, _) in agg.items():
                seg_rows.append(dict(network=ne["name"], variant=vname, edge=eg, peak_pax_M=pmax / 1e6,
                                     flights_day=pmax / SEAT_YR))
            pk = max(v[0] for v in lb.values())
            flights.append(dict(network=ne["name"], variant=vname, peak_seg_pax_M=pk / 1e6,
                                peak_flights_day=pk / SEAT_YR,
                                share_of_cap_10000=pk / SEAT_YR / LID_CAP["base"],
                                share_of_cap_1000=pk / SEAT_YR / LID_CAP["low"]))
    S = pd.DataFrame(seg_rows)
    Sw = S.pivot_table(index=["network", "edge"], columns="variant", values="peak_pax_M").reset_index()
    Sw.to_csv(HERE / "segment_flows.csv", index=False, float_format="%.3f")
    FL = pd.DataFrame(flights)
    FL.to_csv(HERE / "flights_vs_caps.csv", index=False, float_format="%.4g")
    print(Sw.round(1).to_string())
    print(FL.round(3).to_string())

    # ---------------------------------------------------------------- 3d. airports (spine, base, 100 %)
    r = respond(inp, ne_sp, "bl", SCEN["base"], 1.0)
    tot = r["cap_pax"] + r["gross"]
    ap = pd.DataFrame(dict(a=np.r_[code[r["iu"]], code[r["ju"]]], pax=np.r_[tot, tot] / 2,
                           pax_fixed=np.r_[r["cap_pax"], r["cap_pax"]] / 2))
    ap = ap.groupby("a").sum().sort_values("pax", ascending=False)
    ap["overture_deps_day"] = ap.pax / SEAT_YR            # departures/day (each pax-trip departs once)
    ap["narrowbody_deps_day_same_pax"] = ap.pax_fixed / (150 * 0.85 * 365)
    ap["overture_movements_day"] = 2 * ap.overture_deps_day
    MOV_K = {"LGA": 3400, "ORD": 2400, "LAX": 1700, "DEN": 1800, "LAS": 1600, "SFO": 1300, "BOS": 1100,
             "DCA": 2400, "SEA": 1150, "PHX": 1200, "PHL": 900, "MSP": 1000}   # metro movements/day 2024 [K][U]
    ap["metro_movements_day_2024_K"] = [MOV_K.get(c, np.nan) for c in ap.index]
    ap["added_share"] = (ap.overture_movements_day - 2 * ap.narrowbody_deps_day_same_pax) / ap.metro_movements_day_2024_K
    ap.head(15).to_csv(HERE / "airports.csv", float_format="%.3g")
    print(ap.head(12).round(2).to_string())

    # ---------------------------------------------------------------- 3e. emissions flag (spine, 100 %)
    em = []
    for vname in ("low", "base", "high"):
        sc = SCEN[vname]
        for cmp in ("bl", "sub"):
            r = respond(inp, ne_sp, cmp, sc, 1.0)
            dkm = inp["base"]["d"][r["sel"]] * bt.PATH
            ind_pkm = (r["gross"] * sc["n"] * dkm).sum()
            co2_ind = ind_pkm * CO2_SUB * CO2_MULT[vname] / 1e3       # t
            sw_pkm = (r["cap_pax"] * dkm).sum() if cmp == "sub" else 0.0
            co2_sw = sw_pkm * CO2_SUB * (CO2_MULT[vname] - 1) / 1e3
            ben = (r["ex_h"].sum() * V + r["cred_usd"].sum() + r["ind_usd"].sum()) * AVAIL
            em.append(dict(scenario=vname, vs=cmp, induced_net_pax_km_bn=ind_pkm / 1e9, co2_induced_Mt=co2_ind / 1e6,
                           co2_switch_to_supersonic_Mt=co2_sw / 1e6, carbon_cost_busd=(co2_ind + co2_sw) * SCC / 1e9,
                           benefit_busd=ben / 1e9, carbon_cost_over_benefit=(co2_ind + co2_sw) * SCC / ben))
    EM = pd.DataFrame(em)
    EM.to_csv(HERE / "emissions_flag.csv", index=False, float_format="%.4g")
    print(EM.round(3).to_string())

    # ---------------------------------------------------------------- 3f. fare premium -> endogenous capture
    # lognormal value of time (median $60/h, sigma 0.75, mean ~$80/h): choose the shield over subsonic if
    # VOT * saving > fare premium.  Pax-weighted over spine-captured pairs.
    from scipy.stats import norm
    rr = respond(inp, ne_sp, "sub", SCEN["base"], 1.0)
    w = rr["cap_pax"]
    fs = inp["fare"][rr["sel"]]
    caprows = []
    for phi in (1.0, 1.2, 1.5, 2.0, 3.0, 5.0):
        prem = (phi - 1) * fs
        vstar = np.where(rr["s"] > 0, prem / np.maximum(rr["s"], 1e-6) * 60, np.inf)
        share = 1 - norm.cdf((np.log(np.maximum(vstar, 1e-9)) - np.log(60)) / 0.75)
        share = np.where(prem <= 0, 1.0, share)
        caprows.append(dict(phi=phi, mean_fare_premium_usd=float(np.average(prem, weights=w)),
                            willing_share_vs_subsonic=float(np.average(share, weights=w))))
    CAP = pd.DataFrame(caprows)
    CAP.to_csv(HERE / "fare_premium_capture.csv", index=False, float_format="%.4g")
    print(CAP.round(3).to_string())

    # ---------------------------------------------------------------- 3g. network size: does the knee move?
    curve = R[(R.mach == 2.0) & (R.vs == "bl") & (R.capture == 1.0)
              & (R.network.str.startswith("greedy step") | (R.network == "spine"))].copy()
    curve["step"] = curve.network.map(lambda s: 0 if s == "spine" else int(s.split()[-1]))
    curve = curve.sort_values(["scenario", "step"])
    # B/C-max step, and uplift by network size (does induced demand favour bigger networks?)
    curve.to_csv(HERE / "curve_induced.csv", index=False, float_format="%.5g",
                 columns=["step", "km", "scenario", "pax_h_fixed_M", "pax_h_equiv_total_M", "bcr_fixed",
                          "bcr_plus_induced", "bcr", "uplift", "gross_induced_M"])
    kn = []
    for vname, c in curve.groupby("scenario"):
        c = c.reset_index(drop=True)
        at = lambda km: c.iloc[(c.km - km).abs().idxmin()]
        kn.append(dict(scenario=vname, bcr_max_km=c.km[c.bcr.idxmax()], bcr_spine=c.bcr[0],
                       bcr_5421=at(5421).bcr, bcr_6989=at(6989).bcr, bcr_12381=at(12381).bcr,
                       uplift_spine=c.uplift[0], uplift_5421=at(5421).uplift, uplift_6989=at(6989).uplift,
                       uplift_12381=at(12381).uplift))
    KN = pd.DataFrame(kn)
    KN.to_csv(HERE / "knee_shift.csv", index=False, float_format="%.4g")
    print(KN.round(4).to_string())

    # ---------------------------------------------------------------- 4. cost-down hook
    cd = NM.parent / "economics" / "cost-down"
    hook = sorted(cd.glob("*.csv")) if cd.exists() else []
    print("cost-down CSVs:", [h.name for h in hook])
    if (cd / "scenarios.csv").exists():      # sibling cost-reduction study: compose cost factor x demand uplift
        cdf = pd.read_csv(cd / "scenarios.csv")
        jrows = []
        for o in cdf.itertuples():
            M = 2.0 if str(o.design).startswith("M2.0") else 1.7
            for vname in ("fixed", "low", "base", "high", "base + life-avg growth x1.4"):
                row = {"design": o.design, "cost_scenario": o.scenario, "demand_scenario": vname}
                for cmp in ("bl", "sub"):
                    for cap in (0.1, 1.0):
                        q = R[(R.network == "spine") & (R.mach == M) & (R.vs == cmp) & (R.capture == cap)
                              & (R.scenario == ("base" if vname == "fixed" else vname))].iloc[0]
                        up = 1.0 if vname == "fixed" else q.bcr / q.bcr_fixed * (q.bcr_fixed / R[(R.network == "spine") & (R.mach == M) & (R.vs == cmp) & (R.capture == cap) & (R.scenario == "base")].bcr_fixed.iloc[0])
                        row[f"bcr_{cmp}_{int(cap * 100)}"] = getattr(o, f"bcr_{cmp}_{int(cap * 100)}") * up
                jrows.append(row)
        J = pd.DataFrame(jrows)
        J.to_csv(HERE / "joint_costdown.csv", index=False, float_format="%.4g")
        print(J.round(4).to_string())

    # ---------------------------------------------------------------- 5. figure
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    C1, C2, C3, C4, INK, INK2, GRID = "#2a78d6", "#eb6834", "#1baf7a", "#6250d6", "#0b0b0b", "#52514e", "#e4e3df"
    plt.rcParams.update({"font.size": 9, "axes.edgecolor": INK2, "axes.labelcolor": INK, "xtick.color": INK2,
                         "ytick.color": INK2, "font.family": "DejaVu Sans"})
    fig, (a1, a2, a3) = plt.subplots(1, 3, figsize=(13.5, 4.3), dpi=150)
    # (a) induced % vs distance, spine pairs
    pb = P[P.t_spine.notna() & (P.save_bl > 0)]
    for col, lab, cl in (("induced_bl_base", "vs Boomless, base", C1), ("induced_bl_high", "vs Boomless, high", C2),
                         ("induced_sub_base", "vs subsonic, base (fare x3)", C3)):
        a1.scatter(pb.km, pb[col] / pb.pax_yr * 100, s=np.clip(pb.pax_yr / 4e4, 2, 60), color=cl, alpha=0.6,
                   label=lab, edgecolors="none")
    a1.set_xlabel("Great-circle distance, km")
    a1.set_ylabel("Gross induced trips, % of captured trips")
    a1.set_title("(a) Induced demand per spine pair (100 % capture)", loc="left", fontsize=9.5)
    a1.legend(fontsize=7.5, frameon=False, loc="center right")
    a1.annotate("high case: same-day return newly\nfeasible (100 min access, <= 4.5 h)", (2950, 45), xytext=(3350, 52),
                fontsize=7.5, color=INK2, arrowprops=dict(arrowstyle="-", color=INK2, lw=0.8))
    # (b) B/C uplift by scenario, spine
    s2 = sp[(sp.capture == 1.0)]
    order = ["low", "base", "high", "TTB bound", "base + life-avg growth x1.4"]
    xs = np.arange(len(order))
    for k_, (cmp, cl, lab) in enumerate((("bl", C1, "vs Boomless"), ("sub", C3, "vs subsonic"))):
        vals = [s2[(s2.scenario == o) & (s2.vs == cmp)].bcr.iloc[0] for o in order]
        fixed = s2[(s2.scenario == "base") & (s2.vs == cmp)].bcr_fixed.iloc[0]
        a2.bar(xs + (k_ - 0.5) * 0.38, vals, 0.36, color=cl, label=lab)
        a2.axhline(fixed, color=cl, lw=1, ls=(0, (4, 3)))
    a2.axhline(1.0, color=INK, lw=0.8)
    a2.set_yscale("log")
    a2.set_ylim(0.02, 1.3)
    a2.text(len(order) - 0.6, 1.05, "break-even", ha="right", fontsize=7.5, color=INK)
    a2.set_xticks(xs)
    a2.set_xticklabels(["low", "base", "high", "TTB\nbound", "base +\ngrowth x1.4"], fontsize=8)
    a2.set_ylabel("B/C, spine M2.0, 100 % capture (log)")
    a2.set_title("(b) B/C with induced demand (dashed: fixed demand)", loc="left", fontsize=9.5)
    a2.legend(fontsize=7.5, frameon=False, loc="upper left", bbox_to_anchor=(0.0, 0.93))
    # (c) B/C vs network size
    for vname, cl in (("low", C4), ("base", C1), ("high", C2)):
        c = curve[curve.scenario == vname]
        a3.plot(c.km / 1e3, c.bcr, "-o", ms=3, color=cl, lw=1.8, label=f"induced, {vname}")
    c = curve[curve.scenario == "base"]
    a3.plot(c.km / 1e3, c.bcr_fixed, "--", color=INK2, lw=1.4, label="fixed demand (coverage.md)")
    a3.axvspan(5.4, 7.0, color=GRID, alpha=0.6, lw=0)
    a3.text(6.2, a3.get_ylim()[1] if False else 0.046, "knee\n5.4-7.0k km", ha="center", fontsize=7.5, color=INK2)
    a3.set_xlabel("Total shield length, 1,000 km")
    a3.set_ylabel("B/C vs Boomless, 100 % capture")
    a3.set_title("(c) Greedy growth from the spine", loc="left", fontsize=9.5)
    a3.legend(fontsize=7.5, frameon=False)
    for a in (a1, a2, a3):
        a.grid(color=GRID, lw=0.6)
        a.spines[["top", "right"]].set_visible(False)
    fig.tight_layout()
    fig.savefig(FIG / "induced_demand.png", facecolor="white")
    print("figure written")
