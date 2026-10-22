#!/usr/bin/env python3
"""Higher-Mach study: rerun the spine (and candidate networks) at Mach 1.7-4 with each concept's ramp
acceleration, shield thickness/width/altitude and cost; induced demand with demand/ conventions.

Run:  python3 hm_network.py     (one process; caches every network evaluation in .cache/ as it goes)
Writes: mach_designs.csv, network_by_mach.csv, bcr_by_mach.csv, emissions_by_mach.csv

Reads read-only: demand/.cache/inputs.pkl (DOT Table 1 + gravity fill demand, Boomless/subsonic baselines,
fares, time zones), coverage/build_coverage.py (evaluator), demand/build_demand.py (respond: induced demand,
day-return / after-work credits, rule of half), economics/cost-down (cost per km, two tether rows).
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
CACHE.mkdir(exist_ok=True)
sys.dont_write_bytecode = True
sys.path.insert(0, str(NM / "coverage"))
sys.path.insert(0, str(NM / "demand"))
import build_coverage as bc  # noqa: E402
import build_demand as bd  # noqa: E402

import hm_aircraft as ha  # noqa: E402
import hm_cost as hc  # noqa: E402
import hm_shield as hs  # noqa: E402

bt, cm = bc.bt, bc.cm
V = bd.V                      # $86/h
AVAIL = bd.AVAIL              # 0.95
LF = cm.LOAD_FACTOR           # 0.75
FUEL_USD_KG = 0.80            # jet fuel ~$2.4/gal [K]
CO2_PER_KG_FUEL = 3.16
SCC = bd.SCC                  # $190/t
ICCT_SCALE = 5.0 / 3.0        # demand.md base: Overture 5x subsonic per pax (ICCT); this model gives 3.0x -> scale for 'high'
FUEL_SUB_PAX_KM = 0.090 / CO2_PER_KG_FUEL     # kg fuel per pax-km, subsonic (demand.md 90 g CO2)


# ------------------------------------------------------------------------------------- Mach designs
def mach_designs():
    """One row per Mach: aircraft (weight, ramps, fuel), shield (d, width, altitude, gas) and cost."""
    rows = []
    ref = ha.size(1.7, "fixed", kT=1.22)                      # Overture-class, domestic-sized: fuel/fare reference
    spec = [
        ("repo M2.0 (55 m, 1.75 km)", 2.0, None),
        ("M1.7 Overture", 1.7, None),
        ("M2.0 adaptive", 2.0, "adaptive"),
        ("A: M2.2 Al-limit", 2.2, "adaptive"),
        ("M2.5 Ti adaptive", 2.5, "adaptive"),
        ("B: M3.0 Ti adaptive", 3.0, "adaptive"),
        ("C: M3.5 TBCC", 3.5, "tbcc"),
        ("D: M4.0 TBCC (+40 K gas)", 4.0, "tbcc"),
    ]
    for name, M, eng in spec:
        if name.startswith("repo"):
            ac = dict(ref, M=2.0, W_mid_t=147.0, h_cruise_km=20.25)
            a, b = bt.A_ACC, bt.A_DEC
            d, w, H, Tg, by = 55.0, 1.75, 20e3, None, "repo (Concorde 147 t)"
            fuel_ac = ha.size(2.0, "adaptive", kT=1.40)
        elif name.startswith("M1.7"):
            ac = ha.size(1.7, "fixed", kT=1.22, R_km=7870.0)
            ac["W_mid_t"] = 62.0
            a, b = 0.23, 0.46
            H, Tg = 20e3, None
            fuel_ac = ref
        else:
            kT = 1.40 if eng == "adaptive" else 1.60
            ac = ha.size(M, eng, kT=kT)
            a = ha.ramp_accel(ac)
            b = min(2 * a, 0.9)
            H = ac["h_cruise_km"] * 1e3 - hs.H_AC
            Tg = hs.atm(H)["T"] + 40.0 if M >= 3.95 else None
            fuel_ac = ac
        W_N = ac["W_mid_t"] * 1e3 * 9.81
        if not name.startswith("repo"):
            d, by = hs.d_min(W_N, M, H, Tg)
            w = hs.width_km(M, d, W_N, H, Tg)
        lift = hs.lift_per_m(H, Tg)
        T_eq = d * lift / hs.lift_per_m(20e3)
        costs = {lvl: hc.annual(T_eq, w, lvl, two_rows=True) for lvl in ("base", "plausible", "aggressive")}
        costs1 = {lvl: hc.annual(T_eq, w, lvl, two_rows=False) for lvl in ("base", "plausible", "aggressive")}
        # fare multiple relative to Overture-class: 60 % capital/crew (weight per seat / block speed) + 40 % fuel [U]
        fuel_r = ha.trip_fuel(fuel_ac, 4093) / ha.trip_fuel(ref, 4093)
        sp_r = (M / 1.7) ** 0.45            # block-speed ratio on a ~4,000 km trip (hm_time: 198 -> 158 min at M3)
        fare_r = 0.6 * (fuel_ac["t_per_seat"] / ref["t_per_seat"]) / sp_r + 0.4 * fuel_r
        if name.startswith("M1.7") or name.startswith("repo"):
            fare_r = 1.0
        g = hs.gas(Tg or hs.atm(H)["T"], 0.0, hs.atm(H)["p"])
        rows.append(dict(design=name, M=M, W_mid_t=ac["W_mid_t"], a_acc=a, a_dec=b, H_s_km=H / 1e3,
                         gas="+40 K (daytime/heated only)" if Tg else "ambient", q_kPa=0.7 * hs.atm(H + 250)["p"] * M ** 2 / 1e3,
                         ceiling_mach=g["c"] / hs.A_SOUND, d_m=d, set_by=by, width_km=w, T_eq20_m=T_eq,
                         h2_t_km=g["rho"] * d * w * 1e3, lid_kg_m2=hs.lid(d, H, Tg),
                         **{f"cost_{k}_musd_km_yr": v["net_musd_km_yr"] for k, v in costs.items()},
                         **{f"cost_{k}_onerow": v["net_musd_km_yr"] for k, v in costs1.items()},
                         capex_base_usd_m2=costs["base"]["capex_usd_m2"],
                         fuel_ratio_vs_overture=fuel_r, fare_ratio=fare_r,
                         _fuel_ac=fuel_ac, _ref=ref))
    return rows


def fuel_per_pax(ac, km):
    """kg fuel per passenger trip of great-circle km (flown x1.03) at LF 0.75, vectorised."""
    L = np.asarray(km) * 1.03
    fa = ac["f_accel"]
    f = ac["W0_t"] * 1e3 * (fa + (1 - fa) * (1 - np.exp(-L / ac["RF_km"]))) * 0.97
    return f / (ha.SEATS * LF)


# ------------------------------------------------------------------------------------- networks
def networks():
    cc = pd.read_csv(NM / "coverage" / "coverage_curve.csv")
    cc = cc[(cc.path == "greedy from spine") & (cc.mach == 2.0)].set_index("step")
    ed = lambda s: [tuple(e.split("-")) for e in s.split()]
    return {
        "spine": bc.SPINE,
        "spine+LAX": bc.SCEN["spine+LAX"],
        "knee 5,421 km": ed(cc.loc[3, "edges"]),
        "shoulder 6,989 km": ed(cc.loc[6, "edges"]),
        "south only LAX-PHX-DFW-ATL-WAS": bc.SCEN["south only LAX-PHX-DFW-ATL-WAS"],
        "trunk.md M2.5 south spine": ed("LAS-DEN DEN-MCI MCI-IND IND-DCA IND-ATL ATL-JAX"),
        "two-spine": bc.SCEN["two-spine (spine+LAX+south LAX-PHX-DFW-ATL-WAS)"],
    }


def net_eval(dfx, base, T, edges, M, sel, a, b):
    cfg = dict(bc.CFG, a_acc=a, a_dec=b)
    iu, ju = base["iu"][sel], base["ju"][sel]
    tr = bc.make_trunk(dfx, edges)
    r = bc.eval_on(tr, None, M, cfg, (iu, ju))
    t = r["t"]
    pax = T[iu, ju]
    s_bl = base["t_bl"][sel] - t
    s_bl = np.where(s_bl >= bc.CUT_MIN, s_bl, 0.0)
    s_sub = base["t_sub"][sel] - t
    s_sub = np.where(s_sub >= bc.CUT_MIN, s_sub, 0.0)
    return dict(tr=tr, sel=sel, t=t, pax=pax, s_bl=s_bl, s_sub=s_sub, on=r["on"], off=r["off"], km=tr.km,
                B_bl=float((pax * s_bl).sum() / 60), B_sub=float((pax * s_sub).sum() / 60), mach=M)


def evaluate_all(inp, designs):
    f = CACHE / "net_evals.pkl"
    out = pickle.loads(f.read_bytes()) if f.exists() else {}
    dfx = bc.ext_df(inp["df"], inp["grid"])
    base = inp["base"]
    d = base["d"]
    nets = networks()
    for ds in designs:
        M, a, b = ds["M"], ds["a_acc"], ds["a_dec"]
        ub = np.minimum(*(bt.t_trunk(0.0, d, 0.0, M, fd, 5.0, a, b)[0] for fd in ("sub", "bl")))
        sel = np.where(base["t_sub"] - np.minimum.reduce([base["t_ub_M1.7"], base["t_ub_M2.0"], ub]) >= bc.CUT_MIN)[0]
        names = list(nets) if not ds["design"].startswith("repo") else ["spine"]
        for nn in names:
            key = f"{ds['design']}|{nn}"
            if key in out:
                continue
            ne = net_eval(dfx, base, inp["T"], nets[nn], M, sel, a, b)
            ne["tr"] = None                                   # keep the cache small; trunk rebuilt if needed
            out[key] = dict(ne, name=nn, design=ds["design"])
            f.write_bytes(pickle.dumps(out))
            print(f"  {key}: {ne['km']:.0f} km, {ne['B_bl'] / 1e6:.1f} M pax-h vs BL, {ne['B_sub'] / 1e6:.1f} vs sub",
                  flush=True)
    return out


def peak_load(inp, ne, edges, s_key="s_bl"):
    dfx = bc.ext_df(inp["df"], inp["grid"])
    tr = bc.make_trunk(dfx, edges)
    ev = dict(pax=ne["pax"], sav=ne[s_key], on=ne["on"], off=ne["off"], use=ne[s_key] > 0)
    lb = bt.segment_loads(tr, ev, "sav", "use")
    return max([v[0] for v in lb.values()] or [0])


# ------------------------------------------------------------------------------------- economics
def bcr_rows(inp, ds, ne, scen_name):
    sc = bd.SCEN[scen_name]
    fs = inp["fare"][ne["sel"]]
    toll = (ds["fare_ratio"] - 1.0) * sc["phi"] * fs           # fares tied to fuel burn / weight per seat
    rows = []
    dkm = inp["base"]["d"][ne["sel"]]
    f_M = fuel_per_pax(ds["_fuel_ac"], dkm)
    f_ref = fuel_per_pax(ds["_ref"], dkm)
    f_sub = FUEL_SUB_PAX_KM * dkm * 1.03
    for cmp in ("bl", "sub"):
        for cap in (0.1, 1.0):
            r = bd.respond(inp, ne, cmp, sc, cap, toll=toll)
            ex_usd = r["ex_h"].sum() * V
            ben_time = ex_usd * AVAIL
            ben = (ex_usd + r["cred_usd"].sum() + r["ind_usd"].sum()) * AVAIL
            trips = r["cap_pax"] + r["net"]                       # captured + net-new trips on the shield aircraft
            f_base = f_ref if cmp == "bl" else f_sub
            d_fuel = ((r["cap_pax"] * (f_M - f_base)) + r["net"] * f_M).sum()     # kg/yr extra vs reference world
            fuel_usd = d_fuel * FUEL_USD_KG * AVAIL
            co2_t = d_fuel * CO2_PER_KG_FUEL / 1e3 * AVAIL
            for lvl in ("base", "plausible", "aggressive"):
                cost = ds[f"cost_{lvl}_musd_km_yr"] * 1e6 * ne["km"]
                cost1 = ds[f"cost_{lvl}_onerow"] * 1e6 * ne["km"]
                rows.append(dict(design=ds["design"], M=ds["M"], network=ne["name"], km=ne["km"], demand=scen_name,
                                 vs=cmp, capture=cap, cost=lvl, captured_pax_M=r["cap_pax"].sum() / 1e6,
                                 net_induced_M=r["net"].sum() / 1e6, pax_h_M=r["ex_h"].sum() / 1e6,
                                 day_return_pax_M=(r["cap_pax"] * r["day"]).sum() / 1e6,
                                 benefit_busd=ben / 1e9, extra_fuel_busd=fuel_usd / 1e9, extra_co2_Mt=co2_t / 1e6,
                                 carbon_busd=co2_t * SCC / 1e9, cost_busd=cost / 1e9,
                                 bcr_time_onerow=ben_time / cost1, bcr_time=ben_time / cost, bcr=ben / cost,
                                 bcr_net_fuel=(ben - fuel_usd) / cost,
                                 bcr_net_fuel_carbon=(ben - fuel_usd - co2_t * SCC) / cost,
                                 bcr_net_fuel_carbon_icct=(ben - ICCT_SCALE * (fuel_usd + co2_t * SCC)) / cost))
    return rows


if __name__ == "__main__":
    pd.set_option("display.width", 260)
    inp = pickle.loads((NM / "demand" / ".cache" / "inputs.pkl").read_bytes())
    designs = mach_designs()
    D = pd.DataFrame([{k: v for k, v in r.items() if not k.startswith("_")} for r in designs])
    D.to_csv(HERE / "mach_designs.csv", index=False, float_format="%.4g")
    print(D.round(3).to_string())
    evals = evaluate_all(inp, designs)
    nets = networks()
    rows, nrows = [], []
    for ds in designs:
        for nn in (nets if not ds["design"].startswith("repo") else ["spine"]):
            ne = dict(evals[f"{ds['design']}|{nn}"])
            for scen in ("base", "high"):
                rows += bcr_rows(inp, ds, ne, scen)
    R = pd.DataFrame(rows)
    R.to_csv(HERE / "bcr_by_mach.csv", index=False, float_format="%.5g")
    # network choice per Mach: best B/C vs Boomless, base demand, 100 % capture, base cost
    sel = R[(R.demand == "base") & (R.vs == "bl") & (R.capture == 1.0) & (R.cost == "base")]
    print(sel[["design", "network", "km", "captured_pax_M", "pax_h_M", "net_induced_M", "day_return_pax_M", "bcr_time",
               "bcr", "bcr_net_fuel", "bcr_net_fuel_carbon"]].round(4).to_string())
    for ds in designs:
        s = sel[sel.design == ds["design"]]
        best = s.loc[s.bcr.idxmax()]
        ne = evals[f"{ds['design']}|{best.network}"]
        pk = peak_load(inp, ne, nets[best.network])
        nrows.append(dict(design=ds["design"], best_network=best.network, km=best.km, bcr=best.bcr,
                          spine_bcr=float(s[s.network == "spine"].bcr.iloc[0]),
                          peak_seg_pax_M=pk / 1e6, peak_flights_day=pk / (ha.SEATS * LF * 365)))
    N = pd.DataFrame(nrows)
    N.to_csv(HERE / "network_by_mach.csv", index=False, float_format="%.4g")
    print(N.round(4).to_string())
