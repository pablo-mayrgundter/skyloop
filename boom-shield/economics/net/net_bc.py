#!/usr/bin/env python3
"""Canonical gross and net-of-fuel/CO2 benefit/cost for BoomShield (reconciles highmach, global, growth, integration).

Run:  python3 net_bc.py        (light CPU, ~5 s; reads caches, re-runs nothing heavy)
Writes next to this file:
  canonical.csv          gross / net B/C by Mach x cost case x capture x comparator, US spine (single year)
  lifecycle.csv          the same, life-cycle variant for openings 2035 / 2045 / 2050 (growth_model, separate column)
  crossing_years.csv     first opening year (2030-2100) with life-cycle B/C >= 1, gross and net
  corridors.csv          global corridors under the canonical convention (rescaled from global's corridors/greedy.csv)
  costs.csv              $/m2 and $M/km/yr per design: integrated re-run vs the ratio-scaled value
  reconciliation.csv     bridge from each study's published net figure to the canonical one
  assumptions.csv        the assumption set used here

CANONICAL CONVENTION
  * Single-year steady state at 2025 traffic and $86/h x 0.95 availability; 2025 dollars; shield cost annualised
    (7 %, 25 yr, cost_model) net of dual use; base induced demand (demand.md), capture 30 % and 100 % (linear).
  * Gross = time value + day-return credits + induced-demand consumer surplus (highmach / demand convention).
  * Net = gross - (extra fuel x $0.80/kg + extra CO2 x $190/t), CO2 3.16 kg/kg fuel, x 0.95 availability.
    Extra fuel is WHOLE-TRIP and INCREMENTAL vs the comparator: every captured trip burns the shield aircraft's
    fuel over its great-circle O&D km x 1.03 (feeders, ramps and shield km alike) minus what the comparator would
    burn on the same trip (Boomless: the domestic-sized Overture-class aircraft; subsonic: 90 g CO2/pax-km); every
    induced (net-new) trip carries its full fuel and CO2, because it does not exist in the comparator world.
  * One aircraft per Mach, used for shield sizing, ramps, fares AND fuel (highmach designs): M1.7 Overture (62 t),
    M2.0 adaptive (49 t), M2.2 design A (50 t), M3.0 design B (63 t).  The Concorde-class M2.0 row (shield for
    147 t, Overture ramps) is kept only as the legacy reference used by integration and growth.
  * Cost cases, all with two tether rows: base / plausible / aggressive (cost-down via highmach's hm_cost) and
    integrated = integration's plant + learning model RE-RUN on each design's own shield (aggressive levers,
    central), not the M2.0 integrated/aggressive ratio applied to every Mach.

Imports read-only: network-maps/highmach (hm_network, hm_cost; .cache/net_evals.pkl), network-maps/demand
(.cache/inputs.pkl, build_demand.respond), economics/integration (integrated_eval), economics/growth (evaluate,
metrics), network-maps/global (corridors.csv).  Nothing is written outside this folder.
"""
from __future__ import annotations

import copy
import importlib.util
import pickle
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.dont_write_bytecode = True
HERE = Path(__file__).resolve().parent
ECON = HERE.parent
BS = ECON.parent
NM = BS / "network-maps"
sys.path.insert(0, str(NM / "highmach"))
sys.path.insert(0, str(ECON / "integration"))
import hm_network as hn   # noqa: E402  (pulls in build_coverage, build_demand, hm_cost -> cost_down -> cost_model)
import hm_cost as hc      # noqa: E402
import integration as ig  # noqa: E402  (same cost_down / cost_model module objects as hm_cost)

_spec = importlib.util.spec_from_file_location("growth_model", ECON / "growth" / "growth_model.py")
gm = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(gm)
_energy = gm.energy


def _energy_flat(case, y):
    """'flat': canonical prices held constant ($0.80/kg, $190/t, no SAF, no efficiency gain) - isolates the
    life-cycle effect from growth's central fuel/SAF/carbon path (in-memory only; growth_model.py untouched)."""
    return (FUEL_USD_KG, SCC, 1.0, 1.0) if case == "flat" else _energy(case, y)


gm.energy = _energy_flat

bd = hn.bd
VOT, AVAIL = hn.V, hn.AVAIL                  # $86/h, 0.95
FUEL_USD_KG, CO2_KG, SCC = 0.80, 3.16, 190.0  # canonical (= highmach / global); growth's central path starts at $200/t
FUEL_SUB_PAX_KM = hn.FUEL_SUB_PAX_KM
CAPTURES = (0.3, 1.0)
LEVELS = ("base", "plausible", "aggressive", "integrated")
DESIGNS = {  # label -> highmach design name
    "M1.7 Overture (62 t)": "M1.7 Overture",
    "M2.0 adaptive (49 t)": "M2.0 adaptive",
    "M2.2 design A (50 t)": "A: M2.2 Al-limit",
    "M2.5 Ti (57 t)": "M2.5 Ti adaptive",
    "M3.0 design B (63 t)": "B: M3.0 Ti adaptive",
    "M3.5 TBCC (94 t)": "C: M3.5 TBCC",
    "M2.0 Concorde-class shield (147 t; legacy)": "repo M2.0 (55 m, 1.75 km)",
}
OPENINGS = (2035, 2045, 2050)
CROSS_YEARS = tuple(range(2030, 2101, 5))


# ------------------------------------------------------------------------------------------- costs
def integrated_cost(T_eq, width_km, key):
    """integration.integrated_eval (aggressive levers, central plants, two rows) on an arbitrary layer."""
    hc.cm.DESIGNS[key] = {"T": float(T_eq), "W": (width_km, width_km, width_km)}
    orig = ig.lever_set

    def lever_set(level):
        S = orig(level)
        S["T"] = dict(S["T"])
        S["T"][key] = float(T_eq)
        return S
    ig.lever_set = lever_set
    try:
        r, unit, rows, q = ig.integrated_eval(key, "aggressive", "central")
    finally:
        ig.lever_set = orig
    return dict(capex=r["capex_musd_km"], opex=r["opex_musd_km_yr"], dual=r["dual_musd_km_yr"],
                net=r["net_musd_km_yr"], usd_m2=r["usd_m2"])


def int_ratio_published():
    I = pd.read_csv(ECON / "integration" / "results.csv")
    m = I[I.design.str.startswith("M2.0")].set_index("scenario").net_musd_km_yr
    return m["integrated, aggressive levers, central"] / m["cost-down aggressive + 2nd tether row"]


# ------------------------------------------------------------------------------------------- benefits and fuel
def annual_flows(inp, ds, ne, cmp, cap, demand="base"):
    """Annual $ flows on the network at capture `cap` (highmach bcr_rows, split by trip type)."""
    sc = bd.SCEN[demand]
    sel = ne["sel"]
    toll = (ds["fare_ratio"] - 1.0) * sc["phi"] * inp["fare"][sel]
    r = bd.respond(inp, ne, cmp, sc, cap, toll=toll)
    dkm = inp["base"]["d"][sel]
    f_M = hn.fuel_per_pax(ds["_fuel_ac"], dkm)
    f_base = hn.fuel_per_pax(ds["_ref"], dkm) if cmp == "bl" else FUEL_SUB_PAX_KM * dkm * 1.03
    kg_cap = float((r["cap_pax"] * (f_M - f_base)).sum()) * AVAIL      # captured trips, incremental
    kg_ind = float((r["net"] * f_M).sum()) * AVAIL                      # induced trips, full fuel
    usd_kg = FUEL_USD_KG + CO2_KG / 1e3 * SCC
    ex_usd = float(r["ex_h"].sum()) * VOT
    return dict(time_busd=ex_usd * AVAIL / 1e9,
                gross_busd=(ex_usd + r["cred_usd"].sum() + r["ind_usd"].sum()) * AVAIL / 1e9,
                fuel_kg_cap_bn=kg_cap / 1e9, fuel_kg_ind_bn=kg_ind / 1e9,
                pen_cap_busd=kg_cap * usd_kg / 1e9, pen_ind_busd=kg_ind * usd_kg / 1e9,
                extra_fuel_busd=(kg_cap + kg_ind) * FUEL_USD_KG / 1e9,
                extra_co2_Mt=(kg_cap + kg_ind) * CO2_KG / 1e9,
                pax_h_M=float(r["ex_h"].sum()) / 1e6, induced_M=float(r["net"].sum()) / 1e6)


# ------------------------------------------------------------------------------------------- main
def main():
    pd.set_option("display.width", 250)
    inp = pickle.loads((NM / "demand" / ".cache" / "inputs.pkl").read_bytes())
    evals = pickle.loads((NM / "highmach" / ".cache" / "net_evals.pkl").read_bytes())
    hd = {d["design"]: d for d in hn.mach_designs()}
    nbm = pd.read_csv(NM / "highmach" / "network_by_mach.csv").set_index("design")
    ratio = int_ratio_published()

    rows, cost_tab, flows = [], {}, {}
    for lab, name in DESIGNS.items():
        ds, ne = hd[name], evals[f"{name}|spine"]
        km = ne["km"]
        ic = integrated_cost(ds["T_eq20_m"], ds["width_km"], f"NET {name}")
        costs = {lvl: ds[f"cost_{lvl}_musd_km_yr"] for lvl in ("base", "plausible", "aggressive")}
        costs["integrated"] = ic["net"]
        cost_tab[lab] = dict(ic, ratio_scaled=ds["cost_aggressive_musd_km_yr"] * ratio, **{f"cost_{k}": v for k, v in costs.items()})
        for cmp in ("bl", "sub"):
            for dem, cap in [("base", c) for c in CAPTURES] + [("high", 1.0)]:
                f = annual_flows(inp, ds, ne, cmp, cap, dem)
                if dem == "base":
                    flows[(lab, cmp, cap)] = f
                for lvl in LEVELS:
                    C = costs[lvl] * km / 1e3                  # $B/yr
                    rows.append(dict(design=lab, M=ds["M"], demand=dem, vs=cmp, capture=cap, cost=lvl,
                                     cost_musd_km_yr=round(costs[lvl], 3), cost_busd_yr=C,
                                     gross_busd=f["gross_busd"], fuel_co2_busd=f["pen_cap_busd"] + f["pen_ind_busd"],
                                     extra_co2_Mt=f["extra_co2_Mt"], pax_h_M=f["pax_h_M"],
                                     bc_gross=f["gross_busd"] / C,
                                     bc_net=(f["gross_busd"] - f["pen_cap_busd"] - f["pen_ind_busd"]) / C,
                                     bc_net_captured_only=(f["gross_busd"] - f["pen_cap_busd"]) / C))
    T = pd.DataFrame(rows)

    # ---- life-cycle variant (growth_model.evaluate / metrics; base world, central energy path, 7 %, 25 yr)
    lc, cross = [], []
    for lab, name in DESIGNS.items():
        c = cost_tab[lab]
        cost = {"integrated": dict(capex=c["capex"], opex=c["opex"], dual=c["dual"], net=c["net"])}
        flights = float(nbm.loc[name, "peak_flights_day"]) * 4716 / 4500     # base vs fixed demand (flights_vs_caps)
        for cmp in ("bl", "sub"):
            f = flows[(lab, cmp, 1.0)]
            net = dict(benefit_busd=f["gross_busd"], extra_fuel_busd=f["extra_fuel_busd"], extra_co2_Mt=f["extra_co2_Mt"],
                       km=evals[f"{name}|spine"]["km"], cost_scale=1.0, flights=flights, build=10)
            for cap in CAPTURES:
                for oy in OPENINGS:
                    m = gm.metrics(gm.evaluate(cost, net, oy, "base", "central", "integrated", cap), oy)
                    m0 = gm.metrics(gm.evaluate(cost, net, oy, "base", "flat", "integrated", cap), oy)
                    lc.append(dict(design=lab, vs=cmp, capture=cap, cost="integrated", open_year=oy,
                                   bc_gross_lifecycle=m["bc_gross"], bc_net_lifecycle=m["bc_net"],
                                   bc_net_lifecycle_flat_prices=m0["bc_net"]))
            for cap in CAPTURES:
                g, n = [], []
                for oy in CROSS_YEARS:
                    m = gm.metrics(gm.evaluate(cost, net, oy, "base", "central", "integrated", cap), oy)
                    g.append(m["bc_gross"]); n.append(m["bc_net"])
                first = lambda v: next((y for y, x in zip(CROSS_YEARS, v) if x >= 1.0), None)
                cross.append(dict(design=lab, vs=cmp, capture=cap, cost="integrated", world="base", energy="central",
                                  first_open_gross_ge1=first(g), first_open_net_ge1=first(n),
                                  max_net_to_2100=max(n)))
    L = pd.DataFrame(lc)
    X = pd.DataFrame(cross)

    # ---- top global corridors: rescale global's corridors.csv to the canonical convention
    #      (integrated cost re-run per design; induced-trip fuel added at the US-spine ratio for that design)
    C = pd.read_csv(NM / "global" / "corridors.csv")
    md = pd.read_csv(NM / "highmach" / "mach_designs.csv").set_index("design")
    gl = {"M2.2": ("M2.2 design A (50 t)", "A: M2.2 Al-limit"), "M2.0": ("M2.0 adaptive (49 t)", "M2.0 adaptive")}
    cr = []
    Gy = pd.read_csv(NM / "global" / "greedy.csv")
    Gy = Gy.assign(corridor="greedy " + Gy.scope + ": " + Gy.edges.str.slice(0, 60), design="M2.2")
    for _, q in pd.concat([C, Gy], ignore_index=True).iterrows():
        corr, dk = q.corridor, q.design
        if dk in gl:
            lab, name = gl[dk]
            c_glob = float(md.loc[name, "cost_aggressive_musd_km_yr"]) * ratio * max(q.km_overland, 1.0) / 1e3   # $B/yr
            c_can = cost_tab[lab]["net"] * max(q.km_overland, 1.0) / 1e3
            out = dict(corridor=corr, design=lab, km_overland=q.km_overland)
            for cmp in ("bl", "sub"):
                ben = q[f"bc_{cmp}_int"] * c_glob
                pen = (q[f"bc_{cmp}_int"] - q[f"bc_{cmp}_int_net"]) * c_glob
                us = flows[(lab, cmp, 1.0)]
                pen_can = pen * (1 + us["pen_ind_busd"] / us["pen_cap_busd"])
                out.update({f"global_gross_{cmp}": q[f"bc_{cmp}_int"], f"global_net_{cmp}": q[f"bc_{cmp}_int_net"],
                            f"gross_{cmp}": ben / c_can, f"net_{cmp}": (ben - pen_can) / c_can})
            cr.append(out)
    CR = pd.DataFrame(cr)

    # ---- reconciliation bridge (US spine, integrated, 100 % capture, vs Boomless)
    G = C[C.corridor == "US north spine (reference)"].set_index("design")
    sm = pd.read_csv(NM / "highmach" / "summary_by_mach.csv").set_index("design")
    hl = pd.read_csv(ECON / "growth" / "roi_headline.csv")
    br = []

    def add(study, design, step, value, note=""):
        br.append(dict(study=study, design=design, step=step, bc_net_vs_bl=round(float(value), 3), note=note))
    for dk, lab, name in (("M2.2", "M2.2 design A (50 t)", "A: M2.2 Al-limit"), ("M2.0", "M2.0 adaptive (49 t)", "M2.0 adaptive")):
        f = flows[(lab, "bl", 1.0)]
        c_sc = cost_tab[lab]["ratio_scaled"] * 3613.6 / 1e3
        c_can = cost_tab[lab]["net"] * 3613.6 / 1e3
        add("global", lab, "published (ranking.csv)", G.loc[dk, "bc_bl_int_net"],
            "global routing; integrated = aggressive x 0.552; captured trips' fuel only")
        add("global", lab, "highmach routing/demand, same cost and fuel rule", (f["gross_busd"] - f["pen_cap_busd"]) / c_sc,
            "residual = global's re-implemented routing and uplift (~1-2 %)")
        add("global", lab, "+ fuel and CO2 of induced trips", (f["gross_busd"] - f["pen_cap_busd"] - f["pen_ind_busd"]) / c_sc,
            "= highmach aggressive net / 0.552")
        add("highmach", lab, "published aggressive net / 0.552 (STATUS 'x1.8')", sm.loc[name, "bl_aggr_100_net_fuel_co2"] / ratio,
            "ratio-scaled integrated cost")
        add("canonical", lab, "integrated cost re-run on this shield", (f["gross_busd"] - f["pen_cap_busd"] - f["pen_ind_busd"]) / c_can,
            f"{cost_tab[lab]['ratio_scaled']:.2f} -> {cost_tab[lab]['net']:.2f} $M/km/yr")
    # growth: Concorde-class shield (147 t, Overture ramps) with the 49 t aircraft's fuel; life-cycle
    lab = "M2.0 Concorde-class shield (147 t; legacy)"
    f = flows[(lab, "bl", 1.0)]
    g35 = hl[(hl.design == "M2.0") & (hl.vs == "bl") & (hl.cost == "integrated") & (hl.capture == 1.0) & (hl.open_year == 2035)]
    add("growth", lab, "published life-cycle, 2035 opening (roi_headline.csv)", g35.bc_net.iloc[0], "central energy path, growth, IDC, ramp-up")
    lcl = L[(L.design == lab) & (L.vs == "bl") & (L.capture == 1.0) & (L.open_year == 2035)]
    add("growth", lab, "same, re-computed here (re-run integrated cost)", lcl.bc_net_lifecycle.iloc[0], "")
    add("growth", lab, "life-cycle 2035 at flat canonical prices", lcl.bc_net_lifecycle_flat_prices.iloc[0],
        "$0.80/kg, $190/t, no SAF or efficiency path")
    c_can = cost_tab[lab]["net"] * 3613.6 / 1e3
    add("growth", lab, "single year 2025, canonical prices", (f["gross_busd"] - f["pen_cap_busd"] - f["pen_ind_busd"]) / c_can,
        "= integration row net; shield for 147 t, fuel of 49 t")
    lab = "M2.0 adaptive (49 t)"
    f = flows[(lab, "bl", 1.0)]
    add("canonical", lab, "one 49 t aircraft for shield, ramps and fuel", (f["gross_busd"] - f["pen_cap_busd"] - f["pen_ind_busd"]) / (cost_tab[lab]["net"] * 3.6136),
        "single year 2025, canonical prices")
    BR = pd.DataFrame(br)

    A = pd.DataFrame([
        ("reference year / dollars", "2025 traffic, 2025 $; single-year steady state (life-cycle only in lifecycle.csv)"),
        ("value of time", f"${VOT:.0f}/h x {AVAIL} availability; day-return credits and induced surplus included (base demand)"),
        ("fuel price", f"${FUEL_USD_KG}/kg (jet ~$2.4/gal)"), ("CO2", f"{CO2_KG} kg/kg fuel at ${SCC:.0f}/t (EPA 2023, 2 %)"),
        ("km counted", "whole O&D trip, great circle x 1.03 (feeders, ramps, shield), both modes"),
        ("Boomless fuel", "Overture-class fixed-cycle aircraft, domestic-sized (hm_aircraft ref), design-point burn"),
        ("subsonic fuel", "90 g CO2 per pax-km"), ("seats / load factor", f"{hn.ha.SEATS} / {hn.LF}"),
        ("induced trips", "full fuel and CO2 (no such trip in the comparator world)"),
        ("aircraft per Mach", "highmach designs: M1.7 62 t, M2.0 49 t, M2.2 50 t, M3.0 63 t (shield, ramps, fares and fuel)"),
        ("cost", "two tether rows; annualised 7 %/25 yr net of dual use; integrated = integration plant model re-run per shield"),
        ("capture", "30 % and 100 % (linear)"), ("network", "US spine LAS-DEN-OMA-CHI-CLE-NYC, 3,614 km")],
        columns=["item", "canonical value"])

    T.to_csv(HERE / "canonical.csv", index=False, float_format="%.4g")
    L.to_csv(HERE / "lifecycle.csv", index=False, float_format="%.4g")
    CR.to_csv(HERE / "corridors.csv", index=False, float_format="%.4g")
    X.to_csv(HERE / "crossing_years.csv", index=False, float_format="%.4g")
    BR.to_csv(HERE / "reconciliation.csv", index=False)
    A.to_csv(HERE / "assumptions.csv", index=False)
    costs = pd.DataFrame(cost_tab).T
    costs.to_csv(HERE / "costs.csv", float_format="%.4g")

    # ---- console summary
    print(f"published integrated/aggressive ratio (M2.0 Concorde-class): {ratio:.3f}")
    print(costs[["usd_m2", "net", "ratio_scaled", "cost_aggressive"]].round(3).to_string())
    p = T[(T.capture == 1.0) & (T.demand == "base")].pivot_table(index=["design", "cost"], columns="vs", values=["bc_gross", "bc_net"])
    print(p.round(3).to_string())
    print(L[L.capture == 1.0].round(3).to_string())
    print(T[(T.demand == "high") & (T.cost == "integrated")][["design", "vs", "bc_gross", "bc_net"]].round(3).to_string())
    print(X.to_string())
    print(CR.round(3).to_string())
    print(BR.to_string())


if __name__ == "__main__":
    main()
