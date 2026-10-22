#!/usr/bin/env python3
"""BoomShield unit-cost reduction study (screening level).

Run:  python3 cost_down.py
Writes (next to this file): levers.csv, scenarios.csv, waterfall.csv, pv_dualuse.csv, learning.csv,
scale_check.csv, bcr_new.csv, cost_down.png.  Needs physics.py outputs only for the write-up.

Imports ../cost_model.py READ-ONLY (sys.dont_write_bytecode; parameters are overridden in memory inside
a context manager and restored).  The base case reproduces cost_model's base annual cost per km exactly.
Spine benefits come from ../../network-maps/coverage/scenarios.csv (no fatigue cap, per lid-fatigue.md).
"""
from __future__ import annotations

import contextlib
import copy
import csv
import math
import sys
from pathlib import Path

sys.dont_write_bytecode = True
HERE = Path(__file__).resolve().parent
ECON = HERE.parent
BS = ECON.parent
sys.path.insert(0, str(ECON))
import cost_model as cm  # noqa: E402

cm.load_operations()
M20 = "M2.0 (Concorde-class, psf baseline)"
M17 = "M1.7 (Overture-class, lift-set)"
R, N = 0.07, 25.0          # finance held fixed across scenarios (cost_model base)
CRF = cm.crf(R, N)
COVERAGE = BS / "network-maps" / "coverage" / "scenarios.csv"
DEMAND_DIR = BS / "network-maps" / "demand"

# ----------------------------------------------------------------------------------------------
# Scenario parameter sets.  Every value is a single number; [U] = unverified engineering estimate.
# ----------------------------------------------------------------------------------------------
BASE = {
    # geometry / mass
    "T": {M20: 55.0, M17: 36.0}, "holddown": 0.55, "struct_kg": None,   # None = whole lid at film price
    "film_usd_kg": 4.0, "ballast_usd_kg": 4.0,
    # envelope unit costs
    "skin_usd_m2": 2.0, "uhmwpe_usd_kg": 30.0, "floor_usd_kg": 15.0, "partition_usd_m2": 1.0,
    "ties_usd_m2": 0.5, "seams_usd_m2": 6.0, "shuttle_usd_m2": 1.0, "h2_usd_kg": 4.0, "h2_deliv_usd_kg": 1.5,
    # field
    "tether_f": 1.0, "spacing_km": 3.0, "station_musd": 25.0, "controls_usd_cell": 150e3, "install_usd_m2": 6.0,
    "eng_frac": 0.20, "cont_frac": 0.50,
    # opex
    "life_bulk": 10.0, "life_skin": 10.0, "life_floor": 10.0, "skin_robot_usd_m2": 0.0,
    "h2_makeup": 0.10, "tether_life": 15.0, "om_frac": 0.03, "inspection_musd": 0.15,
    "lid_power_kwh_m2": 1.0, "station_mw": 0.5, "fte_site": 8.0, "central_musd": 0.05, "insurance": 0.01,
    # dual use
    "pv_frac": 0.0, "pv": None, "hosting_musd_km": 0.0, "two_rows": False,
}

PV_PLAUS = {"H": 2680.0, "eta": 0.15, "pr": 0.80, "price": 35.0, "module_usd_w": 0.35, "bos_usd_w": 0.15,
            "tether_cond_usd_w": 0.05, "grid_usd_w": 0.15, "om_usd_kw_yr": 8.0, "life": 20.0, "mass_kg_m2": 1.0,
            "superheat_K": 20.0, "ground_lcoe": 38.0}
PV_AGGR = {"H": 2680.0, "eta": 0.18, "pr": 0.82, "price": 40.0, "module_usd_w": 0.20, "bos_usd_w": 0.10,
           "tether_cond_usd_w": 0.03, "grid_usd_w": 0.12, "om_usd_kw_yr": 6.0, "life": 25.0, "mass_kg_m2": 0.8,
           "superheat_K": 15.0, "ground_lcoe": 28.0}

LEVERS = [
    # (key, label, plausible overrides, aggressive overrides)
    ("mass", "1. Mass: cheap ballast ply; M1.7 re-sized to psf (21 m)",
     {"T": {M20: 55.0, M17: 21.0}, "holddown": 0.62, "struct_kg": 0.35, "ballast_usd_kg": 1.2},
     {"T": {M20: 54.6, M17: 20.4}, "holddown": 0.25, "struct_kg": 0.30, "ballast_usd_kg": 0.8}),
    ("materials", "2. Materials & roll-to-roll manufacturing",
     {"film_usd_kg": 2.5, "skin_usd_m2": 1.5, "uhmwpe_usd_kg": 22.0, "floor_usd_kg": 4.0, "partition_usd_m2": 0.6,
      "ties_usd_m2": 0.35, "seams_usd_m2": 3.0, "shuttle_usd_m2": 0.8, "controls_usd_cell": 75e3,
      "install_usd_m2": 3.5, "h2_usd_kg": 3.0, "h2_deliv_usd_kg": 1.0},
     {"film_usd_kg": 1.8, "skin_usd_m2": 1.0, "uhmwpe_usd_kg": 15.0, "floor_usd_kg": 2.5, "partition_usd_m2": 0.4,
      "ties_usd_m2": 0.25, "seams_usd_m2": 1.5, "shuttle_usd_m2": 0.6, "controls_usd_cell": 40e3,
      "install_usd_m2": 2.0, "h2_usd_kg": 2.0, "h2_deliv_usd_kg": 0.5}),
    ("noak", "3. Nth-of-a-kind engineering & contingency",
     {"eng_frac": 0.12, "cont_frac": 0.30}, {"eng_frac": 0.08, "cont_frac": 0.20}),
    ("life", "4. Life: replaceable UV skin, robots",
     {"life_bulk": 20.0, "life_skin": 8.0, "life_floor": 20.0, "skin_robot_usd_m2": 0.5, "h2_makeup": 0.06,
      "tether_life": 20.0},
     {"life_bulk": 30.0, "life_skin": 10.0, "life_floor": 30.0, "skin_robot_usd_m2": 0.3, "h2_makeup": 0.03,
      "tether_life": 25.0}),
    ("ops", "5. Operations: unstaffed stations, insurance, rim rounding",
     {"fte_site": 1.0, "insurance": 0.006, "tether_f": 0.60, "spacing_km": 4.0, "station_musd": 15.0,
      "inspection_musd": 0.08, "om_frac": 0.025, "lid_power_kwh_m2": 0.0, "station_mw": 0.3},
     {"fte_site": 0.3, "insurance": 0.0035, "tether_f": 0.40, "spacing_km": 5.0, "station_musd": 10.0,
      "inspection_musd": 0.05, "om_frac": 0.02, "lid_power_kwh_m2": 0.0, "station_mw": 0.2}),
    ("dual", "6. Dual-use: PV (market-capped) + hosting",
     {"pv_frac": 0.05, "pv": PV_PLAUS, "hosting_musd_km": 0.012},
     {"pv_frac": 0.15, "pv": PV_AGGR, "hosting_musd_km": 0.055}),
]


@contextlib.contextmanager
def overrides(S, design):
    """Temporarily set cost_model module state from scenario S (restored on exit)."""
    saved_P = copy.deepcopy(cm.P)
    saved_hd = cm.HOLDDOWN
    saved_T = cm.DESIGNS[design]["T"]
    try:
        def put(k, val):
            cm.P[k] = ((val, val, val),) + tuple(cm.P[k][1:])
        put("lid_pe_usd_per_kg", S["film_usd_kg"])
        put("lid_skin_usd_m2", S["skin_usd_m2"])
        put("uhmwpe_usd_per_kg", S["uhmwpe_usd_kg"])
        put("floor_film_usd_per_kg", S["floor_usd_kg"])
        put("partition_usd_m2_film", S["partition_usd_m2"])
        put("ties_usd_m2", S["ties_usd_m2"])
        put("seams_fab_usd_m2", S["seams_usd_m2"])
        put("shuttle_usd_m2", S["shuttle_usd_m2"])
        put("h2_usd_per_kg", S["h2_usd_kg"])
        put("h2_delivery_usd_per_kg", S["h2_deliv_usd_kg"])
        put("spacing_km", S["spacing_km"])
        put("station_musd", S["station_musd"])
        put("controls_usd_per_cell", S["controls_usd_cell"])
        put("install_usd_m2", S["install_usd_m2"])
        cm.HOLDDOWN = S["holddown"]
        cm.DESIGNS[design]["T"] = S["T"][design]
        yield
    finally:
        cm.P.clear()
        cm.P.update(saved_P)
        cm.HOLDDOWN = saved_hd
        cm.DESIGNS[design]["T"] = saved_T


def pv_economics(pv, tether_usd_m2):
    """Per m2 of PV-covered lid: revenue, annualised cost and surplus ($/m2/yr)."""
    wp = pv["eta"] * 1000.0                                 # W_p per m2
    e_kwh = pv["H"] * pv["eta"] * pv["pr"]                  # kWh/m2/yr delivered at the ground converter
    rev = e_kwh * pv["price"] / 1000.0
    # competitive view: long-run capture price cannot exceed the LCOE of ground PV of the same era,
    # so the rent the shield can keep is (min(price, ground LCOE) - shield LCOE) x energy
    rev_comp = e_kwh * min(pv["price"], pv["ground_lcoe"]) / 1000.0
    capex_w = pv["module_usd_w"] + pv["bos_usd_w"] + pv["tether_cond_usd_w"] + pv["grid_usd_w"]
    crf_pv = cm.crf(R, pv["life"])
    ann_pv = capex_w * wp * crf_pv + pv["om_usd_kw_yr"] * wp / 1000.0
    # thermal: a dark PV lid adds superheat -> daily lift swing ~0.46 %/K of lift; tethers absorb it.
    # tether mass ~ total tension; +20 K adds ~0.4 kg/m2 of swing on a 4.5 kg/m2 lift layer, ~+55 % of hold-down.
    thermal = tether_usd_m2 * 0.55 * (pv["superheat_K"] / 20.0) * CRF
    ballast_credit = pv["mass_kg_m2"] * 1.0 * CRF            # PV mass replaces ~$1/kg ballast ply
    margin_price_taker = rev - ann_pv - thermal + ballast_credit
    surplus = rev_comp - ann_pv - thermal + ballast_credit
    return {"Wp_m2": wp, "kWh_m2_yr": e_kwh, "revenue": rev, "capex_usd_w": capex_w, "capex_usd_m2": capex_w * wp,
            "annualised_pv": ann_pv, "lcoe_usd_mwh": ann_pv / e_kwh * 1000.0, "thermal_penalty": thermal,
            "ballast_credit": ballast_credit, "margin_price_taker": margin_price_taker, "surplus": surplus}


def evaluate(design, S):
    with overrides(S, design):
        c = cm.capex(design, "base")
        lines = dict(c["lines_usd_m2"])
        A = c["A_m2_per_km"]
        for k in ("Engineering / PM", "Contingency", "Land and easements"):
            lines.pop(k)
        # lid: split structural ply (film price) from ballast ply (filled compound) - mass unchanged
        lid_bulk = max(c["lid_kg_m2"] - cm.v("tape_kg_m2", "base") - 0.05, 0)
        if S["struct_kg"] is not None:
            s_kg = min(S["struct_kg"], lid_bulk)
            b_kg = lid_bulk - s_kg
            if S["pv_frac"] > 0 and S["pv"]:
                b_kg = max(b_kg - S["pv_frac"] * S["pv"]["mass_kg_m2"], 0)   # PV replaces ballast (credited in PV)
            lines["Lid mass (PE sheet/laminate)"] = s_kg * S["film_usd_kg"] + b_kg * S["ballast_usd_kg"]
        lines["Tethers"] *= S["tether_f"]
        lines["Anchors"] *= (0.5 + 0.5 * S["tether_f"])   # block scales with load; hardware/site works do not
        if S.get("two_rows"):   # operations: moorings for both crosswind directions double the count
            for k in ("Tethers", "Anchors", "Winch/tension stations"):
                lines[k] *= 2.0
        direct = sum(lines.values())
        eng = direct * S["eng_frac"]
        cont = (direct + eng) * S["cont_frac"]
        land = cm.land_rights("base")["bs_cost_usd_km"] / A
        usd_m2 = direct + eng + cont + land
        # ---- opex (own form of cost_model.opex with split lives) ----
        s = S["spacing_km"]
        h2_price = S["h2_usd_kg"] + S["h2_deliv_usd_kg"]
        o = {}
        o["H2 make-up"] = c["h2_t_km"] * 1e3 * S["h2_makeup"] * h2_price
        bulk_keys = ["Lid mass (PE sheet/laminate)", "Load tapes / scrim (UHMWPE)", "Partitions and bellows",
                     "Soft ties", "Seams and fabrication", "Launch / installation"]
        bulk = sum(lines[k] for k in bulk_keys)
        skin = lines["Lid UV/weathering skin"]
        floor = lines["Floor film"]
        o["Envelope replacement"] = A * (bulk / S["life_bulk"] + (skin + S["skin_robot_usd_m2"]) / S["life_skin"]
                                         + floor / S["life_floor"])
        o["Tether replacement"] = lines["Tethers"] * A / S["tether_life"]
        hw = sum(lines[k] for k in ["Anchors", "Winch/tension stations", "Controls, sensors, transponders",
                                    "Diurnal gas shuttle"])
        o["O&M stations/controls"] = hw * A * S["om_frac"]
        o["Inspection"] = S["inspection_musd"] * 1e6
        o["Power"] = (S["lid_power_kwh_m2"] * A + S["station_mw"] * 8760e3 / s) * cm.v("power_usd_kwh", "base")
        nsite = 2.0 if S.get("two_rows") else 1.0
        o["Easement rent"] = nsite * cm.v("easement_rent_usd_site_yr", "base") / s
        o["Staffing"] = nsite * S["fte_site"] * cm.v("fte_usd", "base") / s + S["central_musd"] * 1e6
        o["Insurance"] = usd_m2 * A * S["insurance"]
        opex = sum(o.values())
        annual = usd_m2 * A * CRF + opex
        # ---- dual use ----
        dual = S["hosting_musd_km"] * 1e6
        pvd = None
        if S["pv_frac"] > 0 and S["pv"]:
            pvd = pv_economics(S["pv"], lines["Tethers"])
            dual += max(pvd["surplus"], 0.0) * S["pv_frac"] * A
        return {"design": design, "W_km": c["W_km"], "T": c["T"], "lid_kg_m2": c["lid_kg_m2"], "A": A,
                "lines": lines, "eng": eng, "cont": cont, "land": land, "usd_m2": usd_m2,
                "capex_musd_km": usd_m2 * A / 1e6, "opex": o, "opex_musd_km_yr": opex / 1e6,
                "annual_musd_km_yr": annual / 1e6, "dual_musd_km_yr": dual / 1e6,
                "net_musd_km_yr": (annual - dual) / 1e6, "net_usd_m2_yr": (annual - dual) / A, "pv": pvd,
                "h2_t_km": c["h2_t_km"], "tether_t_km": c["tether_t_km"] * S["tether_f"]}


def scenario(level):
    """level: 'base', 'plausible', 'aggressive'. Returns cumulative list [(label, S)]."""
    S = copy.deepcopy(BASE)
    steps = [("Base (cost_model)", copy.deepcopy(S))]
    if level == "base":
        return steps
    for key, label, pl, ag in LEVERS:
        S.update(copy.deepcopy(pl if level == "plausible" else ag))
        steps.append((label, copy.deepcopy(S)))
    return steps


def single_lever(level, key):
    S = copy.deepcopy(BASE)
    for k, label, pl, ag in LEVERS:
        if k == key:
            S.update(copy.deepcopy(pl if level == "plausible" else ag))
            if key == "dual":            # PV needs the ballast split to credit mass; keep base otherwise
                pass
    return S


# ----------------------------------------------------------------------------------------------
def learning_table():
    """Wright's law: cumulative-average cost over the spine build relative to the pilot unit cost."""
    modules = 3614 * 1.75e6 / (250 * 1.75e3)   # 250 m x full-width modules on the spine
    rows = []
    for lr in (0.10, 0.15, 0.20):
        b = -math.log2(1 - lr)
        for x0_frac in (0.01, 0.03):
            x0 = modules * x0_frac
            Nn = modules
            avg = (Nn ** (1 - b) - x0 ** (1 - b)) / ((1 - b) * (Nn - x0)) / x0 ** (-b)
            last = (Nn / x0) ** (-b)
            rows.append({"learning_rate": lr, "pilot_share_of_spine": x0_frac, "modules": round(modules),
                         "doublings_after_pilot": round(math.log2(Nn / x0), 1),
                         "cum_avg_vs_pilot": round(avg, 3), "last_unit_vs_pilot": round(last, 3)})
    return rows


def scale_rows():
    area_km2 = 3614 * 1.75
    out = []
    for label, kgm2 in (("M2.0 lid 3.9 kg/m2", 3.9), ("M2.0 polymer if 60 % filled ballast", 0.35 + 0.4 * 3.4),
                        ("M1.7 re-sized lid 1.0 kg/m2", 1.0), ("M1.7 base lid 2.3 kg/m2", 2.33)):
        mt = area_km2 * 1e6 * kgm2 / 1e9
        out.append({"case": label, "spine_area_km2": round(area_km2), "Mt_total": round(mt, 1),
                    "Mt_per_yr_10yr_build": round(mt / 10, 2),
                    "share_world_PE_110Mt": f"{mt / 10 / 110:.1%}",
                    "share_world_ag_film_6Mt": f"{mt / 10 / 6:.0%}",
                    "area_vs_China_mulch_200000km2_yr": f"{area_km2 / 10 / 200000:.2%}"})
    return out


def coverage_refs():
    ref = {}
    with open(COVERAGE, newline="") as f:
        for r in csv.DictReader(f):
            if r["scenario"] == "spine":
                m = float(r["mach"])
                ref[m] = {k: float(r[k]) for k in ("cost_musd_km_yr", "bcr_bl_cap100_nofatigue", "bcr_sub_cap100_nofatigue",
                                                    "bcr_bl_cap10_nofatigue", "bcr_sub_cap10_nofatigue", "trunk_km")}
    return ref


def joint_with_demand(out_rows):
    """Combine cost scenarios with the sibling induced-demand study (network-maps/demand/scenario_results.csv):
    its 'uplift' multiplies spine B/C at fixed cost; benefits scale linearly, so the factors compose."""
    f = DEMAND_DIR / "scenario_results.csv"
    if not f.exists():
        return []
    up, fixed = {}, {}
    with open(f, newline="") as fh:
        rows_d = [r for r in csv.DictReader(fh) if r["network"] == "spine" and r["capture"] == "1"]
    for r in rows_d:   # reference: fixed-demand B/C of the 'base' row (equals coverage.md's spine)
        if r["scenario"] == "base":
            fixed[(float(r["mach"]), r["vs"])] = float(r["bcr_fixed"])
    for r in rows_d:
        k = (float(r["mach"]), r["vs"])
        up[(k[0], r["scenario"], r["vs"])] = float(r["bcr"]) / fixed[k]
    rows = []
    for o in out_rows:
        mach = 2.0 if o["design"].startswith("M2.0") else 1.7
        for dsc in ("base", "high", "base + life-avg growth x1.4"):
            if (mach, dsc, "bl") not in up:
                continue
            ub, us = up[(mach, dsc, "bl")], up[(mach, dsc, "sub")]
            rows.append({"design": o["design"], "cost_scenario": o["scenario"], "demand_scenario": dsc,
                         "uplift_vs_bl": round(ub, 3), "uplift_vs_sub": round(us, 3),
                         "bcr_bl_10": round(o["bcr_bl_10"] * ub, 4), "bcr_bl_100": round(o["bcr_bl_100"] * ub, 3),
                         "bcr_sub_10": round(o["bcr_sub_10"] * us, 4), "bcr_sub_100": round(o["bcr_sub_100"] * us, 3)})
    return rows


def demand_hook():
    """Induced-demand sibling study. Returns (multiplier on benefits, note)."""
    if not DEMAND_DIR.exists():
        return None, "network-maps/demand/ absent: no joint scenario (hook: multiply benefits by its uplift)"
    files = sorted(p.name for p in DEMAND_DIR.glob("*"))
    return None, f"network-maps/demand/ present ({', '.join(files[:8])}); read manually for the joint scenario"


def main():
    ref = coverage_refs()
    out_rows, wf_rows, lever_rows, pv_rows = [], [], [], []
    results = {}
    for design, mach in ((M20, 2.0), (M17, 1.7)):
        for level in ("plausible", "aggressive"):
            prev = None
            for label, S in scenario(level):
                r = evaluate(design, S)
                if prev is None and level == "plausible":
                    results[(design, "base")] = r
                wf_rows.append({"design": design, "level": level, "step": label, "usd_m2": round(r["usd_m2"], 2),
                                "capex_musd_km": round(r["capex_musd_km"], 1), "opex_musd_km_yr": round(r["opex_musd_km_yr"], 2),
                                "annual_musd_km_yr": round(r["annual_musd_km_yr"], 2),
                                "dual_musd_km_yr": round(r["dual_musd_km_yr"], 3),
                                "net_musd_km_yr": round(r["net_musd_km_yr"], 2),
                                "net_usd_m2_yr": round(r["net_usd_m2_yr"], 2),
                                "step_saving_musd_km_yr": round(prev["net_musd_km_yr"] - r["net_musd_km_yr"], 2) if prev else "",
                                "T_m": r["T"], "lid_kg_m2": round(r["lid_kg_m2"], 2)})
                prev = r
            results[(design, level)] = prev
            if prev["pv"]:
                pv_rows.append({"design": design, "level": level, **{k: round(v, 3) for k, v in prev["pv"].items()},
                                "pv_frac": [s for s in scenario(level)][-1][1]["pv_frac"]})
        # single levers applied alone to base (for the lever table)
        b = results[(design, "base")]
        for key, label, pl, ag in LEVERS:
            row = {"design": design, "lever": label}
            for level, ov in (("plausible", pl), ("aggressive", ag)):
                S = copy.deepcopy(BASE)
                S.update(copy.deepcopy(ov))
                if key == "dual" and S["struct_kg"] is None:
                    S["struct_kg"] = None
                r = evaluate(design, S)
                row[f"{level}_usd_m2"] = round(r["usd_m2"], 1)
                row[f"{level}_net_musd_km_yr"] = round(r["net_musd_km_yr"], 2)
                row[f"{level}_factor"] = round(b["net_musd_km_yr"] / r["net_musd_km_yr"], 2)
            lever_rows.append(row)

    # full-lid PV upper bound (ignores market saturation) for the write-up
    for name, pv in (("plausible", PV_PLAUS), ("aggressive", PV_AGGR), ("low price $25", {**PV_PLAUS, "price": 25.0})):
        d = pv_economics(pv, results[(M20, "base")]["lines"]["Tethers"])
        pv_rows.append({"design": "per m2 of PV (any design)", "level": name, **{k: round(v, 3) for k, v in d.items()},
                        "pv_frac": 1.0})

    # scenarios + B/C
    for design, mach in ((M20, 2.0), (M17, 1.7)):
        rf = ref[mach]
        for level in ("base", "plausible", "aggressive"):
            r = results[(design, level)]
            scale = rf["cost_musd_km_yr"] / r["net_musd_km_yr"]
            out_rows.append({"design": design, "scenario": level, "T_m": r["T"], "lid_kg_m2": round(r["lid_kg_m2"], 2),
                             "capex_usd_m2": round(r["usd_m2"], 1), "capex_musd_km": round(r["capex_musd_km"], 1),
                             "opex_musd_km_yr": round(r["opex_musd_km_yr"], 2),
                             "annual_musd_km_yr": round(r["annual_musd_km_yr"], 2),
                             "dual_use_musd_km_yr": round(r["dual_musd_km_yr"], 3),
                             "net_musd_km_yr": round(r["net_musd_km_yr"], 2), "net_usd_m2_yr": round(r["net_usd_m2_yr"], 2),
                             "factor_vs_base": round(rf["cost_musd_km_yr"] / r["net_musd_km_yr"], 2),
                             "bcr_bl_10": round(rf["bcr_bl_cap10_nofatigue"] * scale, 4),
                             "bcr_bl_100": round(rf["bcr_bl_cap100_nofatigue"] * scale, 3),
                             "bcr_sub_10": round(rf["bcr_sub_cap10_nofatigue"] * scale, 4),
                             "bcr_sub_100": round(rf["bcr_sub_cap100_nofatigue"] * scale, 3),
                             "benefit_bl_100_musd_km_yr": round(rf["bcr_bl_cap100_nofatigue"] * rf["cost_musd_km_yr"], 2),
                             "benefit_sub_100_musd_km_yr": round(rf["bcr_sub_cap100_nofatigue"] * rf["cost_musd_km_yr"], 2)})

    # sensitivities on the scenarios: second tether row; aggressive without dual use; PV price-taker
    sens = []
    for design, mach in ((M20, 2.0), (M17, 1.7)):
        rf = ref[mach]
        for level in ("base", "plausible", "aggressive"):
            S = copy.deepcopy(scenario(level)[-1][1])
            variants = [("two tether rows (both crosswind directions)", {"two_rows": True})]
            if level == "aggressive":
                variants.append(("no dual use", {"pv_frac": 0.0, "hosting_musd_km": 0.0}))
            for vlabel, ov in variants:
                S2 = copy.deepcopy(S)
                S2.update(ov)
                r = evaluate(design, S2)
                sc = rf["cost_musd_km_yr"] / r["net_musd_km_yr"]
                sens.append({"design": design, "scenario": level, "variant": vlabel, "capex_usd_m2": round(r["usd_m2"], 1),
                             "net_musd_km_yr": round(r["net_musd_km_yr"], 2),
                             "bcr_bl_100": round(rf["bcr_bl_cap100_nofatigue"] * sc, 3),
                             "bcr_sub_100": round(rf["bcr_sub_cap100_nofatigue"] * sc, 3)})
            if level == "aggressive":
                r = results[(design, level)]
                pt = r["pv"]["margin_price_taker"] * S["pv_frac"] * r["A"] / 1e6 + S["hosting_musd_km"]
                net = r["annual_musd_km_yr"] - pt
                sc = rf["cost_musd_km_yr"] / net
                sens.append({"design": design, "scenario": level, "variant": "PV as price-taker at $40/MWh (no ground-PV cap)",
                             "capex_usd_m2": round(r["usd_m2"], 1), "net_musd_km_yr": round(net, 2),
                             "bcr_bl_100": round(rf["bcr_bl_cap100_nofatigue"] * sc, 3),
                             "bcr_sub_100": round(rf["bcr_sub_cap100_nofatigue"] * sc, 3)})
    for r in sens:
        print("SENS", r)

    def dump(name, rows):
        with open(HERE / name, "w", newline="") as f:
            w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
            w.writeheader()
            w.writerows(rows)
    dump("scenarios.csv", out_rows)
    dump("waterfall.csv", wf_rows)
    dump("levers.csv", lever_rows)
    dump("pv_dualuse.csv", pv_rows)
    dump("learning.csv", learning_table())
    dump("scale_check.csv", scale_rows())
    dump("sensitivities.csv", sens)
    joint = joint_with_demand(out_rows)
    if joint:
        dump("joint_demand.csv", joint)
        for r in joint:
            print("JOINT", r)
    cl = []
    for design in (M20, M17):
        for level in ("base", "plausible", "aggressive"):
            r = results[(design, level)]
            L = r["lines"]
            cats = {"Lid & envelope": sum(L[k] for k in ["Lid mass (PE sheet/laminate)", "Lid UV/weathering skin",
                                                          "Load tapes / scrim (UHMWPE)", "Floor film", "Partitions and bellows",
                                                          "Soft ties", "Seams and fabrication"]),
                    "Tethers, anchors, stations": L["Tethers"] + L["Anchors"] + L["Winch/tension stations"],
                    "H2, controls, shuttle, install": L["H2 first fill"] + L["Controls, sensors, transponders"]
                    + L["Diurnal gas shuttle"] + L["Launch / installation"],
                    "Engineering, contingency, land": r["eng"] + r["cont"] + r["land"]}
            for k, val in cats.items():
                cl.append({"design": design, "scenario": level, "category": k, "usd_m2": round(val, 2)})
    dump("capex_categories.csv", cl)

    # base reproduction check
    for design in (M20, M17):
        cx = cm.capex(design, "base")
        a = cm.annual_cost_km(cx, cm.opex(cx, "base"), "base") / 1e6
        print(f"check {design}: cost_model {a:.3f}  this {results[(design, 'base')]['annual_musd_km_yr']:.3f} M$/km/yr")
    for r in out_rows:
        print(r)
    for r in lever_rows:
        print(r)
    for r in pv_rows:
        print(r)
    for r in learning_table():
        print(r)
    for r in scale_rows():
        print(r)
    print(demand_hook()[1])
    # line detail for the write-up
    for key in [(M20, "base"), (M20, "plausible"), (M20, "aggressive"), (M17, "aggressive")]:
        r = results[key]
        print(key, {k: round(v, 2) for k, v in r["lines"].items()}, "eng", round(r["eng"], 2), "cont", round(r["cont"], 2))
        print("   opex", {k: round(v / 1e6, 3) for k, v in r["opex"].items()})
    try:
        import make_figure
        make_figure.make(wf_rows, out_rows)
    except ImportError:
        pass
    return results


if __name__ == "__main__":
    main()
