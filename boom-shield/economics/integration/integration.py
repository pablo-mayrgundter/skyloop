#!/usr/bin/env python3
"""BoomShield vertical-integration study (screening level).

Run:  python3 integration.py        (light CPU, < 5 s)
Writes next to this file: bom.csv, plants.csv, learning.csv, results.csv, bcr.csv, integration.png.

Imports ../cost_model.py and ../cost-down/cost_down.py READ-ONLY (sys.dont_write_bytecode; module state is
overridden in memory inside context managers and restored).  The integrated scenario replaces cost-down's lever 2
("materials & roll-to-roll manufacturing") by bottom-up unit costs from owned plants plus Wright's-law learning,
and keeps cost-down's other levers (mass, NOAK overheads, life, operations, dual-use) at plausible or aggressive
settings.  All scenarios here carry the two-row tether field (operations: ~2,700 tethers per 4,000 km).

[U] = unverified engineering estimate.  Web figures are search summaries (sources.csv), not read at source.
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
sys.path.insert(0, str(ECON / "cost-down"))
import cost_model as cm  # noqa: E402
import cost_down as cd   # noqa: E402

cm.load_operations()
M20, M17 = cd.M20, cd.M17
R = 0.07
SPINE_KM = 3614.0
BUILD_YR = 10.0
DEMAND = BS / "network-maps" / "demand" / "scenario_results.csv"
HIGHMACH = BS / "network-maps" / "highmach"


def crf(r, n):
    return r * (1 + r) ** n / ((1 + r) ** n - 1)


def cap_charge(capex, n, salvage):
    """Annual capital charge for a plant used n years, resold at `salvage` x capex at the end."""
    return (capex - salvage * capex / (1 + R) ** n) * crf(R, n)


def wright_avg(lr, n_units, pilot_share=0.01):
    """Cumulative-average cost over units x0..N relative to the cost of unit x0 (end of the pilot), as cost-down."""
    b = -math.log2(1 - lr)
    x0 = max(n_units * pilot_share, 1.0)
    N = n_units
    if N <= x0:
        return 1.0
    return (N ** (1 - b) - x0 ** (1 - b)) / ((1 - b) * (N - x0)) / x0 ** (-b)


# ------------------------------------------------------------------------------------------------
# Integrated-plant parameters: (low cost, central, high cost).  [U] unless a source id is given.
# ------------------------------------------------------------------------------------------------
P = {
    # resin: ethane-indexed long-term offtake (cash cost + capital fee), not an owned cracker
    "resin_usd_kg": (0.70, 0.80, 0.95),          # S17: ethylene cash ~$0.31/kg; LLDPE spot $0.85-0.95 (econ [1])
    "cracker_capex_usd_t": (4000, 4500, 6000),   # S5/S6: $8.5 B / 2 Mt (GTP); Monaca $6-10 B / 1.6 Mt
    "cracker_cash_usd_kg": (0.45, 0.52, 0.60),   # ethylene cash + PE conversion [U]
    "filler_usd_kg": (0.10, 0.15, 0.25),         # S20: GCC $120/t FOB Gulf; $0.22-0.39/kg regional
    "filler_frac": 0.60,                         # CaCO3/ATH mass fraction in ballast ply (cost-down §1)
    # extrusion
    "film_line_tph": (2.2, 2.0, 1.5),            # S7: 2.1 t/h at 18.5 m; >2 t/h at 24 m
    "film_line_capex_musd": (8.0, 12.0, 18.0),   # incl. winders, building share, QC [U] (no new-line price found)
    "sheet_line_tph": (5.0, 4.0, 3.0),           # direct-compounding flat-die sheet, 8-10 m [U]
    "sheet_line_capex_musd": (10.0, 15.0, 22.0), # incl. twin-screw, gravimetric filler feed [U]
    "hours_yr": 7000.0,                          # 80 % OEE
    "fte_line": (12, 15, 20), "fte_usd": (75e3, 85e3, 100e3),
    "kwh_kg": 0.35, "power_usd_kwh": 0.07, "maint_frac": 0.03, "scrap": 0.03, "additive_usd_kg": 0.05,
    "plant_salvage": (0.30, 0.20, 0.10),
    # UV skin: fluoropolymer resin bought, film cast in-house (12-25 um)
    "skin_usd_m2": (0.35, 0.55, 0.90),           # S25: PVF resin ~$10-11/kg -> ~$0.24/m2 resin at 12.5 um
    # UHMWPE fibre (gel spinning)
    "fibre_capex_usd_tpy": (12e3, 20e3, 35e3),   # S8/S9: RMB 82-94 k per t/yr in China = $11.5-13 k; West x2-3 [U]
    "uhmwpe_resin_usd_kg": (1.8, 2.3, 3.0),      # [U]
    "solvent_usd_kg": (0.3, 0.5, 0.8),           # make-up + recovery [U]
    "fibre_kwh_kg": (8.0, 12.0, 18.0),           # spinning, extraction, multi-stage drawing, distillation [U]
    "fibre_fte_kt": (10.0, 15.0, 30.0),          # FTE per kt/yr [U]
    "fibre_salvage": (0.30, 0.20, 0.05),
    "rope_conv_usd_kg": (0.8, 1.25, 2.0),        # braiding+jacket (tethers) / UD tape (tapes) [U]
    # H2: bought under a build-own-operate SMR contract
    "h2_usd_kg": (1.2, 1.6, 2.0), "h2_deliv_usd_kg": (0.3, 0.5, 0.8),   # S22 grey $1.1-2.4
    # module factory (seams, lay-up, partitions, ties, fold/pack)
    "seam_m_per_m2": 0.29,                       # 9 m laminated web + 16 m floor/skin strips + partition joins
    "seam_usd_m_pilot": (0.6, 1.0, 1.6),         # automated dual-track weld + air-channel test [U]
    "lr_seam": (0.15, 0.10, 0.07),
    "fab_h_m2_pilot": (0.004, 0.006, 0.010),     # lay-up, ties, partitions, fold, QA; field geomembrane ~0.013 [U]
    "fab_usd_h": (50.0, 55.0, 70.0),
    "lr_module": (0.20, 0.15, 0.10),
    "factory_capex_musd": (450.0, 700.0, 1100.0),   # ~6 lines of 300 x 80 m halls + SPMT handling [U]
    "ties_usd_m2": (0.07, 0.10, 0.15), "partition_fab_usd_m2_wall": (0.03, 0.05, 0.08),
    "partition_film_kg_m2_wall": 0.046,         # 50 um
    # launch bases
    "n_bases": 12, "base_capex_musd": (120.0, 170.0, 300.0), "robots_capex_musd": (60.0, 100.0, 200.0),
    "launch_pilot_usd_m2": (0.35, 0.55, 0.90),   # launch, tow, zip, tether raise crews at pilot rate [U]
    "lr_launch": (0.20, 0.15, 0.10),
    # tether stations (standard, unstaffed): pilot unit cost and learnable share
    "station_pilot_musd": (7.0, 9.0, 12.0), "station_learn_share": 0.70, "lr_station": (0.15, 0.12, 0.08),
    "station_factory_capex_musd": (80.0, 150.0, 250.0),
    "anchor_learn_share": 0.50, "lr_anchor": (0.10, 0.08, 0.05),
    # controls: own design, EMS-built
    "cell_pilot_usd": (40e3, 50e3, 70e3), "cell_learn_share": 0.80, "lr_cell": (0.22, 0.18, 0.12),
}
IDX = {"low": 0, "central": 1, "high": 2}


def p(k, case):
    v = P[k]
    return v[IDX[case]] if isinstance(v, tuple) else v


# Prices the project would pay buying at its own scale = cost-down's lever-2 settings
BUY = {
    "plausible": {"film": 2.5, "ballast": 1.2, "uhmwpe": 22.0, "floor": 4.0, "skin": 1.5, "seams": 3.0,
                  "partition": 0.6, "ties": 0.35, "install": 3.5, "station": 15.0, "cell": 75e3, "h2": 3.0 + 1.0},
    "aggressive": {"film": 1.8, "ballast": 0.8, "uhmwpe": 15.0, "floor": 2.5, "skin": 1.0, "seams": 1.5,
                   "partition": 0.4, "ties": 0.25, "install": 2.0, "station": 10.0, "cell": 40e3, "h2": 2.0 + 0.5},
}


# ------------------------------------------------------------------------------------------------
def lever_set(level):
    """cost-down's full stack at `level` (plausible/aggressive) with two tether rows."""
    S = copy.deepcopy(cd.scenario(level)[-1][1])
    S["two_rows"] = True
    return S


def quantities(design, S, km=SPINE_KM):
    """Physical bill of materials for `km` of corridor under cost-down scenario S."""
    r = cd.evaluate(design, S)
    A = r["A"]
    lid = r["lid_kg_m2"]
    tape = cm.v("tape_kg_m2", "base")
    bulk = max(lid - tape - 0.05, 0)
    s_kg = min(S["struct_kg"], bulk) if S["struct_kg"] is not None else bulk
    b_kg = bulk - s_kg
    if S["pv_frac"] > 0 and S["pv"]:
        b_kg = max(b_kg - S["pv_frac"] * S["pv"]["mass_kg_m2"], 0)
    T = r["T"]
    wall = T * 1.5 / 100 + T * 1.5 / 250
    area = A * km                                     # m2
    rows = 2.0 if S.get("two_rows") else 1.0
    teth_t_km = r["tether_t_km"] * rows
    n_sites = rows * km / S["spacing_km"]
    q = {
        "area_m2": area, "A": A, "T": T, "lid_kg_m2": lid, "struct_kg_m2": s_kg, "ballast_kg_m2": b_kg,
        "floor_kg_m2": cm.FLOOR_KG, "partition_kg_m2": wall * P["partition_film_kg_m2_wall"], "wall_m2_m2": wall,
        "tape_kg_m2": tape, "tether_t_km": teth_t_km, "h2_t_km": r["h2_t_km"] * cm.v("h2_fill_loss", "base"),
        "film_t": area * (s_kg + cm.FLOOR_KG + wall * P["partition_film_kg_m2_wall"]) / 1e3,
        "ballast_t": area * b_kg / 1e3,
        "uhmwpe_t": area * tape / 1e3 + teth_t_km * km,
        "tether_fibre_t": teth_t_km * km, "tape_fibre_t": area * tape / 1e3,
        "h2_t": r["h2_t_km"] * cm.v("h2_fill_loss", "base") * km,
        "modules": km * 1000 / 250.0,                # one module = 250 m x full width
        "sites": n_sites, "cells": area / 62500.0,
        "seam_m": area * P["seam_m_per_m2"],
    }
    q["polymer_t"] = q["film_t"] + q["ballast_t"] * (1 - P["filler_frac"]) + area * 0.021 / 1e3   # + skin
    return q, r


# ------------------------------------------------------------------------------------------------
def plants(q, case="central", build_yr=BUILD_YR, learning=True):
    """Integrated plant models.  Returns (unit costs for cost_down, plant rows)."""
    n = build_yr
    hrs = P["hours_yr"]
    wavg = wright_avg if learning else (lambda lr, N, pilot_share=0.01: 1.0)
    sal = p("plant_salvage", case)
    rows = []
    resin = p("resin_usd_kg", case)
    fte = p("fte_usd", case)
    pw = P["kwh_kg"] * P["power_usd_kwh"]

    # 1. Resin: owned cracker shown for comparison only (central plan = offtake)
    poly_tpy = q["polymer_t"] / n
    cr_capex = p("cracker_capex_usd_t", case) * poly_tpy
    cr_full = p("cracker_cash_usd_kg", case) + cap_charge(cr_capex, 25.0, 0.0) / (poly_tpy * 1e3)
    rows.append({"plant": "PE resin: ethane-indexed offtake (owned cracker shown)", "make_or_buy": "buy (contract)",
                 "capacity": f"{poly_tpy/1e6:.2f} Mt/yr polymer", "n_units": 1, "capex_musd": 0.0,
                 "alt_capex_musd": round(cr_capex / 1e6), "unit": "$/kg", "integrated": resin,
                 "alt_owned_full_cost": round(cr_full, 2), "buy_plausible": 0.92, "buy_aggressive": 0.90,
                 "volume_per_yr": poly_tpy * 1e3,
                 "post_build": "cracker would be a merchant PE plant (sellable, 25-40 yr life): no stranding, no edge"})

    # 2. Film lines (structural ply, floor, partitions)
    film_tpy = q["film_t"] / n
    film_line_tpy = p("film_line_tph", case) * hrs
    n_film = math.ceil(film_tpy / film_line_tpy)
    film_capex = n_film * p("film_line_capex_musd", case) * 1e6
    film_conv = (cap_charge(film_capex, n, sal) + n_film * (p("fte_line", case) * fte)
                 + film_capex * P["maint_frac"]) / (film_tpy * 1e3) + pw
    film_cost = resin * (1 + P["scrap"]) + P["additive_usd_kg"] + film_conv
    rows.append({"plant": "Blown film lines (structural ply, floor, partitions)", "make_or_buy": "make",
                 "capacity": f"{film_tpy/1e3:.0f} kt/yr", "n_units": n_film, "capex_musd": film_capex / 1e6,
                 "unit": "$/kg", "integrated": film_cost, "buy_plausible": BUY["plausible"]["film"],
                 "buy_aggressive": BUY["aggressive"]["film"], "volume_per_yr": film_tpy * 1e3,
                 "post_build": "skins, network arms; then ag-film/geomembrane market (~6 Mt/yr ag film)"})

    # 3. Ballast sheet lines (direct compounding + extrusion lamination of the structural ply)
    bal_tpy = q["ballast_t"] / n
    sheet_tpy = p("sheet_line_tph", case) * hrs
    n_sheet = math.ceil(bal_tpy / sheet_tpy)
    sheet_capex = n_sheet * p("sheet_line_capex_musd", case) * 1e6
    sheet_conv = (cap_charge(sheet_capex, n, sal) + n_sheet * p("fte_line", case) * fte
                  + sheet_capex * P["maint_frac"]) / (bal_tpy * 1e3) + 0.30 * P["power_usd_kwh"]
    f = P["filler_frac"]
    bal_cost = ((1 - f) * resin + f * p("filler_usd_kg", case) + 0.02) * (1 + P["scrap"]) + sheet_conv
    rows.append({"plant": "Filled-PE ballast sheet lines (in-line compounding + lamination)", "make_or_buy": "make",
                 "capacity": f"{bal_tpy/1e6:.2f} Mt/yr", "n_units": n_sheet, "capex_musd": sheet_capex / 1e6,
                 "unit": "$/kg", "integrated": bal_cost, "buy_plausible": BUY["plausible"]["ballast"],
                 "buy_aggressive": BUY["aggressive"]["ballast"], "volume_per_yr": bal_tpy * 1e3,
                 "post_build": "little: replacement ballast only at 30-yr life; resell lines (geomembrane)"})

    # 4. UV skin film (bought fluoropolymer resin)
    skin = p("skin_usd_m2", case)
    rows.append({"plant": "UV skin casting (fluoropolymer resin bought)", "make_or_buy": "make film, buy resin",
                 "capacity": f"{q['area_m2']/n/1e6:.0f} M m2/yr (~13 kt/yr resin)", "n_units": 3,
                 "capex_musd": 3 * 15.0, "unit": "$/m2", "integrated": skin,
                 "buy_plausible": BUY["plausible"]["skin"], "buy_aggressive": BUY["aggressive"]["skin"],
                 "volume_per_yr": q["area_m2"] / n,
                 "post_build": "re-skinning every 8-10 yr keeps it ~full: the one plant with a lasting job"})

    # 5. UHMWPE fibre + rope/tape conversion
    fib_tpy = q["uhmwpe_t"] / n
    fib_capex = p("fibre_capex_usd_tpy", case) * fib_tpy
    fib_cost = (p("uhmwpe_resin_usd_kg", case) * 1.05 + p("solvent_usd_kg", case)
                + p("fibre_kwh_kg", case) * P["power_usd_kwh"]
                + p("fibre_fte_kt", case) * fib_tpy / 1e3 * fte / (fib_tpy * 1e3)
                + (cap_charge(fib_capex, n, p("fibre_salvage", case)) + fib_capex * P["maint_frac"]) / (fib_tpy * 1e3))
    rope_cost = fib_cost + p("rope_conv_usd_kg", case)
    conv_capex = 0.25e9 * (fib_tpy / 100e3)          # braiders, tape lines, buildings [U]
    rows.append({"plant": "UHMWPE gel-spun fibre (+ braiding, UD tape)", "make_or_buy": "make",
                 "capacity": f"{fib_tpy/1e3:.0f} kt/yr fibre", "n_units": math.ceil(fib_tpy / 1e3),
                 "capex_musd": (fib_capex + conv_capex) / 1e6, "unit": "$/kg", "integrated": rope_cost,
                 "fibre_only": fib_cost, "buy_plausible": BUY["plausible"]["uhmwpe"],
                 "buy_aggressive": BUY["aggressive"]["uhmwpe"], "volume_per_yr": fib_tpy * 1e3,
                 "post_build": "tether replacement ~1/25 per yr + arms; else doubles world supply -> stranded"})

    # 6. H2 (bought, BOO contract) - own SMR shown for comparison
    h2 = p("h2_usd_kg", case) + p("h2_deliv_usd_kg", case)
    h2_tpy = q["h2_t"] / n
    smr_capex = 3500.0 * h2_tpy                       # ~$0.75 B for ~0.2 Mt/yr [U]
    rows.append({"plant": "H2: BOO SMR/ATR contract at the bases (own SMR shown)", "make_or_buy": "buy (contract)",
                 "capacity": f"{h2_tpy/1e3:.0f} kt/yr fill (+make-up)", "n_units": 1, "capex_musd": 0.0,
                 "alt_capex_musd": round(smr_capex / 1e6), "unit": "$/kg delivered", "integrated": h2,
                 "buy_plausible": BUY["plausible"]["h2"], "buy_aggressive": BUY["aggressive"]["h2"],
                 "volume_per_yr": h2_tpy * 1e3,
                 "post_build": "make-up 3-6 %/yr; electrolysis 2.4 GW would cost ~$3-4/kg: do not build"})

    # 7. Module factory: seams, lay-up, partitions, ties, fold and pack
    seam = p("seam_usd_m_pilot", case) * P["seam_m_per_m2"] * wavg(p("lr_seam", case), q["seam_m"])
    lab = p("fab_h_m2_pilot", case) * p("fab_usd_h", case) * wavg(p("lr_module", case), q["modules"])
    fac_capex = p("factory_capex_musd", case) * 1e6
    fac = (cap_charge(fac_capex, n, sal) + fac_capex * P["maint_frac"]) / (q["area_m2"] / n)
    seams_fab = seam + lab + fac + 0.03
    partition = film_cost * P["partition_film_kg_m2_wall"] + p("partition_fab_usd_m2_wall", case)
    rows.append({"plant": "Module factory (seams, lay-up, ties, fold/pack)", "make_or_buy": "make",
                 "capacity": f"{q['modules']/n:.0f} modules/yr ({q['seam_m']/n/1e6:.0f} M seam-m/yr)",
                 "n_units": 6, "capex_musd": fac_capex / 1e6, "unit": "$/m2", "integrated": seams_fab,
                 "buy_plausible": BUY["plausible"]["seams"], "buy_aggressive": BUY["aggressive"]["seams"],
                 "volume_per_yr": q["area_m2"] / n,
                 "post_build": "spare modules, arms; no outside market for 1,500 t membrane modules"})

    # 8. Launch bases and reusable ground handling
    lb_capex = (P["n_bases"] * p("base_capex_musd", case) + p("robots_capex_musd", case)) * 1e6
    launch = (p("launch_pilot_usd_m2", case) * wavg(p("lr_launch", case), q["modules"])
              + (cap_charge(lb_capex, n, 0.15) + lb_capex * P["maint_frac"]) / (q["area_m2"] / n))
    rows.append({"plant": f"Launch bases ({P['n_bases']}) + tugs, zip robots", "make_or_buy": "make (own)",
                 "capacity": f"{q['modules']/n/100:.0f} modules per window day (100 days/yr)",
                 "n_units": P["n_bases"], "capex_musd": lb_capex / 1e6, "unit": "$/m2", "integrated": launch,
                 "buy_plausible": BUY["plausible"]["install"], "buy_aggressive": BUY["aggressive"]["install"],
                 "volume_per_yr": q["area_m2"] / n,
                 "post_build": "replacement modules and re-skinning; bases stay as maintenance depots"})

    # 9. Tether stations: standard, unstaffed, factory-built kits
    st_learn = wavg(p("lr_station", case), q["sites"])
    sh = P["station_learn_share"]
    st_fac = p("station_factory_capex_musd", case) * 1e6
    station = (p("station_pilot_musd", case) * (1 - sh + sh * st_learn)
               + (cap_charge(st_fac, n, sal) + st_fac * P["maint_frac"]) / (q["sites"] / n) / 1e6)
    rows.append({"plant": "Tether-station kit factory (winch, brake, H2 skid, power)", "make_or_buy": "make (assemble)",
                 "capacity": f"{q['sites']/n:.0f} stations/yr", "n_units": 1, "capex_musd": st_fac / 1e6,
                 "unit": "$M/station", "integrated": station, "buy_plausible": BUY["plausible"]["station"],
                 "buy_aggressive": BUY["aggressive"]["station"], "volume_per_yr": q["sites"] / n,
                 "post_build": "spares; mooring/offshore-wind winch market"})
    an = wavg(p("lr_anchor", case), q["sites"])
    anchor_factor = 1 - P["anchor_learn_share"] + P["anchor_learn_share"] * an

    # 10. Controls per 250 x 250 m cell: own design, contract manufacture
    cl = wavg(p("lr_cell", case), q["cells"])
    cell = p("cell_pilot_usd", case) * (1 - P["cell_learn_share"] + P["cell_learn_share"] * cl)

    unit = {"film_usd_kg": film_cost, "ballast_usd_kg": bal_cost, "floor_usd_kg": film_cost, "skin_usd_m2": skin,
            "uhmwpe_usd_kg": rope_cost, "seams_usd_m2": seams_fab, "partition_usd_m2": partition,
            "ties_usd_m2": p("ties_usd_m2", case), "install_usd_m2": launch, "station_musd": station,
            "controls_usd_cell": cell, "h2_usd_kg": p("h2_usd_kg", case), "h2_deliv_usd_kg": p("h2_deliv_usd_kg", case),
            "_anchor_factor": anchor_factor, "_fibre_usd_kg": fib_cost}
    return unit, rows


@contextlib.contextmanager
def anchor_scale(f):
    saved = copy.deepcopy(cm.P)
    try:
        for k in ("concrete_usd_m3", "anchor_hw_musd"):
            lo, b, hi = cm.P[k][0]
            cm.P[k] = ((lo, b * f, hi),) + tuple(cm.P[k][1:])
        yield
    finally:
        cm.P.clear()
        cm.P.update(saved)


def integrated_eval(design, level, case="central", km=SPINE_KM, build_yr=None, learning=True):
    """cost-down stack at `level` with lever 2 replaced by integrated unit costs."""
    S = lever_set(level)
    q, _ = quantities(design, S, km=km)                # network: same plants (spine build rate) run longer,
    by = build_yr if build_yr is not None else BUILD_YR * km / SPINE_KM   # over more units (more learning)
    unit, rows = plants(q, case, build_yr=by, learning=learning)
    S2 = copy.deepcopy(S)
    for k, v in unit.items():
        if not k.startswith("_"):
            S2[k] = v
    with anchor_scale(unit["_anchor_factor"]):
        r = cd.evaluate(design, S2)
    return r, unit, rows, q


# ------------------------------------------------------------------------------------------------
def benefits():
    """Benefit per km/yr at 100 % capture (fixed and base demand) for spine, knee and shoulder networks, M2.0/M1.7.
    joint_costdown.csv convention: B/C = (B/C at cost_model base) x (base annual cost / new net annual cost)."""
    out = {}
    with open(DEMAND, newline="") as fh:
        for r in csv.DictReader(fh):
            if r["capture"] != "1" or r["scenario"] != "base":
                continue
            km = float(r["km"])
            name = "spine" if abs(km - 3613.6) < 1 else ("knee 5,421 km" if abs(km - 5421.4) < 1 else
                                                         ("shoulder 6,989 km" if abs(km - 6989) < 3 else None))
            if name is None:
                continue
            mach = float(r["mach"])
            ref = 23.34 if mach == 2.0 else 17.27            # cost_model base $M/km/yr used by the demand study
            out[(name, mach, r["vs"], "fixed")] = float(r["bcr_fixed"]) * ref
            out[(name, mach, r["vs"], "base")] = float(r["bcr"]) * ref
            out[(name, mach, "km")] = km
    return out


def main():
    out_bom, out_plants, out_learn, out_res, out_bcr = [], [], [], [], []

    # ---- 1. BOM per km, ranked (M2.0) ----
    ref = {}
    for level in ("base", "plausible", "aggressive"):
        S = copy.deepcopy(cd.scenario(level)[-1][1]) if level != "base" else copy.deepcopy(cd.BASE)
        S["two_rows"] = True
        r = cd.evaluate(M20, S)
        ref[level] = r
    q, _ = quantities(M20, lever_set("aggressive"), km=1.0)
    qty = {"Lid mass (PE sheet/laminate)": f"{q['lid_kg_m2']*q['A']/1e3:,.0f} t/km lid ({q['struct_kg_m2']:.2f} film + {q['ballast_kg_m2']:.2f} ballast kg/m2)",
           "Lid UV/weathering skin": f"{q['A']/1e6:.2f} km2/km (~37 t fluoropolymer)",
           "Load tapes / scrim (UHMWPE)": f"{q['tape_fibre_t']:,.0f} t/km UHMWPE",
           "Floor film": f"{q['floor_kg_m2']*q['A']/1e3:,.0f} t/km film",
           "Partitions and bellows": f"{q['wall_m2_m2']*q['A']/1e6:.1f} km2/km of wall",
           "Soft ties": f"~{q['A']/100:,.0f} ties/km",
           "Seams and fabrication": f"{q['seam_m']/1e3:,.0f} km of seam/km; 4 modules/km",
           "Diurnal gas shuttle": "rim reservoirs, pumps",
           "H2 first fill": f"{q['h2_t']:,.0f} t/km H2",
           "Tethers": f"{q['tether_fibre_t']:,.0f} t/km UHMWPE (2 rows, rounded rim; bare: ~370)",
           "Anchors": f"{q['sites']:.2f} sites/km, gravity blocks",
           "Winch/tension stations": f"{q['sites']:.2f} stations/km",
           "Controls, sensors, transponders": f"{q['cells']:.0f} cells/km",
           "Launch / installation": "4 modules/km",
           }
    lines = list(ref["base"]["lines"].keys())
    tot = {lv: sum(ref[lv]["lines"].values()) + ref[lv]["eng"] + ref[lv]["cont"] + ref[lv]["land"] for lv in ref}
    make = {"Lid mass (PE sheet/laminate)": "MAKE sheet/film; BUY resin (offtake)",
            "Lid UV/weathering skin": "MAKE film; BUY fluoropolymer resin",
            "Load tapes / scrim (UHMWPE)": "MAKE (fibre plant)", "Floor film": "MAKE",
            "Partitions and bellows": "MAKE", "Soft ties": "MAKE (simple moulded pads)",
            "Seams and fabrication": "MAKE (core competence)", "Diurnal gas shuttle": "BUY (pumps, valves)",
            "H2 first fill": "BUY (BOO contract)", "Tethers": "MAKE fibre + braid",
            "Anchors": "BUY civil works (standard design, framework contracts)",
            "Winch/tension stations": "MAKE kit (assemble); BUY winches/drives",
            "Controls, sensors, transponders": "DESIGN in-house; BUY (EMS)",
            "Launch / installation": "MAKE (own bases, crews)"}
    for k in lines:
        row = {"line": k, "quantity_per_km_aggressive_geometry": qty.get(k, ""), "make_or_buy": make.get(k, "")}
        for lv in ("base", "plausible", "aggressive"):
            row[f"{lv}_musd_km"] = round(ref[lv]["lines"][k] * q["A"] / 1e6, 2)
            row[f"{lv}_share"] = round(ref[lv]["lines"][k] / tot[lv], 3)
        out_bom.append(row)
    for name, key in (("Engineering / PM", "eng"), ("Contingency", "cont"), ("Land and easements", "land")):
        row = {"line": name, "quantity_per_km_aggressive_geometry": "", "make_or_buy": "own team" if key == "eng" else ""}
        for lv in ("base", "plausible", "aggressive"):
            row[f"{lv}_musd_km"] = round(ref[lv][key] * q["A"] / 1e6, 2)
            row[f"{lv}_share"] = round(ref[lv][key] / tot[lv], 3)
        out_bom.append(row)
    out_bom.sort(key=lambda x: -x["plausible_share"])

    # ---- 2. plants (M2.0 spine, aggressive levers, central) ----
    r_c, unit_c, rows_c, q_sp = integrated_eval(M20, "aggressive", "central")
    _, unit_lo, rows_lo, _ = integrated_eval(M20, "aggressive", "low")
    _, unit_hi, rows_hi, _ = integrated_eval(M20, "aggressive", "high")
    for rc, rl, rh in zip(rows_c, rows_lo, rows_hi):
        vol = rc["volume_per_yr"]
        u = rc["unit"]
        scale = 1e6 if u == "$M/station" else 1.0
        row = {k: (round(v, 3) if isinstance(v, float) else v) for k, v in rc.items()}
        row["integrated_low"] = round(rl["integrated"], 3)
        row["integrated_high"] = round(rh["integrated"], 3)
        row["capex_low_musd"] = round(rl["capex_musd"])
        row["capex_high_musd"] = round(rh["capex_musd"])
        row["margin_vs_plausible_musd_yr"] = round((rc["buy_plausible"] - rc["integrated"]) * vol * scale / 1e6)
        row["margin_vs_aggressive_musd_yr"] = round((rc["buy_aggressive"] - rc["integrated"]) * vol * scale / 1e6)
        cx = rc["capex_musd"]
        # simple payback on cash savings vs the aggressive buy price (capital charge added back, approx.)
        cash = row["margin_vs_aggressive_musd_yr"] + cx * crf(R, BUILD_YR)
        row["payback_yr_vs_aggressive"] = (round(cx / cash, 1) if cash > 0 else "never") if cx > 0 else ""
        row["capex_musd"] = round(cx)
        out_plants.append(row)

    # ---- 3. learning ----
    analog = {"modules": ("shipbuilding 15-20 %; aircraft 15-20 %; Liberty ships >50 % cut", "lr_module", q_sp["modules"]),
              "seam-metres": ("mature welding process: 7-15 %", "lr_seam", q_sp["seam_m"]),
              "tether stations": ("wind turbines 7-15 % (installed); 14-29 % manufacturing", "lr_station", q_sp["sites"]),
              "control cells": ("PV modules ~20 % (10-53 %); electronics", "lr_cell", q_sp["cells"]),
              "launches": ("as modules; launch-crew learning", "lr_launch", q_sp["modules"]),
              "anchors (civil)": ("construction: 5-10 %", "lr_anchor", q_sp["sites"])}
    for item, (src, key, N) in analog.items():
        for case in ("low", "central", "high"):
            lr = p(key, case)
            out_learn.append({"item": item, "units_spine": round(N), "case": case, "learning_rate": lr,
                              "doublings_after_1pct_pilot": round(math.log2(100), 1),
                              "cum_avg_vs_pilot": round(wright_avg(lr, N), 3),
                              "last_unit_vs_pilot": round((1 - lr) ** math.log2(100), 3), "analogue": src})

    # ---- 4. results: per m2, per km, annual ----
    def res_row(design, label, r, extra=""):
        return {"design": design, "scenario": label, "capex_usd_m2": round(r["usd_m2"], 2),
                "capex_musd_km": round(r["capex_musd_km"], 1), "opex_musd_km_yr": round(r["opex_musd_km_yr"], 2),
                "annual_musd_km_yr": round(r["annual_musd_km_yr"], 2), "dual_use_musd_km_yr": round(r["dual_musd_km_yr"], 2),
                "net_musd_km_yr": round(r["net_musd_km_yr"], 2), "net_usd_m2_yr": round(r["net_usd_m2_yr"], 2),
                "note": extra}

    results = {}
    for design in (M20, M17):
        for level in ("plausible", "aggressive"):
            for tr in (False, True):
                S = copy.deepcopy(cd.scenario(level)[-1][1])
                S["two_rows"] = tr
                r = cd.evaluate(design, S)
                lab = f"cost-down {level}" + (" + 2nd tether row" if tr else " (1 row, as published)")
                out_res.append(res_row(design, lab, r))
                results[(design, lab)] = r
        b = cd.evaluate(design, dict(copy.deepcopy(cd.BASE), two_rows=True))
        out_res.insert(len(out_res) - 4, res_row(design, "base + 2nd tether row", b))
        results[(design, "base + 2nd tether row")] = b
        for level in ("plausible", "aggressive"):
            for case in ("low", "central", "high"):
                r, unit, rows, q = integrated_eval(design, level, case)
                lab = f"integrated, {level} levers, {case}"
                out_res.append(res_row(design, lab, r, "2 tether rows; lever 2 replaced by owned plants + learning"))
                results[(design, lab)] = r
        r, *_ = integrated_eval(design, "aggressive", "central", learning=False)
        out_res.append(res_row(design, "integrated, aggressive levers, central, NO learning", r,
                               "pilot-rate unit costs throughout"))
        results[(design, "integrated no learning")] = r
        # bare tethers (no rim rounding) sensitivity on the central integrated aggressive case
        S = lever_set("aggressive")
        S["tether_f"] = 1.0
        S["spacing_km"] = 3.0
        q, _ = quantities(design, S)
        unit, _ = plants(q, "central")
        S2 = copy.deepcopy(S)
        S2.update({k: v for k, v in unit.items() if not k.startswith("_")})
        with anchor_scale(unit["_anchor_factor"]):
            r = cd.evaluate(design, S2)
        out_res.append(res_row(design, "integrated, aggressive levers, central, BARE tethers (no rim rounding)", r,
                               "fibre plant x1.9"))
        results[(design, "integrated bare")] = r

    # ---- 5. B/C: spine and networks, vs Boomless and subsonic, 10 % / 100 % capture ----
    ben = benefits()
    for design, mach in ((M20, 2.0), (M17, 1.7)):
        for net in ("spine", "knee 5,421 km", "shoulder 6,989 km"):
            if (net, mach, "km") not in ben:     # the demand study ran networks at M2.0 only
                continue
            km = ben[(net, mach, "km")]
            cases = [("base + 2nd tether row", results[(design, "base + 2nd tether row")]),
                     ("cost-down plausible + 2nd row", results[(design, "cost-down plausible + 2nd tether row")]),
                     ("cost-down aggressive + 2nd row", results[(design, "cost-down aggressive + 2nd tether row")]),
                     ("cost-down aggressive (1 row, as published)", results[(design, "cost-down aggressive (1 row, as published)")])]
            for level in ("plausible", "aggressive"):
                r, *_ = integrated_eval(design, level, "central", km=km)
                cases.append((f"integrated, {level} levers (central)", r))
            r, *_ = integrated_eval(design, "aggressive", "low", km=km)
            cases.append(("integrated, aggressive levers (low cost)", r))
            for label, r in cases:
                for dem in ("fixed", "base"):
                    bb = ben[(net, mach, "bl", dem)]
                    bs = ben[(net, mach, "sub", dem)]
                    c = r["net_musd_km_yr"]
                    out_bcr.append({"design": design, "network": net, "km": round(km), "cost_case": label,
                                    "demand": dem, "net_musd_km_yr": round(c, 2),
                                    "bcr_bl_10": round(0.1 * bb / c, 4), "bcr_bl_100": round(bb / c, 3),
                                    "bcr_sub_10": round(0.1 * bs / c, 4), "bcr_sub_100": round(bs / c, 3)})
    hm = highmach_hook()
    if hm:
        out_bcr.extend(hm)

    # ---- write ----
    def w(name, rows):
        if not rows:
            return
        keys = []
        for r in rows:
            for k in r:
                if k not in keys:
                    keys.append(k)
        with open(HERE / name, "w", newline="") as fh:
            wr = csv.DictWriter(fh, fieldnames=keys)
            wr.writeheader()
            wr.writerows(rows)

    w("bom.csv", out_bom)
    w("plants.csv", out_plants)
    w("learning.csv", out_learn)
    w("results.csv", out_res)
    w("bcr.csv", out_bcr)

    # console summary
    print("quantities, spine M2.0 aggressive geometry:",
          {k: (f"{v:.3g}") for k, v in q_sp.items()})
    print("unit costs central:", {k: round(v, 3) for k, v in unit_c.items()})
    print("unit costs low:", {k: round(v, 3) for k, v in unit_lo.items()})
    print("unit costs high:", {k: round(v, 3) for k, v in unit_hi.items()})
    print("integrated M2.0 aggressive central lines:", {k: round(v, 2) for k, v in r_c["lines"].items()},
          "eng", round(r_c["eng"], 2), "cont", round(r_c["cont"], 2))
    print("opex:", {k: round(v / 1e6, 3) for k, v in r_c["opex"].items()})
    for r in out_res:
        print(r)
    try:
        import make_figure_int
        make_figure_int.make(HERE, results, out_bcr, out_learn, out_plants)
    except ImportError as e:
        print("figure skipped:", e)
    return results


def highmach_hook():
    """Sibling study on higher-Mach aircraft (network-maps/highmach/).  If it has written a B/C table with a
    'bcr_' column for the spine at base cost, rescale it by the integrated cost; otherwise return []."""
    if not HIGHMACH.exists():
        return []
    rows = []
    for f in sorted(HIGHMACH.glob("*.csv")):
        if f.name == "evidence.csv":
            continue
        with open(f, newline="") as fh:
            rd = csv.DictReader(fh)
            cols = rd.fieldnames or []
            if not any(c.startswith("bcr") for c in cols):
                continue
            rows.append({"design": "highmach (sibling)", "network": f.name, "km": "", "cost_case": "HOOK: read manually",
                         "demand": "", "net_musd_km_yr": "", "bcr_bl_10": "", "bcr_bl_100": "", "bcr_sub_10": "",
                         "bcr_sub_100": ""})
    return rows


if __name__ == "__main__":
    main()
