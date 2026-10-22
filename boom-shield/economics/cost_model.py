#!/usr/bin/env python3
"""BoomShield first-cut economics (screening level).

Run:  python3 cost_model.py          (writes CSVs + economics.png next to this file)

Everything is per km of corridor unless stated. Three cost cases (low / base / high)
take every parameter at its optimistic / central / pessimistic value together.
Inputs from three sibling studies plug in at the top (PLUG-INS); when a file is
missing the model runs on the placeholder values and reports which were used.

Dollar basis: 2025 USD (no escalation of unit costs between 2023-2025 sources).
"""
from __future__ import annotations

import csv
import math
import os
import re
from pathlib import Path

HERE = Path(__file__).resolve().parent
BS = HERE.parent
CASES = ("low", "base", "high")
ACRE_M2 = 4046.86

# --------------------------------------------------------------------------------------
# PLUG-INS.  Paths to sibling outputs and the placeholder values used when absent.
# --------------------------------------------------------------------------------------
ROUTES_CSV = BS / "network-maps" / "routes" / "routes.csv"
TETHERS_MD = BS / "operations" / "tethers-build-failure.md"
FOCUS_MD = BS / "sim" / "results" / "focus_booms.md"

# Tether field (operations/tethers-build-failure.md).  Placeholder = altitude-trade §6
# scaled to the psf layer by STATUS §6/§7: 17 MN, 331 t per tether at 3 km spacing for
# 30 m x 1 km; ~26 MN for 54 m x 1.5 km  ->  mass ∝ tension  ->  ~506 t per tether.
TETHER = {
    "source": "placeholder (altitude-trade §6 scaled to 26 MN per STATUS §7)",
    "tension_MN": 26.0,            # top tension per tether, Mach-2 layer, 3 km spacing
    "mass_t_per_tether": 506.0,    # bare UHMWPE
    "faired": False,               # bare baseline; faired variant below
    "faired_tension_factor": 0.90, # fairing removes most of ~3 MN tether drag (alt-trade §6)
    "fairing_cost_factor": (0.3, 0.6, 1.2),  # fairing cost as fraction of bare tether cost
    "anchor_footprint_ha": (2.0, 4.0, 8.0),  # hard area per site (anchor, winch hall, substation, access)
    "exclusion_radius_m": (300.0, 500.0, 1000.0),  # recoil / drop exclusion easement per site
}
# Acceleration / superboom rules (sim/results/focus_booms.md).  Placeholder: the climb
# and acceleration ramps at each end need shield but give little time saving, and the
# focus zone needs extra width.  Expressed as extra shielded km per corridor end.
FOCUS = {
    "source": "placeholder (no focus_booms.md yet)",
    "ramp_km_per_end": 250.0,      # km of shield at each end flown below cruise Mach
    "ramp_width_factor": 1.3,      # wider shield under the acceleration focus
}

OPS_DIR = BS / "operations"


def load_operations():
    """Pull tether mass/tension, fairing mass, anchor block and site area from the
    operations agent's CSVs (tether_sizing.csv, failure_landuse.csv) when present."""
    msgs = []
    f = OPS_DIR / "tether_sizing.csv"
    if f.exists():
        with open(f, newline="") as fh:
            rows = [r for r in csv.DictReader(fh) if r["layer"].startswith("54m") and float(r["rim_CD"]) == 1.0]
        by = {r["tether"]: r for r in rows}
        if "bare" in by:
            TETHER["mass_t_per_tether"] = float(by["bare"]["total_t"])
            TETHER["tension_MN"] = float(by["bare"]["T_top_MN"])
            msgs.append(f"bare tether {TETHER['mass_t_per_tether']:.0f} t, {TETHER['tension_MN']:.1f} MN (54 m x 1.5 km, rim C_D 1)")
        if "rigid_faired" in by:
            TETHER["faired_fibre_t"] = float(by["rigid_faired"]["fibre_t"])
            TETHER["fairing_t"] = float(by["rigid_faired"]["fairing_t"])
            msgs.append(f"faired: fibre {TETHER['faired_fibre_t']:.0f} t + fairing {TETHER['fairing_t']:.0f} t")
    f = OPS_DIR / "failure_landuse.csv"
    if f.exists():
        with open(f, newline="") as fh:
            q = {r["quantity"]: float(r["value"]) for r in csv.DictReader(fh)}
        if "gravity_block_m3" in q:
            TETHER["gravity_block_m3"] = q["gravity_block_m3"]
            msgs.append(f"gravity anchor {q['gravity_block_m3']:.0f} m3 ({q.get('gravity_block_t', 0):.0f} t)")
        if "site_ha_r100" in q:
            h = q["site_ha_r100"]
            TETHER["anchor_footprint_ha"] = (h * 0.75, h, h * 2)
            msgs.append(f"site {h:.1f} ha (r 100 m)")
    if msgs:
        TETHER["source"] = "operations/*.csv: " + "; ".join(msgs)
    return msgs


# --------------------------------------------------------------------------------------
# Physical design points (STATUS §6, altitude-trade §6-7, structure-concept §2)
# --------------------------------------------------------------------------------------
RHO_H2_20KM = 5529 * 2.016e-3 / (8.314 * 216.65)  # 0.00619 kg/m3 at 20 km, 216.65 K
LIFT_H2_20KM = 0.0827                             # kg lift per m3 (structure-concept Table 2)
HOLDDOWN = 0.55                                   # kg/m2 tether pretension (altitude-trade §6)
FLOOR_KG = 0.10                                   # kg/m2 floor film

DESIGNS = {
    # name: thickness m, width km (low/base/high), side-load scale for tethers
    "M2.0 (Concorde-class, psf baseline)": {"T": 55.0, "W": (1.5, 1.75, 2.0)},
    "M1.7 (Overture-class, lift-set)": {"T": 36.0, "W": (1.5, 1.5, 2.0)},
}
BASE_DESIGN = "M2.0 (Concorde-class, psf baseline)"

# --------------------------------------------------------------------------------------
# Unit costs.  (low, base, high), unit, source tag.  [U] = unverified / engineering guess.
# --------------------------------------------------------------------------------------
P = {
    # Envelope and lid
    "lid_pe_usd_per_kg": ((2.0, 4.0, 10.0), "$/kg",
        "LLDPE resin $0.85-0.95/kg H2 2025 (S&P Platts, IMARC, ChemOrbis) + conversion to "
        "UV/low-temp sheet & cross-ply [U]"),
    "lid_skin_usd_m2": ((1.0, 2.0, 5.0), "$/m2", "UV/weathering top layer [U]"),
    "tape_kg_m2": ((0.10, 0.10, 0.15), "kg/m2", "tether load-spreading tapes, structure-concept §2"),
    "uhmwpe_usd_per_kg": ((15.0, 30.0, 60.0), "$/kg",
        "UHMWPE fibre/rope; no bulk quote retrieved [U]. Retail 3 mm SK78 ~$2.7/m (Binks) "
        "is ~10x higher and not representative"),
    "floor_film_usd_per_kg": ((8.0, 15.0, 30.0), "$/kg", "balloon-grade PE film [U]"),
    "partition_usd_m2_film": ((0.5, 1.0, 2.0), "$/m2 of wall", "slack bellows/compartment film, fabricated [U]"),
    "ties_usd_m2": ((0.2, 0.5, 1.5), "$/m2", "soft ties K~10 Pa/m with dampers [U]"),
    "seams_fab_usd_m2": ((3.0, 6.0, 15.0), "$/m2",
        "seams, gores, edges, tie/tape attachments. Reality check: NASA 40 MCF balloon "
        "~$120k (Air & Space, date unknown) over ~5.3 ha = ~$2.3/m2 all-in for 20 um film"),
    "shuttle_usd_m2": ((0.5, 1.0, 3.0), "$/m2", "diurnal gas shuttle, rim reservoirs ~0.2 kg/m2, pumps (structure §3) [U]"),
    # Hydrogen
    "h2_usd_per_kg": ((1.5, 4.0, 8.0), "$/kg",
        "grey $1.11-2.35, green $3.74-11.70 (BNEF Dec 2024 via BNN Bloomberg/pv mag); "
        "low = grey, base = cheap green, high = mid green"),
    "h2_delivery_usd_per_kg": ((0.5, 1.5, 3.0), "$/kg", "compression, transport, fill to 20 km [U]"),
    "h2_fill_loss": ((1.05, 1.10, 1.20), "x", "venting and fill losses [U]"),
    # Tethers, anchors, stations
    "spacing_km": ((5.0, 3.0, 2.0), "km", "tether spacing 2-5 km (structure §3; ~1,300 tethers/4,000 km)"),
    "anchor_musd_at_3km": ((3.0, 8.0, 20.0), "$M/site",
        "26 MN inclined ground anchor; scales with load (∝ spacing) [U]"),
    "concrete_usd_m3": ((800.0, 1500.0, 3000.0), "$/m3", "reinforced mass concrete in place incl. excavation [U]"),
    "anchor_hw_musd": ((1.0, 2.0, 5.0), "$M/site", "anchor head, swivel, load cells [U]"),
    "fairing_usd_per_kg": ((10.0, 25.0, 60.0), "$/kg", "rigid fairing shells on the tether [U]"),
    "station_musd": ((10.0, 25.0, 60.0), "$M/site",
        "winch/tension control (largest marine winches ~6 MN), power, building, access, "
        "security, Part 101 lighting [U]"),
    # Controls, installation, overheads
    "controls_usd_per_cell": ((50e3, 150e3, 400e3), "$/cell",
        "per 250 m x 250 m H2 compartment: leak detection, valves, sensors, ADS-B, comms [U]"),
    "install_usd_m2": ((2.0, 6.0, 20.0), "$/m2",
        "launch/assembly of a 1.5-2 km sheet at 20 km: no precedent [U]"),
    "engineering_frac": ((0.10, 0.20, 0.30), "x direct", "design, PM, certification [U]"),
    "contingency_frac": ((0.20, 0.50, 1.00), "x (direct+eng)", "first-of-a-kind, AACE Class 5 style -50/+100 % [U]"),
    # Operating
    "h2_makeup_frac_yr": ((0.05, 0.10, 0.20), "/yr", "structure-concept §4 (holes, flex-crack pinholes)"),
    "film_life_yr": ((20.0, 10.0, 5.0), "yr", "UV + 1e5-1e6 slack cycles at 217 K; no coupon data [U]"),
    "tether_life_yr": ((20.0, 15.0, 10.0), "yr", "UHMWPE creep at sustained ~1.2 GPa [U]"),
    "replace_frac": ((0.8, 1.0, 1.2), "x", "replacement cost / original envelope cost [U]"),
    "om_frac_hw_yr": ((0.02, 0.03, 0.05), "/yr", "O&M on stations, anchors, controls [U]"),
    "inspection_musd_km_yr": ((0.05, 0.15, 0.5), "$M/km/yr", "drones/aircraft, envelope survey [U]"),
    "power_kwh_m2_yr": ((0.5, 1.0, 2.0), "kWh/m2/yr", "gas shuttle ~0.7 W/m2 x 4 h/day (structure §3)"),
    "power_usd_kwh": ((0.06, 0.08, 0.12), "$/kWh", "[U]"),
    "station_power_mw": ((0.2, 0.5, 1.0), "MW/site", "winches, heaters, lighting [U]"),
    "fte_per_site": ((4.0, 8.0, 12.0), "FTE", "24/7 site cover [U]"),
    "fte_usd": ((150e3, 180e3, 220e3), "$/FTE/yr", "loaded [U]"),
    "central_ops_musd_km_yr": ((0.02, 0.05, 0.1), "$M/km/yr", "ATC liaison, ops centre [U]"),
    "insurance_frac_yr": ((0.005, 0.01, 0.02), "/yr of capex", "H2 inventory 14-27 kt TNT/km next to traffic [U]"),
    "availability": ((0.98, 0.95, 0.90), "-", "weather/maintenance downtime [U]"),
    # Land
    "land_usd_acre": ((1920.0, 4350.0, 9830.0), "$/acre",
        "USDA NASS Land Values 2025: pasture $1,920, farm real estate $4,350, Pacific cropland $9,830"),
    "fee_multiplier": ((1.5, 3.0, 6.0), "x", "severance, whole-parcel buys, legal (CHSRA experience) [U]"),
    "easement_frac": ((0.3, 0.5, 0.9), "x fee value", "restricted-use easement [U]"),
    "easement_rent_usd_site_yr": ((20e3, 50e3, 150e3), "$/site/yr", "annual payment, wind-turbine style [U]"),
    # Finance
    "discount": ((0.05, 0.07, 0.08), "/yr", "task brief 5-8 %"),
    "asset_life_yr": ((30.0, 25.0, 20.0), "yr", "task brief 20-30 yr"),
}


def v(name, case):
    return P[name][0][CASES.index(case)]


def crf(r, n):
    return r * (1 + r) ** n / ((1 + r) ** n - 1)


# --------------------------------------------------------------------------------------
# Capital cost
# --------------------------------------------------------------------------------------
def design_mass(T):
    lift = LIFT_H2_20KM * T
    lid = lift - HOLDDOWN - FLOOR_KG
    return {"lift": lift, "lid": lid, "h2_kg_m2": RHO_H2_20KM * T}


def tether_t_per_km(T, case):
    """Tether mass per km of corridor; mass ∝ tension ∝ spacing, so per-km mass is
    spacing-independent.  Side load scales ~ (C_D T + 2 W C_f); 18 % of tension is
    tether self-drag (altitude-trade §6: ~3 MN of 17 MN)."""
    base = TETHER["mass_t_per_tether"] / 3.0
    sl = (T + 6.0) / (55.0 + 6.0)
    m = base * (0.82 * sl + 0.18)
    if TETHER["faired"]:
        m *= TETHER["faired_tension_factor"]
    return m


def capex(design=BASE_DESIGN, case="base", width_km=None, usd_m2_scale=1.0, faired=None):
    d = DESIGNS[design]
    T = d["T"]
    W = width_km if width_km is not None else d["W"][CASES.index(case)]
    A = W * 1e6  # m2 per km
    m = design_mass(T)
    s = v("spacing_km", case)
    lines = {}
    # envelope (per m2)
    lid_bulk = max(m["lid"] - v("tape_kg_m2", case) - 0.05, 0)
    lines["Lid mass (PE sheet/laminate)"] = lid_bulk * v("lid_pe_usd_per_kg", case)
    lines["Lid UV/weathering skin"] = v("lid_skin_usd_m2", case)
    lines["Load tapes / scrim (UHMWPE)"] = v("tape_kg_m2", case) * v("uhmwpe_usd_per_kg", case)
    lines["Floor film"] = FLOOR_KG * v("floor_film_usd_per_kg", case)
    wall_ratio = T * 1.5 / 100 + T * 1.5 / 250  # bellows every 100 m + partitions every 250 m
    lines["Partitions and bellows"] = wall_ratio * v("partition_usd_m2_film", case)
    lines["Soft ties"] = v("ties_usd_m2", case)
    lines["Seams and fabrication"] = v("seams_fab_usd_m2", case)
    lines["Diurnal gas shuttle"] = v("shuttle_usd_m2", case)
    h2 = m["h2_kg_m2"] * v("h2_fill_loss", case)
    lines["H2 first fill"] = h2 * (v("h2_usd_per_kg", case) + v("h2_delivery_usd_per_kg", case))
    # tethers etc. (per km -> per m2)
    tt = tether_t_per_km(T, case)
    teth_cost = tt * 1e3 * v("uhmwpe_usd_per_kg", case)
    use_faired = TETHER["faired"] if faired is None else faired
    if use_faired:
        if "fairing_t" in TETHER:   # operations plug-in: fibre + fairing masses
            k = tt / TETHER["mass_t_per_tether"]
            teth_cost = k * 1e3 * (TETHER["faired_fibre_t"] * v("uhmwpe_usd_per_kg", case)
                                   + TETHER["fairing_t"] * v("fairing_usd_per_kg", case))
        else:
            teth_cost *= 1 + TETHER["fairing_cost_factor"][CASES.index(case)]
    lines["Tethers"] = teth_cost / A
    if "gravity_block_m3" in TETHER:  # operations plug-in: concrete gravity block per site at 3 km
        anchor = TETHER["gravity_block_m3"] * v("concrete_usd_m3", case) + v("anchor_hw_musd", case) * 1e6
    else:
        anchor = v("anchor_musd_at_3km", case) * 1e6
    lines["Anchors"] = anchor / 3.0 / A
    lines["Winch/tension stations"] = v("station_musd", case) * 1e6 / s / A
    lines["Controls, sensors, transponders"] = v("controls_usd_per_cell", case) / 62500.0
    lines["Launch / installation"] = v("install_usd_m2", case)
    direct = sum(lines.values()) * usd_m2_scale
    for k in lines:
        lines[k] *= usd_m2_scale
    eng = direct * v("engineering_frac", case)
    cont = (direct + eng) * v("contingency_frac", case)
    lines["Engineering / PM"] = eng
    lines["Contingency"] = cont
    land = land_rights(case)["bs_cost_usd_km"] / A
    lines["Land and easements"] = land
    total_m2 = sum(lines.values())
    return {"design": design, "case": case, "W_km": W, "A_m2_per_km": A, "T": T,
            "lines_usd_m2": lines, "usd_m2": total_m2, "musd_km": total_m2 * A / 1e6,
            "h2_t_km": m["h2_kg_m2"] * A / 1e3, "tether_t_km": tt, "lid_kg_m2": m["lid"],
            "lid_t_km": m["lid"] * A / 1e3, "spacing_km": s}


# --------------------------------------------------------------------------------------
# Operating cost
# --------------------------------------------------------------------------------------
def opex(cx, case="base", film_life=None):
    A = cx["A_m2_per_km"]
    L = cx["lines_usd_m2"]
    s = cx["spacing_km"]
    life = film_life if film_life is not None else v("film_life_yr", case)
    env_keys = ["Lid mass (PE sheet/laminate)", "Lid UV/weathering skin", "Load tapes / scrim (UHMWPE)",
                "Floor film", "Partitions and bellows", "Soft ties", "Seams and fabrication",
                "Launch / installation"]
    env_m2 = sum(L[k] for k in env_keys)
    hw_m2 = sum(L[k] for k in ["Anchors", "Winch/tension stations", "Controls, sensors, transponders",
                               "Diurnal gas shuttle"])
    o = {}
    h2_price = v("h2_usd_per_kg", case) + v("h2_delivery_usd_per_kg", case)
    o["H2 make-up"] = cx["h2_t_km"] * 1e3 * v("h2_makeup_frac_yr", case) * h2_price
    o["Envelope replacement (UV, fatigue)"] = env_m2 * A * v("replace_frac", case) / life
    o["Tether replacement"] = L["Tethers"] * A / v("tether_life_yr", case)
    o["O&M on stations/controls"] = hw_m2 * A * v("om_frac_hw_yr", case)
    o["Inspection"] = v("inspection_musd_km_yr", case) * 1e6
    o["Power"] = (v("power_kwh_m2_yr", case) * A + v("station_power_mw", case) * 8760e3 / s) * v("power_usd_kwh", case)
    o["Staffing"] = v("fte_per_site", case) * v("fte_usd", case) / s + v("central_ops_musd_km_yr", case) * 1e6
    o["Insurance"] = cx["musd_km"] * 1e6 * v("insurance_frac_yr", case)
    o["Easement rent"] = v("easement_rent_usd_site_yr", case) / s
    return o


# --------------------------------------------------------------------------------------
# Land and rights
# --------------------------------------------------------------------------------------
HSR = [
    # project, $M/km total, basis, source
    ("California HSR Merced-Bakersfield (275 km)", 30e3 / 275, 33e3 / 275,
     "$29.8-33.0 B for 171 mi (CHSRA CEO report Feb 2023); IG estimate ~$48 B (~$175M/km)"),
    ("Brightline West (351 km, I-15 median)", 21.5e3 / 351, 21.5e3 / 351,
     "$21.5 B for 218 mi (Railway Age citing Bloomberg, 2025)"),
    ("Spain AVE network average", 17.7 * 1.08, 18.5,
     "EUR 17.7 M/km / ~$18.5 M/km (Global Railway Review)"),
    ("HS2 (UK)", 175.0, 260.0, "~$175 M/km (2023 report); CNN ~$416 M/mile"),
]
HSR_ROW = {
    "ca_cv_row_busd": 1.5, "ca_cv_km": 119 * 1.609,
    "row_width_m": (15.2, 30.5, 36.6),  # 50 ft envelope (FL HSR guidelines), 100 ft, 120 ft (Peninsula)
    "ca_cv_parcels": 2295,
}


def land_rights(case="base"):
    s = v("spacing_km", case)
    i = CASES.index(case)
    hard_ha = TETHER["anchor_footprint_ha"][i]
    r = TETHER["exclusion_radius_m"][i]
    excl_ha = math.pi * r ** 2 / 1e4
    hard_ac_km = hard_ha * 1e4 / ACRE_M2 / s
    excl_ac_km = max(excl_ha - hard_ha, 0) * 1e4 / ACRE_M2 / s
    price = v("land_usd_acre", case)
    cost = hard_ac_km * price * v("fee_multiplier", case) + excl_ac_km * price * v("easement_frac", case)
    hsr_w = HSR_ROW["row_width_m"][i]
    hsr_ac_km = hsr_w * 1000 / ACRE_M2
    hsr_row_cost = HSR_ROW["ca_cv_row_busd"] * 1e9 / HSR_ROW["ca_cv_km"]
    return {"case": case, "bs_hard_acres_km": hard_ac_km, "bs_excl_acres_km": excl_ac_km,
            "bs_cost_usd_km": cost, "hsr_row_acres_km": hsr_ac_km, "hsr_row_cost_usd_km": hsr_row_cost,
            "hsr_parcels_per_km": HSR_ROW["ca_cv_parcels"] / HSR_ROW["ca_cv_km"],
            "bs_sites_per_km": 1 / s}


# --------------------------------------------------------------------------------------
# Benefits
# --------------------------------------------------------------------------------------
A_SOUND = 295.1 * 3.6  # km/h, isothermal 11-20 km
VTTS_2015 = {"personal": 36.10, "business": 63.20}  # DOT 2016 guidance, intercity air, 2015 $
CPI_2015_2025 = 1.36  # CPI-U 2015 avg 237.0 -> 2025 ~322 [approx, unverified]
VTTS = {"low": VTTS_2015["personal"] * CPI_2015_2025,          # ~$49/h
        "base": VTTS_2015["business"] * CPI_2015_2025,         # ~$86/h
        "high": 200.0}                                         # premium-cabin traveller [U]
SEATS = 65            # Overture 64-80 (Boom); 65 used
LOAD_FACTOR = 0.75
FLIGHTS_CAP = 100     # 1e5-1e6 lid cycles over 30 yr -> ~9-90/day (STATUS §7); cap 100/day
FARE_ONEWAY = 2500.0  # Boom target ~$5,000 round trip business (Denver Gazette / OMAAT)
MARKET_PAX_PLACEHOLDER = 10e6  # pax/yr crossing the busiest trunk km (placeholder until routes.csv)
CAPTURE = {"low": 0.05, "base": 0.15, "high": 0.30}


def hours_per_km(mach):
    return 1.0 / (mach * A_SOUND)


def dt_per_km(m_from, m_to):
    return hours_per_km(m_from) - hours_per_km(m_to)


def pax_from_flights(f):
    return f * 2 * SEATS * LOAD_FACTOR * 365 / 2  # f = one-way departures per day each direction


def annual_cost_km(cx, ox, case="base", r=None, n=None):
    r = v("discount", case) if r is None else r
    n = v("asset_life_yr", case) if n is None else n
    return cx["musd_km"] * 1e6 * crf(r, n) + sum(ox.values())


# --------------------------------------------------------------------------------------
# routes.csv loader (flexible on column names)
# --------------------------------------------------------------------------------------
def _find(cols, *keys, avoid=()):
    for c in cols:
        lc = c.lower()
        if all(k in lc for k in keys) and not any(a in lc for a in avoid):
            return c
    return None


def load_routes():
    if not ROUTES_CSV.exists():
        return None, "routes.csv absent: placeholder demand used"
    with open(ROUTES_CSV, newline="") as f:
        rows = list(csv.DictReader(f))
    if not rows:
        return None, "routes.csv empty"
    cols = list(rows[0].keys())
    m = {
        "id": _find(cols, "route") or _find(cols, "id") or cols[0],
        "km": _find(cols, "overland", "km") or _find(cols, "km"),
        "seats": _find(cols, "seat", avoid=("km",)),
        "pax": _find(cols, "pax", avoid=("hour", "km")) or _find(cols, "passeng", avoid=("hour", "km")),
        "ph_sub": _find(cols, "hour", "sub") or _find(cols, "hours_saved", avoid=("incr", "boomless", "1.15", "115")),
        "ph_inc": _find(cols, "hour", "incr") or _find(cols, "hour", "boomless") or _find(cols, "hour", "115"),
    }
    return (rows, m), f"routes.csv read; column map {m}"


def fnum(x):
    try:
        return float(str(x).replace(",", ""))
    except (TypeError, ValueError):
        return float("nan")


# --------------------------------------------------------------------------------------
# Shadow of an opaque strip at 20 km (land-rights note)
# --------------------------------------------------------------------------------------
def shadow_hours(lat=37.0, H=20.0, W=1.75):
    """Max hours/yr any ground point is shaded by an opaque E-W strip of width W (km) at H."""
    import numpy as np
    y = np.arange(-5, 60, 0.25)  # km north of strip south edge
    hrs = np.zeros_like(y)
    maxday = np.zeros_like(y)
    phi = math.radians(lat)
    for doy in range(0, 365, 3):
        dec = math.radians(-23.44 * math.cos(2 * math.pi * (doy + 10) / 365))
        h = np.radians(np.arange(-180, 180, 1.0) * 1.0)  # hour angle, 4-min steps
        sin_alt = np.sin(phi) * np.sin(dec) + np.cos(phi) * np.cos(dec) * np.cos(h)
        up = sin_alt > 0.05
        alt = np.arcsin(np.clip(sin_alt, -1, 1))
        cos_az_n = (np.sin(dec) - sin_alt * np.sin(phi)) / (np.cos(alt) * np.cos(phi))  # cos of azimuth from N
        # shadow displaced away from the sun; N-S component of horizontal offset
        off = -H / np.tan(alt) * cos_az_n  # km north
        day = np.zeros_like(y)
        for o in off[up]:
            day += (y >= o) & (y <= o + W)
        day *= 4 / 60  # hours
        hrs += day * 3
        maxday = np.maximum(maxday, day)
    return float(hrs.max()), float(maxday.max()), float(y[hrs.argmax()])


# --------------------------------------------------------------------------------------
# Main
# --------------------------------------------------------------------------------------
def main():
    notes = []
    load_operations()
    routes, msg = load_routes()
    notes.append(msg)
    notes.append(f"tethers: {TETHER['source']}" + (" | tethers-build-failure.md ABSENT" if not TETHERS_MD.exists() else " | tethers-build-failure.md present"))
    notes.append(f"focus booms: {FOCUS['source']}" + ("" if not FOCUS_MD.exists() else " | focus_booms.md present"))

    # ---- capex ------------------------------------------------------------------
    cap = {(d, c): capex(d, c) for d in DESIGNS for c in CASES}
    keys = list(cap[(BASE_DESIGN, "base")]["lines_usd_m2"].keys())
    with open(HERE / "capex_breakdown.csv", "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["design", "line", *[f"{c}_usd_m2" for c in CASES], *[f"{c}_musd_km" for c in CASES]])
        for d in DESIGNS:
            for k in keys + ["TOTAL"]:
                vals = [cap[(d, c)]["usd_m2"] if k == "TOTAL" else cap[(d, c)]["lines_usd_m2"][k] for c in CASES]
                mk = [x * cap[(d, c)]["A_m2_per_km"] / 1e6 for x, c in zip(vals, CASES)]
                w.writerow([d, k, *[f"{x:.2f}" for x in vals], *[f"{x:.2f}" for x in mk]])
    # faired variant
    fair = {c: capex(BASE_DESIGN, c, faired=True) for c in CASES}

    # ---- opex -------------------------------------------------------------------
    opx = {(d, c): opex(cap[(d, c)], c) for d in DESIGNS for c in CASES}
    with open(HERE / "opex_breakdown.csv", "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["design", "line", *[f"{c}_musd_km_yr" for c in CASES]])
        for d in DESIGNS:
            for k in list(opx[(d, "base")].keys()) + ["TOTAL"]:
                vals = [sum(opx[(d, c)].values()) if k == "TOTAL" else opx[(d, c)][k] for c in CASES]
                w.writerow([d, k, *[f"{x / 1e6:.3f}" for x in vals]])

    # ---- annual cost and break-even (throughput-based, route-length independent) ---
    rows = []
    for d in DESIGNS:
        mach = 2.0 if d.startswith("M2.0") else 1.7
        for c in CASES:
            cx, ox = cap[(d, c)], opx[(d, c)]
            ann = annual_cost_km(cx, ox, c)
            av = v("availability", c)
            for ref, mref in [("subsonic M0.85", 0.85), ("Boomless M1.15", 1.15), ("Boomless M1.3", 1.3)]:
                val_pkm = dt_per_km(mref, mach) * VTTS["base"] * av
                be_pax = ann / val_pkm
                rows.append({"design": d, "case": c, "vs": ref, "capex_musd_km": cx["musd_km"],
                             "usd_m2": cx["usd_m2"], "opex_musd_km_yr": sum(ox.values()) / 1e6,
                             "annual_cost_musd_km_yr": ann / 1e6,
                             "min_saved_per_1000km": dt_per_km(mref, mach) * 60 * 1000,
                             "value_usd_per_pax_km": val_pkm,
                             "breakeven_pax_yr_vtts": be_pax,
                             "breakeven_flights_day_each_way": be_pax / (SEATS * LOAD_FACTOR * 365 * 2),
                             "breakeven_pax_yr_fullfare": ann / (FARE_ONEWAY / 4000.0 * av)})
    with open(HERE / "annual_breakeven.csv", "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        w.writeheader()
        for r in rows:
            w.writerow({k: (f"{x:.4g}" if isinstance(x, float) else x) for k, x in r.items()})

    # ---- corridors ----------------------------------------------------------------
    cx_b, ox_b = cap[(BASE_DESIGN, "base")], opx[(BASE_DESIGN, "base")]
    cx_17, ox_17 = cap[("M1.7 (Overture-class, lift-set)", "base")], opx[("M1.7 (Overture-class, lift-set)", "base")]
    ann_17 = annual_cost_km(cx_17, ox_17, "base")
    corr = []
    seg_best = None  # pax-h per corridor km per year (vs boomless, vs subsonic) on the densest trunk km
    if routes:
        data, m = routes
        av = v("availability", "base")
        for r in data:
            km = fnum(r.get("overland_km", r.get(m["km"])))
            shield = fnum(r.get("shield_km", "nan"))
            if math.isnan(shield):
                shield = km
            ph_sub = fnum(r.get("pax_hours_saved_vs_subsonic", r.get(m["ph_sub"]) if m["ph_sub"] else "nan"))
            ph_inc = fnum(r.get("pax_hours_saved_vs_boomless", r.get(m["ph_inc"]) if m["ph_inc"] else "nan"))
            seats = fnum(r.get("seats_equiv", "nan"))
            capf = fnum(r.get("capture", "nan"))
            lf = fnum(r.get("load_factor", "nan"))
            pax = seats * capf * lf if not any(math.isnan(x) for x in (seats, capf, lf)) else fnum(r.get(m["pax"]))
            if not shield or math.isnan(shield) or shield <= 0:
                continue
            # ramps under the acceleration focus need extra width (FOCUS plug-in)
            shield_eff = shield + 2 * FOCUS["ramp_km_per_end"] * (FOCUS["ramp_width_factor"] - 1)
            cost = ann_17 * shield_eff
            ben_inc = ph_inc * VTTS["base"] * av
            ben_sub = ph_sub * VTTS["base"] * av
            corr.append({"route": f"{r.get('origin', '')}-{r.get('dest', '')}", "rank": r.get("rank"),
                         "overland_km": km, "shield_km": shield,
                         "captured_pax_yr": pax, "capture": capf,
                         "pax_h_saved_vs_sub": ph_sub, "pax_h_saved_vs_boomless": ph_inc,
                         "benefit_vs_sub_musd_yr": ben_sub / 1e6, "benefit_vs_boomless_musd_yr": ben_inc / 1e6,
                         "benefit_vs_boomless_usd_per_shield_km_yr": ben_inc / shield,
                         "annual_cost_musd_yr": cost / 1e6,
                         "annual_cost_musd_per_shield_km_yr": cost / shield / 1e6,
                         "bcr_vs_boomless": ben_inc / cost, "bcr_vs_sub": ben_sub / cost,
                         "breakeven_capture_vs_boomless": capf * cost / ben_inc if ben_inc > 0 else float("nan"),
                         "cost_per_captured_trip_usd": cost / pax if pax and pax > 0 else float("nan")})
        corr.sort(key=lambda r: -r["bcr_vs_boomless"] if not math.isnan(r["bcr_vs_boomless"]) else 0)
        segf = ROUTES_CSV.parent / "corridor_segments.csv"
        if segf.exists():
            with open(segf, newline="") as f:
                segs = list(csv.DictReader(f))
            with open(HERE / "segment_economics.csv", "w", newline="") as f:
                w = csv.writer(f)
                w.writerow(["segment", "type", "km", "pax_h_per_km_vs_sub", "pax_h_per_km_vs_boomless",
                            "benefit_vs_boomless_usd_per_km_yr", "benefit_vs_sub_usd_per_km_yr",
                            "annual_cost_usd_per_km_yr_M17_base", "bcr_vs_boomless", "bcr_vs_sub"])
                for sg in segs:
                    bl = fnum(sg["pax_h_per_km_vs_bl"]); sb = fnum(sg["pax_h_per_km_vs_sub"])
                    w.writerow([f"{sg['from_iata']}-{sg['to_iata']}", sg["segment_type"], sg["distance_km"], sb, bl,
                                f"{bl * VTTS['base'] * av:.0f}", f"{sb * VTTS['base'] * av:.0f}", f"{ann_17:.0f}",
                                f"{bl * VTTS['base'] * av / ann_17:.2e}", f"{sb * VTTS['base'] * av / ann_17:.2e}"])
                    if seg_best is None or bl > seg_best[0]:
                        seg_best = (bl, sb, f"{sg['from_iata']}-{sg['to_iata']}")
            notes.append(f"corridor_segments.csv read; densest trunk km {seg_best[2]}: "
                         f"{seg_best[0]:.0f} pax-h/km/yr vs Boomless, {seg_best[1]:.0f} vs subsonic (10 % capture)")
    else:
        # placeholder: transcontinental trunk LAX-DEN-ORD-JFK (~4,000 km), M1.7 design.
        # Cost cases at base demand, plus one row with every assumption favourable.
        cap_pax = pax_from_flights(FLIGHTS_CAP) * 2
        scen = [(c, c, MARKET_PAX_PLACEHOLDER * CAPTURE["base"], VTTS["base"]) for c in CASES]
        scen.append(("all-favourable", "low", cap_pax, VTTS["high"]))
        for label, c, pax, vt in scen:
            cx, ox = cap[("M1.7 (Overture-class, lift-set)", c)], opx[("M1.7 (Overture-class, lift-set)", c)]
            ann = annual_cost_km(cx, ox, c)
            pax = min(pax, cap_pax)
            L = 4000.0
            cruise = L - 2 * FOCUS["ramp_km_per_end"]
            h_inc = dt_per_km(1.15, 1.7) * cruise
            h_sub = dt_per_km(0.85, 1.7) * cruise
            av = v("availability", c)
            ben_inc = pax * h_inc * vt * av
            ben_sub = pax * h_sub * vt * av
            cost = ann * (L + 2 * FOCUS["ramp_km_per_end"] * (FOCUS["ramp_width_factor"] - 1))
            corr.append({"route": f"PLACEHOLDER trunk LAX-JFK 4,000 km, cost {label}", "overland_km": L, "pax_yr": pax,
                         "pax_h_saved_vs_sub": pax * h_sub, "pax_h_saved_vs_boomless": pax * h_inc,
                         "benefit_vs_sub_musd_yr": ben_sub / 1e6, "benefit_vs_boomless_musd_yr": ben_inc / 1e6,
                         "benefit_vs_boomless_usd_per_km_yr": ben_inc / L,
                         "annual_cost_musd_yr": cost / 1e6, "annual_cost_usd_per_km_yr": cost / L,
                         "bcr_vs_boomless": ben_inc / cost, "bcr_vs_sub": ben_sub / cost,
                         "cost_per_trip_usd": cost / pax})
    if corr:
        with open(HERE / "corridor_economics.csv", "w", newline="") as f:
            w = csv.DictWriter(f, fieldnames=list(corr[0].keys()))
            w.writeheader()
            for r in corr:
                w.writerow({k: (f"{x:.4g}" if isinstance(x, float) else x) for k, x in r.items()})

    # ---- land ---------------------------------------------------------------------
    land = [land_rights(c) for c in CASES]
    with open(HERE / "land_rights.csv", "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(land[0].keys()))
        w.writeheader()
        for r in land:
            w.writerow({k: (f"{x:.4g}" if isinstance(x, float) else x) for k, x in r.items()})
    shade = shadow_hours()

    # ---- tornado: BCR of M1.7 shield vs Boomless M1.15, per corridor km ------------
    def bcr(usd_m2_scale=1.0, width=None, life=None, r=None, n=None, pax=None, vtts=None, design="M1.7 (Overture-class, lift-set)"):
        cx = capex(design, "base", width_km=width, usd_m2_scale=usd_m2_scale)
        ox = opex(cx, "base", film_life=life)
        ann = annual_cost_km(cx, ox, "base", r=r, n=n)
        vt = vtts if vtts is not None else VTTS["base"]
        if seg_best is not None:   # routes plug-in: pax-h per km on the densest trunk km at 10 % capture
            ph = seg_best[0] * (pax / 0.10 if pax is not None else 1.0)
            return ph * vt * v("availability", "base") / ann
        p = MARKET_PAX_PLACEHOLDER * (pax if pax is not None else CAPTURE["base"])
        return p * dt_per_km(1.15, 1.7) * vt * v("availability", "base") / ann

    base_bcr = bcr()
    sens = [
        ("All-in cost per m² (×0.25 / ×4)", bcr(usd_m2_scale=0.25), bcr(usd_m2_scale=4.0)),
        ("Width (1.25 / 2.0 km)", bcr(width=1.25), bcr(width=2.0)),
        ("Film life (20 / 5 yr)", bcr(life=20), bcr(life=5)),
        ("Discount rate (5 / 8 %)", bcr(r=0.05), bcr(r=0.08)),
        ("Asset life (30 / 20 yr)", bcr(n=30), bcr(n=20)),
        ("Capture fraction (30 / 5 %)", bcr(pax=0.30), bcr(pax=0.05)),
        ("Value of time (\\$200 / \\$49 per h)", bcr(vtts=200), bcr(vtts=VTTS["low"])),
    ]
    with open(HERE / "tornado.csv", "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["parameter", "bcr_favourable", "bcr_unfavourable", "bcr_base"])
        for name, a, b in sens:
            w.writerow([name, f"{a:.3e}", f"{b:.3e}", f"{base_bcr:.3e}"])

    # ---- required $/m2 for break-even vs throughput ---------------------------------
    pax_grid = [1e5 * 10 ** (i / 10) for i in range(0, 41)]  # 1e5 .. 1e9
    def req_usd_m2(pax, mref, design="M1.7 (Overture-class, lift-set)"):
        # Solve scale s so that annual cost = benefit; annual cost ~ linear in s except fixed opex
        cx1 = capex(design, "base", usd_m2_scale=1.0)
        a1 = annual_cost_km(cx1, opex(cx1, "base"), "base")
        cx0 = capex(design, "base", usd_m2_scale=1e-9)
        a0 = annual_cost_km(cx0, opex(cx0, "base"), "base")
        ben = pax * dt_per_km(mref, 1.7) * VTTS["base"] * v("availability", "base")
        if ben <= a0:
            return float("nan")
        s = (ben - a0) / (a1 - a0)
        return s * cx1["usd_m2"]
    curve = [(p, req_usd_m2(p, 0.85), req_usd_m2(p, 1.15)) for p in pax_grid]
    with open(HERE / "breakeven_curve.csv", "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["pax_per_yr_through_km", "req_allin_usd_m2_vs_subsonic", "req_allin_usd_m2_vs_boomless_M1.15"])
        for p, a, b in curve:
            w.writerow([f"{p:.4g}", f"{a:.4g}", f"{b:.4g}"])

    make_figure(cap, sens, base_bcr, curve)

    # ---- console summary ------------------------------------------------------------
    print("\n".join(notes))
    for d in DESIGNS:
        for c in CASES:
            cx = cap[(d, c)]
            print(f"{d:40s} {c:5s} W={cx['W_km']:.2f} km  ${cx['usd_m2']:7.1f}/m2  ${cx['musd_km']:7.1f}M/km  "
                  f"opex ${sum(opx[(d, c)].values()) / 1e6:6.1f}M/km/yr  annual ${annual_cost_km(cx, opx[(d, c)], c) / 1e6:6.1f}M/km/yr  "
                  f"H2 {cx['h2_t_km']:.0f} t/km  lid {cx['lid_kg_m2']:.2f} kg/m2 = {cx['lid_t_km']:.0f} t/km  tether {cx['tether_t_km']:.0f} t/km")
    for c in CASES:
        print(f"faired tethers {c}: ${fair[c]['usd_m2']:.1f}/m2 vs bare ${cap[(BASE_DESIGN, c)]['usd_m2']:.1f}/m2")
    print("base BCR (M1.7 vs Boomless M1.15, densest segment or placeholder):", f"{base_bcr:.2e}")
    print("shadow (opaque, 37N, E-W): max h/yr, max h/day, at km north:", [round(x, 1) for x in shade])
    for r in land:
        print(r)
    for r in corr:
        print(r)


def make_figure(cap, sens, base_bcr, curve):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    BLUE, ORANGE, AQUA, YELLOW, MAG, GREEN, VIOLET, RED = (
        "#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4", "#008300", "#6250d6", "#e34948")
    INK, INK2, GRID = "#0b0b0b", "#52514e", "#e4e3df"
    plt.rcParams.update({"font.size": 9, "axes.edgecolor": INK2, "axes.labelcolor": INK, "xtick.color": INK2,
                         "ytick.color": INK2, "axes.spines.top": False, "axes.spines.right": False})
    fig, axs = plt.subplots(1, 3, figsize=(15, 4.8), gridspec_kw={"width_ratios": [1.0, 1.25, 1.1]})
    fig.patch.set_facecolor("#fcfcfb")

    # (a) stacked capex by group, M2 baseline
    groups = {
        "Lid and envelope films": ["Lid mass (PE sheet/laminate)", "Lid UV/weathering skin", "Load tapes / scrim (UHMWPE)",
                                   "Floor film", "Partitions and bellows", "Soft ties"],
        "Seams, fabrication, install": ["Seams and fabrication", "Launch / installation"],
        "H₂ first fill": ["H2 first fill"],
        "Tethers, anchors, stations": ["Tethers", "Anchors", "Winch/tension stations"],
        "Controls, shuttle, land": ["Controls, sensors, transponders", "Diurnal gas shuttle", "Land and easements"],
        "Engineering + contingency": ["Engineering / PM", "Contingency"],
    }
    colors = [BLUE, ORANGE, AQUA, YELLOW, VIOLET, INK2]
    ax = axs[0]
    labels = []
    for j, d in enumerate(DESIGNS):
        for i, c in enumerate(CASES):
            x = j * 3.6 + i
            labels.append((x, f"{c}\n{'M2.0' if j == 0 else 'M1.7'}"))
            bottom = 0
            L = cap[(d, c)]["lines_usd_m2"]
            for (g, ks), col in zip(groups.items(), colors):
                h = sum(L[k] for k in ks)
                ax.bar(x, h, bottom=bottom, color=col, width=0.8, edgecolor="#fcfcfb", linewidth=1.0,
                       label=g if (i == 0 and j == 0) else None)
                bottom += h
            ax.text(x, bottom * 1.02, f"\\${bottom:.0f}", ha="center", va="bottom", fontsize=8, color=INK)
    ax.axhline(10, color=RED, lw=1.2, ls="--")
    ax.set_xlim(-0.6, 7.0)
    ax.text(6.1, 13, "\\$10/m²\nhandoff", color=INK, fontsize=7, ha="left")
    ax.set_xticks([x for x, _ in labels])
    ax.set_xticklabels([t for _, t in labels], fontsize=8)
    ax.set_ylabel("All-in capital cost, $ per m² of planform")
    ax.set_title("(a) Capital cost by case", loc="left", color=INK, fontsize=10)
    ax.grid(axis="y", color=GRID, lw=0.6)
    ax.set_axisbelow(True)
    ax.legend(fontsize=7, frameon=False, loc="upper center", bbox_to_anchor=(0.5, -0.13), ncol=2)

    # (b) tornado on log BCR
    ax = axs[1]
    order = sorted(sens, key=lambda s: abs(math.log10(s[1]) - math.log10(s[2])))
    for i, (name, a, b) in enumerate(order):
        ax.barh(i, a - base_bcr, left=base_bcr, color=AQUA, height=0.6)
        ax.barh(i, b - base_bcr, left=base_bcr, color=ORANGE, height=0.6)
    ax.set_xscale("log")
    ax.axvline(base_bcr, color=INK, lw=1)
    ax.axvline(1.0, color=RED, lw=1.2, ls="--")
    ax.text(0.9, len(order) - 0.6, "break-even", color=INK, fontsize=8, ha="right")
    ax.set_yticks(range(len(order)))
    ax.set_yticklabels([o[0] for o in order], fontsize=8)
    ax.set_xlabel("Benefit/cost ratio, M1.7 shield vs Boomless M1.15 (log)")
    ax.set_title(f"(b) Sensitivity around base BCR {base_bcr:.1e}", loc="left", color=INK, fontsize=10)
    ax.grid(axis="x", color=GRID, lw=0.6, which="both")
    ax.set_axisbelow(True)
    ax.set_xlim(1e-4, 2)
    from matplotlib.patches import Patch
    ax.legend(handles=[Patch(color=AQUA, label="favourable end"), Patch(color=ORANGE, label="unfavourable end")],
              fontsize=7, frameon=False, loc="lower right")

    # (c) required all-in $/m2 vs throughput
    ax = axs[2]
    p = [c[0] for c in curve]
    ax.plot(p, [c[1] for c in curve], color=BLUE, lw=2, label="vs subsonic M0.85")
    ax.plot(p, [c[2] for c in curve], color=ORANGE, lw=2, label="vs Boomless M1.15")
    ax.set_xscale("log")
    ax.set_yscale("log")
    base_m2 = cap[("M1.7 (Overture-class, lift-set)", "base")]["usd_m2"]
    low_m2 = cap[("M1.7 (Overture-class, lift-set)", "low")]["usd_m2"]
    ax.axhline(base_m2, color=INK2, lw=1, ls=":")
    ax.text(1.2e5, base_m2 * 1.1, f"base estimate \\${base_m2:.0f}/m²", fontsize=8, color=INK)
    ax.axhline(low_m2, color=INK2, lw=1, ls=":")
    ax.text(1.2e5, low_m2 * 1.1, f"low estimate \\${low_m2:.0f}/m²", fontsize=8, color=INK)
    ax.axhline(10, color=RED, lw=1.2, ls="--")
    ax.text(1.2e5, 11, "\\$10/m²", fontsize=8, color=INK)
    cap_pax = pax_from_flights(FLIGHTS_CAP) * 2
    ax.axvline(cap_pax, color=GREEN, lw=1.2)
    ax.text(cap_pax * 0.9, 400, "lid-fatigue cap\n~100 flights/day\neach way", fontsize=7, color=INK, ha="right", va="top")
    ax.axvline(1.7e8, color=VIOLET, lw=1.2)
    ax.text(1.7e8 * 0.9, 400, "Tokaido\nShinkansen\n~170 M/yr [U]", fontsize=7, color=INK, ha="right", va="top")
    ax.set_ylim(0.1, 1000)
    ax.set_xlim(1e5, 1e9)
    ax.set_xlabel("Passengers per year crossing each corridor km")
    ax.set_ylabel("All-in \\$/m² at which BCR = 1 (VTTS \\$86/h)")
    ax.set_title("(c) Break-even cost vs traffic, M1.7 design", loc="left", color=INK, fontsize=10)
    ax.grid(color=GRID, lw=0.6, which="major")
    ax.legend(fontsize=7, frameon=False, loc="lower right")
    fig.tight_layout()
    fig.savefig(HERE / "economics.png", dpi=140, facecolor=fig.get_facecolor())


if __name__ == "__main__":
    main()
