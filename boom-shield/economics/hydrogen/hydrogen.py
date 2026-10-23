#!/usr/bin/env python3
"""Cheap-hydrogen scenarios for BoomShield, Skyloop and Hyperloop (screening level).

Run:  python3 -I boom-shield/economics/hydrogen/hydrogen.py      (light CPU, ~15 s; reads caches, re-runs nothing heavy)

Scenarios (company claims are UNVERIFIED search summaries; see assumptions.csv):
  Baseline        current repo: shield H2 $1.6/kg grey SMR + $0.5 delivery (integration), electrolytic $5-7/kg
                  (growth T6), canonical fuel $0.80/kg + CO2 $190/t, ground PV $28/MWh (cost-down aggressive)
  H2 $2/kg        Terraform demo today (< $2.50/kg reported Mar 2024, "< $2" claimed); DAC $250/t; PV $17/MWh
  H2 $1/kg        Terraform plan / DOE Hydrogen Shot target; DAC $100/t; PV $10/MWh (Terraform whitepaper)

Channels traced:
  1. direct H2 cost: shield fill + make-up (integration plant model re-run), Skyloop envelope gas
  1b. cheap solar erodes the shield's PV dual-use credit (cost-down market cap = ground PV LCOE)
  2. e-SAF bottom-up (H2 + DAC CO2 + conversion) and a carbon-neutral fuel convention: fuel at the e-SAF
     price, no CO2 charge on that fuel (non-CO2 effects uncounted, as in the canonical convention)
  3. firmed solar electricity for the tubes' traction (24 h service needs storage, not daytime prices)

Writes next to this file: esaf_cost.csv, electricity.csv, direct_h2.csv, boomshield_net.csv, tubes_net.csv,
breakeven.csv, assumptions.csv.

Imports read-only (nothing is written outside this folder; module state is patched in memory and restored):
economics/net/net_bc.py (integrated_cost, DESIGNS, canonical constants; through it highmach, integration,
cost-down, cost_model, growth_model), net/canonical.csv + costs.csv + lifecycle.csv, hyperloop/hyperloop.py
(bc_rows; through it skyloop-review/system/sky_model.py lifecycle, capex_level, run_cost), skyloop-review/
system/networks.csv + lifecycle.csv, hyperloop/hl_points.csv + hl_econ.csv.
[U] = engineering assumption.
"""
from __future__ import annotations

import contextlib
import importlib.util
import sys
from pathlib import Path

sys.dont_write_bytecode = True
import numpy as np   # noqa: E402
import pandas as pd  # noqa: E402

HERE = Path(__file__).resolve().parent
ECON = HERE.parent
BS = ECON.parent
ROOT = BS.parent


def _load(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    mod = importlib.util.module_from_spec(spec)
    sys.modules[name] = mod
    spec.loader.exec_module(mod)
    return mod


nb = _load("net_bc", ECON / "net" / "net_bc.py")          # canonical convention (gm.energy gains 'flat' in memory)
gm, ig = nb.gm, nb.ig
cd, cm = ig.cd, ig.cm
hl = _load("hyperloop", ROOT / "hyperloop" / "hyperloop.py")   # also loads sky_model (hl.sm)
sm = hl.sm

FUEL0, CO2_KG, SCC0 = nb.FUEL_USD_KG, nb.CO2_KG, nb.SCC
ALLIN0 = FUEL0 + CO2_KG * SCC0 / 1e3          # $1.40/kg: canonical fossil fuel + CO2 charge
CASES = ("low", "central", "high")
SCEN = ["Baseline (repo)", "H2 $2/kg (Terraform demo)", "H2 $1/kg (Terraform/DOE target)"]

# ------------------------------------------------------------------------------------------- scenario inputs
# At-plant H2 for e-SAF (low/central/high), DAC CO2 $/t, shield H2 at plant, ground PV LCOE $/MWh.
S = {
    SCEN[0]: dict(h2=(5.0, 6.0, 7.0), dac=(250.0, 500.0, 1000.0), shield_h2=None, pv=None, solar=None),
    SCEN[1]: dict(h2=(1.5, 2.0, 2.5), dac=(150.0, 250.0, 400.0), shield_h2=2.0, pv=17.0, solar=17.0),
    SCEN[2]: dict(h2=(0.8, 1.0, 1.5), dac=(50.0, 100.0, 200.0), shield_h2=1.0, pv=10.0, solar=10.0),
}
SHIELD_DELIV = (0.3, 0.5, 1.5)   # $/kg compression, pipeline, lift to 20 km: integration low/central, cost_model base

# e-SAF (power-to-liquid jet via RWGS + Fischer-Tropsch or methanol-to-jet).  Both routes share the
# stoichiometric floor CO2 + 3 H2 -> (CH2) + 2 H2O: 0.432 kg H2 and 3.14 kg CO2 per kg of liquid product.
H2_STOICH, CO2_STOICH = 3 * 2.016 / 14.027, 44.01 / 14.027
LHV_KWH_KG = 43.0 / 3.6
ESAF = {  # low-cost / central / high-cost  [U]
    "eta_h2": (0.85, 0.75, 0.65),          # H2 to liquid product (light gases recycled or burnt)
    "eta_c": (0.92, 0.82, 0.72),           # carbon to liquid product
    "buffer_usd_kg_h2": (0.15, 0.40, 0.80),  # storage + compression to run synthesis steadily on solar H2
    "capex_usd_kw": (800.0, 1500.0, 2800.0),  # RWGS/FT or MeOH/MTJ + upgrading, per kW of fuel (LHV)
    "cf": (0.90, 0.85, 0.70),              # synthesis capacity factor
    "om_frac": (0.03, 0.04, 0.05),         # fixed O&M per yr of capex
    "utilities_usd_kg": (0.03, 0.07, 0.15),  # electricity/heat for RWGS, compression, catalyst
    "distribution_usd_kg": (0.03, 0.05, 0.08),  # pipeline/truck to the airport
}
R, LIFE = 0.07, 25

# Firmed solar for 24 h tube service [U]
FIRM = {"overbuild": (1.2, 1.3, 1.6),     # extra PV for winter and cloudy days
        "shift": (0.5, 0.6, 0.7),         # share of traction energy drawn outside solar hours (18 h service)
        "batt_usd_kwh": (120.0, 180.0, 250.0), "batt_life": (20, 15, 12), "batt_om": (0.015, 0.02, 0.025),
        "rte": (0.90, 0.88, 0.85)}

CV = ("H2 only (PV credit kept)", "H2 + cheap-solar PV erosion")
SC_CV = [(SCEN[0], CV[0])] + [(sc, cv) for sc in SCEN[1:] for cv in CV]
DES = ["M1.7 Overture (62 t)", "M2.0 adaptive (49 t)", "M2.2 design A (50 t)"]
OPEN = (2035, 2045, 2050)


def crf(r, n):
    return r * (1 + r) ** n / ((1 + r) ** n - 1)


# ------------------------------------------------------------------------------------------- channel 2: e-SAF
def esaf_price(h2, dac, i):
    """$/kg e-SAF delivered to the airport: H2 at plant h2 $/kg, DAC CO2 dac $/t, parameter set i (0/1/2)."""
    e = {k: v[i] for k, v in ESAF.items()}
    kg_h2 = H2_STOICH / e["eta_h2"]
    kg_co2 = CO2_STOICH / e["eta_c"]
    kg_yr_per_kw = 8760 * e["cf"] / LHV_KWH_KG
    conv = e["capex_usd_kw"] / kg_yr_per_kw * (crf(R, LIFE) + e["om_frac"])
    parts = dict(kg_h2_per_kg=kg_h2, kg_co2_per_kg=kg_co2, h2_usd=kg_h2 * (h2 + e["buffer_usd_kg_h2"]),
                 co2_usd=kg_co2 * dac / 1e3, conversion_usd=conv, utilities_usd=e["utilities_usd_kg"],
                 distribution_usd=e["distribution_usd_kg"])
    parts["esaf_usd_kg"] = sum(parts[k] for k in ("h2_usd", "co2_usd", "conversion_usd", "utilities_usd",
                                                  "distribution_usd"))
    parts["abatement_usd_t"] = (parts["esaf_usd_kg"] - FUEL0) / CO2_KG * 1e3
    return parts


def esaf_table():
    rows = []
    for sc in SCEN:
        for i, case in enumerate(CASES):
            h2, dac = S[sc]["h2"][i], S[sc]["dac"][i]
            rows.append(dict(scenario=sc, case=case, h2_usd_kg=h2, dac_usd_t=dac, **esaf_price(h2, dac, i)))
    # H2 sweep at the central parameter set and $250 / $100 DAC (for the break-even read-across)
    for dac in (250.0, 100.0):
        for h2 in (0.0, 0.5, 1.0, 1.5, 2.0, 3.0, 4.0, 6.0):
            rows.append(dict(scenario=f"sweep, DAC ${dac:.0f}/t", case="central", h2_usd_kg=h2, dac_usd_t=dac,
                             **esaf_price(h2, dac, 1)))
    rows.append(dict(scenario="floor: free H2, DAC $50/t", case="low", h2_usd_kg=0.0, dac_usd_t=50.0, **esaf_price(0.0, 50.0, 0)))
    return pd.DataFrame(rows)


def esaf_h2_for(price, dac, i=1):
    """At-plant H2 $/kg giving e-SAF = price (parameter set i, DAC dac; linear in h2).  NaN if below the
    free-H2 floor."""
    p0, p1 = esaf_price(0.0, dac, i)["esaf_usd_kg"], esaf_price(1.0, dac, i)["esaf_usd_kg"]
    h = (price - p0) / (p1 - p0)
    return h if h >= 0 else np.nan


# ------------------------------------------------------------------------------------------- channel 3: power
def firmed(solar, i):
    f = {k: v[i] for k, v in FIRM.items()}
    stor = f["batt_usd_kwh"] * (crf(R, f["batt_life"]) + f["batt_om"]) / (365 * 0.9) * 1e3   # $/MWh shifted
    return solar * f["overbuild"] + f["shift"] * (stor + solar * (1 / f["rte"] - 1)), stor


def electricity_table():
    rows = []
    for sc in SCEN[1:]:
        sol = S[sc]["solar"]
        for i, case in enumerate(CASES):
            p, stor = firmed(sol, i)
            rows.append(dict(scenario=sc, case=case, daytime_solar_usd_mwh=sol, overbuild=FIRM["overbuild"][i],
                             shift_share=FIRM["shift"][i], storage_usd_per_mwh_shifted=stor, firmed_24h_usd_mwh=p))
    for y in (2035, 2050, 2060):
        _, scc, _, _ = gm.energy("central", y)
        rows.append(dict(scenario=f"repo grid {y} (sky_model)", case="central", firmed_24h_usd_mwh=sm.P["elec_usd_kwh"] * 1e3,
                         grid_co2_charge_usd_mwh=sm.grid_ci("central", y) * scc,
                         all_in_usd_mwh=sm.P["elec_usd_kwh"] * 1e3 + sm.grid_ci("central", y) * scc))
    return pd.DataFrame(rows)


# ------------------------------------------------------------------------------------------- channel 1: shield cost
@contextlib.contextmanager
def shield_prices(h2=None, deliv=None, pv=None, pv_capex_follows=True):
    """In-memory override of integration's H2 price/delivery and cost-down's aggressive PV market cap."""
    sP, sPV = dict(ig.P), dict(cd.PV_AGGR)
    try:
        if h2 is not None:
            ig.P["h2_usd_kg"] = (h2, h2, h2)
        if deliv is not None:
            ig.P["h2_deliv_usd_kg"] = (deliv, deliv, deliv)
        if pv is not None:
            k = pv / sPV["ground_lcoe"]
            cd.PV_AGGR["ground_lcoe"] = pv
            if pv_capex_follows:     # same-technology world: shield modules and BOS get cheaper too
                cd.PV_AGGR["module_usd_w"] = sPV["module_usd_w"] * k
                cd.PV_AGGR["bos_usd_w"] = sPV["bos_usd_w"] * k
        yield
    finally:
        ig.P.clear(); ig.P.update(sP)
        cd.PV_AGGR.clear(); cd.PV_AGGR.update(sPV)


HD = {d["design"]: d for d in nb.hn.mach_designs()}


def shield_cost(lab, **kw):
    ds = HD[nb.DESIGNS[lab]]
    with shield_prices(**kw):
        return nb.integrated_cost(ds["T_eq20_m"], ds["width_km"], f"H2 {lab}")


def direct_table():
    rows, costs = [], {}
    pub = pd.read_csv(ECON / "net" / "costs.csv", index_col=0)
    deliv_c = SHIELD_DELIV[1]
    for lab in DES:
        ds = HD[nb.DESIGNS[lab]]
        base = shield_cost(lab)
        assert abs(base["net"] - pub.loc[lab, "net"]) < 2e-3, (lab, base["net"], pub.loc[lab, "net"])
        zero = shield_cost(lab, h2=0.0, deliv=0.0)
        h2_t_km = cm.RHO_H2_20KM * ds["T_eq20_m"] * ds["width_km"] * 1e3 * cm.v("h2_fill_loss", "base")
        mk = ig.lever_set("aggressive")["h2_makeup"]
        for sc in SCEN:
            p = S[sc]
            variants = [("H2 + PV (scenario)", dict(h2=p["shield_h2"], deliv=deliv_c if p["shield_h2"] else None, pv=p["pv"])),
                        ("H2 only", dict(h2=p["shield_h2"], deliv=deliv_c if p["shield_h2"] else None)),
                        ("PV only, shield PV capex follows", dict(pv=p["pv"])),
                        ("PV only, shield PV capex fixed", dict(pv=p["pv"], pv_capex_follows=False))]
            if p["shield_h2"]:
                variants += [(f"H2 only, delivery ${d}/kg", dict(h2=p["shield_h2"], deliv=d)) for d in SHIELD_DELIV if d != deliv_c]
            for var, kw in variants:
                c = shield_cost(lab, **kw) if (sc != SCEN[0]) else base
                h2_deliv = (kw.get("h2") if kw.get("h2") is not None else ig.P["h2_usd_kg"][1]) + \
                    (kw.get("deliv") if kw.get("deliv") is not None else ig.P["h2_deliv_usd_kg"][1])
                rows.append(dict(design=lab, scenario=sc, variant=var, h2_delivered_usd_kg=h2_deliv,
                                 ground_pv_usd_mwh=kw.get("pv") or cd.PV_AGGR["ground_lcoe"],
                                 h2_t_km=h2_t_km, makeup_frac_yr=mk,
                                 capex_musd_km=c["capex"], opex_musd_km_yr=c["opex"], dual_musd_km_yr=c["dual"],
                                 net_musd_km_yr=c["net"], d_net_vs_baseline=c["net"] - base["net"],
                                 h2_share_capex_baseline=(base["capex"] - zero["capex"]) / base["capex"],
                                 h2_share_net_cost_baseline=(base["net"] - zero["net"]) / base["net"]))
                if var == "H2 + PV (scenario)":
                    costs[(lab, sc, CV[1])] = c
                if var == "H2 only" or sc == SCEN[0]:
                    costs[(lab, sc, CV[0])] = c
                if sc == SCEN[0]:
                    costs[(lab, sc, CV[1])] = c
                    break
    # cost_model base case (green $4 + $1.5 delivered) for the growth.md "< 4 % of capex" claim
    cx = cm.capex(cm.BASE_DESIGN, "base")
    L = cx["lines_usd_m2"]
    mark = (1 + cm.v("engineering_frac", "base")) * (1 + cm.v("contingency_frac", "base"))
    rows.append(dict(design=cm.BASE_DESIGN, scenario="cost_model base ($4 + $1.5 delivered)", variant="share check",
                     h2_delivered_usd_kg=cm.v("h2_usd_per_kg", "base") + cm.v("h2_delivery_usd_per_kg", "base"),
                     h2_t_km=cx["h2_t_km"] * cm.v("h2_fill_loss", "base"), capex_musd_km=cx["musd_km"],
                     h2_share_capex_baseline=L["H2 first fill"] * mark / cx["usd_m2"]))
    # Skyloop floating tube: envelope gas (physics R 43 / 82 / 88 m) against the run cost
    for i, (lvl, Rm) in enumerate(zip(("low", "base", "high"), (43.0, 82.0, 88.0))):
        t_km = np.pi * Rm ** 2 * 1e3 * cm.RHO_H2_20KM / 1e3 * cm.v("h2_fill_loss", "base")
        run = sm.run_cost(lvl)
        mark = (1 + sm.SOFT[i]) * (1 + sm.CONTING[i])
        for sc, price in ((SCEN[0] + " (cost_model base $5.5 delivered)", 5.5),
                          (SCEN[1], 2.0 + deliv_c), (SCEN[2], 1.0 + deliv_c)):
            rows.append(dict(design=f"Skyloop envelope R {Rm:.0f} m ({lvl})", scenario=sc, variant="gas fill",
                             h2_delivered_usd_kg=price, h2_t_km=t_km, capex_musd_km=run,
                             gas_musd_km=t_km * price * mark / 1e3, gas_share_run_capex=t_km * price * mark / 1e3 / run,
                             makeup_10pct_musd_km_yr=0.1 * t_km * price / 1e3,
                             om_musd_km_yr=run * sm.OM_FRAC[i]))
    return pd.DataFrame(rows), costs


# ------------------------------------------------------------------------------------------- energy paths
_ENERGY0 = gm.energy


def _make_energy(base_energy):
    def energy(case, y):
        """'esaf:P': all fuel at P $/kg (real, flat), no CO2 charge on fuel (co2 factor 0), SCC path kept for the
        grid, central efficiency gain.  Anything else -> the original path."""
        if isinstance(case, str) and case.startswith("esaf:"):
            _, scc, _, eff = base_energy("central", y)
            return float(case[5:]), scc, 0.0, eff
        return base_energy(case, y)
    return energy


gm.energy = _make_energy(_ENERGY0)
sm.gm.energy = _make_energy(sm.gm.energy) if sm.gm is not gm else gm.energy


# ------------------------------------------------------------------------------------------- BoomShield net
def bs_inputs():
    T = pd.read_csv(ECON / "net" / "canonical.csv")
    nbm = pd.read_csv(BS / "network-maps" / "highmach" / "network_by_mach.csv").set_index("design")
    out = {}
    for lab in DES:
        for vs in ("bl", "sub"):
            for cap in (0.3, 1.0):
                r = T[(T.design == lab) & (T.vs == vs) & (T.demand == "base") & (T.capture == cap)
                      & (T.cost == "integrated")].iloc[0]
                kg = r.fuel_co2_busd / ALLIN0                      # bn kg extra fuel (captured + induced)
                out[(lab, vs, cap)] = dict(gross=r.gross_busd, kg=kg, co2=r.extra_co2_Mt,
                                           km=r.cost_busd_yr / r.cost_musd_km_yr * 1e3, C=r.cost_busd_yr,
                                           bc_net_pub=r.bc_net,
                                           flights=float(nbm.loc[nb.DESIGNS[lab], "peak_flights_day"]) * 4716 / 4500)
    return out


def bs_lifecycle(f1, cost, oy, energy, cap):
    """growth_model life-cycle B/C (base world, integrated cost) for flows f1 at 100 % capture."""
    net = dict(benefit_busd=f1["gross"], extra_fuel_busd=f1["kg"] * FUEL0, extra_co2_Mt=f1["co2"], km=f1["km"],
               cost_scale=1.0, flights=f1["flights"], build=10)
    c = {"integrated": dict(capex=cost["capex"], opex=cost["opex"], dual=cost["dual"], net=cost["net"])}
    return gm.metrics(gm.evaluate(c, net, oy, "base", energy, "integrated", cap), oy)


def boomshield(costs, E):
    F = bs_inputs()
    LCpub = pd.read_csv(ECON / "net" / "lifecycle.csv")
    rows, be = [], []
    esaf_c = {sc: E[(E.scenario == sc)].set_index("case").esaf_usd_kg for sc in SCEN}
    for lab in DES:
        for vs in ("bl", "sub"):
            for cap in (0.3, 1.0):
                f, f1 = F[(lab, vs, cap)], F[(lab, vs, 1.0)]
                for sc, cv in SC_CV:
                    c = costs[(lab, sc, cv)]
                    C = c["net"] * f["km"] / 1e3
                    fuels = [("canonical fossil $0.80 + $190/t", None, ALLIN0)] + \
                            [(f"e-SAF {k}", esaf_c[sc][k], esaf_c[sc][k]) for k in CASES]
                    for conv, P, allin in fuels:
                        r = dict(design=lab, vs=vs, capture=cap, scenario=sc, cost_variant=cv, fuel_convention=conv,
                                 fuel_all_in_usd_kg=allin, cost_musd_km_yr=c["net"], basis="2025 year",
                                 bc_gross=f["gross"] / C, bc_net=(f["gross"] - f["kg"] * allin) / C)
                        if sc == SCEN[0] and P is None:
                            assert abs(r["bc_net"] - f["bc_net_pub"]) < 2e-3, (lab, vs, cap, r["bc_net"], f["bc_net_pub"])
                        rows.append(r)
                        for oy in OPEN:
                            m = bs_lifecycle(f1, c, oy, "central" if P is None else f"esaf:{P}", cap)
                            if sc == SCEN[0] and P is None:
                                pub = LCpub[(LCpub.design == lab) & (LCpub.vs == vs) & (LCpub.capture == cap)
                                            & (LCpub.open_year == oy)].bc_net_lifecycle.iloc[0]
                                assert abs(m["bc_net"] - pub) < 2e-3, (lab, vs, cap, oy, m["bc_net"], pub)
                            rows.append(dict(r, basis=f"life-cycle {oy} opening",
                                             fuel_convention=conv if P is not None else "canonical (growth central path)",
                                             fuel_all_in_usd_kg=allin if P is not None else np.nan,
                                             bc_gross=m["bc_gross"], bc_net=m["bc_net"]))
                    # break-even e-SAF prices (linear in P)
                    tgt = dict(bl=1.0, sub=0.0)[vs]
                    sy0 = f["gross"] / C
                    slope = f["kg"] / C
                    pe = (sy0 - tgt) / slope if sy0 > tgt else np.nan
                    be.append(dict(concept="BoomShield", case=lab, vs=vs, capture=cap, scenario_costs=f"{sc}; {cv}",
                                   basis="2025 year", target=f"net B/C = {tgt:g}", esaf_breakeven_usd_kg=pe,
                                   net_at_free_fuel=sy0, note="" if sy0 > tgt else "unreachable: gross below target"))
                    for oy in OPEN:
                        m0 = bs_lifecycle(f1, c, oy, "esaf:0.0", cap)["bc_net"]
                        m1 = bs_lifecycle(f1, c, oy, "esaf:1.0", cap)["bc_net"]
                        pe = (m0 - tgt) / (m0 - m1) if m0 > tgt else np.nan
                        be.append(dict(concept="BoomShield", case=lab, vs=vs, capture=cap, scenario_costs=f"{sc}; {cv}",
                                       basis=f"life-cycle {oy} opening", target=f"net B/C = {tgt:g}",
                                       esaf_breakeven_usd_kg=pe, net_at_free_fuel=m0,
                                       note="" if m0 > tgt else "unreachable: gross below target"))
    return pd.DataFrame(rows), be


# ------------------------------------------------------------------------------------------- tubes
@contextlib.contextmanager
def power(usd_mwh=None, clean=False):
    sP, sg = sm.P["elec_usd_kwh"], sm.grid_ci
    try:
        if usd_mwh is not None:
            sm.P["elec_usd_kwh"] = usd_mwh / 1e3
        if clean:
            sm.grid_ci = lambda case, y: 0.0
        yield
    finally:
        sm.P["elec_usd_kwh"] = sP
        sm.grid_ci = sg


def skyloop_rows(energy, oy=2050):
    """sky_model.stage_lifecycle for spine+LAX all-stop, M4 (system.md headline), low/base cost."""
    N = pd.read_csv(ROOT / "skyloop-review" / "system" / "networks.csv")
    out = []
    for r in N[(N.network == "spine+LAX all-stop") & (N.M == 4.0) & N.vs.isin(["sub", "bl"])].itertuples():
        ann = dict(benefit_busd=r.benefit_busd, fuel_saved_Mt=r.fuel_saved_Mt, kwh_TWh=r.kwh_TWh,
                   pax_M=r.captured_pax_M + r.net_induced_M, pax_km_B=r.pax_km_B, peak=r.peak_seg_pax_M * 1e6)
        build = 10 if r.km > 2000 else 7
        cap_pax = r.bores_movblock * r.cap_pair_movblock_M * 1e6
        for capture in (0.3, 1.0):
            for lvl in ("low", "base"):
                lc = sm.lifecycle(ann, sm.capex_level(r, lvl), sm.OM_FRAC[sm.LVL[lvl]], oy, capture, build=build,
                                  cap_pax=cap_pax, energy_case=energy)
                out.append(dict(concept="Skyloop floating tube", case="spine+LAX all-stop M4", vs=r.vs,
                                capture=capture, cost=lvl, opening=oy, bc_gross=lc["bc_gross"], bc_net=lc["bc_net"],
                                pv_ben=lc["pv_ben"], pv_cost=lc["pv_cost"], pv_energy=lc["pv_energy"]))
    return out


def hyperloop_rows(energy, oy=2050):
    """hyperloop.bc_rows for spine+LAX all-stop, elevated full speed, M3 / M4 (README headline)."""
    P = pd.read_csv(ROOT / "hyperloop" / "hl_points.csv")
    out = []
    sel = P[(P.network == "spine+LAX all-stop") & (P.setting == "elevated") & (P.strategy == "full speed")
            & P.Mach.isin([3.0, 4.0])]
    old = sm.gm.energy
    for o in sel.to_dict("records"):
        capex = {lvl: dict(total=o[f"capex_total_{lvl}_busd"], tube=o[f"capex_tube_{lvl}_busd"]) for lvl in hl.LV}
        # bc_rows reads sm.gm.energy(energy_case='central'); route the 'central' call to the requested case
        sm.gm.energy = (lambda e: (lambda case, y: old(e, y)))(energy)
        try:
            rs = hl.bc_rows(o, capex)
        finally:
            sm.gm.energy = old
        for r in rs:
            if r["opening"] != oy or r["vs"] not in ("sub", "bl"):
                continue
            for lvl in ("low", "base"):
                out.append(dict(concept="Hyperloop elevated", case=f"spine+LAX all-stop M{o['Mach']:g}", vs=r["vs"],
                                capture=r["capture"], cost=lvl, opening=oy, bc_gross=r[f"bcg_{lvl}"],
                                bc_net=r[f"bcn_{lvl}"]))
    return out


def tube_eval(energy):
    return skyloop_rows(energy) + hyperloop_rows(energy)


def tubes(E, El):
    rows = []
    esaf_c = {sc: E[(E.scenario == sc)].set_index("case").esaf_usd_kg for sc in SCEN}
    firm = {sc: El[El.scenario == sc].set_index("case").firmed_24h_usd_mwh for sc in SCEN[1:]}
    # published reproduction
    base = pd.DataFrame(tube_eval("central"))
    sl = pd.read_csv(ROOT / "skyloop-review" / "system" / "lifecycle.csv")
    hlE = pd.read_csv(ROOT / "hyperloop" / "hl_econ.csv")
    for r in base.itertuples():
        if r.concept.startswith("Skyloop"):
            pub = sl[(sl.network == "spine+LAX all-stop") & (sl.M == 4) & (sl.vs == r.vs) & (sl.opening == 2050)
                     & (sl.capture == r.capture) & (sl.cost == r.cost)].bc_net.iloc[0]
        else:
            M = float(r.case.split("M")[-1])
            pub = hlE[(hlE.network == "spine+LAX all-stop") & (hlE.setting == "elevated") & (hlE.strategy == "full speed")
                      & (hlE.Mach == M) & (hlE.vs == r.vs) & (hlE.opening == 2050) & (hlE.capture == r.capture)][f"bcn_{r.cost}"].iloc[0]
        assert abs(r.bc_net - pub) < 2e-3, (r, pub)
    for d in base.to_dict("records"):
        rows.append(dict(d, scenario=SCEN[0], fuel_convention="canonical (growth central path)", electricity="repo grid $70/MWh + CO2"))
    for sc in SCEN:
        for k in CASES:
            P = esaf_c[sc][k]
            for d in tube_eval(f"esaf:{P}"):
                rows.append(dict(d, scenario=sc, fuel_convention=f"e-SAF {k}", fuel_all_in_usd_kg=P,
                                 electricity="repo grid $70/MWh + CO2"))
        if sc == SCEN[0]:
            continue
        P = esaf_c[sc]["central"]
        for k in CASES:
            with power(firm[sc][k], clean=True):
                for d in tube_eval(f"esaf:{P}"):
                    rows.append(dict(d, scenario=sc, fuel_convention="e-SAF central", fuel_all_in_usd_kg=P,
                                     electricity=f"firmed clean solar {k} ${firm[sc][k]:.0f}/MWh"))
        with power(S[sc]["solar"], clean=True):
            for d in tube_eval(f"esaf:{P}"):
                rows.append(dict(d, scenario=sc, fuel_convention="e-SAF central", fuel_all_in_usd_kg=P,
                                 electricity=f"daytime solar ${S[sc]['solar']:.0f}/MWh (bound, not 24 h)"))
    T = pd.DataFrame(rows)
    # break-even: net = gross (fuel credit = electricity cost) and low-cost net vs subsonic = 1, repo electricity
    be = []
    a = pd.DataFrame(tube_eval("esaf:0.0"))
    b = pd.DataFrame(tube_eval("esaf:1.0"))
    for (i, r0), (_, r1) in zip(a.iterrows(), b.iterrows()):
        if r0.capture != 1.0:
            continue
        slope = r1.bc_net - r0.bc_net
        key = dict(concept=r0.concept, case=r0.case, vs=r0.vs, capture=r0.capture, scenario_costs=f"{r0.cost} cost",
                   basis="life-cycle 2050 opening")
        be.append(dict(key, target="net = gross (fuel credit = electricity cost)",
                       esaf_breakeven_usd_kg=(r0.bc_gross - r0.bc_net) / slope, net_at_free_fuel=r0.bc_net,
                       note="repo grid electricity"))
        if r0.vs == "sub":
            ok = r0.bc_net < 1
            be.append(dict(key, target="net B/C = 1", esaf_breakeven_usd_kg=(1 - r0.bc_net) / slope if ok else np.nan,
                           net_at_free_fuel=r0.bc_net, note="repo grid electricity; above this price net >= 1" if ok
                           else "net >= 1 at any fuel price"))
    return T, be


# ------------------------------------------------------------------------------------------- main
def main():
    pd.set_option("display.width", 250)
    pd.set_option("display.max_columns", 30)
    E = esaf_table()
    El = electricity_table()
    D, costs = direct_table()
    B, be_bs = boomshield(costs, E)
    T, be_t = tubes(E, El)
    BE = pd.DataFrame(be_bs + be_t)
    BE["h2_needed_central_dac100"] = [esaf_h2_for(p, 100.0, 1) if p == p else np.nan for p in BE.esaf_breakeven_usd_kg]
    BE["h2_needed_lowcost_dac50"] = [esaf_h2_for(p, 50.0, 0) if p == p else np.nan for p in BE.esaf_breakeven_usd_kg]

    A = pd.DataFrame([
        ("scenario: Baseline", "repo as is: integration shield H2 $1.6/kg grey SMR + $0.5 delivery; cost_model base $4 + $1.5; "
         "growth T6 electrolytic $5-7/kg (G21); canonical fuel $0.80/kg + $190/t; ground PV $28/MWh (cost-down aggressive)"),
        ("scenario: H2 $2/kg", "Terraform Industries demo (company claim, unverified): sunlight->H2 < $2/kg now; Mar 2024 demo "
         "< $2.50/kg H2 and < $250/t DAC CO2; Rosamond PV ~$17/MWh (Jul 2026 column)"),
        ("scenario: H2 $1/kg", "Terraform plan '< $1/kg in the coming years' and DOE Hydrogen Shot $1/kg by 2031; DAC $100/t "
         "(DOE Carbon Negative Shot); PV $10/MWh (Terraform whitepaper assumption)"),
        ("source", "humanprogress.org/solar-powered-electrolyzer-slashes-hydrogen-production-costs (search summary, undated); "
         "critique (Liebreich-style): e-methane stays well above Henry Hub"),
        ("shield H2 delivery", "$0.5/kg central (integration), $0.3-1.5 range: compression, pipeline to bases, lift to 20 km"),
        ("e-SAF route", "RWGS + Fischer-Tropsch or methanol-to-jet; floor 0.432 kg H2 + 3.14 kg CO2 per kg; "
         "H2 eff 0.85/0.75/0.65, C eff 0.92/0.82/0.72; conversion $800/1,500/2,800 per kW fuel at CF 0.90/0.85/0.70, "
         "7 %/25 yr, O&M 3-5 %; H2 buffer $0.15/0.40/0.80 per kg; utilities $0.03-0.15/kg; distribution $0.03-0.08/kg [U]"),
        ("carbon-neutral convention", "extra (or avoided) fuel priced at the e-SAF price, flat in real terms; no CO2 charge on it; "
         "non-CO2 effects (contrails, NOx, stratospheric H2O) and H2 leakage still uncounted; grid electricity keeps its CO2 "
         "charge unless the row says 'clean'"),
        ("firmed solar", "daytime PV x overbuild 1.2/1.3/1.6 + share shifted 0.5/0.6/0.7 x battery ($120/180/250 per kWh, "
         "15-20 yr, 1 cycle/day) + round-trip loss [U]"),
        ("PV dual-use", "cost-down aggressive: 15 % of the lid carries PV, credit = (min(market, ground PV LCOE) - shield PV LCOE) "
         "x energy; 'capex follows' scales shield module + BOS by the ground LCOE ratio"),
        ("life-cycle", "growth_model.evaluate / metrics, base world, integrated cost, central efficiency path, 7 %, 25 yr"),
        ("tubes", "sky_model.lifecycle (Skyloop spine+LAX all-stop M4) and hyperloop.bc_rows (elevated full speed M3/M4), 2050 opening"),
    ], columns=["item", "value"])

    for name, df in (("esaf_cost", E), ("electricity", El), ("direct_h2", D), ("boomshield_net", B),
                     ("tubes_net", T), ("breakeven", BE), ("assumptions", A)):
        df.to_csv(HERE / f"{name}.csv", index=False, float_format="%.4g")

    # ---- console summary
    print(f"canonical all-in fossil fuel: ${ALLIN0:.3f}/kg")
    print(E[["scenario", "case", "h2_usd_kg", "dac_usd_t", "h2_usd", "co2_usd", "conversion_usd", "esaf_usd_kg",
             "abatement_usd_t"]].round(3).to_string())
    print("H2 at plant for e-SAF = $1.40/kg (central / low-cost conversion):",
          {d: (round(esaf_h2_for(ALLIN0, d), 2), round(esaf_h2_for(ALLIN0, d, 0), 2)) for d in (250.0, 100.0, 50.0)})
    print(El.round(1).to_string())
    print(D[D.variant.isin(["H2 + PV (scenario)", "H2 only", "PV only, shield PV capex follows",
                            "PV only, shield PV capex fixed", "share check", "gas fill"])
            | D.variant.str.startswith("H2 only, del")].round(4).to_string())
    k = B[(B.capture == 1.0) & B.design.isin(DES[1:])]
    print(k.pivot_table(index=["design", "vs", "basis"], columns=["scenario", "cost_variant", "fuel_convention"],
                        values="bc_net").round(3).T.to_string())
    print(BE.round(3).to_string())
    t = T[(T.capture == 1.0) & (T.vs == "sub")]
    print(t.pivot_table(index=["concept", "case", "cost"], columns=["scenario", "fuel_convention", "electricity"],
                        values="bc_net", aggfunc="first").round(3).T.to_string())


if __name__ == "__main__":
    main()
