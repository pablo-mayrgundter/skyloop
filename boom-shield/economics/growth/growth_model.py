#!/usr/bin/env python3
"""BoomShield traffic growth, value-of-time growth and timing (screening level).

Run:  python3 growth_model.py      (writes CSVs here; make_figure.py draws growth.png)

Reads (read-only):
  economics/cost_model.py                    discount / asset life defaults, VTTS
  economics/integration/results.csv          capex, opex, dual-use per km, per cost case (two tether rows)
  economics/integration/bcr.csv              integrated B/C for the 5,421 / 6,989 km networks
  network-maps/highmach/bcr_by_mach.csv      2025 benefit, extra fuel, extra CO2 (spine; M2.0 and M2.2)
  network-maps/demand/flights_vs_caps.csv    peak flights/day at 2025 traffic

Everything is real 2025 USD.  Benefits are the sibling studies' annual benefits at 2025 traffic and
$86/h, scaled by a traffic index, a VOT index and a ramp-up; costs are the integration study's capex
(spread over the build) + opex - dual-use over a 25-yr operating life.  [U] marks assumptions.
"""
from __future__ import annotations

import csv
import importlib.util
import math
from pathlib import Path

HERE = Path(__file__).resolve().parent
ECON = HERE.parent
BS = ECON.parent

# -------- read-only import of cost_model (discount rate, asset life) --------
_spec = importlib.util.spec_from_file_location("cost_model", ECON / "cost_model.py")
cm = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(cm)
R_BASE = cm.v("discount", "base")          # 0.07
LIFE = int(cm.v("asset_life_yr", "base"))  # 25

M20 = "M2.0 (Concorde-class, psf baseline)"
HM20 = "repo M2.0 (55 m, 1.75 km)"
HM22 = "A: M2.2 Al-limit"

# ---------------------------------------------------------------------------------------------
# 1. Traffic growth: US domestic long-haul O&D pax (> 1,500 km).  Rate g0 to 2045, then a linear
#    taper to g_inf by 2075 (maturity, slower population growth; CBO).  Forecasts exclude shocks;
#    the rates below already carry a ~0.3-0.4 pp shock/over-forecast haircut on FAA RPM growth.
# ---------------------------------------------------------------------------------------------
TRAFFIC = {   # g0, g_inf  (> 1,500 km pax);  > 2,500 km adds DL_2500
    "low":  (0.013, 0.005),   # FAA pessimistic domestic RPM 1.9 % less haircut; Japan/EU-style saturation
    "base": (0.023, 0.012),   # FAA FY2026-46 domestic RPM 2.7 % less haircut; = 2000-19 RPM history (~2.3 %)
    "high": (0.030, 0.020),   # FAA optimistic ~3.0 %; Boeing NA 2.8 %; premium-led
}
DL_2500 = 0.003               # trip length creep: > 2,500 km grows ~0.3 pp/yr faster [U] (RPM - enplanement gap)
PREMIUM_DL = 0.015            # premium long-haul grows ~1.5 pp/yr faster than all long-haul [U]
VOT = {"low": 0.005, "base": 0.010, "high": 0.015}   # real VTTS growth/yr: elasticity 0.5-1.0 x income 1.0-1.5 %
WORLDS = {"low": ("low", "low"), "base": ("base", "base"), "high": ("high", "high")}  # (traffic, VOT)

# Fuel / CO2 / SAF trajectories (real 2025 $).  Extra fuel is costed at the blended fuel price and
# its CO2 at the SC-CO2 net of SAF life-cycle savings.  Fleet efficiency shrinks the extra fuel.
ENERGY = {
    #            fossil $/kg 2025->2050, SCC $/t 2025->2050, SAF share 2030/2040/2050/2060, SAF $/kg 2025->2050, eff/yr
    "benign":    dict(pf=(0.80, 0.60), scc=(51, 75), saf=(0.00, 0.03, 0.10, 0.15), psaf=(2.0, 1.4), eff=0.010),
    "central":   dict(pf=(0.80, 0.95), scc=(200, 308), saf=(0.02, 0.20, 0.50, 0.70), psaf=(2.0, 1.6), eff=0.008),
    "stringent": dict(pf=(0.80, 1.40), scc=(340, 480), saf=(0.06, 0.34, 0.70, 1.00), psaf=(2.5, 2.0), eff=0.005),
}
SAF_CUT = 0.70                # life-cycle CO2 cut of SAF [U]
CO2_PER_KG = 3.16
FUEL_2025 = 0.80              # $/kg in highmach's extra_fuel_busd

# Cost cases (M2.0, two tether rows) from integration/results.csv
COST_KEYS = {
    "base": "base + 2nd tether row",
    "aggressive": "cost-down aggressive + 2nd tether row",
    "integrated": "integrated, aggressive levers, central",
    "integrated-plausible": "integrated, plausible levers, central",
}
COST_LEARN = {"low": 0.0, "base": 0.01, "high": 0.02}  # exogenous real decline/yr after 2030 [U]
COST_FLOOR_NET = 1.48         # $M/km/yr: integrated low-cost end; learning never goes below it

RAMP_YEARS, RAMP_START = 5, 0.5   # demand/capture ramp 50 % -> 100 % over 5 yr (demand.md caveat) [U]
CAP_BASE, CAP_LOW = 10_000, 1_000 # lid fatigue caps, flights/day on the busiest km


def lerp_year(y, y0, v0, y1, v1):
    if y <= y0:
        return v0
    if y >= y1:
        return v1
    return v0 + (v1 - v0) * (y - y0) / (y1 - y0)


def growth_rate(case, y, extra=0.0):
    g0, gi = TRAFFIC[case]
    return lerp_year(y, 2045, g0, 2075, gi) + extra


def traffic_index(case, year, extra=0.0):
    """Index (2025 = 1) of long-haul O&D pax."""
    x = 1.0
    for y in range(2026, int(year) + 1):
        x *= 1.0 + growth_rate(case, y, extra)
    return x


def vot_index(case, year):
    return (1.0 + VOT[case]) ** (year - 2025)


def energy(case, y):
    e = ENERGY[case]
    pf = lerp_year(y, 2025, e["pf"][0], 2050, e["pf"][1])
    scc = lerp_year(y, 2025, e["scc"][0], 2050, e["scc"][1]) if y <= 2050 else e["scc"][1] * 1.015 ** (y - 2050)
    s30, s40, s50, s60 = e["saf"]
    s = (lerp_year(y, 2025, 0.0, 2030, s30) if y < 2030 else lerp_year(y, 2030, s30, 2040, s40) if y < 2040
         else lerp_year(y, 2040, s40, 2050, s50) if y < 2050 else lerp_year(y, 2050, s50, 2060, s60))
    psaf = lerp_year(y, 2025, e["psaf"][0], 2050, e["psaf"][1])
    fuel_price = (1 - s) * pf + s * psaf
    co2_factor = (1 - s) + s * (1 - SAF_CUT)
    return fuel_price, scc, co2_factor, (1 - e["eff"]) ** (y - 2025)


# ---------------------------------------------------------------------------------------------
# Inputs from sibling studies
# ---------------------------------------------------------------------------------------------
def read_csv(p):
    with open(p, newline="") as f:
        return list(csv.DictReader(f))


def load_inputs():
    res = {(r["design"], r["scenario"]): r for r in read_csv(ECON / "integration" / "results.csv")}
    cost = {}
    for k, key in COST_KEYS.items():
        r = res[(M20, key)]
        cost[k] = dict(capex=float(r["capex_musd_km"]), opex=float(r["opex_musd_km_yr"]),
                       dual=float(r["dual_use_musd_km_yr"]), net=float(r["net_musd_km_yr"]))
    hm = read_csv(BS / "network-maps" / "highmach" / "bcr_by_mach.csv")
    sm = {r["design"]: r for r in read_csv(BS / "network-maps" / "highmach" / "summary_by_mach.csv")}

    def hrow(des, net, dem, vs):
        for r in hm:
            if (r["design"] == des and r["network"] == net and r["demand"] == dem and r["vs"] == vs
                    and float(r["capture"]) == 1.0 and r["cost"] == "aggressive"):
                return {k: float(r[k]) for k in ("benefit_busd", "extra_fuel_busd", "extra_co2_Mt", "km")}
        raise KeyError((des, net, dem, vs))

    bcr_int = read_csv(ECON / "integration" / "bcr.csv")

    def int_row(net, dem="base"):
        for r in bcr_int:
            if (r["design"] == M20 and r["network"] == net and r["demand"] == dem
                    and r["cost_case"] == "integrated, aggressive levers (central)"):
                return r
        raise KeyError(net)

    fl = {(r["network"], r["variant"]): float(r["peak_flights_day"])
          for r in read_csv(BS / "network-maps" / "demand" / "flights_vs_caps.csv")}

    nets = {}
    for dem in ("base", "high"):
        for vs in ("bl", "sub"):
            # M2.0 spine straight from highmach (same numbers as integration/demand)
            nets[("M2.0", "spine", dem, vs)] = dict(hrow(HM20, "spine", dem, vs), cost_scale=1.0,
                                                    flights=fl[("spine", dem)], build=10)
            # M2.2 sweet spot: cost scaled by highmach's aggressive cost ratio
            r22 = hrow(HM22, "spine", dem, vs)
            scale22 = float(sm[HM22]["cost_aggr"]) / float(sm[HM20]["cost_aggr"])
            fl22 = fl[("spine", dem)] * 4961 / 4500   # highmach network_by_mach peak flights ratio
            nets[("M2.2", "spine", dem, vs)] = dict(r22, cost_scale=scale22, flights=fl22, build=10)
            # M2.0 knee / shoulder: benefit from integration bcr (B/C x cost), fuel/CO2 scaled by M2.2 ratios
            sp22 = hrow(HM22, "spine", dem, vs)
            for net, km, flk in (("knee 5,421 km", 5421.4, None), ("shoulder 6,989 km", 6988.8, "greedy step 6")):
                ir = int_row(net, dem if dem == "base" else "base")
                ben_base = float(ir["bcr_bl_100" if vs == "bl" else "bcr_sub_100"]) * float(ir["net_musd_km_yr"]) * km / 1e3
                if dem == "high":  # integration has base demand only; apply highmach's M2.2 high/base uplift
                    n22b, n22h = hrow(HM22, net, "base", vs), hrow(HM22, net, "high", vs)
                    ben_base *= n22h["benefit_busd"] / n22b["benefit_busd"]
                n22 = hrow(HM22, net, dem, vs)
                sp20 = hrow(HM20, "spine", dem, vs)
                f = n22["extra_fuel_busd"] / sp22["extra_fuel_busd"]
                c = n22["extra_co2_Mt"] / sp22["extra_co2_Mt"]
                flights = fl[(flk, dem)] if flk else fl[("spine", dem)] * 1.2  # knee: between spine and shoulder [U]
                nets[("M2.0", net.split()[0], dem, vs)] = dict(
                    benefit_busd=ben_base, extra_fuel_busd=sp20["extra_fuel_busd"] * f,
                    extra_co2_Mt=sp20["extra_co2_Mt"] * c, km=km,
                    cost_scale=float(ir["net_musd_km_yr"]) / cost["integrated"]["net"], int_only_scale=True,
                    flights=flights, build=15)
    return cost, nets



# ---------------------------------------------------------------------------------------------
# 2. Life-cycle evaluation for one opening year
# ---------------------------------------------------------------------------------------------
def evaluate(cost, net, open_year, world="base", energy_case="central", cost_case="integrated",
             capture=1.0, learn="base", r=R_BASE, cap=CAP_BASE, premium=False):
    tcase, vcase = WORLDS[world]
    c = cost[cost_case]
    scale = net["cost_scale"]
    if net.get("int_only_scale") and not cost_case.startswith("integrated"):
        scale = 1.0           # knee/shoulder plant-sharing discount applies only to integrated cases
    km = net["km"]
    build = net["build"]
    mid_build = open_year - build / 2
    lf = (1 - COST_LEARN[learn]) ** max(0.0, mid_build - 2030)
    lf = max(lf, min(1.0, COST_FLOOR_NET / (c["net"] * scale)))
    capex_total = c["capex"] * scale * lf * km / 1e3          # $B
    opex_net = (c["opex"] - c["dual"]) * scale * lf * km / 1e3  # $B/yr
    flows = []   # (year, cost, gross_benefit, fuel, co2)
    for y in range(open_year - build, open_year):
        flows.append((y, capex_total / build, 0.0, 0.0, 0.0, 0.0))
    for i, y in enumerate(range(open_year, open_year + LIFE)):
        ramp = min(1.0, RAMP_START + (1 - RAMP_START) * i / RAMP_YEARS)
        T = traffic_index(tcase, y)
        k = capture
        if premium:   # capture share drifts up with the premium segment's faster growth
            k = min(1.0, capture * (1 + PREMIUM_DL) ** (y - 2025))
        flights = net["flights"] * k * T * ramp
        capf = min(1.0, cap / flights) if flights > 0 else 1.0
        vol = k * ramp * T * capf
        ben = net["benefit_busd"] * vol * vot_index(vcase, y)
        pfuel, scc, co2f, eff = energy(energy_case, y)
        fuel_kg_b = net["extra_fuel_busd"] / FUEL_2025 * vol * eff        # bn kg
        fuel = fuel_kg_b * pfuel
        co2 = net["extra_co2_Mt"] * vol * eff * co2f * scc / 1e3          # $B
        flows.append((y, opex_net, ben, fuel, co2, flights))
    return flows


def pv(flows, idx, r, base_year):
    return sum(f[idx] / (1 + r) ** (f[0] - base_year) for f in flows)


def irr(cash):
    """cash: list of (year, value).  Bisection on [-0.5, 1.0]; None if no sign change."""
    y0 = cash[0][0]

    def npv(rate):
        return sum(v / (1 + rate) ** (y - y0) for y, v in cash)
    lo, hi = -0.5, 1.0
    flo, fhi = npv(lo), npv(hi)
    if flo * fhi > 0:
        return None
    for _ in range(200):
        mid = (lo + hi) / 2
        fm = npv(mid)
        if flo * fm <= 0:
            hi, fhi = mid, fm
        else:
            lo, flo = mid, fm
    return (lo + hi) / 2


def metrics(flows, open_year, r=R_BASE):
    C = pv(flows, 1, r, open_year)
    B = pv(flows, 2, r, open_year)
    F = pv(flows, 3, r, open_year)
    E = pv(flows, 4, r, open_year)
    disc26 = (1 + r) ** (open_year - 2026)
    gross_cash = [(f[0], f[2] - f[1]) for f in flows]
    net_cash = [(f[0], f[2] - f[3] - f[4] - f[1]) for f in flows]
    peak = max(f[5] for f in flows)
    return dict(bc_gross=B / C, bc_net=(B - F - E) / C,
                npv26_gross=(B - C) / disc26, npv26_net=(B - F - E - C) / disc26,
                irr_gross=irr(gross_cash), irr_net=irr(net_cash),
                pv_cost=C / disc26, peak_flights_last=flows[-1][5], peak_flights_max=peak)


# ---------------------------------------------------------------------------------------------
# 3. Outputs
# ---------------------------------------------------------------------------------------------
def write(name, rows):
    if not rows:
        return
    with open(HERE / name, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        w.writeheader()
        w.writerows(rows)


def years_to(mult, case, extra=0.0, horizon=2400):
    x, y = 1.0, 2025
    while x < mult and y < horizon:
        y += 1
        x *= 1 + growth_rate(case, y, extra)
    return y if x >= mult else None


def fmt(v, d=3):
    return None if v is None else round(v, d)


def main():
    cost, nets = load_inputs()

    # --- traffic multiples ---
    rows = []
    for seg, extra in ((">1,500 km", 0.0), (">2,500 km", DL_2500), ("premium >1,500 km", PREMIUM_DL)):
        for case in ("low", "base", "high"):
            r = {"segment": seg, "case": case, "g_2026_45_pct": round(100 * (TRAFFIC[case][0] + extra), 2),
                 "g_after_2075_pct": round(100 * (TRAFFIC[case][1] + extra), 2)}
            for m in (1.5, 2, 3, 5, 10):
                r[f"year_{m}x"] = years_to(m, case, extra)
            for y in (2035, 2045, 2050, 2060, 2075):
                r[f"index_{y}"] = round(traffic_index(case, y, extra), 2)
            rows.append(r)
    write("traffic_multiples.csv", rows)
    const = []
    for g in (0.01, 0.015, 0.02, 0.025, 0.03, 0.04):
        const.append({"constant_growth_pct": 100 * g, **{f"years_to_{m}x": round(math.log(m) / math.log(1 + g), 1)
                                                         for m in (1.5, 2, 3, 5, 10)}})
    write("constant_rate_multiples.csv", const)

    # --- trajectories (for the figure) ---
    traj = []
    for y in range(2025, 2101):
        r = {"year": y}
        for case in ("low", "base", "high"):
            r[f"traffic_{case}"] = round(traffic_index(case, y), 4)
            r[f"vot_{case}"] = round(vot_index(case, y), 4)
        for ec in ENERGY:
            pfuel, scc, co2f, eff = energy(ec, y)
            r[f"fuel_usd_kg_{ec}"] = round(pfuel, 3)
            r[f"scc_{ec}"] = round(scc, 1)
            r[f"allin_usd_per_kg_extra_fuel_{ec}"] = round(pfuel + CO2_PER_KG * co2f * scc / 1e3, 3)
        traj.append(r)
    write("trajectories.csv", traj)

    # --- full grid ---
    grid = []
    opening = list(range(2030, 2101, 5))
    for (des, net, dem, vs), nd in nets.items():
        for world in WORLDS:
            for ec in ENERGY:
                for cc in COST_KEYS:
                    for capt in (0.3, 1.0):
                        for oy in opening:
                            m = metrics(evaluate(cost, nd, oy, world, ec, cc, capt, learn=world), oy)
                            grid.append({"design": des, "network": net, "demand": dem, "vs": vs, "world": world,
                                         "energy": ec, "cost": cc, "capture": capt, "open_year": oy,
                                         **{k: fmt(v, 4) for k, v in m.items()}})
    write("roi_grid.csv", grid)

    # --- headline table: opening years x cost, spine M2.0, base world, central energy ---
    head = []
    for des in ("M2.0", "M2.2"):
        for vs in ("bl", "sub"):
            for cc in ("base", "aggressive", "integrated"):
                for capt in (0.3, 1.0):
                    for oy in (2035, 2040, 2045, 2050):
                        m = metrics(evaluate(cost, nets[(des, "spine", "base", vs)], oy, "base", "central", cc, capt), oy)
                        head.append({"design": des, "vs": vs, "cost": cc, "capture": capt, "open_year": oy,
                                     "bc_gross": fmt(m["bc_gross"]), "bc_net": fmt(m["bc_net"]),
                                     "npv26_gross_busd": fmt(m["npv26_gross"], 1), "npv26_net_busd": fmt(m["npv26_net"], 1),
                                     "irr_gross": fmt(m["irr_gross"]), "irr_net": fmt(m["irr_net"])})
    write("roi_headline.csv", head)

    # --- static check: no growth, no ramp, 2025 prices (compare with sibling annualised B/C) ---
    static = []
    for cc in ("base", "aggressive", "integrated"):
        nd = nets[("M2.0", "spine", "base", "bl")]
        c = cost[cc]
        ann = (c["capex"] * cm.crf(R_BASE, LIFE) + c["opex"] - c["dual"]) * nd["km"] / 1e3
        static.append({"cost": cc, "sibling_bc_bl_100": round(nd["benefit_busd"] / ann, 3)})
    write("static_check.csv", static)

    # --- crossing years ---
    cross = []
    keyset = sorted({(g["design"], g["network"], g["demand"], g["vs"], g["world"], g["energy"], g["cost"], g["capture"])
                     for g in grid})
    idx = {}
    for g in grid:
        idx.setdefault((g["design"], g["network"], g["demand"], g["vs"], g["world"], g["energy"], g["cost"], g["capture"]),
                       []).append(g)
    for k in keyset:
        rows_k = sorted(idx[k], key=lambda g: g["open_year"])
        def first(col):
            for g in rows_k:
                if g[col] is not None and g[col] >= 1.0:
                    return g["open_year"]
            return None
        cross.append(dict(zip(("design", "network", "demand", "vs", "world", "energy", "cost", "capture"), k),
                          cross_gross=first("bc_gross"), cross_net=first("bc_net"),
                          bc_gross_2045=next(g["bc_gross"] for g in rows_k if g["open_year"] == 2045),
                          bc_net_2045=next(g["bc_net"] for g in rows_k if g["open_year"] == 2045),
                          bc_gross_2070=next(g["bc_gross"] for g in rows_k if g["open_year"] == 2070),
                          bc_net_2070=next(g["bc_net"] for g in rows_k if g["open_year"] == 2070)))
    write("crossing_years.csv", cross)

    # --- lid cap and airports vs time ---
    caps = []
    for net, base_f in (("spine", nets[("M2.0", "spine", "base", "bl")]["flights"]),
                        ("shoulder 6,989 km", nets[("M2.0", "shoulder", "base", "bl")]["flights"])):
        for case in ("low", "base", "high"):
            for capt in (0.3, 1.0):
                def yr(level):
                    for y in range(2025, 2201):
                        if base_f * capt * traffic_index(case, y) >= level:
                            return y
                    return None
                caps.append({"network": net, "traffic": case, "capture": capt,
                             "flights_2025": round(base_f * capt), "flights_2050": round(base_f * capt * traffic_index(case, 2050)),
                             "flights_2075": round(base_f * capt * traffic_index(case, 2075)),
                             "year_reach_1000": yr(CAP_LOW), "year_reach_10000": yr(CAP_BASE),
                             "LAX_supersonic_movements_2050": round(1160 * capt * traffic_index(case, 2050) * (1 if net == "spine" else 1.0))})
    write("caps_vs_time.csv", caps)

    # --- sensitivities at a 2045 opening, spine M2.0 integrated, base world, central energy ---
    sens = []
    nd_bl = nets[("M2.0", "spine", "base", "bl")]
    nd_sub = nets[("M2.0", "spine", "base", "sub")]
    def s(label, **kw):
        a = dict(open_year=2045, world="base", energy_case="central", cost_case="integrated", capture=1.0)
        a.update(kw)
        oy = a.pop("open_year")
        r = a.pop("r", R_BASE)
        mb = metrics(evaluate(cost, nd_bl, oy, r=r, **a), oy, r)
        ms = metrics(evaluate(cost, nd_sub, oy, r=r, **a), oy, r)
        sens.append({"case": label, "bl_gross": fmt(mb["bc_gross"]), "bl_net": fmt(mb["bc_net"]),
                     "sub_gross": fmt(ms["bc_gross"]), "sub_net": fmt(ms["bc_net"])})
    s("central (2045, integrated, 100 %, 7 %)")
    s("discount 3 %", r=0.03)
    s("low world (traffic + VOT)", world="low")
    s("high world", world="high")
    s("benign energy", energy_case="benign")
    s("stringent energy", energy_case="stringent")
    s("no exogenous cost learning", learn="low")
    s("fatigue cap 1,000/day", cap=CAP_LOW)
    s("30 % capture", capture=0.3)
    s("30 % capture, premium-led drift", capture=0.3, premium=True)
    s("integrated-plausible cost ($17.6/m2)", cost_case="integrated-plausible")
    nd_hi = nets[("M2.0", "spine", "high", "bl")]
    m = metrics(evaluate(cost, nd_hi, 2045), 2045)
    m2 = metrics(evaluate(cost, nets[("M2.0", "spine", "high", "sub")], 2045), 2045)
    sens.append({"case": "high induced demand", "bl_gross": fmt(m["bc_gross"]), "bl_net": fmt(m["bc_net"]),
                 "sub_gross": fmt(m2["bc_gross"]), "sub_net": fmt(m2["bc_net"])})
    write("sensitivities_2045.csv", sens)
    print("wrote CSVs to", HERE)


if __name__ == "__main__":
    main()
