#!/usr/bin/env python3
"""Mach summary table for highmach.md: spine, 100 % / 10 % capture, base and high demand, three cost levels,
plus the vertical-integration sibling (economics/integration/results.csv) as a cost ratio.
Writes summary_by_mach.csv, emissions_by_mach.csv."""
from __future__ import annotations

from pathlib import Path

import pandas as pd

HERE = Path(__file__).resolve().parent
INT = HERE.parents[1] / "economics" / "integration" / "results.csv"

R = pd.read_csv(HERE / "bcr_by_mach.csv")
D = pd.read_csv(HERE / "mach_designs.csv").set_index("design")
T = pd.read_csv(HERE / "door_to_door.csv")
t = T[(T.pair == "NYC-LAX") & (T.overhead == "base (130 min)")].set_index("mode").door_to_door_h

# integration: integrated aggressive central vs cost-down aggressive + 2nd row, M2.0 (applied as a ratio) [scaling]
int_ratio, int_note = None, "economics/integration/ absent"
if INT.exists():
    I = pd.read_csv(INT)
    m = I[I.design.str.startswith("M2.0")].set_index("scenario").net_musd_km_yr
    int_ratio = m["integrated, aggressive levers, central"] / m["cost-down aggressive + 2nd tether row"]
    int_note = f"integrated/aggressive cost ratio {int_ratio:.3f} (M2.0), applied to every Mach"

MODE = {"M1.7 Overture": "shield M1.7", "M2.0 adaptive": "shield M2.0", "A: M2.2 Al-limit": "shield M2.2",
        "M2.5 Ti adaptive": "shield M2.5", "B: M3.0 Ti adaptive": "shield M3.0", "C: M3.5 TBCC": "shield M3.5",
        "D: M4.0 TBCC (+40 K gas)": "shield M4.0", "repo M2.0 (55 m, 1.75 km)": "shield M2.0"}
rows = []
for des in D.index:
    s = R[(R.design == des) & (R.network == "spine")]
    g = lambda dem, vs, cap, cost, col="bcr": float(s[(s.demand == dem) & (s.vs == vs) & (s.capture == cap)
                                                      & (s.cost == cost)][col].iloc[0])
    row = dict(design=des, M=D.loc[des, "M"], nyc_lax_d2d_h=t[MODE[des]], W_mid_t=D.loc[des, "W_mid_t"],
               H_s_km=D.loc[des, "H_s_km"], d_m=D.loc[des, "d_m"], width_km=D.loc[des, "width_km"],
               cost_base=D.loc[des, "cost_base_musd_km_yr"], cost_aggr=D.loc[des, "cost_aggressive_musd_km_yr"],
               pax_h_M=g("base", "bl", 1.0, "base", "pax_h_M"), day_return_pax_M=g("base", "bl", 1.0, "base", "day_return_pax_M"),
               net_induced_M=g("base", "bl", 1.0, "base", "net_induced_M"))
    for vs in ("bl", "sub"):
        for lvl in ("base", "plausible", "aggressive"):
            row[f"{vs}_{lvl}_100"] = g("base", vs, 1.0, lvl)
        row[f"{vs}_aggr_10"] = g("base", vs, 0.1, "aggressive")
        row[f"{vs}_aggr_high_100"] = g("high", vs, 1.0, "aggressive")
        row[f"{vs}_aggr_100_onerow"] = g("base", vs, 1.0, "aggressive", "bcr_time_onerow") * \
            g("base", vs, 1.0, "aggressive") / g("base", vs, 1.0, "aggressive", "bcr_time")
        row[f"{vs}_aggr_100_net_fuel_co2"] = g("base", vs, 1.0, "aggressive", "bcr_net_fuel_carbon")
        if int_ratio:
            row[f"{vs}_integrated_100"] = row[f"{vs}_aggressive_100"] / int_ratio
            row[f"{vs}_integrated_high_100"] = row[f"{vs}_aggr_high_100"] / int_ratio
    # break-even value of time net of extra fuel (vs Boomless, base cost irrelevant): VOT* = 86 x fuel / time value
    q = s[(s.demand == "base") & (s.vs == "bl") & (s.capture == 1.0) & (s.cost == "base")].iloc[0]
    time_val = q.bcr_time * q.cost_busd
    row["vot_breakeven_fuel_bl"] = 86.0 * q.extra_fuel_busd / time_val if time_val > 0 else float("nan")
    row["vot_breakeven_fuel_co2_bl"] = 86.0 * (q.extra_fuel_busd + q.carbon_busd) / time_val
    row["extra_co2_vs_bl_Mt"] = q.extra_co2_Mt
    qs = s[(s.demand == "base") & (s.vs == "sub") & (s.capture == 1.0) & (s.cost == "base")].iloc[0]
    row["extra_co2_vs_sub_Mt"] = qs.extra_co2_Mt
    # per seat-km vs subsonic: domestic-sized Overture-class reference is 2.17x in this model (as-built Overture 3.0x,
    # ICCT 5x); 'icct' scales the model by 5/3
    ref_ratio = 2.165
    row["co2_per_seat_km_vs_subsonic_model"] = ref_ratio * D.loc[des, "fuel_ratio_vs_overture"]
    row["co2_per_seat_km_vs_subsonic_icct"] = ref_ratio * D.loc[des, "fuel_ratio_vs_overture"] * 5 / 3
    rows.append(row)
S = pd.DataFrame(rows)
S.to_csv(HERE / "summary_by_mach.csv", index=False, float_format="%.4g")
pd.set_option("display.width", 250)
print(int_note)
print(S.T.round(4).to_string())
