"""BoomShield lid fatigue, step 4: trunk B/C with revised lid caps.

Re-uses build_trunk.py's cached per-sub-segment loads (network-maps/trunk/.cache/
ev_greedy_M{M}_M{M}.pkl: lb = {segment: [pax/yr, pax-h/yr vs Boomless]}, ls = same vs
subsonic, same-aircraft base-demand runs) and applies the same capacity cap as
build_trunk.summarise():  captured pax on every km <= cap flights/day x 65 seats x LF x 365.
Costs: cost_per_km_by_mach.csv (base case, cost_model.py).  Value of time and availability
are imported from economics/cost_model.py.  No source file is modified.
Writes trunk_bcr_revised.csv.
"""
from __future__ import annotations

import os
import pickle
import sys

import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
BS = os.path.dirname(HERE)
TRUNK = os.path.join(BS, "network-maps", "trunk")
sys.path.insert(0, os.path.join(BS, "economics"))
import cost_model as cm  # noqa: E402

CAPS = {"none": None, "old 100/day": 100, "revised low 1,000/day": 1000,
        "revised base 10,000/day": 10000, "revised high (no fatigue cap)": None}
CAPTURES = [0.10, 1.00]
MACHS = [1.4, 1.7, 2.0, 2.5, 3.2]
SEATS, LF = cm.SEATS, cm.LOAD_FACTOR


def capped(loads, c, cap):
    ben = 0.0
    for p, ph in loads.values():
        f = 1.0 if cap is None else min(1.0, cap * SEATS * LF * 365 / max(p * c, 1))
        ben += ph * c * f
    return ben


def main():
    cost = pd.read_csv(os.path.join(TRUNK, "cost_per_km_by_mach.csv"))
    summ = pd.read_csv(os.path.join(TRUNK, "trunk_summary.csv"))
    vt = cm.VTTS["base"] * cm.v("availability", "base")
    rows = []
    for M in MACHS:
        d = pickle.loads(open(os.path.join(TRUNK, ".cache", f"ev_greedy_M{M}_M{M}.pkl"), "rb").read())
        G = d["G"]
        km = sum(dd["length"] for nb in G.values() for dd in nb.values()) / 2   # dict-of-dicts, each edge twice
        ann = float(cost[(cost.mach == M) & (cost.case == "base")].annual_musd_km_yr.iloc[0]) * 1e6
        ref = summ[(summ.tree == f"greedy_M{M}") & (summ.mach == M) & (summ.variant == "same_aircraft")
                   & (summ.demand == "base")].iloc[0]
        busiest = max(v[0] for v in d["lb"].values()) / (SEATS * LF * 365)
        for cname, cap in CAPS.items():
            for c in CAPTURES:
                bl = capped(d["lb"], c, cap) * vt / (ann * km)
                sub = capped(d["ls"], c, cap) * vt / (ann * km)
                rows.append(dict(mach=M, trunk_km=round(km), cap=cname, cap_flights_day=cap or float("nan"),
                                 capture=c, bcr_vs_boomless=bl, bcr_vs_subsonic=sub,
                                 busiest_flights_day_at_capture=busiest * c,
                                 check_build_trunk_uncapped=ref[f"bcr_vs_bl_cap{int(c * 100)}"],
                                 check_build_trunk_cap100=ref[f"bcr_vs_bl_cap{int(c * 100)}_lidcapped"]))
    df = pd.DataFrame(rows)
    df.to_csv(os.path.join(HERE, "trunk_bcr_revised.csv"), index=False, float_format="%.4g")
    with pd.option_context("display.width", 220, "display.max_rows", 100):
        print(df.to_string(float_format=lambda v: f"{v:.4g}"))


if __name__ == "__main__":
    main()
