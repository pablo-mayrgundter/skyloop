#!/usr/bin/env python3
"""Sensitivities for the headline networks (M4, 2050 opening, 100 % capture, base cost): ground overhead,
acceleration, headway rule (brick-wall -> more bores), energy path, grid, discount rate, life, airline-cost credit.
Run after sky_model.py.  Writes sensitivities.csv."""
from __future__ import annotations

import math

import numpy as np
import pandas as pd

import sky_model as sm

HEAD = {"spine+LAX all-stop": sm.network_defs()["spine+LAX all-stop"],
        "NYC-LAX direct": sm.network_defs()["NYC-LAX direct"],
        "NYC-MIA direct": sm.network_defs()["NYC-MIA direct"]}


def one(inp, name, lines, cmp, M=4.0, O_sky=None, a=None, rule="moving block 10 s", seats=None, lvl="base",
        oy=2050, cap=1.0, energy="central", grid="central", r=0.07, life=25, resource=False, world="base",
        align_gw_per_1000km=0.0):
    old_a, old_seats = sm.P["a_acc"], sm.P["seats"]
    if a:
        sm.P["a_acc"] = a
    if seats:
        sm.P["seats"] = seats
    try:
        net = sm.Net(name, lines, M, inp["df"])
        ev = sm.evaluate(inp, net, cmp, O_sky=O_sky)
        capyr = sm.cap_pair_yr(M, rule)[0]
        bores = int(min(12, max(1, math.ceil(ev["peak"] * 1.3 / capyr))))
        cx = sm.net_capex(net, lvl, bores=bores)
        pods = sm.pods_capex(ev, net, lvl)
        capex = cx["total_ex_pods"] + pods
        t = ev["tot"]
        ann = dict(benefit_busd=t["benefit_busd"], fuel_saved_Mt=t["fuel_saved_Mt"], kwh_TWh=t["kwh_TWh"],
                   pax_M=t["captured_pax_M"] + t["net_induced_M"], pax_km_B=t["pax_km_B"], peak=ev["peak"])
        # active alignment thrust (physics: 0.4-2 GW per 1,000 km, ideal) as extra opex at grid price
        extra = align_gw_per_1000km * net.km() / 1000 * 8.76 * sm.P["elec_usd_kwh"]          # $B/yr
        lc = sm.lifecycle(ann, capex, sm.OM_FRAC[sm.LVL[lvl]] + extra / capex, oy, cap, world=world, energy_case=energy,
                          grid_case=grid, r=r, life=life, build=10 if net.km() > 2000 else 7,
                          cap_pax=bores * capyr, resource=resource)
        return dict(bores=bores, capex_busd=capex, benefit_busd=t["benefit_busd"], bc_gross=lc["bc_gross"],
                    bc_net=lc["bc_net"])
    finally:
        sm.P["a_acc"], sm.P["seats"] = old_a, old_seats


if __name__ == "__main__":
    inp = sm.load_inputs()
    cases = [("central", {}),
             ("ground overhead 70 min (walk-on)", dict(O_sky=(30, 20, 20))),
             ("ground overhead 130 min (airport-like)", dict(O_sky=(40, 55, 35))),
             ("20 km station by lift (+49 min)", dict(O_sky=(35 + 24, 40, 25 + 24))),
             ("accel 0.15 g", dict(a=0.15 * sm.G)),
             ("brick-wall headway, 50-seat pods", dict(rule="brick-wall")),
             ("brick-wall headway, 150-seat pods", dict(rule="brick-wall", seats=150)),
             ("30 % capture", dict(cap=0.3)),
             ("low cost", dict(lvl="low")),
             ("high cost", dict(lvl="high")),
             ("opening 2045", dict(oy=2045)),
             ("opening 2060", dict(oy=2060)),
             ("benign energy (cheap fuel, low SCC)", dict(energy="benign")),
             ("stringent energy", dict(energy="stringent")),
             ("grid stalls at ~0.3 kg/kWh", dict(grid="stalled")),
             ("clean grid 0.02 kg/kWh", dict(grid="clean")),
             ("discount 3 %", dict(r=0.03)),
             ("life 50 yr", dict(life=50)),
             ("+ avoided airline non-fuel cost - Skyloop ops", dict(resource=True)),
             ("high world", dict(world="high")),
             ("active alignment thrust 1 GW per 1,000 km (physics)", dict(align_gw_per_1000km=1.0))]
    rows = []
    for name, lines in HEAD.items():
        for cmp in ("sub", "bl", "m22"):
            for lab, kw in cases:
                res = one(inp, name, lines, cmp, **kw)
                rows.append(dict(network=name, vs=cmp, case=lab, **res))
                print(name, cmp, lab, {k: round(v, 3) for k, v in res.items()}, flush=True)
            pd.DataFrame(rows).to_csv(sm.HERE / "sensitivities.csv", index=False, float_format="%.4g")
    for M in (2.0, 3.0, 6.0):
        for name, lines in HEAD.items():
            for cmp in ("sub", "bl"):
                for lab, kw in (("brick-wall headway, 50-seat pods", dict(rule="brick-wall")),
                                ("central", {})):
                    res = one(inp, name, lines, cmp, M=M, **kw)
                    rows.append(dict(network=name, vs=cmp, case=f"M{M:g}: {lab}", **res))
    pd.DataFrame(rows).to_csv(sm.HERE / "sensitivities.csv", index=False, float_format="%.4g")
