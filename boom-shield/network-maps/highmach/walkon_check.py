"""Spine B/C vs Boomless (base cost, two rows, 100 % capture) with demand's base scenario but 100 / 70 min
ground overhead (applies to both worlds; only threshold credits and induced response change)."""
import pickle, sys
import pandas as pd
sys.dont_write_bytecode = True
import hm_network as hn
inp = pickle.loads((hn.NM / "demand" / ".cache" / "inputs.pkl").read_bytes())
ev = pickle.loads((hn.CACHE / "net_evals.pkl").read_bytes())
designs = {d["design"]: d for d in hn.mach_designs()}
rows = []
for des in ("M1.7 Overture", "M2.0 adaptive", "A: M2.2 Al-limit", "B: M3.0 Ti adaptive", "C: M3.5 TBCC"):
    ds, ne = designs[des], dict(ev[f"{des}|spine"])
    for O in (130, 100, 70):
        sc = dict(hn.bd.SCEN["base"], O=O)
        hn.bd.SCEN["walk"] = sc
        r = [x for x in hn.bcr_rows(inp, ds, ne, "walk") if x["vs"] == "bl" and x["capture"] == 1.0]
        for x in r:
            rows.append(dict(design=des, overhead_min=O, cost=x["cost"], bcr=x["bcr"], day_return_pax_M=x["day_return_pax_M"],
                             bcr_net_fuel_carbon=x["bcr_net_fuel_carbon"]))
W = pd.DataFrame(rows)
W.to_csv(hn.HERE / "walkon_check.csv", index=False, float_format="%.4g")
print(W[W.cost.isin(["base", "aggressive"])].pivot_table(index=["design", "overhead_min"], columns="cost", values=["bcr", "day_return_pax_M"]).round(3))
