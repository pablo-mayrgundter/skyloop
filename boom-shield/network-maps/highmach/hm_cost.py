#!/usr/bin/env python3
"""Shield cost per km for an arbitrary layer (thickness, width, altitude) via cost-down's evaluator (read-only).

cost_model sizes the lid and H2 at 20 km (lid = 0.0827 T - hold-down - floor; H2 = rho_g(20 km) T). A layer
of thickness d at altitude H is passed as the 20 km-equivalent thickness T_eq = d * lift(H) / lift(20 km):
lid mass, H2 mass and the tether side load (~ rho * T) all scale with density, so they come out right;
partition film (~T) is slightly understated above 20 km.
"""
from __future__ import annotations

import copy
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
CD = HERE.parents[1] / "economics" / "cost-down"
sys.dont_write_bytecode = True
sys.path.insert(0, str(CD))
import cost_down as cdn  # noqa: E402

cm = cdn.cm
REF20 = cdn.M20


def annual(T_eq, width_km, level="base", two_rows=True, key=None):
    """Net $M/km/yr (and $/m2 capex) for the cost-down scenario `level` (base / plausible / aggressive)."""
    key = key or f"HM T{T_eq:.1f} W{width_km:.3f}"
    cm.DESIGNS[key] = {"T": float(T_eq), "W": (width_km, width_km, width_km)}
    S = copy.deepcopy(cdn.scenario(level)[-1][1])
    S["T"] = dict(S["T"])
    S["T"][key] = float(T_eq)
    S["two_rows"] = two_rows
    r = cdn.evaluate(key, S)
    return dict(net_musd_km_yr=r["net_musd_km_yr"], capex_usd_m2=r["usd_m2"], capex_musd_km=r["capex_musd_km"],
                h2_t_km=r["h2_t_km"], lid_kg_m2=r["lid_kg_m2"])


if __name__ == "__main__":
    # reproduce cost-down: M2.0 55 m x 1.75 km
    for lvl in ("base", "plausible", "aggressive"):
        print(lvl, "one row", round(annual(55.0, 1.75, lvl, False)["net_musd_km_yr"], 2),
              "two rows", round(annual(55.0, 1.75, lvl, True)["net_musd_km_yr"], 2))
