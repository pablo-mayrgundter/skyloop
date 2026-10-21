"""BoomShield infographic data driver.

  python boom-shield/infographic/make_data.py --run [--procs 3] [--cache DIR]
      Boom Overture (M1.7) preset through the run_3d.py pipeline (same numerics:
      R_m 3 km all azimuths, near-track re-run at 8 km; f_max 100 Hz, df 1/6 Hz,
      50 kHz ground propagation), unshielded and W 1 km H2 30 m slab, h 250 m.
      Writes <cache>/overture_summary.json. Does not touch results/3d_cache.
  python boom-shield/infographic/make_data.py --build [--cache DIR]
      Assembles boom-shield/infographic/data.json from the repo results, the Overture
      summary and the cited reference tables below.
"""
from __future__ import annotations

import argparse
import dataclasses
import json
import math
import os
import sys
import tempfile

HERE = os.path.dirname(os.path.abspath(__file__))
SIM = os.path.join(os.path.dirname(HERE), "sim")
sys.path.insert(0, SIM)

import run_3d as R3  # noqa: E402
from boomsim.source import Aircraft  # noqa: E402

DEFAULT_CACHE = os.environ.get("OVERTURE_CACHE", os.path.join(tempfile.gettempdir(), "boomshield_overture_cache"))

# Boom Overture preset (estimates, see data_notes.md):
#  length 201 ft = 61.3 m (Boom spec sheet); cruise M 1.7 (Boom);
#  MTOW ~77 t (pre-2022 published figure, unverified for the refined design);
#  mid-cruise mass ~0.81 MTOW (same fraction as the repo's Concorde preset 150/185 t)
#  -> 62 t -> W = 0.61 MN; volume ~300 m^3 (Concorde 450 m^3 scaled for a ~2.7 m
#  fuselage and ~40 t fuel); lift on the wing, ~0.40 L ... 0.90 L.
OVERTURE = Aircraft("Overture", length=61.3, mach=1.7, altitude=18300.0,
                    weight=62000.0 * 9.81, volume=300.0,
                    lift_start=0.40, lift_end=0.90)

_aircraft0 = R3.aircraft


def aircraft(name, mach, h):
    if name == "Overture":
        return dataclasses.replace(OVERTURE, mach=float(mach), altitude=R3.H_S + float(h))
    return _aircraft0(name, mach, h)


R3.aircraft = aircraft


def overture_configs():
    base = dict(group="baseline", acn="Overture", mach=1.7, h=250, W=None, fill=None, T=0.0,
                profile="slab", fringe=0.0, cells=0.0, tieK=0.0, ties="line", R=R3.R_MAIN,
                phis=None)
    C = [dict(base, key="O1.7_h250_none", label="Overture M1.7 h250: unshielded")]
    sh = dict(base, group="aircraft", W=1000.0, fill="hydrogen", T=30.0,
              key="O1.7_h250_W1000_H230", label="Overture M1.7: W 1000 m, H2 30 m")
    C.append(sh)
    C.append(dict(sh, key=sh["key"] + "_R8000", group="far", R=R3.R_FAR, phis=R3.PHI_SENS_FAR,
                  label=sh["label"] + " [R_m 8 km]"))
    for c in C:
        c["base"] = "O1.7_h250_none"
    return C


def run(procs, cache):
    R3.CACHE = cache
    R3.WORKERS = procs
    C = overture_configs()
    for c in C:
        print(f"== {c['key']}", flush=True)
        R3.ground(c, procs)
    summ = summarize_overture(cache)
    print(json.dumps(summ, indent=1))


def _clean(v):
    if isinstance(v, float) and (math.isnan(v) or math.isinf(v)):
        return None
    if hasattr(v, "item"):
        return _clean(v.item())
    if isinstance(v, dict):
        return {str(k): _clean(x) for k, x in v.items()}
    return v


def summarize_overture(cache):
    R3.CACHE = cache
    C = {c["key"]: c for c in overture_configs()}
    base = R3.carpet_rows(C["O1.7_h250_none"])
    out = {}
    for k in ("O1.7_h250_none", "O1.7_h250_W1000_H230"):
        rows = R3.carpet_rows(C[k])
        if rows is None:
            continue
        s = R3.summarize(C[k], rows, base if k != "O1.7_h250_none" else None)
        s["per_phi"] = [dict(phi=r["phi"], y_g_km=r["y_g"] / 1e3, pldb=R3._m(r, "pldb"),
                             pldb_proxy=R3._m(r, "pldb_proxy")) for r in rows]
        out[k] = {kk: _clean(v) for kk, v in s.items()}
    ac = aircraft("Overture", 1.7, 250)
    out["preset"] = dict(name=ac.name, length_m=ac.length, mach=ac.mach, altitude_m=ac.altitude,
                         weight_N=ac.weight, volume_m3=ac.volume, lift_start=ac.lift_start,
                         lift_end=ac.lift_end)
    with open(os.path.join(cache, "overture_summary.json"), "w") as fh:
        json.dump(out, fh, indent=1)
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--run", action="store_true")
    ap.add_argument("--build", action="store_true")
    ap.add_argument("--procs", type=int, default=3)
    ap.add_argument("--cache", default=DEFAULT_CACHE)
    a = ap.parse_args()
    os.makedirs(a.cache, exist_ok=True)
    if a.run:
        run(a.procs, a.cache)
    if a.build:
        build(a.cache)


if __name__ == "__main__":
    main()
