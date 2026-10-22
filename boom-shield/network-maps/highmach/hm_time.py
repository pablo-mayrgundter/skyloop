#!/usr/bin/env python3
"""Higher-Mach study: door-to-door time by component, and the infinite-Mach asymptote.

Run:  python3 hm_time.py   -> door_to_door.csv, d2d_sensitivity.csv, day_return.csv

Block-time model = trunk/build_trunk.py (routes.md constants; imported read-only) decomposed into phases, on
each pair's own great-circle shield, Boomless feeder start (M1.15 after climb).  Access + processing + egress
= demand.md's overhead O (base 130 min, split 40 / 55 / 35 [U]).  Totals are checked against bt.t_trunk.
"""
from __future__ import annotations

import math
import sys
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
NM = HERE.parent
sys.dont_write_bytecode = True
sys.path.insert(0, str(NM / "trunk"))
sys.path.insert(0, str(NM / "demand"))
import build_trunk as bt  # noqa: E402

import hm_aircraft as ha  # noqa: E402

PAIRS = {"NYC-LAX": ((40.64, -73.78), (33.94, -118.41), -3),
         "DEN-NYC": ((39.86, -104.67), (40.64, -73.78), +2),
         "CHI-LAX": ((41.98, -87.90), (33.94, -118.41), -2)}
OVERHEAD = {"base (130 min)": (40, 55, 35), "fast (100 min)": (35, 40, 25), "walk-on (70 min)": (30, 20, 20)}
MACHS = [1.7, 2.0, 2.2, 2.5, 3.0, 3.5, 4.0]


def gc_km(p, q):
    return float(bt.gc(p[0], p[1], q[0], q[1]))


def accel_for(M, kT=None):
    if M <= 1.7 + 1e-9:
        return 0.23, 0.46, 20.25
    eng = "adaptive" if M <= 3.05 else "tbcc"
    kT = kT or (1.40 if eng == "adaptive" else 1.60)
    ac = ha.size(M, eng, kT=kT)
    a = ha.ramp_accel(ac)
    return a, min(2 * a, 0.9), ac["h_cruise_km"]


def phases_conv(d_gc, mode):
    L = d_gc * bt.PATH
    if mode == "sub":
        cruise = (L - 490) / bt.V085 * 60
        return dict(taxi=bt.TAXI, climb=18.0, accel=0.0, cruise=cruise, decel=0.0, descent=18.0, approach=10.0, merge=0.0)
    cruise = (L - 695) / bt.V115 * 60
    return dict(taxi=bt.TAXI, climb=18.0 + 7.0, accel=0.0, cruise=cruise, decel=0.0, descent=18.0 + 4.0, approach=10.0, merge=0.0)


def phases_shield(d_gc, M, a, b, h_km=20.25, merge=5.0):
    """Own-GC shield, Boomless feeder start (as bt.t_trunk(0, L, 0, M, 'bl'): a straight shield, no 1.03
    path factor, as coverage.md's own-great-circle upper bound). M = inf -> bang-bang."""
    L = d_gc
    s1 = bt.CLIMB[1] + bt.ACC_BL[1]
    s_fdesc = L - bt.APPR[1] - (bt.DESC_SUB[1] + bt.DEC_BL[1])
    A = (s_fdesc - s1 - bt.FL600_EXTRA[1]) * 1e3
    v0 = bt.V_BL_START
    k = 1 / (2 * a) + 1 / (2 * b)
    VM = M * bt.A_MS if np.isfinite(M) else 1e9
    need = (VM ** 2 - v0 ** 2) * k
    Vp = VM if A >= need else math.sqrt(v0 ** 2 + A / k)
    cruise_m = max(A - (Vp ** 2 - v0 ** 2) * k, 0)
    extra_desc = 0.5 * max(h_km - 20.25, 0)          # min, descent from above FL600+ [U]
    return dict(taxi=bt.TAXI, climb=18.0 + 7.0, accel=(Vp - v0) / a / 60, cruise=cruise_m / Vp / 60,
                decel=(Vp - v0) / b / 60, descent=bt.FL600_EXTRA[0] + 18.0 + 4.0 + extra_desc, approach=10.0,
                merge=merge, peak_mach=Vp / bt.A_MS, ramp_km=(Vp ** 2 - v0 ** 2) * k / 1e3)


def floor_phases(d_gc):
    """Zero cruise time, instantaneous ramps: what is left at infinite Mach and infinite thrust."""
    return dict(taxi=bt.TAXI, climb=18.0 + 7.0, accel=0.0, cruise=0.0, decel=0.0, descent=bt.FL600_EXTRA[0] + 22.0,
                approach=10.0, merge=5.0)


BLOCK = ["taxi", "climb", "accel", "cruise", "decel", "descent", "approach", "merge"]


def total(ph, O):
    return sum(ph[k] for k in BLOCK) + sum(O)


def day_ok(D_h, dz, W=5.0, Dmax=4.0):
    """demand.md day-return rule: leave home 06:00, home by 23:00; >= W h inside 08:30-18:00 at destination."""
    start = max(6.0 + D_h + dz, 8.5)
    end = min(23.0 - D_h + dz, 18.0)
    return (end - start) >= W and D_h <= Dmax


if __name__ == "__main__":
    pd.set_option("display.width", 250)
    rows = []
    for pair, (p, q, dz) in PAIRS.items():
        d = gc_km(p, q)
        modes = [("subsonic", phases_conv(d, "sub"), None), ("Boomless M1.15", phases_conv(d, "bl"), None)]
        for M in MACHS:
            a, b, h = accel_for(M)
            ph = phases_shield(d, M, a, b, h)
            # check against the repo's t_trunk (same constants)
            chk = float(bt.t_trunk(0.0, d, 0.0, M, "bl", 5.0, a, b)[0]) + 0.5 * max(h - 20.25, 0)
            assert abs(sum(ph[k] for k in BLOCK) - chk) < 0.05, (pair, M, sum(ph[k] for k in BLOCK), chk)
            modes.append((f"shield M{M}", ph, (a, b)))
        a4, b4, h4 = accel_for(4.0)
        modes.append(("shield M-inf, M4 thrust", phases_shield(d, np.inf, a4, b4, h4), (a4, b4)))
        modes.append(("shield M-inf, Overture thrust", phases_shield(d, np.inf, 0.23, 0.46), (0.23, 0.46)))
        modes.append(("floor: no cruise, no ramps", floor_phases(d), None))
        for name, ph, ab in modes:
            for oname, O in OVERHEAD.items():
                tot = total(ph, O)
                D_h = tot / 60
                rows.append(dict(pair=pair, gc_km=round(d), mode=name, overhead=oname, access=O[0], processing=O[1],
                                 egress=O[2], **{k: ph[k] for k in BLOCK}, block_min=sum(ph[k] for k in BLOCK),
                                 door_to_door_h=D_h, peak_mach=ph.get("peak_mach", np.nan),
                                 ramp_km=ph.get("ramp_km", np.nan), a_acc=ab[0] if ab else np.nan,
                                 day_return_out=day_ok(D_h, dz), day_return_back=day_ok(D_h, -dz),
                                 day_return_both_high=day_ok(D_h, dz, 4.0, 4.5) and day_ok(D_h, -dz, 4.0, 4.5)))
    T = pd.DataFrame(rows)
    T.to_csv(HERE / "door_to_door.csv", index=False, float_format="%.3f")
    print(T[T.overhead == "base (130 min)"][["pair", "mode", "taxi", "climb", "accel", "cruise", "decel", "descent",
                                              "block_min", "door_to_door_h", "peak_mach", "a_acc",
                                              "day_return_out", "day_return_back"]].round(2).to_string())

    # ---------------- thrust and ground-processing sensitivity, NYC-LAX
    p, q, dz = PAIRS["NYC-LAX"]
    d = gc_km(p, q)
    sens = []
    for M in MACHS + [np.inf]:
        Mm = 4.0 if not np.isfinite(M) else M
        a0, b0, h = accel_for(Mm)
        for lab, f in (("Overture thrust (0.23 m/s2)", None), ("concept thrust", 1.0), ("2x concept thrust", 2.0)):
            a, b = (0.23, 0.46) if f is None else (a0 * f, min(2 * a0 * f, 1.5))
            ph = phases_shield(d, M, a, b, h)
            for oname, O in OVERHEAD.items():
                sens.append(dict(M=M, thrust=lab, a_acc=a, overhead=oname, block_min=sum(ph[k] for k in BLOCK),
                                 d2d_h=total(ph, O) / 60, peak_mach=ph["peak_mach"]))
    SN = pd.DataFrame(sens)
    SN.to_csv(HERE / "d2d_sensitivity.csv", index=False, float_format="%.3f")
    print(SN.pivot_table(index="M", columns=["thrust", "overhead"], values="d2d_h").round(2).to_string())
