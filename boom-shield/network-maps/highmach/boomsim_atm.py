"""Thin read-only wrapper around sim/boomsim atmosphere and ray gain (US Std 1976)."""
from __future__ import annotations

import sys
from functools import lru_cache
from pathlib import Path

SIM = Path(__file__).resolve().parents[2] / "sim"
sys.dont_write_bytecode = True
if str(SIM) not in sys.path:
    sys.path.insert(0, str(SIM))
from boomsim.atmosphere import std_atm  # noqa: E402
from boomsim import propagate as PR  # noqa: E402


def atm(z):
    a = std_atm(float(z))
    return dict(p=float(a.p), rho=float(a.rho), T=float(a.T), c=float(a.c))


def pressure_alt(p):
    lo, hi = 0.0, 40e3
    for _ in range(60):
        mid = 0.5 * (lo + hi)
        if atm(mid)["p"] > p:
            lo = mid
        else:
            hi = mid
    return 0.5 * (lo + hi)


@lru_cache(maxsize=None)
def ray_gain(H_s, V, h_ac=250.0):
    """psf_sweep's R(H): undertrack linear amplitude ratio shield -> ground, x 1.9 ground reflection."""
    return 1.9 * float(PR.RayPath(H_s + h_ac, V, z_min=0.0).amplitude_ratio(H_s, 0.0))
