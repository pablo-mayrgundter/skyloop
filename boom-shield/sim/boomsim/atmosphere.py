"""US Standard Atmosphere 1976 (0-86 km) plus an assumed humidity profile.

Public API
----------
std_atm(z) -> Atm            vectorised; z = geometric altitude [m]
molar_humidity(z, ...) -> h  molar concentration of water vapour [%] (ISO 9613-1 sense)
saturation_ratio(T)          p_sat/p_r over liquid water, ISO 9613-1 Eq. (B.1)/(B.2)
RH_PROFILE                   default (z [m], RH [%]) table used by molar_humidity

Standard atmosphere
-------------------
U.S. Standard Atmosphere 1976 (NOAA-S/T 76-1562): geopotential altitude
H = r0 z / (r0 + z) with r0 = 6 356 766 m; piecewise-linear molecular-scale
temperature in H with the standard layer bases and lapse rates (0, 11, 20, 32, 47,
51, 71 km; -6.5, 0, +1.0, +2.8, 0, -2.8, -2.0 K/km) up to 84.852 km geopotential;
hydrostatic pressure with g0 = 9.80665 m/s^2, R* = 8.31432 J/(mol K),
M0 = 28.9644 g/mol; sound speed c = sqrt(gamma R T) with gamma = 1.4.
(Above 86 km geometric the 1976 model changes form; it is never needed here.)

Humidity
--------
The relative-humidity profile is an *assumption* (the ANSI S1.26-2014 Annex C
profile used for the AIAA Sonic Boom Prediction Workshop optional cases was not
accessible while writing this; the default below follows its qualitative shape:
~60 % at the ground rising to ~70 % in the lower/mid troposphere, falling through
the upper troposphere to a few per cent in the lower stratosphere; above 20 km it
is tapered further so the implied water-vapour mixing ratio stays small). It is converted to a
molar concentration with the ISO 9613-1 Annex B formula

    h = h_r (p_sat/p_r) / (p_a/p_r),
    log10(p_sat/p_r) = -6.8346 (T01/T)^1.261 + 4.6151,  T01 = 273.16 K.

The ISO saturation formula is used outside its nominal validity range in the cold
stratosphere; there the humidity is tiny either way (h ~ 1e-3 %), which only moves
the O2 relaxation frequency, whose contribution is negligible at altitude.
Override with ``molar_humidity(z, rh_profile=(z_tab, rh_tab))``.
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np

# --- constants (US Std Atmosphere 1976) --------------------------------------
G0 = 9.80665              # m/s^2
R_STAR = 8.31432          # J/(mol K)  (1976 value)
M0 = 28.9644e-3           # kg/mol
R_AIR = R_STAR / M0       # 287.0531 J/(kg K)
GAMMA = 1.4
R_EARTH = 6356766.0       # m, effective earth radius for geopotential
P_SL = 101325.0           # Pa
T_SL = 288.15             # K

# geopotential base altitudes [m] and lapse rates [K/m]
_H_BASE = np.array([0.0, 11000.0, 20000.0, 32000.0, 47000.0, 51000.0, 71000.0, 84852.0])
_LAPSE = np.array([-6.5e-3, 0.0, 1.0e-3, 2.8e-3, 0.0, -2.8e-3, -2.0e-3])


def _layer_bases():
    T = [T_SL]
    p = [P_SL]
    for i in range(len(_LAPSE)):
        dH = _H_BASE[i + 1] - _H_BASE[i]
        Tb, pb, L = T[-1], p[-1], _LAPSE[i]
        if L == 0.0:
            T.append(Tb)
            p.append(pb * np.exp(-G0 * M0 * dH / (R_STAR * Tb)))
        else:
            Tn = Tb + L * dH
            T.append(Tn)
            p.append(pb * (Tb / Tn) ** (G0 * M0 / (R_STAR * L)))
    return np.array(T), np.array(p)


_T_BASE, _P_BASE = _layer_bases()


@dataclass
class Atm:
    """Ambient state. Fields are floats or arrays matching the input ``z``."""
    z: np.ndarray
    T: np.ndarray
    p: np.ndarray
    rho: np.ndarray
    c: np.ndarray


def geopotential(z):
    z = np.asarray(z, dtype=float)
    return R_EARTH * z / (R_EARTH + z)


def std_atm(z) -> Atm:
    """US Standard Atmosphere 1976 at geometric altitude(s) ``z`` [m] (0 ... 86 km).

    Scalars in give 0-d numpy floats out; arrays give arrays of the same shape.
    Slightly negative z (down to -1 km) is extrapolated with the tropospheric layer.
    """
    z = np.asarray(z, dtype=float)
    if np.any(z > 86000.0) or np.any(z < -1000.0):
        raise ValueError("std_atm valid for -1 km <= z <= 86 km")
    H = geopotential(z)
    i = np.clip(np.searchsorted(_H_BASE, H, side="right") - 1, 0, len(_LAPSE) - 1)
    Hb, Tb, pb, L = _H_BASE[i], _T_BASE[i], _P_BASE[i], _LAPSE[i]
    T = Tb + L * (H - Hb)
    with np.errstate(divide="ignore", invalid="ignore"):
        p_grad = pb * (Tb / T) ** (G0 * M0 / (R_STAR * np.where(L == 0, 1.0, L)))
    p_iso = pb * np.exp(-G0 * M0 * (H - Hb) / (R_STAR * Tb))
    p = np.where(L == 0.0, p_iso, p_grad)
    rho = p / (R_AIR * T)
    c = np.sqrt(GAMMA * R_AIR * T)
    return Atm(z=z, T=T, p=p, rho=rho, c=c)


# --- humidity -----------------------------------------------------------------
# Assumed relative-humidity profile (z [m], RH [%]); see module docstring.
RH_PROFILE = (
    np.array([0.0, 1000.0, 2000.0, 3000.0, 5000.0, 7000.0, 9000.0, 11000.0,
              13000.0, 16000.0, 20000.0, 30000.0, 47000.0, 86000.0]),
    np.array([60.0, 63.0, 66.0, 69.0, 69.0, 60.0, 45.0, 25.0,
              10.0, 5.0, 3.0, 1.0, 0.05, 0.05]),
)

P_REF = 101325.0       # ISO 9613-1 reference pressure
T01 = 273.16           # triple point (ISO 9613-1 Annex B)


def saturation_ratio(T):
    """p_sat/p_r over water, ISO 9613-1 Annex B: 10^(-6.8346 (T01/T)^1.261 + 4.6151)."""
    T = np.asarray(T, dtype=float)
    return 10.0 ** (-6.8346 * (T01 / T) ** 1.261 + 4.6151)


def relative_humidity(z, rh_profile=None):
    """Relative humidity [%] at altitude z from the (piecewise-linear) profile."""
    zt, rh = RH_PROFILE if rh_profile is None else rh_profile
    return np.interp(np.asarray(z, dtype=float), zt, rh)


def molar_humidity(z, rh_profile=None, atm=None):
    """Molar concentration of water vapour h [%] at altitude z.

    h = RH * (p_sat/p_r) / (p/p_r) (ISO 9613-1 Annex B), with RH from
    ``rh_profile`` (default ``RH_PROFILE``) and T, p from ``atm`` (default std_atm(z)).
    """
    a = std_atm(z) if atm is None else atm
    return relative_humidity(z, rh_profile) * saturation_ratio(a.T) / (a.p / P_REF)
