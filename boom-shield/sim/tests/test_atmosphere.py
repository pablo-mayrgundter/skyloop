import os
import sys

import numpy as np
import pytest

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))

from boomsim.atmosphere import (R_AIR, molar_humidity, relative_humidity,  # noqa: E402
                                saturation_ratio, std_atm)

# US Standard Atmosphere 1976 tabulated values (geometric altitude)
TABLE = [
    # z [m], T [K], p [Pa], rho [kg/m^3], c [m/s]
    (0.0, 288.150, 101325.0, 1.2250, 340.294),
    (5000.0, 255.676, 54048.0, 0.73643, 320.545),
    (11000.0, 216.774, 22700.0, 0.36480, 295.154),
    (16800.0, 216.650, 9131.0, 0.14683, 295.070),   # ~ interpolated table
    (20000.0, 216.650, 5529.3, 0.088910, 295.070),
    (32000.0, 228.490, 889.06, 0.013555, 303.025),
    (47000.0, 270.650, 115.85, 0.0014975, 329.799),
]


@pytest.mark.parametrize("z,T,p,rho,c", TABLE)
def test_std_atm_table(z, T, p, rho, c):
    a = std_atm(z)
    assert float(a.T) == pytest.approx(T, abs=1.2 if z == 47000.0 else 0.02)
    assert float(a.p) == pytest.approx(p, rel=2e-3)
    assert float(a.rho) == pytest.approx(rho, rel=5e-3)
    assert float(a.c) == pytest.approx(c, rel=2e-3)


def test_std_atm_vectorised_and_monotone():
    z = np.linspace(0, 47000, 471)
    a = std_atm(z)
    assert a.T.shape == z.shape and a.p.shape == z.shape
    assert np.all(np.diff(a.p) < 0)
    assert np.all(np.diff(a.rho) < 0)
    # ideal gas & sound speed consistency
    assert R_AIR == pytest.approx(287.0531, rel=1e-6)
    assert np.allclose(a.p / (a.rho * a.T), R_AIR, rtol=1e-12)
    assert np.allclose(a.c, np.sqrt(1.4 * R_AIR * a.T), rtol=1e-12)
    # isothermal layer 11-20 km geopotential
    assert np.allclose(std_atm([12000.0, 18000.0]).T, 216.65)


def test_saturation_and_humidity():
    # ISO 9613-1 saturation: ~2.34 kPa at 20 C, 611 Pa at 0.01 C (triple point)
    assert saturation_ratio(293.15) * 101325 == pytest.approx(2339, rel=0.01)
    assert saturation_ratio(273.16) * 101325 == pytest.approx(611.7, rel=0.01)
    # ground: RH 60 % at 15 C, 1 atm -> h ~ 0.6 * 1705/101325 * 100 ~ 1.01 %
    assert float(molar_humidity(0.0)) == pytest.approx(1.01, rel=0.03)
    z = np.array([0.0, 5000.0, 11000.0, 17000.0])
    h = molar_humidity(z)
    assert np.all(np.diff(h) < 0)
    assert h[-1] < 0.01                     # dry stratosphere
    assert float(molar_humidity(25000.0)) < 0.01
    z = np.append(z, 25000.0)
    rh = relative_humidity(z)
    assert 55 <= rh[0] <= 70 and rh[-1] < 5


def test_humidity_override():
    prof = (np.array([0.0, 50000.0]), np.array([100.0, 100.0]))
    h100 = molar_humidity(0.0, rh_profile=prof)
    assert float(h100) == pytest.approx(float(molar_humidity(0.0)) * 100 / 60, rel=1e-6)
