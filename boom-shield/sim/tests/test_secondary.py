"""Tests for boomsim.secondary (upper atmosphere, absorption, rays, march, binning)."""
import dataclasses
import math

import numpy as np
import pytest

from boomsim import secondary as S
from boomsim.atmosphere import std_atm
from boomsim.propagate import WindProfile, ground_signature, iso_absorption
from boomsim.source import CONCORDE, SR71

AC = dataclasses.replace(CONCORDE, mach=2.0, altitude=20250.0)


# ------------------------------------------------------------------ atmosphere
@pytest.mark.parametrize("zk", sorted(S.USSA76_TABLE))
def test_upper_atmosphere_matches_ussa76_table(zk):
    T, p, rho, M = S.USSA76_TABLE[zk]
    st = S.upper_atm().state(np.array(zk * 1e3))
    assert float(st["T"]) == pytest.approx(T, abs=0.1)
    assert float(st["p"]) == pytest.approx(p, rel=0.015)
    assert float(st["rho"]) == pytest.approx(rho, rel=0.015)


def test_upper_atmosphere_equals_std_atm_below_86km():
    ua = S.upper_atm()
    z = np.array([0.0, 11e3, 20e3, 47e3, 80e3])
    a = std_atm(z)
    st = ua.state(z)
    assert np.allclose(st["T"], a.T, rtol=1e-6)
    assert np.allclose(st["p"], a.p, rtol=1e-4)
    assert np.allclose(st["c"], a.c, rtol=1e-6)


def test_ussa76_temperature_segments():
    T = S.ussa76_temperature_above86(np.array([88e3, 100e3, 115e3, 120e3, 1000e3]))
    assert T[0] == pytest.approx(186.8673)
    assert T[1] == pytest.approx(195.08, abs=0.01)
    assert T[2] == pytest.approx(300.0)
    assert T[3] == pytest.approx(360.0)
    assert T[4] == pytest.approx(1000.0, abs=1.0)


def test_sound_speed_rises_in_thermosphere():
    ua = S.upper_atm()
    c = ua.interp("c", np.array([100e3, 120e3, 150e3]))
    assert c[0] < 290 and 400 < c[1] < 410 and 560 < c[2] < 575


# ------------------------------------------------------------------ absorption
def test_classical_absorption_matches_iso_at_ground():
    st = S.upper_atm().state(np.array(0.0))
    a = S.classical_coefficient(st)
    # ISO 9613-1 classical+rotational: 1.84e-11 (p_r/p) (T/293.15)^0.5
    iso = 1.84e-11 * (101325.0 / float(st["p"])) * (float(st["T"]) / 293.15) ** 0.5
    assert float(a) == pytest.approx(iso, rel=0.02)


def test_absorption_exponent_total_close_to_iso_low_altitude():
    st = S.upper_atm().state(np.array(1000.0))
    f = np.array([100.0, 1000.0, 5000.0])
    ln = S.absorption_exponent(f, st, 1.0)
    iso = iso_absorption(f, st["T"], st["p"], st["h"])
    assert np.allclose(-ln.real, iso, rtol=0.03)


def test_absorption_scales_inversely_with_pressure():
    ua = S.upper_atm()
    a100 = float(S.classical_coefficient(ua.state(np.array(100e3))))
    a120 = float(S.classical_coefficient(ua.state(np.array(120e3))))
    # alpha ~ mu/(p c): p drops 12.6x, T/c rise compensate a little
    assert 8 < a120 / a100 < 20


# ------------------------------------------------------------------ rays
def test_launch_slowness_on_mach_cone_and_dispersion_relation():
    med = S.Medium(WindProfile.uniform(20.0, 5.0))
    phi = np.radians([0, 45, 90, 135, 180])
    sx, sy, sz = S.launch_slowness(AC.V, AC.altitude, phi, med)
    c, _, u, _, v, _ = med.fields(np.full(phi.size, AC.altitude))
    s = np.sqrt(sx ** 2 + sy ** 2 + sz ** 2)
    assert np.allclose(c * s + u * sx + v * sy, 1.0, atol=1e-12)
    # ground speed: crab into the cross wind, V_g = sqrt(V^2 - v^2) + u
    Vg = math.sqrt(AC.V ** 2 - 25.0) + 20.0
    assert np.allclose(sx, 1.0 / Vg, rtol=1e-12)
    assert sz[0] < 0 and sz[-1] > 0


def _uniform_medium(c0=300.0):
    ua = S.UpperAtmosphere(z_max=60e3, dz=50.0)
    ua.c = np.full_like(ua.c, c0)
    return S.Medium(None, ua=ua)


def test_straight_rays_and_cylindrical_tube_in_uniform_air():
    med = _uniform_medium()
    path = S.ray_path(AC, math.radians(150.0), med, "direct", t_max=60.0, z_top=50e3)
    assert path is None          # straight ray never returns
    tr = S.trace([math.radians(150.0)], AC.V, AC.altitude, med, dt=0.5, t_max=30, record=True,
                 z_top=50e3)
    R = tr["rec"][:, :, 0]
    d = np.diff(R[:, :3], axis=0)
    dirs = d / np.linalg.norm(d, axis=1)[:, None]
    assert np.allclose(dirs, dirs[0], atol=1e-9)
    # ray speed = c
    assert np.allclose(np.linalg.norm(d, axis=1) / 0.5, 300.0, rtol=1e-9)


def test_tube_flux_area_grows_linearly_in_uniform_air():
    med = _uniform_medium()
    dphi = 1e-4
    tr = S.trace(np.radians(150.0) + np.array([-dphi, 0.0, dphi]), AC.V, AC.altitude, med,
                 dt=0.5, t_max=40, record=True, z_top=55e3)
    R = tr["rec"]
    rp = (R[:, 0:3, 2] - R[:, 0:3, 0]) / (2 * dphi)
    v = np.diff(R[:, 0:3, 1], axis=0)[-1] / 0.5
    Q = -rp[:, 2] * v[1] + rp[:, 1] * v[2]
    t = np.arange(R.shape[0]) * 0.5
    ratio = Q[10:] / t[10:]
    assert np.allclose(ratio, ratio[0], rtol=1e-6)


def test_turning_height_matches_snell():
    med = S.Medium()
    phi = math.radians(135.0)
    sx, sy, _ = S.launch_slowness(AC.V, AC.altitude, [phi], med)
    v_tr = 1.0 / math.hypot(sx[0], sy[0])
    tr = S.trace([phi], AC.V, AC.altitude, med, dt=0.2, max_ground=1)
    apex = tr["apex"][0, 0]
    z = med.ua.z
    zt = z[(z > 90e3) & (med.c >= v_tr)][0]
    assert apex == pytest.approx(zt, abs=300.0)


def test_calm_has_no_stratospheric_returns_and_thermospheric_carpet():
    med = S.Medium()
    tr, rows = S.fan(AC, med, np.radians(np.arange(90.0, 150.0, 2.0)), dt=0.5)
    direct = [r for r in rows if r["family"] == "direct"]
    assert direct and all(r["apex"] > 100e3 for r in direct)
    # rays between ~90 and ~120 deg are trapped (no landing) in the calm USSA76
    assert np.sum(tr["status"] == 2) > 5


def test_tailwind_creates_stratospheric_returns():
    med = S.Medium(S.upper_wind("winter", 0.0))
    tr, rows = S.fan(AC, med, np.radians(np.arange(90.0, 130.0, 1.0)), dt=0.5)
    strat = [r for r in rows if r["family"] == "direct" and r["apex"] < 70e3]
    assert len(strat) > 3
    assert all(40e3 < r["apex"] < 70e3 for r in strat)


def test_shield_mirror_and_family_classification():
    med = S.Medium()
    sh = S.Shield(z=20000.0, phi_max=math.atan2(500, 250), mode="reflect")
    tr, rows = S.fan(AC, med, np.radians([30.0, 45.0, 70.0]), shield=sh, dt=0.5)
    fams = {round(math.degrees(r["phi"])): r["family"] for r in rows}
    assert fams.get(30) == "shield" and fams.get(45) == "shield"
    assert tr["kind"][2, 0] == 0 or tr["n_hits"][2] == 0 or fams.get(70) != "shield"


# ------------------------------------------------------------------ waveform
def test_hilbert_shift_twice_inverts():
    t = np.linspace(-1, 1, 2048)
    p = np.exp(-(t / 0.05) ** 2) * np.sin(40 * t)
    q = S.hilbert_shift(S.hilbert_shift(p))
    assert np.allclose(q, -p, atol=1e-6 * np.max(np.abs(p)) + 1e-3 * np.max(np.abs(p)))


def test_primary_march_matches_ground_signature():
    ac = dataclasses.replace(SR71, mach=3.2, altitude=20250.0)
    t1, p1 = ground_signature(ac, fs=20e3)
    t2, p2 = S.primary_signature(ac, 0.0, fs=20e3, n_span=10 ** 7, ds_lin=20.0)
    assert np.max(p2) == pytest.approx(np.max(p1), rel=0.02)
    assert np.min(p2) == pytest.approx(np.min(p1), rel=0.03)


def test_shield_operator_energy_split():
    from boomsim import shield as SH
    a = std_atm(20000.0)
    g = SH.gas("hydrogen", float(a.T), float(a.p))
    st = SH.Stack([SH.Membrane(1.24), SH.GasLayer(g, 30.0), SH.Membrane(1.24)])
    t = np.arange(4096) / 4000.0
    p = np.exp(-((t - 0.5) / 0.02) ** 2)
    _, pr = S.shield_operator(st, a, 590.0, "R")(t, p)
    _, pt = S.shield_operator(st, a, 590.0, "T")(t, p)
    e = np.sum(p ** 2)
    assert (np.sum(pr ** 2) + np.sum(pt ** 2)) / e == pytest.approx(1.0, rel=1e-3)


# ------------------------------------------------------------------ binning
class _R:
    def __init__(self, E):
        self.E_h = E
        self.t = np.arange(256) * 1e-2
        self.h = np.zeros(256)
        self.h[100] = math.sqrt(E / 1e-2)


def test_binning_conserves_energy():
    med = S.Medium()
    dphi = math.radians(0.1)
    phis = np.radians(np.arange(130.0, 140.0, 0.1))
    rows = [dict(phi=q, Y=200e3 + 5e3 * (q - phis[0]) / dphi * 0.1, apex=115e3, vz=-150.0,
                 sx=1 / 590.0, sy=1e-3, t=1000.0, X=1e5) for q in phis]
    marched = {phis[0]: _R(2.0), phis[-1]: _R(2.0)}
    bins = S.bin_ground(rows, marched, med, dphi, dY=3000.0)
    tot = sum(b["energy"] for b in bins) * 3000.0 / S.GROUND_FACTOR ** 2
    st = S.upper_atm().state(np.array(0.0))
    c0, _, u, _, v, _ = med.fields(np.zeros(1))
    G = float(st["rho"]) * float(st["c"]) ** 2 / 150.0
    assert tot == pytest.approx(G * 2.0 * dphi * phis.size, rel=1e-9)


def test_branches_split_at_fold():
    phi = np.radians(np.arange(0, 10, 0.5))
    Y = (phi - phi[10]) ** 2
    lab = S.branches(phi, Y, np.full(phi.size, 110e3))
    assert len(np.unique(lab)) == 2
