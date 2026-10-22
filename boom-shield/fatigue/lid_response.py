"""BoomShield lid fatigue, step 1: boom load and lid response per pass.

Reads the undertrack incident waveform at the shield plane (h = 250 m, shield top 20 km)
cached by sim/run_psf_sweep.py (nonlinear near field -> 20 km), and drives the 1-D
floating-lid transfer-matrix model (sim/boomsim/shield.py) with it.  From the top-face
reflection R and transmission T we recover, per pass:

  p_top   = (1 + R) P            pressure on the lid's upper surface (reflection doubling)
  v_top   = (1 - R) P / Z_a      lid normal velocity (Z_a = rho_a omega / k_z, normal impedance)
  p_gas   = p_top - j w m v_top  gas pressure just under the lid (net load = lid inertia)
  v_bot   = T P / Z_a            floor velocity;  e = w_bot - w_top  thickness change (ties)
  u_gas   = p_gas / (rho_g V)    horizontal gas velocity under the lid (evanescent field,
                                  k_x = omega / V); its time integral is the sideways
                                  excursion of a slack transverse partition near the lid

Lid strains (travelling pattern w(x - V t), so d/dx = -(1/V) d/dt):
  slope  = v_top / V;  curvature kappa = a_top / V^2
  bending strain    eps_b = kappa * t_lid / 2
  membrane strain   eps_m ~ slope^2 / 2  (upper bound: all the extra arc length taken
                    as stretch, i.e. in-plane sliding fully restrained)
Writes lid_response.csv, lid_waveforms.csv.
"""
from __future__ import annotations

import math
import os
import pickle
import sys

import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
BS = os.path.dirname(HERE)
sys.path.insert(0, os.path.join(BS, "sim"))
from boomsim import shield as S  # noqa: E402
from boomsim.atmosphere import std_atm  # noqa: E402

H_S = 20e3
FLOOR = 0.1
CASES = [
    # aircraft, Mach, H2 d [m], lid areal mass [kg/m2] (None = lift - floor), tie K [Pa/m]
    ("Overture", 1.7, 36.0, None, 10.0),
    ("Overture", 1.7, 36.0, None, 0.0),
    ("Concorde", 2.0, 55.0, None, 10.0),
    ("Concorde", 2.0, 55.0, None, 0.0),
    ("Concorde", 2.0, 55.0, 2.3, 10.0),
    ("Concorde", 2.0, 30.0, None, 10.0),
]
RHO_LID = 930.0        # PE, kg/m3 -> lid thickness t = m / rho
E_COLD = (1.0e9, 3.0e9)  # PE laminate modulus at 217 K, low/high [Pa] (assumed; ~0.2-0.5 GPa at 23 C)
TIE_PITCH = 10.0       # m, discrete soft-tie spacing (structure §5 tests 5-25 m)


def incident(ac, mach):
    key = f"above_{ac}_M{mach:g}_H{H_S / 1e3:g}"
    with open(os.path.join(BS, "sim", "results", "psf_cache", f"{key}.above.pkl"), "rb") as fh:
        d = pickle.load(fh)
    p = np.asarray(d["p"], float)
    t = d["t0"] + np.arange(p.size) * d["dt"]
    return t, p


def n_wave_metrics(t, p):
    pk, i_pk = p.max(), int(np.argmax(p))
    tr, i_tr = p.min(), int(np.argmin(p))
    big = np.nonzero(np.abs(p) > 0.05 * pk)[0]
    # duration: front shock (first exceedance) to rear shock (last big positive jump)
    dur = (t[i_tr] - t[big[0]]) + 0.0
    dp = np.diff(p)
    i_rear = i_tr + int(np.argmax(dp[i_tr:i_tr + int(0.05 / (t[1] - t[0]))]))
    dur = t[i_rear] - t[big[0]]
    pos = np.clip(p, 0, None).sum() * (t[1] - t[0])
    neg = np.clip(p, None, 0).sum() * (t[1] - t[0])
    return dict(p_inc_pos_pa=pk, p_inc_neg_pa=tr, duration_s=dur, impulse_pos_pa_s=pos,
                impulse_neg_pa_s=neg)


def response(ac, mach, d, lid, K):
    atm = std_atm(H_S)
    rho_a, c_a = float(atm.rho), float(atm.c)
    g = S.gas("hydrogen", float(atm.T), float(atm.p))
    lift = (rho_a - g.rho) * d
    m_lid = (lift - FLOOR) if lid is None else lid
    st = S.Stack([S.Membrane(m_lid), S.GasLayer(g, d, tie_K=K), S.Membrane(FLOOR)])
    V = mach * c_a
    t, p = incident(ac, mach)
    met = n_wave_metrics(t, p)
    dt = t[1] - t[0]
    pad = 8.0
    n0 = int(pad / dt)
    buf = np.concatenate([np.zeros(n0), p, np.zeros(n0)])
    n = 1 << int(np.ceil(np.log2(buf.size)))
    buf = np.concatenate([buf, np.zeros(n - buf.size)])
    P = np.fft.rfft(buf)
    f = np.fft.rfftfreq(n, dt)
    f[0] = 1e-6
    w = 2 * np.pi * f
    T, R = S.transfer(st, f, V, atm, reference="interface")
    kz = np.sqrt((w / c_a) ** 2 - (w / V) ** 2)
    Za = rho_a * w / kz                      # real, frequency independent
    jw = 1j * w
    spec = {
        "p_top": (1 + R) * P,
        "v_top": (1 - R) * P / Za,
        "p_gas": (1 + R) * P - jw * m_lid * (1 - R) * P / Za,
        "v_bot": T * P / Za,
    }
    spec["w_top"] = spec["v_top"] / jw
    spec["a_top"] = spec["v_top"] * jw
    spec["w_bot"] = spec["v_bot"] / jw
    spec["e"] = spec["w_bot"] - spec["w_top"]
    spec["u_gas"] = spec["p_gas"] / (g.rho * V)
    spec["x_gas"] = spec["u_gas"] / jw
    for k in ("w_top", "w_bot", "e", "x_gas"):
        spec[k][0] = 0.0                     # no static offset
    # low-pass the acceleration at 2 kHz (shock rise ~ ms; avoids FFT ringing dominating)
    lp = 1 / (1 + (f / 2000.0) ** 4)
    spec["a_top"] = spec["a_top"] * lp
    tt = (np.arange(n) - n0) * dt + t[0]
    ts = {k: np.fft.irfft(v, n) for k, v in spec.items()}
    win = (tt > t[0] - 2.0) & (tt < t[-1] + 3.0)
    ts = {k: v[win] for k, v in ts.items()}
    tt = tt[win]

    def amp(x):
        return float(np.max(np.abs(x)))

    t_lid = m_lid / RHO_LID
    slope = amp(ts["v_top"]) / V
    curv = amp(ts["a_top"]) / V ** 2
    eps_m = 0.5 * slope ** 2
    eps_b = curv * t_lid / 2
    e_pp = float(ts["e"].max() - ts["e"].min())
    row = dict(aircraft=ac, mach=mach, V=V, h2_m=d, lid_kg_m2=round(m_lid, 2), tie_K=K, **met,
               p_top_max_pa=amp(ts["p_top"]), p_top_pos_pa=float(ts["p_top"].max()),
               p_top_neg_pa=float(ts["p_top"].min()),
               p_gas_max_pa=amp(ts["p_gas"]), net_lid_load_pa=amp(ts["p_top"] - ts["p_gas"]),
               v_top_max=amp(ts["v_top"]), w_top_pp_m=float(ts["w_top"].max() - ts["w_top"].min()),
               w_bot_pp_m=float(ts["w_bot"].max() - ts["w_bot"].min()),
               thickness_change_pp_m=e_pp,
               a_top_max=amp(ts["a_top"]), slope=slope, curvature_1_m=curv,
               t_lid_mm=t_lid * 1e3, eps_membrane=eps_m, eps_bending=eps_b,
               stress_lid_kpa_lo=(eps_m + eps_b) * E_COLD[0] / 1e3,
               stress_lid_kpa_hi=(eps_m + eps_b) * E_COLD[1] / 1e3,
               tension_change_n_m_hi=(eps_m + eps_b) * E_COLD[1] * t_lid,
               tie_force_pp_n=K * e_pp * TIE_PITCH ** 2,
               u_gas_max=amp(ts["u_gas"]), partition_excursion_pp_m=float(ts["x_gas"].max() - ts["x_gas"].min()),
               lift_kg_m2=lift)
    return row, tt, ts


def main():
    rows, wf = [], []
    for ac, M, d, lid, K in CASES:
        r, tt, ts = response(ac, M, d, lid, K)
        rows.append(r)
        print({k: (round(v, 4) if isinstance(v, float) else v) for k, v in r.items()})
        if lid is None and K == 10.0 and d in (36.0, 55.0):
            sel = slice(None, None, 25)      # 2 kHz output
            wf.append(pd.DataFrame(dict(case=f"{ac} M{M:g} H2 {d:g} m", t_s=tt[sel],
                                        p_inc=np.interp(tt[sel], *incident(ac, M), left=0, right=0),
                                        p_top=ts["p_top"][sel], p_gas=ts["p_gas"][sel],
                                        w_top=ts["w_top"][sel], w_bot=ts["w_bot"][sel], e=ts["e"][sel],
                                        x_gas=ts["x_gas"][sel])))
    pd.DataFrame(rows).to_csv(os.path.join(HERE, "lid_response.csv"), index=False, float_format="%.4g")
    pd.concat(wf).to_csv(os.path.join(HERE, "lid_waveforms.csv"), index=False, float_format="%.5g")


if __name__ == "__main__":
    main()
