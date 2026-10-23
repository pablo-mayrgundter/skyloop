#!/usr/bin/env python3
"""Skyloop internal-flight design map (screening level).

Winged pods fly (aerodynamic lift) over an actively steered guideway / motor strip that supplies thrust,
inside a wide (20-80 m) low-pressure bore of an H2-lifted, tensioned tube at ~20 km, with ground portals.
Maps every discipline over Mach (2-12) x tube pressure (5-1,000 Pa) x bore (20/40/80 m) x pod class
(3.5-47 t, 5-150 seats) x wing loading, and marks the binding constraint at each point.

Run:   python3 skyloop-review/flight/flight.py              (~3-6 min, one process; stages resume)
       python3 skyloop-review/flight/flight.py figs         (figures only, from the CSVs)
Out:   flight_stroke.csv     guideway stroke needed vs Mach, bore, line tension, bridle residual (stability.py ride)
       flight_grid.csv       best wing loading per (pod, bore, Mach, pressure): state, margins, status, binding
                             (> 2 MB, gitignored)
       flight_band.csv       pressure band per pod x bore x Mach, and highest feasible Mach
       flight_ceilings.csv   speed ceilings vs Mach: trip length x comfort accel, curvature, power, KE, LSM
       flight_gasdyn.csv     shock angle, wall-reflection distance, Kantrowitz, wake Pi vs Mach x bore
       flight_points.csv     best operating points (highest Mach per bore x trip, energy- and capacity-optimal)
       flight_econ.csv       B/C hook for the best point per Mach x bore (spine+LAX and NYC-LAX direct)
       flight_levers.csv     what widens the band
       flight_maps.png, flight_summary.png

Reuses read-only: heating/heating.py (pod_state, tube_side, sink_fraction, thresholds),
stability/stability.py (line, ride, tension_mass, tube_props, economics), regime/regime.py (pod, q_hole,
beta_kantrowitz, float_cost), system/sky_model.py (headway_s, ride_min, run_cost, conventions) and
system/*.csv (lifecycle, breakeven, networks).  Mach convention: V = M x 295.07 m/s.  (A) = assumption.
"""
from __future__ import annotations

import csv
import math
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
REV = os.path.dirname(HERE)
sys.dont_write_bytecode = True
for _d in ("heating", "stability"):
    sys.path.insert(0, os.path.join(REV, _d))
import heating as H  # noqa: E402  (read-only)
import stability as S  # noqa: E402  (read-only; imports regime, physics, sky_model)

R, sp, sm = S.R, S.sp, S.sm
G, R_E = 9.81, 6.371e6
C_REF = H.C_REF
RHO_AMB, P_AMB, T_AMB = S.RHO, S.P_AMB, S.T_AMB
LF = sm.P["lf"]
PEAK_FACTOR = 2 / sm.P["peak_share"] * 365            # pax/yr per tube pair per (pax/h per direction)

# ----------------------------------------------------------------------------- grids
M_GRID = np.round(np.arange(2.0, 12.001, 0.25), 2)
P_GRID = np.geomspace(5.0, 1000.0, 30)
D_LIST = [20.0, 40.0, 80.0]
PODS = [3500.0, 8000.0, 17000.0, 30000.0, 47000.0]   # 5 / 20 / 50 / 93 / 150 seats (regime.pod, 300 kg/seat)
WS_LIST = [50.0, 100.0, 200.0, 400.0, 800.0, 1600.0]
M_STROKE = np.round(np.arange(2.0, 12.001, 0.5), 2)
TRIPS = [1000e3, 2500e3, 4000e3]
ACCS = [0.1 * G, 0.2 * G, 0.3 * G]

# ----------------------------------------------------------------------------- base assumptions (A)
BASE = dict(
    T_line=0.1e9,        # N, line tension of the "tensioned envelope" (stability: 0.04-0.37 GN for M1-4)
    eps=S.EPS_BRIDLE,    # bridle residual (stability: 10 %)
    clear=1.0,           # m, wing tip / body to wall clearance at full stroke
    AR=2.0,              # wing aspect ratio (heating: delta, sweep 63 deg)
    rule="moving block 10 s",
    a_acc=0.3 * G,       # launch / brake (system, physics)
    L_trip=4000e3,       # NYC-LAX
    f_mix=1.0,           # gas-to-wall mixing multiplier (heating)
    T_lim_le=1700.0,     # C/SiC (heating)
    leak=1.0,            # x (1 cm2/km at D 3 m, scaled with circumference; stability bigbore)
    imp=3.0,             # strip heating amplification by shock impingement (heating caveat: 3-10x)
    g_mag=0.15,          # m, magnetic gap a flying pod holds over the strip
    tau_max=2.0,         # m, largest LSM pole pitch (iron/copper per pole) ; also <= L_pod / 6
    sig0=40e3,           # Pa, LSM thrust shear at zero gap (A; long-stator LSM class)
    w_strip=2.0,         # m, strip width (heating)
    f_ok=1000.0, f_marg=3000.0,          # Hz, LSM fundamental (Transrapid ~270 Hz)
    P_pod_ok=200e6, P_pod_marg=800e6,    # W, peak traction power per pod (EMALS: 484 MJ in 2-3 s ~ 0.2 GW)
    P_stn_ok=1e9, P_stn_marg=5e9,        # W, launch power per station (physics: 0.6-5 GW)
    g_breach_ok=5.0, g_breach_marg=20.0, # g, meeting 20 km ambient in a breach (regime)
    E_ok=340.0, E_marg=680.0,            # Wh/pax-km: jet-fuel parity / 2x (energy is costed in the B/C hook)
    N_ok=120.0, N_marg=30.0,             # pods/h per bore the tube gas must absorb (heating)
    demand=81e6, capture=1.0,
    eta_d=0.85, eta_r=0.85,              # system
    hotel_w_seat=400.0,                  # system
    act_kw_km=2.6,                       # actuated guideway power (stability: 2-3 MW per 1,000 km)
)

CHECKS = ["lift", "stroke", "trip", "power", "motor", "breach", "blockage", "reflection", "crossover",
          "gas_dT", "throughput", "skin_T", "le_T", "heat_sink", "film_T", "guideway_T", "energy"]
GROUP = {"lift": "flight", "stroke": "flight", "trip": "ceiling", "power": "ceiling", "motor": "ceiling",
         "breach": "ceiling", "blockage": "gas dynamics", "reflection": "gas dynamics", "crossover": "gas dynamics",
         "gas_dT": "gas dynamics", "throughput": "heating", "skin_T": "heating", "le_T": "heating",
         "heat_sink": "heating", "film_T": "heating", "guideway_T": "heating", "energy": "energy"}
SHORT = {"lift": "lift", "stroke": "steering stroke", "trip": "trip length", "power": "launch power",
         "motor": "motor strip", "breach": "breach decel", "blockage": "Kantrowitz", "reflection": "shock on pod",
         "crossover": "crossovers", "gas_dT": "gas ΔT / wake", "throughput": "gas heat → pods/h",
         "skin_T": "skin T", "le_T": "nose/LE T", "heat_sink": "heat sink", "film_T": "film T",
         "guideway_T": "strip T", "energy": "energy"}


# ----------------------------------------------------------------------------- io
def write_csv(name, rows, mode="w"):
    path = os.path.join(HERE, name)
    keys = []
    for r in rows:
        keys += [k for k in r if k not in keys]
    new = mode == "w" or not os.path.exists(path)
    with open(path, mode, newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=keys, extrasaction="ignore")
        if new:
            w.writeheader()
        for r in rows:
            w.writerow({k: (f"{v:.4g}" if isinstance(v, float) else v) for k, v in r.items()})
    print(f"wrote {name} ({len(rows)} rows)", flush=True)


def read_csv(name):
    path = os.path.join(HERE, name)
    if not os.path.exists(path):
        return []
    out = []
    for r in csv.DictReader(open(path)):
        rr = {}
        for k, v in r.items():
            try:
                rr[k] = float(v)
            except (TypeError, ValueError):
                rr[k] = v
        out.append(rr)
    return out


# ----------------------------------------------------------------------------- guideway stroke (stability.ride)
_STROKE = {}


def stroke_table(T, eps):
    """Stroke (4 sigma + scallop) the steered guideway needs, lateral and vertical, vs Mach, per bore."""
    key = (round(T / 1e6, 3), round(eps, 4))
    if key in _STROKE:
        return _STROKE[key]
    tab = {}
    rows = []
    for D in D_LIST:
        L = S.line(D=D, payload_extra=20.0 + S.tension_mass(T), eps=eps, T=T)
        sy, sz, apk = [], [], []
        for M in M_STROKE:
            rd, _ = S.ride(L, float(M), "decoupled", 1e-3)          # tiny s_avail -> least stroke meeting the ride
            sy.append(rd["y"]["stroke"])
            sz.append(rd["z"]["stroke"])
            apk.append(max(rd["y"]["a_pk"], rd["z"]["a_pk"]))
            rows.append(dict(T_GN=T / 1e9, eps=eps, D_bore=D, Mach=float(M), stroke_y_m=sy[-1], stroke_z_m=sz[-1],
                             a_pk_g=apk[-1] / G, env_R_m=L["Re"], sig_y_tube_m=rd["y"]["sig_y"],
                             tension_kg_m=S.tension_mass(T)))
        tab[D] = dict(y=np.array(sy), z=np.array(sz), R=L["Re"], L=L)
    _STROKE[key] = (tab, rows)
    return _STROKE[key]


def s_req(tab, D, M):
    return float(np.interp(M, M_STROKE, tab[D]["y"])), float(np.interp(M, M_STROKE, tab[D]["z"]))


# ----------------------------------------------------------------------------- one point
def hd(F, m):
    """heating.py parameter dict for this pod and lever set."""
    return dict(H.BASE, m=m, AR=F["AR"], L_trip=F["L_trip"], a_acc=F["a_acc"], f_mix=F["f_mix"],
                T_lim_le=F["T_lim_le"], N_req=F["N_ok"], N_marg=F["N_marg"])


def _st(mg, tol, k, m, t):
    mg[k], tol[k] = m, t


def evaluate(st, D, m, F, tab, d):
    pd_ = R.pod(m)
    n, pax = pd_["n"], pd_["n"] * LF
    r = dict(pod_t=m / 1e3, seats=n, D_bore=D, Mach=st["Mach"], p_Pa=st["p_Pa"], WS_Pa=st["WS_Pa"],
             q_Pa=st["q_Pa"], span_m=st["span_m"], S_m2=st["S_m2"], relief=st["relief"])
    mg, tol = {}, {}
    if not st.get("lift_ok"):
        for k in CHECKS:
            _st(mg, tol, k, -9.0 if k == "lift" else 0.0, 0.3)
        r.update(alpha_deg=90.0)
        return finish(r, mg, tol)
    V, M, p = st["V"], st["Mach"], st["p_Pa"]
    Dr = st["drag_N"]
    A_b = math.pi * D * D / 4
    r.update(alpha_deg=st["alpha_deg"], CL=st["CL"], LD=st["LD"], drag_kN=Dr / 1e3)
    # 1. lift (heating: lower-surface incidence; curvature relief inside pod_state)
    _st(mg, tol, "lift", math.log10(d["alpha_le_ok"] / st["alpha_deg"]), math.log10(d["alpha_le_marg"] / d["alpha_le_ok"]))
    # 2. steering stroke vs span (stability ride/stroke; pod follows the steered strip)
    sy, sz = s_req(tab, D, M)
    b = st["span_m"]
    av_y = (D - b) / 2 - F["clear"]
    av_z = (D - (pd_["d"] + 1.5)) / 2 - F["clear"]                  # body + strip depth (A)
    m_y = math.log10(max(av_y, 1e-3) / sy)
    m_z = math.log10(max(av_z, 1e-3) / sz)
    _st(mg, tol, "stroke", min(m_y, m_z), math.log10(1.5))
    r.update(stroke_req_y_m=sy, stroke_av_y_m=av_y, stroke_req_z_m=sz)
    # 3. trip length: the cruise Mach must be reachable at a_acc (time-optimal profile = bang-bang at the cap)
    _st(mg, tol, "trip", math.log10(F["a_acc"] * F["L_trip"] / V ** 2), 1e-6)
    # capacity first (it sets traffic for power, pumping, strip and film)
    h_rule = sm.headway_s(M, F["rule"])
    ts0 = H.tube_side(st, D, d, N=F["N_ok"])
    N_heat = ts0["N_max_per_h"]
    N_max = min(3600.0 / h_rule, N_heat)
    cap = N_max * n * LF * PEAK_FACTOR
    dem = F["demand"] * F["capture"]
    n_pairs = max(1.0, dem / cap)
    N_dem = dem / (n * LF * PEAK_FACTOR)                            # pods/h per direction for the whole demand
    N_ops = min(N_max, N_dem / n_pairs)
    r.update(headway_s=h_rule, N_heat_h=N_heat, N_max_h=N_max, cap_pair_Mpax=cap / 1e6, n_pairs=n_pairs,
             N_ops_h=N_ops)
    # 4. launch power: per pod (peak, at the end of the launch) and per station (launches at N_ops)
    KE = 0.5 * m * V * V
    P_pod = (Dr + m * F["a_acc"]) * V / F["eta_d"]
    P_stn = KE / F["eta_d"] * N_ops / 3600.0
    z1 = math.log10(F["P_pod_ok"] / P_pod), math.log10(F["P_pod_marg"] / F["P_pod_ok"])
    z2 = math.log10(F["P_stn_ok"] / P_stn), math.log10(F["P_stn_marg"] / F["P_stn_ok"])
    _st(mg, tol, "power", *(z1 if z1[0] / z1[1] < z2[0] / z2[1] else z2))
    r.update(P_pod_MW=P_pod / 1e6, P_station_GW=P_stn / 1e9, KE_GJ=KE / 1e9, KE_tTNT=KE / 4.184e9)
    # 5. motor strip: LSM fundamental frequency and thrust shear across the magnetic gap
    tau = min(F["tau_max"], pd_["L"] / 6)
    f = V / (2 * tau)
    sig_req = (Dr + m * F["a_acc"]) / (F["w_strip"] * 0.8 * pd_["L"])
    sig_av = F["sig0"] * math.exp(-math.pi * F["g_mag"] / tau)
    z1 = math.log10(F["f_ok"] / f), math.log10(F["f_marg"] / F["f_ok"])
    z2 = math.log10(sig_av / sig_req), math.log10(2)
    _st(mg, tol, "motor", *(z1 if z1[0] / z1[1] < z2[0] / z2[1] else z2))
    r.update(lsm_Hz=f, tau_m=tau, shear_req_kPa=sig_req / 1e3, shear_av_kPa=sig_av / 1e3)
    # 6. breach: drag deceleration in 20 km ambient with the wing unloaded; lift overload if it is not
    q_amb = 0.5 * RHO_AMB * V * V
    a_br = q_amb * (0.5 * pd_["A"] + 0.006 * st["S_m2"]) / (m * G)
    _st(mg, tol, "breach", math.log10(F["g_breach_ok"] / a_br), math.log10(F["g_breach_marg"] / F["g_breach_ok"]))
    r.update(breach_decel_g=a_br, breach_lift_overload=P_AMB / p, breach_q_kPa=q_amb / 1e3,
             breach_nose_MWm2=H.q_stag(RHO_AMB, V, H.freestream(M, p)[3], d["R_n"], 1000.0) / 1e6)
    # 7. Kantrowitz blockage (frontal: body + wing thickness 4 % of chord + incidence)
    A_fr = pd_["A"] + b * 0.04 * st["chord_m"] + st["S_m2"] * math.sin(math.radians(st["alpha_deg"]))
    beta = A_fr / A_b
    bK = R.beta_kantrowitz(M)
    _st(mg, tol, "blockage", math.log10(0.5 * bK / beta), math.log10(2))
    r.update(beta=beta, beta_kantrowitz=bK)
    # 8. wall-reflected bow shock returning onto the pod (typical offset = half the 4-sigma stroke)
    d_w = max(0.5, (D - b) / 2 - sy / 2)
    x_r = 2 * d_w * math.sqrt(M * M - 1)
    L_ref = max(pd_["L"], 2 * st["chord_m"])
    _st(mg, tol, "reflection", math.log10(x_r / L_ref), math.log10(2))
    r.update(mach_angle_deg=math.degrees(math.asin(1 / M)), x_reflect_m=x_r, L_pod_m=pd_["L"])
    # 9. self-pumping: crossover spacing so that pods do not push gas into closed sections (regime)
    L_seg = p * A_b * V * h_rule / Dr / 1e3
    _st(mg, tol, "crossover", math.log10(L_seg / 1.0), math.log10(4.0))
    r.update(crossover_km=L_seg)
    # 10. tube gas per pass = wake seen by the next pod (Pi = D/pA) (heating)
    _st(mg, tol, "gas_dT", math.log10(d["dT_pass_ok"] / ts0["dT_pass"]), math.log10(d["dT_pass_marg"] / d["dT_pass_ok"]))
    _st(mg, tol, "throughput", math.log10(N_heat / F["N_ok"]), math.log10(F["N_ok"] / F["N_marg"]))
    r.update(Pi=ts0["Pi"], dT_pass_K=ts0["dT_pass"], tau_relax_s=ts0["tau_relax_s"])
    # 11. pod heating (heating.evaluate conventions)
    Tsk = max(st["T_wing_lower"], st["T_body"])
    _st(mg, tol, "skin_T", max(4 * math.log10(d["T_lim_skin"] / Tsk), -0.99 * math.log10(2)), math.log10(2))
    Tle = max(st["T_nose"], st["T_le"])
    _st(mg, tol, "le_T", max(4 * math.log10(d["T_lim_le"] / Tle), -0.99 * math.log10(2)), math.log10(2))
    fs, _, _ = H.sink_fraction(st, d)
    if math.isnan(fs):
        _st(mg, tol, "heat_sink", 0.0, 0.3)
    else:
        _st(mg, tol, "heat_sink", math.log10(d["f_sink_ok"] / max(fs, 1e-6)), math.log10(d["f_sink_marg"] / d["f_sink_ok"]))
    r.update(T_skin_K=Tsk, T_le_K=Tle, T0_K=st["T0"], f_sink=fs, Q_in_kW=st["Q_in_kW"])
    # 12. film and strip at operating traffic, strip with shock-impingement amplification
    ts = H.tube_side(st, D, d, N=max(N_ops, 1e-3))
    E_conv = st["q_wing_lower_kWm2"] * 1e3 * 2 * st["chord_m"] / V
    q_strip = N_ops / 3600 * (ts["E_strip_pass_Jm2"] + (F["imp"] - 1) * E_conv)
    T_strip = (H.T_ENV ** 4 + q_strip / (2 * H.EPS * H.SIGMA)) ** 0.25
    Tf = max(ts["T_film_sun"], ts["T_film_tube"])
    _st(mg, tol, "film_T", 4 * math.log10(d["T_film_lim"] / Tf), 4 * math.log10(d["T_film_marg"] / d["T_film_lim"]))
    _st(mg, tol, "guideway_T", 4 * math.log10(d["T_strip_lim"] / T_strip), 4 * math.log10(453.0 / d["T_strip_lim"]))
    r.update(T_strip_K=T_strip, T_film_K=Tf)
    # 13. energy per pax-km at L_trip: kinetic loss + drag through the strip + hotel + pumping + actuators
    L = F["L_trip"]
    e_ke = KE * (1 / F["eta_d"] - F["eta_d"] * F["eta_r"]) / (pax * L) / 3.6
    e_dr = Dr / F["eta_d"] / pax / 3.6
    e_ht = F["hotel_w_seat"] * n / pax / V / 3.6
    leak = F["leak"] * D / 3.0                                       # cm2 per km (stability bigbore)
    P_pump = 10e3 * leak * R.q_hole(P_AMB, T_AMB, 250.0) / p          # W per km (regime: 10 kW per m3/s)
    paxh = max(N_ops * pax, 1e-6)
    e_pu = P_pump / paxh
    e_ac = F["act_kw_km"] * 1e3 / paxh
    E = e_ke + e_dr + e_ht + e_pu + e_ac
    _st(mg, tol, "energy", math.log10(F["E_ok"] / E), math.log10(F["E_marg"] / F["E_ok"]))
    r.update(E_Wh_paxkm=E, E_kinetic=e_ke, E_drag=e_dr, E_pump=e_pu, E_other=e_ht + e_ac, pump_kW_km=P_pump / 1e3,
             d2d_NYC_LAX_h=100 / 60 + float(sm.ride_min(3974.4 * sm.P["path"], M, F["a_acc"])) / 60)
    return finish(r, mg, tol)


def finish(r, mg, tol):
    stt, worst, wz, nf = "pass", None, 1e9, 0
    for k in CHECKS:
        s = "pass" if mg[k] >= 0 else ("marginal" if mg[k] >= -tol[k] else "fail")
        r[f"m_{k}"] = round(mg[k], 4)
        if s == "fail":
            stt, nf = "fail", nf + 1
        elif s == "marginal" and stt == "pass":
            stt = "marginal"
        z = mg[k] / tol[k]
        if z < wz:
            wz, worst = z, k
    # without the trip check (ceilings for other trip lengths are applied afterwards)
    snt, wnt, znt = "pass", None, 1e9
    for k in CHECKS:
        if k == "trip":
            continue
        s = "pass" if mg[k] >= 0 else ("marginal" if mg[k] >= -tol[k] else "fail")
        if s == "fail":
            snt = "fail"
        elif s == "marginal" and snt == "pass":
            snt = "marginal"
        if mg[k] / tol[k] < znt:
            znt, wnt = mg[k] / tol[k], k
    r.update(status=stt, binding=worst, score=wz, nfail=nf, status_nt=snt, binding_nt=wnt, score_nt=znt)
    return r


RANK = {"pass": 2, "marginal": 1, "fail": 0}


def better(a, b, key="status"):
    """True if row a beats row b (status, fewer fails, larger worst margin)."""
    if b is None:
        return True
    sk = "score" if key == "status" else "score_nt"
    return (RANK[a[key]], -a["nfail"], a[sk]) > (RANK[b[key]], -b["nfail"], b[sk])


def run_grid(F, pods=PODS, d_list=D_LIST, m_grid=M_GRID, p_grid=P_GRID, ws_list=WS_LIST, key="status",
             out=None, done=()):
    tab, _ = stroke_table(F["T_line"], F["eps"])
    rows = []
    for m in pods:
        if m in done:
            continue
        d = hd(F, m)
        prow = []
        for M in m_grid:
            for p in p_grid:
                best = {D: None for D in d_list}
                for WS in ws_list:
                    st = H.pod_state(float(M), float(p), WS, d)
                    for D in d_list:
                        r = evaluate(st, D, m, F, tab, d)
                        if better(r, best[D], key):
                            best[D] = r
                prow += [best[D] for D in d_list]
        rows += prow
        if out:
            write_csv(out, prow, mode="a")
        print(f"  pod {m/1e3:g} t: {sum(r['status'] == 'pass' for r in prow)} pass, "
              f"{sum(r['status'] == 'marginal' for r in prow)} marginal of {len(prow)}", flush=True)
    return rows


# ----------------------------------------------------------------------------- bands and points
def bands(rows, key="status"):
    out = []
    bk = "binding" if key == "status" else "binding_nt"
    for m in sorted(set(r["pod_t"] for r in rows)):
        for D in sorted(set(r["D_bore"] for r in rows)):
            sub = [r for r in rows if r["pod_t"] == m and r["D_bore"] == D]
            for M in sorted(set(r["Mach"] for r in sub)):
                s2 = sorted([r for r in sub if r["Mach"] == M], key=lambda r: r["p_Pa"])
                ok = [i for i, r in enumerate(s2) if r[key] == "pass"]
                mk = [i for i, r in enumerate(s2) if r[key] in ("pass", "marginal")]
                row = dict(pod_t=m, seats=s2[0]["seats"], D_bore=D, Mach=M)
                if ok:
                    lo, hi = ok[0], ok[-1]
                    row.update(p_pass_lo=s2[lo]["p_Pa"], p_pass_hi=s2[hi]["p_Pa"],
                               floor_by=s2[lo - 1][bk] if lo > 0 else "grid edge",
                               ceiling_by=s2[hi + 1][bk] if hi < len(s2) - 1 else "grid edge")
                else:
                    best = max(s2, key=lambda r: (RANK[r[key]], r["score"]))
                    row.update(p_pass_lo=float("nan"), p_pass_hi=float("nan"), floor_by=best[bk], ceiling_by="")
                row.update(p_marg_lo=s2[mk[0]]["p_Pa"] if mk else float("nan"),
                           p_marg_hi=s2[mk[-1]]["p_Pa"] if mk else float("nan"))
                out.append(row)
    return out


def max_mach(rows, key="status", ok=("pass",)):
    Ms = [r["Mach"] for r in rows if r[key] in ok]
    return max(Ms) if Ms else float("nan")


def top_binding(rows, Mtop, key="status"):
    """Binding constraint of the best cell one Mach step above Mtop."""
    bk = "binding" if key == "status" else "binding_nt"
    if math.isnan(Mtop):
        cand = rows
    else:
        nxt = [r["Mach"] for r in rows if r["Mach"] > Mtop + 1e-9]
        if not nxt:
            return "grid edge"
        cand = [r for r in rows if r["Mach"] == min(nxt)]
    best = max(cand, key=lambda r: (RANK[r[key]], -r["nfail"], r["score"]))
    return best[bk]


# ----------------------------------------------------------------------------- economics hook
SYS = os.path.join(REV, "system")
NETS = {"spine+LAX all-stop": 4074.0, "NYC-LAX direct": 4045.0}
CMP = {"sub": 7.291, "bl": 6.27, "m22": 5.08}                       # NYC-LAX door to door, h (system d2d.csv)
LSM_EFF, GUIDE_EFF = R.LSM_EFF, R.GUIDE_EFF


def _sysdata():
    import pandas as pd
    L = pd.read_csv(os.path.join(SYS, "lifecycle.csv"))
    B = pd.read_csv(os.path.join(SYS, "breakeven.csv"))
    N = pd.read_csv(os.path.join(SYS, "networks.csv"))
    L = L[L.opening == 2050]
    B = B[B.opening == 2050]
    out = {}
    k_run = sm.run_cost("low") / sm.run_cost("base")
    for net in NETS:
        for vs in CMP:
            for cap in (0.3, 1.0):
                for M in (2.0, 3.0, 4.0, 6.0):
                    lb = L[(L.network == net) & (L.vs == vs) & (L.capture == cap) & (L.M == M) & (L.cost == "base")].iloc[0]
                    ll = L[(L.network == net) & (L.vs == vs) & (L.capture == cap) & (L.M == M) & (L.cost == "low")].iloc[0]
                    bg = B[(B.network == net) & (B.vs == vs) & (B.capture == cap) & (B.M == M) & (B.kind == "gross")].iloc[0]
                    nr = N[(N.network == net) & (N.vs == vs) & (N.M == M)].iloc[0]
                    k = bg.scale_vs_base
                    C1 = lb.pv_cost_busd
                    C0 = (lb.pv_ben_busd - k * C1) / (1 - k)
                    dC = C1 - C0                                               # PV of the tube capex at base
                    C0l = ll.pv_cost_busd - dC * k_run
                    out[(net, vs, cap, M)] = dict(ben=lb.pv_ben_busd, en=lb.pv_energy_busd, C0=C0, dC=dC, C0l=C0l,
                                                  c_base=bg.base_avg_tube_usd_m_km, paxkm=nr.pax_km_B,
                                                  wh_sys=nr.wh_per_pax_km, run=nr.run_km, ramp=nr.ramp_km,
                                                  bc_g=lb.bc_gross, bc_n=lb.bc_net)
    return out


_SD = None


def sysd():
    global _SD
    if _SD is None:
        _SD = _sysdata()
    return _SD


def d2d_h(M, a=0.3 * G):
    return 100 / 60 + float(sm.ride_min(3974.4 * sm.P["path"], M, a)) / 60


_ADD = {}


def adders(D, T, eps):
    """$M/km added to the floating run+ramp by bore wall, tension, envelope/tether growth and actuators
    (stability.economics conventions), at base and low cost."""
    key = (D, round(T / 1e6, 3), round(eps, 4))
    if key not in _ADD:
        tab, _ = stroke_table(T, eps)
        L = tab[D]["L"]
        add_kg = {"tension": S.tension_mass(T), "wall": max(0.0, S.tube_props(D)[0] - S.tube_props(S.D_BASE)[0])}
        _, _, cb, cl = S.economics(4.0, L, add_kg, S.ACT_USD_M_KM, 0.0, 0.0, D)
        f4 = R.float_cost(4.0)
        _ADD[key] = (cb - f4, cl - f4 * 0.28)
    return _ADD[key]


def tube_cost(net, M, D, F):
    """Average $M/km of a tube pair for this network at Mach M (run + ramps + flight adders), base and low."""
    km = NETS[net]
    if net == "NYC-LAX direct":
        ramp = 2 * min(km / 2, max(sm.accel_km(M), sm.P["ramp_min_km"]))
    else:
        sd = sysd()
        fr = np.interp(min(M, 6.0), [2, 3, 4, 6], [sd[(net, "sub", 1.0, x)]["ramp"] / km for x in (2.0, 3.0, 4.0, 6.0)])
        ramp = km * (1.0 if M >= 6 else fr)
    run = km - ramp
    a_b, a_l = adders(D, F["T_line"], F["eps"])
    out = {}
    for lvl, add in (("base", a_b), ("low", a_l)):
        i = sm.LVL[lvl]
        c_run = sm.run_cost(lvl)
        c_ramp = c_run * sm.RAMP_MULT[i] + sm.LSM[i] * (1 + sm.SOFT[i]) * (1 + sm.CONTING[i])
        strip = LSM_EFF[lvl] - GUIDE_EFF[lvl]                          # full-length motor strip, no levitation track
        out[lvl] = (run * (c_run + strip) + ramp * c_ramp) / km + add
    return out


def bc_hook(r, net, F):
    """B/C (2050 opening) for one operating point: gross/net x 30/100 % capture x comparator x cost level."""
    sd = sysd()
    M = r["Mach"]
    Mx = min(M, 6.0)
    tc = tube_cost(net, M, r["D_bore"], F)
    cap_pair = r["cap_pair_Mpax"] * 1e6
    out = dict(net=net, tube_base_usd_m_km=tc["base"], tube_low_usd_m_km=tc["low"])
    for vs in CMP:
        for cap in (0.3, 1.0):
            v = {k: float(np.interp(Mx, [2, 3, 4, 6], [sd[(net, vs, cap, x)][k] for x in (2.0, 3.0, 4.0, 6.0)]))
                 for k in ("ben", "en", "C0", "dC", "C0l", "c_base", "paxkm", "wh_sys")}
            ben = v["ben"]
            e_sys = v["wh_sys"]
            if net == "NYC-LAX direct" and M > 6:
                sv6, svM = CMP[vs] - d2d_h(6.0), CMP[vs] - d2d_h(M)
                ben *= svM / sv6
                e_sys = float(sm.kwh_per_pax(4045, M)[0]) * 1e3 / 4045
            # extra electricity of flying vs the system's maglev pods: drag + pumping + actuators - EDS drag
            e_mag = sm.P["kg_per_seat"] / LF * G / sm.P["LD_mag"] / sm.P["eta_drive"] / 3.6
            dE = r["E_drag"] + r["E_pump"] + r["E_other"] - e_mag
            if net == "NYC-LAX direct" and M > 6:
                dE += e_sys - v["wh_sys"]
            dPV = dE * v["paxkm"] * cap * R.USD_KWH * R.PV / 1e3          # $B (no traffic growth: low side)
            pairs = max(1.0, F["demand"] * cap / cap_pair)
            k_run = sm.run_cost("low") / sm.run_cost("base")
            for lvl in ("base", "low"):
                C0 = v["C0"] if lvl == "base" else v["C0l"]
                dC = v["dC"] * (1.0 if lvl == "base" else k_run)
                cb = v["c_base"] * (1.0 if lvl == "base" else k_run)
                C = C0 + dC * tc[lvl] * pairs / cb
                tag = f"{vs}_{int(cap * 100)}_{lvl}"
                out[f"bcg_{tag}"] = ben / C
                out[f"bcn_{tag}"] = (ben + v["en"] - dPV) / C
            out[f"pairs_{int(cap * 100)}"] = pairs
    return out


# ----------------------------------------------------------------------------- stages
def stage_stroke():
    rows = []
    for T in (0.01e9, 0.1e9, 0.3e9, 1.0e9):
        rows += stroke_table(T, S.EPS_BRIDLE)[1]
    rows += stroke_table(0.1e9, 0.03)[1]
    write_csv("flight_stroke.csv", rows)


def stage_grid(F=BASE):
    path = os.path.join(HERE, "flight_grid.csv")
    done = set()
    if os.path.exists(path) and "--fresh" not in sys.argv:
        done = {r["pod_t"] * 1e3 for r in read_csv("flight_grid.csv")}
        done = {m for m in PODS if any(abs(m - x) < 1 for x in done)}
    elif os.path.exists(path):
        os.remove(path)
    run_grid(F, out="flight_grid.csv", done=done)
    return read_csv("flight_grid.csv")


def stage_band(rows):
    write_csv("flight_band.csv", bands(rows))


def stage_ceilings():
    rows = []
    tab, _ = stroke_table(BASE["T_line"], BASE["eps"])
    for M in np.arange(2.0, 12.01, 1.0):
        V = M * C_REF
        r = dict(Mach=M, V_ms=V, felt_weight=1 - V * V / (G * R_E),
                 R_min_curve_km=V * V / (G * math.sqrt(1.1 ** 2 - 1)) / 1e3,       # coordinated turn, n <= 1.1
                 transition_km=V * math.sqrt(1.1 ** 2 - 1) * G / (0.05 * G) / 1e3,   # 0.05 g/s roll-in jerk
                 launch_km_03g=V * V / (2 * 0.3 * G) / 1e3,
                 lsm_Hz_tau2m=V / 4.0)
        for m in (3500.0, 17000.0, 47000.0):
            r[f"KE_tTNT_{m/1e3:g}t"] = 0.5 * m * V * V / 4.184e9
            r[f"P_launch_MW_{m/1e3:g}t"] = m * 0.3 * G * V / 0.85 / 1e6
        for a in ACCS:
            r[f"d2d_NYC_LAX_h_a{a/G:.1f}g"] = d2d_h(M, a)
            for L in TRIPS:
                r[f"reach_{L/1e6:g}Mm_a{a/G:.1f}g"] = int(V * V <= a * L)
        for D in D_LIST:
            r[f"stroke_y_D{D:g}"] = s_req(tab, D, M)[0]
        rows.append(r)
    # time-optimal Mach cap per trip length x acceleration
    for a in ACCS:
        for L in TRIPS:
            Vm = math.sqrt(a * L)
            rows.append(dict(Mach=Vm / C_REF, kind="trip cap", a_g=a / G, trip_km=L / 1e3,
                             t_min_bangbang_min=2 * math.sqrt(L / a) / 60,
                             t_at_M4_min=float(sm.ride_min(L / 1e3, 4.0, a)),
                             t_at_M6_min=float(sm.ride_min(L / 1e3, 6.0, a))))
    write_csv("flight_ceilings.csv", rows)


def stage_gasdyn():
    rows = []
    for M in (2.0, 3.0, 4.0, 6.0, 8.0, 10.0, 12.0):
        for D in D_LIST:
            bK = R.beta_kantrowitz(M)
            pd_ = R.pod(17000.0)
            rows.append(dict(Mach=M, D_bore=D, mach_angle_deg=math.degrees(math.asin(1 / M)),
                             x_wall_hit_centre_m=D / 2 * math.sqrt(M * M - 1),
                             x_reflect_back_centre_m=D * math.sqrt(M * M - 1),
                             x_reflect_over_D=math.sqrt(M * M - 1), beta_kantrowitz=bK,
                             beta_17t_body=pd_["A"] / (math.pi * D * D / 4),
                             wave_freq_Hz=M * C_REF / (D * math.sqrt(M * M - 1)),
                             wake_du_over_V_Pi025=0.25 / (1.4 * M * M)))
    write_csv("flight_gasdyn.csv", rows)


def stage_points(rows):
    out = []
    ok = [r for r in rows if r["status_nt"] == "pass"]
    okm = [r for r in rows if r["status_nt"] in ("pass", "marginal")]
    # highest feasible Mach per bore x trip length (0.3 g), any pod, and the binding constraint above it
    for D in D_LIST:
        for L in TRIPS:
            Mcap = math.sqrt(0.3 * G * L) / C_REF
            for lvl, pool in (("pass", ok), ("marginal", okm)):
                c = [r for r in pool if r["D_bore"] == D and r["Mach"] <= Mcap + 1e-9]
                if not c:
                    out.append(dict(kind=f"highest Mach ({lvl})", D_bore=D, trip_km=L / 1e3, Mach=float("nan")))
                    continue
                Mt = max(r["Mach"] for r in c)
                cc = [r for r in c if r["Mach"] == Mt]
                b = max(cc, key=lambda r: (r["cap_pair_Mpax"], -r["E_Wh_paxkm"]))
                allD = [r for r in rows if r["D_bore"] == D]
                above = "trip length" if Mt + 0.25 > Mcap else top_binding(allD, Mt, "status_nt")
                out.append(dict(kind=f"highest Mach ({lvl})", trip_km=L / 1e3, bound_above=above, **pick(b)))
    for D in D_LIST + [None]:
        c = [r for r in ok if D is None or r["D_bore"] == D]
        if not c:
            continue
        e = min(c, key=lambda r: r["E_Wh_paxkm"])
        out.append(dict(kind="energy-optimal", trip_km=4000.0, **pick(e)))
        k = max(c, key=lambda r: (r["cap_pair_Mpax"], -r["E_Wh_paxkm"]))
        out.append(dict(kind="capacity-optimal", trip_km=4000.0, **pick(k)))
    write_csv("flight_points.csv", out)
    return out


PICK = ["pod_t", "seats", "D_bore", "Mach", "p_Pa", "WS_Pa", "span_m", "q_Pa", "CL", "LD", "alpha_deg",
        "drag_kN", "stroke_req_y_m", "stroke_av_y_m", "T_le_K", "T_skin_K", "dT_pass_K", "N_heat_h", "N_max_h",
        "cap_pair_Mpax", "n_pairs", "E_Wh_paxkm", "E_drag", "E_kinetic", "E_pump", "P_pod_MW", "P_station_GW",
        "lsm_Hz", "KE_tTNT", "breach_decel_g", "breach_lift_overload", "x_reflect_m", "beta", "crossover_km",
        "d2d_NYC_LAX_h", "status", "status_nt", "binding_nt"]


def pick(r):
    return {k: r.get(k) for k in PICK}


def stage_econ(rows, F=BASE):
    out = []
    cand = [r for r in rows if r["status_nt"] in ("pass", "marginal")]
    for net in NETS:
        for D in D_LIST:
            for M in sorted(set(r["Mach"] for r in cand)):
                if abs(M * 2 - round(M * 2)) > 1e-9:
                    continue                                   # every 0.5 Mach
                c = [r for r in cand if r["D_bore"] == D and r["Mach"] == M]
                if not c:
                    continue
                best, bh = None, None
                for r in c:
                    h = bc_hook(r, net, F)
                    key = h["bcn_sub_100_base"]
                    if best is None or key > best:
                        best, bh, br = key, h, r
                row = dict(Mach=M, D_bore=D, pod_t=br["pod_t"], seats=br["seats"], p_Pa=br["p_Pa"],
                           WS_Pa=br["WS_Pa"], status=br["status_nt"], E_Wh_paxkm=br["E_Wh_paxkm"],
                           cap_pair_Mpax=br["cap_pair_Mpax"], d2d_NYC_LAX_h=br["d2d_NYC_LAX_h"],
                           reach_4000km=int((M * C_REF) ** 2 <= 0.3 * G * 4000e3))
                row.update(bh)
                out.append(row)
    write_csv("flight_econ.csv", out)
    return out


LEVERS = [
    ("base (0.1 GN line, AR 2, moving block 10 s, 0.3 g)", {}),
    ("untensioned line (0.01 GN)", dict(T_line=0.01e9)),
    ("line tension 0.3 GN", dict(T_line=0.3e9)),
    ("line tension 1 GN", dict(T_line=1.0e9)),
    ("bridle residual 3 % (vs 10 %)", dict(eps=0.03)),
    ("low-AR wing (AR 1: span / 1.4)", dict(AR=1.0)),
    ("tip clearance 0.5 m (vs 1 m)", dict(clear=0.5)),
    ("gas stirring x3", dict(f_mix=3.0)),
    ("UHTC leading edge 2,200 K", dict(T_lim_le=2200.0)),
    ("brick-wall headway (0.5 g)", dict(rule="brick-wall")),
    ("2,500 km trip", dict(L_trip=2500e3)),
    ("comfort accel 0.15 g", dict(a_acc=0.15 * G)),
    ("strip impingement x10", dict(imp=10.0)),
    ("leaks x0.1 (welded liner)", dict(leak=0.1)),
    ("launch power cap 100 MW per pod (vs 200)", dict(P_pod_ok=100e6, P_pod_marg=400e6)),
    ("energy pass at 200 Wh/pax-km (vs 340)", dict(E_ok=200.0, E_marg=400.0)),
    ("30 % capture (24 M pax/yr)", dict(capture=0.3)),
]


def stage_levers():
    path = os.path.join(HERE, "flight_levers.csv")
    have = {r["lever"] for r in read_csv("flight_levers.csv")} if os.path.exists(path) and "--fresh" not in sys.argv else set()
    if not have and os.path.exists(path):
        os.remove(path)
    mg = np.round(np.arange(2.0, 12.001, 0.5), 2)
    pg = np.geomspace(5.0, 1000.0, 20)
    for name, mod in LEVERS:
        if name in have:
            continue
        F = dict(BASE, **mod)
        g = run_grid(F, pods=[17000.0, 47000.0], m_grid=mg, p_grid=pg)
        r = dict(lever=name)
        for D in D_LIST:
            for m in (17.0, 47.0):
                sub = [x for x in g if x["D_bore"] == D and x["pod_t"] == m]
                Mp = max_mach(sub)
                Mm = max_mach(sub, ok=("pass", "marginal"))
                r[f"D{D:g}_{m:g}t_Mmax_pass"] = Mp
                r[f"D{D:g}_{m:g}t_Mmax_marg"] = Mm
                r[f"D{D:g}_{m:g}t_bound"] = top_binding(sub, Mp)
                s6 = sorted([x for x in sub if x["Mach"] == 6.0 and x["status"] == "pass"], key=lambda x: x["p_Pa"])
                r[f"D{D:g}_{m:g}t_M6_p"] = f"{s6[0]['p_Pa']:.0f}-{s6[-1]['p_Pa']:.0f}" if s6 else "none"
                ps = [x for x in sub if x["status"] == "pass"]
                r[f"D{D:g}_{m:g}t_minE"] = min((x["E_Wh_paxkm"] for x in ps), default=float("nan"))
                r[f"D{D:g}_{m:g}t_maxcap"] = max((x["cap_pair_Mpax"] for x in ps), default=float("nan"))
        write_csv("flight_levers.csv", [r], mode="a")
        print(f"  lever {name}: " + ", ".join(f"D{D:g}/47t M{r[f'D{D:g}_47t_Mmax_pass']}" for D in D_LIST), flush=True)


# ----------------------------------------------------------------------------- figures
INK, INK2, GRIDC, SURF = "#0b0b0b", "#52514e", "#e4e3df", "#fcfcfb"
PASS, MARG, FAIL = "#2a78d6", "#a9c8ef", "#ecebe7"
S1, S2, S3, S4 = "#2a78d6", "#eb6834", "#1baf7a", "#eda100"


def _style(ax, title):
    ax.set_facecolor(SURF)
    ax.set_title(title, loc="left", fontsize=9, color=INK, fontweight="bold")
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    for s in ("left", "bottom"):
        ax.spines[s].set_color(INK2)
    ax.tick_params(colors=INK2, labelsize=7.5)
    ax.xaxis.label.set_color(INK2)
    ax.yaxis.label.set_color(INK2)
    ax.grid(color=GRIDC, lw=0.6)
    ax.set_axisbelow(True)


def _edges(v, log=False):
    v = np.asarray(v, float)
    if log:
        lv = np.log10(v)
        e = np.concatenate([[lv[0] - (lv[1] - lv[0]) / 2], (lv[1:] + lv[:-1]) / 2, [lv[-1] + (lv[-1] - lv[-2]) / 2]])
        return 10 ** e
    return np.concatenate([[v[0] - (v[1] - v[0]) / 2], (v[1:] + v[:-1]) / 2, [v[-1] + (v[-1] - v[-2]) / 2]])


def fig_maps():
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.colors import ListedColormap
    from matplotlib.lines import Line2D
    from matplotlib.patches import Patch
    rows = read_csv("flight_grid.csv")
    pods = sorted(set(r["pod_t"] for r in rows))
    fig, axs = plt.subplots(len(D_LIST), len(pods), figsize=(3.3 * len(pods), 2.75 * len(D_LIST) + 0.9),
                            sharex=True, sharey=True)
    fig.patch.set_facecolor(SURF)
    cmap = ListedColormap([FAIL, MARG, PASS])
    for i, D in enumerate(D_LIST):
        for j, m in enumerate(pods):
            ax = axs[i, j]
            sub = [r for r in rows if r["D_bore"] == D and r["pod_t"] == m]
            xs = np.array(sorted(set(r["Mach"] for r in sub)))
            ys = np.array(sorted(set(r["p_Pa"] for r in sub)))
            ix = {x: k for k, x in enumerate(xs)}
            iy = {y: k for k, y in enumerate(ys)}
            Z = np.zeros((len(ys), len(xs)))
            B = np.empty((len(ys), len(xs)), object)
            E = np.full((len(ys), len(xs)), np.nan)
            CP = np.full((len(ys), len(xs)), np.nan)
            for r in sub:
                a, b_ = iy[r["p_Pa"]], ix[r["Mach"]]
                Z[a, b_] = RANK[r["status_nt"]]
                B[a, b_] = r["binding_nt"]
                E[a, b_] = r.get("E_Wh_paxkm", np.nan) if r["alpha_deg"] < 89 else np.nan
                CP[a, b_] = r.get("cap_pair_Mpax", np.nan) if r["alpha_deg"] < 89 else np.nan
            _style(ax, f"{m:g} t / {int(sub[0]['seats'])} seats · bore {D:g} m")
            ax.pcolormesh(_edges(xs), _edges(ys, True), Z, cmap=cmap, vmin=-0.5, vmax=2.5, shading="flat")
            codes = {k: n for n, k in enumerate(sorted(set(B.ravel()) - {None}))}
            C = np.vectorize(lambda b: codes.get(b, -1))(B).astype(float)
            C[Z == 2] = np.nan
            for bnd, n_ in codes.items():
                mask = C == n_
                if mask.sum() < 0.04 * mask.size:
                    continue
                jj, ii = np.nonzero(mask)
                k = np.argmin((jj - np.median(jj)) ** 2 + (ii - np.median(ii)) ** 2)
                ax.text(xs[ii[k]], ys[jj[k]], SHORT.get(bnd, bnd), fontsize=6.3, color=INK, ha="center", va="center",
                        bbox=dict(boxstyle="round,pad=0.12", fc=SURF, ec="none", alpha=0.85))
            E[Z == 0] = np.nan
            with np.errstate(invalid="ignore"):
                if np.isfinite(E).any():
                    cs = ax.contour(xs, ys, E, levels=[150, 250], colors=S2, linewidths=0.8, linestyles=["-", "--"])
                    ax.clabel(cs, fmt=lambda v: f"{v:.0f} Wh", fontsize=5.5)
                if np.isfinite(CP).any() and np.nanmax(CP) > 81:
                    ax.contour(xs, ys, CP, levels=[81], colors=S3, linewidths=1.0)
            for L, ls in ((1000e3, ":"), (2500e3, "-.")):
                ax.axvline(math.sqrt(0.3 * G * L) / C_REF, color=INK2, lw=0.8, ls=ls)
            ax.set_yscale("log")
            if i == len(D_LIST) - 1:
                ax.set_xlabel("Mach (V / 295 m/s)", fontsize=8)
            if j == 0:
                ax.set_ylabel("tube pressure, Pa", fontsize=8)
    handles = [Patch(fc=PASS, label="passes every check"), Patch(fc=MARG, label="marginal (worst check within its band)"),
               Patch(fc=FAIL, label="fails (label = binding check)"),
               Line2D([], [], color=S2, lw=0.8, label="150 / 250 Wh per pax-km, where not failing (jet fuel ~340)"),
               Line2D([], [], color=S3, lw=1.0, label="one tube pair carries 81 M pax/yr (above-left)"),
               Line2D([], [], color=INK2, ls=":", label="trip cap 1,000 km at 0.3 g"),
               Line2D([], [], color=INK2, ls="-.", label="trip cap 2,500 km (4,000 km: M11.6)")]
    fig.legend(handles=handles, loc="lower center", ncol=4, fontsize=7.5, frameon=False)
    fig.suptitle("Internal flight: feasible Mach × tube pressure per bore and pod class (best wing loading per cell; "
                 "0.1 GN line, moving block 10 s, 4,000 km)", fontsize=10, color=INK, x=0.01, ha="left")
    fig.tight_layout(rect=(0, 0.06, 1, 0.97))
    fig.savefig(os.path.join(HERE, "flight_maps.png"), dpi=130, facecolor=SURF)
    plt.close(fig)
    print("wrote flight_maps.png")


def fig_summary():
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    st = read_csv("flight_stroke.csv")
    ec = read_csv("flight_econ.csv")
    band = read_csv("flight_band.csv")
    fig, axs = plt.subplots(1, 4, figsize=(16, 3.9))
    fig.patch.set_facecolor(SURF)
    cols = {20.0: S2, 40.0: S4, 80.0: S1}
    # (a) stroke needed vs Mach and tension, D 40
    ax = axs[0]
    _style(ax, "a. Steering stroke needed (lateral, 4σ), 40 m bore")
    for T, ls in ((0.01, ":"), (0.1, "-"), (0.3, "--"), (1.0, "-.")):
        s = [r for r in st if r["T_GN"] == T and r["D_bore"] == 40.0 and abs(r["eps"] - 0.1) < 1e-6]
        ax.plot([r["Mach"] for r in s], [r["stroke_y_m"] for r in s], color=S1, ls=ls, lw=1.6, label=f"{T:g} GN")
    for b, c in ((5.0, INK2), (20.0, S2)):
        ax.axhline((40 - b) / 2 - 1, color=c, lw=0.9, ls=":")
        ax.text(12, (40 - b) / 2 - 1, f" span {b:g} m", fontsize=6.5, color=c, va="bottom", ha="right")
    ax.set_xlabel("Mach")
    ax.set_ylabel("stroke, m")
    ax.legend(fontsize=7, frameon=False, title="line tension", title_fontsize=7)
    # (b) highest passing / marginal Mach per pod and bore (4,000 km)
    ax = axs[1]
    _style(ax, "b. Highest feasible Mach (filled: pass, open: marginal)")
    pods = sorted(set(r["pod_t"] for r in band))
    for k, D in enumerate(D_LIST):
        for j_, m in enumerate(pods):
            s = [r for r in band if r["D_bore"] == D and r["pod_t"] == m]
            mp = [r["Mach"] for r in s if r["p_pass_lo"] == r["p_pass_lo"]]
            mm = [r["Mach"] for r in s if r["p_marg_lo"] == r["p_marg_lo"]]
            x = j_ + (k - 1) * 0.22
            if mm:
                ax.plot([x], [max(mm)], "o", mfc="none", mec=cols[D], ms=7, mew=1.6)
            if mp:
                ax.plot([x], [max(mp)], "o", color=cols[D], ms=7)
            if not mm:
                ax.plot([x], [2.0], "x", color=cols[D], ms=6)
    for D in D_LIST:
        ax.plot([], [], "o", color=cols[D], label=f"bore {D:g} m")
    ax.set_xticks(range(len(pods)))
    ax.set_xticklabels([f"{m:g} t" for m in pods], fontsize=7.5)
    ax.set_ylim(1.5, 12.5)
    ax.set_ylabel("Mach")
    ax.text(0, 1.7, "x = none", fontsize=6.5, color=INK2)
    ax.legend(fontsize=7, frameon=False, loc="upper right")
    # (c) energy along the B/C-best point
    ax = axs[2]
    _style(ax, "c. Energy at the B/C-best point (4,000 km)")
    for D in D_LIST:
        s = sorted([r for r in ec if r["D_bore"] == D and r["net"] == "NYC-LAX direct"], key=lambda r: r["Mach"])
        if not s:
            continue
        ax.plot([r["Mach"] for r in s], [r["E_Wh_paxkm"] for r in s], color=cols[D], lw=1.6, marker="o", ms=3,
                label=f"bore {D:g} m")
    ax.axhline(340, color=INK2, ls="--", lw=0.9)
    ax.text(2.1, 345, "jet fuel ~340", fontsize=6.5, color=INK2, va="bottom")
    ax.set_xlabel("Mach")
    ax.set_ylabel("Wh per pax-km")
    ax.legend(fontsize=7, frameon=False)
    # (d) B/C net vs subsonic, 100 % capture, base and low cost
    ax = axs[3]
    _style(ax, "d. B/C net vs subsonic, 100 %, NYC–LAX direct")
    for D in D_LIST:
        s = sorted([r for r in ec if r["D_bore"] == D and r["net"] == "NYC-LAX direct"], key=lambda r: r["Mach"])
        if not s:
            continue
        ax.plot([r["Mach"] for r in s], [r["bcn_sub_100_base"] for r in s], color=cols[D], lw=1.6, label=f"{D:g} m base")
        ax.plot([r["Mach"] for r in s], [r["bcn_sub_100_low"] for r in s], color=cols[D], lw=1.2, ls="--",
                label=f"{D:g} m low cost")
    ax.axhline(1.0, color=INK, lw=0.8)
    ax.set_xlabel("Mach")
    ax.set_ylabel("B/C")
    ax.legend(fontsize=6.5, frameon=False, ncol=2)
    fig.tight_layout()
    fig.savefig(os.path.join(HERE, "flight_summary.png"), dpi=130, facecolor=SURF)
    plt.close(fig)
    print("wrote flight_summary.png")


# ----------------------------------------------------------------------------- main
def main():
    args = [a for a in sys.argv[1:] if not a.startswith("--")]
    stages = args or ["stroke", "grid", "band", "ceilings", "gasdyn", "points", "econ", "levers", "figs"]
    rows = None
    if "stroke" in stages:
        stage_stroke()
    if "grid" in stages:
        rows = stage_grid()
    if any(s in stages for s in ("band", "points", "econ")) and rows is None:
        rows = read_csv("flight_grid.csv")
    if "band" in stages:
        stage_band(rows)
    if "ceilings" in stages:
        stage_ceilings()
    if "gasdyn" in stages:
        stage_gasdyn()
    if "points" in stages:
        stage_points(rows)
    if "econ" in stages:
        stage_econ(rows)
    if "levers" in stages:
        stage_levers()
    if "figs" in stages:
        fig_maps()
        fig_summary()


if __name__ == "__main__":
    main()
