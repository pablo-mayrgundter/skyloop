"""Run matrix: raft-of-tubes vs slab (and soft-tie slab) cross-plane BEM, Concorde M2, h = 250 m,
W = 1 km, shield top at 20 km. Output: results/raft_bem.csv (per config, f, phi) and
results/raft_bem_summary.csv. One process; set OMP_NUM_THREADS to limit BLAS threads.

Metric: |p_total / p_incident| on the matching circle R = 3 km (the 3-D study's R_m),
phi = 0..50 deg (symmetric geometry), and on a near-field line 100 m below the shield
(|y| <= 300 m). Linear incident field (exact slender-body), so ratios are independent of
the aging fit."""
import csv
import os
import sys
import time

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import raft_bem as RB  # noqa: E402

OUT = os.path.join(HERE, "results")
os.makedirs(OUT, exist_ok=True)
FREQS = [2.0, 3.0, 5.0, 7.0, 10.0, 15.0, 20.0, 30.0]
PHIS = np.arange(0.0, 50.1, 1.0)
R_OBS = 3000.0
H = 250.0
W = 1000.0
AC = RB.aircraft(2.0, H)
LIFT = RB.RHO - RB.H2.rho   # kg/m^3


class Mount:
    """Top face with areal mass m and a foundation spring Kf [Pa/m] to an external fixed
    reference (hold-down net), not to the bottom face: Z = j w m + Kf/(j w)."""

    def __init__(self, m, Kf):
        self.m, self.Kf = float(m), float(Kf)

    def impedance(self, omega, sx):
        omega = np.asarray(omega, float)
        return 1j * omega * self.m + self.Kf / (1j * omega)

CONFIGS = [
    # key, kind, params
    ("slab_H2_30", "slab", dict(T=30.0)),
    ("slab_H2_60", "slab", dict(T=60.0)),
    ("slab_H2_20", "slab", dict(T=20.0)),
    ("raft_D30_limp", "raft", dict(D=30.0, b=1.0, dp=0.0, Et=0.0)),
    ("raft_D30_dp100", "raft", dict(D=30.0, b=1.0, dp=100.0)),
    ("raft_D30_dp300", "raft", dict(D=30.0, b=1.0, dp=300.0)),
    ("raft_D30_dp1000", "raft", dict(D=30.0, b=1.0, dp=1000.0)),
    ("raft_D30_dp3000", "raft", dict(D=30.0, b=1.0, dp=3000.0)),
    ("raft_D30_dp300_b4", "raft", dict(D=30.0, b=4.0, dp=300.0)),
    ("raft_D30_dp300_gap2", "raft", dict(D=30.0, gap=2.0, dp=300.0)),
    ("raft_D30_dp300_gap5", "raft", dict(D=30.0, gap=5.0, dp=300.0)),
    ("raft_D30_dp300_gap15", "raft", dict(D=30.0, gap=15.0, dp=300.0)),
    ("raft_D30_dp300_Et1e6", "raft", dict(D=30.0, b=1.0, dp=300.0, Et=1e6)),
    ("raft_D30_dp300_light", "raft", dict(D=30.0, b=1.0, dp=300.0, m_s=0.2)),
    ("raft_D20_dp300", "raft", dict(D=20.0, b=1.0, dp=300.0)),
    ("raft_D60_dp300", "raft", dict(D=60.0, b=1.0, dp=300.0)),
    ("raft_D60_dp1000", "raft", dict(D=60.0, b=1.0, dp=1000.0)),
    # (d) soft-tie zero-pressure slab: heavy top / light bottom, smeared ties
    ("slab_H2_30_heavytop_K0", "slab", dict(T=30.0, m_top=LIFT * 30 - 0.1, m_bot=0.1, tie_K=0.0)),
    ("slab_H2_30_heavytop_K3", "slab", dict(T=30.0, m_top=LIFT * 30 - 0.1, m_bot=0.1, tie_K=3.0)),
    ("slab_H2_30_heavytop_K10", "slab", dict(T=30.0, m_top=LIFT * 30 - 0.1, m_bot=0.1, tie_K=10.0)),
    ("slab_H2_30_heavytop_K30", "slab", dict(T=30.0, m_top=LIFT * 30 - 0.1, m_bot=0.1, tie_K=30.0)),
    ("slab_H2_30_heavytop_K100", "slab", dict(T=30.0, m_top=LIFT * 30 - 0.1, m_bot=0.1, tie_K=100.0)),
    ("slab_H2_30_K10", "slab", dict(T=30.0, tie_K=10.0)),
    ("slab_H2_30_K30", "slab", dict(T=30.0, tie_K=30.0)),
    ("slab_H2_30_K100", "slab", dict(T=30.0, tie_K=100.0)),
    ("slab_H2_30_K300", "slab", dict(T=30.0, tie_K=300.0)),
    # second pass: isolate the hoop-stretch effect, realistic (non lift-closed) masses
    ("raft_D30_dp0_Et1e5", "raft", dict(D=30.0, b=1.0, dp=0.0, Et=1e5)),
    ("raft_D30_dp300_Et1e4", "raft", dict(D=30.0, b=1.0, dp=300.0, Et=1e4)),
    ("raft_D30_dp1000_light", "raft", dict(D=30.0, b=1.0, dp=1000.0, m_s=0.2)),
    ("slab_H2_30_light_K10", "slab", dict(T=30.0, m_top=0.15, m_bot=0.05, tie_K=10.0)),
    ("slab_H2_30_light_hold10", "slab", dict(T=30.0, m_top=Mount(0.15, 10.0), m_bot=0.05, tie_K=3.0)),
    ("slab_H2_30_light_hold100", "slab", dict(T=30.0, m_top=Mount(0.15, 100.0), m_bot=0.05, tie_K=3.0)),
    ("slab_H2_20_heavytop_K10", "slab", dict(T=20.0, m_top=LIFT * 20 - 0.1, m_bot=0.1, tie_K=10.0)),
    # breathing-stiffness sweep (Et/R^2 = 111 ... 1778 Pa/m)
    ("raft_D30_dp300_Et3e5", "raft", dict(D=30.0, b=1.0, dp=300.0, Et=3e5)),
    ("raft_D20_dp300_Et4e4", "raft", dict(D=20.0, b=1.0, dp=300.0, Et=4e4)),
    ("raft_D60_dp300_Et4e5", "raft", dict(D=60.0, b=1.0, dp=300.0, Et=4e5)),
    ("raft_D60_dp300_Et1.6e6", "raft", dict(D=60.0, b=1.0, dp=300.0, Et=1.6e6)),
]


def points(depth_bottom):
    ph = np.radians(PHIS)
    circ = R_OBS * np.stack([np.sin(ph), -np.cos(ph)], 1)
    yl = np.arange(-300.0, 300.1, 10.0)
    line = np.stack([yl, np.full(yl.size, -depth_bottom - 100.0)], 1)
    return circ, line


def main(only=None):
    rows_path = os.path.join(OUT, "raft_bem.csv")
    done = set()
    if os.path.exists(rows_path):
        with open(rows_path) as fh:
            for r in csv.DictReader(fh):
                done.add((r["config"], float(r["f"])))
    new = not os.path.exists(rows_path)
    fh = open(rows_path, "a", newline="")
    wr = csv.writer(fh)
    if new:
        wr.writerow(["config", "f", "where", "phi_or_y", "dB_re_incident"])
    for key, kind, par in CONFIGS:
        if only and key not in only:
            continue
        for f in FREQS:
            if (key, f) in done:
                continue
            t = time.time()
            pp = dict(par)
            if kind == "slab":
                T = pp.pop("T")
                circ, line = points(H + T)
                X = np.vstack([circ, line])
                p, p0, _ = RB.solve_slab(AC, f, X, W, H, T, **pp)
                npan = -1
            else:
                D = pp.pop("D")
                circ, line = points(H + D)
                X = np.vstack([circ, line])
                p, p0, geo = RB.solve_raft(AC, f, X, D, W, H, **pp)
                npan = len(geo["P"])
            dB = 20 * np.log10(np.abs(p / p0))
            nc = circ.shape[0]
            for i in range(nc):
                wr.writerow([key, f, "circle", PHIS[i], round(float(dB[i]), 3)])
            for i in range(line.shape[0]):
                wr.writerow([key, f, "line", line[i, 0], round(float(dB[nc + i]), 3)])
            fh.flush()
            print(f"{key:30s} f={f:5.1f} panels={npan:5d} {time.time() - t:6.1f}s  "
                  f"undertrack {dB[0]:6.1f}  circle-mean {10*np.log10(np.mean(10**(dB[:nc]/10))):6.1f}", flush=True)
    fh.close()


if __name__ == "__main__":
    main(sys.argv[1:] or None)
