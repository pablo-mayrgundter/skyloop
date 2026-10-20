"""Validation of the propagation chain (source -> augmented Burgers -> ground).

Runs Concorde and an SR-71-like case through ``ground_signature`` and compares with
published ground booms and with two far-field estimates:
  * Carlson simplified method (NASA TP-1122, 1978) with K_p = K_t = 1 and the shape
    factor computed from our own F-function: K_S = 2^{1/4} g (g+1)^{-1/2} sqrt(I)/l^{3/4},
    I = max_y int_0^y F dy (this is Whitham's far-field rule written in Carlson's form);
  * the stratified weak-shock (equal-area) asymptote using the same ray tube:
    age X = int k(z) G(z)/G(z0) ds, P = sqrt(2 I_t / X) G(0)/G(z0), T = 2 sqrt(2 I_t X).

Outputs results/validation_propagation.png and results/validation_propagation.csv.
Usage:  python boom-shield/sim/validate_propagation.py
"""
from __future__ import annotations

import csv
import os
import sys
import time
from dataclasses import replace

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)

import matplotlib  # noqa: E402

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

from boomsim.atmosphere import std_atm  # noqa: E402
from boomsim.propagate import RayPath, ground_signature  # noqa: E402
from boomsim.source import CONCORDE, SR71, f_function, near_field  # noqa: E402

PSF = 47.880259  # Pa per lb/ft^2
FS = 50e3
OUT = os.path.join(HERE, "results")

# Published ground booms (NASA Dryden fact sheet FS-2003-11-016 "Sonic Booms"):
#   Concorde  M 2.0, 52 000 ft: 1.94 psf;   SR-71  M 3.0, 80 000 ft: 0.9 psf.
PUBLISHED = {"Concorde": 1.94, "SR-71": 0.9}


def metrics(t, p):
    dt = t[1] - t[0]
    pk = p.max()
    i_pk = int(np.argmax(p))
    seg = p[: i_pk + 1]
    i10 = np.flatnonzero(seg < 0.1 * pk)[-1]
    i90 = np.flatnonzero(seg < 0.9 * pk)[-1]
    i_mn = int(np.argmin(p))
    i_end = i_mn + np.flatnonzero(p[i_mn:] > 0.1 * p.min())[0]
    return dict(peak_pa=pk, peak_psf=pk / PSF, min_pa=p.min(),
                rise_ms=(i90 - i10) * dt * 1e3, duration_ms=(i_end - i10) * dt * 1e3,
                t_front=t[i10])


def carlson(ac, reflection=1.9):
    """Carlson TP-1122 level-flight estimate with K_p = K_t = 1 (SI units)."""
    y = np.linspace(0, 3 * ac.length, 30001)
    F = f_function(ac, y)
    I = np.max(np.cumsum(F) * (y[1] - y[0]))
    ks = 2 ** 0.25 * 1.4 / np.sqrt(2.4) * np.sqrt(I) / ac.length ** 0.75
    a_v, a_g = std_atm(ac.altitude), std_atm(0.0)
    h = ac.altitude
    dp = reflection * np.sqrt(float(a_v.p) * float(a_g.p)) * ac.beta ** 0.25 \
        * h ** -0.75 * ac.length ** 0.75 * ks
    dt = 3.42 / float(a_v.c) * ac.mach / ac.beta ** 0.75 * h ** 0.25 * ac.length ** 0.75 * ks
    return dict(K_S=ks, peak_pa=dp, duration_ms=dt * 1e3)


def weak_shock_asymptote(ac, r0, reflection=1.9, fs=FS):
    """Equal-area far-field N-wave through the stratified ray tube (lossless)."""
    t, p = near_field(ac, r0, fs)
    dt = t[1] - t[0]
    I_t = np.max(np.cumsum(p) * dt)                      # Pa s, first positive lobe
    ray = RayPath(ac.altitude, ac.V)
    z0 = ac.altitude - r0
    z = np.linspace(z0, 0.0, 4001)
    a = ray.atm(z)
    k = 1.2 / (a.rho * a.c ** 3)
    g = ray.amplitude_ratio(z0, z)                       # G(z)/G(z0)
    s = ray.s(z)
    X = np.trapezoid(k * g, s)
    P = np.sqrt(2 * I_t / X) * g[-1] * reflection
    return dict(peak_pa=P, duration_ms=2 * np.sqrt(2 * I_t * X) * 1e3, age=X)


def run(ac, label=None, **kw):
    st = {}
    t0 = time.time()
    t, p = ground_signature(ac, fs=kw.pop("fs", FS), stats=st, **kw)
    m = metrics(t, p)
    m.update(case=label or ac.name, runtime_s=time.time() - t0, **st)
    print(f"{m['case']:<32s} peak {m['peak_pa']:6.1f} Pa ({m['peak_psf']:.2f} psf) "
          f"min {m['min_pa']:6.1f}  rise {m['rise_ms']:.2f} ms  dur {m['duration_ms']:.0f} ms"
          f"  [{m['runtime_s']:.1f} s, {st['nonlinear_steps']} NL steps]")
    return t, p, m


def main():
    os.makedirs(OUT, exist_ok=True)
    rows, sigs = [], {}
    for ac in (CONCORDE, SR71):
        r0 = 1.5 * ac.length
        t, p, m = run(ac)
        sigs[ac.name] = (t, p)
        c = carlson(ac)
        w = weak_shock_asymptote(ac, r0)
        m.update(published_psf=PUBLISHED[ac.name], carlson_KS=c["K_S"],
                 carlson_peak_pa=c["peak_pa"], carlson_dur_ms=c["duration_ms"],
                 weakshock_peak_pa=w["peak_pa"], weakshock_dur_ms=w["duration_ms"])
        print(f"   Carlson (K_p=K_t=1, K_S={c['K_S']:.3f}): {c['peak_pa']:.1f} Pa "
              f"({c['peak_pa'] / PSF:.2f} psf), {c['duration_ms']:.0f} ms;  stratified "
              f"weak-shock asymptote: {w['peak_pa']:.1f} Pa, {w['duration_ms']:.0f} ms")
        rows.append(m)
        # sensitivities
        rows.append(run(ac, f"{ac.name} lossless", lossless=True)[2])
        rows.append(run(ac, f"{ac.name} r0=1.0L", r0=1.0 * ac.length)[2])
        rows.append(run(ac, f"{ac.name} r0=3.0L", r0=3.0 * ac.length)[2])
        rows.append(run(ac, f"{ac.name} no dispersion", dispersion=False)[2])
    rows.append(run(replace(CONCORDE, altitude=15850.0), "Concorde 52 kft (15.85 km)")[2])
    rows.append(run(CONCORDE, "Concorde fs=100 kHz", fs=100e3)[2])

    keys = ["case", "peak_pa", "peak_psf", "min_pa", "rise_ms", "duration_ms", "runtime_s",
            "nonlinear_steps", "linear_steps", "n_samples", "published_psf", "carlson_KS",
            "carlson_peak_pa", "carlson_dur_ms", "weakshock_peak_pa", "weakshock_dur_ms"]
    with open(os.path.join(OUT, "validation_propagation.csv"), "w", newline="") as fh:
        wr = csv.DictWriter(fh, fieldnames=keys, extrasaction="ignore")
        wr.writeheader()
        for r in rows:
            wr.writerow({k: (f"{r[k]:.4g}" if isinstance(r.get(k), (float, np.floating))
                             else r.get(k, "")) for k in keys})

    # ---------------------------------------------------------------- figure
    col = {"Concorde": "#2a78d6", "SR-71": "#eb6834"}
    ink, muted, grid = "#0b0b0b", "#52514e", "#e4e3df"
    plt.rcParams.update({"font.size": 9, "axes.edgecolor": muted, "axes.labelcolor": ink,
                         "xtick.color": muted, "ytick.color": muted,
                         "axes.spines.top": False, "axes.spines.right": False})
    fig, axs = plt.subplots(2, 2, figsize=(11, 7.2), facecolor="#fcfcfb")
    for ax in axs.flat:
        ax.set_facecolor("#fcfcfb")
        ax.grid(True, color=grid, lw=0.6)
        ax.axhline(0, color=muted, lw=0.6)
    ax = axs[0, 0]
    for ac in (CONCORDE, SR71):
        t, p = near_field(ac, 1.5 * ac.length, FS)
        ax.plot(t * 1e3, p, color=col[ac.name], lw=1.6,
                label=f"{ac.name}  (r0 = 1.5 L = {1.5 * ac.length:.0f} m)")
    ax.set_xlim(-10, 140)
    ax.set_xlabel("time  [ms]")
    ax.set_ylabel("overpressure  [Pa]")
    ax.set_title("Near field (Whitham, linear)\nspikes = Sears-Haack tail-closure log singularity",
                 loc="left", color=ink)
    ax.legend(frameon=False, loc="upper right")

    for ax, name in ((axs[0, 1], "Concorde"), (axs[1, 0], "SR-71")):
        t, p = sigs[name]
        m = next(r for r in rows if r["case"] == name)
        tt = (t - m["t_front"]) * 1e3
        ax.plot(tt, p, color=col[name], lw=1.6, label="this model (ground, x1.9)")
        ax.axhline(PUBLISHED[name] * PSF, color=muted, lw=1.0, ls="--")
        ax.text(0.99, PUBLISHED[name] * PSF, f"published {PUBLISHED[name]} psf",
                transform=ax.get_yaxis_transform(), ha="right", va="bottom", color=muted)
        ac = CONCORDE if name == "Concorde" else SR71
        ax.set_xlim(-40, m["duration_ms"] + 60)
        ax.set_xlabel("time from bow shock  [ms]")
        ax.set_ylabel("overpressure  [Pa]")
        ax.set_title(f"{name} at ground: M {ac.mach:.1f}, {ac.altitude / 1e3:.1f} km\n"
                     f"{m['peak_pa']:.0f} Pa ({m['peak_psf']:.2f} psf), "
                     f"{m['duration_ms']:.0f} ms", loc="left", color=ink)

    ax = axs[1, 1]
    for name in ("Concorde", "SR-71"):
        t, p = sigs[name]
        m = next(r for r in rows if r["case"] == name)
        ax.plot((t - m["t_front"]) * 1e3, p / m["peak_pa"], color=col[name], lw=1.6,
                label=f"{name}: 10-90 % rise {m['rise_ms']:.2f} ms")
    ax.set_xlim(-0.5, 3.0)
    ax.set_ylim(-0.05, 1.08)
    ax.set_xlabel("time from 10 % point  [ms]")
    ax.set_ylabel("p / peak")
    ax.set_title("Bow-shock structure at ground\n(O2/N2 relaxation-limited, no turbulence)",
                 loc="left", color=ink)
    ax.legend(frameon=False, loc="lower right")
    fig.tight_layout()
    fig.savefig(os.path.join(OUT, "validation_propagation.png"), dpi=150,
                facecolor=fig.get_facecolor())
    print("wrote", os.path.join(OUT, "validation_propagation.png"))


if __name__ == "__main__":
    main()
