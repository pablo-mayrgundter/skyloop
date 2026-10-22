"""Exploratory study: plane-wave transmission of a boom through candidate shield stacks
at 20 km, for several aircraft Mach numbers. Not a full simulation (no propagation).

Question: a helium/hydrogen layer is "subsonic" relative to the boom trace speed V when
V < c_gas, so the boom only tunnels through it evanescently (|T| ~ exp(-kappa d),
kappa = omega sqrt(1/V^2 - 1/c_gas^2)). How much does that beat the plain mass law of
the same lift-limited areal mass at boom frequencies (~1-20 Hz)?

Every stack closes its own lift budget (slab model, lift_margin = 0 - optimistic).

Run:  python boom-shield/sim/explore_shield.py
"""
import csv
import os
import sys
from types import SimpleNamespace

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from boomsim.shield import GasLayer, Membrane, Stack, apply, gas, lift_margin, transfer  # noqa: E402

RESULTS = os.path.join(HERE, "results")

T20, P20 = 216.65, 5474.9
_air = gas("air", T20, P20)
ATM = SimpleNamespace(T=T20, p=P20, rho=_air.rho, c=_air.c)

MACHS = [1.4, 2.0, 2.5, 3.2]
DS = [15.0, 30.0, 60.0, 120.0]
FILLS = {
    "helium": gas("helium", T20, P20),
    "hydrogen": gas("hydrogen", T20, P20),
    "air@0.1p": gas("air", T20, 0.1 * P20),
}
F = np.logspace(np.log10(0.5), 3, 600)
PROBE_F = [2.0, 5.0, 10.0, 20.0]
DURATIONS = [0.1, 0.3]
FS = 20000.0


def budget(g, d):
    return (ATM.rho - g.rho) * d


def cases():
    """(case, fill, d, stack, total membrane mass)."""
    m1 = budget(FILLS["helium"], 15.0)
    yield "single membrane (15 m He budget)", "-", 0.0, Stack([Membrane(m1)]), m1
    for name, g in FILLS.items():
        for d in DS:
            b = budget(g, d)
            st = Stack([Membrane(b / 2), GasLayer(g, d), Membrane(b / 2)])
            assert abs(lift_margin(st, ATM)) < 1e-12
            yield f"membrane/{name}/membrane", name, d, st, b


def nwave(dur, pad=3.0):
    t = np.arange(-pad, dur + pad, 1 / FS)
    p = np.where((t >= 0) & (t < dur), 1.0 - 2.0 * t / dur, 0.0)
    return t, p


def pos_impulse(t, p):
    """Integral over the contiguous positive lobe containing the peak."""
    i = int(np.argmax(p))
    if p[i] <= 0:
        return 0.0
    lo = i
    while lo > 0 and p[lo - 1] > 0:
        lo -= 1
    hi = i
    while hi < p.size - 1 and p[hi + 1] > 0:
        hi += 1
    return float(np.trapezoid(p[lo:hi + 1], t[lo:hi + 1]))


def main():
    os.makedirs(RESULTS, exist_ok=True)
    waves = {dur: nwave(dur) for dur in DURATIONS}
    ref_imp = {dur: pos_impulse(*waves[dur]) for dur in DURATIONS}
    rows, curves = [], {}
    for mach in MACHS:
        V = mach * ATM.c
        for case, fill, d, st, mass in cases():
            T, _ = transfer(st, F, V, ATM)
            curves[(mach, case, d)] = np.abs(T)
            ml = Stack([Membrane(mass)])  # plain mass law, same total areal mass
            row = dict(mach=mach, V=round(V, 1), case=case, fill=fill, d_m=d,
                       membrane_mass_kg_m2=round(mass, 4),
                       regime=("-" if fill == "-" else
                               "evanescent" if V < FILLS[fill].c else "propagating"))
            Tp, _ = transfer(st, PROBE_F, V, ATM)
            Tm, _ = transfer(ml, PROBE_F, V, ATM)
            for fp, a, b in zip(PROBE_F, np.abs(Tp), np.abs(Tm)):
                row[f"T_{fp:g}Hz_dB"] = round(20 * np.log10(max(a, 1e-300)), 2)
                row[f"masslaw_{fp:g}Hz_dB"] = round(20 * np.log10(b), 2)
            for dur in DURATIONS:
                t, p = waves[dur]
                out = apply(st, t, p, V, ATM)
                outm = apply(ml, t, p, V, ATM)
                row[f"peak_{dur:g}s"] = round(float(out.max()), 4)
                row[f"imp_{dur:g}s"] = round(pos_impulse(t, out) / ref_imp[dur], 4)
                row[f"masslaw_peak_{dur:g}s"] = round(float(outm.max()), 4)
                row[f"masslaw_imp_{dur:g}s"] = round(pos_impulse(t, outm) / ref_imp[dur], 4)
            rows.append(row)

    csv_path = os.path.join(RESULTS, "shield_transfer_summary.csv")
    with open(csv_path, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
        w.writeheader()
        w.writerows(rows)

    # ---- console table
    print(f"ambient 20 km: rho={ATM.rho:.4f} c={ATM.c:.1f}; c_He={FILLS['helium'].c:.0f} "
          f"c_H2={FILLS['hydrogen'].c:.0f} m/s")
    hdr = (f"{'M':>4} {'case':<34}{'d':>5} {'m':>6} {'regime':>11} {'|T|10Hz dB':>10} {'ML':>7}"
           f" {'pk.1':>6} {'ML':>5} {'I.1':>5} {'pk.3':>6} {'ML':>5} {'I.3':>5}")
    print(hdr)
    for r in rows:
        print(f"{r['mach']:>4} {r['case']:<34}{r['d_m']:>5.0f} {r['membrane_mass_kg_m2']:>6.2f}"
              f" {r['regime']:>11} {r['T_10Hz_dB']:>10.1f} {r['masslaw_10Hz_dB']:>7.1f}"
              f" {r['peak_0.1s']:>6.3f} {r['masslaw_peak_0.1s']:>5.2f} {r['imp_0.1s']:>5.2f}"
              f" {r['peak_0.3s']:>6.3f} {r['masslaw_peak_0.3s']:>5.2f} {r['imp_0.3s']:>5.2f}")
    print(f"wrote {csv_path}")

    # ---- figure: rows = Mach, cols = fill gas; d as an ordinal blue ramp
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    ramp = ["#86b6ef", "#3987e5", "#1c5cab", "#0d366b"]
    ink, ink2, grid = "#0b0b0b", "#52514e", "#e4e3df"
    plt.rcParams.update({"font.size": 9, "axes.edgecolor": ink2, "axes.labelcolor": ink,
                         "xtick.color": ink2, "ytick.color": ink2})
    fig, axes = plt.subplots(len(MACHS), len(FILLS), figsize=(13, 12), sharex=True, sharey=True)
    fig.patch.set_facecolor("#fcfcfb")
    single = "single membrane (15 m He budget)"
    for i, mach in enumerate(MACHS):
        V = mach * ATM.c
        for j, (name, g) in enumerate(FILLS.items()):
            ax = axes[i, j]
            ax.set_facecolor("#fcfcfb")
            ax.axvspan(1, 20, color=grid, alpha=0.5, lw=0)
            for k, d in enumerate(DS):
                y = 20 * np.log10(np.maximum(curves[(mach, f"membrane/{name}/membrane", d)], 1e-30))
                ax.semilogx(F, y, color=ramp[k], lw=2, label=f"{d:g} m")
            ax.semilogx(F, 20 * np.log10(curves[(mach, single, 0.0)]), color=ink2, lw=1.5,
                        ls="--", label="single membrane\n1.14 kg/m²")
            regime = "evanescent" if V < g.c else "propagating"
            ax.set_title(f"M {mach}  (V={V:.0f} m/s) · {name} · {regime}", fontsize=9, color=ink)
            ax.set_ylim(-80, 3)
            ax.grid(True, which="major", color=grid, lw=0.6)
            for s in ("top", "right"):
                ax.spines[s].set_visible(False)
            if i == len(MACHS) - 1:
                ax.set_xlabel("frequency [Hz]")
            if j == 0:
                ax.set_ylabel("|T| [dB]")
    h, lab = axes[0, 0].get_legend_handles_labels()
    lab = [f"d = {x}" if x.endswith(" m") else x.replace("\n", " ") for x in lab]
    fig.legend(h, lab, loc="lower center", ncol=5, frameon=False, fontsize=9)
    fig.suptitle("Shield pressure transmission at 20 km, lift-closed membrane/gas/membrane stacks "
                 "(shaded: 1–20 Hz boom band; lossless, so propagating layers show sharp resonances)",
                 color=ink)
    fig.tight_layout(rect=(0, 0.03, 1, 0.97))
    png = os.path.join(RESULTS, "shield_transfer.png")
    fig.savefig(png, dpi=130, facecolor=fig.get_facecolor())
    print(f"wrote {png}")


if __name__ == "__main__":
    main()
