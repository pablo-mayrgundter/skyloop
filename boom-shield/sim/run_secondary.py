"""Secondary sonic booms with and without BoomShield (boomsim.secondary).

Cases: Concorde-class (M 2.0) and SR-71-class (M 3.2) at 20.25 km, shield top at 20 km
(H2 30 m, 1.24 kg/m^2 faces, 1 km wide: every ray with |phi| < atan(500/250) = 63.4 deg
hits it).  Atmospheres: USSA76 (calm) and representative middle-atmosphere winds.

Families (NPRM, proposed 14 CFR 1.1: "direct" = over-the-top, refracted back by the upper
atmosphere; "indirect" = first reflected by the surface, then refracted back):
  unshielded: direct (all rays that leave upward, incl. those beyond the lateral cut-off),
              indirect (primary boom reflected by the ground, x0.9)
  shielded:   direct (|phi| > 63.4 deg only), shield (downward field reflected by the shield
              top with R(f)), transmit (field transmitted with T(f), reflected by the ground
              or turning below 20 km)

Outputs: results/secondary_booms.csv (record = 'bin' rows: binned ground levels; 'ray'
rows: per-march diagnostics), results/secondary_booms.png; numbers for secondary_booms.md
are printed.  Single process (about 1-2 h).  Usage:
    python run_secondary.py [--quick] [--cache DIR]
"""
from __future__ import annotations

import argparse
import csv
import dataclasses
import math
import os
import pickle
import sys
import tempfile
import time

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from boomsim import secondary as S                      # noqa: E402
from boomsim import shield as SH                        # noqa: E402
from boomsim.atmosphere import std_atm                  # noqa: E402
from boomsim.loudness import metrics                    # noqa: E402
from boomsim.source import CONCORDE, SR71               # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
RES = os.path.join(HERE, "results")
Z_A = 20250.0
Z_S = 20000.0
W_SHIELD = 1000.0
PHI_S = math.atan2(0.5 * W_SHIELD, Z_A - Z_S)
FACE = 1.24
LIMIT_PSF = 0.11
DPHI_FINE = 0.1          # deg (fan for landing geometry / binning)
DPHI_MARCH = 2.5         # deg (marched rays; plus branch end points)
APEX_MAX = 140e3         # rays turning higher are not marched (absorbed; see md)
DY = 3000.0
FS = 10e3

AIRCRAFT = {
    "Concorde M2.0": dataclasses.replace(CONCORDE, mach=2.0, altitude=Z_A),
    "SR-71 M3.2": dataclasses.replace(SR71, mach=3.2, altitude=Z_A),
}
ATMOS = {
    "calm": lambda: None,
    "winter tail": lambda: S.upper_wind("winter", 0.0),
    "winter head": lambda: S.upper_wind("winter", 180.0),
    "winter cross": lambda: S.upper_wind("winter", 90.0),
    "summer tail": lambda: S.upper_wind("summer", 180.0),
}
RUNS = [("Concorde M2.0", a) for a in ATMOS] + [("SR-71 M3.2", a) for a in
                                                 ("calm", "winter tail", "winter cross")]
FAMILY_LABEL = {"direct": "direct", "indirect": "indirect (ground-reflected primary)",
                "shield": "shield-reflected", "transmit": "shield-transmitted (indirect)"}


def stack():
    a = std_atm(Z_S)
    g = SH.gas("hydrogen", float(a.T), float(a.p))
    return SH.Stack([SH.Membrane(FACE), SH.GasLayer(g, 30.0), SH.Membrane(FACE)])


def march_phis(rows):
    """Marched azimuths for a family: a 2.5 deg grid inside each branch plus its ends."""
    rows = [r for r in rows if r["apex"] <= APEX_MAX]
    if not rows:
        return []
    phi = np.array([r["phi"] for r in rows])
    lab = S.branches(phi, np.array([r["Y"] for r in rows]), np.array([r["apex"] for r in rows]))
    out = set()
    d = math.radians(DPHI_FINE)
    for b in np.unique(lab):
        pb = np.sort(phi[lab == b])
        lo, hi = pb[0], pb[-1]
        out.add(round(math.degrees(lo + (0.5 * d if hi - lo > d else 0.0)), 3))
        out.add(round(math.degrees(hi - (0.5 * d if hi - lo > d else 0.0)), 3))
        g = np.arange(math.ceil(math.degrees(lo) / DPHI_MARCH) * DPHI_MARCH,
                      math.degrees(hi), DPHI_MARCH)
        for q in g:
            out.add(round(float(q), 3))
    return sorted(math.radians(q) for q in out)


def run_case(acn, atmn, cache, quick=False):
    ac = AIRCRAFT[acn]
    wind = ATMOS[atmn]()
    med = S.Medium(wind)
    two_sided = wind is not None and np.any(np.abs(med.v) > 0)
    step = DPHI_FINE * (4 if quick else 1)
    ph = np.arange(0.0, 180.0 + 1e-9, step)
    if two_sided:
        ph = np.concatenate([-ph[:0:-1], ph])
    phis = np.radians(ph)
    t0 = time.time()
    trU, rowsU = S.fan(ac, med, phis)
    sh = S.Shield(z=Z_S, phi_max=PHI_S, mode="reflect")
    inS = np.abs(phis) < PHI_S
    trS, rowsS = S.fan(ac, med, phis[inS], shield=sh)
    phi_c = [S.lateral_cutoff_phi(ac, med, side=+1)]
    if two_sided:
        phi_c.append(-S.lateral_cutoff_phi(ac, med, side=-1))
    print(f"  fans {time.time() - t0:.0f} s; lateral cut-off {np.degrees(phi_c)} deg", flush=True)
    fams = {
        "U:direct": ("direct", [r for r in rowsU if r["family"] == "direct"]),
        "U:indirect": ("indirect", [r for r in rowsU if r["family"] == "indirect"]),
        "S:direct": ("direct", [r for r in rowsU if r["family"] == "direct"
                                and abs(r["phi"]) >= PHI_S]),
        "S:shield": ("shield", [r for r in rowsS if r["family"] == "shield"]),
        "S:transmit": ("transmit", [r for r in rowsU if r["family"] in ("direct", "indirect")
                                    and abs(r["phi"]) < PHI_S]),
    }
    st = stack()
    marched = {}
    ray_rows = []
    for key in ("U:direct", "U:indirect", "S:shield", "S:transmit"):
        fam, rows = fams[key]
        mlist = march_phis(rows)
        if quick:
            mlist = mlist[::3]
        res_f = {}
        for q in mlist:
            ck = os.path.join(cache, f"{acn}|{atmn}|{fam}|{q:.6f}.pkl".replace(" ", "_"))
            if os.path.exists(ck):
                with open(ck, "rb") as fh:
                    out = pickle.load(fh)
            else:
                t1 = time.time()
                try:
                    path, res = S.march_family(ac, q, med, fam, shield=S.Shield(Z_S, PHI_S),
                                               stack=st, fs=FS)
                except Exception as e:          # pragma: no cover - logged, skipped
                    print(f"    march {fam} {math.degrees(q):.2f} failed: {e}", flush=True)
                    path, res = None, None
                out = None
                if res is not None:
                    out = dict(res=res, landing=path.landing, caustics=len(path.caustics),
                               dt=time.time() - t1)
                with open(ck, "wb") as fh:
                    pickle.dump(out, fh)
            if out is None:
                continue
            res = out["res"]
            res_f[q] = res
            ray_rows.append(dict(record="ray", aircraft=acn, atmosphere=atmn,
                                 config="unshielded" if key.startswith("U") else "shielded",
                                 family=fam, phi_deg=round(math.degrees(q), 3),
                                 Y_km=round(out["landing"]["y"] / 1e3, 2),
                                 apex_km=round(res.apex / 1e3, 2), z_ref_km=round(res.z_ref / 1e3, 2),
                                 E_h=res.E_h, max_dp_over_p=round(res.max_dp_over_p, 4),
                                 n_caustic_up=res.n_caustic_up, n_caustic_down=res.n_caustic_down,
                                 delay_min=round((out["landing"]["t"] - out["landing"]["x"]
                                                  / ac.V) / 60.0, 2)))
        marched[fam] = res_f
        print(f"  {key}: {len(rows)} fan landings, {len(res_f)} marched "
              f"({time.time() - t0:.0f} s)", flush=True)
    marched_S_direct = {q: r for q, r in marched["direct"].items() if abs(q) >= PHI_S}
    bins = {}
    for key, (fam, rows) in fams.items():
        mm = marched_S_direct if key == "S:direct" else marched.get(fam, {})
        rows_b = [r for r in rows if r["apex"] <= APEX_MAX + 5e3]
        bins[key] = {dy: S.bin_ground(rows_b, mm, med, math.radians(step), dY=dy)
                     for dy in (1500.0, DY, 6000.0)}
    return dict(bins=bins, ray_rows=ray_rows, phi_c=phi_c, trU=trU, rowsU=rowsU, rowsS=rowsS,
                med=med)


def bin_metrics(b):
    m = metrics(b["t"], b["p"])
    return dict(peak_pa=float(np.max(b["p"])), peak_abs_pa=float(np.max(np.abs(b["p"]))),
                pldb=m["pldb"], csel=m["csel"], asel=m["asel"])


def summarise(acn, atmn, out, csv_rows):
    summ = {}
    for key, by_dy in out["bins"].items():
        cfg = "unshielded" if key.startswith("U") else "shielded"
        fam = key.split(":")[1]
        best = None
        for dy, bl in by_dy.items():
            # coherent upper bound per Y: sum of branch peaks in the same bin
            byY = {}
            for b in bl:
                byY.setdefault(round(b["Y"]), []).append(b)
            for Yk, lst in byY.items():
                mets = [bin_metrics(b) for b in lst]
                k = int(np.argmax([m["peak_abs_pa"] for m in mets]))
                b, m = lst[k], mets[k]
                coh = float(sum(mm["peak_abs_pa"] for mm in mets))
                e_tot = float(sum(bb["energy"] for bb in lst))
                if dy == DY:
                    csv_rows.append(dict(record="bin", aircraft=acn, atmosphere=atmn, config=cfg,
                                         family=fam, Y_km=round(b["Y"] / 1e3, 2),
                                         peak_pa=round(m["peak_pa"], 5),
                                         peak_abs_pa=round(m["peak_abs_pa"], 5),
                                         peak_psf=round(m["peak_abs_pa"] / S.PSF, 6),
                                         coherent_peak_psf=round(coh / S.PSF, 6),
                                         pldb=round(m["pldb"], 2) if np.isfinite(m["pldb"]) else "",
                                         csel=round(m["csel"], 2), asel=round(m["asel"], 2),
                                         energy_pa2s=e_tot, n_branches=len(lst),
                                         apex_km=round(b["apex"] / 1e3, 1),
                                         delay_min=round(b["delay_s"] / 60.0, 2),
                                         phi_deg=round(b["phi_deg"], 2)))
                rec = (m["peak_abs_pa"], coh, m["pldb"], m["csel"], b["Y"], b["apex"],
                       b["delay_s"], b)
                if dy == DY and (best is None or rec[0] > best[0]):
                    best = rec
                summ.setdefault((key, dy), [0.0, 0.0, -np.inf, -np.inf])
                s = summ[(key, dy)]
                s[0] = max(s[0], m["peak_abs_pa"])
                s[1] = max(s[1], coh)
                s[2] = max(s[2], m["pldb"] if np.isfinite(m["pldb"]) else -np.inf)
                s[3] = max(s[3], m["csel"])
        summ[(key, "best")] = best
    return summ


def plot(all_out, all_summ, path):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(2, 2, figsize=(13, 9.5))
    # (a) ray paths
    a = ax[0, 0]
    ac = AIRCRAFT["Concorde M2.0"]
    for atmn, col in (("calm", "C0"), ("winter tail", "C3")):
        med = all_out[("Concorde M2.0", atmn)]["med"]
        for q in (100, 112, 125, 132, 140, 150, 165):
            tr = S.trace([math.radians(q)], ac.V, ac.altitude, med, dt=1.0, max_ground=1,
                         record=True)
            R = tr["rec"][:, :, 0]
            a.plot(R[:, 1] / 1e3, R[:, 2] / 1e3, color=col, lw=0.8, alpha=0.8)
        a.plot([], [], color=col, label=f"{atmn} (upward phi 100-165 deg)")
    a.set_xlim(-20, 450)
    a.set_ylim(0, 170)
    a.set_xlabel("lateral range Y [km]")
    a.set_ylabel("altitude [km]")
    a.set_title("(a) Concorde M2 at 20.25 km: upward rays (cross-plane projection)")
    a.legend(fontsize=8)
    # (b) ground peak vs Y, Concorde, two atmospheres
    for col_i, atmn in enumerate(("calm", "winter tail")):
        a = ax[0, 1] if col_i == 0 else ax[1, 0]
        out = all_out[("Concorde M2.0", atmn)]
        sty = {"U:direct": ("C0", "-", "unshielded: direct"),
               "U:indirect": ("C1", "-", "unshielded: indirect"),
               "S:direct": ("C0", ":", "shielded: direct"),
               "S:shield": ("C3", "-", "shielded: shield-reflected"),
               "S:transmit": ("C2", "--", "shielded: transmitted-indirect")}
        for key, (c, ls, lab) in sty.items():
            bl = out["bins"][key][DY]
            if not bl:
                continue
            byY = {}
            for b in bl:
                pk = float(np.max(np.abs(b["p"])))
                byY[b["Y"]] = max(byY.get(b["Y"], 0.0), pk)
            Ys = np.array(sorted(byY))
            P = np.array([byY[y] for y in Ys]) / S.PSF
            a.semilogy(Ys / 1e3, P, color=c, ls=ls, marker=".", ms=3, lw=1, label=lab)
        a.axhline(LIMIT_PSF, color="k", lw=1)
        a.text(5, LIMIT_PSF * 1.15, "NPRM limit 0.11 psf", fontsize=8)
        a.set_xlabel("lateral range Y [km]")
        a.set_ylabel("secondary-boom peak at ground [psf]")
        a.set_title(f"({'b' if col_i == 0 else 'c'}) Concorde M2, {atmn}: binned ground peak")
        a.set_ylim(1e-5, 1.0)
        a.grid(True, which="both", alpha=0.3)
        a.legend(fontsize=7, loc="lower right")
    # (d) waveforms at the worst bins (Concorde winter tail)
    a = ax[1, 1]
    summ = all_summ[("Concorde M2.0", "winter tail")]
    for key, c in (("U:direct", "C0"), ("U:indirect", "C1"), ("S:shield", "C3"),
                   ("S:transmit", "C2")):
        best = summ.get((key, "best"))
        if best is None:
            continue
        b = best[-1]
        t = b["t"] - b["t"][np.argmax(np.abs(b["p"]))]
        a.plot(t, b["p"], color=c, lw=1,
               label=f"{key} Y={b['Y'] / 1e3:.0f} km, apex {b['apex'] / 1e3:.0f} km")
    a.set_xlim(-8, 8)
    a.set_xlabel("time re peak [s]")
    a.set_ylabel("p at ground [Pa] (x1.9 incl.)")
    a.set_title("(d) worst-bin waveforms, Concorde M2, winter tail")
    a.legend(fontsize=7)
    a.grid(alpha=0.3)
    fig.tight_layout()
    fig.savefig(path, dpi=130)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--quick", action="store_true")
    ap.add_argument("--cache", default=os.path.join(tempfile.gettempdir(), "boomsim_secondary"))
    ap.add_argument("--only", default=None, help="aircraft|atmosphere to run (debug)")
    args = ap.parse_args()
    os.makedirs(args.cache, exist_ok=True)
    runs = RUNS
    if args.only:
        runs = [r for r in RUNS if "|".join(r) == args.only]
    all_out, all_summ, csv_rows = {}, {}, []
    for acn, atmn in runs:
        print(f"== {acn} / {atmn}", flush=True)
        out = run_case(acn, atmn, args.cache, quick=args.quick)
        all_out[(acn, atmn)] = out
        all_summ[(acn, atmn)] = summarise(acn, atmn, out, csv_rows)
        csv_rows.extend(out["ray_rows"])
        s = all_summ[(acn, atmn)]
        for key in ("U:direct", "U:indirect", "S:direct", "S:shield", "S:transmit"):
            v = s.get((key, DY))
            b = s.get((key, "best"))
            if v is None:
                print(f"   {key:11s}: no returns")
                continue
            sens = [s.get((key, dy), [np.nan])[0] / S.PSF for dy in (1500.0, 6000.0)]
            print(f"   {key:11s}: max peak {v[0] / S.PSF:.4f} psf ({v[0]:.3f} Pa) "
                  f"coherent {v[1] / S.PSF:.4f} psf, PLdB {v[2]:.1f}, CSEL {v[3]:.1f}; "
                  f"at Y {b[4] / 1e3:.0f} km apex {b[5] / 1e3:.0f} km delay {b[6] / 60:.1f} min; "
                  f"dY 1.5/6 km: {sens[0]:.4f}/{sens[1]:.4f} psf", flush=True)
    if args.only is None and not args.quick:
        keys = sorted({k for r in csv_rows for k in r}, key=lambda k: (k != "record", k))
        order = ["record", "aircraft", "atmosphere", "config", "family", "Y_km", "phi_deg",
                 "peak_pa", "peak_abs_pa", "peak_psf", "coherent_peak_psf", "pldb", "csel",
                 "asel", "energy_pa2s", "n_branches", "apex_km", "z_ref_km", "delay_min",
                 "E_h", "max_dp_over_p", "n_caustic_up", "n_caustic_down"]
        keys = order + [k for k in keys if k not in order]
        with open(os.path.join(RES, "secondary_booms.csv"), "w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=keys)
            w.writeheader()
            for r in csv_rows:
                w.writerow(r)
        plot(all_out, all_summ, os.path.join(RES, "secondary_booms.png"))
    with open(os.path.join(args.cache, "summary.pkl"), "wb") as fh:
        pickle.dump({k: {kk: (vv[:-1] if kk[1] == "best" and vv is not None else vv)
                         for kk, vv in v.items()} for k, v in all_summ.items()}, fh)


if __name__ == "__main__":
    main()
