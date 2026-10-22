"""BoomShield lid fatigue, step 2-3: strain-life by material and detail -> allowable passes/day.

Inputs: lid_response.csv (lid_response.py).  Everything else is an explicit assumption
below; [K] = handbook-class value (author's background knowledge, not a cited source),
[S#] = search-result summary in sources_notes.md (unverified).

Strain-life model (one curve per material family, at 217 K):
  eps_e   strain amplitude survived for 1e7 cycles at room temperature          [K]
  eps_e217 = eps_e * k_T(case) * k_age(case)   cold (-56 C) and end-of-life UV/ozone knock-downs
  eps >= eps_e217 : Coffin-Manson-type line through (eps_e217, 1e7) and the flex anchor
                    (eps_fold, N_fold) taken from Gelbo / MIT-fold / D2176 class data
  eps <  eps_e217 : N = 1e7 (eps_e217/eps)^m, m = 10: no hard endurance limit (conservative;
                    PE does show one [S8])
Miner sum per detail: SF * (n_boom/N_boom + n_weather/N_weather) <= 1 over the design life.
Allowable passes/day = largest n_boom per day meeting that.
Writes fatigue_details.csv (strain per detail), allowable_traffic.csv (passes/day family),
cycle_census.csv (weather vs boom cycles per year).
"""
from __future__ import annotations

import math
import os

import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
CASES = ["low", "base", "high"]          # pessimistic / central / optimistic
SF_LIFE = 10.0                           # scatter factor on life for polymers [K]
LIFE_YR = {"30 yr": 30.0, "20 yr": 20.0}
N_MAX = 1e15
CAP_SHOW = 1e6                           # report >= 1e6 passes/day as "unlimited by fatigue"
# 1-D linear response ignores local 3-D concentrations (edge diffraction, focusing at ramps,
# heavier aircraft, nonlinear slack/taut snatch): multiply stress-derived strains by this
STRAIN_MULT = {"low": 3.0, "base": 1.5, "high": 1.0}
GEOMETRIC = ("Slack partition", "Rolling-lobe", "Joint flap, light")   # strain set by bend radius

# --------------------------------------------------------------------- materials [K]
# eps_e: 1e7-cycle strain amplitude at 23 C; (eps_fold, N_fold): fold/flex anchor;
# k_T: -56 C knock-down on eps_e (low, base, high); k_age: end-of-life UV/ozone knock-down
# with a UV-blocking top skin (low, base, high).
MATERIALS = {
    "LLDPE balloon film": dict(eps_e=0.005, eps_fold=0.20, N_fold=3e3, k_T=(0.35, 0.6, 1.0),
                               k_age=(0.5, 0.8, 1.0),
                               note="PE endurance ~25-35 % of yield [K][S8]; Gelbo pinholes from ~1e3-1e4 full flexes [S2]; "
                                    "balloon PE turns brittle/notch-sensitive near -60..-90 C [S9]"),
    "PE/PET laminate": dict(eps_e=0.005, eps_fold=0.20, N_fold=3e2, k_T=(0.5, 0.8, 1.0),
                            k_age=(0.5, 0.8, 1.0),
                            note="PET skin: no true endurance limit [S8]; poor Gelbo flex (pinholes in 1e2-1e3) [K]; UV -50 % [S12]"),
    "UHMWPE scrim (Dyneema)": dict(eps_e=0.010, eps_fold=0.10, N_fold=3e4, k_T=(0.8, 1.0, 1.0),
                                   k_age=(0.6, 0.8, 1.0),
                                   note="fibre fail strain ~3.5 %; tension fatigue good [S7]; ductile in cold [K]; creep is the static limit"),
    "Vectran scrim": dict(eps_e=0.008, eps_fold=0.10, N_fold=1e5, k_T=(0.8, 1.0, 1.0),
                          k_age=(0.4, 0.7, 1.0),
                          note="D2176 yarn 115k cycles vs aramid 1.4-15k; crease x100 -0.8 % [S11]; UV-sensitive, needs cover"),
    "Aramid scrim": dict(eps_e=0.008, eps_fold=0.10, N_fold=3e3, k_T=(0.8, 1.0, 1.0),
                         k_age=(0.3, 0.6, 0.9),
                         note="kink-band fold damage (crease x100 -23 %) [S11][S10]; UV-poor"),
    "Coated fabric (PU/PVF on polyester)": dict(eps_e=0.004, eps_fold=0.10, N_fold=1e4, k_T=(0.25, 0.5, 0.8),
                                                k_age=(0.4, 0.7, 0.9),
                                                note="TARS envelopes last 5-6 yr at ground [S4]; PU coating stiff/cracks near Tg (-40..-60 C) [K]"),
}


def eps_endurance(mat, case):
    i = CASES.index(case)
    m = MATERIALS[mat]
    return m["eps_e"] * m["k_T"][i] * m["k_age"][i]


def cycles_to_failure(eps, mat, case, m_hc=10.0):
    if eps <= 0:
        return N_MAX
    ee = eps_endurance(mat, case)
    M = MATERIALS[mat]
    if eps >= ee:
        k = math.log(1e7 / M["N_fold"]) / math.log(M["eps_fold"] / ee)
        return max(1.0, 1e7 * (eps / ee) ** (-k))
    return min(N_MAX, 1e7 * (ee / eps) ** m_hc)


# ---------------------------------------------------------------------- details
E_COLD = 1.0e9          # Pa, PE laminate at 217 K, low end (higher strain for a given stress)


def detail_strains(r):
    """Strain amplitude per pass (and per weather cycle) at each detail, for one lid_response row.
    Returns {detail: dict(low=, base=, high=, w_eps=(...), w_n_yr=(...), note)} where w_* are
    the weather-cycle strain and count per year at that detail (low, base, high)."""
    t = r.t_lid_mm / 1e3
    eps_field = r.eps_membrane + r.eps_bending
    F = r.tie_force_pp_n / 2 if r.tie_K > 0 else 0.0          # tie force amplitude
    F = max(F, 1.0)
    pad = lambda D, kt: kt * F / (math.pi * D * t) / E_COLD
    d = {}
    d["Lid field (sheet)"] = dict(v=(eps_field,) * 3, w=(0, 0, 0), n=(0, 0, 0),
                                  note="membrane 0.5 slope^2 + bending kappa t/2")
    d["Heat-sealed lap seam"] = dict(v=(4 * eps_field, 3 * eps_field, 2 * eps_field), w=(0, 0, 0), n=(0, 0, 0),
                                     note="Kt 2-4 at the stiffness step")
    d["Tie fitting, 50 mm grommet"] = dict(v=(pad(0.05, 3.0),) * 3,
                                           w=(pad(0.05, 3.0) * 25, pad(0.05, 3.0) * 13, pad(0.05, 3.0) * 13),
                                           n=(365, 365, 365),
                                           note="Kt 3; daily breathing puts 0-27 Pa/tie (13-25x the boom swing)")
    d["Tie fitting, 0.5 m spread pad"] = dict(v=(pad(0.5, 2.0),) * 3,
                                              w=(pad(0.5, 2.0) * 25, pad(0.5, 2.0) * 13, pad(0.5, 2.0) * 13),
                                              n=(365, 365, 365), note="Kt 2, doubler/load-spreading patch")
    # slack partition / bellows: 50 um film; crease radius sets strain t/2R; the gas carries the
    # partition top 0.4-1 m sideways each pass while the hem is held by the lid.
    tp = 50e-6
    R = (2e-3, 10e-3, 50e-3)
    d["Slack partition hem/crease (50 um)"] = dict(v=tuple(tp / (2 * x) for x in R),
                                                   w=tuple(tp / (2 * x) for x in (1e-3, 4e-3, 20e-3)),
                                                   n=(365 + 2e4, 365 + 1e4, 365),
                                                   note="crease R 2/10/50 mm per pass; breathing folds tighter (R 1/4/20 mm), "
                                                        "365/yr plus gravity-wave strokes 1-2e4/yr (low/base)")
    d["Rolling-lobe partition / rim (R 0.5 m)"] = dict(v=(tp / 1.0,) * 3, w=(tp / 1.0,) * 3, n=(365 + 2e4,) * 3,
                                                       note="pre-formed lobe fixes the bend radius (air-spring style)")
    # joint flap lighter than the lid: the lid is mass-controlled (net load = inertia), so a light
    # flap over a gap L bulges against its own tension: sag = p L^2 / 8N; hinge strain ~ t_f/(2R)
    pf, L, N = r.net_lid_load_pa, 1.0, 200.0
    sag = min(pf * L ** 2 / (8 * N), 0.5)
    tf = 0.5e-3
    d["Joint flap, light (0.5 kg/m2, 1 m)"] = dict(v=(tf / (2 * 5e-3), tf / (2 * 10e-3), tf / (2 * 25e-3)),
                                                   w=(0, 0, 0), n=(0, 0, 0),
                                                   note=f"sag ~{sag:.2f} m per pass -> hinge fold R 5/10/25 mm")
    d["Joint flap, mass-matched"] = dict(v=(3 * eps_field,) * 3, w=(0, 0, 0), n=(0, 0, 0),
                                         note="same areal mass as the lid: moves with it; seam-like")
    d["Tether/tendon attachment"] = dict(v=(1e-5,) * 3, w=(2e-3, 1.5e-3, 1e-3), n=(1e3, 3e2, 1e2),
                                         note="boom ~0.1 % of 26 MN working load; wind/gust cycles dominate")
    return d


def life_rows(resp):
    rows, allow = [], []
    for _, r in resp.iterrows():
        tag = f"{r.aircraft} M{r.mach:g}, H2 {r.h2_m:g} m, lid {r.lid_kg_m2:g} kg/m2"
        for det, dd in detail_strains(r).items():
            for mat in MATERIALS:
                # material/detail pairs that make sense
                if det.startswith("Slack partition") or det.startswith("Rolling-lobe"):
                    if mat not in ("LLDPE balloon film", "PE/PET laminate"):
                        continue
                if det.startswith("Tether"):
                    if mat not in ("UHMWPE scrim (Dyneema)", "Vectran scrim", "Aramid scrim"):
                        continue
                for ci, case in enumerate(CASES):
                    sm = 1.0 if det.startswith(GEOMETRIC) else STRAIN_MULT[case]
                    eps = dd["v"][ci] * sm
                    Nb = cycles_to_failure(eps, mat, case)
                    Nw = cycles_to_failure(dd["w"][ci] * sm, mat, case) if dd["n"][ci] else N_MAX
                    rows.append(dict(load_case=tag, detail=det, material=mat, case=case, eps_boom=eps,
                                     N_boom=Nb, eps_weather=dd["w"][ci] * sm, weather_cycles_yr=dd["n"][ci],
                                     N_weather=Nw, note=dd["note"]))
                    for lname, L in LIFE_YR.items():
                        wd = SF_LIFE * dd["n"][ci] * L / Nw          # weather damage share
                        spare = 1.0 - wd
                        perday = 0.0 if spare <= 0 else spare * Nb / (SF_LIFE * 365 * L)
                        allow.append(dict(load_case=tag, detail=det, material=mat, case=case, life=lname,
                                          eps_boom=eps, N_boom=Nb, weather_damage=wd,
                                          passes_per_day=min(perday, CAP_SHOW),
                                          passes_per_day_boom_only=min(Nb / (SF_LIFE * 365 * L), CAP_SHOW)))
    return pd.DataFrame(rows), pd.DataFrame(allow)


DESIGNS = {
    "A. as sketched": ["Lid field (sheet)", "Heat-sealed lap seam", "Joint flap, light (0.5 kg/m2, 1 m)",
                       "Tie fitting, 50 mm grommet", "Slack partition hem/crease (50 um)", "Tether/tendon attachment"],
    "B. matched flaps, grommet ties, slack partitions": ["Lid field (sheet)", "Heat-sealed lap seam",
                       "Joint flap, mass-matched", "Tie fitting, 50 mm grommet",
                       "Slack partition hem/crease (50 um)", "Tether/tendon attachment"],
    "C. fatigue-aware": ["Lid field (sheet)", "Heat-sealed lap seam", "Joint flap, mass-matched",
                         "Tie fitting, 0.5 m spread pad", "Rolling-lobe partition / rim (R 0.5 m)",
                         "Tether/tendon attachment"],
}
LID_MATERIALS = ["LLDPE balloon film", "PE/PET laminate", "UHMWPE scrim (Dyneema)", "Vectran scrim",
                 "Aramid scrim", "Coated fabric (PU/PVF on polyester)"]


def design_caps(allow):
    """System cap = min over the design's details; lid material for lid details, LLDPE film for
    partitions, UHMWPE for tether attachments."""
    out = []
    for (lc, life), a in allow.groupby(["load_case", "life"]):
        for dname, dets in DESIGNS.items():
            for lm in LID_MATERIALS:
                for case in CASES:
                    vals, vals_b, lim = [], [], None
                    for det in dets:
                        mat = ("LLDPE balloon film" if det.startswith(("Slack partition", "Rolling-lobe"))
                               else "UHMWPE scrim (Dyneema)" if det.startswith("Tether") else lm)
                        x = a[(a.detail == det) & (a.material == mat) & (a.case == case)]
                        vals.append((float(x.passes_per_day.iloc[0]), det))
                        vals_b.append(float(x.passes_per_day_boom_only.iloc[0]))
                    v, lim = min(vals)
                    out.append(dict(load_case=lc, life=life, design=dname, lid_material=lm, case=case,
                                    cap_passes_day=v, limiting_detail=lim if v < CAP_SHOW else "none (>=1e6)",
                                    cap_boom_only=min(vals_b)))
    return pd.DataFrame(out)


def cycle_census():
    """Large-cycle sources per year at a point on the lid (screening estimates)."""
    return pd.DataFrame([
        dict(source="Boom pass, 100 flights/day", cycles_yr=3.65e4, lid_load_pa="300-650 (doubled)",
             thickness_stroke_m="0.1-0.3", note="one reversal per pass; lid strain ~1e-5"),
        dict(source="Boom pass, 4,400 flights/day (trunk, full capture)", cycles_yr=1.6e6,
             lid_load_pa="300-650", thickness_stroke_m="0.1-0.3", note=""),
        dict(source="Diurnal breathing", cycles_yr=365, lid_load_pa="0-27 through ties",
             thickness_stroke_m="~5 (55 m layer)", note="largest stroke: partitions, ties, rim fold fully"),
        dict(source="Gravity / mountain waves (periods 4-60 min)", cycles_yr="1e4-5e4",
             lid_load_pa="few", thickness_stroke_m="0.1-0.5 (est.)",
             note="T' ~ 1 K per 100 m displacement -> 0.5 % volume; strong over the Rockies [S16]"),
        dict(source="Gusts / turbulence (quasi-steady lift 5-22 Pa)", cycles_yr="1e5-1e6", lid_load_pa="5-22 over km",
             thickness_stroke_m="<0.05", note="added mass ~70 kg/m2 makes the layer ride gusts (structure §3)"),
        dict(source="Synoptic wind events (side force 0.5-2 kN/m)", cycles_yr="50-100", lid_load_pa="-",
             thickness_stroke_m="-", note="tethers, rim tendons; not a lid-face cycle"),
    ])


def main():
    resp = pd.read_csv(os.path.join(HERE, "lid_response.csv"))
    resp = resp[resp.tie_K > 0].reset_index(drop=True)
    det, allow = life_rows(resp)
    det.to_csv(os.path.join(HERE, "fatigue_details.csv"), index=False, float_format="%.3g")
    allow.to_csv(os.path.join(HERE, "allowable_traffic.csv"), index=False, float_format="%.3g")
    cycle_census().to_csv(os.path.join(HERE, "cycle_census.csv"), index=False)
    dc = design_caps(allow)
    dc.to_csv(os.path.join(HERE, "design_caps.csv"), index=False, float_format="%.3g")
    with pd.option_context("display.width", 250, "display.max_rows", 400, "display.max_colwidth", 40):
        for lc in dc.load_case.unique():
            x = dc[(dc.load_case == lc) & (dc.life == "30 yr")]
            print("==", lc)
            print(x.pivot_table(index=["design", "lid_material"], columns="case", values="cap_passes_day")[CASES]
                  .map(lambda v: f"{v:.3g}"))
            print(x.pivot_table(index=["design", "lid_material"], columns="case", values="cap_boom_only")[CASES]
                  .map(lambda v: f"{v:.3g}"))
            print(x[x.case != "high"].groupby(["design", "case"]).limiting_detail.agg(lambda s: sorted(set(s))))
    return
    base = allow[(allow.life == "30 yr") & allow.load_case.str.startswith("Concorde M2, H2 55 m, lid 4.45")]
    pv = base.pivot_table(index=["detail", "material"], columns="case", values="passes_per_day")[CASES]
    with pd.option_context("display.width", 200, "display.max_rows", 200):
        print(pv.map(lambda x: f"{x:.3g}"))
        for lc in allow.load_case.unique():
            a = allow[(allow.life == "30 yr") & (allow.load_case == lc)]
            print(lc)
            print(a.pivot_table(index="detail", columns="case", values="passes_per_day", aggfunc="min")[CASES]
                  .map(lambda x: f"{x:.3g}"))


if __name__ == "__main__":
    main()
