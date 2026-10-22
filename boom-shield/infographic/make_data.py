"""BoomShield infographic data driver.

  python boom-shield/infographic/make_data.py --run [--procs 3] [--cache DIR]
      Boom Overture (M1.7) preset through the run_3d.py pipeline (same numerics:
      R_m 3 km all azimuths, near-track re-run at 8 km; f_max 100 Hz, df 1/6 Hz,
      50 kHz ground propagation), unshielded and W 1 km H2 30 m slab, h 250 m.
      Writes <cache>/overture_summary.json. Does not touch results/3d_cache.
  python boom-shield/infographic/make_data.py --build [--cache DIR]
      Assembles boom-shield/infographic/data.json from the repo results, the Overture
      summary and the cited reference tables below.
"""
from __future__ import annotations

import argparse
import dataclasses
import json
import math
import os
import sys
import tempfile

HERE = os.path.dirname(os.path.abspath(__file__))
SIM = os.path.join(os.path.dirname(HERE), "sim")
sys.path.insert(0, SIM)

import run_3d as R3  # noqa: E402
from boomsim.source import Aircraft  # noqa: E402

DEFAULT_CACHE = os.environ.get("OVERTURE_CACHE", os.path.join(tempfile.gettempdir(), "boomshield_overture_cache"))

# Boom Overture preset (estimates, see data_notes.md):
#  length 201 ft = 61.3 m (Boom spec sheet); cruise M 1.7 (Boom);
#  MTOW ~77 t (pre-2022 published figure, unverified for the refined design);
#  mid-cruise mass ~0.81 MTOW (same fraction as the repo's Concorde preset 150/185 t)
#  -> 62 t -> W = 0.61 MN; volume ~300 m^3 (Concorde 450 m^3 scaled for a ~2.7 m
#  fuselage and ~40 t fuel); lift on the wing, ~0.40 L ... 0.90 L.
OVERTURE = Aircraft("Overture", length=61.3, mach=1.7, altitude=18300.0,
                    weight=62000.0 * 9.81, volume=300.0,
                    lift_start=0.40, lift_end=0.90)

_aircraft0 = R3.aircraft


def aircraft(name, mach, h):
    if name == "Overture":
        return dataclasses.replace(OVERTURE, mach=float(mach), altitude=R3.H_S + float(h))
    return _aircraft0(name, mach, h)


R3.aircraft = aircraft


def overture_configs():
    base = dict(group="baseline", acn="Overture", mach=1.7, h=250, W=None, fill=None, T=0.0,
                profile="slab", fringe=0.0, cells=0.0, tieK=0.0, ties="line", R=R3.R_MAIN,
                phis=None)
    C = [dict(base, key="O1.7_h250_none", label="Overture M1.7 h250: unshielded")]
    sh = dict(base, group="aircraft", W=1000.0, fill="hydrogen", T=30.0,
              key="O1.7_h250_W1000_H230", label="Overture M1.7: W 1000 m, H2 30 m")
    C.append(sh)
    C.append(dict(sh, key=sh["key"] + "_R8000", group="far", R=R3.R_FAR, phis=R3.PHI_SENS_FAR,
                  label=sh["label"] + " [R_m 8 km]"))
    for c in C:
        c["base"] = "O1.7_h250_none"
    return C


def run(procs, cache):
    R3.CACHE = cache
    R3.WORKERS = procs
    C = overture_configs()
    for c in C:
        print(f"== {c['key']}", flush=True)
        R3.ground(c, procs)
    summ = summarize_overture(cache)
    summ.update(cruise_points(cache))
    with open(os.path.join(cache, "overture_summary.json"), "w") as fh:
        json.dump(summ, fh, indent=1)
    print(json.dumps({k: v for k, v in summ.items() if k != "O1.7_h250_W1000_H230"}, indent=1)[:4000])


def cruise_points(cache):
    """Unshielded undertrack levels at each type's normal cruise altitude (direct
    near-field path, propagate.ground_carpet, phi = 0 only): Overture M1.7 at 60,000 ft,
    SR-71 M3.0 at 80,000 ft (validation preset). One process."""
    from boomsim import propagate as PR
    from boomsim.source import SR71
    out = {}
    for key, ac in (("cruise_18_3km", dataclasses.replace(OVERTURE, mach=1.7, altitude=18288.0)),
                    ("sr71_cruise_24_4km", SR71)):
        p = os.path.join(cache, key + ".json")
        if os.path.exists(p):
            out[key] = json.load(open(p))
            continue
        res = PR.ground_carpet(ac, [0.0], start="near_field", procs=1, keep=False)
        d = dict(altitude_m=ac.altitude, mach=ac.mach, pldb_raw=float(res["pldb"][0]),
                 pldb_proxy=float(res["pldb_proxy"][0]), peak_pa=float(res["peak_pa"][0]),
                 peak_psf=float(res["peak_pa"][0]) / 47.88)
        json.dump(d, open(p, "w"))
        out[key] = d
    return out


def _clean(v):
    if isinstance(v, float) and (math.isnan(v) or math.isinf(v)):
        return None
    if hasattr(v, "item"):
        return _clean(v.item())
    if isinstance(v, dict):
        return {str(k): _clean(x) for k, x in v.items()}
    return v


def summarize_overture(cache):
    R3.CACHE = cache
    C = {c["key"]: c for c in overture_configs()}
    base = R3.carpet_rows(C["O1.7_h250_none"])
    out = {}
    for k in ("O1.7_h250_none", "O1.7_h250_W1000_H230"):
        rows = R3.carpet_rows(C[k])
        if rows is None:
            continue
        s = R3.summarize(C[k], rows, base if k != "O1.7_h250_none" else None)
        s["per_phi"] = [dict(phi=r["phi"], y_g_km=r["y_g"] / 1e3, pldb=R3._m(r, "pldb"),
                             pldb_proxy=R3._m(r, "pldb_proxy")) for r in rows]
        out[k] = {kk: _clean(v) for kk, v in s.items()}
    ac = aircraft("Overture", 1.7, 250)
    out["preset"] = dict(name=ac.name, length_m=ac.length, mach=ac.mach, altitude_m=ac.altitude,
                         weight_N=ac.weight, volume_m3=ac.volume, lift_start=ac.lift_start,
                         lift_end=ac.lift_end)
    return out



# ============================================================================ build
REPO = os.path.dirname(os.path.dirname(HERE))
RES = os.path.join(SIM, "results")
OUT = os.path.join(HERE, "data.json")

# Source registry. status: "repo" (file in this repo), "official" (agency/manufacturer
# page seen at search-result level), "secondary" (press/third party), "unverified" = the
# figure was read from a search-engine summary/snippet; the page itself could not be
# fetched from this sandbox (WebFetch/curl blocked), so it is not independently checked.
SOURCES = {
    "r3d_sum": dict(ref="boom-shield/sim/results/3d_summary.md", status="repo"),
    "r3d_csv": dict(ref="boom-shield/sim/results/3d_carpet_summary.csv", status="repo"),
    "r3d_res": dict(ref="boom-shield/sim/results/3d_results.csv", status="repo"),
    "r3d_py": dict(ref="boom-shield/sim/run_3d.py", status="repo"),
    "src_py": dict(ref="boom-shield/sim/boomsim/source.py", status="repo"),
    "wt": dict(ref="boom-shield/sim/results/wind_turbulence.md", status="repo"),
    "at": dict(ref="boom-shield/sim/results/alongtrack_summary.md", status="repo"),
    "valc": dict(ref="boom-shield/sim/results/validation_carpet.csv", status="repo"),
    "valp": dict(ref="boom-shield/sim/results/validation_propagation.csv", status="repo"),
    "struct": dict(ref="boom-shield/structure/structure-concept.md", status="repo"),
    "vac": dict(ref="boom-shield/vacuum/vacuum-balloon-review.md", status="repo"),
    "handoff": dict(ref="boom-shield/handoff/2026-10-20-boomshield-handoff.md", status="repo"),
    "addendum": dict(ref="boom-shield/handoff/2026-10-20-review-addendum.md", status="repo"),
    "ovr_run": dict(ref="boom-shield/infographic/make_data.py --run (this study; run_3d.py pipeline)",
                    status="repo"),
    "ovr_preset": dict(ref="boom-shield/infographic/make_data.py (OVERTURE preset; estimates)",
                       status="repo"),
    # web
    "boom_specs": dict(ref="https://boomsupersonic.com/press-release/refined-overture-northrop-grumman",
                       status="official, unverified (search-result level)",
                       note="Overture: length 201 ft, span 106 ft, cruise M1.7 over water, 4,250 nm, 64-80 pax"),
    "fi_overture": dict(ref="https://flightplan.forecastinternational.com/2025/08/25/boom-supersonic-overture-airliner-program-outlook/",
                        status="secondary, unverified", note="cruise M1.7, altitudes up to 60,000 ft"),
    "mtow77": dict(ref="https://aerocorner.com/aircraft/boom-overture/ ; https://knaviation.net/boom-overture/",
                   status="secondary, unverified",
                   note="MTOW ~77 t (170,000 lb); figure dates from the pre-2022 design, Boom publishes none for the refined design"),
    "boomless_pr": dict(ref="https://boomsupersonic.com/press-release/boomless-cruise-announcement",
                        status="official, unverified (search-result level)",
                        note="Feb 2025: XB-1 M1.12 with no boom reaching ground microphones; Overture up to M1.3 overland without audible boom"),
    "boomless_faq": dict(ref="https://boomsupersonic.com/boomless-cruise",
                         status="official, unverified (search-result level)",
                         note="typical Boomless Cruise M1.1-1.2, up to 1.3; at M1.7 over water Overture makes a boom"),
    "x59_target_ntrs": dict(ref="https://ntrs.nasa.gov/citations/20205002229",
                            status="official, unverified (search-result level)",
                            note="Doebler & Rathsam, 'How loud is X-59's shaped sonic boom?', POMA 36, 040005 (2019): design target 75 dB Stevens PL; 'acoustic thermometer' of familiar impulsive sounds"),
    "x59_carpet": dict(ref="https://www.nasa.gov/centers-and-facilities/armstrong/nasa-instruments-will-listen-for-supersonic-x-59s-quiet-thump/",
                       status="official, unverified (search-result level)",
                       note="X-59 thump ~75 PLdB; Concorde 'over 100' PLdB (L. Cliatt)"),
    "x59_cruise": dict(ref="https://www.nasa.gov/blogs/quesst/2026/06/12/x-59-blog-061226/",
                       status="official, unverified (search-result level)",
                       note="target cruise M1.4 (~925 mph) at 55,000 ft reached 12 Jun 2026"),
    "x59_len": dict(ref="https://ntrs.nasa.gov/api/citations/20220009074/downloads/QuesstMissOvuaX59Pres_050822r.pdf",
                    status="official, unverified (search-result level)", note="length 99 ft, span 29.5 ft"),
    "rathsam25": dict(ref="https://ntrs.nasa.gov/api/citations/20250001048/downloads/Rathsam_RAeS_Loughborough_Branch_STRIVES.pdf",
                      status="official, unverified (search-result level)",
                      note="X-59 thump ~ distant thunder or car door slam from across the street; Concorde 105 dB PL ~ nearby thunder or car door slam from inside the car"),
    "nasa_podcast": dict(ref="https://www.nasa.gov/podcasts/houston-we-have-a-podcast/the-quesst-for-low-sonic-booms/",
                         status="official, unverified (search-result level)",
                         note="L. Ozoroski: 'more like a car door closing, you know, across the street at your neighbor's house'"),
    "phys_today": dict(ref="https://physicstoday.aip.org/news/nasa-unveils-a-supersonic-plane-with-a-quiet-boom",
                       status="secondary, unverified",
                       note="X-59 boom perceived about as loud as a car door slamming 6 m away"),
    "nasa_fs016": dict(ref="https://www.nasa.gov/centers/dryden/pdf/120274main_FS-016-DFRC.pdf",
                       status="official, unverified (table values via secondary copies)",
                       note="Concorde 1.94 psf at M2, 52,000 ft; SR-71 0.9 psf at M3, 80,000 ft"),
    "niosh_meter": dict(ref="https://stacks.cdc.gov/view/cdc/78630/cdc_78630_DS1.pdf",
                        status="official (CDC/NIOSH), unverified (search-result level)",
                        note="whisper 30 dB, normal conversation 60 dB, hand drill 98, chain saw 110; REL 85 dBA"),
    "nidcd": dict(ref="https://www.nidcd.nih.gov/health/noise-induced-hearing-loss",
                  status="official (NIH/NIDCD), unverified (search-result level; values via NIDCD-derived copies)",
                  note="conversation 60-70 dBA; lawn mowers 80-100 dBA; motorcycles ~95; sirens 110-129 dBA; concerts/sports 94-110"),
    "cdc_listen": dict(ref="https://www.cdc.gov/nceh/hearing_loss/infographic/",
                       status="official (CDC), unverified (search-result level)",
                       note="leaf blower 90 dB; siren 120 dB; motorcycle 95 dB"),
    "caltrans_tens": dict(ref="https://dot.ca.gov/-/media/dot-media/programs/environmental-analysis/documents/env/tens-sep2013-a11y.pdf",
                          status="state agency; values read from local reproductions (e.g. https://weblink.cityofpleasantonca.gov/WebLink/0/doc/261418/Page23.aspx), unverified",
                          note="jet fly-over at 300 m (1,000 ft) 110 dBA; gas lawn mower at 1 m 100 dBA; diesel truck at 15 m (50 mph) 90 dBA"),
    "ne_dot": dict(ref="https://dot.nebraska.gov/media/eyablfmc/noise-study-7-2013.pdf",
                   status="state agency, unverified (search-result level)",
                   note="jet flyover at 985 ft 105 dBA; diesel truck at 50 ft 85 dBA"),
    "wsdot": dict(ref="https://appstest.wsdot.wa.gov/sites/default/files/2021-11/I5-JBLM-Vicinity-FINALMay2017EA-Section45Noise.pdf",
                  status="state agency, unverified (search-result level)",
                  note="highway traffic typically 70 dBA at 50 ft"),
    "spokane": dict(ref="https://www.spokanecounty.gov/DocumentCenter/View/7680",
                    status="local agency, unverified (search-result level)",
                    note="light auto traffic at 50 ft 60 dBA; freeway traffic at 50 ft 70 dBA"),
    "cpuc_eir": dict(ref="https://ia.cpuc.ca.gov/environment/info/aspen/westofdevers/deir/d13_noise.pdf",
                     status="state agency EIR, unverified (search-result level)",
                     note="jet takeoff at 200 ft 120 dBA; rock music concert at 50 ft 110 dBA"),
    "dog_uk": dict(ref="https://moderngov.stockton.gov.uk/Data/Planning%20Committee/201907031330/Agenda/$att36456.doc.pdf ; https://cdn.centralcoast.nsw.gov.au/sites/default/files/2023-07/addendum_to_noise_report_-_renl417-03f02_technical_memo_r0_1.pdf",
                   status="local-authority acoustic reports, unverified (search-result level)",
                   note="dog barking at ~1 m: LA5 75 dBA (one dog) to LAeq 95.9 dB (two dogs)"),
    "faa_basics": dict(ref="https://www.faa.gov/sites/faa.gov/files/air_traffic/environmental_issues/ared_documentation/Appendix_E_Basics_of_Noise_1.pdf",
                       status="official (FAA), unverified (search-result level)",
                       note="+10 dB heard as ~twice as loud; 5 dB change noticeable; everyday sounds 30-100 dB"),
}


def V(v, src, tag="simulated", **kw):
    """One number with its provenance."""
    d = dict(value=v, src=src if isinstance(src, list) else [src], tag=tag)
    d.update(kw)
    return d


def r1(x):
    return None if x is None or (isinstance(x, float) and (math.isnan(x) or math.isinf(x))) else round(float(x), 1)


def load_summary_csv():
    import csv
    with open(os.path.join(RES, "3d_carpet_summary.csv")) as fh:
        rows = {r["key"]: r for r in csv.DictReader(fh)}

    def f(key, col):
        v = rows[key][col]
        return float(v) if v not in ("", None) else float("nan")
    return f


def sim_entry(f, base, sh, far, src, cutoff_km, ovr=None):
    """Absolute levels (proxy = headline, raw also given) and headline deltas from the
    3-D summary rows (or the Overture summary dict, same column names)."""
    g = (lambda k, c: ovr[k].get(c)) if ovr else f
    un_ut_p, un_ut_r = g(base, "ut_pldb_proxy"), g(base, "ut_pldb")
    un_mx_p, un_mx_r = g(base, "max_pldb_proxy"), g(base, "max_pldb")
    if ovr:
        fr = ovr[sh]
        ut8_p, ut8_r, d8 = fr.get("ut_pldb_proxy_R8"), fr.get("ut_pldb_R8"), fr.get("d_ut_R8")
    else:
        ut8_p, ut8_r = f(far, "ut_pldb_proxy"), f(far, "ut_pldb")
        d8 = f(far, "d_ut")
    return dict(
        unshielded=dict(
            undertrack_pldb=V(r1(un_ut_p), src, metric="proxy (headline)", raw_pldb=r1(un_ut_r)),
            carpet_max_pldb=V(r1(un_mx_p), src, metric="proxy (headline)", raw_pldb=r1(un_mx_r),
                              note="unshielded carpet maximum is undertrack")),
        shielded=dict(
            undertrack_pldb=V(r1(ut8_p), src, metric="proxy (headline)", raw_pldb=r1(ut8_r),
                              R_m="8 km near-track re-run"),
            undertrack_pldb_R3=V(r1(g(sh, "ut_pldb_proxy")), src, metric="proxy",
                                 raw_pldb=r1(g(sh, "ut_pldb")), R_m="3 km (2-7 dB optimistic near track)"),
            carpet_max_pldb=V(r1(g(sh, "max_pldb_proxy")), src, metric="proxy (headline)",
                              raw_pldb=r1(g(sh, "max_pldb")), R_m="3 km (all azimuths)",
                              at_y_km=r1(g(sh, "y_at_max_km")),
                              note=("at the outermost computed azimuth; the carpet edge beyond it was not computed"
                                    if abs(g(sh, "y_at_max_km") - g(sh, "y_max_km")) < 0.05 else None))),
        delta_headline_db=dict(
            undertrack=V(r1(d8), src, R_m="8 km", note="less favourable of raw and proxy",
                         raw=r1(ut8_r - un_ut_r)),
            carpet_max=V(r1(g(sh, "d_max")), src, R_m="3 km", raw=r1(g(sh, "d_max_raw"))),
            turbulence_paired_note=V("headline deltas conservative by ~5.6-8.3 dB at W 1 km (Concorde M2); paired turbulent delta ~ raw delta",
                                     "wt", tag="simulated")),
        shielded_strip=V(dict(half_width_km_ge10dB=r1(g(sh, "strip_half_km")),
                              computed_to_km=r1(g(sh, "y_max_km")), carpet_half_width_km=cutoff_km),
                         src, note="contiguous half-width from the track with headline delta <= -10 dB; '>=' when it reaches the last computed azimuth"),
    )


def build(cache):
    f = load_summary_csv()
    ovr_path = os.path.join(cache, "overture_summary.json")
    ovr = json.load(open(ovr_path)) if os.path.exists(ovr_path) else None
    from boomsim import propagate as PR

    def cutoff(acn, mach):
        ac = aircraft(acn, mach, 250)
        return r1(PR.lateral_cutoff(ac.altitude, ac.V)[1] / 1e3)

    alt_note = "20.25 km = shield top 20 km + h 250 m (repo convention), not the type's normal cruise altitude"

    # ---------------------------------------------------------------- aircraft
    con = dict(name="Concorde", mach=V(2.0, "r3d_sum", "simulated"),
               length_m=V(61.66, "src_py", "published"),
               altitude_used_km=V(20.25, "r3d_sum", "simulated", note=alt_note),
               shield_case="W 1 km, H2 30 m slab, h 250 m",
               **sim_entry(f, "C2_h250_none", "C2_h250_W1000_H230", "C2_h250_W1000_H230_R8000",
                           ["r3d_csv", "r3d_sum"], cutoff("Concorde", 2.0)),
               cruise_reference=dict(
                   measured_pldb=V(105, ["rathsam25", "x59_carpet"], "published",
                                   note="NASA: Concorde ~105 PL ('over 100 PLdB')"),
                   published_overpressure_psf=V(1.94, ["nasa_fs016", "valp"], "published",
                                                note="M2 at 52,000 ft (15.8 km)"),
                   simulated_at_16_8_km=V(dict(pldb_raw=109.7, pldb_proxy=98.9, peak_psf=1.97),
                                          ["valc", "valp"], "simulated",
                                          note="unshielded undertrack at 16.8 km, M2; measured ~105 PL lies between raw and proxy")),
               note="Simulated with the repo's Concorde preset (150 t mid-cruise); measured Concorde booms were ~105 PLdB.")
    sr = dict(name="SR-71", mach=V(3.2, "r3d_sum", "simulated"),
              length_m=V(32.74, "src_py", "published"),
              altitude_used_km=V(20.25, "r3d_sum", "simulated", note=alt_note),
              shield_case="W 1 km, H2 30 m slab, h 250 m",
              **sim_entry(f, "S3.2_h250_none", "S3.2_h250_W1000_H230", "S3.2_h250_W1000_H230_R8000",
                          ["r3d_csv", "r3d_sum"], cutoff("SR-71", 3.2)),
              cruise_reference=dict(
                  published_overpressure_psf=V(0.9, ["nasa_fs016", "valp"], "published",
                                               note="M3 at 80,000 ft (24.4 km); repo model gives 0.99 psf"),
                  measured_pldb=V(None, [], "published", note="no published PLdB found"),
                  simulated_at_24_4_km=V({k: r1(v) for k, v in ovr["sr71_cruise_24_4km"].items()} if ovr and "sr71_cruise_24_4km" in ovr else None,
                                         "ovr_run", "simulated",
                                         note="unshielded undertrack, M3.0 at 80,000 ft (repo SR-71 validation preset), direct near-field path")),
              note="Highest Mach simulated; widest carpet, so the far-edge shielding (about -21 dB) is weakest.")

    if ovr and "O1.7_h250_W1000_H230" in ovr:
        ov_sim = sim_entry(None, "O1.7_h250_none", "O1.7_h250_W1000_H230", None, ["ovr_run"],
                           cutoff("Overture", 1.7), ovr=ovr)
    else:
        ov_sim = dict(pending="run make_data.py --run first")
    pre = ovr["preset"] if ovr else {}
    ov = dict(name="Boom Overture", mach=V(1.7, "boom_specs", "published", note="supersonic cruise over water"),
              length_m=V(61.3, "boom_specs", "published", note="201 ft"),
              mtow_t=V(77, "mtow77", "published", note="pre-2022 design figure; refined design unpublished -> estimate"),
              preset=V(dict(weight_N=pre.get("weight_N"), mid_cruise_mass_t=62, volume_m3=pre.get("volume_m3"),
                            lift_start=pre.get("lift_start"), lift_end=pre.get("lift_end")),
                       "ovr_preset", "extrapolated",
                       note="estimates: mass 0.81 x MTOW (Concorde preset ratio), volume scaled from Concorde; boom level scales with weight, so treat absolute PLdB as +/- several dB"),
              altitude_used_km=V(20.25, "ovr_run", "simulated", note=alt_note),
              normal_cruise_altitude_km=V(18.3, "fi_overture", "published", note="up to 60,000 ft"),
              shield_case="W 1 km, H2 30 m slab, h 250 m",
              **ov_sim,
              boomless_cruise=V(dict(mach_typical="1.1-1.2", mach_max=1.3, xb1_demo_mach=1.12),
                                ["boomless_pr", "boomless_faq"], "published",
                                note="Mach cut-off: the boom refracts upward and does not reach the ground; at M1.7 over water Overture does boom"),
              note="Same pipeline and numerics as Concorde/SR-71 (R_m 3 km all azimuths, 8 km near track); Overture's real shape is not modelled (Sears-Haack + lift).")
    if ovr and "O1.7_h250_none" in ovr and "cruise_18_3km" in ovr:
        ov["cruise_reference"] = dict(simulated_at_18_3_km=V({k: r1(v) for k, v in ovr["cruise_18_3km"].items()},
                                                             "ovr_run", "simulated",
                                                             note="unshielded undertrack at 60,000 ft, direct near-field path; the modelled front is not yet fully steepened there (~13 ms rise vs ~1 ms at 20.25 km), so it comes out ~2 dB quieter than at 20.25 km; indicative only"),
                                      published_pldb=V(None, [], "published", note="Boom publishes no Mach 1.7 boom level"))

    c14 = "C1.4_h250_W1000_He30"
    d_ut14, d_mx14 = f(c14 + "_R8000", "d_ut"), f(c14, "d_max")
    d_ut14_raw = f(c14 + "_R8000", "d_ut_raw")
    x59 = dict(name="NASA X-59", mach=V(1.4, "x59_cruise", "published"),
               length_m=V(30.2, "x59_len", "published", note="99 ft"),
               altitude_used_km=V(16.8, "x59_cruise", "published", note="55,000 ft design cruise; the published 75 PLdB applies there"),
               unshielded=dict(
                   undertrack_pldb=V(75, ["x59_target_ntrs", "x59_carpet"], "published",
                                     note="design target (shaped 'thump'); flight measurements pending as of Sep 2026"),
                   carpet_max_pldb=V(75, ["x59_target_ntrs"], "published", note="design target; taken as the carpet maximum")),
               shielded=dict(
                   undertrack_pldb=V(round(75 + d_ut14), ["x59_target_ntrs", "r3d_csv"], "extrapolated",
                                     note=f"75 + simulated Concorde M1.4 delta {d_ut14:.1f} dB (R_m 8 km, He 30 m, W 1 km)"),
                   carpet_max_pldb=V(round(75 + d_mx14), ["x59_target_ntrs", "r3d_csv"], "extrapolated",
                                     note=f"75 + Concorde M1.4 carpet-max delta {d_mx14:.1f} dB")),
               delta_headline_db=dict(
                   undertrack=V(r1(d_ut14), "r3d_csv", "extrapolated", raw=r1(d_ut14_raw)),
                   carpet_max=V(r1(d_mx14), "r3d_csv", "extrapolated")),
               shielded_strip=V(dict(half_width_km_ge10dB=r1(f(c14, "strip_half_km")),
                                     computed_to_km=r1(f(c14, "y_max_km")), carpet_half_width_km=cutoff("Concorde", 1.4)),
                                "r3d_csv", "extrapolated", note="Concorde M1.4 carpet"),
               caveat="Our Whitham source model cannot represent a shaped low-boom signature. The shield's dB reduction for a shaped thump may differ (its leak is the >5 Hz edge-diffracted part); He vs H2 differ by <~2 dB; the X-59 at 20.25 km would already be somewhat quieter than at 16.8 km.",
               note="Already quiet by shaping; a shield would take it toward inaudible, but with low marginal value.")
    aircraft_list = [con, sr, x59, ov]

    # ---------------------------------------------------------------- refs
    boom_refs = [
        dict(label="Concorde boom (measured)", pldb=V(105, ["rathsam25", "x59_carpet"], "published"),
             nasa_analogy=V("comparable to nearby thunder or a car door slam heard from inside the car", "rathsam25", "published")),
        dict(label="X-59 design target ('sonic thump')", pldb=V(75, ["x59_target_ntrs", "x59_carpet"], "published"),
             nasa_analogy=V("comparable to distant thunder, or a car door slam from across the street; 'more like a car door closing ... across the street at your neighbor's house'",
                            ["rathsam25", "nasa_podcast"], "published"),
             secondary_analogy=V("about as loud as a car door slamming 6 m away", "phys_today", "published")),
        dict(label="F/A-18 dive manoeuvre simulating a quiet boom (NASA CarpetDIEM ground test)",
             pldb=V(67, "x59_carpet", "published", note="not an X-59 measurement")),
    ]
    everyday = [
        dict(label="Whisper / quiet room", dba=V(30, "niosh_meter", "published"), distance="~1 m"),
        dict(label="Normal conversation", dba=V("60 (60-70)", ["niosh_meter", "nidcd"], "published"), distance="~1 m"),
        dict(label="Light car traffic", dba=V(60, "spokane", "published"), distance="15 m (50 ft)"),
        dict(label="Freeway / highway traffic", dba=V(70, ["wsdot", "spokane"], "published"), distance="15 m (50 ft)"),
        dict(label="Diesel truck passing (~80 km/h)", dba=V("85-90", ["ne_dot", "caltrans_tens"], "published"), distance="15 m (50 ft)"),
        dict(label="Dog barking", dba=V("75-96", "dog_uk", "published"), distance="~1 m"),
        dict(label="Gas lawn mower / leaf blower", dba=V("80-100 (mower); 90 (leaf blower)", ["nidcd", "cdc_listen", "caltrans_tens"], "published"), distance="at the operator (~1 m)"),
        dict(label="Motorcycle", dba=V(95, ["cdc_listen", "nidcd"], "published"), distance="rider / nearby"),
        dict(label="Jet airliner fly-over", dba=V("105-110", ["ne_dot", "caltrans_tens"], "published"), distance="300 m (985-1,000 ft)"),
        dict(label="Rock concert", dba=V(110, ["cpuc_eir", "nidcd"], "published"), distance="15 m (50 ft) / front rows"),
        dict(label="Emergency siren", dba=V("110-129 (120 typical)", ["nidcd", "cdc_listen"], "published"), distance="standing nearby (distance not stated)"),
        dict(label="Jet takeoff", dba=V(120, "cpuc_eir", "published"), distance="60 m (200 ft)"),
    ]
    scale_note = dict(
        caveat="PLdB (Stevens perceived level of a short impulsive boom) and dBA (A-weighted level of a steady or slowly varying sound) are different scales; placing booms among dBA sounds is indicative only.",
        recommended_method=V("Put booms on the PL scale and use NASA's impulsive-sound analogies as the bridge (Doebler & Rathsam's 'acoustic thermometer' of familiar impulsive sounds rated in PL): ~105 PL = nearby thunder / car door slam from inside the car; ~75 PL = distant thunder / car door slam across the street. Show dBA everyday sounds on a separate, parallel axis, not as equal numbers. No published general PLdB<->dBA conversion for booms was found (published work correlates PL and A-weighted SEL but gives no fixed offset).",
                             ["x59_target_ntrs", "rathsam25"], "published"),
        db_rule=V("+10 dB is heard as roughly twice as loud; a 5 dB change is noticeable", "faa_basics", "published"),
    )

    shield_case = dict(
        description="Reference BoomShield: a 1 km-wide, 30 m-thick hydrogen layer ('floating-lid' pillow) floating with its top at 20 km; aircraft fly 250 m above it (20.25 km).",
        width_m=V(1000, ["r3d_sum", "struct"], "simulated"),
        gas=V("H2, 30 m thick", "struct", "simulated"),
        h_m=V(250, "r3d_sum", "simulated"),
        shield_top_km=V(20, "r3d_sum", "simulated"),
        architecture=V("floating lid: system mass on the top face (~2.4 kg/m2 lid), 0.1 kg/m2 floor film, soft ties K ~10 Pa/m, slack partitions",
                       "struct", "simulated"),
        acoustic_model_note="Acoustics simulated as a lift-closed H2 30 m slab (3-D carpets); the floating-lid detail was checked in 2-D BEM (-28.4 dB circle, -33.4 undertrack, re unshielded) and matches the slab.",
        caveats=[
            V("PLdB values use the conservative rise-time proxy; paired turbulence runs say the true reduction is ~5-8 dB larger at W 1 km.", "wt"),
            V("Near-track values use R_m 8 km and are still ~1.5 dB optimistic; R_m 3 km values are 2-7 dB optimistic.", "r3d_sum"),
            V("Along-track joints must be sealed: a 1 m open gap every 250 m costs ~9 dB undertrack (-24.5 vs -33.7).", "at"),
            V("Stiff ties or a stiff (load-bearing) skin can make the layer nearly transparent; faces must stay soft.", ["struct", "r3d_sum"]),
            V("Wind shifts and stretches the carpet by km to tens of km (deltas barely change); siting needs wind climatology.", "wt"),
            V("Concept-level 2-D/3-D linear+nonlinear model; no hardware, no flight test; H2 fire hazard is unresolved.", ["struct", "r3d_sum"]),
        ])

    cross_sections = [
        dict(name="Floating-lid pillow (H2 30 m)", verdict="baseline",
             reason="Best modelled acoustics with no static load on soft elements; mass on the lid makes the floor a stable H2-over-air interface.",
             numbers=V(dict(circle_dB=-28.4, undertrack_dB=-33.4, near_line_dB=-38.2, lift_kg_m2=2.46, face_load_kN_m=0.2, tie_K_Pa_m=10),
                       "struct", "simulated")),
        dict(name="Touching compliant-hoop tube raft (D ~60 m)", verdict="alternate",
             reason="Matches the 60 m slab with no ties, but the hoop must be both strong (3-9 kN/m) and soft (E*t <~1e5 N/m, 3-9 % strain): not standard balloon construction.",
             numbers=V(dict(circle_dB="-27.1/-27.9", undertrack_dB="-36.3/-37.8", dp_Pa="100-300", lift_kg_m2=3.90, net_kg_m2=2.58),
                       "struct", "simulated")),
        dict(name="Lens with cross-width struts", verdict="rejected",
             reason="Zero-pressure face tension ~200 kN/m needs 1 km CFRP compression struts: 9-130 kg/m2 of structure for 2-4 kg/m2 of lift.",
             numbers=V(dict(structure_kg_m2="9.1 (zero-pressure) / 70-130 (superpressure)", lift_kg_m2="1.95-3.9"), "struct", "simulated")),
        dict(name="Drop-stitch pillow (superpressure)", verdict="rejected",
             reason="The internal ties that hold it flat couple the faces 12-80x too stiffly (~1,200-2,400 Pa/m), short-circuiting the gas layer.",
             numbers=V(dict(face_coupling_Pa_m="1,210-2,420", smeared_K300_circle_dB=-18.3, vs_slab_dB=-26.0), "struct", "simulated")),
        dict(name="Vacuum / inverted balloon", verdict="rejected",
             reason="Nothing floats (structure 1.23-1.44x lift) and the rings/struts needed to hold the vacuum open are acoustic studs.",
             numbers=V(dict(structure_to_lift="1.23-1.44", raft_D30_dPLdB="-0.3 to -1.9 (M1.4-2.5)", lens_dPLdB="-2.6 to -6.7"),
                       "vac", "simulated")),
        dict(name="Raft with air gaps", verdict="rejected",
             reason="Air slits leak directly: a 2 m slit costs ~15 dB at the 3 km circle.",
             numbers=V(dict(gap_m=[2, 5, 15], circle_dB=[-12.3, -8.4, -4.2], touching_circle_dB=-27.0), "struct", "simulated")),
        dict(name="Narrow modules, 200-300 m wide", verdict="acoustically insufficient",
             reason="Edge diffraction dominates: W 250 m gives only -8 to -9 dB undertrack and -3 dB at the carpet maximum, and is 0-2 dB louder than no shield beyond ~18 km.",
             numbers=V(dict(W250_undertrack_dB="-8.3/-9.4", W250_max_dB="-3.1/-3.3", W500_undertrack_dB=-21, W500_max_dB=-8.4, W1000_undertrack_dB=-33.7),
                       ["r3d_sum", "handoff"], "simulated")),
    ]

    mooring = [
        dict(topic="Tether drag through the jet stream", text="A 5 cm tether crossing ~8 km of 60-80 m/s jet-stream wind carries ~350-630 kN of drag, more than the whole envelope side load; tether drag likely sets station-keeping.",
             numbers=V(dict(drag_kN="350-630", tether_cm=5, jet_m_s="60-80"), "addendum")),
        dict(topic="Tether drag, sized tethers", text="Sized UHMWPE tethers (6-10 cm) see ~0.6-1.0 MN of jet-stream drag each.",
             numbers=V(dict(drag_MN="0.6-1.0", diameter_cm="6-10"), "struct")),
        dict(topic="Module side load", text="Illustrative 250 m module (15 m thick, 30 m/s): ~600 N/m, ~150 kN per module.",
             numbers=V(dict(N_per_m=600, kN_per_module=150, U_m_s=30), "handoff")),
        dict(topic="Shield side force (1 km, 30 m)", text="0.50 / 1.12 / 1.99 kN per metre of corridor at 20 / 30 / 40 m/s; winds at 20 km average 10-20 m/s, winter 20-30, peaks ~50.",
             numbers=V(dict(U_m_s=[20, 30, 40], kN_per_m=[0.50, 1.12, 1.99]), "struct")),
        dict(topic="Tether loads and field", text="Tethers every 2-5 km carry 2.2-5.6 MN each at 30 m/s (90-220 t each); a 4,000 km corridor needs ~800-2,000 tethers through the jet stream.",
             numbers=V(dict(spacing_km="2-5", load_MN="2.2-5.6", mass_t="90-220", count="800-2,000"), "struct")),
        dict(topic="Propulsion instead of tethers", text="Excluded: ~67 kW per metre of corridor, ~270 GW for 4,000 km.",
             numbers=V(dict(kW_per_m=67, GW_4000km=270), "struct")),
        dict(topic="Daily lift swing (superheat)", text="Solar superheat of ~20 K raises lift ~10 % every day; constant-volume pressure swing 255 / 510 / 995 Pa for 10 / 20 / 39 K, so the layer must breathe 5-10 % daily (tethers, gas shuttle or coatings).",
             numbers=V(dict(lift_swing_pct=10, dp_Pa=[255, 510, 995], superheat_K=[10, 20, 39], tether_vertical_kN_m="0-2.4"), "struct")),
        dict(topic="Boom loading on the shield", text="Each pass hits the top face with 200-600 Pa (10-25x the lid's weight); 2e5-1e6 cycles over 30 years; flex-crack pinholes are the main leak driver.",
             numbers=V(dict(pulse_Pa="200-600", cycles="2e5-1e6"), "struct")),
        dict(topic="Hydrogen safety and supply", text="4,000 km x 1 km x 30 m holds 0.74 Mt H2 (0.8 % of a year's world supply; as helium, 45 years of output) and ~5 kt TNT-equivalent per km; at 55 mbar the layer is marginal for spark ignition; needs compartments every <=250 m and leak detection.",
             numbers=V(dict(H2_Mt=0.74, world_supply_pct=0.8, helium_years=45, kt_TNT_per_km=5, compartment_m=250, replenish_pct_yr="5-20"), "struct")),
        dict(topic="Wind moves the carpet", text="Jet-stream wind shifts or stretches the boom carpet by km to tens of km (cut-off 35-115 km) and a crosswind moves the shadow ~2.4 km downwind; the shield's dB reduction changes by <=0.4 dB undertrack.",
             numbers=V(dict(undertrack_delta_change_dB=0.4, crosswind_shift_km=2.4, cutoff_km="34.9-114.7"), "wt")),
    ]

    # headline
    ovl = ov_sim.get("delta_headline_db", {}) if isinstance(ov_sim, dict) else {}
    ovtxt = ""
    if ovl:
        ovtxt = " (Boom Overture, M1.7, simulated with an estimated preset: %s dB undertrack, %s dB at the carpet maximum)" % (
            ovl['undertrack']['value'], ovl['carpet_max']['value'])
    findings = [
        V("A 1 km-wide, 30 m hydrogen shield cuts a Concorde boom by about 34 dB under the track and about 25 dB at the loudest point of the carpet (conservative), taking it from ~98 PLdB (conservative proxy scale; ~108 raw) to ~64 PLdB under the track and at most ~73 PLdB anywhere in the computed carpet.", ["r3d_sum", "r3d_csv"]),
        V("Width, not gas or thickness, sets the result: 250 m wide gives only -8 to -9 dB and 500 m about -21 dB undertrack, because sound diffracts around the edges.", "r3d_sum"),
        V("The result holds from Mach 1.4 to 3.2: about -33 to -36 dB under the track and -21 to -29 dB at the loudest point for Concorde, Overture and the SR-71" + ovtxt + ".", ["r3d_csv", "ovr_run"]),
        V("These are conservative: with realistic turbulence the reduction is about 5-8 dB larger, and wind barely changes it (<=0.4 dB undertrack) though it moves where the carpet lands.", "wt"),
        V("Only a soft structure works: anything stiff that holds the faces apart (struts, stiff ties, drop-stitch, vacuum frames) or any open gap leaks the boom through; the floating-lid pillow is the baseline.", ["struct", "vac", "at"]),
        V("The hard problems are not acoustic: mooring a 1 km layer at 20 km against the jet stream (350-630 kN drag per tether, 800-2,000 tethers), a ~10 % daily lift swing, and 0.74 Mt of hydrogen next to air traffic.", ["addendum", "struct"]),
        V("Shaping (X-59, ~75 PLdB) and Boom's Mach-cutoff 'Boomless Cruise' (M1.1-1.3) work only at low Mach; a shield's niche is faster flight (Overture's Mach 1.7 and especially Mach 2+), where neither applies.", ["x59_target_ntrs", "boomless_faq", "addendum"]),
    ]

    data = dict(
        meta=dict(title="BoomShield infographic data", generated_by="boom-shield/infographic/make_data.py",
                  conventions=dict(
                      pldb="Absolute PLdB values are the rise-time 'proxy' (conservative headline); raw values are given alongside. Deltas are the less favourable of raw and proxy.",
                      altitude="Simulated aircraft fly at 20.25 km (shield top 20 km + h 250 m).",
                      R_m="Undertrack shielded values from R_m 8 km re-runs; carpet maxima and strips from R_m 3 km (all azimuths).",
                      tags="simulated = this repo's model; published = cited source; extrapolated = derived by us from published + simulated numbers"),
                  sources=SOURCES),
        shield_case=shield_case,
        aircraft=aircraft_list,
        boom_refs=boom_refs,
        everyday_refs=dict(items=everyday, **scale_note),
        cross_sections=cross_sections,
        mooring=mooring,
        headline_findings=findings,
    )
    with open(OUT, "w") as fh:
        json.dump(data, fh, indent=1, ensure_ascii=False)
    print("wrote", OUT)
    return data


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--run", action="store_true")
    ap.add_argument("--build", action="store_true")
    ap.add_argument("--procs", type=int, default=3)
    ap.add_argument("--cache", default=DEFAULT_CACHE)
    a = ap.parse_args()
    os.makedirs(a.cache, exist_ok=True)
    if a.run:
        run(a.procs, a.cache)
    elif a.build:          # re-summarise from the cache (no new runs except the cheap cruise points)
        summ = summarize_overture(a.cache)
        summ.update(cruise_points(a.cache))
        with open(os.path.join(a.cache, "overture_summary.json"), "w") as fh:
            json.dump(summ, fh, indent=1)
    if a.build:
        build(a.cache)


if __name__ == "__main__":
    main()
