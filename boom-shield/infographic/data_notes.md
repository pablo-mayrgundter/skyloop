# BoomShield infographic data: derivation notes

`data.json` is built by `make_data.py --build`. Every number has a `src`, which points into `meta.sources`, and a `tag` (simulated, published or extrapolated).

## Conventions

- Aircraft fly at 20.25 km: the shield top at 20 km plus h = 250 m.
- **Absolute PLdB** uses the conservative rise-time proxy; raw values are given alongside.
- **Deltas** take the less favourable of raw and proxy. Paired turbulence runs give 5–8 dB more reduction (`wind_turbulence.md`).
- **Shielded values:**
  - The undertrack value comes from R_m = 8 km.
  - Carpet maxima and the ≥10 dB strip come from R_m = 3 km.
  - If a maximum lies at the outermost computed azimuth, it is flagged: the carpet beyond it was not computed.
- **Shield case:** a 1 km, H₂ 30 m slab, standing in for the floating-lid pillow.

## Aircraft

- **Concorde and SR-71** values are read from `3d_carpet_summary.csv`.
  - At 16.8 km the model gives 1.97 psf, against 1.94 published.
  - Measured ~105 PLdB lies between the model's raw (109.7) and proxy (98.9) values.
- **Overture** uses my own preset; `source.py` is unchanged.
  - Inputs: 61.3 m, M 1.7, 62 t mid-cruise (0.81 × 77 t, the ratio the Concorde preset uses), 300 m³.
  - **The 77 t MTOW predates the 2022 redesign, so absolute levels are estimates.** Deltas depend less on weight than absolute levels do.
  - It was run through the `run_3d.py` functions with identical numerics: 26 azimuths at R_m 3 km, plus 0/4/12/24° at 8 km. The cache is outside `3d_cache`.
  - The 18.3 km cruise point is indicative only: the modelled front there is not yet steep (13 ms rise).
- **X-59:** the source model cannot represent shaped booms.
  - Unshielded: NASA's 75 PL target.
  - Shielded: 75 plus the Concorde M1.4 deltas (He 30 m, which agrees with H₂ to within ~2 dB): −33.4 dB undertrack and −28.8 dB at the maximum. That gives ~42 and ~46 PL.

## References and scales

- **All web figures are marked unverified.** Page fetches were blocked here, so they come from search-result summaries.
- **Re-check before publishing:** the Overture MTOW, siren and dog-bark levels, and the Caltrans values, which were read from local reproductions.
- **The dBA list** mixes CDC/NIOSH, NIDCD, FAA and state DOT figures.
- **PL vs dBA:** no published fixed conversion for booms was found. Recommended approach:
  - Place booms with NASA's impulsive analogies: 105 PL ≈ nearby thunder, or a car door slam heard from inside the car; 75 PL ≈ distant thunder, or a car door across the street.
  - Show dBA on a parallel axis.

## Repo studies

`cross_sections` comes from `structure-concept.md` (2-D BEM), `vacuum-balloon-review.md` and `3d_summary.md`. `mooring` comes from the review addendum, the handoff and `structure-concept.md` §3–4.
