# 3-D finite-width corridor study: ground carpets (v1)

From `run_3d.py` (about 8 h, 4 cores). Data: `3d_results.csv` (per φ, with y_g), `3d_carpet_summary.csv`. Figures: `3d_carpets.png`, `3d_signatures.png`.

## Bottom line

1. **Width sets the result, not opacity.** At each width, He 30, H₂ 30 and He 60 m agree to within about 2 dB. The leak that matters is edge diffraction above 5 Hz.
2. **W = 1 km** (Concorde M2, h 250 m): about **−33 undertrack** and **−25 at the carpet maximum** (outer carpet, y ≈ 39 km). The whole computed half-carpet is at least 10 dB quieter.
3. **W = 500 m:** −21 undertrack, −8 at the maximum, a 16–17 km protected strip. **W = 250 m:** −8 to −9 undertrack, −3 at the maximum. Beyond 18 km off-track the 250 m shield is **0–2 dB louder** than no shield.
4. Gate 0 gave −∞ for these stacks. Finite width costs 25–35 dB at 1 km.

## Setup and numerics

- Concorde M 2.0, h = 250 m, shield top at 20 km. Lift-closed faces of 1.15 kg/m² (He 30 m) and 1.24 kg/m² (H₂ 30 m). The lens uses its mean depth: 1.85 kg/m² faces, 6 m rim.
- Cross-plane BEM (lens: BEM to 25 Hz, PO above) with the aged incident field, then `propagate_from_circle` (age_from = aging_radius) at 50 kHz. Unshielded runs use the identical path.
- f_max = 100 Hz and df = 1/6 Hz everywhere. Circle signatures are cropped at 10⁻³ of peak; a 10⁻⁴ crop changes nothing.
- **R_m = 3 km** for all widths. Near-track azimuths (0, 4, 12, 24°) were re-run at **8 km**, the largest R_m inside the isothermal layer.
- φ = 0–6° in 1° steps, 2° to 30°, 3° to 54°: y_g = 0–39 km of the 47.7 km cutoff half-width, about 375 m per degree near track. Only φ ≥ 0 (symmetric). Edge projections: 26.6° (250 m), 45° (500 m), beyond cutoff (1 km).

## Headline (Concorde M2, h 250 m; Δ = less favourable of raw and proxy)

| stack | W | ΔPLdB undertrack, R_m 3 km / **8 km** | ΔPLdB max over carpet | protected half-width (≥10 dB) | < 5 Hz energy fraction at circle, φ = 0 |
|---|---|---|---|---|---|
| He 30 m | 250 m | −8.6 / **−8.3** | −3.1 | 0 (6.2 km in patches) | 0.53 |
| H₂ 30 m | 250 m | −9.7 / **−9.4** | −3.3 | 0 (6.2 km in patches) | 0.44 |
| He 30 m | 500 m | −22.5 / **−21.0** | −8.4 | 17.2 km | 0.90 |
| H₂ 30 m | 500 m | −22.6 / **−21.1** | −8.4 | 16.2 km | 0.83 |
| He 30 m | 1 km | −38.9 / **−33.9** | −25.5 | ≥ 39 km | 0.94 |
| H₂ 30 m | 1 km | −38.5 / **−33.7** | −25.0 | ≥ 39 km | 0.96 |
| He 60 m slab | 1 km | −40.6 / **−33.8** | −23.5 | ≥ 39 km | 0.96 |
| He 60 m lens | 1 km | −33.5 / **−31.8** | −22.7 | ≥ 39 km | 0.94 |
| He 30 m + 150 m fringe | 1 km | −∞ / **−∞** | −33.4 | ≥ 39 km | 0.95 |
| H₂ 30 m + 150 m fringe | 1 km | −∞ / **−∞** | −32.7 | ≥ 39 km | 0.98 |
| He 30 m, partitions every 100 m | 1 km | −38.9 / **−33.9** | −25.6 | ≥ 39 km | 0.95 |
| H₂ 30 m, 25 m line ties, K = 30 | 1 km | −38.6 / **−33.8** | −24.5 | ≥ 39 km | 0.96 |
| H₂ 30 m, 25 m line ties, K = 300 | 1 km | −0.2 / **−0.3** | +0.1 | 0 | 0.35 |

Unshielded: 108.4 / 97.7 PLdB (raw / proxy). The unshielded circle pipeline matches the direct near-field path to within 0.06 dB at every φ, and is R_m-independent.

## R_m sensitivity (W = 1 km He 30 m; ground PLdB raw/proxy)

| φ (y_g) | 2 km | 3 km | 4 km | 6 km | 8 km | 12 km* | 18 km* |
|---|---|---|---|---|---|---|---|
| 0 (0) | 57.0/54.9 | 61.4/58.8 | 63.8/60.9 | 66.1/62.8 | 67.2/63.8 | 68.2/64.7 | 68.9/65.2 |
| 12° (4.6 km) | 56.2/54.1 | 61.0/58.6 | 63.2/60.5 | 65.2/62.2 | 66.1/62.9 | 66.9/63.6 | 67.4/64.1 |
| 24° (9.7 km) | 68.7/65.0 | 71.0/67.0 | 71.6/67.7 | 72.1/68.2 | 72.4/68.4 | – | – |
| 54° (39 km) | 75.3/71.9 | 75.7/72.2 | 75.9/72.5 | – | – | – | – |

\*The circle extends below 11 km; these values are indicative only.

R_m = 3W is not enough near the track: the shadow keeps filling in with distance, but the ray start freezes the pattern at R_m. The 3 km near-track values are 2–7 dB optimistic. The 8 km values are still about 1.5 dB low.

## Edge treatments, cell walls, ties

- **150 m fringe:** undertrack is inaudible and the carpet maximum is 8 dB lower. **Treat this as qualitative only.** The levels lie below the BEM floor (about −52 dB), and the x-averaged model drops the cones scattered by real serrations.
- **Limp partitions every 100 m:** no effect (±0.3 dB). The leak is edge diffraction, not sloshing.
- **Line ties at K = 300 (25 m pitch) make the shield transparent.** The 1-D smeared K = 300 still gave −9.4. A circle check (H₂ 30 m, φ = 0) shows the discreteness effect:

| 5–10 Hz band, dB re unshielded | no ties | line 25 m | 12.5 m | 5 m | smeared |
|---|---|---|---|---|---|
| | −29.9 | −2.7 | −14.0 | −18.3 | −19.0 |

  Ties spaced wider than the gas's lateral decay length (about 20 m at 5 Hz) leak far more than their smeared equivalent. This is unexpected and needs an independent check.

## Other aircraft (W = 1 km)

| case | ΔPLdB undertrack, R_m 3 / **8 km** | max over carpet | protected strip |
|---|---|---|---|
| SR-71 M3.2, H₂ 30 m | −40.7 / **−35.4** | −20.8 | ≥ 42 of 54 km |
| Concorde M1.4, He 30 m | −37.1 / **−33.4** | −28.8 | ≥ 24.5 of 33 km |
| Concorde M2, h 100 m, He 30 m | −45.8 / **−40.8** | −34.6 | ≥ 39 km |

## Which band sets ground loudness

The < 5 Hz band holds 80–98 % of the circle energy at W ≥ 500 m: the transmitted, smoothed pulse (11–25 Pa peak, 150–250 ms rise). **Propagated alone it is inaudible in every W = 1 km case** and undertrack at 500 m. The > 5 Hz part alone reproduces the full PLdB to within 0.5–2 dB. That part is small edge-diffracted pulses, 3–8 Pa, arriving 0.1–1.1 s late. At W = 250 m the edge wave is a 25 Pa shock, and off-track both bands matter.

## Comparison with Gate 0 and the review

- **Gate 0:** −∞, against −33 / −25 here at 1 km. Ground peaks barely change (H₂ 30 m: 11.5 vs 11.8 Pa); only the audible band leaks.
- **Review floors (mapped in Gate 0):** −4 to −17 at 250 m and −13 to −23 at 1 km.
  - 250 m: inside that range undertrack, at the weak end for the carpet maximum.
  - 1 km: better than predicted. The leak is concentrated above 5 Hz in weak edge pulses, so the flat-leak mapping overstates it.
- **Fragility:**
  - The two edge waves of the symmetric strip add up on the track (an Arago-type spot). This makes y = 0 5–10 dB louder than φ = 1–2°.
  - Stiff 25 m-pitch ties erase the benefit.

## Not included

- Shield variation along x: transverse walls, joints, finite length.
- Wind, V/T drift and turbulence (only the rise-time proxy is used).
- Carpet-edge diffraction; φ beyond 54° was not computed.
- Ray start when R_m is not ≫ W. The near-track values remain R_m-limited.
- BEM floor of about −52 dB. Near-track high bands at 1 km sit at −36 to −44 dB re incident; the fringe cases are below the floor.
- Rigid walls, face tension, rim mass, flanking other than ties; asymmetric shields.
