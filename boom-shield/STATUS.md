# BoomShield — concept status

**Skyloop · 22 October 2026 · status after the first simulation round and the 0.11 psf check.**
Screening-level: acoustic numbers come from a 1-D/2-D model chain (linear shield physics, standard atmosphere, unshaped Concorde-class source) and are not design guarantees.

## 1. Bottom line

One acoustic mechanism holds up through the full source-to-ground chain, and one structural architecture, the floating-lid pillow, does not undo it. In **loudness** terms, a hydrogen layer **1 km wide and 30 m thick at 20 km**, with a Concorde-class aircraft at Mach 2 flying 250 m above, cuts about **−34 dB undertrack** and **−25 dB at the carpet maximum** (conservative scoring): roughly **58–73 PLdB** across the computed carpet against ~98 unshielded.

The FAA's proposed overland standard is not a loudness limit but a **peak overpressure of 0.11 psf (5.27 Pa)** (§4). It binds harder. The 30 m layer leaves an inaudible infrasonic swell of **0.25 psf at Mach 2**, which fails. Passing undertrack needs **H₂ ~55–80 m depending on Mach** (Overture-class only ~20 m, so lift sets its layer), and passing over the whole carpet needs **~1.5–2 km of width**. The thicker layer lifts more than the structure weighs, and the spare covers tether hold-down. The price is hydrogen: ~2–4 Mt instead of 0.74 Mt.

Caveats: 2-D cross-plane models; the psf results are a 1-D screen with 3-D width checks for Concorde M2 only; nothing structural, operational, safety-related or economic is closed. **The concept has passed an acoustic screen. It has not been shown to be buildable.**

## 2. Concept and evidence

BoomShield is a buoyant layer at about 20 km beneath a supersonic corridor: hydrogen between two thin membranes, with aircraft holding a track a few hundred metres above. The chat-derived surrogates ([README](README.md), [handoff](handoff/2026-10-20-boomshield-handoff.md)) are replaced by the `boomsim` package ([DESIGN](sim/DESIGN.md), [DESIGN-3D](sim/DESIGN-3D.md)): Whitham near field → nonlinear augmented-Burgers rays through US Standard Atmosphere 1976 → 1-D transfer matrix, 2-D cross-plane BEM and along-track Floquet orders at the shield → off-track rays to the ground carpet, with wind and turbulence → Mark VII PLdB, CSEL and peak psf, always against an unshielded run with identical numerics.

| Anchor | Model | Reference |
|---|---|---|
| Concorde, Mach 2, ground peak | 1.97 psf | 1.94 psf (NASA) |
| SR-71, Mach 3, ground peak | 0.99 psf | ~0.9 psf |
| Concorde carpet half-width (16.8 km) | 42.8 km | ~44 km |
| Mark VII vs PyLdB; BEM vs Whitham | ≤ 0.012 dB; 1–3.5 % | — |
| Secondary-boom delay | 8–13 min | 10–15 min |

## 3. Key physics findings

**3.1 The evanescent gas layer beats the mass law.** The boom pattern sweeps the shield at the aircraft speed V. Below the gas sound speed (H₂ to ~Mach 3.8, He to ~Mach 2.9) it cannot propagate inside the layer and tunnels through, decaying as exp(−κd), κ = ω√(1/V² − 1/c_g²). Verified independently. A limp mass-law sheet of the same face mass gives only −0.5 to −4.4 dB, because buoyancy and mass-law opacity both scale with air density.

**3.2 Width sets loudness; thickness and weight set peak psf.** In PLdB, He 30 m, H₂ 30 m and He 60 m agree within ~2 dB at each width: 250 m gives −8 / −3 dB (undertrack / carpet maximum), 500 m −21 / −8, 1 km −34 / −25 (SR-71 M3.2 −35 / −21; Concorde M1.4 −33 / −29). Audible loudness comes from 3–8 Pa edge-diffracted pulses above 5 Hz. The sub-5 Hz transmitted swell holds 80–98 % of the energy, is inaudible, and sets the **peak**, which only thickness, lower altitude and lighter aircraft reduce (§4).

**3.3 Symmetric edges focus on the centreline** (5–10 dB); asymmetric strips are worse off-track. A graded fringe looked ~8 dB better, but below the BEM floor.

**3.4 Ties, and anything that locks the volume, are acoustic studs.** For H₂ with tensioned faces at pitch ≤ 25 m, tie stiffness K ≤ 150 Pa/m at Mach 2 (40 at Mach 1.4); halve at 50 m pitch. The earlier "ties make it transparent" result was a mesh bug ([line_ties_check](sim/results/line_ties_check.md)). A skin with hoop stiffness E·t ≥ 3 × 10⁵ N/m is nearly transparent, so every viable architecture is quasi-zero-pressure and breathes 5–10 % per day. Derived for 30 m; recheck at the psf thickness.

**3.5 Re-steepening.** Marginal stacks re-shock on the descent (15 m He: −63 dB linear, −9.7 dB nonlinear). Score at the ground.

**3.6 Wind and turbulence do not change the conclusions.** Undertrack ΔPLdB moves ≤ 0.4 dB in realistic winds, but the carpet shifts kilometres to tens of kilometres. Paired turbulence runs make the headline 5.6–8.3 dB conservative. 251 automated tests.

## 4. The proposed FAA 0.11 psf standard

**What the NPRM says.** FAA, "Enabling Supersonic Overland Flight", 91 FR 40470 (2 July 2026), **docket FAA-2026-6935**. Proposed § 91.817(a)(2)(i) caps surface overpressure at **0.11 psf** for primary booms and for **direct and indirect secondary booms**. The FR text itself was not fetched (unreachable here); the rule and the § 1.1 definitions come from Holland & Knight (9 Jul 2026) and DLA Piper summaries. Treat the details as unverified.

**How our cases compare** (shield top 20 km, aircraft 250 m above; ground peak, psf):

| Aircraft | Unshielded | H₂ 30 m, 1 km | Min H₂ for 0.11 psf undertrack |
|---|---|---|---|
| Overture M1.7 (est. 62 t) | 1.12 | 0.065 (1-D): passes | 20 m |
| Concorde M1.4 ᶜ | 1.73 | 0.10 (1-D): marginal | 28 m |
| Concorde M2.0 | 1.77 | **0.24–0.25: fails** | **54 m** |
| Concorde M2.5 | 1.74 | 0.42 (1-D): fails | 77 m |
| SR-71 M3.2 ᵠ | 1.08 | 0.29: fails | 71 m |

ᶜ C_L > 0.5; ᵠ q > 35 kPa at this altitude. Helium needs about twice the thickness. 3-D undertrack values are 2–8 % below the 1-D screen. In PLdB every passing case, and 289 failing ones, is inaudible: **psf binds, not loudness.**

**Minimum H₂ thickness (m, 1-D undertrack, floating lid)** ([psf_sweep](sim/results/psf_sweep.md)):

| H_s | Overture M1.7 | Concorde M2.0 | Concorde M2.5 | SR-71 M3.2 |
|---|---|---|---|---|
| 18 km | 18 | 49 | 70 | 64ᵠ |
| 20 km | 20 | **54** | 77 | 71ᵠ |
| 22 km | 23 | 60 | 86 | 79 |

**Width.** In 3-D (Concorde M2, H₂ 60 m), **1 km wide fails off-track: 11.7 Pa at y = 39 km**, from edge diffraction of the ~1 Hz swell on rays that cross only 150–250 m inside the edge. **1.5 and 2 km pass everywhere (≤ 4.2 Pa)**. Only this aircraft has been run; plan on ~1.5–2 km.

**Weight × thickness budget.** There is no hard floor: the swell carries the aircraft's lift, and the peak scales as **W · d^−(1.2…1.5)**, rising with altitude through the ray gain R(H). One fitted constant matches 699 cases within ±1.8 dB. Maximum aircraft weight at 20 km, H₂ floating lid: **66 t per 30 m at Mach 2** (174 t at 60 m, 346 t at 100 m); 104 t at Mach 1.7; 39 t at Mach 2.5. Weight × Mach is the market limit. Spare lift belongs in the lid: a fixed 2.4 kg/m² lid needs 61 m at Mach 2.

**Secondary booms** ([secondary_booms](sim/results/secondary_booms.md)). The shield does not make them worse: it moves the ground reflection that returns as an *indirect* secondary up to 20 km, and shielded worst bins are **0.54–1.02×** unshielded in seven cases. But unshaped Concorde- and SR-71-class aircraft at 20 km exceed 0.11 psf in secondary booms whenever winds form a stratospheric duct (**0.14–0.37 psf**, shield or not); calm air and headwinds give ≤ 0.02 psf. The model reads ~×2 high against Concorde observations, and winds change levels 15–40×. This is an aircraft and routing problem; an upturned shield edge might help (needs a BEM check).

**Altitude** ([altitude-trade](altitude/altitude-trade.md) §7). Lower helps psf, but weakly (~10 % less thickness per 2 km), while hydrogen per km, weather, tethers and airspace worsen. **18–19 km does not now look better; 19–22 km stands, 20 km for Mach 2.** At 20 km, 54 m lifts ~2 kg/m² more than the 2.46 kg/m² structure: 0.55 kg/m² closes the tether hold-down shortfall and the rest stays on the lid (~+4 % in peak, ~+2 m of H₂, estimated).

**Docket point.** A peak rule penalises an inaudible 0.4–2.7 s swell; building rattle is the real issue. Worth a comment.

## 5. What is ruled out

| Option | Why it fails | Source |
|---|---|---|
| Narrow 200–300 m modules | −8 to −9 dB undertrack; louder than no shield beyond 18 km | [3d_summary](sim/results/3d_summary.md) |
| 1 km width under the psf limit | Edge diffraction: 11.7 Pa off-track even with 60 m H₂ | [psf_sweep](sim/results/psf_sweep.md) |
| Mass-law-only shields | −0.5 to −4.4 dB; buoyancy cannot buy mass | [gate0](sim/results/gate0_summary.md) |
| Helium | ~2× the thickness for psf; 30 m × 1 km alone is ~45 years of world output; Mach ceiling ~2.9 | [structure §4](structure/structure-concept.md) |
| Vacuum / inverted balloons | Structure 1.23–1.44× lift; spacers are studs | [vacuum review](vacuum/vacuum-balloon-review.md) |
| Drop-stitch, struts, stiff skins, gapped rafts | Studs or leaks: ties 12–80× too stiff, E·t ≥ 3 × 10⁵ N/m, 2 m slit costs 15 dB | [structure](structure/structure-concept.md) |
| Open joints, coarse serrations | 1 m gap per 250 m: −24.5 dB instead of −33.7; Λ 10–30 m up to 6 dB worse | [alongtrack](sim/results/alongtrack_summary.md) |

## 6. Baseline: the floating-lid pillow, re-sized for 0.11 psf

| Parameter | Value |
|---|---|
| Layer | H₂, **~55 m at Mach 2; ~70–80 m at Mach 2.5–3.2; ~30–36 m for Mach 1.7** (set by lift, not psf) |
| Width / altitude | **~1.5–2 km**; top at **19–22 km** (20 km for Mach 2); aircraft h ≈ 250 m above |
| Lift at 20 km, 55 m | ≈ 4.5 kg/m²: 2.46 structure + 0.55 tether hold-down + ~1.5 extra lid mass |
| Lid / floor | Lid carries the mass (~3.8 kg/m²); floor film ~0.1 kg/m² |
| Ties / partitions | Soft, K ≈ 10 Pa/m (3–30; > 0.8 against lid Rayleigh–Taylor); slack bellows ~100 m; H₂ compartments ≤ 250 m |
| Pressure | Quasi-zero; breathing stroke ~5 m/day at 55 m (scaled) |
| Hydrogen | ~2.0 Mt (Mach 2, 1.5 km) to ~3.8 Mt (Mach 2.5, 2 km) per 4,000 km; ~14–27 kt TNT per km |
| Acoustics | Concorde M2, 60 m × 1.5 km: 4.2 Pa (0.088 psf) carpet maximum. PLdB −33 / −25 dB at 30 m × 1 km |

**Why thicker helps the structure.** At 20 km the old 30 m layer lifted exactly the 2.46 kg/m² system and had nothing left to pretension the tethers (+0.55 kg/m² needed). The psf layer lifts ~2 kg/m² spare, which closes that gap. Diverting the whole surplus to tethers instead would cost ~16 % in peak, so keep the rest on the lid.

**Requirements:** ties within §3.4; absorb the ~10 % daily lift swing without stiffening the volume; membrane-closed joints (≲ 1 m per 250 m); graded edges; limp faces. Alternate: a raft of compliant-hoop tubes (D ≈ 60 m).

## 7. Corrections to the original handoff

- **Width:** "higher Mach → narrower shield" is wrong; W ≳ 4h for PLdB, ~1.5–2 km for psf.
- **Buoyancy:** the evanescent layer is the lift body (≥ 75 % of displaced-air lift).
- **X-59 scaling** does not predict loudness; unshaped sources need 20 dB+, not ~8.
- **Boom load ≫ dead load:** 200–600 Pa per pass, 10⁵–10⁶ cycles.
- **Tethers:** 17 MN each at 20 km (~26 MN for the psf layer), not 150 kN; propulsion is excluded (~67 kW/m).
- **Skyloop's inverted balloon is not "entirely in tension"** ([07-sections](../docs/07-sections.md)): an R = 40 m ring carries 11 MN of compression.
- **The README's "Main result"** is superseded.

## 8. Open risks and next steps (ranked)

| # | Item | Why it matters | Next step |
|---|---|---|---|
| 1 | H₂ safety near traffic | 2–4 Mt, ~14–27 kt TNT per km; leaked H₂ rises toward flight level | Hazard analysis; compartments; leak detection that closes sections |
| 2 | Regulatory fit | 0.11 psf penalises an inaudible swell; secondary booms break it without a shield | Docket comment with a building-response metric; route on wind climatology |
| 3 | Fatigue, 10⁵–10⁶ cycles | Faces go slack every pass; pinholes drive 5–20 %/yr top-up | 217 K coupon tests; slack-lid model under 200–600 Pa |
| 4 | Daily lift swing and breathing | ~10 % swing; stroke grows with thickness | Tethers, gas shuttle or coatings, without stiffening |
| 5 | Aerodynamics of a 1.5–2 km surface | Gusts, mountain waves, flutter, lid Rayleigh–Taylor | Coupled aero-structural model |
| 6 | Tether field | ~1,300 tethers through the jet stream at ~26 MN; Class A obstruction | Layout and station-keeping study |
| 7 | Economics | Inventory and area up 2–5×; no cost model | Cost per movement (handoff §9) |
| 8 | psf in 3-D beyond Concorde M2 | Width and thickness rest on one aircraft | 3-D width runs at psf thickness for each aircraft (~10 min each at 30 Hz) |
| 9 | Ties and lid at the new thickness | Tie limits set for 30 m; lid mass affects peak | Floating-lid 3-D carpet at ~55 m × 1.5 km with soft ties |
| 10 | Edges | Graded edges (PLdB upside); upturned edges (secondary booms) | Graded-sheet model; BEM edge-leak check |
| 11 | Shaped source; numerical floors | N-waves only; near-track ~1.5 dB optimistic; BEM floor ~−52 dB | Shaped near field; stratified start; better solver |

## 9. Evidence map

| Claim | Supporting file |
|---|---|
| Original concept and surrogates | [README](README.md), [handoff](handoff/2026-10-20-boomshield-handoff.md) |
| Mass law vs buoyancy; tether drag; boom load | [review addendum](handoff/2026-10-20-review-addendum.md) |
| Validation anchors | `sim/results/validation_*.csv`, `sim/tests/` |
| Gate 0; mechanism review; BEM validation | [gate0](sim/results/gate0_summary.md), [review_evanescent](sim/results/review_evanescent.md), [crossplane](sim/results/crossplane_validation.md) |
| Width sets loudness; −34 / −25 dB | [3d_summary](sim/results/3d_summary.md) |
| Tie rule; mesh-bug correction | [line_ties_check](sim/results/line_ties_check.md) |
| Joints, serrations; wind and turbulence | [alongtrack](sim/results/alongtrack_summary.md), [wind_turbulence](sim/results/wind_turbulence.md) |
| 0.11 psf: min thickness, width, weight budget | [psf_sweep](sim/results/psf_sweep.md), `psf_sweep.csv` |
| Secondary booms; NPRM definitions | [secondary_booms](sim/results/secondary_booms.md) |
| 19–22 km; spare lift closes hold-down | [altitude-trade](altitude/altitude-trade.md) |
| Vacuum branch closed | [vacuum-balloon-review](vacuum/vacuum-balloon-review.md) |
| Floating-lid baseline, locked-volume rule, H₂ safety | [structure-concept](structure/structure-concept.md) |
