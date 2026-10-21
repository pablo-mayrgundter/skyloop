# BoomShield — concept status

**Skyloop · 21 October 2026 · status after the first simulation round.**
Screening-level: acoustic numbers come from a 1-D/2-D model chain (linear shield physics, standard atmosphere, unshaped Concorde-class source) and are not design guarantees.

## 1. Bottom line

The concept now has one acoustic mechanism that holds up through the full source-to-ground chain, and one structural architecture that does not undo it. The configuration is a hydrogen layer about **1 km wide and 30 m thick at 20 km altitude**, with a Concorde-class aircraft at Mach 2 flying 250 m above it. It cuts ground loudness by about **−34 dB undertrack** and about **−25 dB at the loudest point of the carpet**. Both figures use conservative scoring: the less favourable of raw Mark VII PLdB and a rise-time proxy. On the proxy scale the whole computed carpet comes out at roughly **58–73 PLdB**, against about **98 PLdB unshielded** (108 raw). That is at or below the **~75 PLdB** X-59 targets for its own shaped boom.

The benefit comes from width, not opacity. A 500 m shield gives −21 dB undertrack and −8 dB at the carpet maximum. A 250 m shield gives −8 and −3. The main caveats:

- The model is 2-D in the cross-plane and invariant along track. Along-track effects are added only by physical optics.
- Near-track values are still about 1.5 dB optimistic because of where the rays start.
- Faces are assumed lift-closed, which is optimistic.
- Nothing structural, operational, safety-related or economic has been closed.

The concept has passed an acoustic screen. It has not been shown to be buildable.

## 2. What the concept is, and how the evidence was built

BoomShield is a buoyant layer flown at about 20 km beneath a supersonic corridor, so the downward boom has to cross it before reaching the ground. In its current form the layer is hydrogen between two thin membranes, and aircraft hold a fixed track a few hundred metres above it. The earlier chat-derived surrogates ([README](README.md), [handoff](handoff/2026-10-20-boomshield-handoff.md)) have been replaced by a simulation package, `boomsim` ([sim/DESIGN.md](sim/DESIGN.md), [sim/DESIGN-3D.md](sim/DESIGN-3D.md)).

**Simulation chain:**

1. **Source:** Whitham F-function (Sears–Haack volume plus lift) near-field signature.
2. **Propagation:** nonlinear augmented-Burgers along rays through US Standard Atmosphere 1976, with absorption; ×1.9 ground reflection.
3. **Shield:** 1-D layered transfer matrix (Gate 0); then a 2-D cross-plane boundary-element model (BEM) with gas cells, membranes and ties, physical optics above a split frequency; then along-track Floquet orders for serrations and joints.
4. **Off-track rays:** signatures on a matching circle of radius R_m (3–8 km) launched down off-track rays to the ground carpet; wind via moving-medium rays, turbulence via a phase screen.
5. **Scoring:** Stevens Mark VII PLdB plus CSEL, always against an unshielded run through identical numerics.

**Validation anchors**

| Anchor | Model | Reference |
|---|---|---|
| Concorde, Mach 2, ground peak | 1.97 psf | 1.94 psf (NASA) |
| SR-71, Mach 3, ground peak | 0.99 psf | ~0.9 psf |
| Concorde carpet half-width (16.8 km) | 42.8 km | ~44 km (1 mile per 1,000 ft rule) |
| Mark VII PLdB vs PyLdB | ≤ 0.012 dB on N-waves; 0.003 dB on PyLdB's shaped reference signature | — |
| Cross-plane incident field vs Whitham | 1–3.5 % max error | — |
| BEM edge-diffraction floor vs analytic formula | within ~2 dB | — |
| Unshielded circle pipeline vs direct path | ≤ 0.06 dB | — |
| Discrete line ties: BEM vs independent Floquet model | ~1 dB | — |
| Automated tests | 225 (pytest) | — |

## 3. Key physics findings, in order of importance

**3.1 The evanescent "subsonic gas layer" beats the mass law.** On any horizontal plane the boom pattern moves at the aircraft speed V. If V is below the gas's sound speed, the boom cannot travel as a wave inside the layer. It can only tunnel through, decaying as exp(−κd) with κ = ω√(1/V² − 1/c_g²), and the impedance mismatch at the faces adds further loss.

At 216.65 K this holds for helium below about Mach 2.9 (c = 866 m/s) and hydrogen below about Mach 3.8 (c ≈ 1,120 m/s). It is the known boom-into-ocean problem; an adversarial review verified it independently to 0.01 dB.

A limp membrane relying on the mass law gives only −0.5 to −4.4 dB for the same face mass: buoyancy and mass-law opacity both scale with air density, so lift cannot buy opacity. At normal incidence the gas stack is an ordinary double wall; the gain comes from the trace-speed condition.

In 1-D, layers of 30 m or more were inaudible (ΔPLdB = −∞), so real loudness is set by leaks. Purity matters: 5 % air drops c_He to 751 m/s and halves the benefit of 15 m of helium.

**3.2 Width sets ground loudness, not opacity.** At each width, He 30 m, H₂ 30 m and He 60 m agree to within about 2 dB.

| Width (Concorde M2, h 250 m) | ΔPLdB undertrack (R_m 8 km) | ΔPLdB at carpet maximum | Half-width at least 10 dB quieter |
|---|---|---|---|
| 250 m | −8 to −9 | −3 (0–2 dB *louder* beyond 18 km) | 0 (patches only) |
| 500 m | −21 | −8 | 16–17 km |
| 1 km | −34 | −25 | ≥ 39 km (whole computed carpet) |
| 1 km, SR-71 M3.2 (H₂) / Concorde M1.4 / h 100 m | −35 / −33 / −41 | −21 / −29 / −35 | ≥ 42 / ≥ 24.5 / ≥ 39 km |

The sub-5 Hz transmitted pulse holds 80–98 % of the energy but is inaudible; PLdB is set by 3–8 Pa edge-diffracted pulses above 5 Hz. Finite width costs 25–35 dB relative to 1-D. Rule of thumb: W ≳ 4h.

**3.3 Symmetric edges focus on the centreline.** The two edge waves of a symmetric strip add on the track, an Arago-type spot. This makes y = 0 5–10 dB louder than φ = 1–2°.

Making the strip asymmetric is not a fix. It barely changes the track and raises the worst off-track azimuth by 4–9 dB in the 5–100 Hz bands. A 150 m graded-opacity fringe made the undertrack inaudible and lowered the carpet maximum by 8 dB in the x-averaged model. That result is qualitative only, because it lies below the BEM floor.

**3.4 Ties, and anything that locks the volume, are acoustic studs.** Anything coupling the faces more stiffly than the gas bypasses the mechanism, like a stud in a partition wall. Design rules:

- **Ties.** With tensioned faces (N ≥ 10⁴ N/m) at pitch P ≤ 25 m, tie stiffness must satisfy K ≤ 150 Pa/m (H₂, Mach 2) or 40 Pa/m (H₂, Mach 1.4). Halve these limits at 50 m pitch.
- **Tendons.** Tendons spanning the 30 m layer at 1 m along-track spacing and 25 m pitch need EA ≤ 1.1 × 10⁵ N each.
- **The earlier "ties make the shield transparent" result (−0.3 dB) was a mesh bug.** Each tie acted on one boundary element whose width scaled with frequency, which put a spurious resonance in the 5–10 Hz band. With a realistic tie width of 0.5 m, discrete ties leak *less* than smeared ties. An independent Floquet model agrees to about 1 dB. The smeared model is about right or slightly pessimistic.
- **Locked volume.** The structure study extends the rule from members between the faces to anything that holds the volume fixed. A skin with hoop-stretch stiffness E·t ≥ 3 × 10⁵ N/m makes the shield nearly transparent, so E·t must stay ≲ 10⁵ N/m. Superpressure up to 1 kPa is harmless in itself.
- **Consequence.** Every viable architecture is quasi-zero-pressure and breathes 5–10 % per day. Partitions must be slack; limp partitions every 100 m changed the result by ±0.3 dB.

**3.5 Re-steepening.** Marginal stacks re-shock during the 20 km descent: 15 m of helium at Mach 2 gives −63 dB with a linear lower leg but only −9.7 dB with nonlinear propagation. Opaque stacks do not; what they transmit is too weak and slow (11–25 Pa, 150–250 ms rise). The scattered high-frequency leak gains +1.9 to +3.4 dB from nonlinearity. Shield-plane peak ratios understate ground peaks, so scoring must be at the ground, in PLdB.

**3.6 Wind and turbulence do not change the conclusions.**

- **Wind.** In realistic jet-stream and stratospheric-wind profiles, the undertrack ΔPLdB moves ≤ 0.4 dB and the carpet maximum moves −2.8 to +1.9 dB. The carpet and its shadow shift by kilometres to tens of kilometres, so siting needs wind climatology.
- **Turbulence.** The median level stays within about 1 dB, with a P10–P90 spread of about ±3 dB. The paired ΔPLdB tracks the raw value, so the headline numbers are conservative by 5.6–8.3 dB at 1 km.

## 4. What is ruled out, and why

| Option | Why it fails | Source |
|---|---|---|
| Narrow 200–300 m modules as the acoustic surface | −8 to −9 dB undertrack, −3 dB at carpet maximum, louder than no shield beyond 18 km | [3d_summary](sim/results/3d_summary.md) |
| Mass-law-only (limp, non-evanescent) shields | −0.5 to −4.4 dB; buoyancy cannot buy mass | [gate0](sim/results/gate0_summary.md), [addendum §1](handoff/2026-10-20-review-addendum.md) |
| Helium | 4,000 km × 1 km × 30 m is about 45 years of world helium output; Mach ceiling ~2.9; purity-sensitive | [structure §4](structure/structure-concept.md) |
| Vacuum / inverted balloons | Structure weighs 1.23–1.44× the lift at any size; the members holding the faces apart are studs, so buildable sizes fall back to the mass law | [vacuum review](vacuum/vacuum-balloon-review.md) |
| Drop-stitch at superpressure | Smeared ties of 1,210–2,420 Pa/m, 12–80× too stiff | [structure §2](structure/structure-concept.md) |
| Struts or lenses spanning the width | 9–130 kg/m² of structure for 2–4 kg/m² of lift; strutted lens acoustics −2.6 to −6.7 dB | [structure §2](structure/structure-concept.md), [vacuum review](vacuum/vacuum-balloon-review.md) |
| Stiff-skin superpressure rafts; rafts with air gaps | E·t ≥ 3 × 10⁵ N/m is nearly transparent; a 2 m air slit costs 15 dB | [structure §1](structure/structure-concept.md) |
| Open module joints | A 1 m gap every 250 m: −24.5 dB instead of −33.7; ~1 dB loss needs a gap ≲ 0.1 m | [alongtrack](sim/results/alongtrack_summary.md) |
| Coarse edge serrations (Λ ≈ 10–300 m) | Act as a grating: Λ 10–30 m is up to 6 dB worse than a sharp edge at 23–35 km off-track; full benefit only for Λ ≲ 1 m | [alongtrack](sim/results/alongtrack_summary.md) |

## 5. Current baseline: the floating-lid pillow

The structure study's recommended baseline ([structure-concept](structure/structure-concept.md), architecture d):

| Parameter | Value |
|---|---|
| Layer | H₂, 30 m thick, ~1 km wide, top at ~20 km; aircraft h ≈ 250 m above |
| Lift | ≈ 2.46 kg/m² (H₂ lift density 0.0827 kg/m³) |
| Top face ("lid") | Carries the system mass, ≈ lift (2.4 kg/m² modelled; 1.2–2.4 to test) |
| Floor | Light film, ~0.1 kg/m²; a gravitationally stable H₂-over-air interface |
| Ties | Soft, K ≈ 10 Pa/m (3–30); K > 0.8 Pa/m needed against lid Rayleigh–Taylor instability |
| Partitions | Slack bellows (~100 m); H₂ compartments ≤ 250 m |
| Face tension | ~0.2 kN/m; no compression member, no static tie load |
| Pressure | Quasi-zero; layer thickens ~2.7 m per day |
| Acoustics (5–30 Hz, 3 km circle / undertrack / 100 m below) | −28.4 / −33.4 / −38.2 dB re unshielded (H₂ slab reference: −26.0 / −32.0 / −36.6) |

**Requirements it must meet:**

- **Tie stiffness.** Within the limits in §3.4. Daily breathing should push ≤ 27 Pa through the ties at K = 10, or ≤ 8 Pa at K = 3.
- **Breathing volume.** The ~10 % daily lift swing must be absorbed (tethers, a solar gas shuttle into rim tubes, or low-superheat coatings) without stiffening the volume.
- **Joints.** Membrane-closed. The gap must be ≲ 1 m at 250 m module length or ≲ 2 m at 500 m; open gaps must be ≲ 0.1 m.
- **Edges.** No coarse serrations. Use a graded-property edge, or periodic Λ ≲ 1 m.
- **Faces.** Limp, with in-plane wave speed well away from V (avoiding coincidence). Spare lift goes into face mass, since light faces lose 3–9 dB.

The alternate is a touching raft of compliant-hoop tubes, D ≈ 60 m. It matches the 60 m slab acoustically, but it needs a hoop that is both strong and soft, which is not standard balloon construction.

## 6. Corrections to the original handoff

- **Mach-angle width formula.** The handoff already retired W ≈ 2h/√(M²−1). The addendum showed that in uniform air a ray at azimuth φ lands at y = h·tan φ whatever the Mach number, and that a real atmosphere widens the carpet with Mach. "Higher Mach → narrower shield" is wrong. The simulations give W ≳ 4h.
- **Buoyancy closure.** A 15 m helium envelope gives 0.90 kg/m² at 20 km against a 4 kg/m² target. The addendum added that a limp mass-law shield cannot close lift against acoustics (it needs T ≳ 0.47·r·λ); the evanescent mechanism escapes this. And any evanescent gas has density ≤ γp/V² = 0.022 kg/m³, so it lifts ≥ 75 % of the displaced air: the acoustic layer is the lift body and cannot be hung as a neutral sheet.
- **X-59 scaling.** Matching W/h ratios and scaling X-59 amplitudes does not predict loudness. The opacity cases assumed shaped sources. Unshaped sources model at 103 PLdB raw for SR-71 at Mach 3.2 and 108 raw for Concorde at Mach 2. Needed reductions are therefore 20 dB or more, not ~8. The ~75 PLdB X-59 figure is a reference point, not an adopted acceptance threshold.
- **Boom load ≫ dead load.** The top face sees 200–600 Pa every pass, 10–25× the lid's lift head, and goes slack each time: 2 × 10⁵ to 10⁶ cycles over 30 years at 20–100 flights a day. Fatigue, not catenary dead load, sizes the faces.
- **Tether drag.** The drag on a tether crossing the jet stream (350–630 kN; 0.6–1.0 MN scaled to the baseline) exceeds the handoff's 150 kN module side load. Tethers every 2–5 km would carry 2.2–5.6 MN each. Station-keeping by propulsion is excluded: it would need about 67 kW per metre of corridor.
- **The Skyloop inverted balloon is not "entirely in tension".** [`../docs/07-sections.md`](../docs/07-sections.md) claims so, with ~1 t rings every 50 m. The vacuum review finds an R = 40 m ring carries 11 MN of hoop compression and sizes to about 10.5 t. This affects the main Skyloop docs, not only BoomShield.
- **The README's "Main result"** comes from the earlier surrogate; this note supersedes it.

## 7. Open risks and next steps (ranked)

Items 1–6 could close the concept. Items 7–12 refine the acoustic numbers.

| # | Item | Why it matters | Next step |
|---|---|---|---|
| 1 | H₂ safety near traffic | ~5 kt TNT-equivalent of stored energy per km; at 55 mbar the layer is marginal for spark ignition and unsafe near engine exhaust; leaked H₂ rises toward flight level | Hazard analysis before hardware; compartments, leak detection that closes sections to traffic |
| 2 | Fatigue under 10⁵–10⁶ boom cycles | Faces go slack on every pass; flex-crack pinholes drive leakage (5–20 %/yr replenishment) | Coupon tests at 217 K; nonlinear slack-lid model under a 200–600 Pa pulse; tie damping |
| 3 | Daily lift swing (~10 %) | Zero-pressure lift does not self-limit with altitude; stiffening the volume kills the acoustics | Compare tethers, solar gas shuttle, low-superheat coatings |
| 4 | Aerodynamics of a 1 km surface | Gusts, mountain waves, roll, flutter of slack faces; lid Rayleigh–Taylor stability | Coupled aero-structural model |
| 5 | Tether field | 800–2,000 tethers through the jet stream, each carrying MN-class loads; clearance for traffic | Station-keeping and tether-layout study |
| 6 | Economics | No cost model; at 1 km width, every $10/m² is $10 M/km; first H₂ fill $1–4 B | Cost per useful movement, per the handoff's §9 |
| 7 | Module joints | Open gaps leak under the whole carpet; membrane joints work to 1–2 m | Joint detail, end domes, bellows partitions in the 3-D model |
| 8 | Graded-property edges | Largest acoustic upside (~8 dB at carpet maximum), but unverified | Model a sheet graded in y; needs a lower BEM floor |
| 9 | Shaped-aircraft (X-59-type) source | All results use unshaped N-waves; the leak budget for a shaped signature is unknown | Add a shaped near-field source |
| 10 | R_m convergence near track | 3 km values are 2–7 dB optimistic; 8 km values are still ~1.5 dB low | Start rays from a stratified near field |
| 11 | BEM floor (~−52 dB) | Hides the fringe cases and the high-band near-track leak | Raise numerical precision or switch solvers |
| 12 | Reflected cone, wind-shifted shadow | \|R\| ≈ 1 sends the boom back up, reaching flight level 200–700 m behind the nose; the carpet moves km with wind | In-trail traffic rules; siting climatology |

## 8. Evidence map

| Claim | Supporting file |
|---|---|
| Original concept, surrogate results, Mach-angle width rule | [README](README.md), [handoff](handoff/2026-10-20-boomshield-handoff.md) |
| Mass law couples to buoyancy; width independent of Mach; tether drag; boom load; helium inventory | [review addendum](handoff/2026-10-20-review-addendum.md), [quick_checks.py](handoff/quick_checks.py) |
| Model contracts | [sim/DESIGN.md](sim/DESIGN.md), [sim/DESIGN-3D.md](sim/DESIGN-3D.md) |
| Concorde 1.97 psf, carpet half-width 42.8 km | `sim/results/validation_propagation.csv`, `sim/results/validation_carpet.csv` |
| PyLdB match | `sim/tests/test_loudness.py`, [THIRD_PARTY_NOTICES](sim/THIRD_PARTY_NOTICES.md) |
| 1-D Gate 0: evanescent stacks −∞, mass law −0.5 to −4.4, re-steepening, purity | [gate0_summary](sim/results/gate0_summary.md) |
| Mechanism verified; tie, cell-wall and width floors; W ≳ 4h | [review_evanescent](sim/results/review_evanescent.md) |
| BEM vs Whitham and analytic diffraction; asymmetric and fringe edges | [crossplane_validation](sim/results/crossplane_validation.md) |
| Width sets loudness; −34 / −25 dB at 1 km; centreline focus; R_m sensitivity | [3d_summary](sim/results/3d_summary.md), `sim/results/3d_carpet_summary.csv` |
| "Ties make it transparent" was a mesh bug; tie design rule | [line_ties_check](sim/results/line_ties_check.md) |
| Serrations, module joints | [alongtrack_summary](sim/results/alongtrack_summary.md) |
| Wind and turbulence do not change conclusions; headline conservative | [wind_turbulence](sim/results/wind_turbulence.md) |
| Vacuum branch closed; anything holding the faces apart is a stud; invballoon ring correction | [vacuum-balloon-review](vacuum/vacuum-balloon-review.md) |
| Floating-lid baseline, locked-volume rule, H₂ inventory and safety, loads | [structure-concept](structure/structure-concept.md) |
