# Boom Shield

BoomShield is a concept for a long, flat, tethered layer of hydrogen floating at ~20 km under a supersonic flight corridor, so that the sonic boom of the aircraft flying just above it decays in the gas instead of reaching the ground. This directory holds the screening studies (simulation only, with unverified web data): acoustics, structure, operations, network, economics, higher Mach, growth and global corridors. Start with **[STATUS.md](STATUS.md)** and the illustrated **[explainer page](infographic/index.html)**.

**Phase-end verdict (October 2026).** The acoustics pass a screen: ~55 m of H₂, 1.5–2 km wide, meets the FAA's proposed 0.11 psf limit at Mach 2, and buildability is open rather than refuted. The economics do not close: with owned plants (~$14/m²) and a light M2.0–2.2 aircraft the 3,614 km US spine (the best corridor in the world) reaches a gross benefit/cost of ~0.9 against Boomless Cruise, but net of the extra fuel and CO₂ it stays well below 1 against Boomless and negative against subsonic, so the next step is a re-run around 2029–31 once the FAA rule and Overture's real economics are known.

## Contents

- **[STATUS.md](STATUS.md): current concept status (start here).**
- **[infographic/index.html](infographic/index.html): illustrated explainer (scene, aircraft comparison, loudness scale, cross-sections, mooring, network and economics).**

### Acoustics and simulation

- [sim/DESIGN.md](sim/DESIGN.md): design contract for the Gate 0 boom-through-shield-to-ground simulation.
- [sim/DESIGN-3D.md](sim/DESIGN-3D.md): design contract for the 3-D finite-width corridor model (edge diffraction, cell walls, ties).
- [sim/results/gate0_summary.md](sim/results/gate0_summary.md): Gate 0 results: ground PLdB with lift-closed membrane/gas/membrane shields in full nonlinear propagation (1-D, infinite shield).
- [sim/results/review_evanescent.md](sim/results/review_evanescent.md): adversarial review of the evanescent gas-layer shielding mechanism and its realistic floors.
- [sim/results/crossplane_validation.md](sim/results/crossplane_validation.md): validation of the cross-plane solver against Whitham and analytic checks.
- [sim/results/3d_summary.md](sim/results/3d_summary.md): 3-D finite-width corridor study: ground carpets vs shield width, stack, edge treatment, cell walls, ties, aircraft.
- [sim/results/alongtrack_summary.md](sim/results/alongtrack_summary.md): along-track (3-D) study of serrated edges and module joints.
- [sim/results/line_ties_check.md](sim/results/line_ties_check.md): line-tie bug finding, independent Floquet check, and tie-stiffness design rule.
- [sim/results/wind_turbulence.md](sim/results/wind_turbulence.md): wind (moving-medium rays) and boundary-layer turbulence effects on shielded vs unshielded ground carpets.
- [sim/results/psf_sweep.md](sim/results/psf_sweep.md): sweep of gas, thickness, shield altitude and aircraft against the proposed FAA 0.11 psf overland limit.
- [sim/results/secondary_booms.md](sim/results/secondary_booms.md): secondary-boom check (upper-atmosphere returns) with and without the shield against 0.11 psf.
- [sim/results/focus_booms.md](sim/results/focus_booms.md): focused booms (superbooms) from acceleration and turns, with and without the shield.
- [sim/results/onramp_shield.md](sim/results/onramp_shield.md): on-ramp shield geometry and merge rules for feeder traffic.

### Structure and operations

- [structure/structure-concept.md](structure/structure-concept.md): structural concept study for a ~1 km-wide H2 shield with acoustically decoupled faces; recommends the "floating-lid" pillow baseline.
- [vacuum/vacuum-balloon-review.md](vacuum/vacuum-balloon-review.md): screening review of the inverted (vacuum) balloon option: acoustics, structural closure, pumping, failure modes.
- [altitude/altitude-trade.md](altitude/altitude-trade.md): operating-altitude trade study (12-25 km): buoyancy, weather, winds, airspace, aircraft penalty, tethers.
- [fatigue/lid-fatigue.md](fatigue/lid-fatigue.md): lid fatigue under repeated boom loading; retires the ~100 flights/day cap (revised cap 1,000 / 10,000 / none). Literature notes: [fatigue/sources_notes.md](fatigue/sources_notes.md).
- [operations/tethers-build-failure.md](operations/tethers-build-failure.md): faired and buoyant tethers, construction sequence, H2 make-up, failure modes, land use.

### Network, traffic and aircraft

- [network-maps/america-north/README.md](network-maps/america-north/README.md): first-pass North America corridor framing (trunks and ramps).
- [network-maps/routes/routes.md](network-maps/routes/routes.md): route ranking (overland km x seats) and Overture time-savings model.
- [network-maps/trunk/trunk.md](network-maps/trunk/trunk.md): shared cross-country trunk (gravity-model O&D, Mach sweep 1.4-3.2), traffic per km and benefit/cost.
- [network-maps/coverage/coverage.md](network-maps/coverage/coverage.md): total US long-haul traffic (DOT Table 1 based), addressable share, spine coverage, network-size curve.
- [network-maps/demand/demand.md](network-maps/demand/demand.md): induced demand from shorter trips (elasticities, day-trip thresholds, rule-of-half), revised B/C and side effects. Its [joint_costdown.csv](network-maps/demand/joint_costdown.csv) is one tether row; the two-row B/C is in integration's `bcr.csv`.
- [network-maps/highmach/highmach.md](network-maps/highmach/highmach.md): higher-Mach aircraft (M1.7-4.0, adaptive-cycle/TBCC), door-to-door decomposition, shield vs Mach, B/C gross and net of fuel/CO2.
- [network-maps/global/global.md](network-maps/global/global.md): worldwide long-spine corridor ranking vs the US spine; first mover.

### Economics and timing

- [economics/economics.md](economics/economics.md): first-cut cost, land-rights and benefit model.
- [economics/cost-down/cost-down.md](economics/cost-down/cost-down.md): unit-cost reduction levers and dual-use revenue; scenario B/C (its own joint-demand section is superseded).
- [economics/integration/integration.md](economics/integration/integration.md): vertical integration of the main cost drivers (make vs buy, plants, learning curves) and two-tether-row B/C.
- [economics/growth/growth.md](economics/growth/growth.md): traffic-growth projections, ROI by opening year (gross and net of fuel/CO2), decision triggers and a staged pilot plan.

### Origin and handoff (historical)

Copied from the full [ChatG thread here](https://chatgpt.com/share/69e16de1-e3ac-83e8-9bbd-eb2f741c8bd0).

- [boom-shield-summary.md](boom-shield-summary.md): narrative summary of the findings from the acoustic analysis thread.
- [hypersonic_corridor_points.csv](hypersonic_corridor_points.csv): reference Mach / altitude / density / dynamic-pressure points for SR-71, X-51A, X-43A, and a speculative SR-72 corridor.
- [shield_opacity_required_points.csv](shield_opacity_required_points.csv): required direct-wave opacity versus speed for X-59, SR-71-like, and SR-72-like cases.
- [shield_spacing_envelope_points.csv](shield_spacing_envelope_points.csv): first toy spacing-envelope results using a simple edge-leakage law.
- [shield_spacing_fresnel_points.csv](shield_spacing_fresnel_points.csv): improved spacing-envelope results using a Fresnel / knife-edge diffraction surrogate.
- [segment-analysis/boom-shield-segment-analysis-thread-notes.md](segment-analysis/boom-shield-segment-analysis-thread-notes.md): thread notes merging the Skyloop segment logic with the shield's structural and material thread.
- [handoff/2026-10-20-boomshield-handoff.md](handoff/2026-10-20-boomshield-handoff.md): concept and design handoff (October 2026 working draft).
- [handoff/2026-10-20-review-addendum.md](handoff/2026-10-20-review-addendum.md): review of the handoff with independent checks and the proposed Gate 0.
- [handoff/quick_checks.py](handoff/quick_checks.py): back-of-envelope buoyancy, mass-law, inventory, and load checks behind the review.

## Main result (historical)

> **Superseded by [STATUS.md](STATUS.md).** This was the first-pass result from the original thread; the simulation chain, the 0.11 psf check and the network and economics studies replace it.

The analysis suggests that a **Width W=1 km-wide boom shield** with aircraft flying on the order of **Elevation E=250 m above the shield** may be acoustically plausible for suppressing the **direct centerline boom wave**, especially for **SR-71-class** overland flight and potentially for **optimistic SR-72-class** scenarios.

The central viability condition is:

`E_min,ops < E_max,acoustic`

where:

- `E_max,acoustic` is the largest aircraft-to-shield spacing allowed by the acoustic shielding model.
- `E_min,ops` is the minimum safe spacing required by clearance, gusts, navigation error, and shield motion.

That means the next design thread should focus on shield material, buoyancy, tethering, structural motion, and control margins.

<img width="640" alt="image" src="https://github.com/user-attachments/assets/ea8d2b0b-5731-4b88-b698-0fdf8acaab59" />
<img width="640" alt="image" src="https://github.com/user-attachments/assets/e04fadad-679b-497b-a106-292406f3ee56" />
<img width="640" alt="image" src="https://github.com/user-attachments/assets/c29d7f8f-a720-4306-b50d-39c54013f077" />
