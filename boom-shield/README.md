# Boom Shield

This directory collects the first-pass acoustic analysis for a bouyant sonic-boom shield concept, to reflect the main wave of sonic booms from overflying aircraft away from their ground tracks.

## Contents

- **[STATUS.md](STATUS.md): current concept status (start here).**
- **[infographic/index.html](infographic/index.html): illustrated explainer (scene, aircraft comparison, loudness scale, cross-sections, mooring).**

Copied from the full [ChatG thread here](https://chatgpt.com/share/69e16de1-e3ac-83e8-9bbd-eb2f741c8bd0).

- [boom-shield-summary.md](boom-shield-summary.md): narrative summary of the findings from the acoustic analysis thread.
- [hypersonic_corridor_points.csv](hypersonic_corridor_points.csv): reference Mach / altitude / density / dynamic-pressure points for SR-71, X-51A, X-43A, and a speculative SR-72 corridor.
- [shield_opacity_required_points.csv](shield_opacity_required_points.csv): required direct-wave opacity versus speed for X-59, SR-71-like, and SR-72-like cases.
- [shield_spacing_envelope_points.csv](shield_spacing_envelope_points.csv): first toy spacing-envelope results using a simple edge-leakage law.
- [shield_spacing_fresnel_points.csv](shield_spacing_fresnel_points.csv): improved spacing-envelope results using a Fresnel / knife-edge diffraction surrogate.
- [handoff/2026-10-20-boomshield-handoff.md](handoff/2026-10-20-boomshield-handoff.md): concept and design handoff (October 2026 working draft).
- [handoff/2026-10-20-review-addendum.md](handoff/2026-10-20-review-addendum.md): review of the handoff with independent checks and the proposed Gate 0.
- [handoff/quick_checks.py](handoff/quick_checks.py): back-of-envelope buoyancy, mass-law, inventory, and load checks behind the review.
- [sim/DESIGN.md](sim/DESIGN.md): design contract for the Gate 0 boom-through-shield-to-ground simulation.
- [sim/results/gate0_summary.md](sim/results/gate0_summary.md): Gate 0 results: ground PLdB with lift-closed membrane/gas/membrane shields in full nonlinear propagation (1-D, infinite shield).
- [sim/results/review_evanescent.md](sim/results/review_evanescent.md): adversarial review of the evanescent gas-layer shielding mechanism and its realistic floors.
- [sim/results/3d_summary.md](sim/results/3d_summary.md): 3-D finite-width corridor study: ground carpets vs shield width, stack, edge treatment, cell walls, ties, aircraft.
- [sim/results/alongtrack_summary.md](sim/results/alongtrack_summary.md): along-track (3-D) study of serrated edges and module joints.
- [sim/results/wind_turbulence.md](sim/results/wind_turbulence.md): wind (moving-medium rays) and boundary-layer turbulence effects on shielded vs unshielded ground carpets.
- [sim/results/psf_sweep.md](sim/results/psf_sweep.md): sweep of gas, thickness, shield altitude and aircraft against the proposed FAA 0.11 psf overland limit.
- [sim/results/secondary_booms.md](sim/results/secondary_booms.md): secondary-boom check (upper-atmosphere returns) with and without the shield against 0.11 psf.
- [sim/results/line_ties_check.md](sim/results/line_ties_check.md): line-tie bug finding, independent Floquet check, and tie-stiffness design rule.
- [vacuum/vacuum-balloon-review.md](vacuum/vacuum-balloon-review.md): screening review of the inverted (vacuum) balloon option: acoustics, structural closure, pumping, failure modes.
- [structure/structure-concept.md](structure/structure-concept.md): structural concept study for a ~1 km-wide H2 shield with acoustically decoupled faces; recommends the "floating-lid" pillow baseline.
- [altitude/altitude-trade.md](altitude/altitude-trade.md): operating-altitude trade study (12-25 km): buoyancy, weather, winds, airspace, aircraft penalty, tethers.
- [network-maps/routes/routes.md](network-maps/routes/routes.md): route ranking (overland km x seats) and Overture time-savings model.
- [network-maps/trunk/trunk.md](network-maps/trunk/trunk.md): shared cross-country trunk (gravity-model O&D, Mach sweep 1.4-3.2), traffic per km and benefit/cost.
- [economics/economics.md](economics/economics.md): first-cut cost, land-rights and benefit model.
- [operations/tethers-build-failure.md](operations/tethers-build-failure.md): faired and buoyant tethers, construction sequence, H2 make-up, failure modes, land use.
- [fatigue/lid-fatigue.md](fatigue/lid-fatigue.md): lid fatigue under repeated boom loading; revised traffic cap and trunk B/C.

## Main result

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
