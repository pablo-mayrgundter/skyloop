# BoomShield
## Concept and design handoff

**Skyloop project · 20 October 2026 · Working draft**

**Document structure:** One-page executive brief, followed by ten design sections intended for approximately ten pages, then references. Actual pagination depends on document formatting.

**Evidence status:** Concept-stage synthesis. No integrated acoustic, structural, buoyancy, or operational demonstration has been established.

---

# Executive brief

**BoomShield is a proposed buoyant acoustic corridor beneath supersonic aircraft.** A segmented surface would intercept and reshape downward-propagating pressure waves, seeking to reduce ground-level sonic boom while aircraft follow a controlled flight path above it. The central question is whether useful ground-noise reduction can be achieved with sufficiently little material, buoyant volume, and station-keeping infrastructure.

The concept inherits Skyloop’s modular envelopes, tension-supported structures, and distributed stabilization. Its optimization target changes from supporting a transport tube to maximizing acoustically effective area per unit mass. The guiding phrase from the discussions is **“extreme areal parsimony.”**

**Architecture under study.** Thin membranes, high-strength tendons, sparse ribs, and replaceable acoustic panels form a slightly curved deck. Gas cells, evacuated structures, or a hybrid arrangement provide buoyancy. Modules connect through compliant joints and acoustic overlaps. Local support stations resist corridor drift, avoiding a continuous structure that accumulates loads over its entire length.

**Two design families remain open.**

| Family | Parameters from prior work | Purpose |
|---|---|---|
| Acoustic reference | 1 km width; approximately 250 m aircraft clearance | Explore direct-wave attenuation and edge leakage |
| Structural module | 200–300 m width; 12–18 m thickness; 250 m module length | Explore lightweight construction and suspension |
| Initial flight study | Mach 1.6–2.5; 75–150 m clearance | Earlier proposed development regime |
| Higher-speed sensitivity | Mach 3.2 and a notional Mach 6 case | Test sensitivity to aircraft signature and shielding |

These are **study cases, not a unified design**. In particular, the narrower structural modules cannot inherit the acoustic performance estimated for the 1 km shield.

**What the earlier work established.** It produced a useful architecture, parameterized attenuation estimates, a structural load-path concept, and a preliminary North American network scaffold. Its strongest organizing condition is an operating window:

**Minimum operational clearance < maximum acoustically useful clearance.**

**What this handoff adds.** Three checks materially change how the earlier results should be used:

1. **Geometry:** The earlier Mach-angle expression does not establish required cross-track shield width. A moving-source, three-dimensional calculation must replace it.
2. **Buoyancy:** At approximately 20 km altitude, a 15 m-thick elliptical helium envelope provides only about **0.90 kg/m²** of ideal lift capacity for everything except the lifting gas. The proposed **4 kg/m²** system cannot float on that volume alone.
3. **Acoustics:** Matching width-to-clearance ratios and scaling pressure amplitudes from the X-59 do not demonstrate equivalent ground loudness. Full waveforms, diffraction, atmospheric propagation, and structural reradiation remain unresolved.

**Recommended next decision.** Fund a bounded feasibility study before selecting a corridor or full-scale prototype. It should deliver a validated finite-width acoustic model, a complete mass-and-lift budget, a dynamic clearance envelope, and a cost model tied to useful traffic. Preserve gas, hybrid, and vacuum options, but use gas-supported construction as the clearest comparison baseline.

**Handoff priority:** Determine whether any configuration simultaneously closes acoustics, buoyancy, and clearance. That intersection—not an isolated favorable calculation—is the first credible feasibility result.

---

# Design document

## 1. Purpose, provenance, and concept evolution

BoomShield addresses the acoustic constraint on overland supersonic transportation by introducing infrastructure between aircraft and the ground. The intended benefit is a shared corridor that serves multiple aircraft and origin–destination pairs. Unlike Skyloop’s enclosed transport tube, BoomShield leaves aircraft in open air and changes the propagation environment beneath them.

The retrieved March 2026 discussion described a porous ribbon near 20 km altitude, potentially spanning major corridors with two directions of travel. Candidate mechanisms included scattering, refraction, disruption of coherent wavefronts, and optional thermal gradients. Aerogel-supported balloon grids or planks were suggested as a construction direction, including possible insulation and maintenance functions.

The April acoustic work reframed the idea around suppressing the direct centerline wave with a wide horizontal shield. It introduced material transmission, edge leakage, aircraft-to-shield spacing, and a comparison against the X-59 low-boom design reference. The structural discussions then emphasized flattened buoyant envelopes, tension-supported decks, local stabilization, and acoustic performance per unit areal mass.

A subsequent network sketch treated the system as flight infrastructure: shared high-speed trunks, slower connecting ramps, and climb/acceleration paths constrained by aircraft dynamic pressure.

This handoff draws primarily on the current repository’s acoustic summary, segment-analysis notes, associated CSV datasets, and North America network README [P1–P4]. Retrieved conversation context supplements those files. It does not establish that every chat in the project folder has been recovered.

Three evidence categories apply throughout:

- **Inherited concept:** An idea or parameter recorded in the previous work.
- **Prior model result:** An output of the earlier simplified calculations.
- **Handoff assessment:** A new interpretation, calculation, correction, or proposed next step.

The distinction matters because several earlier conclusions were stronger than their models support. Those conclusions remain part of the project history, but should not become requirements or performance guarantees without further analysis.

The immediate purpose is a technical handoff for acoustics, lightweight structures, flight operations, and systems engineering. It is not yet a construction specification or investment-grade business case.

## 2. System definition and requirements

The system boundary includes the acoustic surface, support structure, buoyancy system, joints, stabilization, sensing, power, maintenance equipment, support stations, and ground interfaces. Aircraft remain separate vehicles, but their navigation, wake, source signature, and emergency behavior are essential design inputs.

A working mission statement is:

**Reduce sonic-boom exposure across a specified ground region by operating compatible aircraft above a maintained buoyant corridor, at an acceptable infrastructure cost and availability.**

“Specified ground region” must be explicit. Reducing pressure directly beneath the centerline while creating strong lateral peaks would not establish useful corridor-wide suppression.

The next study should resolve the following requirements rather than assume them:

| Requirement | Quantity to define |
|---|---|
| Acoustic outcome | Ground waveform and loudness limits; protected width; off-axis limits |
| Aircraft compatibility | Source signatures, dimensions, mass, Mach, altitude, maneuver envelope |
| Geometry | Shield width, thickness, curvature, module length, and clearance |
| Structural behavior | Maximum deflection, gust response, fatigue life, damaged-state behavior |
| Buoyancy | Complete mass budget, lift reserve, gas or vacuum performance |
| Operations | Weather availability, lane separation, entry/exit and abort procedures |
| Maintenance | Inspection interval, replacement time, leakage and repair limits |
| Economics | Installed cost per route-km and annual cost per useful movement |

NASA’s X-59 offers a relevant reference because its design objective is approximately 75 PLdB in a standard atmosphere. That is a reference for shaped-waveform research, not an adopted BoomShield acceptance threshold. [R1]

For every acoustic test, preserve the complete pressure-time history. Report positive and negative peaks, duration, rise time, an energy-related measure, and an appropriate perceived-level metric. State the receiver location and atmosphere.

The inherited objective of “boom reduction per kg/m²” remains valuable, but should be accompanied by protected area and cost. Otherwise, an extremely light treatment could score well while producing little useful attenuation.

The decision framework should therefore compare **acoustic benefit, supported mass, operating availability, and lifecycle cost together**. A design cannot compensate for negative lift margin by scoring well on acoustic efficiency.

## 3. Geometry and the operating envelope

Use a common notation across all subsequent work:

| Symbol | Definition |
|---|---|
| H | Shield altitude above the reference ground |
| h | Aircraft clearance above the relevant shield surface |
| W | Total acoustically active width |
| T | Maximum buoyant-envelope thickness |
| Lₘ | Module length along the corridor |
| M | Aircraft Mach number |
| σ | Supported non-gas system mass divided by shield planform area |
| λ | Acoustic wavelength or explicitly defined spectral scale |

The original studies used both *E* and *h* for clearance. This handoff uses *h*. Ground altitude, envelope centerline, deck elevation, and uppermost obstruction should remain separate in the geometry model.

**Correction to the inherited width estimate.** The segment notes used:

W_req ≈ 2h / √(M² − 1)

and inferred that higher Mach directly reduces required shield width. The Mach angle relation, sin μ = 1/M, is valid for the idealized Mach cone. Its application to a cross-track shield width was not established.

For a point-source Mach cone in uniform air, let x be downstream distance behind the aircraft and y the cross-track coordinate. At a horizontal plane h below the aircraft:

y² + h² = x² / (M² − 1).

The centerline intersection begins at x = h√(M² − 1). The cross-track intersection then changes with downstream position. Consequently, h and M alone do not specify a single total width that intercepts the complete relevant wavefield.

The older width tables should be retained as historical calculations, but their geometric-coverage ratio should not govern design acceptance. Replacement analysis must include aircraft length, distributed sources, the receiver region, finite shield length, and atmospheric refraction.

The operational window remains useful:

h_min,ops < h < h_max,acoustic.

A preliminary clearance budget can combine hard clearance, aircraft motion, shield motion, bias, and fault allowance. If statistical errors are combined, account for their covariance; independence is an assumption, not a default fact.

The upper bound must come from a validated acoustic model. The lower bound must come from dynamic relative-motion analysis. Neither is established by the current notes.

## 4. Acoustic model and interpretation of prior results

The acoustic summary explored a 1 km-wide shield at approximately 250 m clearance. Its surrogate combined intrinsic material transmission with edge leakage:

A_total ≈ A_mat + K√(λ_eff h)/W.

Here A represents pressure-amplitude fractions and K is a model coefficient. This is a screening approximation. Real transmitted, diffracted, and reradiated fields have frequency-dependent amplitudes and phases.

The saved dataset gives the following representative results for W = 1 km, A_mat = 0.05, and λ_eff = 30 m:

| Prior model case | Required amplitude fraction | Reported maximum clearance |
|---|---:|---:|
| Mach 3.2, optimistic source | 0.4011 | 5,069.5 m |
| Mach 3.2, source penalty ×2 | 0.2006 | 932.1 m |
| Mach 6, optimistic source | 0.1631 | 525.8 m |
| Mach 6, source penalty ×2 | 0.0815 | 40.9 m |

These values reproduce the earlier CSV, not observed performance. The repository’s proportionality expression omits the full calibration needed to reconstruct every number; that calibration should be made explicit before extending the model.

Three limitations are especially consequential.

First, the width-to-clearance analogy between an open-air case and a shielded case preserves only one dimensionless ratio. Acoustic similarity also depends on wavelength, source geometry, impedance, and receiver distance.

Second, the opacity table treats differences in perceived level as equivalent pressure-amplitude attenuation. The relation −20 log₁₀(A) describes an amplitude ratio; it does not by itself calculate perceived loudness for a waveform whose spectrum and duration have changed.

Third, edge diffraction is not the only residual contribution. Gaps, finite ends, off-axis paths, panel vibration, and reflected fields can matter.

The replacement workflow should propagate a specified aircraft near-field signature through the shield interaction and then through representative atmospheres to a ground receiver grid. NASA’s sonic-boom prediction work provides relevant propagation benchmarks and demonstrates the importance of atmospheric conditions. [R2]

The practical value of the old model is sensitivity identification: performance depends strongly on transmission floor, spectral content, clearance, and aircraft shaping. It does not yet establish a viable operating envelope.

## 5. Acoustic surface and material architecture

The design discussions converge on a surface that reshapes a wavefront using minimal mass. Candidate actions include reflection, scattering, temporal spreading, controlled transmission, and dissipation. Each mechanism needs a measurable energy and waveform accounting.

A useful test hierarchy is:

| Surface concept | Intended benefit | Key uncertainty |
|---|---|---|
| Continuous tensioned membrane | Low mass; simple construction | Transmission and vibration at relevant frequencies |
| Curved or faceted panels | Redirect outgoing waves | Lateral focusing and reflected aircraft exposure |
| Staggered, overlapping layers | Reduce direct gaps; add path diversity | Added mass and multiple scattering |
| Porous/resonant treatment | Combine airflow with acoustic filtering | Bandwidth, size, and pressure dependence |
| Serrated or irregular edges | Redistribute edge radiation | Whether maximum ground exposure improves |

“Wavefront mangler” captures the design intent, but does not identify a proven mechanism. Spatially redistributing energy can reduce a local peak without eliminating the disturbance. Atmospheric propagation may also alter or reconcentrate the modified waveform.

Relevant acoustic metasurface research has experimentally demonstrated ventilated barriers using combined dissipation and interference. One example reports more than 90% incident-energy blocking over 650–2,000 Hz in a 53 mm structure. Those frequency and laboratory conditions are not a demonstration of lightweight stratospheric sonic-boom shielding. [R4]

Inherited material candidates include polyethylene or PET films, aramid tendons, sparse foam or composite ribs, and selective damping layers. They remain candidates: no qualified material stack or life model was recovered.

Separate the mass of the acoustically active deck from the complete airborne system. The earlier 0.2–0.5 kg/m² aspiration and the 2–8 kg/m² structural comparison range describe different levels of ambition and potentially different boundaries. They must not be interchanged.

Aerogel remains a possible local insulation, spacer, or damping material. Filling a buoyant volume with even a very low-density solid can consume its lift. Retain the original construction idea, but evaluate aerogel by measured subsystem benefit and mass rather than assuming it provides economical large-area support.

The first material experiments should measure transmission, reflection, vibration, and damage across the actual pulse spectrum and incidence angles.

## 6. Buoyancy and mass closure

**New handoff calculation:** The proposed shallow envelopes do not support the earlier total areal-mass targets at stratospheric altitude without additional volume.

For an elliptical envelope of width W, thickness T, and length Lₘ:

V = πWTLₘ/4.

Dividing ideal buoyant capacity by planform area WLₘ gives:

σ_capacity = (ρ_air − ρ_gas)πT/4.

This capacity is available for envelope, deck, structure, controls, tether burden, and reserve. It already subtracts lifting-gas mass. It does not subtract any hardware mass.

Using rounded standard-atmosphere densities and helium at equal pressure and temperature yields:

| Altitude | Air density | Helium capacity, T = 15 m | Ideal vacuum upper bound |
|---|---:|---:|---:|
| 18 km | 0.121 kg/m³ | 1.23 kg/m² | 1.42 kg/m² |
| 20 km | 0.089 kg/m³ | 0.90 kg/m² | 1.05 kg/m² |
| 25 km | 0.040 kg/m³ | 0.41 kg/m² | 0.47 kg/m² |

These are handoff calculations using the standard-atmosphere reference [R5], not complete balloon designs.

Consider the inherited S2 module: 250 m wide, 15 m thick, and 250 m long. Its planform is 62,500 m² and volume approximately 736,000 m³. At 20 km, its ideal helium capacity is approximately **56 tonnes**. At σ = 4 kg/m², its supported system mass would be **250 tonnes**.

The resulting shortfall exists before adding reserve. Ideal vacuum increases capacity only to approximately 65 tonnes, before accounting for the structure needed to resist external pressure.

At the same altitude, helium support for 4 kg/m² would require roughly 66 m maximum elliptical thickness under these idealized assumptions. Adding reserve or extra uncounted components increases that requirement.

This exposes a geometric trade: holding WT constant preserves volume per unit length, but widening the envelope spreads that lift over more planform area. A wider, thinner section does not preserve lift capacity per square metre.

Available responses are lower mass, greater thickness, separate buoyancy bodies, lower operating altitude, or a revised architecture. None has yet been selected. The next mass budget must include every airborne subsystem and the vertical forces imposed by station-keeping.

## 7. Structural load paths, support options, and durability

The structural concept remains a segmented, tension-supported acoustic raft. End frames maintain local geometry; catenary nets distribute deck load; lightweight fairings manage envelope shape; compliant connections limit load transfer between modules.

The inherited load path is:

**Deck → suspension net → frames and buoyancy structure → local stabilization interfaces.**

Tension-dominated support is attractive because it avoids distributing bending stiffness across an enormous area. It does not eliminate compression or stability requirements in frames, joints, and evacuated structures.

The earlier suspension estimate used:

H ≈ qs²/(8f),

where q is load per unit span, s is span, and f is sag. For a one-metre-wide strip with σ = 4 kg/m², q is about 39 N/m. With s = 125 m and f = 15 m, horizontal tension is approximately 5.1 kN.

This is a reasonable static screening calculation for that idealized strip. It does not cover gusts, pretension, fatigue, load concentration, frame buckling, acoustic impulses, or damaged conditions.

Three support branches remain:

| Branch | Advantage | Principal unresolved issue |
|---|---|---|
| Gas-dominant | Straightforward buoyancy accounting | Gas volume, containment, leakage, and drag |
| Evacuated structure | Avoids lifting-gas inventory | External-pressure stability and structural mass |
| Hybrid | Potential separation of shape, lift, and trim functions | Complexity and ambiguous allocation of loads |

The original notes favored hybrid support. This handoff preserves that preference as an inherited hypothesis, not a demonstrated optimum.

A low-pressure flexible shell is not automatically stable because its material is strong in tension. External pressure tends to collapse it. An evacuated design needs a defined compression-bearing or otherwise stabilizing architecture, with imperfections and buckling included.

NASA scientific balloons provide relevant precedent for thin polyethylene envelopes and distributed load-bearing members, but free-flying balloons do not demonstrate a stationary, shallow acoustic corridor. [R6]

Durability testing should cover seam strength, tear propagation, creep, thermal cycling, UV exposure, gas permeability, repeated pulses, and repair quality. Replaceable panels and isolated buoyancy cells are valuable only if their interfaces do not dominate mass or create major acoustic leakage.

## 8. Tethering, station-keeping, and aircraft operations

The discussions identify three tether functions: internal suspension, shape control, and external position-keeping. They should remain separate in the mass model and fault analysis.

The proposed architectural rule is sound as a starting point: stabilize each module or small cluster locally rather than passing all lateral force through an entire corridor. However, “local” at approximately 20 km altitude still implies a substantial three-dimensional support system.

A first crosswind estimate per unit corridor length is:

F′ ≈ ½ρU²C_D T.

For an illustrative 20 km case with ρ = 0.089 kg/m³, U = 30 m/s, C_D = 1, and T = 15 m, the side-force estimate is approximately 600 N/m, or 150 kN across a 250 m module. This is an illustrative load case, not a site wind specification.

It excludes lift on the broad surface, gust dynamics, support-body drag, and tether drag. Ground-connected tethers traverse denser lower air and cannot be evaluated using only the density at the shield.

A station-keeping model needs the wind profile along each tether, cable weight and taper, tension limits, anchor geometry, aerodynamic forces, and the vertical component imposed on the buoyancy system. Purely airborne positioning trades some ground infrastructure for propulsion, power, and control mass.

Aircraft operations introduce a separate clearance problem. The shield’s dynamic surface—not a nominal altitude—is the reference obstruction. Corridor control must estimate shape, motion, and uncertainty sufficiently far ahead for an aircraft to respond.

The two-direction concept also needs explicit lane separation, wake analysis, merge rules, and failed-segment routing. A missing panel or displaced module cannot be treated only as a maintenance issue if it changes the acoustic or collision envelope.

Initial operational cases should include normal passage, rejected entry, wind-limit closure, one-cell loss, one-tether loss, and communication loss. Each case needs a defined outcome for both the infrastructure and approaching traffic.

The required operating margin is a result of these analyses, not a discretionary percentage appended to the acoustic spacing estimate.

## 9. Corridor network and economic model

The North America scaffold treats BoomShield as a network of shared trunks and ramps. Its starting hypothesis is LAX–DEN–ORD–JFK, with feeder concepts involving SEA, SFO, DFW, ATL, MIA, YYZ, and MEX [P4].

This is a sketch for evaluation, not a traffic-optimized route. The saved README explicitly states that a recent-year demand dataset is still needed.

The earlier stage-length screening was:

| Great-circle distance | Original screening category |
|---|---|
| Below 700 km | Weak case for dedicated Mach infrastructure |
| 700–1,500 km | Borderline or regional connector |
| 1,500–3,000 km | Stronger corridor candidate |
| Above 3,000 km | Strong trunk candidate |

These thresholds are placeholders. A useful network model must include climb, acceleration, route deviation, airport access, congestion, weather availability, and feasible aircraft performance. Cruise speed alone overstates passenger time savings.

Aircraft dynamic pressure should constrain altitude and speed along the route:

q_aircraft = ½ρV².

This is distinct from the wind-loading calculation for the shield. A corridor spline must also respect turn radius, climb capability, shield construction geometry, and permissible entry/exit behavior.

A first infrastructure cost equation is:

C/km = 1,000Wc_area + C_buoyancy/km + C_stations/km + C_installation/km + C_controls/km.

Avoid double counting: if c_area already includes buoyancy hardware or controls, remove those separate terms. Annual cost should include replacement, inspection, gas replenishment, power, staffing, land access, and downtime.

Scale alone makes unit cost important. A 250 m-wide corridor has 250,000 m² of planform per kilometre. Every hypothetical $10/m² contributes $2.5 million/km before other categories. At 1 km width, the same allowance contributes $10 million/km. These are arithmetic sensitivities, not market quotations.

Compare network options using annualized cost per useful movement and cost per passenger-hour saved. Shared trunks can improve utilization, but impose detours, merges, and common-mode outage exposure.

No credible corridor price or return estimate follows from the current evidence. The handoff requirement is a transparent model whose inputs can later be replaced with supplier estimates and measured maintenance data.

## 10. Validation plan and handoff decisions

Development should proceed through evidence gates that expose failure early.

**Gate 1 — Reproduce and correct the analytical baseline.** Preserve the original CSVs. Reconstruct their equations and coefficients, document units, distinguish pressure-level and perceived-level quantities, and retire the unsupported cross-track width rule. Create one versioned parameter registry shared by acoustic and structural models.

**Gate 2 — Establish acoustic benefit.** Compare an unshielded source against rigid, compliant, porous, and segmented surfaces. Include finite width, finite length, lateral receivers, joints, and reflected fields. Validate simple cases against analytical solutions or published diffraction experiments before relying on a large simulation. NASA’s building-scale sonic-boom diffraction work is relevant methodological precedent, not direct validation of an airborne shield. [R3]

**Gate 3 — Close mass and lift.** Produce a subsystem budget with supplier-based or measured component masses, envelope volume, atmosphere, gas state, reserve, and tether forces. Reject configurations that remain negatively buoyant under their stated assumptions.

**Gate 4 — Demonstrate coupled dynamics.** Combine flexible-structure motion, station-keeping, gusts, and aircraft-relative navigation. Test representative faults. Determine whether an operational clearance interval remains after acoustic and motion uncertainty are included.

**Gate 5 — Instrumented demonstration.** Begin with material coupons and controlled pulse tests. Progress to a suspended finite array and then a buoyant structural demonstrator. Scaling must preserve the relevant acoustic and structural dimensionless parameters; a convenient small prototype is not automatically representative.

The receiver array should measure both the intended quiet region and displaced energy. NASA’s SR-71 measurement work is useful for sensor and waveform considerations, but the cited experiment covered Mach 1.25–1.60; it does not validate the project’s Mach 3.2 reference case. [R7]

| Workstream | Required handoff output | Decision enabled |
|---|---|---|
| Acoustics | Validated ground-signature maps | Does the surface provide useful suppression? |
| Structures/buoyancy | Complete mass and stability model | Can it support itself with reserve? |
| Controls/operations | Dynamic clearance and fault envelope | Can compatible aircraft use it reliably? |
| Network/economics | Demand-weighted lifecycle model | Is a shared corridor worth building? |

No owners, budget, schedule, or acceptance thresholds have been assigned in the recovered material. Those are the next program decisions.

The inherited concept worth carrying forward is a modular, very light, controlled acoustic surface. The next deliverable should establish a credible intersection of attenuation, lift, and clearance before committing to its full-scale geometry.

---

# References and related work

## Project source material

The repository files below were read for this handoff. Links use the repository snapshot surfaced during retrieval.

**[P1] Skyloop — Boom Shield README and acoustic summary.**  
Defines the 1 km shield, nominal 250 m clearance, direct-wave suppression framing, and acoustic/operational spacing window.

[README](https://github.com/pablo-mayrgundter/skyloop/blob/cd579d5532cab87c1c1c6041ef5e403fe3762cd1/boom-shield/README.md) · [Acoustic summary](https://github.com/pablo-mayrgundter/skyloop/blob/cd579d5532cab87c1c1c6041ef5e403fe3762cd1/boom-shield/boom-shield-summary.md)

**[P2] Skyloop — Boom Shield segment analysis: thread notes.**  
Primary source for module geometry, areal-mass targets, catenary suspension, support branches, materials, and local tethering.

[Segment-analysis notes](https://github.com/pablo-mayrgundter/skyloop/blob/cd579d5532cab87c1c1c6041ef5e403fe3762cd1/boom-shield/segment-analysis/boom-shield-segment-analysis-thread-notes.md)

**[P3] Skyloop — Acoustic and flight-reference datasets.**  
Original surrogate outputs. These are project calculations rather than experimental validation.

[Fresnel spacing estimates](https://github.com/pablo-mayrgundter/skyloop/blob/cd579d5532cab87c1c1c6041ef5e403fe3762cd1/boom-shield/shield_spacing_fresnel_points.csv) · [Required opacity estimates](https://github.com/pablo-mayrgundter/skyloop/blob/cd579d5532cab87c1c1c6041ef5e403fe3762cd1/boom-shield/shield_opacity_required_points.csv) · [Flight-reference points](https://github.com/pablo-mayrgundter/skyloop/blob/cd579d5532cab87c1c1c6041ef5e403fe3762cd1/boom-shield/hypersonic_corridor_points.csv)

**[P4] Skyloop — North America network maps README.**  
Source for shared trunks, ramps, stage-length screening, and the preliminary airport network. Demand weighting remains unfinished.

[Network methodology](https://github.com/pablo-mayrgundter/skyloop/blob/cd579d5532cab87c1c1c6041ef5e403fe3762cd1/boom-shield/network-maps/america-north/README.md)

**Supplementary provenance:** Retrieved discussions dated March 27 and April 16–17, 2026 supplied the early porous-ribbon framing, extreme-areal-parsimony preference, and aerogel/support explorations. Where those recollections overlap with repository content, this handoff uses the repository as the primary technical record.

## External related work

**[R1] NASA. “Simulations of X-59 Low-Booms Propagated Through Measured Atmospheric Profiles in Galveston, Texas.”** NTRS document 20205005688.  
Relevant to the 75 PLdB design reference and the need to evaluate actual atmospheric profiles. It does not establish a universal permissible loudness or a BoomShield performance target.  
[Source](https://ntrs.nasa.gov/citations/20205005688)

**[R2] Rallabhandi, S. K. “Propagation Analysis of the 3rd Sonic Boom Prediction Workshop Cases using sBOOM.”** AIAA SciTech Forum, 2021. NTRS document 20205010207.  
Relevant to propagation-model benchmarking and changes in boom footprints under nonstandard atmospheres.  
[Source](https://ntrs.nasa.gov/citations/20205010207)

**[R3] Barger, R. L.; Beasley, W. D.; Brooks, J. D. “Laboratory Investigation of Diffraction and Reflection of Sonic Booms by Buildings.”** NASA TN D-5830, 1970.  
Relevant precedent for transient-wave diffraction, reflection, and experimental validation. Building-scale geometry differs materially from the proposed airborne corridor.  
[Source](https://ntrs.nasa.gov/citations/19700021539)

**[R4] Dong, R.; Mao, D.; Wang, X.; Li, Y. “Ultrabroadband Acoustic Ventilation Barriers via Hybrid-Functional Metasurfaces.”** *Physical Review Applied* 15, 024044, 2021. DOI: 10.1103/PhysRevApplied.15.024044.  
Relevant to combined dissipation and interference in ventilated barriers. Its reported frequency band and laboratory conditions should not be extrapolated directly to sonic-boom shielding.  
[Paper](https://doi.org/10.1103/PhysRevApplied.15.024044) · [Open manuscript](https://arxiv.org/abs/2010.08160)

**[R5] NOAA / NASA / USAF. *U.S. Standard Atmosphere, 1976*.** NOAA-S/T-76-1562; NASA TM X-74335, 1976.  
Baseline atmospheric reference for buoyancy and aerodynamic calculations. Operational design requires location- and season-specific conditions beyond a standard atmosphere.  
[Source](https://ntrs.nasa.gov/citations/19770009539)

**[R6] NASA Balloon Program. “Types of Scientific Balloons”; “‘Seal Team 6’: Constructing a Super Pressure Balloon.”**  
Relevant to thin-film envelopes, load-bearing members, fabrication, and balloon architecture. These precedents do not demonstrate stationary corridor operation.  
[Balloon types](https://www.nasa.gov/scientificballoons/types-of-balloons/) · [Construction](https://science.nasa.gov/blogs/super-pressure-balloon/2017/04/05/seal-team-6-constructing-a-super-pressure-balloon/)

**[R7] Norris, S. R.; Haering, E. A., Jr.; Murray, J. E. “Ground-Based Sensors for the SR-71 Sonic Boom Propagation Experiment.”** Published in *The 1995 NASA High-Speed Research Program Sonic Boom Workshop*, 1996. NTRS document 19960055059.  
Relevant to measurement methods, N-wave signatures, and atmospheric variability. The reported flight conditions were Mach 1.25–1.60 at 30,000–48,000 ft.  
[Source](https://ntrs.nasa.gov/citations/19960055059)

**Reference interpretation:** These sources support component physics, methods, and adjacent technologies. None of the reviewed sources demonstrates the integrated BoomShield concept. Novelty, patentability, and freedom to operate have not been assessed.
