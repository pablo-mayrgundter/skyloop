# BoomShield — concept status

**Skyloop · 22 October 2026 · status after the acoustics, 0.11 psf, superboom, on-ramp, fatigue, operations, network, economics and induced-demand rounds.**
Everything here is screening level: simulation only (no hardware, no flight test), linear shield physics, unshaped sources, and web figures taken from search-result summaries that were not verified at the source. Unit costs marked [U] in the studies are engineering estimates.

## 1. Bottom line

- **The acoustics pass a screen, and nothing found since undoes them.** At 20 km, a hydrogen layer about **55 m thick and 1.5–2 km wide** meets the FAA's proposed **0.11 psf** limit under a Concorde-class aircraft at Mach 2. Superbooms from acceleration vanish over it, on-ramps need almost no extra shield, and lid fatigue is a detail-design requirement, not a traffic cap.
- **Buildability is open, not refuted.** Sealed 250 m modules can be launched and joined. The tether field is twice the size first assumed (two rows, ~2,700 tethers), and the worst failure is a bay rising into the flight level about 100 s after losing its hold-down.
- **The best network is a 3,614 km spine, LAS–DEN–OMA–CHI–CLE–NYC.** It carries ~81 M of the ~217 M addressable passengers a year and **~37 %** of the time the shield could save.
- **It does not pay against Boomless Cruise under any combination studied.** The best case is a B/C of **0.55**: aggressive cost-down, high induced demand and 100 % capture. Against subsonic flight it exceeds 1 only with aggressive cost-down at ~100 % capture, which implies fares near subsonic parity. Its carbon cost then exceeds its time benefit.
- **Three things decide the case:** whether Boomless Cruise is allowed and reliable overland, whether supersonic fares can reach parity, and whether all-in cost reaches ~$20/m².

## 2. What passes

**Evidence chain.** The `boomsim` package ([DESIGN](sim/DESIGN.md), [DESIGN-3D](sim/DESIGN-3D.md)) runs Whitham near field → nonlinear augmented-Burgers rays (US Standard Atmosphere 1976) → 1-D transfer matrix, 2-D cross-plane BEM and along-track Floquet orders at the shield → ground carpet with wind and turbulence → PLdB and peak psf, always against an unshielded run with identical numerics. Anchors: Concorde M2 1.97 psf (NASA 1.94); SR-71 M3 0.99 psf (~0.9); Concorde carpet half-width 42.8 km (~44); secondary-boom delay 8–13 min (10–15).

**2.1 Loudness.** The boom sweeps the shield at aircraft speed V, slower than sound in H₂ (to ~Mach 3.8), so it cannot propagate in the layer and decays as exp(−κd). A 1 km × 30 m H₂ layer cuts **−34 dB undertrack and −25 dB at the carpet maximum** (Concorde M2: ~98 → 58–73 PLdB, conservative proxy). Width sets loudness (250 m: only −8 dB); wind moves the carpet but changes undertrack ΔPLdB ≤ 0.4 dB. Anything that locks the volume (stiff ties, struts, stiff skins, open gaps) acts as an acoustic stud, so every viable structure is soft and quasi-zero-pressure.

**2.2 The 0.11 psf limit** (NPRM 91 FR 40470, 2 July 2026, docket FAA-2026-6935; rule text read only via law-firm summaries, **unverified**). It binds harder than loudness: the 30 m layer leaves an inaudible 0.24–0.25 psf swell at Mach 2. The peak scales as W·d^−(1.2…1.5); one fitted constant matches 699 cases within ±1.8 dB ([psf_sweep](sim/results/psf_sweep.md)).

| Aircraft (shield top 20 km) | Unshielded psf | Min H₂ for 0.11 psf, undertrack |
|---|---|---|
| Overture M1.7 (est. 62 t) | 1.12 | 20 m (lift sets ~30–36 m) |
| Concorde M2.0 | 1.77 | **54 m** |
| Concorde M2.5 | 1.74 | 77 m |
| SR-71 M3.2 | 1.08 | 71 m |

Width: in 3-D (Concorde M2, 60 m), 1 km fails off-track (11.7 Pa at 39 km); **1.5 and 2 km pass everywhere (≤ 4.2 Pa)**. Only this aircraft has been run in 3-D. Altitude: 19–22 km stands, 20 km for Mach 2 ([altitude-trade](altitude/altitude-trade.md)).

**2.3 Superbooms with the shield** ([focus_booms](sim/results/focus_booms.md)). Any acceleration through the cut-off Mach focuses on the ground: unshielded **4.6–5.0 psf** for Overture-class (~40× the limit), 8 psf for Concorde-class. Over 55 m of H₂ it falls to **0.011–0.017 psf** (Concorde 0.031–0.043), because the rays cross the layer < 0.5 km ahead of the aircraft, before they converge. Deceleration forms no caustic; steady turns focus only above 18.7° bank at M1.7. Unshielded values are lower bounds (the NKZ solver reads ~1.7× below measured F-16 foci).

**2.4 On-ramps** ([onramp_shield](sim/results/onramp_shield.md)). **Merge first, then accelerate:** below the day's cut-off (M 1.08–1.23) no ray reaches the ground, so feeders turn, climb and merge without a shield. Each ramp needs only a ~10 km apron (~20 km on the whole spine, +0.6 %). Concorde-class aircraft need the shield flared to 1.7–1.85 km near cut-off (+10 % area); Overture-class needs 1.1–1.25 km. A shield over only the acceleration zone does not work: the steady boom above cut-off is 10–16× the limit. Widths are wavelength-scaled, not 3-D runs.

**2.5 Secondary booms** ([secondary_booms](sim/results/secondary_booms.md)). The shield does not make them worse (worst bins 0.54–1.02× unshielded). But unshaped Mach 2–3 aircraft at 20 km exceed 0.11 psf in secondary booms when winds form a stratospheric duct (0.14–0.37 psf, shield or not; model ~×2 high vs Concorde observations). That is an aircraft and routing problem.

**2.6 Lid fatigue** ([lid-fatigue](fatigue/lid-fatigue.md)). The floating lid is mass-controlled: it carries the reflection-doubled boom (300–650 Pa) by inertia, with strain **~1–6 × 10⁻⁵, 40–200× below cold PE endurance**. Fatigue bites only at bad details (light flaps, slack partition creases, ties on grommets), which weather breaks anyway. With fatigue-aware details (mass-matched flaps, spread tie pads, rolling-lobe partitions, UV top skin) the **cap is ≥ 10,000 flights/day (base; low 1,000, high none)**. **The earlier ~100 flights/day cap is retired.** Spine peak load is 4,500/day (fixed demand) to 6,600/day (base demand + growth).

**Baseline.** Floating-lid H₂ pillow at 20 km, ~55 m at Mach 2 (36 m lift-set for Mach 1.7), 1.5–2 km wide, lid ~3.9 kg/m², soft ties K ≈ 10 Pa/m, 250 m compartments, ~2–4 Mt of H₂ per 4,000 km. Ruled out: narrow modules, mass-law shields, helium, vacuum balloons, drop-stitch/struts/stiff skins, open joints ([3d_summary](sim/results/3d_summary.md), [structure](structure/structure-concept.md), [vacuum](vacuum/vacuum-balloon-review.md)).

## 3. Engineering open items

From [tethers-build-failure](operations/tethers-build-failure.md) unless noted.

- **Tethers.** A bare tether for the psf layer carries **28.4 MN** (21.5 cm, 555 t), and hold-down is **0.62 kg/m², not 0.55**. Seasonal crosswind reversal needs **two anchor rows, ~2,700 tethers** (the one-row ~1,300 figure is retired). Fairings cut tension only ~10 %; **rounding the rim** (C_D 1.0 → 0.3, assumed) brings it to 11.3 MN and 221 t, but needs an acoustic edge check. Self-buoyant tethers are rejected (jet-stream gust loads ~70× bare-rope drag).
- **Construction.** Sealed 250 m × 1.5 km modules inflate on the ground as a 3.5 m mattress and rise in ~1.4 h. Joints close with membrane flaps (~56,000 km of robotic zips), and tethers are raised per 3 km bay. The pace is **8–16 modules per launch-window day for 10 years** (50–150 window days a year, estimated). Film manufacturing (~24 Mt of lid) sets the pace, not launch.
- **H₂ make-up.** 0.1–0.8 Mt/yr (5–20 %/yr). Electrolysis on the lid fails for lack of water (5 ppmv), so the plan is **an all-dielectric H₂ hose on every ~10th tether**: 5–42 TWh/yr of ground electricity.
- **Failure modes.** The top hazard is **a bay rising into the flight level ~100 s after its hold-down is lost**. A bare tether falling from 20 km is lethal (17 kJ/m) and drifts 9–17 km. One lost compartment is survivable, because the bay stays up. Residual risk stays **high** for compartment fire and aircraft–tether strikes. Mitigations: twin ropes, anti-cascade stations, firebreak joints, drogues, automatic traffic stops.
- **Land use and rights.** Two rows take **3–20 ha/km** of point sites (HSR: ~1.4–8 ha/km of continuous strip), with no severance but fall strips up to 24 km downwind. Land costs $0.01–4 M/km, < 1 % of cost ([economics](economics/economics.md) §3). The binding consent is the FAA's: the shield is a moored balloon under 14 CFR 101, so it needs a Part 73 restricted area. An opaque lid could shade farms up to ~380 h/yr (untested takings risk).

## 4. Network and traffic

- **Routes** ([routes](network-maps/routes/routes.md)). Overture's climb, acceleration and descent take ~1,125 km, so a shield buys **nothing below ~1,100 km**. Coast to coast it adds ~51 min over Boomless Cruise, which already saves ~62 min against subsonic. DEN–ORD is the densest shield-km among single routes.
- **Trunk** ([trunk](network-maps/trunk/trunk.md)). A greedy builder finds a single **3,614 km spine, LAS–DEN–OMA–CHI–CLE–NYC**. Feeders reach it at Boomless speed in the same aircraft. Hub transfers (B/C 0.005) and the hand-drawn fans (B/C halved) both lose. Mach barely matters: M1.7 → M3.2 moves B/C only 0.039 → 0.045.
- **Coverage** ([coverage](network-maps/coverage/coverage.md), DOT Table 1 + gravity fill, mirrored data). US domestic air is **~650 M one-way O&D trips/yr** (600–700 M). At Mach 2, **~217 M pax/yr are addressable** (≥ 10 min saved on their own great circle), worth 122 M pax-h/yr against Boomless. The spine carries **81 M pax/yr** (38 % of addressable pax, **~37 % of addressable pax-h**). Its busiest km (OMA–DEN) sees ~80 M pax/yr. The misses are mostly southern: Florida–West and Texas–Northeast get 0 %. Every added arm lowers B/C. The knee lies at 5,400–7,000 km (53 % of pax-h, B/C 0.032), and coverage plateaus at ~61 %. trunk.md's gravity model ran 1.22× high; spine users barely change (79 → 81 M).

## 5. Economics

**Base** ([economics](economics/economics.md)). The M2.0 layer costs **$90/m², $157 M/km, $23.3 M/km/yr** (M1.7: $77.5/m², $17.3 M/km/yr). That is about California HSR per km, and 9× the handoff's $10/m². The pairwise best (DEN–ORD, 10 % capture) scores B/C 1.5 × 10⁻³. *Economics.md's "~100 flights/day fatigue cap" and "free shield still loses" arguments rest on the retired cap.*

**Cost-down** ([cost-down](economics/cost-down/cost-down.md)). Plausible levers cut cost **×2.9** ($36/m², $8.2 M/km/yr); aggressive ones **×7.0** ($20/m², $3.4 M/km/yr net of $0.87 M of PV rent). Materials do most of it; the M2 lid mass can be cheap filled ballast. Break-even against Boomless needs ≤ $1.0 M/km/yr, below the aggressive opex alone. Not in the joint table below: the two-row tether correction (+17–19 % capex) drops the aggressive M2.0 B/C vs subsonic from 1.11 to 0.92 (fixed demand).

**Induced demand** ([demand](network-maps/demand/demand.md)). The shield cuts door-to-door time only 9–16 % against Boomless. Base case: **+4.8 % trips (2.3 M net new a year), B/C uplift ×1.066 vs Boomless and ×1.046 vs subsonic**. High case: ×1.78 / ×1.45, driven by day-return thresholds. The rule of half caps the gain, and the spine's knee does not move. Against subsonic, Boom's ×3 business fares induce nothing, and only ~10 % would switch at ×1.5 fares. **So 100 % capture against subsonic implies fare parity.**

**Joint table: spine, M2.0, B/C** (from [joint_costdown.csv](network-maps/demand/joint_costdown.csv); demand scenario named in each row; cost-down's own `joint_demand.csv` is superseded)

| Cost | Demand scenario | vs Boomless, 10 % / 100 % capture | vs subsonic, 10 % / 100 % capture |
|---|---|---|---|
| Base | fixed | 0.0044 / 0.044 | 0.016 / 0.16 |
| Base | **base** | 0.0047 / 0.047 | 0.017 / 0.17 |
| Base | high | 0.0078 / 0.078 | 0.023 / 0.23 |
| Plausible (÷2.85) | **base** | 0.013 / 0.13 | 0.047 / 0.47 |
| Plausible | high | 0.022 / 0.22 | 0.066 / 0.66 |
| Aggressive (÷7.0) | **base** | 0.033 / 0.33 | 0.12 / **1.16** |
| Aggressive | high | 0.055 / **0.55** | 0.16 / **1.61** |
| Aggressive | base + growth ×1.4 | 0.046 / 0.46 | 0.16 / 1.63 |

M1.7 (Overture-class) is similar: aggressive cost with base demand gives 0.34 vs Boomless and 1.60 vs subsonic at 100 % capture (0.54 / 2.18 with high demand).

**Reading it.** Against Boomless, no row reaches 1; the best is 0.55 (aggressive cost, high demand, 100 % capture). Against subsonic, B/C > 1 needs aggressive cost *and* ~100 % capture (1.16 base demand, 1.61 high demand). At 10 % capture every row stays ≤ 0.16 (aggressive, high demand, vs subsonic).

## 6. Emissions flag

Rough, and not in any B/C above ([demand](network-maps/demand/demand.md) §4; CO₂ 90 g/pax-km subsonic, supersonic 3–7×, $190/t). Base-demand induced trips add **3.6 Mt CO₂/yr**, costing **$0.69 B against $0.11 B** of induced consumer surplus. Moving ~200 M trips from subsonic to supersonic adds **100–300 Mt/yr**, comparable to all of today's US domestic aviation (~130 Mt). Its carbon cost is **1.4–3× the time benefit**. So the only cases above B/C 1 (vs subsonic) would not survive carbon pricing. Non-CO₂ effects at FL600 and stratospheric H₂ leakage are not counted.

## 7. Decision drivers

1. **The regulatory baseline: is Boomless Cruise allowed, and reliable?** If it is, it is the comparator and the shield's best B/C is 0.55; there is no transport case. If it is banned or unavailable on most days, the comparator becomes subsonic (B/C 0.16–0.17 at base cost, base demand, 100 % capture). Today §91.817 still bans overland supersonic flight, and the 0.11 psf NPRM is not final.
2. **Capture and fare parity.** Every 100 % figure assumes the shield's aircraft carries all the traffic on its pairs. Against subsonic, that requires fares near parity. At 10 % capture, B/C falls 10×.
3. **Cost per m².** $90 (base) → $36 (plausible) → $20 (aggressive). Only the aggressive end crosses 1, against subsonic, and only with full capture. The two-row tether correction pulls it back under 1 at fixed demand.

## 8. Next steps (ranked)

| # | Step | Why |
|---|---|---|
| 1 | Establish the regulatory baseline: Boomless status under §91.817 and the NPRM; docket comment on peak vs building-response metrics | Picks the comparator; a 3–4× swing in B/C |
| 2 | Make capture a function of fare (endogenous), incl. Boom's pricing | 100 % capture is the main assumption in every B/C > 0.1 |
| 3 | Supplier quotes (sheet, filled ballast, UHMWPE); 217 K coupons of seam, tie pad, lobe hinge and filled ply to 10⁸ cycles | Cost per m² and the fatigue cap both rest on them |
| 4 | H₂ hazard analysis: bay uplift, compartment fire, tether strike; twin-rope and traffic-stop design | Highest residual risks |
| 5 | 3-D acoustic runs: rounded rim, M1.2 width near cut-off, psf-thickness width for each aircraft, 55 m × 1.5 km floating lid with soft ties | Tether lever, on-ramp flare, width beyond Concorde M2 |
| 6 | Emissions: non-CO₂ at FL600 and H₂ leakage; carbon-priced B/C | Could reverse the only B/C > 1 cases |
| 7 | Demand: DB1B instead of mirrors and gravity fill; winds, merge sequencing, airport capacity (LAX +45 % movements) | Coverage and spine loads |
| 8 | Aero-structural model of a 1.5–2 km surface; daily lift-swing scheme | Gusts, mountain waves, ~10 % lift swing |

## 9. Evidence map

| Claim | File |
|---|---|
| Original concept, surrogates; review and Gate 0 | [README](README.md), [handoff](handoff/2026-10-20-boomshield-handoff.md), [review addendum](handoff/2026-10-20-review-addendum.md) |
| Mechanism, validation, width sets loudness | [gate0](sim/results/gate0_summary.md), [review_evanescent](sim/results/review_evanescent.md), [crossplane](sim/results/crossplane_validation.md), [3d_summary](sim/results/3d_summary.md) |
| Ties, joints, wind and turbulence | [line_ties_check](sim/results/line_ties_check.md), [alongtrack](sim/results/alongtrack_summary.md), [wind_turbulence](sim/results/wind_turbulence.md) |
| 0.11 psf thickness, width, weight budget | [psf_sweep](sim/results/psf_sweep.md) |
| Secondary booms; superbooms; on-ramps | [secondary_booms](sim/results/secondary_booms.md), [focus_booms](sim/results/focus_booms.md), [onramp_shield](sim/results/onramp_shield.md) |
| Altitude; structure; vacuum branch closed | [altitude-trade](altitude/altitude-trade.md), [structure-concept](structure/structure-concept.md), [vacuum review](vacuum/vacuum-balloon-review.md) |
| Lid fatigue; revised cap | [lid-fatigue](fatigue/lid-fatigue.md) |
| Tethers, build, H₂ make-up, failure, land use | [tethers-build-failure](operations/tethers-build-failure.md) |
| Routes, trunk, coverage | [routes](network-maps/routes/routes.md), [trunk](network-maps/trunk/trunk.md), [coverage](network-maps/coverage/coverage.md) |
| Base cost, land rights; cost-down | [economics](economics/economics.md), [cost-down](economics/cost-down/cost-down.md) |
| Induced demand, emissions, **joint B/C** | [demand](network-maps/demand/demand.md), [joint_costdown.csv](network-maps/demand/joint_costdown.csv) |

**Superseded statements.** Retired: the ~100 flights/day lid-fatigue cap (economics.md, trunk.md, onramp_shield.md rule 7); one tether row of ~1,300 tethers at ~26 MN (now two rows, ~2,700, 28.4 MN); hold-down 0.55 kg/m² (now 0.62); cost-down's `joint_demand.csv` and the 0.61 / 1.9 figures in its verdict (preliminary uplifts); the README's "Main result" and the handoff's 10⁵–10⁶-cycle lid limit.
