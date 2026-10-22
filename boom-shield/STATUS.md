# BoomShield — concept status (end of phase 1)

**Skyloop · 22 October 2026 · final status for this phase, after the acoustics, 0.11 psf, superboom, on-ramp, fatigue, operations, network, economics, induced-demand, integration, higher-Mach, growth, global and net-reconciliation rounds.**
Everything here is screening level: simulation only (no hardware, no flight test), linear shield physics, unshaped sources, and web figures taken from search-result summaries that were not verified at the source. Unit costs marked [U] in the studies are engineering estimates.

**How B/C is quoted.** Every benefit/cost ratio below is labelled **gross** or **net**. *Gross* = value of time saved ($86/h × 95 % availability) over the annualised cost of the shield. *Net* = the same, minus the extra fuel ($0.80/kg) and CO₂ ($190/t) of whole trips against the comparator, induced trips at full fuel. One aircraft per Mach sets shield, ramps and fuel; integrated cost is re-run per shield. Canonical values: [net-reconciliation](economics/net/net-reconciliation.md). Unless a row says otherwise: single 2025 year, spine LAS–DEN–OMA–CHI–CLE–NYC, two tether rows, base induced demand, 100 % capture (B/C scales linearly with capture).

## 1. Bottom line

- **The acoustics pass a screen, and nothing found since undoes them.** At 20 km, a hydrogen layer about **55 m thick and 1.5–2 km wide** meets the FAA's proposed **0.11 psf** limit under a Concorde-class aircraft at Mach 2; a ≤ 50 t M2.0–2.2 aircraft needs only 29–35 m × ~1.25 km. Superbooms vanish over it, on-ramps need almost no extra shield, and lid fatigue is a detail-design requirement.
- **Buildability is open, not refuted.** Sealed 250 m modules can be launched and joined; two anchor rows (~2,700 tethers) hold them; the worst failure is a bay rising into the flight level ~100 s after losing its hold-down.
- **Cost can reach ~$14/m² on paper, by owning the plants.** Vertical integration gives **$14.3/m²** ($10.7–21.4), a **$91 B** spine and **$7 B** of plants, against $106/m² bought at base prices.
- **Gross, the shield can reach B/C ≈ 1 against Boomless Cruise, but only with every lever pulled.** Integrated cost alone gives **0.49 gross** (Concorde-class shield); add a light M2.0–2.2 aircraft, **0.79–0.87**; M3, 0.98; M2.2 with high induced demand, 1.47; a ~2050 opening, 1.06 (1.48 with the light M2.0 jet). All at 100 % capture.
- **Net of fuel and CO₂ it does not pay.** In a 2025 year, net B/C against Boomless is **0.28** with the light M2.0 jet, 0.04 at M2.2 and 0.41 at M1.7 (integrated, 100 % capture). Life-cycle, it stays **≤ 0.75** for openings to 2050 (M2.0; 0.85 at M1.7) and reaches 1 only for a ~2060–65 opening (2075 at M2.2; base world). Against subsonic, gross B/C is 1.75–2.7, but **net is negative at every Mach, in every corridor and for every opening this century**.
- **Mach is the weakest lever; the sweet spot is M2.0–2.2** (M2.0 net of fuel). Ground time dominates the door-to-door trip, and above ~M2.3 the extra fuel costs more than the time it buys.
- **The US spine is the best corridor in the world and the US the natural first mover.** Europe–Gulf reaches 55–65 % of it.
- **Timing comes from events, not traffic.** Traffic growth lifts B/C ~3.5 % a year. Re-run the study ~2029–31 (FAA final rule, Boomless status, Overture's real economics); a pilot gate falls ~2038 and the earliest defensible spine opening is ~2055.

## 2. What passes: acoustics

**Evidence chain.** `boomsim` ([DESIGN](sim/DESIGN.md), [DESIGN-3D](sim/DESIGN-3D.md)): Whitham near field → nonlinear rays → 1-D/2-D/along-track shield models → ground carpet with wind and turbulence, always against an identical unshielded run. Anchors: Concorde M2 1.97 psf (NASA 1.94); SR-71 M3 0.99 psf (~0.9).

**Loudness.** The boom sweeps the shield slower than sound in H₂, so it decays across the layer. A 1 km × 30 m H₂ layer cuts **−34 dB undertrack and −25 dB at the carpet maximum** (Concorde M2: ~98 → 58–73 PLdB, conservative proxy). Width sets loudness (250 m: only −8 dB). Anything stiff or open leaks the boom, so every viable structure is soft.

**The 0.11 psf limit** (NPRM 91 FR 40470, 2 July 2026, docket FAA-2026-6935; rule text read only via law-firm summaries, **unverified**) binds harder than loudness ([psf_sweep](sim/results/psf_sweep.md)):

| Aircraft (shield top 20 km) | Unshielded psf | Min H₂ for 0.11 psf, undertrack |
|---|---|---|
| Overture M1.7 (est. 62 t) | 1.12 | 20 m (lift sets ~30–36 m) |
| Concorde M2.0 | 1.77 | **54 m** |
| Concorde M2.5 | 1.74 | 77 m |
| SR-71 M3.2 | 1.08 | 71 m |
| Light M2.2, 50 t (highmach, τ-fit) | — | 35 m |

In 3-D (Concorde M2, 60 m) **1.5 and 2 km pass everywhere**; 1 km fails off-track. Only this aircraft has been run in 3-D.

**Superbooms** ([focus_booms](sim/results/focus_booms.md)): unshielded 4.6–5.0 psf (Overture class), **0.011–0.017 psf** over 55 m of H₂. **On-ramps** ([onramp_shield](sim/results/onramp_shield.md)): merge below cut-off, then accelerate; ~10 km aprons, +0.6 % shield. **Secondary booms** ([secondary_booms](sim/results/secondary_booms.md)): not worsened, but ducted returns from unshaped M2–3 aircraft exceed 0.11 psf with or without a shield. **Lid fatigue** ([lid-fatigue](fatigue/lid-fatigue.md)): strain 40–200× below endurance; with fatigue-aware details the cap is **≥ 10,000 flights/day** (low 1,000). Growth brings it back: at 100 % capture the spine reaches 10,000/day by ~2060.

**Baseline.** Floating-lid H₂ pillow at 20 km, 1.5–2 km wide, soft ties, 250 m compartments, ~2–4 Mt of H₂ per 4,000 km ([structure](structure/structure-concept.md), [altitude-trade](altitude/altitude-trade.md)). Ruled out: narrow modules, mass-law shields, helium, vacuum balloons ([vacuum](vacuum/vacuum-balloon-review.md)), struts, stiff skins, open joints.

## 3. Engineering open items

From [tethers-build-failure](operations/tethers-build-failure.md) unless noted.

- **Tethers.** 28.4 MN per bare tether (21.5 cm, 555 t); hold-down 0.62 kg/m²; **two anchor rows, ~2,700 tethers**. A rounded rim (assumed C_D 0.3) brings it to 11.3 MN and 221 t, but needs an acoustic edge check.
- **Construction.** 250 m × 1.5 km modules inflate on the ground and rise in ~1.4 h; 8–16 per launch-window day for 10 years. Film supply sets the pace.
- **H₂ make-up.** 0.1–0.8 Mt/yr via a dielectric hose on every ~10th tether.
- **Failure modes.** Top hazard: bay uplift into the flight level ~100 s after hold-down loss. Residual risk stays **high** for compartment fire and aircraft–tether strikes.
- **Land and consent.** 3–20 ha/km of point sites; land < 1 % of cost. The binding consent is the FAA's (moored balloon under 14 CFR 101, Part 73 restricted area).
- **Supply.** The spine's UHMWPE (116 kt/yr) is about **the whole world's fibre capacity**, and China has controlled its fibre exports since July 2024 ([integration](economics/integration/integration.md)).

## 4. Network and traffic

- **Routes and trunk** ([routes](network-maps/routes/routes.md), [trunk](network-maps/trunk/trunk.md)). Overture's climb, acceleration and descent take ~1,125 km, so a shield buys nothing below ~1,100 km. A greedy builder finds one **3,614 km spine, LAS–DEN–OMA–CHI–CLE–NYC**, with feeders at Boomless speed in the same aircraft.
- **Coverage** ([coverage](network-maps/coverage/coverage.md); DOT Table 1 + gravity fill, mirrored data). ~650 M one-way US domestic O&D trips/yr; **~217 M addressable** at Mach 2; the spine carries **81 M** (~37 % of addressable pax-h). Every added arm lowers B/C.
- **Induced demand** ([demand](network-maps/demand/demand.md)). Base: +4.8 % trips, B/C ×1.066 vs Boomless; high: ×1.78. Against subsonic, Boom's ×3 business fares induce nothing and only ~10 % switch at ×1.5 fares, **so 100 % capture against subsonic implies fare parity.**

## 5. Economics: base, cost-down, demand

- **Base** ([economics](economics/economics.md)): $106/m² with two tether rows, $27.3 M/km/yr net. Gross B/C 0.040 vs Boomless, 0.14 vs subsonic.
- **Cost-down** ([cost-down](economics/cost-down/cost-down.md)): plausible ×2.5 ($42/m²), aggressive ×4.6 ($23/m², two rows; $20 with one). Materials do most of it.
- **The second tether row** (+17–19 % capex) pulls the aggressive case below 1 against subsonic: **0.97 gross** (1.16 with one row, as published).
- **Against Boomless,** break-even needs ≤ $1.10 M/km/yr net cost (gross, base demand); even the integrated low-cost case ($1.48) misses.

## 6. Vertical integration

([integration](economics/integration/integration.md)) Replaces bought materials with owned plants plus Wright's-law learning; every row carries the second tether row.

**Make:** filled sheet and blown film, the module factory, UHMWPE fibre and rope (there is no market to buy 116 kt/yr from), launch bases, station kits. **Buy on long contracts:** PE resin (ethane-indexed offtake; an owned cracker has no edge), H₂ (< 4 % of capex), fluoropolymer, electronics, civil works, CaCO₃.

| Spine M2.0 · cost case | $/m² | Net $M/km/yr | Spine $B | Gross B/C vs Boomless | Gross B/C vs subsonic |
|---|---|---|---|---|---|
| Base | 106 | 27.3 | 672 | 0.040 | 0.14 |
| Cost-down plausible | 42 | 9.6 | 268 | 0.115 | 0.41 |
| Cost-down aggressive | 23 | 4.0 | 147 | 0.27 | 0.97 |
| Integrated, plausible levers | 17.6 (13–26) | 4.1 | 111 | 0.27 | 0.94 |
| **Integrated, aggressive levers** | **14.3** (10.7–21.4) | **2.2** | **91** | **0.49** | **1.75** |
| Integrated, low end | 10.7 | 1.5 | — | 0.74 | 2.63 |
| M1.7 integrated, aggressive | 11.8 | 1.5 | 64 | 0.51 | 2.40 |

- **Plants:** **$7.1 B** ($4.6–12.5 B), ~8 % of spine capex, saving ~$3.5 B a year during the build (payback ~1.6 yr). Of the $8.8/m² saved against bought-aggressive, ~$6.2 is make-vs-buy and ~$2.6 learning.
- **Net of fuel and CO₂:** the integrated row is **0.12** vs Boomless and −2.8 vs subsonic, but it pairs a 147 t shield with a 49 t jet's fuel. With the 49 t M2.0 jet throughout: gross 0.79 / 2.67, net **0.28** / −4.0 ([net-reconciliation](economics/net/net-reconciliation.md)).
- **Risks:** seven first-of-a-kind plants at once; common-mode defects across 14,000 modules; stranded plants after the build (salvage 5–30 %); open-access rules if the owner is tied to an airline.

## 7. Higher Mach

([highmach](network-maps/highmach/highmach.md)) **Why the gain is small: ground time dominates.** At M2, NYC–LAX is 5.2 h door to door, of which only 1.3 h is cruise; 3.6 h is fixed (130 min ground, 88 min taxi, climb, descent, approach). Infinite Mach gives 4.5 h; Boomless 6.3 h. Ramps set the asymptote, and a 70 min walk-on overhead saves 60 min, more than M2 → M∞, lifting B/C ×1.3–1.4 at every Mach.

**Table: spine, base demand, 100 % capture** (two rows; integrated cost re-run per shield; canonical net, [net-reconciliation](economics/net/net-reconciliation.md))

| Design | NYC–LAX h | H₂ m | Gross vs BL: aggressive / integrated | Net vs BL: aggressive / integrated | Vs subsonic, integrated: gross / net | Break-even VOT net of fuel |
|---|---|---|---|---|---|---|
| Overture M1.7 | 5.46 | 25 | 0.28 / 0.50 | 0.23 / 0.41 | 2.35 / −3.52 | $10/h |
| Concorde-class M2.0 (147 t) | 5.20 | 55 | 0.27 / 0.49 | 0.07 / 0.12 | 1.75 / −2.82 | $39/h |
| **Light M2.0, adaptive (49 t)** | 5.20 | 29 | 0.45 / **0.79** | 0.16 / **0.28** | 2.67 / −4.01 | $34/h |
| **A: M2.2 Al-limit (50 t)** | 5.08 | 35 | 0.50 / **0.87** | 0.02 / **0.04** | 2.69 / −4.41 | $50/h |
| M2.5 Ti | 4.94 | 47 | 0.53 / 0.91 | −0.53 / −0.91 | 2.60 / −5.90 | $105/h |
| **B: M3.0 Ti** | 4.80 | 73 | 0.57 / **0.98** | −1.37 / −2.37 | 2.57 / −8.16 | $177/h |
| C: M3.5 TBCC | 4.68 | 128 | 0.56 / 0.97 | −3.85 / −6.64 | 2.34 / −14.5 | $414/h |

- **Most of the gain is weight and thrust, not Mach:** a 50 t aircraft cuts the layer from 55 m to 29–35 m, and adaptive-cycle thrust doubles ramp acceleration. M2.0 → M2.2 adds ~10 %; M4 adds 24 %.
- **Fuel undoes speed above ~M2.3.** Fuel per seat-km is ×1.16 (M2.2) and ×1.74 (M3) an Overture's; the study uses $86/h of time value. **Net vs subsonic is negative at every Mach.**
- **Engine reality.** The F-47/NGAD-class adaptive-cycle engines (AETP/NGAP) claim +25 % efficiency, but subsonically; in supercruise the cycle is turbojet-like (~3 % SFC gain). The civil gain is thrust margin. The F-47 flies on interim engines in 2028, NGAP integrates ~2030, and export controls push a civil derivative to the **late 2030s–2040s**.
- **Shield ceiling.** V = c_H₂ gives **M3.80 for pure H₂ but M3.37 at 98 % purity** (2.94 at 95 %). Warming the gas to M4 needs ~1,900 TWh/yr on the spine; solar superheat works only by day. Aluminium airframes stop at ~M2.2 (titanium from ~M2.4).
- **Sweet spot: M2.0–2.2** (M2.0 best net of fuel), ≤ 50 t, aluminium, adaptive-cycle-class engines, 29–35 m × ~1.25 km of H₂ at 20 km, fast ground processing.

## 8. Growth and timing

([growth](economics/growth/growth.md)) **Traffic.** FAA and OEM forecasts give 2.4–2.8 %/yr for US domestic; history since 2000 is ~1.3 %/yr enplanements, ~2.4 %/yr RPMs, with each shock costing 3–6 years. US trips > 1,500 km (base) reach **1.5× ~2043, 2× ~2057, 3× ~2086, 10× ~2187** (high case ~2125).

**Life-cycle B/C by opening year** (spine, integrated cost, 100 % capture, base world, central energy; 7 %, 25 yr, 5-yr ramp-up; canonical, [net-reconciliation](economics/net/net-reconciliation.md)):

| Opening · aircraft | Gross vs BL | Net vs BL | Gross vs subsonic | Net vs subsonic |
|---|---|---|---|---|
| 2035 · M2.0, 49 t | 0.99 | 0.40 | 3.34 | −4.27 |
| 2045 · M2.0, 49 t | 1.34 | 0.63 | 4.53 | −4.64 |
| 2050 · M2.0, 49 t | **1.48** | **0.75** | 5.00 | −4.48 |
| 2050 · M2.2, 50 t | 1.67 | 0.46 | 5.17 | −5.18 |
| 2050 · Concorde-class shield (growth.md) | 1.06 | 0.46 | 3.76 | −3.66 |

- Delay lifts B/C ~3.5 %/yr (the light jets exist only from the late 2030s). Net vs Boomless reaches 1 only for a **~2065 opening** (M2.0; 2075 at M2.2 or with the Concorde-class shield). At 30 % capture net reaches 1 only ~2100; base cost never passes.
- **Capture is capped by airports:** at 100 % capture LAX would see ~2,000 supersonic movements a day by 2050, more than all its movements today. Realistic capture is nearer 30 % unless the aircraft grows.
- **Triggers:** T1 FAA final rule and Boomless status (2027–28); T2 Overture's real fares and capture (2030–32); T3 a ≤ 50 t M2.0–2.2 airliner with a civil adaptive-cycle engine (2035–40); T4 film and fibre at scale, ≤ $20/m² (2028–31); T5 a km-scale structure held ≥ 1 yr at ≥ 18 km; T6 H₂ safety precedent; T7 SAF ≤ $1.5/kg or a stable carbon price (sets the sign of net B/C vs subsonic); T8 long-haul ≥ 1.3× 2025 (~2037).

**Staged plan with gates** ([U] costs)

| Stage | Dates | Cost | Gate |
|---|---|---|---|
| 0. Desk studies (§12 steps); full re-run 2029–31 | 2026–29 | $1–3 M | **G-A ~2029:** T1 resolved; if Boomless is legal and reliable, pause to ~2035 unless T4 is strong |
| 1. Lab acoustics (shock tube / ballistic range), fatigue coupons, supplier pilot lines | 2029–33 | $10–40 M | **G-B ~2032:** measured ΔdB/psf within 3 dB of boomsim; T4 met; T2 shows fares and capture |
| 2. Subscale tethered module, a year at 18–20 km, supersonic fly-under | 2033–38 | $0.3–1 B | **G-C ~2038:** ≥ 95 % availability; uplift and fire tests; T3, T5 met |
| 3. On-ramp pilot, 100–200 km at Denver | 2039–45 | $3–8 B | **G-D ~2045:** cost ≤ $20/m², net B/C ≥ 1 on re-run → spine decision, **opening ~2055** |

## 9. Global corridors

([global](network-maps/global/global.md), rescaled to the canonical convention in [net-reconciliation](economics/net/net-reconciliation.md); M2.2 aircraft unless noted, integrated cost, 100 % capture; non-US demand is a gravity model, ±2×)

| Corridor | Overland km | Gross vs BL | Net vs BL: M2.2 / M2.0 | Gross vs subsonic |
|---|---|---|---|---|
| **US north spine** | 3,614 | **0.87** | 0.04 / **0.28** | 2.7 |
| US south LAX–PHX–DFW–ATL–WAS | 4,043 | 0.74 | −0.05 / 0.17 | 2.4 |
| Europe–Levant–Gulf (greedy) | 4,939 | 0.55 | −0.02 / – | 2.2 |
| Europe–Gulf LON–FRA–IST–DXB | 5,299 | 0.46 | 0.00 / 0.14 | 2.0 |
| LON–IST–DXB–BOM–HKG | 9,708 | 0.36 | 0.00 / 0.10 | 1.7 |
| Russia; China X | 2,838; 3,851 | 0.17 | < 0 / ≤ 0.02 | 1.1; 1.8 |
| India, Brazil, Africa, Australia, Andes | — | ≤ 0.05 | ≤ 0 / ≤ 0.01 | 0.1–1.2 |

- **No corridor matches the US spine;** Europe–Gulf is ~55–65 % of it. LON–IST–DXB–BOM–HKG saves the most time (65 M pax-h) but needs 9,700 km of shield and many trips exceed a 5,000 km range (0.19 capped).
- **China is big but short** (19 M trips > 2,500 km vs 109 M in the US). In the tropics Boomless works on 84–95 % of days. **Net vs subsonic is negative in every corridor.**
- **First mover: the US,** the only jurisdiction with a draft numeric boom limit a shield can be engineered to meet. Over nearly empty land, a regulator could simply allow unshielded booms, which competes with a shield on the Asian corridors.

## 10. Emissions

Net B/C now carries fuel and CO₂ (§§6–9). Extra CO₂ on the spine is 6 Mt/yr vs Boomless at M2.0 (11 at M2.2, 53 at M3); against subsonic 84–400 Mt/yr, comparable to all US domestic aviation (~130 Mt). **Non-CO₂ effects** (NOx at 22–26 km, stratospheric water vapour) **and H₂ leakage are not counted** and would lower every net figure.

## 11. Decision drivers

1. **The regulatory baseline.** If Boomless Cruise is legal and reliable overland, it is the comparator and the case rests on **net B/C vs Boomless: 0.28 in a 2025 year (light M2.0 jet), ≤ 0.75 life-cycle for openings to 2050, ~1 only by ~2065**. If not, the comparator is subsonic: gross 1.75–2.7, but net negative. Today §91.817 still bans overland supersonic flight; the 0.11 psf NPRM is not final (final rule targeted mid-2027).
2. **Capture and fares.** Every headline is 100 % capture; 30 % is realistic given airports and fares. Against subsonic, full capture needs fares near parity.
3. **Cost per m².** $106 → $42 → $23 bought → **$14 integrated**. Only integration gets gross vs Boomless near 0.5 with a Concorde-class aircraft; supplier quotes and coupons must confirm it.
4. **A light M2.0–2.2 aircraft.** ≤ 50 t with adaptive-cycle-class thrust nearly doubles gross B/C (0.49 → 0.79–0.87), more than doubles net at M2.0 (0.12 → 0.28) and thins the shield; it does not exist before the late 2030s.
5. **Ground processing.** A 70 min walk-on overhead adds 30–40 % at every Mach, more than any Mach step, and enables day returns.
6. **Fuel and CO₂ price.** It sets the gap between gross and net. With cheap fuel and no carbon price, net vs Boomless rises to 0.66 (2045 opening, Concorde-class shield); stringent prices drop it to 0.11, and keep net vs subsonic negative.

## 12. Next steps (ranked)

| # | Step | Why |
|---|---|---|
| 1 | Regulatory baseline: Boomless status under §91.817 and the final rule; docket comment on peak vs building-response metrics | Picks the comparator (T1, gate G-A) |
| 2 | Endogenous capture: fares by Mach (esp. M2.2), Boom's pricing, airport capacity | 100 % capture underlies every B/C near 1 |
| 3 | Supplier quotes (filled sheet, film, UHMWPE rope) and 217 K coupons to 10⁸ cycles | Tests the $14/m² integrated cost and the fatigue cap (T4) |
| 4 | ~~Reconcile the net-B/C conventions~~ **done** ([net-reconciliation](economics/net/net-reconciliation.md)); next, add non-CO₂ effects and H₂ leakage | Net B/C now decides the case |
| 5 | H₂ hazard analysis: bay uplift, compartment fire, tether strike | Highest residual risks |
| 6 | 3-D acoustic runs: a 50 t M2.2 aircraft (35 m × 1.25 km), rounded rim, on-ramp width near cut-off | The sweet-spot design and the tether lever |
| 7 | Propulsion/airframe sizing check of the ≤ 50 t M2.2 concept | Half the higher-Mach gain rests on it |
| 8 | Demand: DB1B instead of mirrors; MIDT/Sabre for Europe–Gulf | Spine loads; the only non-US corridor worth more study |
| 9 | Aero-structural model of a 1.5–2 km surface; lift-swing scheme | Gusts, mountain waves |
| 10 | Full re-run in 2029–31 when T1 and T2 resolve | Timing comes from events |

## 13. Evidence map

| Claim | File |
|---|---|
| Original concept, review, Gate 0 | [README](README.md), [handoff](handoff/2026-10-20-boomshield-handoff.md), [review addendum](handoff/2026-10-20-review-addendum.md) |
| Mechanism, validation, width | [gate0](sim/results/gate0_summary.md), [review_evanescent](sim/results/review_evanescent.md), [crossplane](sim/results/crossplane_validation.md), [3d_summary](sim/results/3d_summary.md) |
| Ties, joints, wind and turbulence | [line_ties_check](sim/results/line_ties_check.md), [alongtrack](sim/results/alongtrack_summary.md), [wind_turbulence](sim/results/wind_turbulence.md) |
| 0.11 psf thickness and width | [psf_sweep](sim/results/psf_sweep.md) |
| Secondary booms; superbooms; on-ramps | [secondary_booms](sim/results/secondary_booms.md), [focus_booms](sim/results/focus_booms.md), [onramp_shield](sim/results/onramp_shield.md) |
| Altitude; structure; vacuum closed | [altitude-trade](altitude/altitude-trade.md), [structure-concept](structure/structure-concept.md), [vacuum review](vacuum/vacuum-balloon-review.md) |
| Lid fatigue | [lid-fatigue](fatigue/lid-fatigue.md) |
| Tethers, build, H₂, failure, land | [tethers-build-failure](operations/tethers-build-failure.md) |
| Routes, trunk, coverage; induced demand | [routes](network-maps/routes/routes.md), [trunk](network-maps/trunk/trunk.md), [coverage](network-maps/coverage/coverage.md), [demand](network-maps/demand/demand.md) |
| Base cost, cost-down | [economics](economics/economics.md), [cost-down](economics/cost-down/cost-down.md) |
| Make/buy, plants, learning, **two-row B/C** | [integration](economics/integration/integration.md) (`bcr.csv`, `results.csv`) |
| Mach sweep, engines, ceiling | [highmach](network-maps/highmach/highmach.md) (`summary_by_mach.csv`) |
| **Canonical gross and net B/C** (all Mach, cost cases, capture, life-cycle, corridors) | [net-reconciliation](economics/net/net-reconciliation.md) (`net_bc.py`, `canonical.csv`) |
| Traffic growth, opening year, triggers, staged plan | [growth](economics/growth/growth.md) (`roi_headline.csv`, `crossing_years.csv`) |
| World corridor ranking; first mover | [global](network-maps/global/global.md) (`ranking.csv`) |

**Superseded statements.**
- The ~100 flights/day fatigue cap (economics.md, trunk.md, onramp_shield.md); one tether row (~1,300 at ~26 MN); hold-down 0.55 kg/m²; the README's original "Main result"; the handoff's 10⁵–10⁶-cycle lid limit.
- **One-row B/C figures** in demand/joint_costdown.csv and the previous STATUS (1.16 / 1.61 vs subsonic, best 0.55 vs Boomless): use integration's `bcr.csv`. Cost-down's `joint_demand.csv` was already superseded.
- **"No case reaches 1 against Boomless":** true net of fuel and CO₂; gross, ~1–1.5 is reachable with integrated cost, a light fast aircraft and 100 % capture.
- **"Emissions are not in any B/C":** net B/C is now computed, canonically in [net-reconciliation](economics/net/net-reconciliation.md).
- **Earlier net B/C** (highmach "×1.8", global 0.12 / 0.40 without induced-trip fuel, growth's 147 t-shield rows as the M2.0 case) and **ratio-scaled integrated gross** (0.90 → 0.87, 1.03 → 0.98, 1.53 → 1.47): use net-reconciliation.
- **trunk.md's "Mach barely matters" (0.039 → 0.045):** with light adaptive-cycle aircraft gross B/C doubles M1.7 → M3 (highmach), though Mach stays the weakest lever.
- **integration.md's "highmach hook: no B/C table yet"** and **growth.md §4's "global study not yet available":** both studies now exist; global names the US as first mover, not a fast-growing foreign market.
- **"All-in cost must reach ~$20/m²":** integration reaches $14.3/m² on paper.
