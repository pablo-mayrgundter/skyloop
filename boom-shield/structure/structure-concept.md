# BoomShield — structural concept study

**Skyloop · 21 October 2026 · screening study.** Setup: a 1 km wide H₂ layer at 20 km (216.65 K, 5529 Pa, ρ 0.0889 kg/m³) under a Concorde-class aircraft at Mach 2 (V 590 m/s), h = 250 m. All values are per m² of planform.

**Files** (`boom-shield/structure/`): `raft_bem.py` (tube-raft BEM on `boomsim.bem2d`); `run_raft_bem.py` and `analyse_bem.py` → `results/raft_bem*.csv`, `figures/raft_vs_slab.png`; `budget.py` → `results/budget_*.csv`.

## Verdict

1. **A locked volume is a stud; superpressure itself is not.** Raft results:
   - Superpressure up to 1 kPa (ovalisation 3Δp/R ≤ 200 Pa/m) is harmless.
   - The skin's hoop-stretch ("breathing") stiffness must stay at **E·t ≲ 1×10⁵ N/m**, i.e. Et/R² ≈ 100–450 Pa/m. At E·t = 3×10⁵–10⁶ N/m (any load-bearing laminate) the shield is nearly transparent.
   - Holding volume against the 255–1,000 Pa diurnal superheat swing needs exactly that stiffness, so **every viable architecture is quasi-zero-pressure and breathes 5–10 % daily**.
   - The vacuum review's rule extends from "anything holding the faces apart" to "anything holding the volume fixed".
2. **Recommended baseline: the floating-lid pillow (architecture d).**
   - H₂ layer 30 m thick.
   - The system mass sits on the top face ("lid" ≈ lift); the floor is a light film.
   - Soft ties K ≈ 10 Pa/m; slack (bellows) partitions.
   - Best acoustics modelled; the soft elements carry almost no static load and the faces about 0.2 kN/m.
3. **Alternate: a touching raft of compliant-hoop tubes**, D ≈ 60 m, Δp 100–300 Pa. It equals the 60 m slab and needs no ties, but its hoop must be strong (3–9 kN/m) and soft (3–9 % strain), which is not standard balloon construction.
4. **Rejected:**
   - (b) drop-stitch at superpressure: ties 12–80× too stiff.
   - (c) lens with cross-width struts: 9–130 kg/m² of structure for 2–4 kg/m² of lift.
   - Rafts with air gaps: a 2 m slit costs 15 dB.
   - Conventional stiff-skin superpressure rafts.
5. **Hydrogen is forced.**
   - A 4,000 km × 1 km × 30 m corridor holds 0.74 Mt of H₂, 0.8 % of one year's world supply.
   - In helium it would be 45 years of world production.

## 1. Acoustics: cross-plane BEM, raft vs slab

**Model.**
- Each tube is its own K₀ gas domain.
- Touching tubes share a flattened double-skin contact strip, b = 1–4 m, with no gas.
- The skin has mass, along-track tension Δp·R/2, and a curved-membrane hoop operator N(w″ + w/R²) − κ(Et/C)∮κw (N = Δp·R) via `Problem.ties`. The second term is hoop stretch with free sliding.
- **Checks:** single tensioned tube vs exact mode series 0.12 %; contour energy balance 0.6 %.
- **Metric:** |p/p_inc| power-averaged over 5–30 Hz (2–3 Hz in brackets) at three places: the 3 km matching circle (φ 0–50°), undertrack (φ ≤ 3°), and a line 100 m below the shield (|y| ≤ 300 m). The near line exposes local leakage that edge diffraction hides at 3 km.
- Lift-closed faces unless noted. Linear incident field.

**Table 1 — dB re unshielded**

| configuration | circle | undertrack | near line |
|---|---|---|---|
| slab H₂ 30 m (reference) | −26.0 (−15.7) | −32.0 | −36.6 |
| slab H₂ 60 m | −24.7 (−21.7) | −35.5 | −45.6 |
| raft D 30, limp (Δp 0, Et 0) | −24.0 (−9.5) | −23.2 | −25.1 |
| raft D 30, Δp 300, Et 10⁴ / **10⁵** / 3×10⁵ / 10⁶ N/m | −24.9 / **−27.0 (−19.9)** / −11.5 / **−4.2** | −24.5 / −26.9 / −10.3 / −3.1 | −26.3 / −29.4 / −11.1 / −3.6 |
| raft D 30, Et 10⁵, Δp 0 / 1,000 / 3,000 | −26.8 / −27.2 / **−10.3** | −26.9 / −26.4 / −23.9 | −29.2 / −29.1 / −11.3 |
| raft D 30, contact strip 4 m | −26.9 | −26.7 | −29.4 |
| raft D 30, air gaps 2 / 5 / 15 m | **−12.3 / −8.4 / −4.2** | −9.7 / −6.3 / −2.7 | −11.5 / −7.5 / −3.3 |
| raft D 30, skin 0.2 kg/m² (not lift-closed) | −25.7 | −23.7 | −26.0 |
| raft D 20, Et 4×10⁴ / 10⁵ | −27.6 / −17.3 | −29.1 / −15.7 | −30.5 / −17.0 |
| **raft D 60, Et 10⁵, Δp 300 / 1,000** | **−27.1 / −27.9** (−23) | **−36.3 / −37.8** | **−41.2 / −41.5** |
| raft D 60, Et 4×10⁵ / 1.6×10⁶ | −24.8 / −11.0 | −23.7 / −9.1 | −27.4 / −10.1 |
| **floating lid T 30 (lid 2.38, floor 0.1 kg/m²), K 0 / 10 / 30** | −28.1 / **−28.4** / −28.7 (−17) | −32.2 / **−33.4** / −36.3 | −36.4 / **−38.2** / −41.7 |
| floating lid, K 100 | −27.5 (−15.7) | −37.0 | −32.1 |
| light lid 0.15 / floor 0.05, K 10 (or K_f 10–100 hold-down to a net) | −27.8 (−26.9) | −24.6 (−23.6) | −29.3 (−28.0) |
| floating lid T 20, K 10 | −26.3 | −25.3 | −29.2 |
| lift-closed slab, smeared ties K 300 | **−18.3** (−2.6) | −20.4 | −18.9 |

**Reading the table.**
- **Edge diffraction still sets the 3 km level.** Every good case lands at −24 to −29 dB.
- **Breathing stiffness:** allowable hoop E·t ≈ 10⁵ N/m, roughly independent of D (4×10⁴ at D 20). Moderate stiffness, or ties of 10–30 Pa/m, gains 3–10 dB at 2–3 Hz by offsetting the gas's negative dynamic stiffness (review §3.1).
- **Contact cusps are harmless and air slits are fatal.** D 60 matches the 60 m slab, while D 30 leaks 6–9 dB more locally.
- **Face mass matters:** a light lid or light skins lose 3–9 dB undertrack. **Spend the excess lift on face mass** (durable laminate, thermal control, photovoltaics), not on tether pull.

![BEM transmission](figures/raft_vs_slab.png)

## 2. Architectures: mass and feasibility

Assumed: HAA-class laminate 0.2 kg/m² per 100 kN/m ultimate at safety factor 4; minimum film 0.03 kg/m²; seams +15 %; controls 0.05 and tether load-spreading 0.10 kg/m²; reserve 15 %; superheat 20 K (hot case 39 K).

Volume factor: a circle or ellipse holds π/4 = 0.785 of the slab volume, so a D-tube raft lifts like a slab of 0.785 D.

**Table 2 — kg/m², H₂ lift density 0.0827 kg/m³**

| architecture | lift | structure | net after reserve | face coupling [Pa/m] | status |
|---|---|---|---|---|---|
| (a) raft D 30 (Δp 660; hot 1,150) | 1.95 | 0.45 (0.66) | 1.21 (1.00) | oval 132–229; breathing Et/R² ≤ 450 | ✔ only with compliant hoop |
| (a) raft D 60 | 3.90 | 0.73 (1.15) | 2.58 (2.16) | oval 66–115; breathing ≤ 110 | ✔ alternate |
| (b) drop-stitch superpressure, T 30 / 60 | 2.48 / 4.96 | 0.39 / 0.63 | 1.72 / 3.59 | **2,420 / 1,210** = Δp/(ε_w T) | ✘ |
| (c) lens T 30–60, zero-pressure / superpressure | 1.95–3.9 | **9.1 / 70–130** | < 0 | ~0 | ✘ |
| (d) floating lid T 30 | 2.46 | 0.28 + lid mass ≈ lift | ≈ reserve | 3–30 (chosen) | ✔ baseline |

**(a) Raft.** Hoop tension holds Δp, so nothing spans the faces. The tubes must touch, laced at both cusps with horizontal links. The hoop needs lobed or fluted skin on along-track tendons with E·t ≲ 10⁵ N/m; a stiff hoop of 3×10⁵–10⁶ N/m fails at D 30 and D 60.

**(b) Drop-stitch.** 0.5 m pillowing needs ties every √(4·T·f) ≈ 8 m. Smeared tie stiffness Δp/(ε_w·T) is independent of spacing; soft ties at K = 30 Pa/m would need 73 % working strain and a 17 m day–night stroke. **Failure confirmed.**

**(c) Lens.** Zero-pressure face tension N = Δρ·g·W²/4 ≈ 200 kN/m, independent of T (2.8–5.5 MN/m with superpressure). The 2N rim compression needs CFRP struts 1 km long, 14–33 m in radius, every 2 km. Narrow lenses side by side are a tube raft.

**(d) Floating lid.**
- **No evanescent gas is neutral:** ρ_g ≤ γp/V² = 0.022 kg/m³ at Mach 2, so it lifts ≥ 75 % of the displaced air. A hung acoustic layer cannot hang; the layer is the lift body.
- **Put the mass on the lid.** The floor becomes a gravitationally stable H₂-over-air interface (+Δρ·g = 0.81 Pa/m), and the rim bulge (≈ 180 N/m) self-equilibrates the faces: **no compression member, no static tie load.**
- The lid is Rayleigh–Taylor unstable for λ > 2π·√(N/Δρg) ≈ 94 m. Ties with K > 0.8 Pa/m stabilise it.
- Partitions must be slack (bellows), or they become line ties.

## 3. Loads

**Diurnal.**
- Constant-volume Δp swing = p·ΔT/T: 255 / 510 / 995 Pa for 10 / 20 / 39 K. 39 K was measured on a He superpressure balloon (abstract, unverified). Reported NASA superpressure-balloon flight Δp is ≤ 180 Pa (stratocat, unverified).
- **Zero-pressure lift does not self-limit with altitude.** Lift = m_gas·(T_g·M_air/(T_a·M_gas) − 1) rises about 10 % for 20 K regardless of height, so the layer climbs until it is volume-limited.
- **Options:** tethers (0 → 2.4 kN/m vertical swing at T 30, ≥ 10 % night reserve); a solar gas shuttle pumping 9 % of the gas into rim reservoir tubes at ~5 kPa (~0.7 W/m² over 4 h, ~0.2 kg/m², outside the acoustic span); or coatings holding superheat under 10 K (unverified).
- With soft ties (K 10), daily thickening of 2.7 m puts 0–27 Pa through the ties. K 3 keeps that to 8 Pa.

**Wind** (20 km mean 10–20 m/s, winter 20–30 m/s, peaks about 50 m/s; airship literature, abstracts):

| U [m/s] | side force [kN/m of corridor] | quasi-steady gust lift, w′ = 5 % U [Pa] |
|---|---|---|
| 20 | 0.50 | 5.6 |
| 30 | 1.12 | 12.6 |
| 40 | 1.99 | 22.3 |

The gust lift on the 1 km chord is an upper bound: added mass ρπW/4 ≈ 70 kg/m² makes the layer ride gusts. Mountain waves and flutter remain open.

**Station-keeping.**
- **Propulsion is excluded:** F·U/η ≈ 67 kW per metre of corridor (~270 GW for 4,000 km).
- **Tethers** every 2–5 km at 30 m/s carry 2.2–5.6 MN each: UHMWPE 6–10 cm, 90–220 t, 0.16 kg/m² of downforce, plus 0.6–1.0 MN of jet-stream drag on the tether (addendum, scaled).

**Boom transient.**
- With |R| ≈ 1 the top face sees 200–600 Pa (10–25× the lid's lift head, 1–6× the raft Δp): skins go slack every pass (0.1–0.3 m, review §3.7) and soft ties need damping.
- 20–100 flights a day for 30 years is 2×10⁵–10⁶ cycles. **Flex-crack pinholes are the main leak driver.**

## 4. Hydrogen

| quantity | value |
|---|---|
| inventory, 4,000 km × 1 km × 30 m slab (raft D 30) | 0.74 (0.58) Mt; 8.3 (6.5) × 10⁹ Nm³ |
| share of world supply (~97 Mt in 2023, IEA) | 0.8 % (0.6 %) of one year |
| as helium (~185 × 10⁶ m³/yr, 2024, secondary, unverified) | **45 (35) years** of world output |
| permeation, 25 µm LDPE (7 Barrer at 23 °C, E_a 35 kJ/mol, both assumed) | 36 %/yr warm; ~0.2 %/yr at 217 K |
| replenishment (holes and flex cracks dominate) | 5–20 %/yr = 0.04–0.15 Mt/yr; first fill $1–4 B at $2–6/kg |
| stored combustion energy | 60–90 PJ, ≈ 5 kt TNT-equivalent per km of corridor |

**Safety flag.**
- The H₂–air lower flammability limit stays near 4 % down to about 50 mbar.
- Reported minimum ignition pressures: below 50 mbar with a hot wire, about 100 mbar with a spark (abstracts, unverified). At 55 mbar the layer is marginal for sparks and not safe against engine exhaust.
- Leaked H₂ rises toward the flight level; debris falling into the layer is a fire-spread case.
- **Needs:** compartments every ≤ 250 m, leak detection that closes sections to traffic, and a hazard analysis before hardware.

## 5. Next iteration: parameters to test

1. **Floating lid (baseline):** discrete soft ties at 5 / 10 / 25 m pitch with K 3–30 Pa/m (discreteness check; needs the `bem2d` line-tie fix another agent is making); bellows partitions every 100 m; lid 1.2–2.4 kg/m²; T 20 / 30 / 40 m; **full 3-D PLdB carpet** with and without a 150 m fringe; nonlinear slack lid under a 200–600 Pa pulse.
2. **Raft (alternate):** D 40 / 60 / 80 m; E·t 3×10⁴–2×10⁵ N/m; Δp 100–1,000 Pa; contact gap 0 / 0.3 / 1 m; cusp lacing; end domes every 250–500 m; 3-D carpet.
3. **Both:** ±9 % thickness breathing, tether point loads, wind-driven roll.

## 6. Open risks

1. **Diurnal lift swing.** Tethers, a gas shuttle, or coatings must absorb about 10 % of lift every day.
2. **Soft elements.** Tie springs (floating lid) or strong-and-soft hoops (raft) must survive 10⁵–10⁶ slack/taut cycles at 217 K.
3. **H₂ fire next to traffic.**
4. **1 km-chord aerodynamics:** gusts, mountain waves, roll, flutter of slack faces.
5. **Leakage growth** from film fatigue.
6. **Tether field:** 800–2,000 tethers through the jet stream.
7. **Model limits:** 2-D, x-invariant and linear; ratios at the matching circle, not ground PLdB; idealised raft contact; static RT estimate.

**Sources** (search-result level only; page fetches were blocked):
- IEA *Global Hydrogen Review 2024*: <https://www.iea.org/reports/global-hydrogen-review-2024/hydrogen-demand>
- USGS helium (world total via secondary press, unverified): <https://pubs.usgs.gov/periodicals/mcs2026/mcs2026-helium.pdf>
- NASA superpressure balloon:
  - <https://stratocat.com.ar/fichas/2022/WNK-20220527.htm>
  - <https://aerostar.com/news/raven-s-aerostar-subsidiary-announces-successful-flight-of-nasa-prototype-super-pressure-balloon-in>
- Superheat 39 K and 649 Pa (abstract): <https://onlinelibrary.wiley.com/doi/abs/10.1155/2018/7806036>
- Stratospheric wind minimum: <https://arxiv.org/pdf/1402.6706>, <https://ams.confex.com/ams/pdfpapers/128256.pdf>
- H₂ flammability at reduced pressure: <https://arxiv.org/pdf/2602.10279>; Le, Nayak & Mannan, AIChE 2012: <https://www.aiche.org/conferences/videos/conference-presentations/flammability-limits-hydrogen-and-its-mixtures-hydrocarbons-air-atmospheric-and-sub-atmospheric>
- SPICE tethered balloon, Davidson et al., *Phil. Trans. R. Soc. A* 370 (2012): cited from the addendum, not re-verified.
- Materials values marked "assumed" are handbook-class guesses.
