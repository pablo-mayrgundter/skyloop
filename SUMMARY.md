# Three concepts compared: BoomShield, Skyloop, Hyperloop (our version)

**Skyloop · 23 October 2026 · screening-level desk studies.** Each study maps the rough feasible operating regime across disciplines (acoustics or tube physics, structure, network, energy, cost, benefit/cost), not a validated design point. Simulation only, no hardware; web figures are unverified search summaries; unit costs are engineering estimates. Every B/C below assumes **100 % capture** unless stated; realistic capture is nearer 30 %.

All three use BoomShield's demand and value-of-time model (US spine LAS–DEN–OMA–CHI–CLE–NYC, tubes adding LAX; $86/h; 7 %, 25 yr), so their B/C can be compared once the basis is labelled. **Gross** counts time saved; **net** also counts the fuel and CO₂ against the comparator (extra for aircraft, avoided for electric tubes).

Infographics: [BoomShield](https://claude.ai/artifact/1GASTWSYperFtVLyLY29VJ) · [Skyloop](https://claude.ai/artifact/Rvmp18jJLdP7pAPWvuiHHK) · [Hyperloop (our version)](https://claude.ai/artifact/QbRskhNLJEzpye6SzGC3fB)

## Where each stands

**BoomShield** (tethered H₂ layer at 20 km that absorbs sonic booms; [STATUS](boom-shield/STATUS.md))
- **Passes:** acoustics pass a screen. A ≤ 50 t M2.0–2.2 aircraft needs **29–35 m × ~1.25 km** of H₂ to meet the proposed **0.11 psf** limit. Buildability is open, not refuted.
- **Binds:** fuel. Net of fuel and CO₂ it is **negative against subsonic at every Mach and corridor**. Capture is capped by airports. The spine's fibre is about the whole world's UHMWPE capacity.
- **B/C (2025 year, integrated cost):** vs Boomless **0.87 gross / 0.04 net** at M2.2, **0.79 / 0.28** at M2.0. Life-cycle net reaches 1 only for a **~2065 opening**.
- **Verdict:** keep at desk level; re-run ~2029–31 when the FAA rule and Overture's economics are known.

**Skyloop** (floating near-vacuum tube at 20 km; [review](skyloop-review/README.md))
- **Passes:** the time case. A ground portal and launch ramp give **2.7 h NYC–LAX at M4**. Electric drive at **21–85 Wh/pax-km** makes net B/C exceed gross.
- **Binds:** alignment. The tube moves **8–485 m** in gusts against **0.14–1.3 m** allowed. No feasible point at 20 km with a fixed guideway ([regime](skyloop-review/regime/regime.md)); a tensioned line plus a steered guideway in a 40 m bore passes screens on paper but stays uneconomic ([hybrid](skyloop-review/hybrid/hybrid.md)), and wings only add cost.
- **B/C (2050 opening):** net vs subsonic **0.47 base, 1.67 low cost** ([system](skyloop-review/system/system.md)).
- **Verdict:** do not pursue the stratospheric tube.

**Hyperloop, our version** (elevated ground tube, compliant-gap SC-EDS maglev; [README](hyperloop/README.md))
- **Passes:** an active 0.3 Hz / 0.3 m gap rides a realistic pylon guideway to **M13–20**. Plains legs need almost no tunnel. No physics veto to ~M6.
- **Binds:** breach deceleration (**pass to M3, marginal to M6**), then right-of-way above ~M4, then cost. Mountains force **690 km of tunnel at M4** and lift the route average to **$120 M/km**.
- **B/C (2050 opening, M3–4):** net vs subsonic **0.50–0.53 base, 1.27–1.36 low cost**. At 30 % capture, ≤ 0.63.
- **Verdict:** keep as the ground-vactrain watching brief at **M3–4**. It is a cost and certification problem, not a Mach problem.

## Comparison (NYC–LAX)

**Basis.** Tube columns: life-cycle, **2050 opening**, M4, spine+LAX. BoomShield: M2.2, 50 t aircraft, integrated cost, with the STATUS **2025-year** canonical value and the **2050 life-cycle** value shown separately. The two BoomShield bases differ by ~2× gross; compare only the 2050 values with the tubes. "Low" is each study's own low-cost case (Skyloop run $19 M/km; Hyperloop $51 M/km; BoomShield's integrated cost is already its low-cost case).

| | BoomShield M2.2 | Skyloop, floating 20 km (M4) | Hyperloop, elevated (M4) |
|---|---|---|---|
| Door to door (subsonic 7.29 h, Boomless 6.27 h) | 5.08 h | 2.73 h, if buildable | **2.73 h** nonstop, 3.62 h all-stop |
| Physics status · binding | acoustics pass a screen · net of fuel, capture | **alignment veto** (8–485 m motion) | **passes to M3, marginal to M6** · breach, then right-of-way |
| Infrastructure, $M/km | 25 integrated – 186 base ¹ | 66 run (19–187) + 157 ramp (40–623) | 80 standard, **120 route average** (51–300) |
| Corridor capex | $91 B spine ¹ ($672 B at base) | $569 B spine+LAX ($172 B low) | **$552 B** spine+LAX ($233–1,344 B) |
| Energy, Wh/pax-km | ~1,400 (jet fuel) | 21–85 (electric) | 15 nonstop · 59 all-stop |
| B/C vs subsonic, gross / net | 2025 yr: 2.69 / −4.41 · 2050 LC: 5.17 / −5.18 | base 0.36 / 0.47 · low 1.28 / 1.67 | base **0.37 / 0.50** · low 0.94 / 1.27 |
| B/C vs Boomless, gross / net | 2025 yr: 0.87 / 0.04 · 2050 LC: 1.67 / 0.46 | base 0.28 / 0.61 · low 1.09 / 2.17 | base 0.29 / 0.64 · low 0.72 / 1.62 ² |
| Pays below (net, 100 %) | does not pay net in any 2025-year cost case; net vs Boomless = 1 only for a ~2065 (M2.0) or ~2075 (M2.2) opening | **$40–55 M/km** all-in; ~$11–14 M at 30 % | **$52 M/km** vs subsonic ($34 M gross), $71 M vs Boomless; ~$10 M at 30 % |

¹ Integrated cost for the Concorde-class (55 m) shield, $14.3/m²; the thinner M2.2 shield costs less ($1.62 vs $2.22 M/km/yr, [net-reconciliation](boom-shield/economics/net/net-reconciliation.md)).
² Gross from `hyperloop/hl_econ.csv` (M4, 2050, 100 %); the README quotes net only.

## Cross-cutting lessons

- **Ground processing beats Mach.** At M2 only 1.3 h of BoomShield's 5.2 h NYC–LAX is cruise; a 70 min walk-on saves more than M2 → M∞. The tubes win by boarding at ground level, not by speed: Skyloop M2 → M6 saves 70 min, and Hyperloop B/C peaks at **M3–4**, because energy grows as V² and all-stop legs never reach cruise.
- **Net of fuel flips sign between aircraft and electric tubes.** BoomShield loses most of its gross value to fuel and CO₂; vs subsonic, net is negative everywhere. The tubes avoid jet fuel, so net is **1.3–2.2× gross**. The fuel and carbon price sets the gap in both directions.
- **Capture and moving-block certification dominate.** Every headline is 100 % capture; at 30 % divide by about three. For the tubes, brick-wall braking cuts B/C from ~0.5 to **0.04**. Only an uncertified **10 s moving block** carries the spine in one tube pair, and it puts 24 pods within stopping distance of a breach at M4.
- **Cost per km is the gate.** BoomShield needs $106 → **$14/m²** by owning the plants. The tubes break even near **$40–55 M/km**, against $120 M/km priced with mountains. No vactrain has been built at that cost; Hyperloop One's own Bay Area estimate was $52–75 M/km.
- **Physics moves the binding constraint, then cost takes over.** The compliant gap turns alignment from fatal into a non-issue on the ground; breach then caps Mach near M3–6. Nothing rescues the floating tube.

## What to watch, and gates

**BoomShield** (timing comes from events, not traffic)
- T1: FAA final rule (targeted mid-2027) and Boomless Cruise status, 2027–28 → **G-A ~2029**.
- **Full re-run 2029–31**, once T1 and T2 (Overture's real fares and capture, 2030–32) resolve; T4 film and fibre at ≤ $20/m², 2028–31.
- Later gates: lab acoustics G-B ~2032, subscale module G-C ~2038, spine decision G-D ~2045, earliest opening ~2055.

**Hyperloop (our version)** (the study sets no dates)
- A **priced elevated km ≤ $50 M including mountains**.
- A regulator's **path to ~10 s headways** (relative-braking operation).
- A **breach-mitigation demonstration**: sector valves ≤ 1 s with early braking (frontier M3 / 6 → M4 / 11).
- External data point: CASIC T-Flight's 60 km results (~2030, system.md).

**Skyloop:** no gates; the stratospheric tube is stopped.

## Links

- BoomShield: [STATUS](boom-shield/STATUS.md), [README](boom-shield/README.md), [net reconciliation](boom-shield/economics/net/net-reconciliation.md)
- Skyloop re-review: [README](skyloop-review/README.md), [system](skyloop-review/system/system.md), [regime](skyloop-review/regime/regime.md), [hybrid](skyloop-review/hybrid/hybrid.md)
- Hyperloop (our version): [README](hyperloop/README.md)
