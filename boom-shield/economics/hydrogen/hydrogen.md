# Cheap hydrogen: what it changes for BoomShield, Skyloop and Hyperloop

**Skyloop · 23 October 2026 · screening level.** `hydrogen.py` re-uses the canonical net code (`net/net_bc.py`), the integration and cost-down cost models, growth_model's life cycle, `sky_model.lifecycle` and `hyperloop.bc_rows`, all read-only, and writes the CSVs here. It reproduces the published canonical, Skyloop and Hyperloop net B/C to < 0.002 before changing anything. **Canonical numbers are unchanged; these are scenarios beside them.**

**The claim tested.** Terraform Industries (Casey Handmer) says its sunlight → H₂ demo is below **$2/kg** now and plans **< $1/kg** "in the coming years". Its March 2024 demo reported < $2.50/kg H₂ and < $250/t DAC CO₂. A July 2026 column puts its Rosamond PV at ~$17/MWh against a whitepaper assumption of $10/MWh. A Liebreich-style critique argues e-methane stays well above Henry Hub. All are **company claims or search summaries, unverified** (humanprogress.org summary, undated).

## Bottom line

- **Direct H₂ cost is immaterial.** H₂ is **~4 % of shield capex** (3.8–4.4 % integrated, 4.1 % at cost_model base). The canonical integrated case already buys grey H₂ at **$2.1/kg delivered**. Green H₂ at $2 adds **+0.9 %** to net cost; at $1, **−1.4 %**. Skyloop's envelope gas is ≤ 2 % of its run cost.
- **Cheap e-SAF helps BoomShield only below the canonical $1.40/kg.** That is fossil fuel at $0.80 plus CO₂ at $190/t. Central e-SAF is **$2.76/kg at H₂ $2** and **$1.61/kg at H₂ $1**, so carbon-neutral fuel makes net B/C *worse* than canonical at $2. At $1 it is slightly worse on the 2025 basis and better for a 2050 opening, where growth's central path already reaches $1.91/kg.
- **BoomShield turns positive against subsonic only at e-SAF ≤ ~$0.55/kg (2025 year) or ≤ ~$1.0/kg (2050 opening).** That needs the low corner of the $1 scenario: H₂ ≤ ~$1, DAC ≤ $50/t and cheap conversion. Net vs Boomless reaches 1 for a 2050 opening at **≤ $1.30/kg (M2.0)** or **≤ $1.10/kg (M2.2)**; in a 2025 year never, because gross is below 1.
- **The material surprise: cheap solar erases the shield's PV dual-use credit.** At $17/MWh ground PV it falls from **$0.62 to $0.16 M/km/yr**; at $10, to hosting only. Net shield cost rises **+30–36 %** and gross B/C vs Boomless drops **0.87 → 0.65–0.67** (M2.2). This is the same cheap-solar world that makes H₂ cheap.
- **The tubes' net advantage survives.** It disappears only at e-SAF **≤ $0.21–0.53/kg**, below the e-SAF floor (**$0.45/kg** even at free H₂). At $1.61 their net B/C falls ~6 % (Hyperloop M4 **0.50 → 0.47** base, **1.27 → 1.19** low); at $2.76 it rises ~12 %.
- **Solar traction power is a small lever.** Firmed for 24 h service, cheap PV costs **~$56–66/MWh** ($33–110), not $10–17. Against the repo's $70 plus grid CO₂, it adds only **+0.01–0.02 (base) and +0.02–0.06 (low)** to tube net B/C.

## Scenarios (`assumptions.csv`)

| | Baseline (repo) | H₂ $2/kg (Terraform demo) | H₂ $1/kg (Terraform / DOE target) |
|---|---|---|---|
| H₂ at plant, e-SAF | $6 ($5–7; growth T6, G21) | $2 ($1.5–2.5) | $1 ($0.8–1.5) |
| Shield H₂, delivered | $2.1 grey SMR (integration) | $2 + $0.5 = $2.5 | $1 + $0.5 = $1.5 |
| DAC CO₂ | $500/t ($250–1,000) | $250/t ($150–400) | $100/t ($50–200) |
| Ground PV | $28/MWh (cost-down aggressive) | $17/MWh (Rosamond) | $10/MWh (whitepaper) |
| **e-SAF, central (range)** | **$6.0/kg** (3.7–10.5) | **$2.76/kg** (1.55–4.90) | **$1.61/kg** (0.85–3.37) |
| Implied abatement cost | $1,650/t | $620/t | $256/t |

## 1. Direct H₂ cost (`direct_h2.csv`)

**BoomShield.** The integration plant model is re-run on each shield with the scenario's H₂ price. Inventory is 225–293 t/km; make-up is 3 %/yr at aggressive levers.

| Spine, integrated cost | M2.0 (49 t) | M2.2 (50 t) |
|---|---|---|
| H₂ share of capex / of net annual cost | 3.8 % / 4.8 % | 4.4 % / 5.4 % |
| Net cost, baseline | $1.54 M/km/yr | $1.62 M/km/yr |
| H₂ $2 + $0.5 delivery | +$0.014 (+0.9 %) | +$0.017 (+1.0 %) |
| H₂ $1 + $0.5 delivery | −$0.021 (−1.4 %) | −$0.025 (−1.5 %) |
| H₂ $1, delivery $0.3 / $1.5 | −$0.028 / +$0.014 | −$0.034 / +$0.017 |

- **Delivered price is what counts, and at 20 km delivery is a large part of it.** At $1/kg at plant, delivery of $0.3–1.5 (compression, pipeline, lift up a tether hose) is 23–60 % of the delivered price. At $1.5 delivery, $1 green H₂ costs more than today's $2.1 grey.
- **growth.md's "< 4 % of capex" is about right: ~4 %.** Safety, not price, is still the H₂ issue (STATUS T6).
- **Skyloop envelopes** (R 43–88 m, 40–166 t/km): the gas is **$1.3 M/km** at $5.5 delivered, 1.9 % of the $66 M/km run; $0.35 M/km at $1.5. Make-up at 10 %/yr is ≤ $0.08 M/km/yr against $1.7 M of O&M. **Immaterial,** and the floating tube is stopped anyway.

## 1b. Cheap solar and the shield's PV credit (`direct_h2.csv`, PV variants)

Cost-down's aggressive lever puts PV on 15 % of the lid. Its credit is capped at the ground PV cost of the same era. The shield PV costs ~$20/MWh, against $28 on the ground today.

| Ground PV (with the scenario's H₂) | Dual-use credit, $M/km/yr | Net cost M2.2 | Gross vs BL, M2.2 / M2.0 |
|---|---|---|---|
| $28/MWh (baseline) | 0.62 | 1.62 | 0.87 / 0.79 |
| $17/MWh (shield PV cost falls in step) | 0.16 | 2.11 | 0.67 / 0.60 |
| $17/MWh (shield PV cost fixed) | 0.055 | 2.19 | — |
| $10/MWh | 0.055 (hosting only) | 2.17 | 0.65 / 0.58 |

- The credit is gone below **~$15/MWh** (shield PV cost falls in step) or **~$20/MWh** (fixed). Cheap ground PV also floods daytime prices, so the cap is credible.
- **This is the largest single effect of the Terraform world on BoomShield**, larger than the H₂ itself. It is a cost, not a fuel, effect, so it hits gross and net alike.

## 2. e-SAF and a carbon-neutral fuel convention (`esaf_cost.csv`, `boomshield_net.csv`, `tubes_net.csv`, `breakeven.csv`)

**e-SAF bottom-up.** RWGS + Fischer–Tropsch or methanol-to-jet share a floor of **0.43 kg H₂ and 3.14 kg CO₂ per kg** of liquid fuel. With 75 % H₂ and 82 % carbon efficiency: 0.58 kg H₂ and 3.8 kg CO₂ per kg. Conversion is $1,500/kW of fuel at 85 % capacity factor, **$0.30/kg**; H₂ buffer storage $0.40/kg H₂, so synthesis runs steadily on daytime H₂.

| e-SAF, central, $/kg | H₂ | DAC CO₂ | Conversion | Other | **Total** |
|---|---|---|---|---|---|
| H₂ $6, DAC $500 | 3.68 | 1.91 | 0.30 | 0.12 | **6.02** |
| H₂ $2, DAC $250 | 1.38 | 0.96 | 0.30 | 0.12 | **2.76** |
| H₂ $1, DAC $100 | 0.81 | 0.38 | 0.30 | 0.12 | **1.61** |
| Floor: free H₂, DAC $50, low-cost plant | 0.08 | 0.17 | 0.14 | 0.06 | **0.45** |

- e-SAF matches the canonical **$1.40/kg** at H₂ **$0.64** (central plant, DAC $100) or **$1.54** (low-cost plant). At DAC $250 the central plant cannot reach it at any H₂ price. **DAC carbon matters as much as the H₂.** The critique of e-methane applies: even at $1/kg H₂, e-SAF is ~2× fossil jet.

**Convention.** Extra (or avoided) fuel is priced at the e-SAF price, flat in real terms, with **no CO₂ charge** on it. Non-CO₂ effects and H₂ leakage stay uncounted, as now. Grid electricity keeps its CO₂ charge unless a row says "clean".

**BoomShield, spine, integrated cost, 100 % capture, net B/C** (H₂-only costs; PV credit kept)

| Fuel | M2.0 vs BL: 2025 yr / LC 2050 | M2.2 vs BL | M2.0 vs subsonic | M2.2 vs subsonic |
|---|---|---|---|---|
| Canonical ($0.80 + $190/t; LC central path) | **0.28 / 0.75** | **0.04 / 0.46** | −4.0 / −4.5 | −4.4 / −5.2 |
| e-SAF $6.0 (baseline H₂) | −1.42 / −0.73 | −2.69 / −2.00 | −26 / −24 | −28 / −26 |
| e-SAF $2.76 (H₂ $2) | −0.22 / 0.46 | −0.76 / −0.01 | −10.4 / −8.2 | −11.2 / −9.2 |
| **e-SAF $1.61 (H₂ $1)** | **0.20 / 0.89** | **−0.09 / 0.69** | −5.1 / −2.7 | −5.6 / −3.2 |
| e-SAF $0.85 ($1 low corner) | 0.48 / 1.16 | 0.37 / 1.15 | −1.4 / **+0.92** | −1.7 / **+0.72** |

With the cheap-solar PV erosion as well (H₂ $1, e-SAF $1.61): M2.0 vs BL 0.15 / 0.79, M2.2 −0.06 / 0.59.

**Break-even e-SAF prices** (100 % capture; H₂ at plant needed with a central plant and DAC $100, or a low-cost plant and DAC $50)

| Target | Basis | M2.0 (49 t) | M2.2 (50 t) | H₂ needed (central / low-cost) |
|---|---|---|---|---|
| Net vs subsonic > 0 | 2025 year | **$0.56/kg** | **$0.53/kg** | none / ≤ $0.16–0.22 |
| | LC, 2045 opening | $0.96 | $0.91 | none / ≤ $0.91–1.01 |
| | LC, 2050 opening | **$1.05** | **$0.99** | ≤ $0.02 / ≤ $1.07–1.18 |
| Net vs Boomless = 1 | 2025 year | never (gross 0.79) | never (gross 0.87) | — |
| | LC, 2045 opening | $0.93 | $0.87 | none / ≤ $0.83–0.96 |
| | LC, 2050 opening | **$1.30** | **$1.10** | ≤ $0.11–0.47 / ≤ $1.28–1.69 |

- With the PV erosion, net vs Boomless = 1 for a 2050 opening needs **$0.97–1.03** (M2.0) or **$0.81–0.85** (M2.2). The subsonic break-even does not depend on cost.
- At 30 % capture net vs Boomless never reaches 1. The subsonic break-even is unchanged, since gross and fuel scale together.
- **STATUS T7 ("SAF ≤ $1.5/kg sets the sign vs subsonic") is too loose.** At $1.5/kg net vs subsonic is still about −2 to −3 for a 2050 opening. The sign flips at **≤ ~$1.0/kg** carbon-neutral fuel, and a carbon price alone never flips it.

**Tubes, spine+LAX, 2050 opening, 100 % capture, net vs subsonic, base / low cost**

| Fuel · electricity | Skyloop M4 | Hyperloop M3 | Hyperloop M4 |
|---|---|---|---|
| Gross (unchanged) | 0.36 / 1.28 | 0.39 / 0.98 | 0.37 / 0.94 |
| Canonical (growth central path, grid $70 + CO₂) | **0.47 / 1.67** | **0.53 / 1.36** | **0.50 / 1.27** |
| e-SAF $2.76 · grid | 0.53 / 1.88 | 0.60 / 1.52 | 0.56 / 1.42 |
| e-SAF $2.76 · firmed clean solar $66 | 0.55 / 1.93 | 0.61 / 1.54 | 0.57 / 1.45 |
| **e-SAF $1.61 · grid** | **0.44 / 1.57** | **0.50 / 1.28** | **0.47 / 1.19** |
| e-SAF $1.61 · firmed clean solar $56 | 0.46 / 1.64 | 0.51 / 1.30 | 0.48 / 1.23 |
| e-SAF $0.85 · grid | 0.39 / 1.37 | 0.44 / 1.12 | 0.41 / 1.04 |

- **Net equals gross** (avoided fuel = electricity cost) at e-SAF **$0.53** (Skyloop), **$0.21** (Hyperloop M3) and **$0.36/kg** (M4) against subsonic, and $0.08–0.19 against Boomless. All but Skyloop's sit below the $0.45 floor, so **the tubes' net > gross advantage does not disappear at any H₂ price.**
- **Hyperloop low-cost net vs subsonic stays ≥ 1 for e-SAF ≥ $0.65/kg** (M4) or $0.30 (M3), so in every scenario. Skyloop's low case is ≥ 1 at any fuel price.
- **Expensive carbon-neutral fuel helps the tubes.** A mandate at $2.76 lifts them ~12 %; cheap e-SAF at $1.61 trims ~6 %.

## 3. Solar electricity for the tubes (`electricity.csv`, `tubes_net.csv`)

- **Daytime prices do not apply to 24 h service.** Pods run ~18 h a day, so 50–70 % of traction energy must come from storage. Firmed PV = PV × 1.2–1.6 overbuild + the shifted share × battery cost ($120–250/kWh, one cycle a day, $40–115 per MWh shifted).
- **Firmed: $56/MWh at $10 PV and $66/MWh at $17 PV** (range $33–110). The repo uses $70/MWh plus grid CO₂ at the SCC: **$109/MWh all-in in 2050** ($131 in 2035).
- Firmed clean solar instead of grid adds **+0.01–0.02 (base) and +0.02–0.06 (low)** to tube net B/C. The daytime-price bound adds only ~0.01–0.03 more. **Immaterial:** electricity is ~$2 B/yr against ~$14 B of benefit on the spine+LAX.
- Launch peaks (1–5 GW per station) need storage or grid capacity whatever the energy price; the station cost already carries the power plant.

## 4. H₂-fuelled supersonic aircraft (no model)

**Not relevant to this screening.** Liquid H₂ has ~2.8× jet fuel's energy per kg but needs ~4× the volume, so a supersonic LH₂ airliner is longer and fatter. That adds wave drag and boom-relevant volume, offsetting its lower weight; the shield is sized by weight and length, so the net effect on the layer is unclear. LH₂ at $1/kg H₂ plus ~$1–1.5/kg liquefaction and airport delivery costs about $0.9–1.25 per kg-of-jet-fuel-equivalent energy, no cheaper than cheap e-SAF. It would also need new aircraft, engines and airport fuel systems, none before the 2050s. And burning H₂ emits **~2.6× more water per unit of energy** than kerosene, at 18–20 km where stratospheric water vapour persists. That is the non-CO₂ term this study already leaves out, and it gets worse, not better. If cheap H₂ matters to supersonic flight, it is through e-SAF.

## Caveats

- Terraform's figures are company claims via a search summary; the $1/kg plan rests on $10/MWh PV that Rosamond (~$17) has not shown. Solar-only electrolysis runs at a 25–30 % capacity factor.
- e-SAF plant costs, efficiencies, buffer storage and firming parameters are engineering assumptions [U]; no plant has run at scale.
- e-SAF is held flat in real terms over each life cycle; the canonical life cycle uses growth's central path (all-in $1.43 → 1.91/kg, 2025–50).
- The PV erosion uses cost-down's market cap rule unchanged; a different dual-use rule would change its size, not its sign.

## Files

- `hydrogen.py` (run: `python3 -I boom-shield/economics/hydrogen/hydrogen.py`, ~10 s)
- `esaf_cost.csv` (bottom-up e-SAF by scenario, plus an H₂ sweep), `electricity.csv` (firmed solar vs repo grid)
- `direct_h2.csv` (shield H₂ and PV variants, delivery sensitivity, Skyloop envelope gas)
- `boomshield_net.csv` (M1.7 / M2.0 / M2.2 × BL / sub × 30 / 100 % × 2025 year and LC 2035 / 2045 / 2050 × fuel convention × cost variant)
- `tubes_net.csv` (Skyloop M4, Hyperloop M3 / M4; fuel × electricity variants), `breakeven.csv`, `assumptions.csv`
