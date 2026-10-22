# BoomShield net B/C: one convention for every study

**Skyloop · 22 October 2026 · screening level.** `net_bc.py` re-uses highmach, demand, integration and growth code read-only and writes the CSVs here. **These numbers supersede every net B/C quoted elsewhere.**

## Canonical convention

- **Basis:** a single year in steady state at 2025 traffic, in 2025 dollars. The shield cost is annualised (7 %, 25 yr) net of dual use, with two tether rows. Base induced demand; capture 100 % (30 % = 0.3×).
- **Gross:** time saved at $86/h × 0.95, plus day-return credits and induced-demand surplus.
- **Net:** gross minus the extra fuel ($0.80/kg) and CO₂ (3.16 kg/kg at $190/t) over **whole trips** (great circle × 1.03: feeders, ramps and shield km alike), incremental to the comparator. Boomless means a domestic-sized Overture-class aircraft; subsonic means 90 g CO₂/pax-km. **Induced trips carry their full fuel.** These are highmach's and global's prices, and growth's 2025 starting point.
- **One aircraft per Mach** for shield, ramps, fares and fuel: M1.7 62 t, M2.0 49 t, M2.2 50 t, M3.0 63 t. The Concorde-class M2.0 row (shield for 147 t, fuel of the 49 t jet) is kept only as the legacy reference.
- **Integrated cost is re-run** with integration's plant model on each design's own shield, not scaled by the M2.0 ratio of 0.552.
- **Life-cycle** results are a separate column: growth_model, base world, central energy path.

## Assumptions compared

| Item | highmach | global | growth | Canonical |
|---|---|---|---|---|
| Basis | 1 yr, 2025 | 1 yr, 2025 | Life-cycle NPV, opening 2035–2100, traffic + VOT growth | 1 yr, 2025 |
| Fuel / CO₂ | $0.80/kg, $190/t | same | $1.43 → 1.91/kg all-in (SAF, SC-CO₂ $200 → 308/t, −0.8 %/yr fuel) | $0.80, $190 |
| Km counted | Whole trip ×1.03 | same | inherits highmach | same |
| Induced trips' fuel | Counted | **Omitted** (benefit counted) | Counted | Counted |
| Baseline fuel | Overture-class; 90 g/pax-km | same | same | same |
| Seats, LF | 65, 0.75 | same | same | same |
| Aircraft, M2.0 | 49 t (also 147 t row) | 49 t | **147 t shield + 49 t fuel** | 49 t |
| Integrated cost | Aggressive × 0.552, all Mach | same | M2.0 147 t; M2.2 × 0.70 | Re-run per shield |
| Tether rows | 2 | 2 | 2 | 2 |
| Capture | 100 %, 10 % | 100 %, 30 % | 100 %, 30 % | 100 %, 30 % |

## What caused each discrepancy

Bridge, spine, integrated, 100 %, net vs Boomless (`reconciliation.csv`):

| Step | M2.2 | M2.0 49 t | Type |
|---|---|---|---|
| global, published | 0.123 | 0.395 | — |
| highmach routing instead of global's re-implementation | 0.110 | 0.387 | method, ±1–2 % |
| + fuel and CO₂ of induced trips | 0.040 | 0.283 | **error in global** |
| = highmach aggressive ÷ 0.552 | 0.040 | 0.283 | — |
| Integrated cost re-run per shield (M2.2 +3.8 %) | **0.039** | **0.275** | approximation |

| Step (M2.0) | Net vs BL | Type |
|---|---|---|
| growth, published (2035 opening) | 0.196 | — |
| same at flat canonical prices | 0.240 | method (energy path) |
| single year 2025 | 0.123 | method (life-cycle) |
| 49 t aircraft for shield and ramps too | **0.275** | **inconsistent mass in growth** |

- **Global's 0.12 / 0.40** omit the fuel and CO₂ of the 1.1–1.7 M induced trips a year, whose benefit they include through highmach's uplift. That is worth −0.07 at M2.2 and −0.10 at M2.0.
- **highmach's +0.16 / +0.02** are at aggressive cost. Dividing by 0.552 is right (0.28 / 0.04), but the ratio understates the narrower shields' integrated cost by 3–5 %.
- **Growth's 0.20** combines three things: a life-cycle basis (+0.12 against a single year), the central energy path (−0.04 in 2035), and a mixed aircraft (a 147 t shield at $2.22 M/km/yr with the 49 t jet's fuel). That aircraft's shield costs 45 % more and saves 9 % less time than the 49 t design's.
- **Gross values agree** within routing noise (M2.2 0.90 / 0.90; M2.0 0.81 / 0.80; 147 t row 0.49–0.50). The re-run cost lowers them 0–5 %.

## Canonical table

**US spine, integrated cost, base demand, 100 % capture**

| Design | $M/km/yr | Gross vs BL | Net vs BL | Gross vs sub | Net vs sub | Extra CO₂ vs BL, Mt/yr | Life-cycle net vs BL, opening 2035 / 2045 / 2050 |
|---|---|---|---|---|---|---|---|
| M1.7 Overture, 62 t | 1.50 | 0.50 | 0.41 | 2.35 | −3.52 | 1.1 | 0.52 / 0.73 / 0.85 |
| **M2.0 adaptive, 49 t** | 1.54 | **0.79** | **0.28** | **2.67** | **−4.01** | 6.4 | 0.40 / 0.63 / 0.75 |
| **M2.2 design A, 50 t** | 1.62 | **0.87** | **0.04** | **2.69** | **−4.41** | 11.0 | 0.14 / 0.34 / 0.46 |
| M2.5 Ti, 57 t | 1.82 | 0.91 | −0.91 | 2.60 | −5.90 | 27 | — |
| M3.0 design B, 63 t | 1.96 | 0.98 | −2.37 | 2.57 | −8.16 | 53 | — |
| M2.0 147 t shield (legacy) | 2.22 | 0.49 | 0.12 | 1.75 | −2.82 | 6.7 | 0.20 / 0.36 / 0.46 |

- **Other cost cases**, gross / net vs BL. M2.0 (49 t): base 0.070 / 0.024, plausible 0.19 / 0.065, aggressive 0.45 / 0.16. M2.2: 0.077 / 0.003, 0.21 / 0.009, 0.50 / 0.022.
- **High demand, M2.2:** gross 1.47, net −0.07.
- **Life-cycle crossings** (base world, 100 %), in `crossing_years.csv`. Gross reaches 1 for a 2040 opening (M2.0) or 2035 (M2.2). Net reaches 1 at **~2065 (M2.0), 2060 (M1.7), 2075 (M2.2)**. At 30 % capture, not before 2100. The 49–50 t jets are not expected before the late 2030s.
- **At M1.7, net ≈ gross minus induced-trip fuel.** No study has an M1.15 off-design fuel term for Boomless.

**Corridors** (M2.2 / M2.0; spine direct, the others rescaled from global's `corridors.csv`)

| Corridor | Gross vs BL | Net vs BL | Net vs sub |
|---|---|---|---|
| US north spine | 0.87 / 0.79 | 0.04 / 0.28 | −4.4 / −4.0 |
| US south | 0.74 / 0.65 | −0.05 / 0.17 | −4.5 / −4.0 |
| Europe–Gulf LON–IST–DXB | 0.46 / 0.41 | 0.00 / 0.14 | −3.3 / −2.9 |
| Europe–Levant–Gulf (greedy) | 0.55 / – | −0.02 / – | −3.5 / – |
| LON–IST–DXB–BOM–HKG | 0.36 / 0.32 | 0.00 / 0.10 | −2.6 / −2.3 |

Net vs subsonic is negative everywhere. In the other corridors, net vs Boomless runs from −0.17 (China, at M2.2) to +0.05.

## Superseded

- **global.md / ranking.csv:**
  - the net column (0.12 / 0.40 and every corridor's) and "M2.0: much better net": use `corridors.csv` here;
  - the integrated gross values, which are ~4 % high (0.90 → 0.87).
- **highmach.md:**
  - the Table 4 "integrated" columns (ratio-scaled): 0.90 → 0.87, 1.03 → 0.98, high demand 1.53 → 1.47;
  - "net figures scale by ~1.8". The aggressive-cost net values stand.
- **growth.md:**
  - the M2.0 rows of Tables 4–5 are the 147 t legacy shield;
  - the canonical M2.0 life-cycle net is 0.40 / 0.63 / 0.75, and net = 1 comes at ~2065, not 2075;
  - the M2.2 rows are ~4 % high.
- **integration.md, STATUS:** "net ~0.1–0.2" = 0.12 (147 t) or 0.28 (49 t). STATUS and the infographic are updated.

Non-CO₂ effects and H₂ leakage are not counted; they would lower every net figure.
