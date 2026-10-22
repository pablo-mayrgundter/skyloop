# BoomShield — the inverted (vacuum) balloon option

**Skyloop · 20 October 2026 · screening review.** 20 km US-1976 via `boomsim.atmosphere` (Gate-0 convention): 216.65 K, 5529 Pa, ρ = 0.0889 kg/m³, c = 295 m/s. 1-D, smeared screening numbers.

## Verdict

- **Acoustics.** A layer *with nothing holding it open* would be the best isolator modelled so far: it beats H₂ 30 m at every Mach, with no cut-off. But whatever holds the layer open is itself the stud. It couples the faces at roughly Δp/R, i.e. 10²–10⁵ Pa/m, against a tolerance of 1–300 Pa/m.
  - Raft of D ≈ 30 m tubes: ΔPLdB −0.3 to −1.9 (Mach 1.4–2.5) and −0.5 to −8 (Mach 3.2). That equals an air-filled double wall or the plain mass law.
  - Lens with struts: −2.6 to −6.7.
  - Only single tubes with D ≥ 100 m stay useful (−11 to −∞). But an H₂ tube of that size also floats with ~6 kg/m² to spare and is also opaque.
- **Structure.** Nothing floats (10–250 m, any architecture): structure is 1.23–1.44× lift; even ideal, buckling-free members take 0.64–0.97. The ratio is scale-invariant.
- **Overall vs He/H₂.** Worse on both counts. A same-shape He envelope nets 0.6–0.85 of the vacuum's *gross* lift; the vacuum nets −0.23 to −0.44. **Close the pure-vacuum branch** (§9 option B), keeping its lesson: anything that holds the faces apart against a static load is a stud.

## The sound-speed worry

The worry is correct but beside the point.

- c depends only on T, so air at 5.5 Pa still has c = 295 m/s < V. There is no evanescence.
- The near-vacuum works by two other effects:
  - its impedance is 0.026 rayl (against 26 rayl);
  - its gas spring is γp_int/(d·cos²θ) ≈ 0.3–0.5 Pa/m, which puts the mass–air–mass resonance at 0.05–0.3 Hz.
- Ideally this beats the He/H₂ mechanism. What fails is the structure (§C).

## A. Acoustic mechanism

**Gas checks at 5.5 Pa:** mean free path 0.9 mm (Kn ≈ 3×10⁻⁵, continuum holds); absorption 0.02 dB/100 m at 10 Hz, 2 dB/100 m at 100 Hz (negligible); thermal boundary layer 8 cm. Cavity resonances f_n = n·c/(2d·cosθ) (5–14 Hz) are very narrow (Q ≈ ωm/Z_gas).

**Shield plane, lift-closed faces (Gate-0 convention), ideal 0.1 s N-wave: peak ratio / positive-impulse ratio**

| stack | M 1.4 | M 2.0 | M 2.5 | M 3.2 | M 2.0, 0.3 s wave |
|---|---|---|---|---|---|
| He 15 m (Gate 0) | 0.054 / 0.34 | 0.187 / 0.57 | 0.328 / 0.72 | 0.657 / 0.89 | 0.535 / 0.90 |
| He 30 m | 0.015 / 0.18 | 0.064 / 0.34 | 0.147 / 0.47 | 0.460 / 0.88 | 0.305 / 0.71 |
| H₂ 30 m | 0.005 / 0.11 | 0.022 / 0.22 | 0.051 / 0.31 | 0.129 / 0.45 | 0.144 / 0.52 |
| H₂ 60 m | 0.001 / 0.05 | 0.006 / 0.11 | 0.014 / 0.16 | 0.043 / 0.25 | 0.047 / 0.31 |
| vacuum 10⁻¹, 30 m, *no structure* | 0.095 / 0.21 | 0.079 / 0.19 | 0.075 / 0.18 | 0.073 / 0.18 | 0.168 / 0.36 |
| vacuum 10⁻², 30 m, *no structure* | 0.010 / 0.02 | 0.008 / 0.02 | 0.008 / 0.02 | 0.007 / 0.02 | 0.018 / 0.04 |
| vacuum 10⁻³, 15–60 m, *no structure* | ≤ 0.002 / ≤ 0.004 | ≤ 0.001 / ≤ 0.003 | ≤ 0.001 / ≤ 0.003 | ≤ 0.002 / ≤ 0.003 | ≤ 0.004 / ≤ 0.009 |
| **raft D 30 m, ring K = 575–1569 Pa/m** | 0.31–0.49 / 0.74–0.96 | 0.31–0.47 / 0.76–0.94 | 0.31–0.47 / 0.76–0.93 | 0.31–0.47 / 0.76–0.93 | – |
| **lens T 30 m, struts K = 7.7×10⁴ Pa/m** | 0.43 / 0.75 | 0.38 / 0.69 | 0.37 / 0.67 | 0.36 / 0.66 | – |
| single tube D 100 m, ring K = 166–454 | 0.05–0.11 / 0.26–0.42 | 0.05–0.10 / 0.26–0.41 | 0.05–0.10 / 0.26–0.41 | 0.05–0.10 / 0.26–0.40 | – |
| single tube D 250 m, ring K = 63–173 | 0.008–0.018 / 0.10–0.17 | same | same | same | – |
| air-filled double wall, same faces as raft | 0.47 / 0.87 | 0.42 / 0.82 | 0.40 / 0.80 | 0.39 / 0.79 | – |
| mass law, same total mass | 0.43 / 0.81 | 0.38 / 0.75 | 0.37 / 0.74 | 0.36 / 0.73 | – |

**Ground ΔPLdB** (h = 250 m; headline = less favourable of raw and rise-proxy, as in Gate 0; the He/H₂ rows reproduce Gate 0 exactly)

| stack | Concorde M 1.4 | M 2.0 | M 2.5 | SR-71 M 3.2 | ΔCSEL, M 2.0 |
|---|---|---|---|---|---|
| He 15 m | −46.7 | −8.1 | −3.6 | −1.1 | −6.7 |
| H₂ 30 m | −∞ | −∞ | −∞ | −77.4 | −37.3 |
| vacuum 10⁻³, 15–60 m, no structure | −74 to −∞ | −68 to −∞ | −67 to −80 | −94 to −∞ | −41 to −62 |
| vacuum 10⁻², 30 m, no structure | −49.9 | −62.6 | −68.6 | −69.2 | −31.3 |
| vacuum 10⁻¹, 30 m, no structure | −8.5 | −15.8 | −20.4 | −37.9 | −10.6 |
| **raft D 30, rings K 575 / 1569** | −1.1 / −0.3 | −1.4 / −0.4 | −1.9 / −0.5 | −8.1 / −0.5 | −1.2 / −0.4 |
| **lens T 30, struts** | −2.6 | −4.2 | −4.8 | −6.7 | −3.8 |
| single tube D 100, K 166 / 454 | −88 / −11.0 | −∞ / −31.4 | −∞ / −45.0 | −∞ / −∞ | −33.4 / −17.3 |
| single tube D 250, K 63 / 173 | −∞ / −∞ | −∞ / −∞ | −∞ / −∞ | −91 / −91 | −59.3 / −48.9 |

Ideal-vacuum residuals below ~35 PLdB are numerical floor (cropped cavity ringing). Unlike evanescent He, the ideal vacuum also removes impulse (sub-Hz f₀).

## B. Membrane tension and coincidence

- **Tension.** Concave skin: N ≈ Δp·R_curv. Along-track UHMWPE tendons between rings (pitch s, sag 1/10) give N_xx = 1.25·Δp·s ≈ 0.7–8×10⁴ N/m.
- **Smeared face: no issue.** Designs give c_m = √(N_xx/m) = 110–175 m/s (0.12–0.42 V), so |m_eff|/m ≥ 0.82 (≤ 1 dB).
- **Without ties it never matters:** at Z_m → 0 the 1000:1 gas mismatch still holds ≤ −45 dB (`acoustics_coincidence.png`).
- **With ring coupling, c_m ≈ 0.6–1.0 V costs 3–4 dB** (peak 0.27 → 0.40–0.46).
- **c_m ≳ 2V (stiffness-controlled faces) helps**: peak 0.026 at M 1.4. But a fully stressed fibre caps c_m at √(σ_w/(SF·ρ)) ≈ 740 m/s (UHMWPE), so this works only below ~M 1.5.
- **Real risk: bare skin bays between rings.** Film and tendons alone give c_m = 310–620 m/s (M 1.05–2.1), so local coincidence can appear above V/(2s) ≈ 70–200 Hz, in the PLdB band. The smeared model cannot resolve this.

**Verdict on B:** a design lever (keep c_m away from V; tendons cross-track), not a fatal band.

## C. Load path: stability implies stiffness

A tension net cannot resist uniform external pressure. The compression members it needs are the acoustic short-circuit.

| architecture | compression members | top–bottom coupling K | value |
|---|---|---|---|
| (i)/(iii) tube or raft | rings (hoop Δp·R·s; q_cr = 3EI/R³), stringers (Δp·R/2 + tendon pull) | ring ovalisation ≈ (1.1–3)(F−1)Δp/R, F = q_cr/q ≈ 2.2–2.4 | D 10: 1.6–4.2 kPa/m; D 30: 0.58–1.6; D 100: 0.17–0.45; D 250: 0.06–0.17 |
| (ii) flat lens | vertical struts (Δp·a²) + in-plane compression grids (2N + Δp·T per unit width) | K = Δp/(ε·T) at working strain ε | T 10–100 m: 2.3×10⁵–2.6×10⁴ Pa/m |
| tolerance | – | He: 35–1000 Pa/m loses 6 dB. Vacuum: K = 30 already takes peak 0.001 → 0.023 | ≲ 10–100 Pa/m |

This cannot be engineered away, for two reasons:

1. **Every collapse mode needs a minimum stiffness**, ≥ (n²−1)(F−1)Δp/R. So K_ring ≳ (1–3)·Δp/R, 1–4× the *ambient* air column it replaced (γp/d ≈ 0.9·Δp/R): acoustically a double wall filled with something stiffer than air. The baselines confirm it (raft ≈ air-filled double wall ≈ mass law).
2. **Soft mounts cannot carry the static load.** A spring K carrying Δp deflects Δp/K: 55 m at 100 Pa/m. The only escape is quasi-zero stiffness (F → 1, at buckling), the same fragile tuning Gate 0 flagged for ties.

**Spoked/tensegrity rims** raise buckling capacity (an elastic foundation), but they cannot carry uniform pressure, they pre-compress the rim, and every spoke is a strut-class stud. They help buckling and hurt acoustics.

**The key trade.**
- **Confirmed:** lenses need struts, and the struts kill the acoustics.
- **Partly refuted:** tubes avoid *internal* members, but their rings are studs of 0.6–1.6 kPa/m at D = 30 m.
- **Refuted:** tubes do not cost planform (a raft gives ρ·πD/4 per m², like a lens of mean depth πD/4). They cost *height*: K ∝ 1/D and face mass ∝ D, so f₀ ∝ 1/D and tolerable coupling needs D ≥ 100 m, which still does not float.

## D. Mass closure (kg per m² of planform; 20 km, p_ratio 10⁻³, SF 1.5)

Assumptions: CFRP tube or 3-chord lattice (lighter of the two), 0.9 GPa usable, SP-8007 knock-down, L/500 and R/300 imperfections, q_cr ≥ 2q; UHMWPE at 0.8 GPa sustained; 25 µm film; joints +15 %/+10 %; He envelope 60 g/m².

| case | s / a [m] | lift (vacuum gross) | skin | compression | total | **structure / lift** | ideal floor | He net lift |
|---|---|---|---|---|---|---|---|---|
| tube D 10 | 1 | 0.70 | 0.22 | 0.78 | 1.00 | **1.44** | 0.97 | 0.41 |
| tube D 30 | 3 | 2.09 | 0.31 | 2.34 | 2.65 | **1.27** | 0.80 | 1.61 |
| tube D 100 | 12 | 6.98 | 0.72 | 7.88 | 8.60 | **1.23** | 0.77 | 5.82 |
| tube D 250 | 12 | 17.4 | 0.72 | 21.8 | 22.6 | **1.29** | 0.64 | 14.8 |
| lens T 15 (flat, rims free) | 3 | 1.33 | 0.20 | 1.65 | 1.85 | **1.39** | 0.89 | 1.03 |
| lens T 30 | 6 | 2.67 | 0.28 | 3.30 | 3.59 | **1.35** | 0.85 | 2.17 |
| lens T 60 | 12 | 5.33 | 0.46 | 6.60 | 7.06 | **1.33** | 0.83 | 4.47 |
| tube D 30, sea level (check) | 1 | 28.8 | 4.1 | 28.9 | 33.0 | 1.14 | 0.92 | – |

Full table (D 15–60, T 10/100) in `structure_closure.csv`.

- **Strength floor** (cylinder): M/M_lift ≥ 3·R_air·T·ρ_s/σ_s = 0.32 for CFRP without SF, independent of altitude. Matching He's net lift needs ρ_s/σ_s 4–6× better than CFRP.
- **Claim (3) is backwards for buckling-limited members.** Load falls 18× but lift falls 14×, and lightly loaded members are less efficient: sea level → 20 km takes a monocoque CFRP sphere from 10× to 33× its lift, and a ring-stiffened tube from 1.14 to 1.27.
- **Claim (5)** holds for gross lift (0.999 vs 0.86/0.93), not net. Pressure ratio barely affects structure.
- **`docs/07-sections.md`** says the invballoon is "entirely in tension", with ~1 t rings every 50 m. At R = 40 m that ring carries 11 MN of hoop compression and sizes to 10.5 t (7.2 t on strength alone).

## E. Pumping and leakage

Per km of corridor (W = 250 m raft of D 30 tubes: 7.9×10⁵ m² skin, 5.9×10⁶ m³). Rule of thumb: ~10 kW and 250 kg per m³/s of pump speed.

| source | pump speed at 10⁻³ (5.5 Pa) | at 10⁻¹ |
|---|---|---|
| permeation, 25 µm PET (≈0.025/0.005 Barrer O₂/N₂, 23 °C) | 0.17 m³/s | 0.002 m³/s |
| permeation, LDPE / with AlOx or EVOH barrier | 26 / 0.002–0.004 m³/s | 0.26 / ~0 |
| outgassing (10⁻⁶ Pa·m³/s/m²) | 0.14 m³/s | ~0 |
| **holes, 1 cm² total per km (choked)** | **10 m³/s ≈ 100 kW, 2.6 t of pumps** | 0.10 m³/s |
| holes, 10 cm² per km | 100 m³/s ≈ 1 MW, 26 t | 1 m³/s |

- **Holes set the pump size**, not permeation: one 1 mm² defect per 100 m² is 79 cm²/km.
- **1/1000 is the wrong target; use 10⁻¹–10⁻²**: 10–100× less pumping for ≤ 10 % less lift. Lift equals He's at p_ratio 0.14; the ideal 10⁻¹ layer already matches He 30 m at M 2; with real structure the ties dominate anyway.
- **Pump-down** at 10 m³/s/km: 16 days to 10⁻¹, 47 days to 10⁻³.

## F. Failure modes and operations

- **Buckling** is a snap-through and imperfection-sensitive. As with subsea pipelines (propagation ≪ initiation pressure; Kyriakides & Corona 2007, from memory), a collapse may run along a tube, so modules need arrestors. A 30 × 250 m module stores pV ≈ 1 GJ (0.23 t TNT).
- **Puncture** is structurally benign (Δp falls) but drains lift: 10 % of a module's lift in ~3 min through 1 m², ~20 days through 1 cm². Descent raises Δp, so a Δp-limited inflow valve is a natural fail-safe.
- **Thermal.** Lift is immune to superheat (an advantage). There is no internal convection; a 60 K film swing moves tendon N by ~10 % and c_m by ±5 %.
- **Start-up.** It cannot be evacuated on the ground (18× load). A constant-Δp ascent gives lift Δp/(R·T), 0.75 of design at the ground, but needs structure < 0.75 × lift. Otherwise: lift by gas, pump at altitude.

## G. Literature

Verified to search-abstract/catalogue level only (page fetches were blocked).

- **Akhmeteli & Gavrilin**, US 2006/0038062 A1, "Layered shell vacuum balloons": homogeneous shells buckle even in diamond; proposes a honeycomb sandwich.
- **Akhmeteli & Gavrilin**, *Eng* 2(4) 480–491 (2021), doi:10.3390/eng2040030 (arXiv:1903.05171): sea-level B₄C/Al-honeycomb sphere, R 2.5 m, W/B ≈ 0.9; not built.
- **Metlen**, AFIT-ENY-13-J-02 (2013); AIAA 2013-1881: Be isogrid sphere W/B 0.79; CFRP geodesic frame with skin 0.94; twin cylinders infeasible.
- **Cranston**, AFIT PhD (2016); US 10,843,783 (2020): icosahedral frame with skin.
- **Rapport**, US 10,625,842: tensegrity "vacuustat"; reports that an earlier tensegrity claim failed analysis.
- **None found** at stratospheric or corridor scale. Every design is a metre-scale sea-level sphere closing marginally; a sphere's strength floor is half a cylinder's.

## Recommendation and next 3-D tests

1. **Do not carry the vacuum branch forward**: it fails on mass at every size, and at buildable sizes its acoustics fall to the mass law.
2. **Carry its rule into the He/H₂ design**: static loads ride on gas pressure, not members between the faces; tie budget K ≲ 10–100 Pa/m. A vacuum "assist" pocket would bring ties with it.
3. **If the 3-D model revisits vacuum, test only what could change the acoustic verdict:**
   - (a) a periodic ring-and-skin model, to replace the smeared-K bracket;
   - (b) skin-bay coincidence at 70–200 Hz;
   - (c) a single D ≥ 100 m tube as a scatterer, not a slab;
   - (d) buckle propagation and arrestor spacing.

None of these fixes the mass closure.

---

**Files.**
- `vacuum_structure.py` → `structure_closure.{csv,png,txt}`
- `vacuum_acoustics.py` → `acoustics_{shield_plane,ties_tension,baselines,ground,ground_large}.csv`, `acoustics_transfer.png`, `acoustics_coincidence.png`
- `vacuum_ops.py` → `ops_pumping.csv`

**Caveats.**
- 1-D plane waves with smeared ties; the ring coupling is bracketed, not solved.
- Faces are lift-closed (optimistic).
- Permeabilities and pump ratings are handbook-class values; verify before use.
