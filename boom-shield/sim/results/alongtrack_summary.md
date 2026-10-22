# Along-track variation: serrated edges and module joints (3-D)

`run_alongtrack.py` / `boomsim/alongtrack.py`. Data `alongtrack_results.csv`, `alongtrack_budget.csv`; figures `alongtrack_{carpets,orders,joints}.png`. Concorde M2, h 250 m, W 1 km, H₂ 30 m slab. Δ = less favourable of raw/proxy ΔPLdB vs unshielded.

## Bottom line

1. **Serrations act as a grating.** The x-averaged model drops τ(1−τ) of the aperture field (28 % of the n = 0 energy). Once orders k_x,n = ω/V + 2πn/Λ reach the ground, up to 45 % of it lands at y ≈ 20–37 km, 20–35 dB above the x-averaged residual.
2. **Λ = 10–30 m is worse than a sharp edge** at 23–35 km (up to +6 dB). At Λ = 100–300 m the near track falls from −∞ to −38…−43 (sharp −34). **The full benefit survives only for periodic Λ ≲ 1 m**; Λ = 2 m loses up to 12 dB at 27 km. Random phase fills k_x ≈ 0 at every frequency, so it does not help at short Λ.
3. **Open joints are unacceptable:** 1 m every 250 m gives −24.5 instead of −33.7 undertrack; ~1 dB needs g ≲ 0.1 m. **Membrane-only joints** are fine to g ≈ 1 m (L 250) or 2 m (L 500).

## Method

- **3-D PO + Floquet orders.** Aperture field p_inc + χ(x,y)·D(y), with D crossplane's PO difference field and χ the plan-view mask. Exact Fourier coefficients c_n(y) give x-plane-wave orders, each an exact 2-D problem (κ_n² = k₀² − k_x,n²); c₀ = τ(y) is the x-averaged model.
- **Shapes** share the raised-cosine τ(y), so n = 0 is the cached BEM fringe run: one-sided sawtooth, symmetric teeth, and random-phase (symmetric teeth randomly shifted per period, 8-period supercell; two realisations agree within 1 dB).
- **Ground.** Each order's angular spectrum is carried as ray-tube energy flux along its own slowness (k_x,n, k_y)/ω. It includes the stratified Jacobian, absorption, reflection and anchor. The x-averaged ESD Σ_n|P_n|² is added in 1/3-octave bands to the specular ground signature (cached BEM circle spectra through the unchanged nonlinear `propagate_from_circle`). For joints, the coherent n = 0 leak (g/L)(D_gap − D_stack) is added on the circle first. Scattered fields go to 400 Hz; 1 kHz changes Δ by ≤ 0.3 dB.
- **Re-steepening.** The loudest scattered band (63–79 Hz) was taken at ground amplitude, compressed into one 70 ms burst (conservative). Run back along its ray through `propagate()`, it gains **+1.9 to +3.4 dB** nonlinear vs linear (Σ ≤ 1.2). Scattered levels may be ~3 dB optimistic; no verdict changes.

## Validation

| check | result |
|---|---|
| PO D(y) vs `crossplane.po_solve` | 5e-13 |
| n = 0 of a sawtooth (exact c₀) vs crossplane fringe PO, 12/30/60 Hz | ≤ 3e-13 |
| Sum of orders vs brute-force 3-D Rayleigh integral (toy, test) | < 2 % |
| Flux mapping vs linear `propagate_from_circle` (unshielded, 10–60 Hz, φ 0–39°) | ≤ 0.03 dB |
| Ray tables vs `RayPath` (offset, Jacobian) | 1e-6 |
| Re-propagated cached spectra vs `3d_results.csv` (sharp, fringe) | identical (0.00 dB) |
| Sharp edge, PO vs BEM on the circle (5–100 Hz) | within 0.5–5 dB (PO louder near track) |
| x-averaged fringe, PO vs BEM ≥ 20 Hz | PO 10–50 dB lower: the BEM n = 0 values are floor-limited, so the fringe residual used here is an upper bound |

## Q1. Serrated edges (w_f 150 m): Δ PLdB vs ground offset

R_m 3 km (first column 8 km). First order reaches the ground above f₁ = 1/(Λ(1/V + 1/c_g)) = 216, 72, 22, 7, 2 Hz for Λ = 1…100 m (actual onset ~1.4× higher: the edge field has lateral slowness).

| case | 0 km (8 km R_m) | 9.7 | 18.3 | 23.5 | 26.9 | 31.5 | 39.3 |
|---|---|---|---|---|---|---|---|
| sharp edge | −33.7 | −29.4 | −23.8 | −21.7 | −20.2 | −19.0 | −17.0 |
| x-averaged fringe (n = 0) | −∞ | −63.3 | −51.9 | −49.5 | −47.4 | −40.1 | −24.7 |
| sawtooth Λ 1 m | −∞ | −63.3 | −51.9 | −49.5 | −47.4 | −39.9 | −24.7 |
| sawtooth Λ 2 m (to 1 kHz) | −∞ | −63.3 | −51.9 | −49.1 | −35.1 | −32.8 | −24.5 |
| sawtooth Λ 3 m | −∞ | −63.3 | −51.9 | −47.0 | −23.9 | −25.6 | −23.6 |
| sawtooth Λ 10 m | −∞ | −63.3 | −51.5 | **−17.2** | **−13.8** | −19.4 | −22.5 |
| sawtooth Λ 30 m | −∞ | −62.5 | −24.8 | −20.9 | **−16.7** | **−15.7** | −17.8 |
| sawtooth Λ 100 m | −40.2 | −34.5 | −28.1 | −24.9 | −20.9 | −19.5 | −16.5 |
| sawtooth Λ 300 m | −38.5 | −35.1 | −29.9 | −26.7 | −24.1 | −21.5 | −18.6 |
| random Λ 1 m | −∞ | −63.3 | −51.8 | −47.6 | −34.6 | −34.2 | −24.6 |
| random Λ 3 m | −∞ | −63.2 | −51.6 | −41.3 | −24.5 | −24.9 | −23.5 |
| random Λ 10 m | −∞ | −62.8 | −49.9 | −28.1 | **−15.6** | **−16.6** | −21.6 |
| random Λ 30 m | −108.8 | −61.9 | −40.2 | −21.9 | −20.0 | **−16.8** | **−15.8** |
| random Λ 100 m | −76.4 | −39.3 | −30.1 | −26.9 | −24.9 | −22.6 | −18.9 |
| random Λ 300 m | −42.7 | −38.4 | −31.9 | −28.7 | −26.7 | −24.4 | −20.6 |

Bold: louder than sharp. Symmetric teeth match random phase within ~1 dB for Λ ≥ 10 m (CSV); the one-sided sawtooth is up to 11 dB worse (its order phase steers energy inward).

- **Where the orders go.** Orders with k_x,n ≈ 0 leave nearly vertically along-track. They keep the edge field's outward lateral slowness, which the specular wave would carry beyond cut-off, and land at 20–37 km. Some reach 50–57 km, beyond the 47.7 km unshielded carpet edge, at up to 64 PLdB.
- **Long Λ (100–300 m)** gives many weak near-specular orders. They fill the near track and leave the outer carpet 0–4 dB better than the sharp edge.
- **Fringe depth.** A deeper fringe scatters more: the scattered-only maximum for w_f = 75 / 150 / 300 m is 75 / 80 / 83 PLdB at Λ = 10 m. Depth does not open the window.

**Design window: periodic Λ ≲ 1 m** (no ground-reaching order below ~300 Hz), any w_f; a 2 m pitch is a marginal fallback. Such teeth are sub-metre slivers of a 30 m stack. PO overstates their opacity, since the gas field's lateral decay length is 2–20 m. A real apodised edge is better built as a sheet with properties graded in y than as plan-view teeth. The 10–300 m pitches a modular build would produce do not deliver the fringe benefit.

## Q2. Module joints (sharp edge, gap across the full width)

| joint | undertrack, 8 km R_m | 0 km (3 km) | 9.7 km | 14.4 km | 31.5 km |
|---|---|---|---|---|---|
| none (sharp) | −33.7 | −38.5 | −29.4 | −26.1 | −19.0 |
| open, g 0.1 / 0.25 / 0.5 m, L 250 m | – | −37.2 / −34.3 / −30.1 | – | – | −19.0 |
| open, g 1 m, L 250 / 500 m | −24.5 / −27.1 | −24.9 / −28.0 | −24.2 / −26.1 | −23.1 / −24.4 | −18.8 / −18.9 |
| open, g 5 m, L 250 / 500 m | −13.9 / −17.2 | −14.0 / −17.3 | −14.5 / −17.7 | −14.9 / −17.9 | −16.8 / −17.7 |
| open, g 20 m, L 250 / 500 m | −6.5 / −9.8 | −6.5 / −9.8 | −7.2 / −10.4 | −7.8 / −10.9 | −10.5 / −13.2 |
| membranes, g 1 m, L 250 / 500 m | −33.4 / −33.5 | −37.8 / −38.2 | −29.3 / −29.3 | −26.0 / −26.1 | −19.0 / −19.0 |
| membranes, g 2 m, L 250 / 500 m | – | −36.6 / −37.5 | – | – | −19.0 / −19.1 |
| membranes, g 5 m, L 250 / 500 m | −31.0 / −32.1 | −33.6 / −35.5 | −28.6 / −29.0 | −25.5 / −25.8 | −19.0 / −19.0 |
| membranes, g 20 m, L 250 / 500 m | −26.0 / −28.6 | −27.1 / −30.2 | −25.9 / −27.4 | −23.2 / −24.4 | −17.9 / −18.8 |

- **Where it lands:** the incident passes from beneath the shield, so the leak covers **the whole carpet, peaking under the track**.
- **Scattered orders dominate:** n = 0 alone gives −20.6 at worst (g 20 m, L 250 m) vs −6.5 with orders. Order energy grows ∝ g/L for g ≳ λ, ∝ (g/L)² for g ≪ λ.
- **Allowable gap.** Open: ≲ 0.1 m for ~1 dB (0.25 m already costs 3–4 dB). Membranes only: ≲ 1 m (L 250 m) or 2 m (L 500 m) for ≤ 1 dB; 5 m costs 3–5 dB at R_m 3 km (2.7 dB at 8 km).

## Limitations

PO teeth and slots (narrow teeth less opaque than assumed; a 30 m-deep slot treated as thin). Launch points ±75 m (fringe), ±500 m (joint). Energies are x-averaged; few-order cases (Λ ≤ 10 m) can vary ±3 dB along x. The BEM n = 0 fringe is floor-limited (upper bound), and near-track R_m 3 km values are optimistic. No wind or turbulence; upward and ducted orders are not followed.
