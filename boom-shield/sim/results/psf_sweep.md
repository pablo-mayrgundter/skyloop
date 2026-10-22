# BoomShield against the proposed 0.11 psf overland limit

From `run_psf_sweep.py` (`--phase validate|coarse|refine|offtrack|extra|3d|report`; cached in `psf_cache/`). Data: `psf_sweep.csv`; figure: `psf_sweep.png`. The limit is a ground peak of **≤ 0.11 psf = 5.27 Pa** (NPRM, docket FAA-2026-6935).

## Verdict

- **Undertrack: reachable for every aircraft, given enough H₂.** At 20 km the floating lid needs about **20 m of H₂ for Overture M1.7, 28 m for Concorde M1.4, 54 m for Concorde M2, 71 m for SR-71 M3.2 and 77 m for Concorde M2.5**. Helium needs about twice that. The **30 m baseline fails at Mach 2 (12.0 Pa)**.
- **Whole carpet:** in new 3-D runs (Concorde M2, H₂ 60 m), **W = 1 km fails off-track (11.7 Pa at 39 km)**, while **W = 1.5 or 2 km pass everywhere (≤ 4.2 Pa)**.
- **No hard floor.** The peak scales as W·d^−(1.2…1.5), so the limit is a weight × thickness budget: about **66 t per 30 m of H₂ at Mach 2**.
- **PLdB and psf have come apart.** In 1-D every passing case is inaudible, and so are 289 failing ones. psf binds.

## Setup

- **Chain.** The Gate-0 1-D undertrack chain, with the aircraft at H_s + 250 m and preset weights (Overture-like preset from `make_data.py`).
- **Stacks.** Floating lid: lift − 0.1 / gas / 0.1 kg/m². Fixed lid: 2.4 / gas / 0.1. Untied.
- **Grid.** 1,020 cases: He/H₂ × d 15–200 m × H_s 14–25 km × 5 aircraft points × 2 lids, plus baselines.
- **Trimmed.** He at Mach 3.2; 3-D checks for Concorde M2 only.
- **Numerics.** 5 kHz below the shield for smooth waves (equal to 50 kHz within 0.01 Pa). Gate 0 is reproduced exactly.

Unshielded peaks: Concorde 1.55–2.15 psf, Overture 0.96–1.13, SR-71 0.96–1.44. The flags below mark unrealistic points:

- **ᶜ C_L > 0.5:** e.g. Concorde M1.4 at ≥ 20 km.
- **ᵠ q > 35 kPa:** e.g. SR-71 at ≤ 20 km.

## Minimum thickness for ≤ 0.11 psf (1-D undertrack, m, floating lid)

| aircraft | gas | 14 km | 16 | 18 | 20 | 22 | 25 km |
|---|---|---|---|---|---|---|---|
| Concorde M1.4 | H₂ | 21 | 23 | 25 | 28ᶜ | 32ᶜ | 37ᶜ |
| | He | 39 | 42 | 46 | 51ᶜ | 57ᶜ | 67ᶜ |
| Overture M1.7 | H₂ | 16 | 16 | 18 | 20 | 23 | 27ᶜ |
| | He | 28 | 30 | 33 | 37 | 41 | 49ᶜ |
| Concorde M2.0 | H₂ | 41ᵠ | 44 | 49 | **54** | 60 | 71ᶜ |
| | He | 73ᵠ | 79 | 87 | 96 | 107 | 126ᶜ |
| Concorde M2.5 | H₂ | 58ᵠ | 63ᵠ | 70 | 77 | 86 | 103 |
| | He | 117ᵠ | 127ᵠ | 139 | 154 | 171 | >200 |
| SR-71 M3.2 | H₂ | 58ᵠ | 60ᵠ | 64ᵠ | 71ᵠ | 79 | 94 |

Log-interpolated, and slightly conservative for thick layers. **Fixed 2.4 kg/m² lid:** equal at d ≤ 30 m, but it does not float there. It is worse when thick (Concorde M2 needs 61 m; SR-71 > 100 m). **Put spare lift into the lid.**

## Floor physics

**The layer is a low-pass filter.** It passes |T|⁻² ≈ 1 + (ωτ)², with

τ² = a²d(d + m/ρ_g) + ¼[Z_a a²d/ρ_g − (m + ρ_g d)/Z_a]², where a² = 1/V² − 1/c_g², m = face mass and Z_a = ρ_a/√(1/c² − 1/V²).

This matches the exact |T|. For H₂ 30 m at Mach 2, τ = 0.20 s; τ ∝ d.

**The k → 0 part is the lift.** The y-integrated incident pressure, which is what the undertrack far field sees, is exactly (ρV²/2β)·A_e′: the lift distribution. A floating shield must pass it (T(0) = 1) and can only spread it over ℓ = Vτ. The ground hears a lifting line of weight W and length ℓ:

**p ≈ C · 1.9 R(H) · βW / (2π√(2βh) ℓ^{3/2})**, with R(H) the ray gain to the ground (0.62 at 14 km, 0.84 at 20 km, 1.13 at 25 km).

One constant, **C = 1.50**, fits 699 cases within **±1.8 dB rms**. The peak is linear in W, and the volume contributes 2 %.

**No irreducible floor.** The unremovable quasi-static ground footprint, W spread over about πβH², is only **0.1–1.2 mPa**. Maximum weight for 0.11 psf (20 km, H₂, floating lid):

| d (H₂) | M1.4 | M1.7 | M2.0 | M2.5 | M3.2 (SR-71 shape) |
|---|---|---|---|---|---|
| 30 m | 162 t | 104 t | **66 t** | 39 t | 17 t |
| 60 m | 402 t | 261 t | 174 t | 104 t | 43 t |
| 100 m | 752 t | 497 t | 346 t | 215 t | 91 t |
| 200 m | 1,660 t | 1,120 t | 818 t | 541 t | 241 t |

(Concorde shape except for the M1.7 and M3.2 columns.)

H₂ beats He because its larger impedance contrast lowers the corner frequency.

## Altitude dependence

At fixed d the peak **rises 2.0–2.25×** from 14 to 25 km. That is R(H): a βW-set pressure gains impedance on the way down from thin air. The minimum d grows 1.7× while the lift per metre falls 5.7×, so **lower is better**, down to the q and traffic limits.

## 1-D screen vs 3-D (Concorde unless noted, h 250 m, H_s 20 km)

| case | 1-D | 3-D undertrack | 3-D carpet maximum |
|---|---|---|---|
| M2 H₂ 30, W 1 km (existing) | 11.8 | 11.3 | 11.5 (track) |
| M2 He 30, W 1 km | 25.8 | 25.4 | 25.5 (track) |
| M2 He 60, W 1 km | 9.6 | 9.2 | **12.5 @ 39 km** |
| SR-71 M3.2 H₂ 30, W 1 km | 14.4 | 14.1 | 13.9 (track) |
| M1.4 He 30, W 1 km | 11.2 | 10.7 | 10.8 (track) |
| **M2 H₂ 60 floating, W 1 km (new)** | 4.54 | 4.20 | **11.7 @ 39 km** |
| **same, W 1.5 km (new)** | 4.54 | 4.20 | **4.20** (track); 4.0 @ 39 km |
| **same, W 2 km (new)** | 4.54 | 4.19 | **4.19** (track); ≤ 1.7 beyond 14 km |

**Undertrack, the claim holds.** The 1-D peak is 2–8 % conservative.

**Over the carpet, it does not hold for thick layers on a 1 km shield.**

- The outer carpet climbs to 10–13 Pa whatever the stack.
- An infinite shield on the same off-track rays gives only 0.4–0.8 Pa (unshielded off-track values match 3-D: 40.8 vs 40.7 Pa).
- So the excess is edge diffraction of the ~1 Hz swell: those rays cross only 150–250 m inside the edge, about one cross-plane wavelength.
- 1.5 km of width removes it.

**Checks**

- **Band limit.** The new 3-D runs are peak-only at f ≤ 30 Hz. A 30 Hz rerun of He 60 m matches the 100 Hz peaks within 4 %.
- **Near-field record.** Extending it from 3 L to 12 L lowers thick-layer peaks by up to 29 % at 200 m, so the tables are conservative.
- **Padding.** No effect.
- **Soft ties.** K = 10 Pa/m helps (3.7 vs 4.5 Pa), but only through tuning.

## PLdB versus psf

- **1-D.** The residual of every passing case is a 1–5 Pa swell with 0.4–2.7 s rise, at PLdB −∞.
- **3-D.** PLdB is set by 3–8 Pa edge pulses above 5 Hz, not by the peak.
- **What each metric responds to.** Width and edges buy PLdB. Thickness, lower altitude and lighter aircraft buy psf.
- **The rule.** A peak rule penalises an inaudible infrasonic pulse, where building rattle is the real issue. Worth a docket comment.

## Recommendations

1. **Re-baseline for psf.**
   - 150 t at Mach 2 and 20 km: **H₂ ≈ 60 m and W ≈ 1.5 km**.
   - Overture-class: H₂ 20–30 m.
2. **Treat weight × Mach as the market limit.** Heavy aircraft at Mach ≥ 2.5 need ≥ 80–100 m.
3. **Prefer 16–18 km** if q and traffic allow.
4. **Use H₂ with a floating lid.**
5. **Next:** 3-D width runs for the other aircraft (about 10 min each at 30 Hz), graded edges, and a building-response metric.
