# Secondary sonic booms with BoomShield

Code: `boomsim/secondary.py` and `run_secondary.py`. Tests: `tests/test_secondary.py`. Data: `secondary_booms.csv`, with record = bin or ray.

## Bottom line

1. **The shield does not make secondary booms worse.** Unshielded, the downward boom reflects off the ground and climbs anyway, returning as an NPRM *indirect* secondary boom. The shield moves that reflection to 20 km, onto rays with the same slowness. In seven cases, the shielded worst bin is **0.54–1.02×** the unshielded one, and exceedances shrink to the inner carpet edge (150–230 km), where direct booms already land.
2. **The limit is broken without any shield.** Unshaped Concorde- or SR-71-class aircraft at 20 km exceed 0.11 psf whenever winds form a stratospheric duct: **0.14–0.37 psf** in four cases. In calm air or a headwind, only thermospheric returns land, at ≤ 0.02 psf. This is an aircraft and routing problem.

## NPRM definitions

Source: FAA, "Enabling Supersonic Overland Flight", 91 FR 40470 (2 July 2026), docket FAA-2026-6935. Proposed § 91.817(a)(2)(i): surface overpressure ≤ 0.11 psf for primary booms and for **direct and indirect secondary** booms. The § 1.1 definitions below are paraphrased from Holland & Knight (9 Jul 2026) and DLA Piper, because the FR text was unreachable here.

- **Primary boom:** goes straight from the aircraft to the surface.
- **Secondary boom:** refracted or reflected so that it travels farther.
- **Direct secondary:** "over the top", refracted back to the surface by the upper atmosphere.
- **Indirect secondary:** reflected by the surface first, then refracted back.

Applied here:
- A shield-reflected boom fits neither definition literally. It is scored as a secondary boom.
- Cut-off rays are counted as direct.
- Terrain reflections are not modelled.

## Method

- **Atmosphere.** USSA76 to 86 km. From 86 to 200 km: USSA76 temperature formulas, tabulated molar mass, hydrostatic pressure, and γ and c with atomic O included (c = 405 m/s at 120 km, 567 m/s at 150 km).
- **Winds.** 40°N CIRA-86-like profiles (Fleming et al. 1990): winter westerly 64 m/s at 60 km, summer easterly −60 m/s at 65 km. No tides.
- **Rays.** 3-D Hamiltonian rays in time, in a moving stratified medium, launched on the Mach cone (φ = 0–180°, 0.1° steps). The ground and the shield top act as mirrors.
- **Amplitude.** Blokhintzev invariant with tube area Q = (e_x × ∂r/∂φ)·v_ray, which stays finite at turning points. A caustic is applied as a Hilbert transform.
- **Waveform.** Augmented Burgers (`propagate.burgers_step`) from the Whitham near field at each φ; the upward lobe uses A_v − A_L.
  - Absorption: Sutherland & Bass (2004) classical plus rotational, ∝ 1/p (0.43 dB/km at 1 Hz at 100 km), plus ISO relaxation.
  - Shield: `shield.transfer` R(f) or T(f) at the ray's trace speed, for H₂ 30 m with 1.24 kg/m² faces. All rays with |φ| < 63.4° hit the strip.
  - Ground reflection ×0.9; linear below apex − 15 km. Δp/p peaks at 0.05–0.17.
- **Ground levels.** Energy per 3 km lateral bin and branch, p² ∝ (ρc²Ω/|v_z|) Σ E_h Δφ/ΔY, which stays finite at caustics. Ground factor ×1.9; peak = max|p|.

## Validation

| Check | Model | Reference |
|---|---|---|
| p, ρ at 86–200 km | ≤ 1 % | USSA76 tables |
| Ground classical absorption | ≤ 1 % | ISO 9613-1 |
| Primary boom vs `ground_signature` (Concorde, SR-71; φ = 0°, 40°) | peak ≤ 0.5 %, PLdB ≤ 0.1 dB | — |
| Turning heights | 41–51 km; 111–140 km | 40–50; 100–130 km (Balachandran, Donn & Rind, *Science* 197, 1977) |
| Delay after overflight | 8–13 min (stratospheric); 18–31 min (thermospheric) | "10–15 min" (Rogers & Maglieri, *Acoustics Today*, 2015) |
| Range from track | 150–390 km; 110–550 km | 165–1000 km (Balachandran 1977) |
| Amplitude | 5–18 Pa (stratospheric); 0.1–1 Pa (thermospheric) | "1–10 Pa" (Rogers & Maglieri); "~0.2 psf" (Falcon-9 JASA-EL 2026, citing Concorde); 0.1 Pa far range (Liszka 1978)* |

\*Abstract only, unverified.

Stratospheric levels sit at the top of the observed range. The model has no turbulence or multipath smearing, so treat it as good to about ×2.

## Results (worst 3 km bin, psf)

| Aircraft, atmosphere | U direct | U indirect | S direct | S shield-reflected | S transmitted | **Max U → S** | PLdB; CSEL U → S |
|---|---|---|---|---|---|---|---|
| Concorde M2, calm | 0.008 | 0.014 | 0.008 | 0.006 | 0.002 | **0.014 → 0.008** | –; 16 → 8 |
| Concorde, winter tailwind | 0.19 | 0.23 | 0.19 | 0.19 | 0.018 | **0.23 → 0.19** | 65 → 61; 84 → 82 |
| Concorde, winter headwind | 0.009 | 0.020 | 0.009 | 0.017 | 0.003 | **0.020 → 0.017** | –; 20 → 17 |
| Concorde, winter crosswind | 0.18 | 0.20 | 0.18 | 0.205 | 0.016 | **0.20 → 0.205** | 65 → 62; 82 → 82 |
| Concorde, summer tailwind | 0.37 | 0.34 | 0.37 | 0.31 | 0.025 | **0.37 → 0.37** | 73 → 70; 91 → 88 |
| SR-71 M3.2, calm | 0.005 | 0.003 | 0.005 | 0.002 | 0.001 | **0.005 → 0.005** | –; 9 → 11 |
| SR-71, winter tailwind | 0.14 | 0.11 | 0.14 | 0.11 | 0.004 | **0.14 → 0.14** | 58 → 58; 81 → 80 |

U = unshielded, S = shielded. "–" means no Mark VII band above threshold (sub-5 Hz rumble). Exceedance zone: winter tailwind 170–386 → 170–209 km; summer 155–254 → 155–164 km.

**Why the shield adds nothing** (one ray, φ = 55°, winter tailwind; normalised ground energy):

| Case | Normalised ground energy |
|---|---|
| Perfect mirror at 20 km | 1037 |
| Unshielded ground reflection | 1016 |
| Real shield | 199 |

- The fresher wave reflected at 20 km loses more to shocks: −6.4 dB, against −3.3 dB for the aged wave.
- The real shield reflects 99.9 % of the incident energy (|R| ≥ 0.86 above 0.3 Hz). Its phase, however, swings from 150° at 1 Hz to 25° at 30 Hz (face mass and gas stiffness), so the reflected pulse is dispersed and less of it survives the duct.
- That 7 dB margin is model-specific. The robust result is the bound **shield-reflected ≤ perfect mirror ≈ unshielded indirect**.

## Uncertainty

- **Winds dominate:** a duct changes levels 15–40×. Tides and gravity waves can open ducts in any season. The profiles are representative only.
- **Caustics and binning.** Bins of 1.5 or 6 km move peaks by +10–45 % or −5–30 %. Adding branches coherently gives up to +75 %.
- **Not modelled:**
  - scattering into the gap between the primary and secondary carpets;
  - finite-strip diffraction of the reflected field below 2 Hz; the 1-D R likely overstates reflection there;
  - shield motion relative to the air;
  - shaped aircraft.
- **Trimmed after container restarts:** SR-71 crosswind; half the azimuths of the transmitted family when more than 24 were needed.

## Mitigations

- **Aircraft and route.** Route on wind climatology, avoiding ducting azimuths (eastbound in winter, westbound in summer). Shape the upward boom too. Cruise higher or slower.
- **Shield edges.** Only reflections at |y| ≈ 310–400 m (φ ≈ 51–58°) feed the stratospheric duct. A 5–10° upturn of the outer 150–200 m would steepen them by 10–20°, to trace speeds of at least 370 m/s. The modelled ducts accept only up to about 355 m/s, so those rays would become thermospheric returns at ≤ 0.02 psf. This needs a BEM check of the edge leak.
- **Absorptive top layer.** It acts mainly above a few hertz, while the ducted peak comes from 0.3–3 Hz content: low value.
- **Avoid a cleaner low-frequency mirror** (heavier or stiffer stacks). That would raise the reflected secondary toward the unshielded indirect level.

![secondary booms](secondary_booms.png)
