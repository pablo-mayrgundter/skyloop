# BoomShield handoff — review addendum

**Skyloop project · 20 October 2026 · Review of [2026-10-20-boomshield-handoff.md](2026-10-20-boomshield-handoff.md)**

This addendum reviews the BoomShield concept-and-design handoff against the repository record and adds independent back-of-envelope checks. The handoff's own corrections hold up: its buoyancy (≈0.89 kg/m² for a 15 m helium envelope at 20 km; ≈67 m thickness to float 4 kg/m²; ≈56 t of capacity vs 250 t of mass for S2), catenary (≈5.1 kN/m), and crosswind (≈600 N/m) numbers reproduce. Its separation of inherited concept, prior model result, and handoff assessment is the right discipline for chat-derived work.

Where this review pushes further: the handoff treats acoustics, buoyancy, and clearance as three independent gates. Two of them are already coupled by first-order physics, and a one-page calculation suggests they may not close at the envelope thicknesses being considered.

All numbers below come from [quick_checks.py](quick_checks.py) (numpy only). They are rough screening estimates, not designs.

---

## 1. The mass law couples acoustics to buoyancy

A free-floating deck has nothing rigid to push against, so at boom frequencies (~5–30 Hz) its transmission is governed by areal mass. At 20 km the characteristic air impedance is low (ρc ≈ 26 rayl, about 16× lower than at sea level), which helps a light sheet. But buoyant capacity also scales with ρ, so the two cancel:

- mass needed for amplitude transmission A at frequency f (limp, normal incidence): m ≈ r·ρλ/π, with r = √(1/A² − 1) and λ = c/f;
- lift available from an elliptical gas envelope of thickness T: m_avail ≈ (ρ − ρ_gas)·πT/4 ≈ 0.86ρ·πT/4.

Setting these equal gives **T ≳ 0.47·r·λ**, independent of altitude to first order: a buoyant limp shield must be of order one to several boom wavelengths thick. For a 10 Hz component at 20 km (λ ≈ 30 m): A = 0.4 needs T ≈ 32 m; A = 0.2 needs T ≈ 68 m; A = 0.05 needs T ≈ 280 m — and that assumes zero structural mass.

Time-domain check, 0.1 s N-wave at 20 km, **all** lift spent on acoustically active mass (zero structure):

| Envelope T | Lift capacity | Peak passed, single limp sheet | Peak passed, double wall (He gap = T) |
|---:|---:|---:|---:|
| 15 m | 0.89 kg/m² | 0.66 | 0.54 |
| 30 m | 1.8 kg/m² | 0.40 | 0.27 |
| 60 m | 3.6 kg/m² | 0.25 | 0.14 |
| 100 m | 6.0 kg/m² | 0.17 | 0.09 |
| 200 m | 11.9 kg/m² | 0.10 | 0.03 |

The acoustic thread assumed A_mat ≈ 0.02–0.05 while the structural thread aims at sub-kg/m² decks. On this check those are inconsistent: A_mat ≈ 0.05 needs a >100–200 m envelope before any structure is added.

Caveats:

- **Peak is not loudness.** The mass law preferentially removes high-frequency (shock rise) content, which PLdB weights heavily. But the transmitted wave then travels ~20 km to the ground and will partly re-steepen. Whether the low-pass-filtered wave re-shocks is the single most valuable calculation to run next.
- **Shaping does not escape the mass law.** Canting, curving, or "mangling" the surface redistributes the reflected field; it does not reduce how much is transmitted. At oblique incidence the effective mass term is m·cos θ, which is worse, not better.

## 2. Cross-track width is independent of Mach in uniform air

The handoff correctly retires W_req ≈ 2h/√(M² − 1), but understates the correction. In uniform, windless air a ray launched at azimuth φ from straight down meets a horizontal plane a distance h below at

**y = h·tan φ,**

independent of Mach number. Mach sets only the along-track offset. So "higher Mach → narrower shield" is wrong, and in a real atmosphere the ground carpet widens with Mach, which makes the required width larger, not smaller. This also undercuts the segment notes' "sweet regime above X-59" conclusion.

For reference, the 1 km / 250 m acoustic case (W/h = 4) intercepts rays out to φ ≈ ±63°; a 250 m shield at h = 100 m covers about ±51°.

## 3. The source cases do not bracket an unshaped N-wave

The opacity table's shape penalties Cs = 1–3 all assume some degree of X-59-class shaping. A real SR-71 produces an ordinary N-wave of roughly 1 psf at the ground. As a rough estimate, that is mid-to-high 90s PLdB, not the ~83 PLdB used for the "SR-71-like optimistic" case. That implies about 20 dB or more of required reduction, not ~8 dB. Cs ≈ 3 is closer to optimistic than pessimistic for an unshaped aircraft.

## 4. The boom load, not dead load, sizes the deck

A rough Carlson-method estimate (no ground reflection) for 250 m below a Mach 3.2, ~33 m aircraft gives ~100–200 Pa of overpressure (K_p = 0.06–0.13). Deck dead load is 9 Pa at 0.9 kg/m² and 39 Pa at 4 kg/m². The deck therefore sees a pressure pulse 3–20× its own weight on every passage. The catenary dead-load calculation is not the governing case; per-passage fatigue under impulsive loading is.

## 5. Helium inventory

For a 4,000 km corridor at 20 km altitude:

| Section | Gas volume at altitude | At STP | Years of world helium production (~170 million m³/yr) |
|---|---:|---:|---:|
| W = 250 m, T = 15 m | 1.2 × 10¹⁰ m³ | 0.8 × 10⁹ m³ | ≈ 5 |
| W = 1 km, T = 15 m | 4.7 × 10¹⁰ m³ | 3.2 × 10⁹ m³ | ≈ 19 |
| W = 250 m, T = 66 m | 5.2 × 10¹⁰ m³ | 3.5 × 10⁹ m³ | ≈ 21 |

This likely forces hydrogen — at 75–250 m beneath supersonic aircraft — for any gas-supported branch at corridor scale.

## 6. Tether drag through the jet stream

A 5 cm tether crossing ~8 km of jet-stream winds (ρ ≈ 0.41 kg/m³, U = 60–80 m/s, C_D ≈ 1.2) carries roughly 350–630 kN of drag. That exceeds the handoff's illustrative 150 kN side load on a whole 250 m module. Tether drag, not envelope drag, likely sets the station-keeping design. The closest precedent is the 20 km tethered-balloon analysis for the SPICE geoengineering project (Davidson, Burgoyne, Hunt & Causier, "Lifting options for stratospheric aerosol geoengineering: advantages of tethered balloon systems," *Phil. Trans. R. Soc. A* 370, 2012).

## 7. The competitive and regulatory landscape has moved

As of my knowledge (to be verified before external use):

- Boom Supersonic's XB-1 has demonstrated Mach-cutoff ("boomless") supersonic flight, where the boom refracts away before reaching the ground — at low supersonic Mach.
- NASA's X-59 is flying.
- A 2025 executive order directed the FAA toward a noise-based standard for overland supersonic flight in place of the outright ban.

That leaves BoomShield's distinctive niche at roughly Mach 2 and above, where neither Mach cutoff nor current shaping reaches — and where points 1–4 are hardest. Sections 1–2 of the handoff should state this explicitly.

---

## Proposed Gate 0 (ahead of the handoff's Gates 1–5)

A closure test, done before anything else:

**σ_acoustic(A_req) ≤ σ_lift(T, H)**

evaluated with the mass law (and layered-stack transmission) at the shield, followed by nonlinear propagation of the transmitted wave to the ground and a perceived-level metric there — not peak amplitude. This is days of work, not a programme. If it fails at every thickness, the floating-limp-shield family is closed, and candidate escapes include an evacuated double wall (vacuum is a near-perfect isolator, but the shell must carry ~5.5 kPa of external pressure) or a much lower-Mach aircraft (which changes the wavelength to be blocked).

## Next: simulation

The Gate 0 model is specified in [`../sim/DESIGN.md`](../sim/DESIGN.md): Whitham F-function near field → augmented-Burgers propagation along the undertrack ray → layered transfer-matrix shield → propagation to ground → Stevens Mark VII PLdB.

One new hypothesis to test there: **a helium or hydrogen layer is "subsonic" relative to the aircraft's trace speed.** On any horizontal surface the boom's pattern moves at the aircraft speed V. If V is below the gas sound speed (helium ≈ 866 m/s, hydrogen ≈ 1,120 m/s at 216.65 K — i.e. below about Mach 2.9 or Mach 3.8 respectively at 20 km), the boom cannot propagate as a wave inside the layer and can only tunnel through it evanescently, decaying roughly as exp(−κd) with κ = ω·√(1/V² − 1/c_gas²). That is a possible mechanism beyond the plain mass law, and it gets stronger as the lifting layer gets thicker — the same direction the lift budget pushes.
