# Skyloop re-review (October 2026)

A second look at the original Skyloop concept (`docs/`): a buoyant near-vacuum tube at ~20 km carrying pods at Mach 2–6, revisited with the physics and economics machinery built for [BoomShield](../boom-shield/STATUS.md). Screening-level desk study; web sources were checked only at search-summary level.

**Verdict: do not pursue the stratospheric tube. Keep watching ground vactrains.**

## Physics and engineering ([physics/physics.md](physics/physics.md))

Ranked problems:

1. **Alignment (likely fatal).** A floating tube moves ~8–485 m in 1–5 m/s gusts; at Mach 6 a pod tolerates ~0.14 m over a 6 km wavelength (1.3 m at Mach 2) for a 0.05 g ride. Tethers would need ~170× the stiffness, or active thrust 0.4–2 GW per 1,000 km.
2. **Headway (likely fatal for the business case).** Stopping from Mach 2–6 at 1 g takes 18–160 km; safe-braking spacing gives 120–360 passengers/h per tube, not the 12,000 claimed.
3. **Thermal bow (hard).** Cross-tube temperature differences must stay under ~0.24 K at Mach 6.
4. **Buoyancy (solvable, large).** Vacuum lift is 0.3–2.5 kg/m against 90–410 kg/m of tube and guideway; the vacuum "invballoon" nets −250 kg/m. Separate H₂ envelopes of 43–82 m radius are needed.
5. **Energy and geometry (hard).** 7–42 Wh per passenger-km (the docs claim 0.5–2), 0.6–5 GW per station, curve radii 355–3,200 km, launch runs 59–533 km.
6. **Safety (hard).** A Mach 6 pod carries ~1.1 t TNT equivalent; breach, H₂ next to vacuum, evacuation from 20 km and lightning.
7. **Vacuum (solvable).** 1–3 Pa is needed; leak pumping is modest.

The original docs' claims that fail: the invballoon is not "entirely in tension" and its lift is wrong; the tube's first mode is ~4 Hz, not 10; the intended pressure is 10⁻³ of *ambient* (~5.5 Pa at 20 km), not 10⁻³ atm, and that is fine for drag. But "gliders" cannot glide at 5.5 Pa: lift needs ~25–1,400 Pa in a ≤ 8 m bore, and heating closes gliding above ~M3. So the route needs a levitation track ([regime](regime/regime.md)).

## Feasible operating regime ([regime/regime.md](regime/regime.md))

The map spans tube pressure 0.1 Pa–1 kPa, Mach 1–8, pod size 2–30 t, bore, setting and headway rule, for both maglev and glide.

- **Floating at 20 km: no feasible point.** Alignment fails at every Mach (26 m motion vs 0.14–5 m allowed), with glide or maglev alike.
- **Glide works almost nowhere.** At 10⁻³ of ambient the lift is 10–100× short. The only window is Mach ~1.8–3 at 200–1,000 Pa, at 130–260 Wh/pax-km, near jet-fuel parity.
- **Technically open, ground maglev:** Mach 1–2 elevated or 1–2.8 buried, at ~6 Pa–1 kPa, with ≥ 33-seat pods on a 10 s moving block. It is marginal to M4 (alignment, curve radius) or M5.6 (breach deceleration). Leaks set the pressure floor and heating the ceiling, not drag.
- **B/C ≥ 1 only for an elevated tube at ~$40 M/km with a certified moving block** (M1.4–4: 1.0–1.4 at low cost, ≤ 0.75 at base). Buried never pays (≤ 0.51).
- **To open more:** active suspension or 5 mm/km alignment (M2 → 5–8), straighter right-of-way, moving-block certification and cheaper tubes. Nothing rescues the floating tube.

## Floating-tube stability ([stability/stability.md](stability/stability.md))

A spectral study of a tethered tube at 20 km, covering gusts, tethers, keel, taut line, active trim and big bores, at Mach 1–4.

- **Steady loads on discrete tethers bind first, not gusts.** Wind changes add ~3 kN/m and the daily lift swing ~2 kN/m, 10–50× the km-scale gust loads. They need load-spreading bridles.
- **Tethers are stiffness-capped.** Any fibre layout gives at most ~0.074 N/m per kg of fibre, and bracing does not help.
- **Tension beats a keel.** With axial tension T, ride acceleration is V²f/T at every wavelength. 0.1–0.4 GN of CFRP rod (0.1–0.6 t/m) does what 0.3–4 t/m of keel does.
- **A big bore helps only as stroke.** Its wall gives little EI, and it reaches 65 % self-lift only at 80 m. An actuated guideway needs a 40–50 m bore on bare tethers, or ~20 m with 0.1–0.3 GN of line tension. Flying pods cost 155–270 Wh/pax-km and lack propulsion.
- **The best combination is only marginal.** Tension plus an actuated guideway in a 20 m bore reaches M1–4 at R ≈ 100–115 m and $105–170 M/km, but gives B/C 0.24–0.42 at base cost (0.9–1.45 at the low-cost end). Evacuation, breach, headway and line strain energy remain. The verdict is unchanged.

## System and economics ([system/system.md](system/system.md))

- **The time advantage comes from boarding at ground level, not from Mach.** NYC–LAX door to door: 3.6 h at M2, 2.7 h at M4, 2.5 h at M6, against 5.1 h for BoomShield M2.2 and 7.3 h subsonic. Lifts to a 20 km station or docking aircraft at altitude are worse.
- **Capacity:** a tube pair carries ~9 M passengers/yr at M4 with safe braking spacing; the spine needs ~81 M. Switching at speed is impractical.
- **Cost:** ~$66 M/km for the stratospheric run plus ~$157 M/km for ramps; ~$569 B for the spine plus LAX at M4. The docs' $5–10 M/km is far too low.
- **Energy is the one clear advantage:** 21–85 Wh per passenger-km of electricity against ~340 Wh of jet fuel, so net of fuel and CO₂ the case improves (unlike BoomShield).
- **B/C (2050 opening, 100 % capture), gross / net of fuel and CO₂:** base cost 0.36 / 0.47 vs subsonic and 0.28 / 0.61 vs Boomless; low cost 1.28 / 1.67 and 1.09 / 2.17. It pays only below ~$40–55 M/km all-in. A ground vactrain on the same network at ~$62 M/km scores higher (0.68 / 0.89 vs subsonic).
- **Hybrids with BoomShield:** none of the three studied works. The lesson to carry back: ground processing beats Mach.
