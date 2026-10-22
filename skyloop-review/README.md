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

The original docs' claims that fail: the invballoon is not "entirely in tension" and its lift is wrong; the tube's first mode is ~4 Hz, not 10; 10⁻³ atm is far too high a pressure; "gliders" cannot glide in vacuum, so the whole route needs a levitation track.

## System and economics ([system/system.md](system/system.md))

- **The time advantage comes from boarding at ground level, not from Mach.** NYC–LAX door to door: 3.6 h at M2, 2.7 h at M4, 2.5 h at M6, against 5.1 h for BoomShield M2.2 and 7.3 h subsonic. Lifts to a 20 km station or docking aircraft at altitude are worse.
- **Capacity:** a tube pair carries ~9 M passengers/yr at M4 with safe braking spacing; the spine needs ~81 M. Switching at speed is impractical.
- **Cost:** ~$66 M/km for the stratospheric run plus ~$157 M/km for ramps; ~$569 B for the spine plus LAX at M4. The docs' $5–10 M/km is far too low.
- **Energy is the one clear advantage:** 21–85 Wh per passenger-km of electricity against ~340 Wh of jet fuel, so net of fuel and CO₂ the case improves (unlike BoomShield).
- **B/C (2050 opening, 100 % capture), gross / net of fuel and CO₂:** base cost 0.36 / 0.47 vs subsonic and 0.28 / 0.61 vs Boomless; low cost 1.28 / 1.67 and 1.09 / 2.17. It pays only below ~$40–55 M/km all-in. A ground vactrain on the same network at ~$62 M/km scores higher (0.68 / 0.89 vs subsonic).
- **Hybrids with BoomShield:** none of the three studied works. The lesson to carry back: ground processing beats Mach.
