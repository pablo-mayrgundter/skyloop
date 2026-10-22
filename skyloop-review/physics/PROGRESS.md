# physics review — working notes (agent scratch, keep until physics.md final)
- 2026-10-22: read README, docs 01-08, BoomShield vacuum/altitude/tethers/structure/STATUS.
- Plan: skyloop_physics.py -> CSVs (tube_mass, drag_vs_pressure, kantrowitz, accel, curves, alignment, pumping) + physics.png; then physics.md.
- Key early numbers (hand calcs, to be confirmed by script):
  - tube own vacuum lift at 20 km: 0.28-2.5 kg/m for D 2-6 m vs 100-400 kg/m tube+track -> no self-buoyancy.
  - accel length at 0.3 g: M2 59 km, M4 237 km, M6 533 km.
  - alignment: allowable amplitude A = a_max (lambda/2 pi v)^2 ~ 3.6 cm at 3 km, M6, 0.05 g; tether-top gust motion ~20 m.
  - brick-wall braking at 1 g: M6 160 km, 180 s.
- Script skyloop_physics.py runs; CSVs written (tube_mass_budget, envelope_sensitivity, invballoon_check, kantrowitz,
  drag_vs_pressure, pumping_per_km, propulsion, energy_per_paxkm, ride_geometry, alignment_allowable,
  floating_tube_motion, thermal_bow). Key: H2 envelope R 82 m (baseline, holddown dominates 1485 of 1756 kg/m lift);
  C_D0.3+V/H0.75 -> R 43 m; untethered 30 m. Gust motion 8-41 m tether top, allowable 0.14 m (M6) / 1.3 m (M2) at 6 km.
  Brick-wall capacity M6 1 g: 40 pods/h = 120 pax/h. Energy M4 3000 km 8-29 Wh/paxkm.
- Next: make_figure.py -> physics.png; physics.md.
