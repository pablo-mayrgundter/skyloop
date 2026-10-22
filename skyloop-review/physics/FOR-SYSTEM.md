# Physics -> system handoff (numbers you can reuse; source CSVs in this folder)
Conventions: V = M x 295.1 m/s (20 km c); pod 3,000 kg, 3 pax, 1.0 m^2 frontal; tube gas 250 K.
- Accel/decel length (propulsion.csv): 0.3 g -> M2 59 km / 201 s, M4 237 km / 401 s, M6 533 km / 602 s.
  0.5 g -> 36 / 142 / 320 km. Peak power per pod at 0.3 g: 5.2 / 10.4 / 15.6 MW.
- KE per pod 0.52 / 2.09 / 4.70 GJ (145 / 580 / 1306 kWh); per pax 48 / 194 / 435 kWh.
  Station average launch power at docs' 4,000 pods/h: 0.58 / 2.3 / 5.2 GW per direction (gross, before regen).
- Energy per pax-km (energy_per_paxkm.csv; regen 0.8, EDS L/D 1000, 1 Pa): M4 1,000 km 42, 3,000 km 16, 9,000 km 7 Wh;
  M6 3,000 km 32, 9,000 km 13 Wh; M2 1,000 km 12 Wh. Docs claim 0.5-2.
- SAFETY HEADWAY (ride_geometry.csv): brick-wall stopping at 1 g: M2 18 km, M4 71 km, M6 160 km ->
  max 120 / 60 / 40 pods/h = 360 / 180 / 120 pax/h per tube at 3 pax/pod. At 0.5 g half that.
  Docs' 4,000 pods/h needs ~0.9 s headway (relative braking only; a pod failure then collides with many).
- Tube + lift mass at 20 km, D 3 m (tube_mass_budget.csv, envelope_sensitivity.csv): payload 125-315 kg/m
  (central 184: wall 35, guideway 60, LSM avg 71 for M4/1000 km, util 18); H2 envelope R 43-82 m, lift 0.5-1.8 t/m,
  H2 37-131 t/km, tether 45 MN / 0.9 kt every 3 km at baseline (C_D 0.6, V/H 1.5).
- Route geometry (ride_geometry.csv): horizontal R at 0.1 g: 355 / 1,421 / 3,196 km; banked 0.5 g: 71 / 284 / 639 km.
- Pumping (pumping_per_km.csv): 1 cm^2 holes/km at 20 km, 10 Pa: 6.5 m^3/s, ~65 kW/km; at ground 16x more.
