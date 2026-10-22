# Focused sonic booms (superbooms) and BoomShield — DRAFT

*Draft written while the runs complete; numbers marked TBD are pending.*

Code: `boomsim/focus.py`, `run_focus.py`; tests `tests/test_focus.py`. Data: `focus_booms.csv`; figure `focus_booms.png`.

## Bottom line (preliminary)

- Every acceleration through the Mach cut-off focuses on the ground. Overture-like, 20.25 km, 0.02–0.15 g: focus at emission Mach 1.16–1.28, **~4 psf unshielded** (3–4× the steady boom at that Mach), 109–110 PLdB.
- Over the shield (H₂ 55 m) the same focus is **~0.01 psf**: the layer removes the shock before the rays converge, and the remaining infrasonic swell is too long to focus.
- Deceleration forms no caustic. Turns at Mach 1.7 focus only above ~15–20° of bank, at the inside carpet edge.
