# Focused sonic booms (superbooms) and BoomShield

*DRAFT, numbers being finalised.* Code: `boomsim/focus.py`, `run_focus.py`; tests: `tests/test_focus.py`. Data: `focus_booms.csv`; figure: `focus_booms.png`.

## Physics

An aircraft in steady flight lays down a ray pattern that does not change. If it accelerates, turns or pushes over, rays emitted at different times cross. Their envelope is a **caustic**, and where it touches the ground a **focus boom** lands. Here geometric acoustics fails, because the ray-tube area goes to zero.

During acceleration the Mach cone narrows, so later rays are steeper and overtake earlier ones. In a real atmosphere sound is faster near the ground, so the first rays that reach it (Mach cut-off, M_cut = 1.153 at 20.25 km) arrive grazing. Every acceleration through M_cut therefore puts a crescent-shaped focus line on the ground, ahead of the aircraft. A gentle acceleration does not avoid it; it only moves the focus Mach closer to M_cut.
