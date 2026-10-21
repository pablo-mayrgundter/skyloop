# Line ties vs smeared ties: independent check

Script: `check_lineties.py`. Data: `line_ties_check.csv`, `3d_cache/lineties_*.pkl`. Figure: `line_ties_check.png`. Tests: `tests/test_lineties.py`.

## Verdict: a mesh artefact

The 3-D finding that 25 m line ties at "K = 300" make the shield transparent is a mesh artefact, not physics.

- **Cause:** each line tie acts on one boundary element. The BEM faces are limp, so a tie is a spring between two face strips of mass m·b. Its coupling and its strip resonance, f_s ≈ (1/2π)·√(2 k_line / (m_eff b)), depend on b, which in the 3-D runs was the element size, λ⊥/10 ≈ 34/f m at Mach 2.
- **Effect:** at 25 m pitch (k_line = 7500 N/m²) f_s ≈ 7 Hz, inside 5–10 Hz; at 12.5 m ≈ 3.5 Hz; at 5 m several ties share an element (≈ smeared). This reproduces the circle-check ordering.
- **With a physical width (b = 0.5 m):** 25 m ties leak **less** than smeared ties in 5–10 Hz.

## Audit

- **Units, sign, panels: correct.** k_line [N/m per m of line] acts as the pressure k_line·e/L on the tie element; K = k_line/P matches PO and `run_3d.py`. The sign matches `GasLayer._tie_layer_scaled`, and the same-index top and bottom panels are coupled. The BEM's smeared and tie-free results match 1-D to 0.1 dB.
- **The defect:** the single-element application. On tension-free faces the answer moves with the mesh, by up to 30 dB at fixed f.

## Independent Floquet–Bloch model

The model is an infinite membrane / gas / membrane stack under the boom plane wave (k_x = ω/V, k_y = 0), with ties at pitch P. Each order n has k_y = 2πn/P and its own evanescent gas and air k_z. The tie force, uniform over b, couples all orders and is eliminated in closed form. Tension N and loss η enter through `shield.Membrane`.

Checks (in the tests): n = 0 alone reproduces `GasLayer(tie_K)` to 1e-15, and b = P gives the smeared result exactly. With N > 0, P → 0 recovers it too. With N = 0 it does not, unless b/P → 1: the correction −2K(P/b − 1)/(ω²m) does not depend on P.

**BEM vs Floquet** (H₂ 30 m, Mach 2, P = 25 m, K = 300; plane wave on a 1 km slab; |T₀|² in dB):

| f [Hz] | none | smeared | BEM, 3-D mesh (b = 34/f) / Floquet | BEM panel 3 m / Fl | panel 0.75 m / Fl | **fixed: b = 1.5 m, two meshes / Fl** | N = 1e4, b = 1.5 / Fl |
|---|---|---|---|---|---|---|---|
| 4 | −21.1 | −9.4 | −4.3 / −5.6 | −16.1 / −16.2 | −20.0 / −19.9 | −19.0, −18.9 / −18.7 | −10.7 / −10.9 |
| 6 | −30.4 | −17.4 | −2.2 / −4.5 | −18.0 / −18.3 | −27.3 / −27.3 | −25.4, −24.9 / −24.3 | −34.8 / −34.3 |
| 7 | −34.8 | −20.8 | **+0.9 / −0.8** | −16.9 / −17.5 | −30.1 / −30.1 | −27.0, −26.5 / −25.9 | −27.5 / −28.4 |
| 9 | −43.2 | −26.3 | −12.5 / −11.1 | **−2.4 / −5.4** | −33.9 / −33.9 | −26.6, −26.6 / −26.9 | −29.6 / −28.8 |
| 12 | −53.7 | −32.9 | −24.8 / −24.4 | −25.7 / −25.3 | −35.8 / −36.0 | −19.8, −19.7 / −21.2 | −31.1 / −29.5 |

Given the BEM's element width, Floquet matches the BEM to about 1 dB, resonance shift included. The tensioned column adds a finite-difference N·w_yy panel coupling (0.75 m panels): within 1.6 dB of Floquet at 2–12 Hz.

**Circle-check analogue**, mean |T₀|² over 5–10 Hz in dB (3-D circle values in brackets):

| no ties | line 25 m | 12.5 m | 5 m | smeared | 25 m, b 0.5, N 0 | 25 m, b 0.5, N 1e4 |
|---|---|---|---|---|---|---|
| −32.9 (−29.9) | −4.2 (−2.7) | −13.8 (−14.0) | −19.6 (−18.3) | −19.6 (−19.0) | −30.5 | −23.4 |

## Physics

- **Limp faces:** a narrow tie is not a stud that short-circuits the stack. The face cannot spread the force, so the tie moves only its strip and leaks mainly near f_s.
- **Tensioned faces (any pressurised envelope):** the spans act in series with the tie, K_eff ≈ K / (1 + K P² / 6N) statically, and span modes add narrow peaks (figure, middle). At 25 m and K = 300 the 5–30 Hz loss is 7.7 dB at N = 10⁴, against 13.0 dB smeared. At N = 10⁵, P ≤ 25 m acts smeared; at P = 50 m span modes enter the band (16 dB).

So the smeared model is not an underestimate. For realistic tension and P ≤ 25 m it is about right or slightly pessimistic.

## Design rule (≤ 3 dB loss of tie-free 5–30 Hz TL)

The criterion is band power weighted by 1/f², with b = 0.5 m and η = 0.05. The table gives the allowable K = k_line/P in Pa/m.

| stack | smeared | N = 10⁴, P = 5 / 12.5 / 25 / 50 m | N = 10⁵, same | N = 0 (limp), same |
|---|---|---|---|---|
| H₂ 30 m, M 2 | 150 | 159 / 144 / 162 / 81 | 149 / 154 / 179 / 90 | 53 / 24 / 12 / 6.3 |
| H₂ 30 m, M 1.4 | 40 | 41 / 44 / 43 / 30 | 39 / 39 / 40 / 30 | 28 / 18 / 11 / 6.0 |
| He 30 m, M 2 | 397 | 465 / 425 / 1700 / ∞ | 397 / 435 / 680 / 185 | 59 / 24 / 12 / 6.5 |
| He 30 m, M 1.4 | 94 | 97 / 99 / 100 / 62 | 90 / 92 / 100 / 59 | 43 / 21 / 11 / 5.9 |

∞ means the face span compliance alone keeps K_eff below the limit.

- **Tensioned faces (N ≥ 10⁴ N/m) at P ≤ 25 m:** K ≤ K_s, the smeared limit: 150 (H₂, M 2), 40 (H₂, M 1.4), 400 (He, M 2) and 90 (He, M 1.4) Pa/m. At P = 50 m, halve it.
- **Tendon lines:** for tendons of axial stiffness EA, length d = 30 m and spacing s_x, k_line = EA/(d·s_x) ≤ K_s·P. For H₂ at Mach 2, P = 25 m and s_x = 1 m: EA ≤ 1.1 × 10⁵ N per tendon.
- **Limp faces:** keep f_s below the band, k_line ≤ 2π² f_lo² m b (≈ 300 N/m² for b = 0.5 m), at any pitch.
- **Gate-0 K = 300:** it fails the rule both smeared and discrete, but it is not transparent.

## Fix (`boomsim/crossplane.py`)

- `line_tie_width`, or a per-tie `(y, k_line, b)`, gives each tie a physical width, resolved with nodes at y ± b/2; force uniform over the strip, elongation the strip mean. Results are mesh-independent (0.5 dB) and match Floquet (1 dB). The uniform mesh spreads by overlap and warns if b is below a panel.
- Without b: unchanged behaviour plus a RuntimeWarning. The PO path accepts 3-tuples.
- 12 new tests pass; the crossplane and shield suites still pass.

**Recommendation:** rerun the `ties300` and `ties30` cases with smeared ties (valid for N ≥ 10⁴ and P ≤ 25 m) or with `line_tie_width` set, and withdraw the −0.3 dB headline. Limp-face BEM ties understate the leak of tensioned faces.
