# Adversarial review: the evanescent gas-layer mechanism

**Verdict: holds, with caveats.** The math is right, and so is the physics for an ideal layer. The headline numbers (60 m → ≤1 % peak; −100 to −180 dB) assume a layer that is infinite, uniform and tie-free. Realistic geometry floors them at roughly −15 to −45 dB. A 15 m layer suffers less, but its value is unproven until it is measured as ground PLdB.

Scratch scripts: `/tmp/claude-0/-home-user-skyloop/85fc76a2-9709-5cfa-b462-9308def9aea4/scratchpad/review/`

## 1. Math (verified)
- An independent boundary-value solve matches `transfer()` to 0.01 dB (M 1.4–2, He/H₂, 15–60 m, 2–20 Hz, down to −180 dB).
- The branch choice, impedances and scaling are correct.
- A Floquet model with no walls reproduces the transfer-matrix result exactly.
- |T| falls well below exp(−κd), for example −28.5 vs −17.4 dB (15 m He, M 1.4, 10 Hz). The extra loss comes from impedance mismatch: the gas impedance is ρ_g/a ≈ 6 rayl against 37 rayl for air.

## 2. Physics (sound for the ideal case)
- This is the established sonic-boom-into-ocean problem (Sawyers 1968; Cheng & Lee, JFM 2004).
- The gain really is evanescent. At normal incidence the same He stack beats the mass law by only 0–11 dB (an ordinary double wall). At M 1.4–2 it gains 6–47 dB.
- A steady source has ω = V·k_x exactly, so nothing leaks through an infinite layer.
- Components with k_y ≠ 0 decay faster: κ is ×1.5 at 30° azimuth and ×2 at 45°.
- The lift-carrying k = 0 part must pass and does (T(0) = 1). At the ground it is ~W/(20 km)², i.e. mPa, so it is harmless.
- The real low-k leak is the passband: below ~2–5 Hz a 15 m layer is nearly transparent.
- Engine noise and manoeuvres have ω ≠ V·k_x, so the model does not cover them.

## 3. Caveats, ranked by how much they could erode the effect

**1. Ties or a stiff deck can erase the effect.** I modelled smeared tie stiffness K between the two faces. At 10 Hz, 6 dB of the benefit is lost at:

| Case | K [Pa/m] |
|---|---|
| 60 m He, M 2 | 35 |
| 30 m H₂, M 2 | 105 |
| 15 m He, M 1.4 | 275 |
| 15 m He, M 2 | 955 |

- The scale is set by ρ_g/(a²d).
- At K ≥ 10³–10⁴ Pa/m transmission returns to the mass law or worse (a mass–spring resonance).
- For reference, Dyneema tendons on a 10 m grid give ≈ 200 Pa/m. This is the building-acoustics stud short-circuit.
- Sparse line ties at 250 m pitch cost only 0–2 dB.

**2. Finite cross-track width sets a floor independent of opacity.** Two-edge diffraction, undertrack:

| W / h | 2 Hz | 10 Hz | 20 Hz |
|---|---|---|---|
| 250 m / 100 m | −8 dB | −15 dB | −18 dB |
| 250 m / 250 m | −6 dB | −12 dB | −15 dB |
| 1 km / 250 m | −15 dB | −22 dB | −25 dB |

**3. Cell walls and lens rims cap thick layers.** I modelled periodic cells whose partitions block lateral sloshing, with L = 250 m.
- 15 m He, M 2, 10 Hz: total transmitted power is −13 dB, against −18 dB ideal.
- 60 m He: power is floored at −33 dB (10 Hz) and −44 dB (20 Hz), against −67 and −119 dB ideal. At 2 Hz it is −3 dB against −14 dB, because sloshing needs L ≫ V/2f.
- In an elliptic lens the thin rim leaks. For 60 m, M 1.4, 10 Hz: −95 dB becomes −59 dB (coherent) or −45 dB (power).
- This is the wavy-ocean effect: any modulation scatters energy into the propagating cone.

**4. Gas purity near cutoff.**
- 5 % air drops c_He to 751 m/s and c_H₂ to 867 m/s.
- At M 2.5, He with 5 % air keeps only 36 % of κ; with 10 % air the layer propagates.
- Superheat (+20 K) helps; night cooling and wind (±5 % in V) matter only near cutoff.

**5. The metric, and re-steepening on the way down.**
- Peak at the shield overstates the benefit. The impulse ratio is 0.34 where the peak ratio is 0.054 (15 m He, M 1.4).
- Far-field N-wave strength scales as √(positive F-area), so a fully re-shocked wave would be only ~−5 dB.
- Over the 20 km descent, a transmitted 10–20 Pa peak advances 12–25 ms, against 19–32 ms of smoothing: partial re-steepening. Peaks of 50 Pa or more re-shock.
- The 0.3 s test wave is wrong at the shield: the real signature there is ~0.05–0.1 s with multiple shocks.

**6. Face coincidence.**
- If the faces' in-plane wave speed is ≈ 1.1–1.4 V, they resonate against the gas's negative compliance. At 10 Hz that gives −7 dB instead of −18 dB.
- Limp 5 kN/m membranes (c_m ≈ 95 m/s) are safe. Stiff sandwich panels may not be.

**7. Linearity is adequate.**
- p/γp ≤ 0.05.
- Lateral gas velocity is ≈ 20–40 m/s (u/V ≤ 0.1). Displacements are 0.1–0.3 m.
- However, the boom is 10–30× the ≈ 11 Pa superpressure, so the faces go slack and wrinkle on every pass. That is a fatigue issue and adds some scatter.

**8. Reflection.**
- |R| ≈ 1, and the reflected cone reaches flight level 200–700 m behind the nose. That affects in-trail traffic and doubles the load on the top face.

## 4. Next fidelity step
1. **Move the 2-D cell model into `boomsim`.** Promote the scratch periodic-cell model, with walls, rims, ties, tension and loss. Add a test that the no-wall case reproduces the transfer-matrix result.
2. **Model the corridor in 3-D.** Integrate over k_y with a finite-width corridor, and compute both undertrack and lateral ground signatures.
3. **Score at the ground.** Feed the near-field signature through the shield, then Burgers propagation to the ground. Score PLdB and impulse, not peak.
4. **Adopt design rules:**
   - K ≪ ρ_g/(a²d)
   - L ≫ V/2f_min
   - purity ≥ 98 %
   - limp faces
   - W ≳ 4h
5. **Report gains against fair baselines** (a double wall; the air@0.1p case), not the single-membrane mass law alone.
