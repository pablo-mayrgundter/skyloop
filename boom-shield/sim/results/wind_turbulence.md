# Wind and turbulence on the 3-D BoomShield carpets

Concorde M2, h = 250 m. Cases: unshielded, W = 1 km H₂ 30 m, W = 500 m He 30 m. All use the cached R_m = 3 km circle signatures from `run_3d.py` (no cross-plane re-run). Driver: `run_wind_turbulence.py`. Data and figure: `wind_turbulence.csv`, `.png`.

Δ = shielded − unshielded PLdB in the same atmosphere. "Headline" is the less favourable of the raw and Gate-0-proxy Δ, as in the 3-D study.

## Does either change the 3-D conclusions? No.

1. **Wind.**
   - In the same wind, undertrack Δ moves ≤ 0.4 dB, Δ at ±24° ≤ 0.3 dB, and the carpet maximum −2.8 to +1.9 dB. Width still sets the result.
   - The carpet stretches, shrinks or shifts by km to tens of km, and the shadow moves with it. Corridor siting therefore needs wind climatology.
   - The station-kept shield sees ±15–25 m/s of trace-speed change. This alters only the inaudible pulse below 5 Hz: ground PLdB changes ≤ 0.01 dB, even at ±60 m/s.
   - No cross-plane re-runs are needed for along-track wind.
2. **Turbulence.**
   - Median PLdB stays within about 1 dB of the deterministic raw value (P10–P90 about ±3 dB). The Gate-0 proxy sits 2.7–10.7 dB below raw, outside that spread.
   - The paired turbulent Δ is close to the *raw* Δ: −46.8 dB undertrack at 1 km, against a headline of −38.5. The 3-D headline was conservative by 5.6–8.3 dB (1 km, all φ; 500 m, φ ≤ 24°), and by about 1 dB at the 500 m carpet edge.
   - The R_m caveat (8 km is about +5 dB near track) still applies.
3. **Leak pulses vs N-wave.**
   - The 1 km leak pulses (3–8 Pa) are not shocks: they rise over milliseconds.
   - Turbulence barely spikes them: > 20 Hz peak ×0.95 vs ×1.38 for the N-wave. Their max slope drops 11 % vs 69 %.
   - Their PLdB scatters as much as the N-wave's (the multipath coda adds high frequencies to both), but in step with it, so the paired Δ spread is only ±1.5 dB.

## Methods

**Wind** (`propagate.py`; optional `wind=` on RayPath, propagate, ground_signature, propagate_from_circle, ground_carpet, lateral_cutoff).
- Medium and aircraft: stratified w = (u(z), v(z)); the aircraft crabs to hold its track at a fixed air Mach number, V_g = √(V² − v_a²) + u_a.
- Rays: normals on the Mach cone about the air velocity, s = n/(c_a + w·n), conserved s_h (s_x = 1/V_g), Ω = 1 − w·s_h, s_z² = Ω²/c² − s_h², ray velocity c n + w.
- Amplitude: Blokhintzev invariant p²|v_ray|A/(ρc²Ω), i.e. p² s_z ∂Y/∂φ / (ρΩ²).
- Propagation: ground-frame retarded time. Burgers over ds_e = Ωc dt with local values; absorption at Ωf. Air-frame start signatures have t × Ω.
- Default: `wind=None` leaves outputs bit-identical.

**Wind validation** (`tests/test_wind.py`, 40 tests).
- Uniform wind is an exact Galilean shift: geometry; the linear absorbing march with t × Ω to 10⁻⁹; nonlinear < 1 %.
- Jacobian vs neighbouring rays: 10⁻⁶. Pierce's form from traced tubes: 10⁻⁵. Calm vs windless: 10⁻⁷.

**Wind profiles.**
- Winter jet: 65 m/s at 10.5–12 km, 15 m/s at ≥ 20 km. Sources: SKYbrary "Jet Stream" (typical 50 m/s, up to 100 m/s); Holton & Hakim 2013.
- QBO easterly: 25 m/s at 20 km, 30 m/s at 26 km. Sources: Baldwin et al., *Rev. Geophys.* 39 (2001); Pahlavan et al., QJRMS 2021.

**Shield frame.** `shield.transfer` was evaluated at V_g; only the gas-transmitted part changes, P += P_inc·(T(V_g) − T(V)).

**Turbulence** (`turbulence.py`).
- Model: von Kármán phase screen for a 1 km planetary boundary layer (PBL), placed 500 m above the ground.
- Presets: weak / moderate / strong = σ_T 0.2 / 0.5 / 1 K and σ_w 0.3 / 1 / 1.7 m/s.
- Propagation: Fresnel–Kirchhoff to the ground. Impulse response h = A''(t)/(2πcd), where A(t) is the screen area with delay < t, computed exactly on a 0.5 m triangulation. This captures spikes (caustics), rounding (multipath) and diffraction.
- Tests: lens focal gain to 0.3 %, a tilt is a pure delay, identity without turbulence.
- Monte-Carlo: 200 / 100 / 100 kernels, the same for shielded and unshielded.
  - (A) kernel applied at the ground;
  - (B) kernel applied at 500 m, then the nonlinear march down (re-steepening).

## Results: wind

| case | V_g m/s | cut-off +y (φ, km) | cut-off −y | y_g at 24° / 48° / −48°, km |
|---|---|---|---|---|
| calm | 590 | 54.9°, 47.7 | −47.7 | 9.7 / 26.9 / −26.9 |
| jet headwind | 575 | 52.5°, 34.9 | −34.9 | 9.2 / 25.0 / −25.0 |
| jet tailwind | 605 | 57.1°, 88.1* | −88.1* | 10.3 / 30.5 / −30.5 |
| jet crosswind (+y) | 590 | 55.9°, 114.7* | −51.6°, −28.4 | 12.8 / 36.8 / −21.1 |
| QBO headwind | 565 | 50.5°, 41.7 | −41.7 | 10.2 / 31.6 / −31.6 |

\*The grazing ray is sensitive to the lowest kilometre: with a uniform 10 m/s boundary layer the tailwind cut-off is 63 km. Unshielded PLdB there is below 86. Crosswind moves the undertrack ray, and the shadow, 2.4 km downwind.

| ΔPLdB (headline) | calm | jet head | jet tail | jet cross | QBO head |
|---|---|---|---|---|---|
| **1 km**, undertrack | −38.5 | −38.6 | −38.3 | −38.5 | −38.9 |
| 1 km, max over φ ≤ 48° | −25.9 | −24.5 | −28.7 | −23.9 | −27.4 |
| 1 km, φ = ±48° | −20.3 | −20.1 | −21.1 | −21.5 / −20.1 | −21.3 |
| **500 m**, undertrack | −22.5 | −22.3 | −22.5 | −22.5 | −22.6 |
| 500 m, max (φ = 0, 16, 24, 48°) | −9.0 | −7.5 | −10.9 | −7.1 | −9.5 |
| 500 m, φ = ±48° | −3.4 | −3.1 | −3.3 | −3.4 / −3.3 | −3.3 |
| unshielded undertrack PLdB | 108.4 | 108.4 | 108.3 | 108.4 | 108.0 |

**Shield frame.** The evanescence margin c_gas/V_g is 1.72–2.11 (H₂) and 1.33–1.63 (He) for V_g = 530–650 m/s. |T|² in dB at V_g = V − 15 / V / V + 15 m/s:

| band | H₂ 30 m | He 30 m |
|---|---|---|
| 0.2–2 Hz | −4.2 / −4.0 / −3.8 | −2.0 / −1.9 / −1.7 |
| 2–5 Hz | −16.9 / −16.3 / −15.8 | −10.6 / −10.1 / −9.6 |
| 5–10 Hz | −33.8 / −33.1 / −32.3 | −25.3 / −24.4 / −23.5 |
| 10–20 Hz | −57.4 / −56.2 / −55.1 | −46.5 / −45.1 / −43.8 |

| ground, φ = 0 | dV −60 | −15 | +15 | +60 |
|---|---|---|---|---|
| 1 km: transmitted-pulse peak, Pa (calm 11.5) | 8.3 | 10.7 | 12.4 | 15.2 |
| 1 km: ΔPLdB / ΔCSEL | 0.00 / −0.19 | 0.00 / −0.05 | 0.00 / +0.06 | 0.00 / +0.29 |
| 500 m: ΔPLdB (peak 25.2 → 23.4 / 27.0 Pa at ∓15) | – | −0.00 | +0.01 | – |

Scaling the whole field by T(V_g)/T(V) gives ±8 dB, but that is unphysical: the audible edge field never crosses the gas.

A crosswind (crossflow Mach 0.05–0.08) would make the two edge waves asymmetric, which only a convected cross-plane solver can resolve. I estimate the effect at about 1 dB.

## Results: turbulence (moderate; Δ paired)

| config, φ | det. raw / proxy | PLdB A: median [P10, P90] | B | Δ det. raw / proxy | Δ A | Δ B |
|---|---|---|---|---|---|---|
| unshielded 0° | 108.4 / 97.7 | 107.9 [105.4, 111.2] | 107.4 [103.9, 111.5] | – | – | – |
| unshielded 48° | 100.6 / 92.1 | 100.2 [97.9, 103.5] | 99.6 [96.1, 103.5] | – | – | – |
| 1 km 0° | 61.9 / 59.2 | 60.9 [58.0, 64.5] | 61.1 [58.0, 64.4] | −46.5 / −38.5 | −46.8 [−48.4, −45.3] | −46.2 [−48.3, −44.4] |
| 1 km 24° | 71.1 / 67.1 | 70.2 [67.2, 73.6] | 70.2 [67.3, 73.5] | −35.7 / −29.4 | −36.2 [−37.4, −35.0] | −35.6 [−37.4, −34.0] |
| 1 km 48° | 75.2 / 71.8 | 74.2 [71.5, 77.6] | 74.3 [71.4, 77.4] | −25.5 / −20.3 | −25.9 [−27.1, −24.9] | −25.3 [−26.9, −24.1] |
| 500 m 0° | 80.0 / 75.2 | 79.1 [76.7, 82.5] | 79.2 [76.6, 82.4] | −28.4 / −22.5 | −28.6 [−29.5, −27.8] | −28.2 [−29.7, −26.7] |
| 500 m 48° | 96.4 / 88.7 | 95.8 [93.5, 99.0] | 95.3 [92.2, 99.1] | −4.3 / −3.4 | −4.4 [−4.5, −4.2] | −4.3 [−4.8, −3.9] |

Weak and strong presets give Δ within 0.9 dB of these.

| moderate A, ratio to deterministic | N-wave, φ = 0 | 1 km leak, φ = 0 | 500 m edge shock (25 Pa), φ = 48° |
|---|---|---|---|
| peak | 1.27 [1.04, 1.64] | 1.00 | 1.11 [0.98, 1.33] |
| > 20 Hz peak | 1.38 [0.90, 2.16] | 0.95 [0.80, 1.13] | 0.98 [0.81, 1.45] |
| max slope | 0.31 | 0.89 | 0.50 |

| unshielded calibration (weak → moderate) | model | White Sands, low → moderate turbulence (Willshire & DeVilbiss, NASA CP-3172, 1992) |
|---|---|---|
| mean rise time | 1.3 → 3.1 ms (×2.4) | ×4 |
| rise-time std | ×2.2 | ×3 |
| peak std | ×1.9 | +40 % |

The moderate median rise is 2.2 ms [0.7, 6.9], inside the 2–5 ms typical of measured booms. Peaks are biased high (mean ×1.31 in A, ×1.17 in B), whereas flight data are roughly log-normal about 1 (Garrick & Maglieri, NASA TN D-4588, 1968). The model's single screen without saturation over-spikes.

**Proxy.** The 1 ms Gaussian proxy matches the turbulent rise time (2.6 ms), but it deletes everything above about 150 Hz. Turbulent booms keep that energy in spikes and a scattered coda (figure, panel f). So the proxy quietens the N-wave by 10.7 dB but the leak by only 2.7 dB, and the headline credits the unshielded boom with a quietening it does not get.

**Caveats.**
- Calm air is assumed inside R_m.
- One screen at normal incidence is used for all φ.
- Spiking is probably overstated.
- The paired Δ assumes one atmosphere for both booms; with independent atmospheres the spread would be about √2 × 2.3 dB.
- Cache: `$WT_CACHE` (default `results/wind_turbulence_cache/`).
