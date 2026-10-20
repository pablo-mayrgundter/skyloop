# boomsim — design contract

First-cut 1D simulation of a sonic boom passing through a BoomShield and on to the ground.
Goal: answer "Gate 0" — can any buoyant shield stack that closes its own lift budget
reduce *ground* loudness (PLdB) meaningfully, once 20 km of nonlinear propagation is included?

Scope of v0: undertrack ray only (φ = 0), steady level flight, no wind, horizontally
infinite shield (no edge diffraction), linear shield interaction, nonlinear atmosphere.

## Conventions (all modules)

- SI units everywhere. Altitudes `z` in metres above ground.
- A signal is a pair of 1-D numpy arrays `(t, p)`: uniformly sampled time `t` [s] and
  overpressure `p` [Pa]. Sample rate `fs = 1/(t[1]-t[0])`. Functions that change
  the record length return a new `t`.
- Frequency domain uses `numpy.fft.rfft`/`irfft` (forward kernel e^{-jωt}); physical
  harmonic fields are written e^{+j(ωt − k·x)}. A pure delay Δ is `X·exp(-jωΔ)`.
- Aircraft speed `V` [m/s] sets the horizontal trace velocity of the boom on any
  horizontal surface (horizontal slowness 1/V is conserved along the ray, no wind).

## Package layout (`boom-shield/sim/`)

```
boomsim/
  atmosphere.py   US Standard Atmosphere 1976 + humidity profile
  source.py       Aircraft definition, equivalent area, Whitham F-function, near-field signature
  propagate.py    Ray-tube geometry + augmented Burgers propagation (nonlinear, absorbing)
  shield.py       Layered transfer-matrix shield model + lift budget
  loudness.py     Stevens Mark VII PLdB + other waveform metrics
tests/            pytest; run with `python -m pytest boom-shield/sim/tests`
results/          generated CSV/PNG/MD from run scripts
run_initial.py    Gate-0 study driver (written after modules land)
```

## Module interfaces

### atmosphere.py
- `std_atm(z) -> Atm` (dataclass, vectorised OK): `T` [K], `p` [Pa], `rho` [kg/m³], `c` [m/s], `z`.
- `molar_humidity(z) -> h` [% molar concentration of water vapour], from an assumed
  relative-humidity profile (documented in the docstring, overridable).

### source.py
- `@dataclass Aircraft(name, length, mach, altitude, weight, volume, lift_start=0.4, lift_end=1.0)`
  Weight [N] supported by lift; volume [m³] distributed as a Sears–Haack body.
  Property `V` = mach × c(altitude).
- `equivalent_area(ac, x) -> A_e(x)` volume + lift contribution
  `A_L(x) = β/(ρV²) ∫₀ˣ l(ξ) dξ`, l(ξ) smooth and zero at both ends of the lifting segment.
- `f_function(ac, y) -> F(y)` Whitham F-function `(1/2π)∫ A_e''(ξ)/√(y−ξ) dξ`.
- `near_field(ac, r0, fs) -> (t, p)` linear signature at distance r0 below the axis:
  `Δp = p∞ γ M² F(y) / √(2βr0)`, `t = y/V`, zero-padded.

### propagate.py
- `class RayPath(altitude, V)`: undertrack ray from an aircraft at `altitude` with speed `V`.
  Gives arc length `s(z)`, ray-tube area `A(z)` (up to a constant), local Atm.
- `propagate(t, p, ray, z_from, z_to, **opts) -> (t, p)`: augmented Burgers from altitude
  `z_from` down to `z_to` (z_to < z_from). Linear-only and lossless switches in opts.
- `ground_signature(ac, r0=None, insert=None, fs=..., reflection=1.9) -> (t, p)`:
  near field at r0 → propagate to ground. If `insert=(z_ins, fn)`, propagate to `z_ins`,
  apply `p = fn(t, p)` there, then continue to the ground. Ground reflection factor applied last.

### shield.py
- Gas helpers: `gas(kind, T, p)` for `"air"`, `"helium"`, `"hydrogen"`; partial-pressure air allowed.
- Elements: `Membrane(m, N=0.0, eta=0.0)`, `GasLayer(gas, d)`, `Screen(R)` (flow-resistive sheet).
- `Stack(layers)`; ambient air both sides.
- `transfer(stack, f, V, atm) -> (T, R)` complex pressure transmission/reflection at
  frequencies `f` for a wave with horizontal trace velocity `V`.
- `apply(stack, t, p, V, atm) -> p_out` transmitted waveform in the air below.
- `lift_margin(stack, atm) -> σ_net` [kg/m²]: buoyant lift of gas layers minus all element
  masses (positive = floats with that much spare per m²).

### loudness.py
- `pldb(t, p) -> float` Stevens Mark VII perceived level (Shepherd & Sullivan, NASA TP-3134).
- `metrics(t, p) -> dict`: `peak_pa`, `peak_psf`, `rise_time_ms`, `duration_ms`,
  `pos_impulse_pa_s`, `pldb`, plus `asel`/`csel` if implemented.

## Validation anchors
- Concorde cruise (~Mach 2, ~16–18 km): ground N-wave ≈ 1.9–2.0 psf (~93 Pa), ~300 ms.
- Shield model: single limp membrane → mass law; lossless stack conserves energy;
  evanescent gas layer → tunnelling ~ exp(−κd).
- PLdB: reproduce published Mark VII values for reference N-waves to ±0.5 dB.
