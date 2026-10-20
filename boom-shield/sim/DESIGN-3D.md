# boomsim 3-D — finite-width corridor model (design contract)

Gate 0 (`DESIGN.md`) used a 1-D undertrack ray and a horizontally infinite shield. The ideal
evanescent gas layers drove ground PLdB below the Mark VII contours, so the realistic answer is
set by what that model left out: finite corridor width (edge diffraction), cell walls and lens
rims, ties, and off-track rays. This model adds those and scores the whole ground carpet.

## Key reduction: steady flight ⇒ a family of 2-D cross-plane problems

Aircraft in steady level flight along x at speed V; shield invariant along x (v1). In the
ground frame every field ∝ e^{j(ωt − k_x x)} with k_x = ω/V, so for each frequency ω the
problem is 2-D in the cross-plane (y lateral, z vertical):

- air:  ∇²_{yz} p + k_⊥² p = 0,  k_⊥ = ω·sqrt(1/c² − 1/V²) = k_x·β  (propagating)
- gas with c_g > V:  ∇²_{yz} p − κ² p = 0,  κ = ω·sqrt(1/V² − 1/c_g²)  (modified Helmholtz)
- limp membrane (areal mass m) between media: normal velocity continuous,
  p_above − p_below = jωm·v_n (tension terms optional); normal velocity
  v_n = −(1/(jωρ)) ∂p/∂n on each side.
- smeared ties K [Pa/m] between top and bottom faces at the same y (optional, local).
- incident field: exact linear supersonic slender-body solution in the cross-plane — volume
  source ∝ H₀⁽²⁾(k_⊥ r), lift dipole ∝ cos φ · H₁⁽²⁾(k_⊥ r) — with strengths from the
  equivalent-area distributions in `source.py`, normalised so its large-r asymptote equals
  Whitham's near field (`source.near_field`) at φ = 0.

Same time/FFT conventions as `DESIGN.md` (rfft, physical e^{+jωt}).

## Solver (crossplane.py)

Hybrid in frequency, with an explicit overlap check:
- f ≤ f_split: 2-D multi-domain boundary element method (exterior air with Hankel kernel;
  each gas cell its own domain with K₀ kernel; membranes and partitions as thin coupling
  sheets; slab or lens cross-section; internal partitions; ties).
- f > f_split: Kirchhoff/physical optics + angular spectrum: infinite-stack transmission
  T(k_x, k_y) from `shield.transfer` inside the strip, free field outside, propagated with
  exact 2-D angular spectrum.
- Must agree with each other near f_split, and converge in f_max (near-field band limit).

Output: time signatures on a matching circle of radius R_m centred on the flight axis,
for a set of downward azimuths φ (φ = 0 undertrack), in uniform air at the shield
altitude's ambient values. Unshielded runs go through the same path.

## Off-track propagation (propagate.py / source.py extensions)

- `RayPath(altitude, V, phi=0.0, ...)`: off-track ray (no wind), horizontal slowness
  (1/V, s_y) fixed by the launch azimuth; ray-tube area from neighbouring rays; lateral
  cutoff detection; ground lateral offset y_g(φ).
- Start from a point on the ray at slant distance R from the axis, applying the WKB density
  correction for the uniform-air assumption used inside R.
- Source with azimuth: equivalent area A_e(x, φ) = A_v(x) + A_L(x)·cos φ.
- `ground_carpet(...)`: ground signatures and metrics vs y_g for a list of φ.

## Scoring

- PLdB (raw and rise-time proxy, headline = less favourable, as in Gate 0) and CSEL vs y_g.
- Carpet metrics: undertrack PLdB, maximum PLdB over the carpet, carpet half-width.
- Compare against unshielded with identical numerics.

## Known exclusions of v1

Shield variation along x (module joints, transverse cell walls), wind, turbulence,
aircraft manoeuvre, structure-borne paths other than smeared/line ties, nonlinearity inside
R_m (quantified by an unshielded linear-vs-nonlinear check).
