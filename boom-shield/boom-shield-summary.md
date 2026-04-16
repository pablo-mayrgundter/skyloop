# Boom Shield Acoustic Analysis Summary

This directory captures the first-pass acoustic analysis for the Skyloop boom-shield concept.

## Scope

The goal was to test whether a buoyant horizontal shield might suppress the **direct sonic-boom wave** enough to make faster overland supersonic or hypersonic corridors plausible, starting from NASA X-59 low-boom performance as a calibration point.

## Main assumptions

- Calibrate the open-air baseline to **X-59** at about **Mach 1.4**, **55 kft / 16 km**, and about **75 PLdB** on the ground undertrack.
- Preserve the geometric similarity ratio between the open-air reference and the shielded analog:
  - open-air reference: **W = 64 km**, **E = 16 km**
  - shield analog: **W = 1 km**, **E = 250 m**
  - so **W/E = 4** and **(W/2)/E = 2** in both cases.
- First focus only on the **central direct wave**. Side-lobes, edge diffraction, reradiation, and structural motion are acknowledged but deferred or simplified.

## Key variables

- **W**: shield width
- **E**: aircraft height above shield
- **A_req**: required transmitted direct-wave pressure amplitude fraction
- **Π_req**: required shield opacity in dB, with `Π_req = -20 log10(A_req)`
- **A_mat**: intrinsic direct-wave transmission floor of the shield
- **λ_eff**: effective acoustic wavelength controlling diffraction leakage
- **E_max,acoustic**: max spacing still allowed by the direct-wave suppression model
- **E_min,ops**: minimum spacing required by flight-control / shield-motion safety margins
- **Γ = E_max,acoustic / E_min,ops**: viability margin ratio

## What the analysis found

### 1. Direct-wave suppression framing works well

The useful design question is not the absolute boom from scratch, but:

> how much direct-wave attenuation must the shield provide, relative to the X-59 low-boom baseline?

This led to a compact surrogate model for the open-case boom growth with Mach and altitude, anchored at X-59.

### 2. SR-71-class looks plausible on the direct-wave-only model

Using a simple scaling law anchored at X-59:

- **SR-71-like** at **Mach 3.2 / 85 kft** needs about:
  - **7.9 dB** opacity in an optimistic shape case (`Cs = 1`), so `A_req ≈ 0.40`
  - **14 dB** opacity in a rougher source-shape case (`Cs = 2`), so `A_req ≈ 0.20`

This suggests the boom shield concept may have a real design window for SR-71-class overland operation.

### 3. SR-72-class may be plausible, but only in the optimistic regime

For a notional **SR-72-like** case at **Mach 6 / 90 kft**:

- optimistic source-shape case (`Cs = 1`):
  - required opacity about **15.8 dB**
  - `A_req ≈ 0.16`
- rougher source-shape case (`Cs = 2`):
  - required opacity about **21.8 dB**
  - `A_req ≈ 0.08`
- still rougher (`Cs = 3`):
  - required opacity about **25.3 dB**
  - `A_req ≈ 0.054`

So Mach-6 overland looks sensitive to source shaping, diffraction wavelength, and shield transmission floor.

### 4. Width helps, but spacing E is the more sensitive parameter once width is already large

After the shield width is already “big enough” to intercept the main lobe, the bigger question becomes the aircraft height above the shield.

This matches the design intuition that **250 m spacing** may matter more than pushing width beyond 1 km.

### 5. A Fresnel / knife-edge surrogate gave a better spacing model

Replacing a crude leakage law with a diffraction-style surrogate:

- `A_edge ~ sqrt(λ_eff E) / W`

made the trade much clearer.

With **W = 1 km**, **A_mat = 0.05**, and **λ_eff = 30 m**:

- **SR-71-like optimistic:** `E_max ≈ 5.1 km`
- **SR-72-like optimistic:** `E_max ≈ 526 m`
- **SR-72-like with shape penalty:** `E_max ≈ 41 m`

That suggests:

- **250 m** spacing looks very safe for SR-71-like direct-wave suppression
- **250 m** spacing looks plausible for optimistic SR-72-like cases
- **250 m** may fail if the Mach-6 source signature is materially worse than the X-59-like optimistic assumption

### 6. The real master inequality is a window problem

The concept becomes viable only if:

`E_min,ops < E_max,acoustic`

where:

- `E_max,acoustic` comes from the direct-wave shielding model
- `E_min,ops` comes from safety clearance, aircraft motion, shield motion, sag, navigation error, and gust response

That means the next practical design thread should focus on material, buoyancy, tethering, structural modes, and control stability, because those govern both the transmission floor `A_mat` and the operational spacing floor `E_min,ops`.

## Cautions

- These are **surrogate models**, not high-fidelity sonic-boom propagation claims.
- Side-lobes were not solved here.
- Edge diffraction and reradiation are likely the next acoustic bottlenecks after direct-wave suppression.
- The key unknowns going forward are:
  - source shaping quality at higher Mach
  - effective acoustic wavelength `λ_eff`
  - intrinsic shield transmission `A_mat`
  - shield vertical motion / stability

## Suggested next thread

Proceed to the shield design itself:

- material architecture
- buoyancy and mass-per-area
- rigid segment spacing
- tethering and crosswind response
- structural modes / heave / sag
- operational clearance and control margins

Those are the ingredients needed to evaluate whether the spacing window stays open for SR-71-class and SR-72-class use cases.
