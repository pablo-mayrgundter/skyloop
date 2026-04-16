# Boom Shield

This directory collects the first-pass acoustic analysis for a bouyant sonic-boom shield concept, to reflect the main wave of sonic booms from overflying aircraft away from their ground tracks.

## Contents

Copied from the full [ChatG thread here](https://chatgpt.com/share/69e16de1-e3ac-83e8-9bbd-eb2f741c8bd0).

- [boom-shield-summary.md](boom-shield-summary.md): narrative summary of the findings from the acoustic analysis thread.
- [hypersonic_corridor_points.csv](hypersonic_corridor_points.csv): reference Mach / altitude / density / dynamic-pressure points for SR-71, X-51A, X-43A, and a speculative SR-72 corridor.
- [shield_opacity_required_points.csv](shield_opacity_required_points.csv): required direct-wave opacity versus speed for X-59, SR-71-like, and SR-72-like cases.
- [shield_spacing_envelope_points.csv](shield_spacing_envelope_points.csv): first toy spacing-envelope results using a simple edge-leakage law.
- [shield_spacing_fresnel_points.csv](shield_spacing_fresnel_points.csv): improved spacing-envelope results using a Fresnel / knife-edge diffraction surrogate.

## Main result

The analysis suggests that a **Width W=1 km-wide boom shield** with aircraft flying on the order of **Elevation E=250 m above the shield** may be acoustically plausible for suppressing the **direct centerline boom wave**, especially for **SR-71-class** overland flight and potentially for **optimistic SR-72-class** scenarios.

The central viability condition is:

`E_min,ops < E_max,acoustic`

where:

- `E_max,acoustic` is the largest aircraft-to-shield spacing allowed by the acoustic shielding model.
- `E_min,ops` is the minimum safe spacing required by clearance, gusts, navigation error, and shield motion.

That means the next design thread should focus on shield material, buoyancy, tethering, structural motion, and control margins.

<img width="640" alt="image" src="https://github.com/user-attachments/assets/ea8d2b0b-5731-4b88-b698-0fdf8acaab59" />
<img width="640" alt="image" src="https://github.com/user-attachments/assets/e04fadad-679b-497b-a106-292406f3ee56" />
<img width="640" alt="image" src="https://github.com/user-attachments/assets/c29d7f8f-a720-4306-b50d-39c54013f077" />
