# Boom Shield segment analysis — thread notes

Context: these notes capture the tube-to-envelope-to-Boom-Shield pivot from one thread, to be combined later with material from the other Boom Shield thread.

## 1. Structural pivot

Original logic:
- optimize enclosed tube volume
- envelope wraps tubes
- tubes become structural spines
- track loads make everything stiff and heavy

Boom Shield logic:
- optimize planform area per mass
- use buoyant envelope mainly to hold altitude and shape
- use tensioned ring-to-ring suspension
- let the shield deck be a hung acoustic truss / panel system

Primary load path becomes:

`deck load -> radial/catenary tethers -> rings -> outer envelope / global tethers`

This is favorable because it keeps most members in tension.

## 2. Cross-section flattening

Earlier section was roughly:
- width ~120 m
- height ~30 m

Boom Shield candidate family:
- width W = 150 to 300 m
- thickness T = 10 to 20 m

Approximate outer family with ellipse:

`A_env ~= pi W T / 4`

A useful constant-buoyancy-class family keeps `W*T ~= 3600 m^2`:

| Width W | Thickness T | W*T |
|---:|---:|---:|
| 120 m | 30 m | 3600 |
| 180 m | 20 m | 3600 |
| 240 m | 15 m | 3600 |
| 300 m | 12 m | 3600 |

Implication: flattening can preserve rough buoyancy class while greatly improving shield width.

## 3. Module concept

Likely repeated module lengths:
- 100 m
- 250 m
- 500 m

Best first nominal point: `L_m = 250 m`

Each module has:
- ring frame at each end
- catenary support net between rings
- hung shield deck / truss
- outer membrane fairing
- light inter-module joints / overlap lips

The deck should not carry the whole corridor longitudinally. The module is closer to a buoyant, tensioned acoustic raft than a bridge girder.

## 4. Boom interception geometry

Use first-order footprint model:

`sin(mu) = 1/M`

`W_req ~= 2 h tan(mu) = 2 h / sqrt(M^2 - 1)`

Define geometric coverage parameter:

`Lambda = W / W_req = W / (2 h tan(mu))`

Interpretation:
- `Lambda < 1`: lots of side spill
- `Lambda ~ 1`: decent first-pass interception
- `Lambda >> 1`: strong geometric shadow beneath

### Widths required for Lambda = 1 at h = 100 m

| Mach | mu | W_req |
|---:|---:|---:|
| 1.2 | 56.4 deg | 302 m |
| 1.4 | 45.6 deg | 204 m |
| 1.6 | 38.7 deg | 160 m |
| 2.0 | 30.0 deg | 115 m |
| 2.5 | 23.6 deg | 87 m |
| 3.0 | 19.5 deg | 71 m |
| 3.2 | 18.2 deg | 66 m |

Main result: for fixed aircraft-shield gap `h`, higher Mach narrows the Mach cone and reduces required shield width.

### Widths required for selected shield gaps

| Mach | h = 75 m | h = 100 m | h = 150 m | h = 250 m |
|---:|---:|---:|---:|---:|
| 1.2 | 226 | 302 | 452 | 754 |
| 1.4 | 153 | 204 | 306 | 510 |
| 1.6 | 120 | 160 | 240 | 400 |
| 2.0 | 87 | 115 | 173 | 289 |
| 2.5 | 65 | 87 | 131 | 218 |
| 3.0 | 53 | 71 | 106 | 177 |
| 3.2 | 50 | 66 | 99 | 165 |

This suggests a promising regime around:
- Mach 1.6 to 2.5
- aircraft-shield gap h = 75 to 150 m
- shield width W = 200 to 300 m

### Example Lambda values for h = 100 m

| Mach | 100 m shield | 200 m shield | 300 m shield |
|---:|---:|---:|---:|
| 1.2 | 0.33 | 0.66 | 0.99 |
| 1.4 | 0.49 | 0.98 | 1.47 |
| 1.6 | 0.63 | 1.25 | 1.88 |
| 2.0 | 0.87 | 1.73 | 2.60 |
| 2.5 | 1.15 | 2.30 | 3.45 |
| 3.0 | 1.41 | 2.82 | 4.23 |
| 3.2 | 1.52 | 3.03 | 4.55 |

Implication:
- 100 m shield gets serious only above roughly Mach 2
- 200 m shield is useful from roughly Mach 1.4 upward
- 300 m shield catches much of the first-pass footprint even around Mach 1.2 to 1.4

## 5. Areal density framing

Main mass metric:

`sigma = system areal density [kg/m^2]`

Per meter of corridor:

`m' = sigma * W`

### Mass per meter of corridor

For `W = 200 m`:

| sigma | m' |
|---:|---:|
| 2 kg/m^2 | 400 kg/m |
| 4 kg/m^2 | 800 kg/m |
| 6 kg/m^2 | 1200 kg/m |
| 8 kg/m^2 | 1600 kg/m |

For `W = 300 m`:

| sigma | m' |
|---:|---:|
| 2 kg/m^2 | 600 kg/m |
| 4 kg/m^2 | 1200 kg/m |
| 6 kg/m^2 | 1800 kg/m |
| 8 kg/m^2 | 2400 kg/m |

Working rubric:
- green: 2 to 4 kg/m^2
- yellow: 4 to 8 kg/m^2
- red: > 8 kg/m^2

## 6. Tethering recap

There are three tether systems:

### A. Internal suspension tethers
Carry:
- dead load of deck
- local acoustic loads
- local gust and vibration loads

Use suspension-bridge approximation for span `s`, sag `f`, distributed line load `q`:

`H ~= q s^2 / (8 f)`

With areal density `sigma`, per meter of corridor use:

`q ~= sigma g`

Example for `sigma = 4 kg/m^2`, `s = 125 m`, `f = 15 m`:

`q ~= 39 N/m^2`

`H ~= 39 * 125^2 / (8 * 15) ~= 5.1 kN/m`

This suggests dead-load support tensions are manageable.

### B. Shape-control tethers
Carry:
- aerodynamic shape loads
- membrane tension redistribution
- transient off-axis disturbance

These likely require redundancy, segmentation, and possibly active trimming.

### C. Ground / station tethers
Carry:
- corridor position-keeping
- net crosswind loads
- thermal drift / alignment correction
- fault capture

Important architectural rule:
- do not try to transmit all lateral loads longitudinally through the entire corridor
- each module or small cluster should be laterally stabilized mostly locally

Preferred load path:

`module loads -> local station tethers -> ground`

not:

`module loads -> whole corridor spine`

## 7. Crosswind load note

Crude side-force per unit length:

`F'_perp ~= 0.5 rho U_perp^2 C_D T`

For a constant-buoyancy family, flattening reduces thickness `T` and therefore reduces crosswind side projected area. So flattening helps both:
- boom interception geometry
- crosswind side force

## 8. Structural interpretation of the catenary ring-to-ring concept

The catenary-ring concept says:
- rings take hoop load and preserve geometry
- catenary tethers carry longitudinal support in tension
- central truss/deck hangs between ring frames
- radial symmetry gives a clean way to survive off-axis disturbance without heavy compression members

For Boom Shield, this evolves toward:
- strong end rings
- upper cable net / catenary fans
- hung acoustic deck below, not perfectly flat
- slight camber or shallow V may help both structure and side-leakage redirection
- outer membrane serves buoyancy and aerodynamic shaping, not heavy primary deck load

A useful picture is a lenticular suspension bridge in section.

## 9. Nominal candidate configurations to compare

### Flight cases

| Case | Mach | gap h |
|---|---:|---:|
| A | 1.4 | 100 m |
| B | 1.6 | 100 m |
| C | 2.0 | 100 m |
| D | 2.5 | 100 m |
| E | 3.0 | 100 m |
| F | 2.0 | 150 m |

### Shield cases

| Case | Width W | Thickness T | Module L_m | sigma |
|---|---:|---:|---:|---:|
| S1 | 200 m | 18 m | 250 m | 4 kg/m^2 |
| S2 | 250 m | 15 m | 250 m | 4 kg/m^2 |
| S3 | 300 m | 12 m | 250 m | 4 kg/m^2 |
| S4 | 250 m | 15 m | 250 m | 8 kg/m^2 |

## 10. Early conclusions from this thread

- Flattening is the right move.
- Higher Mach helps shield width strongly by simple footprint geometry.
- 200 to 300 m width is the right first shield band to study.
- Internal suspension looks manageable in static dead-load terms.
- Ground/station tethering is likely the real positioning challenge.
- The sweet regime for Boom Shield may sit above X-59 territory, not below it, because shield width becomes easier with Mach even though aircraft-side penalties rise.
- The concept likely lives or dies on areal density and acoustic effectiveness, not on obvious static structural impossibility.

## 11. Recommended next synthesis step

Combine this file with notes from the parallel Boom Shield thread, then produce:
1. one merged `Lambda(M, h, W)` design-space matrix
2. one areal-mass budget by subsystem
3. one tethering/load-path summary sketch
4. one first cost rubric based on width, areal density, module length, and ground anchor spacing
