# North America network maps

This folder starts a first-pass **North America boom-shield corridor solve** for the Skyloop / boom-shield concept.

The framing from the thread is:

- the shield is not just an acoustic surface, but **flight infrastructure**
- aircraft fly a designed path a small distance `E` above the shield
- the corridor is shaped by **q-limited climb / acceleration splines**
- the network should look more like **freeways with ramps and trunks** than point-to-point free flight
- high-speed trunks can be shared across many origin-destination pairs, with slower merge / egress ramps at hubs

## What is in this folder

- `README.md` — this summary and methodology
- `north_america_seed_airports.csv` — seed airport set for the first solve
- `north_america_candidate_pairs.csv` — all pairwise distances for the seed set, with a first-pass Mach-regime classification by stage length
- `network_flow_scaffold.py` — scaffold for turning a recent-year traffic table into a shared-corridor graph solve

## Status

This commit sets up the **geometry and solve scaffold**.

The traffic-weighted solve still needs a recent reference-year commercial passenger dataset. The intended next input is a table like:

- `origin`
- `dest`
- `annual_passengers`
- `reference_year`

A good source would be a recent BTS / DB1B / T-100 or equivalent North America commercial O&D dataset. Once dropped in, the scaffold script can:

1. join O&D demand to the airport seed set
2. compute distance-based Mach-regime benefit
3. threshold out short-haul pairs with weak supersonic benefit
4. identify shared trunk corridors by aggregated flow
5. output candidate shield mainlines and regional ramps

## First-pass design logic from the thread

The working network intuition is:

- **local / short-haul pairs**: likely not worth dedicated Mach-shield infrastructure
- **regional connectors**: plausible in roughly the Mach 1.2–1.8 range
- **mid-continent corridors**: strong candidates for Mach 1.5–2.5
- **long-haul trunks**: strongest candidates for Mach 2.5–3, especially where multiple city-pairs share a mainline

This is consistent with the shield concept because the aircraft can remain only a few hundred meters above the shield while the corridor itself climbs and accelerates along a q-managed spline.

## Pair classification rule used here

This initial pass uses only great-circle stage length:

- `< 700 km`: below threshold for dedicated Mach corridor
- `700–1500 km`: borderline / regional connector
- `1500–3000 km`: strong corridor candidate
- `> 3000 km`: very strong trunk candidate

That is only a placeholder rule. It should later be replaced by a utility model that combines:

- annual passengers
n- total passenger-hours saved
- route overlap / shared-trunk benefit
- capex per shield-km
- merge / egress complexity
- plausible cruise Mach by corridor class

## Next suggested files

Once recent-year O&D data is available, add:

- `north_america_recent_year_od.csv`
- `north_america_top_pairs_by_volume.csv`
- `north_america_shared_trunks.csv`
- `north_america_corridor_graph.geojson`
- simple maps for western, central, and eastern trunk candidates
