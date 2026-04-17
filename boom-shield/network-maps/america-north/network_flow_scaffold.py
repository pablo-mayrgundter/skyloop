#!/usr/bin/env python3
"""
Scaffold for a first North America boom-shield network solve.

Expected input:
  - north_america_seed_airports.csv
  - north_america_recent_year_od.csv   (to be added later)

Expected OD schema:
  origin,dest,annual_passengers,reference_year

This script:
  1. loads the airport seed set
  2. loads recent-year O&D traffic if present
  3. computes great-circle distance
  4. classifies each pair by likely Mach regime
  5. scores pairs by a simple travel-time-saved surrogate
  6. emits a first-pass prioritized table
"""

from __future__ import annotations

import math
from pathlib import Path
import pandas as pd


ROOT = Path(__file__).resolve().parent
AIRPORTS_CSV = ROOT / "north_america_seed_airports.csv"
OD_CSV = ROOT / "north_america_recent_year_od.csv"
OUT_CSV = ROOT / "north_america_pairs_scored.csv"


def haversine_km(lat1: float, lon1: float, lat2: float, lon2: float) -> float:
    r_km = 6371.0
    p = math.pi / 180.0
    dlat = (lat2 - lat1) * p
    dlon = (lon2 - lon1) * p
    a = (
        math.sin(dlat / 2.0) ** 2
        + math.cos(lat1 * p) * math.cos(lat2 * p) * math.sin(dlon / 2.0) ** 2
    )
    return 2.0 * r_km * math.asin(math.sqrt(a))


def classify_regime(distance_km: float) -> tuple[str, str, float]:
    """
    Returns:
      (regime_class, regime_note, nominal_mach)
    """
    if distance_km < 700:
        return (
            "local",
            "Below threshold for dedicated Mach corridor",
            0.8,
        )
    if distance_km < 1500:
        return (
            "connector",
            "Borderline benefit; Mach 1.2-1.8 regional connector",
            1.5,
        )
    if distance_km < 3000:
        return (
            "corridor",
            "Strong benefit; Mach 1.5-2.5 corridor candidate",
            2.0,
        )
    return (
        "trunk",
        "Very strong benefit; Mach 2.5-3 trunk candidate",
        3.0,
    )


def subsonic_block_hours(distance_km: float) -> float:
    # crude block-speed surrogate for screening, not schedule quality
    return distance_km / 780.0 + 0.75


def mach_corridor_block_hours(distance_km: float, nominal_mach: float) -> float:
    # crude corridor-speed surrogate with ramp / terminal overhead
    speed_kmh = 1060.0 * nominal_mach
    overhead = 0.5 if nominal_mach < 2.0 else 0.7
    return distance_km / speed_kmh + overhead


def main() -> None:
    airports = pd.read_csv(AIRPORTS_CSV)
    od = pd.read_csv(OD_CSV) if OD_CSV.exists() else pd.DataFrame(
        columns=["origin", "dest", "annual_passengers", "reference_year"]
    )

    merged = (
        od.merge(
            airports.add_prefix("o_"),
            left_on="origin",
            right_on="o_iata",
            how="left",
        )
        .merge(
            airports.add_prefix("d_"),
            left_on="dest",
            right_on="d_iata",
            how="left",
        )
        .copy()
    )

    if merged.empty:
        print("No OD traffic file found yet. Add north_america_recent_year_od.csv and rerun.")
        return

    merged["distance_km"] = merged.apply(
        lambda r: haversine_km(
            float(r["o_lat_deg"]),
            float(r["o_lon_deg"]),
            float(r["d_lat_deg"]),
            float(r["d_lon_deg"]),
        ),
        axis=1,
    )

    regimes = merged["distance_km"].apply(classify_regime)
    merged["mach_regime_class"] = regimes.apply(lambda x: x[0])
    merged["mach_regime_note"] = regimes.apply(lambda x: x[1])
    merged["nominal_mach"] = regimes.apply(lambda x: x[2])

    merged["subsonic_block_hours"] = merged["distance_km"].apply(subsonic_block_hours)
    merged["mach_corridor_block_hours"] = merged.apply(
        lambda r: mach_corridor_block_hours(r["distance_km"], r["nominal_mach"]),
        axis=1,
    )
    merged["hours_saved_per_passenger"] = (
        merged["subsonic_block_hours"] - merged["mach_corridor_block_hours"]
    ).clip(lower=0.0)
    merged["annual_passenger_hours_saved"] = (
        merged["annual_passengers"] * merged["hours_saved_per_passenger"]
    )

    cols = [
        "origin",
        "dest",
        "reference_year",
        "annual_passengers",
        "distance_km",
        "mach_regime_class",
        "mach_regime_note",
        "nominal_mach",
        "subsonic_block_hours",
        "mach_corridor_block_hours",
        "hours_saved_per_passenger",
        "annual_passenger_hours_saved",
    ]
    merged = merged[cols].sort_values(
        ["annual_passenger_hours_saved", "annual_passengers", "distance_km"],
        ascending=[False, False, False],
    )

    merged.to_csv(OUT_CSV, index=False)
    print(f"Wrote {OUT_CSV}")


if __name__ == "__main__":
    main()
