"""Gate-0 rise-time proxy and scored loudness metrics (loudness-agnostic helper).

The propagation model gives ground rise times (~0.7-0.9 ms for Concorde) shorter than
measured booms, whose rise times are lengthened by turbulence (typically 2-5 ms).  Gate 0
therefore also scored every signature after smoothing with a Gaussian kernel of standard
deviation ``RISE_SIGMA_S`` = 1 ms (an ideal step then has a 10-90 % rise of 2.56 ms),
and took as headline the less favourable of the raw and proxy *differences* versus the
unshielded reference.  This module holds that proxy so the Gate-0 driver and the 3-D
carpet code use one definition (numerics identical to the Gate-0 study).

Public API
----------
rise_time_proxy(t, p, sigma_s=RISE_SIGMA_S) -> p_smooth
scored_metrics(t, p, taper_s=METRIC_TAPER_S, sigma_s=RISE_SIGMA_S) -> dict
    loudness.metrics of the raw signature plus ``pldb_proxy``, ``peak_pa_proxy``,
    ``rise_ms_proxy`` (10-90 %) and ``csel_proxy`` of the smoothed one.
headline_delta(shielded, reference) -> max(dPLdB raw, dPLdB proxy)
"""
from __future__ import annotations

import numpy as np
from scipy.ndimage import gaussian_filter1d

from .loudness import metrics

RISE_SIGMA_S = 1.0e-3     # Gaussian kernel std -> ideal step 10-90 % = 2.56 ms
METRIC_TAPER_S = 0.05     # half-Hann taper at both record ends for the loudness metrics


def rise_time_proxy(t, p, sigma_s: float = RISE_SIGMA_S):
    """Gaussian-smoothed signature (std ``sigma_s`` seconds, zero outside the record)."""
    fs = 1.0 / float(t[1] - t[0])
    return gaussian_filter1d(np.asarray(p, dtype=float), sigma_s * fs, mode="constant")


def scored_metrics(t, p, taper_s: float = METRIC_TAPER_S, sigma_s: float = RISE_SIGMA_S):
    """Raw loudness.metrics plus the rise-time-proxy PLdB (Gate-0 conventions)."""
    m = dict(metrics(t, p, taper_s=taper_s))
    ms = metrics(t, rise_time_proxy(t, p, sigma_s), taper_s=taper_s)
    m["pldb_proxy"] = ms["pldb"]
    m["peak_pa_proxy"] = ms["peak_pa"]
    m["rise_ms_proxy"] = ms["rise_time_10_90_ms"]
    m["csel_proxy"] = ms["csel"]
    return m


def headline_delta(shielded: dict, reference: dict) -> float:
    """Gate-0 headline dPLdB: the less favourable (larger) of the raw and proxy deltas."""
    return max(shielded["pldb"] - reference["pldb"],
               shielded["pldb_proxy"] - reference["pldb_proxy"])
