"""Shared, validated analysis metrics for the Boltzmann Brains probes.

Every metric in this package is validated against synthetic data with a known
ground truth in ``analysis/tests/`` *before* it is trusted on a simulation.
See ``docs/24_m1_implementation_spec.md`` §1.

Design rules:
  - numpy/scipy only; no simulator framework, no sklearn.
  - pure functions where possible; explicit ``rng`` for anything stochastic.
  - every function that can silently lie about its result says so in its
    docstring under **Caveats**.
"""

__all__ = [
    "dimensionality",
    "avalanche",
    "spectral",
    "dynamics",
    "rotation",
    "attractors",
    "decode",
]
