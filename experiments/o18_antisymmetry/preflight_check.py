"""O18 pre-flight — run and pass BEFORE any sweep number is interpreted.

`docs/25` + `docs/working_rules.md`: dynamic range at the experiment's real N/P, control
potency, empirical null through the whole pipeline, name-content. The spec
(§5) already records that T1's two instruments are pinned across the entire
region of interest, so the primary metric here is deliberately neither of them.
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))

from analysis import preflight as pf                                  # noqa: E402
from experiments.o18_antisymmetry.harness import frac_fixed_only  # noqa: E402

N, P = 500, 50          # the experiment's real parameters (alpha = 0.1)


def frac_fixed(cond: str, lam: float, seed: int) -> float:
    return frac_fixed_only(cond, lam, N=N, P=P, seed=seed)


def main() -> int:
    checks = []

    # -- A: can the primary instrument move, at the REAL N and P? ------------
    # `observed` is the operating point the probe will actually claim at:
    # lambda = 0.30, where the sqrt(alpha) prediction says the structured knee
    # sits. If frac_fixed is pinned THERE, the headline cannot be measured.
    obs = float(np.mean([frac_fixed("structured", 0.30, s) for s in range(3)]))
    checks.append(pf.dynamic_range(
        measure=lambda rng: frac_fixed("structured", 0.30, int(rng.integers(1e6))),
        null=lambda rng: frac_fixed("structured", 0.0, int(rng.integers(1e6))),
        saturated=lambda rng: frac_fixed("structured", 8.0, int(rng.integers(1e6))),
        observed=obs, n_seeds=3,
        name="dynamic_range[frac_fixed @ lambda=0.30, N=500, P=50]"))

    # -- B: controls must be able to produce the outcome --------------------
    # "Plant the effect": a large lambda IS the effect. Each control must move
    # frac_fixed away from the lambda=0 baseline, or it cannot control anything.
    base = float(np.mean([frac_fixed("structured", 0.0, s) for s in range(3)]))
    for cond in ("rank_matched_random", "shuffled_order", "random"):
        checks.append(pf.control_potency(
            measure=lambda c: float(np.mean([frac_fixed(c, 4.0, s) for s in range(3)])),
            apply_control=lambda c: c,
            planted=cond, baseline=base, min_effect=0.1,
            name=f"control_potency[{cond}]"))

    # -- C: the empirical null, through THIS pipeline -----------------------
    # The random-A condition re-run at the same N/P/seeds and the same code
    # path. Structured excesses are compared to this, never to T1's table.
    for lam in (0.20, 0.30):
        checks.append(pf.empirical_null(
            pipeline=lambda rng, lam=lam: frac_fixed("random", lam, int(rng.integers(1e6))),
            n_seeds=5, name=f"empirical_null[random A, lambda={lam}]"))

    # -- D: name-content ----------------------------------------------------
    checks.append(pf.depends_on(
        fn=lambda **kw: frac_fixed("structured", kw["lam"], 0),
        kwargs={"lam": 0.0}, varying="lam", values=[0.0, 0.3, 1.0, 4.0],
        name="depends_on[lambda]"))

    print(pf.report(checks))
    return 0 if all(bool(c) for c in checks) else 1


if __name__ == "__main__":
    raise SystemExit(main())
