"""O18 sweep runner — resumable, seed-pure, one JSONL per (arm, seed).

Spec: ``docs/42_o18_antisymmetry_spec.md``. Usage::

    python3 experiments/o18_antisymmetry/run.py --seeds 0 1 2 3 4
    python3 experiments/o18_antisymmetry/run.py --arm load --seeds 0 1 2 3 4

Each ``(condition, lambda, seed)`` cell is written as it completes and skipped
on re-run, so a sweep split across calls is bit-identical to one run in a
single call (`docs/working_rules.md`: chunk long work into resumable foreground steps).

**Patterns depend only on ``seed``**, so at a given seed all four conditions
see the *same* stored patterns. The design is therefore **paired**, and the
analysis reports per-seed paired differences — which cancels the large
pattern-set variance visible in the pre-flight null (0.7 ± 0.17).
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))

from experiments.o18_antisymmetry.harness import (  # noqa: E402
    CONDITIONS, run_cell, cell_dict,
)

OUT = ROOT / "results" / "o18_antisymmetry"

# Finer at the low end than the spec's first draft: the pre-flight found the
# structured condition already at frac_fixed 0.03 by lambda = 0.30, so a grid
# that starts at 0.05 would miss the knee entirely.
LAMBDAS = [0.0, 0.01, 0.02, 0.03, 0.05, 0.075, 0.10, 0.15, 0.20, 0.30, 0.50, 1.0, 1.5]

ARMS = {
    # main: alpha = 0.1, matching T1's random-A check
    "main": [(500, 50)],
    # P3: does the knee move with load?
    "load": [(500, 25), (500, 50), (500, 100)],
}

COST = dict(max_sweeps=25, n_flip_levels=4, n_basin_probes=3, viol_sweeps=2)


def run(arm: str, seeds: list[int]) -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    for N, P in ARMS[arm]:
        for seed in seeds:
            path = OUT / f"cells_{arm}_N{N}_P{P}_seed{seed}.jsonl"
            done = set()
            if path.exists():
                for line in path.read_text().splitlines():
                    if line.strip():
                        r = json.loads(line)
                        done.add((r["condition"], r["lam"]))
            with path.open("a") as fh:
                for cond in CONDITIONS:
                    for lam in LAMBDAS:
                        if (cond, lam) in done:
                            continue
                        c = run_cell(cond, lam, N=N, P=P, seed=seed, **COST)
                        fh.write(json.dumps(cell_dict(c)) + "\n")
                        fh.flush()
                        print(f"{arm} N{N} P{P} s{seed} {cond:>20} "
                              f"lam={lam:<5} frac_fixed={c.frac_fixed:.3f} "
                              f"noncnv={c.frac_nonconverged:.2f} "
                              f"basin={c.basin_radius:.3f} p_next={c.p_next:+.3f}",
                              flush=True)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--arm", default="main", choices=sorted(ARMS))
    ap.add_argument("--seeds", type=int, nargs="+", default=[0])
    a = ap.parse_args()
    run(a.arm, a.seeds)


if __name__ == "__main__":
    main()
