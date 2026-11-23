"""P3 / the sqrt(alpha) law across loads — the arm that makes `alpha` non-constant.

`scan_table` failed the main arm for `N`, `P`, `alpha` constant. That is a real
flag, not a nuisance: the main arm cannot see whether the sqrt(alpha) ratio P1
predicts is a *law* or a coincidence at one load. This arm varies the load and
tests the prediction across it.

It uses `frac_fixed_only` (a single matmul, no dynamics), so a grid an order of
magnitude finer than the main arm costs seconds. The primary metric is
identical -- `test_frac_fixed_only_agrees_with_run_cell` pins them together.
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))
from experiments.o18_antisymmetry.harness import frac_fixed_curve, CONDITIONS  # noqa: E402

OUT = ROOT / "results" / "o18_antisymmetry"
LAMS = np.round(np.concatenate([np.arange(0.0, 0.4, 0.01),
                                np.arange(0.4, 1.61, 0.05)]), 4).tolist()
LOADS = [(500, 15), (500, 25), (500, 50), (500, 75), (500, 100)]   # alpha .03 -> .20
SEEDS = list(range(5))


def knee(lams, vals, frac=0.5):
    base = vals[0]
    thr = frac * base
    for i in range(1, len(vals)):
        if vals[i] < thr:
            x0, x1, y0, y1 = lams[i - 1], lams[i], vals[i - 1], vals[i]
            return x0 + (y0 - thr) * (x1 - x0) / (y0 - y1) if y0 != y1 else x1
    return float("nan")


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    rows = []
    path = OUT / "load_sweep.jsonl"
    for N, P in LOADS:
        for seed in SEEDS:
            for cond in CONDITIONS:
                vals = frac_fixed_curve(cond, LAMS, N, P, seed)
                rows.append({"N": N, "P": P, "alpha": P / N, "seed": seed,
                             "condition": cond, "knee": knee(LAMS, vals),
                             "base": vals[0]})
    path.write_text("\n".join(json.dumps(r) for r in rows) + "\n")

    print(f"{'alpha':>7} {'sqrt(a)':>8} | "
          + " ".join(f"{c[:11]:>14}" for c in CONDITIONS) + " | ratio str/rnd")
    for N, P in LOADS:
        a = P / N
        sub = [r for r in rows if r["P"] == P]
        cells, means = [], {}
        for c in CONDITIONS:
            k = np.array([r["knee"] for r in sub if r["condition"] == c], float)
            means[c] = np.nanmean(k)
            cells.append(f"{np.nanmean(k):>7.3f}±{np.nanstd(k, ddof=1):<6.3f}")
        ratio = means["structured"] / means["random"]
        print(f"{a:>7.3f} {np.sqrt(a):>8.3f} | " + " ".join(cells)
              + f" | {ratio:>7.3f}")


if __name__ == "__main__":
    main()
