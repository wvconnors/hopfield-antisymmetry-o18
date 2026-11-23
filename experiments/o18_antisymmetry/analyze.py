"""O18 analysis. Runs `scan_table` BEFORE printing anything interpretable.

Design is **paired**: at a given seed every condition sees the same stored
patterns, so differences are reported per-seed and then averaged, which
removes the large pattern-set variance (pre-flight null 0.7 ± 0.17).
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))
from analysis import preflight as pf   # noqa: E402

OUT = ROOT / "results" / "o18_antisymmetry"
CONDS = ["structured", "shuffled_order", "rank_matched_random", "random"]


def load(arm="main", N=500, P=50):
    rows = []
    for p in sorted(OUT.glob(f"cells_{arm}_N{N}_P{P}_seed*.jsonl")):
        rows += [json.loads(l) for l in p.read_text().splitlines() if l.strip()]
    return rows


def knee(lams, vals, frac=0.5):
    """First lambda at which the metric falls below `frac` of its lambda=0 value.

    Linearly interpolated. Reported per seed so it carries a spread.
    """
    base = vals[0]
    thr = frac * base
    for i in range(1, len(vals)):
        if vals[i] < thr:
            x0, x1, y0, y1 = lams[i - 1], lams[i], vals[i - 1], vals[i]
            return x0 + (y0 - thr) * (x1 - x0) / (y0 - y1) if y0 != y1 else x1
    return float("nan")


def main():
    rows = load()
    lams = sorted({r["lam"] for r in rows})
    seeds = sorted({r["seed"] for r in rows})

    # ---- read our own output adversarially, before any prose --------------
    print(pf.report([pf.scan_table(rows)]))
    print()

    def grid(metric):
        g = {}
        for c in CONDS:
            g[c] = np.array([[next(r[metric] for r in rows
                                   if r["condition"] == c and r["lam"] == l and r["seed"] == s)
                              for l in lams] for s in seeds])
        return g

    ff = grid("frac_fixed")

    print(f"=== frac_fixed (PRIMARY), mean ± SD over {len(seeds)} seeds, N=500 P=50 α=0.1 ===")
    print(f"{'lam':>6} " + " ".join(f"{c[:18]:>19}" for c in CONDS))
    for i, l in enumerate(lams):
        cells = [f"{ff[c][:, i].mean():>9.3f} ±{ff[c][:, i].std(ddof=1):<8.3f}" for c in CONDS]
        print(f"{l:>6} " + " ".join(cells))

    print(f"\n=== PAIRED difference vs the `random` empirical null (per seed, then averaged) ===")
    print("negative = breaks EARLIER than isotropic random antisymmetry")
    print(f"{'lam':>6} " + " ".join(f"{c[:18]:>19}" for c in CONDS[:3]))
    for i, l in enumerate(lams):
        cells = []
        for c in CONDS[:3]:
            d = ff[c][:, i] - ff["random"][:, i]
            m, s = d.mean(), d.std(ddof=1)
            mark = "*" if abs(m) > s else " "
            cells.append(f"{m:>+9.3f} ±{s:<7.3f}{mark}")
        print(f"{l:>6} " + " ".join(cells))
    print("  * = |difference| exceeds its across-seed SD (the resolvability bar)")

    print(f"\n=== KNEE: lambda at which frac_fixed falls to half its lambda=0 value ===")
    print(f"{'condition':>20}  {'mean':>7} {'SD':>7}   per-seed")
    knees = {}
    for c in CONDS:
        k = np.array([knee(lams, ff[c][j]) for j in range(len(seeds))])
        knees[c] = k
        print(f"{c:>20}  {np.nanmean(k):>7.3f} {np.nanstd(k, ddof=1):>7.3f}   "
              + " ".join(f"{x:.3f}" for x in k))

    alpha = 0.1
    ratio = np.nanmean(knees["structured"]) / np.nanmean(knees["random"])
    print(f"\nP1 (pre-registered): structured knee / random knee ≈ sqrt(alpha) = {np.sqrt(alpha):.3f}")
    print(f"P1 observed        : {ratio:.3f}")

    print(f"\n=== CONTROL: does the sequence structure matter at all? ===")
    print("structured vs its rank/subspace-matched controls, PAIRED, at each lambda")
    for i, l in enumerate(lams):
        out = []
        for c in ("shuffled_order", "rank_matched_random"):
            d = ff["structured"][:, i] - ff[c][:, i]
            m, s = d.mean(), d.std(ddof=1)
            out.append(f"{c[:11]:>11} {m:>+7.3f}±{s:<6.3f}{'*' if abs(m) > s else ' '}")
        print(f"  lam={l:<6} " + "  ".join(out))

    print(f"\n=== P2 (mechanism): p_next, overlap with xi^(mu+1) after settling ===")
    pn = grid("p_next")
    for i, l in enumerate(lams):
        print(f"  lam={l:<6} " + "  ".join(
            f"{c[:11]:>11} {pn[c][:, i].mean():>+6.3f}±{pn[c][:, i].std(ddof=1):<5.3f}"
            for c in ("structured", "random")))

    print(f"\n=== secondary: energy_violations (spec §5 records this is floored) ===")
    ev = grid("energy_violations")
    for i, l in enumerate(lams):
        print(f"  lam={l:<6} " + "  ".join(
            f"{c[:11]:>11} {ev[c][:, i].mean():>8.2f}" for c in CONDS))


if __name__ == "__main__":
    main()
