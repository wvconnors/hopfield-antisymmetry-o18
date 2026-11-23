"""O18 adversarial-review checks — the four corrections, re-runnable.

Every claim in ``notes/o18_antisymmetry.md`` §7 that overturned a V1 conclusion
is reproduced here, so the corrections are artifacts rather than assertions.
V1's own numbers are not re-derived (``analyze.py`` already reproduces them
bit-for-bit against the committed ``analysis.txt``); what is derived here is
what V1 could not see.

    PYTHONPATH=. python3 experiments/o18_antisymmetry/review_checks.py

Runtime ~2 min. Checks, in the order the write-up uses them:

1. **units** — field delivered to the stored patterns per unit ``‖A‖_F``, and
   ``knee × field`` across conditions. The headline correction: the 3× knee
   separation is ``√(N/P)``, and ``knee × field`` is flat.
2. **relabel** — ``shuffled_order`` is an exact relabelling, so its null was an
   algebraic identity (Type B).
3. **sign** — the dynamics applies ``−A`` (``h = S @ W[:, i]`` is ``Wᵀs``), so
   ``p_next`` read the lobe the field pushes away from; ``p_prev`` mirrored.
   Plus the reachability arithmetic that made P2 dead on arrival.
4. **law** — the load law decomposed into unit conversion × residual, and the
   per-seed spread the leave-one-out range understated.
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from experiments.p1_unlearning.hopfield import random_patterns  # noqa: E402
from experiments.o18_antisymmetry import harness as H  # noqa: E402

RESULTS = Path(__file__).resolve().parents[2] / "results" / "o18_antisymmetry"

# λ grid fine enough that interpolation bias is ~0.001 (V1's coarse grid
# inflated the *random* knee — the denominator of every ratio — by ~5 %).
LAMS = np.concatenate([np.arange(0, 0.4, 0.0025), np.arange(0.4, 6.01, 0.02)])


def knee(lams, ff, frac: float = 0.5) -> float:
    """λ at which ``frac_fixed`` reaches ``frac`` × its λ=0 value, interpolated."""
    ff = np.asarray(ff, dtype=float)
    target = frac * ff[0]
    for i in range(1, len(ff)):
        if ff[i] <= target:
            x0, x1, y0, y1 = lams[i - 1], lams[i], ff[i - 1], ff[i]
            return float(x0 + (y0 - target) * (x1 - x0) / (y0 - y1))
    return float("nan")


def field_per_fro(A: np.ndarray, patterns: np.ndarray) -> float:
    """``RMS_{i,μ}(A ξ^μ)_i / ‖A‖_F`` — field delivered per unit Frobenius norm.

    This is the quantity ``sign(Ws)`` responds to, and the one ``λ`` was
    assumed to be a proxy for. It is not: it varies with the operator's support.
    """
    return float(np.sqrt(((A @ patterns.T) ** 2).mean())) / float(np.linalg.norm(A))


def operator(cond: str, patterns: np.ndarray, seed: int) -> np.ndarray:
    """The same operator ``build``/``frac_fixed_curve`` construct, unnormalised."""
    rng = np.random.default_rng(seed + 10_000)
    if cond == "structured":
        return H.stdp_sequence_A(patterns)
    if cond == "shuffled_order":
        return H.stdp_sequence_A(patterns, order=rng.permutation(patterns.shape[0]))
    if cond == "random":
        return H.random_A(patterns.shape[1], rng)
    if cond == "rank_matched_random":
        return H.rank_matched_random_A(patterns, rng)
    raise ValueError(cond)


def fit_loglog(x, y) -> tuple[float, float]:
    b, c = np.polyfit(np.log(x), np.log(y), 1)
    return float(b), float(np.exp(c))


# --------------------------------------------------------------------------
# 1. the units: knee × field is flat, and the 3× is √(N/P)
# --------------------------------------------------------------------------

def check_units(N: int = 500, P: int = 50, seeds: int = 5) -> None:
    print("\n=== 1. units — the headline correction ===")
    print(f"{'condition':22s}{'rank':>6}{'field/||A||F':>14}{'knee':>10}{'sd':>8}{'knee*field':>12}")
    for cond in H.CONDITIONS[:1] + tuple(c for c in H.CONDITIONS if c != "structured"):
        ks, fs, rk = [], [], []
        for seed in range(seeds):
            ks.append(knee(LAMS, H.frac_fixed_curve(cond, LAMS, N, P, seed)))
            pats = random_patterns(P, N, np.random.default_rng(seed))
            A = operator(cond, pats, seed)
            fs.append(field_per_fro(A, pats))
            if seed == 0:
                rk.append(int(np.linalg.matrix_rank(A)))
        m, sd, f = float(np.mean(ks)), float(np.std(ks)), float(np.mean(fs))
        print(f"{cond:22s}{rk[0]:6d}{f:14.4f}{m:10.4f}{sd:8.4f}{m * f:12.4f}")
    print(f"  sqrt(N/P) = {np.sqrt(N / P):.3f}   "
          "<- the whole knee separation; knee*field is flat to ~3%")


# --------------------------------------------------------------------------
# 2. shuffled_order: the null was an algebraic identity (Type B)
# --------------------------------------------------------------------------

def check_relabel(N: int = 500, P: int = 50, seed: int = 0) -> None:
    print("\n=== 2. shuffled_order is an exact relabelling ===")
    pats = random_patterns(P, N, np.random.default_rng(seed))
    sig = np.random.default_rng(seed + 7).permutation(P)
    d_A = float(np.abs(H.stdp_sequence_A(pats, order=sig) - H.stdp_sequence_A(pats[sig])).max())
    S1 = (pats.T @ pats) / N
    S2 = (pats[sig].T @ pats[sig]) / N
    print(f"  max|A(pat, order=sigma) - A(pat[sigma])| = {d_A:.3e}")
    print(f"  max|S(pat) - S(pat[sigma])|              = {float(np.abs(S1 - S2).max()):.3e}")
    print("  -> with i.i.d. patterns the two conditions have the SAME joint law;")
    print("     every permutation-invariant statistic matches by construction.")


# --------------------------------------------------------------------------
# 3. the sign convention, and P2's reachability
# --------------------------------------------------------------------------

def check_sign(N: int = 500, P: int = 50, seeds: int = 8) -> None:
    print("\n=== 3. sign convention — p_next read the wrong lobe ===")
    print(f"{'lam':>6}{'p_next (recorded)':>22}{'p_prev (the drive)':>22}{'p_self':>10}{'frac_fixed':>12}")
    for lam in (0.10, 0.15, 0.20, 0.30):
        nx, pv, sf, ff = [], [], [], []
        for seed in range(seeds):
            pats = random_patterns(P, N, np.random.default_rng(seed))
            net, _ = H.build(pats, "structured", lam,
                             np.random.default_rng(seed + 10_000), compute_rank=False)
            final, _, _ = net.converge(pats.copy(), np.random.default_rng(seed + 999),
                                       max_sweeps=25)
            ov = (final @ pats.T) / N
            nx.append(np.diag(np.roll(ov, -1, axis=1)).mean())
            pv.append(np.diag(np.roll(ov, 1, axis=1)).mean())
            sf.append(np.diag(ov).mean())
            ff.append(net.is_fixed_point(pats).mean())
        f = lambda v: (float(np.mean(v)), float(np.std(v, ddof=1)))  # noqa: E731
        (a, asd), (b, bsd) = f(nx), f(pv)
        print(f"{lam:6.2f}   {a:+.4f} +/- {asd:.4f}      {b:+.4f} +/- {bsd:.4f}"
              f"{np.mean(sf):10.4f}{np.mean(ff):12.4f}")

    pats = random_patterns(P, N, np.random.default_rng(0))
    S = (pats.T @ pats) / N
    np.fill_diagonal(S, 0.0)
    A = H.stdp_sequence_A(pats)
    A = A * (np.linalg.norm(S) / np.linalg.norm(A))
    self_t = float(np.mean([(S @ pats[m]) @ pats[m] / N for m in range(P)]))
    chain = float(np.mean([(A @ pats[m]) @ pats[(m + 1) % P] / N for m in range(P)]))
    print(f"\n  self term {self_t:.4f} ; chain term {chain:.4f} "
          f"-> lambda needed = {self_t / abs(chain):.2f}")
    print("  fixed points are gone by lambda ~ 0.3, so P2 was unreachable at every")
    print("  lambda where the probe had a landscape to measure.")


# --------------------------------------------------------------------------
# 4. the load law, decomposed
# --------------------------------------------------------------------------

def check_law(N: int = 500, seeds: int = 5) -> None:
    print("\n=== 4. the load law decomposes into unit conversion x residual ===")
    rows = [r for r in json.load(open(RESULTS / "knee_law.json")) if r["base"] > 0.5]
    a, fro, fld = [], [], []
    print(f"{'alpha':>7}{'kneeStr':>9}{'kneeRnd':>9}{'froRatio':>10}"
          f"{'1/fieldRatio':>14}{'fieldMatched':>14}")
    for r in rows:
        P = r["P"]
        fs, fr = [], []
        for seed in range(seeds):
            pats = random_patterns(P, N, np.random.default_rng(seed))
            fs.append(field_per_fro(operator("structured", pats, seed), pats))
            fr.append(field_per_fro(operator("random", pats, seed), pats))
        fs, fr = float(np.mean(fs)), float(np.mean(fr))
        ks = float(np.mean(r["knee_structured"]))
        kr = float(np.mean(r["knee_random"]))
        a.append(r["alpha"]); fro.append(ks / kr); fld.append((ks * fs) / (kr * fr))
        print(f"{r['alpha']:7.3f}{ks:9.3f}{kr:9.3f}{ks / kr:10.4f}"
              f"{fr / fs:14.4f}{(ks * fs) / (kr * fr):14.4f}")
    a = np.array(a)
    bf, cf = fit_loglog(a, np.array(fro))
    bm, cm = fit_loglog(a, np.array(fld))
    conv_b = bf - bm
    print(f"\n  Frobenius knee ratio : {cf:.4f} * alpha^{bf:+.4f}   (V1's 'rejection' of sqrt-alpha)")
    print(f"  unit conversion      :          alpha^{conv_b:+.4f}   <- this IS sqrt(alpha), as pre-registered")
    print(f"  field-matched residual: {cm:.4f} * alpha^{bm:+.4f}   <- and it runs the OTHER way")

    per_seed = []
    for s in range(seeds):
        y = np.array([r["knee_structured"][s] / r["knee_random"][s] for r in rows])
        per_seed.append(fit_loglog(a, y)[0])
    per_seed = np.array(per_seed)
    print(f"\n  per-seed exponents {np.round(per_seed, 3).tolist()}")
    print(f"  -> {per_seed.mean():.3f} +/- {per_seed.std(ddof=1):.3f}  (V1 quoted a leave-one-out")
    print("     range over LOADS, +/-0.007, which understates this ~10x; 0.5 still excluded ~4 sigma)")


if __name__ == "__main__":
    check_units()
    check_relabel()
    check_sign()
    check_law()
