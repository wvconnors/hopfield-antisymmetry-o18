"""O18 — structured vs random antisymmetry in the Hopfield energy landscape.

Spec: ``docs/42_o18_antisymmetry_spec.md``.

The validated ``experiments/p1_unlearning/hopfield.py`` is used **unmodified**.
This module only builds antisymmetric operators and installs them on the net's
``W``; every dynamical routine (async ``converge``, ``is_fixed_point``,
``energy``) is probe 1's.

One convenience is load-bearing enough to be tested rather than assumed
(``test_harness.py::test_energy_ignores_antisymmetric_part``): for any
antisymmetric ``A``, ``sᵀAs = 0`` identically, so ``HopfieldNet.energy``
returns the **symmetric-part** energy even after ``A`` is installed. That is
exactly what "energy-increasing update" must mean here.
"""

from __future__ import annotations

import sys
from dataclasses import dataclass, asdict
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from experiments.p1_unlearning.hopfield import HopfieldNet, random_patterns  # noqa: E402

__all__ = [
    "antisym", "stdp_sequence_A", "random_A", "rank_matched_random_A",
    "build", "Cell", "run_cell", "CONDITIONS",
]

CONDITIONS = ("structured", "random", "rank_matched_random", "shuffled_order")


# --------------------------------------------------------------------------
# operators
# --------------------------------------------------------------------------

def antisym(M: np.ndarray) -> np.ndarray:
    """Antisymmetric part ``(M − Mᵀ)/2``, diagonal zeroed (it is zero anyway)."""
    A = 0.5 * (M - M.T)
    np.fill_diagonal(A, 0.0)
    return A


def stdp_sequence_A(patterns: np.ndarray, order: np.ndarray | None = None) -> np.ndarray:
    """Antisymmetric coupling from an STDP kernel applied to a stored *sequence*.

    The Sompolinsky–Kanter construction (PRL 57:2861): consecutive patterns in
    the sequence are chained ``ξ^{μ+1} ξ^{μᵀ}``, then antisymmetrised. This is
    the discrete-time reduction of a temporally-antisymmetric learning window
    (C1, #223): pre-before-post potentiates the forward pair, post-before-pre
    depresses the reverse one, which is exactly ``M − Mᵀ``.

    ``order`` permutes the temporal order without touching the patterns — the
    shuffled-order control (spec §6.2).
    """
    P, N = patterns.shape
    idx = np.arange(P) if order is None else np.asarray(order)
    xi = patterns[idx]
    # forward chain, wrapped: a cycle, so no pattern is privileged as an endpoint
    M = (np.roll(xi, -1, axis=0).T @ xi) / N
    return antisym(M)


def random_A(N: int, rng: np.random.Generator) -> np.ndarray:
    """Isotropic random antisymmetric operator — the null case (T1's condition)."""
    return antisym(rng.standard_normal((N, N)))


def rank_matched_random_A(patterns: np.ndarray, rng: np.random.Generator) -> np.ndarray:
    """THE control (spec §6.1): random, but confined to the pattern subspace.

    Same rank and same row space as :func:`stdp_sequence_A`, with the sequence
    structure destroyed — built by chaining the patterns to *random Gaussian
    combinations of the patterns* rather than to their successors. If this
    breaks descent where the structured operator does, the effect is about low
    rank and subspace alignment, not about sequence structure or STDP.
    """
    P, N = patterns.shape
    G = rng.standard_normal((P, P))
    G /= np.linalg.norm(G, axis=1, keepdims=True) * np.sqrt(P)
    M = ((G @ patterns).T @ patterns) / N
    return antisym(M)


def build(patterns: np.ndarray, condition: str, lam: float,
          rng: np.random.Generator, compute_rank: bool = True) -> tuple[HopfieldNet, dict]:
    """Hopfield net with ``W = S + λ·Â``, ``‖Â‖_F = ‖S‖_F`` (asserted).

    ``λ`` is therefore exactly ``‖A‖/‖S‖`` — the swept ratio. Norm matching *is*
    the measurement (spec §6.3), so it is asserted per cell rather than trusted.
    """
    net = HopfieldNet(patterns)
    net.assert_valid()                     # S is symmetric & hollow before we touch it
    S = net.W.copy()
    s_norm = float(np.linalg.norm(S))

    if condition == "structured":
        A = stdp_sequence_A(patterns)
    elif condition == "shuffled_order":
        A = stdp_sequence_A(patterns, order=rng.permutation(patterns.shape[0]))
    elif condition == "random":
        A = random_A(patterns.shape[1], rng)
    elif condition == "rank_matched_random":
        A = rank_matched_random_A(patterns, rng)
    else:
        raise ValueError(f"unknown condition {condition!r}")

    a_norm = float(np.linalg.norm(A))
    assert a_norm > 0, f"{condition}: degenerate A"
    A = A * (s_norm / a_norm)              # matched Frobenius norm
    assert np.max(np.abs(A + A.T)) < 1e-9, "A lost antisymmetry"

    net.W = S + lam * A
    ratio = float(np.linalg.norm(lam * A) / s_norm)
    assert abs(ratio - lam) < 1e-9, f"norm ratio {ratio} != lambda {lam}"
    # rank is an SVD and dominates the cost of a lambda sweep; it is diagnostic
    # only, so the sweep path turns it off.
    rank = int(np.linalg.matrix_rank(A)) if compute_rank else -1
    return net, {"rank_A": rank, "ratio": ratio}


# --------------------------------------------------------------------------
# metrics
# --------------------------------------------------------------------------

def _energy_violations(net: HopfieldNet, Sym: np.ndarray, S0: np.ndarray,
                       rng: np.random.Generator, n_sweeps: int) -> float:
    """Mean count of async updates that *increase* the symmetric-part energy.

    Secondary metric only — the spec (§5) records that it is floored at 0
    across the whole region of interest, which is why it is not primary.

    ``E = −½ sᵀ Sym s`` with ``Sym`` symmetric and hollow, so flipping neuron
    ``i`` by ``δ = s_new − s_old`` gives **exactly** ``ΔE = −δ·(Sym s)_i``
    (the ``S_ii`` term vanishes). That is ``O(N)`` per update instead of the
    ``O(N²)`` of re-evaluating the quadratic form, and
    ``test_energy_violations_matches_bruteforce`` pins the two together — an
    optimisation that silently changed the quantity would be a name-content
    failure of exactly the kind `docs/25` catalogues.
    """
    S = np.array(np.atleast_2d(S0), dtype=float)
    viol = np.zeros(S.shape[0])
    for _ in range(n_sweeps):
        for i in rng.permutation(net.N):
            h = S @ net.W[:, i]
            new = np.where(h >= 0, 1.0, -1.0)
            delta = new - S[:, i]
            dE = -delta * (S @ Sym[:, i])
            viol += (dE > 1e-9)
            S[:, i] = new
    return float(viol.mean())


@dataclass
class Cell:
    condition: str
    lam: float
    N: int
    P: int
    seed: int
    alpha: float
    frac_fixed: float          # PRIMARY
    frac_nonconverged: float
    basin_radius: float
    mean_overlap: float
    p_next: float
    energy_violations: float
    rank_A: int


def frac_fixed_curve(condition: str, lams, N: int, P: int, seed: int) -> list[float]:
    """``frac_fixed`` over a whole λ grid, building ``S`` and ``A`` **once**.

    ``A`` does not depend on λ, so rebuilding it per point (and re-deriving its
    rank, an SVD) made a fine grid cost minutes instead of seconds. The λ=0
    entry is checked against a fresh :func:`build` so the fast path cannot
    drift from the slow one silently.
    """
    rng = np.random.default_rng(seed)
    patterns = random_patterns(P, N, rng)
    base, _ = build(patterns, condition, 0.0,
                    np.random.default_rng(seed + 10_000), compute_rank=False)
    S = base.W.copy()
    unit, _ = build(patterns, condition, 1.0,
                    np.random.default_rng(seed + 10_000), compute_rank=False)
    A = unit.W - S                       # already norm-matched to ‖S‖
    out = []
    for lam in lams:
        base.W = S + lam * A
        out.append(float(base.is_fixed_point(patterns).mean()))
    return out


def frac_fixed_only(condition: str, lam: float, N: int, P: int, seed: int) -> float:
    """Primary metric alone, with no dynamics — the cheap path for pre-flight.

    ``is_fixed_point`` is a single matmul, so this is milliseconds where the
    full cell is seconds. Used only where the other metrics are not needed;
    ``run_cell`` remains the single definition of the metric itself.
    """
    rng = np.random.default_rng(seed)
    patterns = random_patterns(P, N, rng)
    net, _ = build(patterns, condition, lam,
                   np.random.default_rng(seed + 10_000), compute_rank=False)
    return float(net.is_fixed_point(patterns).mean())


def run_cell(condition: str, lam: float, N: int, P: int, seed: int,
             max_sweeps: int = 60, n_flip_levels: int = 10,
             n_basin_probes: int = 5, viol_sweeps: int = 4) -> Cell:
    rng = np.random.default_rng(seed)
    patterns = random_patterns(P, N, rng)
    net, info = build(patterns, condition, lam, np.random.default_rng(seed + 10_000))
    Sym = (patterns.T @ patterns) / N
    np.fill_diagonal(Sym, 0.0)

    # --- primary: do the stored patterns survive as exact fixed points? ----
    frac_fixed = float(net.is_fixed_point(patterns).mean())

    # --- does the dynamics settle at all, and where does it land? ----------
    final, _, converged = net.converge(patterns.copy(),
                                       np.random.default_rng(seed + 20_000),
                                       max_sweeps=max_sweeps)
    ov = (final @ patterns.T) / N                      # (P, P)
    mean_overlap = float(np.abs(np.diag(ov)).mean())
    p_next = float(np.diag(np.roll(ov, -1, axis=1)).mean())   # overlap with xi^{mu+1}
    # per-state non-convergence, not the net-wide flag
    still_moving = ~net.is_fixed_point(final)
    frac_nonconverged = float(still_moving.mean()) if converged is not None else 1.0

    # --- basin radius: probe 1's metric, kept commensurable -----------------
    brng = np.random.default_rng(seed + 30_000)
    radii = []
    for p in range(min(P, n_basin_probes)):
        ok_frac = []
        for k in range(1, n_flip_levels + 1):
            n_flip = int(round(N * k / (2 * n_flip_levels)))   # up to 50 % corruption
            probes = np.repeat(patterns[p:p + 1], 8, axis=0)
            for r in range(probes.shape[0]):
                flip = brng.choice(N, size=n_flip, replace=False)
                probes[r, flip] *= -1
            out, _, _ = net.converge(probes, brng, max_sweeps=max_sweeps)
            ok_frac.append(float(((out @ patterns[p]) / N > 0.95).mean()))
        radii.append(float(np.mean(ok_frac)))
    basin_radius = float(np.mean(radii))

    viol = _energy_violations(net, Sym, patterns[:8].copy(),
                              np.random.default_rng(seed + 40_000), viol_sweeps)

    return Cell(condition=condition, lam=lam, N=N, P=P, seed=seed,
                alpha=P / N, frac_fixed=frac_fixed,
                frac_nonconverged=frac_nonconverged, basin_radius=basin_radius,
                mean_overlap=mean_overlap, p_next=p_next,
                energy_violations=viol, rank_A=info["rank_A"])


def cell_dict(c: Cell) -> dict:
    return asdict(c)
