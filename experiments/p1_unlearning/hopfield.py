"""A minimal Hopfield network with the four unlearning conditions.

Pure numpy. The one non-obvious implementation choice is that convergence is
**vectorised over initial states, not over neurons**: updates stay strictly
asynchronous (one neuron at a time, fresh random order each sweep, which is
what guarantees the energy decreases monotonically), while many independent
states are carried through that same update in parallel. Converging 2000 states
therefore costs the same number of Python-level iterations as converging one.

Synchronous updates would be simpler and faster still, and are rejected: they
admit period-2 limit cycles that are not fixed points, which would contaminate
exactly the "spurious attractor" count this probe measures.
"""

from __future__ import annotations

import numpy as np

__all__ = ["HopfieldNet", "random_patterns"]


def random_patterns(n_patterns: int, N: int, rng: np.random.Generator) -> np.ndarray:
    """``(n_patterns, N)`` array of unbiased ±1 patterns."""
    return rng.choice(np.array([-1.0, 1.0]), size=(n_patterns, N))


class HopfieldNet:
    """Binary (±1) Hopfield network with Hebbian storage and anti-Hebbian unlearning.

    ``W`` is kept symmetric with zero diagonal at all times — the condition under
    which asynchronous dynamics are guaranteed to descend the energy
    ``E = −½ sᵀWs`` and terminate at a fixed point. Every operation here
    preserves it, and :meth:`assert_valid` checks it.
    """

    def __init__(self, patterns: np.ndarray):
        self.patterns = np.asarray(patterns, dtype=float)
        self.n_patterns, self.N = self.patterns.shape
        self.W = (self.patterns.T @ self.patterns) / self.N
        np.fill_diagonal(self.W, 0.0)
        self.W0_norm = float(np.linalg.norm(self.W))

    # -- properties ---------------------------------------------------------

    @property
    def load(self) -> float:
        """α = P/N. The Hopfield capacity limit is α_c ≈ 0.138."""
        return self.n_patterns / self.N

    def frobenius(self) -> float:
        return float(np.linalg.norm(self.W))

    def assert_valid(self, tol: float = 1e-9) -> None:
        assert np.max(np.abs(self.W - self.W.T)) < tol, "W lost symmetry"
        assert np.max(np.abs(np.diag(self.W))) < tol, "W gained a diagonal"

    def energy(self, S: np.ndarray) -> np.ndarray:
        """Energy of each row of ``S``: ``E = −½ sᵀWs``."""
        S = np.atleast_2d(np.asarray(S, dtype=float))
        return -0.5 * np.einsum("ij,jk,ik->i", S, self.W, S)

    # -- dynamics -----------------------------------------------------------

    def converge(self, S: np.ndarray, rng: np.random.Generator,
                 max_sweeps: int = 100):
        """Run asynchronous dynamics to a fixed point.

        ``S`` is ``(n_states, N)``; returns ``(fixed_points, n_sweeps, all_converged)``.

        All states share the neuron order within a sweep (a fresh random
        permutation each sweep). They are independent, so a shared order costs
        nothing in correctness and buys the vectorisation.
        """
        S = np.array(np.atleast_2d(S), dtype=float)
        sweeps = 0
        for sweeps in range(1, max_sweeps + 1):
            changed = False
            for i in rng.permutation(self.N):
                h = S @ self.W[:, i]
                new = np.where(h >= 0, 1.0, -1.0)
                if np.any(new != S[:, i]):
                    S[:, i] = new
                    changed = True
            if not changed:
                return S, sweeps, True
        return S, sweeps, False

    def is_fixed_point(self, S: np.ndarray) -> np.ndarray:
        """Boolean per row: is every neuron already aligned with its field?"""
        S = np.atleast_2d(np.asarray(S, dtype=float))
        H = S @ self.W
        return np.all(np.where(H >= 0, 1.0, -1.0) == S, axis=1)

    # -- plasticity ---------------------------------------------------------

    def anti_hebbian(self, S: np.ndarray, epsilon: float) -> None:
        """``W ← W − (ε/N) Σ_b s_b s_bᵀ`` over the rows of ``S`` (diagonal re-zeroed).

        This is the negative / free-running phase of energy-based learning:
        the same term as the Boltzmann machine's ``−⟨s_i s_j⟩_model`` and as
        Hopfield–Feinstein–Palmer unlearning (fragment 15, I9).
        """
        S = np.atleast_2d(np.asarray(S, dtype=float))
        self.W -= (epsilon / self.N) * (S.T @ S)
        np.fill_diagonal(self.W, 0.0)

    def decay(self, delta: float) -> None:
        """``W ← (1−δ)W`` — the boring alternative explanation, made runnable."""
        self.W *= (1.0 - delta)

    def decay_to_norm(self, target_norm: float) -> float:
        """Scale ``W`` to a target Frobenius norm; returns the δ applied.

        Used by the weight-decay control to track the sleep condition's norm
        trajectory *exactly*, so that "shrinking the weights" is ruled out as
        the explanation rather than merely gestured at.
        """
        cur = self.frobenius()
        if cur <= 0:
            return 0.0
        delta = 1.0 - target_norm / cur
        self.decay(delta)
        return float(delta)
