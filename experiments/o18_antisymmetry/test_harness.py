"""Ground-truth tests for the O18 harness.

Every test here is a **name-content** test in the sense of `docs/25`: it fails
if the named thing computed something else. Probe 2's f-I test would have
passed with `run()` deleted; these are written so they cannot.
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from experiments.o18_antisymmetry.harness import (  # noqa: E402
    antisym, stdp_sequence_A, random_A, rank_matched_random_A, build, run_cell,
)
from experiments.p1_unlearning.hopfield import random_patterns  # noqa: E402


def test_operators_are_antisymmetric():
    rng = np.random.default_rng(0)
    pats = random_patterns(10, 60, rng)
    for A in (stdp_sequence_A(pats), random_A(60, rng), rank_matched_random_A(pats, rng)):
        assert np.max(np.abs(A + A.T)) < 1e-12
        assert np.max(np.abs(np.diag(A))) < 1e-12


def test_energy_ignores_antisymmetric_part():
    """sᵀAs = 0 identically, so `energy` reports the SYMMETRIC-part energy.

    The spec leans on this to define "energy-increasing update". Asserting it
    rather than assuming it is the whole point.
    """
    rng = np.random.default_rng(1)
    pats = random_patterns(8, 50, rng)
    net0, _ = build(pats, "structured", 0.0, np.random.default_rng(2))
    states = random_patterns(20, 50, rng)
    e_sym = net0.energy(states)
    for lam in (0.3, 1.0, 5.0):
        net, _ = build(pats, "structured", lam, np.random.default_rng(2))
        assert np.allclose(net.energy(states), e_sym, atol=1e-9)


def test_norm_ratio_is_exactly_lambda():
    """λ IS ‖A‖/‖S‖ by construction. If this drifts, P1 is manufactured."""
    rng = np.random.default_rng(3)
    pats = random_patterns(20, 100, rng)
    for cond in ("structured", "random", "rank_matched_random", "shuffled_order"):
        for lam in (0.0, 0.15, 0.5, 1.0):
            net, info = build(pats, cond, lam, np.random.default_rng(4))
            S = (pats.T @ pats) / 100
            np.fill_diagonal(S, 0.0)
            A = net.W - S
            got = np.linalg.norm(A) / np.linalg.norm(S)
            assert abs(got - lam) < 1e-8, (cond, lam, got)


def test_structured_is_low_rank_and_random_is_not():
    """The rank/alignment premise of the §3 prediction, checked not assumed."""
    rng = np.random.default_rng(5)
    P, N = 20, 200
    pats = random_patterns(P, N, rng)
    assert np.linalg.matrix_rank(stdp_sequence_A(pats)) <= 2 * P
    assert np.linalg.matrix_rank(rank_matched_random_A(pats, rng)) <= 2 * P
    assert np.linalg.matrix_rank(random_A(N, rng)) > 4 * P


def test_structured_field_points_at_the_next_pattern():
    """A_seq applied to ξ^μ must produce a field aligned with ξ^{μ+1}.

    This is the mechanism the whole probe rests on. A test that only checked
    antisymmetry would pass on an operator that chained the WRONG patterns.
    """
    rng = np.random.default_rng(6)
    P, N = 10, 2000
    pats = random_patterns(P, N, rng)
    A = stdp_sequence_A(pats)
    xtalk = np.sqrt(P / N)          # the analytic crosstalk scale
    for mu in range(P):
        h = A @ pats[mu]
        nxt = float(h @ pats[(mu + 1) % P]) / N
        prv = float(h @ pats[(mu - 1) % P]) / N
        others = [float(h @ pats[k]) / N for k in range(P)
                  if k not in ((mu + 1) % P, (mu - 1) % P)]
        # antisymmetrisation halves the chain coefficient: 1/2 in the large-N
        # limit, minus O(1/sqrt(N)) crosstalk.
        assert nxt > 0.35, (mu, nxt)
        # the STDP window's defining SHAPE: potentiate forward, depress backward,
        # balanced. An operator that chained the wrong pairs, or that lost the
        # depression lobe, passes the antisymmetry test and fails this one.
        assert prv < -0.35, (mu, prv)
        assert abs(nxt + prv) < 0.15 * nxt, (mu, nxt, prv)
        # off-chain overlaps must sit at the analytic crosstalk scale sqrt(P/N),
        # not merely "smaller than the signal" -- a fixed ratio would silently
        # become untestable as N grows.
        assert max(abs(o) for o in others) < 4 * xtalk, (mu, max(others), xtalk)
        assert nxt > 4 * max(abs(o) for o in others), (mu, nxt, max(others))


def test_shuffled_order_chains_different_pairs():
    """The control must actually differ from the thing it controls for."""
    rng = np.random.default_rng(7)
    pats = random_patterns(12, 80, rng)
    a = stdp_sequence_A(pats)
    b = stdp_sequence_A(pats, order=np.random.default_rng(8).permutation(12))
    assert np.linalg.norm(a - b) / np.linalg.norm(a) > 0.5


def test_lambda_zero_recovers_probe1_symmetric_net():
    """λ=0 must be the untouched symmetric Hopfield net, exactly."""
    rng = np.random.default_rng(9)
    pats = random_patterns(15, 150, rng)
    for cond in ("structured", "random", "rank_matched_random"):
        net, _ = build(pats, cond, 0.0, np.random.default_rng(10))
        net.assert_valid()


def test_run_cell_actually_calls_the_dynamics():
    """A metric that cannot see the integrator is not a measurement.

    Breaking `converge` must break the cell. Probe 2's f-I test failed exactly
    this and nobody noticed for a session.
    """
    import experiments.p1_unlearning.hopfield as hop

    good = run_cell("structured", 0.5, N=120, P=12, seed=0,
                    max_sweeps=20, n_flip_levels=3, n_basin_probes=2, viol_sweeps=1)

    orig = hop.HopfieldNet.converge
    hop.HopfieldNet.converge = lambda self, S, rng, max_sweeps=100: (
        np.array(np.atleast_2d(S), dtype=float), 0, True)
    try:
        broken = run_cell("structured", 0.5, N=120, P=12, seed=0,
                          max_sweeps=20, n_flip_levels=3, n_basin_probes=2,
                          viol_sweeps=1)
    finally:
        hop.HopfieldNet.converge = orig

    assert good.basin_radius != broken.basin_radius


def test_energy_violations_matches_bruteforce():
    """The O(N) delta-energy shortcut must count EXACTLY what re-evaluating
    ``E = −½ sᵀ Sym s`` would count. An optimisation that quietly changed the
    quantity is the name-content failure `docs/25` catalogues."""
    from experiments.o18_antisymmetry.harness import _energy_violations

    rng = np.random.default_rng(11)
    N, P = 80, 8
    pats = random_patterns(P, N, rng)
    Sym = (pats.T @ pats) / N
    np.fill_diagonal(Sym, 0.0)
    net, _ = build(pats, "structured", 0.6, np.random.default_rng(12))

    fast = _energy_violations(net, Sym, pats[:4].copy(),
                              np.random.default_rng(13), n_sweeps=3)

    # brute force: recompute the quadratic form on both sides of every update
    S = np.array(pats[:4], dtype=float)
    viol = np.zeros(S.shape[0])
    r = np.random.default_rng(13)
    for _ in range(3):
        for i in r.permutation(N):
            before = -0.5 * np.einsum("ij,jk,ik->i", S, Sym, S)
            S[:, i] = np.where(S @ net.W[:, i] >= 0, 1.0, -1.0)
            after = -0.5 * np.einsum("ij,jk,ik->i", S, Sym, S)
            viol += (after > before + 1e-9)
    assert fast == pytest.approx(float(viol.mean()))


def test_frac_fixed_only_agrees_with_run_cell():
    """The cheap pre-flight path must compute the SAME number as the real one."""
    from experiments.o18_antisymmetry.harness import frac_fixed_only
    for cond in ("structured", "random"):
        for lam in (0.0, 0.1, 0.5):
            full = run_cell(cond, lam, N=150, P=15, seed=3, max_sweeps=10,
                            n_flip_levels=2, n_basin_probes=1, viol_sweeps=1)
            assert frac_fixed_only(cond, lam, 150, 15, 3) == pytest.approx(full.frac_fixed)


@pytest.mark.parametrize("cond", ["structured", "random"])
def test_frac_fixed_has_dynamic_range(cond):
    """Primary metric must move between the null and the saturated end.

    NOTE the λ=0 baseline is ~0.90, not 1.0: at α=0.1 the plain Hopfield net is
    close enough to α_c≈0.138 that a couple of stored patterns are already not
    exact fixed points. **That is the null floor the sweep must be read
    against** — an excess measured against 1.0 would be measuring the load, not
    the antisymmetry.
    """
    lo = run_cell(cond, 0.0, N=200, P=20, seed=0, max_sweeps=20,
                  n_flip_levels=2, n_basin_probes=1, viol_sweeps=1)
    hi = run_cell(cond, 8.0, N=200, P=20, seed=0, max_sweeps=20,
                  n_flip_levels=2, n_basin_probes=1, viol_sweeps=1)
    assert lo.frac_fixed >= 0.85
    assert hi.frac_fixed < lo.frac_fixed - 0.2


def test_frac_fixed_curve_agrees_with_the_slow_path():
    """The lambda-sweep fast path must reproduce `frac_fixed_only` exactly.

    It caches S and A across lambda; a caching bug would shift the whole curve
    and, with it, every knee -- i.e. the entire P1 result.
    """
    from experiments.o18_antisymmetry.harness import frac_fixed_curve, frac_fixed_only
    lams = [0.0, 0.05, 0.17, 0.5, 1.3]
    for cond in ("structured", "random", "rank_matched_random", "shuffled_order"):
        fast = frac_fixed_curve(cond, lams, 200, 20, seed=5)
        slow = [frac_fixed_only(cond, l, 200, 20, 5) for l in lams]
        assert fast == pytest.approx(slow), cond
