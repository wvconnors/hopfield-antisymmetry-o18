"""Pre-flight checks — run these on a metric *before* believing what it says.

Every check here exists because a specific claim in M1 had to be withdrawn. The
taxonomy and the evidence are in `docs/25_probe_protocol.md`; this module is the
mechanised form, because a checklist that is written and not run is worth
nothing.

    A. the instrument cannot move   -> dynamic_range, quantisation
    B. the control cannot fail      -> control_potency, tolerance_is_discriminative
    C. name does not match content  -> observes, depends_on
    D. the null is not null         -> empirical_null

Plus :func:`scan_table`, the cheapest check of all: read your own output for
columns that are constant, exactly 0, exactly 1, NaN or suspiciously flat. Both
of probe 3's fatal findings were sitting in the printed output, unexamined.

Everything returns a :class:`Check` rather than raising, so a probe can report
its pre-flight in `summary.json` alongside its results. ``Check.__bool__`` is the
pass/fail, so ``assert preflight.dynamic_range(...)`` reads naturally in tests.
"""

from __future__ import annotations

from dataclasses import dataclass, field, asdict
from typing import Callable, Sequence

import numpy as np

__all__ = [
    "Check", "substrate_predictability", "dynamic_range", "quantisation",
    "control_potency", "tolerance_is_discriminative", "empirical_null",
    "observes", "depends_on", "scan_table", "report",
]


@dataclass
class Check:
    name: str
    passed: bool
    detail: str
    values: dict = field(default_factory=dict)

    def __bool__(self) -> bool:
        return self.passed

    def to_dict(self) -> dict:
        return asdict(self)

    def __repr__(self) -> str:  # pragma: no cover - cosmetic
        mark = "PASS" if self.passed else "FAIL"
        return f"[{mark}] {self.name}: {self.detail}"


# --------------------------------------------------------------------------
# 0 — measure the substrate before measuring instruments on it
# --------------------------------------------------------------------------

def _blocked_folds(n: int, k: int):
    """Contiguous blocks, NOT interleaved.

    Interleaved k-fold leaks across neighbouring time bins whenever the data is
    autocorrelated — and neural activity always is. Blocked folds are the
    minimum honest choice for a time series.
    """
    edges = np.linspace(0, n, k + 1).astype(int)
    for i in range(k):
        test = np.arange(edges[i], edges[i + 1])
        train = np.concatenate([np.arange(0, edges[i]), np.arange(edges[i + 1], n)])
        yield train, test


def _ridge_r2(X: np.ndarray, Y: np.ndarray, alphas=(1e-2, 1e0, 1e2, 1e4),
              n_folds: int = 5) -> float:
    """Out-of-sample R² of the best-of-`alphas` linear map X -> Y, blocked CV."""
    n = X.shape[0]
    best = -np.inf
    for a in alphas:
        num = den = 0.0
        for tr, te in _blocked_folds(n, n_folds):
            Xt = np.hstack([X[tr], np.ones((tr.size, 1))])
            Xe = np.hstack([X[te], np.ones((te.size, 1))])
            G = Xt.T @ Xt + a * np.eye(Xt.shape[1])
            G[-1, -1] -= a
            W = np.linalg.solve(G, Xt.T @ Y[tr])
            resid = Y[te] - Xe @ W
            num += float(np.sum(resid ** 2))
            den += float(np.sum((Y[te] - Y[tr].mean(axis=0)) ** 2))
        best = max(best, 1.0 - num / den if den > 0 else -np.inf)
    return float(best)


def substrate_predictability(X: np.ndarray, n_lags: int = 1, n_folds: int = 5,
                             n_surrogates: int = 3, seed: int = 0,
                             min_excess: float = 0.005,
                             name: str = "substrate_predictability") -> Check:
    """**Run this FIRST, before building any instrument.** Is there anything here?

    Fits the best linear predictor of ``x_{t+1}`` from ``n_lags`` of history
    (blocked CV) and compares its out-of-sample R² against rate-matched
    independent surrogates. If the real data is not measurably more predictable
    than surrogates that share only its firing rates, **the substrate has no
    temporal structure to encode**, and any probe asking "which representation
    predicts best" is void before it starts.

    ``X`` is ``(n_units, n_timebins)``.

    Why this exists
    ---------------
    Probe 4 was fully specified — encoders, matched code size, three falsifiers,
    a mandatory preflight — against cached Brunel activity in which predicting
    ``x_{t+1}`` from the **entire observed state** scores *worse* than predicting
    from nothing. 512 of 12 500 neurons were recorded (4 %), so most of the
    causal state was unobserved and each neuron received fresh external Poisson
    noise every step. The probe's own falsifier would have fired by
    construction. This check takes ~30 s and would have stopped it before the
    spec was written.

    The general lesson: **measure the substrate before measuring instruments on
    it.** Every other check in this module validates an instrument. This one
    asks whether the thing being measured exists.
    """
    X = np.asarray(X, dtype=float)
    if X.ndim != 2 or X.shape[1] < 4 * n_folds:
        raise ValueError("X must be (n_units, n_timebins) with enough bins")
    rng = np.random.default_rng(seed)

    def score(M: np.ndarray) -> float:
        T = M.shape[1]
        feats = np.concatenate([M[:, k:T - n_lags + k].T for k in range(n_lags)],
                               axis=1)
        target = M[:, n_lags:].T
        return _ridge_r2(feats, target, n_folds=n_folds)

    real = score(X)
    rates = X.mean(axis=1, keepdims=True)
    surr = []
    for _ in range(n_surrogates):
        # rate-matched, temporally independent: preserves per-unit mean only
        surr.append(score((rng.random(X.shape) < rates).astype(float)))
    s_m = float(np.mean(surr))
    s_sd = float(np.std(surr, ddof=1)) if len(surr) > 1 else 0.0

    excess = real - s_m
    ok = excess > max(min_excess, 3 * s_sd)
    return Check(name, ok,
                 f"linear-predictability R² {real:.5f} vs rate-matched surrogate "
                 f"{s_m:.5f}±{s_sd:.5f} (excess {excess:+.5f})"
                 + ("" if ok else " — SUBSTRATE HAS NO PREDICTABLE STRUCTURE; "
                    "any probe measuring 'which code predicts best' is void"),
                 {"real_r2": real, "surrogate_r2": s_m, "surrogate_sd": s_sd,
                  "excess": excess, "n_lags": n_lags,
                  "n_units": int(X.shape[0]), "n_bins": int(X.shape[1])})


# --------------------------------------------------------------------------
# A — can the instrument move?
# --------------------------------------------------------------------------

def dynamic_range(measure: Callable[[np.random.Generator], float],
                  null: Callable[[np.random.Generator], float],
                  saturated: Callable[[np.random.Generator], float],
                  observed: float, n_seeds: int = 5, n_sigma: float = 3.0,
                  name: str = "dynamic_range") -> Check:
    """Is ``observed`` distinguishable from both the null and saturated endpoints?

    ``null`` and ``saturated`` take an rng and return the metric computed on
    data where the answer is known to be "nothing" and "everything"
    respectively. **They must be evaluated at the experiment's real N, T and
    parameters** — that is the whole point; a metric with plenty of range at
    N=100 can be pinned at N=500.

    Fails if ``observed`` sits within ``n_sigma`` of either endpoint, or if the
    two endpoints are themselves indistinguishable (in which case the metric has
    no range at all and nothing it reports can mean anything).

    Caught in M1: ``spurious_total`` saturated at 1.000 with genuine basin
    volume floored at 0.000; probe 3's separability numerator pinned by its own
    normalisation.
    """
    rng = np.random.default_rng(0)
    lo = np.array([float(null(rng)) for _ in range(n_seeds)])
    hi = np.array([float(saturated(rng)) for _ in range(n_seeds)])
    obs = float(observed)

    lo_m, lo_s = float(lo.mean()), float(lo.std(ddof=1) if n_seeds > 1 else 0.0)
    hi_m, hi_s = float(hi.mean()), float(hi.std(ddof=1) if n_seeds > 1 else 0.0)
    spread = abs(hi_m - lo_m)
    noise = max(lo_s, hi_s, 1e-12)

    if spread < n_sigma * noise:
        return Check(name, False,
                     f"metric has NO range: null {lo_m:.4g}±{lo_s:.2g} vs "
                     f"saturated {hi_m:.4g}±{hi_s:.2g}",
                     {"null": lo_m, "saturated": hi_m, "observed": obs})

    d_lo = abs(obs - lo_m) / max(lo_s, spread * 1e-6)
    d_hi = abs(obs - hi_m) / max(hi_s, spread * 1e-6)
    ok = d_lo > n_sigma and d_hi > n_sigma
    where = ("pinned at the NULL endpoint" if d_lo <= n_sigma
             else "pinned at the SATURATED endpoint" if d_hi <= n_sigma
             else "has room to move")
    return Check(name, ok,
                 f"observed {obs:.4g} vs null {lo_m:.4g}±{lo_s:.2g} / "
                 f"saturated {hi_m:.4g}±{hi_s:.2g} — {where}",
                 {"null": lo_m, "null_sd": lo_s, "saturated": hi_m,
                  "saturated_sd": hi_s, "observed": obs,
                  "sigma_from_null": d_lo, "sigma_from_saturated": d_hi})


def quantisation(measure: Callable[[np.ndarray], float],
                 inputs: Sequence[np.ndarray], min_distinct: int = 3,
                 name: str = "quantisation") -> Check:
    """Does the metric return distinct values for inputs that genuinely differ?

    A statistic that collapses a wide range of inputs onto one value cannot
    detect anything in that range, however well it behaves elsewhere.

    Caught in M1: an autocorrelation statistic defined as "lags until r falls
    below 1/e", which returns the integer 1 for *every* series with lag-1
    r < 0.368 — including strongly anti-correlated ones. It mapped [−1, 0.368]
    to a single value and destroyed a real result (lag-1 r ran −0.012 → +0.612
    across cluster sizes).
    """
    vals = [float(measure(x)) for x in inputs]
    finite = [v for v in vals if np.isfinite(v)]
    n_distinct = len({round(v, 12) for v in finite})
    ok = n_distinct >= min_distinct
    return Check(name, ok,
                 f"{n_distinct} distinct values from {len(inputs)} distinct "
                 f"inputs: {[round(v, 4) for v in vals]}",
                 {"values": vals, "n_distinct": n_distinct})


# --------------------------------------------------------------------------
# B — can the control fail?
# --------------------------------------------------------------------------

def control_potency(measure: Callable[[object], float],
                    apply_control: Callable[[object], object],
                    planted: object, baseline: float | None = None,
                    min_effect: float = 0.0, name: str = "control_potency") -> Check:
    """Can the control produce the outcome *at all*?

    ``planted`` is data in which the effect is present by construction. Applies
    the control and checks the metric actually moves. **A control that cannot
    produce the outcome is not a control** — passing it is guaranteed and
    carries no evidence.

    Caught in M1: weight decay under Hopfield's ``sign(Ws)`` scale invariance
    (provably cannot move a single fixed point); a "matched" shuffled control
    delivering 7.5× less plasticity than the condition it was matched to.
    """
    before = float(measure(planted)) if baseline is None else float(baseline)
    after = float(measure(apply_control(planted)))
    delta = abs(after - before)
    ok = delta > max(min_effect, 1e-12)
    return Check(name, ok,
                 f"control moved the metric by {delta:.4g} "
                 f"({before:.4g} -> {after:.4g})"
                 + ("" if ok else " — INERT, so it rules nothing out"),
                 {"before": before, "after": after, "delta": delta})


def tolerance_is_discriminative(tolerance: float, observed: Sequence[float],
                                name: str = "tolerance") -> Check:
    """Is the pass tolerance smaller than the values being compared?

    If a test asserts ``|a − b| < tol`` and ``tol`` exceeds the typical
    magnitude of ``a`` and ``b``, the test passes for a degenerate implementation
    that returns zero.

    Caught in M1: a Brian2 cross-check asserting ``|CV_ours − CV_brian| < 0.15``
    against measured CVs of 0.09–0.14. A simulator emitting perfectly clock-like
    spike trains would have passed.
    """
    mags = np.abs(np.asarray([v for v in observed if np.isfinite(v)], dtype=float))
    typical = float(np.median(mags)) if mags.size else 0.0
    ok = tolerance < typical
    return Check(name, ok,
                 f"tolerance {tolerance:.4g} vs typical magnitude {typical:.4g}"
                 + ("" if ok else " — a degenerate (all-zero) result would pass"),
                 {"tolerance": float(tolerance), "typical": typical})


# --------------------------------------------------------------------------
# C — does the name match the content?
# --------------------------------------------------------------------------

def observes(target: object, attr: str, run: Callable[[], object],
             name: str = "observes") -> Check:
    """Does ``run`` actually exercise ``target.attr``?

    Replaces the attribute with something that raises, then checks ``run``
    fails. If it still succeeds, ``run`` never touched the thing it is named
    after. This is a mutation test in miniature and it is the single highest-
    yield check here.

    Caught in M1: an f–I test that built a network, silenced it, then hand-rolled
    the integration loop *inside the test body* — it would have passed with
    ``BrunelNet.run`` deleted, while its docstring claimed it validated the
    integrator's leak, threshold, reset and refractory handling.
    """
    original = getattr(target, attr)

    def _boom(*a, **k):
        raise AssertionError(f"{attr} should have been called")

    setattr(target, attr, _boom)
    try:
        run()
        reached = True
    except Exception:
        reached = False
    finally:
        setattr(target, attr, original)

    return Check(name, not reached,
                 f"{getattr(target, '__name__', target)}.{attr} is "
                 + ("exercised" if not reached
                    else "NEVER CALLED — the check does not test what it names"),
                 {"attribute": attr})


def depends_on(fn: Callable, kwargs: dict, varying: str,
               values: Sequence, name: str = "depends_on") -> Check:
    """Does ``fn``'s output actually change when ``varying`` changes?

    A function whose name promises to account for some quantity must respond to
    it. If the output is constant across ``values``, the quantity is not in the
    computation whatever the docstring says.

    Caught in M1: ``mean_input_mv`` documented as the total mean input
    "including recurrent" but omitting the recurrent term entirely — it was
    identically ``ν_ratio · θ``, so the recorded field read 40.0 mV where the
    true value was 21.1.
    """
    outs = []
    for v in values:
        outs.append(float(fn(**{**kwargs, varying: v})))
    n_distinct = len({round(o, 12) for o in outs})
    ok = n_distinct > 1
    return Check(name, ok,
                 f"varying {varying!r} over {list(values)} gave "
                 f"{[round(o, 4) for o in outs]}"
                 + ("" if ok else f" — output IGNORES {varying!r}"),
                 {"values": list(map(float, values)), "outputs": outs})


# --------------------------------------------------------------------------
# D — is the null actually null?
# --------------------------------------------------------------------------

def empirical_null(pipeline: Callable[[np.random.Generator], float],
                   n_seeds: int = 8, analytic: float | None = None,
                   n_sigma: float = 3.0, name: str = "empirical_null") -> Check:
    """Run the **whole pipeline** on no-effect data and report mean ± SD.

    Pass ``analytic`` to test whether the estimator is unbiased. Fails when the
    measured null differs from the analytic value by more than ``n_sigma`` —
    which means excesses must be judged against the *measured* null, not the
    ideal one.

    Caught in M1: greedy max-correlation pairing has a winner's curse (it selects
    the largest sample correlations, which among independent units average
    r ≈ 0.07), putting the pure-noise variance exponent at 1.033 ± 0.002 —
    16σ from the analytic 1.0. Random pairing recovers 1.000.
    """
    rng = np.random.default_rng(0)
    vals = np.array([float(pipeline(rng)) for _ in range(n_seeds)])
    m = float(vals.mean())
    sd = float(vals.std(ddof=1)) if n_seeds > 1 else 0.0
    out = {"null_mean": m, "null_sd": sd, "n_seeds": n_seeds}
    if analytic is None:
        return Check(name, True, f"empirical null {m:.4g} ± {sd:.2g}", out)
    z = abs(m - analytic) / max(sd, 1e-12)
    out.update({"analytic": float(analytic), "z": z})
    ok = z <= n_sigma
    return Check(name, ok,
                 f"empirical null {m:.4g} ± {sd:.2g} vs analytic "
                 f"{analytic:.4g} — {z:.1f}σ"
                 + ("" if ok else "; judge excesses against the MEASURED null"),
                 out)


# --------------------------------------------------------------------------
# the cheapest check — read your own output
# --------------------------------------------------------------------------

def scan_table(rows: Sequence[dict], flat_cv: float = 0.02,
               name: str = "scan_table") -> Check:
    """Flag columns that are constant, all-0, all-1, all-NaN or suspiciously flat.

    Run this on the results table *before* writing prose about it. In probe 3
    both fatal findings were visible here and unexamined: a numerator column
    reading ``2.143, 2.137, 2.124, 2.106…`` (pinned by a normalisation), and
    ``alpha_autocorr = 0.00`` in four of five conditions with ``r² = NaN``.

    ``flat_cv`` is the coefficient of variation below which a column counts as
    suspiciously flat.
    """
    if not rows:
        return Check(name, True, "no rows", {})
    keys = [k for k in rows[0]
            if all(isinstance(r.get(k), (int, float, np.floating, np.integer))
                   and not isinstance(r.get(k), bool) for r in rows)]
    flags = {}
    for k in keys:
        col = np.array([float(r[k]) for r in rows], dtype=float)
        if np.all(np.isnan(col)):
            flags[k] = "all NaN"
        elif np.any(np.isnan(col)):
            flags[k] = f"{int(np.isnan(col).sum())}/{col.size} NaN"
        elif np.allclose(col, 0.0):
            flags[k] = "exactly 0 everywhere"
        elif np.allclose(col, 1.0):
            flags[k] = "exactly 1 everywhere"
        elif col.std() == 0.0:
            flags[k] = f"constant at {col[0]:.4g}"
        elif abs(col.mean()) > 0 and col.std() / abs(col.mean()) < flat_cv:
            flags[k] = (f"suspiciously flat (CV "
                        f"{col.std() / abs(col.mean()):.3f}, mean {col.mean():.4g})")
    ok = not flags
    detail = "no degenerate columns" if ok else "; ".join(
        f"{k}: {v}" for k, v in flags.items())
    return Check(name, ok, detail, {"flags": flags})


def report(checks: Sequence[Check]) -> str:
    """Format checks for a log or a `summary.json` preamble."""
    lines = [repr(c) for c in checks]
    n_fail = sum(1 for c in checks if not c)
    lines.append(f"-- {len(checks) - n_fail}/{len(checks)} passed"
                 + ("" if not n_fail else f"; {n_fail} FAILED — do not interpret"))
    return "\n".join(lines)
