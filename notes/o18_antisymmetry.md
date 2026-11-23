# O18 — how much antisymmetry can the energy landscape survive?

*Spec `docs/42_o18_antisymmetry_spec.md` · code `experiments/o18_antisymmetry/` · results `results/o18_antisymmetry/` · **mode: intuition-build, not contribution** (§0 of the spec, and §8 below).*

> **Version 2.** V1 reported that pattern-aligned antisymmetry destroys fixed points at **≈ ⅓** the perturbation isotropic antisymmetry needs, that the pre-registered `√α` law was **rejected** in favour of `α^0.216`, and that this warranted a scoped downgrade of the energy function's biological-support rating. Adversarial review — three passes, findings re-verified independently before any of this was written — showed that **the ratio is a unit conversion**, that **`√α` is the correct exponent for the only part of the effect that exists**, that **one of the two controls could not have failed**, and that **the P2 sequence-transition arm measured the wrong lobe of its own operator**. §7 keeps the failure record, because how the errors happened is the most transferable part of the probe.
>
> **What survives unchanged:** the arithmetic (every number V1 quoted from the main arm reproduces bit-for-bit from the recorded cells), the negative conclusion that biological *structure* is not what matters, and the reproduction's competence. What changes is nearly every sentence that was headed for a register.

---

## 1. Verdict

**In the units the probe swept — `λ = ‖A‖_F/‖S‖_F` — a pattern-aligned antisymmetric operator destroys the stored patterns' fixed points at ≈ ⅓ the λ an isotropic one needs. That entire factor is the Frobenius norm mis-measuring what the operator does: at matched `‖A‖_F`, an operator confined to the P-dimensional pattern subspace delivers `√(N/P)` times more field to the stored patterns. Matched on delivered field, all four conditions — isotropic random included — break at the same place.**

| condition | rank | field/‖A‖_F | knee (λ, 5 seeds) | **knee × field** |
|---|---|---|---|---|
| structured (STDP on a sequence) | 48 | 0.1477 | 0.1648 ± 0.0169 | **0.0243** |
| shuffled_order | 48 | 0.1483 | 0.1546 ± 0.0096 | **0.0229** |
| rank_matched_random | 50 | 0.1484 | 0.1640 ± 0.0060 | **0.0243** |
| isotropic random | 500 | 0.0448 | 0.5167 ± 0.0425 | **0.0231** |

*(N=500, P=50, α=0.1; knee = λ at which `frac_fixed` reaches half its λ=0 value, interpolated on a Δλ=0.0025 grid. `field/‖A‖_F` is `RMS_{i,μ}(Aξ^μ)_i / ‖A‖_F`. Measured field ratio aligned:isotropic = 3.30; `√(N/P)` = 3.16.)*

`knee × field` is flat to **±3 %** across conditions that differ in rank by a factor of ten. There is no residual alignment effect at α = 0.1. The measured phenomenon is one sentence, and it is a statement about units, not about the landscape.

## 2. Verdict on the pre-registered predictions

- **P1 (direction) — ✅ confirmed in Frobenius units, and the effect is a normalisation.** Knee 0.165 (structured) vs 0.517 (isotropic); ratio 0.32. Real, reproducible, and fully accounted for by §1. Restated: *the landscape's tolerance for antisymmetry, expressed in `‖A‖_F`, depends on where the operator's support lies; expressed in delivered field, it does not.*
- **P1 (the number) — the rejection is withdrawn; `√α` was right about the part that exists.** The Frobenius knee ratio fits `0.518 · α^0.216` over seven loads inside capacity — that fit reproduces exactly. But it **decomposes**:

  | quantity | fitted law |
  |---|---|
  | unit conversion (isotropic field ÷ aligned field, per unit ‖A‖_F) | `0.92 · α^{+0.48}` — i.e. `√α`, the spec's own derivation |
  | residual, matched on delivered field | `0.567 · α^{−0.26}` |
  | product = Frobenius knee ratio | `0.518 · α^{+0.216}` |

  So `0.216 = 0.48 − 0.26`. The pre-registered exponent 0.5 **is** the exponent of the unit conversion, to within the fit's error. The residual runs the *other way*: matched on delivered field, aligned antisymmetry is **less** damaging than isotropic, and increasingly so as load falls (field-matched knee ratio 1.04 at α = 0.10 → 1.60 at α = 0.02). V1 called P1's number "decisively rejected"; the truth is that the spec's mechanism was right and its bookkeeping was incomplete.
- **P1's uncertainty was quoted ~10× too tight.** V1 offered a leave-one-out range over *loads* (0.210–0.224) as if it were an interval. Fitting each seed separately — the seeds are where the noise is — gives per-seed exponents `[0.324, 0.148, 0.220, 0.141, 0.238]`, **0.214 ± 0.075**. A seed bootstrap agrees (95 % CI ≈ [0.15, 0.28]). The Frobenius exponent is excluded from 0.5 at ≈ 4σ, not the ≈ 40σ the LOO range implies. Quote it as **≈ 0.21 ± 0.08**.
- **P2 (sequence transition) — ✗ WITHDRAWN, twice over.** See §5. It measured the wrong direction, and it was arithmetically unreachable in this design regardless.
- **P3 (knee falls with load) — ✅ confirmed**, and it is the one prediction that needs no restatement: the structured knee runs 0.504 → 0.165 over α 0.02 → 0.10, fitting `0.038 · α^{−0.68}`.

## 3. The controls: one fired, one could not

**`rank_matched_random` fired, and it is the probe's real result.** A random operator confined to the pattern subspace is indistinguishable from the STDP-on-a-sequence operator at every λ. So the STDP window and the sequence contribute nothing beyond the support. **The biological specificity the probe was reaching for is not there** — that conclusion is unchanged and survived every attack.

But two caveats the review attached to it, both material:

- **It is not rank-matched.** Rank is **50 vs 48** in every recorded cell (a P-cycle's `R − Rᵀ` has two zero eigenvalues; `G − Gᵀ` does not), spectral norm differs by 17 %, and the participation ratio of `σ²` is 23.6 vs 30.0. What it matches, to 0.05 %, is `field/‖A‖_F` — the one quantity §1 shows the knee depends on. It is a **field-matched** control that has been named for two properties it does not control. `rank_A` is a recorded column in every JSONL row and `analyze.py` never prints it, which is why nobody looked.
- Consequently V1's conclusion — *"what matters is that the perturbation lives in the pattern subspace and is low rank"* — is **not supported by this control**, because rank never varied independently of delivered field. What is supported is *what matters is delivered field*, which is simpler and stronger. Deciding whether rank matters at all needs rank-`r` truncations at constant `RMS(Aξ)`; that is cheap and is the obvious follow-on.

**`shuffled_order` could not have failed. Withdraw it as evidence.**

```
max| stdp_sequence_A(pat, order=σ) − stdp_sequence_A(pat[σ]) | = 0.0
max| S(pat) − S(pat[σ]) |                                     = 0.0
```

`S = ξᵀξ/N` is invariant under permuting the pattern rows and `frac_fixed` averages over the same *set*, so the `shuffled_order` cell **is** the `structured` cell computed on `patterns[σ]`. With i.i.d. patterns, `patterns[σ]` is an exact re-draw from the same distribution: the joint law of `(S, A)` is identical between the two conditions, and every permutation-invariant statistic has the same distribution *by construction*. Its expected difference from `structured` is exactly zero before the run. This is the repo's own **control-can't-fail** type in its purest form.

The deeper reason is worth keeping: with exchangeable patterns there is no temporal structure for an order-shuffle to destroy — "order" is a labelling of interchangeable objects. A potent order control needs either sequentially correlated patterns or a statistic that is not permutation-invariant. `p_next` was the one such statistic in the probe, and §5 is what happened to it.

## 4. What "identical" is worth at five seeds

V1 wrote *"Knees: 0.167 · 0.162 · 0.167. **Identical.** … not one resolvable cell."* The design's minimum detectable difference on the knee is **≈ 22–25 % of the knee** (80 % power, n = 5), and at λ = 0.15 — inside the knee region — MDE on `frac_fixed` is ≈ 39 % relative. At 40 seeds the structured − rank-matched gap is +0.011 ± 0.004 (t = 2.7): resolvable, and an order of magnitude too small to rescue P1. Three of the "every λ" rows (λ ≥ 0.5) have **both** arms floored at `frac_fixed = 0`, so they could not have shown a difference and should not be counted as evidence of none.

**Correct claim: the conditions agree to within ~10 % (95 % CI, 40 seeds) against a predicted 3× separation.** Not "identical". The pairing shares patterns and seeds but not the operator draw, so the residual noise is the control's own randomness — averaging several operator draws per seed would cut it as `1/√n` at negligible cost (`is_fixed_point` is one matmul).

## 5. P2 is withdrawn: the metric pointed at the lobe the dynamics pushes away from

`hopfield.converge` and `is_fixed_point` compute the field as

```python
h = S @ self.W[:, i]      #  h_i = Σ_j s_j W[j,i]   ==   (Wᵀ s)_i
```

For probe 1 this is harmless — `W = Wᵀ`. Here `W = S + λA` with `A` antisymmetric, so `Wᵀ = S − λA` and **the dynamics applies −A**. `stdp_sequence_A` chains `Aξ^μ → +ξ^{μ+1}`; the dynamics therefore drives toward `ξ^{μ−1}`. `p_next` measured the overlap with `ξ^{μ+1}` — the direction the field pushes *away* from, which is exactly why it came back slightly negative.

Mirroring it (8 seeds, structured, N=500, P=50):

| λ | `p_next` (as recorded) | `p_prev` (the direction the dynamics chains) | `p_self` | `frac_fixed` |
|---|---|---|---|---|
| 0.10 | −0.0023 ± 0.0043 | +0.0026 ± 0.0037 | 0.996 | 0.630 |
| 0.15 | −0.0099 ± 0.0090 | +0.0080 ± 0.0052 | 0.976 | 0.458 |
| 0.20 | −0.0165 ± 0.0099 | **+0.0203 ± 0.0085** | 0.872 | 0.283 |
| 0.30 | −0.0050 ± 0.0194 | −0.0006 ± 0.0263 | 0.287 | 0.043 |

So a directional sequence drift **does** occur — and it is the one measurement in the probe where temporal order is resolvable (the aligned controls sit at zero). It is also two orders below a transition.

**But P2 was unreachable anyway, and that is the more useful finding.** From numbers already in spec §3, with no run at all:

```
self term   (S ξ^μ)·ξ^μ     / N = 0.9945
chain term  (A ξ^μ)·ξ^{μ+1} / N = 0.7003
→ λ needed for the chain to overcome the self term = 1.42
→ fixed points are entirely gone by λ ≈ 0.3
```

**The design could not have produced a sequence transition at any λ where it had a landscape to measure.** Spec §8 said P2 was "most likely to fail for an implementation reason" (Sompolinsky–Kanter needs a *slow* asymmetric synapse; these are instantaneous) and shipped it anyway. Naming a risk is not the same as computing that the arm is dead — this is `docs/working_rules.md`'s probe-5 lesson (*check the design can actually produce the phenomenon's shape*) recurring in the probe whose own §7 is about that catalogue.

**Correction:** P2 is withdrawn from the verdict list. Nothing was learned about sequence transition. The reachability calculation belongs in the pre-flight so the next spec catches it before the run, and the harness needs a name–content test that drives the **dynamics** and checks which neighbour the state moves toward — the existing `test_structured_field_points_at_the_next_pattern` asserts on `A @ pats[mu]`, never touches `converge`, `is_fixed_point` or `W`, and **would pass with `converge` deleted**. That is verbatim the probe-2 f–I failure, Type C.

## 6. What it does to the registers

*Every V1 bullet here is replaced. Three of the four did not follow from the results even before §1.*

- **O18 — answered, and the answer is "the question was in the wrong units."** Pattern-aligned antisymmetry is **not** intrinsically more corrosive to the energy landscape than isotropic antisymmetry. Per unit of field delivered to the stored patterns it is marginally *less* corrosive, by a factor rising from 1.04 at α = 0.10 to 1.60 at α = 0.02. The 3× that appears at matched `‖A‖_F` is `√(N/P)`. **Move O18 to closed, with the finding stated in field units and the Frobenius statement kept as the unit-conversion note.**
- **T1 — confirmed, not corrected.** T1 reported its table for random `A` and explicitly refused to generalise: *"`A` here is random — the near-null case."* This probe **reproduces T1's number on an independent code path** (0.80 energy violations at λ = 0.5 vs T1's 0.8) and adds the row T1 named as untested. V1's *"T1's correction softened it slightly too far"* is withdrawn: nothing here shows T1 said too little about what T1 measured, and it sits badly beside this probe's own conclusion that the operator class, not the structure, is what differs. **T1 now carries both rows and the units caveat.**
- **Doc 27 §3's biological-support rating for energy — no downgrade.** Proposed wording: *strong for attractor dynamics; for the energy **function**, robust to antisymmetry in proportion to the field the antisymmetric component delivers to the stored patterns, with no additional penalty for the pattern-aligned, low-rank operators biology's learning rules generate. In `‖A‖_F` units the tolerance falls as roughly `α^{−0.7}` with load and is `√(N/P)` times smaller for aligned than for isotropic operators — a bookkeeping fact about the norm, not a fragility of the formalism.* V1's *"tolerance shrinking as α^0.22"* was wrong in both form and sign — `α^{+0.22}` is the **ratio**, and it *increases* with α; the tolerance itself is `0.038 · α^{−0.68}`.
- **Probe 1's scope limit is NOT quantified by this probe — withdraw the claim.** V1 wrote *"probe 1's results carry to λ ≲ 0.1 at α = 0.1."* Probe 1 runs at **α ∈ {0.12, 0.14, 0.16, 0.18}** with its headline config at **0.16** (`sweep.py`, `config.yaml`); α = 0.1 appears nowhere in it. O18's own inclusion rule stops at α = 0.12 because the λ=0 baseline `frac_fixed` collapses to 0.139 at α = 0.15 and 0.006 at α = 0.20 — i.e. **O18 has no substrate at probe 1's operating point** and cannot bound it by extrapolating a law fitted over α ≤ 0.10. Worse, probe 1's *post-unlearning* networks have `frac_fixed = 1.000` in every cell, so they are not the Hebbian substrate O18 measured at all. Honest statement: **O18 bounds the Hebbian substrate at α ≤ 0.12; probe 1's scope limit remains unquantified, and quantifying it means re-running the load arm on post-unlearning networks — where the baseline is 1.000 and the measurement is actually possible.**
- **Prior art — the library gap is real** and unchanged. We hold nothing on asymmetric or sequence-storing attractor networks. One correction to V1's framing: PRE 52:5261 counts the **total** number of fixed points (overwhelmingly spurious), while T1 and O18 measure whether **stored patterns** remain fixed points. Its retrieval half is what T1 hand-re-derived; its fixed-point-count half is a different quantity we did not measure. Three acquisitions listed in spec §9.

## 7. How V1's conclusions failed

Five passes on this probe — the author's own self-review, then three adversarial subagents under the standing gate, then independent re-verification of each load-bearing finding before it was written here. The self-review found real things (§7.1). It did not find any of the four below, and the reason is instructive: **it audited the numbers and never audited the units.**

**7.1 What the self-review did catch, and it was right to.** (i) The α ≥ 0.15 loads were measuring nothing — baseline `frac_fixed` 0.139 and 0.006, above capacity — and were in the first fit; excluding them moved the exponent 0.235 → 0.216. (ii) The λ grid truncated the random knee at low load, returning NaN at α = 0.03; extended to λ ≤ 6.0. (iii) The knee threshold was an unexamined free parameter; re-fitting at 25–90 % gives 0.252 → 0.197, so P1's Frobenius exponent is not a threshold artefact. All three reproduce.

**7.2 `λ` was named "how much antisymmetry" and its content was "how much field, times a subspace factor."** Every headline followed from comparing operators at matched `‖A‖_F` and reading the result as a property of the landscape. The check that would have caught it — *does the swept parameter mean the same thing in every condition?* — is not in `docs/25`, and is the genuinely new thing this probe produced (see §7.5). Note the shape: this is **name ≠ content one level up**, at the swept parameter rather than at a function or a recorded field.

**7.3 A control whose null was an algebraic identity.** §3. The repo has written down *a control that cannot produce the outcome is not a control* twice, and violated it here a third time, in a form the mechanised `control_potency` check cannot see — because `preflight_check.py` calls it with `apply_control=lambda c: c` (the identity) evaluated at λ = 4.0, where every condition is floored at `frac_fixed = 0`. All four conditions return byte-identical "PASS", **including `structured`, which is not a control at all**. The check as invoked tests that λ matters, which `depends_on(λ)` already tests. It should plant a difference *between conditions* and confirm the n=5 design resolves it, inside the dynamic range.

**7.4 Prose numbers with no recorded counterpart — the rule broken again, in three places.**
  - *"the fitted law gives 0.315 against √0.1 = 0.316 … agreeing to three decimal places."* 0.31501 vs 0.31623 agree to **two**. And the single-load arm's actually-observed ratio was **0.307** — a 3 % gap described as a 0.1 % one.
  - *"disagreeing by 60 % at α = 0.02."* The prediction is 36.5 % below the truth, or the truth 57.4 % above the prediction. Neither is 60 %.
  - *"The first load fit spanned α ∈ [0.03, 0.20] and returned an exponent of 0.239."* α = 0.03 returns NaN for `random` in all five seeds — §6.2 of V1 says so itself, contradicting §6.1 about the same fit — and the reconstructable fit over [0.05, 0.20] gives **0.2347**. No recorded artifact produces 0.239.
  - *"α = 0.10 because T1 was, and T1 was because probe 1 was."* T1 → α = 0.1 is correct. Probe 1 → **false**: α ∈ {0.12…0.18}. A design-provenance claim written from memory rather than read out of probe 1's files — the same species as the error the probe-5 review caught, one level up from a number.

**7.5 §4's "fifth failure type" is demoted to a sharpening of Type A — and its content is rewritten.** V1 proposed *confirmation at a crossing point* as a fifth entry in `docs/25`. It does not clear the register's own inclusion bar:
  - its detector is Type A's detector (`scan_table` flagging `alpha: constant at 0.1`);
  - it is internal-validity, the same axis as A–D, unlike §5b's external-validity gate, which is the precedent V1 invoked;
  - every other type is anchored to a committed error with a withdrawn claim; this one was populated by a counterfactual ("it *would* have been called confirmed"), and `docs/25` §6 explicitly forbids adding structure on that basis;
  - and the crossing itself is not robust — re-deriving it at each knee threshold in V1's own sensitivity table puts it anywhere from α = 0.082 to α = 0.135. It cannot be a disarmed free parameter in §6 and load-bearing in §4.
  - Most of all, the story was wrong: `√α` was **not** a functional form that happened to cross the truth at one point. It was the correct exponent for the unit conversion, sitting inside a product V1 never factored.

  **What goes into `docs/25` instead — two lines, both Type A sharpenings, both paid for:**
  > *When `scan_table` flags the parameter your prediction is a function of, "constant by design" is not an exemption — it means the prediction is untested.* (Probe 5 dismissed the analogous flag; taking it literally here is the entire reason the load arm exists.)
  >
  > *Before comparing conditions at "matched X", check that X means the same thing in each — that the swept parameter is the quantity the metric responds to, not a proxy that varies with the condition.* (Type A, at the experiment rather than the instrument.)

**7.6 Smaller, recorded rather than folded in silently.**
  - The coarse λ grid inflates the *random* knee — the denominator of every ratio — by ~5 % (0.544 ± 0.071 coarse vs 0.517 ± 0.043 fine). One-sided, and it lands only on the denominator. V1 checked knee-*threshold* sensitivity but not knee-*grid* sensitivity.
  - *"resolvable from λ = 0.05 upward"* is false upward: λ = 1.0 and 1.5 read `+0.000 ± 0.000` with both arms dead. Say "λ = 0.05 to 0.5".
  - The control table is headed "at every λ" and shows 5 of 13 rows. None of the omitted rows is resolvable, so nothing changes — but `docs/25` §4 lists exactly this as a scan target, and probe 3's fatal finding was this pattern.
  - *"energy violations appearing by λ ≈ 0.02"* is **two violating updates in one of five seeds**, compared against T1's 0.5, which is merely T1's first sampled grid point above 0.2 under a different protocol (4000 updates at N=200 vs 1000 at N=500). Within this probe's own protocol the isotropic first-nonzero is λ = 0.15. Spec §5 files this metric as "reported, not interpreted alone"; V1's register bullet interpreted it alone, at the floor. Withdrawn. And *"agrees independently"* is wrong in principle: it is computed on the same `S + λA` under the same Frobenius convention, and under field matching its 74× separation collapses to ~1.5×.
  - T1's *"perfect retrieval to 0.5"* is 0.998 in T1's table. Write ≥ 0.998.
  - `frac_nonconverged` and `basin_radius` are pre-registered in spec §5 and appear in **no** table and **no** sentence. `frac_nonconverged` as recorded is also an artefact of `max_sweeps=25` at λ = 0.15–0.2 (0.008 → 0.000 and 0.080 → 0.000 at 200 sweeps). Report them at `max_sweeps ≥ 200` or strike them from the spec with the reason — do not leave a pre-registered metric silently unreported.
  - `observes` is pre-registered in spec §7 for the energy-violation counter and is never called. The "7/7" pre-flight is 1 × `dynamic_range` + 3 × vacuous `control_potency` + 2 × `empirical_null` + 1 × `depends_on`, **all seven on `frac_fixed`**. `p_next` — the metric that turned out inverted — was never range-checked at all.
  - **Reproducibility gap:** `results/o18_antisymmetry/knee_law.json` carries the headline law and **nothing in the repo produces it**. `run_load.py` still holds the pre-§6 configuration (5 loads, λ ceiling 1.6) and writes `load_sweep.jsonl` — the superseded, withdrawn run — so re-running the committed script *overwrites the results with the wrong numbers*. No artifact carries a config hash or git SHA, which `docs/working_rules.md` requires. The knees themselves were reproduced independently in review, so this is provenance, not correctness — but the repo's headline result currently violates its own `(config, seed) → artifacts` rule.
  - **Mutation testing** (review, verified): 7 of 9 mutants survive all 13 harness tests, including `p_next := 0.0`, `frac_nonconverged := 0.0`, `converge`'s `max_sweeps` forced to 1, and — the serious one — `rank_matched_random_A` returning `stdp_sequence_A(patterns)` **itself**. The probe's whole finding is "the rank-matched control matches structured", and a control that literally *is* the structured operator passes every test. `shuffled_order` has a differentiation test; this one does not. The "13 ground-truth tests" badge overstates what the suite constrains.
  - **A latent fragility, not a live failure:** `converge` computes the field per neuron (`S @ W[:, i]`) while `is_fixed_point` uses one matmul (`S @ W`). The two reduction orders can disagree on exactly-zero fields. Under the repo's environment (numpy 2.2.6) `test_converge_returns_actual_fixed_points` passes and no zero fields occur; under numpy 2.4.4 the same test fails on 4/50 states. Nothing here depends on it — at λ = 0 an exact integer recomputation gives the same `frac_fixed` to four decimals — but the tie-breaking should be shared rather than duplicated.

## 8. Novelty: none claimed

Mode was **intuition-build**, declared in spec §0 before the build. Retrieved before building: Sompolinsky & Kanter 1986 (PRL 57:2861) for the sequence construction; [Phys. Rev. E 52:5261 (1995)](https://doi.org/10.1103/PhysRevE.52.5261), Singh, Zhang & Dasgupta, *Fixed points in a Hopfield model with random asymmetric interactions*, for random-asymmetric fixed points; and [arXiv:2605.24611](https://arxiv.org/abs/2605.24611), *Beyond Fixed Points: Superpolynomial Capacity of Asymmetric Hopfield Networks*, for current SOTA. **Both citations were re-verified against the sources in review; neither is fabricated or misattributed.** Two scope notes: the 1995 result counts total fixed points, not stored-pattern stability (§6), and the 2026 paper's construction is *synchronous* with limit-cycle attractors while this probe is asynchronous — a SOTA pointer, not a comparison class.

## 9. What the review did not settle

- **N-dependence.** Everything is at N = 500. Whether `α^{−0.26}` (the field-matched residual) is a finite-size number is untested; the crosstalk scale `√(P/N)` enters both operators, and a second N would separate α-dependence from N-dependence. **This is the largest open threat to §2's decomposition.**
- **Does rank matter at all?** No condition varied rank at fixed delivered field. Rank-`r` truncations of an aligned operator at constant `RMS(Aξ)`, r = 5…50, would decide it — cheap, and the obvious next run.
- **The `p_prev` drift** is established at λ = 0.2, α = 0.1 (+0.020 ± 0.009, and zero in every aligned control). Not mapped over λ or α, and not tested against the delayed-synapse variant spec §8 names.
- **Probe 1's actual operating point** (α = 0.12–0.18, and post-unlearning networks with `frac_fixed = 1.000`) is unmeasured — see §6.
- A formal variance comparison between conditions under the shared-pattern pairing; and whether `frac_fixed`'s 0.78 baseline at α = 0.1 matches the analytic single-pattern stability prediction.
