# Probe protocol — how to not fool ourselves

*Written after M1, where **nine** claims across three probes had to be withdrawn or corrected by adversarial review. The short operational version is in `working_rules.md` (auto-loaded each session); this file carries the evidence, so the rules can be argued with rather than obeyed. Living.*

---

## 1. Why this exists

M1 produced a harness and three probes. It also produced this record:

| probe | headline | what failed |
|---|---|---|
| 1 — sleep as unlearning | **reversed in both directions** | refutation drawn from a saturated metric; all three controls incapable of failing |
| 2 — Brunel testbed | **stood**; 3 of 5 validation bullets defective | a test that never called the integrator; a function that didn't compute its name; a tolerance larger than the values compared; a non-reproducible sweep |
| 3 — phenomenological RG | **stood and strengthened** | O12 result was a normalisation artifact; a "negative" was a broken statistic hiding a real positive |

Two things to take from that table. First, **the primary results got more robust over the three probes** — probe 3's headline survived a determined attack including a confound I hadn't thought to check. The trend is good. Second, **the apparatus around the results did not improve at all**, because I never diagnosed why it kept failing. This document is that diagnosis.

**None of the nine was caught by tests.** 133 tests passed throughout. This is not an argument for more tests; it is an argument that these failures are of a kind tests do not see.

---

## 2. The four failure types, with the real cases

### Type A — the instrument cannot move

The metric has no dynamic range in the region where the answer lives, so it reports the same value whatever the truth is.

- **`spurious_total` saturated at 1.000** (probe 1). Defined as `1 − genuine`, and genuine basin volume is exponentially small in N, so it read 0.000 in every condition at every checkpoint. With 300 samples and zero hits, a true change from 10⁻⁶ to 10⁻² is invisible. I drew a *strong negative conclusion* from it — and I had already documented the flooring on the genuine side one paragraph earlier.
- **Separability pinned by its own normalisation** (probe 3). Features were scaled by the across-regime standard deviation, then between-regime distances measured on that same scale — which fixes the numerator by construction (2.61 ± 0.05 on pure random features).
- **An autocorrelation statistic mapping [−1, 0.368] to the integer 1** (probe 3). "Lags until autocorrelation falls below 1/e" returns 1 for any anti-correlated series. It destroyed a real result: lag-1 autocorrelation actually runs −0.012 → **+0.612** across cluster sizes in one regime.

**Why it kept happening:** every metric was validated in the *positive* direction only. Participation ratio recovers *d*; the power-law fitter recovers α; RG recovers exponents 1 and 2. None was checked for **range at the experiment's actual N, T and parameters**.

**Two sharpenings, both paid for by O18 — read "instrument" as "experiment", not only "metric":**

- **When `scan_table` flags the parameter your prediction is a function of, "constant by design" is not an exemption — it means the prediction is untested.** O18 pre-registered a law in α and ran a single-load arm; `scan_table` duly returned `alpha: constant at 0.1`, which reads like bookkeeping (it *is* constant by design in a one-load arm) and which probe 5 had dismissed on exactly that reasoning. Taking it literally and building the load arm is the only reason the law was tested at all. A point prediction confirmed at one operating point is not evidence for the function it was derived from.
- **Before comparing conditions at "matched X", check that X means the same thing in each.** The swept parameter must be the quantity the metric responds to, not a proxy that varies with the condition. O18 swept `λ = ‖A‖_F/‖S‖_F` across operators with ranks 48 to 500, concluded that pattern-aligned antisymmetry is 3× more corrosive to the energy landscape, and was measuring `√(N/P)`: a subspace-confined operator delivers that much more *field* at matched Frobenius norm, and field is what `sign(Ws)` responds to. Matched on field, all four conditions broke at the same place. **Nothing in the numbers was wrong; the units were the finding.** Neither the test suite, `scan_table`, nor the author's own audit could see it, because every one of them checks numbers against numbers.

### Type B — the control cannot fail

The comparison is structurally incapable of producing the effect, so passing it is guaranteed and carries no evidence.

- **Weight decay under scale invariance** (probe 1). Hopfield dynamics depend on `sign(Ws)`, invariant under `W → cW`. Uniform decay cannot move a single fixed point however it is tuned. It "ruled out weight normalisation" *necessarily*, not evidentially.
- **A "matched" shuffled control delivering 7.5× less plasticity** (probe 1). Matched on *update count*; random rank-1 updates cancel incoherently, so ‖ΔW‖_F was 0.63 against sleep's 4.41.
- **A CV tolerance of 0.15 against CVs of 0.09–0.14** (probe 2). A simulator emitting perfectly clock-like spike trains would have passed.

- **A shuffled-order control that was an exact relabelling** (O18). `stdp_sequence_A(patterns, order=σ)` is *identically* `stdp_sequence_A(patterns[σ])` (`max|Δ| = 0.0`), and `S = ξᵀξ/N` is permutation-invariant, so with i.i.d. patterns the "order-destroyed" condition is a re-draw of the structured condition from the same distribution. Every permutation-invariant statistic has the same law in both **by construction**; the expected difference was exactly zero before the run. **The general form: with exchangeable items there is no order to destroy** — "order" is a labelling of interchangeable objects, and a potent order control needs either sequentially correlated items or a statistic that is not permutation-invariant.

**The tell:** I identified the weight-decay vacuity myself, wrote it up as a lesson — and then failed to apply it to the shuffled control in the same probe. Three probes later O18 shipped a shuffled control whose null was an algebraic identity, and `pf.control_potency` passed it — because it was invoked with the identity function as the control, evaluated at a λ where every condition was floored. **A mechanised check invoked wrongly is worse than no check: it launders the omission as a PASS.**

### Type C — name ≠ content

The artifact is called something it isn't, so downstream reasoning is about a fiction.

- **An f–I test that never called the integrator** (probe 2). It built a network, silenced it, then hand-rolled the integration loop *inside the test body*. It would have passed with `BrunelNet.run` deleted. The physics was right; the test was theatre — and it read as rigour.
- **`mean_input_mv` omitting the recurrent term** (probe 2). Identically `ν_ratio · θ`. The records said 40.0 mV where the true value was 21.1, while the write-up quoted 21.4 computed by hand off-record *and invited readers to check the raw numbers*.

### Type D — the null isn't null

The baseline carries an uncorrected bias, so "excess over baseline" is measured against the wrong zero.

- **Greedy pairing's winner's curse** (probe 3). Max-correlation pairing selects the largest sample correlations, which among independent units average r ≈ 0.07. The pure-noise exponent is **1.033 ± 0.002** — 16σ from the analytic 1.0. Random pairing recovers 1.000.

---

## 3. The pre-flight — four checks before interpreting any number

> **Mechanised in `analysis/preflight.py`.** Each check below names the function that implements it. `analysis/tests/test_preflight.py` is the module's acceptance test: it reproduces all nine failures below and asserts preflight catches each one — so if the module ever stops working, that suite says so. Prose alone would have decayed within a session.

Each check targets one type. Together they mechanically catch seven of the nine.

**0. Substrate check — added after probe 4** *(`pf.substrate_predictability`)*
Before building any instrument, ask whether the substrate has the structure the probe intends to measure. Fit the best linear predictor of the next state from history (**blocked** CV — interleaved folds leak across adjacent bins of any autocorrelated series) and compare against rate-matched surrogates.
*Would have caught:* probe 4, whose primary regime scores R² = 0.00001 — indistinguishable from independent noise, because only 512 of 12 500 neurons (4 %) were recorded and each receives fresh external Poisson drive every step. The probe's own falsifier would have fired by construction. **Cost: 1 second.** Every other check here validates an instrument; this is the only one that asks whether the thing being measured exists.

**1. Dynamic-range check** *(Type A — `pf.dynamic_range`, `pf.quantisation`)*
Compute the metric under the null condition and the saturated condition **at the experiment's real N, T and parameters**. Report both endpoints beside the observed value. If the observation sits within noise of either, the metric is dead — do not interpret it, and say so.
*Would have caught:* `spurious_total` (endpoints 1.000/1.000), separability (numerator pinned), the autocorrelation statistic (constant 1).

**2. Control-potency check** *(Type B — `pf.control_potency`, `pf.tolerance_is_discriminative`)*
Before crediting a control, show it *can* produce the outcome: plant the effect and verify the control reproduces it, or at minimum verify it moves the outcome variable at all. Where a control's inertness is provable — as with scale invariance — say so and stop treating it as evidence.
*Would have caught:* weight decay, the under-powered shuffle, the CV tolerance.

**3. Empirical-null calibration** *(Type D — `pf.empirical_null`)*
Run the **entire pipeline** — not just the metric — on no-effect synthetic data at the real N and T, several seeds. Report mean ± SD. Interpret excesses against that, never against the analytic ideal.
*Would have caught:* the greedy-pairing bias.

**4. Name–content test** *(Type C — `pf.observes`, `pf.depends_on`)*
For any function or recorded field whose name asserts a meaning, write one test that fails if it computed something else. For a test, the question is sharper: **could this test pass if the thing it names were deleted?** If yes, it is not testing that thing.
*Would have caught:* the f–I test, `mean_input_mv`.

---

## 4. The cheapest check of all — read your own output adversarially *(`pf.scan_table`)*

Before writing any prose, scan every results table for:

- columns that are **constant or near-constant** across conditions;
- values of **exactly 0, exactly 1, or NaN** repeated across conditions;
- **suspiciously round agreement** between a real condition and its control;
- rows **omitted** from the summary table that are present in the data.

In probe 3 both fatal findings were visible in the printed output and went unexamined: a flat `2.143, 2.137, 2.124, 2.106…` column, and `alpha_autocorr = 0.00` in four of five conditions with `r² = NaN`. The summary table also silently dropped the one row where the flat column jumps.

And one absolute: **never quote a number in prose that differs from what the code recorded.** If they differ, one of them is wrong and you don't yet know which.

---

## 5. Statistical hygiene that came up repeatedly

- **Fit quality is not evidence.** r² ≈ 0.99 appeared in all ten of probe 3's conditions *including the surrogates*. A power law fitted to a handful of points over two decades is nearly unavoidable. Report the comparison against a surrogate, not the fit.
- **Independence assumptions in SEs.** Probe 1's retrieval trials were clustered within patterns, giving a design effect of ~1.6 that the naive `√(0.25/n)` ignored.
- **Thresholds are conventions.** Report the sensitivity table, not the single value. Probe 3's "16 of 56" ranged 8–21 under ±30 % threshold moves; the defensible claim was the topology, not the count.
- **Selection on the dependent variable.** Probe 1's α = 0.16 was chosen *because it gave the outcome variable headroom*. That is legitimate exploration and illegitimate pre-registration.
- **Convergence in the nuisance parameters.** Probe 3's exponents were not converged in N; the ordering was. Say which.

---

## 5b. A fifth failure, on a different axis — the result is real but already known

The four types above are all *internal*: the apparatus can only say one thing. A later session found a fifth that none of them catch, because it is *external*: **the result is real, and already in the literature, and I called it novel from recall.** A VSA-cleanup-capacity probe (probe 5) was graded "plausibly novel" while the Frady 2018b capacity theory and the Hopfield-unlearning results sat in our own library unread — the same recall-over-retrieval failure as citing from memory, one level up.

The fix is mostly **upstream of this document** and lives in the auto-loaded rules (`working_rules.md` "external validity — the prior-art gate") and the routing table (`../START_HERE.md`): load the thread's **landscape** (`docs/directions/*_D<n>_landscape.md`) before choosing a direction, and treat **novelty as provisional until retrieved, never claimed from recall**. Direction-setting is where prior art does its real work; a probe-time check is too late to redirect it.

What belongs *here*, at probe time, is only the last two narrow steps:
- **Declare the probe's mode before building.** *Intuition-build* (reproduction is fine, novelty off the table) or *contribution* (prior art mandatory, landscape current). Probe 5 was a fine intuition-build mis-labelled as a contribution candidate; naming the mode is what stops the mis-grade — and, unlike the internal gates, costs one sentence.
- **If contribution: one targeted prior-art confirmation** that the *specific* result isn't already published (library + web, in the math / ML / neuro vocabularies — cross-field prior art hides under other names). This is the narrowest layer and only works because the landscape read already happened upstream.

This is not the "more ceremony" §6 warns against: it is a *different axis* (is-it-known vs is-it-real), it was genuinely absent, and it cost real work — the exact justification the four internal checks had.

## 6. Gates, and what not to add

**Keep:**
- **Adversarial review before any register entry moves.** Brief a subagent to find the work wrong, hand it the claims explicitly, require it to run code, and verify its findings before acting. Three passes cost ~45 minutes of compute and caught more than 133 tests.
- **Prefer an analytic or theoretical cross-check to more unit tests.** Mean-field agreement across probe 2's phase diagram (rates to 3–8 %) constrained the integrator far more than every single-neuron test combined — and was the check I skipped in v1.
- **Ground-truth validation in both directions** before a metric touches real data.

**Do not add:** more pre-registration ceremony, a methods manual, or additional review layers. The failures were not caused by insufficient structure — this repo had a spec, registers, and a stated pre-registration discipline throughout, and produced all nine anyway. The missing step was interrogating numbers, which is four checks and a careful read, not a process.

---

## 7. The honest summary

The recurring error is not carelessness about *results*. It is that I built **measuring instruments and controls without measuring them**, then read their output for confirmation rather than for range. Every one of the nine is a case where the apparatus could only say one thing, and I heard it say something.

The corresponding habit is small: *before believing a number, find out what it would say if there were nothing there.*
