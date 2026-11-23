# Working rules for this repo

*Auto-loaded each session. Keep it short — it is the part that actually gets read. Orientation: **`START_HERE.md`** → `NEXT_STEPS.md`. Full protocol and the evidence behind these rules: **`docs/25_probe_protocol.md`**.*

## The one thing that has gone wrong repeatedly

Across M1's three probes, **nine** claims had to be withdrawn or corrected by adversarial review. Not one was caught by the test suite — 133 tests passed the whole time. Every one was an instrument, control, or label problem, and they fall into four repeating types:

| type | example |
|---|---|
| **Instrument can't move** | `spurious_total` saturated at 1.000; separability pinned by its own normalisation; an autocorrelation statistic that mapped [−1, 0.368] to the integer 1 |
| **Control can't fail** | weight decay is provably vacuous under `sign(Ws)` scale invariance; a "matched" shuffled control delivering 7.5× less plasticity; a CV tolerance larger than the CVs compared |
| **Name ≠ content** | an f–I test that never called the integrator (would have passed with `run()` deleted); `mean_input_mv` that omitted the recurrent term |
| **Null isn't null** | greedy pairing's winner's curse put the pure-noise exponent at 1.033 ± 0.002, 16σ off the analytic 1.0 |

## Before interpreting any number — four checks

**These are mechanised in `analysis/preflight.py`. Use it — a checklist that is written and not run is worth nothing.**

```python
from analysis import preflight as pf
pf.substrate_predictability(X)                            # 0: is there ANYTHING here?
pf.dynamic_range(...)   pf.quantisation(...)              # A: instrument can't move
pf.control_potency(...) pf.tolerance_is_discriminative(...)  # B: control can't fail
pf.observes(...)        pf.depends_on(...)                # C: name != content
pf.empirical_null(...)                                    # D: null isn't null
pf.scan_table(rows)                                       # read your own output
print(pf.report(checks))
```

`analysis/tests/test_preflight.py` re-finds all nine M1 failures as its acceptance test.


0. **Measure the substrate before measuring instruments on it.** Every other check validates an *instrument*; this one asks whether the thing being measured exists. Probe 4 was fully specified — five encoders, matched code size, three falsifiers — against data where predicting the next bin from the *entire observed state* scores worse than predicting from nothing (R² = 0.00001; only 4 % of the network was recorded). The check takes **1 second** and would have stopped it before the spec was written.
1. **Dynamic range.** Compute the metric under the *null* and *saturated* conditions **at the experiment's real N, T and parameters**. If the observed value sits within noise of either endpoint, it is dead — do not interpret it. Validating a metric only in the positive direction (does PR recover *d*?) is what let three saturated metrics through.
2. **Control potency.** Show each control *can* produce the effect — plant the effect and check the control reproduces it, or at minimum that it moves the outcome variable at all. **A control that cannot produce the outcome is not a control.**
3. **Empirical null.** Run the *whole pipeline* on no-effect synthetic data at real N/T, several seeds, and report mean ± SD. Compare excesses to that, never to the analytic ideal.
4. **Name–content.** For any function or recorded field whose name asserts a meaning, one test that fails if it computed something else.

## Before writing prose — read your own output adversarially

Scan every results table for columns that are **constant, exactly 0, exactly 1, NaN, or suspiciously flat**. In probe 3 both fatal giveaways were in the printed output and unexamined: a flat `2.14, 2.14, 2.12…` column, and `alpha_autocorr = 0.00` in four of five conditions with `r² = NaN`. Also: **never quote a number in prose that differs from what the code recorded** (probe 2 quoted 21.4 mV while the JSON said 40.0).

## Standing gates

- **Adversarial review before any register entry moves.** Brief a subagent to find the work wrong, give it the claims explicitly, and require it to run code. Three passes cost ~45 min of compute and caught more than the entire suite. Verify its findings yourself before acting.
- **Ground-truth validation before a metric touches real data** (`analysis/tests/`), in *both* directions — positive recovery and degenerate/null behaviour.
- **Every probe states a prediction, a falsifier, and a control that could produce the same result for a boring reason** before it runs.
- **If a probe motivates itself by an empirical phenomenon, cite the row in `docs/directions/CONSTRAINTS.md`** — and check the design can actually produce the phenomenon's *shape*, not just beat a baseline. Probe 5's recency arm was motivated by the serial-position curve and could not produce a U (no primacy mechanism existed in the design): a Type A failure hiding in a spec. If there is no matching row, the motivation is folk psychology.
- **Exploratory ≠ pre-registered.** If the operating point or the primary metric was chosen after seeing data, the run is exploratory. Say so. Documenting a forking path does not convert it into a pre-registered one.

### External validity — the prior-art gate (added after a session mis-graded probe 5 as novel from recall)
The gates above secure *internal* validity ("is it real?"). These secure *external* validity ("is it new / already known?"), which had no gate and cost us: a VSA-cleanup-capacity probe was called "plausibly novel" while the Frady 2018b capacity theory and Hopfield-unlearning results sat in our own library, unread.

- **Novelty is provisional until retrieved — never claim "novel" from the model's recall.** Recall is unreliable for *is-this-known* (the same failure as citing from memory) but *useful for pointing at which fields to search*. Until a landscape + targeted search (library **and** web, in ≥2 vocabularies — the math, ML, and neuro phrasings, because cross-field prior art hides under other names) is done, a result is "unknown-provisional," not new.
- **Before mining a D-thread, load its landscape** — `docs/directions/*_D<n>_landscape.md` (state of the field, retrieved). Direction-setting reads the landscape *first*; it is not a probe-time afterthought. `21_D2_landscape.md` is the prototype. If a thread has no landscape doc, writing one (from the surveys we hold + a web scan) is the first task, not the probe.
- **Name the probe's mode before building it.** *Intuition-build* — reproduction is fine, novelty is off the table. *Contribution* — prior art is mandatory and the landscape must be current. Probe 5 was a fine intuition-build mis-labelled as a contribution candidate; the mode names the difference and stops the mis-grade.

## Recording results

- Corrections go in the write-up **as corrections**, not silent edits — how a claim failed is usually more transferable than the claim.
- **But the *lesson* is timeless; the *session* is not.** Keep "X was claimed, it's wrong because Y" — drop "(Session 19)". Living/reference docs (primers, currency docs, the registers, START_HERE) state the current understanding and the *why* of any correction, with **no session watermarks**: the incremental "what happened when" is exactly what `journal.md` is for, and git recovers the date. **Frozen docs are the exception** — they *are* session records (marked Frozen in `docs/README.md`), so they keep their identity. A section titled `(added Session N)` is cruft; a correction that reads `⚠ corrects the earlier over-claim that…` is not.
- Refutation is a result. `ISOMORPHISMS.md` and `OPPORTUNITIES.md` are designed to hold demotions.
- Runs are `(config, seed) → artifacts` with config hash + git SHA, resumable, and seed-pure so a run split across calls is bit-identical.

## Using the library — search first, do not cite from memory

**The corpus is 218 records. You have not read them; the model's training has. Those are different, and conflating them is a documented failure here.** Doc 27 §6a "rediscovered" five sources we already held — sparse coding, reservoir, polychronization — by the model *recalling* they existed, which is retrieval from priors wearing a library costume. They had sat mis-tagged and unsearched for sessions. The library exists to cite what the model would **not** recall, and to quote what a paper **says** rather than what it is remembered to say.

**So before any doc makes an empirical or citation claim, reach the corpus through tooling — never from recall:**

```
python3 library/scripts/lib.py find "sparse code"        # keyword/tag; zero setup, always works
python3 library/scripts/lib.py find "binding" --thread D2 # filter by thread/section
python3 library/scripts/lib.py semantic "criticality vs reduction"  # fuzzy; needs fastembed (vendored)
python3 library/scripts/lib.py read 119                   # grep-ready text path + head
```

- **`find` is the floor** — pure stdlib over stored keywords/tldr/synopsis, so it works in a cold session with nothing installed. Use it every time.
- **Then `grep library/text/<slug>.txt` and READ before citing.** A `#NN` in prose must point at something the text actually says — the tldr is a signpost, not a source. (Same discipline as "never quote a number that differs from what the code recorded.")
- **After writing a doc, run `python3 library/scripts/check_citations.py`** — it flags sources cited beside a thread the index doesn't carry (retag them, citing the specific claim), and reports the ~22 % of the library never cited anywhere. Retagging is per-doc first, that script as the audit. See `library/PIPELINE.md` Procedure D.

## Environment notes

- Bash calls time out at **45 s** and each runs in a fresh PID namespace — background jobs do not survive. Chunk long work into resumable foreground steps with artifact-based progress.
- numpy-first for simulation; cross-check against an independent implementation (Brian2) rather than trusting a hand-rolled integrator alone.
- **Prefer an analytic or theoretical cross-check over more unit tests.** Mean-field agreement across probe 2's whole phase diagram constrained the integrator far more than any single-neuron test, and it was the check I initially skipped.
