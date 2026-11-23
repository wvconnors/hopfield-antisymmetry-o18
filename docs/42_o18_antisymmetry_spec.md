# O18 — how much *structured* antisymmetry can the energy landscape survive?

*Pre-registered spec. Written before any probe code. Motivated by **C1** (STDP is temporally antisymmetric) via tension **T1**; opportunity **O18**; threads **D2** (energy as a currency) and **D6**. Machinery: `experiments/p1_unlearning/hopfield.py`, unmodified.*

---

## 0. Mode declaration — **intuition-build, not contribution**

`docs/working_rules.md` requires naming the mode before building. **This is an intuition-build with an internal-scope question. Novelty is off the table**, and the prior-art gate below says why in detail. What the probe is *for* is two internal jobs the literature does not do for us:

1. **Bound probe 1's scope.** Probe 1 lives entirely inside `W = Wᵀ` (`assert_valid` checks symmetry every step). Its results are results about the formalism. How far that formalism reaches into biologically-plausible weight matrices has never been measured *here*.
2. **Grade doc 27 §3's biological-support rating for the energy currency**, which T1 argues conflates "attractor dynamics are well supported" with "the symmetric-weight energy *function* is well supported."

Both are questions about **our** claims. Neither is a claim about the field.

## 1. Prior-art gate — retrieved, and it changes the design

**Library: verified gap.** `lib.py find` on `asymmetric`, `temporal association`, `sequence`, `STDP`, `limit cycle` returns **nothing** on asymmetric or sequence-storing attractor networks. The only STDP source is #223 (Bi & Poo). Given that the energy landscape is a scored D2 currency and its two contract slots (`compare`, `evolve`) both depend on symmetry, **this is a real collection gap, not a tagging one** — recorded in §7.

**Web scan (two vocabularies — statmech and ML):**

- **Sompolinsky & Kanter (1986), *Temporal Association in Asymmetric Neural Networks*, PRL 57:2861.** The canonical construction: asymmetric connections built from a stored *sequence* make the network recall time sequences and cycles of patterns. In their model the asymmetric term acts through synapses with a **slow dynamic response** — the delay is load-bearing for clean sequential retrieval, not decoration.
- **Phys. Rev. E 52:5261 (1995), *Fixed points in a Hopfield model with random asymmetric interactions*.** Adding a random antisymmetric part causes an **exponential decrease in the number of fixed points**, but at small relative strength and low load it does **not** substantially degrade retrieval quality.
- **The symmetry-parameter transition literature.** Symmetric nets have fixed points or period-2 cycles; as asymmetry grows there is a **critical symmetry parameter** past which typical limit-cycle length grows exponentially and the dynamics become chaotic.
- **Current SOTA:** arXiv 2605.24611 (2026), *Beyond Fixed Points: Superpolynomial Capacity of Asymmetric Hopfield Networks* — asymmetry as a capacity *resource*, not only a hazard.

**Two consequences, both design changes rather than footnotes.**

**(a) T1's "partial answer" was a re-derivation, and it should be labelled one.** The random-`A` table in `CONSTRAINTS.md` T1 — descent survives to ‖A‖/‖S‖ ≈ 0.2 with zero violations, retrieval perfect to 0.5 — is a hand-measured instance of the 1995 PRE result. It was correct, and it was presented as if we had found it. That is the same novel-from-recall failure the external-validity gate exists to catch, one level down: not claiming novelty out loud, but not checking either.

**(b) The probe must *predict against* the theory, not measure blind.** This is probe 5's closing lesson applied at the start instead of the end. §3 therefore states a prediction with a **number derived from the structure of the two operators**, not a direction.

## 2. The question, stated so it can fail

T1 left this: *does **structured** antisymmetry break descent at a ratio where **random** antisymmetry does not?* Random `A` is the null case; STDP-generated `A` is correlated with the stored patterns and their order, and probe 1 established that **structured** spectral operations are what move this system (eigenvalue flattening was a 99.4 %-diagonal operation).

## 3. Pre-registered prediction — with a number

Write `W = S + λ·Â`, `S` the Hebbian symmetric part, `Â` antisymmetric and normalised so `‖Â‖_F = ‖S‖_F`, so `λ = ‖A‖/‖S‖` is the swept ratio.

The two operators differ in **rank and alignment**, and that is computable:

- **Structured** `A_seq ∝ Σ_μ (ξ^{μ+1}ξ^{μᵀ} − ξ^{μ}ξ^{μ+1ᵀ})` has rank ≤ 2P and its row space **is** the pattern subspace. Acting on a state near `ξ^μ` it produces a **coherent O(1) field along `ξ^{μ+1}`** — a push toward the next stored pattern.
- **Random** `A_rand` is isotropic over all `N²` directions. Its component along any particular pattern is `O(1/√N)` — **incoherent**.

At **matched Frobenius norm**, the structured operator concentrates its power on ≈ `P` directions where the random one spreads it over ≈ `N`. So its per-direction amplitude *where the symmetric part's signal lives* is larger by ≈ `√(N/P) = 1/√α`.

> **P1 (primary, quantitative).** Structured antisymmetry degrades retrieval and descent at a ratio **smaller than random by a factor ≈ √α**. At `α = 0.1` that is **≈ 0.32**: where the random-`A` knee sits near `λ ≈ 1.0`, the structured knee should sit near **`λ ≈ 0.3`**, and structured `A` should produce **non-zero energy-increasing updates at `λ = 0.2`, where random produces exactly zero**.
>
> **P2 (mechanism, structured only).** Where structured `A` destroys a fixed point, the state should not wander — it should land on or near **`ξ^{μ+1}`**. Sequence transition, not noise. (Secondary: random `A` cannot produce this, so it is a mechanism readout, not a fair head-to-head.)
>
> **P3 (regime).** The `λ` at which fixed points are lost should **fall as load `α` rises** — more patterns, more coherent transition drive, and less margin.

## 4. Falsifiers (stated before the run)

- **F1 — the sharpened question answers "no".** Structured and random knees coincide within the across-seed SD at matched `λ`. Then "structure matters" does **not** transfer from probe 1 to this operator, O18's sharpening was wrong, and T1's random-`A` table is the whole answer. *Honourable exit; cheap; run first.*
- **F2 — the factor is wrong in the informative direction.** The knees separate but not by ≈ `√α` (e.g. structured collapses at `λ < 0.05`, or the ratio does not track `√α` across loads). Then the rank/alignment account in §3 is incomplete and P3's load sweep is what diagnoses it.
- **F3 — no knee anywhere.** If neither condition shows a knee inside `λ ∈ [0, 2]` at the chosen `N`, the instrument is dead and nothing is interpreted (see §5).

## 5. Instruments — and the one that is already known to be pinned

> **⚠ The T1 table's two columns are both dead in the region of interest.** Energy-increasing updates read `0.0, 0.0, 0.0, 0.0` for `λ ≤ 0.2` (floor) and overlap reads `1.000, 1.000, 1.000, 1.000` (ceiling). **Both instruments are pinned exactly where P1 says the interesting difference lives.** This is a Type A failure sitting in a register, and it is why the primary metric below is neither of them.

| metric | why | range check |
|---|---|---|
| **`frac_fixed`** — fraction of stored patterns that remain exact fixed points | **primary.** Graded on [0,1], moves in both conditions, and it is the quantity Sompolinsky–Kanter and the PRE result both speak to | null (λ=0) ≈ 1.0, saturated (λ→∞) ≈ 0; must show intermediate values |
| **`frac_noncoverged`** — states not settled within `max_sweeps` | the fixed-point → limit-cycle transition itself | 0 at λ=0; must rise |
| **`basin_radius`** — probe 1's metric | graded, and makes the two probes commensurable | must not saturate at the real N, P |
| `energy_violations` | secondary **only**, given the floor above | reported, not interpreted alone |
| `p_next` — overlap with `ξ^{μ+1}` after escape | P2's mechanism readout | structured only |

**A convenience worth an explicit test, not an assumption:** for any antisymmetric `A`, `sᵀAs = 0` identically, so `HopfieldNet.energy` returns the **symmetric-part** energy even after `A` is added. That is exactly what "energy-increasing update" should mean here — and it is a name–content claim, so §6 tests it rather than trusting it.

## 6. Controls — each able to produce the outcome

Per the standing rule (*a control that cannot produce the outcome is not a control*), both of these are **expected to have a real chance of firing**.

1. **Rank-matched random `A` (the load-bearing control).** A random antisymmetric operator confined to the *same rank / same pattern subspace* as `A_seq`. If it breaks descent at the same `λ` as the STDP-sequence operator, then the effect is about **low rank and subspace alignment** — not about sequence structure, temporal order, or STDP at all. **This is the boring explanation for P1, and it is a likely one.**
2. **Shuffled sequence order.** Same STDP kernel, same patterns, temporal order randomly permuted. Separates *order* from *pairing structure*. If the shuffle behaves like the ordered version, "STDP applied to a sequence" reduces to "any pattern-aligned antisymmetric coupling."
3. **Norm-matching is the measurement, so it must be exact.** `λ` is defined by Frobenius ratio and asserted per cell — a mismatch here would manufacture P1 outright.

## 7. Pre-flight (`analysis/preflight.py`) — run before interpreting anything

`dynamic_range` on `frac_fixed` at the real `N, P` under λ=0 and λ→∞ · `control_potency` for both controls · `empirical_null`: **re-run the whole pipeline on the random-`A` condition at the same N/P/seeds** and compare structured excesses to *that*, never to the T1 table (different code path) · `depends_on(λ)` · `observes` for the energy-violation counter · `scan_table` on every results table before prose.

## 8. Design

`N = 500`, `P = 50` (`α = 0.1`, matching T1's check), sequence = the P patterns in storage order, `λ ∈ {0, 0.05, 0.1, 0.15, 0.2, 0.3, 0.5, 0.7, 1.0, 1.5}`, **5 seeds**, conditions `{structured, random, rank_matched_random, shuffled_order}`. P3 adds `α ∈ {0.05, 0.10, 0.20}`. Pure numpy; seed-pure and resumable per `(config, seed)`.

**Scope limit, stated now:** this is a **rate/binary** model with instantaneous synapses. Sompolinsky–Kanter's sequential retrieval needs a *slow* asymmetric synapse; without delay we should expect degraded fixed points and possibly cycles, but not clean sequence replay. **P2 is therefore the prediction most likely to fail for an implementation reason rather than a scientific one**, and it is filed as secondary for exactly that reason.

## 9. Library acquisitions this generated

- **Sompolinsky & Kanter (1986)**, PRL 57:2861 — the construction this probe builds.
- **Phys. Rev. E 52:5261 (1995)** — random asymmetric fixed points; the paper T1's table re-derived.
- **arXiv 2605.24611 (2026)** — asymmetry as capacity resource; the current-SOTA check any forward claim must clear.
