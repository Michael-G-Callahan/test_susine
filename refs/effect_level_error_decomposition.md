# Effect-Level Error Decomposition for AUPRC

**Date:** 2026-04-08
**Status:** Conceptual framework established; blending experiment design unresolved

---

## 1) Motivation

AUPRC is our primary performance metric, scored on the **top 8 causal variants by |beta|** out of 23 total causals, with L=10 effects available. We want to decompose AUPRC error into mechanistically distinct sources at the **effect level**, so we can diagnose *why* methods differ and *which intervention* (exploration, annotation, model specification) addresses which source.

---

## 2) Three Sources of Error

### Missing
A graded causal has no effect component that "owns" it. The model's budget of L=10 effects was consumed by other targets (ungraded weak causals, noncausal LD proxies, noise), leaving a top-8 causal with low PIP. This is a **recall** loss.

### Hallucinated
An effect component concentrates mass on noncausal variants. Note: with 23 true causals in the region, an effect locking onto an ungraded weak causal is correct model behavior but unrewarded by the top-8 grading. This is a **precision** loss (inflates PIPs of noncausal variants).

### Blended
An effect covers a graded causal but dilutes alpha mass across LD neighbors. The causal's PIP is suppressed (recall loss) and noncausal neighbors' PIPs are inflated (precision loss). The effect "found the right LD block" but lacks resolution within it.

Each source implies a different intervention: missing -> better exploration/more restarts; hallucinated -> better model specification or budget allocation; blended -> better within-block resolution (annotations, sharper priors).

---

## 3) Per-Effect Diagnostic Axes

### Axis 1: Decisiveness (k_eff)

The effective number of variants in an effect's alpha vector, measured as the exponential of Shannon entropy:

    k_eff_l = exp(-sum_j alpha_lj * log(alpha_lj))

- k_eff ~ 1: one-hot, fully decisive
- k_eff ~ 5-10: blended but possibly still useful
- k_eff > 100: effectively diffuse / unconverged

This is superior to CS size because it measures the actual concentration of posterior mass, not the size of a thresholded set. It is also superior to CS purity because purity measures worst-case pairwise LD within a set, which can be low (e.g., 0.5) even when the effect is doing its job (small CS covering the causal).

### Axis 2: Accuracy Ratio

For each effect, compare its peak alpha on a graded causal to its overall peak:

    accuracy_l = max_{j in graded}(alpha_lj) / max_j(alpha_lj)

- 1.0: the effect's single most confident pick is a graded causal
- Near 0: the effect is confidently targeting something else (ungraded causal, noncausal, or proxy)

### Why These Two Axes, Not CS Size + Coverage

1. **No filtering step.** CS-based metrics require purity filtering (selection bias: low-purity effects are dropped, hiding information about blended-but-partially-correct effects).
2. **Continuous, not binary.** Coverage is 0/1; accuracy_ratio gives the degree of correctness.
3. **Separates concentration from correctness.** CS size conflates blending with LD structure; k_eff directly measures posterior concentration regardless of LD.
4. **LD accounted for.** Mass on LD proxies inflates k_eff, so blending from LD is captured in the decisiveness axis.

### Interpreting Joint Shifts Across Methods

When comparing method A vs method B (e.g., vanilla SuSiE vs SuSiNE), plot the joint distribution of (log k_eff, accuracy_ratio) across all effects x all datasets, colored by method. Shifts in this joint distribution are interpretable:

- **k_eff drops, accuracy stable:** same basins selected, but effects are sharper within basins (resolution/sharpening gain from annotations)
- **accuracy rises, k_eff stable:** different/better basins selected, same spread (exploration gain)
- **diagonal shift:** both
- **more points leaving the diffuse region:** more effects converging (capacity gain)

This is analyzed as a **joint distribution**, not as independent marginal summaries. The two axes interact (accuracy_ratio is only meaningful for decisive effects), which is visible in the scatter rather than hidden by marginal averaging.

Model-level AUPRC and PVE serve as confirmatory metrics: they catch failure modes the effect-level scatter cannot distinguish (e.g., redundant effects on the same LD block both look good per-effect but waste model capacity).

---

## 4) Oracle Experiments

The goal is to design interventions that surgically remove one error source, leaving others intact. By comparing AUPRC under different oracles, we can attribute improvement to specific error sources.

### Oracle 1: Truth Warm-Start (removes missingness + hallucination)

Already implemented as `warm_method = "truth_warm"` in `run_model.R:1660-1681`.

**What it does:** Initializes the top 10 effects (by |beta|) as one-hot alpha vectors on the 10 strongest causal variants. This covers all 8 graded causals plus 2 ungraded.

**What it removes at initialization:**
- Missingness: all 8 graded causals have a dedicated effect
- Hallucination: no effect starts on a noncausal variant

**What it preserves:**
- Blending: IBSS dynamics may spread alpha mass across LD neighbors during fitting
- Effect drift: effects may migrate away from their oracle assignments during convergence (pulled by ungraded causals or LD proxies)

**Interpretation of AUPRC gap (truth_warm vs perfect):** How much does the IBSS optimizer degrade an oracle initialization? This is a mix of blending and effect-stealing, not blending alone.

**Interpretation of AUPRC gap (vanilla vs truth_warm):** The combined cost of the model not starting in the right basin -- missingness and hallucination at initialization that the optimizer cannot recover from.

### Oracle 2: Blending Removal (UNRESOLVED)

We want an intervention that removes blending while preserving the model's basin selection (including its mistakes). This turns out to be harder than expected.

---

## 5) The Blending Experiment Impasse

### Why post-hoc interventions fail

Three post-hoc approaches were considered and rejected:

**Threshold-based sharpening** (for each effect, if max graded alpha > threshold, sharpen onto that causal): Creates an arbitrary threshold that is hard to justify. Too loose sharpens effects barely touching a causal; too strict misses blended-but-correct effects.

**CS-membership sharpening** (sharpen onto the graded causal if it falls in the effect's 95% CS): For a near-diffuse effect (alpha ~ 1/p everywhere), 95% of all variants are in the 95% CS, so a graded causal is almost always "in the CS" even when the effect carries no real information. The 95% CS criterion does not distinguish real coverage from chance inclusion.

**Peak-sharpening** (replace each effect's alpha with one-hot on its max-alpha variant): Trivially produces near-perfect AUPRC. With 10 effects and 8 graded causals, even moderate accuracy means most graded causals get PIP=1 after sharpening. The oracle does all the work and there is nothing left to measure.

### The core difficulty

Any post-hoc sharpening oracle that decides *which variant to sharpen onto* smuggles in the information we are trying to diagnose. If you sharpen onto the true causal, you have injected the answer. If you sharpen onto the model's peak, you have either done nothing (if it was already decisive) or amplified whatever mistake it was making.

### During-fitting alternatives (not yet evaluated)

Two approaches were discussed but not yet committed to:

**Oracle prior inclusion weights within LD blocks:** For each graded causal, upweight pi_j relative to its LD neighbors (e.g., 10x within r^2 > 0.5 blocks), but do not change relative weights across blocks. This gives within-block resolution without telling the model which blocks to target. Closest to what a real annotation does. Limitation: reduces blending but does not eliminate it; hard to claim this is a clean "remove blending" oracle.

**Oracle LD-pruned design matrix:** For each graded causal, remove all variants with r^2 > threshold from X, leaving only the causal to represent each block. Blending is impossible by construction. Limitation: changes the problem dimensionality and difficulty, so AUPRC comparisons are not apples-to-apples with the full-p problem.

### What's needed

A blending removal experiment that:
1. Does not inject the answer post-hoc
2. Does not change the problem structure (keep the same X and p)
3. Specifically targets within-block resolution without helping basin selection
4. Has a clean interpretation for the AUPRC delta

This remains an open design question.

---

## 6) Planned Analyses (Once Experiments Are Designed)

### Per-method effect-level scatter
For each method, plot (log k_eff, accuracy_ratio) for all non-diffuse effects across all datasets. Compare joint distributions across methods to identify whether gains come from sharpening, basin selection, or both.

### Oracle 2x2 (pending blending oracle)
Run the four combinations of (basin selection: vanilla vs truth_warm) x (blending: preserved vs oracle-removed) for each method. Attribute AUPRC improvement to basin selection vs blending vs interaction.

### Annotation channel attribution
Compare single-fit SuSiE vs single-fit SuSiNE (no ensembling) on the effect-level scatter. Determine whether annotation-driven improvement manifests as:
- More effects becoming decisive (capacity gain)
- Decisive effects becoming more accurate (exploration/basin selection gain)
- Similar accuracy but smaller k_eff (resolution/sharpening gain)

### Tier-stratified analysis
The existing `compute_cs_power_by_tier()` function already computes CS power at different causal tiers (top 3, 5, 7, 9, 11, 23). The effect-level decomposition should be cross-referenced with tier to see whether missingness/blending/hallucination disproportionately affect weaker graded causals (ranks 6-8 by |beta|).
