# SuSiNE robustness analyses — implementation spec

**Three** self-contained analyses for the manuscript revision. Written for analysts working in the SuSiNE codebase who may not have the full paper context. Each section gives the **construction**, the **intuition** (so you can sanity-check your own implementation), and the **exact endpoints to visualize**.

Priorities, in order: (1) **correct annotation scale/ingestion** — triple-checked, see §1.2; (2) **lightweight compute** — reuse fits wherever a quantity doesn't depend on the perturbation; (3) minimal logging — only what the figures need.

The three analyses: **(1)** associated-null contamination (self-gating stress test), **(2)** sparse-architecture sensitivity, **(3)** diffuse-architecture do-no-harm check. Analyses 2 and 3 share machinery (regenerate phenotypes on the same genotypes, same C-CS-vs-baseline comparison). If compute is constrained, **Analysis 3 is the first to cut** — its claim is the most modest.

---

## Background you need (60-second version)

SuSiNE adds a **signed prior mean** `µ0 = c·a` to SuSiE, where `a` is a per-variant signed annotation (effect-size-scale, e.g. from AlphaGenome) and `c` is a scale swept by the ensemble. The per-effect Bayes factor contains an **alignment term ∝ `c·(1−κ)·z_j·ã_j`** that *rewards* agreement between the annotation sign and the marginal association sign `z_j`. This is the "self-gating" property: annotations that disagree with the data are down-weighted automatically.

The reviewer concern we are addressing: self-gating assumes annotation error is **independent** of association error. If a non-causal variant *both* looks associated (high `|z|`, e.g. an LD tag of a true causal) *and* receives an annotation that agrees in sign, the alignment term **rewards the wrong variant** and self-gating cannot save you. **Analysis 1** is an adversarial stress test that injects exactly this, with one knob. **Analyses 2 and 3** (sparse and diffuse architectures) address a separate reviewer concern that our gains may be specific to the oligogenic architecture — testing the sparse tail (few strong causals) and the polygenic tail (many weak causals) respectively.

> Framing for the paper (use this language): Analysis 1 is an **adversarial associated-null contamination stress test**, NOT a generative model of how S2F models behave. It deliberately violates the independence assumption in the most dangerous way to find the boundary of self-gating.

---

## Analysis 1 — Associated-null directional contamination

### 1.1 Construction

Use the **standard oligogenic architecture** and the **central annotation setting (ϕₐ=0.3, νₐ=0.95)**, exactly as in the main simulation study. 600 datasets (150 genotype matrices × 4 phenotype seeds), reused from the main study.

For each dataset:

1. Generate the clean annotation vector `a_orig` with the **current** generator (no changes to the generator).
2. Compute the **frozen marginal signed z-scores** once, before any model fitting:
   ```
   z_j = b_hat_j / se(b_hat_j)     # simple per-variant marginal regression of y on x_j
   ```
   `z` is computed ONCE and never recomputed. Do **not** couple the contamination to any model-derived (posterior) effect estimate — that would be a moving target as IBSS iterates. Marginal, pre-fit, frozen.
3. Let `C` = the true causal set (known in simulation). **Convex-blend the non-causal entries toward the association direction, then renorm the null block to its original RMS.** Causals untouched.
   ```
   # null block only (j not in C); all RMS taken over the null block
   rms_a0   = rms(a_orig[null])                          # original null-annotation RMS
   z_scaled = z[null] * rms_a0 / rms(z[null])            # association direction, matched to null RMS
   blend    = (1 - lambda) * a_orig[null] + lambda * z_scaled
   a[null]  = blend * rms_a0 / rms(blend)                # RENORM: rms(a[null]) == rms_a0 for EVERY lambda
   a[C]     = a_orig[C]                                  # UNCHANGED, byte-for-byte, all lambda
   ```
   The renorm is the key requirement: `rms(a[null])` is **identical across the whole λ sweep**, so λ rotates the null annotations from random toward association-aligned while holding their magnitude — and hence the effective νₐ — fixed. λ is a *dose of alignment*, not a dose of added variance.

**λ grid:** `{0, 0.1, 0.2, 0.3, 0.5, 0.7, 1.0}`. Add `1.5` only if 1.0 does not produce a clear effect (we want to be sure we bracket any crossing).

**Why convex-blend with a renorm (not addition):** the explicit requirement is that the **null-block RMS is constant across λ**. Pure addition would inflate the null variance as λ grows, silently drifting the effective νₐ upward and confounding "aligned vs. random" with "more vs. less null noise." The convex blend rotates direction; the renorm pins magnitude. So λ isolates the one variable we care about — *alignment* — at fixed null-annotation strength. (Aside: a raw convex blend of two near-orthogonal equal-RMS vectors sags ~30% in RMS at λ=0.5; the renorm removes that sag. The cost is that λ is *monotone*, not perfectly *linear*, in alignment angle — off by a few degrees mid-sweep. Claim "monotone dial at fixed magnitude," never "linear.")

**Why nulls-only:** if we blended causals too, the helpful alignment at causals would partly cancel the harmful alignment at nulls, and a flat result would be uninterpretable. Nulls-only makes λ a **pure-harm dial**: any drop is a genuine self-gating failure.

**Why `z` (not raw `b_hat`) as the direction:** the alignment term literally uses `z_j ã_j`, and the dangerous nulls are the statistically-convincing high-`|z|` LD tags. We scale `z` to the null-block RMS so it blends on equal footing with the random annotations. (For this architecture `z` and the standardized marginal effect are nearly proportional, so the choice barely matters — but use `z`.)

**λ=1 anchor:** the null annotations are **fully replaced by the association direction**, at the same RMS the original random null annotations had — i.e. νₐ is held fixed and only the *direction* of the null annotations has rotated from random to association-aligned. λ=0 is exactly baseline; λ=1 is "nulls are now purely aligned signal." (This displacement is deliberate and arguably the more faithful question: what if the annotation model's null predictions are systematically association-aligned rather than noise?)

### 1.2 Annotation scale & ingestion — TRIPLE-CHECK THIS

This is the highest-risk part. The contaminated `a` must enter the fit through the **identical** path as the clean `a`, at the **identical** scale.

Required checks, all must pass:

1. **Trace ingestion in the susine repo.** Find the argument where the annotation enters the prior mean (the `µ0 = c·a` construction inside the single-effect/SER update). Confirm the simulation harness passes `a` straight into `µ0 = c·a` with `c` taken from the grid, with **no hidden per-locus rescaling** (the real-data pipeline has a `c_base` rescaling step from Eqs. 19–21 of the manuscript — confirm the *simulation* path does NOT apply that, or if it does, that contamination is applied at the matching point).
2. **Apply contamination at the same point** in the pipeline that the clean `a` is constructed, then let it flow through unchanged.
3. **MASTER CHECK — λ=0 reproduces baseline bit-for-bit.** Running this analysis at λ=0 must reproduce the existing baseline C-CS AUPRC/recall numbers to numerical tolerance. If it doesn't, ingestion is wrong — stop and fix before doing anything else. This single check catches scale, pipeline, and indexing errors at once.
4. **Fixed-RMS invariant, per locus:** assert `rms(a[null])` equals `rms(a_orig[null])` at *every* λ (to numerical tolerance). This is the defining property of the construction; if it fails, the renorm is computing RMS over the wrong index set or at the wrong step. Also confirm `a[C]` is byte-for-byte identical to `a_orig[C]` at every λ.

### 1.3 Positive-control arm (engagement check) — REQUIRED

Run the identical construction but on the **causal block**, leaving nulls unchanged — the mirror image:
```
# causal block only (j in C); all RMS over the causal block
rms_a0c   = rms(a_orig[C])
z_scaledc = z[C] * rms_a0c / rms(z[C])
blendc    = (1 - lambda) * a_orig[C] + lambda * z_scaledc
a[C]      = blendc * rms_a0c / rms(blendc)        # causal RMS held fixed across lambda
a[null]   = a_orig[null]                          # UNCHANGED
```
Same λ grid, same metrics. (Causal RMS held fixed for the same reason — isolate alignment from magnitude.)

**Purpose / what it must show:** with 23 causal variants, marginal `z` is a badly LD-smeared shadow of the joint signal, so it is *a priori possible* the null-arm contamination simply never engages the joint fit — in which case a flat null-arm curve would be meaningless (non-engagement), not robustness. The positive control is an **engagement check**: rotating causal annotations onto their (mostly correct) association direction should **materially move the reporting metrics, generally toward improvement**, confirming the model is responsive to this kind of perturbation.

> Do NOT gate on strict monotonicity. The same effect that motivates the null arm — `z` being LD-smeared or sign-cancelled under oligogenic effects — can make even the causal arm non-monotone or partially flat at some λ, *even when ingestion is perfect*. A monotonicity gate would therefore reject a correctly-working pipeline for the very reason the experiment is interesting. The bar is "metrics materially change (ideally improve)," not "monotone increase."

**Decision rule:** if the positive-control arm moves materially but the null arm is flat → that flatness is **real robustness**. If even the positive-control arm is flat at λ=1 (where the causal annotations are *fully* rotated onto the association direction) → engagement or ingestion is suspect; investigate before trusting the null arm. Do not report the null-arm result until the positive control has demonstrated engagement.

### 1.4 Shuffled-z control (OPTIONAL — trigger only if the null arm declines)

Only if the null arm shows a decline worth explaining, run the **identical null-arm construction** but with `z` **permuted across the non-causal variants** before it enters the blend (`z[null] → z[perm(null)]`). Same magnitude distribution, variant-level alignment destroyed; null-block RMS still held fixed by the renorm. Run at the 2–3 λ values nearest the crossing only (not the full grid). If harm appears under aligned `z` but vanishes under shuffled `z`, the harm is **alignment-specific**, not a generic consequence of rotating the null annotations. If the null arm is flat, skip this entirely.

### 1.5 Diagnostic: observable contamination signature

For each dataset and each λ, compute:
```
N0 = { j : anchorPIP_j < 0.01 }          # anchor = the annotation-agnostic (c=0) fit's PIP
diagnostic(lambda) = cor( a_lambda[N0], z[N0] )    # signed Pearson over N0
```
- `anchorPIP` comes from the **c=0 / SuSiE-equivalent** fit, which never sees `a`. Therefore `N0` is **automatically λ-invariant** (no freezing needed) and is **observable in real data** (it's just the SuSiE baseline's low-PIP variants). This is the whole reason it transfers sim→real.
- Use the **contaminated** `a_lambda` and the **frozen** `z`.

**Intuition / what it buys us:** the takeaway of Analysis 1 is NOT a bare λ threshold (λ is unobservable in practice, so "don't exceed λ*" is unfalsifiable advice). The takeaway is a **dose–response we can map to an observable**: as λ rises, this diagnostic rises mechanically, and we then check where **real data** sits on that scale (§1.8). Under clean annotations the diagnostic ≈ 0 (null-vs-noise); contamination lifts it. The real question the figure answers: "do real loci sit in the λ range where SuSiNE still beats baseline?"

### 1.6 Ensemble grid & compute reuse

- **Method:** C-CS, **16-run (4×4) grid**, cluster-weight (credible-shift) aggregation — exactly the recommended deployable workflow. The 4×4 **must include the corners of the full 8×8 grid** (so it spans the full range, not a central subgrid). Use the existing Fig-8F 4×4 if it already includes the corners; otherwise use:
  ```
  c    in {0.00, 0.43, 1.07, 1.50}        # includes c=0 (annotation-free) and c_max
  σ0^2 in {0.01, 0.10, 0.40, 1.00}        # includes both extremes
  ```
  (Confirm against the repo's existing subgrid definition — see open questions.)
- **Comparator:** the zero-prior-mean SuSiE-equivalent baseline.

**Compute reuse (do this — it's most of the savings):**
1. **The baseline does not depend on `a`, hence not on λ.** Compute it ONCE; draw it as a single horizontal reference across the λ axis. Do not recompute per λ.
2. **The c=0 members of the C-CS grid do not use `a`.** Compute those 4 fits ONCE per dataset and reuse across all λ. Only the **12 fits with c>0** are recomputed per λ (with the contaminated `a`). Aggregate all 16 each time.
3. The positive-control arm reuses the same c=0 fits too.

### 1.7 Metrics (only these)

- Pooled **AUPRC** (top-8 causal-recovery convention from the main paper).
- **Recall @ precision = 0.75.**
- **Fallback:** if a method's PR curve never reaches 0.75 precision, report `recall@0.75 = 0`, plus `recall@0.50` and `max_achievable_precision`, and flag the cell. (Likely needed only at extreme λ if performance collapses.)

### 1.8 Endpoints to visualize (be exact)

**Figure 1A — main result.** x = λ; y = AUPRC. Three series:
- SuSiE baseline (horizontal reference line; λ-invariant).
- C-CS **null arm** (the stress test).
- C-CS **positive-control arm** (engagement proof).
Annotate the crossing λ\* where the null-arm line drops below baseline, **or** state "no crossing in range." Companion panel with y = recall@75% (same three series).

> The headline sentence is the *regime*, not the threshold: either (a) C-CS stays ≥ baseline across the whole sweep (self-gating holds — annotations stop helping but never hurt), or (b) C-CS crosses below baseline at λ\* (boundary found). The positive-control line is what licenses interpreting a flat null arm as robustness.

**Figure 1B — observable diagnostic.** x = λ; y = `cor(a,z) | anchorPIP<0.01` (sim mean ± band over 600 datasets). This is the dose–response that ties the unobservable λ to a measurable quantity.

**Figure 1B overlay / real-data tie-in (§ below).** On Figure 1B, overlay the **real GTEx-locus diagnostic values** (computed identically) as a rug or points on the y-axis. This places real loci on the **same observable alignment scale** as the simulation. Defensible takeaway: "real loci sit at diagnostic ≤ X, within the range where the simulation shows SuSiNE at or above baseline." Do **not** assign real loci a precise "implied λ" — the diagnostic→λ map is specific to *our* contamination construction (convex rotation toward `z`), and real annotation error need not take that form. The honest claim is co-location on the observable scale, not an inferred contamination dose.

### 1.9 Real-data quick check (lightweight addendum)

For each of the 20 GTEx Lung loci already in the study, compute the **same diagnostic**, using observable quantities only:
```
anchor   = the susie_rss annotation-agnostic baseline fit (already computed)
N0       = { variants with anchor PIP < 0.01 }
z_j      = GTEx marginal z (already available)
a_j      = the AlphaGenome annotation (already available)
real_diagnostic(locus) = cor( a[N0], z[N0] )
```
One number per locus. Report the range, and place them on Figure 1B. No new fitting required — this reuses existing fits and inputs.

---

## Analysis 2 — Sparse-architecture sensitivity

### 2.1 Construction

Separate, simpler check. Reviewer concern: gains may be specific to the 23-causal oligogenic architecture; our own theory (manuscript §3.2) warns the gain "need not stay positive" as the causal fraction shrinks. So test sparse loci.

- **Genotypes:** reuse the same 150 matrices. **Regenerate phenotypes** under a sparse architecture (4 seeds each → 600 datasets per setting).
- **Causal count k ∈ {1, 3, 5}.** (k=1 is the default — but verify the generator first; see note.)
- **Per-variant PVE ∈ {1%, 3%, 5%, 10%}.** 3×4 = 12 settings. `10%` is a deliberate **high-signal stress point** (label it as such); `{1%, 3%, 5%}` span the realistic range. Per-variant PVE (not total `h²`) is the knob **on purpose**: it holds the signal *per causal variant* fixed while only `k` changes, so the "gain vs. sparsity" read is not confounded by per-variant signal weakening. Report total `h² ≈ k × per-variant` as a derived column for reference, but it is not the knob here.
- **Annotation:** central setting (ϕₐ=0.3, νₐ=0.95), generated fresh per dataset with the current generator (no contamination — this analysis is unrelated to Analysis 1).
- **Method/comparator:** C-CS 16-run grid (same as §1.6) vs SuSiE-equivalent baseline.

**k=1 note — verify one implementation fact, then proceed.** There are two distinct "ϕₐ": (i) a *generation parameter* that sets the signal/noise mixing weight when sampling `a` — well-defined at k=1; and (ii) a *measured statistic* `Corr(a_T, β_T)` computed on generated data — undefined at k=1 (correlation over one point). k=1 is fine **unless** the generator calibrates ϕₐ by measuring the realized correlation *per-locus* on that locus's causal variants. Check which of these the current generator does:
- **(a) Sets the mixing weight analytically from the target ϕₐ and samples once** → k=1 fine, proceed.
- **(b) Calibrates by pooling realized correlation across all 600 datasets** (one causal per dataset → 600 points), in the style of the AlphaGenome calibration → k=1 fine, proceed.
- **(c) Calibrates per-locus on that locus's k causal variants** → k=1 breaks (and k=2 is noisy). In this case either switch that locus's draw to an explicit per-variant SNR form (`a_causal = sqrt(ϕₐ)·sign(β)·scale + sqrt(1−ϕₐ)·noise`), or fall back to k∈{2,3,5}.

Expectation is (a) or (b), so k=1 should work directly — but confirm before running, and record which case applies in the methods.

### 2.2 Metrics & positive set

- AUPRC and recall@75% (same fallback rule as §1.7).
- **Positive set = ALL causal variants** (not top-8 — there are only k of them). Do not reuse the top-8 convention here.
- Expect undefined recall@75% cells at the weak end (k low, 1% PVE); apply the fallback and flag.

### 2.3 Endpoints to visualize

**Figure 2 — sparse sweep.** x = k (1, 3, 5); y = AUPRC; one line per method (baseline vs C-CS); faceted by per-variant PVE (4 panels: 1%, 3%, 5%, 10%). Companion: same layout for recall@75% (or fallback). Headline: does the C-CS − baseline gap stay positive, shrink, or reverse as k drops. Match the narrative to whatever actually happens — a shrinking or reversing gap at high sparsity is an honest, expected-by-theory result, not a failure.

---

## Analysis 3 — Diffuse-architecture sensitivity (do-no-harm check)

### 3.1 Framing — read this before coding

Sparse (Analysis 2) tests one tail; this tests the **opposite** tail (many weak causals), which is the second half of the "your gains are oligogenic-specific" reviewer concern.

**Scope the claim correctly.** With ~100 causals and L=10, the fit is heavily under-specified — this regime primarily stresses **SuSiE's own sparsity assumption**, not the µ₀ channel. Do **not** expect or claim a performance *gain* here; both methods may be weak. The defensible target is a **do-no-harm** result: *the signed prior does not make things materially worse than the annotation-agnostic baseline in the polygenic regime.* Pre-commit to that framing so a flat or null result reads as "bounded harm," not "method fails." If compute is tight, this is the **first arm to cut** — its claim is the most modest of the three analyses.

### 3.2 Construction

- **Genotypes:** reuse the same 150 matrices; regenerate phenotypes under the diffuse architecture (4 seeds → 600 datasets per setting).
- **Causal count K = 100** (fixed).
- **Total locus h² ∈ {0.05, 0.15, 0.25}.** 3 settings. Total `h²` **is** the right knob here (unlike sparse): K is fixed, so nothing else is varying for it to confound, and it is the natural "how much signal does this polygenic locus carry" dial. Distribute effects across the 100 causals per the diffuse architecture (e.g. equal or mild-tiered effect sizes), then calibrate residual variance to hit the target `h²` exactly (same `σ² = Var(η)(1−h²)/h²` mechanism as the main study).
- **Annotation:** central setting (ϕₐ=0.3, νₐ=0.95), fresh per dataset, no contamination.
- **Method/comparator:** C-CS 16-run grid vs SuSiE-equivalent baseline.

### 3.3 Metrics & positive set

- AUPRC and recall@75% (same fallback rule as §1.7) — expect frequent fallback at low `h²`, since precision may never reach 0.75 for either method. Flag those cells.
- **Positive set = top-8 causal by |β|** (matches the oligogenic convention; with 100 causals, "all causals" is not a meaningful recovery target).

### 3.4 Endpoints to visualize

**Figure 3 — diffuse sweep.** x = total `h²` (0.05, 0.15, 0.25); y = AUPRC; one line per method (baseline vs C-CS). Companion: recall@75% (or fallback). Headline (do-no-harm framing): "C-CS tracks baseline within noise across the polygenic regime — the signed prior does not degrade performance when the architecture is far more polygenic than designed." If C-CS instead drops materially below baseline at some `h²`, report it plainly as a scope limitation.

---

## Master verification checklist (must pass before reporting anything)

- [ ] **λ=0 reproduces existing baseline C-CS numbers bit-for-bit** (Analysis 1).
- [ ] `rms(a[null])` constant across all λ, and `a[C]` unchanged across all λ (Analysis 1 §1.2 item 4).
- [ ] Annotation ingestion path traced in susine repo; no hidden rescaling in the sim path (or matched if present).
- [ ] Positive-control arm moves the metrics materially (engagement check — NOT a monotonicity gate) before any null-arm result is trusted.
- [ ] Diagnostic computed on **anchor (c=0) PIP<0.01**, with **frozen z** and **contaminated a**.
- [ ] c=0 fits and the baseline computed once and reused across λ.
- [ ] Sparse positive set = all causals; diffuse positive set = top-8 by |β|; recall@75% fallback wired for both.

## Open questions for the authors (resolve before coding)

1. **4×4 subgrid:** does the repo's existing Fig-8F 4×4 already include the 8×8 corners? If not, use the grid in §1.6 — confirm the exact `c` / `σ0²` values.
2. **Sparse k=1:** default is to include it. Confirm the generator is case (a) or (b) in §2.1 (analytic mixing weight, or pooled-across-600 calibration) so no per-locus correlation is computed at k=1. Only if it's case (c) does k=1 need the per-variant SNR form or a k=2 floor.
3. **Shuffled-z control:** run it if the null arm declines (recommended), or skip regardless?
