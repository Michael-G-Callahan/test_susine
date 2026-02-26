# SuSiNE Paper: Analysis Completion — Status, Decisions & TODO

## Context

This document inventories every piece of remaining analysis work (short of paper writing) to complete the SuSiNE paper per the study plan. It is the master planning artifact. Organization:

1. **Critical Decisions** — strategic choices that gate downstream work
2. **SuSiNE Package Updates** — model-level features in `susine/`
3. **Simulation Harness Updates** — `test_susine/` infrastructure
4. **Ensemble / Aggregation Implementation** — weighting and combination methods
5. **HPC Run Plan** — pilot, full study grid, compute accounting
6. **Consolidated TODO** — prioritized action items

Key references:

- `refs/susine_study_plan_high_level.md` (Sections 0-9)
- `refs/susie_susine_background.md` (SuSiE/SuSiE-2/SuSiNE technical background)
- `refs/pip_ensemble_methods.md` (ensemble weighting theory + methods)
- `vignettes/susie_pathology.ipynb` (SuSiE failure mode demonstrations)

---

# CRITICAL DECISIONS

Every decision here blocks or shapes downstream implementation. Marked with status: **OPEN**, **RESOLVED**, or **DEFERRED**.

---

## D1. Package architecture: extend susine vs port into susieR fork vs hybrid

**Status: RESOLVED — Plan A (two-package strategy)**

**Discovery (2025-07-07):** The local susieR fork (v0.12.40) was stale. The upstream stephenslab/susieR has been rewritten as v2.0.0 with SuSiE-ash, SuSiE-inf, PIP convergence, refinement, NIG priors, stochastic LD, and a full S3 class architecture. This fundamentally changed the analysis. See `refs/d1_package_architecture_analysis.md` (revised) and `refs/susieR_2.0_inventory.md` for full details.

**Decision:** Use upstream susieR 2.0 as-is for SuSiE baselines. Keep susine for SuSiNE-specific features. Harness (test_susine) wraps both packages.

| Component | Package | Notes |
| --------- | ------- | ----- |
| SuSiE baseline | susieR 2.0 (upstream) | `susie()` with `refine`, `convergence_method`, `unmappable_effects` |
| SuSiE-ash baseline | susieR 2.0 | `susie(unmappable_effects = "ash")` — individual data only |
| SuSiE-inf baseline | susieR 2.0 | `susie(unmappable_effects = "inf")` — individual + SS |
| SuSiNE (directional prior) | susine | `susine()` with c-grid, sigma_0^2 grid, restarts |
| BFS refinement | test_susine (harness) | Model-agnostic; works with both `susine()` and `susie()` |
| Aggregation (Method A) | test_susine (harness) | Model-agnostic; operates on flat bag of (PIP, ELBO) |
| Alpha/PIP convergence | susine (add) | For parity with susieR 2.0's `convergence_method = "pip"` |

**Rationale (Plan A over Plan B):**
- Most "missing" SuSiE-2 features are now built into susieR 2.0 — zero porting needed
- Plan B (port SuSiNE into susieR 2.0) would require modifying 15+ S3 method implementations across 3 data classes in an unfamiliar architecture (~500–700 lines), with a permanent upstream-tracking burden
- Plan A requires only ~260–370 lines of new code in codebases you control
- SuSiE-ash/inf become free comparison arms (see D3)
- The harness is the true user interface; whether it calls `susine()` or `susie()` is an internal detail

**Immediate action:** `remotes::install_github("stephenslab/susieR")` to replace the stale local fork.

---

## D2. Use case catalog refactoring: separate model spec from exploration method

**Status: RESOLVED — replace old catalog entirely**

**Decision:** Replace the current 14-entry `use_case_catalog()` with a fully factored design. The old catalog conflated model specification with exploration method; the new design separates them into independent axes with an explicit validity constraint.

### Architecture

The pipeline factorizes as:

```text
(prior_spec  x  exploration_method)  →  fits  →  aggregation_method
```

**Interaction structure:**

| Crossing | Interaction? | Notes |
| -------- | ------------ | ----- |
| Exploration × Prior spec | **Yes** | c-grid only valid for SuSiNE specs; tau-grid only for functional-pi specs; sigma_0^2 grid and random restarts apply to all; refinement applies to all but branches differently depending on the prior (different priors → different CSs → different refinement paths) |
| Aggregation × Prior spec | **No** | Aggregator sees a flat bag of (PIP, ELBO) pairs regardless of what prior generated them |
| Aggregation × Exploration | **No** | Same reason — aggregation is grid-agnostic by design (D6) |

This means exploration methods carry a **validity constraint** against prior specs (enforced programmatically), but aggregation is fully independent of both upstream axes.

### Prior specs (3 core)

| ID | Label | Prior mean | Prior variance | Inclusion prior (pi) |
| -- | ----- | ---------- | -------------- | -------------------- |
| `susie` | Vanilla SuSiE | mu = 0 | fixed sigma_0^2 | uniform (1/p) |
| `susie_fpi` | SuSiE + functional pi | mu = 0 | fixed sigma_0^2 | softmax(\|a_j\| / tau) (see D5) |
| `susine` | SuSiNE | mu_0 = c * a | fixed sigma_0^2 | uniform (1/p) |

### Exploration methods (with validity constraints)

| ID | Label | Applies to | Notes |
| -- | ----- | ---------- | ----- |
| `single` | Single fit (baseline) | all prior specs | K = 1 |
| `restart` | K random restarts (Dirichlet alpha init) | all prior specs | |
| `sigma_grid` | sigma_0^2 grid | all prior specs | 3–5 values |
| `c_grid` | c-grid | `susine` only | 5–7 values of c |
| `tau_grid` | tau-grid (softmax temperature) | `susie_fpi` only | see D15 |
| `refine` | SuSiE-2 refinement | all prior specs | cost is variable (depends on # CSs found); behavior differs across prior specs |
| `combined` | Combinations (e.g., restart + c_grid) | validity inherited from components | cross-product of component grids |

### Aggregation methods (fully independent, applied post-hoc)

| ID | Label | Notes |
| -- | ----- | ----- |
| `max_elbo` | Max ELBO (single best) | Baseline; ignores all other fits |
| `uniform` | Uniform weighting | Ignores ELBO; frequency-biased |
| `elbo_softmax` | ELBO softmax | Simple BMA approximation; still frequency-biased |
| `cluster_weight` | Cluster-then-weight (Method A) | Importance-corrects for optimizer frequency via JSD clustering (see `pip_ensemble_methods.md`) |

### Migration plan

- **Replace entirely.** The old 14-entry catalog is retired. No backward-compatibility shim.
- `use_case_catalog()` in `R/use_cases.R` will be replaced by three lookup tables: `prior_spec_catalog()`, `exploration_catalog()`, `aggregation_catalog()`, plus a `valid_exploration_for_prior()` constraint function.
- The run-control pipeline (`make_job_config()` etc.) will accept `prior_spec`, `exploration_method`, and `aggregation_method` as independent arguments, crossing them as requested and filtering by validity.
- Old output directories that reference the retired use-case IDs (a_i, b_ii, c_ii, etc.) remain readable but will not be re-generated.

### Retired use cases (for reference)

The following 14 entries from the old `use_case_catalog()` are retired:

| Old ID | Disposition |
| ------ | ---------- |
| a_i | → `susie` + `single` |
| a_i_restart | → `susie` + `restart` |
| a_ii | → `susie` + `single` with `estimate_prior_variance = TRUE` (EB sigma baseline; see D13) |
| a_iii | → dropped (EB mu not included; see D13) |
| a_iv | → dropped (EB mu & sigma not included; see D13) |
| b_i | → dropped (functional sigma replaced by functional pi; see D5) |
| b_ii | → `susine` + `single` (functional mu is now just the `susine` prior spec) |
| b_iii | → dropped (EB sigma; see D13) |
| b_iv | → dropped (functional sigma; see D5) |
| b_v | → `susine` + `single` with auto mu scale (subsumed into c-grid exploration) |
| b_vi | → dropped (EB sigma; see D13) |
| c_i | → dropped (annealing dropped from paper; see D12) |
| c_ii | → `susine` + `c_grid` + `elbo_softmax` (or other aggregation) |

---

## D3. SuSiE-ash / SuSiE-inf: implement or skip?

**Status: RESOLVED — include as baseline comparison arms (zero implementation cost)**

**Discovery (2025-07-07):** Both SuSiE-ash and SuSiE-inf are fully implemented in upstream susieR 2.0. No custom implementation needed — just pass `unmappable_effects = "ash"` or `"inf"` to `susie()`.

**Decision:** Include SuSiE-ash and SuSiE-inf as **baseline comparison arms** in the simulation study. This strengthens the paper by showing:
1. The multimodality/instability problem is not unique to vanilla SuSiE (failure-mode comparison)
2. Directional priors (SuSiNE) and background-effect models (ash/inf) address complementary problems
3. Ensemble methods work across all paradigms

**Constraints:**
- SuSiE-ash requires **individual-level data only** (not SS or RSS). Our simulations use individual-level data, so this works.
- SuSiE-inf works with **individual + SS** data, but **not RSS (λ>0)**
- Refinement (`refine = TRUE`) is **incompatible** with unmappable effects
- When using ash/inf, convergence is forced to PIP-based (not ELBO)

**Phase D real-data note:** Real data may use RSS where ash/inf are unavailable. Document as a limitation.

**Original recommendation was "skip" based on the (incorrect) belief that implementation would cost 500+ lines each. The actual cost is zero.**

---

## D4. Refinement: include as exploration competitor

**Status: RESOLVED — yes, include**

Study plan Phase B2 lists refinement as exploration mechanism #3. SuSiE-2's `refine=TRUE` explores alternative basins by zeroing out CS SNPs and refitting. This is the existing "SuSiE answer" to multimodality — omitting it would invite reviewer criticism.

**Decision:** Include refinement as an exploration method in the study. It enters the exploration axis alongside random restarts, sigma_0^2 grid, and c-grid. Compute cost is variable (depends on number of CSs found), which needs to be accounted for in the budget-K fairness comparison.

**Implementation path:** Tied to D1 (package architecture). Options:
- Use susieR's built-in `refine = TRUE` for SuSiE baselines
- Port refinement logic into susine for SuSiNE variants
- Hybrid: susieR for SuSiE, custom for SuSiNE

This does NOT need to be resolved before the pilot — the pilot can use susieR's refinement directly for the SuSiE baseline arm.

---

## D5. Functional prior inclusions (pi) as a model spec baseline

**Status: RESOLVED**

**Decision:** Yes. Replace functional prior *variances* (current b_i, b_iv use cases) with functional prior *inclusions*. The pi vector will be driven by softmax on |a_j|:

```
pi_j = softmax(|a_j| / tau)
```

where tau is a temperature parameter that needs tuning/testing.

**Implications:**

- Drop use cases b_i, b_iv (functional sigma) from the catalog
- Add new model spec: "SuSiE + functional pi" (mu=0, fixed sigma_0^2, pi from softmax(|a_j|))
- Strengthens the paper: shows mean shifts (SuSiNE) add value beyond what pi-only approaches provide
- Need to implement softmax(|a_j|) pi generation in `simulate_data.R` or a new helper
- Need to test sensitivity to tau (temperature) — perhaps as part of pilot

**Sub-question:** What values of tau to test? This affects how peaked vs diffuse the pi distribution is. A very low tau concentrates pi on the highest-annotation SNPs; a high tau approaches uniform.

---

## D6. Aggregation design: flatten-based methods

**Status: RESOLVED**

**Decision:** All aggregation methods operate on a **flattened** pool of fits. For a given (dataset, model_spec), all fits — across random restarts, prior grid points, and any other exploration axes — are collected into a single pool and aggregated without regard to the grid structure that generated them. Each fit's ELBO is evaluated under its own prior ("diagonal" evaluation only; no cross-prior re-scoring).

**Four aggregation methods (increasing sophistication):**

| Method | Formula | Notes |
| ------ | ------- | ----- |
| Max ELBO | Select single fit with highest ELBO | Baseline; ignores all other fits |
| Uniform | PIP_ens = (1/K) Σ PIP_i | Ignores ELBO; frequency-biased |
| ELBO softmax | w_i ∝ exp(ELBO_i); PIP_ens = Σ w_i · PIP_i | Simple BMA approximation; still frequency-biased |
| Cluster-then-weight (Method A) | Cluster by JSD, then w_m ∝ exp(ELBO_m) / f_m | Importance-corrects for optimizer frequency; diversity-aware |

Method A is the primary contribution. The first three are baselines/ablations that isolate what the importance correction buys. Max ELBO and uniform are degenerate cases of softmax at temperature → 0 and → ∞ respectively.

**Why flatten (not two-tier):** Keeping aggregation grid-agnostic simplifies the pipeline — the aggregator receives a bag of (PIP, ELBO) pairs regardless of how they were generated (restarts, c-grid, sigma_0^2 grid, or combinations). The grid structure matters at the *exploration* stage (what fits to generate), not at the *aggregation* stage. This also makes the Method A comparison clean: flatten+softmax vs flatten+cluster+importance-corrected softmax directly isolates the value of correcting for optimizer frequency bias.

**Architecture:** Aggregation is post-hoc, operating on saved per-fit PIP vectors. This requires saving all individual fit outputs (PIP, ELBO, prior spec metadata) during HPC runs, not just within-use-case aggregates.

---

## D7. Global ensembles: pool fits across model specs?

**Status: RESOLVED — yes, implement**

**Decision:** For a given (X, y) dataset, pool ALL fits from ALL model specs (SuSiE, SuSiE+pi, SuSiNE across all exploration points) and aggregate. Study plan Phase C mentions a "pooled modes" experiment.

**Pros:** Maximizes diversity; directly tests whether weighting schemes can identify good modes regardless of origin model.

**Cons:** ELBO comparability across different priors is questionable (different priors → different ELBO scales); implementation complexity; could dilute well-specified model's signal.

**Mitigation:** The ELBO-comparability concern is handled by the cluster-then-weight (Method A) aggregation, which uses JSD-based clustering rather than raw ELBO magnitudes for grouping. Within-cluster representative selection still uses ELBO, but only among fits with similar PIP solutions. Additionally, post-hoc analysis can compare within-model-spec aggregation vs global pooling to quantify any dilution effect.

**Implementation:** The flatten-based aggregation architecture (D6) already supports this — the aggregator is agnostic to the origin of fits. Just widen the bag to include fits from all model specs for a given dataset. Tag each fit with its model spec for diagnostic breakdowns.

---

## D8. Compute budget definition and fairness

**Status: RESOLVED — fixed K (number of fits)**

**Decision:** Budget is defined as a fixed number of model fits K. Walltime varies too much across HPC nodes to be a reliable budget metric.

**Fairness subtlety:** If SuSiNE searches a (c, sigma_0^2) grid of size 5x3 = 15, while SuSiE uses 15 random restarts, is that "fair"? SuSiNE gets 15 distinct prior specs; SuSiE gets 15 noisy samples from one prior spec. The compute-confound pilot (D10) is designed to address this.

**Wall time recording:** Record wall time per fit in all HPC outputs so that post-hoc analyses can verify that per-fit cost is comparable across methods and report any systematic differences.

**Specific K values:** Defer to the run control workbook. Refactor `run_controls.R` / `make_job_config()` if K is not already cleanly parameterized as a top-level budget argument.

---

## D9. Z-score metric selection

**Status: RESOLVED — current four candidates sufficient**

**Decision:** The four implemented z-score metrics (`z_topk_ratio`, `z_max_abs`, `z_count_abs_gt_3`, `z_eff_signals`) are sufficient. No need to implement optional kurtosis or entropy metrics. Primary z-metric selection remains a post-hoc analysis step after the pilot.

---

## D10. Pilot study scope

**Status: RESOLVED — run a pilot for benchmarking exploration + aggregation**

**Decision:** Run a pilot study on a smaller set of datasets (5–10, spanning M1 difficulty range) with the primary purpose of benchmarking exploration methods and aggregation methods. This is a substantive benchmarking run, not just a debugging smoke test — pipeline bugs should be caught in local testing before submitting to HPC.

**Pilot scope:**

- 5–10 datasets selected by M1 spread
- Full exploration grid: random restarts, sigma_0^2 grid, c-grid (SuSiNE), tau-grid (functional pi), refinement
- All aggregation methods applied post-hoc
- Compute-confound comparison: restarts vs c-grid diversity at matched K
- Secondary outputs: compute credit estimation for the full study, z-metric selection data

**Study plan Section 5 suggests:** fixed sigma_0^2 = 0.2, K=12 restarts for SuSiE vs K=12 c-values for SuSiNE.

---

## D11. Real-data Phase D: in scope now or later?

**Status: RESOLVED — plan now, execute after simulation insights**

**Decision:** Plan the Phase D implementation now (pipeline + data prep), but execution follows naturally from simulation insights. The workflow is relatively straightforward once we have the simulation results to guide locus selection and method choice.

**Phase D steps (from study plan):**

- Screen a large set of loci by difficulty metrics
- Select 3–4 additional loci beyond CBX8 by (M1, Z_primary) spread
- Run ensemble workflows on each
- Minimal validation: summary-stat matches individual-level on at least one locus

**Note:** Phase D may require LD regularization (D1 dependency) for real LD matrices.

---

## D12. Annealing: keep as exploration method or drop?

**Status: RESOLVED — dropped from paper**

**Decision:** Annealing is dropped entirely from the simulation study. The susine package retains its annealing implementation (`c_i` use case) for potential future use, but the paper will not include annealing as an exploration method. This reduces study complexity and keeps the exploration axis focused on the mechanisms central to the paper's argument: random restarts, prior grids (sigma_0^2, c, tau), and optionally refinement.

---

## D13. EB prior updates: include as a model spec axis or drop?

**Status: RESOLVED — keep EB prior variance as a use case (for prior variance only)**

**Decision:** Include EB prior variance estimation (`estimate_prior_variance = TRUE`) as a use case in the study, but only for prior variance — not for prior mean. EB mu (a_iii) and EB both (a_iv) are dropped. This serves as the "what practitioners currently do" baseline to contrast against fixed-grid ensembling.

**Rationale:** EB sigma is the SuSiE/susieR default. Omitting it entirely would be odd and invite reviewer questions. Including it — even if only to show it doesn't perform well for ensemble weighting — strengthens the argument for fixed sigma_0^2 grids. The pathology vignette already showed EB distorts ELBO-softmax weights (0.21/0.79 instead of 0.50/0.50 for symmetric basins).

**Pilot plan:** Include EB prior variance in the pilot study (D10) to quantify the distortion on pilot data and confirm it doesn't perform well for aggregation.

---

## D14. Annotation quality / error model

**Status: RESOLVED — current controls sufficient, no sign flips or LD-correlated noise**

**Decision:** Keep annotation quality controls as-is (`annotation_r2`, `inflate_match`, `gamma_shrink` in `simulate_priors()`). No need for explicit sign-flip scenarios or LD-correlated annotation noise.

**Paper note:** The manuscript should briefly explain why we don't model sign flips or LD-correlated noise — the `annotation_r2` parameter already captures annotation quality degradation in a smooth, interpretable way, and the study's focus is on the ensemble/exploration framework rather than annotation error modeling.

---

## D15. Softmax temperature grid for functional pi

**Status: RESOLVED — tau-grid is an exploration knob for functional-pi specs**

With pi_j = softmax(|a_j| / tau), the temperature tau controls how peaked the inclusion prior is. **tau is designated as an exploration knob** — a grid over tau values is an exploration method for the `susie_fpi` model spec, analogous to the c-grid for SuSiNE and the sigma_0^2 grid for all models.

**Implementation needs:**

1. Sensible tau grid values (low tau = peaked on high-annotation SNPs; high tau ≈ uniform) — determine empirically from pilot
2. Testing on pilot data to understand sensitivity and useful range
3. Integration into the exploration method axis of the refactored pipeline

**Note:** tau-grid applies only to model specs that use functional pi (`susie_fpi`). It is irrelevant for vanilla SuSiE (uniform pi) and SuSiNE with uniform pi.

---

## D16. Number of seeds / replicates per phenotype setting

**Status: RESOLVED — defer specifics; add seed management review to pre-pilot checklist**

**Decision:** Seed management is complex and touches many parts of the pipeline (phenotype generation, initialization, restart seeds, etc.). Rather than locking in a number now, add a **seed management review** to the TODO list immediately before kicking off pilot HPC runs. This review should verify reproducibility, confirm seed propagation through the pipeline, and set the number of replicates based on pilot compute estimates.

---

## D17. Specific grid values for sigma_0^2 and c

**Status: RESOLVED — defer to run control workbook**

**Decision:** Specific grid values for sigma_0^2 and c are run-level parameters, not high-level architectural decisions. Set them as inputs in the run control workbook when configuring each study (pilot, full). No need to lock them in here.

---

# IMPLEMENTATION STATUS

---

## 1. SuSiNE Package Updates (susine/)

### 1.1 Refinement — D1 RESOLVED, D4 RESOLVED

- **Status:** NOT IMPLEMENTED
- **Implementation path (per D1 resolution):**
  - **SuSiE baselines**: Use susieR 2.0's built-in `refine = TRUE` as one exploration arm
  - **All models (SuSiE + SuSiNE)**: BFS refinement tree in test_susine harness (model-agnostic, budget-controlled, saves all fits for ensembling)
- **Effort:** Low–Med (~80–120 lines) for harness-level BFS. Zero for susieR's built-in.
- **Files:** `test_susine/R/run_model.R` (or new `refinement.R`); depends on `get_credible_set()` in `evaluation_metrics.R`

### 1.2 Alpha-based convergence criterion

- **Status:** NOT IMPLEMENTED (ELBO-only)
- **What:** Stop when max |delta alpha| < tol. Study plan calls for this for cross-method comparability.
- **Effort:** Low. Add `convergence_method` param, track max alpha change per iteration.
- **Files:** `susine/R/susine.R`, `susine_ss.R`.

### 1.3 SuSiE-ash / SuSiE-inf — D3 RESOLVED

- **Status:** Available in upstream susieR 2.0 — no implementation needed.
- **Usage:** `susie(X, y, unmappable_effects = "ash")` or `susie(X, y, unmappable_effects = "inf")`
- **Constraints:** ash = individual data only; inf = individual + SS; neither works with RSS (λ>0)
- **Role:** Baseline comparison arms to show instability is not unique to vanilla SuSiE, and that directional priors (SuSiNE) address a complementary problem to background-effect models.

### 1.4 LD regularization — blocked on D11

- **Status:** BASIC (symmetrization only). No spectral shrinkage.
- **Priority:** LOW for simulations (true genotype matrices). MEDIUM for real-data Phase D.
- **Files:** `susine/R/susine_rss.R`, `initialize.R`.

### 1.5 Aggregation method expansion in susine package

- **Status:** `aggregate_susine_fits()` supports `"equal"` and `"elbo"` (softmax with temperature).
- **Missing:** Cluster-then-weight (Method A). See Section 3.

### 1.6 Functional pi support in susine

- **Status:** NEEDS VERIFICATION. susine accepts `pi` parameter but need to confirm it flows through correctly to SER updates and doesn't get overwritten by uniform defaults.
- **Action:** Test with non-uniform pi vector.

### 1.7 Computational speedups in susine IBSS loop

- **Status:** NOT IMPLEMENTED — profiled 2025-07-07, pure-R fixes, no compiled code needed.
- **Priority:** MEDIUM. Reduces per-fit wall time by ~40–50%, directly relevant to the full study's ~2.7M model fits.
- **Context:** The SuSiE 2.0 paper's "5x speed improvement" is RSS-specific (regularized LD eigendecomposition reuse) and does not apply to individual-level data. However, profiling susine's own IBSS loop revealed it performs ~7L O(np) matrix-vector products per iteration when ~3L–4L suffice. The fixes below are independent and can be applied incrementally.

#### 1.7.1 Eliminate redundant `compute_Xty` in `BF()` (highest impact)

- **What:** `SER()` calls `compute_Xty(X, r)` at `susine/R/SER.R` ~line 24 to get `Xty`, which feeds into the posterior mean/variance updates. Then `BF()` (called from within SER at ~line 37, defined at ~line 75) calls `compute_Xty(X, r)` **again** with the same `X` and `r` to compute log Bayes factors. This is the identical O(np) matrix-vector product.
- **Fix:** Pass the already-computed `Xty` vector from `SER()` into `BF()` as an argument instead of recomputing it. This eliminates L redundant O(np) mat-vec products per IBSS iteration.
- **Impact:** ~25% wall-clock reduction (single largest win). Eliminates L of the ~7L mat-vec products per iteration.
- **Files:** `susine/R/SER.R` — modify `BF()` signature to accept `Xty`; update the call site in `SER()`.
- **Test:** Verify ELBO and PIP output is bit-identical before/after on a test case.

#### 1.7.2 Compute ERSS once per effect, reuse for ELBO and variance update

- **What:** Expected residual sum of squares (ERSS) is computed in `SER_ERSS()` (`susine/R/ERSS.R`) for the variance/sigma update step, and a functionally equivalent quantity is computed again during the ELBO calculation. Both involve O(np) work.
- **Fix:** Compute ERSS once in the SER update, cache it, and pass the cached value to both the ELBO computation and the variance update.
- **Impact:** Eliminates ~L redundant O(np) products per iteration.
- **Files:** `susine/R/ERSS.R`, `susine/R/SER.R`, `susine/R/susine.R` (where ELBO is assembled).
- **Test:** Confirm ELBO values unchanged.

#### 1.7.3 Eliminate redundant `Xb` computation in `SER_ERSS()`

- **What:** `SER_ERSS()` computes `X %*% b_hat` (the fitted values for the current effect) to get the expected residual. But `SER()` already computes the posterior mean vector `mu1` and alpha weights, from which `Xb = X %*% (alpha * mu1)` can be derived. The same product (or a closely related one) is recomputed inside `SER_ERSS()`.
- **Fix:** Have `SER()` return the intermediate `Xb` (or `X %*% (alpha * mu1)`) and pass it into `SER_ERSS()` instead of recomputing.
- **Impact:** Eliminates ~L O(np) products per iteration.
- **Files:** `susine/R/SER.R`, `susine/R/ERSS.R`.

#### 1.7.4 Carry residual forward across outer iterations (incremental update)

- **What:** At the start of each outer IBSS iteration, `susine.R` ~line 166–167 recomputes the full residual `r = y - X %*% fitted` from scratch. Within each iteration, the inner loop already does incremental add-back/subtract for each effect ℓ (adds back ℓ's contribution, fits SER on residual, subtracts new contribution). But between iterations, the accumulated residual is thrown away and rebuilt.
- **Fix:** After the inner loop completes, the residual `r` already reflects the current fit. Carry it forward to the next iteration instead of recomputing. Only recompute from scratch every N iterations (e.g., every 5–10) as a numerical stability check.
- **Impact:** Eliminates 1 O(np) mat-vec product per outer iteration (minor compared to 1.7.1–1.7.3, but free).
- **Files:** `susine/R/susine.R`.
- **Caution:** Floating-point drift over many iterations. Add an optional `recompute_residual_every` parameter (default ~10) that does a full recompute periodically.

#### 1.7.5 Precompute `diag(X'X)` if not already cached (defensive)

- **What:** `susine` caches `diag(X'X)` as `attr(X, "d")` in `standardize_x()`. Verify this is consistently used everywhere and never recomputed.
- **Fix:** Audit all callsites in `SER.R`, `BF.R`, `ERSS.R` to confirm they read from `attr(X, "d")` rather than recomputing column norms.
- **Impact:** Likely already handled, but worth verifying.
- **Files:** `susine/R/SER.R`, `susine/R/ERSS.R`, `susine/R/matrix_computation.R`.

#### Summary: expected reduction

| Source | Current cost (per IBSS iter) | After fix |
|--------|-----------------------------|-----------|
| SER: compute_Xty | L × O(np) | L × O(np) (kept) |
| BF: compute_Xty (redundant) | L × O(np) | **0** (1.7.1) |
| ERSS + ELBO (redundant) | ~L × O(np) | **0** (1.7.2) |
| SER_ERSS: Xb (redundant) | ~L × O(np) | **0** (1.7.3) |
| Residual rebuild | 1 × O(np) | **0** (1.7.4, amortized) |
| **Total** | **~(4L+1) × O(np)** | **~L × O(np)** + small overhead |

For typical L=10, n=600, p=1000: ~41 → ~10 mat-vec products per iteration. Expect ~40–50% wall-clock improvement depending on iteration count and overhead.

---

## 2. Simulation Harness Updates (test_susine/)

### 2.1 Oligogenic phenotype architecture

- **Status:** PARTIALLY SUPPORTED
- **What exists:** `simulate_effect_sizes(p, p_star, effect_sd, seed)` draws Normal(0, effect_sd). p_star grid supports any value.
- **What's missing:** Heavy-tailed / mixture effect distributions (study plan says optional). No explicit h^2 parameterization (uses `noise_fraction = 1 - h^2`).
- **Priority:** Current Normal effects + p_star grid likely sufficient for main text. Heavy-tailed is supplementary.
- **File:** `test_susine/R/simulate_data.R`

### 2.2 Functional pi generation — NEW, from D5

- **Status:** NOT IMPLEMENTED
- **What:** Generate pi_j = softmax(|a_j| / tau) from annotation vector a.
- **Where:** New helper in `simulate_data.R` or `simulate_priors()`.
- **Needs:** Temperature parameter tau; testing sensitivity; integration with run pipeline.
- **File:** `test_susine/R/simulate_data.R`

### 2.3 Exploration methods — current inventory

| Method | Code Status | Pipeline Status | Notes |
| ------ | ----------- | --------------- | ----- |
| Single fit (baseline) | WORKING | In catalog as `extra_compute = "none"` | Default path |
| Random restarts (Dirichlet alpha) | WORKING | In catalog as `extra_compute = "restart_alpha"` | n_inits from restart settings |
| Model averaging (multi-init) | WORKING | In catalog as `extra_compute = "model_avg"` | Averages PIPs across inits |
| **sigma_0^2 grid** | **NOT IMPL** | **NOT IMPL** | Central to Phase B; needs job-config-level grid |
| **c-grid (annotation scale)** | **NOT IMPL** | **NOT IMPL** | Central to Phase B; mu_0 currently pre-computed, not gridded |
| **tau-grid (functional pi temp)** | **NOT IMPL** | **NOT IMPL** | D2/D15: softmax temperature as exploration knob for functional-pi models |
| **Refinement** | **NOT IMPL** | **NOT IMPL** | D1+D4 RESOLVED: susieR 2.0 built-in `refine=TRUE` for SuSiE baselines; BFS in harness for all models |
| ~~Annealing/tempering~~ | WORKING (in susine) | ~~In catalog~~ | **D12 RESOLVED: dropped from paper** |

**Critical gaps:** sigma_0^2 grid and c-grid are the two exploration mechanisms central to the paper's argument but are not wired as exploration axes.

### 2.4 Use case catalog refactoring (D2 RESOLVED)

- **Current:** 14 use cases mixing model spec + exploration.
- **Needed:** Replace with `prior_spec_catalog()` x `exploration_catalog()` x `aggregation_catalog()` + `valid_exploration_for_prior()` constraint function.
- **Files:** `use_cases.R`, `run_controls.R`, `run_model.R`.
- **Effort:** Medium-high. Touches the core pipeline.

### 2.5 Dataset difficulty metrics — ALL IMPLEMENTED

All metrics from study plan Section 4.2 are implemented in `dataset_metrics.R`:

**LD-only (X-matrix):**

- M1 (mid-LD energy) — primary
- M2, M1_dist, ECS, High_LD_mass, Top5_ev_frac, NSI
- Block structure metrics (tau=0.5, 0.8)
- Spectral metrics (participation ratio, stable rank, gap)
- Row concentration (Herfindahl)
- Submodularity bounds

**Z-score-only:**

- z_topk_ratio (concentration, k=10)
- z_max_abs (peak magnitude)
- z_count_abs_gt_3 (tail count)
- z_eff_signals (effective signal count)

**Action:** Verify all are computed and *saved* during HPC runs (not just available as functions). Check that z-score metrics are computed from the simulated y/X data during `run_task()`.

**D9 RESOLVED:** Current four z-score metrics are sufficient. No additional metrics needed.

### 2.6 Multimodality metrics — ALL IMPLEMENTED

In `run_model.R` `compute_multimodal_metrics()`:

- mean_jsd, median_jsd, max_jsd (pairwise Jensen-Shannon divergence)
- jaccard_top10 (top-10 PIP overlap)
- mean_pip_var (per-SNP PIP variance across fits)
- n_clusters (hierarchical clustering at JSD threshold)

**Minor issue:** Uses average linkage; `pip_ensemble_methods.md` recommends complete linkage. Consider switching.

**Action:** Verify these are computed for ALL multi-fit groups (not just restart_alpha use cases).

---

## 3. Ensemble / Aggregation Methods

### 3.1 Currently implemented

| Method | Where | Status |
| ------ | ----- | ------ |
| Max ELBO (single best) | `run_model.R` `aggregate_use_case_pips()` | WORKING |
| Uniform / mean | `run_model.R` + `aggregate_fits.R` | WORKING |
| ELBO softmax | `run_model.R` + `aggregate_fits.R` | WORKING; temperature in susine but NOT exposed in test_susine job config |

### 3.2 Not yet implemented

| Method | Priority | Effort | Description |
| ------ | -------- | ------ | ----------- |
| **Cluster-then-weight (Method A)** | HIGH | Low-med | Clustering infra exists (JSD + hclust). Need: importance correction w_m ~ exp(ELBO_m) / f_m, representative selection per cluster, ESS diagnostic. Ref: `pip_ensemble_methods.md` |
| **Temperature-exposed softmax** | MEDIUM | Trivial | Wire temperature param from job config through to aggregation. |
| **KDE importance sampling (Method B)** | LOW | Medium | Only for K=100+ ensembles. Skip unless we run large restart counts. |
| **Global cross-model-spec pooling** | HIGH (D7 RESOLVED) | Medium | Pool all fits across model specs for a dataset. |

### 3.3 Aggregation architecture (RESOLVED per D6)

All aggregation is **post-hoc** and **grid-agnostic**. For a given (dataset, model_spec), the aggregator receives a flat bag of (PIP vector, ELBO) pairs from all fits — regardless of whether they came from restarts, prior grids, or combinations. No two-tier or grid-aware aggregation logic is needed.

**Implication:** We must save per-fit PIP vectors and ELBOs (not just within-use-case aggregates) from HPC runs. This may already happen with `verbose_file_output` but needs verification.

**Four methods to implement (see D6 for details):** Max ELBO, Uniform, ELBO softmax, Cluster-then-weight (Method A). The first three are already implemented; Method A is the primary new implementation need.

---

## 4. HPC Run Plan

### 4.1 Pilot (D10 RESOLVED)

**Recommended:** 5-10 datasets spanning M1 difficulty, full exploration grid, serves as pipeline validation + compute-confound pilot.

**Compute-confound pilot spec (study plan Section 5):**

- 5-10 datasets selected by M1 spread (and optionally max|z|)
- Fixed sigma_0^2 = 0.2
- SuSiE: K=12 random restarts (Dirichlet alpha)
- SuSiNE: K=12 c-values (e.g., c in {0, 0.1, 0.2, ..., 1.1})
- Compare: mean/median/max JSD, cluster count, PIP variance, PIP correlation heatmaps
- Report: fraction of datasets with D_c-grid > D_restart, mean ratio

### 4.2 Full study design — blocked on D1, D2

**Genotype matrices:** Up to 150 loci (from first paper simulations).

**Phenotype grid:**

- p_star: {1, 3, 5, 10} — sparse to oligogenic
- noise_fraction: {0.5, 0.8, 0.95} — h^2 = {0.5, 0.2, 0.05}
- effect_sd: scaled to achieve target h^2 given X (or fixed and let h^2 vary?)
- seeds: 5-10 per setting (D16)

**Annotation grid:**

- annotation_r2: {0.0, 0.1, 0.3, 0.5} — no annotation to strong
- Note: annotation_r2 = 0.0 means SuSiNE has no useful directional info (important control)

**Model specs (5 arms, per D2/D3/D5):**

1. **Vanilla SuSiE** — mu=0, fixed sigma_0^2, uniform pi (via susieR 2.0)
2. **SuSiE + functional pi** — mu=0, fixed sigma_0^2, pi from softmax(|a_j| / tau) (via susieR 2.0)
3. **SuSiNE** — mu_0 = c*a, fixed sigma_0^2, uniform pi (via susine)
4. **SuSiE-ash** — vanilla SuSiE + Mr.ASH background (via susieR 2.0, `unmappable_effects = "ash"`, individual data only)
5. **SuSiE-inf** — vanilla SuSiE + infinitesimal background (via susieR 2.0, `unmappable_effects = "inf"`, individual + SS)

**Exploration methods (per model spec, matched budget K; composable):**

- Single fit (control)
- K random restarts
- sigma_0^2 grid (D17)
- c-grid for SuSiNE (D17)
- tau-grid for functional-pi models (D15)
- Refinement (D4 RESOLVED: included)
- Combinations (e.g., random restarts + c-grid at matched total K)

**Aggregation (post-hoc, flatten all fits per (dataset, model_spec) — see D6):**

- Max ELBO, Uniform, ELBO softmax, Cluster-then-weight (Method A)

**Rough compute estimate:**

- 150 loci x 4 p_star x 3 noise x 4 annot_r2 x 5 seeds = 36,000 datasets
- 3 model specs x ~5 exploration configs x ~5 fits each = ~75 fits per dataset
- Total: ~2.7M model fits
- Per fit: ~1-5 seconds (n=600, p=1000) → ~750-3750 CPU-hours
- With SLURM array parallelism: feasible in days

### 4.3 Output storage requirements

Per HPC run, we need to save:

- Per-fit: PIP vector (length p), ELBO, sigma_2, convergence info
- Per-group: aggregated PIPs (per aggregation method), multimodality metrics
- Per-dataset: difficulty metrics (M1, z-score metrics)
- Per-dataset-model-spec: all individual fit PIPs (for post-hoc aggregation flexibility)

**Disk estimate:** ~2.7M fits x p=1000 floats x 8 bytes = ~21 GB for PIPs alone. Plus metadata. Manageable.

---

## 5. Consolidated TODO List

### Phase 0: Decisions (resolve these first)

- [x] D1: Package architecture — **RESOLVED: Plan A (two-package: upstream susieR 2.0 + susine)**
- [x] D2: Use case refactoring approach — **RESOLVED: replace old catalog with prior_spec × exploration × aggregation (see D2 section)**
- [x] D3: SuSiE-ash/inf — **RESOLVED: include as baseline comparison arms (free in susieR 2.0)**
- [x] D4: Refinement — **RESOLVED: yes, include as exploration competitor**
- [x] D5: Functional pi — **RESOLVED: yes, softmax(|a_j| / tau)**
- [x] D6: Aggregation design — **RESOLVED: flatten + softmax / Method A**
- [x] D7: Global ensembles — **RESOLVED: yes, implement cross-model-spec pooling**
- [x] D8: Compute budget — **RESOLVED: fixed K (number of fits); record wall time per fit**
- [x] D9: Z-score metrics — **RESOLVED: current four candidates sufficient**
- [x] D10: Pilot scope — **RESOLVED: benchmarking exploration + aggregation on 5–10 datasets**
- [x] D11: Real-data Phase D — **RESOLVED: plan now, execute after simulation insights**
- [x] D12: Annealing — **RESOLVED: dropped from paper**
- [x] D13: EB variants — **RESOLVED: keep EB prior variance as practitioner-default baseline**
- [x] D14: Annotation error model — **RESOLVED: current controls sufficient, no sign flips**
- [x] D15: Softmax temperature — **RESOLVED: tau-grid is exploration knob for functional-pi specs**
- [x] D16: Seeds — **RESOLVED: seed management review before pilot HPC runs**
- [x] D17: Grid values — **RESOLVED: defer to run control workbook**

### Phase 1: Foundation (after decisions)

- [ ] Refactor use_case_catalog into model_spec x exploration x aggregation (D2)
- [ ] Implement sigma_0^2 grid exploration in harness
- [ ] Implement c-grid exploration for SuSiNE in harness
- [ ] Implement functional pi generation: softmax(|a_j| / tau) (D5)
- [ ] Implement tau-grid exploration for functional-pi models (D15)
- [ ] Wire softmax temperature through job config to aggregation
- [ ] Implement flatten-based post-hoc aggregation pipeline (D6)
- [ ] Verify susine accepts and correctly uses non-uniform pi
- [ ] Add wall time recording per fit to HPC outputs (D8)
- [ ] Verify K is cleanly parameterized as top-level budget in run controls (D8)

### Phase 2: Model features

- [ ] Implement refinement as exploration method (D4; implementation path tied to D1)
- [ ] Implement alpha-based convergence in susine
- [ ] Implement cluster-then-weight aggregation (Method A)
- [ ] Confirm oligogenic architecture coverage (p_star grid + noise_fraction)
- [ ] Test functional pi sensitivity to tau (D15)
- [ ] **susine speedup 1.7.1:** Eliminate redundant `compute_Xty` in `BF()` — pass `Xty` from `SER()` (see §1.7.1)
- [ ] **susine speedup 1.7.2:** Compute ERSS once per effect, reuse for ELBO + variance update (see §1.7.2)
- [ ] **susine speedup 1.7.3:** Eliminate redundant `Xb` in `SER_ERSS()` — pass from `SER()` (see §1.7.3)
- [ ] **susine speedup 1.7.4:** Carry residual forward across outer iterations (see §1.7.4)
- [ ] **susine speedup 1.7.5:** Audit `diag(X'X)` caching — verify no redundant column-norm recomputation (see §1.7.5)

### Phase 3: Pipeline validation

- [ ] Verify all dataset difficulty metrics computed + saved during HPC runs
- [ ] Verify all multimodality metrics computed for all multi-fit groups
- [ ] Verify per-fit PIP vectors are saved (needed for post-hoc aggregation)
- [ ] Consider complete vs average linkage for JSD clustering
- [ ] Run 2-locus quick pipeline smoke test

### Phase 4: Pilot study

- [ ] Seed management review: verify reproducibility + seed propagation through pipeline (D16)
- [ ] Select 5-10 datasets by M1 spread
- [ ] Include EB prior variance as practitioner-default baseline in pilot (D13)
- [ ] Run compute-confound pilot (restarts vs c-grid at matched K)
- [ ] Analyze diversity metrics: JSD, clusters, PIP variance
- [ ] Estimate full-study compute cost from pilot timing
- [ ] Validate aggregation methods on pilot data
- [ ] Test global cross-model-spec pooling on pilot data (D7)

### Phase 5: Full simulation study

- [ ] Finalize full study grid (model specs x exploration x aggregation x phenotype)
- [ ] Generate and submit SLURM job array
- [ ] Post-hoc aggregation and metric computation
- [ ] Z-metric selection protocol: choose primary z-metric (D9)
- [ ] Failure-mode stratification by (M1, Z_primary) bins
- [ ] Generate paper exhibits (power-FDR curves, AUPRC boxplots, failure map, diversity scatter)

### Phase 6: Real data (D11 RESOLVED — plan now, execute after simulation)

- [ ] Screen loci by difficulty metrics
- [ ] Select 3-4 additional loci by (M1, Z_primary) spread
- [ ] Run ensemble workflows on selected loci
- [ ] Summary-stat vs individual-level validation
- [ ] LD regularization for real LD matrices (if needed)

### Deferred

- ~~SuSiE-ash / SuSiE-inf (D3)~~ → **RESOLVED: included as baseline comparison arms via susieR 2.0 (zero cost). See D3, Phase 4.2 arm 4-5.**
- [ ] KDE importance sampling Method B (only for K=100+)
- [ ] Heavy-tailed effect distributions
