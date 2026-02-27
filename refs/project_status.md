# SuSiNE Project Status

**Last updated:** 2026-02-27
**Supersedes:** `analysis_completion_status.md` (which remains as archival reference for resolved decisions D1-D17)

---

## Agent Implementation Protocol

Before implementing any task item (S* or T*), an agent must:

1. **Read the relevant source code** in all applicable repos (test_susine, susine, susieR) to understand the current state — function signatures, data flow, naming conventions, existing patterns.
2. **Assess clarity**: decide whether the task description, combined with the code you've read, gives enough information to implement within reasonable confidence of the author's intent.
3. **If unclear, stop and ask.** Do not guess at ambiguous design choices. Flag the specific ambiguity as a question and wait for guidance before writing code.
4. **For the coupled refactor T1+T2+T12**: these are one logical unit. Do not implement them independently — read all three descriptions, understand the current pipeline end-to-end, and plan the refactor as a single coherent change.

---

## Section 1: High-Level Requirements

### 1.1 Pilot Study

**Purpose:** Benchmark exploration methods and aggregation methods on a small dataset set. Validate the full pipeline before committing to a large HPC run.

**Scope:**
- 5-10 genotype matrices, sampled to stratify by M1 (LD complexity)
- L = 10 for all models
- All use cases (Section 2 below), all exploration methods, all aggregation methods
- Outputs: exploration diversity comparison, aggregation method comparison, use case comparison, compute cost estimates

**Key outputs:**
1. How do exploration methods compare at multimodality metrics? At ensemble performance?
2. How does the aggregation method compare at ensemble performance?
3. How does the use case compare at individual-level model performance? At ensemble performance?
4. Compute credit estimation for the full study
5. Z-metric selection data (which z-score metric best predicts instability)

### 1.2 Full Simulation Study

**Purpose:** Comprehensive evaluation of SuSiNE's ensemble framework across the full design grid.

**Scope:**
- Up to 150 genotype matrices (from first-paper simulations)
- Phenotype grid: p_star x y_noise x (architecture?) x annotation settings x seeds
- All use cases, explorations, and aggregations
- Failure-mode stratification by (M1, Z_primary) bins
- Paper exhibits: power-FDR curves, AUPRC boxplots, failure map, diversity scatter

### 1.3 Real Data Analysis (Phase D)

**Purpose:** Apply simulation-learned workflows to real eQTL loci.

**Scope:**
- Screen a large set of loci by difficulty metrics (M1, Z_primary)
- Select 3-4 additional loci beyond CBX8 by difficulty spread
- Run ensemble workflows on selected loci
- Minimal validation: summary-stat implementation matches individual-level on at least one locus
- LD regularization may be needed for real LD matrices

**Timing:** Plan now, execute after simulation insights guide locus selection.

---

## Section 2: Pilot Study Detailed Requirements

### 2.1 Use Cases

All models use L = 10.

| ID | Label | Package | Prior mean | Prior variance | Inclusion prior (pi) | Notes |
|----|-------|---------|------------|----------------|---------------------|-------|
| 1 | SuSiE vanilla | susieR | mu = 0 | fixed sigma_0^2 | uniform | Primary baseline |
| 2 | SuSiE-ash | susieR | mu = 0 | fixed sigma_0^2 | uniform | `unmappable_effects = "ash"`, individual data only |
| 3 | SuSiE-inf | susieR | mu = 0 | fixed sigma_0^2 | uniform | `unmappable_effects = "inf"` |
| 4 | SuSiE vanilla | susine | mu = 0 | fixed sigma_0^2 | uniform | Verification baseline — should match #1 |
| 5 | SuSiE EB prior var | susieR | mu = 0 | EB estimated | uniform | `estimate_prior_variance = TRUE`; practitioner default |
| 6 | SuSiE EB prior var | susine | mu = 0 | EB estimated | uniform | Verification — should match #5 |
| 7 | SuSiE functional pi | susieR | mu = 0 | fixed sigma_0^2 | softmax(\|a_j\| / tau) | `prior_weights` from annotation |
| 8 | SuSiE functional pi | susine | mu = 0 | fixed sigma_0^2 | softmax(\|a_j\| / tau) | Verification — should match #7 |
| 9 | SuSiNE functional mu_0 | susine | mu_0 = c * a | fixed sigma_0^2 | uniform | Central method |

**Notes:**
- Use cases 4, 6, 8 exist to verify that the susine package matches susieR on equivalent specifications. If verification passes on a small subset, these can potentially be dropped from the full pilot to reduce compute.
- SuSiE-ash/inf (use cases 2, 3): susieR forces PIP-based convergence (not ELBO) when `unmappable_effects != "none"`. Use default convergence settings; no special handling needed.

### 2.2 Exploration Methods

All exploration methods run at equal budget K (number of fits).

**Pilot default: K = 12** (per study plan Section 5). Adjustable in run control workbook.

| Method | Applies to | Description |
|--------|-----------|-------------|
| Random starts | All non-functional prior use cases (1-6) | 1 default initialization + (K-1) Dirichlet warm starts. Each fit labeled as "default" or "warm". |
| c-grid | SuSiNE (#9) | K values of c for scaling prior means (mu_0 = c * a). Pilot default: c in {0, 0.1, 0.2, ..., 1.1}. Implemented via `annotation_scales` in `make_job_config()`. |
| tau-grid | Functional pi use cases (#7, #8) | K values of tau for softmax temperature on \|a_j\|. Pilot default: TBD empirically — needs sensible range where low tau concentrates pi on high-annotation SNPs and high tau approaches uniform. |

**Default sigma_0^2 for fixed-variance use cases:** `scaled_prior_variance = 0.2` (susieR default, meaning `0.2 * var(y)`). Implemented via `sigma_0_2_scalars` in `make_job_config()`.

### 2.3 Aggregation Methods

Applied post-hoc to each (dataset, use_case, exploration) group, **and** to the overall pool of all fits per dataset.

| Method | Description |
|--------|-------------|
| Max ELBO | Select single fit with highest ELBO |
| Uniform | Equal-weight average of all PIPs |
| ELBO softmax | w_i proportional to exp(ELBO_i); flatten all fits first |
| Cluster-then-ELBO-softmax | JSD-based clustering, importance-corrected weights (Method A from pip_ensemble_methods.md) |

### 2.4 Genotype Matrices

- 5-10 matrices from the existing collection
- Selected by stratified sampling on M1 (LD complexity metric)

### 2.5 Phenotype Simulations

**Sparse architecture:**
- Grid over p_star (number of causal variants): e.g., {1, 3, 5}
- Grid over y_noise (noise fraction / residual variance): e.g., {0.5, 0.8, 0.95}
- Full interaction (Cartesian product) of p_star x y_noise

**Oligogenic architecture:** Three-component partition per SuSiE-2 paper (sparse + oligogenic + polygenic tiers, each with own heritability budget). See T3 for implementation details. Whether to include in pilot is a design decision — adds complexity but directly relevant to SuSiE-ash/inf comparison arms.

### 2.6 Annotation Simulations

Controlled by existing parameters in `simulate_priors()`:

- Grid over `annotation_r2` (PVE_causal): e.g., {0.0, 0.1, 0.3, 0.5}. Fraction of causal effect variance explained by the annotation on causal indices. `annotation_r2 = 0` means annotation is pure noise (important control for SuSiNE).
- Grid over `inflate_match` (causal-to-noncausal annotation variance ratio): e.g., {0.0, 0.5, 1.0}. Noncausal annotation variance as a multiple of theoretical causal annotation variance.
- Annotation generation:
  - **Causal indices:** real beta + gaussian noise, with noise calibrated by `annotation_r2`
  - **Noncausal indices:** gaussian noise, with variance set by `inflate_match` relative to causal annotation variance
- Interact `annotation_r2` x `inflate_match` grids
- **Note:** `gamma_shrink` is deprecated (was for functional prior variances, dropped per D5). Do not include in new run configurations.

### 2.7 Metrics to Collect

**Model metrics (per fit, truth-aware):**
- Confusion PIP extract (TP, FP, TN, FN at various thresholds)
- Power @ FDR = 0.1, 0.2, 0.5
- AUPRC (average precision)

**Credible set metrics (per fit, truth-aware):**
- Purity (minimum absolute correlation within CS)
- Power (fraction of causal variants covered by any CS)
- Coverage (fraction of CSs containing a causal variant)
- Size (number of variants per CS)

**Dataset metrics (per dataset, truth-agnostic):**
- Z-score metrics: z_topk_ratio, z_max_abs, z_count_abs_gt_3, z_eff_signals
- LD metrics: M1, count of |r| > 0.95 off-diagonal (expected singleton bandits proxy), related transformations

**Multimodal metrics (per multi-fit group, truth-agnostic):**
- mean_jsd, median_jsd, max_jsd (pairwise Jensen-Shannon divergence)
- jaccard_top10 (top-10 PIP overlap)
- mean_pip_var (per-SNP PIP variance across fits)
- n_clusters (hierarchical clustering at JSD threshold)

---

## Section 3: Repository Update TODOs

### 3.1 susine Package Updates

| # | Task | Priority | Effort | Details |
|---|------|----------|--------|---------|
| S1 | Computational speedups | MEDIUM | Medium | 5 specific optimizations identified (see 3.1.1 below). ~40-50% wall-clock improvement expected. |
| S2 | Alpha/PIP convergence option | HIGH | Low | Add `convergence_method = c("elbo", "alpha")` parameter. Set default to match susieR 2.0's default. Track `max(abs(alpha_new - alpha_old))` per iteration. Files: `susine.R`, `susine_ss.R`. |
| S3 | LD regularization for SS/RSS | LOW (sim) / MED (real data) | Low-Med | Basic PSD projection + spectral shrinkage. Only needed for Phase D real data. Files: `susine_rss.R`, `initialize.R`. |
| S4 | Confirm baseline matching to susieR | HIGH | Low | Run susine and susieR on identical inputs (mu=0, fixed sigma_0^2, uniform pi, same init) and verify PIP/ELBO agreement. This validates use cases #1/#4, #5/#6, #7/#8 equivalence. |
| S5 | Verify non-uniform pi flows correctly | HIGH | Low | Test that `prior_inclusion_weights` parameter in susine is properly used in SER updates and not overwritten by defaults. Needed for functional pi use cases. |

#### 3.1.1 Computational Speedups (S1 detail)

| Sub-task | Description | Impact |
|----------|-------------|--------|
| S1a | Eliminate redundant `compute_Xty` in `BF()` — pass `Xty` from `SER()` | ~25% wall-clock (largest single win) |
| S1b | Compute ERSS once per effect, reuse for ELBO + variance update | Eliminates ~L redundant O(np) products |
| S1c | Eliminate redundant `Xb` in `SER_ERSS()` — pass from `SER()` | Eliminates ~L redundant O(np) products |
| S1d | Carry residual forward across outer iterations (incremental update) | Minor per-iteration savings |
| S1e | Audit `diag(X'X)` caching — verify no redundant column-norm recomputation | Likely already handled; verify |

### 3.2 test_susine Package Updates

| # | Task | Priority | Effort | Details |
|---|------|----------|--------|---------|
| T1 | Wire susieR use cases into harness | HIGH | Medium | The harness currently only calls `susine()`. Need a dispatch layer that calls `susieR::susie()` for use cases 1-3, 5, 7. **Key argument mapping:** susine's `prior_inclusion_weights` = susieR's `prior_weights`; susine takes raw `sigma_0_2` while susieR uses `scaled_prior_variance` (fraction of `var(y)`, default 0.2); susine's `init_alpha` = susieR's `model_init` (but format differs — susieR expects a full susie fit object); susine's `eb_prior_mode` maps to susieR's `estimate_prior_variance`. Read both packages' function signatures before implementing. |
| T2 | Refactor use case catalog | HIGH | Medium | Replace the 14-entry `use_case_catalog()` with a factored design. **New structure:** three independent lookup tables — `prior_spec_catalog()` (rows: susie, susie_fpi, susine, susie_ash, susie_inf, susie_eb, susine_eb, susie_fpi_susine, etc. — one per use case in Section 2.1), `exploration_catalog()` (rows: single, restart, c_grid, tau_grid, refine), `aggregation_catalog()` (rows: max_elbo, uniform, elbo_softmax, cluster_weight). Plus a `valid_exploration_for_prior()` constraint function (c_grid only for susine specs, tau_grid only for fpi specs, restarts for all). Run-control pipeline accepts these as independent arguments and crosses them, filtering by validity. Old 14-entry catalog is retired with no backward-compatibility shim. |
| T3 | Implement oligogenic architecture | MEDIUM | Medium | Per SuSiE-2 paper: three-tier effect partition (sparse: k_S large N(0, sigma_S^2) effects; oligogenic: k_O moderate effects from 2-component Gaussian mixture; polygenic: k_P small effects OR infinitesimal on all remaining variants). Each tier scaled to its h^2 target. New `simulate_effect_sizes_oligogenic()` or `architecture` flag. Annotation simulation should reflect tiered structure. |
| T4 | Implement refinement exploration mode | MEDIUM | Low-Med | BFS refinement in the test_susine harness — model-agnostic, budget-controlled, works with both `susine()` and `susieR::susie()` outputs. We do NOT use susieR's built-in `refine = TRUE`; all refinement is harness-level. Algorithm: from a fitted model, extract CSs, zero out CS SNPs in pi, refit, repeat as BFS tree up to budget K. Saves all fits for ensembling. ~80-120 lines. See BFS algorithm spec in `d1_package_architecture_analysis.md` Section A1. |
| T5 | Implement tau-grid exploration for functional pi | HIGH | Low | Analogous to existing c-grid (`annotation_scales`). Add `tau_grid` parameter to `make_job_config()` that expands functional-pi runs across tau values. |
| T6 | Implement cluster-then-ELBO-softmax aggregation (Method A) | HIGH | Low-Med | Clustering infra exists (JSD + hclust). Need: importance-corrected weights w_m ~ exp(ELBO_m) / f_m, representative selection per cluster, ESS diagnostic. |
| T7 | Label default vs warm start fits | LOW | Low | Tag each fit in the output with `init_type = "default"` or `"warm"`. Ensure the first fit in each restart group uses default initialization (no init_alpha). |
| T8 | Pass warm start K through run control workbook | HIGH | Low | K (budget = number of fits per exploration group) should be a top-level parameter in the run control workbook, flowing through `make_job_config()`. |
| T9 | Deprecate `gamma_shrink` in `simulate_priors()` | LOW | Low | `gamma_shrink` was used to construct functional prior *variances* from annotations (old b_i/b_iv use cases, dropped per D5). Add deprecation warning; remove from new run configurations. Existing `annotation_r2` and `inflate_match` already cover the notes' PVE_causal and causal-to-noncausal variance ratio specs. |
| T10 | Add overall per-dataset aggregation option | MEDIUM | Low | Pool all fits across use cases for a given dataset and apply aggregation methods. Wire as a run control option. |
| T11 | Allow exploration composability | MEDIUM | Medium | Each (exploration, prior_spec) pair is independently valid. Allow interacting explorations (e.g., 5 random starts x 3 c-values = 15 fits). Current grid expansion partially supports this but not fully exposed. |
| T12 | Separate exploration/aggregation control from use cases | HIGH | Medium | Use cases define model specs only. Explorations and aggregations are specified separately in the run control workbook. Related to T2. |
| T13 | Verify all metrics are collected and saved | HIGH | Low | Audit that model metrics, CS metrics, dataset metrics, and multimodal metrics from Section 2.7 are all computed and written to output during HPC runs. |
| T14 | Add high-LD-count metric | LOW | Low | Count of \|r\| > 0.95 off-diagonal in LD matrix (expected singleton bandits). Consider transformation to "expected number of singleton bandits per column." |
| T15 | Ensure per-fit PIP vectors are saved | HIGH | Low | Needed for post-hoc aggregation flexibility. Verify that individual fit PIPs and ELBOs are written, not just within-group aggregates. |
| T16 | Wire softmax temperature through job config to aggregation | MEDIUM | Low | ELBO softmax temperature exists in susine's `aggregate_susine_fits()` but is NOT exposed in test_susine job config. Need to thread it through `make_job_config()` and the aggregation call path. |
| T17 | Add wall time recording per fit | HIGH | Low | Record wall time per model fit in HPC outputs (per D8). Needed to verify per-fit cost is comparable across methods and for compute budgeting. |
| T18 | Seed management review | HIGH | Low | Before pilot HPC runs: verify reproducibility, confirm seed propagation through the full pipeline (phenotype generation, initialization, restart seeds, annotation draws). Set number of replicates. Per D16. |
| T19 | Consider complete vs average linkage for JSD clustering | LOW | Low | Multimodality metrics currently use average linkage in `compute_multimodal_metrics()`. `pip_ensemble_methods.md` recommends complete linkage for Method A (ensures all fits within a cluster are mutually close). Evaluate whether to switch. |
| T20 | Run 2-locus pipeline smoke test | HIGH | Low | Before pilot: run the full pipeline end-to-end on 2 matrices with minimal grid settings to catch integration bugs. Not a benchmarking run — just validation that all pieces connect. |

### 3.3 Dependency: Implementation Order

```
Phase 0 (prerequisite):
  S4 (confirm susine matches susieR) — gates all verification use cases
  S2 (alpha convergence) — for cross-method comparability
  S5 (verify non-uniform pi) — gates functional pi use cases

Phase 1 (core harness):
  T1 + T2 + T12 — ONE COUPLED REFACTOR. Do not implement independently.
  T5 (tau-grid) + T8 (K through run control) + T16 (softmax temp) + T17 (wall time) — wiring
  T6 (Method A aggregation) + T19 (linkage consideration) — aggregation
  T13 (verify metrics) + T15 (per-fit PIPs)

Phase 2 (enhancements):
  T3 (oligogenic) — meaningful new simulation function
  T4 (refinement) — nice to have for pilot
  T7 (label fits) — quality of life
  T9 (deprecate gamma_shrink) — cleanup
  T10 (overall aggregation) + T11 (composability)
  S1 (speedups) — reduces wall-clock for full study

Phase 3 (pre-pilot validation):
  T18 (seed management review) — must pass before submitting HPC jobs
  T20 (2-locus smoke test) — must pass before pilot
  T14 (LD count metric) — low priority
  S3 (LD regularization) — only for Phase D
```

---

## Appendix: Resolved Decisions (Summary)

All 17 critical decisions (D1-D17) from analysis_completion_status.md are resolved. Key outcomes:

| Decision | Resolution |
|----------|-----------|
| D1 Package architecture | Two-package: upstream susieR 2.0 + susine |
| D2 Use case refactoring | Replace old catalog with prior_spec x exploration x aggregation |
| D3 SuSiE-ash/inf | Include as baseline comparison arms (free in susieR 2.0) |
| D4 Refinement | Include as exploration competitor |
| D5 Functional pi | softmax(\|a_j\| / tau) with tau-grid exploration |
| D6 Aggregation design | Flatten-based, grid-agnostic. Four methods: max ELBO, uniform, ELBO softmax, cluster-then-weight |
| D7 Global ensembles | Yes, pool across model specs per dataset |
| D8 Compute budget | Fixed K (number of fits); record wall time per fit |
| D9 Z-score metrics | Current four candidates sufficient |
| D10 Pilot scope | Benchmarking exploration + aggregation on 5-10 datasets |
| D11 Real-data Phase D | Plan now, execute after simulation insights |
| D12 Annealing | Dropped from paper |
| D13 EB variants | Keep EB prior variance as practitioner-default baseline |
| D14 Annotation error model | Current controls sufficient, no sign flips |
| D15 Softmax temperature | tau-grid is exploration knob for functional-pi specs |
| D16 Seeds | Seed management review before pilot HPC runs |
| D17 Grid values | Defer to run control workbook |
