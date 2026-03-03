# SuSiNE Project Status

**Last updated:** 2026-03-02
**Supersedes:** `analysis_completion_status.md` (which remains as archival reference for resolved decisions D1-D17)

---

## Planning Artifacts

- **2026-02-27 detailed refactor scope (planning only):** `test_susine_refactor_scope_plan.md`
- Notes: This artifact expands T1-T20 into file-level scope, dependency constraints, and recommended execution order. No code implementation changes are included in that plan.

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
- mean_jsd, median_jsd, max_jsd (pairwise Jensen-Shannon style divergence on raw model-wide PIP vectors; no /L normalization in current code)
- jaccard_top10 (top-10 PIP overlap)
- mean_pip_var (per-SNP PIP variance across fits)
- n_clusters (hierarchical clustering at JSD threshold; threshold is interpreted on the same raw JSD scale above)

---

## Section 3: Repository Update TODOs

### 3.1 susine Package Updates

| # | Task | Priority | Effort | Status | Details |
|---|------|----------|--------|--------|---------|
| S1 | Computational speedups | MEDIUM | Medium | **COMPLETED (2026-02-27)** | Implemented S1a-S1d in `susine` core paths (cached `Xty` into `BF`, cached `Xb` into `SER_ERSS`, single ERSS pass per iteration for ELBO+sigma updates, residual carry-forward). |
| S2 | Alpha/PIP convergence option | HIGH | Low | **COMPLETED (2026-02-27)** | Added `convergence_method = c("elbo","alpha")` to `susine`, `susine_ss`, and `susine_rss`; tracks `model_fit$alpha_diff` and supports alpha-delta stopping. |
| S3 | LD regularization for SS/RSS | LOW (sim) / MED (real data) | Low-Med | **COMPLETED (2026-02-27)** | Added `regularize_ld_matrix()` (PSD projection + shrinkage) and wired RSS options: `ld_regularization`, `ld_shrinkage`, `ld_min_eig`. |
| S4 | Confirm baseline matching to susieR | HIGH | Low | **COMPLETED (2026-02-27)** | Baseline check script passes with fixed prior variance 0.2, EB prior updates off, residual variance estimated in both engines. |
| S5 | Verify non-uniform pi flows correctly | HIGH | Low | **COMPLETED (2026-02-27)** | Added explicit test confirming `prior_inclusion_weights` propagate through SER alpha updates (no silent overwrite by uniform defaults). |

#### 3.1.1 Computational Speedups (S1 detail)

| Sub-task | Description | Impact |
|----------|-------------|--------|
| S1a | Eliminate redundant `compute_Xty` in `BF()` - pass `Xty` from `SER()` | ~25% wall-clock (largest single win) |
| S1b | Compute ERSS once per effect, reuse for ELBO + variance update | Eliminates ~L redundant O(np) products |
| S1c | Eliminate redundant `Xb` in `SER_ERSS()` - pass from `SER()` | Eliminates ~L redundant O(np) products |
| S1d | Carry residual forward across outer iterations (incremental update) | Minor per-iteration savings |
| S1e | Audit `diag(X'X)` caching - verify no redundant column-norm recomputation | Verified existing `diag(X'X)` handling in SS/RSS paths; no additional change needed. |

#### 3.1.2 Completion Log (2026-02-27)

- Implemented code changes in `susine/R`: `SER.R`, `ERSS.R`, `ERSS_ss.R`, `susine.R`, `susine_ss.R`, `susine_rss.R`, `initialize.R`, `finalize.R`.
- Added focused status tests: `susine/tests/testthat/test-susine-status-updates.R`.
- Validation run: `devtools::test("susine")` -> **PASS 27, FAIL 0**.
- Baseline equivalence rerun: `test_susine/vignettes/one_off_validations/test_susine_susieR_baseline_match.R` -> **6/6 checks passed**; PIP corr `0.99983733`, max |PIP diff| `0.008644`, relative ELBO diff `0.000147`.
### 3.2 test_susine Package Updates

| # | Task | Priority | Effort | Status | Details |
|---|------|----------|--------|--------|---------|
| T1 | Wire susieR use cases into harness | HIGH | Medium | COMPLETED (2026-02-27) | Added backend dispatch in run_use_case() for susine and susieR; mapped fixed/EB prior variance settings; added compatibility guard for older susieR formals. |
| T2 | Refactor use case catalog | HIGH | Medium | COMPLETED (2026-02-27) | Added prior_spec_catalog(), exploration_catalog(), aggregation_catalog(), and valid_exploration_for_prior(); run controls now use factored catalogs. |
| T3 | Implement oligogenic architecture | MEDIUM | Medium | COMPLETED (2026-03-02) | Added `simulate_effect_sizes_oligogenic()` and wired `architecture_grid`/`architecture` through run controls and simulation data generation. |
| T4 | Implement refinement exploration mode | MEDIUM | Low-Med | COMPLETED (2026-03-02) | Implemented harness-level BFS refinement exploration with per-step prior-masking and `refine_step` run-table expansion. |
| T5 | Implement tau-grid exploration for functional pi | HIGH | Low | COMPLETED (2026-02-27) | Added tau_grid_values to make_job_config() / run-table expansion and validity-constrained crossing. |
| T6 | Implement cluster-then-ELBO-softmax aggregation (Method A) | HIGH | Low-Med | COMPLETED (2026-02-27) | Implemented cluster_weight aggregation with JSD clustering, representative selection, importance correction, and ESS output. |
| T7 | Label default vs warm start fits | LOW | Low | COMPLETED (2026-02-27) | Added init_type propagation (default/warm) through run table and output artifacts. |
| T8 | Pass warm start K through run control workbook | HIGH | Low | COMPLETED (2026-02-27) | K now drives exploration group expansion semantics in separate/intersect modes and is threaded in workbook usage. |
| T9 | Deprecate gamma_shrink in simulate_priors() | LOW | Low | COMPLETED (2026-02-27) | gamma_shrink now emits deprecation warning and is ignored in simulation generation path. |
| T10 | Add overall per-dataset aggregation option | MEDIUM | Low | COMPLETED (2026-02-27) | Added global per-dataset pooling with configurable overall_aggregation_methods and include_overall_pool. |
| T11 | Allow exploration composability | MEDIUM | Medium | COMPLETED (2026-02-27) | Exploration axes now compose by dataset x use_case x exploration group, with separate/intersect expansion and validity filtering. |
| T12 | Separate exploration/aggregation control from use cases | HIGH | Medium | COMPLETED (2026-02-27) | Use cases now define prior specs only; explorations and aggregations are independent run-control inputs. |
| T13 | Verify all metrics are collected and saved | HIGH | Low | COMPLETED (2026-03-02) | Added `validate_metrics_coverage()` automated schema/file checks and validated model/effect/dataset/confusion/multimodal/validation/snp outputs in strict smoke tests. |
| T14 | Add high-LD-count metric | LOW | Low | COMPLETED (2026-02-27) | Added high_ld_count_095 and high_ld_count_095_per_snp dataset metrics. |
| T15 | Ensure per-fit PIP vectors are saved | HIGH | Low | COMPLETED (2026-02-27) | Per-fit SNP/PIP parquet persistence and downstream aggregated snps_dataset generation verified in buffered aggregation smoke path. |
| T16 | Wire softmax temperature through job config to aggregation | MEDIUM | Low | COMPLETED (2026-02-27) | Added softmax_temperature to job config and wired through elbo_softmax/cluster_weight aggregation calls. |
| T17 | Add wall time recording per fit | HIGH | Low | COMPLETED (2026-02-27) | Added per-fit wall-clock timing (wall_time_sec) into fit metadata and output tables. |
| T18 | Seed management review | HIGH | Low | COMPLETED (2026-03-02) | Added `seed_management_report()` and validated seed propagation/uniqueness checks in the two-locus smoke workflow. |
| T19 | Consider complete vs average linkage for JSD clustering | LOW | Low | COMPLETED (2026-03-02) | Switched multimodal cluster counting to complete linkage for consistency with Method A clustering assumptions. |
| T20 | Run 2-locus pipeline smoke test | HIGH | Low | COMPLETED (2026-03-02) | Added `run_two_locus_smoke_test()` and executed strict pass (2 bundles, staged aggregation, metrics/seed checks). |

#### 3.2.1 Functional Validation Snapshot (2026-03-02)

- Restart exploration smoke check: K=4, restart produced expected split default=1, warm=3, and run_task() completed.
- Intersect exploration smoke check: exploration_methods = c("single","restart"), exploration_mode = "intersect", K=2 completed after fixing a list-column expansion bug in run controls.
- Buffered staging + aggregation smoke check: with verbose_file_output = FALSE, index_staging_outputs() indexed 6 files and aggregate_staging_outputs() wrote expected outputs (model_metrics.csv, effect_metrics.csv, validation.csv, dataset_metrics.csv, confusion_bins.csv, snps_dataset/...).
- Refine exploration smoke check: exploration_methods = "refine", K=3 completed with BFS prior-masking flow and per-step outputs.
- Oligogenic architecture smoke check: architecture_grid = "oligogenic" completed end-to-end through run_task().
- Two-locus strict smoke check: `run_two_locus_smoke_test()` passed with required multimodal metrics and seed management checks.
- Known packaging blocker: devtools::check() is currently blocked by malformed package name in DESCRIPTION (pre-existing project issue).

#### 3.2.2 Remaining Work for test_susine

1. No open T3-T20 implementation items remain.
2. Pre-existing package metadata blocker remains: malformed package name in DESCRIPTION prevents `devtools::check()`.

### 3.3 Dependency: Implementation Order

Phase 0 (prerequisite) - COMPLETE:
  S4 (confirm susine matches susieR) - COMPLETE 2026-02-27
  S2 (alpha convergence) - COMPLETE 2026-02-27
  S5 (verify non-uniform pi) - COMPLETE 2026-02-27

Phase 1 (core harness) - COMPLETE:
  T1 + T2 + T12 - COMPLETE 2026-02-27
  T5 + T8 + T16 + T17 - COMPLETE 2026-02-27
  T6 - COMPLETE 2026-02-27
  T13 - COMPLETE 2026-03-02
  T15 - COMPLETE 2026-02-27

Phase 2 (enhancements) - COMPLETE:
  T7 - COMPLETE 2026-02-27
  T9 - COMPLETE 2026-02-27
  T10 + T11 - COMPLETE 2026-02-27
  T3 - COMPLETE 2026-03-02
  T4 - COMPLETE 2026-03-02
  T19 - COMPLETE 2026-03-02
  S1 (speedups) - COMPLETE 2026-02-27

Phase 3 (pre-pilot validation) - COMPLETE:
  T18 (seed management review) - COMPLETE 2026-03-02
  T20 (2-locus smoke test) - COMPLETE 2026-03-02
  T14 (LD count metric) - COMPLETE 2026-02-27
  S3 (LD regularization) - COMPLETE 2026-02-27 (available for Phase D)

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


