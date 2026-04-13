# Compute Allocation Experiment: Ensemble Fine-Mapping Scaling Laws

## Context

**Question:** Given a compute budget of K model fits instead of 1, what is the optimal way to allocate those fits across 4 levers for ensemble fine-mapping? How does the answer change as K scales from 10 to 100?

**4 levers:**

1. **Random restarts** - same model, different alpha initializations via `init_alpha` only (Dirichlet-random). First restart = baseline (default init).
2. **Refinements** - use the actual SuSiE-2 / `susieR` refinement procedure (`refine = TRUE` semantics): for each CS, zero prior weights on that CS and refit, then refit from the converged model under the original prior weights, keeping the best improvement. This is **not** `prior_refit`, and it is **not** the current harness BFS notion unless the harness is explicitly reworked to match upstream behavior.
3. **Annotation scale grid (c_grid)** - sweep c in [0, 1.5] for mu_0 = c*a, fixed sigma_0_2=0.2. Required whenever annotations are used without EB (no principled default c).
4. **Prior variance grid (sigma_0_2_grid)** - sweep sigma_0_2 in logspace(0.01, 1.0). Default when not gridded: 0.2. All sigma_0_2 grids hard-code 0.2 as one of the grid values. *(New exploration axis - needs implementation.)*

**Design:** Run full-sized grids for each lever and their interactions. Scaling curves come from subsampling these grids at smaller sizes post-hoc. Interaction grids use the Cartesian product of lever dimensions (e.g., 8 restart seeds x 8 refine steps = 64 actual runs where each run uses a specific restart seed AND refine depth).

---

## Architecture and Dataset Grid

**Architecture:** susie2_oligogenic (three-tier)

- 3 sparse (beta ~ N(0,1), 50% h2) + 5 oligogenic (mixture, 35% h2) + 15 polygenic (beta ~ N(0,0.08), 15% h2) = 23 causal variants, h2=0.25 fixed

| Parameter | Values | Count |
|-----------|--------|-------|
| Genotype matrices | Full 150-locus set | 150 |
| Seeds | 3 per matrix | 3 |
| Architecture | susie2_oligogenic | 1 |
| Annotation R-squared | 0 (null), 0.2 (moderate), 0.5 (strong) | 3 |

**Base datasets:** 150 x 3 = **450**

---

## Method Catalog

### Additional standalone baselines

Add **three additional standalone baseline runs** warm-started on truth:

1. `susine` with `init_alpha` warm start on the true top 10 causal variants
2. annotated EB `susine` with `init_alpha` warm start on the true top 10 causal variants
3. fixed-grid `susine` with `c = 0.5`, `sigma_0_2 = 0.2`, and `init_alpha` warm start on the true top 10 causal variants

Shared settings:

- one variant per effect
- no ensemble aggregation applied

These are standalone benchmarks / upper-context baselines, not part of the 4-lever scaling grids above.

For annotated standalone baselines, include one run per dataset per annotation R-squared level.

### Group A: Baseline SuSiE (no annotations, R-squared-invariant)

3 levers: restarts (R), refinements (F), sigma_0_2 grid (S).

| Spec | Grid | Runs | Purpose |
|------|------|------|---------|
| A-R | 64 restarts | 64 | Pure restart scaling |
| A-F | 32 refine steps x 2 | 64 | Pure refinement scaling |
| A-S | 64 sigma_0_2 values | 64 | Pure variance grid scaling |
| A-RF | 8R x 8F | 64 | Restart x refinement interaction |
| A-RS | 8R x 8S | 64 | Restart x variance interaction |
| A-FS | 8F x 8S | 64 | Refinement x variance interaction |
| A-RFS | 4R x 4F x 4S | 64 | Full 3-way interaction |
| **Total (gross)** | | **448** | |
| **Total (dedup)** | | **~364** | Nested grids share baseline + edge runs |

### Group B: EB SuSiNE (annotations via EB, per R-squared level)

EB estimates c and sigma_0_2 internally. 2 levers: restarts (R), refinements (F).

| Spec | Grid | Runs | Purpose |
|------|------|------|---------|
| B-R | 64 restarts | 64 | EB with restart diversity |
| B-F | 32 refine steps x 2 | 64 | EB with refinement search |
| B-RF | 8R x 8F | 64 | Restart x refinement interaction |
| **Total (gross)** | | **192** | |
| **Total (dedup)** | | **~176** | Shared baseline + r=0 refine chain |

### Group C: Grid SuSiNE (annotations, no EB, per R-squared level)

No default c, so c_grid is always required. 4 levers: c_grid (C), sigma_0_2 grid (S), refinements (F), restarts (R).

| Spec | Grid | Runs | Purpose |
|------|------|------|---------|
| C-C | 64 c values | 64 | Pure c_grid scaling |
| C-CS | 8C x 8S | 64 | Joint c x sigma_0_2 surface |
| C-CSR | 8C x 8S search + exact default-prior refit | 64 | c x sigma search followed by warm refit under default priors |
| C-CSFR | 3C x 3S x 3F x 3R | 81 | Full 4-way interaction (coarse) |
| C-CFR | 4C x 4F x 4R | 64 | c x refine x restart (sigma=0.2) |
| **Total (gross)** | | **273** | |
| **Total (dedup)** | | **~228** | Nested c grids share points at c boundaries |

---

## Compute Budget

| Group | Runs/dataset (gross) | Runs/dataset (dedup) | x Datasets | x R-sq levels | Gross | Dedup |
|-------|---------------------|---------------------|------------|---------------|-------|-------|
| A (baseline) | 448 | ~364 | 450 | 1 | 201,600 | ~163,800 |
| B (EB, x3 R-sq) | 192 | ~176 | 450 | 3 | 259,200 | ~237,600 |
| C (grid, x3 R-sq) | 273 | ~228 | 450 | 3 | 368,550 | ~307,800 |
| **Total** | | | | | **829,350** | **~709,200** |

**71% of 1M budget** after dedup. Margin for follow-up susie2_sparse experiment or reruns.

Deduplicated run counts are still important and should be tracked explicitly for HPC accounting. The key distinction is:

- **gross / requested runs** = the conceptual experiment design
- **deduplicated / realized executions** = the actual number of model calls after cache reuse

Both numbers should be reported.

---

## Grid Value Specifications

### c_grid values: seq(0, 1.5, length.out=N)

| N | Formula | Raw values | Rounded |
|---|---------|-----------|---------|
| 3 | seq(0, 1.5, 3), step=0.75 | 0, 0.75, 1.5 | 0, 0.75, 1.5 |
| 4 | seq(0, 1.5, 4), step=0.5 | 0, 0.5, 1.0, 1.5 | 0, 0.5, 1.0, 1.5 |
| 8 | seq(0, 1.5, 8), step=1.5/7=0.2143 | 0, 0.2143, 0.4286, 0.6429, 0.8571, 1.0714, 1.2857, 1.5 | 0, 0.21, 0.43, 0.64, 0.86, 1.07, 1.29, 1.5 |
| 64 | seq(0, 1.5, 64), step=1.5/63=0.0238 | 0, 0.0238, 0.0476, ..., 1.4762, 1.5 | 64 evenly-spaced values, step ~0.024 |

### sigma_0_2 grid values: logspace with 0.2 hard-coded

Raw formula: 10^seq(-2, 0, length.out=N). Problem: 0.2 = 10^(-0.699) does not land on an evenly-spaced log grid. Solution: use hand-picked grids for small N that include 0.2 as an exact value. For N=64, replace the nearest logspace point with 0.2.

| N | Raw logspace 10^seq(-2,0,N) | Proposed (with 0.2 ensured) |
|---|----|----|
| 3 | 0.01, 0.1, 1.0 | **0.05, 0.2, 1.0** |
| 4 | 0.01, 0.0464, 0.2154, 1.0 | **0.02, 0.1, 0.2, 0.5** |
| 8 | 0.01, 0.0193, 0.0373, 0.0720, 0.1389, 0.2682, 0.5179, 1.0 | **0.01, 0.03, 0.07, 0.1, 0.2, 0.4, 0.7, 1.0** |
| 64 | 0.01, 0.0112, ..., 0.8913, 1.0 (step in log10 = 2/63 = 0.0317) | 10^seq(-2, 0, 64) with nearest point to 0.2 snapped to exactly 0.2 (point 41 at 10^(-0.683)=0.2075 becomes 0.2) |

**Default when not gridded:** sigma_0_2 = 0.2

### Restart grid

N restart seeds. Index 0 = baseline (default init), indices 1..N-1 = Dirichlet-random `init_alpha` warm starts only (`alpha_concentration=0.001`).

### Refinement grid

N refinement ensemble size (including the baseline as member 1). The refinement procedure produces a **tree** of models via BFS:

**Per-node refinement procedure (2 model fits → 1 stored result):**

1. Take the parent node's converged fit. Examine its first L effects.
2. For each effect whose 95% CS has purity > 0.95:
   - **Perturb:** Zero prior inclusion weights for all SNPs in that CS. Set remaining weights uniform. Fit susine to convergence. (1 model fit, intermediate — not stored.)
   - **Refit:** Take the perturbed model's final alpha as `init_alpha`. Restore original (all-uniform) prior inclusion weights. Fit susine to convergence. (1 model fit — stored as a child node.)
3. Each eligible effect produces one child node in the tree.

**Tree construction (BFS to fill ensemble):**

- Root = baseline model (1 run, stored).
- Level 1: one child per eligible effect of the root. Each costs 2 fits.
- If more models needed: go to the first child node, repeat the procedure on its effects to produce grandchildren.
- Continue BFS until ensemble has N stored models (including root).

**Cost accounting:**

- N stored models = 1 baseline + (N-1) perturb-refit pairs = 1 + 2(N-1) = **2N-1 total model fits**.
- For "8 refine": 1 + 14 = 15 model fits, 8 stored models.
- For "32 refine": 1 + 62 = 63 model fits, 32 stored models.
- The "64 runs" in the plan tables counts model fits, not stored models. So "32 refine steps x 2 = 64" means 32 perturb-refit pairs = 64 fits + 1 shared baseline = 65 fits total, producing 33 stored models.

**Note on budget table reconciliation:** The plan table says A-F = 64 runs. This is 32 perturb-refit pairs (64 fits). The baseline is shared with other specs (dedup). So A-F produces 32 child models + 1 shared baseline = 33 ensemble members from 64 dedicated fits.

**Current harness status (T4 BFS in run_model.R:830-913):**

The existing BFS refinement already implements the tree/queue structure with `blocked_idx` and CS extraction. The key gap: **it only does the perturb step (blocked weights, single fit), not the two-step perturb + refit-with-restored-weights.** Specifically:
- Current: 1 fit per node (perturb with blocked weights, store that result)
- Target: 2 fits per node (perturb with blocked weights → take alpha → refit with uniform weights, store only the refit)

The refactoring is localized: inside the BFS loop (run_model.R:857-905), after the perturb fit, add a refit call with restored weights and init_alpha from the perturb.

### Interaction grids

Cartesian product of lever dimensions. E.g., 8R x 8F = 8 restart seeds, each refined to N=8 depth (7 perturb-refit pairs + unrefined base = 8 stored models per start chain). Each start gets its own refinement tree.

---

## Implementation Notes

- **One job config for everything:** this experiment must run from a single job config. That will require a rework of run-table generation so one config can contain heterogeneous specs with different exploration axes and axis sizes.
- **Deduplication logic:** test_susine already has run deduplication via `execution_cache_key()` in run_model.R (17-component key). This is useful and should be retained for HPC run accounting. Need to verify it correctly identifies shared runs across specs (e.g., baseline model shared between A-R, A-F, A-RF), and consolidated outputs should report both gross requested runs and deduplicated realized executions.
- **Aggregation methods to keep:** cluster_weight (unnormalized categorical JSD, 0.15 distance, complete linkage), cluster_weight_050 (same but 0.50 distance), elbo_softmax, uniform, max_elbo. Drop all other aggregation variants (multi-threshold bjsd variants, other jsd thresholds) from the harness to reduce output bloat.
- **Runtime tracking:** Verify that per-model wall-clock runtime is saved in model_metrics.csv. If not, add it. This is needed for the cost-aware scaling analysis (AUPRC per compute-second, not just per model-fit).
- **Oracle handling:** do not add oracle rows to the run table. Oracle should be added only in consolidated post-hoc outputs (at minimum, confusion-bin-derived outputs) with `agg_method = "oracle"`.

---

## Post-Hoc Scaling Analysis (No Additional Model Fits)

### Analysis 1: Single-Lever Scaling Curves (subsampling 64x grids)

For each pure 64-run spec (A-R, A-F, A-S, B-R, B-F, C-C), subsample to ensemble sizes n in {4, 8, 16, 32, 64}:

- **Restarts:** 50 random subsets of size n, aggregate PIPs, evaluate AUPRC. Report mean +/- SE.
- **Refinements:** prefix of first n/2 BFS steps (= n runs).
- **Grids** (c, sigma_0_2): take n evenly-spaced grid points from the 64 (structured subsampling preserving coverage).

**Primary metric:** AUPRC. **Secondary:** Power at FDR <= 0.05, CS coverage, mean CS size.

### Analysis 2: Interaction Scaling (subsampling 8x8 grids)

For each 8x8 interaction spec, subsample to 4x4 (take every other point on each axis):

- A-RF (8x8 -> 4x4): restart-refine interaction at half resolution
- A-RS (8x8 -> 4x4): restart-sigma interaction at half resolution
- A-FS (8x8 -> 4x4): refine-sigma interaction at half resolution
- C-CS (8x8 -> 4x4): c-sigma joint surface at half resolution

### Analysis 3: Compare Allocation Strategies at Fixed Budget

At a fixed budget B (e.g., B=64), compare:

- Pure lever: all 64 on restarts vs all 64 on c_grid vs all 64 on refinements
- Interaction: 8x8 restart-refine vs 8x8 c-sigma vs 4x4x4 three-way
- Key question: at budget 64, is it better to do 64 restarts, 8x8 restart-refine, or 4x4x4 three-way?

### Analysis 4: Aggregation Method Comparison

For each ensemble, compare max_elbo, uniform, elbo_softmax, cluster_weight.

Key hypothesis: cluster_weight dominates for multimodal ensembles (restarts, refinements) where PIPs concentrate on distinct basins. max_elbo or elbo_softmax may suffice for grid searches that explore a smooth hyperparameter surface.

### Analysis 5: Annotation Quality Stratification

Stratify all results by R-squared:

- **R-squared = 0:** Null annotations. EB and grid methods should degrade gracefully (EB finds c~0, c_grid includes c=0). Establishes that annotations don't hurt when uninformative.
- **R-squared = 0.2:** Moderate. Most realistic. Where we expect the biggest divergence between strategies.
- **R-squared = 0.5:** Strong annotations. EB likely dominates since it can adapt c precisely. Grid search still useful but c_grid quickly converges.

### Final Recommendation Output

The scaling analysis produces a **recommendation table**: for each (annotation quality, budget range), the optimal allocation strategy. E.g.:

- Budget 1-4: single fit with EB (if annotations available) or vanilla
- Budget 4-16: restarts dominate (search landscape exploration)
- Budget 16-64: interaction grids (e.g., 4x4x4 c-refine-restart) provide best marginal returns
- Budget 64+: diminishing returns, concentrate on the strongest lever

---

## Setup and Settings

### New Workbooks

Create three clean new workbooks for this experiment (NOT reusing existing pilot workbooks):

1. **run_control workbook** - job config construction, run_table generation, SLURM submission
2. **collect_results workbook** - post-hoc aggregation, subsampling analysis, metric computation
3. **visualize_results workbook** - all figure generation (see Figures section below)

These can edit susine/ and test_susine/ R packages as needed but the workbook code stays separate.

### PIP Breaks

Set PIP bucket breaks to variable resolution:

- 0.02 steps from 1.00 down to 0.10 (46 breaks)
- 0.01 steps from 0.10 down to 0.05 (5 breaks)
- 0.005 steps from 0.05 down to 0.00 (10 breaks)

### Other Settings

- **inflate_match = 1** for all annotations (noncausal annotation variance = causal annotation variance)
- **Aggregation methods (5 total):** cluster_weight (JSD 0.15), cluster_weight_050 (JSD 0.50), elbo_softmax, uniform, max_elbo
- **Oracle aggregation:** do not implement as a scheduled run-table method. Add it only during consolidation / post-hoc analysis as a truth-aware ceiling.
- **Additional standalone baselines:** include all three truth-warm-start standalone baselines in any figure/table where the corresponding ordinary single-fit baseline is shown.

### Grading Convention

For susie2_oligogenic (23 causal), grade using **top 8 as causal** (3 sparse + 5 oligogenic effects). Filter to the 8 largest absolute effect sizes per dataset. The 15 polygenic effects are too small to detect and should not count against FDR/power.

---

## Figures

### Category 1: SuSiE2-Style Plots (individual model metrics, max-ELBO selection)

These use credible sets and effect sizes, which are only available from individual model fits (not ensemble-aggregated PIPs). For each ensemble, select the single best model by max-ELBO.

**1a. CS Power vs. top-n-as-causal**

- Y: CS power (fraction of true causal variants covered by a CS)
- X: top n variants treated as causal (sweep n)
- Color by ensemble method
- Only full ensembles (no subsampled versions)

**1b. CS FDR vs. top-n-as-causal**

- Y: CS FDR (fraction of CSs not containing a true causal)
- X: top n as causal
- Color by ensemble method, full ensembles only

**1c. Pearson R-squared (effect size prediction)**

- Y: R-squared between estimated and true effect sizes
- Color by ensemble method, full ensembles only

**1d. TPR vs FPR curve (FPR 0 to 0.1)**

- Grade using top 8 as causal
- Uses PIPs, so can use all aggregation methods
- Facet by aggregation method (including oracle), color by ensemble method

**1e. PIP calibration plots**

- Bin PIPs in 0.10-width bins (10 bins: [0,0.1), [0.1,0.2), ..., [0.9,1.0])
- Plot observed fraction causal vs. mean PIP per bin
- Facet by aggregation method, color by ensemble method

### Category 2: Ensemble Performance Plots (PIP-based, all aggregation methods)

**2a. Precision-Recall curves**

- Grade using top 8 as causal
- Facet by aggregation method (including oracle), color by ensemble method

**2b. Multimodality metrics boxplots**

- One PNG, faceted by multimodality metric (jaccard_top10, mean_jsd, median_jsd, max_jsd, mean_pip_var, n_clusters)
- Each plot: boxplots colored by ensemble method

**2c. AUPRC bar chart by method (cf. warm_start_pilot_dual auprc_barplot_by_method.png)**

- X-axis groupings: ensemble method
- Color by aggregation method (including oracle)
- Horizontal red dashed line: baseline susine model at sigma_0_2=0.2 (single fit, no ensemble)
- Also include the standalone truth-warm-start baselines anywhere the corresponding baseline comparison is shown
- This is the primary "which ensemble strategy wins?" figure

### Category 3: Dataset Difficulty Plots

**3a. Dataset metrics histograms**

- Histograms of M1, z-based metrics (z_topk_ratio, z_max_abs, z_count_abs_gt_3, z_eff_signals), LD-based metrics (high_ld_count_095)
- Correlation matrix of all dataset metrics

**3b. Dataset difficulty scatter (M1 as primary predictor)**

- Take baseline SuSiE only (sigma_0_2=0.2, single fit) as the reference performance
- Scatter: AUPRC (y) vs M1 (x), one point per dataset
- Find one z-based metric + one LD-based metric that, combined with M1, maximizes explained AUPRC variation (via multiple regression, report R-squared)
- Additional scatter plots: same AUPRC (y) vs the two supplementary metrics (faceted horizontally)
- Regression lines from the combined 3-predictor model overlaid on all three scatter panels
- Then overlay points (different colors) for:
  - Best no-annotation ensemble (single best ensemble method x aggregation combo)
  - Best with-annotation ensemble (single best combo across all annotation ensembles)

### Category 4: Ensemble Scaling Plots

**4a. 64-grid scaling: AUPRC**

- Y: AUPRC, X: ensemble size {4, 8, 16, 32, 64}
- Color by aggregation method, facet by ensemble method (one panel per spec: A-R, A-F, A-S, B-R, B-F, C-C)

**4b. 64-grid scaling: TPR at FPR=0.05**

- Same as 4a but Y = TPR at fixed FPR=0.05

**4c. 8x8 interaction scaling: AUPRC**

- Y: AUPRC, X: grid resolution {4x4, 8x8}
- Color by aggregation method, facet by interaction spec (A-RF, A-RS, A-FS, C-CS)

**4d. 8x8 interaction scaling: TPR at FPR=0.05**

- Same as 4c but Y = TPR at fixed FPR=0.05

### Figure Output Organization

All figures saved under `slurm_output/{job_name}/{job_id}/aggregate/figures/` with subfolders:

```text
figures/
  susie2_plots/          # Category 1: CS power, CS FDR, R2, TPR-FPR, PIP calibration
  ensemble_performance/  # Category 2: PR curves, multimodality, AUPRC bars
  dataset_difficulty/    # Category 3: histograms, correlations, M1 scatters
  ensemble_scaling/      # Category 4: 64-scaling, 8x8-scaling
```

---

## Implementation Plan

### Step 1: Implement sigma_0_2_grid exploration axis

**Pattern:** Mirror c_grid implementation.

**Files to modify:**

1. **test_susine/R/use_cases.R**
   - Add "sigma_0_2_grid" row to exploration_catalog()
   - Update valid_exploration_for_prior(): valid for fixed-variance prior specs (not EB variance)

2. **test_susine/R/run_controls.R**
   - Add sigma_0_2_grid_values param to make_job_config(), default 10^seq(-2, 0, length.out=K) with 0.2 ensured
   - In build_exploration_groups(): expand K rows with sigma_0_2_scalar column
   - Add to make_default_grid()

3. **test_susine/R/run_model.R**
   - In run_use_case(): when sigma_0_2_grid exploration, use per-run sigma_0_2_scalar (column already exists in run_table)

### Step 2: Implement interaction grid support

Need to support Cartesian products of exploration axes (e.g., restart x refine, c_grid x sigma_0_2).

**Check:** Does exploration_mode="intersect" in build_exploration_groups() already support this? From the codebase exploration, intersect mode uses `tidyr::crossing` to create Cartesian products. If so, just need to wire it up for the new sigma_0_2 axis.

**But:** one job config for the whole experiment will require more than that. The run-table builder needs to support heterogeneous per-spec axis definitions and sizes within a single config, rather than one global `K` / `exploration_methods` setting for the whole job.

**Files to modify:**

- test_susine/R/run_controls.R - extend or redesign run-table generation for multi-axis Cartesian products and heterogeneous per-spec grids in one job config
- test_susine/R/run_model.R - ensure the harness can read and use the richer run table cleanly

### Step 3: Verify infrastructure

- annotation_r2 gridding works with {0, 0.2, 0.5}
- architecture_grid = "susie2_oligogenic" works (h2-based, no y_noise/p_star)
- actual refinement dispatch matches upstream `susieR` semantics where intended
- execution cache (`execution_cache_key` in run_model.R) correctly deduplicates shared runs across specs

### Step 4: Build one job config

The goal is **one job config for everything**, not one per group. This will require the run table to support heterogeneous specs such as A-R, A-F, A-RF, and C-CSFR in the same job artifact.

This same job config should also include all three standalone truth-warm-start baseline row types.

### Step 5: Write post-hoc scaling analysis script

New test_susine/analysis/scaling_analysis.R:

1. Read per-run PIP vectors from snps_dataset/ parquet output
2. Group by (dataset, spec, annotation_r2)
3. For pure-lever specs: subsample at {4, 8, 16, 32, 64}
4. For interaction specs: subsample by reducing dimension sizes (8x8 -> 4x4)
5. For each subsample: aggregate PIPs via all 5 methods, evaluate AUPRC + secondary metrics
6. For fixed-budget comparison: line up all specs that cost the same budget
7. Output: scaling_curves.csv, interaction_scaling.csv, budget_comparison.csv, aggregation_comparison.csv

### Step 6: Add oracle only in consolidated outputs

Do **not** add oracle to the run table. Do **not** schedule oracle model fits.

Instead, during consolidation / analysis, add `agg_method = "oracle"` rows to the consolidated outputs by selecting the single best model within each ensemble group using truth-aware AUPRC. This is a ceiling and should live only in post-hoc outputs.

### Step 7: Figure generation

See the detailed Figures section above. Implement in the visualize_results workbook.

---

## Verification

1. **Smoke test:** 2 matrices x 1 seed x 1 R-sq x multiple spec types in a single job config. Verify PIP storage and subsample analysis end-to-end.
2. **Single-config check:** confirm the one job config can contain heterogeneous specs with different exploration axes / sizes.
3. **Budget check:** `nrow(run_table)` matches gross expectation, and deduplicated realized executions match the expected cache-reuse savings.
4. **sigma_0_2 grid includes 0.2:** verify 0.2 is always one of the grid values at any grid size.
5. **Refinement semantics:** verify the refinement path used for this study matches actual upstream SuSiE-2 refinement semantics where intended.
6. **Oracle discipline:** verify oracle rows appear only in consolidated post-hoc outputs, never in the run table.
7. **R-sq=0 sanity:** with null annotations, EB and grid methods should not outperform vanilla due to an implementation artifact.
