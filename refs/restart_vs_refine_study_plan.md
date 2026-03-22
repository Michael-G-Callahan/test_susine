# Restart vs. Refinement Tradeoff Simulation Study

## Context

SuSiE's posterior landscape is multimodal under LD. Given a fixed compute budget of K model fits per locus, how should we allocate between random restarts (global exploration via Dirichlet-sampled pi) and refinements (local exploration via BFS credible-set zeroing)? This study answers that question empirically, then transfers the finding to the annotated SuSiNE case.

**Key insight**: every "fit" (restart or refine) is a full IBSS run — equal compute cost. The only question is the allocation ratio.

---

## 1. Study Design

### 1.1 Datasets

- **150 X matrices** from `sampled_simulated_genotypes/`, **scenario_1** only (n=600)
- **Architecture**: oligogenic only (3-tier: sparse/oligo/poly; p_star and y_noise fixed by architecture)
- **Seeds**: 3 phenotype seeds per matrix (distinct; seeds change which variants fall into which significance tiers)
- **Total**: 150 x 3 = **450 datasets**
- **L = 10** for all fits

### 1.2 Prior Configurations (4 total)

| Config  | sigma_0_2  | prior_update_method | Label     |
|---------|------------|---------------------|-----------|
| Fixed 1 | 0.1        | "none"              | sigma_0.1 |
| Fixed 2 | 0.2        | "none"              | sigma_0.2 |
| Fixed 3 | 0.4        | "none"              | sigma_0.4 |
| EB      | NA (learned) | "var"             | eb_var    |

All use susine backend with zero prior mean (no annotations) and uniform baseline pi.

### 1.3 Arms (6 total, each = 64 fits per prior config)

| Arm | R starts | F refines/start | Total fits | Start % |
|-----|----------|-----------------|------------|---------|
| 1   | 64       | 0               | 64         | 100%    |
| 2   | 32       | 1               | 64         | 50%     |
| 3   | 16       | 3               | 64         | 25%     |
| 4   | 8        | 7               | 64         | 12.5%   |
| 5   | 4        | 15              | 64         | 6.25%   |
| 6   | 1        | 63              | 64         | 1.6%    |

**Start 1 is always the baseline** (uniform pi). Starts 2+ use Dirichlet(alpha=0.01)-sampled pi.

### 1.4 Refinement Mechanics

- BFS order: from a converged parent, refine each purity-filtered CS (zero out CS variants' pi, re-fit from scratch)
- Depth 1 = refine each CS of the parent (up to L=10 children)
- Depth 2+ = refine CSs of depth-1 fits, etc.
- "F refines per start" = F BFS steps consumed from the queue for that start's tree

---

## 2. Smart Reuse Strategy

### 2.1 Key Observation

Arms share starts (same deterministic seeds: `phenotype_seed + restart_id`). Start k in arm 2 produces the identical fit as start k in arm 4. Similarly, the BFS refinement tree from start k is deterministic — refine step 3 from start k is the same fit regardless of which arm requested it.

### 2.2 Maximum Refine Depth Per Start (Union Across Arms)

| Start IDs | Max refines needed | Source arm(s) |
|-----------|--------------------|---------------|
| 1         | 63                 | arm 6         |
| 2-4       | 15                 | arm 5         |
| 5-8       | 7                  | arm 4         |
| 9-16      | 3                  | arm 3         |
| 17-32     | 1                  | arm 2         |
| 33-64     | 0                  | arm 1 only    |

### 2.3 Unique Fits Per Prior Config

| Component                          | Count |
|------------------------------------|-------|
| 64 starts (root fits)              | 64    |
| Start 1: 63 refine steps          | 63    |
| Starts 2-4: 15 refine steps each  | 45    |
| Starts 5-8: 7 refine steps each   | 28    |
| Starts 9-16: 3 refine steps each  | 24    |
| Starts 17-32: 1 refine step each  | 16    |
| **Total unique fits per prior config** | **240** |

### 2.4 Compute Budget Comparison

|                                      | Fits/prior | x 4 configs | x 450 datasets | Total       |
|--------------------------------------|-----------|-------------|----------------|-------------|
| **Without reuse** (6 separate arms)  | 384       | 1,536       | 691,200        | **691,200** |
| **With reuse** (single super-job)    | 240       | 960         | 432,000        | **432,000** |
| **Savings**                          | 144 (37.5%) | 576       | 259,200        | **259,200 fits saved** |

### 2.5 Arm Reconstruction (Post-Hoc Filtering)

Each arm is a filter on (restart_id, bfs_step) from the 240-fit pool:

- **Arm 1** (64R, 0F): restart_id 1-64, bfs_step = 0 only -> 64 fits
- **Arm 2** (32R, 1F): restart_id 1-32, bfs_step 0-1 -> 64 fits
- **Arm 3** (16R, 3F): restart_id 1-16, bfs_step 0-3 -> 64 fits
- **Arm 4** (8R, 7F): restart_id 1-8, bfs_step 0-7 -> 64 fits
- **Arm 5** (4R, 15F): restart_id 1-4, bfs_step 0-15 -> 64 fits
- **Arm 6** (1R, 63F): restart_id 1, bfs_step 0-63 -> 64 fits

---

## 3. Aggregation Plan

### 3.1 Aggregation Methods (all 4)

1. **max_elbo** -- PIP from highest-ELBO fit
2. **uniform** -- equal-weight average of all PIPs
3. **elbo_softmax** -- softmax(ELBO)-weighted average
4. **cluster_weight** -- JSD-clustered, ELBO-weighted (Method A)

### 3.2 Aggregation Scopes

**Per-arm, per-prior-config**: For each of 6 arms x 4 prior configs, aggregate that arm's 64 fits.

**Per-arm, cross-sigma** (non-EB only): For each arm, pool 3 x 64 = 192 fits across sigma_0.1, sigma_0.2, sigma_0.4. Aggregate all 192 together. (EB excluded -- it's a separate comparison axis.)

### 3.3 Subset Aggregation (Arm 4: 8R x 7F)

For each prefix (r, f) where r in {1,...,8} and f in {0,...,7}:

- Take first r starts, first f refine steps per start
- Subset size = r x (1 + f) fits
- Apply all 4 aggregation methods
- Scopes:
  - Per prior config: 4 x (8 x 8) = **256 subsets**
  - Cross-sigma: 1 x (8 x 8) = **64 subsets**
  - **Total: 320 subset aggregation analyses per dataset**

This produces a performance-vs-compute surface over (r, f) showing the optimal restart/refine allocation at every budget level from 1 to 64.

---

## 4. Implementation Plan

### 4.1 Job Structure: Single Super-Job with Smart Reuse

Run ONE job config with all 240 unique fits per prior config per dataset. Arms are reconstructed post-hoc via filtering.

**Why reuse is essentially free**: The existing codebase does NOT support combining restarts + refines -- `"separate"` mode runs them independently, `"intersect"` mode Cartesian-crosses them (wrong). A new `"restart_refine"` exploration method is needed **regardless** of whether we use reuse. Adding variable per-start refine budgets is ~5 lines on top of that. So we get the 37.5% compute savings for free.

### 4.2 Required Code Changes

#### A. `R/use_cases.R` -- Add exploration method (~5 lines)

- Add `"restart_refine"` to `exploration_catalog()`
- Update `valid_exploration_for_prior()` to accept it for all prior specs

#### B. `R/run_controls.R` -- Run table construction (~40 lines)

- New helper: `build_restart_refine_schedule(refine_schedule)` generates run-table rows:
  - `restart_id` (1-64), `run_type` ("default"/"warm"), `max_refines` (per-start budget)
  - Pre-allocates `bfs_step` rows: 0 (root) + 1..max_refines per start
- Extend `axis_table_for_method()` for `"restart_refine"` -- calls the schedule builder
- New param in `make_job_config()`: `restart_refine_schedule` (tibble with restart_id, max_refines columns)

#### C. `R/run_model.R` -- Execution (~50 lines modified)

- New branch in `execute_dataset_bundle()` for `"restart_refine"` groups (lines 825-903)
- Logic: split group_rows by `restart_id`, for each restart:
  1. Run root fit (uniform or Dirichlet pi depending on restart_id)
  2. Initialize fresh BFS queue from root's CSs
  3. Consume up to `max_refines` BFS steps, each zeroing a CS in the parent's pi
- Warm-start pi generation already works (`run_use_case()` lines 1232-1254 handle both `restart_id` and `blocked_idx`)
- Track metadata: `restart_id`, `bfs_step`, `bfs_depth` in `meta_override`

#### D. `R/restart_refine_analysis.R` -- New file (~150 lines)

Post-hoc analysis functions (needed regardless of reuse):

- `define_arm_schedule()` -- returns the 6-arm definition table
- `reconstruct_arm(results, arm_id)` -- filters fit results to a specific arm's 64 fits
- `subset_aggregation(results, max_r, max_f, agg_methods)` -- runs (r,f) prefix grid
- `cross_sigma_aggregation(results, arm_id, sigma_values)` -- pools across fixed-sigma configs

#### E. Vignettes (3 notebooks)

1. **`run_control_workbook_restart_vs_refine.Rmd`** -- Job config, schedule definition, write artifacts, SLURM commands
2. **`collect_results_restart_vs_refine.Rmd`** -- Aggregate flush files, reconstruct arms, run subset analysis
3. **`visualize_restart_vs_refine.Rmd`** -- Figures (see Section 5)

### 4.3 Memory Optimization (target: 8GB HPC nodes)

**Current memory problem**: `execution_cache` in `execute_dataset_bundle()` (line 649, `run_model.R`) stores **full fit objects** (alpha L x p + mu L x p + mu2 L x p + eval_res + data_bundle reference) for every unique fit. With 960 fits per task, this cache alone uses ~225-400MB and produces **zero cache hits** in this study (all fits are unique by construction).

**What we actually need to keep in memory per task:**

| Item                          | Size (p=1000, L=10)  | When needed              |
|-------------------------------|----------------------|--------------------------|
| X matrix                      | ~5 MB                | Entire task              |
| XtX (if precomputed)          | ~8 MB                | Entire task              |
| Current fit object            | ~300 KB              | During BFS step only     |
| Current fit's CS indices      | ~1 KB                | To spawn BFS children    |
| Accumulated PIPs (all 960)    | ~8 MB                | For aggregation at end   |
| Accumulated ELBOs (all 960)   | ~8 KB                | For aggregation at end   |
| BFS queue (blocked_idx sets)  | ~100 KB              | Per restart's BFS tree   |
| Buffer context (flush lists)  | ~2-5 MB              | Flushed periodically     |
| **Total working set**         | **~25 MB**           |                          |

**Implementation changes for memory efficiency:**

1. **Disable execution cache for `restart_refine` groups.** No cache hits possible -- every fit has a unique (restart_id, blocked_idx) combination. Skip the `assign(cache_key, ...)` call entirely (or wrap in a flag). This saves ~225-400 MB.

2. **Release fit objects immediately after use.** In the BFS loop, after extracting `out$eval_res$effects_filtered$indices` for CS children and after `write_run_outputs()` flushes metrics to disk, set `out <- NULL` and call `gc()` every N iterations. The PIPs and ELBOs are already accumulated in `primary_pips_by_group`/`primary_elbos_by_group`.

3. **Flush buffers after each restart_id's BFS tree.** Currently buffers are flushed at task-level intervals. For restart_refine groups, flush after each restart's tree completes (every 1-64 fits depending on restart_id). This prevents buffer accumulation.

4. **Don't store alpha between restarts.** Refinement only needs the CURRENT fit's CSs (to determine `blocked_idx` for children). The parent's alpha is NOT needed -- `blocked_idx` is stored in the BFS queue node. So each restart's entire BFS tree can run with only ~300KB of fit object in memory at any time.

5. **What to flush to disk per fit:**
   - model_metrics CSV: ELBO, AUPRC, power, CS metrics (already implemented)
   - PIPs: keep in memory for aggregation (8 bytes x p per fit)
   - alpha matrices: **do NOT save to disk** unless explicitly requested. They are ~80KB per fit x 960 = 75MB of I/O that's not needed for any analysis in this study.

**Estimated peak memory with optimizations:** ~30-40 MB working set + R overhead (~500MB) = well within 8GB.

### 4.4 SLURM Sizing

- 450 SLURM array tasks (1 dataset per task)
- Per task: 240 fits x 4 prior configs = 960 IBSS fits
- Estimated wall time: ~30-80 min per task (L=10, p~500-1200)
- Request: 3 hours, 4 GB memory (reduced from 6 GB thanks to memory optimizations)

---

## 5. Analysis & Figures

### 5.1 Primary Comparison (6 arms at equal budget)

- **AUPRC vs. start fraction**: 6-point curve, faceted by agg_method and prior_config
- **Power @ FDR threshold** vs. start fraction
- **Diversity metrics** (mean JSD, n_clusters, mean_pip_var) vs. start fraction
- **Best single fit quality** (max-ELBO fit's AUPRC) vs. start fraction

### 5.2 Arm 4 Performance Surface (subset analysis)

- **Heatmap**: AUPRC over (r, f) grid, faceted by agg_method
- **Iso-budget contours**: performance at constant total fits (e.g., budget=16: (8,1) vs (4,3) vs (2,7))
- **Diminishing returns**: per-strategy performance vs. cumulative fits

### 5.3 Cross-Sigma Pooling

- Per-sigma vs. cross-sigma aggregation quality comparison
- Does pooling across sigma values help beyond per-sigma aggregation?

### 5.4 EB vs. Fixed Comparison

- EB (prior_update_method="var") arm performance vs. each fixed-sigma arm
- Does EB's variance learning help or hurt aggregation? (Per A0 vignette: EB distorts ELBO-softmax weights)

### 5.5 Dataset Stratification

- All above stratified by M1 (LD complexity) and z-score difficulty metrics
- Key question: does the optimal restart/refine ratio depend on dataset difficulty?

---

## 6. Verification Plan

1. **Smoke test**: Run on 3 matrices, 1 seed, 1 prior config. Verify 240 fits produced, arm reconstruction yields 64 fits each, all aggregation methods run.
2. **Seed determinism**: Verify start k produces identical fit across different arm reconstructions (same PIPs, same ELBO).
3. **BFS correctness**: Verify refine step ordering matches BFS (breadth-first, not depth-first). Check that start 1's first 7 refines match arm 4's refines for start 1.
4. **Budget accounting**: For each arm, verify exactly 64 fits after filtering. For each (r,f) subset, verify r*(1+f) fits.
5. **Full run**: 150 matrices, 3 seeds, all 4 prior configs. Check aggregation outputs and metrics.

---

## 7. Resolved Design Decisions

- **Scenario**: scenario_1 only (n=600)
- **Architecture**: oligogenic (p_star/y_noise fixed by architecture, no grid)
- **BFS counting**: 1 refine = 1 BFS step = 1 fit
- **BFS ordering**: keep default (effect index order from `effects_filtered$indices`)
- **Reuse**: Yes -- needed anyway since restart+refine combination requires new code regardless
- **Seeds**: `seeds = c(1, 2, 3)` as base seeds (standard convention)
- **Dirichlet concentration**: `alpha_concentration = 0.01` (sparse draws -- most mass on few variants)
  - **NOTE**: Michael may run a small pilot study first to calibrate this value before the full run. The job config should make this easy to change.

---

## 8. Key Files to Modify

| File | Changes |
|------|---------|
| `test_susine/R/use_cases.R` | Add `"restart_refine"` to catalogs |
| `test_susine/R/run_controls.R` | `build_restart_refine_schedule()`, extend `axis_table_for_method()` |
| `test_susine/R/run_model.R` | New BFS-per-restart branch in `execute_dataset_bundle()`, memory opts |
| `test_susine/R/restart_refine_analysis.R` | **New file**: arm reconstruction, subset aggregation, cross-sigma |
| `test_susine/vignettes/simulation pipeline/` | 3 new Rmd notebooks |
