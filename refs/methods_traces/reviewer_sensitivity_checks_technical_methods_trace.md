# Reviewer Sensitivity Checks Technical Trace

Date written: 2026-06-24
Last updated: 2026-06-25

This document traces the reviewer-requested sensitivity checks added after the
main simulation and real-data analyses:

- simulated annotation-contamination sensitivity;
- simulated genetic-architecture sensitivity;
- real-data annotation/z alignment diagnostic among low-SuSiE-PIP variants.

The trace is implementation-facing. Its purpose is to preserve enough detail to
audit exactly how the paper-facing sensitivity figures and supporting CSVs were
generated, where the code lives, and what output artifacts should be treated as
the paper receipts.

## 1. Scope, provenance, and caveats

Main manuscript context:

- `../Writings/second draft/susine_second_draft.tex`

Paper artifact root used for this trace:

- `../Writings/plots/`

Final copied annotation-sensitivity artifacts:

- `../Writings/plots/annotation_sensitivity/annotation_sensitivity_metrics.png`
- `../Writings/plots/annotation_sensitivity/annotation_z_correlations.png`
- `../Writings/plots/annotation_sensitivity/rms_invariant_by_lambda.png`
- `../Writings/plots/annotation_sensitivity/metric_summary.csv`
- `../Writings/plots/annotation_sensitivity/baseline_summary.csv`
- `../Writings/plots/annotation_sensitivity/annotation_gate.csv`
- `../Writings/plots/annotation_sensitivity/annotation_gate_overall.csv`
- `../Writings/plots/annotation_sensitivity/prior_gate.csv`
- `../Writings/plots/annotation_sensitivity/lambda0_gate.csv`
- `../Writings/plots/annotation_sensitivity/annotation_diagnostics_full.csv`
- `../Writings/plots/annotation_sensitivity/*_plot_data.csv`

Final copied architecture-sensitivity artifacts:

- `../Writings/plots/architecture sensitivity/sparse_architecture_metrics.png`
- `../Writings/plots/architecture sensitivity/sparse_architecture_by_total_h2.png`
- `../Writings/plots/architecture sensitivity/diffuse_architecture_metrics.png`
- `../Writings/plots/architecture sensitivity/metric_summary.csv`
- `../Writings/plots/architecture sensitivity/baseline_summary.csv`
- `../Writings/plots/architecture sensitivity/prior_gate.csv`
- `../Writings/plots/architecture sensitivity/confusion_presence.csv`
- `../Writings/plots/architecture sensitivity/*_plot_data.csv`

Real-data low-PIP annotation diagnostic workbook:

- `vignettes/real data pipeline/real_data_annotation_low_pip_diagnostic.Rmd`

Expected real-data low-PIP copied artifacts after rendering that workbook:

- `../Writings/plots/real_data_case_study/paper_real_data_annotation_low_pip_alignment.csv`
- `../Writings/plots/real_data_case_study/paper_real_data_annotation_low_pip_alignment_summary.csv`
- `../Writings/plots/real_data_case_study/paper_real_data_annotation_low_pip_alignment.png`

Trace-time caveats:

- The copied paper artifact directory contains the annotation and architecture
  sensitivity outputs listed above.
- As of the 2026-06-25 update, the real-data low-PIP diagnostic artifacts
  (`paper_real_data_annotation_low_pip_alignment.csv`,
  `paper_real_data_annotation_low_pip_alignment_summary.csv`,
  `paper_real_data_annotation_low_pip_alignment.png`) are present under
  `../Writings/plots/real_data_case_study/`. The diagnostic remains a one-off
  post-hoc workbook; at original trace time (2026-06-24) these copied artifacts
  were not yet present in the local plot directory.
- The paper-copy directories preserve the final figures and key CSV receipts.
  They do not by themselves preserve the raw SLURM parent job IDs unless
  `arm_index.csv` / `job_index.csv` or the original `output/` job trees are also
  retained.
- The fit RDS objects are intentionally not required for these paper artifacts.
  The simulations write metrics, confusion bins, diagnostics, and plot data; the
  final paper figures are generated from those aggregated tables.

## 2. Scientific purpose

The reviewer checks probe whether the main SuSiNE conclusion depends on two
specific assumptions.

1. Annotation contamination:
   - Does performance remain stable when the simulated annotation vector is
     artificially aligned with the marginal phenotype association statistic
     `z` in variants that should be null?
   - Does the same construction behave as a positive control when applied to
     causal variants?

2. Genetic architecture:
   - Does the annotation-informed C-CS procedure remain useful under very sparse
     architectures with only 1, 3, or 5 causal variants?
   - Does it remain useful under a diffuse architecture with 100 causal variants
     and varying total local heritability?

3. Real-data diagnostic:
   - In the GTEx Lung case study, are annotations spuriously aligned with
     marginal association statistics among variants that the SuSiE anchor assigns
     very low PIP?

The simulation checks use the same sampled genotype-matrix substrate as the main
simulation work: 150 genotype matrices, four phenotype seeds per matrix, and
therefore 600 simulated datasets per grid setting.

## 3. Shared implementation path

The sensitivity jobs use the standard simulation pipeline:

1. Run-control workbook builds a `make_job_config()` object.
2. The generated SLURM array calls the standard task runner.
3. `R/run_model.R::run_task()` executes each assigned dataset bundle and writes
   staged outputs.
4. A collect workbook calls `R/collect_results.R::aggregate_staging_outputs()`.
5. A visualization workbook reads the combined outputs and saves figures plus
   plot-data CSVs.
6. Selected files were copied into `../Writings/plots/`.

Core shared code:

- `R/run_controls.R`
  - adds annotation-contamination columns to run tables;
  - includes those columns in `group_key`;
  - carries `write_prior_diagnostics` and `buffer_flush_interval` into the job
    config.
- `R/run_model.R`
  - `run_task()` honors `job_config$job$buffer_flush_interval`;
  - `annotation_contamination_spec()` parses arm/lambda/shuffle settings;
  - `annotation_contamination_cache_key()` prevents contaminated annotations
    from being reused across incompatible runs;
  - `apply_annotation_contamination_for_run()` applies the transform before the
    fixed-c grid consumes `mu_0`;
  - prior diagnostics are written only when `write_prior_diagnostics = TRUE`.
- `R/simulate_data.R`
  - `simulate_priors()` generates the annotation vector;
  - `compute_marginal_z_scores()` freezes marginal `z` once per dataset bundle;
  - `contaminate_annotation_vector()` constructs the contaminated annotation;
  - `annotation_contamination_diagnostics()` emits RMS, unchanged-block, and
    correlation diagnostics.
- `R/collect_results.R`
  - propagates annotation-contamination columns through aggregated outputs.

Relevant tests:

- `tests/testthat/test-annotation-contamination.R`
  - block invariants for null and causal arms;
  - RNG preservation for optional shuffled-z control;
  - end-to-end lambda-zero master check against a clean C-CS run;
  - annotation settings in C-CS group keys;
  - recall-at-precision fallback behavior.
- `tests/testthat/test-architecture-sensitivity.R`
  - sparse grid preserves `p_star` and per-causal heritability;
  - diffuse grid preserves total heritability and `diffuse_k`;
  - diffuse data generation uses requested causal count and total heritability;
  - sparse sign-prior mode avoids constant-beta annotation degeneracy.

## 4. Annotation-contamination sensitivity

### 4.1 Run-control workbooks

The two arms are launched as separate SLURM jobs:

- `vignettes/simulation pipeline/run_control_workbook_annotation_contamination_null.Rmd`
- `vignettes/simulation pipeline/run_control_workbook_annotation_contamination_causal.Rmd`

Both use:

- `n_matrices_target <- 150L`
- `seeds <- 1:4`
- `L_grid <- 10L`
- `architecture_grid <- "susie2_oligogenic"`
- `annotation_r2 <- 0.3`
- `inflate_match <- 0.95`
- `lambda_grid <- c(0, 0.1, 0.2, 0.3, 0.5, 0.7, 1.0)`
- `c_grid_4 <- c(0.00, 0.43, 1.07, 1.50)`
- `sigma_grid_4 <- c(0.01, 0.10, 0.40, 1.00)`
- `max_iter <- 100L`
- `bundles_per_task <- 4L`
- `aggregation_methods <- "cluster_weight_credible"`
- `include_overall_pool <- FALSE`
- `credible_set_rho <- 0.95`
- `purity_threshold <- 0.50`
- `write_confusion_bins <- TRUE`
- `write_prior_diagnostics <- TRUE`

The null arm sets:

- `job_name <- "annotation_contamination_null"`
- `annotation_arm <- "null"`

The positive-control causal arm sets:

- `job_name <- "annotation_contamination_causal_positive"`
- `annotation_arm <- "causal"`

Each job fits two specs:

- `C-CS-annotation-contamination`
  - `use_case_ids = "susine_functional_mu"`
  - `exploration_methods = c("c_grid", "sigma_0_2_grid")`
  - `exploration_mode = "intersect"`
  - `K = 16`
  - fixed 4 x 4 `c` by `sigma_0^2` grid.
- `baseline-single`
  - `use_case_ids = "susine_vanilla"`
  - `exploration_methods = "single"`
  - `K = 1`.

### 4.2 Annotation construction

The original simulated annotation is the usual `mu_0` generated by
`simulate_priors()`. Before fitting contaminated C-CS runs, the pipeline computes
one frozen vector of marginal z-scores from the simulated phenotype:

```text
z_j = simple-regression t statistic for y ~ X_j
```

For each requested arm and contamination strength `lambda`, the selected block
is:

- null arm: all non-causal variants;
- causal arm: all causal variants;
- lambda zero or arm `none`: identity transform.

For the active block `B`, the code rescales `z_B` to match the RMS of the
original annotation block:

```text
z_scaled = z_B * rms(a_B) / rms(z_B)
blend    = (1 - lambda) * a_B + lambda * z_scaled
a'_B     = blend * rms(a_B) / rms(blend)
```

The inactive block is copied byte-for-byte from the original annotation. The RMS
renormalization keeps annotation scale fixed, so the sensitivity check tests
directional `z` alignment rather than simply increasing prior-mean magnitude.

The implementation hard-fails if:

- active-block RMS is not preserved to tolerance;
- the causal block changes in the null arm;
- the null block changes in the causal arm.

The diagnostic also reports:

- `cor_a_z_null`;
- `cor_a_z_causal`;
- `cor_a_z_anchor_low_pip`;
- `n_anchor_low_pip`;
- active and inactive block RMS values;
- active and inactive block relative RMS errors;
- unchanged-block flags.

`cor_a_z_anchor_low_pip` uses the deterministic SuSiE anchor from the C-CS grid:
among `c = 0` members, the smallest `sigma_0^2` is selected and variants with
anchor `PIP < 0.01` define the low-PIP set.

### 4.3 Compute reuse

The cache key intentionally avoids recomputing fits that do not depend on the
annotation:

- annotation-free baseline runs carry `NA` contamination columns and are not
  multiplied across lambda;
- `c = 0` C-CS members have annotation-free prior means and are reusable across
  lambda values;
- `c > 0` C-CS members use the contaminated prior mean and are keyed by
  arm/lambda/shuffle settings.

The fit path uses fixed `c` values. The prior gate verifies that the real-data
style `mu_0_scale_factor`/`c_l` estimation path did not silently rescale the
contaminated annotations.

### 4.4 Collection and gates

Collection workbook:

- `vignettes/simulation pipeline/collect_results_workbook_annotation_contamination.Rmd`

It writes combined outputs under:

- `output/annotation_contamination_sensitivity/combined/`

The final copied CSV receipts under `../Writings/plots/annotation_sensitivity/`
show:

- `metric_summary.csv`: 14 rows;
- `annotation_diagnostics_full.csv`: 8400 rows;
- `prior_gate.csv`: `n_rows = 1344000`, `max_abs_c_l_minus_1 = 0`,
  `all_c_l_one = TRUE`;
- `lambda0_gate.csv`: `n_lambda0_rows = 2`, `auprc_range = 0`,
  `recall75_range = 0`;
- `annotation_gate.csv`: 600 diagnostic rows per arm/lambda, all
  `all_untouched_block_ok = TRUE`, and active-block RMS relative errors on the
  order of machine precision.

The main pooled metrics are:

- pooled AUPRC from confusion bins;
- recall at precision >= 0.75;
- annotation-free SuSiE-equivalent baseline for each arm/job.

### 4.5 Visualization outputs

Visualization workbook:

- `vignettes/simulation pipeline/visualize_results_workbook_annotation_contamination.Rmd`

It reads:

- `output/annotation_contamination_sensitivity/combined/metric_summary.csv`
- `output/annotation_contamination_sensitivity/combined/baseline_summary.csv`
- `output/annotation_contamination_sensitivity/combined/annotation_diagnostics_full.csv`
- `output/annotation_contamination_sensitivity/combined/annotation_gate.csv`
- `output/annotation_contamination_sensitivity/combined/lambda0_gate.csv`
- `output/annotation_contamination_sensitivity/combined/prior_gate.csv`

It saves figures under:

- `output/annotation_contamination_sensitivity/figures/`

The paper-copy directory contains:

- `annotation_sensitivity_metrics.png`
  - two-panel AUPRC and recall-at-precision figure;
  - dashed horizontal lines are the annotation-free baseline.
- `annotation_z_correlations.png`
  - mean `cor(annotation, z)` curves for null variants, causal variants, and
    anchor low-PIP variants.
- `rms_invariant_by_lambda.png`
  - construction invariant check on a log scale.
- `metric_summary_plot_data.csv`
- `baseline_summary_plot_data.csv`
- `annotation_z_diagnostic_plot_data.csv`
- `construction_gate_plot_data.csv`

## 5. Genetic-architecture sensitivity

### 5.1 Sparse architecture run-control

Run-control workbook:

- `vignettes/simulation pipeline/run_control_workbook_architecture_sparse.Rmd`

Key settings:

- `job_name <- "architecture_sparse_sensitivity"`
- `n_matrices_target <- 150L`
- `seeds <- 1:4`
- `L_grid <- 10L`
- `architecture_grid <- "susie2_sparse"`
- `p_star_grid <- c(1L, 3L, 5L)`
- `h2_snp_per_causal_grid <- c(0.01, 0.03, 0.05, 0.10)`
- `annotation_r2 <- 0.3`
- `inflate_match <- 0.95`
- `c_grid_4 <- c(0.00, 0.43, 1.07, 1.50)`
- `sigma_grid_4 <- c(0.01, 0.10, 0.40, 1.00)`
- `max_iter <- 100L`
- `bundles_per_task <- 48L`
- `buffer_flush_interval <- 16L`
- `write_confusion_bins <- TRUE`
- `write_prior_diagnostics <- TRUE`

The sparse settings collapse to 150 SLURM array tasks, with 48 dataset bundles
per task and three buffer flushes per task. This was added to satisfy cluster
array-size limits while keeping the full grid.

Sparse positive set:

- all causal variants are positives;
- `p_star <= 5`, so this is a strict sparse-recovery check.

Sparse annotation construction:

- sparse effects use constant positive effects;
- therefore the annotation generator uses `annotation_signal_mode = "sign"`
  rather than the default centered-effect mode;
- causal annotations are generated as a per-variant SNR mixture:

```text
mu_0[causal] = causal_scale *
  (sqrt(annotation_r2) * sign(beta) / rms(sign(beta)) +
   sqrt(1 - annotation_r2) * uncentered_rms_noise)
```

This avoids the `k = 1` oracle edge case and the `k = 2` centered-noise
degeneracy that would occur if constant positive effects were treated as a
centered effect vector.

### 5.2 Diffuse architecture run-control

Run-control workbook:

- `vignettes/simulation pipeline/run_control_workbook_architecture_diffuse.Rmd`

Key settings:

- `job_name <- "architecture_diffuse_sensitivity"`
- `n_matrices_target <- 150L`
- `seeds <- 1:4`
- `L_grid <- 10L`
- `architecture_grid <- "susie2_diffuse"`
- `h2_total_grid <- c(0.05, 0.15, 0.25)`
- `diffuse_k_grid <- 100L`
- `annotation_r2 <- 0.3`
- `inflate_match <- 0.95`
- `c_grid_4 <- c(0.00, 0.43, 1.07, 1.50)`
- `sigma_grid_4 <- c(0.01, 0.10, 0.40, 1.00)`
- `max_iter <- 100L`
- `bundles_per_task <- 12L`
- `buffer_flush_interval <- 4L`
- `write_confusion_bins <- TRUE`
- `write_prior_diagnostics <- TRUE`

The diffuse settings also collapse to 150 SLURM array tasks, with 12 dataset
bundles per task and three buffer flushes per task.

Diffuse positive set:

- the top eight causal variants by absolute effect are treated as positives for
  precision-recall metrics;
- this avoids making the metric dominated by many extremely weak causal effects.

Diffuse annotation construction:

- diffuse effects retain the default `annotation_signal_mode = "effect"`, so
  the annotation is correlated with heterogeneous causal effect sizes.

### 5.3 Fitted specs

Both architecture jobs fit:

- `C-CS-architecture-sparse` or `C-CS-architecture-diffuse`
  - `use_case_ids = "susine_functional_mu"`
  - fixed 4 x 4 `c` by `sigma_0^2` C-CS grid;
  - aggregation method `cluster_weight_credible`.
- `baseline-single`
  - `use_case_ids = "susine_vanilla"`
  - one annotation-free baseline fit per dataset.

### 5.4 Collection and gates

Collection workbook:

- `vignettes/simulation pipeline/collect_results_workbook_architecture_sensitivity.Rmd`

It writes combined outputs under:

- `output/architecture_sensitivity/combined/`

The final copied CSV receipts under
`../Writings/plots/architecture sensitivity/` show:

- `metric_summary.csv`: 15 rows;
- `baseline_summary.csv`: 15 rows;
- `prior_gate.csv`:
  - diffuse C-CS: `n_rows = 288000`,
    `max_abs_c_l_minus_1 = 0`, `all_c_l_one = TRUE`;
  - sparse C-CS: `n_rows = 1152000`,
    `max_abs_c_l_minus_1 = 0`, `all_c_l_one = TRUE`;
- `confusion_presence.csv`: confirms confusion-bin rows for every architecture
  setting and for both C-CS and baseline specs.

The main pooled metrics are:

- pooled AUPRC;
- recall at precision >= 0.75;
- annotation-free baseline metrics aligned to each architecture setting.

### 5.5 Visualization outputs

Visualization workbook:

- `vignettes/simulation pipeline/visualize_results_workbook_architecture_sensitivity.Rmd`

It reads:

- `output/architecture_sensitivity/combined/metric_summary.csv`
- `output/architecture_sensitivity/combined/baseline_summary.csv`

It saves figures under:

- `output/architecture_sensitivity/figures/`

The paper-copy directory contains:

- `sparse_architecture_metrics.png`
  - sparse AUPRC and recall at precision >= 0.75 by causal count and per-causal
    heritability.
- `sparse_architecture_by_total_h2.png`
  - sparse metrics plotted by approximate total locus heritability
    `p_star * h2_snp_per_causal`.
- `diffuse_architecture_metrics.png`
  - diffuse AUPRC and recall at precision >= 0.75 by total local heritability.
- `architecture_metric_plot_data.csv`
- `sparse_architecture_plot_data.csv`
- `diffuse_architecture_plot_data.csv`

## 6. Real-data low-PIP annotation/z diagnostic

One-off workbook:

- `vignettes/real data pipeline/real_data_annotation_low_pip_diagnostic.Rmd`

This workbook runs after the real-data ensemble collect/visualization pipeline.
It does not refit any real-data model.

Default source inputs:

- `output/slurm_output/real_data_ensemble_geometric_n20/52906940/aggregated/paper_real_data_ensemble_summary.csv`
- `output/slurm_output/real_data_ensemble_geometric_n20/52906940/aggregated/variant_posteriors_dataset/`

Environment-variable overrides are supported through:

- `RD_PREVIEW_AGG_DIR`
- `RD_PREVIEW_FIG_DIR`
- `RD_PREVIEW_PAPER_DIR`

For each locus, the workbook:

1. reads `susie_anchor_run_id` from `paper_real_data_ensemble_summary.csv`;
2. filters `variant_posteriors_dataset` to that `locus_id` and anchor `run_id`;
3. keeps variants with anchor `pip < 0.01`;
4. computes:

```text
cor_a_z_anchor_low_pip = cor(annotation_a, z_score)
```

It also records:

- `n_variants`;
- `n_anchor_low_pip`;
- `n_finite_low_pip`;
- `low_pip_fraction`;
- `cor_a_z_anchor_low_pip_p_value`;
- `abs_cor_a_z_anchor_low_pip`;
- `mean_abs_z_anchor_low_pip`;
- `mean_abs_a_anchor_low_pip`.

Expected outputs:

- `paper_real_data_annotation_low_pip_alignment.csv`
- `paper_real_data_annotation_low_pip_alignment_summary.csv`
- `paper_real_data_annotation_low_pip_alignment.png`

Copied paper receipts under `../Writings/plots/real_data_case_study/` (as of
2026-06-25) show:

- `paper_real_data_annotation_low_pip_alignment.csv`: 20 locus rows, with columns
  `locus_id`, `gene_name`, `susie_anchor_run_id`, `n_variants`,
  `n_anchor_low_pip`, `n_finite_low_pip`, `low_pip_fraction`,
  `cor_a_z_anchor_low_pip`, `cor_a_z_anchor_low_pip_p_value`,
  `abs_cor_a_z_anchor_low_pip`, `mean_abs_z_anchor_low_pip`,
  `mean_abs_a_anchor_low_pip`;
- `paper_real_data_annotation_low_pip_alignment_summary.csv`: `n_loci = 20`,
  `n_loci_with_finite_cor = 20`, `median_cor` ~= 0.0124, `max_abs_cor` ~= 0.0807,
  `median_abs_cor` ~= 0.0174, `median_n_anchor_low_pip = 6275.5`,
  `min_n_anchor_low_pip = 1580`. The near-zero correlations among low-PIP
  variants are the intended negative-control reading: annotations are not
  spuriously aligned with marginal `z` where the SuSiE anchor assigns very low
  PIP.

Default save location:

- `output/slurm_output/real_data_ensemble_geometric_n20/52906940/figures/real_data_ensemble/overall/`

If `../Writings/plots/real_data_case_study/` exists, the workbook copies those
three files there as paper-local artifacts.

## 7. How to regenerate the paper-facing outputs

### 7.1 Annotation sensitivity

1. Render both run-control workbooks after setting any cluster-specific
   submission details:
   - `run_control_workbook_annotation_contamination_null.Rmd`
   - `run_control_workbook_annotation_contamination_causal.Rmd`
2. Submit the generated SLURM jobs.
3. Put the returned parent job IDs into:
   - `collect_results_workbook_annotation_contamination.Rmd`
4. Render the collect workbook.
5. Render:
   - `visualize_results_workbook_annotation_contamination.Rmd`
6. Copy the required files from:
   - `output/annotation_contamination_sensitivity/combined/`
   - `output/annotation_contamination_sensitivity/figures/`
   into:
   - `../Writings/plots/annotation_sensitivity/`

### 7.2 Architecture sensitivity

1. Render both run-control workbooks:
   - `run_control_workbook_architecture_sparse.Rmd`
   - `run_control_workbook_architecture_diffuse.Rmd`
2. Submit the generated SLURM jobs.
3. Put the returned parent job IDs into:
   - `collect_results_workbook_architecture_sensitivity.Rmd`
4. Render the collect workbook.
5. Render:
   - `visualize_results_workbook_architecture_sensitivity.Rmd`
6. Copy the required files from:
   - `output/architecture_sensitivity/combined/`
   - `output/architecture_sensitivity/figures/`
   into:
   - `../Writings/plots/architecture sensitivity/`

### 7.3 Real-data low-PIP diagnostic

1. Ensure the real-data ensemble aggregated outputs exist for
   `real_data_ensemble_geometric_n20/52906940`, or set `RD_PREVIEW_AGG_DIR` to a
   compatible aggregated directory.
2. Render:
   - `vignettes/real data pipeline/real_data_annotation_low_pip_diagnostic.Rmd`
3. Confirm the three files named `paper_real_data_annotation_low_pip_alignment*`
   exist under either:
   - the real-data job `overall/` figure directory; or
   - `../Writings/plots/real_data_case_study/`.

## 8. Audit checklist

For annotation contamination:

- Confirm `write_prior_diagnostics = TRUE` in both annotation run-control
  workbooks.
- Confirm `prior_gate.csv` has `all_c_l_one = TRUE`.
- Confirm `lambda0_gate.csv` has zero AUPRC and recall ranges.
- Confirm `annotation_gate.csv` has `all_untouched_block_ok = TRUE` for every
  arm/lambda.
- Confirm `annotation_diagnostics_full.csv` has expected arms and lambda grid.
- Confirm figures were made from `metric_summary_plot_data.csv`,
  `baseline_summary_plot_data.csv`, `annotation_z_diagnostic_plot_data.csv`, and
  `construction_gate_plot_data.csv`.

For architecture sensitivity:

- Confirm `write_prior_diagnostics = TRUE` in both architecture run-control
  workbooks.
- Confirm `prior_gate.csv` has `all_c_l_one = TRUE` for sparse and diffuse.
- Confirm `confusion_presence.csv` has rows for every architecture setting and
  both C-CS and baseline specs.
- Confirm sparse plot data covers:
  - `p_star` in `{1, 3, 5}`;
  - `h2_snp_per_causal` in `{0.01, 0.03, 0.05, 0.10}`.
- Confirm diffuse plot data covers:
  - `h2_total` in `{0.05, 0.15, 0.25}`;
  - `diffuse_k = 100`.

For the real-data diagnostic:

- Confirm `paper_real_data_ensemble_summary.csv` has `locus_id` and
  `susie_anchor_run_id`.
- Confirm `variant_posteriors_dataset` has `locus_id`, `run_id`, `variant_id`,
  `pip`, `z_score`, and `annotation_a`.
- Confirm the copied output files named
  `paper_real_data_annotation_low_pip_alignment*` exist in the paper plot
  directory before citing them as paper artifacts.

