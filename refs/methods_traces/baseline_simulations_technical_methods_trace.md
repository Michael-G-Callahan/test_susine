# Baseline Simulation Study Technical Trace

Date written: 2026-05-20

This document traces the baseline single-fit simulation study as implemented in
the `test_susine` codebase and as represented in the current manuscript draft.
It is implementation-facing: the goal is to preserve enough detail to audit how
the paper-facing baseline figures and quoted numbers were produced.

The baseline HPC output directory for job `baseline_sims_screen/51509956` was
not present in this local checkout when this trace was written. The trace is
therefore based on source workbooks, helper code, the manuscript draft, final
copied plot files, and the copied paper-side `top_settings_by_family.csv`.

## 1. Scope, provenance, and caveats

Main files traced:

- Manuscript context:
  - `../Writings/second draft/susine_second_draft.tex`
- Run-control workbook:
  - `vignettes/simulation pipeline/run_control_workbook_baseline_sims.Rmd`
- Collection workbook:
  - `vignettes/simulation pipeline/collect_results_workbook_baseline_sims.Rmd`
- Visualization workbook:
  - `vignettes/simulation pipeline/visualize_results_workbook_baseline_sims.Rmd`
- Final copied paper artifacts:
  - `../Writings/plots/baseline_sims/fig_baseline_combined.png`
  - `../Writings/plots/baseline_sims/pip_calibration_top_settings.png`
  - `../Writings/plots/baseline_sims/auprc_by_annotation.png`
  - `../Writings/plots/baseline_sims/auprc_non_annotation_methods.png`
  - `../Writings/plots/baseline_sims/tau_sweep.png`
  - `../Writings/plots/baseline_sims/top_settings_by_family.csv`
- Core helper code:
  - `R/run_controls.R`
  - `R/run_model.R`
  - `R/simulate_data.R`
  - `R/evaluation_metrics.R`
  - `R/heritability.R`
  - `R/visualize_results.R`
  - `R/collect_results.R`
  - `R/use_cases.R`
- Relevant tests:
  - `tests/testthat/test-baseline-sims-screening.R`

Important caveats:

- The expected local generated directories do not exist:
  - `output/slurm_output/baseline_sims_screen/51509956`
  - `output/run_history/baseline_sims_screen/51509956/job_config.json`
- The exact HPC consolidated inputs used by the visualization workbook were
  therefore not locally re-read. If exact regenerated tables are needed, the
  narrow artifact request should be only:
  `output/slurm_output/baseline_sims_screen/51509956/consolidated/`.
- The final paper plot directory does contain copied outputs, including
  `top_settings_by_family.csv`, which preserves the paper-facing selected
  settings and numeric AUPRC/TPR summaries.
- Git status could not be checked because this checkout triggered Git's
  safe-directory ownership guard. I did not mutate global Git config to bypass
  that guard.

## 2. Manuscript context

The manuscript frames the baseline simulations as the first realistic
simulation screen after the toy pathology examples. The baseline study is used
to answer which single-fit SuSiE-family specifications are worth carrying
forward.

The manuscript states that the shared simulation substrate consists of:

- 150 genotype matrices;
- 4 phenotype seeds per matrix;
- 600 total simulated datasets;
- `n = 600` individuals;
- approximately `p = 1000` SNPs per locus;
- SuSiE 2.0-style oligogenic architecture;
- 3 sparse/strong effects, 5 oligogenic/moderate effects, and 15 polygenic
  effects;
- 23 total causal variants per locus;
- `L = 10` fitted single-effect components;
- maximum 100 IBSS iterations;
- credible-set threshold `rho = 0.95`;
- purity threshold `0.50` for filtered credible-set summaries.

The manuscript baseline results section uses this study to support three
decisions:

- annotation-agnostic SuSiE-like methods are tightly clustered;
- the signed prior-mean channel (`mu_0 = c * a`) is stronger than the unsigned
  functional-pi channel under matched signed annotations;
- later ensemble work should carry forward `susine_functional_mu` and
  `susine_eb_clamped_scale_var_nonneg`, while dropping SuSiE-ash, SuSiE-inf,
  and functional-pi from the main ensemble screen.

The manuscript baseline figures are:

- Fig. baseline AUPRC:
  `../Writings/plots/baseline_sims/fig_baseline_combined.png`
- Fig. baseline calibration:
  `../Writings/plots/baseline_sims/pip_calibration_top_settings.png`

## 3. End-to-end workflow

The baseline simulation study has four stages.

1. Build the job configuration.
   - Workbook:
     `vignettes/simulation pipeline/run_control_workbook_baseline_sims.Rmd`
   - Key helper:
     `R/run_controls.R` (`make_job_config`)
   - Outputs expected under `output/run_history/baseline_sims_screen/<job id>/`.

2. Execute one task per dataset bundle batch on HPC.
   - Script:
     `inst/scripts/run_task.R`
   - Key helper:
     `R/run_model.R` (`run_task`, `execute_dataset_bundle`)
   - Each dataset bundle regenerates the genotype-derived data bundle, phenotype,
     simulated annotations, fits, model/effect metrics, tier metrics, and
     confusion bins.

3. Aggregate staged task outputs.
   - Workbook:
     `vignettes/simulation pipeline/collect_results_workbook_baseline_sims.Rmd`
   - Key helper:
     `R/collect_results.R` (`aggregate_staging_outputs`)
   - Writes consolidated CSVs expected under:
     `output/slurm_output/baseline_sims_screen/51509956/consolidated/`.

4. Render analysis and paper-facing plots.
   - Workbook:
     `vignettes/simulation pipeline/visualize_results_workbook_baseline_sims.Rmd`
   - Reads consolidated CSVs, constructs paper plots under the job figure
     directory, and the final selected plot files were copied to:
     `../Writings/plots/baseline_sims/`.

## 4. Run-control settings

The baseline run-control workbook sets:

- `job_name <- "baseline_sims_screen"`
- `n_matrices_target <- 150L`
- `seeds <- 1:4`
- `L_grid <- 10L`
- `y_noise_grid <- NA_real_`
- `p_star_grid <- NA_integer_`
- `architecture_grid <- "susie2_oligogenic"`
- `annotation_r2_levels <- c(0, 0.3, 0.5)`
- `inflate_match_levels <- c(1, 0.9, 0.8)`
- `sigma_0_2_default <- 0.2`
- `max_iter <- 100L`
- `credible_set_rho <- 0.95`
- `purity_threshold <- 0.50`
- `verbose_file_output <- FALSE`
- `write_snps_parquet <- FALSE`
- `write_confusion_bins <- TRUE`
- `write_tier_cs_metrics <- TRUE`
- `write_scaling_confusion_bins <- FALSE`
- `z_top_k <- 10L`
- `jsd_threshold <- 0.15`

The annotation settings map to manuscript notation as:

- `annotation_r2` = `phi_a`;
- `inflate_match` = `nu_a`.

The study evaluates the `3 x 3` grid:

| `phi_a` | `nu_a` |
|---:|---:|
| 0.0 | 1.0 |
| 0.0 | 0.9 |
| 0.0 | 0.8 |
| 0.3 | 1.0 |
| 0.3 | 0.9 |
| 0.3 | 0.8 |
| 0.5 | 1.0 |
| 0.5 | 0.9 |
| 0.5 | 0.8 |

The PIP confusion-bin threshold grid is:

- `0.00` to `0.05` by `0.005`;
- `0.05` to `0.10` by `0.01`;
- `0.10` to `1.00` by `0.02`.

This grid is used by `compute_confusion_bins()` and then pooled by the
collection workbook to compute AUPRC and TPR@FPR=0.05.

## 5. Baseline specification catalog

The run-control workbook defines eight baseline specs:

| Spec | Use case | Exploration | Grid values |
|---|---|---|---|
| `susie_vanilla_sigma` | `susie_vanilla` | `sigma_0_2_grid` | `0.1, 0.2, 0.4` |
| `susine_vanilla_sigma` | `susine_vanilla` | `sigma_0_2_grid` | `0.1, 0.2, 0.4` |
| `susie_eb_single` | `susie_eb` | `single` | default EB |
| `susie_inf_single` | `susie_inf` | `single` | default |
| `susie_ash_fixed_single` | `susie_ash_fixed` | `single` | fixed sigma default |
| `susine_functional_mu_grid` | `susine_functional_mu` | `c_grid x sigma_0_2_grid` | `c = seq(0, 1.5, length.out = 8)`, `sigma_0_2 = 0.01, 0.03, 0.1, 0.2, 0.4` |
| `susine_eb_clamped_scale_var_nonneg_single` | `susine_eb_clamped_scale_var_nonneg` | `single` | EB estimates nonnegative scale and variance |
| `susie_functional_pi_tau` | `susie_functional_pi` | `tau_grid` | `exp(seq(log(0.01), log(10), length.out = 10))` |

The annotation-agnostic specs have `annotation_r2 = NA` and
`inflate_match = NA`. Annotation-consuming specs are crossed with the full
`3 x 3` annotation grid.

The test `tests/testthat/test-baseline-sims-screening.R` checks that:

- vanilla specs expand over only the sigma grid and not the annotation grid;
- functional-mu specs expand over the Cartesian product of annotation quality,
  `c`, and `sigma_0_2`;
- aggregation preserves both filtered and unfiltered effect metrics with the
  setting-context columns needed downstream.

## 6. Matrix and dataset construction

The run-control workbook builds the matrix catalog with:

- `test_susine:::build_data_matrix_catalog(requested_scenarios = "scenario_1")`
- source summary path:
  `data/sampled_simulated_genotypes/scenario_sampling_summary.csv`

It then tries to use:

- `data/sampled_simulated_genotypes/scenario_sampling_summary_with_m1.csv`

If that file is absent, the workbook computes M1 values with
`test_susine::build_x_metrics_table()` and writes the M1-augmented summary.

Matrix selection is M1-stratified:

1. filter rows with non-missing `M1`;
2. sort by `M1`;
3. take `round(seq(1, nrow(catalog), length.out = n_target))`;
4. keep the corresponding 150 matrices.

With `seeds <- 1:4`, the job config creates 600 dataset bundles.

In the locally available 5-arm derived `dataset_bundles.csv`, the 600 baseline
dataset keys use `phenotype_seed` values such as `10008`, `20016`, `30024`,
and `40032` for matrix 1. These are derived internal seeds, not the literal
`1:4` values from the run-control workbook.

## 7. Oligogenic phenotype generator

The relevant implementation is `simulate_effect_sizes_susie2_oligogenic()` in
`R/simulate_data.R`.

The current code defines three fixed tiers:

| Tier | Count | Draw | Energy fraction in code |
|---|---:|---|---:|
| sparse | 3 | `N(0, 1)` | 0.50 |
| oligogenic | 5 | mixture, 35% `sd = 0.8`, 65% `sd = 0.25` | 0.35 |
| polygenic | 15 | `N(0, 0.08^2)` | 0.15 |

The generator samples 23 non-overlapping causal variants, assigns tier labels,
and rescales each tier to its target L2-energy fraction. The combined causal
index is the sorted union of all tier indices.

For `susie2_oligogenic`, `generate_simulation_data()` and
`generate_data_for_bundle()` set the target phenotype heritability with
`simulate_phenotype_h2()`. The code path uses `h2_total = 0.25` for this
architecture.

## 8. Annotation simulation

Annotations are generated by `simulate_priors()` in `R/simulate_data.R`.

At causal variants:

- `annotation_r2` controls the squared correlation between the simulated signed
  annotation and the true causal effect vector;
- the causal annotation vector is built by mixing a standardized causal-effect
  signal with orthogonal noise to target the requested correlation;
- the code preserves a causal annotation scale based on causal effect variance.

At non-causal variants:

- `inflate_match` controls non-causal annotation variance relative to causal
  annotation variance;
- non-causal annotations are centered noise with that variance.

The resulting `mu_0` vector is the signed annotation vector used by
`susine_functional_mu` after multiplication by the fitted or fixed scale `c`.
The same annotation vector is converted to unsigned inclusion weights for
functional-pi methods through the use-case-specific prior construction.

## 9. Model execution and metrics

The HPC task script dispatches to `run_task()` in `R/run_model.R`.

For each dataset bundle:

1. `generate_data_for_bundle()` regenerates `X`, `y`, `beta`, `causal_idx`,
   `mu_0`, `sigma_0_2`, and annotation metadata.
2. `execute_dataset_bundle()` groups run rows by use case.
3. `run_use_case()` dispatches to either `susieR` or local `susine`, depending
   on the use-case catalog.
4. `evaluate_model()` computes model-level and effect-level metrics.
5. `write_run_outputs()` writes or buffers model metrics, effect metrics,
   unfiltered effect metrics, tier credible-set metrics, and confusion bins.

Important metric definitions from `R/evaluation_metrics.R`:

- credible set: smallest descending-alpha prefix with cumulative alpha at least
  `rho`;
- purity: minimum absolute pairwise genotype correlation inside the credible
  set;
- coverage: whether a credible set contains at least one true causal variant;
- model power: fraction of true causals captured by at least one credible set;
- combined PIP: `1 - prod_l(1 - alpha_lj)`;
- effect diffuseness:
  - `effect_pip_entropy`;
  - `effect_pip_entropy_core95`;
  - `effect_k_eff_signal = exp(effect_pip_entropy)`;
  - `effect_k_eff_signal_core95 = exp(effect_pip_entropy_core95)`;
- effect accuracy:
  - `accuracy_ratio = max alpha on true causal variants / max alpha on any variant`.

`evaluate_model()` also emits local-genetic-variance (PVE/heritability) columns
on each model-metric row. Two estimands are emitted:

- legacy `hg2` (still emitted) = `var(fitted_y) / var(y)`, clipped to `[0, 1]`,
  via `estimate_hg2()` in `R/evaluation_metrics.R`;
- a corrected local-genetic-variance decomposition (the now-reported estimand),
  via `hg2_components()` with its engine in `R/heritability.R`. It emits
  `hg2_postmean = var(E[Xβ|y]) / var(y)`, `hg2_uncertainty` (a within-fit
  posterior-variance correction), and
  `hg2_expected_pve = hg2_postmean + hg2_uncertainty = E[var(Xβ)|y] / var(y)`
  for a single fit. An additional `hg2_between_fit` term appears only under
  multi-fit/ensemble aggregation. Provenance: commits `4e0721b` / `236820b`
  (2026-06-09).

`evaluate_model()` produces both unfiltered and purity-filtered views. The
filtered view keeps credible sets with `purity >= 0.50`.

## 10. Confusion bins and top-8 convention

The confusion-bin implementation is `compute_confusion_bins()` in
`R/run_model.R`.

The baseline study uses a top-8 causal mask for PIP-threshold classification:

- within each dataset, causal variants are ordered by `abs(beta)`;
- the top 8 causal variants are treated as positives;
- the `causal_mask` argument replaces the full causal label vector inside
  `compute_confusion_bins()`;
- positives are the top-8 mask; negatives are true non-causal variants only.
  Variants that are truly causal but outside the top-8 mask (the weak /
  polygenic-tier causals) are DROPPED from the confusion-bin table entirely --
  neither positive nor negative. The code does
  `scored <- as.integer(seq_len(n) %in% causal_mask)` and then
  `keep <- !(causal == 1L & scored == 0L)`, removing those weak causals so
  AUPRC is not penalized for them.

This convention corresponds to the sparse plus oligogenic tiers under the
23-causal oligogenic architecture. The 5-arm drift PR regeneration notebook is
intended to match this same convention.

This is a point to watch when writing prose. The source code path for
confusion-bin AUPRC uses the top-8 mask as the positive set, true non-causal
variants as the negative set, and DROPS the weak/polygenic-tier causals (true
causals outside the top-8 mask) from the table entirely. Manuscript captions
that describe the polygenic-tier causals as excluded from negatives are
therefore correct: those weak causals are neither positives nor negatives.

## 11. Collection workflow

The collection workbook uses:

- `job_name <- "baseline_sims_screen"`
- `parent_job_id <- "51509956"`
- `output_root <- here("output")`

It expects:

- job directory:
  `output/slurm_output/baseline_sims_screen/51509956`
- job config:
  `output/run_history/baseline_sims_screen/51509956/job_config.json`

It calls:

```r
test_susine::aggregate_staging_outputs(
  job_name = job_name,
  parent_job_id = parent_job_id,
  output_root = output_root,
  validate = TRUE,
  output_dir = aggregated_dir,
  write_snps_dataset = FALSE
)
```

The collection workbook then reads the aggregated files:

- `model_metrics.csv`
- `effect_metrics.csv`
- `effect_metrics_unfiltered.csv`
- `confusion_bins.csv`
- `dataset_metrics.csv`
- `tier_cs_metrics.csv`
- `validation.csv`

It derives baseline-specific identifiers:

- `method_family`
- `annotation_label`
- `setting_label`

The label logic is:

- annotation label:
  - `no_annotation` when `annotation_r2` or `inflate_match` is missing;
  - otherwise `r2=<annotation_r2> | inflate=<inflate_match>`.
- setting label:
  - sigma-based labels for vanilla, SuSiE-inf, and SuSiE-ash;
  - tau and sigma labels for functional-pi;
  - c and sigma labels for functional-mu;
  - method family only for EB clamped SuSiNE.

It computes:

- per-dataset AUPRC grouped by dataset, run, method, annotation, and tuning
  setting;
- pooled overall AUPRC grouped by method, annotation, and tuning setting;
- pooled-by-annotation AUPRC;
- corresponding TPR@FPR=0.05 tables using pooled confusion bins.

The expected consolidated exports are:

- `model_metrics_full.csv`
- `effect_metrics_full.csv`
- `effect_metrics_unfiltered_full.csv`
- `confusion_bins_full.csv`
- `dataset_metrics_full.csv`
- `tier_cs_metrics_full.csv`
- `validation_full.csv`
- `auprc_per_dataset.csv`
- `auprc_pooled_overall.csv`
- `auprc_pooled_by_annotation.csv`
- `tpr05_per_dataset.csv`
- `tpr05_pooled_overall.csv`
- `tpr05_pooled_by_annotation.csv`
- `screening_summary.csv`

These consolidated files were not present locally at trace time.

## 12. Visualization workflow

The visualization workbook reads the consolidated CSVs and writes figures under:

- `output/slurm_output/baseline_sims_screen/51509956/figures/`

It creates these figure subdirectories:

- `susie2_plots`
- `baseline_performance`
- `baseline_performance/confirmation`
- `dataset_difficulty`
- `effect_diagnostics`

The paper-facing plots are generated as follows.

### 12.1 PIP calibration

Output in job figure directory:

- `susie2_plots/pip_calibration_top_settings.png`

Final copied paper file:

- `../Writings/plots/baseline_sims/pip_calibration_top_settings.png`

The workbook selects:

- annotation-consuming methods at `annotation_r2 = 0.3`, `inflate_match = 0.9`;
- the best AUPRC row within each annotation-consuming method family;
- non-annotation methods `susie_inf`, `susie_ash_fixed`, and `susie_eb` at their
  best AUPRC rows;
- `susie_vanilla` fixed at `sigma_0_2 = 0.2`.

It bins PIP thresholds into 0.10-wide bins, computes observed causal fraction
per dataset and bin, then averages across datasets with standard-error
intervals.

### 12.2 Non-annotation methods

Output in job figure directory:

- `baseline_performance/auprc_non_annotation_methods.png`
- `baseline_performance/tpr05_non_annotation_methods.png`

Final copied paper file:

- `../Writings/plots/baseline_sims/auprc_non_annotation_methods.png`

The AUPRC panel includes:

- `susie_vanilla`
- `susie_inf`
- `susie_ash_fixed`
- `susie_eb`

### 12.3 Annotation response

Output in job figure directory:

- `baseline_performance/auprc_by_annotation.png`
- `baseline_performance/tpr05_annotation_consuming.png`

Final copied paper file:

- `../Writings/plots/baseline_sims/auprc_by_annotation.png`

The panel focuses on `susine_functional_mu` across the c and sigma grid, with
reference lines for:

- best SuSiE vanilla;
- best SuSiE functional-pi within each annotation facet;
- SuSiNE-EB within each annotation facet.

### 12.4 Tau sweep

Output in job figure directory:

- `baseline_performance/tau_sweep.png`

Final copied paper file:

- `../Writings/plots/baseline_sims/tau_sweep.png`

This panel shows functional-pi AUPRC versus tau on a log scale, faceted by
`inflate_match` and colored by `annotation_r2`, with reference lines for best
vanilla and best functional-mu.

### 12.5 Combined baseline figure

Output in job figure directory:

- `baseline_performance/fig_baseline_combined.png`

Final copied paper file:

- `../Writings/plots/baseline_sims/fig_baseline_combined.png`

The combined figure uses `patchwork`:

- panel (a): non-annotation methods;
- panel (b): functional-pi tau sweep;
- panel (c): functional-mu annotation response.

### 12.6 Top settings table

Output in job figure directory:

- `baseline_performance/top_settings_by_family.csv`
- `baseline_performance/top_settings_global_top5.csv`

Final copied paper file:

- `../Writings/plots/baseline_sims/top_settings_by_family.csv`

The table chooses:

- one best row per method family for annotation-agnostic methods;
- one best row per `(method_family, annotation_label)` for annotation-consuming
  methods.

## 13. Final copied baseline artifacts

Local final copied paper files under `../Writings/plots/baseline_sims`:

| File | Size | Last write time |
|---|---:|---|
| `auprc_by_annotation.png` | 329,923 | 2026-04-26 18:37:58 |
| `auprc_non_annotation_methods.png` | 45,423 | 2026-04-26 18:38:00 |
| `fig_baseline_combined.png` | 580,093 | 2026-04-27 11:12:04 |
| `pip_calibration_top_settings.png` | 196,084 | 2026-04-26 18:38:53 |
| `tau_sweep.png` | 134,877 | 2026-04-26 18:38:03 |
| `top_settings_by_family.csv` | 3,675 | 2026-04-26 18:38:05 |

## 14. Paper-facing selected settings and numbers

The final copied `top_settings_by_family.csv` contains the following
paper-facing rows:

```csv
method_family,annotation_label,setting_label,AUPRC,tpr_fpr05
susie_ash_fixed,no_annotation,susie_ash_fixed | sigma=0.2,0.234942,0.515417
susie_eb,no_annotation,susie_eb,0.228315,0.485208
susie_functional_pi,r2=0.0 | inflate=0.8,susie_functional_pi | tau=4.6400,0.245036,0.531458
susie_functional_pi,r2=0.0 | inflate=0.9,susie_functional_pi | tau=4.6400,0.243008,0.530625
susie_functional_pi,r2=0.0 | inflate=1.0,susie_functional_pi | tau=10.0000,0.243364,0.531042
susie_functional_pi,r2=0.3 | inflate=0.8,susie_functional_pi | tau=0.4640,0.271397,0.545208
susie_functional_pi,r2=0.3 | inflate=0.9,susie_functional_pi | tau=0.4640,0.265318,0.538958
susie_functional_pi,r2=0.3 | inflate=1.0,susie_functional_pi | tau=1.0000,0.256406,0.536250
susie_functional_pi,r2=0.5 | inflate=0.8,susie_functional_pi | tau=0.4640,0.293128,0.558750
susie_functional_pi,r2=0.5 | inflate=0.9,susie_functional_pi | tau=0.4640,0.283158,0.542708
susie_functional_pi,r2=0.5 | inflate=1.0,susie_functional_pi | tau=0.4640,0.273789,0.536250
susie_inf,no_annotation,susie_inf | sigma=0.2,0.247328,0.550208
susie_vanilla,no_annotation,susie_vanilla | sigma=0.1,0.243601,0.533750
susine_eb_clamped_scale_var_nonneg,r2=0.0 | inflate=0.8,susine_eb_clamped_scale_var_nonneg,0.219573,0.490417
susine_eb_clamped_scale_var_nonneg,r2=0.0 | inflate=0.9,susine_eb_clamped_scale_var_nonneg,0.213792,0.486667
susine_eb_clamped_scale_var_nonneg,r2=0.0 | inflate=1.0,susine_eb_clamped_scale_var_nonneg,0.220685,0.488542
susine_eb_clamped_scale_var_nonneg,r2=0.3 | inflate=0.8,susine_eb_clamped_scale_var_nonneg,0.258100,0.506667
susine_eb_clamped_scale_var_nonneg,r2=0.3 | inflate=0.9,susine_eb_clamped_scale_var_nonneg,0.261627,0.507917
susine_eb_clamped_scale_var_nonneg,r2=0.3 | inflate=1.0,susine_eb_clamped_scale_var_nonneg,0.257765,0.505625
susine_eb_clamped_scale_var_nonneg,r2=0.5 | inflate=0.8,susine_eb_clamped_scale_var_nonneg,0.267929,0.508125
susine_eb_clamped_scale_var_nonneg,r2=0.5 | inflate=0.9,susine_eb_clamped_scale_var_nonneg,0.266709,0.508542
susine_eb_clamped_scale_var_nonneg,r2=0.5 | inflate=1.0,susine_eb_clamped_scale_var_nonneg,0.264265,0.510208
susine_functional_mu,r2=0.0 | inflate=0.8,susine_functional_mu | c=0.000 | sigma=0.01,0.246172,0.546250
susine_functional_mu,r2=0.0 | inflate=0.9,susine_functional_mu | c=0.000 | sigma=0.01,0.246172,0.546250
susine_functional_mu,r2=0.0 | inflate=1.0,susine_functional_mu | c=0.000 | sigma=0.01,0.246172,0.546250
susine_functional_mu,r2=0.3 | inflate=0.8,susine_functional_mu | c=0.643 | sigma=0.01,0.312642,0.562708
susine_functional_mu,r2=0.3 | inflate=0.9,susine_functional_mu | c=0.643 | sigma=0.01,0.313759,0.568542
susine_functional_mu,r2=0.3 | inflate=1.0,susine_functional_mu | c=0.857 | sigma=0.01,0.308138,0.558333
susine_functional_mu,r2=0.5 | inflate=0.8,susine_functional_mu | c=0.643 | sigma=0.01,0.339148,0.572500
susine_functional_mu,r2=0.5 | inflate=0.9,susine_functional_mu | c=0.643 | sigma=0.01,0.336866,0.570625
susine_functional_mu,r2=0.5 | inflate=1.0,susine_functional_mu | c=0.857 | sigma=0.01,0.334232,0.564167
susine_vanilla,no_annotation,susine_vanilla | sigma=0.1,0.243665,0.535000
```

These rows support the manuscript statements that:

- SuSiE-inf is the strongest annotation-agnostic baseline by AUPRC
  (`0.247328`);
- SuSiE vanilla is essentially tied with SuSiNE vanilla around `0.244`;
- at the AlphaGenome-calibration point `r2=0.3 | inflate=0.9`,
  `susine_functional_mu` reaches `0.313759`;
- at that same annotation point, functional-pi reaches `0.265318`;
- at the high-quality drift-analysis point `r2=0.5 | inflate=0.9`,
  `susine_functional_mu` reaches `0.336866`, functional-pi reaches
  `0.283158`, and SuSiNE-EB reaches `0.266709`.

## 15. Diagnostics and non-paper plots

The visualization workbook also writes diagnostic plots that were not copied to
the final `Writings/plots/baseline_sims` directory:

- TPR/FPR overlays for top settings;
- credible-set power and FDR by top-N tier;
- SuSiE versus SuSiNE vanilla confirmation plots;
- TPR@FPR annotation response;
- dataset difficulty versus `M1` and z-score metrics;
- effect-diagnostic plots for diffuseness and accuracy.

These plots are useful for audit and sanity checking but are not the final
paper figures listed in the manuscript.

## 16. Reproducibility checklist

To fully reproduce the baseline paper figures from raw outputs, the required
artifacts are:

1. HPC run history:
   - `output/run_history/baseline_sims_screen/51509956/job_config.json`
   - associated run tables if not embedded in the JSON.
2. HPC staged outputs:
   - `output/slurm_output/baseline_sims_screen/51509956/task-*`
3. Consolidated outputs:
   - `output/slurm_output/baseline_sims_screen/51509956/consolidated/*.csv`
4. Source code at the same revision used to generate the final plots.
5. Rendering dependencies:
   - `ggplot2`
   - `patchwork`
   - `viridisLite`
   - `readr`
   - `dplyr`
   - `tidyr`
   - local `susine` and `test_susine` source packages.

The minimal artifact request for auditing final numbers is the consolidated
directory only. The staged outputs are only needed if re-aggregation itself is
being audited.
