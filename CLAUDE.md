# SuSiNE Project — test_susine Simulation Harness

## Project Overview

SuSiNE (Sum of Single Non-central Effects) extends SuSiE's Bayesian fine-mapping
model by allowing nonzero prior means on effect sizes: mu_0 = c * a, where a is a
signed annotation vector (e.g., from Borzoi) and c is an unknown global scale treated
as an exploration knob. Rather than estimating c via EB, we fit a grid of (c, sigma_0^2)
values and aggregate via ELBO-softmax BMA. Target journal: Bioinformatics.
Author: Michael Callahan (mgc5166@psu.edu, Penn State). Advisor: Xiang Zhu (Calico).

## Key Context Documents (load on demand)

@refs/susine_study_plan_high_level.md
@refs/susie_susine_background.md

## Multi-Repo Map

All repos live under `../` relative to this one:

| Repo | Role | Language |
|------|------|----------|
| **test_susine** (this repo) | Simulation harness: run-control grids, SLURM jobs, metrics, figures | R package |
| **susine** (`../susine`) | Core model: susine(), susine_ss(), susine_rss(), SER, ERSS, prior updates, aggregation | R package |
| **susieR** (`../susieR`) | Upstream SuSiE implementation (forked from stephenslab). Suggests dependency | R package |
| **eQTL_annotations_for_susine** (`../eQTL_annotations_for_susine`) | Borzoi annotation extraction and data prep for real-data studies | Python |

Dependency chain: susieR (upstream) -> susine (extends it) -> test_susine (imports susine).

## Build / Test / Load Commands

```r
# Load for interactive development (order matters):
devtools::load_all("../susine")   # susine first
devtools::load_all(".")           # then test_susine

# Regenerate docs and check:
devtools::document(".")
devtools::check(".")
devtools::check("../susine")

# Run susine tests:
testthat::test_dir("../susine/tests/testthat")

# Lint (if .lintr exists):
lintr::lint_package(".")
```

## Key File Index — test_susine/R/

- **run_controls.R**: Job config builders (make_job_config, prior_quality_grid, make_run_tables, assign_task_ids_by_bundle, write_job_artifacts, render_slurm_script). Largest/most complex file.
- **run_model.R**: Task execution (run_task, execute_dataset_bundle, execute_single_run, run_use_case). Handles restarts, model averaging, annealing, aggregation, multimodal metrics.
- **simulate_data.R**: Data generation (generate_simulation_data, simulate_effect_sizes, simulate_phenotype, simulate_priors, standardize_x).
- **evaluation_metrics.R**: Model evaluation (evaluate_model, get_credible_set, cs_purity_min_abs, auprc_average_precision). Coverage, power, FDR, AUPRC.
- **use_cases.R**: use_case_catalog() defining 14 model configurations (a_i through c_ii). Add new use cases here.
- **collect_results.R**: Post-run aggregation (aggregate_staging_outputs, index_staging_outputs, validate_staging_outputs, build_confusion_matrix).
- **dataset_metrics.R**: LD metrics (M1, M2, ECS, NSI, block structure, spectral), z-score metrics, GAM screening.
- **visualize_results.R**: Plotting (power-FDR curves, AUPRC boxplots, CS property plots, X-metric scatter).
- **real_case_slurm.R**: Real-data SLURM jobs (build_real_case_job_config, run_real_case_task).
- **run_accounting.R**: Human-readable run count summaries.
- **utils.R**: %||% null-coalescing, ensure_dir, timestamp_utc, resolve_flag.
- **utils-pipe.R**: Magrittr pipe re-export (auto-generated, do not edit).
- **x_matrix_metrics.R**: Deprecated stub (moved to dataset_metrics.R).

## Key File Index — susine/R/

- **susine.R**: Main IBSS loop (individual-level data) with annealing, warm starts, prior updates.
- **susine_ss.R, susine_rss.R**: Summary-stat variants.
- **SER.R, SER_ss.R**: Single-effect regression updates.
- **ERSS.R, ERSS_ss.R**: Expected residual sum of squares.
- **initialize.R**: Prior/effect/model initialization (supports init_alpha for warm starts).
- **finalize.R**: Output formatting.
- **update_priors.R, update_priors_ss.R**: Empirical Bayes prior updates (mean + variance).
- **aggregate_fits.R**: Multi-fit aggregation methods (softmax_elbo, mean, max_elbo).
- **matrix_computation.R**: Xb, Xty, column stats.
- **univariate_regression.R**: Marginal regression.
- **E2E_testing.R**: End-to-end test helpers.

## Vignettes (test_susine/vignettes/)

- **simulation pipeline/**: run_control_workbook.Rmd, collect_results_workbook.Rmd, visualize_results_workbook_paper_exhibits.Rmd (+ pilot variants).
- **real data pipeline/**: get_genotype_matrices.Rmd, susine_real_case_study.Rmd, real_data_model_clustering.Rmd.
- **one_off_validations/**: Ad hoc validation notebooks.

## SLURM Entry Points

- **inst/scripts/run_task.R**: SLURM array entry point. Loads local dev packages, parses CLI args, calls run_task().
- **inst/scripts/run_real_case_task.R**: Real-data SLURM entry point.

## R Code Style Conventions

1. Assignment: `<-` (not `=`).
2. Pipe: `%>%` (magrittr), not native `|>`.
3. Tidy eval: `.data$col` pronoun in dplyr verbs.
4. NULL coalescing: `%||%` (defined in utils.R).
5. Naming: `snake_case` for all functions and variables.
6. roxygen: Markdown syntax. `@keywords internal` for non-exported helpers.
7. Section headers: `# Section name -----` (trailing dashes to ~col 70).
8. Data frames: Use `tibble::tibble()`, not base `data.frame()`.
9. Imports: Prefer explicit `dplyr::filter()` in package code.
10. Error handling: `stop()` with informative messages; `tryCatch` for I/O.

## Important Patterns

- **Job config pipeline**: make_job_config() -> write_job_artifacts() -> SLURM submission -> run_task() per array element -> aggregate_staging_outputs() post-hoc.
- **Use-case catalog**: All model variants defined in use_case_catalog(). 14 entries: SuSiE baseline (a_i), restarts (a_i_restart), EB variants (a_ii-a_iv), functional priors (b_i-b_vi), annealing (c_i), model averaging (c_ii). Add new use cases by adding rows.
- **Dataset bundle pattern**: X,y generated once per (matrix, p_star, y_noise, seed), then shared across all use_case runs. Priors cached by annotation settings.
- **Buffered I/O**: verbose_file_output=FALSE streams metrics to task-level flush files (CSV + Parquet) to reduce HPC filesystem pressure.

## Paper Status

- **Rough draft**: refs/SuSiNE - rough draft.tex (dated 12/19/25). Needs major revision per updated study plan.
- **Key claims to preserve**: M1 metric explains 42% of AUPRC variation; directional priors work even at overall R2~0.002; CBX8 case study shows 4 posterior basins and 50% more candidate variants.
- **New content needed**: Ensemble framework (Phase B), failure-mode characterization (cross-cutting pillar), compute-confound pilot (Section 5 of plan), multi-locus real data (Phase D), explicit testing of hypotheses H1-H6.
- **Active investigations**: tmp_check_scen3_minimax.R, tmp_s3_beta_sensitivity.R, tmp_check_scenario4_init.R.
