# Baseline Drift Analysis Technical Trace

Date written: 2026-05-20

June 15, 2026 audit note: the main baseline simulation trace now points to the
rerun job `baseline_sims_screen/53522724`, but this drift-analysis trace still
documents the 5-arm effect-matching artifacts generated from the historical
baseline source `baseline_sims_screen/51509956`. The drift PNGs currently in
`../Writings/plots/drift_analysis/` have May 2026 write times and do not appear
to have been regenerated with the June 15 baseline rerun. If the drift figures
are rerun, their pooled PR AUPRCs should follow the current production
step-average-precision convention (`sum(precision_k * delta_recall_k)`, no
trapezoidal interpolation).

This document traces the high-quality-regime 5-arm effect-matching and drift
analysis built on top of the baseline simulation study. It is source-first and
implementation-facing. The goal is to make clear how the final drift-analysis
paper figures were intended to be generated, which cached artifacts exist
locally, and which figure inputs were generated on HPC but are not present in
this checkout.

## 1. Scope, provenance, and caveats

Main files traced:

- Shared 5-arm configuration:
  - `vignettes/one_off_validations/_5arm_cfg.R`
- Refit export:
  - `vignettes/one_off_validations/effect_matching_refit_export_5arm.Rmd`
- Heavy drift analysis:
  - `vignettes/one_off_validations/effect_matching_drift_analysis_5arm.Rmd`
- Final cached plotting workbooks:
  - `vignettes/one_off_validations/effect_matching_pr_pve_quality_figure_5arm.Rmd`
  - `vignettes/one_off_validations/effect_matching_drift_paper_figure_5arm.Rmd`
- Final copied paper plots:
  - `../Writings/plots/drift_analysis/paper_figure_pr_pve_effect_quality.png`
  - `../Writings/plots/drift_analysis/paper_figure_truth_drift_basis_divergence.png`
- Local study directory:
  - `output/effect_matching_study/baseline_sims_5arm_effect_matching/`
- Relevant helpers:
  - `R/run_model.R`
  - `R/evaluation_metrics.R`
  - `R/heritability.R`
  - `R/simulate_data.R`
  - `R/run_controls.R`

Important caveats:

- The local study directory contains the main 5-arm refit records and several
  drift-derived CSVs, but it does not contain all cached inputs expected by the
  final lightweight plotting workbooks.
- Missing locally at trace time:
  - `selected_5arm_settings.csv`
  - `pip_records_5arm.rds`
  - `pr_curves_pooled_top8_5arm.csv`
  - `pr_auc_pooled_top8_5arm.csv`
  - `effect_quality_density_shift_5arm.csv`
  - `causal_truth_5arm.rds`
  - `truth_fitted_y_5arm_top23.rds`
  - `truth_drift_summary_5arm_top23.csv`
  - `truth_drift_assignments_5arm_top23.csv`
- The local `drift_summary_5arm.csv` has `weighted_chord_drift`, but not the
  newer `weighted_r2_drift` column used by the final paper-figure workbook.
  Therefore the local drift cache appears older than the May 20 final paper
  PNGs.
- The final paper PNGs in `../Writings/plots/drift_analysis` are newer than the
  local 5-arm derived CSVs and appear to have been rendered from HPC or synced
  caches not present locally.
- This trace distinguishes source-code provenance from locally available
  artifact provenance.

## 2. Scientific purpose

The drift analysis asks how the signed prior-mean channel (`mu_0`) and the
unsigned functional-pi channel change a SuSiE-family fit at a mechanistic
level.

The 5-arm chain is designed to separate three possible mechanisms:

- within-basin PIP reweighting: annotated fit improves variant ranking, but a
  vanilla refit initialized from the annotated solution loses the gain;
- durable basin movement: annotated fit moves IBSS into a new basin, and a
  vanilla refit initialized from that basin keeps the gain after the prior is
  removed;
- unstable prior-dependent movement: annotated fit moves to a different basin,
  but vanilla snaps back or degrades when the prior is removed.

The manuscript places this analysis in the baseline-results section, after the
main baseline AUPRC comparison, focusing on the high-quality annotation regime:

- `phi_a = 0.5`
- `nu_a = 0.9`

## 3. Shared 5-arm configuration

The shared configuration lives in `_5arm_cfg.R`.

Study-level constants:

- `study_name_5arm <- "baseline_sims_5arm_effect_matching"`
- `n_matrices_target_5arm <- 150L`
- `seeds_5arm <- 1:4`
- `L_value_5arm <- 10L`
- `y_noise_value_5arm <- NA_real_`
- `p_star_value_5arm <- NA_integer_`
- `architecture_value_5arm <- "susie2_oligogenic"`
- `data_scenario_5arm <- "scenario_1"`
- `max_iter_5arm <- 100L`
- `credible_set_rho_5arm <- 0.95`
- `purity_threshold_5arm <- 0.50`

Baseline dependency:

- `baseline_job_name_5arm <- "baseline_sims_screen"`
- `baseline_parent_job_id_5arm <- "51509956"`

Default/manual high-quality treatment settings in `_5arm_cfg.R`:

- `annotation_r2_value_5arm <- 0.5`
- `inflate_match_value_5arm <- 0.9`
- `c_value_mu_5arm <- 0.6`
- `tau_value_pi_5arm <- 0.464`
- `sigma_0_2_vanilla_value_5arm <- 0.2`
- `sigma_0_2_mu_value_5arm <- 0.01`
- `sigma_0_2_pi_value_5arm <- 0.2`

However, the local generated `selected_run_rows.csv` for the 5-arm study uses
these three cold settings:

- `susie_vanilla`: `no_annotation`;
  `susie_vanilla | sigma=0.01`.
- `susine_functional_mu`: `r2=0.5 | inflate=0.9`;
  `susine_functional_mu | c=0.6 | sigma=0.01`.
- `susie_functional_pi`: `r2=0.5 | inflate=0.9`;
  `susie_functional_pi | tau=0.464 | sigma=0.01`.

This local generated state matches the manuscript drift caption's statement
that all five arms use `sigma_0_2 = 0.01`. It does not match the baseline
top-settings CSV's best annotation-agnostic SuSiE row, which is
`susie_vanilla | sigma=0.1`. The drift trace should therefore treat the
5-arm selected-run CSV and manuscript drift caption as the source of truth for
the drift study, not the baseline top-settings row for vanilla.

The five method labels are:

| Role | Method label |
|---|---|
| cold baseline | `susie_vanilla` |
| cold signed prior-mean treatment | `susine_functional_mu` |
| warm vanilla from signed treatment | `susie_vanilla_warm_from_mu` |
| cold unsigned functional-pi treatment | `susie_functional_pi` |
| warm vanilla from pi treatment | `susie_vanilla_warm_from_pi` |

## 4. How 5-arm settings are resolved

`_5arm_cfg.R` defines two paths.

The preferred path:

1. read the baseline consolidated `screening_summary.csv`;
2. pick the best row for:
   - `susie_vanilla`, annotation-agnostic;
   - `susine_functional_mu` at `annotation_r2 = 0.5`, `inflate_match = 0.9`;
   - `susie_functional_pi` at `annotation_r2 = 0.5`, `inflate_match = 0.9`;
3. cache the result to `selected_5arm_settings.csv`.

The fallback path:

1. use the manually declared constants in `_5arm_cfg.R`;
2. build a fresh job config if the baseline job config is unavailable.

The local generated `selected_run_rows.csv` indicates that the run actually
used a third effective state:

- high-quality treatment settings `c = 0.6`, `tau = 0.464`;
- `sigma_0_2 = 0.01` for all cold arms;
- 600 dataset bundles;
- 3 cold selected run rows per dataset.

Because `selected_5arm_settings.csv` is not local, the exact setting-resolution
path for the final run cannot be verified locally from the cache alone.

## 5. Refit export workflow

The export workbook is:

- `vignettes/one_off_validations/effect_matching_refit_export_5arm.Rmd`

It writes to:

- `output/effect_matching_study/baseline_sims_5arm_effect_matching/`

The intended files are:

- `effect_refit_records.rds`
- `effect_refit_metadata.csv`
- `pip_records_5arm.rds`
- `selected_run_rows.csv`
- `dataset_bundles.csv`
- `refit_progress_log.csv`
- `selected_5arm_settings.csv`

The local directory contains:

| File | Size | Last write time |
|---|---:|---|
| `effect_refit_records.rds` | 153,654,516 | 2026-04-27 11:41:40 |
| `effect_refit_metadata.csv` | 89,883,337 | 2026-04-27 11:41:49 |
| `selected_run_rows.csv` | 586,915 | 2026-04-27 11:43:00 |
| `dataset_bundles.csv` | 23,541 | 2026-04-27 11:43:00 |
| `refit_progress_log.csv` | 52,439 | 2026-04-27 11:41:49 |

`pip_records_5arm.rds` and `selected_5arm_settings.csv` are not present in the
local directory.

The local `dataset_bundles.csv` has 600 rows. The local
`effect_refit_metadata.csv` has 6,000 effect rows for each of the five methods:

| Method | Effect rows | Datasets |
|---|---:|---:|
| `susie_vanilla` | 6,000 | 600 |
| `susine_functional_mu` | 6,000 | 600 |
| `susie_vanilla_warm_from_mu` | 6,000 | 600 |
| `susie_functional_pi` | 6,000 | 600 |
| `susie_vanilla_warm_from_pi` | 6,000 | 600 |

This is exactly `600 datasets x L=10 effects` per arm.

## 6. Cold fit order

For each dataset, the refit export workbook fits arms in this order:

1. cold `susie_vanilla`;
2. cold `susine_functional_mu`;
3. `susie_vanilla_warm_from_mu`;
4. cold `susie_functional_pi`;
5. `susie_vanilla_warm_from_pi`.

This order is intentional. The warm refit is run immediately after its source
annotated fit so the source fit object is still in memory.

Cold fits are run through:

- `fit_one_method_with_fit()`;
- `test_susine:::lookup_use_case()`;
- `test_susine:::run_use_case()`;
- `test_susine:::evaluate_model()`.

The fit object is retained for warm-start construction.

## 7. Warm-refit mechanics

Warm refits are implemented by `fit_warm_refit()` in the export workbook.

The source fit is converted to `init_effect_fits` by
`build_init_effect_fits_robust()`, which handles three schemas:

- local `susine` fits with `effect_fits$alpha`, `b_hat`, and `b_2_hat`;
- wrapped `susieR` fits with alpha in `effect_fits` and `mu`, `mu2` in `raw`;
- raw `susieR` fits with top-level `alpha`, `mu`, and `mu2`.

The warm fit calls the local `susine` function with:

- `L = 10`;
- same `X` and `y`;
- `mu_0 = 0`;
- `sigma_0_2 = sigma_0_2_vanilla_value`;
- uniform inclusion weights `rep(1 / p, p)`;
- `prior_update_method = "none"`;
- `max_iter = 100`;
- `init_effect_fits` from the source annotated fit.

The warm run row is synthetic:

- `use_case_id` is reset to `susie_vanilla`;
- `method_family` is set to either `susie_vanilla_warm_from_mu` or
  `susie_vanilla_warm_from_pi`;
- annotation fields are set to `NA`;
- `setting_label` records the warm source and vanilla sigma.

This design removes the annotation prior while preserving the starting basin.

## 8. Exported effect records

Each effect row records:

- dataset and matrix identifiers;
- model method;
- run id;
- effect index;
- setting and annotation labels;
- credible-set size, purity, and coverage;
- effect diffuseness:
  - `effect_pip_entropy`;
  - `effect_pip_entropy_core95`;
  - `effect_k_eff_signal`;
  - `effect_k_eff_signal_core95`;
- effect accuracy:
  - `accuracy_ratio`;
- per-effect fitted-y vector:
  - `fitted_y_effect`;
- per-effect explained-variance fraction:
  - `var_y_explained_frac`;
- model-level AUPRC and TPR@FPR=0.05 copied onto the row.

The model-level metrics copied onto each row come from `evaluate_model()`, which
also emits two local-genetic-variance (PVE/heritability) estimands: the legacy
`hg2` (still emitted) = `var(fitted_y) / var(y)`, clipped to `[0, 1]`, via
`estimate_hg2()` in `R/evaluation_metrics.R`; and the corrected
local-genetic-variance decomposition (the now-reported estimand) via
`hg2_components()` with its engine in `R/heritability.R`, emitting
`hg2_postmean = var(E[Xβ|y]) / var(y)`, `hg2_uncertainty` (a within-fit
posterior-variance correction), and
`hg2_expected_pve = hg2_postmean + hg2_uncertainty = E[var(Xβ)|y] / var(y)` for a
single fit (an additional `hg2_between_fit` term appears only under
multi-fit/ensemble aggregation). Provenance: commits `4e0721b` / `236820b`
(2026-06-09). Note this model-level PVE is distinct from the per-effect
`var_y_explained_frac` above.

The export workbook reconstructs per-effect fitted-y vectors from the effect
coefficient matrix and the original `X`. It explicitly checks reconstructed
total fitted-y against `model_fit$fitted_y` when available, after centering and
scale correction.

## 9. Drift-analysis workflow

The heavy drift notebook is:

- `vignettes/one_off_validations/effect_matching_drift_analysis_5arm.Rmd`

It reads:

- `effect_refit_records.rds`
- `selected_5arm_settings.csv` in the intended workflow

In the local cache, `selected_5arm_settings.csv` is absent. The existing derived
CSV outputs therefore likely came from an earlier run where the settings file
was present or the local cache was later partially cleaned.

The notebook defines these method pairs:

| Pair id | Chain | A | B |
|---|---|---|---|
| `vanilla__mu` | `mu` | `susie_vanilla` | `susine_functional_mu` |
| `mu__vanilla_warm_from_mu` | `mu` | `susine_functional_mu` | `susie_vanilla_warm_from_mu` |
| `vanilla__vanilla_warm_from_mu` | `mu` | `susie_vanilla` | `susie_vanilla_warm_from_mu` |
| `vanilla__pi` | `pi` | `susie_vanilla` | `susie_functional_pi` |
| `pi__vanilla_warm_from_pi` | `pi` | `susie_functional_pi` | `susie_vanilla_warm_from_pi` |
| `vanilla__vanilla_warm_from_pi` | `pi` | `susie_vanilla` | `susie_vanilla_warm_from_pi` |

## 10. Effect drift definition

For one dataset and one pair of fits A and B:

1. Stack effect-specific fitted-y vectors into an `n x L` matrix for A and for
   B.
2. Compute the `L x L` correlation matrix between A effects and B effects.
3. Effects with near-zero standard deviation get correlation 0.
4. If any correlations are negative, shift the full matrix before assignment
   because `clue::solve_LSAP(maximum = TRUE)` requires nonnegative weights.
5. Solve the maximum-weight one-to-one assignment with
   `clue::solve_LSAP()`.
6. Record each assigned edge correlation and each side's
   `var_y_explained_frac`.

The summary function records:

- `n_edges`;
- `mean_corr`;
- `median_corr`;
- `min_corr`;
- `n_below_cutoff`, with cutoff `0.9`;
- `n_active_pairs`, with active threshold `var_y_explained_frac > 0.005` on
  both sides;
- `n_active_below_cutoff`;
- `mean_corr_active`;
- `weighted_drift = 1 - weighted mean corr`;
- `weighted_chord_drift = weighted mean sqrt(2 * (1 - corr))`;
- in newer source, `weighted_r2_drift = weighted mean (1 - corr^2)`.

Weights are `min(var_y_explained_A, var_y_explained_B)`, so inactive matched
effects contribute little.

The local `drift_summary_5arm.csv` has 600 rows per pair. Its mean chord drift
by pair is:

| Pair id | Mean weighted chord drift | Approx. median |
|---|---:|---:|
| `vanilla__mu` | 0.063451 | 0.033676 |
| `mu__vanilla_warm_from_mu` | 0.043903 | 0.036383 |
| `vanilla__vanilla_warm_from_mu` | 0.035749 | 0.012414 |
| `vanilla__pi` | 0.034536 | 0.016945 |
| `pi__vanilla_warm_from_pi` | 0.028730 | 0.023050 |
| `vanilla__vanilla_warm_from_pi` | 0.025969 | 0.013121 |

This local summary supports the qualitative claim that the mu chain moves the
effect basis more than the pi chain, but it is not the exact cache used by the
final May 20 paper figure because it lacks `weighted_r2_drift`.

## 11. AUPRC trajectory

The heavy drift notebook builds `auprc_per_dataset_5arm.csv` from the model
AUPRC values stored in the effect records.

It computes:

- `auprc_delta_mu = auprc_susine_functional_mu - auprc_susie_vanilla`;
- `auprc_delta_pi = auprc_susie_functional_pi - auprc_susie_vanilla`;
- `auprc_delta_warm_mu_vs_vanilla`;
- `auprc_delta_warm_pi_vs_vanilla`;
- `auprc_round_trip_mu = warm_from_mu - cold_mu`;
- `auprc_round_trip_pi = warm_from_pi - cold_pi`.

The local `auprc_per_dataset_5arm.csv` has 600 rows. Mean per-dataset AUPRCs
from the local file are:

| Method | Mean per-dataset AUPRC |
|---|---:|
| cold `susie_vanilla` | 0.146921 |
| cold `susine_functional_mu` | 0.184921 |
| warm `susie_vanilla_warm_from_mu` | 0.146644 |
| cold `susie_functional_pi` | 0.156264 |
| warm `susie_vanilla_warm_from_pi` | 0.145353 |

These local per-dataset means should not be confused with pooled PR-curve
AUPRCs in the final paper figure. The final figure legend uses pooled PR AUPRC
from `pr_auc_pooled_top8_5arm.csv`, which is not local.

## 12. PIP JSD

The heavy notebook computes model-wide PIP Jensen-Shannon divergence from
`pip_records_5arm.rds`.

The definition mirrors `R/run_model.R::js_distance`:

- treat raw model-wide PIP vectors as categorical mass vectors after adding a
  small epsilon;
- normalize each vector;
- compute Jensen-Shannon divergence through the midpoint distribution;
- report a scalar distance/divergence for each pair.

The intended output is:

- `pip_drift_summary_5arm.csv`

The local `pip_drift_summary_5arm.csv` exists, but all 3,600 `pip_jsd` values
are `NA`. This is consistent with the missing local `pip_records_5arm.rds` and
indicates that the local PIP-JSD cache is not usable for the final paper figure.
The final May 20 figure that includes PIP-JSD panels must have used a different
or regenerated cache.

## 13. Effect-quality density shifts

The heavy notebook builds pair rows for effect-quality density shifts:

- output:
  `effect_quality_pair_rows_5arm.csv`

The effect-quality axes are:

- x-axis:
  `effect_pip_entropy_core95`, displayed as
  `K_l = exp(H(alpha_l^(95)))`;
- y-axis:
  `accuracy_ratio = max causal alpha / max any alpha`;
- weight:
  `pve_weight = max(var_y_explained_frac, 0)`.

The pair-row table stacks source and destination effects for these comparisons:

- cold vanilla to annotated treatment;
- annotated treatment to warm vanilla;
- cold vanilla to warm vanilla;
- separately for the mu and pi chains.

The density-shift grid is produced by `weighted_kde2d_df()`:

- resolve x/y plotting limits from trimmed quantiles;
- build a fixed grid with `density_grid_n = 140`;
- use bandwidth `c(0.45, 0.10)`;
- estimate source and destination PVE-weighted 2D KDEs;
- plot `destination - source`.

The local directory has `effect_quality_pair_rows_5arm.csv`, but not
`effect_quality_density_shift_5arm.csv`. The final lightweight plotting
workbooks can rebuild the density-shift grid from pair rows if needed.

## 14. Top-8 PR scoring

The heavy notebook regenerates causal truth from `dataset_bundles.csv` by
calling `generate_data_for_bundle()` for each dataset.

For pooled PR curves:

- `pr_top_k <- 8L`;
- causals are ranked by `abs(beta)` within each dataset;
- the top 8 are positives;
- true non-causal variants are negatives; non-top causal variants (true causals
  outside the top 8) are DROPPED entirely -- neither positive nor negative.

This matches the baseline-sims confusion-bin AUPRC convention in
`compute_confusion_bins()` (`R/run_model.R`), which does
`scored <- as.integer(seq_len(n) %in% causal_mask)` then
`keep <- !(causal == 1L & scored == 0L)` so weak causals are removed and AUPRC
is not penalized for them.

Outputs intended by the heavy notebook:

- `causal_truth_5arm.rds`
- `pr_curves_pooled_top8_5arm.csv`
- `pr_auc_pooled_top8_5arm.csv`
- `drift_analysis_plots/pr_curves_pooled_top8.png`

These PR caches are not present locally. The final paper PR plot therefore
cannot be locally audited from cached CSVs without syncing those files.

## 15. Truth-basin drift

The later section of `effect_matching_drift_analysis_5arm.Rmd` computes drift
to a synthetic truth basis.

Purpose:

- compare each fit's explanatory basis to the fitted-y basis generated by the
  true causal effects;
- separate movement toward truth from movement that is merely different.

The final source setting is:

- `truth_drift_top_k <- 23L`

This means the synthetic truth side uses all 23 true causals in the
oligogenic architecture. The fit side has `L = 10`, so the code pads rows to
`target_dim = max(truth_drift_top_k, records_L)`. The extra truth components
can be matched to zero-padded fit components, but their weights do not dominate
because the weighted drift uses the minimum explained variance on the two sides.

Intended outputs:

- `truth_fitted_y_5arm_top23.rds`
- `truth_drift_assignments_5arm_top23.csv`
- `truth_drift_summary_5arm_top23.csv`

These are not present locally. The final copied plot
`paper_figure_truth_drift_basis_divergence.png` depends on this newer truth
drift cache.

## 16. Heavy drift diagnostic plots

The heavy drift notebook writes diagnostic plots under:

- `output/effect_matching_study/baseline_sims_5arm_effect_matching/drift_analysis_plots/`

Local diagnostic files include:

- `drift_profile_histograms.png`
- `chain_auprc_trajectory.png`
- `headline_drift_vs_warm_auprc.png`
- `unpaired_count_vs_auprc_delta.png`
- `mu_vs_pi_drift_comparison.png`
- `basin_and_pip_drift_summary.png`
- `effect_quality_density_matched_unmatched_mu_5arm.png`
- `effect_quality_density_matched_unmatched_pi_5arm.png`

These local diagnostic PNGs were written on 2026-04-27. They are useful for
early analysis but are not the final paper figures copied under
`../Writings/plots/drift_analysis`.

## 17. Final paper plotting workbooks

Two lightweight plotting workbooks build the final drift-analysis paper plots
from cached CSV/RDS inputs.

### 17.1 PR, PVE, and effect-quality figure

Workbook:

- `vignettes/one_off_validations/effect_matching_pr_pve_quality_figure_5arm.Rmd`

Final local copied paper file:

- `../Writings/plots/drift_analysis/paper_figure_pr_pve_effect_quality.png`

The workbook expects:

- `pr_curves_pooled_top8_5arm.csv`
- `pr_auc_pooled_top8_5arm.csv` if available;
- `effect_quality_pair_rows_5arm.csv`;
- `effect_quality_density_shift_5arm.csv`, or it rebuilds that cache from
  pair rows.

It builds:

- panel A: pooled top-8 PR curves for the five arms;
- panel B: per-effect quality density shifts.

The local final copied PNG exists and was written:

- 2026-05-20 15:09:48

The exact PR-curve CSV and PR-AUC CSV used for that figure are not local.

### 17.2 Truth drift and basis-divergence figure

Workbook:

- `vignettes/one_off_validations/effect_matching_drift_paper_figure_5arm.Rmd`

Final local copied paper file:

- `../Writings/plots/drift_analysis/paper_figure_truth_drift_basis_divergence.png`

The workbook expects:

- `drift_summary_5arm.csv`;
- `auprc_per_dataset_5arm.csv`;
- `pip_drift_summary_5arm.csv`;
- `pr_curves_pooled_top8_5arm.csv`;
- `pr_auc_pooled_top8_5arm.csv` if available;
- `effect_quality_pair_rows_5arm.csv`;
- `effect_quality_density_shift_5arm.csv`;
- `pip_records_5arm.rds`;
- `causal_truth_5arm.rds`;
- truth-drift caches when building truth-basis panels.

The source uses `weighted_r2_drift` for the paper-facing basin-drift panels,
not the older local `weighted_chord_drift` column.

The local final copied PNG exists and was written:

- 2026-05-20 15:09:51

The exact truth-drift and PIP-JSD caches used for that figure are not local.

## 18. Final copied drift artifacts

Local final copied paper files under `../Writings/plots/drift_analysis`:

| File | Size | Last write time |
|---|---:|---|
| `paper_figure_pr_pve_effect_quality.png` | 797,995 | 2026-05-20 15:09:48 |
| `paper_figure_truth_drift_basis_divergence.png` | 1,167,756 | 2026-05-20 15:09:51 |

These timestamps are substantially newer than the local derived CSVs in
`output/effect_matching_study/baseline_sims_5arm_effect_matching`, whose visible
derived outputs were last written on 2026-04-27.

## 19. Source-code consistency notes

### 19.1 Sigma choice in the 5-arm drift study

There are three relevant sources:

- `_5arm_cfg.R` manual fallback defaults:
  - vanilla `sigma_0_2 = 0.2`;
  - mu `sigma_0_2 = 0.01`;
  - pi `sigma_0_2 = 0.2`.
- baseline paper `top_settings_by_family.csv`:
  - best `susie_vanilla` row is `sigma=0.1`;
  - best high-quality mu and pi rows use `sigma=0.01`.
- local 5-arm `selected_run_rows.csv` and manuscript drift caption:
  - all five arms use `sigma_0_2 = 0.01`.

For the actual drift analysis, use the generated selected-run rows and the
manuscript caption as the operational source: all five arms are treated as
`sigma_0_2 = 0.01`.

### 19.2 Top-8 positives versus polygenic causals

Both baseline confusion bins and the 5-arm PR regeneration code use top-8
causal variants by `abs(beta)` as positives. In source code, true non-causal
variants are the negatives, and variants that are truly causal but outside the
top-8 set (including members of the 15-variant polygenic causal tier) are
DROPPED from the table entirely -- neither positive nor negative -- so AUPRC is
not penalized for them.

Manuscript prose that says polygenic-tier causals are excluded from negatives is
therefore correct: those weak causals are dropped, not scored as negatives.

### 19.3 Chord drift versus `1 - r^2` drift

The heavy drift source now defines both:

- chord drift: `sqrt(2 * (1 - corr))`;
- basis drift: `1 - corr^2`.

The local `drift_summary_5arm.csv` only contains chord drift. The final
paper-figure source expects `weighted_r2_drift`, so the final paper figure was
not rendered from this exact local `drift_summary_5arm.csv`.

### 19.4 PIP JSD local cache

The local `pip_drift_summary_5arm.csv` has all `pip_jsd` values as `NA`. The
source paper figure expects finite PIP-JSD values. Therefore any final figure
panel using PIP JSD must have been rendered from a cache that is not present
locally.

## 20. Minimal artifact requests for full audit

To audit the final drift figures without downloading the full HPC tree, request
only these files from the HPC study directory:

- `selected_5arm_settings.csv`
- `pip_records_5arm.rds`
- `pr_curves_pooled_top8_5arm.csv`
- `pr_auc_pooled_top8_5arm.csv`
- `effect_quality_density_shift_5arm.csv`
- `causal_truth_5arm.rds`
- `truth_fitted_y_5arm_top23.rds`
- `truth_drift_summary_5arm_top23.csv`
- `truth_drift_assignments_5arm_top23.csv`
- the version of `drift_summary_5arm.csv` containing `weighted_r2_drift`
- the version of `pip_drift_summary_5arm.csv` containing finite `pip_jsd`

The heavy `effect_refit_records.rds` is already local and does not need to be
requested again unless cache regeneration is required.

## 21. Reproducibility checklist

To rerun the 5-arm drift workflow from scratch:

1. Restore or rerun the baseline consolidated outputs for
   `baseline_sims_screen/51509956`.
2. Confirm the intended 5-arm selected settings, especially the all-arm
   `sigma_0_2 = 0.01` convention.
3. Render `effect_matching_refit_export_5arm.Rmd` to regenerate effect records
   and PIP records.
4. Render `effect_matching_drift_analysis_5arm.Rmd` to regenerate:
   - drift assignments and summaries;
   - AUPRC trajectories;
   - PIP JSD;
   - effect-quality density shifts;
   - top-8 PR curves and AUPRCs;
   - top-23 truth-basin drift.
5. Render the two lightweight paper plotting workbooks.
6. Copy final PNGs to `../Writings/plots/drift_analysis/`.
