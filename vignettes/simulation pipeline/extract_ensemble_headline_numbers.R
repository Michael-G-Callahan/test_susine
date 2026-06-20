#!/usr/bin/env Rscript
# ------------------------------------------------------------------------------
# Standalone extraction of ensemble headline numbers (2026-06-20 results triage).
#
# FAST PATH: reads the cached paper plot-data RDS that
# prepare_results_workbook_ensemble_scaling_paper.Rmd already produces, then
# writes the headline CSVs WITHOUT re-running the full prep workbook -- no model
# fits, no 2000-rep bootstrap, no scaling recompute. Every quantity below is a
# pure reshape of objects already inside the RDS, so this runs in seconds.
#
# It writes exactly the same CSVs as Section 6 of the prep workbook. Use the
# workbook section on a full rerun; use this script to refresh the headline CSVs
# on their own.
#
# Usage (from the test_susine repo root on the HPC):
#   Rscript "vignettes/simulation pipeline/extract_ensemble_headline_numbers.R"
# or, in an interactive R session with the working dir at the repo root:
#   source("vignettes/simulation pipeline/extract_ensemble_headline_numbers.R")
#
# Requires only: dplyr, tibble, readr (no devtools::load_all of the package).
# ------------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(dplyr)
  library(tibble)
  library(readr)
})

`%||%` <- function(x, y) if (is.null(x) || length(x) == 0L) y else x

## ---- Config: must match the prep workbook ------------------------------------
job_name      <- "ensemble_scaling_full"
parent_job_id <- "53547760"

# Resolve the output root the same way the workbook does (here::here("output")),
# falling back to <cwd>/output when the 'here' project root is unavailable.
output_root <- tryCatch(
  here::here("output"),
  error = function(e) file.path(getwd(), "output")
)

job_dir       <- file.path(output_root, "slurm_output", job_name, parent_job_id)
plot_data_dir <- file.path(job_dir, "figures", "paper_ensemble_scaling", "plot_data")
plot_data_rds <- file.path(plot_data_dir, "paper_ensemble_plot_data.rds")

if (!file.exists(plot_data_rds)) {
  stop(
    "Cached plot-data RDS not found:\n  ", plot_data_rds, "\n",
    "Run prepare_results_workbook_ensemble_scaling_paper.Rmd once to produce it, ",
    "then this script can refresh the headline CSVs quickly.",
    call. = FALSE
  )
}

pd <- readRDS(plot_data_rds)
cat("Loaded plot-data RDS:", plot_data_rds, "\n")

## ---- Pull cached objects (already computed by the prep workbook) --------------
annotation_r2_focus  <- pd$config$annotation_r2_focus %||% 0.3
baseline_auprc       <- pd$refs$baseline_auprc
truth_warm_vanilla   <- pd$refs$truth_warm_vanilla
baseline_compute_sec <- pd$refs$baseline_compute_sec

heatmap_data            <- pd$data$heatmap_data            %||% tibble()
ccs_by_r2               <- pd$data$ccs_by_r2               %||% tibble()
pr_curve_counts         <- pd$data$pr_curve_counts         %||% tibble()
ccs_compute_coordinates <- pd$diagnostics$ccs_compute_coordinates %||% tibble()

stopifnot(
  is.finite(baseline_auprc),
  nrow(heatmap_data) > 0L,
  nrow(ccs_by_r2) > 0L,
  nrow(pr_curve_counts) > 0L
)

## ---- Minimal inlined helpers (copied from the prep workbook) ------------------
approx_equal <- function(x, y, tol = 1e-8) {
  is.finite(x) & is.finite(y) & abs(x - y) <= tol
}
near_focus <- function(x, focus = annotation_r2_focus) {
  approx_equal(as.numeric(x), focus)
}
format_r2_tag <- function(x) {
  gsub("\\.", "p", format(x, trim = TRUE, scientific = FALSE))
}
pct_of_base <- function(delta) {
  if (is.finite(baseline_auprc) && baseline_auprc != 0) {
    100 * delta / baseline_auprc
  } else {
    NA_real_
  }
}

# Identical to prepare_results_workbook_ensemble_scaling_paper.Rmd so the
# operating points match exactly.
pick_precision_operating_points <- function(curve_counts, target_precision = 0.50) {
  if (!nrow(curve_counts)) return(tibble::tibble())
  curve_counts %>%
    dplyr::group_by(.data$method, .data$annotation_r2) %>%
    dplyr::group_modify(function(.x, .y) {
      candidates <- .x %>%
        dplyr::filter(is.finite(.data$precision), .data$precision >= target_precision)
      if (nrow(candidates) > 0L) {
        candidates %>%
          dplyr::arrange(dplyr::desc(.data$recall), .data$pip_threshold) %>%
          dplyr::slice(1L) %>%
          dplyr::mutate(operating_rule = "max recall with precision >= target")
      } else {
        .x %>%
          dplyr::mutate(precision_abs_diff = abs(.data$precision - target_precision)) %>%
          dplyr::arrange(.data$precision_abs_diff, dplyr::desc(.data$recall)) %>%
          dplyr::slice(1L) %>%
          dplyr::select(-dplyr::any_of("precision_abs_diff")) %>%
          dplyr::mutate(operating_rule = "closest observed precision to target")
      }
    }) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(target_precision = target_precision) %>%
    dplyr::select(
      method, annotation_r2, operating_rule, target_precision,
      pip_threshold, precision, recall, cum_tp, cum_fp, called_variants, total_p
    )
}

r2_focus_tag <- format_r2_tag(annotation_r2_focus)

## ---- (a) AUPRC gap attribution at the central annotation regime ---------------
ccs_focus_anchor <- function(method_id) {
  v <- heatmap_data %>%
    dplyr::filter(.data$spec_name == "C-CS", .data$agg_method == method_id) %>%
    dplyr::pull("AUPRC")
  if (length(v)) v[[1]] else NA_real_
}

a_baseline <- baseline_auprc
a_max_elbo <- ccs_focus_anchor("max_elbo")
a_cluster  <- ccs_focus_anchor("cluster_weight_credible")
a_oracle   <- ccs_focus_anchor("oracle")
a_ceiling  <- 1.0

auprc_gap_attribution <- tibble::tribble(
  ~gap_id,                 ~from_label,            ~to_label,               ~from_auprc, ~to_auprc,  ~interpretation,
  "exploration_maxelbo",   "baseline SuSiE",       "C-CS max-ELBO",         a_baseline,  a_max_elbo, "Exploration benefit you keep when committing to the single best-ELBO fit.",
  "exploration_oracle",    "baseline SuSiE",       "C-CS oracle",           a_baseline,  a_oracle,   "Total exploration value if the best fit per locus could be selected (oracle).",
  "aggregation_modelspec", "C-CS max-ELBO",        "C-CS cluster softmax",  a_max_elbo,  a_cluster,  "Deployable aggregation gain; how far ELBO-softmax cluster mixing pushes past the single best-ELBO fit (model-specification barrier closed).",
  "elbo_suboptimality",    "C-CS max-ELBO",        "C-CS oracle",           a_max_elbo,  a_oracle,   "How far max-ELBO selection falls short of the best-quality fit, i.e. ELBO does not always rank fit quality.",
  "aggregation_shortfall", "C-CS cluster softmax", "C-CS oracle",           a_cluster,   a_oracle,   "Remaining headroom of deployable cluster softmax relative to the oracle.",
  "residual_to_ceiling",   "C-CS oracle",          "perfect (AUPRC = 1)",   a_oracle,    a_ceiling,  "Residual ceiling: LD-induced ambiguity inherent to the data plus basins never explored."
) %>%
  dplyr::mutate(
    annotation_r2         = annotation_r2_focus,
    spec_name             = "C-CS",
    delta_auprc           = .data$to_auprc - .data$from_auprc,
    delta_pct_of_baseline = pct_of_base(.data$delta_auprc)
  ) %>%
  dplyr::select(spec_name, annotation_r2, gap_id, from_label, to_label,
                from_auprc, to_auprc, delta_auprc, delta_pct_of_baseline,
                interpretation)

## ---- (b) C-CS AUPRC vs baseline across annotation quality ---------------------
truth_warm_gain <- truth_warm_vanilla - a_baseline

ccs_auprc_by_annotation_vs_baseline <- ccs_by_r2 %>%
  dplyr::filter(.data$agg_method %in% c("max_elbo", "cluster_weight_credible", "oracle")) %>%
  dplyr::transmute(
    spec_name                = .data$spec_name,
    annotation_r2            = as.numeric(.data$annotation_r2),
    agg_method               = .data$agg_method,
    AUPRC                    = .data$AUPRC,
    baseline_auprc           = a_baseline,
    delta_auprc              = .data$AUPRC - a_baseline,
    delta_pct                = pct_of_base(.data$AUPRC - a_baseline),
    truth_warm_vanilla_auprc = truth_warm_vanilla,
    truth_warm_gain          = truth_warm_gain,
    ensemble_minus_truthwarm = (.data$AUPRC - a_baseline) - truth_warm_gain
  ) %>%
  dplyr::arrange(.data$agg_method, .data$annotation_r2)

## ---- (c) Recall at fixed 75% precision, ensemble vs baseline ------------------
pr_precision75_points <- pick_precision_operating_points(
  pr_curve_counts,
  target_precision = 0.75
)

pr75_baseline <- pr_precision75_points %>%
  dplyr::filter(.data$method == "SuSiE baseline") %>%
  dplyr::select(
    annotation_r2,
    baseline_pip_threshold = .data$pip_threshold,
    baseline_precision     = .data$precision,
    baseline_recall        = .data$recall,
    baseline_tp            = .data$cum_tp,
    baseline_called        = .data$called_variants,
    baseline_total_p       = .data$total_p,
    baseline_rule          = .data$operating_rule
  )

pr_precision75_ensemble_vs_baseline <- pr_precision75_points %>%
  dplyr::filter(.data$method != "SuSiE baseline") %>%
  dplyr::select(
    method, annotation_r2,
    ensemble_pip_threshold = .data$pip_threshold,
    ensemble_precision     = .data$precision,
    ensemble_recall        = .data$recall,
    ensemble_tp            = .data$cum_tp,
    ensemble_called        = .data$called_variants,
    ensemble_total_p       = .data$total_p,
    ensemble_rule          = .data$operating_rule
  ) %>%
  dplyr::inner_join(pr75_baseline, by = "annotation_r2") %>%
  dplyr::mutate(
    target_precision     = 0.75,
    recall_diff          = .data$ensemble_recall - .data$baseline_recall,
    recall_pct_change    = dplyr::if_else(
      .data$baseline_recall > 0,
      100 * (.data$ensemble_recall - .data$baseline_recall) / .data$baseline_recall,
      NA_real_
    ),
    causal_tp_pct_change = dplyr::if_else(
      .data$baseline_tp > 0,
      100 * (.data$ensemble_tp - .data$baseline_tp) / .data$baseline_tp,
      NA_real_
    )
  ) %>%
  dplyr::arrange(.data$method, .data$annotation_r2)

## ---- (d) Compute multiplier: 16-run (4x4) and full (8x8) C-CS vs baseline -----
ccs_compute_multiplier_summary <- if (nrow(ccs_compute_coordinates) > 0L &&
                                      is.finite(baseline_compute_sec)) {
  ccs_compute_coordinates %>%
    dplyr::filter(.data$subscale_label %in% c("4x4", "8x8")) %>%
    dplyr::transmute(
      spec_name,
      annotation_r2,
      subscale_label,
      n_members            = dplyr::recode(.data$subscale_label, "4x4" = 16L, "8x8" = 64L),
      ensemble_wall_sec    = .data$x_model_fit_wall_sec,
      ensemble_compute_min = .data$x_compute_min,
      AUPRC                = .data$y_AUPRC,
      AUPRC_se,
      baseline_wall_sec    = baseline_compute_sec,
      compute_multiplier   = .data$x_model_fit_wall_sec / baseline_compute_sec,
      baseline_auprc       = a_baseline,
      delta_auprc          = .data$y_AUPRC - a_baseline,
      delta_pct            = pct_of_base(.data$y_AUPRC - a_baseline),
      n_datasets_compute
    ) %>%
    dplyr::arrange(.data$n_members)
} else {
  message("Compute multiplier summary skipped: no scaling-derived C-CS compute ",
          "coordinates (or baseline wall time) in the cached RDS.")
  tibble::tibble()
}

## ---- (e) C-CSR durability at the central annotation regime --------------------
ccsr_durability <- heatmap_data %>%
  dplyr::filter(.data$spec_name == "C-CSR") %>%
  dplyr::transmute(
    spec_name,
    annotation_r2,
    agg_method,
    AUPRC,
    baseline_auprc = a_baseline,
    delta_auprc,
    delta_pct,
    is_primary_aggregation = .data$agg_method == "cluster_weight_credible"
  ) %>%
  dplyr::arrange(dplyr::desc(.data$is_primary_aggregation), .data$agg_method)

## ---- (f) Compact abstract roll-up ---------------------------------------------
ccs_cluster_at <- function(r2_val, col) {
  v <- ccs_auprc_by_annotation_vs_baseline %>%
    dplyr::filter(.data$agg_method == "cluster_weight_credible",
                  near_focus(.data$annotation_r2, r2_val)) %>%
    dplyr::pull(col)
  if (length(v)) v[[1]] else NA_real_
}
get_mult <- function(members, col) {
  if (!nrow(ccs_compute_multiplier_summary)) return(NA_real_)
  v <- ccs_compute_multiplier_summary %>%
    dplyr::filter(.data$n_members == members) %>%
    dplyr::pull(col)
  if (length(v)) v[[1]] else NA_real_
}
pr75_change <- function(r2_val) {
  v <- pr_precision75_ensemble_vs_baseline %>%
    dplyr::filter(near_focus(.data$annotation_r2, r2_val)) %>%
    dplyr::pull("recall_pct_change")
  if (length(v)) v[[1]] else NA_real_
}
ccsr_primary <- ccsr_durability %>% dplyr::filter(.data$is_primary_aggregation)

paper_ensemble_headline_numbers <- tibble::tribble(
  ~metric,                                   ~value,                                          ~unit,
  "baseline_susie_auprc",                    a_baseline,                                      "AUPRC",
  "truth_warm_vanilla_auprc",                truth_warm_vanilla,                              "AUPRC",
  "truth_warm_gain_over_baseline",           truth_warm_gain,                                 "AUPRC delta",
  "ccs_cluster_auprc_r2_0",                  ccs_cluster_at(0,   "AUPRC"),                    "AUPRC",
  "ccs_cluster_delta_r2_0",                  ccs_cluster_at(0,   "delta_auprc"),              "AUPRC delta",
  "ccs_cluster_delta_pct_r2_0",              ccs_cluster_at(0,   "delta_pct"),                "percent",
  "ccs_cluster_auprc_r2_0p3",                ccs_cluster_at(0.3, "AUPRC"),                    "AUPRC",
  "ccs_cluster_delta_r2_0p3",                ccs_cluster_at(0.3, "delta_auprc"),              "AUPRC delta",
  "ccs_cluster_delta_pct_r2_0p3",            ccs_cluster_at(0.3, "delta_pct"),                "percent",
  "ccs_cluster_gain_minus_truthwarm_r2_0p3", ccs_cluster_at(0.3, "ensemble_minus_truthwarm"), "AUPRC delta",
  "ccs_cluster_auprc_r2_0p5",                ccs_cluster_at(0.5, "AUPRC"),                    "AUPRC",
  "ccs_cluster_delta_r2_0p5",                ccs_cluster_at(0.5, "delta_auprc"),              "AUPRC delta",
  "ccs_cluster_delta_pct_r2_0p5",            ccs_cluster_at(0.5, "delta_pct"),                "percent",
  "pr75_recall_pct_change_r2_0",             pr75_change(0),                                  "percent",
  "pr75_recall_pct_change_r2_0p3",           pr75_change(0.3),                                "percent",
  "pr75_recall_pct_change_r2_0p5",           pr75_change(0.5),                                "percent",
  "baseline_wall_sec",                       baseline_compute_sec,                            "seconds",
  "ccs_16run_wall_sec",                      get_mult(16L, "ensemble_wall_sec"),              "seconds",
  "ccs_16run_compute_multiplier",            get_mult(16L, "compute_multiplier"),             "x baseline",
  "ccs_16run_auprc",                         get_mult(16L, "AUPRC"),                          "AUPRC",
  "ccs_64run_wall_sec",                      get_mult(64L, "ensemble_wall_sec"),              "seconds",
  "ccs_64run_compute_multiplier",            get_mult(64L, "compute_multiplier"),             "x baseline",
  "ccs_64run_auprc",                         get_mult(64L, "AUPRC"),                          "AUPRC",
  "ccsr_durable_auprc_r2_0p3",               ccsr_primary$AUPRC[1]      %||% NA_real_,        "AUPRC",
  "ccsr_durable_delta_r2_0p3",               ccsr_primary$delta_auprc[1] %||% NA_real_,       "AUPRC delta",
  "ccsr_durable_delta_pct_r2_0p3",           ccsr_primary$delta_pct[1]  %||% NA_real_,        "percent"
)

## ---- Write the delineated CSVs ------------------------------------------------
dir.create(plot_data_dir, recursive = TRUE, showWarnings = FALSE)
headline_csv <- function(df, name) {
  path <- file.path(plot_data_dir, name)
  readr::write_csv(df, path)
  cat("Saved:", path, "(", nrow(df), "rows )\n")
}

headline_csv(auprc_gap_attribution,
             paste0("paper_auprc_gap_attribution_r2_", r2_focus_tag, ".csv"))
headline_csv(ccs_auprc_by_annotation_vs_baseline,
             "paper_ccs_auprc_by_annotation_vs_baseline.csv")
headline_csv(pr_precision75_points,
             "paper_pr_precision75_operating_points.csv")
headline_csv(pr_precision75_ensemble_vs_baseline,
             "paper_pr_precision75_ensemble_vs_baseline.csv")
headline_csv(ccs_compute_multiplier_summary,
             "paper_ccs_compute_multiplier_summary.csv")
headline_csv(ccsr_durability,
             paste0("paper_ccsr_durability_r2_", r2_focus_tag, ".csv"))
headline_csv(paper_ensemble_headline_numbers,
             "paper_ensemble_headline_numbers.csv")

cat("\nHeadline-number extraction complete.\n")
print(paper_ensemble_headline_numbers, n = Inf)
