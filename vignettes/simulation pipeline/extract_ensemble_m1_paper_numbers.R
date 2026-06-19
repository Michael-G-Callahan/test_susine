#!/usr/bin/env Rscript

# =============================================================================
# Ensemble-scaling headline numbers + paired-bootstrap uncertainty (M1 writeup).
#
# This script regenerates the manuscript-facing C-CS ensemble headline for the
# Results sentence in `Writings/second draft/susine_second_draft.tex`:
#
#     C-CS with cluster-weight aggregation, phi_a = 0.3, vs the zero-prior-mean
#     SuSiE-equivalent baseline: Delta AUPRC = +0.066 (point) and the C-CS
#     truth-aware oracle = +0.089.
#
# It MUST match the figure generator, which is the paper-prep workbook
#   vignettes/simulation pipeline/prepare_results_workbook_ensemble_scaling_paper.Rmd
# That workbook defines the authoritative slice:
#   * job:          ensemble_scaling_full / 53547760   (NOT 51250228; see below)
#   * ensemble:     spec_name == "C-CS"
#   * annotation:   annotation_r2 (phi_a) == 0.3
#   * aggregation:  agg_method == "cluster_weight_credible"  (credible-shift,
#                   complete-linkage cut 0.05, Method B). This is the manuscript's
#                   "cluster softmax" primary rule.
#   * baseline:     spec_name == "baseline-single", use_case_id == "susine_vanilla",
#                   is.na(agg_method)  -- the job's own cold SuSiE-equivalent fit.
#   * scoring:      pooled AUPRC over top-8 causal confusion bins
#                   (auprc_from_pooled_bins), pooling per-dataset bucket counts.
#   * oracle:       best individual C-CS member per dataset by AUPRC, then pool
#                   (compute_ccs_resolution_oracle, "8x8" resolution = full grid).
#
# -----------------------------------------------------------------------------
# OBSOLETE PATHS -- DO NOT REUSE:
#   * The old bootstrap artifact targeting `agg_method == "cluster_weight"`
#     (legacy JSD-0.15 + 1/frequency, "Method C") at `annotation_r2 = 0.2`
#     with observed_delta ~= 0.0585 is STALE. It uses a different aggregator AND
#     a different annotation slice than the current figure. The ~+0.059 value in
#     the older technical trace is that legacy Method-C number, not the paper's
#     Method-B `cluster_weight_credible` value (+0.066).
#   * The default parent_job_id was 51250228 (an old exploratory job, also used
#     by the deprecated visualize_results_workbook_ensemble_scaling.Rmd). The
#     paper figures use 53547760. The default below is now 53547760, and the
#     script fails loudly if the requested slice/method is absent rather than
#     silently substituting the legacy method.
# =============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(tibble)
  library(tidyr)
})

`%||%` <- function(x, y) {
  if (is.null(x) || length(x) == 0L) y else x
}

parse_args <- function() {
  args <- commandArgs(trailingOnly = TRUE)
  out <- list(
    job_name = "ensemble_scaling_full",
    # Paper figure job (prepare/visualize paper workbooks use 53547760).
    parent_job_id = "53547760",
    output_root = NULL,
    bootstrap_B = 2000L,
    bootstrap_seed = 20260505L,
    annotation_r2_focus = 0.3,
    # Current primary cluster-weight method (Method B, credible-shift cut 0.05).
    primary_agg_method = "cluster_weight_credible",
    # Legacy method name kept ONLY for detection/auditing; never used to score.
    legacy_agg_method = "cluster_weight",
    # Hard-fail if the recomputed (bins) delta disagrees with the precomputed
    # figure delta (auprc_pooled_agg_by_r2.csv) by more than this. A large gap
    # means a wrong slice/method, not rounding.
    figure_match_tol = 5e-3,
    # Allow running off the phi_a = 0.3 focus on purpose (debugging only).
    allow_nonfocus = FALSE,
    out_dir = NULL
  )
  for (arg in args) {
    if (!grepl("^--", arg)) next
    kv <- strsplit(sub("^--", "", arg), "=", fixed = TRUE)[[1]]
    key <- kv[[1]]
    val <- if (length(kv) > 1L) paste(kv[-1L], collapse = "=") else TRUE
    if (!key %in% names(out)) stop("Unknown argument: --", key, call. = FALSE)
    out[[key]] <- val
  }
  # Interactive (RStudio) override: define a list `m1_args_override` in the
  # global environment BEFORE source()-ing this file to set any argument above
  # without command-line flags, e.g.
  #   m1_args_override <- list(output_root = "/storage/.../output")
  if (exists("m1_args_override", envir = globalenv(), inherits = FALSE)) {
    ov <- get("m1_args_override", envir = globalenv())
    for (key in names(ov)) {
      if (!key %in% names(out)) stop("Unknown m1_args_override key: ", key, call. = FALSE)
      out[[key]] <- ov[[key]]
    }
  }
  out$bootstrap_B <- as.integer(out$bootstrap_B)
  out$bootstrap_seed <- as.integer(out$bootstrap_seed)
  out$annotation_r2_focus <- as.numeric(out$annotation_r2_focus)
  out$figure_match_tol <- as.numeric(out$figure_match_tol)
  out$allow_nonfocus <- isTRUE(as.logical(out$allow_nonfocus))
  out
}

script_path <- function() {
  file_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  if (length(file_arg)) {
    path <- sub("^--file=", "", file_arg[[1]])
    path <- gsub("~\\+~", " ", path, fixed = FALSE)
    if (file.exists(path)) {
      return(normalizePath(path, winslash = "/", mustWork = TRUE))
    }
  }
  if (!is.null(sys.frames()[[1]]$ofile) && file.exists(sys.frames()[[1]]$ofile)) {
    return(normalizePath(sys.frames()[[1]]$ofile, winslash = "/", mustWork = TRUE))
  }
  NA_character_
}

repo_root_from_script <- function() {
  sp <- script_path()
  if (!is.na(sp)) return(normalizePath(file.path(dirname(sp), "..", ".."), winslash = "/", mustWork = TRUE))
  normalizePath(getwd(), winslash = "/", mustWork = TRUE)
}

# --- AUPRC helpers (byte-for-byte equivalent to the package internals
#     test_susine:::compute_auprc_single (visualize_results.R) and
#     test_susine:::auprc_from_pooled_bins (collect_results.R), kept inline so
#     this script is standalone and produces the same number the figure does) ---

compute_auprc_single <- function(precision, recall) {
  df <- data.frame(precision = precision, recall = recall) |>
    dplyr::filter(!is.na(.data$precision), !is.na(.data$recall)) |>
    dplyr::arrange(.data$recall)
  if (nrow(df) < 2L) return(NA_real_)
  delta_recall <- c(df$recall[[1]], diff(df$recall))
  sum(df$precision * delta_recall)
}

# Pool bucket counts across datasets (group by pip_threshold) -> pooled AUPRC.
auprc_from_bins <- function(bins) {
  if (is.null(bins) || nrow(bins) == 0L) return(NA_real_)
  b <- bins |>
    dplyr::group_by(.data$pip_threshold) |>
    dplyr::summarise(
      n_causal_at_bucket = sum(.data$n_causal_at_bucket, na.rm = TRUE),
      n_noncausal_at_bucket = sum(.data$n_noncausal_at_bucket, na.rm = TRUE),
      .groups = "drop"
    ) |>
    dplyr::arrange(dplyr::desc(.data$pip_threshold))
  cum_tp <- cumsum(b$n_causal_at_bucket)
  cum_fp <- cumsum(b$n_noncausal_at_bucket)
  total_p <- sum(b$n_causal_at_bucket)
  if (total_p == 0) return(NA_real_)
  denom <- cum_tp + cum_fp
  precision <- ifelse(denom > 0, cum_tp / denom, NA_real_)
  recall <- cum_tp / total_p
  compute_auprc_single(precision, recall)
}

# Weighted pooled AUPRC: each unit's bucket counts scaled by bootstrap multiplicity.
weighted_auprc <- function(bins, unit_weights, unit_col_chr) {
  joined <- bins |>
    dplyr::inner_join(unit_weights, by = unit_col_chr) |>
    dplyr::group_by(.data$pip_threshold) |>
    dplyr::summarise(
      n_causal_at_bucket = sum(.data$n_causal_at_bucket * .data$boot_weight, na.rm = TRUE),
      n_noncausal_at_bucket = sum(.data$n_noncausal_at_bucket * .data$boot_weight, na.rm = TRUE),
      .groups = "drop"
    )
  auprc_from_bins(joined)
}

# Mirrors test_susine:::normalize_agg_method_id: only canonicalizes the legacy
# JSD-0.50 alias. It does NOT touch cluster_weight / cluster_weight_credible.
normalize_agg_method_id <- function(x) {
  x <- as.character(x)
  x[x %in% c("cluster_weight_050", "cluster_weight_jsd_050")] <- "cluster_weight_jsd_050"
  x
}
normalize_agg_method_df <- function(df, col = "agg_method") {
  if (is.null(df) || !nrow(df) || !col %in% names(df)) return(df)
  df[[col]] <- normalize_agg_method_id(df[[col]])
  df
}

# Mirrors paper_select_axis_subset() in the paper-prep workbook: evenly spaced
# subset for c/sigma axes (restart/refine take the prefix). For the full 8x8
# resolution on an 8-point grid this returns all values (no subsetting).
paper_select_axis_subset <- function(vals, n_keep, axis_col) {
  vals <- sort(unique(vals))
  vals <- vals[!is.na(vals)]
  if (!length(vals)) return(vals)
  n_keep <- min(as.integer(n_keep), length(vals))
  if (n_keep >= length(vals)) return(vals)
  if (axis_col %in% c("restart_id", "refine_step")) {
    return(vals[seq_len(n_keep)])
  }
  vals[unique(round(seq(1, length(vals), length.out = n_keep)))]
}

# Paired dataset-bundle bootstrap on pooled confusion bins.
#   For each replicate: resample units (dataset bundles) with replacement,
#   weight each unit's bucket counts by its multiplicity, recompute pooled
#   AUPRC for baseline and treatment, take the difference.
bootstrap_delta <- function(baseline_bins, ensemble_bins, unit_col, B, seed) {
  unit_col_chr <- as.character(unit_col)[[1]]
  shared_units <- intersect(unique(baseline_bins[[unit_col_chr]]), unique(ensemble_bins[[unit_col_chr]]))
  baseline_bins <- baseline_bins |> dplyr::filter(.data[[unit_col_chr]] %in% shared_units)
  ensemble_bins <- ensemble_bins |> dplyr::filter(.data[[unit_col_chr]] %in% shared_units)

  baseline_auprc <- auprc_from_bins(baseline_bins)
  treatment_auprc <- auprc_from_bins(ensemble_bins)
  observed <- treatment_auprc - baseline_auprc

  boot <- rep(NA_real_, B)
  if (length(shared_units) > 1L) {
    set.seed(seed)
    for (i in seq_len(B)) {
      sampled <- sample(shared_units, size = length(shared_units), replace = TRUE)
      weights <- tibble::tibble(unit_value = sampled)
      names(weights)[[1]] <- unit_col_chr
      weights <- weights |>
        dplyr::count(.data[[unit_col_chr]], name = "boot_weight")
      boot[[i]] <- weighted_auprc(ensemble_bins, weights, unit_col_chr) -
        weighted_auprc(baseline_bins, weights, unit_col_chr)
    }
  }
  ci <- stats::quantile(boot, c(0.025, 0.975), na.rm = TRUE)
  # Sign / two-sided tail probabilities (kept SEPARATE from the CI per spec).
  pr_delta_le_0 <- mean(boot <= 0, na.rm = TRUE)
  p_two_sided <- 2 * min(pr_delta_le_0, mean(boot >= 0, na.rm = TRUE))
  list(
    baseline_auprc = baseline_auprc,
    treatment_auprc = treatment_auprc,
    observed_delta = observed,
    ci_low = unname(ci[[1]]),
    ci_high = unname(ci[[2]]),
    pr_delta_le_0 = pr_delta_le_0,
    p_two_sided_bootstrap_tail = min(1, p_two_sided),
    boot = boot,
    shared_units = length(shared_units)
  )
}

near <- function(x, y, tol = 1e-8) is.finite(x) & abs(x - y) <= tol

ensure_column <- function(df, name, value = NA) {
  if (!name %in% names(df)) df[[name]] <- value
  df
}

# ---------------------------------------------------------------------------
# Setup / file discovery
# ---------------------------------------------------------------------------

args <- parse_args()
repo_root <- repo_root_from_script()
if (is.null(args$output_root)) args$output_root <- file.path(repo_root, "output")

job_dir <- file.path(args$output_root, "slurm_output", args$job_name, args$parent_job_id)
consolidated_dir <- file.path(job_dir, "consolidated")
plot_data_rds <- file.path(
  job_dir, "figures", "paper_ensemble_scaling", "plot_data", "paper_ensemble_plot_data.rds"
)
run_history_dir <- file.path(args$output_root, "run_history", args$job_name, args$parent_job_id)

out_dir <- args$out_dir %||% file.path(job_dir, "figures", "paper_ensemble_scaling", "m1_paper_audit")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# Guard: stay on the manuscript primary annotation slice unless told otherwise.
if (!isTRUE(args$allow_nonfocus) && !near(args$annotation_r2_focus, 0.3, tol = 1e-8)) {
  stop(
    "annotation_r2_focus = ", args$annotation_r2_focus,
    " but the manuscript headline slice is phi_a = 0.3. ",
    "Pass --allow_nonfocus=TRUE to override deliberately.",
    call. = FALSE
  )
}

required <- c(
  "auprc_pooled_aggregated.csv",
  "auprc_pooled_agg_by_r2.csv",
  "auprc_pooled_individual_overall.csv",
  "auprc_pooled_individual.csv",
  "confusion_bins_full.csv",
  "model_metrics_full.csv"
)
optional <- c("hg2_by_agg.csv", "scaling_bins_pooled.csv")

manifest <- tibble::tibble(
  file = c(required, optional, "paper_ensemble_plot_data.rds",
           "dataset_bundles.csv", "run_table.csv", "job_config.json"),
  source = c(rep(consolidated_dir, length(required) + length(optional)),
             dirname(plot_data_rds),
             rep(run_history_dir, 3L))
) |>
  dplyr::mutate(
    path = file.path(.data$source, .data$file),
    exists = file.exists(.data$path),
    size_bytes = ifelse(.data$exists, file.info(.data$path)$size, NA_real_),
    mtime = ifelse(.data$exists, as.character(file.info(.data$path)$mtime), NA_character_)
  )
readr::write_csv(manifest, file.path(out_dir, "m1_file_manifest.csv"))

missing_required <- manifest |> dplyr::filter(.data$file %in% required, !.data$exists)
if (nrow(missing_required)) {
  stop(
    "Missing required consolidated files for job ", args$job_name, "/", args$parent_job_id,
    ". This script must run where the consolidated CSVs live (the ROAR Collab ",
    "output tree). See: ", file.path(out_dir, "m1_file_manifest.csv"),
    call. = FALSE
  )
}

read_consolidated <- function(filename) {
  readr::read_csv(file.path(consolidated_dir, filename), show_col_types = FALSE)
}

auprc_pooled_agg <- read_consolidated("auprc_pooled_aggregated.csv") |> normalize_agg_method_df()
auprc_pooled_agg_by_r2 <- read_consolidated("auprc_pooled_agg_by_r2.csv") |> normalize_agg_method_df()
auprc_pooled_individual_overall <- read_consolidated("auprc_pooled_individual_overall.csv")
auprc_pooled_individual <- read_consolidated("auprc_pooled_individual.csv")
confusion_bins <- read_consolidated("confusion_bins_full.csv") |> normalize_agg_method_df()
model_metrics <- read_consolidated("model_metrics_full.csv") |> normalize_agg_method_df()

auprc_pooled_agg <- ensure_column(auprc_pooled_agg, "annotation_r2", NA_real_)
auprc_pooled_agg_by_r2 <- ensure_column(auprc_pooled_agg_by_r2, "annotation_r2", NA_real_)
auprc_pooled_individual_overall <- ensure_column(auprc_pooled_individual_overall, "annotation_r2", NA_real_)
auprc_pooled_individual <- ensure_column(auprc_pooled_individual, "annotation_r2", NA_real_)
confusion_bins <- ensure_column(confusion_bins, "annotation_r2", NA_real_)
confusion_bins <- ensure_column(confusion_bins, "agg_method", NA_character_)
confusion_bins <- ensure_column(confusion_bins, "explore_method", NA_character_)

primary_agg <- args$primary_agg_method
legacy_agg <- args$legacy_agg_method

# ---------------------------------------------------------------------------
# Defensive aggregation-method resolution
# ---------------------------------------------------------------------------

ccs_focus_agg <- auprc_pooled_agg_by_r2 |>
  dplyr::filter(
    .data$spec_name == "C-CS",
    near(as.numeric(.data$annotation_r2), args$annotation_r2_focus)
  )
available_aggs <- sort(unique(stats::na.omit(ccs_focus_agg$agg_method)))

# (1) Fail if the requested (current primary) aggregation method is absent.
if (!primary_agg %in% available_aggs) {
  stop(
    "Requested primary aggregation '", primary_agg, "' not found for C-CS at phi_a = ",
    args$annotation_r2_focus, ". Available: ", paste(available_aggs, collapse = ", "), ". ",
    "The legacy '", legacy_agg, "' (JSD-0.15 Method C) is OBSOLETE and will NOT be ",
    "substituted automatically. Re-collect the job with the scheduled cluster-weight ",
    "methods, or pass --primary_agg_method explicitly.",
    call. = FALSE
  )
}

primary_ccs_auprc <- ccs_focus_agg |>
  dplyr::filter(.data$agg_method == primary_agg) |>
  dplyr::pull(.data$AUPRC) |>
  dplyr::first()

# (2) If the legacy method is ALSO present and disagrees, surface it loudly. We
#     have an explicit selection (primary_agg_method), so we proceed -- but we
#     record the legacy number in the audit and refuse to let it be confused
#     with the headline. (Per spec: never silently pick one.)
legacy_ccs_auprc <- NA_real_
agg_conflict <- FALSE
if (legacy_agg %in% available_aggs) {
  legacy_ccs_auprc <- ccs_focus_agg |>
    dplyr::filter(.data$agg_method == legacy_agg) |>
    dplyr::pull(.data$AUPRC) |>
    dplyr::first()
  agg_conflict <- is.finite(legacy_ccs_auprc) &&
    !near(legacy_ccs_auprc, primary_ccs_auprc, tol = 1e-9)
  if (agg_conflict) {
    message(
      "NOTE: legacy '", legacy_agg, "' AUPRC (", round(legacy_ccs_auprc, 4),
      ") differs from primary '", primary_agg, "' AUPRC (", round(primary_ccs_auprc, 4),
      "). Using the explicitly-selected primary method. Legacy value is recorded ",
      "in the audit for transparency and is NOT the manuscript headline."
    )
  }
}

# ---------------------------------------------------------------------------
# Baseline (cold SuSiE-equivalent), from precomputed pooled-individual overall
# ---------------------------------------------------------------------------

baseline_auprc <- auprc_pooled_individual_overall |>
  dplyr::filter(
    .data$spec_name == "baseline-single",
    .data$use_case_id == "susine_vanilla",
    is.na(.data$annotation_r2)
  ) |>
  dplyr::pull(.data$AUPRC) |>
  dplyr::first()

if (!is.finite(baseline_auprc)) {
  stop("Baseline AUPRC (baseline-single / susine_vanilla, NA annotation) not found.",
       call. = FALSE)
}

figure_delta_cluster <- primary_ccs_auprc - baseline_auprc

# ---------------------------------------------------------------------------
# Confusion-bin slices: baseline, C-CS cluster-weight (primary), oracle members
# ---------------------------------------------------------------------------

baseline_bins <- confusion_bins |>
  dplyr::filter(
    .data$spec_name == "baseline-single",
    .data$use_case_id == "susine_vanilla",
    is.na(.data$agg_method)
  )

ccs_cluster_bins <- confusion_bins |>
  dplyr::filter(
    .data$spec_name == "C-CS",
    near(as.numeric(.data$annotation_r2), args$annotation_r2_focus),
    .data$agg_method == primary_agg
  )

if (!nrow(ccs_cluster_bins)) {
  stop("No C-CS '", primary_agg, "' confusion-bin rows at phi_a = ",
       args$annotation_r2_focus, ".", call. = FALSE)
}

# --- Oracle: best individual C-CS member per dataset by AUPRC, then pool.
#     Ports compute_ccs_resolution_oracle() from the paper-prep workbook.
#     The full C-CS grid is 8 c-values x 8 sigma-values, so the "8x8" resolution
#     is the full ensemble and is the manuscript headline oracle (+0.089). "4x4"
#     is reported too for the resolution diagnostic.
compute_ccs_oracle_member_bins <- function(conf_bins, focus, resolution) {
  res_sizes <- as.integer(strsplit(resolution, "x", fixed = TRUE)[[1]])
  ccs_rows <- conf_bins |>
    dplyr::filter(
      .data$spec_name == "C-CS",
      near(as.numeric(.data$annotation_r2), focus),
      is.na(.data$agg_method) | .data$explore_method != "aggregation"
    ) |>
    dplyr::mutate(
      c_value = suppressWarnings(as.numeric(.data$c_value)),
      sigma_0_2_scalar_num = suppressWarnings(as.numeric(.data$sigma_0_2_scalar))
    )
  if (!nrow(ccs_rows)) return(tibble::tibble())

  selected <- ccs_rows |>
    dplyr::group_by(.data$dataset_bundle_id) |>
    dplyr::group_modify(function(.x, .y) {
      c_vals <- sort(unique(stats::na.omit(.x$c_value)))
      s_vals <- sort(unique(stats::na.omit(.x$sigma_0_2_scalar_num)))
      keep_c <- paper_select_axis_subset(c_vals, res_sizes[[1]], "c_value")
      keep_s <- paper_select_axis_subset(s_vals, res_sizes[[2]], "sigma_0_2_scalar_num")
      .x |>
        dplyr::filter(.data$c_value %in% keep_c, .data$sigma_0_2_scalar_num %in% keep_s)
    }) |>
    dplyr::ungroup()

  # Best member (run_id) per dataset by per-member AUPRC.
  best <- selected |>
    dplyr::group_by(.data$dataset_bundle_id, .data$run_id) |>
    dplyr::summarise(AUPRC = auprc_from_bins(dplyr::pick(dplyr::everything())),
                     .groups = "drop") |>
    dplyr::group_by(.data$dataset_bundle_id) |>
    dplyr::slice_max(.data$AUPRC, n = 1L, with_ties = FALSE) |>
    dplyr::ungroup() |>
    dplyr::select("dataset_bundle_id", "run_id")

  # Keep that member's bins (one row-set per dataset), tagged for pooling.
  selected |>
    dplyr::inner_join(best, by = c("dataset_bundle_id", "run_id")) |>
    dplyr::transmute(
      dataset_bundle_id = .data$dataset_bundle_id,
      run_id = .data$run_id,
      pip_threshold = .data$pip_threshold,
      n_causal_at_bucket = .data$n_causal_at_bucket,
      n_noncausal_at_bucket = .data$n_noncausal_at_bucket
    )
}

oracle_bins_8x8 <- compute_ccs_oracle_member_bins(confusion_bins, args$annotation_r2_focus, "8x8")
oracle_bins_4x4 <- compute_ccs_oracle_member_bins(confusion_bins, args$annotation_r2_focus, "4x4")

oracle_auprc_8x8 <- auprc_from_bins(oracle_bins_8x8)
oracle_auprc_4x4 <- auprc_from_bins(oracle_bins_4x4)
figure_delta_oracle <- oracle_auprc_8x8 - baseline_auprc

# ---------------------------------------------------------------------------
# Paired dataset-bundle bootstrap
# ---------------------------------------------------------------------------

cluster_boot <- bootstrap_delta(
  baseline_bins, ccs_cluster_bins,
  unit_col = "dataset_bundle_id",
  B = args$bootstrap_B, seed = args$bootstrap_seed
)

oracle_boot <- bootstrap_delta(
  baseline_bins, oracle_bins_8x8,
  unit_col = "dataset_bundle_id",
  B = args$bootstrap_B, seed = args$bootstrap_seed
)

# (3) Consistency: the bins-recomputed observed delta must agree with the
#     precomputed figure delta (auprc_pooled_agg_by_r2.csv) for the SAME slice.
#     A meaningful gap means a wrong method/slice/baseline -- fail loudly.
cluster_figure_gap <- abs(cluster_boot$observed_delta - figure_delta_cluster)
if (is.finite(cluster_figure_gap) && cluster_figure_gap > args$figure_match_tol) {
  stop(
    "Bootstrap observed C-CS cluster-weight delta (", round(cluster_boot$observed_delta, 4),
    ", recomputed from confusion bins over ", cluster_boot$shared_units, " shared bundles) ",
    "disagrees with the figure delta (", round(figure_delta_cluster, 4),
    ", from auprc_pooled_agg_by_r2.csv) by ", round(cluster_figure_gap, 4),
    " > tol ", args$figure_match_tol, ". This indicates a slice/method/baseline mismatch. ",
    "Inspect m1_audit_table.csv before trusting any number.",
    call. = FALSE
  )
}

# ---------------------------------------------------------------------------
# Audit table (the thing the user reads to decide the manuscript numbers)
# ---------------------------------------------------------------------------

audit_table <- tibble::tibble(
  quantity = c(
    "baseline_auprc (SuSiE-equiv, figure)",
    "ccs_cluster_weight_auprc (figure, agg_by_r2)",
    "ccs_cluster_weight_auprc (bins, shared bundles)",
    "ccs_oracle_8x8_auprc (bins)",
    "ccs_oracle_4x4_auprc (bins)",
    "ccs_legacy_cluster_weight_auprc (OBSOLETE)",
    "delta_cluster_weight (figure)",
    "delta_cluster_weight (bins, observed)",
    "delta_oracle_8x8 (figure)",
    "delta_oracle_8x8 (bins, observed)",
    "delta_pct_cluster_weight (figure)",
    "delta_pct_oracle_8x8 (figure)"
  ),
  value = c(
    baseline_auprc,
    primary_ccs_auprc,
    cluster_boot$treatment_auprc,
    oracle_auprc_8x8,
    oracle_auprc_4x4,
    legacy_ccs_auprc,
    figure_delta_cluster,
    cluster_boot$observed_delta,
    figure_delta_oracle,
    oracle_boot$observed_delta,
    100 * figure_delta_cluster / baseline_auprc,
    100 * figure_delta_oracle / baseline_auprc
  )
)
readr::write_csv(audit_table, file.path(out_dir, "m1_audit_table.csv"))

# Back-compat headline CSV.
headline <- tibble::tibble(
  quantity = c(
    "baseline_auprc", "ccs_cluster_weight_auprc", "ccs_oracle_8x8_auprc",
    "delta_auprc_cluster_weight", "delta_pct_cluster_weight",
    "delta_auprc_oracle_8x8", "delta_pct_oracle_8x8",
    "agg_method", "annotation_r2_focus", "bootstrap_B", "bootstrap_seed", "parent_job_id"
  ),
  value = c(
    baseline_auprc, primary_ccs_auprc, oracle_auprc_8x8,
    figure_delta_cluster, 100 * figure_delta_cluster / baseline_auprc,
    figure_delta_oracle, 100 * figure_delta_oracle / baseline_auprc,
    primary_agg, args$annotation_r2_focus, args$bootstrap_B, args$bootstrap_seed,
    args$parent_job_id
  )
)
readr::write_csv(headline, file.path(out_dir, "m1_headline_numbers.csv"))

# ---------------------------------------------------------------------------
# Machine-readable bootstrap summary + replicates
# ---------------------------------------------------------------------------

boot_summary <- tibble::tibble(
  comparison = c("C-CS cluster-weight vs baseline", "C-CS oracle (8x8) vs baseline"),
  phi_a = args$annotation_r2_focus,
  ensemble = "C-CS",
  agg_method = c(primary_agg, "oracle_8x8"),
  baseline_auprc = c(cluster_boot$baseline_auprc, oracle_boot$baseline_auprc),
  treatment_auprc = c(cluster_boot$treatment_auprc, oracle_boot$treatment_auprc),
  observed_delta = c(cluster_boot$observed_delta, oracle_boot$observed_delta),
  ci_low = c(cluster_boot$ci_low, oracle_boot$ci_low),
  ci_high = c(cluster_boot$ci_high, oracle_boot$ci_high),
  B = args$bootstrap_B,
  seed = args$bootstrap_seed,
  bootstrap_unit = "dataset_bundle_id",
  shared_units = c(cluster_boot$shared_units, oracle_boot$shared_units),
  pr_delta_le_0 = c(cluster_boot$pr_delta_le_0, oracle_boot$pr_delta_le_0),
  p_two_sided_bootstrap_tail = c(cluster_boot$p_two_sided_bootstrap_tail,
                                 oracle_boot$p_two_sided_bootstrap_tail),
  figure_delta = c(figure_delta_cluster, figure_delta_oracle),
  figure_match_ok = c(
    !(is.finite(cluster_figure_gap) && cluster_figure_gap > args$figure_match_tol),
    TRUE
  )
)
readr::write_csv(boot_summary, file.path(out_dir, "m1_bootstrap_summary.csv"))

boot_replicates <- dplyr::bind_rows(
  tibble::tibble(comparison = "C-CS cluster-weight vs baseline",
                 rep = seq_along(cluster_boot$boot), delta_auprc = cluster_boot$boot),
  tibble::tibble(comparison = "C-CS oracle (8x8) vs baseline",
                 rep = seq_along(oracle_boot$boot), delta_auprc = oracle_boot$boot)
)
readr::write_csv(boot_replicates, file.path(out_dir, "m1_bootstrap_replicates.csv"))

# ---------------------------------------------------------------------------
# RDS cross-check (paper-prep workbook caches its own bootstrap refs)
# ---------------------------------------------------------------------------

rds_observed_delta <- NA_real_
rds_ci_low <- NA_real_
rds_ci_high <- NA_real_
if (file.exists(plot_data_rds)) {
  paper_plot_data <- readRDS(plot_data_rds)
  refs <- paper_plot_data$refs %||% list()
  rds_observed_delta <- refs$observed_delta %||% NA_real_
  rds_ci <- refs$delta_ci %||% c(NA_real_, NA_real_)
  rds_ci_low <- rds_ci[[1]]
  rds_ci_high <- rds_ci[[2]]
  readr::write_csv(
    tibble::tibble(
      quantity = c("rds_observed_delta", "rds_delta_ci_low", "rds_delta_ci_high",
                   "rds_baseline_auprc", "rds_shared_boot_n"),
      value = c(rds_observed_delta, rds_ci_low, rds_ci_high,
                refs$baseline_auprc %||% NA_real_, refs$shared_boot_n %||% NA_real_)
    ),
    file.path(out_dir, "m1_paper_plot_rds_refs.csv")
  )
}

# ---------------------------------------------------------------------------
# Markdown audit note
# ---------------------------------------------------------------------------

fmt <- function(x, d = 4) if (is.finite(x)) formatC(x, digits = d, format = "f") else "NA"
cluster_confirmed <- is.finite(figure_delta_cluster) && near(round(figure_delta_cluster, 3), 0.066, tol = 5e-4)
oracle_confirmed <- is.finite(figure_delta_oracle) && near(round(figure_delta_oracle, 3), 0.089, tol = 5e-4)

note <- c(
  "# M1 Ensemble-Scaling Headline Audit",
  "",
  sprintf("Job: `%s/%s`  (paper figure job; the old default `51250228` was stale)", args$job_name, args$parent_job_id),
  sprintf("Consolidated dir: `%s`", consolidated_dir),
  "",
  "## Slice (matches prepare_results_workbook_ensemble_scaling_paper.Rmd)",
  "",
  sprintf("- ensemble: **C-CS**, phi_a = **%.1f**", args$annotation_r2_focus),
  sprintf("- aggregation: **%s** (current primary, Method B / credible-shift cut 0.05)", primary_agg),
  "- baseline: **baseline-single / susine_vanilla**, NA agg (cold SuSiE-equivalent)",
  "- scoring: pooled top-8 causal AUPRC (auprc_from_pooled_bins)",
  "- oracle: best individual C-CS member per dataset by AUPRC, then pooled (8x8 = full grid)",
  "",
  "## Input files used",
  paste0("- `", file.path(consolidated_dir, required), "`"),
  if (file.exists(plot_data_rds)) sprintf("- `%s` (RDS cross-check)", plot_data_rds) else "- RDS plot-data not present (no cross-check)",
  "",
  "## Observed headline values (figure extraction)",
  "",
  sprintf("| quantity | AUPRC | delta vs baseline | %% |"),
  "|---|---:|---:|---:|",
  sprintf("| baseline (SuSiE-equiv) | %s | - | - |", fmt(baseline_auprc)),
  sprintf("| C-CS %s | %s | %s | %s%% |", primary_agg, fmt(primary_ccs_auprc),
          fmt(figure_delta_cluster), fmt(100 * figure_delta_cluster / baseline_auprc, 1)),
  sprintf("| C-CS oracle (8x8) | %s | %s | %s%% |", fmt(oracle_auprc_8x8),
          fmt(figure_delta_oracle), fmt(100 * figure_delta_oracle / baseline_auprc, 1)),
  sprintf("| C-CS oracle (4x4) | %s | %s | - |", fmt(oracle_auprc_4x4),
          fmt(oracle_auprc_4x4 - baseline_auprc)),
  if (is.finite(legacy_ccs_auprc))
    sprintf("| C-CS %s (OBSOLETE legacy) | %s | %s | - |", legacy_agg, fmt(legacy_ccs_auprc),
            fmt(legacy_ccs_auprc - baseline_auprc))
  else "| legacy cluster_weight | (absent in this job) | - | - |",
  "",
  "## Paired dataset-bundle bootstrap (B = %d, seed = %d, unit = dataset_bundle_id)",
  "",
  sprintf("- C-CS cluster-weight: observed delta %s, 95%% CI [%s, %s], shared bundles %d, Pr(delta<=0) = %s",
          fmt(cluster_boot$observed_delta), fmt(cluster_boot$ci_low), fmt(cluster_boot$ci_high),
          cluster_boot$shared_units, fmt(cluster_boot$pr_delta_le_0, 4)),
  sprintf("- C-CS oracle (8x8): observed delta %s, 95%% CI [%s, %s], shared bundles %d, Pr(delta<=0) = %s",
          fmt(oracle_boot$observed_delta), fmt(oracle_boot$ci_low), fmt(oracle_boot$ci_high),
          oracle_boot$shared_units, fmt(oracle_boot$pr_delta_le_0, 4)),
  sprintf("- Figure-vs-bins consistency (cluster-weight): gap = %s (tol %s) -> %s",
          fmt(cluster_figure_gap), args$figure_match_tol,
          if (is.finite(cluster_figure_gap) && cluster_figure_gap <= args$figure_match_tol) "PASS" else "FAIL"),
  if (file.exists(plot_data_rds))
    sprintf("- RDS cross-check: workbook observed delta %s, CI [%s, %s]",
            fmt(rds_observed_delta), fmt(rds_ci_low), fmt(rds_ci_high))
  else "- RDS cross-check: not available",
  "",
  "## Manuscript verdict",
  "",
  sprintf("- Current text says Delta AUPRC = +0.066 (cluster-weight). This run gives **%s**. -> **%s**",
          fmt(figure_delta_cluster, 3), if (cluster_confirmed) "CONFIRMED" else "UPDATE the manuscript number"),
  sprintf("- Current text says oracle Delta AUPRC = +0.089. This run gives **%s**. -> **%s**",
          fmt(figure_delta_oracle, 3), if (oracle_confirmed) "CONFIRMED" else "UPDATE the manuscript number"),
  "",
  "## Obsolete paths explicitly rejected",
  "",
  "- legacy `cluster_weight` (JSD-0.15 + 1/frequency, Method C): NOT the paper aggregator.",
  "- `annotation_r2 = 0.2` slice with observed_delta ~= 0.0585: stale; the paper slice is phi_a = 0.3.",
  "- default job `51250228`: stale exploratory job; paper figures use `53547760`.",
  if (agg_conflict)
    sprintf("- NOTE: legacy `%s` is present in this job (AUPRC %s) and differs from the primary; it was recorded but NOT used.", legacy_agg, fmt(legacy_ccs_auprc))
  else "- legacy `cluster_weight` not present in this job's C-CS focus slice."
)
# inject B/seed into the bootstrap header line
note <- gsub("B = %d, seed = %d",
             sprintf("B = %d, seed = %d", args$bootstrap_B, args$bootstrap_seed),
             note, fixed = TRUE)
writeLines(note, file.path(out_dir, "m1_audit_note.md"))

cat("Wrote M1 paper audit bundle to:\n", out_dir, "\n", sep = "")
cat("Cluster-weight delta (figure):", fmt(figure_delta_cluster, 4),
    " 95% CI [", fmt(cluster_boot$ci_low), ", ", fmt(cluster_boot$ci_high), "]\n", sep = "")
cat("Oracle 8x8 delta (figure):    ", fmt(figure_delta_oracle, 4),
    " 95% CI [", fmt(oracle_boot$ci_low), ", ", fmt(oracle_boot$ci_high), "]\n", sep = "")
