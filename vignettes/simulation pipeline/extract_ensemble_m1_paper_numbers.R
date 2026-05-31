#!/usr/bin/env Rscript

# Extract the small set of ensemble-scaling numbers needed for the manuscript
# M1 writeup. Run on the HPC from anywhere after the ensemble collection and
# paper-prep workbooks have produced their outputs.

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
})

`%||%` <- function(x, y) {
  if (is.null(x) || length(x) == 0L) y else x
}

parse_args <- function() {
  args <- commandArgs(trailingOnly = TRUE)
  out <- list(
    job_name = "ensemble_scaling_full",
    parent_job_id = "51250228",
    output_root = NULL,
    bootstrap_B = 2000L,
    bootstrap_seed = 20260505L,
    annotation_r2_focus = 0.3,
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
  out$bootstrap_B <- as.integer(out$bootstrap_B)
  out$bootstrap_seed <- as.integer(out$bootstrap_seed)
  out$annotation_r2_focus <- as.numeric(out$annotation_r2_focus)
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

compute_auprc_single <- function(precision, recall) {
  df <- data.frame(precision = precision, recall = recall) |>
    dplyr::filter(!is.na(.data$precision), !is.na(.data$recall)) |>
    dplyr::arrange(.data$recall)
  if (nrow(df) < 2L) return(NA_real_)
  delta_recall <- c(df$recall[[1]], diff(df$recall))
  sum(df$precision * delta_recall)
}

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
  if (total_p == 0L) return(NA_real_)
  denom <- cum_tp + cum_fp
  precision <- ifelse(denom > 0, cum_tp / denom, NA_real_)
  recall <- cum_tp / total_p
  compute_auprc_single(precision, recall)
}

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

bootstrap_delta <- function(baseline_bins, ensemble_bins, unit_col, B, seed) {
  unit_col_chr <- as.character(unit_col)[[1]]
  shared_units <- intersect(unique(baseline_bins[[unit_col_chr]]), unique(ensemble_bins[[unit_col_chr]]))
  baseline_bins <- baseline_bins |> dplyr::filter(.data[[unit_col_chr]] %in% shared_units)
  ensemble_bins <- ensemble_bins |> dplyr::filter(.data[[unit_col_chr]] %in% shared_units)

  observed <- auprc_from_bins(ensemble_bins) - auprc_from_bins(baseline_bins)

  set.seed(seed)
  boot <- vector("numeric", B)
  for (i in seq_len(B)) {
    sampled <- sample(shared_units, size = length(shared_units), replace = TRUE)
    weights <- tibble::tibble(unit_value = sampled)
    names(weights)[[1]] <- unit_col_chr
    weights <- weights |>
      dplyr::count(.data[[unit_col_chr]], name = "boot_weight")
    boot[[i]] <- weighted_auprc(ensemble_bins, weights, unit_col_chr) -
      weighted_auprc(baseline_bins, weights, unit_col_chr)
  }
  ci <- stats::quantile(boot, c(0.025, 0.975), na.rm = TRUE)
  p_two_sided <- 2 * min(mean(boot <= 0, na.rm = TRUE), mean(boot >= 0, na.rm = TRUE))
  list(
    observed_delta = observed,
    ci_low = unname(ci[[1]]),
    ci_high = unname(ci[[2]]),
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
  stop("Missing required consolidated files. See: ", file.path(out_dir, "m1_file_manifest.csv"), call. = FALSE)
}

read_consolidated <- function(filename) {
  readr::read_csv(file.path(consolidated_dir, filename), show_col_types = FALSE)
}

auprc_pooled_agg <- read_consolidated("auprc_pooled_aggregated.csv")
auprc_pooled_agg_by_r2 <- read_consolidated("auprc_pooled_agg_by_r2.csv")
auprc_pooled_individual_overall <- read_consolidated("auprc_pooled_individual_overall.csv")
auprc_pooled_individual <- read_consolidated("auprc_pooled_individual.csv")
confusion_bins <- read_consolidated("confusion_bins_full.csv")
model_metrics <- read_consolidated("model_metrics_full.csv")

auprc_pooled_agg <- ensure_column(auprc_pooled_agg, "annotation_r2", NA_real_)
auprc_pooled_agg_by_r2 <- ensure_column(auprc_pooled_agg_by_r2, "annotation_r2", NA_real_)
auprc_pooled_individual_overall <- ensure_column(auprc_pooled_individual_overall, "annotation_r2", NA_real_)
auprc_pooled_individual <- ensure_column(auprc_pooled_individual, "annotation_r2", NA_real_)
confusion_bins <- ensure_column(confusion_bins, "annotation_r2", NA_real_)
confusion_bins <- ensure_column(confusion_bins, "agg_method", NA_character_)

baseline_auprc <- auprc_pooled_individual_overall |>
  dplyr::filter(
    .data$spec_name == "baseline-single",
    .data$use_case_id == "susine_vanilla"
  ) |>
  dplyr::pull(.data$AUPRC) |>
  dplyr::first()

ccs_auprc <- auprc_pooled_agg_by_r2 |>
  dplyr::filter(
    .data$spec_name == "C-CS",
    near(as.numeric(.data$annotation_r2), args$annotation_r2_focus),
    .data$agg_method == "cluster_weight"
  ) |>
  dplyr::pull(.data$AUPRC) |>
  dplyr::first()

headline <- tibble::tibble(
  quantity = c(
    "baseline_auprc",
    "ccs_cluster_weight_auprc",
    "delta_auprc",
    "delta_pct",
    "annotation_r2_focus",
    "bootstrap_B",
    "bootstrap_seed"
  ),
  value = c(
    baseline_auprc,
    ccs_auprc,
    ccs_auprc - baseline_auprc,
    100 * (ccs_auprc - baseline_auprc) / baseline_auprc,
    args$annotation_r2_focus,
    args$bootstrap_B,
    args$bootstrap_seed
  )
)
readr::write_csv(headline, file.path(out_dir, "m1_headline_numbers.csv"))

heatmap_focus <- auprc_pooled_agg_by_r2 |>
  dplyr::filter(
    (.data$spec_name %in% c("A-R", "A-F", "A-S", "A-RF", "A-RS", "A-FS", "A-RFS") & is.na(.data$annotation_r2)) |
      (!.data$spec_name %in% c("A-R", "A-F", "A-S", "A-RF", "A-RS", "A-FS", "A-RFS") &
         near(as.numeric(.data$annotation_r2), args$annotation_r2_focus))
  ) |>
  dplyr::mutate(
    baseline_auprc = baseline_auprc,
    delta_auprc = .data$AUPRC - baseline_auprc,
    delta_pct = 100 * .data$delta_auprc / baseline_auprc
  )
readr::write_csv(heatmap_focus, file.path(out_dir, "m1_heatmap_focus_delta_auprc.csv"))

ccs_by_r2 <- auprc_pooled_agg_by_r2 |>
  dplyr::filter(.data$spec_name == "C-CS") |>
  dplyr::mutate(
    baseline_auprc = baseline_auprc,
    delta_auprc = .data$AUPRC - baseline_auprc,
    delta_pct = 100 * .data$delta_auprc / baseline_auprc
  )
readr::write_csv(ccs_by_r2, file.path(out_dir, "m1_ccs_by_annotation_r2.csv"))

baseline_bins <- confusion_bins |>
  dplyr::filter(
    .data$spec_name == "baseline-single",
    .data$use_case_id == "susine_vanilla",
    is.na(.data$agg_method)
  )
ccs_bins <- confusion_bins |>
  dplyr::filter(
    .data$spec_name == "C-CS",
    near(as.numeric(.data$annotation_r2), args$annotation_r2_focus),
    .data$agg_method == "cluster_weight"
  )

shared_dataset_ids <- intersect(unique(baseline_bins$dataset_bundle_id), unique(ccs_bins$dataset_bundle_id))
baseline_bins_ds <- baseline_bins |> dplyr::filter(.data$dataset_bundle_id %in% shared_dataset_ids)
ccs_bins_ds <- ccs_bins |> dplyr::filter(.data$dataset_bundle_id %in% shared_dataset_ids)

dataset_boot <- bootstrap_delta(
  baseline_bins_ds, ccs_bins_ds,
  unit_col = "dataset_bundle_id",
  B = args$bootstrap_B,
  seed = args$bootstrap_seed
)

boot_summary <- tibble::tibble(
  resampling_unit = "dataset_bundle_id",
  B = args$bootstrap_B,
  seed = args$bootstrap_seed,
  shared_units = dataset_boot$shared_units,
  shared_datasets = length(shared_dataset_ids),
  observed_delta = dataset_boot$observed_delta,
  ci_low = dataset_boot$ci_low,
  ci_high = dataset_boot$ci_high,
  p_two_sided_bootstrap_tail = dataset_boot$p_two_sided_bootstrap_tail
)
readr::write_csv(
  tibble::tibble(rep = seq_along(dataset_boot$boot), delta_auprc = dataset_boot$boot),
  file.path(out_dir, "m1_bootstrap_dataset_level_replicates.csv")
)

dataset_bundles_path <- file.path(run_history_dir, "dataset_bundles.csv")
if (file.exists(dataset_bundles_path)) {
  dataset_bundles <- readr::read_csv(dataset_bundles_path, show_col_types = FALSE) |>
    dplyr::select(.data$dataset_bundle_id, .data$matrix_id) |>
    dplyr::distinct()
  baseline_bins_matrix <- baseline_bins_ds |>
    dplyr::left_join(dataset_bundles, by = "dataset_bundle_id") |>
    dplyr::filter(!is.na(.data$matrix_id))
  ccs_bins_matrix <- ccs_bins_ds |>
    dplyr::left_join(dataset_bundles, by = "dataset_bundle_id") |>
    dplyr::filter(!is.na(.data$matrix_id))

  if (nrow(baseline_bins_matrix) && nrow(ccs_bins_matrix)) {
    matrix_boot <- bootstrap_delta(
      baseline_bins_matrix, ccs_bins_matrix,
      unit_col = "matrix_id",
      B = args$bootstrap_B,
      seed = args$bootstrap_seed
    )
    boot_summary <- dplyr::bind_rows(
      boot_summary,
      tibble::tibble(
        resampling_unit = "matrix_id_cluster",
        B = args$bootstrap_B,
        seed = args$bootstrap_seed,
        shared_units = matrix_boot$shared_units,
        shared_datasets = length(shared_dataset_ids),
        observed_delta = matrix_boot$observed_delta,
        ci_low = matrix_boot$ci_low,
        ci_high = matrix_boot$ci_high,
        p_two_sided_bootstrap_tail = matrix_boot$p_two_sided_bootstrap_tail
      )
    )
    readr::write_csv(
      tibble::tibble(rep = seq_along(matrix_boot$boot), delta_auprc = matrix_boot$boot),
      file.path(out_dir, "m1_bootstrap_matrix_cluster_replicates.csv")
    )
  }
}
readr::write_csv(boot_summary, file.path(out_dir, "m1_bootstrap_summary.csv"))

if (file.exists(plot_data_rds)) {
  paper_plot_data <- readRDS(plot_data_rds)
  refs <- paper_plot_data$refs %||% list()
  rds_refs <- tibble::tibble(
    quantity = c("rds_observed_delta", "rds_delta_ci_low", "rds_delta_ci_high", "rds_shared_boot_n"),
    value = c(
      refs$observed_delta %||% NA_real_,
      (refs$delta_ci %||% c(NA_real_, NA_real_))[[1]],
      (refs$delta_ci %||% c(NA_real_, NA_real_))[[2]],
      refs$shared_boot_n %||% NA_real_
    )
  )
  readr::write_csv(rds_refs, file.path(out_dir, "m1_paper_plot_rds_refs.csv"))
  if (!is.null(paper_plot_data$diagnostics$pr_precision50_points)) {
    readr::write_csv(
      paper_plot_data$diagnostics$pr_precision50_points,
      file.path(out_dir, "m1_pr_precision50_operating_points.csv")
    )
  }
  if (!is.null(paper_plot_data$diagnostics$pr_same_recall_points)) {
    readr::write_csv(
      paper_plot_data$diagnostics$pr_same_recall_points,
      file.path(out_dir, "m1_pr_same_recall_operating_points.csv")
    )
  }
}

methods_lines <- c(
  "# M1 Methods Text Stub",
  "",
  sprintf(
    "For the ensemble-scaling headline comparison, we compared the pooled AUPRC of the C-CS SuSiNE ensemble with cluster-weight aggregation at phi_a = %.1f against the zero-prior-mean SuSiE-equivalent baseline fit.",
    args$annotation_r2_focus
  ),
  sprintf(
    "The observed pooled AUPRC gain was %.4f (%.1f%% relative to the baseline AUPRC of %.4f).",
    ccs_auprc - baseline_auprc,
    100 * (ccs_auprc - baseline_auprc) / baseline_auprc,
    baseline_auprc
  ),
  sprintf(
    "The original paper-prep workflow used a paired bootstrap over dataset_bundle_id with B = %d resamples and RNG seed %d; its percentile 95%% CI was [%.4f, %.4f].",
    args$bootstrap_B,
    args$bootstrap_seed,
    boot_summary$ci_low[boot_summary$resampling_unit == "dataset_bundle_id"],
    boot_summary$ci_high[boot_summary$resampling_unit == "dataset_bundle_id"]
  )
)
if ("matrix_id_cluster" %in% boot_summary$resampling_unit) {
  methods_lines <- c(
    methods_lines,
    sprintf(
      "Because each genotype matrix contributes multiple phenotype seeds, we also computed a clustered paired bootstrap that resamples matrix_id units and carries all phenotype replicates for each sampled matrix; its percentile 95%% CI was [%.4f, %.4f].",
      boot_summary$ci_low[boot_summary$resampling_unit == "matrix_id_cluster"],
      boot_summary$ci_high[boot_summary$resampling_unit == "matrix_id_cluster"]
    )
  )
} else {
  methods_lines <- c(
    methods_lines,
    "A matrix-clustered bootstrap was not computed because dataset_bundles.csv with matrix_id was unavailable in run_history."
  )
}
writeLines(methods_lines, file.path(out_dir, "m1_methods_text_stub.md"))

cat("Wrote M1 paper audit bundle to:\n", out_dir, "\n", sep = "")
