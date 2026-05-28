#!/usr/bin/env Rscript

# Export real-data barrier diagnostics for the paper.
#
# Run from the test_susine repo root on HPC after the real-data ensemble has
# been collected. This script does not refit any models. It reads the collected
# aggregated outputs, identifies the nontrivially weighted ensemble runs, and
# exports:
#   1. per-run weight/ELBO/prior/cluster details;
#   2. pairwise PIP JSD and yhat-drift among those runs;
#   3. per-locus summaries useful for optimizer- and model-specification-barrier
#      prose.
#
# Example:
#   Rscript "vignettes/real data pipeline/export_real_data_barrier_details.R"
#
# Optional environment overrides:
#   REAL_DATA_JOB_NAME
#   REAL_DATA_PARENT_JOB_ID
#   REAL_DATA_OUTPUT_ROOT
#   REAL_DATA_EXPORT_DIR
#   REAL_DATA_LOCI              comma-separated locus_ids or gene names
#   REAL_DATA_WEIGHT_THRESHOLD  default 0.01
#   REAL_DATA_TOP_N             default 8

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(tidyr)
  library(purrr)
  library(tibble)
})

require_pkg <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop("Package '", pkg, "' is required. Install/load it on HPC and rerun.")
  }
}

require_pkg("arrow")
require_pkg("jsonlite")
require_pkg("test_susine")

repo_root <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
job_name <- Sys.getenv("REAL_DATA_JOB_NAME", "real_data_ensemble_geometric_n20")
parent_job_id <- Sys.getenv("REAL_DATA_PARENT_JOB_ID", "52906940")
output_root <- Sys.getenv("REAL_DATA_OUTPUT_ROOT", file.path(repo_root, "output"))

aggregated_dir <- file.path(output_root, "slurm_output", job_name, parent_job_id, "aggregated")
if (!dir.exists(aggregated_dir)) {
  stop("Aggregated output directory not found: ", aggregated_dir,
       "\nRun collect_results_workbook_real_data_ensemble.Rmd first, or set REAL_DATA_OUTPUT_ROOT.")
}

export_dir <- Sys.getenv(
  "REAL_DATA_EXPORT_DIR",
  file.path(aggregated_dir, "barrier_details_for_paper")
)
dir.create(export_dir, recursive = TRUE, showWarnings = FALSE)

weight_threshold <- as.numeric(Sys.getenv("REAL_DATA_WEIGHT_THRESHOLD", "0.01"))
top_n <- as.integer(Sys.getenv("REAL_DATA_TOP_N", "8"))

default_loci <- c(
  "ydjc_chr22_lung", "arsa_chr22_lung", "rrp7a_chr22_lung",
  "rnf14_chr5_lung", "serf1a_chr5_lung", "tmtc1_chr12_lung",
  "lgals9_chr17_lung", "znf280b_chr22_lung"
)
loci_arg <- Sys.getenv("REAL_DATA_LOCI", "")
requested_tokens <- if (nzchar(loci_arg)) {
  trimws(strsplit(loci_arg, ",", fixed = TRUE)[[1]])
} else {
  default_loci
}

read_csv_required <- function(name) {
  path <- file.path(aggregated_dir, name)
  if (!file.exists(path)) stop("Missing required aggregated file: ", path)
  readr::read_csv(path, show_col_types = FALSE)
}

job_config_path <- file.path(output_root, "run_history", job_name, parent_job_id, "job_config.json")
if (!file.exists(job_config_path)) {
  temp_job_config <- file.path(output_root, "temp", job_name, "job_config.json")
  if (file.exists(temp_job_config)) job_config_path <- temp_job_config
}
if (!file.exists(job_config_path)) {
  stop("Could not find job_config.json under run_history or temp for job ", job_name)
}
job_config <- jsonlite::read_json(job_config_path, simplifyVector = TRUE)

run_manifest <- read_csv_required("run_manifest.csv")
run_metrics <- read_csv_required("run_metrics_full.csv")
agg_weights <- read_csv_required("aggregation_weights_cluster_weight.csv")
cluster_membership <- read_csv_required("cluster_membership.csv")
cluster_summary <- read_csv_required("cluster_summary.csv")
paper_summary <- read_csv_required("paper_real_data_ensemble_summary.csv")

run_comparisons_path <- file.path(aggregated_dir, "run_comparisons.csv")
run_comparisons <- if (file.exists(run_comparisons_path)) {
  readr::read_csv(run_comparisons_path, show_col_types = FALSE)
} else {
  tibble()
}

pairwise_pip_path <- file.path(aggregated_dir, "pairwise_pip_jsd.parquet")
pairwise_pip <- if (file.exists(pairwise_pip_path)) {
  arrow::read_parquet(pairwise_pip_path) %>% tibble::as_tibble()
} else {
  warning("Missing pairwise_pip_jsd.parquet; pairwise PIP JSD will be NA.")
  tibble(locus_id = character(), run_id_i = integer(), run_id_j = integer(),
         jsd = numeric(), same_cluster = logical())
}

fit_file_index_path <- file.path(aggregated_dir, "fit_file_index.csv")
fit_file_index <- if (file.exists(fit_file_index_path)) {
  readr::read_csv(fit_file_index_path, show_col_types = FALSE)
} else {
  tibble()
}

resolve_fit_path <- function(path) {
  if (is.null(path) || length(path) == 0L || is.na(path) || !nzchar(path)) return(NA_character_)
  if (file.exists(path)) return(normalizePath(path, winslash = "/", mustWork = TRUE))
  path2 <- file.path(repo_root, path)
  if (file.exists(path2)) return(normalizePath(path2, winslash = "/", mustWork = TRUE))
  path3 <- file.path(output_root, path)
  if (file.exists(path3)) return(normalizePath(path3, winslash = "/", mustWork = TRUE))
  path
}

metric_fit_cols <- intersect(c("run_id", "fit_rds_path"), names(run_metrics))
fit_paths <- if (length(metric_fit_cols) == 2L) {
  run_metrics %>% select(all_of(c("run_id", "fit_rds_path")))
} else if (all(c("run_id", "path") %in% names(fit_file_index))) {
  fit_file_index %>% transmute(run_id, fit_rds_path = path)
} else {
  tibble(run_id = integer(), fit_rds_path = character())
}

fit_paths <- fit_paths %>%
  mutate(fit_rds_path = vapply(.data$fit_rds_path, resolve_fit_path, character(1)))

locus_lookup <- paper_summary %>% distinct(locus_id, gene_name)
selected_loci <- locus_lookup %>%
  filter(tolower(.data$locus_id) %in% tolower(requested_tokens) |
           tolower(.data$gene_name) %in% tolower(requested_tokens)) %>%
  pull(.data$locus_id) %>%
  unique()

if (!length(selected_loci)) {
  stop("None of REAL_DATA_LOCI/default loci were found in paper_real_data_ensemble_summary.csv")
}

message("Selected loci: ", paste(selected_loci, collapse = ", "))

run_annotation <- agg_weights %>%
  filter(.data$locus_id %in% selected_loci) %>%
  left_join(
    run_metrics %>%
      select(any_of(c(
        "run_id", "elbo_final", "sigma_2_final", "pve_postmean_std",
        "h2_proxy_std", "max_pip", "pip_mass_top1", "pip_mass_top10"
      ))),
    by = "run_id"
  ) %>%
  left_join(fit_paths, by = "run_id") %>%
  left_join(
    cluster_membership %>%
      select(any_of(c(
        "locus_id", "run_id", "cluster_size", "cluster_freq",
        "is_representative", "representative_run_id",
        "representative_elbo", "jsd_to_representative"
      ))),
    by = c("locus_id", "run_id")
  ) %>%
  left_join(
    cluster_summary %>%
      select(all_of(c("locus_id", "cluster_id", "cluster_weight"))),
    by = c("locus_id", "cluster_id")
  )

top_weighted_runs <- run_annotation %>%
  group_by(.data$locus_id) %>%
  arrange(desc(.data$agg_weight_run), desc(.data$elbo_final), .data$run_id, .by_group = TRUE) %>%
  mutate(
    locus_max_elbo = max(.data$elbo_final, na.rm = TRUE),
    delta_elbo_from_locus_max = .data$elbo_final - .data$locus_max_elbo,
    weight_rank = row_number()
  ) %>%
  filter(.data$agg_weight_run >= weight_threshold | .data$weight_rank <= top_n) %>%
  ungroup() %>%
  arrange(.data$locus_id, .data$weight_rank)

comparison_wide <- run_comparisons %>%
  filter(.data$run_id %in% top_weighted_runs$run_id) %>%
  select(any_of(c(
    "run_id", "target_run_id", "comparison_target", "jsd_pip", "l1_pip",
    "top10_overlap", "top20_overlap", "delta_elbo", "delta_pve_postmean_std"
  ))) %>%
  pivot_wider(
    id_cols = "run_id",
    names_from = "comparison_target",
    values_from = any_of(c(
      "jsd_pip", "l1_pip", "top10_overlap", "top20_overlap",
      "delta_elbo", "delta_pve_postmean_std"
    )),
    names_glue = "{.value}_vs_{comparison_target}"
  )

top_weighted_runs <- top_weighted_runs %>%
  left_join(comparison_wide, by = "run_id")

readr::write_csv(
  top_weighted_runs,
  file.path(export_dir, "barrier_top_weighted_runs.csv")
)

safe_internal <- function(name) {
  get(name, envir = asNamespace("test_susine"), inherits = FALSE)
}

load_bundle <- safe_internal("load_real_data_locus_bundle")
cache_fit <- safe_internal("real_data_basin_r2_cache")
yhat_pair <- safe_internal("real_data_basis_drift_var_y_pair_cached")
basin_pair <- safe_internal("real_data_basin_r2_drift_pair_cached")

manifest_path <- job_config$job$manifest_path
job_repo_root <- job_config$paths$repo_root
if (is.null(job_repo_root) || is.na(job_repo_root) || !nzchar(job_repo_root)) {
  job_repo_root <- repo_root
}

fit_cache_env <- new.env(parent = emptyenv())
get_fit_cache <- function(run_id, fit_path, R) {
  key <- as.character(run_id)
  if (exists(key, envir = fit_cache_env, inherits = FALSE)) {
    return(get(key, envir = fit_cache_env, inherits = FALSE))
  }
  if (is.na(fit_path) || !file.exists(fit_path)) {
    warning("Missing fit RDS for run_id=", run_id, ": ", fit_path)
    assign(key, NULL, envir = fit_cache_env)
    return(NULL)
  }
  fit <- tryCatch(readRDS(fit_path), error = function(e) {
    warning("Failed to read fit RDS for run_id=", run_id, ": ", conditionMessage(e))
    NULL
  })
  fc <- if (is.null(fit)) NULL else cache_fit(fit, R)
  assign(key, fc, envir = fit_cache_env)
  fc
}

pairwise_lookup <- pairwise_pip %>%
  mutate(
    run_lo = pmin(.data$run_id_i, .data$run_id_j),
    run_hi = pmax(.data$run_id_i, .data$run_id_j)
  ) %>%
  transmute(
    locus_id = .data$locus_id,
    run_lo = .data$run_lo,
    run_hi = .data$run_hi,
    pairwise_pip_jsd = .data$jsd,
    same_cluster = .data$same_cluster
  )

pairwise_rows <- map_dfr(selected_loci, function(loc) {
  loc_runs <- top_weighted_runs %>%
    filter(.data$locus_id == !!loc) %>%
    arrange(desc(.data$agg_weight_run), .data$run_id)

  if (nrow(loc_runs) < 2L) return(tibble())

  bundle <- tryCatch(
    load_bundle(
      locus_id = loc,
      manifest_path = manifest_path,
      repo_root = job_repo_root
    ),
    error = function(e) {
      warning("Failed to load locus bundle for ", loc, ": ", conditionMessage(e))
      NULL
    }
  )

  pairs <- utils::combn(seq_len(nrow(loc_runs)), 2L, simplify = FALSE)
  map_dfr(pairs, function(idx) {
    a <- loc_runs[idx[[1]], , drop = FALSE]
    b <- loc_runs[idx[[2]], , drop = FALSE]
    run_lo <- min(a$run_id[[1]], b$run_id[[1]])
    run_hi <- max(a$run_id[[1]], b$run_id[[1]])
    pip_row <- pairwise_lookup %>%
      filter(.data$locus_id == !!loc, .data$run_lo == !!run_lo, .data$run_hi == !!run_hi) %>%
      slice_head(n = 1L)

    cache_a <- cache_b <- NULL
    if (!is.null(bundle)) {
      cache_a <- get_fit_cache(a$run_id[[1]], a$fit_rds_path[[1]], bundle$R)
      cache_b <- get_fit_cache(b$run_id[[1]], b$fit_rds_path[[1]], bundle$R)
    }
    basin <- if (!is.null(cache_a) && !is.null(cache_b)) basin_pair(cache_a, cache_b) else list(drift = NA_real_, n_matched = NA_integer_)

    tibble(
      locus_id = loc,
      gene_name = a$gene_name[[1]],
      run_id_a = a$run_id[[1]],
      run_id_b = b$run_id[[1]],
      c_value_a = a$c_value[[1]],
      c_value_b = b$c_value[[1]],
      sigma_0_2_scalar_a = a$sigma_0_2_scalar[[1]],
      sigma_0_2_scalar_b = b$sigma_0_2_scalar[[1]],
      agg_weight_a = a$agg_weight_run[[1]],
      agg_weight_b = b$agg_weight_run[[1]],
      cluster_id_a = a$cluster_id[[1]],
      cluster_id_b = b$cluster_id[[1]],
      elbo_a = a$elbo_final[[1]],
      elbo_b = b$elbo_final[[1]],
      abs_delta_elbo = abs(a$elbo_final[[1]] - b$elbo_final[[1]]),
      pairwise_pip_jsd = if (nrow(pip_row)) pip_row$pairwise_pip_jsd[[1]] else NA_real_,
      same_cluster = if (nrow(pip_row)) pip_row$same_cluster[[1]] else NA,
      pairwise_yhat_drift_var_y = yhat_pair(cache_a, cache_b),
      pairwise_hungarian_basin_r2_drift = basin$drift,
      n_hungarian_matched = basin$n_matched
    )
  })
})

readr::write_csv(
  pairwise_rows,
  file.path(export_dir, "barrier_pairwise_top_runs.csv")
)

locus_barrier_summary <- top_weighted_runs %>%
  group_by(.data$locus_id, .data$gene_name) %>%
  summarise(
    n_exported_runs = n(),
    n_runs_weight_ge_threshold = sum(.data$agg_weight_run >= weight_threshold, na.rm = TRUE),
    weight_ge_threshold = sum(.data$agg_weight_run[.data$agg_weight_run >= weight_threshold], na.rm = TRUE),
    n_clusters_weight_ge_threshold = n_distinct(.data$cluster_id[.data$agg_weight_run >= weight_threshold]),
    max_run_weight = max(.data$agg_weight_run, na.rm = TRUE),
    top3_weight = sum(head(.data$agg_weight_run[order(-.data$agg_weight_run)], 3L), na.rm = TRUE),
    c0_weight_exported = sum(.data$agg_weight_run[.data$c_value == 0], na.rm = TRUE),
    off_c0_weight_exported = sum(.data$agg_weight_run[.data$c_value > 0], na.rm = TRUE),
    elbo_range_exported = max(.data$elbo_final, na.rm = TRUE) - min(.data$elbo_final, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  left_join(
    pairwise_rows %>%
      group_by(.data$locus_id) %>%
      summarise(
        max_pairwise_pip_jsd_exported = max(.data$pairwise_pip_jsd, na.rm = TRUE),
        max_pairwise_yhat_drift_var_y_exported = max(.data$pairwise_yhat_drift_var_y, na.rm = TRUE),
        max_pairwise_abs_delta_elbo_exported = max(.data$abs_delta_elbo, na.rm = TRUE),
        n_pairwise_between_clusters = sum(!.data$same_cluster, na.rm = TRUE),
        .groups = "drop"
      ),
    by = "locus_id"
  ) %>%
  left_join(
    paper_summary %>%
      select(any_of(c(
        "locus_id", "highest_weight_source_agg_weight_run",
        "ensemble_agg_weight_c0_total", "ensemble_agg_weight_off_c0_total",
        "highest_weight_source_c_value", "highest_weight_source_sigma_0_2_scalar",
        "delta_elbo_refit_vs_susie_anchor", "delta_elbo_source_vs_baseline_c0",
        "pve_susie_anchor", "pve_susine_ensemble", "pve_highest_weight_refit",
        "jsd_susie_anchor_vs_susine_ensemble", "jsd_susie_anchor_vs_refit",
        "jsd_susine_ensemble_vs_refit"
      ))),
    by = "locus_id"
  ) %>%
  arrange(desc(.data$ensemble_agg_weight_off_c0_total), .data$locus_id)

readr::write_csv(
  locus_barrier_summary,
  file.path(export_dir, "barrier_locus_summary.csv")
)

readr::write_csv(
  cluster_summary %>% filter(.data$locus_id %in% selected_loci),
  file.path(export_dir, "barrier_cluster_summary_selected_loci.csv")
)

message("Wrote barrier detail exports to: ", export_dir)
message("  - barrier_top_weighted_runs.csv")
message("  - barrier_pairwise_top_runs.csv")
message("  - barrier_locus_summary.csv")
message("  - barrier_cluster_summary_selected_loci.csv")
