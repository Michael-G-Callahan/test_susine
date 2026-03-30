# Task/flush staging collection helpers -----------------------------------

#' Resolve standard directory paths for a completed job.
#'
#' @param job_name Job identifier.
#' @param parent_job_id Parent SLURM array/job id.
#' @param output_root Root directory for outputs.
#' @return Named list with paths: run_history_dir, results_dir, staging_dir, aggregated_dir, snps_dataset_dir
#' @export
resolve_job_paths <- function(job_name,
                              parent_job_id,
                              output_root = "output") {

  run_history_dir <- file.path(output_root, "run_history", job_name, parent_job_id)
  results_dir <- file.path(output_root, "slurm_output", job_name, parent_job_id)
  staging_dir <- results_dir
  aggregated_dir <- file.path(results_dir, "aggregated")
  snps_dataset_dir <- file.path(aggregated_dir, "snps_dataset")

list(
    run_history_dir = run_history_dir,
    results_dir = results_dir,
    staging_dir = staging_dir,
    aggregated_dir = aggregated_dir,
    snps_dataset_dir = snps_dataset_dir
  )
}

#' Return the staging directory used by task-level flush outputs.
#' @keywords internal
staging_base_dir <- function(job_name,
                             parent_job_id,
                             output_root = "output") {
  file.path(output_root, "slurm_output", job_name, parent_job_id)
}

#' Index all task/flush files written during a job.
#'
#' @param job_name Job identifier.
#' @param parent_job_id Parent SLURM array/job id.
#' @param output_root Root directory for outputs.
#' @return Tibble with columns `task_id`, `flush_id`, `type`, and `path`.
#' @export
index_staging_outputs <- function(job_name,
                                  parent_job_id,
                                  output_root = "output") {
  staging_dir <- staging_base_dir(job_name, parent_job_id, output_root)
  if (!dir.exists(staging_dir)) {
    stop("Staging directory not found: ", staging_dir)
  }
  task_dirs <- list.dirs(staging_dir, recursive = FALSE, full.names = TRUE)
  if (!length(task_dirs)) {
    return(tibble::tibble(task_id = integer(), flush_id = character(), type = character(), path = character()))
  }
  parse_type <- function(fname) {
    if (grepl("_model_metrics\\.csv$", fname)) return("model_metrics")
    if (grepl("_effect_metrics\\.csv$", fname)) return("effect_metrics")
    if (grepl("_validation\\.csv$", fname)) return("validation")
    if (grepl("_confusion_bins\\.csv$", fname)) return("confusion_bins")
    if (grepl("_dataset_metrics\\.csv$", fname)) return("dataset_metrics")
    if (grepl("_multimodal_metrics\\.csv$", fname)) return("multimodal_metrics")
    if (grepl("_prior_diagnostics\\.csv$", fname)) return("prior_diagnostics")
    if (grepl("_tier_cs_metrics\\.csv$", fname)) return("tier_cs_metrics")
    if (grepl("_scaling_bins\\.csv$", fname)) return("scaling_bins")
    if (grepl("_hg2_by_agg\\.csv$", fname)) return("hg2_by_agg")
    if (grepl("_snps\\.parquet$", fname)) return("snps")
    NA_character_
  }
  purrr::map_dfr(task_dirs, function(td) {
    task_lbl <- basename(td)
    task_id <- suppressWarnings(as.integer(sub("^task-", "", task_lbl)))
    files <- list.files(td, pattern = "^flush-[0-9]+_", full.names = TRUE)
    if (!length(files)) return(NULL)
    tibble::tibble(
      task_id = task_id,
      flush_id = sub("^flush-([0-9]+)_.*", "\\1", basename(files)),
      type = vapply(basename(files), parse_type, character(1)),
      path = files
    ) %>%
      dplyr::filter(!is.na(.data$type))
  }) -> out

  if (is.null(out) || !nrow(out) || !all(c("task_id", "flush_id", "type", "path") %in% names(out))) {
    return(tibble::tibble(
      task_id = integer(),
      flush_id = character(),
      type = character(),
      path = character()
    ))
  }

  out %>% dplyr::arrange(.data$task_id, .data$flush_id, .data$type)
}

#' Validate that staging files are readable.
#'
#' @param file_index Output of [index_staging_outputs()].
#' @return Tibble with columns `task_id`, `flush_id`, `type`, `path`, `ok`, `error`.
#' @export
validate_staging_outputs <- function(file_index) {
  if (is.null(file_index) || !nrow(file_index)) {
    return(tibble::tibble(task_id = integer(), flush_id = character(), type = character(), path = character(), ok = logical(), error = character()))
  }
  file_index <- dplyr::select(file_index, "task_id", "flush_id", "type", "path")
  purrr::pmap_dfr(file_index, function(task_id, flush_id, type, path) {
    reader <- switch(
      type,
      model_metrics = function(p) readr::read_csv(p, show_col_types = FALSE),
      effect_metrics = function(p) readr::read_csv(p, show_col_types = FALSE),
      validation = function(p) readr::read_csv(p, show_col_types = FALSE),
      confusion_bins = function(p) readr::read_csv(p, show_col_types = FALSE),
      dataset_metrics = function(p) readr::read_csv(p, show_col_types = FALSE),
      multimodal_metrics = function(p) readr::read_csv(p, show_col_types = FALSE),
      prior_diagnostics = function(p) readr::read_csv(p, show_col_types = FALSE),
      tier_cs_metrics = function(p) readr::read_csv(p, show_col_types = FALSE),
      snps = function(p) arrow::read_parquet(p),
      function(p) stop("Unknown type: ", type)
    )
    res <- tryCatch({
      reader(path)
      list(ok = TRUE, error = NA_character_)
    }, error = function(e) {
      list(ok = FALSE, error = conditionMessage(e))
    })
    tibble::tibble(
      task_id = task_id,
      flush_id = flush_id,
      type = type,
      path = path,
      ok = res$ok,
      error = res$error
    )
  })
}

#' Aggregate task/flush staging files into job-level outputs.
#'
#' @param job_name Job identifier.
#' @param parent_job_id Parent SLURM array/job id.
#' @param output_root Root directory for outputs.
#' @param validate When TRUE, read each staging file first and fail fast on errors.
#' @param output_dir Optional destination; defaults to `combined/aggregated`.
#' @param write_snps_dataset When FALSE, skip writing the aggregated snps_dataset parquet output.
#' @return List with aggregated paths and validation summary.
#' @export
aggregate_staging_outputs <- function(job_name,
                                      parent_job_id,
                                      output_root = "output",
                                      validate = TRUE,
                                      output_dir = NULL,
                                      run_table_path = NULL,
                                      write_snps_dataset = TRUE) {
  idx <- index_staging_outputs(job_name, parent_job_id, output_root)
  validation <- NULL
  if (validate) {
    validation <- validate_staging_outputs(idx)
    if (nrow(validation) && any(!validation$ok, na.rm = TRUE)) {
      bad <- dplyr::filter(validation, !ok)
      stop("Staging validation failed for ", nrow(bad), " file(s); fix before aggregating.")
    }
  }
  results_dir <- file.path(output_root, "slurm_output", job_name, parent_job_id)
  if (is.null(output_dir)) {
  output_dir <- file.path(results_dir, "aggregated")
  }
  ensure_dir(output_dir)

  # Resolve dimension table paths: check temp first, then run_history

  if (is.null(run_table_path)) {
    temp_path <- file.path(output_root, "temp", job_name, "run_table.csv")
    history_path <- file.path(output_root, "run_history", job_name, parent_job_id, "run_table.csv")
    if (file.exists(temp_path)) {
      run_table_path <- temp_path
    } else if (file.exists(history_path)) {
      run_table_path <- history_path
    }
  }

  dataset_bundles_path <- NULL
  temp_db_path <- file.path(output_root, "temp", job_name, "dataset_bundles.csv")
  history_db_path <- file.path(output_root, "run_history", job_name, parent_job_id, "dataset_bundles.csv")
  if (file.exists(temp_db_path)) {
    dataset_bundles_path <- temp_db_path
  } else if (file.exists(history_db_path)) {
    dataset_bundles_path <- history_db_path
  }

  current_mem_gb <- function() {
    rss_gb <- NA_real_
    if (.Platform$OS.type != "windows") {
      rss <- suppressWarnings(tryCatch(
        system("ps -o rss= -p $$", intern = TRUE),
        error = function(e) NA
      ))
      if (length(rss) && !is.na(rss[1])) {
        rss_gb <- as.numeric(rss[1]) / 1024 / 1024
      }
    } else if ("memory.size" %in% ls("package:utils")) {
      rss_gb <- tryCatch(utils::memory.size() / 1024, error = function(e) NA_real_)
    }
    rss_gb
  }
  log_progress <- function(msg) {
    mem <- current_mem_gb()
    if (is.finite(mem)) {
      message(sprintf("%s | mem=%.2f GB", msg, mem))
    } else {
      message(msg)
    }
  }

  # Load dimension tables for join-enrichment of slim flush files.
  # Flush files carry only run_id or dataset_bundle_id as join keys;
  # dimension columns are re-attached here so aggregated CSVs stay
  # analysis-ready.
  run_lookup <- NULL
  if (!is.null(run_table_path) && file.exists(run_table_path)) {
    run_lookup <- readr::read_csv(run_table_path, show_col_types = FALSE) %>%
      dplyr::mutate(run_id = as.integer(.data$run_id)) %>%
      dplyr::distinct(.data$run_id, .keep_all = TRUE)
    log_progress(sprintf("Loaded run_lookup (%d rows) for enrichment", nrow(run_lookup)))
  }

  bundle_lookup <- NULL
  if (!is.null(dataset_bundles_path) && file.exists(dataset_bundles_path)) {
    bundle_lookup <- readr::read_csv(dataset_bundles_path, show_col_types = FALSE) %>%
      dplyr::mutate(dataset_bundle_id = as.integer(.data$dataset_bundle_id)) %>%
      dplyr::distinct(.data$dataset_bundle_id, .keep_all = TRUE)
    log_progress(sprintf("Loaded bundle_lookup (%d rows) for enrichment", nrow(bundle_lookup)))
  }

  # Helper: enrich a data frame by joining dimension columns that were
  # stripped from flush files. Only adds columns specified by `cols`;
  # skips columns already present in df.
  enrich_from_run_table <- function(df, cols = NULL) {
    if (is.null(run_lookup) || !"run_id" %in% names(df)) return(df)
    available <- setdiff(names(run_lookup), c("run_id", names(df)))
    if (!is.null(cols)) available <- intersect(available, cols)
    if (!length(available)) return(df)
    lookup_slim <- dplyr::select(run_lookup, dplyr::all_of(c("run_id", available)))
    dplyr::left_join(df, lookup_slim, by = "run_id")
  }

  enrich_from_bundle_table <- function(df, cols = NULL) {
    if (is.null(bundle_lookup) || !"dataset_bundle_id" %in% names(df)) return(df)
    available <- setdiff(names(bundle_lookup), c("dataset_bundle_id", names(df)))
    if (!is.null(cols)) available <- intersect(available, cols)
    if (!length(available)) return(df)
    lookup_slim <- dplyr::select(bundle_lookup, dplyr::all_of(c("dataset_bundle_id", available)))
    dplyr::left_join(df, lookup_slim, by = "dataset_bundle_id")
  }

  read_csv_safe <- function(paths, idx_tbl) {
    log_progress(sprintf("Reading %d CSV file(s)...", length(paths)))
    n <- length(paths)
    pb <- if (n) utils::txtProgressBar(min = 0, max = n, style = 3, initial = 0, width = NA) else NULL
    on.exit(if (!is.null(pb)) close(pb), add = TRUE)
    out <- vector("list", n)
    for (i in seq_along(paths)) {
      p <- paths[[i]]
      df <- readr::read_csv(p, show_col_types = FALSE)
      meta <- dplyr::filter(idx_tbl, .data$path == !!p)
      if (!"task_id" %in% names(df)) df$task_id <- meta$task_id[[1]]
      if (!"flush_id" %in% names(df)) df$flush_id <- meta$flush_id[[1]]
      out[[i]] <- df
      if (!is.null(pb)) utils::setTxtProgressBar(pb, i)
    }
    dplyr::bind_rows(out)
  }

  model_files <- dplyr::filter(idx, .data$type == "model_metrics")$path
  effect_files <- dplyr::filter(idx, .data$type == "effect_metrics")$path
  validation_files <- dplyr::filter(idx, .data$type == "validation")$path
  confusion_files <- dplyr::filter(idx, .data$type == "confusion_bins")$path
  dataset_files <- dplyr::filter(idx, .data$type == "dataset_metrics")$path
  multimodal_files <- dplyr::filter(idx, .data$type == "multimodal_metrics")$path
  prior_diag_files   <- dplyr::filter(idx, .data$type == "prior_diagnostics")$path
  tier_cs_files      <- dplyr::filter(idx, .data$type == "tier_cs_metrics")$path
  scaling_bin_files <- dplyr::filter(idx, .data$type == "scaling_bins")$path
  hg2_by_agg_files  <- dplyr::filter(idx, .data$type == "hg2_by_agg")$path
  snp_files <- dplyr::filter(idx, .data$type == "snps")$path

  # Column lists for join-enrichment — match the original aggregated schemas.
  run_enrich_core <- c("use_case_id", "phenotype_seed", "dataset_bundle_id",
                       "architecture", "refine_step", "run_type")
  run_enrich_full <- c(run_enrich_core, "group_key", "L", "annotation_r2",
                       "inflate_match", "sigma_0_2_scalar", "c_value", "tau_value",
                       "matrix_id", "y_noise", "p_star", "exploration_mode",
                       "exploration_methods", "exploration_group",
                       "alpha_concentration", "warm_method")
  bundle_enrich_cols <- c("matrix_id", "architecture", "y_noise", "p_star",
                          "phenotype_seed")

  if (length(model_files)) {
    model_tbl <- read_csv_safe(model_files, idx) %>%
      enrich_from_run_table(cols = run_enrich_full)
    readr::write_csv(model_tbl, file.path(output_dir, "model_metrics.csv"))
    log_progress("Wrote model_metrics.csv")
    rm(model_tbl, model_files)
    gc()
  }
  if (length(effect_files)) {
    effect_tbl <- read_csv_safe(effect_files, idx) %>%
      enrich_from_run_table(cols = run_enrich_core)
    readr::write_csv(effect_tbl, file.path(output_dir, "effect_metrics.csv"))
    log_progress("Wrote effect_metrics.csv")
    rm(effect_tbl, effect_files)
    gc()
  }
  if (length(validation_files)) {
    validation_tbl <- read_csv_safe(validation_files, idx)
    readr::write_csv(validation_tbl, file.path(output_dir, "validation.csv"))
    log_progress("Wrote validation.csv")
    rm(validation_tbl, validation_files)
    gc()
  }
  if (length(confusion_files)) {
    conf_tbl <- read_csv_safe(confusion_files, idx) %>%
      enrich_from_run_table(cols = run_enrich_full)
    readr::write_csv(conf_tbl, file.path(output_dir, "confusion_bins.csv"))
    log_progress("Wrote confusion_bins.csv")
    rm(conf_tbl, confusion_files)
    gc()
  }
  if (length(dataset_files)) {
    dataset_tbl <- read_csv_safe(dataset_files, idx) %>%
      enrich_from_bundle_table(cols = bundle_enrich_cols)
    readr::write_csv(dataset_tbl, file.path(output_dir, "dataset_metrics.csv"))
    log_progress("Wrote dataset_metrics.csv")
    rm(dataset_tbl, dataset_files)
    gc()
  }
  if (length(multimodal_files)) {
    multimodal_tbl <- read_csv_safe(multimodal_files, idx) %>%
      enrich_from_bundle_table(cols = bundle_enrich_cols)
    readr::write_csv(multimodal_tbl, file.path(output_dir, "multimodal_metrics.csv"))
    log_progress("Wrote multimodal_metrics.csv")
    rm(multimodal_tbl, multimodal_files)
    gc()
  }
  if (length(prior_diag_files)) {
    prior_diag_tbl <- read_csv_safe(prior_diag_files, idx) %>%
      enrich_from_run_table(cols = run_enrich_core)
    readr::write_csv(prior_diag_tbl, file.path(output_dir, "prior_diagnostics.csv"))
    log_progress("Wrote prior_diagnostics.csv")
    rm(prior_diag_tbl, prior_diag_files)
    gc()
  }
  if (length(tier_cs_files)) {
    tier_cs_tbl <- read_csv_safe(tier_cs_files, idx) %>%
      enrich_from_run_table(cols = run_enrich_core)
    readr::write_csv(tier_cs_tbl, file.path(output_dir, "tier_cs_metrics.csv"))
    log_progress("Wrote tier_cs_metrics.csv")
    rm(tier_cs_tbl, tier_cs_files)
    gc()
  }
  if (length(scaling_bin_files)) {
    scaling_bins_tbl <- read_csv_safe(scaling_bin_files, idx)
    readr::write_csv(scaling_bins_tbl, file.path(output_dir, "scaling_bins.csv"))
    log_progress("Wrote scaling_bins.csv")
    rm(scaling_bins_tbl, scaling_bin_files)
    gc()
  }
  if (length(hg2_by_agg_files)) {
    hg2_by_agg_tbl <- read_csv_safe(hg2_by_agg_files, idx) %>%
      enrich_from_bundle_table(cols = bundle_enrich_cols)
    readr::write_csv(hg2_by_agg_tbl, file.path(output_dir, "hg2_by_agg.csv"))
    log_progress("Wrote hg2_by_agg.csv")
    rm(hg2_by_agg_tbl, hg2_by_agg_files)
    gc()
  }
  if (length(snp_files) && isTRUE(write_snps_dataset)) {
    log_progress(sprintf("Writing SNP dataset from %d parquet fragment(s)...", length(snp_files)))
    snp_dataset <- arrow::open_dataset(snp_files, format = "parquet")

    # Enrich slim parquet with run-table dimension columns.
    # Slim parquet carries only run_id as join key; stripped cols
    # (use_case_id, task_id, dataset_bundle_id, architecture, refine_step,
    # run_type) are re-attached here. use_case_id is required for partitioning.
    snp_enrich_cols <- c("use_case_id", "task_id", "dataset_bundle_id",
                         "architecture", "refine_step", "run_type")
    if (!is.null(run_lookup)) {
      snp_cols <- names(snp_dataset)
      needed_cols <- setdiff(intersect(snp_enrich_cols, names(run_lookup)), snp_cols)
      if (length(needed_cols)) {
        snp_lookup <- dplyr::select(run_lookup, dplyr::all_of(c("run_id", needed_cols))) %>%
          dplyr::distinct()
        snp_lookup_arrow <- arrow::Table$create(snp_lookup)
        log_progress("Joining run-table columns into SNP dataset...")
        snp_dataset <- snp_dataset %>%
          dplyr::left_join(snp_lookup_arrow, by = "run_id")
        rm(snp_lookup, snp_lookup_arrow)
        gc()
      }
    } else if (!"use_case_id" %in% names(snp_dataset)) {
      stop(
        "SNP parquet files lack 'use_case_id' and no run_table found to map it. ",
        "Provide run_table_path or ensure run_table.csv exists in temp/ or run_history/."
      )
    }

    snp_out_dir <- file.path(output_dir, "snps_dataset")
    if (dir.exists(snp_out_dir)) {
      unlink(snp_out_dir, recursive = TRUE)
    }
    arrow::write_dataset(
      snp_dataset,
      path = snp_out_dir,
      format = "parquet",
      partitioning = "use_case_id"
    )
    log_progress("Finished writing snps_dataset (partitioned by use_case_id)")
    rm(snp_dataset)
    gc()
  } else if (length(snp_files) && !isTRUE(write_snps_dataset)) {
    log_progress("Skipping snps_dataset parquet aggregation (write_snps_dataset = FALSE)")
  }

  list(
    output_dir = output_dir,
    validation = validation,
    file_index = idx,
    snp_files = snp_files
  )
}

#' Compute AUPRC rank of each run within its dataset bundle.
#'
#' Ranks all runs sharing the same `dataset_bundle_id` by AUPRC (descending).
#' Uses only purity-filtered rows from the model metrics table.
#'
#' @param model_metrics_tbl Tibble of aggregated model metrics (must contain
#'   `run_id`, `dataset_bundle_id`, `AUPRC`, and `filtering` columns).
#' @return Tibble with columns `run_id`, `auprc_rank`, `n_competitors`.
#' @export
compute_model_rank <- function(model_metrics_tbl) {
  if (is.null(model_metrics_tbl) || !nrow(model_metrics_tbl)) {
    return(tibble::tibble(
      run_id = integer(),
      auprc_rank = integer(),
      n_competitors = integer()
    ))
  }
  model_metrics_tbl %>%
    dplyr::filter(.data$filtering == "purity_filtered") %>%
    dplyr::group_by(.data$dataset_bundle_id) %>%
    dplyr::mutate(
      auprc_rank = rank(-.data$AUPRC, ties.method = "min"),
      n_competitors = dplyr::n()
    ) %>%
    dplyr::ungroup() %>%
    dplyr::select("run_id", "auprc_rank", "n_competitors")
}

# ---- Scaling analysis helpers -----

#' Aggregate a p x N PIP matrix to a single length-p PIP vector.
#'
#' @param pip_mat Numeric matrix (p rows × N columns).
#' @param elbos Length-N ELBO vector.
#' @param method One of "uniform", "max_elbo", "elbo_softmax",
#'   "cluster_weight", "cluster_weight_050".
#' @param jsd_threshold JSD cutoff for cluster_weight (default 0.15).
#' @keywords internal
aggregate_pip_matrix <- function(pip_mat, elbos, method, jsd_threshold = 0.15) {
  N <- ncol(pip_mat)
  if (N == 1L) return(as.vector(pip_mat))
  switch(method,
    uniform = rowMeans(pip_mat),
    max_elbo = as.vector(pip_mat[, which.max(elbos), drop = TRUE]),
    elbo_softmax = {
      w <- exp(elbos - max(elbos)); w <- w / sum(w)
      as.vector(pip_mat %*% w)
    },
    cluster_weight     = .aggregate_cluster_weight(pip_mat, elbos, jsd_threshold),
    cluster_weight_050 = .aggregate_cluster_weight(pip_mat, elbos, 0.50),
    stop("Unknown aggregation method: ", method)
  )
}

.aggregate_cluster_weight <- function(pip_mat, elbos, jsd_threshold) {
  N <- ncol(pip_mat); eps <- 1e-10
  H <- apply(pip_mat + eps, 2L, function(pv) -sum(pv * log(pv)))
  jsd <- matrix(0, N, N)
  for (i in seq_len(N - 1L)) {
    for (j in seq(i + 1L, N)) {
      mix        <- 0.5 * (pip_mat[, i] + pip_mat[, j]) + eps
      jsd[i, j]  <- jsd[j, i] <- max(-sum(mix * log(mix)) - 0.5 * (H[i] + H[j]), 0)
    }
  }
  clusters <- cutree(hclust(as.dist(jsd), method = "complete"), h = jsd_threshold)
  K <- max(clusters)
  w <- numeric(N)
  for (k in seq_len(K)) {
    idx    <- which(clusters == k)
    ew     <- exp(elbos[idx] - max(elbos[idx]))
    w[idx] <- ew / sum(ew) / K
  }
  as.vector(pip_mat %*% w)
}

#' Select run_ids for a subsample of size n_ens from run_meta.
#'
#' @param run_meta Tibble with columns run_id and lever columns
#'   (restart_id, refine_step, c_value, sigma_0_2_scalar).
#' @param n_ens Number of fits to select.
#' @param lever_type One of "restart", "refine", "sigma_0_2", "c_grid".
#' @param rep_i Random-seed index (used for restart lever; default 1L).
#' @keywords internal
select_subsample_run_ids <- function(run_meta, n_ens, lever_type, rep_i = 1L) {
  n_sel <- min(n_ens, nrow(run_meta))
  switch(lever_type,
    restart = {
      run_meta %>%
        dplyr::arrange(.data$restart_id, .data$run_id) %>%
        dplyr::slice_head(n = n_sel) %>%
        dplyr::pull(run_id)
    },
    refine = {
      run_meta %>%
        dplyr::arrange(.data$refine_step) %>%
        dplyr::slice_head(n = n_sel) %>%
        dplyr::pull(run_id)
    },
    sigma_0_2 = {
      vals <- sort(unique(run_meta$sigma_0_2_scalar))
      keep <- vals[round(seq(1, length(vals), length.out = n_sel))]
      run_meta %>%
        dplyr::filter(.data$sigma_0_2_scalar %in% keep) %>%
        dplyr::pull(run_id)
    },
    c_grid = {
      vals <- sort(unique(run_meta$c_value))
      keep <- vals[round(seq(1, length(vals), length.out = n_sel))]
      run_meta %>%
        dplyr::filter(.data$c_value %in% keep) %>%
        dplyr::pull(run_id)
    },
    stop("Unknown lever_type: ", lever_type)
  )
}

#' Build pre-aggregated confusion matrix from SNPs dataset.
#'
#' Creates a compact summary table with per-bucket causal/non-causal counts
#' at each PIP threshold, grouped by the specified variables. Stores per-bucket
#' counts (NOT cumulative) to allow flexible downstream aggregation.
#'
#' @param snps_dataset_path Path to the aggregated snps_dataset folder (partitioned parquet). Optional when `snp_files` supplied.
#' @param snp_files Optional character vector of parquet file paths (e.g., staging flush outputs) to stream over.
#' @param run_table Data frame with run metadata (must contain run_id and grouping columns).
#' @param grouping_vars Character vector of column names to group by.
#' @param pip_bucket_width Numeric bucket width for PIP binning (default 0.01).
#' @param use_cases_path Optional path to use_cases.csv for label lookup.
#' @param output_path Path where the CSV will be written.
#' @param stream When TRUE, aggregates per file and combines results (recommended when using staging `snp_files`).
#' @return The confusion matrix tibble (invisibly).
#' @export
build_confusion_matrix <- function(snps_dataset_path,
                                   snp_files = NULL,
                                   run_table,
                                   grouping_vars,
                                   pip_bucket_width = 0.01,
                                   use_cases_path = NULL,
                                   output_path,
                                   stream = FALSE,
                                   confusion_bins_path = NULL,
                                   confusion_bins_files = NULL) {
  add_compat_aliases <- function(df) {
    if (is.null(df) || !nrow(df)) return(df)
    if (!"explore_method" %in% names(df) && "variant_type" %in% names(df)) {
      df$explore_method <- df$variant_type
    }
    if (!"variant_type" %in% names(df) && "explore_method" %in% names(df)) {
      df$variant_type <- df$explore_method
    }
    if (!"run_type" %in% names(df) && "init_type" %in% names(df)) {
      df$run_type <- df$init_type
    }
    if (!"init_type" %in% names(df) && "run_type" %in% names(df)) {
      df$init_type <- df$run_type
    }
    if (!"exploration_group" %in% names(df) && "group_key" %in% names(df)) {
      df$exploration_group <- df$group_key
    }
    if (!"group_key" %in% names(df) && "exploration_group" %in% names(df)) {
      df$group_key <- df$exploration_group
    }
    df
  }

  run_table <- add_compat_aliases(run_table)
  if ("run_id" %in% names(run_table)) {
    run_table <- dplyr::mutate(run_table, run_id = as.integer(.data$run_id))
  }

  has_conf_bins <- (!is.null(confusion_bins_path) && file.exists(confusion_bins_path)) ||
    (!is.null(confusion_bins_files) && length(confusion_bins_files))
  if (has_conf_bins) {
    conf_tbl <- if (!is.null(confusion_bins_path) && file.exists(confusion_bins_path)) {
      readr::read_csv(confusion_bins_path, show_col_types = FALSE)
    } else {
      dplyr::bind_rows(lapply(confusion_bins_files, readr::read_csv, show_col_types = FALSE))
    }
    conf_tbl <- add_compat_aliases(conf_tbl)
    if ("run_id" %in% names(conf_tbl)) {
      conf_tbl <- dplyr::mutate(conf_tbl, run_id = as.integer(.data$run_id))
    }

    if ("run_id" %in% names(conf_tbl)) {
      missing_grouping <- setdiff(grouping_vars, names(conf_tbl))
      join_cols <- intersect(missing_grouping, names(run_table))
      if (length(join_cols)) {
        run_join <- run_table %>%
          dplyr::select(.data$run_id, dplyr::all_of(join_cols)) %>%
          dplyr::distinct()
        conf_tbl <- conf_tbl %>%
          dplyr::left_join(run_join, by = "run_id")
        conf_tbl <- add_compat_aliases(conf_tbl)
      }
    }

    available_grouping <- intersect(grouping_vars, names(conf_tbl))
    if (!length(available_grouping)) {
      stop("None of the specified grouping_vars found in confusion_bins data.")
    }
    conf_agg <- conf_tbl %>%
      dplyr::group_by(dplyr::across(dplyr::all_of(available_grouping)), .data$pip_threshold) %>%
      dplyr::summarise(
        n_causal_at_bucket = sum(.data$n_causal_at_bucket, na.rm = TRUE),
        n_noncausal_at_bucket = sum(.data$n_noncausal_at_bucket, na.rm = TRUE),
        .groups = "drop"
      )
    run_id_cols <- intersect(c("run_id", "explore_method", "variant_id", "agg_method"), names(conf_tbl))
    if (!length(run_id_cols)) {
      stop("Cannot compute n_runs from confusion_bins: missing run identity columns.")
    }
    n_runs_per_group <- conf_tbl %>%
      dplyr::distinct(dplyr::across(dplyr::all_of(c(available_grouping, run_id_cols)))) %>%
      dplyr::count(dplyr::across(dplyr::all_of(available_grouping)), name = "n_runs")
    confusion_matrix_table <- conf_agg %>%
      dplyr::left_join(n_runs_per_group, by = available_grouping)
    readr::write_csv(confusion_matrix_table, output_path)
    message("Saved confusion matrix to: ", output_path)
    return(invisible(confusion_matrix_table))
  }

  dataset_available <- !is.null(snps_dataset_path) && dir.exists(snps_dataset_path)
  if (!dataset_available && (is.null(snp_files) || !length(snp_files))) {
    stop("No SNP sources provided: supply snps_dataset_path or snp_files.")
  }
  if (!dataset_available) {
    stream <- TRUE
  }
  if (dataset_available && isTRUE(stream) && (is.null(snp_files) || !length(snp_files))) {
    snp_files <- list.files(snps_dataset_path, pattern = "\\.parquet$", recursive = TRUE, full.names = TRUE)
  }

  # Open SNPs dataset (or derive schema from the first file)
  if (dataset_available) {
    snps_ds <- arrow::open_dataset(snps_dataset_path)
    snps_cols <- names(snps_ds)
  } else {
    snps_cols <- names(arrow::open_dataset(snp_files[[1]], format = "parquet"))
  }

  # Prepare run metadata for grouping
  available_grouping <- intersect(grouping_vars, names(run_table))
  if (length(available_grouping) == 0) {
    stop("None of the specified grouping_vars found in run_table")
  }

  run_map <- run_table %>%
    dplyr::select(run_id, dplyr::all_of(available_grouping)) %>%
    dplyr::mutate(run_id = as.integer(run_id))

  # Optionally add labels
  include_label <- FALSE
  if (!is.null(use_cases_path) && file.exists(use_cases_path) && "use_case_id" %in% available_grouping) {
    use_case_labels <- readr::read_csv(
      use_cases_path,
      show_col_types = FALSE,
      col_select = c(use_case_id, label),
      col_types = readr::cols(
        use_case_id = readr::col_character(),
        label = readr::col_character()
      )
    )
    run_map <- dplyr::left_join(run_map, use_case_labels, by = "use_case_id")
    include_label <- TRUE
  }

  # Count runs per grouping combination
  n_runs_per_group <- run_map %>%
    dplyr::group_by(dplyr::across(dplyr::all_of(available_grouping))) %>%
    dplyr::summarise(n_runs = dplyr::n(), .groups = "drop")

  # Prepare join table (exclude columns already in snps to avoid collision)
  join_keep_cols <- c("run_id", if (include_label) "label" else NULL)
  join_cols_to_exclude <- intersect(setdiff(names(run_map), join_keep_cols), snps_cols)
  if (length(join_cols_to_exclude)) {
    run_map_for_join <- dplyr::select(run_map, -dplyr::all_of(join_cols_to_exclude))
  } else {
    run_map_for_join <- run_map
  }
  run_map_ds <- arrow::InMemoryDataset$create(run_map_for_join)

  # Build aggregation grouping
  agg_group_vars <- available_grouping
  if (include_label) agg_group_vars <- c(agg_group_vars, "label")

  message("Aggregating confusion matrix counts...")
  message("  Grouping by: ", paste(agg_group_vars, collapse = ", "))
  if (dataset_available) {
    message("  Source: ", snps_dataset_path, if (isTRUE(stream)) " (streaming files)" else " (dataset)")
  } else {
    message("  Source: staging parquet files (streaming)")
  }

  aggregate_single_file <- function(path) {
    arrow::open_dataset(path, format = "parquet") %>%
      dplyr::left_join(run_map_ds, by = "run_id") %>%
      dplyr::mutate(
        pip_bucket = pmax(0, pmin(1, floor(pip / pip_bucket_width) * pip_bucket_width))
      ) %>%
      dplyr::group_by(dplyr::across(dplyr::all_of(agg_group_vars)), pip_bucket) %>%
      dplyr::summarise(
        n_snps = dplyr::n(),
        n_causal = sum(causal == 1, na.rm = TRUE),
        n_noncausal = sum(causal == 0, na.rm = TRUE),
        .groups = "drop"
      ) %>%
      dplyr::collect()
  }

  if (isTRUE(stream)) {
    if (is.null(snp_files) || !length(snp_files)) {
      stop("Streaming requested but no snp_files supplied.")
    }
    pb <- utils::txtProgressBar(min = 0, max = length(snp_files), style = 3, initial = 0, width = NA)
    on.exit(if (!is.null(pb)) close(pb), add = TRUE)
    agg_list <- vector("list", length(snp_files))
    for (i in seq_along(snp_files)) {
      message(sprintf("  Streaming %d/%d: %s", i, length(snp_files), basename(snp_files[[i]])))
      agg_list[[i]] <- aggregate_single_file(snp_files[[i]])
      if (!is.null(pb)) utils::setTxtProgressBar(pb, i)
      gc()
    }
    confusion_agg <- dplyr::bind_rows(agg_list)
    rm(agg_list)
    if (nrow(confusion_agg)) {
      confusion_agg <- confusion_agg %>%
        dplyr::group_by(dplyr::across(dplyr::all_of(agg_group_vars)), pip_bucket) %>%
        dplyr::summarise(
          n_snps = sum(n_snps, na.rm = TRUE),
          n_causal = sum(n_causal, na.rm = TRUE),
          n_noncausal = sum(n_noncausal, na.rm = TRUE),
          .groups = "drop"
        )
    }
  } else {
    # Non-streaming: rely on Arrow to aggregate lazily over the dataset
    has_use_case_id <- "use_case_id" %in% snps_cols

    if (has_use_case_id) {
      snps_ds_with_ucid <- snps_ds
    } else {
      subfolders <- list.dirs(snps_dataset_path, recursive = FALSE, full.names = FALSE)
      hive_partitions <- grep("^use_case_id=", subfolders, value = TRUE)

      if (length(hive_partitions) > 0) {
        snps_ds_with_ucid <- arrow::open_dataset(snps_dataset_path, partitioning = "use_case_id")
      } else {
        run_lookup <- run_table %>%
          dplyr::select(run_id, use_case_id) %>%
          dplyr::mutate(run_id = as.integer(run_id)) %>%
          dplyr::distinct()
        run_lookup_ds <- arrow::InMemoryDataset$create(run_lookup)
        snps_ds_with_ucid <- snps_ds %>%
          dplyr::left_join(run_lookup_ds, by = "run_id")
        rm(run_lookup, run_lookup_ds)
      }
    }

    confusion_agg <- snps_ds_with_ucid %>%
      dplyr::left_join(run_map_ds, by = "run_id") %>%
      dplyr::mutate(
        pip_bucket = pmax(0, pmin(1, floor(pip / pip_bucket_width) * pip_bucket_width))
      ) %>%
      dplyr::group_by(dplyr::across(dplyr::all_of(agg_group_vars)), pip_bucket) %>%
      dplyr::summarise(
        n_snps = dplyr::n(),
        n_causal = sum(causal == 1, na.rm = TRUE),
        n_noncausal = sum(causal == 0, na.rm = TRUE),
        .groups = "drop"
      ) %>%
      dplyr::collect()
    rm(snps_ds_with_ucid)
  }

  message("Collected ", nrow(confusion_agg), " rows from Arrow")

  # Build final table with n_runs
  confusion_matrix_table <- confusion_agg %>%
    dplyr::left_join(n_runs_per_group, by = available_grouping) %>%
    dplyr::rename(
      n_causal_at_bucket = n_causal,
      n_noncausal_at_bucket = n_noncausal,
      pip_threshold = pip_bucket
    )

  # Reorder columns nicely
  final_cols <- c(available_grouping,
                  if (include_label) "label" else NULL,
                  "pip_threshold", "n_causal_at_bucket", "n_noncausal_at_bucket", "n_runs")
  confusion_matrix_table <- dplyr::select(confusion_matrix_table, dplyr::all_of(final_cols))

  # Save
  readr::write_csv(confusion_matrix_table, output_path)
  message("Saved confusion matrix to: ", output_path)

  result <- confusion_matrix_table
  rm(confusion_agg, run_map, run_map_for_join, run_map_ds)
  gc()

  invisible(result)
}


# Ensemble scaling analysis helpers ----------------------------------------

# Compute AUPRC from a pooled confusion-bins tibble via trapezoidal PR-curve.
# bins: tibble with pip_threshold, n_causal_at_bucket, n_noncausal_at_bucket
# Returns a scalar or NA_real_.
#' @keywords internal
auprc_from_pooled_bins <- function(bins) {
  if (is.null(bins) || nrow(bins) == 0L) return(NA_real_)
  b       <- bins %>% dplyr::arrange(dplyr::desc(.data$pip_threshold))
  cum_tp  <- cumsum(b$n_causal_at_bucket)
  cum_fp  <- cumsum(b$n_noncausal_at_bucket)
  total_p <- sum(b$n_causal_at_bucket)
  if (total_p == 0L) return(NA_real_)
  keep      <- (cum_tp + cum_fp) > 0L
  precision <- (cum_tp / (cum_tp + cum_fp))[keep]
  recall    <- (cum_tp / total_p)[keep]
  if (length(recall) < 2L) return(NA_real_)
  sum(diff(recall) * 0.5 * (precision[-length(precision)] + precision[-1L]))
}

# Aggregate a p x N PIP matrix to a length-p vector using the specified method.
# pip_mat: numeric matrix, p rows (SNPs), N columns (fits)
# elbos:   length-N numeric (NAs treated as -Inf for max_elbo / softmax)
# method:  "max_elbo", "uniform", "elbo_softmax",
#          "cluster_weight" (JSD threshold 0.15), or "cluster_weight_050" (0.50)
#' @keywords internal
aggregate_pip_matrix <- function(pip_mat, elbos, method, jsd_threshold = 0.15) {
  N <- ncol(pip_mat)
  if (N == 1L) return(as.vector(pip_mat))
  switch(method,
    uniform = rowMeans(pip_mat, na.rm = TRUE),
    max_elbo = {
      best <- which.max(ifelse(is.na(elbos), -Inf, elbos))
      pip_mat[, best, drop = TRUE]
    },
    elbo_softmax = {
      e <- ifelse(is.na(elbos), -Inf, elbos)
      w <- exp(e - max(e)); w <- w / sum(w)
      as.vector(pip_mat %*% w)
    },
    cluster_weight     = .aggregate_cluster_weight(pip_mat, elbos, jsd_threshold),
    cluster_weight_050 = .aggregate_cluster_weight(pip_mat, elbos, 0.50),
    stop("Unknown aggregation method: ", method)
  )
}

# Internal JSD-based cluster-then-ELBO-softmax aggregation.
#' @keywords internal
.aggregate_cluster_weight <- function(pip_mat, elbos, jsd_threshold) {
  N   <- ncol(pip_mat)
  eps <- 1e-10
  H   <- apply(pip_mat + eps, 2L, function(pv) -sum(pv * log(pv)))
  jsd <- matrix(0, N, N)
  for (i in seq_len(N - 1L)) {
    for (j in seq(i + 1L, N)) {
      mix        <- 0.5 * (pip_mat[, i] + pip_mat[, j]) + eps
      jsd[i, j]  <- jsd[j, i] <- max(-sum(mix * log(mix)) - 0.5 * (H[i] + H[j]), 0)
    }
  }
  clusters <- cutree(hclust(as.dist(jsd), method = "complete"), h = jsd_threshold)
  K        <- max(clusters)
  w        <- numeric(N)
  e        <- ifelse(is.na(elbos), -Inf, elbos)
  for (k in seq_len(K)) {
    idx    <- which(clusters == k)
    ew     <- exp(e[idx] - max(e[idx]))
    w[idx] <- ew / sum(ew) / K
  }
  as.vector(pip_mat %*% w)
}

# Select n_ens run_ids from run_meta using the appropriate lever strategy.
# run_meta: tibble with columns run_id, restart_id, refine_step, c_value, sigma_0_2_scalar
# lever_type: "restart", "refine", "sigma_0_2", or "c_grid"
# rep_i: replicate index (seeds RNG for "restart"; ignored for deterministic levers)
#' @keywords internal
select_subsample_run_ids <- function(run_meta, n_ens, lever_type, rep_i = 1L) {
  n_sel <- min(as.integer(n_ens), nrow(run_meta))
  switch(lever_type,
    restart = {
      run_meta %>%
        dplyr::arrange(.data$restart_id, .data$run_id) %>%
        dplyr::slice_head(n = n_sel) %>%
        dplyr::pull(run_id)
    },
    refine = {
      run_meta %>%
        dplyr::arrange(.data$refine_step) %>%
        dplyr::slice_head(n = n_sel) %>%
        dplyr::pull(run_id)
    },
    sigma_0_2 = {
      vals <- sort(unique(run_meta$sigma_0_2_scalar))
      keep <- vals[unique(round(seq(1, length(vals), length.out = n_sel)))]
      run_meta %>%
        dplyr::filter(.data$sigma_0_2_scalar %in% keep) %>%
        dplyr::pull(run_id)
    },
    c_grid = {
      vals <- sort(unique(run_meta$c_value))
      keep <- vals[unique(round(seq(1, length(vals), length.out = n_sel)))]
      run_meta %>%
        dplyr::filter(.data$c_value %in% keep) %>%
        dplyr::pull(run_id)
    },
    stop("Unknown lever_type: ", lever_type)
  )
}

# Accumulate new_bins into pooled_bins[[method]], summing bucket counts.
#' @keywords internal
.accumulate_bins <- function(pooled_bins, method, new_bins) {
  if (is.null(pooled_bins[[method]])) {
    pooled_bins[[method]] <- new_bins
  } else {
    pooled_bins[[method]] <- dplyr::bind_rows(pooled_bins[[method]], new_bins) %>%
      dplyr::group_by(.data$pip_threshold) %>%
      dplyr::summarise(
        n_causal_at_bucket    = sum(.data$n_causal_at_bucket),
        n_noncausal_at_bucket = sum(.data$n_noncausal_at_bucket),
        .groups = "drop"
      )
  }
  pooled_bins
}

# Compute AUPRC + TPR@0.05 for each agg_method from a named bins list.
# Returns tibble(agg_method, AUPRC, tpr05).
#' @keywords internal
.metrics_from_bins_list <- function(pooled_bins, agg_methods) {
  purrr::map_dfr(agg_methods, function(method) {
    pb <- pooled_bins[[method]]
    tibble::tibble(
      agg_method = method,
      AUPRC      = auprc_from_pooled_bins(pb),
      tpr05      = if (!is.null(pb) && nrow(pb) > 0L)
        tpr_at_fpr_threshold(pb) else NA_real_
    )
  })
}

# Build a p x N PIP matrix, causal_vec, and elbos vector from a dataset's
# snps tibble and a vector of selected run_ids. Returns a named list.
#' @keywords internal
.build_pip_mat <- function(ds_snps, sel_ids, elbo_info) {
  run_order <- sort(sel_ids)
  ds_sub    <- ds_snps %>% dplyr::filter(.data$run_id %in% run_order)

  pip_wide <- ds_sub %>%
    dplyr::select(snp_index, run_id, pip) %>%
    tidyr::pivot_wider(
      id_cols     = snp_index,
      names_from  = run_id,
      values_from = pip,
      values_fill = 0,
      names_sort  = TRUE
    ) %>%
    dplyr::arrange(.data$snp_index)

  pip_mat <- as.matrix(pip_wide[, as.character(run_order), drop = FALSE])

  causal_vec <- ds_sub %>%
    dplyr::distinct(.data$snp_index, .data$causal) %>%
    dplyr::arrange(.data$snp_index) %>%
    dplyr::pull(.data$causal)

  elbos <- elbo_info %>%
    dplyr::filter(.data$run_id %in% run_order) %>%
    dplyr::arrange(.data$run_id) %>%
    dplyr::pull(.data$elbo_final)

  list(pip_mat = pip_mat, causal_vec = causal_vec, elbos = elbos)
}

#' Compute post-hoc subscaling curves for pure-lever ensemble specs.
#'
#' For each spec in spec_catalog, subsamples from the full 64-run grid at
#' ensemble sizes n_subsample_sizes, aggregates PIPs, and returns pooled
#' AUPRC + TPR@FPR=0.05 as a function of ensemble size.
#'
#' @param snps_indiv Per-run per-SNP tibble (from snps parquet, individual fits only)
#'   with columns: run_id, snp_index, pip, causal, spec_name, dataset_bundle_id,
#'   annotation_r2, restart_id, refine_step, c_value, sigma_0_2_scalar.
#' @param elbo_info Tibble with columns run_id, elbo_final.
#' @param spec_catalog Tibble with columns spec_name, lever_type, n_random_subsets.
#' @param pip_breaks Numeric vector of PIP bucket break points.
#' @param n_subsample_sizes Integer vector of ensemble sizes to evaluate.
#' @return Tibble: spec_name, annotation_r2 (NA = overall pooled), n_ensemble,
#'   agg_method, AUPRC_mean, AUPRC_se, tpr05_mean, tpr05_se.
#' @export
compute_scaling_curves <- function(snps_indiv, elbo_info, spec_catalog,
                                   pip_breaks,
                                   n_subsample_sizes = c(4L, 8L, 16L, 32L, 64L)) {
  agg_methods <- c("max_elbo", "uniform", "elbo_softmax",
                   "cluster_weight", "cluster_weight_050")

  purrr::map_dfr(seq_len(nrow(spec_catalog)), function(si) {
    spec_nm    <- spec_catalog$spec_name[si]
    lever_type <- spec_catalog$lever_type[si]
    n_reps     <- spec_catalog$n_random_subsets[si]

    snps_spec <- snps_indiv %>% dplyr::filter(.data$spec_name == spec_nm)
    datasets  <- unique(snps_spec$dataset_bundle_id)
    r2_levels <- sort(unique(snps_spec$annotation_r2))

    run_meta_list <- lapply(datasets, function(bid) {
      snps_spec %>%
        dplyr::filter(.data$dataset_bundle_id == bid) %>%
        dplyr::distinct(run_id, restart_id, refine_step, c_value, sigma_0_2_scalar) %>%
        dplyr::left_join(elbo_info, by = "run_id")
    })

    purrr::map_dfr(as.integer(n_subsample_sizes), function(n_ens) {
      rep_rows <- purrr::map_dfr(seq_len(n_reps), function(rep_i) {
        overall_bins <- stats::setNames(
          vector("list", length(agg_methods)), agg_methods)
        r2_bins <- stats::setNames(
          lapply(r2_levels, function(x) {
            stats::setNames(vector("list", length(agg_methods)), agg_methods)
          }),
          as.character(r2_levels)
        )

        for (bid_i in seq_along(datasets)) {
          bid      <- datasets[bid_i]
          run_meta <- run_meta_list[[bid_i]]
          ds_snps  <- snps_spec %>% dplyr::filter(.data$dataset_bundle_id == bid)
          r2_val   <- unique(ds_snps$annotation_r2)[[1L]]
          r2_key   <- as.character(r2_val)

          sel_ids <- select_subsample_run_ids(run_meta, n_ens, lever_type, rep_i)
          if (length(sel_ids) == 0L) next

          pm <- .build_pip_mat(ds_snps, sel_ids, elbo_info)

          for (method in agg_methods) {
            agg_pip      <- aggregate_pip_matrix(pm$pip_mat, pm$elbos, method)
            bins         <- compute_bins_from_pip_vec(agg_pip, pm$causal_vec, pip_breaks)
            overall_bins <- .accumulate_bins(overall_bins, method, bins)
            r2_bins[[r2_key]] <- .accumulate_bins(r2_bins[[r2_key]], method, bins)
          }
        }

        overall_m <- .metrics_from_bins_list(overall_bins, agg_methods) %>%
          dplyr::mutate(annotation_r2 = NA_real_)
        r2_m <- purrr::map_dfr(as.character(r2_levels), function(r2k) {
          .metrics_from_bins_list(r2_bins[[r2k]], agg_methods) %>%
            dplyr::mutate(annotation_r2 = as.numeric(r2k))
        })
        dplyr::bind_rows(overall_m, r2_m) %>% dplyr::mutate(rep_i = rep_i)
      })

      rep_rows %>%
        dplyr::group_by(.data$agg_method, .data$annotation_r2) %>%
        dplyr::summarise(
          AUPRC_mean = mean(.data$AUPRC, na.rm = TRUE),
          AUPRC_se   = if (dplyr::n() > 1L)
            stats::sd(.data$AUPRC, na.rm = TRUE) / sqrt(sum(!is.na(.data$AUPRC)))
          else NA_real_,
          tpr05_mean = mean(.data$tpr05, na.rm = TRUE),
          tpr05_se   = if (dplyr::n() > 1L)
            stats::sd(.data$tpr05, na.rm = TRUE) / sqrt(sum(!is.na(.data$tpr05)))
          else NA_real_,
          .groups = "drop"
        ) %>%
        dplyr::mutate(spec_name = spec_nm, n_ensemble = n_ens)
    })
  }) %>%
    dplyr::select(spec_name, annotation_r2, n_ensemble, agg_method,
                  AUPRC_mean, AUPRC_se, tpr05_mean, tpr05_se)
}

# Map a 2-lever interaction spec name to its two axis column names.
#' @keywords internal
.get_interaction_axes <- function(spec_name) {
  switch(spec_name,
    "A-RF" = c("restart_id",  "refine_step"),
    "A-RS" = c("restart_id",  "sigma_0_2_scalar"),
    "A-FS" = c("refine_step", "sigma_0_2_scalar"),
    "C-CS" = c("c_value",     "sigma_0_2_scalar"),
    stop("Unknown interaction spec: ", spec_name)
  )
}

#' Compute AUPRC and TPR@0.05 for 2-lever interaction specs at 4x4 vs 8x8 resolution.
#'
#' For each spec, selects R evenly-spaced values on each of the two grid axes,
#' producing R x R = R^2 fits, then aggregates and evaluates.
#'
#' @param snps_indiv Same format as compute_scaling_curves().
#' @param elbo_info Tibble with columns run_id, elbo_final.
#' @param interaction_spec_names Character vector, e.g. c("A-RF","A-RS","A-FS","C-CS").
#' @param pip_breaks Numeric vector of PIP bucket break points.
#' @return Tibble: spec_name, annotation_r2 (NA = overall), resolution ("4x4"/"8x8"),
#'   agg_method, AUPRC, tpr05.
#' @export
compute_interaction_scaling <- function(snps_indiv, elbo_info,
                                        interaction_spec_names, pip_breaks) {
  agg_methods <- c("max_elbo", "uniform", "elbo_softmax",
                   "cluster_weight", "cluster_weight_050")
  resolutions <- c("4x4", "8x8")

  purrr::map_dfr(interaction_spec_names, function(spec_nm) {
    axes      <- .get_interaction_axes(spec_nm)
    snps_spec <- snps_indiv %>% dplyr::filter(.data$spec_name == spec_nm)
    datasets  <- unique(snps_spec$dataset_bundle_id)
    r2_levels <- sort(unique(snps_spec$annotation_r2))

    run_meta_list <- lapply(datasets, function(bid) {
      snps_spec %>%
        dplyr::filter(.data$dataset_bundle_id == bid) %>%
        dplyr::distinct(run_id, dplyr::all_of(axes)) %>%
        dplyr::left_join(elbo_info, by = "run_id")
    })

    purrr::map_dfr(resolutions, function(res) {
      R <- as.integer(sub("x.*", "", res))

      overall_bins <- stats::setNames(
        vector("list", length(agg_methods)), agg_methods)
      r2_bins <- stats::setNames(
        lapply(r2_levels, function(x) {
          stats::setNames(vector("list", length(agg_methods)), agg_methods)
        }),
        as.character(r2_levels)
      )

      for (bid_i in seq_along(datasets)) {
        bid      <- datasets[bid_i]
        run_meta <- run_meta_list[[bid_i]]
        ds_snps  <- snps_spec %>% dplyr::filter(.data$dataset_bundle_id == bid)
        r2_val   <- unique(ds_snps$annotation_r2)[[1L]]
        r2_key   <- as.character(r2_val)

        ax1_vals <- sort(unique(run_meta[[axes[[1L]]]]))
        ax2_vals <- sort(unique(run_meta[[axes[[2L]]]]))

        keep1 <- ax1_vals[unique(round(
          seq(1, length(ax1_vals), length.out = min(R, length(ax1_vals)))))]
        keep2 <- ax2_vals[unique(round(
          seq(1, length(ax2_vals), length.out = min(R, length(ax2_vals)))))]

        sel_ids <- run_meta %>%
          dplyr::filter(
            .data[[axes[[1L]]]] %in% keep1,
            .data[[axes[[2L]]]] %in% keep2
          ) %>%
          dplyr::pull(run_id)

        if (length(sel_ids) == 0L) next

        pm <- .build_pip_mat(ds_snps, sel_ids, elbo_info)

        for (method in agg_methods) {
          agg_pip      <- aggregate_pip_matrix(pm$pip_mat, pm$elbos, method)
          bins         <- compute_bins_from_pip_vec(agg_pip, pm$causal_vec, pip_breaks)
          overall_bins <- .accumulate_bins(overall_bins, method, bins)
          r2_bins[[r2_key]] <- .accumulate_bins(r2_bins[[r2_key]], method, bins)
        }
      }

      overall_m <- .metrics_from_bins_list(overall_bins, agg_methods) %>%
        dplyr::mutate(spec_name = spec_nm, resolution = res, annotation_r2 = NA_real_)
      r2_m <- purrr::map_dfr(as.character(r2_levels), function(r2k) {
        .metrics_from_bins_list(r2_bins[[r2k]], agg_methods) %>%
          dplyr::mutate(spec_name = spec_nm, resolution = res,
                        annotation_r2 = as.numeric(r2k))
      })
      dplyr::bind_rows(overall_m, r2_m)
    })
  }) %>%
    dplyr::select(spec_name, annotation_r2, resolution, agg_method, AUPRC, tpr05)
}
