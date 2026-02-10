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
  }) %>%
    dplyr::arrange(.data$task_id, .data$flush_id, .data$type)
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

  # Resolve run_table_path: check temp first, then run_history

if (is.null(run_table_path)) {
    temp_path <- file.path(output_root, "temp", job_name, "run_table.csv")
    history_path <- file.path(output_root, "run_history", job_name, parent_job_id, "run_table.csv")
    if (file.exists(temp_path)) {
      run_table_path <- temp_path
    } else if (file.exists(history_path)) {
      run_table_path <- history_path
    }
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
  snp_files <- dplyr::filter(idx, .data$type == "snps")$path

  if (length(model_files)) {
    model_tbl <- read_csv_safe(model_files, idx)
    readr::write_csv(model_tbl, file.path(output_dir, "model_metrics.csv"))
    log_progress("Wrote model_metrics.csv")
    rm(model_tbl, model_files)
    gc()
  }
  if (length(effect_files)) {
    effect_tbl <- read_csv_safe(effect_files, idx)
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
    conf_tbl <- read_csv_safe(confusion_files, idx)
    readr::write_csv(conf_tbl, file.path(output_dir, "confusion_bins.csv"))
    log_progress("Wrote confusion_bins.csv")
    rm(conf_tbl, confusion_files)
    gc()
  }
  if (length(dataset_files)) {
    dataset_tbl <- read_csv_safe(dataset_files, idx)
    readr::write_csv(dataset_tbl, file.path(output_dir, "dataset_metrics.csv"))
    log_progress("Wrote dataset_metrics.csv")
    rm(dataset_tbl, dataset_files)
    gc()
  }
  if (length(multimodal_files)) {
    multimodal_tbl <- read_csv_safe(multimodal_files, idx)
    readr::write_csv(multimodal_tbl, file.path(output_dir, "multimodal_metrics.csv"))
    log_progress("Wrote multimodal_metrics.csv")
    rm(multimodal_tbl, multimodal_files)
    gc()
  }
  if (length(snp_files) && isTRUE(write_snps_dataset)) {
    log_progress(sprintf("Writing SNP dataset from %d parquet fragment(s)...", length(snp_files)))
    snp_dataset <- arrow::open_dataset(snp_files, format = "parquet")

    # Check if use_case_id is already in the dataset
    has_use_case_id <- "use_case_id" %in% names(snp_dataset)

    if (!has_use_case_id) {
      # Join with run_table to get use_case_id from run_id
      if (is.null(run_table_path) || !file.exists(run_table_path)) {
        stop(
          "SNP parquet files lack 'use_case_id' and no run_table found to map it. ",
          "Provide run_table_path or ensure run_table.csv exists in temp/ or run_history/."
        )
      }
      log_progress("Joining use_case_id from run_table...")
      run_lookup <- readr::read_csv(run_table_path, show_col_types = FALSE) %>%
        dplyr::select(run_id, use_case_id) %>%
        dplyr::mutate(run_id = as.integer(run_id)) %>%
        dplyr::distinct()
      # Convert to Arrow Table for efficient join
      run_lookup_arrow <- arrow::Table$create(run_lookup)
      # Perform the join lazily via dplyr
      snp_dataset <- snp_dataset %>%
        dplyr::left_join(run_lookup_arrow, by = "run_id")
      rm(run_lookup, run_lookup_arrow)
      gc()
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
  has_conf_bins <- (!is.null(confusion_bins_path) && file.exists(confusion_bins_path)) ||
    (!is.null(confusion_bins_files) && length(confusion_bins_files))
  if (has_conf_bins) {
    conf_tbl <- if (!is.null(confusion_bins_path) && file.exists(confusion_bins_path)) {
      readr::read_csv(confusion_bins_path, show_col_types = FALSE)
    } else {
      dplyr::bind_rows(lapply(confusion_bins_files, readr::read_csv, show_col_types = FALSE))
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
    n_runs_per_group <- conf_tbl %>%
      dplyr::group_by(dplyr::across(dplyr::all_of(available_grouping))) %>%
      dplyr::summarise(n_runs = dplyr::n_distinct(.data$run_id, .data$variant_type, .data$variant_id, .data$agg_method), .groups = "drop")
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
