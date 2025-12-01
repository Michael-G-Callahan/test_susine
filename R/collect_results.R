# Task/flush staging collection helpers -----------------------------------

#' Return the staging directory used by task-level flush outputs.
#' @keywords internal
staging_base_dir <- function(job_name,
                             parent_job_id,
                             output_root = "output") {
  file.path(output_root, "slurm_output", job_name, parent_job_id, "combined", "staging")
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
#' @return List with aggregated paths and validation summary.
#' @export
aggregate_staging_outputs <- function(job_name,
                                      parent_job_id,
                                      output_root = "output",
                                      validate = TRUE,
                                      output_dir = NULL) {
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
    output_dir <- file.path(results_dir, "combined", "aggregated")
  }
  ensure_dir(output_dir)

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
    purrr::map_dfr(paths, function(p) {
      df <- readr::read_csv(p, show_col_types = FALSE)
      meta <- dplyr::filter(idx_tbl, .data$path == !!p)
      if (!"task_id" %in% names(df)) df$task_id <- meta$task_id[[1]]
      if (!"flush_id" %in% names(df)) df$flush_id <- meta$flush_id[[1]]
      df
    })
  }

  model_files <- dplyr::filter(idx, .data$type == "model_metrics")$path
  effect_files <- dplyr::filter(idx, .data$type == "effect_metrics")$path
  validation_files <- dplyr::filter(idx, .data$type == "validation")$path
  snp_files <- dplyr::filter(idx, .data$type == "snps")$path

  if (length(model_files)) {
    model_tbl <- read_csv_safe(model_files, idx)
    readr::write_csv(model_tbl, file.path(output_dir, "model_metrics.csv"))
    log_progress("Wrote model_metrics.csv")
  }
  if (length(effect_files)) {
    effect_tbl <- read_csv_safe(effect_files, idx)
    readr::write_csv(effect_tbl, file.path(output_dir, "effect_metrics.csv"))
    log_progress("Wrote effect_metrics.csv")
  }
  if (length(validation_files)) {
    validation_tbl <- read_csv_safe(validation_files, idx)
    readr::write_csv(validation_tbl, file.path(output_dir, "validation.csv"))
    log_progress("Wrote validation.csv")
  }
  if (length(snp_files)) {
    log_progress(sprintf("Writing SNP dataset from %d parquet fragment(s)...", length(snp_files)))
    snp_dataset <- arrow::open_dataset(snp_files, format = "parquet")
    snp_out_dir <- file.path(output_dir, "snps_dataset")
    if (dir.exists(snp_out_dir)) {
      unlink(snp_out_dir, recursive = TRUE)
    }
    arrow::write_dataset(snp_dataset, path = snp_out_dir, format = "parquet")
    log_progress("Finished writing snps_dataset")
  }

  list(
    output_dir = output_dir,
    validation = validation,
    file_index = idx
  )
}
