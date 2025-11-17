# Shard-level collection helpers -------------------------------------------

#' Format a shard directory name given a zero-based index.
#' @keywords internal
format_shard_dir <- function(shard_idx) {
  sprintf("shard-%03d", as.integer(shard_idx))
}

#' Infer the shard size recorded for a completed job.
#' @keywords internal
infer_job_shard_size <- function(job_name,
                                 parent_job_id,
                                 output_root = "output",
                                 default = 1000L) {
  cfg_path <- file.path(output_root, "run_history", job_name, parent_job_id, "job_config.json")
  if (!file.exists(cfg_path)) {
    return(as.integer(default))
  }
  cfg <- load_job_config(cfg_path)
  size <- cfg$job$slurm$shard_size_output %||% default
  as.integer(size %||% default)
}

#' Compute the zero-based shard index for a run_id.
#' @keywords internal
shard_index_from_run <- function(run_id, shard_size) {
  run_id <- as.integer(run_id)
  if (!is.finite(shard_size) || shard_size <= 0L) {
    return(0L)
  }
  (run_id - 1L) %/% as.integer(shard_size)
}

#' Return the on-disk directory for a run given its run_id.
#' @keywords internal
run_output_dir <- function(run_id,
                           results_dir,
                           shard_size) {
  run_subdir <- sprintf("run-%05d", as.integer(run_id))
  if (is.finite(shard_size) && shard_size > 0L) {
    shard_dir <- format_shard_dir(shard_index_from_run(run_id, shard_size))
    return(file.path(results_dir, shard_dir, run_subdir))
  }
  file.path(results_dir, run_subdir)
}

#' Build a manifest describing which run_ids live in each shard.
#'
#' @param job_name Completed job slug.
#' @param parent_job_id Parent SLURM array/job id.
#' @param output_root Root directory that houses `run_history` and `slurm_output`.
#' @param shard_size Optional shard size override. Defaults to the value stored
#'   in the job configuration (or 1000 if unavailable).
#'
#' @return Tibble with one row per shard containing the run counts and ranges.
#'   The returned tibble has an attribute `shard_size` containing the size
#'   used when computing indices.
#' @export
build_shard_manifest <- function(job_name,
                                 parent_job_id,
                                 output_root = "output",
                                 shard_size = NULL) {
  stopifnot(length(job_name) == 1L, nzchar(job_name))
  stopifnot(length(parent_job_id) == 1L, nzchar(parent_job_id))

  run_table_path <- file.path(output_root, "run_history", job_name, parent_job_id, "run_table.csv")
  if (!file.exists(run_table_path)) {
    stop("Run table not found: ", run_table_path)
  }
  run_table <- readr::read_csv(run_table_path, show_col_types = FALSE)
  if (!nrow(run_table)) {
    stop("Run table is empty for job '", job_name, "' (parent ", parent_job_id, ").")
  }

  shard_size <- shard_size %||% infer_job_shard_size(job_name, parent_job_id, output_root)
  if (!is.finite(shard_size) || shard_size <= 0L) {
    stop("Shard collection requires a positive shard_size_output. Found: ", shard_size)
  }

  manifest <- run_table %>%
    dplyr::mutate(
      shard_idx = shard_index_from_run(run_id, shard_size)
    ) %>%
    dplyr::group_by(shard_idx) %>%
    dplyr::summarise(
      run_count = dplyr::n(),
      run_id_min = min(.data$run_id),
      run_id_max = max(.data$run_id),
      .groups = "drop"
    ) %>%
    dplyr::arrange(.data$shard_idx) %>%
    dplyr::mutate(shard_dir = format_shard_dir(.data$shard_idx)) %>%
    dplyr::relocate(.data$shard_dir, .after = .data$shard_idx)

  attr(manifest, "shard_size") <- as.integer(shard_size)
  manifest
}

#' Validate that a run's expected artifacts exist and appear well-formed.
#' @keywords internal
validate_single_run <- function(run_row,
                                results_dir,
                                shard_size) {
  run_id <- as.integer(run_row$run_id)
  shard_idx <- shard_index_from_run(run_id, shard_size)
  run_dir <- run_output_dir(run_id, results_dir, shard_size)

  issues <- character(0)
  expect_file <- function(path, label) {
    if (!file.exists(path)) {
      issues <<- c(issues, paste0("missing: ", label))
      return(FALSE)
    }
    TRUE
  }

  # model_metrics
  model_metrics_path <- file.path(run_dir, "model_metrics.csv")
  if (expect_file(model_metrics_path, "model_metrics.csv")) {
    mm <- tryCatch(
      readr::read_csv(model_metrics_path, show_col_types = FALSE),
      error = function(e) {
        issues <<- c(issues, sprintf("model_metrics.csv read error: %s", conditionMessage(e)))
        NULL
      }
    )
    if (!is.null(mm) && nrow(mm) != 2) {
      issues <- c(
        issues,
        sprintf("model_metrics.csv has %d rows (expected 2)", nrow(mm))
      )
    }
  }

  # effect_metrics (optional)
  effect_metrics_path <- file.path(run_dir, "effect_metrics.csv")
  if (file.exists(effect_metrics_path)) {
    exp_rows <- as.integer(run_row$L)
    em <- tryCatch(
      readr::read_csv(effect_metrics_path, show_col_types = FALSE),
      error = function(e) {
        issues <<- c(issues, sprintf("effect_metrics.csv read error: %s", conditionMessage(e)))
        NULL
      }
    )
    if (!is.null(em) && "filtering" %in% names(em)) {
      unfiltered_rows <- sum(em$filtering == "unfiltered", na.rm = TRUE)
      filtered_rows <- sum(em$filtering == "purity_filtered", na.rm = TRUE)
      if (!is.na(exp_rows) && unfiltered_rows != exp_rows) {
        issues <- c(
          issues,
          sprintf(
            "effect_metrics.csv has %d unfiltered rows (expected %d)",
            unfiltered_rows,
            exp_rows
          )
        )
      }
      if (!is.na(exp_rows) && filtered_rows > exp_rows) {
        issues <- c(
          issues,
          sprintf(
            "effect_metrics.csv has %d purity_filtered rows (expected 0-%d)",
            filtered_rows,
            exp_rows
          )
        )
      }
    }
  }

  # snp-level outputs (prefer Parquet; allow legacy CSV pair)
  snps_path <- file.path(run_dir, "snps.parquet")
  if (!file.exists(snps_path)) {
    pip_ok <- expect_file(file.path(run_dir, "pip.csv"), "pip.csv")
    truth_ok <- expect_file(file.path(run_dir, "truth.csv"), "truth.csv")
    if (!pip_ok || !truth_ok) {
      issues <- c(issues, "Missing snps.parquet and/or legacy pip/truth CSV outputs")
    }
  }

  # fit file
  expect_file(file.path(run_dir, "fit.rds"), "fit.rds")

  tibble::tibble(
    run_id = run_id,
    shard_idx = shard_idx,
    has_issues = length(issues) > 0L,
    issues = if (length(issues)) paste(issues, collapse = "; ") else NA_character_
  )
}

#' Read a run's effect metrics if present.
#' @keywords internal
read_run_effect_metrics <- function(run_id,
                                    results_dir,
                                    shard_size) {
  run_dir <- run_output_dir(run_id, results_dir, shard_size)
  effect_path <- file.path(run_dir, "effect_metrics.csv")
  if (!file.exists(effect_path)) {
    return(NULL)
  }
  tryCatch(
    readr::read_csv(effect_path, show_col_types = FALSE) %>%
      dplyr::select(-dplyr::any_of("indices")),
    error = function(e) {
      warning(sprintf("Failed to read effect_metrics for run %s: %s", run_id, conditionMessage(e)))
      NULL
    }
  )
}

#' Read a run's model metrics if present.
#' @keywords internal
read_run_model_metrics <- function(run_id,
                                   results_dir,
                                   shard_size) {
  run_dir <- run_output_dir(run_id, results_dir, shard_size)
  model_path <- file.path(run_dir, "model_metrics.csv")
  if (!file.exists(model_path)) {
    return(NULL)
  }
  tryCatch(
    readr::read_csv(model_path, show_col_types = FALSE),
    error = function(e) {
      warning(sprintf("Failed to read model_metrics for run %s: %s", run_id, conditionMessage(e)))
      NULL
    }
  )
}

#' Read SNP-level outputs, preferring the compact Parquet artifact.
#' @keywords internal
read_run_snp_metrics <- function(run_id,
                                 results_dir,
                                 shard_size) {
  run_dir <- run_output_dir(run_id, results_dir, shard_size)
  snp_path <- file.path(run_dir, "snps.parquet")
  if (file.exists(snp_path)) {
    return(arrow::read_parquet(snp_path))
  }

  pip_path <- file.path(run_dir, "pip.csv")
  truth_path <- file.path(run_dir, "truth.csv")
  if (!file.exists(pip_path) || !file.exists(truth_path)) {
    return(NULL)
  }

  pip_tbl <- tryCatch(
    readr::read_csv(pip_path, show_col_types = FALSE),
    error = function(e) {
      warning(sprintf("Failed to read pip.csv for run %s: %s", run_id, conditionMessage(e)))
      NULL
    }
  )
  truth_tbl <- tryCatch(
    readr::read_csv(truth_path, show_col_types = FALSE),
    error = function(e) {
      warning(sprintf("Failed to read truth.csv for run %s: %s", run_id, conditionMessage(e)))
      NULL
    }
  )
  if (is.null(pip_tbl) || is.null(truth_tbl)) {
    return(NULL)
  }
  dplyr::left_join(pip_tbl, truth_tbl, by = "snp_index")
}

#' Collect validation + metrics for a single shard into consolidated CSV files.
#'
#' @param job_name Completed job slug.
#' @param parent_job_id Parent SLURM array/job id.
#' @param shard_index Zero-based shard index to process.
#' @param output_root Root directory containing `run_history` and `slurm_output`.
#' @param shard_size Optional shard size override (defaults to value in job config).
#' @param force When TRUE, overwrite existing shard-level CSVs.
#' @param quiet Suppress progress messages.
#'
#' @return Invisibly returns a list summarising written files.
#' @export
collect_shard_metrics <- function(job_name,
                                  parent_job_id,
                                  shard_index,
                                  output_root = "output",
                                  shard_size = NULL,
                                  force = FALSE,
                                  quiet = FALSE) {
  stopifnot(length(shard_index) == 1L)
  shard_index <- as.integer(shard_index)

  run_history_dir <- file.path(output_root, "run_history", job_name, parent_job_id)
  results_dir <- file.path(output_root, "slurm_output", job_name, parent_job_id)
  if (!dir.exists(run_history_dir)) {
    stop("Run history directory not found: ", run_history_dir)
  }
  if (!dir.exists(results_dir)) {
    stop("Results directory not found: ", results_dir)
  }

  shard_size <- shard_size %||% infer_job_shard_size(job_name, parent_job_id, output_root)
  if (!is.finite(shard_size) || shard_size <= 0L) {
    stop("Shard collection requires a positive shard_size_output. Found: ", shard_size)
  }

  run_table_path <- file.path(run_history_dir, "run_table.csv")
  run_table <- readr::read_csv(run_table_path, show_col_types = FALSE) %>%
    dplyr::mutate(shard_idx = shard_index_from_run(.data$run_id, shard_size))
  shard_runs <- run_table %>%
    dplyr::filter(.data$shard_idx == !!shard_index) %>%
    dplyr::arrange(.data$run_id)

  if (!nrow(shard_runs)) {
    stop(
      sprintf(
        "Shard %03d has no runs (job=%s parent=%s).",
        shard_index,
        job_name,
        parent_job_id
      )
    )
  }

  if (!quiet) {
    message(sprintf("Collecting shard %03d (%d runs)...", shard_index, nrow(shard_runs)))
  }

  shard_dir_name <- format_shard_dir(shard_index)
  shard_dir <- file.path(results_dir, shard_dir_name)
  if (!dir.exists(shard_dir)) {
    stop("Expected shard directory not found: ", shard_dir)
  }

  shards_combined_dir <- file.path(results_dir, "combined", "shards")
  ensure_dir(shards_combined_dir)

  validation_path <- file.path(shards_combined_dir, sprintf("validation_%s.csv", shard_dir_name))
  model_path <- file.path(shards_combined_dir, sprintf("model_metrics_%s.csv", shard_dir_name))
  effect_path <- file.path(shards_combined_dir, sprintf("effect_metrics_%s.csv", shard_dir_name))
  snps_path <- file.path(shards_combined_dir, sprintf("snps_%s.parquet", shard_dir_name))

  if (!force && file.exists(validation_path) && file.exists(model_path) && file.exists(snps_path)) {
    if (!quiet) {
      message(sprintf("Shard %03d already collected; skipping (force = FALSE).", shard_index))
    }
    return(invisible(list(
      shard_index = shard_index,
      skipped = TRUE,
      validation_path = validation_path,
      model_path = model_path,
      effect_path = effect_path,
      snps_path = snps_path
    )))
  }

  validation_list <- vector("list", nrow(shard_runs))
  model_list <- vector("list", nrow(shard_runs))
  effect_list <- vector("list", nrow(shard_runs))
  snps_list <- vector("list", nrow(shard_runs))

  for (i in seq_len(nrow(shard_runs))) {
    run_row <- shard_runs[i, , drop = FALSE]
      validation_list[[i]] <- validate_single_run(run_row, results_dir, shard_size)
      model_list[[i]] <- read_run_model_metrics(run_row$run_id, results_dir, shard_size)
      effect_list[[i]] <- read_run_effect_metrics(run_row$run_id, results_dir, shard_size)
      snps_list[[i]] <- read_run_snp_metrics(run_row$run_id, results_dir, shard_size)
    }

  validation_df <- dplyr::bind_rows(validation_list)
  readr::write_csv(validation_df, validation_path)

  compact_models <- purrr::compact(model_list)
  if (!length(compact_models)) {
    stop(sprintf("Shard %03d has no readable model_metrics.csv files.", shard_index))
  }
  model_df <- dplyr::bind_rows(compact_models)
  readr::write_csv(model_df, model_path)

  compact_effects <- purrr::compact(effect_list)
  if (length(compact_effects)) {
    effect_df <- dplyr::bind_rows(compact_effects)
    readr::write_csv(effect_df, effect_path)
  } else if (file.exists(effect_path) && force) {
    unlink(effect_path)
  }

  compact_snps <- purrr::compact(snps_list)
  if (!length(compact_snps)) {
    stop(sprintf("Shard %03d has no readable SNP-level outputs.", shard_index))
  }
  snps_df <- dplyr::bind_rows(compact_snps)
  arrow::write_parquet(snps_df, snps_path, compression = "zstd")

  if (!quiet) {
    issues <- sum(validation_df$has_issues)
    message(sprintf(
      "Shard %03d complete: %d runs (%d with issues).",
      shard_index,
      nrow(shard_runs),
      issues
    ))
  }

  invisible(list(
    shard_index = shard_index,
    skipped = FALSE,
    validation_path = validation_path,
    model_path = model_path,
    effect_path = if (length(compact_effects)) effect_path else NA_character_,
    snps_path = snps_path
  ))
}

#' Write a SLURM script that runs shard-level collection tasks.
#'
#' @param job_name Completed job slug.
#' @param parent_job_id Parent SLURM array/job id.
#' @param shard_indices Integer vector of shard indices to include.
#' @param output_root Root path with `slurm_output` and `slurm_scripts`.
#' @param collect_script Path to `collect_shard.R`.
#' @param job_label Name for the SLURM job (defaults to `<job_name>_collect`).
#' @param time,mem,cpus_per_task,partition,email Standard SLURM resource knobs.
#' @param HPC When TRUE, include basic module setup for cluster environments.
#' @param quiet Flag forwarded to each shard-collection task.
#' @param force Whether to add `--force TRUE` when invoking shard tasks.
#'
#' @return Invisibly returns the path to the rendered script.
#' @export
write_shard_collection_job <- function(job_name,
                                       parent_job_id,
                                       shard_indices,
                                       output_root = "output",
                                       collect_script,
                                       job_label = paste0(job_name, "_collect"),
                                       time = "01:00:00",
                                       mem = "4G",
                                       cpus_per_task = 1,
                                       partition = NULL,
                                       email = NULL,
                                       HPC = FALSE,
                                       quiet = TRUE,
                                       force = FALSE) {
  stopifnot(length(job_name) == 1L, nzchar(job_name))
  stopifnot(length(parent_job_id) == 1L, nzchar(parent_job_id))
  stopifnot(length(collect_script) == 1L, nzchar(collect_script))
  shard_indices <- sort(unique(as.integer(shard_indices)))
  if (!length(shard_indices)) {
    stop("No shards supplied to write_shard_collection_job().")
  }
  if (!file.exists(collect_script)) {
    stop("collect_script not found: ", collect_script)
  }

  script_dir <- file.path(output_root, "slurm_scripts")
  ensure_dir(script_dir)
  prints_dir <- file.path(output_root, "slurm_prints", job_name, parent_job_id, "collect_shards")
  ensure_dir(prints_dir)

  script_path <- file.path(script_dir, paste0(job_label, ".slurm"))
  script_lines <- render_shard_collection_slurm(
    job_name = job_name,
    parent_job_id = parent_job_id,
    shard_indices = shard_indices,
    output_root = output_root,
    collect_script = collect_script,
    job_label = job_label,
    time = time,
    mem = mem,
    cpus_per_task = cpus_per_task,
    partition = partition,
    email = email,
    HPC = HPC,
    quiet = quiet,
    force = force
  )
  writeLines(script_lines, con = script_path)
  invisible(script_path)
}

#' Render the shard-collection SLURM script.
#' @keywords internal
render_shard_collection_slurm <- function(job_name,
                                          parent_job_id,
                                          shard_indices,
                                          output_root,
                                          collect_script,
                                          job_label,
                                          time,
                                          mem,
                                          cpus_per_task,
                                          partition,
                                          email,
                                          HPC,
                                          quiet,
                                          force) {
  array_range <- sprintf("0-%d", length(shard_indices) - 1L)
  shard_vec <- paste(shard_indices, collapse = " ")

  hpc_setup <- if (isTRUE(HPC)) {
    c(
      "module load r",
      "",
      'export R_LIBS_USER="${R_LIBS_USER:-$HOME/R/x86_64-pc-linux-gnu-library/4.3}"',
      ""
    )
  } else {
    NULL
  }

  email_lines <- NULL
  if (!is.null(email) && nzchar(email)) {
    email_lines <- c(
      sprintf("#SBATCH --mail-user=%s", email),
      "#SBATCH --mail-type=BEGIN,END,FAIL"
    )
  }

  partition_line <- if (!is.null(partition) && nzchar(partition)) {
    sprintf("#SBATCH --partition=%s", partition)
  } else {
    NULL
  }

  c(
    "#!/bin/bash",
    sprintf("#SBATCH --job-name=%s", job_label),
    sprintf("#SBATCH --array=%s", array_range),
    sprintf("#SBATCH --time=%s", time),
    sprintf("#SBATCH --mem=%s", mem),
    sprintf("#SBATCH --cpus-per-task=%s", cpus_per_task),
    partition_line,
    email_lines,
    "#SBATCH --output=/dev/null",
    "#SBATCH --error=/dev/null",
    "",
    "set -euo pipefail",
    "",
    sprintf('JOB_NAME="%s"', job_name),
    sprintf('PARENT_ID="%s"', parent_job_id),
    sprintf('OUTPUT_ROOT="%s"', normalizePath(output_root, winslash = "/", mustWork = FALSE)),
    sprintf('SLURM_PRINTS_BASE="%s"', normalizePath(file.path(output_root, "slurm_prints"), winslash = "/", mustWork = FALSE)),
    sprintf('COLLECT_SCRIPT="%s"', normalizePath(collect_script, winslash = "/", mustWork = FALSE)),
    sprintf('FORCE_FLAG=%d', as.integer(isTRUE(force))),
    sprintf('QUIET_FLAG=%d', as.integer(isTRUE(quiet))),
    sprintf('SHARD_IDS=(%s)', shard_vec),
    "",
    'PRINTS_DIR="${SLURM_PRINTS_BASE}/${JOB_NAME}/${PARENT_ID}/collect_shards"',
    'mkdir -p "${PRINTS_DIR}"',
    "",
    hpc_setup,
    "",
    'TASK_SLOT="${SLURM_ARRAY_TASK_ID}"',
    'SHARD_IDX="${SHARD_IDS[$TASK_SLOT]}"',
    'LOG_LABEL=$(printf "shard-%03d" "${SHARD_IDX}")',
    'exec >"${PRINTS_DIR}/${LOG_LABEL}.out" 2>"${PRINTS_DIR}/${LOG_LABEL}.err"',
    'echo "[$(date -Is)] Collecting ${LOG_LABEL} for job ${JOB_NAME} (parent ${PARENT_ID})"',
    "",
    'CMD=(Rscript "$COLLECT_SCRIPT" --job-name "$JOB_NAME" --parent-id "$PARENT_ID" --shard-index "${SHARD_IDX}" --output-root "$OUTPUT_ROOT")',
    'if [ "${FORCE_FLAG}" -eq 1 ]; then',
    '  CMD+=("--force" "TRUE")',
    'fi',
    'if [ "${QUIET_FLAG}" -eq 1 ]; then',
    '  CMD+=("--quiet" "TRUE")',
    'fi',
    '"${CMD[@]}"',
    'echo "[$(date -Is)] Completed ${LOG_LABEL}"'
  )
}
