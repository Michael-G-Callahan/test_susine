# Real-case SLURM helpers ---------------------------------------------------

#' Build a real-case job config for SLURM array execution.
#'
#' @param job_name Character scalar.
#' @param gene_name Character scalar.
#' @param param_grid Data frame with columns `sigma_0_2` and `annotation_scale`.
#' @param data_path Path to .RData file containing R, z, a, variant_map, n_sample.
#' @param output_root Base output directory (default "output").
#' @param task_max Maximum SLURM array size.
#' @param email Notification email address.
#' @param time SLURM time limit.
#' @param mem SLURM memory limit.
#' @param shuffle_seed Optional seed for reproducible shuffling.
#' @return Job config list with tables and paths.
#' @export
build_real_case_job_config <- function(job_name,
                                       gene_name,
                                       param_grid,
                                       data_path,
                                       L,
                                       output_root = "output",
                                       task_max = 200,
                                       email = "mgc5166@psu.edu",
                                       time = "02:00:00",
                                       mem = "4G",
                                       shuffle_seed = NULL) {
  if (!is.data.frame(param_grid) || !nrow(param_grid)) {
    stop("param_grid must be a non-empty data.frame.")
  }
  if (!all(c("sigma_0_2", "annotation_scale") %in% names(param_grid))) {
    stop("param_grid must include sigma_0_2 and annotation_scale columns.")
  }
  if (length(L) != 1L || is.na(L)) {
    stop("L must be a single non-missing value.")
  }

  gene_label <- tolower(gene_name)
  output_root <- normalizePath(output_root, winslash = "/", mustWork = FALSE)
  job_root <- file.path(output_root, "real_data_analysis", gene_label)
  job_control_dir <- file.path(job_root, "job_control", job_name)

  data_path <- normalizePath(data_path, winslash = "/", mustWork = FALSE)

  runs <- tibble::as_tibble(param_grid) %>%
    dplyr::mutate(
      run_id = dplyr::row_number(),
      sigma_0_2 = as.numeric(.data$sigma_0_2),
      annotation_scale = as.numeric(.data$annotation_scale),
      L = as.integer(L)
    )

  task_info <- assign_real_case_tasks(runs, task_max = task_max, shuffle_seed = shuffle_seed)

  list(
    job = list(
      name = job_name,
      gene_name = gene_name,
      created_at = timestamp_utc(),
      L = as.integer(L),
      task_max = as.integer(task_max),
      slurm = list(
        time = time,
        mem = mem,
        cpus_per_task = 1,
        email = email
      ),
      data_path = data_path
    ),
    paths = list(
      output_root = output_root,
      job_root = job_root,
      job_control_dir = job_control_dir,
      temp_dir = file.path(job_control_dir, "temp"),
      slurm_scripts_dir = file.path(job_control_dir, "slurm_scripts"),
      slurm_prints_dir = file.path(job_control_dir, "slurm_prints")
    ),
    tables = list(
      runs = task_info$runs,
      tasks = task_info$tasks
    )
  )
}

#' Assign tasks to real-case runs with optional shuffling.
#' @keywords internal
assign_real_case_tasks <- function(runs, task_max = 200, shuffle_seed = NULL) {
  n_models <- nrow(runs)
  if (n_models < 1) {
    stop("No models available to assign.")
  }
  n_tasks <- min(as.integer(task_max), n_models)
  runs <- dplyr::mutate(runs, run_id = dplyr::row_number())

  if (!is.null(shuffle_seed)) {
    set.seed(shuffle_seed)
  }
  runs <- runs %>%
    dplyr::mutate(random_id = sample.int(n_models)) %>%
    dplyr::arrange(.data$random_id)

  if (n_models <= n_tasks) {
    runs <- runs %>%
      dplyr::mutate(task_id = dplyr::row_number())
  } else {
    base <- as.integer(round(n_models / n_tasks))
    if (base < 1L) {
      base <- 1L
    }
    if (base * (n_tasks - 1L) >= n_models) {
      base <- as.integer(floor(n_models / n_tasks))
      if (base < 1L) {
        base <- 1L
      }
    }
    remainder <- n_models - base * (n_tasks - 1L)
    sizes <- c(rep(base, n_tasks - 1L), remainder)
    task_id <- rep(seq_len(n_tasks), times = sizes)
    runs$task_id <- task_id
  }

  tasks <- runs %>%
    dplyr::group_by(.data$task_id) %>%
    dplyr::summarise(runs_per_task = dplyr::n(), .groups = "drop") %>%
    dplyr::arrange(.data$task_id)

  list(
    runs = dplyr::arrange(runs, .data$run_id),
    tasks = tasks
  )
}

#' Write real-case job artifacts to disk.
#'
#' @param job_config Output of [build_real_case_job_config()].
#' @param run_task_script Path to task runner script.
#' @return List of artifact paths.
#' @export
write_real_case_job_artifacts <- function(job_config, run_task_script) {
  paths <- job_config$paths
  ensure_dir(paths$job_root)
  ensure_dir(paths$job_control_dir)
  ensure_dir(paths$slurm_scripts_dir)
  ensure_dir(paths$slurm_prints_dir)

  unlink(paths$temp_dir, recursive = TRUE)
  ensure_dir(paths$temp_dir)

  job_json_path <- file.path(paths$temp_dir, "job_config.json")
  job_config_json <- job_config
  job_config_json$tables$runs <- NULL
  jsonlite::write_json(
    job_config_json,
    path = job_json_path,
    auto_unbox = TRUE,
    digits = NA,
    pretty = TRUE
  )

  run_table_path <- file.path(paths$temp_dir, "run_table.csv")
  readr::write_csv(job_config$tables$runs, run_table_path)

  task_path <- file.path(paths$temp_dir, "task_table.csv")
  readr::write_csv(job_config$tables$tasks, task_path)

  slurm_path <- file.path(paths$slurm_scripts_dir, paste0(job_config$job$name, ".slurm"))
  writeLines(
    render_real_case_slurm_script(job_config, run_task_script = run_task_script),
    con = slurm_path
  )

  list(
    job_config = job_json_path,
    run_table = run_table_path,
    tasks = task_path,
    slurm_script = slurm_path
  )
}

#' Render SLURM submission script for real-case runs.
#' @keywords internal
render_real_case_slurm_script <- function(job_config, run_task_script) {
  job <- job_config$job
  paths <- job_config$paths
  tasks <- job_config$tables$tasks
  n_tasks <- nrow(tasks)
  slurm <- job$slurm

  script <- c(
    "#!/bin/bash",
    sprintf("#SBATCH --job-name=%s", job$name),
    sprintf("#SBATCH --array=1-%d", n_tasks),
    sprintf("#SBATCH --time=%s", slurm$time),
    sprintf("#SBATCH --mem=%s", slurm$mem),
    sprintf("#SBATCH --cpus-per-task=%s", slurm$cpus_per_task),
    sprintf("#SBATCH --mail-user=%s", slurm$email),
    "#SBATCH --mail-type=BEGIN,END,FAIL",
    "#SBATCH --output=/dev/null",
    "#SBATCH --error=/dev/null",
    "",
    "set -euo pipefail",
    "",
    sprintf('OUTPUT_BASE="%s"', normalizePath(paths$job_root, winslash = "/", mustWork = FALSE)),
    sprintf('CONFIG_PATH="%s"', normalizePath(file.path(paths$temp_dir, "job_config.json"), winslash = "/", mustWork = FALSE)),
    sprintf('RUN_TASK_SCRIPT="%s"', normalizePath(run_task_script, winslash = "/", mustWork = FALSE)),
    sprintf('SLURM_PRINTS_BASE="%s"', normalizePath(paths$slurm_prints_dir, winslash = "/", mustWork = FALSE)),
    "",
    "# --- identifiers ---",
    'JOBNAME="${SLURM_JOB_NAME}"',
    'PARENT_ID="${SLURM_ARRAY_JOB_ID:-$SLURM_JOB_ID}"',
    'TASK_ID="${SLURM_ARRAY_TASK_ID:-0}"',
    'PRINTS_DIR="${SLURM_PRINTS_BASE}/${JOBNAME}/${PARENT_ID}"',
    'mkdir -p "${PRINTS_DIR}"',
    "",
    'exec >"${PRINTS_DIR}/${TASK_ID}.out" 2>"${PRINTS_DIR}/${TASK_ID}.err"',
    'echo "[$(date -Is)] Starting task ${TASK_ID} for job ${JOBNAME} (parent ${PARENT_ID})"',
    "",
    'Rscript "$RUN_TASK_SCRIPT" \\',
    '  --job-name "$JOBNAME" \\',
    '  --task-id "$TASK_ID" \\',
    '  --output-base "$OUTPUT_BASE" \\',
    '  --config-path "$CONFIG_PATH" \\',
    '  --job-id "$PARENT_ID"',
    "",
    'echo "[$(date -Is)] Completed task ${TASK_ID}"'
  )
  script[!is.na(script)]
}

#' Run a real-case task from a job config.
#'
#' @param job_name Character scalar.
#' @param task_id Integer SLURM array index (1-based).
#' @param output_base Root output path under real_data_analysis/<gene>.
#' @param config_path Path to job_config.json.
#' @param job_id Job identifier (e.g., SLURM array job id).
#' @param quiet Suppress console output when TRUE.
#' @return Invisible list of fit summaries.
#' @export
run_real_case_task <- function(job_name,
                               task_id,
                               output_base,
                               config_path,
                               job_id = "local",
                               quiet = FALSE) {
  if (is.null(config_path) || !file.exists(config_path)) {
    stop("config_path not found: ", config_path)
  }
  cfg <- jsonlite::read_json(config_path, simplifyVector = TRUE)
  runs_tbl <- cfg$tables$runs

  if (is.null(runs_tbl)) {
    run_tbl_path <- file.path(dirname(config_path), "run_table.csv")
    if (!file.exists(run_tbl_path)) {
      stop("Run table not found: ", run_tbl_path)
    }
    runs_tbl <- readr::read_csv(run_tbl_path, show_col_types = FALSE)
  } else {
    runs_tbl <- tibble::as_tibble(runs_tbl)
  }

  task_id <- as.integer(task_id)
  if (is.na(task_id)) {
    stop("task_id must be an integer.")
  }

  task_runs <- dplyr::filter(runs_tbl, .data$task_id == !!task_id)
  if (!nrow(task_runs)) {
    stop(sprintf("Task %s not found in job config.", task_id))
  }

  data_path <- cfg$job$data_path
  if (is.null(data_path) || !file.exists(data_path)) {
    stop("Data inputs not found: ", data_path)
  }
  load(data_path)

  output_dir <- file.path(output_base, job_id)
  ensure_dir(output_dir)

  fit_summaries <- vector("list", nrow(task_runs))
  total_runs <- nrow(task_runs)

  for (i in seq_len(total_runs)) {
    run_row <- task_runs[i, , drop = FALSE]
    sigma_0_2 <- as.numeric(run_row$sigma_0_2)
    annotation_scale <- as.numeric(run_row$annotation_scale)
    L_val <- as.integer(run_row$L %||% cfg$job$L)

    fit <- susine::susine_rss(
      L = L_val,
      z = z,
      R = R,
      n = n_sample,
      sigma_0_2 = sigma_0_2,
      mu_0 = annotation_scale * a,
      prior_update_method = "none"
    )

    fit_summaries[[i]] <- list(
      settings = list(
        gene_name = cfg$job$gene_name,
        sigma_0_2 = sigma_0_2,
        annotation_scale = annotation_scale,
        L = L_val,
        n_sample = n_sample
      ),
      PIPs = fit$model_fit$PIPs,
      elbo_last = tail(fit$model_fit$elbo, 1),
      sigma_2_last = tail(fit$model_fit$sigma_2, 1),
      n_iter = length(fit$model_fit$elbo),
      b_hat_sum = colSums(fit$effect_fits$b_hat)
    )

    if (!quiet) {
      message(sprintf("Task %d: %d/%d runs completed.", task_id, i, total_runs))
    }
  }

  out_path <- file.path(output_dir, sprintf("task-%03d_fit_summaries.rds", task_id))
  saveRDS(fit_summaries, out_path)

  invisible(fit_summaries)
}
