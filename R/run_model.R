# Task execution and modelling ----------------------------------------------

#' Run all simulation jobs assigned to a SLURM array task.
#'
#' @param job_name Character scalar.
#' @param task_id Integer SLURM array index (1-based).
#' @param job_root Root directory that houses `run_history` and outputs.
#' @param config_path Optional explicit path to `job_config.json`.
#' @param quiet Suppress console output when TRUE.
#'
#' @return Invisible tibble summarising run outcomes.
#' @export
run_task <- function(job_name,
                     task_id,
                     job_root = "output",
                     config_path = NULL,
                     quiet = FALSE) {
  cfg_path <- config_path %||% file.path(job_root, "run_history", job_name, "job_config.json")
  job_config <- load_job_config(cfg_path)
  runs_tbl <- job_config$tables$runs
  task_id <- as.integer(task_id)
  task_runs <- dplyr::filter(runs_tbl, task_id == !!task_id)
  if (!nrow(task_runs)) {
    stop(sprintf("Task %s not found in job config.", task_id))
  }
  verbose_output <- isTRUE(job_config$job$verbose_file_output)
  if (!verbose_output) {
    job_config$.__shard_buffers__ <- new.env(parent = emptyenv())
    job_config$.__buffer_base__ <- NULL
    on.exit(flush_shard_buffers(job_config), add = TRUE)
  }

  if (!quiet) {
    message(sprintf("Running job '%s' task %s with %d run(s).", job_name, task_id, nrow(task_runs)))
  }

  results <- purrr::map_dfr(seq_len(nrow(task_runs)), function(i) {
    run_row <- task_runs[i, , drop = FALSE]
    execute_single_run(run_row, job_config, quiet = quiet)
  })

  invisible(results)
}

#' Load a job configuration JSON file and coerce tibbles.
#' @keywords internal
load_job_config <- function(path) {
  cfg <- jsonlite::read_json(path, simplifyVector = TRUE)
  runs_tbl <- tibble::as_tibble(cfg$tables$runs)
  if (!"matrix_id" %in% names(runs_tbl)) {
    alt_cols <- c("matrix_id.x", "matrix_id.y")
    for (nm in alt_cols) {
      if (nm %in% names(runs_tbl)) {
        runs_tbl$matrix_id <- runs_tbl[[nm]]
        break
      }
    }
  }
  runs_tbl <- dplyr::select(runs_tbl, -dplyr::any_of(c("matrix_id.x", "matrix_id.y")))
  cfg$tables$runs <- runs_tbl
  cfg$tables$scenarios <- tibble::as_tibble(cfg$tables$scenarios)
  cfg$tables$tasks <- tibble::as_tibble(cfg$tables$tasks)
  cfg$tables$use_cases <- tibble::as_tibble(cfg$tables$use_cases)

  repo_override <- Sys.getenv("SUSINE_REPO_ROOT", unset = "")
  if (nzchar(repo_override) && dir.exists(repo_override)) {
    cfg$paths$repo_root <- normalizePath(repo_override, winslash = "/", mustWork = TRUE)
  } else {
    stored_root <- cfg$paths$repo_root
    if (is.null(stored_root) || !dir.exists(stored_root)) {
      cfg$paths$repo_root <- normalizePath(getwd(), winslash = "/", mustWork = FALSE)
    }
  }

  cfg
}

#' Retrieve metadata row for a given use case ID.
#' @keywords internal
lookup_use_case <- function(job_config, use_case_id) {
  dplyr::filter(job_config$tables$use_cases, use_case_id == !!use_case_id)
}

#' Execute a single run row: simulate data, fit model, and persist outputs.
#' @keywords internal
execute_single_run <- function(run_row, job_config, quiet = FALSE) {
  run_row <- tibble::as_tibble(run_row)
  use_case <- lookup_use_case(job_config, run_row$use_case_id)
  if (!nrow(use_case)) {
    stop(sprintf("Unknown use_case_id '%s'.", run_row$use_case_id))
  }

  data_bundle <- generate_data_for_run(run_row, job_config)

  model_result <- run_use_case(
    use_case = use_case,
    run_row = run_row,
    data_bundle = data_bundle,
    job_config = job_config
  )

  eval_res <- evaluate_model(
    fit = model_result$fit,
    X = data_bundle$X,
    y = data_bundle$y,
    causal_idx = data_bundle$causal_idx,
    rho = job_config$job$credible_set_rho %||% 0.95,
    purity_threshold = job_config$job$purity_threshold %||% 0.5,
    compute_curves = FALSE
  )

  write_run_outputs(
    run_row = run_row,
    job_config = job_config,
    evaluation = eval_res,
    data_bundle = data_bundle,
    model_result = model_result
  )

  if (!quiet) {
    summary_msg <- sprintf(
      "run_id=%s use_case=%s seed=%s L=%s | power=%.3f size=%.2f purity=%.2f",
      run_row$run_id,
      run_row$use_case_id,
      run_row$seed,
      run_row$L,
      eval_res$model_filtered$power,
      eval_res$model_filtered$mean_size,
      eval_res$model_filtered$mean_purity
    )
    message(summary_msg)
  }

  tibble::tibble(
    run_id = run_row$run_id,
    task_id = run_row$task_id,
    use_case_id = run_row$use_case_id,
    seed = run_row$seed,
    power = eval_res$model_filtered$power,
    mean_size = eval_res$model_filtered$mean_size,
    mean_purity = eval_res$model_filtered$mean_purity
  )
}

#' Generate data bundle for a given run_row based on scenario.
#' @keywords internal
generate_data_for_run <- function(run_row, job_config) {
  scenario <- run_row$data_scenario
  if (scenario %in% c("simulation_n3", "sim_n3")) {
    spec <- list(
      seed = as.integer(run_row$seed),
      p_star = as.integer(run_row$p_star),
      y_noise = as.numeric(run_row$y_noise),
      annotation_r2 = as.numeric(run_row$annotation_r2),
      inflate_match = as.numeric(run_row$inflate_match),
      gamma_shrink = as.numeric(run_row$gamma_shrink),
      effect_sd = run_row[["effect_sd"]] %||% 1
    )
    return(generate_simulation_data(spec))
  }

  matrix_row <- dplyr::filter(
    job_config$tables$data_matrices,
    .data$matrix_id == !!run_row$matrix_id,
    .data$data_scenario == !!scenario
  )
  if (!nrow(matrix_row)) {
    stop("No matrix metadata found for scenario '", scenario, "' (matrix_id=", run_row$matrix_id, ").")
  }
  repo_root <- job_config$paths$repo_root %||% getwd()
  X_mat <- load_sampled_matrix(matrix_row, repo_root = repo_root)
  spec <- list(
    seed = as.integer(run_row$seed),
    p_star = as.integer(run_row$p_star),
    y_noise = as.numeric(run_row$y_noise),
    annotation_r2 = as.numeric(run_row$annotation_r2),
    inflate_match = as.numeric(run_row$inflate_match),
    gamma_shrink = as.numeric(run_row$gamma_shrink),
    effect_sd = run_row[["effect_sd"]] %||% 1
  )
  generate_simulation_data(spec, base_X = X_mat)
}

#' Fit a use case with the susine backend, handling optional annealing/model averaging.
#' @keywords internal
run_use_case <- function(use_case, run_row, data_bundle, job_config) {
  L <- as.integer(run_row$L)
  mu_strategy <- use_case$mu_strategy[[1]]
  sigma_strategy <- use_case$sigma_strategy[[1]]

  mu_0 <- if (mu_strategy %in% c("functional", "eb_mu")) data_bundle$mu_0 else 0
  sigma_0_2 <- if (sigma_strategy %in% c("functional", "eb_sigma")) data_bundle$sigma_0_2 else NULL
  pi_weights <- NULL

  extra <- use_case$extra_compute[[1]]
  if (length(extra) && is.na(extra)) {
    extra <- NULL
  }
  anneal <- job_config$job$compute$anneal
  model_avg <- job_config$job$compute$model_average

  base_args <- list(
    L = L,
    X = data_bundle$X,
    y = data_bundle$y,
    mu_0 = mu_0,
    sigma_0_2 = sigma_0_2,
    prior_inclusion_weights = pi_weights,
    prior_update_method = use_case$prior_update_method[[1]],
    verbose = FALSE
  )
  if (resolve_flag(use_case$auto_scale_mu[[1]], FALSE)) {
    base_args$auto_scale_mu_0 <- TRUE
  }
  if (resolve_flag(use_case$auto_scale_sigma[[1]], FALSE)) {
    base_args$auto_scale_sigma_0_2 <- TRUE
  }

  if (identical(extra, "anneal")) {
    base_args$anneal_start_T <- anneal$anneal_start_T %||% 5
    base_args$anneal_schedule_type <- anneal$anneal_schedule_type %||% "geometric"
    base_args$anneal_burn_in <- anneal$anneal_burn_in %||% 5
  }

  if (!identical(extra, "model_avg")) {
    fit <- do.call(susine::susine, base_args)
    fit$settings$L <- fit$settings$L %||% L
    fit$metadata <- list(
      use_case_id = use_case$use_case_id[[1]],
      label = use_case$label[[1]]
    )
    return(list(fit = fit, extra = NULL))
  }

  n_inits <- model_avg$n_inits %||% 5
  init_sd <- model_avg$init_sd %||% 0.05
  subfits <- vector("list", n_inits)
  pips_mat <- matrix(0, nrow = n_inits, ncol = length(data_bundle$beta))
  coef_mat <- matrix(0, nrow = n_inits, ncol = length(data_bundle$beta))

  for (i in seq_len(n_inits)) {
    base_args$init_random <- TRUE
    base_args$init_seed <- as.integer(run_row$seed) + i
    base_args$init_random_sd <- init_sd
    subfits[[i]] <- do.call(susine::susine, base_args)
    pips_mat[i, ] <- subfits[[i]]$model_fit$PIPs
    coef_mat[i, ] <- subfits[[i]]$model_fit$coef
  }

  agg_fit <- subfits[[1]]
  agg_fit$model_fit$PIPs <- colMeans(pips_mat)
  agg_fit$model_fit$coef <- colMeans(coef_mat)
  agg_fit$settings$L <- agg_fit$settings$L %||% L
  agg_fit$model_fit$std_coef <- colMeans(vapply(subfits, function(sf) sf$model_fit$std_coef, numeric(ncol(coef_mat))))
  agg_fit$model_fit$fitted_y <- as.vector(data_bundle$X %*% agg_fit$model_fit$coef)
  agg_fit$metadata <- list(
    use_case_id = use_case$use_case_id[[1]],
    label = use_case$label[[1]],
    model_averaging = list(
      n_inits = n_inits,
      init_sd = init_sd,
      seeds = as.integer(run_row$seed) + seq_len(n_inits)
    )
  )

  list(fit = agg_fit, extra = subfits)
}

#' Persist metrics, truth tables, and model fits for a run.
#' @keywords internal
#' Persist metrics, truth tables, and model fits for a run.
#' @keywords internal
write_run_outputs <- function(run_row,
                              job_config,
                              evaluation,
                              data_bundle,
                              model_result) {

  # Compute run-based output directory with sharding by run_id
  run_id_val <- as.integer(run_row$run_id)
  shard_size <- job_config$job$slurm$shard_size_output %||% 1000L
  verbose_output <- isTRUE(job_config$job$verbose_file_output)
  
  # Construct base directory: slurm_output/<job_name>/<parent_id>
  env_parent_id <- Sys.getenv("SUSINE_PARENT_ID", unset = "")
  env_job_name <- Sys.getenv("SUSINE_JOB_NAME", unset = "")
  
  if (nzchar(env_parent_id) && nzchar(env_job_name)) {
    # SLURM environment available
    base_output <- file.path(
      job_config$paths$slurm_output_dir,
      env_job_name,
      env_parent_id
    )
  } else {
    # Fallback for local testing
    base_output <- file.path(job_config$paths$slurm_output_dir, "local_test")
  }
  
  # Compute shard directory
  if (shard_size > 0) {
    shard_idx <- (run_id_val - 1L) %/% as.integer(shard_size)
    shard_dir <- sprintf("shard-%03d", shard_idx)
    run_dir <- file.path(base_output, shard_dir, sprintf("run-%05d", run_id_val))
  } else {
    run_dir <- file.path(base_output, sprintf("run-%05d", run_id_val))
  }
  
  if (verbose_output) {
    ensure_dir(run_dir)
  } else {
    shard_dir <- dirname(run_dir)
    ensure_dir(shard_dir)
  }

  model_metrics <- dplyr::mutate(
    dplyr::bind_rows(
      dplyr::mutate(evaluation$model_unfiltered, filtering = "unfiltered"),
      dplyr::mutate(evaluation$model_filtered, filtering = "purity_filtered")
    ),
    run_id = run_row$run_id,
    task_id = run_row$task_id,
    use_case_id = run_row$use_case_id,
    seed = run_row$seed
  )
  if (verbose_output) {
    readr::write_csv(model_metrics, file.path(run_dir, "model_metrics.csv"))
  }

  format_effects <- function(effects_df, label) {
    if (is.null(effects_df) || !nrow(effects_df)) return(NULL)
    indices <- vapply(effects_df$indices, function(idx) paste(idx, collapse = " "), character(1))
    dplyr::mutate(
      dplyr::select(effects_df, -indices),
      indices = indices,
      filtering = label,
      run_id = run_row$run_id,
      task_id = run_row$task_id,
      use_case_id = run_row$use_case_id,
      seed = run_row$seed
    )
  }

  effect_metrics <- dplyr::bind_rows(
    format_effects(evaluation$effects_unfiltered, "unfiltered"),
    format_effects(evaluation$effects_filtered, "purity_filtered")
  )
  if (!is.null(effect_metrics) && nrow(effect_metrics) && verbose_output) {
    readr::write_csv(effect_metrics, file.path(run_dir, "effect_metrics.csv"))
  }

  snp_tbl <- build_snp_table(
    run_row = run_row,
    data_bundle = data_bundle,
    evaluation = evaluation
  )
  if (verbose_output) {
    write_snps_parquet(run_dir, snp_tbl)
  }

  if (verbose_output && should_write_legacy_snp_csv(job_config)) {
    pip_tbl <- tibble::tibble(
      snp_index = seq_along(evaluation$combined_pip),
      pip = evaluation$combined_pip,
      run_id = run_row$run_id,
      task_id = run_row$task_id,
      use_case_id = run_row$use_case_id,
      seed = run_row$seed
    )
    readr::write_csv(pip_tbl, file.path(run_dir, "pip.csv"))

    truth_tbl <- tibble::tibble(
      snp_index = seq_along(data_bundle$beta),
      beta = data_bundle$beta,
      mu_0 = data_bundle$mu_0,
      sigma_0_2 = data_bundle$sigma_0_2,
      prior_weight = data_bundle$prior_inclusion_weights,
      causal = as.integer(seq_along(data_bundle$beta) %in% data_bundle$causal_idx)
    )
    readr::write_csv(truth_tbl, file.path(run_dir, "truth.csv"))
  }

  if (verbose_output) {
    saveRDS(model_result$fit, file.path(run_dir, "fit.rds"))
    if (!is.null(model_result$extra)) {
      saveRDS(model_result$extra, file.path(run_dir, "subfits.rds"))
    }
  } else {
    if (is.null(job_config$.__buffer_base__)) {
      job_config$.__buffer_base__ <- base_output
    }
    buffer_shard_outputs(
      job_config = job_config,
      shard_label = sprintf("shard-%03d", shard_idx),
      model_metrics = model_metrics,
      effect_metrics = effect_metrics,
      snp_tbl = snp_tbl,
      validation_row = build_validation_row(run_row, shard_idx)
    )
  }
}

build_snp_table <- function(run_row,
                            data_bundle,
                            evaluation) {
  pip <- evaluation$combined_pip
  beta <- data_bundle$beta
  n <- length(pip)
  if (length(beta) != n) {
    stop("combined PIP length does not match beta length for run ", run_row$run_id)
  }

  safe_vec <- function(x, fill = NA_real_) {
    if (is.null(x)) {
      rep(fill, n)
    } else {
      x
    }
  }

  base_tbl <- tibble::tibble(
    snp_index = seq_len(n),
    pip = pip,
    beta = beta,
    mu_0 = safe_vec(data_bundle$mu_0),
    sigma_0_2 = safe_vec(data_bundle$sigma_0_2),
    prior_weight = safe_vec(data_bundle$prior_inclusion_weights),
    causal = as.integer(seq_len(n) %in% data_bundle$causal_idx)
  )

  metadata_cols <- c(
    "run_id",
    "task_id",
    "use_case_id",
    "seed",
    "data_scenario",
    "matrix_id",
    "L",
    "p_star",
    "y_noise",
    "annotation_r2",
    "inflate_match",
    "gamma_shrink"
  )
  present <- intersect(metadata_cols, names(run_row))
  if (!length(present)) {
    return(base_tbl)
  }
  meta_row <- tibble::as_tibble(run_row[, present, drop = FALSE])
  meta_rep <- meta_row[rep(1, n), , drop = FALSE]
  dplyr::bind_cols(meta_rep, base_tbl)
}

write_snps_parquet <- function(run_dir, snp_tbl) {
  snp_path <- file.path(run_dir, "snps.parquet")
  arrow::write_parquet(snp_tbl, snp_path, compression = "zstd")
}

build_validation_row <- function(run_row, shard_idx) {
  tibble::tibble(
    run_id = run_row$run_id,
    shard_idx = shard_idx,
    has_issues = FALSE,
    issues = NA_character_
  )
}

buffer_shard_outputs <- function(job_config,
                                 shard_label,
                                 model_metrics,
                                 effect_metrics,
                                 snp_tbl,
                                 validation_row) {
  buf_env <- job_config$.__shard_buffers__
  if (is.null(buf_env)) {
    buf_env <- new.env(parent = emptyenv())
    job_config$.__shard_buffers__ <- buf_env
  }
  buf <- buf_env[[shard_label]]
  if (is.null(buf)) {
    buf <- list(
      model = list(),
      effect = list(),
      snps = list(),
      validation = list()
    )
  }
  if (!is.null(model_metrics) && nrow(model_metrics)) {
    buf$model[[length(buf$model) + 1L]] <- model_metrics
  }
  if (!is.null(effect_metrics) && nrow(effect_metrics)) {
    buf$effect[[length(buf$effect) + 1L]] <- effect_metrics
  }
  if (!is.null(snp_tbl) && nrow(snp_tbl)) {
    buf$snps[[length(buf$snps) + 1L]] <- snp_tbl
  }
  if (!is.null(validation_row) && nrow(validation_row)) {
    buf$validation[[length(buf$validation) + 1L]] <- validation_row
  }
  buf_env[[shard_label]] <- buf
}

flush_shard_buffers <- function(job_config) {
  buf_env <- job_config$.__shard_buffers__
  if (is.null(buf_env) || !length(ls(buf_env))) {
    return(invisible(NULL))
  }
  base_output <- job_config$.__buffer_base__
  if (is.null(base_output)) {
    warning("Buffered shard outputs found but base output directory is missing.")
    return(invisible(NULL))
  }
  shard_labels <- ls(buf_env, all.names = TRUE)
  for (label in shard_labels) {
    buf <- buf_env[[label]]
    model_df <- if (length(buf$model)) dplyr::bind_rows(buf$model) else NULL
    effect_df <- if (length(buf$effect)) dplyr::bind_rows(buf$effect) else NULL
    snp_df <- if (length(buf$snps)) dplyr::bind_rows(buf$snps) else NULL
    validation_df <- if (length(buf$validation)) dplyr::bind_rows(buf$validation) else NULL
    write_shard_outputs(
      base_output = base_output,
      shard_label = label,
      model_metrics = model_df,
      effect_metrics = effect_df,
      snp_tbl = snp_df,
      validation_row = validation_df
    )
    rm(list = label, envir = buf_env)
  }
  job_config$.__shard_buffers__ <- NULL
  job_config$.__buffer_base__ <- NULL
}

write_shard_outputs <- function(base_output,
                                shard_label,
                                model_metrics,
                                effect_metrics,
                                snp_tbl,
                                validation_row) {
  combined_dir <- file.path(base_output, "combined", "shards")
  ensure_dir(combined_dir)
  if (!is.null(model_metrics) && nrow(model_metrics)) {
    append_shard_csv(
      combined_dir,
      shard_label,
      "model_metrics",
      model_metrics
    )
  }
  if (!is.null(effect_metrics) && nrow(effect_metrics)) {
    append_shard_csv(
      combined_dir,
      shard_label,
      "effect_metrics",
      effect_metrics
    )
  }
  if (!is.null(snp_tbl) && nrow(snp_tbl)) {
    append_shard_parquet(
      combined_dir,
      shard_label,
      "snps",
      snp_tbl
    )
  }
  if (!is.null(validation_row) && nrow(validation_row)) {
    append_shard_csv(
      combined_dir,
      shard_label,
      "validation",
      validation_row
    )
  }
}

append_shard_csv <- function(dir, shard_label, prefix, data) {
  if (!nrow(data)) return()
  ensure_dir(dir)
  shard_file <- file.path(dir, sprintf("%s_%s.csv", prefix, shard_label))
  lock_file <- paste0(shard_file, ".lock")
  lock <- filelock::lock(lock_file, exclusive = TRUE, timeout = 60000)
  on.exit(filelock::unlock(lock), add = TRUE)
  exists <- file.exists(shard_file)
  readr::write_csv(
    data,
    shard_file,
    append = exists,
    col_names = !exists
  )
}

append_shard_parquet <- function(dir, shard_label, prefix, data) {
  ensure_dir(dir)
  shard_file <- file.path(dir, sprintf("%s_%s.parquet", prefix, shard_label))
  lock_file <- paste0(shard_file, ".lock")
  lock <- filelock::lock(lock_file, exclusive = TRUE, timeout = 60000)
  on.exit(filelock::unlock(lock), add = TRUE)
  if (file.exists(shard_file)) {
    existing <- arrow::read_parquet(shard_file)
    data <- dplyr::bind_rows(existing, data)
  }
  arrow::write_parquet(data, shard_file, compression = "zstd")
}

should_write_legacy_snp_csv <- function(job_config) {
  isTRUE(job_config$job$write_legacy_snp_csv)
}

resolve_matrix_absolute_path <- function(path, repo_root) {
  if (is.null(path) || is.na(path) || !nzchar(path)) {
    return(NA_character_)
  }
  sanitize <- function(p) gsub("\\\\", "/", p)
  repo_norm <- NA_character_
  if (!is.null(repo_root) && nzchar(repo_root) && dir.exists(repo_root)) {
    repo_norm <- sanitize(normalizePath(repo_root, winslash = "/", mustWork = TRUE))
  }
  path_clean <- sanitize(path)
  candidates <- c(path_clean)

  if (!is.na(repo_norm)) {
    candidates <- c(candidates, file.path(repo_norm, path_clean))
    stripped_drive <- sub("^[A-Za-z]:", "", path_clean)
    stripped_drive <- sub("^/+", "", stripped_drive)
    candidates <- c(candidates, file.path(repo_norm, stripped_drive))
    repo_base <- basename(repo_norm)
    pattern <- paste0(".*/", repo_base, "/")
    if (grepl(pattern, path_clean)) {
      trimmed <- sub(pattern, "", path_clean)
      candidates <- c(candidates, file.path(repo_norm, trimmed))
    }
  }

  candidates <- unique(candidates)
  for (cand in candidates) {
    if (nzchar(cand) && file.exists(cand)) {
      return(normalizePath(cand, winslash = "/", mustWork = TRUE))
    }
  }
  stop("Matrix path not found: ", path)
}

load_sampled_matrix <- function(matrix_row, repo_root) {
  matrix_path <- matrix_row$matrix_path[[1]]
  abs_matrix <- resolve_matrix_absolute_path(matrix_path, repo_root)
  if (is.na(abs_matrix)) {
    stop("Matrix path is missing for scenario ", matrix_row$data_scenario[[1]])
  }
  col_spec <- readr::cols(.default = readr::col_double())
  col_spec$cols$participant_ID <- readr::col_character()
  tbl <- readr::read_tsv(abs_matrix, col_types = col_spec, progress = FALSE)
  if ("participant_ID" %in% names(tbl)) {
    tbl <- dplyr::select(tbl, -participant_ID)
  }
  X <- as.matrix(tbl)
  storage.mode(X) <- "numeric"
  X
}
