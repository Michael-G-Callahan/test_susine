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
  job_config <- load_job_config(cfg_path, task_id_filter = as.integer(task_id))
  runs_tbl <- job_config$tables$runs
  task_id <- as.integer(task_id)
  task_runs <- dplyr::filter(runs_tbl, task_id == !!task_id)
  if (!nrow(task_runs)) {
    stop(sprintf("Task %s not found in job config.", task_id))
  }
  verbose_output <- isTRUE(job_config$job$verbose_file_output)
  buffer_ctx <- NULL
  if (!verbose_output) {
    buffer_ctx <- new.env(parent = emptyenv())
    buffer_ctx$buffers <- list(
      model = list(),
      effect = list(),
      snps = list(),
      confusion = list(),
      dataset_metrics = list(),
      multimodal = list(),
      validation = list()
    )
    buffer_ctx$base_output <- determine_base_output(job_config)
    buffer_ctx$task_id <- task_id
    buffer_ctx$flush_index <- 0L
  }

  if (!quiet) {
    message(sprintf("Running job '%s' task %s with %d run(s).", job_name, task_id, nrow(task_runs)))
  }

  if (!"dataset_bundle_id" %in% names(task_runs)) {
    stop("run_table is missing dataset_bundle_id; regenerate job config.")
  }
  bundle_ids <- unique(task_runs$dataset_bundle_id)
  results <- purrr::map_dfr(seq_along(bundle_ids), function(i) {
    bundle_id <- bundle_ids[i]
    bundle_runs <- dplyr::filter(task_runs, dataset_bundle_id == !!bundle_id)
    res <- execute_dataset_bundle(bundle_runs, job_config, quiet = quiet, buffer_ctx = buffer_ctx)
    if (!is.null(buffer_ctx)) {
      flush_task_buffers(buffer_ctx)
    }
    res
  })

  if (!is.null(buffer_ctx)) {
    flush_task_buffers(buffer_ctx)
  }

  invisible(results)
}

#' Load a job configuration JSON file and coerce tibbles.
#' When task_id_filter is provided, the runs table is filtered to that task,
#' using the on-disk run_table.csv alongside the config when available.
#' @keywords internal
load_job_config <- function(path, task_id_filter = NULL) {
  cfg <- jsonlite::read_json(path, simplifyVector = TRUE)

  runs_tbl <- cfg$tables$runs

  if (is.null(runs_tbl)) {
    run_tbl_path <- file.path(dirname(path), "run_table.csv")
    if (!file.exists(run_tbl_path)) {
      stop("Run table not found: ", run_tbl_path)
    }
    runs_tbl <- readr::read_csv(run_tbl_path, show_col_types = FALSE)
  } else {
    runs_tbl <- tibble::as_tibble(runs_tbl)
  }

  if ("sigma_0_2_scalar" %in% names(runs_tbl)) {
    runs_tbl$sigma_0_2_scalar <- as.character(runs_tbl$sigma_0_2_scalar)
  }
  if (!"architecture" %in% names(runs_tbl)) {
    runs_tbl$architecture <- "sparse"
  }
  if (!"refine_step" %in% names(runs_tbl)) {
    runs_tbl$refine_step <- NA_integer_
  }
  if (!"prior_spec_id" %in% names(runs_tbl) && "use_case_id" %in% names(runs_tbl)) {
    runs_tbl$prior_spec_id <- runs_tbl$use_case_id
  }

  if (!is.null(task_id_filter)) {
    runs_tbl <- dplyr::filter(runs_tbl, task_id == !!as.integer(task_id_filter))
  }

  if (!"phenotype_seed" %in% names(runs_tbl) && "seed" %in% names(runs_tbl)) {
    runs_tbl <- dplyr::rename(runs_tbl, phenotype_seed = .data$seed)
  }

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
  if (!is.null(cfg$tables$dataset_bundles)) {
    cfg$tables$dataset_bundles <- tibble::as_tibble(cfg$tables$dataset_bundles)
  } else if (!is.null(cfg$tables$scenarios)) {
    cfg$tables$dataset_bundles <- tibble::as_tibble(cfg$tables$scenarios)
  } else {
    bundles_path <- file.path(dirname(path), "dataset_bundles.csv")
    if (file.exists(bundles_path)) {
      cfg$tables$dataset_bundles <- readr::read_csv(bundles_path, show_col_types = FALSE)
    }
  }
  if (!is.null(cfg$tables$data_matrices)) {
    cfg$tables$data_matrices <- tibble::as_tibble(cfg$tables$data_matrices)
  }
  cfg$tables$tasks <- tibble::as_tibble(cfg$tables$tasks)
  cfg$tables$use_cases <- tibble::as_tibble(cfg$tables$use_cases)
  if (!"prior_spec_id" %in% names(cfg$tables$use_cases) && "use_case_id" %in% names(cfg$tables$use_cases)) {
    cfg$tables$use_cases$prior_spec_id <- cfg$tables$use_cases$use_case_id
  }

  if (!is.null(cfg$tables$dataset_bundles)) {
    if (!"phenotype_seed" %in% names(cfg$tables$dataset_bundles) && "seed" %in% names(cfg$tables$dataset_bundles)) {
      cfg$tables$dataset_bundles <- dplyr::rename(cfg$tables$dataset_bundles, phenotype_seed = .data$seed)
    }
    if (!"architecture" %in% names(cfg$tables$dataset_bundles)) {
      cfg$tables$dataset_bundles$architecture <- "sparse"
    }
  }

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
  dplyr::filter(
    job_config$tables$use_cases,
    .data$use_case_id == !!use_case_id | .data$prior_spec_id == !!use_case_id
  )
}

# Compute confusion bins for a PIP vector.
compute_confusion_bins <- function(pip, causal, pip_bucket_width = 0.01) {
  n <- length(pip)
  if (length(causal) != n) {
    stop("pip and causal lengths differ.")
  }
  bins <- pmax(0, pmin(1, floor(pip / pip_bucket_width) * pip_bucket_width))
  bins <- round(bins, 2)
  tibble::tibble(
    pip_threshold = bins,
    causal = as.integer(causal)
  ) %>%
    dplyr::group_by(.data$pip_threshold) %>%
    dplyr::summarise(
      n_causal_at_bucket = sum(.data$causal == 1, na.rm = TRUE),
      n_noncausal_at_bucket = sum(.data$causal == 0, na.rm = TRUE),
      .groups = "drop"
    )
}

prior_cache_key <- function(annotation_r2, inflate_match, annotation_seed = NA_integer_) {
  sprintf("r2=%s|inflate=%s|seed=%s",
          format(annotation_r2, digits = 6, scientific = FALSE),
          format(inflate_match, digits = 6, scientific = FALSE),
          format(annotation_seed, digits = 10, scientific = FALSE))
}

get_priors_cached <- function(cache_env,
                              beta,
                              annotation_r2,
                              inflate_match,
                              base_sigma2,
                              effect_sd,
                              annotation_seed = NA_integer_) {
  key <- prior_cache_key(annotation_r2, inflate_match, annotation_seed)
  if (exists(key, envir = cache_env, inherits = FALSE)) {
    return(get(key, envir = cache_env, inherits = FALSE))
  }
  priors <- simulate_priors(
    beta = beta,
    annotation_r2 = annotation_r2,
    inflate_match = inflate_match,
    base_sigma2 = base_sigma2,
    effect_sd = effect_sd,
    seed = annotation_seed
  )
  assign(key, priors, envir = cache_env)
  priors
}

normalize_blocked_idx <- function(blocked_idx, p) {
  idx <- unique(as.integer(blocked_idx))
  idx <- idx[is.finite(idx) & idx >= 1L & idx <= as.integer(p)]
  sort(idx)
}

resolve_primary_init_type <- function(run_row) {
  restart_id <- suppressWarnings(as.integer(run_row$restart_id %||% NA_integer_))
  init_type <- as.character(run_row$init_type %||% NA_character_)
  if (is.na(init_type) || !nzchar(init_type)) {
    init_type <- if (!is.na(restart_id) && restart_id > 1L) "warm" else "default"
  }
  init_type
}

resolve_restart_seed <- function(run_row, restart_id = NA_integer_) {
  restart_seed <- suppressWarnings(as.integer(run_row$restart_seed %||% NA_integer_))
  if (is.na(restart_seed)) {
    seed_bump <- ifelse(is.na(restart_id), 1L, restart_id)
    restart_seed <- as.integer(run_row$phenotype_seed %||% 1L) + seed_bump
  }
  as.integer(restart_seed)
}

infer_primary_variant_meta <- function(run_row, wall_time_sec = NA_real_) {
  restart_id <- suppressWarnings(as.integer(run_row$restart_id %||% NA_integer_))
  refine_step <- suppressWarnings(as.integer(run_row$refine_step %||% NA_integer_))
  init_type <- resolve_primary_init_type(run_row)
  is_refine <- !is.na(refine_step)
  list(
    variant_type = if (is_refine) "refine" else if (!is.na(restart_id)) "restart" else "base",
    variant_id = if (is_refine) refine_step else if (!is.na(restart_id)) restart_id else NA_integer_,
    init_type = init_type,
    is_restart = !is.na(restart_id),
    is_agg = FALSE,
    wall_time_sec = as.numeric(wall_time_sec)
  )
}

execution_cache_key <- function(use_case,
                                run_row,
                                data_bundle,
                                job_config,
                                blocked_idx = integer(0)) {
  use_case <- tibble::as_tibble(use_case)
  p <- ncol(data_bundle$X)
  L <- as.integer(run_row$L)
  backend <- as.character(use_case$backend[[1]] %||% NA_character_)
  prior_mean_strategy <- as.character(use_case$prior_mean_strategy[[1]] %||% NA_character_)
  prior_variance_strategy <- as.character(use_case$prior_variance_strategy[[1]] %||% NA_character_)
  inclusion_prior_strategy <- as.character(use_case$inclusion_prior_strategy[[1]] %||% NA_character_)
  unmappable_effects <- as.character(use_case$unmappable_effects[[1]] %||% "none")

  fmt_num <- function(x) {
    x <- as.numeric(x)[1]
    if (!is.finite(x)) return("NA")
    format(x, digits = 16, scientific = FALSE, trim = TRUE)
  }
  fmt_chr <- function(x) {
    x <- as.character(x)[1]
    if (is.na(x) || !nzchar(x)) return("NA")
    x
  }
  fmt_int <- function(x) {
    x <- as.integer(x)[1]
    if (!is.finite(x)) return("NA")
    as.character(x)
  }

  init_type <- resolve_primary_init_type(run_row)
  restart_id <- suppressWarnings(as.integer(run_row$restart_id %||% NA_integer_))
  warm_seed <- if (identical(init_type, "warm")) resolve_restart_seed(run_row, restart_id) else NA_integer_
  alpha_conc <- job_config$job$compute$restart$alpha_concentration %||% 1

  annotation_active <- prior_mean_strategy == "functional_mu" || inclusion_prior_strategy == "functional_pi"
  ann_r2_key <- if (annotation_active) fmt_num(run_row$annotation_r2) else "NA"
  inflate_key <- if (annotation_active) fmt_num(run_row$inflate_match) else "NA"
  ann_seed_key <- if (annotation_active) fmt_int(run_row$annotation_seed) else "NA"

  c_value <- as.numeric(run_row$c_value %||% NA_real_)
  if (!is.finite(c_value)) c_value <- 1
  c_key <- if (prior_mean_strategy == "functional_mu") fmt_num(c_value) else "NA"

  tau_value <- as.numeric(run_row$tau_value %||% NA_real_)
  if (!is.finite(tau_value) || tau_value <= 0) tau_value <- 1
  tau_key <- if (inclusion_prior_strategy == "functional_pi") fmt_num(tau_value) else "NA"

  var_y <- stats::var(data_bundle$y)
  sigma_scalar <- parse_sigma_0_2_scalar(run_row[["sigma_0_2_scalar"]], var_y, L)
  sigma_scalar_fixed <- if (prior_variance_strategy == "fixed") sigma_scalar %||% 0.2 else sigma_scalar
  sigma_key <- if (identical(prior_variance_strategy, "fixed")) fmt_num(sigma_scalar_fixed) else "eb"

  blocked_norm <- normalize_blocked_idx(blocked_idx, p)
  blocked_key <- if (length(blocked_norm)) paste(blocked_norm, collapse = ",") else "none"

  paste(
    paste0("bundle=", fmt_int(run_row$dataset_bundle_id)),
    paste0("use_case=", fmt_chr(run_row$use_case_id)),
    paste0("backend=", fmt_chr(backend)),
    paste0("pm=", fmt_chr(prior_mean_strategy)),
    paste0("pv=", fmt_chr(prior_variance_strategy)),
    paste0("pi=", fmt_chr(inclusion_prior_strategy)),
    paste0("unmap=", fmt_chr(unmappable_effects)),
    paste0("L=", fmt_int(L)),
    paste0("init=", fmt_chr(init_type)),
    paste0("warm_seed=", fmt_int(warm_seed)),
    paste0("alpha_conc=", if (identical(init_type, "warm")) fmt_num(alpha_conc) else "NA"),
    paste0("ann_r2=", ann_r2_key),
    paste0("inflate=", inflate_key),
    paste0("ann_seed=", ann_seed_key),
    paste0("c=", c_key),
    paste0("tau=", tau_key),
    paste0("sigma=", sigma_key),
    paste0("blocked=", blocked_key),
    sep = "|"
  )
}

#' Execute a single run row: simulate data, fit model, and persist outputs.
#' @keywords internal
execute_single_run <- function(run_row, job_config, quiet = FALSE, buffer_ctx = NULL) {
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
  fit <- model_result$fits[[1]]
  meta <- model_result$fit_meta[[1]]

  eval_res <- evaluate_model(
    fit = fit,
    X = data_bundle$X,
    y = data_bundle$y,
    causal_idx = data_bundle$causal_idx,
    rho = job_config$job$credible_set_rho %||% 0.95,
    purity_threshold = job_config$job$purity_threshold %||% 0.5,
    compute_curves = FALSE
  )

  run_data_bundle <- data_bundle
  if (!is.null(model_result$mu_0)) {
    run_data_bundle$mu_0 <- model_result$mu_0
  }
  if (!is.null(model_result$sigma_0_2)) {
    run_data_bundle$sigma_0_2 <- model_result$sigma_0_2
  }

  write_run_outputs(
    run_row = run_row,
    job_config = job_config,
    evaluation = eval_res,
    data_bundle = run_data_bundle,
    model_result = list(fit = fit, extra = NULL),
    buffer_ctx = buffer_ctx,
    variant_meta = meta
  )

  if (!quiet) {
    summary_msg <- sprintf(
      "run_id=%s use_case=%s phenotype_seed=%s L=%s | power=%.3f size=%.2f purity=%.2f",
      run_row$run_id,
      run_row$use_case_id,
      run_row$phenotype_seed,
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
    phenotype_seed = run_row$phenotype_seed,
    power = eval_res$model_filtered$power,
    mean_size = eval_res$model_filtered$mean_size,
    mean_purity = eval_res$model_filtered$mean_purity
  )
}

#' Generate data bundle for a given run_row based on scenario.
#' @keywords internal
generate_data_for_run <- function(run_row, job_config) {
  matrix_row <- dplyr::filter(
    job_config$tables$data_matrices,
    .data$matrix_id == !!run_row$matrix_id
  )
  if (!nrow(matrix_row)) {
    stop("No matrix metadata found (matrix_id=", run_row$matrix_id, ").")
  }
  scenario <- matrix_row$data_scenario[[1]]

  if (scenario %in% c("simulation_n3", "sim_n3")) {
    spec <- list(
      phenotype_seed = as.integer(run_row$phenotype_seed),
      p_star = as.integer(run_row$p_star),
      y_noise = as.numeric(run_row$y_noise),
      annotation_r2 = as.numeric(run_row$annotation_r2),
      inflate_match = as.numeric(run_row$inflate_match),
      effect_sd = run_row[["effect_sd"]] %||% 1,
      architecture = as.character(run_row$architecture %||% "sparse")
    )
    return(generate_simulation_data(spec))
  }
  repo_root <- job_config$paths$repo_root %||% getwd()
  X_mat <- load_sampled_matrix(matrix_row, repo_root = repo_root)
  spec <- list(
    phenotype_seed = as.integer(run_row$phenotype_seed),
    p_star = as.integer(run_row$p_star),
    y_noise = as.numeric(run_row$y_noise),
    annotation_r2 = as.numeric(run_row$annotation_r2),
    inflate_match = as.numeric(run_row$inflate_match),
    effect_sd = run_row[["effect_sd"]] %||% 1,
    architecture = as.character(run_row$architecture %||% "sparse")
  )
  generate_simulation_data(spec, base_X = X_mat)
}

#' Generate data bundle for a dataset bundle row (shared across model runs).
#' @keywords internal
generate_data_for_bundle <- function(bundle_row, job_config) {
  bundle_row <- tibble::as_tibble(bundle_row)
  matrix_row <- dplyr::filter(
    job_config$tables$data_matrices,
    .data$matrix_id == !!bundle_row$matrix_id
  )
  if (!nrow(matrix_row)) {
    stop("No matrix metadata found (matrix_id=", bundle_row$matrix_id, ").")
  }
  scenario <- matrix_row$data_scenario[[1]]
  effect_sd <- bundle_row[["effect_sd"]] %||% 1
  architecture <- as.character(bundle_row$architecture %||% "sparse")
  seed <- as.integer(bundle_row$phenotype_seed %||% 1L)
  p_star <- as.integer(bundle_row$p_star)
  y_noise <- as.numeric(bundle_row$y_noise)

  if (scenario %in% c("simulation_n3", "sim_n3")) {
    data_env <- new.env(parent = emptyenv())
    utils::data("SuSiE_N3_X", package = "test_susine", envir = data_env)
    X <- get("SuSiE_N3_X", envir = data_env)
  } else {
    repo_root <- job_config$paths$repo_root %||% getwd()
    X <- load_sampled_matrix(matrix_row, repo_root = repo_root)
  }

  X <- as.matrix(X)
  effects <- if (identical(architecture, "oligogenic")) {
    simulate_effect_sizes_oligogenic(
      p = ncol(X),
      p_star = p_star,
      effect_sd = effect_sd,
      seed = seed
    )
  } else {
    simulate_effect_sizes(
      p = ncol(X),
      p_star = p_star,
      effect_sd = effect_sd,
      seed = seed
    )
  }
  phenotype <- simulate_phenotype(
    X = X,
    beta = effects$beta,
    noise_fraction = y_noise,
    seed = seed + 1L
  )

  list(
    X = X,
    y = phenotype$y,
    beta = effects$beta,
    sigma2 = phenotype$sigma2,
    causal_idx = effects$causal_idx,
    effect_sd = effect_sd,
    data_scenario = scenario,
    matrix_id = bundle_row$matrix_id,
    architecture = architecture,
    phenotype_seed = seed,
    effect_tier = effects$effect_tier %||% rep("sparse", ncol(X)),
    priors_cache = new.env(parent = emptyenv())
  )
}

#' Execute all model runs for a dataset bundle.
#' @keywords internal
execute_dataset_bundle <- function(bundle_runs, job_config, quiet = FALSE, buffer_ctx = NULL) {
  bundle_id <- unique(bundle_runs$dataset_bundle_id)
  if (length(bundle_id) != 1) {
    stop("execute_dataset_bundle expects a single dataset_bundle_id.")
  }
  bundles_tbl <- job_config$tables$dataset_bundles
  bundle_row <- dplyr::filter(bundles_tbl, dataset_bundle_id == !!bundle_id)
  if (!nrow(bundle_row)) {
    stop("Dataset bundle not found: ", bundle_id)
  }

  data_bundle <- generate_data_for_bundle(bundle_row, job_config)
  bundle_row <- dplyr::mutate(bundle_row, data_scenario = data_bundle$data_scenario)

  # Dataset-level metrics (computed once)
  z_top_k <- job_config$job$metrics$z_top_k %||% 10L
  ds_metrics <- compute_dataset_metrics(data_bundle$X, data_bundle$y, top_k = z_top_k) %>%
    dplyr::mutate(
      dataset_bundle_id = bundle_id,
      data_scenario = data_bundle$data_scenario,
      matrix_id = data_bundle$matrix_id,
      architecture = data_bundle$architecture %||% (bundle_row$architecture %||% "sparse"),
      y_noise = bundle_row$y_noise,
      p_star = bundle_row$p_star,
      phenotype_seed = bundle_row$phenotype_seed
    )
  write_dataset_metrics(bundle_row = bundle_row, job_config = job_config, dataset_metrics = ds_metrics, buffer_ctx = buffer_ctx)

  summary_rows <- list()
  global_primary_pips <- list()
  global_primary_elbos <- c()
  global_run_template <- NULL
  execution_cache <- new.env(parent = emptyenv())
  execution_cache_hits <- 0L
  execution_cache_misses <- 0L

  use_case_ids <- unique(bundle_runs$use_case_id)
  for (uc_id in use_case_ids) {
    use_case <- lookup_use_case(job_config, uc_id)
    uc_runs <- dplyr::filter(bundle_runs, use_case_id == !!uc_id)

    primary_pips_by_group <- list()
    primary_elbos_by_group <- list()
    group_run_map <- list()

    run_and_record <- function(run_row,
                               group_key,
                               blocked_idx = integer(0),
                               meta_override = list()) {
      cache_key <- execution_cache_key(
        use_case = use_case,
        run_row = run_row,
        data_bundle = data_bundle,
        job_config = job_config,
        blocked_idx = blocked_idx
      )

      cache_hit <- exists(cache_key, envir = execution_cache, inherits = FALSE)
      if (cache_hit) {
        cached <- get(cache_key, envir = execution_cache, inherits = FALSE)
        fit <- cached$fit
        eval_res <- cached$evaluation
        run_data_bundle <- cached$run_data_bundle
        meta <- infer_primary_variant_meta(run_row, wall_time_sec = 0)
        meta$model_call_executed <- FALSE
        meta$cache_source_run_id <- as.integer(cached$source_run_id)
        execution_cache_hits <<- execution_cache_hits + 1L
      } else {
        model_result <- run_use_case(
          use_case = use_case,
          run_row = run_row,
          data_bundle = data_bundle,
          job_config = job_config,
          blocked_idx = blocked_idx
        )
        fit <- model_result$fits[[1]]
        meta <- model_result$fit_meta[[1]]
        eval_res <- evaluate_model(
          fit = fit,
          X = data_bundle$X,
          y = data_bundle$y,
          causal_idx = data_bundle$causal_idx,
          rho = job_config$job$credible_set_rho %||% 0.95,
          purity_threshold = job_config$job$purity_threshold %||% 0.5,
          compute_curves = FALSE
        )

        run_data_bundle <- data_bundle
        if (!is.null(model_result$mu_0)) {
          run_data_bundle$mu_0 <- model_result$mu_0
        }
        if (!is.null(model_result$sigma_0_2)) {
          run_data_bundle$sigma_0_2 <- model_result$sigma_0_2
        }
        assign(
          cache_key,
          list(
            fit = fit,
            evaluation = eval_res,
            run_data_bundle = run_data_bundle,
            source_run_id = as.integer(run_row$run_id)
          ),
          envir = execution_cache
        )
        meta$model_call_executed <- TRUE
        meta$cache_source_run_id <- as.integer(run_row$run_id)
        execution_cache_misses <<- execution_cache_misses + 1L
      }
      if (length(meta_override)) {
        meta[names(meta_override)] <- meta_override
      }
      meta$init_type <- as.character(meta$init_type %||% resolve_primary_init_type(run_row))
      meta$model_call_executed <- isTRUE(meta$model_call_executed)
      meta$cache_source_run_id <- as.integer(meta$cache_source_run_id %||% run_row$run_id)
      write_run_outputs(
        run_row = run_row,
        job_config = job_config,
        evaluation = eval_res,
        data_bundle = run_data_bundle,
        model_result = list(fit = fit, extra = NULL),
        buffer_ctx = buffer_ctx,
        variant_meta = meta
      )

      if (!isTRUE(meta$is_agg)) {
        elbo_val <- suppressWarnings(as.numeric(tail(fit$model_fit$elbo, 1)))
        if (!length(elbo_val)) {
          elbo_val <- NA_real_
        }
        primary_pips_by_group[[group_key]] <<- c(primary_pips_by_group[[group_key]], list(fit$model_fit$PIPs))
        primary_elbos_by_group[[group_key]] <<- c(primary_elbos_by_group[[group_key]], elbo_val)
        global_primary_pips[[length(global_primary_pips) + 1L]] <<- fit$model_fit$PIPs
        global_primary_elbos <<- c(global_primary_elbos, elbo_val)
        if (is.null(global_run_template)) {
          global_run_template <<- run_row
        }
      }

      summary_rows[[length(summary_rows) + 1L]] <<- tibble::tibble(
        run_id = run_row$run_id,
        dataset_bundle_id = bundle_id,
        use_case_id = run_row$use_case_id,
        phenotype_seed = run_row$phenotype_seed,
        variant_type = meta$variant_type,
        variant_id = meta$variant_id,
        model_call_executed = as.logical(meta$model_call_executed %||% NA),
        cache_source_run_id = as.integer(meta$cache_source_run_id %||% NA_integer_),
        power = eval_res$model_filtered$power,
        mean_size = eval_res$model_filtered$mean_size,
        mean_purity = eval_res$model_filtered$mean_purity
      )

      list(fit = fit, eval_res = eval_res, meta = meta)
    }

    uc_runs <- uc_runs %>%
      dplyr::mutate(
        .group_key = dplyr::coalesce(
          .data$group_key,
          paste0(
            "L=", .data$L,
            "|r2=", dplyr::coalesce(as.character(.data$annotation_r2), "NA"),
            "|inflate=", dplyr::coalesce(as.character(.data$inflate_match), "NA")
          )
        )
      )
    grouped_idx <- split(seq_len(nrow(uc_runs)), uc_runs$.group_key)
    for (gk in names(grouped_idx)) {
      group_rows <- uc_runs[grouped_idx[[gk]], , drop = FALSE]
      if (is.null(group_run_map[[gk]])) {
        group_run_map[[gk]] <- group_rows[1, , drop = FALSE]
      }

      method_raw <- as.character(group_rows$exploration_methods[[1]] %||% "")
      method_ids <- unique(strsplit(method_raw, "x", fixed = TRUE)[[1]])
      is_refine_group <- "refine" %in% method_ids

      if (!is_refine_group) {
        for (ii in seq_len(nrow(group_rows))) {
          run_and_record(
            run_row = group_rows[ii, , drop = FALSE],
            group_key = gk
          )
        }
        next
      }

      group_rows <- group_rows %>%
        dplyr::arrange(dplyr::coalesce(.data$refine_step, .data$run_id), .data$run_id)

      queue <- list(list(blocked = integer(0)))
      seen_block_sets <- new.env(parent = emptyenv())
      assign(".root", TRUE, envir = seen_block_sets)
      processed <- 0L

      for (ii in seq_len(nrow(group_rows))) {
        if (!length(queue)) {
          break
        }
        node <- queue[[1]]
        if (length(queue) == 1L) {
          queue <- list()
        } else {
          queue <- queue[-1]
        }
        run_row <- group_rows[ii, , drop = FALSE]
        step_id <- as.integer(run_row$refine_step %||% ii)
        out <- run_and_record(
          run_row = run_row,
          group_key = gk,
          blocked_idx = node$blocked,
          meta_override = list(
            variant_type = "refine",
            variant_id = step_id,
            init_type = ifelse(step_id <= 1L, "default", "warm"),
            is_restart = FALSE
          )
        )
        processed <- processed + 1L

        cs_sets <- out$eval_res$effects_filtered$indices
        if (is.null(cs_sets) || !length(cs_sets)) {
          cs_sets <- out$eval_res$effects_unfiltered$indices
        }
        if (is.null(cs_sets) || !length(cs_sets)) {
          next
        }
        for (cs in cs_sets) {
          cs <- as.integer(cs)
          cs <- cs[is.finite(cs)]
          if (!length(cs)) {
            next
          }
          child_blocked <- sort(unique(c(node$blocked, cs)))
          if (length(child_blocked) >= ncol(data_bundle$X)) {
            next
          }
          key <- paste(child_blocked, collapse = ",")
          if (!exists(key, envir = seen_block_sets, inherits = FALSE)) {
            assign(key, TRUE, envir = seen_block_sets)
            queue[[length(queue) + 1L]] <- list(blocked = child_blocked)
          }
        }
      }

      if (processed < nrow(group_rows)) {
        warning(
          "Refine BFS queue exhausted before reaching requested K for group '",
          gk, "' (processed ", processed, " of ", nrow(group_rows), ")."
        )
      }
    }

    # Aggregation across fits (per use case)
    agg_methods <- job_config$job$compute$aggregation_methods %||% character(0)
    softmax_temperature <- as.numeric(job_config$job$compute$softmax_temperature %||% 1)
    group_keys <- names(primary_pips_by_group)
    if (length(group_keys) > 0 && length(agg_methods)) {
      for (group_key in group_keys) {
        group_pips <- primary_pips_by_group[[group_key]]
        group_elbos <- primary_elbos_by_group[[group_key]]
        if (length(group_pips) > 1) {
          agg_results <- aggregate_use_case_pips(
            pip_list = group_pips,
            elbo_vec = group_elbos,
            methods = agg_methods,
            softmax_temperature = softmax_temperature,
            jsd_threshold = job_config$job$metrics$jsd_threshold %||% 0.171661
          )
          for (m in names(agg_results)) {
            agg_pip <- agg_results[[m]]
            agg_meta <- list(
              variant_type = "aggregation",
              variant_id = m,
              agg_method = m,
              agg_ess = as.numeric(attr(agg_pip, "ess") %||% NA_real_),
              is_restart = FALSE,
              is_agg = TRUE
            )
            conf_bins <- compute_confusion_bins(
              agg_pip,
              as.integer(seq_along(data_bundle$beta) %in% data_bundle$causal_idx),
              pip_bucket_width = job_config$job$metrics$pip_bucket_width %||% 0.01
            )
            agg_run_row <- group_run_map[[group_key]] %||% uc_runs[1, , drop = FALSE]
            agg_run_row$restart_id <- NA_integer_
            agg_run_row$init_type <- NA_character_
            agg_run_row$c_value <- NA_real_
            agg_run_row$tau_value <- NA_real_
            agg_run_row$sigma_0_2_scalar <- NA_character_
            write_confusion_bins(
              run_row = agg_run_row,
              job_config = job_config,
              conf_bins = conf_bins,
              buffer_ctx = buffer_ctx,
              variant_meta = agg_meta
            )
          }
        }
      }
    }

    for (group_key in names(primary_pips_by_group)) {
      group_pips <- primary_pips_by_group[[group_key]]
      if (length(group_pips) > 1) {
        mm <- compute_multimodal_metrics(group_pips, jsd_threshold = job_config$job$metrics$jsd_threshold %||% 0.171661)
        write_multimodal_metrics(
          bundle_row = bundle_row,
          use_case_id = uc_id,
          group_label = paste0("model_grid|", group_key),
          metrics = mm,
          job_config = job_config,
          buffer_ctx = buffer_ctx
        )
      }
    }

  }

  # Overall per-dataset aggregation across all fits/use-cases.
  include_overall <- isTRUE(job_config$job$compute$include_overall_pool %||% FALSE)
  overall_methods <- job_config$job$compute$overall_aggregation_methods %||% character(0)
  if (include_overall && length(global_primary_pips) > 1 && length(overall_methods) > 0) {
    softmax_temperature <- as.numeric(job_config$job$compute$softmax_temperature %||% 1)
    overall_agg <- aggregate_use_case_pips(
      pip_list = global_primary_pips,
      elbo_vec = global_primary_elbos,
      methods = overall_methods,
      softmax_temperature = softmax_temperature,
      jsd_threshold = job_config$job$metrics$jsd_threshold %||% 0.171661
    )
    global_row <- global_run_template %||% bundle_runs[1, , drop = FALSE]
    global_row$use_case_id <- "global_pool"
    global_row$prior_spec_id <- "global_pool"
    global_row$restart_id <- NA_integer_
    global_row$init_type <- NA_character_
    global_row$c_value <- NA_real_
    global_row$tau_value <- NA_real_
    global_row$sigma_0_2_scalar <- NA_character_
    for (m in names(overall_agg)) {
      agg_pip <- overall_agg[[m]]
      conf_bins <- compute_confusion_bins(
        agg_pip,
        as.integer(seq_along(data_bundle$beta) %in% data_bundle$causal_idx),
        pip_bucket_width = job_config$job$metrics$pip_bucket_width %||% 0.01
      )
      write_confusion_bins(
        run_row = global_row,
        job_config = job_config,
        conf_bins = conf_bins,
        buffer_ctx = buffer_ctx,
        variant_meta = list(
          variant_type = "overall_aggregation",
          variant_id = m,
          agg_method = m,
          agg_ess = as.numeric(attr(agg_pip, "ess") %||% NA_real_),
          is_restart = FALSE,
          is_agg = TRUE
        )
      )
    }

    mm_global <- compute_multimodal_metrics(
      global_primary_pips,
      jsd_threshold = job_config$job$metrics$jsd_threshold %||% 0.171661
    )
    write_multimodal_metrics(
      bundle_row = bundle_row,
      use_case_id = "global_pool",
      group_label = "overall_pool",
      metrics = mm_global,
      job_config = job_config,
      buffer_ctx = buffer_ctx
    )
  }

  if (!quiet && execution_cache_hits > 0L) {
    message(
      sprintf(
        "Execution cache reused %d fit(s) in bundle %s (computed %d unique fit(s)).",
        execution_cache_hits,
        as.integer(bundle_id),
        execution_cache_misses
      )
    )
  }

  if (length(summary_rows)) {
    return(dplyr::bind_rows(summary_rows))
  }
  tibble::tibble()
}

#' Parse sigma_0_2_scalar string and compute the actual value.
#'
#' Supported formats:
#' - Numeric (e.g., 0.1, 0.2, 0.4): interpreted as multiplier of var_y
#' - String with /L suffix (e.g., "1/L", "0.5/L", "2/L"): interpreted as multiplier/L of var_y
#' - NA or NULL: returns NULL (use susine default)
#'
#' @param scalar_spec The sigma_0_2_scalar value from run_row (character or numeric)
#' @param var_y Variance of y
#' @param L Number of effects
#' @return Scalar value to pass to susine, or NULL for default behavior
#' @keywords internal
parse_sigma_0_2_scalar <- function(scalar_spec, var_y, L) {
  if (is.null(scalar_spec) || (length(scalar_spec) == 1 && is.na(scalar_spec))) {
    return(NULL)
  }
  
  spec_str <- as.character(scalar_spec)
  
  # Check for /L suffix
  if (grepl("/L$", spec_str, ignore.case = TRUE)) {
    # Extract the multiplier before /L
    multiplier <- as.numeric(sub("/L$", "", spec_str, ignore.case = TRUE))
    if (is.na(multiplier)) {
      warning("Could not parse sigma_0_2_scalar: ", scalar_spec, ". Using default.")
      return(NULL)
    }
    return(multiplier / L)
  }
  
  # Plain numeric multiplier
  multiplier <- as.numeric(spec_str)
  if (is.na(multiplier)) {
    warning("Could not parse sigma_0_2_scalar: ", scalar_spec, ". Using default.")
    return(NULL)
  }
  return(multiplier)
}

dirichlet_matrix <- function(n, alpha_vec) {
  k <- length(alpha_vec)
  if (k < 1) {
    stop("alpha_vec must have length >= 1")
  }
  draws <- matrix(stats::rgamma(n * k, shape = alpha_vec, rate = 1), nrow = n, byrow = TRUE)
  draws / rowSums(draws)
}

normalize_susier_fit <- function(fit_raw, X, L) {
  alpha <- as.matrix(fit_raw$alpha)
  if (!nrow(alpha)) {
    stop("susieR fit has empty alpha matrix.")
  }
  if (is.null(L) || !is.finite(L)) {
    L <- nrow(alpha)
  }
  pip <- fit_raw$pip
  if (is.null(pip)) {
    pip <- 1 - apply(1 - alpha, 2, prod)
  }
  mu <- fit_raw$mu
  beta <- rep(0, ncol(X))
  if (!is.null(mu)) {
    beta <- colSums(alpha * as.matrix(mu))
  }
  sigma2 <- fit_raw$sigma2 %||% fit_raw$residual_variance %||% NA_real_
  list(
    settings = list(L = L),
    effect_fits = list(alpha = alpha),
    model_fit = list(
      PIPs = as.numeric(pip),
      elbo = fit_raw$elbo %||% NA_real_,
      sigma_2 = sigma2,
      coef = as.numeric(beta),
      std_coef = as.numeric(beta),
      fitted_y = as.vector(X %*% beta)
    ),
    raw = fit_raw
  )
}

#' Fit one run row under a prior spec (susine or susieR backend).
#' @keywords internal
run_use_case <- function(use_case, run_row, data_bundle, job_config, blocked_idx = integer(0)) {
  L <- as.integer(run_row$L)
  backend <- as.character(use_case$backend[[1]])
  prior_mean_strategy <- as.character(use_case$prior_mean_strategy[[1]])
  prior_variance_strategy <- as.character(use_case$prior_variance_strategy[[1]])
  inclusion_prior_strategy <- as.character(use_case$inclusion_prior_strategy[[1]])
  unmappable_effects <- as.character(use_case$unmappable_effects[[1]] %||% "none")

  annotation_r2 <- as.numeric(run_row$annotation_r2 %||% NA_real_)
  inflate_match <- as.numeric(run_row$inflate_match %||% NA_real_)
  c_value <- as.numeric(run_row$c_value %||% NA_real_)
  tau_value <- as.numeric(run_row$tau_value %||% NA_real_)
  if (!is.finite(c_value)) c_value <- 1
  if (!is.finite(tau_value) || tau_value <= 0) tau_value <- 1

  p <- ncol(data_bundle$X)
  var_y <- stats::var(data_bundle$y)
  restart_cfg <- job_config$job$compute$restart %||% list()
  alpha_conc <- restart_cfg$alpha_concentration %||% 1

  annotation_vec <- NULL
  if (prior_mean_strategy == "functional_mu" || inclusion_prior_strategy == "functional_pi") {
    annotation_seed <- run_row$annotation_seed %||% NA_integer_
    priors <- get_priors_cached(
      cache_env = data_bundle$priors_cache,
      beta = data_bundle$beta,
      annotation_r2 = annotation_r2,
      inflate_match = inflate_match,
      base_sigma2 = stats::var(data_bundle$y),
      effect_sd = data_bundle$effect_sd,
      annotation_seed = annotation_seed
    )
    annotation_vec <- priors$mu_0
  }

  base_prior_weights <- rep(1 / p, p)
  if (inclusion_prior_strategy == "functional_pi") {
    if (is.null(annotation_vec)) {
      stop("Functional pi strategy requires annotation vector.")
    }
    scaled <- abs(as.numeric(annotation_vec)) / tau_value
    scaled <- scaled - max(scaled, na.rm = TRUE)
    pw <- exp(scaled)
    pw <- pw / sum(pw)
    base_prior_weights <- as.numeric(pw)
  }

  restart_id <- suppressWarnings(as.integer(run_row$restart_id %||% NA_integer_))
  refine_step <- suppressWarnings(as.integer(run_row$refine_step %||% NA_integer_))
  init_type <- resolve_primary_init_type(run_row)

  prior_weights <- base_prior_weights
  if (identical(init_type, "warm")) {
    restart_seed <- resolve_restart_seed(run_row, restart_id)
    set.seed(restart_seed)
    prior_weights <- as.numeric(dirichlet_matrix(1L, rep(alpha_conc, p)))
  }
  blocked_idx <- normalize_blocked_idx(blocked_idx, p)
  if (length(blocked_idx)) {
    prior_weights[blocked_idx] <- 0
    if (sum(prior_weights) <= 0) {
      unblocked <- setdiff(seq_len(p), blocked_idx)
      if (!length(unblocked)) {
        prior_weights <- rep(1 / p, p)
      } else {
        prior_weights <- numeric(p)
        prior_weights[unblocked] <- 1 / length(unblocked)
      }
    } else {
      prior_weights <- prior_weights / sum(prior_weights)
    }
  }

  mu_0 <- 0
  if (prior_mean_strategy == "functional_mu") {
    if (is.null(annotation_vec)) {
      stop("Functional mu_0 strategy requires annotation vector.")
    }
    mu_0 <- as.numeric(c_value) * as.numeric(annotation_vec)
  }

  sigma_scalar <- parse_sigma_0_2_scalar(run_row[["sigma_0_2_scalar"]], var_y, L)
  sigma_scalar_fixed <- if (prior_variance_strategy == "fixed") sigma_scalar %||% 0.2 else sigma_scalar

  fit <- NULL
  wall_time_sec <- NA_real_
  started <- proc.time()[["elapsed"]]
  if (backend == "susine") {
    args <- list(
      L = L,
      X = data_bundle$X,
      y = data_bundle$y,
      mu_0 = mu_0,
      prior_inclusion_weights = prior_weights,
      verbose = FALSE
    )
    if (prior_variance_strategy == "fixed") {
      args$sigma_0_2 <- as.numeric(sigma_scalar_fixed) * var_y
      args$prior_update_method <- "none"
    } else {
      args$prior_update_method <- "var"
    }
    fit <- do.call(susine::susine, args)
    fit$settings$L <- fit$settings$L %||% L
  } else if (backend == "susieR") {
    if (!requireNamespace("susieR", quietly = TRUE)) {
      stop("susieR backend requested but package 'susieR' is not installed.")
    }
    susie_formals <- names(formals(susieR::susie))
    args <- list(
      X = data_bundle$X,
      y = data_bundle$y,
      L = L,
      prior_weights = prior_weights,
      estimate_prior_variance = identical(prior_variance_strategy, "eb"),
      estimate_residual_variance = TRUE,
      verbose = FALSE
    )
    if (prior_variance_strategy == "fixed") {
      args$scaled_prior_variance <- as.numeric(sigma_scalar_fixed)
    }
    if ("convergence_method" %in% susie_formals) {
      args$convergence_method <- "elbo"
    }
    if ("unmappable_effects" %in% susie_formals) {
      args$unmappable_effects <- unmappable_effects %||% "none"
    } else if (!is.null(unmappable_effects) && unmappable_effects %in% c("ash", "inf")) {
      stop("Current susieR installation does not support `unmappable_effects`; install susieR >= 2.0.")
    }
    args <- args[intersect(names(args), susie_formals)]
    fit_raw <- do.call(susieR::susie, args)
    fit <- normalize_susier_fit(fit_raw, X = data_bundle$X, L = L)
  } else {
    stop("Unknown backend: ", backend)
  }
  wall_time_sec <- proc.time()[["elapsed"]] - started

  sigma_0_2_truth <- NULL
  if (prior_variance_strategy == "fixed") {
    if (backend == "susine") {
      sigma_0_2_truth <- rep(as.numeric(sigma_scalar_fixed) * var_y, p)
    } else {
      sigma_0_2_truth <- rep(as.numeric(sigma_scalar_fixed), p)
    }
  }

  fits <- list(fit)
  fit_meta <- list(infer_primary_variant_meta(run_row, wall_time_sec = wall_time_sec))

  list(
    fits = fits,
    fit_meta = fit_meta,
    mu_0 = if (length(mu_0) == 1L) rep(mu_0, p) else mu_0,
    sigma_0_2 = sigma_0_2_truth
  )
}

softmax_weights <- function(x, temperature = 1) {
  x <- as.numeric(x)
  temperature <- as.numeric(temperature %||% 1)
  if (!is.finite(temperature) || temperature <= 0) {
    temperature <- 1
  }
  x <- x / temperature
  x <- x - max(x, na.rm = TRUE)
  w <- exp(x)
  w / sum(w)
}

aggregate_use_case_pips <- function(pip_list,
                                    elbo_vec,
                                    methods = c("elbo_softmax"),
                                    softmax_temperature = 1,
                                    jsd_threshold = 0.171661) {
  methods <- unique(as.character(methods))
  methods[methods == "softmax_elbo"] <- "elbo_softmax"
  methods[methods == "mean"] <- "uniform"
  res <- list()
  if (!length(pip_list)) return(res)
  pips_mat <- do.call(rbind, lapply(pip_list, as.numeric))
  if ("uniform" %in% methods) {
    res$uniform <- colMeans(pips_mat)
  }
  if ("max_elbo" %in% methods && length(elbo_vec)) {
    idx <- which.max(elbo_vec)
    res$max_elbo <- pips_mat[idx, ]
  }
  if ("elbo_softmax" %in% methods && length(elbo_vec)) {
    w <- softmax_weights(elbo_vec, temperature = softmax_temperature)
    res$elbo_softmax <- as.numeric(crossprod(w, pips_mat))
  }
  if ("cluster_weight" %in% methods && length(elbo_vec) > 1 && nrow(pips_mat) > 1) {
    n <- nrow(pips_mat)
    jsd_mat <- matrix(0, n, n)
    for (i in 1:(n - 1)) {
      for (j in (i + 1):n) {
        d <- js_distance(pips_mat[i, ], pips_mat[j, ])
        jsd_mat[i, j] <- d
        jsd_mat[j, i] <- d
      }
    }
    hc <- stats::hclust(stats::as.dist(jsd_mat), method = "complete")
    clusters <- stats::cutree(hc, h = jsd_threshold)
    cl_levels <- sort(unique(clusters))
    rep_idx <- integer(length(cl_levels))
    freq <- numeric(length(cl_levels))
    elbo_rep <- numeric(length(cl_levels))
    for (k in seq_along(cl_levels)) {
      cid <- cl_levels[[k]]
      idxs <- which(clusters == cid)
      freq[k] <- length(idxs) / n
      pick <- idxs[which.max(elbo_vec[idxs])]
      rep_idx[k] <- pick
      elbo_rep[k] <- elbo_vec[pick]
    }
    w_rep <- softmax_weights(elbo_rep, temperature = softmax_temperature) / freq
    w_rep <- w_rep / sum(w_rep)
    ens <- as.numeric(crossprod(w_rep, pips_mat[rep_idx, , drop = FALSE]))
    attr(ens, "ess") <- as.numeric(1 / sum(w_rep^2))
    res$cluster_weight <- ens
  }
  res
}

js_distance <- function(p, q, eps = 1e-12) {
  p <- as.numeric(p); q <- as.numeric(q)
  p <- p + eps; q <- q + eps
  m <- 0.5 * (p + q)
  kl <- function(a, b) sum(a * log(a / b))
  0.5 * kl(p, m) + 0.5 * kl(q, m)
}

compute_multimodal_metrics <- function(pip_list, jsd_threshold = 0.171661, top_k = 10L) {
  n <- length(pip_list)
  if (n < 2) {
    return(tibble::tibble(
      mean_jsd = NA_real_,
      median_jsd = NA_real_,
      max_jsd = NA_real_,
      jaccard_top10 = NA_real_,
      mean_pip_var = NA_real_,
      n_clusters = NA_integer_
    ))
  }
  pips_mat <- do.call(rbind, lapply(pip_list, as.numeric))
  jsd_vals <- c()
  for (i in 1:(n - 1)) {
    for (j in (i + 1):n) {
      jsd_vals <- c(jsd_vals, js_distance(pips_mat[i, ], pips_mat[j, ]))
    }
  }
  # Jaccard top-k (mean across pairs)
  k <- as.integer(top_k)
  top_sets <- lapply(seq_len(n), function(i) {
    ord <- order(pips_mat[i, ], decreasing = TRUE)
    ord[seq_len(min(k, length(ord)))]
  })
  pair_jaccard <- c()
  for (i in 1:(n - 1)) {
    for (j in (i + 1):n) {
      inter <- length(intersect(top_sets[[i]], top_sets[[j]]))
      union <- length(union(top_sets[[i]], top_sets[[j]]))
      pair_jaccard <- c(pair_jaccard, if (union > 0) inter / union else NA_real_)
    }
  }
  # PIP variance
  pip_var <- mean(apply(pips_mat, 2, stats::var, na.rm = TRUE))
  # Cluster count by threshold (complete-linkage on JSD)
  jsd_mat <- matrix(0, n, n)
  idx <- 1L
  for (i in 1:(n - 1)) {
    for (j in (i + 1):n) {
      jsd_mat[i, j] <- jsd_vals[idx]
      jsd_mat[j, i] <- jsd_vals[idx]
      idx <- idx + 1L
    }
  }
  hc <- stats::hclust(stats::as.dist(jsd_mat), method = "complete")
  n_clusters <- length(unique(stats::cutree(hc, h = jsd_threshold)))

  tibble::tibble(
    mean_jsd = mean(jsd_vals),
    median_jsd = stats::median(jsd_vals),
    max_jsd = max(jsd_vals),
    jaccard_top10 = mean(pair_jaccard, na.rm = TRUE),
    mean_pip_var = pip_var,
    n_clusters = n_clusters
  )
}

write_dataset_metrics <- function(bundle_row, job_config, dataset_metrics, buffer_ctx = NULL) {
  if (is.null(dataset_metrics) || !nrow(dataset_metrics)) return(invisible(NULL))
  verbose_output <- isTRUE(job_config$job$verbose_file_output)
  if (verbose_output) {
    base_output <- determine_base_output(job_config)
    run_dir <- file.path(base_output, "dataset_metrics")
    ensure_dir(run_dir)
    out_path <- file.path(run_dir, sprintf("dataset_metrics_bundle_%05d.csv", bundle_row$dataset_bundle_id))
    readr::write_csv(dataset_metrics, out_path)
  } else {
    buffer_ctx$buffers$dataset_metrics[[length(buffer_ctx$buffers$dataset_metrics) + 1L]] <- dataset_metrics
  }
}

write_multimodal_metrics <- function(bundle_row, use_case_id, group_label, metrics, job_config, buffer_ctx = NULL) {
  if (is.null(metrics) || !nrow(metrics)) return(invisible(NULL))
  out <- metrics %>%
    dplyr::mutate(
      dataset_bundle_id = bundle_row$dataset_bundle_id,
      data_scenario = bundle_row$data_scenario,
      matrix_id = bundle_row$matrix_id,
      y_noise = bundle_row$y_noise,
      p_star = bundle_row$p_star,
      phenotype_seed = bundle_row$phenotype_seed,
      use_case_id = use_case_id,
      group_label = group_label
    )
  verbose_output <- isTRUE(job_config$job$verbose_file_output)
  if (verbose_output) {
    base_output <- determine_base_output(job_config)
    run_dir <- file.path(base_output, "multimodal_metrics")
    ensure_dir(run_dir)
    out_path <- file.path(run_dir, sprintf("multimodal_metrics_bundle_%05d.csv", bundle_row$dataset_bundle_id))
    readr::write_csv(out, out_path)
  } else {
    buffer_ctx$buffers$multimodal[[length(buffer_ctx$buffers$multimodal) + 1L]] <- out
  }
}

write_confusion_bins <- function(run_row, job_config, conf_bins, buffer_ctx = NULL, variant_meta = NULL) {
  if (is.null(conf_bins) || !nrow(conf_bins)) return(invisible(NULL))
  is_agg <- isTRUE(variant_meta$is_agg %||% FALSE)
  meta <- tibble::tibble(
    run_id = if (is_agg) NA_integer_ else run_row$run_id,
    dataset_bundle_id = run_row$dataset_bundle_id,
    architecture = as.character(run_row$architecture %||% NA_character_),
    use_case_id = run_row$use_case_id,
    prior_spec_id = run_row$prior_spec_id %||% run_row$use_case_id,
    group_key = run_row$group_key %||% NA_character_,
    L = run_row$L,
    annotation_r2 = run_row$annotation_r2,
    inflate_match = run_row$inflate_match,
    sigma_0_2_scalar = run_row$sigma_0_2_scalar %||% NA_character_,
    refine_step = as.integer(run_row$refine_step %||% NA_integer_),
    c_value = run_row$c_value %||% NA_real_,
    tau_value = run_row$tau_value %||% NA_real_,
    init_type = as.character(variant_meta$init_type %||% run_row$init_type %||% NA_character_),
    variant_type = as.character(variant_meta$variant_type %||% "base"),
    variant_id = as.character(variant_meta$variant_id %||% NA_character_),
    agg_method = as.character(variant_meta$agg_method %||% NA_character_),
    agg_ess = as.numeric(variant_meta$agg_ess %||% NA_real_)
  )
  out <- dplyr::bind_cols(meta, conf_bins)
  verbose_output <- isTRUE(job_config$job$verbose_file_output)
  if (verbose_output) {
    base_output <- determine_base_output(job_config)
    run_dir <- file.path(base_output, "confusion_bins")
    ensure_dir(run_dir)
    readr::write_csv(out, file.path(run_dir, sprintf("confusion_bins_run_%05d.csv", run_row$run_id)))
  } else {
    buffer_ctx$buffers$confusion[[length(buffer_ctx$buffers$confusion) + 1L]] <- out
  }
}

#' Persist metrics, truth tables, and model fits for a run.
#' @keywords internal
#' Persist metrics, truth tables, and model fits for a run.
#' @keywords internal
write_run_outputs <- function(run_row,
                              job_config,
                              evaluation,
                              data_bundle,
                              model_result,
                              buffer_ctx = NULL,
                              variant_meta = NULL) {

  verbose_output <- isTRUE(job_config$job$verbose_file_output)
  base_output <- determine_base_output(job_config)
  run_id_val <- as.integer(run_row$run_id)
  run_dir <- file.path(base_output, sprintf("run-%05d", run_id_val))
  if (verbose_output) {
    ensure_dir(run_dir)
  }

  if (is.null(variant_meta)) {
    variant_meta <- list(variant_type = "base", variant_id = NA, agg_method = NA)
  }
  meta_tbl <- tibble::tibble(
    dataset_bundle_id = run_row$dataset_bundle_id,
    architecture = as.character(run_row$architecture %||% NA_character_),
    prior_spec_id = run_row$prior_spec_id %||% run_row$use_case_id,
    refine_step = as.integer(run_row$refine_step %||% NA_integer_),
    init_type = as.character(variant_meta$init_type %||% run_row$init_type %||% NA_character_),
    variant_type = as.character(variant_meta$variant_type %||% "base"),
    variant_id = as.character(variant_meta$variant_id %||% NA_character_),
    agg_method = as.character(variant_meta$agg_method %||% NA_character_),
    wall_time_sec = as.numeric(variant_meta$wall_time_sec %||% NA_real_),
    model_call_executed = as.logical(variant_meta$model_call_executed %||% NA),
    cache_source_run_id = as.integer(variant_meta$cache_source_run_id %||% NA_integer_)
  )

  elbo_vals <- evaluation$traces$elbo
  elbo_final <- NA_real_
  if (!is.null(elbo_vals) && length(elbo_vals)) {
    elbo_non_na <- elbo_vals[!is.na(elbo_vals)]
    if (length(elbo_non_na)) {
      elbo_final <- tail(elbo_non_na, 1)
    }
  }
  model_metrics <- dplyr::mutate(
    dplyr::bind_rows(
      dplyr::mutate(evaluation$model_unfiltered, filtering = "unfiltered"),
      dplyr::mutate(evaluation$model_filtered, filtering = "purity_filtered")
    ),
    run_id = run_row$run_id,
    task_id = run_row$task_id,
    use_case_id = run_row$use_case_id,
    phenotype_seed = run_row$phenotype_seed,
    elbo_final = elbo_final
  ) %>% dplyr::bind_cols(meta_tbl)
  if (verbose_output) {
    readr::write_csv(model_metrics, file.path(run_dir, "model_metrics.csv"))
  }

  effect_metrics <- NULL
  if (!is.null(evaluation$effects_filtered) && nrow(evaluation$effects_filtered)) {
    effect_metrics <- dplyr::mutate(
      dplyr::select(evaluation$effects_filtered, -dplyr::any_of("indices")),
      run_id = run_row$run_id,
      task_id = run_row$task_id,
      use_case_id = run_row$use_case_id,
      phenotype_seed = run_row$phenotype_seed
    ) %>% dplyr::bind_cols(meta_tbl)
  }
  if (!is.null(effect_metrics) && nrow(effect_metrics) && verbose_output) {
    readr::write_csv(effect_metrics, file.path(run_dir, "effect_metrics.csv"))
  }

  snp_tbl <- NULL
  should_write_snps <- isTRUE(job_config$job$write_snps_parquet)
  should_write_legacy <- verbose_output && should_write_legacy_snp_csv(job_config)
  if (should_write_snps || should_write_legacy) {
    snp_tbl <- build_snp_table(
      run_row = run_row,
      data_bundle = data_bundle,
      evaluation = evaluation
    ) %>% dplyr::bind_cols(meta_tbl)
    if (verbose_output && should_write_snps) {
      write_snps_parquet(run_dir, snp_tbl)
    }
  }

  conf_bins <- NULL
  if (isTRUE(job_config$job$write_confusion_bins)) {
    conf_bins <- compute_confusion_bins(
      evaluation$combined_pip,
      as.integer(seq_along(data_bundle$beta) %in% data_bundle$causal_idx),
      pip_bucket_width = job_config$job$metrics$pip_bucket_width %||% 0.01
    )
    write_confusion_bins(
      run_row = run_row,
      job_config = job_config,
      conf_bins = conf_bins,
      buffer_ctx = buffer_ctx,
      variant_meta = variant_meta
    )
  }

  if (should_write_legacy) {
    pip_tbl <- tibble::tibble(
      snp_index = seq_along(evaluation$combined_pip),
      pip = evaluation$combined_pip,
      run_id = run_row$run_id,
      task_id = run_row$task_id,
      use_case_id = run_row$use_case_id,
      phenotype_seed = run_row$phenotype_seed
    )
    readr::write_csv(pip_tbl, file.path(run_dir, "pip.csv"))

    truth_tbl <- tibble::tibble(
      snp_index = seq_along(data_bundle$beta),
      beta = data_bundle$beta,
      mu_0 = data_bundle$mu_0,
      sigma_0_2 = data_bundle$sigma_0_2,
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
    buffer_task_outputs(
      buffer_ctx = buffer_ctx,
      model_metrics = model_metrics,
      effect_metrics = effect_metrics,
      snp_tbl = if (isTRUE(job_config$job$write_snps_parquet)) snp_tbl else NULL,
      validation_row = build_validation_row(run_row)
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

  tibble::tibble(
    run_id = run_row$run_id,
    snp_index = seq_len(n),
    pip = pip,
    beta = beta,
    mu_0 = safe_vec(data_bundle$mu_0),
    sigma_0_2 = safe_vec(data_bundle$sigma_0_2),
    causal = as.integer(seq_len(n) %in% data_bundle$causal_idx)
  )
}

write_snps_parquet <- function(run_dir, snp_tbl) {
  snp_path <- file.path(run_dir, "snps.parquet")
  write_compact_snp_parquet(snp_tbl, snp_path)
}

determine_base_output <- function(job_config) {
  env_parent_id <- Sys.getenv("SUSINE_PARENT_ID", unset = "")
  env_job_name <- Sys.getenv("SUSINE_JOB_NAME", unset = "")
  if (nzchar(env_parent_id) && nzchar(env_job_name)) {
    file.path(
      job_config$paths$slurm_output_dir,
      env_job_name,
      env_parent_id
    )
  } else {
    file.path(job_config$paths$slurm_output_dir, "local_test")
  }
}

build_validation_row <- function(run_row) {
  tibble::tibble(
    run_id = run_row$run_id,
    task_id = run_row$task_id,
    has_issues = FALSE,
    issues = NA_character_
  )
}

buffer_task_outputs <- function(buffer_ctx,
                                model_metrics,
                                effect_metrics,
                                snp_tbl,
                                validation_row,
                                confusion_bins = NULL,
                                dataset_metrics = NULL,
                                multimodal_metrics = NULL) {
  if (is.null(buffer_ctx)) {
    stop("Buffer context is missing for streamed outputs.")
  }
  if (!is.null(model_metrics) && nrow(model_metrics)) {
    buffer_ctx$buffers$model[[length(buffer_ctx$buffers$model) + 1L]] <- model_metrics
  }
  if (!is.null(effect_metrics) && nrow(effect_metrics)) {
    buffer_ctx$buffers$effect[[length(buffer_ctx$buffers$effect) + 1L]] <- effect_metrics
  }
  if (!is.null(snp_tbl) && nrow(snp_tbl)) {
    buffer_ctx$buffers$snps[[length(buffer_ctx$buffers$snps) + 1L]] <- snp_tbl
  }
  if (!is.null(confusion_bins) && nrow(confusion_bins)) {
    buffer_ctx$buffers$confusion[[length(buffer_ctx$buffers$confusion) + 1L]] <- confusion_bins
  }
  if (!is.null(dataset_metrics) && nrow(dataset_metrics)) {
    buffer_ctx$buffers$dataset_metrics[[length(buffer_ctx$buffers$dataset_metrics) + 1L]] <- dataset_metrics
  }
  if (!is.null(multimodal_metrics) && nrow(multimodal_metrics)) {
    buffer_ctx$buffers$multimodal[[length(buffer_ctx$buffers$multimodal) + 1L]] <- multimodal_metrics
  }
  if (!is.null(validation_row) && nrow(validation_row)) {
    buffer_ctx$buffers$validation[[length(buffer_ctx$buffers$validation) + 1L]] <- validation_row
  }
}

flush_task_buffers <- function(buffer_ctx) {
  buf <- buffer_ctx$buffers
  if (is.null(buf) ||
      !length(buf$model) && !length(buf$effect) && !length(buf$snps) &&
      !length(buf$confusion) && !length(buf$dataset_metrics) && !length(buf$multimodal) &&
      !length(buf$validation)) {
    return(invisible(NULL))
  }
  base_output <- buffer_ctx$base_output
  if (is.null(base_output)) {
    warning("Buffered outputs found but base output directory is missing.")
    return(invisible(NULL))
  }
  staging_dir <- file.path(
    base_output,
    sprintf("task-%03d", as.integer(buffer_ctx$task_id))
  )
  ensure_dir(staging_dir)
  flush_label <- sprintf("flush-%03d", as.integer(buffer_ctx$flush_index))

    write_flush_outputs(
      staging_dir = staging_dir,
      flush_label = flush_label,
      task_id = buffer_ctx$task_id,
      model_metrics = if (length(buf$model)) dplyr::bind_rows(buf$model) else NULL,
      effect_metrics = if (length(buf$effect)) dplyr::bind_rows(buf$effect) else NULL,
      snp_tbl = if (length(buf$snps)) dplyr::bind_rows(buf$snps) else NULL,
      validation_row = if (length(buf$validation)) dplyr::bind_rows(buf$validation) else NULL,
      confusion_bins = if (length(buf$confusion)) dplyr::bind_rows(buf$confusion) else NULL,
      dataset_metrics = if (length(buf$dataset_metrics)) dplyr::bind_rows(buf$dataset_metrics) else NULL,
      multimodal_metrics = if (length(buf$multimodal)) dplyr::bind_rows(buf$multimodal) else NULL
    )

  buffer_ctx$buffers <- list(
    model = list(),
    effect = list(),
    snps = list(),
    confusion = list(),
    dataset_metrics = list(),
    multimodal = list(),
    validation = list()
  )
  buffer_ctx$flush_index <- buffer_ctx$flush_index + 1L
}

write_flush_outputs <- function(staging_dir,
                                flush_label,
                                task_id,
                                model_metrics,
                                effect_metrics,
                                snp_tbl,
                                validation_row,
                                confusion_bins = NULL,
                                dataset_metrics = NULL,
                                multimodal_metrics = NULL) {
  if (!is.null(model_metrics) && nrow(model_metrics)) {
    model_metrics <- dplyr::mutate(model_metrics, task_id = task_id, flush_id = flush_label)
    readr::write_csv(model_metrics, file.path(staging_dir, sprintf("%s_model_metrics.csv", flush_label)))
  }
  if (!is.null(effect_metrics) && nrow(effect_metrics)) {
    effect_metrics <- dplyr::mutate(effect_metrics, task_id = task_id, flush_id = flush_label)
    readr::write_csv(effect_metrics, file.path(staging_dir, sprintf("%s_effect_metrics.csv", flush_label)))
  }
  if (!is.null(snp_tbl) && nrow(snp_tbl)) {
    snp_tbl <- dplyr::mutate(snp_tbl, task_id = task_id, flush_id = flush_label)
    write_compact_snp_parquet(snp_tbl, file.path(staging_dir, sprintf("%s_snps.parquet", flush_label)))
  }
  if (!is.null(validation_row) && nrow(validation_row)) {
    validation_row <- dplyr::mutate(validation_row, task_id = task_id, flush_id = flush_label)
    readr::write_csv(validation_row, file.path(staging_dir, sprintf("%s_validation.csv", flush_label)))
  }
  if (!is.null(confusion_bins) && nrow(confusion_bins)) {
    confusion_bins <- dplyr::mutate(confusion_bins, task_id = task_id, flush_id = flush_label)
    readr::write_csv(confusion_bins, file.path(staging_dir, sprintf("%s_confusion_bins.csv", flush_label)))
  }
  if (!is.null(dataset_metrics) && nrow(dataset_metrics)) {
    dataset_metrics <- dplyr::mutate(dataset_metrics, task_id = task_id, flush_id = flush_label)
    readr::write_csv(dataset_metrics, file.path(staging_dir, sprintf("%s_dataset_metrics.csv", flush_label)))
  }
  if (!is.null(multimodal_metrics) && nrow(multimodal_metrics)) {
    multimodal_metrics <- dplyr::mutate(multimodal_metrics, task_id = task_id, flush_id = flush_label)
    readr::write_csv(multimodal_metrics, file.path(staging_dir, sprintf("%s_multimodal_metrics.csv", flush_label)))
  }
}

# Write SNP parquet with reduced precision for numeric columns to cut size.
write_compact_snp_parquet <- function(snp_tbl, path) {
  # Ensure task/flush metadata exist even if caller did not add them.
  if (!"task_id" %in% names(snp_tbl)) {
    snp_tbl$task_id <- NA_integer_
  }
  if (!"flush_id" %in% names(snp_tbl)) {
    snp_tbl$flush_id <- NA_character_
  }
  snp_tbl <- dplyr::mutate(
    snp_tbl,
    run_id = as.integer(run_id),
    task_id = as.integer(task_id),
    dataset_bundle_id = as.integer(dataset_bundle_id %||% NA_integer_),
    architecture = as.character(architecture %||% NA_character_),
    prior_spec_id = as.character(prior_spec_id %||% NA_character_),
    refine_step = as.integer(refine_step %||% NA_integer_),
    init_type = as.character(init_type %||% NA_character_),
    snp_index = as.integer(snp_index),
    pip = as.numeric(pip),
    beta = as.numeric(beta),
    mu_0 = as.numeric(mu_0),
    sigma_0_2 = as.numeric(sigma_0_2),
    causal = as.integer(causal),
    flush_id = as.character(flush_id),
    variant_type = as.character(variant_type %||% NA_character_),
    variant_id = as.character(variant_id %||% NA_character_),
    agg_method = as.character(agg_method %||% NA_character_)
  )
  # Fixed schema that includes task/flush metadata and keeps floats compact.
  snp_schema <- arrow::schema(
    run_id = arrow::int32(),
    task_id = arrow::int32(),
    dataset_bundle_id = arrow::int32(),
    architecture = arrow::utf8(),
    prior_spec_id = arrow::utf8(),
    refine_step = arrow::int32(),
    init_type = arrow::utf8(),
    snp_index = arrow::int32(),
    pip = arrow::float32(),
    beta = arrow::float32(),
    mu_0 = arrow::float32(),
    sigma_0_2 = arrow::float32(),
    causal = arrow::int8(),
    flush_id = arrow::utf8(),
    variant_type = arrow::utf8(),
    variant_id = arrow::utf8(),
    agg_method = arrow::utf8()
  )
  # Align columns to schema order; drop any extras to avoid schema mismatch.
  snp_tbl <- snp_tbl[, names(snp_schema), drop = FALSE]
  snp_table <- arrow::Table$create(snp_tbl, schema = snp_schema)
  arrow::write_parquet(snp_table, path, compression = "zstd")
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



