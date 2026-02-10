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

  if (!is.null(cfg$tables$dataset_bundles)) {
    if (!"phenotype_seed" %in% names(cfg$tables$dataset_bundles) && "seed" %in% names(cfg$tables$dataset_bundles)) {
      cfg$tables$dataset_bundles <- dplyr::rename(cfg$tables$dataset_bundles, phenotype_seed = .data$seed)
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
  dplyr::filter(job_config$tables$use_cases, use_case_id == !!use_case_id)
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

prior_cache_key <- function(annotation_r2, inflate_match, gamma_shrink, annotation_seed = NA_integer_) {
  sprintf("r2=%s|inflate=%s|gamma=%s|seed=%s",
          format(annotation_r2, digits = 6, scientific = FALSE),
          format(inflate_match, digits = 6, scientific = FALSE),
          format(gamma_shrink, digits = 6, scientific = FALSE),
          format(annotation_seed, digits = 10, scientific = FALSE))
}

get_priors_cached <- function(cache_env,
                              beta,
                              annotation_r2,
                              inflate_match,
                              gamma_shrink,
                              base_sigma2,
                              effect_sd,
                              annotation_seed = NA_integer_) {
  key <- prior_cache_key(annotation_r2, inflate_match, gamma_shrink, annotation_seed)
  if (exists(key, envir = cache_env, inherits = FALSE)) {
    return(get(key, envir = cache_env, inherits = FALSE))
  }
  priors <- simulate_priors(
    beta = beta,
    annotation_r2 = annotation_r2,
    inflate_match = inflate_match,
    gamma_shrink = gamma_shrink,
    base_sigma2 = base_sigma2,
    effect_sd = effect_sd,
    seed = annotation_seed
  )
  assign(key, priors, envir = cache_env)
  priors
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

  eval_res <- evaluate_model(
    fit = model_result$fit,
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
    model_result = model_result,
    buffer_ctx = buffer_ctx
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
      gamma_shrink = as.numeric(run_row$gamma_shrink),
      effect_sd = run_row[["effect_sd"]] %||% 1
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
    gamma_shrink = as.numeric(run_row$gamma_shrink),
    effect_sd = run_row[["effect_sd"]] %||% 1
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
  effects <- simulate_effect_sizes(
    p = ncol(X),
    p_star = p_star,
    effect_sd = effect_sd,
    seed = seed
  )
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
    phenotype_seed = seed,
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
      y_noise = bundle_row$y_noise,
      p_star = bundle_row$p_star,
      phenotype_seed = bundle_row$phenotype_seed
    )
  write_dataset_metrics(bundle_row = bundle_row, job_config = job_config, dataset_metrics = ds_metrics, buffer_ctx = buffer_ctx)

  summary_rows <- list()

  use_case_ids <- unique(bundle_runs$use_case_id)
  for (uc_id in use_case_ids) {
    use_case <- lookup_use_case(job_config, uc_id)
    uc_runs <- dplyr::filter(bundle_runs, use_case_id == !!uc_id)

    primary_pips_by_group <- list()
    primary_elbos_by_group <- list()
    group_run_map <- list()

    for (i in seq_len(nrow(uc_runs))) {
      run_row <- uc_runs[i, , drop = FALSE]
      group_key <- run_row$group_key %||% paste0(
        "L=", run_row$L,
        "|r2=", run_row$annotation_r2 %||% "NA",
        "|inflate=", run_row$inflate_match %||% "NA"
      )
      if (is.null(group_run_map[[group_key]])) {
        group_run_map[[group_key]] <- run_row
      }
      model_result <- run_use_case(
        use_case = use_case,
        run_row = run_row,
        data_bundle = data_bundle,
        job_config = job_config
      )

      # Base fit always present
      fit_list <- model_result$fits
      fit_meta <- model_result$fit_meta

      for (k in seq_along(fit_list)) {
        fit <- fit_list[[k]]
        meta <- fit_meta[[k]]
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

        if (!isTRUE(meta$is_agg)) {
          primary_pips_by_group[[group_key]] <- c(primary_pips_by_group[[group_key]], list(fit$model_fit$PIPs))
          primary_elbos_by_group[[group_key]] <- c(primary_elbos_by_group[[group_key]], tail(fit$model_fit$elbo, 1))
        }

        summary_rows[[length(summary_rows) + 1L]] <- tibble::tibble(
          run_id = run_row$run_id,
          dataset_bundle_id = bundle_id,
          use_case_id = run_row$use_case_id,
          phenotype_seed = run_row$phenotype_seed,
          variant_type = meta$variant_type,
          variant_id = meta$variant_id,
          power = eval_res$model_filtered$power,
          mean_size = eval_res$model_filtered$mean_size,
          mean_purity = eval_res$model_filtered$mean_purity
        )

      }
    }

    # Aggregation across fits (per use case)
    agg_methods <- job_config$job$compute$aggregation_methods %||% character(0)
    group_keys <- names(primary_pips_by_group)
    if (length(group_keys) > 0 && length(agg_methods)) {
      for (group_key in group_keys) {
        group_pips <- primary_pips_by_group[[group_key]]
        group_elbos <- primary_elbos_by_group[[group_key]]
        if (length(group_pips) > 1) {
          agg_results <- aggregate_use_case_pips(group_pips, group_elbos, methods = agg_methods)
          for (m in names(agg_results)) {
            agg_pip <- agg_results[[m]]
            agg_meta <- list(
              variant_type = "aggregation",
              variant_id = m,
              agg_method = m,
              is_restart = FALSE,
              is_agg = TRUE
            )
            conf_bins <- compute_confusion_bins(
              agg_pip,
              as.integer(seq_along(data_bundle$beta) %in% data_bundle$causal_idx),
              pip_bucket_width = job_config$job$metrics$pip_bucket_width %||% 0.01
            )
            agg_run_row <- group_run_map[[group_key]] %||% uc_runs[1, , drop = FALSE]
            agg_run_row$gamma_shrink <- NA_real_
            agg_run_row$annotation_scale <- NA_real_
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

#' Fit a use case with the susine backend, handling optional annealing/model averaging.
#' @keywords internal
run_use_case <- function(use_case, run_row, data_bundle, job_config) {
  L <- as.integer(run_row$L)
  mu_strategy <- use_case$mu_strategy[[1]]
  sigma_strategy <- use_case$sigma_strategy[[1]]

  annotation_scale <- NA_real_
  if ("annotation_scale" %in% names(run_row)) {
    annotation_scale <- suppressWarnings(as.numeric(run_row$annotation_scale))
  }

  annotation_r2 <- as.numeric(run_row$annotation_r2 %||% NA_real_)
  inflate_match <- as.numeric(run_row$inflate_match %||% NA_real_)
  gamma_shrink <- as.numeric(run_row$gamma_shrink %||% NA_real_)

  uses_mu_ann <- mu_strategy %in% c("functional", "eb_mu")
  uses_sigma_ann <- sigma_strategy %in% c("functional", "eb_sigma")

  mu_0 <- 0
  sigma_0_2 <- NULL

  if (uses_mu_ann || uses_sigma_ann) {
    annotation_seed <- run_row$annotation_seed %||% NA_integer_
    priors <- get_priors_cached(
      cache_env = data_bundle$priors_cache,
      beta = data_bundle$beta,
      annotation_r2 = annotation_r2,
      inflate_match = inflate_match,
      gamma_shrink = gamma_shrink,
      base_sigma2 = stats::var(data_bundle$y),
      effect_sd = data_bundle$effect_sd,
      annotation_seed = annotation_seed
    )
    if (uses_mu_ann) {
      mu_0 <- priors$mu_0
      if (mu_strategy == "functional") {
        scale_val <- annotation_scale
        if (is.null(scale_val) || length(scale_val) == 0 || is.na(scale_val) || !is.finite(scale_val)) {
          scale_val <- 1
        }
        mu_0 <- mu_0 * scale_val
      }
    }
    if (uses_sigma_ann) {
      sigma_0_2 <- priors$sigma_0_2
    }
  }

  if (sigma_strategy == "naive" && !is.null(run_row[["sigma_0_2_scalar"]])) {
    var_y <- var(data_bundle$y)
    sigma_0_2 <- parse_sigma_0_2_scalar(run_row[["sigma_0_2_scalar"]], var_y, L)
  } else if (!uses_sigma_ann) {
    sigma_0_2 <- NULL
  }

  extra <- use_case$extra_compute[[1]]
  if (length(extra) && is.na(extra)) {
    extra <- NULL
  }
  anneal <- job_config$job$compute$anneal
  model_avg <- job_config$job$compute$model_average
  restart_cfg <- job_config$job$compute$restart

  base_args <- list(
    L = L,
    X = data_bundle$X,
    y = data_bundle$y,
    mu_0 = mu_0,
    sigma_0_2 = sigma_0_2,
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

  fits <- list()
  fit_meta <- list()

  if (identical(extra, "restart_alpha")) {
    alpha_conc <- restart_cfg$alpha_concentration %||% 1
    p <- ncol(data_bundle$X)
    alpha_vec <- rep(alpha_conc, p)
    restart_id <- suppressWarnings(as.integer(run_row$restart_id %||% NA_integer_))
    restart_seed <- suppressWarnings(as.integer(run_row$restart_seed %||% NA_integer_))

    if (!is.na(restart_id)) {
      if (!is.na(restart_seed)) {
        set.seed(restart_seed)
      } else {
        set.seed(as.integer(run_row$phenotype_seed) + restart_id)
      }
      init_alpha <- dirichlet_matrix(L, alpha_vec)
      base_args$init_alpha <- init_alpha
      fit <- do.call(susine::susine, base_args)
      fit$settings$L <- fit$settings$L %||% L
      fits[[1]] <- fit
      fit_meta[[1]] <- list(
        variant_type = "restart",
        variant_id = restart_id,
        is_restart = TRUE,
        is_agg = FALSE
      )
      return(list(fits = fits, fit_meta = fit_meta, mu_0 = mu_0, sigma_0_2 = sigma_0_2))
    }

    # Backward-compatible fallback (single run row drives multiple restarts)
    n_inits <- restart_cfg$n_inits %||% 5
    for (i in seq_len(n_inits)) {
      set.seed(as.integer(run_row$phenotype_seed) + i)
      init_alpha <- dirichlet_matrix(L, alpha_vec)
      base_args$init_alpha <- init_alpha
      fit <- do.call(susine::susine, base_args)
      fit$settings$L <- fit$settings$L %||% L
      fits[[length(fits) + 1L]] <- fit
      fit_meta[[length(fit_meta) + 1L]] <- list(
        variant_type = "restart",
        variant_id = i,
        is_restart = TRUE,
        is_agg = FALSE
      )
    }
    return(list(fits = fits, fit_meta = fit_meta, mu_0 = mu_0, sigma_0_2 = sigma_0_2))
  }

  if (!identical(extra, "model_avg")) {
    fit <- do.call(susine::susine, base_args)
    fit$settings$L <- fit$settings$L %||% L
    fits[[1]] <- fit
    fit_meta[[1]] <- list(
      variant_type = "base",
      variant_id = NA_integer_,
      is_restart = FALSE,
      is_agg = FALSE
    )
    return(list(fits = fits, fit_meta = fit_meta, mu_0 = mu_0, sigma_0_2 = sigma_0_2))
  }

  n_inits <- model_avg$n_inits %||% 5
  init_sd <- model_avg$init_sd %||% 0.05
  subfits <- vector("list", n_inits)
  pips_mat <- matrix(0, nrow = n_inits, ncol = length(data_bundle$beta))
  coef_mat <- matrix(0, nrow = n_inits, ncol = length(data_bundle$beta))

  for (i in seq_len(n_inits)) {
    base_args$init_random <- TRUE
    base_args$init_seed <- as.integer(run_row$phenotype_seed) + i
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

  fits[[1]] <- agg_fit
  fit_meta[[1]] <- list(
    variant_type = "model_avg",
    variant_id = NA_integer_,
    is_restart = FALSE,
    is_agg = FALSE
  )

  list(fits = fits, fit_meta = fit_meta, mu_0 = mu_0, sigma_0_2 = sigma_0_2)
}

softmax_weights <- function(x) {
  x <- as.numeric(x)
  x <- x - max(x, na.rm = TRUE)
  w <- exp(x)
  w / sum(w)
}

aggregate_use_case_pips <- function(pip_list, elbo_vec, methods = c("softmax_elbo")) {
  methods <- unique(methods)
  res <- list()
  if (!length(pip_list)) return(res)
  pips_mat <- do.call(rbind, lapply(pip_list, as.numeric))
  if ("mean" %in% methods) {
    res$mean <- colMeans(pips_mat)
  }
  if ("max_elbo" %in% methods && length(elbo_vec)) {
    idx <- which.max(elbo_vec)
    res$max_elbo <- pips_mat[idx, ]
  }
  if ("softmax_elbo" %in% methods && length(elbo_vec)) {
    w <- softmax_weights(elbo_vec)
    res$softmax_elbo <- as.numeric(crossprod(w, pips_mat))
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
  hc <- stats::hclust(stats::as.dist(jsd_mat), method = "average")
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
    use_case_id = run_row$use_case_id,
    group_key = run_row$group_key %||% NA_character_,
    L = run_row$L,
    annotation_r2 = run_row$annotation_r2,
    inflate_match = run_row$inflate_match,
    gamma_shrink = run_row$gamma_shrink,
    annotation_scale = run_row$annotation_scale %||% NA_real_,
    sigma_0_2_scalar = run_row$sigma_0_2_scalar %||% NA_character_,
    variant_type = as.character(variant_meta$variant_type %||% "base"),
    variant_id = as.character(variant_meta$variant_id %||% NA_character_),
    agg_method = as.character(variant_meta$agg_method %||% NA_character_)
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
    variant_type = as.character(variant_meta$variant_type %||% "base"),
    variant_id = as.character(variant_meta$variant_id %||% NA_character_),
    agg_method = as.character(variant_meta$agg_method %||% NA_character_)
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

