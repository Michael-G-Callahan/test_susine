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
      tier_cs_metrics = list(),
      snps = list(),
      confusion = list(),
      dataset_metrics = list(),
      multimodal = list(),
      refine_depth = list(),
      scaling_bins = list(),
      validation = list(),
      prior_diagnostics = list(),
      hg2_by_agg = list()
    )
    buffer_ctx$base_output <- determine_base_output(job_config)
    buffer_ctx$task_id <- task_id
    buffer_ctx$flush_index <- 0L
    # Reruns should start from a clean task staging directory to avoid
    # stale/duplicate flush files and overwrite-lock surprises.
    prepare_task_staging_dir(buffer_ctx$base_output, task_id)
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

  if (!verbose_output && !is.null(buffer_ctx)) {
    staging_dir <- file.path(
      buffer_ctx$base_output,
      sprintf("task-%03d", as.integer(task_id))
    )
    if (!dir.exists(staging_dir)) {
      warning("Expected staging directory not found after task flush: ", staging_dir)
    } else {
      staged_files <- list.files(staging_dir, pattern = "^flush-[0-9]+_", full.names = FALSE)
      if (!length(staged_files)) {
        warning("No flush files found in staging directory after task flush: ", staging_dir)
      } else {
        type_labels <- sub("^flush-[0-9]+_", "", staged_files)
        type_labels <- sub("\\.(csv|parquet)$", "", type_labels)
        type_counts <- table(type_labels)
        if (!quiet) {
          message(
            sprintf(
              "Task %s staged file counts: %s",
              task_id,
              paste0(names(type_counts), "=", as.integer(type_counts), collapse = ", ")
            )
          )
        }
        has_validation <- any(grepl("_validation\\.csv$", staged_files))
        if (!has_validation) {
          warning(
            sprintf(
              "No validation flush CSV detected for task %s in %s. This will prevent aggregated validation.csv from being produced.",
              task_id,
              staging_dir
            )
          )
        }
      }
    }
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
compute_confusion_bins <- function(pip, causal, pip_bucket_width = 0.01,
                                   causal_mask = NULL,
                                   pip_breaks = NULL) {
  n <- length(pip)
  if (length(causal) != n) {
    stop("pip and causal lengths differ.")
  }
  if (!is.null(causal_mask)) {
    causal <- as.integer(seq_len(n) %in% causal_mask)
  }

  # Build bin breaks: use pip_breaks if provided, else uniform grid
  if (is.null(pip_breaks)) {
    pip_breaks <- seq(0, 1, by = pip_bucket_width)
  }
  # findInterval: bin index 1 = [break[1], break[2]), ..., rightmost closed
  bin_idx <- findInterval(pip, pip_breaks, rightmost.closed = TRUE, left.open = FALSE)
  # Map each variant to the lower edge of its bin
  bins <- pip_breaks[pmax(1L, bin_idx)]

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

#' Build variable-width PIP bin breaks optimized for TPR/FPR resolution.
#' Fine bins where FPR action occurs (low PIP), coarser where it doesn't.
#' @return Numeric vector of break points from 0 to 1.
#' @export
pip_breaks_variable <- function() {
  c(
    seq(0.000, 0.0275, by = 0.0025),   # 12 bins: [0.00, 0.03)
    seq(0.030, 0.095,  by = 0.005),     # 14 bins: [0.03, 0.10)
    seq(0.10,  0.19,   by = 0.01),      # 10 bins: [0.10, 0.20)
    seq(0.20,  0.90,   by = 0.05),      # 15 bins: [0.20, 0.95)
    seq(0.95,  1.00,   by = 0.01)       #  5 bins: [0.95, 1.00]
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

resolve_run_type <- function(run_row) {
  restart_id <- suppressWarnings(as.integer(run_row$restart_id %||% NA_integer_))
  run_type <- as.character(run_row$run_type %||% run_row$init_type %||% NA_character_)
  if (is.na(run_type) || !nzchar(run_type)) {
    run_type <- if (!is.na(restart_id) && restart_id > 1L) "warm" else "default"
  }
  run_type
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
  run_type <- resolve_run_type(run_row)
  is_refine <- !is.na(refine_step)
  list(
    explore_method = if (is_refine) "refine" else if (!is.na(restart_id)) "restart" else "base",
    variant_id = if (is_refine) refine_step else if (!is.na(restart_id)) restart_id else NA_integer_,
    run_type = run_type,
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

  run_type <- resolve_run_type(run_row)
  restart_id <- suppressWarnings(as.integer(run_row$restart_id %||% NA_integer_))
  warm_seed <- if (identical(run_type, "warm")) resolve_restart_seed(run_row, restart_id) else NA_integer_
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
    paste0("init=", fmt_chr(run_type)),
    paste0("warm_seed=", fmt_int(warm_seed)),
    paste0("alpha_conc=", if (identical(run_type, "warm")) fmt_num(alpha_conc) else "NA"),
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
  } else if (identical(architecture, "susie2_sparse")) {
    simulate_effect_sizes_susie2_sparse(
      p = ncol(X),
      p_star = p_star,
      seed = seed
    )
  } else if (identical(architecture, "susie2_oligogenic")) {
    simulate_effect_sizes_susie2_oligogenic(
      p = ncol(X),
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

  # SuSiE 2.0 architectures use h2-calibrated noise; others use noise_fraction
  phenotype <- if (architecture %in% c("susie2_sparse", "susie2_oligogenic")) {
    simulate_phenotype_h2(
      X = X,
      beta = effects$beta,
      h2_snp_per_causal = if (identical(architecture, "susie2_sparse")) {
        as.numeric(bundle_row[["h2_snp_per_causal"]] %||% 0.03)
      } else {
        NA_real_
      },
      h2_total = if (identical(architecture, "susie2_oligogenic")) {
        as.numeric(bundle_row[["h2_total"]] %||% 0.25)
      } else {
        NA_real_
      },
      seed = seed + 1L
    )
  } else {
    simulate_phenotype(
      X = X,
      beta = effects$beta,
      noise_fraction = y_noise,
      seed = seed + 1L
    )
  }

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

  # Top-8 causal mask: score TP/FP only on the 8 largest-effect causal variants
  top8_causal_mask <- if (length(data_bundle$causal_idx) > 0L && !is.null(data_bundle$beta)) {
    n_top <- min(8L, length(data_bundle$causal_idx))
    data_bundle$causal_idx[
      order(abs(data_bundle$beta[data_bundle$causal_idx]), decreasing = TRUE)
    ][seq_len(n_top)]
  } else {
    data_bundle$causal_idx
  }

  # Dataset-level metrics (computed once)
  z_top_k <- job_config$job$metrics$z_top_k %||% 10L
  # Slim dataset_metrics: only dataset_bundle_id + data_scenario + metrics.
  # Bundle-level columns (matrix_id, architecture, y_noise, p_star,
  # phenotype_seed) are stripped — re-attached during aggregation via
  # join on dataset_bundle_id to dataset_bundles.
  ds_metrics <- compute_dataset_metrics(data_bundle$X, data_bundle$y, top_k = z_top_k) %>%
    dplyr::mutate(
      dataset_bundle_id = bundle_id,
      data_scenario = data_bundle$data_scenario
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
    primary_run_meta_by_group <- list()
    fitted_y_by_group <- list()
    group_run_map <- list()

    run_and_record <- function(run_row,
                               group_key,
                               blocked_idx = integer(0),
                               meta_override = list(),
                               init_alpha_override = NULL) {
      cache_key <- execution_cache_key(
        use_case = use_case,
        run_row = run_row,
        data_bundle = data_bundle,
        job_config = job_config,
        blocked_idx = blocked_idx
      )

      # Skip cache for refinement refit steps — each has a unique
      # init_alpha_override that the cache key cannot distinguish.
      cache_hit <- is.null(init_alpha_override) &&
        exists(cache_key, envir = execution_cache, inherits = FALSE)
      if (cache_hit) {
        cached <- get(cache_key, envir = execution_cache, inherits = FALSE)
        fit <- cached$fit
        eval_res <- cached$evaluation
        run_data_bundle <- cached$run_data_bundle
        meta <- infer_primary_variant_meta(
          run_row, wall_time_sec = as.numeric(cached$wall_time_sec %||% 0)
        )
        meta$model_call_executed <- FALSE
        meta$cache_source_run_id <- as.integer(cached$source_run_id)
        execution_cache_hits <<- execution_cache_hits + 1L
      } else {
        model_result <- run_use_case(
          use_case = use_case,
          run_row = run_row,
          data_bundle = data_bundle,
          job_config = job_config,
          blocked_idx = blocked_idx,
          init_alpha_override = init_alpha_override
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
            source_run_id = as.integer(run_row$run_id),
            wall_time_sec = as.numeric(meta$wall_time_sec %||% 0)
          ),
          envir = execution_cache
        )
        meta$model_call_executed <- TRUE
        meta$cache_source_run_id <- as.integer(run_row$run_id)
        execution_cache_misses <<- execution_cache_misses + 1L
      }
      if (length(meta_override)) {
        # Add perturb wall time to refit wall time for refine steps so the
        # recorded wall_time_sec reflects the true cost of the full step.
        if (!is.null(meta_override$wall_time_perturb_sec)) {
          meta$wall_time_sec <- as.numeric(meta$wall_time_sec %||% 0) +
            as.numeric(meta_override$wall_time_perturb_sec)
          meta_override$wall_time_perturb_sec <- NULL
        }
        meta[names(meta_override)] <- meta_override
      }
      meta$run_type <- as.character(meta$run_type %||% resolve_run_type(run_row))
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

      # Tier-stratified CS power for susie2_oligogenic architectures
      if (!is.null(buffer_ctx) && isTRUE(job_config$job$write_tier_cs_metrics) &&
          identical(as.character(run_row$architecture), "susie2_oligogenic") &&
          !isTRUE(meta$is_agg) &&
          !is.null(data_bundle$beta) && length(data_bundle$causal_idx) > 0L) {
        tier_cs <- compute_cs_power_by_tier(
          effects_filtered = eval_res$effects_filtered,
          causal_idx       = data_bundle$causal_idx,
          beta             = data_bundle$beta
        ) %>%
          dplyr::mutate(
            run_id        = as.integer(run_row$run_id),
            filtering     = "purity_filtered",
            explore_method = as.character(meta$explore_method %||% "base"),
            agg_method    = as.character(meta$agg_method %||% NA_character_)
          )
        buffer_ctx$buffers$tier_cs_metrics[[
          length(buffer_ctx$buffers$tier_cs_metrics) + 1L
        ]] <- tier_cs
      }

      # Extract and buffer prior diagnostics for susine fits
      if (!is.null(buffer_ctx) && isTRUE(job_config$job$write_prior_diagnostics)) {
        prior_diag <- extract_prior_diagnostics(
          fit = fit,
          run_id = as.integer(run_row$run_id),
          backend = as.character(use_case$backend[[1]])
        )
        if (!is.null(prior_diag) && nrow(prior_diag)) {
          buffer_ctx$buffers$prior_diagnostics[[
            length(buffer_ctx$buffers$prior_diagnostics) + 1L
          ]] <- prior_diag
        }
      }

      if (!isTRUE(meta$is_agg)) {
        elbo_val <- suppressWarnings(as.numeric(tail(fit$model_fit$elbo, 1)))
        if (!length(elbo_val)) {
          elbo_val <- NA_real_
        }
        primary_pips_by_group[[group_key]] <<- c(primary_pips_by_group[[group_key]], list(fit$model_fit$PIPs))
        primary_elbos_by_group[[group_key]] <<- c(primary_elbos_by_group[[group_key]], elbo_val)
        fy_val <- tryCatch(as.numeric(fit$model_fit$fitted_y), error = function(e) NULL)
        if (is.null(fy_val)) {
          b_val <- tryCatch(
            fit$model_fit$std_coef %||% fit$model_fit$coef, error = function(e) NULL
          )
          if (!is.null(b_val)) fy_val <- as.numeric(data_bundle$X %*% b_val)
        }
        fitted_y_by_group[[group_key]] <<- c(fitted_y_by_group[[group_key]], list(fy_val))
        primary_run_meta_by_group[[group_key]] <<- c(
          primary_run_meta_by_group[[group_key]],
          list(list(
            run_id           = as.integer(run_row$run_id),
            spec_name        = as.character(run_row$spec_name %||% NA_character_),
            restart_id       = as.integer(run_row$restart_id %||% NA_integer_),
            refine_step      = as.integer(run_row$refine_step %||% NA_integer_),
            c_value          = as.numeric(run_row$c_value %||% NA_real_),
            sigma_0_2_scalar = as.numeric(run_row$sigma_0_2_scalar %||% NA_real_),
            annotation_r2    = as.numeric(run_row$annotation_r2 %||% NA_real_)
          ))
        )
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
        explore_method = meta$explore_method,
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
    refine_depth_by_group <- list()
    for (gk in names(grouped_idx)) {
      group_rows <- uc_runs[grouped_idx[[gk]], , drop = FALSE]
      if (is.null(group_run_map[[gk]])) {
        grr <- group_rows[1, , drop = FALSE]
        # Store planned axis sizes for scaling analysis (before BFS may truncate)
        grr$.planned_n_restart <- length(unique(group_rows$restart_id))
        grr$.planned_n_refine  <- length(unique(group_rows$refine_step))
        grr$.planned_n_c       <- length(unique(group_rows$c_value))
        grr$.planned_n_sigma   <- length(unique(group_rows$sigma_0_2_scalar))
        group_run_map[[gk]] <- grr
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

      # Split refine group into independent per-parent BFS trees.
      # Parent = unique combination of non-refine exploration axes
      # (restart_id, c_value, sigma_0_2_scalar, tau_value).
      parent_key_cols <- intersect(
        c("restart_id", "c_value", "sigma_0_2_scalar", "tau_value"),
        names(group_rows)
      )
      parent_key_strs <- do.call(paste, c(
        lapply(parent_key_cols, function(col) {
          paste0(col, "=", dplyr::coalesce(as.character(group_rows[[col]]), "NA"))
        }),
        sep = "|"
      ))
      parent_split <- split(seq_len(nrow(group_rows)), parent_key_strs)

      # Two-step refinement BFS per parent: perturb (blocked weights) -> refit
      # (restored weights + init_alpha from perturb). Root node (step 1) is a
      # baseline fit with no blocking. Each non-root node costs 2 model fits
      # (perturb + refit); only the refit is stored in ensemble PIPs.
      refine_purity_threshold <- as.numeric(
        job_config$job$compute$refine$purity_threshold %||% 0.95
      )

      refine_depth_records <- list()

      for (pk in names(parent_split)) {
        parent_rows <- group_rows[parent_split[[pk]], , drop = FALSE] %>%
          dplyr::arrange(.data$refine_step, .data$run_id)
        n_requested <- nrow(parent_rows)

        queue <- list(list(blocked = integer(0), parent_fit = NULL))
        seen_block_sets <- new.env(parent = emptyenv())
        assign(".root", TRUE, envir = seen_block_sets)
        processed <- 0L

        for (ii in seq_len(n_requested)) {
          if (!length(queue)) {
            break
          }
          node <- queue[[1]]
          if (length(queue) == 1L) {
            queue <- list()
          } else {
            queue <- queue[-1]
          }
          run_row <- parent_rows[ii, , drop = FALSE]
          step_id <- as.integer(run_row$refine_step %||% ii)
          is_root <- (step_id <= 1L && length(node$blocked) == 0L)

          if (is_root) {
            # Root: baseline fit (no blocking, no warm start) -----
            out <- run_and_record(
              run_row = run_row,
              group_key = gk,
              blocked_idx = integer(0),
              meta_override = list(
                explore_method = "refine",
                variant_id = step_id,
                run_type = "default",
                is_restart = FALSE
              )
            )
            processed <- processed + 1L
            stored_fit <- out$fit
            stored_eval <- out$eval_res
          } else {
            # Non-root: two-step perturb + refit -----
            # Step 1 (perturb): fit with blocked prior weights.  This is an
            # intermediate model whose alpha we extract for the warm start.
            # We call run_use_case directly so the perturb is NOT recorded
            # in ensemble PIPs or written to outputs.
            perturb_started <- proc.time()[["elapsed"]]
            perturb_result <- run_use_case(
              use_case = use_case,
              run_row = run_row,
              data_bundle = data_bundle,
              job_config = job_config,
              blocked_idx = node$blocked
            )
            perturb_wall_sec <- proc.time()[["elapsed"]] - perturb_started
            perturb_fit <- perturb_result$fits[[1]]
            perturb_alpha <- perturb_fit$effect_fits$alpha

            # Step 2 (refit): restore uniform prior weights, warm-start
            # from perturbed model's alpha.  The init_alpha_override
            # bypasses normal warm-start logic and injects perturbed
            # alpha directly.  This is stored as the actual ensemble member.
            # wall_time includes both perturb + refit (the true cost of
            # one refine step).
            out <- run_and_record(
              run_row = run_row,
              group_key = gk,
              blocked_idx = integer(0),
              meta_override = list(
                explore_method = "refine",
                variant_id = step_id,
                run_type = "warm",
                is_restart = FALSE,
                wall_time_perturb_sec = perturb_wall_sec
              ),
              init_alpha_override = perturb_alpha
            )
            processed <- processed + 1L
            stored_fit <- out$fit
            stored_eval <- out$eval_res
          }

          # Extract eligible CSs from the stored (refit/root) model for
          # branching.  Only CSs with purity above threshold are eligible.
          cs_sets <- list()
          eff_filt <- stored_eval$effects_filtered
          if (!is.null(eff_filt) && !is.null(eff_filt$indices) &&
              !is.null(eff_filt$purity)) {
            for (ci in seq_along(eff_filt$indices)) {
              pur <- as.numeric(eff_filt$purity[[ci]])
              if (is.finite(pur) && pur >= refine_purity_threshold) {
                cs_sets[[length(cs_sets) + 1L]] <- as.integer(eff_filt$indices[[ci]])
              }
            }
          }
          if (!length(cs_sets)) {
            next
          }
          for (cs in cs_sets) {
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
              queue[[length(queue) + 1L]] <- list(
                blocked = child_blocked,
                parent_fit = stored_fit
              )
            }
          }
        }

        refine_depth_records[[pk]] <- list(
          parent_key         = pk,
          n_refine_requested = n_requested,
          n_refine_processed = processed
        )

        if (processed < n_requested) {
          warning(
            "Refine BFS queue exhausted for parent '", pk, "' in group '", gk,
            "' (processed ", processed, " of ", n_requested, " steps)."
          )
        }
      }

      refine_depth_by_group[[gk]] <- refine_depth_records
    }

    # Aggregation across fits (per use case)
    agg_methods <- job_config$job$compute$aggregation_methods %||% character(0)
    softmax_temperature <- as.numeric(job_config$job$compute$softmax_temperature %||% 1)
    group_keys <- names(primary_pips_by_group)
    need_scaling_bins <- !is.null(buffer_ctx) &&
      isTRUE(job_config$job$compute$write_scaling_confusion_bins)
    full_group_cache <- list()
    if (length(group_keys) > 0 && length(agg_methods)) {
      for (group_key in group_keys) {
        group_pips <- primary_pips_by_group[[group_key]]
        group_elbos <- primary_elbos_by_group[[group_key]]
        if (length(group_pips) > 1) {
          group_cache <- prepare_pip_similarity_cache(
            do.call(rbind, lapply(group_pips, as.numeric))
          )
          pip_mat_cols <- t(group_cache$pips_mat)
          agg_methods_full <- unique(c(
            agg_methods,
            if (need_scaling_bins) "cluster_weight_050" else character(0)
          ))
          agg_results <- aggregate_use_case_pips(
            pip_list = group_pips,
            elbo_vec = group_elbos,
            methods = agg_methods_full,
            softmax_temperature = softmax_temperature,
            jsd_threshold = job_config$job$metrics$jsd_threshold %||% 0.15,
            pip_cache = group_cache
          )
          full_group_bins <- list()
          if (need_scaling_bins) {
            scaling_methods <- c("uniform", "max_elbo", "elbo_softmax",
                                 "cluster_weight", "cluster_weight_050")
            for (sm in scaling_methods) {
              agg_pip_scaling <- aggregate_pip_matrix(
                pip_mat_cols,
                group_elbos,
                sm,
                hc = group_cache$hc %||% NULL
              )
              full_group_bins[[sm]] <- compute_confusion_bins(
                agg_pip_scaling,
                as.integer(seq_along(data_bundle$beta) %in% data_bundle$causal_idx),
                pip_bucket_width = job_config$job$metrics$pip_bucket_width %||% 0.01,
                pip_breaks = job_config$job$metrics$pip_breaks,
                causal_mask = top8_causal_mask
              )
            }
          }
          for (m in names(agg_results)) {
            agg_pip <- agg_results[[m]]
            conf_bins <- compute_confusion_bins(
              agg_pip,
              as.integer(seq_along(data_bundle$beta) %in% data_bundle$causal_idx),
              pip_bucket_width = job_config$job$metrics$pip_bucket_width %||% 0.01,
              pip_breaks = job_config$job$metrics$pip_breaks,
              causal_mask = top8_causal_mask
            )
            if (!m %in% agg_methods) next
            agg_meta <- list(
              explore_method = "aggregation",
              variant_id = m,
              agg_method = m,
              agg_ess = as.numeric(attr(agg_pip, "ess") %||% NA_real_),
              is_restart = FALSE,
              is_agg = TRUE
            )
            agg_run_row <- group_run_map[[group_key]] %||% uc_runs[1, , drop = FALSE]
            agg_run_row$restart_id <- NA_integer_
            agg_run_row$run_type <- NA_character_
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
          full_group_cache[[group_key]] <- list(
            pip_cache = group_cache,
            agg_results = agg_results,
            bins = full_group_bins
          )
        }
      }
    }

    # Scaling confusion bins (subsample-N ensemble analysis stored as bins).
    if (need_scaling_bins) {
      sc_n_sizes     <- as.integer(job_config$job$compute$scaling_n_ens_sizes %||% c(4L, 8L, 16L, 32L, 64L))
      sc_restart_reps <- as.integer(job_config$job$compute$scaling_restart_reps %||% 50L)
      sc_pip_breaks  <- job_config$job$metrics$pip_breaks
      sc_pip_bw      <- as.numeric(job_config$job$metrics$pip_bucket_width %||% 0.01)
      causal_vec_sc  <- as.integer(seq_along(data_bundle$beta) %in% data_bundle$causal_idx)
      for (group_key in names(primary_pips_by_group)) {
        pip_list  <- primary_pips_by_group[[group_key]]
        elbo_vec  <- primary_elbos_by_group[[group_key]]
        meta_list <- primary_run_meta_by_group[[group_key]]
        if (length(pip_list) < 2L || length(meta_list) < 2L) next
        sc_bins <- compute_scaling_confusion_bins_for_group(
          pip_list         = pip_list,
          elbo_vec         = elbo_vec,
          meta_list        = meta_list,
          causal_vec       = causal_vec_sc,
          causal_mask      = top8_causal_mask,
          group_run_row    = group_run_map[[group_key]] %||% uc_runs[1, , drop = FALSE],
          n_ens_sizes      = sc_n_sizes,
          n_restart_reps   = sc_restart_reps,
          pip_breaks       = sc_pip_breaks,
          pip_bw           = sc_pip_bw,
          dataset_bundle_id = bundle_id,
          full_agg_cache    = full_group_cache[[group_key]]
        )
        if (!is.null(sc_bins) && nrow(sc_bins)) {
          buffer_ctx$buffers$scaling_bins[[
            length(buffer_ctx$buffers$scaling_bins) + 1L
          ]] <- sc_bins
        }
      }
    }

    # hg2 by aggregation method (one scalar per agg_method per group, full K fits).
    if (!is.null(buffer_ctx)) {
      for (group_key in names(primary_pips_by_group)) {
        pip_list    <- primary_pips_by_group[[group_key]]
        elbo_vec    <- primary_elbos_by_group[[group_key]]
        fy_list     <- fitted_y_by_group[[group_key]]
        grr         <- group_run_map[[group_key]] %||% uc_runs[1, , drop = FALSE]
        if (length(pip_list) >= 2L && !is.null(fy_list)) {
          hg2_tbl <- compute_hg2_by_agg(
            pip_list          = pip_list,
            fitted_y_list     = fy_list,
            elbo_vec          = elbo_vec,
            y                 = data_bundle$y,
            group_run_row     = grr,
            dataset_bundle_id = bundle_id,
            pip_cache         = full_group_cache[[group_key]]$pip_cache %||% NULL
          )
          if (!is.null(hg2_tbl) && nrow(hg2_tbl) > 0L) {
            buffer_ctx$buffers$hg2_by_agg[[
              length(buffer_ctx$buffers$hg2_by_agg) + 1L
            ]] <- hg2_tbl
          }
        }
      }
    }

    for (group_key in names(primary_pips_by_group)) {
      group_pips <- primary_pips_by_group[[group_key]]
      if (length(group_pips) > 1) {
        mm <- compute_multimodal_metrics(
          group_pips,
          jsd_threshold = job_config$job$metrics$jsd_threshold %||% 0.15,
          pip_cache = full_group_cache[[group_key]]$pip_cache %||% NULL
        )
        # Extract explore method from group_key (e.g., "L=10|...|explore=separate:restart" -> "restart")
        gk_explore <- sub(".*explore=[^:]*:", "", group_key)
        write_multimodal_metrics(
          bundle_row = bundle_row,
          use_case_id = uc_id,
          group_label = paste0("model_grid|", group_key),
          explore_method = gk_explore,
          metrics = mm,
          job_config = job_config,
          buffer_ctx = buffer_ctx
        )
      }
    }

    # Emit per-parent refine depth records for all refine groups.
    for (group_key in names(refine_depth_by_group)) {
      rd_records <- refine_depth_by_group[[group_key]]
      if (!is.null(rd_records) && length(rd_records) > 0L) {
        rd_tbl <- dplyr::bind_rows(lapply(rd_records, function(r) {
          tibble::tibble(
            parent_key         = r$parent_key,
            n_refine_requested = as.integer(r$n_refine_requested),
            n_refine_processed = as.integer(r$n_refine_processed)
          )
        })) %>%
          dplyr::mutate(
            dataset_bundle_id = bundle_row$dataset_bundle_id,
            use_case_id       = uc_id,
            group_label       = paste0("model_grid|", group_key)
          )
        write_refine_depth(rd_tbl, job_config, buffer_ctx)
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
      jsd_threshold = job_config$job$metrics$jsd_threshold %||% 0.15,
      jsd_threshold = job_config$job$metrics$jsd_threshold %||% 0.15
    )
    global_row <- global_run_template %||% bundle_runs[1, , drop = FALSE]
    global_row$use_case_id <- "global_pool"
    global_row$prior_spec_id <- "global_pool"
    global_row$restart_id <- NA_integer_
    global_row$run_type <- NA_character_
    global_row$c_value <- NA_real_
    global_row$tau_value <- NA_real_
    global_row$sigma_0_2_scalar <- NA_character_
    for (m in names(overall_agg)) {
      agg_pip <- overall_agg[[m]]
      conf_bins <- compute_confusion_bins(
        agg_pip,
        as.integer(seq_along(data_bundle$beta) %in% data_bundle$causal_idx),
        pip_bucket_width = job_config$job$metrics$pip_bucket_width %||% 0.01,
        pip_breaks = job_config$job$metrics$pip_breaks,
        causal_mask = top8_causal_mask
      )
      write_confusion_bins(
        run_row = global_row,
        job_config = job_config,
        conf_bins = conf_bins,
        buffer_ctx = buffer_ctx,
        variant_meta = list(
          explore_method = "overall_aggregation",
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
      jsd_threshold = job_config$job$metrics$jsd_threshold %||% 0.15
    )
    write_multimodal_metrics(
      bundle_row = bundle_row,
      use_case_id = "global_pool",
      group_label = "global_pool|all_use_cases",
      explore_method = "all",
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
  # With very small shape parameters, rgamma can return all zeros in a row.
  # Resample any all-zero rows (up to 100 attempts) to avoid 0/0 = NaN.
  rs <- rowSums(draws)
  zero_rows <- which(rs == 0)
  attempts <- 0L
  while (length(zero_rows) > 0 && attempts < 100L) {
    draws[zero_rows, ] <- matrix(
      stats::rgamma(length(zero_rows) * k, shape = alpha_vec, rate = 1),
      nrow = length(zero_rows), byrow = TRUE
    )
    rs <- rowSums(draws)
    zero_rows <- which(rs == 0)
    attempts <- attempts + 1L
  }
  if (length(zero_rows) > 0) {
    # Fallback: uniform for any remaining all-zero rows
    draws[zero_rows, ] <- 1
    rs <- rowSums(draws)
  }
  draws / rs
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

#' Extract per-effect prior diagnostics from a converged fit.
#'
#' Returns a tibble with one row per effect (L rows), recording the final
#' learned prior parameters: scale factor c_l, sigma_0_2_l, and tau2_l.
#' Only meaningful for susine backend fits.
#'
#' @param fit A susine fit object (from `susine::susine()`).
#' @param run_id Integer run identifier.
#' @param backend Character scalar: `"susine"` or `"susieR"`.
#' @return A tibble with columns `run_id`, `effect_l`, `c_l`, `mu_0_l`,
#'   `sigma_0_2_l`, `tau2_l`, or `NULL` if not a susine fit.
#' @keywords internal
extract_prior_diagnostics <- function(fit, run_id, backend) {
  if (backend != "susine") return(NULL)
  priors <- fit$priors
  if (is.null(priors)) return(NULL)

  L <- length(priors$mu_0_scale_factor %||% numeric(0))
  if (L == 0L) return(NULL)

  c_l <- as.numeric(priors$mu_0_scale_factor)
  tau2_l <- as.numeric(priors$mu_0_scale_tau2 %||% rep(NA_real_, L))
  # sigma_0_2 is L x p; take column 1 as representative (all p identical per effect)
  sigma_0_2_l <- if (!is.null(priors$sigma_0_2) && is.matrix(priors$sigma_0_2)) {
    as.numeric(priors$sigma_0_2[, 1])
  } else {
    rep(NA_real_, L)
  }
  # mu_0 is L x p; take column 1 as representative (for "mean"/"mu_var" methods
  # all p values are identical; for annotation methods it's SNP-specific so
  # column 1 is just a sample).
  mu_0_l <- if (!is.null(priors$mu_0) && is.matrix(priors$mu_0)) {
    as.numeric(priors$mu_0[, 1])
  } else {
    rep(NA_real_, L)
  }

  tibble::tibble(
    run_id = as.integer(run_id),
    effect_l = seq_len(L),
    c_l = c_l,
    mu_0_l = mu_0_l,
    sigma_0_2_l = sigma_0_2_l,
    tau2_l = tau2_l
  )
}

#' Fit one run row under a prior spec (susine or susieR backend).
#' @keywords internal
run_use_case <- function(use_case, run_row, data_bundle, job_config,
                         blocked_idx = integer(0),
                         init_alpha_override = NULL) {
  L <- as.integer(run_row$L)
  max_iter <- as.integer(job_config$job$max_iter %||% 100L)
  tol <- as.numeric(job_config$job$tol %||% 1e-3)
  backend <- as.character(use_case$backend[[1]])
  prior_mean_strategy <- as.character(use_case$prior_mean_strategy[[1]])
  prior_variance_strategy <- as.character(use_case$prior_variance_strategy[[1]])
  inclusion_prior_strategy <- as.character(use_case$inclusion_prior_strategy[[1]])
  unmappable_effects <- as.character(use_case$unmappable_effects[[1]] %||% "none")
  eb_method <- use_case$eb_method[[1]] %||% NA_character_
  c_nonneg_flag <- isTRUE(use_case$c_nonneg[[1]])

  annotation_r2 <- as.numeric(run_row$annotation_r2 %||% NA_real_)
  inflate_match <- as.numeric(run_row$inflate_match %||% NA_real_)
  c_value <- as.numeric(run_row$c_value %||% NA_real_)
  tau_value <- as.numeric(run_row$tau_value %||% NA_real_)
  if (!is.finite(c_value)) c_value <- 1
  if (!is.finite(tau_value) || tau_value <= 0) tau_value <- 1

  p <- ncol(data_bundle$X)
  var_y <- stats::var(data_bundle$y)
  restart_cfg <- job_config$job$compute$restart %||% list()
  alpha_conc <- as.numeric(run_row$alpha_concentration %||%
    restart_cfg$alpha_concentration %||% 1)

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
  run_type <- resolve_run_type(run_row)

  prior_weights <- base_prior_weights
  warm_init_alpha <- NULL
  warm_method <- as.character(run_row$warm_method %||% "init_alpha")
  if (identical(run_type, "warm")) {
    restart_seed <- resolve_restart_seed(run_row, restart_id)
    set.seed(restart_seed)
    if (identical(warm_method, "init_alpha")) {
      if (alpha_conc == 0) {
        # One-hot: each effect picks one random SNP
        warm_init_alpha <- matrix(0, nrow = L, ncol = p)
        for (l in seq_len(L)) {
          warm_init_alpha[l, sample.int(p, 1L)] <- 1
        }
      } else {
        warm_init_alpha <- dirichlet_matrix(L, rep(alpha_conc, p))
      }
    } else if (identical(warm_method, "prior_refit")) {
      # Step 1: Dirichlet prior_weights for exploration; step 2 refits after convergence
      prior_weights <- as.numeric(dirichlet_matrix(1L, rep(alpha_conc, p)))
    } else if (identical(warm_method, "truth_warm")) {
      # Warm-start from true causal variants: one effect per causal SNP,
      # each effect concentrates all mass on a single true causal.
      causal_idx <- data_bundle$causal_idx
      if (is.null(causal_idx) || !length(causal_idx)) {
        warning("truth_warm requested but no causal_idx in data_bundle. Falling back to default init.")
      } else {
        # Use top L causal variants by absolute effect size
        beta <- data_bundle$beta
        if (!is.null(beta) && length(beta) == p) {
          causal_order <- causal_idx[order(abs(beta[causal_idx]), decreasing = TRUE)]
        } else {
          causal_order <- causal_idx
        }
        n_warm <- min(L, length(causal_order))
        warm_init_alpha <- matrix(1 / p, nrow = L, ncol = p)
        for (l in seq_len(n_warm)) {
          warm_init_alpha[l, ] <- 0
          warm_init_alpha[l, causal_order[l]] <- 1
        }
      }
    }
  }
  # Override init_alpha if provided directly (e.g., from refinement refit step)
  if (!is.null(init_alpha_override)) {
    warm_init_alpha <- init_alpha_override
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
    if (!is.null(warm_init_alpha)) {
      warm_init_alpha[, blocked_idx] <- 0
      warm_init_alpha <- warm_init_alpha / rowSums(warm_init_alpha)
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
      max_iter = max_iter,
      tol = tol,
      verbose = FALSE
    )
    if (!is.null(warm_init_alpha)) {
      args$init_alpha <- warm_init_alpha
    }
    if (prior_variance_strategy == "fixed") {
      # sigma_0_2 is a proportion of var(y); susine::initialize_priors multiplies
      # by var_y internally, so pass the raw scalar — do NOT pre-multiply.
      args$sigma_0_2 <- as.numeric(sigma_scalar_fixed)
      args$prior_update_method <- "none"
    } else {
      args$prior_update_method <- if (!is.na(eb_method)) eb_method else "var"
      args$c_nonneg <- c_nonneg_flag
      # For "scale" and "mean" EB: sigma is fixed at user-specified value,
      # only c (or scalar mu) is optimized.
      if (!is.na(eb_method) && eb_method %in% c("scale", "mean") &&
          !is.null(sigma_scalar_fixed)) {
        args$sigma_0_2 <- as.numeric(sigma_scalar_fixed)
      }
    }
    # scale_var_outer: inner IBSS does "var" only; outer loop re-estimates c
    if (!is.na(eb_method) && eb_method == "scale_var_outer") {
      n_outer <- 5L
      current_c <- c_value
      for (outer_iter in seq_len(n_outer)) {
        args$mu_0 <- as.numeric(current_c) * as.numeric(annotation_vec)
        fit <- do.call(susine::susine, args)
        # Re-estimate c via neg_log_lik_scale optimization (reuses susine internals)
        new_c <- tryCatch({
          c_lower <- if (c_nonneg_flag) 0 else -100
          sigma_0_2_fixed <- fit$priors$sigma_0_2[1, ]
          opt <- stats::optim(
            par = current_c,
            fn = susine::neg_log_lik_scale,
            annotation = as.numeric(annotation_vec),
            X = data_bundle$X,
            y = data_bundle$y - fit$model_fit$fitted_y + data_bundle$X %*% (fit$priors$mu_0[1, ] * fit$effect_fits$alpha[1, ]),
            sigma_2 = fit$model_fit$sigma_2[length(fit$model_fit$sigma_2)],
            sigma_0_2_fixed = sigma_0_2_fixed,
            prior_inclusion_weights = fit$priors$prior_inclusion_weights[1, ],
            method = "Brent",
            lower = c_lower,
            upper = 100
          )$par
          as.numeric(opt)
        }, error = function(e) current_c)
        if (abs(new_c - current_c) < 1e-4) break
        current_c <- new_c
      }
    } else {
      fit <- do.call(susine::susine, args)
    }
    # prior_refit step 2: refit with uniform pi + converged alpha as init
    if (identical(run_type, "warm") && identical(warm_method, "prior_refit") &&
        !is.null(fit)) {
      refit_alpha <- fit$effect_fits$alpha
      refit_alpha <- refit_alpha / rowSums(refit_alpha)
      args$prior_inclusion_weights <- base_prior_weights
      args$init_alpha <- refit_alpha
      fit <- do.call(susine::susine, args)
    }
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
      max_iter = max_iter,
      tol = tol,
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

prepare_pip_similarity_cache <- function(pips_mat) {
  n <- nrow(pips_mat)
  if (n < 2L) {
    return(list(
      pips_mat = pips_mat,
      jsd_vals = numeric(0),
      jsd_mat = matrix(0, n, n),
      hc = NULL
    ))
  }

  jsd_vals <- numeric(n * (n - 1L) / 2L)
  jsd_mat <- matrix(0, n, n)
  idx <- 1L
  for (i in seq_len(n - 1L)) {
    for (j in seq.int(i + 1L, n)) {
      d <- js_distance(pips_mat[i, ], pips_mat[j, ])
      jsd_vals[[idx]] <- d
      jsd_mat[i, j] <- d
      jsd_mat[j, i] <- d
      idx <- idx + 1L
    }
  }

  list(
    pips_mat = pips_mat,
    jsd_vals = jsd_vals,
    jsd_mat = jsd_mat,
    hc = stats::hclust(stats::as.dist(jsd_mat), method = "complete")
  )
}

aggregate_use_case_pips <- function(pip_list,
                                    elbo_vec,
                                    methods = c("elbo_softmax"),
                                    softmax_temperature = 1,
                                    jsd_threshold = 0.15,
                                    pip_cache = NULL) {
  methods <- unique(as.character(methods))
  methods[methods == "softmax_elbo"] <- "elbo_softmax"
  methods[methods == "mean"] <- "uniform"
  res <- list()
  if (!length(pip_list)) return(res)
  pips_mat <- pip_cache$pips_mat %||% do.call(rbind, lapply(pip_list, as.numeric))
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
  if (any(methods %in% c("cluster_weight", "cluster_weight_050")) &&
      length(elbo_vec) > 1 && nrow(pips_mat) > 1) {
    hc_jsd <- pip_cache$hc %||% prepare_pip_similarity_cache(pips_mat)$hc

    if ("cluster_weight" %in% methods) {
      # cluster_weight: categorical JSD at jsd_threshold (default 0.15) -----
      res$cluster_weight <- .cluster_weight_from_hc(
        hc_jsd, jsd_threshold, elbo_vec, pips_mat, softmax_temperature
      )
    }

    if ("cluster_weight_050" %in% methods) {
      # cluster_weight_050: categorical JSD at 0.50 threshold -----
      res$cluster_weight_050 <- .cluster_weight_from_hc(
        hc_jsd, 0.50, elbo_vec, pips_mat, softmax_temperature
      )
    }

    if ("cluster_weight_jsd_050" %in% methods) {
      # Backward-compatible alias.
      res$cluster_weight_jsd_050 <- .cluster_weight_from_hc(
        hc_jsd, 0.50, elbo_vec, pips_mat, softmax_temperature
      )
    }
  }
  if ("cluster_weight_jsd_050" %in% methods &&
      !("cluster_weight_jsd_050" %in% names(res)) &&
      length(elbo_vec) > 1 && nrow(pips_mat) > 1) {
    hc_jsd <- pip_cache$hc %||% prepare_pip_similarity_cache(pips_mat)$hc
    res$cluster_weight_jsd_050 <- .cluster_weight_from_hc(
      hc_jsd, 0.50, elbo_vec, pips_mat, softmax_temperature
    )
  }
  if ("cluster_weight" %in% methods && !(length(elbo_vec) > 1 && nrow(pips_mat) > 1)) {
    res$cluster_weight <- as.numeric(crossprod(softmax_weights(elbo_vec, temperature = softmax_temperature), pips_mat))
  }
  if ("cluster_weight_050" %in% methods && !(length(elbo_vec) > 1 && nrow(pips_mat) > 1)) {
    res$cluster_weight_050 <- as.numeric(crossprod(softmax_weights(elbo_vec, temperature = softmax_temperature), pips_mat))
  }
  if ("cluster_weight_jsd_050" %in% methods && !(length(elbo_vec) > 1 && nrow(pips_mat) > 1)) {
    res$cluster_weight_jsd_050 <- as.numeric(crossprod(softmax_weights(elbo_vec, temperature = softmax_temperature), pips_mat))
  }
  res
}

#' Compute cluster-then-weight representative indices and weights.
#'
#' Returns a list with rep_idx (which rows are cluster representatives) and
#' w_rep (their normalized importance weights). Used by both PIP aggregation
#' and hg2 aggregation.
#' @keywords internal
.cluster_weights_from_hc <- function(hc, threshold, elbo_vec,
                                     n_fits,
                                     softmax_temperature = 1) {
  clusters <- stats::cutree(hc, h = threshold)
  cl_levels <- sort(unique(clusters))
  rep_idx <- integer(length(cl_levels))
  freq <- numeric(length(cl_levels))
  elbo_rep <- numeric(length(cl_levels))
  for (k in seq_along(cl_levels)) {
    cid <- cl_levels[[k]]
    idxs <- which(clusters == cid)
    freq[k] <- length(idxs) / n_fits
    pick <- idxs[which.max(elbo_vec[idxs])]
    rep_idx[k] <- pick
    elbo_rep[k] <- elbo_vec[pick]
  }
  w_rep <- softmax_weights(elbo_rep, temperature = softmax_temperature) / freq
  w_rep <- w_rep / sum(w_rep)
  list(rep_idx = rep_idx, w_rep = w_rep,
       ess = as.numeric(1 / sum(w_rep^2)))
}

#' Apply cluster-then-weight aggregation from a pre-computed hclust object.
#' @keywords internal
.cluster_weight_from_hc <- function(hc, threshold, elbo_vec, pips_mat,
                                    softmax_temperature = 1) {
  cw <- .cluster_weights_from_hc(hc, threshold, elbo_vec,
                                  n_fits = nrow(pips_mat),
                                  softmax_temperature = softmax_temperature)
  ens <- as.numeric(crossprod(cw$w_rep, pips_mat[cw$rep_idx, , drop = FALSE]))
  attr(ens, "ess") <- cw$ess
  ens
}

js_distance <- function(p, q, eps = 1e-12) {
  p <- as.numeric(p); q <- as.numeric(q)
  p <- p + eps; q <- q + eps
  m <- 0.5 * (p + q)
  kl <- function(a, b) sum(a * log(a / b))
  0.5 * kl(p, m) + 0.5 * kl(q, m)
}

compute_multimodal_metrics <- function(pip_list,
                                       jsd_threshold = 0.15,
                                       top_k = 10L,
                                       pip_cache = NULL) {
  n <- length(pip_list)
  if (n < 2) {
    return(tibble::tibble(
      mean_jsd = NA_real_,
      median_jsd = NA_real_,
      max_jsd = NA_real_,
      jaccard_top10 = NA_real_,
      mean_pip_var = NA_real_,
      n_clusters = NA_integer_,
      n_clusters_jsd_050 = NA_integer_
    ))
  }
  pips_mat <- pip_cache$pips_mat %||% do.call(rbind, lapply(pip_list, as.numeric))
  pip_cache <- pip_cache %||% prepare_pip_similarity_cache(pips_mat)
  jsd_vals <- pip_cache$jsd_vals

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
  hc_jsd <- pip_cache$hc

  tibble::tibble(
    mean_jsd = mean(jsd_vals),
    median_jsd = stats::median(jsd_vals),
    max_jsd = max(jsd_vals),
    jaccard_top10 = mean(pair_jaccard, na.rm = TRUE),
    mean_pip_var = pip_var,
    n_clusters = length(unique(stats::cutree(hc_jsd, h = jsd_threshold))),
    n_clusters_jsd_050 = length(unique(stats::cutree(hc_jsd, h = 0.50)))
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

write_multimodal_metrics <- function(bundle_row, use_case_id, group_label, explore_method = NA_character_,
                                     metrics, job_config, buffer_ctx = NULL) {
  if (is.null(metrics) || !nrow(metrics)) return(invisible(NULL))
  # Slim multimodal_metrics: only dataset_bundle_id + use_case_id +
  # group_label + explore_method + metrics. Bundle-level columns
  # (data_scenario, matrix_id, y_noise, p_star, phenotype_seed) are
  # stripped — re-attached during aggregation.
  out <- metrics %>%
    dplyr::mutate(
      dataset_bundle_id = bundle_row$dataset_bundle_id,
      use_case_id       = use_case_id,
      group_label       = group_label,
      explore_method    = as.character(explore_method)
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

write_refine_depth <- function(rd_tbl, job_config, buffer_ctx = NULL) {
  if (is.null(rd_tbl) || !nrow(rd_tbl)) return(invisible(NULL))
  verbose_output <- isTRUE(job_config$job$verbose_file_output)
  if (verbose_output) {
    base_output <- determine_base_output(job_config)
    run_dir <- file.path(base_output, "refine_depth")
    ensure_dir(run_dir)
    out_path <- file.path(run_dir, sprintf("refine_depth_%s.csv",
                          paste0(rd_tbl$dataset_bundle_id[[1]], "_", rd_tbl$use_case_id[[1]])))
    readr::write_csv(rd_tbl, out_path)
  } else {
    buffer_ctx$buffers$refine_depth[[length(buffer_ctx$buffers$refine_depth) + 1L]] <- rd_tbl
  }
}

write_confusion_bins <- function(run_row, job_config, conf_bins, buffer_ctx = NULL, variant_meta = NULL) {
  if (is.null(conf_bins) || !nrow(conf_bins)) return(invisible(NULL))
  # Slim confusion_bins: run_id always set (agg_method discriminates aggregated

  # rows). Run-table columns (dataset_bundle_id, architecture, use_case_id,
  # group_key, L, annotation_r2, inflate_match, sigma_0_2_scalar, refine_step,
  # c_value, tau_value, run_type) are stripped — re-attached during aggregation
  # via join on run_id.
  meta <- tibble::tibble(
    run_id = as.integer(run_row$run_id),
    explore_method = as.character(variant_meta$explore_method %||% "base"),
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
    variant_meta <- list(explore_method = "base", variant_id = NA, agg_method = NA)
  }
  # Runtime-derived columns NOT recoverable from run_table via run_id.
  # Run-table columns (task_id, use_case_id, phenotype_seed, dataset_bundle_id,
  # architecture, refine_step, run_type) are stripped here and re-attached
  # during aggregation via join on run_id.
  meta_tbl <- tibble::tibble(
    explore_method = as.character(variant_meta$explore_method %||% "base"),
    variant_id = as.character(variant_meta$variant_id %||% NA_character_),
    agg_method = as.character(variant_meta$agg_method %||% NA_character_),
    wall_time_sec = as.numeric(variant_meta$wall_time_sec %||% NA_real_),
    model_call_executed = as.logical(variant_meta$model_call_executed %||% NA),
    cache_source_run_id = as.integer(variant_meta$cache_source_run_id %||% NA_integer_)
  )
  # Slim version for effect_metrics — drop run-level accounting fields
  effect_meta_tbl <- dplyr::select(
    meta_tbl,
    -dplyr::any_of(c("wall_time_sec", "model_call_executed", "cache_source_run_id"))
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
    elbo_final = elbo_final
  ) %>% dplyr::bind_cols(meta_tbl)
  if (verbose_output) {
    readr::write_csv(model_metrics, file.path(run_dir, "model_metrics.csv"))
  }

  effect_metrics <- NULL
  if (!is.null(evaluation$effects_filtered) && nrow(evaluation$effects_filtered)) {
    effect_metrics <- dplyr::mutate(
      dplyr::select(evaluation$effects_filtered, -dplyr::any_of("indices")),
      run_id = run_row$run_id
    ) %>% dplyr::bind_cols(effect_meta_tbl)
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
    # Top-8 causal mask: score only the 8 largest-effect causal variants
    top8_mask_local <- if (length(data_bundle$causal_idx) > 0L && !is.null(data_bundle$beta)) {
      n_top <- min(8L, length(data_bundle$causal_idx))
      data_bundle$causal_idx[
        order(abs(data_bundle$beta[data_bundle$causal_idx]), decreasing = TRUE)
      ][seq_len(n_top)]
    } else {
      data_bundle$causal_idx
    }
    conf_bins <- compute_confusion_bins(
      evaluation$combined_pip,
      as.integer(seq_along(data_bundle$beta) %in% data_bundle$causal_idx),
      pip_bucket_width = job_config$job$metrics$pip_bucket_width %||% 0.01,
      pip_breaks = job_config$job$metrics$pip_breaks,
      causal_mask = top8_mask_local
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

prepare_task_staging_dir <- function(base_output, task_id) {
  staging_dir <- file.path(
    base_output,
    sprintf("task-%03d", as.integer(task_id))
  )
  ensure_dir(staging_dir)
  stale_flush_files <- list.files(
    staging_dir,
    pattern = "^flush-[0-9]+_",
    full.names = TRUE
  )
  if (!length(stale_flush_files)) {
    return(invisible(staging_dir))
  }
  rm_ok <- file.remove(stale_flush_files)
  if (!all(rm_ok)) {
    failed <- stale_flush_files[!rm_ok]
    stop(
      "Unable to clear existing flush files for task ",
      as.integer(task_id),
      ". One or more files are likely open/locked by another process. ",
      "Close those files and rerun.\n",
      paste(failed, collapse = "\n")
    )
  }
  invisible(staging_dir)
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

# Compute scaling confusion bins for one group (one spec × one dataset).
# Called from execute_dataset_bundle when write_scaling_confusion_bins = TRUE.
# Returns a tibble with per-(n_ensemble, agg_method, rep) pooled confusion bins
# including spec_name, dataset_bundle_id, annotation_r2, resolution columns.
#' @keywords internal
compute_scaling_confusion_bins_for_group <- function(pip_list,
                                                      elbo_vec,
                                                      meta_list,
                                                      causal_vec,
                                                      causal_mask = NULL,
                                                      group_run_row,
                                                      n_ens_sizes = c(4L, 8L, 16L, 32L, 64L),
                                                      n_restart_reps = 50L,
                                                      pip_breaks = NULL,
                                                      pip_bw = 0.01,
                                                      dataset_bundle_id = NA_character_,
                                                      full_agg_cache = NULL) {
  N <- length(pip_list)
  pip_mat <- do.call(cbind, pip_list)  # p x N

  run_meta <- dplyr::bind_rows(lapply(meta_list, tibble::as_tibble))

  spec_nm   <- as.character(run_meta$spec_name[[1]] %||% NA_character_)
  r2_val    <- as.numeric(run_meta$annotation_r2[[1]] %||% NA_real_)

  expl_raw  <- as.character(group_run_row$exploration_methods[[1]] %||% "")
  method_ids <- unique(trimws(strsplit(expl_raw, "x", fixed = TRUE)[[1]]))
  is_interaction <- length(method_ids) > 1L

  agg_methods <- c("uniform", "max_elbo", "elbo_softmax", "cluster_weight", "cluster_weight_050")

  results <- list()

  select_axis_subset <- function(vals, n_keep, axis_col) {
    vals <- sort(unique(vals))
    if (!length(vals)) return(vals)
    n_keep <- min(as.integer(n_keep), length(vals))
    if (n_keep >= length(vals)) return(vals)
    if (axis_col %in% c("restart_id", "refine_step")) {
      return(vals[seq_len(n_keep)])
    }
    vals[unique(round(seq(1, length(vals), length.out = n_keep)))]
  }

  compute_bins_for_subset <- function(sel_idx, agg_method) {
    is_full_subset <- length(sel_idx) == N && identical(sort(sel_idx), seq_len(N))
    if (is_full_subset && !is.null(full_agg_cache$bins[[agg_method]])) {
      return(full_agg_cache$bins[[agg_method]])
    }

    sub_pip_mat <- pip_mat[, sel_idx, drop = FALSE]
    sub_elbos   <- elbo_vec[sel_idx]
    subset_cache <- NULL
    if (agg_method %in% c("cluster_weight", "cluster_weight_050")) {
      subset_cache <- prepare_pip_similarity_cache(t(sub_pip_mat))
    }
    agg_pip <- tryCatch(
      aggregate_pip_matrix(sub_pip_mat, sub_elbos, agg_method,
                           hc = subset_cache$hc %||% NULL),
      error = function(e) NULL
    )
    if (is.null(agg_pip)) return(NULL)
    compute_confusion_bins(
      agg_pip, causal_vec,
      pip_bucket_width = pip_bw, pip_breaks = pip_breaks,
      causal_mask = causal_mask
    )
  }

  if (is_interaction) {
    # Skip 3+ axis specs — only 2-axis interactions are analyzed for subscaling.
    if (length(method_ids) > 2L) return(dplyr::bind_rows(results))

    # 2-axis interaction: proportional halving.  Full grid = n1 x n2;
    # half grid = ceil(n1/2) x ceil(n2/2).  Labels use PLANNED sizes
    # (from the run table) so all datasets get identical labels even
    # when BFS truncates some refine trees early.
    axis_info <- list(
      if ("restart"        %in% method_ids) list(col = "restart_id",       planned_col = ".planned_n_restart") else NULL,
      if ("refine"         %in% method_ids) list(col = "refine_step",      planned_col = ".planned_n_refine")  else NULL,
      if ("c_grid"         %in% method_ids) list(col = "c_value",          planned_col = ".planned_n_c")       else NULL,
      if ("sigma_0_2_grid" %in% method_ids) list(col = "sigma_0_2_scalar", planned_col = ".planned_n_sigma")   else NULL
    )
    axis_info <- Filter(Negate(is.null), axis_info)

    # Use planned sizes from the run table (falls back to actual if unavailable)
    axis_planned_sizes <- vapply(axis_info, function(ax) {
      planned <- as.integer(group_run_row[[ax$planned_col]] %||% NA_integer_)
      if (is.na(planned) || planned < 1L) {
        as.integer(length(unique(run_meta[[ax$col]])))
      } else {
        planned
      }
    }, integer(1L))

    # Sort axes by descending planned size so labels are consistent
    # (larger axis first, e.g. "8x4" not "4x8")
    sort_order <- order(axis_planned_sizes, decreasing = TRUE)
    axis_info <- axis_info[sort_order]
    axis_planned_sizes <- axis_planned_sizes[sort_order]

    # Two resolution levels: full and half (proportional per axis)
    resolutions <- list(
      half = as.integer(pmax(1L, ceiling(axis_planned_sizes / 2))),
      full = axis_planned_sizes
    )

    for (res_name in names(resolutions)) {
      res_sizes <- resolutions[[res_name]]
      sub_meta <- run_meta
      if (!identical(res_name, "full")) {
        for (ai in seq_along(axis_info)) {
          col <- axis_info[[ai]]$col
          if (!col %in% names(sub_meta)) next
          keep_vals <- select_axis_subset(sub_meta[[col]], res_sizes[ai], col)
          sub_meta <- dplyr::filter(sub_meta, .data[[col]] %in% keep_vals)
        }
      }
      # Label uses planned sizes so all datasets get identical labels
      resolution_label <- paste(res_sizes, collapse = "x")

      sel_idx <- match(sub_meta$run_id, run_meta$run_id)
      sel_idx <- sel_idx[!is.na(sel_idx)]
      if (!length(sel_idx)) next
      for (am in agg_methods) {
        bins <- compute_bins_for_subset(sel_idx, am)
        if (!is.null(bins) && nrow(bins)) {
          results[[length(results) + 1L]] <- dplyr::mutate(bins,
            spec_name         = spec_nm,
            annotation_r2     = r2_val,
            dataset_bundle_id = as.character(dataset_bundle_id),
            n_ensemble        = as.integer(length(sel_idx)),
            resolution        = resolution_label,
            agg_method        = am,
            rep               = 1L
          )
        }
      }
    }
  } else {
    # Pure-lever specs: subsample N fits across valid sizes.
    lever_type <- if ("restart"        %in% method_ids) "restart"   else
                  if ("refine"         %in% method_ids) "refine"    else
                  if ("sigma_0_2_grid" %in% method_ids) "sigma_0_2" else
                  if ("c_grid"         %in% method_ids) "c_grid"    else NA_character_
    if (is.na(lever_type)) return(NULL)

    n_reps <- 1L
    valid_sizes <- n_ens_sizes[n_ens_sizes <= N]

    for (n_ens in valid_sizes) {
      for (rep_i in seq_len(n_reps)) {
        sel_ids <- select_subsample_run_ids(run_meta, n_ens, lever_type, rep_i)
        sel_idx <- match(sel_ids, run_meta$run_id)
        sel_idx <- sel_idx[!is.na(sel_idx)]
        if (!length(sel_idx)) next
        for (am in agg_methods) {
          bins <- compute_bins_for_subset(sel_idx, am)
          if (!is.null(bins) && nrow(bins)) {
            results[[length(results) + 1L]] <- dplyr::mutate(bins,
              spec_name         = spec_nm,
              annotation_r2     = r2_val,
              dataset_bundle_id = as.character(dataset_bundle_id),
              n_ensemble        = as.integer(n_ens),
              resolution        = NA_character_,
              agg_method        = am,
              rep               = as.integer(rep_i)
            )
          }
        }
      }
    }
  }

  if (!length(results)) return(NULL)
  dplyr::bind_rows(results)
}

#' Compute estimated heritability per aggregation method for one ensemble group.
#'
#' Uses the full pool of K fits (no subsampling). Weights are derived from PIP
#' vectors (same JSD-based clustering used in confusion bins), then applied to
#' fitted_y vectors to produce an exact aggregated hg2 scalar per method.
#'
#' @keywords internal
compute_hg2_by_agg <- function(pip_list,
                                fitted_y_list,
                                elbo_vec,
                                y,
                                group_run_row,
                                dataset_bundle_id = NA_character_,
                                pip_cache = NULL) {
  valid <- !vapply(fitted_y_list, is.null, logical(1L))
  pip_list      <- pip_list[valid]
  fitted_y_list <- fitted_y_list[valid]
  elbo_vec      <- elbo_vec[valid]
  K <- length(fitted_y_list)
  if (K == 0L || is.null(y)) return(NULL)
  vy <- stats::var(as.numeric(y))
  if (!is.finite(vy) || vy <= 0) return(NULL)

  # K x n matrix of fitted values; K x p matrix of PIPs (for JSD clustering)
  fy_mat  <- do.call(rbind, lapply(fitted_y_list, as.numeric))  # K x n
  pip_mat <- pip_cache$pips_mat %||% do.call(rbind, lapply(pip_list, as.numeric))       # K x p

  # JSD-based hierarchical clustering on PIP vectors (same as confusion bins)
  pip_cache <- pip_cache %||% if (K > 1L) prepare_pip_similarity_cache(pip_mat) else NULL
  hc <- pip_cache$hc %||% NULL

  agg_methods <- c("uniform", "max_elbo", "elbo_softmax",
                   "cluster_weight", "cluster_weight_050")
  rows <- lapply(agg_methods, function(am) {
    if (am %in% c("cluster_weight", "cluster_weight_050") && !is.null(hc)) {
      thr <- if (am == "cluster_weight") 0.15 else 0.50
      cw <- .cluster_weights_from_hc(hc, thr, elbo_vec,
                                      n_fits = K)
      fy_agg <- as.numeric(crossprod(cw$w_rep,
                                      fy_mat[cw$rep_idx, , drop = FALSE]))
    } else {
      w <- switch(am,
        "uniform"           = rep(1 / K, K),
        "max_elbo"          = {
          idx <- which.max(elbo_vec)
          w2 <- rep(0, K); w2[idx] <- 1; w2
        },
        "elbo_softmax"      = softmax_weights(elbo_vec),
        # fallback for cluster methods when hc is NULL
        "cluster_weight"    = softmax_weights(elbo_vec),
        "cluster_weight_050"= softmax_weights(elbo_vec)
      )
      fy_agg <- as.numeric(crossprod(w, fy_mat))
    }
    hg2 <- max(0, min(1, stats::var(fy_agg) / vy))
    tibble::tibble(agg_method = am, hg2 = hg2, n_fits = K)
  })

  dplyr::bind_rows(rows) %>%
    dplyr::mutate(
      dataset_bundle_id = as.character(dataset_bundle_id),
      use_case_id       = as.character(group_run_row$use_case_id %||% NA_character_),
      group_label       = as.character(
        if ("group_label" %in% names(group_run_row)) group_run_row$group_label
        else paste0("model_grid|", group_run_row$group_key %||% NA_character_)
      ),
      explore_method    = as.character(
        if ("explore_method" %in% names(group_run_row)) group_run_row$explore_method
        else NA_character_
      ),
      variant_id        = as.character(
        if ("variant_id" %in% names(group_run_row)) group_run_row$variant_id
        else NA_character_
      )
    )
}

flush_task_buffers <- function(buffer_ctx) {
  buf <- buffer_ctx$buffers
  if (is.null(buf) ||
      !length(buf$model) && !length(buf$effect) && !length(buf$snps) &&
      !length(buf$confusion) && !length(buf$dataset_metrics) && !length(buf$multimodal) &&
      !length(buf$refine_depth) && !length(buf$scaling_bins) &&
      !length(buf$validation) && !length(buf$prior_diagnostics) &&
      !length(buf$tier_cs_metrics) && !length(buf$hg2_by_agg)) {
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

  model_tbl <- if (length(buf$model)) dplyr::bind_rows(buf$model) else NULL
  effect_tbl <- if (length(buf$effect)) dplyr::bind_rows(buf$effect) else NULL
  snp_tbl <- if (length(buf$snps)) dplyr::bind_rows(buf$snps) else NULL
  validation_tbl <- if (length(buf$validation)) dplyr::bind_rows(buf$validation) else NULL
  confusion_tbl <- if (length(buf$confusion)) dplyr::bind_rows(buf$confusion) else NULL
  dataset_tbl <- if (length(buf$dataset_metrics)) dplyr::bind_rows(buf$dataset_metrics) else NULL
  multimodal_tbl <- if (length(buf$multimodal)) dplyr::bind_rows(buf$multimodal) else NULL
  refine_depth_tbl   <- if (length(buf$refine_depth))      dplyr::bind_rows(buf$refine_depth)      else NULL
  scaling_bins_tbl   <- if (length(buf$scaling_bins))      dplyr::bind_rows(buf$scaling_bins)      else NULL
  prior_diag_tbl     <- if (length(buf$prior_diagnostics)) dplyr::bind_rows(buf$prior_diagnostics) else NULL
  tier_cs_tbl        <- if (length(buf$tier_cs_metrics))   dplyr::bind_rows(buf$tier_cs_metrics)   else NULL
  hg2_by_agg_tbl     <- if (length(buf$hg2_by_agg))        dplyr::bind_rows(buf$hg2_by_agg)        else NULL

  write_flush_outputs(
    staging_dir = staging_dir,
    flush_label = flush_label,
    task_id = buffer_ctx$task_id,
    model_metrics = model_tbl,
    effect_metrics = effect_tbl,
    snp_tbl = snp_tbl,
    validation_row = validation_tbl,
    confusion_bins = confusion_tbl,
    dataset_metrics = dataset_tbl,
    multimodal_metrics = multimodal_tbl,
    refine_depth = refine_depth_tbl,
    scaling_bins = scaling_bins_tbl,
    prior_diagnostics = prior_diag_tbl,
    tier_cs_metrics = tier_cs_tbl,
    hg2_by_agg = hg2_by_agg_tbl
  )

  if (is.null(validation_tbl) || !nrow(validation_tbl)) {
    warning(
      sprintf(
        "Flush %s has no validation rows; %s_validation.csv was not written.",
        flush_label,
        flush_label
      )
    )
  }

  buffer_ctx$buffers <- list(
    model = list(),
    effect = list(),
    snps = list(),
    confusion = list(),
    dataset_metrics = list(),
    multimodal = list(),
    refine_depth = list(),
    scaling_bins = list(),
    validation = list(),
    prior_diagnostics = list(),
    tier_cs_metrics = list(),
    hg2_by_agg = list()
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
                                multimodal_metrics = NULL,
                                refine_depth = NULL,
                                scaling_bins = NULL,
                                prior_diagnostics = NULL,
                                tier_cs_metrics = NULL,
                                hg2_by_agg = NULL) {
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
  if (!is.null(refine_depth) && nrow(refine_depth)) {
    refine_depth <- dplyr::mutate(refine_depth, task_id = task_id, flush_id = flush_label)
    readr::write_csv(refine_depth, file.path(staging_dir, sprintf("%s_refine_depth.csv", flush_label)))
  }
  if (!is.null(scaling_bins) && nrow(scaling_bins)) {
    scaling_bins <- dplyr::mutate(scaling_bins, task_id = task_id, flush_id = flush_label)
    readr::write_csv(scaling_bins, file.path(staging_dir, sprintf("%s_scaling_bins.csv", flush_label)))
  }
  if (!is.null(prior_diagnostics) && nrow(prior_diagnostics)) {
    prior_diagnostics <- dplyr::mutate(prior_diagnostics, task_id = task_id, flush_id = flush_label)
    readr::write_csv(prior_diagnostics, file.path(staging_dir, sprintf("%s_prior_diagnostics.csv", flush_label)))
  }
  if (!is.null(tier_cs_metrics) && nrow(tier_cs_metrics)) {
    tier_cs_metrics <- dplyr::mutate(tier_cs_metrics, task_id = task_id, flush_id = flush_label)
    readr::write_csv(tier_cs_metrics, file.path(staging_dir, sprintf("%s_tier_cs_metrics.csv", flush_label)))
  }
  if (!is.null(hg2_by_agg) && nrow(hg2_by_agg)) {
    hg2_by_agg <- dplyr::mutate(hg2_by_agg, task_id = task_id, flush_id = flush_label)
    readr::write_csv(hg2_by_agg, file.path(staging_dir, sprintf("%s_hg2_by_agg.csv", flush_label)))
  }
}

# Write SNP parquet with reduced precision for numeric columns to cut size.
write_compact_snp_parquet <- function(snp_tbl, path) {
  # Ensure flush_id exists even if caller did not add it.
  if (!"flush_id" %in% names(snp_tbl)) {
    snp_tbl$flush_id <- NA_character_
  }
  # Slim parquet: run-table columns (task_id, dataset_bundle_id, architecture,
  # refine_step, run_type) stripped — recoverable via run_id join.
  snp_tbl <- dplyr::mutate(
    snp_tbl,
    run_id = as.integer(run_id),
    snp_index = as.integer(snp_index),
    pip = as.numeric(pip),
    beta = as.numeric(beta),
    mu_0 = as.numeric(mu_0),
    sigma_0_2 = as.numeric(sigma_0_2),
    causal = as.integer(causal),
    flush_id = as.character(flush_id),
    explore_method = as.character(explore_method %||% variant_type %||% NA_character_),
    variant_id = as.character(variant_id %||% NA_character_),
    agg_method = as.character(agg_method %||% NA_character_)
  )
  # Fixed schema — slim (run-table cols stripped).
  snp_schema <- arrow::schema(
    run_id = arrow::int32(),
    snp_index = arrow::int32(),
    pip = arrow::float32(),
    beta = arrow::float32(),
    mu_0 = arrow::float32(),
    sigma_0_2 = arrow::float32(),
    causal = arrow::int8(),
    flush_id = arrow::utf8(),
    explore_method = arrow::utf8(),
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
