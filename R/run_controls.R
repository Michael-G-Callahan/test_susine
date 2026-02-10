# Run-control builders ------------------------------------------------------

#' Build a tidy grid of prior-quality settings.
#'
#' @param annotation_r2_levels Numeric vector of target annotation R^2 values for causal SNPs.
#' @param inflate_match_levels Numeric vector indicating how strongly to inflate non-causal annotations.
#' @param gamma_shrink_levels Numeric vector of shrinkage slopes for variance-informed runs.
#' @return Tibble with columns `annotation_r2`, `inflate_match`, `gamma_shrink`, and
#'   `prior_quality_id`.
#' @export
prior_quality_grid <- function(annotation_r2_levels = c(0.25, 0.5, 0.75),
                               inflate_match_levels = c(0, 0.5, 1),
                               gamma_shrink_levels = c(0.4, 0.8)) {
  ann_vals <- unique(annotation_r2_levels)
  inflate_vals <- unique(inflate_match_levels)
  gamma_vals <- unique(gamma_shrink_levels)
  if (!any(is.na(gamma_vals))) {
    gamma_vals <- c(gamma_vals, NA_real_)
  }
  grid <- tidyr::expand_grid(
    annotation_r2 = ann_vals,
    inflate_match = inflate_vals,
    gamma_shrink = gamma_vals
  )
  dplyr::mutate(grid, prior_quality_id = dplyr::row_number())
}

#' Convert an absolute path to one relative to the repository root.
#' @keywords internal
relativize_path <- function(path, root) {
  if (is.null(path)) {
    return(NA_character_)
  }
  root_norm <- normalizePath(root, winslash = "/", mustWork = TRUE)
  purrr::map_chr(path, function(p) {
    if (is.na(p) || !nzchar(p)) {
      return(NA_character_)
    }
    path_norm <- normalizePath(p, winslash = "/", mustWork = FALSE)
    path_norm <- gsub("\\\\", "/", path_norm)
    prefix <- paste0(root_norm, "/")
    if (startsWith(path_norm, prefix)) {
      substring(path_norm, nchar(prefix) + 1L)
    } else {
      repo_base <- basename(root_norm)
      pattern <- paste0("/", repo_base, "/")
      match_idx <- regexpr(pattern, path_norm, fixed = TRUE)
      if (match_idx > 0) {
        start <- match_idx + nchar(pattern)
        substring(path_norm, start)
      } else {
        path_norm
      }
    }
  })
}

#' Build a catalog of available design matrices for the requested scenarios.
#'
#' @param requested_scenarios Character vector.
#' @param repo_root Root of the repository (used to relativize paths).
#' @param summary_path Optional CSV produced by the sampling workbook that lists
#'   `matrix_path`/`manifest_path` pairs for scenarios 1-4.
#' @return Tibble with one row per matrix, including `matrix_id`, `data_scenario`,
#'   `matrix_path`, and metadata such as participant/snps counts.
#' @keywords internal
build_data_matrix_catalog <- function(requested_scenarios,
                                      repo_root,
                                      summary_path = NULL) {
  repo_root <- normalizePath(repo_root, winslash = "/", mustWork = TRUE)
  scenarios <- unique(requested_scenarios)
  summary_tbl <- NULL
  if (!is.null(summary_path) && file.exists(summary_path)) {
    summary_tbl <- readr::read_csv(summary_path, show_col_types = FALSE)
  }

  make_entry <- function(scenario) {
    scenario_lower <- tolower(scenario)
    if (scenario_lower %in% c("simulation_n3", "sim_n3")) {
      tibble::tibble(
        data_scenario = scenario,
        dataset_label = "simulation_n3",
        participant_count = NA_integer_,
        snps_post = NA_integer_,
        snp_set = NA_character_,
        matrix_path = NA_character_,
        manifest_path = NA_character_,
        source = "simulation"
      )
    } else {
      if (is.null(summary_tbl) || !nrow(summary_tbl)) {
        stop(
          "No sampled scenario summary found. ",
          "Provide `summary_path` when calling make_job_config() ",
          "or run the sampling workbook first."
        )
      }
      subset <- dplyr::filter(summary_tbl, .data$scenario == !!scenario)
      if (!nrow(subset)) {
        stop("Scenario '", scenario, "' is not present in ", summary_path, ".")
      }
      subset %>%
        dplyr::transmute(
          data_scenario = .data$scenario,
          dataset_label = .data$gene,
          participant_count = as.integer(.data$participant_count),
          snps_post = as.integer(.data$snps_post),
          snp_set = .data$snp_set,
          matrix_path = relativize_path(.data$matrix_path, repo_root),
          manifest_path = relativize_path(.data$manifest_path, repo_root),
          source = "sampled"
        )
    }
  }

  catalog <- purrr::map_dfr(scenarios, make_entry) %>%
    dplyr::group_by(.data$data_scenario) %>%
    dplyr::mutate(matrix_index = dplyr::row_number()) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(matrix_id = dplyr::row_number())

  catalog
}

ensure_repo_root <- function(path) {
  current <- normalizePath(path, winslash = "/", mustWork = TRUE)
  sentinels <- c(".here", "DESCRIPTION")
  repeat {
    if (any(file.exists(file.path(current, sentinels)))) {
      return(current)
    }
    parent <- dirname(current)
    if (identical(parent, current)) {
      return(current)
    }
    current <- parent
  }
}

#' Construct the full run table for a job.
#'
#' @param use_case_ids Character vector of use-case identifiers.
#' @param L_grid Integer vector of SuSiE/SuSiNE L values.
#' @param y_noise_grid Numeric vector of noise fractions (0-1).
#' @param prior_quality Tibble from [prior_quality_grid()].
#' @param p_star_grid Integer vector for number of causal SNPs.
#' @param seeds Integer vector of RNG seeds for phenotype simulation.
#' @param data_scenarios Character vector naming the data sources.
#' @param pair_L_p_star Logical; when TRUE (and `grid_mode == "full"`), the grid pairs matching
#'   `L` and `p_star` values instead of taking their Cartesian product.
#' @return List with elements `dataset_bundles`, `runs`, and `tasks`.
#' @keywords internal
make_run_tables <- function(use_case_ids,
                            L_grid,
                            y_noise_grid,
                            prior_quality,
                            p_star_grid,
                            seeds,
                            data_scenarios,
                            grid_mode = c("full", "minimal", "intersect"),
                            pair_L_p_star = FALSE,
                            data_matrix_catalog = NULL) {
  if (!is.data.frame(prior_quality) ||
      !all(c("annotation_r2", "inflate_match", "gamma_shrink") %in% names(prior_quality))) {
    stop("prior_quality must contain annotation_r2, inflate_match, and gamma_shrink columns.")
  }
  if (is.null(data_matrix_catalog) || !nrow(data_matrix_catalog)) {
    stop("data_matrix_catalog must be supplied with one row per available matrix.")
  }
  grid_mode <- match.arg(grid_mode)
  use_cases <- resolve_use_cases(use_case_ids)
  if (!nrow(use_cases)) {
    stop("No valid use cases selected.")
  }

  use_cases <- use_cases %>%
    dplyr::mutate(
      uses_mu_annotation = .data$mu_strategy %in% c("functional", "eb_mu"),
      uses_sigma_annotation = .data$sigma_strategy %in% c("functional", "eb_sigma"),
      uses_gamma_shrink = .data$sigma_strategy %in% c("functional", "eb_sigma")
    )

  matrix_catalog_subset <- dplyr::filter(data_matrix_catalog, .data$data_scenario %in% unique(data_scenarios))
  if (!nrow(matrix_catalog_subset)) {
    stop("No matrix metadata found for the requested scenarios.")
  }

  seed_values <- unique(seeds)
  if (!length(seed_values)) {
    stop("At least one seed must be supplied.")
  }

  build_dataset_bundles <- function() {
    scenario_vals <- unique(data_scenarios)
    y_vals <- unique(y_noise_grid)
    p_vals <- unique(p_star_grid)
    candidate_matrices <- matrix_catalog_subset %>%
      dplyr::arrange(.data$data_scenario, .data$matrix_id) %>%
      dplyr::pull(.data$matrix_id)

    if (grid_mode == "minimal") {
      # Single bundle per scenario (first matrix, first y_noise, first p_star, first seed).
      matrix_ids <- matrix_catalog_subset %>%
        dplyr::group_by(.data$data_scenario) %>%
        dplyr::slice_head(n = 1) %>%
        dplyr::ungroup() %>%
        dplyr::pull(.data$matrix_id)
      bundles <- tibble::tibble(
        matrix_id = matrix_ids,
        y_noise = y_vals[1],
        p_star = p_vals[1],
        phenotype_seed = seed_values[1]
      )
      return(bundles)
    }

    if (grid_mode == "intersect") {
      values <- list(
        y_noise = y_vals,
        p_star = p_vals,
        phenotype_seed = seed_values,
        matrix_id = candidate_matrices
      )
      lengths <- vapply(values, length, integer(1))
      n_rows <- max(lengths)
      bundles <- tibble::tibble(
        y_noise = rep_len(values$y_noise, n_rows),
        p_star = rep_len(values$p_star, n_rows),
        phenotype_seed = rep_len(values$phenotype_seed, n_rows),
        matrix_id = rep_len(values$matrix_id, n_rows)
      )
      return(bundles)
    }

    # Full grid
    bundles <- tidyr::expand_grid(
      matrix_id = candidate_matrices,
      y_noise = y_vals,
      p_star = p_vals,
      phenotype_seed = seed_values
    )
    bundles
  }

  dataset_bundles <- build_dataset_bundles() %>%
    dplyr::mutate(
      dataset_bundle_id = dplyr::row_number(),
      phenotype_seed = dplyr::row_number()
    )

  if (!nrow(dataset_bundles)) {
    stop("No dataset bundles created.")
  }

  unique_L <- unique(L_grid)
  unique_p <- unique(p_star_grid)
  if (isTRUE(pair_L_p_star) && !setequal(unique_L, unique_p)) {
    stop("When pair_L_p_star = TRUE, L_grid and p_star_grid must have the same unique values.")
  }

  L_by_bundle <- dataset_bundles %>%
    dplyr::distinct(dataset_bundle_id, p_star) %>%
    dplyr::mutate(
      L_vals = if (isTRUE(pair_L_p_star)) list(.data$p_star) else list(unique_L)
    ) %>%
    tidyr::unnest(cols = "L_vals") %>%
    dplyr::rename(L = L_vals) %>%
    dplyr::select(.data$dataset_bundle_id, .data$L)

  annotation_grid_for_use_case <- function(uc_row) {
    if (!isTRUE(uc_row$uses_mu_annotation) && !isTRUE(uc_row$uses_sigma_annotation)) {
      return(tibble::tibble(annotation_r2 = NA_real_, inflate_match = NA_real_, gamma_shrink = NA_real_))
    }
    grid <- prior_quality %>% dplyr::select(annotation_r2, inflate_match, gamma_shrink)
    if (!isTRUE(uc_row$uses_gamma_shrink)) {
      grid <- grid %>%
        dplyr::mutate(gamma_shrink = NA_real_) %>%
        dplyr::distinct()
    }
    grid
  }

  runs <- purrr::map_dfr(seq_len(nrow(use_cases)), function(i) {
    uc <- use_cases[i, ]
    ann_grid <- annotation_grid_for_use_case(uc)
    base <- dataset_bundles %>%
      dplyr::inner_join(L_by_bundle, by = "dataset_bundle_id") %>%
      dplyr::mutate(use_case_id = uc$use_case_id)
    tidyr::crossing(base, ann_grid)
  }) %>%
    dplyr::arrange(.data$dataset_bundle_id, .data$use_case_id, .data$L, .data$phenotype_seed) %>%
    dplyr::mutate(run_id = dplyr::row_number())

  list(
    dataset_bundles = dataset_bundles,
    data_matrices = data_matrix_catalog,
    runs = runs,
    use_cases = use_cases
  )
}

#' Attach task identifiers by dataset bundle.
#'
#' Ensures all runs for a given dataset_bundle_id are assigned to the same task.
#'
#' @param runs Run tibble containing a dataset bundle column.
#' @param bundle_id_col Column name containing the dataset bundle id.
#' @param bundles_per_task Integer number of bundles assigned to each task.
#' @param shuffle_seed Optional seed for reproducible shuffling.
#' @param shuffle Logical; when FALSE, assign bundles to tasks in id order.
#' @return List with `runs` (task_id assigned) and `tasks` summary.
#' @keywords internal
assign_task_ids_by_bundle <- function(runs,
                                      bundle_id_col = "dataset_bundle_id",
                                      bundles_per_task = 1,
                                      shuffle_seed = NULL,
                                      shuffle = TRUE) {
  if (bundles_per_task < 1) {
    stop("bundles_per_task must be >= 1")
  }
  if (!bundle_id_col %in% names(runs)) {
    stop("bundle_id_col not found in runs: ", bundle_id_col)
  }

  bundles <- runs %>%
    dplyr::distinct(.data[[bundle_id_col]]) %>%
    dplyr::rename(bundle_id = !!bundle_id_col)

  if (isTRUE(shuffle)) {
    if (!is.null(shuffle_seed)) {
      set.seed(shuffle_seed)
    }
    bundles <- bundles %>%
      dplyr::mutate(shuffled_order = sample(n())) %>%
      dplyr::arrange(.data$shuffled_order)
  } else {
    bundles <- bundles %>% dplyr::arrange(.data$bundle_id)
  }

  bundles <- bundles %>%
    dplyr::mutate(
      task_id = ((dplyr::row_number() - 1L) %/% as.integer(bundles_per_task)) + 1L
    )

  runs_out <- runs %>%
    dplyr::left_join(bundles %>% dplyr::select(.data$bundle_id, .data$task_id),
                     by = setNames("bundle_id", bundle_id_col))

  tasks_summary <- bundles %>%
    dplyr::group_by(.data$task_id) %>%
    dplyr::summarise(bundles_per_task = dplyr::n(), .groups = "drop") %>%
    dplyr::arrange(.data$task_id)

  list(runs = runs_out, tasks = tasks_summary)
}

#' Create an in-memory job configuration.
#'
#' @param job_name Character scalar.
#' @param use_case_ids Character vector of selected use cases.
#' @param L_grid Integer vector of L values.
#' @param y_noise_grid Numeric vector of noise fractions.
#' @param prior_quality Tibble with prior noise settings.
#' @param p_star_grid Integer vector.
#' @param seeds Integer vector of RNG seeds for phenotype simulation.
#' @param data_scenarios Character vector.
#' @param task_unit One of "dataset" or "run"; controls how tasks are assigned.
#' @param bundles_per_task Integer number of dataset bundles assigned to each task.
#' @param runs_per_task Integer number of runs assigned to each SLURM task (used when task_unit = "run").
#' @param email Notification email address.
#' @param output_root Root directory for outputs (default `output`).
#' @param credible_set_rho Credible set cumulative PIP threshold.
#' @param purity_threshold Minimum purity to keep CS in filtered summary.
#' @param verbose_file_output When FALSE, avoid per-run artifacts and stream
#'   metrics into task/flush staging files (reduces filesystem load).
#' @param write_legacy_snp_csv Logical; when TRUE, retain the legacy per-run
#'   `pip.csv` and `truth.csv` files alongside the new compact SNP table.
#' @param sigma_0_2_scalars Character or numeric vector of sigma_0_2 prior variance
#'   scalars to test. Only applies to use cases with `sigma_strategy == "naive"`.
#'   Supported formats:
#'   - Numeric values (e.g., 0.1, 0.2, 0.4): multiplier of var(y)
#'   - Strings with /L suffix (e.g., "1/L", "0.5/L", "2/L"): multiplier/L of var(y)
#'   When NULL (default), uses susine's internal default (1/L * var_y).
#' @param annotation_scales Numeric vector of multipliers applied to functional
#'   prior means (`mu_0`). Only use cases with `mu_strategy == "functional"`
#'   receive these values; all other runs are assigned `NA`. When NULL or empty,
#'   defaults to a single value of 1 (no scaling).
#' @param anneal_settings Named list for tempering runs. Set to NULL to omit
#'   annealing config from the job config.
#' @param model_average_settings Named list for model averaging runs.
#' @param pair_L_p_star Logical; when TRUE, the full grid pairs `L` and `p_star` values 1-1
#'   instead of expanding their Cartesian product.
#' @param task_assignment_seed Optional integer seed to make task assignment reproducible.
#'
#' @return Job configuration list ready to serialize.
#' @export
make_job_config <- function(job_name,
                            HPC = FALSE,
                            time = "02:59:00",
                            mem = "2G",
                            use_case_ids,
                            L_grid,
                            y_noise_grid,
                            prior_quality,
                            p_star_grid,
                            seeds,
                            data_scenarios = "simulation_n3",
                            repo_root = ".",
                            sampled_scenario_summary = NULL,
                            data_matrix_catalog = NULL,
                            task_unit = c("dataset", "run"),
                            bundles_per_task = 1,
                            runs_per_task = 150,
                            email = "mgc5166@psu.edu",
                            output_root = "output",
                            credible_set_rho = 0.95,
                            purity_threshold = 0.5,
                            write_legacy_snp_csv = FALSE,
                            write_snps_parquet = TRUE,
                            write_confusion_bins = TRUE,
                            verbose_file_output = TRUE,
                            buffer_flush_interval = 300L,
                            task_assignment_seed = NULL,
                            grid_mode = c("full", "minimal", "intersect"),
                            pair_L_p_star = FALSE,
                            sigma_0_2_scalars = NULL,
                            annotation_scales = NULL,
                            anneal_settings = list(
                              anneal_start_T = 5,
                              anneal_schedule_type = "geometric",
                              anneal_burn_in = 5
                            ),
                            model_average_settings = list(
                              n_inits = 5,
                              init_sd = 0.05
                            ),
                            restart_settings = list(
                              n_inits = 5,
                              alpha_concentration = 1
                            ),
                            aggregation_methods = c("softmax_elbo"),
                            metrics_settings = list(
                              pip_bucket_width = 0.01,
                              z_top_k = 10,
                              jsd_threshold = 0.171661
                            )) {

  grid_mode <- match.arg(grid_mode)
  task_unit <- match.arg(task_unit)
  if (task_unit != "dataset") {
    stop("task_unit must be 'dataset' (run-based task assignment is disabled).")
  }
  repo_root <- normalizePath(repo_root, winslash = "/", mustWork = TRUE)
  repo_root <- ensure_repo_root(repo_root)
  if (is.null(data_matrix_catalog)) {
    data_matrix_catalog <- build_data_matrix_catalog(
      requested_scenarios = unique(data_scenarios),
      repo_root = repo_root,
      summary_path = sampled_scenario_summary
    )
  }

  tables <- make_run_tables(
    use_case_ids = use_case_ids,
    L_grid = L_grid,
    y_noise_grid = y_noise_grid,
    prior_quality = prior_quality,
    p_star_grid = p_star_grid,
    seeds = seeds,
    data_scenarios = data_scenarios,
    grid_mode = grid_mode,
    pair_L_p_star = pair_L_p_star,
    data_matrix_catalog = data_matrix_catalog
  )

  use_cases <- tables$use_cases

  # Expand annotation scales across runs that use functional prior means
  functional_mu_ids <- use_cases$use_case_id[use_cases$mu_strategy == "functional"]
  scale_values <- annotation_scales
  if (is.null(scale_values) || length(scale_values) == 0) {
    scale_values <- 1
  }
  scale_values <- as.numeric(scale_values)

  runs_functional <- tables$runs %>%
    dplyr::filter(.data$use_case_id %in% functional_mu_ids)
  runs_non_functional <- tables$runs %>%
    dplyr::filter(!(.data$use_case_id %in% functional_mu_ids)) %>%
    dplyr::mutate(annotation_scale = NA_real_)

  if (nrow(runs_functional) > 0) {
    runs_functional <- runs_functional %>%
      tidyr::crossing(annotation_scale = scale_values)
  } else {
    runs_functional <- dplyr::mutate(runs_functional, annotation_scale = numeric(0))
  }

  tables$runs <- dplyr::bind_rows(runs_functional, runs_non_functional) %>%
    dplyr::arrange(.data$dataset_bundle_id, .data$use_case_id, .data$L, .data$phenotype_seed, .data$annotation_scale)

  # Expand sigma_0_2_scalars across runs with naive sigma strategy
  # Only applies to use cases with sigma_strategy == "naive"
  if (!is.null(sigma_0_2_scalars) && length(sigma_0_2_scalars) > 0) {
    naive_sigma_ids <- use_cases$use_case_id[use_cases$sigma_strategy == "naive"]
    
    # Split runs into those that need expansion and those that don't
    runs_to_expand <- tables$runs %>%
      dplyr::filter(.data$use_case_id %in% naive_sigma_ids)
    runs_no_expand <- tables$runs %>%
      dplyr::filter(!(.data$use_case_id %in% naive_sigma_ids)) %>%
      dplyr::mutate(sigma_0_2_scalar = NA_character_)
    
    if (nrow(runs_to_expand) > 0) {
      # Expand runs with each sigma_0_2_scalar value
      runs_expanded <- runs_to_expand %>%
        tidyr::crossing(sigma_0_2_scalar = as.character(sigma_0_2_scalars))
      
        tables$runs <- dplyr::bind_rows(runs_expanded, runs_no_expand) %>%
          dplyr::arrange(
            .data$dataset_bundle_id,
            .data$use_case_id,
            .data$L,
            .data$phenotype_seed,
            .data$annotation_scale,
            .data$sigma_0_2_scalar
          )
    } else {
      tables$runs <- runs_no_expand
    }
  } else {
    # No sigma_0_2_scalars specified - add column with NA
    tables$runs <- tables$runs %>%
      dplyr::mutate(sigma_0_2_scalar = NA_character_)
  }

  # Attach annotation seeds (unique per dataset bundle x annotation setting)
  runs <- tables$runs
  ann_key <- ifelse(
    is.na(runs$annotation_r2) & is.na(runs$inflate_match),
    NA_character_,
    paste0("r2=", runs$annotation_r2, "|inflate=", runs$inflate_match)
  )
  ann_lookup <- tibble::tibble(
    dataset_bundle_id = runs$dataset_bundle_id,
    ann_key = ann_key
  ) %>%
    dplyr::distinct() %>%
    dplyr::filter(!is.na(.data$ann_key)) %>%
    dplyr::mutate(annotation_seed = dplyr::row_number())
  runs <- runs %>%
    dplyr::mutate(ann_key = ann_key) %>%
    dplyr::left_join(ann_lookup, by = c("dataset_bundle_id", "ann_key")) %>%
    dplyr::select(-.data$ann_key)

  # Expand restart runs (one row per restart_id)
  restart_ids <- use_cases$use_case_id[use_cases$extra_compute == "restart_alpha"]
  if (length(restart_ids)) {
    n_inits <- restart_settings$n_inits %||% 5L
    restart_runs <- runs %>%
      dplyr::filter(.data$use_case_id %in% restart_ids) %>%
      tidyr::crossing(restart_id = seq_len(n_inits))
    non_restart_runs <- runs %>%
      dplyr::filter(!(.data$use_case_id %in% restart_ids)) %>%
      dplyr::mutate(restart_id = NA_integer_)
    runs <- dplyr::bind_rows(restart_runs, non_restart_runs)
  } else {
    runs <- runs %>% dplyr::mutate(restart_id = NA_integer_)
  }

  format_group_val <- function(x) {
    ifelse(is.na(x), "NA", format(x, digits = 6, scientific = FALSE))
  }

  runs <- runs %>%
    dplyr::mutate(
      group_key = paste0(
        "L=", as.integer(.data$L),
        "|r2=", format_group_val(.data$annotation_r2),
        "|inflate=", format_group_val(.data$inflate_match)
      )
    ) %>%
    dplyr::arrange(
      .data$dataset_bundle_id,
      .data$use_case_id,
      .data$L,
      .data$phenotype_seed,
      .data$annotation_scale,
      .data$sigma_0_2_scalar,
      .data$restart_id
    ) %>%
    dplyr::mutate(.row_id = dplyr::row_number())

  restart_lookup <- runs %>%
    dplyr::filter(!is.na(.data$restart_id)) %>%
    dplyr::mutate(restart_seed = dplyr::row_number()) %>%
    dplyr::select(.data$.row_id, .data$restart_seed)

  runs <- runs %>%
    dplyr::left_join(restart_lookup, by = ".row_id") %>%
    dplyr::mutate(run_id = dplyr::row_number()) %>%
    dplyr::select(-.data$.row_id)

  tables$runs <- runs

  shuffle_runs <- TRUE
  runs_tasks <- assign_task_ids_by_bundle(
    tables$runs,
    bundle_id_col = "dataset_bundle_id",
    bundles_per_task = bundles_per_task,
    shuffle_seed = task_assignment_seed,
    shuffle = shuffle_runs
  )

    compute_list <- list(
      model_average = model_average_settings,
      restart = restart_settings,
      aggregation_methods = aggregation_methods
    )
    if (!is.null(anneal_settings)) {
      compute_list$anneal <- anneal_settings
    }

    list(
      job = list(
      name = job_name,
      HPC = HPC,
      email = email,
      created_at = timestamp_utc(),
        task_unit = task_unit,
        bundles_per_task = bundles_per_task,
        runs_per_task = runs_per_task,
        credible_set_rho = credible_set_rho,
        purity_threshold = purity_threshold,
        verbose_file_output = isTRUE(verbose_file_output),
        write_legacy_snp_csv = isTRUE(write_legacy_snp_csv),
        write_snps_parquet = isTRUE(write_snps_parquet),
        write_confusion_bins = isTRUE(write_confusion_bins),
        buffer_flush_interval = as.integer(buffer_flush_interval),
        task_assignment_seed = task_assignment_seed,
          compute = compute_list,
        metrics = metrics_settings,
      slurm = list(
        time = time,
        mem = mem,
        cpus_per_task = 1,
        partition = NULL
      )
    ),
    paths = list(
      repo_root = repo_root,
      output_root = output_root,
      run_history_dir = file.path(output_root, "run_history", job_name),
      temp_dir = file.path(output_root, "temp", job_name),
      slurm_output_dir = file.path(output_root, "slurm_output"),
      slurm_prints_dir = file.path(output_root, "slurm_prints"),
      slurm_scripts_dir = file.path(output_root, "slurm_scripts")
    ),
    tables = list(
      dataset_bundles = tables$dataset_bundles,
      data_matrices = tables$data_matrices,
      runs = runs_tasks$runs,
      tasks = runs_tasks$tasks,
      use_cases = tables$use_cases
    )
  )
}
#' Write job configuration and scripts to disk.
#'
#' @param job_config Output of [make_job_config()].
#' @param run_task_script Path to the general task-runner script.
#' @return List with paths to the serialized artifacts.
#' @export
write_job_artifacts <- function(job_config,
                                run_task_script) {
  paths <- job_config$paths
  ensure_dir(paths$output_root)
  ensure_dir(paths$slurm_output_dir)
  ensure_dir(paths$slurm_prints_dir)
  ensure_dir(paths$slurm_scripts_dir)

  # Clear temp directory and create fresh
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

  dataset_path <- file.path(paths$temp_dir, "dataset_bundles.csv")
  readr::write_csv(job_config$tables$dataset_bundles, dataset_path)

  data_matrix_path <- file.path(paths$temp_dir, "data_matrices.csv")
  readr::write_csv(job_config$tables$data_matrices, data_matrix_path)

  use_case_path <- file.path(paths$temp_dir, "use_cases.csv")
  readr::write_csv(job_config$tables$use_cases, use_case_path)

  task_path <- file.path(paths$temp_dir, "task_table.csv")
  readr::write_csv(job_config$tables$tasks, task_path)

  slurm_path <- file.path(paths$slurm_scripts_dir, paste0(job_config$job$name, ".slurm"))
  writeLines(
    render_slurm_script(job_config, run_task_script = run_task_script),
    con = slurm_path
  )

  list(
    job_config = job_json_path,
    run_table = run_table_path,
    dataset_bundles = dataset_path,
    use_cases = use_case_path,
    tasks = task_path,
    slurm_script = slurm_path
  )
}

#' Render a SLURM submission script string.
#'
#' @param job_config Job configuration list.
#' @param run_task_script Path to `run_task.R`.
#' @return Character vector of script lines.
#' @keywords internal
render_slurm_script <- function(job_config, run_task_script) {
  job    <- job_config$job
  paths  <- job_config$paths
  tasks  <- job_config$tables$tasks
  n_tasks <- nrow(tasks)
  slurm  <- job$slurm

  partition_line <- if (!is.null(slurm$partition)) {
    paste0("#SBATCH --partition=", slurm$partition)
  } else {
    NULL
  }

  hpc_setup <- if (isTRUE(job$HPC)) {
    c(
      "module load r",
      "",
      'export R_LIBS_USER="/storage/home/mgc5166/R/x86_64-pc-linux-gnu-library/4.3"',
      ""
    )
  } else {
    NULL
  }

  script <- c(
    "#!/bin/bash",
    sprintf("#SBATCH --job-name=%s", job$name),
    sprintf("#SBATCH --array=1-%d", n_tasks),
    sprintf("#SBATCH --time=%s", slurm$time),
    sprintf("#SBATCH --mem=%s", slurm$mem),
    sprintf("#SBATCH --cpus-per-task=%s", slurm$cpus_per_task),
    sprintf("#SBATCH --mail-user=%s", job$email),
    "#SBATCH --mail-type=BEGIN,END,FAIL",
    partition_line,
    "# Log to /dev/null; we'll redirect ourselves below.",
    "#SBATCH --output=/dev/null",
    "#SBATCH --error=/dev/null",
    "",
    "set -euo pipefail",
    "",
    sprintf('JOB_ROOT="%s"', normalizePath(paths$output_root, winslash = "/", mustWork = FALSE)),
    sprintf('CONFIG_PATH="%s"', normalizePath(file.path(paths$temp_dir, "job_config.json"), winslash = "/", mustWork = FALSE)),
    sprintf('RUN_TASK_SCRIPT="%s"', normalizePath(run_task_script, winslash = "/", mustWork = FALSE)),
    sprintf('SLURM_PRINTS_BASE="%s"', normalizePath(paths$slurm_prints_dir, winslash = "/", mustWork = FALSE)),
    sprintf('SLURM_OUTPUT_BASE="%s"', normalizePath(paths$slurm_output_dir, winslash = "/", mustWork = FALSE)),
    sprintf('TEMP_DIR="%s"', normalizePath(paths$temp_dir, winslash = "/", mustWork = FALSE)),
    sprintf('RUN_HISTORY_BASE="%s"', normalizePath(file.path(paths$output_root, "run_history"), winslash = "/", mustWork = FALSE)),
    "",
    "# --- identifiers (stable across Slurm versions) ---",
    'JOBNAME="${SLURM_JOB_NAME}"',
    'PARENT_ID="${SLURM_ARRAY_JOB_ID:-$SLURM_JOB_ID}"   # parent array ID if array, else job id',
    'TASK_ID="${SLURM_ARRAY_TASK_ID:-0}"',
    'PRINTS_DIR="${SLURM_PRINTS_BASE}/${JOBNAME}/${PARENT_ID}"',
    'mkdir -p "${PRINTS_DIR}"',
    "",
    "# Redirect logs (after dirs exist) -- one file per task",
    'exec >"${PRINTS_DIR}/${TASK_ID}.out" 2>"${PRINTS_DIR}/${TASK_ID}.err"',
    'echo "[$(date -Is)] Starting task ${TASK_ID} for job ${JOBNAME} (parent ${PARENT_ID})"',
    "",
    "# --- export job info for R to compute run-specific output dirs ---",
    'export SUSINE_JOB_NAME="${JOBNAME}"',
    'export SUSINE_PARENT_ID="${PARENT_ID}"',
    "# --- site/module setup ---",
    hpc_setup,
    "",
    "# --- task 1 copies run_history ---",
    'if [ "${TASK_ID}" = "1" ]; then',
    '  FINAL_HISTORY_DIR="${RUN_HISTORY_BASE}/${JOBNAME}/${PARENT_ID}"',
    '  mkdir -p "${FINAL_HISTORY_DIR}"',
    '  cp "${TEMP_DIR}"/* "${FINAL_HISTORY_DIR}/"',
    '  echo "[$(date -Is)] Task 1: copied run_history from temp to ${FINAL_HISTORY_DIR}"',
    'fi',
    "",
    "# --- run task ---",
    'Rscript "$RUN_TASK_SCRIPT" \\',
    '  --job-name "$JOBNAME" \\',
    '  --task-id "$TASK_ID" \\',
    '  --job-root "$JOB_ROOT" \\',
    '  --config-path "$CONFIG_PATH"',
    "",
    'echo "[$(date -Is)] Completed task ${TASK_ID}"'
  )
  script[!is.na(script)]
}

#' Summarise expected run counts per use case and scenario.
#'
#' @param job_config Job configuration list.
#' @return Tibble summarising runs.
#' @export
summarise_job_config <- function(job_config) {
  runs <- job_config$tables$runs
  group_cols <- intersect(
    c(
      "data_scenario",
      "matrix_id",
      "use_case_id",
      "L",
      "y_noise",
      "p_star",
      "annotation_r2",
      "inflate_match",
        "gamma_shrink",
        "annotation_scale",
        "sigma_0_2_scalar",
        "restart_id"
      ),
    names(runs)
  )
  runs %>%
    dplyr::group_by(dplyr::across(dplyr::all_of(group_cols))) %>%
    dplyr::summarise(n_runs = dplyr::n(), .groups = "drop")
}
