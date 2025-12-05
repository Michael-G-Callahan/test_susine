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
#' @param seeds Integer vector of RNG seeds.
#' @param data_scenarios Character vector naming the data sources.
#' @param pair_L_p_star Logical; when TRUE (and `grid_mode == "full"`), the grid pairs matching
#'   `L` and `p_star` values instead of taking their Cartesian product.
#' @return List with elements `scenarios`, `runs`, and `tasks`.
#' @keywords internal
make_run_tables <- function(use_case_ids,
                            L_grid,
                            y_noise_grid,
                            prior_quality,
                            p_star_grid,
                            seeds,
                            data_scenarios,
                            grid_mode = c("full", "minimal"),
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
  unique_L <- unique(L_grid)
  unique_p <- unique(p_star_grid)

  build_full_grid <- function() {
    lp_pairs <- if (isTRUE(pair_L_p_star)) {
      if (!setequal(unique_L, unique_p)) {
        stop("When pair_L_p_star = TRUE, L_grid and p_star_grid must have the same unique values.")
      }
      tibble::tibble(L = unique_L, p_star = unique_L)
    } else {
      tidyr::expand_grid(L = unique_L, p_star = unique_p)
    }

    base_grid <- tidyr::expand_grid(
      data_scenario = unique(data_scenarios),
      y_noise = unique(y_noise_grid),
      prior_quality_id = prior_quality$prior_quality_id
    )

    base_grid %>%
      dplyr::mutate(.join_key = 1L) %>%
      dplyr::inner_join(
        lp_pairs %>% dplyr::mutate(.join_key = 1L),
        by = ".join_key"
      ) %>%
      dplyr::select(-.join_key) %>%
      dplyr::select(data_scenario, L, y_noise, p_star, prior_quality_id, dplyr::everything()) %>%
      dplyr::left_join(prior_quality, by = "prior_quality_id") %>%
      dplyr::arrange(data_scenario, L, y_noise, p_star, prior_quality_id) %>%
      dplyr::mutate(scenario_id = dplyr::row_number())
  }

  build_minimal_grid <- function() {
    ann_vals <- unique(prior_quality$annotation_r2)
    inflate_vals <- unique(prior_quality$inflate_match)
    gamma_vals <- unique(prior_quality$gamma_shrink)
    if (!length(ann_vals)) ann_vals <- NA_real_
    if (!length(inflate_vals)) inflate_vals <- NA_real_
    if (!length(gamma_vals)) gamma_vals <- NA_real_

    candidate_matrices <- data_matrix_catalog %>%
      dplyr::filter(.data$data_scenario %in% unique(data_scenarios)) %>%
      dplyr::arrange(.data$data_scenario, .data$matrix_id) %>%
      dplyr::pull(.data$matrix_id)
    if (!length(candidate_matrices)) {
      stop("No matrix metadata found for the requested scenarios.")
    }

    values <- list(
      data_scenario = unique(data_scenarios),
      L = unique(L_grid),
      y_noise = unique(y_noise_grid),
      p_star = unique(p_star_grid),
      annotation_r2 = ann_vals,
      inflate_match = inflate_vals,
      gamma_shrink = gamma_vals,
      matrix_id = unique(candidate_matrices)
    )
    lengths <- vapply(values, length, integer(1))
    n_rows <- max(c(lengths, length(candidate_matrices)))

    tibble::tibble(
      data_scenario = rep_len(values$data_scenario, n_rows),
      L = rep_len(values$L, n_rows),
      y_noise = rep_len(values$y_noise, n_rows),
      p_star = rep_len(values$p_star, n_rows),
      annotation_r2 = rep_len(values$annotation_r2, n_rows),
      inflate_match = rep_len(values$inflate_match, n_rows),
      gamma_shrink = rep_len(values$gamma_shrink, n_rows),
      matrix_id = rep_len(candidate_matrices, n_rows),
      prior_quality_id = seq_len(n_rows)
    ) %>%
      dplyr::mutate(scenario_id = seq_along(data_scenario))
  }

  scenarios <- if (grid_mode == "minimal") {
    build_minimal_grid()
  } else {
    build_full_grid()
  }

  matrix_catalog_subset <- dplyr::filter(data_matrix_catalog, .data$data_scenario %in% unique(data_scenarios))

  scenario_matrix_map <- if ("matrix_id" %in% names(scenarios)) {
    scenarios %>%
      dplyr::select(scenario_id, data_scenario, matrix_id) %>%
      dplyr::inner_join(
        matrix_catalog_subset,
        by = c("data_scenario", "matrix_id")
      )
  } else {
    scenarios %>%
      dplyr::select(scenario_id, data_scenario) %>%
      dplyr::inner_join(
        matrix_catalog_subset,
        by = "data_scenario"
      )
  }

  if (!nrow(scenario_matrix_map)) {
    stop("No matrix metadata found for the requested scenarios.")
  }

  seed_values <- unique(seeds)
  if (!length(seed_values)) {
    stop("At least one seed must be supplied.")
  }

  # Offset seeds by matrix_id so each dataset has its own unique seed set.

  # E.g., if seeds = 1:5, matrix_id 1 gets seeds 1-5, matrix_id 2 gets 6-10, etc.
  n_seeds <- length(seed_values)
  offset_seeds_for_matrix <- function(base_seeds, mat_id) {
    base_seeds + (mat_id - 1L) * n_seeds
  }

  seed_first <- seed_values[1]

  runs <- if (grid_mode == "minimal") {
    scenario_matrix_map %>%
      dplyr::mutate(
        use_case_id = rep_len(use_cases$use_case_id, dplyr::n()),
        seed = offset_seeds_for_matrix(seed_first, matrix_id)
      )
  } else {
    # Expand grid per scenario_matrix_index, then compute offset seeds per matrix_id
    tidyr::expand_grid(
      scenario_matrix_index = seq_len(nrow(scenario_matrix_map)),
      use_case_id = use_cases$use_case_id,
      seed_base = seed_values
    ) %>%
      dplyr::left_join(
        scenario_matrix_map %>%
          dplyr::mutate(scenario_matrix_index = dplyr::row_number()),
        by = "scenario_matrix_index"
      ) %>%
      dplyr::mutate(
        seed = offset_seeds_for_matrix(seed_base, matrix_id)
      ) %>%
      dplyr::select(-scenario_matrix_index, -seed_base)
  }

  scenario_join <- dplyr::select(
    scenarios,
    dplyr::all_of(setdiff(names(scenarios), c("data_scenario", "matrix_id")))
  )

  runs <- runs %>%
    dplyr::left_join(scenario_join, by = "scenario_id") %>%
    dplyr::left_join(
      dplyr::select(use_cases, use_case_id, requires_prior_quality, mu_strategy, sigma_strategy),
      by = "use_case_id"
    )

  runs <- runs %>%
    dplyr::mutate(
      needs_annotation = requires_prior_quality,
      needs_gamma = requires_prior_quality & sigma_strategy == "functional"
    ) %>%
    {\(df) {
      if (grid_mode == "minimal") {
        gamma_pool <- prior_quality$gamma_shrink
        gamma_pool <- gamma_pool[!is.na(gamma_pool)]
        gamma_default <- if (length(gamma_pool)) gamma_pool[1] else 0
        df %>%
          dplyr::mutate(
            annotation_r2 = dplyr::if_else(needs_annotation, annotation_r2, NA_real_),
            inflate_match = dplyr::if_else(needs_annotation, inflate_match, NA_real_),
            gamma_shrink = dplyr::case_when(
              needs_gamma ~ dplyr::coalesce(gamma_shrink, gamma_default),
              needs_annotation & !needs_gamma ~ NA_real_,
              TRUE ~ NA_real_
            )
          )
      } else {
        df
      }
    }}() %>%
    dplyr::filter(
      !needs_annotation |
        (needs_gamma & !is.na(gamma_shrink)) |
        (!needs_gamma & is.na(gamma_shrink))
    )

  default_scenarios <- scenarios %>%
    dplyr::group_by(data_scenario, L, y_noise, p_star) %>%
    dplyr::summarise(default_scenario_id = dplyr::first(scenario_id), .groups = "drop")

  runs <- runs %>%
    dplyr::left_join(default_scenarios,
      by = c("data_scenario", "L", "y_noise", "p_star")
    ) %>%
    dplyr::mutate(
      scenario_id = dplyr::if_else(
        requires_prior_quality,
        scenario_id,
        default_scenario_id
      ),
      prior_quality_id = dplyr::if_else(
        requires_prior_quality,
        prior_quality_id,
        as.integer(NA)
      ),
      annotation_r2 = dplyr::if_else(
        requires_prior_quality,
        annotation_r2,
        as.numeric(NA)
      ),
      inflate_match = dplyr::if_else(
        requires_prior_quality,
        inflate_match,
        as.numeric(NA)
      ),
      gamma_shrink = dplyr::if_else(
        needs_gamma,
        gamma_shrink,
        as.numeric(NA)
      )
    ) %>%
    dplyr::select(-default_scenario_id, -needs_annotation, -needs_gamma, -mu_strategy, -sigma_strategy) %>%
    dplyr::distinct() %>%
    dplyr::arrange(scenario_id, use_case_id, seed) %>%
    dplyr::mutate(run_id = dplyr::row_number())

  scenarios <- dplyr::semi_join(scenarios, runs, by = "scenario_id")

  list(
    scenarios = scenarios,
    data_matrices = data_matrix_catalog,
    runs = runs,
    use_cases = use_cases
  )
}

#' Attach task identifiers to the run table.
#'
#' @param runs Run tibble from [make_run_tables()].
#' @param runs_per_task Integer number of runs assigned to each SLURM task.
#' @param shuffle_seed Optional seed for reproducible shuffling (default: NULL for random).
#' @param shuffle Logical; when FALSE, assign runs to tasks in run_id order (no shuffling).
#' @return Tibble with `task_id` and `shuffled_order` columns added, and supporting task summary.
#' @keywords internal
assign_task_ids <- function(runs,
                            runs_per_task,
                            shuffle_seed = NULL,
                            shuffle = TRUE) {
  if (runs_per_task < 1) {
    stop("runs_per_task must be >= 1")
  }

  n_runs <- nrow(runs)

  if (isTRUE(shuffle)) {
    if (!is.null(shuffle_seed)) {
      set.seed(shuffle_seed)
    }
    shuffled_idx <- sample(n_runs)
    runs_shuffled <- runs %>%
      dplyr::mutate(shuffled_order = shuffled_idx) %>%
      dplyr::arrange(shuffled_order) %>%
      dplyr::mutate(task_id = ((dplyr::row_number() - 1L) %/% as.integer(runs_per_task)) + 1L) %>%
      dplyr::arrange(run_id)
  } else {
    runs_shuffled <- runs %>%
      dplyr::mutate(shuffled_order = seq_len(n_runs)) %>%
      dplyr::arrange(run_id) %>%
      dplyr::mutate(task_id = ((dplyr::row_number() - 1L) %/% as.integer(runs_per_task)) + 1L)
  }

  tasks_summary <- runs_shuffled %>%
    dplyr::group_by(task_id) %>%
    dplyr::summarise(runs_per_task = dplyr::n(), .groups = "drop") %>%
    dplyr::arrange(task_id)

  list(runs = runs_shuffled, tasks = tasks_summary)
}

#' Create an in-memory job configuration.
#'
#' @param job_name Character scalar.
#' @param use_case_ids Character vector of selected use cases.
#' @param L_grid Integer vector of L values.
#' @param y_noise_grid Numeric vector of noise fractions.
#' @param prior_quality Tibble with prior noise settings.
#' @param p_star_grid Integer vector.
#' @param seeds Integer vector of RNG seeds.
#' @param data_scenarios Character vector.
#' @param runs_per_task Integer number of runs assigned to each SLURM task.
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
#' @param anneal_settings Named list for tempering runs.
#' @param model_average_settings Named list for model averaging runs.
#' @param pair_L_p_star Logical; when TRUE, the full grid pairs `L` and `p_star` values 1-1
#'   instead of expanding their Cartesian product.
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
                            runs_per_task = 150,
                            email = "mgc5166@psu.edu",
                            output_root = "output",
                            credible_set_rho = 0.95,
                            purity_threshold = 0.5,
                            write_legacy_snp_csv = FALSE,
                            verbose_file_output = TRUE,
                            buffer_flush_interval = 300L,
                            grid_mode = c("full", "minimal"),
                            pair_L_p_star = FALSE,
                            sigma_0_2_scalars = NULL,
                            anneal_settings = list(
                              anneal_start_T = 5,
                              anneal_schedule_type = "geometric",
                              anneal_burn_in = 5
                            ),
                            model_average_settings = list(
                              n_inits = 5,
                              init_sd = 0.05
                            )) {

  grid_mode <- match.arg(grid_mode)
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

  # Expand sigma_0_2_scalars across runs with naive sigma strategy
  # Only applies to use cases with sigma_strategy == "naive"
  if (!is.null(sigma_0_2_scalars) && length(sigma_0_2_scalars) > 0) {
    use_cases <- tables$use_cases
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
        dplyr::arrange(.data$scenario_id, .data$use_case_id, .data$seed, .data$sigma_0_2_scalar) %>%
        dplyr::mutate(run_id = dplyr::row_number())
    } else {
      tables$runs <- runs_no_expand
    }
  } else {
    # No sigma_0_2_scalars specified - add column with NA
    tables$runs <- tables$runs %>%
      dplyr::mutate(sigma_0_2_scalar = NA_character_)
  }

  shuffle_runs <- TRUE
  runs_tasks <- assign_task_ids(
    tables$runs,
    runs_per_task = runs_per_task,
    shuffle_seed = NULL,
    shuffle = shuffle_runs
  )

  list(
    job = list(
      name = job_name,
      HPC = HPC,
      email = email,
      created_at = timestamp_utc(),
      runs_per_task = runs_per_task,
      credible_set_rho = credible_set_rho,
      purity_threshold = purity_threshold,
      verbose_file_output = isTRUE(verbose_file_output),
      write_legacy_snp_csv = isTRUE(write_legacy_snp_csv),
      buffer_flush_interval = as.integer(buffer_flush_interval),
      compute = list(
        anneal = anneal_settings,
        model_average = model_average_settings
      ),
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
      scenarios = tables$scenarios,
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

  scenario_path <- file.path(paths$temp_dir, "scenario_table.csv")
  readr::write_csv(job_config$tables$scenarios, scenario_path)

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
    scenarios = scenario_path,
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
      "sigma_0_2_scalar"
    ),
    names(runs)
  )
  runs %>%
    dplyr::group_by(dplyr::across(dplyr::all_of(group_cols))) %>%
    dplyr::summarise(n_runs = dplyr::n(), .groups = "drop")
}
