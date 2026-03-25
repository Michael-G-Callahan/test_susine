# Run-control builders ------------------------------------------------------

#' Build a tidy grid of prior-quality settings.
#'
#' @param annotation_r2_levels Numeric vector of target annotation R^2 values for causal SNPs.
#' @param inflate_match_levels Numeric vector indicating how strongly to inflate non-causal annotations.
#' @return Tibble with columns `annotation_r2`, `inflate_match`, and
#'   `prior_quality_id`.
#' @export
prior_quality_grid <- function(annotation_r2_levels = c(0.25, 0.5, 0.75),
                               inflate_match_levels = c(0, 0.5, 1)) {
  ann_vals <- unique(annotation_r2_levels)
  inflate_vals <- unique(inflate_match_levels)
  grid <- tidyr::expand_grid(
    annotation_r2 = ann_vals,
    inflate_match = inflate_vals
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
                            exploration_methods = c("single"),
                            exploration_mode = c("separate", "intersect"),
                            K = 1L,
                            L_grid,
                            y_noise_grid,
                            prior_quality,
                            p_star_grid,
                            seeds,
                            architecture_grid = c("sparse"),
                            data_scenarios,
                            grid_mode = c("full", "minimal", "intersect"),
                            pair_L_p_star = FALSE,
                            data_matrix_catalog = NULL,
                            c_grid_values = NULL,
                            tau_grid_values = NULL,
                            sigma_0_2_grid_values = NULL,
                            restart_settings = list(),
                            refine_settings = list()) {
  if (!is.data.frame(prior_quality) ||
      !all(c("annotation_r2", "inflate_match") %in% names(prior_quality))) {
    stop("prior_quality must contain annotation_r2 and inflate_match columns.")
  }
  if (is.null(data_matrix_catalog) || !nrow(data_matrix_catalog)) {
    stop("data_matrix_catalog must be supplied with one row per available matrix.")
  }
  grid_mode <- match.arg(grid_mode)
  exploration_mode <- match.arg(exploration_mode)
  K <- as.integer(K)
  if (!is.finite(K) || K < 1L) {
    stop("K must be a positive integer.")
  }
  exploration_methods <- unique(as.character(exploration_methods))
  valid_exploration_ids <- exploration_catalog()$exploration_id
  bad_methods <- setdiff(exploration_methods, valid_exploration_ids)
  if (length(bad_methods)) {
    stop("Unknown exploration method(s): ", paste(bad_methods, collapse = ", "))
  }

  prior_specs <- resolve_use_cases(use_case_ids)
  if (!nrow(prior_specs)) {
    stop("No valid prior specs selected.")
  }

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
    arch_vals <- unique(as.character(architecture_grid %||% "sparse"))
    candidate_matrices <- matrix_catalog_subset %>%
      dplyr::arrange(.data$data_scenario, .data$matrix_id) %>%
      dplyr::pull(.data$matrix_id)

    if (grid_mode == "minimal") {
      # Covering design: each setting value appears at least once via rep_len cycling.
      values <- list(
        y_noise = y_vals,
        p_star = p_vals,
        phenotype_seed = seed_values,
        architecture = arch_vals,
        matrix_id = candidate_matrices
      )
      lengths <- vapply(values, length, integer(1))
      n_rows <- max(lengths)
      bundles <- tibble::tibble(
        y_noise = rep_len(values$y_noise, n_rows),
        p_star = rep_len(values$p_star, n_rows),
        phenotype_seed = rep_len(values$phenotype_seed, n_rows),
        architecture = rep_len(values$architecture, n_rows),
        matrix_id = rep_len(values$matrix_id, n_rows)
      )
      return(bundles)
    }

    if (grid_mode == "intersect") {
      values <- list(
        y_noise = y_vals,
        p_star = p_vals,
        phenotype_seed = seed_values,
        architecture = arch_vals,
        matrix_id = candidate_matrices
      )
      lengths <- vapply(values, length, integer(1))
      n_rows <- max(lengths)
      bundles <- tibble::tibble(
        y_noise = rep_len(values$y_noise, n_rows),
        p_star = rep_len(values$p_star, n_rows),
        phenotype_seed = rep_len(values$phenotype_seed, n_rows),
        architecture = rep_len(values$architecture, n_rows),
        matrix_id = rep_len(values$matrix_id, n_rows)
      )
      return(bundles)
    }

    # Full grid
    bundles <- tidyr::expand_grid(
      matrix_id = candidate_matrices,
      y_noise = y_vals,
      p_star = p_vals,
      phenotype_seed = seed_values,
      architecture = arch_vals
    )
    bundles
  }

  dataset_bundles <- build_dataset_bundles()

  # h2-based architectures: y_noise and p_star are not meaningful parameters.
  # Set to NA and collapse duplicate rows to avoid redundant datasets.
  h2_archs <- c("susie2_sparse", "susie2_oligogenic")
  if (any(dataset_bundles$architecture %in% h2_archs)) {
    dataset_bundles <- dataset_bundles %>%
      dplyr::mutate(
        y_noise = dplyr::if_else(.data$architecture %in% h2_archs, NA_real_, .data$y_noise)
      )
    # susie2_oligogenic has fixed tier counts — p_star not applicable
    dataset_bundles <- dataset_bundles %>%
      dplyr::mutate(
        p_star = dplyr::if_else(.data$architecture == "susie2_oligogenic", NA_integer_, .data$p_star)
      ) %>%
      dplyr::distinct()
  }

  dataset_bundles <- dataset_bundles %>%
    dplyr::mutate(
      dataset_bundle_id = dplyr::row_number(),
      base_seed = phenotype_seed,
      phenotype_seed = as.integer(phenotype_seed * 10007L + dataset_bundle_id)
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

  annotation_grid_for_spec <- function(spec_row) {
    if (!isTRUE(spec_row$supports_annotation)) {
      return(tibble::tibble(annotation_r2 = NA_real_, inflate_match = NA_real_))
    }
    dplyr::distinct(dplyr::select(prior_quality, .data$annotation_r2, .data$inflate_match))
  }

  make_default_grid <- function(what, k) {
    if (what == "c_grid") {
      return(seq(0, 1.1, length.out = k))
    }
    if (what == "tau_grid") {
      return(seq(0.25, 3, length.out = k))
    }
    if (what == "sigma_0_2_grid") {
      vals <- 10^seq(-2, 0, length.out = k)
      # Snap nearest value to exactly 0.2 (required by experiment plan)
      nearest_idx <- which.min(abs(vals - 0.2))
      vals[nearest_idx] <- 0.2
      return(sort(unique(vals)))
    }
    stop("Unsupported default grid type: ", what)
  }

  axis_table_for_method <- function(method,
                                    mode,
                                    k,
                                    c_vals,
                                    tau_vals,
                                    sigma_0_2_vals = NULL,
                                    restart_n = NULL,
                                    refine_n = NULL,
                                    allow_infer_restart = FALSE,
                                    other_prod = 1L) {
    if (method == "single") {
      return(tibble::tibble())
    }
    if (method == "restart") {
      n_restart <- restart_n %||% NA_integer_
      if (!is.finite(n_restart)) {
        if (mode == "separate") {
          n_restart <- k
        } else if (isTRUE(allow_infer_restart)) {
          if (k %% as.integer(other_prod) != 0L) {
            stop("Cannot infer restart count: K is not divisible by non-restart axis product.")
          }
          n_restart <- as.integer(k / as.integer(other_prod))
        } else {
          stop("restart_settings$n_inits must be set for intersect mode when restart is included.")
        }
      }
      n_restart <- as.integer(n_restart)
      if (mode == "separate" && n_restart != as.integer(k)) {
        stop("For separate-mode restart exploration, restart_settings$n_inits must equal K.")
      }
      if (n_restart < 1L) {
        stop("restart count must be >= 1.")
      }
      return(tibble::tibble(
        restart_id = seq_len(n_restart),
        run_type = dplyr::if_else(seq_len(n_restart) == 1L, "default", "warm")
      ))
    }
    if (method == "c_grid") {
      vals <- c_vals
      if (is.null(vals) || !length(vals)) {
        vals <- make_default_grid("c_grid", if (mode == "separate") k else 1L)
      }
      vals <- as.numeric(vals)
      if (mode == "separate" && length(vals) != k) {
        stop("c_grid values must have length K for separate mode.")
      }
      return(tibble::tibble(c_value = vals))
    }
    if (method == "tau_grid") {
      vals <- tau_vals
      if (is.null(vals) || !length(vals)) {
        vals <- make_default_grid("tau_grid", if (mode == "separate") k else 1L)
      }
      vals <- as.numeric(vals)
      if (mode == "separate" && length(vals) != k) {
        stop("tau_grid values must have length K for separate mode.")
      }
      return(tibble::tibble(tau_value = vals))
    }
    if (method == "sigma_0_2_grid") {
      vals <- sigma_0_2_vals
      if (is.null(vals) || !length(vals)) {
        vals <- make_default_grid("sigma_0_2_grid", if (mode == "separate") k else 1L)
      }
      vals <- as.numeric(vals)
      if (mode == "separate" && length(vals) != k) {
        stop("sigma_0_2_grid values must have length K for separate mode.")
      }
      return(tibble::tibble(sigma_0_2_scalar = vals))
    }
    if (method == "refine") {
      n_refine <- refine_n %||% NA_integer_
      if (!is.finite(n_refine)) {
        if (mode == "separate") {
          n_refine <- k
        } else if (isTRUE(allow_infer_restart)) {
          if (k %% as.integer(other_prod) != 0L) {
            stop("Cannot infer refine count: K is not divisible by non-refine axis product.")
          }
          n_refine <- as.integer(k / as.integer(other_prod))
        } else {
          stop("refine_settings$n_steps must be set for intersect mode when refine is included.")
        }
      }
      n_refine <- as.integer(n_refine)
      if (mode == "separate" && n_refine != as.integer(k)) {
        stop("For separate-mode refine exploration, refine_settings$n_steps must equal K.")
      }
      if (n_refine < 1L) {
        stop("refine count must be >= 1.")
      }
      return(tibble::tibble(
        refine_step = seq_len(n_refine),
        run_type = rep("default", n_refine)
      ))
    }
    stop("Unsupported exploration method: ", method)
  }

  build_exploration_groups <- function(spec_row) {
    validity <- valid_exploration_for_prior(
      prior_spec_ids = spec_row$prior_spec_id,
      exploration_ids = exploration_methods
    )
    methods_for_spec <- validity$exploration_id[validity$valid]
    if (!length(methods_for_spec)) {
      stop("No valid exploration methods for prior spec: ", spec_row$prior_spec_id)
    }

    if (exploration_mode == "separate") {
      group_rows <- purrr::map_dfr(methods_for_spec, function(m) {
        axis_tbl <- axis_table_for_method(
          method = m,
          mode = "separate",
          k = K,
          c_vals = c_grid_values,
          tau_vals = tau_grid_values,
          sigma_0_2_vals = sigma_0_2_grid_values,
          restart_n = restart_settings$n_inits %||% NULL,
          refine_n = refine_settings$n_steps %||% NULL
        )
        if (!nrow(axis_tbl)) {
          axis_tbl <- tibble::tibble(.axis_stub = 1L)
        }
        axis_tbl %>%
          dplyr::mutate(
            exploration_mode = "separate",
            exploration_methods = m
          )
      })
      return(group_rows %>%
        dplyr::group_by(.data$exploration_mode, .data$exploration_methods) %>%
        dplyr::mutate(.exploration_row = dplyr::row_number()) %>%
        dplyr::ungroup())
    }

    # Intersect mode
    methods_intersect <- methods_for_spec
    non_restart_prod <- 1L
    if ("c_grid" %in% methods_intersect) {
      c_vals <- c_grid_values
      if (is.null(c_vals) || !length(c_vals)) {
        stop("c_grid_values must be provided for intersect mode when c_grid is used.")
      }
      non_restart_prod <- non_restart_prod * length(c_vals)
    }
    if ("tau_grid" %in% methods_intersect) {
      t_vals <- tau_grid_values
      if (is.null(t_vals) || !length(t_vals)) {
        stop("tau_grid_values must be provided for intersect mode when tau_grid is used.")
      }
      non_restart_prod <- non_restart_prod * length(t_vals)
    }
    if ("sigma_0_2_grid" %in% methods_intersect) {
      s_vals <- sigma_0_2_grid_values
      if (is.null(s_vals) || !length(s_vals)) {
        stop("sigma_0_2_grid_values must be provided for intersect mode when sigma_0_2_grid is used.")
      }
      non_restart_prod <- non_restart_prod * length(s_vals)
    }
    if ("single" %in% methods_intersect && length(methods_intersect) > 1L) {
      methods_intersect <- setdiff(methods_intersect, "single")
    }

    axis_list <- list()
    for (m in methods_intersect) {
      axis_list[[m]] <- axis_table_for_method(
        method = m,
        mode = "intersect",
        k = K,
        c_vals = c_grid_values,
        tau_vals = tau_grid_values,
        sigma_0_2_vals = sigma_0_2_grid_values,
        restart_n = restart_settings$n_inits %||% NULL,
        refine_n = refine_settings$n_steps %||% NULL,
        allow_infer_restart = TRUE,
        other_prod = non_restart_prod
      )
    }
    if (!length(axis_list)) {
      axis_tbl <- tibble::tibble(.axis_stub = 1L)
    } else {
      # Use unnamed inputs so crossing splices columns instead of creating list-cols.
      axis_tbl <- do.call(tidyr::crossing, unname(axis_list))
    }
    if (nrow(axis_tbl) != K) {
      stop(
        "Intersect exploration produced ", nrow(axis_tbl),
        " fits but K = ", K, ". Provide axis sizes whose product equals K."
      )
    }
    axis_tbl %>%
      dplyr::mutate(
        exploration_mode = "intersect",
        exploration_methods = paste(sort(methods_intersect), collapse = "x"),
        .exploration_row = dplyr::row_number()
      )
  }

  # Minimal-mode cycling state for annotation quality -------------------------
  if (grid_mode == "minimal") {
    ann_full_grid <- dplyr::distinct(
      dplyr::select(prior_quality, .data$annotation_r2, .data$inflate_match)
    )
    minimal_ann_idx <- 0L
  }

  runs <- purrr::map_dfr(seq_len(nrow(prior_specs)), function(i) {
    spec <- prior_specs[i, ]

    if (grid_mode == "minimal") {
      # 1 annotation combo per annotation-supporting spec, cycling through combos
      if (!isTRUE(spec$supports_annotation)) {
        ann_grid <- tibble::tibble(
          annotation_r2 = NA_real_, inflate_match = NA_real_
        )
      } else {
        minimal_ann_idx <<- minimal_ann_idx + 1L
        row_idx <- ((minimal_ann_idx - 1L) %% nrow(ann_full_grid)) + 1L
        ann_grid <- ann_full_grid[row_idx, , drop = FALSE]
      }
    } else {
      ann_grid <- annotation_grid_for_spec(spec)
    }

    exploration_grid <- build_exploration_groups(spec)

    if (grid_mode == "minimal") {
      # 1 bundle per use case, cycling through bundles
      bundle_idx <- ((i - 1L) %% nrow(dataset_bundles)) + 1L
      base <- dataset_bundles[bundle_idx, , drop = FALSE] %>%
        dplyr::inner_join(L_by_bundle, by = "dataset_bundle_id") %>%
        dplyr::mutate(
          prior_spec_id = spec$prior_spec_id,
          use_case_id = spec$prior_spec_id
        )
    } else {
      base <- dataset_bundles %>%
        dplyr::inner_join(L_by_bundle, by = "dataset_bundle_id") %>%
        dplyr::mutate(
          prior_spec_id = spec$prior_spec_id,
          use_case_id = spec$prior_spec_id
        )
    }
    tidyr::crossing(base, ann_grid, exploration_grid)
  }) %>%
    dplyr::mutate(
      group_key = paste0(
        "L=", as.integer(.data$L),
        "|r2=", ifelse(is.na(.data$annotation_r2), "NA", trimws(format(.data$annotation_r2, digits = 6, scientific = FALSE))),
        "|inflate=", ifelse(is.na(.data$inflate_match), "NA", trimws(format(.data$inflate_match, digits = 6, scientific = FALSE))),
        "|explore=", .data$exploration_mode, ":", .data$exploration_methods
      )
    ) %>%
    dplyr::arrange(
      .data$dataset_bundle_id,
      .data$prior_spec_id,
      .data$L,
      .data$phenotype_seed,
      .data$group_key,
      .data$.exploration_row
    ) %>%
    dplyr::mutate(run_id = dplyr::row_number()) %>%
    dplyr::select(-dplyr::any_of(c(".axis_stub", ".exploration_row")))

  list(
    dataset_bundles = dataset_bundles,
    data_matrices = data_matrix_catalog,
    runs = runs,
    use_cases = prior_specs
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
    n_bundle <- nrow(bundles)
    bundles <- bundles %>%
      dplyr::mutate(shuffled_order = sample.int(n_bundle)) %>%
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
#' @param use_case_ids Character vector of selected prior spec ids (or aliases).
#' @param exploration_methods Character vector of exploration ids.
#' @param exploration_mode One of `separate` or `intersect`.
#' @param K Integer fit budget per dataset x use_case x exploration_group.
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
#'   scalars to test for fixed-variance prior specs.
#'   Supported formats:
#'   - Numeric values (e.g., 0.1, 0.2, 0.4): multiplier of var(y)
#'   - Strings with /L suffix (e.g., "1/L", "0.5/L", "2/L"): multiplier/L of var(y)
#'   When NULL, defaults to `0.2`.
#' @param c_grid_values Numeric vector of c values for c-grid exploration.
#' @param tau_grid_values Numeric vector of tau values for tau-grid exploration.
#' @param restart_settings List controlling restart exploration.
#' @param aggregation_methods Character vector of per-group aggregation methods.
#' @param overall_aggregation_methods Character vector of global-pool aggregation methods.
#' @param include_overall_pool Logical; include full per-dataset global pool aggregation.
#' @param max_iter Maximum number of IBSS iterations for susie/susine fits
#'   (default 100). Lower values speed up runs at the cost of convergence.
#' @param softmax_temperature Numeric temperature for ELBO softmax weighting.
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
                            exploration_methods = c("single"),
                            exploration_mode = c("separate", "intersect"),
                            K = 1L,
                            L_grid,
                            y_noise_grid,
                            prior_quality,
                            p_star_grid,
                            seeds,
                            architecture_grid = c("sparse"),
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
                            write_tier_cs_metrics = FALSE,
                            write_prior_diagnostics = FALSE,
                            verbose_file_output = TRUE,
                            buffer_flush_interval = 300L,
                            task_assignment_seed = NULL,
                            grid_mode = c("full", "minimal", "intersect"),
                            pair_L_p_star = FALSE,
                            sigma_0_2_scalars = NULL,
                            c_grid_values = NULL,
                            tau_grid_values = NULL,
                            sigma_0_2_grid_values = NULL,
                            restart_settings = list(
                              n_inits = NULL,
                              alpha_concentration = 1
                            ),
                            alpha_concentration_grid = NULL,
                            warm_start_grid = NULL,
                            refine_settings = list(
                              n_steps = NULL,
                              cs_source = "filtered"
                            ),
                            max_iter = 100L,
                            aggregation_methods = c("elbo_softmax"),
                            overall_aggregation_methods = NULL,
                            exploration_specs = NULL,
                            include_overall_pool = TRUE,
                            softmax_temperature = 1,
                            metrics_settings = list(
                              pip_bucket_width = 0.01,
                              z_top_k = 10,
                              jsd_threshold = 0.15
                            )) {

  grid_mode <- match.arg(grid_mode)
  exploration_mode <- match.arg(exploration_mode)
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

  # Build run tables: either from a single global config or from

  # heterogeneous exploration_specs (list of per-spec overrides).
  if (!is.null(exploration_specs) && length(exploration_specs) > 0) {
    # Heterogeneous mode: one make_run_tables call per spec, then merge.
    all_runs <- list()
    shared_tables <- NULL
    for (si in seq_along(exploration_specs)) {
      spec <- exploration_specs[[si]]
      spec_name <- spec$name %||% paste0("spec_", si)
      spec_tables <- make_run_tables(
        use_case_ids = spec$use_case_ids %||% use_case_ids,
        exploration_methods = spec$exploration_methods %||% exploration_methods,
        exploration_mode = spec$exploration_mode %||% exploration_mode,
        K = spec$K %||% K,
        L_grid = L_grid,
        y_noise_grid = y_noise_grid,
        prior_quality = prior_quality,
        p_star_grid = p_star_grid,
        seeds = seeds,
        architecture_grid = architecture_grid,
        data_scenarios = data_scenarios,
        grid_mode = grid_mode,
        pair_L_p_star = pair_L_p_star,
        data_matrix_catalog = data_matrix_catalog,
        c_grid_values = spec$c_grid_values %||% c_grid_values,
        tau_grid_values = spec$tau_grid_values %||% tau_grid_values,
        sigma_0_2_grid_values = spec$sigma_0_2_grid_values %||% sigma_0_2_grid_values,
        restart_settings = spec$restart_settings %||% restart_settings,
        refine_settings = spec$refine_settings %||% refine_settings
      )
      spec_runs <- spec_tables$runs %>%
        dplyr::mutate(spec_name = spec_name)
      all_runs[[si]] <- spec_runs
      if (is.null(shared_tables)) {
        shared_tables <- spec_tables
      } else {
        # Union use_cases across specs
        shared_tables$use_cases <- dplyr::bind_rows(
          shared_tables$use_cases, spec_tables$use_cases
        ) %>% dplyr::distinct(.data$prior_spec_id, .keep_all = TRUE)
      }
    }
    tables <- shared_tables
    tables$runs <- dplyr::bind_rows(all_runs)
  } else {
    # Single-config mode (original behavior)
    tables <- make_run_tables(
      use_case_ids = use_case_ids,
      exploration_methods = exploration_methods,
      exploration_mode = exploration_mode,
      K = K,
      L_grid = L_grid,
      y_noise_grid = y_noise_grid,
      prior_quality = prior_quality,
      p_star_grid = p_star_grid,
      seeds = seeds,
      architecture_grid = architecture_grid,
      data_scenarios = data_scenarios,
      grid_mode = grid_mode,
      pair_L_p_star = pair_L_p_star,
      data_matrix_catalog = data_matrix_catalog,
      c_grid_values = c_grid_values,
      tau_grid_values = tau_grid_values,
      sigma_0_2_grid_values = sigma_0_2_grid_values,
      restart_settings = restart_settings,
      refine_settings = refine_settings
    )
  }

  prior_specs <- tables$use_cases
  runs <- tables$runs

  if (is.null(sigma_0_2_scalars) || !length(sigma_0_2_scalars)) {
    sigma_0_2_scalars <- 0.2
  }

  # Specs that need sigma_0_2_scalar expansion:
  # 1. Fixed-variance specs (prior_variance_strategy == "fixed")
  # 2. "scale" EB method specs (sigma fixed, c optimized — sigma_0_2_scalar sets the fixed sigma)
  # 3. "mean" EB method specs (sigma fixed, scalar mu optimized)
  needs_sigma_ids <- prior_specs$prior_spec_id[
    prior_specs$prior_variance_strategy == "fixed" |
      (!is.na(prior_specs$eb_method) & prior_specs$eb_method %in% c("scale", "mean"))
  ]
  runs_fixed <- runs %>% dplyr::filter(.data$prior_spec_id %in% needs_sigma_ids)
  runs_non_fixed <- runs %>%
    dplyr::filter(!(.data$prior_spec_id %in% needs_sigma_ids)) %>%
    dplyr::mutate(sigma_0_2_scalar = NA_character_)

  # Runs from sigma_0_2_grid exploration already have sigma_0_2_scalar set;

  # only expand sigma_0_2_scalars for runs without a pre-set value.
  has_sigma <- "sigma_0_2_scalar" %in% names(runs_fixed) &
    !all(is.na(runs_fixed$sigma_0_2_scalar))
  if (has_sigma) {
    runs_fixed_preset <- runs_fixed %>%
      dplyr::filter(!is.na(.data$sigma_0_2_scalar))
    runs_fixed_needs <- runs_fixed %>%
      dplyr::filter(is.na(.data$sigma_0_2_scalar)) %>%
      dplyr::select(-.data$sigma_0_2_scalar)
  } else {
    runs_fixed_preset <- runs_fixed[0, , drop = FALSE]
    runs_fixed_needs <- runs_fixed
  }
  if (nrow(runs_fixed_needs)) {
    if (grid_mode == "minimal") {
      sigma_vals <- as.character(sigma_0_2_scalars)
      spec_sigma_map <- runs_fixed_needs %>%
        dplyr::distinct(.data$prior_spec_id) %>%
        dplyr::mutate(
          .sigma_idx = ((dplyr::row_number() - 1L) %% length(sigma_vals)) + 1L,
          sigma_0_2_scalar = sigma_vals[.sigma_idx]
        ) %>%
        dplyr::select(.data$prior_spec_id, .data$sigma_0_2_scalar)
      runs_fixed_needs <- runs_fixed_needs %>%
        dplyr::inner_join(spec_sigma_map, by = "prior_spec_id")
    } else {
      runs_fixed_needs <- runs_fixed_needs %>%
        tidyr::crossing(sigma_0_2_scalar = as.character(sigma_0_2_scalars))
    }
  } else {
    runs_fixed_needs <- dplyr::mutate(runs_fixed_needs, sigma_0_2_scalar = character(0))
  }
  runs_fixed_preset <- dplyr::mutate(runs_fixed_preset,
                                      sigma_0_2_scalar = as.character(.data$sigma_0_2_scalar))
  runs <- dplyr::bind_rows(runs_fixed_preset, runs_fixed_needs, runs_non_fixed)

  # Warm-start grid expansion: cross restart runs with (warm_method, alpha) pairs
  ws_grid <- warm_start_grid
  if (is.null(ws_grid) && !is.null(alpha_concentration_grid) &&
      length(alpha_concentration_grid) > 0) {
    # Backward compat: old-style alpha-only grid defaults to init_alpha method
    ws_grid <- tibble::tibble(
      warm_method = "init_alpha",
      alpha_concentration = as.numeric(alpha_concentration_grid)
    )
  }
  if (!is.null(ws_grid) && nrow(ws_grid) > 0) {
    if (!"restart_id" %in% names(runs)) runs$restart_id <- NA_integer_
    has_restart <- !is.na(runs$restart_id)
    runs_restart <- runs %>% dplyr::filter(has_restart)
    runs_other <- runs %>%
      dplyr::filter(!has_restart) %>%
      dplyr::mutate(warm_method = NA_character_, alpha_concentration = NA_real_)
    if (nrow(runs_restart)) {
      runs_restart <- runs_restart %>%
        tidyr::crossing(ws_grid)
      # Update group_key to include method + alpha for separate aggregation
      runs_restart <- runs_restart %>%
        dplyr::mutate(
          group_key = paste0(.data$group_key,
                             "|wm=", .data$warm_method,
                             "|alpha=", .data$alpha_concentration)
        )
    } else {
      runs_restart <- dplyr::mutate(runs_restart,
                                    warm_method = character(0),
                                    alpha_concentration = numeric(0))
    }
    runs <- dplyr::bind_rows(runs_restart, runs_other)
  } else {
    runs$warm_method <- NA_character_
    runs$alpha_concentration <- NA_real_
  }

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

  if (!"restart_id" %in% names(runs)) {
    runs$restart_id <- NA_integer_
  }
  if (!"run_type" %in% names(runs)) {
    runs$run_type <- if ("init_type" %in% names(runs)) runs$init_type else NA_character_
  }
  if (!"c_value" %in% names(runs)) {
    runs$c_value <- NA_real_
  }
  if (!"tau_value" %in% names(runs)) {
    runs$tau_value <- NA_real_
  }
  if (!"refine_step" %in% names(runs)) {
    runs$refine_step <- NA_integer_
  }
  if (!"architecture" %in% names(runs)) {
    runs$architecture <- "sparse"
  }

  runs <- runs %>% dplyr::mutate(.run_key = dplyr::row_number())
  restart_lookup <- runs %>%
    dplyr::filter(!is.na(.data$restart_id)) %>%
    dplyr::mutate(restart_seed = dplyr::row_number()) %>%
    dplyr::select(.data$.run_key, .data$restart_seed)

  runs <- runs %>%
    dplyr::left_join(restart_lookup, by = ".run_key") %>%
    dplyr::arrange(
      .data$dataset_bundle_id,
      .data$prior_spec_id,
      .data$L,
      .data$phenotype_seed,
      .data$group_key,
      dplyr::coalesce(.data$restart_id, 0L),
      dplyr::coalesce(.data$refine_step, 0L),
      .data$c_value,
      .data$tau_value
    ) %>%
    dplyr::mutate(run_id = dplyr::row_number()) %>%
    dplyr::select(-.data$.run_key)

  tables$runs <- runs

  # Minimal mode: force all bundles into a single task
  effective_bundles_per_task <- if (grid_mode == "minimal") {
    nrow(tables$dataset_bundles)
  } else {
    bundles_per_task
  }

  runs_tasks <- assign_task_ids_by_bundle(
    tables$runs,
    bundle_id_col = "dataset_bundle_id",
    bundles_per_task = effective_bundles_per_task,
    shuffle_seed = task_assignment_seed,
    shuffle = TRUE
  )

  if (is.null(overall_aggregation_methods)) {
    overall_aggregation_methods <- aggregation_methods
  }
  compute_list <- list(
    restart = restart_settings,
    refine = refine_settings,
    K = as.integer(K),
    exploration_methods = exploration_methods,
    exploration_mode = exploration_mode,
    exploration_specs = exploration_specs,
    aggregation_methods = aggregation_methods,
    overall_aggregation_methods = overall_aggregation_methods,
    include_overall_pool = isTRUE(include_overall_pool),
    softmax_temperature = as.numeric(softmax_temperature)
  )

  list(
    job = list(
      name = job_name,
      HPC = HPC,
      email = email,
      created_at = timestamp_utc(),
      task_unit = task_unit,
      n_tasks = nrow(runs_tasks$tasks),
      bundles_per_task = bundles_per_task,
      runs_per_task = runs_per_task,
      credible_set_rho = credible_set_rho,
      purity_threshold = purity_threshold,
      max_iter = as.integer(max_iter),
      verbose_file_output = isTRUE(verbose_file_output),
      write_legacy_snp_csv = isTRUE(write_legacy_snp_csv),
      write_snps_parquet = isTRUE(write_snps_parquet),
      write_confusion_bins = isTRUE(write_confusion_bins),
      write_tier_cs_metrics = isTRUE(write_tier_cs_metrics),
      write_prior_diagnostics = isTRUE(write_prior_diagnostics),
      buffer_flush_interval = as.integer(buffer_flush_interval),
      task_assignment_seed = task_assignment_seed,
      compute = compute_list,
      metrics = metrics_settings,
      slurm = list(
        time = time,
        mem = mem,
        cpus_per_task = 1,
        partition = NULL,
        account = NULL
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
      use_cases = prior_specs
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

  account_line <- if (!is.null(slurm$account)) {
    paste0("#SBATCH --account=", slurm$account)
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
    account_line,
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
    'export SUSINE_DEV="1"',
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
      "prior_spec_id",
      "L",
      "y_noise",
      "p_star",
      "architecture",
      "annotation_r2",
      "inflate_match",
      "sigma_0_2_scalar",
      "restart_id",
      "refine_step",
      "run_type",
      "c_value",
      "tau_value",
      "exploration_mode",
      "exploration_methods",
      "group_key"
    ),
    names(runs)
  )
  runs %>%
    dplyr::group_by(dplyr::across(dplyr::all_of(group_cols))) %>%
    dplyr::summarise(n_runs = dplyr::n(), .groups = "drop")
}
