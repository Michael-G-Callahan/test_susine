# Shared config builder for the 5-arm effect-matching refit study.
# Sourced by both effect_matching_refit_export_5arm.Rmd (cold + warm fits)
# and effect_matching_drift_analysis_5arm.Rmd (PR curve regeneration of
# causal_idx via generate_data_for_bundle()), so the dataset bundles, run
# rows, and use_case grid stay aligned across the two workbooks.
#
# Caller is expected to have already loaded test_susine via
# devtools::load_all(), and to have dplyr / tidyr / readr / tibble available.

study_name_5arm <- "baseline_sims_5arm_effect_matching"

# Genotype / phenotype grid (must match across both workbooks so
# dataset_bundle_id keys align between effect_refit_records.rds and
# pip_records_5arm.rds).
n_matrices_target_5arm <- 150L
seeds_5arm             <- 1:4
L_value_5arm           <- 10L
y_noise_value_5arm     <- NA_real_
p_star_value_5arm      <- NA_integer_
architecture_value_5arm <- "susie2_oligogenic"
data_scenario_5arm     <- "scenario_1"

# Baseline screening result used to lock the 5-arm treatment settings.
baseline_job_name_5arm <- "baseline_sims_screen"
baseline_parent_job_id_5arm <- "51509956"

# Treatment-arm annotation regime. The scalar settings below are only a
# manual fallback for exploratory use; production 5-arm reruns should select
# from the locked baseline screening_summary.csv.
annotation_r2_value_5arm <- 0.5
inflate_match_value_5arm <- 0.9
c_value_mu_5arm <- 0.6
tau_value_pi_5arm <- 0.464
sigma_0_2_vanilla_value_5arm <- 0.2
sigma_0_2_mu_value_5arm <- 0.01
sigma_0_2_pi_value_5arm <- 0.2
max_iter_5arm        <- 100L
credible_set_rho_5arm <- 0.95
purity_threshold_5arm <- 0.50

cold_method_levels_5arm <- c(
  "susie_vanilla", "susine_functional_mu", "susie_functional_pi"
)
warm_method_levels_5arm <- c(
  "susie_vanilla_warm_from_mu", "susie_vanilla_warm_from_pi"
)
all_method_levels_5arm <- c(cold_method_levels_5arm, warm_method_levels_5arm)

fmt_num_5arm <- function(x, digits = 3) {
  format(signif(as.numeric(x), digits), trim = TRUE, scientific = FALSE)
}

annotation_label_5arm <- function(annotation_r2, inflate_match) {
  ifelse(
    is.na(annotation_r2) | is.na(inflate_match),
    "no_annotation",
    paste0("r2=", fmt_num_5arm(annotation_r2, 2),
           " | inflate=", fmt_num_5arm(inflate_match, 2))
  )
}

settings_label_5arm <- function(method_family, c_value, tau_value,
                                sigma_0_2_scalar) {
  dplyr::case_when(
    method_family == "susie_vanilla" ~
      paste0("susie_vanilla | sigma=", sigma_0_2_scalar),
    method_family == "susine_functional_mu" ~
      paste0("susine_functional_mu | c=", fmt_num_5arm(c_value, 3),
             " | sigma=", sigma_0_2_scalar),
    method_family == "susie_functional_pi" ~
      paste0("susie_functional_pi | tau=", fmt_num_5arm(tau_value, 3),
             " | sigma=", sigma_0_2_scalar),
    TRUE ~ method_family
  )
}

manual_5arm_settings <- function() {
  tibble::tibble(
    method_family = cold_method_levels_5arm,
    annotation_r2 = c(NA_real_, annotation_r2_value_5arm,
                      annotation_r2_value_5arm),
    inflate_match = c(NA_real_, inflate_match_value_5arm,
                      inflate_match_value_5arm),
    c_value = c(NA_real_, c_value_mu_5arm, NA_real_),
    tau_value = c(NA_real_, NA_real_, tau_value_pi_5arm),
    sigma_0_2_scalar = as.character(c(
      sigma_0_2_vanilla_value_5arm,
      sigma_0_2_mu_value_5arm,
      sigma_0_2_pi_value_5arm
    )),
    AUPRC = NA_real_,
    setting_source = "manual_fallback"
  ) %>%
    dplyr::mutate(
      annotation_label = annotation_label_5arm(.data$annotation_r2,
                                               .data$inflate_match),
      setting_label = settings_label_5arm(
        .data$method_family, .data$c_value, .data$tau_value,
        .data$sigma_0_2_scalar
      )
    )
}

pick_best_5arm_settings <- function(screening_summary,
                                    annotation_r2 = annotation_r2_value_5arm,
                                    inflate_match = inflate_match_value_5arm) {
  required_cols <- c(
    "method_family", "annotation_r2", "inflate_match",
    "sigma_0_2_scalar", "c_value", "tau_value", "AUPRC"
  )
  missing_cols <- setdiff(required_cols, names(screening_summary))
  if (length(missing_cols)) {
    stop("Baseline screening summary is missing required columns: ",
         paste(missing_cols, collapse = ", "))
  }

  pick_one <- function(method, annotated) {
    rows <- screening_summary %>%
      dplyr::filter(.data$method_family == !!method)
    if (annotated) {
      rows <- rows %>%
        dplyr::filter(
          abs(as.numeric(.data$annotation_r2) - .env$annotation_r2) < 1e-8,
          abs(as.numeric(.data$inflate_match) - .env$inflate_match) < 1e-8
        )
    }
    rows <- rows %>%
      dplyr::filter(is.finite(.data$AUPRC)) %>%
      dplyr::arrange(dplyr::desc(.data$AUPRC))
    if (!nrow(rows)) {
      stop("No baseline-screening row found for ", method,
           if (annotated) {
             paste0(" at annotation_r2=", annotation_r2,
                    ", inflate_match=", inflate_match)
           } else {
             ""
           })
    }
    rows[1, , drop = FALSE]
  }

  dplyr::bind_rows(
    pick_one("susie_vanilla", annotated = FALSE),
    pick_one("susine_functional_mu", annotated = TRUE),
    pick_one("susie_functional_pi", annotated = TRUE)
  ) %>%
    dplyr::mutate(
      method_family = as.character(.data$method_family),
      annotation_r2 = as.numeric(.data$annotation_r2),
      inflate_match = as.numeric(.data$inflate_match),
      c_value = as.numeric(.data$c_value),
      tau_value = as.numeric(.data$tau_value),
      sigma_0_2_scalar = as.character(.data$sigma_0_2_scalar),
      annotation_label = annotation_label_5arm(.data$annotation_r2,
                                               .data$inflate_match),
      setting_label = settings_label_5arm(
        .data$method_family, .data$c_value, .data$tau_value,
        .data$sigma_0_2_scalar
      ),
      setting_source = "baseline_screening_summary"
    ) %>%
    dplyr::select(dplyr::all_of(c(
      "method_family", "setting_label", "annotation_label",
      "annotation_r2", "inflate_match", "sigma_0_2_scalar",
      "c_value", "tau_value", "AUPRC", "setting_source"
    )))
}

resolve_5arm_settings <- function(output_root,
                                  settings = NULL,
                                  baseline_screening_csv = NULL,
                                  settings_cache_csv = NULL,
                                  prefer_settings_cache = TRUE,
                                  allow_manual_fallback = FALSE) {
  if (!is.null(settings)) {
    return(settings)
  }

  if (isTRUE(prefer_settings_cache) &&
      !is.null(settings_cache_csv) &&
      file.exists(settings_cache_csv)) {
    return(readr::read_csv(settings_cache_csv, show_col_types = FALSE))
  }

  if (is.null(baseline_screening_csv)) {
    baseline_screening_csv <- file.path(
      output_root, "slurm_output", baseline_job_name_5arm,
      baseline_parent_job_id_5arm, "consolidated", "screening_summary.csv"
    )
  }
  if (file.exists(baseline_screening_csv)) {
    out <- pick_best_5arm_settings(
      readr::read_csv(baseline_screening_csv, show_col_types = FALSE)
    )
    if (!is.null(settings_cache_csv)) {
      dir.create(dirname(settings_cache_csv), recursive = TRUE,
                 showWarnings = FALSE)
      readr::write_csv(out, settings_cache_csv)
    }
    return(out)
  }

  if (!isTRUE(prefer_settings_cache) &&
      !is.null(settings_cache_csv) &&
      file.exists(settings_cache_csv)) {
    return(readr::read_csv(settings_cache_csv, show_col_types = FALSE))
  }

  if (isTRUE(allow_manual_fallback)) {
    warning(
      "Baseline screening summary not found; using manual fallback 5-arm ",
      "settings. This is not recommended for production reruns.\nMissing: ",
      baseline_screening_csv
    )
    return(manual_5arm_settings())
  }

  stop(
    "Cannot resolve 5-arm settings. Expected locked baseline summary at:\n",
    baseline_screening_csv,
    "\nRun collect_results_workbook_baseline_sims.Rmd first, or provide ",
    "settings_cache_csv/settings explicitly."
  )
}

filter_baseline_runs_for_5arm <- function(base_runs, settings) {
  required_cols <- c(
    "dataset_bundle_id", "use_case_id", "annotation_r2", "inflate_match",
    "sigma_0_2_scalar", "c_value", "tau_value", "annotation_seed"
  )
  missing_cols <- setdiff(required_cols, names(base_runs))
  if (length(missing_cols)) {
    stop("Baseline run table is missing required columns: ",
         paste(missing_cols, collapse = ", "))
  }

  pick_method_runs <- function(method) {
    setting <- settings %>%
      dplyr::filter(.data$method_family == !!method)
    if (nrow(setting) != 1L) {
      stop("Expected exactly one selected setting for ", method,
           "; found ", nrow(setting), ".")
    }

    rows <- base_runs %>%
      dplyr::filter(.data$use_case_id == !!method) %>%
      dplyr::filter(
        abs(as.numeric(.data$sigma_0_2_scalar) -
              as.numeric(setting$sigma_0_2_scalar[[1]])) < 1e-8
      )

    if (method == "susine_functional_mu") {
      rows <- rows %>%
        dplyr::filter(
          abs(as.numeric(.data$annotation_r2) -
                as.numeric(setting$annotation_r2[[1]])) < 1e-8,
          abs(as.numeric(.data$inflate_match) -
                as.numeric(setting$inflate_match[[1]])) < 1e-8,
          abs(as.numeric(.data$c_value) -
                as.numeric(setting$c_value[[1]])) < 1e-6
        )
    } else if (method == "susie_functional_pi") {
      rows <- rows %>%
        dplyr::filter(
          abs(as.numeric(.data$annotation_r2) -
                as.numeric(setting$annotation_r2[[1]])) < 1e-8,
          abs(as.numeric(.data$inflate_match) -
                as.numeric(setting$inflate_match[[1]])) < 1e-8,
          abs(as.numeric(.data$tau_value) -
                as.numeric(setting$tau_value[[1]])) < 1e-6
        )
    }

    rows <- rows %>%
      dplyr::mutate(
        method_family = .data$use_case_id,
        annotation_label = annotation_label_5arm(.data$annotation_r2,
                                                 .data$inflate_match),
        setting_label = settings_label_5arm(
          .data$use_case_id, .data$c_value, .data$tau_value,
          .data$sigma_0_2_scalar
        )
      )

    if (!nrow(rows)) {
      stop("No baseline run rows matched selected 5-arm setting for ", method,
           ".")
    }
    rows
  }

  selected_runs <- dplyr::bind_rows(lapply(
    cold_method_levels_5arm, pick_method_runs
  )) %>%
    dplyr::mutate(dataset_bundle_id = as.character(.data$dataset_bundle_id)) %>%
    dplyr::arrange(.data$dataset_bundle_id, .data$use_case_id)

  expected_n <- dplyr::n_distinct(selected_runs$dataset_bundle_id) *
    length(cold_method_levels_5arm)
  if (nrow(selected_runs) != expected_n) {
    stop(
      "Selected baseline run rows are not one row per dataset per cold method. ",
      "Rows=", nrow(selected_runs), "; expected=", expected_n, "."
    )
  }

  seed_check <- selected_runs %>%
    dplyr::filter(.data$use_case_id %in%
                    c("susine_functional_mu", "susie_functional_pi")) %>%
    dplyr::count(.data$dataset_bundle_id, .data$use_case_id,
                 name = "n_rows") %>%
    dplyr::filter(.data$n_rows != 1L)
  if (nrow(seed_check)) {
    stop("Selected baseline annotation rows are not unique by dataset/method.")
  }

  selected_runs
}

# Build the job_config + selected run rows for the 5-arm refit grid.
#
# Args:
#   here_root: project root (output of `here::here()`).
#   output_root: directory under which the study output lives (one level
#     above study_dir).
#   study_name: study identifier; defaults to study_name_5arm.
#
# Returns a list with cfg, selected_runs, data_matrix_catalog.
build_5arm_cfg <- function(here_root, output_root,
                           study_name = study_name_5arm,
                           settings = NULL,
                           baseline_screening_csv = NULL,
                           baseline_cfg_path = NULL,
                           settings_cache_csv = NULL,
                           prefer_settings_cache = TRUE,
                           allow_manual_fallback = FALSE) {
  settings <- resolve_5arm_settings(
    output_root = output_root,
    settings = settings,
    baseline_screening_csv = baseline_screening_csv,
    settings_cache_csv = settings_cache_csv,
    prefer_settings_cache = prefer_settings_cache,
    allow_manual_fallback = allow_manual_fallback
  )

  if (is.null(baseline_cfg_path)) {
    baseline_cfg_path <- file.path(
      output_root, "run_history", baseline_job_name_5arm,
      baseline_parent_job_id_5arm, "job_config.json"
    )
  }

  if (file.exists(baseline_cfg_path)) {
    base_cfg <- test_susine:::load_job_config(baseline_cfg_path)
    base_runs <- base_cfg$tables$runs %>%
      tibble::as_tibble()
    selected_runs <- filter_baseline_runs_for_5arm(base_runs, settings)
    dataset_ids <- sort(unique(as.character(selected_runs$dataset_bundle_id)))

    cfg <- base_cfg
    cfg$job$job_name <- study_name
    cfg$job$output_root <- output_root
    cfg$job$HPC <- FALSE
    cfg$job$aggregation_methods <- character(0)
    cfg$job$include_overall_pool <- FALSE
    cfg$job$write_confusion_bins <- FALSE
    cfg$job$write_snps_parquet <- FALSE
    cfg$job$write_tier_cs_metrics <- TRUE
    cfg$job$write_scaling_confusion_bins <- FALSE
    cfg$five_arm_settings <- settings
    cfg$five_arm_baseline_cfg_path <- baseline_cfg_path

    cfg$tables$dataset_bundles <- cfg$tables$dataset_bundles %>%
      tibble::as_tibble() %>%
      dplyr::mutate(dataset_bundle_id = as.character(.data$dataset_bundle_id)) %>%
      dplyr::filter(.data$dataset_bundle_id %in% dataset_ids) %>%
      dplyr::arrange(as.numeric(.data$dataset_bundle_id))
    cfg$tables$runs <- selected_runs

    data_matrix_catalog <- NULL
    if (!is.null(cfg$job$data_matrix_catalog)) {
      data_matrix_catalog <- cfg$job$data_matrix_catalog
    } else if (!is.null(cfg$tables$data_matrix_catalog)) {
      data_matrix_catalog <- cfg$tables$data_matrix_catalog
    } else {
      matrix_ids <- unique(cfg$tables$dataset_bundles$matrix_id)
      data_matrix_catalog <- test_susine:::build_data_matrix_catalog(
        requested_scenarios = data_scenario_5arm,
        repo_root = here_root,
        summary_path = file.path(
          here_root, "data", "sampled_simulated_genotypes",
          "scenario_sampling_summary.csv"
        )
      ) %>%
        dplyr::filter(.data$matrix_id %in% matrix_ids)
    }

    return(list(
      cfg = cfg,
      selected_runs = selected_runs,
      data_matrix_catalog = data_matrix_catalog,
      settings = settings
    ))
  }

  if (!isTRUE(allow_manual_fallback)) {
    stop(
      "Cannot copy baseline run rows because baseline job_config.json was ",
      "not found:\n", baseline_cfg_path,
      "\nProvide baseline_cfg_path or rerun/sync the baseline run history."
    )
  }

  get_setting <- function(method, col) {
    val <- settings[[col]][match(method, settings$method_family)]
    if (!length(val) || is.na(val[[1]])) {
      return(NA)
    }
    val[[1]]
  }

  sigma_vanilla <- as.numeric(get_setting("susie_vanilla",
                                          "sigma_0_2_scalar"))
  sigma_mu <- as.numeric(get_setting("susine_functional_mu",
                                     "sigma_0_2_scalar"))
  sigma_pi <- as.numeric(get_setting("susie_functional_pi",
                                     "sigma_0_2_scalar"))
  c_mu <- as.numeric(get_setting("susine_functional_mu", "c_value"))
  tau_pi <- as.numeric(get_setting("susie_functional_pi", "tau_value"))
  annotation_r2_mu <- as.numeric(get_setting("susine_functional_mu",
                                             "annotation_r2"))
  annotation_r2_pi <- as.numeric(get_setting("susie_functional_pi",
                                             "annotation_r2"))
  inflate_match_mu <- as.numeric(get_setting("susine_functional_mu",
                                             "inflate_match"))
  inflate_match_pi <- as.numeric(get_setting("susie_functional_pi",
                                             "inflate_match"))

  if (!isTRUE(all.equal(annotation_r2_mu, annotation_r2_pi, tolerance = 1e-8)) ||
      !isTRUE(all.equal(inflate_match_mu, inflate_match_pi, tolerance = 1e-8))) {
    stop(
      "Selected functional-mu and functional-pi settings must use the same ",
      "annotation regime. Got mu=(r2=", annotation_r2_mu,
      ", inflate=", inflate_match_mu, ") and pi=(r2=", annotation_r2_pi,
      ", inflate=", inflate_match_pi, ")."
    )
  }

  data_matrix_catalog_full <- test_susine:::build_data_matrix_catalog(
    requested_scenarios = data_scenario_5arm,
    repo_root = here_root,
    summary_path = file.path(
      here_root, "data", "sampled_simulated_genotypes",
      "scenario_sampling_summary.csv"
    )
  )

  m1_summary_path <- file.path(
    here_root, "data", "sampled_simulated_genotypes",
    "scenario_sampling_summary_with_m1.csv"
  )
  if (file.exists(m1_summary_path)) {
    catalog_with_m1 <- readr::read_csv(m1_summary_path, show_col_types = FALSE)
  } else {
    m1_tbl <- test_susine::build_x_metrics_table(
      matrix_paths_df = data_matrix_catalog_full %>%
        dplyr::select(.data$matrix_id, .data$matrix_path),
      base_dir = here_root
    ) %>%
      dplyr::select(.data$matrix_id, .data$M1)
    catalog_with_m1 <- data_matrix_catalog_full %>%
      dplyr::left_join(m1_tbl, by = "matrix_id")
    readr::write_csv(catalog_with_m1, m1_summary_path)
  }

  if (!"matrix_id" %in% names(catalog_with_m1)) {
    if ("gene" %in% names(catalog_with_m1)) {
      catalog_with_m1 <- dplyr::rename(catalog_with_m1,
                                       dataset_label = .data$gene)
    }
    join_col <- intersect(
      c("dataset_label", "matrix_path"),
      names(catalog_with_m1)
    )[[1L]]
    catalog_with_m1 <- catalog_with_m1 %>%
      dplyr::left_join(
        data_matrix_catalog_full %>%
          dplyr::select(.data$matrix_id, dplyr::all_of(join_col)),
        by = join_col
      )
  }

  catalog_with_m1 <- catalog_with_m1 %>%
    dplyr::filter(!is.na(.data$M1)) %>%
    dplyr::arrange(.data$M1)

  n_target <- min(n_matrices_target_5arm, nrow(catalog_with_m1))
  idx <- unique(round(seq(1, nrow(catalog_with_m1), length.out = n_target)))
  selected_matrices <- catalog_with_m1[idx, , drop = FALSE]

  data_matrix_catalog <- data_matrix_catalog_full %>%
    dplyr::semi_join(
      selected_matrices %>% dplyr::select(.data$matrix_id),
      by = "matrix_id"
    )

  spec_susie_vanilla <- list(
    name = "effect_match_susie_vanilla",
    use_case_ids = "susie_vanilla",
    exploration_methods = "sigma_0_2_grid",
    exploration_mode = "separate",
    K = 1L,
    sigma_0_2_grid_values = sigma_vanilla
  )

  spec_susine_functional_mu <- list(
    name = "effect_match_susine_functional_mu",
    use_case_ids = "susine_functional_mu",
    exploration_methods = c("c_grid", "sigma_0_2_grid"),
    exploration_mode = "intersect",
    K = 1L,
    c_grid_values = c_mu,
    sigma_0_2_grid_values = sigma_mu
  )

  spec_susie_functional_pi <- list(
    name = "effect_match_susie_functional_pi",
    use_case_ids = "susie_functional_pi",
    exploration_methods = c("tau_grid", "sigma_0_2_grid"),
    exploration_mode = "intersect",
    K = 1L,
    tau_grid_values = tau_pi,
    sigma_0_2_grid_values = sigma_pi
  )

  cfg <- test_susine::make_job_config(
    job_name = study_name,
    HPC = FALSE,
    time = "02:00:00",
    mem = "8G",
    output_root = output_root,
    grid_mode = "full",
    use_case_ids = cold_method_levels_5arm,
    exploration_specs = list(
      spec_susie_vanilla,
      spec_susine_functional_mu,
      spec_susie_functional_pi
    ),
    L_grid = L_value_5arm,
    y_noise_grid = y_noise_value_5arm,
    p_star_grid = p_star_value_5arm,
    architecture_grid = architecture_value_5arm,
    seeds = seeds_5arm,
    prior_quality = tibble::tibble(
      annotation_r2 = annotation_r2_mu,
      inflate_match = inflate_match_mu
    ),
    sigma_0_2_scalars = unique(c(sigma_vanilla, sigma_mu, sigma_pi)),
    data_scenarios = data_scenario_5arm,
    repo_root = here_root,
    data_matrix_catalog = data_matrix_catalog,
    sampled_scenario_summary = file.path(
      here_root, "data", "sampled_simulated_genotypes",
      "scenario_sampling_summary.csv"
    ),
    task_unit = "dataset",
    bundles_per_task = 1L,
    task_assignment_seed = 42L,
    max_iter = max_iter_5arm,
    aggregation_methods = character(0),
    include_overall_pool = FALSE,
    credible_set_rho = credible_set_rho_5arm,
    purity_threshold = purity_threshold_5arm,
    verbose_file_output = FALSE,
    write_snps_parquet = FALSE,
    write_confusion_bins = FALSE,
    write_tier_cs_metrics = TRUE,
    write_scaling_confusion_bins = FALSE,
    metrics_settings = list(
      z_top_k = 10L,
      jsd_threshold = 0.15
    )
  )
  cfg$five_arm_settings <- settings

  run_table <- cfg$tables$runs %>%
    tibble::as_tibble() %>%
    dplyr::mutate(
      method_family = .data$use_case_id,
      annotation_label = dplyr::if_else(
        is.na(.data$annotation_r2) | is.na(.data$inflate_match),
        "no_annotation",
        annotation_label_5arm(.data$annotation_r2, .data$inflate_match)
      ),
      setting_label = settings_label_5arm(
        .data$use_case_id, .data$c_value, .data$tau_value,
        .data$sigma_0_2_scalar
      )
    )
  cfg$tables$runs <- run_table

  selected_runs <- run_table %>%
    dplyr::filter(.data$use_case_id %in% cold_method_levels_5arm) %>%
    dplyr::arrange(.data$dataset_bundle_id, .data$use_case_id)

  list(
    cfg = cfg,
    selected_runs = selected_runs,
    data_matrix_catalog = data_matrix_catalog,
    settings = settings
  )
}
