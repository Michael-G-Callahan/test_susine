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

# Treatment-arm settings (high-quality annotation calibration point).
annotation_r2_value_5arm <- 0.5
inflate_match_value_5arm <- 0.9
c_value_mu_5arm    <- 0.6
tau_value_pi_5arm  <- 0.464
sigma_0_2_value_5arm <- 0.01
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
                           study_name = study_name_5arm) {
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
    sigma_0_2_grid_values = sigma_0_2_value_5arm
  )

  spec_susine_functional_mu <- list(
    name = "effect_match_susine_functional_mu",
    use_case_ids = "susine_functional_mu",
    exploration_methods = c("c_grid", "sigma_0_2_grid"),
    exploration_mode = "intersect",
    K = 1L,
    c_grid_values = c_value_mu_5arm,
    sigma_0_2_grid_values = sigma_0_2_value_5arm
  )

  spec_susie_functional_pi <- list(
    name = "effect_match_susie_functional_pi",
    use_case_ids = "susie_functional_pi",
    exploration_methods = c("tau_grid", "sigma_0_2_grid"),
    exploration_mode = "intersect",
    K = 1L,
    tau_grid_values = tau_value_pi_5arm,
    sigma_0_2_grid_values = sigma_0_2_value_5arm
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
      annotation_r2 = annotation_r2_value_5arm,
      inflate_match = inflate_match_value_5arm
    ),
    sigma_0_2_scalars = sigma_0_2_value_5arm,
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

  run_table <- cfg$tables$runs %>%
    tibble::as_tibble() %>%
    dplyr::mutate(
      method_family = .data$use_case_id,
      annotation_label = dplyr::if_else(
        is.na(.data$annotation_r2) | is.na(.data$inflate_match),
        "no_annotation",
        paste0("r2=", .data$annotation_r2,
               " | inflate=", .data$inflate_match)
      ),
      setting_label = dplyr::case_when(
        .data$use_case_id == "susie_vanilla" ~
          paste0("susie_vanilla | sigma=", .data$sigma_0_2_scalar),
        .data$use_case_id == "susine_functional_mu" ~
          paste0("susine_functional_mu | c=", .data$c_value,
                 " | sigma=", .data$sigma_0_2_scalar),
        .data$use_case_id == "susie_functional_pi" ~
          paste0("susie_functional_pi | tau=", .data$tau_value,
                 " | sigma=", .data$sigma_0_2_scalar),
        TRUE ~ .data$use_case_id
      )
    )
  cfg$tables$runs <- run_table

  selected_runs <- run_table %>%
    dplyr::filter(.data$use_case_id %in% cold_method_levels_5arm) %>%
    dplyr::arrange(.data$dataset_bundle_id, .data$use_case_id)

  list(
    cfg = cfg,
    selected_runs = selected_runs,
    data_matrix_catalog = data_matrix_catalog
  )
}