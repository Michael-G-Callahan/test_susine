make_baseline_screen_catalog <- function(n_rows = 1L) {
  tibble::tibble(
    data_scenario = rep("simulation_n3", n_rows),
    dataset_label = paste0("baseline_screen_", seq_len(n_rows)),
    participant_count = rep(NA_integer_, n_rows),
    snps_post = rep(NA_integer_, n_rows),
    snp_set = rep(NA_character_, n_rows),
    matrix_path = rep(NA_character_, n_rows),
    manifest_path = rep(NA_character_, n_rows),
    source = rep("simulation", n_rows),
    matrix_index = seq_len(n_rows),
    matrix_id = seq_len(n_rows)
  )
}

test_that("baseline screening specs expand tuning grids and annotation cross-products", {
  prior_quality <- tidyr::expand_grid(
    annotation_r2 = c(0, 0.3),
    inflate_match = c(1, 0.8)
  )

  cfg <- test_susine::make_job_config(
    job_name = "unit_baseline_screen",
    use_case_ids = c("susie_vanilla", "susine_functional_mu"),
    exploration_specs = list(
      list(
        name = "susie_vanilla_sigma",
        use_case_ids = "susie_vanilla",
        exploration_methods = "sigma_0_2_grid",
        exploration_mode = "separate",
        K = 3L,
        sigma_0_2_grid_values = c(0.1, 0.2, 0.4)
      ),
      list(
        name = "susine_functional_mu_grid",
        use_case_ids = "susine_functional_mu",
        exploration_methods = c("c_grid", "sigma_0_2_grid"),
        exploration_mode = "intersect",
        K = 4L,
        c_grid_values = c(0, 1),
        sigma_0_2_grid_values = c(0.1, 0.2)
      )
    ),
    L_grid = 2L,
    y_noise_grid = 0.8,
    prior_quality = prior_quality,
    p_star_grid = 2L,
    seeds = 1L,
    architecture_grid = "sparse",
    data_scenarios = "simulation_n3",
    repo_root = ".",
    data_matrix_catalog = make_baseline_screen_catalog(),
    task_unit = "dataset",
    bundles_per_task = 1L,
    output_root = tempdir(),
    verbose_file_output = FALSE,
    write_snps_parquet = FALSE,
    write_confusion_bins = TRUE,
    write_tier_cs_metrics = TRUE,
    include_overall_pool = FALSE,
    aggregation_methods = character(0),
    write_scaling_confusion_bins = FALSE
  )

  runs <- cfg$tables$runs

  expect_setequal(
    unique(runs$spec_name),
    c("susie_vanilla_sigma", "susine_functional_mu_grid")
  )

  vanilla_runs <- runs %>% dplyr::filter(.data$spec_name == "susie_vanilla_sigma")
  expect_equal(nrow(vanilla_runs), 3L)
  expect_true(all(is.na(vanilla_runs$annotation_r2)))
  expect_true(all(is.na(vanilla_runs$inflate_match)))
  expect_equal(sort(unique(as.numeric(vanilla_runs$sigma_0_2_scalar))), c(0.1, 0.2, 0.4))

  functional_mu_runs <- runs %>% dplyr::filter(.data$spec_name == "susine_functional_mu_grid")
  expect_equal(nrow(functional_mu_runs), 16L)
  expect_equal(sort(unique(functional_mu_runs$c_value)), c(0, 1))
  expect_equal(sort(unique(as.numeric(functional_mu_runs$sigma_0_2_scalar))), c(0.1, 0.2))
  expect_equal(sort(unique(functional_mu_runs$annotation_r2)), c(0, 0.3))
  expect_equal(sort(unique(functional_mu_runs$inflate_match)), c(0.8, 1))
})

test_that("aggregate_staging_outputs preserves filtered and unfiltered effect metrics", {
  root <- file.path(tempdir(), "baseline-screen-agg")
  job_name <- "unit_baseline_collect"
  parent_job_id <- "12345"

  staging_dir <- file.path(root, "slurm_output", job_name, parent_job_id, "task-001")
  history_dir <- file.path(root, "run_history", job_name, parent_job_id)
  dir.create(staging_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(history_dir, recursive = TRUE, showWarnings = FALSE)

  run_lookup <- tibble::tibble(
    run_id = 1L,
    use_case_id = "susie_vanilla",
    spec_name = "susie_vanilla_sigma",
    phenotype_seed = 1L,
    dataset_bundle_id = 11L,
    architecture = "sparse",
    refine_step = NA_integer_,
    run_type = "cold",
    group_key = "base",
    L = 10L,
    annotation_r2 = NA_real_,
    inflate_match = NA_real_,
    sigma_0_2_scalar = "0.2",
    c_value = NA_real_,
    tau_value = 0.464,
    matrix_id = 1L,
    y_noise = 0.8,
    p_star = 2L,
    exploration_mode = "separate",
    exploration_methods = "sigma_0_2_grid",
    exploration_group = "susie_vanilla_sigma",
    alpha_concentration = NA_real_,
    warm_method = NA_character_
  )
  readr::write_csv(run_lookup, file.path(history_dir, "run_table.csv"))

  effect_row <- tibble::tibble(
    run_id = 1L,
    effect = 1L,
    size = 2,
    purity = 0.9,
    coverage = 1L,
    effect_pip_entropy = 0.5,
    effect_pip_entropy_core95 = 0.4,
    effect_k_eff_signal = 1.6,
    effect_k_eff_signal_core95 = 1.5,
    tail_inflation_ratio = 1.067,
    tail_inflation_log = 0.1,
    accuracy_ratio = 0.8,
    tau_value = NA_real_,
    task_id = 1L,
    flush_id = "flush-001"
  )

  readr::write_csv(effect_row, file.path(staging_dir, "flush-001_effect_metrics.csv"))
  readr::write_csv(effect_row, file.path(staging_dir, "flush-001_effect_metrics_unfiltered.csv"))

  out_dir <- file.path(root, "aggregated-test")
  test_susine::aggregate_staging_outputs(
    job_name = job_name,
    parent_job_id = parent_job_id,
    output_root = root,
    validate = TRUE,
    output_dir = out_dir,
    write_snps_dataset = FALSE
  )

  filtered_path <- file.path(out_dir, "effect_metrics.csv")
  unfiltered_path <- file.path(out_dir, "effect_metrics_unfiltered.csv")

  expect_true(file.exists(filtered_path))
  expect_true(file.exists(unfiltered_path))

  filtered_tbl <- readr::read_csv(filtered_path, show_col_types = FALSE)
  unfiltered_tbl <- readr::read_csv(unfiltered_path, show_col_types = FALSE)

  expect_equal(nrow(filtered_tbl), 1L)
  expect_equal(nrow(unfiltered_tbl), 1L)
  expect_true("use_case_id" %in% names(filtered_tbl))
  expect_true("use_case_id" %in% names(unfiltered_tbl))

  expected_context <- c(
    "spec_name", "annotation_r2", "inflate_match",
    "sigma_0_2_scalar", "c_value", "tau_value"
  )
  expect_true(all(expected_context %in% names(filtered_tbl)))
  expect_true(all(expected_context %in% names(unfiltered_tbl)))
  expect_equal(filtered_tbl$tau_value, 0.464, tolerance = 1e-12)
  expect_equal(unfiltered_tbl$tau_value, 0.464, tolerance = 1e-12)
})
