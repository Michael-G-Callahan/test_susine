make_cs_grid_refit_catalog <- function(n_rows = 1L) {
  tibble::tibble(
    data_scenario = rep("simulation_n3", n_rows),
    dataset_label = paste0("simulation_n3_", seq_len(n_rows)),
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

make_test_buffer_ctx <- function() {
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
  buffer_ctx$base_output <- tempdir()
  buffer_ctx$task_id <- 1L
  buffer_ctx$flush_index <- 0L
  buffer_ctx
}

make_cs_grid_refit_run_tables <- function(...) {
  test_susine:::make_run_tables(
    use_case_ids = "susine_functional_mu",
    exploration_methods = "cs_grid_refit",
    exploration_mode = "separate",
    K = 2L,
    L_grid = 2L,
    y_noise_grid = 0.8,
    prior_quality = test_susine::prior_quality_grid(c(0.4), c(1)),
    p_star_grid = 2L,
    seeds = 1L,
    architecture_grid = "sparse",
    data_scenarios = "simulation_n3",
    data_matrix_catalog = make_cs_grid_refit_catalog(),
    ...,
    c_grid_values = c(0.1, 0.9),
    sigma_0_2_grid_values = c(0.05, 0.2)
  )
}

test_that("cs_grid_refit expands a c by sigma grid with warm refit rows", {
  tables <- make_cs_grid_refit_run_tables()

  expect_equal(nrow(tables$runs), 4L)
  expect_true(all(tables$runs$exploration_methods == "cs_grid_refit"))
  expect_true(all(tables$runs$run_type == "warm"))
  expect_true(all(tables$runs$warm_method == "cs_grid_refit"))
  expect_equal(sort(unique(tables$runs$c_value)), c(0.1, 0.9))
  expect_equal(sort(unique(as.numeric(tables$runs$sigma_0_2_scalar))), c(0.05, 0.2))
})

test_that("cs_grid_refit rejects mixed exploration and intersect mode", {
  common_args <- list(
    use_case_ids = "susine_functional_mu",
    K = 2L,
    L_grid = 2L,
    y_noise_grid = 0.8,
    prior_quality = test_susine::prior_quality_grid(c(0.4), c(1)),
    p_star_grid = 2L,
    seeds = 1L,
    architecture_grid = "sparse",
    data_scenarios = "simulation_n3",
    data_matrix_catalog = make_cs_grid_refit_catalog()
  )

  expect_error(
    do.call(
      test_susine:::make_run_tables,
      c(
        common_args,
        list(
          exploration_methods = c("cs_grid_refit", "restart"),
          exploration_mode = "separate"
        )
      )
    ),
    "standalone exploration method"
  )

  expect_error(
    do.call(
      test_susine:::make_run_tables,
      c(
        common_args,
        list(
          exploration_methods = "cs_grid_refit",
          exploration_mode = "intersect"
        )
      )
    ),
    "only supported with exploration_mode = 'separate'"
  )
})

test_that("cs_grid_refit stores the default-prior exact refit, not the search fit", {
  cfg <- test_susine::make_job_config(
    job_name = "unit_cs_grid_refit",
    use_case_ids = "susine_functional_mu",
    exploration_methods = "cs_grid_refit",
    exploration_mode = "separate",
    K = 1L,
    L_grid = 2L,
    y_noise_grid = 0.8,
    prior_quality = test_susine::prior_quality_grid(c(0.4), c(1)),
    p_star_grid = 2L,
    seeds = 1L,
    architecture_grid = "sparse",
    data_scenarios = "simulation_n3",
    repo_root = ".",
    data_matrix_catalog = make_cs_grid_refit_catalog(),
    task_unit = "dataset",
    bundles_per_task = 1L,
    runs_per_task = 10L,
    output_root = tempdir(),
    c_grid_values = 0.7,
    sigma_0_2_grid_values = 0.6,
    write_snps_parquet = FALSE,
    write_confusion_bins = FALSE,
    write_prior_diagnostics = TRUE,
    verbose_file_output = FALSE,
    include_overall_pool = FALSE,
    aggregation_methods = "uniform"
  )

  bundle_runs <- cfg$tables$runs
  run_row <- bundle_runs[1, , drop = FALSE]
  data_bundle <- test_susine:::generate_data_for_run(run_row, cfg)
  data_bundle$priors_cache <- new.env(parent = emptyenv())
  use_case <- test_susine:::lookup_use_case(cfg, run_row$use_case_id)

  source_result <- test_susine:::run_use_case(
    use_case = use_case,
    run_row = run_row,
    data_bundle = data_bundle,
    job_config = cfg
  )
  source_fit <- source_result$fits[[1]]
  exact_refit <- test_susine:::run_default_prior_refit(
    source_fit = source_fit,
    run_row = run_row,
    data_bundle = data_bundle,
    job_config = cfg
  )

  buffer_ctx <- make_test_buffer_ctx()
  test_susine:::execute_dataset_bundle(
    bundle_runs = bundle_runs,
    job_config = cfg,
    quiet = TRUE,
    buffer_ctx = buffer_ctx
  )

  model_tbl <- dplyr::bind_rows(buffer_ctx$buffers$model)
  prior_tbl <- dplyr::bind_rows(buffer_ctx$buffers$prior_diagnostics)
  stored_elbo <- unique(model_tbl$elbo_final)
  exact_elbo <- tail(exact_refit$fits[[1]]$model_fit$elbo, 1)

  expect_equal(length(stored_elbo), 1L)
  expect_equal(stored_elbo, exact_elbo, tolerance = 1e-8)
  expect_true(all(model_tbl$explore_method == "cs_grid_refit"))
  expect_true(all(model_tbl$model_call_executed))
  expect_true(all(abs(prior_tbl$mu_0_l) < 1e-12))
  expect_true(all(abs(prior_tbl$sigma_0_2_l / stats::var(data_bundle$y) - 0.2) < 1e-8))
  expect_false(isTRUE(all.equal(
    prior_tbl$sigma_0_2_l,
    rep(as.numeric(source_fit$priors$sigma_0_2[1, 1]), nrow(prior_tbl)),
    tolerance = 1e-8
  )))
})

test_that("cs_grid_refit scaling is treated as a c by sigma interaction", {
  pip_list <- list(
    c(0.9, 0.1, 0.05),
    c(0.8, 0.2, 0.05),
    c(0.4, 0.7, 0.05),
    c(0.3, 0.8, 0.05)
  )
  elbo_vec <- c(4, 3, 2, 1)
  meta_list <- list(
    list(run_id = 1L, spec_name = "C-CSR", annotation_r2 = 0.4, c_value = 0.1, sigma_0_2_scalar = 0.05),
    list(run_id = 2L, spec_name = "C-CSR", annotation_r2 = 0.4, c_value = 0.1, sigma_0_2_scalar = 0.2),
    list(run_id = 3L, spec_name = "C-CSR", annotation_r2 = 0.4, c_value = 0.9, sigma_0_2_scalar = 0.05),
    list(run_id = 4L, spec_name = "C-CSR", annotation_r2 = 0.4, c_value = 0.9, sigma_0_2_scalar = 0.2)
  )
  group_run_row <- tibble::tibble(
    exploration_methods = "cs_grid_refit",
    .planned_n_c = 2L,
    .planned_n_sigma = 2L
  )

  out <- test_susine:::compute_scaling_confusion_bins_for_group(
    pip_list = pip_list,
    elbo_vec = elbo_vec,
    meta_list = meta_list,
    causal_vec = c(1L, 0L, 0L),
    group_run_row = group_run_row,
    pip_breaks = c(0, 0.05, 0.1, 0.2, 0.5, 1),
    dataset_bundle_id = 1L
  )

  expect_false(is.null(out))
  expect_setequal(unique(out$resolution), c("1x1", "2x2"))
})

test_that("heterogeneous job configs preserve C-CSR as its own spec_name", {
  cfg <- test_susine::make_job_config(
    job_name = "unit_csr_spec_catalog",
    use_case_ids = "susine_vanilla",
    exploration_methods = "single",
    exploration_mode = "separate",
    K = 1L,
    L_grid = 2L,
    y_noise_grid = 0.8,
    prior_quality = test_susine::prior_quality_grid(c(0.4), c(1)),
    p_star_grid = 2L,
    seeds = 1L,
    architecture_grid = "sparse",
    data_scenarios = "simulation_n3",
    repo_root = ".",
    data_matrix_catalog = make_cs_grid_refit_catalog(),
    task_unit = "dataset",
    bundles_per_task = 1L,
    runs_per_task = 20L,
    output_root = tempdir(),
    write_snps_parquet = FALSE,
    write_confusion_bins = FALSE,
    verbose_file_output = FALSE,
    include_overall_pool = FALSE,
    aggregation_methods = "uniform",
    exploration_specs = list(
      list(
        name = "A-F",
        use_case_ids = "susine_vanilla",
        exploration_methods = "refine",
        exploration_mode = "separate",
        K = 2L,
        refine_settings = list(
          n_steps = 2L,
          cs_source = "filtered",
          purity_threshold = 0.95
        )
      ),
      list(
        name = "C-CSR",
        use_case_ids = "susine_functional_mu",
        exploration_methods = "cs_grid_refit",
        exploration_mode = "separate",
        K = 4L,
        c_grid_values = c(0.1, 0.9),
        sigma_0_2_grid_values = c(0.05, 0.2)
      )
    )
  )

  runs <- cfg$tables$runs
  csr_runs <- runs %>% dplyr::filter(.data$spec_name == "C-CSR")

  expect_setequal(unique(runs$spec_name), c("A-F", "C-CSR"))
  expect_equal(nrow(csr_runs), 4L)
  expect_true(all(csr_runs$exploration_methods == "cs_grid_refit"))
  expect_true(all(csr_runs$run_type == "warm"))
  expect_equal(sort(unique(csr_runs$c_value)), c(0.1, 0.9))
  expect_equal(sort(unique(as.numeric(csr_runs$sigma_0_2_scalar))), c(0.05, 0.2))
})
