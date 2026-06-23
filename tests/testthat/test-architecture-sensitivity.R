make_architecture_test_catalog <- function(n_rows = 1L) {
  tibble::tibble(
    data_scenario = rep("simulation_n3", n_rows),
    dataset_label = paste0("architecture_test_", seq_len(n_rows)),
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

test_that("sparse architecture sensitivity grid preserves k and per-causal h2", {
  cfg <- test_susine::make_job_config(
    job_name = "unit_architecture_sparse",
    use_case_ids = "susine_functional_mu",
    exploration_methods = "single",
    exploration_mode = "separate",
    K = 1L,
    L_grid = 2L,
    y_noise_grid = NA_real_,
    prior_quality = test_susine::prior_quality_grid(c(0.3), c(0.95)),
    p_star_grid = c(1L, 3L, 5L),
    h2_snp_per_causal_grid = c(0.01, 0.03),
    h2_total_grid = c(0.05, 0.25),
    diffuse_k_grid = c(50L, 100L),
    seeds = 1:2,
    architecture_grid = "susie2_sparse",
    data_scenarios = "simulation_n3",
    repo_root = ".",
    data_matrix_catalog = make_architecture_test_catalog(2L),
    task_unit = "dataset",
    bundles_per_task = 1L,
    output_root = tempdir(),
    write_snps_parquet = FALSE,
    write_confusion_bins = FALSE,
    verbose_file_output = FALSE,
    include_overall_pool = FALSE,
    aggregation_methods = character(0)
  )

  bundles <- cfg$tables$dataset_bundles
  expect_equal(nrow(bundles), 2L * 2L * 3L * 2L)
  expect_setequal(bundles$p_star, c(1L, 3L, 5L))
  expect_setequal(bundles$h2_snp_per_causal, c(0.01, 0.03))
  expect_true(all(is.na(bundles$h2_total)))
  expect_true(all(is.na(bundles$diffuse_k)))
  expect_true(all(is.na(bundles$y_noise)))
})

test_that("diffuse architecture sensitivity grid keeps only total h2 and diffuse k", {
  cfg <- test_susine::make_job_config(
    job_name = "unit_architecture_diffuse",
    use_case_ids = "susine_functional_mu",
    exploration_methods = "single",
    exploration_mode = "separate",
    K = 1L,
    L_grid = 2L,
    y_noise_grid = NA_real_,
    prior_quality = test_susine::prior_quality_grid(c(0.3), c(0.95)),
    p_star_grid = c(1L, 3L, 5L),
    h2_snp_per_causal_grid = c(0.01, 0.03),
    h2_total_grid = c(0.05, 0.15, 0.25),
    diffuse_k_grid = c(100L),
    seeds = 1:2,
    architecture_grid = "susie2_diffuse",
    data_scenarios = "simulation_n3",
    repo_root = ".",
    data_matrix_catalog = make_architecture_test_catalog(2L),
    task_unit = "dataset",
    bundles_per_task = 1L,
    output_root = tempdir(),
    write_snps_parquet = FALSE,
    write_confusion_bins = FALSE,
    verbose_file_output = FALSE,
    include_overall_pool = FALSE,
    aggregation_methods = character(0)
  )

  bundles <- cfg$tables$dataset_bundles
  expect_equal(nrow(bundles), 2L * 2L * 3L)
  expect_true(all(is.na(bundles$p_star)))
  expect_true(all(is.na(bundles$h2_snp_per_causal)))
  expect_setequal(bundles$h2_total, c(0.05, 0.15, 0.25))
  expect_true(all(bundles$diffuse_k == 100L))
  expect_true(all(is.na(bundles$y_noise)))
})

test_that("diffuse architecture data generation uses requested causal count and total h2", {
  cfg <- test_susine::make_job_config(
    job_name = "unit_architecture_diffuse_data",
    use_case_ids = "susine_functional_mu",
    exploration_methods = "single",
    exploration_mode = "separate",
    K = 1L,
    L_grid = 2L,
    y_noise_grid = NA_real_,
    prior_quality = test_susine::prior_quality_grid(c(0.3), c(0.95)),
    p_star_grid = NA_integer_,
    h2_total_grid = 0.15,
    diffuse_k_grid = 100L,
    seeds = 1L,
    architecture_grid = "susie2_diffuse",
    data_scenarios = "simulation_n3",
    repo_root = ".",
    data_matrix_catalog = make_architecture_test_catalog(1L),
    task_unit = "dataset",
    bundles_per_task = 1L,
    output_root = tempdir(),
    write_snps_parquet = FALSE,
    write_confusion_bins = FALSE,
    verbose_file_output = FALSE,
    include_overall_pool = FALSE,
    aggregation_methods = character(0)
  )

  bundle <- cfg$tables$dataset_bundles[1, ]
  data_bundle <- test_susine:::generate_data_for_bundle(bundle, cfg)
  expected_k <- min(100L, ncol(data_bundle$X))

  expect_identical(data_bundle$architecture, "susie2_diffuse")
  expect_equal(length(data_bundle$causal_idx), expected_k)
  expect_equal(sum(data_bundle$effect_tier == "diffuse"), expected_k)
  expect_equal(sum(data_bundle$beta != 0), expected_k)
  expect_true(is.finite(data_bundle$sigma2))
  expect_gt(data_bundle$sigma2, 0)
})
test_that("sparse sign-prior mode avoids constant-beta annotation degeneracy", {
  beta_single <- c(1, 0, 0, 0)
  prior_single <- test_susine:::simulate_priors(
    beta = beta_single,
    annotation_r2 = 0.3,
    inflate_match = 0.95,
    effect_sd = 1,
    annotation_signal_mode = "sign",
    seed = 3L
  )
  expect_false(isTRUE(all.equal(prior_single$mu_0[1], beta_single[1], tolerance = 0)))

  beta_two <- c(1, 1, 0, 0)
  expect_no_error(
    prior_two <- test_susine:::simulate_priors(
      beta = beta_two,
      annotation_r2 = 0.3,
      inflate_match = 0.95,
      effect_sd = 1,
      annotation_signal_mode = "sign",
      seed = 4L
    )
  )
  expect_true(all(is.finite(prior_two$mu_0)))
  expect_gt(stats::sd(prior_two$mu_0[beta_two != 0]), 0)
})