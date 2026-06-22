make_annotation_test_catalog <- function(n_rows = 1L) {
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

make_annotation_buffer_ctx <- function() {
  buffer_ctx <- new.env(parent = emptyenv())
  buffer_ctx$buffers <- list(
    model = list(),
    effect = list(),
    effect_unfiltered = list(),
    tier_cs_metrics = list(),
    snps = list(),
    confusion = list(),
    dataset_metrics = list(),
    multimodal = list(),
    refine_depth = list(),
    scaling_bins = list(),
    validation = list(),
    prior_diagnostics = list(),
    annotation_diagnostics = list(),
    hg2_by_agg = list()
  )
  buffer_ctx$base_output <- tempdir()
  buffer_ctx$task_id <- 1L
  buffer_ctx$flush_index <- 0L
  buffer_ctx
}

test_that("annotation contamination preserves required block invariants", {
  set.seed(1)
  a <- rnorm(20)
  z <- rnorm(20)
  causal_idx <- c(2L, 5L, 9L)

  clean <- test_susine:::contaminate_annotation_vector(
    a_orig = a,
    z = z,
    causal_idx = causal_idx,
    arm = "null",
    lambda = 0
  )
  expect_equal(clean$mu_0, a, tolerance = 0)

  null_arm <- test_susine:::contaminate_annotation_vector(
    a_orig = a,
    z = z,
    causal_idx = causal_idx,
    arm = "null",
    lambda = 0.7
  )
  null_idx <- setdiff(seq_along(a), causal_idx)
  expect_equal(null_arm$mu_0[causal_idx], a[causal_idx], tolerance = 0)
  expect_equal(
    test_susine:::annotation_rms(null_arm$mu_0[null_idx]),
    test_susine:::annotation_rms(a[null_idx]),
    tolerance = 1e-12
  )
  expect_true(null_arm$diagnostics$causal_unchanged)
  expect_lt(null_arm$diagnostics$null_rms_rel_error, 1e-12)

  causal_arm <- test_susine:::contaminate_annotation_vector(
    a_orig = a,
    z = z,
    causal_idx = causal_idx,
    arm = "causal",
    lambda = 0.7
  )
  expect_equal(causal_arm$mu_0[null_idx], a[null_idx], tolerance = 0)
  expect_equal(
    test_susine:::annotation_rms(causal_arm$mu_0[causal_idx]),
    test_susine:::annotation_rms(a[causal_idx]),
    tolerance = 1e-12
  )
  expect_true(causal_arm$diagnostics$null_unchanged)
  expect_lt(causal_arm$diagnostics$causal_rms_rel_error, 1e-12)
})

test_that("shuffle_z preserves the caller RNG stream", {
  set.seed(11)
  before <- .Random.seed
  a <- stats::rnorm(30)
  z <- stats::rnorm(30)
  causal_idx <- c(3L, 8L, 20L)
  set.seed(1234)
  expected <- {
    stats::runif(4)
  }
  set.seed(1234)
  test_susine:::contaminate_annotation_vector(
    a_orig = a,
    z = z,
    causal_idx = causal_idx,
    arm = "null",
    lambda = 0.5,
    shuffle_z = TRUE,
    seed = 99L
  )
  observed <- stats::runif(4)
  expect_equal(observed, expected, tolerance = 0)
  expect_false(identical(.Random.seed, before))
})

test_that("lambda-zero annotation contamination reproduces clean C-CS metrics end to end", {
  prior_quality <- tidyr::expand_grid(
    annotation_r2 = 0.3,
    inflate_match = 0.95,
    annotation_contamination_arm = c("none", "null"),
    annotation_contamination_lambda = 0,
    annotation_contamination_shuffle_z = FALSE
  )

  cfg <- test_susine::make_job_config(
    job_name = "unit_annotation_contamination_master_check",
    use_case_ids = "susine_functional_mu",
    exploration_specs = list(
      list(
        name = "C-CS-master",
        use_case_ids = "susine_functional_mu",
        exploration_methods = c("c_grid", "sigma_0_2_grid"),
        exploration_mode = "intersect",
        K = 4L,
        c_grid_values = c(0, 0.43),
        sigma_0_2_grid_values = c(0.1, 1)
      )
    ),
    L_grid = 2L,
    y_noise_grid = NA_real_,
    prior_quality = prior_quality,
    p_star_grid = NA_integer_,
    seeds = 1L,
    architecture_grid = "susie2_oligogenic",
    data_scenarios = "simulation_n3",
    repo_root = ".",
    data_matrix_catalog = make_annotation_test_catalog(),
    task_unit = "dataset",
    bundles_per_task = 1L,
    output_root = tempdir(),
    write_snps_parquet = FALSE,
    write_confusion_bins = TRUE,
    write_prior_diagnostics = TRUE,
    verbose_file_output = FALSE,
    include_overall_pool = FALSE,
    aggregation_methods = "cluster_weight_credible",
    max_iter = 40L
  )

  buffer_ctx <- make_annotation_buffer_ctx()
  test_susine:::execute_dataset_bundle(
    bundle_runs = cfg$tables$runs,
    job_config = cfg,
    quiet = TRUE,
    buffer_ctx = buffer_ctx
  )

  run_lookup <- cfg$tables$runs %>%
    dplyr::select(
      run_id,
      spec_name,
      annotation_contamination_arm,
      annotation_contamination_lambda
    )
  confusion <- dplyr::bind_rows(buffer_ctx$buffers$confusion) %>%
    dplyr::left_join(run_lookup, by = "run_id") %>%
    dplyr::filter(
      .data$spec_name == "C-CS-master",
      .data$explore_method == "aggregation",
      .data$agg_method == "cluster_weight_credible"
    )
  expect_setequal(unique(confusion$annotation_contamination_arm), c("none", "null"))

  auprc <- test_susine::compute_auprc_from_confusion(
    confusion,
    group_vars = "annotation_contamination_arm"
  )
  recall <- test_susine::compute_recall_at_precision_from_confusion(
    confusion,
    group_vars = "annotation_contamination_arm",
    precision_threshold = 0.75
  )

  clean_auprc <- auprc$AUPRC[auprc$annotation_contamination_arm == "none"]
  null_auprc <- auprc$AUPRC[auprc$annotation_contamination_arm == "null"]
  clean_recall <- recall$recall_at_precision[recall$annotation_contamination_arm == "none"]
  null_recall <- recall$recall_at_precision[recall$annotation_contamination_arm == "null"]
  expect_equal(null_auprc, clean_auprc, tolerance = 1e-12)
  expect_equal(null_recall, clean_recall, tolerance = 1e-12)

  prior_diag <- dplyr::bind_rows(buffer_ctx$buffers$prior_diagnostics) %>%
    dplyr::left_join(cfg$tables$runs, by = "run_id") %>%
    dplyr::filter(
      .data$spec_name == "C-CS-master",
      .data$annotation_contamination_arm == "null"
    )
  expect_gt(nrow(prior_diag), 0L)
  expect_true(all(abs(prior_diag$c_l - 1) < 1e-12))

  ann_diag <- dplyr::bind_rows(buffer_ctx$buffers$annotation_diagnostics)
  expect_gt(nrow(ann_diag), 0L)
  expect_true(all(ann_diag$anchor_sigma_0_2_scalar == min(c(0.1, 1))))
})

test_that("annotation contamination settings enter C-CS group keys", {
  catalog <- tibble::tibble(
    data_scenario = "simulation_n3",
    dataset_label = "unit",
    participant_count = NA_integer_,
    snps_post = NA_integer_,
    snp_set = NA_character_,
    matrix_path = NA_character_,
    manifest_path = NA_character_,
    source = "simulation",
    matrix_index = 1L,
    matrix_id = 1L
  )

  prior_quality <- tidyr::expand_grid(
    annotation_r2 = 0.3,
    inflate_match = 0.95,
    annotation_contamination_arm = "null",
    annotation_contamination_lambda = c(0, 0.5),
    annotation_contamination_shuffle_z = FALSE
  )

  cfg <- test_susine::make_job_config(
    job_name = "unit_annotation_contamination",
    use_case_ids = c("susine_vanilla", "susine_functional_mu"),
    exploration_specs = list(
      list(
        name = "C-CS-test",
        use_case_ids = "susine_functional_mu",
        exploration_methods = c("c_grid", "sigma_0_2_grid"),
        exploration_mode = "intersect",
        K = 4L,
        c_grid_values = c(0, 1),
        sigma_0_2_grid_values = c(0.1, 1)
      ),
      list(
        name = "baseline-single",
        use_case_ids = "susine_vanilla",
        exploration_methods = "single",
        exploration_mode = "separate",
        K = 1L
      )
    ),
    L_grid = 2L,
    y_noise_grid = NA_real_,
    prior_quality = prior_quality,
    p_star_grid = NA_integer_,
    seeds = 1L,
    architecture_grid = "susie2_oligogenic",
    data_scenarios = "simulation_n3",
    repo_root = ".",
    data_matrix_catalog = catalog,
    task_unit = "dataset",
    bundles_per_task = 1L,
    output_root = tempdir(),
    write_snps_parquet = FALSE,
    write_confusion_bins = TRUE,
    verbose_file_output = FALSE,
    include_overall_pool = FALSE,
    aggregation_methods = "cluster_weight_credible"
  )

  ccs_runs <- cfg$tables$runs %>% dplyr::filter(.data$spec_name == "C-CS-test")
  expect_equal(sort(unique(ccs_runs$annotation_contamination_lambda)), c(0, 0.5))
  expect_equal(dplyr::n_distinct(ccs_runs$group_key), 2L)
  expect_true(any(grepl("ann_lambda=0.5", ccs_runs$group_key, fixed = TRUE)))

  baseline_runs <- cfg$tables$runs %>% dplyr::filter(.data$spec_name == "baseline-single")
  expect_equal(nrow(baseline_runs), 1L)
  expect_true(all(is.na(baseline_runs$annotation_contamination_arm)))
})

test_that("recall at precision reports fallback when threshold is unreachable", {
  bins <- tibble::tibble(
    metric_group = "g1",
    pip_threshold = c(0.9, 0.5, 0.1),
    n_causal_at_bucket = c(0L, 1L, 1L),
    n_noncausal_at_bucket = c(1L, 1L, 0L)
  )

  out <- test_susine::compute_recall_at_precision_from_confusion(
    bins,
    group_vars = "metric_group",
    precision_threshold = 0.75
  )
  expect_equal(out$recall_at_precision, 0)
  expect_false(out$precision_threshold_reached)
  expect_equal(out$max_achievable_precision, 0.5, tolerance = 1e-12)
})