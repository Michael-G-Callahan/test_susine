test_that("aggregation aliases normalize to the canonical jsd_050 id", {
  expect_equal(
    test_susine:::normalize_agg_method_id(c("cluster_weight", "cluster_weight_050", "cluster_weight_jsd_050")),
    c("cluster_weight", "cluster_weight_jsd_050", "cluster_weight_jsd_050")
  )
})

test_that("pooled AUPRC uses step-function average precision", {
  bins <- tibble::tibble(
    metric_group = "g1",
    pip_threshold = c(0.9, 0.8, 0.1),
    n_causal_at_bucket = c(1, 0, 1),
    n_noncausal_at_bucket = c(0, 2, 1)
  )

  step_ap <- test_susine:::auprc_from_pooled_bins(bins)
  precision <- c(1, 1 / 3, 0.5)
  recall <- c(0.5, 0.5, 1.0)
  trapz_ap <- sum(diff(c(0, recall)) * (c(1, head(precision, -1)) + precision) / 2)

  expect_equal(step_ap, 0.70, tolerance = 1e-12)
  expect_false(isTRUE(all.equal(step_ap, trapz_ap, tolerance = 1e-12)))

  pooled <- test_susine:::compute_auprc_from_pooled_confusion(
    pooled_bins = bins,
    group_vars = "metric_group"
  )
  expect_equal(pooled$AUPRC[[1]], 0.70, tolerance = 1e-12)
})

test_that("TPR@FPR helper matches pooled-confusion implementation", {
  bins <- tibble::tibble(
    metric_group = "g1",
    pip_threshold = c(0.9, 0.6, 0.2),
    n_causal_at_bucket = c(1, 1, 0),
    n_noncausal_at_bucket = c(0, 0, 20)
  )

  pooled <- test_susine:::compute_tpr05_from_pooled_confusion(
    pooled_bins = bins,
    group_vars = "metric_group"
  )

  expect_equal(
    test_susine:::tpr_at_fpr_threshold(bins, fpr_threshold = 0.05),
    pooled$tpr_fpr05[[1]],
    tolerance = 1e-12
  )
})

test_that("single-fit aggregation degenerates to the retained fit for every method", {
  pip <- c(0.8, 0.2, 0.05)

  out <- test_susine:::aggregate_use_case_pips(
    pip_list = list(pip),
    elbo_vec = 7,
    methods = c("uniform", "max_elbo", "elbo_softmax", "cluster_weight", "cluster_weight_jsd_050")
  )

  expect_setequal(
    names(out),
    c("uniform", "max_elbo", "elbo_softmax", "cluster_weight", "cluster_weight_jsd_050")
  )
  for (method in names(out)) {
    expect_equal(as.numeric(out[[method]]), pip, tolerance = 1e-12)
  }
})

test_that("collect backfills missing aggregated confusion for single-fit pure-refine groups", {
  confusion_individual <- tibble::tibble(
    run_id = 101L,
    explore_method = "refine",
    variant_id = 1L,
    agg_method = NA_character_,
    pip_threshold = c(0.9, 0.2),
    n_causal_at_bucket = c(1L, 0L),
    n_noncausal_at_bucket = c(0L, 3L),
    dataset_bundle_id = 7L,
    spec_name = "A-F",
    use_case_id = "susine_vanilla",
    annotation_r2 = NA_real_,
    group_key = "L=10|r2=NA|inflate=NA|explore=separate:refine",
    c_value = NA_real_,
    sigma_0_2_scalar = NA_character_
  )

  run_info <- tibble::tibble(
    run_id = 101L,
    task_id = 5L,
    dataset_bundle_id = 7L,
    spec_name = "A-F",
    use_case_id = "susine_vanilla",
    prior_spec_id = "susine_vanilla",
    exploration_methods = "refine",
    sigma_0_2_scalar = NA_character_,
    warm_method = NA_character_,
    c_value = NA_real_,
    refine_step = 1L,
    restart_id = NA_integer_,
    annotation_r2 = NA_real_,
    inflate_match = NA_real_,
    group_key = "L=10|r2=NA|inflate=NA|explore=separate:refine"
  )

  out <- test_susine:::backfill_single_fit_refine_agg_confusion(
    confusion_individual = confusion_individual,
    confusion_agg = confusion_individual[0, ],
    run_info = run_info,
    agg_methods = c("uniform", "max_elbo", "cluster_weight")
  )

  expect_equal(sort(unique(out$confusion_agg$agg_method)), c("cluster_weight", "max_elbo", "uniform"))
  expect_equal(nrow(out$repaired_groups), 3L)
  expect_equal(sum(out$repair_summary$n_backfilled_groups), 3L)
  expect_true(all(out$confusion_agg$explore_method == "aggregation"))
  expect_setequal(unique(out$confusion_agg$variant_id), c("uniform", "max_elbo", "cluster_weight"))
})

test_that("shared plot table preserves the global vanilla baseline row", {
  auprc_tbl <- test_susine:::build_shared_plot_metric_table(
    agg_overall = tibble::tibble(
      spec_name = "A-F",
      use_case_id = "susine_vanilla",
      agg_method = "cluster_weight_050",
      AUPRC = 0.20
    ),
    agg_by_r2 = tibble::tibble(
      spec_name = character(),
      use_case_id = character(),
      agg_method = character(),
      annotation_r2 = numeric(),
      AUPRC = numeric()
    ),
    individual_overall = tibble::tibble(
      spec_name = c("baseline-single", "truth-warm"),
      use_case_id = c("susine_vanilla", "susine_vanilla"),
      AUPRC = c(0.17, 0.24)
    ),
    individual_by_r2 = tibble::tibble(
      spec_name = character(),
      use_case_id = character(),
      annotation_r2 = numeric(),
      AUPRC = numeric()
    ),
    metric_col = "AUPRC",
    label = "test build"
  )

  expect_true(any(
    auprc_tbl$spec_name == "baseline-single" &
      auprc_tbl$series_type == "baseline_single" &
      auprc_tbl$use_case_id == "susine_vanilla" &
      is.na(auprc_tbl$annotation_r2) &
      abs(auprc_tbl$AUPRC - 0.17) < 1e-12
  ))
  expect_true(any(
    auprc_tbl$spec_name == "A-F" &
      auprc_tbl$agg_method == "cluster_weight_jsd_050"
  ))
})

test_that("terminal scaling overwrite replaces pure and interaction endpoints exactly", {
  summary_df <- tibble::tibble(
    spec_name = c("A-F", "A-F", "B-F", "B-F", "C-CS", "C-CS"),
    annotation_r2 = c(NA_real_, NA_real_, NA_real_, NA_real_, 0.5, 0.5),
    n_ensemble = c(16L, 32L, 16L, 32L, NA_integer_, NA_integer_),
    resolution = c(NA_character_, NA_character_, NA_character_, NA_character_, "4x4", "8x8"),
    agg_method = c("cluster_weight_jsd_050", "cluster_weight_jsd_050", "cluster_weight", "cluster_weight", "cluster_weight", "cluster_weight"),
    AUPRC_mean = c(0.23, 0.39, 0.34, 0.26, 0.24, 0.25),
    AUPRC_se = c(0.01, 0.01, 0.01, 0.01, NA_real_, NA_real_),
    tpr05_mean = c(0.41, 0.70, 0.44, 0.30, 0.35, 0.36),
    tpr05_se = c(0.01, 0.01, 0.01, 0.01, NA_real_, NA_real_),
    n_reps = c(5L, 5L, 5L, 5L, 1L, 1L)
  )

  targets <- tibble::tibble(
    spec_name = c("A-F", "B-F", "C-CS"),
    annotation_r2 = c(NA_real_, NA_real_, 0.5),
    agg_method = c("cluster_weight_jsd_050", "cluster_weight", "cluster_weight"),
    AUPRC = c(0.17, 0.27, 0.26),
    tpr_fpr05 = c(0.42, 0.31, 0.37)
  )

  out <- test_susine:::overwrite_scaling_terminal_metrics(summary_df, targets)

  expect_equal(out$AUPRC_mean[out$spec_name == "A-F" & out$n_ensemble == 32L], 0.17, tolerance = 1e-12)
  expect_equal(out$AUPRC_mean[out$spec_name == "B-F" & out$n_ensemble == 32L], 0.27, tolerance = 1e-12)
  expect_equal(out$AUPRC_mean[out$spec_name == "C-CS" & out$resolution == "8x8"], 0.26, tolerance = 1e-12)
  expect_true(is.na(out$AUPRC_se[out$spec_name == "A-F" & out$n_ensemble == 32L]))
  expect_equal(out$n_reps[out$spec_name == "A-F" & out$n_ensemble == 32L], 1L)
})
