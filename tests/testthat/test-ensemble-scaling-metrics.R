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
    agg_methods = c("uniform", "max_elbo", "cluster_weight", "cluster_weight_jsd_050")
  )

  expect_equal(
    sort(unique(out$confusion_agg$agg_method)),
    c("cluster_weight", "cluster_weight_jsd_050", "max_elbo", "uniform")
  )
  expect_equal(nrow(out$repaired_groups), 4L)
  expect_equal(sum(out$repair_summary$n_backfilled_groups), 4L)
  expect_true(all(out$confusion_agg$explore_method == "aggregation"))
  expect_setequal(
    unique(out$confusion_agg$variant_id),
    c("uniform", "max_elbo", "cluster_weight", "cluster_weight_jsd_050")
  )
})

test_that("collect backfills missing aggregated confusion from terminal scaling bins", {
  scaling_bins_raw <- tibble::tibble(
    dataset_bundle_id = c(9L, 9L, 9L, 9L),
    spec_name = c("A-RFS", "A-RFS", "A-RFS", "A-RFS"),
    annotation_r2 = c(NA_real_, NA_real_, NA_real_, NA_real_),
    n_ensemble = c(16L, 16L, 32L, 32L),
    resolution = c(NA_character_, NA_character_, NA_character_, NA_character_),
    agg_method = c("cluster_weight_jsd_050", "cluster_weight_jsd_050", "cluster_weight_jsd_050", "cluster_weight_jsd_050"),
    rep = c(1L, 1L, 1L, 1L),
    pip_threshold = c(0.9, 0.2, 0.9, 0.2),
    n_causal_at_bucket = c(1L, 0L, 2L, 0L),
    n_noncausal_at_bucket = c(0L, 2L, 0L, 3L)
  )

  run_info <- tibble::tibble(
    run_id = 201L,
    task_id = 12L,
    dataset_bundle_id = 9L,
    spec_name = "A-RFS",
    use_case_id = "susine_vanilla",
    prior_spec_id = "susine_vanilla",
    exploration_methods = "restartxrefine",
    sigma_0_2_scalar = NA_character_,
    warm_method = NA_character_,
    c_value = NA_real_,
    refine_step = NA_integer_,
    restart_id = NA_integer_,
    annotation_r2 = NA_real_,
    inflate_match = NA_real_,
    group_key = "L=10|r2=NA|inflate=NA|explore=intersect:restartxrefine"
  )

  out <- test_susine:::backfill_terminal_scaling_agg_confusion(
    scaling_bins_raw = scaling_bins_raw,
    confusion_agg = tibble::tibble(
      run_id = integer(),
      explore_method = character(),
      variant_id = character(),
      agg_method = character(),
      pip_threshold = numeric(),
      n_causal_at_bucket = integer(),
      n_noncausal_at_bucket = integer(),
      dataset_bundle_id = integer(),
      spec_name = character(),
      use_case_id = character(),
      annotation_r2 = numeric(),
      group_key = character(),
      c_value = numeric(),
      sigma_0_2_scalar = character()
    ),
    run_info = run_info,
    agg_methods = c("uniform", "max_elbo", "elbo_softmax", "cluster_weight", "cluster_weight_jsd_050")
  )

  expect_equal(nrow(out$repaired_groups), 1L)
  expect_equal(out$repaired_groups$agg_method[[1]], "cluster_weight_jsd_050")
  expect_equal(sort(unique(out$confusion_agg$pip_threshold)), c(0.2, 0.9))
  expect_true(all(out$confusion_agg$agg_method == "cluster_weight_jsd_050"))
  expect_true(all(out$confusion_agg$explore_method == "aggregation"))
})

test_that("collect backfills missing cluster_weight_jsd_050 from staged SNPs", {
  skip_if_not_installed("arrow")

  snp_tbl <- tibble::tibble(
    run_id = c(301L, 301L, 301L, 302L, 302L, 302L),
    snp_index = c(1L, 2L, 3L, 1L, 2L, 3L),
    pip = c(0.9, 0.2, 0.1, 0.8, 0.3, 0.05),
    causal = c(1L, 0L, 0L, 1L, 0L, 0L)
  )
  pq_path <- file.path(tempdir(), "test-jsd050-snps.parquet")
  arrow::write_parquet(snp_tbl, pq_path)

  confusion_individual <- tibble::tibble(
    run_id = c(301L, 301L, 302L, 302L),
    explore_method = c("refine", "refine", "refine", "refine"),
    variant_id = c("0", "0", "1", "1"),
    agg_method = NA_character_,
    pip_threshold = c(0.9, 0.2, 0.8, 0.3),
    n_causal_at_bucket = c(1L, 0L, 1L, 0L),
    n_noncausal_at_bucket = c(0L, 2L, 0L, 2L),
    dataset_bundle_id = c(17L, 17L, 17L, 17L),
    spec_name = c("A-F", "A-F", "A-F", "A-F"),
    use_case_id = c("susine_vanilla", "susine_vanilla", "susine_vanilla", "susine_vanilla"),
    annotation_r2 = c(NA_real_, NA_real_, NA_real_, NA_real_),
    group_key = c(
      "L=10|r2=NA|inflate=NA|explore=separate:refine",
      "L=10|r2=NA|inflate=NA|explore=separate:refine",
      "L=10|r2=NA|inflate=NA|explore=separate:refine",
      "L=10|r2=NA|inflate=NA|explore=separate:refine"
    ),
    c_value = c(NA_real_, NA_real_, NA_real_, NA_real_),
    sigma_0_2_scalar = c(NA_character_, NA_character_, NA_character_, NA_character_)
  )

  run_info <- tibble::tibble(
    run_id = c(301L, 302L),
    task_id = c(21L, 21L),
    dataset_bundle_id = c(17L, 17L),
    spec_name = c("A-F", "A-F"),
    use_case_id = c("susine_vanilla", "susine_vanilla"),
    prior_spec_id = c("susine_vanilla", "susine_vanilla"),
    exploration_methods = c("refine", "refine"),
    sigma_0_2_scalar = c(NA_character_, NA_character_),
    warm_method = c(NA_character_, NA_character_),
    c_value = c(NA_real_, NA_real_),
    refine_step = c(0L, 1L),
    restart_id = c(NA_integer_, NA_integer_),
    annotation_r2 = c(NA_real_, NA_real_),
    inflate_match = c(NA_real_, NA_real_),
    group_key = c(
      "L=10|r2=NA|inflate=NA|explore=separate:refine",
      "L=10|r2=NA|inflate=NA|explore=separate:refine"
    )
  )

  model_metrics <- tibble::tibble(
    run_id = c(301L, 302L),
    agg_method = c(NA_character_, NA_character_),
    elbo_final = c(10, 8)
  )

  out <- test_susine:::backfill_missing_agg_confusion_from_snps(
    snp_files = pq_path,
    confusion_individual = confusion_individual,
    confusion_agg = confusion_individual[0, ],
    run_info = run_info,
    model_metrics = model_metrics,
    agg_methods = "cluster_weight_jsd_050",
    pip_breaks = c(0, 0.05, 0.1, 0.2, 0.3, 0.8, 0.9, 1)
  )

  expect_equal(nrow(out$repaired_groups), 1L)
  expect_equal(out$repaired_groups$agg_method[[1]], "cluster_weight_jsd_050")
  expect_true(all(out$confusion_agg$agg_method == "cluster_weight_jsd_050"))
  expect_true(all(out$confusion_agg$explore_method == "aggregation"))
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

test_that("per-r2 plot filtering keeps only matching annotation-dependent rows", {
  plot_tbl <- tibble::tibble(
    spec_name = c("A-F", "B-F", "B-F", "baseline-single",
                  "truth-warm", "truth-warm", "truth-warm"),
    use_case_id = c("susine_vanilla", "susine_eb_clamped_scale_var_nonneg",
                    "susine_eb_clamped_scale_var_nonneg", "susine_vanilla",
                    "susine_vanilla", "susine_eb_clamped_scale_var_nonneg",
                    "susine_eb_clamped_scale_var_nonneg"),
    annotation_r2 = c(NA_real_, NA_real_, 0.5, NA_real_, NA_real_, NA_real_, 0.5),
    agg_method = c("uniform", "uniform", "uniform", NA_character_,
                   NA_character_, NA_character_, NA_character_),
    AUPRC = c(0.11, 0.21, 0.31, 0.09, 0.15, 0.25, 0.35),
    series_type = c("aggregated", "aggregated", "aggregated",
                    "baseline_single", "truth_warm", "truth_warm", "truth_warm")
  )

  out <- test_susine:::filter_shared_plot_metric_rows(plot_tbl, r2_filter = 0.5)

  expect_equal(sort(unique(out$spec_name)), c("A-F", "B-F"))
  expect_false(any(out$spec_name == "baseline-single"))
  expect_false(any(out$spec_name == "truth-warm"))
  expect_equal(out$AUPRC[out$spec_name == "A-F"], 0.11, tolerance = 1e-12)
  expect_equal(out$AUPRC[out$spec_name == "B-F"], 0.31, tolerance = 1e-12)
  expect_false(any(out$spec_name == "B-F" & is.na(out$annotation_r2)))
})

test_that("per-r2 truth-warm filtering keeps vanilla plus the matching annotation row", {
  plot_tbl <- tibble::tibble(
    spec_name = c("truth-warm", "truth-warm", "truth-warm"),
    use_case_id = c("susine_vanilla", "susine_eb_clamped_scale_var_nonneg",
                    "susine_eb_clamped_scale_var_nonneg"),
    annotation_r2 = c(NA_real_, NA_real_, 0.5),
    agg_method = c(NA_character_, NA_character_, NA_character_),
    AUPRC = c(0.15, 0.25, 0.35),
    series_type = c("truth_warm", "truth_warm", "truth_warm")
  )

  out <- test_susine:::filter_truth_warm_reference_rows(plot_tbl, r2_filter = 0.5)

  expect_equal(sort(out$use_case_id), c("susine_eb_clamped_scale_var_nonneg", "susine_vanilla"))
  expect_true(any(out$use_case_id == "susine_vanilla" & is.na(out$annotation_r2)))
  expect_true(any(out$use_case_id == "susine_eb_clamped_scale_var_nonneg" &
                    out$annotation_r2 == 0.5))
  expect_false(any(out$use_case_id == "susine_eb_clamped_scale_var_nonneg" &
                     is.na(out$annotation_r2)))
})

test_that("truth-warm reference filtering also accepts extracted rows without series_type", {
  refs_tbl <- tibble::tibble(
    use_case_id = c("susine_vanilla", "susine_eb_clamped_scale_var_nonneg",
                    "susine_eb_clamped_scale_var_nonneg"),
    annotation_r2 = c(NA_real_, NA_real_, 0.5),
    AUPRC = c(0.15, 0.25, 0.35)
  )

  overall <- test_susine:::filter_truth_warm_reference_rows(refs_tbl)
  by_r2 <- test_susine:::filter_truth_warm_reference_rows(refs_tbl, r2_filter = 0.5)

  expect_equal(sort(overall$use_case_id), c("susine_eb_clamped_scale_var_nonneg", "susine_vanilla"))
  expect_equal(sort(by_r2$use_case_id), c("susine_eb_clamped_scale_var_nonneg", "susine_vanilla"))
  expect_true(any(by_r2$use_case_id == "susine_eb_clamped_scale_var_nonneg" &
                    by_r2$annotation_r2 == 0.5))
})

test_that("pure refine scaling uses planned thresholds even when BFS truncates the tree", {
  out <- test_susine:::compute_scaling_confusion_bins_for_group(
    pip_list = list(
      c(0.9, 0.1, 0.05),
      c(0.8, 0.2, 0.05),
      c(0.7, 0.3, 0.05)
    ),
    elbo_vec = c(3, 2, 1),
    meta_list = list(
      list(run_id = 1L, spec_name = "A-F", annotation_r2 = NA_real_, refine_step = 1L),
      list(run_id = 2L, spec_name = "A-F", annotation_r2 = NA_real_, refine_step = 2L),
      list(run_id = 3L, spec_name = "A-F", annotation_r2 = NA_real_, refine_step = 3L)
    ),
    causal_vec = c(1L, 0L, 0L),
    group_run_row = tibble::tibble(
      exploration_methods = "refine",
      .planned_n_refine = 8L
    ),
    n_ens_sizes = c(4L, 8L),
    pip_breaks = c(0, 0.05, 0.1, 0.2, 0.5, 1),
    dataset_bundle_id = 1L
  )

  expect_false(is.null(out))
  expect_setequal(unique(out$n_ensemble), c(4L, 8L))
})

test_that("collect backfills missing pure-refine scaling thresholds from full aggregated bins", {
  scaling_bins_raw <- tibble::tibble(
    spec_name = character(),
    annotation_r2 = numeric(),
    dataset_bundle_id = integer(),
    n_ensemble = integer(),
    resolution = character(),
    agg_method = character(),
    rep = integer(),
    pip_threshold = numeric(),
    n_causal_at_bucket = integer(),
    n_noncausal_at_bucket = integer()
  )

  confusion_agg <- tibble::tibble(
    dataset_bundle_id = c(1L, 1L),
    spec_name = c("A-F", "A-F"),
    use_case_id = c("susine_vanilla", "susine_vanilla"),
    annotation_r2 = c(NA_real_, NA_real_),
    group_key = c("g1", "g1"),
    explore_method = c("aggregation", "aggregation"),
    agg_method = c("uniform", "uniform"),
    pip_threshold = c(0.9, 0.2),
    n_causal_at_bucket = c(1L, 0L),
    n_noncausal_at_bucket = c(0L, 2L)
  )

  model_metrics <- tibble::tibble(
    dataset_bundle_id = c(1L, 1L, 1L),
    spec_name = c("A-F", "A-F", "A-F"),
    use_case_id = c("susine_vanilla", "susine_vanilla", "susine_vanilla"),
    annotation_r2 = c(NA_real_, NA_real_, NA_real_),
    group_key = c("g1", "g1", "g1"),
    explore_method = c("refine", "refine", "refine"),
    agg_method = c(NA_character_, NA_character_, NA_character_),
    run_id = c(11L, 12L, 13L)
  )

  run_info <- tibble::tibble(
    dataset_bundle_id = rep(1L, 8),
    spec_name = rep("A-F", 8),
    use_case_id = rep("susine_vanilla", 8),
    annotation_r2 = rep(NA_real_, 8),
    group_key = rep("g1", 8),
    exploration_methods = rep("refine", 8),
    refine_step = 1:8
  )

  out <- test_susine:::backfill_pure_refine_scaling_bins(
    scaling_bins_raw = scaling_bins_raw,
    confusion_agg = confusion_agg,
    model_metrics = model_metrics,
    run_info = run_info,
    scaling_n_ens_sizes = c(4L, 8L)
  )

  expect_equal(sort(unique(out$scaling_bins_raw$n_ensemble)), c(4L, 8L))
  expect_equal(nrow(out$repaired_groups), 2L)
  expect_equal(out$repair_summary$n_backfilled_groups[[1]], 2L)
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

test_that("spec baseline mapping follows A/B/C model families", {
  expect_equal(
    test_susine:::spec_single_fit_use_case(c("A-F", "B-R", "C-CSR", "truth-warm")),
    c("susine_vanilla", "susine_eb_clamped_scale_var_nonneg", "susine_functional_mu", NA_character_)
  )
})

test_that("scaling summaries gain exact n=1 single-fit anchors per spec family", {
  summary_df <- tibble::tibble(
    spec_name = c("A-F", "B-F", "C-C"),
    annotation_r2 = c(NA_real_, 0.5, 0.5),
    n_ensemble = c(4L, 4L, 4L),
    resolution = c(NA_character_, NA_character_, NA_character_),
    agg_method = c("uniform", "uniform", "cluster_weight"),
    AUPRC_mean = c(0.20, 0.30, 0.31),
    AUPRC_se = c(0.01, 0.01, 0.01),
    tpr05_mean = c(0.40, 0.50, 0.51),
    tpr05_se = c(0.01, 0.01, 0.01),
    n_reps = c(5L, 5L, 5L)
  )

  auprc_plot_tbl <- tibble::tibble(
    spec_name = "baseline-single",
    use_case_id = c("susine_vanilla", "susine_eb_clamped_scale_var_nonneg",
                    "susine_functional_mu", "susine_functional_mu"),
    annotation_r2 = c(NA_real_, 0.5, NA_real_, 0.5),
    agg_method = NA_character_,
    AUPRC = c(0.17, 0.27, 0.18, 0.28),
    series_type = "baseline_single"
  )

  tpr05_plot_tbl <- tibble::tibble(
    spec_name = "baseline-single",
    use_case_id = c("susine_vanilla", "susine_eb_clamped_scale_var_nonneg",
                    "susine_functional_mu", "susine_functional_mu"),
    annotation_r2 = c(NA_real_, 0.5, NA_real_, 0.5),
    agg_method = NA_character_,
    tpr_fpr05 = c(0.37, 0.47, 0.38, 0.48),
    series_type = "baseline_single"
  )

  out <- test_susine:::append_scaling_single_fit_points(
    summary_df,
    auprc_plot_table = auprc_plot_tbl,
    tpr05_plot_table = tpr05_plot_tbl
  )

  expect_true(any(out$spec_name == "A-F" & out$n_ensemble == 1L &
                    is.na(out$annotation_r2) & abs(out$AUPRC_mean - 0.17) < 1e-12))
  expect_true(any(out$spec_name == "B-F" & out$n_ensemble == 1L &
                    out$annotation_r2 == 0.5 & abs(out$AUPRC_mean - 0.27) < 1e-12))
  expect_true(any(out$spec_name == "C-C" & out$n_ensemble == 1L &
                    out$annotation_r2 == 0.5 & abs(out$AUPRC_mean - 0.28) < 1e-12))
  expect_true(any(out$spec_name == "C-C" & out$n_ensemble == 1L &
                    out$annotation_r2 == 0.5 & abs(out$tpr05_mean - 0.48) < 1e-12))
  expect_true(all(is.na(out$AUPRC_se[out$n_ensemble == 1L])))
  expect_true(all(out$n_reps[out$n_ensemble == 1L] == 1L))
})
