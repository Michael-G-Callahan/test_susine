# Regression test for the execution-cache write-poisoning bug.
#
# A non-root refinement refit is produced with an `init_alpha_override`
# (a perturbed warm start) but is recorded with `blocked_idx = integer(0)` and a
# run_row whose run_type resolves to "default". Its execution-cache key is
# therefore identical to a cold susine_vanilla fit's key. The read guard refuses
# to *serve* such a fit from cache, but the write guard must also refuse to
# *store* it — otherwise a later cold request (e.g. `baseline-single`, or the
# sigma base member of a grid spec) gets a cache hit and receives the REFINED
# fit instead of a genuine cold fit. See R/run_model.R (cache write guard).

make_poison_catalog <- function() {
  tibble::tibble(
    data_scenario = "simulation_n3",
    dataset_label = "simulation_n3_1",
    participant_count = NA_integer_,
    snps_post = NA_integer_,
    snp_set = NA_character_,
    matrix_path = NA_character_,
    manifest_path = NA_character_,
    source = "simulation",
    matrix_index = 1L,
    matrix_id = 1L
  )
}

make_poison_buffer_ctx <- function() {
  buffer_ctx <- new.env(parent = emptyenv())
  buffer_ctx$buffers <- list(
    model = list(), effect = list(), effect_unfiltered = list(),
    tier_cs_metrics = list(), snps = list(), confusion = list(),
    dataset_metrics = list(), multimodal = list(), refine_depth = list(),
    scaling_bins = list(), validation = list(), prior_diagnostics = list(),
    hg2_by_agg = list()
  )
  buffer_ctx$base_output <- tempdir()
  buffer_ctx$task_id <- 1L
  buffer_ctx$flush_index <- 0L
  buffer_ctx
}

make_poison_cfg <- function() {
  test_susine::make_job_config(
    job_name = "unit_refine_cache_poison",
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
    sigma_0_2_scalars = "0.2",
    data_scenarios = "simulation_n3",
    repo_root = ".",
    data_matrix_catalog = make_poison_catalog(),
    task_unit = "dataset",
    bundles_per_task = 1L,
    runs_per_task = 50L,
    output_root = tempdir(),
    write_snps_parquet = FALSE,
    write_confusion_bins = FALSE,
    verbose_file_output = FALSE,
    include_overall_pool = FALSE,
    aggregation_methods = "uniform",
    # Refine spec listed FIRST so its refits run before the cold baseline-single
    # request and can poison the shared cold cache key. purity_threshold = 0 so
    # the BFS reliably branches into non-root (override-init) refits.
    exploration_specs = list(
      list(
        name = "A-F",
        use_case_ids = "susine_vanilla",
        exploration_methods = "refine",
        exploration_mode = "separate",
        K = 4L,
        refine_settings = list(
          n_steps = 4L, cs_source = "filtered", purity_threshold = 0
        )
      ),
      # A small sigma grid spec so the runs table carries a sigma_0_2_scalar
      # column (mirrors the real ensemble's A-S; its sigma=0.2 base member is
      # also a cold cache victim).
      list(
        name = "A-S",
        use_case_ids = "susine_vanilla",
        exploration_methods = "sigma_0_2_grid",
        exploration_mode = "separate",
        K = 1L,
        sigma_0_2_grid_values = 0.2
      ),
      list(
        name = "baseline-single",
        use_case_ids = "susine_vanilla",
        exploration_methods = "single",
        exploration_mode = "separate",
        K = 1L
      )
    )
  )
}

test_that("non-root refine refits do not poison the cold baseline cache key", {
  cfg <- make_poison_cfg()
  runs <- cfg$tables$runs

  # Sanity: all specs share a single dataset bundle.
  expect_setequal(unique(runs$spec_name), c("A-F", "A-S", "baseline-single"))
  expect_equal(dplyr::n_distinct(runs$dataset_bundle_id), 1L)

  # Independently compute the genuine COLD susine_vanilla fit for this bundle.
  # Build the data bundle exactly as execute_dataset_bundle does (from the
  # bundle's annotation settings, not the vanilla run's NA), so X,y match.
  baseline_row <- runs %>%
    dplyr::filter(.data$spec_name == "baseline-single") %>%
    dplyr::slice(1)
  bundle_id <- unique(runs$dataset_bundle_id)
  bundle_row <- dplyr::filter(cfg$tables$dataset_bundles,
                              .data$dataset_bundle_id == !!bundle_id)
  data_bundle <- test_susine:::generate_data_for_bundle(bundle_row, cfg)
  if (is.null(data_bundle$priors_cache)) {
    data_bundle$priors_cache <- new.env(parent = emptyenv())
  }
  use_case <- test_susine:::lookup_use_case(cfg, baseline_row$use_case_id)
  cold <- test_susine:::run_use_case(
    use_case = use_case, run_row = baseline_row,
    data_bundle = data_bundle, job_config = cfg
  )
  cold_elbo <- tail(cold$fits[[1]]$model_fit$elbo, 1)

  # Run the full bundle (refine + baseline-single share one execution cache).
  buffer_ctx <- make_poison_buffer_ctx()
  test_susine:::execute_dataset_bundle(
    bundle_runs = runs, job_config = cfg, quiet = TRUE, buffer_ctx = buffer_ctx
  )
  model_tbl <- dplyr::bind_rows(buffer_ctx$buffers$model)

  # The raw model buffer keys rows by run_id / explore_method (spec_name is
  # joined in only at consolidation), so identify rows that way.

  # Scenario validity: refine must have produced non-root (override-init)
  # refits, i.e. the cache-write path for poisoning was actually exercised.
  nonroot <- model_tbl %>%
    dplyr::filter(.data$explore_method == "refine", .data$variant_id >= 2L)
  expect_gt(nrow(nonroot), 0L)

  # Core invariant: the recorded baseline-single fit is the COLD fit, not a
  # refined refit served from a poisoned cache key.
  bs_elbo <- model_tbl %>%
    dplyr::filter(.data$run_id == !!baseline_row$run_id) %>%
    dplyr::pull(.data$elbo_final) %>%
    unique()
  expect_equal(length(bs_elbo), 1L)
  expect_equal(bs_elbo, cold_elbo, tolerance = 1e-8)
})
