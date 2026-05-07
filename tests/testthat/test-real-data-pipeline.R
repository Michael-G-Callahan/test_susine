real_data_test_workspace <- function(tag = "real-data-pipeline") {
  path <- file.path(tempdir(), sprintf("%s-%s", tag, as.integer(stats::runif(1, 1, 1e9))))
  dir.create(path, recursive = TRUE, showWarnings = FALSE)
  normalizePath(path, winslash = "/", mustWork = TRUE)
}

write_real_data_test_source <- function(source_root,
                                        locus_id = "testgene_chr1_lung",
                                        gene_name = "TESTGENE") {
  z_dir <- file.path(source_root, "output", "z_score", locus_id)
  ld_dir <- file.path(source_root, "output", "ld", locus_id)
  mu0_dir <- file.path(
    source_root,
    "output", "susine_mu0",
    "representative_gene_sample_n10_annotation_selection",
    "per_locus_annotations"
  )
  mu0_root <- dirname(mu0_dir)
  alpha_dir <- file.path(source_root, "output", "annotation", "alphagenome")
  dir.create(z_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(ld_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(mu0_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(alpha_dir, recursive = TRUE, showWarnings = FALSE)

  variant_ids <- c(
    "chr1_100_A_G_b38",
    "chr1_200_C_T_b38",
    "chr1_300_G_A_b38"
  )
  z_tbl <- tibble::tibble(
    gene_id = "ENSGTEST0001",
    variant_id = variant_ids,
    tss_distance = c(-1000, -900, -800),
    af = c(0.2, 0.3, 0.4),
    ma_samples = c(50, 60, 70),
    ma_count = c(60, 72, 84),
    pval_nominal = c(1e-5, 5e-3, 0.2),
    slope = c(0.20, 0.08, 0.01),
    slope_se = c(0.05, 0.04, 0.05),
    sample_size = c(200, 200, 200),
    snp_pos = c(100, 200, 300),
    valid_tss_centered = c(TRUE, TRUE, TRUE),
    valid_snp_centered = c(TRUE, TRUE, TRUE),
    valid_any = c(TRUE, TRUE, TRUE),
    z_score = c(4.0, 2.0, 0.2)
  )
  master_tbl <- dplyr::mutate(z_tbl, ld_included = TRUE, ld_matrix_index = 0:2)
  order_tbl <- tibble::tibble(
    index = 0:2,
    id = variant_ids,
    chrom = "chr1",
    pos = c(100L, 200L, 300L),
    ref = c("A", "C", "G"),
    alt = c("G", "T", "A")
  )
  annotation_tbl <- tibble::tibble(
    locus_id = locus_id,
    gene_name = gene_name,
    variant_id = variant_ids,
    ld_matrix_index = 0:2,
    annotation_missing = FALSE,
    source_scoring_mode = "unit_test",
    raw_score = c(0.6, 0.2, -0.1),
    quantile_score = c(0.8, 0.4, -0.2),
    q_star = c(0.8, 0.4, -0.2),
    a_raw = c(1.2, 0.5, -0.3),
    a_clip = c(1.2, 0.5, -0.3),
    annotation_a = c(1.2, 0.5, -0.3),
    beta_hat_std = c(0.20, 0.08, 0.01),
    beta_hat_slope = c(0.20, 0.08, 0.01),
    mu0 = c(0.30, 0.125, -0.075),
    var_y_hat_from_slope = c(1.0, 1.0, 1.0),
    var_y_hat_from_se = c(1.0, 1.0, 1.0),
    baseline_c_l = 0.25
  )
  ld_long_tbl <- tibble::tibble(
    snp_index_1 = c(0L, 0L, 1L),
    snp_index_2 = c(1L, 2L, 2L),
    r = c(0.25, 0.10, 0.30)
  )
  mu0_summary_tbl <- tibble::tibble(
    locus_id = locus_id,
    gene_name = gene_name,
    baseline_c_l = 0.25
  )
  selection_tbl <- tibble::tibble(
    locus_id = locus_id,
    gene_name = gene_name,
    selected = TRUE
  )

  readr::write_csv(z_tbl, file.path(z_dir, paste0(gene_name, "_GTEx_z_scores.csv")))
  readr::write_csv(master_tbl, file.path(ld_dir, paste0(gene_name, "_phase1_master_variants.csv")))
  readr::write_tsv(order_tbl, file.path(ld_dir, paste0(gene_name, "_LD_variant_order.tsv")))
  arrow::write_parquet(ld_long_tbl, file.path(ld_dir, paste0(gene_name, "_phase1_LD_R_long.parquet")))
  readr::write_csv(annotation_tbl, file.path(mu0_dir, paste0(gene_name, "_mu0_variant_annotations.csv")))
  readr::write_csv(mu0_summary_tbl, file.path(mu0_root, "mu0_locus_summary.csv"))
  readr::write_csv(
    selection_tbl,
    file.path(alpha_dir, "representative_gene_sample_n10_annotation_selection_summary.csv")
  )

  list(
    locus_id = locus_id,
    gene_name = gene_name,
    variant_ids = variant_ids
  )
}

test_that("real-data input sync copies files and builds a manifest", {
  skip_if_not_installed("arrow")

  workspace <- real_data_test_workspace("real-data-sync")
  source_root <- file.path(workspace, "source")
  dest_root <- file.path(workspace, "study")
  fixture <- write_real_data_test_source(source_root)

  synced <- test_susine::sync_real_data_inputs(
    source_repo_root = source_root,
    dest_root = dest_root
  )

  expect_true(file.exists(file.path(dest_root, "locus_manifest.csv")))
  expect_true(file.exists(file.path(
    dest_root, fixture$locus_id, paste0(fixture$gene_name, "_GTEx_z_scores.csv")
  )))
  expect_true(file.exists(file.path(
    dest_root, fixture$locus_id, paste0(fixture$gene_name, "_phase1_LD_R_long.parquet")
  )))
  expect_equal(nrow(synced$manifest), 1L)
  expect_true(all(c(
    "dataset_bundle_id", "locus_id", "gene_name", "gtex_tissue", "gtex_chrom",
    "n_variants", "n_sample_median", "baseline_c_l", "z_scores_path",
    "ld_r_long_path", "ld_order_path", "master_variants_path", "annotation_path"
  ) %in% names(synced$manifest)))
})

test_that("real-data locus bundle aligns annotation, z scores, and LD order", {
  skip_if_not_installed("arrow")

  workspace <- real_data_test_workspace("real-data-bundle")
  source_root <- file.path(workspace, "source")
  dest_root <- file.path(workspace, "study")
  fixture <- write_real_data_test_source(source_root)

  test_susine::sync_real_data_inputs(
    source_repo_root = source_root,
    dest_root = dest_root
  )
  manifest_path <- file.path(dest_root, "locus_manifest.csv")

  bundle <- test_susine:::load_real_data_locus_bundle(
    locus_id = fixture$locus_id,
    manifest_path = manifest_path,
    repo_root = getwd()
  )

  expect_identical(bundle$variant_map$variant_id, fixture$variant_ids)
  expect_equal(bundle$variant_map$ld_matrix_index, 0:2)
  expect_equal(diag(bundle$R), c(1, 1, 1), tolerance = 1e-12)
  expect_equal(bundle$R, t(bundle$R), tolerance = 1e-12)
  expect_length(bundle$a, 3L)
  expect_length(bundle$z, 3L)
  expect_length(bundle$beta_hat_std, 3L)
})

test_that("real-data job config expands to the expected run grid", {
  skip_if_not_installed("arrow")

  workspace <- real_data_test_workspace("real-data-config")
  source_root <- file.path(workspace, "source")
  dest_root <- file.path(workspace, "study")
  fixture <- write_real_data_test_source(source_root)

  test_susine::sync_real_data_inputs(
    source_repo_root = source_root,
    dest_root = dest_root
  )
  manifest_path <- file.path(dest_root, "locus_manifest.csv")

  full_cfg <- test_susine::build_real_data_job_config(
    job_name = "real-data-grid-test",
    manifest_path = manifest_path,
    output_root = file.path(workspace, "output")
  )
  expect_equal(nrow(full_cfg$tables$runs), 65L)
  expect_equal(sum(full_cfg$tables$runs$run_family == "functional_grid"), 64L)
  expect_equal(sum(full_cfg$tables$runs$run_family == "susie_anchor"), 1L)
  expect_equal(nrow(full_cfg$tables$tasks), 1L)
  expect_equal(unique(full_cfg$tables$tasks$runs_per_task), 65L)
  expect_equal(
    sort(unique(dplyr::filter(full_cfg$tables$runs, .data$run_family == "functional_grid")$flush_group)),
    1:8
  )
  expect_equal(
    unique(dplyr::filter(full_cfg$tables$runs, .data$run_family == "susie_anchor")$flush_group),
    9L
  )

  small_cfg <- test_susine::build_real_data_job_config(
    job_name = "real-data-small-grid-test",
    manifest_path = manifest_path,
    output_root = file.path(workspace, "output-small"),
    c_grid_values = c(0, 0.5),
    sigma_grid_values = c(0.01, 0.2)
  )
  expect_equal(nrow(small_cfg$tables$runs), 5L)
  expect_equal(unique(small_cfg$tables$tasks$runs_per_task), 5L)
  expect_equal(
    unique(dplyr::filter(small_cfg$tables$runs, .data$run_family == "susie_anchor")$flush_group),
    3L
  )
})

test_that("real-data task runner and collector produce truth-free ensemble outputs", {
  skip_if_not_installed("arrow")
  skip_if_not_installed("susieR")

  workspace <- real_data_test_workspace("real-data-run")
  source_root <- file.path(workspace, "source")
  dest_root <- file.path(workspace, "study")
  output_root <- file.path(workspace, "output")
  fixture <- write_real_data_test_source(source_root)

  test_susine::sync_real_data_inputs(
    source_repo_root = source_root,
    dest_root = dest_root
  )
  manifest_path <- file.path(dest_root, "locus_manifest.csv")

  job_cfg <- test_susine::build_real_data_job_config(
    job_name = "real-data-smoke",
    manifest_path = manifest_path,
    output_root = output_root,
    L = 2L,
    c_grid_values = c(0, 0.5),
    sigma_grid_values = c(0.01, 0.2),
    max_iter = 20L,
    tol = 1e-3
  )
  artifacts <- test_susine::write_real_data_job_artifacts(
    job_cfg,
    run_task_script = file.path(
      test_susine:::ensure_repo_root(getwd()),
      "inst", "scripts", "run_real_data_task.R"
    )
  )

  old_job_name <- Sys.getenv("SUSINE_JOB_NAME", unset = "")
  old_parent_id <- Sys.getenv("SUSINE_PARENT_ID", unset = "")
  old_dev <- Sys.getenv("SUSINE_DEV", unset = "")
  on.exit({
    do.call(Sys.setenv, list(
      SUSINE_JOB_NAME = old_job_name,
      SUSINE_PARENT_ID = old_parent_id,
      SUSINE_DEV = old_dev
    ))
  }, add = TRUE)
  Sys.setenv(
    SUSINE_JOB_NAME = job_cfg$job$name,
    SUSINE_PARENT_ID = "test-parent-001",
    SUSINE_DEV = "1"
  )

  test_susine::run_real_data_task(
    job_name = job_cfg$job$name,
    task_id = 1L,
    config_path = artifacts$job_config
  )

  task_dir <- file.path(
    output_root, "slurm_output", job_cfg$job$name, "test-parent-001", "task-001"
  )
  expect_true(file.exists(file.path(task_dir, "flush-001_run_metrics.csv")))
  expect_true(file.exists(file.path(task_dir, "flush-001_variant_posteriors.parquet")))
  expect_true(file.exists(file.path(task_dir, "flush-001_effect_posteriors.parquet")))
  expect_true(file.exists(file.path(task_dir, "flush-001_credible_set_membership.parquet")))
  expect_true(file.exists(file.path(task_dir, "flush-001_dataset_metrics.csv")))
  expect_true(file.exists(file.path(task_dir, "flush-001_fit_file_index.csv")))

  collected <- test_susine::collect_real_data_results(
    job_name = job_cfg$job$name,
    parent_job_id = "test-parent-001",
    output_root = output_root,
    validate = TRUE
  )
  validation <- test_susine::validate_real_data_outputs(collected$output_dir)
  expect_true(isTRUE(validation$ok))

  agg_weights <- readr::read_csv(
    file.path(collected$output_dir, "aggregation_weights_cluster_weight.csv"),
    show_col_types = FALSE
  )
  expect_equal(sum(agg_weights$agg_weight_run), 1, tolerance = 1e-8)

  comparisons <- readr::read_csv(
    file.path(collected$output_dir, "run_comparisons.csv"),
    show_col_types = FALSE
  )
  c0_rows <- dplyr::filter(comparisons, .data$comparison_target == "same_sigma_c0") %>%
    dplyr::left_join(
      readr::read_csv(file.path(collected$output_dir, "run_manifest.csv"), show_col_types = FALSE) %>%
        dplyr::select("run_id", "c_value"),
      by = "run_id"
    ) %>%
    dplyr::filter(abs(.data$c_value) < 1e-12)
  expect_true(nrow(c0_rows) >= 1L)
  expect_equal(c0_rows$jsd_pip, rep(0, nrow(c0_rows)), tolerance = 1e-12)

  schema_tables <- list(
    readr::read_csv(file.path(collected$output_dir, "run_metrics_full.csv"), show_col_types = FALSE),
    readr::read_csv(file.path(collected$output_dir, "dataset_metrics.csv"), show_col_types = FALSE),
    readr::read_csv(file.path(collected$output_dir, "functional_grid_summary.csv"), show_col_types = FALSE)
  )
  truth_only_cols <- c("AUPRC", "power", "coverage", "causal", "cs_power", "cs_fdr")
  expect_false(any(vapply(schema_tables, function(tbl) any(truth_only_cols %in% names(tbl)), logical(1))))

  aggregated_pips <- arrow::read_parquet(
    file.path(collected$output_dir, "aggregated_variant_pips_cluster_weight.parquet")
  ) %>%
    tibble::as_tibble()
  expect_equal(nrow(aggregated_pips), 3L)
  expect_true(all(aggregated_pips$locus_id == fixture$locus_id))
})
