# Quality checks and smoke validation ---------------------------------------

read_csv_if_exists <- function(path) {
  if (!file.exists(path)) {
    return(NULL)
  }
  readr::read_csv(path, show_col_types = FALSE)
}

missing_required_cols <- function(df, required_cols) {
  if (is.null(df)) {
    return(required_cols)
  }
  setdiff(required_cols, names(df))
}

#' Validate required aggregated metric outputs and schemas.
#'
#' This provides a concrete checklist for Section 2.7 metric persistence.
#'
#' @param aggregated_dir Directory produced by aggregate_staging_outputs().
#' @param require_multimodal Logical; enforce multimodal_metrics.csv presence.
#' @return List with `ok` and per-file `checks` tibble.
#' @export
validate_metrics_coverage <- function(aggregated_dir, require_multimodal = TRUE) {
  required <- list(
    model_metrics = c("run_id", "use_case_id", "filtering", "power", "mean_size", "mean_purity", "AUPRC", "cross_entropy", "hg2", "elbo_final"),
    effect_metrics = c("run_id", "use_case_id", "effect", "size", "purity", "coverage"),
    dataset_metrics = c("dataset_bundle_id", "M1", "z_topk_ratio", "z_max_abs", "z_count_abs_gt_3", "z_eff_signals", "high_ld_count_095", "high_ld_count_095_per_snp", "high_ld_frac_095"),
    confusion_bins = c("run_id", "pip_threshold", "n_causal_at_bucket", "n_noncausal_at_bucket"),
    validation = c("run_id", "task_id", "has_issues"),
    multimodal_metrics = c("use_case_id", "group_label", "explore_method", "mean_jsd", "median_jsd", "max_jsd", "jaccard_top10", "mean_pip_var", "n_clusters")
  )

  files <- c(
    model_metrics = file.path(aggregated_dir, "model_metrics.csv"),
    effect_metrics = file.path(aggregated_dir, "effect_metrics.csv"),
    dataset_metrics = file.path(aggregated_dir, "dataset_metrics.csv"),
    confusion_bins = file.path(aggregated_dir, "confusion_bins.csv"),
    validation = file.path(aggregated_dir, "validation.csv"),
    multimodal_metrics = file.path(aggregated_dir, "multimodal_metrics.csv")
  )

  checks <- purrr::map_dfr(names(files), function(k) {
    f <- files[[k]]
    is_required <- if (k == "multimodal_metrics") isTRUE(require_multimodal) else TRUE
    df <- read_csv_if_exists(f)
    missing_cols <- missing_required_cols(df, required[[k]])
    tibble::tibble(
      table = k,
      path = f,
      required = is_required,
      exists = !is.null(df),
      missing_columns = if (length(missing_cols)) paste(missing_cols, collapse = ", ") else NA_character_,
      ok = if (!is_required && is.null(df)) TRUE else (!is.null(df) && length(missing_cols) == 0L)
    )
  })

  list(
    ok = all(checks$ok),
    checks = checks
  )
}

#' Review seed propagation and uniqueness in a job config.
#'
#' @param job_config Job configuration list returned by make_job_config().
#' @return Tibble with pass/fail checks.
#' @export
seed_management_report <- function(job_config) {
  runs <- tibble::as_tibble(job_config$tables$runs)
  bundles <- tibble::as_tibble(job_config$tables$dataset_bundles)

  bundle_seed_ok <- bundles %>%
    dplyr::group_by(.data$dataset_bundle_id) %>%
    dplyr::summarise(n = dplyr::n_distinct(.data$phenotype_seed), .groups = "drop") %>%
    dplyr::summarise(ok = all(.data$n == 1L), detail = paste0("max distinct phenotype_seed per bundle=", max(.data$n))) 

  ann_rows <- runs %>% dplyr::filter(!is.na(.data$annotation_seed))
  ann_expected <- ann_rows %>%
    dplyr::distinct(.data$dataset_bundle_id, .data$annotation_r2, .data$inflate_match) %>%
    nrow()
  ann_observed <- dplyr::n_distinct(ann_rows$annotation_seed, na.rm = TRUE)
  ann_ok <- ann_observed == ann_expected

  restart_rows <- runs %>% dplyr::filter(!is.na(.data$restart_id))
  restart_observed <- dplyr::n_distinct(restart_rows$restart_seed, na.rm = TRUE)
  restart_expected <- nrow(restart_rows)
  restart_ok <- restart_observed == restart_expected

  run_id_ok <- dplyr::n_distinct(runs$run_id) == nrow(runs)

  n_bundles <- nrow(bundles)
  n_unique_seeds <- dplyr::n_distinct(bundles$phenotype_seed)
  seed_globally_unique_ok <- n_unique_seeds == n_bundles

  tibble::tibble(
    check = c(
      "bundle_has_single_phenotype_seed",
      "phenotype_seed_globally_unique",
      "annotation_seed_unique_per_bundle_annotation_setting",
      "restart_seed_unique_per_restart_row",
      "run_id_unique"
    ),
    pass = c(
      isTRUE(bundle_seed_ok$ok[[1]]),
      seed_globally_unique_ok,
      ann_ok,
      restart_ok,
      run_id_ok
    ),
    detail = c(
      as.character(bundle_seed_ok$detail[[1]]),
      paste0("unique phenotype_seeds=", n_unique_seeds, "; bundles=", n_bundles),
      paste0("observed=", ann_observed, "; expected=", ann_expected),
      paste0("observed=", restart_observed, "; expected=", restart_expected),
      paste0("observed unique run_id=", dplyr::n_distinct(runs$run_id), "; rows=", nrow(runs))
    )
  )
}

#' Run an end-to-end two-locus local smoke test.
#'
#' Uses duplicated simulation_n3 matrix metadata (2 matrix_id rows) to enforce
#' a two-bundle run, executes all tasks, aggregates outputs, and validates metric
#' schemas + seed checks.
#'
#' @param job_name Job name for temporary outputs.
#' @param parent_job_id Parent id for local staged output directory.
#' @param output_root Output root.
#' @param repo_root Repository root.
#' @param strict Logical; stop() on failed checks.
#' @return List with config, aggregation, metric checks, and seed checks.
#' @export
run_two_locus_smoke_test <- function(job_name = "tmp_two_locus_smoke",
                                     parent_job_id = "local_test",
                                     output_root = "output",
                                     repo_root = ".",
                                     strict = TRUE) {
  pq <- prior_quality_grid(c(0.3), c(1))
  data_catalog <- tibble::tibble(
    data_scenario = c("simulation_n3", "simulation_n3"),
    dataset_label = c("simulation_n3_A", "simulation_n3_B"),
    participant_count = c(NA_integer_, NA_integer_),
    snps_post = c(NA_integer_, NA_integer_),
    snp_set = c(NA_character_, NA_character_),
    matrix_path = c(NA_character_, NA_character_),
    manifest_path = c(NA_character_, NA_character_),
    source = c("simulation", "simulation"),
    matrix_index = c(1L, 2L),
    matrix_id = c(1L, 2L)
  )

  cfg <- make_job_config(
    job_name = job_name,
    use_case_ids = c("susie_vanilla"),
    exploration_methods = c("restart"),
    exploration_mode = "separate",
    K = 2L,
    L_grid = 5L,
    y_noise_grid = 0.8,
    prior_quality = pq,
    p_star_grid = 2L,
    seeds = 1L,
    architecture_grid = c("sparse"),
    data_scenarios = "simulation_n3",
    repo_root = repo_root,
    data_matrix_catalog = data_catalog,
    task_unit = "dataset",
    bundles_per_task = 1L,
    sigma_0_2_scalars = "0.2",
    restart_settings = list(n_inits = 2L, alpha_concentration = 0.1),
    verbose_file_output = FALSE
  )
  artifacts <- write_job_artifacts(cfg, "inst/scripts/run_task.R")

  old_job <- Sys.getenv("SUSINE_JOB_NAME", unset = "")
  old_parent <- Sys.getenv("SUSINE_PARENT_ID", unset = "")
  on.exit({
    Sys.setenv(SUSINE_JOB_NAME = old_job, SUSINE_PARENT_ID = old_parent)
  }, add = TRUE)
  Sys.setenv(SUSINE_JOB_NAME = job_name, SUSINE_PARENT_ID = parent_job_id)

  task_ids <- sort(unique(cfg$tables$tasks$task_id))
  for (tid in task_ids) {
    run_task(
      job_name = job_name,
      task_id = as.integer(tid),
      job_root = output_root,
      config_path = artifacts$job_config,
      quiet = TRUE
    )
  }

  agg <- aggregate_staging_outputs(
    job_name = job_name,
    parent_job_id = parent_job_id,
    output_root = output_root,
    validate = TRUE
  )
  metric_check <- validate_metrics_coverage(agg$output_dir, require_multimodal = TRUE)
  seed_check <- seed_management_report(cfg)

  if (isTRUE(strict)) {
    if (!isTRUE(metric_check$ok)) {
      stop("Two-locus smoke failed metrics coverage checks.")
    }
    if (!all(seed_check$pass)) {
      stop("Two-locus smoke failed seed management checks.")
    }
  }

  list(
    config = cfg,
    artifacts = artifacts,
    tasks = task_ids,
    aggregation = agg,
    metric_check = metric_check,
    seed_check = seed_check
  )
}
