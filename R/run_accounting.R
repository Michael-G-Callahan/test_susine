# Run accounting helpers ----------------------------------------------------

#' Build human-readable run accounting lines for a job config.
#'
#' @param job_config Job configuration list produced by make_job_config().
#' @param seeds_per_matrix Integer count of phenotype seeds per matrix.
#' @param use_case_ids Optional vector of use_case_id values to include.
#' @return Character vector of lines suitable for cat().
#' @export
run_accounting_lines <- function(job_config,
                                 seeds_per_matrix,
                                 use_case_ids = NULL) {
  runs_tbl <- job_config$tables$runs
  bundles_tbl <- job_config$tables$dataset_bundles

  n_bundles <- nrow(bundles_tbl)
  n_matrices <- dplyr::n_distinct(bundles_tbl$matrix_id)
  n_p_star <- dplyr::n_distinct(bundles_tbl$p_star)
  n_y_noise <- dplyr::n_distinct(bundles_tbl$y_noise)
  n_seeds_per_matrix <- as.integer(seeds_per_matrix %||% 0L)
  n_seed_values_total <- dplyr::n_distinct(bundles_tbl$phenotype_seed, na.rm = TRUE)

  bundles_formula <- sprintf(
    "(= %d loci x %d p_star x %d y_noise x %d seeds per matrix)",
    n_matrices, n_p_star, n_y_noise, n_seeds_per_matrix
  )

  n_tasks <- nrow(job_config$tables$tasks)
  bundles_per_task <- job_config$job$bundles_per_task %||% 1L
  tasks_formula <- sprintf("(= ceiling(%d / %d))", n_bundles, bundles_per_task)

  use_cases_tbl <- job_config$tables$use_cases
  if (!is.null(use_case_ids)) {
    use_cases_tbl <- dplyr::filter(use_cases_tbl, use_case_id %in% use_case_ids)
    runs_tbl <- dplyr::filter(runs_tbl, use_case_id %in% use_case_ids)
  }
  n_use_cases <- nrow(use_cases_tbl)

  count_or_one <- function(x) {
    val <- dplyr::n_distinct(x[!is.na(x)])
    ifelse(val < 1, 1L, val)
  }

  required_cols <- c("annotation_r2", "inflate_match", "sigma_0_2_scalar", "c_value", "tau_value", "restart_id")
  for (nm in required_cols) {
    if (!nm %in% names(runs_tbl)) {
      runs_tbl[[nm]] <- NA
    }
  }

  uc_counts <- runs_tbl %>%
    dplyr::group_by(use_case_id) %>%
    dplyr::summarise(
      n_L = dplyr::n_distinct(L),
      n_ann_r2 = count_or_one(annotation_r2),
      n_inflate = count_or_one(inflate_match),
      n_sigma = count_or_one(sigma_0_2_scalar),
      n_c = count_or_one(c_value),
      n_tau = count_or_one(tau_value),
      n_restart = count_or_one(restart_id),
      .groups = "drop"
    )

  group_counts <- runs_tbl %>%
    dplyr::group_by(use_case_id, dataset_bundle_id) %>%
    dplyr::summarise(n_groups = dplyr::n_distinct(group_key), .groups = "drop") %>%
    dplyr::group_by(use_case_id) %>%
    dplyr::summarise(min_groups = min(n_groups), max_groups = max(n_groups), .groups = "drop")

  run_counts <- runs_tbl %>%
    dplyr::group_by(use_case_id, dataset_bundle_id, group_key) %>%
    dplyr::summarise(n_runs = dplyr::n(), .groups = "drop") %>%
    dplyr::group_by(use_case_id) %>%
    dplyr::summarise(min_runs = min(n_runs), max_runs = max(n_runs), .groups = "drop")

  uc_counts <- uc_counts %>%
    dplyr::left_join(group_counts, by = "use_case_id") %>%
    dplyr::left_join(run_counts, by = "use_case_id")

  lines <- c(
    sprintf("Number of dataset bundles: %d %s", n_bundles, bundles_formula),
    sprintf("Distinct phenotype_seed values: %d", n_seed_values_total)
  )

  restart_seeds_tbl <- runs_tbl %>% dplyr::filter(!is.na(restart_id))
  n_restart_seeds_total <- dplyr::n_distinct(restart_seeds_tbl$restart_seed, na.rm = TRUE)
  n_restart_ids <- if (nrow(restart_seeds_tbl)) dplyr::n_distinct(restart_seeds_tbl$restart_id) else 0L
  if (n_restart_seeds_total > 0) {
    expected_restart <- nrow(restart_seeds_tbl)
    lines <- c(lines, sprintf("Distinct restart_seed values: %d (= %d restart runs)",
                              n_restart_seeds_total, expected_restart))
    if (n_restart_seeds_total != expected_restart) {
      lines <- c(lines, sprintf("  WARNING: restart_seed distinct count mismatch (%d vs %d)",
                                n_restart_seeds_total, expected_restart))
    }
  } else {
    lines <- c(lines, "Distinct restart_seed values: 0 (no restart rows in runs table)")
  }

  annotation_seeds_tbl <- runs_tbl %>% dplyr::filter(!is.na(annotation_seed))
  n_annotation_seeds_total <- dplyr::n_distinct(annotation_seeds_tbl$annotation_seed, na.rm = TRUE)
  n_ann_r2_vals <- if (nrow(annotation_seeds_tbl)) dplyr::n_distinct(annotation_seeds_tbl$annotation_r2) else 0L
  n_inflate_vals <- if (nrow(annotation_seeds_tbl)) dplyr::n_distinct(annotation_seeds_tbl$inflate_match) else 0L
  if (n_annotation_seeds_total > 0) {
    expected_annotation <- annotation_seeds_tbl %>%
      dplyr::distinct(dataset_bundle_id, annotation_r2, inflate_match) %>%
      nrow()
    lines <- c(lines, sprintf("Distinct annotation_seed values: %d (= %d bundle x annotation settings)",
                              n_annotation_seeds_total, expected_annotation))
    if (n_annotation_seeds_total != expected_annotation) {
      lines <- c(lines, sprintf("  WARNING: annotation_seed distinct count mismatch (%d vs %d)",
                                n_annotation_seeds_total, expected_annotation))
    }
  } else {
    lines <- c(lines, "Distinct annotation_seed values: 0 (no annotation rows in runs table)")
  }

  lines <- c(
    lines,
    sprintf("Number of tasks: %d %s", n_tasks, tasks_formula),
    sprintf("# of use cases: %d", n_use_cases)
  )

  if (nrow(uc_counts)) {
    for (i in seq_len(nrow(uc_counts))) {
      row <- uc_counts[i, ]
      uc <- row$use_case_id
      groups_line <- sprintf(
        "Groups per bundle (%s): observed %d-%d",
        uc, row$min_groups, row$max_groups
      )
      lines <- c(lines, groups_line)

      runs_line <- sprintf(
        "Runs per group (%s): observed %d-%d (= %d c-values x %d tau-values x %d sigma-values x %d restarts)",
        uc, row$min_runs, row$max_runs, row$n_c, row$n_tau, row$n_sigma, row$n_restart
      )
      lines <- c(lines, runs_line)
    }
  }

  total_runs <- nrow(runs_tbl)
  per_uc_bundle <- runs_tbl %>%
    dplyr::group_by(use_case_id, dataset_bundle_id) %>%
    dplyr::summarise(n_runs = dplyr::n(), .groups = "drop") %>%
    dplyr::group_by(use_case_id) %>%
    dplyr::summarise(
      min_bundle_runs = min(n_runs),
      max_bundle_runs = max(n_runs),
      .groups = "drop"
    )
  total_runs_per_bundle <- runs_tbl %>%
    dplyr::group_by(dataset_bundle_id) %>%
    dplyr::summarise(n_runs = dplyr::n(), .groups = "drop") %>%
    dplyr::summarise(total = dplyr::first(n_runs), .groups = "drop") %>%
    dplyr::pull(total)
  if (!length(total_runs_per_bundle)) {
    total_runs_per_bundle <- 0L
  }
  per_uc_bundle_terms <- per_uc_bundle %>%
    dplyr::mutate(term = sprintf("%s: %d-%d per bundle", use_case_id, min_bundle_runs, max_bundle_runs)) %>%
    dplyr::pull(term)
  lines <- c(lines, sprintf(
    "Total runs per bundle: %d (observed; %s)",
    total_runs_per_bundle,
    paste(per_uc_bundle_terms, collapse = " + ")
  ))
  per_uc_terms <- runs_tbl %>%
    dplyr::group_by(use_case_id) %>%
    dplyr::summarise(total = dplyr::n(), .groups = "drop") %>%
    dplyr::mutate(term = sprintf("%s: %d runs", use_case_id, total)) %>%
    pull(term)
  lines <- c(lines, sprintf(
    "Total runs: %d (= %s)",
    total_runs,
    paste(per_uc_terms, collapse = " + ")
  ))

  lines
}
