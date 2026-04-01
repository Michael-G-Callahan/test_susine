# Internal helpers for ensemble-scaling metric standardization.

#' @keywords internal
normalize_agg_method_id <- function(x) {
  x <- as.character(x)
  x[x %in% c("cluster_weight_050", "cluster_weight_jsd_050")] <- "cluster_weight_jsd_050"
  x
}

#' @keywords internal
normalize_agg_method_df <- function(df, col = "agg_method") {
  if (is.null(df) || !nrow(df) || !col %in% names(df)) return(df)
  df[[col]] <- normalize_agg_method_id(df[[col]])
  df
}

#' @keywords internal
append_overall_scaling_bins <- function(scaling_bins) {
  scaling_bins <- normalize_agg_method_df(scaling_bins)
  if (is.null(scaling_bins) || !nrow(scaling_bins) ||
      !"annotation_r2" %in% names(scaling_bins)) {
    return(scaling_bins)
  }

  has_r2 <- scaling_bins %>% dplyr::filter(!is.na(.data$annotation_r2))
  if (!nrow(has_r2)) return(scaling_bins)

  overall_bins <- has_r2 %>%
    dplyr::group_by(
      .data$spec_name,
      .data$n_ensemble,
      .data$resolution,
      .data$agg_method,
      .data$rep,
      .data$pip_threshold
    ) %>%
    dplyr::summarise(
      n_causal_at_bucket = sum(.data$n_causal_at_bucket, na.rm = TRUE),
      n_noncausal_at_bucket = sum(.data$n_noncausal_at_bucket, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    dplyr::mutate(annotation_r2 = NA_real_)

  if ("n_datasets" %in% names(scaling_bins)) {
    overall_bins <- overall_bins %>%
      dplyr::left_join(
        has_r2 %>%
          dplyr::group_by(
            .data$spec_name,
            .data$n_ensemble,
            .data$resolution,
            .data$agg_method,
            .data$rep,
            .data$pip_threshold
          ) %>%
          dplyr::summarise(n_datasets = sum(.data$n_datasets, na.rm = TRUE), .groups = "drop"),
        by = c("spec_name", "n_ensemble", "resolution", "agg_method", "rep", "pip_threshold")
      )
  }

  dplyr::bind_rows(scaling_bins, overall_bins)
}

#' @keywords internal
collapse_plot_metric_table <- function(df, group_vars, metric_col, label) {
  if (is.null(df) || !nrow(df)) return(tibble::tibble())

  metric_sym <- rlang::sym(metric_col)
  raw <- df %>%
    dplyr::filter(!is.na(.data$spec_name))

  dup_summary <- raw %>%
    dplyr::group_by(dplyr::across(dplyr::all_of(group_vars))) %>%
    dplyr::summarise(
      n_rows = dplyr::n(),
      n_metric = sum(is.finite(!!metric_sym)),
      metric_min = suppressWarnings(min(!!metric_sym, na.rm = TRUE)),
      metric_max = suppressWarnings(max(!!metric_sym, na.rm = TRUE)),
      .groups = "drop"
    ) %>%
    dplyr::mutate(
      metric_span = dplyr::if_else(
        .data$n_metric > 0L,
        .data$metric_max - .data$metric_min,
        0
      )
    )

  conflicting <- dup_summary %>%
    dplyr::filter(.data$n_rows > 1L, is.finite(.data$metric_span), .data$metric_span > 1e-12)

  if (nrow(conflicting) > 0L) {
    message(
      sprintf(
        "%s had %d duplicate plot-key rows with non-identical %s values; collapsing by mean.",
        label, nrow(conflicting), metric_col
      )
    )
  }

  raw %>%
    dplyr::group_by(dplyr::across(dplyr::all_of(group_vars))) %>%
    dplyr::summarise(
      !!metric_col := mean(!!metric_sym, na.rm = TRUE),
      .groups = "drop"
    )
}

#' @keywords internal
build_shared_plot_metric_table <- function(agg_overall,
                                           agg_by_r2,
                                           individual_overall,
                                           individual_by_r2,
                                           metric_col,
                                           label) {
  metric_group_vars <- c("spec_name", "use_case_id", "annotation_r2", "agg_method", "series_type")
  keep_cols <- c("spec_name", "use_case_id", "annotation_r2", "agg_method", metric_col, "series_type")

  aggregated_overall <- normalize_agg_method_df(agg_overall) %>%
    dplyr::mutate(annotation_r2 = NA_real_, series_type = "aggregated") %>%
    dplyr::select(dplyr::any_of(keep_cols))
  aggregated_by_r2 <- normalize_agg_method_df(agg_by_r2) %>%
    dplyr::mutate(series_type = "aggregated") %>%
    dplyr::select(dplyr::any_of(keep_cols))
  individual_overall <- individual_overall %>%
    dplyr::filter(.data$spec_name %in% c("baseline-single", "truth-warm")) %>%
    dplyr::mutate(
      annotation_r2 = NA_real_,
      agg_method = NA_character_,
      series_type = dplyr::if_else(
        .data$spec_name == "baseline-single",
        "baseline_single",
        "truth_warm"
      )
    ) %>%
    dplyr::select(dplyr::any_of(keep_cols))
  individual_by_r2 <- individual_by_r2 %>%
    dplyr::filter(.data$spec_name %in% c("baseline-single", "truth-warm")) %>%
    dplyr::mutate(
      agg_method = NA_character_,
      series_type = dplyr::if_else(
        .data$spec_name == "baseline-single",
        "baseline_single",
        "truth_warm"
      )
    ) %>%
    dplyr::select(dplyr::any_of(keep_cols))

  collapse_plot_metric_table(
    dplyr::bind_rows(
      aggregated_overall,
      aggregated_by_r2,
      individual_overall,
      individual_by_r2
    ),
    group_vars = metric_group_vars,
    metric_col = metric_col,
    label = label
  )
}

#' @keywords internal
backfill_single_fit_refine_agg_confusion <- function(confusion_individual,
                                                     confusion_agg,
                                                     run_info,
                                                     agg_methods) {
  agg_methods <- unique(normalize_agg_method_id(agg_methods))
  confusion_agg <- normalize_agg_method_df(confusion_agg)

  empty_result <- list(
    confusion_agg = confusion_agg,
    repaired_groups = tibble::tibble(),
    repair_summary = tibble::tibble()
  )

  if (is.null(confusion_individual) || !nrow(confusion_individual) ||
      is.null(run_info) || !nrow(run_info) || !length(agg_methods)) {
    return(empty_result)
  }

  group_cols <- c("spec_name", "use_case_id", "dataset_bundle_id", "annotation_r2", "group_key")
  task_lookup <- run_info %>%
    dplyr::filter(!is.na(.data$task_id)) %>%
    dplyr::distinct(dplyr::across(dplyr::all_of(c(group_cols, "task_id"))))

  refine_groups <- run_info %>%
    dplyr::filter(
      !is.na(.data$spec_name),
      !.data$spec_name %in% c("baseline-single", "truth-warm"),
      !is.na(.data$use_case_id),
      !is.na(.data$dataset_bundle_id),
      !is.na(.data$group_key),
      as.character(.data$exploration_methods) == "refine"
    ) %>%
    dplyr::distinct(dplyr::across(dplyr::all_of(group_cols)))

  if (!nrow(refine_groups)) {
    return(empty_result)
  }

  single_fit_sources <- confusion_individual %>%
    dplyr::distinct(
      dplyr::across(dplyr::all_of(group_cols)),
      .data$run_id,
      .data$variant_id
    ) %>%
    dplyr::group_by(dplyr::across(dplyr::all_of(group_cols))) %>%
    dplyr::summarise(
      source_run_id = dplyr::first(.data$run_id),
      source_variant_id = dplyr::first(.data$variant_id),
      n_individual_fits = dplyr::n(),
      .groups = "drop"
    ) %>%
    dplyr::filter(.data$n_individual_fits == 1L)

  eligible_groups <- refine_groups %>%
    dplyr::inner_join(single_fit_sources, by = group_cols)

  if (!nrow(eligible_groups)) {
    return(empty_result)
  }

  existing_agg <- confusion_agg %>%
    dplyr::distinct(dplyr::across(dplyr::all_of(c(group_cols, "agg_method"))))

  missing_agg_groups <- eligible_groups %>%
    tidyr::crossing(agg_method = agg_methods) %>%
    dplyr::anti_join(existing_agg, by = c(group_cols, "agg_method")) %>%
    dplyr::rename(fill_agg_method = .data$agg_method) %>%
    dplyr::left_join(task_lookup, by = group_cols)

  if (!nrow(missing_agg_groups)) {
    return(empty_result)
  }

  repaired_confusion <- confusion_individual %>%
    dplyr::inner_join(
      missing_agg_groups,
      by = setNames(
        c(group_cols, "source_run_id", "source_variant_id"),
        c(group_cols, "run_id", "variant_id")
      ),
      relationship = "many-to-many"
    ) %>%
    dplyr::mutate(
      explore_method = "aggregation",
      variant_id = .data$fill_agg_method,
      agg_method = .data$fill_agg_method
    ) %>%
    dplyr::select(dplyr::all_of(names(confusion_individual)))

  if ("variant_id" %in% names(confusion_agg)) {
    confusion_agg$variant_id <- as.character(confusion_agg$variant_id)
    repaired_confusion$variant_id <- as.character(repaired_confusion$variant_id)
  }

  repair_summary <- missing_agg_groups %>%
    dplyr::count(
      .data$task_id,
      .data$spec_name,
      .data$use_case_id,
      .data$fill_agg_method,
      name = "n_backfilled_groups",
      sort = TRUE
    ) %>%
    dplyr::rename(agg_method = .data$fill_agg_method)

  repaired_groups <- missing_agg_groups %>%
    dplyr::rename(agg_method = .data$fill_agg_method) %>%
    dplyr::arrange(.data$task_id, .data$spec_name, .data$dataset_bundle_id, .data$agg_method)

  list(
    confusion_agg = dplyr::bind_rows(confusion_agg, repaired_confusion),
    repaired_groups = repaired_groups,
    repair_summary = repair_summary
  )
}

#' @keywords internal
backfill_terminal_scaling_agg_confusion <- function(scaling_bins_raw,
                                                    confusion_agg,
                                                    run_info,
                                                    agg_methods) {
  agg_methods <- unique(normalize_agg_method_id(agg_methods))
  confusion_agg <- normalize_agg_method_df(confusion_agg)

  empty_result <- list(
    confusion_agg = confusion_agg,
    repaired_groups = tibble::tibble(),
    repair_summary = tibble::tibble()
  )

  if (is.null(scaling_bins_raw) || !nrow(scaling_bins_raw) ||
      is.null(run_info) || !nrow(run_info) || !length(agg_methods)) {
    return(empty_result)
  }

  scaling_bins_raw <- normalize_agg_method_df(scaling_bins_raw)
  scaling_methods <- intersect(unique(na.omit(scaling_bins_raw$agg_method)), agg_methods)
  if (!length(scaling_methods)) {
    return(empty_result)
  }

  mapping_base <- run_info %>%
    dplyr::filter(
      !is.na(.data$spec_name),
      !.data$spec_name %in% c("baseline-single", "truth-warm"),
      !is.na(.data$dataset_bundle_id),
      !is.na(.data$use_case_id),
      !is.na(.data$group_key)
    ) %>%
    dplyr::group_by(.data$dataset_bundle_id, .data$spec_name, .data$annotation_r2) %>%
    dplyr::summarise(
      n_group_rows = dplyr::n_distinct(paste(.data$use_case_id, .data$group_key, sep = "|")),
      .groups = "drop"
    )

  ambiguous_mappings <- mapping_base %>%
    dplyr::filter(.data$n_group_rows > 1L)
  if (nrow(ambiguous_mappings) > 0L) {
    stop(
      "Terminal scaling backfill found ambiguous run-info mappings for: ",
      paste(
        sprintf(
          "%s/%s/%s",
          ambiguous_mappings$dataset_bundle_id,
          ambiguous_mappings$spec_name,
          ifelse(is.na(ambiguous_mappings$annotation_r2), "overall", ambiguous_mappings$annotation_r2)
        ),
        collapse = ", "
      )
    )
  }

  mapping <- run_info %>%
    dplyr::filter(
      !is.na(.data$spec_name),
      !.data$spec_name %in% c("baseline-single", "truth-warm"),
      !is.na(.data$dataset_bundle_id),
      !is.na(.data$use_case_id),
      !is.na(.data$group_key)
    ) %>%
    dplyr::group_by(.data$dataset_bundle_id, .data$spec_name, .data$annotation_r2) %>%
    dplyr::summarise(
      run_id = min(.data$run_id, na.rm = TRUE),
      task_id = dplyr::first(.data$task_id),
      use_case_id = dplyr::first(.data$use_case_id),
      group_key = dplyr::first(.data$group_key),
      c_value = dplyr::first(.data$c_value),
      sigma_0_2_scalar = dplyr::first(.data$sigma_0_2_scalar),
      .groups = "drop"
    )

  terminal_pure <- scaling_bins_raw %>%
    dplyr::filter(is.na(.data$resolution), .data$agg_method %in% scaling_methods) %>%
    dplyr::group_by(.data$dataset_bundle_id, .data$spec_name, .data$annotation_r2, .data$agg_method) %>%
    dplyr::filter(.data$n_ensemble == max(.data$n_ensemble, na.rm = TRUE)) %>%
    dplyr::ungroup()

  terminal_inter <- scaling_bins_raw %>%
    dplyr::filter(!is.na(.data$resolution), .data$agg_method %in% scaling_methods) %>%
    dplyr::mutate(.resolution_area = parse_resolution_area(.data$resolution)) %>%
    dplyr::filter(is.finite(.data$.resolution_area)) %>%
    dplyr::group_by(.data$dataset_bundle_id, .data$spec_name, .data$annotation_r2, .data$agg_method) %>%
    dplyr::slice_max(order_by = .data$.resolution_area, n = 1L, with_ties = TRUE) %>%
    dplyr::ungroup() %>%
    dplyr::select(-.data$.resolution_area)

  terminal_bins <- dplyr::bind_rows(terminal_pure, terminal_inter)
  if (!nrow(terminal_bins)) {
    return(empty_result)
  }

  terminal_group_keys <- terminal_bins %>%
    dplyr::distinct(.data$dataset_bundle_id, .data$spec_name, .data$annotation_r2, .data$agg_method)

  existing_agg <- confusion_agg %>%
    dplyr::distinct(.data$dataset_bundle_id, .data$spec_name, .data$annotation_r2, .data$agg_method)

  missing_groups <- terminal_group_keys %>%
    dplyr::anti_join(
      existing_agg,
      by = c("dataset_bundle_id", "spec_name", "annotation_r2", "agg_method")
    ) %>%
    dplyr::left_join(
      mapping,
      by = c("dataset_bundle_id", "spec_name", "annotation_r2"),
      relationship = "many-to-one"
    )

  if (!nrow(missing_groups)) {
    return(empty_result)
  }

  repaired_confusion <- terminal_bins %>%
    dplyr::inner_join(
      missing_groups,
      by = c("dataset_bundle_id", "spec_name", "annotation_r2", "agg_method"),
      relationship = "many-to-many"
    ) %>%
    dplyr::mutate(
      explore_method = "aggregation",
      variant_id = .data$agg_method
    ) %>%
    dplyr::select(
      dplyr::all_of(c(
        "run_id", "explore_method", "variant_id", "agg_method",
        "pip_threshold", "n_causal_at_bucket", "n_noncausal_at_bucket",
        "dataset_bundle_id", "spec_name", "use_case_id", "annotation_r2",
        "group_key", "c_value", "sigma_0_2_scalar"
      ))
    )

  if ("variant_id" %in% names(confusion_agg)) {
    confusion_agg$variant_id <- as.character(confusion_agg$variant_id)
    repaired_confusion$variant_id <- as.character(repaired_confusion$variant_id)
  }

  repaired_groups <- missing_groups %>%
    dplyr::arrange(.data$task_id, .data$spec_name, .data$dataset_bundle_id, .data$agg_method)

  repair_summary <- repaired_groups %>%
    dplyr::count(
      .data$task_id,
      .data$spec_name,
      .data$use_case_id,
      .data$agg_method,
      name = "n_backfilled_groups",
      sort = TRUE
    )

  list(
    confusion_agg = dplyr::bind_rows(confusion_agg, repaired_confusion),
    repaired_groups = repaired_groups,
    repair_summary = repair_summary
  )
}

#' @keywords internal
parse_resolution_area <- function(x) {
  parts <- strsplit(as.character(x), "x", fixed = TRUE)
  vapply(parts, function(p) {
    nums <- suppressWarnings(as.numeric(p))
    if (length(nums) < 1L || any(!is.finite(nums))) return(NA_real_)
    prod(nums)
  }, numeric(1L))
}

#' @keywords internal
overwrite_scaling_terminal_metrics <- function(summary_df, full_targets_df) {
  if (is.null(summary_df) || !nrow(summary_df)) return(summary_df)

  summary_df <- normalize_agg_method_df(summary_df)
  full_targets_df <- normalize_agg_method_df(full_targets_df) %>%
    dplyr::rename(AUPRC_full = .data$AUPRC, tpr05_full = .data$tpr_fpr05)

  pure_terminal_flags <- summary_df %>%
    dplyr::filter(is.na(.data$resolution)) %>%
    dplyr::group_by(.data$spec_name, .data$agg_method, .data$annotation_r2) %>%
    dplyr::mutate(
      .terminal = .data$n_ensemble == max(.data$n_ensemble, na.rm = TRUE)
    ) %>%
    dplyr::ungroup()

  inter_terminal_flags <- summary_df %>%
    dplyr::filter(!is.na(.data$resolution)) %>%
    dplyr::mutate(.resolution_area = parse_resolution_area(.data$resolution)) %>%
    dplyr::filter(is.finite(.data$.resolution_area)) %>%
    dplyr::group_by(.data$spec_name, .data$agg_method, .data$annotation_r2) %>%
    dplyr::mutate(
      .terminal = dplyr::min_rank(dplyr::desc(.data$.resolution_area)) == 1L
    ) %>%
    dplyr::ungroup() %>%
    dplyr::select(-.data$.resolution_area)

  terminal_flags <- dplyr::bind_rows(pure_terminal_flags, inter_terminal_flags)

  missing_targets <- terminal_flags %>%
    dplyr::filter(.data$.terminal) %>%
    dplyr::anti_join(
      full_targets_df %>%
        dplyr::select(.data$spec_name, .data$agg_method, .data$annotation_r2),
      by = c("spec_name", "agg_method", "annotation_r2")
    )

  if (nrow(missing_targets) > 0L) {
    stop(
      "Missing pooled aggregated target rows for terminal scaling points: ",
      paste(
        sprintf(
          "%s/%s/%s",
          missing_targets$spec_name,
          missing_targets$agg_method,
          ifelse(is.na(missing_targets$annotation_r2), "overall", missing_targets$annotation_r2)
        ),
        collapse = ", "
      )
    )
  }

  terminal_flags %>%
    dplyr::left_join(
      full_targets_df,
      by = c("spec_name", "agg_method", "annotation_r2")
    ) %>%
    dplyr::mutate(
      AUPRC_mean = dplyr::if_else(.data$.terminal, .data$AUPRC_full, .data$AUPRC_mean),
      tpr05_mean = dplyr::if_else(.data$.terminal, .data$tpr05_full, .data$tpr05_mean),
      AUPRC_se = dplyr::if_else(.data$.terminal, NA_real_, .data$AUPRC_se),
      tpr05_se = dplyr::if_else(.data$.terminal, NA_real_, .data$tpr05_se),
      n_reps = dplyr::if_else(.data$.terminal, 1L, .data$n_reps)
    ) %>%
    dplyr::select(
      -.data$AUPRC_full,
      -.data$tpr05_full,
      -.data$.terminal
    )
}
