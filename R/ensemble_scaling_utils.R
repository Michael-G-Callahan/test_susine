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
annotation_r2_join_key <- function(x) {
  x_chr <- as.character(x)
  x_chr[is.na(x)] <- "__NA__"
  x_chr
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
backfill_pure_refine_scaling_bins <- function(scaling_bins_raw,
                                              confusion_agg,
                                              model_metrics,
                                              run_info,
                                              scaling_n_ens_sizes = c(4L, 8L, 16L, 32L, 64L)) {
  if (is.null(scaling_bins_raw) ||
      is.null(confusion_agg) || !nrow(confusion_agg) ||
      is.null(model_metrics) || !nrow(model_metrics) ||
      is.null(run_info) || !nrow(run_info)) {
    return(list(
      scaling_bins_raw = scaling_bins_raw,
      repaired_groups = tibble::tibble(),
      repair_summary = tibble::tibble()
    ))
  }

  scaling_bins_raw <- normalize_agg_method_df(scaling_bins_raw)
  confusion_agg <- normalize_agg_method_df(confusion_agg)
  model_metrics <- normalize_agg_method_df(model_metrics)

  pure_refine_plans <- run_info %>%
    dplyr::filter(.data$exploration_methods == "refine") %>%
    dplyr::group_by(
      .data$dataset_bundle_id,
      .data$spec_name,
      .data$use_case_id,
      .data$annotation_r2,
      .data$group_key
    ) %>%
    dplyr::summarise(
      planned_n = max(.data$refine_step, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    dplyr::filter(is.finite(.data$planned_n), .data$planned_n >= 1L) %>%
    dplyr::mutate(annotation_r2_key = annotation_r2_join_key(.data$annotation_r2))

  if (!nrow(pure_refine_plans)) {
    return(list(
      scaling_bins_raw = scaling_bins_raw,
      repaired_groups = tibble::tibble(),
      repair_summary = tibble::tibble()
    ))
  }

  actual_refine_counts <- model_metrics %>%
    dplyr::filter(
      .data$explore_method == "refine",
      is.na(.data$agg_method)
    ) %>%
    dplyr::group_by(
      .data$dataset_bundle_id,
      .data$spec_name,
      .data$use_case_id,
      .data$annotation_r2,
      .data$group_key
    ) %>%
    dplyr::summarise(
      actual_n = dplyr::n_distinct(.data$run_id),
      .groups = "drop"
    ) %>%
    dplyr::mutate(annotation_r2_key = annotation_r2_join_key(.data$annotation_r2))

  planned_thresholds <- pure_refine_plans %>%
    dplyr::left_join(
      actual_refine_counts,
      by = c("dataset_bundle_id", "spec_name", "use_case_id", "annotation_r2_key", "group_key")
    ) %>%
    dplyr::select(-tidyselect::any_of("annotation_r2.y")) %>%
    dplyr::rename(annotation_r2 = .data$annotation_r2.x) %>%
    dplyr::mutate(actual_n = dplyr::coalesce(.data$actual_n, 0L)) %>%
    dplyr::rowwise() %>%
    dplyr::mutate(
      requested_thresholds = list(
        as.integer(scaling_n_ens_sizes[scaling_n_ens_sizes <= .data$planned_n])
      )
    ) %>%
    tidyr::unnest_longer(.data$requested_thresholds, values_to = "n_ensemble") %>%
    dplyr::ungroup() %>%
    dplyr::filter(.data$n_ensemble > .data$actual_n)

  if (!nrow(planned_thresholds)) {
    return(list(
      scaling_bins_raw = scaling_bins_raw,
      repaired_groups = tibble::tibble(),
      repair_summary = tibble::tibble()
    ))
  }

  full_refine_bins <- confusion_agg %>%
    dplyr::filter(.data$explore_method == "aggregation") %>%
    dplyr::group_by(
      .data$dataset_bundle_id,
      .data$spec_name,
      .data$use_case_id,
      .data$annotation_r2,
      .data$group_key,
      .data$agg_method,
      .data$pip_threshold
    ) %>%
    dplyr::summarise(
      n_causal_at_bucket = sum(.data$n_causal_at_bucket, na.rm = TRUE),
      n_noncausal_at_bucket = sum(.data$n_noncausal_at_bucket, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    dplyr::mutate(annotation_r2_key = annotation_r2_join_key(.data$annotation_r2))

  backfilled_rows <- planned_thresholds %>%
    dplyr::inner_join(
      full_refine_bins,
      by = c("dataset_bundle_id", "spec_name", "use_case_id", "annotation_r2_key", "group_key")
    ) %>%
    dplyr::select(-tidyselect::any_of("annotation_r2.y")) %>%
    dplyr::rename(annotation_r2 = .data$annotation_r2.x) %>%
    dplyr::mutate(
      resolution = NA_character_,
      rep = 1L
    ) %>%
    dplyr::select(
      .data$spec_name,
      .data$annotation_r2,
      .data$dataset_bundle_id,
      .data$n_ensemble,
      .data$resolution,
      .data$agg_method,
      .data$rep,
      .data$pip_threshold,
      .data$n_causal_at_bucket,
      .data$n_noncausal_at_bucket
    )

  if (!nrow(backfilled_rows)) {
    return(list(
      scaling_bins_raw = scaling_bins_raw,
      repaired_groups = tibble::tibble(),
      repair_summary = tibble::tibble()
    ))
  }

  existing_keys <- scaling_bins_raw %>%
    dplyr::mutate(annotation_r2_key = annotation_r2_join_key(.data$annotation_r2)) %>%
    dplyr::select(
      .data$spec_name,
      .data$annotation_r2_key,
      .data$dataset_bundle_id,
      .data$n_ensemble,
      .data$resolution,
      .data$agg_method,
      .data$rep,
      .data$pip_threshold
    ) %>%
    dplyr::distinct()

  backfilled_rows <- backfilled_rows %>%
    dplyr::mutate(annotation_r2_key = annotation_r2_join_key(.data$annotation_r2)) %>%
    dplyr::anti_join(
      existing_keys,
      by = c("spec_name", "annotation_r2_key", "dataset_bundle_id",
             "n_ensemble", "resolution", "agg_method", "rep", "pip_threshold")
    ) %>%
    dplyr::select(-.data$annotation_r2_key)

  repaired_groups <- planned_thresholds %>%
    dplyr::semi_join(
      backfilled_rows %>%
        dplyr::select(
          .data$spec_name,
          .data$annotation_r2,
          .data$dataset_bundle_id,
          .data$n_ensemble
        ) %>%
        dplyr::distinct(),
      by = c("spec_name", "annotation_r2", "dataset_bundle_id", "n_ensemble")
    ) %>%
    dplyr::arrange(.data$spec_name, .data$annotation_r2, .data$dataset_bundle_id, .data$n_ensemble)

  repair_summary <- repaired_groups %>%
    dplyr::group_by(.data$spec_name, .data$annotation_r2) %>%
    dplyr::summarise(
      n_backfilled_groups = dplyr::n(),
      .groups = "drop"
    )

  list(
    scaling_bins_raw = dplyr::bind_rows(scaling_bins_raw, backfilled_rows),
    repaired_groups = repaired_groups,
    repair_summary = repair_summary
  )
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
filter_shared_plot_metric_rows <- function(df, r2_filter = NULL) {
  if (is.null(df) || !nrow(df)) return(tibble::tibble())

  agg_rows <- df %>%
    dplyr::filter(
      .data$series_type == "aggregated",
      !.data$spec_name %in% c("baseline-single", "truth-warm")
    )

  if (is.null(r2_filter)) {
    return(agg_rows %>% dplyr::filter(is.na(.data$annotation_r2)))
  }

  dependent_keys <- agg_rows %>%
    dplyr::filter(!is.na(.data$annotation_r2)) %>%
    dplyr::distinct(.data$spec_name, .data$use_case_id) %>%
    dplyr::mutate(.annotation_dependent = TRUE)

  agg_rows %>%
    dplyr::left_join(dependent_keys, by = c("spec_name", "use_case_id")) %>%
    dplyr::mutate(
      .annotation_dependent = dplyr::coalesce(.data$.annotation_dependent, FALSE)
    ) %>%
    dplyr::filter(
      (!.data$.annotation_dependent & is.na(.data$annotation_r2)) |
        (.data$.annotation_dependent & .data$annotation_r2 == r2_filter)
    ) %>%
    dplyr::select(-.data$.annotation_dependent)
}

#' @keywords internal
filter_truth_warm_reference_rows <- function(df, r2_filter = NULL) {
  if (is.null(df) || !nrow(df)) return(tibble::tibble())

  truth_rows <- if ("series_type" %in% names(df)) {
    df %>%
      dplyr::filter(.data$series_type == "truth_warm")
  } else {
    tibble::as_tibble(df)
  }

  if (is.null(r2_filter)) {
    return(truth_rows %>% dplyr::filter(is.na(.data$annotation_r2)))
  }

  dependent_use_cases <- truth_rows %>%
    dplyr::filter(!is.na(.data$annotation_r2)) %>%
    dplyr::distinct(.data$use_case_id) %>%
    dplyr::mutate(.annotation_dependent = TRUE)

  truth_rows %>%
    dplyr::left_join(dependent_use_cases, by = "use_case_id") %>%
    dplyr::mutate(
      .annotation_dependent = dplyr::coalesce(.data$.annotation_dependent, FALSE)
    ) %>%
    dplyr::filter(
      (!.data$.annotation_dependent & is.na(.data$annotation_r2)) |
        (.data$.annotation_dependent & .data$annotation_r2 == r2_filter)
    ) %>%
    dplyr::select(-.data$.annotation_dependent)
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
backfill_missing_agg_confusion_from_snps <- function(snp_files,
                                                     confusion_individual,
                                                     confusion_agg,
                                                     run_info,
                                                     model_metrics,
                                                     agg_methods,
                                                     pip_breaks) {
  agg_methods <- unique(normalize_agg_method_id(agg_methods))
  confusion_agg <- normalize_agg_method_df(confusion_agg)

  empty_result <- list(
    confusion_agg = confusion_agg,
    repaired_groups = tibble::tibble(),
    repair_summary = tibble::tibble()
  )

  if (is.null(snp_files) || !length(snp_files) ||
      is.null(confusion_individual) || !nrow(confusion_individual) ||
      is.null(run_info) || !nrow(run_info) ||
      is.null(model_metrics) || !nrow(model_metrics) ||
      !length(agg_methods) || is.null(pip_breaks) || !length(pip_breaks)) {
    return(empty_result)
  }

  group_cols <- c("spec_name", "use_case_id", "dataset_bundle_id", "annotation_r2", "group_key")

  existing_agg <- confusion_agg %>%
    dplyr::distinct(dplyr::across(dplyr::all_of(c(group_cols, "agg_method"))))

  individual_runs <- confusion_individual %>%
    dplyr::distinct(dplyr::across(dplyr::all_of(group_cols)), .data$run_id) %>%
    dplyr::group_by(dplyr::across(dplyr::all_of(group_cols))) %>%
    dplyr::summarise(
      run_ids = list(sort(unique(.data$run_id))),
      n_individual_fits = dplyr::n_distinct(.data$run_id),
      .groups = "drop"
    ) %>%
    dplyr::filter(.data$n_individual_fits > 1L)

  missing_groups <- individual_runs %>%
    tidyr::crossing(agg_method = agg_methods) %>%
    dplyr::anti_join(existing_agg, by = c(group_cols, "agg_method"))

  if (!nrow(missing_groups)) {
    return(empty_result)
  }

  group_meta <- run_info %>%
    dplyr::filter(
      !is.na(.data$spec_name),
      !.data$spec_name %in% c("baseline-single", "truth-warm"),
      !is.na(.data$use_case_id),
      !is.na(.data$dataset_bundle_id),
      !is.na(.data$group_key)
    ) %>%
    dplyr::group_by(dplyr::across(dplyr::all_of(group_cols))) %>%
    dplyr::summarise(
      run_id = min(.data$run_id, na.rm = TRUE),
      task_id = dplyr::first(.data$task_id),
      c_value = dplyr::first(.data$c_value),
      sigma_0_2_scalar = dplyr::first(.data$sigma_0_2_scalar),
      .groups = "drop"
    )

  missing_groups <- missing_groups %>%
    dplyr::left_join(group_meta, by = group_cols)

  target_run_rows <- missing_groups %>%
    dplyr::select(dplyr::all_of(group_cols), .data$agg_method, .data$run_ids) %>%
    tidyr::unnest_longer(.data$run_ids, values_to = "run_id") %>%
    dplyr::mutate(run_id = as.integer(.data$run_id))

  target_run_ids <- sort(unique(target_run_rows$run_id))
  if (!length(target_run_ids)) {
    return(empty_result)
  }

  snps_needed <- arrow::open_dataset(snp_files, format = "parquet") %>%
    dplyr::filter(.data$run_id %in% target_run_ids) %>%
    dplyr::select(.data$run_id, .data$snp_index, .data$pip, .data$causal) %>%
    dplyr::collect()

  if (!nrow(snps_needed)) {
    return(empty_result)
  }

  elbo_info <- model_metrics %>%
    dplyr::filter(!is.na(.data$run_id), is.na(.data$agg_method)) %>%
    dplyr::distinct(.data$run_id, .data$elbo_final)

  target_snps <- target_run_rows %>%
    dplyr::left_join(snps_needed, by = "run_id") %>%
    dplyr::left_join(elbo_info, by = "run_id")

  if (!nrow(target_snps)) {
    return(empty_result)
  }

  repaired_list <- vector("list", nrow(missing_groups))
  repaired_group_rows <- vector("list", nrow(missing_groups))

  for (i in seq_len(nrow(missing_groups))) {
    group_row <- missing_groups[i, , drop = FALSE]
    ds_snps <- target_snps %>%
      dplyr::filter(
        .data$spec_name == group_row$spec_name[[1]],
        .data$use_case_id == group_row$use_case_id[[1]],
        .data$dataset_bundle_id == group_row$dataset_bundle_id[[1]],
        ((is.na(.data$annotation_r2) & is.na(group_row$annotation_r2[[1]])) |
           .data$annotation_r2 == group_row$annotation_r2[[1]]),
        .data$group_key == group_row$group_key[[1]],
        .data$agg_method == group_row$agg_method[[1]]
      )

    if (!nrow(ds_snps)) next

    sel_ids <- sort(unique(ds_snps$run_id))
    pm <- .build_pip_mat(ds_snps, sel_ids, ds_snps %>% dplyr::distinct(.data$run_id, .data$elbo_final))
    agg_pip <- aggregate_pip_matrix(pm$pip_mat, pm$elbos, group_row$agg_method[[1]])
    bins <- compute_bins_from_pip_vec(agg_pip, pm$causal_vec, pip_breaks)

    repaired_list[[i]] <- bins %>%
      dplyr::mutate(
        run_id = group_row$run_id[[1]],
        explore_method = "aggregation",
        variant_id = group_row$agg_method[[1]],
        agg_method = group_row$agg_method[[1]],
        dataset_bundle_id = group_row$dataset_bundle_id[[1]],
        spec_name = group_row$spec_name[[1]],
        use_case_id = group_row$use_case_id[[1]],
        annotation_r2 = group_row$annotation_r2[[1]],
        group_key = group_row$group_key[[1]],
        c_value = group_row$c_value[[1]],
        sigma_0_2_scalar = group_row$sigma_0_2_scalar[[1]]
      ) %>%
      dplyr::select(dplyr::all_of(c(
        "run_id", "explore_method", "variant_id", "agg_method",
        "pip_threshold", "n_causal_at_bucket", "n_noncausal_at_bucket",
        "dataset_bundle_id", "spec_name", "use_case_id", "annotation_r2",
        "group_key", "c_value", "sigma_0_2_scalar"
      )))

    repaired_group_rows[[i]] <- group_row %>%
      dplyr::select(.data$task_id, dplyr::all_of(group_cols), .data$agg_method, .data$n_individual_fits)
  }

  repaired_confusion <- dplyr::bind_rows(repaired_list)
  repaired_groups <- dplyr::bind_rows(repaired_group_rows) %>%
    dplyr::arrange(.data$task_id, .data$spec_name, .data$dataset_bundle_id, .data$agg_method)

  if (!nrow(repaired_confusion)) {
    return(empty_result)
  }

  if ("variant_id" %in% names(confusion_agg)) {
    confusion_agg$variant_id <- as.character(confusion_agg$variant_id)
    repaired_confusion$variant_id <- as.character(repaired_confusion$variant_id)
  }

  repair_summary <- repaired_groups %>%
    dplyr::count(
      .data$task_id, .data$spec_name, .data$use_case_id, .data$agg_method,
      name = "n_backfilled_groups", sort = TRUE
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
spec_single_fit_use_case <- function(spec_name) {
  spec_name <- as.character(spec_name)
  dplyr::case_when(
    startsWith(spec_name, "A-") ~ "susine_vanilla",
    startsWith(spec_name, "B-") ~ "susine_eb_clamped_scale_var_nonneg",
    startsWith(spec_name, "C-") ~ "susine_functional_mu",
    TRUE ~ NA_character_
  )
}

#' @keywords internal
append_scaling_single_fit_points <- function(summary_df,
                                             auprc_plot_table,
                                             tpr05_plot_table) {
  if (is.null(summary_df) || !nrow(summary_df)) return(summary_df)
  if (is.null(auprc_plot_table) || !nrow(auprc_plot_table) ||
      is.null(tpr05_plot_table) || !nrow(tpr05_plot_table)) {
    return(summary_df)
  }

  pure_keys <- summary_df %>%
    dplyr::filter(is.na(.data$resolution)) %>%
    dplyr::distinct(.data$spec_name, .data$annotation_r2, .data$agg_method) %>%
    dplyr::mutate(
      baseline_use_case_id = spec_single_fit_use_case(.data$spec_name)
    ) %>%
    dplyr::filter(!is.na(.data$baseline_use_case_id))

  if (!nrow(pure_keys)) {
    return(summary_df)
  }

  baseline_auprc <- auprc_plot_table %>%
    dplyr::filter(
      .data$series_type == "baseline_single",
      .data$spec_name == "baseline-single"
    ) %>%
    dplyr::select(
      use_case_id = .data$use_case_id,
      annotation_r2 = .data$annotation_r2,
      AUPRC = .data$AUPRC
    )

  baseline_tpr05 <- tpr05_plot_table %>%
    dplyr::filter(
      .data$series_type == "baseline_single",
      .data$spec_name == "baseline-single"
    ) %>%
    dplyr::select(
      use_case_id = .data$use_case_id,
      annotation_r2 = .data$annotation_r2,
      tpr_fpr05 = .data$tpr_fpr05
    )

  dependent_use_cases <- dplyr::bind_rows(
    baseline_auprc %>% dplyr::select(.data$use_case_id, .data$annotation_r2),
    baseline_tpr05 %>% dplyr::select(.data$use_case_id, .data$annotation_r2)
  ) %>%
    dplyr::filter(!is.na(.data$annotation_r2)) %>%
    dplyr::distinct(.data$use_case_id) %>%
    dplyr::mutate(.annotation_dependent = TRUE)

  single_fit_points <- pure_keys %>%
    dplyr::left_join(
      dependent_use_cases,
      by = c("baseline_use_case_id" = "use_case_id")
    ) %>%
    dplyr::mutate(
      .annotation_dependent = dplyr::coalesce(.data$.annotation_dependent, FALSE),
      baseline_annotation_r2 = dplyr::if_else(
        .data$.annotation_dependent,
        .data$annotation_r2,
        NA_real_
      )
    ) %>%
    dplyr::left_join(
      baseline_auprc,
      by = c(
        "baseline_use_case_id" = "use_case_id",
        "baseline_annotation_r2" = "annotation_r2"
      )
    ) %>%
    dplyr::left_join(
      baseline_tpr05,
      by = c(
        "baseline_use_case_id" = "use_case_id",
        "baseline_annotation_r2" = "annotation_r2"
      )
    ) %>%
    dplyr::transmute(
      spec_name = .data$spec_name,
      annotation_r2 = .data$annotation_r2,
      n_ensemble = 1L,
      resolution = NA_character_,
      agg_method = .data$agg_method,
      AUPRC_mean = .data$AUPRC,
      AUPRC_se = NA_real_,
      tpr05_mean = .data$tpr_fpr05,
      tpr05_se = NA_real_,
      n_reps = 1L
    ) %>%
    dplyr::filter(
      is.finite(.data$AUPRC_mean) | is.finite(.data$tpr05_mean)
    )

  if (!nrow(single_fit_points)) {
    return(summary_df)
  }

  summary_wo_n1 <- summary_df %>%
    dplyr::filter(!(is.na(.data$resolution) & .data$n_ensemble == 1L))

  dplyr::bind_rows(summary_wo_n1, single_fit_points) %>%
    dplyr::arrange(
      .data$spec_name,
      .data$annotation_r2,
      .data$agg_method,
      .data$n_ensemble,
      .data$resolution
    )
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
