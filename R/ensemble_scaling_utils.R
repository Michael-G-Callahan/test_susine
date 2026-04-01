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
    dplyr::group_by(.data$spec_name, .data$agg_method, .data$annotation_r2) %>%
    dplyr::mutate(
      .terminal = .data$.resolution_area == max(.data$.resolution_area, na.rm = TRUE)
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
