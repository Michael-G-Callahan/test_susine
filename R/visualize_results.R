# visualize_results.R
# Helper functions for generating plots from simulation results
# ------------------------------------------------------------------------------

#' Get use case labels from catalog
#'
#' @param use_case_ids Character vector of use_case_id values to look up
#' @return Named character vector mapping use_case_id to label
#' @export
get_use_case_labels <- function(use_case_ids = NULL) {

  catalog <- use_case_catalog()
  if (!is.null(use_case_ids)) {
    catalog <- dplyr::filter(catalog, use_case_id %in% use_case_ids)
  }
  labels <- stats::setNames(catalog$label, catalog$use_case_id)
  labels
}

#' Compute Power vs FDR curves from pre-aggregated confusion matrix
#'
#' @param confusion_df Data frame with columns: grouping vars, pip_threshold,
#'   n_causal_at_bucket, n_noncausal_at_bucket
#' @param group_vars Character vector of columns to group by (e.g., "use_case_id")
#' @return Data frame with columns: group_vars, pip_threshold, TP, FP, power, fdr, precision
#' @export
compute_power_fdr_curves <- function(confusion_df, group_vars = "use_case_id") {

  # First aggregate per-bucket counts across any extra dimensions
confusion_df %>%
    dplyr::group_by(dplyr::across(dplyr::all_of(group_vars)), pip_threshold) %>%
    dplyr::summarise(
      n_causal_at_bucket = sum(n_causal_at_bucket, na.rm = TRUE),
      n_noncausal_at_bucket = sum(n_noncausal_at_bucket, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    # Compute cumulative TP/FP per group (descending pip)
    dplyr::group_by(dplyr::across(dplyr::all_of(group_vars))) %>%
    dplyr::arrange(dplyr::desc(pip_threshold), .by_group = TRUE) %>%
    dplyr::mutate(
      TP = cumsum(n_causal_at_bucket),
      FP = cumsum(n_noncausal_at_bucket),
      total_causal = sum(n_causal_at_bucket),
      total_noncausal = sum(n_noncausal_at_bucket),
      FN = total_causal - TP,
      TN = total_noncausal - FP,
      precision = dplyr::if_else(TP + FP > 0, TP / (TP + FP), NA_real_),
      recall = dplyr::if_else(total_causal > 0, TP / total_causal, NA_real_),
      power = recall,
      fdr = 1 - precision
    ) %>%
    dplyr::ungroup() %>%
    dplyr::filter(!is.na(precision), !is.na(recall))
}

#' Compute AUPRC from pre-aggregated confusion matrix
#'
#' Integrates the precision-recall curve using trapezoidal rule,
#' anchored at (recall=0, precision=1).
#'
#' @param confusion_df Data frame with per-bucket counts
#' @param group_vars Character vector of columns to group by for AUPRC calculation
#' @return Data frame with group_vars and AUPRC column
#' @export
compute_auprc_from_confusion <- function(confusion_df, group_vars) {

  # First get power/fdr curves
  curves <- compute_power_fdr_curves(confusion_df, group_vars)

  # Compute AUPRC per group
  curves %>%
    dplyr::group_by(dplyr::across(dplyr::all_of(group_vars))) %>%
    dplyr::summarise(
      AUPRC = compute_auprc_single(precision, recall),
      .groups = "drop"
    )
}

#' Compute AUPRC for a single precision-recall curve
#'
#' @param precision Numeric vector of precision values
#' @param recall Numeric vector of recall values
#' @return Single AUPRC value
#' @keywords internal
compute_auprc_single <- function(precision, recall) {
  df <- data.frame(precision = precision, recall = recall) %>%
    dplyr::filter(!is.na(precision), !is.na(recall)) %>%
    dplyr::arrange(recall)

  if (nrow(df) < 2) return(NA_real_)

 # Anchor at recall=0 with precision=1 (standard AUPRC convention)
  if (min(df$recall) > 0) {
    df <- dplyr::bind_rows(
      data.frame(recall = 0, precision = 1),
      df
    )
  }

  # Trapezoidal integration
  auprc <- sum(diff(df$recall) * (utils::head(df$precision, -1) + utils::tail(df$precision, -1)) / 2)
  auprc
}

#' Create Power vs FDR plot
#'
#' @param curve_data Output from compute_power_fdr_curves()
#' @param color_var Column name to use for color aesthetic
#' @param color_labels Named vector mapping color_var values to labels
#' @param title Plot title
#' @param subtitle Plot subtitle
#' @return ggplot object
#' @export
plot_power_fdr <- function(curve_data,
                           color_var = "use_case_id",
                           color_labels = NULL,
                           title = "Power vs FDR",
                           subtitle = NULL) {

  # Apply labels if provided
  if (!is.null(color_labels)) {
    curve_data <- curve_data %>%
      dplyr::mutate(
        color_label = dplyr::coalesce(color_labels[.data[[color_var]]], .data[[color_var]])
      )
    color_aes <- "color_label"
  } else {
    color_aes <- color_var
  }

  ggplot2::ggplot(curve_data, ggplot2::aes(x = fdr, y = power, color = .data[[color_aes]])) +
    ggplot2::geom_line(linewidth = 0.8) +
    ggplot2::geom_point(size = 0.5, alpha = 0.5) +
    ggplot2::labs(
      title = title,
      subtitle = subtitle,
      x = "False Discovery Rate (1 - Precision)",
      y = "Power (Recall)",
      color = "Model"
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(legend.position = "bottom") +
    ggplot2::guides(color = ggplot2::guide_legend(nrow = 2))
}

#' Create AUPRC plot with y_noise on x-axis and facets
#'
#' @param auprc_data Data frame with AUPRC and grouping columns
#' @param color_var Column for color aesthetic
#' @param color_labels Named vector for color labels
#' @param linetype_var Column for secondary grouping aesthetic (uses shape + linetype)
#' @param linetype_label Legend title for the secondary grouping variable
#' @param facet_row Column for facet rows (e.g., inflate_match)
#' @param facet_col Column for facet columns (e.g., annotation_r2)
#' @param title Plot title
#' @param subtitle Plot subtitle
#' @return ggplot object
#' @export
plot_auprc_informed <- function(auprc_data,
                                color_var = "use_case_id",
                                color_labels = NULL,
                                linetype_var = "gamma_shrink",
                                linetype_label = NULL,
                                facet_row = "inflate_match",
                                facet_col = "annotation_r2",
                                title = "AUPRC by Noise Level",
                                subtitle = NULL) {

  # Auto-generate legend title from variable name if not provided

if (is.null(linetype_label) && !is.null(linetype_var)) {
    # Convert snake_case to Title Case
    linetype_label <- gsub("_", " ", linetype_var)
    linetype_label <- tools::toTitleCase(linetype_label)
  }

  # Apply labels if provided
  if (!is.null(color_labels)) {
    auprc_data <- auprc_data %>%
      dplyr::mutate(
        color_label = dplyr::coalesce(color_labels[.data[[color_var]]], .data[[color_var]])
      )
    color_aes <- "color_label"
  } else {
    color_aes <- color_var
  }

  # Format secondary grouping variable (for shape/linetype)
  if (!is.null(linetype_var) && linetype_var %in% names(auprc_data)) {
    lt_values <- auprc_data[[linetype_var]]
    # Check if values are numeric or can be coerced to numeric
    if (is.numeric(lt_values)) {
      auprc_data <- auprc_data %>%
        dplyr::mutate(
          group_label = dplyr::if_else(
            is.na(.data[[linetype_var]]),
            "default",
            sprintf("%.2f", .data[[linetype_var]])
          )
        )
    } else {
      # Character values - use as-is
      auprc_data <- auprc_data %>%
        dplyr::mutate(
          group_label = dplyr::if_else(
            is.na(.data[[linetype_var]]) | .data[[linetype_var]] == "",
            "default",
            as.character(.data[[linetype_var]])
          )
        )
    }
    group_aes <- "group_label"
  } else {
    group_aes <- NULL
  }

  p <- ggplot2::ggplot(auprc_data, ggplot2::aes(x = y_noise, y = AUPRC, color = .data[[color_aes]]))

  if (!is.null(group_aes)) {
    # Use both shape and linetype for better visual distinction
    p <- p + ggplot2::aes(shape = .data[[group_aes]], linetype = .data[[group_aes]])
  }

  p <- p +
    ggplot2::geom_line(linewidth = 0.8) +
    ggplot2::geom_point(size = 2.5) +
    ggplot2::labs(
      title = title,
      subtitle = subtitle,
      x = "y_noise (Noise Level)",
      y = "AUPRC",
      color = "Model",
      shape = linetype_label,
      linetype = linetype_label
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(legend.position = "right") +
    # Use manual scales for more distinct shapes
    ggplot2::scale_shape_manual(values = c(16, 17, 15, 3, 4, 8, 1, 2, 0, 5, 6, 7))

  # Add facets if columns exist and have non-NA values
  if (!is.null(facet_row) && !is.null(facet_col) &&
      facet_row %in% names(auprc_data) && facet_col %in% names(auprc_data)) {
    p <- p + ggplot2::facet_grid(
      stats::reformulate(facet_col, facet_row),
      labeller = ggplot2::labeller(
        .rows = function(x) paste0(facet_row, ": ", x),
        .cols = function(x) paste0(facet_col, ": ", x)
      )
    )
  }

  p
}

#' Create simple AUPRC plot with one grouping var on x-axis
#'
#' For use cases with NA for the x-axis variable, draws horizontal lines.
#'
#' @param auprc_data Data frame with AUPRC and grouping columns
#' @param x_var Column for x-axis
#' @param color_var Column for color aesthetic
#' @param color_labels Named vector for color labels
#' @param title Plot title
#' @return ggplot object
#' @export
plot_auprc_simple <- function(auprc_data,
                              x_var,
                              color_var = "use_case_id",
                              color_labels = NULL,
                              title = NULL) {

  # Apply labels
  if (!is.null(color_labels)) {
    auprc_data <- auprc_data %>%
      dplyr::mutate(
        color_label = dplyr::coalesce(color_labels[.data[[color_var]]], .data[[color_var]])
      )
    color_aes <- "color_label"
  } else {
    color_aes <- color_var
  }

  # Split data into those with valid x values and those with NA
  data_valid <- auprc_data %>% dplyr::filter(!is.na(.data[[x_var]]))
  data_na <- auprc_data %>% dplyr::filter(is.na(.data[[x_var]]))

  # Get x-axis range for horizontal lines
  if (nrow(data_valid) > 0) {
    x_range <- range(data_valid[[x_var]], na.rm = TRUE)
  } else {
    x_range <- c(0, 1)
  }

  p <- ggplot2::ggplot() +
    ggplot2::theme_minimal() +
    ggplot2::labs(
      title = title %||% paste("AUPRC by", x_var),
      x = x_var,
      y = "AUPRC",
      color = "Model"
    ) +
    ggplot2::theme(legend.position = "bottom")

  # Add lines for valid data
  if (nrow(data_valid) > 0) {
    p <- p +
      ggplot2::geom_line(
        data = data_valid,
        ggplot2::aes(x = .data[[x_var]], y = AUPRC, color = .data[[color_aes]]),
        linewidth = 0.8
      ) +
      ggplot2::geom_point(
        data = data_valid,
        ggplot2::aes(x = .data[[x_var]], y = AUPRC, color = .data[[color_aes]]),
        size = 2
      )
  }

  # Add horizontal lines for NA x values (span full range)
  if (nrow(data_na) > 0) {
    # Get mean AUPRC per use_case for horizontal lines
    na_summary <- data_na %>%
      dplyr::group_by(.data[[color_aes]]) %>%
      dplyr::summarise(AUPRC = mean(AUPRC, na.rm = TRUE), .groups = "drop")

    for (i in seq_len(nrow(na_summary))) {
      p <- p +
        ggplot2::geom_hline(
          yintercept = na_summary$AUPRC[i],
          color = scales::hue_pal()(nrow(na_summary))[i],
          linetype = "dashed",
          linewidth = 0.8
        )
    }
  }

  p
}

#' Create CS metrics boxplots
#'
#' Creates a multi-panel figure with boxplots of CS metrics by p_star and use_case.
#'
#' @param effect_metrics Data frame with effect-level metrics
#' @param model_metrics Data frame with model-level metrics
#' @param use_case_filter Character vector of use_case_ids to include
#' @param color_labels Named vector mapping use_case_id to labels
#' @return ggplot object (cowplot combined)
#' @export
plot_cs_metrics_boxplots <- function(effect_metrics,
                                     model_metrics,
                                     use_case_filter,
                                     color_labels = NULL) {

  # Filter to purity_filtered and selected use cases
  effects <- effect_metrics %>%
    dplyr::filter(filtering == "purity_filtered",
                  use_case_id %in% use_case_filter)

  models <- model_metrics %>%
    dplyr::filter(filtering == "purity_filtered",
                  use_case_id %in% use_case_filter)

  # Need p_star - join from a source that has it
  # For now, assume it's in model_metrics or we need to get it from run_table
  # If not present, we'll need to handle that

  # Apply labels
  if (!is.null(color_labels)) {
    effects <- effects %>%
      dplyr::mutate(label = dplyr::coalesce(color_labels[use_case_id], use_case_id))
    models <- models %>%
      dplyr::mutate(label = dplyr::coalesce(color_labels[use_case_id], use_case_id))
  } else {
    effects$label <- effects$use_case_id
    models$label <- models$use_case_id
  }

  # Coverage boxplot (from effects)
  p_coverage <- ggplot2::ggplot(effects, ggplot2::aes(x = factor(p_star), y = coverage, fill = label)) +
    ggplot2::geom_boxplot(position = ggplot2::position_dodge(width = 0.8), outlier.size = 0.5) +
    ggplot2::labs(title = "Coverage", x = "p* (# causal effects)", y = "Coverage", fill = "Model") +
    ggplot2::theme_minimal() +
    ggplot2::theme(legend.position = "none")

  # Size boxplot (from effects)
  p_size <- ggplot2::ggplot(effects, ggplot2::aes(x = factor(p_star), y = size, fill = label)) +
    ggplot2::geom_boxplot(position = ggplot2::position_dodge(width = 0.8), outlier.size = 0.5) +
    ggplot2::labs(title = "CS Size", x = "p* (# causal effects)", y = "Size", fill = "Model") +
    ggplot2::theme_minimal() +
    ggplot2::theme(legend.position = "none")

  # Purity boxplot (from effects)
  p_purity <- ggplot2::ggplot(effects, ggplot2::aes(x = factor(p_star), y = purity, fill = label)) +
    ggplot2::geom_boxplot(position = ggplot2::position_dodge(width = 0.8), outlier.size = 0.5) +
    ggplot2::labs(title = "Purity", x = "p* (# causal effects)", y = "Purity", fill = "Model") +
    ggplot2::theme_minimal() +
    ggplot2::theme(legend.position = "none")

  # Power boxplot (from models)
  p_power <- ggplot2::ggplot(models, ggplot2::aes(x = factor(p_star), y = power, fill = label)) +
    ggplot2::geom_boxplot(position = ggplot2::position_dodge(width = 0.8), outlier.size = 0.5) +
    ggplot2::labs(title = "Power", x = "p* (# causal effects)", y = "Power", fill = "Model") +
    ggplot2::theme_minimal() +
    ggplot2::theme(legend.position = "bottom")

  # Extract legend from p_power for shared legend
  legend <- cowplot::get_legend(p_power)

  # Remove legends from individual plots
  p_coverage <- p_coverage + ggplot2::theme(legend.position = "none")
  p_size <- p_size + ggplot2::theme(legend.position = "none")
  p_purity <- p_purity + ggplot2::theme(legend.position = "none")
  p_power_no_legend <- p_power + ggplot2::theme(legend.position = "none")

  # Combine with cowplot
  plot_grid <- cowplot::plot_grid(
    p_coverage, p_size,
    p_purity, p_power_no_legend,
    ncol = 2, nrow = 2
  )

  # Add title
  title <- cowplot::ggdraw() +
    cowplot::draw_label(
      "Credible Set Metrics by Number of Causal Effects",
      fontface = "bold", hjust = 0.5
    )

  # Combine title, plots, and legend
  combined <- cowplot::plot_grid(
    title, plot_grid, legend,
    ncol = 1, rel_heights = c(0.05, 0.85, 0.1)
  )

  combined
}

#' Null-coalescing operator
#' @keywords internal
`%||%` <- function(a, b) if (is.null(a)) b else a
