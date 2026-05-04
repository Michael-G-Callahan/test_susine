## annotation_calibration_figure.R
##
## Regenerates the three annotation-quality calibration heatmaps for the
## SuSiNE manuscript (equivalent causality auROC, signed Spearman,
## unsigned Spearman) with journal-ready styling and highlight overlays
## marking (i) the baseline simulation grid used downstream and (ii) the
## working calibration point anchored against AlphaGenome.
##
## Re-uses the calibration logic from
## vignettes/one_off_validations/annotation_quality_alphagenome_calibration.Rmd.
## The pooled summary table is cached at
## vignettes/one_off_validations/annotation_calibration_cache.rds so the
## expensive annotation simulation only has to run once.
##
## Outputs (overwrites):
##   Writings/plots/annotation_calibration_top_row.png  (panels a + b side by side)
##   Writings/plots/annotation_calibration_panel_c.png  (panel c alone)
##
## The LaTeX figure block places these in a 2x2 layout: the top row image
## spans both columns; panel c sits in the bottom-left and the caption
## fills the bottom-right slot via a minipage.

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(tibble)
  library(tidyr)
  library(purrr)
  library(scales)
  library(viridisLite)
  library(devtools)
  library(patchwork)
})

pkg_root <- "c:/Users/mgcal/OneDrive/Documents/School/Research/Genetics/cTWAS Generalization/Code/test_susine"
plots_dir <- file.path(
  "c:/Users/mgcal/OneDrive/Documents/School/Research/Genetics",
  "cTWAS Generalization/Writings/plots/annotation_calibration"
)
cache_path <- file.path(pkg_root, "vignettes/one_off_validations/annotation_calibration_cache.rds")

dir.create(plots_dir, recursive = TRUE, showWarnings = FALSE)

devtools::load_all(pkg_root, quiet = TRUE)

## --------------------------------------------------------------------- ##
## Options (mirror the Rmd)                                              ##
## --------------------------------------------------------------------- ##

p                   <- 1000L
n_beta_draws        <- 600L
p_score_causality   <- 3L
p_score_rank        <- 23L
annotation_r2_grid  <- seq(0, 1, by = 0.1)
inflate_match_grid  <- seq(0, 1, by = 0.1)
effect_sd           <- 1

## Baseline simulation grid (carried into Section 2.3 of the paper)
baseline_phi_vals <- c(0.0, 0.3, 0.5)
baseline_nu_vals  <- c(0.8, 0.9, 1.0)

## Working calibration point against AlphaGenome
working_phi_vals <- c(0.3)
working_nu_vals  <- c(0.9, 1.0)

fmt_num <- function(x) formatC(x, format = "f", digits = 2)
fmt_axis <- function(x) formatC(x, format = "f", digits = 1)

## --------------------------------------------------------------------- ##
## Calibration helpers (copied verbatim from the Rmd)                    ##
## --------------------------------------------------------------------- ##

`%||%` <- function(x, y) if (is.null(x)) y else x

top_k_causal_idx <- function(beta, k) {
  beta <- as.numeric(beta)
  causal_idx <- which(beta != 0)
  if (!length(causal_idx)) return(integer(0))
  n_top <- min(as.integer(k), length(causal_idx))
  causal_idx[order(abs(beta[causal_idx]), decreasing = TRUE)[seq_len(n_top)]]
}

simulate_priors_zero_safe <- function(beta, annotation_r2, inflate_match,
                                      base_sigma2 = NULL, effect_sd = NULL,
                                      seed = NULL) {
  if (!is.null(seed) && is.finite(seed)) set.seed(as.integer(seed))
  p <- length(beta)
  causal_idx <- which(beta != 0)
  noncausal_idx <- setdiff(seq_len(p), causal_idx)
  if (is.null(effect_sd) || !is.finite(effect_sd) || effect_sd <= 0)
    effect_sd <- sqrt(max(stats::var(beta[causal_idx]), stats::var(beta), 1))
  if (is.null(annotation_r2) || is.na(annotation_r2)) annotation_r2 <- 0
  if (is.null(inflate_match) || is.na(inflate_match))  inflate_match  <- 0
  annotation_r2 <- pmin(pmax(annotation_r2, 0), 1)
  inflate_match <- max(inflate_match, 0)
  causal_var <- stats::var(beta[causal_idx])
  if (is.na(causal_var) || causal_var == 0) causal_var <- effect_sd^2
  mu_0 <- numeric(p)
  if (length(causal_idx)) {
    causal_beta <- as.numeric(beta[causal_idx])
    causal_mean <- mean(causal_beta)
    causal_centered <- causal_beta - causal_mean
    causal_scale <- sqrt(mean(causal_centered^2))
    if (!is.finite(causal_scale) || causal_scale <= 0) causal_scale <- sqrt(causal_var)
    if (length(causal_idx) == 1L || !is.finite(causal_scale) || causal_scale <= 0) {
      mu_0[causal_idx] <- causal_beta
    } else {
      signal_unit <- causal_centered / causal_scale
      noise_unit  <- draw_centered_noise(length(causal_idx), target_scale = 1,
                                         orthogonal_to = signal_unit)
      causal_annotation_centered <- causal_scale * (
        sqrt(annotation_r2) * signal_unit + sqrt(1 - annotation_r2) * noise_unit)
      mu_0[causal_idx] <- causal_mean + causal_annotation_centered
    }
  }
  theoretical_causal_mu0_var <- if (length(causal_idx) > 1L)
    mean((mu_0[causal_idx] - mean(mu_0[causal_idx]))^2) else causal_var
  if (length(noncausal_idx)) {
    if (inflate_match <= 0) {
      mu_0[noncausal_idx] <- 0
    } else {
      noncausal_scale <- sqrt(inflate_match) * sqrt(theoretical_causal_mu0_var)
      mu_0[noncausal_idx] <- draw_centered_noise(length(noncausal_idx),
                                                 target_scale = noncausal_scale)
    }
  }
  list(mu_0 = mu_0)
}

pooled_spearman <- function(annotation_list, beta_list, positive_idx_list, use_abs = FALSE) {
  pooled_a <- unlist(purrr::map2(annotation_list, positive_idx_list,
                                 ~ as.numeric(.x)[.y]), use.names = FALSE)
  pooled_b <- unlist(purrr::map2(beta_list, positive_idx_list,
                                 ~ as.numeric(.x)[.y]), use.names = FALSE)
  if (use_abs) { pooled_a <- abs(pooled_a); pooled_b <- abs(pooled_b) }
  if (length(pooled_a) < 2L || length(unique(pooled_a)) < 2L ||
      length(unique(pooled_b)) < 2L) return(NA_real_)
  suppressWarnings(stats::cor(pooled_a, pooled_b, method = "spearman"))
}

pooled_average_precision <- function(annotation_list, positive_idx_list, null_idx_list) {
  pos <- unlist(purrr::map2(annotation_list, positive_idx_list,
                            ~ abs(as.numeric(.x)[.y])), use.names = FALSE)
  neg <- unlist(purrr::map2(annotation_list, null_idx_list,
                            ~ abs(as.numeric(.x)[.y])), use.names = FALSE)
  scores <- c(pos, neg)
  labels <- c(rep.int(1L, length(pos)), rep.int(0L, length(neg)))
  if (!length(scores) || !sum(labels)) return(NA_real_)
  auprc_average_precision(scores = scores, labels = labels)
}

ap_from_dprime <- function(dprime, prevalence, lo = -8, hi = 12, n = 4000L) {
  thr <- seq(lo, hi, length.out = n + 1L)
  tpr <- 1 - stats::pnorm(thr, mean = dprime, sd = 1)
  fpr <- 1 - stats::pnorm(thr, mean = 0,      sd = 1)
  prec <- prevalence * tpr / (prevalence * tpr + (1 - prevalence) * fpr)
  prec[is.nan(prec)] <- 1
  sum(0.5 * (prec[-1] + prec[-length(prec)]) * abs(diff(tpr)))
}

build_ap_to_balanced_auroc <- function(prevalence, d_max = 6,
                                       d_grid_length = 1201L,
                                       integration_points = 4000L) {
  d_grid <- seq(0, d_max, length.out = d_grid_length)
  ap_grid <- vapply(d_grid, ap_from_dprime, numeric(1),
                    prevalence = prevalence, n = integration_points)
  ap_grid <- cummax(ap_grid)
  function(observed_ap) {
    if (!is.finite(observed_ap) || observed_ap <= prevalence) return(0.5)
    clamped_ap <- min(observed_ap, max(ap_grid))
    dprime <- stats::approx(x = ap_grid, y = d_grid, xout = clamped_ap,
                            rule = 2, ties = "ordered")$y
    stats::pnorm(dprime / sqrt(2))
  }
}

## --------------------------------------------------------------------- ##
## Run simulation (with cache)                                           ##
## --------------------------------------------------------------------- ##

if (file.exists(cache_path)) {
  cat("Loading cached pooled summary from", cache_path, "\n")
  pooled_summary_tbl <- readRDS(cache_path)
} else {
  cat("Running calibration simulation (600 loci x 121 grid points) ...\n")
  t0 <- Sys.time()

  beta_tbl <- tibble(beta_id   = seq_len(n_beta_draws),
                     beta_seed = 10000L + seq_len(n_beta_draws)) %>%
    mutate(effects = purrr::map(.data$beta_seed,
                                ~ simulate_effect_sizes_susie2_oligogenic(p = p, seed = .x)),
           beta        = purrr::map(.data$effects, "beta"),
           effect_tier = purrr::map(.data$effects, "effect_tier"),
           causal_idx  = purrr::map(.data$effects, "causal_idx"),
           positive_idx_causality = purrr::map(.data$beta, top_k_causal_idx, k = p_score_causality),
           positive_idx_rank      = purrr::map(.data$beta, top_k_causal_idx, k = p_score_rank),
           null_idx               = purrr::map(.data$effect_tier, ~ which(.x == "none")),
           n_scored_causality     = purrr::map_int(.data$positive_idx_causality, length),
           n_null                 = purrr::map_int(.data$null_idx, length)) %>%
    select(-.data$effects)

  settings_tbl <- expand_grid(annotation_r2 = annotation_r2_grid,
                              inflate_match  = inflate_match_grid) %>%
    arrange(.data$annotation_r2, .data$inflate_match) %>%
    mutate(setting_id      = row_number(),
           annotation_r2_f = factor(fmt_num(.data$annotation_r2),
                                    levels = fmt_num(annotation_r2_grid)),
           inflate_match_f = factor(fmt_num(.data$inflate_match),
                                    levels = fmt_num(inflate_match_grid)))

  annotation_tbl <- expand_grid(beta_id    = beta_tbl$beta_id,
                                setting_id = settings_tbl$setting_id) %>%
    left_join(beta_tbl,     by = "beta_id") %>%
    left_join(settings_tbl, by = "setting_id") %>%
    mutate(annotation_seed = 500000L + .data$beta_id * 1000L + .data$setting_id,
           annotation = purrr::pmap(list(.data$beta, .data$annotation_r2,
                                         .data$inflate_match, .data$annotation_seed),
             function(beta, annotation_r2, inflate_match, seed) {
               simulate_priors_zero_safe(beta, annotation_r2, inflate_match,
                                         base_sigma2 = effect_sd^2,
                                         effect_sd = effect_sd, seed = seed)$mu_0
             })) %>%
    select(.data$beta_id, .data$setting_id, .data$annotation_r2, .data$inflate_match,
           .data$annotation_r2_f, .data$inflate_match_f, .data$beta, .data$annotation,
           .data$positive_idx_causality, .data$positive_idx_rank,
           .data$null_idx, .data$n_scored_causality, .data$n_null)

  scored_per_draw <- annotation_tbl$n_scored_causality[[1]]
  null_per_draw   <- annotation_tbl$n_null[[1]]
  scored_prev     <- scored_per_draw / (scored_per_draw + null_per_draw)
  ap_to_auroc     <- build_ap_to_balanced_auroc(prevalence = scored_prev)

  pooled_summary_tbl <- annotation_tbl %>%
    group_by(.data$setting_id, .data$annotation_r2, .data$inflate_match,
             .data$annotation_r2_f, .data$inflate_match_f) %>%
    group_modify(function(.x, .y) {
      ap <- pooled_average_precision(annotation_list = .x$annotation,
                                     positive_idx_list = .x$positive_idx_causality,
                                     null_idx_list = .x$null_idx)
      tibble(pooled_auprc = ap,
             equivalent_causality_auroc = ap_to_auroc(ap),
             signed_spearman = pooled_spearman(.x$annotation, .x$beta,
                                               .x$positive_idx_rank, use_abs = FALSE),
             unsigned_spearman = pooled_spearman(.x$annotation, .x$beta,
                                                 .x$positive_idx_rank, use_abs = TRUE))
    }) %>% ungroup()

  saveRDS(pooled_summary_tbl, cache_path)
  cat(sprintf("Cached pooled summary to %s (%.1f min)\n",
              cache_path, as.numeric(difftime(Sys.time(), t0, units = "mins"))))
}

## --------------------------------------------------------------------- ##
## Plot helper with highlight overlays                                   ##
## --------------------------------------------------------------------- ##

## Grid index lookup: factor position -> integer (1..11)
phi_idx <- function(v) match(fmt_num(v), fmt_num(annotation_r2_grid))
nu_idx  <- function(v) match(fmt_num(v), fmt_num(inflate_match_grid))

baseline_rects <- expand_grid(phi = baseline_phi_vals, nu = baseline_nu_vals) %>%
  mutate(x = phi_idx(.data$phi), y = nu_idx(.data$nu),
         xmin = .data$x - 0.5, xmax = .data$x + 0.5,
         ymin = .data$y - 0.5, ymax = .data$y + 0.5)

working_rects <- expand_grid(phi = working_phi_vals, nu = working_nu_vals) %>%
  mutate(x = phi_idx(.data$phi), y = nu_idx(.data$nu),
         xmin = .data$x - 0.5, xmax = .data$x + 0.5,
         ymin = .data$y - 0.5, ymax = .data$y + 0.5)

ideal_text_color <- function(fill_hex) {
  rgb <- grDevices::col2rgb(fill_hex)
  lum <- 0.2126 * rgb[1, ] + 0.7152 * rgb[2, ] + 0.0722 * rgb[3, ]
  ifelse(lum < 140, "white", "black")
}

plot_calibration_heatmap <- function(tbl, value_col, fill_label,
                                     palette_option, limits, digits = 2,
                                     panel_tag, show_x = TRUE, show_y = TRUE) {
  values <- tbl[[value_col]]
  lims   <- if (is.null(limits)) range(values, na.rm = TRUE) else limits
  palette_values <- viridisLite::viridis(256, option = palette_option, end = 0.98)
  fill_mapper    <- scales::col_numeric(palette = palette_values,
                                        domain = lims, na.color = "#bdbdbd")
  plot_tbl <- tbl %>%
    mutate(metric_value = .data[[value_col]],
           fill_hex    = fill_mapper(.data$metric_value),
           text_color  = ideal_text_color(.data$fill_hex),
           label       = sprintf(paste0("%.", digits, "f"), .data$metric_value),
           x = as.integer(.data$annotation_r2_f),
           y = as.integer(.data$inflate_match_f))

  g <- ggplot(plot_tbl, aes(x = .data$x, y = .data$y)) +
    geom_tile(aes(fill = .data$metric_value),
              colour = "white", linewidth = 0.3) +
    geom_text(aes(label = .data$label, colour = .data$text_color),
              size = 2.3, fontface = "bold") +
    geom_rect(data = baseline_rects,
              aes(xmin = .data$xmin, xmax = .data$xmax,
                  ymin = .data$ymin, ymax = .data$ymax),
              inherit.aes = FALSE,
              colour = "black", fill = NA, linewidth = 0.9, linetype = "solid") +
    geom_rect(data = working_rects,
              aes(xmin = .data$xmin, xmax = .data$xmax,
                  ymin = .data$ymin, ymax = .data$ymax),
              inherit.aes = FALSE,
              colour = "#D7261E", fill = NA, linewidth = 1.6) +
    scale_colour_identity() +
    scale_fill_gradientn(colours = palette_values, limits = lims,
                         oob = scales::squish, name = fill_label,
                         guide = guide_colourbar(barwidth = unit(2.6, "mm"),
                                                 barheight = unit(34, "mm"),
                                                 ticks.colour = "grey30",
                                                 title.position = "top",
                                                 title.hjust = 0)) +
    scale_x_continuous(breaks = seq_along(annotation_r2_grid),
                       labels = fmt_axis(annotation_r2_grid), expand = c(0, 0)) +
    scale_y_continuous(breaks = seq_along(inflate_match_grid),
                       labels = fmt_axis(inflate_match_grid), expand = c(0, 0)) +
    coord_equal(clip = "off") +
    labs(x = if (show_x) expression(paste("Annotation accuracy  ", phi[a])) else NULL,
         y = if (show_y) expression(paste("Non-causal inflation  ", nu[a])) else NULL,
         title = panel_tag) +
    theme_minimal(base_size = 9) +
    theme(panel.grid    = element_blank(),
          plot.title    = element_text(size = 10, face = "bold", hjust = 0,
                                       margin = margin(0, 0, 3, 0)),
          axis.title    = element_text(size = 9),
          axis.text     = element_text(size = 8, colour = "grey20"),
          axis.ticks    = element_blank(),
          legend.title  = element_text(size = 8),
          legend.text   = element_text(size = 7.5),
          legend.margin = margin(0, 0, 0, 3),
          plot.margin   = margin(2, 4, 2, 2))
  if (!show_x) g <- g + theme(axis.text.x = element_blank())
  if (!show_y) g <- g + theme(axis.text.y = element_blank())
  g
}

## --------------------------------------------------------------------- ##
## Write the three PNGs                                                   ##
## --------------------------------------------------------------------- ##

p_a <- plot_calibration_heatmap(pooled_summary_tbl,
                                value_col = "equivalent_causality_auroc",
                                fill_label = "auROC (eq.)",
                                palette_option = "D", limits = c(0.5, 1.0),
                                panel_tag = "(a) Equivalent causality auROC",
                                show_x = TRUE, show_y = TRUE)

p_b <- plot_calibration_heatmap(pooled_summary_tbl,
                                value_col = "signed_spearman",
                                fill_label = expression(rho[s]),
                                palette_option = "C", limits = c(-1.0, 1.0),
                                panel_tag = "(b) Signed Spearman",
                                show_x = TRUE, show_y = FALSE)

p_c <- plot_calibration_heatmap(pooled_summary_tbl,
                                value_col = "unsigned_spearman",
                                fill_label = expression(rho[u]),
                                palette_option = "C", limits = c(0.0, 1.0),
                                panel_tag = "(c) Unsigned Spearman",
                                show_x = TRUE, show_y = TRUE)

top_row <- p_a | p_b

out_a <- file.path(plots_dir, "annotation_calibration_panel_a.png")
ggsave(out_a, p_a, width = 3.5, height = 3.6, dpi = 300, bg = "white")
cat("Wrote", out_a, "\n")

out_b <- file.path(plots_dir, "annotation_calibration_panel_b.png")
ggsave(out_b, p_b, width = 3.5, height = 3.6, dpi = 300, bg = "white")
cat("Wrote", out_b, "\n")

out_top <- file.path(plots_dir, "annotation_calibration_top_row.png")
ggsave(out_top, top_row, width = 7.0, height = 3.6, dpi = 300, bg = "white")
cat("Wrote", out_top, "\n")

out_c <- file.path(plots_dir, "annotation_calibration_panel_c.png")
ggsave(out_c, p_c, width = 3.5, height = 3.6, dpi = 300, bg = "white")
cat("Wrote", out_c, "\n")

cat("\nDone.\n")
