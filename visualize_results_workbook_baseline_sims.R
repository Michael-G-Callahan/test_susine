## ----setup, message=FALSE-----------------------------------------------------
here::i_am("vignettes/simulation pipeline/visualize_results_workbook_baseline_sims.Rmd")
library(here)
library(devtools)

susine_path <- here("..", "susine")
if (file.exists(file.path(susine_path, "DESCRIPTION"))) {
  options(test_susine.local_susine_path = normalizePath(susine_path, winslash = "/", mustWork = TRUE))
  devtools::load_all(susine_path)
}
devtools::load_all(".")

library(dplyr)
library(readr)
library(tidyr)
library(ggplot2)
library(knitr)

theme_set(theme_bw(base_size = 12))


## ----load---------------------------------------------------------------------
job_name <- "baseline_sims_screen"
parent_job_id <- "REPLACE_ME"
output_root <- here("output")

job_dir <- file.path(output_root, "slurm_output", job_name, parent_job_id)
consolidated_dir <- file.path(job_dir, "consolidated")
fig_root <- file.path(job_dir, "figures")

fig_dirs <- c(
  "susie2_plots",
  "baseline_performance",
  "dataset_difficulty",
  "effect_diagnostics"
)
for (d in fig_dirs) {
  dir.create(file.path(fig_root, d), recursive = TRUE, showWarnings = FALSE)
}

read_req <- function(filename) {
  readr::read_csv(file.path(consolidated_dir, filename), show_col_types = FALSE)
}

read_opt <- function(filename) {
  path <- file.path(consolidated_dir, filename)
  if (!file.exists(path)) {
    return(NULL)
  }
  readr::read_csv(path, show_col_types = FALSE)
}

model_metrics <- read_req("model_metrics_full.csv")
effect_metrics <- read_opt("effect_metrics_full.csv")
effect_metrics_unfiltered <- read_opt("effect_metrics_unfiltered_full.csv")
confusion_bins <- read_req("confusion_bins_full.csv")
dataset_metrics <- read_opt("dataset_metrics_full.csv")
tier_cs_metrics <- read_opt("tier_cs_metrics_full.csv")
auprc_per_dataset <- read_req("auprc_per_dataset.csv")
auprc_pooled_overall <- read_req("auprc_pooled_overall.csv")
auprc_pooled_by_annotation <- read_req("auprc_pooled_by_annotation.csv")
tpr05_per_dataset <- read_req("tpr05_per_dataset.csv")
tpr05_pooled_overall <- read_req("tpr05_pooled_overall.csv")
tpr05_pooled_by_annotation <- read_req("tpr05_pooled_by_annotation.csv")
screening_summary <- read_req("screening_summary.csv")

cat("Figure root:", fig_root, "\n")
cat("Pooled AUPRC settings:", nrow(auprc_pooled_overall), "\n")
cat("Unfiltered effect rows:", if (is.null(effect_metrics_unfiltered)) 0 else nrow(effect_metrics_unfiltered), "\n")


## ----helpers------------------------------------------------------------------
fmt_num <- function(x, digits = 3) {
  format(signif(as.numeric(x), digits), trim = TRUE, scientific = FALSE)
}

top_settings_by_family <- screening_summary %>%
  dplyr::group_by(.data$method_family) %>%
  dplyr::slice_max(.data$AUPRC, n = 3L, with_ties = FALSE) %>%
  dplyr::ungroup()

best_settings_by_family <- screening_summary %>%
  dplyr::group_by(.data$method_family) %>%
  dplyr::slice_max(.data$AUPRC, n = 1L, with_ties = FALSE) %>%
  dplyr::ungroup()


## ----tpr-fpr, fig.width=12, fig.height=8--------------------------------------
plot_settings <- top_settings_by_family %>%
  dplyr::select(.data$setting_label, .data$annotation_label)

roc_df <- confusion_bins %>%
  dplyr::semi_join(plot_settings, by = c("setting_label", "annotation_label")) %>%
  dplyr::group_by(
    .data$method_family,
    .data$setting_label,
    .data$annotation_label,
    .data$pip_threshold
  ) %>%
  dplyr::summarise(
    n_causal_at_bucket = sum(.data$n_causal_at_bucket, na.rm = TRUE),
    n_noncausal_at_bucket = sum(.data$n_noncausal_at_bucket, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  dplyr::group_by(.data$method_family, .data$setting_label, .data$annotation_label) %>%
  dplyr::arrange(dplyr::desc(.data$pip_threshold), .by_group = TRUE) %>%
  dplyr::mutate(
    cum_tp = cumsum(.data$n_causal_at_bucket),
    cum_fp = cumsum(.data$n_noncausal_at_bucket),
    total_p = sum(.data$n_causal_at_bucket),
    total_n = sum(.data$n_noncausal_at_bucket),
    tpr = dplyr::if_else(.data$total_p > 0, .data$cum_tp / .data$total_p, NA_real_),
    fpr = dplyr::if_else(.data$total_n > 0, .data$cum_fp / .data$total_n, NA_real_)
  ) %>%
  dplyr::ungroup()

p_roc <- ggplot(roc_df, aes(x = .data$fpr, y = .data$tpr, color = .data$setting_label)) +
  geom_line(linewidth = 0.7) +
  facet_wrap(~method_family, scales = "free_y") +
  labs(
    x = "FPR",
    y = "TPR",
    color = "Setting",
    title = "TPR vs FPR for top baseline settings by method family"
  ) +
  theme(legend.position = "bottom")

ggsave(
  file.path(fig_root, "susie2_plots", "tpr_fpr_top_settings.png"),
  p_roc,
  width = 12,
  height = 8,
  dpi = 150
)


## ----pip-calibration, fig.width=12, fig.height=8------------------------------
cal_df <- confusion_bins %>%
  dplyr::semi_join(plot_settings, by = c("setting_label", "annotation_label")) %>%
  dplyr::mutate(pip_bin = floor(.data$pip_threshold / 0.10) * 0.10) %>%
  dplyr::group_by(
    .data$method_family,
    .data$setting_label,
    .data$annotation_label,
    .data$dataset_bundle_id,
    .data$pip_bin
  ) %>%
  dplyr::summarise(
    n_causal = sum(.data$n_causal_at_bucket, na.rm = TRUE),
    n_noncausal = sum(.data$n_noncausal_at_bucket, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  dplyr::mutate(
    n_total = .data$n_causal + .data$n_noncausal,
    observed_frac = dplyr::if_else(.data$n_total > 0, .data$n_causal / .data$n_total, NA_real_)
  ) %>%
  dplyr::group_by(.data$method_family, .data$setting_label, .data$annotation_label, .data$pip_bin) %>%
  dplyr::summarise(
    mean_frac = mean(.data$observed_frac, na.rm = TRUE),
    se_frac = sd(.data$observed_frac, na.rm = TRUE) / sqrt(sum(!is.na(.data$observed_frac))),
    .groups = "drop"
  ) %>%
  dplyr::mutate(
    ci_lo = pmax(0, .data$mean_frac - 1.96 * .data$se_frac),
    ci_hi = pmin(1, .data$mean_frac + 1.96 * .data$se_frac)
  )

p_cal <- ggplot(cal_df, aes(x = .data$pip_bin + 0.05, y = .data$mean_frac, color = .data$setting_label)) +
  geom_errorbar(aes(ymin = .data$ci_lo, ymax = .data$ci_hi), width = 0.02, alpha = 0.5) +
  geom_point(size = 1.8) +
  geom_line() +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey50") +
  facet_wrap(~method_family) +
  labs(
    x = "Mean PIP bin midpoint",
    y = "Observed fraction causal",
    color = "Setting",
    title = "PIP calibration for top baseline settings"
  ) +
  coord_cartesian(xlim = c(0, 1), ylim = c(0, 1)) +
  theme(legend.position = "bottom")

ggsave(
  file.path(fig_root, "susie2_plots", "pip_calibration_top_settings.png"),
  p_cal,
  width = 12,
  height = 8,
  dpi = 150
)


## ----cs-power-fdr, fig.width=12, fig.height=10--------------------------------
if (!is.null(tier_cs_metrics) && nrow(tier_cs_metrics)) {
  tier_summary <- tier_cs_metrics %>%
    dplyr::semi_join(plot_settings, by = c("setting_label", "annotation_label")) %>%
    dplyr::group_by(.data$method_family, .data$setting_label, .data$n_top) %>%
    dplyr::summarise(
      sum_causal_covered = sum(.data$n_causal_covered, na.rm = TRUE),
      sum_causal_total = sum(.data$n_causal_total, na.rm = TRUE),
      sum_false_cs = sum(.data$n_false_cs, na.rm = TRUE),
      sum_cs_total = sum(.data$n_cs_total, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    dplyr::mutate(
      cs_power = dplyr::if_else(.data$sum_causal_total > 0, .data$sum_causal_covered / .data$sum_causal_total, NA_real_),
      cs_fdr = dplyr::if_else(.data$sum_cs_total > 0, .data$sum_false_cs / .data$sum_cs_total, NA_real_)
    )

  p_power <- ggplot(tier_summary, aes(x = factor(.data$n_top), y = .data$cs_power, color = .data$setting_label, group = .data$setting_label)) +
    geom_point() +
    geom_line() +
    facet_wrap(~method_family) +
    labs(x = "Top-N treated as causal", y = "CS power", color = "Setting",
         title = "CS power for top baseline settings") +
    theme(legend.position = "bottom")

  p_fdr <- ggplot(tier_summary, aes(x = factor(.data$n_top), y = .data$cs_fdr, color = .data$setting_label, group = .data$setting_label)) +
    geom_point() +
    geom_line() +
    geom_hline(yintercept = 0.05, linetype = "dashed", color = "red") +
    facet_wrap(~method_family) +
    labs(x = "Top-N treated as causal", y = "CS FDR", color = "Setting",
         title = "CS FDR for top baseline settings") +
    theme(legend.position = "bottom")

  ggsave(file.path(fig_root, "susie2_plots", "cs_power_top_settings.png"), p_power, width = 12, height = 5, dpi = 150)
  ggsave(file.path(fig_root, "susie2_plots", "cs_fdr_top_settings.png"), p_fdr, width = 12, height = 5, dpi = 150)
}


## ----overall-screening, fig.width=14, fig.height=8----------------------------
p_auprc <- auprc_pooled_overall %>%
  dplyr::left_join(
    screening_summary %>% dplyr::select(.data$setting_label, .data$annotation_label, .data$method_family),
    by = c("setting_label", "annotation_label", "method_family")
  ) %>%
  dplyr::mutate(display_label = paste(.data$setting_label, "[", .data$annotation_label, "]")) %>%
  ggplot(aes(x = reorder(.data$display_label, .data$AUPRC), y = .data$AUPRC, fill = .data$method_family)) +
  geom_col() +
  coord_flip() +
  labs(x = "Setting", y = "Pooled AUPRC", fill = "Method family",
       title = "Overall baseline screening by pooled AUPRC")

p_tpr <- tpr05_pooled_overall %>%
  dplyr::left_join(
    screening_summary %>% dplyr::select(.data$setting_label, .data$annotation_label, .data$method_family),
    by = c("setting_label", "annotation_label", "method_family")
  ) %>%
  dplyr::mutate(display_label = paste(.data$setting_label, "[", .data$annotation_label, "]")) %>%
  ggplot(aes(x = reorder(.data$display_label, .data$tpr_fpr05), y = .data$tpr_fpr05, fill = .data$method_family)) +
  geom_col() +
  coord_flip() +
  labs(x = "Setting", y = "TPR@FPR=0.05", fill = "Method family",
       title = "Overall baseline screening by pooled TPR@FPR=0.05")

ggsave(file.path(fig_root, "baseline_performance", "auprc_overall_all_settings.png"), p_auprc, width = 14, height = 8, dpi = 150)
ggsave(file.path(fig_root, "baseline_performance", "tpr05_overall_all_settings.png"), p_tpr, width = 14, height = 8, dpi = 150)


## ----by-annotation, fig.width=14, fig.height=8--------------------------------
p_auprc_ann <- ggplot(
  auprc_pooled_by_annotation,
  aes(x = reorder(.data$setting_label, .data$AUPRC), y = .data$AUPRC, fill = .data$method_family)
) +
  geom_col() +
  coord_flip() +
  facet_wrap(~annotation_label, scales = "free_y") +
  labs(x = "Setting", y = "Pooled AUPRC", fill = "Method family",
       title = "Baseline screening by annotation setting")

p_tpr_ann <- ggplot(
  tpr05_pooled_by_annotation,
  aes(x = reorder(.data$setting_label, .data$tpr_fpr05), y = .data$tpr_fpr05, fill = .data$method_family)
) +
  geom_col() +
  coord_flip() +
  facet_wrap(~annotation_label, scales = "free_y") +
  labs(x = "Setting", y = "TPR@FPR=0.05", fill = "Method family",
       title = "Baseline TPR@FPR=0.05 by annotation setting")

ggsave(file.path(fig_root, "baseline_performance", "auprc_by_annotation.png"), p_auprc_ann, width = 14, height = 8, dpi = 150)
ggsave(file.path(fig_root, "baseline_performance", "tpr05_by_annotation.png"), p_tpr_ann, width = 14, height = 8, dpi = 150)


## ----family-screening, fig.width=12, fig.height=8-----------------------------
vanilla_sweeps <- screening_summary %>%
  dplyr::filter(.data$method_family %in% c("susie_vanilla", "susine_vanilla", "susie_inf", "susie_ash_fixed")) %>%
  dplyr::mutate(sigma_value = as.numeric(.data$sigma_0_2_scalar))

p_vanilla_sigma <- ggplot(vanilla_sweeps, aes(x = .data$sigma_value, y = .data$AUPRC, color = .data$annotation_label)) +
  geom_point() +
  geom_line() +
  facet_wrap(~method_family, scales = "free_y") +
  labs(x = expression(sigma[0]^2 / var(y)), y = "Pooled AUPRC", color = "Annotation",
       title = "Sigma sweeps for vanilla/fixed-sigma baseline families")

tau_sweeps <- screening_summary %>%
  dplyr::filter(.data$method_family == "susie_functional_pi")

p_tau <- ggplot(tau_sweeps, aes(x = .data$tau_value, y = .data$AUPRC, color = .data$annotation_label)) +
  geom_point() +
  geom_line() +
  scale_x_log10() +
  labs(x = "tau", y = "Pooled AUPRC", color = "Annotation",
       title = "Tau sweep for susie_functional_pi")

mu_grid <- screening_summary %>%
  dplyr::filter(.data$method_family == "susine_functional_mu") %>%
  dplyr::mutate(sigma_value = as.numeric(.data$sigma_0_2_scalar))

p_mu_heat <- ggplot(mu_grid, aes(x = .data$c_value, y = .data$sigma_value, fill = .data$AUPRC)) +
  geom_tile() +
  facet_wrap(~annotation_label) +
  labs(x = "c", y = expression(sigma[0]^2 / var(y)), fill = "Pooled AUPRC",
       title = "c x sigma screening heatmap for susine_functional_mu")

single_point_families <- screening_summary %>%
  dplyr::filter(.data$method_family %in% c("susie_eb", "susine_eb_clamped_scale_var_nonneg"))

p_single <- ggplot(single_point_families, aes(x = .data$method_family, y = .data$AUPRC, fill = .data$annotation_label)) +
  geom_col(position = position_dodge()) +
  labs(x = "Method family", y = "Pooled AUPRC", fill = "Annotation",
       title = "Single-point baseline families")

ggsave(file.path(fig_root, "baseline_performance", "sigma_sweeps.png"), p_vanilla_sigma, width = 12, height = 8, dpi = 150)
ggsave(file.path(fig_root, "baseline_performance", "tau_sweep.png"), p_tau, width = 10, height = 6, dpi = 150)
ggsave(file.path(fig_root, "baseline_performance", "susine_functional_mu_heatmap.png"), p_mu_heat, width = 12, height = 8, dpi = 150)
ggsave(file.path(fig_root, "baseline_performance", "single_point_families.png"), p_single, width = 10, height = 6, dpi = 150)


## ----top-settings-table-------------------------------------------------------
top_table <- screening_summary %>%
  dplyr::group_by(.data$method_family) %>%
  dplyr::slice_max(.data$AUPRC, n = 5L, with_ties = FALSE) %>%
  dplyr::ungroup() %>%
  dplyr::select(
    .data$method_family,
    .data$setting_label,
    .data$annotation_label,
    .data$AUPRC,
    .data$tpr_fpr05
  )

readr::write_csv(top_table, file.path(fig_root, "baseline_performance", "top_settings_by_family.csv"))
kable(top_table, digits = 4, caption = "Top pooled settings by method family")


## ----dataset-difficulty, fig.width=12, fig.height=8---------------------------
if (!is.null(dataset_metrics) && nrow(dataset_metrics)) {
  best_run_ids <- auprc_per_dataset %>%
    dplyr::group_by(.data$method_family, .data$dataset_bundle_id) %>%
    dplyr::slice_max(.data$AUPRC, n = 1L, with_ties = FALSE) %>%
    dplyr::ungroup() %>%
    dplyr::select(.data$dataset_bundle_id, .data$method_family, .data$AUPRC, .data$annotation_label)

  difficulty_df <- dataset_metrics %>%
    dplyr::distinct(.data$dataset_bundle_id, .keep_all = TRUE) %>%
    dplyr::inner_join(best_run_ids, by = "dataset_bundle_id")

  p_m1 <- ggplot(difficulty_df, aes(x = .data$M1, y = .data$AUPRC, color = .data$method_family)) +
    geom_point(alpha = 0.7) +
    geom_smooth(se = FALSE) +
    facet_wrap(~annotation_label) +
    labs(x = "M1", y = "Per-dataset AUPRC", color = "Method family",
         title = "Dataset difficulty vs AUPRC for best per-family settings")

  z_metric <- if ("z_eff_signals" %in% names(difficulty_df)) "z_eff_signals" else "z_topk_ratio"
  p_z <- ggplot(difficulty_df, aes(x = .data[[z_metric]], y = .data$AUPRC, color = .data$method_family)) +
    geom_point(alpha = 0.7) +
    geom_smooth(se = FALSE) +
    facet_wrap(~annotation_label) +
    labs(x = z_metric, y = "Per-dataset AUPRC", color = "Method family",
         title = "Z-score difficulty metric vs AUPRC")

  ggsave(file.path(fig_root, "dataset_difficulty", "difficulty_vs_M1.png"), p_m1, width = 12, height = 8, dpi = 150)
  ggsave(file.path(fig_root, "dataset_difficulty", "difficulty_vs_z_metric.png"), p_z, width = 12, height = 8, dpi = 150)
}


## ----effect-diagnostics, fig.width=12, fig.height=8---------------------------
if (!is.null(effect_metrics_unfiltered) && nrow(effect_metrics_unfiltered)) {
  diag_methods <- c("susie_vanilla", "susine_functional_mu", "susie_functional_pi")
  effect_diag <- effect_metrics_unfiltered %>%
    dplyr::filter(.data$method_family %in% diag_methods)

  p_main <- ggplot(
    effect_diag,
    aes(x = .data$effect_k_eff_signal_core95, y = .data$accuracy_ratio, color = .data$method_family)
  ) +
    geom_point(alpha = 0.35, size = 1.2) +
    facet_wrap(~annotation_label) +
    labs(
      x = "effect_k_eff_signal_core95",
      y = "accuracy_ratio",
      color = "Method family",
      title = "Effect diffuseness vs causal-centering accuracy"
    )

  p_tail <- ggplot(
    effect_diag,
    aes(x = .data$tail_inflation_ratio, y = .data$accuracy_ratio, color = .data$method_family)
  ) +
    geom_point(alpha = 0.35, size = 1.2) +
    facet_wrap(~annotation_label) +
    labs(
      x = "tail_inflation_ratio",
      y = "accuracy_ratio",
      color = "Method family",
      title = "Tail inflation ratio vs accuracy ratio"
    )

  p_tail_log <- ggplot(
    effect_diag,
    aes(x = .data$tail_inflation_log, y = .data$accuracy_ratio, color = .data$method_family)
  ) +
    geom_point(alpha = 0.35, size = 1.2) +
    facet_wrap(~annotation_label) +
    labs(
      x = "tail_inflation_log",
      y = "accuracy_ratio",
      color = "Method family",
      title = "Tail inflation log-ratio vs accuracy ratio"
    )

  p_entropy <- ggplot(
    effect_diag,
    aes(x = .data$effect_pip_entropy_core95, y = .data$accuracy_ratio, color = .data$method_family)
  ) +
    geom_point(alpha = 0.35, size = 1.2) +
    facet_wrap(~annotation_label) +
    labs(
      x = "effect_pip_entropy_core95",
      y = "accuracy_ratio",
      color = "Method family",
      title = "Core-95 entropy vs accuracy ratio"
    )

  p_bin2d <- ggplot(
    effect_diag,
    aes(x = .data$effect_k_eff_signal_core95, y = .data$accuracy_ratio)
  ) +
    geom_bin2d() +
    facet_grid(rows = vars(method_family), cols = vars(annotation_label)) +
    scale_fill_viridis_c() +
    labs(
      x = "effect_k_eff_signal_core95",
      y = "accuracy_ratio",
      fill = "Count",
      title = "Effect diffuseness vs accuracy ratio density"
    )

  ggsave(file.path(fig_root, "effect_diagnostics", "effect_keff_core95_vs_accuracy_ratio.png"), p_main, width = 12, height = 8, dpi = 150)
  ggsave(file.path(fig_root, "effect_diagnostics", "tail_inflation_ratio_vs_accuracy_ratio.png"), p_tail, width = 12, height = 8, dpi = 150)
  ggsave(file.path(fig_root, "effect_diagnostics", "tail_inflation_log_vs_accuracy_ratio.png"), p_tail_log, width = 12, height = 8, dpi = 150)
  ggsave(file.path(fig_root, "effect_diagnostics", "effect_entropy_core95_vs_accuracy_ratio.png"), p_entropy, width = 12, height = 8, dpi = 150)
  ggsave(file.path(fig_root, "effect_diagnostics", "effect_keff_core95_vs_accuracy_ratio_bin2d.png"), p_bin2d, width = 14, height = 10, dpi = 150)
}

