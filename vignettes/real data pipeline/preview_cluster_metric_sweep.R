# Preview sweep: cluster metric x threshold for the real-data ensemble.
#
# For each (metric, threshold) it re-collects the real-data ensemble from the
# saved raw_fits using Method-B aggregation (frequency-free, within-cluster
# ELBO-softmax) and the chosen clustering distance, writes a full aggregated/
# set into a preview dir, then re-renders the paper figures into a parallel
# preview figure dir. The published 52906940 outputs and Writings/plots are
# never touched (separate dirs + env-gated redirects; production defaults are
# unchanged because the options/args below are opt-in).
#
# Run on HPC from the test_susine repo root, e.g.:
#   Rscript "vignettes/real data pipeline/preview_cluster_metric_sweep.R"
# Paste the final cluster-count + PVE tables back for calibration.

suppressMessages({
  library(here)
  library(dplyr)
  library(readr)
  library(tibble)
  library(rmarkdown)
})
here::i_am("vignettes/real data pipeline/preview_cluster_metric_sweep.R")

devtools::load_all(here("..", "susine"))
devtools::load_all(here())

job_name <- "real_data_ensemble_geometric_n20"
parent   <- "52906940"
out_root <- here("output")
viz_rmd  <- here("vignettes", "real data pipeline",
                 "visualize_results_workbook_real_data_ensemble.Rmd")

# First-pass thresholds per metric (recalibrate from the cluster-count table).
#  cosine        d = 1 - cos: 0.01/0.05/0.10 = direction agreement >= 99/95/90%
#  l2            raw Euclidean; a full confident swap ~ sqrt(2) ~ 1.41
#  bernoulli_jsd sum of per-variant binary JSDs; full swap ~ 2*ln2 ~ 1.39
grid <- list(
  cosine        = c(0.01, 0.05, 0.10),
  l2            = c(0.30, 0.60, 0.90),
  bernoulli_jsd = c(0.30, 0.60, 1.00)
)

# Method B everywhere in the sweep (frequency-free cluster aggregation).
options(rd_preview.cluster_aggregation = "method_b")

cluster_rows <- list()
pve_rows     <- list()

for (metric in names(grid)) {
  for (thr in grid[[metric]]) {
    tag      <- sprintf("%s_t%s", metric, format(thr, trim = TRUE))
    base_dir <- file.path(out_root, "slurm_output", job_name, parent, "preview", tag)
    agg_dir  <- file.path(base_dir, "aggregated")
    fig_dir  <- file.path(base_dir, "figures")
    dir.create(agg_dir, recursive = TRUE, showWarnings = FALSE)
    dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

    message(sprintf("[%s] collecting -> %s", tag, agg_dir))
    options(rd_preview.cluster_metric = metric)
    collect_real_data_results(
      job_name, parent,
      output_root       = out_root,
      validate          = FALSE,
      output_dir        = agg_dir,
      cluster_threshold = thr
    )

    # Cluster-count diagnostic (per-locus, then summarized).
    mm_path <- file.path(agg_dir, "multimodal_metrics.csv")
    if (file.exists(mm_path)) {
      mm <- readr::read_csv(mm_path, show_col_types = FALSE)
      cluster_rows[[tag]] <- tibble::tibble(
        metric = metric, threshold = thr,
        mean_n_clusters   = mean(mm$n_clusters, na.rm = TRUE),
        median_n_clusters = stats::median(mm$n_clusters, na.rm = TRUE),
        max_n_clusters    = max(mm$n_clusters, na.rm = TRUE)
      )
    }

    # PVE diagnostic (anchor vs ensemble vs warm-refit, averaged over loci).
    ps_path <- file.path(agg_dir, "paper_real_data_ensemble_summary.csv")
    if (file.exists(ps_path)) {
      ps <- readr::read_csv(ps_path, show_col_types = FALSE)
      pick <- function(nm) if (nm %in% names(ps)) mean(ps[[nm]], na.rm = TRUE) else NA_real_
      pve_rows[[tag]] <- tibble::tibble(
        metric = metric, threshold = thr,
        pve_anchor   = pick("pve_susie_anchor"),
        pve_ensemble = pick("pve_susine_ensemble"),
        pve_refit    = pick("pve_highest_weight_refit")
      )
    }

    # Re-render the paper figures into the preview fig dir (env-gated redirect).
    Sys.setenv(RD_PREVIEW_AGG_DIR   = agg_dir,
               RD_PREVIEW_FIG_DIR   = fig_dir,
               RD_PREVIEW_PAPER_DIR = fig_dir)
    message(sprintf("[%s] rendering figures -> %s", tag, fig_dir))
    ok <- tryCatch({
      rmarkdown::render(viz_rmd, output_dir = fig_dir,
                        output_file = sprintf("viz_%s.html", tag), quiet = TRUE)
      TRUE
    }, error = function(e) { message("  render failed: ", conditionMessage(e)); FALSE })
    Sys.unsetenv(c("RD_PREVIEW_AGG_DIR", "RD_PREVIEW_FIG_DIR", "RD_PREVIEW_PAPER_DIR"))
  }
}

options(rd_preview.cluster_metric = NULL, rd_preview.cluster_aggregation = NULL)

cat("\n=== cluster-count by (metric, threshold) ===\n")
print(as.data.frame(dplyr::bind_rows(cluster_rows)), row.names = FALSE)
cat("\n=== mean PVE proxy by (metric, threshold) ===\n")
print(as.data.frame(dplyr::bind_rows(pve_rows)), row.names = FALSE)
cat(sprintf(
  "\nFigures + figures HTML under: output/slurm_output/%s/%s/preview/<tag>/figures/\n",
  job_name, parent))
