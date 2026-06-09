#!/usr/bin/env Rscript
# Downstream recompute of the corrected local-genetic-variance ("heritability")
# estimand for the real-data case studies. NO re-fitting required: the real-data
# pipeline persists every fit (saveRDS ... compress="xz"), so this script reloads
# the saved fits, pulls the per-locus LD matrix R, and computes
#
#   hg2_postmean     = m' R m / var(y)            (var(y) = 1, standardized RSS)
#   hg2_uncertainty  = tr(R Cov(b|y)) / var(y)    (within-fit correction)
#   hg2_expected_pve = hg2_postmean + hg2_uncertainty   (= local PVE / cis-h2_g)
#
# and, per locus, the fair ensemble aggregate
#   sum_k w_k * hg2_expected_pve_k
# using the SAME cluster weights the pipeline uses (.cluster_weights_from_hc with
# the job's jsd_threshold / softmax_temperature). See
# refs/decisions/heritability_estimand_decision_2026-06-09.md (sec 5.5).
#
# Run on the HPC where output/ and data/ live. Example:
#   Rscript inst/scripts/recompute_real_data_hg2.R \
#       --job-name real_data_ensemble_full \
#       --parent-job-id 53999999 \
#       --output-root output
#
# Outputs (under output/slurm_output/<job>/<id>/):
#   hg2_corrected/hg2_corrected_by_run.csv
#   hg2_corrected/hg2_corrected_by_locus.csv
#
# NOTE: this driver cannot be unit-tested locally (needs the persisted fits +
# LD matrices). Sanity-check the outputs after the first run:
#   - per fit: hg2_expected_pve >= hg2_postmean (= the existing pve_postmean_std);
#   - per locus: the ensemble expected PVE is a convex combination of the
#     positively-weighted per-fit values, so it must lie within
#     [hg2_expected_pve_pos_min, hg2_expected_pve_pos_max]. It can be ABOVE or
#     BELOW the max-ELBO member -- do NOT expect ensemble >= max-ELBO.
#   - n_missing_components == 0 (any positively-weighted fit lacking second
#     moments forces the ensemble aggregate to NA rather than a silent drop).

suppressMessages({
  if (file.exists("../susine/DESCRIPTION")) devtools::load_all("../susine", quiet = TRUE)
  devtools::load_all(".", quiet = TRUE)
})

# --- args -------------------------------------------------------------------
parse_args <- function(args) {
  out <- list(job_name = NULL, parent_job_id = NULL, output_root = "output")
  i <- 1L
  while (i <= length(args)) {
    key <- args[[i]]
    val <- if (i < length(args)) args[[i + 1L]] else NA_character_
    switch(key,
      "--job-name"      = { out$job_name <- val; i <- i + 2L },
      "--parent-job-id" = { out$parent_job_id <- val; i <- i + 2L },
      "--output-root"   = { out$output_root <- val; i <- i + 2L },
      { i <- i + 1L }
    )
  }
  if (is.null(out$job_name) || is.null(out$parent_job_id)) {
    stop("Usage: recompute_real_data_hg2.R --job-name X --parent-job-id Y [--output-root output]")
  }
  out
}
opt <- parse_args(commandArgs(trailingOnly = TRUE))

results_dir <- file.path(opt$output_root, "slurm_output", opt$job_name, opt$parent_job_id)
history_dir <- file.path(opt$output_root, "run_history", opt$job_name, opt$parent_job_id)
temp_dir    <- file.path(opt$output_root, "temp", opt$job_name)
out_dir     <- file.path(results_dir, "hg2_corrected")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

pick_existing <- function(...) {
  for (p in c(...)) if (file.exists(p)) return(p)
  stop("None of the candidate paths exist: ", paste(c(...), collapse = ", "))
}
job_config_path <- pick_existing(file.path(temp_dir, "job_config.json"),
                                 file.path(history_dir, "job_config.json"))
run_manifest_path <- pick_existing(file.path(temp_dir, "run_manifest.csv"),
                                   file.path(history_dir, "run_manifest.csv"))

cfg <- jsonlite::read_json(job_config_path, simplifyVector = TRUE)
jsd_threshold <- as.numeric(cfg$job$jsd_threshold %||% 0.15)
softmax_temperature <- as.numeric(cfg$job$softmax_temperature %||% 1)
manifest_path <- cfg$job$manifest_path
repo_root <- cfg$paths$repo_root %||% getwd()

run_manifest <- readr::read_csv(run_manifest_path, show_col_types = FALSE)

# Per-task run_metrics carry run_id -> fit_rds_path + elbo_final.
rm_files <- Sys.glob(file.path(results_dir, "task-*", "*run_metrics.csv"))
if (!length(rm_files)) {
  rm_files <- list.files(results_dir, pattern = "run_metrics.*\\.csv$",
                         recursive = TRUE, full.names = TRUE)
}
if (!length(rm_files)) stop("No per-task run_metrics CSVs found under ", results_dir)
run_metrics <- dplyr::bind_rows(lapply(rm_files, readr::read_csv, show_col_types = FALSE)) %>%
  dplyr::distinct(.data$run_id, .keep_all = TRUE)

# Join manifest metadata (locus_id, run_family, backend, c_value) onto fit paths.
meta_cols <- intersect(
  c("run_id", "locus_id", "run_family", "backend", "c_value", "sigma_0_2_scalar"),
  names(run_manifest)
)
runs <- run_metrics %>%
  dplyr::select(dplyr::any_of(c("run_id", "fit_rds_path", "elbo_final"))) %>%
  dplyr::left_join(dplyr::select(run_manifest, dplyr::all_of(meta_cols)), by = "run_id") %>%
  dplyr::filter(!is.na(.data$fit_rds_path), nzchar(.data$fit_rds_path))

cat(sprintf("Loaded %d runs across %d loci.\n",
            nrow(runs), dplyr::n_distinct(runs$locus_id)))

# --- per-fit components -----------------------------------------------------
bundle_cache <- new.env(parent = emptyenv())
get_R <- function(locus_id) {
  key <- as.character(locus_id)
  if (exists(key, envir = bundle_cache, inherits = FALSE)) {
    return(get(key, envir = bundle_cache, inherits = FALSE))
  }
  b <- tryCatch(
    test_susine:::load_real_data_locus_bundle(
      locus_id = locus_id, manifest_path = manifest_path, repo_root = repo_root
    ),
    error = function(e) NULL
  )
  R <- if (is.null(b)) NULL else b$R
  assign(key, R, envir = bundle_cache)
  R
}

components_for_fit <- function(fit, R) {
  # Saved susine_rss fits carry unconditional effect_fits$b_hat/b_2_hat; susieR
  # anchors carry conditional mu/mu2 + alpha. Build a uniform fit shell.
  ef <- fit$effect_fits
  if (!is.null(ef$b_hat) && !is.null(ef$b_2_hat)) {
    shell <- list(effect_fits = list(b_hat = ef$b_hat, b_2_hat = ef$b_2_hat))
  } else if (!is.null(fit$mu) && !is.null(fit$mu2) && !is.null(fit$alpha)) {
    a <- as.matrix(fit$alpha)
    shell <- list(effect_fits = list(b_hat = a * as.matrix(fit$mu),
                                      b_2_hat = a * as.matrix(fit$mu2)))
  } else {
    return(list(postmean = NA_real_, uncertainty = NA_real_, expected_pve = NA_real_,
                pips = NULL, elbo = NA_real_))
  }
  comp <- test_susine:::hg2_components(shell, R = R, vy = 1)
  comp$pips <- tryCatch(as.numeric(fit$model_fit$PIPs %||% fit$pip), error = function(e) NULL)
  comp$elbo <- tryCatch(as.numeric(tail(fit$model_fit$elbo %||% fit$elbo, 1)), error = function(e) NA_real_)
  comp
}

per_run_rows <- list()
pip_store <- list()   # locus_id -> list(run_id -> pips)
elbo_store <- list()  # locus_id -> named numeric run_id -> elbo

for (i in seq_len(nrow(runs))) {
  row <- runs[i, , drop = FALSE]
  locus_id <- as.character(row$locus_id[[1]])
  R <- get_R(locus_id)
  if (is.null(R)) next
  fit <- tryCatch(readRDS(row$fit_rds_path[[1]]), error = function(e) NULL)
  if (is.null(fit)) next
  comp <- components_for_fit(fit, R)

  per_run_rows[[length(per_run_rows) + 1L]] <- tibble::tibble(
    run_id = as.integer(row$run_id[[1]]),
    locus_id = locus_id,
    run_family = as.character(row$run_family[[1]] %||% NA_character_),
    c_value = as.numeric(row$c_value[[1]] %||% NA_real_),
    sigma_0_2_scalar = as.numeric(row$sigma_0_2_scalar[[1]] %||% NA_real_),
    hg2_postmean = comp$postmean,
    hg2_uncertainty = comp$uncertainty,
    hg2_expected_pve = comp$expected_pve,
    elbo_final = comp$elbo %||% as.numeric(row$elbo_final[[1]] %||% NA_real_)
  )

  if (identical(as.character(row$run_family[[1]]), "functional_grid") && !is.null(comp$pips)) {
    pip_store[[locus_id]] <- c(pip_store[[locus_id]] %||% list(),
                               stats::setNames(list(comp$pips), as.character(row$run_id[[1]])))
    elbo_store[[locus_id]] <- c(elbo_store[[locus_id]] %||% numeric(),
                                stats::setNames(comp$elbo, as.character(row$run_id[[1]])))
  }
}

per_run <- dplyr::bind_rows(per_run_rows)
readr::write_csv(per_run, file.path(out_dir, "hg2_corrected_by_run.csv"))
cat(sprintf("Wrote %d per-run rows.\n", nrow(per_run)))

# --- per-locus fair ensemble aggregate --------------------------------------
locus_rows <- list()
for (locus_id in names(pip_store)) {
  pips <- pip_store[[locus_id]]
  if (length(pips) < 1L) next
  run_ids <- names(pips)
  expected_by_run <- per_run %>%
    dplyr::filter(.data$locus_id == !!locus_id,
                  as.character(.data$run_id) %in% run_ids,
                  .data$run_family == "functional_grid")
  ord <- match(run_ids, as.character(expected_by_run$run_id))
  expected_vec <- expected_by_run$hg2_expected_pve[ord]
  postmean_vec <- expected_by_run$hg2_postmean[ord]
  unc_vec      <- expected_by_run$hg2_uncertainty[ord]
  elbo_vec <- as.numeric(elbo_store[[locus_id]][run_ids])
  K <- length(run_ids)

  w <- tryCatch({
    if (K > 1L) {
      pip_mat <- do.call(rbind, lapply(pips, as.numeric))
      pip_cache <- test_susine:::prepare_pip_similarity_cache(pip_mat)
      cw <- test_susine:::.cluster_weights_from_hc(
        pip_cache$hc, jsd_threshold, elbo_vec, n_fits = K,
        softmax_temperature = softmax_temperature
      )
      cw$w_full
    } else {
      1
    }
  }, error = function(e) rep(1 / K, K))

  # Weighted aggregate over POSITIVELY-weighted fits. Returns NA (not a silently
  # renormalized partial sum) if any positive-weight fit is missing its
  # component, so a bogus value never reaches the report.
  weighted_or_na <- function(w, x) {
    pos <- w > 0
    if (!any(pos)) return(NA_real_)
    if (any(!is.finite(x[pos]))) return(NA_real_)
    sum(w[pos] * x[pos])
  }
  pos <- w > 0
  n_missing <- sum(pos & !is.finite(expected_vec))
  pos_finite <- pos & is.finite(expected_vec)
  expected_pos_min <- if (any(pos_finite)) min(expected_vec[pos_finite]) else NA_real_
  expected_pos_max <- if (any(pos_finite)) max(expected_vec[pos_finite]) else NA_real_

  locus_rows[[length(locus_rows) + 1L]] <- tibble::tibble(
    locus_id = locus_id,
    n_fits = K,
    n_missing_components = n_missing,
    hg2_expected_pve_ensemble = weighted_or_na(w, expected_vec),
    hg2_within_fit_ensemble = weighted_or_na(w, unc_vec),
    hg2_postmean_weighted = weighted_or_na(w, postmean_vec),
    hg2_expected_pve_max_elbo = expected_vec[which.max(elbo_vec)],
    hg2_expected_pve_pos_min = expected_pos_min,
    hg2_expected_pve_pos_max = expected_pos_max
  )
}

per_locus <- dplyr::bind_rows(locus_rows)
readr::write_csv(per_locus, file.path(out_dir, "hg2_corrected_by_locus.csv"))
cat(sprintf("Wrote %d per-locus ensemble rows to %s\n", nrow(per_locus), out_dir))
if (nrow(per_locus)) {
  n_na <- sum(is.na(per_locus$hg2_expected_pve_ensemble))
  n_miss <- sum(per_locus$n_missing_components > 0, na.rm = TRUE)
  if (n_miss) {
    cat(sprintf("WARNING: %d loci had >=1 positively-weighted fit missing second ",
                n_miss),
        "moments; their ensemble aggregate is NA (not silently dropped).\n", sep = "")
  }
  cat(sprintf("Ensemble aggregate NA for %d/%d loci.\n", n_na, nrow(per_locus)))
}
cat("Done. Sanity-check: per fit hg2_expected_pve >= hg2_postmean; per locus ",
    "ensemble in [hg2_expected_pve_pos_min, hg2_expected_pve_pos_max].\n", sep = "")
