#!/usr/bin/env Rscript
# Downstream recompute of the corrected local-genetic-variance ("heritability")
# estimand for the real-data case studies. No re-fitting. For each locus it
# computes, per fit,
#
#   hg2_postmean     = m' R m / var(y)            (var(y) = 1, standardized RSS)
#   hg2_uncertainty  = tr(R Cov(b|y)) / var(y)    (within-fit correction)
#   hg2_expected_pve = hg2_postmean + hg2_uncertainty   (= local PVE / cis-h2_g)
#
# where, by the SER single-effect structure (Cov(b_l) = diag(b2_l) - m_l m_l'),
#   tr(R Cov(b|y)) = sum_l [ sum_j R_jj * alpha_b_2_hat[l,j] - m_l' R m_l ].
#
# Efficiency / RAM (loci here are 2-10k variants -> R is dense up to ~800 MB):
#   * ONE LOCUS AT A TIME: build R once, then free it before the next locus, so
#     peak memory is a single R plus that locus's moment matrices.
#   * NO FIT RELOADS by default: the per-effect unconditional moments
#     (alpha_b_hat, alpha_b_2_hat) and inclusion probs (alpha) are read straight
#     from the persisted effect_posteriors parquet. Falls back to reloading that
#     locus's fits only if the parquet lacks those columns.
#   * Batched BLAS: all per-(run,effect) quadratic forms are a single R %*% B_all.
#
# Per locus the fair ensemble aggregate is sum_k w_k * hg2_expected_pve_k using
# the SAME cluster weights the pipeline uses (.cluster_weights_from_hc with the
# job's jsd_threshold / softmax_temperature). See
# refs/decisions/heritability_estimand_decision_2026-06-09.md (sec 5.5).
#
# Run on the HPC where output/ and data/ live (allocate enough RAM for one R):
#   salloc --cpus-per-task=4 --mem=32G --time=00:30:00
#   Rscript inst/scripts/recompute_real_data_hg2.R \
#       --job-name real_data_ensemble_geometric_n20 \
#       --parent-job-id 52906940 --output-root output
#
# Outputs (under output/slurm_output/<job>/<id>/hg2_corrected/):
#   hg2_corrected_by_run.csv, hg2_corrected_by_locus.csv
#
# Sanity-check after the run:
#   - per fit: hg2_expected_pve >= hg2_postmean (== the existing pve_postmean_std);
#   - per locus: ensemble lies within [hg2_expected_pve_pos_min, ..._pos_max]
#     (it is a convex combination -> may be ABOVE or BELOW the max-ELBO member);
#   - n_missing_components == 0.

suppressMessages({
  if (file.exists("../susine/DESCRIPTION")) devtools::load_all("../susine", quiet = TRUE)
  devtools::load_all(".", quiet = TRUE)
})

# --- args -------------------------------------------------------------------
parse_args <- function(args) {
  out <- list(job_name = NULL, parent_job_id = NULL, output_root = "output",
              cores = NULL)
  i <- 1L
  while (i <= length(args)) {
    key <- args[[i]]
    val <- if (i < length(args)) args[[i + 1L]] else NA_character_
    switch(key,
      "--job-name"      = { out$job_name <- val; i <- i + 2L },
      "--parent-job-id" = { out$parent_job_id <- val; i <- i + 2L },
      "--output-root"   = { out$output_root <- val; i <- i + 2L },
      "--cores"         = { out$cores <- val; i <- i + 2L },
      { i <- i + 1L }
    )
  }
  if (is.null(out$job_name) || is.null(out$parent_job_id)) {
    stop("Usage: recompute_real_data_hg2.R --job-name X --parent-job-id Y ",
         "[--output-root output] [--cores N]")
  }
  out
}
opt <- parse_args(commandArgs(trailingOnly = TRUE))

# Cores are only used by the (rare) per-locus fit-reload fallback.
resolve_cores <- function(cli) {
  if (!is.null(cli) && nzchar(cli)) return(max(1L, as.integer(cli)))
  slurm <- suppressWarnings(as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", "")))
  if (!is.na(slurm) && slurm > 0L) return(slurm)
  avail <- tryCatch(parallel::detectCores(), error = function(e) 1L)
  max(1L, min(4L, avail - 1L, na.rm = TRUE))
}
n_cores <- resolve_cores(opt$cores)
use_fork <- n_cores > 1L && .Platform$OS.type == "unix"

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

meta_cols <- intersect(
  c("run_id", "locus_id", "run_family", "backend", "c_value", "sigma_0_2_scalar"),
  names(run_manifest)
)
runs <- run_metrics %>%
  dplyr::select(dplyr::any_of(c("run_id", "fit_rds_path", "elbo_final"))) %>%
  dplyr::left_join(dplyr::select(run_manifest, dplyr::all_of(meta_cols)), by = "run_id") %>%
  dplyr::mutate(run_id = as.integer(.data$run_id),
                locus_id = as.character(.data$locus_id))

locus_ids <- sort(unique(runs$locus_id[!is.na(runs$locus_id)]))
cat(sprintf("Loaded %d runs across %d loci.\n", nrow(runs), length(locus_ids)))

# --- effect_posteriors parquet source (preferred, no fit reloads) -----------
ep_files <- Sys.glob(file.path(results_dir, "task-*", "*effect_posteriors.parquet"))
ep_ds <- NULL
ep_cols <- character(0)
need_ep <- c("run_id", "locus_id", "effect_l", "ld_matrix_index",
             "alpha", "alpha_b_hat", "alpha_b_2_hat")
if (length(ep_files)) {
  ep_ds <- tryCatch(arrow::open_dataset(ep_files), error = function(e) NULL)
  if (!is.null(ep_ds)) ep_cols <- names(ep_ds$schema)
}
parquet_ok <- !is.null(ep_ds) && all(need_ep %in% ep_cols)
cat(if (parquet_ok) {
  sprintf("Using effect_posteriors parquet (%d files); no fit reloads.\n", length(ep_files))
} else {
  "effect_posteriors parquet unavailable/incomplete; falling back to fit reloads.\n"
})

# --- shared compute: long (run,effect,variant) moments -> per-run components --
# `long` columns: run_id, effect_l, ld_matrix_index, alpha, alpha_b_hat,
# alpha_b_2_hat (one locus). R is the p x p LD matrix; dR = diag(R).
compute_locus <- function(long, R) {
  dR <- diag(R)
  p <- nrow(R)
  long <- long[!is.na(long$ld_matrix_index), , drop = FALSE]
  vorder <- sort(unique(long$ld_matrix_index))
  if (length(vorder) != p) {
    # align to R's variant order; assume R rows follow ascending ld_matrix_index.
    row_of <- match(long$ld_matrix_index, vorder)
  } else {
    row_of <- match(long$ld_matrix_index, vorder)
  }
  # one column per (run_id, effect_l)
  re_key <- paste(long$run_id, long$effect_l, sep = "::")
  re_levels <- unique(re_key)
  col_of <- match(re_key, re_levels)
  C <- length(re_levels)
  col_run <- as.integer(sub("::.*$", "", re_levels))

  B  <- matrix(0, p, C)
  B2 <- matrix(0, p, C)
  A  <- matrix(0, p, C)
  idx <- cbind(row_of, col_of)
  B[idx]  <- long$alpha_b_hat
  B2[idx] <- long$alpha_b_2_hat
  A[idx]  <- long$alpha

  run_ids <- sort(unique(col_run))
  nR <- length(run_ids)
  # G: C x nR indicator mapping (run,effect) columns -> run.
  G <- matrix(0, C, nR)
  G[cbind(seq_len(C), match(col_run, run_ids))] <- 1

  RB   <- R %*% B                    # p x C  (single big BLAS call)
  quad_col <- colSums(B * RB)        # m_{r,l}' R m_{r,l}
  quad_run <- as.numeric(crossprod(G, quad_col))           # sum_l per run
  M    <- B %*% G                    # p x nR  total posterior mean per run
  RM   <- R %*% M
  postmean_run <- colSums(M * RM)    # m_r' R m_r
  S2   <- B2 %*% G                   # p x nR  total second moment per run
  diag_run <- as.numeric(crossprod(dR, S2))                # sum_j R_jj * s2

  unc_run <- pmax(0, diag_run - quad_run)   # vy = 1; floor at 0
  postmean_run <- pmax(0, postmean_run)
  expected_run <- postmean_run + unc_run

  # PIPs per run for clustering: 1 - prod_l (1 - alpha).
  # A holds alpha at (variant, run-effect); recover per (run, variant).
  pip_list <- lapply(run_ids, function(rid) {
    cols <- which(col_run == rid)
    1 - apply(1 - A[, cols, drop = FALSE], 1L, prod)
  })
  names(pip_list) <- as.character(run_ids)

  list(run_ids = run_ids,
       postmean = postmean_run, uncertainty = unc_run, expected = expected_run,
       pip_list = pip_list)
}

# Fallback: build the `long` frame for a locus by reloading its fits.
long_from_fits <- function(runs_locus) {
  rows <- if (use_fork) {
    parallel::mclapply(seq_len(nrow(runs_locus)), function(i) fit_long_row(runs_locus[i, ]),
                       mc.cores = n_cores)
  } else {
    lapply(seq_len(nrow(runs_locus)), function(i) fit_long_row(runs_locus[i, ]))
  }
  rows <- Filter(function(r) is.data.frame(r), rows)
  if (!length(rows)) return(NULL)
  dplyr::bind_rows(rows)
}
# Reload one saved fit and emit the long (effect, variant) moment frame. Handles
# both susine fits (effect_fits$b_hat/b_2_hat, unconditional) and susieR fits
# (conditional mu/mu2 + alpha). Used by the fit fallback AND the warm-refit step.
fit_to_long <- function(fit_path, run_id) {
  if (is.na(fit_path) || !nzchar(fit_path)) return(NULL)
  fit <- tryCatch(readRDS(fit_path), error = function(e) NULL)
  if (is.null(fit)) return(NULL)
  ef <- fit$effect_fits
  if (!is.null(ef$b_hat) && !is.null(ef$b_2_hat)) {
    ab <- as.matrix(ef$b_hat); ab2 <- as.matrix(ef$b_2_hat); al <- as.matrix(ef$alpha)
  } else if (!is.null(fit$mu) && !is.null(fit$mu2) && !is.null(fit$alpha)) {
    al <- as.matrix(fit$alpha); ab <- al * as.matrix(fit$mu); ab2 <- al * as.matrix(fit$mu2)
  } else {
    return(NULL)
  }
  L <- nrow(ab); p <- ncol(ab)
  tibble::tibble(
    run_id = as.integer(run_id),
    effect_l = rep(seq_len(L), each = p),
    ld_matrix_index = rep(seq_len(p), times = L),
    alpha = as.numeric(t(al)),
    alpha_b_hat = as.numeric(t(ab)),
    alpha_b_2_hat = as.numeric(t(ab2))
  )
}
fit_long_row <- function(row) fit_to_long(row$fit_rds_path[[1]], row$run_id[[1]])

# Warm-refit fits live under aggregated/highest_weight_refits/ (one per locus) and
# are NOT in the per-task effect_posteriors parquet, so they are reloaded directly.
refit_summary_path <- file.path(results_dir, "aggregated", "highest_weight_refit_summary.csv")
refit_summary <- if (file.exists(refit_summary_path)) {
  readr::read_csv(refit_summary_path, show_col_types = FALSE)
} else {
  tibble::tibble(locus_id = character(), refit_run_id = integer(),
                 refit_fit_rds_path = character(), refit_elbo_final = numeric())
}
cat(sprintf("Warm-refit fits available for %d loci.\n",
            sum(!is.na(refit_summary$refit_fit_rds_path) &
                  nzchar(as.character(refit_summary$refit_fit_rds_path)))))

# --- per-locus driver -------------------------------------------------------
per_run_rows <- vector("list", length(locus_ids))
locus_rows   <- vector("list", length(locus_ids))
elbo_lookup  <- stats::setNames(as.numeric(run_metrics$elbo_final %||% NA_real_),
                                as.character(run_metrics$run_id))

for (li in seq_along(locus_ids)) {
  lid <- locus_ids[[li]]
  runs_locus <- runs[runs$locus_id == lid & !is.na(runs$locus_id), , drop = FALSE]
  cat(sprintf("[%d/%d] %s: %d runs ... ", li, length(locus_ids), lid, nrow(runs_locus)))

  bundle <- tryCatch(
    test_susine:::load_real_data_locus_bundle(
      locus_id = lid, manifest_path = manifest_path, repo_root = repo_root),
    error = function(e) NULL)
  if (is.null(bundle) || is.null(bundle$R)) { cat("no LD matrix; skipped.\n"); next }
  R <- bundle$R
  rm(bundle)

  long <- if (parquet_ok) {
    tryCatch(
      ep_ds %>%
        dplyr::filter(.data$locus_id == lid) %>%
        dplyr::select(dplyr::all_of(setdiff(need_ep, "locus_id"))) %>%
        dplyr::collect(),
      error = function(e) NULL)
  } else NULL
  if (is.null(long) || !nrow(long)) long <- long_from_fits(runs_locus)
  if (is.null(long) || !nrow(long)) { cat("no moments; skipped.\n"); rm(R); gc(); next }

  comp <- compute_locus(long, R)
  rm(long)

  fam <- runs_locus$run_family[match(comp$run_ids, runs_locus$run_id)]
  cval <- runs_locus$c_value[match(comp$run_ids, runs_locus$run_id)]
  sval <- runs_locus$sigma_0_2_scalar[match(comp$run_ids, runs_locus$run_id)]
  elbo <- elbo_lookup[as.character(comp$run_ids)]

  main_row <- tibble::tibble(
    run_id = comp$run_ids, locus_id = lid,
    run_family = as.character(fam), c_value = as.numeric(cval),
    sigma_0_2_scalar = as.numeric(sval),
    hg2_postmean = comp$postmean, hg2_uncertainty = comp$uncertainty,
    hg2_expected_pve = comp$expected, elbo_final = as.numeric(elbo))

  # Warm-refit fit (reloaded; uses the same R before it is freed).
  refit_row <- NULL
  rr <- refit_summary[as.character(refit_summary$locus_id) == lid, , drop = FALSE]
  if (nrow(rr) == 1L && !is.na(rr$refit_fit_rds_path[[1]]) &&
      nzchar(as.character(rr$refit_fit_rds_path[[1]]))) {
    rlong <- fit_to_long(as.character(rr$refit_fit_rds_path[[1]]),
                         as.integer(rr$refit_run_id[[1]]))
    if (!is.null(rlong) && nrow(rlong)) {
      rcomp <- tryCatch(compute_locus(rlong, R), error = function(e) NULL)
      if (!is.null(rcomp)) {
        refit_row <- tibble::tibble(
          run_id = as.integer(rr$refit_run_id[[1]]), locus_id = lid,
          run_family = "highest_weight_refit", c_value = NA_real_,
          sigma_0_2_scalar = NA_real_,
          hg2_postmean = rcomp$postmean[[1]], hg2_uncertainty = rcomp$uncertainty[[1]],
          hg2_expected_pve = rcomp$expected[[1]],
          elbo_final = as.numeric(rr$refit_elbo_final[[1]] %||% NA_real_))
      }
    }
    rm(rlong)
  }
  rm(R); gc(verbose = FALSE)

  per_run_rows[[li]] <- dplyr::bind_rows(main_row, refit_row)

  # Fair ensemble over functional_grid members.
  is_fun <- !is.na(fam) & fam == "functional_grid"
  fun_ids <- comp$run_ids[is_fun]
  if (length(fun_ids) >= 1L) {
    sel <- match(fun_ids, comp$run_ids)
    expected_vec <- comp$expected[sel]
    postmean_vec <- comp$postmean[sel]
    unc_vec      <- comp$uncertainty[sel]
    elbo_vec <- as.numeric(elbo_lookup[as.character(fun_ids)])
    K <- length(fun_ids)
    w <- tryCatch({
      if (K > 1L) {
        pip_mat <- do.call(rbind, comp$pip_list[as.character(fun_ids)])
        pc <- test_susine:::prepare_pip_similarity_cache(pip_mat)
        test_susine:::.cluster_weights_from_hc(
          pc$hc, jsd_threshold, elbo_vec, n_fits = K,
          softmax_temperature = softmax_temperature)$w_full
      } else 1
    }, error = function(e) rep(1 / K, K))

    weighted_or_na <- function(w, x) {
      pos <- w > 0
      if (!any(pos)) return(NA_real_)
      if (any(!is.finite(x[pos]))) return(NA_real_)
      sum(w[pos] * x[pos])
    }
    pos <- w > 0
    posf <- pos & is.finite(expected_vec)
    locus_rows[[li]] <- tibble::tibble(
      locus_id = lid, n_fits = K,
      n_missing_components = sum(pos & !is.finite(expected_vec)),
      hg2_expected_pve_ensemble = weighted_or_na(w, expected_vec),
      hg2_within_fit_ensemble = weighted_or_na(w, unc_vec),
      hg2_postmean_weighted = weighted_or_na(w, postmean_vec),
      hg2_expected_pve_max_elbo = expected_vec[which.max(elbo_vec)],
      hg2_expected_pve_pos_min = if (any(posf)) min(expected_vec[posf]) else NA_real_,
      hg2_expected_pve_pos_max = if (any(posf)) max(expected_vec[posf]) else NA_real_)
  }
  cat("ok.\n")
}

per_run <- dplyr::bind_rows(per_run_rows)
per_locus <- dplyr::bind_rows(locus_rows)
readr::write_csv(per_run, file.path(out_dir, "hg2_corrected_by_run.csv"))
readr::write_csv(per_locus, file.path(out_dir, "hg2_corrected_by_locus.csv"))
cat(sprintf("Wrote %d per-run rows and %d per-locus rows to %s\n",
            nrow(per_run), nrow(per_locus), out_dir))
if (nrow(per_locus)) {
  n_miss <- sum(per_locus$n_missing_components > 0, na.rm = TRUE)
  if (n_miss) cat(sprintf("WARNING: %d loci had >=1 positively-weighted fit missing moments (ensemble NA).\n", n_miss))
}
cat("Done. Sanity-check: per fit hg2_expected_pve >= hg2_postmean; per locus ",
    "ensemble in [hg2_expected_pve_pos_min, hg2_expected_pve_pos_max].\n", sep = "")
