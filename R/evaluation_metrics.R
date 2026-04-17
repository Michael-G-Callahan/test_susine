# evaluation_metrics.R
# Comprehensive evaluation for SuSiE/SuSiNE-like model fits.
# Produces BOTH purity-FILTERED (primary, as in SuSiE) and UNFILTERED metrics.
#
# Inputs:
#   fit : model object with
#         - fit$effect_fits$alpha  (L x p matrix of per-effect PIPs α_lj)
#         - fit$settings$L         (integer)
#         - optional: fit$model_fit$fitted_y, fit$model_fit$sigma_2, fit$model_fit$coef
#   X   : n x p matrix (columns align with α / variables)
#   y   : length-n response vector (only needed for h_g^2)
#   causal_idx : integer indices of true causal variables (for simulation-based metrics)
#
# Key definitions (SuSiE-style):
#   - ρ-level CS: sort α decreasing, take smallest prefix with cumulative sum ≥ ρ.
#   - Purity: min absolute correlation among all pairs inside a CS.
#   - Coverage (per-effect): 1 if CS contains ≥ 1 true causal, else 0.
#   - Power (per-model): fraction of true causals captured by ≥ 1 CS.
#   - Combined PIP per variable: 1 - prod_l (1 - α_lj).
#
# References:
#   - Credible sets & purity filtering: SuSiE §3.2.2. :contentReference[oaicite:1]{index=1}
#   - Power/FDR vs PIP thresholds (PR-curve relation). :contentReference[oaicite:2]{index=2}

# ---------- helpers ----------

# ρ-level credible set by cumulative PIP (SuSiE definition)
get_credible_set <- function(alpha, rho = 0.95) {
  if (!is.numeric(alpha)) stop("alpha must be numeric PIPs")
  if (!is.finite(rho) || rho <= 0 || rho >= 1) stop("rho must be in (0,1)")

  a <- alpha
  a[!is.finite(a)] <- 0
  if (sum(a) <= 0) return(integer(0))  # nothing to select

  # order by decreasing PIP; break ties by original index for determinism
  o  <- order(-a, seq_along(a))
  cs <- cumsum(a[o])
  k  <- which(cs >= rho)[1]
  if (is.na(k)) integer(0) else o[seq_len(k)]
}

# Purity = min absolute pairwise correlation inside a CS (size-1 -> 1, empty -> NA)
cs_purity_min_abs <- function(X, idx, cache = NULL, cor_mat = NULL) {
  k <- length(idx)
  if (k == 0) return(NA_real_)
  if (k == 1) return(1.0)
  key <- NULL
  if (!is.null(cache)) {
    key <- paste(idx, collapse = ",")
    if (exists(key, envir = cache, inherits = FALSE)) {
      return(get(key, envir = cache, inherits = FALSE))
    }
  }
  if (!is.null(cor_mat)) {
    sub_cor <- cor_mat[idx, idx, drop = FALSE]
  } else {
    sub_cor <- suppressWarnings(cor(X[, idx, drop = FALSE], use = "pairwise.complete.obs"))
  }
  purity <- min(abs(sub_cor[upper.tri(sub_cor)]))
  if (!is.null(cache) && length(key)) {
    assign(key, purity, envir = cache)
  }
  purity
}

# Overlap rate among a list of index vectors (fraction of CS pairs that overlap)
overlap_rate_from_sets <- function(list_of_sets) {
  m <- length(list_of_sets)
  if (m < 2) return(NA_real_)
  pairs <- utils::combn(m, 2)
  overlaps <- apply(pairs, 2, function(col) {
    i <- list_of_sets[[col[1]]]
    j <- list_of_sets[[col[2]]]
    length(intersect(i, j)) > 0
  })
  mean(overlaps)
}

# Normalize a numeric vector onto the probability simplex.
normalize_prob_vec <- function(x) {
  if (!is.numeric(x)) stop("x must be numeric")
  x <- as.numeric(x)
  x[!is.finite(x)] <- 0
  x[x < 0] <- 0
  s <- sum(x)
  if (!is.finite(s) || s <= 0) {
    return(rep(0, length(x)))
  }
  x / s
}

prob_entropy <- function(prob) {
  p <- normalize_prob_vec(prob)
  if (!any(p > 0)) return(NA_real_)
  -sum(p[p > 0] * log(p[p > 0]))
}

prob_k_eff <- function(prob) {
  h <- prob_entropy(prob)
  if (!is.finite(h)) return(NA_real_)
  exp(h)
}

prob_entropy_core <- function(prob, rho = 0.95) {
  p <- normalize_prob_vec(prob)
  idx <- get_credible_set(p, rho = rho)
  if (!length(idx)) return(NA_real_)
  prob_entropy(p[idx])
}

prob_k_eff_core <- function(prob, rho = 0.95) {
  h <- prob_entropy_core(prob, rho = rho)
  if (!is.finite(h)) return(NA_real_)
  exp(h)
}

effect_accuracy_ratio <- function(alpha_vec, causal_idx) {
  p <- normalize_prob_vec(alpha_vec)
  if (!length(causal_idx)) return(NA_real_)

  valid_causal <- unique(as.integer(causal_idx[is.finite(causal_idx)]))
  valid_causal <- valid_causal[valid_causal >= 1L & valid_causal <= length(p)]
  if (!length(valid_causal)) return(NA_real_)

  max_any <- max(p, na.rm = TRUE)
  if (!is.finite(max_any) || max_any <= 0) return(NA_real_)

  max_causal <- max(p[valid_causal], na.rm = TRUE)
  if (!is.finite(max_causal)) return(NA_real_)

  max_causal / max_any
}

# Per-effect metrics given one alpha vector
cs_metrics_one_effect <- function(alpha_vec, X, causal_idx, rho = 0.95,
                                  purity_cache = NULL, cor_mat = NULL) {
  alpha_prob <- normalize_prob_vec(alpha_vec)
  cs <- get_credible_set(alpha_prob, rho)
  size <- length(cs)
  purity <- cs_purity_min_abs(X, cs, cache = purity_cache, cor_mat = cor_mat)
  coverage <- as.integer(size > 0 && any(cs %in% causal_idx))
  effect_pip_entropy <- prob_entropy(alpha_prob)
  effect_pip_entropy_core95 <- prob_entropy_core(alpha_prob, rho = 0.95)
  effect_k_eff_signal <- prob_k_eff(alpha_prob)
  effect_k_eff_signal_core95 <- prob_k_eff_core(alpha_prob, rho = 0.95)
  tail_inflation_ratio <- if (is.finite(effect_k_eff_signal) &&
                              is.finite(effect_k_eff_signal_core95) &&
                              effect_k_eff_signal_core95 > 0) {
    effect_k_eff_signal / effect_k_eff_signal_core95
  } else {
    NA_real_
  }
  tail_inflation_log <- if (is.finite(effect_pip_entropy) &&
                            is.finite(effect_pip_entropy_core95)) {
    effect_pip_entropy - effect_pip_entropy_core95
  } else {
    NA_real_
  }
  accuracy_ratio <- effect_accuracy_ratio(alpha_prob, causal_idx)
  list(
    indices = cs,
    size = size,
    purity = purity,
    coverage = coverage,
    effect_pip_entropy = effect_pip_entropy,
    effect_pip_entropy_core95 = effect_pip_entropy_core95,
    effect_k_eff_signal = effect_k_eff_signal,
    effect_k_eff_signal_core95 = effect_k_eff_signal_core95,
    tail_inflation_ratio = tail_inflation_ratio,
    tail_inflation_log = tail_inflation_log,
    accuracy_ratio = accuracy_ratio
  )
}

# ---------- classification-style metrics on combined PIPs ----------

# Average Precision (AUPRC) for binary labels using scores in [0,1]
# Implemented as mean precision at each positive's rank.
auprc_average_precision <- function(scores, labels) {
  if (!is.numeric(scores) || !is.numeric(labels))
    stop("scores and labels must be numeric")
  if (length(scores) != length(labels))
    stop("lengths differ")
  labels <- as.integer(labels > 0)
  P <- sum(labels)
  if (P == 0L) return(NA_real_)            # undefined if no positives
  if (P == length(labels)) return(1.0)     # all positive -> AP = 1

  o <- order(scores, decreasing = TRUE)
  y <- labels[o]
  tp <- cumsum(y)
  fp <- cumsum(1L - y)
  pos_idx <- which(y == 1L)
  mean(tp[pos_idx] / (tp[pos_idx] + fp[pos_idx]))
}

# TPR at specified FPR thresholds — derived from a sorted PIP score list.
# Returns a named list with one entry per threshold, suitable for bind_cols.
compute_tpr_at_fpr <- function(pip, causal, fpr_thresholds = c(0.05, 0.1, 0.2, 0.5)) {
  causal      <- as.integer(causal > 0)
  n_causal    <- sum(causal)
  n_noncausal <- length(causal) - n_causal
  col_names   <- paste0("tpr_at_fpr_",
                        formatC(fpr_thresholds * 100L, format = "d", flag = "0", width = 2))
  if (n_causal == 0L || n_noncausal == 0L) {
    return(purrr::set_names(as.list(rep(NA_real_, length(fpr_thresholds))), col_names))
  }
  ord      <- order(pip, decreasing = TRUE)
  causal_o <- causal[ord]
  fpr_vec  <- cumsum(1L - causal_o) / n_noncausal
  tpr_vec  <- cumsum(causal_o) / n_causal
  vals <- purrr::map_dbl(fpr_thresholds, function(thr) {
    idx <- which(fpr_vec >= thr)[1L]
    if (is.na(idx)) NA_real_ else tpr_vec[idx]
  })
  purrr::set_names(as.list(vals), col_names)
}

# CS power and FDR considering only the top-n causal variants by |beta|.
# Used for tier-stratified evaluation of susie2_oligogenic architectures.
# effects_filtered: output of evaluate_model()$effects_filtered (has $indices list-column).
compute_cs_power_by_tier <- function(effects_filtered, causal_idx, beta,
                                     n_top_vec = c(3L, 5L, 7L, 9L, 11L, 23L)) {
  if (length(causal_idx) == 0L) {
    return(tibble::tibble(
      n_top           = as.integer(n_top_vec),
      cs_power        = NA_real_,
      cs_fdr          = NA_real_,
      n_causal_covered = NA_integer_,
      n_causal_total  = NA_integer_,
      n_false_cs      = NA_integer_,
      n_cs_total      = NA_integer_
    ))
  }
  causal_order <- causal_idx[order(abs(beta[causal_idx]), decreasing = TRUE)]
  cs_indices   <- effects_filtered$indices   # list-column
  union_cs     <- sort(unique(unlist(cs_indices)))
  n_cs_total   <- length(cs_indices)
  purrr::map_dfr(n_top_vec, function(n_top) {
    top_causal       <- head(causal_order, n_top)
    n_causal_total   <- length(top_causal)
    n_causal_covered <- if (length(union_cs) > 0L)
      length(intersect(union_cs, top_causal))
    else 0L
    cs_power <- n_causal_covered / max(n_causal_total, 1L)
    if (n_cs_total > 0L) {
      has_causal <- vapply(cs_indices, function(idx) {
        length(intersect(idx, top_causal)) > 0L
      }, logical(1L))
      n_false_cs <- sum(!has_causal)
      cs_fdr     <- n_false_cs / n_cs_total
    } else {
      n_false_cs <- NA_integer_
      cs_fdr     <- NA_real_
    }
    tibble::tibble(
      n_top            = as.integer(n_top),
      cs_power         = cs_power,
      cs_fdr           = cs_fdr,
      n_causal_covered = as.integer(n_causal_covered),
      n_causal_total   = as.integer(n_causal_total),
      n_false_cs       = as.integer(n_false_cs),
      n_cs_total       = as.integer(n_cs_total)
    )
  })
}

# Cross-entropy (log loss) for binary labels and predicted probs
cross_entropy_loss <- function(scores, labels, eps = 1e-12) {
  if (length(scores) != length(labels))
    stop("lengths differ")
  p <- pmin(pmax(scores, eps), 1 - eps)
  y <- as.integer(labels > 0)
  -mean(y * log(p) + (1 - y) * log(1 - p))
}

# SNP-heritability / variance explained estimate
# Prefers fitted_y; else uses residual variance; else uses coef.
estimate_hg2 <- function(fit, y, X) {
  if (missing(y) || is.null(y)) return(NA_real_)
  vy <- stats::var(y)
  if (!is.finite(vy) || vy <= 0) return(NA_real_)

  # 1) fitted_y available?
  fy <- tryCatch(fit$model_fit$fitted_y, error = function(e) NULL)
  if (!is.null(fy)) {
    h2 <- stats::var(as.numeric(fy)) / vy
    return(max(0, min(1, h2)))
  }

  # 2) residual variance path available?
  s2 <- tryCatch(fit$model_fit$sigma_2, error = function(e) NULL)
  if (!is.null(s2) && length(s2)) {
    s2_last <- tail(s2[is.finite(s2)], 1)
    if (length(s2_last)) {
      h2 <- 1 - s2_last / vy
      return(max(0, min(1, h2)))
    }
  }

  # 3) posterior mean coefficients available?
  b <- tryCatch(fit$model_fit$coef, error = function(e) NULL)
  if (!is.null(b)) {
    mu <- as.numeric(X %*% b)
    h2 <- stats::var(mu) / vy
    return(max(0, min(1, h2)))
  }

  NA_real_
}

# ---------- main evaluator ----------

evaluate_model <- function(fit, X, y = NULL, causal_idx = integer(0),
                           rho = 0.95, purity_threshold = 0.5,
                           compute_curves = TRUE) {
  stopifnot(is.matrix(X))
  L <- fit$settings$L
  A <- fit$effect_fits$alpha
  if (is.list(A)) {
    A <- do.call(rbind, A)
  }
  A <- as.matrix(A)
  if (!nrow(A) || !ncol(A)) {
    stop("effect_fits$alpha has invalid dimensions.")
  }
  if (nrow(A) != L) {
    L <- nrow(A)
  }

  purity_cache <- new.env(parent = emptyenv())
  cor_mat <- suppressWarnings(cor(X, use = "pairwise.complete.obs"))

  # Per-effect CS & metrics (UNFILTERED)
  eff <- vector("list", L)
  for (l in seq_len(L)) {
    eff[[l]] <- cs_metrics_one_effect(
      alpha_vec = A[l, ],
      X = X,
      causal_idx = causal_idx,
      rho = rho,
      purity_cache = purity_cache,
      cor_mat = cor_mat
    )
  }

  effects_unfiltered <- data.frame(
    effect                      = seq_len(L),
    size                        = sapply(eff, `[[`, "size"),
    purity                      = sapply(eff, `[[`, "purity"),
    coverage                    = sapply(eff, `[[`, "coverage"),
    effect_pip_entropy          = sapply(eff, `[[`, "effect_pip_entropy"),
    effect_pip_entropy_core95   = sapply(eff, `[[`, "effect_pip_entropy_core95"),
    effect_k_eff_signal         = sapply(eff, `[[`, "effect_k_eff_signal"),
    effect_k_eff_signal_core95  = sapply(eff, `[[`, "effect_k_eff_signal_core95"),
    tail_inflation_ratio        = sapply(eff, `[[`, "tail_inflation_ratio"),
    tail_inflation_log          = sapply(eff, `[[`, "tail_inflation_log"),
    accuracy_ratio              = sapply(eff, `[[`, "accuracy_ratio"),
    stringsAsFactors = FALSE
  )
  effects_unfiltered$indices <- I(lapply(eff, `[[`, "indices"))

  # FILTERED view (keep CS with purity >= threshold)
  keep_idx <- which(!is.na(effects_unfiltered$purity) &
                      effects_unfiltered$purity >= purity_threshold)
  effects_filtered <- effects_unfiltered[keep_idx, , drop = FALSE]

  # Unions for power
  union_unfiltered <- sort(unique(unlist(effects_unfiltered$indices)))
  union_filtered   <- sort(unique(unlist(effects_filtered$indices)))

  # Power (fraction of true causals hit at least once)
  power_unfiltered <- if (length(causal_idx) > 0)
    length(intersect(union_unfiltered, causal_idx)) / length(causal_idx) else NA_real_
  power_filtered <- if (length(causal_idx) > 0)
    length(intersect(union_filtered, causal_idx)) / length(causal_idx) else NA_real_

  # Overlap rate
  overlap_unfiltered <- overlap_rate_from_sets(effects_unfiltered$indices)
  overlap_filtered   <- overlap_rate_from_sets(effects_filtered$indices)

  # Effective #effects
  L_eff_unfiltered <- sum(effects_unfiltered$size > 0, na.rm = TRUE)
  L_eff_filtered   <- nrow(effects_filtered)

  # Variable-level PIPs (combined across effects)
  combined_pip <- fit$model_fit$PIPs

  # -------- classification-style metrics vs ground truth (if provided) --------
  labels <- rep(0L, ncol(X))
  if (length(causal_idx) > 0) labels[causal_idx] <- 1L

  auprc     <- if (sum(labels) > 0) auprc_average_precision(combined_pip, labels) else NA_real_
  xent      <- if (length(labels) == length(combined_pip))
                 cross_entropy_loss(combined_pip, labels) else NA_real_
  tpr_vals  <- compute_tpr_at_fpr(combined_pip, labels)

  # -------- h_g^2 --------
  hg2 <- estimate_hg2(fit, y, X)

  # Optional traces if present
  elbo <- tryCatch(fit$model_fit$elbo, error = function(e) NULL)
  sigma2_path <- tryCatch(fit$model_fit$sigma_2, error = function(e) NULL)

  # Model-level summaries
  model_unfiltered <- dplyr::bind_cols(
    data.frame(
      L_nominal      = L,
      L_effective    = L_eff_unfiltered,
      mean_size      = mean(effects_unfiltered$size, na.rm = TRUE),
      mean_purity    = mean(effects_unfiltered$purity, na.rm = TRUE),
      mean_coverage  = mean(effects_unfiltered$coverage, na.rm = TRUE),
      power          = power_unfiltered,
      overlap_rate   = overlap_unfiltered,
      AUPRC          = auprc,
      cross_entropy  = xent,
      hg2            = hg2,
      stringsAsFactors = FALSE
    ),
    tibble::as_tibble(tpr_vals)
  )

  model_filtered <- dplyr::bind_cols(
    data.frame(
      L_nominal      = L,
      L_effective    = L_eff_filtered,
      mean_size      = mean(effects_filtered$size, na.rm = TRUE),
      mean_purity    = mean(effects_filtered$purity, na.rm = TRUE),
      mean_coverage  = mean(effects_filtered$coverage, na.rm = TRUE),
      power          = power_filtered,
      overlap_rate   = overlap_filtered,
      AUPRC          = auprc,     # classification metrics don't change with filtering
      cross_entropy  = xent,
      hg2            = hg2,
      stringsAsFactors = FALSE
    ),
    tibble::as_tibble(tpr_vals)
  )

  out <- list(
    params = list(rho = rho, purity_threshold = purity_threshold),
    # per-effect
    effects_unfiltered = effects_unfiltered,
    effects_filtered   = effects_filtered,
    # unions
    cs_union_indices_unfiltered = union_unfiltered,
    cs_union_indices_filtered   = union_filtered,
    # model summaries
    model_unfiltered = model_unfiltered,
    model_filtered   = model_filtered,
    # extras
    combined_pip = combined_pip,
    traces = list(elbo = elbo, sigma2 = sigma2_path)
  )

  if (!compute_curves) return(out)

  # -------- optional curves (useful for regression tests & calibration) --------
  # ρ-sweep coverage curve (empirical coverage vs rho)
  rho_grid <- seq(0.75, 0.99, by = 0.02)
  coverage_curve <- if (length(causal_idx) > 0) {
    sapply(rho_grid, function(rh) {
      cov_vec <- sapply(seq_len(L), function(l) {
        cs_l <- get_credible_set(A[l, ], rh)
        as.integer(length(cs_l) > 0 && any(cs_l %in% causal_idx))
      })
      mean(cov_vec, na.rm = TRUE)
    })
  } else rep(NA_real_, length(rho_grid))

  # PIP calibration (bin by combined PIP)
  if (length(causal_idx) > 0) {
    brks <- seq(0, 1, length.out = 21)  # 20 bins
    bins <- cut(combined_pip, brks, include.lowest = TRUE)
    pip_calibration <- aggregate(
      data.frame(pip = combined_pip,
                 is_causal = seq_along(combined_pip) %in% causal_idx),
      by = list(bin = bins),
      FUN = mean
    )
  } else {
    pip_calibration <- NULL
  }

  out$curves <- list(
    rho = rho_grid,
    coverage = coverage_curve,
    pip_calibration = pip_calibration
  )
  out
}

# ---------- convenience printer ----------

print_model_summary <- function(res) {
  cat("Credible-set evaluation (rho =", res$params$rho,
      ", purity_threshold =", res$params$purity_threshold, ")\n\n")

  cat("== Unfiltered ==\n")
  print(res$model_unfiltered, row.names = FALSE)

  cat("\n== Purity-filtered ==\n")
  print(res$model_filtered, row.names = FALSE)

  cat("\n#CS kept (filtered):", nrow(res$effects_filtered),
      " | #vars in CS union (filtered):", length(res$cs_union_indices_filtered), "\n")
}

# ---------- scaling analysis helpers ----------

# Compute confusion bins from a pre-aggregated PIP vector.
# Matches the schema of flush confusion_bins files:
#   pip_threshold, n_causal_at_bucket, n_noncausal_at_bucket
# Used by compute_scaling_curves() and compute_interaction_scaling() in collect_results.R.
#' @keywords internal
compute_bins_from_pip_vec <- function(pip_vec, causal_vec, pip_breaks) {
  bucket_idx <- findInterval(pip_vec, sort(pip_breaks), rightmost.closed = TRUE)
  tibble::tibble(
    bucket = bucket_idx,
    pip    = pip_vec,
    causal = as.integer(causal_vec)
  ) %>%
    dplyr::group_by(.data$bucket) %>%
    dplyr::summarise(
      pip_threshold         = min(.data$pip),
      n_causal_at_bucket    = sum(.data$causal),
      n_noncausal_at_bucket = dplyr::n() - sum(.data$causal),
      .groups = "drop"
    ) %>%
    dplyr::arrange(dplyr::desc(.data$pip_threshold))
}

# TPR at a fixed FPR threshold from a pooled confusion-bins tibble.
# pooled_bins must have columns: pip_threshold, n_causal_at_bucket, n_noncausal_at_bucket.
# Mirrors the compute_tpr05_from_confusion() logic in the collect workbook.
#' @keywords internal
tpr_at_fpr_threshold <- function(pooled_bins, fpr_threshold = 0.05) {
  b       <- pooled_bins %>% dplyr::arrange(dplyr::desc(.data$pip_threshold))
  cum_tp  <- cumsum(b$n_causal_at_bucket)
  cum_fp  <- cumsum(b$n_noncausal_at_bucket)
  total_p <- sum(b$n_causal_at_bucket)
  total_n <- sum(b$n_noncausal_at_bucket)
  if (total_p == 0L || total_n == 0L) return(NA_real_)
  tpr <- cum_tp / total_p
  fpr <- cum_fp / total_n
  idx <- which(fpr <= fpr_threshold)
  if (length(idx) == 0L) return(NA_real_)
  max(tpr[idx])
}

# ---------- example (comment or adapt) ----------
# res <- evaluate_model(fit = fit_a_i, X = X, y = y, causal_idx = causal_idx,
#                       rho = 0.95, purity_threshold = 0.5, compute_curves = TRUE)
# print_model_summary(res)
# head(res$combined_pip)
