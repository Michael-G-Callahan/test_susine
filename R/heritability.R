# Local genetic-variance ("heritability") estimands ---------------------------
#
# Implements the fair estimand from
# refs/decisions/heritability_estimand_decision_2026-06-09.md:
#
#   E[var(Xb) | y] / var(y)  =  m' Sigma m / var(y)            (posterior-mean part)
#                            +  tr(Sigma Cov(b|y)) / var(y)    (within-fit correction)
#
# where m = colSums(b_hat) is the posterior-mean coefficient vector and, under
# the SuSiE/SuSiNE mean-field SER posterior, the per-effect covariance is
#   Cov(b_l) = diag(b_2_hat[l, ]) - m_l m_l',   m_l = b_hat[l, ]
# (off-diagonals of E[b_l b_l'] vanish because each single effect has exactly one
# active variable). Summing over the L independent effects gives
#   tr(Sigma Cov(b|y)) = sum_l [ sum_j Sigma_jj b_2_hat[l, j] - m_l' Sigma m_l ].
#
# `b_hat` / `b_2_hat` are the UNCONDITIONAL per-effect moments stored by susine
# (`SER.R`: b_hat = alpha * mu_1, b_2_hat = alpha * (sigma_1_2 + mu_1^2)). For a
# susieR fit, build them as alpha * mu and alpha * mu2.
#
# Denominator convention: the plotted metric is stats::var(fitted_y)/stats::var(y)
# (the (n-1) cancels). The Gram matrix Sigma must reproduce var(Xb) under that
# same convention, i.e. centered crossprod(X)/(n-1) for individual data (achieved
# here by evaluating quadratic forms via stats::var(X %*% v)), or the LD matrix R
# with var(y) = 1 for standardized RSS.

#' Local-genetic-variance components from per-effect posterior moments.
#'
#' Low-level engine: given the unconditional moment matrices and a quadratic-form
#' closure for the Gram matrix, returns the posterior-mean part, the within-fit
#' correction, and their sum, each as a fraction of `vy`.
#'
#' Components are floored at 0 (the correction is non-negative in exact
#' arithmetic). `expected_pve` is NOT capped at 1 so that the additive
#' decomposition (postmean + uncertainty [+ between-fit]) holds exactly; clamp at
#' the display layer if desired.
#'
#' @param b_hat L x p matrix of unconditional per-effect posterior means.
#' @param b_2_hat L x p matrix of unconditional per-effect posterior second
#'   moments.
#' @param quad_form function(v) returning `v' Sigma v` for a length-p vector `v`.
#' @param gram_diag length-p vector `diag(Sigma)`.
#' @param vy scalar `var(y)` under the matching convention.
#' @return list(postmean, uncertainty, expected_pve).
#' @keywords internal
hg2_components_from_moments <- function(b_hat, b_2_hat, quad_form, gram_diag, vy) {
  na_out <- list(postmean = NA_real_, uncertainty = NA_real_,
                 expected_pve = NA_real_)
  if (is.null(b_hat) || is.null(b_2_hat)) {
    return(na_out)
  }
  b_hat <- as.matrix(b_hat)
  b_2_hat <- as.matrix(b_2_hat)
  if (!nrow(b_hat) || !ncol(b_hat) ||
      !identical(dim(b_hat), dim(b_2_hat)) ||
      length(gram_diag) != ncol(b_hat)) {
    return(na_out)
  }
  if (!is.finite(vy) || vy <= 0) {
    return(na_out)
  }

  m <- colSums(b_hat)
  postmean <- quad_form(m) / vy

  corr <- 0
  for (l in seq_len(nrow(b_hat))) {
    diag_term <- sum(gram_diag * b_2_hat[l, ])
    quad_term <- quad_form(b_hat[l, ])
    corr <- corr + (diag_term - quad_term)
  }
  uncertainty <- corr / vy

  if (!is.finite(postmean) || !is.finite(uncertainty)) {
    return(na_out)
  }
  postmean <- max(0, postmean)
  uncertainty <- max(0, uncertainty)
  list(
    postmean = postmean,
    uncertainty = uncertainty,
    expected_pve = postmean + uncertainty
  )
}

#' Local-genetic-variance components for a single susine/susieR fit.
#'
#' Dispatches between individual-level data (supply `X`, `y`) and standardized
#' RSS (supply `R`; `var(y)` defaults to 1). Reads the unconditional moment
#' matrices from `fit$effect_fits$b_hat` / `$b_2_hat`; returns all-NA if absent.
#'
#' @param fit A finalized susine/susieR-normalized fit.
#' @param X n x p design matrix (individual-level). Provide `X` or `R`.
#' @param y Response vector (length n), individual-level path.
#' @param R p x p LD matrix (standardized RSS path). Provide `X` or `R`.
#' @param vy Optional override for `var(y)` (defaults: `stats::var(y)` for the
#'   individual path, `1` for RSS).
#' @return list(postmean, uncertainty, expected_pve), fractions of `var(y)`.
#' @keywords internal
hg2_components <- function(fit, X = NULL, y = NULL, R = NULL, vy = NULL) {
  ef <- fit$effect_fits
  b_hat <- ef$b_hat
  b_2_hat <- ef$b_2_hat
  if (is.null(b_hat) || is.null(b_2_hat)) {
    return(list(postmean = NA_real_, uncertainty = NA_real_,
                expected_pve = NA_real_))
  }

  if (!is.null(X)) {
    # susine fits coefficients (b_hat) on the COLUMN-STANDARDIZED scale: it
    # standardizes X internally (center = colMeans, scale = colSds) and
    #   fitted_y = compute_Xb(X, b) = ((X - center) / scale) %*% b.
    # So the quadratic forms must use the standardized design; a naive X %*% b on
    # the raw genotype matrix is wrong by the per-column scaling. Prefer explicit
    # scaling attributes if present, else standardize to unit (n-1) variance to
    # match susine (idempotent if X is already standardized).
    ctr <- attr(X, "scaled:center")
    scl <- attr(X, "scaled:scale")
    Xm <- as.matrix(X)
    if (is.null(ctr) || is.null(scl) ||
        length(ctr) != ncol(Xm) || length(scl) != ncol(Xm)) {
      ctr <- colMeans(Xm)
      scl <- apply(Xm, 2L, stats::sd)
    }
    scl[!is.finite(scl) | scl == 0] <- 1
    Xm <- sweep(sweep(Xm, 2L, ctr, "-"), 2L, scl, "/")
    if (is.null(vy)) {
      vy <- stats::var(as.numeric(y))
    }
    # var(X %*% v) = v' [crossprod(Xc)/(n-1)] v, matching stats::var(y)'s (n-1).
    quad_form <- function(v) stats::var(as.numeric(Xm %*% v))
    gram_diag <- apply(Xm, 2L, stats::var)
  } else if (!is.null(R)) {
    Rm <- as.matrix(R)
    if (is.null(vy)) {
      vy <- 1
    }
    quad_form <- function(v) as.numeric(crossprod(v, Rm %*% v))
    gram_diag <- diag(Rm)
  } else {
    stop("hg2_components: provide either X (individual-level) or R (RSS).")
  }

  hg2_components_from_moments(b_hat, b_2_hat, quad_form, gram_diag, vy)
}

#' Within-fit uncertainty correction only (light path for ensemble threading).
#'
#' Returns the scalar `tr(Sigma Cov(b|y)) / vy`; `postmean` and `between-fit`
#' terms are recovered downstream from the fitted-y vectors, so only this scalar
#' must be carried per fit. NA if the fit lacks `b_2_hat`.
#'
#' @inheritParams hg2_components
#' @return Scalar uncertainty component (fraction of `var(y)`), or NA.
#' @keywords internal
hg2_uncertainty_scalar <- function(fit, X = NULL, y = NULL, R = NULL, vy = NULL) {
  hg2_components(fit, X = X, y = y, R = R, vy = vy)$uncertainty
}
