# Simulation data generation helpers ----------------------------------------

#' Column-wise centering and scaling of a numeric matrix.
#'
#' This mirrors the behaviour of `scale()` but keeps the input type and avoids
#' copying when possible. Columns with zero variance are left unchanged.
#'
#' @param X Numeric matrix.
#' @param center Logical; subtract column means when TRUE.
#' @param scale Logical; divide by column standard deviations when TRUE.
#'
#' @return Matrix with standardized columns.
#' @keywords internal
standardize_x <- function(X, center = TRUE, scale = TRUE) {
  stopifnot(is.matrix(X))
  if (!center && !scale) {
    return(X)
  }
  cm <- if (center) matrixStats::colMeans2(X) else rep(0, ncol(X))
  if (scale) {
    csd <- matrixStats::colSds(X)
    csd[csd == 0] <- 1
  } else {
    csd <- rep(1, ncol(X))
  }
  sweep(sweep(X, 2, cm, `-`), 2, csd, `/`)
}

#' Simulate sparse effect sizes for p SNPs.
#'
#' @param p Number of SNPs.
#' @param p_star Number of causal SNPs to activate.
#' @param effect_sd Standard deviation of true effects.
#' @param seed Optional seed for reproducibility.
#'
#' @return List with `beta` vector and `causal_idx` integer indices.
#' @keywords internal
simulate_effect_sizes <- function(p,
                                  p_star,
                                  effect_sd = 1,
                                  seed = NULL) {
  stopifnot(p > 0, p_star >= 0)
  if (!is.null(seed)) {
    set.seed(seed)
  }
  causal_idx <- if (p_star > 0) {
    sort(sample.int(p, size = min(p_star, p), replace = FALSE))
  } else integer(0)
  beta <- numeric(p)
  if (length(causal_idx)) {
    beta[causal_idx] <- stats::rnorm(length(causal_idx), mean = 0, sd = effect_sd)
  }
  list(beta = beta, causal_idx = causal_idx)
}

#' Construct mu_0 and sigma_0_2 priors from effect sizes using annotation knobs.
#'
#' @param beta Numeric vector of true effects.
#' @param annotation_r2 Target fraction of causal effect variance captured by the annotation.
#' @param inflate_match Scalar controlling non-causal annotation variance relative to the
#'   theoretical causal annotation variance. inflate_match = 0 -> no non-causal annotations;
#'   inflate_match = 1 -> non-causal variance equals theoretical causal variance;
#'   inflate_match = 2 -> non-causal variance is 2x the theoretical causal variance.
#' @param gamma_shrink Optional shrinkage slope used when converting annotations to prior variances.
#' @param base_sigma2 Optional baseline prior variance.
#' @param effect_sd Nominal standard deviation used for causal effects (fallback when beta variance is zero).
#'
#' @return List with `mu_0`, `sigma_0_2`, and `observed_r2`.
#' @keywords internal
simulate_priors <- function(beta,
                            annotation_r2,
                            inflate_match,
                            gamma_shrink = NA_real_,
                            base_sigma2 = NULL,
                            effect_sd = NULL) {
  p <- length(beta)
  causal_idx <- which(beta != 0)
  noncausal_idx <- setdiff(seq_len(p), causal_idx)

  if (is.null(effect_sd) || !is.finite(effect_sd) || effect_sd <= 0) {
    effect_sd <- sqrt(max(stats::var(beta[causal_idx]), stats::var(beta), 1))
  }
  # Handle NULL and NA for annotation_r2 and inflate_match
  if (is.null(annotation_r2) || is.na(annotation_r2)) annotation_r2 <- 0
  if (is.null(inflate_match) || is.na(inflate_match)) inflate_match <- 0
  annotation_r2 <- pmin(pmax(annotation_r2, 0), 1)
  inflate_match <- max(inflate_match, 0)

  causal_var <- stats::var(beta[causal_idx])
  if (is.na(causal_var) || causal_var == 0) {
    causal_var <- effect_sd^2
  }

  if (annotation_r2 <= 0) {
    noise_var <- causal_var
  } else if (annotation_r2 >= 1) {
    noise_var <- 0
  } else {
    noise_var <- causal_var * (1 - annotation_r2) / annotation_r2
  }

  # Compute theoretical causal annotation variance (independent of sample size)
  # Var(mu_0[causal]) = causal_var + noise_var = causal_var / annotation_r2
  theoretical_causal_mu0_var <- if (annotation_r2 <= 0) {
    causal_var
  } else if (annotation_r2 >= 1) {
    causal_var
  } else {
    causal_var / annotation_r2
  }

  mu_0 <- numeric(p)
  if (length(causal_idx)) {
    mu_0[causal_idx] <- beta[causal_idx] + stats::rnorm(length(causal_idx), sd = sqrt(noise_var))
  }
  if (length(noncausal_idx)) {
    # Non-causal variance scales with inflate_match relative to theoretical causal variance
    # inflate_match = 0: variance = 0 (no annotation signal)
    # inflate_match = 1: variance = theoretical_causal_mu0_var (same as causal)
    # inflate_match = 2: variance = 2 * theoretical_causal_mu0_var (twice causal)
    noncausal_var <- inflate_match * theoretical_causal_mu0_var
    mu_0[noncausal_idx] <- stats::rnorm(length(noncausal_idx), sd = sqrt(noncausal_var))
  }

  observed_r2 <- NA_real_
  if (length(causal_idx) > 1 && stats::sd(beta[causal_idx]) > 0 && stats::sd(mu_0[causal_idx]) > 0) {
    observed_r2 <- stats::cor(mu_0[causal_idx], beta[causal_idx])^2
  }

  base_sigma2 <- base_sigma2 %||% causal_var
  if (is.na(base_sigma2) || base_sigma2 <= 0) {
    base_sigma2 <- effect_sd^2
  }

  if (!is.null(gamma_shrink) && is.finite(gamma_shrink)) {
    scale <- if (effect_sd > 0) effect_sd else sqrt(causal_var)
    sigma_0_2 <- base_sigma2 * exp(-gamma_shrink * abs(mu_0) / scale)
  } else {
    sigma_0_2 <- rep(base_sigma2, p)
  }

  list(
    mu_0 = mu_0,
    sigma_0_2 = sigma_0_2,
    observed_r2 = observed_r2
  )
}

#' Simulate phenotype with a target noise fraction.
#'
#' @param X Design matrix.
#' @param beta Effect sizes.
#' @param noise_fraction Fraction of variance attributed to noise (between 0 and 1).
#' @param seed Optional seed.
#'
#' @return List with `y` and implied residual variance `sigma2`.
#' @keywords internal
simulate_phenotype <- function(X, beta, noise_fraction, seed = NULL) {
  if (!is.null(seed)) {
    set.seed(seed)
  }
  signal <- as.vector(X %*% beta)
  signal_var <- stats::var(signal)
  if (is.na(signal_var) || signal_var == 0) {
    signal_var <- max(stats::var(beta), 1e-3)
  }
  noise_fraction <- max(min(noise_fraction, 0.999), 0)
  if (noise_fraction == 1) {
    sigma2 <- 1
    y <- stats::rnorm(nrow(X), mean = 0, sd = 1)
  } else {
    sigma2 <- signal_var * (noise_fraction / (1 - noise_fraction))
    y <- signal + stats::rnorm(nrow(X), mean = 0, sd = sqrt(sigma2))
  }
  list(y = as.numeric(y), sigma2 = sigma2)
}

#' Generate a single simulation dataset for a run specification.
#'
#' @param spec Named list or data.frame row containing simulation controls.
#'   Required fields: `seed`, `p_star`, `y_noise`. Optional prior controls:
#'   `annotation_r2`, `inflate_match`, `gamma_shrink`. Additional optional inputs:
#'   `effect_sd`, `standardize_X`.
#' @param base_X Optional matrix to reuse instead of loading the package data.
#'
#' @return List with design matrix, phenotype, ground-truth beta, priors, and
#'   metadata required by the modelling pipeline.
#' @export
generate_simulation_data <- function(spec,
                                     base_X = NULL) {
  if (is.null(base_X)) {
    data_env <- new.env(parent = emptyenv())
    utils::data("SuSiE_N3_X", package = "test_susine", envir = data_env)
    base_X <- get("SuSiE_N3_X", envir = data_env)
  }
  X <- as.matrix(base_X)
  if (resolve_flag(spec$standardize_X, FALSE)) {
    X <- standardize_x(X)
  }

  seed <- spec$seed %||% 1L
  effects <- simulate_effect_sizes(
    p = ncol(X),
    p_star = spec$p_star %||% 5L,
    effect_sd = spec$effect_sd %||% 1,
    seed = seed
  )
  phenotype <- simulate_phenotype(
    X = X,
    beta = effects$beta,
    noise_fraction = spec$y_noise %||% 0.5,
    seed = seed + 1L
  )
  priors <- simulate_priors(
    beta = effects$beta,
    annotation_r2 = spec$annotation_r2 %||% 0,
    inflate_match = spec$inflate_match %||% 0,
    gamma_shrink = spec$gamma_shrink,
    base_sigma2 = stats::var(phenotype$y),
    effect_sd = spec$effect_sd %||% 1
  )

  list(
    X = X,
    y = phenotype$y,
    beta = effects$beta,
    sigma2 = phenotype$sigma2,
    mu_0 = priors$mu_0,
    sigma_0_2 = priors$sigma_0_2,
    causal_idx = effects$causal_idx,
    annotation_r2_observed = priors$observed_r2,
    annotation_r2_target = spec$annotation_r2 %||% NA_real_,
    inflate_match = spec$inflate_match %||% NA_real_,
    gamma_shrink = spec$gamma_shrink %||% NA_real_
  )
}
