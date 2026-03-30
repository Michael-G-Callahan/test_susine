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

scale_tier_to_energy <- function(beta, idx, target_energy) {
  if (!length(idx) || !is.finite(target_energy) || target_energy <= 0) {
    return(beta)
  }
  current <- sum(beta[idx]^2)
  if (!is.finite(current) || current <= 0) {
    return(beta)
  }
  beta[idx] <- beta[idx] * sqrt(target_energy / current)
  beta
}

#' Simulate oligogenic effect sizes with sparse/oligogenic/polygenic tiers.
#'
#' The generator allocates effects across three tiers and rescales each tier to
#' a target L2-energy fraction. This keeps tier strength interpretable while
#' preserving random draw variability within each tier.
#'
#' @param p Number of SNPs.
#' @param p_star Baseline sparse-causal count (used to size tiers).
#' @param effect_sd Baseline effect SD.
#' @param seed Optional seed.
#' @param tier_h2 Named numeric vector of target tier fractions.
#'
#' @return List with `beta`, `causal_idx`, and `effect_tier` labels.
#' @keywords internal
simulate_effect_sizes_oligogenic <- function(p,
                                             p_star,
                                             effect_sd = 1,
                                             seed = NULL,
                                             tier_h2 = c(sparse = 0.6, oligogenic = 0.3, polygenic = 0.1)) {
  stopifnot(p > 0, p_star >= 0)
  if (!is.null(seed)) {
    set.seed(seed)
  }
  all_idx <- seq_len(p)
  p_star <- as.integer(max(1L, min(as.integer(p_star), p)))
  remaining <- all_idx

  sparse_idx <- sort(sample(remaining, size = p_star, replace = FALSE))
  remaining <- setdiff(remaining, sparse_idx)

  n_oligo <- min(length(remaining), max(1L, as.integer(ceiling(p_star / 2))))
  oligo_idx <- if (n_oligo > 0L) sort(sample(remaining, size = n_oligo, replace = FALSE)) else integer(0)
  remaining <- setdiff(remaining, oligo_idx)

  n_poly <- min(length(remaining), max(1L, as.integer(ceiling(1.5 * p_star))))
  poly_idx <- if (n_poly > 0L) sort(sample(remaining, size = n_poly, replace = FALSE)) else integer(0)

  beta <- numeric(p)
  if (length(sparse_idx)) {
    beta[sparse_idx] <- stats::rnorm(length(sparse_idx), mean = 0, sd = effect_sd)
  }
  if (length(oligo_idx)) {
    mix <- stats::rbinom(length(oligo_idx), size = 1, prob = 0.35)
    sd_vec <- ifelse(mix == 1L, effect_sd * 0.8, effect_sd * 0.25)
    beta[oligo_idx] <- stats::rnorm(length(oligo_idx), mean = 0, sd = sd_vec)
  }
  if (length(poly_idx)) {
    beta[poly_idx] <- stats::rnorm(length(poly_idx), mean = 0, sd = effect_sd * 0.08)
  }

  frac <- as.numeric(tier_h2)
  names(frac) <- names(tier_h2)
  frac <- frac / sum(frac)
  total_energy <- sum(beta^2)
  if (!is.finite(total_energy) || total_energy <= 0) {
    total_energy <- effect_sd^2 * max(length(sparse_idx), 1L)
  }

  beta <- scale_tier_to_energy(beta, sparse_idx, frac["sparse"] * total_energy)
  beta <- scale_tier_to_energy(beta, oligo_idx, frac["oligogenic"] * total_energy)
  beta <- scale_tier_to_energy(beta, poly_idx, frac["polygenic"] * total_energy)

  effect_tier <- rep("none", p)
  effect_tier[sparse_idx] <- "sparse"
  effect_tier[oligo_idx] <- "oligogenic"
  effect_tier[poly_idx] <- "polygenic"
  causal_idx <- sort(c(sparse_idx, oligo_idx, poly_idx))

  list(beta = beta, causal_idx = causal_idx, effect_tier = effect_tier)
}

# SuSiE 2.0 architectures ------------------------------------------------

#' Simulate sparse effect sizes matching the SuSiE 2.0 paper.
#'
#' All causal variants receive fixed effect beta = 1 (not random). Noise is
#' calibrated via h2 rather than noise_fraction — see [simulate_phenotype_h2()].
#'
#' @param p Number of SNPs.
#' @param p_star Number of causal SNPs.
#' @param seed Optional seed for reproducibility.
#'
#' @return List with `beta` and `causal_idx`.
#' @keywords internal
simulate_effect_sizes_susie2_sparse <- function(p, p_star, seed = NULL) {
  stopifnot(p > 0, p_star >= 0)
  if (!is.null(seed)) set.seed(seed)
  causal_idx <- if (p_star > 0) {
    sort(sample.int(p, size = min(p_star, p), replace = FALSE))
  } else {
    integer(0)
  }
  beta <- numeric(p)
  if (length(causal_idx)) beta[causal_idx] <- 1
  list(beta = beta, causal_idx = causal_idx)
}

#' Simulate oligogenic effect sizes matching the SuSiE 2.0 paper.
#'
#' Three tiers with fixed variant counts: sparse (3), oligogenic (5), and
#' polygenic (15). Each tier's L2 energy is rescaled to the target heritability
#' fraction. The oligogenic tier uses a two-component Gaussian mixture.
#'
#' @param p Number of SNPs.
#' @param seed Optional seed for reproducibility.
#' @param tier_config Named list of tier specifications. Each element is a list
#'   with `k` (number of variants) and `h2_frac` (fraction of total h2).
#'
#' @return List with `beta`, `causal_idx`, and `effect_tier`.
#' @keywords internal
simulate_effect_sizes_susie2_oligogenic <- function(
    p,
    seed = NULL,
    tier_config = list(
      sparse     = list(k = 3L,  h2_frac = 0.50),
      oligogenic = list(k = 5L,  h2_frac = 0.35),
      polygenic  = list(k = 15L, h2_frac = 0.15)
    )) {
  stopifnot(p > 0)
  if (!is.null(seed)) set.seed(seed)

  n_sparse <- tier_config$sparse$k
  n_oligo  <- tier_config$oligogenic$k
  n_poly   <- tier_config$polygenic$k
  n_total  <- n_sparse + n_oligo + n_poly
  stopifnot(n_total <= p)

  # Sample non-overlapping causal indices
  all_causal <- sort(sample.int(p, size = n_total, replace = FALSE))
  sparse_idx <- all_causal[seq_len(n_sparse)]
  oligo_idx  <- all_causal[n_sparse + seq_len(n_oligo)]
  poly_idx   <- all_causal[n_sparse + n_oligo + seq_len(n_poly)]

  beta <- numeric(p)

  # Sparse: large effects ~ N(0, 1)
  beta[sparse_idx] <- stats::rnorm(n_sparse, mean = 0, sd = 1)

  # Oligogenic: 2-component Gaussian mixture (35% large / 65% small)
  if (n_oligo > 0L) {
    mix <- stats::rbinom(n_oligo, size = 1, prob = 0.35)
    sd_vec <- ifelse(mix == 1L, 0.8, 0.25)
    beta[oligo_idx] <- stats::rnorm(n_oligo, mean = 0, sd = sd_vec)
  }

  # Polygenic: small effects ~ N(0, 0.08^2)
  if (n_poly > 0L) {
    beta[poly_idx] <- stats::rnorm(n_poly, mean = 0, sd = 0.08)
  }

  # Rescale each tier to target energy fractions
  frac_sparse <- tier_config$sparse$h2_frac
  frac_oligo  <- tier_config$oligogenic$h2_frac
  frac_poly   <- tier_config$polygenic$h2_frac
  total_frac  <- frac_sparse + frac_oligo + frac_poly
  frac_sparse <- frac_sparse / total_frac
  frac_oligo  <- frac_oligo / total_frac
  frac_poly   <- frac_poly / total_frac

  total_energy <- sum(beta^2)
  if (!is.finite(total_energy) || total_energy <= 0) total_energy <- 1

  beta <- scale_tier_to_energy(beta, sparse_idx, frac_sparse * total_energy)
  beta <- scale_tier_to_energy(beta, oligo_idx, frac_oligo * total_energy)
  beta <- scale_tier_to_energy(beta, poly_idx, frac_poly * total_energy)

  effect_tier <- rep("none", p)
  effect_tier[sparse_idx] <- "sparse"
  effect_tier[oligo_idx]  <- "oligogenic"
  effect_tier[poly_idx]   <- "polygenic"
  causal_idx <- sort(c(sparse_idx, oligo_idx, poly_idx))

  list(beta = beta, causal_idx = causal_idx, effect_tier = effect_tier)
}

#' Simulate phenotype with heritability-calibrated noise (SuSiE 2.0 style).
#'
#' Instead of specifying a noise fraction, the residual variance is calibrated
#' to achieve a target total heritability (h2). Exactly one of
#' `h2_snp_per_causal` or `h2_total` must be non-NA.
#'
#' @param X Design matrix (n x p).
#' @param beta True effect sizes (length p).
#' @param h2_snp_per_causal Per-SNP heritability; total h2 = h2_snp * n_causal.
#' @param h2_total Total heritability directly (for oligogenic architectures).
#' @param seed Optional seed.
#'
#' @return List with `y`, `sigma2`, and `h2_realized`.
#' @keywords internal
simulate_phenotype_h2 <- function(X, beta,
                                  h2_snp_per_causal = NA_real_,
                                  h2_total = NA_real_,
                                  seed = NULL) {
  if (!is.null(seed)) set.seed(seed)

  has_snp <- !is.null(h2_snp_per_causal) && is.finite(h2_snp_per_causal)
  has_total <- !is.null(h2_total) && is.finite(h2_total)
  stopifnot(has_snp || has_total, !(has_snp && has_total))

  h2 <- if (has_snp) {
    h2_snp_per_causal * sum(beta != 0)
  } else {
    h2_total
  }
  h2 <- max(min(h2, 0.999), 0.001)

  signal <- as.vector(X %*% beta)
  var_signal <- stats::var(signal)
  if (is.na(var_signal) || var_signal <= 0) var_signal <- 1e-3

  sigma2 <- var_signal * (1 - h2) / h2
  y <- signal + stats::rnorm(nrow(X), mean = 0, sd = sqrt(sigma2))

  list(
    y = as.numeric(y),
    sigma2 = sigma2,
    h2_realized = var_signal / (var_signal + sigma2)
  )
}

#' Construct mu_0 and sigma_0_2 priors from effect sizes using annotation knobs.
#'
#' @param beta Numeric vector of true effects.
#' @param annotation_r2 Target squared correlation between the causal annotation
#'   values and the true causal effects. The construction preserves the causal
#'   annotation scale, so `annotation_r2 = 1` yields an exact match to the true
#'   causal effects and `annotation_r2 = 0` yields centered noise at the same
#'   scale.
#' @param inflate_match Scalar controlling the centered second moment of
#'   non-causal annotation values relative to the causal annotation scale.
#'   `inflate_match = 0` -> no non-causal annotations; `inflate_match = 1` ->
#'   non-causal annotations have the same centered scale as causal annotations;
#'   `inflate_match = 2` -> non-causal annotations have sqrt(2)x the causal
#'   scale.
#' @param gamma_shrink Deprecated. Ignored when provided.
#' @param base_sigma2 Optional baseline prior variance.
#' @param effect_sd Nominal standard deviation used for causal effects (fallback when beta variance is zero).
#' @param seed Optional seed for reproducibility of annotation draws.
#'
#' @return List with `mu_0`, `sigma_0_2`, and `observed_r2`.
#' @keywords internal
center_and_scale <- function(x, target_scale = 1, orthogonal_to = NULL) {
  x <- as.numeric(x)
  n <- length(x)
  if (!n) {
    return(numeric(0))
  }

  x <- x - mean(x)
  if (!is.null(orthogonal_to)) {
    ref <- as.numeric(orthogonal_to)
    if (length(ref) == n) {
      denom <- sum(ref^2)
      if (is.finite(denom) && denom > 0) {
        x <- x - sum(x * ref) / denom * ref
      }
    }
  }

  rms_x <- sqrt(mean(x^2))
  if (!is.finite(rms_x) || rms_x <= 0) {
    return(rep(0, n))
  }
  x / rms_x * target_scale
}

draw_centered_noise <- function(n, target_scale = 1, orthogonal_to = NULL) {
  if (n <= 0) {
    return(numeric(0))
  }
  if (n == 1L) {
    return(0)
  }

  for (attempt in seq_len(8L)) {
    noise <- center_and_scale(
      stats::rnorm(n),
      target_scale = target_scale,
      orthogonal_to = orthogonal_to
    )
    if (sqrt(mean(noise^2)) > 0) {
      return(noise)
    }
  }

  stop("Unable to generate non-degenerate centered annotation noise.")
}

simulate_priors <- function(beta,
                            annotation_r2,
                            inflate_match,
                            gamma_shrink = NA_real_,
                            base_sigma2 = NULL,
                            effect_sd = NULL,
                            seed = NULL) {
  if (!is.null(seed) && is.finite(seed)) {
    set.seed(as.integer(seed))
  }
  if (!is.null(gamma_shrink) && is.finite(gamma_shrink)) {
    warning("`gamma_shrink` is deprecated and ignored.")
  }
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

  mu_0 <- numeric(p)
  if (length(causal_idx)) {
    causal_beta <- as.numeric(beta[causal_idx])
    causal_mean <- mean(causal_beta)
    causal_centered <- causal_beta - causal_mean
    causal_scale <- sqrt(mean(causal_centered^2))
    if (!is.finite(causal_scale) || causal_scale <= 0) {
      causal_scale <- sqrt(causal_var)
    }

    if (length(causal_idx) == 1L || !is.finite(causal_scale) || causal_scale <= 0) {
      mu_0[causal_idx] <- causal_beta
    } else {
      signal_unit <- causal_centered / causal_scale
      noise_unit <- draw_centered_noise(
        length(causal_idx),
        target_scale = 1,
        orthogonal_to = signal_unit
      )
      causal_annotation_centered <- causal_scale * (
        sqrt(annotation_r2) * signal_unit +
          sqrt(1 - annotation_r2) * noise_unit
      )
      mu_0[causal_idx] <- causal_mean + causal_annotation_centered
    }
  }

  theoretical_causal_mu0_var <- if (length(causal_idx) > 1L) {
    mean((mu_0[causal_idx] - mean(mu_0[causal_idx]))^2)
  } else {
    causal_var
  }

  if (length(noncausal_idx)) {
    noncausal_scale <- sqrt(inflate_match) * sqrt(theoretical_causal_mu0_var)
    mu_0[noncausal_idx] <- draw_centered_noise(
      length(noncausal_idx),
      target_scale = noncausal_scale
    )
  }

  observed_r2 <- NA_real_
  if (length(causal_idx) > 1 && stats::sd(beta[causal_idx]) > 0 && stats::sd(mu_0[causal_idx]) > 0) {
    observed_r2 <- stats::cor(mu_0[causal_idx], beta[causal_idx])^2
  }

  base_sigma2 <- base_sigma2 %||% causal_var
  if (is.na(base_sigma2) || base_sigma2 <= 0) {
    base_sigma2 <- effect_sd^2
  }

  sigma_0_2 <- rep(base_sigma2, p)

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
#'   Required fields: `phenotype_seed`, `p_star`, `y_noise`. Optional prior controls:
#'   `annotation_r2`, `inflate_match`. Additional optional inputs:
#'   `effect_sd`, `standardize_X`, `architecture`.
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

  seed <- spec$phenotype_seed %||% spec$seed %||% 1L
  architecture <- as.character(spec$architecture %||% "sparse")
  effects <- if (identical(architecture, "oligogenic")) {
    simulate_effect_sizes_oligogenic(
      p = ncol(X),
      p_star = spec$p_star %||% 5L,
      effect_sd = spec$effect_sd %||% 1,
      seed = seed
    )
  } else if (identical(architecture, "susie2_sparse")) {
    simulate_effect_sizes_susie2_sparse(
      p = ncol(X),
      p_star = spec$p_star %||% 3L,
      seed = seed
    )
  } else if (identical(architecture, "susie2_oligogenic")) {
    simulate_effect_sizes_susie2_oligogenic(
      p = ncol(X),
      seed = seed
    )
  } else {
    simulate_effect_sizes(
      p = ncol(X),
      p_star = spec$p_star %||% 5L,
      effect_sd = spec$effect_sd %||% 1,
      seed = seed
    )
  }

  # SuSiE 2.0 architectures use h2-calibrated noise; others use noise_fraction
  phenotype <- if (architecture %in% c("susie2_sparse", "susie2_oligogenic")) {
    simulate_phenotype_h2(
      X = X,
      beta = effects$beta,
      h2_snp_per_causal = if (identical(architecture, "susie2_sparse")) {
        spec$h2_snp_per_causal %||% 0.03
      } else {
        NA_real_
      },
      h2_total = if (identical(architecture, "susie2_oligogenic")) {
        spec$h2_total %||% 0.25
      } else {
        NA_real_
      },
      seed = seed + 1L
    )
  } else {
    simulate_phenotype(
      X = X,
      beta = effects$beta,
      noise_fraction = spec$y_noise %||% 0.5,
      seed = seed + 1L
    )
  }
  priors <- simulate_priors(
    beta = effects$beta,
    annotation_r2 = spec$annotation_r2 %||% 0,
    inflate_match = spec$inflate_match %||% 0,
    gamma_shrink = spec$gamma_shrink,
    base_sigma2 = stats::var(phenotype$y),
    effect_sd = spec$effect_sd %||% 1,
    seed = spec$annotation_seed %||% NULL
  )

  list(
    X = X,
    y = phenotype$y,
    beta = effects$beta,
    sigma2 = phenotype$sigma2,
    mu_0 = priors$mu_0,
    sigma_0_2 = priors$sigma_0_2,
    causal_idx = effects$causal_idx,
    effect_tier = effects$effect_tier %||% rep("sparse", ncol(X)),
    architecture = architecture,
    annotation_r2_observed = priors$observed_r2,
    annotation_r2_target = spec$annotation_r2 %||% NA_real_,
    inflate_match = spec$inflate_match %||% NA_real_,
    h2_realized = phenotype$h2_realized %||% NA_real_
  )
}
