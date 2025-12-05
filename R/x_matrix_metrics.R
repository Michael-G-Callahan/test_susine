# x_matrix_metrics.R
# Functions for computing X-matrix difficulty metrics and relating them to AUPRC
# Based on theory linking LD structure to fine-mapping performance
# ------------------------------------------------------------------------------

#' @importFrom stats cor quantile sd
#' @importFrom Matrix Matrix
NULL

# ==============================================================================
# SECTION 1: Core X-matrix metric computations
# ==============================================================================

#' Compute mid-energy M1 (linear)
#' Higher values indicate more mass in ambiguous mid-LD range
#' @param R Correlation matrix
#' @return Scalar M1 value
#' @keywords internal
mid_energy_M1 <- function(R) {

  A <- abs(R)
  diag(A) <- 0
  ut <- A[upper.tri(A)]
  mean(ut * (1 - ut)) * 2
}

#' Compute mid-energy M2 (quadratic)
#' Sharper emphasis on moderate-high mid-LD
#' @param R Correlation matrix
#' @return Scalar M2 value
#' @keywords internal
mid_energy_M2 <- function(R) {
  A2 <- R * R
  diag(A2) <- 0
  ut <- A2[upper.tri(A2)]
  mean(ut * (1 - ut)) * 2
}

#' Compute distance-weighted mid-energy
#' Detects non-local "plaid" LD patterns
#' @param R Correlation matrix
#' @return Scalar M1_dist value
#' @keywords internal
mid_energy_dist_weighted <- function(R) {
  p <- ncol(R)
  A <- abs(R)
  diag(A) <- 0
  
  # Distance weights: phi(d) = d/(d+1)
  dists <- outer(1:p, 1:p, function(i, j) abs(i - j))
  phi <- dists / (dists + 1)
  diag(phi) <- 0
  
  num <- sum(A * (1 - A) * phi) / 2

  den <- p * (p - 1) / 2
  num / den
}

#' Compute eigenvalue-based summaries
#' @param R Correlation matrix
#' @return Tibble with spectral metrics
#' @keywords internal
eigs_summaries <- function(R) {
  ev <- eigen(R, symmetric = TRUE, only.values = TRUE)$values
  p <- length(ev)
  lam1 <- max(ev)
  
  # Participation ratio (effective rank normalized)
  PR <- (sum(ev)^2) / sum(ev^2) / p
  
  # Spectral entropy (normalized)
  lam_norm <- ev / sum(ev)
  lam_norm <- pmax(lam_norm, .Machine$double.eps)  
  H <- -sum(lam_norm * log(lam_norm))
  spec_entropy_norm <- exp(H) / p
  
  # Stable rank
  stable_rank <- sum(ev^2) / (lam1^2)
  
  # Spectral gap
  ev_sorted <- sort(ev, decreasing = TRUE)
  spectral_gap <- ev_sorted[1] - ev_sorted[2]
  
  tibble::tibble(
    lambda_max_frac = lam1 / p,
    participation_ratio = PR,
    spectral_entropy_norm = spec_entropy_norm,
    stable_rank = stable_rank,
    spectral_gap = spectral_gap
  )
}

#' Compute mutual coherence (max off-diagonal |r|)
#' @param R Correlation matrix
#' @return Scalar mu value
#' @keywords internal
mutual_coherence <- function(R) {
  A <- abs(R)
  diag(A) <- 0
  max(A)
}

#' Compute row concentration (Herfindahl index per row)
#' High = mass on few neighbors (blocky); Low = diffuse
#' @param R Correlation matrix
#' @return Tibble with mean and q90 of row concentration
#' @keywords internal
row_concentration <- function(R) {
  A <- abs(R)
  diag(A) <- 0
  row_sums <- rowSums(A)
  row_sums <- pmax(row_sums, .Machine$double.eps)
  w <- A / row_sums
  H <- rowSums(w^2)
  tibble::tibble(
    row_conc_mean = mean(H),
    row_conc_q90 = unname(stats::quantile(H, 0.9))
  )
}

#' Compute Toeplitz (AR1) fit error
#' Low = near-banded (simple); High = complex structure
#' @param R Correlation matrix
#' @return Tibble with rho_hat and relative Frobenius error
#' @keywords internal
toeplitz_fit_error <- function(R) {
  p <- ncol(R)
  
  # Fit AR(1) by matching first subdiagonal
  sub1 <- R[row(R) == col(R) + 1]
  rho_hat <- mean(sub1)
  
  # Build Toeplitz matrix
  T_rho <- outer(1:p, 1:p, function(a, b) rho_hat^(abs(a - b)))
  
  # Relative error
  err <- norm(R - T_rho, type = "F") / norm(R, type = "F")
  
  tibble::tibble(
    ar1_rho_hat = rho_hat,
    toeplitz_rel_F_error = err
  )
}

# ==============================================================================
# SECTION 1b: Axis A metrics - Basin Size / Resolution Difficulty
# "How big are the LD neighborhoods that SuSiE will find?"
# ==============================================================================

#' Compute Effective Cluster Size (ECS)
#'
#' Measures average effective neighborhood size using Herfindahl inverse.
#' Uses alpha=2 weighting to emphasize strong correlations.
#' High ECS = big/strong clusters = IBSS finds them but CSs are large.
#'
#' @param R Correlation matrix
#' @param alpha Weighting exponent (default 2, emphasizes strong correlations)
#' @return Scalar ECS value (median across SNPs)
#' @keywords internal
effective_cluster_size <- function(R, alpha = 2) {
  A <- abs(R)
  diag(A) <- 0
  
  # Compute row-wise weights: w_ij = A_ij^alpha / sum_k A_ik^alpha
  A_alpha <- A^alpha
  row_sums <- rowSums(A_alpha)
  row_sums <- pmax(row_sums, .Machine$double.eps)
  W <- A_alpha / row_sums
  
  # ECS_i = 1 / sum_j w_ij^2 (inverse Herfindahl)
  ecs_i <- 1 / rowSums(W^2)
  
  # Return median (robust to outliers)
  median(ecs_i)
}

#' Compute High-LD Mass
#'
#' Measures mass in strong-LD pairs using gamma=2 to emphasize tail.
#' High value = big tight blocks.
#'
#' @param R Correlation matrix
#' @param gamma Exponent for emphasizing strong correlations (default 2)
#' @return Scalar high-LD mass value
#' @keywords internal
high_ld_mass <- function(R, gamma = 2) {
  A2 <- R * R  # r^2
  A2g <- A2^gamma  # (r^2)^gamma
  diag(A2g) <- 0
  
  # Mean over upper triangle
  mean(A2g[upper.tri(A2g)]) * 2
}

#' Compute top-k eigenvalue fraction
#'
#' Fraction of variance captured by top k eigenvalues.
#' High value = few dominant directions = block structure.
#'
#' @param R Correlation matrix
#' @param k Number of top eigenvalues (default 5)
#' @return Scalar fraction
#' @keywords internal
top_k_ev_frac <- function(R, k = 5) {
  ev <- eigen(R, symmetric = TRUE, only.values = TRUE)$values
  ev <- pmax(ev, 0)  # numerical stability
  ev_sorted <- sort(ev, decreasing = TRUE)
  k <- min(k, length(ev))
  
  sum(ev_sorted[1:k]) / sum(ev)
}

#' Compute block structure metrics via thresholded connected components
#'
#' Thresholds R at tau, finds connected components, returns block statistics.
#' Directly measures "how blocky" the LD structure is.
#'
#' @param R Correlation matrix
#' @param tau Threshold for |r| to define edges
#' @param suffix Suffix for column names (e.g., "_tau50" or "_tau80")
#' @return Tibble with number of blocks, average/max block size, block entropy
#' @keywords internal
block_structure_metrics <- function(R, tau, suffix) {
  p <- ncol(R)
  
  # Create adjacency matrix at threshold tau
  A <- abs(R) >= tau
  diag(A) <- FALSE
  
  # Find connected components via simple BFS
  visited <- rep(FALSE, p)
  block_sizes <- integer(0)
  
  for (start in 1:p) {
    if (visited[start]) next
    
    # BFS from this node
    queue <- start
    component_size <- 0
    while (length(queue) > 0) {
      node <- queue[1]
      queue <- queue[-1]
      if (visited[node]) next
      visited[node] <- TRUE
      component_size <- component_size + 1
      neighbors <- which(A[node, ] & !visited)
      queue <- c(queue, neighbors)
    }
    block_sizes <- c(block_sizes, component_size)
  }
  
  n_blocks <- length(block_sizes)
  avg_block_size <- mean(block_sizes)
  max_block_size <- max(block_sizes)
  
  # Block entropy: high = many similarly-sized blocks; low = one dominant block
  block_probs <- block_sizes / sum(block_sizes)
  block_entropy <- -sum(block_probs * log(block_probs + 1e-10))
  # Normalize by max possible entropy (uniform blocks)
  max_entropy <- log(n_blocks + 1e-10)
  block_entropy_norm <- if (max_entropy > 0) block_entropy / max_entropy else 1
  
  # Build tibble with suffixed names
  tbl <- tibble::tibble(
    a = n_blocks,
    b = avg_block_size,
    c = max_block_size,
    d = max_block_size / p,
    e = block_entropy_norm
  )
  names(tbl) <- paste0(c("n_blocks", "avg_block_size", "max_block_size", 
                         "max_block_frac", "block_entropy_norm"), suffix)
  tbl
}

# ==============================================================================
# SECTION 1c: Axis B metrics - Identifiability / Factorization Mismatch
# "How interchangeable are different SNP neighborhoods?"
# ==============================================================================

#' Compute Neighborhood Similarity Index (NSI)
#'
#' Measures how similar LD "fingerprints" are across SNPs.
#' L2-normalizes rows of |R|, then computes max cosine similarity per row.
#' High NSI = many SNPs have similar correlation neighborhoods = blur problem.
#'
#' @param R Correlation matrix
#' @return Scalar NSI value (mean of max cosine similarities)
#' @keywords internal
neighborhood_similarity_index <- function(R) {
  A <- abs(R)
  diag(A) <- 0
  
  # L2-normalize rows
  row_norms <- sqrt(rowSums(A^2))
  row_norms[row_norms == 0] <- 1  # avoid division by zero
  B <- A / row_norms
  
  # Compute cosine similarity matrix
  S <- B %*% t(B)
  diag(S) <- -Inf  # exclude self-similarity
  
  # Mean of max similarity per row
  mean(apply(S, 1, max))
}

#' Compute off-diagonal means
#' @param R Correlation matrix
#' @return Tibble with mean |r| and mean r^2
#' @keywords internal
offdiag_means <- function(R) {
  A <- abs(R)
  diag(A) <- 0
  A2 <- R * R
  diag(A2) <- 0
  
  tibble::tibble(
    mean_abs_r = mean(A[upper.tri(A)]) * 2,
    mean_r2 = mean(A2[upper.tri(A2)]) * 2
  )
}

#' Compute submodularity bounds via Monte Carlo sparse eigenvalues
#' @param R Correlation matrix
#' @param k Number of causal effects (uses s = 2k)
#' @param n_samples Number of random subsets to sample
#' @param seed Random seed for reproducibility
#' @return Tibble with lambda_min bounds and gamma estimates
#' @keywords internal
submodularity_bounds <- function(R, k, n_samples = 300, seed = 1L) {
  stopifnot(k >= 1)
  p <- ncol(R)
  s <- min(2 * k, p)
  
  set.seed(seed)
  lam_min_mc <- Inf
  
  for (t in seq_len(n_samples)) {
    S <- sort(sample.int(p, s))
    lam_min <- min(eigen(R[S, S, drop = FALSE], symmetric = TRUE, only.values = TRUE)$values)
    if (lam_min < lam_min_mc) lam_min_mc <- lam_min
  }
  
  # Mutual coherence
  A <- abs(R)
  diag(A) <- 0
  mu <- max(A)
  
  # Gershgorin bound: lambda_min >= 1 - (s-1) * mu
  gamma_lb_gersh <- max(0, 1 - (s - 1) * mu)
  
  tibble::tibble(
    k_ref = k,
    s_for_bounds = s,
    mu = mu,
    lambda_min_2k_mc = lam_min_mc,
    gamma_lb_mc = lam_min_mc,
    gamma_lb_gersh = gamma_lb_gersh
  )
}

#' Compute all X-matrix metrics for a single correlation matrix
#'
#' Computes metrics organized by theoretical axes:
#' - Axis A (Basin Size): ecs, high_ld_mass, top5_ev_frac, block structure
#' - Axis B (Identifiability): nsi, M1, M2, M1_dist, lambda_min_2k_mc
#' - Plus spectral and structural descriptors
#'
#' @param R Correlation matrix (p x p)
#' @param k_ref Reference number of causal effects for submodularity bounds
#' @param n_samples Number of MC samples for sparse eigenvalue estimation
#' @return Tibble with all metrics
#' @export
compute_x_metrics_from_R <- function(R, k_ref = 3, n_samples = 200) {
  # Ensure valid correlation matrix
  R[is.na(R)] <- 0
  diag(R) <- 1
  
  # Base metrics (Axis B: mid-energy / blur)
  base_tbl <- tibble::tibble(
    M1 = mid_energy_M1(R),
    M2 = mid_energy_M2(R),
    M1_dist = mid_energy_dist_weighted(R)
  ) %>%
    dplyr::bind_cols(eigs_summaries(R)) %>%
    dplyr::bind_cols(offdiag_means(R)) %>%
    dplyr::bind_cols(row_concentration(R)) %>%
    dplyr::bind_cols(toeplitz_fit_error(R))
  
  # Axis A metrics: Basin size / resolution difficulty
  axis_a_tbl <- tibble::tibble(
    ecs = effective_cluster_size(R, alpha = 2),
    high_ld_mass = high_ld_mass(R, gamma = 2),
    top5_ev_frac = top_k_ev_frac(R, k = 5)
  ) %>%
    dplyr::bind_cols(block_structure_metrics(R, tau = 0.5, suffix = "_tau50")) %>%
    dplyr::bind_cols(block_structure_metrics(R, tau = 0.8, suffix = "_tau80"))
  
  # Axis B metrics: Identifiability / factorization mismatch
  axis_b_tbl <- tibble::tibble(
    nsi = neighborhood_similarity_index(R)
  )
  
  # Submodularity bounds (Axis B: theoretical identifiability limit)
  subm_tbl <- submodularity_bounds(R, k = k_ref, n_samples = n_samples)
  
  dplyr::bind_cols(base_tbl, axis_a_tbl, axis_b_tbl, subm_tbl)
}

#' Compute X-matrix metrics from a genotype matrix file
#'
#' @param x_path Path to file containing X matrix (samples x SNPs).
#'   Supports CSV, TSV, and gzipped versions (.gz)
#' @param k_ref Reference number of causal effects
#' @param n_samples Number of MC samples for sparse eigenvalue
#' @return Tibble with all metrics
#' @export
compute_x_metrics_from_file <- function(x_path, k_ref = 3, n_samples = 200) {
  # Detect file type and read accordingly
  if (grepl("\\.tsv\\.gz$|\\.tsv$", x_path, ignore.case = TRUE)) {
    X_df <- suppressWarnings(readr::read_tsv(x_path, show_col_types = FALSE))
  } else if (grepl("\\.csv\\.gz$|\\.csv$", x_path, ignore.case = TRUE)) {
    X_df <- suppressWarnings(readr::read_csv(x_path, show_col_types = FALSE))
  } else {
    # Try delim auto-detect
    X_df <- suppressWarnings(readr::read_delim(x_path, show_col_types = FALSE, delim = "\t"))
  }
  
  # Check if first column looks like sample IDs (non-numeric)
  first_col <- X_df[[1]]
  if (is.character(first_col) && all(is.na(suppressWarnings(as.numeric(first_col))))) {
    # First column is sample IDs - drop it
    X_df <- X_df[, -1, drop = FALSE]
  }
  
  # Convert to numeric matrix, suppressing coercion warnings from headers
  X <- suppressWarnings(apply(as.matrix(X_df), 2, as.numeric))
  
  # Check for issues
  n_na <- sum(is.na(X))
  if (n_na > 0) {
    pct_na <- 100 * n_na / length(X)
    if (pct_na > 5) {
      warning("High NA rate (", round(pct_na, 1), "%) in matrix from: ", basename(x_path))
    }
  }
  
  R <- suppressWarnings(stats::cor(X, use = "pairwise.complete.obs"))
  
  compute_x_metrics_from_R(R, k_ref = k_ref, n_samples = n_samples)
}

#' Build metrics table for multiple datasets
#'
#' @param matrix_paths_df Data frame with columns `matrix_id` and `matrix_path`
#' @param k_ref Reference k for submodularity bounds
#' @param n_samples MC samples per dataset
#' @param base_dir Optional base directory to prepend to paths (e.g., here())
#' @return Tibble with matrix_id and all metrics
#' @export
build_x_metrics_table <- function(matrix_paths_df, k_ref = 3, n_samples = 200, base_dir = NULL) {
  stopifnot("matrix_id" %in% names(matrix_paths_df))
  stopifnot("matrix_path" %in% names(matrix_paths_df))
  
  n_matrices <- nrow(matrix_paths_df)
  message("Computing X-matrix metrics for ", n_matrices, " datasets...")
  
  purrr::map_dfr(seq_len(n_matrices), function(i) {
    row <- matrix_paths_df[i, ]
    id <- row$matrix_id
    path <- row$matrix_path
    
    # Prepend base_dir if provided
    if (!is.null(base_dir)) {
      path <- file.path(base_dir, path)
    }
    
    if (!file.exists(path)) {
      warning("Missing X file for matrix_id=", id, ": ", path)
      return(tibble::tibble(matrix_id = id, .missing_file = TRUE))
    }
    
    tryCatch({
      m <- compute_x_metrics_from_file(path, k_ref = k_ref, n_samples = n_samples)
      m %>% dplyr::mutate(matrix_id = id, .missing_file = FALSE)
    }, error = function(e) {
      warning("Error processing matrix_id=", id, ": ", e$message)
      tibble::tibble(matrix_id = id, .missing_file = TRUE)
    })
  }, .progress = TRUE)
}

# ==============================================================================
# SECTION 2: AUPRC aggregation schemes
# ==============================================================================

#' Aggregate AUPRC by Dataset × Use Case
#' 
#' Computes single AUPRC per (matrix_id, use_case_id) by pooling
#' all seeds and settings within each combination.
#'
#' @param confusion_df Pre-aggregated confusion matrix data
#' @return Tibble with matrix_id, use_case_id, and AUPRC
#' @export
aggregate_auprc_by_dataset_usecase <- function(confusion_df) {
  confusion_df %>%
    dplyr::group_by(matrix_id, use_case_id, pip_threshold) %>%
    dplyr::summarise(
      n_causal_at_bucket = sum(n_causal_at_bucket, na.rm = TRUE),
      n_noncausal_at_bucket = sum(n_noncausal_at_bucket, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    # Compute AUPRC per (matrix_id, use_case_id)
    dplyr::group_by(matrix_id, use_case_id) %>%
    dplyr::arrange(dplyr::desc(pip_threshold), .by_group = TRUE) %>%
    dplyr::mutate(
      TP = cumsum(n_causal_at_bucket),
      FP = cumsum(n_noncausal_at_bucket),
      total_causal = sum(n_causal_at_bucket),
      precision = dplyr::if_else(TP + FP > 0, TP / (TP + FP), NA_real_),
      recall = dplyr::if_else(total_causal > 0, TP / total_causal, NA_real_)
    ) %>%
    dplyr::summarise(
      AUPRC = compute_auprc_single_internal(precision, recall),
      n_runs = dplyr::n_distinct(pip_threshold),  # proxy for data volume
      .groups = "drop"
    )
}

#' Aggregate AUPRC by Dataset × Use Case × Key Settings
#'
#' Keeps important settings separate to detect interactions.
#'
#' @param confusion_df Pre-aggregated confusion matrix with setting columns
#' @param setting_vars Character vector of setting columns to preserve
#' @return Tibble with matrix_id, use_case_id, settings, and AUPRC
#' @export
aggregate_auprc_by_dataset_usecase_settings <- function(confusion_df,
                                                         setting_vars = c("p_star", "y_noise")) {
  group_vars <- c("matrix_id", "use_case_id", setting_vars, "pip_threshold")
  
  confusion_df %>%
    dplyr::group_by(dplyr::across(dplyr::all_of(group_vars))) %>%
    dplyr::summarise(
      n_causal_at_bucket = sum(n_causal_at_bucket, na.rm = TRUE),
      n_noncausal_at_bucket = sum(n_noncausal_at_bucket, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    dplyr::group_by(dplyr::across(dplyr::all_of(c("matrix_id", "use_case_id", setting_vars)))) %>%
    dplyr::arrange(dplyr::desc(pip_threshold), .by_group = TRUE) %>%
    dplyr::mutate(
      TP = cumsum(n_causal_at_bucket),
      FP = cumsum(n_noncausal_at_bucket),
      total_causal = sum(n_causal_at_bucket),
      precision = dplyr::if_else(TP + FP > 0, TP / (TP + FP), NA_real_),
      recall = dplyr::if_else(total_causal > 0, TP / total_causal, NA_real_)
    ) %>%
    dplyr::summarise(
      AUPRC = compute_auprc_single_internal(precision, recall),
      .groups = "drop"
    )
}

#' Compute residualized AUPRC (Simpson's paradox correction)
#'
#' Fits a baseline model explaining AUPRC from settings, then extracts
#' residuals representing "unexplained" performance variation.
#'
#' @param auprc_df Data frame with AUPRC and setting columns
#' @param setting_vars Columns to use as predictors in baseline model
#' @return Original df with added AUPRC_residual column
#' @export
compute_residualized_auprc <- function(auprc_df,
                                        setting_vars = c("p_star", "y_noise", "use_case_id")) {
  # Build formula: AUPRC ~ setting1 + setting2 + ...
  # Use factors for categorical variables
  df_model <- auprc_df %>%
    dplyr::mutate(dplyr::across(dplyr::where(is.character), as.factor))
  
  formula_str <- paste("AUPRC ~", paste(setting_vars, collapse = " + "))
  
  fit <- stats::lm(stats::as.formula(formula_str), data = df_model)
  
  auprc_df %>%
    dplyr::mutate(
      AUPRC_predicted = stats::predict(fit, newdata = df_model),
      AUPRC_residual = AUPRC - AUPRC_predicted
    )
}

#' Compute within-setting AUPRC ranks
#'
#' For each setting combination, rank datasets by AUPRC.
#' This removes setting confounds entirely.
#'
#' @param auprc_df Data frame with AUPRC and setting columns
#' @param setting_vars Columns defining settings to group within
#' @return Original df with AUPRC_rank and AUPRC_zscore columns
#' @export
compute_within_setting_ranks <- function(auprc_df,
                                          setting_vars = c("p_star", "y_noise", "use_case_id")) {
  auprc_df %>%
    dplyr::group_by(dplyr::across(dplyr::all_of(setting_vars))) %>%
    dplyr::mutate(
      AUPRC_rank = rank(AUPRC) / dplyr::n(),  # percentile rank
      AUPRC_zscore = (AUPRC - mean(AUPRC, na.rm = TRUE)) / 
                     pmax(stats::sd(AUPRC, na.rm = TRUE), 0.001)
    ) %>%
    dplyr::ungroup()
}

#' Internal AUPRC computation (single curve)
#' @keywords internal
compute_auprc_single_internal <- function(precision, recall) {
  df <- data.frame(precision = precision, recall = recall) %>%
    dplyr::filter(!is.na(precision), !is.na(recall)) %>%
    dplyr::arrange(recall)
  
  if (nrow(df) < 2) return(NA_real_)
  
  # Anchor at (recall=0, precision=1)
  if (min(df$recall) > 0) {
    df <- dplyr::bind_rows(
      data.frame(recall = 0, precision = 1),
      df
    )
  }
  
  # Trapezoidal integration
  sum(diff(df$recall) * (utils::head(df$precision, -1) + utils::tail(df$precision, -1)) / 2)
}

# ==============================================================================
# SECTION 3: Plotting X-metrics vs AUPRC
# ==============================================================================

#' Plot X-metric vs AUPRC with smooth fit and R²
#'
#' @param df_joined Data frame with X metrics and AUPRC
#' @param metric_name Column name of metric to plot on x-axis
#' @param auprc_var Column name for AUPRC (allows residuals, ranks, etc.)
#' @param color_var Optional column for color aesthetic
#' @param smoother Either "gam" or "loess"
#' @return ggplot object
#' @export
plot_x_metric_vs_auprc <- function(df_joined,
                                    metric_name,
                                    auprc_var = "AUPRC",
                                    color_var = NULL,
                                    smoother = c("gam", "loess")) {
  smoother <- match.arg(smoother)
  
  dfm <- df_joined %>%
    dplyr::filter(!is.na(.data[[metric_name]]), !is.na(.data[[auprc_var]]))
  
  if (nrow(dfm) < 10) {
    warning("Not enough data points for metric: ", metric_name)
    return(NULL)
  }
  
  # Check for sufficient unique x-values for smoothing
  n_unique <- length(unique(dfm[[metric_name]]))
  use_smooth <- n_unique >= 6
  
  if (!use_smooth) {
    message("Metric '", metric_name, "' has only ", n_unique, 
            " unique values - using linear fit instead of ", smoother)
  }
  
 # Compute R² for annotation
  r2 <- NA_real_
  yhat <- NULL
  
  tryCatch({
    if (use_smooth && smoother == "gam") {
      # Adaptive k based on unique values
      k_use <- min(6, n_unique - 1)
      fit <- mgcv::gam(
        stats::as.formula(paste(auprc_var, "~ s(", metric_name, ", k =", k_use, ")")),
        data = dfm, method = "REML"
      )
      yhat <- stats::predict(fit, newdata = dfm)
    } else if (use_smooth && smoother == "loess") {
      fit <- stats::loess(
        stats::as.formula(paste(auprc_var, "~", metric_name)),
        data = dfm, span = 0.75
      )
      yhat <- stats::predict(fit)
    } else {
      # Fall back to linear model
      fit <- stats::lm(
        stats::as.formula(paste(auprc_var, "~", metric_name)),
        data = dfm
      )
      yhat <- stats::predict(fit)
    }
    
    r2 <- 1 - sum((dfm[[auprc_var]] - yhat)^2, na.rm = TRUE) / 
              sum((dfm[[auprc_var]] - mean(dfm[[auprc_var]], na.rm = TRUE))^2, na.rm = TRUE)
  }, error = function(e) {
    warning("Could not fit model for ", metric_name, ": ", e$message)
  })
  
  r2_lab <- if (is.na(r2)) "R² = NA" else sprintf("R² = %.3f", r2)
  
  # Build plot
  p <- ggplot2::ggplot(dfm, ggplot2::aes(x = .data[[metric_name]], y = .data[[auprc_var]]))
  
  if (!is.null(color_var)) {
    p <- p + ggplot2::geom_point(ggplot2::aes(color = .data[[color_var]]), alpha = 0.5, size = 1.5)
  } else {
    p <- p + ggplot2::geom_point(alpha = 0.5, size = 1.5)
  }
  
  # Add smooth or linear fit
  if (use_smooth) {
    if (smoother == "gam") {
      k_use <- min(6, n_unique - 1)
      p <- p + ggplot2::geom_smooth(
        method = "gam", 
        formula = as.formula(paste("y ~ s(x, k =", k_use, ")")), 
        se = TRUE, color = "blue"
      )
    } else {
      p <- p + ggplot2::geom_smooth(method = "loess", se = TRUE, color = "blue")
    }
  } else {
    # Linear fit for low-variance metrics
    p <- p + ggplot2::geom_smooth(method = "lm", se = TRUE, color = "blue")
  }
  
  p <- p +
    ggplot2::annotate("label", x = Inf, y = Inf, label = r2_lab,
                      vjust = 1.2, hjust = 1.05, size = 4) +
    ggplot2::labs(
      title = paste0(metric_name, " vs ", auprc_var),
      x = metric_name,
      y = auprc_var
    ) +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(
      panel.background = ggplot2::element_rect(fill = "white", color = NA),
      plot.background = ggplot2::element_rect(fill = "white", color = NA)
    )
  
  p
}

#' Generate all X-metric vs AUPRC plots
#'
#' @param df_joined Data frame with X metrics and AUPRC columns
#' @param metric_names Character vector of metric columns to plot
#' @param auprc_var AUPRC column name
#' @param out_dir Output directory for plots
#' @param color_var Optional column for coloring points
#' @param smoother "gam" or "loess"
#' @return Invisible tibble of R² values per metric
#' @export
generate_x_metric_plots <- function(df_joined,
                                     metric_names,
                                     auprc_var = "AUPRC",
                                     out_dir = "plots_x_metrics",
                                     color_var = NULL,
                                     smoother = "gam") {
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  
  r2_results <- purrr::map_dfr(metric_names, function(m) {
    # Skip metrics with zero variance (like mu = 1 for all)
    vals <- df_joined[[m]]
    vals <- vals[!is.na(vals)]
    if (length(unique(vals)) < 2) {
      message("Skipping ", m, " - no variance (all values identical)")
      return(tibble::tibble(metric = m, r2 = NA_real_, skipped = TRUE))
    }
    
    p <- tryCatch({
      plot_x_metric_vs_auprc(df_joined, m, auprc_var, color_var, smoother)
    }, error = function(e) {
      warning("Error plotting ", m, ": ", e$message)
      NULL
    })
    
    if (is.null(p)) {
      return(tibble::tibble(metric = m, r2 = NA_real_, skipped = TRUE))
    }
    
    # Save plot
    fname <- file.path(out_dir, paste0(auprc_var, "_vs_", m, ".png"))
    ggplot2::ggsave(fname, p, width = 7, height = 5, dpi = 150, bg = "white")
    message("Saved: ", fname)
    
    # Recompute R² for return value (with error handling)
    dfm <- df_joined %>%
      dplyr::filter(!is.na(.data[[m]]), !is.na(.data[[auprc_var]]))
    
    r2 <- tryCatch({
      n_unique <- length(unique(dfm[[m]]))
      
      if (n_unique >= 6 && smoother == "gam") {
        k_use <- min(6, n_unique - 1)
        fit <- mgcv::gam(
          stats::as.formula(paste(auprc_var, "~ s(", m, ", k =", k_use, ")")),
          data = dfm, method = "REML"
        )
        yhat <- stats::predict(fit, newdata = dfm)
      } else if (n_unique >= 6 && smoother == "loess") {
        fit <- stats::loess(
          stats::as.formula(paste(auprc_var, "~", m)),
          data = dfm, span = 0.75
        )
        yhat <- stats::predict(fit)
      } else {
        # Linear model fallback
        fit <- stats::lm(stats::as.formula(paste(auprc_var, "~", m)), data = dfm)
        yhat <- stats::predict(fit)
      }
      
      1 - sum((dfm[[auprc_var]] - yhat)^2, na.rm = TRUE) / 
          sum((dfm[[auprc_var]] - mean(dfm[[auprc_var]], na.rm = TRUE))^2, na.rm = TRUE)
    }, error = function(e) NA_real_)
    
    tibble::tibble(metric = m, r2 = r2, skipped = FALSE)
  })
  
  # Print summary sorted by R²
  message("\n=== R² Summary (sorted) ===")
  r2_sorted <- r2_results %>% dplyr::arrange(dplyr::desc(r2))
  print(r2_sorted, n = Inf)
  
  invisible(r2_sorted)
}

#' Get default list of X-metric names
#' @return Character vector of metric column names
#' @export
get_x_metric_names <- function() {
  c(
    # Axis A: Basin size / resolution difficulty
    "ecs", "high_ld_mass", "top5_ev_frac",
    "n_blocks_tau50", "avg_block_size_tau50", "max_block_size_tau50",
    "max_block_frac_tau50", "block_entropy_norm_tau50",
    "n_blocks_tau80", "avg_block_size_tau80", "max_block_size_tau80",
    "max_block_frac_tau80", "block_entropy_norm_tau80",
    # Axis B: Identifiability / blur
    "nsi", "M1", "M2", "M1_dist",
    # Spectral metrics (mixed)
    "lambda_max_frac", "participation_ratio", "spectral_entropy_norm",
    "stable_rank", "spectral_gap",
    # Simple summaries
    "mean_abs_r", "mean_r2", "row_conc_mean", "row_conc_q90",
    # Structure metrics
    "ar1_rho_hat", "toeplitz_rel_F_error",
    # Submodularity/coherence (Axis B theoretical)
    "mu", "lambda_min_2k_mc", "gamma_lb_mc", "gamma_lb_gersh"
  )
}

#' Get Axis A metric names (Basin Size / Resolution Difficulty)
#' 
#' Returns metrics that measure how big/strong LD neighborhoods are.
#' High Axis A = SuSiE finds right basin but CS is large.
#' 
#' @return Character vector of Axis A metric names
#' @export
get_axis_a_metrics <- function() {
  c(
    "ecs", "high_ld_mass", "top5_ev_frac",
    "n_blocks_tau50", "avg_block_size_tau50", "max_block_size_tau50",
    "max_block_frac_tau50", "block_entropy_norm_tau50",
    "n_blocks_tau80", "avg_block_size_tau80", "max_block_size_tau80",
    "max_block_frac_tau80", "block_entropy_norm_tau80",
    # These also relate to block structure
    "row_conc_mean", "row_conc_q90", "lambda_max_frac"
  )
}

#' Get Axis B metric names (Identifiability / Factorization Mismatch)
#' 
#' Returns metrics that measure how interchangeable SNP neighborhoods are.
#' High Axis B = many SNPs look alike, SuSiE's factorization may be wrong shape.
#' 
#' @return Character vector of Axis B metric names
#' @export
get_axis_b_metrics <- function() {
  c(
    "nsi", "M1", "M2", "M1_dist",
    "spectral_entropy_norm", "participation_ratio",
    "lambda_min_2k_mc", "gamma_lb_mc",
    "toeplitz_rel_F_error", "mean_r2"
  )
}

#' Validate computed X-matrix metrics
#' 
#' Checks that metrics are in expected ranges and reports summary statistics.
#' Use this to verify the computation worked correctly.
#'
#' @param x_metrics Tibble of X-matrix metrics from build_x_metrics_table
#' @return Invisible TRUE if all checks pass, prints diagnostic info
#' @export
validate_x_metrics <- function(x_metrics) {
  cat("\n=== X-Matrix Metrics Validation ===\n\n")
  
  # Check for missing files
  n_missing <- sum(x_metrics$.missing_file, na.rm = TRUE)
  n_total <- nrow(x_metrics)
  cat(sprintf("Datasets: %d total, %d missing (%.1f%%)\n", 
              n_total, n_missing, 100 * n_missing / n_total))
  
  if (n_missing == n_total) {
    warning("All files missing! Check paths.")
    return(invisible(FALSE))
  }
  
  # Filter to valid rows
  valid <- x_metrics %>% dplyr::filter(!.missing_file | is.na(.missing_file))
  cat(sprintf("Valid datasets: %d\n\n", nrow(valid)))
  
  # Define expected ranges for key metrics
  checks <- list(
    M1 = c(0, 0.5),           # Mid-energy: 0 to 0.5 (max at r=0.5)
    M2 = c(0, 0.25),          # Quadratic mid-energy: 0 to 0.25
    mu = c(0, 1),             # Mutual coherence: 0 to 1
    lambda_max_frac = c(0, 1), # Largest eigenvalue fraction
    participation_ratio = c(0, 1),
    mean_abs_r = c(0, 1),
    ar1_rho_hat = c(-1, 1)    # AR(1) coefficient
  )
  
  cat("Metric ranges (should be within expected bounds):\n")
  all_ok <- TRUE
  for (metric in names(checks)) {
    if (!metric %in% names(valid)) next
    
    vals <- valid[[metric]]
    vals <- vals[!is.na(vals)]
    if (length(vals) == 0) next
    
    lo <- checks[[metric]][1]
    hi <- checks[[metric]][2]
    actual_lo <- min(vals)
    actual_hi <- max(vals)
    
    in_range <- actual_lo >= lo - 0.01 && actual_hi <= hi + 0.01
    status <- if (in_range) "OK" else "WARNING"
    if (!in_range) all_ok <- FALSE
    
    cat(sprintf("  %-20s: [%6.3f, %6.3f]  expected [%.1f, %.1f]  %s\n",
                metric, actual_lo, actual_hi, lo, hi, status))
  }
  
  cat("\nSummary statistics for key metrics:\n")
  summary_metrics <- c("M1", "mu", "lambda_max_frac", "mean_abs_r", "ar1_rho_hat")
  for (metric in summary_metrics) {
    if (!metric %in% names(valid)) next
    vals <- valid[[metric]]
    cat(sprintf("  %-20s: mean=%.3f, sd=%.3f, median=%.3f\n",
                metric, mean(vals, na.rm = TRUE), sd(vals, na.rm = TRUE), 
                median(vals, na.rm = TRUE)))
  }
  
  # Check for suspicious patterns
  cat("\nDiagnostics:\n")
  
  # All identical values (suggests computation issue)
  for (metric in c("M1", "M2", "mu")) {
    if (!metric %in% names(valid)) next
    vals <- valid[[metric]]
    if (length(unique(vals[!is.na(vals)])) == 1) {
      cat(sprintf("  WARNING: %s has identical values for all datasets\n", metric))
      all_ok <- FALSE
    }
  }
  
  # High NA rate
  for (metric in get_x_metric_names()) {
    if (!metric %in% names(valid)) next
    na_rate <- mean(is.na(valid[[metric]]))
    if (na_rate > 0.1) {
      cat(sprintf("  WARNING: %s has %.1f%% NA values\n", metric, 100 * na_rate))
    }
  }
  
  if (all_ok) {
    cat("\n✓ All validation checks passed!\n")
  } else {
    cat("\n⚠ Some checks had warnings - review above\n")
  }
  
  invisible(all_ok)
}

#' Test all 1-way and 2-way GAM models for X-metrics vs AUPRC
#'
#' Systematically fits GAMs with all single metrics and all pairs of metrics
#' to identify which metrics (and interactions) best predict AUPRC.
#'
#' @param df_joined Data frame with X metrics and AUPRC columns
#' @param metric_names Character vector of metric columns to test
#' @param auprc_var AUPRC column name (can be raw, residualized, or ranked)
#' @param out_path Optional path to save CSV results
#' @param k_smooth Basis dimension for smooth terms (default 5)
#' @param cross_axis_only If TRUE, only test 2-way models where one metric is
#'   from Axis A (basin size) and one from Axis B (identifiability). This is
#'   much more efficient and theoretically motivated than all pairs.
#' @return Tibble with model diagnostics for all combinations
#' @export
test_gam_combinations <- function(df_joined,
                                   metric_names = NULL,
                                   auprc_var = "AUPRC",
                                   out_path = NULL,
                                   k_smooth = 5,
                                   cross_axis_only = FALSE) {
  
  if (is.null(metric_names)) {
    metric_names <- get_x_metric_names()
  }
  
  # Filter to metrics present in data with sufficient variance
  metric_names <- metric_names[metric_names %in% names(df_joined)]
  
  valid_metrics <- sapply(metric_names, function(m) {
    vals <- df_joined[[m]]
    vals <- vals[!is.na(vals)]
    length(unique(vals)) >= k_smooth
  })
  metric_names <- metric_names[valid_metrics]
  
  # For cross-axis mode, determine which metrics belong to which axis
  if (cross_axis_only) {
    axis_a <- get_axis_a_metrics()
    axis_b <- get_axis_b_metrics()
    metrics_a <- intersect(metric_names, axis_a)
    metrics_b <- intersect(metric_names, axis_b)
    
    # Generate only cross-axis pairs
    pairs <- expand.grid(m1 = metrics_a, m2 = metrics_b, stringsAsFactors = FALSE)
    pairs <- lapply(seq_len(nrow(pairs)), function(i) c(pairs$m1[i], pairs$m2[i]))
    
    message("Testing ", length(metric_names), " metrics with sufficient variance")
    message("  Axis A metrics: ", length(metrics_a))
    message("  Axis B metrics: ", length(metrics_b))
    message("Fitting 1-way GAMs: ", length(metric_names), " models")
    message("Fitting 2-way GAMs (cross-axis only): ", length(pairs), " pairs")
  } else {
    pairs <- utils::combn(metric_names, 2, simplify = FALSE)
    
    message("Testing ", length(metric_names), " metrics with sufficient variance")
    message("Fitting 1-way GAMs: ", length(metric_names), " models")
    message("Fitting 2-way GAMs (all pairs): ", length(pairs), " pairs")
  }
  
  # Prepare clean data
  df_clean <- df_joined %>%
    dplyr::select(dplyr::all_of(c(auprc_var, metric_names))) %>%
    tidyr::drop_na()
  
  n_obs <- nrow(df_clean)
  message("Using ", n_obs, " complete observations")
  
  if (n_obs < 20) {
    warning("Too few observations for reliable GAM fitting")
    return(NULL)
  }
  
  results <- list()
  
  # ===== 1-WAY MODELS =====
  message("\n--- Fitting 1-way GAMs ---")
  
  for (m1 in metric_names) {
    res <- fit_gam_safe(df_clean, auprc_var, m1, NULL, k_smooth)
    if (!is.null(res)) {
      results[[length(results) + 1]] <- res
    }
  }
  
  # ===== 2-WAY MODELS (additive: s(m1) + s(m2)) =====
  if (cross_axis_only) {
    message("\n--- Fitting 2-way additive GAMs (cross-axis) ---")
  } else {
    message("\n--- Fitting 2-way additive GAMs ---")
  }
  
  for (pair in pairs) {
    m1 <- pair[1]
    m2 <- pair[2]
    res <- fit_gam_safe(df_clean, auprc_var, m1, m2, k_smooth, interaction = FALSE)
    if (!is.null(res)) {
      results[[length(results) + 1]] <- res
    }
  }
  
  # ===== 2-WAY MODELS WITH INTERACTION: s(m1, m2) =====
  if (cross_axis_only) {
    message("\n--- Fitting 2-way interaction GAMs (cross-axis) ---")
  } else {
    message("\n--- Fitting 2-way interaction GAMs ---")
  }
  
  for (pair in pairs) {
    m1 <- pair[1]
    m2 <- pair[2]
    res <- fit_gam_safe(df_clean, auprc_var, m1, m2, k_smooth, interaction = TRUE)
    if (!is.null(res)) {
      results[[length(results) + 1]] <- res
    }
  }
  
  # Combine results
  results_df <- dplyr::bind_rows(results) %>%
    dplyr::arrange(dplyr::desc(r2_adj))
  
  # Print top models
  message("\n=== Top 20 Models by Adjusted R² ===")
  print(results_df %>% head(20), n = 20)
  
  # Save if path provided
 if (!is.null(out_path)) {
    readr::write_csv(results_df, out_path)
    message("\nSaved results to: ", out_path)
  }
  
  invisible(results_df)
}

#' Safely fit a GAM and extract diagnostics
#' @keywords internal
fit_gam_safe <- function(df, y_var, m1, m2 = NULL, k = 5, interaction = FALSE) {
  
  tryCatch({
    # Build formula
    if (is.null(m2)) {
      # 1-way model
      formula_str <- sprintf("%s ~ s(%s, k = %d)", y_var, m1, k)
      model_type <- "1-way"
      metric1 <- m1
      metric2 <- NA_character_
    } else if (!interaction) {
      # 2-way additive
      formula_str <- sprintf("%s ~ s(%s, k = %d) + s(%s, k = %d)", y_var, m1, k, m2, k)
      model_type <- "2-way_additive"
      metric1 <- m1
      metric2 <- m2
    } else {
      # 2-way with tensor product interaction
      formula_str <- sprintf("%s ~ te(%s, %s, k = %d)", y_var, m1, m2, k)
      model_type <- "2-way_interaction"
      metric1 <- m1
      metric2 <- m2
    }
    
    # Fit GAM
    fit <- mgcv::gam(
      stats::as.formula(formula_str),
      data = df,
      method = "REML"
    )
    
    # Extract diagnostics
    summ <- summary(fit)
    
    # Deviance explained
    dev_expl <- summ$dev.expl
    
    # R² and adjusted R²
    r2 <- summ$r.sq
    n <- nrow(df)
    p_eff <- sum(fit$edf)  # effective degrees of freedom
    r2_adj <- 1 - (1 - r2) * (n - 1) / (n - p_eff - 1)
    
    # AIC/BIC
    aic <- stats::AIC(fit)
    bic <- stats::BIC(fit)
    
    # GCV score
    gcv <- fit$gcv.ubre
    
    # Smooth term p-values and edf
    smooth_table <- summ$s.table
    
    if (is.null(m2)) {
      # 1-way
      edf1 <- smooth_table[1, "edf"]
      pval1 <- smooth_table[1, "p-value"]
      edf2 <- NA_real_
      pval2 <- NA_real_
    } else if (!interaction) {
      # 2-way additive
      edf1 <- smooth_table[1, "edf"]
      pval1 <- smooth_table[1, "p-value"]
      edf2 <- smooth_table[2, "edf"]
      pval2 <- smooth_table[2, "p-value"]
    } else {
      # 2-way interaction (tensor)
      edf1 <- smooth_table[1, "edf"]
      pval1 <- smooth_table[1, "p-value"]
      edf2 <- NA_real_  # tensor combines them
      pval2 <- NA_real_
    }
    
    # RMSE
    rmse <- sqrt(mean(fit$residuals^2))
    
    tibble::tibble(
      model_type = model_type,
      metric1 = metric1,
      metric2 = metric2,
      formula = formula_str,
      n_obs = n,
      r2 = r2,
      r2_adj = r2_adj,
      dev_expl = dev_expl,
      aic = aic,
      bic = bic,
      gcv = gcv,
      rmse = rmse,
      edf_total = p_eff,
      edf1 = edf1,
      pval1 = pval1,
      edf2 = edf2,
      pval2 = pval2
    )
    
  }, error = function(e) {
    # Return NULL on error (will be filtered out)
    warning("Failed to fit: ", m1, if (!is.null(m2)) paste0(" + ", m2), " - ", e$message)
    NULL
  })
}

#' Compare 1-way vs 2-way models to find meaningful interactions
#'
#' For each pair of metrics, compares the 2-way interaction model to
#' the sum of individual 1-way models to identify synergistic effects.
#'
#' @param gam_results Output from test_gam_combinations
#' @return Tibble with interaction strength for each pair
#' @export
find_metric_interactions <- function(gam_results) {
  
  # Get 1-way results
  one_way <- gam_results %>%
    dplyr::filter(model_type == "1-way") %>%
    dplyr::select(metric = metric1, r2_1way = r2_adj)
  
  # Get 2-way additive results
  two_way_add <- gam_results %>%
    dplyr::filter(model_type == "2-way_additive") %>%
    dplyr::select(metric1, metric2, r2_additive = r2_adj, aic_additive = aic)
  
  # Get 2-way interaction results
  two_way_int <- gam_results %>%
    dplyr::filter(model_type == "2-way_interaction") %>%
    dplyr::select(metric1, metric2, r2_interaction = r2_adj, aic_interaction = aic)
  
  # Join 1-way R² for each metric in the pair
  interactions <- two_way_add %>%
    dplyr::left_join(two_way_int, by = c("metric1", "metric2")) %>%
    dplyr::left_join(one_way %>% dplyr::rename(metric1 = metric, r2_m1 = r2_1way), by = "metric1") %>%
    dplyr::left_join(one_way %>% dplyr::rename(metric2 = metric, r2_m2 = r2_1way), by = "metric2") %>%
    dplyr::mutate(
      # Expected R² if metrics were independent
      r2_expected = r2_m1 + r2_m2 - r2_m1 * r2_m2,
      
      # Synergy: how much better is additive model than expected?
      synergy_additive = r2_additive - pmax(r2_m1, r2_m2),
      
      # Interaction effect: how much better is interaction vs additive?
      interaction_effect = r2_interaction - r2_additive,
      
      # AIC improvement from interaction
      aic_improvement = aic_additive - aic_interaction
    ) %>%
    dplyr::arrange(dplyr::desc(interaction_effect))
  
  message("\n=== Metric Pairs with Strongest Interactions ===")
  message("(interaction_effect > 0 means tensor product helps over additive)\n")
  
  top_int <- interactions %>%
    dplyr::filter(interaction_effect > 0.01) %>%
    dplyr::select(metric1, metric2, r2_m1, r2_m2, r2_additive, r2_interaction, 
                  interaction_effect, aic_improvement)
  
  if (nrow(top_int) > 0) {
    print(top_int, n = 20)
  } else {
    message("No strong interactions found (all interaction_effect < 0.01)")
  }
  
  invisible(interactions)
}
