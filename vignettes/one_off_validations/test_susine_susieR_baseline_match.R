# =============================================================================
# test_susine_susieR_baseline_match.R
#
# Informal unit test for fast regression checks after susine edits.
#
# Goal:
# - Generate a simple synthetic X matrix and y vector.
# - Fit vanilla SuSiE in both engines:
#   * susieR::susie
#   * susine::susine
# - Use the same fixed prior variance (0.2 * var(y)).
# - Keep ALL prior EB updates off.
# - Keep residual variance estimation on.
#
# Exit code:
# - 0 if all checks pass
# - 1 if one or more checks fail
#
# Usage from Code/ root:
#   Rscript test_susine/vignettes/one_off_validations/test_susine_susieR_baseline_match.R
# Usage from test_susine/:
#   Rscript vignettes/one_off_validations/test_susine_susieR_baseline_match.R
# =============================================================================

# --- Setup -------------------------------------------------------------------

try_load <- function(pkg, path) {
  # Prefer local repo code when available so this checks in-worktree edits.
  if (file.exists(file.path(path, "DESCRIPTION"))) {
    if (!requireNamespace("devtools", quietly = TRUE)) {
      stop(sprintf(
        "Local path for '%s' exists (%s) but 'devtools' is unavailable.",
        pkg, path
      ))
    }
    message(sprintf("Loading local %s via devtools::load_all('%s')", pkg, path))
    devtools::load_all(path, quiet = TRUE)
    return(invisible(TRUE))
  }
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop(sprintf("Package '%s' is not installed and local path was not found: %s", pkg, path))
  }
  invisible(TRUE)
}

wd <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
if (grepl("/test_susine$", wd)) {
  susine_path <- "../susine"
  susieR_path <- "../susieR"
} else {
  susine_path <- "susine"
  susieR_path <- "susieR"
}

try_load("susine", susine_path)
try_load("susieR", susieR_path)

# --- Synthetic data -----------------------------------------------------------

set.seed(20260227)

n <- 300L
p <- 150L
L <- 10L
scaled_prior_variance <- 0.2
tol <- 1e-5
max_iter <- 200L

# Create a basic correlated design matrix to mimic LD structure.
X <- matrix(rnorm(n * p), nrow = n, ncol = p)
for (j in 2:p) {
  X[, j] <- 0.7 * X[, j - 1] + sqrt(1 - 0.7^2) * X[, j]
}

causal_idx <- sort(sample(seq_len(p), 4))
beta_true <- rep(0, p)
beta_true[causal_idx] <- rnorm(length(causal_idx), mean = 0, sd = 0.8)

signal <- drop(X %*% beta_true)
target_pve <- 0.35
noise_var <- var(signal) * (1 - target_pve) / target_pve
y <- signal + rnorm(n, sd = sqrt(noise_var))

cat("=== susine vs susieR vanilla baseline check ===\n")
cat(sprintf(
  "Data: n=%d, p=%d, L=%d, causal=%s, target_pve=%.2f\n",
  n, p, L, paste(causal_idx, collapse = ","), target_pve
))
cat(sprintf(
  "Settings: scaled_prior_variance=%.3f, estimate_prior_variance=FALSE, prior_update_method='none'\n",
  scaled_prior_variance
))

# --- Fits --------------------------------------------------------------------

susie_formals <- names(formals(susieR::susie))
susie_args <- list(
  X = X,
  y = y,
  L = L,
  scaled_prior_variance = scaled_prior_variance,
  estimate_prior_variance = FALSE,
  estimate_residual_variance = TRUE,
  standardize = TRUE,
  intercept = TRUE,
  tol = tol,
  max_iter = max_iter
)
if ("convergence_method" %in% susie_formals) {
  susie_args$convergence_method <- "elbo"
}
if ("refine" %in% susie_formals) {
  susie_args$refine <- FALSE
}

fit_susieR <- do.call(susieR::susie, susie_args)

fit_susine <- susine::susine(
  L = L,
  X = X,
  y = y,
  mu_0 = 0,
  sigma_0_2 = scaled_prior_variance,
  prior_inclusion_weights = NULL,
  prior_update_method = "none",
  auto_scale_mu_0 = FALSE,
  auto_scale_sigma_0_2 = FALSE,
  tol = tol,
  max_iter = max_iter
)

pip_susieR <- fit_susieR$pip
pip_susine <- fit_susine$model_fit$PIPs
sigma2_susieR <- fit_susieR$sigma2
sigma2_susine <- tail(fit_susine$model_fit$sigma_2, 1)
elbo_susieR <- tail(fit_susieR$elbo, 1)
elbo_susine <- tail(fit_susine$model_fit$elbo, 1)

# --- Checks ------------------------------------------------------------------

n_checks <- 0L
n_pass <- 0L

check <- function(name, condition, detail = "") {
  n_checks <<- n_checks + 1L
  if (isTRUE(condition)) {
    n_pass <<- n_pass + 1L
    cat(sprintf("[PASS] %s%s\n", name, if (nzchar(detail)) paste0(" - ", detail) else ""))
  } else {
    cat(sprintf("[FAIL] %s%s\n", name, if (nzchar(detail)) paste0(" - ", detail) else ""))
  }
}

pip_abs_diff <- abs(pip_susine - pip_susieR)
max_pip_diff <- max(pip_abs_diff)
mean_pip_diff <- mean(pip_abs_diff)
pip_cor <- cor(pip_susine, pip_susieR)
top10_overlap <- length(intersect(order(pip_susine, decreasing = TRUE)[1:10],
                                  order(pip_susieR, decreasing = TRUE)[1:10]))
sigma2_rel_diff <- abs(sigma2_susine - sigma2_susieR) / sigma2_susieR
elbo_rel_diff <- abs(elbo_susine - elbo_susieR) / abs(elbo_susieR)

cat("\n--- Metrics ---\n")
cat(sprintf("Max |PIP diff|: %.6f\n", max_pip_diff))
cat(sprintf("Mean |PIP diff|: %.6f\n", mean_pip_diff))
cat(sprintf("PIP correlation: %.8f\n", pip_cor))
cat(sprintf("Top-10 PIP overlap: %d/10\n", top10_overlap))
cat(sprintf("Relative sigma2 diff: %.6f\n", sigma2_rel_diff))
cat(sprintf("Relative ELBO diff: %.6f\n", elbo_rel_diff))

check("max PIP diff < 0.02", max_pip_diff < 0.02, sprintf("%.6f", max_pip_diff))
check("mean PIP diff < 0.003", mean_pip_diff < 0.003, sprintf("%.6f", mean_pip_diff))
check("PIP correlation > 0.995", pip_cor > 0.995, sprintf("%.8f", pip_cor))
check("top-10 overlap >= 7", top10_overlap >= 7, sprintf("%d/10", top10_overlap))
check("relative sigma2 diff < 0.05", sigma2_rel_diff < 0.05, sprintf("%.6f", sigma2_rel_diff))
check("relative ELBO diff < 0.01", elbo_rel_diff < 0.01, sprintf("%.6f", elbo_rel_diff))

cat(sprintf("\n=== %d/%d checks passed ===\n", n_pass, n_checks))

if (n_pass == n_checks) {
  cat("ALL CHECKS PASSED.\n")
} else {
  cat("CHECKS FAILED. Review metrics above.\n")
}

if (!interactive()) {
  quit(status = if (n_pass == n_checks) 0L else 1L, save = "no")
}
