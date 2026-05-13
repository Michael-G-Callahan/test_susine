test_that("susieR normalization preserves raw and standardized coefficient scales", {
  skip_if_not_installed("susieR")

  set.seed(1101)
  n <- 120
  p <- 12
  X <- matrix(stats::rbinom(n * p, size = 2, prob = 0.18), nrow = n, ncol = p)
  storage.mode(X) <- "double"
  y <- 0.8 * X[, 3] - 0.4 * X[, 8] + stats::rnorm(n, sd = 0.8)

  fit_raw <- susieR::susie(
    X = X,
    y = y,
    L = 5,
    scaled_prior_variance = 0.2,
    estimate_prior_variance = FALSE,
    estimate_residual_variance = TRUE,
    max_iter = 100,
    tol = 1e-4,
    verbose = FALSE
  )

  fit <- test_susine:::normalize_susier_fit(fit_raw, X = X, L = 5)

  expect_equal(fit$model_fit$fitted_y, as.numeric(fit_raw$fitted), tolerance = 1e-8)
  expect_equal(fit$model_fit$coef, as.numeric(stats::coef(fit_raw)[-1]), tolerance = 1e-8)
  expect_equal(
    fit$model_fit$std_coef,
    as.numeric(colSums(fit_raw$alpha * fit_raw$mu)),
    tolerance = 1e-8
  )
})
