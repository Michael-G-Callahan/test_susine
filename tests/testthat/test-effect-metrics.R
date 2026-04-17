test_that("evaluate_model reports effect diffuseness and accuracy metrics", {
  X <- matrix(
    c(
      1, 0, 0, 1,
      0, 1, 1, 0,
      1, 1, 0, 0,
      0, 0, 1, 1,
      1, 0, 1, 0,
      0, 1, 0, 1
    ),
    ncol = 4,
    byrow = TRUE
  )
  X <- scale(X, center = TRUE, scale = TRUE)

  alpha_1 <- c(0.8, 0.1, 0.05, 0.05)
  alpha_2 <- c(0.4, 0.3, 0.2, 0.1)
  alpha <- rbind(alpha_1, alpha_2)

  fit <- list(
    settings = list(L = 2L),
    effect_fits = list(alpha = alpha),
    model_fit = list(
      PIPs = 1 - apply(1 - alpha, 2, prod),
      fitted_y = rep(0, nrow(X)),
      sigma_2 = 1
    )
  )

  res <- test_susine:::evaluate_model(
    fit = fit,
    X = X,
    y = rep(0, nrow(X)),
    causal_idx = 2L,
    purity_threshold = -Inf,
    compute_curves = FALSE
  )

  eff <- res$effects_unfiltered

  expect_true(all(c(
    "effect_pip_entropy",
    "effect_pip_entropy_core95",
    "effect_k_eff_signal",
    "effect_k_eff_signal_core95",
    "tail_inflation_ratio",
    "tail_inflation_log",
    "accuracy_ratio"
  ) %in% names(eff)))
  expect_equal(nrow(res$effects_filtered), 2L)

  expect_equal(eff$size, c(3, 4))
  expect_equal(eff$coverage, c(1L, 1L))
  expect_equal(eff$accuracy_ratio, c(0.1 / 0.8, 0.3 / 0.4), tolerance = 1e-12)

  entropy_1 <- -sum(alpha_1 * log(alpha_1))
  core_1 <- alpha_1[c(1, 2, 3)] / sum(alpha_1[c(1, 2, 3)])
  entropy_core_1 <- -sum(core_1 * log(core_1))

  entropy_2 <- -sum(alpha_2 * log(alpha_2))

  expect_equal(eff$effect_pip_entropy[[1]], entropy_1, tolerance = 1e-12)
  expect_equal(eff$effect_pip_entropy_core95[[1]], entropy_core_1, tolerance = 1e-12)
  expect_equal(eff$effect_k_eff_signal[[1]], exp(entropy_1), tolerance = 1e-12)
  expect_equal(
    eff$effect_k_eff_signal_core95[[1]],
    exp(entropy_core_1),
    tolerance = 1e-12
  )
  expect_equal(
    eff$tail_inflation_ratio[[1]],
    exp(entropy_1 - entropy_core_1),
    tolerance = 1e-12
  )
  expect_equal(
    eff$tail_inflation_log[[1]],
    entropy_1 - entropy_core_1,
    tolerance = 1e-12
  )

  expect_equal(eff$effect_pip_entropy[[2]], entropy_2, tolerance = 1e-12)
  expect_equal(eff$effect_pip_entropy_core95[[2]], entropy_2, tolerance = 1e-12)
  expect_equal(eff$effect_k_eff_signal[[2]], exp(entropy_2), tolerance = 1e-12)
  expect_equal(eff$effect_k_eff_signal_core95[[2]], exp(entropy_2), tolerance = 1e-12)
  expect_equal(eff$tail_inflation_ratio[[2]], 1, tolerance = 1e-12)
  expect_equal(eff$tail_inflation_log[[2]], 0, tolerance = 1e-12)
})
