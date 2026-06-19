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
test_that("purity_filter_confusion_sweep recomputes PIPs after dropping diffuse effects", {
  alpha <- rbind(
    c(0.90, 0.02, 0.02, 0.03, 0.03),
    c(0.01, 0.55, 0.14, 0.15, 0.15),
    c(0.01, 0.01, 0.01, 0.85, 0.12)
  )

  effects <- tibble::tibble(
    purity = c(0.8, 0.8, 0.8),
    effect_k_eff_signal_core95 = c(1.4, 5.0, 1.5)
  )
  causal_vec <- c(1L, 1L, 0L, 0L, 0L)

  bins <- test_susine:::purity_filter_confusion_sweep(
    alpha = alpha,
    effects_unfiltered = effects,
    causal_vec = causal_vec,
    causal_mask = c(1L, 2L),
    pip_breaks = seq(0, 1, by = 0.001),
    purity_thresholds = 0.5,
    keff_core95_max = 3,
    pip_bucket_width = 0.001
  )

  expect_equal(
    unique(bins$n_effects_kept[bins$filter_type == "keff_core95_le_3"]),
    2
  )

  auprc <- test_susine::compute_auprc_from_confusion(
    bins,
    group_vars = "filter_type"
  )
  no_filter <- auprc$AUPRC[auprc$filter_type == "no_filter"]
  keff_filter <- auprc$AUPRC[auprc$filter_type == "keff_core95_le_3"]

  expect_lt(keff_filter, no_filter)

  pip_unfiltered <- test_susine:::combined_pip_from_alpha(alpha)
  pip_filtered <- test_susine:::combined_pip_from_alpha(alpha[c(1, 3), , drop = FALSE])
  expect_lt(pip_filtered[[2]], pip_unfiltered[[2]])
})
