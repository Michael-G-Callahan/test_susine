# Tests for the local-genetic-variance ("heritability") estimand helpers and the
# corrected hg2_by_agg decomposition. See
# refs/decisions/heritability_estimand_decision_2026-06-09.md.

make_fake_fit <- function(alpha, mu, mu2) {
  list(effect_fits = list(
    b_hat = alpha * mu,
    b_2_hat = alpha * mu2
  ))
}

test_that("hg2_components: standardizes X (susine scale), not the raw design", {
  set.seed(1)
  n <- 200; p <- 8; L <- 3
  # Raw design with very different per-column means/scales so raw vs
  # standardized give materially different quadratic forms (regression guard
  # for the X-scaling bug: susine fits coefficients on the standardized scale).
  X <- matrix(rnorm(n * p), n, p)
  X <- sweep(X, 2L, runif(p, 0.5, 5), "*")
  X <- sweep(X, 2L, rnorm(p, 0, 3), "+")
  alpha <- matrix(0, L, p)
  for (l in 1:L) {
    w <- runif(p); w[sample(p, p - 3)] <- 0; alpha[l, ] <- w / sum(w)
  }
  mu <- matrix(rnorm(L * p, 0, 0.5), L, p)
  fit <- make_fake_fit(alpha, mu, mu^2 + 0.3)
  m <- colSums(fit$effect_fits$b_hat)
  y <- as.numeric(scale(X) %*% m) + rnorm(n)

  ci <- hg2_components(fit, X = X, y = y)
  expect_true(is.finite(ci$expected_pve))
  expect_equal(ci$postmean + ci$uncertainty, ci$expected_pve, tolerance = 1e-10)
  expect_gte(ci$uncertainty, 0)
  # postmean must use the STANDARDIZED design (matches var(fitted_y)) ...
  fy_std <- as.numeric(scale(X) %*% m)
  expect_equal(ci$postmean, stats::var(fy_std) / stats::var(y), tolerance = 1e-9)
  # ... and must NOT equal the (wrong) raw-design quadratic form.
  fy_raw <- as.numeric(X %*% m)
  expect_false(isTRUE(all.equal(ci$postmean, stats::var(fy_raw) / stats::var(y))))
})

test_that("hg2_components: RSS(R) path matches the standardized individual path", {
  set.seed(2)
  n <- 200; p <- 8; L <- 3
  X <- matrix(rnorm(n * p), n, p)
  X <- sweep(X, 2L, runif(p, 0.5, 3), "*")
  alpha <- matrix(0, L, p)
  for (l in 1:L) {
    w <- runif(p); w[sample(p, p - 3)] <- 0; alpha[l, ] <- w / sum(w)
  }
  mu <- matrix(rnorm(L * p, 0, 0.5), L, p)
  fit <- make_fake_fit(alpha, mu, mu^2 + 0.3)
  y <- as.numeric(scale(X) %*% colSums(fit$effect_fits$b_hat)) + rnorm(n)

  R <- crossprod(scale(X)) / (n - 1)   # standardized Gram == helper's convention
  vy <- stats::var(y)
  rss <- hg2_components(fit, R = R, vy = vy)
  ind <- hg2_components(fit, X = X, y = y, vy = vy)  # raw X; helper standardizes
  expect_equal(rss$postmean, ind$postmean, tolerance = 1e-8)
  expect_equal(rss$uncertainty, ind$uncertainty, tolerance = 1e-8)
})

test_that("hg2_components: NA-safe when per-effect second moments absent", {
  fit0 <- list(effect_fits = list(alpha = matrix(0.5, 2, 4)))
  ci <- hg2_components(fit0, X = matrix(rnorm(40), 10, 4), y = rnorm(10))
  expect_true(is.na(ci$expected_pve))
})

test_that("compute_hg2_by_agg: corrected columns and additive identity", {
  set.seed(3)
  n <- 250; p <- 12; K <- 6
  X <- matrix(rnorm(n * p), n, p)
  beta <- rep(0, p); beta[c(2, 7)] <- c(0.5, -0.4)
  y <- as.numeric(X %*% beta) + rnorm(n)

  fy_list <- lapply(seq_len(K), function(k) as.numeric(X %*% (beta + rnorm(p, 0, 0.15))))
  pip_list <- lapply(seq_len(K), function(k) runif(p))
  elbo_vec <- sort(rnorm(K), decreasing = TRUE)
  unc_list <- as.list(runif(K, 0.01, 0.05))
  grr <- tibble::tibble(use_case_id = "uc", group_key = "g")

  tbl <- compute_hg2_by_agg(
    pip_list = pip_list, fitted_y_list = fy_list, elbo_vec = elbo_vec,
    y = y, group_run_row = grr, dataset_bundle_id = "b",
    uncertainty_list = unc_list, pip_cache = NULL
  )
  need <- c("hg2", "hg2_weighted_pve", "hg2_postmean", "hg2_uncertainty",
            "hg2_between_fit", "hg2_expected_pve")
  expect_true(all(need %in% names(tbl)))
  expect_equal(
    tbl$hg2_expected_pve,
    tbl$hg2_postmean + tbl$hg2_between_fit + tbl$hg2_uncertainty,
    tolerance = 1e-9
  )
  # max-ELBO puts all weight on one fit -> no between-fit dispersion.
  me <- tbl[tbl$agg_method == "max_elbo", ]
  expect_equal(me$hg2_between_fit, 0, tolerance = 1e-9)
  # corrected estimand is never below the posterior-mean part.
  expect_true(all(tbl$hg2_expected_pve >= tbl$hg2_postmean - 1e-9))
})

test_that("softmax_temperature lives at job$compute and survives JSON round-trip", {
  # Regression for the call-site config-path bug: compute_hg2_by_agg() must read
  # the temperature from the same location make_job_config writes it. The sim
  # schema nests it under job$compute (run_controls.R compute_list); a flat read
  # of job$softmax_temperature would silently fall back to 1.
  skip_if_not_installed("jsonlite")
  catalog <- tibble::tibble(
    data_scenario = "simulation_n3",
    dataset_label = "simulation_n3_1",
    participant_count = NA_integer_, snps_post = NA_integer_,
    snp_set = NA_character_, matrix_path = NA_character_,
    manifest_path = NA_character_, source = "simulation",
    matrix_index = 1L, matrix_id = 1L
  )
  cfg <- test_susine::make_job_config(
    job_name = "unit_hg2_softmax_temp",
    use_case_ids = "susine_functional_mu",
    exploration_methods = "c_grid",
    exploration_mode = "separate",
    K = 1L, L_grid = 2L, y_noise_grid = 0.8,
    prior_quality = test_susine::prior_quality_grid(c(0.4), c(1)),
    p_star_grid = 2L, seeds = 1L, architecture_grid = "sparse",
    data_scenarios = "simulation_n3", repo_root = ".",
    data_matrix_catalog = catalog,
    task_unit = "dataset", bundles_per_task = 1L, runs_per_task = 10L,
    output_root = tempdir(), c_grid_values = 0.7,
    write_snps_parquet = FALSE, write_confusion_bins = FALSE,
    verbose_file_output = FALSE, include_overall_pool = FALSE,
    aggregation_methods = "uniform",
    softmax_temperature = 0.5
  )

  # The exact resolver used at the compute_hg2_by_agg() call site, and the
  # OLD (buggy) flat-only read, written out explicitly.
  resolve_site <- function(job) {
    v <- job$compute$softmax_temperature
    if (is.null(v) || !length(v)) v <- job$softmax_temperature
    if (is.null(v) || !length(v)) v <- 1
    as.numeric(v)
  }
  resolve_flat_only <- function(job) {
    v <- job$softmax_temperature
    if (is.null(v) || !length(v)) v <- 1
    as.numeric(v)
  }

  # In-memory: nested under compute, resolver lands on it, flat-only is wrong (1).
  expect_equal(cfg$job$compute$softmax_temperature, 0.5)
  expect_equal(resolve_site(cfg$job), 0.5)
  expect_equal(resolve_flat_only(cfg$job), 1)  # documents the bug the fix avoids

  # Round-trip through the same JSON write/read the SLURM task uses.
  jpath <- tempfile(fileext = ".json")
  job_json <- cfg
  job_json$tables <- NULL
  jsonlite::write_json(job_json, jpath, auto_unbox = TRUE, digits = NA, pretty = TRUE)
  cfg2 <- jsonlite::read_json(jpath, simplifyVector = TRUE)
  expect_equal(resolve_site(cfg2$job), 0.5)
})

test_that("compute_hg2_by_agg: zero-weight NA does not poison max_elbo", {
  set.seed(7)
  n <- 200; p <- 10; K <- 5
  X <- matrix(rnorm(n * p), n, p)
  y <- as.numeric(X[, 1] * 0.5) + rnorm(n)
  fy_list <- lapply(seq_len(K), function(k) as.numeric(X[, 1] * (0.5 + rnorm(1, 0, 0.1))))
  pip_list <- lapply(seq_len(K), function(k) runif(p))
  elbo_vec <- c(10, 1, 1, 1, 1)  # fit 1 dominates -> max_elbo weight on fit 1 only
  # Only the discarded (zero-weight) fits lack second moments.
  unc_list <- list(0.02, NA_real_, NA_real_, NA_real_, NA_real_)
  grr <- tibble::tibble(use_case_id = "uc", group_key = "g")

  tbl <- compute_hg2_by_agg(
    pip_list = pip_list, fitted_y_list = fy_list, elbo_vec = elbo_vec,
    y = y, group_run_row = grr, dataset_bundle_id = "b",
    uncertainty_list = unc_list, pip_cache = NULL
  )
  me <- tbl[tbl$agg_method == "max_elbo", ]
  # max_elbo weights only fit 1 (unc = 0.02); zero-weight NAs must not poison it.
  expect_true(is.finite(me$hg2_uncertainty))
  expect_equal(me$hg2_uncertainty, 0.02, tolerance = 1e-9)
  expect_true(is.finite(me$hg2_expected_pve))
})

test_that("compute_hg2_by_agg: softmax_temperature changes elbo_softmax estimand", {
  set.seed(8)
  n <- 200; p <- 10; K <- 5
  X <- matrix(rnorm(n * p), n, p)
  y <- as.numeric(X[, 1] * 0.4) + rnorm(n)
  # Distinct per-fit fitted-y AND spread-out ELBOs so the softmax weights (and
  # hence the weighted estimand) genuinely depend on temperature.
  fy_list <- lapply(seq_len(K), function(k) as.numeric(X[, 1] * (0.2 + 0.15 * k)))
  pip_list <- lapply(seq_len(K), function(k) runif(p))
  elbo_vec <- c(5, 2, 0, -2, -5)
  unc_list <- as.list(seq(0.01, 0.05, length.out = K))
  grr <- tibble::tibble(use_case_id = "uc", group_key = "g")

  run_at <- function(temp) {
    tbl <- compute_hg2_by_agg(
      pip_list = pip_list, fitted_y_list = fy_list, elbo_vec = elbo_vec,
      y = y, group_run_row = grr, dataset_bundle_id = "b",
      uncertainty_list = unc_list, softmax_temperature = temp, pip_cache = NULL
    )
    tbl[tbl$agg_method == "elbo_softmax", ]
  }
  cold <- run_at(0.25)   # sharper -> closer to max-ELBO fit
  hot  <- run_at(4.0)    # flatter -> closer to uniform

  # Temperature must move the elbo_softmax expected-PVE and its components.
  expect_false(isTRUE(all.equal(cold$hg2_expected_pve, hot$hg2_expected_pve)))
  expect_false(isTRUE(all.equal(cold$hg2_postmean, hot$hg2_postmean)))
  expect_false(isTRUE(all.equal(cold$hg2_uncertainty, hot$hg2_uncertainty)))
  # max_elbo is temperature-invariant (point mass on the best fit) -> unchanged.
  me_cold <- compute_hg2_by_agg(
    pip_list = pip_list, fitted_y_list = fy_list, elbo_vec = elbo_vec,
    y = y, group_run_row = grr, dataset_bundle_id = "b",
    uncertainty_list = unc_list, softmax_temperature = 0.25, pip_cache = NULL
  )
  me_hot <- compute_hg2_by_agg(
    pip_list = pip_list, fitted_y_list = fy_list, elbo_vec = elbo_vec,
    y = y, group_run_row = grr, dataset_bundle_id = "b",
    uncertainty_list = unc_list, softmax_temperature = 4.0, pip_cache = NULL
  )
  expect_equal(
    me_cold$hg2_expected_pve[me_cold$agg_method == "max_elbo"],
    me_hot$hg2_expected_pve[me_hot$agg_method == "max_elbo"],
    tolerance = 1e-12
  )
})

test_that("compute_hg2_by_agg: NA uncertainty propagates, legacy cols intact", {
  set.seed(4)
  n <- 150; p <- 10; K <- 4
  X <- matrix(rnorm(n * p), n, p)
  y <- as.numeric(X[, 1] * 0.5) + rnorm(n)
  fy_list <- lapply(seq_len(K), function(k) as.numeric(X[, 1] * (0.5 + rnorm(1, 0, 0.1))))
  pip_list <- lapply(seq_len(K), function(k) runif(p))
  grr <- tibble::tibble(use_case_id = "uc", group_key = "g")

  # uncertainty_list = NULL -> within-fit term NA, legacy hg2 still finite.
  tbl <- compute_hg2_by_agg(
    pip_list = pip_list, fitted_y_list = fy_list,
    elbo_vec = sort(rnorm(K), decreasing = TRUE),
    y = y, group_run_row = grr, dataset_bundle_id = "b",
    uncertainty_list = NULL, pip_cache = NULL
  )
  expect_true(all(is.na(tbl$hg2_uncertainty)))
  expect_true(all(is.finite(tbl$hg2)))
  expect_true(all(is.finite(tbl$hg2_postmean)))
})
