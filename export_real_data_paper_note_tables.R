#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(tidyr)
  library(purrr)
  library(arrow)
})

args <- commandArgs(trailingOnly = TRUE)
arg_value <- function(flag, default = NULL) {
  i <- which(args == flag)
  if (!length(i) || i[[1]] == length(args)) return(default)
  args[[i[[1]] + 1L]]
}

repo_root <- arg_value("--repo-root", Sys.getenv("TEST_SUSINE_REPO_ROOT", getwd()))
job_name <- arg_value("--job-name", Sys.getenv("REAL_DATA_JOB_NAME", "real_data_ensemble_geometric_n20"))
parent_job_id <- arg_value("--parent-job-id", Sys.getenv("REAL_DATA_PARENT_JOB_ID", "52906940"))
selected_l2_tokens <- strsplit(
  arg_value("--selected-l2-tokens", Sys.getenv("REAL_DATA_SELECTED_L2_TOKENS", "ydjc,arsa,rrp7a,tmtc1,lgals9,znf280b")),
  ",",
  fixed = TRUE
)[[1]]
selected_l2_tokens <- trimws(selected_l2_tokens)
selected_l2_n_variants <- as.integer(arg_value("--n-variants", Sys.getenv("REAL_DATA_SELECTED_L2_N_VARIANTS", "6")))

repo_root <- normalizePath(repo_root, winslash = "/", mustWork = FALSE)
aggregated_dir <- arg_value(
  "--aggregated-dir",
  file.path(repo_root, "output", "slurm_output", job_name, parent_job_id, "aggregated")
)
figure_dir <- arg_value(
  "--figure-dir",
  file.path(repo_root, "output", "slurm_output", job_name, parent_job_id, "figures", "real_data_ensemble")
)
overall_figure_dir <- file.path(figure_dir, "overall")
out_dir <- arg_value("--out-dir", overall_figure_dir)

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

message("repo_root:      ", repo_root)
message("aggregated_dir: ", aggregated_dir)
message("out_dir:        ", out_dir)

required_files <- c(
  file.path(aggregated_dir, "functional_grid_summary.csv"),
  file.path(aggregated_dir, "paper_real_data_ensemble_summary.csv"),
  file.path(aggregated_dir, "aggregated_variant_pips_cluster_weight.parquet")
)
missing_required <- required_files[!file.exists(required_files)]
if (length(missing_required)) {
  stop("Missing required input(s):\n", paste(missing_required, collapse = "\n"), call. = FALSE)
}
if (!dir.exists(file.path(aggregated_dir, "variant_posteriors_dataset"))) {
  stop("Missing required input: ", file.path(aggregated_dir, "variant_posteriors_dataset"), call. = FALSE)
}

read_csv_if_exists <- function(path, empty) {
  if (file.exists(path)) readr::read_csv(path, show_col_types = FALSE) else empty
}

is_abs_path <- function(path) {
  grepl("^(/|[A-Za-z]:[/\\\\])", path)
}

resolve_path <- function(path) {
  if (length(path) != 1L || is.na(path) || !nzchar(path)) return(NA_character_)
  if (is_abs_path(path)) return(path)
  file.path(repo_root, path)
}

js_distance <- function(p, q, eps = 1e-12) {
  p <- as.numeric(p)
  q <- as.numeric(q)
  p <- p + eps
  q <- q + eps
  m <- 0.5 * (p + q)
  kl <- function(a, b) sum(a * log(a / b))
  0.5 * kl(p, m) + 0.5 * kl(q, m)
}

safe_pip_jsd <- function(x, y) {
  if (is.null(x) || is.null(y) || length(x) != length(y) || !length(x)) return(NA_real_)
  as.numeric(js_distance(x, y))
}

safe_pip_l2 <- function(x, y) {
  if (is.null(x) || is.null(y) || length(x) != length(y) || !length(x)) return(NA_real_)
  sqrt(sum((x - y)^2, na.rm = TRUE))
}

safe_pip_credible_shift <- function(x, y) {
  if (is.null(x) || is.null(y) || length(x) != length(y) || !length(x)) return(NA_real_)
  out <- max(pmax(x, y) * abs(x - y), na.rm = TRUE)
  if (is.infinite(out)) NA_real_ else as.numeric(out)
}

core95_entropy <- function(prob, rho = 0.95) {
  p <- as.numeric(prob)
  p[!is.finite(p) | p < 0] <- 0
  s <- sum(p)
  if (!is.finite(s) || s <= 0) return(NA_real_)
  p <- p / s
  ord <- order(p, decreasing = TRUE)
  cum_p <- cumsum(p[ord])
  keep_n <- which(cum_p >= rho)[[1]]
  if (!length(keep_n) || is.na(keep_n)) keep_n <- length(ord)
  core <- p[ord[seq_len(keep_n)]]
  core <- core / sum(core)
  -sum(core[core > 0] * log(core[core > 0]))
}

functional_grid_summary <- readr::read_csv(
  file.path(aggregated_dir, "functional_grid_summary.csv"),
  show_col_types = FALSE
)
paper_summary <- readr::read_csv(
  file.path(aggregated_dir, "paper_real_data_ensemble_summary.csv"),
  show_col_types = FALSE
)
highest_weight_refit_basin <- read_csv_if_exists(
  file.path(aggregated_dir, "highest_weight_refit_basin_r2_drift.csv"),
  tibble::tibble(locus_id = character())
)
aggregated_variants <- arrow::read_parquet(
  file.path(aggregated_dir, "aggregated_variant_pips_cluster_weight.parquet")
) %>%
  tibble::as_tibble()
highest_weight_refit_variants <- if (file.exists(file.path(aggregated_dir, "highest_weight_refit_variant_posteriors.parquet"))) {
  arrow::read_parquet(file.path(aggregated_dir, "highest_weight_refit_variant_posteriors.parquet")) %>%
    tibble::as_tibble()
} else {
  tibble::tibble(locus_id = character(), variant_id = character(), pip = numeric())
}
variant_posteriors <- arrow::open_dataset(file.path(aggregated_dir, "variant_posteriors_dataset"))
loci <- sort(unique(functional_grid_summary$locus_id))

nearest_reference_run <- function(locus_id, target_c = 0, target_sigma = 0.2) {
  df <- functional_grid_summary %>%
    dplyr::filter(
      .data$locus_id == !!locus_id,
      is.finite(.data$c_value),
      is.finite(.data$sigma_0_2_scalar)
    ) %>%
    dplyr::mutate(
      c_distance = abs(.data$c_value - target_c),
      sigma_distance = abs(log(.data$sigma_0_2_scalar / target_sigma))
    ) %>%
    dplyr::arrange(.data$c_distance, .data$sigma_distance, .data$run_id)
  if (!nrow(df)) return(NA_integer_)
  as.integer(df$run_id[[1]])
}

reference_pip_tbl <- function(locus_id) {
  ref_run_id <- nearest_reference_run(locus_id)
  if (is.na(ref_run_id)) {
    return(tibble::tibble(
      variant_id = character(),
      baseline_pip = numeric(),
      baseline_posterior_mean = numeric(),
      baseline_conditional_effect = numeric()
    ))
  }
  df <- variant_posteriors %>%
    dplyr::filter(.data$locus_id == !!locus_id, .data$run_id == !!ref_run_id) %>%
    dplyr::collect()
  if (!"posterior_mean" %in% names(df)) df$posterior_mean <- NA_real_
  df %>%
    dplyr::transmute(
      .data$variant_id,
      baseline_pip = .data$pip,
      baseline_posterior_mean = .data$posterior_mean,
      baseline_conditional_effect = dplyr::if_else(
        is.finite(.data$pip) & .data$pip >= 0.001,
        .data$posterior_mean / .data$pip,
        0
      )
    )
}

refit_pip_tbl <- function(locus_id) {
  if (!nrow(highest_weight_refit_variants)) {
    return(tibble::tibble(
      variant_id = character(),
      refit_pip = numeric(),
      refit_posterior_mean = numeric(),
      refit_conditional_effect = numeric()
    ))
  }
  df <- highest_weight_refit_variants %>%
    dplyr::filter(.data$locus_id == !!locus_id)
  if (!"posterior_mean" %in% names(df)) df$posterior_mean <- NA_real_
  df %>%
    dplyr::transmute(
      .data$variant_id,
      refit_pip = .data$pip,
      refit_posterior_mean = .data$posterior_mean,
      refit_conditional_effect = dplyr::if_else(
        is.finite(.data$pip) & .data$pip >= 0.001,
        .data$posterior_mean / .data$pip,
        0
      )
    )
}

run_pip_vec <- function(locus_id, run_id) {
  if (length(run_id) != 1L || is.na(run_id)) return(NULL)
  x <- variant_posteriors %>%
    dplyr::filter(.data$locus_id == !!locus_id, .data$run_id == !!as.integer(run_id)) %>%
    dplyr::select(.data$ld_matrix_index, .data$pip) %>%
    dplyr::collect() %>%
    dplyr::arrange(.data$ld_matrix_index) %>%
    dplyr::pull(.data$pip)
  if (!length(x)) return(NULL)
  x
}

run_pip_jsd <- function(locus_id, run_id_a, run_id_b) {
  safe_pip_jsd(run_pip_vec(locus_id, run_id_a), run_pip_vec(locus_id, run_id_b))
}

refit_pip_jsd <- function(locus_id, anchor_run_id) {
  if (length(anchor_run_id) != 1L || is.na(anchor_run_id)) return(NA_real_)
  anchor_pip <- variant_posteriors %>%
    dplyr::filter(.data$locus_id == !!locus_id, .data$run_id == !!as.integer(anchor_run_id)) %>%
    dplyr::select(.data$variant_id, anchor_pip = .data$pip) %>%
    dplyr::collect()
  refit_pip <- refit_pip_tbl(locus_id) %>%
    dplyr::select(.data$variant_id, .data$refit_pip)
  joined <- dplyr::inner_join(anchor_pip, refit_pip, by = "variant_id")
  if (!nrow(joined)) return(NA_real_)
  safe_pip_jsd(joined$anchor_pip, joined$refit_pip)
}

weighted_annotation_tbl <- functional_grid_summary %>%
  dplyr::filter(
    is.finite(.data$c_value),
    is.finite(.data$sigma_0_2_scalar),
    is.finite(.data$agg_weight_run)
  ) %>%
  dplyr::group_by(.data$locus_id, .data$gene_name) %>%
  dplyr::summarise(
    weight_total = sum(.data$agg_weight_run, na.rm = TRUE),
    weighted_c_value = sum(.data$agg_weight_run * .data$c_value, na.rm = TRUE) /
      sum(.data$agg_weight_run, na.rm = TRUE),
    weighted_sigma_0_2_scalar = sum(.data$agg_weight_run * .data$sigma_0_2_scalar, na.rm = TRUE) /
      sum(.data$agg_weight_run, na.rm = TRUE),
    off_c0_weight = sum(.data$agg_weight_run[abs(.data$c_value) > 1e-12], na.rm = TRUE) /
      sum(.data$agg_weight_run, na.rm = TRUE),
    max_model_weight = max(.data$agg_weight_run, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  dplyr::filter(
    is.finite(.data$weight_total),
    .data$weight_total > 0,
    is.finite(.data$weighted_c_value),
    is.finite(.data$weighted_sigma_0_2_scalar)
  )

effect_diffuseness_run_tbl <- paper_summary %>%
  dplyr::select(.data$locus_id, .data$susie_anchor_run_id, .data$highest_weight_source_run_id) %>%
  tidyr::pivot_longer(
    cols = c("susie_anchor_run_id", "highest_weight_source_run_id"),
    names_to = "fit_type",
    values_to = "run_id"
  ) %>%
  dplyr::mutate(
    fit_type = dplyr::recode(
      .data$fit_type,
      susie_anchor_run_id = "SuSiE baseline",
      highest_weight_source_run_id = "Highest-weight model"
    ),
    fit_type = factor(.data$fit_type, levels = c("SuSiE baseline", "Highest-weight model"))
  ) %>%
  dplyr::filter(!is.na(.data$run_id)) %>%
  dplyr::distinct(.data$locus_id, .data$fit_type, .data$run_id)

effect_posteriors_path <- file.path(aggregated_dir, "effect_posteriors_dataset")
effect_diffuseness_tbl <- if (nrow(effect_diffuseness_run_tbl) && dir.exists(effect_posteriors_path)) {
  arrow::open_dataset(effect_posteriors_path) %>%
    dplyr::filter(.data$run_id %in% !!effect_diffuseness_run_tbl$run_id) %>%
    dplyr::select(.data$run_id, .data$locus_id, .data$effect_l, .data$alpha) %>%
    dplyr::collect() %>%
    dplyr::inner_join(effect_diffuseness_run_tbl, by = c("run_id", "locus_id")) %>%
    dplyr::group_by(.data$locus_id, .data$run_id, .data$fit_type, .data$effect_l) %>%
    dplyr::summarise(
      effect_pip_entropy_core95 = core95_entropy(.data$alpha),
      effect_k_eff_core95 = exp(.data$effect_pip_entropy_core95),
      alpha_mass = sum(.data$alpha, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    dplyr::filter(is.finite(.data$effect_k_eff_core95), .data$alpha_mass > 1e-8)
} else {
  warning("effect_posteriors_dataset is missing or no run ids are available; writing empty diffuseness counts.")
  tibble::tibble(
    locus_id = character(),
    run_id = integer(),
    fit_type = factor(levels = c("SuSiE baseline", "Highest-weight model")),
    effect_l = integer(),
    effect_pip_entropy_core95 = numeric(),
    effect_k_eff_core95 = numeric(),
    alpha_mass = numeric()
  )
}

diffuseness_breaks <- c(0, 1.5, 5, 20, 100, 2500, 4000, Inf)
diffuseness_labels <- c("0-1.5", "1.5-5", "5-20", "20-100", "100-2.5k", "2.5k-4k", "4k+")
effect_diffuseness_count_tbl <- effect_diffuseness_tbl %>%
  dplyr::mutate(
    diffuseness_bin = cut(
      .data$effect_k_eff_core95,
      breaks = diffuseness_breaks,
      labels = diffuseness_labels,
      include.lowest = TRUE,
      right = FALSE
    )
  ) %>%
  dplyr::filter(!is.na(.data$diffuseness_bin)) %>%
  dplyr::count(.data$fit_type, .data$diffuseness_bin, name = "n") %>%
  tidyr::complete(
    fit_type = factor(c("SuSiE baseline", "Highest-weight model"),
                      levels = c("SuSiE baseline", "Highest-weight model")),
    diffuseness_bin = factor(diffuseness_labels, levels = diffuseness_labels),
    fill = list(n = 0L)
  ) %>%
  tidyr::pivot_wider(names_from = .data$fit_type, values_from = .data$n) %>%
  dplyr::mutate(
    delta_highest_minus_baseline = .data$`Highest-weight model` - .data$`SuSiE baseline`
  )

run_metrics_full <- read_csv_if_exists(
  file.path(aggregated_dir, "run_metrics_full.csv"),
  tibble::tibble(locus_id = character(), run_id = integer(), fit_rds_path = character())
) %>%
  dplyr::select(any_of(c("locus_id", "run_id", "fit_rds_path")))
aggregation_weights <- read_csv_if_exists(
  file.path(aggregated_dir, "aggregation_weights_cluster_weight.csv"),
  tibble::tibble(locus_id = character(), run_id = integer(), agg_weight_run = numeric())
) %>%
  dplyr::select(any_of(c("locus_id", "run_id", "agg_weight_run")))

floored_conditional_alpha_floor <- 0.0001
floored_conditional_diagnostics <- list()
floored_conditional_cache <- new.env(parent = emptyenv())

log_floored_conditional <- function(locus_id, kind, detail = NA_character_) {
  floored_conditional_diagnostics[[length(floored_conditional_diagnostics) + 1L]] <<-
    tibble::tibble(locus_id = as.character(locus_id), kind = as.character(kind), detail = as.character(detail))
}

floored_conditional_for_fit_object <- function(fit, alpha_floor = floored_conditional_alpha_floor) {
  if (is.null(fit)) return(NULL)
  alpha <- NULL
  mu <- NULL
  if (!is.null(fit$effect_fits$alpha) && !is.null(fit$effect_fits$b_hat)) {
    alpha <- as.matrix(fit$effect_fits$alpha)
    alpha_b_hat <- as.matrix(fit$effect_fits$b_hat)
    if (!all(dim(alpha) == dim(alpha_b_hat))) return(NULL)
    mu <- ifelse(alpha > 0, alpha_b_hat / alpha, 0)
  } else if (!is.null(fit$alpha) && !is.null(fit$mu)) {
    alpha <- as.matrix(fit$alpha)
    mu <- as.matrix(fit$mu)
    if (!all(dim(alpha) == dim(mu))) return(NULL)
  } else {
    return(NULL)
  }
  alpha_clip <- pmax(alpha, alpha_floor)
  list(
    num = colSums(alpha_clip * mu, na.rm = TRUE),
    denom = colSums(alpha_clip, na.rm = TRUE),
    p = ncol(alpha)
  )
}

floored_conditional_for_locus <- function(locus_id) {
  key <- as.character(locus_id)
  if (exists(key, envir = floored_conditional_cache, inherits = FALSE)) {
    return(get(key, envir = floored_conditional_cache, inherits = FALSE))
  }

  baseline_run_id <- nearest_reference_run(locus_id)
  variant_order <- if (!is.na(baseline_run_id)) {
    variant_posteriors %>%
      dplyr::filter(.data$locus_id == !!locus_id, .data$run_id == !!baseline_run_id) %>%
      dplyr::select(.data$variant_id, .data$ld_matrix_index) %>%
      dplyr::collect() %>%
      dplyr::arrange(.data$ld_matrix_index)
  } else {
    tibble::tibble(variant_id = character(), ld_matrix_index = integer())
  }
  p <- nrow(variant_order)
  na_vec <- rep(NA_real_, p)
  empty_result <- tibble::tibble(
    variant_id = variant_order$variant_id,
    baseline_conditional_effect_floored = na_vec,
    aggregated_conditional_effect_floored = na_vec,
    refit_conditional_effect_floored = na_vec
  )
  if (!p) {
    assign(key, empty_result, envir = floored_conditional_cache)
    return(empty_result)
  }

  path_for <- function(rid) {
    if (is.na(rid) || !"fit_rds_path" %in% names(run_metrics_full)) return(NA_character_)
    out <- run_metrics_full %>%
      dplyr::filter(.data$locus_id == !!locus_id, .data$run_id == !!rid) %>%
      dplyr::pull(.data$fit_rds_path)
    if (!length(out)) NA_character_ else resolve_path(out[[1]])
  }

  compute_for_path <- function(path, role, run_id = NA) {
    detail <- if (!is.na(run_id)) paste0("run_id=", run_id, ", path=", path) else paste0("path=", path)
    if (length(path) != 1L || is.na(path) || !nzchar(path) || !file.exists(path)) {
      log_floored_conditional(locus_id, paste0(role, "_fit_unreadable"), detail)
      return(NULL)
    }
    fit <- tryCatch(readRDS(path), error = function(e) {
      log_floored_conditional(locus_id, paste0(role, "_fit_unreadable"), paste(conditionMessage(e), detail))
      NULL
    })
    fc <- floored_conditional_for_fit_object(fit)
    if (is.null(fc) || fc$p != p) {
      log_floored_conditional(locus_id, paste0(role, "_fit_invalid"), detail)
      return(NULL)
    }
    fc
  }

  baseline_fc <- if (!is.na(baseline_run_id)) {
    compute_for_path(path_for(baseline_run_id), role = "baseline", run_id = baseline_run_id)
  } else {
    log_floored_conditional(locus_id, "baseline_run_missing", "nearest_reference_run returned NA")
    NULL
  }

  refit_dir <- file.path(aggregated_dir, "highest_weight_refits", locus_id)
  refit_files <- if (dir.exists(refit_dir)) list.files(refit_dir, pattern = "\\.rds$", full.names = TRUE) else character(0)
  refit_fc <- if (length(refit_files)) {
    compute_for_path(refit_files[[1]], role = "refit")
  } else {
    log_floored_conditional(locus_id, "refit_file_missing", paste0("no .rds files under ", refit_dir))
    NULL
  }

  nominees <- aggregation_weights %>%
    dplyr::filter(.data$locus_id == !!locus_id, .data$agg_weight_run > 0)
  ens_num <- rep(0, p)
  ens_denom <- rep(0, p)
  ens_valid <- FALSE
  if (nrow(nominees)) {
    for (i in seq_len(nrow(nominees))) {
      rid <- nominees$run_id[[i]]
      w <- nominees$agg_weight_run[[i]]
      fc <- compute_for_path(path_for(rid), role = "ensemble_nominee", run_id = rid)
      if (is.null(fc)) next
      ens_num <- ens_num + w * fc$num
      ens_denom <- ens_denom + w * fc$denom
      ens_valid <- TRUE
    }
  } else {
    log_floored_conditional(locus_id, "ensemble_no_nominees", "agg_weight_run > 0 returned 0 rows")
  }

  safe_div <- function(n, d) as.numeric(ifelse(d > 0, n / d, NA_real_))
  result <- tibble::tibble(
    variant_id = variant_order$variant_id,
    baseline_conditional_effect_floored = if (!is.null(baseline_fc)) safe_div(baseline_fc$num, baseline_fc$denom) else na_vec,
    aggregated_conditional_effect_floored = if (ens_valid) safe_div(ens_num, ens_denom) else na_vec,
    refit_conditional_effect_floored = if (!is.null(refit_fc)) safe_div(refit_fc$num, refit_fc$denom) else na_vec
  )
  assign(key, result, envir = floored_conditional_cache)
  result
}

selected_lollipop_variant_tbl <- function(
    locus_id,
    n_variants = 6L,
    selection_metric = c("pip_delta", "posterior_mean_delta")
) {
  n_variants <- as.integer(n_variants)
  selection_metric <- match.arg(selection_metric)
  locus_aggregated_variants <- aggregated_variants %>%
    dplyr::filter(.data$locus_id == !!locus_id)
  if (!"aggregated_posterior_mean" %in% names(locus_aggregated_variants)) {
    locus_aggregated_variants$aggregated_posterior_mean <- NA_real_
  }
  floored_cond_tbl <- floored_conditional_for_locus(locus_id)

  plot_df <- locus_aggregated_variants %>%
    dplyr::left_join(reference_pip_tbl(locus_id), by = "variant_id") %>%
    dplyr::left_join(refit_pip_tbl(locus_id), by = "variant_id") %>%
    dplyr::left_join(floored_cond_tbl, by = "variant_id") %>%
    dplyr::mutate(
      baseline_conditional_effect = .data$baseline_conditional_effect_floored,
      refit_conditional_effect = .data$refit_conditional_effect_floored,
      aggregated_conditional_effect = .data$aggregated_conditional_effect_floored,
      abs_delta_baseline_ensemble = abs(.data$aggregated_pip - .data$baseline_pip),
      abs_delta_posterior_mean = abs(.data$aggregated_posterior_mean - .data$baseline_posterior_mean),
      selection_score = dplyr::case_when(
        .env$selection_metric == "posterior_mean_delta" ~ .data$abs_delta_posterior_mean,
        TRUE ~ .data$abs_delta_baseline_ensemble
      )
    ) %>%
    dplyr::filter(is.finite(.data$selection_score)) %>%
    dplyr::arrange(dplyr::desc(.data$selection_score), .data$ld_matrix_index) %>%
    dplyr::slice_head(n = n_variants)

  if (!nrow(plot_df)) {
    plot_df <- locus_aggregated_variants %>%
      dplyr::arrange(dplyr::desc(.data$aggregated_pip), .data$ld_matrix_index) %>%
      dplyr::slice_head(n = n_variants) %>%
      dplyr::left_join(reference_pip_tbl(locus_id), by = "variant_id") %>%
      dplyr::left_join(refit_pip_tbl(locus_id), by = "variant_id") %>%
      dplyr::left_join(floored_cond_tbl, by = "variant_id") %>%
      dplyr::mutate(
        baseline_conditional_effect = .data$baseline_conditional_effect_floored,
        refit_conditional_effect = .data$refit_conditional_effect_floored,
        aggregated_conditional_effect = .data$aggregated_conditional_effect_floored,
        abs_delta_posterior_mean = abs(.data$aggregated_posterior_mean - .data$baseline_posterior_mean),
        selection_score = .data$abs_delta_posterior_mean
      )
  }

  plot_df %>%
    dplyr::arrange(dplyr::desc(.data$aggregated_pip), .data$ld_matrix_index)
}

resolve_locus_token <- function(token, locus_ids) {
  hit <- locus_ids[grepl(paste0("^", token), locus_ids, ignore.case = TRUE)]
  if (!length(hit)) hit <- locus_ids[grepl(token, locus_ids, ignore.case = TRUE)]
  if (!length(hit)) stop("Could not resolve selected locus token: ", token, call. = FALSE)
  hit[[1]]
}

selected_zoom_l2_loci <- vapply(
  selected_l2_tokens,
  resolve_locus_token,
  character(1),
  locus_ids = loci
)

build_selected_drift_tbl <- function(selected_loci) {
  purrr::map_dfr(selected_loci, function(loc) {
    ref_pip <- reference_pip_tbl(loc)
    run_tbl <- functional_grid_summary %>%
      dplyr::filter(.data$locus_id == !!loc) %>%
      dplyr::arrange(.data$sigma_0_2_scalar, .data$c_value)
    pip_tbl <- variant_posteriors %>%
      dplyr::filter(.data$locus_id == !!loc, .data$run_id %in% run_tbl$run_id) %>%
      dplyr::select(.data$run_id, .data$variant_id, .data$pip) %>%
      dplyr::collect() %>%
      dplyr::inner_join(ref_pip, by = "variant_id") %>%
      dplyr::arrange(.data$variant_id)
    pip_tbl_by_run <- split(pip_tbl, pip_tbl$run_id)
    drift_for_run <- function(rid, kind) {
      x <- pip_tbl_by_run[[as.character(rid)]]
      if (is.null(x) || !nrow(x)) return(NA_real_)
      if (identical(kind, "jsd")) {
        js_distance(x$pip, x$baseline_pip)
      } else if (identical(kind, "l2")) {
        safe_pip_l2(x$pip, x$baseline_pip)
      } else if (identical(kind, "credible")) {
        safe_pip_credible_shift(x$pip, x$baseline_pip)
      } else {
        stop("Unknown drift metric: ", kind, call. = FALSE)
      }
    }
    run_tbl %>%
      dplyr::mutate(
        jsd_from_reference_c0_sigma02 = purrr::map_dbl(.data$run_id, drift_for_run, kind = "jsd"),
        l2_from_reference_c0_sigma02 = purrr::map_dbl(.data$run_id, drift_for_run, kind = "l2"),
        credible_from_reference_c0_sigma02 = purrr::map_dbl(.data$run_id, drift_for_run, kind = "credible"),
        locus_id = loc,
        .before = 1L
      )
  })
}

selected_l2_drift_tbl <- build_selected_drift_tbl(selected_zoom_l2_loci)
selected_l2_lollipop_source_tbl <- purrr::map_dfr(
  selected_zoom_l2_loci,
  ~ selected_lollipop_variant_tbl(.x, n_variants = selected_l2_n_variants, selection_metric = "posterior_mean_delta")
) %>%
  dplyr::select(
    .data$locus_id,
    .data$variant_id,
    .data$ld_matrix_index,
    .data$selection_score,
    .data$baseline_pip,
    .data$aggregated_pip,
    .data$refit_pip,
    .data$baseline_posterior_mean,
    .data$aggregated_posterior_mean,
    .data$refit_posterior_mean,
    .data$baseline_conditional_effect,
    .data$aggregated_conditional_effect,
    .data$refit_conditional_effect,
    .data$abs_delta_posterior_mean
  )

basis_cols <- c(
  "basis_drift_var_y_highest_weight_vs_susie_anchor",
  "basis_drift_var_y_refit_vs_susie_anchor"
)
for (nm in basis_cols) {
  if (!nm %in% names(highest_weight_refit_basin)) highest_weight_refit_basin[[nm]] <- NA_real_
}

drift_source_tbl <- paper_summary %>%
  dplyr::left_join(
    highest_weight_refit_basin %>%
      dplyr::select(.data$locus_id, dplyr::all_of(basis_cols)),
    by = "locus_id"
  ) %>%
  dplyr::mutate(
    baseline_reference_run_id = purrr::map_int(.data$locus_id, nearest_reference_run),
    highest_weight_pip_jsd = purrr::map2_dbl(
      .data$locus_id,
      .data$highest_weight_source_run_id,
      ~ run_pip_jsd(.x, nearest_reference_run(.x), .y)
    ),
    warm_refit_pip_jsd = purrr::map2_dbl(
      .data$locus_id,
      .data$baseline_reference_run_id,
      refit_pip_jsd
    )
  )

drift_jsd_long_tbl <- dplyr::bind_rows(
  drift_source_tbl %>%
    dplyr::transmute(
      .data$locus_id,
      fit_type = "Highest-weight model",
      basis_drift_var_y = .data$basis_drift_var_y_highest_weight_vs_susie_anchor,
      pip_jsd = .data$highest_weight_pip_jsd,
      source_weight = .data$highest_weight_source_agg_weight_run
    ),
  drift_source_tbl %>%
    dplyr::transmute(
      .data$locus_id,
      fit_type = "Warm refit",
      basis_drift_var_y = .data$basis_drift_var_y_refit_vs_susie_anchor,
      pip_jsd = .data$warm_refit_pip_jsd,
      source_weight = .data$highest_weight_source_agg_weight_run
    )
) %>%
  dplyr::arrange(.data$locus_id, .data$fit_type)

readr::write_csv(
  weighted_annotation_tbl,
  file.path(out_dir, "weighted_annotation_tbl_for_paper_notes.csv")
)
readr::write_csv(
  drift_jsd_long_tbl,
  file.path(out_dir, "drift_jsd_long_tbl_for_paper_notes.csv")
)
readr::write_csv(
  effect_diffuseness_count_tbl,
  file.path(out_dir, "effect_diffuseness_count_tbl_for_paper_notes.csv")
)
readr::write_csv(
  selected_l2_drift_tbl,
  file.path(out_dir, "selected_l2_drift_for_paper_notes.csv")
)
readr::write_csv(
  selected_l2_lollipop_source_tbl,
  file.path(out_dir, "selected_l2_lollipop_source_for_paper_notes.csv")
)

diagnostics_tbl <- if (length(floored_conditional_diagnostics)) {
  dplyr::bind_rows(floored_conditional_diagnostics)
} else {
  tibble::tibble(locus_id = character(), kind = character(), detail = character())
}
readr::write_csv(
  diagnostics_tbl,
  file.path(out_dir, "selected_l2_lollipop_floored_conditional_diagnostics.csv")
)

message("Wrote paper-note helper CSVs:")
message("  ", file.path(out_dir, "weighted_annotation_tbl_for_paper_notes.csv"))
message("  ", file.path(out_dir, "drift_jsd_long_tbl_for_paper_notes.csv"))
message("  ", file.path(out_dir, "effect_diffuseness_count_tbl_for_paper_notes.csv"))
message("  ", file.path(out_dir, "selected_l2_drift_for_paper_notes.csv"))
message("  ", file.path(out_dir, "selected_l2_lollipop_source_for_paper_notes.csv"))
message("  ", file.path(out_dir, "selected_l2_lollipop_floored_conditional_diagnostics.csv"))
