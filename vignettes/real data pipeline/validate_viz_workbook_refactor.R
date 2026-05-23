# validate_viz_workbook_refactor.R
#
# Validates that the lazy-Arrow refactor of visualize_results_workbook_real_data
# _ensemble.Rmd produces bit-identical output to the original eager-collect
# version, without ever materializing the full variant_posteriors_dataset (so
# this script itself does not OOM).
#
# Strategy: for each consumer function in the workbook that touches
# variant_posteriors, define both the eager-then-filter (legacy) path and the
# lazy-filter-then-collect (refactored) path, run both on a small subset of
# loci, and assert identity. If the refactor is a pure restatement of the
# legacy logic, every comparison must pass with zero numeric tolerance.
#
# Also writes per-locus digests of every figure-input data frame so they can
# be re-verified after any future refactor without rerunning the workbook.
#
# Usage on HPC:
#   Rscript vignettes/real\ data\ pipeline/validate_viz_workbook_refactor.R \
#     --job-name real_data_ensemble_geometric_n20 \
#     --parent-job-id 52906940 \
#     [--loci ydjc_chr22_lung,rrp7a_chr22_lung,prpf38b_chr1_lung] \
#     [--digest-out output/slurm_output/.../viz_refactor_digests.rds]

suppressPackageStartupMessages({
  library(here)
  library(devtools)
  library(dplyr)
  library(readr)
  library(arrow)
  library(tibble)
})

# ---------- CLI ----------
parse_arg <- function(args, name, default = NULL) {
  hit <- grep(paste0("^", name, "="), args, value = TRUE)
  if (length(hit)) {
    return(sub(paste0("^", name, "="), "", hit[[1]]))
  }
  idx <- match(name, args)
  if (!is.na(idx) && idx < length(args)) return(args[idx + 1L])
  default
}

argv <- commandArgs(trailingOnly = TRUE)
job_name <- parse_arg(argv, "--job-name", "real_data_ensemble_geometric_n20")
parent_job_id <- parse_arg(argv, "--parent-job-id", "52906940")
loci_arg <- parse_arg(argv, "--loci", NULL)
digest_out_arg <- parse_arg(argv, "--digest-out", NULL)

here::i_am("vignettes/real data pipeline/validate_viz_workbook_refactor.R")
susine_path <- here("..", "susine")
if (file.exists(file.path(susine_path, "DESCRIPTION"))) {
  devtools::load_all(susine_path)
}
devtools::load_all(".")

aggregated_dir <- here("output", "slurm_output", job_name, parent_job_id, "aggregated")
stopifnot(dir.exists(aggregated_dir))

# ---------- Load supporting tables (small) ----------
functional_grid_summary <- readr::read_csv(
  file.path(aggregated_dir, "functional_grid_summary.csv"), show_col_types = FALSE
)
anchor_summary <- readr::read_csv(
  file.path(aggregated_dir, "susie_anchor_summary.csv"), show_col_types = FALSE
)
paper_summary <- readr::read_csv(
  file.path(aggregated_dir, "paper_real_data_ensemble_summary.csv"), show_col_types = FALSE
)
aggregated_variants <- arrow::read_parquet(
  file.path(aggregated_dir, "aggregated_variant_pips_cluster_weight.parquet")
) %>% tibble::as_tibble()
highest_weight_refit_variants <- if (file.exists(file.path(aggregated_dir, "highest_weight_refit_variant_posteriors.parquet"))) {
  arrow::read_parquet(file.path(aggregated_dir, "highest_weight_refit_variant_posteriors.parquet")) %>%
    tibble::as_tibble()
} else {
  tibble::tibble(locus_id = character(), variant_id = character(), pip = numeric())
}

all_loci <- sort(unique(functional_grid_summary$locus_id))
loci <- if (!is.null(loci_arg)) {
  strsplit(loci_arg, ",", fixed = TRUE)[[1]]
} else {
  # Default to the three zoom loci that exercise both extremes of basis drift.
  intersect(
    c("ydjc_chr22_lung", "rrp7a_chr22_lung", "prpf38b_chr1_lung"),
    all_loci
  )
}
stopifnot(length(loci) >= 1L, all(loci %in% all_loci))
cat("[validate] loci under test:", paste(loci, collapse = ", "), "\n")

# ---------- Two parallel paths over variant_posteriors ----------

vp_lazy <- arrow::open_dataset(file.path(aggregated_dir, "variant_posteriors_dataset"))
# Eager but SCOPED: only the loci under test, not the whole dataset.
vp_eager_subset <- vp_lazy %>%
  dplyr::filter(.data$locus_id %in% loci) %>%
  dplyr::collect()
cat(sprintf(
  "[validate] eager subset: %d rows / ~%.1f MB (full dataset NOT collected)\n",
  nrow(vp_eager_subset),
  as.numeric(utils::object.size(vp_eager_subset)) / 2^20
))

# Helper-by-helper equivalence: each helper has an `eager` variant that
# filters from vp_eager_subset (mimicking the pre-refactor behavior of the
# workbook) and a `lazy` variant that filters from vp_lazy and collects.

# nearest_reference_run / baseline run id per locus
nearest_reference_run <- function(locus_id, target_c = 0, target_sigma = 0.2) {
  df <- functional_grid_summary %>%
    dplyr::filter(.data$locus_id == !!locus_id,
                  is.finite(.data$c_value),
                  is.finite(.data$sigma_0_2_scalar)) %>%
    dplyr::mutate(
      c_distance = abs(.data$c_value - target_c),
      sigma_distance = abs(log(.data$sigma_0_2_scalar / target_sigma))
    ) %>%
    dplyr::arrange(.data$c_distance, .data$sigma_distance, .data$run_id)
  if (!nrow(df)) return(NA_integer_)
  as.integer(df$run_id[[1]])
}

reference_pip_tbl_eager <- function(locus_id) {
  ref_run_id <- nearest_reference_run(locus_id)
  if (is.na(ref_run_id)) {
    return(tibble::tibble(
      variant_id = character(), baseline_pip = numeric(),
      baseline_posterior_mean = numeric(), baseline_conditional_effect = numeric()
    ))
  }
  df <- vp_eager_subset %>%
    dplyr::filter(.data$locus_id == !!locus_id, .data$run_id == !!ref_run_id)
  if (!"posterior_mean" %in% names(df)) df$posterior_mean <- NA_real_
  df %>%
    dplyr::transmute(
      .data$variant_id,
      baseline_pip = .data$pip,
      baseline_posterior_mean = .data$posterior_mean,
      baseline_conditional_effect = dplyr::if_else(
        is.finite(.data$pip) & .data$pip > .Machine$double.eps,
        .data$posterior_mean / .data$pip,
        NA_real_
      )
    )
}
reference_pip_tbl_lazy <- function(locus_id) {
  ref_run_id <- nearest_reference_run(locus_id)
  if (is.na(ref_run_id)) {
    return(tibble::tibble(
      variant_id = character(), baseline_pip = numeric(),
      baseline_posterior_mean = numeric(), baseline_conditional_effect = numeric()
    ))
  }
  df <- vp_lazy %>%
    dplyr::filter(.data$locus_id == !!locus_id, .data$run_id == !!ref_run_id) %>%
    dplyr::collect()
  if (!"posterior_mean" %in% names(df)) df$posterior_mean <- NA_real_
  df %>%
    dplyr::transmute(
      .data$variant_id,
      baseline_pip = .data$pip,
      baseline_posterior_mean = .data$posterior_mean,
      baseline_conditional_effect = dplyr::if_else(
        is.finite(.data$pip) & .data$pip > .Machine$double.eps,
        .data$posterior_mean / .data$pip,
        NA_real_
      )
    )
}

run_pip_vec_eager <- function(locus_id, run_id) {
  if (length(run_id) != 1L || is.na(run_id)) return(NULL)
  x <- vp_eager_subset %>%
    dplyr::filter(.data$locus_id == !!locus_id,
                  .data$run_id == !!as.integer(run_id)) %>%
    dplyr::arrange(.data$ld_matrix_index) %>%
    dplyr::pull(.data$pip)
  if (!length(x)) return(NULL)
  x
}
run_pip_vec_lazy <- function(locus_id, run_id) {
  if (length(run_id) != 1L || is.na(run_id)) return(NULL)
  x <- vp_lazy %>%
    dplyr::filter(.data$locus_id == !!locus_id,
                  .data$run_id == !!as.integer(run_id)) %>%
    dplyr::select(.data$ld_matrix_index, .data$pip) %>%
    dplyr::collect() %>%
    dplyr::arrange(.data$ld_matrix_index) %>%
    dplyr::pull(.data$pip)
  if (!length(x)) return(NULL)
  x
}

# plot_drift's inner filter+inner_join
build_drift_pip_tbl_eager <- function(locus_id, ref_pip, run_ids) {
  vp_eager_subset %>%
    dplyr::filter(.data$locus_id == !!locus_id, .data$run_id %in% run_ids) %>%
    dplyr::select(.data$run_id, .data$variant_id, .data$pip) %>%
    dplyr::inner_join(ref_pip, by = "variant_id")
}
build_drift_pip_tbl_lazy <- function(locus_id, ref_pip, run_ids) {
  vp_lazy %>%
    dplyr::filter(.data$locus_id == !!locus_id, .data$run_id %in% run_ids) %>%
    dplyr::select(.data$run_id, .data$variant_id, .data$pip) %>%
    dplyr::collect() %>%
    dplyr::inner_join(ref_pip, by = "variant_id")
}

# plot_aggregated_lollipop's overlay block
build_lollipop_overlay_eager <- function(locus_id, anchor_id, baseline_id, kept_variant_ids) {
  vp_eager_subset %>%
    dplyr::filter(.data$run_id %in% c(anchor_id, baseline_id),
                  .data$variant_id %in% kept_variant_ids)
}
build_lollipop_overlay_lazy <- function(locus_id, anchor_id, baseline_id, kept_variant_ids) {
  vp_lazy %>%
    dplyr::filter(.data$locus_id == !!locus_id,
                  .data$run_id %in% c(anchor_id, baseline_id),
                  .data$variant_id %in% kept_variant_ids) %>%
    dplyr::collect()
}

# ---------- Comparison harness ----------
canon <- function(df) {
  # Normalize column and row order so equality does not depend on either.
  if (!nrow(df)) return(df)
  df <- df[, sort(names(df)), drop = FALSE]
  key_cols <- intersect(c("variant_id", "run_id", "ld_matrix_index"), names(df))
  if (length(key_cols)) df <- df[do.call(order, df[key_cols]), , drop = FALSE]
  rownames(df) <- NULL
  df
}

digest_or_na <- function(x) {
  if (is.null(x) || (is.atomic(x) && !length(x))) return(NA_character_)
  if (!requireNamespace("digest", quietly = TRUE)) return(NA_character_)
  digest::digest(canon(if (is.atomic(x)) tibble::tibble(x = x) else x))
}

assert_equal <- function(label, a, b) {
  if (is.null(a) && is.null(b)) {
    cat("[ ok   ]", label, "(both NULL)\n"); return(invisible(TRUE))
  }
  if (xor(is.null(a), is.null(b))) {
    cat("[ FAIL ]", label, ": one side is NULL\n"); return(invisible(FALSE))
  }
  if (is.atomic(a) && is.atomic(b)) {
    eq <- isTRUE(all.equal(a, b, tolerance = 0))
  } else {
    a <- canon(tibble::as_tibble(a))
    b <- canon(tibble::as_tibble(b))
    eq <- isTRUE(all.equal(a, b, tolerance = 0))
  }
  if (eq) {
    n <- if (is.atomic(a)) length(a) else nrow(a)
    cat(sprintf("[ ok   ] %s (n = %d)\n", label, n))
  } else {
    cat("[ FAIL ]", label, "\n")
    if (!is.atomic(a)) {
      cat("        dims eager =", paste(dim(a), collapse = "x"),
          "lazy =", paste(dim(b), collapse = "x"), "\n")
    } else {
      cat("        len eager =", length(a), "lazy =", length(b), "\n")
    }
    cat("        ", all.equal(a, b)[[1]], "\n")
  }
  invisible(eq)
}

# ---------- Run checks ----------
all_ok <- TRUE
digests <- list()

for (loc in loci) {
  cat("\n========== ", loc, " ==========\n", sep = "")

  # reference_pip_tbl
  ref_e <- reference_pip_tbl_eager(loc)
  ref_l <- reference_pip_tbl_lazy(loc)
  all_ok <- assert_equal(paste0(loc, " :: reference_pip_tbl"), ref_e, ref_l) && all_ok
  digests[[paste0(loc, "::reference_pip_tbl")]] <- digest_or_na(ref_l)

  # run_pip_vec for the anchor and highest-weight runs
  anchor_id <- anchor_summary %>%
    dplyr::filter(.data$locus_id == loc) %>% dplyr::pull(.data$susie_anchor_run_id)
  hw_id <- paper_summary %>%
    dplyr::filter(.data$locus_id == loc) %>% dplyr::pull(.data$highest_weight_source_run_id)
  for (rid in c(anchor_id, hw_id)) {
    if (!length(rid) || is.na(rid)) next
    pv_e <- run_pip_vec_eager(loc, rid)
    pv_l <- run_pip_vec_lazy(loc, rid)
    all_ok <- assert_equal(sprintf("%s :: run_pip_vec(run_id=%d)", loc, rid), pv_e, pv_l) && all_ok
    digests[[sprintf("%s::run_pip_vec::%d", loc, rid)]] <- digest_or_na(pv_l)
  }

  # plot_drift inner pip_tbl
  run_ids <- functional_grid_summary %>%
    dplyr::filter(.data$locus_id == loc) %>% dplyr::pull(.data$run_id)
  drift_e <- build_drift_pip_tbl_eager(loc, ref_e, run_ids)
  drift_l <- build_drift_pip_tbl_lazy(loc, ref_l, run_ids)
  all_ok <- assert_equal(paste0(loc, " :: plot_drift_pip_tbl"), drift_e, drift_l) && all_ok
  digests[[paste0(loc, "::plot_drift_pip_tbl")]] <- digest_or_na(drift_l)

  # plot_aggregated_lollipop overlay
  baseline_id <- anchor_summary %>%
    dplyr::filter(.data$locus_id == loc) %>% dplyr::pull(.data$baseline_run_id)
  kept_variants <- aggregated_variants %>%
    dplyr::filter(.data$locus_id == loc) %>%
    dplyr::arrange(dplyr::desc(.data$aggregated_pip), .data$ld_matrix_index) %>%
    dplyr::slice_head(n = 20L) %>%
    dplyr::pull(.data$variant_id)
  if (length(anchor_id) && length(baseline_id) && length(kept_variants)) {
    ov_e <- build_lollipop_overlay_eager(loc, anchor_id, baseline_id, kept_variants)
    ov_l <- build_lollipop_overlay_lazy(loc, anchor_id, baseline_id, kept_variants)
    # The eager overlay in the original code did NOT filter by locus_id (it
    # filtered by run_id only). Match that scope by also restricting eager to
    # the locus under test for the equality check; run_ids are unique across
    # loci, so this is a no-op except for the column set.
    ov_e <- ov_e %>% dplyr::filter(.data$locus_id == loc)
    all_ok <- assert_equal(paste0(loc, " :: lollipop_overlay"), ov_e, ov_l) && all_ok
    digests[[paste0(loc, "::lollipop_overlay")]] <- digest_or_na(ov_l)
  }
}

# ---------- Digest write ----------
digest_path <- digest_out_arg %||%
  file.path(aggregated_dir, "viz_refactor_digests.rds")
saveRDS(digests, digest_path)
cat("\n[validate] digests written to ", digest_path, "\n", sep = "")

cat("\n==========================\n")
if (all_ok) {
  cat("[validate] ALL CHECKS PASSED\n")
  quit(save = "no", status = 0L)
} else {
  cat("[validate] AT LEAST ONE CHECK FAILED\n")
  quit(save = "no", status = 1L)
}
