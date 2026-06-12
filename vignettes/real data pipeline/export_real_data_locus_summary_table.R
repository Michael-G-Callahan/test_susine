#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(purrr)
})

argv <- commandArgs(trailingOnly = TRUE)

arg_value <- function(flag, default = NULL) {
  hit <- which(argv == flag)
  if (!length(hit) || hit[[1]] == length(argv)) return(default)
  argv[[hit[[1]] + 1L]]
}

arg_flag <- function(flag) {
  flag %in% argv
}

script_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
script_path_candidates <- c(
  if (length(script_arg)) sub("^--file=", "", script_arg[[1]]) else character(0),
  file.path("vignettes", "real data pipeline", "export_real_data_locus_summary_table.R")
)
script_path_candidates <- unique(script_path_candidates[nzchar(script_path_candidates)])
script_path_hit <- script_path_candidates[file.exists(script_path_candidates)]
script_path <- if (length(script_path_hit)) {
  normalizePath(script_path_hit[[1]], winslash = "/", mustWork = TRUE)
} else {
  normalizePath(script_path_candidates[[length(script_path_candidates)]],
                winslash = "/", mustWork = FALSE)
}
repo_root <- if (file.exists(file.path(getwd(), "DESCRIPTION")) &&
                 dir.exists(file.path(getwd(), "vignettes"))) {
  normalizePath(getwd(), winslash = "/", mustWork = TRUE)
} else {
  normalizePath(file.path(dirname(script_path), "..", ".."),
                winslash = "/", mustWork = TRUE)
}

job_name <- arg_value("--job-name", "real_data_ensemble_geometric_n20")
parent_job_id <- arg_value("--parent-job-id", "52906940")
aggregated_dir <- arg_value(
  "--aggregated-dir",
  file.path(repo_root, "output", "slurm_output", job_name, parent_job_id, "aggregated")
)
paper_dir <- arg_value(
  "--paper-dir",
  file.path(repo_root, "..", "..", "Writings", "plots", "real_data_case_study")
)
out_dir <- arg_value("--out-dir", paper_dir)

aggregated_dir <- normalizePath(aggregated_dir, winslash = "/", mustWork = FALSE)
paper_dir <- normalizePath(paper_dir, winslash = "/", mustWork = FALSE)
out_dir <- normalizePath(out_dir, winslash = "/", mustWork = FALSE)

summary_candidates <- c(
  file.path(aggregated_dir, "paper_real_data_ensemble_summary.csv"),
  file.path(paper_dir, "paper_real_data_ensemble_summary.csv")
)
summary_hits <- summary_candidates[file.exists(summary_candidates)]
summary_csv <- if (length(summary_hits)) summary_hits[[1]] else NA_character_
if (is.na(summary_csv) || !nzchar(summary_csv)) {
  stop("Could not find paper_real_data_ensemble_summary.csv in aggregated_dir or paper_dir.")
}

paper_summary <- readr::read_csv(summary_csv, show_col_types = FALSE)

# Corrected local-genetic-variance estimand (E[var(Xb)|y]/var(y)) from
# inst/scripts/recompute_real_data_hg2.R. Maps the SuSiE anchor by its run id and
# the SuSiNE ensemble by locus; falls back to the posterior-mean pve_* columns
# below when the recompute outputs are absent. See
# refs/decisions/heritability_estimand_decision_2026-06-09.md (sec 5.5).
hg2_corrected_dir <- file.path(dirname(aggregated_dir), "hg2_corrected")
.read_if <- function(fname) {
  p <- file.path(hg2_corrected_dir, fname)
  if (file.exists(p)) readr::read_csv(p, show_col_types = FALSE) else NULL
}
hg2_by_run <- .read_if("hg2_corrected_by_run.csv")
hg2_by_locus <- .read_if("hg2_corrected_by_locus.csv")

corrected_pve_tbl <- paper_summary %>%
  dplyr::distinct(.data$locus_id, .data$susie_anchor_run_id)
corrected_pve_tbl$pve_susie_anchor_expected <- if (!is.null(hg2_by_run)) {
  hg2_by_run$hg2_expected_pve[match(as.integer(corrected_pve_tbl$susie_anchor_run_id),
                                    as.integer(hg2_by_run$run_id))]
} else NA_real_
corrected_pve_tbl$pve_susine_ensemble_expected <- if (!is.null(hg2_by_locus)) {
  hg2_by_locus$hg2_expected_pve_ensemble[match(as.character(corrected_pve_tbl$locus_id),
                                               as.character(hg2_by_locus$locus_id))]
} else NA_real_
corrected_pve_tbl <- dplyr::select(corrected_pve_tbl, -dplyr::any_of("susie_anchor_run_id"))

weighted_annotation_tbl <- NULL
weighted_notes_csv <- file.path(paper_dir, "weighted_annotation_tbl_for_paper_notes.csv")
functional_grid_csv <- file.path(aggregated_dir, "functional_grid_summary.csv")
if (file.exists(weighted_notes_csv)) {
  weighted_annotation_tbl <- readr::read_csv(weighted_notes_csv, show_col_types = FALSE)
} else if (file.exists(functional_grid_csv)) {
  weighted_annotation_tbl <- readr::read_csv(functional_grid_csv, show_col_types = FALSE) %>%
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
      weighted_sigma_0_2_scalar =
        sum(.data$agg_weight_run * .data$sigma_0_2_scalar, na.rm = TRUE) /
        sum(.data$agg_weight_run, na.rm = TRUE),
      off_c0_weight = sum(.data$agg_weight_run[abs(.data$c_value) > 1e-12], na.rm = TRUE) /
        sum(.data$agg_weight_run, na.rm = TRUE),
      max_model_weight = max(.data$agg_weight_run, na.rm = TRUE),
      .groups = "drop"
    )
} else {
  weighted_annotation_tbl <- tibble::tibble(
    locus_id = character(),
    weighted_c_value = numeric(),
    off_c0_weight = numeric()
  )
}

locus_label <- function(locus_id, gene_name = NULL) {
  out <- if (!is.null(gene_name)) as.character(gene_name) else rep(NA_character_, length(locus_id))
  missing <- is.na(out) | !nzchar(out)
  out[missing] <- toupper(sub("_chr.*$", "", as.character(locus_id[missing])))
  out
}

locus_chromosome <- function(locus_id) {
  x <- as.character(locus_id)
  dplyr::if_else(
    grepl("_chr[^_]+_", x),
    sub("^.*_chr([^_]+)_.*$", "chr\\1", x),
    NA_character_
  )
}

compute_pip_l2_tbl <- function(paper_summary, aggregated_dir) {
  if ("pip_l2_susie_anchor_vs_susine_ensemble" %in% names(paper_summary)) {
    return(paper_summary %>%
      dplyr::transmute(
        .data$locus_id,
        pip_l2_susie_anchor_vs_susine_ensemble =
          .data$pip_l2_susie_anchor_vs_susine_ensemble
      ))
  }

  variant_dataset_dir <- file.path(aggregated_dir, "variant_posteriors_dataset")
  aggregated_variant_file <- file.path(aggregated_dir, "aggregated_variant_pips_cluster_weight.parquet")
  if (!dir.exists(variant_dataset_dir) || !file.exists(aggregated_variant_file)) {
    warning("PIP L2 could not be computed: full parquet outputs are missing.")
    return(tibble::tibble(
      locus_id = paper_summary$locus_id,
      pip_l2_susie_anchor_vs_susine_ensemble = NA_real_
    ))
  }
  if (!requireNamespace("arrow", quietly = TRUE)) {
    warning("PIP L2 could not be computed: R package 'arrow' is unavailable.")
    return(tibble::tibble(
      locus_id = paper_summary$locus_id,
      pip_l2_susie_anchor_vs_susine_ensemble = NA_real_
    ))
  }

  variant_posteriors <- arrow::open_dataset(variant_dataset_dir)
  aggregated_variants <- arrow::read_parquet(aggregated_variant_file) %>%
    tibble::as_tibble()

  purrr::map_dfr(seq_len(nrow(paper_summary)), function(i) {
    loc <- paper_summary$locus_id[[i]]
    anchor_run_id <- paper_summary$susie_anchor_run_id[[i]]
    if (is.na(anchor_run_id)) {
      return(tibble::tibble(
        locus_id = loc,
        pip_l2_susie_anchor_vs_susine_ensemble = NA_real_
      ))
    }
    anchor_pip <- variant_posteriors %>%
      dplyr::filter(.data$locus_id == !!loc,
                    .data$run_id == !!anchor_run_id) %>%
      dplyr::select(.data$variant_id, susie_anchor_pip = .data$pip) %>%
      dplyr::collect()
    ensemble_pip <- aggregated_variants %>%
      dplyr::filter(.data$locus_id == !!loc) %>%
      dplyr::select(.data$variant_id, susine_ensemble_pip = .data$aggregated_pip)
    pip_joined <- dplyr::inner_join(anchor_pip, ensemble_pip, by = "variant_id")
    pip_l2 <- if (nrow(pip_joined)) {
      sqrt(sum((pip_joined$susine_ensemble_pip - pip_joined$susie_anchor_pip)^2,
               na.rm = TRUE))
    } else {
      NA_real_
    }
    tibble::tibble(
      locus_id = loc,
      pip_l2_susie_anchor_vs_susine_ensemble = pip_l2
    )
  })
}

pip_l2_tbl <- compute_pip_l2_tbl(paper_summary, aggregated_dir)

locus_summary_tbl <- paper_summary %>%
  dplyr::select(-dplyr::any_of("pip_l2_susie_anchor_vs_susine_ensemble")) %>%
  dplyr::left_join(
    weighted_annotation_tbl %>%
      dplyr::select(locus_id, weighted_c_value,
                    weighted_off_c0_weight = off_c0_weight),
    by = "locus_id"
  ) %>%
  dplyr::left_join(pip_l2_tbl, by = "locus_id") %>%
  dplyr::left_join(corrected_pve_tbl, by = "locus_id") %>%
  dplyr::transmute(
    locus_id = .data$locus_id,
    locus_name = locus_label(.data$locus_id, .data$gene_name),
    chromosome = locus_chromosome(.data$locus_id),
    n_susie_pip_gt_05 = .data$n_pip_gt_05_susie_anchor,
    n_susine_pip_gt_05 = .data$n_pip_gt_05_susine_ensemble,
    n_agreed_pip_gt_05 = .data$n_pip_gt_05_susie_anchor_and_susine_ensemble,
    # Corrected estimand when available, else posterior-mean point estimate.
    delta_pve_susine_minus_susie =
      dplyr::coalesce(.data$pve_susine_ensemble_expected, .data$pve_susine_ensemble) -
      dplyr::coalesce(.data$pve_susie_anchor_expected, .data$pve_susie_anchor),
    delta_elbo_warm_refit_vs_susie = .data$delta_elbo_refit_vs_susie_anchor,
    total_off_c0_weight = dplyr::coalesce(
      .data$ensemble_agg_weight_off_c0_total,
      .data$weighted_off_c0_weight
    ),
    elbo_weighted_annotation_scale_c = .data$weighted_c_value,
    pip_l2_susine_ensemble_vs_susie = .data$pip_l2_susie_anchor_vs_susine_ensemble
  ) %>%
  dplyr::arrange(dplyr::desc(.data$pip_l2_susine_ensemble_vs_susie),
                 .data$locus_name)

pip_gt05_totals <- locus_summary_tbl %>%
  dplyr::summarise(
    n_loci = dplyr::n(),
    total_susie_pip_gt_05 = sum(.data$n_susie_pip_gt_05, na.rm = TRUE),
    total_susine_pip_gt_05 = sum(.data$n_susine_pip_gt_05, na.rm = TRUE),
    total_agreed_pip_gt_05 = sum(.data$n_agreed_pip_gt_05, na.rm = TRUE),
    .groups = "drop"
  )

if (!dir.exists(out_dir)) {
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
}

locus_out <- file.path(out_dir, "paper_real_data_locus_summary_table.csv")
totals_out <- file.path(out_dir, "paper_real_data_pip_gt05_totals.csv")
readr::write_csv(locus_summary_tbl, locus_out)
readr::write_csv(pip_gt05_totals, totals_out)

cat("Wrote:\n")
cat("  ", locus_out, "\n", sep = "")
cat("  ", totals_out, "\n", sep = "")
cat(sprintf(
  "Totals across %d loci: SuSiE=%d, SuSiNE=%d, agreed=%d PIP > 0.5\n",
  pip_gt05_totals$n_loci,
  pip_gt05_totals$total_susie_pip_gt_05,
  pip_gt05_totals$total_susine_pip_gt_05,
  pip_gt05_totals$total_agreed_pip_gt_05
))

if (arg_flag("--print")) {
  print(locus_summary_tbl, n = Inf)
  print(pip_gt05_totals)
}
