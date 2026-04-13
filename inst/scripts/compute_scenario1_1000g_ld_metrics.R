#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(readr)
  library(Rsamtools)
  library(IRanges)
  library(GenomicRanges)
})

find_pkg_root <- function() {
  if (file.exists("DESCRIPTION")) {
    return(normalizePath(".", winslash = "/", mustWork = TRUE))
  }

  file_arg <- commandArgs(FALSE)
  file_arg <- file_arg[grepl("^--file=", file_arg)]
  if (!length(file_arg)) {
    stop("Unable to determine package root.")
  }

  script_path <- normalizePath(sub("^--file=", "", file_arg[[1]]),
                               winslash = "/",
                               mustWork = TRUE)
  normalizePath(dirname(dirname(dirname(script_path))),
                winslash = "/",
                mustWork = TRUE)
}

pkg_root <- find_pkg_root()
data_dir <- file.path(pkg_root, "data", "sampled_simulated_genotypes")
cache_dir <- file.path(pkg_root, "data", "reference_cache")
dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)

summary_path <- file.path(data_dir, "scenario_sampling_summary.csv")
output_csv <- file.path(data_dir, "scenario_1_1000g_ld_metrics.csv")
selected_csv <- file.path(data_dir, "scenario_1_1000g_ld_representative_genes.csv")
plot_path <- file.path(data_dir, "scenario_1_1000g_ld_metrics_scatter.png")

panel_url <- paste0(
  "https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/",
  "integrated_call_samples_v3.20130502.ALL.panel"
)
panel_cache <- file.path(cache_dir, "integrated_call_samples_v3.20130502.ALL.panel")

vcf_url_for_chrom <- function(chrom) {
  sprintf(
    paste0(
      "https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/",
      "data_collections/1000_genomes_project/release/",
      "20190312_biallelic_SNV_and_INDEL/",
      "ALL.chr%s.shapeit2_integrated_snvindels_v2a_27022019.GRCh38.phased.vcf.gz"
    ),
    chrom
  )
}

rescale01 <- function(x) {
  rng <- range(x, na.rm = TRUE)
  if (!all(is.finite(rng)) || diff(rng) == 0) {
    return(rep(0.5, length(x)))
  }
  (x - rng[[1]]) / diff(rng)
}

mid_energy_M1 <- function(R) {
  A <- abs(R)
  diag(A) <- 0
  ut <- A[upper.tri(A)]
  mean(ut * (1 - ut)) * 2
}

ld_abs_alpha_thr <- function(R, alpha = 4, threshold = 0.5) {
  abs_r <- abs(R[upper.tri(R)])
  abs_r_thr <- ifelse(abs_r >= threshold, abs_r, 0)
  mean(abs_r_thr^alpha, na.rm = TRUE)
}

gt_to_dosage <- function(gt_vec) {
  out <- rep(NA_real_, length(gt_vec))
  out[gt_vec %in% c("0|0", "0/0")] <- 0
  out[gt_vec %in% c("0|1", "1|0", "0/1", "1/0")] <- 1
  out[gt_vec %in% c("1|1", "1/1")] <- 2
  out
}

ensure_panel <- function() {
  if (!file.exists(panel_cache)) {
    download.file(panel_url, panel_cache, mode = "wb", quiet = TRUE)
  }
  readr::read_tsv(
    panel_cache,
    show_col_types = FALSE,
    col_names = c("sample", "pop", "super_pop", "gender")
  )
}

read_vcf_sample_names <- function(vcf_url) {
  con <- gzcon(url(vcf_url, open = "rb"))
  on.exit(close(con), add = TRUE)
  hdr <- readLines(con, n = 500)
  chrom_line <- hdr[grepl("^#CHROM", hdr)]
  if (!length(chrom_line)) {
    stop("Unable to find #CHROM header in ", vcf_url)
  }
  strsplit(chrom_line[[1]], "\t", fixed = TRUE)[[1]][-(1:9)]
}

tabix_cache <- new.env(parent = emptyenv())

get_tabix_file <- function(chrom) {
  key <- as.character(chrom)
  if (!exists(key, envir = tabix_cache, inherits = FALSE)) {
    tf <- TabixFile(vcf_url_for_chrom(chrom))
    open(tf)
    assign(key, tf, envir = tabix_cache)
  }
  get(key, envir = tabix_cache, inherits = FALSE)
}

close_tabix_cache <- function() {
  for (nm in ls(tabix_cache, all.names = TRUE)) {
    tf <- get(nm, envir = tabix_cache, inherits = FALSE)
    try(close(tf), silent = TRUE)
  }
}

make_target_key <- function(pos, allele1, allele2) {
  alleles <- sort(toupper(c(allele1, allele2)))
  paste(pos, alleles[[1]], alleles[[2]], sep = "|")
}

compute_reference_metrics_one <- function(gene_row, eur_field_idx) {
  manifest <- readr::read_tsv(gene_row$manifest_path, show_col_types = FALSE) %>%
    dplyr::transmute(
      manifest_index = dplyr::row_number(),
      snp_id = .data$snp_id,
      bp = as.integer(.data$bp),
      allele_a = toupper(.data$allele_ref),
      allele_b = toupper(.data$allele_alt),
      match_key = make_target_key(.data$bp, .data$allele_ref, .data$allele_alt)
    )

  requested_snps <- nrow(manifest)
  if (requested_snps < 2L) {
    return(tibble::tibble(
      requested_snps = requested_snps,
      matched_snps = 0L,
      match_rate = 0,
      M1 = NA_real_,
      ld_abs_alpha_4_thr_0p5 = NA_real_,
      status = "too_few_requested_snps"
    ))
  }

  match_idx <- setNames(manifest$manifest_index, manifest$match_key)
  chr <- as.character(gene_row$chromosome)
  region <- GenomicRanges::GRanges(
    seqnames = chr,
    ranges = IRanges::IRanges(min(manifest$bp), max(manifest$bp))
  )

  tf <- get_tabix_file(chr)
  records <- scanTabix(tf, param = region)[[1]]

  G <- matrix(NA_real_, nrow = requested_snps, ncol = length(eur_field_idx))

  for (rec in records) {
    fields <- strsplit(rec, "\t", fixed = TRUE)[[1]]
    if (length(fields) < max(9L + eur_field_idx)) {
      next
    }

    alt <- strsplit(fields[[5]], ",", fixed = TRUE)[[1]][[1]]
    key <- make_target_key(as.integer(fields[[2]]), fields[[4]], alt)
    idx <- match_idx[[key]]
    if (is.null(idx) || !is.na(G[idx, 1])) {
      next
    }

    gt_fields <- fields[9L + eur_field_idx]
    G[idx, ] <- gt_to_dosage(gt_fields)
  }

  matched_rows <- which(rowSums(!is.na(G)) > 0)
  matched_snps <- length(matched_rows)
  if (matched_snps < 2L) {
    return(tibble::tibble(
      requested_snps = requested_snps,
      matched_snps = matched_snps,
      match_rate = matched_snps / requested_snps,
      M1 = NA_real_,
      ld_abs_alpha_4_thr_0p5 = NA_real_,
      status = "too_few_matches"
    ))
  }

  G <- G[matched_rows, , drop = FALSE]
  row_means <- rowMeans(G, na.rm = TRUE)
  for (i in seq_len(nrow(G))) {
    na_idx <- is.na(G[i, ])
    if (any(na_idx)) {
      G[i, na_idx] <- row_means[[i]]
    }
  }

  G_centered <- G - rowMeans(G)
  row_sd <- apply(G_centered, 1, stats::sd)
  keep <- is.finite(row_sd) & row_sd > 0
  G_centered <- G_centered[keep, , drop = FALSE]
  row_sd <- row_sd[keep]

  if (nrow(G_centered) < 2L) {
    return(tibble::tibble(
      requested_snps = requested_snps,
      matched_snps = matched_snps,
      match_rate = matched_snps / requested_snps,
      M1 = NA_real_,
      ld_abs_alpha_4_thr_0p5 = NA_real_,
      status = "too_few_variable_snps"
    ))
  }

  X <- G_centered / row_sd
  R <- tcrossprod(X) / max(ncol(X) - 1L, 1L)
  R[is.na(R)] <- 0
  diag(R) <- 1

  tibble::tibble(
    requested_snps = requested_snps,
    matched_snps = matched_snps,
    match_rate = matched_snps / requested_snps,
    M1 = mid_energy_M1(R),
    ld_abs_alpha_4_thr_0p5 = ld_abs_alpha_thr(R, alpha = 4, threshold = 0.5),
    status = "ok"
  )
}

select_representative_genes <- function(metrics_tbl) {
  valid <- metrics_tbl %>%
    dplyr::filter(is.finite(.data$M1), is.finite(.data$ld_abs_alpha_4_thr_0p5))

  if (nrow(valid) < 4L) {
    return(valid)
  }

  valid <- valid %>%
    dplyr::mutate(
      x01 = rescale01(.data$M1),
      y01 = rescale01(.data$ld_abs_alpha_4_thr_0p5)
    )

  targets <- tibble::tibble(
    rep_id = paste0("G", 1:4),
    target_x = c(0.2, 0.2, 0.8, 0.8),
    target_y = c(0.2, 0.8, 0.2, 0.8)
  )

  picked <- integer(0)
  chosen <- vector("list", nrow(targets))
  for (i in seq_len(nrow(targets))) {
    d2 <- (valid$x01 - targets$target_x[[i]])^2 +
      (valid$y01 - targets$target_y[[i]])^2
    if (length(picked)) {
      d2[picked] <- Inf
    }
    idx <- which.min(d2)
    picked <- c(picked, idx)
    chosen[[i]] <- dplyr::bind_cols(targets[i, ], valid[idx, ])
  }

  dplyr::bind_rows(chosen) %>%
    dplyr::select(
      .data$rep_id,
      .data$dataset_label,
      .data$chromosome,
      .data$tss_bp,
      .data$matrix_id,
      .data$requested_snps,
      .data$matched_snps,
      .data$match_rate,
      .data$M1,
      .data$ld_abs_alpha_4_thr_0p5,
      .data$target_x,
      .data$target_y
    )
}

stopifnot(file.exists(summary_path))

message("Loading summary and 1000G EUR sample metadata...")
summary_tbl <- readr::read_csv(summary_path, show_col_types = FALSE) %>%
  dplyr::filter(.data$scenario == "scenario_1") %>%
  dplyr::mutate(matrix_id = dplyr::row_number()) %>%
  dplyr::arrange(.data$matrix_id)

panel_tbl <- ensure_panel()
eur_samples <- panel_tbl %>%
  dplyr::filter(.data$super_pop == "EUR") %>%
  dplyr::pull(.data$sample)

sample_names <- read_vcf_sample_names(vcf_url_for_chrom(summary_tbl$chromosome[[1]]))
eur_field_idx <- which(sample_names %in% eur_samples)

if (!length(eur_field_idx)) {
  stop("No EUR samples from the panel were found in the 1000G VCF header.")
}

message("EUR samples used: ", length(eur_field_idx))

on.exit(close_tabix_cache(), add = TRUE)

metric_rows <- vector("list", nrow(summary_tbl))
for (i in seq_len(nrow(summary_tbl))) {
  gene_row <- summary_tbl[i, , drop = FALSE]
  message(sprintf(
    "[%d/%d] %s",
    i,
    nrow(summary_tbl),
    gene_row$gene[[1]]
  ))

  metric_rows[[i]] <- gene_row %>%
    dplyr::transmute(
      matrix_id = .data$matrix_id,
      dataset_label = .data$gene,
      chromosome = .data$chromosome,
      tss_bp = .data$tss_bp,
      scenario = .data$scenario,
      snp_set = .data$snp_set,
      participant_count = .data$participant_count,
      snps_post = .data$snps_post,
      matrix_path = .data$matrix_path,
      manifest_path = .data$manifest_path
    ) %>%
    dplyr::bind_cols(
      compute_reference_metrics_one(gene_row, eur_field_idx = eur_field_idx)
    )
}

metrics_tbl <- dplyr::bind_rows(metric_rows) %>%
  dplyr::arrange(.data$matrix_id)

selected_tbl <- select_representative_genes(metrics_tbl)

plot_tbl <- metrics_tbl %>%
  dplyr::filter(is.finite(.data$M1), is.finite(.data$ld_abs_alpha_4_thr_0p5))

plot_obj <- ggplot(plot_tbl, aes(x = .data$M1, y = .data$ld_abs_alpha_4_thr_0p5)) +
  geom_point(color = "#4C78A8", alpha = 0.75, size = 2) +
  geom_point(
    data = selected_tbl,
    color = "#D62728",
    fill = "#D62728",
    size = 3
  ) +
  geom_text(
    data = selected_tbl,
    aes(label = .data$rep_id),
    nudge_x = 0.0025,
    nudge_y = 0.00025,
    size = 3.2,
    color = "#D62728",
    check_overlap = TRUE
  ) +
  labs(
    title = "Scenario 1 Gene-Level 1000G LD Metrics",
    subtitle = "EUR 1000 Genomes reference LD using sampled-manifest SNP lists",
    x = "M1",
    y = "ld_abs_alpha_4_thr_0p5"
  ) +
  theme_minimal(base_size = 12)

readr::write_csv(metrics_tbl, output_csv)
readr::write_csv(selected_tbl, selected_csv)
ggplot2::ggsave(plot_path, plot_obj, width = 8, height = 6, dpi = 150, bg = "white")

message("Wrote metrics: ", output_csv)
message("Wrote representatives: ", selected_csv)
message("Wrote plot: ", plot_path)
