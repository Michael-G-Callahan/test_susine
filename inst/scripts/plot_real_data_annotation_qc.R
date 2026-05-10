#!/usr/bin/env Rscript

parse_args <- function(args) {
  out <- list()
  i <- 1L
  while (i <= length(args)) {
    key <- args[[i]]
    if (!startsWith(key, "--")) {
      stop("Unexpected positional argument: ", key)
    }
    key <- sub("^--", "", key)
    if (i == length(args) || startsWith(args[[i + 1L]], "--")) {
      out[[key]] <- TRUE
      i <- i + 1L
    } else {
      out[[key]] <- args[[i + 1L]]
      i <- i + 2L
    }
  }
  out
}

`%||%` <- function(x, y) {
  if (is.null(x) || length(x) == 0L || all(is.na(x))) y else x
}

repo_root <- function(path = getwd()) {
  current <- normalizePath(path, winslash = "/", mustWork = TRUE)
  repeat {
    if (file.exists(file.path(current, "DESCRIPTION")) || file.exists(file.path(current, ".here"))) {
      return(current)
    }
    parent <- dirname(current)
    if (identical(parent, current)) return(normalizePath(path, winslash = "/", mustWork = TRUE))
    current <- parent
  }
}

resolve_path <- function(path, root) {
  if (is.na(path) || !nzchar(path)) return(NA_character_)
  if (grepl("^/", path)) return(normalizePath(path, winslash = "/", mustWork = TRUE))
  normalizePath(file.path(root, path), winslash = "/", mustWork = TRUE)
}

zscale <- function(x) {
  x <- as.numeric(x)
  s <- stats::sd(x, na.rm = TRUE)
  if (!is.finite(s) || s <= 0) return(rep(NA_real_, length(x)))
  (x - mean(x, na.rm = TRUE)) / s
}

rms <- function(x) {
  x <- as.numeric(x)
  sqrt(mean(x^2, na.rm = TRUE))
}

q_safe <- function(x, p) {
  x <- as.numeric(x)
  x <- x[is.finite(x)]
  if (!length(x)) return(NA_real_)
  unname(stats::quantile(x, probs = p, na.rm = TRUE))
}

read_manifest <- function(path) {
  if (!file.exists(path)) stop("Missing manifest: ", path)
  utils::read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
}

read_annotation <- function(path) {
  df <- utils::read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
  required <- c("variant_id", "annotation_a", "beta_hat_std")
  missing <- setdiff(required, names(df))
  if (length(missing)) {
    stop("Annotation file missing columns: ", paste(missing, collapse = ", "), " in ", path)
  }
  df
}

plot_overlay <- function(a, b, title, subtitle) {
  a <- as.numeric(a)
  b <- as.numeric(b)
  lim <- range(c(a, b), na.rm = TRUE)
  if (!all(is.finite(lim)) || diff(lim) == 0) lim <- c(-1, 1)
  hist(
    a,
    breaks = 50,
    freq = FALSE,
    col = grDevices::adjustcolor("#2f6c8f", alpha.f = 0.35),
    border = "white",
    xlim = lim,
    main = title,
    xlab = "",
    ylab = "Density"
  )
  hist(
    b,
    breaks = 50,
    freq = FALSE,
    col = grDevices::adjustcolor("#c84b31", alpha.f = 0.35),
    border = "white",
    add = TRUE
  )
  legend(
    "topright",
    legend = c("annotation_a", "beta_hat_std"),
    fill = grDevices::adjustcolor(c("#2f6c8f", "#c84b31"), alpha.f = 0.35),
    border = NA,
    bty = "n"
  )
  mtext(subtitle, side = 3, line = 0.25, cex = 0.75)
}

plot_scatter <- function(a, b, title, subtitle) {
  plot(
    a,
    b,
    pch = 16,
    cex = 0.35,
    col = grDevices::adjustcolor("#2b2b2b", alpha.f = 0.35),
    xlab = "annotation_a",
    ylab = "beta_hat_std",
    main = title
  )
  abline(h = 0, v = 0, col = "#999999", lty = 3)
  fit <- try(stats::lm(b ~ a), silent = TRUE)
  if (!inherits(fit, "try-error")) abline(fit, col = "#2f6c8f", lwd = 2)
  mtext(subtitle, side = 3, line = 0.25, cex = 0.75)
}

main <- function() {
  args <- parse_args(commandArgs(trailingOnly = TRUE))
  root <- repo_root(getwd())
  manifest_path <- args[["manifest"]] %||%
    file.path(root, "data", "real_case_studies", "geometric_n20_loci", "locus_manifest.csv")
  out_dir <- args[["out-dir"]] %||%
    file.path(root, "output", "real_data_annotation_qc", "geometric_n20")
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  manifest <- read_manifest(manifest_path)
  if (!nrow(manifest)) stop("Manifest has no rows: ", manifest_path)
  if (!"annotation_path" %in% names(manifest)) stop("Manifest missing annotation_path column.")

  pdf_path <- file.path(out_dir, "annotation_beta_hat_overlay.pdf")
  summary_path <- file.path(out_dir, "annotation_beta_hat_summary.csv")

  rows <- vector("list", nrow(manifest))
  grDevices::pdf(pdf_path, width = 11, height = 8.5, onefile = TRUE)
  on.exit(grDevices::dev.off(), add = TRUE)

  for (i in seq_len(nrow(manifest))) {
    locus_id <- manifest$locus_id[[i]]
    gene_name <- manifest$gene_name[[i]]
    annotation_path <- resolve_path(manifest$annotation_path[[i]], root)
    df <- read_annotation(annotation_path)

    a <- as.numeric(df$annotation_a)
    b <- as.numeric(df$beta_hat_std)
    ok <- is.finite(a) & is.finite(b)
    a_z <- zscale(a)
    b_z <- zscale(b)
    pearson <- if (sum(ok) >= 3L) suppressWarnings(stats::cor(a[ok], b[ok], method = "pearson")) else NA_real_
    spearman <- if (sum(ok) >= 3L) suppressWarnings(stats::cor(a[ok], b[ok], method = "spearman")) else NA_real_

    rows[[i]] <- data.frame(
      locus_id = locus_id,
      gene_name = gene_name,
      n_variants = nrow(df),
      annotation_mean = mean(a, na.rm = TRUE),
      annotation_sd = stats::sd(a, na.rm = TRUE),
      annotation_rms = rms(a),
      annotation_q01 = q_safe(a, 0.01),
      annotation_q05 = q_safe(a, 0.05),
      annotation_q50 = q_safe(a, 0.50),
      annotation_q95 = q_safe(a, 0.95),
      annotation_q99 = q_safe(a, 0.99),
      beta_hat_mean = mean(b, na.rm = TRUE),
      beta_hat_sd = stats::sd(b, na.rm = TRUE),
      beta_hat_rms = rms(b),
      beta_hat_q01 = q_safe(b, 0.01),
      beta_hat_q05 = q_safe(b, 0.05),
      beta_hat_q50 = q_safe(b, 0.50),
      beta_hat_q95 = q_safe(b, 0.95),
      beta_hat_q99 = q_safe(b, 0.99),
      pearson_annotation_beta_hat = pearson,
      spearman_annotation_beta_hat = spearman,
      stringsAsFactors = FALSE
    )

    old_par <- par(no.readonly = TRUE)
    par(mfrow = c(2, 2), mar = c(4.2, 4.2, 3.6, 1.2), oma = c(0, 0, 2, 0))
    subtitle <- sprintf("n=%s | cor=%.3f | rho=%.3f", format(nrow(df), big.mark = ","), pearson, spearman)
    plot_overlay(a, b, "Raw distributions", subtitle)
    plot_overlay(a_z, b_z, "Within-locus z-scaled distributions", subtitle)
    plot_scatter(a, b, "Raw annotation vs marginal beta-hat", subtitle)
    plot_scatter(a_z, b_z, "Z-scaled annotation vs marginal beta-hat", subtitle)
    mtext(sprintf("%s (%s)", gene_name, locus_id), outer = TRUE, cex = 1.1, font = 2)
    par(old_par)
  }

  summary_df <- do.call(rbind, rows)
  utils::write.csv(summary_df, summary_path, row.names = FALSE)
  message("Wrote: ", pdf_path)
  message("Wrote: ", summary_path)
}

main()
