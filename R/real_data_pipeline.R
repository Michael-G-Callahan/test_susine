# Real-data RSS ensemble pipeline -------------------------------------------

real_data_required_annotation_cols <- function() {
  c(
    "variant_id",
    "ld_matrix_index",
    "annotation_a",
    "beta_hat_std",
    "baseline_c_l",
    "var_y_hat_from_slope",
    "var_y_hat_from_se"
  )
}

real_data_study_root <- function(repo_root = ensure_repo_root(getwd())) {
  file.path(repo_root, "data", "real_case_studies", "geometric_n20_loci")
}

real_data_default_manifest_path <- function(repo_root = ensure_repo_root(getwd())) {
  file.path(real_data_study_root(repo_root), "locus_manifest.csv")
}

real_data_resolve_path <- function(path, repo_root) {
  if (is.na(path) || !nzchar(path)) {
    return(NA_character_)
  }
  if (grepl("^/", path)) {
    return(normalizePath(path, winslash = "/", mustWork = TRUE))
  }
  normalizePath(file.path(repo_root, path), winslash = "/", mustWork = TRUE)
}

real_data_parse_locus_id <- function(locus_id) {
  m <- regexec("^(.+?)_(chr[^_]+)_(.+)$", locus_id)
  parts <- regmatches(locus_id, m)[[1]]
  if (length(parts) != 4L) {
    stop("Unable to parse locus_id: ", locus_id)
  }
  list(
    gene_key = parts[[2]],
    chrom = parts[[3]],
    tissue = parts[[4]]
  )
}

real_data_source_paths <- function(source_repo_root, locus_id, gene_name) {
  gene_upper <- toupper(gene_name)
  list(
    z_scores = file.path(source_repo_root, "output", "z_score", locus_id,
                         paste0(gene_upper, "_GTEx_z_scores.csv")),
    ld_long = file.path(source_repo_root, "output", "ld", locus_id,
                        paste0(gene_upper, "_phase1_LD_R_long.parquet")),
    ld_order = file.path(source_repo_root, "output", "ld", locus_id,
                         paste0(gene_upper, "_LD_variant_order.tsv")),
    master_variants = file.path(source_repo_root, "output", "ld", locus_id,
                                paste0(gene_upper, "_phase1_master_variants.csv")),
    annotations = file.path(
      source_repo_root,
      "output", "susine_mu0",
      "geometric_n20",
      "per_locus_annotations",
      paste0(gene_upper, "_mu0_variant_annotations.csv")
    )
  )
}

real_data_check_required_files <- function(paths, locus_id) {
  missing <- names(paths)[!vapply(paths, file.exists, logical(1))]
  if (length(missing)) {
    stop(
      "Missing required real-data source file(s) for ", locus_id, ": ",
      paste(missing, collapse = ", ")
    )
  }
  invisible(paths)
}

real_data_first_or_default <- function(x, default = NA) {
  if (length(x) < 1L) {
    return(default)
  }
  x[[1]]
}

real_data_scalar_or_na <- function(x, default = NA) {
  x <- unique(stats::na.omit(x))
  if (!length(x)) {
    return(default)
  }
  x[[1]]
}

real_data_bind_csv_files <- function(files) {
  if (!length(files)) {
    return(tibble::tibble())
  }
  purrr::map_dfr(
    files,
    ~ readr::read_csv(.x, show_col_types = FALSE, progress = FALSE)
  )
}

sync_real_data_inputs <- function(
    source_repo_root = "/storage/work/mgc5166/Annotations/eQTL_annotations_for_susine",
    dest_root = real_data_study_root(),
    loci = NULL,
    source_mu0_name = "geometric_n20",
    source_annotation_summary = NULL
) {
  repo_root <- ensure_repo_root(getwd())
  source_repo_root <- normalizePath(source_repo_root, winslash = "/", mustWork = TRUE)
  dest_root <- normalizePath(dest_root, winslash = "/", mustWork = FALSE)
  ensure_dir(dest_root)

  mu0_summary_path <- file.path(
    source_repo_root,
    "output", "susine_mu0",
    source_mu0_name,
    "mu0_locus_summary.csv"
  )
  selection_summary_path <- source_annotation_summary %||% file.path(
    source_repo_root,
    "output", "annotation", "alphagenome",
    paste0(source_mu0_name, "_annotation_batch_summary.csv")
  )
  if (!file.exists(mu0_summary_path)) {
    stop("Missing mu0_locus_summary.csv at ", mu0_summary_path)
  }

  mu0_summary <- readr::read_csv(mu0_summary_path, show_col_types = FALSE)
  if (!is.null(loci)) {
    mu0_summary <- dplyr::filter(mu0_summary, .data$locus_id %in% loci)
  }
  if (!nrow(mu0_summary)) {
    stop("No loci selected for real-data sync.")
  }

  copied_rows <- purrr::map_dfr(seq_len(nrow(mu0_summary)), function(i) {
    row <- mu0_summary[i, , drop = FALSE]
    locus_id <- row$locus_id[[1]]
    gene_name <- row$gene_name[[1]]
    locus_dir <- file.path(dest_root, locus_id)
    ensure_dir(locus_dir)

    src <- real_data_source_paths(source_repo_root, locus_id, gene_name)
    src$annotations <- file.path(
      source_repo_root,
      "output", "susine_mu0",
      source_mu0_name,
      "per_locus_annotations",
      paste0(gene_name, "_mu0_variant_annotations.csv")
    )
    real_data_check_required_files(src, locus_id)

    dest_files <- list(
      z_scores = file.path(locus_dir, basename(src$z_scores)),
      ld_long = file.path(locus_dir, basename(src$ld_long)),
      ld_order = file.path(locus_dir, basename(src$ld_order)),
      master_variants = file.path(locus_dir, basename(src$master_variants)),
      annotations = file.path(locus_dir, basename(src$annotations))
    )

    ok <- file.copy(
      from = unname(unlist(src)),
      to = unname(unlist(dest_files)),
      overwrite = TRUE
    )
    if (!all(ok)) {
      failed <- names(dest_files)[!ok]
      stop("Failed to copy ", locus_id, " file(s): ", paste(failed, collapse = ", "))
    }

      tibble::tibble(
      locus_id = locus_id,
      gene_name = gene_name,
      z_scores_path = dest_files[["z_scores"]],
      ld_r_long_path = dest_files[["ld_long"]],
      ld_order_path = dest_files[["ld_order"]],
      master_variants_path = dest_files[["master_variants"]],
      annotation_path = dest_files[["annotations"]]
    )
  })

  provenance_files <- c(mu0_summary_path)
  if (file.exists(selection_summary_path)) {
    provenance_files <- c(provenance_files, selection_summary_path)
  }
  prov_ok <- file.copy(
    provenance_files,
    file.path(dest_root, basename(provenance_files)),
    overwrite = TRUE
  )
  if (!all(prov_ok)) {
    stop("Failed to copy one or more study-level provenance files.")
  }

  manifest <- build_real_data_manifest(dest_root = dest_root, repo_root = repo_root)
  list(
    dest_root = dest_root,
    manifest_path = file.path(dest_root, "locus_manifest.csv"),
    manifest = manifest,
    copied = copied_rows
  )
}

build_real_data_manifest <- function(
    dest_root = real_data_study_root(),
    repo_root = ensure_repo_root(getwd())
) {
  dest_root <- normalizePath(dest_root, winslash = "/", mustWork = TRUE)
  repo_root <- normalizePath(repo_root, winslash = "/", mustWork = TRUE)

  locus_dirs <- list.dirs(dest_root, recursive = FALSE, full.names = TRUE)
  if (!length(locus_dirs)) {
    stop("No locus directories found under ", dest_root)
  }

  manifest <- purrr::map_dfr(sort(locus_dirs), function(locus_dir) {
    locus_id <- basename(locus_dir)
    parsed <- real_data_parse_locus_id(locus_id)

    annotation_path <- list.files(
      locus_dir,
      pattern = "_mu0_variant_annotations\\.csv$",
      full.names = TRUE
    )
    z_scores_path <- list.files(
      locus_dir,
      pattern = "_GTEx_z_scores\\.csv$",
      full.names = TRUE
    )
    ld_long_path <- list.files(
      locus_dir,
      pattern = "_phase1_LD_R_long\\.parquet$",
      full.names = TRUE
    )
    ld_order_path <- list.files(
      locus_dir,
      pattern = "_LD_variant_order\\.tsv$",
      full.names = TRUE
    )
    master_variants_path <- list.files(
      locus_dir,
      pattern = "_phase1_master_variants\\.csv$",
      full.names = TRUE
    )

    paths <- c(
      annotation = real_data_first_or_default(annotation_path, NA_character_),
      z_scores = real_data_first_or_default(z_scores_path, NA_character_),
      ld_long = real_data_first_or_default(ld_long_path, NA_character_),
      ld_order = real_data_first_or_default(ld_order_path, NA_character_),
      master_variants = real_data_first_or_default(master_variants_path, NA_character_)
    )
    if (any(!file.exists(paths))) {
      missing <- names(paths)[!file.exists(paths)]
      stop("Missing copied real-data inputs in ", locus_id, ": ", paste(missing, collapse = ", "))
    }

    annotation_tbl <- readr::read_csv(paths[["annotation"]], show_col_types = FALSE)
    z_tbl <- readr::read_csv(paths[["z_scores"]], show_col_types = FALSE)

    tibble::tibble(
      locus_id = locus_id,
      gene_name = if ("gene_name" %in% names(annotation_tbl)) {
        as.character(real_data_scalar_or_na(annotation_tbl$gene_name, toupper(parsed$gene_key)))
      } else {
        toupper(parsed$gene_key)
      },
      gtex_tissue = parsed$tissue,
      gtex_chrom = parsed$chrom,
      n_variants = nrow(annotation_tbl),
      n_sample_median = stats::median(z_tbl$sample_size, na.rm = TRUE),
      baseline_c_l = as.numeric(real_data_scalar_or_na(annotation_tbl$baseline_c_l, NA_real_)),
      z_scores_path = relativize_path(paths[["z_scores"]], repo_root),
      ld_r_long_path = relativize_path(paths[["ld_long"]], repo_root),
      ld_order_path = relativize_path(paths[["ld_order"]], repo_root),
      master_variants_path = relativize_path(paths[["master_variants"]], repo_root),
      annotation_path = relativize_path(paths[["annotation"]], repo_root)
    )
  }) %>%
    dplyr::arrange(.data$locus_id) %>%
    dplyr::mutate(dataset_bundle_id = dplyr::row_number(), .before = 1L)

  manifest_path <- file.path(dest_root, "locus_manifest.csv")
  readr::write_csv(manifest, manifest_path)
  manifest
}

real_data_validate_annotation_columns <- function(annotation_tbl, locus_id) {
  required <- real_data_required_annotation_cols()
  missing <- setdiff(required, names(annotation_tbl))
  if (length(missing)) {
    stop(
      "Annotation file for ", locus_id, " is missing required columns: ",
      paste(missing, collapse = ", ")
    )
  }
  invisible(annotation_tbl)
}

real_data_validate_contiguous_indices <- function(idx, label) {
  idx <- as.integer(idx)
  expected <- seq.int(0L, length(idx) - 1L)
  if (!identical(idx, expected)) {
    stop(label, " indices must be contiguous 0:(p-1).")
  }
  invisible(idx)
}

real_data_reconstruct_ld_matrix <- function(ld_long_tbl, p, locus_id) {
  required <- c("snp_index_1", "snp_index_2", "r")
  missing <- setdiff(required, names(ld_long_tbl))
  if (length(missing)) {
    stop("LD long table for ", locus_id, " is missing columns: ", paste(missing, collapse = ", "))
  }

  idx1 <- as.integer(ld_long_tbl$snp_index_1)
  idx2 <- as.integer(ld_long_tbl$snp_index_2)
  r <- as.numeric(ld_long_tbl$r)

  if (any(idx1 < 0L | idx1 >= p | idx2 < 0L | idx2 >= p, na.rm = TRUE)) {
    stop("LD indices out of range for ", locus_id)
  }
  key <- paste(idx1, idx2, sep = "_")
  if (anyDuplicated(key)) {
    stop("Duplicate LD pairs detected for ", locus_id)
  }

  R <- diag(1, p)
  R[cbind(idx1 + 1L, idx2 + 1L)] <- r
  R[cbind(idx2 + 1L, idx1 + 1L)] <- r

  if (max(abs(R - t(R)), na.rm = TRUE) > 1e-8) {
    stop("Reconstructed LD matrix is not symmetric for ", locus_id)
  }
  if (max(abs(diag(R) - 1), na.rm = TRUE) > 1e-8) {
    stop("Reconstructed LD matrix does not have unit diagonal for ", locus_id)
  }
  R
}

#' @param apply_sample_size_filter Opt-in flag. If TRUE, drops variants with
#'   per-variant effective n < 0.5 * max(n) and subsets z, R, annotations
#'   together. Default FALSE: keep all variants that survived the upstream
#'   prep filter (`sample_size > 50` in `1_get_z_scores.ipynb`). Filtering on
#'   n excises exactly the low-power variants where the directional prior is
#'   most informative, so it is off by default; available for sensitivity
#'   analyses or for loci where the within-locus n spread is so wide that the
#'   single-n assumption in susie_rss becomes badly violated.
load_real_data_locus_bundle <- function(
    locus_id,
    manifest_path = real_data_default_manifest_path(),
    repo_root = ensure_repo_root(getwd()),
    apply_sample_size_filter = FALSE
) {
  repo_root <- normalizePath(repo_root, winslash = "/", mustWork = TRUE)
  manifest_path <- normalizePath(manifest_path, winslash = "/", mustWork = TRUE)
  manifest <- readr::read_csv(manifest_path, show_col_types = FALSE)
  row <- dplyr::filter(manifest, .data$locus_id == !!locus_id)
  if (nrow(row) != 1L) {
    stop("Expected exactly one manifest row for locus_id ", locus_id)
  }

  annotation_path <- real_data_resolve_path(row$annotation_path[[1]], repo_root)
  z_scores_path <- real_data_resolve_path(row$z_scores_path[[1]], repo_root)
  ld_long_path <- real_data_resolve_path(row$ld_r_long_path[[1]], repo_root)
  ld_order_path <- real_data_resolve_path(row$ld_order_path[[1]], repo_root)
  master_variants_path <- real_data_resolve_path(row$master_variants_path[[1]], repo_root)

  annotation_tbl <- readr::read_csv(annotation_path, show_col_types = FALSE) %>%
    dplyr::arrange(.data$ld_matrix_index)
  real_data_validate_annotation_columns(annotation_tbl, locus_id)
  real_data_validate_contiguous_indices(annotation_tbl$ld_matrix_index, "Annotation")

  ld_order_tbl <- readr::read_tsv(ld_order_path, show_col_types = FALSE) %>%
    dplyr::arrange(.data$index)
  real_data_validate_contiguous_indices(ld_order_tbl$index, "LD order")

  if (!identical(as.character(ld_order_tbl$id), as.character(annotation_tbl$variant_id))) {
    stop("LD order and annotation variant_id values do not match exactly for ", locus_id)
  }

  z_tbl <- readr::read_csv(z_scores_path, show_col_types = FALSE)
  if (anyDuplicated(z_tbl$variant_id)) {
    dupes <- unique(z_tbl$variant_id[duplicated(z_tbl$variant_id)])
    stop(
      "Duplicate variant_id rows found in z-score file for ", locus_id, ": ",
      paste(head(dupes, 5L), collapse = ", ")
    )
  }
  z_match <- match(ld_order_tbl$id, z_tbl$variant_id)
  if (any(is.na(z_match))) {
    missing_variants <- ld_order_tbl$id[is.na(z_match)]
    stop(
      "Missing z-score rows for ", locus_id, ": ",
      paste(head(missing_variants, 5L), collapse = ", ")
    )
  }
  z_aligned <- z_tbl[z_match, , drop = FALSE]

  master_tbl <- readr::read_csv(master_variants_path, show_col_types = FALSE) %>%
    dplyr::arrange(.data$ld_matrix_index)
  real_data_validate_contiguous_indices(master_tbl$ld_matrix_index, "Master variants")
  if (!identical(as.character(master_tbl$variant_id), as.character(annotation_tbl$variant_id))) {
    stop("Master variants and annotation variant_id values do not match exactly for ", locus_id)
  }

  ld_long_tbl <- arrow::read_parquet(ld_long_path) %>%
    tibble::as_tibble()
  p <- nrow(ld_order_tbl)
  R <- real_data_reconstruct_ld_matrix(ld_long_tbl, p = p, locus_id = locus_id)

  variant_map <- tibble::tibble(
    variant_id = as.character(annotation_tbl$variant_id),
    ld_matrix_index = as.integer(annotation_tbl$ld_matrix_index),
    chrom = as.character(ld_order_tbl$chrom %||% NA_character_),
    pos = as.integer(ld_order_tbl$pos %||% NA_integer_),
    ref = as.character(ld_order_tbl$ref %||% NA_character_),
    alt = as.character(ld_order_tbl$alt %||% NA_character_),
    z_score = as.numeric(z_aligned$z_score),
    slope = as.numeric(z_aligned$slope %||% NA_real_),
    slope_se = as.numeric(z_aligned$slope_se %||% NA_real_),
    sample_size = as.numeric(z_aligned$sample_size %||% NA_real_),
    af = as.numeric(z_aligned$af %||% NA_real_),
    annotation_a = as.numeric(annotation_tbl$annotation_a),
    annotation_missing = as.logical(annotation_tbl$annotation_missing %||% FALSE),
    beta_hat_std = as.numeric(annotation_tbl$beta_hat_std),
    mu0_source = as.numeric(annotation_tbl$mu0 %||% NA_real_),
    var_y_hat_from_slope = as.numeric(annotation_tbl$var_y_hat_from_slope),
    var_y_hat_from_se = as.numeric(annotation_tbl$var_y_hat_from_se),
    baseline_c_l = as.numeric(annotation_tbl$baseline_c_l)
  )

  # Sample-size filter: drop variants whose effective n is < 0.5 * max(n).
  # This bounds within-locus n spread to 2x so the single-n assumption used by
  # susie_rss / susine_rss is more honest, and pairs with `n = min(...)` below.
  # Skipped when apply_sample_size_filter = FALSE (e.g. for post-hoc analysis
  # of fits trained before the filter was introduced).
  ss_vec <- variant_map$sample_size
  ss_max <- suppressWarnings(max(ss_vec, na.rm = TRUE))
  if (apply_sample_size_filter && is.finite(ss_max) && ss_max > 0) {
    ss_threshold <- 0.5 * ss_max
    keep <- !is.na(ss_vec) & ss_vec >= ss_threshold
  } else {
    ss_threshold <- NA_real_
    keep <- rep(TRUE, nrow(variant_map))
  }
  n_dropped_sample_size <- sum(!keep)
  if (n_dropped_sample_size > 0L) {
    if (sum(keep) < 2L) {
      stop("After sample-size filter (n >= 0.5 * max), fewer than 2 variants remain for ", locus_id)
    }
    variant_map <- variant_map[keep, , drop = FALSE]
    R <- R[keep, keep, drop = FALSE]
    # Re-index 1..p_new so downstream sort/join on ld_matrix_index stays valid
    variant_map$ld_matrix_index <- seq_len(nrow(variant_map))
  }

  baseline_c <- unique(stats::na.omit(variant_map$baseline_c_l))
  if (length(baseline_c) > 1L) {
    stop("Multiple baseline_c_l values found for ", locus_id)
  }

  list(
    dataset_bundle_id = as.integer(row$dataset_bundle_id[[1]]),
    locus_id = locus_id,
    gene_name = row$gene_name[[1]],
    gtex_tissue = row$gtex_tissue[[1]],
    gtex_chrom = row$gtex_chrom[[1]],
    z = as.numeric(variant_map$z_score),
    R = R,
    a = as.numeric(variant_map$annotation_a),
    variant_map = variant_map,
    n_sample = as.numeric(stats::median(variant_map$sample_size, na.rm = TRUE)),
    n_sample_min = as.numeric(min(variant_map$sample_size, na.rm = TRUE)),
    n_sample_max = as.numeric(max(variant_map$sample_size, na.rm = TRUE)),
    n_sample_median = as.numeric(stats::median(variant_map$sample_size, na.rm = TRUE)),
    n_sample_threshold = as.numeric(ss_threshold),
    n_dropped_sample_size = as.integer(n_dropped_sample_size),
    baseline_c_l = as.numeric(real_data_scalar_or_na(baseline_c, NA_real_)),
    beta_hat_std = as.numeric(variant_map$beta_hat_std),
    var_y_hat_from_slope = as.numeric(variant_map$var_y_hat_from_slope),
    var_y_hat_from_se = as.numeric(variant_map$var_y_hat_from_se),
    paths = list(
      annotation = annotation_path,
      z_scores = z_scores_path,
      ld_long = ld_long_path,
      ld_order = ld_order_path,
      master_variants = master_variants_path
    )
  )
}

build_real_data_job_config <- function(
    job_name = "real_data_ensemble_geometric_n20",
    manifest_path = real_data_default_manifest_path(),
    loci = NULL,
    output_root = "output",
    L = 10L,
    c_grid_values = seq(0, 1.5, length.out = 8),
    sigma_grid_values = c(0.01, 0.03, 0.07, 0.1, 0.2, 0.4, 0.7, 1.0),
    max_iter = 100L,
    tol = 1e-5,
    cs_coverage = 0.95,
    cs_min_purity = 0.5,
    jsd_threshold = 0.15,
    softmax_temperature = 1,
    estimate_residual_variance = FALSE,
    email = "mgc5166@psu.edu",
    time = "05:59:59",
    mem = "8G",
    cpus_per_task = 1L,
    partition = NULL,
    account = NULL,
    HPC = TRUE
) {
  repo_root <- ensure_repo_root(getwd())
  manifest_path <- normalizePath(manifest_path, winslash = "/", mustWork = TRUE)
  bundles <- readr::read_csv(manifest_path, show_col_types = FALSE)
  if (!is.null(loci)) {
    bundles <- dplyr::filter(bundles, .data$locus_id %in% loci)
  }
  if (!nrow(bundles)) {
    stop("No real-data loci selected for the job config.")
  }

  bundles <- bundles %>%
    dplyr::arrange(.data$dataset_bundle_id) %>%
    dplyr::mutate(task_id = dplyr::row_number())

  functional_grid <- tidyr::crossing(
    tibble::tibble(
      sigma_index = seq_along(sigma_grid_values),
      sigma_0_2_scalar = as.numeric(sigma_grid_values)
    ),
    tibble::tibble(
      c_index = seq_along(c_grid_values),
      c_value = as.numeric(c_grid_values)
    )
  ) %>%
    dplyr::mutate(
      backend = "susine_rss",
      run_family = "functional_grid",
      flush_group = as.integer(.data$sigma_index)
    )
  anchor_flush_group <- length(sigma_grid_values) + 1L
  runs_per_task <- nrow(functional_grid) + 1L

  runs <- purrr::map_dfr(seq_len(nrow(bundles)), function(i) {
    bundle <- bundles[i, , drop = FALSE]
    functional <- functional_grid %>%
      dplyr::transmute(
        dataset_bundle_id = bundle$dataset_bundle_id[[1]],
        task_id = bundle$task_id[[1]],
        locus_id = bundle$locus_id[[1]],
        gene_name = bundle$gene_name[[1]],
        baseline_c_l = bundle$baseline_c_l[[1]],
        backend = .data$backend,
        run_family = .data$run_family,
        c_value = .data$c_value,
        sigma_0_2_scalar = .data$sigma_0_2_scalar,
        group_key = paste0("functional_grid|", bundle$locus_id[[1]]),
        flush_group = .data$flush_group
      )

    anchor <- tibble::tibble(
      dataset_bundle_id = bundle$dataset_bundle_id[[1]],
      task_id = bundle$task_id[[1]],
      locus_id = bundle$locus_id[[1]],
      gene_name = bundle$gene_name[[1]],
      baseline_c_l = bundle$baseline_c_l[[1]],
      backend = "susie_rss",
      run_family = "susie_anchor",
      c_value = NA_real_,
      sigma_0_2_scalar = 0.2,
      group_key = paste0("susie_anchor|", bundle$locus_id[[1]]),
      flush_group = anchor_flush_group
    )

    dplyr::bind_rows(functional, anchor)
  }) %>%
    dplyr::mutate(run_id = dplyr::row_number(), .before = 1L)

  tasks <- bundles %>%
    dplyr::transmute(
      task_id = .data$task_id,
      dataset_bundle_id = .data$dataset_bundle_id,
      locus_id = .data$locus_id,
      gene_name = .data$gene_name,
      runs_per_task = as.integer(runs_per_task)
    )

  output_root <- normalizePath(output_root, winslash = "/", mustWork = FALSE)
  list(
    job = list(
      name = job_name,
      email = email,
      created_at = timestamp_utc(),
      manifest_path = manifest_path,
      study_root = dirname(manifest_path),
      L = as.integer(L),
      max_iter = as.integer(max_iter),
      tol = as.numeric(tol),
      cs_coverage = as.numeric(cs_coverage),
      cs_min_purity = as.numeric(cs_min_purity),
      jsd_threshold = as.numeric(jsd_threshold),
      softmax_temperature = as.numeric(softmax_temperature),
      estimate_residual_variance = isTRUE(estimate_residual_variance),
      c_grid_values = as.numeric(c_grid_values),
      sigma_grid_values = as.numeric(sigma_grid_values),
      HPC = isTRUE(HPC),
      slurm = list(
        time = time,
        mem = mem,
        cpus_per_task = as.integer(cpus_per_task),
        partition = partition,
        account = account
      )
    ),
    paths = list(
      repo_root = repo_root,
      output_root = output_root,
      temp_dir = file.path(output_root, "temp", job_name),
      slurm_scripts_dir = file.path(output_root, "slurm_scripts"),
      slurm_prints_dir = file.path(output_root, "slurm_prints"),
      slurm_output_dir = file.path(output_root, "slurm_output")
    ),
    tables = list(
      dataset_bundles = bundles,
      runs = runs,
      tasks = tasks
    )
  )
}

write_real_data_job_artifacts <- function(job_config, run_task_script) {
  paths <- job_config$paths
  ensure_dir(paths$temp_dir)
  ensure_dir(paths$slurm_scripts_dir)
  ensure_dir(paths$slurm_prints_dir)
  ensure_dir(paths$slurm_output_dir)

  unlink(paths$temp_dir, recursive = TRUE)
  ensure_dir(paths$temp_dir)

  run_task_script <- normalizePath(run_task_script, winslash = "/", mustWork = TRUE)

  job_config_json <- job_config
  job_config_json$tables <- NULL
  job_config_path <- file.path(paths$temp_dir, "job_config.json")
  jsonlite::write_json(
    job_config_json,
    path = job_config_path,
    auto_unbox = TRUE,
    digits = NA,
    pretty = TRUE
  )

  run_manifest_path <- file.path(paths$temp_dir, "run_manifest.csv")
  dataset_bundles_path <- file.path(paths$temp_dir, "dataset_bundles.csv")
  task_table_path <- file.path(paths$temp_dir, "task_table.csv")
  readr::write_csv(job_config$tables$runs, run_manifest_path)
  readr::write_csv(job_config$tables$dataset_bundles, dataset_bundles_path)
  readr::write_csv(job_config$tables$tasks, task_table_path)

  slurm_path <- file.path(paths$slurm_scripts_dir, paste0(job_config$job$name, ".slurm"))
  writeLines(
    render_real_data_slurm_script(job_config, run_task_script = run_task_script),
    con = slurm_path
  )

  list(
    job_config = job_config_path,
    run_manifest = run_manifest_path,
    dataset_bundles = dataset_bundles_path,
    task_table = task_table_path,
    slurm_script = slurm_path
  )
}

render_real_data_slurm_script <- function(job_config, run_task_script) {
  job <- job_config$job
  paths <- job_config$paths
  tasks <- job_config$tables$tasks
  n_tasks <- nrow(tasks)
  slurm <- job$slurm

  partition_line <- if (!is.null(slurm$partition)) {
    paste0("#SBATCH --partition=", slurm$partition)
  } else {
    NULL
  }

  account_line <- if (!is.null(slurm$account)) {
    paste0("#SBATCH --account=", slurm$account)
  } else {
    NULL
  }

  hpc_setup <- if (isTRUE(job$HPC)) {
    c(
      "module load r",
      "",
      'export R_LIBS_USER="/storage/home/mgc5166/R/x86_64-pc-linux-gnu-library/4.3"',
      ""
    )
  } else {
    NULL
  }

  script <- c(
    "#!/bin/bash",
    sprintf("#SBATCH --job-name=%s", job$name),
    sprintf("#SBATCH --array=1-%d", n_tasks),
    sprintf("#SBATCH --time=%s", slurm$time),
    sprintf("#SBATCH --mem=%s", slurm$mem),
    sprintf("#SBATCH --cpus-per-task=%s", slurm$cpus_per_task),
    sprintf("#SBATCH --mail-user=%s", job$email),
    "#SBATCH --mail-type=BEGIN,END,FAIL",
    partition_line,
    account_line,
    "#SBATCH --output=/dev/null",
    "#SBATCH --error=/dev/null",
    "",
    "set -euo pipefail",
    "",
    sprintf('JOB_ROOT="%s"', normalizePath(paths$output_root, winslash = "/", mustWork = FALSE)),
    sprintf('CONFIG_PATH="%s"', normalizePath(file.path(paths$temp_dir, "job_config.json"), winslash = "/", mustWork = FALSE)),
    sprintf('RUN_TASK_SCRIPT="%s"', normalizePath(run_task_script, winslash = "/", mustWork = FALSE)),
    sprintf('SLURM_PRINTS_BASE="%s"', normalizePath(paths$slurm_prints_dir, winslash = "/", mustWork = FALSE)),
    sprintf('SLURM_OUTPUT_BASE="%s"', normalizePath(paths$slurm_output_dir, winslash = "/", mustWork = FALSE)),
    sprintf('TEMP_DIR="%s"', normalizePath(paths$temp_dir, winslash = "/", mustWork = FALSE)),
    sprintf('RUN_HISTORY_BASE="%s"', normalizePath(file.path(paths$output_root, "run_history"), winslash = "/", mustWork = FALSE)),
    "",
    'JOBNAME="${SLURM_JOB_NAME}"',
    'PARENT_ID="${SLURM_ARRAY_JOB_ID:-$SLURM_JOB_ID}"',
    'TASK_ID="${SLURM_ARRAY_TASK_ID:-0}"',
    'PRINTS_DIR="${SLURM_PRINTS_BASE}/${JOBNAME}/${PARENT_ID}"',
    'mkdir -p "${PRINTS_DIR}"',
    "",
    'exec >"${PRINTS_DIR}/${TASK_ID}.out" 2>"${PRINTS_DIR}/${TASK_ID}.err"',
    'echo "[$(date -Is)] Starting task ${TASK_ID} for job ${JOBNAME} (parent ${PARENT_ID})"',
    "",
    'export SUSINE_JOB_NAME="${JOBNAME}"',
    'export SUSINE_PARENT_ID="${PARENT_ID}"',
    'export SUSINE_DEV="1"',
    hpc_setup,
    'if [ "${TASK_ID}" = "1" ]; then',
    '  FINAL_HISTORY_DIR="${RUN_HISTORY_BASE}/${JOBNAME}/${PARENT_ID}"',
    '  mkdir -p "${FINAL_HISTORY_DIR}"',
    '  cp "${TEMP_DIR}"/* "${FINAL_HISTORY_DIR}/"',
    '  echo "[$(date -Is)] Task 1: copied run_history from temp to ${FINAL_HISTORY_DIR}"',
    'fi',
    "",
    'Rscript "$RUN_TASK_SCRIPT" \\',
    '  --job-name "$JOBNAME" \\',
    '  --task-id "$TASK_ID" \\',
    '  --job-root "$JOB_ROOT" \\',
    '  --config-path "$CONFIG_PATH"',
    "",
    'echo "[$(date -Is)] Completed task ${TASK_ID}"'
  )
  script[!is.na(script)]
}

real_data_read_job_artifacts <- function(config_path) {
  cfg <- jsonlite::fromJSON(config_path, simplifyVector = TRUE)
  temp_dir <- dirname(normalizePath(config_path, winslash = "/", mustWork = TRUE))
  list(
    job_config = cfg,
    run_manifest = readr::read_csv(file.path(temp_dir, "run_manifest.csv"), show_col_types = FALSE),
    dataset_bundles = readr::read_csv(file.path(temp_dir, "dataset_bundles.csv"), show_col_types = FALSE),
    task_table = readr::read_csv(file.path(temp_dir, "task_table.csv"), show_col_types = FALSE)
  )
}

real_data_safe_last <- function(x) {
  x <- as.numeric(x)
  x <- x[is.finite(x)]
  if (!length(x)) NA_real_ else tail(x, 1L)
}

real_data_safe_first <- function(x) {
  x <- as.numeric(x)
  x <- x[is.finite(x)]
  if (!length(x)) NA_real_ else x[[1]]
}

real_data_top_overlap <- function(x, y, k) {
  if (!length(x) || !length(y)) {
    return(NA_real_)
  }
  idx_x <- head(order(x, decreasing = TRUE), min(k, length(x)))
  idx_y <- head(order(y, decreasing = TRUE), min(k, length(y)))
  union_len <- length(union(idx_x, idx_y))
  if (union_len == 0L) {
    return(NA_real_)
  }
  length(intersect(idx_x, idx_y)) / union_len
}

real_data_pve_from_posterior_mean <- function(posterior_mean, R) {
  b <- as.numeric(posterior_mean)
  if (!length(b) || is.null(R) || nrow(R) != length(b) || ncol(R) != length(b)) {
    return(NA_real_)
  }
  val <- as.numeric(crossprod(b, R %*% b))
  if (!is.finite(val)) {
    return(NA_real_)
  }
  pmin(pmax(val, 0), 1)
}

real_data_pip_count <- function(pip, threshold = 0.25) {
  sum(as.numeric(pip) > threshold, na.rm = TRUE)
}

real_data_pip_intersection_count <- function(pip_a, pip_b, threshold = 0.25) {
  a <- as.numeric(pip_a)
  b <- as.numeric(pip_b)
  n <- min(length(a), length(b))
  if (!n) {
    return(NA_integer_)
  }
  sum(a[seq_len(n)] > threshold & b[seq_len(n)] > threshold, na.rm = TRUE)
}

real_data_safe_js_distance <- function(pip_a, pip_b) {
  if (is.null(pip_a) || is.null(pip_b)) {
    return(NA_real_)
  }
  a <- as.numeric(pip_a)
  b <- as.numeric(pip_b)
  if (length(a) != length(b) || !length(a)) {
    return(NA_real_)
  }
  as.numeric(js_distance(a, b))
}

real_data_init_effect_fits_from_susine <- function(fit) {
  ef <- fit$effect_fits %||% NULL
  if (is.null(ef) || is.null(ef$alpha) || is.null(ef$b_hat) || is.null(ef$b_2_hat)) {
    stop("Source SuSiNE fit does not contain alpha, b_hat, and b_2_hat effect_fits.")
  }
  list(
    alpha = as.matrix(ef$alpha),
    b_hat = as.matrix(ef$b_hat),
    b_2_hat = as.matrix(ef$b_2_hat)
  )
}

real_data_run_highest_weight_refit <- function(source_fit, bundle, job_config,
                                               sigma_0_2_scalar = 0.2) {
  if (!requireNamespace("susine", quietly = TRUE)) {
    stop("Package 'susine' is required for highest-weight warm refits.")
  }
  p <- length(bundle$z)
  susine::susine_rss(
    L = as.integer(job_config$job$L),
    z = bundle$z,
    R = bundle$R,
    n = bundle$n_sample,
    mu_0 = 0,
    sigma_0_2 = as.numeric(sigma_0_2_scalar),
    prior_inclusion_weights = rep(1 / p, p),
    prior_update_method = "none",
    estimate_residual_variance = isTRUE(job_config$job$estimate_residual_variance),
    verbose = FALSE,
    convergence_method = "elbo",
    tol = as.numeric(job_config$job$tol %||% 1e-5),
    max_iter = as.integer(job_config$job$max_iter %||% 100L),
    init_effect_fits = real_data_init_effect_fits_from_susine(source_fit)
  )
}

real_data_write_parquet <- function(df, path) {
  ensure_dir(dirname(path))
  arrow::write_parquet(arrow::Table$create(df), path, compression = "zstd")
}

real_data_backend_matrices <- function(fit, backend) {
  if (identical(backend, "susine_rss")) {
    alpha <- as.matrix(fit$effect_fits$alpha)
    alpha_b_hat <- as.matrix(fit$effect_fits$b_hat)
    alpha_b_2_hat <- as.matrix(fit$effect_fits$b_2_hat)
    b_hat <- ifelse(alpha > 0, alpha_b_hat / alpha, 0)
    b_2_hat <- ifelse(alpha > 0, alpha_b_2_hat / alpha, 0)
    pips <- as.numeric(fit$model_fit$PIPs)
    elbo <- as.numeric(fit$model_fit$elbo)
    sigma_2_trace <- as.numeric(fit$model_fit$sigma_2)
    n_iter <- max(length(elbo) - 1L, 0L)
    converged <- if (!is.null(fit$model_fit$alpha_diff) && length(fit$model_fit$alpha_diff)) {
      as.logical(tail(fit$model_fit$alpha_diff, 1L) <= (fit$settings$tol %||% 1e-5))
    } else if (!is.null(fit$settings$max_iter)) {
      n_iter < fit$settings$max_iter
    } else {
      NA
    }
    return(list(
      alpha = alpha,
      b_hat = b_hat,
      b_2_hat = b_2_hat,
      alpha_b_hat = alpha_b_hat,
      alpha_b_2_hat = alpha_b_2_hat,
      pip = pips,
      elbo = elbo,
      sigma_2_trace = sigma_2_trace,
      sigma_2_final = real_data_safe_last(sigma_2_trace),
      n_iter = n_iter,
      converged = converged
    ))
  }

  if (identical(backend, "susie_rss")) {
    alpha <- as.matrix(fit$alpha)
    b_hat <- as.matrix(fit$mu)
    b_2_hat <- as.matrix(fit$mu2)
    alpha_b_hat <- alpha * b_hat
    alpha_b_2_hat <- alpha * b_2_hat
    pips <- fit$pip %||% (1 - apply(1 - alpha, 2L, prod))
    elbo <- as.numeric(fit$elbo %||% numeric())
    sigma_2_trace <- as.numeric(fit$sigma2 %||% fit$sigma_2 %||% numeric())
    n_iter <- as.integer(fit$niter %||% max(length(elbo) - 1L, 0L))
    return(list(
      alpha = alpha,
      b_hat = b_hat,
      b_2_hat = b_2_hat,
      alpha_b_hat = alpha_b_hat,
      alpha_b_2_hat = alpha_b_2_hat,
      pip = as.numeric(pips),
      elbo = elbo,
      sigma_2_trace = sigma_2_trace,
      sigma_2_final = real_data_safe_last(sigma_2_trace),
      n_iter = n_iter,
      converged = as.logical(fit$converged %||% NA)
    ))
  }

  stop("Unsupported backend: ", backend)
}

compute_real_data_cs_summaries <- function(
    alpha,
    variant_map,
    R,
    run_id,
    locus_id,
    gene_name,
    backend,
    cs_coverage = 0.95,
    cs_min_purity = 0.5
) {
  L <- nrow(alpha)
  raw_sets <- list()
  filtered_sets <- list()
  membership_rows <- list()
  effect_rows <- list()

  for (l in seq_len(L)) {
    alpha_l <- as.numeric(alpha[l, ])
    alpha_prob <- normalize_prob_vec(alpha_l)
    cs_raw <- get_credible_set(alpha_prob, rho = cs_coverage)
    purity_raw <- cs_purity_min_abs(NULL, cs_raw, cor_mat = R)
    keep_filtered <- length(cs_raw) > 0L && is.finite(purity_raw) && purity_raw >= cs_min_purity
    cs_filtered <- if (keep_filtered) cs_raw else integer(0)
    purity_filtered <- if (keep_filtered) purity_raw else NA_real_

    raw_sets[[length(raw_sets) + 1L]] <- cs_raw
    filtered_sets[[length(filtered_sets) + 1L]] <- cs_filtered

    add_membership <- function(set_idx, set_type) {
      if (!length(set_idx)) {
        return(invisible(NULL))
      }
      ord <- set_idx[order(alpha_prob[set_idx], decreasing = TRUE)]
      cum_alpha <- cumsum(alpha_prob[ord])
      membership_rows[[length(membership_rows) + 1L]] <<- tibble::tibble(
        run_id = as.integer(run_id),
        locus_id = as.character(locus_id),
        gene_name = as.character(gene_name),
        effect_l = as.integer(l),
        backend = as.character(backend),
        set_type = as.character(set_type),
        member_rank = seq_along(ord),
        variant_id = as.character(variant_map$variant_id[ord]),
        ld_matrix_index = as.integer(variant_map$ld_matrix_index[ord]),
        alpha_value = as.numeric(alpha_prob[ord]),
        cumulative_alpha = as.numeric(cum_alpha)
      )
    }

    add_membership(cs_raw, "raw")
    add_membership(cs_filtered, "filtered")

    top_idx <- which.max(alpha_prob)
    effect_rows[[length(effect_rows) + 1L]] <- tibble::tibble(
      run_id = as.integer(run_id),
      locus_id = as.character(locus_id),
      gene_name = as.character(gene_name),
      effect_l = as.integer(l),
      backend = as.character(backend),
      alpha_mass = sum(alpha_l, na.rm = TRUE),
      alpha_max = max(alpha_prob, na.rm = TRUE),
      alpha_entropy = prob_entropy(alpha_prob),
      alpha_k_eff = prob_k_eff(alpha_prob),
      top_variant_id = as.character(variant_map$variant_id[[top_idx]]),
      top_ld_index = as.integer(variant_map$ld_matrix_index[[top_idx]]),
      cs_size_raw = length(cs_raw),
      cs_size_filtered = length(cs_filtered),
      cs_purity_raw = as.numeric(purity_raw),
      cs_purity_filtered = as.numeric(purity_filtered)
    )
  }

  effect_tbl <- dplyr::bind_rows(effect_rows)
  membership_tbl <- dplyr::bind_rows(membership_rows)

  valid_raw_sets <- raw_sets[vapply(raw_sets, length, integer(1)) > 0L]
  valid_filtered_sets <- filtered_sets[vapply(filtered_sets, length, integer(1)) > 0L]

  list(
    effect_summaries = effect_tbl,
    credible_set_membership = membership_tbl,
    summary = list(
      cs_count_raw = length(valid_raw_sets),
      cs_count_filtered = length(valid_filtered_sets),
      mean_cs_size_raw = if (nrow(effect_tbl)) mean(effect_tbl$cs_size_raw, na.rm = TRUE) else NA_real_,
      mean_cs_size_filtered = if (nrow(effect_tbl)) mean(effect_tbl$cs_size_filtered, na.rm = TRUE) else NA_real_,
      mean_cs_purity_raw = if (nrow(effect_tbl)) mean(effect_tbl$cs_purity_raw, na.rm = TRUE) else NA_real_,
      mean_cs_purity_filtered = if (nrow(effect_tbl)) mean(effect_tbl$cs_purity_filtered, na.rm = TRUE) else NA_real_,
      cs_overlap_rate_raw = overlap_rate_from_sets(valid_raw_sets),
      cs_overlap_rate_filtered = overlap_rate_from_sets(valid_filtered_sets)
    )
  )
}

compute_real_data_run_metrics <- function(
    run_row,
    bundle,
    fit,
    backend_mats,
    cs_info,
    mu_0,
    fit_rds_path,
    wall_time_sec
) {
  pips <- as.numeric(backend_mats$pip)
  elbo <- as.numeric(backend_mats$elbo)
  sigma_2_final <- as.numeric(backend_mats$sigma_2_final)
  posterior_mean <- colSums(backend_mats$alpha_b_hat, na.rm = TRUE)
  pve_postmean_std <- real_data_pve_from_posterior_mean(posterior_mean, bundle$R)
  top_ord <- order(pips, decreasing = TRUE)
  top_mass <- function(k) sum(pips[head(top_ord, min(k, length(top_ord)))], na.rm = TRUE)

  tibble::tibble(
    run_id = as.integer(run_row$run_id),
    dataset_bundle_id = as.integer(run_row$dataset_bundle_id),
    locus_id = as.character(run_row$locus_id),
    gene_name = as.character(run_row$gene_name),
    backend = as.character(run_row$backend),
    run_family = as.character(run_row$run_family),
    c_value = as.numeric(run_row$c_value),
    sigma_0_2_scalar = as.numeric(run_row$sigma_0_2_scalar),
    baseline_c_l = as.numeric(run_row$baseline_c_l),
    annotation_used = identical(run_row$backend[[1]], "susine_rss") &&
      is.finite(run_row$c_value[[1]]) &&
      abs(run_row$c_value[[1]]) > 0,
    elbo_first = real_data_safe_first(elbo),
    elbo_final = real_data_safe_last(elbo),
    elbo_gain = real_data_safe_last(elbo) - real_data_safe_first(elbo),
    n_iter = as.integer(backend_mats$n_iter),
    converged = as.logical(backend_mats$converged),
    sigma_2_final = sigma_2_final,
    h2_proxy_std = pmin(pmax(1 - sigma_2_final, 0), 1),
    pve_postmean_std = as.numeric(pve_postmean_std),
    sum_pip = sum(pips, na.rm = TRUE),
    max_pip = max(pips, na.rm = TRUE),
    pip_entropy = prob_entropy(pips),
    n_snps_pip_ge_0_1 = sum(pips >= 0.1, na.rm = TRUE),
    n_snps_pip_ge_0_5 = sum(pips >= 0.5, na.rm = TRUE),
    n_snps_pip_ge_0_9 = sum(pips >= 0.9, na.rm = TRUE),
    pip_mass_top1 = top_mass(1L),
    pip_mass_top5 = top_mass(5L),
    pip_mass_top10 = top_mass(10L),
    L_nominal = as.integer(fit$settings$L %||% nrow(backend_mats$alpha)),
    L_active = sum(apply(backend_mats$alpha, 1L, function(row) {
      max(row, na.rm = TRUE) > (1 / max(length(row), 1L) + 1e-8)
    })),
    cs_count_raw = as.integer(cs_info$summary$cs_count_raw),
    cs_count_filtered = as.integer(cs_info$summary$cs_count_filtered),
    mean_cs_size_raw = as.numeric(cs_info$summary$mean_cs_size_raw),
    mean_cs_size_filtered = as.numeric(cs_info$summary$mean_cs_size_filtered),
    mean_cs_purity_raw = as.numeric(cs_info$summary$mean_cs_purity_raw),
    mean_cs_purity_filtered = as.numeric(cs_info$summary$mean_cs_purity_filtered),
    cs_overlap_rate_raw = as.numeric(cs_info$summary$cs_overlap_rate_raw),
    cs_overlap_rate_filtered = as.numeric(cs_info$summary$cs_overlap_rate_filtered),
    mu0_mean = mean(mu_0, na.rm = TRUE),
    mu0_sd = stats::sd(mu_0, na.rm = TRUE),
    mu0_l2 = sqrt(sum(mu_0^2, na.rm = TRUE)),
    mu0_corr_z = suppressWarnings(stats::cor(mu_0, bundle$z, use = "pairwise.complete.obs")),
    mu0_corr_beta_hat_std = suppressWarnings(stats::cor(mu_0, bundle$beta_hat_std, use = "pairwise.complete.obs")),
    wall_time_sec = as.numeric(wall_time_sec),
    fit_rds_path = as.character(fit_rds_path),
    error_message = NA_character_
  )
}

# Stub run_metrics row written when a fit errors, so failed runs remain visible
# to downstream aggregators. Schema must match compute_real_data_run_metrics().
real_data_failed_run_metrics_stub <- function(run_row, error_message) {
  tibble::tibble(
    run_id = as.integer(run_row$run_id[[1]]),
    dataset_bundle_id = as.integer(run_row$dataset_bundle_id[[1]]),
    locus_id = as.character(run_row$locus_id[[1]]),
    gene_name = as.character(run_row$gene_name[[1]]),
    backend = as.character(run_row$backend[[1]]),
    run_family = as.character(run_row$run_family[[1]]),
    c_value = as.numeric(run_row$c_value[[1]]),
    sigma_0_2_scalar = as.numeric(run_row$sigma_0_2_scalar[[1]]),
    baseline_c_l = as.numeric(run_row$baseline_c_l[[1]]),
    annotation_used = NA,
    elbo_first = NA_real_,
    elbo_final = NA_real_,
    elbo_gain = NA_real_,
    n_iter = NA_integer_,
    converged = FALSE,
    sigma_2_final = NA_real_,
    h2_proxy_std = NA_real_,
    pve_postmean_std = NA_real_,
    sum_pip = NA_real_,
    max_pip = NA_real_,
    pip_entropy = NA_real_,
    n_snps_pip_ge_0_1 = NA_integer_,
    n_snps_pip_ge_0_5 = NA_integer_,
    n_snps_pip_ge_0_9 = NA_integer_,
    pip_mass_top1 = NA_real_,
    pip_mass_top5 = NA_real_,
    pip_mass_top10 = NA_real_,
    L_nominal = NA_integer_,
    L_active = NA_integer_,
    cs_count_raw = NA_integer_,
    cs_count_filtered = NA_integer_,
    mean_cs_size_raw = NA_real_,
    mean_cs_size_filtered = NA_real_,
    mean_cs_purity_raw = NA_real_,
    mean_cs_purity_filtered = NA_real_,
    cs_overlap_rate_raw = NA_real_,
    cs_overlap_rate_filtered = NA_real_,
    mu0_mean = NA_real_,
    mu0_sd = NA_real_,
    mu0_l2 = NA_real_,
    mu0_corr_z = NA_real_,
    mu0_corr_beta_hat_std = NA_real_,
    wall_time_sec = NA_real_,
    fit_rds_path = NA_character_,
    error_message = as.character(error_message)
  )
}

compute_real_data_dataset_metrics <- function(bundle) {
  hl_count <- high_ld_count(bundle$R, threshold = 0.95)
  z_metrics <- z_score_metrics(bundle$z, top_k = 10L)
  var_y_median_slope <- stats::median(bundle$var_y_hat_from_slope, na.rm = TRUE)
  var_y_median_se <- stats::median(bundle$var_y_hat_from_se, na.rm = TRUE)
  var_y_iqr <- stats::IQR(bundle$var_y_hat_from_slope, na.rm = TRUE)
  ld_z_diag <- real_data_ld_z_consistency(bundle$R, bundle$z, bundle$n_sample)

  tibble::tibble(
    dataset_bundle_id = as.integer(bundle$dataset_bundle_id),
    locus_id = as.character(bundle$locus_id),
    gene_name = as.character(bundle$gene_name),
    n_variants = length(bundle$z),
    n_sample_used = as.numeric(bundle$n_sample),
    n_sample_min = as.numeric(bundle$n_sample_min %||% NA_real_),
    n_sample_max = as.numeric(bundle$n_sample_max %||% NA_real_),
    n_sample_median = as.numeric(bundle$n_sample_median %||% NA_real_),
    n_sample_threshold = as.numeric(bundle$n_sample_threshold %||% NA_real_),
    n_dropped_sample_size = as.integer(bundle$n_dropped_sample_size %||% NA_integer_),
    baseline_c_l = as.numeric(bundle$baseline_c_l),
    M1 = mid_energy_M1(bundle$R),
    high_ld_count_095 = hl_count,
    high_ld_count_095_per_snp = if (length(bundle$z) > 0) (2 * hl_count) / length(bundle$z) else NA_real_,
    high_ld_frac_095 = if (length(bundle$z) > 1) hl_count / (length(bundle$z) * (length(bundle$z) - 1) / 2) else NA_real_,
    z_topk_ratio = z_metrics$z_topk_ratio,
    z_max_abs = z_metrics$z_max_abs,
    z_count_abs_gt_3 = z_metrics$z_count_abs_gt_3,
    z_eff_signals = z_metrics$z_eff_signals,
    annotation_corr_beta_hat_std = suppressWarnings(stats::cor(bundle$a, bundle$beta_hat_std, use = "pairwise.complete.obs")),
    annotation_corr_z = suppressWarnings(stats::cor(bundle$a, bundle$z, use = "pairwise.complete.obs")),
    annotation_mean_abs = mean(abs(bundle$a), na.rm = TRUE),
    annotation_max_abs = max(abs(bundle$a), na.rm = TRUE),
    var_y_hat_slope_median = as.numeric(var_y_median_slope),
    var_y_hat_se_median = as.numeric(var_y_median_se),
    var_y_hat_iqr = as.numeric(var_y_iqr),
    var_y_near1_flag = isTRUE(is.finite(var_y_median_slope) &&
      is.finite(var_y_median_se) &&
      abs(var_y_median_slope - 1) <= 0.1 &&
      abs(var_y_median_se - 1) <= 0.1),
    mom_sigma2_hat = ld_z_diag$mom_sigma2_hat,
    ld_chol_ok = ld_z_diag$ld_chol_ok,
    ld_z_inconsistency_flag = ld_z_diag$ld_z_inconsistency_flag
  )
}

# Diagnostic for (R, z, n) consistency under the susie_rss model.
# - ld_chol_ok: TRUE iff R is positive-definite without regularization.
# - mom_sigma2_hat: 1 - z' (R + lambda I)^-1 z / n. The susie_rss MoM estimator
#   for residual variance; values < 0 indicate the (R, z, n) triple is
#   internally inconsistent (e.g. external LD reference panel, allele flips,
#   wrong n) and susieR will throw "estimated value is negative".
# - ld_z_inconsistency_flag: TRUE if either R is not PSD or mom_sigma2_hat < 0.
real_data_ld_z_consistency <- function(R, z, n, lambda = 1e-4) {
  p <- length(z)
  if (!is.finite(n) || n <= 0 || p < 2L) {
    return(list(
      mom_sigma2_hat = NA_real_,
      ld_chol_ok = NA,
      ld_z_inconsistency_flag = NA
    ))
  }
  ld_chol_ok <- tryCatch({
    suppressWarnings(chol(R))
    TRUE
  }, error = function(e) FALSE)
  ch <- tryCatch(suppressWarnings(chol(R + diag(lambda, p))),
                 error = function(e) NULL)
  mom_sigma2 <- if (is.null(ch)) {
    NA_real_
  } else {
    qf <- as.numeric(crossprod(backsolve(ch, z, transpose = TRUE)))
    1 - qf / n
  }
  list(
    mom_sigma2_hat = as.numeric(mom_sigma2),
    ld_chol_ok = as.logical(ld_chol_ok),
    ld_z_inconsistency_flag = isTRUE(!ld_chol_ok) ||
                              isTRUE(is.finite(mom_sigma2) && mom_sigma2 < 0)
  )
}

# Pull an L x p coefficient matrix b = alpha * (posterior mean | inclusion)
# from a fit. susine stores this directly as fit$effect_fits$b_hat (already
# alpha-weighted; see real_data_backend_matrices()). susieR stores
# fit$alpha (inclusion probs, L x p) and fit$mu (conditional posterior mean,
# L x p) separately, so their product gives the same quantity.
real_data_extract_b_matrix <- function(fit) {
  if (!is.null(fit$effect_fits$b_hat)) {
    return(as.matrix(fit$effect_fits$b_hat))
  }
  if (!is.null(fit$alpha) && !is.null(fit$mu)) {
    alpha <- as.matrix(fit$alpha)
    mu    <- as.matrix(fit$mu)
    if (nrow(alpha) != nrow(mu) || ncol(alpha) != ncol(mu)) return(NULL)
    return(alpha * mu)
  }
  NULL
}

# Weighted (1 - r^2) drift between two fits' per-effect coefficient vectors,
# matched 1-to-1 by max-weight assignment. RSS analog of weighted_r2_drift in
# effect_matching_drift_analysis_5arm.Rmd: in the SS regime where var(y) = 1
# and X is column-standardized, var(X b_l) = b_l' R b_l and
# cor(X b_i, X b_j) = (b_i' R b_j) / sqrt((b_i' R b_i)(b_j' R b_j)).
# Weights: pmin(var_y_explained_a[i], var_y_explained_b[j]).
real_data_basin_r2_drift_pair <- function(b_a, b_b, R) {
  if (is.null(b_a) || is.null(b_b)) return(list(drift = NA_real_, n_matched = 0L))
  L_a <- nrow(b_a); L_b <- nrow(b_b)
  if (L_a == 0L || L_b == 0L) return(list(drift = NA_real_, n_matched = 0L))
  if (!requireNamespace("clue", quietly = TRUE))
    stop("Package 'clue' required for solve_LSAP. install.packages('clue').")

  Rb_a <- tcrossprod(b_a, R)              # L_a x p ; rows are R %*% b_a[i,]
  Rb_b <- tcrossprod(b_b, R)              # L_b x p
  num  <- b_a %*% t(Rb_b)                 # L_a x L_b ; (b_a[i,]' R b_b[j,])
  var_a <- pmax(rowSums(b_a * Rb_a), 0)   # var_y_explained per A-effect
  var_b <- pmax(rowSums(b_b * Rb_b), 0)
  denom <- sqrt(outer(var_a, var_b))
  cor_mat <- ifelse(denom > 0, num / denom, 0)
  cor_mat <- pmin(pmax(cor_mat, -1), 1)

  # solve_LSAP requires non-negative entries when maximum = TRUE. Cosines live
  # in [-1, 1], so shift by +1 (adding a constant to every entry preserves the
  # optimal assignment). Pad to square with 0 so unmatched effects pair off
  # outside the L_a x L_b block.
  L <- max(L_a, L_b)
  M_sq <- matrix(0, L, L)
  M_sq[seq_len(L_a), seq_len(L_b)] <- cor_mat + 1
  perm <- as.integer(clue::solve_LSAP(M_sq, maximum = TRUE))

  w_sum <- 0
  num_sum <- 0
  n_matched <- 0L
  for (i in seq_len(L_a)) {
    j <- perm[i]
    if (j > L_b) next
    w <- min(var_a[i], var_b[j])
    if (!is.finite(w) || w <= 0) next
    cij <- cor_mat[i, j]
    num_sum <- num_sum + w * (1 - cij^2)
    w_sum <- w_sum + w
    n_matched <- n_matched + 1L
  }
  list(
    drift = if (w_sum > 0) num_sum / w_sum else NA_real_,
    n_matched = n_matched
  )
}

real_data_basin_r2_cache <- function(fit, R) {
  b <- real_data_extract_b_matrix(fit)
  if (is.null(b)) return(NULL)
  if (ncol(b) != ncol(R)) return(NULL)
  Rb <- tcrossprod(b, R)
  list(
    b = b,
    Rb = Rb,
    var = pmax(rowSums(b * Rb), 0)
  )
}

real_data_basin_r2_drift_pair_cached <- function(cache_a, cache_b) {
  if (is.null(cache_a) || is.null(cache_b)) return(list(drift = NA_real_, n_matched = 0L))
  L_a <- nrow(cache_a$b); L_b <- nrow(cache_b$b)
  if (L_a == 0L || L_b == 0L) return(list(drift = NA_real_, n_matched = 0L))
  if (!requireNamespace("clue", quietly = TRUE))
    stop("Package 'clue' required for solve_LSAP. install.packages('clue').")

  num <- cache_a$b %*% t(cache_b$Rb)
  denom <- sqrt(outer(cache_a$var, cache_b$var))
  cor_mat <- ifelse(denom > 0, num / denom, 0)
  cor_mat <- pmin(pmax(cor_mat, -1), 1)

  L <- max(L_a, L_b)
  M_sq <- matrix(0, L, L)
  M_sq[seq_len(L_a), seq_len(L_b)] <- cor_mat + 1
  perm <- as.integer(clue::solve_LSAP(M_sq, maximum = TRUE))

  w_sum <- 0
  num_sum <- 0
  n_matched <- 0L
  for (i in seq_len(L_a)) {
    j <- perm[i]
    if (j > L_b) next
    w <- min(cache_a$var[i], cache_b$var[j])
    if (!is.finite(w) || w <= 0) next
    cij <- cor_mat[i, j]
    num_sum <- num_sum + w * (1 - cij^2)
    w_sum <- w_sum + w
    n_matched <- n_matched + 1L
  }
  list(
    drift = if (w_sum > 0) num_sum / w_sum else NA_real_,
    n_matched = n_matched
  )
}

# Compute basin r^2 drift for every (sigma_0_2, c) cell in a locus's
# functional_grid runs, vs the same-sigma c=0 baseline. Returns a long tibble.
# `runs_tbl` must have columns: run_id, c_value, sigma_0_2_scalar, fit_rds_path.
real_data_basin_r2_drift_locus <- function(R, runs_tbl, progress = NULL) {
  rows <- list()
  fit_cache <- vector("list", nrow(runs_tbl))
  for (i in seq_len(nrow(runs_tbl))) {
    path <- runs_tbl$fit_rds_path[[i]]
    if (!is.na(path)) {
      fit <- tryCatch(readRDS(path), error = function(e) NULL)
      fit_cache[[i]] <- real_data_basin_r2_cache(fit, R)
    }
    if (is.function(progress)) {
      progress(stage = "cached fit", index = i, total = nrow(runs_tbl))
    }
  }
  sigmas <- sort(unique(runs_tbl$sigma_0_2_scalar))
  for (sig in sigmas) {
    sig_runs <- runs_tbl %>%
      dplyr::mutate(.drift_cache_idx = dplyr::row_number()) %>%
      dplyr::filter(abs(.data$sigma_0_2_scalar - !!sig) < 1e-12) %>%
      dplyr::arrange(.data$c_value)
    base_row <- sig_runs %>% dplyr::filter(abs(.data$c_value) < 1e-12)
    if (nrow(base_row) != 1L || is.na(base_row$fit_rds_path[[1]])) next
    base_cache <- fit_cache[[base_row$.drift_cache_idx[[1]]]]
    if (is.null(base_cache)) next
    for (k in seq_len(nrow(sig_runs))) {
      row <- sig_runs[k, , drop = FALSE]
      if (is.na(row$fit_rds_path[[1]])) {
        rows[[length(rows) + 1L]] <- tibble::tibble(
          sigma_0_2_scalar = sig, c_value = row$c_value[[1]],
          basin_r2_drift = NA_real_, n_matched_effects = NA_integer_,
          run_id = row$run_id[[1]]
        )
        next
      }
      drift_res <- real_data_basin_r2_drift_pair_cached(
        fit_cache[[row$.drift_cache_idx[[1]]]],
        base_cache
      )
      rows[[length(rows) + 1L]] <- tibble::tibble(
        sigma_0_2_scalar = sig,
        c_value = row$c_value[[1]],
        basin_r2_drift = drift_res$drift,
        n_matched_effects = drift_res$n_matched,
        run_id = row$run_id[[1]]
      )
      if (is.function(progress)) {
        progress(stage = "computed drift", index = length(rows), total = nrow(runs_tbl))
      }
    }
  }
  dplyr::bind_rows(rows)
}

extract_rss_fit_artifacts <- function(
    fit,
    run_row,
    bundle,
    job_config,
    fit_rds_path,
    wall_time_sec
) {
  backend_mats <- real_data_backend_matrices(fit, backend = run_row$backend)
  mu_0 <- if (identical(run_row$backend, "susine_rss")) {
    as.numeric(run_row$c_value) * bundle$a
  } else {
    rep(0, length(bundle$a))
  }

  cs_info <- compute_real_data_cs_summaries(
    alpha = backend_mats$alpha,
    variant_map = bundle$variant_map,
    R = bundle$R,
    run_id = run_row$run_id,
    locus_id = bundle$locus_id,
    gene_name = bundle$gene_name,
    backend = run_row$backend,
    cs_coverage = job_config$job$cs_coverage,
    cs_min_purity = job_config$job$cs_min_purity
  )

  alpha_max <- apply(backend_mats$alpha, 2L, max, na.rm = TRUE)
  alpha_sum <- colSums(backend_mats$alpha, na.rm = TRUE)
  posterior_mean <- colSums(backend_mats$alpha_b_hat, na.rm = TRUE)
  posterior_second_moment <- colSums(backend_mats$alpha_b_2_hat, na.rm = TRUE)
  posterior_var <- pmax(posterior_second_moment - posterior_mean^2, 0)

  variant_tbl <- bundle$variant_map %>%
    dplyr::transmute(
      run_id = as.integer(run_row$run_id),
      locus_id = as.character(bundle$locus_id),
      gene_name = as.character(bundle$gene_name),
      backend = as.character(run_row$backend),
      run_family = as.character(run_row$run_family),
      variant_id = as.character(.data$variant_id),
      ld_matrix_index = as.integer(.data$ld_matrix_index),
      pip = as.numeric(backend_mats$pip),
      posterior_mean = as.numeric(posterior_mean),
      posterior_second_moment = as.numeric(posterior_second_moment),
      posterior_var = as.numeric(posterior_var),
      alpha_max = as.numeric(alpha_max),
      alpha_sum = as.numeric(alpha_sum),
      z_score = as.numeric(.data$z_score),
      beta_hat_std = as.numeric(.data$beta_hat_std),
      annotation_a = as.numeric(.data$annotation_a),
      mu_0 = as.numeric(mu_0),
      baseline_c_l = as.numeric(.data$baseline_c_l)
    )

  effect_tbl <- purrr::map_dfr(seq_len(nrow(backend_mats$alpha)), function(l) {
    tibble::tibble(
      run_id = as.integer(run_row$run_id),
      locus_id = as.character(bundle$locus_id),
      gene_name = as.character(bundle$gene_name),
      effect_l = as.integer(l),
      backend = as.character(run_row$backend),
      variant_id = as.character(bundle$variant_map$variant_id),
      ld_matrix_index = as.integer(bundle$variant_map$ld_matrix_index),
      alpha = as.numeric(backend_mats$alpha[l, ]),
      b_hat = as.numeric(backend_mats$b_hat[l, ]),
      b_2_hat = as.numeric(backend_mats$b_2_hat[l, ]),
      alpha_b_hat = as.numeric(backend_mats$alpha_b_hat[l, ]),
      alpha_b_2_hat = as.numeric(backend_mats$alpha_b_2_hat[l, ]),
      annotation_a = as.numeric(bundle$variant_map$annotation_a),
      mu_0 = as.numeric(mu_0)
    )
  })

  elbo_trace <- tibble::tibble(
    run_id = as.integer(run_row$run_id),
    locus_id = as.character(bundle$locus_id),
    gene_name = as.character(bundle$gene_name),
    backend = as.character(run_row$backend),
    iteration = seq_along(backend_mats$elbo) - 1L,
    elbo = as.numeric(backend_mats$elbo)
  )

  list(
    run_metrics = compute_real_data_run_metrics(
      run_row = run_row,
      bundle = bundle,
      fit = fit,
      backend_mats = backend_mats,
      cs_info = cs_info,
      mu_0 = mu_0,
      fit_rds_path = fit_rds_path,
      wall_time_sec = wall_time_sec
    ),
    effect_summaries = cs_info$effect_summaries,
    credible_set_membership = cs_info$credible_set_membership,
    variant_posteriors = variant_tbl,
    effect_posteriors = effect_tbl,
    elbo_trace = elbo_trace,
    validation = tibble::tibble(
      run_id = as.integer(run_row$run_id),
      task_id = as.integer(run_row$task_id),
      has_issues = FALSE,
      issues = NA_character_
    )
  )
}

real_data_build_fit_index_row <- function(run_row, fit_rds_path) {
  tibble::tibble(
    run_id = as.integer(run_row$run_id),
    task_id = as.integer(run_row$task_id),
    dataset_bundle_id = as.integer(run_row$dataset_bundle_id),
    locus_id = as.character(run_row$locus_id),
    gene_name = as.character(run_row$gene_name),
    backend = as.character(run_row$backend),
    run_family = as.character(run_row$run_family),
    fit_rds_path = as.character(fit_rds_path)
  )
}

real_data_run_susie_anchor <- function(bundle, run_row, job_config) {
  if (!requireNamespace("susieR", quietly = TRUE)) {
    stop("susieR is required for backend = 'susie_rss'.")
  }
  susie_fn <- getExportedValue("susieR", "susie_rss")
  fn_formals <- names(formals(susie_fn))

  args <- list(
    R = bundle$R,
    n = as.numeric(bundle$n_sample),
    L = as.integer(job_config$job$L)
  )
  if ("z" %in% fn_formals) {
    args$z <- bundle$z
  } else if ("zhat" %in% fn_formals) {
    args$zhat <- bundle$z
  } else {
    stop("Unable to identify z-score argument for susieR::susie_rss().")
  }
  if ("scaled_prior_variance" %in% fn_formals) {
    args$scaled_prior_variance <- as.numeric(run_row$sigma_0_2_scalar)
  } else if ("prior_variance" %in% fn_formals) {
    args$prior_variance <- as.numeric(run_row$sigma_0_2_scalar)
  } else {
    stop("Unable to identify prior variance argument for susieR::susie_rss().")
  }
  if ("estimate_prior_variance" %in% fn_formals) {
    args$estimate_prior_variance <- FALSE
  }
  if ("estimate_residual_variance" %in% fn_formals) {
    args$estimate_residual_variance <- isTRUE(job_config$job$estimate_residual_variance)
  }
  if (!isTRUE(job_config$job$estimate_residual_variance) &&
      "residual_variance" %in% fn_formals) {
    args$residual_variance <- 1
  }
  if ("check_prior" %in% fn_formals) {
    args$check_prior <- FALSE
  }
  if ("max_iter" %in% fn_formals) {
    args$max_iter <- as.integer(job_config$job$max_iter)
  }
  if ("tol" %in% fn_formals) {
    args$tol <- as.numeric(job_config$job$tol)
  }
  if ("convergence_method" %in% fn_formals) {
    args$convergence_method <- "elbo"
  }
  if ("verbose" %in% fn_formals) {
    args$verbose <- FALSE
  }
  if ("coverage" %in% fn_formals) {
    args$coverage <- as.numeric(job_config$job$cs_coverage)
  }
  if ("min_abs_corr" %in% fn_formals) {
    args$min_abs_corr <- as.numeric(job_config$job$cs_min_purity)
  }
  do.call(susie_fn, args)
}

run_real_data_task <- function(job_name, task_id, job_root = "output", config_path, quiet = FALSE) {
  artifacts <- real_data_read_job_artifacts(config_path)
  job_config <- artifacts$job_config
  runs <- artifacts$run_manifest
  bundles <- artifacts$dataset_bundles

  task_id <- as.integer(task_id)
  task_runs <- dplyr::filter(runs, .data$task_id == !!task_id) %>%
    dplyr::arrange(.data$flush_group, .data$run_id)
  if (!nrow(task_runs)) {
    stop("No runs found for task_id ", task_id)
  }

  bundle_id <- unique(task_runs$dataset_bundle_id)
  if (length(bundle_id) != 1L) {
    stop("Each real-data task must map to exactly one dataset bundle.")
  }
  bundle_row <- dplyr::filter(bundles, .data$dataset_bundle_id == !!bundle_id)
  if (nrow(bundle_row) != 1L) {
    stop("Missing dataset bundle row for task_id ", task_id)
  }

  bundle <- load_real_data_locus_bundle(
    locus_id = bundle_row$locus_id[[1]],
    manifest_path = job_config$job$manifest_path,
    repo_root = job_config$paths$repo_root %||% ensure_repo_root(getwd())
  )

  base_output <- determine_base_output(job_config)
  task_dir <- prepare_task_staging_dir(base_output, task_id)
  raw_fits_dir <- file.path(task_dir, "raw_fits")
  ensure_dir(raw_fits_dir)

  dataset_metrics_tbl <- compute_real_data_dataset_metrics(bundle)
  fit_index_rows <- list()
  dataset_metrics_written <- FALSE

  flush_groups <- sort(unique(task_runs$flush_group))
  for (flush_group in flush_groups) {
    flush_label <- sprintf("flush-%03d", as.integer(flush_group))
    group_runs <- dplyr::filter(task_runs, .data$flush_group == !!flush_group) %>%
      dplyr::arrange(.data$run_id)

    run_metrics_rows <- list()
    effect_summary_rows <- list()
    membership_rows <- list()
    variant_rows <- list()
    effect_rows <- list()
    elbo_rows <- list()
    validation_rows <- list()

    for (i in seq_len(nrow(group_runs))) {
      run_row <- group_runs[i, , drop = FALSE]
      start_time <- proc.time()[["elapsed"]]

      fit_result <- tryCatch({
        if (identical(run_row$backend[[1]], "susine_rss")) {
          susine::susine_rss(
            L = as.integer(job_config$job$L),
            z = bundle$z,
            R = bundle$R,
            n = bundle$n_sample,
            mu_0 = as.numeric(run_row$c_value[[1]]) * bundle$a,
            sigma_0_2 = as.numeric(run_row$sigma_0_2_scalar[[1]]),
            prior_update_method = "none",
            estimate_residual_variance = isTRUE(job_config$job$estimate_residual_variance),
            verbose = FALSE,
            convergence_method = "elbo",
            tol = as.numeric(job_config$job$tol),
            max_iter = as.integer(job_config$job$max_iter)
          )
        } else {
          real_data_run_susie_anchor(bundle, run_row, job_config)
        }
      }, error = function(e) {
        e
      })

      if (inherits(fit_result, "error")) {
        err_msg <- conditionMessage(fit_result)
        validation_rows[[length(validation_rows) + 1L]] <- tibble::tibble(
          run_id = as.integer(run_row$run_id[[1]]),
          task_id = as.integer(run_row$task_id[[1]]),
          has_issues = TRUE,
          issues = err_msg
        )
        # Emit a stub run_metrics row so failed runs stay visible to aggregators
        # (otherwise downstream joins silently drop the locus).
        run_metrics_rows[[length(run_metrics_rows) + 1L]] <-
          real_data_failed_run_metrics_stub(run_row, err_msg)
        next
      }

      fit_path <- file.path(raw_fits_dir, sprintf("run-%05d_fit.rds", as.integer(run_row$run_id[[1]])))
      extracted <- tryCatch({
        saveRDS(fit_result, fit_path, compress = "xz")
        fit_index_rows[[length(fit_index_rows) + 1L]] <- real_data_build_fit_index_row(run_row, fit_path)

        wall_time <- proc.time()[["elapsed"]] - start_time
        extract_rss_fit_artifacts(
          fit = fit_result,
          run_row = run_row,
          bundle = bundle,
          job_config = job_config,
          fit_rds_path = fit_path,
          wall_time_sec = wall_time
        )
      }, error = function(e) {
        validation_rows[[length(validation_rows) + 1L]] <<- tibble::tibble(
          run_id = as.integer(run_row$run_id[[1]]),
          task_id = as.integer(run_row$task_id[[1]]),
          has_issues = TRUE,
          issues = conditionMessage(e)
        )
        NULL
      })
      if (is.null(extracted)) {
        next
      }

      run_metrics_rows[[length(run_metrics_rows) + 1L]] <- extracted$run_metrics
      effect_summary_rows[[length(effect_summary_rows) + 1L]] <- extracted$effect_summaries
      membership_rows[[length(membership_rows) + 1L]] <- extracted$credible_set_membership
      variant_rows[[length(variant_rows) + 1L]] <- extracted$variant_posteriors
      effect_rows[[length(effect_rows) + 1L]] <- extracted$effect_posteriors
      elbo_rows[[length(elbo_rows) + 1L]] <- extracted$elbo_trace
      validation_rows[[length(validation_rows) + 1L]] <- extracted$validation
    }

    run_metrics_tbl <- dplyr::bind_rows(run_metrics_rows)
    effect_summaries_tbl <- dplyr::bind_rows(effect_summary_rows)
    membership_tbl <- dplyr::bind_rows(membership_rows)
    variant_tbl <- dplyr::bind_rows(variant_rows)
    effect_tbl <- dplyr::bind_rows(effect_rows)
    elbo_tbl <- dplyr::bind_rows(elbo_rows)
    validation_tbl <- dplyr::bind_rows(validation_rows)

    if (nrow(run_metrics_tbl)) {
      readr::write_csv(run_metrics_tbl, file.path(task_dir, sprintf("%s_run_metrics.csv", flush_label)))
    }
    if (nrow(effect_summaries_tbl)) {
      readr::write_csv(effect_summaries_tbl, file.path(task_dir, sprintf("%s_effect_summaries.csv", flush_label)))
    }
    if (nrow(membership_tbl)) {
      real_data_write_parquet(membership_tbl, file.path(task_dir, sprintf("%s_credible_set_membership.parquet", flush_label)))
    }
    if (nrow(variant_tbl)) {
      real_data_write_parquet(variant_tbl, file.path(task_dir, sprintf("%s_variant_posteriors.parquet", flush_label)))
    }
    if (nrow(effect_tbl)) {
      real_data_write_parquet(effect_tbl, file.path(task_dir, sprintf("%s_effect_posteriors.parquet", flush_label)))
    }
    if (nrow(elbo_tbl)) {
      readr::write_csv(elbo_tbl, file.path(task_dir, sprintf("%s_elbo_trace.csv", flush_label)))
    }
    if (nrow(validation_tbl)) {
      readr::write_csv(validation_tbl, file.path(task_dir, sprintf("%s_validation.csv", flush_label)))
    }
    if (!dataset_metrics_written) {
      readr::write_csv(dataset_metrics_tbl, file.path(task_dir, sprintf("%s_dataset_metrics.csv", flush_label)))
      dataset_metrics_written <- TRUE
    }
  }

  fit_index_tbl <- dplyr::bind_rows(fit_index_rows)
  if (nrow(fit_index_tbl)) {
    readr::write_csv(fit_index_tbl, file.path(task_dir, "flush-001_fit_file_index.csv"))
  }

  invisible(
    list(
      task_dir = task_dir,
      raw_fits_dir = raw_fits_dir
    )
  )
}

compute_real_data_comparisons <- function(run_metrics_tbl, pip_lookup, run_rows, anchor_run_id = NA_integer_) {
  run_rows <- dplyr::arrange(run_rows, .data$run_id)
  functional_rows <- dplyr::filter(run_rows, .data$backend == "susine_rss", .data$run_family == "functional_grid")
  if (!nrow(functional_rows)) {
    return(tibble::tibble())
  }

  metric_lookup <- run_metrics_tbl %>%
    dplyr::select(
      .data$run_id,
      .data$elbo_final,
      .data$sigma_2_final,
      .data$h2_proxy_std,
      .data$pve_postmean_std,
      .data$pip_mass_top1,
      .data$pip_mass_top10
    )

  build_row <- function(run_id, target_run_id, comparison_target) {
    p <- pip_lookup[[as.character(run_id)]]
    q <- pip_lookup[[as.character(target_run_id)]]
    if (is.null(p) || is.null(q)) {
      return(NULL)
    }
    run_m <- dplyr::filter(metric_lookup, .data$run_id == !!run_id)
    target_m <- dplyr::filter(metric_lookup, .data$run_id == !!target_run_id)
    if (nrow(run_m) != 1L || nrow(target_m) != 1L) {
      return(NULL)
    }
    tibble::tibble(
      run_id = as.integer(run_id),
      target_run_id = as.integer(target_run_id),
      comparison_target = as.character(comparison_target),
      jsd_pip = as.numeric(js_distance(p, q)),
      l1_pip = sum(abs(p - q), na.rm = TRUE),
      spearman_pip = suppressWarnings(stats::cor(p, q, method = "spearman")),
      top10_overlap = real_data_top_overlap(p, q, 10L),
      top20_overlap = real_data_top_overlap(p, q, 20L),
      delta_elbo = run_m$elbo_final[[1]] - target_m$elbo_final[[1]],
      delta_sigma_2 = run_m$sigma_2_final[[1]] - target_m$sigma_2_final[[1]],
      delta_h2_proxy_std = run_m$h2_proxy_std[[1]] - target_m$h2_proxy_std[[1]],
      delta_pve_postmean_std = run_m$pve_postmean_std[[1]] - target_m$pve_postmean_std[[1]],
      delta_pip_mass_top1 = run_m$pip_mass_top1[[1]] - target_m$pip_mass_top1[[1]],
      delta_pip_mass_top10 = run_m$pip_mass_top10[[1]] - target_m$pip_mass_top10[[1]]
    )
  }

  rows <- list()
  for (i in seq_len(nrow(functional_rows))) {
    run_row <- functional_rows[i, , drop = FALSE]
    baseline_row <- dplyr::filter(
      functional_rows,
      .data$sigma_0_2_scalar == run_row$sigma_0_2_scalar[[1]],
      .data$c_value == 0
    )
    if (nrow(baseline_row) == 1L) {
      rows[[length(rows) + 1L]] <- build_row(
        run_id = run_row$run_id[[1]],
        target_run_id = baseline_row$run_id[[1]],
        comparison_target = "same_sigma_c0"
      )
    }
    if (is.finite(anchor_run_id)) {
      rows[[length(rows) + 1L]] <- build_row(
        run_id = run_row$run_id[[1]],
        target_run_id = anchor_run_id,
        comparison_target = "susie_anchor"
      )
    }
  }
  dplyr::bind_rows(rows)
}

real_data_write_dataset <- function(files, output_dir, partitioning = c("locus_id")) {
  if (!length(files)) {
    return(invisible(NULL))
  }
  unlink(output_dir, recursive = TRUE)
  ds <- arrow::open_dataset(files, format = "parquet")
  arrow::write_dataset(ds, output_dir, format = "parquet", partitioning = partitioning)
  invisible(output_dir)
}

collect_real_data_results <- function(
    job_name,
    parent_job_id,
    output_root = "output",
    validate = TRUE,
    output_dir = NULL
) {
  idx <- index_staging_outputs(job_name, parent_job_id, output_root = output_root)
  validation_index <- NULL
  if (validate) {
    validation_index <- validate_staging_outputs(idx)
    if (nrow(validation_index) && any(!validation_index$ok, na.rm = TRUE)) {
      bad <- dplyr::filter(validation_index, !.data$ok)
      stop("Real-data staging validation failed for ", nrow(bad), " file(s).")
    }
  }

  results_dir <- file.path(output_root, "slurm_output", job_name, parent_job_id)
  if (is.null(output_dir)) {
    output_dir <- file.path(results_dir, "aggregated")
  }
  ensure_dir(output_dir)
  highest_weight_refit_dir <- file.path(output_dir, "highest_weight_refits")
  ensure_dir(highest_weight_refit_dir)

  temp_dir <- file.path(output_root, "temp", job_name)
  history_dir <- file.path(output_root, "run_history", job_name, parent_job_id)
  job_config_path <- if (file.exists(file.path(temp_dir, "job_config.json"))) {
    file.path(temp_dir, "job_config.json")
  } else {
    file.path(history_dir, "job_config.json")
  }
  run_manifest_path <- if (file.exists(file.path(temp_dir, "run_manifest.csv"))) {
    file.path(temp_dir, "run_manifest.csv")
  } else {
    file.path(history_dir, "run_manifest.csv")
  }
  dataset_bundles_path <- if (file.exists(file.path(temp_dir, "dataset_bundles.csv"))) {
    file.path(temp_dir, "dataset_bundles.csv")
  } else {
    file.path(history_dir, "dataset_bundles.csv")
  }

  job_config <- if (file.exists(job_config_path)) {
    jsonlite::read_json(job_config_path, simplifyVector = TRUE)
  } else {
    NULL
  }
  jsd_threshold <- as.numeric(job_config$job$jsd_threshold %||% 0.15)
  softmax_temperature <- as.numeric(job_config$job$softmax_temperature %||% 1)

  run_manifest <- readr::read_csv(run_manifest_path, show_col_types = FALSE)
  dataset_bundles <- readr::read_csv(dataset_bundles_path, show_col_types = FALSE)
  readr::write_csv(run_manifest, file.path(output_dir, "run_manifest.csv"))

  csv_types <- c("run_metrics", "effect_summaries", "dataset_metrics", "validation", "elbo_trace", "fit_file_index")
  csv_out_names <- c(
    run_metrics = "run_metrics_full.csv",
    effect_summaries = "effect_summaries_full.csv",
    dataset_metrics = "dataset_metrics.csv",
    validation = "validation.csv",
    elbo_trace = "elbo_trace_full.csv",
    fit_file_index = "fit_file_index.csv"
  )
  csv_tables <- list()
  for (type in csv_types) {
    files <- dplyr::filter(idx, .data$type == !!type)$path
    tbl <- real_data_bind_csv_files(files)
    if (identical(type, "dataset_metrics") && nrow(tbl)) {
      tbl <- dplyr::distinct(tbl, .data$dataset_bundle_id, .keep_all = TRUE)
    }
    csv_tables[[type]] <- tbl
    if (nrow(tbl)) {
      readr::write_csv(tbl, file.path(output_dir, csv_out_names[[type]]))
    }
  }

  variant_files <- dplyr::filter(idx, .data$type == "variant_posteriors")$path
  effect_files <- dplyr::filter(idx, .data$type == "effect_posteriors")$path
  cs_files <- dplyr::filter(idx, .data$type == "credible_set_membership")$path

  variant_dataset_dir <- file.path(output_dir, "variant_posteriors_dataset")
  effect_dataset_dir <- file.path(output_dir, "effect_posteriors_dataset")
  cs_dataset_dir <- file.path(output_dir, "credible_set_membership_dataset")

  real_data_write_dataset(variant_files, variant_dataset_dir)
  real_data_write_dataset(effect_files, effect_dataset_dir)
  real_data_write_dataset(cs_files, cs_dataset_dir)

  variant_ds <- if (length(variant_files)) arrow::open_dataset(variant_dataset_dir) else NULL
  pairwise_rows <- list()
  multimodal_rows <- list()
  cluster_rows <- list()
  cluster_summary_rows <- list()
  agg_weight_rows <- list()
  aggregated_variant_rows <- list()
  top_variant_rows <- list()
  comparison_rows <- list()
  functional_summary_rows <- list()
  anchor_summary_rows <- list()
  highest_weight_refit_summary_rows <- list()
  highest_weight_refit_variant_rows <- list()
  highest_weight_refit_basin_rows <- list()
  paper_summary_rows <- list()

  run_metrics_tbl <- csv_tables$run_metrics
  if (!"pve_postmean_std" %in% names(run_metrics_tbl)) {
    run_metrics_tbl$pve_postmean_std <- NA_real_
  }
  next_refit_run_id <- max(as.integer(run_manifest$run_id), na.rm = TRUE) + 1L

  for (i in seq_len(nrow(dataset_bundles))) {
    locus_row <- dataset_bundles[i, , drop = FALSE]
    locus_id <- locus_row$locus_id[[1]]
    functional_runs <- run_manifest %>%
      dplyr::filter(
        .data$locus_id == !!locus_id,
        .data$backend == "susine_rss",
        .data$run_family == "functional_grid"
      ) %>%
      dplyr::arrange(.data$run_id)
    anchor_runs <- run_manifest %>%
      dplyr::filter(
        .data$locus_id == !!locus_id,
        .data$backend == "susie_rss",
        .data$run_family == "susie_anchor"
      ) %>%
      dplyr::arrange(.data$run_id)

    if (!nrow(functional_runs) || is.null(variant_ds)) {
      next
    }

    locus_variant_tbl <- variant_ds %>%
      dplyr::filter(.data$locus_id == !!locus_id) %>%
      dplyr::select(
        .data$run_id, .data$backend, .data$run_family, .data$variant_id,
        .data$ld_matrix_index, .data$pip, .data$z_score, .data$beta_hat_std,
        .data$annotation_a, .data$baseline_c_l, .data$posterior_mean
      ) %>%
      dplyr::collect()

    locus_bundle_for_pve <- tryCatch(
      load_real_data_locus_bundle(
        locus_id = locus_id,
        manifest_path = job_config$job$manifest_path,
        repo_root = job_config$paths$repo_root %||% ensure_repo_root(getwd())
      ),
      error = function(e) NULL
    )
    if (!is.null(locus_bundle_for_pve) && "posterior_mean" %in% names(locus_variant_tbl)) {
      pve_by_run <- locus_variant_tbl %>%
        dplyr::arrange(.data$run_id, .data$ld_matrix_index) %>%
        dplyr::group_by(.data$run_id) %>%
        dplyr::summarise(
          pve_postmean_std = real_data_pve_from_posterior_mean(
            .data$posterior_mean,
            locus_bundle_for_pve$R
          ),
          .groups = "drop"
        )
      if (!"pve_postmean_std" %in% names(run_metrics_tbl)) {
        run_metrics_tbl$pve_postmean_std <- NA_real_
      }
      idx_run <- match(pve_by_run$run_id, run_metrics_tbl$run_id)
      keep_idx <- !is.na(idx_run)
      run_metrics_tbl$pve_postmean_std[idx_run[keep_idx]] <- pve_by_run$pve_postmean_std[keep_idx]
    }

    locus_functional_tbl <- locus_variant_tbl %>%
      dplyr::filter(.data$backend == "susine_rss", .data$run_family == "functional_grid") %>%
      dplyr::arrange(.data$run_id, .data$ld_matrix_index)
    pip_lookup_all <- split(
      locus_variant_tbl %>%
        dplyr::arrange(.data$run_id, .data$ld_matrix_index) %>%
        dplyr::pull(.data$pip),
      locus_variant_tbl %>%
        dplyr::arrange(.data$run_id, .data$ld_matrix_index) %>%
        dplyr::pull(.data$run_id)
    )

    pip_split <- split(locus_functional_tbl$pip, locus_functional_tbl$run_id)
    run_ids <- as.integer(names(pip_split))
    functional_runs_present <- functional_runs %>%
      dplyr::filter(.data$run_id %in% run_ids) %>%
      dplyr::arrange(match(.data$run_id, run_ids))
    pip_mat <- do.call(rbind, lapply(pip_split, as.numeric))
    rownames(pip_mat) <- run_ids
    locus_run_metrics_present <- run_metrics_tbl %>%
      dplyr::filter(.data$run_id %in% run_ids) %>%
      dplyr::distinct(.data$run_id, .keep_all = TRUE) %>%
      dplyr::arrange(match(.data$run_id, run_ids)) %>%
      tibble::as_tibble()
    if (nrow(locus_run_metrics_present) != length(run_ids)) {
      stop("Missing run_metrics rows needed to collect functional runs for locus ", locus_id)
    }
    elbo_vec <- locus_run_metrics_present$elbo_final

    pip_cache <- prepare_pip_similarity_cache(pip_mat)
    mm <- compute_multimodal_metrics(
      pip_list = lapply(seq_len(nrow(pip_mat)), function(j) pip_mat[j, ]),
      jsd_threshold = jsd_threshold,
      pip_cache = pip_cache
    ) %>%
      dplyr::mutate(
        dataset_bundle_id = locus_row$dataset_bundle_id[[1]],
        locus_id = locus_id,
        gene_name = locus_row$gene_name[[1]],
        group_key = paste0("functional_grid|", locus_id)
      )
    multimodal_rows[[length(multimodal_rows) + 1L]] <- mm

    if (nrow(pip_mat) > 1L) {
      clusters <- stats::cutree(pip_cache$hc, h = jsd_threshold)
      cw <- .cluster_weights_from_hc(
        pip_cache$hc,
        jsd_threshold,
        elbo_vec,
        n_fits = nrow(pip_mat),
        softmax_temperature = softmax_temperature
      )
      cluster_levels <- sort(unique(clusters))
      weight_by_run <- rep(0, length(run_ids))
      weight_by_run[cw$rep_idx] <- cw$w_rep

      pair_rows <- list()
      for (a in seq_len(nrow(pip_mat) - 1L)) {
        for (b in seq.int(a + 1L, nrow(pip_mat))) {
          pair_rows[[length(pair_rows) + 1L]] <- tibble::tibble(
            locus_id = locus_id,
            gene_name = locus_row$gene_name[[1]],
            run_id_i = as.integer(run_ids[[a]]),
            run_id_j = as.integer(run_ids[[b]]),
            c_i = functional_runs_present$c_value[match(run_ids[[a]], functional_runs_present$run_id)],
            sigma_i = functional_runs_present$sigma_0_2_scalar[match(run_ids[[a]], functional_runs_present$run_id)],
            c_j = functional_runs_present$c_value[match(run_ids[[b]], functional_runs_present$run_id)],
            sigma_j = functional_runs_present$sigma_0_2_scalar[match(run_ids[[b]], functional_runs_present$run_id)],
            jsd = as.numeric(pip_cache$jsd_mat[a, b]),
            same_cluster = as.logical(clusters[[a]] == clusters[[b]])
          )
        }
      }
      pairwise_rows[[length(pairwise_rows) + 1L]] <- dplyr::bind_rows(pair_rows)

      membership_tbl <- purrr::map_dfr(seq_along(run_ids), function(j) {
        cluster_id <- clusters[[j]]
        cluster_members <- which(clusters == cluster_id)
        rep_local_idx <- cw$rep_idx[match(cluster_id, cluster_levels)]
        tibble::tibble(
          locus_id = locus_id,
          gene_name = locus_row$gene_name[[1]],
          run_id = as.integer(run_ids[[j]]),
          cluster_id = as.integer(cluster_id),
          cluster_size = length(cluster_members),
          cluster_freq = length(cluster_members) / length(run_ids),
          is_representative = as.logical(j == rep_local_idx),
          representative_run_id = as.integer(run_ids[[rep_local_idx]]),
          representative_elbo = as.numeric(elbo_vec[[rep_local_idx]]),
          jsd_to_representative = as.numeric(if (j == rep_local_idx) 0 else pip_cache$jsd_mat[j, rep_local_idx])
        )
      })
      cluster_rows[[length(cluster_rows) + 1L]] <- membership_tbl

      agg_weights_locus <- membership_tbl %>%
        dplyr::left_join(
          functional_runs_present %>% dplyr::select(.data$run_id, .data$c_value, .data$sigma_0_2_scalar),
          by = "run_id"
        ) %>%
        dplyr::mutate(
          agg_weight_run = weight_by_run[match(.data$run_id, run_ids)]
        ) %>%
        dplyr::select(
          .data$locus_id, .data$gene_name, .data$run_id, .data$c_value,
          .data$sigma_0_2_scalar, .data$cluster_id, .data$agg_weight_run
        )
      agg_weight_rows[[length(agg_weight_rows) + 1L]] <- agg_weights_locus

      cluster_summary_rows[[length(cluster_summary_rows) + 1L]] <- purrr::map_dfr(seq_along(cluster_levels), function(k) {
        cid <- cluster_levels[[k]]
        rep_idx <- cw$rep_idx[[k]]
        members <- which(clusters == cid)
        tibble::tibble(
          locus_id = locus_id,
          gene_name = locus_row$gene_name[[1]],
          cluster_id = as.integer(cid),
          representative_run_id = as.integer(run_ids[[rep_idx]]),
          cluster_size = length(members),
          cluster_freq = length(members) / length(run_ids),
          representative_elbo = as.numeric(elbo_vec[[rep_idx]]),
          cluster_weight = as.numeric(cw$w_rep[[k]])
        )
      })

      agg_pip <- aggregate_pip_matrix(
        pip_mat = t(pip_mat),
        elbos = elbo_vec,
        method = "cluster_weight",
        jsd_threshold = jsd_threshold,
        hc = pip_cache$hc
      )
    } else {
      membership_tbl <- tibble::tibble(
        locus_id = locus_id,
        gene_name = locus_row$gene_name[[1]],
        run_id = as.integer(run_ids[[1]]),
        cluster_id = 1L,
        cluster_size = 1L,
        cluster_freq = 1,
        is_representative = TRUE,
        representative_run_id = as.integer(run_ids[[1]]),
        representative_elbo = as.numeric(elbo_vec[[1]]),
        jsd_to_representative = 0
      )
      cluster_rows[[length(cluster_rows) + 1L]] <- membership_tbl
      agg_weights_locus <- membership_tbl %>%
        dplyr::left_join(
          functional_runs_present %>% dplyr::select(.data$run_id, .data$c_value, .data$sigma_0_2_scalar),
          by = "run_id"
        ) %>%
        dplyr::mutate(agg_weight_run = 1) %>%
        dplyr::select(
          .data$locus_id, .data$gene_name, .data$run_id, .data$c_value,
          .data$sigma_0_2_scalar, .data$cluster_id, .data$agg_weight_run
        )
      agg_weight_rows[[length(agg_weight_rows) + 1L]] <- agg_weights_locus
      cluster_summary_rows[[length(cluster_summary_rows) + 1L]] <- tibble::tibble(
        locus_id = locus_id,
        gene_name = locus_row$gene_name[[1]],
        cluster_id = 1L,
        representative_run_id = as.integer(run_ids[[1]]),
        cluster_size = 1L,
        cluster_freq = 1,
        representative_elbo = as.numeric(elbo_vec[[1]]),
        cluster_weight = 1
      )
      agg_pip <- as.numeric(pip_mat[1, ])
    }

    weight_vec <- agg_weights_locus$agg_weight_run[match(run_ids, agg_weights_locus$run_id)]
    weight_vec[!is.finite(weight_vec)] <- 0
    pm_split <- split(locus_functional_tbl$posterior_mean, locus_functional_tbl$run_id)
    pm_mat <- do.call(rbind, lapply(as.character(run_ids), function(rid) as.numeric(pm_split[[rid]])))
    agg_posterior_mean <- if (length(weight_vec) == nrow(pm_mat) && sum(weight_vec) > 0) {
      as.numeric(crossprod(weight_vec / sum(weight_vec), pm_mat))
    } else {
      rep(NA_real_, ncol(pm_mat))
    }

    ensemble_pve_postmean_std <- if (!is.null(locus_bundle_for_pve)) {
      real_data_pve_from_posterior_mean(agg_posterior_mean, locus_bundle_for_pve$R)
    } else {
      NA_real_
    }

    base_variant_tbl <- locus_functional_tbl %>%
      dplyr::filter(.data$run_id == run_ids[[1]]) %>%
      dplyr::arrange(.data$ld_matrix_index)
    aggregated_tbl <- base_variant_tbl %>%
      dplyr::transmute(
        locus_id = locus_id,
        gene_name = locus_row$gene_name[[1]],
        variant_id = as.character(.data$variant_id),
        ld_matrix_index = as.integer(.data$ld_matrix_index),
        aggregated_pip = as.numeric(agg_pip),
        aggregated_posterior_mean = as.numeric(agg_posterior_mean),
        z_score = as.numeric(.data$z_score),
        beta_hat_std = as.numeric(.data$beta_hat_std),
        annotation_a = as.numeric(.data$annotation_a),
        baseline_c_l = as.numeric(.data$baseline_c_l)
      )
    aggregated_variant_rows[[length(aggregated_variant_rows) + 1L]] <- aggregated_tbl
    top_variant_rows[[length(top_variant_rows) + 1L]] <- aggregated_tbl %>%
      dplyr::arrange(dplyr::desc(.data$aggregated_pip), .data$ld_matrix_index) %>%
      dplyr::slice_head(n = 20L) %>%
      dplyr::mutate(rank = dplyr::row_number(), .before = 1L)

    highest_weight_source <- agg_weights_locus %>%
      dplyr::left_join(
        locus_run_metrics_present %>%
          dplyr::select(.data$run_id, .data$fit_rds_path, .data$elbo_final,
                        .data$pve_postmean_std, .data$max_pip),
        by = "run_id"
      ) %>%
      dplyr::arrange(
        dplyr::desc(.data$agg_weight_run),
        dplyr::desc(.data$elbo_final),
        .data$run_id
      ) %>%
      dplyr::slice_head(n = 1L)

    refit_fit <- NULL
    source_fit <- NULL
    refit_artifacts <- NULL
    refit_metric <- tibble::tibble()
    refit_variant_tbl <- tibble::tibble()
    refit_pip <- NULL
    refit_error <- NA_character_
    refit_run_id <- next_refit_run_id
    next_refit_run_id <- next_refit_run_id + 1L

    source_fit_path <- if (nrow(highest_weight_source)) highest_weight_source$fit_rds_path[[1]] else NA_character_
    can_refit <- !is.null(locus_bundle_for_pve) &&
      nrow(highest_weight_source) == 1L &&
      is.character(source_fit_path) &&
      !is.na(source_fit_path) &&
      file.exists(source_fit_path)

    if (can_refit) {
      refit_locus_dir <- file.path(highest_weight_refit_dir, locus_id)
      ensure_dir(refit_locus_dir)
      refit_fit_path <- file.path(
        refit_locus_dir,
        paste0("highest_weight_refit_run_", refit_run_id, ".rds")
      )
      source_fit <- tryCatch(readRDS(source_fit_path), error = function(e) {
        refit_error <<- conditionMessage(e)
        NULL
      })
      if (!is.null(source_fit)) {
        refit_run_row <- tibble::tibble(
          run_id = as.integer(refit_run_id),
          task_id = NA_integer_,
          dataset_bundle_id = as.integer(locus_row$dataset_bundle_id[[1]]),
          locus_id = as.character(locus_id),
          gene_name = as.character(locus_row$gene_name[[1]]),
          backend = "susine_rss",
          run_family = "highest_weight_warm_refit",
          c_value = 0,
          sigma_0_2_scalar = 0.2,
          baseline_c_l = as.numeric((locus_row$baseline_c_l %||% NA_real_)[[1]])
        )
        elapsed <- system.time({
          refit_fit <- tryCatch(
            real_data_run_highest_weight_refit(
              source_fit = source_fit,
              bundle = locus_bundle_for_pve,
              job_config = job_config,
              sigma_0_2_scalar = 0.2
            ),
            error = function(e) {
              refit_error <<- conditionMessage(e)
              NULL
            }
          )
        })[["elapsed"]]
        if (!is.null(refit_fit)) {
          saveRDS(refit_fit, refit_fit_path)
          refit_artifacts <- extract_rss_fit_artifacts(
            fit = refit_fit,
            run_row = refit_run_row,
            bundle = locus_bundle_for_pve,
            job_config = job_config,
            fit_rds_path = refit_fit_path,
            wall_time_sec = elapsed
          )
          refit_metric <- refit_artifacts$run_metrics %>%
            dplyr::mutate(
              source_run_id = as.integer(highest_weight_source$run_id[[1]]),
              source_agg_weight_run = as.numeric(highest_weight_source$agg_weight_run[[1]]),
              source_c_value = as.numeric(highest_weight_source$c_value[[1]]),
              source_sigma_0_2_scalar = as.numeric(highest_weight_source$sigma_0_2_scalar[[1]]),
              source_fit_rds_path = as.character(source_fit_path)
            )
          refit_variant_tbl <- refit_artifacts$variant_posteriors %>%
            dplyr::mutate(
              source_run_id = as.integer(highest_weight_source$run_id[[1]]),
              source_agg_weight_run = as.numeric(highest_weight_source$agg_weight_run[[1]])
            )
          refit_pip <- refit_variant_tbl %>%
            dplyr::arrange(.data$ld_matrix_index) %>%
            dplyr::pull(.data$pip)
          highest_weight_refit_variant_rows[[length(highest_weight_refit_variant_rows) + 1L]] <- refit_variant_tbl
        }
      }
    } else if (!nrow(highest_weight_source)) {
      refit_error <- "No highest-weight source run was available."
    } else if (is.null(locus_bundle_for_pve)) {
      refit_error <- "Locus bundle could not be loaded."
    } else {
      refit_error <- paste0("Highest-weight source fit is missing: ", source_fit_path)
    }

    highest_weight_refit_summary_rows[[length(highest_weight_refit_summary_rows) + 1L]] <- tibble::tibble(
      locus_id = locus_id,
      gene_name = locus_row$gene_name[[1]],
      refit_status = if (is.null(refit_fit)) "failed" else "ok",
      refit_error = if (is.null(refit_fit)) refit_error else NA_character_,
      refit_run_id = as.integer(refit_run_id),
      refit_fit_rds_path = if (nrow(refit_metric)) refit_metric$fit_rds_path[[1]] else NA_character_,
      refit_elbo_final = if (nrow(refit_metric)) refit_metric$elbo_final[[1]] else NA_real_,
      refit_pve_postmean_std = if (nrow(refit_metric)) refit_metric$pve_postmean_std[[1]] else NA_real_,
      refit_max_pip = if (nrow(refit_metric)) refit_metric$max_pip[[1]] else NA_real_,
      source_run_id = if (nrow(highest_weight_source)) as.integer(highest_weight_source$run_id[[1]]) else NA_integer_,
      source_agg_weight_run = if (nrow(highest_weight_source)) highest_weight_source$agg_weight_run[[1]] else NA_real_,
      source_c_value = if (nrow(highest_weight_source)) highest_weight_source$c_value[[1]] else NA_real_,
      source_sigma_0_2_scalar = if (nrow(highest_weight_source)) highest_weight_source$sigma_0_2_scalar[[1]] else NA_real_,
      source_elbo_final = if (nrow(highest_weight_source)) highest_weight_source$elbo_final[[1]] else NA_real_,
      source_pve_postmean_std = if (nrow(highest_weight_source)) highest_weight_source$pve_postmean_std[[1]] else NA_real_,
      source_max_pip = if (nrow(highest_weight_source)) highest_weight_source$max_pip[[1]] else NA_real_,
      source_fit_rds_path = source_fit_path
    )

    anchor_run_id <- if (nrow(anchor_runs)) as.integer(anchor_runs$run_id[[1]]) else NA_integer_
    locus_run_metrics <- run_metrics_tbl %>% dplyr::filter(.data$locus_id == !!locus_id)
    comparison_tbl <- compute_real_data_comparisons(
      run_metrics_tbl = locus_run_metrics,
      pip_lookup = pip_lookup_all,
      run_rows = functional_runs_present,
      anchor_run_id = anchor_run_id
    ) %>%
      dplyr::mutate(locus_id = locus_id, gene_name = locus_row$gene_name[[1]], .before = 1L)
    comparison_rows[[length(comparison_rows) + 1L]] <- comparison_tbl

    same_sigma_tbl <- comparison_tbl %>%
      dplyr::filter(.data$comparison_target == "same_sigma_c0") %>%
      dplyr::rename(
        jsd_from_same_sigma_c0 = .data$jsd_pip,
        l1_from_same_sigma_c0 = .data$l1_pip,
        spearman_from_same_sigma_c0 = .data$spearman_pip,
        top10_overlap_same_sigma_c0 = .data$top10_overlap,
        top20_overlap_same_sigma_c0 = .data$top20_overlap,
        delta_elbo_same_sigma_c0 = .data$delta_elbo,
        delta_sigma_2_same_sigma_c0 = .data$delta_sigma_2,
        delta_h2_proxy_std_same_sigma_c0 = .data$delta_h2_proxy_std,
        delta_pve_postmean_std_same_sigma_c0 = .data$delta_pve_postmean_std,
        delta_pip_mass_top1_same_sigma_c0 = .data$delta_pip_mass_top1,
        delta_pip_mass_top10_same_sigma_c0 = .data$delta_pip_mass_top10
      )

    functional_summary_rows[[length(functional_summary_rows) + 1L]] <- functional_runs %>%
      dplyr::filter(.data$run_id %in% run_ids) %>%
      dplyr::arrange(match(.data$run_id, run_ids)) %>%
      dplyr::left_join(
        dplyr::select(
          locus_run_metrics_present,
          .data$run_id,
          dplyr::all_of(setdiff(names(locus_run_metrics_present), c(names(functional_runs_present), "locus_id", "gene_name")))
        ),
        by = "run_id"
      ) %>%
      dplyr::left_join(
        membership_tbl %>% dplyr::select(
          .data$run_id, .data$cluster_id,
          .data$cluster_size, .data$cluster_freq, .data$is_representative
        ),
        by = "run_id"
      ) %>%
      dplyr::left_join(
        agg_weights_locus %>% dplyr::select(
          .data$run_id, .data$agg_weight_run
        ),
        by = "run_id"
      ) %>%
      dplyr::left_join(
        same_sigma_tbl %>% dplyr::select(-.data$comparison_target, -.data$gene_name, -.data$locus_id),
        by = "run_id"
      )

    if (nrow(anchor_runs)) {
      anchor_id <- as.integer(anchor_runs$run_id[[1]])
      anchor_metric <- locus_run_metrics %>% dplyr::filter(.data$run_id == anchor_id)
      baseline_metric <- locus_run_metrics %>%
        dplyr::filter(.data$backend == "susine_rss", .data$run_family == "functional_grid",
                      abs(.data$sigma_0_2_scalar - 0.2) < 1e-12,
                      abs(.data$c_value) < 1e-12)
      # Emit a row whenever the anchor was attempted; mark status so notebooks
      # can distinguish ok vs failed without losing the locus.
      if (nrow(anchor_metric) == 1L) {
        anchor_failed <- isTRUE(is.na(anchor_metric$elbo_final[[1]]))
        anchor_err <- if ("error_message" %in% names(anchor_metric))
          anchor_metric$error_message[[1]] else NA_character_
        base_ok <- nrow(baseline_metric) == 1L &&
          isTRUE(is.finite(baseline_metric$elbo_final[[1]]))
        anchor_pip <- if (!anchor_failed) {
          locus_variant_tbl %>% dplyr::filter(.data$run_id == anchor_id) %>%
            dplyr::arrange(.data$ld_matrix_index) %>% dplyr::pull(.data$pip)
        } else NULL
        base_pip <- if (base_ok) {
          locus_variant_tbl %>%
            dplyr::filter(.data$run_id == baseline_metric$run_id[[1]]) %>%
            dplyr::arrange(.data$ld_matrix_index) %>% dplyr::pull(.data$pip)
        } else NULL
        can_compare <- !is.null(anchor_pip) && !is.null(base_pip) &&
          length(anchor_pip) == length(base_pip) && length(anchor_pip) > 0L

        anchor_summary_rows[[length(anchor_summary_rows) + 1L]] <- tibble::tibble(
          locus_id = locus_id,
          gene_name = locus_row$gene_name[[1]],
          susie_anchor_run_id = anchor_id,
          susie_anchor_status = if (anchor_failed) "failed" else "ok",
          susie_anchor_error = if (anchor_failed) as.character(anchor_err) else NA_character_,
          susie_anchor_sigma_0_2_scalar = anchor_metric$sigma_0_2_scalar[[1]],
          susie_anchor_elbo_final = anchor_metric$elbo_final[[1]],
          susie_anchor_sigma_2_final = anchor_metric$sigma_2_final[[1]],
          susie_anchor_h2_proxy_std = anchor_metric$h2_proxy_std[[1]],
          susie_anchor_pve_postmean_std = anchor_metric$pve_postmean_std[[1]],
          susie_anchor_max_pip = anchor_metric$max_pip[[1]],
          baseline_run_id = if (base_ok) baseline_metric$run_id[[1]] else NA_integer_,
          baseline_elbo_final = if (base_ok) baseline_metric$elbo_final[[1]] else NA_real_,
          baseline_sigma_2_final = if (base_ok) baseline_metric$sigma_2_final[[1]] else NA_real_,
          baseline_h2_proxy_std = if (base_ok) baseline_metric$h2_proxy_std[[1]] else NA_real_,
          baseline_pve_postmean_std = if (base_ok) baseline_metric$pve_postmean_std[[1]] else NA_real_,
          ensemble_pve_postmean_std = ensemble_pve_postmean_std,
          highest_weight_refit_run_id = refit_run_id,
          highest_weight_refit_status = if (is.null(refit_fit)) "failed" else "ok",
          highest_weight_refit_pve_postmean_std = if (nrow(refit_metric)) refit_metric$pve_postmean_std[[1]] else NA_real_,
          highest_weight_refit_max_pip = if (nrow(refit_metric)) refit_metric$max_pip[[1]] else NA_real_,
          highest_weight_source_run_id = if (nrow(highest_weight_source)) highest_weight_source$run_id[[1]] else NA_integer_,
          highest_weight_source_agg_weight_run = if (nrow(highest_weight_source)) highest_weight_source$agg_weight_run[[1]] else NA_real_,
          baseline_max_pip = if (base_ok) baseline_metric$max_pip[[1]] else NA_real_,
          jsd_anchor_vs_baseline = if (can_compare)
            as.numeric(js_distance(anchor_pip, base_pip)) else NA_real_,
          l1_anchor_vs_baseline = if (can_compare)
            sum(abs(anchor_pip - base_pip), na.rm = TRUE) else NA_real_,
          top10_overlap_anchor_vs_baseline = if (can_compare)
            real_data_top_overlap(anchor_pip, base_pip, 10L) else NA_real_
        )

        source_run_id <- if (nrow(highest_weight_source)) as.integer(highest_weight_source$run_id[[1]]) else NA_integer_
        source_pip <- if (!is.na(source_run_id)) {
          pip_lookup_all[[as.character(source_run_id)]]
        } else {
          NULL
        }
        ensemble_weight_c0 <- agg_weights_locus %>%
          dplyr::filter(is.finite(.data$c_value), abs(.data$c_value) < 1e-12) %>%
          dplyr::summarise(weight = sum(.data$agg_weight_run, na.rm = TRUE), .groups = "drop") %>%
          dplyr::pull(.data$weight)
        ensemble_weight_c0 <- if (length(ensemble_weight_c0)) ensemble_weight_c0[[1]] else NA_real_
        ensemble_weight_off_c0 <- if (is.finite(ensemble_weight_c0)) {
          1 - ensemble_weight_c0
        } else {
          NA_real_
        }
        paper_summary_rows[[length(paper_summary_rows) + 1L]] <- tibble::tibble(
          locus_id = locus_id,
          gene_name = locus_row$gene_name[[1]],
          n_variants = length(agg_pip),
          susie_anchor_run_id = anchor_id,
          baseline_c0_run_id = if (base_ok) baseline_metric$run_id[[1]] else NA_integer_,
          highest_weight_source_run_id = source_run_id,
          highest_weight_refit_run_id = refit_run_id,
          highest_weight_source_agg_weight_run = if (nrow(highest_weight_source)) highest_weight_source$agg_weight_run[[1]] else NA_real_,
          ensemble_agg_weight_c0_total = ensemble_weight_c0,
          ensemble_agg_weight_off_c0_total = ensemble_weight_off_c0,
          highest_weight_source_c_value = if (nrow(highest_weight_source)) highest_weight_source$c_value[[1]] else NA_real_,
          highest_weight_source_sigma_0_2_scalar = if (nrow(highest_weight_source)) highest_weight_source$sigma_0_2_scalar[[1]] else NA_real_,
          elbo_susie_anchor = anchor_metric$elbo_final[[1]],
          elbo_baseline_c0 = if (base_ok) baseline_metric$elbo_final[[1]] else NA_real_,
          elbo_highest_weight_source = if (nrow(highest_weight_source)) highest_weight_source$elbo_final[[1]] else NA_real_,
          elbo_highest_weight_refit = if (nrow(refit_metric)) refit_metric$elbo_final[[1]] else NA_real_,
          delta_elbo_refit_vs_susie_anchor = if (nrow(refit_metric))
            refit_metric$elbo_final[[1]] - anchor_metric$elbo_final[[1]] else NA_real_,
          delta_elbo_refit_vs_baseline_c0 = if (nrow(refit_metric) && base_ok)
            refit_metric$elbo_final[[1]] - baseline_metric$elbo_final[[1]] else NA_real_,
          delta_elbo_source_vs_baseline_c0 = if (nrow(highest_weight_source) && base_ok)
            highest_weight_source$elbo_final[[1]] - baseline_metric$elbo_final[[1]] else NA_real_,
          n_pip_gt_025_susie_anchor = if (!is.null(anchor_pip)) real_data_pip_count(anchor_pip) else NA_integer_,
          n_pip_gt_025_baseline_c0 = if (!is.null(base_pip)) real_data_pip_count(base_pip) else NA_integer_,
          n_pip_gt_025_susine_ensemble = real_data_pip_count(agg_pip),
          n_pip_gt_025_highest_weight_source = if (!is.null(source_pip)) real_data_pip_count(source_pip) else NA_integer_,
          n_pip_gt_025_highest_weight_refit = if (!is.null(refit_pip)) real_data_pip_count(refit_pip) else NA_integer_,
          n_pip_gt_025_susie_anchor_and_susine_ensemble = real_data_pip_intersection_count(anchor_pip, agg_pip),
          n_pip_gt_025_susie_anchor_and_refit = real_data_pip_intersection_count(anchor_pip, refit_pip),
          n_pip_gt_025_baseline_c0_and_susine_ensemble = real_data_pip_intersection_count(base_pip, agg_pip),
          n_pip_gt_025_baseline_c0_and_refit = real_data_pip_intersection_count(base_pip, refit_pip),
          n_pip_gt_025_susine_ensemble_and_refit = real_data_pip_intersection_count(agg_pip, refit_pip),
          pve_susie_anchor = anchor_metric$pve_postmean_std[[1]],
          pve_baseline_c0 = if (base_ok) baseline_metric$pve_postmean_std[[1]] else NA_real_,
          pve_susine_ensemble = ensemble_pve_postmean_std,
          pve_highest_weight_source = if (nrow(highest_weight_source)) highest_weight_source$pve_postmean_std[[1]] else NA_real_,
          pve_highest_weight_refit = if (nrow(refit_metric)) refit_metric$pve_postmean_std[[1]] else NA_real_,
          jsd_susie_anchor_vs_susine_ensemble = real_data_safe_js_distance(anchor_pip, agg_pip),
          jsd_susie_anchor_vs_refit = real_data_safe_js_distance(anchor_pip, refit_pip),
          jsd_baseline_c0_vs_susine_ensemble = real_data_safe_js_distance(base_pip, agg_pip),
          jsd_baseline_c0_vs_refit = real_data_safe_js_distance(base_pip, refit_pip),
          jsd_susine_ensemble_vs_refit = real_data_safe_js_distance(agg_pip, refit_pip)
        )

        if (!is.null(locus_bundle_for_pve) && !anchor_failed && nrow(highest_weight_source)) {
          anchor_fit <- tryCatch(readRDS(anchor_metric$fit_rds_path[[1]]), error = function(e) NULL)
          if (is.null(source_fit) && !is.na(source_fit_path) && file.exists(source_fit_path)) {
            source_fit <- tryCatch(readRDS(source_fit_path), error = function(e) NULL)
          }
          anchor_cache <- real_data_basin_r2_cache(anchor_fit, locus_bundle_for_pve$R)
          source_cache <- real_data_basin_r2_cache(source_fit, locus_bundle_for_pve$R)
          refit_cache <- real_data_basin_r2_cache(refit_fit, locus_bundle_for_pve$R)
          source_vs_anchor <- real_data_basin_r2_drift_pair_cached(source_cache, anchor_cache)
          refit_vs_anchor <- real_data_basin_r2_drift_pair_cached(refit_cache, anchor_cache)
          refit_vs_source <- real_data_basin_r2_drift_pair_cached(refit_cache, source_cache)
          highest_weight_refit_basin_rows[[length(highest_weight_refit_basin_rows) + 1L]] <- tibble::tibble(
            locus_id = locus_id,
            gene_name = locus_row$gene_name[[1]],
            susie_anchor_run_id = anchor_id,
            highest_weight_source_run_id = source_run_id,
            highest_weight_refit_run_id = refit_run_id,
            basin_r2_drift_highest_weight_vs_susie_anchor = source_vs_anchor$drift,
            basin_r2_drift_refit_vs_susie_anchor = refit_vs_anchor$drift,
            basin_r2_drift_refit_vs_highest_weight = refit_vs_source$drift,
            n_matched_highest_weight_vs_susie_anchor = source_vs_anchor$n_matched,
            n_matched_refit_vs_susie_anchor = refit_vs_anchor$n_matched,
            n_matched_refit_vs_highest_weight = refit_vs_source$n_matched
          )
        }
      }
    }
  }

  multimodal_tbl <- dplyr::bind_rows(multimodal_rows)
  cluster_tbl <- dplyr::bind_rows(cluster_rows)
  agg_weights_tbl <- dplyr::bind_rows(agg_weight_rows)
  cluster_summary_tbl <- dplyr::bind_rows(cluster_summary_rows)
  functional_summary_tbl <- dplyr::bind_rows(functional_summary_rows)
  anchor_summary_tbl <- dplyr::bind_rows(anchor_summary_rows)
  highest_weight_refit_summary_tbl <- dplyr::bind_rows(highest_weight_refit_summary_rows)
  highest_weight_refit_variants_tbl <- dplyr::bind_rows(highest_weight_refit_variant_rows)
  highest_weight_refit_basin_tbl <- dplyr::bind_rows(highest_weight_refit_basin_rows)
  paper_summary_tbl <- dplyr::bind_rows(paper_summary_rows)
  top_variants_tbl <- dplyr::bind_rows(top_variant_rows)
  comparisons_tbl <- dplyr::bind_rows(comparison_rows)
  aggregated_variants_tbl <- dplyr::bind_rows(aggregated_variant_rows)
  pairwise_tbl <- dplyr::bind_rows(pairwise_rows)

  if (nrow(run_metrics_tbl)) readr::write_csv(run_metrics_tbl, file.path(output_dir, "run_metrics_full.csv"))
  if (nrow(multimodal_tbl)) readr::write_csv(multimodal_tbl, file.path(output_dir, "multimodal_metrics.csv"))
  if (nrow(cluster_tbl)) readr::write_csv(cluster_tbl, file.path(output_dir, "cluster_membership.csv"))
  if (nrow(agg_weights_tbl)) readr::write_csv(agg_weights_tbl, file.path(output_dir, "aggregation_weights_cluster_weight.csv"))
  if (nrow(cluster_summary_tbl)) readr::write_csv(cluster_summary_tbl, file.path(output_dir, "cluster_summary.csv"))
  if (nrow(functional_summary_tbl)) readr::write_csv(functional_summary_tbl, file.path(output_dir, "functional_grid_summary.csv"))
  if (nrow(anchor_summary_tbl)) readr::write_csv(anchor_summary_tbl, file.path(output_dir, "susie_anchor_summary.csv"))
  if (nrow(highest_weight_refit_summary_tbl)) readr::write_csv(highest_weight_refit_summary_tbl, file.path(output_dir, "highest_weight_refit_summary.csv"))
  if (nrow(highest_weight_refit_basin_tbl)) readr::write_csv(highest_weight_refit_basin_tbl, file.path(output_dir, "highest_weight_refit_basin_r2_drift.csv"))
  if (nrow(paper_summary_tbl)) readr::write_csv(paper_summary_tbl, file.path(output_dir, "paper_real_data_ensemble_summary.csv"))
  if (nrow(top_variants_tbl)) readr::write_csv(top_variants_tbl, file.path(output_dir, "top_variants_cluster_weight.csv"))
  if (nrow(comparisons_tbl)) readr::write_csv(comparisons_tbl, file.path(output_dir, "run_comparisons.csv"))
  if (nrow(pairwise_tbl)) real_data_write_parquet(pairwise_tbl, file.path(output_dir, "pairwise_pip_jsd.parquet"))
  if (nrow(aggregated_variants_tbl)) real_data_write_parquet(aggregated_variants_tbl, file.path(output_dir, "aggregated_variant_pips_cluster_weight.parquet"))
  if (nrow(highest_weight_refit_variants_tbl)) real_data_write_parquet(highest_weight_refit_variants_tbl, file.path(output_dir, "highest_weight_refit_variant_posteriors.parquet"))

  metric_inventory <- tibble::tribble(
    ~artifact, ~path, ~format, ~granularity, ~description,
    "run_manifest", file.path(output_dir, "run_manifest.csv"), "csv", "run", "Full job run manifest",
    "run_metrics_full", file.path(output_dir, "run_metrics_full.csv"), "csv", "run", "Per-run real-data fit metrics",
    "effect_summaries_full", file.path(output_dir, "effect_summaries_full.csv"), "csv", "run-effect", "Per-effect alpha and credible-set summaries",
    "dataset_metrics", file.path(output_dir, "dataset_metrics.csv"), "csv", "locus", "Per-locus LD, z, and annotation metrics",
    "validation", file.path(output_dir, "validation.csv"), "csv", "run", "Per-run validation and issues",
    "fit_file_index", file.path(output_dir, "fit_file_index.csv"), "csv", "run", "Index of persisted raw fit RDS files",
    "elbo_trace_full", file.path(output_dir, "elbo_trace_full.csv"), "csv", "run-iteration", "ELBO trace for each run",
    "multimodal_metrics", file.path(output_dir, "multimodal_metrics.csv"), "csv", "locus", "PIP diversity metrics across the functional grid",
    "cluster_membership", file.path(output_dir, "cluster_membership.csv"), "csv", "run", "Per-run JSD cluster assignments",
    "aggregation_weights_cluster_weight", file.path(output_dir, "aggregation_weights_cluster_weight.csv"), "csv", "run", "Cluster-weight aggregation weights per run",
    "cluster_summary", file.path(output_dir, "cluster_summary.csv"), "csv", "locus-cluster", "Representative runs and cluster weights",
    "functional_grid_summary", file.path(output_dir, "functional_grid_summary.csv"), "csv", "run", "Functional grid summary joined to cluster and drift metrics",
    "susie_anchor_summary", file.path(output_dir, "susie_anchor_summary.csv"), "csv", "locus", "susieR anchor vs functional baseline summary",
    "highest_weight_refit_summary", file.path(output_dir, "highest_weight_refit_summary.csv"), "csv", "locus", "Warm baseline refit initialized from the highest-weight SuSiNE ensemble member",
    "highest_weight_refit_basin_r2_drift", file.path(output_dir, "highest_weight_refit_basin_r2_drift.csv"), "csv", "locus", "Basin r2 drift among susie anchor, highest-weight source, and warm refit",
    "paper_real_data_ensemble_summary", file.path(output_dir, "paper_real_data_ensemble_summary.csv"), "csv", "locus", "Paper-facing summary of PIP counts, overlaps, PVE, and JSD",
    "top_variants_cluster_weight", file.path(output_dir, "top_variants_cluster_weight.csv"), "csv", "locus-variant", "Top aggregated variants per locus",
    "run_comparisons", file.path(output_dir, "run_comparisons.csv"), "csv", "run-target", "Per-run comparison metrics against baselines",
    "variant_posteriors_dataset", variant_dataset_dir, "parquet-dataset", "run-variant", "Per-run variant posterior summaries",
    "effect_posteriors_dataset", effect_dataset_dir, "parquet-dataset", "run-effect-variant", "Per-effect posterior summaries",
    "credible_set_membership_dataset", cs_dataset_dir, "parquet-dataset", "run-effect-member", "Credible-set membership rows",
    "pairwise_pip_jsd", file.path(output_dir, "pairwise_pip_jsd.parquet"), "parquet", "run-pair", "Pairwise PIP JSD within each locus",
    "aggregated_variant_pips_cluster_weight", file.path(output_dir, "aggregated_variant_pips_cluster_weight.parquet"), "parquet", "locus-variant", "Cluster-weight aggregated variant PIPs and posterior means",
    "highest_weight_refit_variant_posteriors", file.path(output_dir, "highest_weight_refit_variant_posteriors.parquet"), "parquet", "locus-variant", "Variant posterior summaries for highest-weight warm refits"
  )
  readr::write_csv(metric_inventory, file.path(output_dir, "metric_inventory.csv"))

  list(
    output_dir = output_dir,
    validation = validation_index,
    metric_inventory = metric_inventory
  )
}

validate_real_data_outputs <- function(output_dir) {
  required_csv <- list(
    run_manifest = c("run_id", "task_id", "dataset_bundle_id", "locus_id", "backend", "run_family"),
    run_metrics_full = c("run_id", "locus_id", "backend", "elbo_final", "sigma_2_final", "fit_rds_path"),
    effect_summaries_full = c("run_id", "effect_l", "backend", "top_variant_id", "cs_size_raw"),
    dataset_metrics = c("dataset_bundle_id", "locus_id", "M1", "z_topk_ratio", "annotation_corr_z"),
    validation = c("run_id", "task_id", "has_issues"),
    fit_file_index = c("run_id", "task_id", "fit_rds_path"),
    elbo_trace_full = c("run_id", "iteration", "elbo"),
    multimodal_metrics = c("locus_id", "mean_jsd", "n_clusters"),
    cluster_membership = c("locus_id", "run_id", "cluster_id", "representative_run_id"),
    aggregation_weights_cluster_weight = c("locus_id", "run_id", "agg_weight_run"),
    cluster_summary = c("locus_id", "cluster_id", "representative_run_id", "cluster_weight"),
    functional_grid_summary = c("run_id", "locus_id", "cluster_id", "agg_weight_run"),
    susie_anchor_summary = c("locus_id", "susie_anchor_run_id", "baseline_run_id"),
    top_variants_cluster_weight = c("locus_id", "variant_id", "aggregated_pip"),
    run_comparisons = c("run_id", "target_run_id", "comparison_target", "jsd_pip"),
    metric_inventory = c("artifact", "path", "format", "granularity")
  )

  required_parquet <- list(
    pairwise_pip_jsd = c("locus_id", "run_id_i", "run_id_j", "jsd", "same_cluster"),
    aggregated_variant_pips_cluster_weight = c("locus_id", "variant_id", "aggregated_pip", "aggregated_posterior_mean")
  )

  checks <- list()
  for (nm in names(required_csv)) {
    path <- file.path(output_dir, paste0(nm, ".csv"))
    df <- if (file.exists(path)) readr::read_csv(path, show_col_types = FALSE, n_max = 1L) else NULL
    missing_cols <- if (is.null(df)) required_csv[[nm]] else setdiff(required_csv[[nm]], names(df))
    checks[[length(checks) + 1L]] <- tibble::tibble(
      artifact = nm,
      path = path,
      exists = !is.null(df),
      missing_columns = if (length(missing_cols)) paste(missing_cols, collapse = ", ") else NA_character_,
      ok = !is.null(df) && !length(missing_cols)
    )
  }

  dataset_dirs <- list(
    variant_posteriors_dataset = c("run_id", "locus_id", "variant_id", "pip", "mu_0"),
    effect_posteriors_dataset = c("run_id", "locus_id", "effect_l", "variant_id", "alpha", "b_hat"),
    credible_set_membership_dataset = c("run_id", "locus_id", "effect_l", "set_type", "variant_id")
  )
  for (nm in names(dataset_dirs)) {
    path <- file.path(output_dir, nm)
    df <- if (dir.exists(path)) {
      arrow::open_dataset(path) %>%
        dplyr::slice_head(n = 1L) %>%
        dplyr::collect()
    } else {
      NULL
    }
    missing_cols <- if (is.null(df)) dataset_dirs[[nm]] else setdiff(dataset_dirs[[nm]], names(df))
    checks[[length(checks) + 1L]] <- tibble::tibble(
      artifact = nm,
      path = path,
      exists = !is.null(df),
      missing_columns = if (length(missing_cols)) paste(missing_cols, collapse = ", ") else NA_character_,
      ok = !is.null(df) && !length(missing_cols)
    )
  }

  for (nm in names(required_parquet)) {
    path <- file.path(output_dir, paste0(nm, ".parquet"))
    df <- if (file.exists(path)) arrow::read_parquet(path) %>% tibble::as_tibble() %>% utils::head(1L) else NULL
    missing_cols <- if (is.null(df)) required_parquet[[nm]] else setdiff(required_parquet[[nm]], names(df))
    checks[[length(checks) + 1L]] <- tibble::tibble(
      artifact = nm,
      path = path,
      exists = !is.null(df),
      missing_columns = if (length(missing_cols)) paste(missing_cols, collapse = ", ") else NA_character_,
      ok = !is.null(df) && !length(missing_cols)
    )
  }

  checks_tbl <- dplyr::bind_rows(checks)
  list(
    ok = all(checks_tbl$ok),
    checks = checks_tbl
  )
}
