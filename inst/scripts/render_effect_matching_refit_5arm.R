#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)

arg_value <- function(name, default = NULL) {
  prefix <- paste0("--", name, "=")
  hit <- args[startsWith(args, prefix)]
  if (length(hit)) {
    return(sub(prefix, "", hit[[length(hit)]]))
  }
  flag_idx <- which(args == paste0("--", name))
  if (length(flag_idx) && flag_idx[[length(flag_idx)]] < length(args)) {
    return(args[[flag_idx[[length(flag_idx)]] + 1L]])
  }
  default
}

arg_flag <- function(name, default = FALSE) {
  val <- arg_value(name, NA_character_)
  if (is.na(val)) {
    return(isTRUE(default))
  }
  tolower(val) %in% c("1", "true", "t", "yes", "y")
}

repo_root <- normalizePath(arg_value("repo-root", getwd()), winslash = "/", mustWork = TRUE)
fresh <- arg_flag("fresh", TRUE)
max_datasets <- arg_value("max-datasets", Sys.getenv("EFFECT_MATCHING_MAX_DATASETS", ""))
checkpoint_every <- arg_value("checkpoint-every", Sys.getenv("EFFECT_MATCHING_CHECKPOINT_EVERY", "25"))
shard_total <- as.integer(arg_value("shard-total", Sys.getenv("EFFECT_MATCHING_SHARD_TOTAL", "1")))
shard_index <- as.integer(arg_value("shard-index", Sys.getenv("EFFECT_MATCHING_SHARD_INDEX", "1")))
if (!is.finite(shard_total) || shard_total < 1L) shard_total <- 1L
if (!is.finite(shard_index) || shard_index < 1L || shard_index > shard_total) {
  stop("Invalid shard settings: shard-index=", shard_index,
       ", shard-total=", shard_total)
}
shard_suffix <- if (shard_total > 1L) {
  sprintf("_shard_%03d_of_%03d", shard_index, shard_total)
} else {
  ""
}

setwd(repo_root)
options(bitmapType = "cairo")

study_dir <- file.path(
  repo_root,
  "output",
  "effect_matching_study",
  "baseline_sims_5arm_effect_matching"
)
dir.create(study_dir, recursive = TRUE, showWarnings = FALSE)

if (fresh && shard_total <= 1L) {
  stamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
  backup_dir <- file.path(study_dir, paste0("backup_before_slurm_refit_", stamp))
  stale_files <- file.path(
    study_dir,
    c(
      "effect_refit_records.rds",
      "effect_refit_metadata.csv",
      "pip_records_5arm.rds",
      "refit_progress_log.csv"
    )
  )
  existing <- stale_files[file.exists(stale_files)]
  if (length(existing)) {
    dir.create(backup_dir, recursive = TRUE, showWarnings = FALSE)
    file.copy(existing, backup_dir, overwrite = TRUE)
    unlink(existing)
    message("Backed up stale 5-arm outputs to: ", backup_dir)
  } else {
    message("No stale 5-arm outputs found.")
  }
} else if (fresh && shard_total > 1L) {
  shard_files <- file.path(
    study_dir,
    paste0(
      c("effect_refit_records", "effect_refit_metadata",
        "refit_progress_log"),
      shard_suffix,
      c(".rds", ".csv", ".csv")
    )
  )
  shard_files <- c(
    shard_files,
    file.path(study_dir, paste0("pip_records_5arm", shard_suffix, ".rds"))
  )
  existing <- shard_files[file.exists(shard_files)]
  if (length(existing)) {
    unlink(existing)
    message("Removed stale shard outputs: ", paste(basename(existing), collapse = ", "))
  } else {
    message("No stale shard outputs found for ", shard_suffix, ".")
  }
}

Sys.setenv(
  TEST_SUSINE_REPO_ROOT = repo_root,
  EFFECT_MATCHING_FORCE_REFIT = "false",
  EFFECT_MATCHING_CHECKPOINT_EVERY = checkpoint_every,
  EFFECT_MATCHING_SHARD_TOTAL = as.character(shard_total),
  EFFECT_MATCHING_SHARD_INDEX = as.character(shard_index)
)
if (nzchar(max_datasets)) {
  Sys.setenv(EFFECT_MATCHING_MAX_DATASETS = max_datasets)
}

message("Repo root: ", repo_root)
message("Study dir: ", study_dir)
message("Fresh run: ", fresh)
message("Checkpoint every: ", checkpoint_every)
message("Max datasets: ", Sys.getenv("EFFECT_MATCHING_MAX_DATASETS", unset = "Inf"))
message("Shard: ", shard_index, " of ", shard_total)

if (!requireNamespace("rmarkdown", quietly = TRUE)) {
  stop("The rmarkdown package is required to render the 5-arm export workbook.")
}

rmarkdown::render(
  input = file.path("vignettes", "one_off_validations", "effect_matching_refit_export_5arm.Rmd"),
  output_format = "html_document",
  quiet = FALSE,
  envir = new.env(parent = globalenv())
)
