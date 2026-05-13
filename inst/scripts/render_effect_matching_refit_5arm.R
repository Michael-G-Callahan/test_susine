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

setwd(repo_root)

study_dir <- file.path(
  repo_root,
  "output",
  "effect_matching_study",
  "baseline_sims_5arm_effect_matching"
)
dir.create(study_dir, recursive = TRUE, showWarnings = FALSE)

if (fresh) {
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
}

Sys.setenv(
  EFFECT_MATCHING_FORCE_REFIT = "false",
  EFFECT_MATCHING_CHECKPOINT_EVERY = checkpoint_every
)
if (nzchar(max_datasets)) {
  Sys.setenv(EFFECT_MATCHING_MAX_DATASETS = max_datasets)
}

message("Repo root: ", repo_root)
message("Study dir: ", study_dir)
message("Fresh run: ", fresh)
message("Checkpoint every: ", checkpoint_every)
message("Max datasets: ", Sys.getenv("EFFECT_MATCHING_MAX_DATASETS", unset = "Inf"))

if (!requireNamespace("rmarkdown", quietly = TRUE)) {
  stop("The rmarkdown package is required to render the 5-arm export workbook.")
}

rmarkdown::render(
  input = file.path("vignettes", "one_off_validations", "effect_matching_refit_export_5arm.Rmd"),
  output_format = "html_document",
  quiet = FALSE,
  envir = new.env(parent = globalenv())
)
