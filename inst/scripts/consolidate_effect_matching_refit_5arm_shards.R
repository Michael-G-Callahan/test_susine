#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)

arg_value <- function(name, default = NULL) {
  prefix <- paste0("--", name, "=")
  hit <- args[startsWith(args, prefix)]
  if (length(hit)) return(sub(prefix, "", hit[[length(hit)]]))
  flag_idx <- which(args == paste0("--", name))
  if (length(flag_idx) && flag_idx[[length(flag_idx)]] < length(args)) {
    return(args[[flag_idx[[length(flag_idx)]] + 1L]])
  }
  default
}

repo_root <- normalizePath(arg_value("repo-root", getwd()), winslash = "/", mustWork = TRUE)
shard_total <- as.integer(arg_value("shard-total", Sys.getenv("EFFECT_MATCHING_SHARD_TOTAL", "20")))
if (!is.finite(shard_total) || shard_total < 1L) stop("Invalid shard-total: ", shard_total)

study_dir <- file.path(
  repo_root, "output", "effect_matching_study", "baseline_sims_5arm_effect_matching"
)
if (!dir.exists(study_dir)) stop("Study dir does not exist: ", study_dir)

suffixes <- sprintf("_shard_%03d_of_%03d", seq_len(shard_total), shard_total)
records_files <- file.path(study_dir, paste0("effect_refit_records", suffixes, ".rds"))
metadata_files <- file.path(study_dir, paste0("effect_refit_metadata", suffixes, ".csv"))
progress_files <- file.path(study_dir, paste0("refit_progress_log", suffixes, ".csv"))
pip_files <- file.path(study_dir, paste0("pip_records_5arm", suffixes, ".rds"))

missing <- c(records_files, metadata_files, pip_files)[!file.exists(c(records_files, metadata_files, pip_files))]
if (length(missing)) {
  stop("Missing shard outputs:\n", paste(missing, collapse = "\n"))
}

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
})

message("Consolidating ", shard_total, " shards from: ", study_dir)

records <- bind_rows(lapply(records_files, readRDS))
saveRDS(records, file.path(study_dir, "effect_refit_records.rds"))

metadata <- bind_rows(lapply(metadata_files, readr::read_csv, show_col_types = FALSE))
readr::write_csv(metadata, file.path(study_dir, "effect_refit_metadata.csv"))

pip_records <- bind_rows(lapply(pip_files, readRDS))
saveRDS(pip_records, file.path(study_dir, "pip_records_5arm.rds"))

existing_progress <- progress_files[file.exists(progress_files)]
if (length(existing_progress)) {
  progress <- bind_rows(lapply(existing_progress, readr::read_csv, show_col_types = FALSE))
  readr::write_csv(progress, file.path(study_dir, "refit_progress_log.csv"))
}

record_check <- records %>%
  count(.data$dataset_bundle_id, .data$model_method, name = "n_effects") %>%
  count(.data$dataset_bundle_id, name = "n_methods")

bad <- record_check %>% filter(.data$n_methods != 5L)
if (nrow(bad)) {
  stop("Consolidated records have datasets without all 5 methods. Bad datasets: ",
       paste(head(bad$dataset_bundle_id, 20), collapse = ", "))
}

message("Wrote consolidated effect_refit_records.rds rows: ", nrow(records))
message("Wrote consolidated effect_refit_metadata.csv rows: ", nrow(metadata))
message("Wrote consolidated pip_records_5arm.rds rows: ", nrow(pip_records))
message("Datasets: ", n_distinct(records$dataset_bundle_id))
