#!/usr/bin/env Rscript

sanitize_path <- function(path) {
  if (is.null(path)) return(path)
  gsub("~+~", " ", path, fixed = TRUE)
}

suppressPackageStartupMessages(library(devtools))

reload_local_packages <- function(susine_path, repo_root) {
  if (!requireNamespace("pkgload", quietly = TRUE)) {
    stop("pkgload is required to reload local packages cleanly.")
  }

  if ("package:test_susine" %in% search() || "test_susine" %in% loadedNamespaces()) {
    try(pkgload::unload("test_susine"), silent = TRUE)
  }
  if ("package:susine" %in% search() || "susine" %in% loadedNamespaces()) {
    try(pkgload::unload("susine"), silent = TRUE)
  }

  pkgload::load_all(susine_path, quiet = TRUE, reset = TRUE)
  pkgload::load_all(repo_root, quiet = TRUE, reset = TRUE)
}

find_repo_root <- function() {
  cmd <- commandArgs(trailingOnly = FALSE)
  file_arg <- "--file="
  file_idx <- grep(file_arg, cmd)
  script_dir <- NULL
  if (length(file_idx)) {
    script_path <- sub(file_arg, "", cmd[file_idx[1]])
    script_path <- sanitize_path(script_path)
    script_dir <- dirname(normalizePath(script_path, winslash = "/", mustWork = TRUE))
  } else if (!is.null(sys.frames()[[1]]$ofile)) {
    script_dir <- dirname(normalizePath(sanitize_path(sys.frames()[[1]]$ofile), winslash = "/", mustWork = TRUE))
  }
  if (is.null(script_dir)) {
    script_dir <- getwd()
  }
  normalizePath(file.path(script_dir, "..", ".."), winslash = "/", mustWork = TRUE)
}

repo_root <- find_repo_root()
setwd(repo_root)

use_dev <- Sys.getenv("SUSINE_DEV", unset = "") == "1"
if (use_dev) {
  susine_path <- file.path(dirname(repo_root), "susine")
  if (dir.exists(susine_path) && file.exists(file.path(susine_path, "DESCRIPTION"))) {
    message("Loading local dev susine from: ", susine_path)
    options(test_susine.local_susine_path = normalizePath(susine_path, winslash = "/", mustWork = TRUE))
    reload_local_packages(susine_path, repo_root)
    susine_fn <- getExportedValue("susine", "susine_rss")
    required_susine_rss_args <- c("prior_update_method", "estimate_residual_variance")
    missing_susine_rss_args <- setdiff(required_susine_rss_args, names(formals(susine_fn)))
    if (length(missing_susine_rss_args)) {
      stop(
        "Local dev susine was loaded from ", susine_path,
        " but susine::susine_rss() lacks expected argument(s): ",
        paste(missing_susine_rss_args, collapse = ", "), ". ",
        "That sibling checkout is stale or not the expected branch."
      )
    }
  } else {
    pkgload::load_all(repo_root, quiet = TRUE, reset = TRUE)
  }
} else {
  suppressPackageStartupMessages(library(test_susine))
}

parse_args <- function(args) {
  out <- list()
  i <- 1L
  while (i <= length(args)) {
    arg <- args[[i]]
    if (startsWith(arg, "--")) {
      key_val <- sub("^--", "", arg)
      if (grepl("=", key_val, fixed = TRUE)) {
        kv <- strsplit(key_val, "=", fixed = TRUE)[[1]]
        out[[kv[1]]] <- kv[2]
      } else {
        if (i == length(args)) stop("Missing value for argument ", arg)
        out[[key_val]] <- args[[i + 1L]]
        i <- i + 1L
      }
    }
    i <- i + 1L
  }
  out
}

main <- function() {
  args <- parse_args(commandArgs(trailingOnly = TRUE))
  job_name <- args[["job-name"]]
  task_id <- args[["task-id"]]
  job_root <- args[["job-root"]]
  config_path <- args[["config-path"]]
  quiet_flag <- args[["quiet"]]
  if (is.null(quiet_flag)) quiet_flag <- "FALSE"

  if (is.null(job_name) || is.null(task_id) || is.null(config_path)) {
    stop("Arguments --job-name, --task-id, and --config-path are required.")
  }

  task_id <- as.integer(task_id)
  if (is.na(task_id)) {
    stop("task-id must be an integer.")
  }

  quiet <- tolower(quiet_flag) %in% c("1", "true", "yes")

  tryCatch(
    {
      test_susine::run_real_data_task(
        job_name = job_name,
        task_id = task_id,
        job_root = if (is.null(job_root)) "output" else job_root,
        config_path = config_path,
        quiet = quiet
      )
      invisible(NULL)
    },
    error = function(e) {
      message("Error during task execution: ", conditionMessage(e))
      quit(save = "no", status = 1L)
    }
  )
}

main()
