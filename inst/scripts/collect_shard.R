#!/usr/bin/env Rscript

library(devtools)
devtools::load_all(".")

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

truthy <- function(x) {
  tolower(x) %in% c("1", "true", "yes")
}

main <- function() {
  args <- parse_args(commandArgs(trailingOnly = TRUE))
  job_name <- args[["job-name"]]
  parent_job_id <- args[["parent-id"]]
  if (is.null(parent_job_id)) {
    parent_job_id <- args[["parent_job_id"]]
  }
  shard_index <- args[["shard-index"]]
  output_root <- args[["output-root"]]
  if (is.null(output_root)) {
    output_root <- "output"
  }
  shard_size <- args[["shard-size"]]
  force_flag <- args[["force"]]
  if (is.null(force_flag)) {
    force_flag <- "FALSE"
  }
  quiet_flag <- args[["quiet"]]
  if (is.null(quiet_flag)) {
    quiet_flag <- "FALSE"
  }

  if (is.null(job_name) || is.null(parent_job_id) || is.null(shard_index)) {
    stop("Arguments --job-name, --parent-id, and --shard-index are required.")
  }

  shard_index <- as.integer(shard_index)
  if (is.na(shard_index)) {
    stop("shard-index must be an integer.")
  }

  if (!is.null(shard_size)) {
    shard_size <- as.integer(shard_size)
  }

  force <- truthy(force_flag)
  quiet <- truthy(quiet_flag)

  tryCatch(
    {
      test_susine::collect_shard_metrics(
        job_name = job_name,
        parent_job_id = parent_job_id,
        shard_index = shard_index,
        output_root = output_root,
        shard_size = shard_size,
        force = force,
        quiet = quiet
      )
      invisible(NULL)
    },
    error = function(e) {
      message("Error during shard collection: ", conditionMessage(e))
      quit(save = "no", status = 1L)
    }
  )
}

main()
