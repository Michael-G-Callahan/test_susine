options(warn = 1)
cat('R:', R.version.string, '\n')
suppressPackageStartupMessages({
  library(devtools)
  library(readr)
})

devtools::load_all('../susine', quiet = TRUE)
devtools::load_all('.', quiet = TRUE)

mk_cfg <- function(job_name, methods, mode = 'separate', K = 1L, restart_settings = NULL, p_star = 2L) {
  pq <- prior_quality_grid(c(0.3), c(1))
  make_job_config(
    job_name = job_name,
    use_case_ids = c('susie_vanilla'),
    exploration_methods = methods,
    exploration_mode = mode,
    K = as.integer(K),
    L_grid = 5L,
    y_noise_grid = 0.8,
    prior_quality = pq,
    p_star_grid = as.integer(p_star),
    seeds = 1L,
    data_scenarios = 'simulation_n3',
    repo_root = '.',
    task_unit = 'dataset',
    bundles_per_task = 1L,
    sigma_0_2_scalars = '0.2',
    restart_settings = restart_settings
  )
}

# 1) single
cfg_single <- mk_cfg('tmp_smoke_codex', methods = c('single'), mode = 'separate', K = 1L)
write_job_artifacts(cfg_single, 'inst/scripts/run_task.R')
run_task('tmp_smoke_codex', 1L, job_root = 'output', config_path = 'output/temp/tmp_smoke_codex/job_config.json', quiet = TRUE)
cat('OK single run\n')

# 2) restart K=4
cfg_restart <- mk_cfg('tmp_restart_codex', methods = c('restart'), mode = 'separate', K = 4L,
  restart_settings = list(n_inits = 4L, alpha_concentration = 0.1)
)
print(table(cfg_restart[['tables']][['runs']][['init_type']], useNA = 'ifany'))
write_job_artifacts(cfg_restart, 'inst/scripts/run_task.R')
run_task('tmp_restart_codex', 1L, job_root = 'output', config_path = 'output/temp/tmp_restart_codex/job_config.json', quiet = TRUE)
cat('OK restart run\n')

# 3) intersect valid
cfg_intersect <- mk_cfg('tmp_intersect_codex', methods = c('single', 'restart'), mode = 'intersect', K = 2L,
  restart_settings = list(n_inits = 2L, alpha_concentration = 0.1)
)
cat('intersect rows:', nrow(cfg_intersect[['tables']][['runs']]), '\n')
write_job_artifacts(cfg_intersect, 'inst/scripts/run_task.R')
run_task('tmp_intersect_codex', 1L, job_root = 'output', config_path = 'output/temp/tmp_intersect_codex/job_config.json', quiet = TRUE)
cat('OK intersect run\n')

# 4) aggregation path + index/aggregate
job <- 'tmp_agg_codex'
parent <- 'local_test'
cfg_agg <- mk_cfg(job, methods = c('single'), mode = 'separate', K = 1L)
write_job_artifacts(cfg_agg, 'inst/scripts/run_task.R')
Sys.setenv(SUSINE_JOB_NAME = job, SUSINE_PARENT_ID = parent)
run_task(job, 1L, job_root = 'output', config_path = file.path('output','temp',job,'job_config.json'), quiet = TRUE)
base_dir <- file.path('output','slurm_output',job,parent)
cat('base_dir exists:', dir.exists(base_dir), 'path:', base_dir, '\n')
if (dir.exists(base_dir)) {
  td <- list.dirs(base_dir, recursive = FALSE, full.names = TRUE)
  cat('task dirs:', length(td), '\n')
  if (length(td)) {
    ff <- list.files(td[[1]], full.names = TRUE)
    cat('files in first task dir:', length(ff), '\n')
    print(basename(ff))
  }
}
idx <- index_staging_outputs(job, parent, output_root = 'output')
cat('indexed files:', nrow(idx), '\n')
if (nrow(idx)) print(table(idx[['type']]))
agg <- aggregate_staging_outputs(job, parent, output_root = 'output', validate = TRUE)
cat('agg output dir:', agg[['output_dir']], '\n')
out_files <- if (dir.exists(agg[['output_dir']])) list.files(agg[['output_dir']], recursive = TRUE) else character(0)
cat('aggregated file count:', length(out_files), '\n')
print(out_files)

cat('ALL SMOKE TESTS COMPLETED\n')
