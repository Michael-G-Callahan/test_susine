suppressPackageStartupMessages(library(devtools))
devtools::load_all('../susine', quiet = TRUE)
devtools::load_all('.', quiet = TRUE)
job <- 'tmp_agg_codex2'
parent <- 'local_test'
pq <- prior_quality_grid(c(0.3), c(1))
cfg <- make_job_config(
  job_name = job,
  use_case_ids = c('susie_vanilla'),
  exploration_methods = c('single'),
  exploration_mode = 'separate',
  K = 1L,
  L_grid = 5L,
  y_noise_grid = 0.8,
  prior_quality = pq,
  p_star_grid = 2L,
  seeds = 1L,
  data_scenarios = 'simulation_n3',
  repo_root = '.',
  task_unit = 'dataset',
  bundles_per_task = 1L,
  sigma_0_2_scalars = '0.2'
)
write_job_artifacts(cfg, 'inst/scripts/run_task.R')
Sys.setenv(SUSINE_JOB_NAME = job, SUSINE_PARENT_ID = parent)
run_task(job, 1L, job_root = 'output', config_path = file.path('output','temp',job,'job_config.json'), quiet = TRUE)
idx <- index_staging_outputs(job, parent, output_root = 'output')
cat('indexed:', nrow(idx), '\n')
if (nrow(idx)) print(table(idx[['type']]))
agg <- aggregate_staging_outputs(job, parent, output_root = 'output', validate = TRUE)
cat('outdir:', agg[['output_dir']], '\n')
print(list.files(agg[['output_dir']], recursive = TRUE))
