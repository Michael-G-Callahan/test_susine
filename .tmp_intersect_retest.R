suppressPackageStartupMessages(library(devtools))
devtools::load_all('../susine', quiet = TRUE)
devtools::load_all('.', quiet = TRUE)
pq <- prior_quality_grid(c(0.3), c(1))
cfg <- make_job_config(
  job_name = 'tmp_intersect_codex_fixed',
  use_case_ids = c('susie_vanilla'),
  exploration_methods = c('single','restart'),
  exploration_mode = 'intersect',
  K = 2L,
  L_grid = 5L,
  y_noise_grid = 0.8,
  prior_quality = pq,
  p_star_grid = 2L,
  seeds = 1L,
  data_scenarios = 'simulation_n3',
  repo_root = '.',
  task_unit = 'dataset',
  bundles_per_task = 1L,
  sigma_0_2_scalars = '0.2',
  restart_settings = list(n_inits = 2L, alpha_concentration = 0.1)
)
cat('rows:', nrow(cfg[['tables']][['runs']]), '\n')
print(vapply(cfg[['tables']][['runs']], function(x) is.list(x), logical(1)))
write_job_artifacts(cfg, 'inst/scripts/run_task.R')
run_task('tmp_intersect_codex_fixed', 1L, job_root = 'output', config_path = 'output/temp/tmp_intersect_codex_fixed/job_config.json', quiet = TRUE)
cat('OK intersect run\n')
