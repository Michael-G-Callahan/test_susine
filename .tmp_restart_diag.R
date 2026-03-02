suppressPackageStartupMessages(library(devtools))
devtools::load_all('../susine', quiet = TRUE)
devtools::load_all('.', quiet = TRUE)
pq <- prior_quality_grid(c(0.3), c(1))
cfg <- make_job_config(
  job_name = 'tmp_restart_diag',
  use_case_ids = c('susie_vanilla'),
  exploration_methods = c('restart'),
  exploration_mode = 'separate',
  K = 4L,
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
  restart_settings = list(n_inits = 4L, alpha_concentration = 0.1)
)
r <- cfg[['tables']][['runs']]
print(vapply(r, function(x) paste(class(x), collapse='|'), character(1)))
cat('list cols:', paste(names(r)[vapply(r, is.list, logical(1))], collapse=', '), '\n')
