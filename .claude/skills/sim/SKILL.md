---
name: sim
description: Enter simulation workflow mode. Use when configuring, running, debugging, or analyzing simulation experiments.
user-invocable: true
---

# Simulation Workflow Mode

You are helping configure, run, debug, or analyze simulation experiments for the SuSiNE benchmarking study.

## Research Context (load on demand for harness design decisions)

@refs/susine_study_plan_high_level.md
@refs/susie_susine_background.md

Read these when the user asks about *why* the harness is structured a certain way,
when adding new scenarios or use cases, or when making architectural changes to the
simulation pipeline. They contain the study hypotheses (H1-H6), failure-mode
characterization framework, and the phased study structure that drives what the
harness needs to support.

## Pipeline Overview

```
1. Configure  -->  make_job_config() or run_control_workbook.Rmd
2. Write      -->  write_job_artifacts() produces job_config.json, run_table.csv, .slurm
3. Submit     -->  sbatch output/slurm_scripts/<job>.slurm
4. Execute    -->  Each SLURM array task calls inst/scripts/run_task.R -> run_task()
5. Collect    -->  aggregate_staging_outputs() combines flush files
6. Visualize  -->  visualize_results_workbook_paper_exhibits.Rmd
```

## Key Configuration Parameters (make_job_config)

| Parameter | Description | Example values |
|-----------|-------------|----------------|
| `use_case_ids` | Model configs from use_case_catalog() | `c("a_i", "b_ii", "c_ii")` |
| `L_grid` | Number of single effects | `5`, `10`, `c(5, 10)` |
| `p_star_grid` | True causal SNP count | `c(1, 2, 3, 4, 5, 10, 20)` |
| `y_noise_grid` | Noise fraction (0=no noise) | `c(0.05, 0.1, 0.2, 0.4)` |
| `prior_quality` | Annotation settings | `prior_quality_grid(...)` |
| `seeds` | RNG seeds for phenotype sim | `1:3` |
| `data_scenarios` | Genotype source | `"simulation_n3"` |
| `grid_mode` | How grids combine | `"full"`, `"minimal"`, `"intersect"` |
| `sigma_0_2_scalars` | Prior variance multipliers | `c("0.1", "1/L")` |
| `annotation_scales` | c-grid for functional mu | `c(0, 0.5, 1, 2, 5)` |
| `aggregation_methods` | BMA strategies | `c("softmax_elbo", "mean", "max_elbo")` |
| `restart_settings` | Random init config | `list(n_inits = 5, alpha_concentration = 1)` |

## Use Case Groups (from use_case_catalog())

- **Group 1a** (SuSiE variants): a_i (baseline), a_i_restart (with restarts), a_ii (EB sigma), a_iii (EB mu), a_iv (EB both)
- **Group 1b** (functional priors): b_i (functional sigma, mu=0), b_ii (functional mu), b_iii (functional mu + EB sigma), b_iv-b_vi (auto-scaled variants)
- **Group 1c** (extra compute): c_i (annealing), c_ii (model averaging)

## Output Directory Structure

```
output/
  temp/<job>/                       # Temporary config (overwritten on re-create)
    job_config.json
    run_table.csv
    dataset_bundles.csv
    use_cases.csv
  run_history/<job>/<parent>/       # Immutable config snapshot
  slurm_output/<job>/<parent>/      # Per-task flush files
    task-001/
      flush-000_model_metrics.csv
      flush-000_confusion_bins.csv
      flush-000_snps.parquet
    aggregated/                     # Post-collection output
      model_metrics.csv
      confusion_bins.csv
      snps_dataset/                 # Partitioned by use_case_id
  slurm_scripts/<job>.slurm
  slurm_prints/<job>/<parent>/      # SLURM stdout/stderr
```

## Quick Local Test

```r
devtools::load_all("../susine")
devtools::load_all(".")

cfg <- make_job_config(
  job_name = "debug_run",
  use_case_ids = c("a_i", "b_ii"),
  L_grid = 5, y_noise_grid = 0.2,
  prior_quality = prior_quality_grid(),
  p_star_grid = 5, seeds = 1:2
)
write_job_artifacts(cfg)
run_task("debug_run", task_id = 1,
         config_path = "output/temp/debug_run/job_config.json")
```

## Common Debugging

- **"Task not found"**: Check task_id range vs nrow(tasks) in the config JSON.
- **"Matrix path not found"**: Path normalization issue (Windows vs Linux). Check repo_root in job_config.json.
- **SLURM OOM**: Increase `mem` in make_job_config() or set `verbose_file_output = FALSE` to reduce memory from buffering.
- **Local testing**: Set SUSINE_DEV=1 env var, call run_task() directly in R.
- **Stale package**: Always re-run devtools::load_all() for both susine and test_susine after code changes.
- **Flush file missing columns**: Check that evaluation_metrics.R and run_model.R agree on output schema after edits.