# Ensemble Simulation Study Technical Trace

Date written: 2026-05-08
Last updated: 2026-06-15

This document traces the ensemble-scaling simulation study as implemented in the
`test_susine` codebase. It is intentionally technical and implementation-facing:
the goal is to preserve the details needed to later write a cleaner methods
section.

I did not include Slurm resource settings or cluster submission mechanics except
where they affect the scientific workflow or the shape of the outputs.

## 1. Scope, provenance, and key caveats

The main files traced were:

- Manuscript context:
  - `../Writings/second draft/susine_second_draft.tex`
- Run-control workbook:
  - `test_susine/vignettes/simulation pipeline/run_control_workbook_ensemble_scaling.Rmd`
- Collection workbook:
  - `test_susine/vignettes/simulation pipeline/collect_results_workbook_ensemble_scaling.Rmd`
- Broad exploratory visualization workbook:
  - `test_susine/vignettes/simulation pipeline/visualize_results_workbook_ensemble_scaling.Rmd`
- Paper-plot preparation workbook:
  - `test_susine/vignettes/simulation pipeline/prepare_results_workbook_ensemble_scaling_paper.Rmd`
- Paper-plot rendering workbook:
  - `test_susine/vignettes/simulation pipeline/visualize_results_workbook_ensemble_scaling_paper.Rmd`
- Core R helpers:
  - `test_susine/R/use_cases.R`
  - `test_susine/R/run_controls.R`
  - `test_susine/R/simulate_data.R`
  - `test_susine/R/run_model.R`
  - `test_susine/R/evaluation_metrics.R`
  - `test_susine/R/collect_results.R`
  - `test_susine/R/ensemble_scaling_utils.R`
  - `test_susine/R/visualize_results.R`

Important provenance notes:

- The manuscript draft has placeholders for the ensemble simulation methods and
  results. In `susine_second_draft.tex`, the ensemble methods subsection is
  currently just `\subsection{Simulation study --- ensemble}` and the ensemble
  results subsection is also empty. The baseline simulation methods/results and
  the introductory framing provide the relevant context.
- The collection workbook literal differs from the paper plot workbooks:
  `collect_results_workbook_ensemble_scaling.Rmd` has
  `parent_job_id <- "53520405"`, while both paper-prep and paper-rendering
  workbooks use `parent_job_id <- "53547760"`. The paper figures therefore point
  to consolidated outputs for job `ensemble_scaling_full/53547760`.
- In this local checkout, the full `ensemble_scaling_full/53547760` output
  directory may not be present under `test_susine/output/slurm_output` or
  `output/slurm_output`. The trace below is therefore based on the source code
  and workbook literals rather than re-reading the final generated CSVs.
- The run-control workbook scheduled annotation quality as
  `annotation_r2_levels <- c(0, 0.3, 0.5)` with `inflate_match = 0.95`. This maps
  onto the simulator's two knobs as `annotation_r2` and `inflate_match`, not the
  full manuscript baseline grid of `(phi_a, nu_a)`. The paper-prep workflow
  focuses on `annotation_r2_focus <- 0.3`.
- The final per-variant classification curves in this ensemble workflow use the
  top 8 causal variants by absolute effect size as positives. The remaining
  (weak) causal variants are not counted as positives -- and, critically, they
  are not counted as negatives either: `compute_confusion_bins()` DROPS them
  entirely (`keep <- !(causal == 1L & scored == 0L)`), so the confusion-bin
  scoring used in pooled AUPRC and precision-recall plots is neither rewarded
  nor penalized for them.
- The scheduled cluster-weight methods are now the pilot-locked
  `cluster_weight_credible` (clustering metric `credible_shift` =
  max_j[ max(p_j,q_j) * |p_j - q_j| ], complete-linkage cut at 0.05) and
  `cluster_weight_max` (clustering metric `max_shift` = max_j |p_j - q_j|,
  L-infinity, complete-linkage cut at 0.10). Both aggregate with **Method B**
  (frequency-free): each cluster gets weight proportional to
  `exp(max-ELBO-in-cluster)`, split within the cluster across all member fits by
  ELBO-softmax. All fits contribute -- there is no single nominee and no
  1/frequency correction. Ground truth: `run_model.R` ~line 2005
  (`credible_shift`), ~2009 (`max_shift`), ~2086-2126
  (`.cluster_weights_from_hc(..., aggregation = "method_b")`), ~2143-2146
  (`.cluster_weight_specs`). For these scheduled methods both the full-ensemble
  path and the scaling-bin path now route through Method B identically, so the
  older full-vs-scaling discrepancy is obsolete for them. The legacy JSD-0.15 +
  frequency `cluster_weight` method (Method C) still exists for backward
  compatibility in lower-level helpers, but it is not the current paper
  ensemble-simulation aggregator and is not the current real-data aggregator.
  The paper-prep workbook highlights `cluster_weight_credible` as the primary
  method.

## 2. Manuscript context

The manuscript currently frames SuSiNE as opening a signed prior-mean channel,
`mu_0 = c * a`, where `a` is a directional functional annotation vector and `c`
is a scale. The relevant conceptual point for the ensemble study is that the
manuscript treats two failure modes as distinct but related:

- an optimizer barrier, where IBSS can converge to a bad local basin;
- a model-specification barrier, where one variational SuSiE/SuSiNE fit cannot
  represent OR-of-ANDs posterior uncertainty across distinct multi-variant
  configurations.

The manuscript introduction says the proposed ensembling workflow fits multiple
SuSiNE models under perturbed initializations and priors, then aggregates PIPs
across the ensemble using ELBO/JSD-informed weights. That is exactly what the
ensemble-scaling harness operationalizes.

The baseline simulation section establishes the simulation substrate:

- 150 genotype matrices;
- 4 phenotype seeds per matrix, for 600 simulated datasets;
- `n = 600` individuals;
- approximately `p = 1000` SNPs per locus;
- SuSiE 2.0 oligogenic architecture:
  - 3 strong effects;
  - 5 moderate/oligogenic effects;
  - 15 polygenic-background effects;
  - 23 total causal variants per locus;
  - target total heritability `h2_total = 0.25`;
- `L = 10` fitted single-effect components;
- cap of 100 IBSS iterations;
- credible-set cumulative posterior threshold `rho = 0.95`.

The baseline results section also explains why the ensemble study carries
forward only the `mu_0`-based SuSiNE variants:

- annotation-agnostic SuSiE-like methods were tightly clustered;
- the signed `mu_0` channel beat the unsigned functional-pi channel under matched
  signed annotations;
- the manuscript explicitly says SuSiNE functional-mu, both grid and EB, is
  carried forward, while SuSiE-ash, SuSiE-inf, and functional-pi are discontinued
  for the ensemble study.

The ensemble run-control code implements that decision:

- Group A: annotation-agnostic, zero-prior-mean `susine_vanilla`;
- Group B: EB SuSiNE with signed annotations;
- Group C: fixed-variance grid SuSiNE with signed annotations.

## 3. End-to-end workflow

The study has four main computational stages.

1. Build the heterogeneous ensemble job configuration.
   - Workbook:
     `test_susine/vignettes/simulation pipeline/run_control_workbook_ensemble_scaling.Rmd`
   - Key helper:
     `test_susine/R/run_controls.R:750` (`make_job_config`)
   - Output tables include dataset bundles, matrix catalog rows, use-case
     metadata, run rows, and task rows.

2. Execute all run rows by dataset bundle.
   - Key helper:
     `test_susine/R/run_model.R:13` (`run_task`)
   - For each dataset bundle:
     - regenerate the simulated phenotype;
     - generate annotations on demand for annotation-consuming methods;
     - fit all individual ensemble members;
     - write individual metrics and confusion bins;
     - aggregate PIPs within ensemble groups;
     - write aggregated confusion bins;
     - write scaling-subset confusion bins.

3. Collect staged task outputs into consolidated CSVs.
   - Workbook:
     `test_susine/vignettes/simulation pipeline/collect_results_workbook_ensemble_scaling.Rmd`
   - Key helper modules:
     `test_susine/R/collect_results.R`,
     `test_susine/R/ensemble_scaling_utils.R`,
     `test_susine/R/visualize_results.R`
   - This stage pools confusion bins to compute AUPRC and TPR@FPR=0.05,
     constructs oracle rows, backfills missing aggregation rows where possible,
     and exports the consolidated analysis tables.

4. Prepare and render paper plots.
   - Prep workbook:
     `test_susine/vignettes/simulation pipeline/prepare_results_workbook_ensemble_scaling_paper.Rmd`
   - Rendering workbook:
     `test_susine/vignettes/simulation pipeline/visualize_results_workbook_ensemble_scaling_paper.Rmd`
   - The prep workbook filters, reshapes, bootstraps, and caches plot data in an
     RDS file.
   - The rendering workbook reads that RDS and saves final PNGs.

## 4. Run-control settings

The ensemble run-control workbook sets:

- `job_name <- "ensemble_scaling_full"`
- `n_matrices_target <- 150L`
- `seeds <- 1:4`
- `L_grid <- 10L`
- `architecture_grid <- "susie2_oligogenic"`
- `annotation_r2_levels <- c(0, 0.3, 0.5)`
- `noncausal_inflation_factor <- 0.95` for all annotation levels
- `sigma_0_2_default <- 0.2`
- `max_iter <- 100L`
- `credible_set_rho <- 0.95`
- `purity_threshold <- 0.50`
- `z_top_k <- 10`
- `jsd_threshold <- 0.15`
- `include_overall_pool <- FALSE`
- `softmax_temperature <- 1`
- `verbose_file_output <- FALSE`
- `write_snps_parquet <- FALSE`
- `write_confusion_bins <- TRUE`
- `write_tier_cs_metrics <- TRUE`
- `write_scaling_confusion_bins <- TRUE`

The confusion-bin threshold grid is nonuniform:

- 0 to 0.05 in increments of 0.005;
- 0.05 to 0.10 in increments of 0.01;
- 0.10 to 1.00 in increments of 0.02.

This puts more threshold resolution at low PIP values, which is where FPR changes
quickly in these sparse-positive fine-mapping evaluations.

The scheduled aggregation methods are:

- `max_elbo`
- `uniform`
- `elbo_softmax`
- `cluster_weight_credible` (clustering metric `credible_shift`, cut 0.05, Method B)
- `cluster_weight_max` (clustering metric `max_shift` / L-infinity, cut 0.10, Method B)

The paper-prep workflow now keeps all scheduled methods: `drop_aggs <- character(0)`.
The paper-prep workbook highlights `cluster_weight_credible` as the primary method.

The scaling-subset sizes are:

- `4`
- `8`
- `16`
- `32`
- `64`

The run-control workbook sets `scaling_restart_reps <- 50L`, but the runtime
scaling helper in `run_model.R` currently sets `n_reps <- 1L` for pure-lever
subsets. In practice, the scaling bins emitted by the runtime harness are
deterministic subsets, not 50 independent random restart subsets.

## 5. Genotype matrix selection

The run-control workbook builds a matrix catalog from sampled simulated genotype
matrices:

- helper: `test_susine/R/run_controls.R:60` (`build_data_matrix_catalog`)
- requested scenario: `"scenario_1"`
- summary path:
  `data/sampled_simulated_genotypes/scenario_sampling_summary.csv`

The workbook then either reads or constructs:

- `data/sampled_simulated_genotypes/scenario_sampling_summary_with_m1.csv`

It arranges available matrices by `M1`, takes 150 approximately evenly spaced
rows over that ordering, and then semi-joins the full matrix catalog to those
selected matrix IDs. The intent is not to sample 150 matrices uniformly at
random; it is to stratify coverage over the matrix-level `M1` statistic.

This selected catalog is passed to `make_job_config()` as `data_matrix_catalog`.

## 6. Dataset bundles and phenotype seeds

The core run-table construction happens in:

- `test_susine/R/run_controls.R:146` (`make_run_tables`)

For the ensemble study:

- matrix count is 150;
- base seeds are 1, 2, 3, 4;
- architecture is only `susie2_oligogenic`;
- `y_noise` is not meaningful for heritability-calibrated architectures and is
  collapsed to `NA`;
- `p_star` is not meaningful for `susie2_oligogenic` and is collapsed to `NA`.

This gives 150 x 4 = 600 dataset bundles.

Important seed transformation:

`make_run_tables()` does not use the base seed literally as the phenotype seed.
After assigning `dataset_bundle_id`, it sets:

```r
base_seed = phenotype_seed
phenotype_seed = as.integer(phenotype_seed * 10007L + dataset_bundle_id)
```

So the four requested seeds `1:4` become deterministic, dataset-specific
phenotype seeds after combining the base seed with the dataset bundle ID. This
prevents identical effect draws across different genotype matrices that share
the same base seed.

## 7. Oligogenic effect simulation

The SuSiE 2.0-style oligogenic architecture is implemented in:

- `test_susine/R/simulate_data.R:178`
  (`simulate_effect_sizes_susie2_oligogenic`)

For each locus with `p` SNPs:

- sample 23 non-overlapping causal SNP indices;
- divide them into tiers:
  - sparse tier: 3 SNPs;
  - oligogenic tier: 5 SNPs;
  - polygenic tier: 15 SNPs;
- draw raw effects:
  - sparse tier: `N(0, 1)`;
  - oligogenic tier: two-component Gaussian mixture:
    - 35 percent large, standard deviation 0.8;
    - 65 percent small, standard deviation 0.25;
  - polygenic tier: `N(0, 0.08^2)`;
- rescale tier energy so that the total squared-effect energy fractions are:
  - sparse: 0.50;
  - oligogenic: 0.35;
  - polygenic: 0.15.

The returned object includes:

- `beta`, the true effect vector;
- `causal_idx`, all 23 causal SNP indices;
- `effect_tier`, a per-SNP label: `sparse`, `oligogenic`, `polygenic`, or
  `none`.

The phenotype is generated by:

- `test_susine/R/simulate_data.R:257` (`simulate_phenotype_h2`)

For `susie2_oligogenic`, the target is:

- `h2_total = 0.25`

The phenotype simulation computes:

- signal: `X %*% beta`;
- signal variance: `var(signal)`;
- residual variance:

```r
sigma2 = var_signal * (1 - h2) / h2
```

Then:

```r
y = signal + rnorm(n, 0, sqrt(sigma2))
```

The phenotype noise seed is `phenotype_seed + 1`.

During actual model execution, dataset generation is handled by:

- `test_susine/R/run_model.R:525` (`generate_data_for_bundle`)

That function loads the selected matrix, regenerates effects and phenotype using
the transformed seed, and returns a data bundle containing `X`, `y`, `beta`,
`sigma2`, `causal_idx`, `effect_tier`, `matrix_id`, `architecture`,
`phenotype_seed`, and an empty annotation cache environment.

## 8. Annotation simulation

Annotations are generated on demand only for use cases that need them:

- prior mean strategy `functional_mu`; or
- inclusion prior strategy `functional_pi`.

In this ensemble study only `functional_mu` variants are used.

The annotation generator is:

- `test_susine/R/simulate_data.R:356` (`simulate_priors`)

It takes:

- `beta`;
- `annotation_r2`;
- `inflate_match`;
- `base_sigma2`;
- `effect_sd`;
- `annotation_seed`.

The key construction:

1. Identify causal and noncausal indices from `beta != 0`.
2. On causal positions, center the true causal effects and scale them.
3. Generate a centered noise vector orthogonal to the centered causal-effect
   direction.
4. Construct causal annotations as a signal-noise mixture:

```r
causal_annotation_centered =
  causal_scale * (
    sqrt(annotation_r2) * signal_unit +
    sqrt(1 - annotation_r2) * noise_unit
  )
```

5. Add back the causal mean.
6. On noncausal positions, generate centered noise with RMS:

```r
sqrt(inflate_match) * sqrt(theoretical_causal_mu0_var)
```

7. Return:
   - `mu_0`, the simulated signed annotation vector;
   - `sigma_0_2`, a vector initially set to `base_sigma2`;
   - `observed_r2`, the realized squared correlation between `mu_0` and `beta`
     on causal positions.

For the ensemble study:

- target `annotation_r2` values are 0, 0.3, and 0.5;
- `inflate_match` is always 0.95;
- annotation-consuming specs are expanded over all three annotation levels;
- annotation-agnostic specs get `annotation_r2 = NA` and `inflate_match = NA`.

The run table assigns annotation seeds in:

- `test_susine/R/run_controls.R:997`

For ordinary run-table rows, it builds a distinct mapping from
`(dataset_bundle_id, annotation_r2, inflate_match)` to `annotation_seed`.

Important truth-warm caveat:

The run-control workbook appends truth-warm rows manually after
`make_job_config()`. Those rows set `annotation_r2` and `inflate_match`, but the
appended code does not recompute `annotation_seed`. If the template row used for
truth-warm has `annotation_seed = NA`, annotated truth-warm fits may generate
annotations without a fixed annotation seed. The paper plots mainly use the
truth-warm vanilla reference, so this likely does not affect the plotted
annotation-dependent ensemble curves, but it matters if annotated truth-warm
rows are later interpreted.

## 9. Model use cases carried into the ensemble study

Use cases are defined in:

- `test_susine/R/use_cases.R:9` (`prior_spec_catalog`)

The ensemble run-control workbook uses three use cases.

### 9.1 `susine_vanilla`

Catalog row:

- backend: `susine`;
- prior mean strategy: `zero`;
- prior variance strategy: `fixed`;
- inclusion prior strategy: `uniform`;
- EB method: `NA`;
- supports annotations: false.

Despite the name, this is not the `susieR::susie` backend. It is a SuSiE-like
zero-prior-mean fit through the local `susine` backend. It is used as the
annotation-agnostic baseline family in Group A.

### 9.2 `susine_eb_clamped_scale_var_nonneg`

Catalog row:

- backend: `susine`;
- prior mean strategy: `functional_mu`;
- prior variance strategy: `eb`;
- inclusion prior strategy: `uniform`;
- EB method: `clamped_scale_var`;
- `c_nonneg = TRUE`;
- supports annotations: true.

This is the EB SuSiNE family in Group B. In `run_use_case()`, the model call
sets:

```r
prior_update_method = "clamped_scale_var"
c_nonneg = TRUE
```

and passes `mu_0 = c_value * annotation_vec`. If `c_value` is absent, the code
defaults it to 1 before calling the EB updater; the EB routine then updates the
scale/variance internally.

### 9.3 `susine_functional_mu`

Catalog row:

- backend: `susine`;
- prior mean strategy: `functional_mu`;
- prior variance strategy: `fixed`;
- inclusion prior strategy: `uniform`;
- EB method: `NA`;
- supports annotations: true.

This is the fixed-variance grid SuSiNE family in Group C. It consumes the signed
annotation vector as:

```r
mu_0 = c_value * annotation_vec
```

For fixed-variance `susine` fits, `sigma_0_2` is passed as a scalar multiplier.
The code comment in `run_use_case()` is explicit: `susine::initialize_priors`
multiplies by `var(y)` internally, so the harness passes the raw scalar, such as
`0.2`, rather than pre-multiplying by `var(y)`.

## 10. Exploration and aggregation catalogs

Exploration IDs are defined in:

- `test_susine/R/use_cases.R:65` (`exploration_catalog`)

Relevant ensemble-study exploration modes:

- `single`: one fit;
- `restart`: randomized initialization through `init_alpha`;
- `c_grid`: grid over the signed prior-mean scale `c`;
- `sigma_0_2_grid`: grid over fixed prior-variance scalars;
- `refine`: harness-level perturb-and-refit BFS;
- `cs_grid_refit`: fit a `(c, sigma_0_2)` grid under functional-mu, then refit
  with default priors from the converged effect fits.

Aggregation IDs are defined in:

- `test_susine/R/use_cases.R:82` (`aggregation_catalog`)

The core aggregation methods are:

- `max_elbo`;
- `uniform`;
- `elbo_softmax`;
- `cluster_weight_credible` (clustering metric `credible_shift`, cut 0.05, Method B);
- `cluster_weight_max` (clustering metric `max_shift` / L-infinity, cut 0.10, Method B).

Compatibility rules are in:

- `test_susine/R/use_cases.R:99` (`valid_exploration_for_prior`)

Important compatibility rules:

- `c_grid` applies to functional-mu fixed-variance models and is redundant for
  EB methods that optimize `c`;
- `sigma_0_2_grid` applies to fixed-variance models;
- `cs_grid_refit` is only valid for the `susine_functional_mu` fixed-variance
  setup.

## 11. Grid values and perturbation settings

### 11.1 `c` grids

The run-control workbook defines:

```r
c_grid_64 <- seq(0, 1.5, length.out = 64)
c_grid_8  <- seq(0, 1.5, length.out = 8)
c_grid_4  <- seq(0, 1.5, length.out = 4)
c_grid_3  <- seq(0, 1.5, length.out = 3)
```

So all `c` grids span 0 to 1.5, with different resolution.

### 11.2 `sigma_0_2` grids

The run-control workbook defines:

```r
sigma_grid_64 <- 10^seq(-2, 0, length.out = 64)
```

Then it snaps the nearest value to exactly `0.2` and sorts unique values.

The smaller grids are:

```r
sigma_grid_8 <- c(0.01, 0.03, 0.07, 0.1, 0.2, 0.4, 0.7, 1.0)
sigma_grid_4 <- c(0.02, 0.1, 0.2, 0.5)
sigma_grid_3 <- c(0.05, 0.2, 1.0)
```

The workbook asserts that `0.2` is present in all sigma grids.

### 11.3 Restart settings

Restart lists:

```r
restart_64 <- list(n_inits = 64L, alpha_concentration = 0.001)
restart_8  <- list(n_inits = 8L,  alpha_concentration = 0.001)
restart_4  <- list(n_inits = 4L,  alpha_concentration = 0.001)
restart_3  <- list(n_inits = 3L,  alpha_concentration = 0.001)
```

In `run_use_case()`, restart rows with `run_type = "warm"` generate
`init_alpha` from a Dirichlet distribution:

```r
dirichlet_matrix(L, rep(alpha_concentration, p))
```

With `alpha_concentration = 0.001`, these warm starts are very sparse over
variants for each effect.

The first restart row has `run_type = "default"` and subsequent restart rows are
`run_type = "warm"`; this is assigned in `axis_table_for_method()` in
`test_susine/R/run_controls.R:350`.

### 11.4 Refine settings

Refine lists:

```r
refine_32 <- list(n_steps = 32L, cs_source = "filtered", purity_threshold = 0.95)
refine_8  <- list(n_steps = 8L,  cs_source = "filtered", purity_threshold = 0.95)
refine_4  <- list(n_steps = 4L,  cs_source = "filtered", purity_threshold = 0.95)
refine_3  <- list(n_steps = 3L,  cs_source = "filtered", purity_threshold = 0.95)
refine_2  <- list(n_steps = 2L,  cs_source = "filtered", purity_threshold = 0.95)
```

The general model-evaluation purity threshold is 0.50, but refine branching uses
0.95. That is deliberate: output metrics keep CSs at purity >= 0.50, but refine
only blocks high-purity CSs when generating children.

## 12. Exact ensemble specification table

The run-control workbook constructs 16 non-baseline ensemble specs plus a
standalone baseline-single spec. These are appended into `all_specs`.

The names below are the exact `spec_name` values used downstream.

### 12.1 Group A: annotation-agnostic SuSiE-like ensembles

All Group A specs use:

- use case: `susine_vanilla`;
- annotations: none (`annotation_r2 = NA`);
- fixed prior variance unless explicitly gridded;
- default fixed `sigma_0_2 = 0.2` for non-sigma-grid specs.

| spec | exploration | mode | nominal K | axes |
|---|---:|---:|---:|---|
| `A-R` | `restart` | `separate` | 64 | 64 restarts |
| `A-F` | `refine` | `separate` | 32 | 32 stored refine steps |
| `A-S` | `sigma_0_2_grid` | `separate` | 64 | 64 sigma values |
| `A-RF` | `restart x refine` | `intersect` | 64 | 8 restarts x 4 refine steps = 32 stored rows |
| `A-RS` | `restart x sigma_0_2_grid` | `intersect` | 64 | 8 restarts x 8 sigma values = 64 rows |
| `A-FS` | `refine x sigma_0_2_grid` | `intersect` | 64 | 4 refine steps x 8 sigma values = 32 stored rows |
| `A-RFS` | `restart x refine x sigma_0_2_grid` | `intersect` | 64 | 4 restarts x 2 refine steps x 4 sigma values = 32 stored rows |

For refine-containing specs, stored row count is not the same as model-call
count. Non-root refine steps perform a perturb fit and a refit; only the refit is
stored as an ensemble member.

### 12.2 Group B: EB SuSiNE ensembles

All Group B specs use:

- use case: `susine_eb_clamped_scale_var_nonneg`;
- annotations: `annotation_r2` in 0, 0.3, 0.5 and `inflate_match = 0.95`;
- EB prior update method: `clamped_scale_var`;
- nonnegative annotation scale.

| spec | exploration | mode | nominal K | axes |
|---|---:|---:|---:|---|
| `B-R` | `restart` | `separate` | 64 | 64 restarts |
| `B-F` | `refine` | `separate` | 32 | 32 stored refine steps |
| `B-RF` | `restart x refine` | `intersect` | 64 | 8 restarts x 4 refine steps = 32 stored rows |

Each Group B spec is evaluated separately at all three annotation qualities.

### 12.3 Group C: fixed-variance grid SuSiNE ensembles

All Group C specs use:

- use case: `susine_functional_mu`;
- annotations: `annotation_r2` in 0, 0.3, 0.5 and `inflate_match = 0.95`;
- fixed-variance SuSiNE;
- uniform inclusion prior.

| spec | exploration | mode | nominal K | axes |
|---|---:|---:|---:|---|
| `C-C` | `c_grid` | `separate` | 64 | 64 c values |
| `C-CS` | `c_grid x sigma_0_2_grid` | `intersect` | 64 | 8 c values x 8 sigma values |
| `C-CSR` | `cs_grid_refit` | `separate` | 64 | 8 c values x 8 sigma values, then default-prior refit |
| `C-CSFR` | `c_grid x sigma_0_2_grid x refine x restart` | `intersect` | 81 | 3 c values x 3 sigma values x 3 refine steps x 3 restarts |
| `C-CFR` | `c_grid x refine x restart` | `intersect` | 64 | 4 c values x 2 refine steps x 4 restarts = 32 stored rows |
| `C-CF` | `c_grid x refine` | `intersect` | 64 | 8 c values x 4 refine steps = 32 stored rows |

Each Group C spec is evaluated separately at all three annotation qualities.

### 12.4 Baseline single fits

The run-control workbook appends:

```r
spec_baseline_single <- list(
  name = "baseline-single",
  use_case_ids = c(
    "susine_vanilla",
    "susine_eb_clamped_scale_var_nonneg",
    "susine_functional_mu"
  ),
  exploration_methods = "single",
  exploration_mode = "separate",
  K = 1L
)
```

Consequences:

- `baseline-single / susine_vanilla` has no annotations and fixed
  `sigma_0_2 = 0.2`; this is the paper's SuSiE baseline reference.
- `baseline-single / susine_eb_clamped_scale_var_nonneg` is expanded over the
  three annotation qualities.
- `baseline-single / susine_functional_mu` is expanded over the three annotation
  qualities, uses fixed `sigma_0_2 = 0.2`, and because no `c_value` is set by the
  single-fit exploration, `run_use_case()` defaults `c_value` to 1.

### 12.5 Truth-warm rows

After `make_job_config()`, the run-control workbook manually appends a
`truth-warm` spec. This is not part of the heterogeneous `all_specs` list; it is
added by directly binding rows onto `cfg$tables$runs`.

The truth-warm use cases are:

| use case | `c_value` | `sigma_0_2_scalar` |
|---|---:|---:|
| `susine_vanilla` | `NA` | `"0.2"` |
| `susine_eb_clamped_scale_var_nonneg` | `NA` | `NA` |
| `susine_functional_mu` | `0.5` | `"0.2"` |

For annotation-supporting truth-warm use cases, rows are expanded over
annotation levels 0, 0.3, and 0.5. For `susine_vanilla`, there is a single
annotation-agnostic row per dataset bundle.

In `run_use_case()`, `warm_method = "truth_warm"` creates an `L x p` initial
alpha matrix:

- all rows start uniform at `1 / p`;
- the top `min(L, n_causal)` causal variants by `abs(beta)` are identified;
- the first `n_warm` effects are made one-hot at those top causal variants.

With `L = 10` and 23 causal variants, this means the first 10 effects are
initialized exactly on the 10 largest-effect causal SNPs.

## 13. Scheduled row counts

Given 600 dataset bundles, the scheduled run-table size can be derived from the
spec grid.

Per dataset bundle:

- Group A ensemble rows:
  - `A-R`: 64
  - `A-F`: 32
  - `A-S`: 64
  - `A-RF`: 32
  - `A-RS`: 64
  - `A-FS`: 32
  - `A-RFS`: 32
  - subtotal: 320
- Group B ensemble rows:
  - 128 rows per annotation level;
  - 3 annotation levels;
  - subtotal: 384
- Group C ensemble rows:
  - 337 rows per annotation level;
  - 3 annotation levels;
  - subtotal: 1011
- `baseline-single` rows:
  - 1 annotation-agnostic `susine_vanilla`;
  - 3 EB rows;
  - 3 functional-mu rows;
  - subtotal: 7
- `truth-warm` rows:
  - 1 annotation-agnostic `susine_vanilla`;
  - 3 EB rows;
  - 3 functional-mu rows;
  - subtotal: 7

Total scheduled run rows per dataset bundle:

```text
320 + 384 + 1011 + 7 + 7 = 1729
```

Across 600 dataset bundles:

```text
1729 * 600 = 1,037,400 scheduled run rows
```

This is a scheduled-row count, not an exact model-call count. Actual model calls
can differ because:

- refine rows beyond the root perform two model fits, perturb plus refit, while
  only the refit is stored;
- refine BFS can exhaust its queue early, so fewer stored refine fits may be
  emitted than were scheduled;
- duplicate run configurations can be served from the per-dataset execution
  cache rather than refit.

## 14. Run-table construction mechanics

`make_job_config()` supports a heterogeneous `exploration_specs` list. For the
ensemble study, every spec in `all_specs` is passed into that argument. Internally
it calls `make_run_tables()` once per spec, tags rows with `spec_name`, and binds
all rows together.

Important run-table mechanics:

### 14.1 Annotation expansion

For each use case:

- if `supports_annotation = FALSE`, the annotation grid is one row:
  `annotation_r2 = NA`, `inflate_match = NA`;
- if `supports_annotation = TRUE`, the annotation grid is all distinct rows of
  the provided prior-quality table:
  `(0, 0.95)`, `(0.3, 0.95)`, `(0.5, 0.95)`.

This expansion is controlled by `annotation_grid_for_spec()` in
`make_run_tables()`.

### 14.2 Sigma expansion

After initial run rows are built, `make_job_config()` expands fixed-variance
specs over `sigma_0_2_scalars` unless the row already has a preset
`sigma_0_2_scalar` from `sigma_0_2_grid`.

For the ensemble study:

- global scalar is `0.2`;
- sigma-grid specs already carry their sigma value and are not crossed with the
  global scalar again;
- fixed-variance non-sigma-grid specs get `sigma_0_2_scalar = "0.2"`;
- EB specs get `sigma_0_2_scalar = NA`.

This happens in `test_susine/R/run_controls.R:904`.

### 14.3 Intersect-mode row crossing

For intersect-mode specs, `build_exploration_groups()` builds one axis table per
exploration method, then crosses them with `tidyr::crossing()`.

For refine-containing intersect specs, the code treats each non-root refine row
as roughly twice as expensive as a normal row, because non-root refine performs a
perturb fit and a refit. The validation therefore allows:

- row count equal to nominal `K`; or
- row count times 2 equal to nominal `K`.

There is a special case for an inferred refine count of 3: it is not halved.

### 14.4 Group keys

The `group_key` is used to decide which individual fits belong to the same
ensemble for aggregation. It includes:

- `L`;
- annotation r2;
- `inflate_match`;
- exploration mode;
- exploration-method combination.

For example:

```text
L=10|r2=0.3|inflate=0.95|explore=intersect:c_gridxsigma_0_2_grid
```

For restart warm-start grids, older code can append warm-method and alpha fields
to the group key, but in this workbook restart alpha is supplied directly in
`restart_settings`, not through an additional warm-start grid crossing.

### 14.5 Task assignment

Rows are assigned to tasks by dataset bundle, so all runs for a dataset bundle
stay together. This matters scientifically because aggregation is performed
within a dataset bundle over all ensemble members in a group.

Slurm execution details are otherwise outside the scope of this note.

## 15. Dataset execution loop

Task execution enters through:

- `test_susine/R/run_model.R:13` (`run_task`)

The per-dataset work is in:

- `test_susine/R/run_model.R:622` (`execute_dataset_bundle`)

For each dataset bundle:

1. Regenerate the simulated data bundle with `generate_data_for_bundle()`.
2. Construct the top-8 causal mask:
   - sort `data_bundle$causal_idx` by `abs(beta)` decreasing;
   - keep the first 8;
   - use this mask for confusion-bin scoring.
3. Compute dataset-level metrics once with `compute_dataset_metrics()`, using
   `z_top_k = 10`.
4. Create a per-dataset `execution_cache`.
5. Iterate over use cases and run rows.
6. Store primary PIPs, ELBOs, run metadata, and fitted values by ensemble group.
7. After individual fits, aggregate PIPs within each group and write aggregated
   confusion bins.
8. If scaling bins are enabled, compute deterministic sub-ensemble confusion
   bins from the already-computed PIPs.

Because the execution cache is per dataset bundle, repeated identical model
calls inside the same dataset can be avoided. The cache key includes use-case
settings, annotation settings, `c`, `sigma_0_2`, restart/warm-start fields, and
blocked indices.

Refine refits with direct `init_alpha_override` skip the cache, because the
override is not encoded in the cache key.

## 16. Model fitting details

The model-fitting function is:

- `test_susine/R/run_model.R:1623` (`run_use_case`)

For each run row, it extracts:

- `L`;
- max iterations and tolerance;
- backend (`susine` or `susieR`);
- prior mean strategy;
- prior variance strategy;
- inclusion prior strategy;
- EB method;
- annotation settings;
- `c_value`;
- `tau_value`;
- restart/refine metadata;
- `run_type`;
- `warm_method`.

### 16.1 Annotation use inside model fitting

If the prior mean strategy is `functional_mu`, the function calls
`get_priors_cached()` to generate or retrieve the simulated annotation vector.
Then:

```r
mu_0 = c_value * annotation_vec
```

If `c_value` is missing or non-finite, it defaults to 1.

The inclusion prior remains uniform in this ensemble study. Functional-pi code
exists, but it is not used by the ensemble specs.

### 16.2 Restart warm starts

If `run_type = "warm"` and `warm_method = "init_alpha"`, the code creates a warm
initial alpha matrix.

With `alpha_concentration = 0.001`, each effect's initialization is a sparse
Dirichlet draw over SNPs.

If `alpha_concentration = 0`, the code would use one-hot random SNP choices,
but the ensemble workbook uses 0.001, not 0.

### 16.3 Truth warm starts

If `warm_method = "truth_warm"`, the code initializes the first `L` effects on
the largest true causal variants by absolute effect size, as described above.

This is a truth-aware reference, not a deployable method.

### 16.4 Blocking logic for refine

If `blocked_idx` is supplied:

- prior weights at blocked variants are set to 0 and renormalized;
- if a warm initial alpha exists, the blocked columns are also zeroed and rows
  are renormalized.

The perturb fit in refine uses blocked indices. The stored refit restores the
unblocked prior weights but uses the perturb fit's alpha matrix as the direct
warm start.

### 16.5 Fixed prior variance

For fixed-variance `susine` fits:

```r
args$sigma_0_2 <- as.numeric(sigma_scalar_fixed)
args$prior_update_method <- "none"
```

The raw scalar is passed, not multiplied by `var(y)`.

### 16.6 EB prior variance and scale

For EB `susine` fits:

```r
args$prior_update_method <- eb_method
args$c_nonneg <- c_nonneg_flag
```

For Group B, this means:

```r
prior_update_method = "clamped_scale_var"
c_nonneg = TRUE
```

## 17. Refine algorithm

Refine execution is implemented inside `execute_dataset_bundle()`, starting
around `test_susine/R/run_model.R:956`.

The algorithm is a breadth-first search over blocked credible-set combinations.

For each parent group, where "parent" means the non-refine axes are fixed
(for example one restart, one sigma value, one c value), the code:

1. Starts a queue with the root node:
   - `blocked = integer(0)`;
   - no parent fit.
2. Processes requested refine rows in order.
3. For the root node:
   - fit the model with no blocking;
   - store the root fit as the first ensemble member.
4. For each non-root node:
   - perturb step:
     - fit with prior weights blocked on the node's accumulated CS indices;
     - do not write this perturb fit as an ensemble member;
     - extract `perturb_fit$effect_fits$alpha`;
   - refit step:
     - restore unblocked weights;
     - inject `init_alpha_override = perturb_alpha`;
     - write this refit as the ensemble member;
     - wall time includes perturb plus refit cost.
5. From the stored fit's filtered credible sets, collect branch candidates with
   purity >= 0.95.
6. For each eligible CS, create a child node whose blocked set is the union of
   the parent blocked set and that CS.
7. Avoid duplicate blocked sets through a `seen_block_sets` environment.
8. Continue until the requested row count is processed or the BFS queue is
   exhausted.

The code records refine depth diagnostics:

- `n_refine_requested`;
- `n_refine_processed`;
- parent key;
- group key.

If the queue exhausts early, the number of stored refine fits can be lower than
the scheduled `n_steps`. The collection workflow has specific backfills to
handle cases where a pure-refine group only produced a single fit.

## 18. `C-CSR`: CS-grid refit

The special `cs_grid_refit` mode is used by `C-CSR`.

Run-table construction for `cs_grid_refit` is in:

- `test_susine/R/run_controls.R:451`

It crosses the provided `c_grid_values` and `sigma_0_2_grid_values`, sets:

```r
run_type = "warm"
warm_method = "cs_grid_refit"
exploration_mode = "separate"
exploration_methods = "cs_grid_refit"
```

For `C-CSR`, the grid is 8 c values x 8 sigma values = 64 source settings.

During model execution, the source fit is a normal functional-mu fit at that
`c` and `sigma_0_2`. Then `run_default_prior_refit()` is called:

- file: `test_susine/R/run_model.R:1587`

The refit:

- uses `mu_0 = 0`;
- uses `sigma_0_2 = 0.2`;
- uses uniform prior inclusion weights;
- disables prior updates;
- passes exact converged effect fits from the source fit as
  `init_effect_fits = build_exact_init_effect_fits(source_fit)`.

So `C-CSR` asks: after functional-mu grid exploration finds a basin, does the
default zero-prior model retain that basin when exactly warm-started from the
source effect fits?

## 19. Individual model evaluation

Model evaluation is handled by:

- `test_susine/R/evaluation_metrics.R:315` (`evaluate_model`)

For each fit:

1. Extract the `L x p` alpha matrix.
2. For each effect, construct a credible set at cumulative mass `rho = 0.95`.
3. Compute credible-set metrics:
   - size;
   - purity;
   - coverage;
   - effect PIP entropy;
   - core95 entropy;
   - effective signal count;
   - tail inflation;
   - accuracy ratio.
4. Create unfiltered and purity-filtered effect tables.
5. Filtered effects require purity >= `purity_threshold`, which is 0.50 in this
   ensemble job.
6. Compute model-level summaries:
   - nominal and effective `L`;
   - mean CS size;
   - mean purity;
   - mean coverage;
   - CS power;
   - overlap rate;
   - AUPRC;
   - cross entropy;
   - estimated genetic variance fraction `hg2`;
   - TPR at several FPR thresholds.

Important nuance:

- The `evaluate_model()` AUPRC is computed against all causal indices.
- The confusion-bin outputs used by the collection workflow and paper PR/AUPRC
  plots are computed against the top-8 causal mask.

So the paper's pooled AUPRC is based on confusion bins, not directly on the
per-fit `model_metrics$AUPRC` column.

## 20. Confusion bins and top-8 causal scoring

Confusion bins are built by:

- `test_susine/R/run_model.R:218` (`compute_confusion_bins`)

The function takes:

- a PIP vector;
- a causal vector;
- PIP bucket breaks;
- optional `causal_mask`.

If `causal_mask` is supplied, it builds the scored-positive vector and then
DROPS the weak (true-causal-but-not-masked) variants:

```r
scored <- as.integer(seq_len(n) %in% causal_mask)
keep <- !(causal == 1L & scored == 0L)
pip <- pip[keep]
causal <- scored[keep]
```

In `execute_dataset_bundle()`, the mask is:

- the top 8 true causal SNPs by `abs(beta)`;
- selected from the full 23-causal support.

Therefore:

- positives = 8 largest-effect causal SNPs per dataset;
- true non-causal SNPs are negatives;
- the smaller-effect (weak) causal SNPs are DROPPED entirely -- neither
  positive nor negative -- so AUPRC is not penalized for failing to recover
  them.

The manuscript baseline figure caption for the 5-arm refit chain similarly says
positives are the top 8 causal SNPs per dataset by absolute beta, with the weak
(non-top-8) causals excluded from the negative set. The confusion-bin function
implements exactly that exclusion via the `keep` filter above.

The confusion bins record, for each PIP bucket lower edge:

- `pip_threshold`;
- `n_causal_at_bucket`;
- `n_noncausal_at_bucket`.

The collection code pools these bucket counts across datasets before computing
pooled AUPRC and TPR.

## 21. Primary full-ensemble aggregation

Primary aggregation happens in `execute_dataset_bundle()` after all individual
fits in a group are stored:

- `test_susine/R/run_model.R:1101`
- aggregation helper:
  `test_susine/R/run_model.R:1931` (`aggregate_use_case_pips`)

For each ensemble group, the code has:

- list of PIP vectors;
- vector of final ELBOs;
- run metadata;
- fitted values.

It aggregates PIPs by each scheduled method.

### 21.1 `uniform`

Mean PIP across fits:

```r
colMeans(pips_mat)
```

Here `pips_mat` is `n_fits x p`.

### 21.2 `max_elbo`

Select the PIP vector from the fit with highest final ELBO:

```r
pips_mat[which.max(elbo_vec), ]
```

### 21.3 `elbo_softmax`

Convert ELBOs to softmax weights with temperature 1:

```r
w = exp(elbo - max(elbo))
w = w / sum(w)
```

Then aggregate:

```r
crossprod(w, pips_mat)
```

### 21.4 `cluster_weight_credible` and `cluster_weight_max`

The two scheduled cluster-weight methods each build their own dendrogram with a
single-variant PIP-shift metric and complete-linkage clustering, then aggregate
with Method B (frequency-free):

- `cluster_weight_credible`: clustering metric `credible_shift` =
  max_j[ max(p_j,q_j) * |p_j - q_j| ] (largest credibility-weighted single-variant
  PIP shift), complete-linkage cut at threshold 0.05;
- `cluster_weight_max`: clustering metric `max_shift` = max_j |p_j - q_j|
  (Chebyshev / L-infinity), complete-linkage cut at 0.10.

The per-method (metric, threshold) pairs are locked in
`.cluster_weight_specs` (`run_model.R` ~lines 2143-2146). Each method builds its
dendrogram via `prepare_pip_similarity_cache(pips_mat, metric = spec$metric)`,
then aggregates via `.aggregate_cluster_method()` ->
`.cluster_weights_from_hc(..., aggregation = "method_b")`.

**Method B** (`.cluster_weights_from_hc`, `run_model.R` ~lines 2086-2126):

- cut the dendrogram at the method's threshold;
- each cluster gets weight proportional to `exp(max-ELBO-in-cluster)`
  (softmax over per-cluster max ELBOs);
- that cluster weight is split within the cluster across ALL member fits by
  ELBO-softmax;
- there is no single nominee and no 1/frequency correction; the weight vector
  `w_full` spans all fits.

A degenerate single-fit group falls back to plain ELBO-softmax.

The final aggregated PIP vector gets an `ess` attribute:

```r
1 / sum(w_full^2)
```

The legacy JSD-0.15 + frequency-corrected `cluster_weight` (Method C: max-ELBO
nominee per JSD cluster, weight divided by cluster frequency) still exists in
code for backward compatibility, but it is not used by the current
ensemble-simulation paper workflow or the current real-data paper workflow.

## 22. Scaling-bin aggregation and the cluster-weight discrepancy

Scaling-bin subsets are computed at runtime by:

- `test_susine/R/run_model.R:2454`
  (`compute_scaling_confusion_bins_for_group`)

That function uses:

- `test_susine/R/collect_results.R:830` (`aggregate_pip_matrix`)

For `uniform`, `max_elbo`, and `elbo_softmax`, this is equivalent to the primary
aggregation logic.

For the scheduled `cluster_weight_credible` and `cluster_weight_max` methods,
both the full-ensemble path and the scaling-bin path now route through the same
Method B logic (each method's own metric + threshold, frequency-free cluster
weights split by within-cluster ELBO-softmax). The older full-vs-scaling
cluster-weight discrepancy -- where the scaling helper weighted clusters
uniformly (`1/K`) without representatives or frequency adjustment -- is therefore
obsolete for the scheduled methods; both paths agree at full resolution.

(The legacy JSD-0.15 + frequency `cluster_weight` method is the one for which
the historical discrepancy was defined; current paper workflows avoid it.)

## 23. Scaling-subset construction

Runtime scaling bins are emitted only if `write_scaling_confusion_bins = TRUE`,
which it is in the run-control workbook.

The scaling helper receives all PIPs and ELBOs for an ensemble group and returns
confusion bins for sub-ensembles.

### 23.1 Pure-lever specs

Pure-lever specs have one exploration method:

- pure restart: `A-R`, `B-R`;
- pure refine: `A-F`, `B-F`;
- pure sigma grid: `A-S`;
- pure c grid: `C-C`.

The helper identifies the lever type in priority order:

1. restart;
2. refine;
3. sigma;
4. c grid.

It then evaluates `n_ensemble` values from:

```r
c(4, 8, 16, 32, 64)
```

but only those less than or equal to the planned full size for that lever.

Subset selection:

- restart: first `n` restart IDs in sorted order;
- refine: first `n` refine steps in sorted order;
- sigma grid: evenly spaced values over sorted unique sigma values;
- c grid: evenly spaced values over sorted unique c values.

Although the function signature has `n_restart_reps = 50`, the implementation
sets `n_reps <- 1L`, so only one deterministic subset is emitted per size.

### 23.2 Two-axis interaction specs

For two-axis interaction specs, the helper emits two resolution levels:

- half resolution:
  - each axis size is `ceil(planned_axis_size / 2)`;
- full resolution:
  - each axis size is the planned full axis size.

Axis values are selected as:

- first values for restart/refine axes;
- evenly spaced values for c/sigma axes.

The resolution label is made from planned sizes and sorted by descending planned
axis size, so labels are stable across datasets.

Examples:

- `C-CS` full: `8x8`;
- `C-CS` half: `4x4`;
- `A-RF` full: `8x4`;
- `A-RF` half: `4x2`;
- `A-RS` full: `8x8`;
- `A-FS` full: `8x4` or `4x8` depending on planned-axis sorting.

### 23.3 Three-or-more-axis interaction specs

The runtime scaling helper skips specs with more than two exploration axes:

```r
if (length(method_ids) > 2L) return(dplyr::bind_rows(results))
```

So subscaling bins are not generated for:

- `A-RFS`;
- `C-CSFR`;
- `C-CFR`;

and any other three-or-more-axis spec.

`C-CF` has two axes (`c_grid` and `refine`), so it is eligible.

## 24. Additional runtime diagnostics

### 24.1 Multimodality metrics

The runtime computes multimodality metrics in:

- `test_susine/R/run_model.R:2029`

For groups with multiple fits, it computes:

- pairwise PIP-vector JSD mean/median/max;
- top-`k` Jaccard similarity across fits, with `k = z_top_k = 10`;
- mean PIP variance;
- number of clusters at JSD 0.15;
- number of clusters at JSD 0.50.

These diagnostics are written to `multimodal_metrics`.

### 24.2 Aggregated `hg2`

Aggregated fitted-value heritability is computed in:

- `test_susine/R/run_model.R:2850` (`compute_hg2_by_agg`)

The code now emits a corrected "fair" local-genetic-variance decomposition rather
than only the old posterior-mean point estimate. For each aggregation method it
combines fitted values with the method's aggregation weights `w` and computes:

- `hg2_postmean` = var(weighted posterior-mean fitted y) / var(y);
- `hg2_between_fit` = Σ_k w_k var(fy_k − fy_agg) / var(y) (mixture dispersion);
- `hg2_uncertainty` = Σ_k w_k * within-fit posterior-variance correction
  (`tr(Σ Cov(b|y))/var(y)`, weighted over positively-weighted fits);
- the REPORTED estimand
  `hg2_expected_pve = hg2_postmean + hg2_between_fit + hg2_uncertainty`
  = E[var(Xβ)|y] / var(y).

The additive decomposition is kept exact (not clamped at 1); clamping happens at
the display layer. The legacy point-estimate columns `hg2` (=
var(weighted fitted y)/var(y), the old estimand) and `hg2_weighted_pve` (=
Σ_k w_k var(fy_k)/var(y)) are still emitted for backward compatibility. The
decomposition engine lives in the new file `test_susine/R/heritability.R`
(`hg2_uncertainty_scalar()` supplies the per-fit correction). Provenance commits
`4e0721b` / `236820b` (2026-06-09); see
`refs/decisions/heritability_estimand_decision_2026-06-09.md`.

For cluster-weight aggregation, each scheduled method uses its own metric +
threshold dendrogram and Method B weights (same logic as primary `run_model.R`
cluster-weight PIP aggregation). The paper-prep workbook consumes
`hg2_expected_pve` / `hg2_postmean`.

## 25. Runtime output shape

With `verbose_file_output = FALSE`, the harness avoids per-run file storms and
buffers task outputs into staged flush CSVs.

Important outputs written per task include:

- `model_metrics`;
- `effect_metrics`;
- `confusion_bins`;
- `dataset_metrics`;
- `validation`;
- `prior_diagnostics` if enabled;
- `tier_cs_metrics`;
- `multimodal_metrics`;
- `refine_depth`;
- `hg2_by_agg`;
- `scaling_bins`.

Because `write_snps_parquet = FALSE`, per-SNP PIP parquet output is not written
for this full ensemble job. That limits the ability to reconstruct missing
aggregation rows from SNP-level data after the fact. The collection workflow has
an SNP-based backfill path, but it will not help if no SNP parquet files exist.

Individual fits write model metrics and confusion bins. Aggregated ensembles
write confusion bins but not full model metrics, because there is no fitted model
object for an aggregated PIP vector. Aggregated `hg2` is written separately via
`hg2_by_agg`.

## 26. Collection workflow

The collection workbook is:

- `test_susine/vignettes/simulation pipeline/collect_results_workbook_ensemble_scaling.Rmd`

At the top of the currently inspected file:

```r
job_name      <- "ensemble_scaling_full"
parent_job_id <- "53520405"
output_root   <- here("output")
```

Again, this parent ID differs from the paper-prep and paper-rendering workbooks,
which use `53547760`.

The collection workflow:

1. Reads job configuration from:
   - `output/run_history/<job_name>/<parent_job_id>/job_config.json`
2. Defines the job output directory:
   - `output/slurm_output/<job_name>/<parent_job_id>`
3. Indexes staged task outputs.
4. Validates that required output types are present and readable.
5. Checks validation rows for model-level problems.
6. Merges staged CSVs by output type.
7. Joins run metadata from the job config.
8. Backfills certain missing aggregation/scaling rows.
9. Pools confusion bins into AUPRC/TPR summaries.
10. Constructs oracle rows.
11. Exports consolidated CSVs.

## 27. Collection validation

The collection workbook defines required staged output types:

- `model_metrics`;
- `confusion_bins`;
- `dataset_metrics`;
- `validation`.

It writes diagnostics under:

- `<job_dir>/collect_diagnostics`

The validation stage checks:

- expected task IDs from `cfg$tables$tasks`;
- discovered staged files;
- file readability;
- `validation` rows with `has_issues`.

If required files are missing or unreadable, or if validation reports problems,
the workbook stops before metric consolidation.

## 28. Metadata joins during collection

The collection workbook creates `run_info` from `cfg$tables$runs`, selecting
columns such as:

- `run_id`;
- `task_id`;
- `dataset_bundle_id`;
- `spec_name`;
- `use_case_id`;
- `exploration_methods`;
- `exploration_mode`;
- `sigma_0_2_scalar`;
- `warm_method`;
- `c_value`;
- `refine_step`;
- `restart_id`;
- `annotation_r2`;
- `inflate_match`;
- `group_key`.

It joins `run_info` into:

- `model_metrics`;
- `tier_cs_metrics`;
- individual confusion bins;
- aggregated confusion bins.

Confusion bins are split as:

- individual rows:
  - `agg_method` missing;
- aggregated rows:
  - `agg_method` non-missing.

Then both are normalized and later bound together for export.

## 29. Collection backfills

The collection workflow applies several backfills from
`test_susine/R/ensemble_scaling_utils.R`.

### 29.1 Single-fit refine aggregation backfill

Function:

- `test_susine/R/ensemble_scaling_utils.R:423`
  (`backfill_single_fit_refine_agg_confusion`)

Purpose:

- Some refine groups can produce only one individual fit because the BFS queue
  exhausts before branching.
- Aggregation over one fit is degenerate but still should have rows for each
  aggregation method.

Mechanism:

- find refine groups with exactly one individual fit;
- identify missing aggregation methods;
- copy the individual fit's confusion bins;
- mark them as `explore_method = "aggregation"` and `agg_method = <method>`.

### 29.2 Terminal scaling aggregation backfill

Function:

- `test_susine/R/ensemble_scaling_utils.R:540`
  (`backfill_terminal_scaling_agg_confusion`)

Purpose:

- If a full-group aggregate confusion row is missing but terminal scaling bins
  exist, use the terminal scaling bins to fill the aggregate row.

Mechanism:

- for pure scaling rows, take the maximum `n_ensemble`;
- for interaction scaling rows, take the largest resolution area;
- map those bins back onto the group metadata;
- add missing aggregate confusion rows.

Because terminal scaling bins use `aggregate_pip_matrix()`, this backfill can
import the scaling implementation of cluster-weight into the main aggregate
confusion table for repaired groups.

### 29.3 SNP-based aggregation backfill

Function:

- `test_susine/R/ensemble_scaling_utils.R:700`
  (`backfill_missing_agg_confusion_from_snps`)

Purpose:

- If per-SNP parquet files exist, reconstruct missing aggregate PIP vectors and
  confusion bins.

In the ensemble run-control workbook, `write_snps_parquet = FALSE`, so this path
is expected to do nothing for the full run unless SNP parquet files were produced
by some other version of the job.

### 29.4 Pure-refine scaling-bin backfill

Function:

- `test_susine/R/ensemble_scaling_utils.R:72`
  (`backfill_pure_refine_scaling_bins`)

Purpose:

- Pure-refine scaling curves request multiple `n_ensemble` sizes, but the actual
  BFS may produce fewer stored refine fits than planned.
- This function fills missing pure-refine scaling bins for sizes above the
  actual count by using the full refine aggregate bins.

## 30. AUPRC and TPR computation

The collection workbook computes metrics from pooled confusion bins.

The core pooled-AUPRC implementation is:

- `test_susine/R/collect_results.R:809` (`auprc_from_pooled_bins`)
- `test_susine/R/visualize_results.R:131`
  (`compute_auprc_from_pooled_confusion`)

The algorithm:

1. Sort PIP thresholds descending.
2. Compute cumulative true positives and false positives:

```r
cum_tp = cumsum(n_causal_at_bucket)
cum_fp = cumsum(n_noncausal_at_bucket)
```

3. Compute:

```r
precision = cum_tp / (cum_tp + cum_fp)
recall = cum_tp / total_positives
```

4. Apply the package's average-precision / step-function AUPRC helper:
   `sum(precision_k * delta_recall_k)`, with no `(0, 1)` anchor and no
   trapezoidal interpolation.

TPR@FPR=0.05 is computed by:

- `test_susine/R/visualize_results.R:337`
  (`compute_tpr05_from_pooled_confusion`)

The collection workbook computes both:

- per-dataset metrics;
- pooled metrics where counts are summed across datasets before metric
  computation.

For paper figures, the key AUPRC values come from pooled confusion bins, not from
averaging per-dataset AUPRC values.

## 31. Collection grouping structure

The collection workbook creates:

### 31.1 Per-dataset individual AUPRC

Group variables include:

- `dataset_bundle_id`;
- `run_id`;
- `spec_name`;
- `use_case_id`;
- `annotation_r2`;
- `group_key`;
- `explore_method`;
- `agg_method`;
- `variant_id`.

### 31.2 Per-dataset aggregated AUPRC

Same general per-dataset structure, but `agg_method` identifies the aggregation
method and `explore_method = "aggregation"`.

### 31.3 Pooled aggregated AUPRC

Overall grouping:

- `spec_name`;
- `use_case_id`;
- `agg_method`.

By-annotation grouping:

- `spec_name`;
- `use_case_id`;
- `agg_method`;
- `annotation_r2`.

### 31.4 Pooled individual AUPRC

Overall grouping:

- `spec_name`;
- `use_case_id`.

By-annotation grouping:

- `spec_name`;
- `use_case_id`;
- `annotation_r2`.

The same grouping logic is repeated for TPR@FPR=0.05.

## 32. Oracle construction

The oracle is built in the collection workbook's `oracle` chunk.

It is not a model that was run. It is a truth-aware upper bound constructed from
the already-run individual fits.

Mechanism:

1. Start from per-dataset individual AUPRC rows.
2. For each `(dataset_bundle_id, spec_name, use_case_id, annotation_r2,
   group_key)`, choose the individual run with maximum AUPRC.
3. Label it:

```r
agg_method = "oracle"
explore_method = "aggregation"
```

4. Pull the corresponding individual confusion bins.
5. Pool those selected bins to compute oracle pooled AUPRC and TPR@FPR=0.05.
6. Append oracle rows to the aggregated metric tables.

Because oracle selection is by per-dataset AUPRC, it uses ground truth and is
only an interpretive ceiling for how much exploration could help if one could
choose the best basin after seeing truth.

## 33. Consolidated CSV outputs

The collection workbook exports consolidated files to:

```text
output/slurm_output/ensemble_scaling_full/<parent_job_id>/consolidated
```

The exported files include:

- `auprc_aggregated.csv`
- `auprc_individual.csv`
- `tpr05_aggregated.csv`
- `tpr05_individual.csv`
- `auprc_pooled_aggregated.csv`
- `auprc_pooled_agg_by_r2.csv`
- `auprc_pooled_individual_overall.csv`
- `auprc_pooled_individual.csv`
- `model_metrics_full.csv`
- `confusion_bins_full.csv`
- `dataset_metrics.csv`
- `multimodal_metrics.csv`
- `refine_depth.csv`
- `tier_cs_metrics_full.csv`
- `hg2_by_agg.csv`
- `scaling_bins_pooled.csv`
- `tpr05_pooled_aggregated.csv`
- `tpr05_pooled_agg_by_r2.csv`
- `tpr05_pooled_individual_overall.csv`
- `tpr05_pooled_individual.csv`

Some files are conditional on source data existing, but the paper-prep workflow
expects several of them to be present.

The scaling bins are pooled as:

```r
group_by(
  spec_name,
  annotation_r2,
  n_ensemble,
  resolution,
  agg_method,
  rep,
  pip_threshold
)
summarise(
  n_causal_at_bucket = sum(n_causal_at_bucket),
  n_noncausal_at_bucket = sum(n_noncausal_at_bucket),
  n_datasets = n_distinct(dataset_bundle_id)
)
```

## 34. Broad exploratory visualization workbook

The broad visualization workbook is:

- `test_susine/vignettes/simulation pipeline/visualize_results_workbook_ensemble_scaling.Rmd`

It uses:

```r
parent_job_id <- "51250228"
```

(this exploratory workbook still points at the older job; the paper-prep and
paper-rendering workbooks use `53547760`) and reads the consolidated CSVs for
that job.

It creates a wider set of exploratory plots, organized into categories:

1. SuSiE2-style plots:
   - TPR/FPR curves;
   - PIP calibration;
   - CS power by top-n causal definition;
   - CS FDR by top-n causal definition;
   - `hg2` boxplots.
2. Ensemble performance plots:
   - precision-recall curves;
   - multimodality metrics;
   - AUPRC bar charts;
   - TPR@FPR=0.05 bar charts;
   - C-CS grid landscape over `c` and `sigma_0_2`.
3. Dataset difficulty plots:
   - dataset metrics;
   - difficulty scatter;
   - model-based difficulty analyses;
   - failure maps;
   - refine-depth diagnostics.
4. Ensemble scaling plots:
   - AUPRC scaling curves;
   - TPR scaling curves;
   - two-lever interaction scaling;
   - compute vs AUPRC.
5. Bootstrap AUPRC stability.
6. Summary tables and metric inventory.

This workbook is useful for audit and exploration, but the newer paper-prep and
paper-rendering workbooks are the final paper-figure path.

## 35. Paper plot preparation

The paper-prep workbook is:

- `test_susine/vignettes/simulation pipeline/prepare_results_workbook_ensemble_scaling_paper.Rmd`

It uses:

```r
job_name <- "ensemble_scaling_full"
parent_job_id <- "53547760"
annotation_r2_focus <- 0.3
drop_aggs <- character(0)
bootstrap_B <- 2000L
bootstrap_seed <- 20260505L
true_h2 <- 0.25
```

It reads from:

```text
output/slurm_output/ensemble_scaling_full/53547760/consolidated
```

and writes plot data to:

```text
output/slurm_output/ensemble_scaling_full/53547760/figures/paper_ensemble_scaling/plot_data
```

The primary cached output is:

```text
paper_ensemble_plot_data.rds
```

The prep workbook loads required consolidated files:

- `auprc_pooled_aggregated.csv`
- `auprc_pooled_agg_by_r2.csv`
- `auprc_pooled_individual_overall.csv`
- `auprc_pooled_individual.csv`
- `confusion_bins_full.csv`
- `model_metrics_full.csv`

Optional files:

- `hg2_by_agg.csv`
- `scaling_bins_pooled.csv`

## 36. Paper-prep filtering and labels

The paper-prep workflow defines:

- `spec_order`:
  - `A-R`, `A-F`, `A-S`, `A-RF`, `A-RS`, `A-FS`, `A-RFS`,
  - `B-R`, `B-F`, `B-RF`,
  - `C-C`, `C-CS`, `C-CSR`, `C-CSFR`, `C-CFR`, `C-CF`;
- aggregation order:
  - `max_elbo`;
  - `uniform`;
  - `elbo_softmax`;
  - `cluster_weight_credible` (primary);
  - `cluster_weight_max`;
  - `oracle`.

It defines `is_focus_annotation_row()`:

- Group A rows are kept at `annotation_r2 = NA`;
- Groups B and C are kept at `annotation_r2 = 0.3` for focus plots.

The main paper plots therefore compare:

- annotation-agnostic Group A at its single annotation-agnostic level;
- annotation-consuming Groups B/C at `annotation_r2 = 0.3`;
- baseline SuSiE from `baseline-single / susine_vanilla`.

## 37. Paper-prep baseline and truth-warm references

The prep workbook extracts the baseline metric from:

- `spec_name == "baseline-single"`;
- `use_case_id == "susine_vanilla"`;
- `annotation_r2 = NA`;
- individual pooled metric tables.

This baseline is labelled as the SuSiE baseline in paper plots, even though the
backend is `susine_vanilla` through the local `susine` package.

The truth-warm SuSiE reference is extracted from:

- `spec_name == "truth-warm"`;
- `use_case_id == "susine_vanilla"`;
- `annotation_r2 = NA`.

This reference is used in at least the C-CS-by-annotation plot as a horizontal
truth-warm ceiling.

## 38. Paper-prep heatmap data

`heatmap_data` is built from the shared pooled AUPRC plot table:

- keep focus annotation rows;
- `drop_aggs` is now `character(0)`, so no aggregation methods are dropped;
- keep aggregations in `agg_order`;
- keep specs in `spec_order`;
- compute:

```r
delta_auprc = AUPRC - baseline_auprc
delta_pct = 100 * delta_auprc / baseline_auprc
```

This feeds:

- `paper_delta_auprc_heatmap_r2_0p2.png`

The visualizer highlights the `C-CS / cluster_weight_credible` tile with a red
border.

## 39. Paper-prep C-CS by annotation data

`ccs_by_r2` focuses on:

- `spec_name == "C-CS"`;
- annotation levels 0, 0.3, 0.5;
- `agg_method` in the main aggregation order (`drop_aggs` is now `character(0)`,
  so nothing is dropped).

The paper-rendering workbook further keeps:

- `cluster_weight_credible`;
- `oracle`.

This feeds:

- `paper_ccs_aggregation_by_annotation_r2.png`

The plot includes horizontal reference lines for:

- baseline SuSiE;
- truth-warm SuSiE.

## 40. Paper-prep bootstrap

The bootstrap compares:

- baseline:
  - `baseline-single / susine_vanilla`;
  - individual confusion bins;
- ensemble:
  - `C-CS`;
  - `annotation_r2 = 0.3`;
  - `agg_method = "cluster_weight_credible"`;
  - aggregated confusion bins.

The bootstrap unit is `dataset_bundle_id`.

Procedure:

1. Find shared dataset IDs present in both baseline and C-CS bins.
2. For each of `bootstrap_B = 2000` replicates:
   - sample dataset IDs with replacement;
   - weight each dataset's bins by how many times it was sampled;
   - pool weighted bins;
   - compute AUPRC for baseline;
   - compute AUPRC for C-CS cluster_weight_credible;
   - store the delta.
3. Compute the observed paired delta using the unresampled shared dataset set.
4. Compute the 2.5 percent and 97.5 percent bootstrap quantiles.

This feeds:

- `paper_bootstrap_auprc_r2_0p2.png`

## 41. Paper-prep precision-recall curves

The PR data are built by `build_pr_curve_counts()` from confusion bins.

Rows included:

- baseline:
  - `baseline-single / susine_vanilla`;
  - individual bins;
  - crossed with all annotation r2 values for plotting alongside C-CS;
- C-CS:
  - `spec_name == "C-CS"`;
  - `agg_method == "cluster_weight_credible"`;
  - annotation r2 values 0, 0.3, 0.5.

The prep workbook also computes two operating-point diagnostics:

- `paper_pr_precision50_operating_points.csv`:
  - maximum recall with precision >= 0.50;
- `paper_pr_same_recall_as_baseline_precision50.csv`:
  - non-baseline rows closest to the baseline's 50-percent-precision recall.

This feeds:

- `paper_pr_curves_ccs_cluster_weight_by_r2.png`

## 42. Paper-prep `hg2` data

If `hg2_by_agg.csv` is available, the prep workbook:

- normalizes aggregation-method IDs;
- enriches metadata from `model_metrics`;
- keeps focus annotation rows;
- keeps `agg_method == "cluster_weight_credible"`;
- labels `C-CS` specially;
- adds baseline `hg2` rows from `model_metrics` for
  `baseline-single / susine_vanilla`.

The plotted reference heritability is:

```r
true_h2 <- 0.25
```

This feeds:

- `paper_hg2_cluster_weight_r2_0p2.png`

## 43. Paper-prep PIP calibration data

Calibration compares:

- baseline:
  - `baseline-single / susine_vanilla`;
  - individual bins;
- ensemble:
  - `C-CS`;
  - `annotation_r2 = 0.3`;
  - `agg_method = "cluster_weight_credible"`.

The prep workflow bins by PIP threshold into width-0.10 bins:

```r
pip_bin = floor(pip_threshold / 0.10) * 0.10
```

with a cap at 0.9.

For each method, dataset, and PIP bin, it computes observed causal fraction from
bucket counts. It then summarizes across datasets with standard errors and
approximate 95-percent intervals.

This feeds:

- `paper_pip_calibration_ccs_cluster_weight_r2_0p2.png`

## 44. Paper-prep scaling data

The prep workbook uses `scaling_bins_pooled.csv` if available.

Function:

- `compute_scaling_metrics()`

It:

- applies `drop_aggs` (now `character(0)`, so no aggregation methods are dropped);
- splits pure scaling rows from interaction rows:
  - pure rows have `resolution = NA`;
  - interaction rows have non-missing `resolution`;
- for each group and aggregation method, computes AUPRC from pooled scaling
  bins;
- summarizes mean and SE over `rep`.

Because runtime scaling emits `rep = 1` in the current implementation, most SEs
from runtime scaling are expected to be missing or zero-equivalent.

### 44.1 C-CS 4x4/8x8 cluster-weight and oracle

The prep workbook constructs `ccs_resolution_plot_data` by combining:

1. `ccs_resolution_cluster`:
   - from `scaling_metrics`;
   - `spec_name == "C-CS"`;
   - `annotation_r2 = 0.3`;
   - `agg_method == "cluster_weight_credible"`;
   - resolutions `4x4` and `8x8`;
2. `compute_ccs_resolution_oracle()`:
   - recomputes an oracle from individual C-CS confusion bins;
   - for each dataset and resolution, selects evenly spaced c and sigma values;
   - computes per-run AUPRC;
   - picks the best individual run per dataset;
   - pools the selected bins and computes AUPRC.

This feeds:

- `paper_ccs_4x4_8x8_r2_0p2.png`

The oracle here is a subgrid oracle, not the same as the full collection oracle
over every individual fit in the full C-CS grid.

### 44.2 Compute vs AUPRC

The prep workbook estimates compute for scaling points using:

- `prepare_individual_run_rows()`;
- `select_subscale_run_rows()`;
- `estimate_subscale_compute()`;
- `baseline_compute_minutes()`.

It uses `model_metrics$wall_time_sec` from purity-filtered individual rows,
selects the individual run rows corresponding to the subscale point, sums wall
time within dataset, and averages across datasets.

The plotted compute metric is:

```r
compute_min = mean_total_wall_sec / 60
```

For baseline, it uses the mean wall time of:

- `baseline-single / susine_vanilla`;
- purity-filtered rows.

This feeds:

- `paper_compute_vs_auprc_cluster_weight_r2_0p2.png`

## 45. Paper plot rendering

The final rendering workbook is:

- `test_susine/vignettes/simulation pipeline/visualize_results_workbook_ensemble_scaling_paper.Rmd`

It uses:

```r
job_name <- "ensemble_scaling_full"
parent_job_id <- "53547760"
```

and reads:

```text
figures/paper_ensemble_scaling/plot_data/paper_ensemble_plot_data.rds
```

It saves figures with `save_paper_png()`, using `ragg` if available and Cairo
PNG otherwise.

The rendering workbook's expected output figures are:

- `paper_delta_auprc_heatmap_r2_0p2.png`
- `paper_ccs_aggregation_by_annotation_r2.png`
- `paper_bootstrap_auprc_r2_0p2.png`
- `paper_pr_curves_ccs_cluster_weight_by_r2.png`
- `paper_hg2_cluster_weight_r2_0p2.png`
- `paper_hg2_decomposition_r2_0p2.png`
- `paper_hg2_truth_error_r2_0p2.png` when truth-error data are available
- `paper_pip_calibration_ccs_cluster_weight_r2_0p2.png`
- `paper_ccs_4x4_8x8_r2_0p2.png`
- `paper_compute_vs_auprc_cluster_weight_r2_0p2.png`
- `paper_ensemble_scaling_composite_heatmap.png`
- `paper_ensemble_scaling_composite_performance.png`
- `paper_ensemble_multimodality.png`

Those are the workflow-level figure names expected under the job output tree.
The current local manuscript/writing folder contains copied composite and
supplement PNGs:

```text
../Writings/plots/ensemble_sims/paper_ensemble_scaling_composite_heatmap.png
../Writings/plots/ensemble_sims/paper_ensemble_scaling_composite_performance.png
../Writings/plots/ensemble_sims/paper_ensemble_multimodality.png
```

These three files are the current manuscript-facing ensemble simulation
artifacts copied under `../Writings/plots`.

Observed local figure metadata:

| file | dimensions | size | modified time |
|---|---:|---:|---|
| `paper_ensemble_multimodality.png` | 5400 x 2100 | 257,785 | 2026-06-15 11:17:03 |
| `paper_ensemble_scaling_composite_heatmap.png` | 5400 x 3060 | 812,809 | 2026-06-15 11:17:05 |
| `paper_ensemble_scaling_composite_performance.png` | 5940 x 6210 | 1,100,957 | 2026-06-15 11:17:07 |

The local `paper_ensemble_scaling_composite_heatmap.png` contains the heatmap
version of the main ensemble result, centered on:

- A: delta-AUPRC heatmap over ensemble spec and aggregation method;
- selected supporting performance panels from the same prepared RDS, depending
  on which components were available during rendering.

The local `paper_ensemble_scaling_composite_performance.png` contains the
performance subset:

- B: C-CS AUPRC by annotation accuracy;
- C: PIP calibration;
- D: precision-recall curves by annotation accuracy;
- F: estimated heritability;
- G: AUPRC by compute cost.

The local `paper_ensemble_multimodality.png` is a supplementary multimodality
summary rendered from the same paper workflow. There is still no local
individual-panel PNG folder under `../Writings/plots/ensemble_sims`; the
code-level individual panels are left under the HPC/job-output figure tree or
not copied into the writing directory.

The local composite labels the annotation axis as `phi_a`, while the code field
feeding these plots is `annotation_r2`. In methods text, these should be tied
together explicitly: the plotted `phi_a` is the simulator's target causal
annotation squared-correlation parameter as stored in code as `annotation_r2`.

Visible values in panel A of the local composite provide a useful sanity check
against the code trace:

- the highlighted primary method, `C-CS` with `cluster_weight_credible`, has
  `Delta AUPRC = +0.059 (+25%)` at the focus annotation setting;
- `C-CS` with `oracle` has `Delta AUPRC = +0.090 (+38%)`;
- among non-oracle aggregators, C-CS cluster-weight is the largest visible gain
  in the heatmap;
- the baseline AUPRC implied by `+0.059 (+25%)` is about 0.236, consistent with
  the dashed SuSiE baseline line in the compute-cost panel.

## 46. Figure-by-figure data provenance

### 46.1 Delta-AUPRC heatmap

File:

- `paper_delta_auprc_heatmap_r2_0p2.png`

Data:

- `heatmap_data` from paper-prep RDS.

Source metric:

- pooled AUPRC from consolidated pooled AUPRC tables.

Rows:

- specs in `spec_order`;
- Group A at `annotation_r2 = NA`;
- Groups B/C at `annotation_r2 = 0.3`;
- aggregation methods:
  - max ELBO;
  - uniform;
  - ELBO softmax;
  - cluster_weight_credible (primary);
  - cluster_weight_max;
  - oracle.

Displayed value:

- AUPRC delta relative to `baseline-single / susine_vanilla`.

### 46.2 C-CS aggregation by annotation r2

File:

- `paper_ccs_aggregation_by_annotation_r2.png`

Data:

- `ccs_by_r2`.

Rows:

- `spec_name == "C-CS"`;
- annotation r2 values 0, 0.3, 0.5;
- plotted methods:
  - cluster_weight_credible;
  - oracle.

References:

- baseline SuSiE horizontal line;
- truth-warm SuSiE horizontal line.

### 46.3 Bootstrap delta AUPRC

File:

- `paper_bootstrap_auprc_r2_0p2.png`

Data:

- `boot_plot_data`.

Comparison:

- C-CS cluster_weight_credible at `annotation_r2 = 0.3`;
- minus baseline SuSiE.

Bootstrap:

- 2000 paired dataset resamples;
- dataset IDs sampled with replacement.

### 46.4 PR curves

File:

- `paper_pr_curves_ccs_cluster_weight_by_r2.png`

Data:

- `pr_curve_counts`.

Methods:

- SuSiE baseline;
- SuSiNE C-CS cluster-weight.

Annotation r2 facets/series:

- 0;
- 0.2;
- 0.5.

### 46.5 `hg2` plot

File:

- `paper_hg2_cluster_weight_r2_0p2.png`

Data:

- `hg2_plot_data`.

Methods:

- baseline;
- cluster-weight ensemble series;
- C-CS highlighted.

Reference:

- true `h2 = 0.25`.

### 46.6 PIP calibration

File:

- `paper_pip_calibration_ccs_cluster_weight_r2_0p2.png`

Data:

- `pip_cal_with_se`.

Methods:

- SuSiE baseline;
- SuSiNE C-CS cluster-weight at `annotation_r2 = 0.3`.

### 46.7 C-CS 4x4/8x8 resolution

File:

- `paper_ccs_4x4_8x8_r2_0p2.png`

Data:

- `ccs_resolution_plot_data`.

Series:

- cluster-weight aggregation;
- subgrid oracle.

Resolutions:

- 4x4;
- 8x8.

Reference:

- baseline AUPRC.

### 46.8 Compute vs AUPRC

File:

- `paper_compute_vs_auprc_cluster_weight_r2_0p2.png`

Data:

- `compute_plot_data`.

Rows:

- cluster-weight scaling points;
- focus annotation rows;
- specs in `spec_order`.

X-axis:

- mean compute minutes per dataset, estimated from individual
  `wall_time_sec`.

Y-axis:

- pooled AUPRC from scaling bins.

Reference:

- baseline compute/AUPRC point or line.

## 47. What is scientifically varied?

The ensemble study varies four main exploration levers:

1. Random restarts.
   - Perturb the initial alpha matrix.
   - Meant to test whether optimizer basins are sensitive to initialization.

2. Refine.
   - Iteratively block high-purity credible sets, perturb the fit, then refit
     without the block from the perturbed alpha.
   - Meant to force the algorithm to explore alternative basins or secondary
     causal configurations.

3. Prior-mean scale `c`.
   - Vary the strength of the signed annotation vector in `mu_0 = c * a`.
   - Meant to test annotation-strength calibration and basin movement under
     different annotation influence.

4. Prior variance `sigma_0_2`.
   - Vary fixed prior-effect variance.
   - Meant to test sensitivity to prior shrinkage and effect-size scale.

The specs test these levers singly and in combinations under:

- no annotations (Group A);
- EB annotation scale/variance learning (Group B);
- fixed grid-based annotation scale/variance sweeps (Group C).

The central paper emphasis appears to be C-CS:

- `c_grid x sigma_0_2_grid`;
- 8 x 8 grid;
- functional-mu SuSiNE;
- cluster-weight aggregation;
- focus annotation r2 0.3.

This is the spec highlighted in the heatmap and carried through the PR,
bootstrap, calibration, h2, resolution, and compute figures.

## 48. Important implementation notes for methods writing

### 48.1 "SuSiE baseline" in the ensemble figures

The paper plots use `baseline-single / susine_vanilla` as the SuSiE baseline.
This uses the `susine` backend with zero prior mean, fixed prior variance, and
uniform inclusion weights. If the manuscript calls this "SuSiE", it may be worth
phrasing as "SuSiE-equivalent zero-prior-mean fit implemented through the
SuSiNE code path" unless exact agreement with `susieR` is established elsewhere.

### 48.2 Annotation quality parameterization differs from the baseline draft

The manuscript baseline section discusses a 3 x 3 grid in `(phi_a, nu_a)` and
an AlphaGenome-calibrated regime around `(0.3, 0.9)` or `(0.3, 1.0)`.

The ensemble run uses:

- `annotation_r2` in `{0, 0.3, 0.5}`;
- `inflate_match = 0.95`.

So the ensemble study is a one-dimensional annotation-accuracy sweep at fixed
noncausal annotation scale, not the full two-dimensional baseline annotation
grid.

### 48.3 The paper focus is r2 = 0.3

The paper-prep workflow sets:

```r
annotation_r2_focus <- 0.3
```

Most main figures compare Group B/C at r2 0.3 against Group A and baseline
annotation-agnostic fits.

### 48.4 Top-8 causal scoring should be stated clearly

The confusion-bin AUPRC used for paper curves treats only the top 8 causal SNPs
by absolute effect size as positives. This is a major methods detail.

The data-generating model has 23 causal SNPs, but AUPRC/PR paper plots are not
asking the methods to recover all 23 equally. They score recovery of the largest
8 causal effects. The remaining (weak) causal SNPs are DROPPED from the scoring
entirely -- they are neither positives nor negatives
(`keep <- !(causal == 1L & scored == 0L)` in `compute_confusion_bins()`) -- so
the curves are not penalized for failing to recover them.

### 48.5 Purity thresholds differ by purpose

There are two purity thresholds:

- 0.50 for output model summaries and purity-filtered credible-set metrics;
- 0.95 for refine branching eligibility.

Do not collapse these into one threshold in the methods text.

### 48.6 Refine nominal K is not exact compute K

Refine rows store one ensemble member per row, but non-root refine steps cost two
fits:

- perturb with blocked weights;
- refit with restored weights and perturb alpha warm start.

Therefore nominal `K` in refine specs is a planning label, not always exact
model-call count. Compute plots use observed `wall_time_sec`, which is better
than nominal K for cost.

### 48.7 Scaling curves are deterministic subsets

Although `scaling_restart_reps` is set to 50 in the run-control workbook, the
runtime scaling helper emits one deterministic subset per scaling size. Restart
and refine subsets use the first IDs/steps, while c/sigma subsets use evenly
spaced grid values.

### 48.8 Cluster-weight aggregation now agrees between full and scaling paths

For the scheduled `cluster_weight_credible` and `cluster_weight_max` methods,
both the full-ensemble path and the scaling-bin path route through the same
Method B logic (each method's own metric + threshold, frequency-free cluster
weights split by within-cluster ELBO-softmax, all fits contributing). The older
full-vs-scaling discrepancy -- full path used a max-ELBO representative per JSD
cluster with a 1/frequency correction, while the scaling helper weighted clusters
uniformly with no representative -- is obsolete for the scheduled methods; both
paths agree at full resolution.

That historical discrepancy was defined only for the legacy JSD-0.15 +
frequency `cluster_weight` method (Method C). Current paper workflows use
metric-specific `cluster_weight_credible` / `cluster_weight_max` Method-B
weights instead.

### 48.9 Truth-warm is an upper-bound diagnostic, not a method

Truth-warm rows initialize effects on known causal variants. They should be
described only as a diagnostic ceiling or sanity check.

### 48.10 `C-CSR` asks a different question from C-CS

`C-CSR` is not simply C-CS plus aggregation. It uses functional-mu grid fits to
find basins, then refits with default zero-prior settings from exact converged
effect fits. It tests basin durability after removing the annotation prior.

### 48.11 The collection workbook parent ID should be updated before rerun

If rerunning collection for the paper figures, update:

```r
parent_job_id <- "53520405"
```

in `collect_results_workbook_ensemble_scaling.Rmd`, or parameterize it. The
paper-prep and paper-rendering workbooks already use `53547760`.

## 49. Methods-section distillation

A concise methods version could say:

We evaluated ensemble fine-mapping on the same 600 oligogenic simulated loci used
in the baseline simulation screen. For each of 150 genotype matrices, we
generated four independent phenotypes under the SuSiE 2.0 oligogenic
architecture with 23 causal SNPs partitioned into 3 sparse, 5 oligogenic, and 15
polygenic effects, rescaled to explain 50, 35, and 15 percent of effect energy,
respectively, with total target heritability 0.25. All fits used `L = 10` and at
most 100 IBSS iterations.

We generated signed functional annotations by mixing the true causal effects with
orthogonal Gaussian noise on causal SNPs to achieve target squared correlation
`annotation_r2`, and by drawing centered noncausal annotations with RMS matched
to 95% of the causal annotation centered second moment (`inflate_match = 0.95`).
Ensemble simulations used `annotation_r2` values 0, 0.3, and 0.5; primary paper
plots focused on 0.3.

We compared three model families: an annotation-agnostic zero-prior-mean
SuSiE-equivalent fit through the SuSiNE backend, an empirical-Bayes SuSiNE model
that estimates nonnegative annotation scale and prior variance, and a fixed-
variance SuSiNE functional-mu model with `mu_0 = c * a`. We constructed
ensembles by varying random initializations, refine perturbations, prior-mean
scale `c`, and prior variance `sigma_0_2`, both individually and in selected
combinations. The main fixed-grid SuSiNE ensemble (`C-CS`) crossed 8 values of
`c` between 0 and 1.5 with 8 values of `sigma_0_2` between 0.01 and 1.0.

Within each dataset and ensemble group, we aggregated PIPs using several rules:
selecting the max-ELBO fit, uniformly averaging PIPs, ELBO-softmax averaging,
and two cluster-weighted rules. The primary cluster-weight rule
(`cluster_weight_credible`) clustered PIP vectors by complete linkage on a
single-variant credibility-weighted shift metric,
`credible_shift = max_j[ max(p_j,q_j) * |p_j - q_j| ]`, cut at 0.05; a companion
rule (`cluster_weight_max`) used the unweighted L-infinity shift
`max_j |p_j - q_j|` cut at 0.10. Both aggregated with a frequency-free scheme
(Method B): each cluster received weight proportional to exp(max-ELBO-in-cluster),
split within the cluster across all member fits by ELBO-softmax, so every fit
contributes and no inverse-frequency correction is applied.

Performance curves were computed from PIP confusion bins pooled across datasets.
For classification-style AUPRC/PR analyses, positives were defined as the top 8
true causal SNPs per dataset by absolute effect size; the remaining (weak) causal
SNPs were excluded from scoring entirely (neither positive nor negative). The
oracle aggregator was a truth-aware diagnostic that selected the best individual
ensemble member per dataset by AUPRC before pooling its confusion bins.

## 50. Key code reference index

Run configuration:

- `test_susine/vignettes/simulation pipeline/run_control_workbook_ensemble_scaling.Rmd`
- `test_susine/R/run_controls.R:60` - `build_data_matrix_catalog`
- `test_susine/R/run_controls.R:146` - `make_run_tables`
- `test_susine/R/run_controls.R:750` - `make_job_config`
- `test_susine/R/use_cases.R:9` - `prior_spec_catalog`
- `test_susine/R/use_cases.R:65` - `exploration_catalog`
- `test_susine/R/use_cases.R:82` - `aggregation_catalog`
- `test_susine/R/use_cases.R:99` - `valid_exploration_for_prior`

Simulation:

- `test_susine/R/simulate_data.R:178` - `simulate_effect_sizes_susie2_oligogenic`
- `test_susine/R/simulate_data.R:257` - `simulate_phenotype_h2`
- `test_susine/R/simulate_data.R:356` - `simulate_priors`
- `test_susine/R/run_model.R:525` - `generate_data_for_bundle`

Execution:

- `test_susine/R/run_model.R:13` - `run_task`
- `test_susine/R/run_model.R:622` - `execute_dataset_bundle`
- `test_susine/R/run_model.R:333` - `execution_cache_key`
- `test_susine/R/run_model.R:1587` - `run_default_prior_refit`
- `test_susine/R/run_model.R:1623` - `run_use_case`
- `test_susine/R/run_model.R:1931` - `aggregate_use_case_pips`
- `test_susine/R/run_model.R:2029` - `compute_multimodal_metrics`
- `test_susine/R/run_model.R:2165` - `write_run_outputs`
- `test_susine/R/run_model.R:2454` - `compute_scaling_confusion_bins_for_group`
- `test_susine/R/run_model.R:2086` - `.cluster_weights_from_hc` (Method B)
- `test_susine/R/run_model.R:2143` - `.cluster_weight_specs` (locked-in cluster-weight metrics/thresholds)
- `test_susine/R/run_model.R:2850` - `compute_hg2_by_agg`
- `test_susine/R/heritability.R` - fair-PVE decomposition engine (`hg2_uncertainty_scalar`, etc.)

Metrics:

- `test_susine/R/evaluation_metrics.R:315` - `evaluate_model`
- `test_susine/R/run_model.R:218` - `compute_confusion_bins`
- `test_susine/R/collect_results.R:809` - `auprc_from_pooled_bins`
- `test_susine/R/collect_results.R:830` - `aggregate_pip_matrix`
- `test_susine/R/visualize_results.R:131` - `compute_auprc_from_pooled_confusion`
- `test_susine/R/visualize_results.R:337` - `compute_tpr05_from_pooled_confusion`

Collection and plotting:

- `test_susine/vignettes/simulation pipeline/collect_results_workbook_ensemble_scaling.Rmd`
- `test_susine/R/ensemble_scaling_utils.R:25` - `append_overall_scaling_bins`
- `test_susine/R/ensemble_scaling_utils.R:72` - `backfill_pure_refine_scaling_bins`
- `test_susine/R/ensemble_scaling_utils.R:307` - `build_shared_plot_metric_table`
- `test_susine/R/ensemble_scaling_utils.R:423` - `backfill_single_fit_refine_agg_confusion`
- `test_susine/R/ensemble_scaling_utils.R:540` - `backfill_terminal_scaling_agg_confusion`
- `test_susine/R/ensemble_scaling_utils.R:700` - `backfill_missing_agg_confusion_from_snps`
- `test_susine/vignettes/simulation pipeline/prepare_results_workbook_ensemble_scaling_paper.Rmd`
- `test_susine/vignettes/simulation pipeline/visualize_results_workbook_ensemble_scaling_paper.Rmd`
