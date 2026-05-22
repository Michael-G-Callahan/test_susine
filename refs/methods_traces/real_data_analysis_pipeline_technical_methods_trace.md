# Real Data Analysis Pipeline Technical Trace

Date written: 2026-05-22

This document traces the real-data eQTL case-study pipeline that produced the
paper figures:

- `../Writings/plots/real_data_case_study/paper_overall_annotation_drift_pve.png`
- `../Writings/plots/real_data_case_study/paper_selected_locus_zoom_lollipop_l2_drift.png`

The trace is intentionally implementation-facing. Its goal is to make the full
data provenance and analysis path recoverable from source files, even though the
main generated real-data inputs and SLURM outputs were produced on HPC and are
not present in this local checkout.

## 1. Scope, provenance, and caveats

Main manuscript context:

- `../Writings/second draft/susine_second_draft.tex`

Main downstream fitting and plotting workbooks:

- `test_susine/vignettes/real data pipeline/prepare_real_data_inputs_workbook.Rmd`
- `test_susine/vignettes/real data pipeline/run_control_workbook_real_data_ensemble.Rmd`
- `test_susine/vignettes/real data pipeline/collect_results_workbook_real_data_ensemble.Rmd`
- `test_susine/vignettes/real data pipeline/visualize_results_workbook_real_data_ensemble.Rmd`

Main downstream helper code:

- `test_susine/R/real_data_pipeline.R`
- `test_susine/inst/scripts/run_real_data_task.R`
- `test_susine/R/run_model.R`
- `test_susine/R/collect_results.R`

Main upstream annotation/data-prep repo:

- `eQTL_annotations_for_susine/`
- `eQTL_annotations_for_susine/real_data_geometric_n20_workflow.md`
- `eQTL_annotations_for_susine/scripts/prepare_geometric_n20_mu0.py`
- `eQTL_annotations_for_susine/scripts/run_phase1_task.py`
- `eQTL_annotations_for_susine/scripts/run_annotation_batch.py`
- `eQTL_annotations_for_susine/utils/pipeline_zscore.py`
- `eQTL_annotations_for_susine/utils/pipeline_ld.py`
- `eQTL_annotations_for_susine/utils/pipeline_alphagenome_batch.py`

Important local-state caveats at trace time:

- The local checkout does not contain
  `test_susine/data/real_case_studies/geometric_n20_loci/`.
- The local checkout does not contain
  `test_susine/output/slurm_output/real_data_ensemble_geometric_n20/52906940/`.
- Therefore, the final figures could not be regenerated locally from the
  currently synced artifacts.
- The plotting workbook uses `parent_job_id <- "52906940"` and saves its paper
  figures under:
  `test_susine/output/slurm_output/real_data_ensemble_geometric_n20/52906940/figures/real_data_ensemble/overall/`.
- The paper copies currently live under:
  `../Writings/plots/real_data_case_study/`.
- The final PNGs in the paper directory have local write times on
  2026-05-22, while the plotting workbook was last written on 2026-05-21. This
  is consistent with rendering or copying the final figures after editing the
  visualization workbook.
- The upstream annotation-prep workbook/code references the source repo as
  `/storage/work/mgc5166/Annotations/eQTL_annotations_for_susine`, while the
  user-provided HPC shell example points to
  `/storage/home/mgc5166/work/Annotations/eQTL_annotations_for_susine`. These
  may be equivalent via a filesystem convention or symlink, but that is not
  verifiable locally.
- The manuscript still has placeholder real-data methods/results text. The
  pipeline traced here is more specific than the current draft.

## 2. Scientific purpose

The real-data study applies the SuSiNE RSS workflow to GTEx Lung eQTL summary
statistics at selected protein-coding gene loci. It asks whether signed
AlphaGenome annotations, inserted through the prior-mean channel
`mu_0 = c * a`, change fine-mapping outputs relative to annotation-agnostic
SuSiE/SuSiNE RSS fits.

The final figures focus on four linked diagnostics:

1. how much weight the ensemble assigns to annotation-influenced fits;
2. how much the highest-weight annotation-informed model and its warm baseline
   refit drift from a SuSiE baseline in PIP space and effect-basis space;
3. whether the annotation-informed ensemble or warm refit changes the
   standardized posterior-mean PVE proxy;
4. whether annotation-informed fits sharpen per-effect posterior support.

The selected-locus zoom figure then drills into three loci:

- `rrp7a` -> expected full locus id `rrp7a_chr22_lung`
- `arsa` -> expected full locus id `arsa_chr22_lung`
- `prpf38b` -> expected full locus id `prpf38b_chr1_lung`

The workbook resolves these by matching the tokens against the available
`locus_id` values.

## 3. Upstream data sources

### 3.1 GTEx Lung eQTL summary statistics

The upstream repo stores GTEx inputs under:

```text
eQTL_annotations_for_susine/data/gtex/
```

Local files visible at trace time include:

```text
GTEx_Analysis_v10_QTLs_GTEx_Analysis_v10_eQTL_all_associations_Lung.v10.allpairs.chr1.parquet
GTEx_Analysis_v10_QTLs_GTEx_Analysis_v10_eQTL_all_associations_Lung.v10.allpairs.chr5.parquet
GTEx_Analysis_v10_QTLs_GTEx_Analysis_v10_eQTL_all_associations_Lung.v10.allpairs.chr12.parquet
GTEx_Analysis_v10_QTLs_GTEx_Analysis_v10_eQTL_all_associations_Lung.v10.allpairs.chr17.parquet
GTEx_Analysis_v10_QTLs_GTEx_Analysis_v10_eQTL_all_associations_Lung.v10.allpairs.chr22.parquet
Lung.v10.eQTLs.SuSiE_summary.csv
Lung.v10.eQTLs.SuSiE_summary.parquet
```

The z-score export reads the chromosome-specific GTEx Lung allpairs parquet for
the target gene. The source fields used downstream are:

- `variant_id`
- `slope`
- `slope_se`
- `af`
- `ma_count`
- `gene_id`

The effective sample size is computed as:

```text
sample_size = ma_count / (2 * af)
```

The marginal association statistic is computed as:

```text
z_score = slope / slope_se
```

Default filters from the manifest are:

```text
0.01 <= af <= 0.99
sample_size > 50
```

The z-score stage writes, per locus:

```text
output/z_score/<locus_id>/<GENE>_GTEx_z_scores.csv
output/z_score/<locus_id>/<GENE>_variants.vcf
output/prelim/<locus_id>/<GENE>_phase1_count_funnel.csv
```

The implementation is in:

- `eQTL_annotations_for_susine/utils/pipeline_zscore.py`
- prototype notebook: `eQTL_annotations_for_susine/vignettes/1_get_z_scores.ipynb`

### 3.2 Gene models

Gene metadata are read from GENCODE annotations cached under:

```text
eQTL_annotations_for_susine/data/gtf_cache/
```

If missing, the code downloads GENCODE v44 for hg38 from EBI. Gene shortcuts are
stored under:

```text
eQTL_annotations_for_susine/data/gtf_shortcuts/
```

The z-score and AlphaGenome stages use the gene coordinates, exons, and TSS to
evaluate sequence-window eligibility.

### 3.3 LD reference panel

LD is computed from the 1000 Genomes Project GRCh38 reference panel.

The sample panel is:

```text
integrated_call_samples_v3.20130502.ALL.panel
```

The code keeps EUR samples only, where the third field in the panel is `EUR`.

The genotype source is the EBI 1000 Genomes GRCh38 biallelic SNV/INDEL release:

```text
https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/
1000_genomes_project/release/20190312_biallelic_SNV_and_INDEL/
ALL.chr<chrom>.shapeit2_integrated_snvindels_v2a_27022019.GRCh38.phased.vcf.gz
```

Variants exported from GTEx are matched to 1000 Genomes by:

- chromosome,
- position,
- reference allele,
- alternate allele.

Dosage genotypes for EUR samples are mean-imputed for missing genotypes,
centered, scaled, and used to compute:

```text
R = cor(G_EUR)
```

The implementation rounds the LD matrix to three decimals and stores the upper
triangle in long-form parquet.

The LD stage writes, per locus:

```text
output/ld/<locus_id>/<GENE>_LD_variant_order.tsv
output/ld/<locus_id>/<GENE>_phase1_master_variants.csv
output/ld/<locus_id>/<GENE>_phase1_variant_map.parquet
output/ld/<locus_id>/<GENE>_phase1_z_scores.parquet
output/ld/<locus_id>/<GENE>_phase1_LD_R_long.parquet
output/prelim/<locus_id>/<GENE>_phase1_dataset_metrics.csv
```

The implementation is in:

- `eQTL_annotations_for_susine/utils/pipeline_ld.py`
- prototype notebook: `eQTL_annotations_for_susine/vignettes/2_cbx8_1000g_ld.ipynb`

The canonical downstream variant order is always:

```text
output/ld/<locus_id>/<GENE>_LD_variant_order.tsv
```

The user-provided HPC example for `abca10_chr17_lung` matches this layout:

```text
/storage/home/mgc5166/work/Annotations/eQTL_annotations_for_susine/output/ld/abca10_chr17_lung/
  ABCA10_LD_variant_order.tsv
  ABCA10_phase1_LD_R_long.parquet
  ABCA10_phase1_master_variants.csv
  ABCA10_phase1_variant_map.parquet
  ABCA10_phase1_z_scores.parquet
```

## 4. Locus selection

The current production note describes the selected panel as the geometric N20
real-data workflow. The relevant source note is:

```text
eQTL_annotations_for_susine/real_data_geometric_n20_workflow.md
```

Candidate loci came from GTEx v10 Lung allpairs parquet files. The main
candidate manifest for the geometric workflow is:

```text
eQTL_annotations_for_susine/config/loci_manifest_sample_100_per_chrom.csv
```

The manifest schema is documented in:

```text
eQTL_annotations_for_susine/config/manifest_schema.md
```

The locked sample is described as a predicted-vs-actual geometric maximin
sample in the plane:

```text
x = log1p(predicted pip_k_eff from dataset metrics)
y = log1p(actual baseline SuSiE pip_k_eff)
```

The official geometric sample artifacts are expected at:

```text
output/baseline_susie_screen/selected_geometric_n20/locked_loci.csv
output/baseline_susie_screen/selected_geometric_n20/locked_loci_with_manifest.csv
output/baseline_susie_screen/selected_geometric_n20/annotation_selection_geometric_n20.csv
config/annotation_selection_geometric_n20.csv
```

The locking script is:

```text
eQTL_annotations_for_susine/scripts/lock_geometric_n20_sample.py
```

At trace time, `config/annotation_selection_geometric_n20.csv` was not present
locally, but `loci_manifest_sample_100_per_chrom.csv` did contain the selected
zoom-locus candidates:

- `prpf38b_chr1_lung`
- `arsa_chr22_lung`
- `rrp7a_chr22_lung`

Open provenance issue: the manuscript draft says the real-data study uses
10 protein-coding loci, while the current source workflow and job names are
`geometric_n20`. The final paper figures appear to come from the geometric N20
pipeline unless a later filter reduced the displayed or successful loci.

## 5. AlphaGenome annotation generation

AlphaGenome scoring is run from the upstream repo with:

```text
eQTL_annotations_for_susine/scripts/run_annotation_batch.py
eQTL_annotations_for_susine/utils/pipeline_alphagenome_batch.py
eQTL_annotations_for_susine/utils/alphagenome.py
```

Expected command for the geometric N20 production run:

```bash
python scripts/run_annotation_batch.py \
  --manifest config/loci_manifest_sample_100_per_chrom.csv \
  --selection config/annotation_selection_geometric_n20.csv \
  --summary-path output/annotation/alphagenome/geometric_n20_annotation_batch_summary.csv
```

The API key is read from:

```text
ALPHAGENOME_API_KEY
```

unless a different environment variable is supplied with `--api-key-env-var`.

For each selected locus, the annotation runner reads:

```text
output/ld/<locus_id>/<GENE>_phase1_master_variants.csv
```

and scores variants for the target GTEx/Lung RNA-seq context. The manifest
default sequence length is:

```text
1MB
```

The interval policy is:

- `tss_centered`: preferred when both alleles fit inside the canonical 1 Mb
  TSS-centered interval.
- `midpoint_centered`: fallback when the target-gene exon span plus both alleles
  fits inside a 1 Mb interval but the TSS-centered interval is not usable.
- `ineligible`: variants that fit neither interval are excluded before API
  submission and logged.

The scorer uses:

```python
dna_client.OutputType.RNA_SEQ
variant_scorers.GeneMaskLFCScorer(requested_output=target_output)
```

The runner performs a preflight request to confirm the target gene and target
GTEx/Lung rows are present in AlphaGenome output. It filters returned rows by:

- target gene name or Ensembl gene id,
- data source `gtex`,
- GTEx tissue `Lung`.

For each locus, it writes:

```text
output/annotation/alphagenome/<locus_id>/<GENE>_alphagenome_filtered_scores.parquet
output/annotation/alphagenome/<locus_id>/<GENE>_alphagenome_variant_scores.csv
output/annotation/alphagenome/<locus_id>/<GENE>_alphagenome_variant_scores_histogram.png
output/annotation/alphagenome/<locus_id>/<GENE>_alphagenome_variant_scores_histogram_trimmed.png
output/prelim/<locus_id>/<GENE>_alphagenome_variant_window_eligibility.csv
output/prelim/<locus_id>/<GENE>_alphagenome_interval_genes.csv
```

The compact per-variant AlphaGenome score file contains:

```text
source_variant_id
alphagenome_raw_mean
alphagenome_quantile_mean
alphagenome_track_count
```

The geometric N20 SuSiNE handoff uses `alphagenome_quantile_mean` as the signed
annotation source.

## 6. Annotation transform and `mu_0` handoff

The production handoff script is:

```text
eQTL_annotations_for_susine/scripts/prepare_geometric_n20_mu0.py
```

Expected command:

```bash
python scripts/prepare_geometric_n20_mu0.py \
  --manifest config/loci_manifest_sample_100_per_chrom.csv \
  --selection config/annotation_selection_geometric_n20.csv \
  --output-dir output/susine_mu0/geometric_n20
```

For each selected locus, the script aligns:

- LD order,
- Phase 1 master variants,
- AlphaGenome compact variant scores.

The join key is:

```text
LD order id == AlphaGenome source_variant_id
```

Missing AlphaGenome scores are treated as neutral annotation evidence:

```text
alphagenome_quantile_mean_filled = 0
annotation_missing = TRUE
```

Let `q_j` be the filled signed AlphaGenome quantile score. The transform into
the SuSiNE annotation vector is:

```text
q_star_j = clip(q_j, -1 + 1e-4, 1 - 1e-4)
a_raw_j = qnorm((q_star_j + 1) / 2)
a_clip_j = clip(a_raw_j, -2.5, 2.5)
a_j = a_clip_j / sqrt(mean(a_clip^2))
```

The annotation is RMS-normalized within locus and is not mean-centered. The
reason is semantic: `0` represents no directional AlphaGenome evidence.

The standardized marginal effect calibration is computed from the z-only path.
For variant `j`:

```text
adj_j = (n_j - 1) / (z_j^2 + n_j - 2)
beta_hat_std_j = z_j * sqrt(adj_j / (n_j - 1))
shat2_std_j = 1 / (n_j - 1)
```

The per-locus baseline annotation scale is:

```text
c_rms_l = sqrt(mean(beta_hat_std^2))
c_cap_l = quantile(abs(beta_hat_std), 0.95) / max(abs(a))
baseline_c_l = min(c_rms_l, c_cap_l)
mu0_j = baseline_c_l * a_j
```

The cap prevents the RMS scale from creating a maximum absolute prior mean that
is too large relative to the observed marginal standardized effect scale.

The script writes:

```text
output/susine_mu0/geometric_n20/mu0_locus_summary.csv
output/susine_mu0/geometric_n20/mu0_variant_table.parquet
output/susine_mu0/geometric_n20/mu0_variant_annotations_all.csv
output/susine_mu0/geometric_n20/per_locus_annotations/<GENE>_mu0_variant_annotations.csv
```

Required per-locus annotation columns for downstream `test_susine` are:

```text
variant_id
ld_matrix_index
annotation_a
beta_hat_std
baseline_c_l
var_y_hat_from_slope
var_y_hat_from_se
```

## 7. Sync into `test_susine`

The downstream input-sync workbook is:

```text
test_susine/vignettes/real data pipeline/prepare_real_data_inputs_workbook.Rmd
```

Its configuration is:

```r
source_repo_root <- "/storage/work/mgc5166/Annotations/eQTL_annotations_for_susine"
dest_root <- here("data", "real_case_studies", "geometric_n20_loci")
source_mu0_name <- "geometric_n20"
source_annotation_summary <- NULL
loci <- NULL
```

The workbook calls:

```r
test_susine::sync_real_data_inputs(
  source_repo_root = source_repo_root,
  dest_root = dest_root,
  source_mu0_name = source_mu0_name,
  source_annotation_summary = source_annotation_summary,
  loci = loci
)
```

For each row in:

```text
<source_repo_root>/output/susine_mu0/geometric_n20/mu0_locus_summary.csv
```

the sync step copies:

```text
<source_repo_root>/output/z_score/<locus_id>/<GENE>_GTEx_z_scores.csv
<source_repo_root>/output/ld/<locus_id>/<GENE>_phase1_LD_R_long.parquet
<source_repo_root>/output/ld/<locus_id>/<GENE>_LD_variant_order.tsv
<source_repo_root>/output/ld/<locus_id>/<GENE>_phase1_master_variants.csv
<source_repo_root>/output/susine_mu0/geometric_n20/per_locus_annotations/<GENE>_mu0_variant_annotations.csv
```

into:

```text
test_susine/data/real_case_studies/geometric_n20_loci/<locus_id>/
```

The sync also copies study-level provenance files:

```text
mu0_locus_summary.csv
geometric_n20_annotation_batch_summary.csv
```

if the AlphaGenome batch summary exists.

Then it builds:

```text
test_susine/data/real_case_studies/geometric_n20_loci/locus_manifest.csv
```

The manifest contains, per locus:

- `dataset_bundle_id`
- `locus_id`
- `gene_name`
- `gtex_tissue`
- `gtex_chrom`
- `n_variants`
- `n_sample_median`
- `baseline_c_l`
- relative paths to z-scores, LD long parquet, LD order, master variants, and
  annotations.

The preview chunk loads the first bundle with:

```r
test_susine:::load_real_data_locus_bundle(
  locus_id = preview_locus,
  manifest_path = file.path(dest_root, "locus_manifest.csv"),
  repo_root = here()
)
```

Bundle validation checks include:

- annotation `ld_matrix_index` is contiguous;
- LD order `index` is contiguous;
- LD order IDs exactly match annotation `variant_id`;
- z-score file has one row per LD-order variant;
- master variant table exactly matches annotation `variant_id`;
- LD long parquet reconstructs a symmetric matrix with unit diagonal.

The returned bundle contains:

```text
bundle$z
bundle$R
bundle$a
bundle$variant_map
bundle$n_sample
bundle$n_sample_min
bundle$n_sample_max
bundle$n_sample_median
bundle$baseline_c_l
bundle$beta_hat_std
bundle$var_y_hat_from_slope
bundle$var_y_hat_from_se
```

The sample size passed to RSS fits is the median variant-level effective sample
size within the locus.

## 8. Real-data RSS job configuration

The real-data run-control workbook is:

```text
test_susine/vignettes/real data pipeline/run_control_workbook_real_data_ensemble.Rmd
```

Production settings, after turning `DRY_RUN <- FALSE`, are:

```r
job_name <- "real_data_ensemble_geometric_n20"
manifest_path <- here("data", "real_case_studies", "geometric_n20_loci", "locus_manifest.csv")

L <- 10L
max_iter <- 100L
tol <- 1e-5
cs_coverage <- 0.95
cs_min_purity <- 0.5
jsd_threshold <- 0.15
softmax_temperature <- 1
estimate_residual_variance <- FALSE

c_grid_8 <- seq(0, 1.5, length.out = 8)
sigma_grid_8 <- c(0.01, 0.03, 0.07, 0.1, 0.2, 0.4, 0.7, 1.0)
```

SLURM settings in the workbook:

```r
email <- "mgc5166@psu.edu"
slurm_time <- "05:59:59"
slurm_mem <- "8G"
cpus_per_task <- 1L
slurm_partition <- NULL
slurm_account <- "statsresearch_sc_default"
```

The workbook calls:

```r
test_susine::build_real_data_job_config(...)
test_susine::write_real_data_job_artifacts(
  job_config,
  run_task_script = here("inst", "scripts", "run_real_data_task.R")
)
```

For each locus, the run manifest contains:

- 64 SuSiNE RSS functional-grid fits:
  `8 c values x 8 sigma_0^2 values`;
- 1 `susieR::susie_rss()` anchor fit at `sigma_0^2 = 0.2`.

The functional-grid fits use:

```r
susine::susine_rss(
  L = 10,
  z = bundle$z,
  R = bundle$R,
  n = bundle$n_sample,
  mu_0 = c_value * bundle$a,
  sigma_0_2 = sigma_0_2_scalar,
  prior_update_method = "none",
  estimate_residual_variance = FALSE,
  convergence_method = "elbo",
  tol = 1e-5,
  max_iter = 100
)
```

The anchor fit uses:

```r
susieR::susie_rss(
  z = bundle$z,
  R = bundle$R,
  n = bundle$n_sample,
  L = 10,
  scaled_prior_variance = 0.2,
  estimate_prior_variance = FALSE,
  estimate_residual_variance = FALSE,
  residual_variance = 1,
  check_prior = FALSE,
  max_iter = 100,
  tol = 1e-5,
  convergence_method = "elbo",
  coverage = 0.95,
  min_abs_corr = 0.5
)
```

The generated SLURM script loads R, sets:

```bash
export R_LIBS_USER="/storage/home/mgc5166/R/x86_64-pc-linux-gnu-library/4.3"
export SUSINE_DEV="1"
```

and runs:

```bash
Rscript "$RUN_TASK_SCRIPT" \
  --job-name "$JOBNAME" \
  --task-id "$TASK_ID" \
  --job-root "$JOB_ROOT" \
  --config-path "$CONFIG_PATH"
```

When `SUSINE_DEV=1`, `run_real_data_task.R` attempts to load a sibling local
development checkout:

```text
../susine
```

and validates that `susine::susine_rss()` has the expected arguments
`prior_update_method` and `estimate_residual_variance`.

## 9. Per-task outputs

Each SLURM array task maps to one locus. Within a task, fits are grouped by
`flush_group`, primarily to periodically write compact outputs rather than keep
all 65 fits in memory.

The task runner writes staged outputs under:

```text
test_susine/output/slurm_output/<job_name>/<parent_job_id>/staging/task-<task_id>/
```

or the equivalent HPC output root.

Per-flush outputs include:

```text
flush-XXX_run_metrics.csv
flush-XXX_effect_summaries.csv
flush-XXX_credible_set_membership.parquet
flush-XXX_variant_posteriors.parquet
flush-XXX_effect_posteriors.parquet
flush-XXX_elbo_trace.csv
flush-XXX_validation.csv
flush-XXX_dataset_metrics.csv
flush-001_fit_file_index.csv
raw_fits/run-<run_id>_fit.rds
```

Per-run metrics include:

- final ELBO;
- residual variance trace endpoint;
- `h2_proxy_std = 1 - sigma_2_final`, clipped to `[0, 1]`;
- posterior-mean PVE proxy;
- PIP entropy and PIP threshold counts;
- top-k PIP mass;
- credible-set counts, sizes, purities, and overlap rates;
- annotation correlation diagnostics;
- fit path and wall time.

The posterior-mean PVE proxy is computed as:

```text
b = posterior mean effect vector
PVE_proxy = b' R b
```

clipped to `[0, 1]`.

For SuSiNE RSS fits:

```text
b = colSums(effect_fits$b_hat)
```

where `effect_fits$b_hat` is already alpha-weighted in the current SuSiNE fit
object.

For `susieR` fits, the code constructs alpha-weighted effect summaries from:

```text
alpha * mu
alpha * mu2
```

Credible sets are built per effect by sorting normalized `alpha_l` until the
cumulative mass reaches `rho = 0.95`. Filtered sets require minimum absolute LD
purity at least `0.5`.

The per-locus dataset metrics include:

- `M1`;
- high-LD pair counts at `|r| >= 0.95`;
- z-score concentration metrics;
- annotation correlations with `z` and `beta_hat_std`;
- `var_y_hat` diagnostics;
- an `(R, z, n)` consistency diagnostic.

The `(R, z, n)` consistency diagnostic computes:

```text
mom_sigma2_hat = 1 - z' (R + lambda I)^-1 z / n
```

with `lambda = 1e-4`. It flags a locus if Cholesky decomposition fails or
`mom_sigma2_hat < 0`, because those conditions indicate an internally
inconsistent RSS triple.

## 10. Collection and derived artifacts

The collection workbook is:

```text
test_susine/vignettes/real data pipeline/collect_results_workbook_real_data_ensemble.Rmd
```

Its production configuration is:

```r
job_name <- "real_data_ensemble_geometric_n20"
parent_job_id <- "52906940"
output_root <- here("output")
```

It calls:

```r
test_susine::index_staging_outputs(...)
test_susine::collect_real_data_results(...)
test_susine::validate_real_data_outputs(...)
```

The consolidated output directory is:

```text
test_susine/output/slurm_output/real_data_ensemble_geometric_n20/52906940/aggregated/
```

The collector writes:

```text
run_manifest.csv
run_metrics_full.csv
effect_summaries_full.csv
dataset_metrics.csv
validation.csv
fit_file_index.csv
elbo_trace_full.csv
multimodal_metrics.csv
cluster_membership.csv
aggregation_weights_cluster_weight.csv
cluster_summary.csv
functional_grid_summary.csv
susie_anchor_summary.csv
highest_weight_refit_summary.csv
highest_weight_refit_basin_r2_drift.csv
paper_real_data_ensemble_summary.csv
top_variants_cluster_weight.csv
run_comparisons.csv
variant_posteriors_dataset/
effect_posteriors_dataset/
credible_set_membership_dataset/
pairwise_pip_jsd.parquet
aggregated_variant_pips_cluster_weight.parquet
highest_weight_refit_variant_posteriors.parquet
metric_inventory.csv
```

The collector uses `metric_inventory.csv` as the paper-facing artifact index.

### 10.1 Functional-grid clustering

Within each locus, the collector forms a matrix of PIP vectors across
functional-grid runs. It computes pairwise Jensen-Shannon divergence with:

```text
JSD(p, q) = 0.5 * KL(p || m) + 0.5 * KL(q || m)
m = 0.5 * (p + q)
```

using a small epsilon for numerical stability. The implementation calls this
quantity `js_distance`, although it is the Jensen-Shannon divergence rather
than its square root.

It then runs complete-linkage hierarchical clustering on the JSD matrix and
cuts the tree at:

```text
jsd_threshold = 0.15
```

Cluster summary outputs include:

- cluster id;
- representative run id;
- cluster size;
- cluster frequency;
- representative ELBO;
- cluster weight.

### 10.2 Aggregation weights

The collector constructs `aggregation_weights_cluster_weight.csv` by selecting
the highest-ELBO representative within each JSD cluster and assigning
representative weights using:

```text
w_rep = softmax(ELBO_rep / temperature) / cluster_frequency
w_rep = w_rep / sum(w_rep)
```

with `temperature = 1`.

The displayed annotation-weight panel uses these run-level weights through
`functional_grid_summary$agg_weight_run`.

Implementation note: the PIP vector written to
`aggregated_variant_pips_cluster_weight.parquet` is generated by
`aggregate_pip_matrix(..., method = "cluster_weight")`. In the available source,
that function gives each cluster total weight `1/K` and does within-cluster ELBO
softmax. By contrast, the run-level `agg_weight_run` table uses the
frequency-adjusted representative weighting above. The ensemble posterior mean
and PVE are computed from the run-level `agg_weight_run` table. This means the
paper PIP aggregation and paper PVE/annotation-weight aggregation may not be
using exactly the same cluster-weight formula unless the HPC run used a newer
implementation not present locally.

This should be checked before final manuscript submission.

### 10.3 Warm refit from highest-weight source

For each locus, the collector identifies the highest-weight source run by:

1. descending `agg_weight_run`,
2. descending `elbo_final`,
3. ascending `run_id`.

It then runs a warm baseline refit initialized from that source model:

```r
susine::susine_rss(
  L = 10,
  z = bundle$z,
  R = bundle$R,
  n = bundle$n_sample,
  mu_0 = 0,
  sigma_0_2 = 0.2,
  prior_inclusion_weights = rep(1 / p, p),
  prior_update_method = "none",
  estimate_residual_variance = FALSE,
  convergence_method = "elbo",
  tol = 1e-5,
  max_iter = 100,
  init_effect_fits = source_fit$effect_fits[c("alpha", "b_hat", "b_2_hat")]
)
```

This refit tests whether an annotation-discovered basin remains stable after the
annotation prior is removed.

The refit outputs are written to:

```text
aggregated/highest_weight_refits/<locus_id>/highest_weight_refit_run_<id>.rds
aggregated/highest_weight_refit_summary.csv
aggregated/highest_weight_refit_variant_posteriors.parquet
```

### 10.4 Effect-basis `r^2` drift

Effect-basis drift compares the per-effect coefficient vectors of two fits.

For each effect, define:

```text
b_l = alpha_l * posterior_mean_l
var(X b_l) = b_l' R b_l
cor(X b_i, X b_j) = (b_i' R b_j) / sqrt((b_i' R b_i)(b_j' R b_j))
```

Effects are matched one-to-one by maximum assignment using the Hungarian
algorithm (`clue::solve_LSAP`). The pairwise drift is:

```text
weighted mean of (1 - cor^2)
```

with weights:

```text
min(var_explained_effect_a, var_explained_effect_b)
```

The collector writes three comparisons:

- highest-weight source vs SuSiE anchor;
- warm refit vs SuSiE anchor;
- warm refit vs highest-weight source.

These are stored in:

```text
highest_weight_refit_basin_r2_drift.csv
```

## 11. Paper figure rendering

The visualization workbook is:

```text
test_susine/vignettes/real data pipeline/visualize_results_workbook_real_data_ensemble.Rmd
```

Production configuration:

```r
job_name <- "real_data_ensemble_geometric_n20"
parent_job_id <- "52906940"
force_recompute_basin_r2_drift <- FALSE

aggregated_dir <- here("output", "slurm_output", job_name, parent_job_id, "aggregated")
figure_dir <- here("output", "slurm_output", job_name, parent_job_id, "figures", "real_data_ensemble")
overall_figure_dir <- file.path(figure_dir, "overall")
locus_figure_dir <- file.path(figure_dir, "per_locus")
```

The final paper copies in:

```text
../Writings/plots/real_data_case_study/
```

appear to be copied from:

```text
test_susine/output/slurm_output/real_data_ensemble_geometric_n20/52906940/figures/real_data_ensemble/overall/
```

### 11.1 Overall annotation-drift-PVE figure

Final paper file:

```text
paper_overall_annotation_drift_pve.png
```

Workbook output paths:

```text
overall/paper_overall_annotation_drift_pve_panel_a.png
overall/paper_overall_annotation_drift_pve_panel_b.png
overall/paper_overall_annotation_drift_pve_panel_c.png
overall/paper_overall_annotation_drift_pve_panel_d.png
overall/paper_overall_annotation_drift_pve.png
overall/paper_overall_annotation_drift_pve.pdf
```

The combined figure is rendered with `patchwork` at:

```text
width = 7.52
height = 7.8
dpi = 300
bg = "white"
```

Panel A: annotation influence.

Source table:

```text
functional_grid_summary.csv
```

Data construction:

```r
weighted_c_value =
  sum(agg_weight_run * c_value) / sum(agg_weight_run)

weighted_sigma_0_2_scalar =
  sum(agg_weight_run * sigma_0_2_scalar) / sum(agg_weight_run)

off_c0_weight =
  sum(agg_weight_run[c_value > 0]) / sum(agg_weight_run)

max_model_weight = max(agg_weight_run)
```

Loci are labeled when:

```text
weighted_c_value > 0.1 OR weighted_sigma_0_2_scalar > 0.1
```

Panel B: PIP drift vs effect-basis drift.

Source tables:

```text
highest_weight_refit_basin_r2_drift.csv
paper_real_data_ensemble_summary.csv
variant_posteriors_dataset/
```

The x-axis is effect-basis `r^2` drift from the SuSiE baseline. The y-axis is
PIP JSD from the SuSiE baseline. Each locus contributes up to two points:

- highest-weight model;
- warm refit.

For the highest-weight model, PIP JSD is recomputed directly between:

```text
susie_anchor_run_id
highest_weight_source_run_id
```

For the warm refit, PIP JSD comes from:

```text
paper_real_data_ensemble_summary.csv: jsd_susie_anchor_vs_refit
```

Loci are labeled for the highest-weight-model point when:

```text
basin_r2_drift > 0.1 OR pip_jsd > 0.3
```

Panel C: posterior-mean PVE.

Source tables:

```text
susie_anchor_summary.csv
highest_weight_refit_summary.csv
```

The panel compares:

- SuSiE anchor PVE,
- SuSiNE ensemble PVE,
- warm-refit PVE.

The PVE proxy is:

```text
b' R b
```

where `b` is the posterior mean effect vector on standardized genotype and
standardized phenotype scale.

Loci are labeled for the SuSiNE ensemble point when:

```text
susie_anchor_pve > 0.15 OR comparison_pve > 0.15
```

Panel D: effect diffuseness.

Source tables:

```text
paper_real_data_ensemble_summary.csv
effect_posteriors_dataset/
```

The panel uses only:

- SuSiE baseline run ids;
- highest-weight source run ids.

For each selected run and effect, it computes entropy on the top 95 percent
posterior mass core:

```text
effect_k_eff_core95 = exp(core95_entropy(alpha_l))
```

Effects with negligible alpha mass are dropped:

```text
alpha_mass > 1e-8
```

The histogram bins are:

```text
[0, 1.5)
[1.5, 5)
[5, 20)
[20, 100)
[100, 2500)
[2500, 4000)
[4000, Inf)
```

### 11.2 Selected-locus zoom figure

Final paper file:

```text
paper_selected_locus_zoom_lollipop_l2_drift.png
```

Workbook output paths:

```text
overall/paper_selected_locus_zoom_lollipop_drift.png
overall/paper_selected_locus_zoom_lollipop_drift.pdf
overall/paper_selected_locus_zoom_lollipop_l2_drift.png
overall/paper_selected_locus_zoom_lollipop_l2_drift.pdf
```

The selected-locus figure is rendered with `patchwork` at:

```text
width = 7.52
height = 9.2
dpi = 300
bg = "white"
```

The selected locus tokens are:

```r
selected_locus_tokens <- c("rrp7a", "arsa", "prpf38b")
```

Each token is resolved by first matching `^token` against `locus_id`, then by a
case-insensitive substring fallback.

For each locus, the left panel is a lollipop plot of the 15 variants with the
largest absolute baseline-to-ensemble PIP difference:

```text
abs_delta_baseline_ensemble = abs(aggregated_pip - baseline_pip)
```

The three plotted point types are:

- Baseline: nearest functional-grid run to `c = 0`, `sigma_0^2 = 0.2`;
- SuSiNE ensemble: `aggregated_pip`;
- Warm refit: PIP from `highest_weight_refit_variant_posteriors.parquet`.

The right panel in the final L2 version plots, for each functional-grid run:

```text
x = c_value
y = sqrt(sum((pip_run - baseline_pip)^2))
color = sigma_0_2_scalar
```

The black open circle marks the highest-weight source run for that locus, using
the same selection rule as the warm-refit stage:

```text
desc(agg_weight_run), desc(elbo_final), run_id
```

The figure also renders a JSD version, but the final paper PNG named in this
trace is the L2-drift version.

## 12. Reproduction recipe

This is the practical end-to-end recipe implied by the codebase.

### 12.1 Build upstream GTEx z, LD, and annotations

From:

```text
eQTL_annotations_for_susine/
```

install Python dependencies:

```bash
pip install -r requirements.txt
```

Make sure GTEx Lung allpairs parquet files are present under:

```text
data/gtex/
```

Make sure the geometric N20 selection exists:

```bash
python scripts/lock_geometric_n20_sample.py
```

Run Phase 1 in processing mode for the selected loci. The exact batch command
depends on the SLURM array setup, but the single-task entry point is:

```bash
python scripts/run_phase1_task.py \
  --manifest config/loci_manifest_sample_100_per_chrom.csv \
  --task-id <task_id> \
  --chunk-size 1 \
  --run-mode processing
```

After all selected loci are processed, run AlphaGenome scoring:

```bash
python scripts/run_annotation_batch.py \
  --manifest config/loci_manifest_sample_100_per_chrom.csv \
  --selection config/annotation_selection_geometric_n20.csv \
  --summary-path output/annotation/alphagenome/geometric_n20_annotation_batch_summary.csv
```

Then build SuSiNE annotation handoff files:

```bash
python scripts/prepare_geometric_n20_mu0.py \
  --manifest config/loci_manifest_sample_100_per_chrom.csv \
  --selection config/annotation_selection_geometric_n20.csv \
  --output-dir output/susine_mu0/geometric_n20
```

Expected upstream validation:

- `output/susine_mu0/geometric_n20/mu0_locus_summary.csv` exists.
- `output/susine_mu0/geometric_n20/per_locus_annotations/<GENE>_mu0_variant_annotations.csv`
  exists for every selected locus.
- Each annotation CSV has `variant_id` exactly matching the LD order.
- The AlphaGenome batch summary reports each selected locus as `completed` or
  `skipped_existing`.

### 12.2 Sync into `test_susine`

From:

```text
test_susine/
```

run the input-prep workbook:

```text
vignettes/real data pipeline/prepare_real_data_inputs_workbook.Rmd
```

or call:

```r
test_susine::sync_real_data_inputs(
  source_repo_root = "/storage/work/mgc5166/Annotations/eQTL_annotations_for_susine",
  dest_root = here::here("data", "real_case_studies", "geometric_n20_loci"),
  source_mu0_name = "geometric_n20"
)
```

Expected downstream validation:

- `data/real_case_studies/geometric_n20_loci/locus_manifest.csv` exists.
- `load_real_data_locus_bundle()` succeeds for at least one locus.
- Bundle diagnostics show finite correlations and a plausible sample-size range.

### 12.3 Build and submit the real-data ensemble job

Run:

```text
vignettes/real data pipeline/run_control_workbook_real_data_ensemble.Rmd
```

with:

```r
DRY_RUN <- FALSE
```

This writes SLURM-ready artifacts under:

```text
output/temp/real_data_ensemble_geometric_n20/
output/slurm_scripts/real_data_ensemble_geometric_n20.slurm
```

Submit:

```bash
sbatch output/slurm_scripts/real_data_ensemble_geometric_n20.slurm
```

The run traced here used:

```text
parent_job_id = 52906940
```

### 12.4 Collect the real-data ensemble job

Run:

```text
vignettes/real data pipeline/collect_results_workbook_real_data_ensemble.Rmd
```

with:

```r
job_name <- "real_data_ensemble_geometric_n20"
parent_job_id <- "52906940"
```

Expected validation:

- `validate_real_data_outputs(collected$output_dir)$ok` is `TRUE`.
- `metric_inventory.csv` exists in the aggregated output directory.
- `paper_real_data_ensemble_summary.csv` exists.
- `aggregated_variant_pips_cluster_weight.parquet` exists.
- `highest_weight_refit_variant_posteriors.parquet` exists.
- `highest_weight_refit_basin_r2_drift.csv` exists.

### 12.5 Render paper figures

Run:

```text
vignettes/real data pipeline/visualize_results_workbook_real_data_ensemble.Rmd
```

with:

```r
job_name <- "real_data_ensemble_geometric_n20"
parent_job_id <- "52906940"
force_recompute_basin_r2_drift <- FALSE
```

The workbook writes:

```text
output/slurm_output/real_data_ensemble_geometric_n20/52906940/figures/real_data_ensemble/overall/paper_overall_annotation_drift_pve.png
output/slurm_output/real_data_ensemble_geometric_n20/52906940/figures/real_data_ensemble/overall/paper_selected_locus_zoom_lollipop_l2_drift.png
```

Copy those two files into:

```text
../Writings/plots/real_data_case_study/
```

to match the current paper figure paths.

## 13. Artifact inventory expected for final provenance

For a fully locked reproducibility bundle, archive or record checksums for:

Upstream GTEx and Phase 1 inputs:

- GTEx Lung allpairs parquet files used for the selected chromosomes;
- `config/loci_manifest_sample_100_per_chrom.csv`;
- `config/annotation_selection_geometric_n20.csv`;
- per-locus z-score CSVs;
- per-locus LD order TSVs;
- per-locus LD long parquet files;
- per-locus phase1 master variants CSVs.

Upstream AlphaGenome and annotation handoff:

- `output/annotation/alphagenome/geometric_n20_annotation_batch_summary.csv`;
- per-locus AlphaGenome variant score CSVs;
- per-locus AlphaGenome eligibility CSVs;
- `output/susine_mu0/geometric_n20/mu0_locus_summary.csv`;
- `output/susine_mu0/geometric_n20/mu0_variant_table.parquet`;
- per-locus `*_mu0_variant_annotations.csv` files.

Downstream fitting inputs and outputs:

- `test_susine/data/real_case_studies/geometric_n20_loci/locus_manifest.csv`;
- `output/run_history/real_data_ensemble_geometric_n20/52906940/job_config.json`;
- `output/run_history/real_data_ensemble_geometric_n20/52906940/run_manifest.csv`;
- `output/slurm_output/real_data_ensemble_geometric_n20/52906940/aggregated/metric_inventory.csv`;
- all files listed in the aggregated `metric_inventory.csv`;
- the final two paper PNGs and PDFs, if PDFs were rendered.

## 14. Open questions to resolve

These are the points that could not be verified from the local checkout alone.

1. Confirm whether the real-data study is intended to be described as 10 loci
   or 20 loci. The manuscript says 10, but the pipeline and job name are
   geometric N20.
2. Confirm the final source root used on HPC:
   `/storage/work/mgc5166/Annotations/eQTL_annotations_for_susine` versus
   `/storage/home/mgc5166/work/Annotations/eQTL_annotations_for_susine`.
3. Confirm that `parent_job_id = 52906940` is the exact run used for the two
   final paper PNGs.
4. Confirm whether any selected loci failed the `susieR` anchor fit or were
   excluded from the final plots due to missing outputs.
5. Confirm the AlphaGenome client/model version used by the production run.
   The source code records the API call pattern but not the remote model
   version.
6. Check the cluster-weight implementation discrepancy noted in Section 10.2:
   `agg_weight_run` and ensemble posterior means use one weighting formula,
   while `aggregated_pip` appears to use another in the local source.
7. If the final manuscript will claim exact reproducibility, archive the
   generated `locus_manifest.csv`, `job_config.json`, `run_manifest.csv`, and
   `metric_inventory.csv` for job `52906940` with checksums.
