# Real Data Analysis Pipeline Technical Trace

Date written: 2026-05-22
Last updated: 2026-06-20

This document traces the real-data eQTL case-study pipeline that produced the
paper figures:

- `../Writings/plots/real_data_case_study/paper_overall_annotation_drift_pve.png`
- `../Writings/plots/real_data_case_study/paper_overall_annotation_matched_basis_drift_pve.png`
- `../Writings/plots/real_data_case_study/paper_selected_locus_zoom_lollipop_l2_drift.png`
- `../Writings/plots/real_data_case_study/paper_supplement_real_vs_sim_multimodality.png`
- `../Writings/plots/real_data_case_study/paper_real_data_locus_summary_table.csv`
- `../Writings/plots/real_data_case_study/paper_real_data_pip_gt05_totals.csv`
- `../Writings/plots/real_data_case_study/paper_real_data_matched_basis_drift_summary.csv`
- `../Writings/plots/real_data_case_study/paper_real_data_annotation_weight_labels.csv`

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
- The final paper directory currently contains the two main paper PNGs, one
  supplement PNG, paper-facing CSV summaries, and older paper-note artifacts.
  The main PNGs and generated paper summary CSVs have local write times on
  2026-06-15.
- The production upstream annotation/data-prep source root used for this run is:
  `/storage/home/mgc5166/work/Annotations/eQTL_annotations_for_susine`.
  Some source code and older notebooks still show
  `/storage/work/mgc5166/Annotations/eQTL_annotations_for_susine`; treat that as
  an outdated/default path unless the two locations are verified to be aliases
  on HPC.
- User-confirmed provenance updates, 2026-05-26: the real-data study uses 20
  loci, `52906940` is the parent job id for the run used in the paper, and no
  selected loci failed.
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
3. whether the annotation-informed ensemble or warm refit changes local genetic
   variance, using corrected expected PVE when available and otherwise falling
   back to the posterior-mean PVE proxy;
4. whether annotation-informed fits sharpen per-effect posterior support.

The selected-locus zoom figure (the paper L2-drift figure) then drills into six
loci:

- `ydjc` -> expected full locus id `ydjc_chr22_lung`
- `arsa` -> expected full locus id `arsa_chr22_lung`
- `rrp7a` -> expected full locus id `rrp7a_chr22_lung`
- `tmtc1` -> expected full locus id `tmtc1_chr12_lung`
- `lgals9` -> expected full locus id `lgals9_chr17_lung`
- `znf280b` -> expected full locus id `znf280b_chr22_lung`

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

The paper real-data study uses the geometric N20 selected panel, i.e. 20
protein-coding GTEx Lung eQTL loci. Earlier 10-locus wording in drafts and
legacy notebooks is outdated. The relevant source note is:

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
locally, but the final paper-side `paper_real_data_ensemble_summary.csv`
contains 20 completed loci:

- `arhgef37_chr5_lung`
- `arsa_chr22_lung`
- `cd1d_chr1_lung`
- `erbb3_chr12_lung`
- `lgals9_chr17_lung`
- `nbpf19_chr1_lung`
- `pcdhac2_chr5_lung`
- `prpf38b_chr1_lung`
- `pus7l_chr12_lung`
- `pzp_chr12_lung`
- `rnf14_chr5_lung`
- `rrp7a_chr22_lung`
- `serf1a_chr5_lung`
- `sil1_chr5_lung`
- `snrnp35_chr12_lung`
- `supt4h1_chr17_lung`
- `tmtc1_chr12_lung`
- `trabd_chr22_lung`
- `ydjc_chr22_lung`
- `znf280b_chr22_lung`

The local `loci_manifest_sample_100_per_chrom.csv` also contains the selected
zoom-locus candidates:

- `prpf38b_chr1_lung`
- `arsa_chr22_lung`
- `rrp7a_chr22_lung`

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

The code creates the remote client with:

```python
dna_model = dna_client.create(api_key)
```

No local code path or generated summary file records an AlphaGenome remote model
version. Therefore the exact remote model/version used by the API at production
time is not recoverable from the current repo artifacts.

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
source_repo_root <- "/storage/home/mgc5166/work/Annotations/eQTL_annotations_for_susine"
dest_root <- here("data", "real_case_studies", "geometric_n20_loci")
source_mu0_name <- "geometric_n20"
source_annotation_summary <- NULL
loci <- NULL
```

The default argument in `test_susine/R/real_data_pipeline.R` still points to
`/storage/work/mgc5166/Annotations/eQTL_annotations_for_susine`; override it
with the production path above when reproducing the paper run unless the two
HPC paths are verified aliases.

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

The run-time posterior-mean PVE proxy is computed as:

```text
b = posterior mean effect vector
PVE_proxy = b' R b
```

clipped to `[0, 1]`.

For the current paper PVE exhibit, the visualization workbook first tries to
load a corrected local-genetic-variance estimand from the sibling
`hg2_corrected/` directory, produced by `inst/scripts/recompute_real_data_hg2.R`.
When those files are present, the plotted quantity is:

```text
hg2_expected_pve = E[var(Xb)|y] / var(y)
```

This is the estimand switched in by commit `4e0721b` (see
`R/heritability.R`, `hg2_components()` / `hg2_components_from_moments()`, and
`refs/decisions/heritability_estimand_decision_2026-06-09.md`). It replaces the
old "variance of the posterior-mean predictor" point estimate
`var(E[Xb|y])/var(y)` with the posterior mean of the genetic variance, which
decomposes additively as

```text
E[var(Xb)|y] / var(y)  =  m' R m / var(y)            (postmean part, = hg2_postmean)
                       +  tr(R Cov(b|y)) / var(y)    (within-fit correction, = hg2_uncertainty)
```

with `m = colSums(b_hat)` the unconditional posterior-mean coefficient vector
and, under the SuSiE/SuSiNE single-effect (SER) structure,
`Cov(b_l) = diag(b_2_hat[l, ]) - m_l m_l'`, so
`tr(R Cov(b|y)) = sum_l [ sum_j R_jj b_2_hat[l, j] - m_l' R m_l ]`. For
standardized RSS `var(y) = 1`, the Gram matrix is `R`, and `b_hat`/`b_2_hat` are
the unconditional per-effect moments stored on the susine fit (`alpha * mu_1`,
`alpha * (sigma_1_2 + mu_1^2)`); for a susieR fit they are reconstructed as
`alpha * mu` and `alpha * mu2`. Components are floored at 0 and `expected_pve` is
NOT capped at 1, so the additive decomposition holds exactly.

The recompute is downstream-only (no refitting): the per-effect unconditional
moments are read from the persisted `effect_posteriors` parquet
(`alpha`, `alpha_b_hat`, `alpha_b_2_hat`), falling back to reloading the saved
`.rds` fits only when those columns are absent. The estimand is computed per fit
(SuSiE anchor, every functional-grid fit, and the warm refit reloaded from
`aggregated/highest_weight_refits/`). The **between-fit** term enters at the
ensemble level: the per-locus ensemble PVE
`hg2_expected_pve_ensemble = sum_k w_k * hg2_expected_pve_k` is a convex
combination of the functional-grid fits' expected PVEs, so it may lie above or
below the max-ELBO member (the script also reports
`hg2_expected_pve_pos_min/max` to bound it). As of the 2026-06-16 source fix,
the recompute uses the same credible-shift + Method-B weights as the collector:
`cluster_weight_credible` (`credible_shift` at threshold `0.05`) followed by
`.cluster_weights_from_hc(..., aggregation = "method_b")$w_full`. This matches
the `agg_weight_run` / `aggregated_pip` rule written to
`aggregation_weights_cluster_weight.csv` (Section 10.2). Corrected
`hg2_corrected/` outputs should be regenerated after this source fix before the
PVE numbers are quoted.

The estimand is used consistently for the SuSiE anchor, SuSiNE ensemble, and
warm refit. If `hg2_corrected/` is absent, the workbook (and the export script)
fall back to the posterior-mean proxy below and emit a warning.

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
functional-grid runs. Clustering and aggregation use the locked-in
`cluster_weight_credible` spec (pulled through from the simulation pipeline on
2026-06-11), i.e. the **credible-shift** PIP distance

```text
credible_shift(p, q) = max_j [ max(p_j, q_j) * |p_j - q_j| ]
```

with complete-linkage hierarchical clustering cut at threshold `0.05`
(`.cluster_weight_specs[["cluster_weight_credible"]]`). The cache is built by
`prepare_pip_similarity_cache(pip_mat, metric = "credible_shift")` and stored as
`cred_cache`.

Jensen-Shannon divergence is still computed (`pip_cache <-
prepare_pip_similarity_cache(pip_mat)`, `js_distance`) but only as the
**reported** pairwise-divergence diagnostic that feeds `pairwise_pip_jsd`, the
multimodal metrics, and `jsd_to_representative`; it no longer governs the
clustering. The `jsd_threshold = 0.15` value now applies only to the JSD
multimodal `n_clusters` diagnostic.

Cluster summary outputs include:

- cluster id;
- representative run id (highest-ELBO fit in the cluster);
- cluster size;
- cluster frequency;
- representative ELBO;
- cluster weight (total Method-B weight of the cluster).

### 10.2 Aggregation weights

Both the run-level `agg_weight_run` table and the aggregated PIP vector use
**Method B** (frequency-free) from `.cluster_weights_from_hc(..., aggregation =
"method_b")` in `R/run_model.R`:

1. Cluster fits by credible-shift at threshold `0.05` (complete-linkage hclust
   then cutree on `cred_cache$hc`).
2. Each cluster receives a weight proportional to `exp(max-ELBO-in-cluster)`.
3. That cluster weight is split *within* the cluster across all member fits by
   `softmax(ELBO)` (temperature `1`). Every fit therefore contributes a nonzero
   `w_full` weight; there is no single nominee and no `1 / cluster_frequency`
   correction.

The displayed annotation-weight panel uses these run-level weights through
`functional_grid_summary$agg_weight_run = cw$w_full`. The PIP vector written to
`aggregated_variant_pips_cluster_weight.parquet` is generated by
`aggregate_pip_matrix(..., method = "cluster_weight_credible")`, which builds its
own credible-shift dendrogram + Method-B weights internally via
`.aggregate_cluster_method()` (`R/collect_results.R` -> `R/run_model.R`), so the
run-level weights and the aggregated PIP vector are consistent.

History: before 2026-06-11 the real-data path used the legacy `cluster_weight`
method (JSD at `0.15` + nominee-only `softmax(ELBO_rep) / cluster_frequency`,
"Method C"). Any aggregated PIP / `agg_weight_run` outputs produced before that
date must be regenerated to reflect the credible-shift + Method-B rule.
(The output file *names* still carry the `cluster_weight` stem for downstream
compatibility, but the underlying method is `cluster_weight_credible`.)

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

### 10.4 Effect-basis drift (three distinct metrics)

The collector computes three related per-pair drift scalars and writes them to
`highest_weight_refit_basin_r2_drift.csv`. They are easy to conflate. The current
paper overall figure Panel B uses the third metric.

(1) Hungarian effect-basis `r^2` drift
(`real_data_basin_r2_drift_pair_cached`, `R/real_data_pipeline.R`). For each
effect, define:

```text
b_l = colSums-free per-effect vector alpha_l * mu_l  (one row of the effect basis)
var(X b_l) = b_l' R b_l
cor(X b_i, X b_j) = (b_i' R b_j) / sqrt((b_i' R b_i)(b_j' R b_j))
```

Effects are matched one-to-one by maximum assignment using the Hungarian
algorithm (`clue::solve_LSAP`, with cosines shifted by `+1` and the cost matrix
padded to square with zeros). The pairwise drift is the weighted mean of
`(1 - cor^2)` over matched effects, with weights
`min(var_explained_effect_a, var_explained_effect_b)`. This is stored in the
`basin_r2_drift_*` columns. It is scale-free within the matched effects and is
not the x-axis of the current paper Panel B.

(2) Global fitted-y drift in standardized-phenotype units
(`real_data_basis_drift_var_y_pair_cached`, `R/real_data_pipeline.R`). This is a
single scalar with no per-effect matching:

```text
b = colSums(alpha * mu)        (global posterior-mean effect vector of a fit)
basis_drift_var_y = (b_a - b_b)' R (b_a - b_b) = ||b_a - b_b||_R^2
```

The value is floored at 0. With `var(y) = 1` it is directly a fraction of
`var(y)`: the variance of the difference in the two fits' linear predictors.
Unlike (1), it is not scale-invariant. A locus whose two fits explain very
different amounts of `var(y)` can register larger drift even if the fitted basis
is not the main change. This is stored in the `basis_drift_var_y_*` columns for
audit and backward comparison. It is no longer the paper Panel B x-axis.

(3) Hungarian matched basis drift in standardized-phenotype units
(`real_data_matched_basis_drift_var_y_pair_cached`, `R/real_data_pipeline.R`).
This keeps the per-effect Hungarian matching from (1), but reports the numerator
on the standardized phenotype scale:

```text
matched_basis_drift_var_y = sum_matched min(var_a, var_b) * (1 - cor_ab^2)
matched_basis_var_y       = sum_matched min(var_a, var_b)
matched_basis_drift_per_var_y = matched_basis_drift_var_y / matched_basis_var_y
```

Because the RSS fits use standardized phenotypes, `matched_basis_drift_var_y` is
interpretable as the fraction of `var(y)` carried by matched effect bases that
changed orientation. The assignment score uses `1 + cor^2` for real effect pairs
and 0 for dummy padding, so real pairs are preferred over unmatched dummy slots
even when the matched correlation is near zero. This is the current x-axis of
paper overall Panel B.

All three metrics are written for the same three comparisons:

- highest-weight source vs SuSiE anchor
  (`*_highest_weight_vs_susie_anchor`);
- warm refit vs SuSiE anchor (`*_refit_vs_susie_anchor`);
- warm refit vs highest-weight source (`*_refit_vs_highest_weight`).

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

Previously observed local paper-facing artifacts in that directory. After the 2026-06-20 matched-basis update, rerender the workbook and use the `matched_basis_drift` outputs listed below for the paper figure:

| File | Dimensions / rows | Size | Last write time |
|---|---:|---:|---|
| `paper_overall_annotation_drift_pve.png` | 2256 x 2340 | 301,701 | 2026-06-15 16:00:09 |
| `paper_selected_locus_zoom_lollipop_l2_drift.png` | 2256 x 2520 | 708,743 | 2026-06-15 15:59:53 |
| `paper_supplement_real_vs_sim_multimodality.png` | 2100 x 1500 | 230,792 | 2026-06-15 16:00:23 |
| `paper_real_data_locus_summary_table.csv` | 20 data rows | 2,937 | 2026-06-15 16:02:34 |
| `paper_real_data_pip_gt05_totals.csv` | 1 data row | 87 | 2026-06-15 16:02:36 |
| `paper_real_data_ensemble_summary.csv` | 20 data rows | 11,952 | 2026-05-26 13:10:33 |
| `paper_overall_annotation_matched_basis_drift_pve.png` | produced after rerender | not yet recorded | 2026-06-20 code update |
| `paper_real_data_matched_basis_drift_summary.csv` | produced after rerender | not yet recorded | 2026-06-20 code update |
| `paper_real_data_annotation_weight_labels.csv` | produced after rerender | not yet recorded | 2026-06-20 code update |

The visualization workbook writes the table CSV summaries directly into
`../Writings/plots/real_data_case_study/` when that directory exists. The
matched-basis combined figure is also written directly there when that directory
exists. Older PNG paper copies may have been copied manually from the job output
`overall/` directory.

### 11.1 Overall annotation-drift-PVE figure

Final paper files after the 2026-06-20 matched-basis update:

```text
paper_overall_annotation_matched_basis_drift_pve.png
paper_overall_annotation_matched_basis_drift_pve.pdf
```

Workbook output paths:

```text
overall/paper_overall_annotation_drift_pve_panel_a.png
overall/paper_overall_annotation_drift_pve_panel_b.png
overall/paper_overall_annotation_matched_basis_drift_pve_panel_b.png
overall/paper_overall_annotation_drift_pve_panel_c.png
overall/paper_overall_annotation_drift_pve_panel_d.png
overall/paper_overall_annotation_drift_pve.png
overall/paper_overall_annotation_drift_pve.pdf
overall/paper_overall_annotation_matched_basis_drift_pve.png
overall/paper_overall_annotation_matched_basis_drift_pve.pdf
overall/paper_real_data_annotation_weight_labels.csv
overall/paper_real_data_matched_basis_drift_summary.csv
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
  sum(agg_weight_run[abs(c_value) > 1e-12]) / sum(agg_weight_run)

max_model_weight = max(agg_weight_run)
```

The point fill is `off_c0_weight` (viridis "D", limits `[0, 1]`, legend "Total
c>0 weight"). Loci are labeled when:

```text
off_c0_weight > 0.1
```

Panel B: PIP L2 drift vs effect-basis drift.

Source tables:

```text
highest_weight_refit_basin_r2_drift.csv
paper_real_data_ensemble_summary.csv
variant_posteriors_dataset/
```

The x-axis is the Hungarian matched basis drift in standardized-phenotype units
from the SuSiE baseline: the `matched_basis_drift_var_y_*` scalar of Section
10.4(3). The workbook axis title is `Matched basis drift from SuSiE baseline
(fraction var(y))`. The y-axis is PIP L2 distance from the SuSiE baseline. Each
locus contributes up to two points, joined by a grey path:

- highest-weight SuSiNE fit (point color `#0B6E4F`);
- SuSiE warm refit (point color `#C65D00`).

The x-coordinates come straight from `highest_weight_refit_basin` columns
`matched_basis_drift_var_y_highest_weight_vs_susie_anchor` and
`matched_basis_drift_var_y_refit_vs_susie_anchor`. The old total fitted-y drift
columns are retained in `paper_real_data_matched_basis_drift_summary.csv` as
`total_yhat_drift_var_y` for audit. The y-coordinates (PIP L2) are recomputed
live in the workbook from the variant-posterior parquet, not read from a stored
column:

- highest-weight fit: `run_pip_l2(locus_id, susie_anchor_run_id,
  highest_weight_source_run_id)`;
- warm refit: `refit_pip_l2(locus_id, susie_anchor_run_id)`, which joins the
  anchor variant posteriors to `highest_weight_refit_variant_posteriors.parquet`
  by `variant_id`.

JSD values are still collected in `paper_real_data_ensemble_summary.csv` and
auxiliary diagnostics, but the paper-facing overall panel is the L2 version.

Loci are labeled (only on the highest-weight-fit point) when:

```text
matched_basis_drift_var_y > 0.02 OR pip_l2 > 1.0
```

Panel C: local genetic variance / PVE.

Source tables:

```text
susie_anchor_summary.csv
highest_weight_refit_summary.csv
```

The panel compares:

- SuSiE anchor PVE,
- SuSiNE ensemble PVE,
- warm-refit PVE.

The preferred plotted estimand is:

```text
hg2_expected_pve = E[var(Xb)|y] / var(y)
```

selected per fit via `.use_expected_pve <- have_corrected_pve`: anchor uses
`susie_anchor_pve_expected`, ensemble uses `ensemble_pve_expected` (joined by
locus from `hg2_corrected_by_locus.csv`), warm refit uses `refit_pve_expected`.
The figure never mixes estimands within a panel. The axis title uses `h_g^2`.
When corrected-PVE files are unavailable, the workbook falls back to the
posterior-mean proxy `*_pve_postmean_std` (i.e. `b' R b` with `b` the
posterior-mean effect vector on the standardized genotype / phenotype scale).
For each locus a grey vertical segment connects the ensemble and warm-refit
points at the shared anchor x.

Loci are labeled for the SuSiNE ensemble point when:

```text
label_pve = pmax(susie_anchor_pve, susine_ensemble_pve, warm_refit_pve)
label_pve >= 0.19
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
overall/paper_selected_locus_zoom_lollipop_credible_drift.png
overall/paper_selected_locus_zoom_lollipop_credible_drift.pdf
```

The selected-locus figure is rendered with `patchwork` at:

```text
width = 7.52
height = 9.2
dpi = 300
bg = "white"
```

The selected locus tokens for the paper (L2-drift) figure are:

```r
selected_l2_locus_tokens <- c("ydjc", "arsa", "rrp7a", "tmtc1", "lgals9", "znf280b")
selected_l2_n_variants <- 6L
```

Each token is resolved (`resolve_locus_token`) by first matching `^token`
against `locus_id` (case-insensitive), then by a case-insensitive substring
fallback. (A separate deprecated JSD figure uses
`selected_locus_tokens <- c("rrp7a", "arsa", "prpf38b")`.)

The L2-drift figure is a one-row-per-locus grid; each row is FOUR patchwork
cells (`plot_layout(widths = c(0.07, 1, 1, 1))`): a rotated locus-name strip,
then three data columns. Column titles appear only on the first row.

Column 1 — "PIPs for shifted effects" (`paper_lollipop_l2_plot` ->
`plot_aggregated_lollipop(locus, n_variants = 6, selection_metric =
"posterior_mean_delta")`). The 6 displayed variants are the top variants by
`abs_delta_posterior_mean = abs(aggregated_posterior_mean -
baseline_posterior_mean)` (NOT the PIP delta — note the standalone
`plot_aggregated_lollipop` default is `selection_metric = "pip_delta"`, but the
paper L2 column overrides it to `posterior_mean_delta`). Each variant shows three
points plus a grey segment from baseline to ensemble PIP:

- SuSiNE c=0 baseline: nearest functional-grid run to `c = 0`, `sigma_0^2 = 0.2`
  (`nearest_reference_run`, ties broken by `c_distance`, then
  `abs(log(sigma/0.2))`, then `run_id`);
- SuSiNE ensemble: `aggregated_pip`;
- SuSiE warm refit: PIP from `highest_weight_refit_variant_posteriors.parquet`.

Column 2 — "Fitted effect size" (`paper_signal_contribution_l2_plot`; despite
the legacy function name it plots the signed conditional posterior-mean effect,
not a signal contribution). For the same 6 variants it plots the
floored-alpha-weighted conditional effect `E[b_j | gamma_j = 1]` for baseline,
ensemble, and warm refit. Per fit:
`mu_lj = b_hat[l,j]/alpha[l,j]`, `alpha_lj_clip = max(alpha[l,j], 1e-4)`,
`conditional_j = sum_l(alpha_lj_clip * mu_lj) / sum_l(alpha_lj_clip)`; the
ensemble aggregates `(numerator, denominator)` across cluster-weight nominees
weighted by `agg_weight_run` (Option-B aggregation,
`floored_conditional_for_locus`). The x-axis is a signed-log
(`scales::pseudo_log_trans(sigma = 1, base = 2)`) with hardcoded breaks
`{0, +/-1, +/-5, +/-20}`.

Column 3 — "PIP L2 drift from baseline" (`paper_pip_drift_l2_plot` ->
`paper_pip_drift_metric_plot(drift_col = "l2_from_reference_c0_sigma02")`). For
each functional-grid run it plots:

```text
x = c_value
y = sqrt(sum((pip_run - baseline_pip)^2))   (l2_from_reference_c0_sigma02)
color = factor(sigma_0_2_scalar)
```

where `baseline_pip` is the nearest-`(c=0, sigma_0^2=0.2)` reference run. A black
open circle marks the highest-weight source run for that locus, labeled with its
`agg_weight_run` (formatted `%.2f`), selected by:

```text
desc(agg_weight_run), desc(elbo_final), run_id
```

`build_selected_drift_tbl` also computes JSD and credible-shift drift columns
(`jsd_from_reference_c0_sigma02`, `credible_from_reference_c0_sigma02`) for the
sibling figures. The figure additionally renders a deprecated two-column JSD
version (`selected_locus_zoom_figure`, tokens `rrp7a/arsa/prpf38b`) and a
six-locus credible-shift drift version
(`paper_selected_locus_zoom_lollipop_credible_drift.*`, column 3 uses
`credible_from_reference_c0_sigma02`). The final paper PNG named in this trace is
the six-locus L2-drift version above (saved at `width = 7.52`, `height = 8.4`,
even though the chunk header sets `fig.height = 9.2`).

### 11.3 Real-vs-simulation multimodality supplement

Final paper file:

```text
paper_supplement_real_vs_sim_multimodality.png
```

Workbook chunk: `supplement-real-vs-sim-multimodality` (section "5b. Real vs
Simulation Multimodality Benchmark"). Rendered at `width = 7`, `height = 5`,
`dpi = 300`, `bg = "white"`. Also emits a `.pdf`.

This figure benchmarks the real-data loci's posterior multimodality against an
oligogenic *simulation* ensemble. The simulation source is a separate job,
configured in the same workbook:

```r
simulation_multimodality_job_name <- "ensemble_scaling_full"
simulation_multimodality_parent_job_id <- "53547760"
```

read from
`output/slurm_output/ensemble_scaling_full/53547760/consolidated/multimodal_metrics.csv`
(with `model_metrics_full.csv` joined to recover `spec_name`). The simulation
background is filtered to `architecture == "susie2_oligogenic"` and, when
present, `spec_name == "C-CS"`.

IMPORTANT — both axes use the **credible-shift** multimodal metrics, NOT the JSD
metrics. The figure plots:

```text
x = n_clusters_credible   (number of credible-shift clusters at threshold 0.05)
y = max_credible_dist     (largest pairwise credible-shift distance)
```

Both come from `compute_multimodal_metrics()` in `R/run_model.R`:
`n_clusters_credible = length(unique(cutree(cred_cache$hc, h = 0.05)))` and
`max_credible_dist = max(cred_cache$jsd_vals)` where `cred_cache` is built with
`metric = "credible_shift"` (`credible_shift(p,q) = max_j max(p_j,q_j)|p_j-q_j|`).
The JSD-based fields (`mean_jsd`, `median_jsd`, `max_jsd`, `jaccard_top10`,
`mean_pip_var`, `n_clusters` at `jsd_threshold = 0.15`) are still computed and
power the older `multimodality_summary.png` in-workbook diagnostic, but the paper
supplement does NOT use them. The simulation ensemble is drawn as a `geom_bin2d`
density (viridis "C", sqrt-scaled) with an optional white 2D-density contour; the
20 real-data loci are overlaid as labeled green points (label = `locus_id`
stripped of `_lung`).

The chunk skips gracefully (message only) if the simulation
`multimodal_metrics.csv` is missing or the required columns
(`n_clusters_credible`, `max_credible_dist`) are absent.

### 11.4 Paper locus summary CSVs

Two paper-facing CSVs are written to
`../Writings/plots/real_data_case_study/` (and to the job `overall/` dir):

```text
paper_real_data_locus_summary_table.csv   (one row per locus, 20 rows)
paper_real_data_pip_gt05_totals.csv       (one across-locus totals row)
```

There are TWO generators that produce identical schemas:

1. the visualization workbook, section "8. Paper Locus Summary Table" (chunk
   `paper-real-data-locus-summary-table`), which writes both into
   `overall/` and, when the writing directory exists, into the paper dir;
2. the standalone script
   `vignettes/real data pipeline/export_real_data_locus_summary_table.R`
   (defaults `--job-name real_data_ensemble_geometric_n20`,
   `--parent-job-id 52906940`, `--out-dir` = paper dir). Use this to refresh the
   CSVs without re-rendering the full workbook.

Both read `paper_real_data_ensemble_summary.csv` and join the corrected PVE from
`hg2_corrected/` (anchor by `susie_anchor_run_id` -> `hg2_corrected_by_run.csv`;
ensemble by `locus_id` -> `hg2_corrected_by_locus.csv`). The locus-summary
columns are:

```text
locus_id
locus_name                          (gene_name, else upper-cased locus prefix)
chromosome                          (parsed from locus_id "_chr<X>_")
n_susie_pip_gt_05                   = n_pip_gt_05_susie_anchor
n_susine_pip_gt_05                  = n_pip_gt_05_susine_ensemble
n_agreed_pip_gt_05                  = n_pip_gt_05_susie_anchor_and_susine_ensemble
delta_pve_susine_minus_susie        (corrected expected-PVE delta when available,
                                     else postmean: coalesce(pve_*_expected, pve_*))
delta_elbo_warm_refit_vs_susie      = delta_elbo_refit_vs_susie_anchor
total_off_c0_weight                 = coalesce(ensemble_agg_weight_off_c0_total,
                                               weighted off_c0_weight)
elbo_weighted_annotation_scale_c    = weighted_c_value
pip_l2_susine_ensemble_vs_susie     = pip_l2_susie_anchor_vs_susine_ensemble
```

Rows are sorted by descending `pip_l2_susine_ensemble_vs_susie`, then
`locus_name`. PIP counts use threshold `0.5`
(`real_data_pip_count`/`real_data_pip_intersection_count`,
`R/real_data_pipeline.R`). The export script recomputes
`pip_l2_susie_anchor_vs_susine_ensemble` from the anchor `variant_posteriors`
and the aggregated PIP parquet when the column is not already present in the
summary (anchor vs aggregated-ensemble L2 — distinct from the Panel B
anchor-vs-highest-weight-fit L2). The totals CSV is `n_loci`,
`total_susie_pip_gt_05`, `total_susine_pip_gt_05`, `total_agreed_pip_gt_05`.

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
  source_repo_root = "/storage/home/mgc5166/work/Annotations/eQTL_annotations_for_susine",
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

This parent job id is confirmed as the run used for the paper figures.

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
- `highest_weight_refit_basin_r2_drift.csv` exists and includes `matched_basis_drift_var_y_*` columns.
- No selected loci failed. The final paper-side
  `paper_real_data_ensemble_summary.csv` has 20 rows and non-missing anchor,
  highest-weight-source, and warm-refit run ids for all 20 loci.

### 12.4b Recompute the corrected heritability estimand

The corrected PVE estimand `E[var(Xb)|y]/var(y)` is computed downstream (no
refitting) by:

```bash
Rscript inst/scripts/recompute_real_data_hg2.R \
  --job-name real_data_ensemble_geometric_n20 \
  --parent-job-id 52906940 --output-root output
```

This writes
`output/slurm_output/real_data_ensemble_geometric_n20/52906940/hg2_corrected/{hg2_corrected_by_run.csv,hg2_corrected_by_locus.csv}`.
Run it before rendering figures so Panel C and the locus summary table use the
corrected estimand instead of the posterior-mean fallback. After the 2026-06-16
source fix, the script's ensemble PVE uses the same credible-shift + Method-B
weights as the collector; rerun this script to refresh any older
`hg2_corrected/` outputs.

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
output/slurm_output/real_data_ensemble_geometric_n20/52906940/figures/real_data_ensemble/overall/paper_overall_annotation_matched_basis_drift_pve.png
output/slurm_output/real_data_ensemble_geometric_n20/52906940/figures/real_data_ensemble/overall/paper_overall_annotation_matched_basis_drift_pve.pdf
output/slurm_output/real_data_ensemble_geometric_n20/52906940/figures/real_data_ensemble/overall/paper_overall_annotation_matched_basis_drift_pve_panel_b.png
output/slurm_output/real_data_ensemble_geometric_n20/52906940/figures/real_data_ensemble/overall/paper_real_data_annotation_weight_labels.csv
output/slurm_output/real_data_ensemble_geometric_n20/52906940/figures/real_data_ensemble/overall/paper_real_data_matched_basis_drift_summary.csv
output/slurm_output/real_data_ensemble_geometric_n20/52906940/figures/real_data_ensemble/overall/paper_selected_locus_zoom_lollipop_l2_drift.png
output/slurm_output/real_data_ensemble_geometric_n20/52906940/figures/real_data_ensemble/overall/paper_supplement_real_vs_sim_multimodality.png
output/slurm_output/real_data_ensemble_geometric_n20/52906940/figures/real_data_ensemble/overall/paper_real_data_locus_summary_table.csv
output/slurm_output/real_data_ensemble_geometric_n20/52906940/figures/real_data_ensemble/overall/paper_real_data_pip_gt05_totals.csv
```

Copy the paper PNGs into:

```text
../Writings/plots/real_data_case_study/
```

to match the current paper figure paths. The table CSV summaries are written
there directly by the visualization workbook when the writing directory exists. The matched-basis combined figure is also written there directly when possible.

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
- the final paper PNGs/PDFs, including the matched-basis overall figure, if PDFs were rendered.

## 14. Resolved and remaining provenance notes

Resolved for this trace:

1. The real-data study should be described as 20 loci. Earlier 10-locus claims
   are outdated.
2. The production eQTL data-prep source root is
   `/storage/home/mgc5166/work/Annotations/eQTL_annotations_for_susine`.
3. `parent_job_id = 52906940` is confirmed as the run used for the paper
   figures.
4. No selected loci failed. The paper-side summary has 20 locus rows with
   populated anchor, highest-weight-source, and warm-refit run ids.
5. The cluster-weight ambiguity is resolved in the current source: both
   `agg_weight_run` and `aggregated_pip` (in `real_data_pipeline.R` lines
   ~2124-2235) use the credible-shift, Method-B, frequency-free
   `.cluster_weights_from_hc(..., aggregation = "method_b")` rule described in
   Section 10.2. Every fit can receive nonzero weight; the representative run ids
   remain cluster descriptors, not the only contributors.
6. The PVE/heritability estimand switched (commit `4e0721b`) to
   `E[var(Xb)|y]/var(y)` via `R/heritability.R` (`hg2_components()`), recomputed
   downstream for real data by `inst/scripts/recompute_real_data_hg2.R`. Panel C
   and the locus summary table plot the corrected estimand when `hg2_corrected/`
   is present, falling back to the posterior-mean proxy otherwise (Section 9).
7. The supplement real-vs-sim multimodality figure
   (`paper_supplement_real_vs_sim_multimodality.png`) uses the credible-shift
   multimodal metrics `n_clusters_credible` / `max_credible_dist`, not the JSD
   metrics; benchmarked against simulation job `ensemble_scaling_full` /
   `53547760` (Section 11.3).

Remaining or intentionally deferred:

1. AlphaGenome remote model/version is unknown. The source records the
   AlphaGenome client call pattern and scorer class, but does not record a
   remote model identifier in the generated outputs.
2. Exact reproducibility packaging is deferred to the supplemental information
   package. That package should archive generated `locus_manifest.csv`,
   `job_config.json`, `run_manifest.csv`, `metric_inventory.csv`, and relevant
   upstream annotation artifacts for job `52906940` with checksums.
3. **Resolved 2026-06-16 source fix, rerun still required**:
   `recompute_real_data_hg2.R` was migrated from legacy JSD + Method-C weights
   to the collector's credible-shift + Method-B `cluster_weight_credible` rule.
   Regenerate `hg2_corrected/` before quoting the corrected ensemble-PVE numbers
   or rerendering Panel C.
4. The `ensemble_scaling_full` simulation `parent_job_id = 53547760` used as the
   multimodality background (Section 11.3) is set in the visualization workbook
   but was not independently verified against the simulation trace. FLAG for
   Michael to confirm it matches the simulation figures' job id.
5. The "selected" zoom loci for the L2 figure
   (`ydjc, arsa, rrp7a, tmtc1, lgals9, znf280b`) are hardcoded tokens in the
   workbook; which loci are "selected" is an editorial choice, not derivable
   from outputs. Recorded here verbatim from the source, not independently
   justified.
