# tmp_prove_cache_poisoning.R
# ---------------------------------------------------------------------------
# 100%-confidence, DATA-ONLY proof of the cache-write poisoning diagnosis.
# No model runs, no corrective action. Reads only recorded outputs.
#
# Claim to prove:
#   (1) IDENTITY  : ensemble `baseline-single` metrics == the metrics of the run
#                   it cache-sourced from (cache_source_run_id), to machine eps.
#   (2) SOURCE    : those source runs are NON-ROOT refine refits
#                   (explore_method == "refine", variant_id >= 2).
#   (3) WRONG FIT : ensemble `baseline-single` != baseline-screen cold
#                   susine_vanilla (the un-poisoned cold baseline).
#
#   (1)+(2) => the ensemble "baseline" IS a refined fit.
#   (3)     => it is NOT the cold baseline it claims to be.
#
# Run from the test_susine repo root on ROAR.
# ---------------------------------------------------------------------------

suppressPackageStartupMessages({ library(dplyr); library(readr) })

ens_dir  <- "output/slurm_output/ensemble_scaling_full/53398149/consolidated"
base_dir <- "output/slurm_output/baseline_sims_screen/51509956/consolidated"

e_mm <- read_csv(file.path(ens_dir,  "model_metrics_full.csv"), show_col_types = FALSE, progress = FALSE)
b_mm <- read_csv(file.path(base_dir, "model_metrics_full.csv"), show_col_types = FALSE, progress = FALSE)

# Real metric columns (verified from the prior schema dump)
metrics <- intersect(
  c("AUPRC", "power", "hg2", "elbo_final", "cross_entropy", "mean_size",
    "mean_purity", "mean_coverage", "tpr_at_fpr_05"),
  names(e_mm)
)
sig_eq02 <- function(x) abs(suppressWarnings(as.numeric(x)) - 0.2) < 1e-9

# Ensemble baseline-single rows (susine_vanilla, sigma=0.2)
bs <- e_mm %>%
  filter(spec_name == "baseline-single", use_case_id == "susine_vanilla",
         sig_eq02(sigma_0_2_scalar))

cat("================ (1) IDENTITY: baseline-single == its cache source ================\n")
cat(sprintf("baseline-single rows: %d  (distinct bundles: %d)\n",
            nrow(bs), n_distinct(bs$dataset_bundle_id)))

# Pull the cache-SOURCE run's own recorded metrics by joining run_id.
# (One primary recorded row per run_id; take the model_call_executed==TRUE base row.)
src_pool <- e_mm %>%
  filter(is.na(agg_method) | agg_method == "") %>%
  select(run_id, spec_name, explore_method, variant_id, model_call_executed,
         all_of(metrics)) %>%
  distinct(run_id, .keep_all = TRUE)

bs_src <- bs %>%
  select(dataset_bundle_id, run_id, cache_source_run_id,
         all_of(metrics)) %>%
  left_join(src_pool, by = c("cache_source_run_id" = "run_id"),
            suffix = c("", ".src"))

cat(sprintf("baseline-single rows joined to a cache source: %d / %d\n",
            sum(!is.na(bs_src$spec_name)), nrow(bs_src)))
cat("\nPer-metric |baseline-single - cache_source| (should be ~0 if identical):\n")
for (m in metrics) {
  d <- abs(bs_src[[m]] - bs_src[[paste0(m, ".src")]])
  cat(sprintf("    %-16s max|diff|=%.3e   n_exact_equal=%d/%d\n",
              m, max(d, na.rm = TRUE), sum(d < 1e-12, na.rm = TRUE), sum(!is.na(d))))
}

cat("\n================ (2) SOURCE: what kind of run is the cache source? ================\n")
bs_src %>%
  mutate(is_root = (variant_id == 1L)) %>%
  count(src_spec = spec_name, src_explore = explore_method, is_root,
        model_call_executed, sort = TRUE) %>%
  print(n = 50)
cat(sprintf(
  "\nFraction of baseline-single rows sourced from a NON-ROOT refine refit (explore=='refine' & variant_id>=2): %.1f%%\n",
  100 * mean(bs_src$explore_method == "refine" & bs_src$variant_id >= 2L, na.rm = TRUE)))

cat("\n================ (3) WRONG FIT: baseline-single vs baseline-screen COLD ================\n")
# baseline-screen cold susine_vanilla sigma=0.2 (this job has NO refine specs -> un-poisoned)
bc <- b_mm %>%
  filter(spec_name == "susine_vanilla_sigma", use_case_id == "susine_vanilla",
         sig_eq02(sigma_0_2_scalar))

one_per_bundle <- function(df) df %>% group_by(dataset_bundle_id) %>% slice(1) %>% ungroup()
bs1 <- one_per_bundle(bs); bc1 <- one_per_bundle(bc)

j <- inner_join(
  bs1 %>% select(dataset_bundle_id, all_of(metrics)),
  bc1 %>% select(dataset_bundle_id, all_of(metrics)),
  by = "dataset_bundle_id", suffix = c(".ens", ".base"))
cat(sprintf("matched bundles: %d\n", nrow(j)))
for (m in metrics) {
  d <- abs(j[[paste0(m, ".ens")]] - j[[paste0(m, ".base")]])
  cat(sprintf("    %-16s max|diff|=%.3e   n_exact_equal=%d/%d\n",
              m, max(d, na.rm = TRUE), sum(d < 1e-12, na.rm = TRUE), nrow(j)))
}

cat("\n================ VERDICT KEY ================\n")
cat("If (1) all max|diff|~0, (2) ~100% non-root refine, (3) AUPRC differs on many bundles:\n")
cat("  => ensemble 'baseline-single' IS a refined fit served from the poisoned cold cache key,\n")
cat("     and is NOT the cold baseline. Diagnosis confirmed.\n")
cat("================ DONE ================\n")
