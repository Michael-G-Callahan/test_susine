# tmp_diagnose_ensemble_baseline.R
# ---------------------------------------------------------------------------
# Localize WHY the ensemble-internal "baseline-single" susine_vanilla sigma=0.2
# fit differs from the independent baseline-screen susine_vanilla sigma=0.2 fit.
#
# Logic: a cold susine_vanilla sigma=0.2 fit is deterministic (no restarts) and
# annotation-independent, and X,y depend only on phenotype_seed+architecture
# (identical across both jobs). So per-bundle model metrics (AUPRC etc.) of the
# true cold fit MUST match exactly across sources. Where they DON'T match tells
# us which layer is responsible.
#
# Run from the test_susine repo root on ROAR (paths are relative to output/).
# Prints only schemas, distinct provenance patterns, and max-abs diffs.
# ---------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
})

# ── Paths ──────────────────────────────────────────────────────────────────
ens_dir  <- "output/slurm_output/ensemble_scaling_full/53398149/consolidated"
base_dir <- "output/slurm_output/baseline_sims_screen/51509956/consolidated"
ens_hist <- "output/run_history/ensemble_scaling_full/53398149"

rd <- function(dir, f) {
  p <- file.path(dir, f)
  if (!file.exists(p)) { message("MISSING: ", p); return(NULL) }
  readr::read_csv(p, show_col_types = FALSE, progress = FALSE)
}

e_mm  <- rd(ens_dir,  "model_metrics_full.csv")
b_mm  <- rd(base_dir, "model_metrics_full.csv")
e_cb  <- rd(ens_dir,  "confusion_bins_full.csv")
b_cb  <- rd(base_dir, "confusion_bins_full.csv")
# run_table location can vary; try a couple of names
e_rt  <- rd(ens_hist, "run_table.csv")
if (is.null(e_rt)) e_rt <- rd(ens_dir, "run_table.csv")

cat("\n================ SCHEMAS ================\n")
cat("model_metrics cols:\n"); print(names(e_mm))
cat("\nconfusion_bins cols:\n"); print(names(e_cb))
cat("\nrun_table cols:\n"); if (!is.null(e_rt)) print(names(e_rt)) else cat("(no run_table)\n")

# ── Flexible column helpers ──────────────────────────────────────────────────
pick <- function(df, cands) cands[cands %in% names(df)][1]
sig_num <- function(x) suppressWarnings(as.numeric(as.character(x)))

mm_metric_cols <- function(df) {
  cand <- c("auprc", "auprc_avg_precision", "average_precision",
            "power_at_fdr", "tpr_at_fdr05", "tpr05", "power", "fdr",
            "coverage", "cs_power", "n_cs", "mean_cs_size")
  intersect(cand, names(df))
}
bundle_col <- function(df) pick(df, c("dataset_bundle_id", "bundle_id"))
sigma_col  <- function(df) pick(df, c("sigma_0_2_scalar", "sigma_0_2", "sigma_scalar"))
agg_col    <- function(df) pick(df, c("agg_method", "aggregation_method"))
explore_col<- function(df) pick(df, c("explore_method", "exploration_method", "run_type"))

# ── Filter to cold vanilla sigma=0.2 base rows for a given source ────────────
cold_vanilla <- function(df, spec, use_case = "susine_vanilla") {
  bc <- bundle_col(df); sc <- sigma_col(df); ac <- agg_col(df); ec <- explore_col(df)
  out <- df %>% dplyr::filter(.data$spec_name == spec,
                              .data$use_case_id == use_case)
  if (!is.na(sc)) out <- out %>% dplyr::filter(abs(sig_num(.data[[sc]]) - 0.2) < 1e-9)
  # base (non-aggregate) rows only
  if (!is.na(ac)) out <- out %>% dplyr::filter(is.na(.data[[ac]]) | .data[[ac]] == "" )
  if (!is.na(ec)) out <- out %>% dplyr::filter(is.na(.data[[ec]]) | .data[[ec]] %in% c("base","default","single"))
  out
}

# Reduce to one row per bundle; warn if duplicates disagree on metrics.
one_per_bundle <- function(df, tag) {
  bc <- bundle_col(df)
  mcols <- mm_metric_cols(df)
  if (is.na(bc) || !nrow(df)) { cat(sprintf("[%s] no rows / no bundle col\n", tag)); return(df) }
  dup <- df %>% group_by(.data[[bc]]) %>% summarise(n = n(), .groups = "drop")
  cat(sprintf("[%s] rows=%d  distinct bundles=%d  max copies/bundle=%d\n",
              tag, nrow(df), nrow(dup), max(dup$n)))
  if (length(mcols)) {
    disagree <- df %>% group_by(.data[[bc]]) %>%
      summarise(across(all_of(mcols), ~ diff(range(.x, na.rm = TRUE))), .groups = "drop") %>%
      summarise(across(all_of(mcols), ~ max(.x, na.rm = TRUE)))
    cat(sprintf("[%s] max within-bundle metric spread (should be ~0 if duplicates are identical):\n", tag))
    print(disagree)
  }
  df %>% group_by(.data[[bc]]) %>% slice(1) %>% ungroup()
}

# ── Compare two cold-fit sources per bundle ──────────────────────────────────
compare_sources <- function(x, y, tag) {
  bcx <- bundle_col(x); bcy <- bundle_col(y)
  mcols <- intersect(mm_metric_cols(x), mm_metric_cols(y))
  if (is.na(bcx) || is.na(bcy) || !length(mcols)) { cat(sprintf("[%s] cannot compare\n", tag)); return(invisible()) }
  j <- inner_join(
    x %>% select(all_of(c(bcx, mcols))) %>% rename(.bundle = all_of(bcx)),
    y %>% select(all_of(c(bcy, mcols))) %>% rename(.bundle = all_of(bcy)),
    by = ".bundle", suffix = c(".x", ".y")
  )
  cat(sprintf("\n[%s] matched bundles: %d\n", tag, nrow(j)))
  for (m in mcols) {
    dx <- j[[paste0(m, ".x")]]; dy <- j[[paste0(m, ".y")]]
    d  <- abs(dx - dy)
    cat(sprintf("    %-22s max|diff|=%.3e   n_exact_equal=%d/%d\n",
                m, max(d, na.rm = TRUE), sum(d < 1e-9, na.rm = TRUE), length(d)))
  }
}

cat("\n================ COLD-FIT FINGERPRINT (per-bundle model metrics) ================\n")
ens_base <- one_per_bundle(cold_vanilla(e_mm, "baseline-single"), "ENS baseline-single (susine)")
ens_as   <- one_per_bundle(cold_vanilla(e_mm, "A-S"),             "ENS A-S sigma=0.2 (susine)")
bse_van  <- one_per_bundle(cold_vanilla(b_mm, "susine_vanilla_sigma"), "BASE susine_vanilla (susine)")
bse_sr   <- one_per_bundle(cold_vanilla(b_mm, "susie_vanilla_sigma", use_case = "susie_vanilla"),
                           "BASE susie_vanilla (susieR)")

compare_sources(ens_base, ens_as,  "A : ENS baseline-single   vs ENS A-S(0.2)        [within-job, susine]")
compare_sources(ens_base, bse_van, "B : ENS baseline-single   vs BASE susine_vanilla [cross-job, same backend]")
compare_sources(ens_as,   bse_van, "B': ENS A-S(0.2)          vs BASE susine_vanilla [cross-job, same backend]")
compare_sources(ens_base, bse_sr,  "D : ENS baseline-single   vs BASE susie_vanilla  [BACKEND gap: susine vs susieR]")
compare_sources(bse_van,  bse_sr,  "D': BASE susine_vanilla   vs BASE susie_vanilla  [BACKEND gap within baseline job]")

# ── Provenance of ENS baseline-single rows (what the cache actually served) ──
cat("\n================ C: ENS baseline-single PROVENANCE ================\n")
prov_cols <- intersect(c("model_call_executed", "cache_source_run_id"), names(e_mm))
if (length(prov_cols) == 2 && !is.null(e_rt)) {
  src <- cold_vanilla(e_mm, "baseline-single") %>%
    distinct(model_call_executed, cache_source_run_id)
  rt_keep <- intersect(c("run_id","spec_name","use_case_id","exploration_methods",
                         "run_type","restart_id","refine_step","c_value",
                         "sigma_0_2_scalar","warm_method","alpha_concentration"),
                       names(e_rt))
  joined <- src %>% left_join(e_rt %>% select(all_of(rt_keep)),
                              by = c("cache_source_run_id" = "run_id"))
  cat("Distinct provenance patterns of baseline-single source runs:\n")
  joined %>% count(across(any_of(setdiff(rt_keep, "run_id"))), model_call_executed,
                   sort = TRUE) %>% print(n = 50)
} else {
  cat("Provenance columns or run_table unavailable. Have: ", paste(prov_cols, collapse=", "), "\n")
}

# ── Confusion-bin calibration: ENS baseline-single vs BASE-SCREEN ────────────
cat("\n================ CALIBRATION (pooled observed freq per PIP bin) ================\n")
cb_cold <- function(df, spec, use_case = "susine_vanilla") {
  if (is.null(df)) return(NULL)
  sc <- sigma_col(df); ac <- agg_col(df)
  out <- df %>% dplyr::filter(.data$spec_name == spec, .data$use_case_id == use_case)
  if (!is.na(sc)) out <- out %>% dplyr::filter(abs(sig_num(.data[[sc]]) - 0.2) < 1e-9)
  if (!is.na(ac)) out <- out %>% dplyr::filter(is.na(.data[[ac]]) | .data[[ac]] == "")
  out
}
calib <- function(df, tag) {
  if (is.null(df) || !nrow(df)) { cat(sprintf("[%s] no bins\n", tag)); return(NULL) }
  binc <- pick(df, c("pip_bin","bin_mid","bin_lower","bin","pip_break"))
  totc <- pick(df, c("n_total","n","count","total"))
  posc <- pick(df, c("n_causal","n_true","n_pos","causal"))
  if (any(is.na(c(binc,totc,posc)))) { cat(sprintf("[%s] bin/total/pos col not found; cols=%s\n", tag, paste(names(df),collapse=","))); return(NULL) }
  df %>% group_by(.bin = .data[[binc]]) %>%
    summarise(tot = sum(.data[[totc]], na.rm=TRUE),
              pos = sum(.data[[posc]], na.rm=TRUE), .groups="drop") %>%
    mutate(obs = pos/tot, src = tag)
}
ce  <- calib(cb_cold(e_cb, "baseline-single"),                                  "ENS_susine")
cb  <- calib(cb_cold(b_cb, "susine_vanilla_sigma"),                             "BASE_susine")
csr <- calib(cb_cold(b_cb, "susie_vanilla_sigma", use_case = "susie_vanilla"),  "BASE_susieR(Fig4)")
if (!is.null(ce) && !is.null(cb)) {
  cmp <- ce %>% select(.bin, obs_ens_susine = obs, tot_ens = tot) %>%
    full_join(cb  %>% select(.bin, obs_base_susine = obs, tot_base = tot), by = ".bin")
  if (!is.null(csr)) {
    cmp <- cmp %>% full_join(csr %>% select(.bin, obs_base_susieR = obs, tot_sr = tot), by = ".bin")
  }
  cmp <- cmp %>% arrange(.bin)
  cat("Per-bin observed causal frequency. KEY comparisons:\n")
  cat("  obs_ens_susine  vs obs_base_susine  -> cross-job, SAME backend (should match if no bug)\n")
  cat("  obs_ens_susine  vs obs_base_susieR  -> backend gap (this is what 'poor vs healthy' likely is)\n")
  print(cmp, n = 100)
  cat(sprintf("\nMax |ens_susine - base_susine| (same backend, cross-job): %.4f\n",
              max(abs(cmp$obs_ens_susine - cmp$obs_base_susine), na.rm = TRUE)))
  if ("obs_base_susieR" %in% names(cmp)) {
    cat(sprintf("Max |ens_susine - base_susieR| (BACKEND gap, Fig4): %.4f\n",
                max(abs(cmp$obs_ens_susine - cmp$obs_base_susieR), na.rm = TRUE)))
  }
}

cat("\n================ DONE ================\n")
