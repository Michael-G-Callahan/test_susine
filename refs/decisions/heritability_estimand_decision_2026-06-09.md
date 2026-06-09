# Heritability ("genetic causal variance") estimand — findings & fix plan

Date: 2026-06-09
Status: **Full path implemented (2026-06-09).** Shared helper, fit-time plumbing,
persisted columns, realized-truth persistence, paper workbooks (panel E +
decomposition + truth diagnostics), the real-data recompute script, and the
hotfix minibatch run-control are all written and tested/syntax-checked. **What
remains is execution, not code:** run the hotfix SLURM job + the real-data
recompute on the HPC, then regenerate panel E. See "Implementation status".

## 0. TL;DR

- Panel E of the ensemble-scaling composite figure
  (`output/slurm_output/ensemble_scaling_full/53547760/figures/paper_ensemble_scaling/paper_ensemble_scaling_composite_performance.png`)
  showed **SuSiNE (max-ELBO) ≈ SuSiNE (ensemble), both well below SuSiE**, with
  all three under the true `h2 = 0.25`. The question was whether heritability is
  measured apples-to-apples across the three boxes.
- **It is the same estimand in form** (`var` of the **posterior-mean** genetic
  predictor over `var(y)`), so not a units mismatch — but that estimand is
  **structurally biased downward by an amount that scales with posterior
  diffuseness**, which provides a plausible mechanism for SuSiNE (more diffuse
  posteriors under directional priors) reading lower than SuSiE.
  **Caveat (do not overclaim):** the diagnosis is mathematically sound, but how
  much of the observed gap is this artifact versus genuine model shrinkage /
  prior-variance effects is **not established until the corrected metric is
  actually recomputed**. Treat "the gap is an artifact" as a hypothesis to test,
  not a finding.
- **The real-data case studies use the same estimand** (`b'Rb`, the RSS form of
  `var(Xb)/var(y)`), so they inherit the same bias and the same ensemble
  double-shrinkage.
- **Fix**: switch to `E[var(Xβ) | y] / var(y)` (posterior mean of the genetic
  variance = local PVE / cis-`h2_g`), which adds back the dropped
  `tr(Σ·Cov(β|y))` term and is fair across models. For the ensemble, aggregate as
  the **weighted mean of per-fit variances**, not the variance of the weighted
  mean.
- **Rerun implications**: simulation ensemble **must refit** (raw fits/second
  moments are not persisted in buffered mode); real-data case studies **do not**
  (fits are saved to disk → downstream recompute).

## 1. How heritability is currently computed

### 1.1 Simulation (panel E mapping)

Source: `prepare_results_workbook_ensemble_scaling_paper.Rmd:832-882`.

| Box | Source table | Filter | Column |
|-----|--------------|--------|--------|
| SuSiE | `model_metrics` (baseline-single, `susine_vanilla`) | — | `hg2` |
| SuSiNE (max-ELBO fit) | `hg2_by_agg` (spec C-CS, `annotation_r2=0.3`) | `agg_method == "max_elbo"` | `hg2` |
| SuSiNE (ensemble) | `hg2_by_agg` (spec C-CS, `annotation_r2=0.3`) | `agg_method == "cluster_weight_credible"` | `hg2` |

Two code paths feed these:

- **SuSiE box** → `estimate_hg2(fit, y, X)` (`evaluation_metrics.R:280-311`).
  Priority cascade: (1) `var(fit$model_fit$fitted_y)/var(y)`; (2)
  `1 - sigma_2/var(y)`; (3) `var(X %*% coef)/var(y)`. Branch 1 fires in practice.
- **SuSiNE boxes** → `compute_hg2_by_agg(...)` (`run_model.R:2775-2837`). Emits
  `hg2 = var(fy_agg)/var(y)` (weight-aggregated fitted-y) and
  `hg2_weighted_pve = Σ_k w_k·var(fy_k)/var(y)`. Panel E uses the `hg2` column.

The per-fit `fitted_y` is the posterior-mean linear predictor:
susine via `finalize_model_fit` → `compute_Xb(X, std_coef)`
(`../susine/R/finalize.R:82`); the susieR-dispatch baseline via `fit_raw$fitted`
(`run_model.R:1486-1504`). The intercept difference is irrelevant (`var()` is
invariant to additive constants).

**Net:** all three boxes are `var(E[Xβ|y]) / var(y)` — the **point-estimate**
estimand (call it estimand #1).

### 1.2 Real data case studies

Source: `R/real_data_pipeline.R`.

- **Operative metric:** `pve_postmean_std = b' R b`, with
  `b = colSums(alpha ⊙ mu)` the posterior-mean effect vector and `R` the locus
  LD matrix; clamped to `[0,1]`. See `real_data_pve_from_posterior_mean`
  (`:803-813`), per-fit (`:1063-1064`), ensemble via `agg_posterior_mean`
  (`:2271-2272`).
- **Degenerate proxy:** `h2_proxy_std = 1 - sigma_2_final` (`:1087`). Real-data
  runs fix `estimate_residual_variance = FALSE`, `residual_variance = 1`
  (`:1600-1602`), so `sigma_2_final ≈ 1` ⇒ `h2_proxy_std ≈ 0`. **Ignore it.**

In the standardized RSS regime (`var(y)=1`, X column-standardized),
`var(Xb) = b'Rb` with `R = X'X/n`. So **real-data `pve_postmean_std` is exactly
estimand #1** in summary-stat algebra — same quantity as the simulation `hg2`
box, with `R` in place of `X'X/n`.

## 2. The statistical problem

Target ("genetic causal variance") is a fixed functional of β:
`V_g = var(Xβ) = β'Σβ`, `Σ = X'X/n` (individual) or `R` (RSS), reported as a
fraction of `var(y)` = local PVE / cis-`h2_g`.

Law of total variance:

```
E[ V_g | y ]            =     g(E[β|y])              +   tr(Σ · Cov(β|y))
(posterior mean of V_g)   (estimand #1, current)       (dropped correction ≥ 0)
```

- Estimand #1 is the principled quantity **minus a nonnegative term that grows
  with posterior diffuseness**. SuSiNE's directional priors elevate
  annotation-aligned variants → PIP mass spreads → larger `Cov(β|y)` → larger
  dropped term → lower reading. **This is a plausible driver of the SuSiE >
  SuSiNE gap, to be confirmed by recomputing the corrected metric** — it may
  coexist with real shrinkage / prior-variance differences.
- **max-ELBO ≈ ensemble is expected**, not a bug: when the ELBO-dominant fit
  also carries the cluster weight, the weighted `fy_agg` ≈ the single best fit's
  fitted-y, so the BMA point estimate barely shrinks.
- The **ensemble box double-shrinks**: `var(Σ_k w_k E_k[Xβ])` (variance of the
  weighted-mean predictor) loses the between-fit dispersion, which is why it
  collapses onto the max-ELBO box.

## 3. The fix — recommended estimand

Report `E[var(Xβ) | y] / var(y)`, with the clean decomposition:

```
V_g/var(y)  =  g(E[β|y])/var(y)        "confidently localized"
            +  tr(Σ Cov(β|y))/var(y)    "known-present, unlocalized" (within-fit)
            +  between-fit dispersion    (ensemble only)
```

- **Fair**: posterior mean of the *same* scalar under each model; does not
  penalize a more diffuse posterior. Localization is handled separately by
  AUPRC/CS metrics. `V_g` is (desirably) near-invariant to which SNP in an LD
  block carries the weight.
- **Interpretable**: it *is* local PVE / cis-`h2_g`, the object geneticists
  already use.
- **Cheap & exact, using the moments the repo already stores.** In susine,
  `effect_fits$b_hat` / `b_2_hat` are **unconditional per-effect moments**
  (`SER.R:55`: `b_hat = α⊙μ_1`, `b_2_hat = α⊙(σ²_1 + μ_1²) = E[b²]`), *not*
  conditional `μ/μ2`. Because each single effect has exactly one active variable,
  `E[b_l b_l'] = diag(b_2_hat[l,])` (off-diagonals vanish), so the helper is:

  ```r
  m_l <- b_hat[l, ]
  m   <- colSums(b_hat)                      # = posterior-mean coef
  hg2_postmean   <- t(m) %*% Sigma %*% m
  hg2_correction <- sum_l( sum(diag(Sigma) * b_2_hat[l, ]) - t(m_l) %*% Sigma %*% m_l )
  hg2_expected_pve <- hg2_postmean + hg2_correction
  ```

  This relies on the SER single-effect structure (one active variable per effect)
  and the mean-field independence across the `L` effects (`Cov(β) = Σ_l Cov(b_l)`)
  — i.e. it reports the model's own variational posterior covariance. Do **not**
  apply it to a non-SER fit. (It is the same `E[‖Xb‖²]` content that
  `ERSS.R`/`ERSS_ss.R` compute, but those use a `colSums(X^2)` /
  `diag(X'X)` scaling — see the denominator note below before reusing them.)

- **Fair ensemble form** (headline): `V_g^ens / var(y) = Σ_k w_k · hg2_expected_pve_k`,
  the weighted mean of per-fit corrected variances, using the *same* cluster/ELBO
  weights. This already equals the mixture posterior expectation.
  **Do not add a between-fit term on top of it** — that double-counts. If you want
  the *decomposition*, use the identity (PSD quadratic form):

  ```text
  Σ_k w_k m_k'Σm_k  =  m̄'Σm̄  +  Σ_k w_k (m_k − m̄)'Σ(m_k − m̄),   m̄ = Σ_k w_k m_k
  ```

  so

  ```text
  V_g^ens/var(y) =  m̄'Σm̄                                  (ensemble posterior-mean part)
                 +  Σ_k w_k (m_k − m̄)'Σ(m_k − m̄)          (between-fit dispersion)
                 +  Σ_k w_k · hg2_correction_k             (weighted within-fit corrections)
  ```

  The **old ensemble box is only the first term** (`m̄'Σm̄`), dropping *both* the
  between-fit and the within-fit pieces.

Shortcuts considered and rejected as the primary metric: `1 - σ̂²/var(y)`
(inherits residual-estimator choice differences; no clean ensemble definition —
σ̂² of the BMA ≠ BMA of σ̂²). Keep it only as a sanity cross-check. Stay
**in-sample**.

### 3.1 Denominator convention (must match the plotted ratio)

The live metric is `stats::var(fitted_y) / stats::var(y)` — both use the `(n−1)`
denominator, which **cancels in the ratio**. The corrected helper must keep that
consistency: pick `Σ` so that `m'Σm = var(Xm)` under the *same* convention as the
`var(y)` it divides by. Concretely, use **centered, `(n−1)`-scaled** Gram:

```r
Xc    <- scale(X, center = TRUE, scale = FALSE)
Sigma <- crossprod(Xc) / (nrow(X) - 1)        # individual-data sims
denom <- stats::var(y)
```

Do **not** silently use `X'X / n` unless every denominator is adjusted to match.
For **real-data RSS**, `var(y) = 1` and `Σ = R` (LD correlation matrix) already on
the standardized scale, so `b'Rb` is directly a fraction of `var(y)` with no extra
denominator — keep that path as-is, but route it through the *same* helper with
`Σ = R`, `denom = 1`.

## 4. Persistence — what drives rerun needs

| Pipeline | Raw fits / per-fit `b_hat`,`b_2_hat` persisted? | Can patch downstream-only? |
|----------|-------------------------------|----------------------------|
| **Simulation ensemble** | **No.** Buffered I/O (`verbose_file_output=FALSE`); `saveRDS(...fit.rds)` is gated by the verbose branch (`run_model.R:2293-2299, 2430`). Flush files store only scalar `hg2`/`hg2_weighted_pve` (`:3017-3020`). Per-fit second moments are discarded after flush. | **No → must refit.** |
| **Real data** | **Yes.** `saveRDS(fit_result, fit_path, compress="xz")` (`:1723`), reloaded for drift/refit (`:1401, 2357, 2639`). `R` lives in `data/real_case_studies/...`. | **Yes → no rerun.** |

The corrected estimand is computed at **fit time** inside `compute_hg2_by_agg`
(needs per-fit `b_hat`/`b_2_hat` + the Gram matrix), so any future
ensemble/baseline run picks it up automatically. Only **already-run** simulation
jobs need a refit.

## 5. Plan for repo changes

### 5.1 Shared estimand helper (single source of truth)

- New helper parameterized by the **Gram matrix `Σ`** so individual-data and RSS
  pipelines compute a provably identical quantity. Inputs: a fit's `b_hat`,
  `b_2_hat` (L×p unconditional moments), `Σ`, and `denom` (= `var(y)`, or 1 for
  RSS). Returns the column schema in §5.3:
  - `hg2_postmean   = (m'Σm) / denom`, `m = colSums(b_hat)` (current estimand #1),
  - `hg2_uncertainty = (Σ_l [sum_j Σ_jj·b_2_hat[l,j] − m_l'Σm_l]) / denom`
    (within-fit correction),
  - `hg2_expected_pve = hg2_postmean + hg2_uncertainty`.
- `Σ` must match the plotted ratio's denominator (see §3.1): centered
  `crossprod(Xc)/(n−1)` for sims, `R` for RSS. **Do not** blindly reuse
  `ERSS.R`'s `colSums(X^2)` scaling — it uses a different (uncentered, `n`-style)
  convention.
- Proposed location: `test_susine/R/heritability.R` (both pipelines live in
  test_susine). **Open decision: helper home (test_susine vs an exported susine
  accessor).**

### 5.2 Simulation fit-time path

- Extend `compute_hg2_by_agg` to receive per-fit `b_hat`/`b_2_hat` (or the
  precomputed per-fit scalars `hg2_postmean_k`, `hg2_uncertainty_k`), thread the
  bundle `Σ`, and emit: `hg2_postmean`, `hg2_uncertainty`, `hg2_expected_pve`,
  plus the fair ensemble aggregate `Σ_k w_k·hg2_expected_pve_k` and the
  `hg2_between_fit` dispersion term (for decomposition only — **not** added to the
  headline; see §3). Keep the old `hg2` column unchanged for continuity.
- Route the **SuSiE baseline** through the identical helper. **Gotcha:**
  `normalize_susier_fit()` stores only `effect_fits = list(alpha = alpha)`
  (`run_model.R:1497`) and already computes `mu` (`:1475`). Add
  `b_hat = alpha*mu` and `b_2_hat = alpha*mu2` (pull `fit_raw$mu2`) so the baseline
  forms its correction term identically (replaces `estimate_hg2`'s
  `var(fitted_y)` branch for this purpose).
- (Optional, future-proofing) add a `write_hg2_components` flag that flushes the
  per-fit `hg2_postmean`/`hg2_uncertainty` scalars, so future estimand tweaks on
  *this* data become downstream-only and never need another refit.

### 5.3 Persisted column schema (use in both pipelines)

| Column | Meaning |
|--------|---------|
| `hg2` (legacy) | unchanged old `var(fitted_y)/var(y)`; keep for backward comparability |
| `hg2_postmean` | `var(E[Xb \| y])/var(y) = m'Σm/denom` (= old estimand #1, explicit) |
| `hg2_uncertainty` | within-fit posterior second-moment correction |
| `hg2_expected_pve` | `E[var(Xb) \| y]/var(y)` = `hg2_postmean + hg2_uncertainty` — **the headline** |
| `hg2_between_fit` | ensemble mixture-dispersion component (where relevant) |

Naming rationale: avoid `hg2_fair` (motivational, not descriptive) and bare
`hg2_pve` (the old metric is also PVE-like colloquially). New plots use
`hg2_expected_pve`. *(Minor open call: rename `hg2_uncertainty` → `hg2_within_fit`
to parallel `hg2_between_fit`.)*

### 5.3b Downstream — simulation

- `prepare_results_workbook_ensemble_scaling_paper.Rmd`: swap the panel-E column
  → `hg2_expected_pve`; relabel; add the decomposition stacked-bar exhibit
  (`hg2_postmean` / `hg2_within_fit` / `hg2_between_fit`).
- `visualize_results_workbook_ensemble_scaling_paper.Rmd`: render the new panel +
  decomposition; fix the caption to say "posterior mean of local genetic
  variance (PVE)", not a raw `var(fitted_y)`.
- **Truth line (panel E):** keep the single dashed `0.25` design-target line in
  the main figure (calibrated `h2_total = 0.25`, easier to read). Separately
  compute **dataset-specific realized truth** with the *same* denominator
  convention as the estimator:

  ```r
  true_hg2_realized <- as.numeric(var(X %*% beta) / var(y))
  ```

  Report mean / SD / min / max of `true_hg2_realized` in a diagnostics table, and
  add a supplement panel of truth-centered error
  `hg2_expected_pve - true_hg2_realized`. If realized truth is tightly centered
  near 0.25, the main dashed line stands; if it varies noticeably, the error plot
  becomes the stronger scientific comparison.

### 5.4 Hotfix minibatch ensemble run (only the grid we need)

- Build a **dedicated mini job-config** scoped to exactly what panel E consumes:
  **spec C-CS at `annotation_r2 = 0.3`** + the **`susine_vanilla` baseline-single**
  rows. Refit with the patched code under a new `parent_job_id`. Seed-deterministic
  data ⇒ reproduces the same fits and just emits the fair `hg2` columns.
- Re-run collect → prepare → visualize for **panel E only** against the hotfix
  job id. (Confirm this is acceptable for the composite figure, or splice panel E
  back into the full composite.)
- Do **not** rerun the full annotation × spec × lever grid unless we want the
  fair estimand everywhere for completeness.

### 5.5 Downstream — real data (no rerun)

- Add a post-hoc recompute pass (in `collect_results_workbook_real_data_ensemble.Rmd`
  or a small `inst/scripts` driver): reload each saved `fit_rds_path`, pull `R`
  via `load_real_data_locus_bundle`, run the shared helper with `Σ = R`,
  `denom = 1` to get `hg2_postmean`/`hg2_uncertainty`/`hg2_expected_pve` per fit
  and the fair ensemble aggregate, write new columns alongside `pve_postmean_std`
  (which equals `hg2_postmean`). No SLURM job.
- Going forward, also compute the fair estimand inline in
  `run_real_data_task` so new real-data runs carry it natively.

## 6. Open decisions

1. Helper home: `test_susine/R/heritability.R` vs an exported susine accessor.
2. ~~Column naming~~ **Resolved**: `hg2_expected_pve` headline; schema in §5.3.
   Sub-call still open: `hg2_uncertainty` vs `hg2_within_fit`.
3. ~~Truth line~~ **Resolved**: single `0.25` in main panel; dataset-specific
   `true_hg2_realized` + truth-centered error in diagnostics/supplement (§5.3b).
4. Whether to persist per-fit components now (future-proofing) — recommended.
5. Composite-figure handling for the hotfix (panel-E splice vs full re-render).
6. Keep `hg2` (old) columns for backward comparability — recommended yes.

## 6a. Resolved caveats from review (2026-06-09)

- **Do not overclaim "artifact".** Bias direction is proven; the *share* of the
  Panel-E gap due to it is a hypothesis until recomputed (§0, §2).
- **Moment convention.** Use repo-native unconditional `b_hat`/`b_2_hat`
  (verified `SER.R:55`), not conditional `μ/μ2`; single-effect structure gives the
  diagonal `E[b_l b_l'] = diag(b_2_hat[l,])` (§3, §5.1).
- **Gram denominator.** Match the `stats::var` `(n−1)` convention; use centered
  `crossprod(Xc)/(n−1)`, not `X'X/n`, and not `ERSS`'s `colSums(X^2)` scaling
  (§3.1).
- **Ensemble decomposition.** `Σ_k w_k·hg2_expected_pve_k` is the headline and
  already includes between-fit; the `m̄'Σm̄ + between + within` split is an
  equality, not an additive extra (§3).

## 7. Key file/line index

- `R/evaluation_metrics.R:280-311` — `estimate_hg2` (SuSiE box path).
- `R/run_model.R:2775-2837` — `compute_hg2_by_agg` (SuSiNE boxes).
- `R/run_model.R:2293-2299, 2430` — verbose-gated raw-fit save (sim: off).
- `R/run_model.R:3017-3020` — flush writes scalar `hg2_by_agg` only.
- `R/run_model.R:1461-1504` — susieR→susine adaptor (stores only `alpha`).
- `../susine/R/finalize.R:82` — `fitted_y = compute_Xb(X, std_coef)`.
- `R/real_data_pipeline.R:803-813` — `real_data_pve_from_posterior_mean` (`b'Rb`).
- `R/real_data_pipeline.R:1063-1064, 2271-2272` — per-fit & ensemble PVE.
- `R/real_data_pipeline.R:1087, 1600-1602` — degenerate `h2_proxy_std`.
- `R/real_data_pipeline.R:1723` — fits persisted (`saveRDS ... compress="xz"`).
- `prepare_results_workbook_ensemble_scaling_paper.Rmd:832-882` — panel E data.

## 8. Implementation status (2026-06-09)

**Done (committed to working tree, load_all + tests green):**

- `R/heritability.R` (NEW) — shared engine `hg2_components_from_moments()`,
  dispatcher `hg2_components(fit, X=/R=, y=, vy=)`, light
  `hg2_uncertainty_scalar()`. Gram-parameterized (X for sims, R for RSS),
  matches `stats::var` `(n−1)` via `var(X %*% v)`; components floored at 0,
  `expected_pve` left un-capped to preserve exact additivity.
- `R/run_model.R` — `normalize_susier_fit()` now stores
  `b_hat = alpha*mu`, `b_2_hat = alpha*mu2`; per-fit `hg2_uncertainty_scalar`
  collected into `hg2_unc_by_group` at the fitted-y site and threaded into
  `compute_hg2_by_agg` (new `uncertainty_list` arg). `compute_hg2_by_agg` now
  emits `hg2_postmean`, `hg2_uncertainty`, `hg2_between_fit`, `hg2_expected_pve`
  (legacy `hg2`, `hg2_weighted_pve` unchanged). Between-fit computed directly as
  `Σ_k w_k var(fy_k − fy_agg)/vy`; additive identity verified.
- `R/evaluation_metrics.R` — `evaluate_model` adds `hg2_postmean`,
  `hg2_uncertainty`, `hg2_expected_pve` to model_metrics (both filtered +
  unfiltered) via the shared helper → SuSiE baseline box gets the fair value.
- `tests/testthat/test-heritability.R` (NEW) — 14 assertions: additivity,
  `postmean == var(fy)/var(y)`, RSS≡standardized-individual, NA-safety, the
  `compute_hg2_by_agg` decomposition + max-ELBO between-fit==0 + NA propagation.
- New columns flow through unchanged: `write_flush_outputs` binds the tibble,
  `collect_results.R:425` binds+enriches+writes `hg2_by_agg.csv`, model_metrics
  carried as-is. No schema plumbing needed.

Resolved open decisions: helper home = `test_susine/R/heritability.R` (#1);
per-fit components persisted now (#4, via the new columns). Naming kept as
`hg2_uncertainty` (the `hg2_within_fit` rename sub-call left to MGC).

**Done — path-wide (2026-06-09, second pass):**

- `R/run_model.R` — `execute_dataset_bundle` now persists
  `true_hg2_realized = var(X beta)/var(y)` per dataset into `dataset_metrics`
  (flows through `write_dataset_metrics` → collect → consolidated, verified).
- §5.3b sim paper workbooks **done**:
  - `prepare_results_workbook_ensemble_scaling_paper.Rmd` — loads
    `dataset_metrics.csv`; hg2 chunk now carries
    `hg2_postmean`/`hg2_uncertainty`/`hg2_between_fit`/`hg2_expected_pve` per
    series (graceful legacy fallback when corrected columns absent), builds
    `hg2_decomp_data` + `hg2_truth_diag`, merges realized truth; all added to the
    saved `paper_plot_data$data`.
  - `visualize_results_workbook_ensemble_scaling_paper.Rmd` — panel E plots
    `hg2_expected_pve` (y-label = local genetic variance E[var(Xb)|y]/var(y),
    dashed 0.25 kept); new §7a decomposition stacked-bar; new §7b realized-truth
    diagnostics table + truth-centered error supplement; comparison chunk uses
    `hg2_expected_pve`. (Pre-existing MD025 H1 lint warnings unrelated.)
- §5.5 real-data recompute **done**: `inst/scripts/recompute_real_data_hg2.R`
  (reloads saved fits, `Σ = R`, `vy = 1`, per-run components + per-locus fair
  ensemble via `.cluster_weights_from_hc`). HPC-run, sanity-check after.
- §5.4 hotfix run-control **done**:
  `run_control_workbook_ensemble_scaling_hg2_hotfix.Rmd`
  (`job_name = ensemble_scaling_hg2_hotfix`; C-CS `c_grid_8 × sigma_grid_8` +
  baseline-single `susine_vanilla`; `annotation_r2 = 0.3`; same M1-stratified
  matrices/seeds; scaling off; aggregations `max_elbo` + `cluster_weight_credible`).

**Remaining — execution only (MGC, on HPC):**

1. Knit the hotfix run-control → `write_job_artifacts` → `sbatch`.
2. Collect the hotfix job (`collect_results_workbook_ensemble_scaling.Rmd`, point
   `parent_job_id` at it).
3. Prepare+visualize paper workbooks pointed at the hotfix job for panel E;
   splice panel E (+ decomposition/diagnostics) into the full composite.
4. Run `inst/scripts/recompute_real_data_hg2.R` on the real-data job; sanity-check.
5. Run `devtools::document()` before any formal `R CMD check` (new internal
   roxygen blocks have no `.Rd` yet; not needed for load_all/runtime).

**Note:** pre-existing unrelated failure `test-ensemble-scaling-metrics.R:60`
(`aggregate_use_case_pips` method naming, the cluster-weight rename WIP) — present
on clean `HEAD`, not caused by this change.

## 8a. Review-pass-2 robustness fixes (2026-06-09)

From a second code review; all in tests (18/18) + integration (23/23):

- **Sim aggregation respected non-default softmax temperature.** Threaded
  `softmax_temperature` through `compute_hg2_by_agg()` (signature + the
  `softmax_weights` / `.cluster_weights_from_hc` calls). **Config-path correction
  (review pass 3):** the sim schema nests this under `job$compute$softmax_temperature`
  (`make_job_config` → `compute_list`, `run_controls.R:1081-1091`; read that way at
  `run_model.R:1144,1349`), NOT the flat `job$softmax_temperature` the first patch
  used — which would have silently stayed at 1. Call site now reads
  `job$compute$softmax_temperature %||% job$softmax_temperature %||% 1` (the flat
  fallback covers the real-data schema, which *is* flat — `real_data_pipeline.R:1973`).
  Regression test asserts temperature actually *moves* the `elbo_softmax`
  expected-PVE/postmean/uncertainty while leaving `max_elbo` invariant.
  **Config round-trip test (review pass 3):** builds a real `make_job_config(...,
  softmax_temperature = 0.5)`, asserts it nests at `job$compute`, that the
  call-site resolver returns 0.5 both in-memory and after the JSON write/read the
  SLURM task uses, and that a flat-only read returns 1 — so reverting the call
  site to the flat path is now caught by a test (closes the exact gap the direct
  `compute_hg2_by_agg(..., temperature=)` test missed). 25/25 in
  `test-heritability.R`.
- **Zero-weight NA no longer poisons the within-fit term.** `0 * NA = NA` meant a
  discarded (zero-weight) fit lacking moments could NA-out e.g. `max_elbo`'s
  uncertainty. Now restricted to `w > 0` fits; NA still propagates correctly when
  a *positively*-weighted fit is missing moments. (regression test added.)
- **Real-data recompute no longer silently drops NA positive-weight fits.**
  `inst/scripts/recompute_real_data_hg2.R` replaced `sum(..., na.rm=TRUE)` with
  `weighted_or_na()` → returns NA (not a renormalized partial sum / 0) when any
  positive-weight component is missing; records `n_missing_components` +
  `hg2_expected_pve_pos_{min,max}` for bounds checking.
- **Corrected the recompute sanity note:** the weighted ensemble is a convex
  combination of positively-weighted per-fit expected PVEs, so it lies within
  `[pos_min, pos_max]` and may be **below** the max-ELBO member — not "at or
  above" it (the old note was wrong).
