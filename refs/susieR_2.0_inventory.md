# stephenslab/susieR 2.0 — Detailed Inventory

> **Generated**: 2025-07-17
> **Source**: https://github.com/stephenslab/susieR (main branch)
> **Release**: susieR 2.0.0, date 2025-09-22
> **Target journal**: Genome Biology
> **Author of rewrite**: Alex McCreight (vignette attribution)

---

## 1. R/ Source File Listing

| File | Role | Approx. lines |
|------|------|:---:|
| `susie.R` | User-facing interface: `susie()`, `susie_ss()`, `susie_rss()` | ~700 |
| `susie_workhorse.R` | Main orchestration: initialize → iterate → finalize | ~65 |
| `iterative_bayesian_stepwise_selection.R` | IBSS core: `ibss_initialize()`, `ibss_fit()`, `ibss_finalize()` | med |
| `single_effect_regression.R` | SER: `single_effect_regression()`, `optimize_prior_variance()`, `single_effect_update()` | med |
| `susie_constructors.R` | Constructor pattern: `individual_data_constructor()`, `sufficient_stats_constructor()`, `summary_stats_constructor()`, `rss_lambda_constructor()` | ~850 |
| `generic_methods.R` | S3 generics (~20+): `configure_data`, `compute_residuals`, `compute_ser_statistics`, `SER_posterior_e_loglik`, `calculate_posterior_moments`, `update_fitted_values`, `update_variance_components`, etc. | med |
| `individual_data_methods.R` | Backend methods for class `individual` (X, y) | ~250 |
| `sufficient_stats_methods.R` | Backend methods for class `ss` (XtX, Xty, yty). **Also houses SuSiE-ash masking logic and SuSiE-inf omega-weighted residuals** | ~885 |
| `rss_lambda_methods.R` | Backend methods for class `rss_lambda` (z, R, lambda > 0) | ~300 |
| `susie_utils.R` | Internal utilities: initialization, LBF computation, variance estimation (MoM/MLE/NIG), convergence, ELBO, credible sets | ~1,450 |
| `susie_rss_utils.R` | RSS-specific utils: `compute_suff_stat()`, `estimate_s_rss()`, `compute_Dinv()`, stochastic LD precomputation, allele-switch detection | med |
| `mr.ash.R` | Full Mr.ASH implementation (VEB coordinate-wise); uses `caisa_rcpp` C++ backend | ~350 |
| `refinement.R` | `run_refine()` — two-step iterative procedure to escape local optima | small |
| `susie_get_functions.R` | Accessors: `susie_get_cs()`, `susie_get_pip()`, `susie_get_posterior_mean()`, `susie_get_lfsr()`, `susie_get_posterior_samples()`, etc.  | ~500 |
| `susie_auto.R` | `susie_auto()` — automatic L selection by doubling | ~120 |
| `susie_plot.R` | `susie_plot()`, `susie_plot_iteration()`, `susie_plot_changepoint()` | med |
| `predict.susie.R` | `coef.susie()`, `predict.susie()` — includes theta for unmappable effects | small |
| `summary.susie.R` | `summary.susie()`, `print.summary.susie()` | small |
| `diagnosis_reports.R` | `diagnose_susie_ash_iter()` — diagnostic reporting for SuSiE-ash iteration details | ~60 |
| `RcppExports.R` | Auto-generated: `caisa_rcpp()`, `random_order()` | small |
| `susieR-package.R` | Package-level documentation stub | tiny |
| `example_dataset.R` | Dataset documentation (N2finemapping, unmappable_data, etc.) | small |

**C++ source** (`src/`):
- `caisa_rcpp` — coordinate ascent for Mr.ASH (17 parameters, Armadillo-based)
- `random_order` — generates random coordinate update orders

---

## 2. DESCRIPTION / Dependencies

| Field | Value |
|-------|-------|
| Version | 2.0.0 |
| Compiled code | Yes (Rcpp + RcppArmadillo) |
| **Imports** | Rcpp, RcppArmadillo, matrixStats, Matrix, reshape, stats, utils |
| **Suggests** | ggplot2, Rfast (optional), testthat, knitr, rmarkdown |
| **LinkingTo** | Rcpp, RcppArmadillo |

---

## 3. New Features in 2.0 vs Original susieR (≤ v0.12.x)

### From `announcements.Rmd` release notes:

| Feature | Description | Status in v0.12.x |
|---------|-------------|:--:|
| **SuSiE-ash** | Adaptive shrinkage for unmappable effects via Mr.ASH (Stephens 2017). Handles oligogenic/complex backgrounds through flexible unimodal mixture-of-normals prior on residual effects. | NEW |
| **SuSiE-inf** | Infinitesimal effects model (Cui et al. 2024, Nature Genetics). Models polygenic background as $\tau^2 I$ alongside sparse signals. | NEW |
| **Servin-Stephens NIG prior** | Normal-Inverse-Gamma prior that integrates out $\sigma^2$ analytically; produces $t$-distributed marginals. Improves coverage when $n < 80$. | NEW |
| **PIP convergence** | `convergence_method = "pip"` — convergence based on PIP stability instead of ELBO. Forced when unmappable effects are active. | NEW |
| **MoM residual variance** | `estimate_residual_method = "MoM"` — Method of Moments estimator for $\sigma^2$ (more stable than MLE for SuSiE-inf). | NEW (MoM variant) |
| **Stochastic LD correction** | `stochastic_ld_sample = B` in `susie_rss()` — SER variance inflation correction for random-projection LD sketches. | NEW |
| **Improved refinement** | `run_refine()` — two-step iterative CS escape procedure. Incompatible with unmappable effects. | Enhanced |
| **S3 method dispatch** | Unified architecture: 3 data classes (`individual`, `ss`, `rss_lambda`) with polymorphic backends. | NEW (was monolithic) |
| **Constructor pattern** | `Interface → Constructor → Workhorse → IBSS → Backend` pipeline with immutable Data/Params + mutable Model. | NEW |
| **>1,000 tests** | Comprehensive test suite with reference tests against v1.0 outputs for regression testing. | NEW at this scale |

---

## 4. `susie()` / `susie_ss()` / `susie_rss()` — Complete Parameter Signatures

### `susie(X, y, ...)`

```r
susie(X, y,
      L = min(10, ncol(X)),
      scaled_prior_variance = 0.2,
      residual_variance = NULL,
      prior_weights = NULL,
      null_weight = 0,
      standardize = TRUE,
      intercept = TRUE,
      estimate_residual_variance = TRUE,
      estimate_residual_method = c("MoM", "MLE", "Servin_Stephens"),
      estimate_prior_variance = TRUE,
      estimate_prior_method = c("optim", "EM", "simple"),
      unmappable_effects = c("none", "inf", "ash"),
      check_null_threshold = 0,
      prior_tol = 1e-9,
      residual_variance_upperbound = Inf,
      model_init = NULL,
      coverage = 0.95,
      min_abs_corr = 0.5,
      compute_univariate_zscore = FALSE,
      na.rm = FALSE,
      max_iter = 100,
      tol = 1e-3,
      convergence_method = c("elbo", "pip"),
      verbose = FALSE,
      track_fit = FALSE,
      residual_variance_lowerbound = NULL,
      refine = FALSE,
      n_purity = 100,
      alpha0 = 0.1,
      beta0 = 0.1)
```

### `susie_ss(XtX, Xty, yty, n, ...)`

```r
susie_ss(XtX, Xty, yty, n,
         L = min(10, ncol(XtX)),
         X_colmeans = NA, y_mean = NA,
         maf = NULL, maf_thresh = 0,
         check_input = FALSE,
         r_tol = 1e-8,
         standardize = TRUE,
         scaled_prior_variance = 0.2,
         residual_variance = NULL,
         prior_weights = NULL,
         null_weight = 0,
         model_init = NULL,
         estimate_residual_variance = TRUE,
         estimate_residual_method = c("MoM", "MLE", "Servin_Stephens"),
         residual_variance_lowerbound = 0,
         residual_variance_upperbound = Inf,
         estimate_prior_variance = TRUE,
         estimate_prior_method = c("optim", "EM", "simple"),
         unmappable_effects = c("none", "inf"),          # NO "ash"
         check_null_threshold = 0,
         prior_tol = 1e-9,
         max_iter = 100,
         tol = 1e-3,
         convergence_method = c("elbo", "pip"),
         coverage = 0.95,
         min_abs_corr = 0.5,
         n_purity = 100,
         verbose = FALSE,
         track_fit = FALSE,
         check_prior = FALSE,
         refine = FALSE,
         alpha0 = 0.1,
         beta0 = 0.1)
```

### `susie_rss(z, R, n, ...)`

```r
susie_rss(z = NULL, R = NULL, n = NULL,
          X = NULL,
          bhat = NULL, shat = NULL, var_y = NULL,
          L = min(10, ...),
          lambda = 0,
          maf = NULL,
          maf_thresh = 0,
          z_ld_weight = 0,
          prior_variance = 50,
          scaled_prior_variance = 0.2,
          residual_variance = NULL,
          prior_weights = NULL,
          null_weight = 0,
          standardize = TRUE,
          intercept_value = 0,
          estimate_residual_variance = FALSE,        # default FALSE for RSS
          estimate_residual_method = c("MoM", "MLE", "Servin_Stephens"),
          estimate_prior_variance = TRUE,
          estimate_prior_method = c("optim", "EM", "simple"),
          unmappable_effects = c("none", "inf"),     # NO "ash"
          check_null_threshold = 0,
          prior_tol = 1e-9,
          residual_variance_lowerbound = 0,
          residual_variance_upperbound = Inf,
          model_init = NULL,
          coverage = 0.95,
          min_abs_corr = 0.5,
          max_iter = 100,
          tol = 1e-3,
          convergence_method = c("elbo", "pip"),
          verbose = FALSE,
          track_fit = FALSE,
          check_R = FALSE,
          check_z = FALSE,
          stochastic_ld_sample = NULL,               # NEW: sketch size B
          refine = FALSE,
          n_purity = 100,
          alpha0 = 0.1,
          beta0 = 0.1)
```

**Key new parameters** (not in v0.12.x):

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `estimate_residual_method` | char | `"MoM"` | `"MoM"`, `"MLE"`, or `"Servin_Stephens"` |
| `unmappable_effects` | char | `"none"` | `"none"`, `"inf"`, or `"ash"` (ash only for individual data) |
| `convergence_method` | char | `"elbo"` | `"elbo"` or `"pip"` (pip forced for unmappable) |
| `alpha0` | numeric | 0.1 | NIG prior shape hyperparameter |
| `beta0` | numeric | 0.1 | NIG prior rate hyperparameter |
| `stochastic_ld_sample` | int/NULL | NULL | Sketch size B for stochastic LD correction (RSS only) |

---

## 5. Architecture Overview

### Object Taxonomy

```
Data (immutable, S3-dispatched)       Params (immutable)       Model (mutable)
├─ individual: X, y, n, p             L, max_iter, tol,        alpha [L × p]
├─ ss: XtX, Xty, yty, n, p           convergence_method,      mu    [L × p]
└─ rss_lambda: z, R, λ, n, p         estimate_prior_method,   mu2   [L × p]
   + eigen decomp fields              unmappable_effects,      V     [L]
   + scaling attributes                refine, alpha0, beta0,  sigma2
   + shat2_inflation (stoch LD)        ...                     theta [p] (unmappable)
                                                                tau2  (unmappable)
                                                                Xr / XtXr / Rz
```

### Pipeline Flow

```
User Interface    Constructor              Workhorse        IBSS Core          Backend (S3)
─────────────    ──────────────           ──────────        ─────────          ────────────
susie()      →   individual_data_         susie_          ibss_initialize()   compute_residuals.{class}
                 constructor()            workhorse()  →  ibss_fit()       →  compute_ser_statistics.{class}
susie_ss()   →   sufficient_stats_                     ↓  ibss_finalize()     calculate_posterior_moments.{class}
                 constructor()                         │                      update_fitted_values.{class}
susie_rss()  →   summary_stats_                        │  SER pipeline:       update_variance_components.{class}
                 constructor()                         │  single_effect_      loglik.{class}
                 ├─ λ=0 → ss route                    │  regression()        neg_loglik.{class}
                 └─ λ>0 → rss_lambda_                 │  optimize_prior_     get_cs.{class}
                          constructor()                │  variance()          ...
                                                       ↓
                                           run_refine() (optional)
```

### S3 Generic Methods (~20+ generics in `generic_methods.R`)

Core generics dispatched by data class:

| Generic | Purpose | Classes |
|---------|---------|---------|
| `configure_data` | Pre-processing (eigen decomp for inf/rss) | all 3 |
| `get_var_y` | Variance of y | all 3 |
| `initialize_susie_model` | Create initial model state | all 3 |
| `initialize_fitted` | Initial fitted values | all 3 |
| `validate_prior` | Check prior variance | all 3 |
| `track_ibss_fit` | Track iteration history | all 3 |
| `compute_residuals` | Partial residuals for SER | all 3 |
| `compute_ser_statistics` | betahat, shat2 for SER | all 3 |
| `SER_posterior_e_loglik` | Expected log-likelihood | all 3 |
| `calculate_posterior_moments` | Posterior mean/variance | all 3 |
| `compute_kl` | KL divergence | all 3 |
| `get_ER2` | Expected residual sum of squares | all 3 |
| `Eloglik` | Expected log-likelihood for ELBO | all 3 |
| `loglik` / `neg_loglik` | Log-likelihood / neg for optim | all 3 |
| `update_fitted_values` | Update Xr/XtXr/Rz after SER | all 3 |
| `update_variance_components` | σ² (and τ² for unmappable) | all 3 |
| `update_derived_quantities` | Post-variance-update bookkeeping | all 3 |
| `get_scale_factors` / `get_intercept` | Rescaling for output | all 3 |
| `get_fitted` / `get_cs` | Fitted values / credible sets | all 3 |
| `get_variable_names` / `get_zscore` | Metadata extraction | all 3 |
| `cleanup_model` | Remove internal fields before return | all 3 |

---

## 6. Feature × Data-Type Support Matrix

| Feature | `individual` (X, y) | `ss` (XtX, Xty) | `rss_lambda` (z, R, λ>0) |
|---------|:---:|:---:|:---:|
| Standard SuSiE | ✅ | ✅ | ✅ |
| SuSiE-inf (`unmappable_effects = "inf"`) | ✅ | ✅ | ❌ |
| SuSiE-ash (`unmappable_effects = "ash"`) | ✅ | ❌ (error) | ❌ |
| Servin-Stephens NIG prior | ✅ | ✅ | ❌ |
| `estimate_residual_method = "MoM"` | ✅ | ✅ | N/A |
| `estimate_residual_method = "MLE"` | ✅ | ✅ | ✅ (profile) |
| `estimate_prior_method = "optim"` | ✅ | ✅ | ✅ |
| `estimate_prior_method = "EM"` | ✅ | ✅ | ✅ |
| `estimate_prior_method = "simple"` | ✅ | ✅ | ✅ |
| `convergence_method = "pip"` | ✅ | ✅ | ✅ |
| `convergence_method = "elbo"` | ✅ | ✅ | ✅ |
| Refinement (`refine = TRUE`) | ✅ | ✅ | ✅ |
| Stochastic LD correction | N/A | via `susie_rss(λ=0)` | ❌ (requires λ=0 path) |
| `model_init` warm start | ✅ | ✅ | ✅ |
| `track_fit` | ✅ | ✅ | ✅ |
| `compute_univariate_zscore` | ✅ | ❌ | ❌ |
| `check_input` / `check_R` / `check_z` | — | ✅ | ✅ |
| MAF filtering | — | — | ✅ (rss path) |
| `intercept` estimation | ✅ | ❌ (always FALSE) | N/A |

**Critical constraint**: Refinement is **incompatible** with unmappable effects (`refine = TRUE` + `unmappable_effects != "none"` will error or be ignored).

**Convergence override**: When `unmappable_effects ∈ {"inf", "ash"}`, convergence is forced to PIP-based (not ELBO).

---

## 7. SuSiE-ash — Adaptive Shrinkage for Unmappable Effects

### Concept
SuSiE-ash extends SuSiE with a Mr.ASH (Kim 2020, Stephens 2017) component that absorbs effects which cannot be mapped to individual sparse signals. It provides a "middle ground" between standard SuSiE (which may register false positives from background effects) and SuSiE-inf (which can be over-conservative).

### Mathematical Model
- Sparse effects: Standard SuSiE L single-effect components
- Non-sparse effects: $\theta \sim \text{Mr.ASH}(\pi, \sigma^2 \cdot \sigma_k^2)$ — mixture of normals with adaptive shrinkage
- $\tau^2 = \sum_k \sigma_k^2 \pi_k \cdot \sigma^2$ (implied unmappable variance)
- $\theta$ is fit by alternating between SuSiE updates and Mr.ASH coordinate ascent

### Key Implementation Details

**Data requirement**: Individual-level data only. `sufficient_stats_constructor()` rejects `unmappable_effects = "ash"` with an explicit error: *"Adaptive shrinkage (ash) requires individual-level data."*

**Masking logic** (in `update_variance_components.ss`, lines ~430–760 of `sufficient_stats_methods.R`):
- **Purpose**: Prevent double-counting between SuSiE's sparse effects and Mr.ASH's dense effects.
- SuSiE effects are classified into 3 cases:
  1. **Case 1 (confident)**: High purity, no sentinel collision → protected from Mr.ASH
  2. **Case 2 (borderline)**: Some diffusion or collision → partially protected
  3. **Case 3 (diffuse)**: Low purity or high collision → "donated" to Mr.ASH
- Masking: $\theta_j = 0$ at positions where SuSiE has signal (or in LD neighborhood), based on:
  - Direct PIP threshold
  - Neighborhood PIP threshold (via LD adjacency matrix)
  - Force-mask for positions under CS coverage
- **Unmasking**: Positions can be unmasked after a waiting period if Mr.ASH finds signal there, with a "second chance" mechanism.

**Mr.ASH warm starts**: Each IBSS iteration warm-starts Mr.ASH from the previous $\theta$, $\pi$, $\sigma^2$.

**Final unmasked pass** (`run_final_ash_pass()`): After SuSiE convergence, Mr.ASH is run one final time *without* masking to get accurate $\theta$ estimates for prediction.

**Diagnostic reporting** (`diagnosis_reports.R`): `diagnose_susie_ash_iter()` prints per-iteration Case assignments, masking statistics, Mr.ASH convergence info.

**Model output**: Returns `theta` (p-vector), `tau2` (scalar), along with standard SuSiE output.

### Initialization (in `initialize_susie_model.ss`)
```
model$tau2 = 0
model$theta = rep(0, p)
model$XtX_theta = rep(0, p)
model$masked = rep(FALSE, p)
model$ash_iter = 0
model$ash_pi = NULL
model$diffuse_iter_count = rep(0, L)
model$prev_sentinel = rep(0, L)
model$unmask_candidate_iters = rep(0, p)
model$ever_unmasked = rep(FALSE, p)
model$force_exposed_iter = rep(0, p)
model$ever_diffuse = rep(0, L)
model$second_chance_used = rep(FALSE, p)
model$prev_case = rep(0, L)
```

---

## 8. SuSiE-inf — Infinitesimal Effects Model

### Concept
Based on Cui et al. (2024, Nature Genetics). Decomposes genetic effects into L sparse causal signals plus an infinitesimal polygenic background: $\beta = \sum_{l=1}^L \gamma_l b_l + u$, where $u \sim N(0, \tau^2 I)$.

### Mathematical Framework
- **Covariance**: $\Omega^{-1} = \tau^2 X'X + \sigma^2 I$ (in eigen space: $\omega_k = \tau^2 d_k + \sigma^2$)
- **BLUP for infinitesimal effects**: $\theta = \tau^2 X (X'X \tau^2 + \sigma^2 I)^{-1} (y - X\hat{\beta})$
- **Omega-weighted residuals**: SER uses $X'\Omega^{-1} r$ instead of $X'r$
- **Eigen decomposition**: Pre-computed on `configure_data.ss()` via thin SVD when X available

### Data Support
| Interface | Support |
|-----------|:---:|
| `susie()` | ✅ — converts individual data to SS internally (`convert_individual_to_ss()`) |
| `susie_ss()` | ✅ — native SS operation |
| `susie_rss()` (λ=0) | ✅ — routes through SS constructor |
| `susie_rss()` (λ>0) | ❌ — not supported for `rss_lambda` class |

### Key Functions
- `compute_omega_quantities(data, tau2, sigma2)` → `{omega_var, diagXtOmegaX}` (susie_utils.R)
- `compute_theta_blup(data, model)` → p-vector of infinitesimal BLUP coefficients (susie_utils.R)
- `mom_unmappable(data, params, model, omega, tau2, ...)` → MoM for (σ², τ²) (susie_utils.R)
- `mle_unmappable(data, params, model, omega, ...)` → MLE via L-BFGS-B for (σ², τ²) (susie_utils.R)
- `compute_elbo_inf(alpha, mu, omega, ...)` → ELBO for infinitesimal model (susie_utils.R)

### Variance Estimation Methods
- **MoM** (default): Solves linear system $A \cdot (\sigma^2, \tau^2)' = x$ using eigenvalue moments
- **MLE**: L-BFGS-B optimization of profile likelihood over $(\sigma^2, \tau^2)$

### SER Modifications for SuSiE-inf
- Residuals: $X'\Omega^{-1}(y - X\hat{\beta}_{-l})$ instead of $X'(y - X\hat{\beta}_{-l})$
- `predictor_weights = diag(X'\Omega^{-1} X)` (updates each iteration as τ² changes)
- `residual_variance = 1` (already absorbed into Ω)
- Prior variance optimization on **linear** scale (bounds [0, 1]) vs log scale for standard SuSiE

### Model Output
- `theta`: p-vector of posterior BLUP means for infinitesimal effects
- `tau2`: estimated infinitesimal variance component
- Standard SuSiE output (alpha, mu, pip, sets, etc.)

---

## 9. Stochastic LD Correction

### Problem
When $R$ is approximated by a stochastic sketch $\hat{R} = \frac{1}{B} U U'$ (where $U = X_\text{std}' W$, $W \sim N(0, I_B/n)$), the null variance of $\hat{z}_j$ is inflated:

$$\text{Var}(\hat{z}_j \mid H_0) = \left(1 + \frac{1}{B}\right) \sigma^2_{j,0} + \frac{R_{jj}}{B} \|R\|_F^2 := \tau_j^2$$

### Correction Mechanism
Set `stochastic_ld_sample = B` in `susie_rss()` to activate. The correction:
1. Precomputes per-SNP inflation factors $\tau_j^2 / \sigma_{j,0}^2$ in the constructor
2. Stores as `data$shat2_inflation`
3. During SER, inflates `shat2` by this factor: `shat2 <- shat2 * data$shat2_inflation`
4. This attenuates Bayes factors for SNPs where apparent signal is due to LD noise

### Debiasing
Uses **debiased Wishart moment estimators** (Proposition 4.4) to avoid upward bias:
```
R_frob_sq_db = (B * R_frob_sq - p^2) / (B + 1)
ell_j_db     = (B * ell_j - p) / (B + 1)
sigma2_j0_db = pmax((B^2 * sigma2_j0 - (B+1)(2p*ell_j_db + R_frob_sq_db) - p^2) / (B^2 + 3B + 4), 1)
tau2_j       = (1 + 1/B) * sigma2_j0_db + d_R * R_frob_sq_db / B
```

### Computation Paths
- **From X (sketch matrix)**: When `X = t(U)` is provided, computes from Gram matrix $A = X X' = U'U$
  - $d_R = \text{colSums}(X^2) / B$ (diagonal of $\hat{R}$)
  - Uses $A^2 X'$ for higher-order moments
- **From R (LD matrix)**: When only R is provided, computes directly from R entries

### Output Diagnostics
When active, `fit$stochastic_ld_diagnostics` contains:
- `B`: sketch sample size
- `p`: number of SNPs
- `effective_rank`: $\hat{r} = p^2 / \|R\|_F^2$
- `r_over_B`: $\hat{r}/B$ (values ≤ 0.2 → sketch adequate)
- `per_snp_inflation`: $\tau_j^2 / \sigma_{j,0}^2 - 1$ (values ≤ 0.2 → minimal power loss)

### Benchmark Notebook
`inst/notebooks/stochastic_ld_benchmark.ipynb` — compares 6 methods:
1. Gold standard (in-sample R)
2. Subsample R from B individuals
3. Stochastic sketch, no correction
4. Stochastic sketch, with correction
5. Stochastic sketch + NIG prior, no correction
6. Stochastic sketch + NIG prior, with correction

Tested with $B \in \{2000, 5000, 10000\}$ and block-diagonal LD structures.

---

## 10. Normal-Inverse-Gamma (NIG) / Servin-Stephens Prior

### Concept
Based on Servin & Stephens (2007). Places an Inverse-Gamma prior on $\sigma^2$:

$$\sigma^2 \sim \text{IG}(\alpha_0/2, \beta_0/2)$$

This integrates out $\sigma^2$ analytically, producing $t$-distributed marginals instead of Gaussian. Primary motivation: improved CS coverage when **$n$ is small** (recommended when $n < 80$).

### Activation
```r
susie(X, y, estimate_residual_method = "Servin_Stephens",
      estimate_prior_method = "EM", alpha0 = 0.1, beta0 = 0.1)
```

### Key Implementation Functions (all in `susie_utils.R`)

| Function | Purpose |
|----------|---------|
| `compute_lbf_servin_stephens(x, y, s0, alpha0, beta0)` | Log BF for individual-level data |
| `compute_lbf_NIG(n, xx, xy, yy, sxy, s0, a0, b0, tau)` | Log BF for SS data |
| `compute_stats_NIG(n, xx, xy, yy, sxy, s0, a0, b0, tau)` | Full NIG statistics (lbf + posterior moments + rv) |
| `compute_posterior_moments_NIG(n, xx, xy, yy, sxy, s0, a0, b0, tau)` | Posterior mean, variance, rv |
| `update_prior_variance_NIG_EM(n, xx, xy, yy, sxy, pip, s0, a0, b0, tau)` | EM update for prior variance |
| `compute_kl_NIG(alpha, post_mean, post_mean2, pi, V, a0, b0, a_post, b_post)` | KL divergence |
| `compute_null_loglik_NIG(n, yy, a0, b0, ...)` | Null log-likelihood |
| `compute_marginal_loglik(lbf_model, n, yy, a0, b0, ...)` | Marginal log-likelihood |
| `inv_gamma_factor(a, b)` | $a \log b - \log \Gamma(a)$ helper |
| `get_nig_sufficient_stats(data, model)` | Extracts (yy, sxy, tau) regardless of data type |
| `posterior_mean_servin_stephens(xtx, xty, s0_t)` | Posterior mean via SS |
| `posterior_var_servin_stephens(xtx, xty, yty, n, s0_t)` | Posterior variance via SS |

### Mathematical Details

**Log Bayes Factor**:
$$\text{lbf}_j = -\frac{1}{2}\left[\log(1 + s_0 x_j'x_j / \tau) + (\alpha_0 + n) \log\frac{\beta_0 + \text{RSS}_j}{\beta_0 + y'y}\right]$$

where $r_0 = s_0 / (s_0 + \tau / x_j'x_j)$, $\text{RSS}_j = y'y(1 - r_0 \cdot r_{xy,j}^2)$.

**Posterior moments**:
- $\hat{\beta}_j^{\text{post}} = r_0 \cdot \hat{\beta}_j^{\text{LS}}$
- $\text{Var}(\beta_j | \text{data}) = \frac{b_1}{a_1 - 2} \cdot \frac{r_0 \tau}{x_j' x_j}$
- $\sigma^2_{\text{mode}} = b_1 / (a_1 - 2)$ where $a_1 = \alpha_0 + n$, $b_1 = \beta_0 + \text{RSS}$

**Residual variance**: Not estimated separately — the NIG prior inherently integrates out $\sigma^2$. The code stores per-effect posterior modes in `model$rv[l]` and uses `mean(model$rv)` when updating `sigma2`.

### Parameter Override Behavior
- `Servin_Stephens` **forces** `estimate_residual_variance = TRUE` (warns if FALSE)
- `Servin_Stephens` with `L > 1` **forces** `convergence_method = "pip"` (ELBO not well-defined for NIG with multiple effects)
- `estimate_prior_method` **forced to "EM"** when Servin_Stephens is active

### Support
- `susie()`: ✅ (full support)
- `susie_ss()`: ✅ (via `get_nig_sufficient_stats` which handles both data types)
- `susie_rss()`: ❌ (not exposed in `rss_lambda` backend)

### Validation Script
`inst/code/small_sim.R` — 250 replicates with n=40, using Thyroid/FMO2 real genotypes, comparing standard SuSiE vs NIG SuSiE.

---

## 11. Vignettes

| Vignette | Description |
|----------|-------------|
| `susie_unmappable_effects.Rmd` | *"Fine-mapping with SuSiE-ash and SuSiE-inf"* — n=1000, p=5000, 3 sparse + 5 oligogenic + 15 polygenic effects. Demonstrates standard SuSiE → SuSiE-inf → SuSiE-ash comparison. |
| `finemapping_summary_statistics.Rmd` | RSS fine-mapping with z-scores and LD matrices. Shows equivalence with individual-level data when in-sample LD is used. |
| `small_sample.Rmd` | Servin-Stephens NIG prior for small-sample fine-mapping. |
| `announcements.Rmd` | Release notes for v2.0.0. |
| (plus existing v1 vignettes) | |

---

## 12. Relevance to SuSiNE

### What SuSiNE can learn from / build on:

| susieR 2.0 Feature | SuSiNE Relevance |
|---------------------|-------------------|
| S3 dispatch architecture | SuSiNE currently has a monolithic structure. Could adopt the Data/Params/Model pattern. |
| SuSiE-inf | Direct overlap — SuSiNE's nonzero prior mean ($\mu_0 = c \cdot a$) is orthogonal to infinitesimal effects. Could combine. |
| SuSiE-ash | Different philosophy from SuSiNE's directional priors, but the masking logic for preventing double-counting is instructive. |
| NIG prior | Could be integrated into SuSiNE for small-sample robustness. |
| Stochastic LD | Directly useful for SuSiNE's real-data pipeline (eQTL studies with reference-panel LD). |
| `convergence_method = "pip"` | SuSiNE already uses ELBO-softmax BMA; PIP convergence is a simpler alternative worth considering. |
| Refinement | SuSiNE already has restart-based exploration; `run_refine()` is a lighter alternative. |
| Constructor pattern | Clean separation of concerns that SuSiNE could adopt for the `susine()` / `susine_ss()` / `susine_rss()` interfaces. |
| Reference test suite | SuSiNE should have reference tests against known-good outputs for regression safety. |

### susieR 2.0 features NOT in SuSiNE's scope:
- Mr.ASH coordinate ascent (SuSiNE uses nonzero-mean Gaussian, not mixture-of-normals)
- `rss_lambda` class (SuSiNE only needs λ=0 path where RSS reduces to SS)
