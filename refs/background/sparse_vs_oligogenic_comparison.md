# Sparse vs. Oligogenic Simulation: test_susine vs. SuSiE 2.0

**Date:** 2026-03-05
**Extends:** `simulation_design_susie2_vs_test_susine.md` (2026-03-02, sparse-focused)
**Scope:** Deep comparison of the oligogenic/complex architecture simulation, mapping p_star regimes to SuSiE 2 scenarios, and proposing concrete changes for alignment.

---

## 1. SuSiE 2 Complex Architecture (S.4, lines 872-928)

SuSiE 2 defines three scenarios under a **fixed total heritability** h2_total = 0.25, with 150 replicates each, using UK Biobank LD blocks (n=1000, p=5000).

### 1.1 Three-component partition

All scenarios partition variants into sparse (S), oligogenic (O), and polygenic (P) tiers:

- **Sparse (S):** k=3 variants, effects drawn from N(0, sigma_S^2), scaled to h2_S
- **Oligogenic (O):** 5-10 variants, effects from a **two-component Gaussian mixture**:
  - z_j ~ Categorical(pi_1, pi_2)
  - beta_j | z_j ~ N(0, sigma_{z_j}^2), where sigma_1^2 < sigma_2^2
  - Scaled to h2_O
- **Polygenic (P):** 15 variants (Scenario 1) or all remaining (Scenarios 2-3), effects from N(0, sigma_P^2), scaled to h2_P

### 1.2 Scenario table

| Scenario | Sparse | Oligogenic | Polygenic | h2 split (S/O/P) | Background |
|----------|--------|------------|-----------|-------------------|------------|
| 1 (Fig 1F-I) | 3 | 5 | 15 | 50/35/15 | Finite polygenic |
| 2 (Fig S2 A-D) | 3 | 5 | all remaining | 50/35/15 | Moderate infinitesimal |
| 3 (Fig S2 E-H) | 3 | 10 | all remaining | 50/15/35 | Extensive infinitesimal |

### 1.3 Key design choices

1. **Total causal count is fixed**, not varied as a grid dimension
2. **h2 targeting is explicit**: effect vectors scaled so var(X beta_tier) matches target h2_tier, then residual variance sigma^2 set by Eq. S28
3. **Oligogenic mixture is Gaussian**: two variance components, not a Bernoulli switch on SD
4. **Sparse effects also Normal**: N(0, sigma_S^2), *not* uniform beta=1 (that's the sparse-only setting). The complex setting uses Normal draws for all tiers.
5. **Infinitesimal = all remaining variants**, not a finite count

---

## 2. Current test_susine Implementation

**Source:** `R/simulate_data.R`, `simulate_effect_sizes_oligogenic()` (lines 68-137)

### 2.1 Three-tier architecture

| Tier | Count | Effect distribution | Default h2 fraction |
|------|-------|-------------------|-------------------|
| Sparse | p_star | N(0, effect_sd^2) | 0.60 |
| Oligogenic | ceil(p_star/2) | Bernoulli mixture: 35% at 0.8x SD, 65% at 0.25x SD | 0.30 |
| Polygenic | ceil(1.5 * p_star) | N(0, (0.08 * effect_sd)^2) | 0.10 |

### 2.2 Energy rescaling

After drawing effects, each tier is L2-energy-rescaled via `scale_tier_to_energy()`:
```
beta[tier_idx] *= sqrt(target_energy / current_energy)
```
where `target_energy = tier_frac * total_energy` and `total_energy = sum(beta^2)`.

This preserves within-tier heterogeneity (relative magnitudes) while hitting the target energy fraction.

### 2.3 Phenotype generation

`simulate_phenotype()` uses noise-fraction control:
```
sigma^2 = var(X beta) * (noise_fraction / (1 - noise_fraction))
```
So `h2 = 1 - noise_fraction` (approximately, when X is standardized). The mapping:

| y_noise | Implied h2 |
|---------|-----------|
| 0.50 | 0.50 |
| 0.80 | 0.20 |
| 0.95 | 0.05 |

---

## 3. Feature-by-Feature Comparison

### 3.1 Signal control knob

| | test_susine | SuSiE 2 |
|-|------------|---------|
| Primary knob | `y_noise` (noise fraction) | `h2_total` (total heritability) |
| Relationship | h2 = 1 - y_noise | h2_total = 0.25 (fixed) |
| Tier targeting | Energy fractions (L2 norm) | Heritability fractions (genetic variance) |

**Substantiveness:** Mostly a reparameterization. Energy fractions approximate heritability fractions when X is standardized (columns have unit variance), because `var(X beta_tier) approx sum(beta_tier^2)` when X has orthonormal columns. With real LD structure, this approximation degrades — `var(X beta_tier)` depends on LD among causal SNPs, while `sum(beta_tier^2)` does not. SuSiE 2's h2-based targeting is more precise for controlled simulation.

### 3.2 Tier sizing

| | test_susine | SuSiE 2 |
|-|------------|---------|
| Sparse | p_star (grid parameter) | Fixed: 3 |
| Oligogenic | ceil(p_star/2) | Fixed: 5 (Sc. 1-2) or 10 (Sc. 3) |
| Polygenic | ceil(1.5 * p_star) | Fixed: 15 (Sc. 1) or all remaining (Sc. 2-3) |
| Total causal | 2.5 * p_star (approx) | 23 (Sc. 1) or 8-13 + all remaining |

**Substantiveness:** This is a significant structural difference. In test_susine, increasing p_star scales all tiers proportionally — the "shape" of the architecture is invariant, only the total causal count changes. In SuSiE 2, the architecture is fixed per scenario with specific ecological interpretations (e.g., "3 strong eQTL + 5 moderate + background").

### 3.3 Oligogenic mixture model

| | test_susine | SuSiE 2 |
|-|------------|---------|
| Model | Bernoulli switch on SD | Two-component Gaussian mixture |
| Components | 35% at 0.8x base SD, 65% at 0.25x base SD | z_j ~ Cat(pi_1, pi_2), then N(0, sigma_{z_j}^2) |
| Effect | Creates bimodal |effect| within tier | Creates bimodal |effect| within tier |
| Parameterization | Hardcoded proportions/multipliers | Unspecified pi_1, pi_2, sigma_1^2, sigma_2^2 |

**Substantiveness:** Both produce bimodal effect magnitudes within the oligogenic tier. The current Bernoulli-on-SD approach draws from N(0, 0.8x) or N(0, 0.25x) — i.e., the variance ratio between components is (0.8/0.25)^2 = 10.24. The SuSiE 2 paper doesn't specify exact (pi_1, pi_2, sigma_1^2, sigma_2^2) values, only that sigma_1^2 < sigma_2^2. The qualitative behavior is similar; the exact parameterization matters less than the h2-rescaling that follows.

### 3.4 Polygenic / infinitesimal background

| | test_susine | SuSiE 2 |
|-|------------|---------|
| Finite poly | ceil(1.5 * p_star) SNPs | 15 SNPs (Sc. 1 only) |
| Infinitesimal | Not implemented | All remaining SNPs (Sc. 2-3) |
| Effect size | 0.08x base SD | N(0, sigma_P^2), scaled to h2_P |

**Substantiveness:** This is the most important structural difference. The infinitesimal background (all remaining ~4970 SNPs contributing small effects) creates a very different fitting challenge than a finite polygenic tier of ~8-30 SNPs. The infinitesimal background is specifically what SuSiE-inf was designed to handle. Without an infinitesimal mode, test_susine cannot directly replicate the benchmarks showing SuSiE-inf's advantage.

### 3.5 h2 fraction defaults

| | test_susine | SuSiE 2 Sc. 1 | SuSiE 2 Sc. 3 |
|-|------------|---------------|---------------|
| Sparse | 60% | 50% | 50% |
| Oligogenic | 30% | 35% | 15% |
| Polygenic | 10% | 15% | 35% |

**Substantiveness:** Modest difference. test_susine allocates more to sparse (60% vs 50%) and less to polygenic (10% vs 15%). Easily aligned by changing `tier_h2` defaults.

---

## 4. What p_star Range Maps to SuSiE 2's "Realistic/Oligogenic"?

### 4.1 Total causal count mapping

With test_susine's proportional tier sizing (p_star + ceil(p_star/2) + ceil(1.5*p_star)):

| p_star | Sparse | Oligo | Poly | Total | Closest SuSiE 2 scenario |
|--------|--------|-------|------|-------|--------------------------|
| 1 | 1 | 1 | 2 | 4 | None (too sparse) |
| 3 | 3 | 2 | 5 | 10 | Between sparse (k=3) and Sc. 1 |
| 5 | 5 | 3 | 8 | 16 | Approaching Sc. 1 (23 total) |
| 10 | 10 | 5 | 15 | 30 | **Close to Sc. 1** (3+5+15=23) |
| 15 | 15 | 8 | 23 | 46 | Between Sc. 1 and Sc. 3 |
| 20 | 20 | 10 | 30 | 60 | Beyond any SuSiE 2 finite scenario |

### 4.2 Interpretation

**p_star in {8, 10, 12} is the closest match to SuSiE 2 Scenario 1** (23 total finite causals). At p_star=10 (total=30), the tier breakdown {10, 5, 15} is remarkably close to SuSiE 2's {3, 5, 15} in the oligogenic and polygenic tiers, though with more sparse causals.

However, the **qualitative character differs**:
- SuSiE 2 Sc. 1 has 3 strong + 5 moderate + 15 weak = specific "eQTL-like" architecture
- test_susine at p_star=10 has 10 strong + 5 moderate + 15 weak = more "sparse-heavy" architecture

To truly match SuSiE 2 Sc. 1, test_susine would need `tier_mode="fixed"` with `tier_counts=c(3, 5, 15)`.

### 4.3 What about the sparse setting?

SuSiE 2's sparse setting uses k = {1,2,3,4,5} with uniform beta=1 and h2_snp=0.03. This maps directly to test_susine's `architecture="sparse"` with `p_star={1,2,3,4,5}` — but the effect magnitude distribution differs (Normal vs uniform). The p_star range in the current pilot ({1,3,5,10,20}) covers and exceeds SuSiE 2's sparse range.

---

## 5. Proposed Changes to `simulate_effect_sizes_oligogenic()`

All changes are backward-compatible (defaults reproduce current behavior).

### 5.1 New arguments

```r
simulate_effect_sizes_oligogenic <- function(
  p,
  p_star,
  effect_sd = 1,
  seed = NULL,
  tier_h2 = c(sparse = 0.6, oligogenic = 0.3, polygenic = 0.1),
  # --- NEW ARGUMENTS ---
  tier_mode = c("proportional", "fixed"),
  tier_counts = NULL,              # e.g., c(sparse = 3, oligogenic = 5, polygenic = 15)
  oligo_mixture = c("bernoulli", "gaussian_2comp"),
  oligo_2comp_probs = c(0.5, 0.5), # pi_1, pi_2 for gaussian_2comp
  oligo_2comp_sd_ratio = 3.0,      # sigma_2 / sigma_1 for gaussian_2comp
  sparse_effect_mode = c("normal", "uniform"),
  infinitesimal = FALSE            # if TRUE, polygenic = all remaining SNPs
)
```

### 5.2 Tier sizing logic

```r
tier_mode <- match.arg(tier_mode)
if (tier_mode == "fixed") {
  stopifnot(!is.null(tier_counts), length(tier_counts) == 3)
  n_sparse <- tier_counts[["sparse"]]
  n_oligo  <- tier_counts[["oligogenic"]]
  n_poly   <- if (infinitesimal) p - n_sparse - n_oligo else tier_counts[["polygenic"]]
} else {
  # Current proportional logic
  n_sparse <- p_star
  n_oligo  <- max(1L, as.integer(ceiling(p_star / 2)))
  n_poly   <- if (infinitesimal) p - n_sparse - n_oligo else max(1L, as.integer(ceiling(1.5 * p_star)))
}
```

### 5.3 Oligogenic mixture

```r
oligo_mixture <- match.arg(oligo_mixture)
if (oligo_mixture == "bernoulli") {
  # Current: 35% at 0.8x SD, 65% at 0.25x SD
  mix <- stats::rbinom(n_oligo, size = 1, prob = 0.35)
  sd_vec <- ifelse(mix == 1L, effect_sd * 0.8, effect_sd * 0.25)
  beta[oligo_idx] <- stats::rnorm(n_oligo, mean = 0, sd = sd_vec)
} else if (oligo_mixture == "gaussian_2comp") {
  # SuSiE 2 style: two-component Gaussian
  z <- sample(1:2, size = n_oligo, replace = TRUE, prob = oligo_2comp_probs)
  # Set sigma_1 and sigma_2 so that sigma_2/sigma_1 = oligo_2comp_sd_ratio
  # and the mean variance matches effect_sd^2
  sigma_1 <- effect_sd / sqrt(oligo_2comp_probs[1] + oligo_2comp_probs[2] * oligo_2comp_sd_ratio^2)
  sigma_2 <- sigma_1 * oligo_2comp_sd_ratio
  sd_vec <- ifelse(z == 1L, sigma_1, sigma_2)
  beta[oligo_idx] <- stats::rnorm(n_oligo, mean = 0, sd = sd_vec)
}
```

### 5.4 Sparse effect mode

```r
sparse_effect_mode <- match.arg(sparse_effect_mode)
if (sparse_effect_mode == "normal") {
  beta[sparse_idx] <- stats::rnorm(n_sparse, mean = 0, sd = effect_sd)
} else if (sparse_effect_mode == "uniform") {
  # SuSiE 2 sparse-only convention: all effects = effect_sd, random signs
  signs <- sample(c(-1, 1), size = n_sparse, replace = TRUE)
  beta[sparse_idx] <- signs * effect_sd
}
```

### 5.5 Infinitesimal background

When `infinitesimal = TRUE`, the polygenic tier includes all SNPs not in sparse or oligogenic. With `tier_mode = "fixed"`, this is `p - n_sparse - n_oligo` SNPs (potentially ~4990 for SuSiE 2 dimensions). The effect_sd for the polygenic tier would be very small (calibrated by the h2 rescaling step).

### 5.6 Run control example: SuSiE 2 Scenario 1

```r
job_inputs <- list(
  architecture_grid = "oligogenic",
  oligo_config = list(
    tier_mode = "fixed",
    tier_counts = c(sparse = 3, oligogenic = 5, polygenic = 15),
    tier_h2 = c(sparse = 0.50, oligogenic = 0.35, polygenic = 0.15),
    oligo_mixture = "gaussian_2comp",
    sparse_effect_mode = "normal"  # SuSiE 2 uses Normal for complex setting
  ),
  y_noise = 0.75  # h2_total = 0.25
)
```

### 5.7 Run control example: SuSiE 2 Scenario 2 (moderate infinitesimal)

```r
oligo_config = list(
  tier_mode = "fixed",
  tier_counts = c(sparse = 3, oligogenic = 5, polygenic = NA),  # NA = all remaining
  tier_h2 = c(sparse = 0.50, oligogenic = 0.35, polygenic = 0.15),
  oligo_mixture = "gaussian_2comp",
  infinitesimal = TRUE
)
```

---

## 6. Summary of Substantive vs. Cosmetic Differences

| Difference | Substantive? | Impact on comparability | Proposed fix |
|-----------|-------------|------------------------|-------------|
| y_noise vs h2 targeting | Low | Reparameterization; near-equivalent for standardized X | Document mapping |
| Tier sizing: proportional vs fixed | **High** | Changes architecture shape with p_star | Add `tier_mode="fixed"` |
| Oligo mixture: Bernoulli vs 2-comp Gaussian | Medium | Both produce bimodal, but different tail behavior | Add `oligo_mixture="gaussian_2comp"` |
| No infinitesimal background | **High** | Cannot benchmark SuSiE-inf advantage | Add `infinitesimal=TRUE` |
| h2 fractions: 60/30/10 vs 50/35/15 | Low | Easily changed via `tier_h2` arg | Already configurable |
| Sparse effects: Normal vs uniform | Low | Both reasonable; Normal more realistic | Add `sparse_effect_mode="uniform"` opt-in |
| Matrix dimensions (variable vs 1000x5000) | Medium | Affects regime interpretation | Use appropriate matrix catalog |

---

## 7. Recommendations

1. **For the current pilot study (Phase A-B):** Keep the proportional-sizing, Normal-draw design as primary. It tests the exploration/aggregation machinery across a range of architectures, which is the paper's focus.

2. **For external comparability (optional sensitivity):** Implement the SuSiE-2-compatible mode as a separate architecture option, usable via `oligo_config` in `make_job_config()`. Run a small sensitivity study (not full factorial) to show results are qualitatively similar.

3. **For SuSiE-inf benchmarking:** The infinitesimal background is essential. Without it, the advantage of `unmappable_effects="inf"` cannot be demonstrated. Prioritize implementing `infinitesimal=TRUE`.

4. **p_star interpretation for the paper:** When discussing results, note that p_star=10-15 in the proportional-sizing design corresponds roughly to SuSiE 2 Scenario 1 in total causal count (~23-46), but with a different sparse/oligo/poly balance. The paper should state this mapping explicitly.
