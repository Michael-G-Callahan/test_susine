# Simulation Design Note: SuSiE 2 S4 vs test_susine

**Date:** 2026-03-02  
**Scope:** Compare simulation design in SuSiE 2.0 Supplementary Note S4 with current `test_susine` implementation, and document current project decision.

## 1. Source documents

- SuSiE 2.0 supplemental text: `refs/citations/SuSiE_2.txt` (Section `S.4 Simulation Study Details`)
- Current simulator implementation: `R/simulate_data.R`
- Current matrix catalog examples: `data/sampled_simulated_genotypes/scenario_sampling_summary.csv`

## 2. SuSiE 2.0 S4 simulation design (summary)

### 2.1 Sparse setting

S4 describes a sparse simulation design using UK Biobank-derived LD blocks with:

- nominal dimensions: `n = 1000`, `p = 5000`
- causal count grid: `k = {1,2,3,4,5}`
- per-SNP heritability fixed at `h2_snp = 0.03`
- effect construction: causal effects set uniformly (then noise variance set via heritability equations)
- 150 replicates per scenario

Interpretation: with fixed per-SNP heritability and fixed `k`, causal effect magnitudes are close to uniform (modulo sign/scale conventions), rather than heterogeneous draws.

### 2.2 Complex (oligogenic/polygenic) settings

S4 defines explicit effect partitions and target heritability shares, e.g.:

- Scenario A: `3 sparse (50%) + 5 oligogenic (35%) + 15 polygenic (15%)`
- Scenario B: same sparse/oligogenic core with remaining heritability spread across all remaining variants
- Scenario C: `3 sparse (50%) + 10 oligogenic (15%) + remaining variants (35%)`

Other details:

- total heritability fixed (`h2_total = 0.25`)
- oligogenic component uses a two-component Gaussian mixture for effect magnitudes
- 150 replicates per scenario

## 3. Current test_susine simulation design (summary)

### 3.1 Sparse setting (`simulate_effect_sizes`)

Current sparse generator in `R/simulate_data.R`:

- samples causal indices
- samples causal effects from a Normal distribution (`rnorm`) with configurable `effect_sd`
- leaves non-causal effects at 0

This naturally produces heterogeneous causal effect magnitudes.

### 3.2 Oligogenic setting (`simulate_effect_sizes_oligogenic`)

Current oligogenic generator:

- builds sparse/oligogenic/polygenic tiers as functions of `p_star`
- uses tier-specific random draws (including a Bernoulli-driven oligogenic mixture)
- rescales tiers to target energy fractions (default `0.6/0.3/0.1`)

This is not parameterized to exactly reproduce S4's fixed counts/heritability partitions.

### 3.3 Phenotype/noise parameterization

Current phenotype simulation is noise-fraction based (`y_noise`), not explicitly heritability-targeted by `h2_snp`/`h2_total`.

### 3.4 Matrix dimensions/sources

Current runs use:

- built-in `SuSiE_N3_X` (observed ~`574 x 1001`), or
- sampled scenario matrices with `participant_count` often 600 or 6000 and `snps_post` typically around ~700-1500 (varies by locus/scenario)

So default matrix dimensions differ from S4's `1000 x 5000` design.

## 4. Key differences

1. **Effect-size heterogeneity (sparse):**  
   S4 is close to equal-magnitude causal effects under fixed per-SNP heritability; current repo uses random Normal causal effects.

2. **Primary control knob:**  
   S4 is heritability-driven (`h2_snp`, `h2_total`); current repo is noise-fraction-driven (`y_noise`).

3. **Complex architecture parameterization:**  
   S4 uses fixed count/share templates (e.g., 3/5/15 with 50/35/15); current repo uses `p_star`-driven tier sizes and default energy fractions.

4. **Infinitesimal background scenario:**  
   S4 explicitly includes "all remaining variants nonzero" backgrounds; current repo does not make this the default architecture mode.

5. **Dimensional regime and replicates:**  
   S4 emphasizes `n=1000, p=5000`, 150 replicates; current repo varies dimensions by matrix catalog and uses seed count from run controls.

## 5. Commentary from project discussion

Current project perspective:

- Forcing equal-ish causal magnitudes in sparse settings can feel unnatural and may underrepresent realistic locus-to-locus effect heterogeneity.
- The S4 complex architecture setup is useful for specific method stress tests, but is comparatively prescriptive and heavier to maintain.
- For the present study goals (exploration/aggregation workflow benchmarking), a simpler sparse generator with random Normal causal effects is preferred.

## 6. Decision (as of 2026-03-02)

**Decision:** keep current `test_susine` sparse simulation method as the primary study simulator.

Rationale:

1. preserves natural variation in causal effect magnitudes,
2. avoids over-constraining architecture assumptions in primary benchmarks,
3. keeps run-control complexity manageable.

## 7. Optional future work (not current default)

If needed for external comparability:

1. Add an S4-style compatibility mode as a sensitivity scenario (not default), with:
   - explicit `h2_total` targeting,
   - fixed sparse/oligo/polygenic count templates,
   - optional infinitesimal-background mode.
2. Keep current sparse Normal-draw design as the mainline benchmark path.
