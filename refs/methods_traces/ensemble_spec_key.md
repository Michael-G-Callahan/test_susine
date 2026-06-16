# Ensemble Spec Key

## Letter prefixes

| Prefix | Model family | Description |
|--------|-------------|-------------|
| **A** | SuSiE (vanilla) | Use case `susine_vanilla`: zero prior mean, fixed σ₀² = 0.2, no EB. Runs through the local `susine` backend (not `susieR`), so "SuSiE" here means a SuSiE-equivalent zero-prior-mean fit. |
| **B** | SuSiNE (EB) | Use case `susine_eb_clamped_scale_var_nonneg`: EB scale + variance, clamped & non-negative annotation scale |
| **C** | SuSiNE (grid) | Use case `susine_functional_mu`: directional prior μ₀ = c·a; c is the exploration knob, fixed σ₀² |

## Exploration axis codes

| Code | Axis | Values |
|------|------|--------|
| **R** | Random restarts | Dirichlet-initialized α (concentration = 0.001) |
| **F** | Refine (BFS) | CS-blocking refinement; each run-table row = 2 model fits |
| **S** | σ₀² grid | Log-spaced from 0.01–1.0, 0.2 always included |
| **C** (in C-*) | c grid | Linear from 0 to 1.5 |

## Full spec catalog

| Spec | Model | Exploration axes | Run-table rows | Effective fits | Notes |
|------|-------|-----------------|---------------|----------------|-------|
| **A-R** | SuSiE vanilla | 64 restarts | 64 | 64 | Pure restart lever |
| **A-F** | SuSiE vanilla | 32 refine steps | 32 | 64 | Each step = 2 fits |
| **A-S** | SuSiE vanilla | 64 σ₀² values | 64 | 64 | Log-spaced grid |
| **A-RF** | SuSiE vanilla | 8 restarts × 4 refine | 64 | ~64 | 4 refine/parent × 2 fits |
| **A-RS** | SuSiE vanilla | 8 restarts × 8 σ₀² | 64 | 64 | Full Cartesian |
| **A-FS** | SuSiE vanilla | 8 σ₀² × 4 refine | 64 | ~64 | 4 refine/σ value × 2 fits |
| **A-RFS** | SuSiE vanilla | 4 restarts × 4 σ₀² × 2 refine | 64 | 64 | All three levers |
| **B-R** | SuSiE EB | 64 restarts | 64 | 64 | Like A-R but with EB σ₀² |
| **B-F** | SuSiE EB | 32 refine steps | 32 | 64 | Like A-F but with EB σ₀² |
| **B-RF** | SuSiE EB | 8 restarts × 4 refine | 64 | ~64 | Like A-RF but with EB σ₀² |
| **C-C** | SuSiNE | 64 c values | 64 | 64 | Pure c-grid, σ₀² = 0.2 |
| **C-CS** | SuSiNE | 8 c × 8 σ₀² | 64 | 64 | Joint c × σ₀² grid |
| **C-CSR** | SuSiNE | 8 c × 8 σ₀² search, then exact default-prior refit | 64 | 64 | `cs_grid_refit`: warm refit using exact alpha, b_hat, b_2_hat |
| **C-CF** | SuSiNE | 8 c × 4 refine | 64 | ~64 | c-grid + refinement, σ₀² = 0.2 fixed |
| **C-CFR** | SuSiNE | 4 c × 4 restarts × 2 refine | 64 | 64 | c + restart + refine |
| **C-CSFR** | SuSiNE | 3 c × 3 σ₀² × 3 refine × 3 restarts | 81 | 81 | All four levers |

## Baselines and reference bounds

These appear in figures as horizontal reference lines or separate panels alongside the main specs.

### Lower bound

**baseline-single** — One cold-start fit per model family (A, B, C), no exploration.
The floor: represents vanilla SuSiE/SuSiNE out of the box.

### Oracle upper bounds (truth-aware, not achievable in practice)

**truth-warm** — Single fit warm-started from the true top-10 causal variants (one per
effect). Shows the ceiling for a perfect initializer. Run for all three model families
(vanilla, EB, SuSiNE with c=0.5).

**oracle** (aggregation) — Post-hoc selection of the single best individual run per
dataset by AUPRC, using ground truth. Upper bound on what any aggregation rule could
achieve. Computed from the same confusion bins as the main specs; not a separate SLURM
run.

### Interpretation

```text
baseline-single  ≤  ensemble methods  ≤  oracle  ≤  truth-warm
(no exploration)       (practical)    (best in pool) (perfect init)
```

The gap between **baseline-single** and **oracle** is the total exploitable headroom in the pool.
The gap between **oracle** and any aggregation method (e.g., elbo_softmax) is the cost of not knowing which fit to pick.
The gap between **oracle** and **truth-warm** reflects whether the best discovered basin is as good as the truth-initialized one.

## Aggregation methods

Applied post-hoc to each spec's pool of fits.

| Method (`agg_method` ID) | Figure label | Description |
|--------|--------------|-------------|
| **max_elbo** | max ELBO | Select single best fit by ELBO |
| **uniform** | uniform | Average PIPs equally across all fits |
| **elbo_softmax** | ELBO softmax | ELBO-softmax weighted average of PIPs |
| **cluster_weight_credible** | cluster softmax | Cluster fits by `credible_shift` = max_j[ max(p_j,q_j)·\|p_j−q_j\| ], complete-linkage cut at 0.05; aggregate with Method B (frequency-free). Primary scheduled cluster-weight method. |
| **cluster_weight_max** | cluster (max) | Cluster fits by `max_shift` = max_j \|p_j−q_j\| (L-infinity), complete-linkage cut at 0.10; aggregate with Method B (frequency-free). Dropped from the delta-AUPRC heatmap panel. |
| **oracle** | oracle | Best individual fit by AUPRC (truth-aware; reference only) |

The scheduled method list comes from the run-control workbook
(`aggregation_methods`) and the `.cluster_weight_specs` registry in
`run_model.R`, NOT from `aggregation_catalog()` in `use_cases.R` (which is stale
and still lists only the legacy `cluster_weight` / `cluster_weight_jsd_050`
JSD methods). Figure display labels are from `agg_labels` in the paper-prep
workbook.

**Method B (frequency-free):** each cluster gets weight proportional to
`exp(max-ELBO-in-cluster)`; that weight is split within the cluster across member
fits by ELBO-softmax. All fits contribute (no single nominee), so no 1/frequency
correction is applied. The legacy JSD-0.15 + frequency `cluster_weight` method
(Method C, max-ELBO nominee + inverse-frequency) still exists in code for
backward compatibility, but it is not used by the current ensemble-simulation or
real-data paper workflows.

## Manuscript figure mapping

The three manuscript-facing PNGs (in `../Writings/plots/ensemble_sims/`) are
the only files the rendering workbook writes; `write_component_pngs = FALSE`
suppresses individual panels. They are produced by
`visualize_results_workbook_ensemble_scaling_paper.Rmd` from the prepared RDS.

| Manuscript PNG | Composite panels | Content |
|----------------|------------------|---------|
| `paper_ensemble_scaling_composite_heatmap.png` | A | Delta-AUPRC heatmap (spec × aggregation; `cluster_weight_max` dropped) |
| `paper_ensemble_scaling_composite_performance.png` | B, C, D, E, F | B: C-CS AUPRC vs annotation; C: PIP calibration; D: PR curves; E: heritability; F: AUPRC vs compute |
| `paper_ensemble_multimodality.png` | G | Faceted multimodality boxplots (metric × model family), per-locus |

Panel E (heritability) inputs are spliced from the separate hotfix job
`ensemble_scaling_hg2_hotfix/53588814`; all other panels read the main job
`ensemble_scaling_full/53547760`. Panel G prefers the `n_clusters_credible` /
`max_credible_dist` multimodal columns, falling back to the JSD-era columns when
absent.
