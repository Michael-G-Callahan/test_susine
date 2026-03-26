# Ensemble Spec Key

## Letter prefixes

| Prefix | Model family | Description |
|--------|-------------|-------------|
| **A** | SuSiE (vanilla) | Fixed σ₀² = 0.2, no EB prior variance estimation |
| **B** | SuSiE (EB) | EB prior variance estimation, clamped & non-negative |
| **C** | SuSiNE | Directional prior μ₀ = c·a; c is the exploration knob |

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

| Method | Description |
|--------|-------------|
| **max_elbo** | Select single best fit by ELBO |
| **uniform** | Average PIPs equally across all fits |
| **elbo_softmax** | ELBO-softmax weighted average of PIPs |
| **cluster_weight** | JSD-cluster fits (threshold 0.15), ELBO-softmax within clusters, uniform across clusters |
| **cluster_weight_050** | Same as above but with JSD threshold 0.50 (coarser clustering) |
| **oracle** | Best individual fit by AUPRC (truth-aware; reference only) |
