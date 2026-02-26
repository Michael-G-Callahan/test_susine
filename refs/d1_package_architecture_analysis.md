# D1 Package Architecture Analysis: susine vs susieR

**Date:** 2025-02-25 (original); **Revised:** 2025-07-07
**Purpose:** Concrete implementation plans for D1, updated after discovering that the upstream stephenslab/susieR has undergone a complete rewrite to v2.0.0. The original analysis was based on the stale local fork (v0.12.40) and is superseded.

> **KEY CORRECTION:** The original analysis stated "susieR does NOT have: alpha/PIP convergence, SuSiE-ash, SuSiE-inf." This was **wrong** — it reflected the stale local fork. The upstream susieR 2.0.0 has **all of these features** and more. See `refs/susieR_2.0_inventory.md` for the full inventory.

---

## Package Inventories

### susieR — upstream stephenslab/susieR v2.0.0

- **~22 R source files + Rcpp/C++ backend** (Mr.ASH via `caisa_rcpp`)
- **Complete architectural rewrite** by Alex McCreight
- **S3 class dispatch**: 3 data classes (`individual`, `ss`, `rss_lambda`) with ~20+ generic methods
- **Pipeline**: Interface → Constructor → Workhorse → IBSS → Backend (S3)
- **Deps**: Rcpp, RcppArmadillo, matrixStats, Matrix

**Feature inventory:**

| Feature | Status | Data support |
|---------|--------|--------------|
| Standard SuSiE (IBSS) | ✅ | individual, SS, RSS |
| SuSiE-ash (Mr.ASH background) | ✅ `unmappable_effects = "ash"` | **individual only** — rejects SS/RSS |
| SuSiE-inf (infinitesimal background) | ✅ `unmappable_effects = "inf"` | individual + SS — **not RSS (λ>0)** |
| PIP convergence | ✅ `convergence_method = "pip"` | all 3 |
| Refinement | ✅ `refine = TRUE` → `run_refine()` | all 3 |
| NIG/Servin-Stephens prior | ✅ `estimate_residual_method = "Servin_Stephens"` | individual + SS — not RSS |
| Stochastic LD correction | ✅ `stochastic_ld_sample = B` | RSS (λ=0 path) |
| EB prior variance (optim/EM/simple) | ✅ | all 3 |
| Non-uniform pi (`prior_weights`) | ✅ | all 3 |
| Warm starts (`model_init`) | ✅ | all 3 |
| MoM residual variance | ✅ `estimate_residual_method = "MoM"` | individual + SS |
| `susie_auto()` L-doubling | ✅ | all 3 |
| `susie_trendfilter()` | ✅ | individual |

**Key constraints:**
- Refinement is **incompatible** with unmappable effects (`refine = TRUE` + `unmappable_effects != "none"` errors)
- When `unmappable_effects ∈ {"inf", "ash"}`, convergence is forced to PIP-based (not ELBO)
- SuSiE-ash requires individual-level data (X, y) and uses compiled C++ (caisa_rcpp)
- NIG prior forces `convergence_method = "pip"` when L > 1

**What susieR 2.0 does NOT have:**
- Nonzero prior means (μ₀ ≠ 0)
- Annotation-driven priors
- Multi-fit aggregation (single-fit output only)
- Annealing / temperature schedules

### susine (v0.0.0.9000) — unchanged

- **16 R source files**, ~2,850 lines
- **Zero dependency on susieR** — fully standalone reimplementation
- Implements: `susine()`, `susine_ss()`, `susine_rss()`
- Has: nonzero prior means (mu_0 = c·a), warm starts (`init_alpha`), annealing, EB prior updates (mean/var/both/none with rejection), aggregation (equal + ELBO-softmax), auto-scale mu_0/sigma_0_2, non-uniform pi (`prior_inclusion_weights`)
- Does NOT have: refinement, alpha/PIP convergence, SuSiE-ash/inf, advanced LD regularization, NIG prior, stochastic LD

### Local susieR fork (v0.12.40) — STALE

The local fork is 2+ years behind upstream. It lacks all features added in the 2.0 rewrite. **The local fork should be replaced with upstream susieR 2.0** (installed from GitHub) to serve as the SuSiE baseline.

### Structural relationship

susine is a **clean reimplementation** of the SuSiE IBSS algorithm with SuSiNE extensions. Same algorithmic template (build full residuals → for each L, add-back → SER → subtract) but different internal data structures (data.frame with list-columns vs flat list), different helper functions, and additional machinery (annealing, EB mean, warm starts, auto-scaling). No shared source files. The upstream susieR 2.0 has a much more sophisticated architecture (S3 dispatch, Data/Params/Model objects) that is fundamentally different from both the old susieR and susine.

---

## Plan A (Revised): Two-Package Strategy

**Use upstream susieR 2.0 as-is for SuSiE baselines. Keep susine for SuSiNE extensions. Harness-level BFS refinement for both.**

This is the natural evolution of the original Plan A, but dramatically simplified because upstream susieR 2.0 already has most features we previously thought we'd need to implement.

### What upstream susieR 2.0 provides for free

| Feature | Previously... | Now... |
|---------|---------------|--------|
| Refinement | Needed porting | ✅ Built into susieR: `refine = TRUE` |
| PIP convergence | Needed implementing | ✅ Built into susieR: `convergence_method = "pip"` |
| SuSiE-ash | Would skip (D3) | ✅ Available: `unmappable_effects = "ash"` (individual only) |
| SuSiE-inf | Would skip (D3) | ✅ Available: `unmappable_effects = "inf"` (individual + SS) |
| NIG prior | Not considered | ✅ Available: `estimate_residual_method = "Servin_Stephens"` |
| EB prior variance | Had it (old susieR) | ✅ Enhanced: optim/EM/simple methods |
| Stochastic LD | Not considered | ✅ Available for RSS path |

### What still needs implementing

| Feature | Where | Effort | Notes |
|---------|-------|--------|-------|
| **BFS refinement** | test_susine (harness) | Low–Med (~80–120 lines) | Model-agnostic; works with both `susine()` and `susie()` |
| **Alpha/PIP convergence in susine** | susine package | Low (~30–50 lines) | For cross-method comparability; susieR 2.0 already has this |
| **LD regularization in susine** | susine package | Low–Med (~100 lines) | Only for Phase D real data |
| **Aggregation pipeline** | test_susine (harness) | Med (~150–200 lines) | Method A cluster-then-weight; model-agnostic |

### A1. BFS Refinement in test_susine (unchanged from original)

The harness-level BFS tree exploration strategy remains the right approach:
- All fits are saved for the flatten-then-aggregate pipeline (D6), not greedily discarded
- Budget K is enforced at the harness level (D8)
- Model-agnostic — works with `susine()` and `susie()` alike
- Complements (rather than replaces) susieR's built-in `refine = TRUE`

For **SuSiE baselines**, we can also benchmark susieR's built-in `refine = TRUE` as a separate exploration arm. This provides a direct comparison: susieR's greedy refinement (keeps single best) vs our BFS refinement (saves all fits for ensembling).

**Algorithm:** (unchanged — see below)
```
refine_bfs(model_fn, data, K, cs_coverage, cs_min_purity):
  fit_0 ← model_fn(data, prior_inclusion_weights = uniform)
  fits ← [fit_0]
  queue ← [fit_0]

  while length(fits) < K and queue is not empty:
    parent ← dequeue(queue)
    cs_list ← get_credible_sets(parent, coverage, min_purity)

    for cs in cs_list:
      if length(fits) >= K: break
      pi_new ← parent's prior_inclusion_weights
      pi_new[cs$snp_indices] ← 0
      pi_new ← pi_new / sum(pi_new)
      fit_branch ← model_fn(data, prior_inclusion_weights = pi_new)
      fits ← append(fits, fit_branch)
      queue ← append(queue, fit_branch)

  return fits
```

### A2. Alpha/PIP Convergence in susine

**susieR 2.0 now has this**: `convergence_method = c("elbo", "pip")` with `check_convergence()` utility.

Still implement in susine for parity:
- Add `convergence_method` parameter: `"elbo"` (default) or `"alpha"`
- Track `max(abs(alpha_new - alpha_old))` across all L effects each iteration
- **Effort: Low (~30–50 lines)**

### A3. SuSiE-ash / SuSiE-inf — now available upstream (revises D3)

**Major change from original analysis.** These are fully implemented in susieR 2.0:

- **SuSiE-ash**: Mr.ASH background component with masking/unmasking logic, C++ compiled coordinate ascent. Individual-level data only.
- **SuSiE-inf**: Infinitesimal polygenic background (Cui et al. 2024 Nature Genetics). Individual + SS data. Omega-weighted SER updates, BLUP theta estimation.

**Implications for D3:**
- We can run SuSiE-ash and SuSiE-inf **as baseline comparison arms** with zero implementation effort — just call `susie(X, y, unmappable_effects = "ash")` or `susie(X, y, unmappable_effects = "inf")`
- The study plan's suggestion of a "compact failure-mode comparison vs SuSiE-ash/SuSiE-inf" becomes trivially achievable
- No need to implement these in susine — they are SuSiE extensions that address a different problem (background effects) than SuSiNE (directional annotation priors)

**Constraint**: Both require individual-level data. Our simulations use individual-level data, so this is fine. Real-data Phase D may need to use RSS, where ash/inf are unavailable.

### A4. Advanced LD Regularization

**susieR 2.0 now has**:
- Eigenvalue decomposition and clipping in constructors
- Stochastic LD correction (`stochastic_ld_sample = B`) for RSS with sketch matrices
- NIG prior for small-sample robustness
- Ridge lambda parameter in `rss_lambda` class

For susine: still port basic PSD projection and `estimate_s_rss()` for Phase D real data. Low priority for simulations.

**Effort: Low–Med (~100 lines)**

### Plan A Summary (Revised)

| Feature | Effort | Where | Priority | Notes |
|---------|--------|-------|----------|-------|
| **Install upstream susieR 2.0** | **Trivial** | local env | **IMMEDIATE** | `remotes::install_github("stephenslab/susieR")` |
| BFS refinement | Low–Med | test_susine | HIGH | Harness-level; model-agnostic |
| Alpha/PIP convergence | Low | susine | MEDIUM | For parity with susieR 2.0 |
| SuSiE-ash/inf baselines | **Zero** | — | MEDIUM | Just call susieR 2.0 with `unmappable_effects` |
| LD regularization | Low–Med | susine | LOW (sim) / MED (real) | Phase D only |
| Aggregation (Method A) | Med | test_susine | HIGH | Cluster-then-weight |

**Total new code: ~260–370 lines** (BFS refinement + alpha convergence + LD reg). Plus aggregation pipeline (~150–200 lines) which is needed regardless of plan.

---

## Plan B (Revised): Port SuSiNE Into upstream susieR 2.0

### What changed from the original Plan B

The original Plan B was based on the stale v0.12.40 fork — a monolithic codebase with flat lists, 3 separate SER files, and straightforward (if tedious) modification. **susieR 2.0 is a completely different beast:**

- **S3 dispatch across 3 data classes**: Every SER-related computation goes through generic methods (`compute_ser_statistics`, `calculate_posterior_moments`, `compute_kl`, etc.). Modifying the posterior mean requires changes in **every S3 method implementation** across all 3 data classes.
- **Compiled C++ backend**: Mr.ASH uses `caisa_rcpp`. If nonzero prior means interact with ash (they shouldn't, but the masking logic might), C++ code may need modification.
- **Complex initialization**: The `initialize_susie_model` generic has class-specific implementations that would all need to accept and store `mu_0`.
- **Constructor pipeline**: Data flows through Interface → Constructor → Workhorse → IBSS, with validation at each stage. The `mu_0` parameter would need threading through all 4 stages.

### B1. Nonzero Prior Means (μ₀ = c · a) in susieR 2.0

**Changes needed (much more than original estimate):**

1. **S3 generic methods** (all 3 data class implementations):
   - `compute_ser_statistics.{individual,ss,rss_lambda}` — modified Bayes factor with non-central null
   - `calculate_posterior_moments.{individual,ss,rss_lambda}` — μ₁ formula changes
   - `compute_kl.{individual,ss,rss_lambda}` — KL with nonzero prior mean
   - `SER_posterior_e_loglik.{individual,ss,rss_lambda}` — expected log-likelihood changes
   - `loglik.{individual,ss,rss_lambda}` / `neg_loglik.{individual,ss,rss_lambda}` — for EB optim
   - That's **15+ method implementations** to modify (5 generics × 3 classes)

2. **Prior variance optimization** — `optimize_prior_variance()` and its methods:
   - Currently assumes zero-mean prior in the objective function
   - The optim/EM/simple estimators all need the nonzero-mean correction terms

3. **IBSS loop** — `ibss_fit()`:
   - Pass `mu_0[l,]` to `single_effect_regression()` for each effect l

4. **Constructors** — all 3 constructors need to accept, validate, and store `mu_0`

5. **Interface functions** — `susie()`, `susie_ss()`, `susie_rss()`:
   - Add `prior_mean` parameter, thread through workhorse

6. **Tests**: susieR 2.0 has **>1000 tests**. All must still pass after μ₀ changes.

**Effort: HIGH (~500–700 lines of changes across 15+ files, in an unfamiliar S3 architecture)**

### B2–B5: Other features (unchanged assessment)

- Warm starts: susieR 2.0's `model_init` already covers this. **Low effort**.
- EB prior mean: Not needed (D13).
- Aggregation: Model-agnostic; lives in harness regardless. **Not a Plan B factor**.

### Plan B Summary (Revised)

| Feature | Effort | Notes |
|---------|--------|-------|
| Nonzero prior means in susieR 2.0 | **HIGH** (~500–700 lines) | 15+ S3 method implementations; unfamiliar architecture; >1000 tests to maintain |
| Warm start helper | Low (~50 lines) | `model_init` already exists |
| Aggregation | Lives in harness regardless | Not a differentiator |

**Total for Plan B: ~550–750 lines across 15+ files in a complex S3 architecture you didn't write.**

---

## Head-to-Head Comparison (Revised)

| Dimension | Plan A (two-package) | Plan B (SuSiNE → susieR 2.0) |
|-----------|---------------------|-------------------------------|
| **New code** | ~260–370 lines | ~550–750 lines |
| **Files touched** | 3–5 (mostly new files in test_susine + susine) | 15+ (modifying S3 methods across 3 data classes) |
| **Risk** | Low — adding self-contained features to your own code | **Very high** — modifying the mathematical core of a sophisticated S3 architecture with compiled C++ |
| **Testing** | Existing susine tests + new refinement tests | >1000 existing susieR tests must pass; need new μ₀-specific tests |
| **Upstream drift** | N/A — susine is independent; susieR installed from upstream as-is | **Critical** — any upstream susieR update would conflict with your μ₀ modifications |
| **SuSiE-ash/inf access** | ✅ Free — just call susieR 2.0 | ✅ Also free — but nonzero means + ash interaction is undefined and potentially dangerous |
| **Refinement** | Built-in `refine = TRUE` for SuSiE baselines; BFS in harness for both | Same built-in; BFS still in harness |
| **Harness impact** | Minimal — test_susine already calls `susine()` for SuSiNE, adds `susie()` calls for baselines | **Significant** — retool all SuSiNE calls from `susine()` to `susie()` |
| **Familiarity** | You wrote susine; susieR 2.0 is read-only | susieR 2.0's S3 internals are someone else's complex architecture |
| **"Drop-in" claim** | Two packages, but the harness wraps both — user-facing interface is the harness anyway | Single `susie()` call, but the μ₀ fork drifts from upstream |
| **Time to first pilot** | **Days** — install susieR 2.0, write BFS, wire harness | **Weeks** — understand S3 architecture, modify 15+ methods, pass >1000 tests |

---

## Recommendation: Plan A (Two-Package Strategy)

**Plan A is now a clear winner.** The discovery of susieR 2.0 fundamentally changes the calculus:

1. **Most of the "missing features" are already in susieR 2.0.** We don't need to implement refinement, PIP convergence, ash, inf, NIG priors, or stochastic LD. Just install the upstream package.

2. **Plan B's effort estimate roughly doubled** because the S3 architecture is far more complex than the old monolithic code. Porting nonzero prior means now means modifying 15+ S3 method implementations (5 generics × 3 data classes) instead of editing 3 flat SER files.

3. **Plan B creates a permanent upstream-tracking burden.** susieR 2.0 is under active development. Any fork with μ₀ modifications diverges immediately. With Plan A, we install susieR from upstream and get updates for free.

4. **The harness is the true user interface.** Neither `susine()` nor `susie()` is called directly by practitioners — the harness (test_susine) wraps both. The "drop-in replacement" argument for Plan B is moot.

5. **SuSiE-ash/inf become free comparison arms.** The paper can include a "compact failure-mode comparison vs SuSiE-ash/SuSiE-inf" (study plan Phase A) with zero implementation cost, strengthening the story.

### Concrete implementation path

```
Step 1: Install upstream susieR 2.0
  remotes::install_github("stephenslab/susieR")

Step 2: Add alpha/PIP convergence to susine (~30–50 lines)

Step 3: Implement BFS refinement in test_susine (~80–120 lines)

Step 4: Wire harness to call both packages:
  - susie()   for: SuSiE, SuSiE+refinement, SuSiE-ash, SuSiE-inf baselines
  - susine()  for: SuSiNE with c-grid, sigma_0^2 grid, restarts
  - Both via BFS refinement (harness-level, model-agnostic)

Step 5: Implement aggregation pipeline in test_susine (~150–200 lines)

Step 6: LD regularization in susine (Phase D, later)
```

### What this means for D3

**D3 should be RESOLVED: include SuSiE-ash and SuSiE-inf as baseline comparison arms.**

The original "leaning skip" was based on the (incorrect) belief that implementation would cost 500+ lines each. Since they're built into susieR 2.0, the cost is **zero** — one function argument. The paper benefits from showing that directional priors (SuSiNE) address a complementary problem to background-effect models (ash/inf), and that ensemble methods work across both paradigms.

**Constraint**: ash/inf require individual-level data. Our simulations use individual-level data, so this works. Phase D real data uses RSS where ash/inf are unavailable — note as limitation.
