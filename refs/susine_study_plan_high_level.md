# SuSiNE paper: high-level study plan (handoff)

**Target journal:** Bioinformatics  
**Scope:** eQTL fine-mapping (typical cis window), individual-level genotype/phenotype where available, plus a summary-statistic implementation as secondary.

---

## 0) One-paragraph framing (for anyone new)

SuSiE is a fast variational fine-mapper that approximates a posterior over sparse multi-SNP models using a “sum of single effects” factorization and an IBSS-style coordinate ascent algorithm. This speed comes with a known limitation: **the posterior landscape can be multimodal** (multiple statistically plausible causal configurations under LD), and a single variational optimum cannot represent multiple well-separated modes simultaneously. Additionally, *hyperparameters* (e.g., prior variance) are often uncertain/unstable under empirical Bayes, causing large changes in inferred CS/PIPs. SuSiNE extends SuSiE by allowing **nonzero prior means** for effect sizes (a directional functional prior), using a single signed annotation vector **a** and a global scale **c**:  
\[
\mu_0 = c\,a.
\]
Because **c is not identifiable** from the data in realistic settings, we treat it as an **interpretable exploration knob**—a principled way to traverse posterior basins—then **aggregate** across model fits under a fixed compute budget using BMA-style weighting (ELBO-softmax) and/or robust alternatives. We show (i) why fine-mapping with SuSiE is multimodal and hyperparameter-sensitive, (ii) how ensembling improves SuSiE, (iii) how SuSiNE’s additional knob improves exploration and aggregation beyond SuSiE’s knobs, (iv) how these insights transfer to real eQTL loci, and (v) how **dataset-level failure mode diagnostics** can be built by screening multiple z-score-only metrics and selecting a primary one to pair with LD-only \(M_1\).

---

## 1) Contributions (what the paper should *clearly* claim)

### 1.1 Method contributions
1. **Directional prior means for sparse fine-mapping (SuSiNE):** incorporate signed functional information through a nonzero mean effect-size prior, \(\mu_0=c a\), shared across effects.
2. **Principled workflow for non-identifiable hyperparameters:** treat \(c\) (and \(\sigma_0^2\)) as exploration knobs; avoid unstable EB point estimates.
3. **Budgeted ensembling and aggregation:** fit multiple candidate models per locus, then aggregate PIPs across models using BMA-style weights (primary: ELBO-softmax), with sensitivity analyses to alternative weighting schemes.
4. **Diagnostics for *why* it helps:** quantify (a) posterior diversity / multimodality exploration, (b) quality of discovered modes, and (c) the ability of weighting to recognize good modes.
5. **Failure-mode characterization framework:** systematically study a candidate family of z-score-only difficulty metrics, select the best-performing primary z-metric, then define a two-axis, truth-agnostic map using LD-only **M1** + that primary z-metric to stratify expected SuSiE behavior.

### 1.2 Empirical contributions
1. A simulation suite for eQTL-like loci showing:
   - SuSiE multimodality + instability under hyperparameter uncertainty,
   - ensemble gains for SuSiE,
   - additional gains when allocating budget to SuSiNE’s annotation-scale grid.
2. Real-data eQTL case studies on multiple loci spanning difficulty (LD complexity, z-score concentration/kurtosis), including ensemble outputs and biologically interpretable differences vs vanilla SuSiE.
3. A practical stratification recipe for practitioners: given a new locus, compute M1 + z-score diagnostics first, then decide whether default SuSiE is sufficient or whether heavy multi-start / grid exploration is required.

---

## 2) Core hypotheses (what we need to demonstrate)

### H1 — Multimodality and instability are practical, not pathological edge cases
- Under realistic LD + oligogenic architectures, SuSiE solutions vary substantially across initialization and prior variance settings.
- EB estimation of prior variance (and any EB attempt to infer \(c\)) is unstable due to circularity/non-identifiability.

### H2 — Ensemble workflows can improve SuSiE in a predictable way
- Given a fixed budget \(K\) model fits, SuSiE ensembles (random starts, \(\sigma_0^2\) grids, optionally refinement) yield improved PIP calibration / power-FDR tradeoffs vs a single fit.

### H3 — SuSiNE’s \(c\)-grid provides **structured exploration** that random restarts cannot replicate
- Define diversity metrics \(D\) (e.g., mean pairwise JSD across PIP vectors).
- For matched compute (same \(K\)), expect:
  - \(D_{c\text{-grid}} \gg D_{restart}\) for many datasets, because changing \(c\) changes the objective function (inductive bias), not just initialization noise.

### H4 — Gains are not “just more compute”
- If SuSiNE uses a larger hyperparameter grid than SuSiE, show:
  - at equal budget \(K\), allocating fits to \(c\)-grid yields better (a) diversity, (b) best-mode quality, and/or (c) aggregation quality than allocating all fits to restarts.

### H5 — Directional priors are most valuable when fine-mapping is hardest
- Under easy sparse settings, SuSiNE may offer little improvement even with good annotations.
- Under difficult/ambiguous settings (LD complexity, diffuse background, multiple plausible models), directional priors reduce ambiguity and improve inference.

### H6 — Failure modes are predictable from LD-only + z-score-only diagnostics
- M1 and z-score shape metrics explain a meaningful portion of between-dataset variation in the following operational failure endpoints (matching current code outputs):
  - **Instability / multimodality endpoints** from `compute_multimodal_metrics` (`R/run_model.R`):
    - `mean_jsd`, `median_jsd`, `max_jsd` (pairwise Jensen-Shannon divergence across fit PIP vectors),
    - `jaccard_top10` (mean pairwise Jaccard overlap of top-10 PIP SNP sets),
    - `mean_pip_var` (mean per-SNP variance of PIPs across fits),
    - `n_clusters` (hierarchical clustering count at JSD cutoff `jsd_threshold`, default `0.171661`).
  - **Baseline fit quality endpoints** from `evaluate_model` (`R/evaluation_metrics.R`):
    - `AUPRC` (primary),
    - plus sensitivity endpoints such as `power`, `mean_purity`, `mean_size`, and `cross_entropy`.
  - **Structured-exploration gain endpoints**:
    - between-strategy deltas/ratios for the above metrics (e.g., \( \Delta \)AUPRC and \( \Delta \)mean/max JSD for SuSiNE \(c\)-grid vs restart-only exploration at matched budget).
- A two-axis stratification (M1, Z-metric) gives actionable pre-fit guidance for expected fit quality and workflow choice.

---

## 3) Study structure (high-level)

### Cross-cutting pillar: failure mode characterization (key contribution)
Treat this as a first-class contribution, not a side analysis.

1. Build a two-axis dataset map:
   - x-axis: LD-only metric \(M_1\)
   - y-axis: selected z-score-only metric $Z_{\text{primary}}$
2. Define failure endpoints per dataset:
   - instability across restarts / prior variance choices
   - posterior multimodality (pairwise JSD, cluster count)
   - gap between single-fit and ensemble performance
3. Quantify how endpoints vary across $(M_1, Z_{\text{primary}})$ bins.
4. Use this map to recommend a workflow:
   - low-difficulty region: single-fit or light exploration
   - high-difficulty region: structured exploration (variance grid / \(c\)-grid) + aggregation

### Phase A: Establish SuSiE pathology + motivate ensembling
**Goal:** show that single-fit SuSiE is not a reliable “one-shot inference” tool when the posterior is multimodal and hyperparameters are uncertain.

#### A0) Illustrative Vignette (Pathological Examples)
**Goal:** Create a clean R notebook vignette demonstrating SuSiE failure modes on synthetic eQTL-like data ($X \in \{0,1,2\}$, centered/scaled).

1. **Baseline Success:**
   - $n=600, p=10$ (or small $p$), 3 causal effects.
   - SuSiE correctly identifies all 3 in separate CSs.
2. **Correlated "Bandits" (Ambiguity):**
   - 3 causal variables, each with a highly correlated "bandit" (singleton proxy).
   - Expected: 3 CSs, each containing 2 variables (causal + bandit).
3. **Pathology 1: Identical Subspaces (Multimodality):**
   - **Math/Logic:** Construct two triplets $A=\{x_1, x_2, x_3\}$ and $B=\{x_4, x_5, x_6\}$ such that their column spaces are nearly identical ($\text{span}(A) \approx \text{span}(B)$), but variables $x_i$ and $x_j$ are not individually highly correlated across sets.
   - **Implementation:** Let $A$ be random orthogonal vectors. Let $B = A Q$ where $Q$ is a random orthogonal rotation matrix. Then mixing signals $\beta_A$ generate $y$. The model $y \approx B (Q^T \beta_A)$ explains data equally well.
   - **Expected:** SuSiE converges to *one* of the two sets (A or B) depending on initialization variance/random seed, but rarely mixes them or captures the multimodal uncertainty without multiple fits.
4. **Pathology 2: The "Ultimate" Greedy Trap (Cancellation):**
   - **Math/Logic:** Create a signal where marginal correlations vanish but joint signal is strong.
     - Causal variables $x_1, x_2$ such that $y = \beta(x_1 - x_2) + \epsilon$.
     - If $x_1, x_2$ are highly correlated ($r \approx 1$), then $\text{Var}(x_1 - x_2)$ is small, but if scaled up, they fit $y$.
     - Crucially, marginal $X^T y$ for $x_1, x_2$ can be small.
     - Add a "Decoy" $x_d$ independent of $x_1, x_2$ with moderate marginal correlation with $y$ (just noise alignment or weak true effect).
   - **Expected:** SuSiE's greedy forward selection (L=1 pass) picks $x_d$ because $|\text{cor}(x_d, y)| > |\text{cor}(x_1, y)|$. Subsequent passes fail to swap $x_d$ for $\{x_1, x_2\}$ because swapping one at a time doesn't improve the objective (you need to add both $x_1, x_2$ simultaneously to see the benefit).
   - **Demo:** Show ELBO($x_1, x_2$) $\gg$ ELBO($x_d$).


**Deliverables (main text):**
1. Random-start variability: same dataset, multiple seeds → different solutions.
2. Prior-variance sensitivity: \(\sigma_0^2\) grid → multiple clusters of solutions.
3. EB instability: EB underperforms reasonable fixed values, and is sensitive to initialization/data idiosyncrasies.
4. Link instability to observable dataset difficulty metrics (M1 + z-score metrics).
5. Build a two-axis failure map (M1 on x-axis, chosen z-score metric on y-axis) and label dominant failure mode regions.
6. Quantify how failure mode region changes expected benefit of random restarts, \(\sigma_0^2\) grids, and SuSiNE \(c\)-grids.

**Key conceptual output from A0:** The vignette establishes two orthogonal failure axes: (1) *model-specification barriers* (OR-of-ANDs uncertainty that no single fit can represent → ensembling required) vs (2) *optimizer barriers* (better basin exists but greedy updates miss it → multi-start selection sufficient). It also demonstrates that EB estimation of $\sigma_0^2$ distorts ELBO-softmax ensemble weights even for symmetric basins (0.21/0.79 long-run weights vs 0.47/0.53 with fixed prior variance), motivating the fixed-prior-variance grid approach in Phase B.

**SuSiE-ash / SuSiE-inf instability comparison (supplement, previously optional — now planned):**
- D3 RESOLVED: SuSiE-ash and SuSiE-inf are fully implemented in upstream susieR 2.0 (`unmappable_effects = "ash"` / `"inf"`) at zero implementation cost. Include as baseline comparison arms in Phase B (see B1 below). Within Phase A, include a compact failure-mode comparison showing that multimodality/instability is not unique to vanilla SuSiE.

---

### Phase B: Ensemble study — exploration + aggregation under a fixed budget
**Goal:** quantify how different exploration knobs and aggregation rules trade off performance.

#### B1) Define “use cases” (model families)

**[UPDATED 2025-07-07 per D1/D3 resolution — see `analysis_completion_status.md`]**

5 model arms (factored design per D2):

| Arm | Prior spec | Package | Notes |
| --- | ---------- | ------- | ----- |
| 1. Vanilla SuSiE | mu=0, fixed sigma_0^2, uniform pi | susieR 2.0 | Primary baseline |
| 2. SuSiE + functional pi | mu=0, fixed sigma_0^2, pi=softmax(\|a_j\|/tau) | susieR 2.0 | D5: shows pi-only is different from mean-shifts |
| 3. SuSiNE | mu_0 = c*a, fixed sigma_0^2, uniform pi | susine | Central method |
| 4. SuSiE-ash | vanilla SuSiE + Mr.ASH background | susieR 2.0 | D3: `unmappable_effects = "ash"`, individual data only |
| 5. SuSiE-inf | vanilla SuSiE + infinitesimal background | susieR 2.0 | D3: `unmappable_effects = "inf"`, individual + SS |

*(Central comparison remains Arms 1 vs 3. Arms 4-5 are free baseline comparison arms showing that background-effect models address a complementary problem to directional priors. Arm 2 isolates the value of mean-shifts vs pi-only functional information.)*

#### B2) Exploration mechanisms (how we generate \(K\) candidate fits per locus)
Under a fixed budget \(K\):
1. **Random restarts** (SuSiE): vary initialization of \(\alpha\) (e.g., Dirichlet).
2. **Prior variance grid** (SuSiE): vary \(\sigma_0^2\) over a small grid.
3. **SuSiE refinement** (if used): treat refinement as an exploration sampler with its own cost accounting.
4. **Annotation-scale grid** (SuSiNE): vary \(c\) over a small interpretable grid (optionally include a shrinkage exponent \(\gamma\) only if already in code; otherwise omit to keep scope tight).
5. **Combined grids** (secondary): \((c, \sigma_0^2)\) grid for SuSiNE if needed for robustness, but note the compute fairness issue explicitly.

#### B3) Aggregation mechanisms (how we combine candidate fits)
Primary:
- **ELBO-softmax weighting** within a use case and dataset.

Secondary/sensitivity:
- Uniform weights
- Max-ELBO (select best model)
- Simple PIP average
- Optional: stacking/LOO-style weights if you can produce per-observation predictive densities reliably.

**Caution — EB prior variance estimation distorts ELBO-softmax weights:**
The pathology vignette (A0, Scenario 3) demonstrates that `estimate_prior_variance = TRUE` amplifies small numerical asymmetries between symmetric basins, producing long-run ensemble weights of 0.21/0.79 instead of the expected ~0.50/0.50. Fixing prior variance restores balance. This motivates using a fixed $\sigma_0^2$ grid (B2) rather than per-fit EB estimation when computing ELBO-softmax weights for aggregation.

**Aggregation design decision (RESOLVED — see `analysis_completion_status.md` D6):**
- All aggregation methods operate on a **flattened pool** of fits for a given (dataset, model_spec). Each fit's ELBO is evaluated under its own prior only (no cross-prior re-scoring).
- Primary methods: Max ELBO, Uniform, ELBO softmax, and **Cluster-then-weight (Method A)** — which importance-corrects for optimizer frequency bias via JSD-based clustering.
- Grid structure matters at the *exploration* stage (what fits to generate), not at the *aggregation* stage.

---

### Phase C: Decompose *why* SuSiNE helps (mode discovery vs mode quality vs weighting)
For each dataset, compare ensembles on three axes:

1. **Exploration / diversity:** do we find more distinct basins?
   - pairwise distances between PIP vectors (JSD, correlations)
   - clustering of fits (number of clusters above a JSD cutoff)
   - top-k Jaccard overlap of high-PIP SNPs
2. **Best-mode quality:** among the \(K\) candidate fits, how good is the single best PIP vector?
   - AUPRC (simulation)
   - power at FDR threshold
3. **Weighting / recognition:** if we pool candidate models across methods, which weighting scheme best selects/weights good modes?
   - “pooled modes” experiment: pool all fits, then apply different weighting rules to see who recognizes good modes when they exist.

---

### Phase D: Real data eQTL study (multiple loci)
**Goal:** take the simulation-learned workflow and apply to real loci, emphasizing interpretability and robustness.

Steps:
1. Screen a large set of loci/genes; compute difficulty metrics:
   - LD-only metric M1 (colinearity proxy)
   - z-score-only candidate metrics (concentration / tail-heaviness / effective signal count)
2. Select **3–4 additional loci** representative across the difficulty spread.
3. For each locus:
   - run the chosen ensemble workflows
   - report multimodality metrics and final aggregated PIPs
   - compare SuSiE vs SuSiNE; highlight SNPs where SuSiNE changes prioritization and give biological context.

Also: include a minimal validation that the summary-statistic implementation matches individual-level results on at least one shared dataset (or remove the claim).

---

## 4) Simulation design (high-level spec)

### 4.1 Data generation
- Genotypes: use real LD / genotype matrices if available; n≈600, p≈1000 per locus.
- Phenotype architecture grid:
  - \(p^*\) (# causal), e.g. 1–10 (include oligogenic)
  - total \(h^2\): vary (important robustness dimension)
  - effect sizes: start with Normal for causal effects; optionally include heavy-tail / mixture in sensitivity analyses.
- Annotation model:
  - one signed vector \(a\in\mathbb R^p\)
  - control annotation quality via \(R^2(a, \beta)\) or an equivalent alignment metric
  - include scenarios with imperfect directionality (optional: sign flips) only if needed.

### 4.2 Dataset difficulty metrics (computed per dataset before fitting)
LD-based:
- **Primary LD-only metric (M1):**
  \[
  M_1 = \frac{2}{p(p-1)}\sum_{i<j}|r_{ij}|(1-|r_{ij}|).
  \]
  (Normalization optional if all loci use similar \(p\); keep normalized form for cross-locus comparability.)

z-score-based (candidate family; summary-stat only and truth-agnostic):
- **Top-k concentration ratio(s)** `z_topk_ratio` (default `top_k = 10` in current code):
  \[
  C_k = \frac{\sum_{i\in top\text{-}k} z_i^2}{\sum_{i\notin top\text{-}k} z_i^2}, \quad k\in\{1,5,10\}.
  \]
- **Peak magnitude** `z_max_abs`: \(Z_{\max}=\max_i |z_i|\).
- **Tail count** `z_count_abs_gt_3`: \(N_{>|3|}=\#\{i:|z_i|>3\}\).
- **Effective signal count** `z_eff_signals`:
  \[
  S_{\mathrm{eff}}=\frac{(\sum_i z_i^2)^2}{\sum_i z_i^4}.
  \]
- **Optional shape metrics (if stable):**
  - kurtosis of \(|z|\) or \(z^2\),
  - entropy of normalized \(z_i^2\): \(H_z=-\sum_i p_i\log p_i,\; p_i=z_i^2/\sum_j z_j^2\).

Metric selection protocol (to end with one primary z-metric in main text):
1. Compute all candidate z-metrics for each dataset/locus.
2. Evaluate each metric’s ability to predict instability/diversity endpoints (e.g., mean pairwise JSD, cluster count, seed sensitivity) using simple monotonic fits or rank correlations.
3. Choose one **primary z-metric** (best predictive + most interpretable + stable across simulation/real data), keep others as supplement/sensitivity.
4. Define stratification bins (e.g., terciles) on $(M_1, Z_{\text{primary}})$ and report performance/instability by bin.

### 4.3 Model runs per dataset (skeleton)
For each dataset \((X,y,a)\):
1. Compute difficulty metrics.
2. Run candidate fits for each exploration method with budget \(K\).
3. Aggregate PIPs per use case under each aggregation rule.
4. Compute:
   - performance metrics (simulation truth-aware)
   - multimodality metrics (truth-agnostic)
5. Save everything (raw fits + aggregated outputs).

---

## 5) Addressing the “compute confound” critique (pilot study)

**Question:** Are SuSiNE gains due to directional priors or simply because SuSiNE’s hyperparameter grid creates more distinct fits?

**Pilot design (supplement-level):**
- Select 5–10 datasets spanning difficulty (M1 x $Z_{\text{primary}}$, or M1 x max|z| during pilot screening).
- Fix \(\sigma_0^2\) (e.g., 0.2).
- Run:
  - SuSiE: \(K\) random restarts (e.g., 12) with Dirichlet-initialized \(\alpha\).
  - SuSiNE: \(K\) values of \(c\) (e.g., 0, 0.2, 0.4, 0.6, 0.8, …) at same \(\sigma_0^2\).
- Compare diversity metrics:
  - mean/median/max pairwise JSD
  - cluster count
  - PIP variance
  - PIP correlation matrices (heatmaps)
- Report: fraction of datasets with \(D_{c\text{-grid}} > D_{restart}\) and mean ratio \(D_{c\text{-grid}}/D_{restart}\).

**Interpretation plan:**
- If \(D_{restart} \ll D_{c\text{-grid}}\): argue SuSiNE offers structured exploration not replicable by restarts.
- If \(D_{restart} \approx D_{c\text{-grid}}\): either match compute more explicitly or reframe gains around interpretability/weighting rather than exploration.

---

## 6) Primary outputs (figures/tables checklist)

### Simulations (main)
1. SuSiE pathology: solution variability vs seed and \(\sigma_0^2\); EB instability.
2. Failure map: heatmap/contours of instability and multimodality across $(M_1, Z_{\text{primary}})$.
3. Ensemble gains: performance (AUPRC, power @ FDR) vs ensemble budget \(K\) by exploration method, stratified by failure-map region.
4. Why it helps: diversity metrics vs performance; pooled-modes weighting analysis.

### Pilot (supplement)
- restart diversity vs c-grid diversity:
  - violin/boxplots
  - head-to-head scatter with diagonal
  - correlation heatmaps

### Real data (main)
1. Loci selection: scatter of difficulty metrics (M1 vs $Z_{\text{primary}}$); highlight chosen loci and their failure-map regions.
2. Per-locus comparison: aggregated PIPs and credible sets; multimodality metrics; differences in prioritized SNPs.

### Table (main or supplement)
- Summary across loci/datasets: mean improvements, stability, diversity ratios, compute budget.

---

## 7) Implementation plan (HPC-friendly)
- Treat each dataset \((X,y,a)\) as a node-level unit of work:
  - dataset = (genotype id, phenotype settings, annotation settings, seed)
- Pipeline:
  1. generate dataset
  2. compute difficulty metrics
  3. run all candidate fits (by use case x exploration)
  4. aggregate outputs
  5. compute and save all metrics
- Store:
  - per-fit PIP vectors (and optionally \(\alpha\) matrices)
  - per-fit ELBO (if available)
  - aggregated PIPs for each rule
  - all metrics (difficulty, multimodality, performance)

---

## 8) Scope control (to keep the paper finishable)

**[UPDATED 2025-07-07 per D1/D3/D4/D5 resolution — see `analysis_completion_status.md`]**

1. ~~Do not implement SuSiNE-ash/inf unless absolutely necessary.~~ **D3 RESOLVED:** SuSiE-ash/inf are free baseline arms via susieR 2.0 (`unmappable_effects`). Include them.
2. ~~Keep the "functional SuSiE via π" baseline minimal.~~ **D5 RESOLVED:** Implement functional pi as `softmax(|a_j| / tau)` with tau-grid exploration (D15).
3. Prefer a small, interpretable \(c\)-grid (and small \(\sigma_0^2\)-grid) over a huge Cartesian product. *(unchanged)*
4. ~~Treat SuSiE refinement as an exploration competitor only if it fits cleanly.~~ **D4 RESOLVED:** Refinement is included. susieR 2.0 `refine=TRUE` for SuSiE baselines; harness-level BFS refinement for all models.

---

## 9) Immediate next actions (concrete TODO list)
~~**Create the "SuSiE Pathology" R notebook vignette**~~ — **DONE** (`vignettes/susie_pathology.ipynb`). Introduces AND-of-ORs vs OR-of-ANDs taxonomy and model-specification vs optimizer barrier dichotomy.
2. Confirm Bioinformatics scope/format constraints (figure limits, supplement norms).
3. Finalize simulation harness changes:
   - alpha-based convergence for comparability
   - consistent genotype QC filters (MAF, missingness)
3. Implement flatten-based aggregation pipeline (see `analysis_completion_status.md` D6):
   - save all per-fit PIP vectors + ELBOs
   - post-hoc: flatten all fits per (dataset, model_spec), apply four aggregation methods
4. Implement dataset difficulty metrics (M1 + full z-score metric family).
5. Implement multimodality metrics (JSD, clustering, top-k Jaccard, PIP variance).
6. Run metric-screening analysis to choose one primary z-score metric for main text.
7. Run the compute-confound pilot on 5–10 datasets.
8. Expand real-data study to 3–4 additional loci chosen by difficulty spread on the $(M_1, Z_{\text{primary}})$ plane.
9. Add minimal validation for summary-statistic implementation (or remove claim).

---

## Appendix: glossary (minimal)
- **IBSS:** iterative Bayesian stepwise selection; coordinate ascent over SuSiE effects.
- **PIP:** posterior inclusion probability per SNP.
- **CS:** credible set per effect; intended to contain one causal variant per effect.
- **M1:** LD complexity metric \(\sum_{i<j} |r_{ij}|(1-|r_{ij}|)\).
- **Z-primary metric:** selected z-score-only difficulty metric used with M1 to stratify failure modes.
- **JSD:** Jensen–Shannon divergence between PIP vectors (diversity proxy).
- **Flatten-based aggregation:** pool all fits per (dataset, model_spec); aggregate via Max ELBO, Uniform, ELBO softmax, or Cluster-then-weight (Method A). See `analysis_completion_status.md` D6.




