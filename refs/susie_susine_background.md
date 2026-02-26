# SuSiE / SuSiE-2 / SuSiNE: Background Brief (handoff)

Goal of this note: quickly bring a technically literate reader up to speed on (i) what **vanilla SuSiE** is and why it can fail, (ii) what **SuSiE-ash** and **SuSiE-inf** change, (iii) how existing functional-prior work typically enters through **prior inclusion probabilities** $\pi$, and (iv) what **SuSiNE** adds by allowing **nonzero prior means** for effect sizes (with an unknown scale), plus a few other SuSiE-2 “plumbing” updates that matter for benchmarking.

---

## 1) Vanilla SuSiE

### 1.1 Purpose
**SuSiE (Sum of Single Effects)** is a Bayesian variable selection regression / genetic fine-mapping model designed to infer which variants in a locus are causal (via **posterior inclusion probabilities**, PIPs) while accounting for LD among variants. It provides:
- Variant-level PIPs
- Credible sets (CS) intended to contain one causal variant per effect component
- A fast approximate posterior that scales to typical cis-eQTL loci

In practice, SuSiE is a computationally efficient alternative to exact or near-exact **model-space exploration** methods (e.g., FINEMAP, DAP-G, CAVIAR) that sum/integrate over discrete multi-SNP models but can become expensive as $p$ and the number of causal variants ($p^\*$) increase.

### 1.2 Model specification (core ideas)
**Likelihood (standard linear model):**
$$
y = X b + \varepsilon, \quad \varepsilon \sim \mathcal N(0, \sigma^2 I_n),
$$
where $X\in\mathbb R^{n\times p}$ is a standardized genotype matrix (typically centered/scaled), $b\in\mathbb R^p$ are SNP effect sizes, and $\sigma^2$ is residual variance.

**Sum of single effects decomposition:**
$$
b = \sum_{\ell=1}^L b^{(\ell)},
$$
with each component $b^{(\ell)}$ constrained to be **“single-effect”** (exactly one nonzero entry), i.e.
$$
b^{(\ell)} = \gamma^{(\ell)} \cdot \theta_{\ell},
$$
where $\gamma^{(\ell)}\in\{e_1,\dots,e_p\}$ is a one-hot indicator selecting a SNP, and $\theta_{\ell}$ is that effect’s scalar magnitude. The prior typically takes:
- $\Pr(\gamma^{(\ell)} = e_j) = \pi_j$ (often uniform $\pi_j=1/p$)
- $\theta_{\ell} \sim \mathcal N(0, \sigma_0^2)$

So each effect $\ell$ chooses one SNP (via $\gamma^{(\ell)}$) and gives it a Normal prior on its effect size. The posterior over $\gamma^{(\ell)}$ yields per-effect “alpha” vectors $\alpha^{(\ell)}$ (probabilities over SNPs), and the final **PIP** for SNP $j$ is:
$$
\text{PIP}_j = 1 - \prod_{\ell=1}^L (1 - \alpha^{(\ell)}_j).
$$

### 1.3 Algorithm (IBSS / coordinate ascent variational inference)
SuSiE uses a variational approximation in which each single-effect component has its own variational factors (intuitively: each effect is updated while holding the others fixed). The classic fitting loop is often described as:

1. Initialize $\{\alpha^{(\ell)}\}_{\ell=1}^L$, $\{\mu^{(\ell)}, v^{(\ell)}\}$ (effect-size posterior mean/var for each SNP within each component), and hyperparameters ($\sigma^2$, $\sigma_0^2$, etc.).
2. For $\ell = 1,\dots,L$ repeatedly:
   - Form a residual by subtracting contributions of all other effects (a “leave-one-effect-out” residual).
   - Fit a **single-effect regression** to that residual, producing updated $\alpha^{(\ell)}$ and effect-size posterior parameters.
3. Update hyperparameters (often by empirical Bayes steps) and iterate until convergence.

This is a form of **coordinate ascent** on a variational objective (ELBO). It is fast and scales well.

### 1.4 Where vanilla SuSiE fails (why multimodality matters)

#### (A) Local optima from greedy / coordinate updates
IBSS is a coordinate-ascent procedure. Like many greedy/coordinate algorithms, it can get stuck in **local optima** when:
- The posterior is multi-modal over discrete SNP configurations.
- Different configurations explain the data similarly well.

Two intuitive constructions:

**1) “Pair-only signal” (marginals negligible, joint strong):**  
Construct two causal SNPs whose *joint* inclusion explains substantial variation, but each marginal association is small due to cancellation / correlation patterns. A greedy single-effect update that relies heavily on marginal regression statistics can fail to pick either early, and subsequent components may fill with suboptimal proxies.

**2) “Bandit variables” (sequential confounding):**  
A noncausal SNP $x_b$ highly correlated with a true causal SNP $x_c$ can look best marginally. Once selected, the remaining conditional signal can be weak or misassigned, leading later components to pick a second “bandit” that jointly mimics a causal pair. The sum-of-single-effects structure is robust to some *single* proxy/bandit situations (it can place probability mass across LD proxies within a component), but it can struggle with **bandit pairs / triplets** where a *set* of noncausal SNPs together approximates the signal of the true causal set, yet no single SNP is strongly correlated with the true set as a whole.

In contrast, discrete model-space samplers (FINEMAP / DAP-G / CAVIAR) explicitly enumerate or stochastically traverse multi-SNP models, allowing them (at higher computational cost) to compare evidence across **distinct multi-SNP configurations** more directly.

#### (B) Multi-modality: “one VI solution ≠ multiple basins”
Even when SuSiE’s posterior approximation allows uncertainty *within* a component (via $\alpha^{(\ell)}$), the overall fitted solution is still a **single** variational optimum. If there are multiple well-separated basins (different sets of CS allocations / effect decompositions), a single optimum may reflect only one basin, depending on initialization, hyperparameters, and update order.

This motivates **ensembling / multi-start / hyperparameter grids**:
- Not because “SuSiE is bad,” but because the optimization landscape can contain multiple plausible explanations that a single VI optimum cannot represent simultaneously.

#### (B′) Structural taxonomy: AND-of-ORs vs OR-of-ANDs

SuSiE's variational posterior is a product of $L$ independent categorical distributions (one per effect component). This means a single fit naturally represents **AND-of-ORs** uncertainty — each component independently spreads mass across alternatives (e.g., a causal SNP and its LD proxy). When the true uncertainty has this structure, a single fit suffices.

What a single fit *cannot* represent is **OR-of-ANDs** uncertainty — “either this complete multi-SNP configuration or that one” — because this requires *correlated* choices across components. When two sets of variables span the same column subspace (even with moderate pairwise correlations), each fit commits to one set, creating distinct posterior basins.

This distinction maps to intervention type:

- **Model-specification barrier** (OR-of-ANDs): ensembling across fits is necessary to represent cross-basin uncertainty.
- **Optimizer barrier** (greedy trap): a single better basin exists; multi-start selection or `refine` suffices.

See `vignettes/susie_pathology.ipynb` (Scenarios 2–4) for concrete demonstrations.

#### (C) Mis-specification / “background effects” (hallucinated sparse effects)
Vanilla SuSiE assumes the signal is well represented by a **small number** $L$ of single-SNP effects plus noise. When the true architecture includes many weak effects (polygenic background), SuSiE can:
- Attribute diffuse background signal to a smaller number of sparse components,
- Produce **synthetic associations** or overconfident CSs,
- Inflate CS-level false discoveries under certain complex architectures.

This motivates SuSiE-2 variants that explicitly model “background” differently (next section).

---

## 2) Existing SuSiE extensions: SuSiE-ash and SuSiE-inf

SuSiE-2 introduces alternative prior/likelihood modeling choices aimed at robustness when the “few sparse effects” assumption is strained.

### 2.1 SuSiE-ash (adaptive shrinkage background component)
**Idea:** Keep the sparse SuSiE component for mappable effects, but add an extra **background/unmappable component** with an adaptive-shrinkage (ash) **mixture-of-normals** prior. In the SuSiE-2 formulation this is:
$$
y = X\beta + X\theta + \varepsilon,
$$
where $\beta$ is the sparse SuSiE part and $\theta$ captures diffuse moderate/weak effects via ash. The ash prior can represent:
- A spike near zero (many tiny effects),
- A range of small/moderate/large effects,
- Better-calibrated shrinkage across effect magnitudes.

**What it addresses:**
- Better calibration / FDR control in settings with many weak/moderate effects or heavy-tailed effect distributions.
- Reduces overconfidence that can arise when the prior is too rigid.

**What it does *not* inherently fix:**
- Multi-modality and local optima due to LD and combinatorial configurations can still exist.
- It changes the “shape” of the posterior in each basin but does not guarantee exploration across basins.

**Practical note for benchmarking:** SuSiE-ash is best viewed as changing **background-effect modeling and shrinkage calibration** rather than adding a fundamentally different search/exploration algorithm.

#### 2.1.1 Rough implementation spec (SuSiE-ash)

> **Note (2025-07-07):** SuSiE-ash is now fully implemented in upstream susieR 2.0 via `susie(X, y, unmappable_effects = "ash")`. The spec below remains useful as a conceptual reference but is no longer a porting target. See `refs/susieR_2.0_inventory.md` for the actual implementation details (C++ backend via `caisa_rcpp`, masking/unmasking pattern, individual-data only constraint).

If someone needs a first working implementation, this is the minimum useful spec.

**Model (individual-level):**
$$
y = X\beta + X\theta + \varepsilon,\qquad \varepsilon\sim\mathcal N(0,\sigma^2 I_n)
$$
Sparse part (SuSiE):
$$
\beta=\sum_{\ell=1}^L \beta^{(\ell)},\quad \beta^{(\ell)}=b_\ell\gamma^{(\ell)},\quad
\gamma^{(\ell)}\sim\text{Mult}(1,\pi),\quad b_\ell\sim\mathcal N(0,\sigma_{0\ell}^2).
$$
Background part (ash):
$$
\theta_j\sim\sum_{k=1}^K w_k\,\mathcal N\!\left(0,\sigma^2 s_k^2\right),
$$
where $s_k^2$ is a fixed nonnegative variance grid (include $s_1^2=0$), and $w_k$ are mixture weights.

**Key quantities:**
- Global background variance:
$$
\tau^2=\mathrm{Var}(\theta_j)=\sigma^2\sum_{k=1}^K w_k s_k^2.
$$
- Precision used by sparse-effect updates:
$$
\Omega=(\sigma^2 I_n+\tau^2 XX^\top)^{-1}.
$$

**One outer iteration (GIBSS-style):**
1. Update sparse effects $\{\beta^{(\ell)}\}_{\ell=1}^L$ by IBSS using weighted SER likelihood with $\Omega$.
2. Form posterior mean $\bar\beta=\sum_\ell \alpha^{(\ell)}\odot m_1^{(\ell)}$ and residual $r_\theta=y-X\bar\beta$.
3. Run Mr.ASH-style update on $r_\theta$: for each SNP $j$, update mixture responsibilities $\phi_{jk}$, component posterior means/variances $(\mu_{jk},s_{jk}^2)$, then set $\bar\theta_j=\sum_k\phi_{jk}\mu_{jk}$ and $w_k\leftarrow p^{-1}\sum_j\phi_{jk}$.
4. Update $\sigma^2$ and $\tau^2$ (MoM or EB step), then recompute $\Omega$.
5. Check convergence (`elbo` or max change in $\alpha$), repeat.

**Weighted SER update (per effect $\ell$, per SNP $j$):**
$$
v_{\ell j}=\left(x_j^\top\Omega x_j + \frac{1}{\sigma_{0\ell}^2}\right)^{-1},\qquad
m_{\ell j}=v_{\ell j}\,x_j^\top\Omega r_\ell.
$$
Then compute a log-BF-like score
$$
\log \mathrm{BF}_{\ell j}=\frac{1}{2}\log\!\left(\frac{v_{\ell j}}{\sigma_{0\ell}^2}\right)+\frac{m_{\ell j}^2}{2v_{\ell j}},
$$
and normalize
$$
\alpha_{\ell j}\propto \pi_j \exp(\log \mathrm{BF}_{\ell j}).
$$

**Mr.ASH normal-means update (rough):**
Set $\tilde\theta_j=(x_j^\top x_j)^{-1}x_j^\top r_{\theta,j}$ and $\gamma_j^2=(x_j^\top x_j)^{-1}$. For each mixture component $k$:
$$
\mu_{jk}=\frac{s_k^2}{s_k^2+\gamma_j^2}\tilde\theta_j,\qquad
\phi_{jk}\propto w_k\,\mathcal N\!\left(\tilde\theta_j;0,\sigma^2(s_k^2+\gamma_j^2)\right).
$$
Normalize $\phi_{jk}$ over $k$, then update $w_k\leftarrow p^{-1}\sum_j\phi_{jk}$.

**Minimum pseudo-code sketch:**
```text
init alpha, m1, sigma2, w, variance grid s2
repeat:
  tau2 = sigma2 * sum_k w[k] * s2[k]
  Omega = inv(sigma2 * I + tau2 * X X')

  # sparse update
  for l in 1..L:
    r_l = y - X * sum_{l'!=l} beta_bar[l']
    (alpha[l], m1[l], v1[l]) = SER_weighted(r_l, X, Omega, pi, sigma0_l2)
    beta_bar[l] = alpha[l] .* m1[l]
  beta_bar_total = sum_l beta_bar[l]

  # background ash update
  r_theta = y - X * beta_bar_total
  (theta_bar, w, sigma2) = mr_ash_update(r_theta, X, s2, w, sigma2)

until converged
```

Implementation note: numerically, precomputing an eigendecomposition of $XX^\top$ is useful so repeated $\Omega$ updates are cheap.
Summary-stat note: in an RSS implementation, replace raw $X$ operations with LD/sufficient-stat equivalents (e.g., using $R$, $z$, and sample size), but keep the same outer-loop logic.

### 2.2 SuSiE-inf (infinitesimal/background effect model)
**Idea:** Add an **infinitesimal** (very polygenic) component capturing diffuse background signal, so that the sparse components focus on identifiable “spikes” rather than trying to explain background with sparse effects.

**What it addresses:**
- Settings where there is a genuine polygenic background that otherwise causes SuSiE to “hallucinate” sparse effects.
- Can improve calibration when background effects are the dominant mismatch.

**What it does *not* inherently fix:**
- If sparse signals are still multi-modal due to LD, you can still have multiple basins.
- Depending on implementation details, it may be less competitive on some sparse metrics if the model is too “diffuse” relative to truth.

#### 2.2.1 Rough implementation spec (SuSiE-inf)

> **Note (2025-07-07):** SuSiE-inf is now fully implemented in upstream susieR 2.0 via `susie(X, y, unmappable_effects = "inf")` (individual data) and `susie_ss(unmappable_effects = "inf")` (summary stats). Uses omega-weighted SER updates and BLUP theta. See `refs/susieR_2.0_inventory.md` for details. Incompatible with RSS (λ>0).

SuSiE-inf can be viewed as the same decomposition, but with a single Gaussian background prior:
$$
y=X\beta+X\theta+\varepsilon,\qquad \theta_j\sim\mathcal N(0,\tau^2),\qquad \varepsilon\sim\mathcal N(0,\sigma^2I_n).
$$
So the sparse updates again use:
$$
\Omega=(\sigma^2 I_n+\tau^2 XX^\top)^{-1}.
$$

**Practical rough algorithm:**
1. Initialize $\sigma^2,\tau^2$ and standard SuSiE parameters.
2. Given current $(\sigma^2,\tau^2)$, run one IBSS sweep for sparse effects under weighted SER likelihood with $\Omega$.
3. From current posterior moments of $\beta$, update $(\sigma^2,\tau^2)$ (MoM or EB update).
4. Recompute $\Omega$, iterate until convergence.

A common rough MoM update solves:
$$
\mathbb E_q\!\left[\lVert y-X\beta\rVert^2\right]=\sigma^2 n+\tau^2\operatorname{tr}(X^\top X),
$$
$$
\mathbb E_q\!\left[\lVert X^\top(y-X\beta)\rVert^2\right]
=\sigma^2\operatorname{tr}(X^\top X)+\tau^2\operatorname{tr}\!\big((X^\top X)^2\big),
$$
for $(\sigma^2,\tau^2)$ each outer iteration.

In code, treat this as a 2x2 linear system:
$$
\begin{bmatrix}
n & \operatorname{tr}(X^\top X)\\
\operatorname{tr}(X^\top X) & \operatorname{tr}((X^\top X)^2)
\end{bmatrix}
\begin{bmatrix}
\sigma^2\\
\tau^2
\end{bmatrix}
=
\begin{bmatrix}
\mathbb E_q[\lVert y-X\beta\rVert^2]\\
\mathbb E_q[\lVert X^\top(y-X\beta)\rVert^2]
\end{bmatrix},
$$
then clip to small positive floors for numerical stability.

**Minimum pseudo-code sketch:**
```text
init alpha, m1, sigma2, tau2
repeat:
  Omega = inv(sigma2 * I + tau2 * X X')
  for l in 1..L:
    r_l = y - X * sum_{l'!=l} beta_bar[l']
    (alpha[l], m1[l], v1[l]) = SER_weighted(r_l, X, Omega, pi, sigma0_l2)
    beta_bar[l] = alpha[l] .* m1[l]
  (sigma2, tau2) = update_sigma2_tau2_from_beta_moments()
until converged
```

In practice, SuSiE-inf and SuSiE-ash share most infrastructure; SuSiE-ash replaces the single $\tau^2$ background prior with a learned ash mixture that still feeds into a global $\tau^2$ for the sparse-effect precision update.
Summary-stat note: the same applies in `susie_rss`-style code; only the linear algebra backend changes, not the update ordering.

### 2.3 How to interpret these variants relative to your project
- **ash/inf** primarily address *model mis-specification* (effect-size distribution and background), not the **unknown-scale directional prior mean** problem.
- Your method (SuSiNE) is designed to incorporate **orthogonal information** (signed prior means) and to encourage **better exploration/selection of modes** under scale uncertainty.

---

## 3) SuSiE prior extensions in existing functional fine-mapping pipelines ($\pi$-only)

Many “functional prior” approaches modify SuSiE (or other fine-mappers) by changing the **prior inclusion probabilities** $\pi_j$, i.e.
$$
\Pr(\gamma^{(\ell)} = e_j) = \pi_j \quad \text{instead of uniform}.
$$

### 3.1 PolyFun-style functional priors (conceptual)
A common approach is:
- Use functional annotations (chromatin, conservation, etc.) and external GWAS/eQTL training to estimate per-SNP “prior causality” scores.
- Convert those scores into $\pi_j$ (often via logistic/softmax transformations, regularization, and LD-aware adjustments).
- Run a fine-mapper using these $\pi_j$ as prior inclusion weights.

**Key limitation for your context:** $\pi$-only priors encode **“which SNP is more likely causal”**, but do not directly encode:
- direction of effect (sign),
- expected effect size magnitude in trait units,
- or a directional “logic gate” that can downweight basins inconsistent with functional directionality.

### 3.2 Sequence-to-function model features used as $\pi$-priors
In some GWAS contexts, sequence-to-function models (e.g., predicting regulatory activity changes) are used to derive broad functional scores that are then mapped into $\pi_j$. This can involve additional steps to:
- aggregate across tissues/cell types,
- map predicted regulatory effects to trait relevance,
- normalize/transform into a usable prior inclusion weight.

**Takeaway:** these approaches typically remain $\pi$-only in the fine-mapping model, because $\pi$ is the “standard injection point” for functional information in many fine-mappers.

---

## 4) SuSiNE prior extension: nonzero prior means (directional priors)

### 4.1 What SuSiNE adds
SuSiNE generalizes the per-effect prior on effect size from:
$$
\theta_{\ell} \sim \mathcal N(0,\sigma_0^2)
$$
to a **noncentral Normal**:
$$
\theta_{\ell} \sim \mathcal N(\mu_{0j},\sigma_0^2) \quad \text{(when SNP } j \text{ is selected)}.
$$

In practice, you construct $\mu_{0}\in\mathbb R^p$ from functional annotations that are **directly interpretable as signed effect-direction predictions**, and then introduce an **unknown global scale** $c$ such that:
$$
\mu_{0} = c\, a,
$$
where $a$ is a signed annotation vector (from a sequence-to-function model) and $c$ converts annotation units to phenotype-effect units.

### 4.2 Why nonzero means are “orthogonal” to $\pi$-only priors
- $\pi$-only priors can say “SNP j is more likely causal,” but not “if causal, the effect is likely positive and of roughly this magnitude.”
- A nonzero mean prior acts as a **directional filter / logic gate**:
  - If a basin requires effects opposite to the functional direction, it becomes less plausible under the prior.
  - If a basin aligns with functional direction, it receives a boost.
- This can meaningfully reshape the posterior landscape, not just reweight within a basin.

### 4.3 The key practical complication: scale $c$ is unknown and often non-identifiable
In eQTL contexts, sequence-to-function predictions can be **directionally meaningful**, but the mapping to effect-size units is usually not known a priori. As a result:
- A single EB estimate of $c$ (or $\sigma_0^2$) can be unstable / circular (since estimating $c$ depends on which SNPs are causal, which is what you are trying to infer).
- A **grid over $c$** (and/or $\sigma_0^2$) paired with **ensembling** becomes a principled workflow:
  - explore multiple basins induced by different $c$,
  - then aggregate models under a fixed compute budget.

This is central to the story: the “knob” $c$ is both a prior parameter and a practical exploration axis.

---

## 5) Other SuSiE-2 updates relevant for benchmarking (non-exhaustive)

These are not the core scientific novelty for your project, but they matter for reproducing/benchmarking SuSiE-2 behavior.

### 5.1 Refinement (model branching / exploration)
SuSiE-2 includes an enhanced **refinement** stage (enabled by `refine = TRUE` in `susie`/`susie_ss`) that explicitly explores alternative basins after an initial IBSS fit. The procedure used in `susieR` is:
1. Fit an initial model $s$ and obtain its credible sets (CSs). Suppose there are $K$ CSs.
2. For each CS $k$, set prior weights of SNPs in that CS to zero and refit, producing candidate model $t_k$.
3. For each $t_k$, run another fit initialized from $t_k$ (`\alpha,\mu,\mu^2`) to get $s_k$.
4. If $\max_k \text{ELBO}(s_k) > \text{ELBO}(s)$, replace $s$ with the best $s_k$ and repeat; otherwise stop.

This creates targeted branches away from currently selected CS configurations, and compute cost generally grows with the number of CSs explored per refinement round.

### 5.2 Convergence criterion: flexible ELBO- or PIP-based stopping
SuSiE-2 exposes flexible convergence rules rather than a hard switch away from ELBO. In `susieR`, `convergence_method` supports:
- `"elbo"`: stop when ELBO improvement is below tolerance (default behavior).
- `"pip"`: stop when the maximum absolute change in posterior allocations $\alpha$ is below tolerance.

### 5.3 Regularized LD / numerical stabilization
SuSiE-2 emphasizes numerically stable handling of LD (e.g., regularization) to improve robustness when working from summary LD matrices or noisy estimates.

---

## 6) Summary: how the pieces fit your project narrative

**Vanilla SuSiE**: fast VI fine-mapper; can fail due to (i) local optima/multimodality and (ii) mis-specification under background effects.  
**SuSiE-ash / SuSiE-inf**: add explicit background-effect modeling (ash mixture or infinitesimal), improving calibration under complex architectures; do not inherently solve unknown-scale directional prior means.  
**Functional prior work (PolyFun, etc.)**: usually injects functional information via $\pi$ only (inclusion probabilities).  
**SuSiNE**: injects *signed* functional information through **nonzero prior means**, introducing an unknown scale $c$ that (a) must be handled robustly and (b) can be exploited as an exploration axis; combined with ensembling, this targets multimodality + scale uncertainty in eQTL fine-mapping.

---

## 7) Open questions for the team (optional, to align scope)
These questions are here so a new collaborator/AI can quickly diagnose what remains ambiguous.

1. **Primary output**: are we optimizing for (a) calibrated posterior PIPs/CSs, or (b) robust candidate selection under hyperparameter uncertainty?
2. **Compute budget definition**: will we compare methods under a fixed number of model fits per locus, or fixed walltime? **→ D8 RESOLVED: fixed K (number of fits); record wall time per fit.**
3. **Annotation error model**: will we include explicit sign-flip scenarios, LD-correlated annotation noise, or just R²-controlled Gaussian noise? **→ D14 RESOLVED: current controls sufficient, no sign flips or LD-correlated noise.**
4. **Role of SuSiE-ash/inf**: do we treat them as baseline "engines" (single-fit) or include them in any ensemble comparisons (even via non-ELBO pooling)? **→ D3 RESOLVED: include as baseline comparison arms via susieR 2.0 (`unmappable_effects`). Zero implementation cost. See `analysis_completion_status.md` D3.**


