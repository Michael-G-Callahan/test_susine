# Future Ideas & Intellectual Parking Lot

Ideas that are interesting but out of scope for the current SuSiNE paper. Organized loosely by theme.

---

## 1. IBSS as Mirror Descent / Particle-Based Variational Inference

IBSS coordinate ascent can be viewed through the lens of **mirror descent** on the simplex (the alpha updates are multiplicative / entropic). This opens several directions:

- **Formal connection:** Cast each SER update as a mirror descent step with KL divergence as the Bregman divergence. Characterize convergence rates and basin geometry under this lens.
- **Particle-based flow:** Treat each model fit (PIP vector) as a **particle** on the probability simplex and define a repulsive-attractive flow analogous to **Stein Variational Gradient Descent (SVGD)**. Particles would be attracted toward high-ELBO regions but repelled from each other via a kernel (JSD-based, naturally). This could give a principled multi-modal variational approximation without post-hoc ensembling.
- **Connection to Wasserstein gradient flows:** The simplex geometry + JSD metric already gives us a natural Riemannian structure (Fisher-Rao). A particle flow on this manifold could be related to Wasserstein gradient descent on the space of posteriors.
- **Practical implication:** Instead of "run K independent fits, then aggregate," you'd run K *interacting* fits that jointly explore the posterior landscape. The repulsive kernel prevents mode collapse; the attractive term keeps particles near high-density regions.

**Why it's out of scope:** This is a fundamentally different algorithm, not an analysis of the current method. It's a follow-up paper, not a paragraph in this one.

---

## 2. Method B: KDE Importance Sampling for Large Ensembles

From `refs/pip_ensemble_methods.md`. For K = 100+ restarts, replace discrete clustering (Method A) with continuous density estimation on the PIP simplex:

- Use the **heat kernel** on the simplex (the geometrically principled kernel under Fisher-Rao / JSD):
  $$K_t(\text{PIP}_i, \text{PIP}_j) \propto \exp(-\text{JSD}(\text{PIP}_i/L, \text{PIP}_j/L) / 2t)$$
- Bandwidth $t$ selected by leave-one-out cross-validation on the JSD matrix alone (no external data needed).
- Importance weights: $w_i \propto \exp(\text{ELBO}_i) / \hat{f}(\text{PIP}_i)$ where $\hat{f}$ is the KDE estimate of optimizer discovery frequency.
- The information-geometric foundation (heat kernel as fundamental solution to diffusion on the simplex, Varadhan's approximation) is elegant and publishable on its own.

**Why it's out of scope:** We won't run K = 100+ restarts in the current study (budget is ~5-20 fits per model spec). Method A (cluster-then-weight) covers the small-K regime cleanly. Method B becomes relevant if someone scales this up or applies it to a setting where restarts are cheap.

---

## 3. SuSiE-ash / SuSiE-inf Under Ensemble Framework

SuSiE-ash (adaptive shrinkage background) and SuSiE-inf (infinitesimal background) address model mis-specification from polygenic effects rather than multimodality. But they could interact with ensembling in interesting ways:

- **Do background components reduce multimodality?** If the sparse components no longer need to "hallucinate" diffuse signal, basins might be better separated. An empirical test: run SuSiE-ash with restarts and measure JSD diversity vs vanilla SuSiE with restarts.
- **SuSiNE-ash:** Combine directional priors with adaptive background. The sparse effects get $\mu_0 = ca$; the background gets an ash mixture. This could be the "ultimate" model but is a lot of machinery.
- **Ensemble across model classes:** Pool SuSiE, SuSiE-ash, SuSiE-inf, and SuSiNE fits and let the aggregation sort it out. ELBO comparability across different likelihoods is questionable, but Method A's JSD-based clustering might sidestep this (cluster by PIP similarity, not ELBO).

**Why it's out of scope:** High implementation cost (weighted SER with Omega, iterative background estimation). The paper's message is about exploration/aggregation under prior uncertainty, not background modeling. Per study plan Section 8: "do not implement unless absolutely necessary."

---

## 4. Full-Grid BMA (Cross-Prior ELBO Re-Scoring)

In our D6 discussion, we considered evaluating every fit's ELBO under every prior specification (not just its own). The "full grid" approach would give a $K \times G$ matrix of ELBOs (K fits, G prior specs), then marginalize:

$$w_k \propto \sum_g \exp(\text{ELBO}(q_k; \text{prior}_g)) \cdot p(\text{prior}_g)$$

**Why it's intellectually interesting:**
- This is the "proper" BMA if we treat prior specs as a discrete hyperprior.
- It would let a fit discovered under one prior spec be recognized as good under a different prior spec.
- Connection to thermodynamic integration / path sampling: the grid of priors defines a path through model space.

**Why we killed it:** The variational posterior $q_k$ is optimized for prior spec $k$. Evaluating it under a different prior gives a *loose* ELBO bound — the gap between the bound and the true log marginal likelihood is uncontrolled and potentially large. This makes the off-diagonal entries unreliable, and the resulting weights could be worse than diagonal-only. Would need to re-optimize $q$ under each prior (K x G fits instead of K fits) to do this properly.

**Possible future rehabilitation:** If someone develops a fast "prior-perturbation" correction (e.g., first-order Taylor expansion of the ELBO around the fitted prior), the off-diagonal entries might become usable without full re-optimization.

---

## 5. Annealing / Tempering as Exploration

Dropped from the paper (D12 RESOLVED) but the susine package retains the implementation. Ideas:

- **Simulated annealing on the ELBO:** Temperature schedule $T(t)$ that flattens the ELBO landscape early (encouraging exploration) and sharpens it later (convergence to a mode). Already implemented in susine.
- **Parallel tempering:** Run multiple chains at different temperatures simultaneously, with swap moves. This is the MCMC version of what our grid approach does approximately.
- **Annealing as "soft refinement":** Instead of the hard basin-hopping of `refine = TRUE` (zero out CSs and refit), anneal gradually to let the optimizer drift between basins.

**Why it's out of scope:** Adds another exploration axis to an already complex study. The paper focuses on grid-based exploration (sigma_0^2, c, tau) + restarts, which are simpler to reason about and compare fairly. Annealing's compute cost is harder to account for in the "fixed budget K" framework.

---

## 6. EB Estimation of c (The Road Not Taken)

We treat c as a non-identifiable exploration knob. But what if you could estimate it?

- **Why EB fails:** Estimating c requires knowing which SNPs are causal (to evaluate whether their effects align with $ca$), but which SNPs are causal is what you're trying to infer. Circular.
- **Possible escape hatches:**
  - Two-stage: use a quick-and-dirty method (e.g., marginal regression + thresholding) to get a rough causal set, estimate c from that, then run SuSiNE. Fragile but might work in easy cases.
  - Hierarchical model: put a prior on c and integrate it out analytically (conjugate if the prior is Gaussian and the model is linear). The resulting marginal likelihood would be a function of $a^T \hat{\beta}$ type quantities.
  - Variational EB: add c as a variational parameter. But this is exactly the circularity problem — the ELBO landscape over c may be flat or multimodal.
- **Connection to transfer learning:** If you have c estimated from a training set of loci, you could use it as a prior for new loci. This is essentially what PolyFun does for pi, but for effect size means.

**Why it's out of scope:** The paper's argument is precisely that c *shouldn't* be point-estimated. Showing that grid + aggregate beats EB is part of the story (Phase A pathology vignette). Exploring EB estimation of c would undermine the framing.

---

## 7. Heavy-Tailed / Mixture Effect-Size Distributions

Current simulations use Normal(0, effect_sd) for causal effect sizes. More realistic alternatives:

- **t-distribution** effects (heavier tails, occasional large effects)
- **Spike-and-slab mixture:** most causal effects small, a few large
- **Sign-constrained effects:** all causal effects positive (or matching annotation sign) — tests whether SuSiNE's directional prior is *necessary* when the truth is already directional

**Why it's out of scope:** Study plan says "optionally include heavy-tail / mixture in sensitivity analyses." The Normal case is standard and sufficient for the main claims. Heavy tails would be a supplementary robustness check at most.

---

## 8. Global Cross-Model-Spec Ensembles (D7)

Pool ALL fits from ALL model specs (SuSiE, SuSiE+pi, SuSiNE) for a given dataset and aggregate.

**Why it's interesting:**
- Maximizes diversity — different priors induce genuinely different posterior landscapes.
- Tests whether the aggregation machinery can identify good modes regardless of which model found them.
- A "SuSiE fit that happens to land in a SuSiNE-favorable basin" is still a good fit.

**Why it's problematic:**
- ELBO comparability: different priors produce different ELBO scales. A SuSiNE fit at c=0.5 has a different prior contribution to the ELBO than a SuSiE fit at c=0. Softmax on incomparable ELBOs is meaningless.
- Method A might help: cluster by JSD (PIP similarity, not ELBO), then weight within clusters. But the importance correction still needs comparable ELBOs.
- Possible fix: normalize ELBOs by subtracting the prior's contribution, comparing only likelihood terms. But the variational approximation quality also differs across priors.

**Status:** May appear as a "pooled modes" supplementary analysis per study plan Phase C, but not a primary method.

---

## 9. Annotation Error Models (Sign Flips, LD-Correlated Noise)

Current annotation simulation uses R^2-controlled Gaussian noise. More adversarial scenarios:

- **Sign flips:** A fraction of SNPs have annotations with the wrong sign. Tests SuSiNE's robustness to partially incorrect directional information.
- **LD-correlated annotation noise:** Annotation errors that track LD structure (e.g., a causal SNP's neighbors all get similar annotation values regardless of their true effect). This is realistic for sequence-based predictors that can't distinguish causal from LD-proxy SNPs.
- **Annotation-truth misalignment patterns:** Cases where the annotation correctly predicts *which* SNPs are important but gets the *magnitude* wrong (e.g., the annotation is thresholded / saturated).

**Why it's deferred:** Study plan says "only if needed." The R^2-controlled noise already spans a range from useless (R^2 = 0) to strong (R^2 = 0.5). Sign flips and LD-correlated noise are refinements that matter more for a methods-comparison paper than for establishing the basic ensemble framework.

---

## 10. Spectral Shrinkage for LD Regularization

For real-data analysis (Phase D), noisy LD matrices from reference panels need regularization. SuSiE-2 uses spectral shrinkage. We currently have only basic symmetrization.

- **Spectral shrinkage:** Eigendecompose R, shrink small eigenvalues toward zero (or a floor), reconstruct. Controls the condition number.
- **Banding:** Zero out entries beyond a certain genomic distance. Exploits the block-diagonal-ish structure of LD.
- **Connection to random matrix theory:** The Marchenko-Pastur distribution tells you which eigenvalues are "signal" vs "noise" in the LD matrix. Could inform an adaptive shrinkage threshold.

**Why it's deferred:** Only needed for real data. Simulations use true genotype matrices (no noise in LD). Will implement when we get to Phase D.

---

## 11. Stacking / LOO Weights as Alternative Aggregation

Study plan B3 mentions "stacking/LOO-style weights if you can produce per-observation predictive densities reliably." The idea:

- For each fit $k$, compute leave-one-out predictive log-likelihood for each observation.
- Find stacking weights $w$ that maximize $\sum_i \log(\sum_k w_k \cdot p_k(y_i | y_{-i}))$.
- These weights are asymptotically optimal for predictive performance (Yao et al., 2018).

**Why it's interesting:** Stacking weights are theoretically superior to BMA weights for prediction, especially under model mis-specification. They don't rely on ELBO as a proxy for marginal likelihood.

**Why it's out of scope:** Requires per-observation predictive densities, which SuSiE's variational approximation doesn't directly provide in a clean form. The computational cost of K x n LOO evaluations is also substantial. And our goal is variable selection (PIP quality), not prediction — stacking optimizes the wrong objective.

---

## 12. Submodularity and Greedy Guarantees

The dataset_metrics.R code computes submodularity bounds. The connection:

- Variable selection in linear regression has a submodular structure (diminishing returns) when variables are not too correlated.
- The greedy algorithm (which IBSS approximates) has a (1 - 1/e) approximation guarantee for submodular maximization.
- But LD breaks submodularity: adding a variable correlated with an already-selected one can have *increasing* returns if they form a cancellation pair.

**Why it's interesting:** Could formalize *when* IBSS is guaranteed to find a good solution (low LD, submodular regime) vs when it's expected to fail (high LD, non-submodular regime). This connects M1 to theoretical guarantees, not just empirical correlations.

**Why it's out of scope:** Heavy theory for a Bioinformatics methods paper. Better suited to a statistics/ML venue.
