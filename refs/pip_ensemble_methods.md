# Ensembling Fine-Mapping PIPs Across Multiple SuSiE Restarts

## Motivation

Running SuSiE from multiple initializations or hyperparameter settings produces multiple PIP vectors, each representing a different local optimum of the variational objective. These optima may correspond to genuinely different causal configurations (modes of the posterior), and the frequency with which the optimizer discovers each mode reflects the **basin of attraction size** under optimization dynamics — not the posterior mass at that mode. A principled ensemble should correct for this discrepancy.

---

## Simple (Diversity-Agnostic) Methods

### Single Best ELBO

Select the PIP vector from the fit with the highest ELBO. Discard all others.

- **Pros:** No tuning, no ambiguity.
- **Cons:** Ignores all information from other modes. If the optimizer gets lucky once on a suboptimal mode, you're stuck with it. Assumes the optimization landscape is unimodal or that the best ELBO found is globally best — both often false.

### Uniform Weighting

$$\text{PIP}_{\text{ens}} = \frac{1}{K} \sum_{i=1}^{K} \text{PIP}_i$$

- **Pros:** Simple, hedges across all fits.
- **Cons:** Treats all fits as equally good regardless of ELBO. If 8 of 10 restarts find the same mediocre mode, that mode dominates — this is exactly the optimizer-frequency bias we want to correct.

### ELBO Softmax Weighting

$$w_i \propto \exp(\text{ELBO}_i / \tau)$$
$$\text{PIP}_{\text{ens}} = \sum_i w_i \cdot \text{PIP}_i$$

where $\tau$ is a temperature parameter.

- At $\tau \to 0$: converges to single-best-ELBO.
- At $\tau \to \infty$: converges to uniform weighting.
- **Pros:** Smooth interpolation, easy to implement.
- **Cons:** Still frequency-biased. If one mode is found 80% of the time, 80% of the softmax-weighted terms point to that mode regardless of $\tau$. The temperature controls how much ELBO matters, not how much diversity matters.

---

## Diversity-Aware Methods

### Core Insight: Importance Sampling Correction

Each optimizer restart is a "sample" from a proposal distribution $q$ over modes, where the discovery frequency reflects basin-of-attraction geometry rather than posterior mass. The target distribution is $p_m \propto \exp(\text{ELBO}_m)$, which approximates the posterior mass at mode $m$. The importance weight is:

$$w_m = \frac{p_m}{q_m} \propto \frac{\exp(\text{ELBO}_m)}{f_m}$$

where $f_m$ is the (estimated) discovery frequency of mode $m$. This upweights high-ELBO modes that are rarely found and downweights frequently-found mediocre modes.

**Key modeling choice:** Weights are computed at the **model level**, not per-SNP. The ensemble PIP is obtained by marginalizing:

$$\text{PIP}_{\text{ens}} = \sum_m w_m \cdot \text{PIP}_m \bigg/ \sum_m w_m$$

This preserves the joint correlation structure across SNPs (e.g., LD-driven mutual exclusivity between causal candidates). Per-SNP mixture weights would destroy this structure and could assign high PIP to configurations that no single model supports.

---

### Method A: Cluster-Then-Weight (Recommended for Small Ensembles, K = 5–20)

**When to use:** Small numbers of restarts where nonparametric density estimation is infeasible.

#### Workflow

1. **Run K restarts.** Collect $\{(\text{PIP}_i, \text{ELBO}_i)\}_{i=1}^K$.

2. **Compute pairwise JSD matrix.** Normalize PIP vectors by $L$ (number of single effects) to place them on the simplex:
$$D_{ij} = \text{JSD}(\text{PIP}_i / L,\; \text{PIP}_j / L)$$

3. **Cluster into modes.** Hierarchical clustering with complete linkage on $D$. Complete linkage defines cluster diameter by the maximum internal pairwise distance, ensuring all fits within a cluster are mutually close. Cut the dendrogram at a natural gap (with K = 5–20, this is typically visually obvious).

4. **Select a representative per cluster.** Take the fit with the highest ELBO within each cluster. Call these $\{(\text{PIP}_m^*, E_m)\}_{m=1}^M$.

5. **Compute importance-corrected weights.** Let $n_m$ = number of fits in cluster $m$, $f_m = n_m / K$. Then:
$$w_m \propto \frac{\exp(E_m)}{f_m}$$
Normalize so $\sum_m w_m = 1$.

6. **Ensemble PIP.**
$$\text{PIP}_{\text{ens}} = \sum_{m=1}^M w_m \cdot \text{PIP}_m^*$$

#### Diagnostics

- **Effective sample size (ESS):** $\text{ESS} = 1 / \sum_m w_m^2$. If ESS $\approx 1$, one cluster dominates. Could indicate genuine posterior concentration or insufficient restarts.
- **Representative sensitivity:** Replace best-ELBO representative with within-cluster average PIP. If ensemble PIP changes meaningfully, clusters are too heterogeneous — consider splitting.
- **K sensitivity:** Run at K and 2K restarts. If weights or ensemble PIP shift substantially, more restarts are needed.

#### Note on ELBO Convergence

The importance correction assumes the best ELBO found per cluster is a good proxy for the log marginal likelihood of that mode. If the optimizer under-converges for some modes (e.g., those with more complex correlation structure), consider running additional optimization iterations on each cluster representative before computing weights.

---

### Method B: KDE-Based Importance Sampling (For Large Ensembles, K = 100+)

**When to use:** Large numbers of restarts where you have enough samples to estimate the proposal density nonparametrically. Avoids the discrete clustering step entirely.

#### Information-Geometric Foundation

PIP vectors normalized by $L$ live on the probability simplex. The natural Riemannian metric on the simplex is the Fisher information metric, under which $\sqrt{\text{JSD}}$ is proportional to the geodesic (Fisher-Rao) distance. The **heat kernel** — the fundamental solution to the diffusion equation on the manifold — is the unique geometrically principled kernel. Under Varadhan's approximation:

$$K_t(\text{PIP}_i, \text{PIP}_j) \propto \exp\!\left(-\frac{\text{JSD}(\text{PIP}_i/L,\; \text{PIP}_j/L)}{2t}\right)$$

This is not an arbitrary choice of Gaussian form — it is the leading-order heat kernel on the simplex. The parameter $t$ is diffusion time (bandwidth).

#### Workflow

1. **Run K restarts.** Collect $\{(\text{PIP}_i, \text{ELBO}_i)\}_{i=1}^K$.

2. **Compute pairwise JSD matrix** $D_{ij}$ as above.

3. **Select bandwidth by leave-one-out cross-validation.** For each candidate $t$, compute:
$$\hat{f}_{-i}(\text{PIP}_i) = \frac{1}{K-1} \sum_{j \neq i} K_t(\text{PIP}_i, \text{PIP}_j)$$
$$\text{CV}(t) = \sum_{i=1}^K \log \hat{f}_{-i}(\text{PIP}_i)$$
Maximize $\text{CV}(t)$ over a grid of $t$ values. No genotype data or external information is needed — this operates entirely on the JSD matrix.

4. **Compute proposal density at each fit.**
$$\hat{f}(\text{PIP}_i) = \frac{1}{K-1} \sum_{j \neq i} K_{\hat{t}}(\text{PIP}_i, \text{PIP}_j)$$

5. **Importance weights.**
$$w_i \propto \frac{\exp(\text{ELBO}_i)}{\hat{f}(\text{PIP}_i)}$$

6. **Ensemble PIP.**
$$\text{PIP}_{\text{ens}} = \sum_i w_i \cdot \text{PIP}_i \bigg/ \sum_i w_i$$

#### Sample Size Considerations

- **Dimensionality is not the bottleneck.** JSD collapses each PIP pair to a scalar, so KDE operates in a one-dimensional distance space regardless of the number of SNPs.
- **K must be large enough for both mode discovery and frequency estimation.** Aim for 5–10+ discoveries per mode. If a mode is found only 1–2 times, its density estimate is unreliable.
- **Stability diagnostic:** Run LOO-CV at multiple values of K (e.g., 50, 100, 200). Plot selected $\hat{t}$ vs. K. Convergence indicates sufficient restarts.
- **ESS diagnostic:** Monitor $\text{ESS} = (\sum w_i)^2 / \sum w_i^2$. Low ESS signals that rare high-ELBO modes are receiving extreme weights — the fundamental tension in importance sampling that no kernel choice resolves.

---

## Summary Table

| Method | Diversity-Aware? | Tuning Parameters | Min K | Complexity |
|---|---|---|---|---|
| Single Best ELBO | No | None | 1 | Trivial |
| Uniform Weighting | No | None | Any | Trivial |
| ELBO Softmax | No | Temperature $\tau$ | Any | Trivial |
| Cluster-Then-Weight | Yes | Dendrogram cut (visual) | 5–20 | Low |
| KDE Importance Sampling | Yes | Bandwidth $t$ (LOO-CV, automatic) | 100+ | Moderate |
