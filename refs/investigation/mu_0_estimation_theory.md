# Why EB Works for sigma_0^2 but Not for c: Theory of Prior Parameter Estimation in SuSiNE

**Date:** 2026-03-05
**Author:** Michael Callahan / Claude (adversarial investigation)
**Status:** Draft for internal validation

---

## 1. Problem Statement

SuSiNE extends SuSiE by placing a nonzero prior mean on effect sizes:

    b_j | gamma_j = 1 ~ N(mu_{0j}, sigma_0^2)

where mu_{0j} = c * a_j for a signed annotation vector a and unknown global scale c.

We have three strategies for learning the prior hyperparameters:

| Strategy | sigma_0^2 | c (or mu_0) | When |
|----------|-----------|-------------|------|
| **Iterative EB** | Brent on log(V) per IBSS iter | Brent on mu_0 or L-BFGS-B on (mu_0, log V) | During fitting |
| **One-shot auto_scale** | N/A | WLS from unconditional marginal | Pre-fitting |
| **Grid ensemble** | Fixed or grid | Fixed grid, aggregate via ELBO softmax | Multiple fits |

Empirically, EB for sigma_0^2 works well (pilot_phase_a shows EB outperforming vanilla). Iterative EB for c was abandoned as unreliable. This document investigates **why**.

---

## 2. The SER Marginal Likelihood

### 2.1 Setup

For a single effect regression (SER) on partial residuals y with design X:

- beta_hat_j = X_j^T y / d_j (univariate OLS coefficient)
- shat2_j = sigma^2 / d_j (OLS variance)
- d_j = ||X_j||^2 (column norm)

The log Bayes factor for SNP j (from `SER.R:98-99`) is:

    log BF_j = log N(beta_hat_j; mu_{0j}, sigma_0^2 + shat2_j)
             - log N(beta_hat_j; 0, shat2_j)

Expanding:

    log BF_j = -1/2 log(1 + sigma_0^2/shat2_j)
             + (beta_hat_j - mu_{0j})^2 / (2*(sigma_0^2 + shat2_j))
             ... wait, let's be careful.

Let V_j = sigma_0^2 + shat2_j. Then:

    log BF_j = -1/2 log(V_j / shat2_j)
             - 1/2 (beta_hat_j - mu_{0j})^2 / V_j
             + 1/2 beta_hat_j^2 / shat2_j

    = -1/2 log(1 + sigma_0^2/shat2_j)
      + 1/2 beta_hat_j^2 / shat2_j
      - 1/2 (beta_hat_j - mu_{0j})^2 / V_j

### 2.2 Posterior inclusion

    alpha_j proportional to pi_j * exp(log BF_j)

The SER marginal log-likelihood (from `update_priors.R:193-211`) is:

    log p(y | sigma_0^2, mu_0, pi) = log sum_j [pi_j * exp(log BF_j)] + const

where the constant is sum(log N(y_i; 0, sigma^2)).

---

## 3. How sigma_0^2 Enters the BF: The "Temperature" Interpretation

### 3.1 Constant-shat2 case

When all SNPs have the same precision (shat2_j = shat2 for all j), which occurs when X has equal column norms:

    log BF_j = -1/2 log(1 + sigma_0^2/shat2)                    [constant across j]
             + 1/2 beta_hat_j^2 / shat2                          [data signal]
             - 1/2 (beta_hat_j - mu_{0j})^2 / (sigma_0^2 + shat2) [prior-data fit]

With mu_0 = 0 (the sigma_0^2-only case):

    log BF_j = -1/2 log(1 + sigma_0^2/shat2) + 1/2 beta_hat_j^2 * sigma_0^2 / (shat2 * V)

The j-dependent part is **beta_hat_j^2 * sigma_0^2 / (shat2 * V)**, which is proportional to beta_hat_j^2 regardless of sigma_0^2.

**Key result:** When shat2 is constant and mu_0 = 0, changing sigma_0^2 scales the log BF contrasts but does not change the SNP ranking. The SNP with the largest |beta_hat_j| always has the highest BF.

sigma_0^2 acts as an **inverse temperature** on the softmax:
- Large sigma_0^2 → small scaling → diffuse alpha (all SNPs get similar probability)
- Small sigma_0^2 → large scaling → concentrated alpha (winner-take-all)

Wait — that's backwards. Let's check. The j-dependent log BF term is:

    beta_hat_j^2 * sigma_0^2 / (2 * shat2 * (sigma_0^2 + shat2))

As sigma_0^2 → infinity: this → beta_hat_j^2 / (2 * shat2), which is the maximum contrast.
As sigma_0^2 → 0: this → 0, so all BFs converge to 0 and alpha → uniform (pi).

So **larger sigma_0^2 = sharper alpha** (more concentrated on best SNP), not diffuser. This matches the intuition that a wider prior accommodates more extreme effects, making the model more confident about which SNP is causal.

### 3.2 Variable-shat2 case

When d_j varies (realistic), the ranking CAN change with sigma_0^2 because:

    log BF_j = -1/2 log(1 + sigma_0^2/shat2_j) + 1/2 beta_hat_j^2 * sigma_0^2 / (shat2_j * V_j)

The log(1 + sigma_0^2/shat2_j) term penalizes SNPs with small shat2_j (well-measured SNPs) more heavily when sigma_0^2 is large. This can swap rankings between a well-measured SNP with moderate signal and a poorly-measured SNP with strong signal. But the effect is typically second-order — the beta_hat_j^2 term dominates.

### 3.3 Summary: sigma_0^2 is a "safe" knob

sigma_0^2 primarily controls **contrast** (sharpness of alpha) without changing **ranking** (which SNP is favored). This means:

1. The EB update for sigma_0^2 has a **unimodal** objective: increase V to match the signal strength, decrease V if no signal
2. The ranking is approximately preserved across iterations → the residuals that the next effect sees don't change much → **limited circularity**
3. The accept/reject step in `update_priors.R:174` almost always accepts because the landscape is smooth

---

## 4. How c Enters the BF: The "Ranking Shift" Interpretation

### 4.1 The c-dependent BF

With mu_{0j} = c * a_j, the j-dependent log BF becomes:

    log BF_j = -1/2 log(V_j / shat2_j)
             + 1/2 beta_hat_j^2 / shat2_j
             - 1/2 (beta_hat_j - c * a_j)^2 / V_j

The c-dependent term is:

    -(beta_hat_j - c * a_j)^2 / (2 * V_j)

This is maximized when **c * a_j = beta_hat_j** — i.e., the prior mean matches the observed signal. Different SNPs want different values of c:

    c_j^* = beta_hat_j / a_j   (the "ideal c" for SNP j)

Unless the annotation perfectly predicts the effects (a_j proportional to beta_hat_j), different SNPs disagree on the optimal c.

### 4.2 c changes the ranking

Unlike sigma_0^2, changing c **directly changes which SNPs are favored**:

- At c = 0: ranking determined by |beta_hat_j| (standard SuSiE)
- At c > 0: SNPs where a_j > 0 and beta_hat_j > 0 get a boost; SNPs where a_j > 0 but beta_hat_j < 0 get penalized
- At c < 0: the reverse

This means alpha(c) is a fundamentally different function than alpha(sigma_0^2):
- alpha(sigma_0^2) preserves ranking (approximately) → smooth, unimodal in sigma_0^2
- alpha(c) can have **discontinuous jumps** as c crosses values where rankings swap → potentially multimodal

### 4.3 Two-SNP example

Consider p = 2 with:
- SNP 1: beta_hat_1 = 2.0, a_1 = 1.0 (annotation says "causal, positive")
- SNP 2: beta_hat_2 = 1.8, a_2 = -0.5 (annotation says "mildly anti-causal")
- shat2 = 1.0, sigma_0^2 = 1.0, V = 2.0

**At c = 0:** log BF_1 = f(4.0), log BF_2 = f(3.24) → SNP 1 wins (barely).

**At c = 1:**
- log BF_1: -(2.0 - 1.0)^2 / 4 = -0.25 penalty reduction (prior helps)
- log BF_2: -(1.8 + 0.5)^2 / 4 = -1.3225 (prior hurts badly)
→ SNP 1 wins decisively.

**At c = -1:**
- log BF_1: -(2.0 + 1.0)^2 / 4 = -2.25 (prior hurts badly)
- log BF_2: -(1.8 - 0.5)^2 / 4 = -0.4225 penalty reduction
→ **SNP 2 wins** (ranking swapped!)

The ranking swap at some critical c* creates a **non-smooth landscape** in the marginal likelihood.

---

## 5. Fisher Information Comparison

### 5.1 Fisher information for sigma_0^2

Under the SER mixture model, the expected information for sigma_0^2 is:

    I(sigma_0^2) = sum_j alpha_j * [ 1/(2*V_j^2) - (beta_hat_j - mu_{0j})^2 / (2*V_j^3) + 1/(2*V_j^2) ]

Simplifying at the MLE:

    I(sigma_0^2) = sum_j alpha_j / (2 * V_j^2)

This is **always positive** as long as at least one alpha_j > 0 and V_j < infinity.

### 5.2 Fisher information for c

With mu_{0j} = c * a_j, by the chain rule:

    d(log BF_j)/dc = a_j * (beta_hat_j - c * a_j) / V_j

The score function is:

    S(c) = sum_j alpha_j * a_j * (beta_hat_j - c * a_j) / V_j

The Fisher information:

    I(c) = sum_j alpha_j * a_j^2 / V_j

**This is always non-negative**, but its magnitude depends critically on **sum_j alpha_j * a_j^2**.

### 5.3 When is I(c) large vs. small?

I(c) = sum_j alpha_j * a_j^2 / V_j

This is large when:
- The SNPs receiving high posterior weight (large alpha_j) also have large |a_j| (informative annotation at the causal SNPs)
- V_j is small (precise measurements at annotated SNPs)

This is small when:
- The annotation is noisy (a_j uncorrelated with the true causal pattern)
- The causal SNPs have small |a_j|
- alpha is diffuse (no strong signal) → all alpha_j are small

**Contrast with I(sigma_0^2):** I(sigma_0^2) = sum_j alpha_j / (2*V_j^2) is always substantial when alpha is concentrated, regardless of annotation quality. sigma_0^2 is identifiable from the data alone; c is identifiable only when the annotation is informative.

### 5.4 Effective sample size interpretation

Think of I(c) as the "effective sample size" for estimating c:
- I(sigma_0^2) ~ n_eff / V^2 where n_eff = 1/sum(alpha_j^2) (inverse Herfindahl of alpha)
- I(c) ~ n_eff * <a^2>_alpha / V where <a^2>_alpha = sum(alpha_j * a_j^2) is the annotation-weighted mean squared annotation

When the annotation is uninformative (a_j independent of the causal structure), <a^2>_alpha is just the population variance of a, which is the same regardless of which SNPs are causal. There's no signal for c.

---

## 6. The Circularity Problem

### 6.1 The IBSS fixed-point equation

The IBSS loop iterates:
1. For each effect l = 1, ..., L:
   a. Compute partial residuals: y_resid = y - sum_{l' != l} X b_hat_{l'}
   b. Run SER on y_resid → get alpha_l, b_hat_l
   c. (Optional) EB update of sigma_0^2 and/or c

### 6.2 sigma_0^2 circularity is benign

When we update sigma_0^2 for effect l:
1. The EB update finds sigma_0^2 that maximizes log p(y_resid | sigma_0^2, mu_0, pi)
2. This changes the **sharpness** of alpha_l but (approximately) not the **ranking**
3. The new b_hat_l = alpha_l * mu_1_l changes in magnitude but not in which SNPs dominate
4. Other effects' partial residuals change slightly → their alpha changes slightly
5. The perturbation is **contractive**: small changes in sigma_0^2 → small changes in b_hat → small changes in other effects' residuals

Formally, in the Gaussian compound decision framework (e.g., Robbins 1956), the NPMLE for the prior variance of a normal means problem is consistent and has a unique maximum. The IBSS decomposition doesn't break this because each SER's marginal likelihood in sigma_0^2 is unimodal.

### 6.3 c circularity is dangerous

When we update c for effect l:
1. The EB update finds c that maximizes log p(y_resid | sigma_0^2, c, a, pi)
2. This changes the **ranking** of SNPs — different SNPs get more/less alpha weight
3. The new b_hat_l can shift dramatically (different SNP selected)
4. Other effects' partial residuals change substantially → their alpha changes substantially
5. The perturbation can be **expansive**: a small change in c → large change in b_hat → large change in other effects' residuals → large change in their c estimates

**Positive feedback loop:** If c is slightly too large, the annotation-aligned SNPs get too much alpha, the model "explains" more variance through them, the partial residuals for other effects shrink, and the next round's c estimate remains inflated because the model has committed to the annotation-driven solution.

### 6.4 The IBSS-EB fixed-point: uniqueness

For sigma_0^2: the fixed-point sigma_0^2*(alpha(sigma_0^2)) is typically unique because alpha(sigma_0^2) varies smoothly and the mapping is contractive.

For c: the fixed-point c*(alpha(c)) may not be unique because alpha(c) can have discontinuities (ranking swaps) and the mapping can be expansive. Multiple fixed points correspond to:
- c = 0: the model ignores the annotation (vanilla SuSiE solution)
- c = c_true: the model correctly leverages the annotation
- c = c_spurious: the model incorrectly leverages the annotation (e.g., wrong sign)

The accept/reject safeguard in `update_priors.R:174` prevents moves to worse solutions, but can get stuck at c = 0 (a local minimum when the annotation is informative but the model hasn't yet discovered this).

---

## 7. The One-Shot Estimator: Why It Works Better

### 7.1 The unconditional marginal model

`estimate_mu_0_scale_factor()` (`update_priors.R:307-362`) uses:

    beta_hat_j ~ N(c * a_j, shat2_j + tau^2)

This is a **single-effect unconditional model** — it doesn't condition on which SNP is causal or on any IBSS state. It treats all p univariate regression coefficients as exchangeable observations from a linear regression of beta_hat on a.

### 7.2 Why this avoids circularity

The one-shot estimator:
1. Is computed **before** any IBSS iteration
2. Uses the raw (unresidualised) beta_hat_j
3. Does not depend on alpha (no feedback loop)
4. Is a standard weighted least squares problem with well-understood properties

The WLS estimate (line 330):

    c_hat(tau^2) = sum_j w_j * beta_hat_j * a_j / sum_j w_j * a_j^2

where w_j = 1/(shat2_j + tau^2). This is consistent when the annotation is informative (E[beta_hat_j] = c_true * a_j for causal SNPs) and the nuisance tau^2 absorbs model misspecification.

### 7.3 Limitations

1. **Multi-effect contamination:** When L > 1, the unconditional beta_hat_j reflects contributions from ALL effects, not just the one being annotated. The one-shot estimator finds c that best explains the **aggregate** signal, which may differ from the per-effect c.

2. **LD contamination:** If causal SNPs are in LD, their correlated beta_hat values inflate the apparent signal for c, potentially overestimating it.

3. **Single c for all effects:** The estimator returns one global c, but different effects might want different c values (e.g., effect 1 is well-annotated, effect 5 is not).

---

## 8. The c-Grid Ensemble: Why It's the Most Robust Strategy

### 8.1 Bayesian model averaging interpretation

The c-grid approach treats c as a discrete hyperparameter:

    p(PIPs | data) = sum_k w_k * PIPs_k

where PIPs_k comes from fitting SuSiNE at c = c_k, and w_k is proportional to exp(ELBO_k) (ELBO-softmax weighting).

### 8.2 Advantages over point estimation

| Property | Iterative EB | One-shot | c-grid |
|----------|-------------|----------|--------|
| Circularity | Dangerous (ranking shifts) | None (pre-fit) | None (c is fixed per fit) |
| Identifiability required | Yes, per-effect per-iter | Yes, marginally | No — BMA marginalizes |
| Multimodality in c | Stuck at one mode | Finds one mode | **Explores all modes** |
| Computational cost | 1 fit (with EB overhead) | 1 fit (+ pre-computation) | K fits |
| Can handle L effects wanting different c | No (one c for all) | No (one c for all) | Partially (ELBO picks best global c) |
| Handles annotation noise | Poorly (amplifies noise) | Moderately (WLS averages) | Well (bad c gets low ELBO weight) |

### 8.3 The key insight

The c-grid doesn't need c to be identifiable. Even if the ELBO surface in c is flat (uninformative annotation), the grid gives approximately uniform weights → the ensemble reduces to averaging over c values → the annotation has negligible effect → graceful degradation to vanilla SuSiE.

In contrast, iterative EB on a flat surface produces **noisy, unreliable point estimates** that can destabilize the entire fit.

---

## 9. A Subtlety: The "mean" Mode Estimates a SCALAR mu_0

An important implementation detail (from `update_priors.R:111-132`): the `"mean"` EB mode optimizes a **single scalar** mu_0 applied uniformly to all SNPs:

    proposed_priors$mu_0 = I(list(rep(proposed_mean, p)))

This is NOT the annotation-structured mu_{0j} = c * a_j. It's finding one global shift that best explains the partial residuals. This is even less useful than estimating c, because:

1. A uniform prior mean shifts all BFs equally → ranking preserved → no annotation leverage
2. It's equivalent to shifting the intercept of the regression, which is already absorbed by centering y

The `"both"` mode has the same issue (line 147: `rep(proposed_mean_var[[1]], p)`).

**Neither iterative mode actually estimates c in the annotation-structured sense.** They estimate a uniform shift, which is the wrong parameterization for leveraging annotations.

### 9.1 What a proper iterative c estimator would look like

A correct iterative EB for c would:
1. Parameterize mu_{0j} = c * a_j (not rep(mu_0, p))
2. Optimize c (scalar) by Brent, with objective neg_log_lik evaluated at mu_0 = c * a
3. Still face the circularity issues from Section 6.3, but at least estimate the right parameter

This has NOT been implemented. The current code only has:
- Uniform-shift EB (wrong parameterization)
- One-shot c estimation (right parameterization, but pre-fit only)
- c-grid ensemble (right parameterization, robust)

---

## 10. When Would Iterative EB for c Be Worth Pursuing?

### 10.1 Theoretical conditions

Iterative EB for c could work when:
1. **High annotation quality** (R^2 > 0.3 between a and true beta at causal sites)
2. **Well-separated effects** (each effect's partial residual is clean)
3. **Consistent sign** (all effects want the same sign of c)
4. **Low LD** (reduces cross-effect contamination)

### 10.2 What would need to change in the code

1. Replace `rep(proposed_mean, p)` with `c_hat * mu_0_unscaled` (use the annotation structure)
2. Change optimization to 1D Brent on c_hat, not on mu_0
3. Add per-effect c tracking to detect sign disagreements
4. Consider a hierarchical model: c_l ~ N(c_global, tau_c^2) where each effect has its own c_l

### 10.3 Is it worth it?

Probably not for the current paper. The c-grid ensemble already handles the key use case (exploring directional priors) with:
- Greater robustness (no circularity)
- Natural uncertainty quantification (ELBO weights)
- Compatibility with the exploration/aggregation framework that is the paper's main contribution

Iterative EB for c would be a follow-up paper's contribution: "Can we reduce the K-fold compute cost by estimating c adaptively?"

---

## 11. Summary of Theoretical Findings

### 11.1 Why EB for sigma_0^2 works

1. sigma_0^2 controls **contrast** (temperature) without changing **ranking** → smooth, unimodal objective
2. Fisher information I(sigma_0^2) is always positive → always identifiable from data
3. IBSS fixed-point with sigma_0^2 EB is approximately **contractive** → stable convergence
4. sigma_0^2 = 0 acts as a natural null hypothesis (automatic sparsity) → well-posed

### 11.2 Why EB for c doesn't work (currently)

1. c changes **ranking** of SNPs → non-smooth, potentially multimodal objective
2. Fisher information I(c) depends on annotation quality → may be near-zero for noisy annotations
3. IBSS fixed-point with c EB can be **expansive** → unstable, positive feedback loops
4. Current implementation estimates a uniform shift (wrong parameterization), not annotation-scaled c
5. Even with correct parameterization, circularity creates multiple fixed points

### 11.3 The c-grid ensemble is optimal for the current paper

- Sidesteps all identifiability/circularity issues
- Provides natural BMA uncertainty quantification
- Degrades gracefully when annotation is uninformative
- Fits within the exploration/aggregation framework
- Computational cost is shared with other exploration axes

---

## 12. Predictions for Experimental Validation (D3B)

The theory predicts:

1. **auto_scale_mu_0 should work well** when annotation R^2 is high and L is small. It should degrade (overestimate |c|) when L is large and effects share LD.

2. **The (c, sigma_0_2) likelihood surface** should be:
   - Unimodal in sigma_0^2 at any fixed c (Section 3)
   - Potentially multimodal or flat in c at any fixed sigma_0^2 (Section 4)
   - The sigma_0^2 ridge should be steeper than the c ridge (higher curvature)

3. **Iterative "mean" mode** should be ineffective because it estimates a uniform shift, not c (Section 9). It should converge to mu_0 ≈ 0 in most cases.

4. **Iterative "both" mode** should either:
   - Collapse to (mu_0 = 0, sigma_0^2 = large) — ignoring the annotation
   - Or find a pathological (mu_0, sigma_0^2) that fits noise rather than signal

5. **Oracle c should give the best AUPRC**, with the c-grid's best coming close and the one-shot estimator being intermediate.

These predictions should be tested in the experiments workbook (D3B).
