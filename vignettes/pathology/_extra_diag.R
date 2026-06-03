## _extra_diag.R -- extra diagnostics for the manuscript writeup.
## Reproduces the four scenarios from pathology_figure_purity.R verbatim
## (same seeds/construction) and prints the quantities the figure script
## does not already echo: per-fit PVE for every (scenario, init), the middle
## decoy x4 marginal PVE, and AUPRC under "drop diffuse effects 4&5".

suppressPackageStartupMessages({ library(susieR) })

make_init <- function(idx, betas, L, P) {
  alpha <- mu <- mu2 <- matrix(0, L, P)
  alpha[cbind(seq_along(idx), idx)] <- 1
  mu[cbind(seq_along(idx), idx)]    <- betas
  mu2[cbind(seq_along(idx), idx)]   <- betas^2
  structure(list(alpha = alpha, mu = mu, mu2 = mu2), class = "susie")
}
fit_three <- function(X, y, true_idx, decoy_idx, b_true, L, b_decoy = NULL, ...) {
  P <- ncol(X)
  res_def  <- susie(X, y, L = L, ...)
  res_true <- susie(X, y, L = L, model_init = make_init(true_idx, b_true, L, P), ...)
  if (is.null(b_decoy)) {
    b_decoy <- as.numeric(coef(lm(y ~ X[, decoy_idx]))[-1])
    b_decoy <- ifelse(abs(b_decoy) < 1e-6, sign(b_decoy + 1e-12) * 1e-6, b_decoy)
  }
  res_decoy <- susie(X, y, L = L, model_init = make_init(decoy_idx, b_decoy, L, P), ...)
  list(Default = res_def, True = res_true, Decoy = res_decoy)
}
fit_pve <- function(fit, X, y) {
  beta_hat <- colSums(fit$alpha * fit$mu)
  var(as.numeric(X %*% beta_hat)) / var(y)
}
auprc <- function(scores, true_idx) {
  ord <- order(scores, decreasing = TRUE)
  yv <- as.integer(seq_along(scores) %in% true_idx)[ord]
  tp <- cumsum(yv); fp <- cumsum(1 - yv)
  rec <- tp / length(true_idx); prec <- tp / pmax(tp + fp, 1)
  sum((rec - c(0, rec[-length(rec)])) * prec)
}
## tie-aware average precision (sklearn average_precision_score semantics):
## one operating point per distinct score; tied variants enter together.
ap_tieaware <- function(scores, true_idx) {
  yv <- as.integer(seq_along(scores) %in% true_idx)
  npos <- length(true_idx)
  uq <- sort(unique(scores), decreasing = TRUE)
  prev_rec <- 0; ap <- 0
  for (t in uq) {
    sel <- scores >= t
    tp <- sum(yv[sel]); rec <- tp / npos; prec <- tp / sum(sel)
    ap <- ap + (rec - prev_rec) * prec
    prev_rec <- rec
  }
  ap
}
cs95 <- function(a) { ord <- order(a, decreasing = TRUE); sort(ord[seq_len(which(cumsum(a[ord]) >= 0.95)[1])]) }
purity_of <- function(idx, R) { if (length(idx) < 2) return(1); sub <- abs(R[idx, idx]); diag(sub) <- NA_real_; min(sub, na.rm = TRUE) }
pip_from_effects <- function(alpha, rows) {
  if (!any(rows)) return(rep(0, ncol(alpha)))
  1 - apply(1 - alpha[rows, , drop = FALSE], 2, prod)
}

N <- 600; P <- 50

report <- function(fits, X, y, label) {
  cat(sprintf("[%s] PVE default=%.4f true=%.4f decoy=%.4f\n",
              label, fit_pve(fits$Default, X, y), fit_pve(fits$True, X, y), fit_pve(fits$Decoy, X, y)))
}

## ---- Scenario 1 ---- ##
set.seed(42)
X1 <- apply(matrix(runif(N * P), N, P), 2, function(x) scale(rbinom(N, 2, x)))
b1 <- c(0.8, -0.8, 1.0); g1 <- X1[, c(1,3,5)] %*% b1
y1 <- g1 + rnorm(N, sd = sqrt(var(g1) * (1/0.3 - 1)))
fits1 <- fit_three(X1, y1, c(1,3,5), c(2,4,6), b1, L = 5,
                   estimate_residual_variance = TRUE, estimate_prior_variance = TRUE)
report(fits1, X1, y1, "Scenario 1")

## ---- Scenario 2 ---- ##
make_centered_basis <- function(n, k, seed = NULL) { if (!is.null(seed)) set.seed(seed); M <- scale(matrix(rnorm(n*k), n, k), center=TRUE, scale=FALSE); qr.Q(qr(M)) * sqrt(n-1) }
make_low_ld_pair <- function(z, u, r, a2 = (1+r)/2) {
  a <- sqrt(a2); b <- sqrt(1-a2); rho <- (r - a2)/(1 - a2); rho <- max(-1, min(1, rho))
  u1 <- u[,1]; u2 <- rho*u1 + sqrt(1-rho^2)*u[,2]; cbind(a*z + b*u1, a*z + b*u2)
}
set.seed(42)
r34 <- 0.10
B2 <- make_centered_basis(N, 120, seed = 42)
X2 <- matrix(0, N, P); z2_1 <- B2[,1]; z2_3 <- B2[,2]; z2_5 <- B2[,3]; eps2 <- B2[,4]
X2[,1:2] <- make_low_ld_pair(z2_1, B2[,5:6], r = 0.98)
X2[,3:4] <- make_low_ld_pair(z2_3, B2[,7:8], r = r34)
X2[,5:6] <- make_low_ld_pair(z2_5, B2[,9:10], r = 0.98)
X2[,7:P] <- B2[,11:(10+P-6)]
X2 <- scale(X2)
s2_pve <- c(0.09, 0.12, 0.05)
y2 <- sqrt(s2_pve[1])*z2_1 + sqrt(s2_pve[2])*z2_3 + sqrt(s2_pve[3])*z2_5 + sqrt(1-sum(s2_pve))*eps2
y2 <- as.numeric(scale(y2))
b2 <- as.numeric(coef(lm(y2 ~ X2[, c(1,3,5)]))[-1])
fits2 <- fit_three(X2, y2, c(1,3,5), c(2,4,6), b2, L = 5, max_iter = 1000)
report(fits2, X2, y2, "Scenario 2")
R2 <- cor(X2)
fit2 <- fits2$Default
pur <- vapply(seq_len(nrow(fit2$alpha)), function(l) purity_of(cs95(fit2$alpha[l,]), R2), numeric(1))
## per-effect summary: owner, cs95 members, effect PVE, purity
for (l in seq_len(nrow(fit2$alpha))) {
  cs <- cs95(fit2$alpha[l, ])
  bl <- fit2$alpha[l, ] * fit2$mu[l, ]
  pvel <- var(as.numeric(X2 %*% bl)) / var(y2)
  cat(sprintf("  S2 effect %d: cs={%s} purity=%.2f effectPVE=%.4f\n", l, paste(cs, collapse=","), pur[l], pvel))
}
keep_pure <- pur >= 0.95
auprc_raw   <- auprc(fit2$pip, c(1,3,5))
auprc_pure  <- auprc(pip_from_effects(fit2$alpha, keep_pure), c(1,3,5))
## "drop diffuse effects 4 and 5" = keep effects 1,2,3
keep_no45   <- c(TRUE, TRUE, TRUE, FALSE, FALSE)
auprc_no45  <- auprc(pip_from_effects(fit2$alpha, keep_no45), c(1,3,5))
cat(sprintf("  S2 marginal PVE: x3(causal)=%.4f x4(decoy)=%.4f | r(x3,x4)=%.3f\n",
            cor(X2[,3], y2)^2, cor(X2[,4], y2)^2, cor(X2[,3], X2[,4])))
cat(sprintf("  S2 marginal PVE all six: x1=%.4f x2=%.4f x3=%.4f x4=%.4f x5=%.4f x6=%.4f\n",
            cor(X2[,1],y2)^2, cor(X2[,2],y2)^2, cor(X2[,3],y2)^2, cor(X2[,4],y2)^2, cor(X2[,5],y2)^2, cor(X2[,6],y2)^2))
cat(sprintf("  S2 AUPRC(index-tie) raw=%.4f drop-diffuse(4,5)=%.4f purity0.95=%.4f | keep_pure={%s}\n",
            auprc_raw, auprc_no45, auprc_pure, paste(which(keep_pure), collapse=",")))
cat(sprintf("  S2 AP(tie-aware)    raw=%.4f drop-diffuse(4,5)=%.4f purity0.95=%.4f\n",
            ap_tieaware(fit2$pip, c(1,3,5)),
            ap_tieaware(pip_from_effects(fit2$alpha, keep_no45), c(1,3,5)),
            ap_tieaware(pip_from_effects(fit2$alpha, keep_pure), c(1,3,5))))
cat("  S2 raw PIP[1:6] hi-prec = ", paste(sprintf("%.8f", fit2$pip[1:6]), collapse=" "), "\n")
cat("  S2 raw PIP[7:10]        = ", paste(sprintf("%.2e", fit2$pip[7:10]), collapse=" "), "\n")
## proper tie-aware AP, rounding scores to a tolerance so 50/50 ties group
ap_round <- function(scores, true_idx, dig = 4) ap_tieaware(round(scores, dig), true_idx)
cat(sprintf("  S2 AP(ties@1e-4) raw=%.4f drop45=%.4f purity=%.4f\n",
            ap_round(fit2$pip, c(1,3,5)),
            ap_round(pip_from_effects(fit2$alpha, keep_no45), c(1,3,5)),
            ap_round(pip_from_effects(fit2$alpha, keep_pure), c(1,3,5))))
cat(sprintf("  S2 AP(ties@1e-2) raw=%.4f drop45=%.4f purity=%.4f\n",
            ap_round(fit2$pip, c(1,3,5), 2),
            ap_round(pip_from_effects(fit2$alpha, keep_no45), c(1,3,5), 2),
            ap_round(pip_from_effects(fit2$alpha, keep_pure), c(1,3,5), 2)))
cat("  S2 raw PIP[1:6]   = ", paste(sprintf("%.3f", fit2$pip[1:6]), collapse=" "), "\n")
cat("  S2 pure PIP[1:6]  = ", paste(sprintf("%.3f", pip_from_effects(fit2$alpha, keep_pure)[1:6]), collapse=" "), "\n")
cat("  S2 pure PIP top non-pair = ", paste(sprintf("%.3f", sort(pip_from_effects(fit2$alpha, keep_pure)[7:P], decreasing=TRUE)[1:3]), collapse=" "), "\n")
cat("  S2 no45 PIP[1:6]  = ", paste(sprintf("%.3f", pip_from_effects(fit2$alpha, keep_no45)[1:6]), collapse=" "), "\n")

## ---- Scenario 3 ---- ##
set.seed(42)
A <- qr.Q(qr(scale(matrix(rnorm(N*3), N, 3), center=TRUE, scale=FALSE)))
Q <- diag(3) - (2/3)*matrix(1,3,3); Bm <- A %*% Q
Z0 <- scale(matrix(rnorm(N*(P-6)), N, P-6), center=TRUE, scale=FALSE)
Z0 <- Z0 - A %*% crossprod(A, Z0)
X3 <- scale(cbind(A[,1], Bm[,1], A[,2], Bm[,2], A[,3], Bm[,3], qr.Q(qr(Z0))))
b3 <- c(1, -1.1, 0.5); g3 <- as.numeric(X3[, c(1,3,5)] %*% b3)
set.seed(115)
y3 <- g3 + rnorm(N, sd = sqrt(var(g3) * (1/0.40 - 1)))
fits3 <- fit_three(X3, y3, c(1,3,5), c(2,4,6), b3, L = 5, b_decoy = as.numeric(t(Q) %*% b3), max_iter = 1000)
report(fits3, X3, y3, "Scenario 3")
e3 <- sapply(fits3, susie_get_objective)
cat(sprintf("  S3 ELBOs: default=%.4f true=%.4f decoy=%.4f range=%.4f\n",
            e3["Default"], e3["True"], e3["Decoy"], diff(range(e3))))
w3 <- exp(e3 - max(e3)); w3 <- w3 / sum(w3)
cat(sprintf("  S3 3-way ELBO-softmax (default,true,decoy) = (%.3f, %.3f, %.3f)\n",
            w3["Default"], w3["True"], w3["Decoy"]))
w3td <- exp(e3["True"] - e3["Decoy"]); w3td <- w3td / (1 + w3td)
cat(sprintf("  S3 2-way softmax(true|decoy) = (%.3f, %.3f)\n", w3td, 1 - w3td))

## ---- Scenario 4 ---- ##
make_signal_proxy <- function(signal, gamma, seed) {
  signal <- as.numeric(scale(signal)); set.seed(seed)
  z <- scale(rnorm(length(signal)), center=TRUE, scale=FALSE)
  z <- z - signal * as.numeric(crossprod(signal, z) / crossprod(signal))
  z <- as.numeric(scale(z)); as.numeric(scale(gamma*signal + sqrt(1-gamma^2)*z))
}
set.seed(123)
f <- rnorm(N)
X_true <- scale(sqrt(0.98)*f + sqrt(0.02)*matrix(rnorm(N*3), N, 3))
beta_true <- c(1.5, -2.2, 0.7); g_true <- as.numeric(scale(X_true %*% beta_true))
d1 <- make_signal_proxy(g_true, 0.35, 124); d2 <- make_signal_proxy(g_true, 0.25, 125); d3 <- make_signal_proxy(g_true, 0.15, 126)
y4 <- g_true + rnorm(N, sd = sqrt(var(g_true) * (1/0.98 - 1)))
Z4 <- scale(matrix(rnorm(N*(P-6)), N, P-6), center=TRUE, scale=FALSE)
base4 <- cbind(X_true[,1], d1, X_true[,2], d2, X_true[,3], d3, y4)
Z4 <- Z4 - base4 %*% solve(crossprod(base4), crossprod(base4, Z4))
X4 <- scale(cbind(X_true[,1], d1, X_true[,2], d2, X_true[,3], d3, qr.Q(qr(Z4))))
fits4 <- fit_three(X4, y4, c(1,3,5), c(2,4,6), beta_true, L = 5, max_iter = 300)
report(fits4, X4, y4, "Scenario 4")
cat(sprintf("  S4 ELBOs: default=%.4f true=%.4f decoy=%.4f\n",
            susie_get_objective(fits4$Default), susie_get_objective(fits4$True), susie_get_objective(fits4$Decoy)))
cat(sprintf("  S4 within-causal r(1,3)=%.3f r(1,5)=%.3f r(3,5)=%.3f\n",
            cor(X4[,1],X4[,3]), cor(X4[,1],X4[,5]), cor(X4[,3],X4[,5])))
