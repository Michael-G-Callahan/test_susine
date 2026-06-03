## purity_pathology_test.R
##
## Scratch test for the SuSiNE purity discussion (NOT a manuscript figure).
##
## Goal: duplicate Scenario 2 but retune the SECOND (causal, decoy) pair
## (cols 3-4) into a "low-purity but honest" credible set:
##   - pairs (1,2) and (5,6) stay at r = 0.98, b = +-1  (purity PASSES, honest)
##   - pair (3,4) is a pure-LD decoy at correlation r34 with a weak effect b_mid
##
## A pure-LD decoy cannot reproduce the causal's fitted-y at low r (that needs
## |r| = 1); the decoy only stays in the CS when its z-score is close to the
## causal's, i.e. z3^2 (1 - r34^2) is small. So we sweep (r34, b_mid) and read
## off, per fit, the CS membership, PIPs on {3,4}, purity, effective size K_l,
## accuracy A_l, and the marginal z-scores -- to locate the regime where the
## set is concentrated + accurate + truth-containing yet purity < 0.95.

suppressPackageStartupMessages({
  library(susieR)
})

set.seed(42)
N <- 600
P <- 10
TRUE_CAUSAL <- c(1, 3, 5)

## --- build a Scenario-2-style design with a tunable pair (3,4) ----------- ##
make_design <- function(r34) {
  X <- matrix(0, N, P)
  ## pair (1,2): high-LD honest pair
  X[, 1] <- rnorm(N); X[, 2] <- 0.98 * X[, 1] + sqrt(1 - 0.98^2) * rnorm(N)
  ## pair (3,4): tunable-correlation pure-LD decoy
  X[, 3] <- rnorm(N); X[, 4] <- r34 * X[, 3] + sqrt(1 - r34^2) * rnorm(N)
  ## pair (5,6): high-LD honest pair
  X[, 5] <- rnorm(N); X[, 6] <- 0.98 * X[, 5] + sqrt(1 - 0.98^2) * rnorm(N)
  ## cols 7-10: independent noise
  X[, 7:P] <- rnorm(N * 4)
  scale(X)
}

## --- effect-level diagnostics from the paper ----------------------------- ##
## K_l = exp(H(alpha^(95)))   (effective size of the top-95% of an effect)
## A_l = max_{j in T} alpha_lj / max_j alpha_lj   (accuracy toward truth)
eff_diag <- function(alpha_row, true_idx) {
  ord <- order(alpha_row, decreasing = TRUE)
  cum <- cumsum(alpha_row[ord])
  k95 <- which(cum >= 0.95)[1]
  sel <- ord[seq_len(k95)]
  p95 <- alpha_row[sel] / sum(alpha_row[sel])
  H   <- -sum(p95 * log(p95))
  K   <- exp(H)
  A   <- max(alpha_row[true_idx]) / max(alpha_row)
  list(K = K, A = A, cs95 = sort(sel))
}

purity_of <- function(idx, R) {
  if (length(idx) < 2) return(1)
  sub <- abs(R[idx, idx]); diag(sub) <- NA
  min(sub, na.rm = TRUE)
}

## --- one fit: report the effect that owns the (3,4) pair ----------------- ##
run_one <- function(r34, b_mid, target_pve = 0.20) {
  X <- make_design(r34)
  R <- cor(X)
  b <- c(1, -b_mid, 1)                       # effects on cols 1,3,5
  g <- X[, TRUE_CAUSAL] %*% b
  y <- as.numeric(g + rnorm(N, sd = sqrt(var(as.numeric(g)) * (1 / target_pve - 1))))

  ## per-variant PVE of the middle causal (col 3), empirically
  pve_mid <- (b_mid^2 * var(X[, 3])) / var(y)
  ## marginal z-scores for cols 3 and 4
  z <- function(j) summary(lm(y ~ X[, j]))$coef[2, 3]
  z3 <- z(3); z4 <- z(4)

  fit <- susie(X, y, L = 5)
  alpha <- fit$alpha
  ## the effect that carries the (3,4) pair
  l_mid <- which.max(alpha[, 3] + alpha[, 4])
  d <- eff_diag(alpha[l_mid, ], TRUE_CAUSAL)
  pur <- purity_of(d$cs95, R)

  data.frame(
    r34      = r34,
    b_mid    = b_mid,
    pve_mid  = round(pve_mid, 4),
    z3       = round(z3, 2),
    z4       = round(z4, 2),
    pip3     = round(fit$pip[3], 3),
    pip4     = round(fit$pip[4], 3),
    cs_size  = length(d$cs95),
    cs_mem   = paste(d$cs95, collapse = ","),
    K_l      = round(d$K, 2),
    A_l      = round(d$A, 3),
    purity   = round(pur, 3),
    filtered = pur < 0.95
  )
}

## --- sweep --------------------------------------------------------------- ##
grid <- expand.grid(
  r34   = c(0.3, 0.5, 0.7, 0.9, 0.98),
  b_mid = c(0.30, 0.45, 0.60, 1.00)
)

res <- do.call(rbind, Map(function(r, bm) {
  set.seed(42)                 # fix design+noise across cells for comparability
  run_one(r, bm)
}, grid$r34, grid$b_mid))

res <- res[order(res$b_mid, res$r34), ]
cat("\n===== purity-pathology sweep (Scenario-2 middle pair) =====\n")
print(res, row.names = FALSE)

cat("\nLegend: cs_mem = members of the (3,4)-owning effect's 95% set;",
    "\n        K_l = effective size, A_l = accuracy (1 = truth on top),",
    "\n        purity = min|r| in that set, filtered = purity < 0.95.\n")
