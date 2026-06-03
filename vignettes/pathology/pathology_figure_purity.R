## pathology_figure.R
##
## Generates the duplicated pathology figure for the SuSiNE manuscript with
## Scenario 2 engineered to make variants 3 and 4 low-LD equal-fit alternatives.
## Output:
##
##   Writings/plots/pathology_combined_purity.png  (7.0 x 8.0 in, 300 dpi)

suppressPackageStartupMessages({
  library(susieR)
  library(ggplot2)
  library(patchwork)
  library(tibble)
})

## --------------------------------------------------------------------- ##
## Helpers                                                                ##
## --------------------------------------------------------------------- ##

make_init <- function(idx, betas, L, P) {
  alpha <- mu <- mu2 <- matrix(0, L, P)
  alpha[cbind(seq_along(idx), idx)] <- 1
  mu[cbind(seq_along(idx), idx)]    <- betas
  mu2[cbind(seq_along(idx), idx)]   <- betas^2
  structure(list(alpha = alpha, mu = mu, mu2 = mu2), class = "susie")
}

fit_three <- function(X, y, true_idx, decoy_idx, b_true, L,
                      b_decoy = NULL, ...) {
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

## ggplot theme: compact, journal-style
theme_paper <- function() {
  theme_bw(base_size = 8) +
    theme(panel.grid.minor = element_blank(),
          panel.grid.major = element_line(linewidth = 0.2, colour = "grey92"),
          axis.title       = element_text(size = 8),
          axis.text        = element_text(size = 7),
          plot.title       = element_blank(),
          legend.position  = "none",
          plot.margin      = margin(2, 2, 2, 2))
}

## PIP panel (one per scenario x init)
fit_pve <- function(fit, X, y) {
  beta_hat <- colSums(fit$alpha * fit$mu)
  var(as.numeric(X %*% beta_hat)) / var(y)
}

pip_panel <- function(fit, true_causal, X, y, display_idx = seq_len(10),
                      show_y = FALSE, show_x = FALSE) {
  display_idx <- display_idx[display_idx <= length(fit$pip)]
  df <- tibble(SNP      = display_idx,
               PIP      = fit$pip[display_idx],
               IsCausal = display_idx %in% true_causal)
  stats_label <- sprintf("PVE %.2f\nELBO %.1f", fit_pve(fit, X, y),
                         susie_get_objective(fit))
  g <- ggplot(df, aes(x = SNP, y = PIP, colour = IsCausal)) +
    geom_point(size = 1.4) +
    annotate("text", x = max(display_idx) + 0.25, y = 0.98, label = stats_label,
             hjust = 1, vjust = 1, size = 2.1, colour = "grey20",
             lineheight = 0.9) +
    scale_colour_manual(values = c(`FALSE` = "grey30", `TRUE` = "#D7261E"),
                        labels = c("Non-causal", "Causal"),
                        name   = NULL) +
    scale_x_continuous(breaks = intersect(c(1, 3, 5, 7, 10), display_idx),
                       limits = range(display_idx) + c(-0.5, 0.5)) +
    scale_y_continuous(limits = c(-0.02, 1.02), breaks = c(0, 0.5, 1)) +
    theme_paper()
  if (!show_y) g <- g + theme(axis.title.y = element_blank(),
                              axis.text.y  = element_blank(),
                              axis.ticks.y = element_blank())
  else         g <- g + labs(y = "PIP")
  if (!show_x) g <- g + theme(axis.title.x = element_blank(),
                              axis.text.x  = element_blank(),
                              axis.ticks.x = element_blank())
  else         g <- g + labs(x = "SNP index")
  g
}

## LD correlation panel (one per scenario)
ld_panel <- function(X, display_idx = seq_len(10), show_y = FALSE, show_x = FALSE) {
  display_idx <- display_idx[display_idx <= ncol(X)]
  R  <- cor(X[, display_idx, drop = FALSE])
  p  <- nrow(R)
  df <- tibble(row = rep(display_idx, times = p),
               col = rep(display_idx, each  = p),
               r   = as.numeric(R))
  g <- ggplot(df, aes(x = col, y = row, fill = r)) +
    geom_tile() +
    scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B",
                         midpoint = 0, limits = c(-1, 1), name = "r",
                         breaks = c(-1, 0, 1)) +
    scale_x_continuous(breaks = intersect(c(1, 5, 10), display_idx),
                       limits = range(display_idx) + c(-0.5, 0.5), expand = c(0, 0)) +
    scale_y_reverse(breaks = intersect(c(1, 5, 10), display_idx),
                    limits = rev(range(display_idx) + c(-0.5, 0.5)), expand = c(0, 0)) +
    coord_equal() +
    theme_paper() +
    theme(panel.grid.major = element_blank(),
          panel.border     = element_rect(colour = "grey60", fill = NA, linewidth = 0.3))
  if (!show_y) g <- g + theme(axis.title.y = element_blank(),
                              axis.text.y  = element_blank(),
                              axis.ticks.y = element_blank())
  else         g <- g + labs(y = "SNP index")
  if (!show_x) g <- g + theme(axis.title.x = element_blank(),
                              axis.text.x  = element_blank(),
                              axis.ticks.x = element_blank())
  else         g <- g + labs(x = "SNP index")
  g
}

## --------------------------------------------------------------------- ##
## Scenario 1 --- independent, strong signal                             ##
## --------------------------------------------------------------------- ##

set.seed(42)
N <- 600; P <- 50
X1 <- apply(matrix(runif(N * P), N, P), 2,
            function(x) scale(rbinom(N, 2, x)))
b1 <- c(0.8, -0.8, 1.0)
g1 <- X1[, c(1, 3, 5)] %*% b1
y1 <- g1 + rnorm(N, sd = sqrt(var(g1) * (1/0.3 - 1)))
fits1 <- fit_three(X1, y1, c(1, 3, 5), c(2, 4, 6), b1, L = 5,
                   estimate_residual_variance = TRUE,
                   estimate_prior_variance    = TRUE)

## --------------------------------------------------------------------- ##
## Scenario 2 --- low-purity but useful middle pair                      ##
## --------------------------------------------------------------------- ##

make_centered_basis <- function(n, k, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  M <- scale(matrix(rnorm(n * k), n, k), center = TRUE, scale = FALSE)
  qr.Q(qr(M)) * sqrt(n - 1)
}

make_low_ld_pair <- function(z, u, r, a2 = (1 + r) / 2) {
  a <- sqrt(a2)
  b <- sqrt(1 - a2)
  rho <- (r - a2) / (1 - a2)
  if (rho < -1 - 1e-10 || rho > 1 + 1e-10) {
    stop("low-LD pair construction is infeasible")
  }
  rho <- max(-1, min(1, rho))
  u1 <- u[, 1]
  u2 <- rho * u1 + sqrt(1 - rho^2) * u[, 2]
  cbind(a * z + b * u1, a * z + b * u2)
}

auprc <- function(scores, true_idx) {
  ord <- order(scores, decreasing = TRUE)
  y <- as.integer(seq_along(scores) %in% true_idx)[ord]
  tp <- cumsum(y)
  fp <- cumsum(1 - y)
  rec <- tp / length(true_idx)
  prec <- tp / pmax(tp + fp, 1)
  sum((rec - c(0, rec[-length(rec)])) * prec)
}

cs95 <- function(alpha_row) {
  ord <- order(alpha_row, decreasing = TRUE)
  sort(ord[seq_len(which(cumsum(alpha_row[ord]) >= 0.95)[1])])
}

purity_of <- function(idx, R) {
  if (length(idx) < 2) return(1)
  sub <- abs(R[idx, idx])
  diag(sub) <- NA_real_
  min(sub, na.rm = TRUE)
}

purity_filtered_pip <- function(fit, X, threshold = 0.95) {
  R <- cor(X)
  pur <- vapply(seq_len(nrow(fit$alpha)), function(l) purity_of(cs95(fit$alpha[l, ]), R), numeric(1))
  keep <- pur >= threshold
  pip <- if (any(keep)) 1 - apply(1 - fit$alpha[keep, , drop = FALSE], 2, prod) else rep(0, ncol(X))
  list(pip = pip, keep = keep, purity = pur)
}

set.seed(42)
r34 <- 0.10
B2 <- make_centered_basis(N, 120, seed = 42)
X2 <- matrix(0, N, P)
z2_1 <- B2[, 1]
z2_3 <- B2[, 2]
z2_5 <- B2[, 3]
eps2 <- B2[, 4]
X2[, 1:2] <- make_low_ld_pair(z2_1, B2[, 5:6], r = 0.98)
X2[, 3:4] <- make_low_ld_pair(z2_3, B2[, 7:8], r = r34)
X2[, 5:6] <- make_low_ld_pair(z2_5, B2[, 9:10], r = 0.98)
X2[, 7:P] <- B2[, 11:(10 + P - 6)]
X2 <- scale(X2)

## The middle pair has r=0.10 but identical marginal fitted-y support.
## At this signal level SuSiE keeps variants 3 and 4 in one 50/50 effect.
s2_component_pve <- c(0.09, 0.12, 0.05)
y2 <- sqrt(s2_component_pve[1]) * z2_1 +
  sqrt(s2_component_pve[2]) * z2_3 +
  sqrt(s2_component_pve[3]) * z2_5 +
  sqrt(1 - sum(s2_component_pve)) * eps2
y2 <- as.numeric(scale(y2))
b2 <- as.numeric(coef(lm(y2 ~ X2[, c(1, 3, 5)]))[-1])

fits2 <- fit_three(X2, y2, c(1, 3, 5), c(2, 4, 6), b2, L = 5,
                   max_iter = 1000)

s2_filter <- purity_filtered_pip(fits2$Default, X2)
s2_raw_auprc <- auprc(fits2$Default$pip, c(1, 3, 5))
s2_filtered_auprc <- auprc(s2_filter$pip, c(1, 3, 5))
s2_owner <- which.max(rowSums(fits2$Default$alpha[, c(3, 4), drop = FALSE]))
s2_pair_share <- fits2$Default$pip[3] / sum(fits2$Default$pip[c(3, 4)])
s2_owner_beta <- fits2$Default$alpha[s2_owner, ] * fits2$Default$mu[s2_owner, ]
s2_owner_pve <- var(as.numeric(X2 %*% s2_owner_beta)) / var(y2)
if (!setequal(cs95(fits2$Default$alpha[s2_owner, ]), c(3, 4)) ||
    s2_pair_share < 0.4 || s2_pair_share > 0.6) {
  stop("Scenario 2 purity construction failed: expected a low-purity {3,4} CS within 60/40 PIP split")
}
cat(sprintf(
  "\n[Scenario 2 purity] r12=%.3f r34=%.3f r56=%.3f x3_marginal_PVE=%.4f cs34_effect_PVE=%.4f PIPs=(%.3f,%.3f,%.3f,%.3f,%.3f,%.3f) raw_AUPRC=%.3f purity_AUPRC=%.3f keep=%s purity={%s}\n",
  cor(X2[, 1], X2[, 2]), cor(X2[, 3], X2[, 4]), cor(X2[, 5], X2[, 6]),
  cor(X2[, 3], y2)^2, s2_owner_pve, fits2$Default$pip[1], fits2$Default$pip[2],
  fits2$Default$pip[3], fits2$Default$pip[4], fits2$Default$pip[5], fits2$Default$pip[6],
  s2_raw_auprc, s2_filtered_auprc,
  paste(as.integer(s2_filter$keep), collapse = ""),
  paste(sprintf("%.2f", s2_filter$purity), collapse = ",")
))

## --------------------------------------------------------------------- ##
## Scenario 3 --- OR-of-ANDs (model-specification barrier)               ##
## --------------------------------------------------------------------- ##

set.seed(42)
A <- qr.Q(qr(scale(matrix(rnorm(N * 3), N, 3), center = TRUE, scale = FALSE)))
Q <- diag(3) - (2/3) * matrix(1, 3, 3)
B <- A %*% Q
Z0 <- scale(matrix(rnorm(N * (P - 6)), N, P - 6), center = TRUE, scale = FALSE)
Z0 <- Z0 - A %*% crossprod(A, Z0)
X3 <- scale(cbind(A[, 1], B[, 1], A[, 2], B[, 2], A[, 3], B[, 3], qr.Q(qr(Z0))))
colnames(X3) <- paste0("V", 1:P)
b3 <- c(1, -1.1, 0.5)
g3 <- as.numeric(X3[, c(1, 3, 5)] %*% b3)
target3_pve <- 0.40
set.seed(115)
y3 <- g3 + rnorm(N, sd = sqrt(var(g3) * (1/target3_pve - 1)))
fits3 <- fit_three(X3, y3, c(1, 3, 5), c(2, 4, 6), b3, L = 5,
                   b_decoy = as.numeric(t(Q) %*% b3),
                   max_iter = 1000)
s3_elbo <- sapply(fits3, susie_get_objective)
if (diff(range(s3_elbo)) > 0.4) {
  stop("Scenario 3 construction failed: all-start ELBO range exceeds 0.4")
}
cat(sprintf(
  "\n[Scenario 3] all-start ELBO range=%.4f PIPs default=(%s) true=(%s) decoy=(%s)\n",
  diff(range(s3_elbo)),
  paste(sprintf("%.3f", fits3$Default$pip[1:6]), collapse = ","),
  paste(sprintf("%.3f", fits3$True$pip[1:6]), collapse = ","),
  paste(sprintf("%.3f", fits3$Decoy$pip[1:6]), collapse = ",")
))

## --------------------------------------------------------------------- ##
## Scenario 4 --- optimizer barrier                                      ##
## --------------------------------------------------------------------- ##

make_signal_proxy <- function(signal, gamma, seed) {
  signal <- as.numeric(scale(signal))
  set.seed(seed)
  z <- scale(rnorm(length(signal)), center = TRUE, scale = FALSE)
  z <- z - signal * as.numeric(crossprod(signal, z) / crossprod(signal))
  z <- as.numeric(scale(z))
  as.numeric(scale(gamma * signal + sqrt(1 - gamma^2) * z))
}

set.seed(123)
f <- rnorm(N)
X_true <- scale(sqrt(0.98) * f + sqrt(0.02) * matrix(rnorm(N * 3), N, 3))
beta_true <- c(1.5, -2.2, 0.7)
g_true <- as.numeric(scale(X_true %*% beta_true))
d1 <- make_signal_proxy(g_true, 0.35, seed = 124)
d2 <- make_signal_proxy(g_true, 0.25, seed = 125)
d3 <- make_signal_proxy(g_true, 0.15, seed = 126)
y4 <- g_true + rnorm(N, sd = sqrt(var(g_true) * (1/0.98 - 1)))
Z4 <- scale(matrix(rnorm(N * (P - 6)), N, P - 6), center = TRUE, scale = FALSE)
base4 <- cbind(X_true[, 1], d1, X_true[, 2], d2, X_true[, 3], d3, y4)
Z4 <- Z4 - base4 %*% solve(crossprod(base4), crossprod(base4, Z4))
X4 <- scale(cbind(X_true[, 1], d1, X_true[, 2], d2, X_true[, 3], d3, qr.Q(qr(Z4))))
colnames(X4) <- paste0("V", 1:P)
fits4 <- fit_three(X4, y4, c(1, 3, 5), c(2, 4, 6), beta_true, L = 5,
                   max_iter = 300)
cat(sprintf(
  "\n[Scenario 4] marginal r(y,X1:6)=(%s)\n",
  paste(sprintf("%.3f", cor(X4, y4)[1:6, 1]), collapse = ",")
))

## --------------------------------------------------------------------- ##
## ELBO diagnostics (printed; used in manuscript body)                   ##
## --------------------------------------------------------------------- ##

report_elbos <- function(fits, label) {
  e <- sapply(fits, susie_get_objective)
  w <- plogis(e["True"] - e["Decoy"])
  cat(sprintf("\n[%s] ELBO default=%.4f true=%.4f decoy=%.4f | softmax(A|B)=(%.3f, %.3f)\n",
              label, e["Default"], e["True"], e["Decoy"], w, 1 - w))
}
report_elbos(fits1, "Scenario 1")
report_elbos(fits2, "Scenario 2")
report_elbos(fits3, "Scenario 3")
report_elbos(fits4, "Scenario 4")

## --------------------------------------------------------------------- ##
## Composite figure                                                       ##
## --------------------------------------------------------------------- ##

build_row <- function(fits, X, y, true_idx, is_last) {
  ld  <- ld_panel(X, show_y = TRUE, show_x = is_last)
  p_d <- pip_panel(fits$Default, true_idx, X, y, show_y = TRUE,  show_x = is_last)
  p_t <- pip_panel(fits$True,    true_idx, X, y, show_y = FALSE, show_x = is_last)
  p_y <- pip_panel(fits$Decoy,   true_idx, X, y, show_y = FALSE, show_x = is_last)
  list(ld = ld, default = p_d, true = p_t, decoy = p_y)
}

row1 <- build_row(fits1, X1, y1, c(1, 3, 5), is_last = FALSE)
row2 <- build_row(fits2, X2, y2, c(1, 3, 5), is_last = FALSE)
row3 <- build_row(fits3, X3, y3, c(1, 3, 5), is_last = FALSE)
row4 <- build_row(fits4, X4, y4, c(1, 3, 5), is_last = TRUE)

## Add column titles (only on row 1)
col_title <- function(lbl) ggplot() + theme_void() +
  annotate("text", x = 0.5, y = 0.5, label = lbl, size = 2.8, fontface = "bold") +
  xlim(0, 1) + ylim(0, 1)

row_label <- function(lbl) ggplot() + theme_void() +
  annotate("text", x = 0.5, y = 0.5, label = lbl, size = 2.8,
           fontface = "bold", angle = 90) +
  xlim(0, 1) + ylim(0, 1)

titles <- (col_title("") | col_title("LD correlation") | col_title("Default") |
           col_title("True warm start") | col_title("Decoy warm start")) +
  plot_layout(widths = c(0.12, 1, 1, 1, 1))

assemble_row <- function(row, label) {
  (row_label(label) | row$ld | row$default | row$true | row$decoy) +
    plot_layout(widths = c(0.12, 1, 1, 1, 1))
}

r1 <- assemble_row(row1, "Scenario 1")
r2 <- assemble_row(row2, "Scenario 2")
r3 <- assemble_row(row3, "Scenario 3")
r4 <- assemble_row(row4, "Scenario 4")

combined <- titles / r1 / r2 / r3 / r4 +
  plot_layout(heights = c(0.12, 1, 1, 1, 1.08),
              guides  = "collect") &
  theme(legend.position  = "bottom",
        legend.box       = "horizontal",
        legend.direction = "horizontal",
        legend.key.height = unit(3, "mm"),
        legend.key.width  = unit(6, "mm"),
        legend.text      = element_text(size = 7),
        legend.title     = element_text(size = 7))

out_path <- file.path(
  "c:/Users/mgcal/OneDrive/Documents/School/Research/Genetics",
  "cTWAS Generalization/Writings/plots/pathology_combined_purity.png"
)

ggsave(out_path, combined, width = 7.0, height = 8.0, dpi = 300, bg = "white")
cat("\nWrote", out_path, "\n")


