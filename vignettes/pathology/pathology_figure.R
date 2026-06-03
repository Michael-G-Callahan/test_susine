## pathology_figure.R
##
## Generates the combined pathology figure for the SuSiNE manuscript
## (four scenarios x four panels: LD + three initializations).
##
## Reproduces the four scenarios from susie_pathology.ipynb verbatim
## (same seeds, same generation, same fits) and composites the per-panel
## ggplots with patchwork. Output:
##
##   Writings/plots/pathology_combined.png  (7.0 x 8.0 in, 300 dpi)

suppressPackageStartupMessages({
  library(susieR)
  library(ggplot2)
  library(patchwork)
  library(dplyr)
  library(tidyr)
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
pip_panel <- function(fit, true_causal, show_y = FALSE, show_x = FALSE) {
  df <- tibble(SNP      = seq_along(fit$pip),
               PIP      = fit$pip,
               IsCausal = seq_along(fit$pip) %in% true_causal)
  g <- ggplot(df, aes(x = SNP, y = PIP, colour = IsCausal)) +
    geom_point(size = 1.4) +
    scale_colour_manual(values = c(`FALSE` = "grey30", `TRUE` = "#D7261E"),
                        labels = c("Non-causal", "Causal"),
                        name   = NULL) +
    scale_x_continuous(breaks = c(1, 3, 5, 7, 10), limits = c(0.5, 10.5)) +
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
ld_panel <- function(X, show_y = FALSE, show_x = FALSE) {
  R  <- cor(X)
  p  <- nrow(R)
  df <- tibble(row = rep(seq_len(p), times = p),
               col = rep(seq_len(p), each  = p),
               r   = as.numeric(R))
  g <- ggplot(df, aes(x = col, y = row, fill = r)) +
    geom_tile() +
    scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B",
                         midpoint = 0, limits = c(-1, 1), name = "r",
                         breaks = c(-1, 0, 1)) +
    scale_x_continuous(breaks = c(1, 5, 10), expand = c(0, 0)) +
    scale_y_reverse(breaks = c(1, 5, 10), expand = c(0, 0)) +
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
N <- 600; P <- 10
X1 <- apply(matrix(runif(N * P), N, P), 2,
            function(x) scale(rbinom(N, 2, x)))
b1 <- c(0.8, -0.8, 1.0)
g1 <- X1[, c(1, 3, 5)] %*% b1
y1 <- g1 + rnorm(N, sd = sqrt(var(g1) * (1/0.3 - 1)))
fits1 <- fit_three(X1, y1, c(1, 3, 5), c(2, 4, 6), b1, L = 5,
                   estimate_residual_variance = TRUE,
                   estimate_prior_variance    = TRUE)

## --------------------------------------------------------------------- ##
## Scenario 2 --- within-pair high LD (AND-of-ORs)                       ##
## --------------------------------------------------------------------- ##

set.seed(42)
X2 <- matrix(0, N, P)
for (i in seq(1, 6, 2)) {
  X2[, i]     <- rnorm(N)
  X2[, i + 1] <- 0.98 * X2[, i] + sqrt(1 - 0.98^2) * rnorm(N)
}
X2[, 7:P] <- rnorm(N * 4)
X2 <- scale(X2)
b2 <- c(1, -1, 1)
g2 <- X2[, c(1, 3, 5)] %*% b2
y2 <- g2 + rnorm(N, sd = sqrt(var(g2) * (1/0.2 - 1)))
fits2 <- fit_three(X2, y2, c(1, 3, 5), c(2, 4, 6), b2, L = 5)

## --------------------------------------------------------------------- ##
## Scenario 3 --- OR-of-ANDs (model-specification barrier)               ##
## --------------------------------------------------------------------- ##

set.seed(42)
A <- qr.Q(qr(scale(matrix(rnorm(N * 3), N, 3), center = TRUE, scale = FALSE)))
Q <- diag(3) - (2/3) * matrix(1, 3, 3)
B <- A %*% Q
Z0 <- scale(matrix(rnorm(N * 4), N, 4), center = TRUE, scale = FALSE)
Z0 <- Z0 - A %*% crossprod(A, Z0)
X3 <- scale(cbind(A[, 1], B[, 1], A[, 2], B[, 2], A[, 3], B[, 3], qr.Q(qr(Z0))))
colnames(X3) <- c("A1", "B1", "A2", "B2", "A3", "B3", "V7", "V8", "V9", "V10")
b3 <- c(1, -1, 0.5)
g3 <- X3[, c(1, 3, 5)] %*% b3
y3 <- g3 + rnorm(N, sd = sqrt(var(g3) * (1/0.7 - 1)))
fits3 <- fit_three(X3, y3, c(1, 3, 5), c(2, 4, 6), b3, L = 5,
                   b_decoy = as.numeric(t(Q) %*% b3),
                   estimate_prior_variance = FALSE,
                   scaled_prior_variance   = 0.2)

## --------------------------------------------------------------------- ##
## Scenario 4 --- optimizer barrier                                      ##
## --------------------------------------------------------------------- ##

set.seed(123)
f <- rnorm(N)
X_true <- scale(sqrt(0.8) * f + sqrt(0.2) * matrix(rnorm(N * 3), N, 3))
beta_true <- c(1.2, -1.8, 0.6)
g_true <- as.numeric(scale(X_true %*% beta_true))
d1 <- as.numeric(scale(g_true + 0.50 * rnorm(N)))
d2 <- as.numeric(scale(d1     + 1.00 * rnorm(N)))
d3 <- as.numeric(scale(d2     + 1.30 * rnorm(N)))
X4 <- scale(cbind(X_true[, 1], d1, X_true[, 2], d2, X_true[, 3], d3,
                  matrix(rnorm(N * 4), N, 4)))
colnames(X4) <- paste0("V", 1:P)
y4 <- g_true + rnorm(N, sd = sqrt(var(g_true) * (1/0.70 - 1)))
fits4 <- fit_three(X4, y4, c(1, 3, 5), c(2, 4, 6), beta_true, L = 3,
                   max_iter = 200)

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

build_row <- function(fits, X, true_idx, is_last) {
  ld  <- ld_panel(X, show_y = TRUE, show_x = is_last)
  p_d <- pip_panel(fits$Default, true_idx, show_y = TRUE,  show_x = is_last)
  p_t <- pip_panel(fits$True,    true_idx, show_y = FALSE, show_x = is_last)
  p_y <- pip_panel(fits$Decoy,   true_idx, show_y = FALSE, show_x = is_last)
  list(ld = ld, default = p_d, true = p_t, decoy = p_y)
}

row1 <- build_row(fits1, X1, c(1, 3, 5), is_last = FALSE)
row2 <- build_row(fits2, X2, c(1, 3, 5), is_last = FALSE)
row3 <- build_row(fits3, X3, c(1, 3, 5), is_last = FALSE)
row4 <- build_row(fits4, X4, c(1, 3, 5), is_last = TRUE)

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
  "cTWAS Generalization/Writings/plots/pathology_combined.png"
)

ggsave(out_path, combined, width = 7.0, height = 8.0, dpi = 300, bg = "white")
cat("\nWrote", out_path, "\n")
