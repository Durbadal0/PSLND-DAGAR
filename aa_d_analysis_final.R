# ==============================================================================
# AA-D Stratum Analysis with Model Diagnostics (Final Version)
# ==============================================================================
#
# This script:
#   1. Runs full PS-LND-DAGAR analysis for African American, Distant stage patients
#      with production-level MCMC iterations
#   2. Creates Nelson-Aalen diagnostic plot for Cox-Snell residuals
#   3. Computes CPO (Conditional Predictive Ordinate) diagnostics with single plot
#
# Addressing Reviewer Comment 2: Model diagnostics for log-normal assumption
#
# ==============================================================================

rm(list = ls())
set.seed(12345)

# Load implementation
source("ps_lnd_dagar.R")

# Additional packages for diagnostics
if (!require(survival)) install.packages("survival", repos = "https://cloud.r-project.org")
library(survival)

cat("==============================================================================\n")
cat("AA-D Stratum Analysis with Model Diagnostics (Final Version)\n")
cat("==============================================================================\n\n")

# ==============================================================================
# 1. LOAD AND PREPARE DATA
# ==============================================================================

cat("1. Loading and preparing AA-D data...\n")

# Load FCR data
data <- read.csv("Surv_data2.csv")

# Load Florida county adjacency matrix
fl_adj <- as.matrix(read.csv("W.mat.csv", row.names = 1))

# Filter to AA-D subset (African American, Distant stage)
# Race = 2 (African American), Stage = 3 (Distant)
aa_d <- data[data$Race == 2 & data$Stage == 3, ]

cat("  Total FCR patients:", nrow(data), "\n")
cat("  AA-D subset: N =", nrow(aa_d), "\n")
cat("  Counties with AA-D patients:", length(unique(aa_d$county)), "\n")

# Prepare variables
time <- aa_d$as.numeric.date_diff. / 30.44  # Convert days to months
status <- aa_d$death                         # 1 = death, 0 = censored
Z <- ifelse(aa_d$TX_Delay == 1, 1, 0)        # 1 = short TD, 0 = long TD
county <- aa_d$county

# Covariates
Age <- aa_d$Age
BD <- ifelse(aa_d$BX_Delay == 3, 1, 0)  # 1 = biopsy delay, 0 = no delay
HR <- aa_d$HR_p                          # 1 = positive, 0 = negative
TG <- aa_d$Tgrade                        # 1, 2, 3

X <- cbind(Age, BD, HR, TG)
colnames(X) <- c("Age", "BD", "HR", "TG")

cat("\n  Censoring rate:", round(100 * mean(status == 0), 1), "%\n")
cat("  Short TD (Z=1):", sum(Z == 1), "\n")
cat("  Long TD (Z=0):", sum(Z == 0), "\n")
cat("  Median survival time:", round(median(time), 1), "months\n")
cat("  Time range:", round(min(time), 1), "to", round(max(time), 1), "months\n")

# Get unique counties
unique_counties <- sort(unique(county))
K_data <- length(unique_counties)
cat("  Counties in data:", K_data, "out of 67\n\n")

# ==============================================================================
# 2. FIT PS-LND-DAGAR MODEL (Production-level iterations)
# ==============================================================================

cat("2. Fitting PS-LND-DAGAR model...\n")
cat("   MCMC Configuration:\n")

# MCMC settings - production-level iterations (40k total)
# Stage 1: Propensity Score Model
n_mcmc_ps <- 20000    # Total PS iterations
burn_in_ps <- 10000   # PS burn-in
# Stage 2: Outcome Model
n_mcmc_out <- 40000   # Total outcome iterations
burn_in_out <- 10000  # Outcome burn-in
thin <- 10            # Thinning interval

n_samples_ps <- (n_mcmc_ps - burn_in_ps) / thin
n_samples_out <- (n_mcmc_out - burn_in_out) / thin

cat("   Stage 1 (PS): ", n_mcmc_ps, " iterations, ", burn_in_ps, " burn-in, thin=", thin,
    " -> ", n_samples_ps, " samples\n", sep="")
cat("   Stage 2 (Outcome): ", n_mcmc_out, " iterations, ", burn_in_out, " burn-in, thin=", thin,
    " -> ", n_samples_out, " samples\n", sep="")
cat("   (This will take approximately 45-60 minutes)\n\n")

fit <- ps_lnd_dagar(
  time = time,
  status = status,
  X = X,
  Z = Z,
  cluster_id = county,
  adj_matrix = fl_adj,
  ordering = 1:67,
  n_mcmc_ps = n_mcmc_ps,
  burn_in_ps = burn_in_ps,
  n_mcmc_out = n_mcmc_out,
  burn_in_out = burn_in_out,
  thin = thin,
  num_tree_ps = 50,
  num_tree_out = 50,
  verbose = TRUE
)

cat("\n")

# ==============================================================================
# 3. MODEL SUMMARY
# ==============================================================================

cat("3. Model Summary...\n")

n_samples <- fit$out_fit$n_save
cat("   Posterior samples retained:", n_samples, "\n")
cat("   Spatial rho: mean =", round(mean(fit$out_fit$rho_samples), 3),
    ", 95% CI = [", round(quantile(fit$out_fit$rho_samples, 0.025), 3), ",",
    round(quantile(fit$out_fit$rho_samples, 0.975), 3), "]\n")
cat("   Error sigma: mean =", round(mean(sqrt(fit$out_fit$sigma2_samples)), 3),
    ", 95% CI = [", round(quantile(sqrt(fit$out_fit$sigma2_samples), 0.025), 3), ",",
    round(quantile(sqrt(fit$out_fit$sigma2_samples), 0.975), 3), "]\n")
cat("   Spatial tau_w: mean =", round(mean(fit$out_fit$tau_w_samples), 3), "\n")
cat("   PS equicorrelation rho_0: mean =", round(mean(fit$ps_fit$rho_0_samples), 3), "\n\n")

# ==============================================================================
# 4. NELSON-AALEN DIAGNOSTIC PLOT (Cox-Snell Residuals)
# ==============================================================================

cat("4. Creating Nelson-Aalen diagnostic plot...\n")

# --------------------------------------------------
# Cox-Snell Residuals for Log-Normal AFT Model
# --------------------------------------------------
# If model is correct, H(Y_i | X_i, theta) ~ Exp(1) (censored)
#
# For log-normal:
#   S(t) = 1 - Phi((log(t) - mu) / sigma)
#   H(t) = -log(S(t))
#
# Cox-Snell residual: r_i = H(y_i | mu_i, sigma)
# --------------------------------------------------

N <- length(time)

# Get posterior means for predictions
sigma_mean <- mean(sqrt(fit$out_fit$sigma2_samples))
f_mean <- colMeans(fit$out_fit$f_samples)
W_mean <- colMeans(fit$out_fit$w_samples)

# Map counties to internal indices
county_internal <- fit$orig_to_internal[as.character(county)]

# Compute mu_i = f_i + W_{county_i} for each observation (posterior mean)
mu_hat <- numeric(N)
for (i in 1:N) {
  mu_hat[i] <- f_mean[i] + W_mean[county_internal[i]]
}

# Compute Cox-Snell residuals
# r_i^CS = -log(S(y_i)) = -log(1 - Phi((log(y_i) - mu_i) / sigma))
log_time <- log(pmax(time, 0.001))  # Avoid log(0)
z_scores <- (log_time - mu_hat) / sigma_mean
S_hat <- 1 - pnorm(z_scores)
S_hat <- pmax(S_hat, 1e-10)  # Avoid log(0)
cox_snell_resid <- -log(S_hat)

# Create survival object for the Cox-Snell residuals
surv_resid <- Surv(cox_snell_resid, status)

# Fit Nelson-Aalen estimator to residuals
na_fit <- survfit(surv_resid ~ 1, type = "fleming-harrington")

# Extract Nelson-Aalen cumulative hazard estimate
na_time <- na_fit$time
na_cumhaz <- na_fit$cumhaz

# Create Nelson-Aalen diagnostic plot
pdf("aa_d_nelson_aalen_diagnostic.pdf", width = 8, height = 7)
par(mar = c(5, 5, 4, 2))

plot(na_time, na_cumhaz, type = "s", lwd = 2, col = "blue",
     xlab = expression("Cox-Snell Residual " * (r[i]^{CS})),
     ylab = expression("Nelson-Aalen Cumulative Hazard " * hat(H)(r^{CS})),
     main = "Nelson-Aalen Diagnostic for Log-Normal AFT Model",
     xlim = c(0, max(na_time) * 1.05),
     ylim = c(0, max(na_cumhaz) * 1.05),
     cex.lab = 1.2, cex.main = 1.1)

# Add 45-degree reference line (theoretical Exp(1))
abline(0, 1, col = "red", lwd = 2, lty = 2)

# Add legend
legend("topleft",
       legend = c("Nelson-Aalen estimate", "Theoretical Exp(1): H(r) = r"),
       col = c("blue", "red"), lwd = 2, lty = c(1, 2), bty = "n", cex = 1.0)

dev.off()
cat("   Saved: aa_d_nelson_aalen_diagnostic.pdf\n")

# Create Cox-Snell residual plot (separate)
pdf("aa_d_cox_snell_residuals.pdf", width = 8, height = 7)
par(mar = c(5, 5, 4, 2))

plot(1:N, cox_snell_resid,
     pch = ifelse(status == 1, 16, 1),
     col = ifelse(status == 1, "darkgreen", "red"),
     xlab = "Observation Index",
     ylab = expression("Cox-Snell Residual " * (r[i]^{CS})),
     main = "Cox-Snell Residuals",
     cex = 0.8, cex.lab = 1.2, cex.main = 1.1)

# Add reference line at 1 (mean of Exp(1))
abline(h = 1, col = "gray40", lty = 2, lwd = 1.5)

# Add legend
legend("topright",
       legend = c("Observed (death)", "Censored"),
       pch = c(16, 1), col = c("darkgreen", "red"), bty = "n", cex = 1.0)

dev.off()
cat("   Saved: aa_d_cox_snell_residuals.pdf\n\n")

# ==============================================================================
# 5. CPO (CONDITIONAL PREDICTIVE ORDINATE) DIAGNOSTICS
# ==============================================================================

cat("5. Computing CPO diagnostics...\n")
cat("   (This requires computing likelihoods for all N x M combinations)\n")

# --------------------------------------------------
# CPO Calculation using Harmonic Mean Estimator
# --------------------------------------------------
# CPO_i = [E_{theta|D}(1/f(Y_i|theta))]^{-1}
#
# For log-normal AFT with censoring:
# If delta_i = 1 (observed):
#   f(y_i|theta) = (1/(y_i * sigma * sqrt(2*pi))) * exp(-z_i^2/2)
#   where z_i = (log(y_i) - mu_i) / sigma
#
# If delta_i = 0 (censored):
#   f(y_i|theta) = S(y_i|theta) = 1 - Phi(z_i)
# --------------------------------------------------

# Extract MCMC samples
f_samples <- fit$out_fit$f_samples      # n_samples x N
W_samples <- fit$out_fit$w_samples      # n_samples x K
sigma2_samples <- fit$out_fit$sigma2_samples  # n_samples

# Storage for inverse likelihoods (for harmonic mean)
inv_lik_matrix <- matrix(0, nrow = n_samples, ncol = N)

cat("   Computing likelihoods across", n_samples, "MCMC samples for", N, "observations...\n")

pb_interval <- max(1, floor(n_samples / 10))

for (m in 1:n_samples) {
  if (m %% pb_interval == 0) cat("     Progress:", round(100 * m / n_samples), "%\n")

  sigma_m <- sqrt(sigma2_samples[m])
  f_m <- f_samples[m, ]
  W_m <- W_samples[m, ]

  for (i in 1:N) {
    # Compute mu_i for this MCMC sample
    mu_i <- f_m[i] + W_m[county_internal[i]]
    z_i <- (log_time[i] - mu_i) / sigma_m

    if (status[i] == 1) {
      # Observed event: use log-normal density
      # f(y) = (1/(y*sigma*sqrt(2*pi))) * exp(-z^2/2)
      # log f(y) = -log(y) - log(sigma) - 0.5*log(2*pi) - 0.5*z^2
      log_lik <- -log(time[i]) - log(sigma_m) - 0.5 * log(2 * pi) - 0.5 * z_i^2
    } else {
      # Censored: use survival function S(y) = 1 - Phi(z)
      log_lik <- pnorm(z_i, lower.tail = FALSE, log.p = TRUE)
    }

    # Store inverse likelihood for harmonic mean calculation
    inv_lik_matrix[m, i] <- exp(-log_lik)
  }
}

# Compute CPO as inverse of harmonic mean of likelihoods
# CPO_i = 1 / mean(1/f(y_i|theta^(m)))
harmonic_mean_inv <- colMeans(inv_lik_matrix)
CPO <- 1 / harmonic_mean_inv
log_CPO <- log(CPO)

# Compute LPML (Log Pseudo Marginal Likelihood)
LPML <- sum(log_CPO)

cat("\n   CPO Summary Statistics:\n")
cat("     Mean CPO:    ", sprintf("%.6f", mean(CPO)), "\n")
cat("     Median CPO:  ", sprintf("%.6f", median(CPO)), "\n")
cat("     Min CPO:     ", sprintf("%.6f", min(CPO)), "\n")
cat("     Max CPO:     ", sprintf("%.6f", max(CPO)), "\n")
cat("     Mean log(CPO):", sprintf("%.3f", mean(log_CPO)), "\n")
cat("     LPML:        ", sprintf("%.2f", LPML), "\n\n")

# Create single CPO diagnostic plot (colored by censoring status only)
pdf("aa_d_cpo_diagnostic.pdf", width = 9, height = 7)
par(mar = c(5, 5, 4, 2))

# Plot log(CPO) vs survival time
plot(time, log_CPO,
     pch = ifelse(status == 1, 16, 1),
     col = ifelse(status == 1, "darkblue", "darkorange"),
     xlab = "Survival Time (months)",
     ylab = expression(log(CPO[i])),
     main = "CPO Diagnostic for Log-Normal AFT Model",
     cex = 0.9, cex.lab = 1.2, cex.main = 1.1)

# Add legend
legend("bottomright",
       legend = c("Observed (death)", "Censored"),
       pch = c(16, 1),
       col = c("darkblue", "darkorange"),
       bty = "n", cex = 1.0)

dev.off()
cat("   Saved: aa_d_cpo_diagnostic.pdf\n\n")

# ==============================================================================
# 6. SAVE RESULTS
# ==============================================================================

cat("6. Saving results...\n")

# Save all diagnostic data
diagnostics <- list(
  # Data info
  N = N,
  K = K_data,
  time = time,
  status = status,
  censoring_rate = mean(status == 0),

  # MCMC configuration
  mcmc_config = list(
    n_mcmc_ps = n_mcmc_ps, burn_in_ps = burn_in_ps,
    n_mcmc_out = n_mcmc_out, burn_in_out = burn_in_out,
    thin = thin, n_samples = n_samples
  ),

  # Model fit
  fit = fit,
  sigma_mean = sigma_mean,
  rho_mean = mean(fit$out_fit$rho_samples),
  tau_w_mean = mean(fit$out_fit$tau_w_samples),

  # Cox-Snell residuals
  cox_snell_resid = cox_snell_resid,
  mu_hat = mu_hat,

  # Nelson-Aalen
  na_time = na_time,
  na_cumhaz = na_cumhaz,

  # CPO diagnostics
  CPO = CPO,
  log_CPO = log_CPO,
  LPML = LPML
)

save(diagnostics, file = "aa_d_diagnostics_final.RData")
cat("   Saved: aa_d_diagnostics_final.RData\n\n")

# ==============================================================================
# 7. SUMMARY
# ==============================================================================

cat("==============================================================================\n")
cat("ANALYSIS COMPLETE\n")
cat("==============================================================================\n\n")

cat("MCMC Configuration:\n")
cat("  Stage 1 (PS):     ", n_mcmc_ps, " iterations (", burn_in_ps, " burn-in)\n", sep="")
cat("  Stage 2 (Outcome): ", n_mcmc_out, " iterations (", burn_in_out, " burn-in)\n", sep="")
cat("  Thinning:          ", thin, "\n", sep="")
cat("  Posterior samples: ", n_samples, "\n\n", sep="")

cat("Model Parameter Estimates:\n")
cat("  Spatial rho:       ", round(mean(fit$out_fit$rho_samples), 3),
    " [", round(quantile(fit$out_fit$rho_samples, 0.025), 3), ", ",
    round(quantile(fit$out_fit$rho_samples, 0.975), 3), "]\n", sep="")
cat("  Error sigma:       ", round(sigma_mean, 3),
    " [", round(quantile(sqrt(fit$out_fit$sigma2_samples), 0.025), 3), ", ",
    round(quantile(sqrt(fit$out_fit$sigma2_samples), 0.975), 3), "]\n", sep="")
cat("  Spatial tau_w:     ", round(mean(fit$out_fit$tau_w_samples), 3), "\n", sep="")
cat("  PS rho_0:          ", round(mean(fit$ps_fit$rho_0_samples), 3), "\n\n", sep="")

cat("Diagnostic Summary:\n")
cat("  LPML:              ", round(LPML, 2), "\n", sep="")
cat("  Mean log(CPO):     ", round(mean(log_CPO), 3), "\n\n", sep="")

cat("Output files:\n")
cat("  - aa_d_nelson_aalen_diagnostic.pdf : Nelson-Aalen diagnostic plot\n")
cat("  - aa_d_cox_snell_residuals.pdf     : Cox-Snell residuals plot\n")
cat("  - aa_d_cpo_diagnostic.pdf          : CPO diagnostic plot\n")
cat("  - aa_d_diagnostics_final.RData     : All diagnostic data\n\n")

cat("Interpretation:\n")
cat("  1. Nelson-Aalen Plot: If Cox-Snell residuals follow the 45-degree line,\n")
cat("     the log-normal assumption is reasonable for the observed data.\n\n")

cat("  2. CPO Plot: No systematic pattern suggests good model fit.\n")
cat("     Low CPO values indicate observations that are poorly predicted.\n")
cat("     LPML = ", round(LPML, 2), " (higher is better for model comparison).\n\n", sep="")

cat("  3. For RMST-based causal inference (restricted to t* = 10 years),\n")
cat("     the log-normal tail behavior is less critical.\n\n")

cat("==============================================================================\n")
