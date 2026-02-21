# ==============================================================================
# AA-D Stratum Analysis with Model Diagnostics
# ==============================================================================
#
# This script:
#   1. Runs full PS-LND-DAGAR analysis for African American, Distant stage patients
#   2. Creates Nelson-Aalen/Kaplan-Meier diagnostic plots
#   3. Computes CPO (Conditional Predictive Ordinate) diagnostics
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
cat("AA-D Stratum Analysis with Model Diagnostics\n")
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
# 2. FIT PS-LND-DAGAR MODEL (Long iterations)
# ==============================================================================

cat("2. Fitting PS-LND-DAGAR model with long MCMC...\n")
cat("   (This will take several minutes)\n\n")

# MCMC settings - long iterations for final analysis
n_mcmc_ps <- 3000    # Propensity score MCMC iterations
burn_in_ps <- 1000   # PS burn-in
n_mcmc_out <- 5000   # Outcome model MCMC iterations
burn_in_out <- 2000  # Outcome burn-in
thin <- 5            # Thinning interval

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
cat("   Posterior samples:", n_samples, "\n")
cat("   Spatial rho: mean =", round(mean(fit$out_fit$rho_samples), 3),
    ", 95% CI = [", round(quantile(fit$out_fit$rho_samples, 0.025), 3), ",",
    round(quantile(fit$out_fit$rho_samples, 0.975), 3), "]\n")
cat("   Error sigma: mean =", round(mean(sqrt(fit$out_fit$sigma2_samples)), 3),
    ", 95% CI = [", round(quantile(sqrt(fit$out_fit$sigma2_samples), 0.025), 3), ",",
    round(quantile(sqrt(fit$out_fit$sigma2_samples), 0.975), 3), "]\n")
cat("   PS equicorrelation rho_0: mean =", round(mean(fit$ps_fit$rho_0_samples), 3), "\n\n")

# ==============================================================================
# 4. NELSON-AALEN / KAPLAN-MEIER DIAGNOSTIC PLOT
# ==============================================================================

cat("4. Creating Nelson-Aalen/Kaplan-Meier diagnostic plot...\n")

# The idea: If the log-normal AFT model is correct, then the cumulative hazard
# H(t|X) = -log(S(t|X)) should follow a specific pattern.
#
# For log-normal: H(t) = -log(1 - Phi((log(t) - mu)/sigma))
#
# Transform: H_i = -log(S(Y_i | X_i, theta))
# If model is correct, these should be censored samples from Exp(1)
# => Nelson-Aalen estimator should be approximately linear

N <- length(time)

# Get posterior means for predictions
sigma_mean <- mean(sqrt(fit$out_fit$sigma2_samples))
f_mean <- colMeans(fit$out_fit$f_samples)
W_mean <- colMeans(fit$out_fit$w_samples)

# Map counties to internal indices
county_internal <- fit$orig_to_internal[as.character(county)]

# Compute mu_i = f_i + W_{county_i} for each observation
mu_hat <- numeric(N)
for (i in 1:N) {
  mu_hat[i] <- f_mean[i] + W_mean[county_internal[i]]
}

# Compute cumulative hazard residuals (Cox-Snell residuals)
# For log-normal: S(t) = 1 - Phi((log(t) - mu) / sigma)
# H(t) = -log(S(t))

log_time <- log(pmax(time, 0.001))  # Avoid log(0)
z_scores <- (log_time - mu_hat) / sigma_mean
S_hat <- 1 - pnorm(z_scores)
S_hat <- pmax(S_hat, 1e-10)  # Avoid log(0)
H_residuals <- -log(S_hat)

# Cox-Snell residuals: if model is correct, these should be Exp(1) censored
# Create survival object for the residuals
surv_resid <- Surv(H_residuals, status)

# Fit Nelson-Aalen estimator to residuals
na_fit <- survfit(surv_resid ~ 1, type = "fleming-harrington")

# Extract cumulative hazard from NA estimator
na_time <- na_fit$time
na_cumhaz <- na_fit$cumhaz

# Also compute Kaplan-Meier for visual comparison
km_fit <- survfit(surv_resid ~ 1)

# Create diagnostic plot
pdf("aa_d_na_km_diagnostic.pdf", width = 12, height = 5)
par(mfrow = c(1, 3), mar = c(5, 5, 4, 2))

# Panel 1: Nelson-Aalen cumulative hazard vs. theoretical Exp(1)
plot(na_time, na_cumhaz, type = "s", lwd = 2, col = "blue",
     xlab = "Cox-Snell Residuals (Cumulative Hazard)",
     ylab = "Nelson-Aalen Cumulative Hazard",
     main = "Nelson-Aalen Diagnostic\n(Should be y = x if model correct)",
     xlim = c(0, max(na_time) * 1.1),
     ylim = c(0, max(na_cumhaz) * 1.1))
abline(0, 1, col = "red", lwd = 2, lty = 2)
legend("topleft", legend = c("Nelson-Aalen estimate", "Theoretical Exp(1)"),
       col = c("blue", "red"), lwd = 2, lty = c(1, 2), bty = "n")

# Panel 2: Scatter of Cox-Snell residuals by censoring status
plot(1:N, H_residuals, pch = ifelse(status == 1, 16, 1),
     col = ifelse(status == 1, "darkgreen", "red"),
     xlab = "Observation Index", ylab = "Cox-Snell Residual",
     main = "Cox-Snell Residuals by Censoring Status",
     cex = 0.8)
abline(h = 1, col = "gray", lty = 2)
legend("topright", legend = c("Observed (death)", "Censored"),
       pch = c(16, 1), col = c("darkgreen", "red"), bty = "n")

# Panel 3: Q-Q plot of uncensored residuals vs Exp(1)
uncensored_resid <- H_residuals[status == 1]
n_uncens <- length(uncensored_resid)
theoretical_quantiles <- qexp(ppoints(n_uncens))
empirical_quantiles <- sort(uncensored_resid)

plot(theoretical_quantiles, empirical_quantiles, pch = 16, cex = 0.6,
     xlab = "Theoretical Exp(1) Quantiles",
     ylab = "Empirical Cox-Snell Residual Quantiles",
     main = "Q-Q Plot (Uncensored Observations Only)")
abline(0, 1, col = "red", lwd = 2)

dev.off()
cat("   Saved: aa_d_na_km_diagnostic.pdf\n\n")

# ==============================================================================
# 5. CPO (CONDITIONAL PREDICTIVE ORDINATE) DIAGNOSTICS
# ==============================================================================

cat("5. Computing CPO diagnostics...\n")

# CPO_i = [E_{theta|D}(1/f(Y_i|theta))]^{-1}
# This is the harmonic mean of the likelihoods across MCMC samples
#
# For log-normal AFT with censoring:
# - If delta_i = 1 (observed): f(Y_i|theta) = (1/sigma*Y_i) * phi((log(Y_i) - mu_i)/sigma)
# - If delta_i = 0 (censored): f(Y_i|theta) = S(Y_i|theta) = 1 - Phi((log(Y_i) - mu_i)/sigma)

# Extract MCMC samples
f_samples <- fit$out_fit$f_samples      # n_samples x N
W_samples <- fit$out_fit$w_samples      # n_samples x K
sigma2_samples <- fit$out_fit$sigma2_samples  # n_samples

# Storage for inverse likelihoods
inv_lik_matrix <- matrix(0, nrow = n_samples, ncol = N)

cat("   Computing likelihoods across", n_samples, "MCMC samples...\n")

for (m in 1:n_samples) {
  if (m %% 100 == 0) cat("     Sample", m, "/", n_samples, "\n")

  sigma_m <- sqrt(sigma2_samples[m])
  f_m <- f_samples[m, ]
  W_m <- W_samples[m, ]

  for (i in 1:N) {
    mu_i <- f_m[i] + W_m[county_internal[i]]
    z_i <- (log_time[i] - mu_i) / sigma_m

    if (status[i] == 1) {
      # Observed: log-normal density
      # f(y) = (1/(y*sigma*sqrt(2*pi))) * exp(-z^2/2)
      log_lik <- -log(time[i]) - log(sigma_m) - 0.5 * log(2 * pi) - 0.5 * z_i^2
    } else {
      # Censored: survival function
      # S(y) = 1 - Phi(z)
      log_lik <- pnorm(z_i, lower.tail = FALSE, log.p = TRUE)
    }

    # Store inverse likelihood (for harmonic mean)
    inv_lik_matrix[m, i] <- exp(-log_lik)
  }
}

# Compute CPO as harmonic mean
harmonic_mean_inv <- colMeans(inv_lik_matrix)
CPO <- 1 / harmonic_mean_inv
log_CPO <- log(CPO)

# Compute LPML (Log Pseudo Marginal Likelihood)
LPML <- sum(log_CPO)

cat("\n   CPO Summary:\n")
cat("     Mean CPO:", round(mean(CPO), 6), "\n")
cat("     Median CPO:", round(median(CPO), 6), "\n")
cat("     Min CPO:", round(min(CPO), 6), "\n")
cat("     Max CPO:", round(max(CPO), 6), "\n")
cat("     LPML:", round(LPML, 2), "\n\n")

# Create CPO diagnostic plots
pdf("aa_d_cpo_diagnostic.pdf", width = 12, height = 5)
par(mfrow = c(1, 3), mar = c(5, 5, 4, 2))

# Panel 1: CPO vs survival time
plot(time, CPO, pch = ifelse(status == 1, 16, 1),
     col = ifelse(status == 1, "darkblue", "orange"),
     xlab = "Survival Time (months)", ylab = "CPO",
     main = "CPO vs. Survival Time",
     cex = 0.8, log = "y")
legend("topright", legend = c("Observed", "Censored"),
       pch = c(16, 1), col = c("darkblue", "orange"), bty = "n")

# Panel 2: log(CPO) vs survival time
plot(time, log_CPO, pch = ifelse(status == 1, 16, 1),
     col = ifelse(status == 1, "darkblue", "orange"),
     xlab = "Survival Time (months)", ylab = "log(CPO)",
     main = "log(CPO) vs. Survival Time",
     cex = 0.8)
legend("topright", legend = c("Observed", "Censored"),
       pch = c(16, 1), col = c("darkblue", "orange"), bty = "n")

# Identify potential outliers (low CPO)
threshold <- quantile(log_CPO, 0.05)
outliers <- which(log_CPO < threshold)
if (length(outliers) > 0) {
  points(time[outliers], log_CPO[outliers], pch = 4, col = "red", cex = 1.5)
}

# Panel 3: Histogram of log(CPO)
hist(log_CPO, breaks = 30, col = "lightblue", border = "white",
     main = "Distribution of log(CPO)",
     xlab = "log(CPO)", ylab = "Frequency")
abline(v = mean(log_CPO), col = "red", lwd = 2, lty = 2)
abline(v = median(log_CPO), col = "blue", lwd = 2, lty = 2)
legend("topleft", legend = c(paste("Mean =", round(mean(log_CPO), 2)),
                              paste("Median =", round(median(log_CPO), 2))),
       col = c("red", "blue"), lwd = 2, lty = 2, bty = "n")

dev.off()
cat("   Saved: aa_d_cpo_diagnostic.pdf\n\n")

# ==============================================================================
# 6. COMBINED DIAGNOSTIC FIGURE
# ==============================================================================

cat("6. Creating combined diagnostic figure...\n")

pdf("aa_d_model_diagnostics_combined.pdf", width = 14, height = 10)
par(mfrow = c(2, 3), mar = c(5, 5, 4, 2))

# Panel 1: Nelson-Aalen diagnostic
plot(na_time, na_cumhaz, type = "s", lwd = 2, col = "blue",
     xlab = "Cox-Snell Residuals", ylab = "Nelson-Aalen Cumulative Hazard",
     main = "(A) Nelson-Aalen Diagnostic",
     xlim = c(0, max(na_time) * 1.1),
     ylim = c(0, max(na_cumhaz) * 1.1))
abline(0, 1, col = "red", lwd = 2, lty = 2)
legend("topleft", legend = c("Observed", "Theoretical Exp(1)"),
       col = c("blue", "red"), lwd = 2, lty = c(1, 2), bty = "n", cex = 0.9)

# Panel 2: Q-Q plot
plot(theoretical_quantiles, empirical_quantiles, pch = 16, cex = 0.6,
     xlab = "Theoretical Exp(1) Quantiles",
     ylab = "Empirical Quantiles",
     main = "(B) Q-Q Plot (Uncensored)")
abline(0, 1, col = "red", lwd = 2)

# Panel 3: Cox-Snell residuals by censoring
plot(1:N, H_residuals, pch = ifelse(status == 1, 16, 1),
     col = ifelse(status == 1, "darkgreen", "red"),
     xlab = "Observation Index", ylab = "Cox-Snell Residual",
     main = "(C) Cox-Snell Residuals",
     cex = 0.7)
abline(h = 1, col = "gray", lty = 2)
legend("topright", legend = c("Observed", "Censored"),
       pch = c(16, 1), col = c("darkgreen", "red"), bty = "n", cex = 0.9)

# Panel 4: CPO vs time
plot(time, log_CPO, pch = ifelse(status == 1, 16, 1),
     col = ifelse(status == 1, "darkblue", "orange"),
     xlab = "Survival Time (months)", ylab = "log(CPO)",
     main = "(D) log(CPO) vs. Survival Time",
     cex = 0.7)
legend("topright", legend = c("Observed", "Censored"),
       pch = c(16, 1), col = c("darkblue", "orange"), bty = "n", cex = 0.9)

# Panel 5: Histogram of log(CPO)
hist(log_CPO, breaks = 30, col = "lightblue", border = "white",
     main = "(E) Distribution of log(CPO)",
     xlab = "log(CPO)", ylab = "Frequency")
abline(v = median(log_CPO), col = "red", lwd = 2, lty = 2)

# Panel 6: Model fit summary text
plot.new()
plot.window(xlim = c(0, 1), ylim = c(0, 1))
title(main = "(F) Model Fit Summary")

text(0.1, 0.95, "AA-D Stratum Analysis", adj = 0, font = 2, cex = 1.1)
text(0.1, 0.85, paste("N =", N, "patients,", K_data, "counties"), adj = 0)
text(0.1, 0.75, paste("Censoring rate:", round(100 * mean(status == 0), 1), "%"), adj = 0)

text(0.1, 0.60, "Model Parameters:", adj = 0, font = 2)
text(0.1, 0.50, paste("  Spatial rho =", round(mean(fit$out_fit$rho_samples), 3)), adj = 0)
text(0.1, 0.42, paste("  Error sigma =", round(mean(sqrt(fit$out_fit$sigma2_samples)), 3)), adj = 0)
text(0.1, 0.34, paste("  MCMC samples =", n_samples), adj = 0)

text(0.1, 0.20, "CPO Diagnostics:", adj = 0, font = 2)
text(0.1, 0.10, paste("  LPML =", round(LPML, 2)), adj = 0)
text(0.1, 0.02, paste("  Mean log(CPO) =", round(mean(log_CPO), 2)), adj = 0)

dev.off()
cat("   Saved: aa_d_model_diagnostics_combined.pdf\n\n")

# ==============================================================================
# 7. SAVE RESULTS
# ==============================================================================

cat("7. Saving results...\n")

# Save diagnostic data
diagnostics <- list(
  # Data info
  N = N,
  K = K_data,
  time = time,
  status = status,
  censoring_rate = mean(status == 0),

  # Model fit
  fit = fit,
  sigma_mean = sigma_mean,
  rho_mean = mean(fit$out_fit$rho_samples),

  # Cox-Snell residuals
  H_residuals = H_residuals,
  mu_hat = mu_hat,

  # Nelson-Aalen
  na_time = na_time,
  na_cumhaz = na_cumhaz,

  # CPO
  CPO = CPO,
  log_CPO = log_CPO,
  LPML = LPML
)

save(diagnostics, file = "aa_d_diagnostics.RData")
cat("   Saved: aa_d_diagnostics.RData\n\n")

# ==============================================================================
# 8. SUMMARY
# ==============================================================================

cat("==============================================================================\n")
cat("ANALYSIS COMPLETE\n")
cat("==============================================================================\n\n")

cat("Output files:\n")
cat("  - aa_d_na_km_diagnostic.pdf        : Nelson-Aalen/KM diagnostic plots\n")
cat("  - aa_d_cpo_diagnostic.pdf          : CPO diagnostic plots\n")
cat("  - aa_d_model_diagnostics_combined.pdf : Combined diagnostic figure\n")
cat("  - aa_d_diagnostics.RData           : All diagnostic data\n\n")

cat("Diagnostic Interpretation:\n")
cat("  1. Nelson-Aalen Plot: If the cumulative hazard of Cox-Snell residuals\n")
cat("     follows the 45-degree line, the log-normal assumption is reasonable.\n\n")

cat("  2. CPO Diagnostics:\n")
cat("     - LPML =", round(LPML, 2), "(higher is better)\n")
cat("     - Observations with very low CPO may be outliers or poorly fit\n\n")

cat("  3. For causal inference using RMST (restricted mean), the tail behavior\n")
cat("     of the log-normal is less critical since we integrate up to t*.\n\n")

cat("==============================================================================\n")
