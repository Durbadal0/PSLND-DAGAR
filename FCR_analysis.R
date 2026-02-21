# ==============================================================================
# Florida Cancer Registry (FCR) Data Analysis
# ==============================================================================
#
# This script replicates the data analysis from Section 4 of the paper:
# "Estimation of cluster-specific causal effects on spatially associated
#  survival data using SoftBART"
#
# Analysis focuses on WA-D patients (White/Non-AA, Distant stage)
# N = 1994 patients across 61 Florida counties
#
# Outputs:
#   1. Patient-level CERM estimates (Table 2 in paper)
#   2. County-level ACERM, ACERMT, ACERMU maps (Figures 1-2)
#   3. State-level RATE, RATT, RATU estimates (Table 3)
#   4. SPTE survival curves (Figure 3)
#
# ==============================================================================

rm(list = ls())
set.seed(12345)

# Load implementation
source("ps_lnd_dagar.R")

cat("==============================================================================\n")
cat("Florida Cancer Registry (FCR) Data Analysis\n")
cat("Causal Effects of Treatment Delay on Breast Cancer Survival\n")
cat("==============================================================================\n\n")

# ==============================================================================
# 1. LOAD AND PREPARE DATA
# ==============================================================================

cat("1. Loading and preparing data...\n")

# Load FCR data
data <- read.csv("Surv_data2.csv")

# Load Florida county adjacency matrix
fl_adj <- as.matrix(read.csv("W.mat.csv", row.names = 1))

# Filter to WA-D subset (White/Non-AA, Distant stage)
# Race = 1 (White), Stage = 3 (Distant)
wa_d <- data[data$Race == 1 & data$Stage == 3, ]

cat("  Total FCR patients:", nrow(data), "\n")
cat("  WA-D subset: N =", nrow(wa_d), "\n")
cat("  Counties with WA-D patients:", length(unique(wa_d$county)), "\n")

# Prepare variables
time <- wa_d$as.numeric.date_diff. / 30.44  # Convert days to months
status <- wa_d$death                 # 1 = death, 0 = censored
Z <- ifelse(wa_d$TX_Delay == 1, 1, 0)  # 1 = short TD (<90d), 0 = long TD (>90d)
county <- wa_d$county

# Covariates: Age, BX_Delay, HR_p, Tgrade
# BX_Delay: 1 = short (no delay), 3 = long (delay) -> convert to 0/1
# HR_p: 1 = positive, 0 = negative
# Tgrade: 1, 2, 3
Age <- wa_d$Age
BD <- ifelse(wa_d$BX_Delay == 3, 1, 0)  # 1 = biopsy delay, 0 = no delay
HR <- wa_d$HR_p                          # 1 = positive, 0 = negative
TG <- wa_d$Tgrade                        # 1, 2, 3

X <- cbind(Age, BD, HR, TG)
colnames(X) <- c("Age", "BD", "HR", "TG")

cat("\n  Censoring rate:", round(100 * mean(status == 0), 1), "%\n")
cat("  Short TD (Z=1):", sum(Z == 1), "\n")
cat("  Long TD (Z=0):", sum(Z == 0), "\n")
cat("  Median survival time:", round(median(time), 1), "months\n")

# Get unique counties and create mapping
unique_counties <- sort(unique(county))
K_data <- length(unique_counties)
cat("  Counties in data:", K_data, "out of 67\n\n")

# ==============================================================================
# 2. FIT PS-LND-DAGAR MODEL
# ==============================================================================

cat("2. Fitting PS-LND-DAGAR model...\n")
cat("   (This may take several minutes)\n\n")

# MCMC settings - use more iterations for final analysis
# For testing, reduce these values
n_mcmc_ps <- 2000    # Propensity score MCMC iterations
burn_in_ps <- 500    # PS burn-in
n_mcmc_out <- 3000   # Outcome model MCMC iterations
burn_in_out <- 1000  # Outcome burn-in
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
# 3. PATIENT-LEVEL CERM ESTIMATES (Table 2)
# ==============================================================================

cat("3. Computing patient-level CERM estimates (Table 2)...\n")

# Focus on Broward County (i = 6) as in the paper
broward_id <- 6
broward_idx <- which(county == broward_id)
cat("   Broward County (i=6): n =", length(broward_idx), "patients\n")

# Create grid of covariate combinations
# Age: 55, 65, 75
# BD: 0 (No), 1 (Yes)
# HR: 0 (Negative), 1 (Positive)
# TG: 1, 2, 3
ages <- c(55, 65, 75)
bds <- c(0, 1)
hrs <- c(0, 1)
tgs <- c(1, 2, 3)

# Create all combinations
cov_grid <- expand.grid(Age = ages, BD = bds, HR = hrs, TG = tgs)
cov_grid <- cov_grid[order(cov_grid$BD, cov_grid$HR, cov_grid$TG, cov_grid$Age), ]

# Get mean propensity score for Broward (use as representative)
e_hat_broward <- mean(fit$e_hat[broward_idx])

# Compute CERM for each covariate combination
cat("   Computing CERM for 36 covariate combinations...\n")

cerm_results <- data.frame(
  BD = character(),
  HR = character(),
  TG = integer(),
  Age = integer(),
  CERM_5yr = numeric(),
  CERM_5yr_lower = numeric(),
  CERM_5yr_upper = numeric(),
  CERM_10yr = numeric(),
  CERM_10yr_lower = numeric(),
  CERM_10yr_upper = numeric(),
  stringsAsFactors = FALSE
)

for (i in 1:nrow(cov_grid)) {
  x <- as.numeric(cov_grid[i, ])

  # CERM at 5 years (60 months)
  cerm_5 <- estimate_cerm(t_star = 60, x = x, county_id = broward_id,
                          e_hat = e_hat_broward, fit = fit$out_fit,
                          orig_to_internal = fit$orig_to_internal)

  # CERM at 10 years (120 months)
  cerm_10 <- estimate_cerm(t_star = 120, x = x, county_id = broward_id,
                           e_hat = e_hat_broward, fit = fit$out_fit,
                           orig_to_internal = fit$orig_to_internal)

  cerm_results <- rbind(cerm_results, data.frame(
    BD = ifelse(cov_grid$BD[i] == 0, "No", "Yes"),
    HR = ifelse(cov_grid$HR[i] == 0, "Negative", "Positive"),
    TG = cov_grid$TG[i],
    Age = cov_grid$Age[i],
    CERM_5yr = cerm_5$mean,
    CERM_5yr_lower = cerm_5$lower,
    CERM_5yr_upper = cerm_5$upper,
    CERM_10yr = cerm_10$mean,
    CERM_10yr_lower = cerm_10$lower,
    CERM_10yr_upper = cerm_10$upper,
    stringsAsFactors = FALSE
  ))
}

cat("\n   Table 2: Patient-level CERM estimates (months) for Broward County\n")
cat("   ====================================================================\n")
print(cerm_results, row.names = FALSE)

# Save results
write.csv(cerm_results, "table2_cerm_broward.csv", row.names = FALSE)
cat("\n   Saved to: table2_cerm_broward.csv\n\n")

# ==============================================================================
# 4. COUNTY-LEVEL TREATMENT EFFECTS
# ==============================================================================

cat("4. Computing county-level treatment effects...\n")

# Initialize storage for all 67 counties
county_results <- data.frame(
  county_id = 1:67,
  n_patients = integer(67),
  n_treated = integer(67),
  n_untreated = integer(67),
  ACERM_5yr = numeric(67),
  ACERM_5yr_lower = numeric(67),
  ACERM_5yr_upper = numeric(67),
  ACERM_10yr = numeric(67),
  ACERM_10yr_lower = numeric(67),
  ACERM_10yr_upper = numeric(67),
  ACERMT_5yr = numeric(67),
  ACERMT_5yr_lower = numeric(67),
  ACERMT_5yr_upper = numeric(67),
  ACERMT_10yr = numeric(67),
  ACERMT_10yr_lower = numeric(67),
  ACERMT_10yr_upper = numeric(67),
  ACERMU_5yr = numeric(67),
  ACERMU_5yr_lower = numeric(67),
  ACERMU_5yr_upper = numeric(67),
  ACERMU_10yr = numeric(67),
  ACERMU_10yr_lower = numeric(67),
  ACERMU_10yr_upper = numeric(67),
  stringsAsFactors = FALSE
)

# Set NA for counties without patients
county_results[, 5:22] <- NA

n_draws <- fit$out_fit$n_save

for (cid in unique_counties) {
  cat("   Processing county", cid, "...\n")

  idx <- which(county == cid)
  n_county <- length(idx)
  X_county <- X[idx, , drop = FALSE]
  Z_county <- Z[idx]
  e_hat_county <- fit$e_hat[idx]

  county_results$n_patients[cid] <- n_county
  county_results$n_treated[cid] <- sum(Z_county == 1)
  county_results$n_untreated[cid] <- sum(Z_county == 0)

  # Compute individual CERM for all patients in county
  CERM_5yr_draws <- matrix(0, nrow = n_draws, ncol = n_county)
  CERM_10yr_draws <- matrix(0, nrow = n_draws, ncol = n_county)

  for (j in 1:n_county) {
    cerm_5 <- estimate_cerm(60, X_county[j, ], cid, e_hat_county[j], fit$out_fit,
                             fit$orig_to_internal)
    cerm_10 <- estimate_cerm(120, X_county[j, ], cid, e_hat_county[j], fit$out_fit,
                              fit$orig_to_internal)
    CERM_5yr_draws[, j] <- cerm_5$draws
    CERM_10yr_draws[, j] <- cerm_10$draws
  }

  # ACERM: Average over all patients
  ACERM_5yr_draws <- rowMeans(CERM_5yr_draws)
  ACERM_10yr_draws <- rowMeans(CERM_10yr_draws)

  county_results$ACERM_5yr[cid] <- mean(ACERM_5yr_draws)
  county_results$ACERM_5yr_lower[cid] <- quantile(ACERM_5yr_draws, 0.025)
  county_results$ACERM_5yr_upper[cid] <- quantile(ACERM_5yr_draws, 0.975)
  county_results$ACERM_10yr[cid] <- mean(ACERM_10yr_draws)
  county_results$ACERM_10yr_lower[cid] <- quantile(ACERM_10yr_draws, 0.025)
  county_results$ACERM_10yr_upper[cid] <- quantile(ACERM_10yr_draws, 0.975)

  # ACERMT: Average over treated patients only
  if (sum(Z_county == 1) > 0) {
    treated_idx <- which(Z_county == 1)
    ACERMT_5yr_draws <- rowMeans(CERM_5yr_draws[, treated_idx, drop = FALSE])
    ACERMT_10yr_draws <- rowMeans(CERM_10yr_draws[, treated_idx, drop = FALSE])

    county_results$ACERMT_5yr[cid] <- mean(ACERMT_5yr_draws)
    county_results$ACERMT_5yr_lower[cid] <- quantile(ACERMT_5yr_draws, 0.025)
    county_results$ACERMT_5yr_upper[cid] <- quantile(ACERMT_5yr_draws, 0.975)
    county_results$ACERMT_10yr[cid] <- mean(ACERMT_10yr_draws)
    county_results$ACERMT_10yr_lower[cid] <- quantile(ACERMT_10yr_draws, 0.025)
    county_results$ACERMT_10yr_upper[cid] <- quantile(ACERMT_10yr_draws, 0.975)
  }

  # ACERMU: Average over untreated patients only
  if (sum(Z_county == 0) > 0) {
    untreated_idx <- which(Z_county == 0)
    ACERMU_5yr_draws <- rowMeans(CERM_5yr_draws[, untreated_idx, drop = FALSE])
    ACERMU_10yr_draws <- rowMeans(CERM_10yr_draws[, untreated_idx, drop = FALSE])

    county_results$ACERMU_5yr[cid] <- mean(ACERMU_5yr_draws)
    county_results$ACERMU_5yr_lower[cid] <- quantile(ACERMU_5yr_draws, 0.025)
    county_results$ACERMU_5yr_upper[cid] <- quantile(ACERMU_5yr_draws, 0.975)
    county_results$ACERMU_10yr[cid] <- mean(ACERMU_10yr_draws)
    county_results$ACERMU_10yr_lower[cid] <- quantile(ACERMU_10yr_draws, 0.025)
    county_results$ACERMU_10yr_upper[cid] <- quantile(ACERMU_10yr_draws, 0.975)
  }
}

# Save county results
write.csv(county_results, "county_effects.csv", row.names = FALSE)
cat("\n   Saved to: county_effects.csv\n")

# Print summary for key counties mentioned in paper
cat("\n   Key county results (ACERM in months):\n")
cat("   ----------------------------------------\n")
key_counties <- c(2, 4, 6, 21, 63)  # Baker, Bradford, Broward, Gilchrist, Union
for (cid in key_counties) {
  if (!is.na(county_results$ACERM_10yr[cid])) {
    cat(sprintf("   County %2d: ACERM(10yr) = %5.1f [%5.1f, %5.1f], n = %d\n",
                cid,
                county_results$ACERM_10yr[cid],
                county_results$ACERM_10yr_lower[cid],
                county_results$ACERM_10yr_upper[cid],
                county_results$n_patients[cid]))
  }
}
cat("\n")

# ==============================================================================
# 5. STATE-LEVEL TREATMENT EFFECTS (Table 3)
# ==============================================================================

cat("5. Computing state-level treatment effects (Table 3)...\n")

# Compute CERM for all patients in the state
N <- length(time)
N_1 <- sum(Z == 1)  # Treated
N_0 <- sum(Z == 0)  # Untreated

cat("   Computing CERM for all", N, "patients...\n")

CERM_5yr_all <- matrix(0, nrow = n_draws, ncol = N)
CERM_10yr_all <- matrix(0, nrow = n_draws, ncol = N)

for (i in 1:N) {
  if (i %% 200 == 0) cat("     Patient", i, "/", N, "\n")
  cerm_5 <- estimate_cerm(60, X[i, ], county[i], fit$e_hat[i], fit$out_fit,
                           fit$orig_to_internal)
  cerm_10 <- estimate_cerm(120, X[i, ], county[i], fit$e_hat[i], fit$out_fit,
                            fit$orig_to_internal)
  CERM_5yr_all[, i] <- cerm_5$draws
  CERM_10yr_all[, i] <- cerm_10$draws
}

# RATE: Average over all patients
RATE_5yr_draws <- rowMeans(CERM_5yr_all)
RATE_10yr_draws <- rowMeans(CERM_10yr_all)

# RATT: Average over treated patients
treated_idx <- which(Z == 1)
RATT_5yr_draws <- rowMeans(CERM_5yr_all[, treated_idx, drop = FALSE])
RATT_10yr_draws <- rowMeans(CERM_10yr_all[, treated_idx, drop = FALSE])

# RATU: Average over untreated patients
untreated_idx <- which(Z == 0)
RATU_5yr_draws <- rowMeans(CERM_5yr_all[, untreated_idx, drop = FALSE])
RATU_10yr_draws <- rowMeans(CERM_10yr_all[, untreated_idx, drop = FALSE])

# Create Table 3
table3 <- data.frame(
  Estimand = c("RATE", "RATT", "RATU", "RATE", "RATT", "RATU"),
  Interval_years = c(10, 10, 10, 5, 5, 5),
  Sample_Size = c(N, N_1, N_0, N, N_1, N_0),
  Estimate_months = c(
    mean(RATE_10yr_draws), mean(RATT_10yr_draws), mean(RATU_10yr_draws),
    mean(RATE_5yr_draws), mean(RATT_5yr_draws), mean(RATU_5yr_draws)
  ),
  CI_lower = c(
    quantile(RATE_10yr_draws, 0.025), quantile(RATT_10yr_draws, 0.025),
    quantile(RATU_10yr_draws, 0.025),
    quantile(RATE_5yr_draws, 0.025), quantile(RATT_5yr_draws, 0.025),
    quantile(RATU_5yr_draws, 0.025)
  ),
  CI_upper = c(
    quantile(RATE_10yr_draws, 0.975), quantile(RATT_10yr_draws, 0.975),
    quantile(RATU_10yr_draws, 0.975),
    quantile(RATE_5yr_draws, 0.975), quantile(RATT_5yr_draws, 0.975),
    quantile(RATU_5yr_draws, 0.975)
  ),
  stringsAsFactors = FALSE
)

cat("\n   Table 3: State-level treatment effects (months)\n")
cat("   ================================================\n")
print(table3, row.names = FALSE)

write.csv(table3, "table3_state_effects.csv", row.names = FALSE)
cat("\n   Saved to: table3_state_effects.csv\n\n")

# ==============================================================================
# 6. SPTE SURVIVAL CURVES (Figure 3)
# ==============================================================================

cat("6. Computing SPTE survival curves (Figure 3)...\n")

# Time points for survival curves (0 to 120 months)
t_points <- seq(0, 120, by = 6)
n_t <- length(t_points)

# Storage for survival probabilities
S1_draws <- matrix(0, nrow = n_draws, ncol = n_t)
S0_draws <- matrix(0, nrow = n_draws, ncol = n_t)

# For each MCMC draw, compute average survival probabilities
cat("   Computing survival curves at", n_t, "time points...\n")

for (m in 1:n_draws) {
  if (m %% 50 == 0) cat("     MCMC draw", m, "/", n_draws, "\n")

  sigma_m <- sqrt(fit$out_fit$sigma2_samples[m])

  for (t_idx in 1:n_t) {
    t <- t_points[t_idx]
    if (t == 0) {
      S1_draws[m, t_idx] <- 1
      S0_draws[m, t_idx] <- 1
    } else {
      S1_sum <- 0
      S0_sum <- 0

      for (i in 1:N) {
        # Get b2 predictions
        x_out_1 <- scale_new(matrix(c(X[i, ], 1, fit$e_hat[i]), nrow = 1),
                             fit$out_fit$scaling_params)
        x_out_0 <- scale_new(matrix(c(X[i, ], 0, fit$e_hat[i]), nrow = 1),
                             fit$out_fit$scaling_params)

        b2_1 <- fit$out_fit$forest$predict_iteration(x_out_1, iter = m)
        b2_0 <- fit$out_fit$forest$predict_iteration(x_out_0, iter = m)

        county_idx <- fit$orig_to_internal[as.character(county[i])]
        W_m <- fit$out_fit$w_samples[m, county_idx]

        mu_1 <- b2_1 + W_m
        mu_0 <- b2_0 + W_m

        # S(t) = 1 - Phi((log(t) - mu) / sigma)
        S1_sum <- S1_sum + (1 - pnorm((log(t) - mu_1) / sigma_m))
        S0_sum <- S0_sum + (1 - pnorm((log(t) - mu_0) / sigma_m))
      }

      S1_draws[m, t_idx] <- S1_sum / N
      S0_draws[m, t_idx] <- S0_sum / N
    }
  }
}

# Compute SPTE = S1 - S0
SPTE_draws <- S1_draws - S0_draws

# Summary statistics
S1_mean <- colMeans(S1_draws)
S0_mean <- colMeans(S0_draws)
SPTE_mean <- colMeans(SPTE_draws)
SPTE_lower <- apply(SPTE_draws, 2, quantile, 0.025)
SPTE_upper <- apply(SPTE_draws, 2, quantile, 0.975)

# Create survival curves data frame
survival_curves <- data.frame(
  time_months = t_points,
  S1 = S1_mean,
  S0 = S0_mean,
  SPTE = SPTE_mean,
  SPTE_lower = SPTE_lower,
  SPTE_upper = SPTE_upper
)

write.csv(survival_curves, "survival_curves.csv", row.names = FALSE)
cat("   Saved to: survival_curves.csv\n")

# Find peak SPTE
peak_idx <- which.max(SPTE_mean)
cat("\n   Peak SPTE:", round(SPTE_mean[peak_idx] * 100, 1), "% at",
    t_points[peak_idx], "months\n")
cat("   SPTE at 120 months:", round(SPTE_mean[n_t] * 100, 1), "%\n\n")

# ==============================================================================
# 7. CREATE PLOTS
# ==============================================================================

cat("7. Creating plots...\n")

# Plot 1: SPTE survival curves (Figure 3)
pdf("figure3_spte.pdf", width = 10, height = 6)
par(mfrow = c(1, 2), mar = c(5, 4, 3, 2))

# Left panel: S(t) curves
plot(t_points, S1_mean, type = "l", lty = 2, lwd = 2,
     xlab = "Time (months)", ylab = "Survival Probability",
     main = "Survival Curves", ylim = c(0, 1))
lines(t_points, S0_mean, lty = 3, lwd = 2)
legend("topright", legend = c("S(1)(t) - Short TD", "S(0)(t) - Long TD"),
       lty = c(2, 3), lwd = 2, bty = "n")

# Right panel: SPTE
plot(t_points, SPTE_mean, type = "l", lwd = 2, col = "blue",
     xlab = "Time (months)", ylab = "SPTE = S(1)(t) - S(0)(t)",
     main = "Survival Probability Treatment Effect",
     ylim = c(min(SPTE_lower) - 0.02, max(SPTE_upper) + 0.02))
polygon(c(t_points, rev(t_points)),
        c(SPTE_lower, rev(SPTE_upper)),
        col = rgb(0, 0, 1, 0.2), border = NA)
lines(t_points, SPTE_mean, lwd = 2, col = "blue")
abline(h = 0, lty = 2, col = "gray")

dev.off()
cat("   Saved: figure3_spte.pdf\n")

# Plot 2: County ACERM map (simplified - actual maps would need shapefile)
pdf("figure1_acerm_map.pdf", width = 12, height = 6)
par(mfrow = c(1, 2), mar = c(5, 4, 3, 2))

# 5-year ACERM by county
valid_counties <- which(!is.na(county_results$ACERM_5yr))
barplot(county_results$ACERM_5yr[valid_counties],
        names.arg = valid_counties,
        main = "ACERM at 5 years by County",
        xlab = "County ID", ylab = "ACERM (months)",
        col = "steelblue", las = 2, cex.names = 0.6)

# 10-year ACERM by county
barplot(county_results$ACERM_10yr[valid_counties],
        names.arg = valid_counties,
        main = "ACERM at 10 years by County",
        xlab = "County ID", ylab = "ACERM (months)",
        col = "darkgreen", las = 2, cex.names = 0.6)

dev.off()
cat("   Saved: figure1_acerm_map.pdf\n")

# ==============================================================================
# 8. SUMMARY
# ==============================================================================

cat("\n")
cat("==============================================================================\n")
cat("ANALYSIS COMPLETE\n")
cat("==============================================================================\n")
cat("\nOutput files:\n")
cat("  - table2_cerm_broward.csv    : Patient-level CERM for Broward County\n")
cat("  - county_effects.csv         : County-level ACERM, ACERMT, ACERMU\n")
cat("  - table3_state_effects.csv   : State-level RATE, RATT, RATU\n")
cat("  - survival_curves.csv        : SPTE survival curves data\n")
cat("  - figure3_spte.pdf           : SPTE survival curves plot\n")
cat("  - figure1_acerm_map.pdf      : County ACERM bar plots\n")
cat("\nModel diagnostics:\n")
cat("  - Spatial rho:", round(mean(fit$out_fit$rho_samples), 3), "\n")
cat("  - Error sigma:", round(mean(sqrt(fit$out_fit$sigma2_samples)), 3), "\n")
cat("  - PS equicorrelation rho_0:", round(mean(fit$ps_fit$rho_0_samples), 3), "\n")
cat("\n")

# Save full model fit for later use
save(fit, county_results, table3, survival_curves,
     file = "fcr_analysis_results.RData")
cat("Full results saved to: fcr_analysis_results.RData\n")
cat("==============================================================================\n")
