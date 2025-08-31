# ------------------------------------------------------------------------------
# 0.  Load required libraries + source Functions.R
# ------------------------------------------------------------------------------
library(dplyr)
library(readr)
library(SoftBart)
library(truncnorm)
library(MASS)
library(ggplot2)
library(sf)
library(tigris)

# Source your Functions.R (from Functions.docx)
source("Functions.R")

# ------------------------------------------------------------------------------
# 1.  Load FCR data and preprocess
# ------------------------------------------------------------------------------
fcr_data <- read_csv("Surv_data2.csv")


fcr_modified <- fcr_data %>%
  mutate(
    time   = fcr_data$`as.numeric.date_diff.`,
    status = death,
    county = county,
    Age    = Age,
    BX_Delay = BX_Delay,
    HR_p   = HR_p,
    Tgrade = Tgrade,

    Race_cat  = Race,
    Stage_cat = factor(Stage, levels = c(1,2,3), labels = c("Local","Regional","Distant")),
   =
    Z = ifelse(TX_Delay == 1, 0, 1)
  ) %>%
  filter(!is.na(Age), !is.na(BX_Delay), !is.na(HR_p), !is.na(Tgrade), !is.na(Z))

# (1b) Load adjacency W (67×67) for DAGAR
W_mat <- read_csv("W.mat.csv", col_types = cols(.default = col_double()))
W_mat <- as.matrix(W_mat); storage.mode(W_mat) <- "numeric"

# ------------------------------------------------------------------------------
# 2.  Subset to one stratum: race = 1, stage = “Distant”
# ------------------------------------------------------------------------------

race_val  <- 1
stage_val <- "Distant"

data_stratum <- fcr_modified %>%
  filter(Race_cat == race_val, Stage_cat == stage_val)

cat("Observations in that stratum:", nrow(data_stratum), "\n")
stratum_counties <- sort(unique(data_stratum$county))
cat("Number of counties in stratum:", length(stratum_counties), "\n")

K_stratum <- length(stratum_counties)
W_subset <- matrix(0, nrow = K_stratum, ncol = K_stratum)
for(i in seq_along(stratum_counties)) {
  for(j in seq_along(stratum_counties)) {
    ci <- stratum_counties[i]
    cj <- stratum_counties[j]
    if(ci <= ncol(W_mat) && cj <= nrow(W_mat)) {
      W_subset[i,j] <- W_mat[ci, cj]
    }
  }
}

W_subset <- (W_subset + t(W_subset)) / 2
diag(W_subset) <- 0

# ------------------------------------------------------------------------------
n_iter_ps <- 120000  # # iterations for propensity SBART
n_mcmc   <- 40000
burn_in  <- 10
thin     <- 5

)
rho_init    <- 0.5
rho_prop_sd <- 0.02
a_tau       <- 1; b_tau <- 1
a_sigma     <- 2; b_sigma <- 1

cat("Fitting AFT_mixed_DAGAR_probit_causal() ...\n")
fit_fcr2 <- AFT_mixed_DAGAR_probit_causal2(
  time        = time,
  status      = status,
  X           = X,
  V           = V,
  Z           = Z,
  group       = group,
  X_test      = X,      # we predict on training points
  V_test      = V,
  Z_test      = Z,
  group_test  = group,
  num_tree_ps = 50,
  num_tree_out= 50,
  n_iter_ps   = n_iter_ps,
  n_mcmc      = n_mcmc,
  burn_in     = burn_in,
  thin        = thin,
  rho_init    = rho_init,
  rho_prop_sd = rho_prop_sd,
  a_tau       = a_tau,
  b_tau       = b_tau,
  a_sigma     = a_sigma,
  b_sigma     = b_sigma
)
cat("Model fitting complete.\n")
fit=fit_fcr2
# Extract MCMC outputs
forest_ps   <- fit$probit_forest
forest_out  <- fit$forest_out
w_samples   <- fit$w_samples       # (M × K_stratum) matrix of DAGAR draws
sigma2_samps<- fit$sigma2_samps    # length‐M vector
# Note: M = floor((n_mcmc - burn_in)/thin)

M_total <- length(sigma2_samps)
cat("Number of posterior draws used (M):", M_total, "\n")





forest_ps    <- fit$probit_forest
forest_out   <- fit$forest_out
w_samples    <- fit$w_samples      # M × K matrix
sigma2_samps <- fit$sigma2_samps   # length-M vector

burn_in   <- burn_in   # e.g. 500
thin      <- thin      # e.g. 5
n_mcmc    <- n_mcmc    # e.g. 10000
n_iter_ps <- n_iter_ps # e.g. 15000

# Choose t_star = 5 (in same time‐units as your 'time' variable;
# if time is days, use 5*365)
t_star <- 5*365

# ----- 2. Prepare data dimensions -----
n        <- nrow(X)       # number of patients
sample_iters_out <- seq(burn_in + 1, n_mcmc, by = thin)
Mprime   <- length(sample_iters_out)

# ----- 3. Pre‐allocate matrices to hold all posterior draws -----
CATE_draws_matrix  <- matrix(NA, nrow = Mprime, ncol = n)
CRATE_draws_matrix <- matrix(NA, nrow = Mprime, ncol = n)

# ----- 4. Loop over patients to get draws -----
for (j in seq_len(n)) {
  # patient‐level covariates and county info:
  x_j       <- X[j, ]
  v_j       <- V[j, ]        # if V is a matrix of intercept+features
  county_j  <- group[j]

  # (a) CATE draws for patient j
  res_cate <- estimate_CATE_corrected(
    x            = x_j,
    v_i          = v_j,
    county_id    = county_j,
    forest_ps    = forest_ps,
    forest_out   = forest_out,
    w_samples    = w_samples,
    sigma2_samps = sigma2_samps,
    burn_in      = burn_in,
    thin         = thin,
    n_mcmc       = n_mcmc,
    n_iter_ps    = n_iter_ps
  )
  CATE_draws_matrix[, j] <- res_cate$CATE_draws

  # (b) CRATE@5 draws for patient j
  res_crate <- estimate_CRATE(
    t_star       = t_star,
    x            = x_j,
    v_i          = v_j,
    county_id    = county_j,
    forest_ps    = forest_ps,
    forest_out   = forest_out,
    w_samples    = w_samples,
    sigma2_samps = sigma2_samps,
    burn_in      = burn_in,
    thin         = thin,
    n_mcmc       = n_mcmc,
    n_iter_ps    = n_iter_ps
  )
  CRATE_draws_matrix[, j] <- res_crate$CRATE_draws
}

# ----- 5. Summaries: posterior mean & 95% CI for each patient -----
# CATE
CATE_mean  <- colMeans(CATE_draws_matrix)
CATE_ci    <- apply(CATE_draws_matrix, 2, quantile, probs = c(0.025, 0.975))

# CRATE@5
CRATE_mean <- colMeans(CRATE_draws_matrix)
CRATE_ci   <- apply(CRATE_draws_matrix, 2, quantile, probs = c(0.025, 0.975))

# ----- 6. Combine into a results data.frame -----
results <- data.frame(
  patient   = seq_len(n),
  county    = group,

  CATE_mean   = CATE_mean,
  CATE_lower  = CATE_ci[1, ],
  CATE_upper  = CATE_ci[2, ],

  CRATE_mean  = CRATE_mean,
  CRATE_lower = CRATE_ci[1, ],
  CRATE_upper = CRATE_ci[2, ]
)

# Inspect the first few rows
head(results)




write.csv(
  CATE_draws_matrix,
  file      = "CATE_draws_matrix_WD_2.csv",
  row.names = FALSE
)

write.csv(
  CRATE_draws_matrix,
  file      = "CRATE_draws_matrix_WD_2.csv",
  row.names = FALSE
)

# 2. (Optional) Save the patient‐level summary
write.csv(
  results,
  file      = "CATE_CRATE_summary_WD_2.csv",
  row.names = FALSE
)

message("Saved CATE draws to 'CATE_draws_matrix_WD.csv', CRATE draws to 'CRATE_draws_matrix_WD.csv', and summaries to 'CATE_CRATE_summary_WD.csv'.")






age_min   <- min(data_stratum$Age)
age_max   <- max(data_stratum$Age)

bx_min    <- min(data_stratum$BX_Delay)
bx_max    <- max(data_stratum$BX_Delay)

hr_min    <- min(data_stratum$HR_p)
hr_max    <- max(data_stratum$HR_p)

grade_min <- min(data_stratum$Tgrade)
grade_max <- max(data_stratum$Tgrade)


cov_raw <- expand.grid(
  age          = c(50, 60, 74),
  biopsy_delay = c(0, 1),
  hs_status    = c(0, 1),
  tumor_grade  = c(1, 2, 3)
)


cov_grid <- cov_raw %>%
  mutate(
    age_norm        = (age          - age_min)   / (age_max   - age_min),
    BX_Delay_norm   = (biopsy_delay - bx_min)    / (bx_max    - bx_min),
    HR_p_norm       = (hs_status    - hr_min)    / (hr_max    - hr_min),
    Tgrade_norm     = (tumor_grade  - grade_min) / (grade_max - grade_min)
  )

# ——— 3. Extract the normalized X–matrix for prediction ———
X_grid_norm <- as.matrix(cov_grid[, c("age_norm", "BX_Delay_norm",
                                      "HR_p_norm", "Tgrade_norm")])


Mprime <- length(seq(burn_in + 1, n_mcmc, by = thin))
n_grid <- nrow(X_grid_norm)

# Pre‐allocate
CATE_big_draws    <- matrix(NA, Mprime, n_grid)
CRATE5_big_draws  <- matrix(NA, Mprime, n_grid)
CATE_small_draws  <- matrix(NA, Mprime, n_grid)
CRATE5_small_draws<- matrix(NA, Mprime, n_grid)

for(i in seq_len(n_grid)) {
  x_i <- X_grid_norm[i, ]

  # Big county
  tmp1 <- estimate_CATE_corrected(
    x            = x_i,
    v_i          = v_big,
    county_id    = big_county,
    forest_ps    = forest_ps,
    forest_out   = forest_out,
    w_samples    = w_samples,
    sigma2_samps = sigma2_samps,
    burn_in      = burn_in,
    thin         = thin,
    n_mcmc       = n_mcmc,
    n_iter_ps    = n_iter_ps
  )
  CATE_big_draws[, i]   <- tmp1$CATE_draws

  tmp2 <- estimate_CRATE(
    t_star       = 5*365,
    x            = x_i,
    v_i          = v_big,
    county_id    = big_county,
    forest_ps    = forest_ps,
    forest_out   = forest_out,
    w_samples    = w_samples,
    sigma2_samps = sigma2_samps,
    burn_in      = burn_in,
    thin         = thin,
    n_mcmc       = n_mcmc,
    n_iter_ps    = n_iter_ps
  )
  CRATE5_big_draws[, i] <- tmp2$CRATE_draws

  # Small county
  tmp3 <- estimate_CATE_corrected(
    x            = x_i,
    v_i          = v_small,
    county_id    = small_county,
    forest_ps    = forest_ps,
    forest_out   = forest_out,
    w_samples    = w_samples,
    sigma2_samps = sigma2_samps,
    burn_in      = burn_in,
    thin         = thin,
    n_mcmc       = n_mcmc,
    n_iter_ps    = n_iter_ps
  )
  CATE_small_draws[, i] <- tmp3$CATE_draws

  tmp4 <- estimate_CRATE(
    t_star       = 5*365,
    x            = x_i,
    v_i          = v_small,
    county_id    = small_county,
    forest_ps    = forest_ps,
    forest_out   = forest_out,
    w_samples    = w_samples,
    sigma2_samps = sigma2_samps,
    burn_in      = burn_in,
    thin         = thin,
    n_mcmc       = n_mcmc,
    n_iter_ps    = n_iter_ps
  )
  CRATE5_small_draws[, i] <- tmp4$CRATE_draws
}

summ <- function(mat) {
  mu  <- colMeans(mat)
  ci  <- apply(mat, 2, quantile, c(0.025, 0.975))
  data.frame(mean = mu, lower = ci[1,], upper = ci[2,])
}

res_big_cate   <- summ(CATE_big_draws)
res_big_crate5 <- summ(CRATE5_big_draws)
res_sml_cate   <- summ(CATE_small_draws)
res_sml_crate5 <- summ(CRATE5_small_draws)

results_grid <- cbind(
  cov_raw,
  res_big_cate,     res_big_crate5,
  res_sml_cate,     res_sml_crate5
)
head(results_grid)



# ——— Save full posterior‐draw matrices ———
write.csv(CATE_big_draws,    "CATE_big_draws_gpt_2.csv",    row.names = FALSE)
write.csv(CRATE5_big_draws,  "CRATE5_big_draws_gpt_2.csv",  row.names = FALSE)
write.csv(CATE_small_draws,  "CATE_small_draws_gpt_2.csv",  row.names = FALSE)
write.csv(CRATE5_small_draws,"CRATE5_small_draws_gpt_2.csv",row.names = FALSE)

# ——— Save the summarized grid ———
write.csv(results_grid,      "results_grid_2.csv",      row.names = FALSE)

message("All posterior draws and summaries saved to CSV files.")



# 1. Create draw‐column names (length = Mprime)
draw_names <- paste0("draw", seq_len(Mprime))

CATEBigSamples <- cbind(
  cov_raw,
  setNames(
    as.data.frame(t(CATE_big_draws), check.names = FALSE),
    draw_names
  )
)

CRATEFiveBigSamples <- cbind(
  cov_raw,
  setNames(
    as.data.frame(t(CRATE5_big_draws), check.names = FALSE),
    draw_names
  )
)

CATESmallSamples <- cbind(
  cov_raw,
  setNames(
    as.data.frame(t(CATE_small_draws), check.names = FALSE),
    draw_names
  )
)

CRATEFiveSmallSamples <- cbind(
  cov_raw,
  setNames(
    as.data.frame(t(CRATE5_small_draws), check.names = FALSE),
    draw_names
  )
)

# 3. Save to CSV
write.csv(CATEBigSamples,        "CATEBigSamples.csv",        row.names = FALSE)
write.csv(CRATEFiveBigSamples,   "CRATEFiveBigSamples.csv",   row.names = FALSE)
write.csv(CATESmallSamples,      "CATESmallSamples.csv",      row.names = FALSE)
write.csv(CRATEFiveSmallSamples, "CRATEFiveSmallSamples.csv", row.names = FALSE)

message("Saved four sample matrices with covariate combos (transposed) to CSV.")






summarize_draws <- function(draws_mat) {
  mu  <- colMeans(draws_mat)
  ci  <- apply(draws_mat, 2, quantile, probs = c(0.025, 0.975))
  data.frame(mean = mu, lower = ci[1,], upper = ci[2,])
}

# ——— 2. Compute summaries for each matrix ———
res_big_cate    <- summarize_draws(CATE_big_draws)
res_big_crate5  <- summarize_draws(CRATE5_big_draws)
res_small_cate  <- summarize_draws(CATE_small_draws)
res_small_crate5<- summarize_draws(CRATE5_small_draws)

# ——— 3. Assemble final summary table ———
summary_grid <- data.frame(
  cov_raw,

  # Big county summaries
  CATE_big_mean    = res_big_cate$mean,
  CATE_big_lower   = res_big_cate$lower,
  CATE_big_upper   = res_big_cate$upper,

  CRATE5_big_mean  = res_big_crate5$mean,
  CRATE5_big_lower = res_big_crate5$lower,
  CRATE5_big_upper = res_big_crate5$upper,

  # Small county summaries
  CATE_small_mean    = res_small_cate$mean,
  CATE_small_lower   = res_small_cate$lower,
  CATE_small_upper   = res_small_cate$upper,

  CRATE5_small_mean  = res_small_crate5$mean,
  CRATE5_small_lower = res_small_crate5$lower,
  CRATE5_small_upper = res_small_crate5$upper
)

# ——— 4. View or save ———
print(summary_grid)
# write.csv(summary_grid, "summary_grid.csv", row.names = FALSE)





# 0. Prereqs --------------------------------------------------------------
# CATE_draws_matrix: M' × n
# CRATE_draws_matrix: M' × n
# group: length-n vector of county IDs (1…K)
# Z: length-n treatment assignment (0/1)
Mprime <- nrow(CATE_draws_matrix)
counties <- sort(unique(group))
K       <- length(counties)

# 1. Prepare storage for posterior draws -------------------------------
ATE_county_draws    <- matrix(NA, Mprime, K)
ATT_county_draws    <- matrix(NA, Mprime, K)
ATU_county_draws    <- matrix(NA, Mprime, K)
CRATE_county_draws  <- matrix(NA, Mprime, K)
CRATT_county_draws  <- matrix(NA, Mprime, K)
CRATU_county_draws  <- matrix(NA, Mprime, K)

# 2. Sample sizes per county ---------------------------------------------
N_all        <- integer(K)
N_treated    <- integer(K)
N_untreated  <- integer(K)

# 3. Loop over counties & iterations ------------------------------------
for (k in seq_along(counties)) {
  cid <- counties[k]
  idx <- which(group == cid)

  N_all[k]       <- length(idx)
  idx_treat      <- idx[Z[idx] == 1]
  idx_untreat    <- idx[Z[idx] == 0]
  N_treated[k]   <- length(idx_treat)
  N_untreated[k] <- length(idx_untreat)

  for (m in 1:Mprime) {
    # individualized draws for this iteration
    cate_m  <- CATE_draws_matrix[m, idx]
    crate_m <- CRATE_draws_matrix[m, idx]

    # ATE, ATT, ATU draws
    ATE_county_draws[m, k] <- mean(cate_m)
    ATT_county_draws[m, k] <- if(length(idx_treat)>0) mean(CATE_draws_matrix[m, idx_treat]) else NA
    ATU_county_draws[m, k] <- if(length(idx_untreat)>0) mean(CATE_draws_matrix[m, idx_untreat]) else NA

    # CRATE@5, CRATT@5, CRATU@5 draws
    CRATE_county_draws[m, k] <- mean(crate_m)
    CRATT_county_draws[m, k] <- if(length(idx_treat)>0) mean(CRATE_draws_matrix[m, idx_treat]) else NA
    CRATU_county_draws[m, k] <- if(length(idx_untreat)>0) mean(CRATE_draws_matrix[m, idx_untreat]) else NA
  }
}

# 4. Summarize draws ------------------------------------------------------
summarize_draws <- function(mat) {
  mu  <- colMeans(mat, na.rm = TRUE)
  ci  <- apply(mat, 2, quantile, probs = c(0.025, 0.975), na.rm = TRUE)
  data.frame(mean = mu, lower = ci[1, ], upper = ci[2, ])
}

res_ATE   <- summarize_draws(ATE_county_draws)
res_ATT   <- summarize_draws(ATT_county_draws)
res_ATU   <- summarize_draws(ATU_county_draws)
res_CRATE <- summarize_draws(CRATE_county_draws)
res_CRATT <- summarize_draws(CRATT_county_draws)
res_CRATU <- summarize_draws(CRATU_county_draws)

# 5. Assemble county‐level summary table -------------------------------
county_summary <- data.frame(
  county       = counties,
  N_all        = N_all,
  N_treated    = N_treated,
  N_untreated  = N_untreated,

  ATE_mean     = res_ATE$mean,
  ATE_lower    = res_ATE$lower,
  ATE_upper    = res_ATE$upper,

  ATT_mean     = res_ATT$mean,
  ATT_lower    = res_ATT$lower,
  ATT_upper    = res_ATT$upper,

  ATU_mean     = res_ATU$mean,
  ATU_lower    = res_ATU$lower,
  ATU_upper    = res_ATU$upper,

  CRATE_mean   = res_CRATE$mean,
  CRATE_lower  = res_CRATE$lower,
  CRATE_upper  = res_CRATE$upper,

  CRATT_mean   = res_CRATT$mean,
  CRATT_lower  = res_CRATT$lower,
  CRATT_upper  = res_CRATT$upper,

  CRATU_mean   = res_CRATU$mean,
  CRATU_lower  = res_CRATU$lower,
  CRATU_upper  = res_CRATU$upper
)

# 6. Save to CSV ---------------------------------------------------------
write.csv(
  county_summary,
  file      = "county_level_treatment_effects_2.csv",
  row.names = FALSE
)

message("Saved county‐specific ATE/ATT/ATU and CRATE/CRATT/CRATU with CIs and sample sizes to 'county_level_treatment_effects.csv'")




library(dplyr)

county_summary_fmt <- county_summary %>%
  transmute(
    county,
    N_all,
    N_treated,
    N_untreated,

    ATE = sprintf("%.2f (%.2f, %.2f)",
                  ATE_mean, ATE_lower, ATE_upper),
    ATT = sprintf("%.2f (%.2f, %.2f)",
                  ATT_mean, ATT_lower, ATT_upper),
    ATU = sprintf("%.2f (%.2f, %.2f)",
                  ATU_mean, ATU_lower, ATU_upper),

    CRATE = sprintf("%.2f (%.2f, %.2f)",
                    CRATE_mean, CRATE_lower, CRATE_upper),
    CRATT = sprintf("%.2f (%.2f, %.2f)",
                    CRATT_mean, CRATT_lower, CRATT_upper),
    CRATU = sprintf("%.2f (%.2f, %.2f)",
                    CRATU_mean, CRATU_lower, CRATU_upper)
  )

# inspect
print(county_summary_fmt)

# save
write.csv(
  county_summary_fmt,
  "county_summary_formatted2.csv",
  row.names = FALSE
)





# ─── Libraries ───
library(tigris)    # county shapefiles
library(sf)        # spatial data
library(dplyr)     # data wrangling
library(ggplot2)   # plotting
library(patchwork) # combine plots
options(tigris_use_cache = TRUE)

# ─── 1. Load Florida counties and join your summaries ───
fl_counties <- counties(state = "FL", cb = TRUE, class = "sf")

county_summary <- county_summary %>%
  mutate(GEOID = as.character(county))

map_data <- fl_counties %>%
  left_join(county_summary, by = "GEOID")

# ─── 2. Compute common color‐scale ranges ───
r_ate   <- range(c(map_data$ATE_mean,   map_data$CRATE_mean),   na.rm = TRUE)
r_att   <- range(c(map_data$ATT_mean,   map_data$CRATT_mean),   na.rm = TRUE)
r_atu   <- range(c(map_data$ATU_mean,   map_data$CRATU_mean),   na.rm = TRUE)

# ─── 3. Build individual maps ───
theme_map <- theme_void() +
  theme(legend.position = "right")

# Figure 1: ATE & CRATE
p_ate   <- ggplot(map_data) +
  geom_sf(aes(fill = ATE_mean), color = NA) +
  scale_fill_viridis_c(limits = r_ate, name = "Post. Mean") +
  theme_map

p_crate <- ggplot(map_data) +
  geom_sf(aes(fill = CRATE_mean), color = NA) +
  scale_fill_viridis_c(limits = r_ate, name = "Post. Mean") +
  theme_map

fig1 <- p_ate + p_crate +
  plot_layout(ncol = 2, guides = "collect") &
  theme(legend.position = "right")

# Figure 2: ATT & CRATT
p_att   <- ggplot(map_data) +
  geom_sf(aes(fill = ATT_mean), color = NA) +
  scale_fill_viridis_c(limits = r_att, name = "Post. Mean") +
  theme_map

p_cratt <- ggplot(map_data) +
  geom_sf(aes(fill = CRATT_mean), color = NA) +
  scale_fill_viridis_c(limits = r_att, name = "Post. Mean") +
  theme_map

fig2 <- p_att + p_cratt +
  plot_layout(ncol = 2, guides = "collect") &
  theme(legend.position = "right")

# Figure 3: ATU & CRATU
p_atu   <- ggplot(map_data) +
  geom_sf(aes(fill = ATU_mean), color = NA) +
  scale_fill_viridis_c(limits = r_atu, name = "Post. Mean") +
  theme_map

p_cratu <- ggplot(map_data) +
  geom_sf(aes(fill = CRATU_mean), color = NA) +
  scale_fill_viridis_c(limits = r_atu, name = "Post. Mean") +
  theme_map

fig3 <- p_atu + p_cratu +
  plot_layout(ncol = 2, guides = "collect") &
  theme(legend.position = "right")

# ─── 4. Display the figures ───
fig1
fig2
fig3




ATE_pop_draws  <- rowMeans(CATE_draws_matrix)
ATT_pop_draws  <- apply(CATE_draws_matrix[, Z == 1, drop = FALSE], 1, mean)
ATU_pop_draws  <- apply(CATE_draws_matrix[, Z == 0, drop = FALSE], 1, mean)

RATE_pop_draws <- rowMeans(CRATE_draws_matrix)
RATT_pop_draws <- apply(CRATE_draws_matrix[, Z == 1, drop = FALSE], 1, mean)
RATU_pop_draws <- apply(CRATE_draws_matrix[, Z == 0, drop = FALSE], 1, mean)


summarize_draws <- function(draws) {
  m   <- mean(draws)
  ci  <- quantile(draws, probs = c(0.025, 0.975))
  c(mean = m, lower = ci[1], upper = ci[2])
}


pop_df <- data.frame(
  Estimand = rownames(res_pop),
  mean      = res_pop[, 1],
  lower     = res_pop[, 2],
  upper     = res_pop[, 3],
  stringsAsFactors = FALSE
)

library(dplyr)

pop_summary <- pop_df %>%
  mutate(
    `Estimate (95% CI)` = sprintf(
      "%.2f (%.2f, %.2f)",
      mean, lower, upper
    )
  ) %>%
  select(Estimand, `Estimate (95% CI)`)

print(pop_summary)

# 3.LaTeX
library(knitr)

cat(
  kable(
    pop_summary,
    format   = "latex",
    booktabs = TRUE,
    caption  = "Population‐level treatment effects and 95\\% credible intervals",
    label    = "tab:pop_effects",
    align    = c("l","r")
  )
)
pop_summary <- pop_df %>%
  mutate(
    `Estimate (95% CI)` = sprintf(
      "%.2f (%.2f, %.2f)",
      mean, lower, upper
    )
  ) %>%
  select(Estimand, `Estimate (95% CI)`)

# --- 4. Inspect
print(pop_summary)

# turn into a data.frame and format
library(tibble)
library(dplyr)

pop_summary <- as.data.frame(res_pop) %>%
  rownames_to_column("Estimand") %>%
  mutate(
    `Estimate (95% CI)` = sprintf(
      "%.2f (%.2f, %.2f)",
      mean, lower, upper
    )
  ) %>%
  select(Estimand, `Estimate (95% CI)`)

# view
print(pop_summary)


library(knitr)

latex_table <- kable(
  pop_summary,
  format    = "latex",
  booktabs  = TRUE,
  caption   = "Population‐level treatment effects and 95\\% credible intervals",
  label     = "tab:pop_effects",
  align     = c("l", "r")
)

cat(latex_table)



pop_df_months <- data.frame(
  Estimand = rownames(res_pop),
  mean      = res_pop[, 1] / 30,
  lower     = res_pop[, 2] / 30,
  upper     = res_pop[, 3] / 30,
  stringsAsFactors = FALSE
)
# Round to 2 decimals
library(dplyr)
pop_df_months <- pop_df_months %>%
  mutate(across(c(mean, lower, upper), ~ round(.x, 2)))

pop_table <- pop_df_months %>%
  select(Estimand, mean, lower, upper)


print(pop_table)



library(knitr)


latex_code <- kable(
  pop_table,
  format    = "latex",
  booktabs  = TRUE,
  caption   = "Population‐level treatment effects (in months) and 95\\% credible intervals",
  label     = "tab:pop_effects",
  align     = c("l", "r", "r", "r"),
  digits    = 2
)

# 2. Print it
cat(latex_code)








------------------------------------------------------------------

t_star <- 5 * 365.25
# Get dimensions
n_patients <- nrow(data_stratum)
n_counties <- length(stratum_counties)
M_samples <- length(sigma2_samps)  # Number of posterior samples

cat("Computing CATE and CRATE for", n_patients, "patients across", n_counties, "counties...\n")
cat("Using", M_samples, "posterior samples\n")
cat("Time horizon:", t_star, "days (5 years)\n\n")

# ------------------------------------------------------------------------------
# 1. CATE (Conditional Average Treatment Effect) for each patient
# ------------------------------------------------------------------------------

# Initialize  for CATE results
CATE_results <- list(
  samples = matrix(0, nrow = M_samples, ncol = n_patients),  # M x n matrix
  means = numeric(n_patients),                               # length n vector
  lower_ci = numeric(n_patients),                           # 2.5% quantiles
  upper_ci = numeric(n_patients),                           # 97.5% quantiles
  patient_info = data.frame(
    patient_id = 1:n_patients,
    county_id = data_stratum$county,
    county_index = group,
    Age = data_stratum$Age,
    BX_Delay = data_stratum$BX_Delay,
    HR_p = data_stratum$HR_p,
    Tgrade = data_stratum$Tgrade,
    observed_Z = Z,
    observed_time = time,
    observed_status = status
  )
)

# Compute CATE for each patient
cat("Computing CATE for each patient...\n")
pb <- txtProgressBar(min = 0, max = n_patients, style = 3)

for (j in 1:n_patients) {
  # Extract patient j's covariates
  x_j <- X[j, ]                    # Patient-level covariates (scaled)
  v_j <- V[j, ]                    # County-level covariates
  county_j <- group[j]             # County index for patient j

  # Compute CATE for patient j
  cate_j <- estimate_CATE_corrected(
    x = x_j,
    v_i = v_j,
    county_id = county_j,
    forest_ps = forest_ps,
    forest_out = forest_out,
    w_samples = w_samples,
    sigma2_samps = sigma2_samps,
    burn_in = burn_in,
    thin = thin,
    n_mcmc = n_mcmc,
    n_iter_ps = n_iter_ps
  )

  # Store results
  CATE_results$samples[, j] <- cate_j$CATE_draws
  CATE_results$means[j] <- cate_j$CATE_mean
  CATE_results$lower_ci[j] <- cate_j$CATE_lower
  CATE_results$upper_ci[j] <- cate_j$CATE_upper

  setTxtProgressBar(pb, j)
}
close(pb)

# ------------------------------------------------------------------------------
# 2. CRATE (Conditional Restricted Average Treatment Effect) at 5 years
# ------------------------------------------------------------------------------

# Initialize storage for CRATE results
CRATE_results <- list(
  samples = matrix(0, nrow = M_samples, ncol = n_patients),  # M x n matrix
  means = numeric(n_patients),                               # length n vector
  lower_ci = numeric(n_patients),                           # 2.5% quantiles
  upper_ci = numeric(n_patients),                           # 97.5% quantiles
  t_star = t_star,                                          # Time horizon used
  patient_info = CATE_results$patient_info                  # Same patient info
)

# Compute CRATE at 5 years for each patient
cat("\nComputing CRATE at 5 years for each patient...\n")
pb <- txtProgressBar(min = 0, max = n_patients, style = 3)

for (j in 1:n_patients) {
  # Extract patient j's covariates
  x_j <- X[j, ]                    # Patient-level covariates (scaled)
  v_j <- V[j, ]                    # County-level covariates
  county_j <- group[j]             # County index for patient j

  # Compute CRATE for patient j at time t_star
  crate_j <- estimate_CRATE(
    t_star = t_star,
    x = x_j,
    v_i = v_j,
    county_id = county_j,
    forest_ps = forest_ps,
    forest_out = forest_out,
    w_samples = w_samples,
    sigma2_samps = sigma2_samps,
    burn_in = burn_in,
    thin = thin,
    n_mcmc = n_mcmc,
    n_iter_ps = n_iter_ps
  )

  # Store results
  CRATE_results$samples[, j] <- crate_j$CRATE_draws
  CRATE_results$means[j] <- crate_j$CRATE_mean
  CRATE_results$lower_ci[j] <- crate_j$CRATE_lower
  CRATE_results$upper_ci[j] <- crate_j$CRATE_upper

  setTxtProgressBar(pb, j)
}
close(pb)

# ------------------------------------------------------------------------------
# 3. Summary and Visualization
# ------------------------------------------------------------------------------

summary_df <- data.frame(
  patient_id = 1:n_patients,
  county_id = data_stratum$county,
  county_index = group,
  observed_Z = Z,
  observed_time = time,
  observed_status = status,
  Age = data_stratum$Age,
  BX_Delay = data_stratum$BX_Delay,
  HR_p = data_stratum$HR_p,
  Tgrade = data_stratum$Tgrade,
  CATE_mean = CATE_results$means,
  CATE_lower = CATE_results$lower_ci,
  CATE_upper = CATE_results$upper_ci,
  CRATE_mean = CRATE_results$means,
  CRATE_lower = CRATE_results$lower_ci,
  CRATE_upper = CRATE_results$upper_ci
)

# Print summary statistics
cat("\n" + "="*80 + "\n")
cat("SUMMARY RESULTS\n")
cat("="*80 + "\n")

cat("CATE (Conditional Average Treatment Effect):\n")
cat("  Mean across all patients:", round(mean(CATE_results$means), 3), "\n")
cat("  Range: [", round(min(CATE_results$means), 3), ", ", round(max(CATE_results$means), 3), "]\n")
cat("  SD:", round(sd(CATE_results$means), 3), "\n\n")

cat("CRATE at 5 years (Conditional Restricted Average Treatment Effect):\n")
cat("  Mean across all patients:", round(mean(CRATE_results$means), 3), "\n")
cat("  Range: [", round(min(CRATE_results$means), 3), ", ", round(max(CRATE_results$means), 3), "]\n")
cat("  SD:", round(sd(CRATE_results$means), 3), "\n\n")

# Show first few patients
cat("First 10 patients:\n")
print(head(summary_df, 10))


# Save posterior samples
cat("\nSaving results...\n")

# Save CATE samples (M x n matrix)
write.csv(CATE_results$samples, "CATE_posterior_samples_claudio.csv", row.names = FALSE)

# Save CRATE samples (M x n matrix)
write.csv(CRATE_results$samples, "CRATE_posterior_samples_claudio.csv", row.names = FALSE)

# Save summary table
write.csv(summary_df, "CATE_CRATE_summary_claudio.csv", row.names = FALSE)

# Save complete results as RData
save(CATE_results, CRATE_results, summary_df,
     file = "CATE_CRATE_complete_results_claudio.RData")

cat("Results saved to:\n")
cat("  - CATE_posterior_samples_claudio.csv (", M_samples, " x ", n_patients, " matrix)\n")
cat("  - CRATE_posterior_samples_claudio.csv (", M_samples, " x ", n_patients, " matrix)\n")
cat("  - CATE_CRATE_summary_claudio.csv (summary with means and CIs)\n")
cat("  - CATE_CRATE_complete_results_claudio.RData (complete R objects)\n")


# Example: Get results for patient 1
patient_idx <- 1
cat("\n" + "="*50 + "\n")
cat("EXAMPLE: Results for Patient", patient_idx, "\n")
cat("="*50 + "\n")

cat("Patient characteristics:\n")
cat("  County ID:", summary_df$county_id[patient_idx], "\n")
cat("  Age:", summary_df$Age[patient_idx], "\n")
cat("  BX_Delay:", summary_df$BX_Delay[patient_idx], "\n")
cat("  HR_p:", summary_df$HR_p[patient_idx], "\n")
cat("  Tgrade:", summary_df$Tgrade[patient_idx], "\n")
cat("  Observed treatment (Z):", summary_df$observed_Z[patient_idx], "\n")
cat("  Observed time:", summary_df$observed_time[patient_idx], "\n")
cat("  Observed status:", summary_df$observed_status[patient_idx], "\n\n")

cat("CATE:\n")
cat("  Posterior mean:", round(CATE_results$means[patient_idx], 3), "\n")
cat("  95% CI: [", round(CATE_results$lower_ci[patient_idx], 3), ", ",
    round(CATE_results$upper_ci[patient_idx], 3), "]\n")
cat("  Posterior samples: ", length(CATE_results$samples[, patient_idx]), " draws\n\n")

cat("CRATE at 5 years:\n")
cat("  Posterior mean:", round(CRATE_results$means[patient_idx], 3), "\n")
cat("  95% CI: [", round(CRATE_results$lower_ci[patient_idx], 3), ", ",
    round(CRATE_results$upper_ci[patient_idx], 3), "]\n")
cat("  Posterior samples: ", length(CRATE_results$samples[, patient_idx]), " draws\n")


# Plot distributions of CATE and CRATE means across patients
library(ggplot2)

# CATE distribution
p1 <- ggplot(summary_df, aes(x = CATE_mean)) +
  geom_histogram(bins = 30, fill = "steelblue", alpha = 0.7) +
  labs(title = "Distribution of CATE Posterior Means Across Patients",
       x = "CATE Posterior Mean",
       y = "Count") +
  theme_minimal()

# CRATE distribution
p2 <- ggplot(summary_df, aes(x = CRATE_mean)) +
  geom_histogram(bins = 30, fill = "darkgreen", alpha = 0.7) +
  labs(title = "Distribution of CRATE (5-year) Posterior Means Across Patients",
       x = "CRATE Posterior Mean",
       y = "Count") +
  theme_minimal()

# Print plots
print(p1)
print(p2)

# Save plots
ggsave("CATE_distribution.png", p1, width = 8, height = 6)
ggsave("CRATE_distribution.png", p2, width = 8, height = 6)

cat("\nAnalysis complete!\n")
cat("You now have:\n")
cat("1. CATE_results$samples: ", M_samples, " x ", n_patients, " matrix of CATE posterior samples\n")
cat("2. CRATE_results$samples: ", M_samples, " x ", n_patients, " matrix of CRATE posterior samples\n")
cat("3. summary_df: Data frame with posterior means and 95% CIs for each patient\n")
cat("4. All results saved to CSV and RData files\n")







# ------------------------------------------------------------------------------
# CATE and CRATE Analysis for All Covariate Combinations
# Big County vs Small County Comparison
# ------------------------------------------------------------------------------

# Define time horizon (5 years)
t_star <- 5 * 365.25  # 5 years in days

# ------------------------------------------------------------------------------
# 1. Create all combinations of covariates
# ------------------------------------------------------------------------------

# Define covariate levels (you can adjust these based on your data range)
age_levels <- c(50, 60, 80)  # Young, middle-aged, older
biopsy_delay_levels <- c(0, 1)  # No delay, delay
hs_status_levels <- c(0, 1)     # HR negative, HR positive
tumor_grade_levels <- c(1, 2, 3)  # Low, intermediate, high grade

# Create all combinations
covariate_combinations <- expand.grid(
  Age = age_levels,
  BX_Delay = biopsy_delay_levels,
  HR_p = hs_status_levels,
  Tgrade = tumor_grade_levels
)

n_combinations <- nrow(covariate_combinations)
cat("Total covariate combinations:", n_combinations, "\n")
print(covariate_combinations)

# Scale the combinations to [0,1] using the same scaling as original data
# (Use the min/max from your original training data for consistency)
scale01_with_range <- function(x, min_val, max_val) {
  if(max_val == min_val) return(rep(0, length(x)))
  (x - min_val) / (max_val - min_val)
}

# Get original data ranges for consistent scaling
age_min <- min(data_stratum$Age, na.rm = TRUE)
age_max <- max(data_stratum$Age, na.rm = TRUE)
bx_delay_min <- min(data_stratum$BX_Delay, na.rm = TRUE)
bx_delay_max <- max(data_stratum$BX_Delay, na.rm = TRUE)
hr_min <- min(data_stratum$HR_p, na.rm = TRUE)
hr_max <- max(data_stratum$HR_p, na.rm = TRUE)
tgrade_min <- min(data_stratum$Tgrade, na.rm = TRUE)
tgrade_max <- max(data_stratum$Tgrade, na.rm = TRUE)

# Scale covariate combinations
X_combinations_scaled <- data.frame(
  Age = scale01_with_range(covariate_combinations$Age, age_min, age_max),
  BX_Delay = scale01_with_range(covariate_combinations$BX_Delay, bx_delay_min, bx_delay_max),
  HR_p = scale01_with_range(covariate_combinations$HR_p, hr_min, hr_max),
  Tgrade = scale01_with_range(covariate_combinations$Tgrade, tgrade_min, tgrade_max)
)
X_combinations_scaled <- as.matrix(X_combinations_scaled)

cat("\nScaled covariate combinations:\n")
print(round(X_combinations_scaled, 3))


# Count patients per county in the stratum
county_sizes <- table(data_stratum$county)
cat("\nCounty sizes in stratum:\n")
print(county_sizes)

# Define big and small counties (you can adjust these criteria)
median_size <- median(county_sizes)
big_counties <- names(county_sizes)[county_sizes >= median_size]
small_counties <- names(county_sizes)[county_sizes < median_size]

big_county_id <- as.numeric(big_counties[which.max(county_sizes[big_counties])])  # Largest county
small_county_id <- as.numeric(small_counties[1])  # First small county

# Get their indices in the stratum
big_county_index <- which(stratum_counties == big_county_id)
small_county_index <- which(stratum_counties == small_county_id)

cat("\nSelected counties:\n")
cat("Big county ID:", big_county_id, "(index:", big_county_index, ") with", county_sizes[as.character(big_county_id)], "patients\n")
cat("Small county ID:", small_county_id, "(index:", small_county_index, ") with", county_sizes[as.character(small_county_id)], "patients\n")

# County-level covariates (assuming intercept-only as in original analysis)
v_county <- matrix(1, nrow = 1, ncol = 1)  # Intercept


compute_effects_for_combinations <- function(county_index, county_name) {
  cat("\n", paste(rep("=", 60), collapse=""), "\n")
  cat("Computing effects for", county_name, "(County", stratum_counties[county_index], ")\n")
  cat(paste(rep("=", 60), collapse=""), "\n")


  n_covariate_cols <- 4  # Age, BX_Delay, HR_p, Tgrade

  CATE_results <- list(
    samples_with_covariates = matrix(0, nrow = n_combinations, ncol = n_covariate_cols + M_total),
    samples = matrix(0, nrow = M_total, ncol = n_combinations),
    means = numeric(n_combinations),
    lower_ci = numeric(n_combinations),
    upper_ci = numeric(n_combinations)
  )

  CRATE_results <- list(
    samples_with_covariates = matrix(0, nrow = n_combinations, ncol = n_covariate_cols + M_total),
    samples = matrix(0, nrow = M_total, ncol = n_combinations),
    means = numeric(n_combinations),
    lower_ci = numeric(n_combinations),
    upper_ci = numeric(n_combinations)
  )

  CATE_results$samples_with_covariates[, 1:4] <- as.matrix(covariate_combinations)
  CRATE_results$samples_with_covariates[, 1:4] <- as.matrix(covariate_combinations)

  # Add column names
  covariate_names <- c("Age", "BX_Delay", "HR_p", "Tgrade")
  sample_names <- paste0("Sample_", 1:M_total)
  colnames(CATE_results$samples_with_covariates) <- c(covariate_names, sample_names)
  colnames(CRATE_results$samples_with_covariates) <- c(covariate_names, sample_names)

  # Progress bar
  pb <- txtProgressBar(min = 0, max = n_combinations, style = 3)

  for (k in 1:n_combinations) {

    x_k <- X_combinations_scaled[k, ]

    # Compute CATE for combination k
    cate_k <- estimate_CATE_corrected(
      x = x_k,
      v_i = v_county,
      county_id = county_index,
      forest_ps = forest_ps,
      forest_out = forest_out,
      w_samples = w_samples,
      sigma2_samps = sigma2_samps,
      burn_in = burn_in,
      thin = thin,
      n_mcmc = n_mcmc,
      n_iter_ps = n_iter_ps
    )

    # Compute CRATE for combination k
    crate_k <- estimate_CRATE(
      t_star = t_star,
      x = x_k,
      v_i = v_county,
      county_id = county_index,
      forest_ps = forest_ps,
      forest_out = forest_out,
      w_samples = w_samples,
      sigma2_samps = sigma2_samps,
      burn_in = burn_in,
      thin = thin,
      n_mcmc = n_mcmc,
      n_iter_ps = n_iter_ps
    )

    # Store CATE results
    CATE_results$samples[, k] <- cate_k$CATE_draws
    CATE_results$samples_with_covariates[k, (n_covariate_cols + 1):(n_covariate_cols + M_total)] <- cate_k$CATE_draws
    CATE_results$means[k] <- cate_k$CATE_mean
    CATE_results$lower_ci[k] <- cate_k$CATE_lower
    CATE_results$upper_ci[k] <- cate_k$CATE_upper

    # Store CRATE results
    CRATE_results$samples[, k] <- crate_k$CRATE_draws
    CRATE_results$samples_with_covariates[k, (n_covariate_cols + 1):(n_covariate_cols + M_total)] <- crate_k$CRATE_draws
    CRATE_results$means[k] <- crate_k$CRATE_mean
    CRATE_results$lower_ci[k] <- crate_k$CRATE_lower
    CRATE_results$upper_ci[k] <- crate_k$CRATE_upper

    setTxtProgressBar(pb, k)
  }
  close(pb)

  return(list(CATE = CATE_results, CRATE = CRATE_results))
}


big_county_results <- compute_effects_for_combinations(big_county_index, "Big County")


small_county_results <- compute_effects_for_combinations(small_county_index, "Small County")


# Create summary data frame
create_summary_df <- function(results, county_name, county_id) {
  data.frame(
    County_Type = county_name,
    County_ID = county_id,
    Combination_ID = 1:n_combinations,
    Age = covariate_combinations$Age,
    BX_Delay = covariate_combinations$BX_Delay,
    HR_p = covariate_combinations$HR_p,
    Tgrade = covariate_combinations$Tgrade,
    CATE_mean = results$CATE$means,
    CATE_lower = results$CATE$lower_ci,
    CATE_upper = results$CATE$upper_ci,
    CRATE_mean = results$CRATE$means,
    CRATE_lower = results$CRATE$lower_ci,
    CRATE_upper = results$CRATE$upper_ci,
    CATE_width = results$CATE$upper_ci - results$CATE$lower_ci,
    CRATE_width = results$CRATE$upper_ci - results$CRATE$lower_ci
  )
}

big_county_summary <- create_summary_df(big_county_results, "Big County", big_county_id)
small_county_summary <- create_summary_df(small_county_results, "Small County", small_county_id)


combined_summary <- rbind(big_county_summary, small_county_summary)

# ------------------------------------------------------------------------------
# 7. Display results and comparisons
# ------------------------------------------------------------------------------

cat("\n", paste(rep("=", 80), collapse=""), "\n")
cat("SUMMARY RESULTS FOR ALL COVARIATE COMBINATIONS\n")
cat(paste(rep("=", 80), collapse=""), "\n")

# Big County Results
cat("\nBIG COUNTY (ID:", big_county_id, ") RESULTS:\n")
cat("CATE - Mean:", round(mean(big_county_results$CATE$means), 3),
    ", Range: [", round(min(big_county_results$CATE$means), 3),
    ", ", round(max(big_county_results$CATE$means), 3), "]\n")
cat("CRATE - Mean:", round(mean(big_county_results$CRATE$means), 3),
    ", Range: [", round(min(big_county_results$CRATE$means), 3),
    ", ", round(max(big_county_results$CRATE$means), 3), "]\n")

# Small County Results
cat("\nSMALL COUNTY (ID:", small_county_id, ") RESULTS:\n")
cat("CATE - Mean:", round(mean(small_county_results$CATE$means), 3),
    ", Range: [", round(min(small_county_results$CATE$means), 3),
    ", ", round(max(small_county_results$CATE$means), 3), "]\n")
cat("CRATE - Mean:", round(mean(small_county_results$CRATE$means), 3),
    ", Range: [", round(min(small_county_results$CRATE$means), 3),
    ", ", round(max(small_county_results$CRATE$means), 3), "]\n")


cat("\nFirst 10 combinations for both counties:\n")
print(head(combined_summary[, c("County_Type", "Age", "BX_Delay", "HR_p", "Tgrade",
                                "CATE_mean", "CATE_lower", "CATE_upper",
                                "CRATE_mean", "CRATE_lower", "CRATE_upper")], 20))


cat("\nSaving results...\n")

# Create the named matrices as requested
CATEBigSamples_claudio<- big_county_results$CATE$samples_with_covariates
CRATEFiveBigSamples_claudio <- big_county_results$CRATE$samples_with_covariates
CATESmallSamples_claudio <- small_county_results$CATE$samples_with_covariates
CRATEFiveSmallSamples_claudio <- small_county_results$CRATE$samples_with_covariates

# Save the matrices with covariate combinations included
write.csv(CATEBigSamples_claudio, "CATEBigSamples_claudio.csv", row.names = FALSE)
write.csv(CRATEFiveBigSamples_claudio, "CRATEFiveBigSamples_claudio.csv", row.names = FALSE)
write.csv(CATESmallSamples_claudio, "CATESmallSamples_claudio.csv", row.names = FALSE)
write.csv(CRATEFiveSmallSamples_claudio, "CRATEFiveSmallSamples_claudio.csv", row.names = FALSE)

# Also save the original format (just posterior samples)
write.csv(big_county_results$CATE$samples, "big_county_CATE_samples_only.csv", row.names = FALSE)
write.csv(big_county_results$CRATE$samples, "big_county_CRATE_samples_only.csv", row.names = FALSE)
write.csv(small_county_results$CATE$samples, "small_county_CATE_samples_only.csv", row.names = FALSE)
write.csv(small_county_results$CRATE$samples, "small_county_CRATE_samples_only.csv", row.names = FALSE)

# Save summary tables
write.csv(big_county_summary, "big_county_summary.csv", row.names = FALSE)
write.csv(small_county_summary, "small_county_summary.csv", row.names = FALSE)
write.csv(combined_summary, "combined_county_summary.csv", row.names = FALSE)

# Save covariate combinations
write.csv(covariate_combinations, "covariate_combinations.csv", row.names = FALSE)

# Save complete results including the new matrices
save(big_county_results, small_county_results,
     big_county_summary, small_county_summary, combined_summary,
     covariate_combinations, X_combinations_scaled,
     big_county_id, small_county_id, big_county_index, small_county_index,
     CATEBigSamples, CRATEFiveBigSamples, CATESmallSamples, CRATEFiveSmallSamples,
     file = "complete_combinations_analysis.RData")

cat("Results saved to:\n")
cat("  - CATEBigSamples.csv (", n_combinations, " x ", ncol(CATEBigSamples), " with covariates + samples)\n")
cat("  - CRATEFiveBigSamples.csv (", n_combinations, " x ", ncol(CRATEFiveBigSamples), " with covariates + samples)\n")
cat("  - CATESmallSamples.csv (", n_combinations, " x ", ncol(CATESmallSamples), " with covariates + samples)\n")
cat("  - CRATEFiveSmallSamples.csv (", n_combinations, " x ", ncol(CRATEFiveSmallSamples), " with covariates + samples)\n")
cat("  - combined_county_summary.csv (summary with means and CIs)\n")
cat("  - complete_combinations_analysis.RData (all R objects)\n")

# Display structure of the new matrices
cat("\nMatrix structures:\n")
cat("CATEBigSamples: ", nrow(CATEBigSamples), " combinations x ", ncol(CATEBigSamples), " columns\n")
cat("  - Columns 1-4: Age, BX_Delay, HR_p, Tgrade\n")
cat("  - Columns 5-", ncol(CATEBigSamples), ": Posterior samples (Sample_1 to Sample_", M_total, ")\n")
cat("CRATEFiveBigSamples: Same structure as CATEBigSamples\n")
cat("CATESmallSamples: Same structure as CATEBigSamples\n")
cat("CRATEFiveSmallSamples: Same structure as CATEBigSamples\n")


library(ggplot2)
library(reshape2)

# Prepare data for plotting
plot_data <- combined_summary %>%
  select(County_Type, Combination_ID, Age, BX_Delay, HR_p, Tgrade,
         CATE_mean, CRATE_mean) %>%
  melt(id.vars = c("County_Type", "Combination_ID", "Age", "BX_Delay", "HR_p", "Tgrade"),
       variable.name = "Effect_Type", value.name = "Effect_Value")

# Create combination labels
combined_summary$Combo_Label <- paste0("Age:", combined_summary$Age,
                                       ", BX:", combined_summary$BX_Delay,
                                       ", HR:", combined_summary$HR_p,
                                       ", Grade:", combined_summary$Tgrade)

# CATE comparison plot
p1 <- ggplot(combined_summary, aes(x = factor(Combination_ID), y = CATE_mean,
                                   color = County_Type, group = County_Type)) +
  geom_point(size = 2) +
  geom_line() +
  geom_errorbar(aes(ymin = CATE_lower, ymax = CATE_upper), width = 0.2, alpha = 0.6) +
  labs(title = "CATE Comparison: Big County vs Small County",
       subtitle = "Across All Covariate Combinations",
       x = "Covariate Combination",
       y = "CATE (Posterior Mean ± 95% CI)",
       color = "County Type") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_color_manual(values = c("Big County" = "blue", "Small County" = "red"))

# CRATE comparison plot
p2 <- ggplot(combined_summary, aes(x = factor(Combination_ID), y = CRATE_mean,
                                   color = County_Type, group = County_Type)) +
  geom_point(size = 2) +
  geom_line() +
  geom_errorbar(aes(ymin = CRATE_lower, ymax = CRATE_upper), width = 0.2, alpha = 0.6) +
  labs(title = "CRATE (5-year) Comparison: Big County vs Small County",
       subtitle = "Across All Covariate Combinations",
       x = "Covariate Combination",
       y = "CRATE (Posterior Mean ± 95% CI)",
       color = "County Type") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_color_manual(values = c("Big County" = "blue", "Small County" = "red"))

# Print plots
print(p1)
print(p2)

# Save plots
ggsave("CATE_big_vs_small_county.png", p1, width = 12, height = 8)
ggsave("CRATE_big_vs_small_county.png", p2, width = 12, height = 8)








# Load required library for county data
library(dplyr)

--------------------------------------


florida_counties <- data.frame(
  county_id = 1:67,
  county_name = c(
    "Alachua", "Baker", "Bay", "Bradford", "Brevard", "Broward", "Calhoun",
    "Charlotte", "Citrus", "Clay", "Collier", "Columbia", "DeSoto", "Dixie",
    "Duval", "Escambia", "Flagler", "Franklin", "Gadsden", "Gilchrist",
    "Glades", "Gulf", "Hamilton", "Hardee", "Hendry", "Hernando", "Highlands",
    "Hillsborough", "Holmes", "Indian River", "Jackson", "Jefferson", "Lafayette",
    "Lake", "Lee", "Leon", "Levy", "Liberty", "Madison", "Manatee", "Marion",
    "Martin", "Miami-Dade", "Monroe", "Nassau", "Okaloosa", "Okeechobee",
    "Orange", "Osceola", "Palm Beach", "Pasco", "Pinellas", "Polk", "Putnam",
    "Santa Rosa", "Sarasota", "Seminole", "St. Johns", "St. Lucie", "Sumter",
    "Suwannee", "Taylor", "Union", "Volusia", "Wakulla", "Walton", "Washington"
  )
)

# Display the mapping
cat("Florida County ID to Name Mapping:\n")
print(florida_counties)

# ------------------------------------------------------------------------------
# 2. Add County Names to Summary Data
# ------------------------------------------------------------------------------

# Add county names to the summary data
combined_summary_months_with_names <- combined_summary_months %>%
  left_join(florida_counties, by = c("County_ID" = "county_id"))

big_county_summary_months_with_names <- big_county_summary_months %>%
  left_join(florida_counties, by = c("County_ID" = "county_id"))

small_county_summary_months_with_names <- small_county_summary_months %>%
  left_join(florida_counties, by = c("County_ID" = "county_id"))

# Display which counties we're working with
cat("\nCounties in analysis:\n")
cat("Big County:", big_county_id, "=", florida_counties$county_name[florida_counties$county_id == big_county_id], "\n")
cat("Small County:", small_county_id, "=", florida_counties$county_name[florida_counties$county_id == small_county_id], "\n")

# Show all counties in the stratum
counties_in_stratum <- florida_counties[florida_counties$county_id %in% stratum_counties, ]
cat("\nAll counties in stratum (Race =", race_val, ", Stage =", stage_val, "):\n")
print(counties_in_stratum)

# ------------------------------------------------------------------------------
# 3. Updated LaTeX Table Generation Function
# ------------------------------------------------------------------------------

generate_latex_table_with_names <- function(data, caption, label, table_type = "combined") {

  # Round numerical values to 3 decimal places
  data_rounded <- data
  numeric_cols <- c("CATE_mean", "CATE_lower", "CATE_upper",
                    "CRATE_mean", "CRATE_lower", "CRATE_upper")
  for(col in numeric_cols) {
    if(col %in% colnames(data_rounded)) {
      data_rounded[[col]] <- round(data_rounded[[col]], 3)
    }
  }


  data_rounded$BX_Delay_label <- ifelse(data_rounded$BX_Delay == 1, "Yes", "No")
  data_rounded$HR_p_label <- ifelse(data_rounded$HR_p == 1, "Negative", "Positive")

  # Create formatted CATE and CRATE columns: Mean (Lower, Upper)
  data_rounded$CATE_formatted <- paste0(data_rounded$CATE_mean, " (",
                                        data_rounded$CATE_lower, ", ",
                                        data_rounded$CATE_upper, ")")
  data_rounded$CRATE_formatted <- paste0(data_rounded$CRATE_mean, " (",
                                         data_rounded$CRATE_lower, ", ",
                                         data_rounded$CRATE_upper, ")")

  # Start LaTeX table
  latex_code <- "\\begin{table}[htbp]\n"
  latex_code <- paste0(latex_code, "\\centering\n")
  latex_code <- paste0(latex_code, "\\caption{", caption, "}\n")
  latex_code <- paste0(latex_code, "\\label{", label, "}\n")

  if(table_type == "combined") {
    # Combined table with county names
    latex_code <- paste0(latex_code, "\\begin{tabular}{lcccccc}\n")
    latex_code <- paste0(latex_code, "\\hline\n")
    latex_code <- paste0(latex_code, "County & Age & BX & HR & Grade & CATE & CRATE \\\\\n")
    latex_code <- paste0(latex_code, "Name &  & Delay & Status &  & (Months) & (Months) \\\\\n")
    latex_code <- paste0(latex_code, "\\hline\n")

    for(i in 1:nrow(data_rounded)) {
      # Use county name, handling any potential missing values
      county_display <- ifelse(is.na(data_rounded$county_name[i]),
                               paste0("County ", data_rounded$County_ID[i]),
                               data_rounded$county_name[i])

      latex_code <- paste0(latex_code,
                           county_display, " & ",
                           data_rounded$Age[i], " & ",
                           data_rounded$BX_Delay_label[i], " & ",
                           data_rounded$HR_p_label[i], " & ",
                           data_rounded$Tgrade[i], " & ",
                           data_rounded$CATE_formatted[i], " & ",
                           data_rounded$CRATE_formatted[i], " \\\\\n"
      )
    }
  } else {

    county_name_display <- ifelse(is.na(data_rounded$county_name[1]),
                                  paste0("County ", data_rounded$County_ID[1]),
                                  data_rounded$county_name[1])
    latex_code <- paste0(latex_code, "\\begin{tabular}{cccccc}\n")
    latex_code <- paste0(latex_code, "\\hline\n")
    latex_code <- paste0(latex_code, "Age & BX & HR & Grade & CATE & CRATE \\\\\n")
    latex_code <- paste0(latex_code, " & Delay & Status &  & (Months) & (Months) \\\\\n")
    latex_code <- paste0(latex_code, "\\hline\n")

    for(i in 1:nrow(data_rounded)) {
      latex_code <- paste0(latex_code,
                           data_rounded$Age[i], " & ",
                           data_rounded$BX_Delay_label[i], " & ",
                           data_rounded$HR_p_label[i], " & ",
                           data_rounded$Tgrade[i], " & ",
                           data_rounded$CATE_formatted[i], " & ",
                           data_rounded$CRATE_formatted[i], " \\\\\n"
      )
    }
  }

  latex_code <- paste0(latex_code, "\\hline\n")
  latex_code <- paste0(latex_code, "\\end{tabular}\n")
  latex_code <- paste0(latex_code, "\\end{table}\n")

  return(latex_code)
}


big_county_name <- florida_counties$county_name[florida_counties$county_id == big_county_id]
small_county_name <- florida_counties$county_name[florida_counties$county_id == small_county_id]


latex_combined_names <- generate_latex_table_with_names(
  data = combined_summary_months_with_names,
  caption = "CATE and CRATE Estimates (Months) for All Covariate Combinations by County",
  label = "tab:cate_crate_combined_names",
  table_type = "combined"
)

# Big county only
latex_big_names <- generate_latex_table_with_names(
  data = big_county_summary_months_with_names,
  caption = paste0("CATE and CRATE Estimates (Months) for ", big_county_name, " County"),
  label = "tab:cate_crate_big_county_names",
  table_type = "single"
)

# Small county only
latex_small_names <- generate_latex_table_with_names(
  data = small_county_summary_months_with_names,
  caption = paste0("CATE and CRATE Estimates (Months) for ", small_county_name, " County"),
  label = "tab:cate_crate_small_county_names",
  table_type = "single"
)


summary_stats_names <- data.frame(
  County_Name = c(big_county_name, small_county_name),
  County_ID = c(big_county_id, small_county_id),
  N_Patients = c(sum(data_stratum$county == big_county_id),
                 sum(data_stratum$county == small_county_id)),
  CATE_Mean = c(mean(big_county_summary_months$CATE_mean),
                mean(small_county_summary_months$CATE_mean)),
  CATE_SD = c(sd(big_county_summary_months$CATE_mean),
              sd(small_county_summary_months$CATE_mean)),
  CRATE_Mean = c(mean(big_county_summary_months$CRATE_mean),
                 mean(small_county_summary_months$CRATE_mean)),
  CRATE_SD = c(sd(big_county_summary_months$CRATE_mean),
               sd(small_county_summary_months$CRATE_mean))
)

# Round summary statistics
summary_stats_names[, 4:7] <- round(summary_stats_names[, 4:7], 3)

latex_summary_stats_names <- "\\begin{table}[htbp]\n"
latex_summary_stats_names <- paste0(latex_summary_stats_names, "\\centering\n")
latex_summary_stats_names <- paste0(latex_summary_stats_names, "\\caption{Summary Statistics for CATE and CRATE by County (Months)}\n")
latex_summary_stats_names <- paste0(latex_summary_stats_names, "\\label{tab:summary_stats_names}\n")
latex_summary_stats_names <- paste0(latex_summary_stats_names, "\\begin{tabular}{lcccccc}\n")
latex_summary_stats_names <- paste0(latex_summary_stats_names, "\\hline\n")
latex_summary_stats_names <- paste0(latex_summary_stats_names, "County Name & N Patients & \\multicolumn{2}{c}{CATE (Months)} & \\multicolumn{2}{c}{CRATE (Months)} \\\\\n")
latex_summary_stats_names <- paste0(latex_summary_stats_names, " &  & Mean & SD & Mean & SD \\\\\n")
latex_summary_stats_names <- paste0(latex_summary_stats_names, "\\hline\n")

for(i in 1:nrow(summary_stats_names)) {
  latex_summary_stats_names <- paste0(latex_summary_stats_names,
                                      summary_stats_names$County_Name[i], " & ",
                                      summary_stats_names$N_Patients[i], " & ",
                                      summary_stats_names$CATE_Mean[i], " & ",
                                      summary_stats_names$CATE_SD[i], " & ",
                                      summary_stats_names$CRATE_Mean[i], " & ",
                                      summary_stats_names$CRATE_SD[i], " \\\\\n"
  )
}

latex_summary_stats_names <- paste0(latex_summary_stats_names, "\\hline\n")
latex_summary_stats_names <- paste0(latex_summary_stats_names, "\\end{tabular}\n")
latex_summary_stats_names <- paste0(latex_summary_stats_names, "\\end{table}\n")


cat("\n", paste(rep("=", 80), collapse=""), "\n")
cat("UPDATED RESULTS WITH COUNTY NAMES\n")
cat(paste(rep("=", 80), collapse=""), "\n")

cat("\nCounty Information:\n")
cat("Big County: ", big_county_name, " (ID: ", big_county_id, ") with ",
    sum(data_stratum$county == big_county_id), " patients\n")
cat("Small County: ", small_county_name, " (ID: ", small_county_id, ") with ",
    sum(data_stratum$county == small_county_id), " patients\n")

# Show first few rows with county names
cat("\nFirst 10 combinations with county names and readable labels:\n")
# Add readable labels to the display data
display_data <- combined_summary_months_with_names[, c("county_name", "Age", "BX_Delay", "HR_p", "Tgrade", "CATE_mean", "CRATE_mean")]
display_data$BX_Delay_label <- ifelse(display_data$BX_Delay == 1, "Yes", "No")
display_data$HR_p_label <- ifelse(display_data$HR_p == 1, "Negative", "Positive")
print(head(display_data[, c("county_name", "Age", "BX_Delay_label", "HR_p_label", "Tgrade", "CATE_mean", "CRATE_mean")], 20))



write.csv(combined_summary_months_with_names, "combined_summary_months_with_names.csv", row.names = FALSE)
write.csv(big_county_summary_months_with_names, "big_county_summary_months_with_names.csv", row.names = FALSE)
write.csv(small_county_summary_months_with_names, "small_county_summary_months_with_names.csv", row.names = FALSE)
write.csv(summary_stats_names, "summary_statistics_months_with_names.csv", row.names = FALSE)


writeLines(latex_combined_names, "latex_table_combined_with_names.tex")
writeLines(latex_big_names, "latex_table_big_county_with_names.tex")
writeLines(latex_small_names, "latex_table_small_county_with_names.tex")
writeLines(latex_summary_stats_names, "latex_summary_statistics_with_names.tex")


all_latex_names <- paste(
  "% Combined Table with County Names\n", latex_combined_names, "\n\n",
  "% Big County Table with County Names\n", latex_big_names, "\n\n",
  "% Small County Table with County Names\n", latex_small_names, "\n\n",
  "% Summary Statistics with County Names\n", latex_summary_stats_names, sep = ""
)
writeLines(all_latex_names, "all_latex_tables_with_names.tex")


cat("\n", paste(rep("=", 80), collapse=""), "\n")
cat("UPDATED LATEX TABLE CODE WITH COUNTY NAMES\n")
cat(paste(rep("=", 80), collapse=""), "\n")

cat("\n1. COMBINED TABLE (Both Counties with Names):\n")
cat(paste(rep("-", 50), collapse=""), "\n")
cat(latex_combined_names)

cat("\n2. BIG COUNTY TABLE (", big_county_name, "):\n")
cat(paste(rep("-", 50), collapse=""), "\n")
cat(latex_big_names)

cat("\n3. SMALL COUNTY TABLE (", small_county_name, "):\n")
cat(paste(rep("-", 50), collapse=""), "\n")
cat(latex_small_names)

cat("\n4. SUMMARY STATISTICS TABLE WITH COUNTY NAMES:\n")
cat(paste(rep("-", 50), collapse=""), "\n")
cat(latex_summary_stats_names)


cat("\n", paste(rep("=", 60), collapse=""), "\n")
cat("FILES SAVED WITH COUNTY NAMES\n")
cat(paste(rep("=", 60), collapse=""), "\n")

cat("CSV Files (with county names):\n")
cat("  - combined_summary_months_with_names.csv\n")
cat("  - big_county_summary_months_with_names.csv\n")
cat("  - small_county_summary_months_with_names.csv\n")
cat("  - summary_statistics_months_with_names.csv\n\n")

cat("LaTeX Files (with county names):\n")
cat("  - latex_table_combined_with_names.tex\n")
cat("  - latex_table_big_county_with_names.tex\n")
cat("  - latex_table_small_county_with_names.tex\n")
cat("  - latex_summary_statistics_with_names.tex\n")
cat("  - all_latex_tables_with_names.tex (all tables in one file)\n\n")

cat("Usage in LaTeX document:\n")
cat("\\input{latex_table_combined_with_names.tex}\n")
cat("\\input{latex_summary_statistics_with_names.tex}\n")

cat("\nCounty names successfully added to all tables!\n")
cat("Your tables now show actual Florida county names instead of numbers.\n")





t_star <- 5 * 365.25  # 5 years in days


compute_county_specific_effects <- function(
    individual_CATE_samples,     # M x n matrix of individual CATE samples
    individual_CRATE_samples,    # M x n matrix of individual CRATE samples
    patient_counties,            # vector of county IDs for each patient
    patient_treatments,          # vector of treatment assignments (Z) for each patient
    county_list = NULL           # optional: specific counties to analyze
) {

  # Get unique counties
  if(is.null(county_list)) {
    unique_counties <- sort(unique(patient_counties))
  } else {
    unique_counties <- county_list
  }

  n_counties <- length(unique_counties)
  n_samples <- nrow(individual_CATE_samples)

  cat("Computing county-specific effects for", n_counties, "counties...\n")

  # Initialize storage matrices
  county_results <- list(
    # ATE (Average Treatment Effect) - all patients in county
    ATE_samples = matrix(0, nrow = n_samples, ncol = n_counties),
    ATE_means = numeric(n_counties),
    ATE_lower = numeric(n_counties),
    ATE_upper = numeric(n_counties),

    # ATT (Average Treatment Effect on Treated) - only treated patients
    ATT_samples = matrix(0, nrow = n_samples, ncol = n_counties),
    ATT_means = numeric(n_counties),
    ATT_lower = numeric(n_counties),
    ATT_upper = numeric(n_counties),

    # ATU (Average Treatment Effect on Untreated) - only untreated patients
    ATU_samples = matrix(0, nrow = n_samples, ncol = n_counties),
    ATU_means = numeric(n_counties),
    ATU_lower = numeric(n_counties),
    ATU_upper = numeric(n_counties),

    # CRATE (Conditional Restricted ATE at 5 years) - all patients
    CRATE_samples = matrix(0, nrow = n_samples, ncol = n_counties),
    CRATE_means = numeric(n_counties),
    CRATE_lower = numeric(n_counties),
    CRATE_upper = numeric(n_counties),

    # CRATT (Conditional Restricted ATT at 5 years) - only treated patients
    CRATT_samples = matrix(0, nrow = n_samples, ncol = n_counties),
    CRATT_means = numeric(n_counties),
    CRATT_lower = numeric(n_counties),
    CRATT_upper = numeric(n_counties),

    # CRATU (Conditional Restricted ATU at 5 years) - only untreated patients
    CRATU_samples = matrix(0, nrow = n_samples, ncol = n_counties),
    CRATU_means = numeric(n_counties),
    CRATU_lower = numeric(n_counties),
    CRATU_upper = numeric(n_counties),

    # Sample sizes
    n_total = numeric(n_counties),      # Total patients in county
    n_treated = numeric(n_counties),    # Treated patients in county
    n_untreated = numeric(n_counties),  # Untreated patients in county


    county_ids = unique_counties
  )

  # Add column names
  colnames(county_results$ATE_samples) <- paste0("County_", unique_counties)
  colnames(county_results$ATT_samples) <- paste0("County_", unique_counties)
  colnames(county_results$ATU_samples) <- paste0("County_", unique_counties)
  colnames(county_results$CRATE_samples) <- paste0("County_", unique_counties)
  colnames(county_results$CRATT_samples) <- paste0("County_", unique_counties)
  colnames(county_results$CRATU_samples) <- paste0("County_", unique_counties)

  # Progress bar
  pb <- txtProgressBar(min = 0, max = n_counties, style = 3)

  # Compute effects for each county
  for(c in 1:n_counties) {
    county_id <- unique_counties[c]


    county_patients <- which(patient_counties == county_id)
    treated_patients <- which(patient_counties == county_id & patient_treatments == 1)
    untreated_patients <- which(patient_counties == county_id & patient_treatments == 0)

    # Sample sizes
    county_results$n_total[c] <- length(county_patients)
    county_results$n_treated[c] <- length(treated_patients)
    county_results$n_untreated[c] <- length(untreated_patients)

    # Check if we have patients in each subgroup
    if(length(county_patients) > 0) {
      # ATE: Average over all patients in county
      for(m in 1:n_samples) {
        county_results$ATE_samples[m, c] <- mean(individual_CATE_samples[m, county_patients])
        county_results$CRATE_samples[m, c] <- mean(individual_CRATE_samples[m, county_patients])
      }
      county_results$ATE_means[c] <- mean(county_results$ATE_samples[, c])
      county_results$ATE_lower[c] <- quantile(county_results$ATE_samples[, c], 0.025)
      county_results$ATE_upper[c] <- quantile(county_results$ATE_samples[, c], 0.975)

      county_results$CRATE_means[c] <- mean(county_results$CRATE_samples[, c])
      county_results$CRATE_lower[c] <- quantile(county_results$CRATE_samples[, c], 0.025)
      county_results$CRATE_upper[c] <- quantile(county_results$CRATE_samples[, c], 0.975)
    }

    if(length(treated_patients) > 0) {
      # ATT: Average over treated patients only
      for(m in 1:n_samples) {
        county_results$ATT_samples[m, c] <- mean(individual_CATE_samples[m, treated_patients])
        county_results$CRATT_samples[m, c] <- mean(individual_CRATE_samples[m, treated_patients])
      }
      county_results$ATT_means[c] <- mean(county_results$ATT_samples[, c])
      county_results$ATT_lower[c] <- quantile(county_results$ATT_samples[, c], 0.025)
      county_results$ATT_upper[c] <- quantile(county_results$ATT_samples[, c], 0.975)

      county_results$CRATT_means[c] <- mean(county_results$CRATT_samples[, c])
      county_results$CRATT_lower[c] <- quantile(county_results$CRATT_samples[, c], 0.025)
      county_results$CRATT_upper[c] <- quantile(county_results$CRATT_samples[, c], 0.975)
    } else {

      county_results$ATT_samples[, c] <- NA
      county_results$ATT_means[c] <- NA
      county_results$ATT_lower[c] <- NA
      county_results$ATT_upper[c] <- NA
      county_results$CRATT_samples[, c] <- NA
      county_results$CRATT_means[c] <- NA
      county_results$CRATT_lower[c] <- NA
      county_results$CRATT_upper[c] <- NA
    }

    if(length(untreated_patients) > 0) {
      # ATU: Average over untreated patients only
      for(m in 1:n_samples) {
        county_results$ATU_samples[m, c] <- mean(individual_CATE_samples[m, untreated_patients])
        county_results$CRATU_samples[m, c] <- mean(individual_CRATE_samples[m, untreated_patients])
      }
      county_results$ATU_means[c] <- mean(county_results$ATU_samples[, c])
      county_results$ATU_lower[c] <- quantile(county_results$ATU_samples[, c], 0.025)
      county_results$ATU_upper[c] <- quantile(county_results$ATU_samples[, c], 0.975)

      county_results$CRATU_means[c] <- mean(county_results$CRATU_samples[, c])
      county_results$CRATU_lower[c] <- quantile(county_results$CRATU_samples[, c], 0.025)
      county_results$CRATU_upper[c] <- quantile(county_results$CRATU_samples[, c], 0.975)
    } else {

      county_results$ATU_samples[, c] <- NA
      county_results$ATU_means[c] <- NA
      county_results$ATU_lower[c] <- NA
      county_results$ATU_upper[c] <- NA
      county_results$CRATU_samples[, c] <- NA
      county_results$CRATU_means[c] <- NA
      county_results$CRATU_lower[c] <- NA
      county_results$CRATU_upper[c] <- NA
    }

    setTxtProgressBar(pb, c)
  }
  close(pb)

  return(county_results)
}


cat("Computing county-specific effects from individual patient results...\n")

county_effects <- compute_county_specific_effects(
  individual_CATE_samples = CATE_results$samples,      # M x n matrix
  individual_CRATE_samples = CRATE_results$samples,    # M x n matrix
  patient_counties = data_stratum$county,              # county ID for each patient
  patient_treatments = Z,                              # treatment assignment for each patient
  county_list = stratum_counties                       # counties in our stratum
)


# Convert from days to months
days_per_month <- 30.44

# Create comprehensive summary matrix with formatted credible intervals
county_summary_matrix <- data.frame(
  County_ID = county_effects$county_ids,
  County_Name = florida_counties$county_name[match(county_effects$county_ids, florida_counties$county_id)],
  N_Total = county_effects$n_total,
  N_Treated = county_effects$n_treated,
  N_Untreated = county_effects$n_untreated,

  # Formatted effects with credible intervals (in months)
  ATE = ifelse(is.na(county_effects$ATE_means),
               "NA",
               paste0(round(county_effects$ATE_means / days_per_month, 3),
                      " (", round(county_effects$ATE_lower / days_per_month, 3),
                      ", ", round(county_effects$ATE_upper / days_per_month, 3), ")")),

  ATT = ifelse(is.na(county_effects$ATT_means),
               "NA",
               paste0(round(county_effects$ATT_means / days_per_month, 3),
                      " (", round(county_effects$ATT_lower / days_per_month, 3),
                      ", ", round(county_effects$ATT_upper / days_per_month, 3), ")")),

  ATU = ifelse(is.na(county_effects$ATU_means),
               "NA",
               paste0(round(county_effects$ATU_means / days_per_month, 3),
                      " (", round(county_effects$ATU_lower / days_per_month, 3),
                      ", ", round(county_effects$ATU_upper / days_per_month, 3), ")")),

  CRATE = ifelse(is.na(county_effects$CRATE_means),
                 "NA",
                 paste0(round(county_effects$CRATE_means / days_per_month, 3),
                        " (", round(county_effects$CRATE_lower / days_per_month, 3),
                        ", ", round(county_effects$CRATE_upper / days_per_month, 3), ")")),

  CRATT = ifelse(is.na(county_effects$CRATT_means),
                 "NA",
                 paste0(round(county_effects$CRATT_means / days_per_month, 3),
                        " (", round(county_effects$CRATT_lower / days_per_month, 3),
                        ", ", round(county_effects$CRATT_upper / days_per_month, 3), ")")),

  CRATU = ifelse(is.na(county_effects$CRATU_means),
                 "NA",
                 paste0(round(county_effects$CRATU_means / days_per_month, 3),
                        " (", round(county_effects$CRATU_lower / days_per_month, 3),
                        ", ", round(county_effects$CRATU_upper / days_per_month, 3), ")"))
)

county_summary_matrix_numeric <- data.frame(
  County_ID = county_effects$county_ids,
  County_Name = florida_counties$county_name[match(county_effects$county_ids, florida_counties$county_id)],
  N_Total = county_effects$n_total,
  N_Treated = county_effects$n_treated,
  N_Untreated = county_effects$n_untreated,

  # ATE (in months)
  ATE_Mean = round(county_effects$ATE_means / days_per_month, 3),
  ATE_Lower = round(county_effects$ATE_lower / days_per_month, 3),
  ATE_Upper = round(county_effects$ATE_upper / days_per_month, 3),

  # ATT (in months)
  ATT_Mean = round(county_effects$ATT_means / days_per_month, 3),
  ATT_Lower = round(county_effects$ATT_lower / days_per_month, 3),
  ATT_Upper = round(county_effects$ATT_upper / days_per_month, 3),

  # ATU (in months)
  ATU_Mean = round(county_effects$ATU_means / days_per_month, 3),
  ATU_Lower = round(county_effects$ATU_lower / days_per_month, 3),
  ATU_Upper = round(county_effects$ATU_upper / days_per_month, 3),

  # CRATE at 5 years (in months)
  CRATE_Mean = round(county_effects$CRATE_means / days_per_month, 3),
  CRATE_Lower = round(county_effects$CRATE_lower / days_per_month, 3),
  CRATE_Upper = round(county_effects$CRATE_upper / days_per_month, 3),

  # CRATT at 5 years (in months)
  CRATT_Mean = round(county_effects$CRATT_means / days_per_month, 3),
  CRATT_Lower = round(county_effects$CRATT_lower / days_per_month, 3),
  CRATT_Upper = round(county_effects$CRATT_upper / days_per_month, 3),

  # CRATU at 5 years (in months)
  CRATU_Mean = round(county_effects$CRATU_means / days_per_month, 3),
  CRATU_Lower = round(county_effects$CRATU_lower / days_per_month, 3),
  CRATU_Upper = round(county_effects$CRATU_upper / days_per_month, 3)
)


cat("\n", paste(rep("=", 80), collapse=""), "\n")
cat("COUNTY-SPECIFIC TREATMENT EFFECTS SUMMARY\n")
cat(paste(rep("=", 80), collapse=""), "\n")

cat("Computed for", nrow(county_summary_matrix), "counties in stratum\n")
cat("Race =", race_val, ", Stage =", stage_val, "\n\n")

# Show summary
print(county_summary_matrix)

# Show which counties have sufficient sample sizes
cat("\nCounties with both treated and untreated patients:\n")
balanced_counties <- county_summary_matrix[
  county_summary_matrix$N_Treated > 0 & county_summary_matrix$N_Untreated > 0,
  c("County_Name", "N_Total", "N_Treated", "N_Untreated")
]
print(balanced_counties)

# ------------------------------------------------------------------------------
# 5. Save results
# ------------------------------------------------------------------------------

cat("\nSaving results...\n")

# Save summary matrix (formatted version)
write.csv(county_summary_matrix, "county_specific_effects_summary_formatted.csv", row.names = FALSE)

# Save numerical version for further analysis
write.csv(county_summary_matrix_numeric, "county_specific_effects_summary_numeric.csv", row.names = FALSE)

# Save posterior samples (convert to months)
county_samples_months <- list(
  ATE_samples = county_effects$ATE_samples / days_per_month,
  ATT_samples = county_effects$ATT_samples / days_per_month,
  ATU_samples = county_effects$ATU_samples / days_per_month,
  CRATE_samples = county_effects$CRATE_samples / days_per_month,
  CRATT_samples = county_effects$CRATT_samples / days_per_month,
  CRATU_samples = county_effects$CRATU_samples / days_per_month,
  county_ids = county_effects$county_ids,
  sample_sizes = data.frame(
    County_ID = county_effects$county_ids,
    N_Total = county_effects$n_total,
    N_Treated = county_effects$n_treated,
    N_Untreated = county_effects$n_untreated
  )
)

# Save individual sample matrices
write.csv(county_samples_months$ATE_samples, "county_ATE_samples_months.csv", row.names = FALSE)
write.csv(county_samples_months$ATT_samples, "county_ATT_samples_months.csv", row.names = FALSE)
write.csv(county_samples_months$ATU_samples, "county_ATU_samples_months.csv", row.names = FALSE)
write.csv(county_samples_months$CRATE_samples, "county_CRATE_samples_months.csv", row.names = FALSE)
write.csv(county_samples_months$CRATT_samples, "county_CRATT_samples_months.csv", row.names = FALSE)
write.csv(county_samples_months$CRATU_samples, "county_CRATU_samples_months.csv", row.names = FALSE)

# Save complete results
save(county_effects, county_summary_matrix, county_summary_matrix_numeric, county_samples_months,
     file = "county_specific_effects_complete.RData")

cat("Files saved:\n")
cat("  - county_specific_effects_summary_formatted.csv (formatted with CIs)\n")
cat("  - county_specific_effects_summary_numeric.csv (numerical values)\n")
cat("  - county_ATE_samples_months.csv (", nrow(county_samples_months$ATE_samples), " x ", ncol(county_samples_months$ATE_samples), ")\n")
cat("  - county_ATT_samples_months.csv (", nrow(county_samples_months$ATT_samples), " x ", ncol(county_samples_months$ATT_samples), ")\n")
cat("  - county_ATU_samples_months.csv (", nrow(county_samples_months$ATU_samples), " x ", ncol(county_samples_months$ATU_samples), ")\n")
cat("  - county_CRATE_samples_months.csv (", nrow(county_samples_months$CRATE_samples), " x ", ncol(county_samples_months$CRATE_samples), ")\n")
cat("  - county_CRATT_samples_months.csv (", nrow(county_samples_months$CRATT_samples), " x ", ncol(county_samples_months$CRATT_samples), ")\n")
cat("  - county_CRATU_samples_months.csv (", nrow(county_samples_months$CRATU_samples), " x ", ncol(county_samples_months$CRATU_samples), ")\n")
cat("  - county_specific_effects_complete.RData (all R objects)\n")

# ------------------------------------------------------------------------------
# 7. Generate LaTeX table for county-specific effects
# ------------------------------------------------------------------------------

generate_county_latex_table <- function(data, caption, label) {

  # Start LaTeX table
  latex_code <- "\\begin{table}[htbp]\n"
  latex_code <- paste0(latex_code, "\\centering\n")
  latex_code <- paste0(latex_code, "\\caption{", caption, "}\n")
  latex_code <- paste0(latex_code, "\\label{", label, "}\n")
  latex_code <- paste0(latex_code, "\\begin{tabular}{lcccccccc}\n")
  latex_code <- paste0(latex_code, "\\hline\n")
  latex_code <- paste0(latex_code, "County & \\multicolumn{3}{c}{Sample Size} & \\multicolumn{3}{c}{Treatment Effects (Months)} & \\multicolumn{3}{c}{Restricted Effects (Months)} \\\\\n")
  latex_code <- paste0(latex_code, "Name & Total & Treated & Untreated & ATE & ATT & ATU & CRATE & CRATT & CRATU \\\\\n")
  latex_code <- paste0(latex_code, "\\hline\n")

  for(i in 1:nrow(data)) {
    # Handle county name (replace any special characters for LaTeX)
    county_name <- gsub("&", "\\&", data$County_Name[i])
    county_name <- gsub("_", "\\_", county_name)

    latex_code <- paste0(latex_code,
                         county_name, " & ",
                         data$N_Total[i], " & ",
                         data$N_Treated[i], " & ",
                         data$N_Untreated[i], " & ",
                         data$ATE[i], " & ",
                         data$ATT[i], " & ",
                         data$ATU[i], " & ",
                         data$CRATE[i], " & ",
                         data$CRATT[i], " & ",
                         data$CRATU[i], " \\\\\n"
    )
  }

  latex_code <- paste0(latex_code, "\\hline\n")
  latex_code <- paste0(latex_code, "\\end{tabular}\n")
  latex_code <- paste0(latex_code, "\\end{table}\n")

  return(latex_code)
}

# Generate the LaTeX table
county_latex_table <- generate_county_latex_table(
  data = county_summary_matrix,
  caption = paste0("County-Specific Treatment Effects for Race = ", race_val,
                   " and Stage = ", stage_val, " (Months with 95\\% Credible Intervals)"),
  label = "tab:county_specific_effects"
)

# Generate a more compact version without CRATT and CRATU for better fit
generate_county_latex_table_compact <- function(data, caption, label) {

  # Start LaTeX table
  latex_code <- "\\begin{table}[htbp]\n"
  latex_code <- paste0(latex_code, "\\centering\n")
  latex_code <- paste0(latex_code, "\\caption{", caption, "}\n")
  latex_code <- paste0(latex_code, "\\label{", label, "}\n")
  latex_code <- paste0(latex_code, "\\begin{tabular}{lcccccc}\n")
  latex_code <- paste0(latex_code, "\\hline\n")
  latex_code <- paste0(latex_code, "County & \\multicolumn{3}{c}{Sample Size} & \\multicolumn{3}{c}{Treatment Effects (Months)} \\\\\n")
  latex_code <- paste0(latex_code, "Name & Total & Treated & Untreated & ATE & ATT & ATU \\\\\n")
  latex_code <- paste0(latex_code, "\\hline\n")

  for(i in 1:nrow(data)) {
    # Handle county name (replace any special characters for LaTeX)
    county_name <- gsub("&", "\\&", data$County_Name[i])
    county_name <- gsub("_", "\\_", county_name)

    latex_code <- paste0(latex_code,
                         county_name, " & ",
                         data$N_Total[i], " & ",
                         data$N_Treated[i], " & ",
                         data$N_Untreated[i], " & ",
                         data$ATE[i], " & ",
                         data$ATT[i], " & ",
                         data$ATU[i], " \\\\\n"
    )
  }

  latex_code <- paste0(latex_code, "\\hline\n")
  latex_code <- paste0(latex_code, "\\end{tabular}\n")
  latex_code <- paste0(latex_code, "\\end{table}\n")

  return(latex_code)
}

# Generate restricted effects table separately
generate_county_latex_table_restricted <- function(data, caption, label) {

  # Start LaTeX table
  latex_code <- "\\begin{table}[htbp]\n"
  latex_code <- paste0(latex_code, "\\centering\n")
  latex_code <- paste0(latex_code, "\\caption{", caption, "}\n")
  latex_code <- paste0(latex_code, "\\label{", label, "}\n")
  latex_code <- paste0(latex_code, "\\begin{tabular}{lcccc}\n")
  latex_code <- paste0(latex_code, "\\hline\n")
  latex_code <- paste0(latex_code, "County & \\multicolumn{3}{c}{Restricted Treatment Effects at 5 Years (Months)} \\\\\n")
  latex_code <- paste0(latex_code, "Name & CRATE & CRATT & CRATU \\\\\n")
  latex_code <- paste0(latex_code, "\\hline\n")

  for(i in 1:nrow(data)) {
    # Handle county name (replace any special characters for LaTeX)
    county_name <- gsub("&", "\\&", data$County_Name[i])
    county_name <- gsub("_", "\\_", county_name)

    latex_code <- paste0(latex_code,
                         county_name, " & ",
                         data$CRATE[i], " & ",
                         data$CRATT[i], " & ",
                         data$CRATU[i], " \\\\\n"
    )
  }

  latex_code <- paste0(latex_code, "\\hline\n")
  latex_code <- paste0(latex_code, "\\end{tabular}\n")
  latex_code <- paste0(latex_code, "\\end{table}\n")

  return(latex_code)
}

# Generate compact and restricted tables
county_latex_table_compact <- generate_county_latex_table_compact(
  data = county_summary_matrix,
  caption = paste0("County-Specific Treatment Effects for Race = ", race_val,
                   " and Stage = ", stage_val, " (Months with 95\\% Credible Intervals)"),
  label = "tab:county_effects_compact"
)

county_latex_table_restricted <- generate_county_latex_table_restricted(
  data = county_summary_matrix,
  caption = paste0("County-Specific Restricted Treatment Effects at 5 Years for Race = ", race_val,
                   " and Stage = ", stage_val, " (Months with 95\\% Credible Intervals)"),
  label = "tab:county_effects_restricted"
)

# ------------------------------------------------------------------------------
# 8. Display and Save LaTeX Tables
# ------------------------------------------------------------------------------

cat("\n", paste(rep("=", 80), collapse=""), "\n")
cat("LATEX TABLES FOR COUNTY-SPECIFIC EFFECTS\n")
cat(paste(rep("=", 80), collapse=""), "\n")

cat("\n1. COMPLETE TABLE (All Effects):\n")
cat(paste(rep("-", 50), collapse=""), "\n")
cat(county_latex_table)

cat("\n2. COMPACT TABLE (ATE, ATT, ATU only):\n")
cat(paste(rep("-", 50), collapse=""), "\n")
cat(county_latex_table_compact)

cat("\n3. RESTRICTED EFFECTS TABLE (CRATE, CRATT, CRATU):\n")
cat(paste(rep("-", 50), collapse=""), "\n")
cat(county_latex_table_restricted)

# Save LaTeX tables to files
writeLines(county_latex_table, "county_effects_complete_table.tex")
writeLines(county_latex_table_compact, "county_effects_compact_table.tex")
writeLines(county_latex_table_restricted, "county_effects_restricted_table.tex")

# Save all LaTeX tables to a single file
all_county_latex <- paste(
  "% Complete County Effects Table\n", county_latex_table, "\n\n",
  "% Compact County Effects Table\n", county_latex_table_compact, "\n\n",
  "% Restricted County Effects Table\n", county_latex_table_restricted, sep = ""
)
writeLines(all_county_latex, "all_county_latex_tables.tex")

cat("\nLaTeX files saved:\n")
cat("  - county_effects_complete_table.tex (all effects)\n")
cat("  - county_effects_compact_table.tex (ATE, ATT, ATU only)\n")
cat("  - county_effects_restricted_table.tex (CRATE, CRATT, CRATU only)\n")
cat("  - all_county_latex_tables.tex (all tables combined)\n")

if(nrow(county_summary_matrix) > 0) {
  cat("\n", paste(rep("=", 60), collapse=""), "\n")
  cat("EXAMPLE: Results for", county_summary_matrix$County_Name[1], "County\n")
  cat(paste(rep("=", 60), collapse=""), "\n")

  example_county <- 1
  cat("Sample sizes:\n")
  cat("  Total patients:", county_summary_matrix$N_Total[example_county], "\n")
  cat("  Treated patients:", county_summary_matrix$N_Treated[example_county], "\n")
  cat("  Untreated patients:", county_summary_matrix$N_Untreated[example_county], "\n\n")

  cat("Treatment effects (months):\n")
  cat("  ATE:", county_summary_matrix$ATE[example_county], "\n")

  if(!is.na(county_summary_matrix_numeric$ATT_Mean[example_county])) {
    cat("  ATT:", county_summary_matrix$ATT[example_county], "\n")
  }

  if(!is.na(county_summary_matrix_numeric$ATU_Mean[example_county])) {
    cat("  ATU:", county_summary_matrix$ATU[example_county], "\n")
  }

  cat("  CRATE (5yr):", county_summary_matrix$CRATE[example_county], "\n")
}

cat("\nCounty-specific analysis complete!\n")
cat("You now have ATE, ATT, ATU, CRATE, CRATT, and CRATU for each county\n")
cat("with full posterior samples and credible intervals.\n")






# Load required libraries
library(ggplot2)
library(sf)
library(tigris)
library(dplyr)
library(viridis)
library(cowplot)
library(gridExtra)

# Set options for tigris
options(tigris_use_cache = TRUE)

cat("Creating Florida county maps with legends from scratch...\n")

# ------------------------------------------------------------------------------
# 1. Download and prepare Florida county shapefiles
# ------------------------------------------------------------------------------

cat("Downloading Florida county shapefiles...\n")

# Download Florida counties shapefile
florida_counties_sf <- counties(state = "FL", cb = TRUE, year = 2021)

# Clean county names and add FIPS codes
florida_counties_sf <- florida_counties_sf %>%
  mutate(
    county_name_clean = gsub(" County", "", NAME),
    county_fips = as.numeric(COUNTYFP)
  )


florida_fips_mapping <- data.frame(
  county_id = 1:67,
  county_name = c(
    "Alachua", "Baker", "Bay", "Bradford", "Brevard", "Broward", "Calhoun",
    "Charlotte", "Citrus", "Clay", "Collier", "Columbia", "DeSoto", "Dixie",
    "Duval", "Escambia", "Flagler", "Franklin", "Gadsden", "Gilchrist",
    "Glades", "Gulf", "Hamilton", "Hardee", "Hendry", "Hernando", "Highlands",
    "Hillsborough", "Holmes", "Indian River", "Jackson", "Jefferson", "Lafayette",
    "Lake", "Lee", "Leon", "Levy", "Liberty", "Madison", "Manatee", "Marion",
    "Martin", "Miami-Dade", "Monroe", "Nassau", "Okaloosa", "Okeechobee",
    "Orange", "Osceola", "Palm Beach", "Pasco", "Pinellas", "Polk", "Putnam",
    "Santa Rosa", "Sarasota", "Seminole", "St. Johns", "St. Lucie", "Sumter",
    "Suwannee", "Taylor", "Union", "Volusia", "Wakulla", "Walton", "Washington"
  ),
  fips_code = c(
    1, 3, 5, 7, 9, 11, 13, 15, 17, 19, 21, 23, 27, 29, 31, 33, 35, 37, 39, 41,
    43, 45, 47, 49, 51, 53, 55, 57, 59, 61, 63, 65, 67, 69, 71, 73, 75, 77, 79,
    81, 83, 85, 86, 87, 89, 91, 93, 95, 97, 99, 101, 103, 105, 107, 109, 111,
    113, 115, 117, 119, 121, 123, 125, 127, 129, 131, 133
  )
)

# ------------------------------------------------------------------------------
# 2. Merge treatment effects data with county mapping and shapefile
# ------------------------------------------------------------------------------

# Merge county results with FIPS mapping
county_data_with_fips <- county_summary_matrix_numeric %>%
  left_join(florida_fips_mapping, by = c("County_ID" = "county_id"))

# Merge with shapefile
florida_map_data <- florida_counties_sf %>%
  left_join(county_data_with_fips, by = c("county_fips" = "fips_code"))

cat("Merged data for", nrow(florida_map_data), "counties\n")
cat("Counties with treatment effect data:", sum(!is.na(florida_map_data$ATE_Mean)), "\n")


# Calculate overall range from all treatment effects for consistency
all_effects <- c(
  florida_map_data$ATE_Mean,
  florida_map_data$ATT_Mean,
  florida_map_data$ATU_Mean,
  florida_map_data$CRATE_Mean,
  florida_map_data$CRATT_Mean,
  florida_map_data$CRATU_Mean
)


overall_range <- range(all_effects, na.rm = TRUE)
cat("Using consistent color scale range: [", round(overall_range[1], 3), ", ", round(overall_range[2], 3), "] months\n")

# Create breaks for the color scale
color_breaks <- pretty(overall_range, n = 5)


create_county_map <- function(map_data, effect_column, map_title = "",
                              color_range = overall_range, show_legend = TRUE) {

  p <- ggplot(map_data) +
    geom_sf(aes_string(fill = effect_column),
            color = "white",
            size = 0.2) +
    scale_fill_viridis_c(
      name = "Treatment Effect\n(Months)",
      na.value = "lightgray",
      limits = color_range,
      breaks = color_breaks,
      labels = round(color_breaks, 2),
      option = "viridis"
    ) +
    theme_void() +
    theme(
      plot.background = element_rect(fill = "white", color = NA),
      panel.background = element_rect(fill = "white", color = NA),
      legend.position = if(show_legend) "bottom" else "none",
      legend.title = element_text(size = 12, face = "bold"),
      legend.text = element_text(size = 10),
      legend.key.width = unit(2, "cm"),
      legend.key.height = unit(0.5, "cm"),
      legend.margin = margin(t = 10, b = 5),
      plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
      plot.margin = margin(10, 10, 10, 10)
    )

  if(map_title != "") {
    p <- p + ggtitle(map_title)
  }

  return(p)
}



cat("Creating Figure 1: ATE and CRATE...\n")

# Create individual maps without any titles
map_ATE <- create_county_map(florida_map_data, "ATE_Mean", "", overall_range, FALSE)
map_CRATE <- create_county_map(florida_map_data, "CRATE_Mean", "", overall_range, FALSE)

# Create a separate plot just for the legend
legend_plot <- create_county_map(florida_map_data, "ATE_Mean", "", overall_range, TRUE)
shared_legend <- get_legend(legend_plot)

# Combine maps with legend - NO LABELS OR TITLES
figure1 <- plot_grid(
  plot_grid(map_ATE, map_CRATE, ncol = 2),
  shared_legend,
  ncol = 1,
  rel_heights = c(1, 0.15)
)

# ------------------------------------------------------------------------------
# 6. Create Figure 2: ATT and CRATT (No titles/headings)
# ------------------------------------------------------------------------------

cat("Creating Figure 2: ATT and CRATT...\n")

# Create individual maps without any titles
map_ATT <- create_county_map(florida_map_data, "ATT_Mean", "", overall_range, FALSE)
map_CRATT <- create_county_map(florida_map_data, "CRATT_Mean", "", overall_range, FALSE)

# Combine maps with legend - NO LABELS OR TITLES
figure2 <- plot_grid(
  plot_grid(map_ATT, map_CRATT, ncol = 2),
  shared_legend,
  ncol = 1,
  rel_heights = c(1, 0.15)
)

# ------------------------------------------------------------------------------
# 7. Create Figure 3: ATU and CRATU (No titles/headings)
# ------------------------------------------------------------------------------

cat("Creating Figure 3: ATU and CRATU...\n")

# Create individual maps without any titles
map_ATU <- create_county_map(florida_map_data, "ATU_Mean", "", overall_range, FALSE)
map_CRATU <- create_county_map(florida_map_data, "CRATU_Mean", "", overall_range, FALSE)

# Combine maps with legend - NO LABELS OR TITLES
figure3 <- plot_grid(
  plot_grid(map_ATU, map_CRATU, ncol = 2),
  shared_legend,
  ncol = 1,
  rel_heights = c(1, 0.15)
)

# ------------------------------------------------------------------------------
# 8. Create clean simple versions (No titles/headings)
# ------------------------------------------------------------------------------

cat("Creating clean simple figures with individual legends...\n")

# Figure 1 - ATE with legend, CRATE without legend, no titles
map_ATE_legend <- create_county_map(florida_map_data, "ATE_Mean", "", overall_range, TRUE)
map_CRATE_no_legend <- create_county_map(florida_map_data, "CRATE_Mean", "", overall_range, FALSE)

figure1_simple <- plot_grid(map_ATE_legend, map_CRATE_no_legend, ncol = 2)

# Figure 2 - ATT with legend, CRATT without legend, no titles
map_ATT_legend <- create_county_map(florida_map_data, "ATT_Mean", "", overall_range, TRUE)
map_CRATT_no_legend <- create_county_map(florida_map_data, "CRATT_Mean", "", overall_range, FALSE)

figure2_simple <- plot_grid(map_ATT_legend, map_CRATT_no_legend, ncol = 2)

# Figure 3 - ATU with legend, CRATU without legend, no titles
map_ATU_legend <- create_county_map(florida_map_data, "ATU_Mean", "", overall_range, TRUE)
map_CRATU_no_legend <- create_county_map(florida_map_data, "CRATU_Mean", "", overall_range, FALSE)

figure3_simple <- plot_grid(map_ATU_legend, map_CRATU_no_legend, ncol = 2)

# ------------------------------------------------------------------------------
# 9. Save all figure versions
# ------------------------------------------------------------------------------

cat("Saving figures...\n")


ggsave("Figure1_ATE_CRATE_clean.png", figure1,
       width = 14, height = 8, dpi = 300, bg = "white")
ggsave("Figure1_ATE_CRATE_clean.pdf", figure1,
       width = 14, height = 8, bg = "white")

ggsave("Figure2_ATT_CRATT_clean.png", figure2,
       width = 14, height = 8, dpi = 300, bg = "white")
ggsave("Figure2_ATT_CRATT_clean.pdf", figure2,
       width = 14, height = 8, bg = "white")

ggsave("Figure3_ATU_CRATU_clean.png", figure3,
       width = 14, height = 8, dpi = 300, bg = "white")
ggsave("Figure3_ATU_CRATU_clean.pdf", figure3,
       width = 14, height = 8, bg = "white")

# Save the simple clean figures (legend on left map only)
ggsave("Figure1_ATE_CRATE_simple_clean.png", figure1_simple,
       width = 14, height = 6, dpi = 300, bg = "white")
ggsave("Figure1_ATE_CRATE_simple_clean.pdf", figure1_simple,
       width = 14, height = 6, bg = "white")

ggsave("Figure2_ATT_CRATT_simple_clean.png", figure2_simple,
       width = 14, height = 6, dpi = 300, bg = "white")
ggsave("Figure2_ATT_CRATT_simple_clean.pdf", figure2_simple,
       width = 14, height = 6, bg = "white")

ggsave("Figure3_ATU_CRATU_simple_clean.png", figure3_simple,
       width = 14, height = 6, dpi = 300, bg = "white")
ggsave("Figure3_ATU_CRATU_simple_clean.pdf", figure3_simple,
       width = 14, height = 6, bg = "white")


cat("\n", paste(rep("=", 80), collapse=""), "\n")
cat("FLORIDA COUNTY MAPS WITH LEGENDS - SUMMARY\n")
cat(paste(rep("=", 80), collapse=""), "\n")

cat("Created three figures with county-specific treatment effects:\n\n")

cat("Figure 1: ATE and CRATE\n")
cat("  - Left map: Average Treatment Effect\n")
cat("  - Right map: Conditional Restricted Average Treatment Effect (5 years)\n")
cat("  - No titles or headings, legend only\n")
cat("  - Files: Figure1_ATE_CRATE_clean.png/.pdf\n\n")

cat("Figure 2: ATT and CRATT\n")
cat("  - Left map: Average Treatment Effect on Treated\n")
cat("  - Right map: Conditional Restricted Average Treatment Effect on Treated (5 years)\n")
cat("  - No titles or headings, legend only\n")
cat("  - Files: Figure2_ATT_CRATT_clean.png/.pdf\n\n")

cat("Figure 3: ATU and CRATU\n")
cat("  - Left map: Average Treatment Effect on Untreated\n")
cat("  - Right map: Conditional Restricted Average Treatment Effect on Untreated (5 years)\n")
cat("  - No titles or headings, legend only\n")
cat("  - Files: Figure3_ATU_CRATU_clean.png/.pdf\n\n")

cat("Map specifications:\n")
cat("  - Consistent color scale across all figures: [", round(overall_range[1], 3), ", ", round(overall_range[2], 3), "] months\n")
cat("  - Viridis color palette (blue to yellow)\n")
cat("  - Counties with no data: light gray\n")
cat("  - Counties with data:", sum(!is.na(florida_map_data$ATE_Mean)), "out of 67\n")
cat("  - White background, no grids\n")
cat("  - Resolution: 300 DPI\n")
cat("  - Two versions: shared legend at bottom, and simple with individual legends\n")

# Test that individual maps work
cat("\nTesting individual map creation...\n")
test_map <- create_county_map(florida_map_data, "ATE_Mean", "Test ATE Map", overall_range, TRUE)
ggsave("test_individual_map.png", test_map, width = 8, height = 6, dpi = 300, bg = "white")
cat("Test map saved as test_individual_map.png\n")

# Save all objects
save(figure1, figure2, figure3, figure1_simple, figure2_simple, figure3_simple,
     florida_map_data, overall_range,
     file = "florida_maps_with_legends_complete.RData")

cat("\nAll figures created successfully with guaranteed legends!\n")
cat("If legends still don't appear, check the simple versions which have individual legends.\n")


























# Conversion factor from days to months
days_per_month <- 30.44


cat("Computing population-level treatment effects...\n")

# Check if we have individual patient results
if(!exists("CATE_results") || !exists("CRATE_results")) {
  stop("CATE_results and CRATE_results objects not found. Please run individual patient analysis first.")
}

# Get the individual patient samples (M x n matrices)
individual_CATE_samples <- CATE_results$samples    # M x n matrix
individual_CRATE_samples <- CRATE_results$samples  # M x n matrix

n_patients <- ncol(individual_CATE_samples)
n_samples <- nrow(individual_CATE_samples)

cat("Using", n_patients, "patients and", n_samples, "posterior samples\n")
cat("Patient treatment assignments:\n")
table(Z)


# ATE: Average over all patients for each posterior sample
population_ATE_samples <- rowMeans(individual_CATE_samples)
population_ATE_mean <- mean(population_ATE_samples)
population_ATE_lower <- quantile(population_ATE_samples, 0.025)
population_ATE_upper <- quantile(population_ATE_samples, 0.975)

# RATE: Average over all patients for each posterior sample
population_RATE_samples <- rowMeans(individual_CRATE_samples)
population_RATE_mean <- mean(population_RATE_samples)
population_RATE_lower <- quantile(population_RATE_samples, 0.025)
population_RATE_upper <- quantile(population_RATE_samples, 0.975)

# Get indices of treated patients
treated_indices <- which(Z == 1)
n_treated <- length(treated_indices)

if(n_treated > 0) {
  # ATT: Average over treated patients only for each posterior sample
  population_ATT_samples <- rowMeans(individual_CATE_samples[, treated_indices, drop = FALSE])
  population_ATT_mean <- mean(population_ATT_samples)
  population_ATT_lower <- quantile(population_ATT_samples, 0.025)
  population_ATT_upper <- quantile(population_ATT_samples, 0.975)

  # RATT: Average over treated patients only for each posterior sample
  population_RATT_samples <- rowMeans(individual_CRATE_samples[, treated_indices, drop = FALSE])
  population_RATT_mean <- mean(population_RATT_samples)
  population_RATT_lower <- quantile(population_RATT_samples, 0.025)
  population_RATT_upper <- quantile(population_RATT_samples, 0.975)
} else {
  population_ATT_samples <- rep(NA, n_samples)
  population_ATT_mean <- NA
  population_ATT_lower <- NA
  population_ATT_upper <- NA
  population_RATT_samples <- rep(NA, n_samples)
  population_RATT_mean <- NA
  population_RATT_lower <- NA
  population_RATT_upper <- NA
}


# Get indices of untreated patients
untreated_indices <- which(Z == 0)
n_untreated <- length(untreated_indices)

if(n_untreated > 0) {
  # ATU: Average over untreated patients only for each posterior sample
  population_ATU_samples <- rowMeans(individual_CATE_samples[, untreated_indices, drop = FALSE])
  population_ATU_mean <- mean(population_ATU_samples)
  population_ATU_lower <- quantile(population_ATU_samples, 0.025)
  population_ATU_upper <- quantile(population_ATU_samples, 0.975)

  # RATU: Average over untreated patients only for each posterior sample
  population_RATU_samples <- rowMeans(individual_CRATE_samples[, untreated_indices, drop = FALSE])
  population_RATU_mean <- mean(population_RATU_samples)
  population_RATU_lower <- quantile(population_RATU_samples, 0.025)
  population_RATU_upper <- quantile(population_RATU_samples, 0.975)
} else {
  population_ATU_samples <- rep(NA, n_samples)
  population_ATU_mean <- NA
  population_ATU_lower <- NA
  population_ATU_upper <- NA
  population_RATU_samples <- rep(NA, n_samples)
  population_RATU_mean <- NA
  population_RATU_lower <- NA
  population_RATU_upper <- NA
}

# ------------------------------------------------------------------------------
# 5. Create summary table (in months)
# ------------------------------------------------------------------------------

population_effects_summary <- data.frame(
  Effect = c("ATE", "ATT", "ATU", "RATE", "RATT", "RATU"),
  Estimand = c(
    "Average Treatment Effect",
    "Average Treatment Effect on Treated",
    "Average Treatment Effect on Untreated",
    "Restricted Average Treatment Effect (5 years)",
    "Restricted Average Treatment Effect on Treated (5 years)",
    "Restricted Average Treatment Effect on Untreated (5 years)"
  ),
  Sample_Size = c(
    n_patients, n_treated, n_untreated,
    n_patients, n_treated, n_untreated
  ),
  Mean_Days = c(
    population_ATE_mean, population_ATT_mean, population_ATU_mean,
    population_RATE_mean, population_RATT_mean, population_RATU_mean
  ),
  Lower_Days = c(
    population_ATE_lower, population_ATT_lower, population_ATU_lower,
    population_RATE_lower, population_RATT_lower, population_RATU_lower
  ),
  Upper_Days = c(
    population_ATE_upper, population_ATT_upper, population_ATU_upper,
    population_RATE_upper, population_RATT_upper, population_RATU_upper
  ),
  Mean_Months = round(c(
    population_ATE_mean, population_ATT_mean, population_ATU_mean,
    population_RATE_mean, population_RATT_mean, population_RATU_mean
  ) / days_per_month, 3),
  Lower_Months = round(c(
    population_ATE_lower, population_ATT_lower, population_ATU_lower,
    population_RATE_lower, population_RATT_lower, population_RATU_lower
  ) / days_per_month, 3),
  Upper_Months = round(c(
    population_ATE_upper, population_ATT_upper, population_ATU_upper,
    population_RATE_upper, population_RATT_upper, population_RATU_upper
  ) / days_per_month, 3)
)

# Add formatted column with mean (credible interval)
population_effects_summary$Effect_Formatted <- ifelse(
  is.na(population_effects_summary$Mean_Months),
  "NA",
  paste0(population_effects_summary$Mean_Months, " (",
         population_effects_summary$Lower_Months, ", ",
         population_effects_summary$Upper_Months, ")")
)

# ------------------------------------------------------------------------------
# 6. Display results
# ------------------------------------------------------------------------------

cat("\n", paste(rep("=", 80), collapse=""), "\n")
cat("POPULATION-LEVEL TREATMENT EFFECTS\n")
cat(paste(rep("=", 80), collapse=""), "\n")

cat("Study population: Race =", race_val, ", Stage =", stage_val, "\n")
cat("Total patients:", n_patients, "\n")
cat("Treated patients (Z=1):", n_treated, "\n")
cat("Untreated patients (Z=0):", n_untreated, "\n")
cat("Posterior samples:", n_samples, "\n\n")

cat("Population-Level Treatment Effects (Months with 95% Credible Intervals):\n\n")
for(i in 1:nrow(population_effects_summary)) {
  cat(sprintf("%-6s: %s\n",
              population_effects_summary$Effect[i],
              population_effects_summary$Effect_Formatted[i]))
}

# Show detailed table
cat("\nDetailed Summary Table:\n")
print(population_effects_summary[, c("Effect", "Sample_Size", "Mean_Months", "Lower_Months", "Upper_Months", "Effect_Formatted")])


generate_population_effects_latex <- function(data, caption, label) {

  # Start LaTeX table
  latex_code <- "\\begin{table}[htbp]\n"
  latex_code <- paste0(latex_code, "\\centering\n")
  latex_code <- paste0(latex_code, "\\caption{", caption, "}\n")
  latex_code <- paste0(latex_code, "\\label{", label, "}\n")
  latex_code <- paste0(latex_code, "\\begin{tabular}{lcc}\n")
  latex_code <- paste0(latex_code, "\\hline\n")
  latex_code <- paste0(latex_code, "Treatment Effect & Sample Size & Effect (95\\% CI) \\\\\n")
  latex_code <- paste0(latex_code, " &  & (Months) \\\\\n")
  latex_code <- paste0(latex_code, "\\hline\n")

  for(i in 1:nrow(data)) {
    # Handle NA values
    effect_display <- if(is.na(data$Mean_Months[i])) "NA" else data$Effect_Formatted[i]

    latex_code <- paste0(latex_code,
                         data$Effect[i], " & ",
                         data$Sample_Size[i], " & ",
                         effect_display, " \\\\\n"
    )
  }

  latex_code <- paste0(latex_code, "\\hline\n")

  # Add footnote explaining effects
  latex_code <- paste0(latex_code, "\\end{tabular}\n")
  latex_code <- paste0(latex_code, "\\begin{tablenotes}\n")
  latex_code <- paste0(latex_code, "\\small\n")
  latex_code <- paste0(latex_code, "\\item ATE: Average Treatment Effect (all patients)\n")
  latex_code <- paste0(latex_code, "\\item ATT: Average Treatment Effect on Treated\n")
  latex_code <- paste0(latex_code, "\\item ATU: Average Treatment Effect on Untreated\n")
  latex_code <- paste0(latex_code, "\\item RATE: Restricted Average Treatment Effect at 5 years\n")
  latex_code <- paste0(latex_code, "\\item RATT: Restricted Average Treatment Effect on Treated at 5 years\n")
  latex_code <- paste0(latex_code, "\\item RATU: Restricted Average Treatment Effect on Untreated at 5 years\n")
  latex_code <- paste0(latex_code, "\\end{tablenotes}\n")
  latex_code <- paste0(latex_code, "\\end{table}\n")

  return(latex_code)
}

# Generate LaTeX table
population_latex_table <- generate_population_effects_latex(
  data = population_effects_summary,
  caption = paste0("Population-Level Treatment Effects for Race = ", race_val,
                   " and Stage = ", stage_val, " (Months with 95\\% Credible Intervals)"),
  label = "tab:population_effects"
)

# Alternative simpler LaTeX table without footnotes
generate_population_effects_latex_simple <- function(data, caption, label) {

  latex_code <- "\\begin{table}[htbp]\n"
  latex_code <- paste0(latex_code, "\\centering\n")
  latex_code <- paste0(latex_code, "\\caption{", caption, "}\n")
  latex_code <- paste0(latex_code, "\\label{", label, "}\n")
  latex_code <- paste0(latex_code, "\\begin{tabular}{lcc}\n")
  latex_code <- paste0(latex_code, "\\hline\n")
  latex_code <- paste0(latex_code, "Treatment Effect & Sample Size & Effect (95\\% CI) \\\\\n")
  latex_code <- paste0(latex_code, " &  & (Months) \\\\\n")
  latex_code <- paste0(latex_code, "\\hline\n")

  for(i in 1:nrow(data)) {
    effect_display <- if(is.na(data$Mean_Months[i])) "NA" else data$Effect_Formatted[i]

    latex_code <- paste0(latex_code,
                         data$Effect[i], " & ",
                         data$Sample_Size[i], " & ",
                         effect_display, " \\\\\n"
    )
  }

  latex_code <- paste0(latex_code, "\\hline\n")
  latex_code <- paste0(latex_code, "\\end{tabular}\n")
  latex_code <- paste0(latex_code, "\\end{table}\n")

  return(latex_code)
}

population_latex_table_simple <- generate_population_effects_latex_simple(
  data = population_effects_summary,
  caption = paste0("Population-Level Treatment Effects for Race = ", race_val,
                   " and Stage = ", stage_val, " (Months with 95\\% Credible Intervals)"),
  label = "tab:population_effects_simple"
)

# ------------------------------------------------------------------------------
# 8. Save results
# ------------------------------------------------------------------------------

cat("\nSaving results...\n")

# Save summary table
write.csv(population_effects_summary, "population_effects_summary.csv", row.names = FALSE)

# Save posterior samples
population_samples <- data.frame(
  ATE = population_ATE_samples / days_per_month,
  ATT = if(n_treated > 0) population_ATT_samples / days_per_month else rep(NA, n_samples),
  ATU = if(n_untreated > 0) population_ATU_samples / days_per_month else rep(NA, n_samples),
  RATE = population_RATE_samples / days_per_month,
  RATT = if(n_treated > 0) population_RATT_samples / days_per_month else rep(NA, n_samples),
  RATU = if(n_untreated > 0) population_RATU_samples / days_per_month else rep(NA, n_samples)
)

write.csv(population_samples, "population_effects_posterior_samples.csv", row.names = FALSE)

# Save LaTeX tables
writeLines(population_latex_table, "population_effects_table.tex")
writeLines(population_latex_table_simple, "population_effects_table_simple.tex")

# Save all results
save(population_effects_summary, population_samples,
     population_ATE_samples, population_ATT_samples, population_ATU_samples,
     population_RATE_samples, population_RATT_samples, population_RATU_samples,
     population_latex_table, population_latex_table_simple,
     file = "population_effects_complete.RData")


cat("\n", paste(rep("=", 80), collapse=""), "\n")
cat("LATEX TABLE CODE\n")
cat(paste(rep("=", 80), collapse=""), "\n")

cat("\n1. DETAILED TABLE (with footnotes):\n")
cat(paste(rep("-", 50), collapse=""), "\n")
cat(population_latex_table)

cat("\n2. SIMPLE TABLE (no footnotes):\n")
cat(paste(rep("-", 50), collapse=""), "\n")
cat(population_latex_table_simple)

# ------------------------------------------------------------------------------
# 10. Summary statistics
# ------------------------------------------------------------------------------

cat("\n", paste(rep("=", 60), collapse=""), "\n")
cat("FILES SAVED\n")
cat(paste(rep("=", 60), collapse=""), "\n")

cat("CSV Files:\n")
cat("  - population_effects_summary.csv (main results table)\n")
cat("  - population_effects_posterior_samples.csv (", n_samples, " posterior samples)\n\n")

cat("LaTeX Files:\n")
cat("  - population_effects_table.tex (with footnotes)\n")
cat("  - population_effects_table_simple.tex (clean version)\n\n")

cat("R Data File:\n")
cat("  - population_effects_complete.RData (all objects)\n\n")

cat("Usage in LaTeX document:\n")
cat("\\input{population_effects_table_simple.tex}\n")

cat("\nPopulation-level treatment effects analysis complete!\n")
cat("Results show treatment effects averaged over the entire study population.\n")



####------------------







ATE_out <- estimate_ATE_corrected(
  X            = X,
  V            = V,
  group        = group,
  forest_ps    = forest_ps,
  forest_out   = forest_out,
  w_samples    = w_samples,
  sigma2_samps = sigma2_samps,
  burn_in      = burn_in,
  thin         = thin,
  n_mcmc       = n_mcmc,
  n_iter_ps    = n_iter_ps
)
cat("Marginal ATE (mean, lower, upper): ",
    round(ATE_out$ATE_mean,3), ", ",
    round(ATE_out$ATE_lower,3), ", ",
    round(ATE_out$ATE_upper,3), "\n")

# 5b.  Marginal RATE at t* = 5 yr (5×365 days)
RATE5_out <- estimate_RATE(
  t_star       = 5 * 365,
  X            = X,
  V            = V[,1, drop=FALSE],   # For estimate_RATE, V argument expects a vector or one‐column
  group        = group,
  forest_ps    = forest_ps,
  forest_out   = forest_out,
  w_samples    = w_samples,
  sigma2_samps = sigma2_samps,
  burn_in      = burn_in,
  thin         = thin,
  n_mcmc       = n_mcmc,
  n_iter_ps    = n_iter_ps
)
cat("Marginal RATE@5 yr (mean, lower, upper): ",
    round(RATE5_out$RATE_mean,3), ", ",
    round(RATE5_out$RATE_lower,3), ", ",
    round(RATE5_out$RATE_upper,3), "\n")

# 5c.  Marginal RATE at t* = 10 yr (10×365 days)
RATE10_out <- estimate_RATE(
  t_star       = 10 * 365,
  X            = X,
  V            = V[,1, drop=FALSE],
  group        = group,
  forest_ps    = forest_ps,
  forest_out   = forest_out,
  w_samples    = w_samples,
  sigma2_samps = sigma2_samps,
  burn_in      = burn_in,
  thin         = thin,
  n_mcmc       = n_mcmc,
  n_iter_ps    = n_iter_ps
)
cat("Marginal RATE@10 yr (mean, lower, upper): ",
    round(RATE10_out$RATE_mean,3), ", ",
    round(RATE10_out$RATE_lower,3), ", ",
    round(RATE10_out$RATE_upper,3), "\n")

# (5d)  Combine into a small data.frame and save:
marginal_effects <- data.frame(
  Measure   = c("ATE", "RATE_5yr", "RATE_10yr"),
  Estimate  = c(ATE_out$ATE_mean,
                RATE5_out$RATE_mean,
                RATE10_out$RATE_mean),
  Lower_CI  = c(ATE_out$ATE_lower,
                RATE5_out$RATE_lower,
                RATE10_out$RATE_lower),
  Upper_CI  = c(ATE_out$ATE_upper,
                RATE5_out$RATE_upper,
                RATE10_out$RATE_upper)
)
write.csv(marginal_effects, file = "marginal_effects_stratum.csv", row.names = FALSE)

# ------------------------------------------------------------------------------
# 6.  County‐specific ATE and RATE (for each county in this stratum)
# ------------------------------------------------------------------------------
county_results <- data.frame(
  county         = stratum_counties,
  sample_size    = sapply(stratum_counties, function(cc) sum(data_stratum$county==cc)),
  ATE_est        = NA,
  ATE_lower      = NA,
  ATE_upper      = NA,
  RATE5_est      = NA,
  RATE5_lower    = NA,
  RATE5_upper    = NA,
  RATE10_est     = NA,
  RATE10_lower   = NA,
  RATE10_upper   = NA
)

for(i in seq_along(stratum_counties)) {
  cid     <- stratum_counties[i]
  idx_i   <- which(data_stratum$county == cid)     # which rows belong to county i
  X_i     <- X[idx_i, , drop = FALSE]
  V_i_row <- V[idx_i, , drop = FALSE]              # same “V” repeated for each patient
  g_i     <- rep(i, length(idx_i))                 # group index = i for all patients in that county

  # 6a. County‐specific ATE
  c_ATE_out <- estimate_County_ATE(
    X_i          = X_i,
    v_i          = as.numeric(V[idx_i,1]),  # pulling the one‐column intercept vector
    county_id    = i,
    forest_ps    = forest_ps,
    forest_out   = forest_out,
    w_samples    = w_samples,
    sigma2_samps = sigma2_samps,
    burn_in      = burn_in,
    thin         = thin,
    n_mcmc       = n_mcmc,
    n_iter_ps    = n_iter_ps
  )
  county_results$ATE_est[i]   <- c_ATE_out$County_ATE_mean
  county_results$ATE_lower[i] <- c_ATE_out$County_ATE_lower
  county_results$ATE_upper[i] <- c_ATE_out$County_ATE_upper

  # 6b. County‐specific RATE at 5 yr
  c_RATE5_out <- estimate_County_RATE(
    t_star       = 5 * 365,
    X_i          = X_i,
    v_i          = as.numeric(V[idx_i,1]),
    county_id    = i,
    forest_ps    = forest_ps,
    forest_out   = forest_out,
    w_samples    = w_samples,
    sigma2_samps = sigma2_samps,
    burn_in      = burn_in,
    thin         = thin,
    n_mcmc       = n_mcmc,
    n_iter_ps    = n_iter_ps
  )
  county_results$RATE5_est[i]   <- c_RATE5_out$County_RATE_mean
  county_results$RATE5_lower[i] <- c_RATE5_out$County_RATE_lower
  county_results$RATE5_upper[i] <- c_RATE5_out$County_RATE_upper

print(i)
}

# Save county‐specific tables
write.csv(county_results, file = "county_specific_ATE_RATE_WD.csv", row.names = FALSE)


sorted_counties <- county_results %>% arrange(desc(sample_size))
big_county   <- sorted_counties$county[1]
small_county <- sorted_counties$county[nrow(sorted_counties)]
cat("Biggest county ID:", big_county, " (n=", sorted_counties$sample_size[1], ")\n")
cat("Smallest county ID:", small_county, " (n=", sorted_counties$sample_size[nrow(sorted_counties)], ")\n")


age_vals      <- c(50, 60, 74)
bxdelay_vals  <- c(0, 1)
hrp_vals      <- c(0, 1)
tgrade_vals   <- c(1, 2, 3)
covariate_grid <- expand.grid(
  Age      = age_vals,
  BX_Delay = bxdelay_vals,
  HR_p     = hrp_vals,
  Tgrade   = tgrade_vals
)
# Scale to [0,1] using the same scale01
covariate_grid_scaled <- data.frame(
  Age      = scale01(covariate_grid$Age),
  BX_Delay = scale01(covariate_grid$BX_Delay),
  HR_p     = scale01(covariate_grid$HR_p),
  Tgrade   = scale01(covariate_grid$Tgrade)
)
X_grid <- as.matrix(covariate_grid_scaled)


compute_CATE_CRATE_for_county <- function(county_id, is_big = TRUE) {
  idx_county <- which(data_stratum$county == county_id)
  X_i_all    <- X[idx_county, , drop = FALSE]

  V_i_row    <- rep(1, length(idx_county))

  county_index_in_group <- which(stratum_counties == county_id)


  results_df <- data.frame(
    county        = county_id,
    is_big        = ifelse(is_big, "Big", "Small"),
    Age           = covariate_grid$Age,
    BX_Delay      = covariate_grid$BX_Delay,
    HR_p          = covariate_grid$HR_p,
    Tgrade        = covariate_grid$Tgrade,
    CATE_est      = NA,
    CATE_lower    = NA,
    CATE_upper    = NA,
    CRATE5_est    = NA,
    CRATE5_lower  = NA,
    CRATE5_upper  = NA,
    CRATE10_est   = NA,
    CRATE10_lower = NA,
    CRATE10_upper = NA
  )

  for(j in seq_len(nrow(X_grid))) {
    x_j     <- matrix(X_grid[j, ], nrow = 1)
    v_i_vec <- V_i_row[1]  # always 1

    # 7c(i). CATE for this X
    cate_out <- estimate_CATE_corrected(
      x            = as.numeric(x_j),
      v_i          = v_i_vec,
      county_id    = county_index_in_group,
      forest_ps    = forest_ps,
      forest_out   = forest_out,
      w_samples    = w_samples,
      sigma2_samps = sigma2_samps,
      burn_in      = burn_in,
      thin         = thin,
      n_mcmc       = n_mcmc,
      n_iter_ps    = n_iter_ps
    )
    results_df$CATE_est[j]   <- cate_out$CATE_mean
    results_df$CATE_lower[j] <- cate_out$CATE_lower
    results_df$CATE_upper[j] <- cate_out$CATE_upper

    # 7c(ii). CRATE @5 yr
    crate5_out <- estimate_CRATE(
      t_star       = 5*365,
      x            = as.numeric(x_j),
      v_i          = v_i_vec,
      county_id    = county_index_in_group,
      forest_ps    = forest_ps,
      forest_out   = forest_out,
      w_samples    = w_samples,
      sigma2_samps = sigma2_samps,
      burn_in      = burn_in,
      thin         = thin,
      n_mcmc       = n_mcmc,
      n_iter_ps    = n_iter_ps
    )
    results_df$CRATE5_est[j]   <- crate5_out$CRATE_mean
    results_df$CRATE5_lower[j] <- crate5_out$CRATE_lower
    results_df$CRATE5_upper[j] <- crate5_out$CRATE_upper

    # 7c(iii). CRATE @10 yr
    crate10_out <- estimate_CRATE(
      t_star       = 10*365,
      x            = as.numeric(x_j),
      v_i          = v_i_vec,
      county_id    = county_index_in_group,
      forest_ps    = forest_ps,
      forest_out   = forest_out,
      w_samples    = w_samples,
      sigma2_samps = sigma2_samps,
      burn_in      = burn_in,
      thin         = thin,
      n_mcmc       = n_mcmc,
      n_iter_ps    = n_iter_ps
    )
    results_df$CRATE10_est[j]   <- crate10_out$CRATE_mean
    results_df$CRATE10_lower[j] <- crate10_out$CRATE_lower
    results_df$CRATE10_upper[j] <- crate10_out$CRATE_upper
  }

  return(results_df)
}


big_df   <- compute_CATE_CRATE_for_county(big_county,   is_big = TRUE)

small_df <- compute_CATE_CRATE_for_county(small_county, is_big = FALSE)

# Save each table:
write.csv(big_df,   file = sprintf("CATE_CRATE_grid_county_%s_big.csv",   big_county),   row.names = FALSE)
write.csv(small_df, file = sprintf("CATE_CRATE_grid_county_%s_small.csv", small_county), row.names = FALSE)

# ------------------------------------------------------------------------------

spte_results <- data.frame(
  county       = stratum_counties,
  SPTE_est     = NA,
  SPTE_lower   = NA,
  SPTE_upper   = NA
)
for(i in seq_along(stratum_counties)) {
  cid <- stratum_counties[i]

  # Extract X_i, V_i, group_i for that county
  idx_i   <- which(data_stratum$county == cid)
  X_i     <- X[idx_i, , drop = FALSE]
  V_i_row <- V[idx_i, , drop = FALSE]
  g_i     <- rep(i, length(idx_i))

  cs_spte <- estimate_County_SPTE(
    t            = 5*365,
    X_i          = X_i,
    V_i_row      = as.numeric(V_i_row[,1]),
    group_i      = g_i,
    forest_ps    = forest_ps,
    forest_out   = forest_out,
    w_samples    = w_samples,
    sigma2_samps = sigma2_samps,
    burn_in      = burn_in,
    thin         = thin,
    n_mcmc       = n_mcmc,
    n_iter_ps    = n_iter_ps
  )
  spte_results$SPTE_est[i]   <- cs_spte$County_SPTE_mean
  spte_results$SPTE_lower[i] <- cs_spte$County_SPTE_lower
  spte_results$SPTE_upper[i] <- cs_spte$County_SPTE_upper
}

# Save the SPTE table
write.csv(spte_results, "county_SPTE_stratum.csv", row.names = FALSE)


spte_results <- spte_results %>%
  mutate(
    county_label = as.character(county)
  )

ggplot(spte_results, aes(x = county_label, y = SPTE_est)) +
  geom_point(size = 2) +
  geom_errorbar(aes(ymin = SPTE_lower, ymax = SPTE_upper), width = 0.2) +
  labs(
    x = "County (ID)",
    y = "SPTE @5 yr",
    title = "County‐specific SPTE (stratum: race=1, stage=Distant)"
  ) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)
  )
ggsave("All_county_SPTE_5yr.png", width = 12, height = 6, dpi = 300)

cat("All steps complete. Tables and plots have been saved.\n")






time_years <- seq(0, 12, by = 4/12)
time_days  <- round(time_years * 365)

# Check length
cat("Number of time‐points:", length(time_days), "\n")  # should be 37


spte_marginal <- data.frame(
  t_years    = time_years,
  t_days     = time_days,
  SPTE_mean  = NA_real_,
  SPTE_lower = NA_real_,
  SPTE_upper = NA_real_
)


for(i in seq_along(time_days)) {
  t_i <- time_days[i]

  out_i <- estimate_SPTE(
    t            = t_i,
    X            = X,                 # n×p matrix of patient‐level covariates
    V            = V[,1],             # if V is intercept‐only, V[,1] is a vector of 1’s
    group        = group,             # length‐n vector of county indices (1:K_stratum)
    forest_ps    = forest_ps,
    forest_out   = forest_out,
    w_samples    = w_samples,
    sigma2_samps = sigma2_samps,
    burn_in      = burn_in,
    thin         = thin,
    n_mcmc       = n_mcmc,
    n_iter_ps    = n_iter_ps
  )

  # (C.2) Store mean, lower, upper
  spte_marginal$SPTE_mean[i]  <- out_i$SPTE_mean
  spte_marginal$SPTE_lower[i] <- out_i$SPTE_lower
  spte_marginal$SPTE_upper[i] <- out_i$SPTE_upper
  print(i)
}

# Write marginal SPTE to CSV
write.csv(spte_marginal,
          file = "SPTE_marginal_time_grid.csv",
          row.names = FALSE)


spte_county <- expand.grid(
  county   = stratum_counties,  # e.g. 1:K_stratum
  t_years  = time_years
) %>%
  arrange(county, t_years) %>%
  as.data.frame()

# Add t_days and empty columns
spte_county <- spte_county %>%
  mutate(
    t_days     = round(t_years * 365),
    SPTE_mean  = NA_real_,
    SPTE_lower = NA_real_,
    SPTE_upper = NA_real_
  )


for(j in seq_len(nrow(spte_county))) {
  cid     <- spte_county$county[j]
  t_i     <- spte_county$t_days[j]

  county_idx_in_group <- which(stratum_counties == cid)  # 1…K_stratum

  idx_i <- which(data_stratum$county == cid)
  X_i   <- X[idx_i, , drop = FALSE]     # (n_i × p)


  V_i   <- rep(V[idx_i, 1][1], length(idx_i))  # length‐n_i vector


  g_i   <- rep(county_idx_in_group, length(idx_i))


  cs_out <- estimate_County_SPTE(
    t            = t_i,
    X_i          = X_i,
    V_i_row      = V_i,
    group_i      = g_i,
    forest_ps    = forest_ps,
    forest_out   = forest_out,
    w_samples    = w_samples,
    sigma2_samps = sigma2_samps,
    burn_in      = burn_in,
    thin         = thin,
    n_mcmc       = n_mcmc,
    n_iter_ps    = n_iter_ps
  )

  # (E.2) Store
  spte_county$SPTE_mean[j]  <- cs_out$County_SPTE_mean
  spte_county$SPTE_lower[j] <- cs_out$County_SPTE_lower
  spte_county$SPTE_upper[j] <- cs_out$County_SPTE_upper
  print(j)
}


write.csv(spte_county,
          file = "SPTE_county_time_grid.csv",
          row.names = FALSE)


library(ggplot2)
ggplot(spte_marginal, aes(x = t_years, y = SPTE_mean)) +
  geom_line() +
  geom_ribbon(aes(ymin = SPTE_lower, ymax = SPTE_upper), alpha = 0.2) +
  labs(
    x = "Time (years)",
    y = "Marginal SPTE",
    title = "Marginal SPTE from 0 to 12 years (4‐month grid)"
  ) +
  theme_minimal()

ggsave("SPTE_marginal_time_grid.png", width = 6, height = 4, dpi = 300)


one_county  <- stratum_counties[1]
temp_county <- spte_county %>% filter(county == one_county)

ggplot(temp_county, aes(x = t_years, y = SPTE_mean)) +
  geom_line(color = "steelblue") +
  geom_ribbon(aes(ymin = SPTE_lower, ymax = SPTE_upper), alpha = 0.2, fill = "steelblue") +
  labs(
    x = "Time (years)",
    y = paste0("County ", one_county, " SPTE"),
    title = paste0("County‐specific SPTE (County ", one_county, ")")
  ) +
  theme_minimal()

ggsave(sprintf("County_%s_SPTE_time_grid.png", one_county),
       width = 6, height = 4, dpi = 300)














library(ggplot2)

# Assuming in your R session you have:
#   spte_marginal   # data.frame with t_years, SPTE_mean, SPTE_lower, SPTE_upper
#   ATE_out         # list with ATE_mean, ATE_lower, ATE_upper
#   RATE5_out       # list with RATE_mean, RATE_lower, RATE_upper

# Compose the annotation strings:
ate_label  <- sprintf(
  " ATE: %.2f (%.2f, %.2f)",
  ATE_out$ATE_mean,
  ATE_out$ATE_lower,
  ATE_out$ATE_upper
)
rate_label <- sprintf(
  "RATE (5yr): %.2f (%.2f, %.2f)",
  RATE5_out$RATE_mean,
  RATE5_out$RATE_lower,
  RATE5_out$RATE_upper
)

# Find reasonable y‐positions for the text (just above the top of the SPTE ribbon):
y_max <- max(spte_marginal$SPTE_upper, na.rm = TRUE)
y_min <- min(spte_marginal$SPTE_lower, na.rm = TRUE)
y_range <- y_max - y_min

# Place ATE label near x = 1 year, y = y_max + 0.05*y_range
# Place RATE@5yr label near x = 1 year, y = y_max + 0.10*y_range (slightly above ATE text)
x_pos <- 1
y_ate <- y_max + 0.05 * y_range
y_rate <- y_max + 0.10 * y_range

ggplot(spte_marginal, aes(x = t_years, y = SPTE_mean)) +
  geom_line(size = 1) +
  geom_ribbon(aes(ymin = SPTE_lower, ymax = SPTE_upper), alpha = 0.2) +
  # Add the two annotation texts:
  annotate(
    "text",
    x     = x_pos,
    y     = y_ate,
    label = ate_label,
    hjust = 0,
    size  = 4
  ) +
  annotate(
    "text",
    x     = x_pos,
    y     = y_rate,
    label = rate_label,
    hjust = 0,
    size  = 4
  ) +
  labs(
    x     = "Time (years)",
    y     = "Marginal SPTE",
    title = "Marginal SPTE from 0 to 12 years (4‐month grid)"
  ) +
  theme_minimal() +
  theme(
    plot.background   = element_rect(fill = "white", color = NA),
    panel.grid.major  = element_blank(),
    panel.grid.minor  = element_blank(),
    panel.background  = element_rect(fill = "white", color = NA)
  )




library(ggplot2)
library(dplyr)

# (1) Example data frame: spte_marginal has t_days, SPTE_mean, SPTE_lower, SPTE_upper
#     We'll add t_years = t_days / 365
spte_marginal <- spte_marginal %>%
  mutate(
    t_years = t_days / 365
  )

# (2) Build the annotation strings, converting days to years:
ate_label  <- sprintf(
  " ATE: %.2f (%.2f, %.2f)",
  ATE_out$ATE_mean/365,
  ATE_out$ATE_lower/365,
  ATE_out$ATE_upper/365
)

# If RATE5_out was computed at 1825 days:
rate_label <- sprintf(
  "RATE in %d  years: %.2f (%.2f, %.2f)",
  (5 * 365) / 365,            # 5
  RATE5_out$RATE_mean/365,
  RATE5_out$RATE_lower/365,
  RATE5_out$RATE_upper/365
)

# (3) Find a good vertical position just above SPTE’s maximum:
y_max   <- max(spte_marginal$SPTE_upper, na.rm = TRUE)
y_range <- diff(range(spte_marginal$SPTE_lower, spte_marginal$SPTE_upper, na.rm = TRUE))

# Position the two texts slightly above that band:
x_pos   <- 1             # place at x = 1 year
y_ate   <- y_max + 0.05 * y_range
y_rate  <- y_max + 0.10 * y_range

# (4) Plot:
ggplot(spte_marginal, aes(x = t_years, y = SPTE_mean)) +
  geom_line(size = 1) +
  geom_ribbon(aes(ymin = SPTE_lower, ymax = SPTE_upper), alpha = 0.2) +
  # Annotate ATE text:
  annotate(
    "text",
    x     = x_pos,
    y     = y_ate,
    label = ate_label,
    hjust = 0,
    size  = 4
  ) +
  # Annotate RATE@5yr text:
  annotate(
    "text",
    x     = x_pos,
    y     = y_rate,
    label = rate_label,
    hjust = 0,
    size  = 4
  ) +
  labs(
    x     = "Time (years)",
    y     = "SPTE",
    title = ""
  ) +
  theme_minimal() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_rect(fill = "white", color = NA),
    plot.background  = element_rect(fill = "white", color = NA)
  )






# A vector of all 67 Florida‐county FIPS (in the same order your county IDs run 1:67)
florida_fips <- c(
  "12001","12003","12005","12007","12009","12011","12013","12015","12017","12019",
  "12021","12023","12027","12029","12031","12033","12035","12037","12039","12041",
  "12043","12045","12047","12049","12051","12053","12055","12057","12059","12061",
  "12063","12065","12067","12069","12071","12073","12075","12077","12079","12081",
  "12083","12085","12086","12087","12089","12091","12093","12095","12097","12099",
  "12101","12103","12105","12107","12109","12111","12113","12115","12117","12119",
  "12121","12123","12125","12127","12129"
)

county_results <- county_results %>%
  mutate(
    # If 'county' is numeric 1:67, map to the corresponding FIPS string
    fips = florida_fips[county]
  )

# Now join to fl_counties by GEOID:
county_map <- fl_counties %>%
  left_join(
    county_results,
    by = c("GEOID" = "fips")
  )






library(ggplot2)

county_map <- county_map %>%
  mutate(
    ATE_years = ATE_est / 365,
    RATE5     = RATE5_est/365
  )

# Plot 1: ATE (in years)
ggplot(county_map) +
  geom_sf(aes(fill = ATE_years), color = "white", size = 0.2) +
  scale_fill_viridis_c(
    name    = "ATE (years)",
    option  = "plasma",
    na.value= "grey90"
  ) +
  labs(
    title    = " County ATE ",
    subtitle = ""
  ) +
  theme_void() +
  theme(
    legend.position   = "right",
    plot.title        = element_text(size = 14, face = "bold", hjust = 0.5),
    plot.subtitle     = element_text(size = 11, hjust = 0.5),
    legend.title      = element_text(size = 10),
    legend.text       = element_text(size = 9),
    panel.background  = element_rect(fill = "white", color = NA),
    plot.background   = element_rect(fill = "white", color = NA)
  )

# Plot 2: RATE @ 5 years
ggplot(county_map) +
  geom_sf(aes(fill = RATE5), color = "white", size = 0.2) +
  scale_fill_viridis_c(
    name    = "RATE @ 5 yr",
    option  = "magma",
    na.value= "grey90"
  ) +
  labs(
    title    = "County RATE at 5 Years",
    subtitle = ""
  ) +
  theme_void() +
  theme(
    legend.position   = "right",
    plot.title        = element_text(size = 14, face = "bold", hjust = 0.5),
    plot.subtitle     = element_text(size = 11, hjust = 0.5),
    legend.title      = element_text(size = 10),
    legend.text       = element_text(size = 9),
    panel.background  = element_rect(fill = "white", color = NA),
    plot.background   = element_rect(fill = "white", color = NA)
  )






# -------------------------------------------------------
# 0. Prepare posterior‐iteration indices (once, globally)
# -------------------------------------------------------
sample_iters_out <- seq(from = burn_in + 1, to = n_mcmc, by = thin)
Mprime           <- length(sample_iters_out)

# -------------------------------------------------------
# 1. ATT over entire stratum
# -------------------------------------------------------
ATT_out <- estimate_ATT_corrected(
  X            = X,
  V            = V[,1],      # if V is intercept‐only
  group        = group,
  Z            = Z,
  forest_ps    = forest_ps,
  forest_out   = forest_out,
  w_samples    = w_samples,
  sigma2_samps = sigma2_samps,
  burn_in      = burn_in,
  thin         = thin,
  n_mcmc       = n_mcmc,
  n_iter_ps    = n_iter_ps
)
print(ATT_out$ATT_mean)
print(ATT_out$ATT_lower)
print(ATT_out$ATT_upper)

# -------------------------------------------------------
# 2. ATU over entire stratum
# -------------------------------------------------------
ATU_out <- estimate_ATU_corrected(
  X            = X,
  V            = V[,1],
  group        = group,
  Z            = Z,
  forest_ps    = forest_ps,
  forest_out   = forest_out,
  w_samples    = w_samples,
  sigma2_samps = sigma2_samps,
  burn_in      = burn_in,
  thin         = thin,
  n_mcmc       = n_mcmc,
  n_iter_ps    = n_iter_ps
)
print(ATU_out$ATU_mean)
print(ATU_out$ATU_lower)
print(ATU_out$ATU_upper)

# -------------------------------------------------------
# 3.  County‐specific ATT and ATU for each county in stratum
# -------------------------------------------------------
county_ids <- sort(unique(data_stratum$county))  # e.g., 1:K_stratum
n_counties <- length(county_ids)

county_att_results <- data.frame(
  county         = county_ids,
  sample_size    = NA_integer_,
  County_ATT_est = NA_real_,
  County_ATT_lower = NA_real_,
  County_ATT_upper = NA_real_
)

county_atu_results <- data.frame(
  county         = county_ids,
  sample_size    = NA_integer_,
  County_ATU_est = NA_real_,
  County_ATU_lower = NA_real_,
  County_ATU_upper = NA_real_
)

for (i in seq_along(county_ids)) {
  cid <- county_ids[i]
  idx_i   <- which(data_stratum$county == cid)
  X_i     <- X[idx_i, , drop = FALSE]
  V_i_row <- V[idx_i, 1]    # intercept‐only column, repeated
  Z_i     <- Z[idx_i]
  group_i <- rep(i, length(idx_i))  # county index in {1..K_stratum}

  county_att_results$sample_size[i] <- length(idx_i)
  county_atu_results$sample_size[i] <- length(idx_i)

  # (3a) CountyATT
  cat("Computing CountyATT for county:", cid, "\n")
  cat_att <- estimate_County_ATT(
    X_i          = X_i,
    v_i          = V_i_row,
    Z_i          = Z_i,
    county_id    = i,            # i matches “group indexes” used in fit
    forest_ps    = forest_ps,
    forest_out   = forest_out,
    w_samples    = w_samples,
    sigma2_samps = sigma2_samps,
    burn_in      = burn_in,
    thin         = thin,
    n_mcmc       = n_mcmc,
    n_iter_ps    = n_iter_ps
  )
  county_att_results$County_ATT_est[i]   <- cat_att$County_ATT_mean
  county_att_results$County_ATT_lower[i] <- cat_att$County_ATT_lower
  county_att_results$County_ATT_upper[i] <- cat_att$County_ATT_upper

  # (3b) CountyATU
  cat("Computing CountyATU for county:", cid, "\n")
  cat_atu <- estimate_County_ATU(
    X_i          = X_i,
    v_i          = V_i_row,
    Z_i          = Z_i,
    county_id    = i,
    forest_ps    = forest_ps,
    forest_out   = forest_out,
    w_samples    = w_samples,
    sigma2_samps = sigma2_samps,
    burn_in      = burn_in,
    thin         = thin,
    n_mcmc       = n_mcmc,
    n_iter_ps    = n_iter_ps
  )
  county_atu_results$County_ATU_est[i]   <- cat_atu$County_ATU_mean
  county_atu_results$County_ATU_lower[i] <- cat_atu$County_ATU_lower
  county_atu_results$County_ATU_upper[i] <- cat_atu$County_ATU_upper
}

# Save to CSV if desired:
write.csv(county_att_results, "county_ATT_results.csv", row.names = FALSE)
write.csv(county_atu_results, "county_ATU_results.csv", row.names = FALSE)

















# -------------------------------------------------------
# 0.  (Same as before) Prepare posterior‐iteration indices
# -------------------------------------------------------
sample_iters_out <- seq(from = burn_in + 1, to = n_mcmc, by = thin)
Mprime           <- length(sample_iters_out)

# -------------------------------------------------------
# 1.  (Same as before) Global ATT and ATU
# -------------------------------------------------------
ATT_out <- estimate_ATT_corrected(
  X            = X,
  V            = V[,1],
  group        = group,
  Z            = Z,
  forest_ps    = forest_ps,
  forest_out   = forest_out,
  w_samples    = w_samples,
  sigma2_samps = sigma2_samps,
  burn_in      = burn_in,
  thin         = thin,
  n_mcmc       = n_mcmc,
  n_iter_ps    = n_iter_ps
)

ATU_out <- estimate_ATU_corrected(
  X            = X,
  V            = V[,1],
  group        = group,
  Z            = Z,
  forest_ps    = forest_ps,
  forest_out   = forest_out,
  w_samples    = w_samples,
  sigma2_samps = sigma2_samps,
  burn_in      = burn_in,
  thin         = thin,
  n_mcmc       = n_mcmc,
  n_iter_ps    = n_iter_ps
)

# -------------------------------------------------------
# 2.  County‐specific ATT/ATU with NA if no (un)treated
# -------------------------------------------------------
county_ids <- sort(unique(data_stratum$county))
n_counties <- length(county_ids)

county_att_results <- data.frame(
  county            = county_ids,
  sample_size       = NA_integer_,
  n_treated         = NA_integer_,
  n_untreated       = NA_integer_,
  County_ATT_est    = NA_real_,
  County_ATT_lower  = NA_real_,
  County_ATT_upper  = NA_real_,
  County_ATU_est    = NA_real_,
  County_ATU_lower  = NA_real_,
  County_ATU_upper  = NA_real_
)

for (i in seq_along(county_ids)) {
  cid     <- county_ids[i]
  idx_i   <- which(data_stratum$county == cid)
  X_i     <- X[idx_i, , drop = FALSE]
  V_i_row <- V[idx_i, 1]    # intercept‐only
  Z_i     <- Z[idx_i]

  n_i         <- length(idx_i)
  n1_i        <- sum(Z_i == 1)
  n0_i        <- sum(Z_i == 0)

  county_att_results$sample_size[i]  <- n_i
  county_att_results$n_treated[i]    <- n1_i
  county_att_results$n_untreated[i]  <- n0_i

  # (2a) CountyATT only if n1_i > 0
  if (n1_i > 0) {
    cat("Computing CountyATT for county:", cid, " (n_treated =", n1_i, ")\n")
    cat_att <- estimate_County_ATT(
      X_i          = X_i,
      v_i          = V_i_row,
      Z_i          = Z_i,
      county_id    = i,
      forest_ps    = forest_ps,
      forest_out   = forest_out,
      w_samples    = w_samples,
      sigma2_samps = sigma2_samps,
      burn_in      = burn_in,
      thin         = thin,
      n_mcmc       = n_mcmc,
      n_iter_ps    = n_iter_ps
    )
    county_att_results$County_ATT_est[i]   <- cat_att$County_ATT_mean
    county_att_results$County_ATT_lower[i] <- cat_att$County_ATT_lower
    county_att_results$County_ATT_upper[i] <- cat_att$County_ATT_upper
  } else {
    cat("Skipping CountyATT for county:", cid, "- no treated units.\n")
    county_att_results$County_ATT_est[i]   <- NA
    county_att_results$County_ATT_lower[i] <- NA
    county_att_results$County_ATT_upper[i] <- NA
  }

  # (2b) CountyATU only if n0_i > 0
  if (n0_i > 0) {
    cat("Computing CountyATU for county:", cid, " (n_untreated =", n0_i, ")\n")
    cat_atu <- estimate_County_ATU(
      X_i          = X_i,
      v_i          = V_i_row,
      Z_i          = Z_i,
      county_id    = i,
      forest_ps    = forest_ps,
      forest_out   = forest_out,
      w_samples    = w_samples,
      sigma2_samps = sigma2_samps,
      burn_in      = burn_in,
      thin         = thin,
      n_mcmc       = n_mcmc,
      n_iter_ps    = n_iter_ps
    )
    county_att_results$County_ATU_est[i]   <- cat_atu$County_ATU_mean
    county_att_results$County_ATU_lower[i] <- cat_atu$County_ATU_lower
    county_att_results$County_ATU_upper[i] <- cat_atu$County_ATU_upper
  } else {
    cat("Skipping CountyATU for county:", cid, "- no untreated units.\n")
    county_att_results$County_ATU_est[i]   <- NA
    county_att_results$County_ATU_lower[i] <- NA
    county_att_results$County_ATU_upper[i] <- NA
  }
  print(i)
}

# -------------------------------------------------------
# 3.  Save results
# -------------------------------------------------------
write.csv(county_att_results, "county_ATT_ATU_results.csv", row.names = FALSE)















library(dplyr)
library(sf)
library(tigris)
library(ggplot2)

# (A) If county_att_results does NOT yet have a “fips” column, create it now.
#     Here we assume county = 1 maps to "12001", county = 2 → "12003", …, county = 67 → "12129".
florida_fips <- c(
  "12001","12003","12005","12007","12009","12011","12013","12015","12017","12019",
  "12021","12023","12027","12029","12031","12033","12035","12037","12039","12041",
  "12043","12045","12047","12049","12051","12053","12055","12057","12059","12061",
  "12063","12065","12067","12069","12071","12073","12075","12077","12079","12081",
  "12083","12085","12086","12087","12089","12091","12093","12095","12097","12099",
  "12101","12103","12105","12107","12109","12111","12113","12115","12117","12119",
  "12121","12123","12125","12127","12129"
)

county_att_results <- county_att_results %>%
  mutate(
    fips = florida_fips[county]
  )

# (B) Load Florida county shapefile as an sf object
fl_counties <- counties(state = "FL", cb = TRUE, class = "sf")

# (C) Left‐join ATT/ATU results onto the shapefile by FIPS
county_map <- fl_counties %>%
  left_join(
    county_att_results %>% mutate(fips = as.character(fips)),
    by = c("GEOID" = "fips")
  )

# (D) Convert ATT and ATU from days → years (if not already)
county_map <- county_map %>%
  mutate(
    ATT_years = County_ATT_est / 365,
    ATU_years = County_ATU_est / 365
  )

# ------------------------------------------------------------------------------
# 1.  Choropleth: County_ATT (in years)
# ------------------------------------------------------------------------------
ggplot(county_map) +
  geom_sf(
    aes(fill = ATT_years),
    color = "white",
    size  = 0.2
  ) +
  scale_fill_viridis_c(
    name     = "ATT (years)",
    option   = "plasma",
    na.value = "grey90"
  ) +
  labs(
    title    = "",
    subtitle = "",
    caption  = ""
  ) +
  theme_void() +
  theme(
    legend.position   = "right",
    plot.title        = element_text(size = 14, face = "bold", hjust = 0.5),
    plot.subtitle     = element_text(size = 11, hjust = 0.5),
    legend.title      = element_text(size = 10),
    legend.text       = element_text(size = 9),
    panel.background  = element_rect(fill = "white", color = NA),
    plot.background   = element_rect(fill = "white", color = NA)
  )

# (Optional) Save to file:
# ggsave("Florida_County_ATT_years.png", width = 8, height = 6, dpi = 300)


# ------------------------------------------------------------------------------
# 2.  Choropleth: County_ATU (in years)
# ------------------------------------------------------------------------------
ggplot(county_map) +
  geom_sf(
    aes(fill = ATU_years),
    color = "white",
    size  = 0.2
  ) +
  scale_fill_viridis_c(
    name     = "ATU (years)",
    option   = "magma",
    na.value = "grey90"
  ) +
  labs(
    title    = "",
    subtitle = "",
    caption  = ""
  ) +
  theme_void() +
  theme(
    legend.position   = "right",
    plot.title        = element_text(size = 14, face = "bold", hjust = 0.5),
    plot.subtitle     = element_text(size = 11, hjust = 0.5),
    legend.title      = element_text(size = 10),
    legend.text       = element_text(size = 9),
    panel.background  = element_rect(fill = "white", color = NA),
    plot.background   = element_rect(fill = "white", color = NA)
  )

# (Optional) Save to file:
# ggsave("Florida_County_ATU_years.png", width = 8, height = 6, dpi = 300)


# (1) Extract the vector of posterior draws for rho
rho_draws <- fit$rho_samps

# (2) Compute posterior summary:
rho_mean <- mean(rho_draws)
rho_ci   <- quantile(rho_draws, probs = c(0.025, 0.975))

# (3) Print results
cat("Posterior mean of rho:    ", round(rho_mean, 3), "\n")
cat("95% credible interval:    ",
    paste0("(", round(rho_ci[1], 3), ", ", round(rho_ci[2], 3), ")\n"))
hist(rho_draws,main = "Posterior Distribution of ")
boxplot(rho_draws)



# Assuming you have already run:
#   fit <- AFT_mixed_DAGAR_probit_causal(…)
# and that the posterior draws of rho are stored in:
#   rho_draws <- fit$rho_samps

library(ggplot2)

# Extract rho posterior draws
rho_draws <- fit$rho_samps

# Create a data frame for ggplot
df_rho <- data.frame(rho = rho_draws)

# Assuming `df_rho` exists with a column named “rho”:

library(ggplot2)

ggplot(df_rho, aes(x = rho)) +
  geom_histogram(binwidth = 0.02, fill = "steelblue", color = "white") +
  labs(
    x     = expression(rho),
    y     = "Frequency"
  ) +
  theme_classic()



















# ------------------------------------------------------------------------------
# 1.  Load required packages
# ------------------------------------------------------------------------------
if (!requireNamespace("knitr", quietly = TRUE)) {
  install.packages("knitr")
}
library(knitr)
library(dplyr)
library(readr)
library(stringr)

# ------------------------------------------------------------------------------
# 2.  Read in the CSVs/data frames we created earlier
# ------------------------------------------------------------------------------
# 2.1. Marginal effects: marginal_effects_stratum.csv
marginal_effects <- read_csv("marginal_effects_stratum.csv",
                             col_types = cols(
                               Measure  = col_character(),
                               Estimate = col_double(),
                               Lower_CI = col_double(),
                               Upper_CI = col_double()
                             ))

# 2.2. County‐specific ATE & RATE@5yr: county_specific_ATE_RATE.csv
county_results   <- read_csv("county_specific_ATE_RATE.csv",
                             col_types = cols(
                               county       = col_character(),
                               sample_size  = col_integer(),
                               ATE_est      = col_double(),
                               ATE_lower    = col_double(),
                               ATE_upper    = col_double(),
                               RATE5_est    = col_double(),
                               RATE5_lower  = col_double(),
                               RATE5_upper  = col_double(),
                               RATE10_est   = col_double(),   # we’ll drop this
                               RATE10_lower = col_double(),
                               RATE10_upper = col_double()
                             ))

# 2.3. CATE/CRATE grids (big & small): big and small CSVs
#     (e.g. "CATE_CRATE_grid_county_<big>_big.csv" and "..._small.csv")
#     For simplicity, assume you saved them as:
big_df   <- read_csv("CATE_CRATE_grid_county_6_big.csv",
                     col_types = cols(
                       Age       = col_double(),
                       BX_Delay  = col_integer(),
                       HR_p      = col_double(),
                       Tgrade    = col_integer(),
                       CATE_est  = col_double(),
                       CATE_lower= col_double(),
                       CATE_upper= col_double(),
                       CRATE5_est   = col_double(),
                       CRATE5_lower = col_double(),
                       CRATE5_upper = col_double(),
                       CRATE10_est   = col_double(),  # drop
                       CRATE10_lower = col_double(),
                       CRATE10_upper = col_double()
                     ))

small_df <- read_csv("CATE_CRATE_grid_county_30_small.csv",
                     col_types = cols(
                       Age       = col_double(),
                       BX_Delay  = col_integer(),
                       HR_p      = col_double(),
                       Tgrade    = col_integer(),
                       CATE_est  = col_double(),
                       CATE_lower= col_double(),
                       CATE_upper= col_double(),
                       CRATE5_est   = col_double(),
                       CRATE5_lower = col_double(),
                       CRATE5_upper = col_double(),
                       CRATE10_est   = col_double(),  # drop
                       CRATE10_lower = col_double(),
                       CRATE10_upper = col_double()
                     ))


spte_county <- read_csv("SPTE_county_time_grid.csv",
                        col_types = cols(
                          county       = col_character(),
                          t_years      = col_double(),
                          t_days       = col_double(),
                          SPTE_mean    = col_double(),
                          SPTE_lower   = col_double(),
                          SPTE_upper   = col_double()
                        ))

# Filter only t_years == 5 (i.e. t_days == 5*365)
spte_5yr <- spte_county %>% filter(abs(t_years - 5) < 1e-6)

# 2.5. County‐specific ATT/ATU: county_ATT_ATU_results.csv
county_att_atu <- read_csv("county_ATT_ATU_results.csv",
                           col_types = cols(
                             county            = col_character(),
                             sample_size       = col_integer(),
                             n_treated         = col_integer(),
                             n_untreated       = col_integer(),
                             County_ATT_est    = col_double(),
                             County_ATT_lower  = col_double(),
                             County_ATT_upper  = col_double(),
                             County_ATU_est    = col_double(),
                             County_ATU_lower  = col_double(),
                             County_ATU_upper  = col_double()
                           ))

# ------------------------------------------------------------------------------
# 3.  Transform each data frame so that Estimate±CI appear in one column
# ------------------------------------------------------------------------------
# 3.1. Marginal effects
marginal_latex <- marginal_effects %>%
  mutate(
    `Estimate (95\\% CI)` =
      sprintf("%.2f (%.2f--%.2f)", Estimate, Lower_CI, Upper_CI)
  ) %>%
  select(Measure, `Estimate (95\\% CI)`)

# 3.2. County‐specific ATE & RATE@5yr
county_latex <- county_results %>%
  transmute(
    County       = county,
    `Sample Size`= sample_size,
    `ATE (95\\% CI)`    = sprintf(
      "%.2f (%.2f--%.2f)",
      ATE_est, ATE_lower, ATE_upper
    ),
    `RATE@5yr (95\\% CI)` = sprintf(
      "%.2f (%.2f--%.2f)",
      RATE5_est, RATE5_lower, RATE5_upper
    )
  )

# 3.3. CATE & CRATE@5yr for Big county
big_latex <- big_df %>%
  transmute(
    Age       = Age,
    BX_Delay  = BX_Delay,
    HR_p      = HR_p,
    Tgrade    = Tgrade,
    `CATE (95\\% CI)`   = sprintf(
      "%.2f (%.2f--%.2f)",
      CATE_est, CATE_lower, CATE_upper
    ),
    `CRATE@5yr (95\\% CI)` = sprintf(
      "%.2f (%.2f--%.2f)",
      CRATE5_est, CRATE5_lower, CRATE5_upper
    )
  )

# 3.4. CATE & CRATE@5yr for Small county
small_latex <- small_df %>%
  transmute(
    Age       = Age,
    BX_Delay  = BX_Delay,
    HR_p      = HR_p,
    Tgrade    = Tgrade,
    `CATE (95\\% CI)`   = sprintf(
      "%.2f (%.2f--%.2f)",
      CATE_est, CATE_lower, CATE_upper
    ),
    `CRATE@5yr (95\\% CI)` = sprintf(
      "%.2f (%.2f--%.2f)",
      CRATE5_est, CRATE5_lower, CRATE5_upper
    )
  )

# 3.5. County‐specific SPTE@5yr
spte5_latex <- spte_5yr %>%
  transmute(
    County        = county,
    `SPTE@5yr (95\\% CI)` = sprintf(
      "%.2f (%.2f--%.2f)",
      SPTE_mean, SPTE_lower, SPTE_upper
    )
  )

# 3.6. County‐specific ATT & ATU
county_attatu_latex <- county_att_atu %>%
  transmute(
    County       = county,
    `Sample Size` = sample_size,
    `# Treated`   = n_treated,
    `ATT (95\\% CI)` = if_else(
      is.na(County_ATT_est),
      "--",
      sprintf("%.2f (%.2f--%.2f)",
              County_ATT_est, County_ATT_lower, County_ATT_upper)
    ),
    `ATU (95\\% CI)` = if_else(
      is.na(County_ATU_est),
      "--",
      sprintf("%.2f (%.2f--%.2f)",
              County_ATU_est, County_ATU_lower, County_ATU_upper)
    )
  )

# ------------------------------------------------------------------------------
# 4.  Generate LaTeX code for each table using knitr::kable()
# ------------------------------------------------------------------------------
# 4.1. Marginal Effects
cat("%% ---------------------- Marginal Effects ----------------------\n")
cat(kable(
  marginal_latex,
  format     = "latex",
  booktabs   = TRUE,
  linesep    = "",
  caption    = "Marginal ATE and RATE@5yr (95\\% CI)",
  label      = "tab:marginal_effects"
),
"\n\n")

# 4.2. County‐Specific ATE & RATE@5yr
cat("%% --------------- County‐Specific ATE & RATE@5yr --------------\n")
cat(kable(
  county_latex,
  format     = "latex",
  booktabs   = TRUE,
  linesep    = "",
  caption    = "County‐Specific ATE and RATE@5yr (95\\% CI)",
  label      = "tab:county_ate_rate5"
),
"\n\n")

# 4.3. CATE & CRATE@5yr for Big county
cat("%% --------- CATE & CRATE@5yr for Big County ---------\n")
cat(kable(
  big_latex,
  format     = "latex",
  booktabs   = TRUE,
  linesep    = "",
  caption    = "CATE and CRATE@5yr for Big County (95\\% CI)",
  label      = "tab:CATE_CRATE_big"
),
"\n\n")

# 4.4. CATE & CRATE@5yr for Small county
cat("%% -------- CATE & CRATE@5yr for Small County --------\n")
cat(kable(
  small_latex,
  format     = "latex",
  booktabs   = TRUE,
  linesep    = "",
  caption    = "CATE and CRATE@5yr for Small County (95\\% CI)",
  label      = "tab:CATE_CRATE_small"
),
"\n\n")

# 4.5. County‐Specific SPTE@5yr
cat("%% -------------- County‐Specific SPTE@5yr ---------------\n")
cat(kable(
  spte5_latex,
  format     = "latex",
  booktabs   = TRUE,
  linesep    = "",
  caption    = "County‐Specific SPTE@5yr (95\\% CI)",
  label      = "tab:county_spte"
),
"\n\n")

# 4.6. County‐Specific ATT & ATU
cat("%% ---------- County‐Specific ATT & ATU (Days→Years) ----------\n")
cat(kable(
  county_attatu_latex,
  format     = "latex",
  booktabs   = TRUE,
  linesep    = "",
  caption    = "County‐Specific ATT and ATU (95\\% CI)",
  label      = "tab:county_att_atu"
),
"\n\n")

# ------------------------------------------------------------------------------
# 5.  (Optional) Write each LaTeX table to its own .tex file
# ------------------------------------------------------------------------------
writeLines(
  kable(
    marginal_latex,
    format   = "latex",
    booktabs = TRUE,
    linesep  = "",
    caption  = "Marginal ATE and RATE@5yr (95\\% CI)",
    label    = "tab:marginal_effects"
  ),
  con = "table_marginal_effects.tex"
)

writeLines(
  kable(
    county_latex,
    format   = "latex",
    booktabs = TRUE,
    linesep  = "",
    caption  = "County‐Specific ATE and RATE@5yr (95\\% CI)",
    label    = "tab:county_ate_rate5"
  ),
  con = "table_county_ate_rate5.tex"
)

writeLines(
  kable(
    big_latex,
    format   = "latex",
    booktabs = TRUE,
    linesep  = "",
    caption  = "CATE and CRATE@5yr for Big County (95\\% CI)",
    label    = "tab:CATE_CRATE_big"
  ),
  con = "table_CATE_CRATE_big.tex"
)

writeLines(
  kable(
    small_latex,
    format   = "latex",
    booktabs = TRUE,
    linesep  = "",
    caption  = "CATE and CRATE@5yr for Small County (95\\% CI)",
    label    = "tab:CATE_CRATE_small"
  ),
  con = "table_CATE_CRATE_small.tex"
)

writeLines(
  kable(
    spte5_latex,
    format   = "latex",
    booktabs = TRUE,
    linesep  = "",
    caption  = "County‐Specific SPTE@5yr (95\\% CI)",
    label    = "tab:county_spte"
  ),
  con = "table_county_spte.tex"
)

writeLines(
  kable(
    county_attatu_latex,
    format   = "latex",
    booktabs = TRUE,
    linesep  = "",
    caption  = "County‐Specific ATT and ATU (95\\% CI)",
    label    = "tab:county_att_atu"
  ),
  con = "table_county_att_atu.tex"
)

cat("All LaTeX tables have been printed to the console and written to .tex files.\n")
















