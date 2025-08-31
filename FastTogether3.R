library(SoftBart)
library(BART)
library(truncnorm)
library(MASS)
library(parallel)

n_replications <- 300
n_cores <- 10

cat("Number of cores available:", detectCores(), "\n")
cat("Using", n_cores, "cores for parallel processing\n")

create_dagar_components <- function(W, ordering, rho) {
  K <- nrow(W)
  N_pi <- vector("list", K)
  n_pi <- integer(K)

  for (i in 1:K) {
    node <- ordering[i]
    neighbors <- which(W[node, ] == 1)
    directed_neighbors <- neighbors[
      sapply(neighbors, function(j) which(ordering == j) < i)
    ]
    N_pi[[node]] <- directed_neighbors
    n_pi[node] <- length(directed_neighbors)
  }

  B <- matrix(0, K, K)
  for (i in 1:K) {
    node <- ordering[i]
    if (n_pi[node] > 0) {
      b_ij <- rho / (1 + (n_pi[node] - 1) * rho^2)
      B[node, N_pi[[node]]] <- b_ij
    }
  }

  tau_i <- (1 + (n_pi - 1)*rho^2) / (1 - rho^2)
  F_mat <- diag(tau_i)
  L <- diag(K) - B
  Q <- t(L) %*% F_mat %*% L

  return(list(B = B, F_mat = F_mat, L = L, Q = Q, N_pi = N_pi, n_pi = n_pi, ordering = ordering))
}

b2_func <- function(x_vals, z_val, w_val) {
  x1 <- x_vals[1]; x2 <- x_vals[2]
  x3 <- x_vals[3]; x4 <- x_vals[4]; x5 <- x_vals[5]

  main_effect <- 1.0 * z_val +
    0.7 * x1 - 0.5 * x2 + 0.3 * (z_val * x1 * exp(x2)) +
    0.5 * x3 - 0.2 * x4 + 0.1 * x5 +
    0.4 * x1 * x2 * exp(x2) +
    0.3 * x2 * x4 +
    0.2 * x1 * x3 * x4 +
    0.35 * z_val * x2 * x3 +
    0.25 * z_val * log(x1) * x4 +
    0.15 * z_val * x3 * x4 * x5 / 10 +
    0.5 * sin(pi * x5 / 5) +
    0.3 * sin(pi * x1 * x3) +
    0.2 * z_val * sin(pi * x2 * x4) +
    0.5 * z_val * exp(x1) * log(x2)

  covariate_w_interaction <- 0.3 * x1 * w_val +
    0.2 * x2 * w_val +
    0.15 * x3 * w_val +
    0.1 * (x5 - 5.5) * w_val +
    0.25 * z_val * x1 * w_val

  return(main_effect + covariate_w_interaction)
}

set.seed(2015)
K <- 35
cluster_sizes <- rep(50, K)

cluster_sizes[4] <- 100
cluster_sizes[35] <- 100

cluster_sizes[10] <- 25
cluster_sizes[20] <- 25

total_so_far <- sum(cluster_sizes)
remaining <- 1750 - total_so_far
other_counties <- setdiff(1:K, c(4, 35, 10, 20))
additional_per_county <- remaining / length(other_counties)
for (i in other_counties) {
  cluster_sizes[i] <- cluster_sizes[i] + floor(additional_per_county)
}
cluster_sizes[other_counties[1]] <- cluster_sizes[other_counties[1]] + (1750 - sum(cluster_sizes))

cat("Total sample size:", sum(cluster_sizes), "\n")
cat("County 4 size:", cluster_sizes[4], "\n")
cat("County 35 size:", cluster_sizes[35], "\n")
cat("County 10 size:", cluster_sizes[10], "\n")
cat("County 20 size:", cluster_sizes[20], "\n")

target_counties <- c(4, 35, 10, 20)

t_star_values <- c(2, 3, 4, 5, 15)

start_time <- Sys.time()

run_single_replication <- function(rep_num) {
  library(SoftBart)
  library(BART)
  library(truncnorm)
  library(MASS)

  b2_func_sim <- function(x_vals, z_val, w_val) {
    x1 <- x_vals[1]; x2 <- x_vals[2]
    x3 <- x_vals[3]; x4 <- x_vals[4]; x5 <- x_vals[5]

    main_effect <- 1.0 * z_val +
      0.7 * x1 - 0.5 * x2 + 0.3 * (z_val * x1 * exp(x2)) +
      0.5 * x3 - 0.2 * x4 + 0.1 * x5 +
      0.4 * x1 * x2 * exp(x2) +
      0.3 * x2 * x4 +
      0.2 * x1 * x3 * x4 +
      0.35 * z_val * x2 * x3 +
      0.25 * z_val * log(x1) * x4 +
      0.15 * z_val * x3 * x4 * x5 / 10 +
      0.5 * sin(pi * x5 / 5) +
      0.3 * sin(pi * x1 * x3) +
      0.2 * z_val * sin(pi * x2 * x4) +
      0.5 * z_val * exp(x1) * log(x2)

    covariate_w_interaction <- 0.3 * x1 * w_val +
      0.2 * x2 * w_val +
      0.15 * x3 * w_val +
      0.1 * (x5 - 5.5) * w_val +
      0.25 * z_val * x1 * w_val

    return(main_effect + covariate_w_interaction)
  }

  estimate_CRATE <- function(t_star, x, county_id,
                             e_hat,
                             forest_out,
                             w_samples, sigma2_samps,
                             burn_in, thin, n_mcmc) {

    sample_iters_out <- seq(from = burn_in + 1, to = n_mcmc, by = thin)
    Mprime <- length(sample_iters_out)
    CRATE_draws <- numeric(Mprime)

    for (m in seq_len(Mprime)) {
      x_out_1  <- matrix(c(x, 1, e_hat), nrow = 1)
      x_out_0  <- matrix(c(x, 0, e_hat), nrow = 1)

      iter_out <- sample_iters_out[m]
      b2_1_m   <- forest_out$predict_iteration(x_out_1, iter = iter_out)
      b2_0_m   <- forest_out$predict_iteration(x_out_0, iter = iter_out)
      W_i_m    <- w_samples[m, county_id]
      sigma2_m <- sigma2_samps[m]
      sigma_m  <- sqrt(sigma2_m)

      mu1      <- exp(b2_1_m + W_i_m + sigma2_m/2)
      mu0      <- exp(b2_0_m + W_i_m + sigma2_m/2)
      z1       <- (log(t_star) - b2_1_m - W_i_m) / sigma_m
      z0       <- (log(t_star) - b2_0_m - W_i_m) / sigma_m
      term1    <- mu1 * pnorm(z1-sigma_m)
      term2    <- t_star * (1 - pnorm(z1))
      term3    <- mu0 * pnorm(z0-sigma_m)
      term4    <- t_star * (1 - pnorm(z0))

      CRATE_draws[m] <- (term1 + term2) - (term3 + term4)
    }

    list(
      CRATE_mean  = mean(CRATE_draws),
      CRATE_lower = quantile(CRATE_draws, 0.025),
      CRATE_upper = quantile(CRATE_draws, 0.975),
      CRATE_draws = CRATE_draws
    )
  }

  estimate_County_RATE <- function(t_star, X_i, county_id,
                                   e_hat_vec,
                                   forest_out,
                                   w_samples, sigma2_samps,
                                   burn_in, thin, n_mcmc) {

    n_i <- nrow(X_i)
    sample_iters_out <- seq(from = burn_in + 1, to = n_mcmc, by = thin)
    Mprime <- length(sample_iters_out)
    CRATE_matrix <- matrix(0, nrow = Mprime, ncol = n_i)

    for (j in seq_len(n_i)) {
      x_j   <- X_i[j, ]
      e_hat_j <- e_hat_vec[j]

      tmp_j <- estimate_CRATE(
        t_star        = t_star,
        x             = x_j,
        county_id     = county_id,
        e_hat         = e_hat_j,
        forest_out    = forest_out,
        w_samples     = w_samples,
        sigma2_samps  = sigma2_samps,
        burn_in       = burn_in,
        thin          = thin,
        n_mcmc        = n_mcmc
      )
      CRATE_matrix[, j] <- tmp_j$CRATE_draws
    }

    County_RATE_draws <- rowMeans(CRATE_matrix)

    list(
      County_RATE_mean  = mean(County_RATE_draws),
      County_RATE_lower = quantile(County_RATE_draws, 0.025),
      County_RATE_upper = quantile(County_RATE_draws, 0.975),
      County_RATE_draws = County_RATE_draws
    )
  }

  cat("\n========================================\n")
  cat("Starting Replication", rep_num, "of", n_replications, "\n")
  cat("========================================\n")

  set.seed(2015 + rep_num - 1)

  cat("Generating common data...\n")

  K <- 35
  nrow_grid <- 5; ncol_grid <- 7
  A <- matrix(0, nrow = K, ncol = K)

  for (i in 1:K) {
    r_i <- ceiling(i / ncol_grid)
    c_i <- i - (r_i - 1) * ncol_grid
    neighbors <- list(
      c(r_i - 1, c_i),
      c(r_i + 1, c_i),
      c(r_i, c_i - 1),
      c(r_i, c_i + 1)
    )
    for (nb in neighbors) {
      if (nb[[1]] >= 1 && nb[[1]] <= nrow_grid &&
          nb[[2]] >= 1 && nb[[2]] <= ncol_grid) {
        j <- (nb[[1]] - 1) * ncol_grid + nb[[2]]
        A[i, j] <- 1
      }
    }
  }

  A <- (A + t(A)) > 0
  A <- 1 * A

  rho_true <- 0.7
  tau_w_true <- 2.0
  dagar_true <- create_dagar_components(W = A, ordering = 1:K, rho = rho_true)
  Sigma_W <- solve(t(dagar_true$L) %*% dagar_true$F_mat %*% dagar_true$L) / tau_w_true
  W_true <- as.numeric(mvrnorm(1, mu = rep(0, K), Sigma = Sigma_W))

  n <- sum(cluster_sizes)
  cluster_id <- rep(1:K, times = cluster_sizes)
  X1 <- runif(n, 0, 1)
  X2 <- rbeta(n, 2, 6)
  X3 <- rbinom(n, 1, 0.8)
  X4 <- rbinom(n, 1, 0.7)
  X5 <- sample(1:10, n, replace = TRUE)
  X <- cbind(X1, X2, X3, X4, X5)

  b1_lin <- 0.5 * X1 - 0.3 * X2 + 0.2 * X3 - 0.1 * X4 + 0.05 * X2*(X5 - 5.5)+ X1*X2
  e_true <- pnorm(b1_lin)
  Z <- rbinom(n, 1, prob = e_true)

  b2_true_log <- numeric(n)
  for (i in 1:n) {
    b2_true_log[i] <- b2_func_sim(X[i,], Z[i], W_true[cluster_id[i]])
  }

  sigma2_true <- 0.15
  logT_true <- b2_true_log + rnorm(n, 0, sqrt(sigma2_true))
  T_true <- exp(logT_true)

  C_time <- rexp(n, rate = 0.1)
  time_obs <- pmin(T_true, C_time)
  status <- as.integer(T_true <= C_time)

  V <- matrix(0, nrow = n, ncol = 1)
  V_test <- V

  cat("Computing true CRATE for 4 counties at 3 time points...\n")

  true_values <- expand.grid(
    county = target_counties,
    t_star = t_star_values,
    stringsAsFactors = FALSE
  )
  true_values$true_CRATE <- numeric(nrow(true_values))

  for (row in 1:nrow(true_values)) {
    i <- true_values$county[row]
    t_star <- true_values$t_star[row]

    idx_i <- which(cluster_id == i)
    X_i <- X[idx_i, , drop = FALSE]
    W_i_true <- W_true[i]
    n_i <- nrow(X_i)

    b2_1_i <- numeric(n_i)
    b2_0_i <- numeric(n_i)

    for (j in 1:n_i) {
      x_j <- X_i[j, ]
      b2_1_i[j] <- b2_func_sim(x_j, 1, W_i_true)
      b2_0_i[j] <- b2_func_sim(x_j, 0, W_i_true)
    }

    sigma_m <- sqrt(sigma2_true)
    CRATE_i_j <- numeric(n_i)

    for (j in 1:n_i) {
      z1_j <- (log(t_star) - b2_1_i[j]) / sigma_m
      z0_j <- (log(t_star) - b2_0_i[j]) / sigma_m

      mu1_j <- exp(b2_1_i[j] + sigma2_true / 2)
      mu0_j <- exp(b2_0_i[j] + sigma2_true / 2)

      term1_j <- mu1_j * pnorm(z1_j - sigma_m)
      term2_j <- t_star * (1 - pnorm(z1_j))
      term3_j <- mu0_j * pnorm(z0_j - sigma_m)
      term4_j <- t_star * (1 - pnorm(z0_j))

      CRATE_i_j[j] <- (term1_j + term2_j) - (term3_j + term4_j)
    }

    true_values$true_CRATE[row] <- mean(CRATE_i_j)
  }

  cat("True CRATE values computed:\n")
  print(true_values)

  results_all_methods <- list()

  cat("\n--- METHOD 1: Fitting AFT-DAGAR-SOFTBART model ---\n")

  burn_in_dagar <- 50000
  thin_dagar <- 5
  n_mcmc_dagar <- 90000
  n_iter_ps <-90000

  fit_dagar <- AFT_mixed_DAGAR_probit_causal3(
    time = time_obs,
    status = status,
    X = X,
    Z = Z,
    group = cluster_id,
    X_test = X,
    Z_test = Z,
    group_test = cluster_id,
    num_tree_ps = 20,
    num_tree_out = 20,
    n_iter_ps = n_iter_ps,
    n_mcmc = n_mcmc_dagar,
    burn_in = burn_in_dagar,
    thin = thin_dagar,
    rho_init = 0.5,
    rho_prop_sd = 0.02,
    a_tau = 1,
    b_tau = 1,
    a_sigma = 2,
    b_sigma = 1
  )

  forest_ps <- fit_dagar$probit_forest
  forest_out_dagar <- fit_dagar$forest_out
  w_samples_dagar <- fit_dagar$w_samples
  sigma2_samps_dagar <- fit_dagar$sigma2_samps

  dagar_results <- expand.grid(
    county = target_counties,
    t_star = t_star_values,
    method = "AFT-DAGAR",
    stringsAsFactors = FALSE
  )

  dagar_results$true_CRATE <- true_values$true_CRATE
  dagar_results$est_CRATE <- numeric(nrow(dagar_results))
  dagar_results$CI_lower_CRATE <- numeric(nrow(dagar_results))
  dagar_results$CI_upper_CRATE <- numeric(nrow(dagar_results))
  dagar_results$MAE_CRATE <- numeric(nrow(dagar_results))
  dagar_results$coverage_CRATE <- logical(nrow(dagar_results))
  dagar_results$CI_width_CRATE <- numeric(nrow(dagar_results))

  e_hat_train <- apply(pnorm(fit_dagar$r_train), 2, mean)

  for (row in 1:nrow(dagar_results)) {
    i <- dagar_results$county[row]
    t_star <- dagar_results$t_star[row]

    idx_i <- which(cluster_id == i)
    X_i <- X[idx_i, , drop = FALSE]

    e_hat_vec_i <- e_hat_train[idx_i]

    est_crate <- estimate_County_RATE(
      t_star = t_star,
      X_i = X_i,
      county_id = i,
      e_hat_vec = e_hat_vec_i,
      forest_out = forest_out_dagar,
      w_samples = w_samples_dagar,
      sigma2_samps = sigma2_samps_dagar,
      burn_in = burn_in_dagar,
      thin = thin_dagar,
      n_mcmc = n_mcmc_dagar
    )

    dagar_results$est_CRATE[row] <- est_crate$County_RATE_mean
    dagar_results$CI_lower_CRATE[row] <- est_crate$County_RATE_lower
    dagar_results$CI_upper_CRATE[row] <- est_crate$County_RATE_upper
    dagar_results$MAE_CRATE[row] <- abs(est_crate$County_RATE_mean - dagar_results$true_CRATE[row])
    dagar_results$coverage_CRATE[row] <- (dagar_results$true_CRATE[row] >= est_crate$County_RATE_lower &&
                                            dagar_results$true_CRATE[row] <= est_crate$County_RATE_upper)
    dagar_results$CI_width_CRATE[row] <- est_crate$County_RATE_upper - est_crate$County_RATE_lower
  }

  cat("\n--- METHOD 2: Fitting riAFT (Simplified AFT) model ---\n")

  burn_in_riaft <- 200
  thin_riaft <- 5
  n_mcmc_riaft <- 5000

  fit_riaft <- AFT_mixed_simple_causal(
    time = time_obs,
    status = status,
    X = X,
    V = V,
    Z = Z,
    group = cluster_id,
    X_test = X,
    V_test = V,
    Z_test = Z,
    group_test = cluster_id,
    num_tree_out = 20,
    n_mcmc = n_mcmc_riaft,
    burn_in = burn_in_riaft,
    thin = thin_riaft,
    a_tau = 1,
    b_tau = 1,
    a_sigma = 2,
    b_sigma = 1
  )

  forest_out_riaft <- fit_riaft$forest_out
  w_samples_riaft <- fit_riaft$w_samples
  sigma2_samps_riaft <- fit_riaft$sigma2_samps

  riaft_results <- expand.grid(
    county = target_counties,
    t_star = t_star_values,
    method = "riAFT",
    stringsAsFactors = FALSE
  )

  riaft_results$true_CRATE <- true_values$true_CRATE
  riaft_results$est_CRATE <- numeric(nrow(riaft_results))
  riaft_results$CI_lower_CRATE <- numeric(nrow(riaft_results))
  riaft_results$CI_upper_CRATE <- numeric(nrow(riaft_results))
  riaft_results$MAE_CRATE <- numeric(nrow(riaft_results))
  riaft_results$coverage_CRATE <- logical(nrow(riaft_results))
  riaft_results$CI_width_CRATE <- numeric(nrow(riaft_results))

  for (row in 1:nrow(riaft_results)) {
    i <- riaft_results$county[row]
    t_star <- riaft_results$t_star[row]

    idx_i <- which(cluster_id == i)
    X_i <- X[idx_i, , drop = FALSE]

    est_rate <- estimate_County_RATE_simple(
      t_star = t_star,
      X_i = X_i,
      v_i = c(0),
      county_id = i,
      forest_out = forest_out_riaft,
      w_samples = w_samples_riaft,
      sigma2_samps = sigma2_samps_riaft,
      burn_in = burn_in_riaft,
      thin = thin_riaft,
      n_mcmc = n_mcmc_riaft
    )

    riaft_results$est_CRATE[row] <- est_rate$County_RATE_mean
    riaft_results$CI_lower_CRATE[row] <- est_rate$County_RATE_lower
    riaft_results$CI_upper_CRATE[row] <- est_rate$County_RATE_upper
    riaft_results$MAE_CRATE[row] <- abs(est_rate$County_RATE_mean - riaft_results$true_CRATE[row])
    riaft_results$coverage_CRATE[row] <- (riaft_results$true_CRATE[row] >= est_rate$County_RATE_lower &&
                                            riaft_results$true_CRATE[row] <= est_rate$County_RATE_upper)
    riaft_results$CI_width_CRATE[row] <- est_rate$County_RATE_upper - est_rate$County_RATE_lower
  }

  cat("\n--- METHOD 3: Fitting PS-PH-Frailty model ---\n")

  bayes_fit <- bayesian_weibull_frailty_ps_causal(
    time = time_obs,
    event = status,
    treatment = Z,
    covariates = X,
    cluster = cluster_id,
    n_iter = 2000,
    n_burn = 1000,
    n_iter_ps = 1000,
    n_burn_ps = 500
  )

  propensity_scores <- bayes_fit$propensity_score

  frailty_results <- expand.grid(
    county = target_counties,
    t_star = t_star_values,
    method = "PS-PH-Frailty",
    stringsAsFactors = FALSE
  )

  frailty_results$true_CRATE <- true_values$true_CRATE
  frailty_results$est_CRATE <- numeric(nrow(frailty_results))
  frailty_results$CI_lower_CRATE <- numeric(nrow(frailty_results))
  frailty_results$CI_upper_CRATE <- numeric(nrow(frailty_results))
  frailty_results$MAE_CRATE <- numeric(nrow(frailty_results))
  frailty_results$coverage_CRATE <- logical(nrow(frailty_results))
  frailty_results$CI_width_CRATE <- numeric(nrow(frailty_results))

  for (row in 1:nrow(frailty_results)) {
    i <- frailty_results$county[row]
    t_star <- frailty_results$t_star[row]

    idx_i <- which(cluster_id == i)
    X_i <- X[idx_i, , drop = FALSE]

    est_cnty_RATE <- County_RATE_weibull(
      bayes_fit = bayes_fit,
      tau = t_star,
      x = X_i,
      cluster = rep(i, nrow(X_i)),
      propensity_scores = propensity_scores[idx_i],
      credible_level = 0.95
    )

    frailty_results$est_CRATE[row] <- est_cnty_RATE$County_RATE_mean
    frailty_results$CI_lower_CRATE[row] <- est_cnty_RATE$County_RATE_lower
    frailty_results$CI_upper_CRATE[row] <- est_cnty_RATE$County_RATE_upper
    frailty_results$MAE_CRATE[row] <- abs(est_cnty_RATE$County_RATE_mean - frailty_results$true_CRATE[row])
    frailty_results$coverage_CRATE[row] <- (frailty_results$true_CRATE[row] >= est_cnty_RATE$County_RATE_lower &&
                                              frailty_results$true_CRATE[row] <= est_cnty_RATE$County_RATE_upper)
    frailty_results$CI_width_CRATE[row] <- est_cnty_RATE$County_RATE_upper - est_cnty_RATE$County_RATE_lower
  }

  combined_results <- rbind(dagar_results, riaft_results, frailty_results)

  cat("\nReplication", rep_num, "completed.\n")

  return(combined_results)
}

cat("\n========================================\n")
cat("Running", n_replications, "replications in parallel on", n_cores, "cores\n")
cat("========================================\n")

all_results <- mclapply(
  1:n_replications,
  run_single_replication,
  mc.cores = n_cores,
  mc.preschedule = FALSE
)

errors <- sapply(all_results, function(x) inherits(x, "try-error"))
if (any(errors)) {
  cat("Warning: The following replications encountered errors:\n")
  cat(which(errors), "\n")
  all_results <- all_results[!errors]
}

cat("\n========================================\n")
cat("Aggregating and comparing results\n")
cat("========================================\n")

all_results_df <- do.call(rbind,
                          lapply(seq_along(all_results), function(rep) {
                            df <- all_results[[rep]]
                            df$replication <- rep
                            return(df)
                          }))

summary_results <- aggregate(
  cbind(MAE_CRATE, coverage_CRATE, CI_width_CRATE) ~ method + county + t_star,
  data = all_results_df,
  FUN = function(x) c(mean = mean(x), se = sd(x)/sqrt(length(x)))
)

summary_clean <- data.frame(
  method = rep(rep(c("AFT-DAGAR", "riAFT", "PS-PH-Frailty"), each = length(target_counties)), length(t_star_values)),
  county = rep(rep(target_counties, 3), length(t_star_values)),
  t_star = rep(t_star_values, each = length(target_counties) * 3),
  MAE_CRATE_mean = summary_results[, "MAE_CRATE"][, "mean"],
  MAE_CRATE_se = summary_results[, "MAE_CRATE"][, "se"],
  coverage_CRATE = summary_results[, "coverage_CRATE"][, "mean"] * 100,
  CI_width_CRATE = summary_results[, "CI_width_CRATE"][, "mean"]
)

write.csv(all_results_df, "Sim3_combined_simulation_detailed_results.csv", row.names = FALSE)
write.csv(summary_clean, "Sim3_combined_simulation_summary.csv", row.names = FALSE)

cat("\nDetailed results saved to: Sim3_combined_simulation_detailed_results.csv\n")
cat("Summary results saved to: Sim3_combined_simulation_summary.csv\n")

cat("\n========================================\n")
cat("COMPARISON OF THREE METHODS FOR 4 COUNTIES AT 3 TIME POINTS\n")
cat("========================================\n")

print(summary_clean)

if (n_replications > 1) {
  par(mfrow = c(3, 4))

  for (t_idx in 1:length(t_star_values)) {
    t_star <- t_star_values[t_idx]

    for (county in target_counties) {
      boxplot(MAE_CRATE ~ method,
              data = all_results_df[all_results_df$county == county & all_results_df$t_star == t_star,],
              main = paste("County", county, "- t* =", t_star),
              ylab = "MAE CRATE")
    }
  }
}

end_time <- Sys.time()
total_time <- difftime(end_time, start_time, units = "hours")
cat("\nTotal simulation time:", round(total_time, 2), "hours\n")
cat("Average time per replication:", round(as.numeric(total_time) * 60 / n_replications, 2), "minutes\n")

cat("\n========================================\n")
cat("METHOD COMPARISON SUMMARY\n")
cat("========================================\n")

for (t_star in t_star_values) {
  cat("\n\nTime point t* =", t_star, "\n")
  cat("================\n")

  for (county in target_counties) {
    cat("\nCounty", county, "(n =", cluster_sizes[county], "):\n")
    cat("----------\n")

    dagar_data <- summary_clean[summary_clean$method == "AFT-DAGAR" &
                                  summary_clean$county == county &
                                  summary_clean$t_star == t_star, ]
    riaft_data <- summary_clean[summary_clean$method == "riAFT" &
                                  summary_clean$county == county &
                                  summary_clean$t_star == t_star, ]
    frailty_data <- summary_clean[summary_clean$method == "PS-PH-Frailty" &
                                    summary_clean$county == county &
                                    summary_clean$t_star == t_star, ]

    cat("CRATE Performance:\n")
    cat("  AFT-DAGAR MAE:    ", round(dagar_data$MAE_CRATE_mean, 4),
        "(Coverage:", round(dagar_data$coverage_CRATE, 1), "%)\n")
    cat("  riAFT MAE:        ", round(riaft_data$MAE_CRATE_mean, 4),
        "(Coverage:", round(riaft_data$coverage_CRATE, 1), "%)\n")
    cat("  PS-PH-Frailty MAE:", round(frailty_data$MAE_CRATE_mean, 4),
        "(Coverage:", round(frailty_data$coverage_CRATE, 1), "%)\n")
  }
}

save(all_results, file = "Sim3_three_methods_sim_all_results_list.RData")
cat("\nAll results list saved to: Sim3_three_methods_sim_all_results_list.RData\n")
