# ================================
# Functions.R (Modified with BART pbart)
# ================================
#
# This file contains:
#   1. create_dagar_components()
#   2. AFT_mixed_DAGAR_probit_causal()
#   3. estimate_CRATE()
#   4. estimate_County_RATE()
#
# Uses BART::pbart for propensity scores and SoftBart for outcome model

library(SoftBart)
library(truncnorm)
library(MASS)
library(BART)

# --------------------------------------------------
# 1. create_dagar_components()
# --------------------------------------------------
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
    n_pi[node]   <- length(directed_neighbors)
  }
  B     <- matrix(0, K, K)
  for (i in 1:K) {
    node <- ordering[i]
    if (n_pi[node] > 0) {
      b_ij <- rho / (1 + (n_pi[node] - 1) * rho^2)
      B[node, N_pi[[node]]] <- b_ij
    }
  }
  tau_i  <- (1 + (n_pi - 1)*rho^2) / (1 - rho^2)
  F_mat  <- diag(tau_i)
  L      <- diag(K) - B
  Q      <- t(L) %*% F_mat %*% L
  return(list(B = B, F_mat = F_mat, L = L, Q = Q, N_pi = N_pi, n_pi = n_pi, ordering = ordering))
}

# --------------------------------------------------
# 2. AFT_mixed_DAGAR_probit_causal()
# --------------------------------------------------
AFT_mixed_DAGAR_probit_causal2 <- function(
    time, status, X, Z, group,
    X_test = NULL, Z_test = NULL, group_test = NULL,
    num_tree_ps = 50, num_tree_out = 50,
    n_iter_ps = 200, n_mcmc = 2000, burn_in = 500, thin = 5,
    rho_init = 0.5, rho_prop_sd = 0.02,
    a_tau = 1, b_tau = 1, a_sigma = 2, b_sigma = 1,  a_rho = 2, b_rho = 2
) {
  n <- length(time)
  K <- length(unique(group))
  if (is.null(Z_test))       Z_test    <- Z
  if (is.null(group_test))   group_test <- group
  if (is.null(X_test))       X_test    <- X

  cluster_indices <- split(1:n, group)
  n_i_vec <- sapply(cluster_indices, length)

  # Stage 1: Probit BART using BART::pbart
  X_train   <- X
  X_test_ps <- X_test

  # Run pbart for propensity score model
  probit_fit <- pbart(x.train = X_train,
                      y.train = Z,
                      x.test = X_test_ps,
                      ntree = num_tree_ps,
                      ndpost = n_iter_ps,
                      nskip = 100,
                      printevery = 100)

  # Extract propensity scores (apply pnorm to get probabilities)
  e_hat_train <- apply(pnorm(probit_fit$yhat.train), 2, mean)
  e_hat_test  <- apply(pnorm(probit_fit$yhat.test), 2, mean)

  # Stage 2: Outcome SBART + DAGAR
  y_log <- log(time)
  X_out_train <- cbind(X, Z, e_hat_train)
  X_out_test  <- cbind(X_test, Z_test, e_hat_test)

  hypers_out <- Hypers(X = X_out_train, Y = y_log, k = 2,
                       num_tree = num_tree_out, sigma_hat = sd(y_log))
  opts_out   <- Opts(update_sigma = FALSE, cache_trees = TRUE)
  forest_out <- MakeForest(hypers_out, opts_out)

  f_current      <- forest_out$do_predict(X_out_train)
  f_test_current <- forest_out$do_predict(X_out_test)

  sigma2_current <- 1
  tau_w_current  <- 1
  rho_current    <- rho_init
  W_current      <- rep(0, K)

  n_save <- floor((n_mcmc - burn_in) / thin)
  f_samples      <- matrix(0, nrow = n_save, ncol = n)
  f_test_samples <- matrix(0, nrow = n_save, ncol = nrow(X_out_test))
  w_samples      <- matrix(0, nrow = n_save, ncol = K)
  sigma2_samples <- numeric(n_save)
  tau_w_samples  <- numeric(n_save)
  rho_samples    <- numeric(n_save)

  sample_idx <- 1
  tilde_y    <- numeric(n)

  # Initial imputation for censoring
  mu_init <- f_current + W_current[group]
  for (j in 1:n) {
    if (status[j] == 1) {
      tilde_y[j] <- y_log[j]
    } else {
      tilde_y[j] <- truncnorm::rtruncnorm(
        n    = 1,
        a    = y_log[j],
        b    = Inf,
        mean = mu_init[j],
        sd   = sqrt(sigma2_current)
      )
    }
  }

  # MCMC loop
  for (iter in 1:n_mcmc) {
    # 1. Impute censored log‐time
    mu_i <- f_current + W_current[group]
    for (j in 1:n) {
      if (status[j] == 1) {
        tilde_y[j] <- y_log[j]
      } else {
        tilde_y[j] <- truncnorm::rtruncnorm(
          n    = 1,
          a    = y_log[j],
          b    = Inf,
          mean = mu_i[j],
          sd   = sqrt(sigma2_current)
        )
      }
    }

    # 2. Update outcome SBART
    w_obs      <- W_current[group]
    resid_for_f <- tilde_y - w_obs
    forest_out$do_gibbs(X_out_train, resid_for_f, X_out_test, 1)
    f_current      <- forest_out$do_predict(X_out_train)
    f_test_current <- forest_out$do_predict(X_out_test)

    # 3. Update DAGAR random effects W
    A <- matrix(0, nrow = K, ncol = K)
    for (ii in 1:(K - 1)) {
      A[ii, ii + 1] <- 1; A[ii + 1, ii] <- 1
    }
    dagar     <- create_dagar_components(W = A, ordering = 1:K, rho = rho_current)
    Q_current <- dagar$Q
    b_vec     <- numeric(K)
    for (i in 1:K) {
      idx_i   <- cluster_indices[[i]]
      b_vec[i] <- sum(tilde_y[idx_i] - f_current[idx_i]) / sigma2_current
    }
    Lik_prec <- diag(n_i_vec) / sigma2_current
    Prec_w   <- Lik_prec + tau_w_current * Q_current
    Cov_w    <- solve(Prec_w)
    mu_w     <- Cov_w %*% b_vec
    W_current <- as.numeric(mvrnorm(1, mu = mu_w, Sigma = Cov_w))

    # 4. Update sigma^2
    w_obs    <- W_current[group]
    resid_all <- tilde_y - f_current - w_obs
    shape_s   <- a_sigma + n / 2
    rate_s    <- b_sigma + 0.5 * sum(resid_all^2)
    inv_s     <- rgamma(1, shape = shape_s, rate = rate_s)
    sigma2_current <- 1 / inv_s

    # 5. Update tau_w
    shape_t   <- a_tau + K / 2
    rate_t    <- b_tau + 0.5 * as.numeric(t(W_current) %*% Q_current %*% W_current)
    tau_w_current <- rgamma(1, shape = shape_t, rate = rate_t)

    # 6. Update rho by MH with Beta prior
    rho_prop <- rho_current + rnorm(1, mean = 0, sd = rho_prop_sd)

    if (rho_prop > 0 && rho_prop < 1) {
      dager_prop <- create_dagar_components(W = A, ordering = 1:K, rho = rho_prop)
      Q_prop <- dager_prop$Q

      log_det_prop <- as.numeric(determinant(Q_prop, logarithm = TRUE)$modulus)
      log_det_curr <- as.numeric(determinant(Q_current, logarithm = TRUE)$modulus)

      # Log-likelihood part
      log_lik_prop <- 0.5 * log_det_prop -
        0.5 * tau_w_current * as.numeric(t(W_current) %*% Q_prop %*% W_current)
      log_lik_curr <- 0.5 * log_det_curr -
        0.5 * tau_w_current * as.numeric(t(W_current) %*% Q_current %*% W_current)

      # Add Beta prior log-density
      log_prior_prop <- (a_rho - 1) * log(rho_prop) + (b_rho - 1) * log(1 - rho_prop)
      log_prior_curr <- (a_rho - 1) * log(rho_current) + (b_rho - 1) * log(1 - rho_current)

      # Total log-posterior
      log_target_prop <- log_lik_prop + log_prior_prop
      log_target_curr <- log_lik_curr + log_prior_curr

      if (log(runif(1)) < (log_target_prop - log_target_curr)) {
        rho_current <- rho_prop
        Q_current <- Q_prop
      }
    }
    # 7. Store draws after burn‐in & thinning
    if (iter > burn_in && ((iter - burn_in) %% thin == 0)) {
      f_samples[sample_idx, ]      <- f_current
      f_test_samples[sample_idx, ] <- f_test_current
      w_samples[sample_idx, ]      <- W_current
      sigma2_samples[sample_idx]   <- sigma2_current
      tau_w_samples[sample_idx]    <- tau_w_current
      rho_samples[sample_idx]      <- rho_current
      sample_idx <- sample_idx + 1
    }
  }

  return(list(
    probit_forest = probit_fit,
    r_train       = probit_fit$yhat.train,
    r_test        = probit_fit$yhat.test,
    forest_out    = forest_out,
    f_train       = f_samples,
    f_test        = f_test_samples,
    w_samples     = w_samples,
    sigma2_samps  = sigma2_samples,
    tau_w_samps   = tau_w_samples,
    rho_samps     = rho_samples
  ))
}

# --------------------------------------------------
# 3. estimate_CRATE()
# --------------------------------------------------
estimate_CRATE <- function(t_star, x, county_id,
                           forest_ps, forest_out,
                           w_samples, sigma2_samps,
                           burn_in, thin, n_mcmc, n_iter_ps) {
  sample_iters_out <- seq(from = burn_in + 1, to = n_mcmc, by = thin)
  Mprime <- length(sample_iters_out)
  CRATE_draws <- numeric(Mprime)
  x_ps_mat   <- matrix(x, nrow = 1)

  # Get propensity score posterior draws
  # Since we're using pbart, we need to predict for this specific x
  # predict returns a matrix with posterior draws
  ps_pred <- predict(forest_ps, newdata = x_ps_mat)
  # ps_pred is a matrix where each row is a posterior draw
  # Take the mean of the probabilities across all posterior draws
  e_hat_mean <- mean(pnorm(ps_pred$yhat.test))

  for (m in seq_len(Mprime)) {
    x_out_1  <- matrix(c(x, 1, e_hat_mean), nrow = 1)
    x_out_0  <- matrix(c(x, 0, e_hat_mean), nrow = 1)
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

# --------------------------------------------------
# 4. estimate_County_RATE()
# --------------------------------------------------
estimate_County_RATE <- function(t_star, X_i, county_id,
                                 forest_ps, forest_out,
                                 w_samples, sigma2_samps,
                                 burn_in, thin, n_mcmc, n_iter_ps) {
  n_i <- nrow(X_i)
  sample_iters_out <- seq(from = burn_in + 1, to = n_mcmc, by = thin)
  Mprime <- length(sample_iters_out)
  CRATE_matrix <- matrix(0, nrow = Mprime, ncol = n_i)
  for (j in seq_len(n_i)) {
    x_j   <- X_i[j, ]
    tmp_j <- estimate_CRATE(
      t_star        = t_star,
      x             = x_j,
      county_id     = county_id,
      forest_ps     = forest_ps,
      forest_out    = forest_out,
      w_samples     = w_samples,
      sigma2_samps  = sigma2_samps,
      burn_in       = burn_in,
      thin          = thin,
      n_mcmc        = n_mcmc,
      n_iter_ps     = n_iter_ps
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

# Helper Functions for Simulation 2

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

#  b2 function WITH COVARIATE-RANDOM EFFECT INTERACTION
b2_func <- function(x_vals, z_val, w_val) {
  x1 <- x_vals[1]; x2 <- x_vals[2]
  x3 <- x_vals[3]; x4 <- x_vals[4]; x5 <- x_vals[5]


  main_effect <- 1.0 * z_val +
    0.7 * x1 - 0.5 * x2 + 0.3 * (z_val * x1) +
    0.5 * x3 - 0.2 * x4 + 0.1 * x5 +
    
    0.4 * x1 * x2 * x3 +
    0.3 * x2 * x4 +
    0.2 * x1 * x3 * x4 +
    
    0.35 * z_val * x2 * x3 +
    0.25 * z_val * x1 * x4 +
    0.15 * z_val * x3 * x4 * x5 / 10 +
 
    0.5 * sin(pi * x5 / 5) +
    0.3 * sin(pi * x1 * x3) +
    0.2 * z_val * sin(pi * x2 * x4)

  covariate_w_interaction <- 0.3 * x1 * w_val +
    0.2 * x2 * w_val +
    0.15 * x3 * w_val +
    0.1 * (x5 - 5.5) * w_val / 10 +
    0.25 * z_val * x1 * w_val

  return(main_effect + covariate_w_interaction)
}

# AFT mixed model 
AFT_mixed_simple_causal <- function(
    time, status, X, V, Z, group,
    X_test = NULL, V_test = NULL, Z_test = NULL, group_test = NULL,
    num_tree_out = 50,
    n_mcmc = 2000, burn_in = 500, thin = 5,
    a_tau = 1, b_tau = 1, a_sigma = 2, b_sigma = 1
) {
  n <- length(time)
  K <- length(unique(group))

  if (is.null(Z_test)) Z_test <- Z
  if (is.null(group_test)) group_test <- group
  if (is.null(X_test)) X_test <- X
  if (is.null(V_test)) V_test <- V

  cluster_indices <- split(1:n, group)
  n_i_vec <- sapply(cluster_indices, length)

  # Single stage: Outcome SBART only
  y_log <- log(time)
  X_out_train <- cbind(X, V, Z)
  X_out_test <- cbind(X_test, V_test, Z_test)

  hypers_out <- Hypers(X = X_out_train, Y = y_log, k = 2,
                       num_tree = num_tree_out, sigma_hat = sd(y_log))
  opts_out <- Opts(update_sigma = FALSE, cache_trees = TRUE)
  forest_out <- MakeForest(hypers_out, opts_out)

  f_current <- forest_out$do_predict(X_out_train)
  f_test_current <- forest_out$do_predict(X_out_test)

  sigma2_current <- 1
  tau_w_current <- 1
  W_current <- rep(0, K)

  n_save <- floor((n_mcmc - burn_in) / thin)
  f_samples <- matrix(0, nrow = n_save, ncol = n)
  f_test_samples <- matrix(0, nrow = n_save, ncol = nrow(X_out_test))
  w_samples <- matrix(0, nrow = n_save, ncol = K)
  sigma2_samples <- numeric(n_save)
  tau_w_samples <- numeric(n_save)

  sample_idx <- 1
  tilde_y <- numeric(n)

  # Initial imputation for censoring
  mu_init <- f_current + W_current[group]
  for (j in 1:n) {
    if (status[j] == 1) {
      tilde_y[j] <- y_log[j]
    } else {
      tilde_y[j] <- truncnorm::rtruncnorm(
        n = 1,
        a = y_log[j],
        b = Inf,
        mean = mu_init[j],
        sd = sqrt(sigma2_current)
      )
    }
  }

  # MCMC loop
  for (iter in 1:n_mcmc) {
    # 1. Impute censored log-time
    mu_i <- f_current + W_current[group]
    for (j in 1:n) {
      if (status[j] == 1) {
        tilde_y[j] <- y_log[j]
      } else {
        tilde_y[j] <- truncnorm::rtruncnorm(
          n = 1,
          a = y_log[j],
          b = Inf,
          mean = mu_i[j],
          sd = sqrt(sigma2_current)
        )
      }
    }

    # 2. Update outcome SBART
    w_obs <- W_current[group]
    resid_for_f <- tilde_y - w_obs
    forest_out$do_gibbs(X_out_train, resid_for_f, X_out_test, 1)
    f_current <- forest_out$do_predict(X_out_train)
    f_test_current <- forest_out$do_predict(X_out_test)

    # 3. Update IID normal random effects W
    for (i in 1:K) {
      idx_i <- cluster_indices[[i]]
      sum_resid_i <- sum(tilde_y[idx_i] - f_current[idx_i])
      prec_w_i <- n_i_vec[i] / sigma2_current + tau_w_current
      var_w_i <- 1 / prec_w_i
      mu_w_i <- var_w_i * (sum_resid_i / sigma2_current)
      W_current[i] <- rnorm(1, mean = mu_w_i, sd = sqrt(var_w_i))
    }

    # 4. Update sigma^2
    w_obs <- W_current[group]
    resid_all <- tilde_y - f_current - w_obs
    shape_s <- a_sigma + n / 2
    rate_s <- b_sigma + 0.5 * sum(resid_all^2)
    inv_s <- rgamma(1, shape = shape_s, rate = rate_s)
    sigma2_current <- 1 / inv_s

    # 5. Update tau_w
    shape_t <- a_tau + K / 2
    rate_t <- b_tau + 0.5 * sum(W_current^2)
    tau_w_current <- rgamma(1, shape = shape_t, rate = rate_t)

    # 6. Store draws after burn-in & thinning
    if (iter > burn_in && ((iter - burn_in) %% thin == 0)) {
      f_samples[sample_idx, ] <- f_current
      f_test_samples[sample_idx, ] <- f_test_current
      w_samples[sample_idx, ] <- W_current
      sigma2_samples[sample_idx] <- sigma2_current
      tau_w_samples[sample_idx] <- tau_w_current
      sample_idx <- sample_idx + 1
    }
  }

  return(list(
    forest_out = forest_out,
    f_train = f_samples,
    f_test = f_test_samples,
    w_samples = w_samples,
    sigma2_samps = sigma2_samples,
    tau_w_samps = tau_w_samples
  ))
}

# County-ATE estimation
estimate_County_ATE_simple <- function(X_i, v_i, county_id,
                                       forest_out, w_samples, sigma2_samps,
                                       burn_in, thin, n_mcmc) {
  n_i <- nrow(X_i)
  n_draws <- nrow(w_samples)
  County_ATE_draws <- numeric(n_draws)

  for (m in seq_len(n_draws)) {
    exp1_vec <- numeric(n_i)
    exp0_vec <- numeric(n_i)

    for (j in seq_len(n_i)) {
      x_j <- X_i[j, ]
      x_out_1 <- matrix(c(x_j, v_i, 1), nrow = 1)
      x_out_0 <- matrix(c(x_j, v_i, 0), nrow = 1)

      # Get predictions from BART 
      b2_1_m <- forest_out$predict_iteration(x_out_1, iter = m)
      b2_0_m <- forest_out$predict_iteration(x_out_0, iter = m)

      # Add random effect
      W_i_m <- w_samples[m, county_id]
      sigma2_m <- sigma2_samps[m]

      # E[T | Z=1] = exp(b2(1) + W + sigma2/2)
      # E[T | Z=0] = exp(b2(0) + W + sigma2/2)
      exp1_vec[j] <- exp(b2_1_m + W_i_m + sigma2_m/2)
      exp0_vec[j] <- exp(b2_0_m + W_i_m + sigma2_m/2)
    }

    County_ATE_draws[m] <- mean(exp1_vec - exp0_vec)
  }

  list(
    County_ATE_mean = mean(County_ATE_draws),
    County_ATE_lower = quantile(County_ATE_draws, 0.025),
    County_ATE_upper = quantile(County_ATE_draws, 0.975),
    County_ATE_draws = County_ATE_draws
  )
}

# Simplified County-RATE estimation
estimate_County_RATE_simple <- function(t_star, X_i, v_i, county_id,
                                        forest_out, w_samples, sigma2_samps,
                                        burn_in, thin, n_mcmc) {
  n_i <- nrow(X_i)
  n_draws <- nrow(w_samples)
  County_RATE_draws <- numeric(n_draws)

  for (m in seq_len(n_draws)) {
    CRATE_vec <- numeric(n_i)

    for (j in seq_len(n_i)) {
      x_j <- X_i[j, ]
      x_out_1 <- matrix(c(x_j, v_i, 1), nrow = 1)
      x_out_0 <- matrix(c(x_j, v_i, 0), nrow = 1)

      # Get predictions from BART
      b2_1_m <- forest_out$predict_iteration(x_out_1, iter = m)
      b2_0_m <- forest_out$predict_iteration(x_out_0, iter = m)

      # Add random effect to get full linear predictor
      W_i_m <- w_samples[m, county_id]
      sigma2_m <- sigma2_samps[m]
      sigma_m <- sqrt(sigma2_m)

      # Full linear predictors on log scale
      b2_1_full <- b2_1_m + W_i_m
      b2_0_full <- b2_0_m + W_i_m

      # E[T | Z] for this individual
      mu1 <- exp(b2_1_full + sigma2_m/2)
      mu0 <- exp(b2_0_full + sigma2_m/2)

      # Standardized log(t*)
      z1 <- (log(t_star) - b2_1_full) / sigma_m
      z0 <- (log(t_star) - b2_0_full) / sigma_m

      # E[min(T, t*) | Z=1] = mu1 * Φ(z1 - sigma) + t* * (1 - Φ(z1))
      term1 <- mu1 * pnorm(z1 - sigma_m)
      term2 <- t_star * (1 - pnorm(z1))

      # E[min(T, t*) | Z=0] = mu0 * Φ(z0 - sigma) + t* * (1 - Φ(z0))
      term3 <- mu0 * pnorm(z0 - sigma_m)
      term4 <- t_star * (1 - pnorm(z0))

      CRATE_vec[j] <- (term1 + term2) - (term3 + term4)
    }

    County_RATE_draws[m] <- mean(CRATE_vec)
  }

  list(
    County_RATE_mean = mean(County_RATE_draws),
    County_RATE_lower = quantile(County_RATE_draws, 0.025),
    County_RATE_upper = quantile(County_RATE_draws, 0.975),
    County_RATE_draws = County_RATE_draws
  )
}




# Modified estimate_CRATE function
estimate_CRATE <- function(t_star, x, county_id,
                           e_hat,  
                           forest_out,
                           w_samples, sigma2_samps,
                           burn_in, thin, n_mcmc) {

  sample_iters_out <- seq(from = burn_in + 1, to = n_mcmc, by = thin)
  Mprime <- length(sample_iters_out)
  CRATE_draws <- numeric(Mprime)

  for (m in seq_len(Mprime)) {
    # Use the provided e_hat directly
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

# Modified estimate_County_RATE function
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
    e_hat_j <- e_hat_vec[j]  # Use the j-th element of the PS vector

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



