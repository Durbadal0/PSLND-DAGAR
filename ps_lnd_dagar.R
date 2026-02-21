# ==============================================================================
# PS-LND-DAGAR: Propensity Score Log-Normal DAGAR Method
# ==============================================================================
#
# Implementation of the method from:
# "Estimation of cluster-specific causal effects on spatially associated
#  survival data using SoftBART"
#
# This file contains:
#   1. scale_to_unit(), scale_new() - Covariate scaling for SoftBART
#   2. create_dagar_precision() - DAGAR precision matrix construction
#   3. fit_ps_sbart() - Stage 1: Propensity score with SoftBART + equicorrelation
#   4. fit_outcome_dagar() - Stage 2: Outcome model with SoftBART + DAGAR
#   5. estimate_cerm(), estimate_county_acerm() - Causal effect estimation
#   6. ps_lnd_dagar() - Main wrapper function
#
# ==============================================================================

library(SoftBart)
library(truncnorm)
library(MASS)

# ==============================================================================
# 1. SCALING FUNCTIONS
# ==============================================================================

#' Scale matrix columns to [0, 1]
#'
#' @param X Numeric matrix
#' @return List with scaled matrix and scaling parameters
scale_to_unit <- function(X) {
  X <- as.matrix(X)
  mins <- apply(X, 2, min)
  maxs <- apply(X, 2, max)
  ranges <- maxs - mins
  ranges[ranges == 0] <- 1  # avoid division by zero for constant columns
  X_scaled <- sweep(sweep(X, 2, mins), 2, ranges, "/")
  list(X_scaled = X_scaled, mins = mins, ranges = ranges)
}

#' Apply saved scaling to new data
#'
#' @param x Numeric matrix or vector (single row)
#' @param params Scaling parameters from scale_to_unit()
#' @return Scaled matrix
scale_new <- function(x, params) {
  x <- as.matrix(x)
  if (ncol(x) != length(params$mins)) {
    stop("Dimension mismatch: x has ", ncol(x), " columns but params has ",
         length(params$mins))
  }
  sweep(sweep(x, 2, params$mins), 2, params$ranges, "/")
}

# ==============================================================================
# 2. DAGAR PRECISION MATRIX CONSTRUCTION
# ==============================================================================

#' Create DAGAR precision matrix components
#'
#' @param adj_matrix K x K symmetric binary adjacency matrix
#' @param ordering Permutation of 1:K specifying DAGAR ordering
#' @param rho Spatial correlation parameter in [0, 1)
#' @return List with Q (precision), B, F_mat, L, N_pi, n_pi
create_dagar_precision <- function(adj_matrix, ordering, rho) {
  K <- nrow(adj_matrix)

  # Validate inputs
  if (!is.matrix(adj_matrix) || nrow(adj_matrix) != ncol(adj_matrix)) {
    stop("adj_matrix must be a square matrix")
  }
  if (!all(adj_matrix == t(adj_matrix))) {
    stop("adj_matrix must be symmetric")
  }
  if (rho < 0 || rho >= 1) {
    stop("rho must be in [0, 1)")
  }
  if (length(ordering) != K || !setequal(ordering, 1:K)) {
    stop("ordering must be a permutation of 1:K")
  }

  # Find directed neighbors for each node
  N_pi <- vector("list", K)
  n_pi <- integer(K)

  for (i in 1:K) {
    node <- ordering[i]  # node at position i in ordering
    neighbors <- unname(which(adj_matrix[node, ] == 1))  # unname to avoid issues with named matrices
    # Directed neighbors: those appearing earlier in the ordering
    if (length(neighbors) == 0) {
      N_pi[[node]] <- integer(0)
      n_pi[node] <- 0
    } else {
      filter_mask <- sapply(neighbors, function(j) which(ordering == j) < i)
      directed_neighbors <- neighbors[filter_mask]
      N_pi[[node]] <- directed_neighbors
      n_pi[node] <- length(directed_neighbors)
    }
  }

  # Construct B matrix
  B <- matrix(0, K, K)
  for (i in 1:K) {
    node <- ordering[i]
    if (n_pi[node] > 0) {
      b_ij <- rho / (1 + (n_pi[node] - 1) * rho^2)
      B[node, N_pi[[node]]] <- b_ij
    }
  }

  # Construct F matrix (diagonal)
  tau_i <- (1 + (n_pi - 1) * rho^2) / (1 - rho^2)
  F_mat <- diag(tau_i)

  # Compute L and Q
  L <- diag(K) - B
  Q <- t(L) %*% F_mat %*% L

  list(Q = Q, B = B, F_mat = F_mat, L = L, N_pi = N_pi, n_pi = n_pi,
       ordering = ordering)
}

# ==============================================================================
# 3. STAGE 1: PROPENSITY SCORE MODEL
# ==============================================================================

#' Fit propensity score model using SoftBART with equicorrelation
#'
#' @param X Covariate matrix (n x p)
#' @param Z Binary treatment indicator
#' @param cluster_id Cluster assignments (1 to K)
#' @param n_mcmc Total MCMC iterations
#' @param burn_in Burn-in iterations
#' @param thin Thinning interval
#' @param num_tree Number of trees for SoftBART
#' @param rho_0_init Initial equicorrelation parameter
#' @param rho_0_prop_sd Proposal SD for rho_0 MH updates
#' @param a_rho_0 Beta prior shape1 for rho_0
#' @param b_rho_0 Beta prior shape2 for rho_0
#' @param verbose Print progress
#' @return List with e_hat, ps_samples, rho_0_samples, forest
fit_ps_sbart <- function(X, Z, cluster_id,
                          n_mcmc = 2000, burn_in = 500, thin = 5,
                          num_tree = 50,
                          rho_0_init = 0.3, rho_0_prop_sd = 0.05,
                          a_rho_0 = 2, b_rho_0 = 2,
                          verbose = TRUE) {

  n <- length(Z)
  unique_clusters <- sort(unique(cluster_id))
  K <- length(unique_clusters)
  cluster_indices <- split(1:n, cluster_id)

  # Standardize covariates (for probit, standardization is fine)
  # Handle constant columns to avoid NaN
  X_scaled <- scale(X)
  # Replace NaN with 0 for constant columns
  X_scaled[is.nan(X_scaled)] <- 0
  # Ensure values are in reasonable range
  X_scaled <- pmax(pmin(X_scaled, 10), -10)

  # Initialize SoftBART for probit
  # For probit: sigma = 1 (fixed), so we use update_sigma = FALSE
  y_init <- qnorm(pmax(0.01, pmin(0.99, mean(Z) + 0.1 * (Z - mean(Z)))))
  hypers <- Hypers(X = X_scaled, Y = y_init, k = 2, num_tree = num_tree,
                   sigma_hat = 1)
  opts <- Opts(update_sigma = FALSE, cache_trees = TRUE)
  forest <- MakeForest(hypers, opts, warn = FALSE)
  b1_current <- as.numeric(forest$do_predict(X_scaled))

  # Initialize latent probit variables
  Z_star <- numeric(n)
  for (j in 1:n) {
    if (Z[j] == 1) {
      Z_star[j] <- rtruncnorm(1, a = 0, b = Inf, mean = b1_current[j], sd = 1)
    } else {
      Z_star[j] <- rtruncnorm(1, a = -Inf, b = 0, mean = b1_current[j], sd = 1)
    }
  }

  # Initialize cluster random effects and equicorrelation
  A_star <- rep(0, K)
  names(A_star) <- as.character(unique_clusters)
  rho_0 <- rho_0_init

  # Storage
  n_save <- floor((n_mcmc - burn_in) / thin)
  ps_samples <- matrix(0, nrow = n_save, ncol = n)
  rho_0_samples <- numeric(n_save)
  sample_idx <- 1

  if (verbose) cat("Stage 1: Propensity Score Estimation\n")

  for (iter in 1:n_mcmc) {

    # Step 1: Update Z* (latent probit variables)
    for (i in seq_along(unique_clusters)) {
      clust_id <- unique_clusters[i]
      idx_i <- cluster_indices[[as.character(clust_id)]]
      mean_vec <- b1_current[idx_i] + rho_0 * A_star[i]
      sd_val <- sqrt(1 - rho_0^2)
      for (j in idx_i) {
        j_local <- which(idx_i == j)
        if (Z[j] == 1) {
          Z_star[j] <- rtruncnorm(1, a = 0, b = Inf,
                                   mean = mean_vec[j_local], sd = sd_val)
        } else {
          Z_star[j] <- rtruncnorm(1, a = -Inf, b = 0,
                                   mean = mean_vec[j_local], sd = sd_val)
        }
      }
    }

    # Step 2: Update A* (cluster random effects)
    for (i in seq_along(unique_clusters)) {
      clust_id <- unique_clusters[i]
      idx_i <- cluster_indices[[as.character(clust_id)]]
      n_i <- length(idx_i)
      sum_resid <- sum(Z_star[idx_i] - b1_current[idx_i])
      denom <- n_i * rho_0^2 + 1 - rho_0^2
      mu_i_star <- rho_0 * sum_resid / denom
      sigma_i_star_sq <- (1 - rho_0^2) / denom
      A_star[i] <- rnorm(1, mean = mu_i_star, sd = sqrt(sigma_i_star_sq))
    }

    # Step 3: Update b1 via SoftBART backfitting
    A_star_expanded <- numeric(n)
    for (i in seq_along(unique_clusters)) {
      idx_i <- cluster_indices[[as.character(unique_clusters[i])]]
      A_star_expanded[idx_i] <- A_star[i]
    }
    y_for_sbart <- Z_star - rho_0 * A_star_expanded
    forest$do_gibbs(X_scaled, y_for_sbart, X_scaled, 1)
    b1_current <- as.numeric(forest$do_predict(X_scaled))

    # Step 4: Update rho_0 via Metropolis-Hastings
    rho_0_prop <- rho_0 + rnorm(1, 0, rho_0_prop_sd)
    if (rho_0_prop > 0 && rho_0_prop < 1) {
      log_lik_curr <- 0
      log_lik_prop <- 0
      for (i in seq_along(unique_clusters)) {
        idx_i <- cluster_indices[[as.character(unique_clusters[i])]]
        resid_i <- Z_star[idx_i] - b1_current[idx_i]
        log_lik_curr <- log_lik_curr +
          sum(dnorm(resid_i, mean = rho_0 * A_star[i],
                    sd = sqrt(1 - rho_0^2), log = TRUE)) +
          dnorm(A_star[i], 0, 1, log = TRUE)
        log_lik_prop <- log_lik_prop +
          sum(dnorm(resid_i, mean = rho_0_prop * A_star[i],
                    sd = sqrt(1 - rho_0_prop^2), log = TRUE)) +
          dnorm(A_star[i], 0, 1, log = TRUE)
      }
      log_prior_curr <- (a_rho_0 - 1) * log(rho_0) +
                        (b_rho_0 - 1) * log(1 - rho_0)
      log_prior_prop <- (a_rho_0 - 1) * log(rho_0_prop) +
                        (b_rho_0 - 1) * log(1 - rho_0_prop)
      log_alpha <- (log_lik_prop + log_prior_prop) -
                   (log_lik_curr + log_prior_curr)
      if (log(runif(1)) < log_alpha) {
        rho_0 <- rho_0_prop
      }
    }

    # Store samples
    if (iter > burn_in && ((iter - burn_in) %% thin == 0)) {
      ps_samples[sample_idx, ] <- pnorm(b1_current)
      rho_0_samples[sample_idx] <- rho_0
      sample_idx <- sample_idx + 1
    }

    if (verbose && iter %% 500 == 0) {
      cat("  Iter", iter, "| rho_0 =", round(rho_0, 3), "\n")
    }
  }

  # Posterior estimates
  e_hat <- colMeans(ps_samples)

  if (verbose) {
    cat("Stage 1 complete. Mean PS:", round(mean(e_hat), 3),
        "| rho_0:", round(mean(rho_0_samples), 3), "\n")
  }

  list(
    e_hat = e_hat,
    ps_samples = ps_samples,
    rho_0_samples = rho_0_samples,
    forest = forest,
    X_scaled = X_scaled
  )
}

# ==============================================================================
# 4. STAGE 2: OUTCOME MODEL WITH DAGAR
# ==============================================================================

#' Fit outcome model with SoftBART + DAGAR spatial random effects
#'
#' @param time Observed survival times
#' @param status Censoring indicator (1 = event, 0 = censored)
#' @param X Covariate matrix
#' @param Z Treatment indicator
#' @param e_hat Estimated propensity scores from Stage 1
#' @param cluster_id Cluster assignments (1 to K)
#' @param adj_matrix K x K spatial adjacency matrix
#' @param ordering DAGAR ordering (default: 1:K)
#' @param n_mcmc Total MCMC iterations
#' @param burn_in Burn-in iterations
#' @param thin Thinning interval
#' @param num_tree Number of trees for SoftBART
#' @param rho_init Initial spatial correlation
#' @param rho_prop_sd Proposal SD for rho MH updates
#' @param a_tau, b_tau Gamma prior parameters for tau_w
#' @param a_sigma, b_sigma Inverse-Gamma prior parameters for sigma^2
#' @param a_rho, b_rho Beta prior parameters for rho
#' @param verbose Print progress
#' @return List with all MCMC samples and forest object
fit_outcome_dagar <- function(time, status, X, Z, e_hat, cluster_id,
                               adj_matrix, ordering = NULL,
                               n_mcmc = 2000, burn_in = 500, thin = 5,
                               num_tree = 50,
                               rho_init = 0.5, rho_prop_sd = 0.02,
                               a_tau = 1, b_tau = 1,
                               a_sigma = 2, b_sigma = 1,
                               a_rho = 2, b_rho = 2,
                               verbose = TRUE) {

  n <- length(time)
  unique_clusters <- sort(unique(cluster_id))
  K <- length(unique_clusters)

  # Create cluster index mapping (ensure 1:K)
  cluster_map <- setNames(1:K, as.character(unique_clusters))
  cluster_idx <- cluster_map[as.character(cluster_id)]
  cluster_indices <- split(1:n, cluster_idx)
  n_i_vec <- sapply(cluster_indices, length)

  # Default ordering
  if (is.null(ordering)) {
    ordering <- 1:K
  }

  # Prepare outcome model covariates
  # Handle zero or very small times: add small constant to avoid log(0) = -Inf
  time_adj <- time
  min_nonzero <- min(time[time > 0])
  zero_idx <- time <= 0
  if (any(zero_idx)) {
    time_adj[zero_idx] <- min_nonzero / 10  # Use 1/10 of smallest positive time
    if (verbose) {
      cat("  Note:", sum(zero_idx), "observations with time <= 0 adjusted\n")
    }
  }
  y_log <- log(time_adj)
  X_out <- cbind(X, Z, e_hat)

  # Scale to [0, 1] for SoftBART
  scaling <- scale_to_unit(X_out)
  X_out_scaled <- scaling$X_scaled

  # Initialize SoftBART
  hypers <- Hypers(X = X_out_scaled, Y = y_log, k = 2, num_tree = num_tree,
                   sigma_hat = sd(y_log), normalize_Y = FALSE)
  opts <- Opts(update_sigma = FALSE, cache_trees = TRUE)
  forest <- MakeForest(hypers, opts, warn = FALSE)
  f_current <- as.numeric(forest$do_predict(X_out_scaled))

  # Initialize parameters
  sigma2 <- var(y_log)
  tau_w <- 1
  rho <- rho_init
  W <- rep(0, K)

  # Initial DAGAR precision
  dagar <- create_dagar_precision(adj_matrix, ordering, rho)
  Q <- dagar$Q

  # Initialize imputed times
  tilde_y <- y_log
  mu_init <- f_current + W[cluster_idx]
  censored_idx <- which(status == 0)
  for (j in censored_idx) {
    tilde_y[j] <- rtruncnorm(1, a = y_log[j], b = Inf,
                              mean = mu_init[j], sd = sqrt(sigma2))
  }

  # Storage
  n_save <- floor((n_mcmc - burn_in) / thin)
  f_samples <- matrix(0, nrow = n_save, ncol = n)
  w_samples <- matrix(0, nrow = n_save, ncol = K)
  sigma2_samples <- numeric(n_save)
  tau_w_samples <- numeric(n_save)
  rho_samples <- numeric(n_save)
  sample_idx <- 1

  if (verbose) cat("Stage 2: Outcome Model with DAGAR\n")

  for (iter in 1:n_mcmc) {

    # Step 1: Impute censored log-times
    mu_current <- f_current + W[cluster_idx]
    for (j in censored_idx) {
      tilde_y[j] <- rtruncnorm(1, a = y_log[j], b = Inf,
                                mean = mu_current[j], sd = sqrt(sigma2))
    }

    # Step 2: Update W (DAGAR random effects)
    b_vec <- numeric(K)
    for (i in 1:K) {
      idx_i <- cluster_indices[[i]]
      b_vec[i] <- sum(tilde_y[idx_i] - f_current[idx_i]) / sigma2
    }
    Lik_prec <- diag(n_i_vec) / sigma2
    Prec_W <- Lik_prec + tau_w * Q
    Cov_W <- solve(Prec_W)
    mu_W <- Cov_W %*% b_vec
    W <- as.numeric(mvrnorm(1, mu = mu_W, Sigma = Cov_W))

    # Step 3: Update tau_w
    shape_tau <- a_tau + K / 2
    rate_tau <- b_tau + 0.5 * as.numeric(t(W) %*% Q %*% W)
    tau_w <- rgamma(1, shape = shape_tau, rate = rate_tau)

    # Step 4: Update rho via Metropolis-Hastings
    rho_prop <- rho + rnorm(1, 0, rho_prop_sd)
    if (rho_prop > 0 && rho_prop < 1) {
      dagar_prop <- create_dagar_precision(adj_matrix, ordering, rho_prop)
      Q_prop <- dagar_prop$Q

      log_det_prop <- as.numeric(determinant(Q_prop, logarithm = TRUE)$modulus)
      log_det_curr <- as.numeric(determinant(Q, logarithm = TRUE)$modulus)

      log_lik_prop <- 0.5 * log_det_prop -
                      0.5 * tau_w * as.numeric(t(W) %*% Q_prop %*% W)
      log_lik_curr <- 0.5 * log_det_curr -
                      0.5 * tau_w * as.numeric(t(W) %*% Q %*% W)

      log_prior_prop <- (a_rho - 1) * log(rho_prop) +
                        (b_rho - 1) * log(1 - rho_prop)
      log_prior_curr <- (a_rho - 1) * log(rho) +
                        (b_rho - 1) * log(1 - rho)

      log_alpha <- (log_lik_prop + log_prior_prop) -
                   (log_lik_curr + log_prior_curr)

      if (log(runif(1)) < log_alpha) {
        rho <- rho_prop
        Q <- Q_prop
      }
    }

    # Step 5: Update b2 via SoftBART backfitting
    R <- tilde_y - W[cluster_idx]
    forest$do_gibbs(X_out_scaled, R, X_out_scaled, 1)
    f_current <- as.numeric(forest$do_predict(X_out_scaled))

    # Step 6: Update sigma^2
    resid_all <- tilde_y - f_current - W[cluster_idx]
    shape_sigma <- a_sigma + n / 2
    rate_sigma <- b_sigma + 0.5 * sum(resid_all^2)
    inv_sigma2 <- rgamma(1, shape = shape_sigma, rate = rate_sigma)
    sigma2 <- 1 / inv_sigma2

    # Store samples
    if (iter > burn_in && ((iter - burn_in) %% thin == 0)) {
      f_samples[sample_idx, ] <- f_current
      w_samples[sample_idx, ] <- W
      sigma2_samples[sample_idx] <- sigma2
      tau_w_samples[sample_idx] <- tau_w
      rho_samples[sample_idx] <- rho
      sample_idx <- sample_idx + 1
    }

    if (verbose && iter %% 500 == 0) {
      cat("  Iter", iter, "| sigma =", round(sqrt(sigma2), 3),
          "| rho =", round(rho, 3), "| tau_w =", round(tau_w, 3), "\n")
    }
  }

  if (verbose) {
    cat("Stage 2 complete. Mean sigma:", round(mean(sqrt(sigma2_samples)), 3),
        "| Mean rho:", round(mean(rho_samples), 3), "\n")
  }

  list(
    forest = forest,
    f_samples = f_samples,
    w_samples = w_samples,
    sigma2_samples = sigma2_samples,
    tau_w_samples = tau_w_samples,
    rho_samples = rho_samples,
    scaling_params = scaling,
    cluster_map = cluster_map,
    n_save = n_save
  )
}

# ==============================================================================
# 5. CERM ESTIMATION
# ==============================================================================

#' Compute CERM (Causal Effect on Restricted Mean) for one subject
#'
#' @param t_star Restriction time
#' @param x Covariate vector for the subject
#' @param county_id Cluster/county ID (original, not mapped)
#' @param e_hat Propensity score for this subject
#' @param fit Output from fit_outcome_dagar()
#' @param orig_to_internal Optional mapping from original county IDs to internal indices
#' @return List with CERM mean, 95% CI, and draws
estimate_cerm <- function(t_star, x, county_id, e_hat, fit,
                          orig_to_internal = NULL) {

  n_draws <- fit$n_save
  CERM_draws <- numeric(n_draws)

  # Map county_id to internal index
  if (!is.null(orig_to_internal)) {
    county_idx <- orig_to_internal[as.character(county_id)]
  } else {
    county_idx <- fit$cluster_map[as.character(county_id)]
  }

  # Construct outcome model inputs for Z=1 and Z=0
  x_out_1 <- matrix(c(x, 1, e_hat), nrow = 1)
  x_out_0 <- matrix(c(x, 0, e_hat), nrow = 1)

  # Scale using saved parameters
  x_out_1_scaled <- scale_new(x_out_1, fit$scaling_params)
  x_out_0_scaled <- scale_new(x_out_0, fit$scaling_params)

  for (m in 1:n_draws) {
    # Get b2 predictions
    b2_1 <- fit$forest$predict_iteration(x_out_1_scaled, iter = m)
    b2_0 <- fit$forest$predict_iteration(x_out_0_scaled, iter = m)

    # Add random effect
    W_m <- fit$w_samples[m, county_idx]
    sigma2_m <- fit$sigma2_samples[m]
    sigma_m <- sqrt(sigma2_m)

    # Full linear predictors
    mu_1 <- b2_1 + W_m
    mu_0 <- b2_0 + W_m

    # E[min(T, t*) | Z=z] for log-normal
    # Formula: exp(mu + sigma^2/2) * Phi((log(t*) - mu - sigma^2)/sigma)
    #        + t* * (1 - Phi((log(t*) - mu)/sigma))
    # Use numerically stable computation to avoid overflow

    log_t_star <- log(t_star)
    z_1 <- (log_t_star - mu_1) / sigma_m
    z_0 <- (log_t_star - mu_0) / sigma_m

    # Cap the exponent to avoid overflow (exp(700) ~ Inf)
    exp_arg_1 <- pmin(mu_1 + sigma2_m / 2, 700)
    exp_arg_0 <- pmin(mu_0 + sigma2_m / 2, 700)

    RM_1 <- exp(exp_arg_1) * pnorm(z_1 - sigma_m) +
            t_star * (1 - pnorm(z_1))
    RM_0 <- exp(exp_arg_0) * pnorm(z_0 - sigma_m) +
            t_star * (1 - pnorm(z_0))

    # Cap at t_star (restricted mean cannot exceed t_star)
    RM_1 <- pmin(RM_1, t_star)
    RM_0 <- pmin(RM_0, t_star)

    CERM_draws[m] <- RM_1 - RM_0
  }

  list(
    mean = mean(CERM_draws),
    lower = quantile(CERM_draws, 0.025),
    upper = quantile(CERM_draws, 0.975),
    draws = CERM_draws
  )
}

#' Compute Average CERM for a county (ACERM)
#'
#' @param t_star Restriction time
#' @param X_county Covariate matrix for subjects in county
#' @param county_id County ID
#' @param e_hat_vec Propensity scores for subjects in county
#' @param fit Output from fit_outcome_dagar()
#' @param orig_to_internal Optional mapping from original county IDs to internal indices
#' @return List with county ACERM mean, 95% CI, and draws
estimate_county_acerm <- function(t_star, X_county, county_id, e_hat_vec, fit,
                                   orig_to_internal = NULL) {

  X_county <- as.matrix(X_county)
  n_county <- nrow(X_county)
  n_draws <- fit$n_save

  # Matrix to hold individual CERM draws
  CERM_matrix <- matrix(0, nrow = n_draws, ncol = n_county)

  for (j in 1:n_county) {
    cerm_j <- estimate_cerm(t_star, X_county[j, ], county_id, e_hat_vec[j], fit,
                             orig_to_internal)
    CERM_matrix[, j] <- cerm_j$draws
  }

  # Average across subjects at each MCMC iteration
  ACERM_draws <- rowMeans(CERM_matrix)

  list(
    mean = mean(ACERM_draws),
    lower = quantile(ACERM_draws, 0.025),
    upper = quantile(ACERM_draws, 0.975),
    draws = ACERM_draws
  )
}

# ==============================================================================
# 6. MAIN WRAPPER FUNCTION
# ==============================================================================

#' PS-LND-DAGAR: Full two-stage causal inference method
#'
#' @param time Observed survival times
#' @param status Censoring indicator (1 = event, 0 = censored)
#' @param X Covariate matrix (n x p)
#' @param Z Binary treatment indicator
#' @param cluster_id Cluster/county assignments
#' @param adj_matrix K x K spatial adjacency matrix
#' @param ordering DAGAR ordering (default: 1:K)
#' @param n_mcmc_ps MCMC iterations for propensity score (Stage 1)
#' @param burn_in_ps Burn-in for Stage 1
#' @param n_mcmc_out MCMC iterations for outcome model (Stage 2)
#' @param burn_in_out Burn-in for Stage 2
#' @param thin Thinning interval
#' @param num_tree_ps Trees for propensity score SoftBART
#' @param num_tree_out Trees for outcome SoftBART
#' @param verbose Print progress
#' @param ... Additional arguments passed to fit_ps_sbart and fit_outcome_dagar
#' @return List with Stage 1 and Stage 2 results
ps_lnd_dagar <- function(time, status, X, Z, cluster_id, adj_matrix,
                          ordering = NULL,
                          n_mcmc_ps = 2000, burn_in_ps = 500,
                          n_mcmc_out = 2000, burn_in_out = 500,
                          thin = 5,
                          num_tree_ps = 50, num_tree_out = 50,
                          verbose = TRUE, ...) {

  # Validate inputs
  n <- length(time)
  if (length(status) != n || length(Z) != n || length(cluster_id) != n) {
    stop("time, status, Z, and cluster_id must have the same length")
  }
  X <- as.matrix(X)
  if (nrow(X) != n) {
    stop("X must have same number of rows as length of time")
  }

  K_adj <- nrow(adj_matrix)
  unique_clusters <- sort(unique(cluster_id))
  K <- length(unique_clusters)

  # Check that all cluster IDs are valid (within 1:K_adj)
  if (any(cluster_id < 1) || any(cluster_id > K_adj)) {
    stop("cluster_id values must be between 1 and ", K_adj)
  }

  # If not all clusters have data, subset the adjacency matrix
  if (K < K_adj) {
    if (verbose) {
      cat("Note: Data has", K, "of", K_adj, "clusters with observations\n")
      cat("      Subsetting adjacency matrix to observed clusters\n")
    }
    adj_matrix <- adj_matrix[unique_clusters, unique_clusters]
    # Remove row/column names to avoid issues with named indices
    rownames(adj_matrix) <- NULL
    colnames(adj_matrix) <- NULL
    # Remap cluster IDs to 1:K
    cluster_id_map <- setNames(1:K, as.character(unique_clusters))
    cluster_id_orig <- cluster_id
    cluster_id <- cluster_id_map[as.character(cluster_id)]
  } else {
    cluster_id_orig <- cluster_id
  }

  # Handle ordering: must be 1:K for the (possibly subsetted) adjacency matrix
  if (is.null(ordering)) {
    ordering <- 1:K
  } else {
    # If ordering was provided for the full adjacency matrix, adjust it
    if (length(ordering) == K_adj && K < K_adj) {
      # Use default ordering 1:K for the subsetted matrix
      ordering <- 1:K
      if (verbose) {
        cat("      Ordering adjusted to 1:", K, " for subsetted matrix\n")
      }
    } else if (length(ordering) != K) {
      stop("ordering must be a permutation of 1:", K)
    }
  }

  if (verbose) {
    cat("==================================================\n")
    cat("PS-LND-DAGAR: Causal Inference with Spatial Effects\n")
    cat("==================================================\n")
    cat("N =", n, "| K =", K, "clusters\n")
    cat("Treated:", sum(Z == 1), "| Control:", sum(Z == 0), "\n")
    cat("Censored:", sum(status == 0), "(", round(100*mean(status==0), 1), "%)\n\n")
  }

  # Stage 1: Propensity Score
  ps_fit <- fit_ps_sbart(X, Z, cluster_id,
                          n_mcmc = n_mcmc_ps, burn_in = burn_in_ps, thin = thin,
                          num_tree = num_tree_ps, verbose = verbose, ...)
  e_hat <- ps_fit$e_hat

  if (verbose) cat("\n")

  # Stage 2: Outcome Model with DAGAR
  out_fit <- fit_outcome_dagar(time, status, X, Z, e_hat, cluster_id,
                                adj_matrix, ordering = ordering,
                                n_mcmc = n_mcmc_out, burn_in = burn_in_out,
                                thin = thin, num_tree = num_tree_out,
                                verbose = verbose, ...)

  if (verbose) {
    cat("\n==================================================\n")
    cat("Estimation complete.\n")
    cat("Use estimate_cerm() or estimate_county_acerm() for causal effects.\n")
    cat("==================================================\n")
  }

  # Create mapping from original county IDs to internal IDs
  if (exists("cluster_id_map")) {
    orig_to_internal <- cluster_id_map
  } else {
    orig_to_internal <- setNames(1:K, as.character(1:K))
  }

  list(
    ps_fit = ps_fit,
    out_fit = out_fit,
    e_hat = e_hat,
    n = n,
    K = K,
    cluster_id = cluster_id,
    cluster_id_orig = cluster_id_orig,
    orig_to_internal = orig_to_internal,
    unique_clusters = unique_clusters
  )
}
