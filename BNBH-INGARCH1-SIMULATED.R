# ==============================================================================
#  SIMULATION STUDY: Beta-Negative Binomial Hurdle-INGARCHX(1,1) Model
#  via Two-Phase Adaptive MCMC (Chen & So, 2006)
# ==============================================================================

# ------------------------------
# Load Packages and Set Seed
# ------------------------------
library(mvtnorm)
library(MASS)
library(coda)
library(foreach)       # added for parallel replication
library(doParallel)    # added for parallel replication
set.seed(7)

# ==============================================================================
# SECTION 1: BNB Distribution Utilities
# ==============================================================================

log_dbnb <- function(y, r, gamma_t, phi) {
  lgamma(y + r) - lgamma(y + 1) - lgamma(r) +
    lbeta(phi + r, gamma_t + y) - lbeta(phi, gamma_t)
}

log_fbnb_zero <- function(r, gamma_t, phi) {
  lbeta(phi + r, gamma_t) - lbeta(phi, gamma_t)
}

rbnb <- function(n, r, gamma_t, phi) {
  p <- rbeta(n, phi, gamma_t)
  y <- rnbinom(n, size = r, prob = p)
  return(y)
}

rbnb_truncated <- function(n, r, gamma_t, phi) {
  result <- integer(n)
  for (i in 1:n) {
    repeat {
      val <- rbnb(1, r, gamma_t, phi)
      if (val > 0) {
        result[i] <- val
        break
      }
    }
  }
  return(result)
}

# ==============================================================================
# SECTION 2: Alpha Constraint
# ==============================================================================
check_alpha_constraints <- function(alpha_prop) {
  if (length(alpha_prop) == 3) {
    return(abs(alpha_prop[2]) < 1 && abs(alpha_prop[3]) < 1 &&
             abs(alpha_prop[2] + alpha_prop[3]) < 1)
  } else {
    return(abs(alpha_prop[1]) < 1 && abs(alpha_prop[2]) < 1 &&
             abs(alpha_prop[1] + alpha_prop[2]) < 1)
  }
}

# ==============================================================================
# SECTION 3: Beta Constraint
# ==============================================================================
check_beta_constraints <- function(beta_prop) {
  return(beta_prop[1] > 0 &&
           beta_prop[2] > 0 &&
           beta_prop[3] >= 0 &&
           (beta_prop[2] + beta_prop[3]) < 1)
}

# ==============================================================================
# SECTION 4: Generalized Data Generation Function
# ==============================================================================
BNBH_data <- function(n, alpha, beta, omega, r, phi, X_list, b, b0,
                      pi_init, lambda_init, k) {
  y        <- rep(0, n)
  pi_t     <- rep(pi_init, n)
  lambda_t <- rep(lambda_init, n)
  gamma_1  <- (phi - 1) / r * lambda_init
  if (runif(1) < pi_init) {
    y[1] <- 0
  } else {
    y[1] <- rbnb_truncated(1, r, gamma_1, phi)
  }
  
  for (t in (b0 + 1):n) {
    logit_pi_t_minus_1 <- log(pi_t[t - 1] / (1 - pi_t[t - 1]))
    exponent_pi <- alpha[1] + alpha[2] * y[t - 1] + alpha[3] * logit_pi_t_minus_1
    pi_t[t] <- exp(exponent_pi) / (1 + exp(exponent_pi))
    lambda_t[t] <- beta[1] + beta[2] * y[t - 1] + beta[3] * lambda_t[t - 1]
    if (k >= 1) {
      for (i in 1:k) {
        lambda_t[t] <- lambda_t[t] + omega[i] * X_list[[i]][1, t - b[i]]
      }
    }
    lambda_t[t] <- max(lambda_t[t], 1e-6)
    gamma_t <- (phi - 1) / r * lambda_t[t]
    if (runif(1) < pi_t[t]) {
      y[t] <- 0
    } else {
      y[t] <- rbnb_truncated(1, r, gamma_t, phi)
    }
  }
  return(list(y = y, lambda = lambda_t, pi = pi_t))
}

# ==============================================================================
# SECTION 5: Generalized Log-Likelihood Function
# ==============================================================================
log_likelihood <- function(Y, alpha, beta, omega, r, phi, X_list, b, b0, k) {
  n        <- length(Y)
  log_like <- 0
  pi_t     <- rep(0.3, n)
  lambda_t <- rep(0.1, n)
  
  for (t in (b0 + 1):n) {
    logit_pi_t_minus_1 <- log(pi_t[t - 1] / (1 - pi_t[t - 1]))
    exponent_pi <- alpha[1] + alpha[2] * Y[t - 1] + alpha[3] * logit_pi_t_minus_1
    pi_t[t] <- exp(exponent_pi) / (1 + exp(exponent_pi))
    lambda_t[t] <- beta[1] + beta[2] * Y[t - 1] + beta[3] * lambda_t[t - 1]
    
    if (k >= 1) {
      for (i in 1:k) {
        lambda_t[t] <- lambda_t[t] + omega[i] * X_list[[i]][1, t - b[i]]
      }
    }
    lambda_t[t] <- max(lambda_t[t], 1e-6)
    gamma_t <- (phi - 1) / r * lambda_t[t]
    
    indicator_y_zero <- ifelse(Y[t] == 0, 1, 0)
    
    log_fbnb_y <- lgamma(Y[t] + r) - lgamma(Y[t] + 1) - lgamma(r) +
      lbeta(phi + r, gamma_t + Y[t]) - lbeta(phi, gamma_t)
    log_fbnb_0 <- lbeta(phi + r, gamma_t) - lbeta(phi, gamma_t)
    fbnb_0     <- exp(log_fbnb_0)
    
    log_like <- log_like +
      indicator_y_zero * log(pi_t[t]) +
      (1 - indicator_y_zero) * (log(1 - pi_t[t]) +
                                  log_fbnb_y -
                                  log(1 - fbnb_0))
  }
  return(log_like)
}

# ==============================================================================
# SECTION 6: Generalized MCMC Function (Single Chain)
# ==============================================================================
run_mcmc <- function(Y, X_list, k, N = 20000, burn_in = 8000, b0_mcmc = 3,
                     prior_hyp = list(c1 = 1, c2 = 1,
                                      a1 = 3, a2 = 1,
                                      d1 = 3, d2 = 1)) {
  c1 <- prior_hyp$c1; c2 <- prior_hyp$c2
  a1 <- prior_hyp$a1; a2 <- prior_hyp$a2
  d1 <- prior_hyp$d1; d2 <- prior_hyp$d2
  
  # Step sizes — case-specific tuning
  step_size_alpha <- c(0.04, 0.04, 0.04)
  
  if (k == 0) {
    # Case 1: unchanged
    step_size_beta  <- c(0.12, 0.12, 0.12)
    step_size_omega <- NULL
    step_size_r     <- 3.0
    step_size_phi   <- 2.0
    
  } else if (k == 1) {
    # Case 2
    step_size_beta  <- c(0.10, 0.10, 0.10)
    step_size_omega <- rep(0.10, k)
    step_size_r     <- 2.00
    step_size_phi   <- 2.50
    
  } else if (k == 2) {
    # Case 3: slight adjustment only
    step_size_beta  <- c(0.10, 0.10, 0.10)
    step_size_omega <- rep(0.30, k)
    step_size_r     <- 2.5
    step_size_phi   <- 1.5
  }
  
  alpha_current <- rep(0.1, 3)
  beta_current  <- c(0.1, 0.1, 0.1)
  omega_current <- if (k > 0) rep(0.1, k) else NULL
  r_current     <- 2
  phi_current   <- 3
  b             <- if (k > 0) rep(1, k) else NULL
  
  n_param         <- 3 + 3 + k + 1 + 1 + k
  samples         <- matrix(0, nrow = N - burn_in, ncol = n_param)
  burn_in_samples <- matrix(0, nrow = burn_in,     ncol = n_param)
  
  accept_count       <- 0
  alpha_accept_count <- 0
  beta_accept_count  <- 0
  omega_accept_count <- 0
  r_accept_count     <- 0
  phi_accept_count   <- 0
  
  pack_params <- function() {
    c(alpha_current, beta_current,
      if (k > 0) omega_current else numeric(0),
      r_current, phi_current,
      if (k > 0) b else numeric(0))
  }
  
  idx_alpha <- 1:3
  idx_beta  <- 4:6
  idx_omega <- if (k > 0) (7):(6 + k)           else integer(0)
  idx_r     <- 6 + k + 1
  idx_phi   <- 6 + k + 2
  idx_b     <- if (k > 0) (6 + k + 3):(6 + 2*k + 2) else integer(0)
  
  # ============================================================
  # PHASE 1: Random Walk Metropolis-Hastings (Burn-in)
  # ============================================================
  for (iter in 1:burn_in) {
    alpha_proposal <- c()
    beta_proposal  <- c()
    omega_proposal <- c()
    
    if (iter %% 1000 == 0) cat("  Burn-in Iteration:", iter, "\n")
    
    repeat {
      alpha_proposal[1] <- alpha_current[1] + rnorm(1, 0, step_size_alpha[1])
      alpha_proposal[2] <- alpha_current[2] + rnorm(1, 0, step_size_alpha[2])
      alpha_proposal[3] <- alpha_current[3] + rnorm(1, 0, step_size_alpha[3])
      if (check_alpha_constraints(alpha_proposal)) break
    }
    log_accept_ratio <- log_likelihood(Y, alpha_proposal, beta_current, omega_current,
                                       r_current, phi_current, X_list, b, b0_mcmc, k) -
      log_likelihood(Y, alpha_current, beta_current, omega_current,
                     r_current, phi_current, X_list, b, b0_mcmc, k)
    accept_prob <- min(1, exp(log_accept_ratio)); U <- runif(1)
    if (!is.na(log_accept_ratio) && !is.nan(log_accept_ratio)) {
      if (U < accept_prob) { alpha_current <- alpha_proposal; alpha_accept_count <- alpha_accept_count + 1 }
    }
    
    repeat {
      beta_proposal[1] <- beta_current[1] + rnorm(1, 0, step_size_beta[1])
      beta_proposal[2] <- beta_current[2] + rnorm(1, 0, step_size_beta[2])
      beta_proposal[3] <- beta_current[3] + rnorm(1, 0, step_size_beta[3])
      if (check_beta_constraints(beta_proposal)) break
    }
    log_accept_ratio <- log_likelihood(Y, alpha_current, beta_proposal, omega_current,
                                       r_current, phi_current, X_list, b, b0_mcmc, k) -
      log_likelihood(Y, alpha_current, beta_current, omega_current,
                     r_current, phi_current, X_list, b, b0_mcmc, k)
    accept_prob <- min(1, exp(log_accept_ratio)); U <- runif(1)
    if (!is.na(log_accept_ratio) && !is.nan(log_accept_ratio)) {
      if (U < accept_prob) { beta_current <- beta_proposal; beta_accept_count <- beta_accept_count + 1 }
    }
    
    if (k > 0) {
      repeat {
        for (i in 1:k) omega_proposal[i] <- omega_current[i] + rnorm(1, 0, step_size_omega[i])
        if (all(omega_proposal > 0)) break
      }
      omega_prop_lposterior <- log_likelihood(Y, alpha_current, beta_current, omega_proposal,
                                              r_current, phi_current, X_list, b, b0_mcmc, k) +
        sum(dgamma(omega_proposal, shape = c1, rate = c2, log = TRUE))
      omega_curr_lposterior <- log_likelihood(Y, alpha_current, beta_current, omega_current,
                                              r_current, phi_current, X_list, b, b0_mcmc, k) +
        sum(dgamma(omega_current, shape = c1, rate = c2, log = TRUE))
      log_accept_ratio <- omega_prop_lposterior - omega_curr_lposterior
      accept_prob <- min(1, exp(log_accept_ratio)); U <- runif(1)
      if (!is.na(log_accept_ratio) && !is.nan(log_accept_ratio)) {
        if (U < accept_prob) { omega_current <- omega_proposal; omega_accept_count <- omega_accept_count + 1 }
      }
    }
    
    repeat {
      r_proposal <- r_current[1] + rnorm(1, 0, step_size_r)
      if (r_proposal > 0) break
    }
    r_prop_lposterior <- log_likelihood(Y, alpha_current, beta_current, omega_current,
                                        r_proposal, phi_current, X_list, b, b0_mcmc, k) +
      sum(dgamma(r_proposal, shape = a1, rate = a2, log = TRUE))
    r_curr_lposterior <- log_likelihood(Y, alpha_current, beta_current, omega_current,
                                        r_current, phi_current, X_list, b, b0_mcmc, k) +
      sum(dgamma(r_current, shape = a1, rate = a2, log = TRUE))
    log_accept_ratio <- r_prop_lposterior - r_curr_lposterior
    accept_prob <- min(1, exp(log_accept_ratio)); U <- runif(1)
    if (!is.na(log_accept_ratio) && !is.nan(log_accept_ratio)) {
      if (U < accept_prob) { r_current <- r_proposal; r_accept_count <- r_accept_count + 1 }
    }
    
    repeat {
      phi_proposal <- phi_current + rnorm(1, 0, step_size_phi)
      if (phi_proposal > 2) break
    }
    phi_prop_lposterior <- log_likelihood(Y, alpha_current, beta_current, omega_current,
                                          r_current, phi_proposal, X_list, b, b0_mcmc, k) +
      (d1 - 1) * log(phi_proposal) - d2 * phi_proposal
    phi_curr_lposterior <- log_likelihood(Y, alpha_current, beta_current, omega_current,
                                          r_current, phi_current, X_list, b, b0_mcmc, k) +
      (d1 - 1) * log(phi_current) - d2 * phi_current
    log_accept_ratio <- phi_prop_lposterior - phi_curr_lposterior
    accept_prob <- min(1, exp(log_accept_ratio)); U <- runif(1)
    if (!is.na(log_accept_ratio) && !is.nan(log_accept_ratio)) {
      if (U < accept_prob) { phi_current <- phi_proposal; phi_accept_count <- phi_accept_count + 1 }
    }
    
    if (k > 0) {
      for (bi_idx in 1:k) {
        lik_b <- NULL; lik_bm <- NULL; sum_lik <- 0; prob <- NULL
        b_temp    <- b
        b[bi_idx] <- rep(0, 1)
        lik_b <- sapply(1:b0_mcmc, function(j) {
          b_candidate         <- b_temp
          b_candidate[bi_idx] <- j
          log_likelihood(Y, alpha_current, beta_current, omega_current,
                         r_current, phi_current, X_list, b_candidate, b0_mcmc, k)
        })
        max_lik  <- max(lik_b)
        lik_bm   <- exp(lik_b - max_lik)
        prob     <- lik_bm / sum(lik_bm)
        cum_prob <- cumsum(prob)
        U <- runif(1)
        if (!is.na(cum_prob[1]) && !is.nan(cum_prob[1])) {
          if (U < cum_prob[1]) {
            b[bi_idx] <- 1
          } else {
            B = FALSE; I = 1
            while (B == FALSE) {
              int <- c(cum_prob[I], cum_prob[I + 1])
              if ((U > int[1]) & (U < int[2])) { b[bi_idx] <- I + 1; B = TRUE } else { I = I + 1 }
            }
          }
        }
      }
    }
    burn_in_samples[iter, ] <- pack_params()
  }
  
  cat("\n  Phase 1 Acceptance Rates:\n")
  cat("    alpha:", alpha_accept_count / burn_in, "\n")
  cat("    beta: ", beta_accept_count  / burn_in, "\n")
  if (k > 0) cat("    omega:", omega_accept_count / burn_in, "\n")
  cat("    r:    ", r_accept_count   / burn_in, "\n")
  cat("    phi:  ", phi_accept_count / burn_in, "\n")
  
  mu_alpha  <- colMeans(burn_in_samples[-1:-1000, idx_alpha])
  cov_alpha <- cov(burn_in_samples[-1:-1000, idx_alpha])
  mu_beta   <- colMeans(burn_in_samples[-1:-1000, idx_beta])
  cov_beta  <- cov(burn_in_samples[-1:-1000, idx_beta])
  
  if (k > 0) {
    if (k == 1) {
      mu_omega  <- mean(burn_in_samples[-1:-1000, idx_omega])
      cov_omega <- var(burn_in_samples[-1:-1000,  idx_omega])
    } else {
      mu_omega  <- colMeans(burn_in_samples[-1:-1000, idx_omega])
      cov_omega <- cov(burn_in_samples[-1:-1000,      idx_omega])
    }
  }
  
  mu_r  <- mean(burn_in_samples[-1:-1000, idx_r])
  cov_r <- var(burn_in_samples[-1:-1000,  idx_r])
  mu_phi  <- mean(burn_in_samples[-1:-1000, idx_phi])
  cov_phi <- var(burn_in_samples[-1:-1000,  idx_phi])
  
  log_gaussian_current  <- 0
  log_gaussian_proposal <- 0
  
  alpha_accept_count <- 0; beta_accept_count  <- 0
  omega_accept_count <- 0; r_accept_count     <- 0; phi_accept_count <- 0
  
  # ============================================================
  # PHASE 2: Independent Kernel Metropolis-Hastings
  # ============================================================
  for (iter in 1:(N - burn_in)) {
    if (iter %% 1000 == 0) cat("  Independent-Kernel Iteration:", iter, "\n")
    
    repeat {
      alpha_proposal <- as.numeric(mvrnorm(1, mu = mu_alpha, Sigma = cov_alpha))
      if (check_alpha_constraints(alpha_proposal)) break
    }
    log_gaussian_current  <- mvtnorm::dmvnorm(alpha_current,  mean = mu_alpha, sigma = cov_alpha, log = TRUE)
    log_gaussian_proposal <- mvtnorm::dmvnorm(alpha_proposal, mean = mu_alpha, sigma = cov_alpha, log = TRUE)
    log_accept_ratio <- log_likelihood(Y, alpha_proposal, beta_current, omega_current,
                                       r_current, phi_current, X_list, b, b0_mcmc, k) +
      log_gaussian_current -
      log_likelihood(Y, alpha_current, beta_current, omega_current,
                     r_current, phi_current, X_list, b, b0_mcmc, k) -
      log_gaussian_proposal
    accept_prob <- min(1, exp(log_accept_ratio)); U <- runif(1)
    if (!is.na(log_accept_ratio) && !is.nan(log_accept_ratio)) {
      if (U < accept_prob) { alpha_current <- alpha_proposal; alpha_accept_count <- alpha_accept_count + 1 }
    }
    
    repeat {
      beta_proposal <- as.numeric(mvrnorm(1, mu = mu_beta, Sigma = cov_beta))
      if (check_beta_constraints(beta_proposal)) break
    }
    log_gaussian_current  <- mvtnorm::dmvnorm(beta_current,  mean = mu_beta, sigma = cov_beta, log = TRUE)
    log_gaussian_proposal <- mvtnorm::dmvnorm(beta_proposal, mean = mu_beta, sigma = cov_beta, log = TRUE)
    log_accept_ratio <- log_likelihood(Y, alpha_current, beta_proposal, omega_current,
                                       r_current, phi_current, X_list, b, b0_mcmc, k) +
      log_gaussian_current -
      log_likelihood(Y, alpha_current, beta_current, omega_current,
                     r_current, phi_current, X_list, b, b0_mcmc, k) -
      log_gaussian_proposal
    accept_prob <- min(1, exp(log_accept_ratio)); U <- runif(1)
    if (!is.na(log_accept_ratio) && !is.nan(log_accept_ratio)) {
      if (U < accept_prob) { beta_current <- beta_proposal; beta_accept_count <- beta_accept_count + 1 }
    }
    
    if (k > 0) {
      if (k == 1) {
        repeat { omega_proposal <- rnorm(1, mu_omega, sqrt(cov_omega)); if (omega_proposal > 0) break }
        log_gaussian_current  <- dnorm(omega_current,  mu_omega, sqrt(cov_omega), log = TRUE)
        log_gaussian_proposal <- dnorm(omega_proposal, mu_omega, sqrt(cov_omega), log = TRUE)
      } else {
        repeat { omega_proposal <- as.numeric(mvrnorm(1, mu = mu_omega, Sigma = cov_omega)); if (all(omega_proposal > 0)) break }
        log_gaussian_current  <- dmvnorm(omega_current,  mean = mu_omega, sigma = cov_omega, log = TRUE)
        log_gaussian_proposal <- dmvnorm(omega_proposal, mean = mu_omega, sigma = cov_omega, log = TRUE)
      }
      omega_prop_lposterior <- log_likelihood(Y, alpha_current, beta_current, omega_proposal,
                                              r_current, phi_current, X_list, b, b0_mcmc, k) +
        sum(dgamma(omega_proposal, shape = c1, rate = c2, log = TRUE))
      omega_curr_lposterior <- log_likelihood(Y, alpha_current, beta_current, omega_current,
                                              r_current, phi_current, X_list, b, b0_mcmc, k) +
        sum(dgamma(omega_current, shape = c1, rate = c2, log = TRUE))
      log_accept_ratio <- omega_prop_lposterior + log_gaussian_current -
        omega_curr_lposterior - log_gaussian_proposal
      accept_prob <- min(1, exp(log_accept_ratio)); U <- runif(1)
      if (!is.na(log_accept_ratio) && !is.nan(log_accept_ratio)) {
        if (U < accept_prob) { omega_current <- omega_proposal; omega_accept_count <- omega_accept_count + 1 }
      }
    }
    
    repeat { r_proposal <- rnorm(1, mu_r, sqrt(cov_r)); if (r_proposal > 0) break }
    log_gaussian_current  <- dnorm(r_current,  mu_r, sqrt(cov_r), log = TRUE)
    log_gaussian_proposal <- dnorm(r_proposal, mu_r, sqrt(cov_r), log = TRUE)
    r_prop_lposterior <- log_likelihood(Y, alpha_current, beta_current, omega_current,
                                        r_proposal, phi_current, X_list, b, b0_mcmc, k) +
      sum(dgamma(r_proposal, shape = a1, rate = a2, log = TRUE))
    r_curr_lposterior <- log_likelihood(Y, alpha_current, beta_current, omega_current,
                                        r_current, phi_current, X_list, b, b0_mcmc, k) +
      sum(dgamma(r_current, shape = a1, rate = a2, log = TRUE))
    log_accept_ratio <- r_prop_lposterior + log_gaussian_current - r_curr_lposterior - log_gaussian_proposal
    accept_prob <- min(1, exp(log_accept_ratio)); U <- runif(1)
    if (!is.na(log_accept_ratio) && !is.nan(log_accept_ratio)) {
      if (U < accept_prob) { r_current <- r_proposal; r_accept_count <- r_accept_count + 1 }
    }
    
    repeat { phi_proposal <- rnorm(1, mu_phi, sqrt(cov_phi)); if (phi_proposal > 2) break }
    log_gaussian_current  <- dnorm(phi_current,  mu_phi, sqrt(cov_phi), log = TRUE)
    log_gaussian_proposal <- dnorm(phi_proposal, mu_phi, sqrt(cov_phi), log = TRUE)
    phi_prop_lposterior <- log_likelihood(Y, alpha_current, beta_current, omega_current,
                                          r_current, phi_proposal, X_list, b, b0_mcmc, k) +
      (d1 - 1) * log(phi_proposal) - d2 * phi_proposal
    phi_curr_lposterior <- log_likelihood(Y, alpha_current, beta_current, omega_current,
                                          r_current, phi_current, X_list, b, b0_mcmc, k) +
      (d1 - 1) * log(phi_current) - d2 * phi_current
    log_accept_ratio <- phi_prop_lposterior + log_gaussian_current - phi_curr_lposterior - log_gaussian_proposal
    accept_prob <- min(1, exp(log_accept_ratio)); U <- runif(1)
    if (!is.na(log_accept_ratio) && !is.nan(log_accept_ratio)) {
      if (U < accept_prob) { phi_current <- phi_proposal; phi_accept_count <- phi_accept_count + 1 }
    }
    
    if (k > 0) {
      for (bi_idx in 1:k) {
        lik_b <- NULL; lik_bm <- NULL; sum_lik <- 0; prob <- NULL
        b_temp    <- b
        b[bi_idx] <- rep(0, 1)
        lik_b <- sapply(1:b0_mcmc, function(j) {
          b_candidate         <- b_temp
          b_candidate[bi_idx] <- j
          log_likelihood(Y, alpha_current, beta_current, omega_current,
                         r_current, phi_current, X_list, b_candidate, b0_mcmc, k)
        })
        max_lik  <- max(lik_b)
        lik_bm   <- exp(lik_b - max_lik)
        prob     <- lik_bm / sum(lik_bm)
        cum_prob <- cumsum(prob)
        U <- runif(1)
        if (!is.na(cum_prob[1]) && !is.nan(cum_prob[1])) {
          if (U < cum_prob[1]) {
            b[bi_idx] <- 1
          } else {
            B = FALSE; I = 1
            while (B == FALSE) {
              int <- c(cum_prob[I], cum_prob[I + 1])
              if ((U > int[1]) & (U < int[2])) { b[bi_idx] <- I + 1; B = TRUE } else { I = I + 1 }
            }
          }
        }
      }
    }
    samples[iter, ] <- pack_params()
  }
  
  n_phase2     <- N - burn_in
  n_blocks     <- if (k > 0) 5 else 4
  total_accept <- alpha_accept_count + beta_accept_count + r_accept_count + phi_accept_count
  if (k > 0) total_accept <- total_accept + omega_accept_count
  accept_rate  <- (total_accept / n_blocks) / n_phase2
  
  cat("\n  Phase 2 Acceptance Rates:\n")
  cat("    alpha:", alpha_accept_count / n_phase2, "\n")
  cat("    beta: ", beta_accept_count  / n_phase2, "\n")
  if (k > 0) cat("    omega:", omega_accept_count / n_phase2, "\n")
  cat("    r:    ", r_accept_count   / n_phase2, "\n")
  cat("    phi:  ", phi_accept_count / n_phase2, "\n")
  cat("    Overall:", accept_rate, "\n")
  
  return(list(samples = samples, n_param = n_param, k = k,
              idx_alpha = idx_alpha, idx_beta = idx_beta,
              idx_omega = idx_omega, idx_r = idx_r, idx_phi = idx_phi,
              idx_b = idx_b,
              accept_alpha = alpha_accept_count / n_phase2,
              accept_beta  = beta_accept_count  / n_phase2,
              accept_omega = if (k > 0) omega_accept_count / n_phase2 else NULL,
              accept_r     = r_accept_count   / n_phase2,
              accept_phi   = phi_accept_count / n_phase2))
}

# ==============================================================================
# SECTION 7: Generalized Plotting and Summary Functions
# ==============================================================================
plot_diagnostics <- function(result, true_vals, param_names, case_label) {
  samples <- result$samples
  n_cont  <- length(true_vals)
  for (i in 1:n_cont) {
    par(mfrow = c(1, 1))
    plot(samples[, i], type = "l",
         main = paste0(case_label, " - Traceplot: ", param_names[i]),
         ylab = "Value", xlab = "Iteration")
    abline(h = true_vals[i], col = "red", lwd = 2)
    par(mfrow = c(1, 1))
    acf(samples[, i], main = paste0(case_label, " - ACF: ", param_names[i]))
  }
}

summarize_results <- function(result, true_vals_all, param_names_all, k) {
  samples <- result$samples
  n_cont  <- length(true_vals_all) - k
  mean_cont   <- round(apply(samples[, 1:n_cont, drop = FALSE], 2, mean),   4)
  median_cont <- round(apply(samples[, 1:n_cont, drop = FALSE], 2, median), 4)
  std_cont    <- round(apply(samples[, 1:n_cont, drop = FALSE], 2, sd),     4)
  p025_cont   <- round(apply(samples[, 1:n_cont, drop = FALSE], 2, quantile, probs = 0.025), 4)
  p975_cont   <- round(apply(samples[, 1:n_cont, drop = FALSE], 2, quantile, probs = 0.975), 4)
  mode_fn <- function(x) { ux <- unique(x); ux[which.max(tabulate(match(x, ux)))] }
  if (k > 0) {
    mode_b <- sapply(1:k, function(i) mode_fn(samples[, n_cont + i]))
    summaries <- data.frame(
      Parameter = param_names_all, True = true_vals_all,
      Mean = c(mean_cont, mode_b), Median = c(median_cont, rep("", k)),
      Std = c(std_cont, rep("", k)), P0.025 = c(p025_cont, rep("", k)),
      P0.975 = c(p975_cont, rep("", k)), stringsAsFactors = FALSE)
  } else {
    summaries <- data.frame(
      Parameter = param_names_all, True = true_vals_all,
      Mean = mean_cont, Median = median_cont, Std = std_cont,
      P0.025 = p025_cont, P0.975 = p975_cont, stringsAsFactors = FALSE)
  }
  return(summaries)
}

# ==============================================================================
# SECTION 8: Common Settings
# ==============================================================================
n           <- 1000
alpha_true  <- c(0.30, -0.39, 0.50)
beta_true   <- c(0.50, 0.30, 0.15)
r_true      <- 3
phi_true    <- 4
omega_true_c2 <- 0.12
b_true_c2     <- 1
omega_true_c3 <- c(0.12, 0.13)
b_true_c3     <- c(1, 1)
pi_init     <- 0.3
lambda_init <- 0.1

X1 <- matrix(rnorm(n, 0, 1), 1, ncol = n)
X2 <- matrix(rnorm(n, 0, 1), 1, ncol = n)  # keeping rnorm for now

# ==============================================================================
# SECTION 9: REPLICATION BLOCK — CASE 2 AND CASE 3
# ==============================================================================
prime_seeds <- c(2, 3, 5, 7, 11, 13, 17, 19, 23, 29,
                 31, 37, 41, 43, 47, 53, 59, 61, 67, 71,
                 73, 79, 83, 89, 97, 101, 103, 107, 109, 113,
                 127, 131, 137, 139, 149, 151, 157, 163, 167, 173,
                 179, 181, 191, 193, 197, 199, 211, 223, 227, 229,
                 233, 239, 241, 251, 257, 263, 269, 271, 277, 281,
                 283, 293, 307, 311, 313, 317, 331, 337, 347, 349,
                 353, 359, 367, 373, 379, 383, 389, 397, 401, 409,
                 419, 421, 431, 433, 439, 443, 449, 457, 461, 463,
                 467, 479, 487, 491, 499, 503, 509, 521, 523, 541)

M         <- 100
n_cont_c2 <- 9
n_cont_c3 <- 10

mode_fn <- function(x) {
  ux <- unique(x)
  ux[which.max(tabulate(match(x, ux)))]
}

param_names_c2 <- c("alpha1", "alpha2", "alpha3",
                    "beta1",  "beta2",  "beta3",
                    "omega1", "r",      "phi")

param_names_c3 <- c("alpha1", "alpha2", "alpha3",
                    "beta1",  "beta2",  "beta3",
                    "omega1", "omega2", "r",      "phi")

col_labels_c2 <- c(
  paste0("mean_",   param_names_c2),
  paste0("median_", param_names_c2),
  paste0("sd_",     param_names_c2),
  paste0("p025_",   param_names_c2),
  paste0("p975_",   param_names_c2),
  "mode_b1",
  "accept_alpha", "accept_beta", "accept_omega", "accept_r", "accept_phi"
)

col_labels_c3 <- c(
  paste0("mean_",   param_names_c3),
  paste0("median_", param_names_c3),
  paste0("sd_",     param_names_c3),
  paste0("p025_",   param_names_c3),
  paste0("p975_",   param_names_c3),
  "mode_b1", "mode_b2",
  "accept_alpha", "accept_beta", "accept_omega", "accept_r", "accept_phi"
)

checkpoint_file_c2 <- "checkpoint_case2.csv"
checkpoint_file_c3 <- "checkpoint_case3.csv"

# ==============================================================================
# --- Case 2 Replication ---
# ==============================================================================
if (file.exists(checkpoint_file_c2)) {
  post_means_c2 <- as.matrix(read.csv(checkpoint_file_c2, check.names = FALSE, row.names = NULL)[, -1])
  start_m_c2    <- nrow(post_means_c2) + 1
  cat("Case 2: Resuming from replication", start_m_c2, "\n")
} else {
  post_means_c2 <- matrix(NA, nrow = 0, ncol = n_cont_c2 * 5 + 1 + 5)
  start_m_c2    <- 1
}


# ==============================================================================
# PARALLEL WORKER FUNCTION — Case 2
# ==============================================================================
run_one_rep_c2 <- function(m) {
  set.seed(prime_seeds[m])
  
  X1_rep <- matrix(rnorm(n, 0, 1), 1, ncol = n)
  
  sim_rep <- BNBH_data(n = n, alpha = alpha_true, beta = beta_true,
                       omega = omega_true_c2, r = r_true, phi = phi_true,
                       X_list = list(X1_rep), b = b_true_c2, b0 = 1,
                       pi_init = pi_init, lambda_init = lambda_init, k = 1)
  Y_rep <- sim_rep$y
  res   <- run_mcmc(Y = Y_rep, X_list = list(X1_rep), k = 1)
  
  new_row <- c(
    apply(res$samples[, 1:n_cont_c2], 2, mean),
    apply(res$samples[, 1:n_cont_c2], 2, median),
    apply(res$samples[, 1:n_cont_c2], 2, sd),
    apply(res$samples[, 1:n_cont_c2], 2, quantile, probs = 0.025),
    apply(res$samples[, 1:n_cont_c2], 2, quantile, probs = 0.975),
    mode_fn(res$samples[, res$idx_b[1]]),
    res$accept_alpha, res$accept_beta,
    res$accept_omega, res$accept_r, res$accept_phi
  )
  
  return(new_row)
}

if (start_m_c2 <= M) {
  reps_to_run <- start_m_c2:M
  
  ncores <- max(1L, parallel::detectCores() - 1L)
  cl     <- makeCluster(ncores)
  registerDoParallel(cl)
  
  cat("Case 2: Running replications", start_m_c2, "to", M,
      "on", ncores, "cores...\n")
  t0 <- Sys.time()
  
  results_c2 <- foreach(
    m         = reps_to_run,
    .packages = c("mvtnorm", "MASS"),
    .export   = c("run_one_rep_c2", "run_mcmc", "log_likelihood",
                  "BNBH_data", "check_alpha_constraints",
                  "check_beta_constraints", "mode_fn",
                  "prime_seeds", "n", "n_cont_c2",
                  "alpha_true", "beta_true", "omega_true_c2",
                  "r_true", "phi_true", "b_true_c2",
                  "pi_init", "lambda_init",
                  "param_names_c2", "col_labels_c2",
                  "rbnb", "rbnb_truncated")
  ) %dopar% {
    run_one_rep_c2(m)
  }
  
  stopCluster(cl)
  cat("Case 2 parallel run completed in",
      round(difftime(Sys.time(), t0, units = "mins"), 2), "minutes.\n")
  
  # Collect results into matrix
  new_rows                <- do.call(rbind, results_c2)
  post_means_c2           <- rbind(post_means_c2, new_rows)
  colnames(post_means_c2) <- col_labels_c2
  
  # Save checkpoint once after all replications complete
  write.csv(post_means_c2, checkpoint_file_c2)
  cat("Checkpoint saved: Case 2,", nrow(post_means_c2), "replications.\n")
  
  # Generate plots serially for replication 1
  cat("Generating plots for replication 1...\n")
  set.seed(prime_seeds[1])
  X1_plot  <- matrix(rnorm(n, 0, 1), 1, ncol = n)
  sim_plot <- BNBH_data(n = n, alpha = alpha_true, beta = beta_true,
                        omega = omega_true_c2, r = r_true, phi = phi_true,
                        X_list = list(X1_plot), b = b_true_c2, b0 = 1,
                        pi_init = pi_init, lambda_init = lambda_init, k = 1)
  res_plot     <- run_mcmc(Y = sim_plot$y, X_list = list(X1_plot), k = 1)
  true_vals_c2 <- c(alpha_true, beta_true, omega_true_c2, r_true, phi_true)
  
  pdf("traceplot_case2.pdf")
  par(mfrow = c(3, 3))
  for (i in 1:n_cont_c2) {
    plot(res_plot$samples[, i], type = "l",
         main = param_names_c2[i], ylab = "Value", xlab = "Iteration")
    abline(h = true_vals_c2[i], col = "red", lwd = 2)
  }
  dev.off()
  
  pdf("acf_case2.pdf")
  par(mfrow = c(3, 3))
  for (i in 1:n_cont_c2) {
    acf(res_plot$samples[, i], main = param_names_c2[i])
  }
  dev.off()
  
  cat("Plots saved: traceplot_case2.pdf and acf_case2.pdf\n")
}



# ==============================================================================
# SECTION 10: SIMULATION SUMMARY TABLE
# ==============================================================================
make_summary_table <- function(post_means, param_names, true_vals) {
  n_p      <- length(param_names)
  means_mx <- post_means[, 1:n_p, drop = FALSE]
  col_mean   <- round(colMeans(means_mx), 4)
  col_median <- round(apply(means_mx, 2, median), 4)
  col_sd     <- round(apply(means_mx, 2, sd), 4)
  col_p025   <- round(apply(means_mx, 2, quantile, probs = 0.025), 4)
  col_p975   <- round(apply(means_mx, 2, quantile, probs = 0.975), 4)
  col_bias   <- round(col_mean - true_vals, 4)
  col_rmse   <- round(sqrt(colMeans((means_mx - matrix(true_vals,
                                                       nrow = nrow(means_mx),
                                                       ncol = n_p, byrow = TRUE))^2)), 4)
  data.frame(
    Parameter = param_names, True = true_vals,
    Mean = col_mean, Median = col_median, Std = col_sd,
    P0.025 = col_p025, P0.975 = col_p975,
    Bias = col_bias, RMSE = col_rmse, row.names = NULL)
}

# --- Case 2 ---
true_cont_c2 <- c(alpha_true, beta_true, omega_true_c2, r_true, phi_true)
result_c2    <- make_summary_table(post_means_c2, param_names_c2, true_cont_c2)
mode_b1_c2   <- mode_fn(post_means_c2[, "mode_b1"])
result_c2    <- rbind(result_c2, data.frame(
  Parameter = "b1", True = b_true_c2, Mean = mode_b1_c2,
  Median = "", Std = "", P0.025 = "", P0.975 = "", Bias = "", RMSE = "",
  stringsAsFactors = FALSE))
cat("\n========== Case 2 Summary ==========\n")
print(result_c2)
write.csv(result_c2, "result_case2.csv", row.names = FALSE)

cat("\nCase 2 Average Phase 2 Acceptance Rates across", nrow(post_means_c2), "replications:\n")
cat("  alpha:", round(mean(as.numeric(post_means_c2[, "accept_alpha"]), na.rm = TRUE), 4), "\n")
cat("  beta: ", round(mean(as.numeric(post_means_c2[, "accept_beta"]),  na.rm = TRUE), 4), "\n")
cat("  omega:", round(mean(as.numeric(post_means_c2[, "accept_omega"]), na.rm = TRUE), 4), "\n")
cat("  r:    ", round(mean(as.numeric(post_means_c2[, "accept_r"]),     na.rm = TRUE), 4), "\n")
cat("  phi:  ", round(mean(as.numeric(post_means_c2[, "accept_phi"]),   na.rm = TRUE), 4), "\n")

