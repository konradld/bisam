estimate_bisam <- function(
    data,
    do_constant = FALSE,
    do_individual_fe = FALSE,
    do_time_fe = FALSE,
    y_index = 3,
    i_index = 1,
    t_index = 2,
    do_center_y = FALSE,
    do_scale_y = FALSE,
    do_center_x = FALSE,
    do_scale_x = FALSE,
    Ndraw = 10^4L,
    Nburn = 2000L,
    beta_prior = "g",  # Options: "g" (Zellner's G-Prior), "f" (Hagan's Fractional Prior)
    step_size_prior = c("mom", "imom"),
    step_incl_prior = c("bern", "beta_bern"),
    # Prior hyperparameters
    beta_variance_scale = 1000, # Variance scaling parameter for beta coefficients
    sigma2_shape = NULL,        # Inverse Gamma shape parameter for sigma^2
    sigma2_rate = NULL,         # Inverse Gamma rate parameter for sigma^2
    sigma2_hyper_p = 0.9,       # Probability that the true sigma^2 is smaller than OLS (if sigma2 params are NULL)
    step_incl_prob = 0.5,       # Bernoulli inclusion probability
    step_incl_alpha = 1,        # Beta-Bernoulli alpha parameter
    step_incl_beta = 1,         # Beta-Bernoulli beta parameter
    step_size_scale = 0.348,    # iMOM/MOM scale parameter relative to sigma^2
    # Additional settings
    cred_int = 0.95,
    do_split_Z = TRUE,
    do_cluster_s2 = FALSE,
    do_check_outlier = FALSE,
    outlier_incl_alpha = 1,     # P(outlier) ~ Beta(outlier_incl_alpha, outlier_incl_beta)
    outlier_incl_beta = 10,     # P(outlier) ~ Beta(outlier_incl_alpha, outlier_incl_beta)
    outlier_scale = 10,
    do_sparse_computation = FALSE,
    do_geweke_test = FALSE
) {
  # ============================================================================
  # GEWEKE DIAGNOSTIC MODE SETUP
  # ============================================================================
  if (do_geweke_test) {
    library(extraDistr)
    beta_variance_scale <- 1.234
    Ndraw <- 100000
    sigma2_shape <- 3
    sigma2_rate <- 1
    beta_prior <- "g"
  }
  
  # ============================================================================
  # REQUIRED PACKAGES
  # ============================================================================
  require(mvtnorm)
  library(Matrix)
  library(matrixStats)
  require(dplyr)
  library(mombf)
  library(glmnet)
  
  # ============================================================================
  # INITIAL SETUP
  # ============================================================================
  cl <- match.call()
  start_time <- Sys.time()
  
  # Rename columns for clarity
  colnames(data)[-c(y_index, i_index, t_index)] <- 
    paste0("beta.", colnames(data[, -c(y_index, i_index, t_index), drop = FALSE]))
  
  # Extract indices
  n_ind <- unique(factor(data[, i_index]))
  t_ind <- unique(data[, t_index])
  
  # Extract response and covariates
  y <- as.matrix(data[, y_index])
  y <- scale(y, center = do_center_y, scale = do_scale_y)
  
  X <- as.matrix(data[, -c(y_index, i_index, t_index)])
  X <- scale(X, center = do_center_x, scale = do_scale_x)
  
  # Dimensions
  n <- length(n_ind)
  t <- length(t_ind)
  N <- n * t
  N_idx <- matrix(1:N, ncol = t, byrow = TRUE)
  
  # ============================================================================
  # BUILD DESIGN MATRIX X
  # ============================================================================
  if (do_sparse_computation) library(Matrix); library(RcppEigen)
  orig_p <- ncol(X) # save original nr. of columns
  
  # Add constant if requested
  if (do_constant) {
    X <- cbind(X, "const" = 1)
  }
  
  # Add individual fixed effects
  if (do_individual_fe) {
    IFE <- kronecker(
      if (do_sparse_computation) Diagonal(n) else diag(n),
      matrix(1, t)
    )
    colnames(IFE) <- paste0("ife.", n_ind)
    if (do_constant) IFE <- IFE[, -1]
    X <- cbind(X, IFE)
  }
  
  # Add time fixed effects
  if (do_time_fe) {
    TFE <- kronecker(
      matrix(1, n),
      if (do_sparse_computation) Diagonal(t) else diag(t)
    )
    colnames(TFE) <- paste0("tfe.", unique(t_ind))
    if (do_constant) TFE <- TFE[, -1]
    X <- cbind(X, TFE)
  }
  
  # Handle collinearity when both FE are used
  if (do_time_fe & do_individual_fe & !do_constant) {
    warning("Both time and unit fixed effects used.\nDropping first indiv. FE to avoid perfect colinearity")
    X <- X[, -(orig_p + 1)]
  }
  
  # Precompute cross-product
  XX <- crossprod(X)
  
  p <- ncol(X)
  
  # ============================================================================
  # BUILD STRUCTURAL BREAK MATRIX Z
  # ============================================================================
  z <- lower.tri(matrix(1, nrow = t, ncol = t), diag = TRUE)[, -c(1:2, t)]
  mode(z) <- "integer"
  Z <- kronecker(if (do_sparse_computation) Diagonal(n) else diag(n), z)
  if (!do_sparse_computation) mode(Z) <- "integer"
  
  # Create grids for labeling
  sis_grid <- expand.grid(t_ind[-c(1:2, t)], n_ind)
  iis_grid <- expand.grid(t_ind, n_ind)
  colnames(Z) <- paste0("sis.", sis_grid[, "Var2"], ".", sis_grid[, "Var1"])
  
  r <- ncol(Z)
  
  # ============================================================================
  # PRIOR SPECIFICATION
  # ============================================================================
  
  # --- Sigma^2 Prior ---
  if (is.null(sigma2_shape) | is.null(sigma2_rate)) {
    print("Sigma^2 prior is inadmissable. Using default spec. based on OLS")
    if (do_cluster_s2) {
      sigma2_shape <- sigma2_rate <- numeric(n)
      for (i in 1:n) {
        n_idx <- N_idx[i,]
        X_tmp <- X[n_idx, , drop = FALSE]
        X_tmp <- X_tmp[, colSums(X_tmp != 0) > 1, drop = FALSE] # >1 to drop TFE
        mod_prior <- lm.fit(as.matrix(X_tmp), y[n_idx])
        s2_OLS <- sum(mod_prior$residuals^2) / mod_prior$df.residual
        s2_pars <- inv_gamma_params(shape = 3, s2_OLS, p = sigma2_hyper_p)
        sigma2_shape[i] <- s2_pars$shape
        sigma2_rate[i] <- s2_pars$rate
      }
    } else {
      mod_prior <- lm.fit(as.matrix(X), y)
      s2_OLS <- sum(mod_prior$residuals^2) / mod_prior$df.residual
      s2_pars <- inv_gamma_params(shape = 3, s2_OLS, p = 0.9)
      sigma2_shape <- s2_pars$shape
      sigma2_rate <- s2_pars$rate
    }
  }
  
  # --- Beta Prior ---
  if (beta_prior == "g" || beta_prior == "hs") {
    beta_mean <- numeric(p)
  }
  if (beta_prior == "f") {
    if (exists("mod_prior") & !do_cluster_s2) { 
      beta_mean <- mod_prior$coefficients
    } else {
      beta_mean <- lm.fit(as.matrix(X), y)$coefficients
    }
  }
  if (beta_prior == "flasso") {
    lasso_model <- glmnet(X, y, alpha = 1, intercept = FALSE)
    best_lambda <- cv.glmnet(X, y, alpha = 1)$lambda.min
    beta_mean <- as.array(coef(lasso_model, s = best_lambda))[-1]
  }
  # beta_var_inv <- XX / beta_variance_scale
  
  # --- Structural Break Prior (MOM/iMOM) ---
  if (step_size_prior == "imom") {
    sis_prior_f <- mombf::imomprior(tau = step_size_scale)
  } else if (step_size_prior == "mom") {
    sis_prior_f <- mombf::momprior(tau = step_size_scale)
  } else {
    stop("sis prior not implemented, choose step_size_prior = c('mom','imom')")
  }
  
  # --- Inclusion Prior ---
  if (step_incl_prior == "bern") {
    incl_prior_f <- mombf::modelbinomprior(step_incl_prob)
    cat("Inclusion prior is Bernoulli(step_incl_prob) - 'step_incl_alpha' and 'step_incl_beta' have no meaning\n")
  } else if (step_incl_prior == "beta_bern") {
    incl_prior_f <- mombf::modelbbprior(step_incl_alpha, step_incl_beta)
    cat("Inclusion prior is Beta-Bernoulli(step_incl_alpha, step_incl_beta) - 'step_incl_prob' has no meaning\n")
  }
  
  # ============================================================================
  # STARTING VALUES
  # ============================================================================
  b_i <- beta_mean
  g_i <- numeric(r)
  g_incl_i <- rep(sd(y) / 100, r)  # Always non-zero for rnlp within Gibbs
  s2_i <- rep(1 / rgamma(1, shape = sigma2_shape, rate = sigma2_rate), N)
  sqrt_s2_i <- sqrt(s2_i)
  
  w_i <- logical(r)
  z_cols <- rep(1:ncol(z), n)
  w_1 <- z_cols * w_i
  
  o_i <- numeric(N)
  pip_i <- numeric(r)
  
  # --- Starting Values for Horseshoe Plug-In ---
  hs_lamb_b <- rep(1, p)
  hs_tau_b <- 1
  hs_nu_b <- rep(1, p)
  hs_xi_b <- 1
  
  # ============================================================================
  # STORAGE ARRAYS
  # ============================================================================
  if (Nburn >= Ndraw) {
    stop('The number of burn-in exceeds number of draws')
  }
  
  Nstore <- Ndraw - Nburn
  
  # Main parameter storage
  b_store <- matrix(NA, nrow = Nstore, ncol = p)
  g_store <- matrix(NA, nrow = Nstore, ncol = r)
  s2_store <- matrix(NA, nrow = Nstore, ncol = ifelse(do_cluster_s2, n, 1))
  
  colnames(b_store) <- colnames(X)
  colnames(g_store) <- colnames(Z)
  colnames(s2_store) <- if (do_cluster_s2) paste0('sigma2_', unique(data[, i_index])) else 'sigma2'
  
  # Selection indicator storage
  w_store <- matrix(NA, nrow = Nstore, ncol = r) # for step selection indicator
  pip_store <- matrix(NA, nrow = Nstore, ncol = r) # for step selection probability
  o_store <- matrix(NA, nrow = Nstore, ncol = N) # for outlier indicator
  
  colnames(w_store) <- colnames(Z)
  colnames(pip_store) <- colnames(Z)
  colnames(o_store) <- paste0("iis.", iis_grid[, "Var2"], ".", iis_grid[, "Var1"])
  
  # Pre-compute
  XX_inv <- solve(XX)
  sqrt_outlier_scale <- sqrt(outlier_scale)
  
  nn <- 2 # number of modelSelection iterations
  
  ngdraw <- 1 # number of rnlp_new iterations saved
  ngburn <- 2 # number of rnlp_new iterations burned
  
  # Initialize objects in loop
  y_aug <- y
  Xb_i <- X %*% b_i
  Zg_i <- Z %*% g_i
  
  if(do_cluster_s2) {
    s2_i_unique <- numeric(n)
    residuals_matrix <- matrix(NA, nrow = t, ncol = n)
  } else {
    s2_i_unique <- numeric(1)
  }
  
  if (do_split_Z) {
    batch <- matrix(1:r, ncol = r / n, byrow = TRUE)
  } else {
    batch <- matrix(1:r, ncol = r)
  }
  
  # ============================================================================
  # GIBBS SAMPLER
  # ============================================================================
  pb <- txtProgressBar(min = 0, max = Ndraw, style = 3)
  
  for (i in (1 - Nburn):Nstore) {
    
    if(do_check_outlier){
      # ========================================================================
      # DRAW p(kappa | omikron, y)
      # ========================================================================
      n_outliers <- sum(o_i)
      k_i <- rbeta(1, outlier_incl_alpha + n_outliers, outlier_incl_beta + N - n_outliers)
      
      # ========================================================================
      # DRAW p(omikron | kappa, beta, gamma, sigma, y)
      # ========================================================================
      residuals <- as.vector(y - Xb_i - Zg_i)

      # Compute log probabilities for outlier models
      log_p0 <- dnorm(residuals, 0, sqrt_s2_i, log = TRUE) + log(1 - k_i)
      log_p1 <- dnorm(residuals, 0, sqrt_outlier_scale * sqrt_s2_i, log = TRUE) + log(k_i)
      
      # Direct calculation of log probability
      log_prob_outlier <- log_p1 - matrixStats::rowLogSumExps(cbind(log_p0, log_p1))
      o_i <- rbinom(N, 1, exp(log_prob_outlier)) # works
      
      # Standardize outliers
      outlier_rescaling <- ifelse(o_i == 0, 1, sqrt_outlier_scale)
      if (any(o_i > 0)) {
        y_aug <- y / outlier_rescaling
      } else {
        y_aug <- y
      }
    }
    
    # ==========================================================================
    # DRAW p(sigma^2 | beta, gamma, y)
    # ==========================================================================
    y_tmp <- y_aug - Zg_i
    residuals <- y_tmp - Xb_i
    
    if (do_cluster_s2) {
      residuals_matrix <- matrix(residuals, nrow = t, ncol = n)
      cN <- sigma2_shape+t/2
      CN <- sigma2_rate + 0.5 * colSums(residuals_matrix^2)
      s2_i_unique <- 1 / rgamma(n, shape = cN, rate = CN) # works
      s2_i <- rep(s2_i_unique, each = t)
    } else {
      cN <- sigma2_shape + n * t / 2
      CN <- sigma2_rate + 0.5 * sum(residuals^2)
      s2_i_unique <- 1 / rgamma(1, shape = cN, rate = CN)
      s2_i[] <- s2_i_unique
    }
    sqrt_s2_i <- sqrt(s2_i)
    
    # ==========================================================================
    # DRAW p(beta | sigma^2, gamma, y)
    # ==========================================================================
    if (beta_prior == "g" || beta_prior == "f"|| beta_prior == "flasso") {
      if (do_cluster_s2) {
        XtSX <- crossprod(X / s2_i, X)
        beta_var_inv <- XtSX / beta_variance_scale
        
        BN <- Matrix::solve(crossprod(X * (1 / s2_i + 1 / (beta_variance_scale * s2_i)), X))
        bN <- BN %*% (crossprod(X / s2_i, y_tmp) + beta_var_inv %*% beta_mean)
      } else {
        BN <- s2_i_unique * beta_variance_scale / (1 + beta_variance_scale) * XX_inv
        bN <- 1 / (1 + beta_variance_scale) * 
          (beta_variance_scale * XX_inv %*% crossprod(X, y_tmp) + beta_mean)
      }
    } else {
      stop("For 'beta' only g-prior and fractional-prior is implemented!")
    }
    
    b_i <- try((bN + t(chol(BN)) %*% rnorm(p)), silent = TRUE)
    if (is(b_i, "try-error")) { b_i <- t(rmvnorm(1,bN,BN)) }
    
    Xb_i <- X %*% b_i
    
    # ==========================================================================
    # DRAW p(omega | beta, sigma^2, y) - Selection Indicators
    # ==========================================================================
    y_tmp <- y_aug - Xb_i
    
    # Standardize for model selection
    y_tmp_sd <- y_tmp / sqrt_s2_i
    s2_i_tmp <- 1
    
    for (j in 1:n) {

      p_idx <- batch[j, ]
      p_idx_rand <- p_idx  # Keep order intact (can be randomized: sample(p_idx))
      
      if (do_split_Z) {
        n_idx <- N_idx[j,]
      } else {
        n_idx <- 1:N
      }
      
      # Set initial parameters
      if (i == (1 - Nburn)) {
        initpar <- "auto"
      } else {
        initpar <- g_i[p_idx_rand]
      }
      
      # Set variance prior
      if (do_cluster_s2) {
        var_prior_f <- igprior(sigma2_shape[j], sigma2_rate[j])
      } else {
        var_prior_f <- igprior(sigma2_shape, sigma2_rate)
      }
      
      # Model selection using mombf
      w_i_mod <- mombf::modelSelection(
        y = y_tmp_sd[n_idx],
        x = Z[n_idx, p_idx_rand],
        groups = 1:length(p_idx),
        nknots = 9,
        center = FALSE,
        scale = FALSE,
        enumerate = FALSE,
        includevars = rep(FALSE, length(p_idx)),
        niter = nn,
        thinning = 1,
        burnin = nn - 1,
        family = "normal",
        priorCoef = sis_prior_f,
        priorDelta = incl_prior_f,
        phi = s2_i_tmp,
        deltaini = w_i[p_idx_rand],
        initSearch = 'none',
        method = 'ALA',
        hess = "asymp",
        initpar = initpar,
        adj.overdisp = 'intercept',
        optimMethod = "auto",
        optim_maxit = 10,
        B = 10^5,
        priorVar = igprior(sigma2_shape[j], sigma2_rate[j]),
        priorSkew = momprior(tau = step_size_scale),
        XtXprecomp = TRUE,
        verbose = FALSE
      )
      
      w_i[p_idx_rand] <- as.logical(w_i_mod$postSample)
      pip_i[p_idx_rand] <- w_i_mod$margpp
      
      # ========================================================================
      # DRAW p(gamma | omega, beta, sigma^2, y) - Break Magnitudes
      # ========================================================================
      if (any(w_i[p_idx_rand])) {
        g_draw <- matrix(0, nrow = ngdraw, ncol = t - 2)
        colsel <- which(w_i_mod$postSample == 1)
        
        g_draw[1, c(colsel, t - 2)] <- rnlp_new(
          y = y_tmp[n_idx],
          x = w_i_mod$xstd[, colsel, drop = FALSE],
          priorCoef = w_i_mod$priors$priorCoef,
          priorGroup = w_i_mod$priors$priorGroup,
          priorVar = w_i_mod$priors$priorVar,
          isgroup = rep(FALSE, length(colsel)),
          niter = ngburn + ngdraw,
          burnin = ngburn,
          knownphi = TRUE,
          phi = s2_i_unique,
          use_thinit = TRUE,
          thinit = g_incl_i[p_idx][colsel]
        )
        g_incl_i[p_idx][w_i[p_idx]] <- g_draw[colsel]
      } else {
        g_draw <- matrix(0, nrow = ngdraw, ncol = t - 2)
      }
      g_i[p_idx_rand] <- g_draw[-(t - 2)]  # Remove phi estimate
    }
    
    Zg_i <- Z[, w_i, drop = FALSE] %*% g_i[w_i]
    
    # ==========================================================================
    # STORE DRAWS (After Burn-in)
    # ==========================================================================
    if (i > 0) {
      b_store[i, ] <- as.matrix(b_i)
      g_store[i, ] <- g_i
      w_store[i, ] <- w_i
      o_store[i, ] <- o_i
      s2_store[i, ] <- s2_i_unique
      pip_store[i,] <- pip_i
    }
    
    setTxtProgressBar(pb, (i + Nburn))
  }
  
  # ============================================================================
  # POST-PROCESSING
  # ============================================================================
  timer <- Sys.time() - start_time
  close(pb)
  cat("Finished after ", format(round(timer, 2)), ".\n", sep = "")
  
  # --- Posterior Means ---
  beta_hat <- colMeans(b_store)
  s2_hat <- colMeans(s2_store)
  gamma_hat <- colMeans(g_store)
  omega_hat <- colMeans(w_store)
  p_outl_hat <- colMeans(o_store)
  
  # --- Credible Intervals ---
  estimates <- data.frame(
    MAP = c(beta_hat, gamma_hat, s2_hat),
    CI_Lower = c(
      apply(b_store, 2, quantile, (1 - cred_int) / 2),
      apply(g_store, 2, quantile, (1 - cred_int) / 2),
      apply(s2_store, 2, quantile, (1 - cred_int) / 2)
    ),
    CI_Upper = c(
      apply(b_store, 2, quantile, 1 - (1 - cred_int) / 2),
      apply(g_store, 2, quantile, 1 - (1 - cred_int) / 2),
      apply(s2_store, 2, quantile, 1 - (1 - cred_int) / 2)
    ),
    PIP = c(
      rep(NA, ncol(b_store)),
      omega_hat,
      rep(NA, ncol(s2_store))
    )
  )
  
  # --- Model Fit Statistics ---
  y_pred <- X %*% beta_hat + Z %*% gamma_hat
  
  # Total R-squared
  SST <- sum((y - mean(y))^2)
  SSR <- sum((y - y_pred)^2)
  R2_total <- (SST - SSR) / SST
  
  # Country-specific R-squared
  R2_country <- numeric(n)
  for (i in 1:n) {
    idx <- ((i - 1) * t + 1):(i * t)
    y_i <- y[idx]
    y_pred_i <- y_pred[idx]
    SST_i <- sum((y_i - mean(y_i))^2)
    SSR_i <- sum((y_i - y_pred_i)^2)
    R2_country[i] <- (SST_i - SSR_i) / SST_i
  }
  
  # ============================================================================
  # OUTPUT
  # ============================================================================
  out <- list(
    "draws" = list(
      "beta" = b_store,
      "sigma2" = s2_store,
      "sis" = g_store,
      "omega" = w_store,
      "iis" = o_store
    ),
    "data" = list(
      "y" = y,
      "X" = X,
      "Z" = Z, 
      "original_data" = data
    ),
    "fitted" = y_pred,
    "meta" = list(
      "timer" = timer, 
      "sample" = c("n" = n, "t" = t, "N" = N, "K" = p + r),
      "MCMC" = c("Ndraw" = Ndraw, "Nburn" = Nburn, "Nstore" = Nstore), 
      "call" = cl
    ),
    "coefs" = list(
      "beta" = beta_hat,
      "sigma2" = s2_hat,
      "sis" = gamma_hat,
      "omega" = omega_hat,
      "iis" = p_outl_hat
    ),
    "coef_list" = estimates,
    "R^2" = list(
      "total_R2" = R2_total,
      "country_R2" = R2_country
    )
  )
  
  class(out) <- "ism"
  
  return(out)
}

#===============================================================================
#
#                  Helper Funcions for inv. Gamma spec
#
#===============================================================================

# Single quantile constraint: P(X <= x) = p
inv_gamma_params <- function(shape = 3, x, p = 0.9) {
  require(invgamma)
  # Solve for scale such that P(X <= x) = p
  result <- uniroot(
    # \propto x^(-shape-1) exp(-rate/x) == 1/gamma(shape,rate)
    function(b) invgamma::pinvgamma(x, shape = shape, scale = b) - p,
    interval = c(1e-6, 1e6)
  )
  
  return(list(shape = shape, scale = result$root, rate = 1/result$root))
}

# Two quantile constraints: P(X <= x1) = p1 and P(X <= x2) = p2
inv_gamma_params_dual <- function(x1, p1, x2, p2) {
  require(nleqslv)
  require(invgamma)
  # Solve system of equations
  result <- nleqslv::nleqslv(
    x = c(3, 2),  # initial guess for (shape, scale)
    fn = function(params) {
      shape <- params[1]
      scale <- params[2]
      c(
        invgamma::pinvgamma(x1, shape = shape, scale = scale) - p1,
        invgamma::pinvgamma(x2, shape = shape, scale = scale) - p2
      )
    }
  )
  
  return(list(shape = result$x[1], scale = result$x[2]))
}