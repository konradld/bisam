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
    do_sv = FALSE,               # stochastic volatility: one AR(1) log-variance process per unit over time
    sv_prior_mu_mean = NULL,     # SV level prior mean (NULL -> data-driven log(OLS variance))
    sv_prior_mu_var = 100,       # SV level prior variance
    sv_prior_phi_a = 20,         # SV persistence prior: (phi + 1) / 2 ~ Beta(sv_prior_phi_a, sv_prior_phi_b)
    sv_prior_phi_b = 1.5,
    sv_prior_sigma_shape = 2.5,  # SV vol-of-vol prior: sigma^2_eta ~ InvGamma(shape, rate)
    sv_prior_sigma_rate = 0.025,
    sv_backend = c("internal", "stochvol"), # SV sampler: hand-rolled KSC/FFBS ("internal") or stochvol::update_fast_sv ("stochvol")
    do_check_outlier = FALSE,
    outlier_incl_alpha = 1,     # P(outlier) ~ Beta(outlier_incl_alpha, outlier_incl_beta)
    outlier_incl_beta = 10,     # P(outlier) ~ Beta(outlier_incl_alpha, outlier_incl_beta)
    outlier_scale = 10,
    steps_to_check = "all",
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
  
  if(steps_to_check[1] == "all") {steps_to_check <- 1:r}
  # ============================================================================
  # PRIOR SPECIFICATION
  # ============================================================================
  
  # --- Stochastic volatility compatibility checks ---
  if (do_sv) {
    sv_backend <- match.arg(sv_backend)
    if (sv_backend == "stochvol" && !requireNamespace("stochvol", quietly = TRUE)) {
      stop("sv_backend = 'stochvol' requires the 'stochvol' package. ",
           "Install it (install.packages('stochvol')) or use sv_backend = 'internal'.")
    }
    if (do_cluster_s2) {
      warning("do_sv = TRUE overrides do_cluster_s2: stochastic volatility already provides ",
              "one unit-specific variance process over time. Setting do_cluster_s2 = FALSE.")
      do_cluster_s2 <- FALSE
    }
    if (do_geweke_test) {
      stop("do_geweke_test is only implemented for the constant-variance sampler, not for do_sv.")
    }
    if (t < 20) {
      stop("Stochastic volatility cannot be reliably estimated with T = ", t, 
           " time periods (minimum ~30, recommended 100+). ",
           "Use do_sv = FALSE or acquire more data.")
    }
    if (t < 30) {
      warning("T = ", t, " is below the recommended minimum for SV estimation. ",
              "Results will be heavily prior-dependent. ",
              "Consider using very informative priors and report sensitivity analysis.")
    }
    if (t < 50) {
      message("Note: SV estimation with T = ", t, " will have wide credible intervals. ",
              "Consider T >= 100 for reliable inference.")
    }
  }
  
  # --- Sigma^2 Prior ---
  if (is.null(sigma2_shape) | is.null(sigma2_rate)) {
    print("Sigma^2 prior is inadmissable. Using default spec. based on OLS")
    if(do_sparse_computation) {
      mod_prior <- list()
      mod_prior$coefficients <- MatrixModels:::lm.fit.sparse(X, y)
      mod_prior$residuals <- as.vector(y - X %*% mod_prior$coefficients)
    } else {
      mod_prior <- lm.fit(X, y)
    }
    if (do_cluster_s2) {
      sigma2_shape <- sigma2_rate <- numeric(n)
      for (i in 1:n) {
        n_idx <- N_idx[i,]
        res2 <- mod_prior$residuals[n_idx]^2
        Qtop <- quantile(res2, 0.90)
        s2_OLS <- sum(res2[res2 < Qtop]) / (sum(res2 < Qtop) - p / n * 0.9) # p / n is the relative share of coefs per unit
        s2_pars <- inv_gamma_params(shape = 3, s2_OLS, p = sigma2_hyper_p)
        sigma2_shape[i] <- s2_pars$shape
        sigma2_rate[i] <- s2_pars$rate
      }
    } else {
      res2 <- mod_prior$residuals^2
      Qtop <- quantile(res2, 0.90)
      s2_OLS <- sum(res2[res2 < Qtop]) / (sum(res2 < Qtop) - p * 0.9) # p * 0.9 to correct for outlier exclusion
      s2_pars <- inv_gamma_params(shape = 3, s2_OLS, p = sigma2_hyper_p)
      sigma2_shape <- s2_pars$shape
      sigma2_rate <- s2_pars$rate
    }
  }
  
  # --- Beta Prior ---
  if (beta_prior == "g" || beta_prior == "hs") {
    beta_mean <- numeric(p)
  }
  if (beta_prior == "f" || beta_prior == "f_indep") {
    if (!exists("mod_prior")) {
      if(do_sparse_computation) {
        mod_prior <- list()
        mod_prior$coefficients <- MatrixModels:::lm.fit.sparse(X, y)
        mod_prior$residuals <- as.vector(y - X %*% mod_prior$coefficients)
      } else {
        mod_prior <- lm.fit(X, y)
      }
    }
    res2 <- mod_prior$residuals^2
    select_i <- c(sapply(1:n, \(i) {
      idx <- (i-1)*t + (1:t)
      idx[res2[idx] < quantile(res2[idx], 0.90)]
    }))
    if(do_sparse_computation) {
      beta_mean <- MatrixModels:::lm.fit.sparse(X[select_i,, drop = FALSE], y[select_i,])
    } else {
      beta_mean <- lm.fit(X[select_i,, drop = FALSE], y[select_i,])$coefficients
    }
    
    if (beta_prior == "f_indep") {
      D <- if (do_sparse_computation) Diagonal(p) else diag(p)
      beta_var <- D * beta_variance_scale
      beta_var_inv <- D / beta_variance_scale
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
    incl_prior_f <- mombf::modelbinomprior(step_incl_prob) # change for Jakob!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    # incl_prior_f <- BISAM::modelbinomprior(step_incl_prob) # change for Jakob
    cat("Inclusion prior is Bernoulli(step_incl_prob) - 'step_incl_alpha' and 'step_incl_beta' have no meaning\n")
  } else if (step_incl_prior == "beta_bern") {
    incl_prior_f <- mombf::modelbbprior(step_incl_alpha, step_incl_beta) # change for Jakob !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    # incl_prior_f <- BISAM::modelbbprior(step_incl_alpha, step_incl_beta) # change for Jakob
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
  
  # --- Stochastic volatility state and hyperparameters ---
  if (do_sv) {
    # Data-driven initialisation of the per-unit log-variance level
    if (!exists("mod_prior")) {
      if (do_sparse_computation) {
        mod_prior <- list()
        mod_prior$coefficients <- MatrixModels:::lm.fit.sparse(X, y)
        mod_prior$residuals <- as.vector(y - X %*% mod_prior$coefficients)
      } else {
        mod_prior <- lm.fit(X, y)
      }
    }
    res0_mat  <- matrix(mod_prior$residuals, nrow = t, ncol = n)
    sv_mu     <- log(pmax(colMeans(res0_mat^2), 1e-6)) # per-unit level (log-variance)
    sv_phi    <- rep(0.95, n)                           # persistence
    sv_sigma2 <- rep(0.05, n)                           # vol-of-vol (innovation variance of h)
    sv_h      <- matrix(sv_mu, nrow = t, ncol = n, byrow = TRUE) # latent log-variance paths (t x n)
    s2_i      <- as.vector(exp(sv_h))                   # per-observation variance, unit-major order
    sqrt_s2_i <- sqrt(s2_i)
    if (is.null(sv_prior_mu_mean)){
      sv_prior_mu_mean <- sv_mu
    } else {
      sv_prior_mu_mean <- rep(sv_prior_mu_mean, n)
    } 
    
    sv_priors <- list(
      mu_mean   = sv_prior_mu_mean, mu_var  = sv_prior_mu_var,
      phi_a     = sv_prior_phi_a,   phi_b   = sv_prior_phi_b,
      sig_shape = sv_prior_sigma_shape, sig_rate = sv_prior_sigma_rate
    )

    # --- stochvol backend: prior spec, expert settings and extra latent state ---
    if (sv_backend == "stochvol") {
      # stochvol's fast_sv sampler places a Gamma prior on sigma^2_eta (vol-of-vol).
      # Match its prior mean to the internal InvGamma(shape, rate) mean = rate/(shape-1)
      # so the two backends are comparable; keep stochvol's default Gamma shape of 0.5.
      # Bsigma <- if (sv_prior_sigma_shape > 1) {
      #   sv_prior_sigma_rate / (sv_prior_sigma_shape - 1)
      # } else {
      #   sv_prior_sigma_rate
      # }
      
      sv_prior_spec <- list()
      for(ii in 1:n) {
        if(sv_prior_phi_a == 0 | sv_prior_phi_b == 0 | sv_prior_sigma_shape == 0 | sv_prior_sigma_rate == 0) {
          sv_prior_spec[[ii]] <- stochvol::specify_priors(
            mu     = stochvol::sv_constant(sv_prior_mu_mean[ii]),
            phi    = stochvol::sv_constant(0),
            sigma2 = stochvol::sv_constant(0.0001)
          )
        } else {
          sv_prior_spec[[ii]] <- stochvol::specify_priors(
            mu     = stochvol::sv_normal(mean = sv_prior_mu_mean[ii], sd = sqrt(sv_prior_mu_var)),
            phi    = stochvol::sv_beta(shape1 = sv_prior_phi_a, shape2 = sv_prior_phi_b),
            sigma2 = stochvol::sv_gamma(shape = sv_prior_sigma_shape, rate = sv_prior_sigma_rate) # changed that here: shape = 0.5, rate = 0.5 / Bsigma
          )
        }
      }
      
      sv_expert <- stochvol::get_default_fast_sv()
      sv_h0 <- sv_mu                            # per-unit initial log-variance state h_0
      sv_r  <- matrix(5L, nrow = t, ncol = n)   # mixture-component indicators (resampled every sweep)
    }
  }
  
  
  w_i <- logical(r)
  z_cols <- rep(1:ncol(z), n)
  w_1 <- z_cols * w_i
  
  o_i <- numeric(N)
  pip_i <- numeric(r)
  
  # --- Outlier variance inflation factor ---
  lambda_i <- rep(1, N)
  
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
  s2_store <- matrix(NA, nrow = Nstore, ncol = if (do_sv) N else if (do_cluster_s2) n else 1)
  
  colnames(b_store) <- colnames(X)
  colnames(g_store) <- colnames(Z)
  colnames(s2_store) <- if (do_sv) {
    paste0("sigma2.", iis_grid[, "Var2"], ".", iis_grid[, "Var1"])
  } else if (do_cluster_s2) {
    paste0('sigma2_', unique(data[, i_index]))
  } else {
    'sigma2'
  }
  # Stochastic volatility hyperparameter storage (one AR(1) process per unit)
  if (do_sv) {
    sv_mu_store     <- matrix(NA, nrow = Nstore, ncol = n)
    sv_phi_store    <- matrix(NA, nrow = Nstore, ncol = n)
    sv_sigma2_store <- matrix(NA, nrow = Nstore, ncol = n)
    colnames(sv_mu_store)     <- paste0("sv_mu.", n_ind)
    colnames(sv_phi_store)    <- paste0("sv_phi.", n_ind)
    colnames(sv_sigma2_store) <- paste0("sv_sigma2.", n_ind)
  }
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
  obs_with_steps <- unique(ceiling(steps_to_check / (t - 3)))
  
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
      
      # PICK OUTLIER DISTRIBUTION
      log_p1 <- log_dimom(residuals, gamma0=0, k=1, nu=3, tau=outlier_scale, sigma2=s2_i) + log(k_i)
      # log_p1 <- dnorm(residuals, 0, sqrt_outlier_scale * sqrt_s2_i, log = TRUE) + log(k_i)
      
      # Direct calculation of log probability
      log_prob_outlier <- log_p1 - matrixStats::rowLogSumExps(cbind(log_p0, log_p1))
      o_i <- rbinom(N, 1, exp(log_prob_outlier)) # works
      
      # UPDATE VARIANCE INFLATION FACTORS
      lambda_i <- ifelse(o_i == 0, 1, 2 * outlier_scale) # 2*tau*s2 is the variance of iMom(0,1,3,tau,s2)
      # lambda_i <- ifelse(o_i == 0, 1, outlier_scale)
      
    }
    
    # ==========================================================================
    # DRAW p(sigma^2 | beta, gamma, y)
    # ==========================================================================
    y_tmp <- y - Zg_i
    residuals <- y_tmp - Xb_i
    
    if (do_sv) {
      # One AR(1) log-variance (stochastic volatility) process per unit over time.
      # Feed the sampler the base residuals (outlier inflation lambda_i removed) so
      # the SV process models the base variance. total variance stays s2_i * lambda_i.
      res_for_sv <- residuals / sqrt(lambda_i)
      res_mat    <- matrix(res_for_sv, nrow = t, ncol = n)
      if (sv_backend == "stochvol") {
        sv_upd  <- update_sv_stochvol(res_mat, sv_h, sv_h0, sv_r,
                                      sv_mu, sv_phi, sv_sigma2,
                                      sv_prior_spec, sv_expert)
        sv_h0   <- sv_upd$h0
        sv_r    <- sv_upd$r
      } else {
        sv_upd  <- update_sv(res_mat, sv_h, sv_mu, sv_phi, sv_sigma2, sv_priors)
      }
      sv_h       <- sv_upd$h
      sv_mu      <- sv_upd$mu
      sv_phi     <- sv_upd$phi
      sv_sigma2  <- sv_upd$sigma2
      s2_i       <- as.vector(exp(sv_h)) # base per-observation variance (unit-major order)
    } else if (do_cluster_s2) {
      weighted_res2 <- matrix(residuals^2 / lambda_i, nrow = t, ncol = n)
      cN <- sigma2_shape + t / 2
      CN <- sigma2_rate + 0.5 * colSums(weighted_res2)
      s2_i_unique <- 1 / rgamma(n, shape = cN, rate = CN)
      s2_i <- rep(s2_i_unique, each = t)
    } else {
      cN <- sigma2_shape + N / 2
      CN <- sigma2_rate + 0.5 * sum(residuals^2 / lambda_i)
      s2_i_unique <- 1 / rgamma(1, shape = cN, rate = CN)
      s2_i[] <- s2_i_unique
    }
    sqrt_s2_i <- sqrt(s2_i)
    
    # ==========================================================================
    # DRAW p(beta | sigma^2, gamma, y)
    # ==========================================================================
    if (beta_prior == "g" || beta_prior == "f"|| beta_prior == "flasso" || beta_prior == "f_indep") {
      
      if (do_check_outlier || do_cluster_s2 || do_sv) {
        XtW <- X / (s2_i * lambda_i)
        XtWX <- crossprod(XtW, X)
        XtWy <- crossprod(XtW, y_tmp)
        
        if (beta_prior == "f_indep") {
          BN <- safe_invert(beta_var_inv + XtWX, do_sparse_computation)
        } else {
          # beta_var_inv <- crossprod(X / (beta_variance_scale * s2_i_unique), X) # g-prior unconditional on lambda_i
          beta_var_inv <- XtWX / beta_variance_scale # g-prior conditional on lambda_i
          BN <- safe_invert(beta_var_inv + XtWX, do_sparse_computation)
        }
        bN <- BN %*% (XtWy + beta_var_inv %*% beta_mean)
      }
      else {
        if (beta_prior == "f_indep") {
          XtSX <- XX / s2_i_unique
          XtSy <- crossprod(X, y_tmp) / s2_i_unique
          
          BN <- safe_invert(beta_var_inv + XtSX, do_sparse_computation)
          bN <- BN %*% (XtSy + beta_var_inv %*% beta_mean)
        } else {
          BN <- s2_i_unique * beta_variance_scale / (1 + beta_variance_scale) * XX_inv
          bN <- 1 / (1 + beta_variance_scale) * 
            (beta_variance_scale * XX_inv %*% crossprod(X, y_tmp) + beta_mean)
        }
      }
    } else {
      stop("For 'beta' only g-prior and fractional-prior is implemented!")
    }
    
    b_i <- try((bN + t(chol(BN)) %*% rnorm(p)), silent = TRUE)
    if (is(b_i, "try-error")) { b_i <- t(rmvnorm(1, as.vector(bN), as.matrix(BN))) }
    
    Xb_i <- X %*% b_i
    
    # ==========================================================================
    # DRAW p(omega | beta, sigma^2, y) - Selection Indicators
    # ==========================================================================
    y_tmp <- y - Xb_i
    
    # Standardize for model selection:
    y_tmp_sd <- y_tmp / sqrt(s2_i * lambda_i)
    if (!do_sv) {
      s2_i_tmp <- if (do_cluster_s2) s2_i_unique else rep(s2_i_unique, n)
    }
    
    for (j in obs_with_steps) {
      
      p_idx_full <- batch[j, ]
      p_idx <- p_idx_full[p_idx_full %in% steps_to_check]
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
      if (do_sv) {
        # y and Z are standardized below by the per-observation sqrt-variance, so the
        # working residual variance is ~1: use a unit-scaled inverse-gamma prior (mean 1).
        var_prior_f <- igprior(3, 2)
      } else if (do_cluster_s2) {
        var_prior_f <- igprior(sigma2_shape[j], sigma2_rate[j])
      } else {
        var_prior_f <- igprior(sigma2_shape, sigma2_rate)
      }
      
      if (do_sv) {
        # Full GLS standardization by the time-varying (and outlier-inflated) variance,
        # so the model-selection regression has unit-variance errors (phi = 1 below).
        sqrt_total_var_j <- sqrt((s2_i * lambda_i)[n_idx])
        Z_std_j <- Z[n_idx, p_idx_rand, drop = FALSE] / sqrt_total_var_j
      } else if (do_check_outlier) {
        lambda_j <- lambda_i[n_idx]
        sqrt_total_var_j <- sqrt(lambda_j)
        Z_std_j <- Z[n_idx, p_idx_rand, drop = FALSE] / sqrt_total_var_j
      } else {
        Z_std_j <- Z[n_idx, p_idx_rand, drop = FALSE]
      }
      
      # Model selection using mombf # change for Jakob
      w_i_mod <- mombf::modelSelection(
        y = y_tmp_sd[n_idx],
        x = Z_std_j,
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
        phi = 1,
        deltaini = w_i[p_idx_rand],
        initSearch = 'none',
        method = 'ALA',
        hess = "asymp",
        initpar = initpar,
        adj.overdisp = 'intercept',
        optimMethod = "auto",
        optim_maxit = 10,
        B = 10^5,
        priorVar = var_prior_f,
        priorSkew = momprior(tau = step_size_scale),
        XtXprecomp = TRUE,
        verbose = FALSE
      )
      
      # # Model selection using BISAM # change for Jakob
      # w_i_mod <- BISAM::fast_model_selection(
      #   y = y_tmp_sd[n_idx],
      #   x = Z_std_j,
      #   # groups = 1:length(p_idx),
      #   # nknots = 9,
      #   center = FALSE,
      #   scale = FALSE,
      #   # enumerate = FALSE,
      #   # includevars = rep(FALSE, length(p_idx)),
      #   niter = nn,
      #   thinning = 1,
      #   burnin = nn - 1,
      #   # family = "normal",
      #   # priorCoef = sis_prior_f@priorPars["tau"],
      #   tau = sis_prior_f@priorPars["tau"],
      #   priorDelta = incl_prior_f,
      #   phi = 1,
      #   thinit = g_i[p_idx_rand],
      #   initpar_type = 0,
      #   deltaini = w_i[p_idx_rand],
      #   # initSearch = 'none',
      #   method = 0, #'ALA', !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      #   hesstype = 1, #"asymp", !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      #   # initpar = initpar,
      #   # adj.overdisp = 'intercept',
      #   optimMethod = 2, #"auto", !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      #   optim_maxit = 10,
      #   B = 10^5,
      #   knownphi = 1,
      #   # priorVar = var_prior_f,
      #   priorSkew = BISAM::momprior(tau = step_size_scale)$tau,
      #   computation_strategy = BISAM::STANDARD(), # !!!!!!!!!!!!!!!!!!!!!!!!!!
      #   XtXprecomp = TRUE#,
      #   # verbose = FALSE
      # )
      
      w_i[p_idx_rand] <- as.logical(w_i_mod$postSample) # change for Jakob !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      # w_i[p_idx_rand] <- as.logical(w_i_mod$post_sample) # change for Jakob
      pip_i[p_idx_rand] <- w_i_mod$margpp
      
      # ========================================================================
      # DRAW p(gamma | omega, beta, sigma^2, y) - Break Magnitudes
      # ========================================================================
      if (any(w_i[p_idx_rand])) {
        g_draw <- matrix(0, nrow = ngdraw, ncol = t - 2)
        colsel <- which(w_i[p_idx_full] == TRUE)
        
        # --- Adjust for time-varying / outlier variance inflation ---
        if (do_sv) {
          # GLS transform: divide through by the per-observation total sqrt-variance so
          # the break magnitudes are drawn with known unit variance (phi_rnlp = 1).
          sqrt_var_j <- sqrt((s2_i * lambda_i)[n_idx])
          y_rnlp <- y_tmp[n_idx] / sqrt_var_j
          z_rnlp <- z[, colsel, drop = FALSE] / sqrt_var_j
          phi_rnlp <- 1
        } else if (do_check_outlier) {
          sqrt_lambda_j <- sqrt(lambda_i[n_idx])
          y_rnlp <- y_tmp[n_idx] / sqrt_lambda_j
          z_rnlp <- z[, colsel, drop = FALSE] / sqrt_lambda_j
          phi_rnlp <- s2_i_tmp[j]
        } else {
          y_rnlp <- y_tmp[n_idx]
          z_rnlp <- z[, colsel, drop = FALSE]
          phi_rnlp <- s2_i_tmp[j]
        }
        
        g_draw[1, c(colsel, t - 2)] <- rnlp_new( # the t-2 is for the sigma^2 that is thrown away if sigma is known
          y = y_rnlp,
          x = z_rnlp,
          priorCoef = sis_prior_f,
          priorGroup = sis_prior_f,
          priorVar = var_prior_f,
          isgroup = rep(FALSE, length(colsel)),
          niter = ngburn + ngdraw,
          burnin = ngburn,
          knownphi = TRUE,
          phi = phi_rnlp,
          use_thinit = TRUE,
          thinit = g_incl_i[p_idx_full][colsel]
        )
        g_incl_i[p_idx][w_i[p_idx]] <- g_draw[colsel]
      } else {
        g_draw <- matrix(0, nrow = ngdraw, ncol = t - 2)
      }
      g_i[p_idx_full] <- g_draw[-(t - 2)]  # Remove phi estimate
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
      if (do_sv) {
        s2_store[i, ]        <- s2_i        # time-varying variance per observation
        sv_mu_store[i, ]     <- sv_mu
        sv_phi_store[i, ]    <- sv_phi
        sv_sigma2_store[i, ] <- sv_sigma2
      } else {
        s2_store[i, ] <- s2_i_unique
      }
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
  
  # --- Stochastic volatility output (AR(1) log-variance process per unit) ---
  if (do_sv) {
    out$draws$sv <- list(
      "mu"     = sv_mu_store,
      "phi"    = sv_phi_store,
      "sigma2" = sv_sigma2_store
    )
    out$coefs$sv <- list(
      "mu"     = colMeans(sv_mu_store),
      "phi"    = colMeans(sv_phi_store),
      "sigma2" = colMeans(sv_sigma2_store)
    )
    # Posterior-mean variance path, reshaped to time x unit for convenience
    out$sv_variance_path <- matrix(s2_hat, nrow = t, ncol = n,
                                   dimnames = list(t_ind, n_ind))
  }
  
  class(out) <- "ism"
  
  return(out)
}

#===============================================================================
#
#                  Helper Function for inv. Gamma spec
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

#===============================================================================
#
#                  Helper Function for inverting BN
#
#===============================================================================
safe_invert <- function(m, do_sparse_computation = TRUE) {
  # Convert to dense if requested
  if (!do_sparse_computation && inherits(m, "sparseMatrix")) {
    m <- as.matrix(m)
  }
  
  tryCatch(Matrix::solve(m), error = function(e) {
    tryCatch({
      if (do_sparse_computation && inherits(m, "sparseMatrix")) {
        chol2inv(chol(as(m, "symmetricMatrix")))
      } else {
        chol2inv(chol(as.matrix(m)))
      }
    }, error = function(e2) {
      warning("Falling back to nearPD", if (!do_sparse_computation || !inherits(m, "sparseMatrix")) " (dense computation)" else " (losing sparsity)")
      Matrix::solve(nearPD(if (do_sparse_computation) m else as.matrix(m), corr = FALSE, keepDiag = TRUE)$mat)
    })
  })
}

#===============================================================================
#
#                  Helper Function for computing log density of imom
#
#===============================================================================
log_dimom <- function(x, gamma0, k, nu, tau, sigma2) {
  stopifnot(k > 0, nu > 0, tau > 0, sigma2 > 0)
  
  z      <- x - gamma0
  result <- rep(-Inf, length(z))
  valid  <- z != 0
  
  if (any(valid)) {
    lz        <- log(abs(z[valid]))
    log_scale <- 0.5 * log(tau * sigma2)
    
    const <- log(k) +
      (nu / 2) * log(tau * sigma2) -
      lgamma(nu / (2 * k))
    
    result[valid] <- const -
      (nu + 1) * lz -
      exp(-2 * k * (lz - log_scale))
  }
  result
}

#===============================================================================
#
#     Helper Function for one Gibbs update of the Stochastic Volatility block
#
#===============================================================================
# One full sweep of the stochastic-volatility Gibbs step, run in parallel for
# every unit (column of the residual matrix). Each unit has its own AR(1)
# log-variance process
#
#     r_{i,t} = exp(h_{i,t} / 2) * u_{i,t},   u_{i,t} ~ N(0, 1)
#     h_{i,t} = mu_i + phi_i * (h_{i,t-1} - mu_i) + eta_{i,t},  eta ~ N(0, sigma2_i)
#
# Arguments
#   res_mat : t x n matrix of residuals (time in rows, units in columns)
#   h       : t x n matrix of current log-variance states
#   mu, phi, sigma2 : length-n vectors of current hyperparameters
#   priors  : list(mu_mean, mu_var, phi_a, phi_b, sig_shape, sig_rate)
# Returns list(h, mu, phi, sigma2).
update_sv <- function(res_mat, h, mu, phi, sigma2, priors) {
  t <- nrow(res_mat)
  n <- ncol(res_mat)
  
  # --- Kim, Shephard & Chib (1998) 7-component mixture for log(chi^2_1) ---
  ksc_prob <- c(0.00730, 0.10556, 0.00002, 0.04395, 0.34001, 0.24566, 0.25750)
  ksc_mean <- c(-11.40039, -5.24321, -9.83726, 1.50746, -0.65098, 0.52478, -2.35859)
  ksc_var  <- c(5.79596, 2.61369, 5.17950, 0.16735, 0.64009, 0.34023, 1.26261)
  
  # Log squared residuals (offset avoids log(0) for near-zero residuals)
  ystar <- log(res_mat^2 + 1e-7)
  
  # ---------------------------------------------------------------------------
  # 1. Sample mixture component indicators s given the current states h
  # ---------------------------------------------------------------------------
  maxlog <- matrix(-Inf, t, n)
  logw   <- vector("list", 7)
  for (jc in 1:7) {
    lw <- log(ksc_prob[jc]) - 0.5 * log(ksc_var[jc]) -
      0.5 * (ystar - h - ksc_mean[jc])^2 / ksc_var[jc]
    logw[[jc]] <- lw
    maxlog <- pmax(maxlog, lw)
  }
  wsum  <- matrix(0, t, n)
  wlist <- vector("list", 7)
  for (jc in 1:7) {
    wj <- exp(logw[[jc]] - maxlog)
    wlist[[jc]] <- wj
    wsum <- wsum + wj
  }
  u_draw <- matrix(runif(t * n), t, n) * wsum
  csum   <- matrix(0, t, n)
  chosen <- matrix(1L, t, n)
  found  <- matrix(FALSE, t, n)
  for (jc in 1:7) {
    csum <- csum + wlist[[jc]]
    pick <- (!found) & (u_draw <= csum)
    chosen[pick] <- jc
    found[pick]  <- TRUE
  }
  m_sel <- matrix(ksc_mean[chosen], t, n) # mixture means per observation
  v_sel <- matrix(ksc_var[chosen],  t, n) # mixture variances per observation
  
  # ---------------------------------------------------------------------------
  # 2. Sample the log-variance states h via FFBS (vectorized over units)
  #    Observation:  (ystar - m_sel) = h + N(0, v_sel)
  #    State:        h_t = mu + phi (h_{t-1} - mu) + N(0, sigma2)
  # ---------------------------------------------------------------------------
  yy <- ystar - m_sel
  
  filt_mean <- matrix(0, t, n)
  filt_var  <- matrix(0, t, n)
  pred_mean <- matrix(0, t, n)
  pred_var  <- matrix(0, t, n)
  
  a <- mu                     # predicted state mean at t = 1 (stationary mean)
  P <- sigma2 / (1 - phi^2)   # predicted state variance at t = 1 (stationary var)
  for (tt in 1:t) {
    pred_mean[tt, ] <- a
    pred_var[tt, ]  <- P
    Fv <- P + v_sel[tt, ]
    K  <- P / Fv
    filt_mean[tt, ] <- a + K * (yy[tt, ] - a)
    filt_var[tt, ]  <- P * (1 - K)
    # one-step-ahead prediction for the next time point
    a <- mu + phi * (filt_mean[tt, ] - mu)
    P <- phi^2 * filt_var[tt, ] + sigma2
  }
  
  h_new <- matrix(0, t, n)
  h_new[t, ] <- filt_mean[t, ] + sqrt(filt_var[t, ]) * rnorm(n)
  if (t >= 2) {
    for (tt in (t - 1):1) {
      J      <- filt_var[tt, ] * phi / pred_var[tt + 1, ]
      m_star <- filt_mean[tt, ] + J * (h_new[tt + 1, ] - pred_mean[tt + 1, ])
      P_star <- filt_var[tt, ] * (1 - J * phi)
      h_new[tt, ] <- m_star + sqrt(pmax(P_star, 0)) * rnorm(n)
    }
  }
  h <- h_new
  
  # ---------------------------------------------------------------------------
  # 3. Sample the AR(1) hyperparameters (per unit, vectorized over units)
  # ---------------------------------------------------------------------------
  h1   <- h[1, ]
  hlag <- h[-t, , drop = FALSE]  # h_{t-1}, t = 2..T
  hcur <- h[-1, , drop = FALSE]  # h_t,     t = 2..T
  
  # --- sigma2 (vol-of-vol): inverse-gamma full conditional ---
  hcen_lag <- sweep(hlag, 2, mu)
  hcen_cur <- sweep(hcur, 2, mu)
  eta      <- hcen_cur - sweep(hcen_lag, 2, phi, `*`)
  SSE      <- (1 - phi^2) * (h1 - mu)^2 + colSums(eta^2)
  sigma2   <- 1 / rgamma(n,
                         shape = priors$sig_shape + t / 2,
                         rate  = priors$sig_rate + 0.5 * SSE)
  
  # --- mu (level): Gaussian full conditional ---
  sum_term <- colSums(hcur - sweep(hlag, 2, phi, `*`)) # sum_{t=2}^T (h_t - phi h_{t-1})
  prec <- 1 / priors$mu_var +
    (1 - phi^2) / sigma2 +
    (t - 1) * (1 - phi)^2 / sigma2
  numr <- priors$mu_mean / priors$mu_var +
    (1 - phi^2) / sigma2 * h1 +
    (1 - phi) / sigma2 * sum_term
  mu <- numr / prec + sqrt(1 / prec) * rnorm(n)
  
  # --- phi (persistence): Metropolis-within-Gibbs, Beta prior on (phi + 1)/2 ---
  hcen_lag <- sweep(hlag, 2, mu)
  hcen_cur <- sweep(hcur, 2, mu)
  Sxx <- colSums(hcen_lag^2)
  Sxy <- colSums(hcen_lag * hcen_cur)
  Sxx[Sxx <= 0] <- 1e-8
  mean_prop <- Sxy / Sxx
  sd_prop   <- sqrt(sigma2 / Sxx)
  phi_prop  <- rnorm(n, mean_prop, sd_prop)
  
  h1 <- h[1, ]
  log_target <- function(ph) {
    # stationary-initial-state term + Beta((phi+1)/2) prior; the AR(2..T)
    # likelihood cancels against the Gaussian proposal in the acceptance ratio.
    stat  <- 0.5 * log(1 - ph^2) - 0.5 * (1 - ph^2) * (h1 - mu)^2 / sigma2
    prior <- (priors$phi_a - 1) * log((1 + ph) / 2) +
      (priors$phi_b - 1) * log((1 - ph) / 2)
    stat + prior
  }
  valid   <- phi_prop > -1 & phi_prop < 1
  phi_eval <- ifelse(valid, phi_prop, 0)          # dummy for invalid draws (masked out below)
  log_acc <- log_target(phi_eval) - log_target(phi)
  log_acc[!valid] <- -Inf
  accept  <- valid & (log(runif(n)) < log_acc)
  phi[accept] <- phi_prop[accept]
  
  list(h = h, mu = mu, phi = phi, sigma2 = sigma2)
}

#===============================================================================
#
#     Helper Function: one SV Gibbs sweep via the stochvol package backend
#
#===============================================================================
# Note on parameterisation: stochvol works with sigma (the SD of the log-variance
# innovations); this wrapper takes/returns sigma2 (the variance) to stay
# interchangeable with update_sv(). It also threads the extra latent state
# stochvol needs across sweeps: h0 (initial log-variance) and r (mixture-component
# indicators)
#
# Arguments
#   res_mat    : t x n matrix of residuals (time in rows, units in columns)
#   h          : t x n matrix of current log-variance states
#   h0         : length-n vector of current initial states h_0
#   r          : t x n integer matrix of current mixture indicators
#   mu, phi, sigma2 : length-n vectors of current hyperparameters (sigma2 = variance)
#   prior_spec : stochvol prior object from specify_priors()
#   expert     : stochvol expert settings from get_default_fast_sv()
# Returns list(h, h0, r, mu, phi, sigma2).
update_sv_stochvol <- function(res_mat, h, h0, r, mu, phi, sigma2, prior_spec, expert) {
  t <- nrow(res_mat)
  n <- ncol(res_mat)
  offset <- 1e-7  # guards log(0) for near-zero residuals (matches update_sv)
  
  for (i in 1:n) {
    # log_data2 <- res_mat[, i]^2 + offset)
    y_raw <- res_mat[, i]
    y_raw[abs(y_raw) < offset] <- offset * sign(y_raw[abs(y_raw) < offset] + offset)
    
    para <- list(mu = mu[i], 
                 phi = phi[i], 
                 sigma = sqrt(sigma2[i]), 
                 nu = Inf, # only necessary for t-model
                 rho = 0, # only necessary for leverage model
                 beta = NA,
                 latent0 = h0[i]) 
    
    latent <- h[, i]

    upd <- tryCatch({
      stochvol::svsample_fast_cpp(
        y = y_raw,
        startpara = para,
        startlatent = latent,
        priorspec = prior_spec[[i]],
        fast_sv = expert
      )
    }, error = function(e) {
      # Fallback: use svsample with 1 draw
      warning(sprintf("svsample_fast_cpp failed for unit %d: %s. Using svsample fallback.", i, e$message))
      sv_fit <- stochvol::svsample(
        y_raw,
        draws = 1,
        burnin = 0,
        startpara = list(mu = mu[i], phi = phi[i], sigma = sqrt(sigma2[i])),
        startlatent = latent,
        priorspec = prior_spec[[i]],
        quiet = TRUE
      )
      list(
        latent = as.numeric(sv_fit$latent[[1]]),
        latent0 = as.numeric(sv_fit$latent0[[1]]),
        para = matrix(c(
          as.numeric(sv_fit$para[, "mu"]),
          as.numeric(sv_fit$para[, "phi"]),
          as.numeric(sv_fit$para[, "sigma"])
        ), nrow = 1, dimnames = list(NULL, c("mu", "phi", "sigma")))
      )
    })

    h[, i]    <- upd$latent
    h0[i]     <- upd$latent0
    # r[, i]    <- r[, i] # handled internally by stochvol
    mu[i]     <- upd$para[, "mu"]
    phi[i]    <- upd$para[, "phi"]
    sigma2[i] <- upd$para[, "sigma"]^2
  }
  
  list(h = h, h0 = h0, r = r, mu = mu, phi = phi, sigma2 = sigma2)
}
