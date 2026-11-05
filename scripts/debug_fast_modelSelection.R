# ==============================================================================
# SETUP AND INITIALIZATION
# ==============================================================================

rm(list = ls())
library(bisam)

set.seed(12345)

setup_type <- "sparse" # or "dense"

# ==============================================================================
# SIMULATION PARAMETERS
# ==============================================================================

# Prior specification
PRIOR <- "imom"

# Twostage Estimation?
DO_TWOSTAGE <- FALSE

# Data dimensions
Ni <- 10          # number of sim. observations
Nt <- 30          # number of sim. time periods
NX <- 0           # number of regressors

# Model structure
DO_CONST <- FALSE    # inclusion of a constant
DO_INDIV_FE <- TRUE     # inclusion of indiv. fixed effects
DO_TIME_FE <- TRUE     # inclusion of time fixed effects
DO_OUTLIERS <- FALSE      # inclusion of indicator saturation
DO_STEP_SATURATION <- TRUE      # inclusion of stepshift saturation

# Outlier and break parameters
P_OUTL <- 0.0    # probability of outlier in a Series
P_STEP <- 0.0    # probability of a stepshift in a Series

# Error distribution
ERROR_SD <- 1   # standard deviation of the error

# Outlier characteristics
OUTL_MEAN <- 0   # mean of size of outlier
OUTL_SD <- 0     # variance of size of outlier

# Stepshift characteristics
STEP_MEAN_REL <- 3      # relative mean of size of stepshift in error.sd
STEP_SD <- 0.00         # variance of size of stepshift

# Break positions
POS_OUTL <- 0

if(setup_type == "sparse") {
  POS_STEP <- c(43, 108, 169, 221)
} else if (setup_type == "dense") {
  POS_STEP <- c(7, 22, 43, 80, 108, 115, 127, 144, 169, 190, 200, 221)
} else {
  stop("Simulation setup not available.")
}


POS_STEP_IN_Z <- POS_STEP - 2 * (POS_STEP %/% Nt + 1) - (POS_STEP %/% Nt)
STEP_MEAN_ABS <- STEP_MEAN_REL * ERROR_SD
S2_TRUE <- ERROR_SD^2

# ==============================================================================
# DATA SIMULATION
# ==============================================================================

source("./scripts/00_sim_breaks_for_viz.R")

sim <- contr_sim_breaks(
  n = Ni, 
  t = Nt, 
  nx = NX, 
  iis = DO_OUTLIERS, 
  sis = DO_STEP_SATURATION,
  const = DO_CONST, 
  ife = DO_INDIV_FE, 
  tfe = DO_TIME_FE,
  pos.outl = POS_OUTL, 
  pos.step = POS_STEP,
  outl.mean = OUTL_MEAN, 
  step.mean = STEP_MEAN_ABS,
  error.sd = ERROR_SD
)

# ==============================================================================
# HELPER FUNCTIONS
# ==============================================================================

make_Z <- function(n, t, i_index = 1, t_index = 2) {
  z <- lower.tri(matrix(1, nrow = t, ncol = t), diag = TRUE)[, -c(1:2, t)]
  mode(z) <- "integer"
  Z <- kronecker(diag(n), z)
  mode(Z) <- "integer"
  
  sis_grid <- expand.grid(unique(data[, t_index])[-c(1:2, t)], 
                          unique(data[, i_index]))
  iis_grid <- expand.grid(unique(data[, t_index]), 
                          unique(data[, i_index]))
  colnames(Z) <- paste0("sis.", sis_grid[, "Var2"], ".", sis_grid[, "Var1"])
  
  return(Z)
}

# ==============================================================================
# TRUE MODEL FIT
# ==============================================================================

data <- sim$data
y <- data[, 3]
Z <- make_Z(Ni, Nt)
true_fit <- lm(y ~ Z[, POS_STEP_IN_Z])

df <- cbind(as.data.frame(data), as.data.frame(Z[, POS_STEP_IN_Z]))
names(df) <- c("n","t","y", paste0("x",1:length(POS_STEP)))
if (requireNamespace("fixest", quietly = TRUE)) {
  if (DO_INDIV_FE & DO_TIME_FE) {
    true_fit <- fixest::feols(as.formula(paste0("y~", paste0("x", 1:length(POS_STEP), collapse = "+"), "| n + t")) , data = df)
  } else if (DO_INDIV_FE) {
    true_fit <- fixest::feols(as.formula(paste0("y~", paste0("x", 1:length(POS_STEP), collapse = "+"), "| n")) , data = df)  
  } else if (DO_TIME_FE) {
    true_fit <- fixest::feols(as.formula(paste0("y~", paste0("x", 1:length(POS_STEP), collapse = "+"), "| t")), data = df)
  } else {
    true_fit <- fixest::feols(as.formula(paste0("y~", paste0("x", 1:length(POS_STEP), collapse = "+"))), data = df)
  }
}

# ==============================================================================
# PRIOR HYPERPARAMETERS
# ==============================================================================

if (PRIOR == "imom") {
  TAU <- priorp2g(0.05, 1, nu = 1, prior = "iMom")
} else if (PRIOR == "mom") {
  TAU <- priorp2g(0.05, 1, nu = 1, prior = "normalMom")
} else {
  stop("selected prior not implemented")
}

# ==============================================================================
# MODEL ESTIMATION PARAMETERS
# ==============================================================================

# Data processing
I_INDEX <- 1
T_INDEX <- 2
Y_INDEX <- 3
DO_CENTER_Y <- FALSE
DO_SCALE_Y <- FALSE
DO_CENTER_X <- FALSE
DO_SCALE_X <- FALSE

# MCMC settings
NDRAW <- 2000L
NBURN <- 500L

# Prior settings
BETA_VARIANCE_SCALE <- 10

SIGMA2_SHAPE <- NULL
SIGMA2_RATE <- NULL
SIGMA2_HYPER_P <- 0.90

STEP_INCL_PROB <- 0.5
STEP_INCL_ALPHA <- 1
STEP_INCL_BETA <- 1

# Prior specifications
BETA_PRIOR <- "f"
STEP_SIZE_PRIOR <- PRIOR
STEP_INCL_PRIOR <- "bern"

# Advanced options
DO_SPLIT_Z <- TRUE
DO_CLUSTER_S2 <- FALSE
DO_CHECK_OUTLIER <- TRUE
# Outlier detection options
OUTLIER_INCL_ALPHA <- 1
OUTLIER_INCL_BETA <- 10
OUTLIER_SCALE <- 20
# Set computational strategy
DO_SPARSE_COMPUTATION <- FALSE
# Check model Validity
DO_GEWEKE_TEST <- FALSE
# Integrated likelihood for iMOM break sizes and break-size Gibbs iterations
IMOM_MARGINAL <- "ala_legacy"
RNLP_ITER <- 25L

# ==============================================================================
# RUN MODEL: arguments, exactly as estimate_bisam() receives them
# ==============================================================================

args <- list(
  data = data, do_constant = DO_CONST, do_individual_fe = DO_INDIV_FE, do_time_fe = DO_TIME_FE,
  y_index = Y_INDEX, i_index = I_INDEX, t_index = T_INDEX,
  do_center_y = DO_CENTER_Y, do_scale_y = DO_SCALE_Y, do_center_x = DO_CENTER_X, do_scale_x = DO_SCALE_X,
  Ndraw = NDRAW, Nburn = NBURN, beta_prior = BETA_PRIOR, step_size_prior = STEP_SIZE_PRIOR,
  imom_marginal = IMOM_MARGINAL, rnlp_iter = RNLP_ITER, step_incl_prior = STEP_INCL_PRIOR,
  beta_variance_scale = BETA_VARIANCE_SCALE, sigma2_shape = SIGMA2_SHAPE, sigma2_rate = SIGMA2_RATE,
  sigma2_hyper_p = SIGMA2_HYPER_P, step_incl_prob = STEP_INCL_PROB, step_incl_alpha = STEP_INCL_ALPHA,
  step_incl_beta = STEP_INCL_BETA, step_size_scale = TAU, do_split_Z = DO_SPLIT_Z,
  do_cluster_s2 = DO_CLUSTER_S2, do_check_outlier = DO_CHECK_OUTLIER,
  outlier_incl_alpha = OUTLIER_INCL_ALPHA, outlier_incl_beta = OUTLIER_INCL_BETA,
  outlier_scale = OUTLIER_SCALE, do_sparse_computation = DO_SPARSE_COMPUTATION,
  do_geweke_test = DO_GEWEKE_TEST
)

# Full run with the packaged function (uncomment to compare):
# res <- do.call(estimate_bisam, args)

# Put every argument (given or default) into the workspace, as estimate_bisam() would see it
for (.nm in names(formals(estimate_bisam))) {
  assign(.nm, if (.nm %in% names(args)) args[[.nm]] else eval(formals(estimate_bisam)[[.nm]]))
}
rm(.nm)

# Internal helpers of the bisam package used in the body below
bisam_ms_gibbs     <- bisam:::bisam_ms_gibbs
bisam_rnlp         <- bisam:::bisam_rnlp
sparse_lm_coef     <- bisam:::sparse_lm_coef
inv_gamma_params   <- bisam:::inv_gamma_params
safe_invert        <- bisam:::safe_invert
log_dimom          <- bisam:::log_dimom
update_sv          <- bisam:::update_sv
update_sv_stochvol <- bisam:::update_sv_stochvol

# ==============================================================================
# BODY OF estimate_bisam() (bisam 0.2.0), inlined for line-by-line debugging.
# The Gibbs loop (i) and the unit loop (j) are unrolled: the loop headers are
# replaced by `i <- 1` and `j <- obs_with_steps[1]` and their closing braces are
# commented out, so the code below runs ONE Gibbs iteration for ONE unit. Change
# i / j and re-run the parts you need.
# ==============================================================================
# ============================================================================
# GEWEKE DIAGNOSTIC MODE SETUP
# ============================================================================
if (do_geweke_test) {
  beta_variance_scale <- 1.234
  Ndraw <- 100000
  sigma2_shape <- 3
  sigma2_rate <- 1
  beta_prior <- "g"
}

# ============================================================================
# INITIAL SETUP
# ============================================================================
cl <- as.call(c(quote(estimate_bisam), args))
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
  sv_backend <- match.arg(sv_backend, eval(formals(estimate_bisam)$sv_backend))
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
    mod_prior$coefficients <- sparse_lm_coef(X, y)
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
      mod_prior$coefficients <- sparse_lm_coef(X, y)
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
    beta_mean <- sparse_lm_coef(X[select_i,, drop = FALSE], y[select_i,])
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
  if (!requireNamespace("glmnet", quietly = TRUE)) {
    stop("beta_prior = 'flasso' requires the 'glmnet' package. Install it with install.packages('glmnet').")
  }
  lasso_model <- glmnet::glmnet(X, y, alpha = 1, intercept = FALSE)
  best_lambda <- glmnet::cv.glmnet(X, y, alpha = 1)$lambda.min
  beta_mean <- as.array(stats::coef(lasso_model, s = best_lambda))[-1]
}
# beta_var_inv <- XX / beta_variance_scale

# --- Structural Break Prior (MOM/iMOM) ---
# Passed to the compiled sampler as a code: 0 = product MOM (r = 1), 1 = product iMOM
step_size_prior <- step_size_prior[1]
if (step_size_prior == "imom") {
  sis_prior_code <- 1L
} else if (step_size_prior == "mom") {
  sis_prior_code <- 0L
} else {
  stop("sis prior not implemented, choose step_size_prior = c('mom','imom')")
}
# iMOM integrated likelihood used in the break-indicator Gibbs step (ignored for MOM):
#   "ala_legacy" ALA at the mode found by mombf 3.5.4's coordinate iteration from a ridge
#                start (default; fastest, results as in mombf)
#   "ala"        ALA (posterior mode + BIC-type term, as mombf) at the near-global posterior
#                mode (more accurate integrated likelihood, ~1.5-3.5x slower)
#   "laplace"    Laplace approximation (posterior mode + Hessian) at the near-global mode;
#                closer to the exact integrated likelihood than ALA
imom_marginal <- match.arg(imom_marginal, eval(formals(estimate_bisam)$imom_marginal))
imom_method_code <- switch(imom_marginal, ala_legacy = 0L, ala = 1L, laplace = 2L)

# --- Inclusion Prior ---
# Code for the compiled sampler: 1 = Bernoulli(step_incl_prob), 2 = Beta-Bernoulli(alpha, beta)
step_incl_prior <- step_incl_prior[1]
if (step_incl_prior == "bern") {
  if (step_incl_prob <= 0 || step_incl_prob > 1) stop("step_incl_prob must lie in (0, 1]")
  incl_prior_code <- 1L
  cat("Inclusion prior is Bernoulli(step_incl_prob) - 'step_incl_alpha' and 'step_incl_beta' have no meaning\n")
} else if (step_incl_prior == "beta_bern") {
  incl_prior_code <- 2L
  cat("Inclusion prior is Beta-Bernoulli(step_incl_alpha, step_incl_beta) - 'step_incl_prob' has no meaning\n")
} else {
  stop("step_incl_prior not implemented, choose step_incl_prior = c('bern','beta_bern')")
}

# --- Outlier Prior ---
outlier_prior <- match.arg(outlier_prior, eval(formals(estimate_bisam)$outlier_prior))
if(!outlier_prior %in% c("imom", "norm")) {
  stop("outlier prior not implemented, choose outlier_prior = c('imom','norm')")
}
sqrt_outlier_scale <- sqrt(outlier_scale)


# ============================================================================
# STARTING VALUES
# ============================================================================
b_i <- beta_mean
g_i <- numeric(r)
s2_i <- rep(1 / rgamma(1, shape = sigma2_shape, rate = sigma2_rate), N)
sqrt_s2_i <- sqrt(s2_i)

# --- Stochastic volatility state and hyperparameters ---
if (do_sv) {
  # Data-driven initialisation of the per-unit log-variance level
  if (!exists("mod_prior")) {
    if (do_sparse_computation) {
      mod_prior <- list()
      mod_prior$coefficients <- sparse_lm_coef(X, y)
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
lambda <- ifelse(outlier_prior == "imom", 
                 2 * outlier_scale, # 2*tau*s2 is variance of iMom(0,1,3,tau,s2) 
                 outlier_scale)     # tau*s2 for N(0, tau*s2)
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
# Per-unit Gram matrices X_u'X_u (column u of G_unit = vec of a p x p matrix), used to form
# X'WX cheaply when the weights are constant within units apart from flagged outliers
if (!do_sparse_computation && !do_sv) {
  G_unit <- vapply(seq_len(n), function(u) as.vector(crossprod(X[N_idx[u, ], , drop = FALSE])), numeric(p * p))
}
sqrt_outlier_scale <- sqrt(outlier_scale)

nn <- 2 # number of model-selection Gibbs iterations (burn-in nn - 1)

# Break sizes are drawn with the latent truncation Gibbs sampler (Rossell & Telesca 2017,
# Alg. 2). The break indicators are drawn with the sizes integrated out, so:
#  * if a unit's set of breaks is unchanged, its previous sizes are already a draw from their
#    conditional posterior and one Gibbs iteration (invariant for the new conditional, given
#    the updated beta and sigma^2) is a valid update;
#  * if the set changed, the previous sizes belong to another model, so a converged draw is
#    needed: rnlp_iter iterations, new breaks starting at their conditional posterior mode.
# Earlier code: 3 iterations always, new breaks started from stale values. In full runs with
# many breaks (MOM prior) this inflated sigma^2 by ~5% and shifted break sizes and PIPs.
# With this scheme, rnlp_iter = 25 matched a 200-iteration reference within Monte Carlo
# error; rnlp_iter = 12 still left sigma^2 about 1% too high.
rnlp_iter <- as.integer(rnlp_iter)
if (rnlp_iter < 1) stop("rnlp_iter must be >= 1")

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

# Per-batch column indices and dense break-design blocks (fixed over the sampler)
p_idx_list <- lapply(seq_len(nrow(batch)), function(j) {
  p_idx_full <- batch[j, ]
  p_idx_full[p_idx_full %in% steps_to_check]
})
n_idx_list <- lapply(seq_len(nrow(batch)), function(j) if (do_split_Z) N_idx[j, ] else 1:N)
Z_block <- lapply(seq_len(nrow(batch)), function(j) {
  zb <- as.matrix(Z[n_idx_list[[j]], p_idx_list[[j]], drop = FALSE])
  storage.mode(zb) <- "double"
  zb
})
storage.mode(z) <- "double"

# ============================================================================
# GIBBS SAMPLER
# ============================================================================
pb <- txtProgressBar(min = 0, max = Ndraw, style = 3)

i <- 1  # for (i in (1 - Nburn):Nstore) {
  
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
    
    if(outlier_prior == "imom") {
      log_p1 <- log_dimom(residuals, gamma0 = 0, k = 1, nu = 3, tau = outlier_scale, sigma2 = s2_i) + log(k_i)
    } else if(outlier_prior == "norm") {
      log_p1 <- dnorm(residuals, 0, sqrt_outlier_scale * sqrt_s2_i, log = TRUE) + log(k_i)
    }
    
    # Direct calculation of log probability
    log_prob_outlier <- log_p1 - matrixStats::rowLogSumExps(cbind(log_p0, log_p1))
    o_i <- rbinom(N, 1, exp(log_prob_outlier)) # works
    
    # UPDATE VARIANCE INFLATION FACTORS
    lambda_i <- ifelse(o_i == 0, 1, lambda)
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
  if (!(beta_prior %in% c("g", "f", "flasso", "f_indep"))) {
    stop("For 'beta' only g-prior and fractional-prior is implemented!")
  }
  if (!do_sparse_computation) {
    # Dense path: one Cholesky factor R of the posterior precision (R'R = Prec) gives
    # both the posterior mean (two triangular solves) and the draw bN + R^{-1} z.
    weighted <- do_check_outlier || do_cluster_s2 || do_sv
    if (weighted) {
      w_obs <- 1 / (s2_i * lambda_i)
      if (!do_sv) {
        # weights = 1/s2 per unit (or overall), times 1/lambda for flagged outliers:
        # X'WX = sum_u X_u'X_u / s2_u + outlier-row corrections (per-unit Gram matrices precomputed)
        s2_unit <- if (do_cluster_s2) s2_i_unique else rep(s2_i_unique, n)
        XtWX <- matrix(G_unit %*% (1 / s2_unit), p, p)
        jj <- which(lambda_i != 1)
        if (length(jj) > 0) {
          Xo <- X[jj, , drop = FALSE]
          XtWX <- XtWX + crossprod(Xo * (w_obs[jj] - 1 / s2_i[jj]), Xo)
        }
      } else {
        XtWX <- crossprod(X * sqrt(w_obs))
      }
      XtWy <- crossprod(X, w_obs * y_tmp)
    } else {
      XtWX <- XX / s2_i_unique
      XtWy <- crossprod(X, y_tmp) / s2_i_unique
    }
    if (beta_prior == "f_indep") {
      prior_prec <- beta_var_inv
    } else {
      prior_prec <- XtWX / beta_variance_scale  # g-prior conditional on lambda_i (and s2)
    }
    Prec <- prior_prec + XtWX
    rhs <- XtWy + prior_prec %*% beta_mean
    R_prec <- tryCatch(chol(Prec), error = function(e) NULL)
    if (!is.null(R_prec)) {
      bN <- backsolve(R_prec, backsolve(R_prec, rhs, transpose = TRUE))
      b_i <- bN + backsolve(R_prec, rnorm(p))
    } else {
      BN <- safe_invert(Prec, do_sparse_computation)
      bN <- BN %*% rhs
      b_i <- t(rmvnorm(1, as.vector(bN), as.matrix(BN)))
    }
  } else {
    # Sparse path (unchanged)
    if (do_check_outlier || do_cluster_s2 || do_sv) {
      XtW <- X / (s2_i * lambda_i)
      XtWX <- crossprod(XtW, X)
      XtWy <- crossprod(XtW, y_tmp)
      
      if (beta_prior == "f_indep") {
        BN <- safe_invert(beta_var_inv + XtWX, do_sparse_computation)
      } else {
        beta_var_inv <- XtWX / beta_variance_scale # g-prior conditional on lambda_i
        BN <- safe_invert(beta_var_inv + XtWX, do_sparse_computation)
      }
      bN <- BN %*% (XtWy + beta_var_inv %*% beta_mean)
    } else {
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
    b_i <- try((bN + t(chol(BN)) %*% rnorm(p)), silent = TRUE)
    if (is(b_i, "try-error")) { b_i <- t(rmvnorm(1, as.vector(bN), as.matrix(BN))) }
  }
  
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
  
  j <- obs_with_steps[1]  # for (j in obs_with_steps) {
    
    p_idx_full <- batch[j, ]
    p_idx <- p_idx_list[[j]]
    p_idx_rand <- p_idx  # Keep order intact (can be randomized: sample(p_idx))
    n_idx <- n_idx_list[[j]]
    
    # Note: the residual-variance prior plays no role here. phi is known (= 1 after
    # standardisation) both in the model selection and in the break-size draws.
    if (do_sv) {
      # Full GLS standardization by the time-varying (and outlier-inflated) variance,
      # so the model-selection regression has unit-variance errors (phi = 1 below).
      sqrt_total_var_j <- sqrt((s2_i * lambda_i)[n_idx])
      Z_std_j <- Z_block[[j]] / sqrt_total_var_j
    } else if (do_check_outlier) {
      sqrt_total_var_j <- sqrt(lambda_i[n_idx])
      Z_std_j <- Z_block[[j]] / sqrt_total_var_j
    } else {
      Z_std_j <- Z_block[[j]]
    }
    
    # Model selection: Gibbs sweeps over this unit's break indicators with ALA
    # integrated likelihoods (compiled core, known phi = 1)
    w_i_mod <- bisam_ms_gibbs(
      y = y_tmp_sd[n_idx],
      x = Z_std_j,
      prior = sis_prior_code,
      tau = step_size_scale,
      phi = 1,
      priorDelta = incl_prior_code,
      p_incl = step_incl_prob,
      bb_a = step_incl_alpha,
      bb_b = step_incl_beta,
      niter = nn,
      burnin = nn - 1L,
      deltaini = w_i[p_idx_rand],
      imom_method = imom_method_code
    )
    
    w_i[p_idx_rand] <- w_i_mod$postSample[nrow(w_i_mod$postSample), ]
    pip_i[p_idx_rand] <- w_i_mod$margpp
    
    # ========================================================================
    # DRAW p(gamma | omega, beta, sigma^2, y) - Break Magnitudes
    # ========================================================================
    g_draw <- numeric(length(p_idx_full))
    if (any(w_i[p_idx_rand])) {
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
      
      # Non-local posterior draw of the selected break sizes (known phi). Breaks that stay
      # in the model start from their previous draw; breaks that just entered (previous
      # value 0) start at the maximiser of their conditional posterior given the others.
      g_prev <- g_i[p_idx_full][colsel]
      is_new <- g_prev == 0                       # break just entered (was excluded)
      set_changed <- any((g_i[p_idx_full] != 0) != w_i[p_idx_full])  # a break entered or left
      n_rnlp <- if (set_changed) rnlp_iter else 1L
      g_draw[colsel] <- bisam_rnlp(
        y = as.vector(y_rnlp),
        x = z_rnlp,
        prior = sis_prior_code,
        tau = step_size_scale,
        phi = phi_rnlp,
        niter = n_rnlp,
        burnin = n_rnlp - 1L,
        thinit = ifelse(is_new, 1, g_prev),
        newflag = is_new
      )
    }
    g_i[p_idx_full] <- g_draw
  # }  # end of unit loop (j)
  
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
# }  # end of Gibbs loop (i)

# ============================================================================
# POST-PROCESSING
# ============================================================================
# NOTE: with the loops unrolled only draw i is stored; the other rows of the *_store matrices are NA.
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
# (matrixStats::colQuantiles uses the same default as quantile(), type 7)
ci_probs <- c((1 - cred_int) / 2, 1 - (1 - cred_int) / 2)
ci_b  <- matrixStats::colQuantiles(b_store,  probs = ci_probs, drop = FALSE)
ci_g  <- matrixStats::colQuantiles(g_store,  probs = ci_probs, drop = FALSE)
ci_s2 <- matrixStats::colQuantiles(s2_store, probs = ci_probs, drop = FALSE)
estimates <- data.frame(
  MAP = c(beta_hat, gamma_hat, s2_hat),
  CI_Lower = c(ci_b[, 1], ci_g[, 1], ci_s2[, 1]),
  CI_Upper = c(ci_b[, 2], ci_g[, 2], ci_s2[, 2]),
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


str(out, max.level = 1)
