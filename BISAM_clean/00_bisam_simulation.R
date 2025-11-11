# ==============================================================================
# SETUP AND INITIALIZATION
# ==============================================================================

rm(list = ls())
library(mombf)

set.seed(192837612)

# ==============================================================================
# SIMULATION PARAMETERS
# ==============================================================================

# Prior specification
PRIOR <- "imom"

# Twostage Estimation?
DO_TWOSTAGE <- TRUE

# Data dimensions
Ni <- 10          # number of sim. observations
Nt <- 30          # number of sim. time periods
NX <- 0           # number of regressors

# Model structure
DO_CONST <- FALSE    # inclusion of a constant
DO_INDIV_FE <- TRUE     # inclusion of indiv. fixed effects
DO_TIME_FE <- FALSE     # inclusion of time fixed effects
DO_OUTLIERS <- TRUE      # inclusion of indicator saturation
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
STEP_MEAN_REL <- 2      # relative mean of size of stepshift in error.sd
STEP_SD <- 0.00         # variance of size of stepshift

# Break positions
POS_OUTL <- 0
POS_STEP <- c(41, 108, 169, 196, 221)

POS_STEP_IN_Z <- POS_STEP - 2 * (POS_STEP %/% Nt + 1) - (POS_STEP %/% Nt)
STEP_MEAN_ABS <- STEP_MEAN_REL * ERROR_SD
S2_TRUE <- ERROR_SD^2

# ==============================================================================
# DATA SIMULATION
# ==============================================================================

source("./BISAM_clean/contr_sim_breaks_fun.R")

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
if (DO_INDIV_FE & DO_TIME_FE) {
  true_fit <- fixest::feols(y ~ x1 + x2 + x3 + x4 + x5 | n + t, data = df)
} else if (DO_INDIV_FE) {
  true_fit <- fixest::feols(y ~ x1 + x2 + x3 + x4 + x5 | n, data = df)  
} else if (DO_TIME_FE) {
  true_fit <- fixest::feols(y ~ x1 + x2 + x3 + x4 + x5 | t, data = df)
} else {
  true_fit <- fixest::feols(y ~ x1 + x2 + x3 + x4 + x5, data = df)
}

# ==============================================================================
# PRIOR HYPERPARAMETERS
# ==============================================================================

if (PRIOR == "imom") {
  TAU <- priorp2g(1/12, STEP_MEAN_REL, nu = 1, prior = "iMom")
} else if (PRIOR == "mom") {
  TAU <- priorp2g(1/12, STEP_MEAN_REL, nu = 1, prior = "normalMom")
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
DO_CLUSTER_S2 <- TRUE
DO_CHECK_OUTLIER <- TRUE
# Outlier detection options
OUTLIER_INCL_ALPHA <- 1
OUTLIER_INCL_BETA <- 10
OUTLIER_SCALE <- 10
# Set computational strategy
DO_SPARSE_COMPUTATION <- FALSE
# Check model Validity
DO_GEWEKE_TEST <- FALSE

# ==============================================================================
# RUN MODEL
# ==============================================================================

source("./BISAM_clean/estimate_bisam_fun.R")

data = data
do_constant = DO_CONST
do_individual_fe = DO_INDIV_FE
do_time_fe = DO_TIME_FE
# do_step_saturation = DO_STEP_SATURATION
y_index = Y_INDEX
i_index = I_INDEX
t_index = T_INDEX
do_center_y = DO_CENTER_Y
do_scale_y = DO_SCALE_Y
do_center_x = DO_CENTER_X
do_scale_x = DO_SCALE_X
Ndraw = NDRAW
Nburn = NBURN
beta_prior = BETA_PRIOR
step_size_prior = STEP_SIZE_PRIOR
step_incl_prior = STEP_INCL_PRIOR
beta_variance_scale = BETA_VARIANCE_SCALE
sigma2_shape = SIGMA2_SHAPE
sigma2_rate = SIGMA2_RATE
sigma2_hyper_p = SIGMA2_HYPER_P
step_incl_prob = STEP_INCL_PROB
step_incl_alpha = STEP_INCL_ALPHA
step_incl_beta = STEP_INCL_BETA
step_size_scale = TAU
do_split_Z = DO_SPLIT_Z
do_cluster_s2 = DO_CLUSTER_S2
do_check_outlier = DO_CHECK_OUTLIER
outlier_incl_alpha = OUTLIER_INCL_ALPHA
outlier_incl_beta = OUTLIER_INCL_BETA
outlier_scale = OUTLIER_SCALE
do_sparse_computation = DO_SPARSE_COMPUTATION
do_geweke_test = DO_GEWEKE_TEST

mod <- estimate_bisam(
  data = data,
  do_constant = DO_CONST,
  do_individual_fe = DO_INDIV_FE,
  do_time_fe = DO_TIME_FE,
  # do_step_saturation = DO_STEP_SATURATION,
  y_index = Y_INDEX,
  i_index = I_INDEX,
  t_index = T_INDEX,
  do_center_y = DO_CENTER_Y,
  do_scale_y = DO_SCALE_Y,
  do_center_x = DO_CENTER_X,
  do_scale_x = DO_SCALE_X,
  Ndraw = NDRAW,
  Nburn = NBURN,
  beta_prior = BETA_PRIOR,
  step_size_prior = STEP_SIZE_PRIOR,
  step_incl_prior = STEP_INCL_PRIOR,
  beta_variance_scale = BETA_VARIANCE_SCALE,
  sigma2_shape = SIGMA2_SHAPE,
  sigma2_rate = SIGMA2_RATE,
  sigma2_hyper_p = SIGMA2_HYPER_P,
  step_incl_prob = STEP_INCL_PROB,
  step_incl_alpha = STEP_INCL_ALPHA,
  step_incl_beta = STEP_INCL_BETA,
  step_size_scale = TAU,
  do_split_Z = DO_SPLIT_Z,
  do_cluster_s2 = DO_CLUSTER_S2,
  do_check_outlier = DO_CHECK_OUTLIER,
  outlier_incl_alpha = OUTLIER_INCL_ALPHA,
  outlier_incl_beta = OUTLIER_INCL_BETA,
  outlier_scale = OUTLIER_SCALE,
  do_sparse_computation = DO_SPARSE_COMPUTATION,
  do_geweke_test = DO_GEWEKE_TEST
)

if (DO_TWOSTAGE) {
  pip <- mod$coefs$omega
  pip_q80 <- quantile(pip[pip > 0], 0.5)
  steps_to_check = which(pip > pip_q80)
  
  mod <- estimate_bisam(
    data = data,
    do_constant = DO_CONST,
    do_individual_fe = DO_INDIV_FE,
    do_time_fe = DO_TIME_FE,
    # do_step_saturation = DO_STEP_SATURATION,
    y_index = Y_INDEX,
    i_index = I_INDEX,
    t_index = T_INDEX,
    do_center_y = DO_CENTER_Y,
    do_scale_y = DO_SCALE_Y,
    do_center_x = DO_CENTER_X,
    do_scale_x = DO_SCALE_X,
    Ndraw = NDRAW,
    Nburn = NBURN,
    beta_prior = BETA_PRIOR,
    step_size_prior = STEP_SIZE_PRIOR,
    step_incl_prior = STEP_INCL_PRIOR,
    beta_variance_scale = BETA_VARIANCE_SCALE,
    sigma2_shape = SIGMA2_SHAPE,
    sigma2_rate = SIGMA2_RATE,
    sigma2_hyper_p = SIGMA2_HYPER_P,
    step_incl_prob = STEP_INCL_PROB,
    step_incl_alpha = STEP_INCL_ALPHA,
    step_incl_beta = STEP_INCL_BETA,
    step_size_scale = TAU,
    do_split_Z = DO_SPLIT_Z,
    do_cluster_s2 = DO_CLUSTER_S2,
    do_check_outlier = DO_CHECK_OUTLIER,
    outlier_incl_alpha = OUTLIER_INCL_ALPHA,
    outlier_incl_beta = OUTLIER_INCL_BETA,
    outlier_scale = OUTLIER_SCALE,
    steps_to_check = steps_to_check,
    do_sparse_computation = DO_SPARSE_COMPUTATION,
    do_geweke_test = DO_GEWEKE_TEST
  )
}

# ==============================================================================
# PLOTTING
# ==============================================================================

# Set up plotting parameters
par(
  mfrow = c(2, 1), 
  mar = c(1, 4.5, 1.5, 1),
  oma = c(1, 0, 2, 0),
  cex.axis = 0.9,
  cex.lab = 1,
  las = 1
)

# Color palette
COL_MAIN <- "#2166AC"    # Blue
COL_STEP <- "#B2182B"    # Red
COL_FIT <- "#1B9E77"     # Teal
COL_GRID <- "gray70"

# ---- Plot 1: Omega coefficients ----
plot(
  mod$coefs$omega, 
  type = "l", 
  col = COL_MAIN, 
  lwd = 2,
  ylim = c(0, 1),
  xlab = "",
  ylab = expression(omega),
  main = "",
  xaxt = "n",
  frame.plot = FALSE
)

# Reference lines
abline(h = c(0, 0.25, 0.5, 0.75, 1), lty = 3, col = COL_GRID, lwd = 0.8)
abline(v = 0:Ni * (Nt - 3), lty = 2, col = COL_GRID, lwd = 0.8)
abline(v = POS_STEP_IN_Z, col = COL_STEP, lty = 2, lwd = 1.5)

# Legend
legend(
  "topright", 
  legend = c(expression(omega ~ "coefficients"), "True Breaks"),
  col = c(COL_MAIN, COL_STEP),
  lty = c(1, 2),
  lwd = c(2, 1.5),
  bty = "n",
  cex = 0.85
)

mtext(
  side = 3, 
  line = 0.5, 
  text = bquote("Panel A: " ~ omega ~ "estimates"), 
  adj = 0, 
  font = 2, 
  cex = 0.95
)

# ---- Plot 2: Response variable ----
par(mar = c(3, 4.5, 0.5, 1))

plot(
  data[, 3], 
  cex = 0.6,
  pch = 19, 
  col = adjustcolor(COL_MAIN, alpha.f = 0.4),
  xlab = "",
  ylab = "Response (y)",
  main = "",
  frame.plot = FALSE
)

# Reference lines
abline(v = 0:Ni * Nt, lty = 2, col = COL_GRID, lwd = 0.8)
abline(v = POS_STEP, col = COL_STEP, lty = 2, lwd = 1.5)
abline(h = sim$true.const, lty = 3, col = COL_GRID, lwd = 0.8)
abline(h = sim$true.const + STEP_MEAN_ABS, lty = 3, col = COL_GRID, lwd = 0.8)

# Fitted lines
lines(true_fit$fitted.values, col = COL_FIT, lwd = 2)
lines(mod$fitted, col = COL_STEP, lwd = 2)

# Legend
legend(
  "topright",
  legend = c("Observed", "True Mean", "BMA-fit", "True Breaks"),
  col = c(
    adjustcolor(COL_MAIN, alpha.f = 0.9), 
    COL_FIT, 
    COL_STEP, 
    COL_STEP
  ),
  pch = c(19, NA, NA, NA),
  lty = c(NA, 1, 1, 2),
  lwd = c(NA, 2, 2, 1.5),
  bty = "n",
  cex = 0.85,
  bg = "white"
)

mtext(
  side = 3, 
  line = -0.3, 
  text = bquote("Panel B: Fitted values of" ~ y), 
  adj = 0, 
  font = 2, 
  cex = 0.95
)

# Overall title
title(
  main = bquote(
    "Step size: " ~ .(round(STEP_MEAN_ABS, 2)) ~ "units / " ~ 
      .(round(STEP_MEAN_REL, 2)) ~ sigma ~ "/" ~ tau ~ "=" ~ .(TAU) ~ 
      "/" ~ hat(sigma) ~ "=" ~ .(round(sqrt(mod$coefs$sigma2), 2))
  ),
  outer = TRUE, 
  cex.main = 1.1, 
  font.main = 2
)