# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#                       Method Comparison for Break Detection
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# ==============================================================================
# SETUP AND INITIALIZATION
# ==============================================================================
rm(list = ls())

# Get SLURM array ID
run <- commandArgs(trailingOnly = TRUE)
is_slurm <- if (length(run) > 0) TRUE else FALSE
# code is to be run on a SLURM cluster. For illustration default to setting "1" if run locally 
run_numeric <- if (is_slurm) as.numeric(run) else 1

if(is_slurm) {
  .libPaths("~/R_LIBS")
}

library(Matrix)
library(mombf)

config <- expand.grid(
  sis_prior = c("imom"),
  gets_lvl = c(0.01),
  rel_effect = c(3),
  tau = c(priorp2g(0.01, 1)),
  
  Nt = seq(from = 20, to = 100, by = 5),
  Ni = c(10),
  
  number_reps = 1:100,
  setup = c("sparse"),
  date = "2026-07-15_no-FEs_fixed-breaks_NT",
  stringsAsFactors = FALSE
)
conf <- config[run_numeric,]

# ==============================================================================
# SIMULATION PARAMETERS
# ==============================================================================

# Prior specification
PRIOR <- conf$sis_prior

# Data dimensions
Ni <- conf$Ni          # number of sim. observations
Nt <- conf$Nt          # number of sim. time periods
NX <- 0           # number of regressors

# Model structure
DO_CONST <- TRUE    # inclusion of a constant
DO_INDIV_FE <- FALSE     # inclusion of indiv. fixed effects
DO_TIME_FE <- FALSE     # inclusion of time fixed effects
DO_OUTLIERS <- FALSE      # inclusion of indicator saturation
DO_STEP_SATURATION <- TRUE      # inclusion of stepshift saturation
DO_INDICATOR_SATURATION <- FALSE
# Outlier and break parameters
P_OUTL <- 0.0    # probability of outlier in a Series
P_STEP <- 0.0    # probability of a stepshift in a Series
# Error distribution
ERROR_SD <- 1   # standard deviation of the error
# Outlier characteristics
OUTL_MEAN <- 0   # mean of size of outlier
# Stepshift characteristics
STEP_MEAN_REL <- conf$rel_effect    # relative mean of size of stepshift in error.sd

# Break positions
POS_OUTL <- 0

# Sample random breaks in the N_STEPS observations
if(conf$setup == "sparse") {
  N_STEPS <- c(1:4)
  # N_STEPS <- sample(1:Ni, ceiling(0.01 * Nt * Ni), replace = T)
  
} else if (conf$setup == "dense") {
  N_STEPS <- c(1:8, 2, 4, 6, 8)
  # N_STEPS <- sample(1:Ni, ceiling(0.05 * Nt * Ni), replace = T)
  
} else if (is.numeric(conf$setup)) {
  N_STEPS <- sample(1:Ni, conf$setup, replace = T)
} else {
  stop("Simulation setup not available.")
}

chk <- TRUE
while(chk) {
  POS_STEP_IN_Z <- sapply(N_STEPS, \(x) sample(1:(Nt - 3) + (x - 1) * (Nt - 3), 1))
  if(!any(duplicated(POS_STEP_IN_Z))) chk <- FALSE
}

POS_STEP <- POS_STEP_IN_Z + 2 * (N_STEPS) + ((N_STEPS) - 1)
STEP_MEAN_ABS <- STEP_MEAN_REL * ERROR_SD
S2_TRUE <- ERROR_SD^2

# ==============================================================================
# DATA SIMULATION
# ==============================================================================

source("./R/contr_sim_breaks_fun.R")
source("./R/estimate_bisam_fun.R")
source("./R/pip_window_fun.R")

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

data <- sim$data

# To save results in
results <- list()

# =========================== BISAM ======================================.=====

# Data processing
I_INDEX <- 1
T_INDEX <- 2
Y_INDEX <- 3
DO_CENTER_Y <- FALSE
DO_SCALE_Y <- FALSE
DO_CENTER_X <- FALSE
DO_SCALE_X <- FALSE

# MCMC settings
NDRAW <- 5000L
NBURN <- 1000L

# Prior settings
BETA_VARIANCE_SCALE <- 100

SIGMA2_SHAPE <- NULL
SIGMA2_RATE <- NULL
SIGMA2_HYPER_P <- 0.9

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
# Outlier detection options
OUTLIER_INCL_ALPHA <- 1
OUTLIER_INCL_BETA <- 10
OUTLIER_SCALE <- 10
# Set computational strategy
DO_SPARSE_COMPUTATION <- FALSE
# Check model Validity
DO_GEWEKE_TEST <- FALSE

if(conf$tau == "auto") {
  if (PRIOR == "imom") {
    TAU <- priorp2g(0.05, STEP_MEAN_REL, nu = 1, prior = "iMom")
  } else if (PRIOR == "mom") {
    TAU <- STEP_MEAN_REL^2 / 2
  } else {
    stop("selected prior not implemented")
  }
} else {
  TAU <- as.numeric(conf$tau)
}


# ==============================================================================
# RUN MODEL
# ==============================================================================

results$b_ssvs <- estimate_bisam(
  data = data,
  do_constant = DO_CONST,
  do_individual_fe = DO_INDIV_FE,
  do_time_fe = DO_TIME_FE,
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
  do_check_outlier = DO_INDICATOR_SATURATION,
  outlier_incl_alpha = OUTLIER_INCL_ALPHA,
  outlier_incl_beta = OUTLIER_INCL_BETA,
  outlier_scale = OUTLIER_SCALE,
  do_sparse_computation = DO_SPARSE_COMPUTATION,
  do_geweke_test = DO_GEWEKE_TEST
)

#===============================================================================
# Collect ground truth
#===============================================================================
# The BISAM object does NOT contain the ground truth, so pull it from the
# simulator. All downstream metrics (power, FDR, joint recovery) are computed in
# the analysis script from the saved draws + this truth.

# True break set (names "sis.i.t"), plus sizes / indices / net effect.
tr_breaks <- rownames(sim$tr.idx)

#===============================================================================
# Per-replicate posterior summaries  (ALL computation happens here)
#===============================================================================
# Reduce the Nstore x K inclusion draws to the posterior quantities that
# characterise model-selection consistency (Johnson & Rossell 2010, 2012). The
# analysis script only stacks and plots these numbers. Expectations E[.] and
# probabilities P(.) below are over the MCMC draws of omega.
#
# Headline (strong -> marginal):
#   p_true = P( active set == true set )     posterior mass on the true model  -> 1
#   power  = E[ TP / |true| ]                posterior mean power              -> 1
#   fdr    = E[ FP / |selected| ]            posterior mean false-discovery    -> 0
#
# Diagnostics (graceful weaker convergences that never collapse to 0/1 abruptly):
#   e_fn / e_fp / hamming  expected # missed / spurious / total errors         -> 0
#   p_all_true             P( all true breaks included )   recall-side one-sided-> 1
#   p_no_fp                P( no spurious break )          precision-side       -> 1
#   pip_true / pip_false   mean marginal PIP on true / false candidates      -> 1 / 0

omega_draws <- results$b_ssvs$draws$omega          # Nstore x K, 0/1
true_cols   <- colnames(omega_draws) %in% tr_breaks
n_true      <- length(tr_breaks)

TP_draw   <- rowSums(omega_draws[, true_cols, drop = FALSE])   # true steps included per draw
nsel_draw <- rowSums(omega_draws)                              # total selected per draw
FP_draw   <- nsel_draw - TP_draw                               # spurious steps per draw
FN_draw   <- n_true - TP_draw                                  # missed true steps per draw

pip <- results$b_ssvs$coefs$omega   # marginal PIPs = colMeans(omega_draws), named

metrics <- data.frame(
  # --- headline: strong (joint) + marginal ---
  p_true    = mean(TP_draw == n_true & FP_draw == 0),
  power     = mean(TP_draw / n_true),
  fdr       = mean(ifelse(nsel_draw > 0, FP_draw / nsel_draw, 0)),
  # --- diagnostics: expected error counts (soft Hamming distance) ---
  e_fn      = mean(FN_draw),
  e_fp      = mean(FP_draw),
  hamming   = mean(FP_draw + FN_draw),
  # --- diagnostics: one-sided recovery probabilities ---
  p_all_true = mean(TP_draw == n_true),   # no false negatives in the draw
  p_no_fp    = mean(FP_draw == 0),        # no false positives in the draw
  # --- diagnostics: marginal PIP separation ---
  pip_true  = mean(pip[true_cols]),
  pip_false = mean(pip[!true_cols]),
  # --- extra diagnostics ---
  post_size = mean(nsel_draw),   # posterior mean model size (vs n_true)
  n_true    = n_true,
  stringsAsFactors = FALSE
)

#===============================================================================
# Assemble the object to save
#===============================================================================
# The plot script consumes only $meta + $metrics (all computed above). We also
# keep the raw MCMC draws for any later analysis the summaries don't cover
# (marginal PIP profiles, step magnitudes, model-size distribution). The step
# draws are overwhelmingly zero / FALSE (spike-and-slab), so they are stored as
# sparse matrices to keep files small.

output <- list(
  meta = data.frame(
    prior      = conf$sis_prior,
    tau        = TAU,
    gets_lvl   = conf$gets_lvl,
    Ni         = Ni,
    Nt         = Nt,
    setup      = conf$setup,
    rep        = conf$number_reps,
    rel_effect = conf$rel_effect,
    step_mean  = STEP_MEAN_ABS,
    error_sd   = ERROR_SD,
    Ndraw      = NDRAW,
    Nburn      = NBURN,
    Nstore     = as.integer(results$b_ssvs$meta$MCMC[["Nstore"]]),
    stringsAsFactors = FALSE
  ),

  # ---- per-replicate posterior summaries (what the plot script consumes) ----
  metrics = metrics,

  # ---- ground truth (not available inside the BISAM object) -----------------
  truth = list(
    breaks = tr_breaks,     # names "sis.i.t" of the true steps
    tr_idx = sim$tr.idx     # size / index / rel_net_eff, rownames = breaks
  ),

  # ---- BISAM posterior draws (the "right thing" to keep) --------------------
  bisam = list(
    # Nstore x K Bernoulli inclusion indicators (the draws behind the PIPs)
    omega  = Matrix::Matrix(results$b_ssvs$draws$omega, sparse = TRUE),
    # Nstore x K step magnitudes (spike-and-slab coefficient draws)
    sis    = Matrix::Matrix(results$b_ssvs$draws$sis,   sparse = TRUE),
    # Nstore x 1 error-variance draws (diagnostics / uncertainty)
    sigma2 = results$b_ssvs$draws$sigma2,
    # marginal PIPs = colMeans(omega), kept for convenience
    pip    = results$b_ssvs$coefs$omega
  )
)

#===============================================================================
# Save Results
#===============================================================================

folder_path <- sprintf("./output/simulation/bisam_convergence-%0.2f_bisam_prior-%s_tau-%s/",
                       conf$gets_lvl, conf$sis_prior, conf$tau)

if (!dir.exists(folder_path)) {dir.create(folder_path, recursive = TRUE)}

file_name <- sprintf("breaksize-%0.1fSD_breaknumber-%s_Ni-%s_Nt-%s_rep%0.0f.RDS", conf$rel_effect, conf$setup, conf$Ni, conf$Nt, conf$number_reps)

saveRDS(output, file = paste0(folder_path, file_name))


#===============================================================================
# End of File
#===============================================================================