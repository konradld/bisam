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

library(stringr)
library(gets)
library(getspanel)
library(Matrix)
library(mombf)
library(glmnet)  # Added for adaptive Lasso

config <- expand.grid(
  sis_prior = c("imom"),
  gets_lvl = c(0.05,0.01),
  rel_effect = c(1, 1.5, 2, 3, 6, 10),
  tau = mombf::priorp2g(0.05, 1, nu = 1, prior = "iMom"),
  number_reps = 1:100,
  date = "2025-11-18_sparse",
  stringsAsFactors = FALSE
)
conf <- config[run_numeric,]

# ==============================================================================
# SIMULATION PARAMETERS
# ==============================================================================

# Prior specification
PRIOR <- conf$sis_prior

# Data dimensions
Ni <- 10          # number of sim. observations
Nt <- 30          # number of sim. time periods
NX <- 0           # number of regressors

# Model structure
DO_CONST <- TRUE    # inclusion of a constant
DO_INDIV_FE <- FALSE     # inclusion of indiv. fixed effects
DO_TIME_FE <- FALSE     # inclusion of time fixed effects
DO_OUTLIERS <- FALSE      # inclusion of indicator saturation
DO_STEP_SATURATION <- TRUE      # inclusion of stepshift saturation
DO_INDICATOR_SATURATION = FALSE
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

# Sample random breaks in the first N_STEPS observations
N_STEPS <- c(1:4)
POS_STEP_IN_Z <- sapply(N_STEPS, \(x) sample(1:(Nt - 3) + (x - 1) * (Nt - 3), 1))
POS_STEP <- POS_STEP_IN_Z + 2 * (N_STEPS) + ((N_STEPS) - 1)
STEP_MEAN_ABS <- STEP_MEAN_REL * ERROR_SD
S2_TRUE <- ERROR_SD^2

# ==============================================================================
# DATA SIMULATION
# ==============================================================================

if (is_slurm) {
  source("../code/contr_sim_breaks_fun.R")
  source("../code/estimate_bisam_fun.R")
  source("../code/pip_window_fun.R")
} else {
  source("./Functions/contr_sim_breaks_fun.R")
  source("./Functions/estimate_bisam_fun.R")
  source("./Functions/pip_window_fun.R")
}

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

data<- sim$data
dat <- as.data.frame(data)

# To save results in
results <- list()

# =========================== GETS =======================================.=====

results$gets <- list() # initiate new sub-list

formula  <- paste('y',
                  paste(colnames(dat)[grepl('x\\d+',colnames(data))],collapse = '+'),
                  sep = '~')

if(NX==0 & DO_CONST){formula <- 'y~c'} #------------- careful not always true !!

index = c("n","t")

# Break analysis:
sig_lvl <- conf$gets_lvl

# run model
res_i <- isatpanel(
  data = cbind(dat, c = 1),
  formula = as.formula(formula),
  index = index,
  effect = "none",
  iis = DO_INDICATOR_SATURATION,
  jsis = FALSE,
  fesis = TRUE,
  t.pval = sig_lvl,
  print.searchinfo = FALSE
)

# fix names
gets_breaks <- res_i$isatpanel.result$mean.results %>%
  rownames() %>%
  str_replace_all(c(
    "x" = "beta.x",
    "time" = "tfe.",
    "id" = "ife.",
    "(fesis|sis)" = "sis.",
    "iis" = "iis."
  ))

gets_coefs <- res_i$isatpanel.result$mean.results
rownames(gets_coefs) <- gets_breaks

results$gets[[as.character(sig_lvl)]] <- gets_coefs

# =========================== ADAPTIVE LASSO =================================.=

results$alasso <- list()

# Create design matrix matching BISAM setup
# Extract response variable
y_vec <- data[, 3]  # Y_INDEX = 3

# Create step saturation indicators
n_obs <- nrow(data)
n_units <- Ni
n_time <- Nt

# Initialize design matrix
X_base <- matrix(0, nrow = n_obs, ncol = 0)

# Add constant if needed
if (DO_CONST) {
  X_base <- cbind(X_base, const = 1)
}

# Add regressors if any
if (NX > 0) {
  X_cols <- grep("^x\\d+$", colnames(data), value = TRUE)
  X_base <- cbind(X_base, data[, X_cols, drop = FALSE])
}

# Create step saturation matrix (SIS)
# For each unit, create step indicators for each time period
sis_matrix <- matrix(0, nrow = n_obs, ncol = n_units * (n_time - 1))
sis_names <- character(n_units * (n_time - 1))

col_idx <- 1
for (i in 1:n_units) {
  unit_rows <- which(data[, 1] == i)  # I_INDEX = 1
  
  for (t_start in 2:n_time) {
    # Create step indicator: 0 before t_start, 1 from t_start onwards
    time_vals <- data[unit_rows, 2]  # T_INDEX = 2
    sis_matrix[unit_rows, col_idx] <- ifelse(time_vals >= t_start, 1, 0)
    sis_names[col_idx] <- paste0("sis.", i, ".", t_start)
    col_idx <- col_idx + 1
  }
}

colnames(sis_matrix) <- sis_names

# Combine base and SIS matrices
X_full <- cbind(X_base, sis_matrix)

# Remove constant from penalization
n_base_cols <- ncol(X_base)
sis_start_col <- n_base_cols + 1

# Step 1: Initial estimation (Ridge regression for stability)
# Use ridge to get initial weights
cv_ridge <- cv.glmnet(
  x = X_full[, sis_start_col:ncol(X_full), drop = FALSE],
  y = y_vec,
  alpha = 0,
  standardize = TRUE,
  intercept = DO_CONST
)

ridge_coef <- as.vector(coef(cv_ridge, s = "lambda.min"))
# Remove intercept from coefficients
if (DO_CONST) {
  ridge_coef <- ridge_coef[-1]
}

# Step 2: Create adaptive weights
# Weight = 1 / |ridge coefficient|^gamma, gamma typically = 1
gamma <- 1
adaptive_weights <- 1 / (abs(ridge_coef) + 1e-8)^gamma

# Step 3: Run adaptive Lasso
# Create penalty factor vector (0 for unpenalized, adaptive weights for SIS)
penalty_factors <- rep(0, ncol(X_full))
penalty_factors[sis_start_col:ncol(X_full)] <- adaptive_weights

# Run cross-validated Lasso with adaptive weights
cv_alasso <- cv.glmnet(
  x = X_full,
  y = y_vec,
  alpha = 1,  # Lasso
  standardize = TRUE,
  intercept = FALSE,  # Already included in X_full if needed
  penalty.factor = penalty_factors,
  nfolds = 10
)

# Extract selected breaks (using lambda.1se for more conservative selection)
alasso_coef <- coef(cv_alasso, s = "lambda.1se")
alasso_coef_vec <- as.vector(alasso_coef)[-1]  # Remove intercept

# Get selected SIS indicators
selected_sis_idx <- which(alasso_coef_vec[sis_start_col:ncol(X_full)] != 0)
alasso_breaks <- sis_names[selected_sis_idx]

# Store results
results$alasso$selected_breaks <- alasso_breaks
results$alasso$coefficients <- alasso_coef_vec[sis_start_col:ncol(X_full)]
results$alasso$lambda_min <- cv_alasso$lambda.min
results$alasso$lambda_1se <- cv_alasso$lambda.1se
results$alasso$cv_error <- min(cv_alasso$cvm)

# Also store results using lambda.min (less conservative)
alasso_coef_min <- coef(cv_alasso, s = "lambda.min")
alasso_coef_vec_min <- as.vector(alasso_coef_min)[-1]
selected_sis_idx_min <- which(alasso_coef_vec_min[sis_start_col:ncol(X_full)] != 0)
alasso_breaks_min <- sis_names[selected_sis_idx_min]
results$alasso$selected_breaks_min <- alasso_breaks_min

# =========================== BISAM ======================================.=====

data <- sim$data

# Data processing
I_INDEX <- 1
T_INDEX <- 2
Y_INDEX <- 3
DO_CENTER_Y <- FALSE
DO_SCALE_Y <- FALSE
DO_CENTER_X <- FALSE
DO_SCALE_X <- FALSE

# MCMC settings
NDRAW <- 10000L
NBURN <- 2000L

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
    TAU <- priorp2g(1/12, STEP_MEAN_REL, nu = 1, prior = "iMom")
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
# Start Analysis
#===============================================================================
#helper function
make_sis_names <- function (window_breaks) {
  win_size <- unique(as.numeric(window_breaks$end_time) -
                       as.numeric(window_breaks$start_time) + 1)
  
  if (length(win_size) > 1) {
    stop("incorrect window construction - non uniform window length!")
  }
  
  paste("sis",
        rep(window_breaks$unit, each = win_size),
        mapply(
          \(x1, x2) x1:x2,
          as.numeric(window_breaks$start_time),
          as.numeric(window_breaks$end_time)
        ),
        sep = ".")
}

# comparison ----

tr_breaks     <- rownames(sim$tr.idx)
tr_breaks_index_matrix <- str_match(tr_breaks, "sis\\.(\\d+)\\.(\\d+)")

ssvs_breaks   <- names(results$b_ssvs$coefs$omega[results$b_ssvs$coefs$omega > 0.5])

ssvs_breaks <- pip_window(results$b_ssvs, win_size = 1, op = ">=", pip_threshold = 0.50)
ssvs_breaks <- make_sis_names(ssvs_breaks)

ssvs_breaks_t2 <- pip_window(results$b_ssvs, win_size = 2, op = ">=", pip_threshold = 0.75)
ssvs_breaks_t2 <- make_sis_names(ssvs_breaks_t2)

gets_breaks <- rownames(results$gets[[1]])[grepl("iis.+|sis.+",rownames(results$gets[[1]]))]

# Adaptive Lasso breaks (using lambda.1se)
alasso_breaks <- results$alasso$selected_breaks

# Adaptive Lasso breaks (using lambda.min - less conservative)
alasso_breaks_min <- results$alasso$selected_breaks_min

all_breaks  <- str_sort(unique(c(tr_breaks, ssvs_breaks, gets_breaks, 
                                 alasso_breaks, alasso_breaks_min, 
                                 ssvs_breaks_t2)), numeric = TRUE)

#check immediate neighbourhood
tr_breaks_nbh_ <- unique(c(
  paste(
    "sis",
    tr_breaks_index_matrix[, 2],
    as.numeric(tr_breaks_index_matrix[, 3]) - 1,
    sep = "."
  ),
  paste(
    "sis",
    tr_breaks_index_matrix[, 2],
    as.numeric(tr_breaks_index_matrix[, 3]) + 1,
    sep = "."
  )
))
tr_breaks_nbh  <- tr_breaks_nbh_[!tr_breaks_nbh_%in%tr_breaks]

# Extended comparison matrix to include adaptive Lasso
break_comparison <- matrix(NA, nrow = length(all_breaks), ncol = 24)
colnames(break_comparison) <- c(
  "true", "ssvs", "gets", "alasso", "alasso_min",
  "tr.ssvs", "tr.gets", "tr.alasso", "tr.alasso_min",
  "fp.ssvs", "fp.gets", "fp.alasso", "fp.alasso_min",
  "fn.ssvs", "fn.gets", "fn.alasso", "fn.alasso_min",
  "tr.all3", "fp.all3", "fn.all3",
  "ssvs_1nn_fp", "gets_1nn_fp", "alasso_1nn_fp",
  "ssvs_t3_tr"
)
rownames(break_comparison) <- all_breaks

# which are true/found
break_comparison[,1] <- ifelse(all_breaks %in% tr_breaks, 1, 0)
break_comparison[,2] <- ifelse(all_breaks %in% ssvs_breaks, 1, 0)
break_comparison[,3] <- ifelse(all_breaks %in% gets_breaks, 1, 0)
break_comparison[,4] <- ifelse(all_breaks %in% alasso_breaks, 1, 0)
break_comparison[,5] <- ifelse(all_breaks %in% alasso_breaks_min, 1, 0)

# which are true positive
break_comparison[,6] <- (break_comparison[,2] + break_comparison[,1]) == 2
break_comparison[,7] <- (break_comparison[,3] + break_comparison[,1]) == 2
break_comparison[,8] <- (break_comparison[,4] + break_comparison[,1]) == 2
break_comparison[,9] <- (break_comparison[,5] + break_comparison[,1]) == 2

# which are false positive
break_comparison[,10] <- (break_comparison[,2] - break_comparison[,1]) == 1
break_comparison[,11] <- (break_comparison[,3] - break_comparison[,1]) == 1
break_comparison[,12] <- (break_comparison[,4] - break_comparison[,1]) == 1
break_comparison[,13] <- (break_comparison[,5] - break_comparison[,1]) == 1

# which are false negative
break_comparison[,14] <- (break_comparison[,2] - break_comparison[,1]) == -1
break_comparison[,15] <- (break_comparison[,3] - break_comparison[,1]) == -1
break_comparison[,16] <- (break_comparison[,4] - break_comparison[,1]) == -1
break_comparison[,17] <- (break_comparison[,5] - break_comparison[,1]) == -1

# which tr/fp/fn are in common across all three methods (SSVS, GETS, ALASSO with lambda.1se)
break_comparison[,18] <- (break_comparison[,6] + break_comparison[,7] + break_comparison[,8]) == 3
break_comparison[,19] <- (break_comparison[,10] + break_comparison[,11] + break_comparison[,12]) == 3
break_comparison[,20] <- (break_comparison[,14] + break_comparison[,15] + break_comparison[,16]) == 3

# 1-nearest neighbor false positives
break_comparison[,21] <- (ifelse(all_breaks %in% tr_breaks_nbh, 1, 0) + break_comparison[,10]) == 2
break_comparison[,22] <- (ifelse(all_breaks %in% tr_breaks_nbh, 1, 0) + break_comparison[,11]) == 2
break_comparison[,23] <- (ifelse(all_breaks %in% tr_breaks_nbh, 1, 0) + break_comparison[,12]) == 2

all_t3 <- ifelse(all_breaks %in% ssvs_breaks_t2, 1, 0)
break_comparison[,24] <- (all_t3 + break_comparison[,1]) == 2  


#===============================================================================
# Save Results
#===============================================================================

if (is_slurm) {
  dir.create(sprintf("./results/%s/", conf$date), showWarnings = FALSE)
  
  folder_path <- sprintf("./results/%s/gets_bisam_alasso_comparison_gets-%0.2f_bisam_prior-%s_tau-%s/",
                         conf$date, conf$gets_lvl, conf$sis_prior, conf$tau)
} else {
  folder_path <- sprintf("./Simulations/gets_bisam_alasso_comparison_gets-%0.2f_bisam_prior-%s_tau-%s/",
                         conf$gets_lvl, conf$sis_prior, conf$tau)
}

if (!dir.exists(folder_path)) {dir.create(folder_path)}

file_name <- sprintf("breaksize-%0.1fSD_rep%0.0f.RDS", conf$rel_effect, conf$number_reps)

saveRDS(break_comparison, file = paste0(folder_path, file_name))


#===============================================================================
# End of File
#===============================================================================