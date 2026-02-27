# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#                       Estimation file for CO2 Break Detection
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# clean environment
rm(list=ls())


# Settings ----------------------------------------------------------------

# Get SLURM array ID
run <- commandArgs(trailingOnly = TRUE)
is_slurm <- if (length(run) > 0) TRUE else FALSE
# code is to be run on a SLURM cluster. For illustration default to setting "1" if run locally 
run_numeric <- if (is_slurm) as.numeric(run) else 1

if(is_slurm) {
  .libPaths("~/R_LIBS")
}


# load libraries
library(dplyr)
library(stringr)
library(gets)
library(getspanel)
library(Matrix)
library(mombf)

config <- list(
  ssvs_settings = data.frame(
    sis_prior = c("imom"),
    tau = sapply(c(0.05), priorp2g, q = 1, prior = c("iMom")), 
    beta_prior = c("f"),
    incl_prior = c("bern"),
    stringsAsFactors = FALSE
  ),
  gets_settings = data.frame(
    gets_lvl = c(0.05)
  ),
  date = "2026-02-27"
)

if(is_slurm) {
  dir_res <- sprintf("./results/%s/", 
                     config$date)
  dir_data <- "./CO2DriversEU_dataset_CLEAN.csv"
  source("../code/estimate_bisam_fun.R")
} else {
  dir_res <- sprintf("./output/emissions/%s/", 
                     config$date)
  dir_data <- "./data/CO2DriversEU_dataset_CLEAN.csv"
  source("./functions/estimate_bisam_fun.R")
}

dir.create(dir_res, 
           showWarnings = FALSE)


# Load and prepare data ---------------------------------------------------

data <- read.csv(dir_data)[-1]

# Group specification
EU15   <- c("Austria", "Belgium", "Germany", "Denmark", "Spain", "Finland",
            "France", "United Kingdom", "Ireland", "Italy", "Luxembourg", 
            "Netherlands", "Greece", "Portugal", "Sweden")

data_    <- data[, c('country','year','ltransport.emissions','lgdp','lgdp_sq','lpop')]
dat      <- filter(data_, country %in% EU15, year>=1995)

i_names <- unique(dat$country)
t_names <- unique(dat$year)

n <- length(i_names)
t <- length(t_names)

# BISAM -------------------------------------------------------------------
conf <- config$ssvs_settings[run_numeric, ]

# Data processing
I_INDEX <- 1
T_INDEX <- 2
Y_INDEX <- 3
DO_CONST <- FALSE
DO_INDIV_FE <- TRUE
DO_TIME_FE <- TRUE 
DO_CENTER_Y <- FALSE
DO_SCALE_Y <- FALSE
DO_CENTER_X <- FALSE
DO_SCALE_X <- FALSE

# MCMC settings
NDRAW <- 10000L
NBURN <- 2000L

# Prior settings
BETA_VARIANCE_SCALE <- 10

SIGMA2_SHAPE <- NULL
SIGMA2_RATE <- NULL
SIGMA2_HYPER_P <- 0.9

STEP_INCL_PROB <- 0.5
STEP_INCL_ALPHA <- 1
STEP_INCL_BETA <- 1

DO_INDICATOR_SATURATION <- TRUE
OUTLIER_INCL_ALPHA <- 1
OUTLIER_INCL_BETA <- 10
OUTLIER_SCALE <- 10^3

# Advanced options
DO_SPLIT_Z <- TRUE
DO_CLUSTER_S2 <- TRUE
# Set computational strategy
DO_SPARSE_COMPUTATION <- FALSE
# Check model Validity
DO_GEWEKE_TEST <- FALSE


# Prior specifications
BETA_PRIOR <- conf$beta_prior
STEP_SIZE_PRIOR <- conf$sis_prior
TAU <- conf$tau
STEP_INCL_PRIOR <- conf$incl_prior


ssvs_i <- estimate_bisam(
  data = dat,
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

dir_save <- sprintf(paste0(dir_res, "ssvs_tau-%s_prior-%s_c0-%s_C0-%s_modprior-%s_v0-%s.RDS"), 
                      conf$tau,
                      conf$sis_prior,
                      if(is.null(SIGMA2_SHAPE)){"auto"}else{as.character(SIGMA2_SHAPE)},
                      if(is.null(SIGMA2_RATE)){"auto"}else{as.character(SIGMA2_RATE)},
                      conf$incl_prior, 
                      STEP_INCL_PROB)
  
saveRDS(ssvs_i, dir_save)



# GETS --------------------------------------------------------------------

formula <- "ltransport.emissions ~ lgdp + lgdp_sq + lpop"
index   <- c("country", "year")

p.value <- config$gets_settings$gets_lvl

gets_i <- isatpanel(
  data    = dat,
  formula = as.formula(formula),
  index   = index,
  effect  = "twoways",
  iis     = TRUE,
  jsis     = FALSE,
  fesis   = TRUE, 
  t.pval  = p.value,
  print.searchinfo = FALSE)

dir_save <- sprintf(paste0(dir_res, "gets_%s.RDS"), 
                    p.value)

saveRDS(gets_i, dir_save)


# ==============================================================================
# END OF SCRIPT
# ==============================================================================