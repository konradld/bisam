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
library(stochvol)

config <- list(
  ssvs_settings = expand.grid(
    sis_prior = c("imom"),
    tau = sapply(c(0.05), priorp2g, q = 1, prior = c("iMom")), 
    check_outl = c(FALSE),
    out_scale = c(10),
    beta_prior = c("f"),
    incl_prior = c("bern"),
    
    do_sv = c(TRUE),
    sv_version = c("internal", "stochvol"), 
    sv_mu_mean = c("auto"),
    sv_mu_var = c(10),
    sv_phi_a = c(5),
    sv_phi_b = c(5),
    sv_sigma_shape = c(100, 1000),
    sv_sigma_rate = c(1000, 100000),
    
    stringsAsFactors = FALSE
  ),
  gets_lvl = c(0.05), # just use the setting as done by Koch et al (2022)
  date = "2026-07-20_SV-stochvol"
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
  source("./R/estimate_bisam_fun.R")
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


if(run_numeric <= nrow(config$ssvs_settings)) {
  
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
  
  OUTLIER_INCL_ALPHA <- 1
  OUTLIER_INCL_BETA <- 10
  OUTLIER_SCALE <- conf$out_scale
  
  DO_SV <- conf$do_sv
  
  if(conf$sv_mu_mean == "auto") mean_mu <- NULL else mean_mu <- as.numeric(conf$sv_mu_mean)
  
  SV_PRIOR_MU_MEAN <- mean_mu
  SV_PRIOR_MU_VAR <- conf$sv_mu_var
  SV_PRIOR_PHI_A <- conf$sv_phi_a
  SV_PRIOR_PHI_B <- conf$sv_phi_b
  SV_PRIOR_SIGMA_SHAPE <- conf$sv_sigma_shape
  SV_PRIOR_SIGMA_RATE <- conf$sv_sigma_rate

  SV_BACKEND <- conf$sv_version
  
  
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
  DO_INDICATOR_SATURATION <- conf$check_outl
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
    
    do_sv = DO_SV,
    sv_prior_mu_mean = SV_PRIOR_MU_MEAN, 
    sv_prior_mu_var = SV_PRIOR_MU_VAR,
    sv_prior_phi_a = SV_PRIOR_PHI_A,
    sv_prior_phi_b = SV_PRIOR_PHI_B,
    sv_prior_sigma_shape = SV_PRIOR_SIGMA_SHAPE,
    sv_prior_sigma_rate = SV_PRIOR_SIGMA_RATE,
    sv_backend = SV_BACKEND, 
    
    do_check_outlier = DO_INDICATOR_SATURATION,
    outlier_incl_alpha = OUTLIER_INCL_ALPHA,
    outlier_incl_beta = OUTLIER_INCL_BETA,
    outlier_scale = OUTLIER_SCALE,
    do_sparse_computation = DO_SPARSE_COMPUTATION,
    do_geweke_test = DO_GEWEKE_TEST
  )
  
  dir_save <- sprintf(paste0(dir_res, 
                             "ssvs_sv-%s_mu-mean-%s-var-%s_phi-a-%s-b-%s_sigma-shape-%s-rate-%s.RDS"), 
                      conf$sv_version,
                      conf$sv_mu_mean,
                      conf$sv_mu_var,
                      conf$sv_phi_a,
                      conf$sv_phi_b,
                      conf$sv_sigma_shape,
                      conf$sv_sigma_rate
                      )
  saveRDS(ssvs_i, dir_save)
  
} else {
  
  # GETS --------------------------------------------------------------------
  
  formula <- "ltransport.emissions ~ lgdp + lgdp_sq + lpop"
  index   <- c("country", "year")
  
  p.value <- config$gets_lvl[run_numeric - nrow(config$ssvs_settings)]
  
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
}

