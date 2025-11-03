par_iter_fun <- function(
    Ni=n,
    Nt=t,
    NX=nx,
    DO_OUTLIERS=iis,
    DO_STEP_SATURATION=sis,
    DO_CONST=const,
    DO_INDIV_FE = ife,
    DO_TIME_FE = tfe,
    POS_OUTL = pos.outl,
    POS_STEP = pos.step,
    OUTL_MEAN=outl.mean,
    STEP_MEAN_REL=step.mean,
    ERROR_SD=error.sd,
    PRIOR = prior_fam,
){
  
# ==============================================================================
# SETUP AND INITIALIZATION
# ==============================================================================

rm(list = ls())
library(mombf)

set.seed(192837612)

# ==============================================================================
# SIMULATION PARAMETERS
# ==============================================================================

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
  
  
  data<- sim$data
  dat <- as.data.frame(data)
  
  # To save results in
  results <- list()
  
  # =========================== GETS =======================================.=====
  
  results$gets <- list() # initiate new sub-list
  
  significance_levels <- c(0.05,0.01)
  formula  <- paste('y',
                    paste(colnames(dat)[grepl('x\\d+',colnames(data))],collapse = '+'),
                    sep = '~')
  if(NX==0 & DO_CONST){formula <- 'y~c'} #--------------------------------------------- careful not always true !!
  
  index = c("n","t")
  
  # Break analysis:
  for(sig_lvl in significance_levels){
    
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
    
    results$gets[[as.character(sig_lvl)]]<-gets_coefs
  }
  
  # =========================== Bayesian SSVS ==============================.=====

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
  BETA_VARIANCE_SCALE <- 1000
  
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
  
  if (PRIOR == "imom") {
    TAU <- priorp2g(1/12, STEP_MEAN_REL, nu = 1, prior = "iMom")
  } else if (PRIOR == "mom") {
    TAU <- priorp2g(1/12, STEP_MEAN_REL, nu = 1, prior = "normalMom")
  } else {
    stop("selected prior not implemented")
  }
  # ==============================================================================
  # RUN MODEL
  # ==============================================================================
  
  source("./BISAM_clean/estimate_bisam_fun.R")
  results$b_ssvs_5 <- estimate_bisam(
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
  
  #=============================================================================
  # Start Analysis
  #=============================================================================
  #helper function
  make_sis_names <- function (window_breaks) {
    paste("sis",
          window_breaks$unit,
          mapply(
            \(x1, x2) x1:x2,
            as.numeric(window_breaks$start_time),
            as.numeric(window_breaks$end_time)
          ),
          sep = ".")
  }
  
  # comparison 0.01----
  results$b_ssvs_1 <- results$b_ssvs_5 #------------------------------------------ Only when iterating over Tau!!!
  
  tr_breaks     <- rownames(sim$tr.idx)
  tr_breaks_index_matrix <- str_match(tr_breaks, "sis\\.(\\d+)\\.(\\d+)")
  
  ssvs_breaks   <- names(results$b_ssvs_1$coefs$omega[results$b_ssvs_1$coefs$omega > 0.5])
  source("./BISAM_clean/pip_window_fun.R")
  ssvs_breaks <- pip_window(results$b_ssvs_1, win_size = 1, op = ">=", pip_threshold = 0.50)
  ssvs_breaks <- make_sis_names(ssvs_breaks)
  
  ssvs_breaks_t2 <- pip_window(results$b_ssvs_1, win_size = 2, op = ">=", pip_threshold = 0.75)
  ssvs_breaks_t2 <- make_sis_names(ssvs_breaks_t2)
  
  gets_breaks01 <- rownames(results$gets$`0.01`)[grepl("iis.+|sis.+",rownames(results$gets$`0.01`))]
  all_breaks01  <- str_sort(unique(c(tr_breaks,ssvs_breaks,gets_breaks01,ssvs_breaks_t2)),numeric = T)
  
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
  
  break_comparison01 <- matrix(NA,nrow = length(all_breaks01),ncol = 15)
  colnames(break_comparison01) <- c("true","ssvs","gets",
                                    "tr.ssvs", "tr.gets",
                                    "fp.ssvs","fp.gets",
                                    "fn.ssvs","fn.gets",
                                    "tr.both","fp.both","fn.both","ssvs_1nn_fp","gets_1nn_fp","ssvs_t3_tr")
  rownames(break_comparison01) <- all_breaks01
  
  # which are true/found
  break_comparison01[,1] <- ifelse(all_breaks01%in%tr_breaks,1,0)
  break_comparison01[,2] <- ifelse(all_breaks01%in%ssvs_breaks,1,0)
  break_comparison01[,3] <- ifelse(all_breaks01%in%gets_breaks01,1,0)
  # which are true positive
  break_comparison01[,4] <- (break_comparison01[,2]+break_comparison01[,1])== 2
  break_comparison01[,5] <- (break_comparison01[,3]+break_comparison01[,1])== 2
  # which are false positive
  break_comparison01[,6] <- (break_comparison01[,2]-break_comparison01[,1])== 1
  break_comparison01[,7] <- (break_comparison01[,3]-break_comparison01[,1])== 1
  # which are false negative
  break_comparison01[,8] <- (break_comparison01[,2]-break_comparison01[,1])== -1
  break_comparison01[,9] <- (break_comparison01[,3]-break_comparison01[,1])== -1
  # which tr/fp/fn are in common
  break_comparison01[,10] <- (break_comparison01[,4]+break_comparison01[,5])== 2
  break_comparison01[,11] <- (break_comparison01[,6]+break_comparison01[,7])== 2
  break_comparison01[,12] <- (break_comparison01[,8]+break_comparison01[,9])== 2
  
  break_comparison01[,13] <- (ifelse(all_breaks01%in%tr_breaks_nbh,1,0)+break_comparison01[,6])== 2
  break_comparison01[,14] <- (ifelse(all_breaks01%in%tr_breaks_nbh,1,0)+break_comparison01[,7])== 2
  
  all_t3                  <- ifelse(all_breaks01%in%ssvs_breaks_t2,1,0)
  break_comparison01[,15] <- (all_t3+break_comparison01[,1])== 2  
  
  return(rss01 = break_comparison01)
}
