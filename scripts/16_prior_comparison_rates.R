# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#        What the non-local prior buys: rates of evidence accumulation
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#
# Restores, as an empirical result, the claim the paper states at
# final_revised_paper.tex:268 and then CUT at lines 270-276:
#
#   "Under a pMOM prior of order k, BF_T(1|0) = O_p(T^{-k-1/2}), while under the
#    piMOM prior, log BF_T(1|0) / T^{k/(k+1)} -> c < 0, implying a near-exponential
#    rate of evidence accumulation against spurious breaks.
#    Both priors preserve exponential learning rates when breaks are truly present."
#
# The paper asserts this by citation (Johnson & Rossell 2010) and never shows it:
# Figure 1 plots two tau settings of the SAME iMOM density, and Figures 2-4 fix
# bisam-imom and compare against GETS/ALASSO, never against another prior.
#
# WHY THE BAYES-FACTOR LEVEL. BF_T(1|0) is exactly the object those lines make
# claims about. Computing it directly with mombf::nlpMarginal keeps the result
# clean of the two things that confound the full sampler: the Bernoulli(0.5)
# inclusion prior (no multiplicity control) and the ALA/getthinit patch. It also
# runs in minutes instead of a cluster sweep.
#
# DESIGN. One cross-sectional unit of the paper's SIS design: candidate steps at
# dates 3..(T-1), i.e. the columns of the panel's Z block for a single unit (see
# R/contr_sim_breaks_fun.R). This is the right unit of analysis because the
# sampler factorises over units (paper, line 604), so the per-unit problem is
# "T-3 candidates from T observations".
#
# Output: output/simulation/prior_rate_comparison.RDS
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

suppressPackageStartupMessages(library(mombf))

OUT_FILE <- "./output/simulation/prior_rate_comparison.RDS"
TS       <- c(20, 40, 80, 160, 320, 640)
REPS     <- 200
STEP_SD  <- 3      # true break size in error SDs, as in the paper's headline design
NI       <- 10     # units, for the multiplicity burden log K = log(Ni*(T-3))

# ------------------------------------------------------------------------------
# Prior calibration
# ------------------------------------------------------------------------------
# The paper calibrates tau by inverting P(|gamma| <= gamma_bar | tau) = alpha with
# gamma_bar = sigma and alpha = 0.01 (tex:308, tex:683; Table A.2). Apply that SAME
# principle to every prior so the comparison is like-for-like and follows the
# paper's own stated methodology rather than an arbitrary default.
#
#   iMOM : priorp2g(0.01, 1, nu=1, prior="iMom")      = 3.3174  (the paper's tau)
#   pMOM : priorp2g(0.01, 1, nu=1, prior="normalMom") = 8.7084
#   normal slab N(0, tau*sigma^2): P(|g|<=sigma) = 2*pnorm(1/sqrt(tau)) - 1 = 0.01
#          => tau = (1/qnorm(0.505))^2, a very diffuse slab. That is not a strawman:
#          it is what "only 1% prior mass within one SD" forces a local prior to be.
#
# RATE IS CALIBRATION-INVARIANT. tau moves the intercept of log BF, not its slope
# in T. The Zellner arm is therefore also run at the conventional unit-information
# g = n, and both settings must land on the same -1/2 slope. That is the figure's
# built-in answer to "you just picked a bad tau".
TAU_IMOM <- priorp2g(0.01, 1, nu = 1, prior = "iMom")
TAU_MOM  <- priorp2g(0.01, 1, nu = 1, prior = "normalMom")
TAU_NORM <- (1 / qnorm(0.505))^2

prior_arms <- list(
  list(key = "zellner_g_n",   fam = "local", lab = "Local: Zellner (g = n)",
       f = function(n) zellnerprior(tau = n)),
  list(key = "normalid_cal",  fam = "local", lab = "Local: normal slab (calibrated)",
       f = function(n) normalidprior(tau = TAU_NORM)),
  list(key = "pmom",          fam = "mom",   lab = "pMOM (k = 1)",
       f = function(n) momprior(tau = TAU_MOM)),
  list(key = "pimom",         fam = "imom",  lab = "piMOM (k = 1)",
       f = function(n) imomprior(tau = TAU_IMOM))
)

# ------------------------------------------------------------------------------
# One unit's SIS design: step indicators 1{t >= s}, s = 3..(T-1)
# ------------------------------------------------------------------------------
sis_block <- function(Ti) {
  1 * lower.tri(matrix(1, Ti, Ti), diag = TRUE)[, 3:(Ti - 1), drop = FALSE]
}

# log BF for including candidate j against the null model, on data y.
log_bf <- function(y, x, pr) {
  m0 <- nlpMarginal(sel = integer(0), y = y, x = x, priorCoef = pr,
                    logscale = TRUE, method = "Laplace")
  m1 <- nlpMarginal(sel = 1, y = y, x = x, priorCoef = pr,
                    logscale = TRUE, method = "Laplace")
  m1 - m0
}

# ------------------------------------------------------------------------------
# Simulate
# ------------------------------------------------------------------------------
# H0 arm: pure-noise unit, evidence AGAINST a randomly chosen spurious candidate.
# H1 arm: a true STEP_SD break at a random date, evidence FOR the true candidate.
run_all <- function() {
  res <- expand.grid(T = TS, arm = c("H0", "H1"), key = sapply(prior_arms, `[[`, "key"),
                     stringsAsFactors = FALSE)
  res$mean <- NA_real_; res$se <- NA_real_
  t0 <- Sys.time()
  for (Ti in TS) {
    L <- sis_block(Ti)
    K <- ncol(L)
    for (a in prior_arms) {
      pr <- a$f(Ti)
      v0 <- v1 <- numeric(REPS)
      for (r in seq_len(REPS)) {
        j <- sample(K, 1)
        # --- H0: null is TRUE at candidate j
        y0 <- as.numeric(scale(rnorm(Ti), scale = FALSE))
        v0[r] <- log_bf(y0, L[, j, drop = FALSE], pr)
        # --- H1: a real break of STEP_SD sigma at candidate j
        y1 <- as.numeric(scale(STEP_SD * L[, j] + rnorm(Ti), scale = FALSE))
        v1[r] <- log_bf(y1, L[, j, drop = FALSE], pr)
      }
      i0 <- res$T == Ti & res$arm == "H0" & res$key == a$key
      i1 <- res$T == Ti & res$arm == "H1" & res$key == a$key
      res$mean[i0] <- mean(v0); res$se[i0] <- sd(v0) / sqrt(REPS)
      res$mean[i1] <- mean(v1); res$se[i1] <- sd(v1) / sqrt(REPS)
    }
    cat(sprintf("  T = %4d done (%.1f min)\n", Ti,
                as.numeric(difftime(Sys.time(), t0, units = "mins"))))
  }
  res
}

# ------------------------------------------------------------------------------
# Rate fits — the quantitative claim, reported alongside the figure
# ------------------------------------------------------------------------------
# Theory (Johnson & Rossell 2010, eq. 11; the paper's cut lines 270-276):
#   local      log BF ~ -0.5 * log T          (polynomial)
#   pMOM(k)    log BF ~ -(k + 0.5) * log T    (polynomial)
#   piMOM(k)   log BF / T^{k/(k+1)} -> c < 0  (root-exponential at k = 1)
rate_table <- function(res) {
  out <- do.call(rbind, lapply(prior_arms, function(a) {
    d <- res[res$key == a$key & res$arm == "H0", ]
    d <- d[order(d$T), ]
    f_log <- lm(d$mean ~ log(d$T))
    f_sqrt <- lm(d$mean ~ sqrt(d$T))
    h <- res[res$key == a$key & res$arm == "H1", ]
    h <- h[order(h$T), ]
    f_lin <- lm(h$mean ~ h$T)
    data.frame(
      prior        = a$lab,
      slope_logT   = coef(f_log)[2],  r2_logT  = summary(f_log)$r.squared,
      slope_sqrtT  = coef(f_sqrt)[2], r2_sqrtT = summary(f_sqrt)$r.squared,
      favours      = ifelse(AIC(f_log) < AIC(f_sqrt), "log T (polynomial)",
                            "sqrt(T) (root-exponential)"),
      H1_slope_T   = coef(f_lin)[2],  H1_r2 = summary(f_lin)$r.squared,
      stringsAsFactors = FALSE)
  }))
  rownames(out) <- NULL
  out
}

# ------------------------------------------------------------------------------
if (sys.nframe() == 0L) {
  set.seed(20260716)
  cat(sprintf("tau: iMOM %.4f | pMOM %.4f | normal slab %.1f\n",
              TAU_IMOM, TAU_MOM, TAU_NORM))
  cat(sprintf("simulating %d reps x %d T-values x %d priors x 2 arms ...\n",
              REPS, length(TS), length(prior_arms)))
  res <- run_all()
  tab <- rate_table(res)
  saveRDS(list(res = res, rates = tab, reps = REPS, step_sd = STEP_SD, Ni = NI,
               tau = c(imom = TAU_IMOM, mom = TAU_MOM, normal = TAU_NORM)),
          OUT_FILE)
  cat("\n--- Evidence AGAINST a spurious break (H0 true) ---\n")
  print(format(tab[, c("prior", "slope_logT", "r2_logT", "slope_sqrtT",
                       "r2_sqrtT", "favours")], digits = 3), row.names = FALSE)
  cat("\n--- Evidence FOR a true break (H1 true): slope in T ---\n")
  print(format(tab[, c("prior", "H1_slope_T", "H1_r2")], digits = 3),
        row.names = FALSE)
  cat("\nwrote", OUT_FILE, "\n")
}
