# ==============================================================================
# REPLICATION CODE FOR:
# "On Spike and Slab Empirical Bayes Multiple Testing"
# by Ismaël Castillo and Étienne Roquain (2019)
# 
# This code replicates the numerical experiments from Section 4 of the paper.
# The goal is to compare different Bayesian Multiple Testing (BMT) procedures
# in terms of FDR control, Power, and Hamming Risk.
# ==============================================================================

# ==============================================================================
# MAIN SIMULATION PROCEDURE
# ==============================================================================
# 
# CONTEXT FROM THE PAPER:
# -----------------------
# The paper considers the Gaussian sequence model (Equation 1):
#   X_i = θ_{0,i} + ε_i,  where ε_i ~ N(0,1) i.i.d.
# 
# The goal is to test simultaneously for i = 1,...,n:
#   H_{0,i}: θ_{0,i} = 0  vs  H_{1,i}: θ_{0,i} ≠ 0
#
# The paper studies "spike and slab" prior distributions (Equation 3):
#   Π_{w,γ} = ((1-w)δ_0 + wG)^{⊗n}
# where:
#   - δ_0 is the Dirac mass at 0 (the "spike")
#   - G is a distribution with density γ (the "slab")
#   - w ∈ (0,1) is a hyperparameter (estimated via marginal maximum likelihood)
#
# KEY QUANTITIES:
# ---------------
# 1. ℓ-value (local FDR, Equation 13): 
#    ℓ(x;w,g) = (1-w)φ(x) / [(1-w)φ(x) + wg(x)]
#    This is the posterior probability that the null is true given X_i = x
#
# 2. q-value (Equation 15):
#    q(x;w,g) = (1-w)Φ̄(|x|) / [(1-w)Φ̄(|x|) + w Ḡ(|x|)]
#    This is the posterior probability that the null is true given |X_i| ≥ |x|
#
# 3. g(x) is defined in Equation 10 as the convolution: g(x) = ∫γ(x-u)φ(u)du
#
# ==============================================================================

procedure = function(prior, alpha, n, sn, nbsimu) {
  # ---------------------------------------------------------------------------
  # PARAMETERS:
  # ---------------------------------------------------------------------------
  # prior  : Choice of slab prior ("cauchy" or "laplace")
  # alpha  : Target FDR level (e.g., 0.2)
  # n      : Total number of hypotheses (dimension)
  # sn     : Number of non-null (signal) hypotheses (sparsity parameter)
  # nbsimu : Number of Monte Carlo simulations
  #
  # The paper assumes θ_0 ∈ ℓ_0[s_n], i.e., at most s_n non-zero coordinates
  # (see Equation 2 for the definition of ℓ_0[s_n])
  # ---------------------------------------------------------------------------
  
  m = n  # m is used as the total number of tests
  
  # ===========================================================================
  # DEFINE THE MARGINAL DENSITY g(x) = ∫γ(x-u)φ(u)du
  # ===========================================================================
  # 
  # This is the convolution of the slab prior γ with the standard normal φ
  # (see Equation 10 in the paper)
  #
  # The function g appears in:
  # - The ℓ-value formula (Equation 13)
  # - The q-value formula via Ḡ (Equation 15-16)
  # ===========================================================================
  
  if (prior == "cauchy") {
    # -------------------------------------------------------------------------
    # QUASI-CAUCHY PRIOR (Equations 20-21)
    # -------------------------------------------------------------------------
    # The paper notes that for the true Cauchy prior, the integral g(x) is not
    # explicit. However, one can use the "quasi-Cauchy" prior (see [31] in refs)
    # which has explicit formulas:
    #
    # γ(x) = (2π)^{-1/2} (1 - |x|Φ̄(x)/φ(x))   [Equation 20]
    # g(x) = (2π)^{-1/2} x^{-2} (1 - e^{-x²/2})  [Equation 21]
    #
    # Note: The quasi-Cauchy satisfies the required conditions (17)-(19):
    # - Lipschitz log-density
    # - Heavy tails (κ=2 in Equation 18)
    # - Bounded x²γ(x)
    # -------------------------------------------------------------------------
    g = function(x) {
      res = 1 / (2 * sqrt(2 * pi))  # Limit as x → 0
      if (x != 0) {
        res = (1 / sqrt(2 * pi)) * (1 - exp(-x^2 / 2)) / x^2
      }
      return(res)
    }
  }
  
  if (prior == "laplace") {
    # -------------------------------------------------------------------------
    # LAPLACE PRIOR (Remark 7 in the paper, Equation 89)
    # -------------------------------------------------------------------------
    # The Laplace prior with parameter a > 0:
    #   γ(x) = γ_a(x) = (a/2) e^{-a|x|}
    #
    # From Remark 7, the convolution g(x) has the explicit form:
    #   g(x) = (a/2)e^{a²/2} [e^{-ax}Φ̄(a-x) + e^{ax}Φ̄(a+x)]
    #
    # The code below uses a = 1 (implied by the formula structure)
    # Note: exp(1/2) = e^{a²/2} when a = 1
    # -------------------------------------------------------------------------
    g = function(x) {
      (exp(1/2) / 2) * (exp(-x) * pnorm(x - 1) + exp(x) * (1 - pnorm(x + 1)))
    }
  }
  
  # ===========================================================================
  # DEFINE Ḡ(u) = ∫_u^∞ g(τ) dτ
  # ===========================================================================
  # This is the tail function of g, needed for q-values (Equation 16)
  # The q-value formula uses Ḡ(|x|) in the denominator
  # ===========================================================================
  Gbar = function(u) {
    integrate(Vectorize(function(tau) g(tau)), u, 40)$value
  }
  
  # ===========================================================================
  # INITIALIZE RESULT MATRICES
  # ===========================================================================
  # 
  # The code compares several methods across different signal strengths:
  # - BH: Benjamini-Hochberg procedure (benchmark)
  # - ℓ-value thresholding at various cutoffs t (Algorithm EBayesL, Eq. 30)
  # - SC procedure: Sun-Cai averaging of ℓ-values (Section 5.2)
  # - q-value thresholding at various cutoffs t (Algorithm EBayesq, Eq. 32)
  #
  # Number of methods = 1 (BH) + 3 × length(seqt)
  #   = 1 + 3 × 5 = 16 methods
  # ===========================================================================
  
  nbmethod = 3 * length(seqt) + 1
  nbmoy = length(seqmoy)
  
  # FDR: False Discovery Rate (Definition in Equation 9)
  #      FDR(θ_0, φ) = E_{θ_0}[#{false discoveries} / (#{total discoveries} ∨ 1)]
  FDR = matrix(0, nbmethod, nbmoy)
  
  # Pow: Power (proportion of true signals detected)
  Pow = matrix(0, nbmethod, nbmoy)
  
  # Risk: Hamming Risk = #{false positives} + #{false negatives}
  #       This is related to FDR+FNR (see Corollary 1 and discussion around Eq. 39)
  Risk = matrix(0, nbmethod, nbmoy)
  
  # ===========================================================================
  # MAIN SIMULATION LOOP OVER SIGNAL STRENGTHS
  # ===========================================================================
  i = 1
  for (moy in seqmoy) {
    # -------------------------------------------------------------------------
    # CONSTRUCT THE TRUE PARAMETER VECTOR θ_0
    # -------------------------------------------------------------------------
    # The paper considers:
    # - Constant alternatives: θ_{0,i} = μ for i ∈ S_0 (support set)
    # - θ_0 ∈ ℓ_0[s_n]: sparse vectors with at most s_n non-zeros
    #
    # Here: first (n - s_n) coordinates are 0, last s_n are equal to 'moy'
    # -------------------------------------------------------------------------
    mu = 0
    if (sn > 0) {
      mu = c(rep(0, m - sn), rep(moy, sn))
    }
    
    # =========================================================================
    # MONTE CARLO LOOP
    # =========================================================================
    for (simu in 1:nbsimu) {
      # -----------------------------------------------------------------------
      # GENERATE DATA FROM THE GAUSSIAN SEQUENCE MODEL (Equation 1)
      # X_i = θ_{0,i} + ε_i, where ε_i ~ N(0,1)
      # -----------------------------------------------------------------------
      X = mu + rnorm(m)
      
      # =======================================================================
      # METHOD 1: BENJAMINI-HOCHBERG (BH) PROCEDURE
      # =======================================================================
      # The BH procedure (Reference [7]) is the benchmark for FDR control.
      # 
      # Algorithm:
      # 1. Compute p-values: p_i = 2Φ̄(|X_i|) (two-sided test)
      # 2. Order p-values: p_{(1)} ≤ p_{(2)} ≤ ... ≤ p_{(n)}
      # 3. Find k* = max{k : p_{(k)} ≤ αk/n}
      # 4. Reject H_{0,i} if p_i ≤ αk*/n
      #
      # The paper mentions (Section 1.1) that Abramovich et al. [2] showed
      # BH achieves minimax adaptive risk properties.
      # =======================================================================
      pvalues = 2 * (1 - pnorm(abs(X)))  # Two-sided p-values
      pvaluesort = sort(pvalues)
      set = which(pvaluesort <= alpha * 1:m / m)  # Find valid indices
      
      if (length(set) > 0) {
        BHth = alpha * max(set) / m  # Threshold
        
        # FDR: proportion of false discoveries among all discoveries
        FDR[1, i] = FDR[1, i] + 
          (sum((mu == 0) & (pvalues <= BHth)) / sum(pvalues <= BHth)) / nbsimu
        
        # Power: proportion of true signals discovered
        if (sum(mu > 0)) {
          Pow[1, i] = Pow[1, i] + 
            (sum((mu > 0) & (pvalues <= BHth)) / sum(mu > 0)) / nbsimu
        }
        
        # Hamming Risk: false positives + false negatives
        Risk[1, i] = Risk[1, i] + 
          (sum((mu == 0) & (pvalues <= BHth)) + sum((mu > 0) & (pvalues > BHth))) / nbsimu
      }
      
      # =======================================================================
      # COMPUTE ℓ-VALUES USING EMPIRICAL BAYES
      # =======================================================================
      # 
      # The function ebayesthresh from the EbayesThresh package [31] computes:
      # 1. The marginal maximum likelihood estimate ŵ (Equation 29)
      # 2. The posterior distribution
      #
      # The ℓ-value is computed from Equation 13:
      #   ℓ_i(X) = ℓ(X_i; ŵ, g) = (1-ŵ)φ(X_i) / [(1-ŵ)φ(X_i) + ŵg(X_i)]
      #
      # Key insight from the paper (Section 2.6):
      # The MMLE ŵ is found by maximizing the marginal likelihood L(w),
      # or equivalently finding where the score function S(w) = 0 (Equation 28)
      # =======================================================================
      JS = ebayesthresh(X, verbose = TRUE, prior)
      
      # Compute ℓ-values using the formula from Equation 13
      lvalues = (1 - JS$w) * dnorm(X) / 
        ((1 - JS$w) * dnorm(X) + JS$w * sapply(X, FUN = g))
      
      # =======================================================================
      # METHODS 2-6: ℓ-VALUE THRESHOLDING (EBayesL)
      # =======================================================================
      # Algorithm EBayesL (page 12):
      # Reject H_{0,i} if ℓ̂_i(X) ≤ t
      #
      # THEOREM 1 guarantees (Equation 31):
      #   sup_{θ_0 ∈ ℓ_0[s_n]} FDR(θ_0, φ^{ℓ-val}) ≤ C × log(log n) / (log n)^{κ/2}
      # 
      # This bound goes to 0 as n → ∞, meaning the ℓ-value procedure is
      # conservative - it doesn't "spend" the full FDR budget α.
      # =======================================================================
      method = 2  # Counter for methods (BH is method 1)
      for (t in seqt) {
        if (sum(lvalues <= t) > 0) {
          FDR[method, i] = FDR[method, i] + 
            (sum((mu == 0) & (lvalues <= t)) / sum(lvalues <= t)) / nbsimu
          
          if (sum(mu > 0)) {
            Pow[method, i] = Pow[method, i] + 
              (sum((mu > 0) & (lvalues <= t)) / sum(mu > 0)) / nbsimu
          }
        }
        Risk[method, i] = Risk[method, i] + 
          (sum((mu == 0) & (lvalues <= t)) + sum((mu > 0) & (lvalues > t))) / nbsimu
        
        method = method + 1
      }
      
      # =======================================================================
      # METHODS 7-11: SUN-CAI (SC) PROCEDURE
      # =======================================================================
      # The SC procedure (Section 5.2) averages ℓ-values:
      # 
      # Algorithm:
      # 1. Order ℓ-values: ℓ̂_{(1)} ≤ ℓ̂_{(2)} ≤ ... ≤ ℓ̂_{(n)}
      # 2. Find k̂ = max{k : (1/k) Σ_{k'=1}^k ℓ̂_{(k')} ≤ t}
      # 3. Reject the k̂ smallest ℓ-values
      #
      # The paper notes (Section 16) that SC shows FDR inflation compared to
      # q-values, especially when s_n/n is not small. The convergence to
      # target level t is slow (logarithmic in n/s_n).
      # =======================================================================
      for (t in seqt) {
        lvaluesort = sort(lvalues)
        cummean = cumsum(lvaluesort) / 1:m  # Running average of ordered ℓ-values
        set = which(cummean <= t)
        selected = rep(FALSE, m)
        
        if (length(set) > 0) {
          kSC = max(set)
          selected[(order(lvalues)[1:kSC])] = TRUE  # Select kSC smallest ℓ-values
          
          FDR[method, i] = FDR[method, i] + 
            (sum((mu == 0) & selected) / sum(selected)) / nbsimu
          
          if (sum(mu > 0)) {
            Pow[method, i] = Pow[method, i] + 
              (sum((mu > 0) & selected) / sum(mu > 0)) / nbsimu
          }
        }
        Risk[method, i] = Risk[method, i] + 
          (sum((mu == 0) & selected) + sum((mu > 0) & (!selected))) / nbsimu
        
        method = method + 1
      }
      
      # =======================================================================
      # COMPUTE q-VALUES
      # =======================================================================
      # The q-value (Equation 15):
      #   q(x; w, g) = (1-w)Φ̄(|x|) / [(1-w)Φ̄(|x|) + w Ḡ(|x|)]
      #
      # where Φ̄(s) = P(|ε| ≥ s) = 2Φ̄(s) is the tail probability
      # and Ḡ(s) = ∫_s^∞ g(x)dx
      #
      # Key property (Lemma 10): q_i(X) ≤ ℓ_i(X), so q-values are more liberal
      # =======================================================================
      qvalues = (1 - JS$w) * (1 - pnorm(abs(X))) / 
        ((1 - JS$w) * (1 - pnorm(abs(X))) + JS$w * sapply(abs(X), FUN = Gbar))
      
      # =======================================================================
      # METHODS 12-16: q-VALUE THRESHOLDING (EBayesq)
      # =======================================================================
      # Algorithm EBayesq (page 13):
      # Reject H_{0,i} if q̂_i(X) ≤ t
      #
      # THEOREM 2 guarantees (Equation 35):
      #   sup_{θ_0 ∈ ℓ_0[s_n]} FDR(θ_0, φ^{q-val}) ≤ C × t × log(1/t)
      #
      # THEOREM 3 shows that for strong signals (θ_0 ∈ L_0[s_n], Eq. 37):
      #   lim_n FDR(θ_0, φ^{q-val}) = t  (exactly!)
      #
      # This is the KEY RESULT: q-values achieve exact FDR control at level t
      # for large signals, unlike ℓ-values which are conservative.
      # =======================================================================
      for (t in seqt) {
        if (sum(qvalues <= t) > 0) {
          FDR[method, i] = FDR[method, i] + 
            (sum((mu == 0) & (qvalues <= t)) / sum(qvalues <= t)) / nbsimu
          
          if (sum(mu > 0)) {
            Pow[method, i] = Pow[method, i] + 
              (sum((mu > 0) & (qvalues <= t)) / sum(mu > 0)) / nbsimu
          }
        }
        Risk[method, i] = Risk[method, i] + 
          (sum((mu == 0) & (qvalues <= t)) + sum((mu > 0) & (qvalues > t))) / nbsimu
        
        method = method + 1
      }
    }
    i = i + 1
  }
  
  return(list(FDR = FDR, Pow = Pow, RiskHam = Risk))
}


# ==============================================================================
# FILE NAMING UTILITY
# ==============================================================================
makefilename = function(prior, alpha, n, sn, nbsimu, Risk) {
  filename = paste(Risk, "EBmovemoy", sep = "")
  filename = paste(filename, "_prior", prior, "_alpha", alpha, "_n", n,
                   "_snsurn", sn/n, "_nbsimu", nbsimu, sep = "")
  filename = gsub("\\.", "_", filename)
  return(filename)
}


# ==============================================================================
# MAIN FUNCTION FOR PARALLEL EXECUTION
# ==============================================================================
mainfunction = function(param) {
  prior = param[1]
  alpha = as.numeric(param[2])
  n = as.numeric(param[3])
  sn = ceiling(as.numeric(param[4]) * n)  # Convert proportion to count
  
  res = procedure(prior, alpha, n, sn, nbsimu)
  filename = makefilename(prior, alpha, n, sn, nbsimu, "")
  save(res, file = paste("./Replication_Castillo_Roquain_FDR/", filename, sep = ""))
  
  print("ok")
}


# ==============================================================================
# SIMULATION PARAMETERS
# ==============================================================================

# Signal strength sequence (μ values for non-null hypotheses)
# Paper Figure 1 uses μ on x-axis from 0 to 10
seqmoy = c(seq(0.01, 0.901, 0.1), 1:10)

# Threshold values for ℓ-value, q-value, and SC procedures
# Paper tests t ∈ {0.05, 0.1, 0.2} (see Figure 1)
seqt = c(0.05, 0.1, 0.2, 0.5, 0.75)

# ==============================================================================
# PARAMETERS FOR PARALLEL EXECUTION
# ==============================================================================

# Choice of slab prior (Section 2.3)
# - "cauchy": quasi-Cauchy prior (κ = 2 in Theorem 1)
# - "laplace": Laplace prior (κ = 1 in Theorem 1)
seqprior = c("cauchy", "laplace")

# Target FDR level
# Paper uses α = 0.2 for main simulations
seqalpha = c(0.2)

# Dimension n = 10^5 (paper uses n = 10^4 in Section 4)
n = 10^5
seqn = c(n)

# Sparsity levels s_n/n
# These determine the proportion of non-null hypotheses
# Paper requires s_n ≤ n^υ for υ ∈ (0,1) (polynomial sparsity)
# Values tested: 1%, 3%, 6%, 10%, 20%, 30%, 40%, 50%, 60%
seqsnsurn = c(0.01, 0.03, 0.06, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6)

# Number of Monte Carlo replications
# Paper uses 2000 replications (Section 4)
nbsimu = 50

# ==============================================================================
# CREATE PARAMETER GRID FOR PARALLEL EXECUTION
# ==============================================================================
parameter = list()
for (prior in seqprior) {
  for (alpha in seqalpha) {
    for (n in seqn) {
      for (snsurn in seqsnsurn) {
        parameter = c(parameter, list(c(prior, alpha, n, snsurn)))
      }
    }
  }
}

# ==============================================================================
# LOAD REQUIRED PACKAGES
# ==============================================================================
library(EbayesThresh)  # For empirical Bayes thresholding [31]
library(parallel)       # For parallel computing
library(foreach)        # For parallel loops
library(doParallel)     # Parallel backend

# ==============================================================================
# RUN SIMULATIONS IN PARALLEL
# ==============================================================================
mclapply(parameter, mainfunction, mc.cores = 18)


# ==============================================================================
# PLOTTING CODE - REPLICATES FIGURE 1 FROM THE PAPER
# ==============================================================================
# 
# Figure 1 in the paper shows:
# - FDR as function of signal strength μ
# - Different rows: different sparsity levels s_n/n
# - Different columns: quasi-Cauchy vs Laplace
# - Horizontal line at α shows the target FDR level
#
# Legend:
# - BH: Benjamini-Hochberg (dashed black, lty=2)
# - lval: ℓ-value thresholding (solid, lty=1)
# - SC: Sun-Cai procedure (dotted-dashed, lty=4)
# - qval: q-value thresholding (dashed, lty=3)
# ==============================================================================

for (prior in seqprior) {
  for (alpha in seqalpha) {
    for (n in seqn) {
      for (snsurn in seqsnsurn) {
        sn = ceiling(snsurn * n)
        filename = makefilename(prior, alpha, n, sn, nbsimu, "")
        load(filename)
        
        # =====================================================================
        # PLOT FDR
        # =====================================================================
        # Key observations from Figure 1 and Theorems 1-3:
        # - ℓ-values: FDR → 0 (conservative, Theorem 1)
        # - q-values: FDR → t for large signals (exact control, Theorem 3)
        # - SC: Shows FDR inflation, especially for weaker sparsity
        # =====================================================================
        resFDR = res$FDR
        filename = makefilename(prior, alpha, n, sn, nbsimu, "FDR")
        filename = paste("./Replication_Castillo_Roquain_FDR/", filename, ".pdf", sep = "")
        pdf(filename)
        
        # Plot BH (method 1)
        plot(seqmoy, resFDR[1,], lty = 2, col = 1, xlab = "", ylab = "",
             cex.axis = 2, type = "l", lwd = 3, ylim = c(0, 1),
             xlim = c(0, max(seqmoy) + 2))
        
        # Plot ℓ-value methods (methods 2 to length(seqt)+1)
        for (meth in 2:(length(seqt) + 1)) {
          lines(seqmoy, resFDR[meth,], lty = 1, col = meth, lwd = 3)
        }
        
        # Plot SC methods (methods length(seqt)+2 to 2*length(seqt)+1)
        for (meth in (length(seqt) + 2):(2 * length(seqt) + 1)) {
          lines(seqmoy, resFDR[meth,], lty = 4, col = meth - length(seqt), lwd = 3)
        }
        
        # Plot q-value methods (methods 2*length(seqt)+2 to 3*length(seqt)+1)
        for (meth in (2 * length(seqt) + 2):(3 * length(seqt) + 1)) {
          lines(seqmoy, resFDR[meth,], lty = 3, col = meth - 2 * length(seqt), lwd = 3)
        }
        
        # Horizontal line at target level α
        abline(h = alpha)
        
        # Legend
        lg = c("BH", paste("lval", seqt), paste("SC", seqt), paste("qval", seqt))
        legend("topright",
               col = c(1:(length(seqt) + 1), 2:(length(seqt) + 1), 2:(length(seqt) + 1)),
               lty = c(2, rep(1, length(seqt)), rep(4, length(seqt)), rep(3, length(seqt))),
               legend = lg, bg = "white", cex = 0.7)
        dev.off()
        
        # =====================================================================
        # PLOT POWER
        # =====================================================================
        # Power = proportion of true signals detected
        # Theorem 4 shows: for large signals, FNR → 0, hence Power → 1
        # =====================================================================
        resPow = res$Pow
        filename = makefilename(prior, alpha, n, sn, nbsimu, "Pow")
        filename = paste("./Replication_Castillo_Roquain_FDR/", filename, ".pdf", sep = "")
        pdf(filename)
        
        plot(seqmoy, resPow[1,], lty = 2, col = 1, xlab = "", ylab = "",
             cex.axis = 2, type = "l", lwd = 3, ylim = c(0, 1),
             xlim = c(0, max(seqmoy) + 2))
        
        for (meth in 2:(length(seqt) + 1)) {
          lines(seqmoy, resPow[meth,], lty = 1, col = meth, lwd = 3)
        }
        for (meth in (length(seqt) + 2):(2 * length(seqt) + 1)) {
          lines(seqmoy, resPow[meth,], lty = 4, col = meth - length(seqt), lwd = 3)
        }
        for (meth in (2 * length(seqt) + 2):(3 * length(seqt) + 1)) {
          lines(seqmoy, resPow[meth,], lty = 3, col = meth - 2 * length(seqt), lwd = 3)
        }
        
        abline(h = alpha)
        lg = c("BH", paste("lval", seqt), paste("SC", seqt), paste("qval", seqt))
        legend("topright",
               col = c(1:(length(seqt) + 1), 2:(length(seqt) + 1), 2:(length(seqt) + 1)),
               lty = c(2, rep(1, length(seqt)), rep(4, length(seqt)), rep(3, length(seqt))),
               legend = lg, bg = "white", cex = 0.7)
        dev.off()
        
        # =====================================================================
        # PLOT HAMMING RISK
        # =====================================================================
        # Hamming Risk = #{false positives} + #{false negatives}
        # Related to the classification risk R(θ_0, φ) = FDR + FNR
        # Corollary 1 shows:
        # - For ℓ-values: FDR + FNR → 0 over L_0[s_n;a], a > 1
        # - For q-values: FDR → t, FNR → 0, so R → t
        # Proposition 2: Classification impossible for a < 1
        # =====================================================================
        resRiskHam = res$RiskHam
        filename = makefilename(prior, alpha, n, sn, nbsimu, "RiskHam")
        filename = paste("./Replication_Castillo_Roquain_FDR/", filename, ".pdf", sep = "")
        pdf(filename)
        
        plot(seqmoy, resRiskHam[1,], lty = 2, col = 1, xlab = "", ylab = "",
             cex.axis = 2, type = "l", lwd = 3, ylim = c(0, max(resRiskHam)),
             xlim = c(0, max(seqmoy) + 2))
        
        for (meth in 2:(length(seqt) + 1)) {
          lines(seqmoy, resRiskHam[meth,], lty = 1, col = meth, lwd = 3)
        }
        for (meth in (length(seqt) + 2):(2 * length(seqt) + 1)) {
          lines(seqmoy, resRiskHam[meth,], lty = 4, col = meth - length(seqt), lwd = 3)
        }
        for (meth in (2 * length(seqt) + 2):(3 * length(seqt) + 1)) {
          lines(seqmoy, resRiskHam[meth,], lty = 3, col = meth - 2 * length(seqt), lwd = 3)
        }
        
        abline(h = alpha)
        lg = c("BH", paste("lval", seqt), paste("SC", seqt), paste("qval", seqt))
        legend("topright",
               col = c(1:(length(seqt) + 1), 2:(length(seqt) + 1), 2:(length(seqt) + 1)),
               lty = c(2, rep(1, length(seqt)), rep(4, length(seqt)), rep(3, length(seqt))),
               legend = lg, bg = "white", cex = 0.7)
        dev.off()
      }
    }
  }
}

# ==============================================================================
# END OF REPLICATION CODE
# ==============================================================================