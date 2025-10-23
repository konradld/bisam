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

# # Usage:
# # Find (shape,scale) such that P(X <= 1) = 0.1 and P(X <= 2) = 0.9
# library(nleqslv)
# library(invgamma)
# 
# params <- inv_gamma_params(shape = 3, x = 3, p = 0.9)
# 
# N <- 1e6
# test <- 1/rgamma(N, shape = params$shape, scale = params$scale)
# sum(test > 3) / N
# 
# 
# params <- inv_gamma_params_dual(x1 = 1, p1 = 0.1, x2 = 2, p2 = 0.9)
# params
