#' bisam: Bayesian indicator saturated models
#'
#' Main function: \code{\link{estimate_bisam}}. See also \code{\link{pip_window}},
#' \code{\link{contr_sim_breaks}} and \code{\link{priorp2g}}.
#'
#' The compiled core (break-indicator Gibbs sampler, iMOM posterior mode, non-local posterior
#' draws of break sizes) is derived from mombf 3.5.4 (D. Rossell et al., GPL >= 2); see
#' \code{system.file("IMPLEMENTATION.md", package = "bisam")} for its differences to mombf
#' and how it was validated.
#'
#' @useDynLib bisam, .registration = TRUE
#' @importFrom Rcpp evalCpp
#' @import Matrix
#' @importFrom methods as is
#' @importFrom mvtnorm rmvnorm
#' @importFrom stats dnorm lm.fit nlminb pgamma pnorm qgamma quantile rbeta rbinom rgamma rnorm
#'   runif
#' @importFrom utils setTxtProgressBar txtProgressBar
#' @keywords internal
"_PACKAGE"
