### Prior helpers (verbatim from mombf 3.5.4, R/dmom.R; D. Rossell et al., GPL >= 2)

#' MOM and iMOM prior distribution functions
#'
#' Cumulative distribution and quantile functions of the univariate MOM prior
#' \eqn{\theta^2/(\tau V_1) N(\theta; 0, \tau V_1)} and of the iMOM prior with dispersion
#' \code{tau} (from mombf 3.5.4).
#'
#' @param q Vector of quantiles.
#' @param p Vector of probabilities.
#' @param V1 Scale (e.g. the residual variance).
#' @param tau Prior dispersion.
#' @param nu Degrees of freedom of the iMOM prior.
#' @return \code{pmom}, \code{pimom}: probabilities; \code{qmom}, \code{qimom}: quantiles.
#' @examples
#' pimom(-1, tau = 0.133)
#' qmom(0.05, tau = 0.348)
#' @name nlp_distributions
NULL

#' @rdname nlp_distributions
#' @export
pmom <- function(q,V1=1,tau=1) {

  z <- .5-(pnorm(abs(q)/sqrt(V1*tau)) - abs(q)/sqrt(2*pi*V1*tau) * exp(-.5*q^2/(tau*V1)) - .5)
  return(z*(q<=0)+(1-z)*(q>0))
}


#' @rdname nlp_distributions
#' @export
qmom <- function(p, V1=1, tau=1) {

  e <- function(q) { return((pmom(q, V1, tau)-pneg)^2) }
  ans <- double(length(p))
  for (i in 1:length(p)) {
    pneg <- ifelse(p[i]<=.5, p[i], 1-p[i])
    ans[i] <- nlminb(start=-1, objective=e)$par
    if (p[i]>.5) {
      ans[i] <- -ans[i]
    }
  }
  ans
}


#' @rdname nlp_distributions
#' @export
pimom <- function(q,V1=1,tau=1,nu=1) {

  z <- (.5*pgamma(tau*V1/q^2,nu/2,1))
  return(z*(q<=0)+(1-z)*(q>0))
}


#' @rdname nlp_distributions
#' @export
qimom <- function(p, V1=1, tau=1, nu=1) {

  ans <- double(length(p))
  ans[p<=.5] <- -sqrt(tau*V1/qgamma(2*p[p<=.5], nu/2, 1))
  ans[p>.5] <- sqrt(tau*V1/qgamma(2*(1-p[p>.5]), nu/2, 1))
  ans
}


#' Prior elicitation for MOM / iMOM priors
#'
#' \code{priorp2g} returns the prior dispersion \code{tau} such that the prior probability of
#' \eqn{|\theta| < q} equals \code{priorp}; \code{mode2g} returns the \code{tau} giving a prior
#' mode at \code{prior.mode} (from mombf 3.5.4).
#'
#' @param priorp Prior probability of \eqn{|\theta| < q}.
#' @param q Threshold (e.g. a negligible effect size relative to the residual sd).
#' @param prior.mode Desired prior mode.
#' @param nu Degrees of freedom (iMOM, tMOM).
#' @param dim Dimension (for \code{mode2g}).
#' @param prior \code{"iMom"}, \code{"normalMom"} or \code{"tMom"}.
#' @return The prior dispersion \code{tau}.
#' @examples
#' priorp2g(0.01, 1, prior = "iMom")       # tau used in the BISAM simulations
#' priorp2g(0.01, 1, prior = "normalMom")
#' @name prior_elicitation
NULL

#' @rdname prior_elicitation
#' @export
mode2g <- function(prior.mode,
                   prior=c('iMom', 'normalMom', 'tMom'),
                   nu=1,
                   dim=1) {
  prior <- match.arg(prior)
  switch(EXPR=prior,
         normalMom = {
           prior.mode / 2
         },
         iMom = {
           prior.mode * (nu+dim) / 2
         },
         tMom = {
           if (nu<3) {
             stop('tMom prior must have nu>2 degrees of freedom')
           }
           prior.mode * (nu-2+dim) / (2*nu)
         })
}


#' @rdname prior_elicitation
#' @export
priorp2g <- function(priorp,
                     q,
                     nu=1,
                     prior=c('iMom', 'normalMom','tMom')) {
  prior <- match.arg(prior)
  switch(EXPR=prior,
         normalMom = {
           e <- function(logg) {
             return((1-2*pmom(-abs(q), tau=exp(logg)) - priorp[i])^2)
           }

           ans <- double(length(priorp))
           for (i in 1:length(priorp)) {
             ans[i] <- exp(nlminb(start=0, objective=e)$par)
           }
           ans
         },
         tMom = {
           stop("prior=='tMom' is not currently implemented")
         },
         iMom = {
           p <- (1-priorp)/2
           qgamma(2*p,nu/2,1)*q^2
         })
}
