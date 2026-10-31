### Prior helpers (verbatim from mombf 3.5.4, R/dmom.R)
###
### Product MOM / iMOM cumulative distribution and quantile functions, and the
### prior-elicitation helpers priorp2g() / mode2g() used by the BISAM scripts to
### choose the prior dispersion tau.

pmom <- function(q,V1=1,tau=1) {

  z <- .5-(pnorm(abs(q)/sqrt(V1*tau)) - abs(q)/sqrt(2*pi*V1*tau) * exp(-.5*q^2/(tau*V1)) - .5)
  return(z*(q<=0)+(1-z)*(q>0))
}


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


pimom <- function(q,V1=1,tau=1,nu=1) {

  z <- (.5*pgamma(tau*V1/q^2,nu/2,1))
  return(z*(q<=0)+(1-z)*(q>0))
}


qimom <- function(p, V1=1, tau=1, nu=1) {

  ans <- double(length(p))
  ans[p<=.5] <- -sqrt(tau*V1/qgamma(2*p[p<=.5], nu/2, 1))
  ans[p>.5] <- sqrt(tau*V1/qgamma(2*(1-p[p>.5]), nu/2, 1))
  ans
}


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
