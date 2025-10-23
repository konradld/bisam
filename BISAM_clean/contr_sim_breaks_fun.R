# _________________ Controlled Simulation ______________----
contr_sim_breaks <- function(
    n,
    t,
    nx,
    iis=TRUE,
    sis=TRUE,
    pos.outl, # count in total position of panel setup
    pos.step, # count in total position of panel setup
    const = FALSE,
    ife = FALSE,
    tfe = FALSE,
    outl.mean,
    step.mean,
    error.sd){
  
  require(mvtnorm)
  
  # Create individual and time indicators
  n_<- c(t(matrix(rep(1:n,t), ncol = rep(t,n)))) # unit index
  t_<- c(matrix(rep(1:t,n), ncol = n)) # time index
  
  # Create obs. matrices and corresponding betas
  if(nx > 0){
    X_ <- rmvnorm(n*t,rep(0,nx),diag(nx))
    X  <- X_ # original version will be required later

    b_ <- sample(x = -10:10,size = nx,replace = T)
    b  <- b_
  }else{
    X_ <- NULL
    X  <- NULL
    b_ <- NULL
    b  <- NULL
  }

  if(const) { # Intercept
    X <- cbind(X, "const" = rep(1,n*t))
    b0<- sample(x = -10:10,1)
    b <- c(b,b0)
  }
  if(ife) { # Individual fixed effects
    IFE <- kronecker(diag(n), matrix(1, t))
    if(const) {IFE <- IFE[,-1]}
    X <- cbind(X, IFE)
    bife<-sample(x = -10:10,size = ncol(IFE),replace = T)
    b <- c(b,bife)
  }
  if(tfe) { # Time fixed effects
    TFE <- kronecker(matrix(1, n), diag(t))
    if(const) {TFE <- TFE[,-1]}
    X <- cbind(X, TFE)
    btfe<-sample(x = -10:10,size = ncol(TFE),replace = T)
    b <- c(b,btfe)
  }
  if(tfe & ife & !const) {
    warning("Both time and unit fixed effects used.\n Dropping first indiv. FE dummy to avoid perfect colinearity")
    X <- X[, -(ncol(X_)+1)]
    b <- b[-(ncol(X_)+1)]
  }
  # IIS
  I <- diag(n*t)
  tr_a <- rep(0,n*t)
  tr_a[pos.outl] <- 1
  a <- tr_a*outl.mean
  
  # SIS
  Z <- kronecker(diag(n),lower.tri(matrix(1,nrow=t,ncol=t),diag = T)[,-c(1,t*iis)])
  tr_g <- rep(0,n*(t-1-iis))
  tr_g[(pos.step%/%t)*(t-1-iis)+pos.step%%t-1] <- 1
  g <- tr_g*step.mean
  
  e <- rnorm(n*t,0,error.sd)
  
  if(is.null(X)){
    X = matrix(0,nrow=n*t)
    b = 0
  }
  
  y <- X%*%b + I%*%a + Z%*%g + e
  
  tr.ind <- cbind(which(t(matrix(a,ncol=t,byrow = T))!=0,arr.ind = T),a[a!=0])
  tr.stp <- cbind(which(t(cbind(rep(0,n),matrix(g,ncol=(t-1-iis),byrow = T)))!=0,arr.ind = T),g[g!=0])
  if(length(tr.ind)!=0){
    rownames(tr.ind) <- paste('iis',tr.ind[,2],tr.ind[,1],sep = '.')
  }
  if(length(tr.stp)!=0){
    rownames(tr.stp) <- paste('sis',tr.stp[,2],tr.stp[,1],sep = '.')
  }   
  tr.ind <- cbind(as.matrix(tr.ind),'index'=(tr.ind[,2]-1)*t+tr.ind[,1])
  tr.stp <- cbind(as.matrix(tr.stp),'index'=(tr.stp[,2]-1)*t+tr.stp[,1])
  treated_<- rbind(tr.ind,tr.stp)[,c(3,4),drop=F]
  treated<- cbind(treated_,'rel_net_eff' = (treated_[,1]+e[treated_[,2]])/error.sd)
  
  errors <- c(e,e)
  names(errors)<- c(paste('iis',n_,t_,sep = '.'),paste('sis',n_,t_,sep = '.'))
  # make data pretty
  data      <- cbind(n_, t_, y, X_)
  if(nx>0){
    colnames(data)<- c("n", "t", "y", paste0('x',1:nx))
  }else{
    colnames(data)<- c("n", "t", "y")
  }
  
  
  # valuate exact treatment timing
  colnames(treated)<-c('size','index','rel_net_eff')
  
  # structure output
  sim       <- list()
  sim$data  <- data
  sim$true.b<- b_
  sim$errors<- errors
  if(const){
    sim$true.const <- b0
  }
  if(ife){
    sim$true.ife <- bife
  }
  if(tfe){
    sim$true.tfe <- btfe
  }
  sim$tr.idx<- treated
  
  return(sim)
}

