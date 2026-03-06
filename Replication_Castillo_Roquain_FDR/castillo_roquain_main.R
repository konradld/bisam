# programme pour cluster EB 

procedure=function(prior,alpha,n,sn,nbsimu){
  m=n
  if (prior=="cauchy")  #changement 23102017 quasiCauchy<-Cauchy
    g=function(x){
      res=1/(2*sqrt(2*pi))
      if(x!=0) res=(1/sqrt(2*pi))*(1-exp(-x^2/2))/x^2
      return(res)
    }
  #   plot(seq(-5,5,0.01),sapply(seq(-5,5,0.01),g),type="l",ylim=c(0,0.5))
  #   lines(seq(-5,5,0.01),sapply(seq(-5,5,0.01),function(x)
  #     integrate(Vectorize(function(tau) dnorm(tau)*dcauchy(x-tau)),-Inf,Inf)$value))
  
  if (prior=="laplace") #changement 23102017 laplace<-formule
    g=function(x) (exp(1/2)/2)*(exp(-x)*pnorm(x-1)+exp(x)*(1-pnorm(x+1))) 
  #   plot(seq(-5,5,0.01),sapply(seq(-5,5,0.01),g),type="l")
  #   lines(seq(-5,5,0.01),sapply(seq(-5,5,0.01),function(x)
  #   integrate(Vectorize(function(tau) dnorm(tau)*exp(-abs(x-tau))/2),-Inf,Inf)$value),col=2,lty=2)
  #   
  Gbar=function(u) integrate(Vectorize(function(tau) g(tau)),u,40)$value
  #plot(seq(-5,5,0.01),sapply(seq(-5,5,0.01),Gbar))
  
  nbmethod=3*length(seqt)+1
  nbmoy=length(seqmoy)
  
  FDR=matrix(0,nbmethod,nbmoy)
  Pow=matrix(0,nbmethod,nbmoy)
  Risk=matrix(0,nbmethod,nbmoy)
  
  i=1
  for (moy in seqmoy){
    mu=0
    if(sn>0) mu=c(rep(0,m-sn),rep(moy,sn))
    for (simu in 1:nbsimu){
      X=mu+rnorm(m)
      
      ## method 1 : BH
      pvalues=2*(1-pnorm(abs(X)))
      pvaluesort=sort(pvalues)
      set=which(pvaluesort<=alpha*1:m/m)
      if(length(set)>0){ 
        BHth=alpha*max(set)/m
        FDR[1,i]=FDR[1,i]+(sum((mu==0)&(pvalues<=BHth))/sum(pvalues<=BHth))/nbsimu
        if(sum(mu>0)) Pow[1,i]=Pow[1,i]+(sum((mu>0)&(pvalues<=BHth))/sum(mu>0))/nbsimu
        Risk[1,i]=Risk[1,i]+(sum((mu==0)&(pvalues<=BHth)) + sum((mu>0)&(pvalues>BHth)))/nbsimu
      }
      
      ## Computing lvalues
      JS=ebayesthresh(X,verbose=TRUE,prior)
      lvalues=(1-JS$w)*dnorm(X)/((1-JS$w)*dnorm(X)+ JS$w*sapply(X,FUN=g))
      method=2
      
      for (t in seqt){
        
        ## method thresholding : lvalue <= t
        if(sum(lvalues<=t)>0){
          FDR[method,i]=FDR[method,i]+(sum((mu==0)&(lvalues<=t))/sum(lvalues<=t))/nbsimu
          if(sum(mu>0)) Pow[method,i]=Pow[method,i]+(sum((mu>0)&(lvalues<=t))/sum(mu>0))/nbsimu
        }
        Risk[method,i]=Risk[method,i]+(sum((mu==0)&(lvalues<=t)) + sum((mu>0)&(lvalues>t)))/nbsimu
        method=method+1
      } 
      
      for (t in seqt){
        
        ## SC method t
        lvaluesort=sort(lvalues)
        cummean=cumsum(lvaluesort)/1:m
        set=which(cummean<=t)
        selected=rep(FALSE,m)
        if(length(set)>0){ 
          kSC=max(set)
          selected[(order(lvalues)[1:kSC])]=TRUE
          FDR[method,i]=FDR[method,i]+(sum((mu==0)&selected)/sum(selected))/nbsimu
          if(sum(mu>0)) Pow[method,i]=Pow[method,i]+(sum((mu>0)&(selected))/sum(mu>0))/nbsimu
        }
        Risk[method,i]=Risk[method,i]+(sum((mu==0)&(selected)) + sum((mu>0)&(!selected)))/nbsimu
        method=method+1
      }
      
      qvalues=(1-JS$w)*(1-pnorm(abs(X)))/((1-JS$w)*(1-pnorm(abs(X)))+ JS$w*sapply(abs(X),FUN=Gbar))
      
      for (t in seqt){
        
        ## method thresholding : qvalue <= t
        if(sum(qvalues<=t)>0){ #corrected 15/08/2017
          FDR[method,i]=FDR[method,i]+(sum((mu==0)&(qvalues<=t))/sum(qvalues<=t))/nbsimu
          if(sum(mu>0)) Pow[method,i]=Pow[method,i]+(sum((mu>0)&(qvalues<=t))/sum(mu>0))/nbsimu
        }
        Risk[method,i]=Risk[method,i]+(sum((mu==0)&(qvalues<=t)) + sum((mu>0)&(qvalues>t)))/nbsimu
        method=method+1
      } 
      
      
    }       
    i=i+1
  }
  return(list(FDR=FDR,Pow=Pow,RiskHam=Risk)) 
}


makefilename=function(prior,alpha,n,sn,nbsimu,Risk){
  filename=paste(Risk,"EBmovemoy",sep="")
  filename=paste(filename,"_prior",prior,"_alpha",alpha,"_n",n,"_snsurn",sn/n,"_nbsimu",nbsimu,sep="")
  filename=gsub("\\.","_",filename)
  return(filename)
}

mainfunction=function(param){
  
  prior=param[1]
  alpha=as.numeric(param[2])
  n=as.numeric(param[3])
  sn=ceiling(as.numeric(param[4])*n)
  
  res=procedure(prior,alpha,n,sn,nbsimu)
  filename=makefilename(prior,alpha,n,sn,nbsimu,"")
  save(res,file=paste(filename,sep=""))
  
  print("ok")
  
}

# parametres
seqmoy=c(seq(0.01,0.901,0.1),1:10)
seqt=c(0.05,0.1,0.2,0.5,0.75)

# parametres parallelized 
seqprior=c("cauchy","laplace")
seqalpha=c(0.2)
n=10^5
seqn=c(n)
seqsnsurn=c(0.01,0.03,0.06,0.1,0.2,0.3,0.4,0.5,0.6)


nbsimu=50

#creation de la liste de para
parameter=list()
for (prior in seqprior){
  for (alpha in seqalpha){
    for (n in seqn){
      for (snsurn in seqsnsurn){
        parameter=c(parameter,list(c(prior,alpha,n,snsurn)))
      }}}}

library(EbayesThresh)
library(parallel)
library(foreach)
library(doParallel)

# run
mclapply(parameter,  mainfunction, mc.cores=18)


# Make plots

for (prior in seqprior){
  for(alpha in seqalpha){
    for (n in seqn){
      for (snsurn in seqsnsurn){
        sn=ceiling(snsurn*n)
        filename=makefilename(prior,alpha,n,sn,nbsimu,"")
        load(filename)
        
        # FDR
        resFDR=res$FDR     
        filename=makefilename(prior,alpha,n,sn,nbsimu,"FDR")
        filename=paste(filename,".pdf",sep="")
        pdf(filename)
        plot(seqmoy,resFDR[1,],lty=2,col=1,xlab="",ylab="",cex.axis=2,type="l",lwd=3,ylim=c(0,1),xlim=c(0,max(seqmoy)+2))
        for (meth in 2:(length(seqt)+1)){
          lines(seqmoy,resFDR[meth,],lty=1,col=meth,xlab="",ylab="",cex.axis=2,lwd=3)
        }
        for (meth in (length(seqt)+2):(2*length(seqt)+1)){
          lines(seqmoy,resFDR[meth,],lty=4,col=meth-length(seqt),xlab="",ylab="",cex.axis=2,lwd=3)
        }
        for (meth in (2*length(seqt)+2):(3*length(seqt)+1)){
          lines(seqmoy,resFDR[meth,],lty=3,col=meth-2*length(seqt),xlab="",ylab="",cex.axis=2,lwd=3)
        }
        
        abline(h=alpha)
        lg=c("BH",paste("lval",seqt),paste("SC",seqt),paste("qval",seqt))
        legend("topright",col=c(1:(length(seqt)+1),2:(length(seqt)+1),2:(length(seqt)+1)),
               lty=c(2,(2:(length(seqt)+1))*0+1,(2:(length(seqt)+1))*0+4,(2:(length(seqt)+1))*0+3),
               legend=lg,bg="white",cex=0.7)
        dev.off()
        
        #Pow
        resPow=res$Pow     
        filename=makefilename(prior,alpha,n,sn,nbsimu,"Pow")
        filename=paste(filename,".pdf",sep="")
        pdf(filename)
        plot(seqmoy,resPow[1,],lty=2,col=1,xlab="",ylab="",cex.axis=2,type="l",lwd=3,ylim=c(0,1),xlim=c(0,max(seqmoy)+2))
        for (meth in 2:(length(seqt)+1)){
          lines(seqmoy,resPow[meth,],lty=1,col=meth,xlab="",ylab="",cex.axis=2,lwd=3)
        }
        for (meth in (length(seqt)+2):(2*length(seqt)+1)){
          lines(seqmoy,resPow[meth,],lty=4,col=meth-length(seqt),xlab="",ylab="",cex.axis=2,lwd=3)
        }
        for (meth in (2*length(seqt)+2):(3*length(seqt)+1)){
          lines(seqmoy,resPow[meth,],lty=3,col=meth-2*length(seqt),xlab="",ylab="",cex.axis=2,lwd=3)
        }
        
        abline(h=alpha)
        lg=c("BH",paste("lval",seqt),paste("SC",seqt),paste("qval",seqt))
        legend("topright",col=c(1:(length(seqt)+1),2:(length(seqt)+1),2:(length(seqt)+1)),
               lty=c(2,(2:(length(seqt)+1))*0+1,(2:(length(seqt)+1))*0+4,(2:(length(seqt)+1))*0+3),
               legend=lg,bg="white",cex=0.7)
        dev.off()
        
        #RiskHam
        resRiskHam=res$RiskHam     
        filename=makefilename(prior,alpha,n,sn,nbsimu,"RiskHam")
        filename=paste(filename,".pdf",sep="")
        pdf(filename)
        plot(seqmoy,resRiskHam[1,],lty=2,col=1,xlab="",ylab="",cex.axis=2,type="l",lwd=3,ylim=c(0,max(resRiskHam)),xlim=c(0,max(seqmoy)+2))
        for (meth in 2:(length(seqt)+1)){
          lines(seqmoy,resRiskHam[meth,],lty=1,col=meth,xlab="",ylab="",cex.axis=2,lwd=3)
        }
        for (meth in (length(seqt)+2):(2*length(seqt)+1)){
          lines(seqmoy,resRiskHam[meth,],lty=4,col=meth-length(seqt),xlab="",ylab="",cex.axis=2,lwd=3)
        }    
        for (meth in (2*length(seqt)+2):(3*length(seqt)+1)){
          lines(seqmoy,resRiskHam[meth,],lty=3,col=meth-2*length(seqt),xlab="",ylab="",cex.axis=2,lwd=3)
        }
        
        abline(h=alpha)
        lg=c("BH",paste("lval",seqt),paste("SC",seqt),paste("qval",seqt))
        legend("topright",col=c(1:(length(seqt)+1),2:(length(seqt)+1),2:(length(seqt)+1)),
               lty=c(2,(2:(length(seqt)+1))*0+1,(2:(length(seqt)+1))*0+4,(2:(length(seqt)+1))*0+3),
               legend=lg,bg="white",cex=0.7)
        dev.off()
        
        
      }
    }
  }
}
