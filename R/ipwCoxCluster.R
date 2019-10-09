#' Inference of marginal HR in IPW Cox model based on CSV with clustering
#'
#' Inference of marginal hazard ratios (HR) in inverse probability weighted (IPW) Cox model with clustered data, under both the conventional inverse probability weights and the stabilized weights. Corrected sandwich variance (CSV) estimation method is used for variance estimation of estimated log marginal hazard ratios. 
#'@param data A clustered dataset to be analyzed in the form of R data frame.
#'@param indID A column name indicating the cluster variable. Observations from the same cluster have the same indID value. 
#'@param indA A column name indicating the treatment variable.
#'@param indX A vector of column names indicating the covariates included in the propensity score model.
#'@param indStatus A column name indicating the non-censoring status (1 if observed and 0 if censored).
#'@param indTime A column name indicating the outcome variable, i.e., min(true event time, censoring time).
#'@param ties A character string indicating the method ("efron","breslow",or "exact") to deal with tied events for point estimation; the default is "breslow". For variance estimation, Breslow method is used to deal with tied events. 
#'@param confidence A confidence level between 0 and 1; the default is 0.95 corresponding to a 95 per cent confidence interval.
#'@return A matrix of inference results from inverse probability weighted Cox model with clustered data. The first and the second rows report log marginal hazard ratio estimate and associated corrected sandwich based standard error, marginal hazard ratio estimate and associated normality-based confidence interval, under conventional inverse probability weights and stabilized weights, respectively.
#'@import survival
#'@import stats
#'
#'@examples
#'#simulate a dataset under marginal hazard ratio 1.5 with clustering
#'set.seed(100)
#'n=300
#'oneGene<-function(id){
#'  tau=0.7
#'  alpha=(1/tau-1)/2
#'  u1=runif(1)
#'  u2=runif(1)
#'  af=(1-u2)^(-1/alpha)
#'  T0v1=-log(1-u2)
#'  T0v2=alpha*log((1-af)+af*((1-u1)^(-1/(1+alpha))))
#'  T0=c(T0v1,T0v2)
#'  ZZZ=rep(rnorm(1),2)
#'  X1before= 0.5*(T0+0.2)/(T0+1)+0.3*ZZZ
#'  X1=mean(X1before)
#'  X2= 1/log(1.3*T0+3)-0.3*ZZZ
#'  X3= rbinom(2,1,0.3+0.5/(T0+1))  
#'  A=rbinom(2,1,1/(1+exp(1-2*X1-X2-X3)))
#'  Ttime <- T0*exp(-log(1.5)*A)
#'  rateC=0.5
#'  C <- rexp(2, rate=rateC)
#'  time <- pmin(Ttime, C)
#'  status <- as.numeric(Ttime <= C)
#'  cbind(id=id,time=time,status=status,A=A,X1=X1,X2=X2,X3=X3)
#'} 
#'da=as.data.frame(do.call("rbind",lapply(1:n,oneGene)))
#'head(da)
#'#inference results for marginal hazard ratio 
#'ipwCoxCluster(data=da,indID="id",indA="A",indX=c("X1","X2","X3"),indStatus="status",indTime="time")
#'
#'@export
ipwCoxCluster<-function(data,indID,indA,indX,indStatus,indTime,ties="breslow",confidence=0.95){
  dat=data
  N=nrow(dat)
  dat$id=dat[,indID]
  dat$A=dat[,indA]
  dat$time=dat[,indTime]
  dat$status=dat[,indStatus]
  
  nX=length(indX)+1
  covX0=dat[,indX]
  A=dat$A
  psmd=glm(A~.,family="binomial",data=as.data.frame(cbind(A,covX0)))
  psfit=predict(psmd, type = "response")
  
  dat$wt=dat$A/psfit+(1-dat$A)/(1-psfit)
  prevA=mean(dat$A)
  dat$swt=prevA*dat$A/psfit+(1-prevA)*(1-dat$A)/(1-psfit)
  
  fit <- coxph(Surv(time, status) ~ A, weights=dat$wt,data=dat, ties=ties)
  logHR=fit$coefficients
   
  fits <- coxph(Surv(time, status) ~ A, weights=dat$swt,data=dat, ties=ties)
  logHRs=fits$coefficients
  
  eventid=which(dat$status==1)
  covX=as.matrix(cbind(rep(1,N),covX0))
  
  dgvec=-dat$A*(1-psfit)/psfit+(1-dat$A)*psfit/(1-psfit)
  dgvecS=-prevA*dat$A*(1-psfit)/psfit+(1-prevA)*(1-dat$A)*psfit/(1-psfit)
  WdevR=t(diag(dgvec)%*%covX)
  WdevRS=t(diag(dgvecS)%*%covX)
  WfS=dat$A/psfit-(1-dat$A)/(1-psfit)
  
  A11x=rep(0,N)
  A12x=matrix(0,nX,N)
  A11xS=rep(0,N)
  A12xS=matrix(0,nX,N)
  A13xS=rep(0,N)
  
  for (x in eventid){
    idrs=which(dat$time>=dat$time[x])
    s0x=sum(dat$wt[idrs]*exp(dat$A[idrs]*logHR))
    s1x=sum(dat$wt[idrs]*exp(dat$A[idrs]*logHR)*dat$A[idrs])
    A11x[x]=dat$wt[x]*(s1x/s0x-s1x^2/(s0x^2))
    A12a=(dat$A[x]-s1x/s0x)*WdevR[,x]
    
    s0xS=sum(dat$swt[idrs]*exp(dat$A[idrs]*logHRs))
    s1xS=sum(dat$swt[idrs]*exp(dat$A[idrs]*logHRs)*dat$A[idrs])
    A11xS[x]=dat$swt[x]*(s1xS/s0xS-s1xS^2/(s0xS^2))
    A12aS=(dat$A[x]-s1xS/s0xS)*WdevRS[,x]
    A13aS=(dat$A[x]-s1xS/s0xS)*WfS[x]
    
    if (length(idrs)==1){
      A12c=(WdevR[,idrs]*exp(dat$A[idrs]*logHR))*s1x
      A12d=WdevR[,idrs]*(exp(dat$A[idrs]*logHR)*dat$A[idrs])
      
      A12cS=(WdevRS[,idrs]*exp(dat$A[idrs]*logHRs))*s1xS
      A12dS=WdevRS[,idrs]*(exp(dat$A[idrs]*logHRs)*dat$A[idrs])
      A13cS=(WfS[idrs]*exp(dat$A[idrs]*logHRs))*s1xS
      A13dS=WfS[idrs]*(exp(dat$A[idrs]*logHRs)*dat$A[idrs])
    } else{
      A12c=(WdevR[,idrs]%*%exp(dat$A[idrs]*logHR))*s1x
      A12d=WdevR[,idrs]%*%(exp(dat$A[idrs]*logHR)*dat$A[idrs])
      
      A12cS=(WdevRS[,idrs]%*%exp(dat$A[idrs]*logHRs))*s1xS
      A12dS=WdevRS[,idrs]%*%(exp(dat$A[idrs]*logHRs)*dat$A[idrs])
      A13cS=(WfS[idrs]%*%exp(dat$A[idrs]*logHRs))*s1xS
      A13dS=WfS[idrs]%*%(exp(dat$A[idrs]*logHRs)*dat$A[idrs])
    }
    
    
    A12b=dat$wt[x]*(A12d/s0x-A12c/(s0x^2))
    A12x[,x]=-(A12a-A12b)
    
    A12bS=dat$swt[x]*(A12dS/s0xS-A12cS/(s0xS^2))
    A12xS[,x]=-(A12aS-A12bS)
    A13bS=dat$swt[x]*(A13dS/s0xS-A13cS/(s0xS^2))
    A13xS[x]=-(A13aS-A13bS)
  }
  
  A11A12=c(sum(A11x),apply(A12x,1,sum))
  A11A12A13s=c(sum(A11xS),apply(A12xS,1,sum),sum(A13xS))
  
  sumsquare<-function(u){
    u%*%t(u)
  }
  
  A22mat=apply(covX,1,sumsquare)%*%(psfit*(1-psfit))
  A22=matrix(apply(A22mat,1,sum),nX,nX)
  
  AA=as.matrix(rbind(A11A12,cbind(0,A22)))
  invAA=solve(AA)
  
  AAs=as.matrix(rbind(A11A12A13s,cbind(0,A22,0),c(rep(0,nX+1),N)))
  invAAs=solve(AAs)
  
  eventPa=subset(dat,dat$status==1)
  eventimes=unique(eventPa$time)
  
  RScol5=eventimes
  RScol3=rep(0,length(eventimes))
  RScol4=RScol3
  RScol1=RScol3
  RScol2=RScol3
  
  RScol5s=eventimes
  RScol3s=rep(0,length(eventimes))
  RScol4s=RScol3s
  RScol1s=RScol3s
  RScol2s=RScol3s
  
  
  for (i in 1:length(eventimes)){
    idc=which(dat$time==eventimes[i]&dat$status==1)
    RD1=dat[idc,]
    
    RScol1[i]=sum(RD1$wt*RD1$A)
    RScol2[i]=sum(RD1$wt)
    
    RSind=which(dat$time>=eventimes[i])
    RS1=dat[RSind,]
    
    RScol3[i]=sum(RS1$wt*RS1$A)
    RScol4[i]=sum(RS1$wt*(1-RS1$A))
    
    RScol1s[i]=sum(RD1$swt*RD1$A)
    RScol2s[i]=sum(RD1$swt)
    RScol3s[i]=sum(RS1$swt*RS1$A)
    RScol4s[i]=sum(RS1$swt*(1-RS1$A))   
  }
   
  risksetFull=data.frame(sumAiWi=RScol1,sumWi=RScol2,sumW1=RScol3,sumW0=RScol4,time=RScol5)
  risksetFulls=data.frame(sumAiWi=RScol1s,sumWi=RScol2s,sumW1=RScol3s,sumW0=RScol4s,time=RScol5s)
  
  indevent=which(dat$status==1)
  gg=rep(0,nrow(dat))
  ggS=gg
  
  HRest=exp(logHR)
  HRests=exp(logHRs)
  
  for (i in indevent){
    ind1=which(risksetFull$time==dat$time[i])
    vv=(risksetFull$sumW1*HRest)/(risksetFull$sumW1*HRest+risksetFull$sumW0)
    gg[i]=dat$wt[i]*(dat$A[i]- vv[ind1])
    
    ind1=which(risksetFulls$time==dat$time[i])
    vv=(risksetFulls$sumW1*HRests)/(risksetFulls$sumW1*HRests+risksetFulls$sumW0)
    ggS[i]=dat$swt[i]*(dat$A[i]- vv[ind1])
  }
  
  
  rs11=rep(0,nrow(dat))
  rs12=rs11
  
  rs11s=rep(0,nrow(dat))
  rs12s=rs11s
  
  for (i in 1:nrow(dat)){
    ind2=which(risksetFull$time<=dat$time[i])
    rsred=risksetFull[ind2,]
    rs11[i]=dat$wt[i]*dat$A[i]*exp(log(HRest)*dat$A[i])*sum(rsred$sumWi/(rsred$sumW1*HRest+rsred$sumW0))
    rs12[i]=dat$wt[i]*exp(log(HRest)*dat$A[i])*sum(rsred$sumWi*(rsred$sumW1*HRest)/((rsred$sumW1*HRest+rsred$sumW0)^2))
    
    ind2=which(risksetFulls$time<=dat$time[i])
    rsred=risksetFulls[ind2,]
    rs11s[i]=dat$swt[i]*dat$A[i]*exp(log(HRests)*dat$A[i])*sum(rsred$sumWi/(rsred$sumW1*HRests+rsred$sumW0))
    rs12s[i]=dat$swt[i]*exp(log(HRests)*dat$A[i])*sum(rsred$sumWi*(rsred$sumW1*HRests)/((rsred$sumW1*HRests+rsred$sumW0)^2))
  }
  
  eta=gg-rs11+rs12
  etaS=ggS-rs11s+rs12s
  
  covX=as.matrix(cbind(rep(1,N),covX0))
  
  matpi=diag(dat$A-psfit)%*%covX
   
  sumsquare<-function(u){
    u%*%t(u)
  }
  
  fornow=as.data.frame(cbind(id=dat$id,eta=eta,etaS=etaS,matpi,sigm=dat$A-prevA))
  
  aggreOut=aggregate(.~id, fornow, sum)
  naggr=ncol(aggreOut)
  
  bbmat=aggreOut[,-c(1,3,naggr)]
  bbmatS=aggreOut[,-c(1,2)]
  
  oot=apply(bbmat,1,sumsquare)
  BB=matrix(apply(oot,1,sum),nX+1,nX+1)
  
  propVar=invAA%*%BB%*%t(invAA)
  proposeStdErr=(diag(propVar)^0.5)[1]
  
  ootS=apply(bbmatS,1,sumsquare)
  BBs=matrix(apply(ootS,1,sum),nX+2,nX+2)
  
  propVarS=invAAs%*%BBs%*%t(invAAs)
  proposeStdErrS=(diag(propVarS)^0.5)[1]
  
  lowProp=logHR-qnorm(1-(1-confidence)/2)*proposeStdErr
  upProp=logHR+qnorm(1-(1-confidence)/2)*proposeStdErr
  
  lowPropS=logHRs-qnorm(1-(1-confidence)/2)*proposeStdErrS
  upPropS=logHRs+qnorm(1-(1-confidence)/2)*proposeStdErrS
  
  mtd=c("conventional weights","stabilized weights")
  est=c(logHR,logHRs)
  hrest=exp(est)
  se=c(proposeStdErr,proposeStdErrS)
  low=exp(c(lowProp,lowPropS))
  up=exp(c(upProp,upPropS))
  output=cbind(est,se,hrest,low,up)
  colnames(output)=c("log HR Estimate","Standard Error","HR Estimate",
    paste("HR ", confidence*100,"% CI", "-low", sep =""),
    paste("HR ", confidence*100,"% CI", "-up", sep =""))
  rownames(output)=mtd
  output
}

 
