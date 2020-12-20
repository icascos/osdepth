# OSD-chart depicts an OSD-chart from data matrix samples. Each row corresponds to a sample and 
# the number of columns is the number of observations in a sample
# By default the semiparametric procedure with False Alarm Probability (alpha) equal to 0.0027 is used.
# Alternatively the semiparametric procedure with alpha 0.05 can be selected or the control limit
# is approximated by resampling, and it possible to fix any alpha and number of resamples.
# The OSD EWMA chart is plotted by selecting type EWMA. The default lambda is 0.25
# The OSD EWMA chart plotted is the parametric if theta and lambda are given, otherwise it is the nonparametric.

source("https://raw.githubusercontent.com/icascos/osdepth/master/OSD.R")

require(spc)
osd.chart=function(trials,type="sexp.002",alpha=0.05,B=min(10000,max(1000,20/alpha)),Bew=10000,delta=0.25,newdata=NULL,theta=NULL,lambda=NULL){
  if(alpha<=0|alpha>=1) alpha=0.05
  n=ncol(trials)
  Me=c(0.29675,0.48761,0.58949,0.65388,0.69842,0.73132,0.7568,0.77704,
       0.79377,0.80784,0.81989,0.82994,0.83889,0.84656,0.85375,0.8601,
       0.86584,0.87084,0.87559)[n-1]
  CL=c(0.00099,0.02437,0.07245,0.12436,0.1714,0.21218,0.24739,0.28018,
       0.30647,0.32926,0.34974,0.3695,0.38624,0.40079,0.41618,0.42871,
       0.44117,0.45323,0.46522)[n-1]
  alfa=0.002
  if(type=="sexp.0027") {
  CL=c(0.00137,0.02858,0.07989,0.13491,0.18306,0.22509,0.26083,0.29393,
    0.32069,0.34402,0.36551,0.38391,0.40126,0.41588,0.43054,0.44348,
    0.45611,0.46776,0.47965)[n-1]
    alfa=0.0027
  }
  if(type=="sexp.05") {
  CL=c(0.02541,0.12754,0.22602,0.30371,0.36534,0.41416,0.45329,0.48647,
       0.51411,0.53778,0.55908,0.57694,0.59267,0.60697,0.61969,0.63118,
       0.64155,0.6519,0.66077)[n-1]
       alfa=0.05
       }
  else if(type=="boot"){
    simulations=vector(length=B)
    simulations=osdepth(samples=matrix(sample(as.vector(trials),B*n,replace=TRUE),ncol=n),x=as.vector(trials))
    Me=quantile(simulations,.5)
    CL=quantile(simulations,alpha)
    alfa=alpha
  }
  if(type=="param"){
      theta.h=theta
      if(is.null(theta.h)) theta.h=min(trials)
      lambda.h=lambda
      if(is.null(lambda.h)) lambda.h=mean(trials)-theta.h
      depth=depthsexp.sample(samples=rbind(trials,newdata),theta=theta.h,lambda=lambda.h)
      alfa=0.002
  }
  else if(type=="EWMA"){
    if(delta<=0|delta>=1) delta=0.25
    
    if(is.null(theta)|is.null(lambda)) {
      sdepth=osdepth(samples=matrix(sample(as.vector(trials),Bew*n,replace=TRUE),ncol=n),x=as.vector(trials))
      depth2=osdepth(samples=rbind(trials,newdata),x=as.vector(trials))
      }
    else {
      sdepth=depthsexp.sample(matrix(rexp(Bew*n),ncol=n))
      depth2=depthsexp.sample(samples=rbind(trials,newdata),theta=theta,lambda=lambda)
    } 
    sdepth=c(sdepth,0)
    
    rdata=vector(length=length(depth2))
    for(i in 1:length(rdata)) {rdata[i]=mean(depth2[i]>sdepth)}
    
    s1=-delta*log(rdata)[1]+(1-delta)
    sk=s1
    depth=vapply(-log(rdata)[-1],function(x) sk <<- delta*x+(1-delta)*sk,0)
    depth=c(s1,depth)
    
    CL=sewma.crit(l=delta,L0=1/alpha,df=2)[2]
    alfa=alpha
  }
  else depth=osdepth(samples=rbind(trials,newdata),x=as.vector(trials))
  plot(depth,ylim=c(min(depth,CL),max(depth,CL)),type="l",xlab="sample",ylab="") 
  if(type=="boot") title("OSD chart nonparametric",sub=paste("alpha=",alfa," CL=",round(CL*100000)/100000),ylab="OSD")
  else if (type=="EWMA")  title(paste("OSD EWMA chart, delta=",delta),sub=paste("alpha=",alfa," CL=",round(CL*100000)/100000),ylab="EWMA")
  else if (type=="param")  title("OSD chart parametric",sub=paste("alpha=",alfa," CL=",round(CL*100000)/100000),ylab="OSD")
  else title("OSD chart semiparametric",sub=paste("alpha=",alfa," CL=",round(CL*100000)/100000),ylab="OSD")
  if(type!="EWMA") {
    for(i in 1:length(depth)) {
      if(depth[i]>=CL) {points(i,depth[i])}
      else {points(i,depth[i],pch=8)}
    }
    abline(h=Me,lty=2)
  }
  else {
    for(i in 1:length(depth)) {
      if(depth[i]<=CL) {points(i,depth[i])}
      else {points(i,depth[i],pch=8)}
    }
  }
  abline(h=CL)
  if(!is.null(newdata)) abline(v=nrow(trials)+.5,lty=3)
  text(length(depth),CL+.03,c("CL"))
}
