# univariate zonoid depth of s w.r.t. sample x

zdepth<-function(s,x,ord=F){
  
  if(ord==F){x=sort(x)}
  n=length(x)
  
  i=1
  if((s<x[1])|(s>x[n])){zdepth=0;return(zdepth)}
  if(s==x[1]){while((x[i+1]==s)&(i<n)){i=i+1};zdepth=i/n;return(zdepth)}
  if(s==x[n]){while((x[n-i]==s)&(i<n)){i=i+1};zdepth=i/n;return(zdepth)}
  
  m=mean(x)
  
  if(s==m){zdepth=1;return(zdepth)}
  if(n==1){zdepth=0;return(zdepth)}
  if(s<m){
    suma=x[1]
    while(((suma+x[i+1])/(i+1)<s)&(i<n)){i=i+1;suma=suma+x[i]}
    zdepth=i*(suma/i-x[i+1])/(n*(s-x[i+1]))
    return(zdepth)
  }
  
  if(s>m){
    suma=x[n]
    while(((suma+x[n-i])/(i+1)>s)&(i<(n-1))){i=i+1;suma=suma+x[n-i+1]}
    zdepth=i*(suma/i-x[n-i])/(n*(s-x[n-i]))
    return(zdepth)
  }
  
}


# cdf of shifted exponential with parameters theta and lambda evaluated on x

psexp=function(x,theta=1,lambda=1){if(x<theta){return(0)};return(1-exp(-(x-theta)/lambda))}


# zonoid depth of t w.r.t. exponential with parameter lambda

depthexp=function(t,lambda=1){
  t=t/lambda
  if(t>=1){return(exp(1-t))}
  if(t<1){
    f<-function(d,x){return(((1-d)*log(1-d)/d-x+1)^2)}
    return(optimize(f,c(0,1),tol=0.000001,x=t)$minimum)
  }
}

# origin-scale depth of (t1,t2) w.r.t. shifted exponential with parameters theta and lambda

depthsexp=function(t1,t2,theta=1,lambda=1){
  if(t1<theta){return(0)}
  return((1-psexp(t1,theta=theta,lambda=lambda))*depthexp(t2,lambda=lambda))
}

# origin-scale depth of samples in rows of matrix samples w.r.t. sample x

osdepth=function(samples,x,ord=F){
  if(ord==F){x=sort(x)}
  n=length(x)
  
  osdepth=vector(length=nrow(samples))
  
  for(i in 1:nrow(samples)){
    l=sum(x<min(samples[i,]))
    if(l==n){osdepth[i]=0}
    else {osdepth[i]=(n-l)*zdepth(mean(samples[i,]),x[(l+1):n],ord=T)/n}
  }
  return(osdepth)
}


# n random obsevations of a shifted exponential distribution with parameters theta and lambda

rsexp<-function(n,theta,lambda){return(rexp(n)*lambda+theta)}



# OSD-chart depicts a OSD-chart from data matrix samples. Each row corresponds to a sample and 
# the number of columns is the number of observations in a sample
# By default the semiparametric procedure with False Alarm Probability (alpha) equal to 0.0027 is used.
# Alternatively the semiparametric procedure with alpha 0.05 can be selected or the control limit
# is approximated by resampling, and it possible to fix any alpha and number of resamples.


osd.chart=function(samples,type="sexp.0027",alpha=0.05,B=1000){
  depth=osdepth(samples=samples,x=as.vector(samples))
  n=ncol(samples)
  Me=c(0.29675,0.48761,0.58949,0.65388,0.69842,0.73132,0.7568,0.77704,
       0.79377,0.80784,0.81989,0.82994,0.83889,0.84656,0.85375,0.8601,
       0.86584,0.87084,0.87559)[n-1]
  CL=c(0.00137,0.02858,0.07989,0.13491,0.18306,0.22509,0.26083,0.29393,
    0.32069,0.34402,0.36551,0.38391,0.40126,0.41588,0.43054,0.44348,
    0.45611,0.46776,0.47965)[n-1]
  alfa=0.0027
  if(type=="sexp.05") {
  CL=c(0.02541,0.12754,0.22602,0.30371,0.36534,0.41416,0.45329,0.48647,
       0.51411,0.53778,0.55908,0.57694,0.59267,0.60697,0.61969,0.63118,
       0.64155,0.6519,0.66077)[n-1]
       alfa=0.05
       }
  else if(type=="boot"){
    simulations=vector(length=B)
    simulations=osdepth(samples=matrix(sample(as.vector(samples),B*n,replace=TRUE),ncol=n),x=as.vector(samples))
    Me=quantile(simulations,.5)
    CL=quantile(simulations,alpha)
    alfa=alpha
  }
  plot(depth,ylim=c(min(depth,CL),1),type="l",ylab="OSD",xlab="sample",main=paste("OSD chart, Phase I, alpha=",alfa))
  for(i in 1:length(depth)) {
    if(depth[i]>=CL) {points(i,depth[i])}
    else {points(i,depth[i],pch=8)}
  }
  abline(h=CL)
  abline(h=Me,lty=2)
  text(length(depth),CL+.03,c("CL"))
}

# data for examples. Source: Kao (2010) Normalization of the origin-shifted exponential distribution for control chart construction, Journal of Applied Statistics 37, 1067-1087.

sexp.data=c(51.7, 54.3, 62.4, 52.7, 53.7, 58.9, 54.7, 60.4, 52.1, 52.1,
            54.9, 60.3, 60.9, 52.5, 58.3, 56.6, 56.8, 62.8, 52.5, 55.5,
            55.5, 57.2, 62.5, 52.1, 53.3, 52.2, 52.9, 57.4, 56.2, 52.2,
            57.3, 55.3, 65.9, 56.0, 55.9, 53.2, 55.4, 58.7, 52.5, 54.8,
            52.1, 55.5, 58.8, 52.4, 52.4, 52.2, 74.8, 59.5, 52.7, 52.7,
            51.8, 54.2, 59.5, 51.7, 52.9, 51.0, 53.5, 59.9, 51.9, 51.2,
            52.3, 55.5, 58.0, 79.2, 51.2, 50.4, 50.5, 61.9, 85.4, 50.8,
            53.8, 57.1, 60.9, 56.6, 52.4, 78.6, 58.9, 61.6, 51.0, 51.0,
            52.4, 54.8, 54.0, 51.0, 50.3, 50.7, 51.2, 59.5, 50.2, 51.9,
            52.5, 52.7, 58.6, 50.3, 82.1, 54.7, 51.5, 57.1, 51.8, 54.2,
            52.1, 88.0, 58.5, 50.3, 52.2, 52.5, 51.3, 57.6, 50.8, 51.5,
            50.7, 52.3, 64.5, 51.3, 52.8, 56.8, 54.5, 108, 58.4, 60.0)

# data with abnormal sample

sexp.data.out=c(51.7, 54.3, 62.4, 52.7, 53.7, 58.9, 54.7, 60.4, 52.1, 52.1,
                54.9, 60.3, 60.9, 52.5, 58.3, 56.6, 56.8, 62.8, 52.5, 55.5,
                55.5, 57.2, 62.5, 52.1, 53.3, 52.2, 52.9, 57.4, 56.2, 52.2,
                57.3, 55.3, 65.9, 56.0, 55.9, 53.2, 55.4, 58.7, 52.5, 54.8,
                108,  107,  107,  105,  104,  103,  102,  101,
                52.1, 55.5, 58.8, 52.4, 52.4, 52.2, 74.8, 59.5, 52.7, 52.7,
                51.8, 54.2, 59.5, 51.7, 52.9, 51.0, 53.5, 59.9, 51.9, 51.2,
                52.3, 55.5, 58.0, 79.2, 51.2, 50.4, 50.5, 61.9, 85.4, 50.8,
                53.8, 57.1, 60.9, 56.6, 52.4, 78.6, 58.9, 61.6, 51.0, 51.0,
                52.4, 54.8, 54.0, 51.0, 50.3, 50.7, 51.2, 59.5, 50.2, 51.9,
                52.5, 52.7, 58.6, 50.3, 82.1, 54.7, 51.5, 57.1, 51.8, 54.2,
                52.1, 88.0, 58.5, 50.3, 52.2, 52.5, 51.3, 57.6, 50.8, 51.5,
                50.7, 52.3, 64.5, 51.3, 52.8, 56.8, 54.5, 108, 58.4, 60.0)



# formating the data

samples=matrix(sexp.data,ncol=8,byrow=TRUE)
samples.out=matrix(sexp.data.out,ncol=8,byrow=TRUE)
set.seed(1)
samples.sim=matrix(rsexp(n=160,theta=10,lambda=5),ncol=4,byrow=TRUE)

# Some examples

par(mfrow=c(2,2))

osd.chart(samples)
osd.chart(samples,type="sexp.05")
osd.chart(samples.out)
set.seed(12)
osd.chart(samples.sim,type="boot",alpha=0.04)