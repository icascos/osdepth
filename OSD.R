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

# origin-scale depth of of samples in rows of matrix samples w.r.t. shifted exponential with parameters theta and lambda

depthsexp.sample=function(samples,theta=1,lambda=1){
  depthsexp.sample=vector(length=nrow(samples))
  for(i in 1:nrow(samples)){
    depthsexp.sample[i]=depthsexp(t1=min(samples[i,]),t2=mean(samples[i,])-min(samples[i,]),theta,lambda)
  }
  return(depthsexp.sample)
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