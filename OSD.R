# univariate zonoid depth of s w.r.t. sample x

zdepth <- function(x,s,prob=rep(1/length(s),length(s)),ord=FALSE) {
  prob=prob/sum(prob)
  if(ord==FALSE) {or<-order(s);prob<-prob[or];s<-s[or]}
  n <- length(s)
  if((x<s[1])|(x>s[n])) {return(0)}
  me <- sum(s*prob)
  if(x==sum(s*prob)) {return(1)}
  if(x>me) {x <- -x ; s <- -s[n:1] ; prob <- prob[n:1]}
  i <- 1
  ac <- s[i]*prob[i]
  ap <- prob[i]
  while(x*ap>=ac) {i <- i+1 ; ap <- ap+prob[i] ; ac <- ac+s[i]*prob[i]}
  return((ap*s[i]-ac)/(s[i]-x))
}


# cdf of shifted exponential with parameters theta and lambda evaluated on x

psexp=function(x,theta=0,lambda=1){if(x<theta){return(0)};return(1-exp(-(x-theta)/lambda))}


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

depthsexp=function(t1,t2,theta=0,lambda=1){
  if(t1<theta){return(0)}
  return((1-psexp(t1,theta=theta,lambda=lambda))*depthexp(t2,lambda=lambda))
}

# origin-scale depth of of samples in rows of matrix samples w.r.t. shifted exponential with parameters theta and lambda

depthsexp.sample=function(samples,theta=0,lambda=1){
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
