# osdepth
R functions to compute the origin-scale depth and plot the OSD-chart

```{r echo=T,eval=T,fig.width=7,fig.height=6,message=FALSE,warning=FALSE,results="hide"}
source("https://raw.githubusercontent.com/icascos/osdepth/master/OSD_chart.R")
require(qcc)
military=read.csv("https://raw.githubusercontent.com/icascos/osdepth/master/military.csv",
                  sep=",",header=TRUE,colClasses=c(trial="logical"))
m.trial<- subset(military,trial==TRUE)
m.new <- subset(military,trial==FALSE)
samples.dist.trial=qcc.groups(m.trial$dist,m.trial$sample)
samples.dist.new=qcc.groups(m.new$dist,m.new$sample)
par(mfrow=c(2,2))
osd.chart(samples.dist.trial,newdata=samples.dist.new)
osd.chart(samples.dist.trial,type="param",newdata=samples.dist.new)
osd.chart(samples.dist.trial,type="param",newdata=samples.dist.new,
          theta=162,lambda=836)
set.seed(1)
osd.chart(samples.dist.trial,type="boot",newdata=samples.dist.new,
          alpha=0.002,B=10000)
```
