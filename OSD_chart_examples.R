source("https://raw.githubusercontent.com/icascos/osdepth/master/OSD_chart.R")

# data for examples. Source: Kao (2010) Normalization of the origin-shifted exponential distribution for control chart  construction, Journal of Applied Statistics 37, 1067-1087.

HPM.data <- c(51.7, 54.3, 62.4, 52.7, 53.7, 58.9, 54.7, 60.4, 52.1, 52.1,
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

samples.current <- matrix(HPM.data, ncol = 8, byrow = TRUE)

#require(qcc)
#HPM <- read.csv("https://raw.githubusercontent.com/icascos/osdepth/master/HPM.csv", sep=",", header = TRUE, colClasses = c(trial="logical"))
#samples.current <- qcc.groups(HPM$current[trial], HPM$sample[trial])


theta.h <- min(samples.current)
lambda.h <- mean(samples.current)-theta.h
set.seed(12)
samples.current.new <- matrix(rsexp(n = 80, theta = theta.h+2, lambda = lambda.h+5), ncol = 8)

# Some examples

par(mfrow = c(2,2))

osd.chart(samples.current)
osd.chart(samples.current, type = "param", newdata = samples.current.new)
set.seed(1)
osd.chart(samples.current, type = "boot", alpha = 0.02, newdata = samples.current.new)
set.seed(1)
osd.chart(samples.current, type = "EWMA", alpha = 0.002, newdata = samples.current.new)



require(qcc)

military <- read.csv("https://raw.githubusercontent.com/icascos/osdepth/master/military.csv", sep=",", header=TRUE, colClasses=c(trial="logical"))
m.trial <- subset(military, trial == TRUE)
m.new <- subset(military, trial == FALSE)
samples.dist.trial <- qcc.groups(m.trial$dist, m.trial$sample)
samples.dist.new <- qcc.groups(m.new$dist, m.new$sample)


par(mfrow = c(2,2))
osd.chart(samples.dist.trial, newdata = samples.dist.new)
osd.chart(samples.dist.trial, type = "param", newdata = samples.dist.new, theta = 162, lambda = 836)
set.seed(1)
osd.chart(samples.dist.trial, type = "boot", alpha = 0.002, newdata = samples.dist.new, B = 10000)
set.seed(1)
osd.chart(samples.dist.trial, type = "EWMA", alpha = 0.002, newdata = samples.dist.new, theta = 162, lambda = 836)
