### Chicago Air Pollution Data ###
# For paper: Identifiability Constraints in Generalized Additive Models
# Section 5
# Alex Stringer
# 2022/05

## Global parameters ----

## Setup ----
# Packages
library(mgcv)
library(gamair)
library(tidyverse)
library(parallel)
options(mc.cores = parallel::detectCores())
data(chicago) # In package: gamair

# Paths
resultspath <- tempdir()

## Preprocess data ----
# See Wood (2017) Generalized Additive Models (the book), page 247
# Code taken from there.
lag.sum <- function(a,l0,l1) {
  n<-length(a)
  b<-rep(0,n-l1)
  for (i in 0:(l1-l0)) b <- b + a[(i+1):(n-l1+i)]
  b
}
chicago2 <- chicago[4:5114, ]
chicago2$o3 <- lag.sum(chicago$o3median,0,3)
chicago2$tmp <- lag.sum(chicago$tmpd,0,3)
chicago2$pm10 <- lag.sum(log(chicago$pm10median+40),0,3)
chicago2$so2 <- lag.sum(log(chicago$so2median+10),0,3)

## Fit models ----
# See Wood (2017) Generalized Additive Models (the book), page 248
# Sum to zero
mod_stz <- gam(
  death ~ s(time,bs="cr",k=200) + 
          te(o3,tmp,k=10) +
          # s(o3,k=8) + s(tmp,k=8) +
          s(pm10,bs="cr",k=8),
  family=poisson,
  data = chicago2,
  method = "REML"
)
# Point in middle
mod_pc <- gam(
  death ~ s(time,bs="cr",k=200,pc=mean(chicago2$time,na.rm = TRUE)) + 
    te(o3,tmp,k=10) +
    # s(o3,k=8) + s(tmp,k=8) +
    s(pm10,bs="cr",k=8,pc=mean(chicago2$pm10,na.rm = TRUE)),
  family=poisson,
  data = chicago2,
  method = "REML"
)

## Summarize models ----
# Get the termwise fits and their standard errors
numpred <- 1e03
preddat <- data.frame(
  o3 = 0,
  tmp = 0,
  pm10 = seq(min(chicago2$pm10,na.rm=TRUE),max(chicago2$pm10,na.rm=TRUE),length.out=numpred),
  time = seq(min(chicago2$time),max(chicago2$time),length.out=numpred)
)
pred_time_stz_dat <- data.frame(o3=0,tmp=0,pm10=0,time=mean(chicago2$time,na.rm = TRUE))
pred_pm10_stz_dat <- data.frame(o3=0,tmp=0,pm10=mean(chicago2$pm10,na.rm = TRUE),time=0)
pred_center_time <- predict(mod_stz,newdata = pred_time_stz_dat,type='terms')[1]
pred_center_pm10 <- predict(mod_stz,newdata = pred_pm10_stz_dat,type='terms')[3]



pred_stz <- predict(mod_stz,newdata = preddat,se.fit = TRUE,unconditional = FALSE,type = 'terms')
fit_stz_time <- pred_stz$fit[, 's(time)']
se_stz_time <- pred_stz$se.fit[, 's(time)']
fit_stz_pm10 <- pred_stz$fit[, 's(pm10)']
se_stz_pm10 <- pred_stz$se.fit[, 's(pm10)']

pred_pc <- predict(mod_pc,newdata = preddat,se.fit = TRUE,unconditional = FALSE,type = 'terms')
fit_pc_time <- pred_pc$fit[, 's(time)']
se_pc_time <- pred_pc$se.fit[, 's(time)']
fit_pc_pm10 <- pred_pc$fit[, 's(pm10)']
se_pc_pm10 <- pred_pc$se.fit[, 's(pm10)']

# Time

textsize <- 1.5

pdf(file = file.path(resultspath,'chicago-time-plot.pdf'),width=14,height=7)
plot(preddat$time,fit_stz_time-pred_center_time,type='l',xlab="Time",ylab="f(Time)",cex.axis=textsize,cex.lab=textsize,ylim = c(-.25,.07))
lines(preddat$time,fit_stz_time-pred_center_time - 2*se_stz_time,lty='dashed')
lines(preddat$time,fit_stz_time-pred_center_time + 2*se_stz_time,lty='dashed')

lines(preddat$time,fit_pc_time)
lines(preddat$time,fit_pc_time - 2*se_pc_time,lty='dotted')
lines(preddat$time,fit_pc_time + 2*se_pc_time,lty='dotted')
dev.off()

pdf(file = file.path(resultspath,'chicago-pm10-plot.pdf'),width=14,height=7)
plot(preddat$pm10,fit_stz_pm10-pred_center_pm10,type='l',xlab="PM10",ylab="f(PM10)",cex.axis=textsize,cex.lab=textsize)
lines(preddat$pm10,fit_stz_pm10-pred_center_pm10 - 2*se_stz_pm10,lty='dashed')
lines(preddat$pm10,fit_stz_pm10-pred_center_pm10 + 2*se_stz_pm10,lty='dashed')

lines(preddat$pm10,fit_pc_pm10)
lines(preddat$pm10,fit_pc_pm10 - 2*se_pc_pm10,lty='dotted')
lines(preddat$pm10,fit_pc_pm10 + 2*se_pc_pm10,lty='dotted')
dev.off()

# Average difference in standard errors
se_diff_time <- mean(se_pc_time - se_stz_time)
se_diff_pm10 <- mean(se_pc_pm10 - se_stz_pm10)

# Standard deviations of the log mean
meanpred_stz <- predict(mod_stz,se.fit = TRUE,unconditional = FALSE)
meanpred_pc <- predict(mod_pc,se.fit = TRUE,unconditional = FALSE)
mean(abs(meanpred_stz$fit - meanpred_pc$fit))
mean(meanpred_pc$se.fit - meanpred_stz$se.fit)
# Identical





