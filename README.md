# R package quantoptimR

## Overview

This package accompanies the paper [_Testing Quantile Forecast Optimality_](https://arxiv.org/abs/2302.02747) by Jack Fosten, Daniel Gutknecht and Marc-Oliver Pohle.

The package provides functions that execute tests for optimality of quantile forecasts over possibly multiple forecast horizons
and multiple quantiles based on Mincer-Zarnowitz quantile regressions. The tests lead to an overall decision about the optimality/calibration of a forecasting approach (over all horizons and quantiles), but also give valuable feedback about the horizons and quantiles, where miscalibration is strongest. The basic test checks the null hypothesis of autocalibration. There are two extensions: one test for stronger forms of optimality and one for multivariate forecasts.

There are three main functions:

* `MZ_test` executes the basic Mincer-Zarnowitz test for autocalibration over multiple horizons and quantiles when provided with forecasts and corresponding observations.
* `MZ_test_multi` is an extension of the basic test for multiple variables.
* `MZ_test_augmented` is an extension of the basic test, where additional regressors can be added to the Mincer-Zarnowitz regressions to test a corresponding stronger form of optimality.

## Installation

To install the package from GitHub and load it, run the following `R` code:

```
install.packages("devtools")
library(devtools)
install_github("MarcPohle/quantoptimR")
library(quantoptimR)

```

## Example

This is a simple example to illustrate the execution of the basic MZ test for 
autocalibration. The setting is the same as for the first simulations in the
paper. We generate observations from an AR(1) model with parameter b and 
generate the optimal quantile forecasts for an AR(1) model with parameter
btilde. If b equals btilde, we are under the null. If b does not equal btilde,
we are under the alternative.


```
T <- 100 # length of evaluation sample
H <- 3   # maximum forecast horizon
tau <- c(0.25,0.75) # quantile levels of interest

b <- 0.6 # AR(1) parameter
btilde <- 0.8 # AR(1) parameter used by the forecaster


# simulate series of size T+h from AR(1) process with coefficient 0.6 (no intercept)

eps <- rnorm(T+H)
y <- rep(NA,T+H)
y[1] <- eps[1] # for simplicity, let the first value of y be the first value of epsilon here
for (t in 2:(T+H)){
  y[t] <- b*y[t-1]+eps[t]
}


# forecasts

# generate a list of length(tau) of T x H matrices to save the forecasts
yhat <- lapply(1:length(tau),function(x){matrix(NA,T,H)})

# generate the optimal quantile forecasts for the parameter btilde
for (taucount in 1:length(tau)){
  for (h in 1:H){
    for (t in 1:T){
      yhat[[taucount]][t,h] <- (btilde^h)*y[t+H-h] + sqrt(1-btilde^(2*h))*qnorm(tau[taucount])
    }
  }
}


# observations

y <- y[(H+1):(T+H)]


# apply the basic MZ test

B=1000 # number of bootstrap draws
l=4     # block length for the bootstrap

MZ_test(y,tau,H,yhat,B,l)


```

Test results:

```

# $resultsMZ
#          Stat      90%     95%      99% p-value
# [1,] 433.7722 201.8357 351.499 500.5457    0.02
# 
# $tables
# $tables$cont
#      tau=0.25 tau=0.75
# h=1  39.79224 24.64671
# h=2 102.03140 61.19475
# h=3 124.63756 81.46951
# 
# $tables$alphahat
#       tau=0.25  tau=0.75
# h=1 -0.2863600 0.3717130
# h=2 -0.4013776 0.5329346
# h=3 -0.4614399 0.5750336
# 
# $tables$betahat
#        tau=0.25  tau=0.75
# h=1  0.43793208 0.6709155
# h=2  0.07306417 0.4273498
# h=3 -0.01658682 0.3042763

```
