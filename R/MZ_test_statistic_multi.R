#' @title MZ_test_statistic_multi
#'
#' @description computes the MZ test statistic for multiple time series from an evaluation sample
#'
#'
#' @param y the realizations, a list of length G (the number of time series) with vectors with realizations (call their lengths T)
#' @param tau vector of quantile levels, for which forecasts are provided
#' @param H maximum forecast horizon
#' @param yhat the forecasts, a list of length G containing lists of length length(tau), which contain T x H matrices with forecasts for a specific quantile for horizons 1 to H corresponding to the realizations (e.g. row 1 contains 1- up to H-period-ahead forecasts for time 1, where the first value in the realizations vector is also for time 1)
#'
#' @return a list with two elements;
#' * `Uhat_MZ` value of the MZ test statistic
#' * `mhat` vector of the individual mhats (see paper for their definition)
#'
#'
#' @examples
#' @export
#' @importFrom quantreg rq
#'




MZ_test_statistic_multi <- function(y, tau, H, yhat){

T <- length(y[[1]])
G <- length(y)

#generate empty vector to save mhats
m=matrix(NA,2*length(tau)*H*G,1)

# calculate joint test statistic as sum of individual test statistics and determine vector of mhats

teststat <- 0

for (g in 1:G){

  indresult <- MZ_test_statistic(y=y[[g]], tau=tau, H=H, yhat=yhat[[g]])

  teststat <- teststat + indresult[[1]]

  m[((g-1)*2*H*length(tau)+1):(g*2*H*length(tau))] <- indresult[[3]]


}


list(teststat=teststat, mhat=m)
}
