#' @title MZ_test
#'
#' @description executes the MZ test for a given evaluation sample
#'
#'
#' @param y the realizations, a vector
#' @param tau vector of quantile levels, for which forecasts are provided
#' @param H maximum forecast horizon
#' @param yhat the forecasts, a list of length length(tau), which contains length(y) x H matrices with forecasts for a specific quantile for horizons 1 to H corresponding to the realizations (e.g. row 1 contains 1- up to H-period-ahead forecasts for time 1, where the first value in the realizations vector is also for time 1)
#' @param B number of bootstrap draws for the moving block bootstrap
#' @param l block length for the moving block bootstrap
#'
#' @return a list with two elements;
#' * `resultsMZ` main results of the test, a named one matrix with one row containing five numbers:
#'     * `Stat` value of the MZ test statistic
#'     * `90%/95%/99%` 0.9/0.95/0.99 quantile of bootstrap distribution
#'     * `p-value` p-value of the test
#' * `tables` additional information for each horizon-quantile combination, a list with three data frames:
#'    * `cont` individual contributions to test statistic
#'    * `alphahat` MZ regression intercepts
#'    * `betahat` MZ regression slopes
#'
#'
#' @examples
#' @export
#' @importFrom quantreg rq




MZ_test <- function(y, tau, H, yhat,B,l){

  # calculate test statistic and save tables and mhat

  T <- length(y)

  result <- MZ_test_statistic(y=y, tau=tau, H=H, yhat=yhat)

  teststat <- result[[1]]
  tables <- result[[2]]
  m <- result[[3]]

  # obtain bootstrap distribution of test statistic

  # generate empty vector to save bootstrap test statistics
  bootstats=matrix(NA,B,1)

  # set bootstrap parameters
  bl=floor(T/l)
  Pstar=l*bl

  for (boot in 1:B){

    # draw the start indices for the moving block bootstrap
    I=as.matrix(rep(sample((T-l+1),bl,replace=TRUE),each=l)+rep(seq(0,(l-1),1),bl))

    # draw realizations and forecasts
    ystar=y[I]
    yhatstar=lapply(1:length(tau),function(x){as.matrix(yhat[[x]][I,])})

    # compute bootstrap test statistic

    bootstats[boot] <- Pstar*sum((MZ_test_statistic(y=ystar, tau=tau, H=H, yhat=yhatstar)[[3]] - m)^2)

  }

  # generate results matrix/row vector

  resultsMZ=matrix(NA,1,5)
  colnames(resultsMZ)=c("Stat","90%","95%","99%","p-value")

  # save results

  resultsMZ[,1] <- teststat
  resultsMZ[,2:4] <- quantile(bootstats,probs=c(0.9,0.95,0.99))
  resultsMZ[,5] <- mean(bootstats>teststat)


  list(resultsMZ = resultsMZ, tables=tables)

}
