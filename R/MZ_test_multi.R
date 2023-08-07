#' @title MZ_test_multi
#'
#' @description executes the MZ test for multiple time series for a given evaluation sample
#'
#'
#' @param y the realizations, a list of length G (the number of time series) with vectors with realizations (call their lengths T)
#' @param tau vector of quantile levels, for which forecasts are provided
#' @param H maximum forecast horizon
#' @param yhat the forecasts, a list of length G containing lists of length length(tau), which contain T x H matrices with forecasts for a specific quantile for horizons 1 to H corresponding to the realizations (e.g. row 1 contains 1- up to H-period-ahead forecasts for time 1, where the first value in the realizations vector is also for time 1)
#' @param B number of bootstrap draws for the moving block bootstrap
#' @param l block length for the moving block bootstrap
#'
#' @return a list with two elements;
#' * `resultsMZ` results of the test, a named one matrix with one row containing five numbers:
#'     * `Stat` value of the MZ test statistic
#'     * `90%/95%/99%` 0.9/0.95/0.99 quantile of bootstrap distribution
#'     * `p-value` p-value of the test
#'
#'
#' @examples
#' @export
#' @importFrom quantreg rq




MZ_test_multi <- function(y, tau, H, yhat,B,l){

  # calculate test statistic and save mhat

  T <- length(y[[1]])
  G <- length(y)

  result <- MZ_test_statistic_multi(y=y, tau=tau, H=H, yhat=yhat)

  teststat <- result[[1]]
  m <- result[[2]]

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
    ystar=lapply(1:G,function(x){y[[x]][I]})
    yhatstar=lapply(1:G,function(y){lapply(1:length(tau),function(x){as.matrix(yhat[[y]][[x]][I,])})})


    # compute bootstrap test statistic

    bootstats[boot] <- Pstar*sum((MZ_test_statistic_multi(y=ystar, tau=tau, H=H, yhat=yhatstar)[[2]] - m)^2)

  }

  # generate results matrix/row vector

  resultsMZ=matrix(NA,1,5)
  colnames(resultsMZ)=c("Stat","90%","95%","99%","p-value")

  # save results

  resultsMZ[,1] <- teststat
  resultsMZ[,2:4] <- quantile(bootstats,probs=c(0.9,0.95,0.99))
  resultsMZ[,5] <- mean(bootstats>teststat)


  resultsMZ

}
