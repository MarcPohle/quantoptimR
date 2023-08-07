#' @title MZ_test_augmented
#'
#' @description executes the augmented MZ test for a given evaluation sample and given additional regressors
#'
#'
#' @param y the realizations, a vector
#' @param tau vector of quantile levels, for which forecasts are provided
#' @param H maximum forecast horizon
#' @param yhat the forecasts, a list of length length(tau), which contains length(y) x H matrices with forecasts for a specific quantile for horizons 1 to H corresponding to the realizations (e.g. row 1 contains 1- up to H-period-ahead forecasts for time 1, where the first value in the realizations vector is also for time 1)
#' @param B number of bootstrap draws for the moving block bootstrap
#' @param l block length for the moving block bootstrap
#' @param adreg additional regressors,  a (length(y) + H - 1) X A matrix with additional regressors to include in the MZ regression, A is the number of additional regressors; the columns contain the time series of additional regressors (known at the time of forecasting, i.e. from the forecaster's information set); the time period corresponding to a realization and forecasts at time t is t-h, i.e. this series has to start H periods prior to y and yhat and to end one period earlier
#'
#' @return a list with two elements;
#' * `resultsMZ` main results of the test, a named one matrix with one row containing five numbers:
#'     * `Stat` value of the MZ test statistic
#'     * `90%/95%/99%` 0.9/0.95/0.99 quantile of bootstrap distribution
#'     * `p-value` p-value of the test
#' * `tables` additional information for each horizon-quantile combination, a list with 3+A data frames:
#'    * `cont` individual contributions to test statistic
#'    * `alphahat` MZ regression intercepts
#'    * `betahat` MZ regression slopes
#'    * data frames 4 to 3+A.: regression coefficients of additional regressors
#'
#'
#' @examples
#' @export
#' @importFrom quantreg rq




MZ_test_augmented <- function(y, tau, H, yhat,B,l,adreg){

  T <- length(y)
  A <- ncol(adreg)

  # generate a list with matrices (one for each additional regressors) of the same format as yhat, where each entry of the matrix contains the respective regressor for horizon hcount and time t

  adregmat <- vector(mode = "list", length = A)

  for (a in 1:A){
    mat <- matrix(NA,nrow=T,ncol=H)

    for (hcount in 1:H){
      mat[,hcount] <- adreg[(H-hcount+1):(H-hcount+T),a]
    }

    adregmat[[a]] <- mat
  }


  # calculate test statistic and save tables and mhat



  result <- MZ_test_statistic_augmented(y=y, tau=tau, H=H, yhat=yhat,adregmat=adregmat)

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

    # draw realizations, forecasts and additional regressors
    ystar=y[I]
    yhatstar=lapply(1:length(tau),function(x){as.matrix(yhat[[x]][I,])})
    adregmatstar=lapply(1:A,function(x){as.matrix(adregmat[[x]][I,])})

    # compute bootstrap test statistic

    bootstats[boot] <- Pstar*sum((MZ_test_statistic_augmented(y=ystar, tau=tau, H=H, yhat=yhatstar,adregmat=adregmatstar)[[3]] - m)^2)

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
