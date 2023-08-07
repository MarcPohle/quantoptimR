#' @title MZ_test_statistic
#'
#' @description computes the MZ test statistic from an evaluation sample
#'
#'
#' @param y the realizations, a vector
#' @param tau vector of quantile levels, for which forecasts are provided
#' @param H maximum forecast horizon
#' @param yhat the forecasts, a list of length length(tau), which contains length(y) x H matrices with forecasts for a specific quantile for horizons 1 to H corresponding to the realizations (e.g. row 1 contains 1- up to H-period-ahead forecasts for time 1, where the first value in the realizations vector is also for time 1)
#'
#'
#' @return a list with three elements;
#' * `Uhat_MZ` value of the MZ test statistic
#' * additional information for each horizon-quantile combination, a list with three data frames:
#'    * `cont` individual contributions to test statistic
#'    * `alphahat` MZ regression intercepts
#'    * `betahat` MZ regression slopes
#' * `mhat` vector of the individual mhats (see paper for their definition)
#'
#'
#' @examples
#' @export
#' @importFrom quantreg rq






MZ_test_statistic <- function(y, tau, H, yhat){

  T <- length(y)

  #generate empty vector to save MZ coefficients
  coef=matrix(NA,2*length(tau)*H,1)

  # Save quantile regression coefficients into coef

  for (taucount in 1:length(tau)){

    for (hcount in 1:H){

      coef[((taucount-1)*(2*H)+2*hcount-1):((taucount-1)*(2*H)+2*hcount),]=coef(rq(y~yhat[[taucount]][,hcount],tau=tau[taucount]))

    }
  }

  # Obtain MZ test statistic

  m <- coef-rep(c(0,1),H*length(tau)) # First subtract 1 from the slope coefficient
  Uave=sum((sqrt(T)*m)^2)


  # generate tables with individual contributions to test statistic, intercepts and slopes

  # Extract individual contributions to test stat
  Uaveind=sapply(1:(H*length(tau)),function(x){sum((sqrt(T)*m[(2*x-1):(2*x)])^2)})
  Uaveind=matrix(Uaveind,nrow=H)
  rownames(Uaveind)=paste("h=",1:H,sep="")
  colnames(Uaveind)=paste("tau=",tau,sep="")
  Uaveind <- as.data.frame(Uaveind)

  #Extract intercepts and slopes

  alpha <- coef[seq(from=1,to=length(coef),by=2)]
  alpha=matrix(alpha,nrow=H)
  rownames(alpha)=paste("h=",1:H,sep="")
  colnames(alpha)=paste("tau=",tau,sep="")
  alpha <- as.data.frame(alpha)

  beta <- coef[seq(from=2,to=length(coef),by=2)]
  beta=matrix(beta,nrow=H)
  rownames(beta)=paste("h=",1:H,sep="")
  colnames(beta)=paste("tau=",tau,sep="")
  beta <- as.data.frame(beta)

  list(Uhat_MZ = Uave, list(cont = Uaveind, alphahat = alpha, betahat=beta), mhat=m)
}

