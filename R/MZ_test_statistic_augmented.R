#' @title MZ_test_statistic_augmented
#'
#' @description computes the test statistic of the augmented MZ test from an evaluation sample
#'
#'
#' @param y the realizations, a vector
#' @param tau vector of quantile levels, for which forecasts are provided
#' @param H maximum forecast horizon
#' @param yhat the forecasts, a list of length length(tau), which contains length(y) x H matrices with forecasts for a specific quantile for horizons 1 to H corresponding to the realizations (e.g. row 1 contains 1- up to H-period-ahead forecasts for time 1, where the first value in the realizations vector is also for time 1)
#' @param adregmat additional regressors, a list of length A (the number of additional regressors), with length(y) x H matrices with additional regressors for horizons 1 to H corresponding to the realizations and forecasts; i.e. column hcount contains in every row the corresponding regressors for the MZ regression using the same row in yhat and y; this matrix structure is necessary in the bootstrap, where this function is used repeatedly
#'
#' @return a list with three elements;
#' * `Uhat_MZ` value of the MZ test statistic
#' * additional information for each horizon-quantile combination, a list with 3+A data frames:
#'    * `cont` individual contributions to test statistic
#'    * `alphahat` MZ regression intercepts
#'    * `betahat` MZ regression slopes
#'    * data frames 4 to 3+A.: regression coefficients of additional regressors
#' * `mhat` vector of the individual mhats (see paper for their definition)
#'
#'
#' @examples
#' @export
#' @importFrom quantreg rq






MZ_test_statistic_augmented <- function(y, tau, H, yhat, adregmat){

  T <- length(y)
  A <- ncol(adregmat[[1]])/H
  k <- 2+A # number of regression coefficients in augmented MZ regression


  # generate empty vector to save MZ coefficients
  coef=matrix(NA,k*length(tau)*H,1)

  # Save quantile regression coefficients into coef

  for (taucount in 1:length(tau)){

    for (hcount in 1:H){



      # generate matrix with additional regressors for the current horizon hcount
      adreg_hcount <- matrix(NA,nrow=T,ncol=A)
      for (a in 1:A){
        adreg_hcount[,a] <- adregmat[[a]][,hcount]
      }

      # add the forecasts

      MZ_regressors <- as.data.frame(cbind(yhat[[taucount]][,hcount],adreg_hcount))

      # generate full data matrix for MZ quantile regression by adding realizations

      data <- cbind(data.frame(y=y), MZ_regressors)


      coef[((taucount-1)*(k*H)+k*(hcount-1) +1):((taucount-1)*(k*H)+k*hcount),]=coef(rq(y~.,data=data,tau=tau[taucount]))

    }
  }

  # Obtain MZ test statistic

  m <- coef-rep(c(0,1,rep(0,times=A)),H*length(tau)) # First subtract 1 from the slope coefficient

  Uave=sum((sqrt(T)*m)^2)


  # generate tables with individual contributions to test statistic, intercepts and slopes

  # Extract individual contributions to test stat
  Uaveind=sapply(1:(H*length(tau)),function(x){sum((sqrt(T)*m[(k*(x-1)+1):(k*x)])^2)})
  Uaveind=matrix(Uaveind,nrow=H)
  rownames(Uaveind)=paste("h=",1:H,sep="")
  colnames(Uaveind)=paste("tau=",tau,sep="")
  Uaveind <- as.data.frame(Uaveind)

  #Extract intercepts, regression coefficients of forecasts and of additional regressors

  alpha <- coef[seq(from=1,to=length(coef),by=k)]
  alpha=matrix(alpha,nrow=H)
  rownames(alpha)=paste("h=",1:H,sep="")
  colnames(alpha)=paste("tau=",tau,sep="")
  alpha <- as.data.frame(alpha)

  beta <- coef[seq(from=2,to=length(coef),by=k)]
  beta=matrix(beta,nrow=H)
  rownames(beta)=paste("h=",1:H,sep="")
  colnames(beta)=paste("tau=",tau,sep="")
  beta <- as.data.frame(beta)


  gammas <- vector(mode = "list", length = A)

  for (a in 1:A){
    gamma <- coef[seq(from=2+a,to=length(coef),by=k)]
    gamma=matrix(gamma,nrow=H)
    rownames(gamma)=paste("h=",1:H,sep="")
    colnames(gamma)=paste("tau=",tau,sep="")
    gamma <- as.data.frame(gamma)
    gammas[[a]] <- gamma
  }

 list(Uhat_MZ = Uave, c(list(cont = Uaveind, alphahat = alpha, betahat=beta),gammas), mhat=m)
}

