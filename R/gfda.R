#' GFDA model
#' 
#' This function is an wrapper to fitting a Geostatistical Functional Data 
#' Analysis (GFDA) model to spatio-temporal data.
#'
#' @param training.set A numeric matrix with the response variable in the first column, 
#'              time values in the second column and X, Y coordinates in the third
#'              column. Corresponds to the training data.
#' @param prediction.set Same as training.set, but corresponding to points to be predicted for 
#'              MSPE calculation.
#' @param ssensors Integer corresponding to how many static sensors in the dataset.
#' @param rsensors Integer corresponding to how many roving sensors in the dataset.
#' @param J Number of eigenfunctions in spatio-temporal covariance, defaults to 2.
#' @param spline.df Number of spline basis for the deterministic spline component (see \code{\link{bs}} function).
#' @param fpca.df Number of spline basis for the stochastic spline component (see \code{\link{tfpca}} function).
#' @param verbose Wheter \code{gfda} gives detailed information on optimization step, defaults to \code{TRUE}.
#'
#' @export
#' @return List of four items 
#'   \item{beta.est}{Coefficient estimates}
#'   \item{MSPEKrig}{Mean squared Kriging Prediction error}
#'   \item{spatCov}{Theta parameters for stochastic term}
#'   \item{sigma2s}{Static sensor variability}
#'   \item{sigma2r}{Roving sensor variability}
#'
#' @examples
#' ## require(fda)
#' ## Y <- CanadianWeather$dailyAv
#' ## XY <- CanadianWeather$coordinates
#' # static sensor simulation example:
#' set.seed(1)
#' X <- cbind(rnorm(200), rep(1:50,4), rep(c(2,2,4,4), each = 50), rep(c(2,4,2,4), each = 50))
#' Xpred <- cbind(rnorm(50), 1:50, 3, 3)
#' results <- gfda(X, Xpred, ssensors = 4)
#' 
#' @references
#'  \url{http://www.google.com}
#'
#' @seealso \code{\link{tfpca}}
#' @keywords Spatial Statistics
#' @keywords Functional Data Analysis
gfda <- function(training.set, prediction.set, subtfpca = NULL, ssensors = 6, rsensors = 2, 
                 J = 2, spline.df = NULL, fpca.df = 20, verbose = TRUE){
  
  # NEEDS TO TAKE GENERIC TIME!!
  # ssensors might be useless! They do not have the same # of observations.
  
  nTR <- dim(training.set)[1]  
  nTS <- dim(prediction.set)[1]
  nT <- nTR + nTS
  
  #!# Gui: XTR <- cbind(1, x, y, S(t))
  XTR <- cbind(1, training.set[ ,3:4], splines::bs(training.set[ ,2], df = spline.df))
  YTR <- training.set[ ,1] #!# Gui: YTR=training.set$Leq
  
  Step1 <- .lm.fit(XTR, YTR) # Fastest least squares
  training.set.detrended <- training.set
  training.set.detrended[,1] <- Step1$residuals
  
  timerange <- (diff(range(training.set[ ,2])))
  if(is.null(subtfpca)){
    NoiStDe <- matrix(training.set.detrended[, 1], ncol = ssensors) #!# IMPORTANT
    t.pred <- t.fit <- (unique(training.set[ ,2]) - min(training.set[ ,2])) / timerange #  +2 was here
  } else {
    NoiStDe <- matrix(training.set.detrended[subtfpca, 1], ncol = ssensors-rsensors) #!# IMPORTANT
    t.pred <- t.fit <- (unique(training.set[subtfpca,2]) - min(training.set[subtfpca,2]) )/timerange
  }
  nSt <- length(t.fit)
  Step2 <- tfpca(NoiStDe, J, t.fit, t.pred, fpca.df) 
  lamb.est <- Step2$values
  Phi.est <- Phi.TR <- Phi.TE <- Step2$vectors
  
  ##################################
  # Step3:   Spatial Parameters    #
  ##################################
  
  DTR <- matrix(0, nrow = nTR, ncol = nTR)        
  for(i in 1:nTR) {
    DTR[i, ] <- sqrt((training.set[i,3] - training.set[,3])^2 + (training.set[i,4] - training.set[,4])^2)
  }
  
  Dmax = max(DTR)
  Vmax = var(training.set.detrended[,1])
  theta0 = c(rep(Dmax/2, J), Vmax/2, Vmax/2) # theta_{1:J}, sigma_S, sigma_R
  
  if(verbose) cat("Optimizing... \n")
  
  # Constraint: ui %*% theta - ci >= 0 
  # Check: http://eigen.tuxfamily.org/dox/AsciiQuickReference.txt
  prof.max <- constrOptim(theta0, logProfileCpp, grad = dlogProfileCpp,
                          ui = cbind(c(1,0,0,0,-1,0,0,0),
                                     c(0,1,0,0,0,-1,0,0), 
                                     c(0,0,1,0,0,0,-1,0), 
                                     c(0,0,0,1,0,0,0,-10)), 
                          ci = c(0.001, 0.001, 0.001, 0.001, 
                                 -Dmax, -Dmax, -Vmax, -Vmax), 
                          DTR = DTR, Y = training.set[,1], XTR = XTR, 
                          PhiTime = Phi.est[training.set[,2],], LambEst = lamb.est)
  ### CORRECT THIS: Phi.est[training.set[,2],]
  ### TRY:  Phi.est[match(((prediction.set[ ,2] - min(training.set[ ,2])) / timerange), t.pred), ]
  theta <- prof.max$par
  
  if(verbose) cat("Optimization done. \n")
  
  # PhiTime <- Phi.est[training.set[,2],] # WRITE EVALUATOR HERE!
  PhiTime <- Phi.est[match(((training.set[ ,2] - min(training.set[ ,2])) / timerange), t.fit), ]
  psi.cov <- matrix(0, nrow = nTR, ncol = nTR)
  Gamma1 <- lamb.est[1]*exp(-DTR/theta[1])
  Gamma2 <- lamb.est[2]*exp(-DTR/theta[2])
  for (i in 1:nTR){
    psi.cov[i,] <- Gamma1[i,]*PhiTime[i,1]*PhiTime[,1]+Gamma2[i,]*PhiTime[i,2]*PhiTime[,2]
  }
  diag(psi.cov) <- diag(psi.cov)+theta[3]
  
  beta.est <- solve(crossprod(XTR, solve(psi.cov, XTR)), crossprod(XTR, solve(psi.cov, YTR)))
  resid <- YTR-XTR%*%beta.est
  
  # Starts Kriging
  
  DKrig <- psi.krig <- matrix(0, nrow = nTR, ncol = nTS)        
  
  for(i in 1:nTR) {
    DKrig[i,] <- sqrt((training.set[i,3] - prediction.set[,3])^2 + (training.set[i,4] - prediction.set[,4])^2)
  }
  
  # PhiTimeTE <- Phi.est[prediction.set[,2],] # WRITE EVALUATOR HERE!
  PhiTimeTE <- Phi.est[match(((prediction.set[ ,2] - min(training.set[ ,2])) / timerange), t.pred), ]
  Gamma1TE <- lamb.est[1]*exp(-DKrig/theta[1])
  Gamma2TE <- lamb.est[2]*exp(-DKrig/theta[2])
  for (i in 1:nTR){
    psi.krig[i,] <- Gamma1TE[i,]*PhiTime[i,1]*PhiTimeTE[,1]+Gamma2TE[i,]*PhiTime[i,2]*PhiTimeTE[,2]
  }
  XTE <- cbind(1, prediction.set[,3:4], splines::bs(prediction.set[ ,2], df = spline.df))
  YKrig <- XTE%*%beta.est+crossprod(psi.krig, solve(psi.cov, resid))
  MSPEKrig <- mean((prediction.set[,1]-YKrig)^2)
  
  beta.est <- as.numeric(beta.est)
  names(beta.est) <- c("b0", "bx", "by", paste0("sp", seq_len(length(beta.est)-3))) # 3 basis function, time using defaults
  ret <- list(beta.est = beta.est, MSPEKrig = MSPEKrig, spatCov = theta)
  return(ret)
}