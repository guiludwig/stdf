#' GFDA model
#' 
#' This function is an wrapper to fitting a Geostatistical Functional Data 
#' Analysis (GFDA) model to spatio-temporal data.
#'
#' @param NoiTR A spatio-temporal data model with vertical alignment (TODO: Generalize)
#' @param NoiTE A spatio-temporal data model with vertical alignment, prediction set (TODO: Generalize)
#' @param ssensors How many sensors in the dataset (To be removed)
#' @param verbose Wheter \code{gfda} gives detailed information, defaults to \code{TRUE}
#' @param J Number of eigenfunctions in spatio-temporal covariance, defaults to 2 (TODO: implement)
#'
#' @export
#' @return List of four items 
#'   \item{beta.est}{Coefficient estimates}
#'   \item{fitted.values}{TODO}
#'   \item{MSPEKrig}{Mean squared Kriging Prediction error}
#'   \item{Cov.matrix}{TODO}
#'
#' @examples
#' require(fda)
#' Y <- CanadianWeather$dailyAv
#' XY <- CanadianWeather$coordinates
#' mean(rnorm(20))
#' 
#' @references
#'  \url{http://www.google.com}
#'
#' @seealso \code{\link{median}},
#'   \code{\link{mean}}
#' @keywords Spatial Statistics
#' @keywords Functional Data Analysis
gfda <- function(NoiTR, NoiTE, subtfpca = NULL, ssensors = 6, 
                 excsensors = 2, verbose = TRUE, tbas = 20){
  
  # NEEDS TO TAKE GENERIC TIME!!
  # ssensors might be useless! They do not have the same # of observations.
  
  nTR <- dim(NoiTR)[1]  
  nTS <- dim(NoiTE)[1]
  nT <- nTR + nTS
  
  #!# Can insert knots here
  #!# Gui: XTR <- cbind(1, x, y, S(t))
  XTR <- cbind(1, NoiTR[ ,3:4], splines::bs(NoiTR[ ,2]))
  YTR <- NoiTR[ ,1] #!# Gui: YTR=NoiTR$Leq
  
  Step1 <- lm(YTR~XTR-1)
  NoiTRDe <- NoiTR
  NoiTRDe[,1] <- Step1$res
  
  timerange <- (diff(range(NoiTR[ ,2])))
  J <- 2
  if(is.null(subtfpca)){
    NoiStDe <- matrix(NoiTRDe[, 1], ncol=ssensors) #!# IMPORTANT
    t.pred <- t.fit <- (unique(NoiTR[ ,2]) - min(NoiTR[ ,2])) / timerange #  +2 was here
  } else {
    NoiStDe <- matrix(NoiTRDe[subtfpca, 1], ncol=ssensors-excsensors) #!# IMPORTANT
    t.pred <- t.fit <- (unique(NoiTR[subtfpca,2]) - min(NoiTR[subtfpca,2]) )/timerange
  }
  nSt <- length(t.fit)
  Step2 <- tfpca(NoiStDe, J, t.fit, t.pred, tbas) 
  lamb.est <- Step2$values
  Phi.est <- Phi.TR <- Phi.TE <- Step2$vectors
  
  ##################################
  # Step3:   Spatial Parameters    #
  ##################################
  
  DTR <- matrix(0, nrow = nTR, ncol = nTR)        
  for(i in 1:nTR) {
    DTR[i, ] <- sqrt((NoiTR[i,3] - NoiTR[,3])^2 + (NoiTR[i,4] - NoiTR[,4])^2)
  }
  
  Dmax = max(DTR)
  Vmax = var(NoiTR[,1])
  theta0=c(Dmax/2,Dmax/2,Vmax/2)
  
  if(verbose) cat("Optimizing... \n")
  
  # Constraint: ui %*% theta - ci >= 0 
  # system.time({log.profile(theta0, DTR = DTR, NoiTR = NoiTR, XTR = XTR, Phi.est = Phi.est, lamb.est = lamb.est)})
  prof.max <- constrOptim(theta0, logProfileCpp, grad = dlogProfileCpp,
                          ui = cbind(c(1,0,0,-1,0,0),c(0,1,0,0,-1,0), c(0,0,1,0,0,-1)), 
                          ci = c(0.001, 0.001, 0.001, -Dmax, -Dmax, -Vmax), 
                          DTR = DTR, Y = NoiTR[,1], XTR = XTR, 
                          PhiTime = Phi.est[NoiTR[,2],], LambEst = lamb.est)
  
  theta <- prof.max$par
  
  if(verbose) cat("Optimization done. \n")
  
  # PhiTime <- Phi.est[NoiTR[,2],] # WRITE EVALUATOR HERE!
  PhiTime <- Phi.est[match(((NoiTR[ ,2] - min(NoiTR[ ,2])) / timerange), t.fit), ]
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
    DKrig[i,] <- sqrt((NoiTR[i,3] - NoiTE[,3])^2 + (NoiTR[i,4] - NoiTE[,4])^2)
  }
  
  # PhiTimeTE <- Phi.est[NoiTE[,2],] # WRITE EVALUATOR HERE!
  PhiTimeTE <- Phi.est[match(((NoiTE[ ,2] - min(NoiTR[ ,2])) / timerange), t.pred), ]
  Gamma1TE <- lamb.est[1]*exp(-DKrig/theta[1])
  Gamma2TE <- lamb.est[2]*exp(-DKrig/theta[2])
  for (i in 1:nTR){
    psi.krig[i,] <- Gamma1TE[i,]*PhiTime[i,1]*PhiTimeTE[,1]+Gamma2TE[i,]*PhiTime[i,2]*PhiTimeTE[,2]
  }
  XTE <- cbind(1, NoiTE[,3:4], splines::bs(NoiTE[ ,2]))
  YKrig <- XTE%*%beta.est+crossprod(psi.krig, solve(psi.cov, resid))
  MSPEKrig <- mean((NoiTE[,1]-YKrig)^2)
  
  beta.est <- as.numeric(beta.est)
  names(beta.est) <- c("b0", "bx", "by", "sp1", "sp2", "sp3") # 3 basis function, time using defaults
  ret <- list(beta.est = beta.est, MSPEKrig = MSPEKrig, spatCov = theta)
  return(ret)
}