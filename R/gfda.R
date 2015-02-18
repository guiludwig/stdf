#' GFDA model
#' 
#' This function is an wrapper to fitting a Geostatistical Functional Data 
#' Analysis (GFDA) model to spatio-temporal data.
#'
#' @param training.set A numeric matrix with the response variable in the first column, 
#'                     time values in the second column and X, Y coordinates in the third
#'                     column. Corresponds to the training data.
#' @param prediction.set Same as training.set, but corresponding to points to be predicted for 
#'                       MSPE calculation.
#' @param subtfpca Logical vector with the same length as the number of rows in training.set; tells
#'                 gfda which lines should be included in the calculation of the temporal dependency
#'                 component. Defaults to NULL, which includes all lines.
#' @param ssensors Integer corresponding to how many static sensors in the dataset.
#' @param L Number of eigenfunctions in spatio-temporal covariance, defaults to 2.
#' @param spline.df Number of spline basis for the deterministic spline component (see 
#'                  \code{\link{bs}} function).
#' @param fpca.df Number of spline basis for the stochastic spline component (see \code{\link{tfpca}} 
#'                function).
#' @param homogeneous Whether the variance of static and roving sensors is assumed to be the same or not. 
#'                    Defaults to not.
#' @param verbose Whether \code{gfda} prints information about the optimization step, defaults to 
#'                \code{TRUE}.
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
#' ## Y <- CanadianWeather$dailyAv
#' ## XY <- CanadianWeather$coordinates
#' # =sensor simulation example:
#' library(fda)
#' set.seed(1)
#' static <- cbind(rep(0, 200), rep(1:50,4), rep(c(2,2,4,4), each = 50), rep(c(2,4,2,4), each = 50))
#' roving <- cbind(rep(0, 50), 1:50, seq(2,4, length=50), seq(2,4, length=50))
#' prediction <- cbind(rep(0, 50), 1:50, 3, 3)
#' complete <- rbind(static, roving, prediction)
#' D <- matrix(0, 300, 300)
#' for(i in 1:300) D[i,] <- sqrt((complete[i,3] - complete[,3])^2 + (complete[i,4] - complete[,4])^2)
#' S <- matrix(0, 300, 300)
#' Phi <- cbind(sin(2*pi*complete[,2]/50), cos(2*pi*complete[,2]/50), 2*sin(4*pi*complete[,2]/50))
#' sigmaR <- 2
#' sigmaS <- 1
#' theta <- 1:3
#' for(L in 1:3) S <- S + exp(-D/theta[L])*outer(Phi[,L], Phi[,L])
#' diag(S) <- diag(S) + sigmaS
#' diag(S)[201:250] <- diag(S)[201:250] - sigmaS + sigmaR
#' complete[,1] <- 2 + 2*complete[,3] + 2*complete[,4] + 2*cos(4*pi*complete[,2]/50) + 
#'                 t(chol(S)) %*% rnorm(300)
#' static[,1] <- complete[1:200,1]
#' roving[,1] <- complete[201:250,1]
#' prediction[,1] <- complete[251:300,1]
#' 
#' results <- gfda(static, prediction, ssensors = 4, L = 3)
#' resultsH <- gfda(static, prediction, ssensors = 4, L = 3, homogeneous = TRUE)
#' resultsPlusRoving <- gfda(rbind(static, roving), prediction, 
#'                           subtfpca = c(rep(TRUE,200), rep(FALSE,50)), 
#'                           ssensors = 4, L = 3)
#' # Wrong model:                         
#' resultsPlusRovingH <- gfda(rbind(static, roving), prediction, 
#'                           subtfpca = c(rep(TRUE,200), rep(FALSE,50)), 
#'                           ssensors = 4, L = 3, homogeneous = TRUE)
#' 
#' @references
#'  \url{http://www.google.com}
#'
#' @seealso \code{\link{tfpca}}
#' @keywords Spatial Statistics
#' @keywords Functional Data Analysis
gfda <- function(training.set, prediction.set, subtfpca = NULL, ssensors = 6, 
                 L = 2, spline.df = NULL, fpca.df = 10, homogeneous = FALSE, 
                 verbose = TRUE, ...){
  
  if(!is.null(subtfpca) && sum(subtfpca) == length(subtfpca)) message("Only static sensors, consider setting homogeneous = TRUE or subtfpca = NULL")
  
  nTR <- nrow(training.set)
  nTS <- nrow(prediction.set)
  nT <- nTR + nTS
  
  #!# XTR <- cbind(1, x, y, S(t))
  XTR <- cbind(1, training.set[ ,3:4], splines::bs(training.set[ ,2], df = spline.df))
  YTR <- training.set[ ,1] #!# YTR <- training.set$Leq
  
  Step1 <- .lm.fit(XTR, YTR) # Fastest least squares
  training.set.detrended <- training.set
  training.set.detrended[,1] <- Step1$residuals
  
  if(is.null(subtfpca)){
    NoiStDe <- matrix(training.set.detrended[, 1], ncol = ssensors) #!# IMPORTANT
    t.pred <- unique(prediction.set[ ,2])
    t.fit <- unique(training.set[ ,2])
  } else {
    NoiStDe <- matrix(training.set.detrended[subtfpca, 1], ncol = ssensors) #!# IMPORTANT
    t.pred <- unique(prediction.set[ ,2])
    t.fit <- unique(training.set[subtfpca,2])
  }
  nSt <- length(t.fit)
  Step2 <- tfpca(NoiStDe, L, t.fit, t.pred, fpca.df, ...) 
  lamb.est <- Step2$values
  Phi.est <- Phi.TR <- Phi.TE <- Step2$vectors
  
  ##################################
  # Step3:   Spatial Parameters    #
  ##################################
  
  DTR <- matrix(0, nrow = nTR, ncol = nTR)  
  for(i in 1:nTR) {
    DTR[i, ] <- sqrt((training.set[i,3] - training.set[,3])^2 + (training.set[i,4] - training.set[,4])^2)
  }
  
  Dmax <- max(DTR)
  Vmax <- var(training.set.detrended[,1])
  if(homogeneous){
    theta0 <- c(rep(Dmax/2, L), Vmax/2) # theta_{1:L}, sigma
    UI <- rbind(diag(L+1), -1*diag(L+1))
    CI <- c(rep(0.001, L+1), rep(-Dmax, L), -Vmax)
  } else {
    theta0 <- c(rep(Dmax/2, L), Vmax/2, Vmax/2) # theta_{1:L}, sigma_S, sigma_R
    UI <- rbind(diag(L+2), -1*diag(L+2))
    CI <- c(rep(0.001, L+2), rep(-Dmax, L), -Vmax, -Vmax)
  }
  
  if(verbose) cat("Optimizing... \n")
  
  # Constraint: ui %*% theta - ci >= 0 
  # Check: http://eigen.tuxfamily.org/dox/AsciiQuickReference.txt
  if(is.null(subtfpca)) {
    subsetStatic <- rep(1, nTR)
  } else {
    subsetStatic <- as.numeric(subtfpca)
  }
  
  if(!homogeneous){
    prof.max <- constrOptim(theta0, logProfileCpp, grad = dlogProfileCpp, # grad = dlogProfileCpp,
                            ui = UI, ci = CI, # Constraints
                            DTR = DTR, Y = YTR, XTR = XTR, 
                            subsetStatic = subsetStatic,
                            PhiTime = Phi.est[match(training.set[ ,2], t.fit),], 
                            LambEst = lamb.est)
  } else {
    prof.max <- constrOptim(theta0, logProfileCppH, grad = dlogProfileCppH,
                            ui = UI, ci = CI, # Constraints
                            DTR = DTR, Y = YTR, XTR = XTR, 
                            PhiTime = Phi.est[match(training.set[ ,2], t.fit),], 
                            LambEst = lamb.est)
  }
  
  theta <- prof.max$par
  if(is.null(subtfpca) & !homogeneous) {
    theta <- theta[-length(theta)] 
  }
  
  if(verbose) cat("Optimization done. \n")
  
  PhiTime <- Phi.est[match(training.set[ ,2], t.fit),]
  psi.cov <- matrix(0, nrow = nTR, ncol = nTR)
  for(j in 1:L){
    psi.cov <- psi.cov + lamb.est[j]*exp(-DTR/theta[j])*outer(PhiTime[,j],PhiTime[,j])
  }
  diag(psi.cov) <- diag(psi.cov) + theta[L+1]
  if (!homogeneous & sum(subsetStatic) < length(subsetStatic)){
    diag(psi.cov) <- diag(psi.cov) + (theta[L+2] - theta[L+1])*(1-subsetStatic)
  } 
  
  beta.est <- solve(crossprod(XTR, solve(psi.cov, XTR)), crossprod(XTR, solve(psi.cov, YTR)))
  resid <- YTR-XTR%*%beta.est
  
  # Starts Kriging
  
  DKrig <- psi.krig <- matrix(0, nrow = nTR, ncol = nTS)        
  
  for(i in 1:nTR) {
    DKrig[i,] <- sqrt((training.set[i,3] - prediction.set[,3])^2 + (training.set[i,4] - prediction.set[,4])^2)
  }
  
  PhiTimeTE <- Phi.est[match(prediction.set[ ,2], t.pred), ]
  for(j in 1:L){
    psi.krig <- psi.krig + lamb.est[j]*exp(-DKrig/theta[j])*outer(PhiTime[,j],PhiTimeTE[,j])
  }
  
  XTE <- cbind(1, prediction.set[,3:4], splines::bs(prediction.set[ ,2], df = spline.df))
  YKrig <- XTE%*%beta.est+crossprod(psi.krig, solve(psi.cov, resid))
  MSPEKrig <- mean((prediction.set[,1]-YKrig)^2)
  
  beta.est <- as.numeric(beta.est)
  names(beta.est) <- c("b0", "bx", "by", paste0("sp", seq_len(length(beta.est)-3))) # 3 basis function, time using defaults
  ret <- list(beta.est = beta.est, MSPEKrig = MSPEKrig, spatCov = theta)
  return(ret)
}