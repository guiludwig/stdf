#' STDF model
#' 
#' This function is an wrapper to fitting a Spatio-Temporal Data 
#' Fusion (STDF) model to spatio-temporal data.
#'
#' @param training.set A numeric matrix with the response variable in the 
#'                     first column, time values in the second column and x, y 
#'                     spatial coordinates in the third column. Corresponds to 
#'                     the training data.
#' @param subtfpca Logical vector with the same length as the number of rows in 
#'                 training.set; tells stdf which lines should be included in 
#'                 the calculation of the temporal dependency component. 
#'                 Defaults to NULL, which includes all lines.
#' @param ssensors Integer corresponding to how many static sensors in the 
#'                 dataset.
#' @param L Number of eigenfunctions in spatio-temporal covariance, defaults 
#'                 to 2.
#' @param spline.df Number of spline basis for the deterministic spline 
#'                  component (see \code{\link{bs}} function). Defaults to 
#'                  NULL, which uses 3 spline basis.
#' @param zeta Smoothness parameter for the fpca function. Defaults 
#'                    to 0 (no smoothing). See \code{\link{tfpca}} function.
#' @param fpca.df Number of spline basis for the stochastic spline component 
#'                (see \code{\link{tfpca}} function). Defaults to 20. 
#' @param homogeneous Whether the variance of static and roving sensors is 
#'                    assumed to be the same or not. Defaults to FALSE.
#' @param verbose Whether \code{stdf} prints a message when the optimization 
#'                step is done. Defaults to \code{TRUE}.
#' @param method Either "L-BFGS-B" or "Nelder-Mead", which are passed to the 
#'               optimization function \code{\link{constrOptim}}. 
#' @param ... further arguments passed to or from other methods.
#'
#' @export
#' @return \code{stdf} returns an object of \code{\link{class}} "stdf" containing 
#' at least the following components
#'   \item{beta.est}{Coefficient estimates}
#'   \item{MSPEKrig}{Mean squared Kriging Prediction error}
#'   \item{spatCov}{Theta parameters for stochastic term}
#'   \item{sigma2s}{Static sensor variance estimate}
#'   \item{sigma2r}{Roving sensor variance estimate}
#'
#' @examples
#' # sensor simulation example, placeholder for unit tests but also shows
#' # options available for optimization:
#' set.seed(1)
#' static <- cbind(Y = rep(0, 200), 
#'                 t = rep(1:50,4), 
#'                 sx = rep(c(2,2,4,4), each = 50), 
#'                 sy = rep(c(2,4,2,4), each = 50))
#' roving <- cbind(Y = rep(0, 50), 
#'                 t = 1:50, 
#'                 sx = seq(2,4, length=50), 
#'                 sy = seq(2,4, length=50))
#' prediction <- cbind(Y = rep(0, 50), t = 1:50, sx = 3, sy = 3)
#' complete <- rbind(static, roving, prediction)
#' D <- matrix(0, 300, 300)
#' for(i in 1:300) D[i,] <- sqrt((complete[i,3] - complete[,3])^2 + 
#'                               (complete[i,4] - complete[,4])^2)
#' S <- matrix(0, 300, 300)
#' Phi <- cbind(sin(2*pi*complete[,2]/50), 
#'              cos(2*pi*complete[,2]/50), 
#'              2*sin(4*pi*complete[,2]/50))
#' sigmaR <- 2
#' sigmaS <- 1
#' theta <- 1:3
#' for(L in 1:3) S <- S + exp(-D/theta[L])*outer(Phi[,L], Phi[,L])
#' diag(S) <- diag(S) + sigmaS
#' diag(S)[201:250] <- diag(S)[201:250] - sigmaS + sigmaR
#' complete[,1] <- 2 + 2*complete[,3] + 2*complete[,4] + 
#'                 2*cos(4*pi*complete[,2]/50) + 
#'                 t(chol(S)) %*% rnorm(300)
#' static[,1] <- complete[1:200,1]
#' roving[,1] <- complete[201:250,1]
#' prediction[,1] <- complete[251:300,1]
#' 
#' results <- stdf(static, ssensors = 4, L = 3)
#' resultsNM <- stdf(static, ssensors = 4, L = 3, 
#'                   method = "Nelder-Mead")
#' # Compare results of using "L-BFGS-B" and "Nelder-Mead";
#' # Should be almost the same, bar numberical accuracy issues
#' results
#' resultsNM
#' # How to predict, and some numerical accuracy difference
#' test.A <- predict(results, prediction)
#' test.B <- predict(resultsNM, prediction)
#' all.equal(test.A$fitted, test.B$fitted)
#' # Compare with fitting homogeneous case
#' resultsH <- stdf(static, ssensors = 4, L = 3, 
#'                  homogeneous = TRUE)
#' resultsHNM <- stdf(static, ssensors = 4, L = 3, 
#'                    homogeneous = TRUE, method = "Nelder-Mead")
#'                    
#' # Nelder-Mead is more time-consuming, for example in this case:
#' system.time({resultsPlusRoving <- stdf(rbind(static, roving), 
#'                           subtfpca = c(rep(TRUE,200), rep(FALSE,50)), 
#'                           ssensors = 4, L = 3)})
#' system.time({resultsPlusRovingNM <- stdf(rbind(static, roving), 
#'                             subtfpca = c(rep(TRUE,200), rep(FALSE,50)), 
#'                             ssensors = 4, L = 3, method = "Nelder-Mead")})
#' # Results should be almost the same, bar numerical accuracy issues
#' resultsPlusRoving
#' resultsPlusRovingNM
#' # Forces homogeneous model:                         
#' resultsPlusRovingH <- stdf(rbind(static, roving), 
#'                           subtfpca = c(rep(TRUE,200), rep(FALSE,50)), 
#'                           ssensors = 4, L = 3, homogeneous = TRUE)
#' 
#' # Example with dataset from fda package; this dataset is
#' # huge and takes a long time to run.
#' \dontrun{
#' library(fda)
#' n <- 35 # stations
#' Y <- as.numeric(CanadianWeather$dailyAv[ , ,"Temperature.C"])
#' XY <- cbind(Y, 
#'             rep(1:365, n),  
#'             rep(CanadianWeather$coordinates[, "N.latitude"], each = 365), 
#'             rep(CanadianWeather$coordinates[, "W.longitude"], each = 365))
#' fakePred <- XY[1:3,]
#' model <- stdf(XY, ssensors = n, homogeneous = TRUE,
#'               L = 2, zeta = 1e5)
#' # Creating a hazard map requires the ggplot2 package
#' # TODO
#' }
#'   
#' @author Guilherme Ludwig and Tingjin Chu
#'
#' @references
#' 
#'   Chu, T., Zhu, J. and Wang, H. (2014) On Semiparametric Inference of Geostatistical Models via Local Karhunen-Loeve Expansion. \emph{Journal of the Royal Statistical Society}, 76, 817-832.
#'   
#'   Ludwig, G., Chu, T., Zhu, J., Wang, H. and Koehler, K. (2015) Static and Roving Sensor Data Fusion for Spatio-Temporal Mapping with Application to Occupational Exposure Assessment, \emph{to appear}.
#'
#' @seealso \code{\link{tfpca}}
#' @keywords Spatial Statistics
#' @keywords Functional Data Analysis
stdf <- function(training.set, subtfpca = NULL, ssensors = 6, 
                 L = 2, spline.df = NULL, zeta = 0, fpca.df = 20, 
                 homogeneous = FALSE, matern.nu = 0.5,
                 method = c("L-BFGS-B","Nelder-Mead"),
                 verbose = TRUE, ...){
  
  if(!is.null(subtfpca) && sum(subtfpca) == length(subtfpca)){
    message("Only static sensors, consider setting homogeneous = TRUE or subtfpca = NULL")
  } 
  if(matern.nu <= 0){
    message("Select an appropriate smoothness parameter for the Matern covariance function (nu > 0)")
  } 
  method <- match.arg(method)
  nu <- matern.nu
  
  nTR <- nrow(training.set)
  
  # Sets up spline, i.e. XTR <- cbind(1, x, y, S(t))
  # Potentially include covariates here?
  XTR <- cbind(1, training.set[ ,3:4], splines::bs(training.set[ ,2], df = spline.df))
  #!# Assuming response is in first column of training set
  YTR <- training.set[ ,1] 
  
  Step1 <- .lm.fit(XTR, YTR) # Fastest least squares
  res0 <- Step1$residuals
  
  if(is.null(subtfpca)){
    NoiStDe <- matrix(res0, ncol = ssensors)
    TOTAL.fit <- t.fit <- unique(training.set[ ,2])
  } else {
    NoiStDe <- matrix(res0[subtfpca], ncol = ssensors)
    t.fit <- unique(training.set[subtfpca,2])
    TOTAL.fit <- unique(training.set[ ,2])
  }
  Step2 <- tfpca(MatY = NoiStDe, L = L, t.fit = t.fit, 
                 zeta = zeta, tbas = fpca.df, ...) 
  lamb.est <- Step2$values
  Phi.est <- with(Step2, fda::eval.fd(TOTAL.fit, harmfd)*t(matrix(sqrt(nObs/etan), L, length(TOTAL.fit))))
  
  # Step3:   Spatial Parameters
  
  DTR <- matrix(0, nrow = nTR, ncol = nTR)  
  for(i in 1:nTR) {
    DTR[i, ] <- sqrt((training.set[i,3] - training.set[,3])^2 + (training.set[i,4] - training.set[,4])^2)
  }
  
  Dmax <- max(DTR)
  Vmax <- var(res0)
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
  
  if(method == "Nelder-Mead"){
    if(!homogeneous){
      prof.max <- constrOptim(theta0, logProfileCpp, grad = NULL, 
                              ui = UI, ci = CI, # Constraints
                              DTR = DTR, Y = YTR, XTR = XTR, 
                              subsetStatic = subsetStatic,
                              PhiTime = Phi.est[match(training.set[ ,2], TOTAL.fit),], 
                              LambEst = lamb.est, nu = nu)
    } else {
      prof.max <- constrOptim(theta0, logProfileCppH, grad = NULL,
                              ui = UI, ci = CI, # Constraints
                              DTR = DTR, Y = YTR, XTR = XTR, 
                              PhiTime = Phi.est[match(training.set[ ,2], TOTAL.fit),], 
                              LambEst = lamb.est, nu = nu)
    }
  } else {
    if(!homogeneous){
      prof.max <- constrOptim(theta0, logProfileCpp, grad = dlogProfileCpp, 
                              ui = UI, ci = CI, # Constraints
                              DTR = DTR, Y = YTR, XTR = XTR, 
                              subsetStatic = subsetStatic,
                              PhiTime = Phi.est[match(training.set[ ,2], TOTAL.fit),], 
                              LambEst = lamb.est, nu = nu)
    } else {
      prof.max <- constrOptim(theta0, logProfileCppH, grad = dlogProfileCppH,
                              ui = UI, ci = CI, # Constraints
                              DTR = DTR, Y = YTR, XTR = XTR, 
                              PhiTime = Phi.est[match(training.set[ ,2], TOTAL.fit),], 
                              LambEst = lamb.est, nu = nu)
    }
  } 
  
  theta <- prof.max$par
  if(is.null(subtfpca) & !homogeneous) {
    theta <- theta[-length(theta)] 
  }
  
  if(verbose) cat("Optimization done. \n")
  
  PhiTime <- Phi.est[match(training.set[ ,2], TOTAL.fit),]
  psi.cov <- evalPsi(DTR, L, lamb.est, theta, PhiTime, PhiTimeTE = NULL,
                     homogeneous, subsetStatic, kriging = FALSE)
  
  beta.est <- solve(crossprod(XTR, solve(psi.cov, XTR)), crossprod(XTR, solve(psi.cov, YTR)))
  resid <- YTR-XTR%*%beta.est
  
  beta.est <- as.numeric(beta.est)
  names(beta.est) <- c("b0", "bx", "by",
                       paste0("sp", seq_len(length(beta.est)-3))) 
  ret <- list(beta.est = beta.est, spatCov = theta[1:L])
  ret$sigma2s <- theta[L+1]
  if(!is.na(theta[L+2]))
    ret$sigma2r <- theta[L+2]
  ret$phi <- Phi.est[order(TOTAL.fit),]
  ret$times <- sort(TOTAL.fit)
  ret$tfpca.params <- Step2
  ret$training.set <- training.set
  ret$subtfpca <- subtfpca
  ret$homogeneous <- homogeneous
  ret$residuals <- resid
  ret$spline.df <- spline.df
  class(ret) <- "stdf"
  return(ret)
}