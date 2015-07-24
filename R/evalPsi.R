#' Evaluate covariance matrix for kriging
#'
#' This function is for internal use on the STDF algorithm. It evaluates the 
#' covariance matrix for fitting, and the kriging matrix for prediction.
#'
#' @param DTR Matrix of distances
#' @param L Number of functional principal components.
#' @param lamb.est Vector of eigenvalues.
#' @param theta Vector of spatial dependence parameters.
#' @param homogeneous Logical flag for homogeneous model.
#' @param subsetStatic Logical vector for subset of static sensors
#' @param kriging Logical flag, whether the function is evaluating a covariance
#'                matrix (FALSE) or a kriging cross-covariance matrix (TRUE). 
#'                Defaults to FALSE.
#'
#' @export
#' @return List of one items
#'   \item{psi.cov}{Desired covariance matrix}
#'
#' @examples
#' evalPsi(DTR = matrix(c(1,0,0,1), 2),
#'         L = 1, lamb.est = 1,
#'         theta = 1:2, PhiTime = matrix(1:2, ncol=1), PhiTimeTE = NULL,
#'         homogeneous = TRUE, subsetStatic = rep(1, 2))
#'
#' @seealso \code{\link{stdf}}
#' @keywords Spatial Statistics
evalPsi <- function(DTR, L, lamb.est, theta, PhiTime, PhiTimeTE = NULL,
                    homogeneous, subsetStatic, kriging = FALSE){
  nTR <- dim(DTR)[1]
  nTS <- dim(DTR)[2]
  psi.cov <- matrix(0, nrow = nTR, ncol = nTS)
  if(is.null(PhiTimeTE)){
    PhiTimeTE <- PhiTime
  }
  
  for(j in 1:L){
    psi.cov <- psi.cov + lamb.est[j]*exp(-DTR/theta[j])*outer(PhiTime[,j],PhiTimeTE[,j])
  }
  if(!kriging){
    diag(psi.cov) <- diag(psi.cov) + theta[L+1]
    if (!homogeneous & sum(subsetStatic) < length(subsetStatic)){
      diag(psi.cov) <- diag(psi.cov) + (theta[L+2] - theta[L+1])*(1-subsetStatic)
    } 
  }
  return(psi.cov)
}