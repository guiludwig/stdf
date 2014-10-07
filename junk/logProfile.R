logProfile <- function(theta, DTR, NoiTR, XTR, Phi.est, lamb.est){
  # Compute the correlation matrix.
  Y <- NoiTR[,1]
  # X <- XTR
  N <- length(Y)
  PhiTime <- Phi.est[NoiTR[,2],] # Assumes integer NoiTR
  psi <- matrix(0, nrow = N, ncol = N)
  Gamma1 <- lamb.est[1]*exp(-DTR/theta[1])
  Gamma2 <- lamb.est[2]*exp(-DTR/theta[2])
  psi <- Gamma1*outer(PhiTime[,1],PhiTime[,1])+Gamma2*outer(PhiTime[,2],PhiTime[,2])
  diag(psi) <- diag(psi)+theta[3]
  U <- chol(psi)
  psi.inv <- chol2inv(U)
  # Compute the MLE for beta.
  beta <- solve(crossprod(XTR, psi.inv %*% XTR), crossprod(XTR, psi.inv %*% Y))
  # Compute the MLE for sigma.
  resid <- Y - XTR %*% beta
  Nsigma2 <- crossprod(resid, psi.inv %*% resid) #/N
  # Evaluate -log-profile likelihood.
  log.U.det <- log(sum(diag(U))) # log|Sigma|/2 = log(Det(UU'))/2 = Det(U)
  # Log of the determinant of U.
  prof <- log.U.det + Nsigma2/2
  return(as.numeric(prof))
}
# system.time({logProfile(theta0, DTR = DTR, NoiTR = NoiTR, XTR = XTR, Phi.est = Phi.est, lamb.est = lamb.est)})
# OLD: 0.705
# NEW: 0.357
