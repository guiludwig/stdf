dlogProfile <- function(theta, DTR, NoiTR, XTR, Phi.est, lamb.est){
  # Compute the correlation matrix.
  Y <- NoiTR[,1]
  X <- XTR
  N <- length(Y)
  PhiTime <- Phi.est[NoiTR[,2],]
  psi <- matrix(0, nrow = N, ncol = N)
  Gamma1 <- lamb.est[1]*exp(-DTR/theta[1])
  Gamma2 <- lamb.est[2]*exp(-DTR/theta[2])
  psi <- Gamma1*outer(PhiTime[,1],PhiTime[,1])+Gamma2*outer(PhiTime[,2],PhiTime[,2])
  diag(psi) <- diag(psi)+theta[3]
  U <- chol(psi)
  psi.inv <- chol2inv(U)
  # Compute the MLE for beta.
  beta <- solve(crossprod(X, psi.inv %*% X), crossprod(X, psi.inv %*% Y))
  # Compute the MLE for sigma.
  resid <- Y - X %*% beta
  
  psi1 <- 1*DTR*Gamma1*outer(PhiTime[,1],PhiTime[,1])/(theta[1]^2) # if c = 1/theta
  psi2 <- 1*DTR*Gamma2*outer(PhiTime[,2],PhiTime[,2])/(theta[2]^2) # change to -1
  sigmaRes <- psi.inv %*% resid
  theta1 <- sum(diag(psi.inv%*%psi1)) # -1*crossprod(sigmaRes, psi1%*%sigmaRes) + sum(diag(psi.inv%*%psi1))
  theta2 <- sum(diag(psi.inv%*%psi2)) # -1*crossprod(sigmaRes, psi2%*%sigmaRes) + sum(diag(psi.inv%*%psi2))
  sigma2 <- -1*crossprod(sigmaRes) + sum(diag(psi.inv))
  # Log of the determinant of U.
  prof <- c(theta1, theta2, sigma2)
  return(prof)
}
# system.time({dlogProfile(theta0, DTR = DTR, NoiTR = NoiTR, XTR = XTR, Phi.est = Phi.est, lamb.est = lamb.est)})
# OLD: 0.705
# NEW: 0.357