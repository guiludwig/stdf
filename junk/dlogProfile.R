dlogProfile <- function(theta, DTR, Y, XTR, PhiTime, LambEst){
  # Compute the correlation matrix.
  X <- XTR
  N <- length(Y)
  psi <- matrix(0, nrow = N, ncol = N)
  Gamma1 <- LambEst[1]*exp(-DTR/theta[1])
  Gamma2 <- LambEst[2]*exp(-DTR/theta[2])
  Gamma3 <- LambEst[3]*exp(-DTR/theta[3])
  psi <- Gamma1*outer(PhiTime[,1],PhiTime[,1]) + 
         Gamma2*outer(PhiTime[,2],PhiTime[,2]) + 
         Gamma3*outer(PhiTime[,3],PhiTime[,3])
  diag(psi) <- diag(psi)+theta[4]
  U <- chol(psi)
  psi.inv <- chol2inv(U)
  # Compute the MLE for beta.
  beta <- solve(crossprod(X, psi.inv %*% X), crossprod(X, psi.inv %*% Y))
  # Compute the MLE for sigma.
  resid <- Y - X %*% beta
  
  psi1 <- 1*DTR*Gamma1*outer(PhiTime[,1],PhiTime[,1])/(theta[1]^2) # if c = 1/theta
  psi2 <- 1*DTR*Gamma2*outer(PhiTime[,2],PhiTime[,2])/(theta[2]^2) # change to -1
  psi3 <- 1*DTR*Gamma3*outer(PhiTime[,3],PhiTime[,3])/(theta[3]^2) # change to -1
  sigmaRes <- psi.inv %*% resid
  theta1 <- sum(diag(psi.inv%*%psi1)) -1*crossprod(sigmaRes, psi1%*%sigmaRes)
  theta2 <- sum(diag(psi.inv%*%psi2)) -1*crossprod(sigmaRes, psi2%*%sigmaRes)
  theta3 <- sum(diag(psi.inv%*%psi3)) -1*crossprod(sigmaRes, psi3%*%sigmaRes)
  sigma2 <- -1*crossprod(sigmaRes) + sum(diag(psi.inv))
  # Log of the determinant of U.
  prof <- c(theta1, theta2, theta3, sigma2)
  return(prof)
  # This is a test
}
# system.time({dlogProfile(theta = c(1,1.5,2,1), DTR = DTR, Y = YTR, XTR = XTR, 
#                          PhiTime = Phi.est[match(training.set[ ,2], t.fit),], 
#                          LambEst = lamb.est)})
# OLD: 0.705
# NEW: 0.357
