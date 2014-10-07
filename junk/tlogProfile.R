log.profile=function(theta, DTR, NoiTR, XTR, Phi.est, lamb.est){
  # Compute the correlation matrix.
  Y = NoiTR[,1]
  X = XTR
  N = length(Y)
  PhiTime=Phi.est[NoiTR[,2],]
  psi=matrix(0, nrow = N, ncol = N)
  Gamma1=lamb.est[1]*exp(-DTR/theta[1])
  Gamma2=lamb.est[2]*exp(-DTR/theta[2])
  for (i in 1:N){
    psi[i,]=Gamma1[i,]*PhiTime[i,1]*PhiTime[,1]+Gamma2[i,]*PhiTime[i,2]*PhiTime[,2]
  }
  diag(psi)=diag(psi)+theta[3]
  U <- chol(psi)
  U.inv <- backsolve(U, diag(1, nrow = N))
  psi.inv <- U.inv %*% t(U.inv)
  # Compute the MLE for beta.
  beta <- solve(t(X) %*% psi.inv %*% X) %*% t(X) %*% psi.inv %*% Y
  # Compute the MLE for sigma.
  resid <- Y - X %*% beta
  sigma2 <- t(resid) %*% psi.inv %*% resid/N
  # Evaluate -log-profile likelihood.
  log.U.det <- sum(log(diag(U)))
  # Log of the determinant of U.
  prof <- log.U.det + (N * log(sigma2))/2 + N/2
  return(prof)
}
