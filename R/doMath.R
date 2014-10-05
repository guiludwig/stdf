doMath <- function(NoiTR, NoiTE, testset, trainset, ssensors=6, verbose=TRUE){
  
  NoiS <- matrix(NoiTR[trainset, 1], ncol=ssensors) #!# IMPORTANT
  
  nTR <- dim(NoiTR)[1]  
  nTS <- dim(NoiTE)[1]
  nT <- nTR + nTS
  
  SpT <- bs(1:dim(NoiS)[1]) #!# Can insert knots here
  XTR <- cbind(1,NoiTR[,3:4],SpT[NoiTR[,2],]) #!# Gui: XTR=cbind(1,NoiTR$X, NoiTR$Y, SpT[NoiTR[,2],])
  YTR <- NoiTR[,1] #!# Gui: YTR=NoiTR$Leq
  
  Step1 <- lm(YTR~XTR-1)
  NoiTRDe <- NoiTR
  NoiTRDe[,1] <- Step1$res
  
  nT <- dim(NoiS)[1]
  nSt <- dim(NoiS)[2]
  # teststa <- NoiTR[trainset[1:round(.1*length(trainset))],2]
  # SOMETHING WRONG HERE
  # Has Detrended data per column
  # NoiStDe <- matrix(NoiTRDe[1:(nSt*(nT-length(teststa))),1],ncol=nSt)
  NoiStDe <- matrix(NoiTRDe[trainset, 1], ncol=ssensors) #!# IMPORTANT
  
  J <- 2 # trainset[1:round(.1*length(trainset))]
  # -c(trainset[1:round(.1*length(trainset))])
  t.pred <- t.fit <- unique(NoiTR[,2]/(diff(range(NoiTR[,2]))+2))  # (1:nT)[-teststa]/(nT+1)
  # t.pred <- NoiTE[,2]/(diff(range(NoiTE[,2]))+2) # teststa/(nT+1)
  Step2 <- myfpca(NoiStDe, J, t.fit, t.pred, 20) #!# based on Ramsay's ANOVA PCA
  lamb.est <- Step2$values
  Phi.est <- matrix(0,nT,2)
  Phi.est <- Phi.TR <- Phi.TE <- Step2$vectors
  #Phi.est[teststa,] <- Step2$pred.vec
  
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
  prof.max <- constrOptim(theta0, logProfileCpp, grad = NULL,
                          ui = cbind(c(1,0,0,-1,0,0),c(0,1,0,0,-1,0), c(0,0,1,0,0,-1)), 
                          ci = c(0.001, 0.0001, 0.001, -Dmax, -Dmax, -Vmax), 
                          DTR = DTR, Y = NoiTR[,1], XTR = XTR, 
                          PhiTime = Phi.est[NoiTR[,2],], LambEst = lamb.est)
  
  theta <- prof.max$par
  
  if(verbose) cat("Optimization done. \n")
  
  #### covaraince matrix
  PhiTime <- Phi.est[NoiTR[,2],]
  psi.cov <- matrix(0, nrow = nTR, ncol = nTR)
  Gamma1 <- lamb.est[1]*exp(-DTR/theta[1])
  Gamma2 <- lamb.est[2]*exp(-DTR/theta[2])
  for (i in 1:nTR){
    psi.cov[i,]=Gamma1[i,]*PhiTime[i,1]*PhiTime[,1]+Gamma2[i,]*PhiTime[i,2]*PhiTime[,2]
  }
  diag(psi.cov) <- diag(psi.cov)+theta[3]
  
  beta.est <- solve(crossprod(XTR, solve(psi.cov, XTR)), crossprod(XTR, solve(psi.cov, YTR)))
  resid <- YTR-XTR%*%beta.est
  
  # Starts Kriging, note it does extrapolation!
  
  DKrig <- psi.krig <- matrix(0, nrow = nTR, ncol = length(testset))        
  
  for(i in 1:nTR) {
    DKrig[i,] <- sqrt((NoiTR[i,3] - NoiTE[,3])^2 + (NoiTR[i,4] - NoiTE[,4])^2)
  }
  
  PhiTimeTE <- Phi.est[NoiTE[,2],]
  Gamma1TE <- lamb.est[1]*exp(-DKrig/theta[1])
  Gamma2TE <- lamb.est[2]*exp(-DKrig/theta[2])
  for (i in 1:nTR){
    psi.krig[i,] <- Gamma1TE[i,]*PhiTime[i,1]*PhiTimeTE[,1]+Gamma2TE[i,]*PhiTime[i,2]*PhiTimeTE[,2]
  }
  XTE <- cbind(1,NoiTE[,3:4],SpT[NoiTE[,2],])
  YKrig <- XTE%*%beta.est+crossprod(psi.krig, solve(psi.cov, resid))
  MSPEKrig <- mean((NoiTE[,1]-YKrig)^2)
  
  beta.est <- as.numeric(beta.est)
  names(beta.est) <- c("b0", "bx", "by", "sp1", "sp2", "sp3") # 3 basis function, time using defaults
  ret <- list(beta.est = beta.est, MSPEKrig = MSPEKrig)
  return(ret)
}