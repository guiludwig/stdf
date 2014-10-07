library(MASS)
library(fda)
library(splines)
library(Rcpp)
library(RcppEigen)

load("SimulatedData200.Rdata")
source("./R/myfpca.R")
source("./R/logProfile.R")
source("./R/dlogProfile.R")
source("./R/tlogProfile.R")
sourceCpp("./src/logProfile.cpp")
sourceCpp("./src/dlogProfile.cpp")

CASES <- LETTERS[1:12]
total.cases <- length(CASES)

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

##########################################################################################################

##############################################
################## FITTING ###################
##############################################

############
# BASELINE #
############

totalIter <- total.cases*3*3*B
vecMSPE <- numeric(totalIter)
matBETA <- matrix(0, nrow=totalIter, ncol=6)
Iter <- 1

for(cas in CASES) {
  for(mu.model in c("A", "B", "C")) {
    for(cov.model in c("I", "S", "TS")){
      for(b in 1:B){
        cat("Model: ", cas, mu.model, cov.model, b, "\n", sep="")
        current.case <- get(paste0("case", cas))
        temp <- current.case[c("x","y","t","id")]
        temp$mu <- current.case[[paste0("mu.", mu.model)]]
        temp$Y <- temp$mu + current.case[[paste0(cov.model, sprintf("%003d",b))]]
        
        testset <- which(temp$id == "FakePE")
        NoiTR <- matrix(0, ncol=4, nrow=sum(temp$id != "FakePE"))
        NoiTE <- matrix(0, ncol=4, nrow=sum(temp$id == "FakePE"))
        NoiTR[,1] <- temp$Y[-testset]
        NoiTR[,2] <- temp$t[-testset]
        NoiTR[,3] <- temp$x[-testset]
        NoiTR[,4] <- temp$y[-testset]
        NoiTE[,1] <- temp$Y[testset]
        NoiTE[,2] <- temp$t[testset]
        NoiTE[,3] <- temp$x[testset]
        NoiTE[,4] <- temp$y[testset]
        trainset <- which(temp$id != "FakePE") # grep("Fix*", temp$id)
        ssensors <- sum(unique(as.numeric(gsub(".*\\D(\\d{1,})\\D.*", "\\1" , temp$id[trainset]))))
        print(ssensors)
        tempRes <- doMath(NoiTR, NoiTE, testset, trainset, ssensors=ssensors)
        vecMSPE[Iter] <- tempRes$MSPEKrig
        matBETA[Iter, ] <- tempRes$beta.est
        Iter <- Iter + 1
      }
    }
  }
}

results <- vector("list", length=total.cases*3*3*B)
names(results) <- paste0(rep(CASES, times=1, each=3*3*B), rep(c("A", "B", "C"), times=total.cases, each=3*B),
                         rep(c("I", "S", "TS"), times=total.cases*3, each=B), rep(1:B, times=total.cases*3*3, each=1))
results <- lapply(results, function(z) `<-`(z, list(beta.est = numeric(6), MSPEKrig = 0)))

Iter <- 1
for(cas in CASES) {
  for(mu.model in c("A", "B", "C")) {
    for(cov.model in c("I", "S", "TS")){
      for(b in 1:B){
        results[[paste0(cas, mu.model, cov.model, b)]]$beta.est <- matBETA[Iter, ]
        results[[paste0(cas, mu.model, cov.model, b)]]$MSPEKrig <- vecMSPE[Iter]
        Iter <- Iter + 1
      }
    }
  }
}
save(results, file=paste0("results",B,".Rdata"))