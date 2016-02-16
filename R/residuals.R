residuals.stdf <- function(object, scale = FALSE, ...) {
  if(scale){
    theta <- c(object$spatCov, object$sigma2s)
    if(!is.null(object$sigma2r)) {
      theta <- c(theta, object$sigma2r)
    }    
    tfpca.params <- object$tfpca.params
    training.set <- object$training.set
    nTR <- nrow(training.set)
    DTR <- matrix(0, nrow = nTR, ncol = nTR)  
    for(i in 1:nTR) {
      DTR[i, ] <- sqrt((training.set[i,3] - training.set[,3])^2 + (training.set[i,4] - training.set[,4])^2)
    }
    if(is.null(object$subtfpca)) {
      subsetStatic <- rep(1, nTR)
    } else {
      subsetStatic <- as.numeric(object$subtfpca)
    }
    TOTAL.fit <- unique(training.set[ ,2])
    lamb.est <- tfpca.params$values  
    L <- length(object$spatCov)
    Phi.est <- with(tfpca.params, fda::eval.fd(TOTAL.fit, harmfd)*t(matrix(sqrt(nObs/etan), L, length(TOTAL.fit))))
    PhiTime <- Phi.est[match(training.set[ ,2], TOTAL.fit),]    
    psi.cov <- evalPsi(DTR, L, lamb.est, theta, PhiTime, PhiTimeTE = NULL,
                       object$homogeneous, subsetStatic, kriging = FALSE)   
    RawRes <- as.numeric(object$residuals)
    SVD <- svd(psi.cov)
    # Overkill to make sure p.d.? (d + |d|)/2 sets nonpositive to zero, 
    # add then 10^{-5} to make pd instead of psd
    return(as.numeric(SVD$u%*%diag(1/sqrt((SVD$d + abs(SVD$d))/2 + 0.00001))%*%crossprod(SVD$v, RawRes)))
  } else {
    return(as.numeric(object$residuals))
  }
}