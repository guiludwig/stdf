predict.stdf <- function(model, newdata = NULL, 
                         what = c("fitted", "predict", "loadings")) {
  what <- match.call(what)
  attach(model)
  
  if(!is.null(newdata)){
    nTS <- nrow(newdata)
    t.pred <- unique(newdata[ ,2])
  } else {
    # repeat things?
  }
  
  if(what == "predict"){
    # Starts Kriging
    nTS <- nrow(newdata)
    DKrig <- matrix(0, nrow = nTR, ncol = nTS)        
    for(i in 1:nTR) {
      DKrig[i,] <- sqrt((training.set[i,3] - newdata[,3])^2 + (training.set[i,4] - newdata[,4])^2)
    }
    
    PhiTimeTE <- with(tfpca.params, fda::eval.fd(t.pred, harmfd)*t(matrix(sqrt(nObs/etan), L, length(t.pred))))
    psi.krig <- evalPsi(DKrig, L, lamb.est, theta, PhiTime, PhiTimeTE,
                        homogeneous, subsetStatic, kriging = TRUE)
    
    XTE <- cbind(1, prediction.set[,3:4], splines::bs(prediction.set[ ,2], df = spline.df))
    YKrig <- XTE%*%beta.est+crossprod(psi.krig, solve(psi.cov, resid))
    MSPEKrig <- mean((prediction.set[,1]-YKrig)^2)
    
    ret <- list(fitted = YKrig, MSPEKrig = MSPEKrig)
  } else if(what == "loadings"){
    # Evaluate loadings
    
    static.loadings <- crossprod(Phi.est[order(t.fit),], 
                                 matrix(solve(psi.cov, resid)[as.logical(subsetStatic)],
                                        ncol = ssensors)[order(t.fit),])
    for(ell in 1:L)
      static.loadings[ell,] <- lamb.est[ell]*static.loadings[ell,]
    colnames(static.loadings) <- paste0("curve ",1:ncol(static.loadings))
    static.loadings <- t(static.loadings)
    ret$static.loadings <- static.loadings
  }

  class(ret) <- "stdfPred"
  dettach(model)
  return(ret)
}