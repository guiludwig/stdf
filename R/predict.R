#' Predict method STDF Model Fits
#' 
#' Predicted values for a \code{stdf} object.
predict.stdf <- function(object, newdata = NULL, 
                         what = c("fit", "predict", "loadings"), 
                         se = FALSE, 
                         loadings.range = list(x = c(0,1), 
                                               y = c(0,1), 
                                               mesh = 20), ...) {
  # NOTES:
  # * Make sure spline.df is matching when calculating kriging s.e.
  # * Rethink what's the best way to load data into model; too much 
  #   copying of the original dataset results in inefficiency.
  # * Did I use the right formula for se? Was it the most efficient?
  
  what <- match.arg(what)
  if(what == "predict" & !(is.data.frame(newdata)|is.matrix(newdata))){
    stop("Please provide a \"newdata\" object with entries matching the training.set columns.")
    # Need to write some test to match names
  }
  if(what == "fit" & is.data.frame(newdata)|is.matrix(newdata)){
    what <- "predict"
  } 
  if(what == "loadings" & se){
    stop("Standard error for loadings not implemented yet!")
  }
  
  beta.est <- object$beta.est
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
  if(!L) stop("Spatial parameters are missing!\n")
  Phi.est <- with(tfpca.params, fda::eval.fd(TOTAL.fit, harmfd)*t(matrix(sqrt(nObs/etan), L, length(TOTAL.fit))))
  PhiTime <- Phi.est[match(training.set[ ,2], TOTAL.fit),]
  psi.cov <- evalPsi(DTR, L, lamb.est, theta, PhiTime, PhiTimeTE = NULL,
                     object$homogeneous, subsetStatic, kriging = FALSE)
  if(is.null(newdata)){
    prediction.set <- training.set
  } else {
    prediction.set <- newdata
  }
  nTS <- nrow(prediction.set)
  t.pred <- unique(prediction.set[ ,2])
  
  if(what == "fit" | what == "predict"){
    # Starts Kriging
    DKrig <- matrix(0, nrow = nTR, ncol = nTS)        
    for(i in 1:nTR) {
      DKrig[i,] <- sqrt((training.set[i,3] - prediction.set[,3])^2 + (training.set[i,4] - prediction.set[,4])^2)
    }
    DTE <- matrix(0, nrow = nTS, ncol = nTS)        
    for(i in 1:nTS) {
      DTE[i,] <- sqrt((prediction.set[i,3] - prediction.set[,3])^2 + (prediction.set[i,4] - prediction.set[,4])^2)
    }

    Phi.estTE <- with(tfpca.params, fda::eval.fd(t.pred, harmfd)*t(matrix(sqrt(nObs/etan), L, length(t.pred))))
    PhiTimeTE <- Phi.estTE[match(prediction.set[ ,2], t.pred), ]
    psi.krig <- evalPsi(DKrig, L, lamb.est, theta, PhiTime, PhiTimeTE,
                        object$homogeneous, subsetStatic, kriging = TRUE)
    var.krig <- evalPsi(DTE, L, lamb.est, theta, PhiTimeTE, PhiTimeTE,
                        object$homogeneous, subsetStatic, kriging = TRUE)
    
    XTE <- cbind(1, prediction.set[,3:4], splines:::predict.bs(splines::bs(training.set[ ,2], df = object$spline.df), prediction.set[ ,2]))
    YKrig <- as.numeric(XTE%*%beta.est+crossprod(psi.krig, solve(psi.cov, object$resid)))
    MSPEKrig <- mean((prediction.set[,1]-YKrig)^2)
    
    ret <- list(fit = YKrig, MSPEKrig = MSPEKrig)
    if(se){
      XTR <- cbind(1, object$training.set[ ,3:4], 
                   splines::bs(object$training.set[ ,2], df = object$spline.df))
      # Is this right?
      temp.XTR <- solve(psi.cov, XTR) # Cz^{-1}X in C&W, p148
      # c'Cz^{-1}c + (x_0 - X'Cz^{-1}c)'(X'Cz^{-1}X)^{-1}(x_0 - X'Cz^{-1}c)
      # cat("temp.XTR = ", dim(temp.XTR), "\n")
      # cat("psi.krig = ", dim(psi.krig), "\n")
      # cat("psi.cov  = ", dim(psi.cov), "\n")
      # cat("     XTR = ", dim(XTR), "\n")
      # cat("     XTE = ", dim(XTE), "\n")
      # cat("   X'SX  = ", dim(crossprod(temp.XTR, XTR)), "\n")
      # cat("   c  = ", dim(t(XTE - crossprod(psi.krig, temp.XTR))), "\n")
      SE <- var.krig - crossprod(psi.krig, solve(psi.cov, psi.krig)) +
        (XTE - crossprod(psi.krig, temp.XTR))%*%solve(crossprod(temp.XTR, XTR),
                                                      t(XTE - crossprod(psi.krig, temp.XTR)))
      ret$se.fit <- sqrt(diag(SE))
      # SE <- numeric(nrow(XTE))
      # TEMP <- solve(crossprod(temp.XTR, XTR))
      # for(i in 1:nrow(XTE)){
      #   SE[i] <- as.numeric(crossprod(psi.krig[ , i, drop=FALSE], solve(psi.cov, psi.krig[ , i, drop=FALSE])) +
      #                                   t(XTE[i,] - crossprod(XTR, solve(psi.cov,psi.krig[ , i, drop=FALSE])))%*%TEMP%*%(XTE[i,] - crossprod(XTR, solve(psi.cov,psi.krig[ , i, drop=FALSE]))))
      # }
      # ret$se.fit <- sqrt(SE)
    }
  } else if(what == "loadings"){
    # Evaluate loadings
    new.s <- with(loadings.range, expand.grid(x = seq(x[1], x[2], length.out = mesh), 
                                              y = seq(y[1], y[2], length.out = mesh)))
    nTS <- nrow(new.s)
    ret <- list()
    DKrig <- matrix(0, nrow = nTR, ncol = nTS)        
    for(i in 1:nTR) {
      DKrig[i,] <- sqrt((training.set[i,3] - new.s$x)^2 + 
                          (training.set[i,4] - new.s$y)^2)
    }  
    for(j in 1:L){
      tempPhiTime <- matrix(rep(PhiTime[,j], ncol(DKrig)), byrow = FALSE, 
                            ncol = ncol(DKrig), 
                            nrow = nrow(DKrig))
      tempPsi <- lamb.est[j]*exp(-DKrig/theta[j])*tempPhiTime
      tempLoad <- crossprod(tempPsi, solve(psi.cov, object$resid))
      ret[[j]] <- matrix(as.numeric(tempLoad), nrow = loadings.range$mesh)
    }
    ret$loadings.range <- new.s
  }

  ret$what <- what
  return(ret)
}

# a <- predict(results, what="loadings", loadings.range = list(x = c(1,5), y = c(1,5), mesh = 20))
# image(a[[1]])