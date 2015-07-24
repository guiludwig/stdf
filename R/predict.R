predict.stdf <- function(model, newdata = NULL, 
                         what = c("fit", "predict", "loadings"), 
                         loadings.range = list(x = c(0,1), y = c(0,1), mesh = 20)) {
  what <- match.arg(what)
  if(what == "predict" & !is.data.frame(newdata)){
    stop("Please provide a \"newdata\" object with entries matching the training.set columns.")
  }
  if(what == "fit" & is.data.frame(newdata)){
    what <- "predict"
  }
  
  beta.est <- model$beta.est
  theta <- c(model$spatCov, model$sigma2s)
  if(!is.null(model$sigma2r)) {
    theta <- c(theta, model$sigma2r)
  }
  tfpca.params <- model$tfpca.params
  training.set <- model$training.set
  nTR <- nrow(training.set)
  DTR <- matrix(0, nrow = nTR, ncol = nTR)  
  for(i in 1:nTR) {
    DTR[i, ] <- sqrt((training.set[i,3] - training.set[,3])^2 + (training.set[i,4] - training.set[,4])^2)
  }
  if(is.null(model$subtfpca)) {
    subsetStatic <- rep(1, nTR)
    t.fit <- unique(training.set[ ,2])
  } else {
    subsetStatic <- as.numeric(model$subtfpca)
    t.fit <- unique(training.set[model$subtfpca,2])
  }
  lamb.est <- tfpca.params$values
  Phi.est <- tfpca.params$vectors
  PhiTime <- Phi.est[match(training.set[ ,2], t.fit),]
  L <- length(model$spatCov)
  psi.cov <- evalPsi(DTR, L, lamb.est, theta, PhiTime, PhiTimeTE = NULL,
                     model$homogeneous, subsetStatic, kriging = FALSE)
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

    Phi.estTE <- with(tfpca.params, fda::eval.fd(t.pred, harmfd)*t(matrix(sqrt(nObs/etan), L, length(t.pred))))
    PhiTimeTE <- Phi.estTE[match(prediction.set[ ,2], t.pred), ]
    psi.krig <- evalPsi(DKrig, L, lamb.est, theta, PhiTime, PhiTimeTE,
                        model$homogeneous, subsetStatic, kriging = TRUE)

    XTE <- cbind(1, prediction.set[,3:4], splines::bs(prediction.set[ ,2], df = model$spline.df))
    YKrig <- as.numeric(XTE%*%beta.est+crossprod(psi.krig, solve(psi.cov, model$resid)))
    MSPEKrig <- mean((prediction.set[,1]-YKrig)^2)
    
    ret <- list(fitted = YKrig, MSPEKrig = MSPEKrig)
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
      tempPsi <- lamb.est[j]*exp(-DKrig/theta[j])
      tempLoad <- -1*crossprod(tempPsi, solve(psi.cov, model$resid))
      ret[[j]] <- matrix(as.numeric(tempLoad), nrow = loadings.range$mesh)
    }
    ret$loadings.range <- new.s
  }
  
  class(ret) <- "stdfPred"
  return(ret)
}

# a <- predict(results, what="loadings", loadings.range = list(x = c(1,5), y = c(1,5), mesh = 20))
# image(a[[1]])