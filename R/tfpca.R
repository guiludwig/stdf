#' Temporal Functional Principal Component
#'
#' This function is used on step II of the GFDA algorithm by Chu, Zhu and Wang (2014). It does 
#' not need to be used directly, since it's mostly an wrapper to function calls from J. O. Ramsay 
#' et al.'s \code{fda} package function calls. It behaves similarly to the \code{pca.fd} function, 
#' but rescales the eigenfunctions/eigenvalues to a scale that's needed by the \code{gfda} function.
#'
#' @param MatY The data organized in a matrix, with each time series in a column. Missing 
#' values are currently not handled properly.
#' @param J The number of functional principal components to be returned.
#' @param t.fit Values of the time-domain at which a value is fitted.
#' @param t.pred Values of the time-domain at which a predicted value needs to be produced.
#' @param tbas Number of spline basis functions used in the smooth.
#'
#' @export
#' @return List of three items 
#'   \item{values}{numeric vector of J Eigenvalue estimates}
#'   \item{vectors}{matrix of J columns with Eigenfunction estimates observed at t.fit}
#'   \item{pred.vec}{matrix of J columns with Eigenfunction estimates observed at t.pred}
#'
#' @examples
#' # This example uses fda's CanadianWeather dataset:
#' library(fda)
#' test <- tfpca(CanadianWeather$dailyAv[,,"Temperature.C"], 2, 1:365, 32:60, 20)
#' layout(matrix(1:2, ncol=2))
#' plot(test[,1], type="l", main="First Principal Component")
#' plot(test[,2], type="l", main="Second Principal Component")
#' 
#' @references
#'  Ramsay, J. O. (2006) Functional Data Analysis. New York: Springer.
#'
#' @seealso \code{\link{pca.fd}}, \code{\link{gfda}}
#' @keywords Spatial Statistics
#' @keywords Functional Data Analysis
tfpca <- function(MatY, J, t.fit, t.pred, tbas=10){
  nObs <- dim(MatY)[1]
  rng <- range(c(t.fit, t.pred))
  timebasis <- fda::create.bspline.basis(rng, nbasis=tbas) # Smoothness controlled with # basis
  timefdPar <- fda::fdPar(timebasis)
  daytempfd <- fda::smooth.basis(t.fit, MatY, timebasis)$fd
  pcafd <- fda::pca.fd(daytempfd, nharm = J, timefdPar)
  harmfd <- pcafd[[1]]
  fdmat <- fda::eval.fd(t.fit, harmfd)
  
  etan <- apply(fdmat^2,2,sum)
  values <- pcafd$values[1:J]*etan/nObs # Truncates at J, multiply by |phi|, divide by n
  vectors <- fdmat*t(matrix(sqrt(nObs/etan), J, nObs)) # Vectors at scale of data?
  fd.pre <- fda::eval.fd(t.pred, harmfd)*t(matrix(sqrt(nObs/etan),J,length(t.pred)))
  return(list(values=values, vectors=vectors, pred.vec=fd.pre))
}
