#' Temporal Functional Principal Component
#'
#' This function is used on step II of the GFDA algorithm by Chu, Zhu and Wang (2014). It does 
#' not need to be used directly, since it's mostly an wrapper to function calls from J. O. Ramsay 
#' et al.'s \code{fda} package function calls. It behaves similarly to the \code{pca.fd} function, 
#' but rescales the eigenfunctions/eigenvalues to a scale that's needed by the \code{gfda} function.
#' Notice \code{tfpca} do not handle NA's in the data values. It might be useful do to some
#' kind of nearest neightbor imputation if possible. The \code{fda} function \code{data2fd} does
#' a least-squares fit of the basis function, which can be helpful in preprocessing the data.
#'
#' @param MatY The data organized in a matrix, with each time series in a column. Missing 
#' values are currently not handled properly.
#' @param L The number of functional principal components to be returned.
#' @param t.fit Values of the time-domain at which a value is fitted.
#' @param t.pred Values of the time-domain at which a predicted value needs to be produced.
#' @param tbas Number of spline basis functions used in the smooth.
#'
#' @export
#' @return List of three items 
#'   \item{values}{numeric vector of L Eigenvalue estimates}
#'   \item{vectors}{matrix of L columns with Eigenfunction estimates observed at t.fit}
#'   \item{pred.vec}{matrix of L columns with Eigenfunction estimates observed at t.pred}
#'
#' @examples
#' # This example uses fda's CanadianWeather dataset:
#' library(fda)
#' test <- tfpca(CanadianWeather$dailyAv[,,"Temperature.C"], 2, 1:365, 32:60, 20)
#' layout(matrix(1:2, ncol=2))
#' plot(test$vectors[,1], type="l", main="First Principal Component")
#' plot(test$vectors[,2], type="l", main="Second Principal Component")
#' 
#' @references
#'  Ramsay, J. O. (2006) Functional Data Analysis. New York: Springer.
#'
#' @seealso \code{\link{pca.fd}}, \code{\link{gfda}}
#' @keywords Spatial Statistics
#' @keywords Functional Data Analysis
tfpca <- function(MatY, L, t.fit, t.pred, tbas=10){
  nObs <- nrow(MatY)
  rng <- range(c(t.fit, t.pred))
  timebasis <- fda::create.bspline.basis(rng, nbasis=tbas) # Smoothness controlled with # basis
  timefdPar <- fda::fdPar(timebasis) # Sets parameters getting defaults (lambda = 0 => no additional smoothing)
  smtempfd <- fda::smooth.basis(t.fit, MatY, timebasis)$fd # Preliminary smoothing
  pcafd <- fda::pca.fd(smtempfd, nharm = L, timefdPar)
  harmfd <- pcafd[[1]]
  fdmat <- fda::eval.fd(t.fit, harmfd)
  etan <- apply(fdmat^2,2,sum)
  values <- pcafd$values[1:L]*etan/nObs # Truncates at L, multiply by |phi|, divide by n
  vectors <- fdmat*t(matrix(sqrt(nObs/etan), L, nObs)) # Vectors at scale of data?
  fd.pre <- fda::eval.fd(t.pred, harmfd)*t(matrix(sqrt(nObs/etan),L,length(t.pred)))
  return(list(values=values, vectors=vectors, pred.vec=fd.pre))
}