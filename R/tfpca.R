#' Temporal Functional Principal Component
#'
#' This function is used on step II of the STDF algorithm. It does not need to be used 
#' directly, since it's mostly an wrapper to \code{fda} package function calls. It 
#' behaves similarly to the \code{pca.fd} function, but rescales the 
#' eigenfunctions/eigenvalues in a way that's needed by the \code{stdf} function.
#' Notice \code{tfpca} do not handle NA's in the data values. It might be useful do to some
#' kind of nearest neightbor imputation if possible. The \code{fda} function \code{data2fd} does
#' a least-squares fit of the basis function, which can be helpful in preprocessing the data.
#' The tuning parameter lambda must be passed manually.
#'
#' @param MatY The data organized in a matrix, with each time series in a column. Missing 
#' values are currently not handled properly.
#' @param L The number of functional principal components to be returned.
#' @param t.fit Values of the time-domain at which a value is fitted.
#' @param lambda Smoothness parameter, must be provided. Defaults to 0 (no smoothing).
#' @param tbas Number of spline basis functions used in the smooth.
#'
#' @export
#' @return List of three items 
#'   \item{values}{numeric vector of L Eigenvalue estimates}
#'   \item{vectors}{matrix of L columns with Eigenfunction estimates observed at t.fit}
#'   \item{harmfd}{Internal use. Prediction of the Eigenfunction at new values of t, 
#'                 say under vector t.pred, can be obtained with the next three 
#'                 parameters. The syntax is 
#'                 fda::eval.fd(t.pred, harmfd)*t(matrix(sqrt(nObs/etan), L, length(t.pred)))}
#'   \item{nObs}{Internal use. }
#'   \item{etan}{Internal use. }
#'
#' @examples
#' # This example uses fda's CanadianWeather dataset:
#' library(fda)
#' test <- tfpca(CanadianWeather$dailyAv[,,"Temperature.C"], 2, 1:365, 32:60, 0, 20)
#' test2 <- tfpca(CanadianWeather$dailyAv[,,"Temperature.C"], 2, 1:365, 32:60, 1e5, 20)
#' layout(matrix(1:4, byrow = TRUE, ncol = 2))
#' plot(test$vectors[,1], type="l", main="First Principal Component")
#' plot(test$vectors[,2], type="l", main="Second Principal Component")
#' plot(test2$vectors[,1], type="l", main="Smooth First Principal Component")
#' plot(test2$vectors[,2], type="l", main="Smooth Second Principal Component")
#' 
#' @references
#'  Ramsay, J. O. (2006) \emph{Functional Data Analysis}. New York: Springer.
#'
#' @seealso \code{\link{pca.fd}}, \code{\link{stdf}}
#' @keywords Spatial Statistics
#' @keywords Functional Data Analysis
tfpca <- function(MatY, L, t.fit, lambda = 0, tbas = 20){
  nObs <- nrow(MatY)
  if(nObs < tbas) stop("Number of basis > number of observations, please adjust \"tbas\" parameter\n")
  rng <- range(t.fit)
  timebasis <- fda::create.bspline.basis(rng, nbasis = tbas) # Smoothness controlled with # basis
  timefdPar <- fda::fdPar(timebasis, lambda = lambda) # lambda = 0 => no additional smoothing
  smtempfd <- fda::smooth.basis(t.fit, MatY, timebasis)$fd # Preliminary smoothing
  pcafd <- fda::pca.fd(smtempfd, nharm = L, timefdPar)
  harmfd <- pcafd[[1]]
  fdmat <- fda::eval.fd(t.fit, harmfd)
  etan <- apply(fdmat^2,2,sum)
  values <- pcafd$values[1:L]*etan/nObs # Truncates at L, multiply by |phi|, divide by n
  vectors <- fdmat*t(matrix(sqrt(nObs/etan), L, nObs)) # Vectors at scale of data?
  return(list(values = values,
              vectors = vectors, 
              harmfd = harmfd,
              nObs = nObs,
              etan = etan))
}