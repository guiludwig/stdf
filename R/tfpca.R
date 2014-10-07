#' Temporal Functional Principal Component
#'
#' This function is used on step II of the GFDA algorithm by Chu, Zhu and Wang (2014). It does 
#' not need to be used directly, since it's mostly an wrapper to function calls from J. O. Ramsay 
#' et al.'s \code{fda} package function calls.
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
#' mean(rnorm(20)) # Write some test cases!
#' 
#' @references
#'  Ramsay, J. O. and Silverman, B. (2005) Functional Data Analysis. New York: Springer.
#'
#' @seealso \code{\link{median}},
#'   \code{\link{mean}}
#' @keywords Spatial Statistics
#' @keywords Functional Data Analysis
tfpca <- function(MatY, J, t.fit, t.pred, tbas){
  nt <- dim(MatY)[1]
  timebasis <- fda::create.bspline.basis(c(0,1), nbasis=tbas)
  timefdPar <- fda::fdPar(timebasis)
  daytempfd <- fda::smooth.basis(t.fit, MatY, timebasis)$fd
  pcafd <- fda::pca.fd(daytempfd, nharm=J, timefdPar)
  harmfd <- pcafd[[1]]
  fdmat.all <- fda::eval.fd(c(t.fit,t.pred), harmfd)
  fdmat <- fdmat.all[1:nt,]
  
  etan <- apply(fdmat^2,2,sum)
  values <- pcafd$values[1:J]*etan/nt
  vectors <- fdmat*t(matrix(sqrt(nt/etan),J,nt))
  fd.pre <- fdmat.all[(nt+1):(nt+length(t.pred)),]*t(matrix(sqrt(nt/etan),J,length(t.pred)))
  return(list(values=values, vectors=vectors, pred.vec=fd.pre))
}
