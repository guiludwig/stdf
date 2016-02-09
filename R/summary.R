#' Summarizing STDF Model Fits
#' 
#' \code{summary} method for class \code{"stdf"}.
summary.stdf <- function(object, ...) {
  cat("\nSpatio-Temporal Data Fusion object\n\n", sep = "")
  cat("\nResiduals:\n", sep = "")
  temp <- quantile(residuals(object))
  names(temp) <- c("Min","1Q","Median","3Q","Max")
  print(temp, digits = 4)
  cat("\n")
  print(object, summary = TRUE, ...)
}