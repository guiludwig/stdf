summary.stdf <- function(object, ...) {
  cat("\nResiduals:\n", sep = "")
  temp <- quantile(residuals(object))
  names(temp) <- c("Min","1Q","Median","3Q","Max")
  print(temp, digits = 4)
  print(object, ...)
}