print.stdfPred <- function(x, digits = max(3L, getOption("digits") - 3L), ...) {
  cat("\nSpatio-Temporal Data Fusion object\n\n", sep = "")
  if (length(x$fitted)) {
    cat("Predicted values:\n")
    print.default(format(x$fitted, digits = digits), print.gap = 2L, 
                  quote = FALSE)
  }
  else cat("\tNo values predicted\n")
  cat("\n")
  if (length(x$MSPEKrig)) {
    cat("Mean Squared Prediction Error:\n")
        print.default(format(x$MSPEKrig, digits = digits), print.gap = 2L, 
                      quote = FALSE)
  }
  else cat("\tI could not calculate MSPE\n")
  cat("\n")
  invisible(x)
}
