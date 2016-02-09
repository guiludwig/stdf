print.stdf <- function(x, digits = max(3L, getOption("digits") - 3L), summary = FALSE, ...) {
  if(!summary) {cat("\nSpatio-Temporal Data Fusion object\n\n", sep = "")}
  if (length(x$beta.est)) {
    cat("Deterministic coefficients:\n")
    print.default(format(x$beta.est, digits = digits), print.gap = 2L, 
                  quote = FALSE)
  }
  else cat("No deterministic coefficients\n")
  cat("\n")
  if (length(x$spatCov)) {
    cat("Spatial coefficients:\n")
    temp <- x$spatCov
    names(temp) <- paste0("phi",1:length(temp))
    print.default(format(temp, digits = digits), print.gap = 2L, 
                  quote = FALSE)
    temp <- x$tfpca.params$values
    names(temp) <- paste0("lambda",1:length(temp))
    print.default(format(temp, digits = digits), print.gap = 2L, 
                  quote = FALSE)
  }
  else cat("No deterministic coefficients\n")
  cat("\n")
  if (length(x$sigma2s)) {
    cat("Static instrument variance estimate: ")
    cat(format(x$sigma2s, digits = digits))
  }
  cat("\n")
  if (length(x$sigma2r)) {
    cat("Roving instrument variance estimate: ")
    cat(format(x$sigma2r, digits = digits))
  }
  cat("\n")
  invisible(x)
}