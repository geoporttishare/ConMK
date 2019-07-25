#' Check 1D timeseries
#'
#' For internal use.
#'
#' @param x time series object
#' @param ... ignored
#' @export


check_x <- function(x ,...) {
  x <- cbind(x)
  n <- nrow(x)
  if(n < 3) stop("too few values to do anything useful.")
  y <- x[,1]
  if(ncol(x) == 1){
    time <- 1:n
  }
  else if(ncol(x)==2) {
    time <- as.numeric(x[,2])
  }
  else stop("x should be a n-vector or n x 2 matrix with cbind(values, times)")
  cbind(y, time)
}
