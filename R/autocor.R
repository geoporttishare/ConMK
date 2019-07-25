#' Quick Autocorrelation of Lag 1
#'
#' 1D only, inteded for internal use.
#'
#' @param x values
#' @param ... ignored
#' @param check check input?
#'
#' @export
autocor_lag1 <- function(x, ..., check=TRUE) {
  if(check) x <- check_x(x)[,1]

  c_autocorrelation_1(x)
}
