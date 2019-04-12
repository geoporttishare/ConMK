#' Quick autocorrelation of lag 1
#'
#' @param x values
#'
#' @export
autocor_lag1 <- function(x, ..., check=TRUE) {
  if(check) x <- check_x(x)
  c_autocorrelation_1(x[,1])
}
