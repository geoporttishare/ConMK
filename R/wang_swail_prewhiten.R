#' Prewhitening Raster Stack
#'
#' rasterStack version
#'
#' @param x rasterStack
#' @param ... Passed on to \link{wang_swail_prewhiten_1d}
#'
#' @details Apply Wang&Swail prewhitening pixelwise. The method assumes AR(1) auto-correlation.
#' The resulting stack will have first layer just NA.
#' @seealso \link{wang_swail_prewhiten_1d}
#' @references
#' Wang, X.L. and Swail, V.R. (2000) Changes of Extreme Wave Heights
#' in Northern Hemisphere Oceans and Related Atmoshperic Circulation Regimes, \emph{Journal of Climate}
#' @import raster
#' @export

wang_swail_prewhiten_stack <- function(x, ...) {
  fun <- function(v) {
    if(any(is.na(v))) return(rep(NA, length(v)))
    c(NA, wang_swail_prewhiten_1d(v, ...)$W)
  }
  rW <- calc(x, fun)
  rW
}

#' Prewhitening 1D
#'
#' Simple 1d timeseries: attempt to remove autocorrelation. Assumes AR(1) errors, and 1 step time-grid.
#'
#' @param x either vector or matrix with c(x, time)
#'
#' @details time is ignored, and assumed to be 1-step integers 1:nrow(x).
#'
#' @return list with element W containing the whitened series, one shorter than the original.
#'
#' @references
#' Wang, X.L. and Swail, V.R. (2000) Changes of Extreme Wave Heights
#' in Northern Hemisphere Oceans and Related Atmoshperic Circulation Regimes, \emph{Journal of Climate}
#' @export

wang_swail_prewhiten_1d <- function(x, ..., eps = 1e-4, rho_th = 0.05, itmax=20, useC=FALSE) {
  #
  x <- check_x(x)
  y <- x[,1]
  n <- length(y)
  time <- x[,2]

  if(useC) return( c_wang_swail_prewithen_1d(y, time, eps, rho_th, itmax)  )
  #
  rho0 <- 2
  rho <- 0
  beta0 <- beta <- 0
  it <- 0
  # iterate
  while( (abs(rho - rho0) > eps | abs(beta0 - beta) > eps) & it < itmax){
    beta0 <- beta
    rho0 <- c_autocorrelation_1(y - time * beta0)[3]
    W <- (y[-1] - rho0 * y[-n])/(1-rho0)
    beta0 <- median(c_mann_kendall_test_and_beta(W, time)$X)
    rho <- c_autocorrelation_1(y - time * beta0)[3]
    W <- (y[-1] - rho * y[-n])/(1-rho)
    beta <- median(c_mann_kendall_test_and_beta(W, time)$X)
    it <- it + 1
  }
  #if(it == itmax) warning("itmax reached.")
  list(W=W, rho = rho, beta = beta, iter = it)
}

