#' Contextual Mann-Kendall Trend Test
#'
#' Mann-Kendall trend test for 2D rasterStacks, assuming layers
#' represent timepoints.
#'
#' @param x Raster stack, each layer denoting one timestep. Equidistant steps of time 1 assumed.
#' @param neighbourhood 2:queen, 1:rook, 0:none (classical M-K test) Only 2 and 0 implemented.
#' @param calc_slope Calculate Theil-Sen slope estimate for each series? Default:FALSE
#' @param time Optional numerical vector for the sampling times.
#'
#' @details  Assume that the values of each rasterStack location (cell) over time (layers) are a time-series. This function calculates the Mann-Kendall trend test for each series, and returns the relevant statistics as a new rasterStack.
#'
#'  The **Contextual** version (default; see ref.) simply averages the
#'  Mann-Kendall-test's S-statistic at
#'  each location over neighbours, very much like running 'focal(S, w=matrix(1/9,3,3))'
#'  on the non-contextual/cellwise statistic raster.
#'  However, this function also calculates the adjusted variances (and hence the p-values).
#'
#' NOTE assumes equal interval timeseries if 'time' is missing. Matters only for the slope calculation.
#'
#'
#' @return A rasterStack with layers: 'S' for the Mann-Kendall Statistic (classical or smoothed);
#'  's2' for the variance of the statistic; 'p' for p-value of trend detected;
#' (optional) 'slope' for the Theil-Sen slope estimate.
#'
#' @references
#' Neeti, N. and Eastman J.R. (2011) A Contextual Mann-Kendall Approach for the Assesment of Trend Significance in Image Time Series, \emph{Transactions in GIS}
#' @useDynLib ConMK
#' @import raster
#' @export

contextual_mann_kendall <- function(x, ..., neighbourhood = 2, calc_slope = FALSE,  time) {
  if(!neighbourhood %in% c(0,2)) stop("Only '0' and '2' neighbourhoods implemented.")
  # center first
  if( !canProcessInMemory(x) ) stop("'canProcessInMemory' return FALSE")
  t0 <- Sys.time()
  x <- x - mean(x)
  X <- values(stack(x))
  nr <- nrow(x)
  # check time
  if(missing(time)) time <- 1:ncol(X)
  time <- sort(time)
  if(length(time) != ncol(X)|!is.numeric(time)) stop(paste0("'time' should be a numeric increasing vector of length ", ncol(X)))
  # c_contextual_mann_kendall(v, nr)
  #browser()
  res <- c_contextual_mann_kendall( X , nr, time = time,
                                    neigh = neighbourhood, calc_slope = calc_slope)
  V <- res$S_and_s2
  Sm <- V[,1]
  s2 <- V[,2]
  D <- (Sm > 0) * 1
  Z <- (Sm+(1-2*D))/sqrt(s2)
  p <- 2*(1 - abs(pnorm(abs(Z))))
  r <- raster(x)
  out <- stack(list(S=setValues(r, Sm), s2 = setValues(r, s2), p = setValues(r, p)))
  if(calc_slope)  out$slope <- setValues(r, res$slope)
  attr(out,"timing") <-Sys.time() - t0
  out
}
