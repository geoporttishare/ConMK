#' Contextual Mann-Kendall Trend Test
#'
#' 2D raster version
#'
#' @param x Raster stack, each layer denoting one timestep. Equidistant steps of time 1.
#' @param neighbourhood 2:queen, 1:rook, 0:none (classical M-K test)
#'
#' @references
#' Neeti, N. and Eastman J.R. (2011) A Contextual Mann-Kendall Approach for the Assesment of Trend Significance in Image Time Series, \emph{Transactions in GIS}
#' @useDynLib ConMK
#' @import raster
#' @export

contextual_mann_kendall <- function(x, ..., neighbourhood = 2) {
  # center first
  if( !canProcessInMemory(x) ) stop("'canProcessInMemory' return FALSE")
  x <- x - mean(x)
  X <- values(x)
  nr <- nrow(x)
  # c_contextual_mann_kendall(v, nr)
  out <- c_contextual_mann_kendall( X , nr , neigh = neighbourhood)
  V <- out[[3]]
  Sm <- V[,1]
  s2 <- V[,2]
  D <- (Sm > 0) * 1
  Z <- (Sm+(1-2*D))/sqrt(s2)
  p <- 2*(1 - abs(pnorm(abs(Z))))
  r <- x[[1]]
  stack(list(S=setValues(r, Sm), s2 = setValues(r, s2), p = setValues(r, p)))
}
