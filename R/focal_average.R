#' Raster Focal Statistics
#'
#' Average value in neighbourhood of each cell.
#'
#' @param x raster or rasterStack
#' @param neighbourhood 2: queen
#'
#' @details NA-neighbour values will be ignored and will not increase the denominator.
#'
#' Edge pixels have the correct denominators.
#'
#' The result will not fill in NA cells, so if the focal cell value is NA, it will be NA in the output.
#'
#' The input 'x' will be converted to a rasterStack with 'stack(x)' before calculations.
#'
#' @import raster
#' @export

focal_average <- function(x, ..., neighbourhood = 2) {
  if(!neighbourhood %in% c(2)) stop("Only '2' neighbourhood implemented.")
  # center first
  if( !canProcessInMemory(x) ) stop("'canProcessInMemory' return FALSE")
  t0 <- Sys.time()
  X <- values(stack(x))
  nr <- nrow(x)
  #browser()
  res <- c_focal_average(X, nr, neighbourhood = neighbourhood)
  x[,] <- if(ncol(X) == 1) res[,1] else res
  x
}
