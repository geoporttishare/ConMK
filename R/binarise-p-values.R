#' Adjust p-values in a Raster
#'
#' Simple wrapper for the p.adjust-function.
#'
#' @param x A raster of p-values
#' @param alpha nominal level, like 0.05 or 0.01
#' @param method Type of discovery
#'
#' @details See \link{p.adjust} for the method parameters.
#'
#' @return a raster-object
#'
#'
#'
#' @export

p.adjust_raster <- function(x, method = c("none", "holm")) {
  #
  pv <- x[]
  ok <- !is.na(pv)
  e <- calc(x, p.adjust, method = method)
  o <- lapply(method, function(m) {
                r <- raster(x)
                v <- p.adjust(pv[ok], method = m)
                r[ok] <- v
                r
            })
  out <- stack(o)
  names(out) <- method
  out
}
