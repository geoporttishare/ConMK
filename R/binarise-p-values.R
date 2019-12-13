#' Adjust p-values in a Raster
#'
#' Simple wrapper for the p.adjust-function.
#'
#' @param x A raster of p-values
#' @param method Type of discovery
#'
#' @details See \link{p.adjust} for the method parameters.
#'
#' @return a raster-object
#'
#' @importFrom stats p.adjust
#' @export

p.adjust_raster <- function(x, method = c("none", "fdr")) {
  #
  #e <- calc(x, p.adjust, method = method)
  pv <- x[]
  ok <- !is.na(pv)
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
