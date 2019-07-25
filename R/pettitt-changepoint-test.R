#' Pettitt's Test For a Change-Point
#'
#' Test for a single changepoint in a timeseries.
#'
#' @param x A vector of time series values in order
#' @param check check the data. If you give nx2 matrix then this FALSE saves time.
#'
#' @export
pettitt <- function(x, check = TRUE) {
  if(check){
    x <- check_x(x)
    x <- na.omit(x)
    y <- x[,1]
  }else y <- as.numeric(x)
  r <- rank(y)
  n <- length(r)
  U <- 2 * cumsum(r) - 1:n * (n+1)
  U <- abs(U)
  K <- U[Ki <- which.max(U)]
  L <- n
  p <- min( 2 * exp(-6 * K^2/(L^3 + L^2)), 1)
  data.frame(Ustar = K, p.approx = p, idx = Ki)
}



#' Pettitt's Change-Point Test for Stack
#'
#' Assuming time as layers
#'
#' @param x rasterStack
#' @param ... passed on to \link{pettitt}
#' @details Stack's layers are taken as time. Each cell's timeseries is passed on to \link{pettitt}.
#' @export

pettitt_stack <- function(x, ...) {
  if( !canProcessInMemory(x) ) stop("'canProcessInMemory' return FALSE")
  fun <- function(y, ...) unlist(pettitt(y, check = FALSE))
  X <- values(x)
  # Calculate
  V <- apply(X, 1, fun)
  # compile
  r <- raster(x)
  rl <- lapply(rownames(V), function(i) setValues(r, V[i,]))
  s <- stack(rl)
  names(s) <- rownames(V)
  s
}
