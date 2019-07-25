#' Cox-Stuart Trend Test
#'
#' @param x A vector of time series values in order
#' @param check check the data. If you give nx2 matrix then this FALSE saves time.
#' @param split <= 0.5. Size of the beginning/end of series to use in comparison.
#' @param alternative "two.sided", "less" (downwards) or "greater" (upwards trend)
#' @details
#' The test is done two tailed, no up/down direction.
#' @importFrom stats pbinom
#' @export
cox_stuart <- function(x, check = TRUE, split = .5, alternative = "two.sided") {
  if(check){
    x <- check_x(x)
    y <- x[,1]
  }else y <- as.numeric(x)
  y <- na.omit(y)
  n <- length(y)
  if(n < 3) return(data.frame(S=NA, p.value=NA))
  l <- floor(split * n) # items to use
  v <- sign( y[1:l] - y[(n-l+1):n] )
  v <- v[v != 0]# eliminate ties
  nn <- length(v)
  S <- table(v) # count -1 +1
  if(alternative == "two.sided") p <-  2 * pbinom(min(S), nn, .5) # two tailed
  else if(alternative == "less") p <- pbinom(min(S), nn, .5)
  else if(alternative == "greater")  p <- pbinom(max(S), nn, .5)
  else stop("'alternative' not recognized")
  #browser()
  data.frame(S = max(S), p.value = p)
}



#' Cox-Stuart Trend Test For a Stack
#'
#'
#' @param x rasterStack
#' @param split See \link{cox_stuart}
#' @param ... ignored
#' @details Stack's layers are taken as time. Each cell's timeseries is passed on to \link{cox_stuart}.
#'
#' @export

cox_stuart_stack <- function(x, split = 0.5,  ...) {
  if( !canProcessInMemory(x) ) stop("'canProcessInMemory' return FALSE")
  fun <- function(y, ...) unlist(cox_stuart(y, check = FALSE, split = split))
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
