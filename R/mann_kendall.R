#' Mann-Kendall Trend Test
#'
#' 1D. Assuming time
#'
#' @param x Time series, either a vector of values or n x 2 matrix of value,time pairs
#' @param ... ignored
#' @param calc_slope Estimate linear coefficient as well?
#' @param est_beta calc_slope (legacy)
#' @param check check the data. If you give nx2 matrix then this FALSE saves time.
#'
#' @useDynLib ConMK
#' @importFrom stats na.omit pnorm qnorm
#' @export

mann_kendall <- function(x, ..., est_beta = TRUE, calc_slope = est_beta, check = TRUE) {
  if(check){
    x <- check_x(x)
    x <- na.omit(x)
  }

  y <- x[,1]
  time <- x[,2]
  res <- if(calc_slope) c_mann_kendall_test_and_beta(y, time) else c_mann_kendall_test(y)
  S <- res$S
  s2 <- res$s2
  Z <- (S+ifelse(S<0, 1, -1) * 1)/sqrt(s2)
  p <- 2*(1 - abs(pnorm(abs(Z))))
  out <- data.frame(S=S, s2=s2, Z=Z, p = p)
  #
  if(calc_slope) {
    Nstar <- sqrt(s2) * qnorm(1-0.05/2)
    N <- length(res$X)
    M1 <- (N - Nstar)/2
    M2 <- (N + Nstar)/2
    z <- quantile(res$X, c(M1/N, 0.5, M2/N), na.rm=TRUE)
    out$slope <- z[2]
    out$slope_low5 <- z[1]
    out$slope_high5 <- z[3]
  }
  out
}


#' Mann-Kendall Trend Test for Stack
#'
#' 1D. Assuming time as layers.
#'
#' @param x rasterStack
#' @param ... passed on to \link{contextual_mann_kendall}
#' @param neighbourhood default: 0. See \link{contextual_mann_kendall}
#'
#' @details Just a wrapper for \link{contextual_mann_kendall} with neighbourhood=0 i.e. no
#' smoothing applied.
#' @useDynLib ConMK
#' @export

mann_kendall_stack <- function(x, ..., neighbourhood = 0) {
  contextual_mann_kendall(x, ..., neighbourhood = neighbourhood)
  # if( !canProcessInMemory(x) ) stop("'canProcessInMemory' return FALSE")
  # time <- 1:nlayers(x)
  # fun <- function(y, ...) unlist(mann_kendall(cbind(y, time), check = FALSE, est_beta=est_beta))
  # X <- values(x)
  # # Calculate
  # V <- apply(X, 1, fun)
  # # compile
  # r <- raster(x)
  # rl <- lapply(rownames(V), function(i) setValues(r, V[i,]))
  # s <- stack(rl)
  # names(s) <- rownames(V)
  # s
}
