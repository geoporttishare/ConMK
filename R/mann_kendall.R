#' Mann-Kendall Trend Test
#'
#' 1D. Assuming time
#'
#' @param x Time series, either a vector of values or n x 2 matrix of value,time pairs
#' @param est_beta Estimate linear coefficient as well?
#' @param check check the data. If you give nx2 matrix then this TRUE saves time.
#'
#' @useDynLib ConMK
#' @export

mann_kendall <- function(x, ..., est_beta = TRUE, check = TRUE) {
  if(check){
    x <- check_x(x)
    x <- na.omit(x)
  }
  y <- x[,1]
  time <- x[,2]
  res <- if(est_beta) c_mann_kendall_test_and_beta(y, time) else c_mann_kendall_test(y)
  S <- res$S
  s2 <- res$s2
  Z <- (S+ifelse(S<0, 1, -1) * 1)/sqrt(s2)
  p <- 2*(1 - abs(pnorm(abs(Z))))
  out <- data.frame(S=S, s2=s2, Z=Z, p = p)
  #
  if(est_beta) {
    Nstar <- sqrt(s2) * qnorm(1-0.05/2)
    N <- length(res$X)
    M1 <- (N - Nstar)/2
    M2 <- (N + Nstar)/2
    z <- quantile(res$X, c(M1/N, 0.5, M2/N), na.rm=TRUE)
    out$betahat <- z[2]
    out$beta_low5 <- z[1]
    out$beta_high5 <- z[3]
  }
  out
}


#' Mann-Kendall Trend Test for Stack
#'
#' 1D. Assuming time as layers
#'
#' @param x rasterStack
#' @param est_beta Estimate linear coefficient as well? Takes extra effort
#'
#' @useDynLib ConMK
#' @export

mann_kendall_stack <- function(x, ..., est_beta = TRUE) {
  if( !canProcessInMemory(x) ) stop("'canProcessInMemory' return FALSE")
  time <- 1:nlayers(x)
  fun <- function(y, ...) unlist(mann_kendall(cbind(y, time), check = FALSE, est_beta=est_beta))
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
