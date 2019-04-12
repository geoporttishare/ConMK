# Standard Mann-Kendall trend test

devtools::load_all()

# Test significance using Mann-Kendall:
MK <- function(y) {
  n <- length(y)
  k <- 0
  W <- vector(l=n*(n-1)/2, "numeric")
  for(i in 1:(n-1))for(j in (i+1):n)
    W[k<-k+1] <- sign(y[j]-y[i])
  S <- sum(W)
  s2 <- n*(n-1)*(2*n+5)/18
  Z <- (S+ifelse(S<0, 1, -1) * 1)/sqrt(s2)
  p <- 2*(1 - abs(pnorm(abs(Z))))
  data.frame(S=S, s2=s2, Z=Z, p = p)
}

# Example time 1d series
g <- seq(0, 1, l = 100) * 100
trend <- sin(g/50) # almost monotonic
# AR1:
library(Matrix)
rho <- 0.4
S <- bandSparse(length(g), k=0:1, symmetric = TRUE,
                diagonals = list(rep(1, length(g)), rep(rho, length(g)))) #exp(-.5*as.matrix(dist(g))^2/1^2)
A <- chol(S)
est <- NULL

for(i in 1:100){
  y0 <- c(A %*% rnorm(100))[[1]][,1]
  y <- y0 + trend * 5
  est <- rbind(est, (c(mann_kendall(y)$p, MK(y)$p)))
}
print(all.equal(est[,1], est[,2]))
print(mean(est[,1]<.05)) # should be too high thanks to autocorrelation




################################################
# 2D timeserieses as raster stack
# Stack:
if(0){
r <- stack( readRDS("test_stack.rds") )
p <- stackApply(r, 1, function(...) mann_kendall(...)$p )
plot(p, zlim = c(0,1)) # not much to see

# add trend
r1 <- susiapu::restack_layered(list(a=r))
r2 <- lapply(1:length(r1), function(i) r1[[i]] + i*0.05)
r3 <- susiapu::restack_layered(r2)[[1]]
p3 <- stackApply(r3, 1, function(...) mann_kendall(...)$p )
plot(p3, zlim = c(0,1)) # not much to see
}

