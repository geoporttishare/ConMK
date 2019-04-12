# Dev doodles

# Example time 1d series
g <- seq(0, 1, l = 100) * 100
S <- exp(-.5*as.matrix(dist(g))^2/1^2)
y0 <- c(chol(S) %*% rnorm(100))
trend <- sin(g/50)
y <- y0 + trend*4
#
# Theil-Sen slope estimator
theil_sen <- function(y){
  n <- length(y)
  k <- 0
  W <- vector(l=n*(n-1)/2, "numeric")
  for(i in 1:(n-1))for(j in (i+1):n)
    W[k<-k+1] <- (y[j]-y[i])/(j-i)
  median(W, na.rm=TRUE)
}


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



plot(y)
o <- modifiedmk::tfpwmk(y)
o2 <- MK(y)
print(o)
print(o2)


### raster stack


r <- stack( readRDS("test_stack.rds") )
p <- stackApply(r, 1, sd )
p <- calc(r, sd)
plot(p) # not much to see

## Smoorhing
filter <- focalWeight(r, 5000, "rect")
mask  <- is.na(r[[1]])
ps <- stack(lapply(as.list(r), function(l){
  e <- focal(l, w=filter, pad=T, fun = function(x,...) {i <- !is.na(x); sum(x[i])/sum(filter[i])} )
  #e <- focal(l, w=filter, pad=T )
  e[mask] <- NA
  e
}))
#plot(stack(r[[1]],ps))
plot(ps)
