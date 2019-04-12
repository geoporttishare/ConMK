# Test Theil-Sen slope estimator

devtools::load_all(".")

# Example time 1d series
est <- NULL
for(i in 1:100){

g <- seq(0, 1, l = 100) * 100
S <- exp(-.5*as.matrix(dist(g))^2/1^2)
y0 <- c(chol(S) %*% rnorm(100))
trend <- sin(g/50) # almost monotonic
y <- y0 + trend*4
#
# Theil-Sen slope estimator
theil_sen_b <- function(y){
  n <- length(y)
  k <- 0
  W <- vector(l=n*(n-1)/2, "numeric")
  for(i in 1:(n-1))for(j in (i+1):n)
    W[k<-k+1] <- (y[j]-y[i])/(j-i)
  median(W, na.rm=TRUE)
}

est <- rbind(est, (c(theil_sen(y), theil_sen_b(y))))
}
(all.equal(est[,1], est[,2]))


