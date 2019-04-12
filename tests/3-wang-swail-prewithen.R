# Wann Swail prewhitening testing


devtools::load_all()

# series with autocorrelation, AR(1)

library(Matrix)
g <- 1:10
rho <- 0.1
S <- bandSparse(length(g), k=0:1, symmetric = TRUE,
                diagonals = list(rep(1, length(g)), rep(rho, length(g))))
A <- chol(S)
est <- NULL

if(1){
y <- as.numeric(A %*% rnorm(length(g)))
trend <- sin(g/200)
x <- y + trend * 1
#plot(y, type="l", ylim =c(-1,1)*3)

#rhoo <- autocor_lag1(x)

t0 <- system.time( res <- wang_swail_prewhiten_1d(x)  )
t1 <- system.time( resc <- wang_swail_prewhiten_1d(x, useC=TRUE) )

print(all.equal(res, resc ))

print(rbind(r=unlist(res[-1]), c=unlist(resc[-1])))

print(rbind(r=t0, c=t1))

}


if(0){
# Check false positive if trend is found:


z <- sapply(1:2000, function(i) {
  y <- as.numeric(A %*% rnorm(length(g)))
  # without withening
  m <- mann_kendall(y)$p
  # with whitening
  mw <- mann_kendall(wang_swail_prewhiten_1d(y, itmax = 20))$p
  c(mk=m, mkw=mw)
})

print( rowMeans(z < 0.05) )
}
