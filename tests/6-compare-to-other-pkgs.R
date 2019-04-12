# Compare to other implementations


devtools::load_all(".")

if(!exists("y")){
# a series

library(Matrix)
g <- 1:10
rho <- 0.1
S <- bandSparse(length(g), k=0:1, symmetric = TRUE,
                diagonals = list(rep(1, length(g)), rep(rho, length(g))))
A <- chol(S)
est <- NULL

y <- as.numeric(A %*% rnorm(length(g)))
trend <- sin(g/10)
x <- y + trend * 2
}

v <- y

### mann-kendall
mann_kendall(v)
modifiedmk::mkttest(v)

mann_kendall(wang_swail_prewhiten_1d(v)$W)
modifiedmk::pwmk(v)
modifiedmk::pbmk(v)
modifiedmk::mmkh(v)
modifiedmk::mmky(v)
modifiedmk::tfpwmk(v)
# pretty much the same
