# Covariance matrix of a stack. Will be big...

devtools::load_all()


if(!exists("x")) x<-(s<-readRDS("test_stacks2.rds"))$trend

nh <- 2
mb <- microbenchmark::microbenchmark(times = 5,
   a <- stack_covariance(x, neighbourhood = nh, sparse = FALSE),
   b <- stack_covariance(x, neighbourhood = nh, sparse = TRUE, ver = 1),
   e <- stack_covariance(x, neighbourhood = nh, sparse = TRUE))

print(mb)

print(c(sum(a-b), sum(a-e), sum(b-e)))

k <- which(s$hastrend[]==1)
mis <- (is.na(x[][,1]))

image(Matrix(a[k,k]))#-e[k,k]))

cat("nnzero:", nnzero(b)/length(b), "\n")
# very sparse indeed
