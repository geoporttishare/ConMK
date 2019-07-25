# Pettitt's changepoint tests

devtools::load_all()


x <- (ss <- readRDS("test_stacks2.rds"))$trend
Q <- pettitt_stack(x)

QQ <- split_calc_wrapper(x, 2, 2, buffer = c(0,0), dbg=TRUE, fun = pettitt_stack)

# plot ones
p <- mask(Q$p.approx, is.na(x$layer.1))
k <- which(p[] < 0.01)

par(mfrow=c(4,4))
for(i in k[1:16]) {
  plot(NA, xlim = c(0, nlayers(x)), ylim = range(x[],na.rm=T))
  j <- Q$idx[i]
  v <- x[][i,]
  n <- length(v)
  lines(1:j, v[1:j])
  points(j, v[j], pch = 19, cex = .5)
  lines(j:n, v[j:n], col = 2)
}
