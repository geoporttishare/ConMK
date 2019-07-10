# Dev and test contextual_mann_kendall

devtools::load_all()

library(raster)

# test stack
if(!exists("v")){
o <- readRDS("test_stacks2.rds")
v <- values(s<-o$noise)
r <- s[[1]]
}
# indexing
if(0) {
  r <- s[[1]]
  tr <- t(r)
  d <- dim(tr)
  i <- sample(prod(d[1:2])-1,1)
  a <- neighbour_cells_queen_row_col(i, d[1], d[2])
  b <- neighbour_cells_queen_col_row(i, d[2], d[1])
  e <- adjacent(r,i+1,8)[,2]-1
  print(c(all.equal(sort(a),sort(b)), all.equal(sort(a),sort(e))))
}
# forward indexing
if(0) {
  r <- s[[1]]
  tr <- t(r)
  d <- dim(tr)
  i <- sample(prod(d[1:2])-1,1)
  a <- forward_neighbour_cells_queen_row_col(i, d[1], d[2])
  b <- forward_neighbour_cells_queen_col_row(i, d[2], d[1])
  e <- adjacent(r,i+1,8)[,2]-1;e <- e[e>i]
  print(c(all.equal(sort(a),sort(b)), all.equal(sort(a),sort(e))))
}

# variances and sum of covariances
if(0){
  ss <- s - mean(s)
  v <- values(ss)
  nr <- dim(ss)[1]
  out <- c_contextual_mann_kendall(v, nr)
  V <- out[[1]]
  e <- stackApply(ss, 1, var) #
  #e <- sum(ss^2)/(nlayers(ss)-1)
  s0 <- apply(v[1:4,], 1, var)
  print(cbind(my=V[1:4,1], sp=e[1:4], s0))
  print( all.equal(V[,1], values(e)) )

  ##
  S <- out[[2]]
  Sm <- V[,4]
  print( cbind( S[1:4],Sm[1:4]) )
  plot( stack(setValues(r, S), setValues(r, Sm), focal(setValues(r,S), matrix(1, nc=3,nr=3)/9)) , zlim = c(-1,1)*80)
}

# smoothed statistic and its s2
if(0){
  s<-o$trend
  ss <- s - mean(s)
  v <- values(ss)
  nr <- dim(ss)[1]
  out <- c_contextual_mann_kendall(v, nr)
  V <- out[[3]]
  plot( stack(S <- setValues(r, V[,1]), s2 <- setValues(r, V[,2])) )
  D <- (S<0) * 1
  Z <- (S+(1-2*D))/sqrt(s2)
  p <- calc(Z, function(v) 2*(1 - abs(pnorm(abs(v)))))
  plot(stack(p, p < 0.05))
}

# first version of wrapper
if(0){
  x <- o$trend
  res <- contextual_mann_kendall(x)
  res$p_th <- res$p < 0.05
  plot(res)
}



# check if the non-contextual stack version works
if(1) {
  x <- o$trend
  fun<-function(v,...) if(sum(is.na(v))>1) return(NA) else mann_kendall(v, est_beta = FALSE)$p
  t0 <- system.time(   a <- stackApply(x, 1, fun)  )
  t1 <- system.time(   aa <- mann_kendall_stack(x, est_beta=FALSE)$p )
  t2 <- system.time(   b <- (s <- contextual_mann_kendall(x, neighbourhood = 0))$p )

  (all.equal(a,b))
  print(rbind(t0,t1,t2))


}
