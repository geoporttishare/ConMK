# Generate test stacks


library(raster)


# Stacks 1.*: small
nc <- 60
nr <- 40
nt <- 20  # cut to length
Nx <- nc * nr
N <- nc * nr * nt
r0 <- raster(nrows = nr, ncols = nc, crs = CRS("+init=epsg:3067"), xmn=0, xmx=1,
             ymn = 0, ymx=nr/nc)
r0 <- setValues(r0, rep(0, Nx))
s0 <- stack(rep(list(r0), nt))

# Mask for NA
xy <- coordinates(r0)
R2 <- 0.2^2
c0 <- c(1,0)
themNA <- which(colSums((c0-t(xy))^2) < R2 )
mask <- r0
mask[themNA] <- NA


##
s <- 0.1
# just noise
v1n <- rnorm(N, 0, s)
s1n <- mask(setValues(s0, v1n), mask)

#
# noise with strong AR(1) noise
v1 <- v1n
A <- Matrix::bandSparse(nt, k = 0:1, diagonals = list(rep(1, nt), rep(0.5, nt)), symmetric=TRUE)
U <- Matrix::chol(A)
s1a <- mask(calc(setValues(s0, v1), function(v,...) as.numeric(U%*%v) ), mask)
#acf(s1a[1][1,])

# noise with trendy area in the middle
D <- distanceFromPoints(r0, c(.5,.5 * nr/nc))
them <- which(values(D) < 0.1)
hastrend <- mask(r0 * 0, mask)
hastrend[them] <- 1
s1t <- s1n
for(k in them) s1t[k] <- s1t[k] + 1:nt * .3/nt
#plot(s1t[k][1,])

# AR + trend
s1at <- s1a
for(k in them) s1at[k] <- s1at[k] + 1:nt * .3/nt
#plot(s1at[k][1,])

###
out <- list(noise=s1n, ar=s1a, trend=s1t, artrend=s1at, hastrend = hastrend)


par(mfrow=c(4,1))
for(n in names(out)[-5])
  plot(out[[n]][k][1,], main =n)

#### store
saveRDS(out, "test_stacks2.rds")
test_stacks2 <- out
save(test_stacks2, file="data/test_stacks2.rda", compress = TRUE)



