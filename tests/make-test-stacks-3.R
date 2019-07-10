# Generate test stacks 3: Gaussian RF, larger size


library(raster)
library(RandomFields)

# Stacks 3.*: somewhat larger
nc <- 200
nr <- 200
nt <- 20  # time. cut to length
Nx <- nc * nr
N <- nc * nr * nt
r0 <- raster(nrows = nr, ncols = nc, crs = CRS("+init=epsg:3067"), xmn=0, xmx=1,
             ymn = 0, ymx=nr/nc)
r0 <- setValues(r0, rep(0, Nx))
s0 <- stack(rep(list(r0), nt))

# Mask for NA
# Combine this with excursion set-ting the field to get more ragged NA regions.
xy <- coordinates(r0)
R2 <- 0.4^2
c0 <- c(1,0)
themNA <- which(colSums((c0-t(xy))^2) < R2 )
mask <- r0
mask[themNA] <- NA

########################################################
# Basic field,
RFoptions(spConform=FALSE)
s2 <- 1
mod <- RMgauss(s2, 0.1)
v1 <- 2.0 + RFsimulate(mod, x = xy[,1], y = xy[,2])
mask1 <- mask
mask1[ v1 < 0 ] <- NA
r1 <- mask(setValues(r0, v1), mask1)

s <- 0.3 * sqrt(s2)
v1n <- v1 + rnorm(N, 0, s)
s1n <- mask(setValues(s0, v1n), mask1)

#
# noise with strong AR(1) noise
v1 <- v1n
A <- Matrix::bandSparse(nt, k = 0:1, diagonals = list(rep(1, nt), rep(0.5, nt)), symmetric=TRUE)
U <- Matrix::chol(A)
s1a <- mask(calc(setValues(s0, v1), function(v,...) as.numeric(U%*%v) ), mask1)
#acf(s1a[1][1,])


# noise with trendy area in the middle
D <- distanceFromPoints(r0, c(.5,.5 * nr/nc))
them <- which(values(D) < 0.1)
hastrend <- r0 * 0
hastrend[them] <- 1
hastrend <- mask(hastrend, mask1)
s1t <- s1n
for(k in them) s1t[k] <- s1t[k] + 1:nt * 1/nt
#plot(s1t[k][1,])
#plot(s1t,1:5*4)

# AR + trend
s1at <- s1a
for(k in them) s1at[k] <- s1at[k] + 1:nt * 1/nt
#plot(s1at[k][1,], ylim = c(0,4))
# plot(s1at,1:5*4)

###
out <- list(noise=s1n, ar=s1a, trend=s1t, artrend=s1at, hastrend = hastrend)


par(mfrow=c(4,1))
for(n in names(out)[-5])
  plot(out[[n]][k][1,], main =n)




#### store
saveRDS(out, "test_stacks3.rds")
test_stacks3 <- out
save(test_stacks3, file="data/test_stacks3.rda", compress = TRUE)

