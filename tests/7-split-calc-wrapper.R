# Test the split-calc wrapper

devtools::load_all()

data("test_stacks3")
x <- stack(test_stacks3$trend)
#x <- crop(  x , extent(0, 1, 0, .5))


# direct
if(!exists("a")) a <- contextual_mann_kendall(x)


########################
# Check splitting
if(0){

# divide and conquer
k <-  split_calc_wrapper(x, nx=3, ny=4, buffer = c(1,1), fun = contextual_mann_kendall, dbg=TRUE, bycell=F)
k2 <- split_calc_wrapper(x, nx=3, ny=4, buffer = c(1,1), fun = contextual_mann_kendall, dbg=TRUE, bycell=T)
b <- k2$result

(all.equal(a[],b[]))
par(mfrow=c(2,2))
#plot(a$s2 - b$s2)
#plot(k$result$s2 - k2$result$s2)
#lapply(k$extents, plot, add=T, col = rgb(1,0,0,.5))
#lapply(k2$extents, plot, add=T, col = rgb(0,0,1,.5))
plot(a$s2)
plot(k$result$s2)
plot(k2$result$s2)

print(c(direct=attr(a,"timing"), ex=k$timing$stitch, rc= k2$timing$stitch))
}


########################
# Check snow
if(1){
beginCluster(n=3)
k <-  split_calc_wrapper(x, nx=3, ny=4, buffer = c(1,1), fun = contextual_mann_kendall, dbg=TRUE)
b <- k$result
(all.equal(a[],b[]))
print(c(direct=attr(a,"timing"), ex=k$timing$stitch))
endCluster()
}
