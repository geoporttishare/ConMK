# Wang-Swail prewhitening a stack pixelwise

devtools::load_all()

library(raster)

# test stack
o <- readRDS("test_stacks2.rds")
x <- r <- stack( o$trend )

t0 <- system.time( rW<- wang_swail_prewhiten_stack(r, useC=FALSE) )
t1 <- system.time( rWc<-wang_swail_prewhiten_stack(r, useC=TRUE) )


print(rbind(r=t0, c=t1))

print(all.equal(rW, rWc))



if(1){
## mann_kendall test for the whitened:
fun <- function(v, ...) {
  if(sum(is.na(v))>1) return(NA)
  v <- na.omit(v)
  mann_kendall(v, est_beta = FALSE)$p
}
mk <- stackApply(r, 1, fun)
mkw <- stackApply(rW, 1, fun)
par(mfrow=c(2,1))
plot(mk<.05, main="mk < 0.05")
plot(mkw<.05, main="whitened")
}

# mighty noisy
