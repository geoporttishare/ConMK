# Local averaging of a raster layer, accounting for missing and edge pixels

devtools::load_all()
o <- readRDS("test_stacks2.rds")

x <- o$trend


if(1) {
r0 <- contextual_mann_kendall(x, neighbourhood = 2)
s0 <- contextual_mann_kendall(x, neighbourhood = 0)
S1 <- focal_average(s0$S)
print(  all.equal(r0$S, S1)  )
plot(stack(list(s0$S, r0$S, S1)))
}

if(0){
  y <- focal_average(x)
  plot(stack(x,y), c(1,21))
}

