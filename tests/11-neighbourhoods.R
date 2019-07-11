# Check neighbourhoods
library(raster)
devtools::load_all()

o <- readRDS("test_stacks2.rds")
x <- o$trend[[1]]


x1 <- focal_average(x, neighbourhood = 1)
x2 <- focal_average(x, neighbourhood = 2)

plot(stack(list(x,x1,x2)))


x <- o$trend
m <- contextual_mann_kendall(x, neighbourhood = 0)$S
m1 <- contextual_mann_kendall(x, neighbourhood = 1)$S
m2 <- contextual_mann_kendall(x, neighbourhood = 2)$S

plot(stack(list(m,m1,m2)))
