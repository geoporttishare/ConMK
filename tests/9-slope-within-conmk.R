# Check the slope calculator

devtools::load_all()
o <- readRDS("test_stacks2.rds")

x <- o$trend

t0 <- system.time( r0 <- contextual_mann_kendall(x, neighbourhood = 2, calc_slope = FALSE) )
t1 <- system.time( r  <- contextual_mann_kendall(x, neighbourhood = 2, calc_slope = TRUE ) )

print(rbind(t0, t1))


rs <- stackApply(x, 1, theil_sen)
r0$slope <- rs


print( all.equal(r, r0) )

plot(r0)
