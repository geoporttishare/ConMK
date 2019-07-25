# Cox Stuart trend test

devtools::load_all()


data("maxau", package = "trend")
Q <- maxau[,"s"]
print( trend::cs.test(Q) )
print( cox_stuart(c(Q),split = .5) )




# ok...

x <- (ss <- readRDS("test_stacks2.rds"))$trend

mk <- mann_kendall_stack(x, est_beta = FALSE)
cs <- cox_stuart_stack(x, split = .5)

plot(stack(mk$p <0.05, cs$p.value < 0.05))

e <- split_calc_wrapper(x, 2,2, fun = function(z,...)
  stackApply(z, 1, fun = function(v,...) if(all(is.na(v))) NA else trend::cs.test(v)$p.value))

plot(stack(e$result<0.05, cs$p.value< 0.05 ))
