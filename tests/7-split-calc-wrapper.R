# Test the split-calc wrapper

devtools::load_all()

data("test_stacks3")

x <- stack(test_stacks3$trend, layers = 1:9*2)

# direct
#a <- contextual_mann_kendall(x)

# divide and conquer
k <- split_calc_wrapper(x, nx=3, ny=3, buffer = c(1,1), fun = contextual_mann_kendall, dbg=TRUE)
b <- k$result

all.equal(a[],b[])
plot(a$s2 - b$s2)

