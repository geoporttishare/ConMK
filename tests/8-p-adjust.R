# Test the p.adjust-function

devtools::load_all()

data("test_stacks2")

st <- test_stacks2$artrend

res <- contextual_mann_kendall(st)

x <- res$p

# Holm-Bonferroni by hand
a <- 0.05
pv <- sort(values(x))
e <- a/( length(pv) - 1:length(pv) + 1 )
th <- pv[ 1L+sum(!(pv>e)) ]

xp <- p.adjust_raster(x, method = "holm")

plot(xp<.05)
