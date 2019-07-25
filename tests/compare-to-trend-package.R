# Comparisons to 'trend' package

library(trend)
devtools::load_all()

data("maxau")
Q <- maxau[,"s"]
mk.test(Q)
sens.slope(Q)
mann_kendall(as.numeric(Q))
# ok


csmk.test(nottem)
mann_kendall(c(nottem))

pettitt.test(Q)
pettitt(c(Q))
