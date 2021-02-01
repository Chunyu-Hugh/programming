library(GA)
library(foreach)
library(iterators)

f <- function(x) 1/(x^2+5*x+4)

ga.result <- ga("real-valued"
                ,fitness = f
                ,lower = -5
                , upper = 0
                , nBits = 5
                , popSize = 20
                , pcrossover = 1
                , pmutation = 0.001
                ,maxiter = 200
                
)
summary(ga.result)
plot(ga.result)