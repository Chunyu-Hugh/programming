library(mcga)
f <- function(x) {
  return(-exp(1)
   + 20 * exp(-0.2 * sqrt((x[1] ^ 2 + x[2] ^ 2) / 2))
    + exp((cos(2 * pi * x[1]) + cos(2 * pi * x[2])) / 2)
   )
}

m <- mcga(popsize = 6,
             chsize = 2,
             minval = 0,
             maxval = 4,
           maxiter = 1000,
            crossprob = 1.0,
            mutateprob = 0.02,
            evalFunc = f)